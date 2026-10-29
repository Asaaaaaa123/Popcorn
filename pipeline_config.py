#!/usr/bin/env python3
"""
pipeline_config.py
==================
One configuration object for bone_analysis.py, bone_figures.py and the web application in webapp/.

[APP] Added 2026-10-28, when the pipeline was turned into a configurable web application. DEFAULT_CONFIG reproduces,
value for value, the constants that used to be hard-coded in bone_analysis.py and bone_figures.py (the manuscript
setup: LL-2-1 control and RL-1-5 tumor-bearing tibia of one NSG mouse, 16 um bins), so a run with the defaults is the
manuscript analysis.

A configuration is plain JSON with these sections:
  data_root, samples, reference_sample, bin_size, qc_bin_size      where the Space Ranger outputs are
  qc, clustering, scoring, roi, compartments, uncoupling,
  proximity, moran, sensitivity                                     analysis parameters
  signatures, classes, primary, programs, receptor_panel,
  adipocyte_genes, ubiquitous                                       gene sets (Table 1) and how they are used
  steps, figures, seed, deterministic                               what runs and how it is drawn

Public helpers: default_config(), merge_with_defaults(), validate_config(), effective_config(), load_config(),
save_config(), resolve_bin_dir(), available_bin_sizes(), sample_by_id(), bin_um(), slugify(), Progress.
"""
from __future__ import annotations

import copy
import json
import math
import re
import sys
import time
from pathlib import Path

CONFIG_VERSION = 1

# ============================================================================
# DEFAULTS (the manuscript setup)
# ============================================================================

# Gene signatures (Table 1). Per signature:
#   label        name drawn in figures (bone_figures.lbl(): "_" -> " ", neutrophil -> granulocyte)
#   table_label  name written to Table 1
#   color        figure color (bone_analysis.PALETTE)
#   class        tissue class of Figure 3 (bone_figures.CLASS); "removed" = assigned but not drawn
#   canonical    canonical markers listed in Table 1
#   core         lineage-restricted genes; at least `scoring.min_core` must be detected for a positive call
#   extended     genes scored
#   competes     a bin that scores higher on one of these is not called positive
#   excluded     genes deliberately left out, with the reason (Table 1)
_SIGNATURES = {
    "osteoclast": {
        "label": "osteoclast", "table_label": "osteoclast", "color": "#377EB8", "class": "bone lineage",
        "canonical": ["Atp6v0d2", "Ctsk"],
        "core": ["Acp5", "Ctsk", "Calcr", "Atp6v0d2", "Dcstamp", "Ocstamp"],
        "extended": ["Acp5", "Ctsk", "Calcr", "Mmp9", "Atp6v0d2", "Dcstamp", "Ocstamp", "Nfatc1"],
        "competes": ["erythroid", "neutrophil"],
        "excluded": {"Car2": "abundant in erythroid precursors"}},
    "osteoblast": {
        "label": "osteoblast", "table_label": "osteoblast", "color": "#E41A1C", "class": "bone lineage",
        "canonical": ["Col1a1", "Sp7"],
        "core": ["Runx2", "Sp7", "Alpl", "Ibsp"],
        "extended": ["Runx2", "Sp7", "Alpl", "Ibsp"],
        "competes": ["fibroblast", "chondrocyte"],
        "excluded": {"Col1a1": "pan-fibroblast", "Col1a2": "pan-fibroblast", "Spp1": "major macrophage transcript",
                     "Bglap": "absent from the probe set", "Bglap2": "absent from the probe set"}},
    "osteocyte": {
        "label": "osteocyte", "table_label": "osteocyte", "color": "#984EA3", "class": "bone lineage",
        "canonical": [],
        "core": ["Sost", "Dmp1", "Mepe"],
        "extended": ["Sost", "Dmp1", "Mepe", "Fgf23", "Pdpn"],
        "competes": [],
        "excluded": {"Phex": "not detected in this probe set"}},
    "fibroblast": {
        "label": "fibroblast", "table_label": "fibroblast", "color": "#FB9A99", "class": "stroma",
        "canonical": [],
        "core": ["Pdgfra", "Dcn", "Lum"],
        "extended": ["Pdgfra", "Pdgfrb", "Fap", "Thy1", "Acta2", "Dcn", "Lum", "Col3a1", "Postn", "Fbn1"],
        "competes": [],
        "excluded": {}},
    "macrophage": {
        "label": "macrophage", "table_label": "macrophage", "color": "#DA70D6", "class": "hematopoietic marrow",
        "canonical": [],
        "core": ["Adgre1", "Csf1r", "C1qa", "C1qb", "Mrc1"],
        "extended": ["Adgre1", "Csf1r", "C1qa", "C1qb", "C1qc", "Mrc1", "Cd68", "Trem2", "Ctss", "Lyz2", "Apoe"],
        "competes": ["neutrophil"],
        "excluded": {}},
    # "neutrophil" is the granulocyte signature of the manuscript (marrow granulopoiesis)
    "neutrophil": {
        "label": "granulocyte", "table_label": "granulocyte", "color": "#FF7F00", "class": "hematopoietic marrow",
        "canonical": [],
        "core": ["S100a8", "S100a9", "Retnlg", "Ngp", "Ltf"],
        "extended": ["S100a8", "S100a9", "Retnlg", "Ly6g", "Mpo", "Elane", "Ngp", "Ltf", "Camp", "Lcn2"],
        "competes": [],
        "excluded": {}},
    "erythroid": {
        "label": "erythroid", "table_label": "erythroid", "color": "#A65628", "class": "hematopoietic marrow",
        "canonical": [],
        "core": ["Hba-a1", "Hbb-bs", "Alas2"],
        "extended": ["Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2", "Slc4a1", "Gypa", "Klf1", "Blvrb", "Dmtn"],
        "competes": [],
        "excluded": {}},
    "megakaryocyte": {
        "label": "megakaryocyte", "table_label": "megakaryocyte", "color": "#7FC97F", "class": "hematopoietic marrow",
        "canonical": [],
        "core": ["Pf4", "Itga2b", "Vwf"],
        "extended": ["Pf4", "Itga2b", "Vwf", "Ppbp", "Gp9", "Gp1bb", "Thbs1"],
        "competes": [],
        "excluded": {}},
    "adipocyte": {
        "label": "adipocyte", "table_label": "adipocyte", "color": "#FFD92F", "class": "stroma",
        "canonical": [],
        "core": ["Adipoq", "Plin1", "Cidec"],
        "extended": ["Adipoq", "Plin1", "Cidec", "Lep", "Retn", "Cfd", "Plin4", "Lipe", "Pparg"],
        "competes": ["skeletal_muscle"],
        "excluded": {"Car3": "skeletal muscle marker", "Fabp4": "endothelium and macrophages"}},
    # Endothelium competes in the tissue assignment so that vessel-rich bins are not assigned to fibroblast or
    # osteoblast; like skeletal muscle, nerve and smooth muscle it is not drawn in Figure 3 (class "removed").
    "endothelium": {
        "label": "endothelium", "table_label": "endothelial", "color": "#17BECF", "class": "removed",
        "canonical": [],
        "core": ["Pecam1", "Cdh5", "Emcn"],
        "extended": ["Pecam1", "Cdh5", "Emcn", "Flt1", "Kdr", "Cldn5"],
        "competes": [],
        "excluded": {}},
    "chondrocyte": {
        "label": "chondrocyte", "table_label": "chondrocyte", "color": "#33A02C", "class": "chondrocyte",
        "canonical": [],
        "core": ["Col2a1", "Acan", "Col10a1"],
        "extended": ["Col2a1", "Acan", "Col10a1", "Sox9", "Comp", "Matn3"],
        "competes": [],
        "excluded": {"Mmp13": "also induced in remodeling tissue; circular when used to define the growth-plate "
                              "compartment"}},
    "skeletal_muscle": {
        "label": "skeletal muscle", "table_label": "skeletal muscle", "color": "#8C8C8C", "class": "removed",
        "canonical": [],
        "core": ["Ttn", "Ckm", "Acta1"],
        "extended": ["Ttn", "Tnnt3", "Acta1", "Ckm", "Tpm2", "Mylpf", "Pvalb", "Myh1", "Tnnc2", "Atp2a1", "Car3"],
        "competes": [],
        "excluded": {}},
    "smooth_muscle": {
        "label": "smooth muscle", "table_label": "vascular smooth muscle", "color": "#D2B48C", "class": "removed",
        "canonical": [],
        "core": ["Myh11", "Cnn1"],
        "extended": ["Myh11", "Tagln", "Cnn1", "Myl9", "Acta2"],
        "competes": [],
        "excluded": {}},
    "nerve": {
        "label": "nerve", "table_label": "peripheral nerve", "color": "#000000", "class": "removed",
        "canonical": [],
        "core": ["Mpz", "Plp1"],
        "extended": ["Mpz", "Pmp22", "Plp1", "Sox10", "Mbp"],
        "competes": [],
        "excluded": {}},
}

DEFAULT_CONFIG = {
    "config_version": CONFIG_VERSION,
    "name": "A549 intratibial lesion vs contralateral control",

    # ---- data ----------------------------------------------------------------------------------------------------
    # Folder that contains the sample folders. Relative sample paths are resolved against it.
    "data_root": ".",
    # One entry per Visium HD section. `path` may point at the sample folder, its outs/ folder, outs/binned_outputs/
    # or a square_XXXum/ folder (see resolve_bin_dir). microns_per_pixel: null = read `microns_per_pixel` from the
    # section's spatial/scalefactors_json.json.
    "samples": [
        {"id": "LL-2-1", "label": "control", "condition": "control", "color": "#4C72B0",
         "path": "SR007144_LL-2-1/outs/binned_outputs", "microns_per_pixel": None},
        {"id": "RL-1-5", "label": "tumor-bearing", "condition": "tumor", "color": "#C44E52",
         "path": "SR007144_RL-1-5/outs/binned_outputs", "microns_per_pixel": None},
    ],
    # Fold changes and ratios are every other sample over this one.
    "reference_sample": "LL-2-1",
    # Bin size analysed (Space Ranger folder square_<bin_size>), and the bin size of Supplementary Figure S1 (QC).
    "bin_size": "016um",
    "qc_bin_size": "008um",

    # ---- analysis parameters ------------------------------------------------------------------------------------
    "qc": {"min_genes": 10, "min_counts": 15, "max_pct_mt": 20.0, "sweep_thresholds": [0, 5, 10, 25, 50, 100, 200]},
    "clustering": {
        "target_sum": 10000.0, "n_top_genes": 2000, "scale_max_value": 10.0, "n_pcs": 30, "harmony": True,
        "n_neighbors": 15, "resolution": 0.3, "leiden_n_iterations": 2,
        "n_top_markers": 30, "screen_top": 10, "ubiquitous_fraction": 0.3, "min_domain_bins": 100,
        "screen_genes_per_umi": 0.94, "screen_ubiquitous_flag": 0.2,
        "stability_seeds": [0, 1, 2, 3, 4],
    },
    "scoring": {
        "ctrl_size": 50, "null_percentile": 95.0, "n_perm": 50, "n_perm_resampled": 20, "null_n_sub": 20000,
        "min_core": 2, "n_depth_reps": 5, "pool_k": 8,
        # cap (manuscript): bins above the shallowest section's median UMI are downsampled to it, bins below are kept;
        # quantile: every other section is downsampled bin by bin to the shallowest section's depth distribution.
        "depth_matching": "cap",
    },
    "roi": {
        # Transcriptional region of interest: a bin is excluded when its score for any of these signatures is above
        # that signature's `exclude_percentile` over all bins of all sections (a listed signature that has a positivity
        # call, i.e. is in `primary`, excludes its positive bins instead, as in the original code).
        "exclude": ["skeletal_muscle", "nerve", "smooth_muscle"], "exclude_percentile": 90.0,
        # Hand-drawn region of interest (GeoJSON polygons). masks maps sample id -> GeoJSON path; mask_frame is the
        # image the polygons were drawn on (hires = spatial/tissue_hires_image.png).
        "use_drawn": False, "masks": {}, "mask_frame": "hires",
    },
    # Osteoclast-positive bins above the split_by percentile are the growth-plate compartment; the others are the
    # tumor-associated compartment in sections whose condition is listed in tumor_conditions (elsewhere they are
    # "non_growth_plate_control" and not tested). de_sections: all = growth-plate bins of every section against the
    # tumor-associated bins; tumor = both compartments from tumor sections only.
    "compartments": {"signature": "osteoclast", "split_by": "chondrocyte", "percentile": 90.0,
                     "tumor_conditions": ["tumor"], "de_sections": "all"},
    "uncoupling": {"a": "osteoblast", "b": "osteoclast", "calls": "pooled_dm", "n_perm": 199, "min_pos": 10},
    "proximity": {
        "radius_um": 50.0,
        # false (manuscript): a bin that is itself in the target set is its own nearest target (0 um);
        # true: such a bin is measured to the nearest other target bin. Co-positive bins are counted either way.
        "exclude_self": False,
        "pairs": [
            {"from": "fibroblast", "to": "osteoblast", "mode": "assignment"},
            {"from": "macrophage", "to": "osteoclast", "mode": "positivity"},
        ],
    },
    "moran": {"n_genes": 2000, "n_consensus": 200},
    "sensitivity": {
        "signatures": ["osteoclast", "osteoblast"],
        "competitors": ["erythroid", "neutrophil", "fibroblast", "chondrocyte"],
        "variants": {
            "canonical_two_gene": {"osteoclast": ["Atp6v0d2", "Ctsk"], "osteoblast": ["Col1a1", "Sp7"]},
            "original_study_design": {
                "osteoclast": ["Acp5", "Ctsk", "Calcr", "Mmp9", "Atp6v0d2", "Dcstamp", "Ocstamp", "Nfatc1"],
                "osteoblast": ["Runx2", "Sp7", "Alpl", "Col1a1", "Col1a2", "Spp1", "Bglap"]},
            "table1_signatures": {},
        },
    },

    # ---- gene sets --------------------------------------------------------------------------------------------
    "signatures": _SIGNATURES,
    # Tissue classes drawn in Figure 3a, in legend order. emphasis = drawn with larger markers.
    "classes": [
        {"name": "bone lineage", "color": "#984EA3", "emphasis": True},
        {"name": "chondrocyte", "color": "#33A02C", "emphasis": True},
        {"name": "hematopoietic marrow", "color": "#FF7F00", "emphasis": False},
        {"name": "stroma", "color": "#E6AB02", "emphasis": False},
    ],
    # Signatures with a positivity call (Figures 4g, 6a; Methods), in this order.
    "primary": ["osteoclast", "osteoblast", "osteocyte", "macrophage", "neutrophil", "erythroid", "adipocyte",
                "chondrocyte", "fibroblast", "endothelium", "megakaryocyte"],
    # Signatures competing in the deterministic tissue assignment (Figure 3), in this order (ties go to the first).
    "programs": ["osteoclast", "osteoblast", "osteocyte", "chondrocyte", "macrophage", "neutrophil", "erythroid",
                 "megakaryocyte", "adipocyte", "endothelium", "fibroblast", "skeletal_muscle", "smooth_muscle",
                 "nerve"],
    # Host receptor panel (Methods, "Host receptor panel and proximity analysis"; Figure 6b): gene -> role.
    "receptor_panel": {
        "Pth1r": "PTHrP/PTH receptor (host stroma)", "Tnfsf11": "RANKL (host-produced, PTHrP-inducible)",
        "Tnfrsf11b": "OPG (host-produced decoy)", "Tnfrsf11a": "RANK (osteoclast precursor)", "Csf1": "M-CSF (host)",
        "Csf1r": "M-CSF receptor", "Il6ra": "IL-6R", "Il6st": "gp130", "Il11ra1": "IL-11R",
        "Cxcr2": "IL-8/CXCL receptor", "Socs3": "STAT3 target", "Tgfbr1": "TGF-beta R1", "Tgfbr2": "TGF-beta R2",
        "Serpine1": "TGF-beta target",
    },
    # Supplementary Figure S4.
    "adipocyte_genes": ["Adipoq", "Plin1", "Cidec", "Lep", "Retn", "Cfd", "Plin4", "Lipe", "Pparg"],
    # Transcripts abundant in essentially every cell; a domain headed by these has no interpretable identity.
    "ubiquitous": ["Itm2b", "Tmsb10", "Tmsb4x", "App", "Ptms", "Arglu1", "Actb", "Gapdh", "Eef1a1", "Eef2", "Hspa8",
                   "Ftl1", "Fth1", "Ubc", "Ubb", "Npm1", "Hnrnpa1", "Myl6", "Pabpc1", "Ywhae", "Calm1", "Calm2",
                   "Atp5f1e", "Cox4i1", "Cfl1", "Rack1", "Serf2", "Btf3", "Sumo2", "Uba52", "Nedd8", "Tpt1", "Sept7"],

    # ---- what runs ----------------------------------------------------------------------------------------------
    "steps": {
        "qc_figure_bins": True,      # per-bin QC of `qc_bin_size` for Supplementary Figure S1
        "receptor_panel": True,      # Figure 6b
        "adipocyte_panel": True,     # Supplementary Figure S4
        "depth_matched": True,       # Figures 4g, 6a (right panels)
        "leiden": True,              # domain markers, screen, collapse (Supplementary Figure S5)
        "compartments": True,        # Figure 4h
        "assignment": True,          # Figure 3 (spatially pooled scores, deterministic assignment)
        "stability": True,           # adjusted Rand index across seeds (Methods)
        "uncoupling": True,          # Figure 4i
        "proximity": True,           # "within 50 um" tables
        "moran": True,               # Supplementary Figure S2 (slow)
        "sensitivity": False,        # three gene sets (Methods; slow)
        "project_fine_bins": False,  # project labels onto `qc_bin_size` bins
        "write_h5ad": True,          # joint_<bin>.h5ad
        "figures": True,             # bone_figures.make_figures after the analysis
    },
    "figures": {
        "dpi": 200, "map_dpi": 250,
        "positivity_panels": [
            {"name": "Figure4g_bone_lineage_positivity", "signatures": ["osteoclast", "osteoblast", "osteocyte"]},
            {"name": "Figure6a_myeloid_positivity",
             "signatures": ["macrophage", "fibroblast", "neutrophil", "erythroid", "megakaryocyte"]},
        ],
        "receptor_genes": ["Tnfsf11", "Tnfrsf11b", "Tnfrsf11a", "Pth1r", "Csf1", "Csf1r", "Serpine1"],
        "receptor_shade_first_two": True,
        "composition_signatures": ["osteoclast", "osteoblast", "osteocyte", "chondrocyte", "fibroblast", "macrophage",
                                   "neutrophil", "erythroid", "megakaryocyte", "adipocyte"],
        "moran_top_n": 30,
        # Figure 4a-f: the two signatures mapped (blue, red) and their overlay
        "lineage_maps": {"a": "osteoclast", "b": "osteoblast"},
    },
    "seed": 0,
    # Single-threaded numerics and PYTHONHASHSEED=0, so that a rerun reproduces every number (slower).
    "deterministic": True,
}

# Keys whose value replaces the default wholesale instead of being merged key by key (so entries can be deleted).
_REPLACE = {"samples", "signatures", "classes", "primary", "programs", "receptor_panel", "adipocyte_genes",
            "ubiquitous", "masks", "pairs", "variants", "positivity_panels", "receptor_genes",
            "composition_signatures", "exclude", "stability_seeds", "sweep_thresholds", "tumor_conditions"}

# ============================================================================
# PARAMETER SCHEMA (drives the "Parameters" form of the web application)
# ============================================================================
# field types: int, float, bool, select (options), int_list, str_list, signature (one signature name),
#              signature_list (several signature names)
PARAM_SCHEMA = [
    {"key": "qc", "title": "Quality control", "description": "Bins kept for analysis (Methods, \"Quality control and "
     "filtering\"). The two sections differ about two-fold in capture; a single-cell gene threshold removes bins "
     "preferentially from the lower-capture section.", "fields": [
        {"key": "min_genes", "label": "Minimum genes per bin (more than)", "type": "int", "min": 0, "max": 100000,
         "help": "A bin is kept when it detects more than this many genes. Original pipeline: 200."},
        {"key": "min_counts", "label": "Minimum UMIs per bin (more than)", "type": "int", "min": 0, "max": 10000000,
         "help": "A bin is kept when it has more than this many UMIs."},
        {"key": "max_pct_mt", "label": "Maximum % mitochondrial counts (less than)", "type": "float", "min": 0,
         "max": 100, "step": 0.5, "help": "Mitochondrial genes are those whose name starts with mt- (any case). "
         "Original pipeline: 5%."},
        {"key": "sweep_thresholds", "label": "Retention sweep: minimum gene counts", "type": "int_list",
         "help": "Thresholds of the retention curve in Supplementary Figure S1."},
    ]},
    {"key": "clustering", "title": "Normalization, embedding and Leiden domains", "description": "Used for the "
     "Leiden domains of Supplementary Figure S5 only; scores and tests use the full log-normalized matrix.",
     "fields": [
        {"key": "target_sum", "label": "Normalization target (counts per bin)", "type": "float", "min": 1,
         "help": "normalize_total target before log1p."},
        {"key": "n_top_genes", "label": "Highly variable genes", "type": "int", "min": 50, "max": 20000,
         "help": "Seurat-flavor HVGs computed per section (batch_key = sample)."},
        {"key": "scale_max_value", "label": "Scaling clip value", "type": "float", "min": 1, "max": 1000},
        {"key": "n_pcs", "label": "Principal components used", "type": "int", "min": 2, "max": 50},
        {"key": "harmony", "label": "Harmony integration on the section label", "type": "bool",
         "help": "Skipped automatically when only one section is analysed."},
        {"key": "n_neighbors", "label": "Neighbors in the kNN graph", "type": "int", "min": 2, "max": 200},
        {"key": "resolution", "label": "Leiden resolution (modularity)", "type": "float", "min": 0.01, "max": 10,
         "step": 0.05},
        {"key": "leiden_n_iterations", "label": "Leiden iterations", "type": "int", "min": -1, "max": 100,
         "help": "-1 iterates until stable."},
        {"key": "n_top_markers", "label": "Markers exported per domain", "type": "int", "min": 1, "max": 1000},
        {"key": "screen_top", "label": "Top markers inspected by the domain screen", "type": "int", "min": 1,
         "max": 200},
        {"key": "ubiquitous_fraction", "label": "Merge domains whose top markers are at least this fraction "
         "ubiquitous", "type": "float", "min": 0, "max": 1, "step": 0.05},
        {"key": "min_domain_bins", "label": "Domains smaller than this become 'unassigned'", "type": "int",
         "min": 0},
        {"key": "screen_genes_per_umi", "label": "Shallow-sampling flag: median genes per UMI at least",
         "type": "float", "min": 0, "max": 1, "step": 0.01},
        {"key": "screen_ubiquitous_flag", "label": "Shallow-sampling flag: ubiquitous fraction at least",
         "type": "float", "min": 0, "max": 1, "step": 0.05},
        {"key": "stability_seeds", "label": "Seeds for the clustering-stability check", "type": "int_list"},
    ]},
    {"key": "scoring", "title": "Gene-set scoring and positivity", "description": "A bin is positive for a "
     "signature when its score exceeds the null threshold, enough core genes are detected and no competing "
     "signature scores higher (Methods; Table 1).", "fields": [
        {"key": "ctrl_size", "label": "Control genes per scored gene (score_genes ctrl_size)", "type": "int",
         "min": 1, "max": 1000},
        {"key": "null_percentile", "label": "Null threshold percentile", "type": "float", "min": 50, "max": 100,
         "step": 0.5},
        {"key": "n_perm", "label": "Random gene sets in the null (native depth)", "type": "int", "min": 1,
         "max": 10000},
        {"key": "n_perm_resampled", "label": "Random gene sets in the null (depth-matched and pooled)",
         "type": "int", "min": 1, "max": 10000},
        {"key": "null_n_sub", "label": "Bins sampled for the null", "type": "int", "min": 100},
        {"key": "min_core", "label": "Core genes required for a positive call", "type": "int", "min": 0, "max": 50},
        {"key": "n_depth_reps", "label": "Depth-matching replicates", "type": "int", "min": 1, "max": 100,
         "help": "Raw counts downsampled to the median UMI per bin of the shallowest section."},
        {"key": "pool_k", "label": "Spatial neighbors pooled per bin (Figure 3)", "type": "int", "min": 0,
         "max": 48, "help": "Counts of each bin and its k nearest bins are summed before scoring."},
        {"key": "depth_matching", "label": "Depth matching", "type": "select", "options": ["cap", "quantile"],
         "help": "cap (manuscript): bins above the shallowest section's median UMI per bin are downsampled to it, "
                 "bins below are kept, so deeper sections keep more bins at that depth. quantile: each bin of the "
                 "other sections is downsampled to the shallowest section's depth at the same quantile, so the "
                 "per-bin depth distributions match."},
    ]},
    {"key": "roi", "title": "Region of interest", "description": "Marrow and bone region of interest. The "
     "transcriptional part excludes bins scoring high for the listed signatures; a hand-drawn region (ROI tab) "
     "can be intersected with it.", "fields": [
        {"key": "exclude", "label": "Signatures that exclude a bin", "type": "signature_list",
         "help": "A bin scoring above the exclusion percentile for any of these is outside the region. A signature with "
                 "a positivity call excludes its positive bins instead."},
        {"key": "exclude_percentile", "label": "Exclusion percentile", "type": "float", "min": 0, "max": 100,
         "step": 0.5},
        {"key": "use_drawn", "label": "Intersect with the hand-drawn region", "type": "bool",
         "help": "Sections without a drawn region are not restricted further."},
        {"key": "mask_frame", "label": "Image the regions were drawn on", "type": "select",
         "options": ["hires", "lowres", "fullres"]},
    ]},
    {"key": "compartments", "title": "Osteoclast compartments (Figure 4h)", "fields": [
        {"key": "signature", "label": "Signature split into compartments", "type": "signature"},
        {"key": "split_by", "label": "Split by the score of", "type": "signature"},
        {"key": "percentile", "label": "Growth-plate percentile", "type": "float", "min": 0, "max": 100,
         "step": 0.5},
        {"key": "de_sections", "label": "Sections in the compartment test", "type": "select",
         "options": ["all", "tumor"], "help": "all: growth-plate bins of every section against tumor-associated "
         "bins; tumor: both compartments from tumor sections only (no between-section depth difference)."},
    ]},
    {"key": "uncoupling", "title": "Distance to nearest partner (Figure 4i)", "fields": [
        {"key": "a", "label": "From (signature)", "type": "signature"},
        {"key": "b", "label": "To (signature)", "type": "signature"},
        {"key": "calls", "label": "Positivity calls used", "type": "select",
         "options": ["pooled_dm", "pooled", "native"],
         "help": "pooled_dm: spatially pooled, depth-matched (manuscript)."},
        {"key": "n_perm", "label": "Random-placement replicates", "type": "int", "min": 1, "max": 100000},
        {"key": "min_pos", "label": "Minimum positive bins per side", "type": "int", "min": 1},
    ]},
    {"key": "proximity", "title": "Proximity tables", "fields": [
        {"key": "radius_um", "label": "Radius (um, center to center)", "type": "float", "min": 0, "step": 1},
        {"key": "exclude_self", "label": "Measure target bins to the nearest other target bin", "type": "bool",
         "help": "Off (manuscript): a bin positive for both signatures counts as within the radius (distance 0). "
                 "On: it is measured to the nearest other target bin. Co-positive bins are reported either way."},
    ]},
    {"key": "moran", "title": "Spatially variable genes (Supplementary Figure S2)", "fields": [
        {"key": "n_genes", "label": "Highly variable genes tested", "type": "int", "min": 10, "max": 50000},
        {"key": "n_consensus", "label": "Top genes per section for the consensus", "type": "int", "min": 1},
    ]},
]

STEP_LABELS = {
    "qc_figure_bins": "Per-bin QC for Supplementary Figure S1",
    "receptor_panel": "Host receptor panel (Figure 6b)",
    "adipocyte_panel": "Adipocyte genes (Supplementary Figure S4)",
    "depth_matched": "Depth-matched positivity (Figures 4g, 6a)",
    "leiden": "Leiden domain markers and screen (Supplementary Figure S5)",
    "compartments": "Osteoclast compartments (Figure 4h)",
    "assignment": "Deterministic tissue assignment (Figure 3)",
    "stability": "Leiden seed stability (Methods)",
    "uncoupling": "Osteoblast-to-osteoclast distance (Figure 4i)",
    "proximity": "Proximity tables (within 50 um)",
    "moran": "Moran's I (Supplementary Figure S2, slow)",
    "sensitivity": "Three-gene-set sensitivity (Methods, slow)",
    "project_fine_bins": "Project labels onto the QC bin size",
    "write_h5ad": "Write the joint AnnData object",
    "figures": "Draw figures",
}

# ============================================================================
# HELPERS
# ============================================================================
# \Z, not $: "$" also matches before a trailing newline, which let "abc\n" pass as an id
_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]{0,63}\Z")
_SIG_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,63}\Z")
_BIN_RE = re.compile(r"^(\d{3})um\Z")
_HEX_RE = re.compile(r"^#[0-9A-Fa-f]{6}\Z")


def default_config() -> dict:
    """A deep copy of DEFAULT_CONFIG."""
    return copy.deepcopy(DEFAULT_CONFIG)


def _merge(base, over, key=None):
    if isinstance(base, dict) and isinstance(over, dict) and key not in _REPLACE:
        out = dict(base)
        for k, v in over.items():
            out[k] = _merge(base.get(k), v, k) if k in base else copy.deepcopy(v)
        return out
    return copy.deepcopy(over)


def merge_with_defaults(cfg: dict | None) -> dict:
    """Deep-merge `cfg` over the defaults. Dicts merge key by key, except the _REPLACE keys (samples, signatures,
    gene lists, ...), which replace the default wholesale so that entries can be removed."""
    return _merge(default_config(), cfg or {})


def bin_um(bin_size: str) -> float:
    """'016um' -> 16.0"""
    m = _BIN_RE.match(str(bin_size))
    if not m:
        raise ValueError(f"bin size must look like '016um', got {bin_size!r}")
    return float(int(m.group(1)))


def slugify(text: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "-", str(text)).strip("-")
    return s or "x"


def sample_by_id(cfg: dict, sample_id: str) -> dict:
    for s in cfg["samples"]:
        if s["id"] == sample_id:
            return s
    raise KeyError(f"unknown sample {sample_id!r}")


def sample_ids(cfg: dict) -> list[str]:
    return [s["id"] for s in cfg["samples"]]


def signature_label(cfg: dict, name: str) -> str:
    return cfg["signatures"].get(name, {}).get("label") or name.replace("_", " ")


def _sample_base(cfg: dict, sample: dict) -> Path:
    p = Path(str(sample.get("path") or sample["id"])).expanduser()
    if not p.is_absolute():
        p = Path(str(cfg.get("data_root") or ".")).expanduser() / p
    return p


def _bin_candidates(base: Path, bin_size: str) -> list[Path]:
    name = f"square_{bin_size}"
    c = [base / name, base / "binned_outputs" / name, base / "outs" / "binned_outputs" / name]
    if base.name.startswith("square_"):
        c.insert(0, base.parent / name)
    return c


def resolve_bin_dir(cfg: dict, sample: dict | str, bin_size: str | None = None) -> Path | None:
    """The square_<bin_size>/ folder of a sample (the one that holds filtered_feature_bc_matrix.h5), or None.
    `path` may point at the sample folder, outs/, outs/binned_outputs/ or any square_XXXum/ folder."""
    if isinstance(sample, str):
        sample = sample_by_id(cfg, sample)
    bin_size = bin_size or cfg["bin_size"]
    for d in _bin_candidates(_sample_base(cfg, sample), bin_size):
        if (d / "filtered_feature_bc_matrix.h5").is_file():
            return d
    return None


def available_bin_sizes(cfg: dict, sample: dict | str) -> list[str]:
    """Bin sizes (e.g. ['002um', '008um', '016um']) whose square_XXXum/filtered_feature_bc_matrix.h5 exists."""
    if isinstance(sample, str):
        sample = sample_by_id(cfg, sample)
    base = _sample_base(cfg, sample)
    roots = [base, base / "binned_outputs", base / "outs" / "binned_outputs"]
    if base.name.startswith("square_"):
        roots.insert(0, base.parent)
    found = set()
    for r in roots:
        if r.is_dir():
            for d in r.glob("square_*um"):
                m = re.match(r"^square_(\d{3}um)$", d.name)
                if m and (d / "filtered_feature_bc_matrix.h5").is_file():
                    found.add(m.group(1))
    return sorted(found)


def validate_config(cfg: dict) -> tuple[list[str], list[str]]:
    """Structural checks that need no data. Returns (errors, warnings); a configuration with errors must not run."""
    err, warn = [], []

    def num(sec, key, lo=None, hi=None, integer=False):
        try:
            v = cfg[sec][key]
        except (KeyError, TypeError):
            err.append(f"{sec}.{key} is missing")
            return
        if isinstance(v, bool) or not isinstance(v, (int, float)) or (isinstance(v, float) and not math.isfinite(v)):
            err.append(f"{sec}.{key} must be a number, got {v!r}")
            return
        if integer and int(v) != v:
            err.append(f"{sec}.{key} must be a whole number, got {v!r}")
        if lo is not None and v < lo:
            err.append(f"{sec}.{key} must be >= {lo}, got {v}")
        if hi is not None and v > hi:
            err.append(f"{sec}.{key} must be <= {hi}, got {v}")

    # samples
    samples = cfg.get("samples") or []
    if not samples:
        err.append("at least one sample is required")
    ids = [s.get("id") for s in samples]
    for s in samples:
        sid = s.get("id")
        if not isinstance(sid, str) or not _ID_RE.match(sid):
            err.append(f"sample id {sid!r} must be 1-64 characters of letters, digits, '.', '_' or '-'")
        if not s.get("label"):
            warn.append(f"sample {sid!r} has no label; the id is used")
        c = s.get("color")
        if c and not _HEX_RE.match(str(c)):
            err.append(f"sample {sid!r}: color must be #RRGGBB, got {c!r}")
        mpp = s.get("microns_per_pixel")
        if mpp is not None and (isinstance(mpp, bool) or not isinstance(mpp, (int, float)) or mpp <= 0):
            err.append(f"sample {sid!r}: microns_per_pixel must be a positive number or null")
    if len(set(ids)) != len(ids):
        err.append("sample ids must be unique")
    elif len({i.casefold() for i in ids if isinstance(i, str)}) != len(ids):
        err.append("sample ids must differ by more than letter case (their region files would collide on Windows "
                   "and macOS)")
    if samples and cfg.get("reference_sample") not in ids:
        err.append(f"reference_sample {cfg.get('reference_sample')!r} is not one of the samples {ids}")
    for key in ("bin_size", "qc_bin_size"):
        if not _BIN_RE.match(str(cfg.get(key, ""))):
            err.append(f"{key} must look like '016um', got {cfg.get(key)!r}")
    if (cfg.get("steps") or {}).get("project_fine_bins") and all(_BIN_RE.match(str(cfg.get(k, "")))
                                                                for k in ("bin_size", "qc_bin_size")):
        fine, coarse = bin_um(cfg["qc_bin_size"]), bin_um(cfg["bin_size"])
        if fine <= 0 or fine > coarse or coarse % fine:
            err.append(f"steps.project_fine_bins projects {cfg['bin_size']} labels onto {cfg['qc_bin_size']} bins: "
                       f"qc_bin_size must be finer than bin_size and divide it")
    seeds = (cfg.get("clustering") or {}).get("stability_seeds", [])
    if not isinstance(seeds, list) or not all(isinstance(s, int) and not isinstance(s, bool) and s >= 0 for s in seeds):
        err.append("clustering.stability_seeds must be a list of non-negative whole numbers")
    elif len(set(seeds)) != len(seeds):
        err.append("clustering.stability_seeds lists a seed twice")

    # numbers
    num("qc", "min_genes", 0, integer=True)
    num("qc", "min_counts", 0)
    num("qc", "max_pct_mt", 0, 100)
    num("clustering", "target_sum", 1)
    num("clustering", "n_top_genes", 2, integer=True)
    num("clustering", "n_pcs", 2, integer=True)
    num("clustering", "n_neighbors", 2, integer=True)
    num("clustering", "resolution", 0)
    num("clustering", "ubiquitous_fraction", 0, 1)
    num("clustering", "min_domain_bins", 0, integer=True)
    num("scoring", "ctrl_size", 1, integer=True)
    num("scoring", "null_percentile", 0, 100)
    num("scoring", "n_perm", 1, integer=True)
    num("scoring", "n_perm_resampled", 1, integer=True)
    num("scoring", "null_n_sub", 2, integer=True)
    num("scoring", "min_core", 0, integer=True)
    num("scoring", "n_depth_reps", 1, integer=True)
    num("scoring", "pool_k", 0, integer=True)
    num("roi", "exclude_percentile", 0, 100)
    num("compartments", "percentile", 0, 100)
    num("uncoupling", "n_perm", 1, integer=True)
    num("uncoupling", "min_pos", 1, integer=True)
    num("proximity", "radius_um", 0)
    num("moran", "n_genes", 2, integer=True)
    num("moran", "n_consensus", 1, integer=True)
    try:
        if int(cfg["clustering"]["n_pcs"]) > 50:
            warn.append("clustering.n_pcs > 50: PCA computes 50 components; it will be capped")
    except (KeyError, TypeError, ValueError):
        pass
    if cfg.get("roi", {}).get("mask_frame") not in ("hires", "lowres", "fullres"):
        err.append("roi.mask_frame must be hires, lowres or fullres")
    if cfg.get("uncoupling", {}).get("calls") not in ("pooled_dm", "pooled", "native"):
        err.append("uncoupling.calls must be pooled_dm, pooled or native")
    if cfg.get("scoring", {}).get("depth_matching") not in ("cap", "quantile"):
        err.append("scoring.depth_matching must be cap or quantile")
    if cfg.get("compartments", {}).get("de_sections") not in ("all", "tumor"):
        err.append("compartments.de_sections must be all or tumor")
    if not isinstance(cfg.get("compartments", {}).get("tumor_conditions", []), list):
        err.append("compartments.tumor_conditions must be a list of condition names")
    if not isinstance(cfg.get("proximity", {}).get("exclude_self", False), bool):
        err.append("proximity.exclude_self must be true or false")

    # signatures
    sigs = cfg.get("signatures") or {}
    if not sigs:
        err.append("at least one signature is required")
    for name, d in sigs.items():
        if not _SIG_RE.match(str(name)):
            err.append(f"signature name {name!r} must start with a letter and use letters, digits or '_'")
            continue
        for key in ("core", "extended"):
            if not isinstance(d.get(key), list) or not all(isinstance(g, str) and g for g in d.get(key, [])):
                err.append(f"signature {name}: {key} must be a list of gene names")
        if isinstance(d.get("extended"), list) and len(d["extended"]) < 2:
            warn.append(f"signature {name}: fewer than two scored genes; it cannot be scored")
        missing_core = [g for g in d.get("core", []) if g not in d.get("extended", [])]
        if missing_core:
            warn.append(f"signature {name}: core genes not in the scored list: {', '.join(missing_core)}")
        for c in d.get("competes", []):
            if c not in sigs:
                err.append(f"signature {name} competes with unknown signature {c!r}")
            if c == name:
                err.append(f"signature {name} cannot compete with itself")
        col = d.get("color")
        if col and not _HEX_RE.match(str(col)):
            err.append(f"signature {name}: color must be #RRGGBB, got {col!r}")
        if not isinstance(d.get("excluded", {}), dict):
            err.append(f"signature {name}: excluded must map gene -> reason")
    for key in ("primary", "programs"):
        lst = cfg.get(key) or []
        for n in lst:
            if n not in sigs:
                err.append(f"{key} lists unknown signature {n!r}")
        if len(set(lst)) != len(lst):
            err.append(f"{key} lists a signature twice")
    if not cfg.get("primary"):
        err.append("primary (signatures with a positivity call) must not be empty")
    for n in cfg.get("roi", {}).get("exclude", []):
        if n not in sigs:
            err.append(f"roi.exclude lists unknown signature {n!r}")
        elif n in (cfg.get("primary") or []):
            warn.append(f"roi.exclude: {n} has a positivity call, so its positive bins are excluded and "
                        f"roi.exclude_percentile does not apply to it")
    for sec, keys in (("compartments", ("signature", "split_by")), ("uncoupling", ("a", "b"))):
        for k in keys:
            n = cfg.get(sec, {}).get(k)
            if n not in sigs:
                warn.append(f"{sec}.{k} = {n!r} is not a signature; that step will be skipped")
            elif n not in cfg.get("primary", []) and not (sec == "compartments" and k == "split_by"):
                warn.append(f"{sec}.{k} = {n!r} has no positivity call (not in primary); that step will be skipped")
    for p in cfg.get("proximity", {}).get("pairs", []):
        if p.get("mode") not in ("assignment", "positivity"):
            err.append(f"proximity pair {p}: mode must be assignment or positivity")
        for k in ("from", "to"):
            if p.get(k) not in sigs:
                warn.append(f"proximity pair {p.get('from')} -> {p.get('to')}: {p.get(k)!r} is not a signature; "
                            f"the pair will be skipped")
    for c in cfg.get("classes", []):
        if not c.get("name"):
            err.append("every tissue class needs a name")
        if c.get("color") and not _HEX_RE.match(str(c["color"])):
            err.append(f"class {c.get('name')!r}: color must be #RRGGBB")
    if not isinstance(cfg.get("receptor_panel"), dict):
        err.append("receptor_panel must map gene -> role")
    for k, v in (cfg.get("steps") or {}).items():
        if not isinstance(v, bool):
            err.append(f"steps.{k} must be true or false")
    return err, warn


def effective_config(cfg: dict | None) -> dict:
    """Defaults merged with `cfg`, then checked. Raises ValueError listing every error."""
    full = merge_with_defaults(cfg)
    errors, _ = validate_config(full)
    if errors:
        raise ValueError("invalid configuration:\n  - " + "\n  - ".join(errors))
    return full


def load_config(path) -> dict:
    """Read a JSON configuration and merge it over the defaults (not validated)."""
    return merge_with_defaults(json.loads(Path(path).read_text(encoding="utf-8")))


def save_config(cfg: dict, path) -> None:
    Path(path).write_text(json.dumps(cfg, indent=2), encoding="utf-8")


# ============================================================================
# PROGRESS REPORTING
# ============================================================================
EVENT_PREFIX = "@@EVENT "


class Progress:
    """
    [APP] Step reporting shared by bone_analysis.run_analysis_from_config, bone_figures.make_figures and the web
    worker. Every call prints a human-readable line. With events=True it also prints one machine-readable line
    '@@EVENT {json}' that the web application's job runner parses. Event types:
      plan {steps: [label, ...]}            the steps expected, in order
      step {label, index, total}            a step starts (closes the previous one)
      step_done {label, seconds}            a step finished
      skip {label, reason}                  a planned step will not run
      warning {message}
      done {}                               everything finished
      error {message, traceback}
    """

    def __init__(self, events: bool = False, stream=None):
        self.events = events
        self.stream = stream
        self.labels: list[str] = []
        self.current: str | None = None
        self._t0 = 0.0
        self.warnings: list[str] = []
        self.timings: list[dict] = []
        self.skipped: list[dict] = []

    def _out(self):
        return self.stream or sys.stdout

    def _emit(self, payload: dict) -> None:
        if self.events:
            print(EVENT_PREFIX + json.dumps(payload, default=str), file=self._out(), flush=True)

    def plan(self, labels) -> None:
        self.labels = list(dict.fromkeys(labels))
        self._emit({"type": "plan", "steps": self.labels})

    def _close(self) -> None:
        if self.current is not None:
            secs = round(time.perf_counter() - self._t0, 2)
            self.timings.append({"step": self.current, "seconds": secs})
            self._emit({"type": "step_done", "label": self.current, "seconds": secs})
            print(f"    ({secs:.1f} s)", file=self._out(), flush=True)
            self.current = None

    def step(self, label: str) -> None:
        self._close()
        if label not in self.labels:
            self.labels.append(label)
        idx = self.labels.index(label)
        self.current, self._t0 = label, time.perf_counter()
        print(f"\n=== [{idx + 1}/{len(self.labels)}] {label}", file=self._out(), flush=True)
        self._emit({"type": "step", "label": label, "index": idx, "total": len(self.labels)})

    def skip(self, label: str, reason: str) -> None:
        self._close()
        if label not in self.labels:
            self.labels.append(label)
        self.skipped.append({"step": label, "reason": reason})
        print(f"--- skipped: {label} ({reason})", file=self._out(), flush=True)
        self._emit({"type": "skip", "label": label, "reason": reason})

    def info(self, message: str) -> None:
        print(message, file=self._out(), flush=True)

    def warn(self, message: str) -> None:
        self.warnings.append(message)
        print(f"WARNING: {message}", file=self._out(), flush=True)
        self._emit({"type": "warning", "message": message})

    def done(self) -> None:
        self._close()
        self._emit({"type": "done"})

    def error(self, message: str, traceback_text: str = "") -> None:
        self.current = None
        self._emit({"type": "error", "message": message, "traceback": traceback_text})


if __name__ == "__main__":
    # python pipeline_config.py > config.json   writes the default (manuscript) configuration
    print(json.dumps(default_config(), indent=2))
