#!/usr/bin/env python3
"""
bone_analysis.py
================
Visium HD analysis pipeline: A549 intratibial lesion (right tibia, RL-1-5) against the
contralateral control (left tibia, LL-2-1) of the same NSG mouse, both sections on one capture slide, one CytAssist run,
Visium Mouse Transcriptome Probe Set v2.0 against mm10-2020-A. Human A549 transcripts are not captured: every transcript
is murine.

Reads the Space Ranger binned outputs and writes every table and number the manuscript reports from the transcriptome,
plus the joint AnnData object the figure script draws from. 

PROVENANCE TAGS (every function and every step of run_analysis carries one)
----------------------------------------------------------------------------
  [SHAMIM]            taken unchanged from Shamim Mollah's original pipeline,
                      00_shamim_original/visium_spatial_processing_full_pipeline.py (2024-12-15); "step N" is the step
                      number in that file
  [SHAMIM, MODIFIED]  originates in that pipeline, but its parameters or implementation were changed by Run Zhang;
                      the original setting and the reason for the change are written next to it
  [RUN]               added by Run Zhang; not in the original pipeline

Steps of the original pipeline that are NOT used in the manuscript and are therefore not carried over:
  step 2 violin plot, step 6 PCA variance plot, step 8 UMAP plot, step 9 spatial cluster plot, step 10 marker plot,
  step 12 Enrichr GO/KEGG/Reactome on cluster 2, step 14 image features, step 15 CellPhoneDB export.

Usage (from the project root, the folder that contains SR007144_LL-2-1/ and SR007144_RL-1-5/)
-----
    python Manuscript/code/bone_analysis.py --root . --out Manuscript/code/results
    # optional: --geneset-sensitivity (Methods, three gene sets; slow), --skip-moran (Supplementary Figure S2; slow),
    #           --roi-mask SAMPLE=PATH.geojson (hand-drawn ROI, once it exists; repeat per sample)

Environment: ~/venvs/bone_spatial (scanpy 1.11.5, anndata 0.12.9, harmonypy 0.0.10, squidpy 1.8.1, igraph 1.0.0,
leidenalg 0.12.0, numpy 2.3.5, pandas 2.3.3, scipy 1.16.3, statsmodels 0.14.6, scikit-learn 1.8.0).
"""
from __future__ import annotations

import os

# [RUN] Reproducibility. pynndescent builds its trees in parallel and Harmony seeds its k-means internally, so
# random_state alone does not pin the result. 
if os.environ.get("VISIUM_DETERMINISTIC", "1") == "1":
    for _v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(_v, "1")
    os.environ.setdefault("PYTHONHASHSEED", "0")

import argparse
import copy
import json
import re
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

warnings.filterwarnings("ignore")
sc.settings.verbosity = 2
RANDOM_SEED = 0

# ============================================================================
# CONSTANTS  [RUN]
# ============================================================================
SAMPLES = {
    "LL-2-1": {"condition": "control", "limb": "left", "label": "control"},
    "RL-1-5": {"condition": "tumor", "limb": "right", "label": "tumor-bearing"},
}
SECTION = {"LL-2-1": "control", "RL-1-5": "tumor-bearing"}
# microns per full-resolution pixel, `microns_per_pixel` in each section's scalefactors_json.json
UM_PER_PX = {"LL-2-1": 0.5936042159338343, "RL-1-5": 0.5554917593191749}

# QC thresholds (Methods, "Quality control and filtering"). See qc_filter().
MIN_GENES, MIN_COUNTS, MAX_PCT_MT = 10, 15, 20.0

# Gene signatures (Table 1).
#   core     : lineage-restricted genes; at least two must be detected for a positive call
#   extended : genes scored
#   competes : a bin that scores higher on one of these is not called positive
#   excluded : genes deliberately left out, with the reason (Table 1)
SIG = {
    "osteoclast": {"core": ["Acp5", "Ctsk", "Calcr", "Atp6v0d2", "Dcstamp", "Ocstamp"],
                   "extended": ["Acp5", "Ctsk", "Calcr", "Mmp9", "Atp6v0d2", "Dcstamp", "Ocstamp", "Nfatc1"],
                   "competes": ["erythroid", "neutrophil"],
                   "excluded": {"Car2": "abundant in erythroid precursors"}},
    "osteoblast": {"core": ["Runx2", "Sp7", "Alpl", "Ibsp"],
                   "extended": ["Runx2", "Sp7", "Alpl", "Ibsp"],
                   "competes": ["fibroblast", "chondrocyte"],
                   "excluded": {"Col1a1": "pan-fibroblast", "Col1a2": "pan-fibroblast",
                                "Spp1": "major macrophage transcript",
                                "Bglap": "absent from the probe set", "Bglap2": "absent from the probe set"}},
    "osteocyte": {"core": ["Sost", "Dmp1", "Mepe"], "extended": ["Sost", "Dmp1", "Mepe", "Fgf23", "Pdpn"],
                  "competes": [], "excluded": {"Phex": "not detected in this probe set"}},
    "fibroblast": {"core": ["Pdgfra", "Dcn", "Lum"],
                   "extended": ["Pdgfra", "Pdgfrb", "Fap", "Thy1", "Acta2", "Dcn", "Lum", "Col3a1", "Postn", "Fbn1"],
                   "competes": [], "excluded": {}},
    "macrophage": {"core": ["Adgre1", "Csf1r", "C1qa", "C1qb", "Mrc1"],
                   "extended": ["Adgre1", "Csf1r", "C1qa", "C1qb", "C1qc", "Mrc1", "Cd68", "Trem2", "Ctss", "Lyz2", "Apoe"],
                   "competes": ["neutrophil"], "excluded": {}},
    # "neutrophil" is the granulocyte signature of the manuscript (marrow granulopoiesis)
    "neutrophil": {"core": ["S100a8", "S100a9", "Retnlg", "Ngp", "Ltf"],
                   "extended": ["S100a8", "S100a9", "Retnlg", "Ly6g", "Mpo", "Elane", "Ngp", "Ltf", "Camp", "Lcn2"],
                   "competes": [], "excluded": {}},
    "erythroid": {"core": ["Hba-a1", "Hbb-bs", "Alas2"],
                  "extended": ["Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2", "Slc4a1", "Gypa", "Klf1", "Blvrb", "Dmtn"],
                  "competes": [], "excluded": {}},
    "megakaryocyte": {"core": ["Pf4", "Itga2b", "Vwf"],
                      "extended": ["Pf4", "Itga2b", "Vwf", "Ppbp", "Gp9", "Gp1bb", "Thbs1"], "competes": [], "excluded": {}},
    "adipocyte": {"core": ["Adipoq", "Plin1", "Cidec"],
                  "extended": ["Adipoq", "Plin1", "Cidec", "Lep", "Retn", "Cfd", "Plin4", "Lipe", "Pparg"],
                  "competes": ["skeletal_muscle"],
                  "excluded": {"Car3": "skeletal muscle marker", "Fabp4": "endothelium and macrophages"}},
    # Hypoxia was removed on 2026-10-27: not a cell type, and not part of the question. Endothelium stays as a competitor
    # in the tissue assignment, so that vessel-rich bins are not assigned to fibroblast or osteoblast. Like skeletal
    # muscle, nerve and smooth muscle it is not drawn in Figure 3.
    "endothelium": {"core": ["Pecam1", "Cdh5", "Emcn"], "extended": ["Pecam1", "Cdh5", "Emcn", "Flt1", "Kdr", "Cldn5"],
                    "competes": [], "excluded": {}},
    "chondrocyte": {"core": ["Col2a1", "Acan", "Col10a1"], "extended": ["Col2a1", "Acan", "Col10a1", "Sox9", "Comp", "Matn3"],
                    "competes": [], "excluded": {"Mmp13": "also induced in remodeling tissue; circular when used to "
                                                           "define the growth-plate compartment"}},
    "skeletal_muscle": {"core": ["Ttn", "Ckm", "Acta1"],
                        "extended": ["Ttn", "Tnnt3", "Acta1", "Ckm", "Tpm2", "Mylpf", "Pvalb", "Myh1", "Tnnc2", "Atp2a1", "Car3"],
                        "competes": [], "excluded": {}},
    "smooth_muscle": {"core": ["Myh11", "Cnn1"], "extended": ["Myh11", "Tagln", "Cnn1", "Myl9", "Acta2"],
                      "competes": [], "excluded": {}},
    "nerve": {"core": ["Mpz", "Plp1"], "extended": ["Mpz", "Pmp22", "Plp1", "Sox10", "Mbp"], "competes": [], "excluded": {}},
}
# signatures with a positivity call (Figures 4g, 6a; Methods)
PRIMARY = ["osteoclast", "osteoblast", "osteocyte", "macrophage", "neutrophil", "erythroid", "adipocyte", "chondrocyte",
           "fibroblast", "endothelium", "megakaryocyte"]
# signatures competing in the deterministic tissue assignment (Figure 3)
PROGRAMS = ["osteoclast", "osteoblast", "osteocyte", "chondrocyte", "macrophage", "neutrophil", "erythroid",
            "megakaryocyte", "adipocyte", "endothelium", "fibroblast", "skeletal_muscle", "smooth_muscle", "nerve"]

# Host receptor panel (Methods, "Host receptor panel and proximity analysis"; Figure 6b)
RECEPTOR_PANEL = {
    "Pth1r": "PTHrP/PTH receptor (host stroma)", "Tnfsf11": "RANKL (host-produced, PTHrP-inducible)",
    "Tnfrsf11b": "OPG (host-produced decoy)", "Tnfrsf11a": "RANK (osteoclast precursor)", "Csf1": "M-CSF (host)",
    "Csf1r": "M-CSF receptor", "Il6ra": "IL-6R", "Il6st": "gp130", "Il11ra1": "IL-11R", "Cxcr2": "IL-8/CXCL receptor",
    "Socs3": "STAT3 target", "Tgfbr1": "TGF-beta R1", "Tgfbr2": "TGF-beta R2", "Serpine1": "TGF-beta target",
}
ADIPOCYTE_GENES = ["Adipoq", "Plin1", "Cidec", "Lep", "Retn", "Cfd", "Plin4", "Lipe", "Pparg"]   # Supplementary Figure S4

# Transcripts abundant in essentially every cell; a domain headed by these has no interpretable identity (Methods)
UBIQUITOUS = {
    "Itm2b", "Tmsb10", "Tmsb4x", "App", "Ptms", "Arglu1", "Actb", "Gapdh", "Eef1a1", "Eef2", "Hspa8", "Ftl1", "Fth1",
    "Ubc", "Ubb", "Npm1", "Hnrnpa1", "Myl6", "Pabpc1", "Ywhae", "Calm1", "Calm2", "Atp5f1e", "Cox4i1", "Cfl1", "Rack1",
    "Serf2", "Btf3", "Sumo2", "Uba52", "Nedd8", "Tpt1", "Sept7",
}

PALETTE = {
    "control": "#4C72B0", "tumor": "#C44E52",
    "osteoclast": "#377EB8", "osteoblast": "#E41A1C", "osteocyte": "#984EA3", "macrophage": "#DA70D6",
    "neutrophil": "#FF7F00", "adipocyte": "#FFD92F", "erythroid": "#A65628", "chondrocyte": "#33A02C",
    "fibroblast": "#FB9A99", "skeletal_muscle": "#8C8C8C", "nerve": "#000000", "smooth_muscle": "#D2B48C",
    "endothelium": "#17BECF", "megakaryocyte": "#7FC97F",
}


# ============================================================================
# 1. INPUT AND QUALITY CONTROL
# ============================================================================

def _norm_bc(x):
    """[RUN] Barcode string without the bin suffix or the -1 suffix, so the matrix and the position table match."""
    if pd.isna(x):
        return x
    s = re.sub(r"_[0-9]+_[0-9]+(?:_[0-9]+)*$", "", str(x).strip())
    return re.sub(r"-\d+$", "", s)


def read_visium_hd(path, sample_name):
    """
    [SHAMIM, MODIFIED] step 1, "Load Visium spatial data".
    Original: sc.read_visium(path) and adata.var_names_make_unique().
    Changed because sc.read_visium reads standard Visium spots and does not read the Visium HD binned output
    (square_XXXum/, positions in tissue_positions.parquet). The matrix is read with sc.read_10x_h5 and joined to the
    parquet positions by barcode; bins without valid coordinates are dropped (Methods, "Data import and preprocessing").
    var_names_make_unique() is kept as in the original.
    """
    path, spatial = Path(path), Path(path) / "spatial"
    a = sc.read_10x_h5(path / "filtered_feature_bc_matrix.h5")
    a.var_names_make_unique()                                                   # [SHAMIM] step 1
    a.obs_names = [_norm_bc(x) for x in a.obs_names.astype(str)]
    a.obs_names_make_unique()
    tp = pd.read_parquet(spatial / "tissue_positions.parquet")
    bc = "barcode" if "barcode" in tp.columns else tp.columns[0]
    tp[bc] = tp[bc].astype(str).map(_norm_bc)
    tp = tp[tp[bc].isin(a.obs_names)].set_index(bc).reindex(a.obs_names)
    coords = ("pxl_col_in_fullres", "pxl_row_in_fullres")
    for c in coords:
        tp[c] = pd.to_numeric(tp[c], errors="coerce")
    ok = tp[list(coords)].notna().all(axis=1)
    a = a[ok.values, :].copy()
    tp = tp.loc[ok]
    a.obsm["spatial"] = tp[list(coords)].values                                 # full-resolution pixels
    for c in tp.columns:
        if c not in coords:
            a.obs[c] = tp[c].values                                             # in_tissue, array_row, array_col
    sf = spatial / "scalefactors_json.json"
    a.uns["spatial"] = {sample_name: {"scalefactors": json.loads(sf.read_text()) if sf.exists() else {}}}
    a.obs["sample"] = sample_name
    a.obs["condition"] = SAMPLES[sample_name]["condition"]
    return a


def qc_metrics(a):
    """
    [SHAMIM, MODIFIED] step 2, "Basic quality control".
    Original: sc.pp.calculate_qc_metrics(adata, inplace=True), which does not define mitochondrial genes, so the
    pct_counts_mt used in step 3 did not exist. Mitochondrial genes are now flagged by the mt- prefix, case-insensitive.
    """
    a.var["mt"] = a.var_names.str.lower().str.startswith("mt-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    return a


def qc_sweep(a, name):
    """[RUN] Percentage of bins retained at each minimum gene count (Supplementary Figure S1, right panel)."""
    rows = []
    for t in [0, 5, 10, 25, 50, 100, 200]:
        k = a.obs["n_genes_by_counts"] > t
        rows.append({"sample": name, "min_genes": t, "n_retained": int(k.sum()), "pct_retained": round(100 * k.mean(), 2),
                     "median_genes": float(a.obs.loc[k, "n_genes_by_counts"].median()) if k.any() else np.nan})
    return pd.DataFrame(rows)


def qc_filter(a):
    """
    [SHAMIM, MODIFIED] step 3, "Filter low-quality spots".
    Original: more than 200 genes and less than 5% mitochondrial counts.
    Now: more than 10 genes, more than 15 UMI and less than 20% mitochondrial counts. The two sections differ about
    two-fold in capture, and a gene threshold carried over from single-cell work removes bins preferentially from the
    lower-capture control section (more than 100 genes kept 0.99% of control and 3.85% of tumor-bearing 8 um bins),
    which alone produced an apparent absence of osteoclast signal in the control (Methods; Supplementary Figure S1).
    """
    keep = ((a.obs["n_genes_by_counts"] > MIN_GENES) & (a.obs["total_counts"] > MIN_COUNTS)
            & (a.obs["pct_counts_mt"] < MAX_PCT_MT))
    return a[keep.values].copy()


# ============================================================================
# 2. NORMALIZATION, EMBEDDING AND CLUSTERING (Leiden domains: Supplementary Figure S5 only)
# ============================================================================

def normalize_and_cluster(a, resolution=0.3):
    """
    Normalization, highly variable genes, PCA, Harmony, neighbor graph, Leiden and UMAP on the joint object.
    Each step is tagged below.
    """
    a.layers["counts"] = a.X.copy()                                             # [RUN] raw counts kept
    sc.pp.normalize_total(a, target_sum=1e4)                                    # [SHAMIM] step 4
    sc.pp.log1p(a)                                                              # [SHAMIM] step 4
    a.layers["lognorm"] = a.X.copy()                                            # [RUN] every score and test uses this
    a.raw = a

    # [SHAMIM, MODIFIED] step 5. Original: flavor="seurat", n_top_genes=2000, then adata = adata[:, highly_variable],
    # which discards every other gene. Now computed per section (batch_key="sample"), and only a separate copy is
    # subset and scaled, so that gene-set scoring and differential expression keep the full log-normalized matrix.
    sc.pp.highly_variable_genes(a, flavor="seurat", n_top_genes=2000, batch_key="sample")
    tmp = a[:, a.var.highly_variable].copy()
    sc.pp.scale(tmp, max_value=10)                                              # [SHAMIM] step 6
    sc.tl.pca(tmp, svd_solver="arpack", random_state=RANDOM_SEED)              # [SHAMIM] step 6 (seed added)
    a.obsm["X_pca"] = tmp.obsm["X_pca"]
    del tmp

    # [RUN] Harmony on the section label (residual between-section variation; Methods).
    import harmonypy
    ho = harmonypy.run_harmony(a.obsm["X_pca"], a.obs, ["sample"], random_state=RANDOM_SEED)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != a.n_obs:          # harmonypy 0.0.x returns (n_pcs, n_cells)
        Z = Z.T
    a.obsm["X_pca_harmony"] = Z

    # [SHAMIM, MODIFIED] step 7. Original: sc.pp.neighbors(n_neighbors=10, n_pcs=30); sc.tl.leiden(resolution=0.5).
    # Now 15 neighbors on the Harmony embedding (30 components as before), and Leiden with the modularity objective set
    # explicitly: flavor="igraph" defaults to the CPM objective, on whose scale resolution 0.5 produced 687 clusters of
    # median 23 bins. Resolution 0.3.
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony", random_state=RANDOM_SEED)
    sc.tl.leiden(a, resolution=resolution, key_added="leiden", flavor="igraph", objective_function="modularity",
                 n_iterations=2, random_state=RANDOM_SEED, directed=False)
    sc.tl.umap(a, random_state=RANDOM_SEED)                                     # [SHAMIM] step 8 (embedding only)
    return a


def domain_markers(a, key="leiden", n_top=30):
    """
    [SHAMIM, MODIFIED] steps 10-11, "Differential expression between clusters" and "Export DEGs".
    Original: sc.tl.rank_genes_groups(groupby='leiden', method='wilcoxon') on the HVG-subset, scaled matrix.
    Now on the full log-normalized matrix (use_raw=False, layer copied into X), with detection fractions (pts=True).
    Benjamini-Hochberg adjustment is scanpy's default, as before.
    """
    a.X = a.layers["lognorm"].copy()
    sc.tl.rank_genes_groups(a, groupby=key, method="wilcoxon", use_raw=False, pts=True)
    df = sc.get.rank_genes_groups_df(a, group=None)
    return (df.sort_values(["group", "scores"], ascending=[True, False]).groupby("group").head(n_top)
              .reset_index(drop=True))


def domain_screen(a, de, key="leiden", n_top=10):
    """
    [RUN] Per Leiden domain: the fraction of its top markers that are ubiquitously expressed transcripts (the only
    merging criterion), and median detected genes per UMI. A domain is flagged when shallow sampling (>= 0.94 genes per
    UMI) combines with non-specific markers; flagged domains are not merged (Methods).
    """
    cplx = a.obs["n_genes_by_counts"] / a.obs["total_counts"].clip(lower=1)
    med_cplx = cplx.groupby(a.obs[key], observed=True).median()
    rows = []
    for g, sub in de.groupby("group"):
        top = sub.nlargest(n_top, "scores")["names"].tolist()
        rows.append({"domain": g, "n_bins": int((a.obs[key].astype(str) == str(g)).sum()),
                     "ubiquitous_fraction": round(sum(t in UBIQUITOUS for t in top) / max(len(top), 1), 2),
                     "median_genes_per_umi": round(float(med_cplx.get(g, np.nan)), 3), "top_markers": ",".join(top[:5])})
    df = pd.DataFrame(rows)
    df["de_uninterpretable"] = (df["median_genes_per_umi"] >= 0.94) & (df["ubiquitous_fraction"] >= 0.2)
    return df


def collapse_unresolved(a, screen, key="leiden", ubiq_frac=0.3, label="unresolved"):
    """
    [RUN] Merge domains whose top markers are ubiquitously expressed transcripts into 'unresolved' (Methods;
    Supplementary Figure S5: 50,273 bins, 33.1%). 0.3 is the natural break in this dataset (merged domains 0.3-0.6,
    next value 0.1). The mechanism is not established, and 'unresolved' is not a tumor label.
    """
    hits = screen.loc[screen["ubiquitous_fraction"] >= ubiq_frac, "domain"].astype(str).tolist()
    new = a.obs[key].astype(str)
    a.obs[key] = pd.Categorical(new.where(~new.isin(hits), label))
    print(f"[domains] merged {len(hits)} domains ({int(new.isin(hits).sum()):,} bins) into '{label}'")
    return a


def collapse_fragments(a, key="leiden", min_bins=100, label="unassigned"):
    """[RUN] Domains of fewer than 100 bins (kNN-graph debris) relabeled 'unassigned' (Methods)."""
    vc = a.obs[key].value_counts()
    small = {str(x) for x in vc[vc < min_bins].index}
    if small:
        a.obs[key] = pd.Categorical(a.obs[key].astype(str).where(~a.obs[key].astype(str).isin(small), label))
    return a


def clustering_stability(a, resolution=0.3, seeds=(0, 1, 2, 3, 4), n_neighbors=15, n_pcs=30):
    """
    [RUN] Adjusted Rand index between Leiden runs that differ only by seed (graph and partition re-seeded together).
    Methods and Supplementary Figure S5: "across five seeds the adjusted Rand index was 0.72" (median of the pairs).
    """
    from sklearn.metrics import adjusted_rand_score
    labs = {}
    for sd in seeds:
        sc.pp.neighbors(a, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep="X_pca_harmony", random_state=sd,
                        key_added=f"n{sd}")
        sc.tl.leiden(a, resolution=resolution, key_added=f"_s{sd}", flavor="igraph", objective_function="modularity",
                     n_iterations=2, random_state=sd, directed=False, neighbors_key=f"n{sd}")
        labs[sd] = a.obs[f"_s{sd}"].astype(str).values
    rows = [{"seed_a": i, "seed_b": j, "n_clusters_a": len(set(labs[i])), "n_clusters_b": len(set(labs[j])),
             "ARI": round(adjusted_rand_score(labs[i], labs[j]), 4)}
            for k, i in enumerate(seeds) for j in seeds[k + 1:]]
    for sd in seeds:
        del a.obs[f"_s{sd}"], a.uns[f"n{sd}"], a.obsp[f"n{sd}_distances"], a.obsp[f"n{sd}_connectivities"]
    df = pd.DataFrame(rows)
    print(f"[stability] median ARI {df['ARI'].median():.3f}")
    return df


def map_domains_from_coarser(a_fine, coarse, fine_um, coarse_um, keys=("leiden", "program")):
    """
    [RUN] Project 16 um labels onto 8 um bins by array geometry alone (Methods: at 8 um many bins share identical PCA
    coordinates and the graph fragments, so clustering is done at 16 um). Space Ranger bins nest exactly: an 8 um array
    index maps to the 16 um index by integer division by 2. Fine bins whose parent was removed at QC are 'unassigned'.
    """
    f = int(round(coarse_um / fine_um))
    for key in keys:
        lut = {(s, int(r), int(c)): lab for s, r, c, lab in zip(coarse.obs["sample"].astype(str), coarse.obs["array_row"],
                                                                 coarse.obs["array_col"], coarse.obs[key].astype(str))}
        a_fine.obs[key] = pd.Categorical([lut.get((s, int(r) // f, int(c) // f), "unassigned") for s, r, c in
                                          zip(a_fine.obs["sample"].astype(str), a_fine.obs["array_row"], a_fine.obs["array_col"])])
    return a_fine


# ============================================================================
# 3. GENE-SET SCORING AND POSITIVITY  [RUN]
# ============================================================================

def present(a, genes):
    return [g for g in genes if g in a.var_names]


def score_all(a, sigs=SIG):
    """
    [RUN] Per-bin score of every signature with scanpy score_genes (expression-matched background of 50 genes) on the
    log-normalized matrix, and the number of core genes detected per bin.
    """
    a.X = a.layers["lognorm"].copy()
    cov = []
    for name, d in sigs.items():
        ext, core = present(a, d["extended"]), present(a, d["core"])
        cov.append({"signature": name, "extended_detected": ",".join(ext), "core_detected": ",".join(core),
                    "missing": ",".join(g for g in d["extended"] if g not in a.var_names)})
        if len(ext) >= 2:
            sc.tl.score_genes(a, ext, score_name=f"score_{name}", ctrl_size=50, random_state=RANDOM_SEED)
        else:
            a.obs[f"score_{name}"] = np.nan
        if core:
            X = a[:, core].X
            X = X.toarray() if sp.issparse(X) else np.asarray(X)
            a.obs[f"ncore_{name}"] = (X > 0).sum(axis=1)
        else:
            a.obs[f"ncore_{name}"] = 0
    a.uns["signature_coverage"] = pd.DataFrame(cov)
    return a


def null_threshold(a, name, sigs=SIG, pct=95.0, n_perm=50, n_sub=20000):
    """
    [RUN] 95th percentile of a null of randomly drawn gene sets of the same size, scored the same way on 20,000
    sampled bins. Bins whose null score is exactly 0 (no gene of either set detected) are excluded, otherwise the
    threshold collapses to 0.
    """
    genes = present(a, sigs[name]["extended"])
    if len(genes) < 2:
        return np.nan
    rng = np.random.default_rng(RANDOM_SEED)
    sub = a[rng.choice(a.n_obs, min(n_sub, a.n_obs), replace=False)].copy()
    sub.X = sub.layers["lognorm"].copy()
    nulls = []
    for _ in range(n_perm):
        rg = list(rng.choice(sub.var_names, len(genes), replace=False))
        sc.tl.score_genes(sub, rg, score_name="_n", ctrl_size=50, random_state=int(rng.integers(1e6)))
        v = sub.obs["_n"].values
        nulls.append(v[v != 0])
    nulls = np.concatenate(nulls) if any(len(x) for x in nulls) else np.array([])
    nulls = nulls[nulls != 0]
    return float(np.percentile(nulls, pct)) if nulls.size else 0.0


def call_positive(a, names=PRIMARY, sigs=SIG, min_core=2, n_perm=50):
    """
    [RUN] A bin is positive for a signature when its score exceeds the null threshold, at least two core genes are
    detected, and no designated competing signature scores higher (Methods; Table 1).
    """
    thr = {}
    for n in names:
        key = f"score_{n}"
        if key not in a.obs or a.obs[key].isna().all():
            continue
        t = null_threshold(a, n, sigs=sigs, n_perm=n_perm)
        thr[n] = t
        base = (a.obs[key] > t) & (a.obs[f"ncore_{n}"] >= min_core)
        for comp in sigs[n]["competes"]:
            if f"score_{comp}" in a.obs:
                base &= a.obs[key] >= a.obs[f"score_{comp}"]
        a.obs[f"pos_{n}"] = base.fillna(False)
    a.uns["positivity_thresholds"] = thr
    return a


def positivity_table(a, names=PRIMARY, roi_key=None):
    """[RUN] Positive bins per signature and section, with Wilson 95% intervals (Figures 4g and 6a, left panels)."""
    from statsmodels.stats.proportion import proportion_confint
    obs = a.obs if roi_key is None else a.obs[a.obs[roi_key]]
    rows = []
    for n in names:
        if f"pos_{n}" not in obs:
            continue
        for s, sub in obs.groupby("sample", observed=True):
            x, tot = int(sub[f"pos_{n}"].sum()), len(sub)
            lo, hi = proportion_confint(x, tot, method="wilson")
            rows.append({"signature": n, "sample": s, "roi": roi_key or "whole_section", "n_positive": x, "n_bins": tot,
                         "pct_positive": round(100 * x / tot, 3), "ci_lo": round(100 * lo, 3), "ci_hi": round(100 * hi, 3),
                         "mean_score": float(sub[f"score_{n}"].mean())})
    df = pd.DataFrame(rows)
    w = df.pivot_table(index="signature", columns="sample", values="pct_positive")
    w["fold_tumor_over_control"] = (w["RL-1-5"] / w["LL-2-1"].replace(0, np.nan)).round(2)   # not quoted in the text
    return df.merge(w["fold_tumor_over_control"].reset_index(), on="signature", how="left")


def depth_match_target(a):
    """[RUN] Median UMI per bin of the shallower section, the common depth of every between-section comparison."""
    med = a.obs.groupby("sample", observed=True)["total_counts"].median()
    return int(med.min()), med


def downsampled_counts(a, target, rep):
    """[RUN] One depth-matching replicate: raw counts downsampled to at most `target` UMI per bin, seed 0 + rep."""
    b = ad.AnnData(a.layers["counts"].copy(), obs=a.obs[[]].copy())
    sc.pp.downsample_counts(b, counts_per_cell=target, random_state=RANDOM_SEED + rep)
    return b.X


def depth_matched(a, names=PRIMARY, n_reps=5, roi_key=None):
    """[RUN] Positivity after depth matching, mean and SD over five replicates (Figures 4g and 6a, right panels)."""
    target, _ = depth_match_target(a)
    print(f"[depth-match] downsampling to {target} counts/bin")
    out = []
    for rep in range(n_reps):
        b = a.copy()
        b.uns.pop("log1p", None)
        b.X = downsampled_counts(a, target, rep)
        sc.pp.normalize_total(b, target_sum=1e4)
        sc.pp.log1p(b)
        b.layers["lognorm"] = b.X.copy()
        b = call_positive(score_all(b), names, n_perm=20)
        t = positivity_table(b, names, roi_key=roi_key)
        t["rep"] = rep
        out.append(t)
    df = pd.concat(out, ignore_index=True)
    return (df.groupby(["signature", "sample"], as_index=False)
              .agg(pct_positive_mean=("pct_positive", "mean"), pct_positive_sd=("pct_positive", "std"),
                   mean_score=("mean_score", "mean")))


# ============================================================================
# 4. REGION OF INTEREST  [RUN]
# ============================================================================

def define_marrow_roi(a, exclude=("skeletal_muscle", "nerve", "smooth_muscle"), drawn=None):
    """
    [RUN] Marrow and bone region of interest. Transcriptional part: skeletal muscle, peripheral nerve and vascular
    smooth muscle have no positivity call, so a bin is excluded when its score for any of the three is above that
    signature's 90th percentile over all bins of both sections. When the hand-drawn ROI exists
    (bone_pending_masks.rasterise_drawn_roi), it is intersected with the transcriptional part and is then the only thing
    that excludes the diaphysis. All results in the current manuscript use the transcriptional ROI (drawn=None).
    """
    keep = np.ones(a.n_obs, dtype=bool)
    for e in exclude:
        if f"pos_{e}" in a.obs:
            keep &= ~a.obs[f"pos_{e}"].values
        else:
            keep &= (a.obs[f"score_{e}"] <= np.nanpercentile(a.obs[f"score_{e}"], 90)).values
    a.obs["in_roi_transcriptional"] = keep
    if drawn is not None:
        a.obs["in_roi_drawn"] = np.asarray(drawn, dtype=bool)
        keep = keep & a.obs["in_roi_drawn"].values
    a.obs["in_marrow_roi"] = keep
    return a


# ============================================================================
# 5. DETERMINISTIC TISSUE ASSIGNMENT (Figure 3)  [RUN]
# ============================================================================

def spatially_pooled_scores(a, pool_k=8, depth_match=False):
    """
    [RUN] Raw counts summed over each bin and its eight nearest spatial neighbors before normalization and scoring
    (about ninefold depth, still deterministic). depth_match=True first downsamples every bin to the median of the
    shallower section. Writes pooled_* (native depth) or pooled_dm_* (depth-matched) columns.
    """
    from scipy.spatial import cKDTree
    b = a.copy()
    b.uns.pop("log1p", None)
    b.X = b.layers["counts"].copy()
    if depth_match:
        tgt = int(b.obs.groupby("sample", observed=True)["total_counts"].median().min())
        sc.pp.downsample_counts(b, counts_per_cell=tgt, random_state=RANDOM_SEED)
    X = b.X.tocsr() if sp.issparse(b.X) else sp.csr_matrix(b.X)
    pooled = sp.lil_matrix(X.shape, dtype=np.float32)
    for smp in b.obs["sample"].astype(str).unique():
        m = np.where((b.obs["sample"].astype(str) == smp).values)[0]
        _, idx = cKDTree(b.obsm["spatial"][m]).query(b.obsm["spatial"][m], k=pool_k + 1)
        rows = np.repeat(np.arange(len(m)), pool_k + 1)
        A = sp.csr_matrix((np.ones(len(rows), np.float32), (rows, idx.ravel())), shape=(len(m), len(m)))
        pooled[m] = A @ X[m]
    b.X = pooled.tocsr()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    b.layers["lognorm"] = b.X.copy()
    b = call_positive(score_all(b), PRIMARY, n_perm=20)
    prefix = "pooled_dm_" if depth_match else "pooled_"
    for c in b.obs.columns:
        if c.startswith(("score_", "ncore_", "pos_")):
            a.obs[f"{prefix}{c}"] = b.obs[c].values
    a.uns[f"{prefix}thresholds"] = b.uns.get("positivity_thresholds", {})
    return a


def assign_dominant_program(a, prefix, label_col, programs=PROGRAMS, min_core=2, ambiguous="ambiguous"):
    """
    [RUN] Each bin takes the highest-scoring signature that clears its null threshold with at least two core genes
    detected, otherwise 'ambiguous'. No graph, no seed, no resolution parameter (Methods; Figure 3).
    Thresholds for signatures outside PRIMARY (muscle, nerve, smooth muscle) are 0.
    """
    cols = [f"{prefix}score_{p}" for p in programs if f"{prefix}score_{p}" in a.obs]
    S = a.obs[cols].to_numpy(dtype=float)
    thr_src = a.uns.get(f"{prefix}thresholds", {})
    names = np.array([c.replace(f"{prefix}score_", "") for c in cols])
    thr = np.array([thr_src.get(n, 0.0) for n in names])
    core = np.array([a.obs[f"{prefix}ncore_{n}"].to_numpy() for n in names]).T
    ok = (S > thr) & (core >= min_core)
    best = np.nanargmax(np.where(ok, S, -np.inf), axis=1)
    a.obs[label_col] = pd.Categorical(np.where(ok.any(axis=1), names[best], ambiguous))
    return a


def program_composition(a, col):
    """[RUN] Bins per assigned signature and section, with both denominators and the ambiguous rate (Figure 3b)."""
    ct = pd.crosstab(a.obs["sample"], a.obs[col])
    rows = []
    for smp, r in ct.iterrows():
        tot = r.sum()
        assigned = tot - r.get("ambiguous", 0)
        for prog, n in r.items():
            rows.append({"sample": smp, "program": prog, "n_bins": int(n), "pct_of_all_bins": round(100 * n / tot, 3),
                         "pct_of_assigned_bins": round(100 * n / assigned, 3) if prog != "ambiguous" and assigned else np.nan,
                         "sample_ambiguous_pct": round(100 * r.get("ambiguous", 0) / tot, 2)})
    return pd.DataFrame(rows)


# ============================================================================
# 6. OSTEOCLAST COMPARTMENTS (Figure 4h)  [RUN, with SHAMIM step 10 for the test]
# ============================================================================

def split_osteoclast_compartments(a, pct=90):
    """
    [RUN] Osteoclast-positive bins split by chondrocyte score: above the 90th percentile of all bins = growth-plate
    (chondro-osseous junction) compartment, the rest = tumor-associated (Methods; Figure 4h).
    """
    a.obs["growth_plate"] = a.obs["score_chondrocyte"] > np.nanpercentile(a.obs["score_chondrocyte"], pct)
    a.obs["oc_compartment"] = "negative"
    oc = a.obs["pos_osteoclast"].fillna(False)
    a.obs.loc[oc & a.obs["growth_plate"], "oc_compartment"] = "growth_plate"
    a.obs.loc[oc & ~a.obs["growth_plate"], "oc_compartment"] = "tumor_associated"
    return a


def compare_osteoclast_states(a):
    """[SHAMIM, MODIFIED] step 10 (Wilcoxon rank-sum, BH-adjusted), applied to the two osteoclast compartments
    instead of Leiden clusters, on the full log-normalized matrix."""
    sub = a[a.obs["oc_compartment"] != "negative"].copy()
    sub.X = sub.layers["lognorm"].copy()
    sc.tl.rank_genes_groups(sub, groupby="oc_compartment", method="wilcoxon", use_raw=False)
    return sc.get.rank_genes_groups_df(sub, group=None)


# ============================================================================
# 7. HOST RECEPTOR PANEL AND ADIPOCYTE GENES (Figure 6b, Supplementary Figure S4)  [RUN]
# ============================================================================

def receptor_panel(a, panel=RECEPTOR_PANEL, n_reps=5, control="LL-2-1", tumour="RL-1-5", roi_key="in_marrow_roi"):
    """
    [RUN] Percentage of ROI bins detecting each gene, raw and after depth matching (same downsampling as
    depth_matched(), five replicates), with Wilson intervals on the depth-matched rate. Used for the host receptor panel
    (Figure 6b) and, with panel=ADIPOCYTE_GENES, for Supplementary Figure S4.
    """
    from statsmodels.stats.proportion import proportion_confint
    target, med = depth_match_target(a)
    a.X = a.layers["lognorm"].copy()
    for g in panel:
        if g in a.var_names:
            v = a[:, g].X
            a.obs[f"expr_{g}"] = np.asarray(v.todense() if sp.issparse(v) else v).ravel()
    genes = [g for g in panel if g in a.var_names]
    gi = [a.var_names.get_loc(g) for g in genes]
    smp = a.obs["sample"].astype(str).values.copy()
    samples = list(pd.unique(smp))
    smp[~np.asarray(a.obs[roi_key], bool)] = "__outside_roi__"
    det = {s: np.zeros((n_reps, len(genes))) for s in samples}
    mean = {s: np.zeros((n_reps, len(genes))) for s in samples}
    for rep in range(n_reps):
        X = downsampled_counts(a, target, rep)
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        lib = np.asarray(X.sum(axis=1)).ravel()
        Xg = X[:, gi].toarray()
        ln = np.log1p(Xg / np.maximum(lib, 1)[:, None] * 1e4)
        for s in samples:
            m = smp == s
            det[s][rep] = (Xg[m] > 0).mean(axis=0)
            mean[s][rep] = ln[m].mean(axis=0)
    rows = []
    for j, g in enumerate(genes):
        for s in samples:
            sub = a.obs.loc[smp == s, f"expr_{g}"]
            n = int(len(sub))
            p_dm = float(det[s][:, j].mean())
            npos = int(round(p_dm * n))
            lo, hi = proportion_confint(npos, n, method="wilson")
            rows.append({"gene": g, "role": panel[g], "sample": s, "roi": roi_key, "mean_raw": float(sub.mean()),
                         "pct_detected_raw": 100 * float((sub > 0).mean()), "mean_dm": float(mean[s][:, j].mean()),
                         "pct_detected_dm": 100 * p_dm, "n_positive_dm": npos, "n_bins": n, "ci_lo_dm": 100 * lo,
                         "ci_hi_dm": 100 * hi, "median_umi_per_bin": float(med.get(s, np.nan))})
    df = pd.DataFrame(rows)
    df["fold_dm"] = np.nan
    for g in genes:
        c = df[(df.gene == g) & (df["sample"] == control)]["pct_detected_dm"]
        t = df[(df.gene == g) & (df["sample"] == tumour)]
        if len(c) and len(t) and float(c.iloc[0]) > 0:
            df.loc[t.index, "fold_dm"] = t["pct_detected_dm"].values / float(c.iloc[0])
    return df


# ============================================================================
# 8. SPATIAL ANALYSES ON THE JOINT OBJECT
# ============================================================================

def _array_xy_um(adata, bin_um):
    """[RUN] Bin positions in um from the array grid (adjacent bins are exactly bin_um apart; no pixel calibration)."""
    return np.column_stack([adata.obs["array_col"].to_numpy(float), adata.obs["array_row"].to_numpy(float)]) * bin_um


def uncoupling_metric(adata, a_key="pooled_dm_pos_osteoblast", b_key="pooled_dm_pos_osteoclast", bin_um=16.0,
                      n_perm=199, roi_key="in_marrow_roi", min_pos=10, reference_sample="LL-2-1"):
    """
    [RUN] Figure 4i. Distance from each osteoblast-positive bin to the nearest osteoclast-positive bin (pooled,
    depth-matched calls), within the ROI of each section, against a random-placement null that draws the same numbers of
    positive bins in the same ROI (199 replicates). Primary statistic: observed mean / null mean per limb (< 1 closer than
    chance, > 1 farther); between limbs, the ratio of those ratios. Raw distances are density-dominated and are not
    compared between sections. Returns (summary, distributions); null counts in `distributions` are per replicate.
    """
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(RANDOM_SEED)
    res, dist_rows = [], []

    def _hist(d, s, source, scale=1.0):
        v, c = np.unique(np.round(d, 2), return_counts=True)
        dist_rows.extend({"sample": s, "source": source, "distance_um": float(vi), "count": float(ci) * scale}
                         for vi, ci in zip(v, c))

    for s in adata.obs["sample"].unique():
        m = (adata.obs["sample"] == s).values & np.asarray(adata.obs[roi_key], bool)
        sub = adata[m]
        xy = _array_xy_um(sub, bin_um)
        A = np.asarray(sub.obs[a_key].fillna(False), bool)
        B = np.asarray(sub.obs[b_key].fillna(False), bool)
        nA, nB = int(A.sum()), int(B.sum())
        if nA < min_pos or nB < min_pos:
            res.append({"sample": s, "n_a": nA, "n_b": nB, "status": f"not estimable (n_positive < {min_pos})"})
            continue
        d_obs = cKDTree(xy[B]).query(xy[A])[0]
        mean = float(np.mean(d_obs))
        _hist(d_obs, s, "observed")
        null_mean, pooled = [], []
        for _ in range(n_perm):
            ia = rng.choice(len(xy), nA, replace=False)
            ib = rng.choice(len(xy), nB, replace=False)
            dn = cKDTree(xy[ib]).query(xy[ia])[0]
            null_mean.append(dn.mean())
            pooled.append(dn)
        null_mean = np.array(null_mean)
        _hist(np.concatenate(pooled), s, "null", scale=1.0 / n_perm)
        lo, hi = float(np.percentile(null_mean, 2.5)), float(np.percentile(null_mean, 97.5))
        res.append({"sample": s, "n_a": nA, "n_b": nB, "n_units_roi": int(len(xy)),
                    "ratio_obs_null_mean": mean / float(null_mean.mean()), "ratio_obs_null_lo": mean / hi,
                    "ratio_obs_null_hi": mean / lo, "p_greater_than_null_mean": float((np.sum(null_mean >= mean) + 1) / (n_perm + 1)),
                    "median_nn_um": float(np.median(d_obs)), "mean_nn_um": mean, "null_mean_um": float(null_mean.mean()),
                    "status": "estimable"})
    df = pd.DataFrame(res)
    ref = df[(df["sample"] == reference_sample) & (df["status"] == "estimable")]
    if not ref.empty:
        r0 = ref.iloc[0]
        df["ratio_of_ratios_vs_reference"] = df["ratio_obs_null_mean"] / r0["ratio_obs_null_mean"]
    df["note"] = "raw distances (median_nn_um, mean_nn_um) are density-dominated; do not compare between sections"
    return df, pd.DataFrame(dist_rows)


def proximity_tables(a):
    """
    [RUN] Results and Methods, "within 50 um": fraction of bins of one signature whose nearest bin of the other
    signature lies within 50 um (center to center), against the same fraction for all ROI bins.
      fibroblast-assigned -> osteoblast-assigned bins (deterministic assignment, column 'program'): 41.8/41.3/67.3/5.6%
      macrophage-positive -> osteoclast-positive bins (per-bin positivity): 115 of 405, 28.4% against 17.1%
    Positions are full-resolution pixels times each section's microns per pixel.
    """
    from scipy.spatial import cKDTree
    smp = a.obs["sample"].astype(str).values
    um = np.array([UM_PER_PX[s] for s in smp])[:, None] * np.asarray(a.obsm["spatial"])
    roi = a.obs["in_marrow_roi"].values.astype(bool)
    prog = a.obs["program"].astype(str).values
    fib, mac = [], []
    for s in SAMPLES:
        m = (smp == s) & roi
        P, pr = um[m], prog[m]
        d = cKDTree(P[pr == "osteoblast"]).query(P, k=1)[0]
        fib.append({"section": SECTION[s], "fibroblast_assigned_bins": int((pr == "fibroblast").sum()),
                    "osteoblast_assigned_bins": int((pr == "osteoblast").sum()), "region_of_interest_bins": int(m.sum()),
                    "percent_fibroblast_within_50um_of_osteoblast": round(100 * float((d[pr == "fibroblast"] <= 50).mean()), 1),
                    "percent_all_bins_within_50um_of_osteoblast": round(100 * float((d <= 50).mean()), 1)})
        oc = a.obs["pos_osteoclast"].values.astype(bool)[m]
        mp = a.obs["pos_macrophage"].values.astype(bool)[m]
        d = cKDTree(P[oc]).query(P, k=1)[0]
        mac.append({"section": SECTION[s], "osteoclast_positive_bins": int(oc.sum()), "macrophage_positive_bins": int(mp.sum()),
                    "macrophage_within_50um_of_osteoclast": int((d[mp] <= 50).sum()),
                    "percent_macrophage_within_50um": round(100 * float((d[mp] <= 50).mean()), 1),
                    "region_of_interest_bins": int(m.sum()), "percent_all_bins_within_50um": round(100 * float((d <= 50).mean()), 1)})
    return pd.DataFrame(fib), pd.DataFrame(mac)


def morans_i(adata, n_svg=2000, n_consensus=200):
    """
    [SHAMIM, MODIFIED] step 13, "Spatial statistics".
    Original: sq.gr.spatial_neighbors(adata) and sq.gr.spatial_autocorr(adata, mode="moran") on the whole object.
    Now computed PER SECTION (the two sections' coordinate ranges overlap, so a graph on the joint object would join
    bins on different slides; Methods), with a Delaunay graph on the bin coordinates, for the top 2,000 highly variable
    genes chosen once on the joint object. The combined statement is the set of genes in the top 200 by Moran's I in
    both sections (Supplementary Figure S2: 89 genes).
    """
    import squidpy as sq
    hv = adata.copy()
    hv.X = hv.layers["lognorm"].copy()
    sc.pp.highly_variable_genes(hv, flavor="seurat", n_top_genes=n_svg)
    genes = hv.var_names[hv.var.highly_variable].tolist()
    del hv
    samples = list(adata.obs["sample"].astype(str).unique())
    morans = []
    for s in samples:
        sub = adata[(adata.obs["sample"].astype(str) == s).values].copy()
        sub.X = sub.layers["lognorm"].copy()
        sq.gr.spatial_neighbors(sub, coord_type="generic", delaunay=True)
        # only Moran's I itself is reported, so no permutation P values (they do not change I)
        sq.gr.spatial_autocorr(sub, mode="moran", genes=genes, n_perms=None, n_jobs=1)
        m = sub.uns["moranI"].reset_index().rename(columns={"index": "gene"})
        m.insert(0, "sample", s)
        m["rank_in_sample"] = m["I"].rank(ascending=False, method="min", na_option="keep").astype("Int64")
        morans.append(m)
        del sub
    moran = pd.concat(morans, ignore_index=True)
    top = moran[(moran["rank_in_sample"] <= n_consensus).fillna(False).astype(bool)]
    shared = top.groupby("gene")["sample"].nunique()
    shared = shared[shared == len(samples)].index
    cons = moran[moran["gene"].isin(shared)].pivot_table(index="gene", columns="sample", values="I").add_prefix("moranI_")
    cons["mean_rank"] = moran[moran["gene"].isin(shared)].groupby("gene")["rank_in_sample"].mean()
    return moran, cons.sort_values("mean_rank").reset_index()


def leiden_signature_matrix(a):
    """
    [RUN] Supplementary Figure S5: mean signature score per Leiden domain, z-scored across domains. Each domain is
    named by the signature most enriched among its assigned bins relative to the whole dataset (at least 20 bins).
    """
    def lbl(c):
        return c.replace("_", " ").replace("neutrophil", "granulocyte")
    prog = a.obs["program"].astype(str).values
    ld = a.obs["leiden"].astype(str).values
    ct = pd.crosstab(pd.Series(ld, name="d"), pd.Series(prog, name="p")).drop(columns="ambiguous", errors="ignore")
    enr = (ct.div(ct.sum(1), axis=0) / (ct.sum(0) / ct.values.sum())).where(ct >= 20)
    names = {c: c if c in ("unresolved", "unassigned") or c not in enr.index or enr.loc[c].isna().all()
             else f"{lbl(enr.loc[c].idxmax())}-enriched (domain {c})" for c in pd.unique(ld)}
    cols = [c for c in a.obs.columns if c.startswith("score_") and c.replace("score_", "") in PRIMARY]
    M = a.obs.groupby("leiden", observed=True)[cols].mean()
    M = (M - M.mean()) / M.std()
    M.index = [names[str(i)] for i in M.index]
    M.columns = [c.replace("score_", "").replace("neutrophil", "granulocyte") for c in M.columns]
    un = ld == "unresolved"
    return (M.reset_index().rename(columns={"index": "leiden_domain"}),
            pd.DataFrame([{"unresolved_bins": int(un.sum()), "all_bins": int(len(ld)), "percent": round(100 * un.mean(), 1)}]))


def geneset_sensitivity(a, n_reps=5):
    """
    [RUN] Methods: bone-lineage positivity recomputed under three gene sets with the calling rule held constant
    (canonical two-gene markers; the lists of the original study design; the scoring signatures of Table 1), at native
    depth and after depth matching, within the ROI. Slow (about 30 runs of scoring and null thresholds).
    """
    from statsmodels.stats.proportion import proportion_confint
    oc8 = ["Acp5", "Ctsk", "Calcr", "Mmp9", "Atp6v0d2", "Dcstamp", "Ocstamp", "Nfatc1"]
    ob7 = ["Runx2", "Sp7", "Alpl", "Col1a1", "Col1a2", "Spp1", "Bglap"]
    variants = {"canonical_two_gene": ({"core": ["Atp6v0d2", "Ctsk"], "extended": ["Atp6v0d2", "Ctsk"]},
                                       {"core": ["Col1a1", "Sp7"], "extended": ["Col1a1", "Sp7"]}),
                "original_study_design": ({"core": oc8, "extended": oc8}, {"core": ob7, "extended": ob7}),
                "table1_signatures": (None, None)}
    names, comp = ["osteoclast", "osteoblast"], ["erythroid", "neutrophil", "fibroblast", "chondrocyte"]
    roi = a.obs["in_marrow_roi"].values.astype(bool)
    smp = a.obs["sample"].astype(str).values
    target, _ = depth_match_target(a)
    rows = []
    for vname, (oc, ob) in variants.items():
        sigs = copy.deepcopy(SIG)
        if oc:
            sigs["osteoclast"].update(oc)
        if ob:
            sigs["osteoblast"].update(ob)
        sigs = {k: sigs[k] for k in names + comp}
        b = call_positive(score_all(a.copy(), sigs), names, sigs=sigs, n_perm=20)
        nat = {(n, s): int(b.obs[f"pos_{n}"].values[roi & (smp == s)].sum()) for n in names for s in SAMPLES}
        dm = {k: [] for k in nat}
        for rep in range(n_reps):
            c = a.copy()
            c.uns.pop("log1p", None)
            c.X = downsampled_counts(a, target, rep)
            sc.pp.normalize_total(c, target_sum=1e4)
            sc.pp.log1p(c)
            c.layers["lognorm"] = c.X.copy()
            c = call_positive(score_all(c, sigs), names, sigs=sigs, n_perm=20)
            for (n, s) in dm:
                dm[(n, s)].append(int(c.obs[f"pos_{n}"].values[roi & (smp == s)].sum()))
        for (n, s), k in nat.items():
            N = int((roi & (smp == s)).sum())
            lo, hi = proportion_confint(k, N, method="wilson")
            d = np.array(dm[(n, s)]) / N * 100
            rows.append({"gene_set": vname, "signature": n, "sample": s, "n_units_roi": N, "n_positive_native": k,
                         "pct_native": round(100 * k / N, 4), "ci_lo": round(100 * lo, 4), "ci_hi": round(100 * hi, 4),
                         "pct_depth_matched_mean": round(d.mean(), 4), "pct_depth_matched_sd": round(d.std(ddof=1), 4)})
    return pd.DataFrame(rows)


def table1_signatures():
    """[RUN] Table 1 as a CSV, generated from SIG so that the table and the code cannot diverge."""
    nm = {"neutrophil": "granulocyte", "endothelium": "endothelial", "skeletal_muscle": "skeletal muscle",
          "smooth_muscle": "vascular smooth muscle", "nerve": "peripheral nerve"}
    canon = {"osteoclast": "Atp6v0d2, Ctsk", "osteoblast": "Col1a1, Sp7"}
    return pd.DataFrame([{"signature": nm.get(k, k), "canonical_markers": canon.get(k, ""),
                          "signature_genes": ", ".join(d["extended"]), "genes_counted_toward_two_gene_rule": ", ".join(d["core"]),
                          "excluded_genes": "; ".join(f"{g}: {r}" for g, r in d["excluded"].items())} for k, d in SIG.items()])


def per_bin_table(a):
    """[RUN] One row per 16 um bin (Manuscript/data/per_bin): position in um, ROI flag, assignment, positivity."""
    o = a.obs
    smp = o["sample"].astype(str).values
    um = np.array([UM_PER_PX[s] for s in smp])[:, None] * np.asarray(a.obsm["spatial"])
    group = {"osteoclast": "bone lineage", "osteoblast": "bone lineage", "osteocyte": "bone lineage", "chondrocyte": "chondrocyte",
             "neutrophil": "hematopoietic marrow", "erythroid": "hematopoietic marrow", "megakaryocyte": "hematopoietic marrow",
             "macrophage": "hematopoietic marrow", "fibroblast": "stroma", "adipocyte": "stroma",
             "endothelium": "removed", "skeletal_muscle": "removed", "nerve": "removed", "smooth_muscle": "removed",
             "ambiguous": "unassigned"}
    pb = pd.DataFrame({"barcode": o.index.astype(str), "sample_id": smp, "section": [SECTION[s] for s in smp],
                       "x_um": um[:, 0].round(2), "y_um": um[:, 1].round(2),
                       "in_region_of_interest": o["in_marrow_roi"].values.astype(bool), "total_umi": o["total_counts"].values,
                       "assigned_signature": o["program"].astype(str).str.replace("neutrophil", "granulocyte").values,
                       "tissue_class": [group.get(p, "unassigned") for p in o["program"].astype(str)],
                       "leiden_domain": o["leiden"].astype(str).values})
    for n in ["osteoclast", "osteoblast", "osteocyte", "macrophage", "neutrophil", "erythroid", "megakaryocyte", "adipocyte", "fibroblast"]:
        pb[f"positive_{n.replace('neutrophil', 'granulocyte')}"] = o[f"pos_{n}"].fillna(False).values.astype(bool)
    for n in ["osteoclast", "osteoblast"]:
        pb[f"positive_{n}_pooled_depthmatched"] = o[f"pooled_dm_pos_{n}"].fillna(False).values.astype(bool)
    return pb


# ============================================================================
# MAIN
# ============================================================================

def load_sections(root, bin_str):
    """[SHAMIM, MODIFIED] steps 1-3 for each section; returns the filtered sections and the QC sweep."""
    adatas, sweeps = {}, []
    for name in SAMPLES:
        p = Path(root) / f"SR007144_{name}/outs/binned_outputs/square_{bin_str}"
        a = qc_metrics(read_visium_hd(p, name))
        sweeps.append(qc_sweep(a, name))
        a = qc_filter(a)
        print(f"[{name}] {a.n_obs:,} bins retained at {bin_str}")
        adatas[name] = a
    return adatas, pd.concat(sweeps, ignore_index=True)


def run_analysis(root, out, roi_masks=None, mask_frame="hires", do_moran=True, do_sensitivity=False, project_8um=False):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)

    # ---- 1. input and QC at 16 um -------------------------------------------- [SHAMIM, MODIFIED] steps 1-3
    adatas, sweep = load_sections(root, "016um")
    sweep.to_csv(out / "qc_threshold_sweep_16um.csv", index=False)
    a = ad.concat(adatas, label="sample_key", index_unique="-", uns_merge="unique", join="outer")   # [RUN] joint object
    a.obs["sample"] = a.obs["sample"].astype("category")

    # ---- 2. normalization, embedding, Leiden ---------------------------------- [SHAMIM, MODIFIED] steps 4-8
    a = normalize_and_cluster(a, resolution=0.3)

    # ---- 3. signature scores and positivity ----------------------------------- [RUN]
    a = score_all(a)
    a = call_positive(a)

    # ---- 4. region of interest ------------------------------------------------ [RUN]
    drawn = None
    if roi_masks:
        from bone_pending_masks import rasterise_drawn_roi                       # pending the hand-drawn ROI
        drawn = rasterise_drawn_roi(a, roi_masks, root, "016um", frame=mask_frame)
    a = define_marrow_roi(a, drawn=drawn)

    # ---- 5. host receptor panel (Figure 6b) and adipocyte genes (Supplementary Figure S4) ---------------- [RUN]
    receptor_panel(a).to_csv(out / "Figure6b_receptor_panel_depth_matched.csv", index=False)
    receptor_panel(a, panel={g: "adipocyte" for g in ADIPOCYTE_GENES}).to_csv(
        out / "FigureS4_adipocyte_genes_depth_matched.csv", index=False)

    # ---- 6. positivity, native depth and depth-matched (Figures 4g, 6a) ------ [RUN]
    positivity_table(a, roi_key="in_marrow_roi").to_csv(out / "Figure4g_6a_positivity_native_depth_roi.csv", index=False)
    depth_matched(a, PRIMARY, roi_key="in_marrow_roi").to_csv(out / "Figure4g_6a_positivity_depth_matched.csv", index=False)

    # ---- 7. Leiden domain markers and post-processing ------------------------ [SHAMIM, MODIFIED] steps 10-11; [RUN]
    de = domain_markers(a)
    screen = domain_screen(a, de)
    screen.to_csv(out / "leiden_domain_screen.csv", index=False)
    a = collapse_unresolved(a, screen)
    a = collapse_fragments(a)
    domain_markers(a).to_csv(out / "leiden_domain_markers.csv", index=False)  # [SHAMIM] step 11 (export DEGs)

    # ---- 8. osteoclast compartments (Figure 4h) ------------------------------- [RUN]; test [SHAMIM, MODIFIED] step 10
    a = split_osteoclast_compartments(a)
    compare_osteoclast_states(a).to_csv(out / "Figure4h_osteoclast_compartment_DE.csv", index=False)

    # ---- 9. deterministic tissue assignment (Figure 3) ------------------------ [RUN]
    a = spatially_pooled_scores(a)                                              # native depth -> 'program' (Fig 3a)
    a = assign_dominant_program(a, "pooled_", "program")
    a = spatially_pooled_scores(a, depth_match=True)                            # depth-matched -> 'program_dm' (Fig 3b)
    a = assign_dominant_program(a, "pooled_dm_", "program_dm")
    program_composition(a, "program").to_csv(out / "Figure3a_program_composition_native.csv", index=False)
    program_composition(a, "program_dm").to_csv(out / "Figure3b_program_composition_depth_matched.csv", index=False)

    # ---- 10. Leiden stability (Methods; Supplementary Figure S5) -------------- [RUN]
    clustering_stability(a).to_csv(out / "leiden_seed_stability_adjusted_rand_index.csv", index=False)
    M, un = leiden_signature_matrix(a)
    M.to_csv(out / "FigureS5_domain_by_signature_zscores.csv", index=False)
    un.to_csv(out / "FigureS5_unresolved_bins.csv", index=False)

    # ---- 11. spatial analyses -------------------------------------------------- [RUN]; Moran [SHAMIM, MODIFIED] step 13
    unc, dist = uncoupling_metric(a)
    unc.to_csv(out / "Figure4i_observed_to_null_summary.csv", index=False)
    dist.to_csv(out / "Figure4i_distance_distributions_observed_and_null.csv", index=False)
    fib, mac = proximity_tables(a)
    fib.to_csv(out / "fibroblast_to_osteoblast_proximity.csv", index=False)
    mac.to_csv(out / "macrophage_to_osteoclast_proximity.csv", index=False)
    if do_moran:
        moran, cons = morans_i(a)
        moran.to_csv(out / "morans_I_per_section.csv", index=False)
        cons.to_csv(out / "FigureS2_spatially_variable_genes_morans_I.csv", index=False)

    # ---- 12. Table 1, per-bin table, optional gene-set sensitivity ----------- [RUN]
    table1_signatures().to_csv(out / "Table1_gene_signatures.csv", index=False)
    per_bin_table(a).to_csv(out / "per_bin_signature_calls_16um.csv.gz", index=False, compression="gzip")
    if do_sensitivity:
        geneset_sensitivity(a).to_csv(out / "gene_set_sensitivity_three_gene_sets.csv", index=False)

    a.write(out / "joint_16um.h5ad")                                            # [SHAMIM] "Save final processed data"

    # ---- 13. optional: project the 16 um labels onto 8 um bins (Methods) ------ [RUN]
    if project_8um:
        fine, _ = load_sections(root, "008um")
        f = ad.concat(fine, label="sample_key", index_unique="-", join="outer")
        f.obs["sample"] = f.obs["sample"].astype(str)
        f = map_domains_from_coarser(f, a, 8, 16)
        f.obs[["sample", "array_row", "array_col", "leiden", "program"]].to_csv(out / "labels_8um_projected.csv")
    print(f"done -> {out}")
    return a


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", help="folder containing SR007144_LL-2-1/ and SR007144_RL-1-5/")
    ap.add_argument("--out", default="Manuscript/code/results")
    ap.add_argument("--roi-mask", action="append", default=None,
                    help="hand-drawn ROI as SAMPLE=PATH (.geojson on tissue_hires_image.png); repeat per sample")
    ap.add_argument("--mask-frame", default="hires", choices=["hires", "fullres", "lowres"])
    ap.add_argument("--skip-moran", action="store_true", help="skip Supplementary Figure S2 (slow)")
    ap.add_argument("--geneset-sensitivity", action="store_true", help="run the three-gene-set check of Methods (slow)")
    ap.add_argument("--project-8um", action="store_true", help="write the 16 um labels projected onto 8 um bins")
    args = ap.parse_args()
    masks = None
    if args.roi_mask:
        masks = dict(v.split("=", 1) for v in args.roi_mask)
    run_analysis(args.root, args.out, roi_masks=masks, mask_frame=args.mask_frame, do_moran=not args.skip_moran,
                 do_sensitivity=args.geneset_sensitivity, project_8um=args.project_8um)
