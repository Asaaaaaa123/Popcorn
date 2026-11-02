#!/usr/bin/env python3
"""
bone_analysis.py
================
Visium HD analysis pipeline: A549 intratibial lesion (right tibia, RL-1-5) against the
contralateral control (left tibia, LL-2-1) of the same NSG mouse, both sections on one capture slide, one CytAssist run,
Visium Mouse Transcriptome Probe Set v2.0 against mm10-2020-A. Every feature is a mouse probe, but probes of conserved
genes also capture the human A549 RNA, so the lesion is visible in the data.
[FIX] This paragraph said "Human A549 transcripts are not captured: every transcript is murine". The probe set is mouse
(all 19,059 features are ENSMUSG genes), but mouse probes for conserved genes cross-capture A549 RNA: Krt7 is detected
in 2,523 RL-1-5 bins against 17 LL-2-1 bins, Krt18 in 552 against 13, and Tmsb10 and Itm2b track the lesion. The
lesion is therefore visible inside the region of interest, and the tumor-bearing section's 'unresolved' Leiden domain
is strongly enriched for it (bone_sensitivity/ marks it by pooled Krt7 + Krt18 UMIs).

Reads the Space Ranger binned outputs and writes every table and number the manuscript reports from the transcriptome,
plus the joint AnnData object the figure script draws from.

[APP] Every dataset-specific value (sections, folders, parameters, gene signatures, which steps run) now comes from a
configuration (pipeline_config.py). Its defaults are the manuscript setup above, so the old command line still runs the
manuscript analysis; the web application (python -m webapp) edits and runs configurations.

PROVENANCE TAGS (every function and every step of run_analysis carries one)
----------------------------------------------------------------------------
  [SHAMIM]            taken unchanged from Shamim Mollah's original pipeline,
                      00_shamim_original/visium_spatial_processing_full_pipeline.py (2024-12-15); "step N" is the step
                      number in that file
  [SHAMIM, MODIFIED]  originates in that pipeline, but its parameters or implementation were changed by Run Zhang;
                      the original setting and the reason for the change are written next to it
  [RUN]               added by Run Zhang; not in the original pipeline
  [FIX]               defect corrected on 2026-10-28; the note says what was wrong and what it changed (CHANGES.md)
  [APP]               added on 2026-10-28 for the configurable pipeline and the web application

Steps of the original pipeline that are NOT used in the manuscript and are therefore not carried over:
  step 2 violin plot, step 6 PCA variance plot, step 8 UMAP plot, step 9 spatial cluster plot, step 10 marker plot,
  step 12 Enrichr GO/KEGG/Reactome on cluster 2, step 14 image features, step 15 CellPhoneDB export.

Usage (from the project root, the folder that contains SR007144_LL-2-1/ and SR007144_RL-1-5/)
-----
    python Manuscript/code/bone_analysis.py --root . --out Manuscript/code/results
    # optional: --geneset-sensitivity (Methods, three gene sets), --skip-moran (Supplementary Figure S2; slow),
    #           --roi-mask SAMPLE=PATH.geojson (hand-drawn ROI; repeat per sample), --project-8um
    python bone_analysis.py --config my_config.json --out results    # [APP] any configuration (pipeline_config.py)
    python -m webapp                                                  # [APP] the web application

Environment: ~/venvs/bone_spatial (scanpy 1.11.5, anndata 0.12.9, harmonypy 0.0.10, squidpy 1.8.1, igraph 1.0.0,
leidenalg 0.12.0, numpy 2.3.5, pandas 2.3.3, scipy 1.16.3, statsmodels 0.14.6, scikit-learn 1.8.0). [APP] harmonypy and
squidpy are optional: without harmonypy the graph is built on PCA, without squidpy Moran's I is computed here (same
numbers).
"""
from __future__ import annotations

import os

# [RUN] Reproducibility. pynndescent builds its trees in parallel and Harmony seeds its k-means internally, so
# random_state alone does not pin the result.
# [FIX] These variables take effect only when they are set before numpy, numba and the BLAS libraries load, i.e. when
# this file is the program being run. PYTHONHASHSEED (formerly also set here) is read at interpreter start and never
# applied; no output depends on string hashing. NUMBA_NUM_THREADS is left alone once numba is loaded: numba raises a
# RuntimeError on its next configuration reload when the value changes after its threads started (a notebook or test
# that imports scanpy first). run_analysis_from_config() limits the thread pools at run time instead, and the web
# application sets every variable in the environment of its worker process.
import sys as _sys

if os.environ.get("VISIUM_DETERMINISTIC", "1") == "1":
    for _v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        if not (_v == "NUMBA_NUM_THREADS" and "numba" in _sys.modules):
            os.environ.setdefault(_v, "1")
# [APP] numba (scanpy, pynndescent) caches compiled code next to the package sources; in a deep Windows folder without
# long-path support those file names exceed 260 characters and the import fails. A short per-user cache folder avoids
# it (the same folder the web application gives its worker).
if "numba" not in _sys.modules and "NUMBA_CACHE_DIR" not in os.environ:
    if os.name == "nt":
        _cache = os.path.join(os.environ.get("LOCALAPPDATA") or os.path.expanduser("~"), "BoneSpatialStudio", "numba")
    elif _sys.platform == "darwin":
        _cache = os.path.expanduser("~/Library/Caches/BoneSpatialStudio/numba")
    else:
        _cache = os.path.join(os.environ.get("XDG_CACHE_HOME") or os.path.expanduser("~/.cache"), "BoneSpatialStudio", "numba")
    os.environ["NUMBA_CACHE_DIR"] = _cache

import argparse
import contextlib
import copy
import json
import platform
import re
import sys
import traceback
import warnings
from datetime import datetime
from importlib import metadata
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pipeline_config as pc                                                    # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 2
RANDOM_SEED = 0

# ============================================================================
# CONSTANTS  [RUN]
# ============================================================================
# [APP] The gene signatures (Table 1), the signature lists, the receptor panel, the QC thresholds and the palette are
# defined once, in pipeline_config.DEFAULT_CONFIG, with their comments; the constants below are derived from it and
# describe the manuscript configuration.
SAMPLES = {
    "LL-2-1": {"condition": "control", "limb": "left", "label": "control"},
    "RL-1-5": {"condition": "tumor", "limb": "right", "label": "tumor-bearing"},
}
SECTION = {"LL-2-1": "control", "RL-1-5": "tumor-bearing"}
# microns per full-resolution pixel, `microns_per_pixel` in each section's scalefactors_json.json
# [FIX] Only a fallback now: the value is read from each section's scalefactors_json.json (microns_per_pixel()).
UM_PER_PX = {"LL-2-1": 0.5936042159338343, "RL-1-5": 0.5554917593191749}

_DEFAULT = pc.DEFAULT_CONFIG
# QC thresholds (Methods, "Quality control and filtering"). See qc_filter().
MIN_GENES, MIN_COUNTS, MAX_PCT_MT = (_DEFAULT["qc"][k] for k in ("min_genes", "min_counts", "max_pct_mt"))
# Gene signatures (Table 1): core, extended, competes, excluded (see pipeline_config._SIGNATURES).
SIG = {k: {f: copy.deepcopy(d[f]) for f in ("core", "extended", "competes", "excluded")}
       for k, d in _DEFAULT["signatures"].items()}
# signatures with a positivity call (Figures 4g, 6a; Methods)
PRIMARY = list(_DEFAULT["primary"])
# signatures competing in the deterministic tissue assignment (Figure 3)
PROGRAMS = list(_DEFAULT["programs"])
# Host receptor panel (Methods, "Host receptor panel and proximity analysis"; Figure 6b)
RECEPTOR_PANEL = dict(_DEFAULT["receptor_panel"])
ADIPOCYTE_GENES = list(_DEFAULT["adipocyte_genes"])                              # Supplementary Figure S4
# Transcripts abundant in essentially every cell; a domain headed by these has no interpretable identity (Methods)
UBIQUITOUS = set(_DEFAULT["ubiquitous"])
PALETTE = {"control": "#4C72B0", "tumor": "#C44E52", **{k: d["color"] for k, d in _DEFAULT["signatures"].items()}}

DE_COLUMNS = ["group", "names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
# [FIX] Figure 4h table (compare_osteoclast_states): scanpy's logfoldchanges replaced by a fold change of mean CP10k
COMPARTMENT_DE_COLUMNS = ["group", "names", "scores", "pvals", "pvals_adj", "log2fc_mean_cp10k", "pct_group",
                          "pct_rest", "selection_gene"]
COMPARTMENT_SUMMARY_COLUMNS = ["sample", "threshold", "threshold_is_zero", "pct_bins_above_threshold", "n_growth_plate",
                               "n_tumor_associated", "n_non_growth_plate_control", "n_growth_plate_outside_roi",
                               "n_tumor_associated_outside_roi", "median_umi_growth_plate",
                               "median_umi_tumor_associated"]
UNRESOLVED_COLUMNS = ["sample", "unresolved_bins", "all_bins", "percent", "unresolved_bins_roi", "all_bins_roi",
                      "percent_roi"]
GENESET_SENSITIVITY_COLUMNS = ["gene_set", "signature", "sample", "n_units_roi", "n_positive_native", "pct_native",
                               "ci_lo", "ci_hi", "pct_depth_matched_mean", "pct_depth_matched_sd", "genes_used",
                               "genes_missing"]
_NULL_PROGRESS = pc.Progress()


def _warn(progress, message):
    """[APP] Report a warning once per run (functions are called once per replicate)."""
    progress = progress or _NULL_PROGRESS
    if message not in progress.warnings:
        progress.warn(message)


def _as_bool(values) -> np.ndarray:
    """[APP] Boolean array of a column whose missing values count as False (NaN.astype(bool) would be True)."""
    s = pd.Series(values)
    if s.dtype == bool:
        return s.to_numpy()
    return np.asarray(s.astype("boolean").fillna(False), dtype=bool)


def _bin_label(bin_size: str) -> str:
    """[APP] '016um' -> '16um', as in the original file names."""
    return f"{int(pc.bin_um(bin_size))}um"


# ============================================================================
# 1. INPUT AND QUALITY CONTROL
# ============================================================================

def _norm_bc(x):
    """
    [RUN] Barcode string without the -1 suffix, so the matrix and the position table match.
    [FIX] A first expression also stripped a trailing "_<digits>_<digits>" ("bin suffix"). On a Visium HD barcode
    without the -1 suffix ("s_016um_00012_00034") that removed the bin row and column, so every barcode became
    "s_016um". Visium HD barcodes keep their coordinates; only the -<n> suffix is removed (barcodes ending in -1, as in
    the manuscript data, normalise exactly as before).
    """
    if pd.isna(x):
        return x
    return re.sub(r"-\d+$", "", str(x).strip())


def read_visium_hd(path, sample_name, condition=None, microns_per_pixel=None):
    """
    [SHAMIM, MODIFIED] step 1, "Load Visium spatial data".
    Original: sc.read_visium(path) and adata.var_names_make_unique().
    Changed because sc.read_visium reads standard Visium spots and does not read the Visium HD binned output
    (square_XXXum/, positions in tissue_positions.parquet). The matrix is read with sc.read_10x_h5 and joined to the
    parquet positions by barcode; bins without valid coordinates are dropped (Methods, "Data import and preprocessing").
    var_names_make_unique() is kept as in the original.
    [APP] condition and microns_per_pixel come from the configuration (defaults: SAMPLES, scalefactors_json.json).
    [FIX] Duplicate barcodes in the position table no longer break the join, and a join that places fewer than half of
    the matrix barcodes stops with an explanation instead of silently continuing on a fraction of the bins.
    """
    path, spatial = Path(path), Path(path) / "spatial"
    a = sc.read_10x_h5(path / "filtered_feature_bc_matrix.h5")
    a.var_names_make_unique()                                                   # [SHAMIM] step 1
    a.obs_names = [_norm_bc(x) for x in a.obs_names.astype(str)]
    a.obs_names_make_unique()
    tp = pd.read_parquet(spatial / "tissue_positions.parquet")
    bc = "barcode" if "barcode" in tp.columns else tp.columns[0]
    tp[bc] = tp[bc].astype(str).map(_norm_bc)
    tp = tp[tp[bc].isin(a.obs_names)].drop_duplicates(bc).set_index(bc).reindex(a.obs_names)
    coords = ("pxl_col_in_fullres", "pxl_row_in_fullres")
    for c in coords:
        tp[c] = pd.to_numeric(tp[c], errors="coerce")
    ok = tp[list(coords)].notna().all(axis=1)
    if a.n_obs and ok.mean() < 0.5:
        raise ValueError(f"{sample_name}: only {int(ok.sum()):,} of {a.n_obs:,} barcodes of "
                         f"{path / 'filtered_feature_bc_matrix.h5'} have a position in "
                         f"{spatial / 'tissue_positions.parquet'}; are both files from the same run and bin size?")
    a = a[ok.values, :].copy()
    tp = tp.loc[ok]
    a.obsm["spatial"] = tp[list(coords)].values                                 # full-resolution pixels
    for c in tp.columns:
        if c not in coords:
            a.obs[c] = tp[c].values                                             # in_tissue, array_row, array_col
    for c in ("in_tissue", "array_row", "array_col"):
        if c in a.obs:
            a.obs[c] = a.obs[c].astype(int)
    sf = spatial / "scalefactors_json.json"
    scalefactors = json.loads(sf.read_text()) if sf.exists() else {}
    if microns_per_pixel:
        scalefactors["microns_per_pixel"] = float(microns_per_pixel)
    a.uns["spatial"] = {sample_name: {"scalefactors": scalefactors}}
    a.obs["sample"] = sample_name
    a.obs["condition"] = condition or SAMPLES.get(sample_name, {}).get("condition", "unknown")
    return a


def microns_per_pixel(a, sample) -> float:
    """
    [APP] Microns per full-resolution pixel of a section: the value recorded at load time from its
    scalefactors_json.json (or the configuration's override), else the manuscript constant UM_PER_PX.
    [FIX] Replaces the hard-coded UM_PER_PX lookups, which raised KeyError for any other sample id.
    """
    v = a.uns.get("microns_per_pixel", {}).get(sample)
    if v is None:
        v = a.uns.get("spatial", {}).get(sample, {}).get("scalefactors", {}).get("microns_per_pixel")
    if v is None:
        v = UM_PER_PX.get(sample)
    if not v:
        raise ValueError(f"{sample}: microns_per_pixel is missing from scalefactors_json.json; set it in the "
                         f"sample's configuration")
    return float(v)


def qc_metrics(a):
    """
    [SHAMIM, MODIFIED] step 2, "Basic quality control".
    Original: sc.pp.calculate_qc_metrics(adata, inplace=True), which does not define mitochondrial genes, so the
    pct_counts_mt used in step 3 did not exist. Mitochondrial genes are now flagged by the mt- prefix, case-insensitive.
    """
    a.var["mt"] = a.var_names.str.lower().str.startswith("mt-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    return a


def qc_sweep(a, name, thresholds=(0, 5, 10, 25, 50, 100, 200)):
    """[RUN] Percentage of bins retained at each minimum gene count (Supplementary Figure S1, right panel)."""
    rows = []
    for t in thresholds:
        k = a.obs["n_genes_by_counts"] > t
        rows.append({"sample": name, "min_genes": t, "n_retained": int(k.sum()), "pct_retained": round(100 * k.mean(), 2),
                     "median_genes": float(a.obs.loc[k, "n_genes_by_counts"].median()) if k.any() else np.nan})
    return pd.DataFrame(rows)


def qc_filter(a, min_genes=MIN_GENES, min_counts=MIN_COUNTS, max_pct_mt=MAX_PCT_MT):
    """
    [SHAMIM, MODIFIED] step 3, "Filter low-quality spots".
    Original: more than 200 genes and less than 5% mitochondrial counts.
    Now: more than 10 genes, more than 15 UMI and less than 20% mitochondrial counts. The two sections differ about
    two-fold in capture, and a gene threshold carried over from single-cell work removes bins preferentially from the
    lower-capture control section (more than 100 genes kept 0.99% of control and 3.85% of tumor-bearing 8 um bins),
    which alone produced an apparent absence of osteoclast signal in the control (Methods; Supplementary Figure S1).
    [APP] Thresholds from the configuration (qc.*).
    """
    keep = ((a.obs["n_genes_by_counts"] > min_genes) & (a.obs["total_counts"] > min_counts)
            & (a.obs["pct_counts_mt"] < max_pct_mt))
    return a[keep.values].copy()


# ============================================================================
# 2. NORMALIZATION, EMBEDDING AND CLUSTERING (Leiden domains: Supplementary Figure S5 only)
# ============================================================================

def normalize_and_cluster(a, resolution=0.3, n_top_genes=2000, n_pcs=30, n_neighbors=15, harmony=True,
                          target_sum=1e4, scale_max_value=10, n_iterations=2, seed=RANDOM_SEED, cluster=True,
                          progress=None):
    """
    Normalization, highly variable genes, PCA, Harmony, neighbor graph, Leiden and UMAP on the joint object.
    Each step is tagged below.
    [FIX] No .raw is stored any more. `a.raw = a` froze the native log-normalized matrix, and scanpy's score_genes reads
    .raw whenever it exists; AnnData.copy() carries it, so every score computed on a copy holding a new matrix (depth
    matching, spatial pooling) silently came from the native matrix. layers["lognorm"] holds the same values.
    [APP] Parameters from the configuration (clustering.*). With one section: HVGs without batch_key and no Harmony.
    cluster=False stops after normalization (Leiden step disabled).
    """
    a.layers["counts"] = a.X.copy()                                             # [RUN] raw counts kept
    sc.pp.normalize_total(a, target_sum=target_sum)                             # [SHAMIM] step 4
    sc.pp.log1p(a)                                                              # [SHAMIM] step 4
    a.layers["lognorm"] = a.X.copy()                                            # [RUN] every score and test uses this
    if not cluster:
        return a
    n_sections = a.obs["sample"].nunique()

    # [SHAMIM, MODIFIED] step 5. Original: flavor="seurat", n_top_genes=2000, then adata = adata[:, highly_variable],
    # which discards every other gene. Now computed per section (batch_key="sample"), and only a separate copy is
    # subset and scaled, so that gene-set scoring and differential expression keep the full log-normalized matrix.
    sc.pp.highly_variable_genes(a, flavor="seurat", n_top_genes=n_top_genes,
                                batch_key="sample" if n_sections > 1 else None)
    tmp = a[:, a.var.highly_variable].copy()
    sc.pp.scale(tmp, max_value=scale_max_value)                                 # [SHAMIM] step 6
    n_comps = int(min(50, tmp.n_obs - 1, tmp.n_vars - 1))
    sc.tl.pca(tmp, n_comps=n_comps, svd_solver="arpack", random_state=seed)     # [SHAMIM] step 6 (seed added)
    a.obsm["X_pca"] = tmp.obsm["X_pca"]
    del tmp

    # [RUN] Harmony on the section label (residual between-section variation; Methods).
    rep = "X_pca"
    if harmony and n_sections > 1:
        try:
            import harmonypy
        except ImportError:
            _warn(progress, "harmonypy is not installed: the neighbor graph is built on PCA without Harmony")
        else:
            ho = harmonypy.run_harmony(a.obsm["X_pca"], a.obs, ["sample"], random_state=seed)
            Z = np.asarray(ho.Z_corr)
            if Z.shape[0] != a.n_obs:          # harmonypy 0.0.x returns (n_pcs, n_cells)
                Z = Z.T
            a.obsm["X_pca_harmony"] = Z
            rep = "X_pca_harmony"
    a.uns["neighbors_rep"] = rep

    # [SHAMIM, MODIFIED] step 7. Original: sc.pp.neighbors(n_neighbors=10, n_pcs=30); sc.tl.leiden(resolution=0.5).
    # Now 15 neighbors on the Harmony embedding (30 components as before), and Leiden with the modularity objective set
    # explicitly: flavor="igraph" defaults to the CPM objective, on whose scale resolution 0.5 produced 687 clusters of
    # median 23 bins. Resolution 0.3.
    sc.pp.neighbors(a, n_neighbors=n_neighbors, n_pcs=int(min(n_pcs, a.obsm[rep].shape[1])), use_rep=rep,
                    random_state=seed)
    sc.tl.leiden(a, resolution=resolution, key_added="leiden", flavor="igraph", objective_function="modularity",
                 n_iterations=n_iterations, random_state=seed, directed=False)
    sc.tl.umap(a, random_state=seed)                                            # [SHAMIM] step 8 (embedding only)
    return a


def _testable_groups(labels, what, progress=None):
    """[APP] Groups with at least two bins (rank_genes_groups cannot test smaller ones)."""
    counts = pd.Series(labels).astype(str).value_counts()
    small = counts[counts < 2]
    if len(small):
        _warn(progress, f"{what}: {len(small)} group(s) with fewer than two bins cannot be tested and are left out "
                        f"({', '.join(map(str, small.index[:10]))})")
    return [str(g) for g in counts[counts >= 2].index]


def domain_markers(a, key="leiden", n_top=30, progress=None):
    """
    [SHAMIM, MODIFIED] steps 10-11, "Differential expression between clusters" and "Export DEGs".
    Original: sc.tl.rank_genes_groups(groupby='leiden', method='wilcoxon') on the HVG-subset, scaled matrix.
    Now on the full log-normalized matrix (use_raw=False, layer copied into X), with detection fractions (pts=True).
    Benjamini-Hochberg adjustment is scanpy's default, as before.
    [FIX] A group of fewer than two bins made rank_genes_groups raise; such groups are left out with a warning.
    [APP] Reads layers["lognorm"] directly instead of copying it into X.
    """
    groups = _testable_groups(a.obs[key], f"domain markers ({key})", progress)
    if len(groups) < 2:
        _warn(progress, f"domain markers ({key}): fewer than two testable groups, no markers")
        return pd.DataFrame(columns=DE_COLUMNS + ["pct_nz_group", "pct_nz_reference"])
    all_testable = len(groups) == a.obs[key].astype(str).nunique()
    sc.tl.rank_genes_groups(a, groupby=key, groups="all" if all_testable else groups, method="wilcoxon",
                            use_raw=False, layer="lognorm", pts=True)
    df = sc.get.rank_genes_groups_df(a, group=None)
    if "group" not in df:                                                       # a single tested group
        df.insert(0, "group", groups[0])
    return (df.sort_values(["group", "scores"], ascending=[True, False]).groupby("group").head(n_top)
              .reset_index(drop=True))


def domain_screen(a, de, key="leiden", n_top=10, ubiquitous=UBIQUITOUS, genes_per_umi=0.94, ubiquitous_flag=0.2):
    """
    [RUN] Per Leiden domain: the fraction of its top markers that are ubiquitously expressed transcripts (the only
    merging criterion), and median detected genes per UMI. A domain is flagged when shallow sampling (>= 0.94 genes per
    UMI) combines with non-specific markers; flagged domains are not merged (Methods).
    """
    columns = ["domain", "n_bins", "ubiquitous_fraction", "median_genes_per_umi", "top_markers", "de_uninterpretable"]
    if de.empty:
        return pd.DataFrame(columns=columns)
    cplx = a.obs["n_genes_by_counts"] / a.obs["total_counts"].clip(lower=1)
    med_cplx = cplx.groupby(a.obs[key].astype(str), observed=True).median()
    rows = []
    for g, sub in de.groupby("group", observed=True):
        top = sub.nlargest(n_top, "scores")["names"].tolist()
        rows.append({"domain": g, "n_bins": int((a.obs[key].astype(str) == str(g)).sum()),
                     "ubiquitous_fraction": round(sum(t in ubiquitous for t in top) / max(len(top), 1), 2),
                     "median_genes_per_umi": round(float(med_cplx.get(str(g), np.nan)), 3),
                     "top_markers": ",".join(top[:5])})
    df = pd.DataFrame(rows, columns=columns[:-1])
    df["de_uninterpretable"] = (df["median_genes_per_umi"] >= genes_per_umi) & (df["ubiquitous_fraction"] >= ubiquitous_flag)
    return df


def collapse_unresolved(a, screen, key="leiden", ubiq_frac=0.3, label="unresolved"):
    """
    [RUN] Merge domains whose top markers are ubiquitously expressed transcripts into 'unresolved' (Methods;
    Supplementary Figure S5: 50,273 bins, 33.1%). 0.3 is the natural break in this dataset (merged domains 0.3-0.6,
    next value 0.1). The mechanism is not established.
    [FIX] The note ended "and 'unresolved' is not a tumor label". The label is defined without reference to the tumor,
    but it is not tumor-free: the mouse probes cross-capture A549 RNA (module docstring), and in the tumor-bearing
    section the unresolved domain is strongly enriched for the lesion (cross-captured Krt7 and Krt18; Tmsb10 and Itm2b
    track the lesion). Its share also
    differs between sections (this code on the SR007144 data: 18.3% of the control and 42.3% of the tumor-bearing bins,
    33.3% pooled; leiden_signature_matrix reports it per section).
    """
    hits = screen.loc[screen["ubiquitous_fraction"] >= ubiq_frac, "domain"].astype(str).tolist() if len(screen) else []
    new = a.obs[key].astype(str)
    a.obs[key] = pd.Categorical(new.where(~new.isin(hits), label))
    print(f"[domains] merged {len(hits)} domains ({int(new.isin(hits).sum()):,} bins) into '{label}'")
    return a


def collapse_fragments(a, key="leiden", min_bins=100, label="unassigned"):
    """[RUN] Domains of fewer than 100 bins (kNN-graph debris) relabeled 'unassigned' (Methods)."""
    vc = a.obs[key].astype(str).value_counts()
    small = {str(x) for x in vc[vc < min_bins].index}
    if small:
        a.obs[key] = pd.Categorical(a.obs[key].astype(str).where(~a.obs[key].astype(str).isin(small), label))
    return a


def clustering_stability(a, resolution=0.3, seeds=(0, 1, 2, 3, 4), n_neighbors=15, n_pcs=30, n_iterations=2,
                         use_rep=None):
    """
    [RUN] Adjusted Rand index between Leiden runs that differ only by seed (graph and partition re-seeded together).
    Methods and Supplementary Figure S5: "across five seeds the adjusted Rand index was 0.72" (median of the pairs).
    [FIX] The cleanup now also removes the parameters sc.tl.leiden stores in uns["_s<seed>"] (they were left in the
    h5ad).
    """
    from sklearn.metrics import adjusted_rand_score
    use_rep = use_rep or a.uns.get("neighbors_rep", "X_pca_harmony")
    seeds = list(dict.fromkeys(int(s) for s in seeds))                          # a repeated seed would compare itself
    use_pcs = int(min(n_pcs, a.obsm[use_rep].shape[1]))
    labs = {}
    for sd in seeds:
        sc.pp.neighbors(a, n_neighbors=n_neighbors, n_pcs=use_pcs, use_rep=use_rep, random_state=sd,
                        key_added=f"n{sd}")
        sc.tl.leiden(a, resolution=resolution, key_added=f"_s{sd}", flavor="igraph", objective_function="modularity",
                     n_iterations=n_iterations, random_state=sd, directed=False, neighbors_key=f"n{sd}")
        labs[sd] = a.obs[f"_s{sd}"].astype(str).values
    rows = [{"seed_a": i, "seed_b": j, "n_clusters_a": len(set(labs[i])), "n_clusters_b": len(set(labs[j])),
             "ARI": round(adjusted_rand_score(labs[i], labs[j]), 4)}
            for k, i in enumerate(seeds) for j in seeds[k + 1:]]
    for sd in seeds:
        del a.obs[f"_s{sd}"], a.uns[f"n{sd}"], a.obsp[f"n{sd}_distances"], a.obsp[f"n{sd}_connectivities"]
        a.uns.pop(f"_s{sd}", None)
    df = pd.DataFrame(rows, columns=["seed_a", "seed_b", "n_clusters_a", "n_clusters_b", "ARI"])
    if len(df):
        print(f"[stability] median ARI {df['ARI'].median():.3f}")
    return df


def map_domains_from_coarser(a_fine, coarse, fine_um, coarse_um, keys=("leiden", "program"), missing="no_parent_bin"):
    """
    [RUN] Project 16 um labels onto 8 um bins by array geometry alone (Methods: at 8 um many bins share identical PCA
    coordinates and the graph fragments, so clustering is done at 16 um). Space Ranger bins nest exactly: an 8 um array
    index maps to the 16 um index by integer division by 2.
    [FIX] Fine bins whose parent was removed at QC are now 'no_parent_bin'; they used to be 'unassigned', the label
    collapse_fragments gives to bins of fragment domains, so the two could not be told apart.
    """
    ratio = coarse_um / fine_um
    if ratio < 1 or ratio != int(ratio):
        raise ValueError(f"cannot project {coarse_um:g} um labels onto {fine_um:g} um bins: the fine bin size must "
                         f"divide the coarse one")
    f = int(ratio)
    for key in keys:
        if key not in coarse.obs:
            continue
        lut = {(s, int(r), int(c)): lab for s, r, c, lab in zip(coarse.obs["sample"].astype(str), coarse.obs["array_row"],
                                                                 coarse.obs["array_col"], coarse.obs[key].astype(str))}
        a_fine.obs[key] = pd.Categorical([lut.get((s, int(r) // f, int(c) // f), missing) for s, r, c in
                                          zip(a_fine.obs["sample"].astype(str), a_fine.obs["array_row"], a_fine.obs["array_col"])])
    return a_fine


# ============================================================================
# 3. GENE-SET SCORING AND POSITIVITY  [RUN]
# ============================================================================

def present(a, genes):
    return [g for g in genes if g in a.var_names]


def score_all(a, sigs=SIG, ctrl_size=50, seed=RANDOM_SEED):
    """
    [RUN] Per-bin score of every signature with scanpy score_genes on the log-normalized matrix, and the number of
    core genes detected per bin.
    [FIX] Scores are computed from layers["lognorm"] explicitly. Without layer/use_raw, score_genes read .raw, which
    held the native matrix, so the scores of every depth-matched and spatially pooled copy were the native scores (and
    only the core-gene counts changed).
    [FIX] Wording: the background is not "50 genes". score_genes (n_bins=25) draws ctrl_size control genes from EACH of
    the 25 average-expression bins that contains a signature gene (signature genes excluded), so it holds up to
    ctrl_size x (number of such bins) genes.
    """
    lognorm = a.layers["lognorm"]
    cov = []
    for name, d in sigs.items():
        ext, core = present(a, d["extended"]), present(a, d["core"])
        cov.append({"signature": name, "extended_detected": ",".join(ext), "core_detected": ",".join(core),
                    "missing": ",".join(g for g in d["extended"] if g not in a.var_names)})
        if len(ext) >= 2:
            sc.tl.score_genes(a, ext, score_name=f"score_{name}", ctrl_size=ctrl_size, random_state=seed,
                              use_raw=False, layer="lognorm")
        else:
            a.obs[f"score_{name}"] = np.nan
        if core:
            X = lognorm[:, a.var_names.get_indexer(core)]
            a.obs[f"ncore_{name}"] = np.asarray((X > 0).sum(axis=1)).ravel()
        else:
            a.obs[f"ncore_{name}"] = 0
    a.uns["signature_coverage"] = pd.DataFrame(cov)
    return a


def null_threshold(a, name, sigs=SIG, pct=95.0, n_perm=50, n_sub=20000, ctrl_size=50, seed=RANDOM_SEED):
    """
    [RUN] 95th percentile of a null of randomly drawn gene sets of the same size, scored the same way on 20,000
    sampled bins. Bins whose null score is exactly 0 (no gene of either set detected) are excluded, otherwise the
    threshold collapses to 0.
    [FIX] The null is scored from layers["lognorm"] (it read the native .raw, see score_all, so the thresholds of every
    depth-matched and pooled call were native-depth thresholds), on a light object instead of a full copy. The random
    draws are unchanged, so native thresholds equal the manuscript's.
    """
    genes = present(a, sigs[name]["extended"])
    if len(genes) < 2:
        return np.nan
    rng = np.random.default_rng(seed)
    idx = rng.choice(a.n_obs, min(n_sub, a.n_obs), replace=False)
    sub = ad.AnnData(X=a.layers["lognorm"][idx], var=pd.DataFrame(index=a.var_names.copy()))
    nulls = []
    for _ in range(n_perm):
        rg = list(rng.choice(sub.var_names, len(genes), replace=False))
        sc.tl.score_genes(sub, rg, score_name="_n", ctrl_size=ctrl_size, random_state=int(rng.integers(1e6)),
                          use_raw=False)
        v = sub.obs["_n"].values
        nulls.append(v[v != 0])
    nulls = np.concatenate(nulls) if nulls else np.array([])
    return float(np.percentile(nulls, pct)) if nulls.size else 0.0


def call_positive(a, names=PRIMARY, sigs=SIG, min_core=2, n_perm=50, pct=95.0, n_sub=20000, ctrl_size=50,
                  seed=RANDOM_SEED, progress=None):
    """
    [RUN] A bin is positive for a signature when its score exceeds the null threshold, at least two core genes are
    detected, and no designated competing signature scores higher (Methods; Table 1).
    [FIX] A competitor vetoes a bin only where it scores strictly higher. `score >= competitor` is False whenever the
    competitor cannot be scored (NaN), which silently removed every positive bin of the signature.
    [FIX] A signature that cannot be scored gets an all-False pos_ column and a NaN threshold, with a warning. It used to
    get no column: copies kept a stale pos_ column from the previous run (the gene-set sensitivity then reported the
    main calls for a variant whose genes were absent), and later steps indexing pos_<name> crashed.
    """
    thr = {}
    for n in names:
        key = f"score_{n}"
        if key not in a.obs or a.obs[key].isna().all():
            _warn(progress, f"signature {n}: fewer than two of its scored genes are in the data, so no bin is "
                            f"called positive")
            thr[n] = np.nan
            a.obs[f"pos_{n}"] = False
            continue
        t = null_threshold(a, n, sigs=sigs, pct=pct, n_perm=n_perm, n_sub=n_sub, ctrl_size=ctrl_size, seed=seed)
        thr[n] = t
        score = a.obs[key]
        base = (score > t) & (a.obs[f"ncore_{n}"] >= min_core)
        for comp in sigs[n]["competes"]:
            if f"score_{comp}" in a.obs:
                base &= ~(a.obs[f"score_{comp}"] > score)
        a.obs[f"pos_{n}"] = base.fillna(False).astype(bool)
    a.uns["positivity_thresholds"] = thr
    return a


def positivity_table(a, names=PRIMARY, roi_key=None, reference=None, samples=None, tumor_samples=None):
    """
    [RUN] Positive bins per signature and section, with Wilson 95% intervals (Figures 4g and 6a, left panels).
    [FIX] The section columns were hard-coded (w["RL-1-5"] / w["LL-2-1"]): any other sample id, or a section without
    ROI bins, raised KeyError (in every depth-matched replicate too). Every row now carries fold_vs_reference (its
    percentage over the reference section's); fold_tumor_over_control is kept, computed as before (tumor over control
    percentage, rounded to 2 decimals), when there are exactly two sections and one of them is the tumor section
    (tumor_samples; default: the non-reference one). A section without bins gives NaN instead of a division by zero.
    [FIX] Both folds are computed from the exact rates n_positive / n_bins; only the stored folds are rounded
    (fold_vs_reference 4 decimals, fold_tumor_over_control 2). They were ratios of pct_positive already rounded to 3
    decimals, which moved small rates: on the SR007144 data the osteoclast fold_tumor_over_control was 169.36 instead of
    171.30 (the Figure 4g caption, from exact counts), and the megakaryocyte fold_vs_reference was off by 26%.
    """
    from statsmodels.stats.proportion import proportion_confint
    obs = a.obs if roi_key is None else a.obs[_as_bool(a.obs[roi_key])]
    samples = list(samples) if samples else list(pd.unique(a.obs["sample"].astype(str)))
    smp = obs["sample"].astype(str)
    rows, exact = [], []
    for n in names:
        if f"pos_{n}" not in obs:
            continue
        for s in samples:
            sub = obs[smp == s]
            x, tot = int(_as_bool(sub[f"pos_{n}"]).sum()), len(sub)
            if tot:
                lo, hi = proportion_confint(x, tot, method="wilson")
                rate, lo, hi = 100 * x / tot, 100 * lo, 100 * hi
            else:
                rate = lo = hi = np.nan
            exact.append(rate)
            rows.append({"signature": n, "sample": s, "roi": roi_key or "whole_section", "n_positive": x, "n_bins": tot,
                         "pct_positive": round(rate, 3), "ci_lo": round(lo, 3), "ci_hi": round(hi, 3),
                         "mean_score": float(sub[f"score_{n}"].mean()) if tot else np.nan})
    df = pd.DataFrame(rows, columns=["signature", "sample", "roi", "n_positive", "n_bins", "pct_positive", "ci_lo",
                                     "ci_hi", "mean_score"])
    if df.empty or not samples:
        return df
    rate = pd.Series(exact, index=pd.MultiIndex.from_arrays([df["signature"], df["sample"]]), dtype=float)
    reference = reference if reference in samples else samples[0]
    ref = rate.xs(reference, level=1).replace(0, np.nan)
    df["fold_vs_reference"] = (rate.to_numpy() / df["signature"].map(ref).to_numpy(float)).round(4)
    tumor = [s for s in samples if s in tumor_samples] if tumor_samples is not None else \
        [s for s in samples if s != reference]
    if len(samples) == 2 and len(tumor) == 1:                                   # not quoted in the text
        control = next(s for s in samples if s != tumor[0])
        fold = (rate.xs(tumor[0], level=1) / rate.xs(control, level=1).replace(0, np.nan)).round(2)
        df["fold_tumor_over_control"] = df["signature"].map(fold)
    return df


def depth_match_target(a):
    """[RUN] Median UMI per bin of the shallower section, the common depth of every between-section comparison."""
    med = a.obs.groupby("sample", observed=True)["total_counts"].median()
    return int(med.min()), med


def depth_match_targets(a, mode="cap"):
    """
    [APP] Per-bin UMI targets of the depth matching (configuration scoring.depth_matching).
      cap       (manuscript) every bin above the median UMI per bin of the shallowest section is downsampled to that
                median, in every section including the shallowest; bins below it keep their counts. Deeper sections
                keep more bins at the target than the shallowest one, a residual depth advantage (CHANGES.md).
      quantile  the shallowest section keeps its counts; each bin of every other section is downsampled to the UMI count
                at the same quantile of the shallowest section's per-bin distribution (never above its own total), so the
                per-bin depth distributions match.
    [FIX] Wording. cap lowers both sections (SR007144: the control ROI mean falls from 81.9 to 47.5 UMI per bin), and
    the residual advantage of the deeper section is small (+4.8%). quantile matches whole-section distributions, which
    leaves the tumor-bearing ROI about 15% shallower than the control ROI: conservative for tumor-enriched programs.
    Returns (targets per bin, median target, per-section medians).
    """
    target, med = depth_match_target(a)
    totals = a.obs["total_counts"].to_numpy(dtype=float)
    if mode == "cap":
        return np.full(a.n_obs, target, dtype=np.int64), target, med
    if mode != "quantile":
        raise ValueError(f"depth matching must be cap or quantile, got {mode!r}")
    smp = a.obs["sample"].astype(str).to_numpy()
    shallow = str(med.idxmin())
    ref = np.sort(totals[smp == shallow])
    targets = totals.copy()
    for s in pd.unique(smp):
        if s == shallow:
            continue
        m = np.flatnonzero(smp == s)
        q = (pd.Series(totals[m]).rank(method="first").to_numpy() - 0.5) / len(m)
        targets[m] = np.minimum(totals[m], np.floor(np.quantile(ref, q)))
    return targets.astype(np.int64), target, med


def downsample_counts_matrix(counts, targets, seed):
    """
    [FIX] Depth matching by an independent subsample without replacement of each bin's UMIs (one multivariate
    hypergeometric draw per bin from a single numpy Generator per replicate). sc.pp.downsample_counts re-seeds numba's
    generator with the same seed for every bin, so every bin with the same UMI total kept the same UMI positions: one
    draw shared by the section rather than an independent draw per bin. Bins at or below their target are unchanged.
    """
    X = sp.csr_matrix(counts, dtype=np.float32, copy=True)
    totals = np.asarray(X.sum(axis=1)).ravel()
    targets = np.broadcast_to(np.asarray(targets), totals.shape)
    rng = np.random.default_rng(seed)
    data = np.rint(X.data).astype(np.int64)
    for r in np.flatnonzero(totals > targets):
        s, e = X.indptr[r], X.indptr[r + 1]
        data[s:e] = rng.multivariate_hypergeometric(data[s:e], int(targets[r]))
    X.data = data.astype(np.float32)
    X.eliminate_zeros()
    return X


def downsampled_counts(a, target, rep):
    """[RUN] One depth-matching replicate: raw counts downsampled to at most `target` UMI per bin, seed 0 + rep.
    [FIX] Independent per-bin subsampling (downsample_counts_matrix) instead of sc.pp.downsample_counts."""
    return downsample_counts_matrix(a.layers["counts"], target, RANDOM_SEED + rep)


def _lognorm_object(counts, a, obs_cols, target_sum=1e4):
    """[APP] Light AnnData (no .raw, no other layers) with log1p(normalize_total(counts)) in X and layers["lognorm"]."""
    b = ad.AnnData(X=sp.csr_matrix(counts, dtype=np.float32), obs=a.obs[[c for c in obs_cols if c in a.obs]].copy(),
                   var=pd.DataFrame(index=a.var_names.copy()))
    sc.pp.normalize_total(b, target_sum=target_sum)
    sc.pp.log1p(b)
    b.layers["lognorm"] = b.X
    return b


def depth_matched(a, names=PRIMARY, n_reps=5, roi_key=None, sigs=SIG, mode="cap", min_core=2, n_perm=20, pct=95.0,
                  n_sub=20000, ctrl_size=50, seed=RANDOM_SEED, target_sum=1e4, samples=None, progress=None):
    """
    [RUN] Positivity after depth matching, mean and SD over five replicates (Figures 4g and 6a, right panels).
    [FIX] Each replicate is scored, and its null thresholds computed, from its own downsampled matrix; both used to come
    from the native .raw (see score_all), so only the core-gene count reflected the downsampling. Replicates use
    independent per-bin subsampling (downsample_counts_matrix) and a light object instead of a full copy. The table
    gains n_positive_mean and n_bins, so figures label bars with counts instead of rebuilding them from rounded
    percentages, and the means are taken over unrounded percentages.
    [FIX] The table also carries ci_lo, ci_hi: the Wilson 95% interval (percent) of round(n_positive_mean) positive bins
    out of n_bins, a bin-level interval like the native table's. pct_positive_sd is the spread over the downsampling
    replicates only; it says nothing about the sampling uncertainty of the rate.
    """
    from statsmodels.stats.proportion import proportion_confint
    targets, target, _ = depth_match_targets(a, mode)
    print(f"[depth-match] {mode}: downsampling to {target} counts/bin")
    out = []
    for rep in range(n_reps):
        X = downsample_counts_matrix(a.layers["counts"], targets, seed + rep)
        b = _lognorm_object(X, a, ["sample", roi_key] if roi_key else ["sample"], target_sum)
        b = call_positive(score_all(b, sigs, ctrl_size, seed), names, sigs=sigs, min_core=min_core, n_perm=n_perm,
                          pct=pct, n_sub=n_sub, ctrl_size=ctrl_size, seed=seed, progress=progress)
        t = positivity_table(b, names, roi_key=roi_key, samples=samples)
        t["rep"] = rep
        out.append(t)
    df = pd.concat(out, ignore_index=True)
    df["pct_exact"] = 100 * df["n_positive"] / df["n_bins"].where(df["n_bins"] > 0)
    g = (df.groupby(["signature", "sample"], as_index=False, sort=False)
           .agg(pct_positive_mean=("pct_exact", "mean"), pct_positive_sd=("pct_exact", "std"),
                n_positive_mean=("n_positive", "mean"), n_bins=("n_bins", "first"), mean_score=("mean_score", "mean")))
    ci = [proportion_confint(int(round(k)), int(n), method="wilson") if n > 0 else (np.nan, np.nan)
          for k, n in zip(g["n_positive_mean"], g["n_bins"])]
    g["ci_lo"] = [100 * lo for lo, _ in ci]
    g["ci_hi"] = [100 * hi for _, hi in ci]
    return g.round({"pct_positive_mean": 4, "pct_positive_sd": 4, "n_positive_mean": 2, "ci_lo": 4, "ci_hi": 4})


# ============================================================================
# 4. REGION OF INTEREST  [RUN]
# ============================================================================

def define_marrow_roi(a, exclude=("skeletal_muscle", "nerve", "smooth_muscle"), drawn=None, pct=90.0, progress=None):
    """
    [RUN] Marrow and bone region of interest. Transcriptional part: skeletal muscle, peripheral nerve and vascular
    smooth muscle have no positivity call, so a bin is excluded when its score for any of the three is above that
    signature's 90th percentile over all bins of both sections. When the hand-drawn ROI exists
    (bone_pending_masks.rasterise_drawn_roi), it is intersected with the transcriptional part and is then the only thing
    that excludes the diaphysis. All results in the current manuscript use the transcriptional ROI (drawn=None).
    [FIX] An exclusion signature that cannot be scored (all NaN) made the percentile NaN, and `score <= NaN` is False
    for every bin, so the region of interest came out empty. Such a signature is now skipped with a warning, and only a
    score above the percentile excludes a bin (a NaN score never does).
    """
    keep = np.ones(a.n_obs, dtype=bool)
    for e in exclude:
        if f"pos_{e}" in a.obs:
            keep &= ~_as_bool(a.obs[f"pos_{e}"])
            continue
        score = a.obs[f"score_{e}"].to_numpy(dtype=float) if f"score_{e}" in a.obs else np.full(a.n_obs, np.nan)
        if np.isnan(score).all():
            _warn(progress, f"region of interest: signature {e} cannot be scored, so it excludes no bin")
            continue
        keep &= ~(score > np.nanpercentile(score, pct))
    a.obs["in_roi_transcriptional"] = keep
    if drawn is not None:
        a.obs["in_roi_drawn"] = np.asarray(drawn, dtype=bool)
        keep = keep & a.obs["in_roi_drawn"].values
    a.obs["in_marrow_roi"] = keep
    return a


# ============================================================================
# 5. DETERMINISTIC TISSUE ASSIGNMENT (Figure 3)  [RUN]
# ============================================================================

def spatial_neighbor_matrix(a, k):
    """[APP] Sparse n_obs x n_obs 0/1 matrix: row i marks bin i and its k nearest bins of the same section (pixel
    coordinates), the pooling of spatially_pooled_scores. A section with at most k bins pools all of its bins."""
    from scipy.spatial import cKDTree
    smp = a.obs["sample"].astype(str).to_numpy()
    xy = np.asarray(a.obsm["spatial"], dtype=float)
    rows, cols = [], []
    for s in pd.unique(smp):
        m = np.flatnonzero(smp == s)
        kk = min(k + 1, len(m))
        _, idx = cKDTree(xy[m]).query(xy[m], k=kk)
        idx = np.asarray(idx).reshape(len(m), kk)
        rows.append(np.repeat(m, kk))
        cols.append(m[idx.ravel()])
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    return sp.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(a.n_obs, a.n_obs))


def spatially_pooled_scores(a, pool_k=8, depth_match=False, mode="cap", sigs=SIG, names=PRIMARY, min_core=2,
                            n_perm=20, pct=95.0, n_sub=20000, ctrl_size=50, seed=RANDOM_SEED, target_sum=1e4,
                            progress=None):
    """
    [RUN] Raw counts summed over each bin and its eight nearest spatial neighbors before normalization and scoring
    (about ninefold depth, still deterministic). depth_match=True first downsamples every bin to the median of the
    shallower section. Writes pooled_* (native depth) or pooled_dm_* (depth-matched) columns.
    [FIX] The pooled scores and their null thresholds are now computed from the pooled matrix. They came from the
    native, unpooled .raw (see score_all), so 'program' and 'program_dm' (Figure 3), the pooled calls of Figure 4i and
    the fibroblast proximity table were driven by native per-bin scores, with only the core-gene count pooled.
    [FIX] Pooling is one sparse matrix product. Assigning `A @ X[m]` into a lil_matrix by fancy row index densified each
    section (tens of GB on the real data), and a section with fewer than pool_k + 1 bins crashed the neighbor query.
    """
    counts = a.layers["counts"]
    if depth_match:
        targets, _, _ = depth_match_targets(a, mode)
        counts = downsample_counts_matrix(counts, targets, seed)
    pooled = (spatial_neighbor_matrix(a, pool_k) @ sp.csr_matrix(counts)).tocsr()
    b = _lognorm_object(pooled, a, ["sample"], target_sum)
    b = call_positive(score_all(b, sigs, ctrl_size, seed), names, sigs=sigs, min_core=min_core, n_perm=n_perm, pct=pct,
                      n_sub=n_sub, ctrl_size=ctrl_size, seed=seed, progress=progress)
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
    names = np.array([p for p in programs if f"{prefix}score_{p}" in a.obs and f"{prefix}ncore_{p}" in a.obs])
    if not len(names):
        a.obs[label_col] = pd.Categorical([ambiguous] * a.n_obs)
        return a
    S = a.obs[[f"{prefix}score_{n}" for n in names]].to_numpy(dtype=float)
    thr_src = a.uns.get(f"{prefix}thresholds", {})
    thr = np.array([thr_src.get(n, 0.0) for n in names], dtype=float)
    core = np.array([a.obs[f"{prefix}ncore_{n}"].to_numpy() for n in names]).T
    ok = (S > thr) & (core >= min_core)
    best = np.argmax(np.where(ok, S, -np.inf), axis=1)
    a.obs[label_col] = pd.Categorical(np.where(ok.any(axis=1), names[best], ambiguous))
    return a


def program_composition(a, col, samples=None):
    """[RUN] Bins per assigned signature and section, with both denominators and the ambiguous rate (Figure 3b)."""
    ct = pd.crosstab(a.obs["sample"].astype(str), a.obs[col].astype(str))
    if samples:
        ct = ct.reindex([s for s in samples if s in ct.index])
    rows = []
    for smp, r in ct.iterrows():
        tot = r.sum()
        assigned = tot - r.get("ambiguous", 0)
        for prog, n in r.items():
            rows.append({"sample": smp, "program": prog, "n_bins": int(n), "pct_of_all_bins": round(100 * n / tot, 3),
                         "pct_of_assigned_bins": round(100 * n / assigned, 3) if prog != "ambiguous" and assigned else np.nan,
                         "sample_ambiguous_pct": round(100 * r.get("ambiguous", 0) / tot, 2)})
    return pd.DataFrame(rows, columns=["sample", "program", "n_bins", "pct_of_all_bins", "pct_of_assigned_bins",
                                       "sample_ambiguous_pct"])


# ============================================================================
# 6. OSTEOCLAST COMPARTMENTS (Figure 4h)  [RUN, with SHAMIM step 10 for the test]
# ============================================================================
COMPARTMENTS = ["negative", "growth_plate", "tumor_associated", "non_growth_plate_control"]


def split_osteoclast_compartments(a, pct=90, signature="osteoclast", split_by="chondrocyte",
                                  tumor_conditions=("tumor",), progress=None):
    """
    [RUN] Osteoclast-positive bins split by chondrocyte score: above the 90th percentile of all bins = growth-plate
    (chondro-osseous junction) compartment, the rest = tumor-associated (Methods; Figure 4h).
    [FIX] Only sections whose condition is a tumor condition (compartments.tumor_conditions) hold tumor-associated bins.
    Non-growth-plate osteoclast bins of the contralateral control limb, which has no tumor, were labelled
    'tumor_associated' and pooled into the Figure 4h test; they are now 'non_growth_plate_control' and not tested.
    [FIX] The percentile value is stored in a.uns["compartment_threshold"] (NaN when the score cannot be computed), and
    a warning says when it is 0. The chondrocyte score is zero-inflated (44.5% of the SR007144 bins score exactly 0),
    so its 90th percentile is 0 and 'above the 90th percentile' means a score above 0, which requires at least one UMI
    of a chondrocyte gene (2.1% of bins), not the top 10% of bins.
    """
    split = a.obs[f"score_{split_by}"].to_numpy(dtype=float)
    thr = float(np.nanpercentile(split, pct)) if not np.isnan(split).all() else np.nan
    gp = split > thr if np.isfinite(thr) else np.zeros(a.n_obs, dtype=bool)
    a.uns["compartment_threshold"] = thr
    if thr == 0:
        _warn(progress, f"osteoclast compartments: the {pct:g}th percentile of the {split_by} score is 0 "
                        f"({100 * float(np.mean(split == 0)):.1f}% of bins score exactly 0), so the growth-plate rule "
                        f"is a {split_by} score above 0, i.e. at least one UMI of a {split_by} gene "
                        f"({100 * float(np.mean(gp)):.1f}% of bins)")
    a.obs["growth_plate"] = gp
    oc = _as_bool(a.obs[f"pos_{signature}"])
    tumor = a.obs["condition"].astype(str).isin(list(tumor_conditions)).to_numpy()
    comp = np.full(a.n_obs, "negative", dtype=object)
    comp[oc & gp] = "growth_plate"
    comp[oc & ~gp & tumor] = "tumor_associated"
    comp[oc & ~gp & ~tumor] = "non_growth_plate_control"
    a.obs["oc_compartment"] = pd.Categorical(comp, categories=COMPARTMENTS)
    return a


def compare_osteoclast_states(a, sections="all", tumor_conditions=("tumor",), selection_genes=(), progress=None):
    """
    [SHAMIM, MODIFIED] step 10 (Wilcoxon rank-sum, BH-adjusted), applied to the two osteoclast compartments
    instead of Leiden clusters, on the full log-normalized matrix.
    [FIX] Each compartment needs at least two bins: with 0 or 1 bins rank_genes_groups raised or returned all-zero
    scores; the test is now skipped with a warning and an empty table.
    [APP] sections="tumor" takes both compartments from tumor sections only (no between-section depth difference).
    [FIX] scanpy's logfoldchanges is dropped: it is log2 of the ratio of expm1(mean log-expression) with a 1e-9 floor,
    so a gene undetected in one compartment got |values| of 19-31.6 (26.9% of the rows on the SR007144 data). In its
    place: log2fc_mean_cp10k = log2((mean CP10k in the group + 1) / (mean CP10k in the other compartment + 1)), CP10k =
    expm1(lognorm); pct_group and pct_rest, the percentage of bins detecting the gene in each; and selection_gene, True
    for the genes that define the split (selection_genes: the split_by signature's scored genes), which are higher in
    the growth-plate compartment by construction. Scores and P values are unchanged. Columns: COMPARTMENT_DE_COLUMNS.
    """
    comp = a.obs["oc_compartment"].astype(str)
    mask = comp.isin(["growth_plate", "tumor_associated"]).to_numpy()
    if sections == "tumor":
        mask &= a.obs["condition"].astype(str).isin(list(tumor_conditions)).to_numpy()
    n = comp[mask].value_counts()
    if n.get("growth_plate", 0) < 2 or n.get("tumor_associated", 0) < 2:
        _warn(progress, f"osteoclast compartments: growth-plate {int(n.get('growth_plate', 0))} and tumor-associated "
                        f"{int(n.get('tumor_associated', 0))} bins; at least two of each are needed, no test")
        return pd.DataFrame(columns=COMPARTMENT_DE_COLUMNS)
    X = sp.csr_matrix(a.layers["lognorm"][mask])
    labels = comp[mask].to_numpy()
    order = ["growth_plate", "tumor_associated"]
    sub = ad.AnnData(X=X, var=pd.DataFrame(index=a.var_names.copy()),
                     obs=pd.DataFrame({"oc_compartment": pd.Categorical(labels, categories=order)},
                                      index=a.obs_names[mask]))
    sc.tl.rank_genes_groups(sub, groupby="oc_compartment", method="wilcoxon", use_raw=False)
    df = sc.get.rank_genes_groups_df(sub, group=None)
    groups = [X[labels == g] for g in order]
    mean_cp10k = np.vstack([np.asarray(Xg.astype(np.float64).expm1().mean(axis=0)).ravel() for Xg in groups])
    pct = np.vstack([100 * np.asarray((Xg > 0).sum(axis=0)).ravel() / Xg.shape[0] for Xg in groups])
    own = (df["group"].astype(str).to_numpy() == order[1]).astype(int)                # row of the group, 1 - own: rest
    gi = sub.var_names.get_indexer(df["names"].astype(str))
    df["log2fc_mean_cp10k"] = np.round(np.log2((mean_cp10k[own, gi] + 1) / (mean_cp10k[1 - own, gi] + 1)), 4)
    df["pct_group"] = np.round(pct[own, gi], 3)
    df["pct_rest"] = np.round(pct[1 - own, gi], 3)
    df["selection_gene"] = df["names"].astype(str).isin(set(selection_genes)).to_numpy()
    return df[COMPARTMENT_DE_COLUMNS].reset_index(drop=True)


def compartment_summary(a, roi_key="in_marrow_roi", samples=None):
    """
    [FIX] Figure 4h compartments as computed (Figure4h_compartment_summary.csv): one row for all sections ("all"),
    then one per section, with the threshold actually used (a.uns["compartment_threshold"], the joint percentile of the
    split score; on the SR007144 data it is 0, so 'above the 90th percentile' means at least one UMI of a chondrocyte
    gene), the percentage of bins above it, the bins of each compartment, how many of them lie outside the region of
    interest (the compartments are not restricted to it; 39 of 116 growth-plate and 595 of 2,252 tumor-associated bins
    on the SR007144 data), and the median UMI per bin of each compartment (growth-plate bins are about twice as deep).
    """
    thr = float(a.uns.get("compartment_threshold", np.nan))
    comp = a.obs["oc_compartment"].astype(str).to_numpy()
    gp = _as_bool(a.obs["growth_plate"]) if "growth_plate" in a.obs else np.zeros(a.n_obs, dtype=bool)
    roi = _as_bool(a.obs[roi_key]) if roi_key in a.obs else np.ones(a.n_obs, dtype=bool)
    umi = a.obs["total_counts"].to_numpy(dtype=float)
    smp = a.obs["sample"].astype(str).to_numpy()
    samples = list(samples) if samples else list(pd.unique(smp))
    rows = []
    for s, m in [("all", np.ones(a.n_obs, dtype=bool))] + [(s, smp == s) for s in samples]:
        g, t = m & (comp == "growth_plate"), m & (comp == "tumor_associated")
        rows.append({"sample": s, "threshold": thr, "threshold_is_zero": bool(thr == 0),
                     "pct_bins_above_threshold": round(100 * float(gp[m].mean()), 3) if m.any() else np.nan,
                     "n_growth_plate": int(g.sum()), "n_tumor_associated": int(t.sum()),
                     "n_non_growth_plate_control": int((m & (comp == "non_growth_plate_control")).sum()),
                     "n_growth_plate_outside_roi": int((g & ~roi).sum()),
                     "n_tumor_associated_outside_roi": int((t & ~roi).sum()),
                     "median_umi_growth_plate": float(np.median(umi[g])) if g.any() else np.nan,
                     "median_umi_tumor_associated": float(np.median(umi[t])) if t.any() else np.nan})
    return pd.DataFrame(rows, columns=COMPARTMENT_SUMMARY_COLUMNS)


# ============================================================================
# 7. HOST RECEPTOR PANEL AND ADIPOCYTE GENES (Figure 6b, Supplementary Figure S4)  [RUN]
# ============================================================================

def receptor_panel(a, panel=RECEPTOR_PANEL, n_reps=5, control="LL-2-1", tumour="RL-1-5", roi_key="in_marrow_roi",
                   reference=None, mode="cap", seed=RANDOM_SEED, target_sum=1e4, samples=None, progress=None):
    """
    [RUN] Percentage of ROI bins detecting each gene, raw and after depth matching (same downsampling as
    depth_matched(), five replicates), with Wilson intervals on the depth-matched rate. Used for the host receptor panel
    (Figure 6b) and, with panel=ADIPOCYTE_GENES, for Supplementary Figure S4.
    [FIX] fold_dm is computed for every section against the reference section (the control and tumour ids were
    hard-coded); a section without ROI bins gives NaN instead of crashing; the expression values are no longer added to
    a.obs as expr_<gene> columns (they bloated the h5ad); depth matching uses the per-bin subsampling of
    depth_matched() with the same replicate seeds, so both analyses see identical downsampled matrices.
    """
    from statsmodels.stats.proportion import proportion_confint
    reference = reference or control
    targets, _, med = depth_match_targets(a, mode)
    genes = [g for g in panel if g in a.var_names]
    missing = [g for g in panel if g not in a.var_names]
    if missing:
        _warn(progress, f"panel genes absent from the data and left out: {', '.join(missing)}")
    gi = a.var_names.get_indexer(genes)
    raw = a.layers["lognorm"][:, gi]
    raw = raw.toarray() if sp.issparse(raw) else np.asarray(raw)
    smp = a.obs["sample"].astype(str).to_numpy()
    samples = list(samples) if samples else list(pd.unique(smp))
    roi = _as_bool(a.obs[roi_key])
    masks = {s: (smp == s) & roi for s in samples}
    det = {s: np.full((n_reps, len(genes)), np.nan) for s in samples}
    mean = {s: np.full((n_reps, len(genes)), np.nan) for s in samples}
    for rep in range(n_reps):
        X = downsample_counts_matrix(a.layers["counts"], targets, seed + rep)
        lib = np.asarray(X.sum(axis=1)).ravel()
        Xg = X[:, gi].toarray()
        ln = np.log1p(Xg / np.maximum(lib, 1)[:, None] * target_sum)
        for s, m in masks.items():
            if m.any():
                det[s][rep] = (Xg[m] > 0).mean(axis=0)
                mean[s][rep] = ln[m].mean(axis=0)
    rows = []
    for j, g in enumerate(genes):
        for s in samples:
            m = masks[s]
            n = int(m.sum())
            p_dm = float(np.mean(det[s][:, j])) if n else np.nan
            npos = int(round(p_dm * n)) if n else 0
            lo, hi = proportion_confint(npos, n, method="wilson") if n else (np.nan, np.nan)
            rows.append({"gene": g, "role": panel[g], "sample": s, "roi": roi_key,
                         "mean_raw": float(raw[m, j].mean()) if n else np.nan,
                         "pct_detected_raw": 100 * float((raw[m, j] > 0).mean()) if n else np.nan,
                         "mean_dm": float(np.mean(mean[s][:, j])) if n else np.nan, "pct_detected_dm": 100 * p_dm,
                         "n_positive_dm": npos, "n_bins": n, "ci_lo_dm": 100 * lo, "ci_hi_dm": 100 * hi,
                         "median_umi_per_bin": float(med.get(s, np.nan))})
    df = pd.DataFrame(rows, columns=["gene", "role", "sample", "roi", "mean_raw", "pct_detected_raw", "mean_dm",
                                     "pct_detected_dm", "n_positive_dm", "n_bins", "ci_lo_dm", "ci_hi_dm",
                                     "median_umi_per_bin"])
    df["fold_dm"] = np.nan
    for g in genes:
        c = df[(df.gene == g) & (df["sample"] == reference)]["pct_detected_dm"]
        if not len(c) or not float(c.iloc[0]) > 0:
            continue
        t = df[(df.gene == g) & (df["sample"] != reference)]
        df.loc[t.index, "fold_dm"] = t["pct_detected_dm"].values / float(c.iloc[0])
    return df


# ============================================================================
# 8. SPATIAL ANALYSES ON THE JOINT OBJECT
# ============================================================================

def _array_xy_um(adata, bin_um):
    """[RUN] Bin positions in um from the array grid (adjacent bins are exactly bin_um apart; no pixel calibration)."""
    return np.column_stack([adata.obs["array_col"].to_numpy(float), adata.obs["array_row"].to_numpy(float)]) * bin_um


UNCOUPLING_COLUMNS = ["sample", "n_a", "n_b", "n_units_roi", "ratio_obs_null_mean", "ratio_obs_null_lo",
                      "ratio_obs_null_hi", "p_greater_than_null_mean", "median_nn_um", "mean_nn_um", "null_mean_um",
                      "status", "ratio_of_ratios_vs_reference", "note"]


def uncoupling_metric(adata, a_key="pooled_dm_pos_osteoblast", b_key="pooled_dm_pos_osteoclast", bin_um=16.0,
                      n_perm=199, roi_key="in_marrow_roi", min_pos=10, reference_sample="LL-2-1", seed=RANDOM_SEED,
                      samples=None):
    """
    [RUN] Figure 4i. Distance from each osteoblast-positive bin to the nearest osteoclast-positive bin (pooled,
    depth-matched calls), within the ROI of each section, against a random-placement null that draws the same numbers of
    positive bins in the same ROI (199 replicates). Primary statistic: observed mean / null mean per limb (< 1 closer than
    chance, > 1 farther); between limbs, the ratio of those ratios. Raw distances are density-dominated and are not
    compared between sections. Returns (summary, distributions); null counts in `distributions` are per replicate.
    [APP] Sections, reference and bin size from the configuration (bin_um was fixed at 16).
    [FIX] Every summary column exists even when a section is not estimable, so readers and figures no longer fail on a
    missing column; the ratio of ratios is NaN when the reference section is not estimable.
    """
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(seed)
    res, dist_rows = [], []

    def _hist(d, s, source, scale=1.0):
        v, c = np.unique(np.round(d, 2), return_counts=True)
        dist_rows.extend({"sample": s, "source": source, "distance_um": float(vi), "count": float(ci) * scale}
                         for vi, ci in zip(v, c))

    smp = adata.obs["sample"].astype(str).to_numpy()
    samples = list(samples) if samples else list(pd.unique(smp))
    roi = _as_bool(adata.obs[roi_key])
    xy_all = _array_xy_um(adata, bin_um)
    A_all, B_all = _as_bool(adata.obs[a_key]), _as_bool(adata.obs[b_key])
    for s in samples:
        m = (smp == s) & roi
        xy, A, B = xy_all[m], A_all[m], B_all[m]
        nA, nB = int(A.sum()), int(B.sum())
        if nA < min_pos or nB < min_pos:
            res.append({"sample": s, "n_a": nA, "n_b": nB, "n_units_roi": int(m.sum()),
                        "status": f"not estimable (n_positive < {min_pos})"})
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
    df = pd.DataFrame(res).reindex(columns=UNCOUPLING_COLUMNS)
    ref = df[(df["sample"] == reference_sample) & (df["status"] == "estimable")]
    if not ref.empty:
        df["ratio_of_ratios_vs_reference"] = df["ratio_obs_null_mean"] / ref.iloc[0]["ratio_obs_null_mean"]
    df["note"] = "raw distances (median_nn_um, mean_nn_um) are density-dominated; do not compare between sections"
    return df, pd.DataFrame(dist_rows, columns=["sample", "source", "distance_um", "count"])


def _nearest_target_distance(P, target, exclude_self=False):
    """[APP] Distance from every point to the nearest target point; with exclude_self a target point is measured to the
    nearest OTHER target point. No target point: infinite distances (an empty tree is never queried)."""
    from scipy.spatial import cKDTree
    n_t = int(target.sum())
    if n_t == 0:
        return np.full(len(P), np.inf)
    tree = cKDTree(P[target])
    if not exclude_self:
        return tree.query(P, k=1)[0]
    k = 2 if n_t > 1 else 1
    d = np.asarray(tree.query(P, k=k)[0]).reshape(len(P), k)
    out = d[:, 0].copy()
    out[target] = d[target, 1] if k == 2 else np.inf
    return out


def _pct(mask):
    return round(100 * float(mask.mean()), 1) if mask.size else np.nan


def proximity_table(a, src, dst, mode="assignment", radius_um=50.0, exclude_self=False, roi_key="in_marrow_roi",
                    samples=None, labels=None):
    """
    [RUN] Results and Methods, "within 50 um": fraction of bins of one signature whose nearest bin of the other
    signature lies within 50 um (center to center), against the same fraction for all ROI bins.
      fibroblast-assigned -> osteoblast-assigned bins (deterministic assignment, column 'program'): 41.8/41.3/67.3/5.6%
      macrophage-positive -> osteoclast-positive bins (per-bin positivity): 115 of 405, 28.4% against 17.1%
    Positions are full-resolution pixels times each section's microns per pixel.
    [APP] Any pair of signatures and radius (configuration proximity.pairs; the defaults reproduce the two original
    tables, column names included). mode "assignment" uses the 'program' column, "positivity" the pos_ columns.
    exclude_self=True measures a bin that is itself in the target set to the nearest OTHER target bin (by default a bin
    positive for both signatures is its own nearest target, at 0 um); co-positive bins are counted either way.
    [FIX] Microns per pixel come from each section's scalefactors_json.json (UM_PER_PX was keyed by sample id), and a
    section without target bins no longer queries an empty tree.
    """
    smp = a.obs["sample"].astype(str).to_numpy()
    samples = list(samples) if samples else list(pd.unique(smp))
    labels = labels or {}
    roi = _as_bool(a.obs[roi_key])
    xy = np.asarray(a.obsm["spatial"], dtype=float)
    R = f"{radius_um:g}"
    rows = []
    for s in samples:
        m = (smp == s) & roi
        P = xy[m] * microns_per_pixel(a, s)
        if mode == "assignment":
            prog = a.obs["program"].astype(str).to_numpy()[m]
            fr, to = prog == src, prog == dst
        else:
            fr, to = _as_bool(a.obs[f"pos_{src}"])[m], _as_bool(a.obs[f"pos_{dst}"])[m]
        within = _nearest_target_distance(P, to, exclude_self) <= radius_um
        if mode == "assignment":
            rows.append({"section": labels.get(s, s), f"{src}_assigned_bins": int(fr.sum()),
                         f"{dst}_assigned_bins": int(to.sum()), "region_of_interest_bins": int(m.sum()),
                         f"percent_{src}_within_{R}um_of_{dst}": _pct(within[fr]),
                         f"percent_all_bins_within_{R}um_of_{dst}": _pct(within)})
        else:
            rows.append({"section": labels.get(s, s), f"{dst}_positive_bins": int(to.sum()),
                         f"{src}_positive_bins": int(fr.sum()), f"{src}_within_{R}um_of_{dst}": int(within[fr].sum()),
                         f"percent_{src}_within_{R}um": _pct(within[fr]), "region_of_interest_bins": int(m.sum()),
                         f"percent_all_bins_within_{R}um": _pct(within),
                         f"{src}_also_{dst}_positive_bins": int((fr & to).sum())})
    return pd.DataFrame(rows)


def proximity_tables(a):
    """[RUN] The two tables of the manuscript (see proximity_table)."""
    labels = {s: SECTION.get(s, s) for s in pd.unique(a.obs["sample"].astype(str))}
    return (proximity_table(a, "fibroblast", "osteoblast", "assignment", labels=labels),
            proximity_table(a, "macrophage", "osteoclast", "positivity", labels=labels))


def _morans_i_builtin(sub, genes):
    """
    [APP] squidpy.gr.spatial_neighbors(coord_type="generic", delaunay=True) + spatial_autocorr(mode="moran",
    n_perms=None) without squidpy: the same Delaunay graph (scipy), row-normalised, scanpy.metrics.morans_i (the function
    squidpy calls) and squidpy's analytic normal-approximation p-values with Benjamini-Hochberg correction. Genes with
    zero variance get I = NaN and NaN p-values (squidpy leaves those p-values undefined).
    """
    from scipy import stats
    from scipy.spatial import Delaunay
    from sklearn.preprocessing import normalize
    coords = np.asarray(sub.obsm["spatial"], dtype=float)
    n = coords.shape[0]
    indptr, indices = Delaunay(coords).vertex_neighbor_vertices
    g = sp.csr_matrix((np.ones_like(indices, dtype=np.float32), indices, indptr), shape=(n, n))
    g.eliminate_zeros()
    normalize(g, norm="l1", axis=1, copy=False)
    score = np.asarray(sc.metrics.morans_i(g, sub[:, genes].X.T), dtype=float)
    s0 = g.sum()
    t = g.T + g
    s1 = t.multiply(t).sum() / 2.0
    s2 = float(((np.asarray(g.sum(1)).ravel() + np.asarray(g.sum(0)).ravel()) ** 2).sum())
    v_norm = (n * n * s1 - n * s2 + 3 * s0 * s0) / ((n - 1) * (n + 1) * s0 * s0) - (1.0 / (n - 1)) ** 2
    z = (score - (-1.0 / (n - 1))) / np.sqrt(v_norm)
    p = np.where(z > 0, 1 - stats.norm.cdf(z), stats.norm.cdf(z))
    p[np.isnan(z)] = np.nan
    df = pd.DataFrame({"I": score, "pval_norm": p, "var_norm": v_norm}, index=pd.Index(genes))
    df["pval_norm_fdr_bh"] = _bh_on_finite(p)
    return df.sort_values("I", ascending=False)


def _bh_on_finite(p):
    """[FIX] Benjamini-Hochberg adjusted P values over the defined P values only (NaN stays NaN). statsmodels'
    multipletests returns NaN for every gene as soon as one P value is NaN."""
    from statsmodels.stats.multitest import multipletests
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], alpha=0.05, method="fdr_bh")[1]
    return out


def morans_i(adata, n_svg=2000, n_consensus=200, samples=None, progress=None):
    """
    [SHAMIM, MODIFIED] step 13, "Spatial statistics".
    Original: sq.gr.spatial_neighbors(adata) and sq.gr.spatial_autocorr(adata, mode="moran") on the whole object.
    Now computed PER SECTION (the two sections' coordinate ranges overlap, so a graph on the joint object would join
    bins on different slides; Methods), with a Delaunay graph on the bin coordinates, for the top 2,000 highly variable
    genes chosen once on the joint object. The combined statement is the set of genes in the top 200 by Moran's I in
    both sections (Supplementary Figure S2: 89 genes).
    [APP] squidpy is used when it is installed, otherwise _morans_i_builtin (same graph, statistic and p-values).
    [FIX] The highly variable genes are chosen on a light object holding only the log-normalized matrix; the whole
    joint object (four copies of the matrix) used to be copied for it.
    [FIX] With squidpy, pval_norm_fdr_bh is recomputed over the defined P values only (_bh_on_finite). A gene with zero
    variance in a section has a NaN P value, and squidpy passed it to multipletests, which made the adjusted P value NaN
    for every gene of that section (all RL-1-5 genes on the SR007144 data: 12 HVGs are constant there).
    """
    try:
        import squidpy as sq
    except ImportError:
        sq = None
        _warn(progress, "squidpy is not installed: Moran's I is computed with the built-in implementation")
    hv = ad.AnnData(X=adata.layers["lognorm"], var=pd.DataFrame(index=adata.var_names.copy()))
    sc.pp.highly_variable_genes(hv, flavor="seurat", n_top_genes=min(n_svg, adata.n_vars))
    genes = hv.var_names[hv.var.highly_variable].tolist()
    del hv
    smp = adata.obs["sample"].astype(str).to_numpy()
    samples = list(samples) if samples else list(pd.unique(smp))
    morans = []
    for s in samples:
        m = smp == s
        if m.sum() < 4:
            _warn(progress, f"Moran's I: section {s} has fewer than four bins, skipped")
            continue
        sub = ad.AnnData(X=adata.layers["lognorm"][m], obs=pd.DataFrame(index=adata.obs_names[m]),
                         var=pd.DataFrame(index=adata.var_names.copy()),
                         obsm={"spatial": np.asarray(adata.obsm["spatial"])[m]})
        if sq is not None:
            sq.gr.spatial_neighbors(sub, coord_type="generic", delaunay=True)
            # only Moran's I itself is reported, so no permutation P values (they do not change I)
            sq.gr.spatial_autocorr(sub, mode="moran", genes=genes, n_perms=None, n_jobs=1)
            res = sub.uns["moranI"]
            if "pval_norm" in res:
                res["pval_norm_fdr_bh"] = _bh_on_finite(res["pval_norm"].to_numpy())
        else:
            res = _morans_i_builtin(sub, genes)
        mdf = res.reset_index().rename(columns={"index": "gene"})
        mdf.insert(0, "sample", s)
        mdf["rank_in_sample"] = mdf["I"].rank(ascending=False, method="min", na_option="keep").astype("Int64")
        morans.append(mdf)
        del sub
    if not morans:
        return pd.DataFrame(columns=["sample", "gene", "I", "rank_in_sample"]), pd.DataFrame(columns=["gene", "mean_rank"])
    moran = pd.concat(morans, ignore_index=True)
    done = list(pd.unique(moran["sample"]))
    top = moran[(moran["rank_in_sample"] <= n_consensus).fillna(False).astype(bool)]
    shared = top.groupby("gene")["sample"].nunique()
    shared = shared[shared == len(done)].index
    cons = moran[moran["gene"].isin(shared)].pivot_table(index="gene", columns="sample", values="I").add_prefix("moranI_")
    cons["mean_rank"] = moran[moran["gene"].isin(shared)].groupby("gene")["rank_in_sample"].mean().astype(float)
    return moran, cons.sort_values("mean_rank").reset_index()


def leiden_signature_matrix(a, primary=PRIMARY, sigs=SIG, label_of=None, min_bins=20, samples=None,
                            roi_key="in_marrow_roi"):
    """
    [RUN] Supplementary Figure S5: mean signature score per Leiden domain, z-scored across domains. Each domain is
    named by the signature most enriched among its assigned bins relative to the whole dataset (at least 20 bins).
    [APP] Signature labels from the configuration; works without the tissue assignment (domains keep their numbers).
    [FIX] A signature whose mean score is identical in every domain gets NaN instead of a division by zero.
    [FIX] The unresolved table (UNRESOLVED_COLUMNS) has a first row sample="all" (the pooled share, as before) and one
    row per section, each also within the region of interest (roi_key; NaN without it). The pooled share hid the
    difference between sections (SR007144: 33.3% pooled, 18.3% control and 42.3% tumor-bearing; ROI 18.6% and 47.8%).
    """
    def lbl(c):
        return label_of(c) if label_of else c.replace("_", " ").replace("neutrophil", "granulocyte")
    ld = a.obs["leiden"].astype(str).values
    names = {c: c for c in pd.unique(ld)}
    if "program" in a.obs:
        prog = a.obs["program"].astype(str).values
        ct = pd.crosstab(pd.Series(ld, name="d"), pd.Series(prog, name="p")).drop(columns="ambiguous", errors="ignore")
        if ct.shape[1]:
            enr = (ct.div(ct.sum(1).replace(0, np.nan), axis=0) / (ct.sum(0) / ct.values.sum())).where(ct >= min_bins)
            names = {c: c if c in ("unresolved", "unassigned") or c not in enr.index or enr.loc[c].isna().all()
                     else f"{lbl(enr.loc[c].idxmax())}-enriched (domain {c})" for c in pd.unique(ld)}
    cols = [f"score_{n}" for n in sigs if n in primary and f"score_{n}" in a.obs]
    M = a.obs.groupby("leiden", observed=True)[cols].mean()
    M = (M - M.mean()) / M.std().replace(0, np.nan)
    M.index = [names[str(i)] for i in M.index]
    M.columns = [lbl(c.replace("score_", "")) for c in M.columns]
    un = ld == "unresolved"
    smp = a.obs["sample"].astype(str).to_numpy()
    roi = _as_bool(a.obs[roi_key]) if roi_key in a.obs else None
    samples = list(samples) if samples else list(pd.unique(smp))
    rows = []
    for s, m in [("all", np.ones(len(ld), dtype=bool))] + [(s, smp == s) for s in samples]:
        r = {"sample": s, "unresolved_bins": int(un[m].sum()), "all_bins": int(m.sum()),
             "percent": round(100 * un[m].mean(), 1) if m.any() else np.nan,
             "unresolved_bins_roi": np.nan, "all_bins_roi": np.nan, "percent_roi": np.nan}
        if roi is not None:
            mr = m & roi
            r.update(unresolved_bins_roi=int(un[mr].sum()), all_bins_roi=int(mr.sum()),
                     percent_roi=round(100 * un[mr].mean(), 1) if mr.any() else np.nan)
        rows.append(r)
    return M.reset_index().rename(columns={"index": "leiden_domain"}), pd.DataFrame(rows, columns=UNRESOLVED_COLUMNS)


def geneset_sensitivity(a, n_reps=5, variants=None, names=("osteoclast", "osteoblast"),
                        comp=("erythroid", "neutrophil", "fibroblast", "chondrocyte"), sigs=SIG, roi_key="in_marrow_roi",
                        mode="cap", min_core=2, n_perm=50, n_perm_resampled=20, pct=95.0, n_sub=20000, ctrl_size=50,
                        seed=RANDOM_SEED, target_sum=1e4, samples=None, progress=None):
    """
    [RUN] Methods: bone-lineage positivity recomputed under three gene sets with the calling rule held constant
    (canonical two-gene markers; the lists of the original study design; the scoring signatures of Table 1), at native
    depth and after depth matching, within the ROI (about 30 runs of scoring and null thresholds; 1-2 minutes).
    Every listed gene of a variant is both scored and counted toward the two-gene rule (core = extended = the list, as
    in the original code). The osteoclast list of the original study design equals the scored (extended) list of Table 1,
    so its rows differ from the table1_signatures rows only in that Mmp9 and Nfatc1 also count toward the two-gene rule.
    [FIX] Variant genes absent from the data are reported (warning; Bglap of the original osteoblast list is not in the
    probe set), and every row names the genes used and missing (columns genes_used, genes_missing).
    [FIX] The depth-matched replicates are scored from their own downsampled matrices (they were scored from the native
    .raw, see score_all), with the downsampling and null sizes of the main analysis (native calls n_perm, depth-matched
    calls n_perm_resampled; all calls used 20 before), so the table1_signatures rows equal Figures 4g/6a. A variant that
    cannot be scored reports zero positives with a warning instead of the main calls, and a section without ROI bins
    gives NaN instead of a division by zero.
    [APP] Variants, signatures and competitors from the configuration (sensitivity.*).
    """
    from statsmodels.stats.proportion import proportion_confint
    if variants is None:
        variants = pc.DEFAULT_CONFIG["sensitivity"]["variants"]
    names, comp = list(names), [c for c in comp if c not in names]
    roi = _as_bool(a.obs[roi_key])
    smp = a.obs["sample"].astype(str).to_numpy()
    samples = list(samples) if samples else list(pd.unique(smp))
    targets, _, _ = depth_match_targets(a, mode)
    kw = dict(min_core=min_core, pct=pct, n_sub=n_sub, ctrl_size=ctrl_size, seed=seed, progress=progress)
    rows = []
    for vname, override in variants.items():
        vs = copy.deepcopy(sigs)
        for n, genes in (override or {}).items():
            if n in vs:
                vs[n].update({"core": list(genes), "extended": list(genes)})
        vs = {k: vs[k] for k in names + comp if k in vs}
        used = {n: [g for g in vs[n]["extended"] if g in a.var_names] if n in vs else [] for n in names}
        missing = {n: [g for g in vs[n]["extended"] if g not in a.var_names] if n in vs else [] for n in names}
        for n in names:
            if missing[n]:
                _warn(progress, f"gene-set sensitivity, {vname}: {n} genes absent from the data: "
                                f"{', '.join(missing[n])}")
        b = _lognorm_object(a.layers["counts"], a, ["sample"], target_sum)
        b = call_positive(score_all(b, vs, ctrl_size, seed), names, sigs=vs, n_perm=n_perm, **kw)
        nat = {(n, s): int(b.obs[f"pos_{n}"].values[roi & (smp == s)].sum()) for n in names for s in samples}
        dm = {k: [] for k in nat}
        for rep in range(n_reps):
            X = downsample_counts_matrix(a.layers["counts"], targets, seed + rep)
            c = _lognorm_object(X, a, ["sample"], target_sum)
            c = call_positive(score_all(c, vs, ctrl_size, seed), names, sigs=vs, n_perm=n_perm_resampled, **kw)
            for (n, s) in dm:
                dm[(n, s)].append(int(c.obs[f"pos_{n}"].values[roi & (smp == s)].sum()))
        for (n, s), k in nat.items():
            N = int((roi & (smp == s)).sum())
            lo, hi = proportion_confint(k, N, method="wilson") if N else (np.nan, np.nan)
            d = np.array(dm[(n, s)]) / N * 100 if N else np.full(n_reps, np.nan)
            rows.append({"gene_set": vname, "signature": n, "sample": s, "n_units_roi": N, "n_positive_native": k,
                         "pct_native": round(100 * k / N, 4) if N else np.nan, "ci_lo": round(100 * lo, 4),
                         "ci_hi": round(100 * hi, 4), "pct_depth_matched_mean": round(float(d.mean()), 4),
                         "pct_depth_matched_sd": round(float(d.std(ddof=1)), 4) if n_reps > 1 else np.nan,
                         "genes_used": ",".join(used[n]), "genes_missing": ",".join(missing[n])})
    return pd.DataFrame(rows, columns=GENESET_SENSITIVITY_COLUMNS)


def table1_signatures(signatures=None):
    """[RUN] Table 1 as a CSV, generated from SIG so that the table and the code cannot diverge.
    [APP] From the configuration's signatures (table_label, canonical markers)."""
    signatures = signatures or _DEFAULT["signatures"]
    return pd.DataFrame([{"signature": d.get("table_label") or d.get("label") or k,
                          "canonical_markers": ", ".join(d.get("canonical", [])),
                          "signature_genes": ", ".join(d["extended"]), "genes_counted_toward_two_gene_rule": ", ".join(d["core"]),
                          "excluded_genes": "; ".join(f"{g}: {r}" for g, r in d.get("excluded", {}).items())}
                         for k, d in signatures.items()])


def _slug(cfg, name):
    """[APP] Per-bin column name of a signature: its label with spaces as underscores (neutrophil -> granulocyte)."""
    return pc.signature_label(cfg, name).replace(" ", "_")


def per_bin_table(a, cfg=None):
    """[RUN] One row per 16 um bin (Manuscript/data/per_bin): position in um, ROI flag, assignment, positivity.
    [APP] Signatures, labels, classes and microns per pixel from the configuration and the data (were hard-coded)."""
    cfg = cfg or pc.default_config()
    o = a.obs
    smp = o["sample"].astype(str).to_numpy()
    mpp = {s: microns_per_pixel(a, s) for s in pd.unique(smp)}
    um = np.array([mpp[s] for s in smp])[:, None] * np.asarray(a.obsm["spatial"], dtype=float)
    labels = {s["id"]: s.get("label") or s["id"] for s in cfg["samples"]}
    drawn = {c["name"] for c in cfg["classes"]}
    group = {k: (d.get("class") if d.get("class") in drawn else "removed") for k, d in cfg["signatures"].items()}
    group["ambiguous"] = "unassigned"
    program = o["program"].astype(str) if "program" in o else pd.Series("ambiguous", index=o.index)
    pb = pd.DataFrame({"barcode": o.index.astype(str), "sample_id": smp, "section": [labels.get(s, s) for s in smp],
                       "x_um": um[:, 0].round(2), "y_um": um[:, 1].round(2),
                       "in_region_of_interest": _as_bool(o["in_marrow_roi"]), "total_umi": o["total_counts"].values,
                       "assigned_signature": [p if p == "ambiguous" else _slug(cfg, p) for p in program],
                       "tissue_class": [group.get(p, "removed") for p in program],
                       "leiden_domain": o["leiden"].astype(str).values if "leiden" in o else ""})
    for n in cfg["primary"]:
        if f"pos_{n}" in o:
            pb[f"positive_{_slug(cfg, n)}"] = _as_bool(o[f"pos_{n}"])
    unc = cfg["uncoupling"]
    prefix = {"pooled_dm": "pooled_dm_", "pooled": "pooled_", "native": ""}[unc["calls"]]
    suffix = {"pooled_dm": "_pooled_depthmatched", "pooled": "_pooled", "native": "_native"}[unc["calls"]]
    for n in (unc["b"], unc["a"]):
        if f"{prefix}pos_{n}" in o and prefix:
            pb[f"positive_{_slug(cfg, n)}{suffix}"] = _as_bool(o[f"{prefix}pos_{n}"])
    # [APP] sensitivity analyses: share of depth-matching draws in which a bin is pooled_dm-positive (the column above
    # is draw 0 only), and the cross-captured A549 lesion flag
    for n in (unc["b"], unc["a"]):
        if f"pooled_dm_pos_frac_{n}" in o:
            col = f"positive_fraction_{_slug(cfg, n)}_pooled_depthmatched_draws"
            pb[col] = o[f"pooled_dm_pos_frac_{n}"].to_numpy(float)
    if "in_lesion" in o:
        pb["lesion_umi_pooled"] = o["lesion_umi_pooled"].to_numpy(float)
        pb["in_lesion"] = _as_bool(o["in_lesion"])
    return pb


# ============================================================================
# MAIN
# ============================================================================

def load_sections(root, bin_str, samples=None, qc=None, progress=None):
    """[SHAMIM, MODIFIED] steps 1-3 for each section; returns the filtered sections and the QC sweep.
    [APP] Sections, folders and thresholds from the configuration; root may be a configuration (dict)."""
    cfg = root if isinstance(root, dict) else {**pc.default_config(), "data_root": str(root)}
    qc = qc or cfg["qc"]
    adatas, sweeps, loaded = {}, [], {}
    for s in samples or cfg["samples"]:
        p = pc.resolve_bin_dir(cfg, s, bin_str)
        if p is None:
            raise FileNotFoundError(f"{s['id']}: no square_{bin_str}/filtered_feature_bc_matrix.h5 under "
                                    f"{Path(str(cfg['data_root'])) / str(s.get('path') or s['id'])}")
        a = qc_metrics(read_visium_hd(p, s["id"], s.get("condition"), s.get("microns_per_pixel")))
        loaded[s["id"]] = a.n_obs
        sweeps.append(qc_sweep(a, s["id"], qc.get("sweep_thresholds", (0, 5, 10, 25, 50, 100, 200))))
        a = qc_filter(a, qc["min_genes"], qc["min_counts"], qc["max_pct_mt"])
        print(f"[{s['id']}] {a.n_obs:,} bins retained at {bin_str}")
        if a.n_obs == 0:
            raise ValueError(f"{s['id']}: no bin passes quality control at {bin_str}; check the qc thresholds")
        adatas[s["id"]] = a
    return adatas, pd.concat(sweeps, ignore_index=True), loaded


class _Outputs:
    """[APP] Writes the output tables and records every file in manifest.json."""

    def __init__(self, out: Path):
        self.out = out
        self.files: list[dict] = []

    def add(self, name, kind, title, description, rows=None):
        self.files = [f for f in self.files if f["file"] != name]
        self.files.append({"file": name, "kind": kind, "title": title, "description": description, "rows": rows})

    def csv(self, df, name, title, description, index=False):
        df.to_csv(self.out / name, index=index)
        self.add(name, "table", title, description, int(len(df)))
        return df

    def write(self):
        (self.out / "manifest.json").write_text(json.dumps({"files": self.files}, indent=2), encoding="utf-8")


def _environment():
    env = {"python": platform.python_version()}
    for name, dist in (("scanpy", "scanpy"), ("anndata", "anndata"), ("numpy", "numpy"), ("pandas", "pandas"),
                       ("scipy", "scipy"), ("statsmodels", "statsmodels"), ("sklearn", "scikit-learn"),
                       ("squidpy", "squidpy"), ("harmonypy", "harmonypy"), ("igraph", "igraph"),
                       ("leidenalg", "leidenalg")):
        try:
            env[name] = metadata.version(dist)
        except metadata.PackageNotFoundError:
            env[name] = None
    return env


@contextlib.contextmanager
def _thread_limits(enabled):
    """[APP] Single-threaded numerics while the analysis runs (see the reproducibility note at the top)."""
    if not enabled:
        yield
        return
    with contextlib.ExitStack() as stack:
        try:
            import numba
            previous = numba.get_num_threads()
            numba.set_num_threads(1)
            stack.callback(numba.set_num_threads, previous)
        except (ImportError, ValueError):
            pass
        try:
            from threadpoolctl import threadpool_limits
            stack.enter_context(threadpool_limits(1))
        except ImportError:
            pass
        yield


def _steps(cfg):
    """[APP] (key, progress label, enabled) of every step of run_analysis_from_config, in order."""
    B, Q, st = _bin_label(cfg["bin_size"]), _bin_label(cfg["qc_bin_size"]), cfg["steps"]
    return [
        ("load", f"Load sections and quality control ({B} bins)", True),
        ("qc_figure_bins", f"Per-bin QC at {Q} (Supplementary Figure S1)", st["qc_figure_bins"]),
        ("normalize", "Normalization, embedding and Leiden domains" if st["leiden"] else "Normalization", True),
        ("scores", "Gene-set scores and positivity calls", True),
        ("roi", "Region of interest", True),
        ("receptor_panel", "Host receptor panel (Figure 6b)", st["receptor_panel"]),
        ("adipocyte_panel", "Adipocyte genes (Supplementary Figure S4)", st["adipocyte_panel"]),
        ("positivity", "Positivity at native depth (Figures 4g, 6a)", True),
        ("depth_matched", "Depth-matched positivity (Figures 4g, 6a)", st["depth_matched"]),
        ("leiden", "Leiden domain markers and screen", st["leiden"]),
        ("compartments", "Osteoclast compartments (Figure 4h)", st["compartments"]),
        ("assignment", "Spatially pooled scores and tissue assignment (Figure 3)", st["assignment"]),
        ("stability", "Leiden seed stability (Methods)", st["stability"]),
        ("domain_matrix", "Leiden domains by signature (Supplementary Figure S5)", st["leiden"]),
        ("uncoupling", "Distance to nearest partner (Figure 4i)", st["uncoupling"]),
        ("proximity", "Proximity tables", st["proximity"]),
        ("moran", "Moran's I (Supplementary Figure S2)", st["moran"]),
        ("sensitivity", "Gene-set sensitivity (Methods)", st["sensitivity"]),
        ("sensitivity_analyses", "Sensitivity analyses (second review)", st.get("sensitivity_analyses", False)),
        ("tables", "Table 1 and per-bin calls", True),
        ("write_h5ad", f"Write joint_{B}.h5ad", st["write_h5ad"]),
        ("project_fine_bins", f"Project labels onto {Q} bins", st["project_fine_bins"]),
        ("summary", "Summary", True),
    ]


def planned_steps(cfg) -> list[str]:
    """[APP] Progress labels run_analysis_from_config reports for this configuration, in order (disabled steps are
    reported as skipped)."""
    return [label for _, label, _ in _steps(pc.merge_with_defaults(cfg))]


def run_analysis_from_config(cfg, out, progress=None):
    """[APP] The whole analysis for a configuration (pipeline_config.py); writes every output into `out` (ARCHITECTURE.md
    section 3) and returns the joint AnnData object."""
    cfg = pc.effective_config(cfg)
    progress = progress or pc.Progress()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    pc.save_config(cfg, out / "config.json")
    with _thread_limits(cfg.get("deterministic", True)):
        return _run(cfg, out, progress)


def _resolve_mask_paths(cfg):
    masks = {}
    for sid, p in (cfg["roi"].get("masks") or {}).items():
        path = Path(str(p)).expanduser()
        if not path.is_absolute() and not path.exists():
            path = Path(str(cfg["data_root"])) / path
        if not path.is_file():
            raise FileNotFoundError(f"region of interest of {sid}: {path} not found")
        masks[sid] = path
    return masks


def _run(cfg, out, progress):
    steps = _steps(cfg)
    label = {k: lab for k, lab, _ in steps}
    enabled = {k: en for k, _, en in steps}
    progress.plan([lab for _, lab, _ in steps])
    sco, clu = cfg["scoring"], cfg["clustering"]
    sigs = {k: {"core": list(d["core"]), "extended": list(d["extended"]), "competes": list(d.get("competes") or []),
                "excluded": dict(d.get("excluded") or {})} for k, d in cfg["signatures"].items()}
    primary, programs, seed = cfg["primary"], cfg["programs"], int(cfg["seed"])
    samples = pc.sample_ids(cfg)
    names_by_id = {s["id"]: s.get("label") or s["id"] for s in cfg["samples"]}
    ref = cfg["reference_sample"]
    B, Q = _bin_label(cfg["bin_size"]), _bin_label(cfg["qc_bin_size"])
    call_kw = dict(sigs=sigs, min_core=sco["min_core"], pct=sco["null_percentile"], n_sub=sco["null_n_sub"],
                   ctrl_size=sco["ctrl_size"], seed=seed, progress=progress)
    target_sum = float(clu["target_sum"])
    W = _Outputs(out)
    W.add("config.json", "data", "Configuration", "The effective configuration of this run (pipeline_config.py).")
    done, tables = set(), {}

    def begin(key, needs=()):
        if not enabled[key]:
            progress.skip(label[key], "disabled in the configuration")
            return False
        missing = [label[n] for n in needs if n not in done]
        if missing:
            progress.skip(label[key], "needs " + "; ".join(missing))
            return False
        progress.step(label[key])
        return True

    try:
        # ---- 1. input and QC -------------------------------------------------- [SHAMIM, MODIFIED] steps 1-3
        begin("load")
        adatas, sweep, loaded = load_sections(cfg, cfg["bin_size"], progress=progress)
        W.csv(sweep, f"qc_threshold_sweep_{B}.csv", f"QC retention sweep ({B} bins)",
              "Bins retained at each minimum gene count, before filtering (Supplementary Figure S1).")
        mpp = {}
        for sid, x in adatas.items():
            sf = x.uns["spatial"][sid]["scalefactors"]
            mpp[sid] = float(sf["microns_per_pixel"]) if sf.get("microns_per_pixel") else UM_PER_PX.get(sid)
            if not mpp[sid]:
                raise ValueError(f"{sid}: microns_per_pixel missing from scalefactors_json.json; set it in the sample")
        a = ad.concat(adatas, label="sample_key", index_unique="-", uns_merge="unique", join="outer")  # [RUN] joint
        a.obs["sample"] = pd.Categorical(a.obs["sample"].astype(str), categories=samples)
        a.uns["microns_per_pixel"] = mpp
        del adatas
        done.add("load")

        # ---- 1b. per-bin QC at the QC bin size (Supplementary Figure S1) ------------------- [APP]
        if begin("qc_figure_bins"):
            per_bin, sweeps = [], []
            for s in cfg["samples"]:
                size = cfg["qc_bin_size"]
                d = pc.resolve_bin_dir(cfg, s, size)
                if d is None:
                    _warn(progress, f"{s['id']}: no square_{size} folder; Supplementary Figure S1 uses the "
                                    f"{cfg['bin_size']} bins")
                    size = cfg["bin_size"]
                    d = pc.resolve_bin_dir(cfg, s, size)
                x = qc_metrics(read_visium_hd(d, s["id"], s.get("condition"), s.get("microns_per_pixel")))
                per_bin.append(x.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].assign(sample=s["id"], bin_size=size)
                               [["sample", "bin_size", "n_genes_by_counts", "total_counts", "pct_counts_mt"]])
                sweeps.append(qc_sweep(x, s["id"], cfg["qc"]["sweep_thresholds"]).assign(bin_size=size))
                del x
            W.csv(pd.concat(per_bin, ignore_index=True), f"FigureS1_qc_per_bin_{Q}.csv.gz", f"Per-bin QC ({Q} bins)",
                  "Genes, UMIs and % mitochondrial counts of every bin before filtering (Supplementary Figure S1).")
            W.csv(pd.concat(sweeps, ignore_index=True), f"qc_threshold_sweep_{Q}.csv", f"QC retention sweep ({Q} bins)",
                  "Bins retained at each minimum gene count (Supplementary Figure S1, right panel).")
            done.add("qc_figure_bins")

        # ---- 2. normalization, embedding, Leiden ---------------------------------- [SHAMIM, MODIFIED] steps 4-8
        begin("normalize")
        a = normalize_and_cluster(a, resolution=clu["resolution"], n_top_genes=clu["n_top_genes"], n_pcs=clu["n_pcs"],
                                  n_neighbors=clu["n_neighbors"], harmony=clu["harmony"], target_sum=target_sum,
                                  scale_max_value=clu["scale_max_value"], n_iterations=clu["leiden_n_iterations"],
                                  seed=seed, cluster=enabled["leiden"], progress=progress)
        done.add("normalize")

        # ---- 3. signature scores and positivity ----------------------------------- [RUN]
        begin("scores")
        a = score_all(a, sigs, sco["ctrl_size"], seed)
        for r in a.uns["signature_coverage"].itertuples():
            if r.missing:
                _warn(progress, f"signature {r.signature}: genes not in the data: {r.missing.replace(',', ', ')}")
        a = call_positive(a, primary, n_perm=sco["n_perm"], **call_kw)
        done.add("scores")

        # ---- 4. region of interest ------------------------------------------------ [RUN]
        begin("roi")
        drawn = None
        if cfg["roi"]["use_drawn"]:
            from bone_pending_masks import rasterise_drawn_roi
            drawn, unrestricted = rasterise_drawn_roi(a, _resolve_mask_paths(cfg), cfg["roi"]["mask_frame"])
            for s in unrestricted:
                _warn(progress, f"{s}: no hand-drawn region of interest; only the transcriptional region applies")
        a = define_marrow_roi(a, exclude=cfg["roi"]["exclude"], drawn=drawn, pct=cfg["roi"]["exclude_percentile"],
                              progress=progress)
        roi_n = a.obs.groupby("sample", observed=True)["in_marrow_roi"].sum()
        for s in samples:
            if not roi_n.get(s, 0):
                _warn(progress, f"{s}: no bin in the region of interest")
        done.add("roi")
        mode = sco["depth_matching"]

        # ---- 5. host receptor panel (Figure 6b) and adipocyte genes (Supplementary Figure S4) ---------------- [RUN]
        panel_kw = dict(n_reps=sco["n_depth_reps"], reference=ref, mode=mode, seed=seed, target_sum=target_sum,
                        samples=samples, progress=progress)
        if begin("receptor_panel"):
            tables["receptor"] = W.csv(receptor_panel(a, cfg["receptor_panel"], **panel_kw),
                                       "Figure6b_receptor_panel_depth_matched.csv", "Host receptor panel",
                                       "Percentage of ROI bins detecting each receptor-panel gene, raw and depth-matched "
                                       "(mean of the replicates, Wilson 95% interval); fold_dm vs the reference section.")
            done.add("receptor_panel")
        if begin("adipocyte_panel"):
            W.csv(receptor_panel(a, {g: "adipocyte" for g in cfg["adipocyte_genes"]}, **panel_kw),
                  "FigureS4_adipocyte_genes_depth_matched.csv", "Adipocyte genes",
                  "Detection of the adipocyte genes in ROI bins, raw and depth-matched (Supplementary Figure S4).")
            done.add("adipocyte_panel")

        # ---- 6. positivity, native depth and depth-matched (Figures 4g, 6a) ------ [RUN]
        begin("positivity")
        tumor_ids = [s["id"] for s in cfg["samples"] if s.get("condition") in cfg["compartments"]["tumor_conditions"]]
        tables["native"] = W.csv(positivity_table(a, primary, roi_key="in_marrow_roi", reference=ref, samples=samples,
                                                  tumor_samples=tumor_ids),
                                 "Figure4g_6a_positivity_native_depth_roi.csv", "Positivity, native depth (ROI)",
                                 "Positive bins per signature and section inside the region of interest, Wilson 95% "
                                 "intervals (Figures 4g, 6a, left panels).")
        done.add("positivity")
        if begin("depth_matched"):
            tables["dm"] = W.csv(depth_matched(a, primary, n_reps=sco["n_depth_reps"], roi_key="in_marrow_roi", mode=mode,
                                               n_perm=sco["n_perm_resampled"], target_sum=target_sum, samples=samples,
                                               **call_kw),
                                 "Figure4g_6a_positivity_depth_matched.csv", "Positivity, depth-matched (ROI)",
                                 "Positivity after depth matching: mean and SD over the replicates (Figures 4g, 6a, "
                                 "right panels).")
            done.add("depth_matched")

        # ---- 7. Leiden domain markers and post-processing ------------------------ [SHAMIM, MODIFIED] steps 10-11; [RUN]
        if begin("leiden", ("normalize",)):
            de = domain_markers(a, n_top=clu["n_top_markers"], progress=progress)
            screen = W.csv(domain_screen(a, de, n_top=clu["screen_top"], ubiquitous=set(cfg["ubiquitous"]),
                                         genes_per_umi=clu["screen_genes_per_umi"],
                                         ubiquitous_flag=clu["screen_ubiquitous_flag"]),
                           "leiden_domain_screen.csv", "Leiden domain screen",
                           "Ubiquitous-transcript fraction of each domain's top markers and median genes per UMI.")
            a = collapse_unresolved(a, screen, ubiq_frac=clu["ubiquitous_fraction"])
            a = collapse_fragments(a, min_bins=clu["min_domain_bins"])
            W.csv(domain_markers(a, n_top=clu["n_top_markers"], progress=progress), "leiden_domain_markers.csv",
                  "Leiden domain markers", "Wilcoxon markers of each domain after merging (scanpy rank_genes_groups).")
            tables["unresolved_pct"] = round(100 * float((a.obs["leiden"].astype(str) == "unresolved").mean()), 1)
            done.add("leiden")

        # ---- 8. osteoclast compartments (Figure 4h) ------------------------------- [RUN]; test [SHAMIM, MODIFIED] step 10
        if begin("compartments"):
            cp = cfg["compartments"]
            if f"pos_{cp['signature']}" not in a.obs or f"score_{cp['split_by']}" not in a.obs:
                _warn(progress, f"osteoclast compartments need a positivity call for {cp['signature']} and a score for "
                                f"{cp['split_by']}; no test")
            else:
                a = split_osteoclast_compartments(a, cp["percentile"], cp["signature"], cp["split_by"],
                                                  cp["tumor_conditions"], progress)
                W.csv(compare_osteoclast_states(a, cp["de_sections"], cp["tumor_conditions"],
                                                cfg["signatures"][cp["split_by"]]["extended"], progress),
                      "Figure4h_osteoclast_compartment_DE.csv", "Osteoclast compartment DE",
                      "Wilcoxon rank-sum, growth-plate vs tumor-associated osteoclast-positive bins (Figure 4h), "
                      "with log2 fold change of mean CP10k, % of bins detecting each gene, and selection_gene (genes "
                      "of the split signature, higher in the growth-plate compartment by construction).")
                W.csv(compartment_summary(a, samples=samples), "Figure4h_compartment_summary.csv",
                      "Osteoclast compartments as computed",
                      f"Per section and overall: the {cp['split_by']} score threshold actually used (the "
                      f"{cp['percentile']:g}th percentile of all bins; 0 means at least one UMI of a {cp['split_by']} "
                      f"gene), % of bins above it, bins per compartment, compartment bins outside the region of "
                      f"interest and median UMI per bin of each compartment (Figure 4h).")
            done.add("compartments")

        # ---- 9. deterministic tissue assignment (Figure 3) ------------------------ [RUN]
        if begin("assignment"):
            pool_kw = dict(pool_k=sco["pool_k"], mode=mode, names=primary, n_perm=sco["n_perm_resampled"],
                           target_sum=target_sum, **call_kw)
            a = spatially_pooled_scores(a, **pool_kw)                              # native depth -> 'program' (Fig 3a)
            a = assign_dominant_program(a, "pooled_", "program", programs, sco["min_core"])
            a = spatially_pooled_scores(a, depth_match=True, **pool_kw)            # depth-matched -> 'program_dm' (3b)
            a = assign_dominant_program(a, "pooled_dm_", "program_dm", programs, sco["min_core"])
            tables["comp"] = W.csv(program_composition(a, "program", samples), "Figure3a_program_composition_native.csv",
                                   "Tissue assignment, native depth", "Bins per assigned signature and section (Figure 3a).")
            tables["comp_dm"] = W.csv(program_composition(a, "program_dm", samples),
                                      "Figure3b_program_composition_depth_matched.csv", "Tissue assignment, depth-matched",
                                      "Bins per assigned signature and section after depth matching (Figure 3b).")
            done.add("assignment")

        # ---- 10. Leiden stability (Methods; Supplementary Figure S5) -------------- [RUN]
        if begin("stability", ("leiden",)):
            ari = W.csv(clustering_stability(a, clu["resolution"], clu["stability_seeds"], clu["n_neighbors"],
                                             clu["n_pcs"], clu["leiden_n_iterations"]),
                        "leiden_seed_stability_adjusted_rand_index.csv", "Leiden seed stability",
                        "Adjusted Rand index between Leiden runs that differ only by seed.")
            tables["ari"] = float(ari["ARI"].median()) if len(ari) else None
            done.add("stability")
        if begin("domain_matrix", ("leiden",)):
            M, un = leiden_signature_matrix(a, primary, cfg["signatures"], lambda c: pc.signature_label(cfg, c),
                                            samples=samples)
            W.csv(M, "FigureS5_domain_by_signature_zscores.csv", "Leiden domains by signature",
                  "Mean signature score per Leiden domain, z-scored across domains (Supplementary Figure S5).")
            tables["unresolved"] = W.csv(un, "FigureS5_unresolved_bins.csv", "Unresolved bins",
                                         "Bins in the 'unresolved' domain: all sections (first row, sample 'all'), "
                                         "then each section, in all bins and in the region of interest.")
            done.add("domain_matrix")

        # ---- 11. spatial analyses -------------------------------------------------- [RUN]; Moran [SHAMIM, MODIFIED] step 13
        unc = cfg["uncoupling"]
        prefix = {"pooled_dm": "pooled_dm_", "pooled": "pooled_", "native": ""}[unc["calls"]]
        if begin("uncoupling", ("assignment",) if prefix else ()):
            if f"{prefix}pos_{unc['a']}" not in a.obs or f"{prefix}pos_{unc['b']}" not in a.obs:
                _warn(progress, f"distance to nearest partner: {unc['a']} and {unc['b']} need positivity calls")
            else:
                summ, dist = uncoupling_metric(a, f"{prefix}pos_{unc['a']}", f"{prefix}pos_{unc['b']}",
                                               bin_um=pc.bin_um(cfg["bin_size"]), n_perm=unc["n_perm"],
                                               min_pos=unc["min_pos"], reference_sample=ref, seed=seed, samples=samples)
                tables["unc"] = W.csv(summ, "Figure4i_observed_to_null_summary.csv", "Distance to nearest partner",
                                      "Observed / random-placement mean distance per section, 95% null interval "
                                      "(Figure 4i).")
                W.csv(dist, "Figure4i_distance_distributions_observed_and_null.csv", "Distance distributions",
                      "Distances observed and under the random-placement null (counts per replicate).")
            done.add("uncoupling")
        if begin("proximity"):
            prox = cfg["proximity"]
            tables["prox"] = []
            for pair in prox["pairs"]:
                src, dst, pmode = pair["from"], pair["to"], pair["mode"]
                need = ("program" in a.obs) if pmode == "assignment" else (f"pos_{src}" in a.obs and f"pos_{dst}" in a.obs)
                if not need or src not in cfg["signatures"] or dst not in cfg["signatures"]:
                    _warn(progress, f"proximity {src} -> {dst}: the {pmode} calls are not available; table skipped")
                    continue
                t = W.csv(proximity_table(a, src, dst, pmode, prox["radius_um"], prox["exclude_self"], samples=samples,
                                          labels=names_by_id), f"{src}_to_{dst}_proximity.csv",
                          f"Proximity {src} to {dst}",
                          f"Fraction of {src} bins with a {dst} bin within {prox['radius_um']:g} um ({pmode} calls), "
                          f"against all ROI bins.")
                tables["prox"].append((src, dst, pmode, t))
            done.add("proximity")
        if begin("moran"):
            moran, cons = morans_i(a, cfg["moran"]["n_genes"], cfg["moran"]["n_consensus"], samples, progress)
            W.csv(moran, "morans_I_per_section.csv", "Moran's I per section", "Moran's I of the highly variable genes.")
            W.csv(cons, "FigureS2_spatially_variable_genes_morans_I.csv", "Spatially variable genes",
                  "Genes in the top Moran's I of every section (Supplementary Figure S2).")
            done.add("moran")
        if begin("sensitivity"):
            sen = cfg["sensitivity"]
            W.csv(geneset_sensitivity(a, sco["n_depth_reps"], sen["variants"], sen["signatures"], sen["competitors"],
                                      mode=mode, n_perm=sco["n_perm"], n_perm_resampled=sco["n_perm_resampled"],
                                      target_sum=target_sum, samples=samples, **call_kw),
                  "gene_set_sensitivity_three_gene_sets.csv", "Gene-set sensitivity",
                  "Bone-lineage positivity under three gene sets, native and depth-matched (Methods).")
            done.add("sensitivity")
        # ---- 11b. sensitivity analyses of the second review (bone_sensitivity/) ---- [APP]
        if begin("sensitivity_analyses", ("assignment",)):
            tables["sens"] = _sensitivity_analyses(a, cfg, W, progress, done)
            done.add("sensitivity_analyses")

        # ---- 12. Table 1, per-bin table ------------------------------------------- [RUN]
        begin("tables")
        W.csv(table1_signatures(cfg["signatures"]), "Table1_gene_signatures.csv", "Table 1: gene signatures",
              "Signature genes, genes counted toward the two-gene rule, excluded genes.")
        W.csv(per_bin_table(a, cfg), f"per_bin_signature_calls_{B}.csv.gz", "Per-bin calls",
              "One row per bin: position (um), region of interest, assignment, tissue class, positivity.")
        done.add("tables")
        if begin("write_h5ad"):
            a.write(out / f"joint_{B}.h5ad")                                   # [SHAMIM] "Save final processed data"
            W.add(f"joint_{B}.h5ad", "object", "Joint AnnData object",
                  "All sections, counts and log-normalized layers, scores, calls, ROI, domains and assignments.")
            done.add("write_h5ad")

        # ---- 13. optional: project the labels onto the QC bin size (Methods) ------ [RUN]
        if begin("project_fine_bins"):
            fine, _, _ = load_sections(cfg, cfg["qc_bin_size"], progress=progress)
            f = ad.concat(fine, label="sample_key", index_unique="-", join="outer")
            f.obs["sample"] = f.obs["sample"].astype(str)
            f = map_domains_from_coarser(f, a, pc.bin_um(cfg["qc_bin_size"]), pc.bin_um(cfg["bin_size"]))
            cols = [c for c in ("sample", "array_row", "array_col", "leiden", "program") if c in f.obs]
            W.csv(f.obs[cols], f"labels_{Q}_projected.csv", f"Labels projected onto {Q} bins",
                  f"Leiden domain and assignment of each {Q} bin's {B} parent.", index=True)
            del f, fine
            done.add("project_fine_bins")

        begin("summary")
        summary = _summary(cfg, a, loaded, tables, progress)
        (out / "summary.json").write_text(json.dumps(summary, indent=2, default=_json_default), encoding="utf-8")
        W.add("summary.json", "data", "Summary", "Key numbers of this run.")
        print(f"done -> {out}")
        return a
    finally:
        W.write()


def _sensitivity_analyses(a, cfg, W, progress, done):
    """[APP] The sensitivity analyses of the second review (bone_sensitivity/; CHANGES.md): first the cross-captured A549
    lesion flag (obs lesion_umi_pooled, in_lesion), then each module whose inputs exist. A module that fails becomes a
    warning; the other modules and every manuscript output are unaffected. Returns {file name: table}."""
    from bone_sensitivity import clustering, common, compartments, positivity, proximity, uncoupling
    les = common.settings(cfg, "lesion")
    pooled, flag, found = common.lesion_signal(a, les["genes"], les["pool_k"], les["min_pooled_umi"])
    a.obs["lesion_umi_pooled"] = pooled                # all 0 / False without the genes (the modules then need not
    a.obs["in_lesion"] = flag                          # recompute and re-warn)
    if not found:
        _warn(progress, f"sensitivity analyses: none of the lesion genes ({', '.join(les['genes'])}) is in the data, "
                        f"so no bin is flagged as lesion")
    modules = [("positivity", positivity, ("positivity",)), ("proximity", proximity, ()),
               ("uncoupling", uncoupling, ()), ("compartments", compartments, ()),
               ("clustering", clustering, ("leiden", "domain_matrix"))]
    written = {}
    for name, module, needs in modules:
        missing = [n for n in needs if n not in done]
        if missing:
            _warn(progress, f"sensitivity analysis {name} skipped: needs the step(s) {', '.join(missing)}")
            continue
        print(f"[sensitivity] {name}")
        try:
            outputs = module.run(a, cfg, progress)
        except Exception as exc:                       # an extra analysis must not stop the manuscript outputs
            _warn(progress, f"sensitivity analysis {name} failed: {type(exc).__name__}: {exc}")
            print(traceback.format_exc())
            continue
        for o in outputs:
            written[o.file] = W.csv(o.df, o.file, o.title, o.description)
    return written


def _sensitivity_metrics(cfg, sens, add):
    """[APP] summary.json entries of the sensitivity analyses (group 'Sensitivity analyses')."""
    G = "Sensitivity analyses"
    sl = {s["id"]: s.get("label") or s["id"] for s in cfg["samples"]}
    lbl = lambda n: pc.signature_label(cfg, n)                                   # noqa: E731
    t = sens.get("sensitivity_Figure4i_summary.csv")
    if t is not None and len(t):
        calls = cfg["uncoupling"]["calls"]
        for r in t[t["call_set"] == calls].to_dict("records"):
            if _num(r.get("ratio_mean")) is None:
                continue
            detail = (f"{calls} calls; range {_num(r['ratio_min'], 3)}–{_num(r['ratio_max'], 3)} over "
                      f"{r['n_draws']} draw(s)")
            if _num(r.get("ratio_of_ratios_mean")) is not None and r["sample"] != cfg["reference_sample"]:
                detail += (f"; ratio of ratios vs reference {_num(r['ratio_of_ratios_mean'], 3)} "
                           f"({_num(r['ratio_of_ratios_min'], 3)}–{_num(r['ratio_of_ratios_max'], 3)})")
            add(G, f"Figure 4i, {sl.get(r['sample'], r['sample'])}: observed / null, {r['null_model']}",
                r["ratio_mean"], "", detail)
    t = sens.get("sensitivity_proximity.csv")
    if t is not None and len(t):
        for r in t[t["status"] == "estimable"].to_dict("records"):
            ac = r.get("assignment_column")
            col = f" ({ac})" if isinstance(ac, str) and ac else ""
            add(G, f"{sl.get(r['sample'], r['sample'])}: {lbl(r['source'])} within {float(r['radius_um']):g} µm of "
                   f"{lbl(r['target'])}{col}, minus depth-matched expectation", r["diff_pp"], "pp",
                f"observed {_num(r['pct_source_within'], 1)}% vs depth-matched "
                f"{_num(r['expected_depth_matched_pct'], 1)}% (95% block CI of the difference "
                f"{_num(r['diff_ci_lo'], 1)} to {_num(r['diff_ci_hi'], 1)}); non-target baseline "
                f"{_num(r['pct_nontarget_within'], 1)}%, all ROI bins {_num(r['pct_all_within'], 1)}%")
    t = sens.get("sensitivity_positivity_folds.csv")
    if t is not None and len(t):
        shown = [s for p in cfg["figures"]["positivity_panels"] for s in p.get("signatures", [])]
        for (sig, smp), g in t[t["signature"].isin(shown)].groupby(["signature", "sample"], sort=False):
            ms = g[g["manuscript"].astype(str).str.lower().isin(["true", "1"]) & (g["depth"] == "native")]
            f = pd.to_numeric(g["fold"], errors="coerce")
            if ms.empty or not np.isfinite(f).any():
                continue
            m = ms.iloc[0]
            lo, hi = _num(m.get("block_fold_ci_lo"), 2), _num(m.get("block_fold_ci_hi"), 2)
            ci = (f"; spatial block-bootstrap 95% CI {lo if lo is not None else 'n/a'}–"
                  f"{hi if hi is not None else 'unbounded'}") if lo is not None or hi is not None else ""
            add(G, f"{lbl(sig)}: {sl.get(smp, smp)} / reference fold, native depth", m["fold"], "×",
                f"range {_num(f.min(), 2)}–{_num(f.max(), 2)} over {int(np.isfinite(f).sum())} region, rule and "
                f"depth-matching variants{ci}")
    t = sens.get("sensitivity_Figure4h_diagnostics.csv")
    if t is not None and len(t):
        r = t[t["sample"].astype(str) == "all"]
        if len(r):
            r = r.iloc[0]
            add(G, "Figure 4h: growth-plate threshold of the split score", r["threshold"], "",
                f"{_num(r['pct_bins_above_threshold'], 2)}% of bins above it"
                + ("; the threshold is 0, so the rule is 'at least one UMI of a split gene'"
                   if str(r["threshold_is_zero"]).lower() in ("true", "1") else ""))
    t = sens.get("sensitivity_Figure4h_variants.csv")
    if t is not None and len(t):
        for r in t[t["status"].astype(str) != "not configured"].to_dict("records"):
            add(G, f"Figure 4h, {r['variant']}: genes significant (growth-plate / tumor-associated)",
                r.get("n_significant_growth_plate"), "", f"{r['status']}; {r.get('n_group')} vs {r.get('n_rest')} "
                f"{r.get('units')}; tumor-associated significant: {r.get('n_significant_tumor_associated')}")
    t = sens.get("sensitivity_leiden_harmony.csv")
    if t is not None and len(t):
        for r in t[t["variant"] == "no_harmony"].to_dict("records"):
            add(G, "Leiden without Harmony: adjusted Rand index vs the Harmony partition", r["ari_vs_harmony"], "",
                f"{r['n_clusters']} raw clusters; same-section kNN edges {_num(r['same_section_edge_fraction'], 3)}")


def _json_default(x):
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return None if np.isnan(x) else float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return str(x)


def _num(x, digits=4):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if not np.isfinite(v) else round(v, digits)


def _summary(cfg, a, loaded, tables, progress):
    """[APP] summary.json: the key numbers of the run, for the web application."""
    sl = {s["id"]: s.get("label") or s["id"] for s in cfg["samples"]}
    target, med = depth_match_target(a)
    obs = a.obs
    samples = []
    for s in cfg["samples"]:
        m = (obs["sample"].astype(str) == s["id"]).to_numpy()
        n_qc = int(m.sum())
        samples.append({"id": s["id"], "label": sl[s["id"]], "condition": s.get("condition"), "color": s.get("color"),
                        "n_bins_loaded": int(loaded.get(s["id"], 0)), "n_bins_qc": n_qc,
                        "pct_retained": _num(100 * n_qc / loaded[s["id"]], 2) if loaded.get(s["id"]) else None,
                        "median_umi": _num(obs.loc[m, "total_counts"].median(), 1),
                        "median_genes": _num(obs.loc[m, "n_genes_by_counts"].median(), 1),
                        "n_bins_roi": int(_as_bool(obs["in_marrow_roi"])[m].sum()),
                        "microns_per_pixel": _num(microns_per_pixel(a, s["id"]), 6)})
    metrics = []

    def add(group, lab, value, unit="", detail=""):
        metrics.append({"group": group, "label": lab, "value": _num(value), "unit": unit, "detail": detail})

    for s in samples:
        add("Sections", f"{s['label']}: bins after QC", s["n_bins_qc"], "",
            f"{s['pct_retained']}% of {s['n_bins_loaded']:,} loaded; median {s['median_umi']} UMI per bin")
        add("Sections", f"{s['label']}: bins in the region of interest", s["n_bins_roi"])
    lbl = lambda n: pc.signature_label(cfg, n)                                   # noqa: E731
    if "native" in tables:
        for r in tables["native"].itertuples():
            add("Positivity, native depth (ROI)", f"{lbl(r.signature)}: {sl.get(r.sample, r.sample)}", r.pct_positive,
                "%", f"{r.n_positive:,} of {r.n_bins:,} bins; Wilson 95% CI {r.ci_lo}–{r.ci_hi}%")
    if "dm" in tables:
        for r in tables["dm"].itertuples():
            add("Positivity, depth-matched (ROI)", f"{lbl(r.signature)}: {sl.get(r.sample, r.sample)}",
                r.pct_positive_mean, "%", f"mean ± SD {_num(r.pct_positive_mean)} ± {_num(r.pct_positive_sd)}% over "
                                          f"{cfg['scoring']['n_depth_reps']} replicates; {cfg['scoring']['depth_matching']} matching")
    for key, title in (("comp", "native"), ("comp_dm", "depth-matched")):
        if key in tables and len(tables[key]):
            for s, sub in tables[key].groupby("sample", sort=False):
                add("Tissue assignment", f"{sl.get(s, s)}: ambiguous bins ({title})", sub["sample_ambiguous_pct"].iloc[0], "%")
    if "unc" in tables:
        for r in tables["unc"].itertuples():
            if r.status == "estimable":
                add("Distance to nearest partner", f"{sl.get(r.sample, r.sample)}: observed / null mean",
                    r.ratio_obs_null_mean, "", f"95% null interval {_num(r.ratio_obs_null_lo, 3)}–"
                                               f"{_num(r.ratio_obs_null_hi, 3)}; P(null ≥ observed) = "
                                               f"{_num(r.p_greater_than_null_mean, 3)}")
                if pd.notna(r.ratio_of_ratios_vs_reference) and r.sample != cfg["reference_sample"]:
                    add("Distance to nearest partner", f"{sl.get(r.sample, r.sample)}: ratio of ratios vs reference",
                        r.ratio_of_ratios_vs_reference)
            else:
                add("Distance to nearest partner", f"{sl.get(r.sample, r.sample)}: {r.status}", None)
    for src, dst, pmode, t in tables.get("prox", []):
        R = f"{cfg['proximity']['radius_um']:g}"
        for r in t.to_dict("records"):
            col = f"percent_{src}_within_{R}um_of_{dst}" if pmode == "assignment" else f"percent_{src}_within_{R}um"
            base = f"percent_all_bins_within_{R}um_of_{dst}" if pmode == "assignment" else f"percent_all_bins_within_{R}um"
            add("Proximity", f"{r['section']}: {lbl(src)} within {R} µm of {lbl(dst)}", r.get(col), "%",
                f"against {r.get(base)}% of all ROI bins")
    if "unresolved_pct" in tables:
        add("Leiden domains", "Bins in unresolved domains", tables["unresolved_pct"], "%")
        for r in tables.get("unresolved", pd.DataFrame(columns=["sample"])).to_dict("records"):
            if r["sample"] != "all":
                roi = _num(r["percent_roi"], 1)
                add("Leiden domains", f"{sl.get(r['sample'], r['sample'])}: bins in unresolved domains", r["percent"],
                    "%", f"{int(r['unresolved_bins']):,} of {int(r['all_bins']):,} bins"
                         + (f"; region of interest {roi}%" if roi is not None else ""))
        add("Leiden domains", "Domains", obs["leiden"].nunique())
    if tables.get("ari") is not None:
        add("Leiden domains", "Median adjusted Rand index across seeds", tables["ari"])
    if "receptor" in tables:
        rec = tables["receptor"]
        for r in rec[rec["fold_dm"].notna()].itertuples():
            add("Host receptor panel", f"{r.gene} ({r.role}): {sl.get(r.sample, r.sample)} / reference", r.fold_dm,
                "×", f"{_num(r.pct_detected_dm, 3)}% of ROI bins detecting after depth matching")
    if tables.get("sens"):
        try:
            _sensitivity_metrics(cfg, tables["sens"], add)
        except Exception as exc:                       # a summary line must not lose the run
            _warn(progress, f"summary of the sensitivity analyses incomplete: {type(exc).__name__}: {exc}")
    return {"version": 1, "bin_size": cfg["bin_size"], "reference_sample": cfg["reference_sample"],
            "depth_match_target_umi": int(target), "samples": samples, "metrics": metrics,
            "warnings": list(progress.warnings), "skipped": list(progress.skipped), "timings": list(progress.timings),
            "environment": _environment(), "finished": datetime.now().isoformat(timespec="seconds")}


def run_analysis(root, out, roi_masks=None, mask_frame="hires", do_moran=True, do_sensitivity=False, project_8um=False):
    """[RUN] The manuscript analysis. [APP] A wrapper over run_analysis_from_config with the default configuration."""
    cfg = pc.default_config()
    cfg["data_root"] = str(root)
    cfg["steps"].update(moran=do_moran, sensitivity=do_sensitivity, project_fine_bins=project_8um)
    if roi_masks:
        cfg["roi"].update(use_drawn=True, masks={k: str(v) for k, v in roi_masks.items()}, mask_frame=mask_frame)
    return run_analysis_from_config(cfg, out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=None, help="folder containing SR007144_LL-2-1/ and SR007144_RL-1-5/ "
                                                 "(default: the configuration's data_root, or .)")
    ap.add_argument("--out", default="Manuscript/code/results")
    ap.add_argument("--config", default=None, help="JSON configuration (pipeline_config.py); default: the manuscript")
    ap.add_argument("--roi-mask", action="append", default=None,
                    help="hand-drawn ROI as SAMPLE=PATH (.geojson on tissue_hires_image.png); repeat per sample")
    ap.add_argument("--mask-frame", default=None, choices=["hires", "fullres", "lowres"])
    ap.add_argument("--skip-moran", action="store_true", help="skip Supplementary Figure S2 (slow)")
    ap.add_argument("--geneset-sensitivity", action="store_true", help="run the three-gene-set check of Methods")
    ap.add_argument("--project-8um", action="store_true", help="write the 16 um labels projected onto 8 um bins")
    ap.add_argument("--events", action="store_true", help="also print machine-readable @@EVENT progress lines")
    args = ap.parse_args()
    config = pc.load_config(args.config) if args.config else pc.default_config()
    if args.root is not None or not args.config:
        config["data_root"] = args.root or "."
    if args.skip_moran:
        config["steps"]["moran"] = False
    if args.geneset_sensitivity:
        config["steps"]["sensitivity"] = True
    if args.project_8um:
        config["steps"]["project_fine_bins"] = True
    if args.roi_mask:
        config["roi"].update(use_drawn=True, masks=dict(v.split("=", 1) for v in args.roi_mask))
    if args.mask_frame:
        config["roi"]["mask_frame"] = args.mask_frame
    prog = pc.Progress(events=args.events)
    try:
        run_analysis_from_config(config, args.out, prog)
        prog.done()
    except Exception as exc:
        prog.error(str(exc), traceback.format_exc())
        raise
