"""
Figure 4h (osteoclast compartments) under other compartment definitions and tests (second review F7-F10, F26, G1.2).

The manuscript calls the <signature>-positive bins whose <split_by> score is above the <percentile>th percentile of all
bins the growth-plate compartment and the other positive bins of tumor sections the tumor-associated compartment
(bone_analysis.split_osteoclast_compartments), then compares them gene by gene with a Wilcoxon rank-sum test over bins
(compare_osteoclast_states). The review found that the percentile is exactly 0 when most bins score <= 0, so the rule
becomes 'any split_by gene detected'; that the selected bins are not at the growth plate; that they are deeper than the
tumor-associated bins; that contiguous bins of one section are treated as independent; that the split_by genes are
higher in the growth-plate group by construction; and that many tested bins lie outside the region of interest.
This module writes
  diagnostics  the effective threshold and the bins it selects, per section: depth, region of interest, split_by UMIs;
  regions      an anatomical proxy of the growth plate: 8-connected regions of split_by-positive bins (pooled calls) on
               the array grid with at least min_region_bins bins, and the positive bins within anatomical_distance_um;
  variants, de, top_genes   the test under each variant:
    manuscript        the bins and test of Figure 4h (reproduces Figure4h_osteoclast_compartment_DE.csv)
    roi_only          both groups restricted to the marrow region of interest
    depth_adjusted    the deeper group downsampled bin by bin to the other group's UMI distribution (rank-quantile
                      targets), n_depth_reps replicates, re-log-normalized
    block_pseudobulk  units = spatial blocks x compartment per section (summed counts), not bins; the deeper
                      compartment's units downsampled to the other's unit-UMI distribution (the dense compartment
                      sums many more bins per block), n_depth_reps replicates
    anatomical        growth plate = positive bins near a split_by-dense region
    drawn             growth plate = positive bins inside hand-drawn polygons (growth_plate_masks)
The anatomical and drawn variants take both compartments from tumor sections whatever de_sections says: with
de_sections = all, a control limb's intact growth plate would enter the growth-plate side and the test would compare
sections and conditions instead of compartments within the tumor-bearing section.
selection_gene flags the genes of the split_by extended list: they choose the growth-plate bins.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import pipeline_config as pc
from bone_sensitivity import common

GROUPS = ("growth_plate", "tumor_associated")
VARIANTS = ("manuscript", "roi_only", "depth_adjusted", "block_pseudobulk", "anatomical", "drawn")
ALPHA = 0.05
TOP_N_DE = 50                      # genes per variant and group kept in the DE table (plus manuscript-significant ones)
FIGURE_TOP_N = 8                   # genes per compartment drawn in Figure 4h (bone_figures.DE_TOP_N)
DIAG_CSV = "sensitivity_Figure4h_diagnostics.csv"
REGIONS_CSV = "sensitivity_Figure4h_regions.csv"
VARIANTS_CSV = "sensitivity_Figure4h_variants.csv"
DE_CSV = "sensitivity_Figure4h_de.csv"
TOP_CSV = "sensitivity_Figure4h_top_genes.csv"
DIAG_COLUMNS = ["sample", "threshold", "threshold_is_zero", "pct_bins_above_threshold", "n_bins_with_split_gene_umi",
                "n_growth_plate", "n_tumor_associated", "n_non_growth_plate_control", "n_growth_plate_outside_roi",
                "n_tumor_associated_outside_roi", "median_umi_growth_plate", "median_umi_tumor_associated",
                "n_growth_plate_one_split_umi", "n_growth_plate_without_core_gene", "n_bins", "n_bins_above_threshold",
                "pct_bins_score_zero"]
REGION_COLUMNS = ["sample", "n_regions", "region_sizes", "n_osteoclast_bins_near_regions",
                  "n_osteoclast_bins_near_regions_manuscript_growth_plate", "calls", "n_split_positive_bins",
                  "n_region_bins", "n_osteoclast_positive", "n_manuscript_growth_plate"]
VARIANT_COLUMNS = ["variant", "status", "units", "n_group", "n_rest", "median_umi_group", "median_umi_rest",
                   "n_significant_growth_plate", "n_significant_tumor_associated", "note",
                   "median_umi_group_unadjusted", "median_umi_rest_unadjusted"]
DE_COLUMNS = ["variant", "group", "names", "scores", "pvals_adj", "log2fc_mean_cp10k", "pct_group", "pct_rest",
              "selection_gene", "frac_reps_significant"]


# ---- compartment rule ------------------------------------------------------------------------------------------------

def effective_threshold(score, pct):
    """The manuscript split of the split_by score: (pct-th percentile of all scores, bins strictly above it, whether
    the percentile is exactly 0). When at least pct % of bins score <= 0 the percentile is 0 and the rule is
    'score > 0', i.e. some split_by gene detected above the background genes."""
    s = np.asarray(score, float)
    if not s.size or np.isnan(s).all():
        return np.nan, np.zeros(s.shape, dtype=bool), False
    thr = float(np.nanpercentile(s, pct))
    return thr, s > thr, bool(thr == 0.0)


def manuscript_compartments(a, cp):
    """oc_compartment of every bin by the manuscript rule: bone_analysis.split_osteoclast_compartments on a light
    object holding only the obs columns it reads, so `a` is not modified."""
    import anndata as ad
    import bone_analysis as ba
    b = ad.AnnData(obs=a.obs[[f"score_{cp['split_by']}", f"pos_{cp['signature']}", "condition"]].copy())
    ba.split_osteoclast_compartments(b, cp["percentile"], cp["signature"], cp["split_by"], cp["tumor_conditions"])
    return b.obs["oc_compartment"].astype(str).to_numpy()


def gene_umis(counts, var_names, genes):
    """UMIs of `genes` (those present) summed per bin."""
    idx = pd.Index(var_names).get_indexer([g for g in genes if g in set(var_names)])
    if not len(idx):
        return np.zeros(counts.shape[0])
    return np.asarray(counts[:, idx].sum(axis=1), dtype=float).ravel()


# ---- anatomical proxy ------------------------------------------------------------------------------------------------

def connected_regions(rows, cols, mask, min_bins):
    """8-connected components of the bins in `mask` on the integer array grid (rows, cols of one section). Returns the
    region of every bin (0 = largest kept region, ...; -1 = not in a region of at least min_bins bins) and the sizes of
    the kept regions, largest first."""
    from scipy import ndimage
    rows, cols, mask = np.asarray(rows, np.int64), np.asarray(cols, np.int64), np.asarray(mask, bool)
    label = np.full(len(rows), -1, dtype=np.int64)
    if not mask.any():
        return label, []
    r, c = rows[mask] - rows[mask].min(), cols[mask] - cols[mask].min()
    grid = np.zeros((r.max() + 1, c.max() + 1), dtype=bool)
    grid[r, c] = True
    comp, _ = ndimage.label(grid, structure=np.ones((3, 3), dtype=int))
    per_bin = comp[r, c]
    sizes = np.bincount(per_bin)
    kept = np.flatnonzero(sizes >= max(int(min_bins), 1))
    kept = kept[kept > 0]
    kept = kept[np.argsort(-sizes[kept], kind="stable")]
    relabel = np.full(len(sizes), -1, dtype=np.int64)
    relabel[kept] = np.arange(len(kept))
    label[mask] = relabel[per_bin]
    return label, [int(sizes[k]) for k in kept]


def near_regions(xy_um, region, query, distance_um):
    """Bins of `query` whose nearest `region` bin (one section, positions in um) lies within distance_um."""
    from scipy.spatial import cKDTree
    xy, region, query = np.asarray(xy_um, float), np.asarray(region, bool), np.asarray(query, bool)
    out = np.zeros(len(xy), dtype=bool)
    if region.any() and query.any():
        d, _ = cKDTree(xy[region]).query(xy[query], k=1)
        out[query] = d <= float(distance_um) * (1 + 1e-9)
    return out


def anatomical_growth_plate(a, oc, manuscript_gp, split_by, distance_um, min_bins, bin_um, sections, progress=None):
    """Anatomical proxy of the growth plate. Per section, split_by-dense regions are the 8-connected components (array
    grid) of bins with pooled_pos_<split_by> (pos_<split_by> when the pooled calls are absent) of at least min_bins
    bins; anatomical growth-plate bins are the `oc` bins within distance_um (array um) of a region bin.
    Returns (anatomical growth-plate flag per bin, regions table, calls column used or None)."""
    calls = next((c for c in (f"pooled_pos_{split_by}", f"pos_{split_by}") if c in a.obs), None)
    if calls is None:
        common.warn(progress, f"Figure 4h sensitivity: no {split_by} positivity calls, so no {split_by}-dense region")
    elif calls != f"pooled_pos_{split_by}":
        common.warn(progress, f"Figure 4h sensitivity: pooled_pos_{split_by} is absent (tissue assignment not run); "
                              f"the {split_by}-dense regions use the native calls {calls}")
    pos = common.as_bool(a.obs[calls]) if calls else np.zeros(a.n_obs, dtype=bool)
    smp = a.obs["sample"].astype(str).to_numpy()
    rows_all, cols_all = a.obs["array_row"].to_numpy(), a.obs["array_col"].to_numpy()
    xy = common.array_xy_um(a, bin_um)
    near = np.zeros(a.n_obs, dtype=bool)
    rows = []
    for s in sections:
        m = np.flatnonzero(smp == s)
        label, sizes = connected_regions(rows_all[m], cols_all[m], pos[m], min_bins)
        near[m] = near_regions(xy[m], label >= 0, oc[m], distance_um)
        rows.append({"sample": s, "n_regions": len(sizes), "region_sizes": ";".join(map(str, sizes)),
                     "n_osteoclast_bins_near_regions": int(near[m].sum()),
                     "n_osteoclast_bins_near_regions_manuscript_growth_plate": int((near[m] & manuscript_gp[m]).sum()),
                     "calls": calls or "", "n_split_positive_bins": int(pos[m].sum()),
                     "n_region_bins": int((label >= 0).sum()), "n_osteoclast_positive": int(oc[m].sum()),
                     "n_manuscript_growth_plate": int(manuscript_gp[m].sum())})
    return near, pd.DataFrame(rows, columns=REGION_COLUMNS), calls


def drawn_growth_plate(a, masks, frame, data_root):
    """Bins inside the hand-drawn growth-plate polygons (bone_pending_masks.rasterise_drawn_roi; masks: sample id ->
    GeoJSON path, text or object, relative paths against data_root). A section without polygons has no growth-plate
    bin. Returns (flag per bin, sections with polygons)."""
    from bone_pending_masks import rasterise_drawn_roi
    resolved = {}
    for sid, src in (masks or {}).items():
        if isinstance(src, dict) or (isinstance(src, str) and src.lstrip().startswith("{")):
            resolved[sid] = src
            continue
        path = Path(str(src)).expanduser()
        if not path.is_absolute() and not path.exists():
            path = Path(str(data_root)) / path
        if not path.is_file():
            raise FileNotFoundError(f"growth-plate polygons of {sid}: {path} not found")
        resolved[sid] = path
    inside, unrestricted = rasterise_drawn_roi(a, resolved, frame)
    smp = a.obs["sample"].astype(str).to_numpy()
    inside = inside & ~np.isin(smp, list(unrestricted))
    return inside, [s for s in pd.unique(smp) if s not in set(unrestricted)]


# ---- depth adjustment and pseudobulk ---------------------------------------------------------------------------------

def quantile_targets(totals, reference):
    """UMI target of each bin of a deeper group: the reference group's UMI count at the bin's rank quantile
    ((rank - 0.5) / n, ties broken by order), never above the bin's own total (bone_analysis.depth_match_targets,
    quantile mode). The targets then follow the reference distribution."""
    t, ref = np.asarray(totals, float), np.sort(np.asarray(reference, float))
    if not t.size or not ref.size:
        return t.astype(np.int64)
    q = (pd.Series(t).rank(method="first").to_numpy() - 0.5) / len(t)
    return np.minimum(t, np.floor(np.quantile(ref, q))).astype(np.int64)


def depth_adjusted_targets(totals, labels):
    """Targets of the depth-adjusted variant: the group with the higher median UMI per bin is downsampled to the other
    group's distribution (quantile_targets); the other group keeps its counts. Returns (targets, deeper group)."""
    t, labels = np.asarray(totals, float), np.asarray(labels)
    g, r = (labels == GROUPS[0]), (labels == GROUPS[1])
    deeper = GROUPS[0] if np.median(t[g]) >= np.median(t[r]) else GROUPS[1]
    deep, shallow = (g, r) if deeper == GROUPS[0] else (r, g)
    targets = t.astype(np.int64)
    targets[deep] = quantile_targets(t[deep], t[shallow])
    return targets, deeper


def block_units(sections, xy_um, labels, block_um):
    """Pseudobulk units: per section, square blocks of side block_um (common.spatial_blocks on array um) crossed with
    the compartment. Returns (unit of every bin, compartment of every unit, section of every unit)."""
    sections, labels = np.asarray(sections).astype(str), np.asarray(labels).astype(str)
    xy = np.asarray(xy_um, float)
    block = np.zeros(len(sections), dtype=np.int64)
    for s in pd.unique(sections):
        m = sections == s
        block[m] = common.spatial_blocks(xy[m], block_um)
    unit, uniq = pd.MultiIndex.from_arrays([sections, block, labels]).factorize()
    return unit, np.asarray(uniq.get_level_values(2)), np.asarray(uniq.get_level_values(0))


def sum_units(counts, unit, n_units):
    """Counts summed per unit (sparse indicator product)."""
    ind = sp.csr_matrix((np.ones(len(unit), np.float32), (unit, np.arange(len(unit)))), shape=(n_units, len(unit)))
    return (ind @ sp.csr_matrix(counts, dtype=np.float32)).tocsr()


def log_normalize(counts, target_sum):
    """log1p of the counts scaled to target_sum per row (scanpy normalize_total + log1p), float32 CSR."""
    X = sp.csr_matrix(counts, dtype=np.float32, copy=True)
    tot = np.asarray(X.sum(axis=1), dtype=float).ravel()
    scale = np.divide(float(target_sum), tot, out=np.zeros_like(tot), where=tot > 0)
    X = (sp.diags(scale.astype(np.float32)) @ X).tocsr()
    X.data = np.log1p(X.data)
    return X


# ---- differential expression -----------------------------------------------------------------------------------------

def wilcoxon(X, labels, var_names, obs_names=None):
    """scanpy rank_genes_groups (Wilcoxon rank-sum, Benjamini-Hochberg) of the growth_plate against the
    tumor_associated rows of X, built as bone_analysis.compare_osteoclast_states: both groups, every gene."""
    import anndata as ad
    import scanpy as sc
    labels = np.asarray(labels).astype(str)
    index = pd.Index(obs_names if obs_names is not None else [f"u{i}" for i in range(len(labels))]).astype(str)
    sub = ad.AnnData(X=X, var=pd.DataFrame(index=pd.Index(var_names).copy()),
                     obs=pd.DataFrame({"oc_compartment": pd.Categorical(labels, categories=list(GROUPS))}, index=index))
    sc.tl.rank_genes_groups(sub, groupby="oc_compartment", method="wilcoxon", use_raw=False)
    de = sc.get.rank_genes_groups_df(sub, group=None)
    de["group"] = de["group"].astype(str)
    return de


def group_stats(counts, labels, target_sum):
    """Per compartment: (mean CP10k, percentage of units with >= 1 UMI) of every gene."""
    X = sp.csr_matrix(counts, dtype=np.float64)
    tot = np.asarray(X.sum(axis=1)).ravel()
    N = sp.diags(np.divide(float(target_sum), tot, out=np.zeros_like(tot), where=tot > 0)) @ X
    labels = np.asarray(labels).astype(str)
    out = {}
    for g in GROUPS:
        m = labels == g
        n = max(int(m.sum()), 1)
        out[g] = (np.asarray(N[m].sum(axis=0)).ravel() / n, 100 * np.asarray((X[m] > 0).sum(axis=0)).ravel() / n)
    return out


def de_table(X, counts, labels, var_names, target_sum, selection, obs_names=None):
    """Wilcoxon table of one variant (group, names, scores, pvals_adj) with log2fc_mean_cp10k = log2((mean CP10k of
    the group + 1) / (mean CP10k of the other + 1)), pct_group, pct_rest (% of units detecting the gene) and
    selection_gene."""
    de = wilcoxon(X, labels, var_names, obs_names)
    stats = group_stats(counts, labels, target_sum)
    gi = pd.Index(var_names).get_indexer(de["names"])
    grp = de["group"].to_numpy()
    mean_g, mean_o, pct_g, pct_o = (np.full(len(de), np.nan) for _ in range(4))
    for g, o in (GROUPS, GROUPS[::-1]):
        m = grp == g
        mean_g[m], pct_g[m] = stats[g][0][gi[m]], stats[g][1][gi[m]]
        mean_o[m], pct_o[m] = stats[o][0][gi[m]], stats[o][1][gi[m]]
    return pd.DataFrame({"group": grp, "names": de["names"].astype(str).to_numpy(),
                         "scores": de["scores"].to_numpy(float), "pvals_adj": de["pvals_adj"].to_numpy(float),
                         "log2fc_mean_cp10k": np.log2((mean_g + 1) / (mean_o + 1)), "pct_group": pct_g,
                         "pct_rest": pct_o, "selection_gene": de["names"].astype(str).isin(set(selection)).to_numpy(),
                         "frac_reps_significant": np.nan})


def significant(de):
    """Genes significantly higher in the row's compartment: adjusted P < 0.05 and a positive score."""
    return (de["pvals_adj"] < ALPHA) & (de["scores"] > 0)


def aggregate_replicates(reps):
    """Depth-adjusted replicates -> one table: per gene and group the median of every statistic, and the fraction of
    replicates in which the gene is significantly higher in that group."""
    df = pd.concat(reps, ignore_index=True)
    df["_sig"] = significant(df).astype(float)
    agg = (df.groupby(["group", "names"], sort=False)
             .agg(scores=("scores", "median"), pvals_adj=("pvals_adj", "median"),
                  log2fc_mean_cp10k=("log2fc_mean_cp10k", "median"), pct_group=("pct_group", "median"),
                  pct_rest=("pct_rest", "median"), selection_gene=("selection_gene", "first"),
                  frac_reps_significant=("_sig", "mean"))
             .reset_index())
    return agg.sort_values(["group", "scores"], ascending=[True, False], kind="stable", ignore_index=True)


def block_pseudobulk(counts, sections, xy_um, labels, block_um, var_names, target_sum, selection, n_reps, seed):
    """The block_pseudobulk test. Units = block_units of the tested bins with their counts summed (sum_units). The
    compartment whose units are deeper (higher median UMI per unit; the dense compartment sums many more bins per
    block) is downsampled unit by unit to the other compartment's unit-UMI distribution (depth_adjusted_targets on the
    unit totals, bone_analysis.downsample_counts_matrix, seed + replicate): after CP10k and log1p a shallower unit has
    more zeros, so unbalanced units favour the deeper compartment for every lowly expressed gene. Each of the n_reps
    (>= 1) replicates is re-log-normalized and tested between units (de_table); aggregate_replicates combines them.
    Returns (DE table, or None with fewer than two units in a compartment; info: units, median bins and UMI per unit
    before and after the adjustment per compartment, the deeper compartment)."""
    import bone_analysis as ba
    unit, ulab, usec = block_units(sections, xy_um, labels, block_um)
    U = sum_units(counts, unit, len(ulab))
    uumi = np.asarray(U.sum(axis=1)).ravel()
    nbins = np.bincount(unit, minlength=len(ulab))
    info = {"n": {g: int((ulab == g).sum()) for g in GROUPS}, "n_bins": len(unit),
            "bins_per_unit": {g: _median(nbins[ulab == g]) for g in GROUPS},
            "umi": {g: _median(uumi[ulab == g]) for g in GROUPS}, "target": {g: np.nan for g in GROUPS},
            "deeper": None}
    if min(info["n"].values()) < 2:
        return None, info
    targets, info["deeper"] = depth_adjusted_targets(uumi, ulab)
    info["target"] = {g: _median(targets[ulab == g]) for g in GROUPS}
    names = [f"{s}|{u}|{g}" for u, (s, g) in enumerate(zip(usec, ulab))]
    reps = []
    for rep in range(int(n_reps)):
        Ud = ba.downsample_counts_matrix(U, targets, seed + rep)
        reps.append(de_table(log_normalize(Ud, target_sum), Ud, ulab, var_names, target_sum, selection, names))
        del Ud
    return aggregate_replicates(reps), info


# ---- tables ----------------------------------------------------------------------------------------------------------

def diagnostics(sections, smp, thr, above, is_zero, score, comp, roi, umi, split_umi, core_umi):
    """sensitivity_Figure4h_diagnostics.csv: one row per section and one 'all'."""
    rows = []
    for s in list(sections) + ["all"]:
        m = np.ones(len(smp), dtype=bool) if s == "all" else smp == s
        gp, ta = m & (comp == GROUPS[0]), m & (comp == GROUPS[1])
        n = int(m.sum())
        rows.append({
            "sample": s, "threshold": thr, "threshold_is_zero": is_zero,
            "pct_bins_above_threshold": 100 * above[m].mean() if n else np.nan,
            "n_bins_with_split_gene_umi": int((m & (split_umi >= 1)).sum()),
            "n_growth_plate": int(gp.sum()), "n_tumor_associated": int(ta.sum()),
            "n_non_growth_plate_control": int((m & (comp == "non_growth_plate_control")).sum()),
            "n_growth_plate_outside_roi": int((gp & ~roi).sum()) if roi is not None else np.nan,
            "n_tumor_associated_outside_roi": int((ta & ~roi).sum()) if roi is not None else np.nan,
            "median_umi_growth_plate": float(np.median(umi[gp])) if gp.any() else np.nan,
            "median_umi_tumor_associated": float(np.median(umi[ta])) if ta.any() else np.nan,
            "n_growth_plate_one_split_umi": int((gp & (split_umi == 1)).sum()),
            "n_growth_plate_without_core_gene": int((gp & (core_umi == 0)).sum()),
            "n_bins": n, "n_bins_above_threshold": int(above[m].sum()),
            "pct_bins_score_zero": 100 * float((score[m] == 0).mean()) if n else np.nan})
    return pd.DataFrame(rows, columns=DIAG_COLUMNS)


def _row(variant, status, units, n_group=0, n_rest=0, umi_group=np.nan, umi_rest=np.nan, de=None, note="",
         raw_group=None, raw_rest=None):
    """One variants row; median_umi_* after any depth adjustment, *_unadjusted before it (raw_*, default the same)."""
    sig = significant(de) if de is not None else None
    return {"variant": variant, "status": status, "units": units, "n_group": int(n_group), "n_rest": int(n_rest),
            "median_umi_group": umi_group, "median_umi_rest": umi_rest,
            "n_significant_growth_plate": int((sig & (de["group"] == GROUPS[0])).sum()) if sig is not None else np.nan,
            "n_significant_tumor_associated": int((sig & (de["group"] == GROUPS[1])).sum()) if sig is not None
            else np.nan, "note": note,
            "median_umi_group_unadjusted": umi_group if raw_group is None else raw_group,
            "median_umi_rest_unadjusted": umi_rest if raw_rest is None else raw_rest}


def _median(x):
    return float(np.median(x)) if len(x) else np.nan


def _labels(gp, rest):
    """Compartment label per bin ('' = not tested)."""
    lab = np.full(len(gp), "", dtype=object)
    lab[rest] = GROUPS[1]
    lab[gp] = GROUPS[0]
    return lab


def de_subset(de_all, manuscript, n_top=TOP_N_DE):
    """Rows of every variant for the genes in the top n_top by score of any variant and group, or significant (either
    direction) in the manuscript variant."""
    keep = set()
    for de in de_all.values():
        for _, g in de.groupby("group", sort=False):
            keep |= set(g.nlargest(n_top, "scores")["names"])
    if manuscript is not None:
        keep |= set(manuscript.loc[manuscript["pvals_adj"] < ALPHA, "names"])
    parts = [de[de["names"].isin(keep)].assign(variant=v) for v, de in de_all.items()]
    if not parts:
        return pd.DataFrame(columns=DE_COLUMNS)
    return pd.concat(parts, ignore_index=True).reindex(columns=DE_COLUMNS)


def top_genes(de_all, selection, n_top=FIGURE_TOP_N):
    """The genes of Figure 4h (top n_top per compartment by manuscript score, as bone_figures.figure4h) with each
    variant's score, adjusted P and significance in the same direction (NaN where the variant was not estimable)."""
    cols = ["gene", "group", "selection_gene"] + [f"{k}_{v}" for v in VARIANTS
                                                   for k in ("padj", "significant", "score")]
    man = de_all.get("manuscript")
    if man is None:
        return pd.DataFrame(columns=cols)
    man = man[np.isfinite(man["scores"])]
    rows = pd.concat([man[man["group"] == g].nlargest(n_top, "scores") for g in GROUPS], ignore_index=True)
    out = pd.DataFrame({"gene": rows["names"].to_numpy(), "group": rows["group"].to_numpy(),
                        "selection_gene": rows["names"].isin(set(selection)).to_numpy()})
    for v in VARIANTS:
        de = de_all.get(v)
        if de is None:
            out[f"padj_{v}"], out[f"significant_{v}"], out[f"score_{v}"] = np.nan, None, np.nan
            continue
        d = de.set_index(["group", "names"])
        key = pd.MultiIndex.from_arrays([out["group"], out["gene"]])
        padj, score = d["pvals_adj"].reindex(key).to_numpy(), d["scores"].reindex(key).to_numpy()
        out[f"padj_{v}"], out[f"score_{v}"] = padj, score
        out[f"significant_{v}"] = np.where(np.isfinite(padj), (padj < ALPHA) & (score > 0), None)
    return out.reindex(columns=cols)


# ---- run -------------------------------------------------------------------------------------------------------------

def run(a, cfg, progress=None):
    """Figure 4h sensitivity: diagnostics of the manuscript rule, the anatomical proxy, and the compartment test under
    six variants (module docstring). Adds no obs column."""
    import bone_analysis as ba
    cp, st = cfg["compartments"], common.settings(cfg, "compartments")
    sig, split = cp["signature"], cp["split_by"]
    if f"pos_{sig}" not in a.obs or f"score_{split}" not in a.obs:
        common.warn(progress, f"Figure 4h sensitivity needs a positivity call for {sig} and a score for {split}; "
                              f"skipped")
        return []
    seed, target_sum = int(cfg.get("seed", 0)), float(cfg["clustering"]["target_sum"])
    bin_um, block_um = pc.bin_um(cfg["bin_size"]), float(st["pseudobulk_block_um"])
    dist_um, min_bins, n_reps = float(st["anatomical_distance_um"]), int(st["min_region_bins"]), int(st["n_depth_reps"])
    S, P = pc.signature_label(cfg, sig), pc.signature_label(cfg, split)
    tumor_txt = ", ".join(cp["tumor_conditions"])
    sections = common.sample_ids(a, cfg)
    smp = a.obs["sample"].astype(str).to_numpy()
    tumor = a.obs["condition"].astype(str).isin(list(cp["tumor_conditions"])).to_numpy()
    in_sections = tumor if cp["de_sections"] == "tumor" else np.ones(a.n_obs, dtype=bool)
    oc = common.as_bool(a.obs[f"pos_{sig}"])
    roi = common.as_bool(a.obs["in_marrow_roi"]) if "in_marrow_roi" in a.obs else None
    counts, lognorm, var_names = a.layers["counts"], a.layers["lognorm"], a.var_names
    umi = np.asarray(counts.sum(axis=1), dtype=float).ravel()
    split_sig = cfg["signatures"].get(split, {})
    selection = [g for g in split_sig.get("extended", []) if g in set(var_names)]

    # 1. diagnostics of the manuscript rule
    score = a.obs[f"score_{split}"].to_numpy(float)
    thr, above, is_zero = effective_threshold(score, cp["percentile"])
    comp = manuscript_compartments(a, cp)
    split_umi = gene_umis(counts, var_names, split_sig.get("extended", []))
    core_umi = gene_umis(counts, var_names, split_sig.get("core", []))
    diag = diagnostics(sections, smp, thr, above, is_zero, score, comp, roi, umi, split_umi, core_umi)
    if is_zero:
        common.warn(progress, f"Figure 4h: the {cp['percentile']:g}th percentile of the {P} score is exactly 0, so the "
                              f"growth-plate rule is '{P} score > 0' ({int(above.sum()):,} bins; "
                              f"{int((split_umi >= 1).sum()):,} bins have >= 1 {P} UMI)")
    rule = (f"{S}-positive bins with {P} score > {thr:.4g} (the {cp['percentile']:g}th percentile of all bins"
            + (", exactly 0" if is_zero else "") + ")")

    # 2. anatomical proxy
    gp_m = comp == GROUPS[0]
    near, regions, calls = anatomical_growth_plate(a, oc, gp_m, split, dist_um, min_bins, bin_um, sections, progress)
    tumor_secs = [s for s in sections if tumor[smp == s].any()]

    def per_section(flag, tumor_side):
        return ", ".join(f"{s} {int((flag & (smp == s)).sum()):,}" for s in sections
                         if (s in tumor_secs) == tumor_side and (tumor_side or (flag & (smp == s)).any()))

    if cp["de_sections"] == "tumor":
        where = f"both compartments from {tumor_txt} sections"
    else:
        off = per_section(gp_m, False)
        where = (f"de_sections = all: growth-plate bins of every section ({int((gp_m & ~tumor).sum()):,} of "
                 f"{int(gp_m.sum()):,} from non-{tumor_txt} sections" + (f": {off}" if off else "")
                 + f"), tumor-associated bins of {tumor_txt} sections")

    def within_tumor(gp):
        """Note of a variant that takes both compartments from tumor sections: growth-plate bins per section."""
        off = per_section(gp, False)
        return (f"both compartments from {tumor_txt} sections only (whatever de_sections), so that the contrast stays "
                f"within the tumor-bearing sections; growth-plate bins per {tumor_txt} section: "
                f"{per_section(gp, True) or 'none'}" + (f"; left out, not {tumor_txt} sections: {off}" if off else ""))

    # 3. variants
    de_all, rows = {}, []

    def bins_variant(name, lab, note, X=None):
        m = lab != ""
        labels = lab[m].astype(str)
        n = {g: int((labels == g).sum()) for g in GROUPS}
        if min(n.values()) < 2:
            common.warn(progress, f"Figure 4h sensitivity ({name}): {n[GROUPS[0]]} growth-plate and {n[GROUPS[1]]} "
                                  f"tumor-associated bins; at least two of each are needed, not estimable")
            rows.append(_row(name, "not estimable", "bins", n[GROUPS[0]], n[GROUPS[1]],
                             _median(umi[m][labels == GROUPS[0]]), _median(umi[m][labels == GROUPS[1]]),
                             note=note + "; fewer than two bins in a compartment"))
            return
        de = de_table(lognorm[m] if X is None else X, counts[m], labels, var_names, target_sum, selection,
                      a.obs_names[m])
        de_all[name] = de
        rows.append(_row(name, "estimable", "bins", n[GROUPS[0]], n[GROUPS[1]], _median(umi[m][labels == GROUPS[0]]),
                         _median(umi[m][labels == GROUPS[1]]), de, note))

    lab_m = np.where(in_sections, _labels(gp_m, comp == GROUPS[1]), "")
    bins_variant("manuscript", lab_m, f"growth plate = {rule}; tumor-associated = the other {S}-positive bins of "
                                      f"{tumor_txt} sections; {where}; Wilcoxon rank-sum over bins on the "
                                      f"log-normalized expression (Figure4h_osteoclast_compartment_DE.csv)")
    if roi is None:
        rows.append(_row("roi_only", "not estimable", "bins", note="no in_marrow_roi column"))
    else:
        bins_variant("roi_only", np.where(roi, lab_m, ""),
                     f"the manuscript bins inside the marrow region of interest (in_marrow_roi); "
                     f"{int(((lab_m != '') & ~roi).sum()):,} of {int((lab_m != '').sum()):,} tested bins lie "
                     f"outside it; {where}")

    # depth_adjusted: the deeper group downsampled to the other group's UMI distribution
    m = lab_m != ""
    labels = lab_m[m].astype(str)
    n = {g: int((labels == g).sum()) for g in GROUPS}
    if min(n.values()) < 2 or n_reps < 1:
        rows.append(_row("depth_adjusted", "not estimable", "bins", n[GROUPS[0]], n[GROUPS[1]],
                         note="fewer than two bins in a compartment" if n_reps >= 1 else "n_depth_reps is 0"))
    else:
        X = counts[m]
        targets, deeper = depth_adjusted_targets(umi[m], labels)
        other = GROUPS[1] if deeper == GROUPS[0] else GROUPS[0]
        reps = []
        for rep in range(n_reps):
            Xd = ba.downsample_counts_matrix(X, targets, seed + rep)
            reps.append(de_table(log_normalize(Xd, target_sum), Xd, labels, var_names, target_sum, selection,
                                 a.obs_names[m]))
            del Xd
        de = de_all["depth_adjusted"] = aggregate_replicates(reps)
        del reps, X
        rows.append(_row("depth_adjusted", "estimable", "bins", n[GROUPS[0]], n[GROUPS[1]],
                         _median(targets[labels == GROUPS[0]]), _median(targets[labels == GROUPS[1]]), de,
                         f"the {deeper.replace('_', '-')} bins (median {_median(umi[m][labels == deeper]):g} UMI) "
                         f"downsampled bin by bin to the {other.replace('_', '-')} bins' UMI distribution "
                         f"(rank-quantile targets), {n_reps} replicates (seed {seed} + replicate), re-log-normalized "
                         f"(target_sum {target_sum:g}); per gene the median score and adjusted P over replicates and "
                         f"the fraction of replicates significant in the same direction; {where}",
                         _median(umi[m][labels == GROUPS[0]]), _median(umi[m][labels == GROUPS[1]])))

    # block_pseudobulk: blocks x compartment per section as units, unit depths balanced
    bp_units = f"blocks of {block_um:g} um"
    if min(n.values()) < 2 or n_reps < 1:
        rows.append(_row("block_pseudobulk", "not estimable", bp_units,
                         note="fewer than two bins in a compartment" if n_reps >= 1 else
                         "n_depth_reps is 0, so the unit depths cannot be balanced"))
    else:
        idx = np.flatnonzero(m)
        de, info = block_pseudobulk(counts[idx], smp[idx], common.array_xy_um(a, bin_um)[idx], labels, block_um,
                                    var_names, target_sum, selection, n_reps, seed)
        un, bpu, uu, ut = info["n"], info["bins_per_unit"], info["umi"], info["target"]
        g, t = GROUPS
        note = (f"units = {block_um:g} um array-grid blocks x compartment per section, counts summed over the "
                f"manuscript bins of each unit ({sum(un.values())} units from {len(idx):,} bins; median bins per unit "
                f"{bpu[g]:g} growth plate vs {bpu[t]:g} tumor-associated, median UMI per unit {uu[g]:g} vs "
                f"{uu[t]:g}); {where}")
        if de is None:
            rows.append(_row("block_pseudobulk", "not estimable", bp_units, un[g], un[t], uu[g], uu[t],
                             note=note + "; fewer than two units in a compartment"))
            common.warn(progress, f"Figure 4h sensitivity (block_pseudobulk): {un[g]} growth-plate and {un[t]} "
                                  f"tumor-associated units; not estimable")
        else:
            deep = info["deeper"]
            de_all["block_pseudobulk"] = de
            rows.append(_row("block_pseudobulk", "estimable", bp_units, un[g], un[t], ut[g], ut[t], de,
                             note + f"; the {deep.replace('_', '-')} units downsampled unit by unit to the "
                             f"{(t if deep == g else g).replace('_', '-')} units' UMI distribution (rank-quantile "
                             f"targets; median UMI per unit after: {ut[g]:g} vs {ut[t]:g}), {n_reps} replicates (seed "
                             f"{seed} + replicate), CP10k, log1p; Wilcoxon rank-sum between units; per gene the "
                             f"median score and adjusted P over replicates and the fraction of replicates significant "
                             f"in the same direction", uu[g], uu[t]))

    # anatomical and drawn growth-plate definitions, both compartments from tumor sections
    bins_variant("anatomical", np.where(tumor, _labels(oc & near, oc & ~near), ""),
                 f"growth plate = {S}-positive bins within {dist_um:g} um (array grid) of a {P}-dense region "
                 f"(8-connected {calls or 'no'} bins, >= {min_bins} bins); tumor-associated = the other {S}-positive "
                 f"bins; {within_tumor(oc & near)}")
    masks = st.get("growth_plate_masks") or {}
    if not masks:
        rows.append(_row("drawn", "not configured", "bins",
                         note="no hand-drawn growth plate (sensitivity_analyses.compartments.growth_plate_masks)"))
    else:
        try:
            drawn, with_polygons = drawn_growth_plate(a, masks, cfg["roi"].get("mask_frame", "hires"),
                                                      cfg.get("data_root", "."))
        except (OSError, ValueError) as e:
            common.warn(progress, f"Figure 4h sensitivity (drawn): {e}")
            rows.append(_row("drawn", "not estimable", "bins", note=str(e)))
        else:
            bins_variant("drawn", np.where(tumor, _labels(oc & drawn, oc & ~drawn), ""),
                         f"growth plate = {S}-positive bins inside the hand-drawn growth-plate polygons "
                         f"({', '.join(with_polygons) or 'no section'}; {cfg['roi'].get('mask_frame', 'hires')} "
                         f"frame); tumor-associated = the other {S}-positive bins; {within_tumor(oc & drawn)}")

    variants = pd.DataFrame(rows, columns=VARIANT_COLUMNS)
    variants["variant"] = pd.Categorical(variants["variant"], categories=list(VARIANTS), ordered=True)
    variants = variants.sort_values("variant", kind="stable", ignore_index=True)
    variants["variant"] = variants["variant"].astype(str)
    de_out = de_subset(de_all, de_all.get("manuscript"))
    top = top_genes(de_all, selection)

    unit_txt = "bins (block_pseudobulk: blocks x compartment)"
    return [
        common.Output(DIAG_CSV, "Sensitivity — Figure 4h rule diagnostics",
                      f"The growth-plate rule of Figure 4h per section and over all bins ('all'): threshold = the "
                      f"{cp['percentile']:g}th percentile of the {P} score over all bins and whether it is exactly 0; "
                      f"% of bins above it; bins with >= 1 UMI of a {P} gene ({', '.join(selection) or 'none'}); the "
                      f"compartment sizes, how many lie outside the marrow region of interest, their median UMI per "
                      f"bin; growth-plate bins with a single {P} UMI and without any {P} core gene.", diag),
        common.Output(REGIONS_CSV, "Sensitivity — Figure 4h anatomical growth plate",
                      f"{P}-dense regions per section: 8-connected components on the array grid of {calls or 'no'} "
                      f"bins with at least {min_bins} bins (region_sizes, largest first), and the {S}-positive bins "
                      f"within {dist_um:g} um (array um) of a region bin, of which manuscript growth-plate bins.",
                      regions),
        common.Output(VARIANTS_CSV, "Sensitivity — Figure 4h variants",
                      f"The Figure 4h test under each variant: status, units ({unit_txt}), units per compartment "
                      f"(n_group = growth plate, n_rest = tumor-associated), their median UMI per unit as tested "
                      f"(median_umi_group / median_umi_rest: after the depth adjustment of depth_adjusted and "
                      f"block_pseudobulk) and before any adjustment (*_unadjusted), and the genes significantly higher "
                      f"in each compartment (BH-adjusted P < {ALPHA:g}, positive score). The anatomical and drawn "
                      f"variants take both compartments from {tumor_txt} sections. note says what the variant tests.",
                      variants),
        common.Output(DE_CSV, "Sensitivity — Figure 4h DE per variant",
                      f"Wilcoxon rank-sum (scanpy, BH) growth plate vs tumor-associated per variant, both directions; "
                      f"log2fc_mean_cp10k = log2((mean CP10k of the group + 1) / (mean CP10k of the other + 1)); "
                      f"pct_group / pct_rest = % of units detecting the gene; selection_gene = gene of the {P} "
                      f"extended list (selects the growth-plate bins); depth_adjusted and block_pseudobulk: medians "
                      f"over {n_reps} depth-adjustment replicates and frac_reps_significant (fraction of replicates "
                      f"significant in the same direction). Genes in the top {TOP_N_DE} of any variant and group or "
                      f"significant in the manuscript variant.", de_out),
        common.Output(TOP_CSV, "Sensitivity — Figure 4h genes across variants",
                      f"The genes drawn in Figure 4h (top {FIGURE_TOP_N} per compartment by manuscript score) with "
                      f"each variant's adjusted P, significance in the same direction (adjusted P < {ALPHA:g} and a "
                      f"positive score in the row's compartment; empty = variant not estimable) and score.", top),
    ]
