"""
"Within 50 um" proximity tables against other baselines (second review G1.1, F11).

The manuscript tables (bone_analysis.proximity_table) compare the percentage of source bins whose nearest target bin
lies within the radius with the same percentage over all ROI bins. Two things make that baseline unlike the source
bins: it contains the target bins themselves (0 um from a target unless proximity.exclude_self), and positivity calls
need >= 2 core genes, so the source bins are deep bins, and the chance of a target within the radius rises steeply
with depth. For every configured pair and section this module reports the observed percentage next to
  - all ROI bins (the manuscript baseline) and the largest ratio it allows (100 / baseline),
  - non-target ROI bins and bins that are neither source nor target,
  - a depth-matched expectation: each source bin takes the within-radius rate of the reference bins (ROI bins neither
    source nor target) of its total-UMI stratum; the expectation is the mean over source bins, and diff_pp =
    observed - expected. Strata are cut at the n_strata quantiles of total UMIs over the section's ROI bins AND at
    the same quantiles over its source bins. ROI quantiles alone put most source bins (median ~4x the ROI median) into
    the top ROI stratum, whose UMI range is wide and inside which the within-radius rate keeps rising, so the source
    bins get the rate of shallower reference bins: in a depth-only simulation with continuous, spatially smooth
    log-normal depth that leaves +5 to +10 pp of null difference, against about +/-1 pp with the added source cuts.
    The ROI-quantile match is still reported (_roi_strata columns, the definition of the review's G1.1 numbers), and
    median_umi_matched_reference shows how well each match balances depth,
  - the same for source bins that are not themselves targets (co-positive bins, at 0 um unless exclude_self),
with 95% spatial block-bootstrap intervals of the differences (whole square blocks resampled; stratum labels fixed,
stratum rates recomputed in every resample). Within-radius flags use exactly the positions, region of interest,
distances and exclude_self setting of the manuscript table, so pct_source_within and pct_all_within reproduce it.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from bone_sensitivity import common

FILE = "sensitivity_proximity.csv"
MIN_SOURCE = 5                      # fewer source bins (or no target bin): 'not estimable'
COLUMNS = ["pair", "source", "target", "mode", "assignment_column", "sample", "status", "n_roi", "n_source", "n_target",
           "n_copositive", "radius_um", "pct_source_within", "pct_all_within", "pct_nontarget_within",
           "pct_neither_within", "pct_source_within_excl_copositive", "expected_depth_matched_pct", "diff_pp",
           "diff_ci_lo", "diff_ci_hi", "ratio_depth_matched", "expected_excl_copositive_pct", "diff_excl_copositive_pp",
           "diff_excl_copositive_ci_lo", "diff_excl_copositive_ci_hi", "max_possible_ratio_vs_all", "median_umi_source",
           "median_umi_roi", "exclude_self", "n_depth_strata", "median_umi_matched_reference",
           "expected_roi_strata_pct", "diff_roi_strata_pp", "diff_roi_strata_ci_lo", "diff_roi_strata_ci_hi",
           "median_umi_matched_reference_roi_strata"]


def depth_strata(depth, n_strata, source=None):
    """Depth stratum (0 .. k-1) of every bin: cut at the quantiles 0, 1/n, ..., 1 of `depth` and, with `source` flags,
    also at the same quantiles of the source bins' depths (a refinement: every ROI-quantile stratum is a union of
    strata), duplicate edges merged (right-closed intervals, as pd.qcut), so ties never split."""
    d = np.asarray(depth, float)
    if d.size == 0:
        return np.zeros(0, dtype=np.int64)
    q = np.linspace(0, 1, max(int(n_strata), 1) + 1)
    edges = np.quantile(d, q)
    if source is not None and np.any(source):
        edges = np.r_[edges, np.quantile(d[np.asarray(source, bool)], q)]
    edges = np.unique(edges)
    return np.searchsorted(edges[1:-1], d, side="left").astype(np.int64)


def _borrow(n):
    """Per stratum, the nearest stratum with n > 0 (the shallower one on a tie); None when there is none."""
    have = np.flatnonzero(np.asarray(n) > 0)
    if not have.size:
        return None
    return have[np.argmin(np.abs(np.arange(len(n))[:, None] - have[None, :]), axis=1)]


def stratum_rates(hits, n):
    """Within-radius rate of the reference bins per stratum (hits / n); a stratum without reference bins borrows the
    rate of the nearest stratum that has some (the shallower one on a tie). All NaN without any reference bin."""
    hits, n = np.asarray(hits, float), np.asarray(n, float)
    nearest = _borrow(n)
    if nearest is None:
        return np.full(len(n), np.nan)
    return hits[nearest] / n[nearest]


def matched_reference_median(depth, strata, source, target):
    """Median total UMIs of the reference bins (neither source nor target) weighted as the depth match weights them:
    each reference bin carries source bins / reference bins of its stratum (a stratum without reference bins hands
    its source bins to the stratum it borrows from). It equals the source bins' median under a perfect match and is
    lower when the source bins sit at the deep end of their strata. NaN without source or reference bins."""
    d, strata = np.asarray(depth, float), np.asarray(strata, np.int64)
    src, tgt = np.asarray(source, bool), np.asarray(target, bool)
    ref = ~src & ~tgt
    k = int(strata.max()) + 1 if strata.size else 1
    n_ref, n_src = np.bincount(strata[ref], minlength=k), np.bincount(strata[src], minlength=k)
    nearest = _borrow(n_ref)
    if nearest is None or not n_src.sum():
        return np.nan
    per_bin = np.bincount(nearest, weights=n_src, minlength=k) / np.maximum(n_ref, 1)
    weight = per_bin[strata[ref]]
    order = np.argsort(d[ref], kind="stable")
    order = order[weight[order] > 0]
    v, w = d[ref][order], np.cumsum(weight[order])
    i = int(np.searchsorted(w, w[-1] / 2 * (1 - 1e-12)))
    if np.isclose(w[i], w[-1] / 2) and i + 1 < len(v):           # exactly half the weight below: average, as np.median
        return float((v[i] + v[i + 1]) / 2)
    return float(v[i])


def _codes(strata, source, target):
    """Bin code stratum + k * category (0 target only, 1 reference = neither, 2 co-positive source, 3 other source)."""
    k = int(strata.max()) + 1 if strata.size else 1
    cat = np.where(source, np.where(target, 2, 3), np.where(target, 0, 1))
    return strata + k * cat, k


def depth_matched(within, strata, source, target, idx=None):
    """[observed %, depth-matched expected %, observed % excl. co-positive, expected % excl. co-positive] of the
    source bins among the bins idx (all when None, repeats allowed). Expected = mean over the source bins of their
    stratum's reference rate (stratum_rates); NaN where there is no source bin or no reference bin."""
    within, strata = np.asarray(within, bool), np.asarray(strata, np.int64)
    code, k = _codes(strata, np.asarray(source, bool), np.asarray(target, bool))
    return _depth_matched(code, k, within.astype(float), idx)


def _depth_matched(code, k, within_f, idx=None):
    c = code if idx is None else code[idx]
    w = within_f if idx is None else within_f[idx]
    n = np.bincount(c, minlength=4 * k).reshape(4, k)
    h = np.bincount(c, weights=w, minlength=4 * k).reshape(4, k)
    rate = stratum_rates(h[1], n[1])
    out = []
    for n_s, h_s in ((n[2] + n[3], h[2] + h[3]), (n[3], h[3])):
        tot = n_s.sum()
        if tot == 0:
            out += [np.nan, np.nan]
            continue
        out += [100 * h_s.sum() / tot, 100 * float((n_s * rate).sum()) / tot]
    return np.asarray(out, float)


def _pct(mask_values):
    return 100 * float(mask_values.mean()) if mask_values.size else np.nan


def section_row(P, source, target, depth, radius_um=50.0, n_strata=10, block_um=200.0, n_boot=500, seed=0,
                exclude_self=False):
    """Every metric column of one pair in one section from its ROI bins: positions P (um), source and target flags,
    total UMIs per bin. Within-radius flags: bone_analysis._nearest_target_distance(P, target, exclude_self) <=
    radius_um, as in the manuscript table. Depth match: strata at the ROI and source quantiles (depth_strata with
    source); _roi_strata columns: at the ROI quantiles only."""
    import bone_analysis as ba
    P = np.asarray(P, float)
    src, tgt = np.asarray(source, bool), np.asarray(target, bool)
    depth = np.asarray(depth, float)
    n_src, n_tgt = int(src.sum()), int(tgt.sum())
    within = ba._nearest_target_distance(P, tgt, exclude_self) <= radius_um
    neither, excl = ~src & ~tgt, src & ~tgt
    row = {"status": "estimable" if n_src >= MIN_SOURCE and n_tgt > 0 else "not estimable", "n_roi": int(len(P)),
           "n_source": n_src, "n_target": n_tgt, "n_copositive": int((src & tgt).sum()), "radius_um": float(radius_um),
           "pct_source_within": _pct(within[src]), "pct_all_within": _pct(within),
           "pct_nontarget_within": _pct(within[~tgt]), "pct_neither_within": _pct(within[neither]),
           "pct_source_within_excl_copositive": _pct(within[excl]),
           "median_umi_source": float(np.median(depth[src])) if n_src else np.nan,
           "median_umi_roi": float(np.median(depth)) if len(depth) else np.nan, "exclude_self": bool(exclude_self)}
    row["max_possible_ratio_vs_all"] = float(common.fold(100.0, row["pct_all_within"]))
    nan = dict.fromkeys(["n_depth_strata", "expected_depth_matched_pct", "diff_pp", "diff_ci_lo", "diff_ci_hi",
                         "ratio_depth_matched", "expected_excl_copositive_pct", "diff_excl_copositive_pp",
                         "diff_excl_copositive_ci_lo", "diff_excl_copositive_ci_hi", "median_umi_matched_reference",
                         "expected_roi_strata_pct", "diff_roi_strata_pp", "diff_roi_strata_ci_lo",
                         "diff_roi_strata_ci_hi", "median_umi_matched_reference_roi_strata"], np.nan)
    if row["status"] != "estimable":
        return {**row, **nan}
    strata, roi_strata = depth_strata(depth, n_strata, src), depth_strata(depth, n_strata)
    code, k = _codes(strata, src, tgt)
    code_r, k_r = _codes(roi_strata, src, tgt)
    within_f = within.astype(float)
    obs, exp, obs_x, exp_x = _depth_matched(code, k, within_f)
    exp_r = _depth_matched(code_r, k_r, within_f)[1]
    row.update(nan)
    row.update({"n_depth_strata": int(len(np.unique(strata))), "expected_depth_matched_pct": exp,
                "diff_pp": obs - exp, "ratio_depth_matched": float(common.fold(obs, exp)),
                "expected_excl_copositive_pct": exp_x, "diff_excl_copositive_pp": obs_x - exp_x,
                "median_umi_matched_reference": matched_reference_median(depth, strata, src, tgt),
                "expected_roi_strata_pct": exp_r, "diff_roi_strata_pp": obs - exp_r,
                "median_umi_matched_reference_roi_strata": matched_reference_median(depth, roi_strata, src, tgt)})
    if int(n_boot) > 0:
        def stat(idx):
            o, e, ox, ex = _depth_matched(code, k, within_f, idx)
            return np.array([o - e, ox - ex, o - _depth_matched(code_r, k_r, within_f, idx)[1]])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)       # all-NaN resamples (no co-positive-free source)
            lo, hi, _ = common.block_bootstrap(stat, common.spatial_blocks(P, block_um), n_boot, seed)
        row.update({"diff_ci_lo": float(lo[0]), "diff_ci_hi": float(hi[0]),
                    "diff_excl_copositive_ci_lo": float(lo[1]), "diff_excl_copositive_ci_hi": float(hi[1]),
                    "diff_roi_strata_ci_lo": float(lo[2]), "diff_roi_strata_ci_hi": float(hi[2])})
    return row


def _call_sets(a, src, dst, mode):
    """[(assignment_column, source flags, target flags)] over all bins: 'program' (and 'program_dm' when present)
    for mode assignment, the pos_<signature> columns ('pos') for mode positivity."""
    if mode == "assignment":
        labels = {col: a.obs[col].astype(str).to_numpy() for col in ("program", "program_dm") if col in a.obs}
        return [(col, lab == src, lab == dst) for col, lab in labels.items()]
    if f"pos_{src}" in a.obs and f"pos_{dst}" in a.obs:
        return [("pos", common.as_bool(a.obs[f"pos_{src}"]), common.as_bool(a.obs[f"pos_{dst}"]))]
    return []


def run(a, cfg, progress=None):
    """sensitivity_proximity.csv: every proximity pair x call set x section against the all-ROI, non-target, neither
    and depth-matched baselines, with spatial block-bootstrap intervals of the depth-matched differences."""
    st, top = common.settings(cfg, "proximity"), common.settings(cfg, None)
    prox = cfg["proximity"]
    radius, exclude_self = float(prox["radius_um"]), bool(prox.get("exclude_self", False))
    seed = int(cfg.get("seed", 0))
    smp = a.obs["sample"].astype(str).to_numpy()
    roi = common.as_bool(a.obs["in_marrow_roi"])
    xy = common.pixel_xy_um(a)
    if "total_counts" in a.obs:
        depth = a.obs["total_counts"].to_numpy(float)
    else:
        depth = np.asarray(a.layers["counts"].sum(axis=1), dtype=float).ravel()
    rows = []
    for pair in prox.get("pairs", []):
        src, dst, mode = pair["from"], pair["to"], pair.get("mode", "assignment")
        sets = _call_sets(a, src, dst, mode) if src in cfg["signatures"] and dst in cfg["signatures"] else []
        if not sets:
            common.warn(progress, f"proximity sensitivity {src} -> {dst}: the {mode} calls are not available; skipped")
            continue
        for col, fr, to in sets:
            for s in common.sample_ids(a, cfg):
                m = (smp == s) & roi
                r = section_row(xy[m], fr[m], to[m], depth[m], radius, st["n_strata"], st["block_um"], top["n_boot"],
                                seed, exclude_self)
                rows.append({"pair": f"{src}_to_{dst}", "source": src, "target": dst, "mode": mode,
                             "assignment_column": col, "sample": s, **r})
                if r["status"] != "estimable":
                    common.warn(progress, f"proximity sensitivity {src} -> {dst} ({col}), {s}: {r['n_source']} source "
                                          f"and {r['n_target']} target bins; depth-matched comparison not estimable")
    df = pd.DataFrame(rows).reindex(columns=COLUMNS)
    df["n_depth_strata"] = pd.to_numeric(df["n_depth_strata"]).astype("Int64")
    target_note = ("target bins are measured to the nearest other target bin" if exclude_self else
                   "target bins sit at 0 um, so the all-ROI baseline counts them as within")
    desc = (f"Percentage of source bins (marrow ROI) whose nearest target bin lies within {radius:g} um "
            f"(exclude_self={exclude_self}: {target_note}), against all ROI bins (manuscript baseline), non-target ROI "
            f"bins, bins that are neither source nor target, and the depth-matched expectation (reference = ROI bins "
            f"neither source nor target; total-UMI strata cut at the {st['n_strata']} quantiles over the section's ROI "
            f"bins and at the same quantiles over its source bins); diff_pp = observed - expected, with 95% spatial "
            f"block-bootstrap intervals ({st['block_um']:g} um blocks, {top['n_boot']} resamples). _roi_strata: the "
            f"same with strata at the ROI quantiles only, which leaves source bins at the deep end of the top ROI "
            f"stratum and biases the difference upward when depth drives proximity. median_umi_matched_reference: "
            f"median UMI of the reference bins as weighted by the match (equals median_umi_source when depth is "
            f"balanced). _excl_copositive: source bins that are not themselves targets. max_possible_ratio_vs_all = "
            f"100 / pct_all_within. assignment_column: program / program_dm (assignment pairs) or pos (positivity "
            f"pairs).")
    return [common.Output(FILE, "Sensitivity — proximity baselines", desc, df)]
