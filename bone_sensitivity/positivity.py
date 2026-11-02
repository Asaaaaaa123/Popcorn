"""
Positivity under other regions of interest, calling rules and depth-matching schemes (Figures 4g, 6a; second review
F2-F6, F12, F23, F24).

The manuscript calls a bin positive for a signature when its score exceeds the null threshold (95th percentile of
random gene sets), at least scoring.min_core core genes are detected and no competing signature scores higher, and
counts the calls inside the transcriptional region of interest (in_marrow_roi). The review found that the threshold
almost never decides a call (F3), that the region of interest keeps much skeletal muscle in the control section (F4)
and excludes a bin on one UMI of a smooth-muscle or nerve gene (F5), that the competitor veto acts almost only in the
control section (F23), that depth-matched magnitudes depend on the matching scheme (F6, F24), that positive bins
cluster so binomial intervals are too narrow (F12), and that the A549 lesion lies inside the tumor ROI (F2). The calls
are recounted under
  regions of interest  manuscript; call_rule_exclusion (muscle and nerve bins removed by the positivity rule instead
                       of the percentile rule); assignment_exclusion (bins assigned to them in Figure 3); no_lesion
  calling rules        manuscript; min_core_<k> (k core genes); no_competitor_veto
  depth                native; cap, quantile (whole sections) and quantile_roi (ROI bins only) depth matching
with Wilson and spatial block-bootstrap intervals of the percentages, and Katz and (native manuscript rule) spatial
block-bootstrap intervals of the folds. Region-of-interest masks are defined at native depth and applied unchanged to
the depth-matched calls, as in bone_analysis.depth_matched.
"""
from __future__ import annotations

import gc
import re
import warnings

import numpy as np
import pandas as pd

from bone_sensitivity import common

FOCUS = ("osteoclast", "osteoblast", "macrophage")           # signatures counted in the region-of-interest tables
MUSCLE = ("skeletal_muscle", "smooth_muscle", "nerve")        # assignment shares reported for every ROI
ROI_VARIANTS = ("manuscript", "call_rule_exclusion", "assignment_exclusion", "no_lesion")
SCHEMES = ("cap", "quantile", "quantile_roi")
POSITIVITY_CSV = "sensitivity_positivity.csv"
FOLDS_CSV = "sensitivity_positivity_folds.csv"
COMPOSITION_CSV = "sensitivity_roi_composition.csv"
EXCLUSION_CSV = "sensitivity_roi_exclusion_rules.csv"
DIAGNOSTICS_CSV = "sensitivity_positivity_rule_diagnostics.csv"
DEPTH_CSV = "sensitivity_positivity_depth_targets.csv"
POSITIVITY_COLUMNS = ["depth", "rule", "roi", "signature", "sample", "n_bins", "n_positive", "n_positive_sd",
                      "pct_positive", "ci_lo", "ci_hi", "fold_vs_reference", "block_ci_lo", "block_ci_hi",
                      "n_candidates", "n_vetoed", "manuscript"]
FOLD_COLUMNS = ["depth", "rule", "roi", "signature", "sample", "fold", "n_positive", "n_positive_reference",
                "fold_ci_lo", "fold_ci_hi", "block_fold_ci_lo", "block_fold_ci_hi", "manuscript"]
DIAGNOSTIC_COLUMNS = ["depth", "signature", "sample", "threshold", "n_ncore_ge_min", "n_removed_by_threshold",
                      "n_vetoed", "n_positive"]


def _ba():
    import bone_analysis as ba
    return ba


def _info(progress, message):
    if progress is not None:
        progress.info(message)
    else:
        print(message)


def signature_dicts(cfg):
    """Scoring dictionaries of every configured signature, built as bone_analysis._run builds them."""
    return {k: {"core": list(d["core"]), "extended": list(d["extended"]), "competes": list(d.get("competes") or []),
                "excluded": dict(d.get("excluded") or {})} for k, d in cfg["signatures"].items()}


def rule_parameters(rule, min_core):
    """(minimum core genes, competitor veto) of a calling rule: manuscript (scoring.min_core, veto), min_core_<k>
    (k core genes, veto) or no_competitor_veto (scoring.min_core, no veto)."""
    if rule == "manuscript":
        return int(min_core), True
    if rule == "no_competitor_veto":
        return int(min_core), False
    m = re.fullmatch(r"min_core_(\d+)", str(rule))
    if m:
        return int(m.group(1)), True
    raise ValueError(f"unknown calling rule {rule!r}: manuscript, min_core_<k> or no_competitor_veto")


def katz_interval(x1, n1, x2, n2, z=common.Z95):
    """Ratio of two independent proportions (x1/n1) / (x2/n2) and its Katz 95% interval, exp(log ratio +/- z SE) with
    SE^2 = 1/x1 - 1/n1 + 1/x2 - 1/n2. Bins are treated as independent. Returns (ratio, lo, hi); the interval is NaN
    where a count or a denominator is not positive, the ratio where the reference proportion is not positive."""
    x1, n1, x2, n2 = (np.asarray(v, float) for v in (x1, n1, x2, n2))
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where((n1 > 0) & (n2 > 0) & (x2 > 0), (x1 / n1) / (x2 / n2), np.nan)
        ok = (x1 > 0) & (x2 > 0) & (n1 > 0) & (n2 > 0)
        se = np.sqrt(np.where(ok, 1 / x1 - 1 / n1 + 1 / x2 - 1 / n2, np.nan))
        return ratio, np.where(ok, ratio * np.exp(-z * se), np.nan), np.where(ok, ratio * np.exp(z * se), np.nan)


def quantile_roi_targets(totals, samples, roi, shallow=None):
    """Per-bin UMI targets of the quantile_roi depth matching. The reference distribution is the UMI per bin of the
    region-of-interest bins of the shallowest section (lowest median UMI over all its bins, as in
    bone_analysis.depth_match_target, unless `shallow` names it); the ROI bins of every other section, ranked among
    themselves ((rank - 0.5) / n, ties in order), are downsampled to the UMI count at the same quantile of that
    distribution, never above their own total. Every other bin keeps its counts. The quantile scheme matches whole
    sections, which leaves the tumor ROI shallower than the control ROI (F24); this one matches the bins that are
    counted. Returns (targets, shallowest section), targets None when that section has no ROI bin."""
    totals = np.asarray(totals, dtype=float)
    smp = np.asarray(samples).astype(str)
    roi = np.asarray(roi, dtype=bool)
    if shallow is None:
        shallow = str(pd.Series(totals).groupby(smp, sort=False).median().idxmin())
    ref = np.sort(totals[(smp == shallow) & roi])
    if not ref.size:
        return None, shallow
    targets = totals.copy()
    for s in pd.unique(smp):
        m = np.flatnonzero((smp == s) & roi)
        if s == shallow or not m.size:
            continue
        q = (pd.Series(totals[m]).rank(method="first").to_numpy() - 0.5) / len(m)
        targets[m] = np.minimum(totals[m], np.floor(np.quantile(ref, q)))
    return targets.astype(np.int64), shallow


def null_thresholds(b, names, sigs, pct, n_perm, n_sub, ctrl_size, seed):
    """bone_analysis.null_threshold of each signature in `names`. The null depends on a signature only through the
    number of its scored genes present (the same random sets are drawn for the same size), so it is computed once per
    set size; the values equal those of call_positive."""
    ba = _ba()
    out, cache = {}, {}
    for n in names:
        size = len(ba.present(b, sigs[n]["extended"]))
        if size not in cache:
            cache[size] = ba.null_threshold(b, n, sigs=sigs, pct=pct, n_perm=n_perm, n_sub=n_sub,
                                            ctrl_size=ctrl_size, seed=seed)
        out[n] = cache[size]
    return out


class Calls:
    """What every calling rule needs at one depth, bins x signatures: score above the null threshold (`passes`), core
    genes detected (`ncore`) and a competing signature scoring strictly higher (`vetoed`; a competitor that cannot be
    scored vetoes nothing, as in bone_analysis.call_positive)."""

    def __init__(self, obs, names, sigs, thresholds):
        n = len(obs)

        def column(c):
            return obs[c].to_numpy(dtype=float) if c in obs else np.full(n, np.nan)

        self.names = list(names)
        self.thresholds = np.array([np.nan if thresholds.get(x) is None else float(thresholds[x]) for x in names])
        S = np.column_stack([column(f"score_{x}") for x in names])
        rival = np.full_like(S, np.nan)
        for j, x in enumerate(names):
            for c in sigs[x]["competes"]:
                if f"score_{c}" in obs:
                    rival[:, j] = np.fmax(rival[:, j], column(f"score_{c}"))
        with np.errstate(invalid="ignore"):
            self.passes = S > self.thresholds
            self.vetoed = rival > S
        self.ncore = np.nan_to_num(np.column_stack([column(f"ncore_{x}") for x in names])).astype(np.int32)

    def positive(self, min_core, veto=True):
        """Positive bins under a rule: score above threshold, at least min_core core genes, optionally the veto."""
        pos = self.passes & (self.ncore >= min_core)
        return pos & ~self.vetoed if veto else pos

    def diagnostics(self, idx, min_core):
        """Within the bins idx: bins with at least min_core core genes, those of them at or below the threshold, those
        passing both rules but vetoed, and the positive bins (one value per signature)."""
        core, passes, vetoed = self.ncore[idx] >= min_core, self.passes[idx], self.vetoed[idx]
        return {"n_ncore_ge_min": core.sum(axis=0), "n_removed_by_threshold": (core & ~passes).sum(axis=0),
                "n_vetoed": (core & passes & vetoed).sum(axis=0), "n_positive": (core & passes & ~vetoed).sum(axis=0)}


def block_intervals(pos, blocks, n_boot=500, seed=0):
    """Spatial block-bootstrap 95% interval of the percentage of positive bins (each column of pos) among the given
    bins: whole spatial blocks are resampled with replacement (common.block_bootstrap). Blocks are first reduced to
    their positive and total bins; resampling these sums is the same as resampling the blocks' bins. Returns
    (lo, hi, draws): the percentile interval and the resampled percentages (n_boot x columns; none with fewer than
    two blocks)."""
    pos = np.asarray(pos, dtype=float)
    _, inv = np.unique(np.asarray(blocks), return_inverse=True)
    size = np.bincount(inv).astype(float)
    hits = np.column_stack([np.bincount(inv, weights=pos[:, j], minlength=len(size)) for j in range(pos.shape[1])])

    def pct(ids):
        return 100 * hits[ids].sum(axis=0) / size[ids].sum()

    return common.block_bootstrap(pct, np.arange(len(size)), n_boot=n_boot, seed=seed)


def block_fold_interval(draws, draws_ref, alpha=0.05):
    """Percentile 95% interval of the fold of two sections' percentages from their spatial block-bootstrap draws
    (block_intervals, each section with its own seed so the two are resampled independently; draws x signatures,
    paired by draw). Unlike the Katz interval it allows for the clustering of positive bins (F12). A draw whose
    reference percentage is 0 has an unbounded fold (and is dropped when its own percentage is 0 too); a bound of the
    interval that falls among such draws is NaN (unbounded). Returns (lo, hi), NaN without draws."""
    d = np.asarray(draws, dtype=float)
    r = np.asarray(draws_ref, dtype=float) if draws_ref is not None else np.empty((0,) + d.shape[1:])
    n = min(len(d), len(r))
    if n == 0:
        return np.full(d.shape[1:], np.nan), np.full(d.shape[1:], np.nan)
    d, r = d[:n], r[:n]
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(r > 0, d / r, np.where(d > 0, np.inf, np.nan))
    top = np.max(np.where(np.isfinite(f), f, -np.inf), axis=0)       # largest bounded draw
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lo, hi = np.nanpercentile(np.where(np.isinf(f), np.finfo(float).max, f),
                                  [100 * alpha / 2, 100 * (1 - alpha / 2)], axis=0)
    return np.where(lo > top, np.nan, lo), np.where(hi > top, np.nan, hi)


def lesion_flags(a, cfg, progress=None):
    """In-lesion flag per bin: obs in_lesion when present, otherwise common.lesion_signal with the settings of
    sensitivity_analyses.lesion (not written to obs)."""
    if "in_lesion" in a.obs:
        return common.as_bool(a.obs["in_lesion"])
    les = common.settings(cfg, "lesion")
    _, flag, found = common.lesion_signal(a, tuple(les["genes"]), int(les["pool_k"]), int(les["min_pooled_umi"]))
    if not found:
        common.warn(progress, f"positivity sensitivity: none of the lesion genes ({', '.join(les['genes'])}) is in "
                              f"the data, so the no_lesion region equals the manuscript region")
    return np.asarray(flag, dtype=bool)


def exclusion_rules(a, cfg, sigs, progress=None):
    """Bins excluded from the region of interest by each roi.exclude signature under two rules.
    manuscript: bone_analysis.define_marrow_roi; a signature with a positivity call (pos_ column) excludes its positive
    bins, any other a bin scoring above its roi.exclude_percentile over all bins of all sections.
    call: the positivity rule of bone_analysis.call_positive at native depth (score above the null threshold with
    scoring.n_perm sets, at least scoring.min_core core genes), without competitors; a pos_ column is used as is.
    Returns {signature: {"rule", "percentile", "threshold", "manuscript", "call"}} (percentile: value of the score
    percentile that excludes bins, NaN when the signature is excluded through its positivity call or cannot be
    scored)."""
    sco, pct = cfg["scoring"], float(cfg["roi"]["exclude_percentile"])
    thr_uns = dict(a.uns.get("positivity_thresholds") or {})
    none = np.zeros(a.n_obs, dtype=bool)
    out = {}
    for e in cfg["roi"]["exclude"]:
        if f"pos_{e}" in a.obs:
            pos = common.as_bool(a.obs[f"pos_{e}"])
            out[e] = {"rule": "positivity call", "percentile": np.nan, "manuscript": pos, "call": pos,
                      "threshold": float(thr_uns.get(e, np.nan))}
            continue
        score = a.obs[f"score_{e}"].to_numpy(dtype=float) if f"score_{e}" in a.obs else np.full(a.n_obs, np.nan)
        value = float(np.nanpercentile(score, pct)) if not np.isnan(score).all() else np.nan
        if np.isnan(value) or e not in sigs:
            common.warn(progress, f"positivity sensitivity: region-of-interest signature {e} cannot be scored, so it "
                                  f"excludes no bin")
            out[e] = {"rule": "not scored", "percentile": np.nan, "manuscript": none, "call": none,
                      "threshold": np.nan}
            continue
        t = null_thresholds(a, [e], sigs, sco["null_percentile"], sco["n_perm"], sco["null_n_sub"],
                            sco["ctrl_size"], int(cfg["seed"]))[e]
        ncore = a.obs[f"ncore_{e}"].to_numpy() if f"ncore_{e}" in a.obs else np.zeros(a.n_obs)
        with np.errstate(invalid="ignore"):
            out[e] = {"rule": "percentile", "percentile": value, "manuscript": score > value, "threshold": t,
                      "call": (score > t) & (ncore >= int(sco["min_core"]))}
    return out


def roi_variants(a, cfg, variants, excl, in_lesion, progress=None):
    """Region-of-interest masks by name (the manuscript one always first). call_rule_exclusion and
    assignment_exclusion replace the transcriptional part of the manuscript region and keep the hand-drawn region
    (in_roi_drawn) when there is one; no_lesion is the manuscript region without the A549 lesion."""
    base = common.as_bool(a.obs["in_marrow_roi"])
    drawn = common.as_bool(a.obs["in_roi_drawn"]) if "in_roi_drawn" in a.obs else np.ones(a.n_obs, dtype=bool)
    out = {"manuscript": base}
    for v in variants:
        if v == "manuscript":
            continue
        if v == "call_rule_exclusion":
            keep = drawn.copy()
            for d in excl.values():
                keep &= ~d["call"]
            out[v] = keep
        elif v == "assignment_exclusion":
            if "program" not in a.obs:
                common.warn(progress, "positivity sensitivity: no tissue assignment ('program'), so the "
                                      "assignment_exclusion region is skipped")
                continue
            out[v] = drawn & ~a.obs["program"].astype(str).isin(list(cfg["roi"]["exclude"])).to_numpy()
        elif v == "no_lesion":
            out[v] = base & ~in_lesion
        else:
            common.warn(progress, f"positivity sensitivity: unknown region of interest {v!r} skipped "
                                  f"({', '.join(ROI_VARIANTS)})")
    return out


def _rules(variants, min_core, progress):
    rules = {"manuscript": rule_parameters("manuscript", min_core)}
    for r in variants:
        try:
            rules.setdefault(r, rule_parameters(r, min_core))
        except ValueError as exc:
            common.warn(progress, f"positivity sensitivity: {exc}; skipped")
    return rules


def _count(calls, rules, rois, groups, samples):
    """{(rule, roi, sample): positive bins per signature} for the given rules and regions."""
    out = {}
    for r in rules:
        pos = calls.positive(*rules[r])
        for roi in rois:
            for s in samples:
                out[(r, roi, s)] = pos[groups[(roi, s)]].sum(axis=0)
    return out


def _records(depth, counts, names, groups):
    """Long rows of the positivity table from {(rule, roi, sample): counts per signature or replicates x signatures}."""
    rows = []
    for (rule, roi, s), k in counts.items():
        k = np.asarray(k, dtype=float)
        mean = k if k.ndim == 1 else k.mean(axis=0)
        sd = np.full(len(names), np.nan) if k.ndim == 1 or len(k) < 2 else k.std(axis=0, ddof=1)
        for j, n in enumerate(names):
            rows.append({"depth": depth, "rule": rule, "roi": roi, "signature": n, "sample": s,
                         "n_bins": len(groups[(roi, s)]), "n_positive": mean[j], "n_positive_sd": sd[j]})
    return rows


def _section_blocks(a, smp, block_um, progress):
    """Spatial block of every bin: squares of side block_um anchored at each section's lowest coordinates, so every
    region of interest of a section shares the same grid. None when the coordinates cannot be converted to um."""
    try:
        xy = common.pixel_xy_um(a)
    except (KeyError, TypeError, ValueError) as exc:
        common.warn(progress, f"positivity sensitivity: no microns per pixel ({exc}); no block-bootstrap intervals")
        return None
    blocks = np.zeros(a.n_obs, dtype=np.int64)
    for s in pd.unique(smp):
        m = smp == s
        blocks[m] = common.spatial_blocks(xy[m], block_um)
    return blocks


def _positivity_table(rows, extras, reference, configured):
    df = pd.DataFrame(rows, columns=["depth", "rule", "roi", "signature", "sample", "n_bins", "n_positive",
                                     "n_positive_sd"])
    n = df["n_bins"].to_numpy(dtype=float)
    k = df["n_positive"].to_numpy(dtype=float)
    native = (df["depth"] == "native").to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.where(n > 0, 100 * k / n, np.nan)
    lo, hi = common.wilson(np.where(native, k, np.round(k)), n)
    df["pct_positive"], df["ci_lo"], df["ci_hi"] = pct, lo, hi
    key = ["depth", "rule", "roi", "signature"]
    ref = df[df["sample"] == reference][key + ["pct_positive"]].rename(columns={"pct_positive": "_ref"})
    df = df.merge(ref, on=key, how="left")
    df["fold_vs_reference"] = common.fold(df["pct_positive"], df["_ref"])
    ex = pd.DataFrame(extras, columns=["roi", "signature", "sample", "block_ci_lo", "block_ci_hi", "n_candidates",
                                       "n_vetoed"]).assign(depth="native", rule="manuscript")
    df = df.drop(columns="_ref").merge(ex, on=["depth", "rule", "roi", "signature", "sample"], how="left")
    df["manuscript"] = (df["rule"] == "manuscript") & (df["roi"] == "manuscript") & \
        df["depth"].isin(["native", configured])
    df[["n_candidates", "n_vetoed"]] = df[["n_candidates", "n_vetoed"]].astype("Int64")
    return df[POSITIVITY_COLUMNS].round({"pct_positive": 5, "ci_lo": 5, "ci_hi": 5, "block_ci_lo": 5,
                                         "block_ci_hi": 5})


def _folds_table(pos, reference, block_folds):
    """Folds of every non-reference section over the reference: Katz intervals on every row, the block-bootstrap
    intervals (block_folds rows) on the native manuscript-rule rows; both NaN when a count is 0."""
    key = ["depth", "rule", "roi", "signature"]
    ref = pos[pos["sample"] == reference][key + ["n_positive", "n_bins"]] \
        .rename(columns={"n_positive": "n_positive_reference", "n_bins": "_n_ref"})
    df = pos[pos["sample"] != reference].merge(ref, on=key, how="left")
    df["fold"], df["fold_ci_lo"], df["fold_ci_hi"] = katz_interval(df["n_positive"], df["n_bins"],
                                                                   df["n_positive_reference"], df["_n_ref"])
    bf = pd.DataFrame(block_folds, columns=["roi", "signature", "sample", "block_fold_ci_lo", "block_fold_ci_hi"]) \
        .assign(depth="native", rule="manuscript")
    df = df.merge(bf, on=key + ["sample"], how="left")
    zero = ~((df["n_positive"] > 0) & (df["n_positive_reference"] > 0))
    df.loc[zero, ["block_fold_ci_lo", "block_fold_ci_hi"]] = np.nan
    return df[FOLD_COLUMNS].reset_index(drop=True)


def _diagnostic_rows(depth, diag, thresholds, names, samples):
    """diag: {sample: {column: values per signature}} or a list of them (replicates, averaged)."""
    reps = diag if isinstance(diag, list) else [diag]
    thr = np.mean(np.vstack(thresholds), axis=0)
    rows = []
    for s in samples:
        vals = {c: np.mean([r[s][c] for r in reps], axis=0) for c in DIAGNOSTIC_COLUMNS[4:]}
        for j, n in enumerate(names):
            rows.append({"depth": depth, "signature": n, "sample": s, "threshold": thr[j],
                         **{c: float(v[j]) for c, v in vals.items()}})
    return rows


def run(a, cfg, progress=None):
    """Tables of the positivity sensitivity analysis (module docstring). Adds nothing to `a`."""
    ba = _ba()
    st, top, sco = common.settings(cfg, "positivity"), common.settings(cfg, None), cfg["scoring"]
    sigs, seed, min_core = signature_dicts(cfg), int(cfg["seed"]), int(sco["min_core"])
    target_sum, configured = float(cfg["clustering"]["target_sum"]), str(sco["depth_matching"])
    if "in_marrow_roi" not in a.obs:
        common.warn(progress, "positivity sensitivity: no region of interest (in_marrow_roi); nothing computed")
        return []
    names = [n for n in cfg["primary"]
             if n in sigs and all(f"{p}_{n}" in a.obs for p in ("score", "ncore", "pos"))]
    left_out = [n for n in cfg["primary"] if n not in names]
    if left_out:
        common.warn(progress, f"positivity sensitivity: no scores or calls for {', '.join(left_out)}; left out")
    if not names:
        return []
    samples = common.sample_ids(a, cfg)
    reference = cfg["reference_sample"]
    if reference not in samples:
        common.warn(progress, f"positivity sensitivity: reference section {reference} is not in the data; folds "
                              f"are over {samples[0]}")
        reference = samples[0]
    smp = a.obs["sample"].astype(str).to_numpy()
    in_lesion = lesion_flags(a, cfg, progress)
    excl = exclusion_rules(a, cfg, sigs, progress)
    rois = roi_variants(a, cfg, st["roi_variants"], excl, in_lesion, progress)
    rules = _rules(st["rule_variants"], min_core, progress)
    groups = {(r, s): np.flatnonzero(m & (smp == s)) for r, m in rois.items() for s in samples}

    # native depth: the scores and thresholds of the manuscript run
    thr = dict(a.uns.get("positivity_thresholds") or {})
    missing = [n for n in names if thr.get(n) is None]
    if missing:
        common.warn(progress, f"positivity sensitivity: no stored threshold for {', '.join(missing)}; recomputed")
        thr.update(null_thresholds(a, missing, sigs, sco["null_percentile"], sco["n_perm"], sco["null_n_sub"],
                                   sco["ctrl_size"], seed))
    native = Calls(a.obs, names, sigs, thr)
    manuscript_pos = native.positive(*rules["manuscript"])
    for j, n in enumerate(names):
        diff = int((manuscript_pos[:, j] != common.as_bool(a.obs[f"pos_{n}"])).sum())
        if diff:
            common.warn(progress, f"positivity sensitivity: the manuscript rule differs from pos_{n} in {diff} bins")
    rows = _records("native", _count(native, rules, rois, groups, samples), names, groups)
    blocks = _section_blocks(a, smp, float(top["block_um"]), progress)
    candidates = native.passes & (native.ncore >= min_core)
    extras, draws = [], {}
    for (roi, s), idx in groups.items():
        lo = hi = np.full(len(names), np.nan)
        if blocks is not None and len(idx):             # one seed per section: sections are resampled independently
            lo, hi, draws[(roi, s)] = block_intervals(manuscript_pos[idx], blocks[idx], int(top["n_boot"]),
                                                      seed + samples.index(s))
        cand = candidates[idx]
        vet = (cand & native.vetoed[idx]).sum(axis=0)
        extras += [{"roi": roi, "signature": n, "sample": s, "block_ci_lo": lo[j], "block_ci_hi": hi[j],
                    "n_candidates": int(cand[:, j].sum()), "n_vetoed": int(vet[j])} for j, n in enumerate(names)]
    block_folds = []
    for (roi, s), d in draws.items():
        if s != reference:
            lo, hi = block_fold_interval(d, draws.get((roi, reference)))
            block_folds += [{"roi": roi, "signature": n, "sample": s, "block_fold_ci_lo": lo[j],
                             "block_fold_ci_hi": hi[j]} for j, n in enumerate(names)]
    roi_idx = {s: groups[("manuscript", s)] for s in samples}
    diag_rows = _diagnostic_rows("native", {s: native.diagnostics(roi_idx[s], min_core) for s in samples},
                                 [native.thresholds], names, samples)

    # depth matching: every replicate rescored and its thresholds recomputed, as bone_analysis.depth_matched does
    totals = a.obs["total_counts"].to_numpy(dtype=float)
    roi_base = rois["manuscript"]
    depth_rows = [_depth_row("native", s, totals, totals, roi_base, smp) for s in samples]
    needed = {}
    for n in names:
        needed[n] = sigs[n]
        needed.update({c: sigs[c] for c in sigs[n]["competes"] if c in sigs})
    n_reps = int(sco["n_depth_reps"])
    ran, full_ran, diag_ran = [], [], False
    for scheme in dict.fromkeys(st["depth_schemes"]):
        if scheme not in SCHEMES:
            common.warn(progress, f"positivity sensitivity: unknown depth matching {scheme!r} skipped "
                                  f"({', '.join(SCHEMES)})")
            continue
        if n_reps < 1:
            common.warn(progress, "positivity sensitivity: scoring.n_depth_reps is 0; no depth-matched counts")
            break
        if scheme == "quantile_roi":
            shallow = str(ba.depth_match_target(a)[1].idxmin())
            targets, _ = quantile_roi_targets(totals, smp, roi_base, shallow)
            if targets is None:
                common.warn(progress, f"positivity sensitivity: {shallow} has no ROI bin; quantile_roi skipped")
                continue
        else:
            targets = ba.depth_match_targets(a, scheme)[0]
        depth_rows += [_depth_row(scheme, s, totals, np.minimum(totals, targets), roi_base, smp) for s in samples]
        full = scheme in ("cap", configured)
        run_rules = rules if full else {"manuscript": rules["manuscript"]}
        run_rois = list(rois) if full else ["manuscript"]
        reps, diag, thr_reps = {}, [], []
        for rep in range(n_reps):
            _info(progress, f"positivity sensitivity: {scheme} depth matching, replicate {rep + 1} of {n_reps}")
            X = ba.downsample_counts_matrix(a.layers["counts"], targets, seed + rep)
            b = ba._lognorm_object(X, a, ["sample"], target_sum)
            del X
            ba.score_all(b, needed, sco["ctrl_size"], seed)
            t = null_thresholds(b, names, sigs, sco["null_percentile"], sco["n_perm_resampled"], sco["null_n_sub"],
                                sco["ctrl_size"], seed)
            calls = Calls(b.obs, names, sigs, t)
            del b
            gc.collect()
            for key, k in _count(calls, run_rules, run_rois, groups, samples).items():
                reps.setdefault(key, []).append(k)
            if scheme == configured:
                diag.append({s: calls.diagnostics(roi_idx[s], min_core) for s in samples})
                thr_reps.append(calls.thresholds)
        rows += _records(scheme, {k: np.vstack(v) for k, v in reps.items()}, names, groups)
        ran.append(scheme)
        if full:
            full_ran.append(scheme)
        if diag:
            diag_ran = True
            diag_rows += _diagnostic_rows(scheme, diag, thr_reps, names, samples)

    pos = _positivity_table(rows, extras, reference, configured)
    folds = _folds_table(pos, reference, block_folds)
    composition = _composition(a, rois, samples, smp, in_lesion, manuscript_pos, names, cfg)
    exclusion = _exclusion_table(a, excl, samples, smp, names, cfg)
    diagnostics = pd.DataFrame(diag_rows, columns=DIAGNOSTIC_COLUMNS)
    depth = pd.DataFrame(depth_rows)
    reps_text = f"mean over {n_reps} replicates" if n_reps > 1 else "one replicate"
    depth_text = "native depth"
    if ran:
        partial = [d for d in ran if d not in full_ran]
        depth_text += (f" and after {', '.join(ran)} depth matching, {reps_text}; every rule and region under "
                       f"{', '.join(full_ran) or 'none'}"
                       + (f", the manuscript rule and region only under {', '.join(partial)}" if partial else ""))
    block_text = (f"spatial block-bootstrap 95% intervals of the native manuscript-rule percentages "
                  f"({top['block_um']:g} um blocks, {int(top['n_boot'])} resamples, each section resampled with its "
                  f"own seed)" if blocks is not None else "no block-bootstrap intervals (no microns per pixel)")
    block_fold_text = ("; for the native manuscript-rule rows, the spatial block-bootstrap 95% interval (both "
                       "sections' blocks resampled independently, percentile interval of the ratio; the upper end NaN "
                       "when unbounded, i.e. the reference section has no positive bin in more than 2.5% of "
                       "resamples), which allows for the clustering of positive bins" if blocks is not None else "")
    pct_cut = f"{float(cfg['roi']['exclude_percentile']):g}"
    return [
        common.Output(POSITIVITY_CSV, "Positivity under other regions, rules and depth matching",
                      f"Positive bins per signature, section, region of interest ({', '.join(rois)}), calling rule "
                      f"({', '.join(rules)}) and depth ({depth_text}). Wilson 95% intervals (rounded mean count when "
                      f"depth-matched), {block_text}, and the candidates (score and core-gene rules passed) and "
                      f"competitor vetoes of the native manuscript rule. Regions are defined at native depth.", pos),
        common.Output(FOLDS_CSV, "Positivity folds under other regions, rules and depth matching",
                      f"Percentage positive in each non-reference section over the reference section's ({reference}) "
                      f"for every row of {POSITIVITY_CSV}, with the Katz 95% interval of a ratio of two independent "
                      f"proportions (bins treated as independent; on depth-matched rows from the replicate-mean "
                      f"counts, replicate variation not included){block_fold_text}. Intervals are NaN when a count "
                      f"is 0.", folds),
        common.Output(COMPOSITION_CSV, "Region-of-interest composition",
                      "Bins of each region of interest per section: share of the section, share assigned to skeletal "
                      "muscle, smooth muscle and nerve (Figure 3 assignment; empty without one), share in the A549 "
                      "lesion, bins removed and added against the manuscript region, and native manuscript-rule "
                      "positive bins.", composition),
        common.Output(EXCLUSION_CSV, "Region-of-interest exclusion rules",
                      f"Per exclusion signature and section, the rule that excludes bins from the manuscript region "
                      f"(exclusion_rule): for 'percentile', the {pct_cut}th-percentile score cut over all bins of all "
                      f"sections (percentile_value, and whether it is 0); for 'positivity call', the signature's pos_ "
                      f"calls (no percentile). Then the share of bins excluded, the share of those with a core gene "
                      f"detected, the share the positivity rule would exclude instead (call_threshold: its null "
                      f"threshold), and the native positive bins removed ('any': union of the signatures).",
                      exclusion),
        common.Output(DIAGNOSTICS_CSV, "Which part of the calling rule decides a call",
                      f"Region-of-interest bins per signature and section with at least {min_core} core genes, those "
                      f"of them removed by the null threshold, those vetoed by a competitor and the positive bins, at "
                      f"native depth" + (f" and after {configured} depth matching ({reps_text})" if diag_ran else "")
                      + ".", diagnostics),
        common.Output(DEPTH_CSV, "Depth of the region of interest under each matching scheme",
                      f"UMI per manuscript region-of-interest bin per section at native depth and after each "
                      f"depth-matching scheme ({', '.join(ran) or 'none'}; the per-bin targets), and the share of ROI "
                      f"bins downsampled.", depth),
    ]


def _depth_row(scheme, s, totals, matched, roi, smp):
    m = roi & (smp == s)
    n = int(m.sum())
    return {"depth": scheme, "sample": s, "n_roi_bins": n,
            "median_umi_native": float(np.median(totals[m])) if n else np.nan,
            "median_umi_matched": float(np.median(matched[m])) if n else np.nan,
            "mean_umi_matched": float(np.mean(matched[m])) if n else np.nan,
            "pct_downsampled": float(100 * np.mean(matched[m] < totals[m])) if n else np.nan}


def _composition(a, rois, samples, smp, in_lesion, pos, names, cfg):
    program = a.obs["program"].astype(str).to_numpy() if "program" in a.obs else None
    shares = list(MUSCLE) + [e for e in cfg["roi"]["exclude"] if e not in MUSCLE]
    focus = [n for n in FOCUS if n in names]
    base = rois["manuscript"]
    rows = []
    for roi, mask in rois.items():
        for s in samples:
            sec = smp == s
            m = mask & sec
            n = int(m.sum())
            row = {"roi": roi, "sample": s, "n_bins": n, "pct_of_section": 100 * n / max(int(sec.sum()), 1)}
            for e in shares:
                row[f"pct_{e}"] = 100 * float(np.mean(program[m] == e)) if program is not None and n else np.nan
            row["pct_in_lesion"] = 100 * float(np.mean(in_lesion[m])) if n else np.nan
            row["n_removed_vs_manuscript"] = int((base & ~mask & sec).sum())
            row["n_added_vs_manuscript"] = int((mask & ~base & sec).sum())
            for g in focus:
                row[f"n_pos_{g}"] = int(pos[m, names.index(g)].sum())
            rows.append(row)
    return pd.DataFrame(rows)


def _exclusion_table(a, excl, samples, smp, names, cfg):
    focus = [n for n in FOCUS if n in names]
    pos = {g: common.as_bool(a.obs[f"pos_{g}"]) for g in focus}
    ncore = {e: (a.obs[f"ncore_{e}"].to_numpy() >= 1) if f"ncore_{e}" in a.obs else np.zeros(a.n_obs, dtype=bool)
             for e in excl}
    entries = [(e, d["rule"], d["percentile"], d["threshold"], d["manuscript"], d["call"], d["manuscript"] & ncore[e])
               for e, d in excl.items()]
    if len(entries) > 1:
        union = np.zeros(a.n_obs, dtype=bool)
        call, core = union.copy(), union.copy()
        for e, d in excl.items():
            union |= d["manuscript"]
            call |= d["call"]
            core |= d["manuscript"] & ncore[e]
        entries.append(("any", "union", np.nan, np.nan, union, call, core))
    rows = []
    for e, rule, value, t, man, call, core in entries:
        for s in samples:
            sec = smp == s
            n, x = int(sec.sum()), int((man & sec).sum())
            row = {"signature": e, "exclusion_rule": rule, "percentile_value": value,
                   "percentile_is_zero": bool(value == 0) if np.isfinite(value) else pd.NA,
                   "call_threshold": t, "sample": s, "n_bins": n,
                   "pct_bins_excluded": 100 * x / n if n else np.nan,
                   "pct_excluded_with_core_gene": 100 * int((core & sec).sum()) / x if x else np.nan,
                   "pct_excluded_call_rule": 100 * int((call & sec).sum()) / n if n else np.nan}
            for g in focus:
                k, tot = int((man & sec & pos[g]).sum()), int((sec & pos[g]).sum())
                row[f"n_excluded_positive_{g}"] = k
                row[f"pct_excluded_positive_{g}"] = 100 * k / tot if tot else np.nan
            rows.append(row)
    df = pd.DataFrame(rows)
    if len(df):
        df["percentile_is_zero"] = df["percentile_is_zero"].astype("boolean")    # empty where no percentile applies
    return df
