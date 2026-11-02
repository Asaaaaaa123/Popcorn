"""
Figure 4i under other null models, call sets and depth-matching draws (second review, findings F1, G2.1, G2.2, G4.2).

Figure 4i (bone_analysis.uncoupling_metric) compares the mean distance from each osteoblast-called bin (A) to the
nearest osteoclast-called bin (B) with a null that places BOTH sets uniformly at random over the region of interest
(complete spatial randomness, csr_both). The pooled calls are 3x3 spatial smooths and both lineages sit on bone
surfaces, so a B pattern that is clumped, or confined to part of the ROI, makes observed/null > 1 on its own. The
ratio factorises exactly:
    obs / csr null = (obs / baseline) x (baseline / csr null),   baseline = mean distance from every ROI bin to B,
the first factor measured against the real osteoclast pattern, the second how far the B pattern departs from random
placement in the ROI at any scale (local clumping and confinement to part of the ROI; frac_b_clumped and
nn_self_ratio_b give the local clumping alone). This module repeats the test under nulls that keep the pattern of one
or both sets, for the three call sets (pooled depth-matched, pooled, native), and for several depth-matching draws
(the manuscript's pooled_dm calls are one draw, seed 0), and reports how the per-bin calls change between draws.

Null models (n_perm random configurations each, per section, within the ROI):
  csr_both         A and B placed at random (bone_analysis.uncoupling_metric itself: the manuscript test)
  fixed_b          B kept, nA bins drawn at random from the ROI
  fixed_a          A kept, nB bins drawn at random from the ROI (distance still from each A bin to its nearest B bin)
  habitat_fixed_b  B kept, A drawn from the bone habitat: ROI bins assigned (program_dm for pooled_dm calls, program
                   otherwise) to a signature of the habitat classes, plus every A or B bin. The pool always holds A,
                   so n_habitat >= nA; n_habitat_extra counts the habitat bins outside the calls, and the null is not
                   estimable when it is 0 (A would only be reshuffled over A or B).
  random_labels    the called positions kept; co-positive bins stay both, the A-only and B-only labels are permuted
                   over the other called positions (counts kept)
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import pipeline_config as pc
from bone_sensitivity import common

CALL_PREFIX = {"pooled_dm": "pooled_dm_", "pooled": "pooled_", "native": ""}
ASSIGNMENT = {"pooled_dm": "program_dm", "pooled": "program", "native": "program"}
NULL_MODELS = ("csr_both", "fixed_b", "fixed_a", "habitat_fixed_b", "random_labels")
N_SELF_SETS = 20            # random sets of the same size behind nn_self_ratio_a / _b
CLUMP_FACTOR = 1.5          # the 8 grid neighbours of a bin lie within 1.5 bin widths

NULL_COLUMNS = ["call_set", "draw", "null_model", "sample", "status", "n_a", "n_b", "n_units", "n_habitat",
                "n_habitat_extra", "obs_mean_um", "null_mean_um", "ratio", "ratio_lo", "ratio_hi", "p_greater",
                "p_less", "ratio_of_ratios_vs_reference", "n_perm"]
SUMMARY_COLUMNS = ["call_set", "null_model", "sample", "n_draws", "ratio_mean", "ratio_sd", "ratio_min", "ratio_max",
                   "ratio_of_ratios_mean", "ratio_of_ratios_sd", "ratio_of_ratios_min", "ratio_of_ratios_max",
                   "p_greater_max", "p_less_max", "manuscript"]
DECOMPOSITION_COLUMNS = ["call_set", "draw", "sample", "n_a", "n_b", "obs_mean_um", "baseline_all_roi_to_b_um",
                         "csr_null_mean_um", "factor_a_specific", "factor_b_pattern", "frac_a_clumped",
                         "frac_b_clumped", "nn_self_ratio_a", "nn_self_ratio_b"]
STABILITY_COLUMNS = ["call_set", "signature", "sample", "region", "n_draws", "n_bins", "n_called_draw0",
                     "n_called_min", "n_called_max", "n_ever_called", "n_always_called", "n_unstable",
                     "pct_unstable_of_ever_called", "pct_unstable_of_bins", "pct_changed_vs_draw0_mean",
                     "pct_changed_vs_draw0_max"]


def _info(progress, message):
    (progress.info if progress is not None and hasattr(progress, "info") else print)(message)


# ---- statistical core (numpy; positions in um, boolean masks over the ROI bins of one section) -------------------

def nearest_distance(xy, A, B):
    """Distance from each A point to its nearest B point."""
    return cKDTree(xy[B]).query(xy[A])[0]


def null_means(model, xy, A, B, n_perm, rng, habitat=None):
    """Mean A-to-nearest-B distance of n_perm random configurations under one null model (see the module docstring).
    habitat: boolean mask of the habitat_fixed_b draw pool (A and B bins are added to it)."""
    A, B = np.asarray(A, bool), np.asarray(B, bool)
    n, nA, nB = len(xy), int(A.sum()), int(B.sum())
    out = np.empty(int(n_perm))
    if model in ("fixed_b", "habitat_fixed_b"):
        to_b = cKDTree(xy[B]).query(xy)[0]
        pool = np.arange(n) if model == "fixed_b" else np.flatnonzero(np.asarray(habitat, bool) | A | B)
        for i in range(len(out)):
            out[i] = to_b[rng.choice(pool, nA, replace=False)].mean()
    elif model == "fixed_a":
        xa = xy[A]
        for i in range(len(out)):
            out[i] = cKDTree(xy[rng.choice(n, nB, replace=False)]).query(xa)[0].mean()
    elif model == "csr_both":                       # bone_analysis.uncoupling_metric's draws, in the same order
        for i in range(len(out)):
            ia = rng.choice(n, nA, replace=False)
            ib = rng.choice(n, nB, replace=False)
            out[i] = cKDTree(xy[ib]).query(xy[ia])[0].mean()
    elif model == "random_labels":
        for i in range(len(out)):
            out[i] = nearest_distance(xy, *random_labels_draw(A, B, rng)).mean()
    else:
        raise ValueError(f"unknown null model {model!r}")
    return out


def random_labels_draw(A, B, rng):
    """One random_labels configuration (A', B'): co-positive bins stay in both sets, the A-only and B-only labels are
    permuted over the remaining called positions."""
    A, B = np.asarray(A, bool), np.asarray(B, bool)
    free = np.flatnonzero(A ^ B)
    perm = rng.permutation(free)
    n_a_only = int((A & ~B).sum())
    A2, B2 = A & B, A & B
    A2[perm[:n_a_only]] = True
    B2[perm[n_a_only:]] = True
    return A2, B2


def null_statistics(obs, nulls):
    """ratio = obs / mean null, its 95% null interval (obs over the 97.5th and 2.5th null percentiles) and one-sided
    permutation P values (count + 1) / (n + 1)."""
    nulls = np.asarray(nulls, float)
    lo, hi = np.percentile(nulls, [2.5, 97.5])
    n = len(nulls)
    return {"null_mean_um": float(nulls.mean()), "ratio": float(common.fold(obs, nulls.mean())),
            "ratio_lo": float(common.fold(obs, hi)), "ratio_hi": float(common.fold(obs, lo)),
            "p_greater": float((np.sum(nulls >= obs) + 1) / (n + 1)),
            "p_less": float((np.sum(nulls <= obs) + 1) / (n + 1))}


def fraction_clumped(xy, M, radius):
    """Share of the points of M with another point of M within `radius` (8 grid neighbours: 1.5 bin widths)."""
    P = xy[np.asarray(M, bool)]
    if len(P) < 2:
        return np.nan
    d = cKDTree(P).query(P, k=2)[0][:, 1]
    return float(np.mean(d <= radius))


def _mean_nn_self(P):
    return float(cKDTree(P).query(P, k=2)[0][:, 1].mean())


def nn_self_ratio(xy, M, rng, n_sets=N_SELF_SETS):
    """Mean distance from each point of M to the nearest other point of M, over its mean across n_sets random sets of
    the same size drawn from all points (> 1 dispersed, < 1 clumped relative to the ROI)."""
    M = np.asarray(M, bool)
    k = int(M.sum())
    if k < 2:
        return np.nan
    ref = np.mean([_mean_nn_self(xy[rng.choice(len(xy), k, replace=False)]) for _ in range(int(n_sets))])
    return float(common.fold(_mean_nn_self(xy[M]), ref))


def decomposition(xy, A, B, csr_null_mean, bin_um, rng):
    """obs / csr null = factor_a_specific (obs / baseline) x factor_b_pattern (baseline / csr null), baseline = mean
    distance from every ROI bin to its nearest B bin; plus the clumping of A and B."""
    A, B = np.asarray(A, bool), np.asarray(B, bool)
    row = {"n_a": int(A.sum()), "n_b": int(B.sum())}
    if row["n_a"] and row["n_b"]:
        to_b = cKDTree(xy[B]).query(xy)[0]
        obs, base = float(to_b[A].mean()), float(to_b.mean())
    else:
        obs = base = np.nan
    row.update(obs_mean_um=obs, baseline_all_roi_to_b_um=base, csr_null_mean_um=float(csr_null_mean),
               factor_a_specific=float(common.fold(obs, base)),
               factor_b_pattern=float(common.fold(base, csr_null_mean)),
               frac_a_clumped=fraction_clumped(xy, A, CLUMP_FACTOR * bin_um),
               frac_b_clumped=fraction_clumped(xy, B, CLUMP_FACTOR * bin_um),
               nn_self_ratio_a=nn_self_ratio(xy, A, rng), nn_self_ratio_b=nn_self_ratio(xy, B, rng))
    return row


def call_stability(draws, smp, roi, samples):
    """How the pooled depth-matched calls change between downsampling draws, per signature, section and region (the
    marrow ROI of Figure 4i, and the whole section). n_unstable counts the bins called in some draws but not in all;
    pct_changed_vs_draw0 is, for each later draw, 100 x the bins whose call differs from draw 0 over the bins called in
    draw 0 or that draw (100 x (1 - Jaccard overlap)). draws: {signature: bool array (n_draws, n_obs)}."""
    rows = []
    for sig, D in draws.items():
        D = np.asarray(D, bool)
        for s in samples:
            for region, m in (("roi", (smp == s) & roi), ("section", smp == s)):
                X = D[:, m]
                ever, always = X.any(axis=0), X.all(axis=0)
                n_unstable, n_called = int((ever & ~always).sum()), X.sum(axis=1)
                changed = [float(common.fold(100.0 * np.sum(X[0] ^ X[d]), np.sum(X[0] | X[d])))
                           for d in range(1, len(X))]
                changed = np.asarray(changed, float)
                rows.append({"call_set": "pooled_dm", "signature": sig, "sample": s, "region": region,
                             "n_draws": int(len(X)), "n_bins": int(m.sum()), "n_called_draw0": int(n_called[0]),
                             "n_called_min": int(n_called.min()), "n_called_max": int(n_called.max()),
                             "n_ever_called": int(ever.sum()), "n_always_called": int(always.sum()),
                             "n_unstable": n_unstable,
                             "pct_unstable_of_ever_called": float(common.fold(100.0 * n_unstable, ever.sum())),
                             "pct_unstable_of_bins": float(common.fold(100.0 * n_unstable, m.sum())),
                             "pct_changed_vs_draw0_mean": (float(np.nanmean(changed)) if np.isfinite(changed).any()
                                                           else np.nan),
                             "pct_changed_vs_draw0_max": (float(np.nanmax(changed)) if np.isfinite(changed).any()
                                                          else np.nan)})
    return pd.DataFrame(rows, columns=STABILITY_COLUMNS)


# ---- depth-matching draws ------------------------------------------------------------------------------------------

def scoring_signatures(cfg, names):
    """The scoring dicts of `names` and of their competitors ({core, extended, competes, excluded}, as in _run)."""
    sig = cfg["signatures"]
    keep = []
    for n in names:
        for k in [n, *(sig[n].get("competes") or [])]:
            if k in sig and k not in keep:
                keep.append(k)
    return {k: {"core": list(sig[k]["core"]), "extended": list(sig[k]["extended"]),
                "competes": list(sig[k].get("competes") or []), "excluded": dict(sig[k].get("excluded") or {})}
            for k in keep}


def pooled_dm_call_draws(a, cfg, n_draws, progress=None):
    """Pooled, depth-matched calls of uncoupling.a and .b for n_draws downsampling draws (seed + d), each computed as
    bone_analysis.spatially_pooled_scores(depth_match=True) does: downsample, pool over the bin and its pool_k nearest
    bins, normalize, score, null thresholds, calls. Scores and thresholds are per signature with fresh seeds, so scoring
    only the two signatures and their competitors gives the full run's calls: draw 0 equals pooled_dm_pos_<a>/<b>.
    Returns {signature: bool array (n_draws, n_obs)}."""
    import bone_analysis as ba
    sco, unc = cfg["scoring"], cfg["uncoupling"]
    names = list(dict.fromkeys([unc["a"], unc["b"]]))
    sigs = scoring_signatures(cfg, names)
    seed = int(cfg["seed"])
    targets = ba.depth_match_targets(a, sco["depth_matching"])[0]
    neighbors = ba.spatial_neighbor_matrix(a, sco["pool_k"])
    out = {n: np.zeros((int(n_draws), a.n_obs), dtype=bool) for n in names}
    for d in range(int(n_draws)):
        _info(progress, f"[uncoupling] pooled depth-matched calls, draw {d + 1}/{n_draws} (seed {seed + d})")
        pooled = (neighbors @ ba.downsample_counts_matrix(a.layers["counts"], targets, seed + d)).tocsr()
        b = ba._lognorm_object(pooled, a, ["sample"], float(cfg["clustering"]["target_sum"]))
        del pooled
        b = ba.call_positive(ba.score_all(b, sigs, sco["ctrl_size"], seed), names, sigs=sigs, min_core=sco["min_core"],
                             n_perm=sco["n_perm_resampled"], pct=sco["null_percentile"], n_sub=sco["null_n_sub"],
                             ctrl_size=sco["ctrl_size"], seed=seed, progress=progress)
        for n in names:
            out[n][d] = common.as_bool(b.obs[f"pos_{n}"])
        del b
    return out


# ---- run -------------------------------------------------------------------------------------------------------------

def _call_sets(a, cfg, st, progress):
    """{call_set: [(draw, A over all bins, B over all bins), ...]} and the pooled_dm draws (or None)."""
    unc = cfg["uncoupling"]
    sets, draws = {}, None
    for cs in st["call_sets"]:
        if cs not in CALL_PREFIX:
            common.warn(progress, f"Figure 4i sensitivity: unknown call set {cs!r} skipped")
            continue
        if cs == "pooled_dm" and int(st["n_draws"]) > 0:
            draws = pooled_dm_call_draws(a, cfg, int(st["n_draws"]), progress)
            for n in (unc["a"], unc["b"]):
                col = f"pooled_dm_pos_{n}"
                if col in a.obs and not np.array_equal(draws[n][0], common.as_bool(a.obs[col])):
                    common.warn(progress, f"Figure 4i sensitivity: draw 0 of the pooled depth-matched {n} calls "
                                          f"differs from {col}")
            sets[cs] = [(d, draws[unc["a"]][d], draws[unc["b"]][d]) for d in range(int(st["n_draws"]))]
            continue
        if cs == "pooled_dm":
            common.warn(progress, "Figure 4i sensitivity: n_draws is 0, so the pooled depth-matched calls are skipped")
            continue
        cols = [f"{CALL_PREFIX[cs]}pos_{unc['a']}", f"{CALL_PREFIX[cs]}pos_{unc['b']}"]
        if not all(c in a.obs for c in cols):
            common.warn(progress, f"Figure 4i sensitivity: {cs} calls ({', '.join(cols)}) are not in the data; skipped")
            continue
        sets[cs] = [(0, common.as_bool(a.obs[cols[0]]), common.as_bool(a.obs[cols[1]]))]
    return sets, draws


def _habitat(a, cfg, st, call_set):
    """Bins assigned to a habitat-class signature (None when the assignment column is absent). The pooled_dm draws
    all use program_dm, the assignment of the pipeline's draw (seed 0); the called bins of each draw are added later."""
    col = ASSIGNMENT[call_set]
    if col not in a.obs:
        return None, col
    sigs = [k for k, d in cfg["signatures"].items() if d.get("class") in set(st["habitat_classes"])]
    return a.obs[col].astype(str).isin(sigs).to_numpy(), col


def _csr_with_p_less(light, A, B, xy_all, roi, smp, samples, unc, ref, seed, bin_um):
    """bone_analysis.uncoupling_metric on the light object, and the per-replicate null means of its draws (the same
    generator stream, replayed) for p_less; None per section when the replay does not reproduce its null mean."""
    import bone_analysis as ba
    light.obs["_a"], light.obs["_b"] = A, B
    summ, _ = ba.uncoupling_metric(light, "_a", "_b", bin_um=bin_um, n_perm=int(unc["n_perm"]),
                                   min_pos=int(unc["min_pos"]), reference_sample=ref, seed=seed, samples=samples)
    rng = np.random.default_rng(seed)
    nulls = {}
    for s in samples:
        m = (smp == s) & roi
        if int(A[m].sum()) < int(unc["min_pos"]) or int(B[m].sum()) < int(unc["min_pos"]):
            continue
        nulls[s] = null_means("csr_both", xy_all[m], A[m], B[m], int(unc["n_perm"]), rng)
    for s, v in list(nulls.items()):
        r = summ[summ["sample"] == s]
        if r.empty or not np.isclose(v.mean(), float(r["null_mean_um"].iloc[0]), rtol=1e-10, atol=0):
            nulls[s] = None
    return summ, nulls


def run(a, cfg, progress=None):
    """Figure 4i under the configured null models, call sets and depth-matching draws, and the stability of the pooled
    depth-matched calls between draws; adds obs columns pooled_dm_pos_frac_<a> and pooled_dm_pos_frac_<b> (fraction of
    draws in which each bin is called)."""
    st = common.settings(cfg, "uncoupling")
    unc, seed, ref = cfg["uncoupling"], int(cfg["seed"]), cfg["reference_sample"]
    bin_um = pc.bin_um(cfg["bin_size"])
    min_pos, n_perm = int(unc["min_pos"]), int(st["n_perm"])
    samples = common.sample_ids(a, cfg)
    models = []
    for m in st["null_models"]:
        if m in NULL_MODELS:
            models.append(m)
        else:
            common.warn(progress, f"Figure 4i sensitivity: unknown null model {m!r} skipped")
    sets, draws = _call_sets(a, cfg, st, progress)
    if draws is not None:
        for n in dict.fromkeys([unc["a"], unc["b"]]):
            a.obs[f"pooled_dm_pos_frac_{n}"] = draws[n].mean(axis=0)

    smp = a.obs["sample"].astype(str).to_numpy()
    roi = common.as_bool(a.obs["in_marrow_roi"])
    xy_all = common.array_xy_um(a, bin_um)
    light = ad.AnnData(obs=a.obs[["sample", "array_row", "array_col", "in_marrow_roi"]].copy())
    rows, dec = [], []
    for cs, entries in sets.items():
        habitat_all, habitat_col = _habitat(a, cfg, st, cs)
        for d, A_all, B_all in entries:
            _info(progress, f"[uncoupling] null models: {cs} calls, draw {d}")
            # always run: the decomposition needs the csr null mean
            csr, csr_nulls = _csr_with_p_less(light, A_all, B_all, xy_all, roi, smp, samples, unc, ref, seed, bin_um)
            if any(v is None for v in csr_nulls.values()):
                common.warn(progress, f"Figure 4i sensitivity ({cs}, draw {d}): the csr_both null draws could not be "
                                      f"replayed, so its p_less is missing")
            csr_by = {r["sample"]: r for r in csr.to_dict("records")}
            for s in samples:
                m = (smp == s) & roi
                xy, A, B = xy_all[m], A_all[m], B_all[m]
                nA, nB = int(A.sum()), int(B.sum())
                ok = nA >= min_pos and nB >= min_pos
                obs = float(nearest_distance(xy, A, B).mean()) if nA and nB else np.nan
                c = csr_by.get(s, {})
                dec.append({"call_set": cs, "draw": d, "sample": s,
                            **decomposition(xy, A, B, c.get("null_mean_um", np.nan) if ok else np.nan, bin_um,
                                            np.random.default_rng(seed))})
                base = {"call_set": cs, "draw": d, "sample": s, "n_a": nA, "n_b": nB, "n_units": int(m.sum()),
                        "n_habitat": np.nan, "n_habitat_extra": np.nan}
                for model in models:
                    row = {**base, "null_model": model, "n_perm": int(unc["n_perm"]) if model == "csr_both" else n_perm}
                    habitat = None
                    if model == "habitat_fixed_b" and habitat_all is not None:
                        habitat = habitat_all[m]
                        row["n_habitat"] = int((habitat | A | B).sum())
                        row["n_habitat_extra"] = int((habitat & ~A & ~B).sum())
                    if not ok:
                        row["status"] = f"not estimable (n_positive < {min_pos})"
                    elif model == "habitat_fixed_b" and habitat is None:
                        row["status"] = f"not estimable (no {habitat_col} assignment)"
                    elif model == "habitat_fixed_b" and row["n_habitat_extra"] == 0:
                        row["status"] = "not estimable (no habitat bins beyond the calls)"
                    elif model == "csr_both":
                        nulls = csr_nulls.get(s)
                        row.update(status=c.get("status"), obs_mean_um=c.get("mean_nn_um"),
                                   null_mean_um=c.get("null_mean_um"), ratio=c.get("ratio_obs_null_mean"),
                                   ratio_lo=c.get("ratio_obs_null_lo"), ratio_hi=c.get("ratio_obs_null_hi"),
                                   p_greater=c.get("p_greater_than_null_mean"),
                                   p_less=(np.nan if nulls is None else
                                           float((np.sum(nulls <= c.get("mean_nn_um")) + 1) / (len(nulls) + 1))))
                    else:
                        nulls = null_means(model, xy, A, B, n_perm, np.random.default_rng(seed), habitat)
                        row.update(status="estimable", obs_mean_um=obs, **null_statistics(obs, nulls))
                    rows.append(row)
    df = pd.DataFrame(rows).reindex(columns=NULL_COLUMNS)
    df = _ratio_of_ratios(df, ref)
    dec = pd.DataFrame(dec).reindex(columns=DECOMPOSITION_COLUMNS)
    summary = summarize(df, unc["calls"])
    stability = (call_stability(draws, smp, roi, samples) if draws is not None
                 else pd.DataFrame(columns=STABILITY_COLUMNS))
    a_lbl, b_lbl = pc.signature_label(cfg, unc["a"]), pc.signature_label(cfg, unc["b"])
    return [
        common.Output("sensitivity_Figure4i_null_models.csv", "Figure 4i under other null models",
                      f"Observed / null mean distance from each {a_lbl}-called bin to the nearest {b_lbl}-called bin "
                      f"per call set, depth-matching draw, null model and section (ROI); csr_both is the manuscript "
                      f"test; n_habitat is the habitat_fixed_b draw pool (always holds the calls), n_habitat_extra "
                      f"its bins outside the calls; P values treat bins as independent.", df),
        common.Output("sensitivity_Figure4i_summary.csv", "Figure 4i null models, summary over draws",
                      "Ratio and ratio of ratios vs the reference section per call set, null model and section: mean, "
                      "SD and range over the depth-matching draws; manuscript marks the Figure 4i setup.", summary),
        common.Output("sensitivity_Figure4i_decomposition.csv", "Figure 4i ratio decomposition",
                      f"The csr_both ratio as factor_a_specific (observed / mean distance from every ROI bin to the "
                      f"nearest {b_lbl}-called bin) x factor_b_pattern (that baseline / csr null mean: how far the "
                      f"{b_lbl} pattern departs from random placement in the ROI, by local clumping or by confinement "
                      f"to part of it), and the local clumping of both call sets (frac_*_clumped, nn_self_ratio_*).",
                      dec),
        common.Output("sensitivity_Figure4i_call_stability.csv", "Figure 4i calls between depth-matching draws",
                      f"Pooled depth-matched {a_lbl} and {b_lbl} calls over the downsampling draws, per section, in "
                      f"the marrow ROI and the whole section: bins called in any, all, or some but not all draws "
                      f"(n_unstable), and per later draw the % of bins whose call differs from draw 0 among those "
                      f"called in either (empty when the pooled_dm call set was not run).", stability),
    ]


def _ratio_of_ratios(df, ref):
    """ratio(section) / ratio(reference) per call set, draw and null model, when both are estimable."""
    df = df.copy()
    ok = df["status"].astype(str).eq("estimable")
    key = ["call_set", "draw", "null_model"]
    r = df[ok & df["sample"].eq(ref)].set_index(key)["ratio"]
    idx = pd.MultiIndex.from_frame(df[key])
    ref_ratio = r.reindex(idx).to_numpy(float) if len(r) else np.full(len(df), np.nan)
    df["ratio_of_ratios_vs_reference"] = np.where(ok, common.fold(df["ratio"].to_numpy(float), ref_ratio), np.nan)
    return df


def summarize(df, manuscript_calls):
    """Mean, SD and range over draws of the ratio and the ratio of ratios (estimable rows) per call set, null model and
    section; n_draws counts the draws in which the section was estimable."""
    out = []
    for (cs, model, s), g in df.groupby(["call_set", "null_model", "sample"], sort=False):
        e = g[g["status"].astype(str).eq("estimable")]
        r, rr = e["ratio"].astype(float), e["ratio_of_ratios_vs_reference"].astype(float)
        out.append({"call_set": cs, "null_model": model, "sample": s, "n_draws": int(len(e)),
                    "ratio_mean": r.mean(), "ratio_sd": r.std(), "ratio_min": r.min(), "ratio_max": r.max(),
                    "ratio_of_ratios_mean": rr.mean(), "ratio_of_ratios_sd": rr.std(),
                    "ratio_of_ratios_min": rr.min(), "ratio_of_ratios_max": rr.max(),
                    "p_greater_max": e["p_greater"].astype(float).max(), "p_less_max": e["p_less"].astype(float).max(),
                    "manuscript": bool(cs == manuscript_calls and model == "csr_both")})
    return pd.DataFrame(out, columns=SUMMARY_COLUMNS)
