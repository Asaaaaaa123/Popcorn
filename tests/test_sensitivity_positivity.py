"""Positivity sensitivity (bone_sensitivity.positivity, positivity_figures): statistical core on synthetic data with
known answers, run() on the tiny demo (reproduction of the manuscript calls and depth-matched counts), the figure."""
from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from bone_sensitivity import common
from bone_sensitivity import positivity as P
from bone_sensitivity import positivity_figures as PF


class _Progress:
    def __init__(self):
        self.warnings, self.messages = [], []

    def warn(self, message):
        self.warnings.append(message)

    def info(self, message):
        self.messages.append(message)


# ---- statistical core ---------------------------------------------------------------------------------------------

def test_katz_interval_known_value():
    ratio, lo, hi = P.katz_interval(20, 100, 10, 100)
    se = np.sqrt(1 / 20 - 1 / 100 + 1 / 10 - 1 / 100)
    assert ratio == pytest.approx(2.0)
    assert lo == pytest.approx(2 * np.exp(-common.Z95 * se)) and lo == pytest.approx(0.98652, abs=1e-4)
    assert hi == pytest.approx(2 * np.exp(common.Z95 * se)) and hi == pytest.approx(4.0546, abs=1e-3)
    ratio, lo, hi = P.katz_interval([0, 5, 5], [100, 100, 0], [10, 0, 10], [100, 100, 100])
    assert ratio[0] == 0 and np.isnan(ratio[1]) and np.isnan(ratio[2])
    assert np.isnan(lo).all() and np.isnan(hi).all()


def test_rule_parameters():
    assert P.rule_parameters("manuscript", 2) == (2, True)
    assert P.rule_parameters("min_core_1", 2) == (1, True)
    assert P.rule_parameters("min_core_3", 2) == (3, True)
    assert P.rule_parameters("no_competitor_veto", 2) == (2, False)
    with pytest.raises(ValueError):
        P.rule_parameters("strict", 2)


def _synthetic_obs(n=2000, seed=1):
    rng = np.random.default_rng(seed)
    obs = pd.DataFrame({"score_x": rng.normal(0, 1, n), "ncore_x": rng.integers(0, 5, n),
                        "score_y": rng.normal(0, 1, n), "ncore_y": rng.integers(0, 5, n),
                        "score_c": rng.normal(0, 1, n)})
    obs.loc[obs.index[:50], "score_x"] = np.nan
    obs.loc[obs.index[50:100], "score_c"] = np.nan
    sigs = {"x": {"competes": ["c", "y", "absent"]}, "y": {"competes": []}, "c": {"competes": []}}
    return obs, sigs


def test_calls_follow_the_call_positive_rule():
    obs, sigs = _synthetic_obs()
    calls = P.Calls(obs, ["x", "y"], sigs, {"x": 0.2, "y": None})
    expect = (obs["score_x"] > 0.2) & (obs["ncore_x"] >= 2)
    for comp in ("c", "y"):                                       # as bone_analysis.call_positive, one by one
        expect &= ~(obs[f"score_{comp}"] > obs["score_x"])
    pos = calls.positive(2)
    assert np.array_equal(pos[:, 0], expect.fillna(False).to_numpy())
    assert not pos[:50, 0].any()                                  # a NaN score is never positive
    assert not pos[:, 1].any()                                    # no threshold: nothing positive


def test_rule_variants_are_monotone():
    obs, sigs = _synthetic_obs()
    calls = P.Calls(obs, ["x", "y"], sigs, {"x": 0.0, "y": -0.5})
    man = calls.positive(*P.rule_parameters("manuscript", 2))
    one = calls.positive(*P.rule_parameters("min_core_1", 2))
    three = calls.positive(*P.rule_parameters("min_core_3", 2))
    noveto = calls.positive(*P.rule_parameters("no_competitor_veto", 2))
    assert (one >= man).all() and (man >= three).all() and (noveto >= man).all()
    assert one.sum() > man.sum() > three.sum() and noveto[:, 0].sum() > man[:, 0].sum()
    d = calls.diagnostics(np.arange(len(obs)), 2)
    assert np.array_equal(d["n_positive"], man.sum(axis=0))
    assert np.array_equal(d["n_ncore_ge_min"] - d["n_removed_by_threshold"] - d["n_vetoed"], d["n_positive"])


def test_quantile_roi_targets_match_the_roi_distributions():
    rng = np.random.default_rng(0)
    # shallow section A: ROI bins around 50 UMI, non-ROI bins around 120; deeper section B around 200 everywhere
    a_roi, a_out, b = rng.poisson(50, 600) + 1, rng.poisson(120, 300), rng.poisson(200, 1000) + 60
    totals = np.concatenate([a_roi, a_out, b]).astype(float)
    smp = np.array(["A"] * 900 + ["B"] * 1000)
    roi = np.concatenate([np.ones(600, bool), np.zeros(300, bool), rng.random(1000) < 0.7])
    targets, shallow = P.quantile_roi_targets(totals, smp, roi)
    assert shallow == "A"
    assert np.array_equal(targets[smp == "A"], totals[smp == "A"])            # the reference keeps its counts
    assert np.array_equal(targets[(smp == "B") & ~roi], totals[(smp == "B") & ~roi])
    m = (smp == "B") & roi
    assert (targets[m] <= totals[m]).all()
    n = int(m.sum())
    expect = np.floor(np.quantile(np.sort(a_roi.astype(float)), (np.arange(n) + 0.5) / n))
    assert np.array_equal(np.sort(targets[m]), expect)                           # B's ROI gets A's ROI quantiles
    assert abs(np.median(targets[m]) - np.median(a_roi)) <= 1 and abs(targets[m].mean() - a_roi.mean()) < 1
    order = np.argsort(totals[m], kind="stable")                                 # deeper bins get higher targets
    assert (np.diff(targets[m][order]) >= 0).all()
    # whole-section quantile matching (bone_analysis) leaves B's ROI deeper than A's ROI; quantile_roi does not
    import bone_analysis as ba
    obj = ad.AnnData(obs=pd.DataFrame({"sample": pd.Categorical(smp), "total_counts": totals},
                                      index=[str(i) for i in range(len(smp))]))
    whole = ba.depth_match_targets(obj, "quantile")[0]
    assert whole[m].mean() > 1.3 * a_roi.mean()
    targets_none, _ = P.quantile_roi_targets(totals, smp, np.zeros(len(smp), bool))
    assert targets_none is None


def test_block_intervals_equal_the_point_bootstrap_and_widen_for_clusters():
    xy = np.array([(i, j) for i in range(20) for j in range(20)], float) * 10.0
    blocks = common.spatial_blocks(xy, 50.0)                     # 16 blocks of 25 bins
    clustered = np.zeros((400, 1), bool)
    clustered[np.isin(blocks, [0, 5]), 0] = True                 # two whole blocks positive (12.5%)
    lo, hi, draws = P.block_intervals(clustered, blocks, n_boot=400, seed=3)
    plo, phi, pdraws = common.block_bootstrap(lambda idx: 100 * clustered[idx].mean(axis=0), blocks, n_boot=400,
                                              seed=3)
    assert np.allclose(lo, plo) and np.allclose(hi, phi)          # block sums give the bin-level bootstrap exactly
    assert draws.shape == (400, 1) and np.allclose(draws, pdraws)
    wlo, whi = common.wilson(50, 400)
    assert lo[0] <= 12.5 <= hi[0]
    assert (hi[0] - lo[0]) > 2 * (whi - wlo)


def test_block_fold_interval_known_values():
    d = np.arange(1, 201, dtype=float)[:, None]
    lo, hi = P.block_fold_interval(d, np.full((200, 1), 2.0))
    assert lo[0] == pytest.approx(np.percentile(d[:, 0] / 2, 2.5))
    assert hi[0] == pytest.approx(np.percentile(d[:, 0] / 2, 97.5))
    ref = np.full((200, 1), 2.0)
    ref[:10] = 0                                                  # 5% of draws without a positive reference bin
    lo, hi = P.block_fold_interval(d, ref)
    assert np.isfinite(lo[0]) and np.isnan(hi[0])                 # unbounded above
    ref[3:10] = 2.0                                               # 1.5%: still bounded, by a finite draw
    lo, hi = P.block_fold_interval(d, ref)
    assert np.isfinite(hi[0]) and hi[0] <= 100
    zero = d.copy()
    zero[:10] = 0
    ref[:10] = 0                                                  # 0 / 0 draws are dropped, not unbounded
    lo, hi = P.block_fold_interval(zero, ref)
    assert hi[0] == pytest.approx(np.percentile(d[10:, 0] / 2, 97.5))
    for empty in (np.empty((0, 1)), None):
        lo, hi = P.block_fold_interval(d, empty)
        assert lo.shape == (1,) and np.isnan(lo).all() and np.isnan(hi).all()


def test_block_fold_interval_widens_for_clustered_sections():
    xy = np.array([(i, j) for i in range(20) for j in range(20)], float) * 10.0
    blocks = common.spatial_blocks(xy, 50.0)                     # 16 blocks of 25 bins
    tumor = np.isin(blocks, [0, 5, 10, 15])[:, None]             # 4 whole blocks positive: 25%
    ref = np.zeros((400, 1), bool)
    for b in (1, 6, 11, 12):                                      # 6 positive bins in each of 4 blocks: 6%
        ref[np.flatnonzero(blocks == b)[:6], 0] = True
    _, _, dt = P.block_intervals(tumor, blocks, n_boot=1000, seed=0)
    _, _, dr = P.block_intervals(ref, blocks, n_boot=1000, seed=1)
    lo, hi = P.block_fold_interval(dt, dr)
    fold, klo, khi = P.katz_interval(100, 400, 24, 400)
    assert fold == pytest.approx(25 / 6)
    assert lo[0] < klo and (np.isnan(hi[0]) or hi[0] > khi)       # wider than the bin-level Katz interval
    assert lo[0] <= fold


# ---- run() on the tiny demo -------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def positivity_run(tiny_joint):
    a = tiny_joint["a"].copy()
    obs_before, uns_before = a.obs.copy(), set(a.uns.keys())
    prog = _Progress()
    outs = P.run(a, tiny_joint["cfg"], prog)
    return {"a": a, "outs": {o.file: o for o in outs}, "obs_before": obs_before, "uns_before": uns_before,
            "progress": prog}


def test_run_outputs_and_columns(positivity_run):
    outs = positivity_run["outs"]
    assert set(outs) == {P.POSITIVITY_CSV, P.FOLDS_CSV, P.COMPOSITION_CSV, P.EXCLUSION_CSV, P.DIAGNOSTICS_CSV,
                         P.DEPTH_CSV}
    for o in outs.values():
        assert isinstance(o, common.Output) and o.title and o.description and len(o.df)
    pos = outs[P.POSITIVITY_CSV].df
    assert list(pos.columns) == P.POSITIVITY_COLUMNS
    assert list(outs[P.FOLDS_CSV].df.columns) == P.FOLD_COLUMNS
    assert list(outs[P.DIAGNOSTICS_CSV].df.columns) == P.DIAGNOSTIC_COLUMNS
    comp = outs[P.COMPOSITION_CSV].df
    for c in ["roi", "sample", "n_bins", "pct_of_section", "pct_skeletal_muscle", "pct_smooth_muscle", "pct_nerve",
              "pct_in_lesion", "n_removed_vs_manuscript", "n_added_vs_manuscript", "n_pos_osteoclast",
              "n_pos_osteoblast", "n_pos_macrophage"]:
        assert c in comp.columns
    exc = outs[P.EXCLUSION_CSV].df
    for c in ["signature", "percentile_value", "percentile_is_zero", "sample", "pct_bins_excluded",
              "pct_excluded_with_core_gene", "n_excluded_positive_osteoclast", "n_excluded_positive_osteoblast",
              "n_excluded_positive_macrophage"]:
        assert c in exc.columns
    assert set(pos["depth"]) == {"native", "cap", "quantile", "quantile_roi"}
    assert set(pos["roi"]) == set(P.ROI_VARIANTS) and set(pos["rule"]) == {"manuscript", "min_core_1", "min_core_3",
                                                                           "no_competitor_veto"}
    # every rule x region under native and cap, the manuscript analysis only under the other schemes
    combos = pos.groupby("depth")[["rule", "roi"]].apply(lambda d: len(d.drop_duplicates()))
    assert combos["native"] == combos["cap"] == 16 and combos["quantile"] == combos["quantile_roi"] == 1
    assert pos["manuscript"].sum() == 2 * 2 * pos["signature"].nunique()


def test_run_leaves_the_object_unchanged(positivity_run):
    a, before = positivity_run["a"], positivity_run["obs_before"]
    assert list(a.obs.columns) == list(before.columns)
    pd.testing.assert_frame_equal(a.obs, before)
    assert set(a.uns.keys()) == positivity_run["uns_before"]


def test_manuscript_rule_reproduces_the_calls(tiny_joint, positivity_run):
    a, cfg = tiny_joint["a"], tiny_joint["cfg"]
    names = [n for n in cfg["primary"] if f"pos_{n}" in a.obs]
    calls = P.Calls(a.obs, names, P.signature_dicts(cfg), a.uns["positivity_thresholds"])
    pos = calls.positive(*P.rule_parameters("manuscript", cfg["scoring"]["min_core"]))
    for j, n in enumerate(names):
        assert np.array_equal(pos[:, j], a.obs[f"pos_{n}"].to_numpy(bool)), n
    assert not any("differs from pos_" in w for w in positivity_run["progress"].warnings)
    # native manuscript rows = the manuscript table (Figures 4g, 6a, left)
    t = positivity_run["outs"][P.POSITIVITY_CSV].df
    mine = t[(t["depth"] == "native") & t["manuscript"]].set_index(["signature", "sample"])
    ref = pd.read_csv(tiny_joint["out"] / "Figure4g_6a_positivity_native_depth_roi.csv").set_index(["signature",
                                                                                                    "sample"])
    assert len(mine) == len(ref)
    assert np.array_equal(mine.loc[ref.index, "n_positive"].to_numpy(), ref["n_positive"].to_numpy())
    assert np.array_equal(mine.loc[ref.index, "n_bins"].to_numpy(), ref["n_bins"].to_numpy())
    assert np.allclose(mine.loc[ref.index, "pct_positive"], 100 * ref["n_positive"] / ref["n_bins"])
    assert np.allclose(mine.loc[ref.index, "ci_lo"], ref["ci_lo"], atol=1e-3)
    sub = mine.dropna(subset=["n_candidates"])
    assert len(sub) == len(mine)
    assert np.array_equal((sub["n_candidates"] - sub["n_vetoed"]).astype(float).to_numpy(),
                          sub["n_positive"].to_numpy(float))


def test_configured_depth_matching_reproduces_the_manuscript_table(tiny_joint, positivity_run):
    cfg = tiny_joint["cfg"]
    scheme = cfg["scoring"]["depth_matching"]
    t = positivity_run["outs"][P.POSITIVITY_CSV].df
    mine = t[(t["depth"] == scheme) & t["manuscript"]].set_index(["signature", "sample"])
    ref = pd.read_csv(tiny_joint["out"] / "Figure4g_6a_positivity_depth_matched.csv").set_index(["signature",
                                                                                                  "sample"])
    assert len(mine) == len(ref)
    assert np.allclose(mine.loc[ref.index, "n_positive"], ref["n_positive_mean"], atol=0.006)
    assert np.allclose(mine.loc[ref.index, "pct_positive"], ref["pct_positive_mean"], atol=1e-4)
    diag = positivity_run["outs"][P.DIAGNOSTICS_CSV].df
    d = diag[diag["depth"] == scheme].set_index(["signature", "sample"])
    assert np.allclose(d.loc[ref.index, "n_positive"], ref["n_positive_mean"], atol=0.006)


def test_null_thresholds_equal_the_pipeline(tiny_joint):
    import bone_analysis as ba
    a, cfg = tiny_joint["a"], tiny_joint["cfg"]
    sco, sigs = cfg["scoring"], P.signature_dicts(cfg)
    names = list(cfg["primary"])
    got = P.null_thresholds(a, names, sigs, sco["null_percentile"], 3, 500, sco["ctrl_size"], 7)
    for n in names:
        assert got[n] == ba.null_threshold(a, n, sigs=sigs, pct=sco["null_percentile"], n_perm=3, n_sub=500,
                                           ctrl_size=sco["ctrl_size"], seed=7), n


def test_rule_variants_monotone_on_tiny(positivity_run):
    t = positivity_run["outs"][P.POSITIVITY_CSV].df
    w = t.pivot_table(index=["depth", "roi", "signature", "sample"], columns="rule", values="n_positive")
    w = w.dropna()
    assert len(w)
    assert (w["min_core_1"] >= w["manuscript"] - 1e-9).all()
    assert (w["manuscript"] >= w["min_core_3"] - 1e-9).all()
    assert (w["no_competitor_veto"] >= w["manuscript"] - 1e-9).all()


def test_folds_and_intervals(tiny_joint, positivity_run):
    ref = tiny_joint["cfg"]["reference_sample"]
    t = positivity_run["outs"][P.POSITIVITY_CSV].df
    f = positivity_run["outs"][P.FOLDS_CSV].df
    assert ref not in set(f["sample"])
    key = ["depth", "rule", "roi", "signature"]
    m = f.merge(t[t["sample"] == ref][key + ["pct_positive"]], on=key).merge(
        t[t["sample"] != ref][key + ["sample", "pct_positive"]], on=key + ["sample"], suffixes=("_ref", ""))
    ok = m["pct_positive_ref"] > 0
    assert np.allclose(m.loc[ok, "fold"], m.loc[ok, "pct_positive"] / m.loc[ok, "pct_positive_ref"], rtol=1e-3)
    zero = (f["n_positive"] <= 0) | (f["n_positive_reference"] <= 0)
    assert f.loc[zero, "fold_ci_lo"].isna().all()
    fin = f["fold_ci_lo"].notna()
    assert (f.loc[fin, "fold_ci_lo"] <= f.loc[fin, "fold"] + 1e-12).all()
    assert (f.loc[fin, "fold"] <= f.loc[fin, "fold_ci_hi"] + 1e-12).all()
    # block intervals only on native manuscript-rule rows; wider than Wilson on average there
    has = t["block_ci_lo"].notna()
    assert ((t.loc[has, "depth"] == "native") & (t.loc[has, "rule"] == "manuscript")).all()
    native = t[(t["depth"] == "native") & (t["rule"] == "manuscript")]
    assert native["block_ci_lo"].notna().all()
    dm = t[t["depth"] != "native"]
    assert (dm["n_positive_sd"].notna()).all() and t[t["depth"] == "native"]["n_positive_sd"].isna().all()
    # block-bootstrap fold intervals: native manuscript-rule rows only, NaN when a count is 0
    has = f["block_fold_ci_lo"].notna()
    assert has.any()
    assert ((f.loc[has, "depth"] == "native") & (f.loc[has, "rule"] == "manuscript")).all()
    assert f.loc[zero, ["block_fold_ci_lo", "block_fold_ci_hi"]].isna().all().all()
    both = has & f["block_fold_ci_hi"].notna()
    assert (f.loc[both, "block_fold_ci_lo"] <= f.loc[both, "block_fold_ci_hi"]).all()


def test_block_fold_intervals_reproduce(tiny_joint, positivity_run):
    """The block intervals of a native manuscript row, recomputed from the calls, the section blocks and the
    per-section seeds."""
    a, cfg = tiny_joint["a"], tiny_joint["cfg"]
    top, seed, ref = common.settings(cfg, None), int(cfg["seed"]), cfg["reference_sample"]
    samples = common.sample_ids(a, cfg)
    smp = a.obs["sample"].astype(str).to_numpy()
    blocks = P._section_blocks(a, smp, float(top["block_um"]), None)
    roi = a.obs["in_marrow_roi"].to_numpy(bool)
    t = positivity_run["outs"][P.POSITIVITY_CSV].df
    f = positivity_run["outs"][P.FOLDS_CSV].df
    checked = 0
    for sig in ("osteoclast", "macrophage"):
        pos = a.obs[f"pos_{sig}"].to_numpy(bool)[:, None]
        res = {}
        for s in samples:
            idx = np.flatnonzero(roi & (smp == s))
            res[s] = P.block_intervals(pos[idx], blocks[idx], int(top["n_boot"]), seed + samples.index(s))
            row = t[(t["depth"] == "native") & (t["rule"] == "manuscript") & (t["roi"] == "manuscript")
                    & (t["signature"] == sig) & (t["sample"] == s)].iloc[0]
            assert row["block_ci_lo"] == pytest.approx(res[s][0][0], abs=1e-4)
            assert row["block_ci_hi"] == pytest.approx(res[s][1][0], abs=1e-4)
        for s in samples:
            if s == ref:
                continue
            row = f[(f["depth"] == "native") & (f["rule"] == "manuscript") & (f["roi"] == "manuscript")
                    & (f["signature"] == sig) & (f["sample"] == s)].iloc[0]
            if row["n_positive"] > 0 and row["n_positive_reference"] > 0:
                lo, hi = P.block_fold_interval(res[s][2], res[ref][2])
                assert np.allclose([row["block_fold_ci_lo"], row["block_fold_ci_hi"]], [lo[0], hi[0]], equal_nan=True)
                checked += 1
    assert checked


def test_regions_of_interest(tiny_joint, positivity_run):
    a, cfg = tiny_joint["a"], tiny_joint["cfg"]
    comp = positivity_run["outs"][P.COMPOSITION_CSV].df.set_index(["roi", "sample"])
    roi_n = a.obs.groupby("sample", observed=True)["in_marrow_roi"].sum()
    for s, n in roi_n.items():
        assert comp.loc[("manuscript", s), "n_bins"] == n
        assert comp.loc[("manuscript", s), "n_removed_vs_manuscript"] == 0
    # assignment shares against an independent crosstab of the tissue assignment, and the assignment_exclusion region
    # against the assignment labels (they must match the roi.exclude names, or the region would remove nothing)
    exclude = list(cfg["roi"]["exclude"])
    program = a.obs["program"].astype(str)
    assert set(exclude) & set(program)
    in_roi = a.obs["in_marrow_roi"].to_numpy(bool)
    tab = pd.crosstab(a.obs["sample"].astype(str)[in_roi], program[in_roi], normalize="index") * 100
    for s in roi_n.index:
        for e in ("skeletal_muscle", "smooth_muscle", "nerve"):
            expect = tab.loc[s, e] if e in tab.columns else 0.0
            assert comp.loc[("manuscript", s), f"pct_{e}"] == pytest.approx(expect)
        kept = (a.obs["sample"].astype(str) == s) & ~program.isin(exclude)
        assert comp.loc[("assignment_exclusion", s), "n_bins"] == int(kept.sum())
    shares = [f"pct_{e}" for e in exclude]
    assert comp.xs("manuscript", level="roi")[shares].to_numpy().sum() > 0
    assert (comp.xs("assignment_exclusion", level="roi")[shares] == 0).all().all()
    sigs = P.signature_dicts(cfg)
    excl = P.exclusion_rules(a, cfg, sigs)
    union = np.zeros(a.n_obs, bool)
    for d in excl.values():
        union |= d["manuscript"]
    assert np.array_equal(~union, a.obs["in_roi_transcriptional"].to_numpy(bool))  # define_marrow_roi, recomputed
    # in_lesion and in_roi_drawn are honoured; a missing assignment skips its region with a warning
    b = a.copy()
    rng = np.random.default_rng(0)
    b.obs["in_lesion"] = rng.random(b.n_obs) < 0.2
    b.obs["in_roi_drawn"] = b.obs["sample"].astype(str).to_numpy() == cfg["reference_sample"]
    lesion = P.lesion_flags(b, cfg)
    assert np.array_equal(lesion, b.obs["in_lesion"].to_numpy())
    prog = _Progress()
    del b.obs["program"]
    rois = P.roi_variants(b, cfg, P.ROI_VARIANTS, excl, lesion, prog)
    assert "assignment_exclusion" not in rois and any("program" in w for w in prog.warnings)
    base = b.obs["in_marrow_roi"].to_numpy(bool)
    assert np.array_equal(rois["no_lesion"], base & ~lesion)
    assert not rois["call_rule_exclusion"][~b.obs["in_roi_drawn"].to_numpy()].any()


def test_call_rule_exclusion_equals_call_positive(tiny_joint):
    """The call rule of call_rule_exclusion is bone_analysis.call_positive without competitors at native n_perm."""
    import bone_analysis as ba
    a, cfg = tiny_joint["a"], tiny_joint["cfg"]
    sco, sigs = cfg["scoring"], P.signature_dicts(cfg)
    excl = P.exclusion_rules(a, cfg, sigs)
    checked = 0
    for e, d in excl.items():
        if d["rule"] != "percentile":
            continue
        b = a.copy()
        alone = {k: dict(v) for k, v in sigs.items()}
        alone[e]["competes"] = []
        ba.call_positive(b, [e], sigs=alone, min_core=sco["min_core"], n_perm=sco["n_perm"],
                         pct=sco["null_percentile"], n_sub=sco["null_n_sub"], ctrl_size=sco["ctrl_size"],
                         seed=int(cfg["seed"]))
        assert np.array_equal(d["call"], b.obs[f"pos_{e}"].to_numpy(bool)), e
        assert d["threshold"] == b.uns["positivity_thresholds"][e], e
        checked += 1
        del b
    assert checked == len(cfg["roi"]["exclude"])


def test_quantile_scheme_equals_depth_matched(tiny_joint, positivity_run):
    """The quantile rows (manuscript rule and region) equal bone_analysis.depth_matched(mode="quantile")."""
    import bone_analysis as ba
    a, cfg = tiny_joint["a"], tiny_joint["cfg"]
    sco = cfg["scoring"]
    dm = ba.depth_matched(a.copy(), list(cfg["primary"]), n_reps=sco["n_depth_reps"], roi_key="in_marrow_roi",
                          mode="quantile", n_perm=sco["n_perm_resampled"],
                          target_sum=float(cfg["clustering"]["target_sum"]), samples=common.sample_ids(a, cfg),
                          sigs=P.signature_dicts(cfg), min_core=sco["min_core"], pct=sco["null_percentile"],
                          n_sub=sco["null_n_sub"], ctrl_size=sco["ctrl_size"], seed=int(cfg["seed"]))
    t = positivity_run["outs"][P.POSITIVITY_CSV].df
    mine = t[(t["depth"] == "quantile") & (t["rule"] == "manuscript") & (t["roi"] == "manuscript")]
    m = mine.merge(dm, on=["signature", "sample"])
    assert len(m) == len(mine) == len(dm)
    assert np.allclose(m["n_positive"], m["n_positive_mean"], atol=0.006)


def test_descriptions_name_what_ran(tiny_joint):
    """With no depth-matching scheme and no tissue assignment, the descriptions name only what was computed."""
    import copy
    a, cfg = tiny_joint["a"].copy(), copy.deepcopy(tiny_joint["cfg"])
    del a.obs["program"]
    cfg.setdefault("sensitivity_analyses", {}).setdefault("positivity", {})["depth_schemes"] = []
    outs = {o.file: o for o in P.run(a, cfg, _Progress())}
    pos = outs[P.POSITIVITY_CSV]
    assert set(pos.df["depth"]) == {"native"} and "assignment_exclusion" not in set(pos.df["roi"])
    for word in ("quantile", "depth matching", "assignment_exclusion"):
        assert word not in pos.description, word
    assert "no_lesion" in pos.description and "min_core_3" in pos.description
    assert "depth matching" not in outs[P.DIAGNOSTICS_CSV].description
    assert set(outs[P.DIAGNOSTICS_CSV].df["depth"]) == {"native"}


def test_exclusion_by_a_positivity_call_has_no_percentile(tiny_joint):
    """A roi.exclude signature with pos_ calls is excluded through them (define_marrow_roi): no percentile reported."""
    import copy
    a, cfg = tiny_joint["a"], copy.deepcopy(tiny_joint["cfg"])
    cfg["roi"]["exclude"] = list(cfg["roi"]["exclude"]) + ["osteoclast"]
    sigs = P.signature_dicts(cfg)
    excl = P.exclusion_rules(a, cfg, sigs)
    assert excl["osteoclast"]["rule"] == "positivity call" and np.isnan(excl["osteoclast"]["percentile"])
    assert np.array_equal(excl["osteoclast"]["manuscript"], a.obs["pos_osteoclast"].to_numpy(bool))
    names = [n for n in cfg["primary"] if f"pos_{n}" in a.obs]
    t = P._exclusion_table(a, excl, common.sample_ids(a, cfg), a.obs["sample"].astype(str).to_numpy(), names, cfg)
    by = t.set_index("signature")
    assert by.loc["osteoclast", "percentile_value"].isna().all()
    assert by.loc["osteoclast", "percentile_is_zero"].isna().all() and by.loc["any", "percentile_is_zero"].isna().all()
    pct = t[t["exclusion_rule"] == "percentile"]
    assert pct["percentile_value"].notna().all() and pct["percentile_is_zero"].notna().all()


# ---- figure -------------------------------------------------------------------------------------------------------

def _write(outs, folder):
    folder.mkdir(parents=True, exist_ok=True)
    for o in outs.values():
        o.df.to_csv(folder / o.file, index=False)


def _check_entries(entries, out):
    assert [e["name"] for e in entries] == [n for n, _ in PF.FIGURES]
    for e in entries:
        assert set(e) == {"name", "group", "title", "caption", "png", "pdf", "error"}
        assert e["group"] == common.GROUP and e["error"] is None and e["caption"]
        assert (out / e["png"]).stat().st_size > 0 and (out / e["pdf"]).stat().st_size > 0
        json.dumps(e)


def test_draw(tiny_joint, positivity_run, tmp_path):
    results, out = tmp_path / "results", tmp_path / "figures"
    _write(positivity_run["outs"], results)
    assert set(PF.INPUTS) <= {p.name for p in results.iterdir()}
    entries = PF.draw(results, tiny_joint["cfg"], out)
    _check_entries(entries, out)
    f = positivity_run["outs"][P.FOLDS_CSV].df
    oc = f[(f["signature"] == "osteoclast") & (f["depth"] == "native") & f["manuscript"]]["fold"].iloc[0]
    cap = entries[0]["caption"]
    assert f"osteoclast {oc:.2f}" in cap
    assert "replicate-mean counts" in cap and "block-bootstrap 95% intervals" in cap
    cfg = tiny_joint["cfg"]
    tumor, ref = common.tumor_and_reference(cfg)
    plotted = set(f["sample"])
    assert ref not in plotted
    assert ("tumor-bearing section" in cap) == (plotted <= set(tumor))
    assert PF.draw(tmp_path / "nothing", tiny_joint["cfg"], out) == []


def test_caption_names_non_reference_sections(tiny_joint):
    cfg = tiny_joint["cfg"]
    tumor, ref = common.tumor_and_reference(cfg)
    assert PF._who(cfg, ["not-a-tumor-section"]) == "the non-reference section"
    assert PF._who(cfg, ["a", "b"]) == "each non-reference section"
    if tumor:
        assert PF._who(cfg, tumor[:1]) == "the tumor-bearing section"


def test_draw_survives_non_estimable_sections(tiny_joint, positivity_run, tmp_path):
    f = positivity_run["outs"][P.FOLDS_CSV].df.copy()
    f.loc[f["signature"] == "osteoclast", ["fold", "fold_ci_lo", "fold_ci_hi"]] = np.nan
    f.loc[f["signature"] == "macrophage", ["fold", "fold_ci_lo", "fold_ci_hi"]] = [0.0, np.nan, np.nan]
    f.loc[f["signature"] == "osteoblast", ["fold_ci_lo", "fold_ci_hi"]] = np.nan
    ob = (f["signature"] == "osteoblast") & (f["depth"] == "native") & f["manuscript"]
    f.loc[ob, ["fold", "block_fold_ci_lo", "block_fold_ci_hi"]] = [2.0, 0.5, np.nan]     # unbounded above
    results, out = tmp_path / "results", tmp_path / "figures"
    results.mkdir()
    f.to_csv(results / P.FOLDS_CSV, index=False)
    entries = PF.draw(results, tiny_joint["cfg"], out)
    _check_entries(entries, out)
    assert "osteoclast not estimable" in entries[0]["caption"]
    assert "block bootstrap 0.50–unbounded" in entries[0]["caption"]
    # an older table without the block-bootstrap columns still draws, without them in the caption
    f.drop(columns=["block_fold_ci_lo", "block_fold_ci_hi"]).to_csv(results / P.FOLDS_CSV, index=False)
    entries = PF.draw(results, tiny_joint["cfg"], out)
    _check_entries(entries, out)
    assert "block-bootstrap" not in entries[0]["caption"]
