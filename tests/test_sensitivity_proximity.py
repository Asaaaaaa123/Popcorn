"""Proximity sensitivity (bone_sensitivity/proximity.py and proximity_figures.py): the baselines of the "within
50 um" tables, the depth-matched expectation, its spatial block-bootstrap interval, and the figure."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.ndimage import gaussian_filter

import bone_analysis as ba
from bone_sensitivity import common
from bone_sensitivity import proximity as prox
from bone_sensitivity import proximity_figures as proxfig


def _grid(n=60, step=16.0):
    return np.indices((n, n)).reshape(2, -1).T * step


def _scenario(seed, kind):
    """A 60 x 60 grid of 16 um bins: the left half is deep (Poisson 300 UMI) and holds ten times more targets than
    the shallow right half (Poisson 50). 'depth': source bins are drawn like targets, from depth alone; 'adjacent':
    source bins are drawn among the bins next to a target (4-neighbours); 'copositive': like 'depth', plus half of the
    target bins also called source."""
    rng = np.random.default_rng(seed)
    P = _grid()
    deep = P[:, 0] < P[:, 0].max() / 2
    depth = np.where(deep, rng.poisson(300, len(P)), rng.poisson(50, len(P))).astype(float)
    p_t = np.where(deep, .05, .005)
    target = rng.random(len(P)) < p_t
    if kind == "adjacent":
        source = ~target & (ba._nearest_target_distance(P, target) <= 20) & (rng.random(len(P)) < .3)
    else:
        source = ~target & (rng.random(len(P)) < np.where(deep, .08, .005) / (1 - p_t))
    if kind == "copositive":
        source |= target & (rng.random(len(P)) < .5)
    return P, source, target, depth


def _smooth_depth_null(seed, side=160):
    """Continuous depth with no adjacency (the review's check): log UMI = log 65 + a spatially smooth field + noise
    (log-normal, median ~70), source and target calls logistic in z-scored log depth (slope 2.2), so source bins are
    ~4x deeper than the ROI median and lie mostly in the top ROI decile, where the within-radius rate keeps rising."""
    rng = np.random.default_rng(seed)
    P = _grid(side)
    field = gaussian_filter(rng.normal(size=(side, side)), 6).ravel()
    logd = np.log(65) + 0.7 * field / field.std() + 0.6 * rng.normal(size=side * side)
    z = (logd - logd.mean()) / logd.std()
    target = rng.random(len(P)) < 1 / (1 + np.exp(6.0 - 2.2 * z))
    source = (rng.random(len(P)) < 1 / (1 + np.exp(6.3 - 2.2 * z))) & ~target
    return P, source, target, np.round(np.exp(logd))


# ---- statistical core --------------------------------------------------------------------------------------------

def test_depth_strata_are_quantiles_with_merged_ties():
    s = prox.depth_strata(np.arange(100), 10)
    assert np.array_equal(np.bincount(s), np.full(10, 10))
    assert np.array_equal(s[:10], np.zeros(10)) and np.array_equal(s[-10:], np.full(10, 9))
    assert np.all(np.diff(s) >= 0)
    assert np.array_equal(prox.depth_strata(np.full(7, 3.0), 10), np.zeros(7))
    tied = prox.depth_strata(np.r_[np.zeros(60), np.arange(1, 41)], 10)       # 60% zeros: one stratum, not six
    assert len(np.unique(tied[:60])) == 1 and tied.max() == 4
    assert prox.depth_strata([], 10).size == 0


def test_source_cuts_refine_the_roi_strata_and_spread_the_source_bins():
    depth = np.arange(100, dtype=float)
    source = depth >= 80                                                        # the deepest 20 bins
    roi, both = prox.depth_strata(depth, 10), prox.depth_strata(depth, 10, source)
    assert np.all(np.diff(both) >= 0)
    for k in np.unique(both):                                                   # every ROI stratum is a union
        assert len(np.unique(roi[both == k])) == 1
    assert len(np.unique(roi[source])) == 2 and len(np.unique(both[source])) >= 9
    assert np.array_equal(prox.depth_strata(depth, 10, np.zeros(100, bool)), roi)


def test_matched_reference_median_known_answer():
    # stratum 0: reference 1, 2, 3 and source 4 (weight 1/3 each); stratum 1: reference 10, 11 and source 12, 13;
    # stratum 2: source 20 only, handed to stratum 1, whose reference bins carry 3/2 each. Half of 4 is reached at 10.
    depth = np.array([1, 2, 3, 4, 10, 11, 12, 13, 20], float)
    strata = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2])
    source = np.array([0, 0, 0, 1, 0, 0, 1, 1, 1], bool)
    assert prox.matched_reference_median(depth, strata, source, np.zeros(9, bool)) == 10.0
    # a perfect match returns the source bins' median, exactly as np.median
    depth = np.array([5, 9, 5, 9, 5, 9], float)
    strata, source = np.array([0, 1, 0, 1, 0, 1]), np.array([0, 0, 0, 0, 1, 1], bool)
    assert prox.matched_reference_median(depth, strata, source, np.zeros(6, bool)) == np.median(depth[source]) == 7.0
    assert np.isnan(prox.matched_reference_median(depth, strata, np.zeros(6, bool), np.zeros(6, bool)))


def test_stratum_rates_borrow_the_nearest_stratum():
    rates = prox.stratum_rates([1, 0, 0, 3], [2, 0, 0, 4])
    assert np.allclose(rates, [.5, .5, .75, .75])
    assert np.allclose(prox.stratum_rates([1, 0, 3], [2, 0, 4]), [.5, .5, .75])     # tie: the shallower stratum
    assert np.isnan(prox.stratum_rates([0, 0], [0, 0])).all()


def test_depth_matched_known_answer():
    # stratum 0: reference within 1, 0 (rate 1/2); stratum 1: reference 1, 1, 1, 0 (rate 3/4); stratum 2: no
    # reference (borrows stratum 1). Source: two bins in stratum 0 (within 1, 1), one in 1 (0), one co-positive in 2.
    strata = np.array([0, 0, 1, 1, 1, 1, 0, 0, 1, 2, 2])
    within = np.array([1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1], bool)
    source = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0], bool)
    target = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1], bool)
    obs, exp, obs_x, exp_x = prox.depth_matched(within, strata, source, target)
    assert obs == pytest.approx(75.0) and exp == pytest.approx(100 * (.5 + .5 + .75 + .75) / 4)
    assert obs_x == pytest.approx(100 * 2 / 3) and exp_x == pytest.approx(100 * (.5 + .5 + .75) / 3)
    # a resample (indices with repeats) is weighted by its multiplicities
    idx = np.array([0, 1, 6, 6, 8, 2, 3])
    o, e, _, _ = prox.depth_matched(within, strata, source, target, idx)
    assert o == pytest.approx(100 * 2 / 3) and e == pytest.approx(100 * (.5 + .5 + 1.0) / 3)
    assert np.isnan(prox.depth_matched(within, strata, np.zeros(11, bool), target)).all()


def test_section_row_baselines_match_direct_counts():
    rng = np.random.default_rng(3)
    P = rng.uniform(0, 600, (400, 2))
    target, source = rng.random(400) < .05, rng.random(400) < .1
    depth = rng.poisson(80, 400).astype(float)
    r = prox.section_row(P, source, target, depth, radius_um=40.0, n_boot=50)
    d = np.sqrt(((P[:, None, :] - P[None, target, :]) ** 2).sum(-1)).min(axis=1)
    within = d <= 40.0
    assert r["n_roi"] == 400 and r["n_source"] == source.sum() and r["n_target"] == target.sum()
    assert r["n_copositive"] == (source & target).sum() > 0
    assert r["pct_source_within"] == pytest.approx(100 * within[source].mean())
    assert r["pct_all_within"] == pytest.approx(100 * within.mean())
    assert r["pct_nontarget_within"] == pytest.approx(100 * within[~target].mean())
    assert r["pct_neither_within"] == pytest.approx(100 * within[~source & ~target].mean())
    assert r["pct_source_within_excl_copositive"] == pytest.approx(100 * within[source & ~target].mean())
    # F11: the all-ROI baseline counts the target bins themselves (0 um), the non-target baseline does not
    assert r["pct_all_within"] > r["pct_nontarget_within"]
    assert r["max_possible_ratio_vs_all"] == pytest.approx(100 / r["pct_all_within"])
    assert r["diff_pp"] == pytest.approx(r["pct_source_within"] - r["expected_depth_matched_pct"])
    assert r["ratio_depth_matched"] == pytest.approx(r["pct_source_within"] / r["expected_depth_matched_pct"])
    assert r["median_umi_source"] == np.median(depth[source]) and r["median_umi_roi"] == np.median(depth)
    assert r["diff_ci_lo"] <= r["diff_pp"] <= r["diff_ci_hi"]
    strata = prox.depth_strata(depth, 10, source)
    assert r["n_depth_strata"] == len(np.unique(strata)) > 10
    assert r["expected_depth_matched_pct"] == pytest.approx(prox.depth_matched(within, strata, source, target)[1])
    assert r["median_umi_matched_reference"] == prox.matched_reference_median(depth, strata, source, target)
    roi = prox.depth_strata(depth, 10)
    assert r["expected_roi_strata_pct"] == pytest.approx(prox.depth_matched(within, roi, source, target)[1])
    assert r["diff_roi_strata_pp"] == pytest.approx(r["pct_source_within"] - r["expected_roi_strata_pct"])
    assert r["diff_roi_strata_ci_lo"] <= r["diff_roi_strata_pp"] <= r["diff_roi_strata_ci_hi"]
    assert r["exclude_self"] is False


def test_exclude_self_measures_targets_to_the_nearest_other_target():
    P = _grid(10)
    target = np.zeros(100, bool)
    target[[0, 99]] = True                                                      # far apart: no other target near
    source = np.zeros(100, bool)
    source[[0, 1, 10, 50, 60, 70]] = True                                       # bin 0 is co-positive
    keep = prox.section_row(P, source, target, np.ones(100), n_boot=0)
    drop = prox.section_row(P, source, target, np.ones(100), n_boot=0, exclude_self=True)
    assert keep["exclude_self"] is False and drop["exclude_self"] is True
    assert keep["pct_source_within"] == pytest.approx(50.0) and drop["pct_source_within"] == pytest.approx(100 * 2 / 6)
    assert keep["pct_all_within"] - drop["pct_all_within"] == pytest.approx(2.0)   # the two target bins


def test_not_estimable_sections_and_single_block():
    P = _grid(10)
    depth = np.arange(100, dtype=float)
    target = np.zeros(100, bool)
    target[[5, 55]] = True
    few = np.zeros(100, bool)
    few[[3, 4, 6, 7]] = True                                                    # 16 or 32 um from bin 5
    r = prox.section_row(P, few, target, depth)
    assert r["status"] == "not estimable" and r["n_source"] == 4
    assert r["pct_source_within"] == 100.0 and np.isnan(r["diff_pp"]) and np.isnan(r["diff_ci_lo"])
    many = np.zeros(100, bool)
    many[10:30] = True
    r = prox.section_row(P, many, np.zeros(100, bool), depth)
    assert r["status"] == "not estimable" and r["pct_source_within"] == 0.0
    assert np.isnan(r["expected_depth_matched_pct"])
    r = prox.section_row(P, many, target, depth, block_um=1e6)                  # one block: no interval
    assert r["status"] == "estimable" and np.isfinite(r["diff_pp"]) and np.isnan(r["diff_ci_lo"])


def test_depth_alone_gives_no_difference_while_the_manuscript_baseline_shows_enrichment():
    diffs = []
    for seed in range(5):
        r = prox.section_row(*_scenario(seed, "depth"), n_boot=200, seed=0)
        assert r["status"] == "estimable" and r["n_copositive"] == 0
        assert r["pct_source_within"] > 1.4 * r["pct_all_within"]                # the manuscript comparison
        assert r["median_umi_source"] > 1.5 * r["median_umi_roi"]                # source bins are the deep bins
        assert r["diff_ci_lo"] < 0 < r["diff_ci_hi"]
        diffs.append(r["diff_pp"])
    assert abs(np.mean(diffs)) < 3


def test_continuous_depth_alone_gives_no_difference_and_roi_deciles_alone_do():
    """The review's null: continuous, spatially smooth log-normal depth drives both calls, no adjacency. The match on
    ROI and source quantiles balances depth and its difference centres on 0; ROI deciles alone leave the source bins
    at the deep end of the top decile (matched reference median ~0.8x the source median) and a positive difference."""
    rows = [prox.section_row(*_smooth_depth_null(seed), n_boot=200, seed=0) for seed in range(16)]
    assert all(r["status"] == "estimable" for r in rows)
    assert min(r["pct_source_within"] / r["pct_all_within"] for r in rows) > 1.5      # the manuscript comparison
    assert min(r["median_umi_source"] / r["median_umi_roi"] for r in rows) > 3
    diff, diff_roi = (np.array([r[c] for r in rows]) for c in ("diff_pp", "diff_roi_strata_pp"))
    covers = sum(r["diff_ci_lo"] <= 0 <= r["diff_ci_hi"] for r in rows)
    assert abs(diff.mean()) < 2.5 and covers >= 12
    balance = np.array([r["median_umi_matched_reference"] / r["median_umi_source"] for r in rows])
    assert np.all(np.abs(balance - 1) < .05)
    assert diff_roi.mean() > 4 and diff_roi.mean() > diff.mean() + 3
    assert np.all(np.array([r["median_umi_matched_reference_roi_strata"] / r["median_umi_source"] for r in rows]) < .9)


def test_true_adjacency_still_shows_a_positive_difference():
    for seed in range(5):
        r = prox.section_row(*_scenario(seed, "adjacent"), n_boot=200, seed=0)
        assert r["pct_source_within"] == 100.0
        assert r["diff_pp"] > 20 and r["diff_ci_lo"] > 10
        assert r["ratio_depth_matched"] > 1.3


def test_copositive_bins_inflate_the_observed_percentage_only():
    r = prox.section_row(*_scenario(0, "copositive"), n_boot=200, seed=0)
    assert r["n_copositive"] > 30
    assert r["pct_source_within"] > r["pct_source_within_excl_copositive"]
    assert r["diff_pp"] > r["diff_excl_copositive_pp"] + 5 and r["diff_ci_lo"] > 0
    assert r["diff_excl_copositive_ci_lo"] < 0 < r["diff_excl_copositive_ci_hi"]


# ---- run() on the tiny demo analysis -----------------------------------------------------------------------------

@pytest.fixture(scope="module")
def prox_run(tiny_joint):
    a = tiny_joint["a"].copy()
    before = (list(a.obs.columns), sorted(a.layers.keys()), sorted(a.uns.keys()), a.obs.copy())
    outs = prox.run(a, tiny_joint["cfg"])
    return {"a": a, "before": before, "outs": outs}


def test_run_outputs_and_leaves_the_object_unchanged(prox_run, tiny_joint):
    a, (cols, layers, uns, obs) = prox_run["a"], prox_run["before"]
    assert list(a.obs.columns) == cols and sorted(a.layers.keys()) == layers and sorted(a.uns.keys()) == uns
    pd.testing.assert_frame_equal(a.obs, obs)
    (out,) = prox_run["outs"]
    assert out.file == "sensitivity_proximity.csv" and list(out.df.columns) == prox.COLUMNS
    df, cfg = out.df, tiny_joint["cfg"]
    samples = common.sample_ids(a, cfg)
    expected_sets = []
    for p in cfg["proximity"]["pairs"]:
        assignment = ["program"] + (["program_dm"] if "program_dm" in a.obs else [])
        cols_ = assignment if p["mode"] == "assignment" else ["pos"]
        expected_sets += [(f"{p['from']}_to_{p['to']}", c, s) for c in cols_ for s in samples]
    assert list(zip(df["pair"], df["assignment_column"], df["sample"])) == expected_sets
    assert set(df["status"]) <= {"estimable", "not estimable"}
    small = (df["n_source"] < prox.MIN_SOURCE) | (df["n_target"] == 0)
    assert (df.loc[small, "status"] == "not estimable").all() and (df.loc[~small, "status"] == "estimable").all()
    est = df[df["status"] == "estimable"]
    assert len(est) and est[["expected_depth_matched_pct", "diff_pp", "ratio_depth_matched"]].notna().all().all()
    assert np.allclose(est["diff_pp"], est["pct_source_within"] - est["expected_depth_matched_pct"])
    assert np.allclose(est["diff_roi_strata_pp"], est["pct_source_within"] - est["expected_roi_strata_pct"])
    assert est[["median_umi_matched_reference", "median_umi_matched_reference_roi_strata"]].notna().all().all()
    assert np.allclose(df["max_possible_ratio_vs_all"], 100 / df["pct_all_within"])
    assert df.loc[small, ["diff_pp", "diff_roi_strata_pp", "n_depth_strata"]].isna().all().all()
    assert (df["radius_um"] == cfg["proximity"]["radius_um"]).all()
    assert (df["exclude_self"] == bool(cfg["proximity"]["exclude_self"])).all()


def test_run_reproduces_the_manuscript_proximity_tables(prox_run, tiny_joint):
    """pct_source_within and pct_all_within equal bone_analysis.proximity_table (and the CSVs of the run) after its
    1-decimal rounding; the counts equal its counts."""
    a, cfg, df = prox_run["a"], tiny_joint["cfg"], prox_run["outs"][0].df
    samples, R = common.sample_ids(a, cfg), f"{cfg['proximity']['radius_um']:g}"
    for p in cfg["proximity"]["pairs"]:
        src, dst, mode = p["from"], p["to"], p["mode"]
        table = ba.proximity_table(a, src, dst, mode, cfg["proximity"]["radius_um"], cfg["proximity"]["exclude_self"],
                                   samples=samples)
        saved = pd.read_csv(tiny_joint["out"] / f"{src}_to_{dst}_proximity.csv")
        rows = df[(df["pair"] == f"{src}_to_{dst}") & df["assignment_column"].isin(["program", "pos"])]
        assert list(rows["sample"]) == list(table["section"]) == samples
        if mode == "assignment":
            pct_src, pct_all = f"percent_{src}_within_{R}um_of_{dst}", f"percent_all_bins_within_{R}um_of_{dst}"
            counts = {"n_source": f"{src}_assigned_bins", "n_target": f"{dst}_assigned_bins"}
        else:
            pct_src, pct_all = f"percent_{src}_within_{R}um", f"percent_all_bins_within_{R}um"
            counts = {"n_source": f"{src}_positive_bins", "n_target": f"{dst}_positive_bins",
                      "n_copositive": f"{src}_also_{dst}_positive_bins"}
        counts["n_roi"] = "region_of_interest_bins"
        for reference in (table, saved):
            assert np.array_equal(rows["pct_source_within"].round(1).to_numpy(), reference[pct_src].to_numpy())
            assert np.array_equal(rows["pct_all_within"].round(1).to_numpy(), reference[pct_all].to_numpy())
            for mine, theirs in counts.items():
                assert np.array_equal(rows[mine].to_numpy(), reference[theirs].to_numpy())


def test_run_is_deterministic(prox_run, tiny_joint):
    again = prox.run(tiny_joint["a"].copy(), tiny_joint["cfg"])[0].df
    pd.testing.assert_frame_equal(again, prox_run["outs"][0].df)


# ---- figure -----------------------------------------------------------------------------------------------------

def test_draw_from_the_run_tables(prox_run, tiny_joint, tmp_path):
    cfg = tiny_joint["cfg"]
    df = prox_run["outs"][0].df
    df.to_csv(tmp_path / prox.FILE, index=False)
    assert proxfig.INPUTS == [prox.FILE] and proxfig.FIGURES[0][0] == "Sensitivity_proximity_baselines"
    entries = proxfig.draw(tmp_path, cfg, tmp_path / "fig")
    assert len(entries) == 1
    e = entries[0]
    assert set(e) == {"name", "group", "title", "caption", "png", "pdf", "error"} and e["error"] is None
    assert e["name"] == "Sensitivity_proximity_baselines" and e["group"] == common.GROUP
    assert (tmp_path / "fig" / e["png"]).is_file() and (tmp_path / "fig" / e["pdf"]).is_file()
    r = df[df["status"] == "estimable"].iloc[0]
    for x in (r["pct_source_within"], r["pct_all_within"], r["pct_nontarget_within"], r["expected_depth_matched_pct"],
              r["expected_roi_strata_pct"]):
        assert f"{x:.1f}%" in e["caption"]
    assert f"{r['diff_pp']:+.1f} pp" in e["caption"] and f"{r['diff_roi_strata_pp']:+.1f} pp" in e["caption"]
    assert f"{r['median_umi_matched_reference']:.0f} matched reference" in e["caption"]


@pytest.mark.parametrize("exclude_self", [False, True])
def test_caption_follows_exclude_self(prox_run, tiny_joint, tmp_path, exclude_self):
    """The all-ROI and co-positive wording follows the table's exclude_self column, and the configuration when the
    table has no such column."""
    df = prox_run["outs"][0].df.assign(exclude_self=exclude_self)
    df.to_csv(tmp_path / prox.FILE, index=False)
    caption = proxfig.draw(tmp_path, tiny_joint["cfg"], tmp_path / "fig")[0]["caption"]
    assert ("at 0 um" in caption) is not exclude_self
    assert ("nearest other target bin" in caption) is exclude_self
    df.drop(columns="exclude_self").to_csv(tmp_path / prox.FILE, index=False)
    cfg = {**tiny_joint["cfg"], "proximity": {**tiny_joint["cfg"]["proximity"], "exclude_self": exclude_self}}
    caption = proxfig.draw(tmp_path, cfg, tmp_path / "fig")[0]["caption"]
    assert ("nearest other target bin" in caption) is exclude_self


def test_draw_survives_sections_that_are_not_estimable(tmp_path, tiny_joint):
    cfg = tiny_joint["cfg"]
    base = {**dict.fromkeys(prox.COLUMNS, np.nan), "pair": "a_to_b", "source": "a", "target": "b",
            "mode": "positivity", "assignment_column": "pos", "status": "not estimable", "radius_um": 50.0}
    sid = [s["id"] for s in cfg["samples"]]
    df = pd.DataFrame([{**base, "sample": sid[0], "n_roi": 10, "n_source": 2, "n_target": 0, "n_copositive": 0,
                        "pct_source_within": 0.0, "pct_all_within": 0.0},
                       {**base, "sample": sid[1], "n_roi": 12, "n_source": 0, "n_target": 3}], columns=prox.COLUMNS)
    df.to_csv(tmp_path / prox.FILE, index=False)
    entries = proxfig.draw(tmp_path, cfg, tmp_path / "fig")
    assert len(entries) == 1 and (tmp_path / "fig" / entries[0]["png"]).is_file()
    assert "not estimable" in entries[0]["caption"]
    assert proxfig.draw(tmp_path / "missing", cfg, tmp_path / "fig2") == []
