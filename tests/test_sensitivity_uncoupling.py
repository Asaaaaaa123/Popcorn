"""Figure 4i sensitivity (bone_sensitivity.uncoupling): null models on synthetic patterns with known answers, the
depth-matching draws and the manuscript reproduction on the tiny demo, and the figure."""
from __future__ import annotations

import copy

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.spatial import cKDTree

import pipeline_config as pc
from bone_sensitivity import common
from bone_sensitivity import uncoupling as U
from bone_sensitivity import uncoupling_figures as UF

BIN = 16.0
OUTPUT_FILES = ["sensitivity_Figure4i_null_models.csv", "sensitivity_Figure4i_summary.csv",
                "sensitivity_Figure4i_decomposition.csv", "sensitivity_Figure4i_call_stability.csv"]
CALL_COLUMNS = {"pooled_dm": ("pooled_dm_pos_", "program_dm"), "pooled": ("pooled_pos_", "program"),
                "native": ("pos_", "program")}


def _grid(n=60):
    """Every bin of an n x n array grid, in um."""
    r, c = np.divmod(np.arange(n * n), n)
    return np.column_stack([c, r]).astype(float) * BIN


def _clusters(xy, n_clusters, radius_bins, rng):
    """Union of discs of radius_bins bins around random bin centres."""
    centres = xy[rng.choice(len(xy), n_clusters, replace=False)]
    d = np.min(np.linalg.norm(xy[:, None, :] - centres[None, :, :], axis=2), axis=1)
    return d <= radius_bins * BIN + 1e-9


# ---- statistical core ----------------------------------------------------------------------------------------------

def test_independent_clumped_patterns_inflate_csr_but_not_fixed_b():
    """Two independent clumped patterns: CSR of both calls says 'farther than chance' (ratio > 1) only because the
    clumped B pattern leaves most of the ROI far from any B bin; keeping B fixed gives ratio ~ 1 for every pattern
    pair."""
    xy = _grid()
    csr, fixed_b, fixed_a = [], [], []
    for k in range(4):
        rng = np.random.default_rng(100 + k)
        A, B = _clusters(xy, 40, 1, rng), _clusters(xy, 12, 2, rng)
        obs = U.nearest_distance(xy, A, B).mean()
        csr.append(U.null_statistics(obs, U.null_means("csr_both", xy, A, B, 99, np.random.default_rng(0)))["ratio"])
        fixed_b.append(U.null_statistics(obs, U.null_means("fixed_b", xy, A, B, 99, np.random.default_rng(0)))["ratio"])
        fixed_a.append(U.null_statistics(obs, U.null_means("fixed_a", xy, A, B, 99, np.random.default_rng(0)))["ratio"])
    assert min(csr) > 1.5
    assert min(fixed_a) > 1.5                                   # fixed_a still scatters the clumped B calls
    assert 0.85 < np.mean(fixed_b) < 1.15 and all(0.6 < r < 1.4 for r in fixed_b)


def test_attraction_is_detected_by_the_nulls_that_keep_b():
    """A placed right next to B (a coupled pattern): ratio < 1 with B kept and with the habitat null, and the
    one-sided P value for 'closer than chance' is at its floor."""
    xy = _grid(40)
    rng = np.random.default_rng(1)
    B = _clusters(xy, 8, 1, rng)
    near = U.nearest_distance(xy, np.ones(len(xy), bool), B)
    A = (near > 0) & (near <= 1.5 * BIN)
    obs = U.nearest_distance(xy, A, B).mean()
    for model in ("fixed_b", "habitat_fixed_b"):
        st = U.null_statistics(obs, U.null_means(model, xy, A, B, 49, np.random.default_rng(0),
                                                 habitat=_clusters(xy, 20, 3, rng)))
        assert st["ratio"] < 0.5, model
        assert st["p_less"] == pytest.approx(1 / 50) and st["p_greater"] == pytest.approx(1), model


def test_random_labels_detect_segregation_among_the_called_bins():
    """A and B in separate halves of the grid: permuting the labels over the called positions mixes them, so the
    observed distance is far above the null; with every called bin co-positive the null equals the observation."""
    xy = _grid(40)
    rng = np.random.default_rng(5)
    blobs = _clusters(xy, 16, 1, rng)
    left = xy[:, 0] < 20 * BIN
    A, B = blobs & left, blobs & ~left
    obs = U.nearest_distance(xy, A, B).mean()
    st = U.null_statistics(obs, U.null_means("random_labels", xy, A, B, 49, np.random.default_rng(0)))
    assert st["ratio"] > 2 and st["p_greater"] == pytest.approx(1 / 50)
    same = U.null_means("random_labels", xy, blobs, blobs, 5, np.random.default_rng(0))
    np.testing.assert_array_equal(same, 0.0)


def test_habitat_covering_the_roi_equals_fixed_b():
    xy = _grid(30)
    rng = np.random.default_rng(2)
    A, B = _clusters(xy, 10, 1, rng), _clusters(xy, 5, 1, rng)
    hab = U.null_means("habitat_fixed_b", xy, A, B, 25, np.random.default_rng(7), habitat=np.ones(len(xy), bool))
    np.testing.assert_array_equal(hab, U.null_means("fixed_b", xy, A, B, 25, np.random.default_rng(7)))
    small = np.zeros(len(xy), bool)                             # habitat = the called bins only: A stays inside them
    pool = U.null_means("habitat_fixed_b", xy, A, B, 25, np.random.default_rng(7), habitat=small)
    to_b = U.nearest_distance(xy, np.ones(len(xy), bool), B)
    assert pool.max() <= to_b[A | B].max() + 1e-9


def test_random_labels_keep_counts_positions_and_co_positives():
    rng = np.random.default_rng(3)
    A, B = rng.random(500) < 0.2, rng.random(500) < 0.15
    for seed in range(5):
        A2, B2 = U.random_labels_draw(A, B, np.random.default_rng(seed))
        assert A2.sum() == A.sum() and B2.sum() == B.sum()
        np.testing.assert_array_equal(A2 | B2, A | B)
        np.testing.assert_array_equal(A2 & B2, A & B)
    assert not np.array_equal(U.random_labels_draw(A, B, np.random.default_rng(0))[0], A)


def test_null_statistics_known_values():
    nulls = np.arange(1, 100, dtype=float)                      # mean 50
    st = U.null_statistics(50.0, nulls)
    assert st["ratio"] == pytest.approx(1.0)
    assert st["ratio_lo"] == pytest.approx(50 / np.percentile(nulls, 97.5))
    assert st["ratio_hi"] == pytest.approx(50 / np.percentile(nulls, 2.5))
    assert st["p_greater"] == pytest.approx(51 / 100) and st["p_less"] == pytest.approx(51 / 100)
    assert U.null_statistics(200.0, nulls)["p_greater"] == pytest.approx(1 / 100)
    assert U.null_statistics(0.0, nulls)["p_less"] == pytest.approx(1 / 100)


def test_decomposition_product_is_the_csr_ratio_and_clumping_is_measured():
    xy = _grid(40)
    rng = np.random.default_rng(4)
    A, B = _clusters(xy, 15, 1, rng), _clusters(xy, 6, 2, rng)
    csr = U.null_means("csr_both", xy, A, B, 49, np.random.default_rng(0))
    row = U.decomposition(xy, A, B, csr.mean(), BIN, np.random.default_rng(0))
    obs = U.nearest_distance(xy, A, B).mean()
    assert row["factor_a_specific"] * row["factor_b_pattern"] == pytest.approx(obs / csr.mean(), rel=1e-12)
    assert row["baseline_all_roi_to_b_um"] == pytest.approx(U.nearest_distance(xy, np.ones(len(xy), bool), B).mean())
    assert row["frac_a_clumped"] == 1.0 and row["frac_b_clumped"] == 1.0
    assert row["nn_self_ratio_a"] < 1 and row["nn_self_ratio_b"] < 1
    sparse = np.zeros(len(xy), bool)                            # isolated bins three bins apart: dispersed, unclumped
    sparse[(xy[:, 0] % (3 * BIN) == 0) & (xy[:, 1] % (3 * BIN) == 0)] = True
    assert U.fraction_clumped(xy, sparse, 1.5 * BIN) == 0.0
    assert U.nn_self_ratio(xy, sparse, np.random.default_rng(0)) > 1
    assert np.isnan(U.fraction_clumped(xy, np.eye(1, len(xy), 0, dtype=bool)[0], 1.5 * BIN))


def test_confined_unclumped_b_still_inflates_factor_b_pattern():
    """factor_b_pattern is the departure of B from random placement at any scale: isolated, dispersed B bins confined
    to half of the ROI (no local clumping) still give a factor well above 1, while A (random) gives factor_a ~ 1."""
    xy = _grid(60)
    B = (xy[:, 0] < 30 * BIN) & (xy[:, 0] % (3 * BIN) == 0) & (xy[:, 1] % (3 * BIN) == 0)
    A = np.random.default_rng(8).random(len(xy)) < 0.05
    csr = U.null_means("csr_both", xy, A, B, 49, np.random.default_rng(0))
    row = U.decomposition(xy, A, B, csr.mean(), BIN, np.random.default_rng(0))
    assert row["frac_b_clumped"] == 0.0 and row["nn_self_ratio_b"] > 1
    assert row["factor_b_pattern"] > 3 and 0.8 < row["factor_a_specific"] < 1.2


def test_call_stability_known_counts():
    """3 draws over 6 bins of one section (bins 0-3 in the ROI). In the ROI bin 0 is always called and bin 1 only in
    draw 0, so draws 1 and 2 each differ from draw 0 in 1 of the 2 bins called in either; over the section they differ
    in 3 of 4 and 2 of 3."""
    D = np.array([[1, 1, 0, 0, 0, 0],
                  [1, 0, 0, 0, 1, 1],
                  [1, 0, 0, 0, 0, 1]], bool)
    smp = np.array(["S"] * 6)
    roi = np.array([True, True, True, True, False, False])
    t = U.call_stability({"x": D}, smp, roi, ["S"]).set_index("region")
    assert list(U.call_stability({"x": D}, smp, roi, ["S"]).columns) == U.STABILITY_COLUMNS
    r = t.loc["roi"]
    assert (r["n_draws"], r["n_bins"], r["n_called_draw0"], r["n_called_min"], r["n_called_max"]) == (3, 4, 2, 1, 2)
    assert (r["n_ever_called"], r["n_always_called"], r["n_unstable"]) == (2, 1, 1)
    assert r["pct_unstable_of_ever_called"] == pytest.approx(50) and r["pct_unstable_of_bins"] == pytest.approx(25)
    assert r["pct_changed_vs_draw0_mean"] == pytest.approx(50) and r["pct_changed_vs_draw0_max"] == pytest.approx(50)
    s = t.loc["section"]
    assert (s["n_bins"], s["n_ever_called"], s["n_always_called"], s["n_unstable"]) == (6, 4, 1, 3)
    assert s["pct_changed_vs_draw0_mean"] == pytest.approx((100 * 3 / 4 + 100 * 2 / 3) / 2)
    assert s["pct_changed_vs_draw0_max"] == pytest.approx(75)
    one = U.call_stability({"x": D[:1]}, smp, roi, ["S"])
    assert (one["n_unstable"] == 0).all() and one["pct_changed_vs_draw0_mean"].isna().all()


def _light_grid(n=30):
    """Two sections (C, T) of an n x n grid, every bin in the ROI, as the light object of _csr_with_p_less."""
    r, c = np.divmod(np.arange(n * n), n)
    obs = pd.DataFrame({"sample": np.repeat(["C", "T"], n * n), "array_row": np.tile(r, 2),
                        "array_col": np.tile(c, 2), "in_marrow_roi": True})
    obs.index = obs.index.astype(str)
    return ad.AnnData(obs=obs)


def test_csr_p_less_from_the_replayed_draws():
    """A right next to B in both sections: the observed distance is below every csr_both null mean, so p_less is at
    its floor 1 / (n_perm + 1) and p_greater is 1; A far from B (other half): p_less is 1."""
    light = _light_grid()
    xy = common.array_xy_um(light, BIN)
    smp, roi = light.obs["sample"].to_numpy(), np.ones(light.n_obs, bool)
    unc = {"n_perm": 39, "min_pos": 5}
    B = np.random.default_rng(3).random(light.n_obs) < 0.04
    near = np.concatenate([U.nearest_distance(xy[smp == s], np.ones((smp == s).sum(), bool), B[smp == s])
                           for s in ("C", "T")])
    A = (near > 0) & (near <= 1.5 * BIN)
    summ, nulls = U._csr_with_p_less(light, A, B, xy, roi, smp, ["C", "T"], unc, "C", 0, BIN)
    for s in ("C", "T"):
        assert nulls[s] is not None and len(nulls[s]) == 39
        r = summ.set_index("sample").loc[s]
        assert nulls[s].mean() == pytest.approx(r["null_mean_um"], rel=1e-12)
        assert (np.sum(nulls[s] <= r["mean_nn_um"]) + 1) / 40 == pytest.approx(1 / 40)
        assert r["p_greater_than_null_mean"] == pytest.approx(1)
    left = xy[:, 0] < 15 * BIN
    summ, nulls = U._csr_with_p_less(light, A & left, B & ~left, xy, roi, smp, ["C", "T"], unc, "C", 0, BIN)
    for s in ("C", "T"):
        obs = summ.set_index("sample").loc[s, "mean_nn_um"]
        assert (np.sum(nulls[s] <= obs) + 1) / 40 == pytest.approx(1)


def test_ratio_of_ratios_and_summary():
    rows = []
    for d, (rc, rt) in enumerate([(2.0, 3.0), (2.0, 4.0)]):
        rows += [{"call_set": "pooled_dm", "draw": d, "null_model": "csr_both", "sample": "C", "status": "estimable",
                  "ratio": rc, "p_greater": 0.01, "p_less": 1.0},
                 {"call_set": "pooled_dm", "draw": d, "null_model": "csr_both", "sample": "T", "status": "estimable",
                  "ratio": rt, "p_greater": 0.02, "p_less": 0.9}]
    rows.append({"call_set": "native", "draw": 0, "null_model": "csr_both", "sample": "C",
                 "status": "not estimable (n_positive < 10)", "ratio": np.nan})
    rows.append({"call_set": "native", "draw": 0, "null_model": "csr_both", "sample": "T", "status": "estimable",
                 "ratio": 1.5, "p_greater": 0.5, "p_less": 0.5})
    df = U._ratio_of_ratios(pd.DataFrame(rows).reindex(columns=U.NULL_COLUMNS), "C")
    np.testing.assert_allclose(df["ratio_of_ratios_vs_reference"].to_numpy()[:4], [1, 1.5, 1, 2])
    assert df["ratio_of_ratios_vs_reference"].iloc[4:].isna().all()     # reference not estimable
    s = U.summarize(df, "pooled_dm").set_index(["call_set", "sample"])
    assert list(U.summarize(df, "pooled_dm").columns) == U.SUMMARY_COLUMNS
    t = s.loc[("pooled_dm", "T")]
    assert t["n_draws"] == 2 and t["ratio_mean"] == 3.5 and t["ratio_min"] == 3 and t["ratio_max"] == 4
    assert t["ratio_of_ratios_mean"] == 1.75 and t["p_greater_max"] == 0.02 and bool(t["manuscript"])
    assert s.loc[("native", "C")]["n_draws"] == 0 and not bool(s.loc[("native", "T")]["manuscript"])


# ---- tiny demo -----------------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def tiny_run(tiny_joint, tmp_path_factory):
    """run() once on a copy of the tiny joint object (2 draws, 19 permutations); CSVs written to a folder."""
    a = tiny_joint["a"].copy()
    cfg = copy.deepcopy(tiny_joint["cfg"])
    cfg["sensitivity_analyses"]["uncoupling"].update(n_draws=2, n_perm=19)
    before = {"obs": a.obs.copy(), "layers": sorted(a.layers.keys()), "uns": sorted(a.uns.keys())}
    outputs = U.run(a, cfg)
    res = tmp_path_factory.mktemp("uncoupling_results")
    for o in outputs:
        o.df.to_csv(res / o.file, index=False)
    return {"a": a, "cfg": cfg, "outputs": {o.file: o for o in outputs}, "before": before, "results": res,
            "out": tiny_joint["out"]}


def test_draw_zero_reproduces_the_pipeline_calls(tiny_joint):
    a = tiny_joint["a"]
    cfg = tiny_joint["cfg"]
    draws = U.pooled_dm_call_draws(a, cfg, 2)
    for n in (cfg["uncoupling"]["a"], cfg["uncoupling"]["b"]):
        assert draws[n].shape == (2, a.n_obs)
        np.testing.assert_array_equal(draws[n][0], common.as_bool(a.obs[f"pooled_dm_pos_{n}"]))
        assert not np.array_equal(draws[n][0], draws[n][1])      # another downsampling draw changes some calls


def _tiny_calls(tiny_run):
    """{(call_set, draw): (A, B, assignment column)} over all bins of the tiny demo, read from obs; pooled_dm draw 1
    is recovered from the 2-draw call fractions (2 x fraction - draw-0 call)."""
    a, unc = tiny_run["a"], tiny_run["cfg"]["uncoupling"]
    out = {}
    for cs, (prefix, assign) in CALL_COLUMNS.items():
        out[(cs, 0)] = (common.as_bool(a.obs[prefix + unc["a"]]), common.as_bool(a.obs[prefix + unc["b"]]), assign)
    d1 = []
    for n, d0 in zip((unc["a"], unc["b"]), out[("pooled_dm", 0)][:2]):
        v = np.round(2 * a.obs[f"pooled_dm_pos_frac_{n}"].to_numpy(float)).astype(int) - d0.astype(int)
        assert set(np.unique(v)) <= {0, 1}
        d1.append(v.astype(bool))
    out[("pooled_dm", 1)] = (d1[0], d1[1], "program_dm")
    return out


def test_run_outputs_and_obs_columns(tiny_run):
    outs, a, before = tiny_run["outputs"], tiny_run["a"], tiny_run["before"]
    assert list(outs) == OUTPUT_FILES
    nm = outs["sensitivity_Figure4i_null_models.csv"].df
    assert list(nm.columns) == U.NULL_COLUMNS
    assert list(outs["sensitivity_Figure4i_summary.csv"].df.columns) == U.SUMMARY_COLUMNS
    assert list(outs["sensitivity_Figure4i_decomposition.csv"].df.columns) == U.DECOMPOSITION_COLUMNS
    assert list(outs["sensitivity_Figure4i_call_stability.csv"].df.columns) == U.STABILITY_COLUMNS
    assert set(nm["null_model"]) == set(U.NULL_MODELS) and set(nm["call_set"]) == {"pooled_dm", "pooled", "native"}
    assert sorted(nm.loc[nm["call_set"] == "pooled_dm", "draw"].unique()) == [0, 1]
    assert (nm["status"] == "estimable").all()
    hab = nm[nm["null_model"] == "habitat_fixed_b"]
    assert nm.loc[nm.index.difference(hab.index), ["n_habitat", "n_habitat_extra"]].isna().all().all()
    unc = tiny_run["cfg"]["uncoupling"]
    added = [c for c in a.obs.columns if c not in before["obs"].columns]
    assert added == [f"pooled_dm_pos_frac_{unc['a']}", f"pooled_dm_pos_frac_{unc['b']}"]
    for c in added:
        assert a.obs[c].between(0, 1).all() and set(np.round(a.obs[c] * 2)) <= {0, 1, 2}
    for c in before["obs"].columns:
        pd.testing.assert_series_equal(a.obs[c], before["obs"][c])
    assert sorted(a.layers.keys()) == before["layers"] and sorted(a.uns.keys()) == before["uns"]


def test_habitat_pool_and_csr_p_less_match_an_independent_recomputation(tiny_run):
    """n_habitat / n_habitat_extra from the assignment column of each call set (program_dm for pooled_dm, program
    otherwise) and the habitat classes; csr_both null mean, p_greater and p_less from a replay, written here, of
    bone_analysis.uncoupling_metric's generator stream (one default_rng(seed) over the estimable sections)."""
    a, cfg = tiny_run["a"], tiny_run["cfg"]
    unc, st = cfg["uncoupling"], common.settings(cfg, "uncoupling")
    nm = tiny_run["outputs"]["sensitivity_Figure4i_null_models.csv"].df.set_index(
        ["call_set", "draw", "null_model", "sample"])
    hab_sigs = [k for k, d in cfg["signatures"].items() if d.get("class") in set(st["habitat_classes"])]
    assert hab_sigs
    smp, roi = a.obs["sample"].astype(str).to_numpy(), common.as_bool(a.obs["in_marrow_roi"])
    xy_all = np.column_stack([a.obs["array_col"], a.obs["array_row"]]).astype(float) * pc.bin_um(cfg["bin_size"])
    fig4i = pd.read_csv(tiny_run["out"] / "Figure4i_observed_to_null_summary.csv").set_index("sample")
    n_perm = int(unc["n_perm"])
    for (cs, d), (A_all, B_all, assign) in _tiny_calls(tiny_run).items():
        H_all = a.obs[assign].astype(str).isin(hab_sigs).to_numpy()
        rng = np.random.default_rng(int(cfg["seed"]))
        for s in common.sample_ids(a, cfg):
            m = (smp == s) & roi
            xy, A, B, H = xy_all[m], A_all[m], B_all[m], H_all[m]
            h = nm.loc[(cs, d, "habitat_fixed_b", s)]
            assert (h["n_a"], h["n_b"]) == (A.sum(), B.sum())
            assert h["n_habitat"] == (H | A | B).sum() and h["n_habitat_extra"] == (H & ~A & ~B).sum() > 0
            if A.sum() < unc["min_pos"] or B.sum() < unc["min_pos"]:
                continue
            obs = cKDTree(xy[B]).query(xy[A])[0].mean()
            nulls = []
            for _ in range(n_perm):
                ia = rng.choice(len(xy), int(A.sum()), replace=False)
                ib = rng.choice(len(xy), int(B.sum()), replace=False)
                nulls.append(cKDTree(xy[ib]).query(xy[ia])[0].mean())
            nulls = np.asarray(nulls)
            c = nm.loc[(cs, d, "csr_both", s)]
            assert c["null_mean_um"] == pytest.approx(nulls.mean(), rel=1e-12)
            assert c["p_greater"] == pytest.approx((np.sum(nulls >= obs) + 1) / (n_perm + 1))
            assert c["p_less"] == pytest.approx((np.sum(nulls <= obs) + 1) / (n_perm + 1))
            if (cs, d) == ("pooled_dm", 0):                          # the replay itself reproduces Figure 4i
                assert nulls.mean() == pytest.approx(fig4i.loc[s, "null_mean_um"], rel=1e-12)


def test_call_stability_on_the_tiny_demo(tiny_run):
    """The stability table agrees with the per-bin draw fractions and the draw-0 (pipeline) calls."""
    a, cfg = tiny_run["a"], tiny_run["cfg"]
    unc = cfg["uncoupling"]
    t = tiny_run["outputs"]["sensitivity_Figure4i_call_stability.csv"].df
    samples = common.sample_ids(a, cfg)
    assert len(t) == 2 * len(samples) * 2 and set(t["call_set"]) == {"pooled_dm"} and (t["n_draws"] == 2).all()
    smp, roi = a.obs["sample"].astype(str).to_numpy(), common.as_bool(a.obs["in_marrow_roi"])
    t = t.set_index(["signature", "sample", "region"])
    for n in (unc["a"], unc["b"]):
        frac = a.obs[f"pooled_dm_pos_frac_{n}"].to_numpy(float)
        d0 = common.as_bool(a.obs[f"pooled_dm_pos_{n}"])
        for s in samples:
            for region, m in (("roi", (smp == s) & roi), ("section", smp == s)):
                r = t.loc[(n, s, region)]
                assert r["n_bins"] == m.sum() and r["n_called_draw0"] == (d0 & m).sum()
                assert r["n_unstable"] == ((frac == 0.5) & m).sum() and r["n_always_called"] == ((frac == 1) & m).sum()
                assert r["n_ever_called"] == ((frac > 0) & m).sum()
                assert r["pct_changed_vs_draw0_mean"] == pytest.approx(100 * r["n_unstable"] / r["n_ever_called"])
    assert t["n_unstable"].sum() > 0


def test_habitat_without_bins_beyond_the_calls_is_not_estimable(tiny_joint):
    """With no habitat-class assignment the habitat pool is the called bins only: habitat_fixed_b is not estimable
    (n_habitat_extra 0), the other nulls are; without the pooled_dm call set the stability table is empty."""
    a = tiny_joint["a"].copy()
    cfg = copy.deepcopy(tiny_joint["cfg"])
    cfg["uncoupling"]["calls"] = "native"
    cfg["sensitivity_analyses"]["uncoupling"].update(call_sets=["native"], null_models=["fixed_b", "habitat_fixed_b"],
                                                     habitat_classes=["no such class"], n_perm=9)
    outs = U.run(a, cfg)
    assert [o.file for o in outs] == OUTPUT_FILES
    nm = outs[0].df
    hab, fb = nm[nm["null_model"] == "habitat_fixed_b"], nm[nm["null_model"] == "fixed_b"]
    assert (hab["status"] == "not estimable (no habitat bins beyond the calls)").all() and hab["ratio"].isna().all()
    assert (hab["n_habitat_extra"] == 0).all()
    smp, roi = a.obs["sample"].astype(str).to_numpy(), common.as_bool(a.obs["in_marrow_roi"])
    calls = np.logical_or(*(common.as_bool(a.obs[f"pos_{cfg['uncoupling'][k]}"]) for k in ("a", "b")))
    for _, r in hab.iterrows():
        assert r["n_habitat"] == (calls & (smp == r["sample"]) & roi).sum()
    assert (fb["status"] == "estimable").all()
    assert outs[3].df.empty and list(outs[3].df.columns) == U.STABILITY_COLUMNS


def test_csr_both_draw_zero_reproduces_figure_4i(tiny_run):
    nm = tiny_run["outputs"]["sensitivity_Figure4i_null_models.csv"].df
    ref = pd.read_csv(tiny_run["out"] / "Figure4i_observed_to_null_summary.csv").set_index("sample")
    m = nm[(nm["call_set"] == "pooled_dm") & (nm["draw"] == 0) & (nm["null_model"] == "csr_both")].set_index("sample")
    assert list(m.index) == list(ref.index)
    for mine, theirs in (("n_a", "n_a"), ("n_b", "n_b"), ("n_units", "n_units_roi"), ("obs_mean_um", "mean_nn_um"),
                         ("null_mean_um", "null_mean_um"), ("ratio", "ratio_obs_null_mean"),
                         ("ratio_lo", "ratio_obs_null_lo"), ("ratio_hi", "ratio_obs_null_hi"),
                         ("p_greater", "p_greater_than_null_mean"),
                         ("ratio_of_ratios_vs_reference", "ratio_of_ratios_vs_reference")):
        np.testing.assert_allclose(m[mine].astype(float), ref.loc[m.index, theirs].astype(float), rtol=1e-12,
                                   err_msg=mine)
    summ = tiny_run["outputs"]["sensitivity_Figure4i_summary.csv"].df
    man = summ[summ["manuscript"]]
    assert set(man["call_set"]) == {"pooled_dm"} and set(man["null_model"]) == {"csr_both"}
    assert (man["n_draws"] == 2).all()


def test_decomposition_on_the_tiny_demo(tiny_run):
    nm = tiny_run["outputs"]["sensitivity_Figure4i_null_models.csv"].df
    dec = tiny_run["outputs"]["sensitivity_Figure4i_decomposition.csv"].df
    csr = nm[nm["null_model"] == "csr_both"].set_index(["call_set", "draw", "sample"])["ratio"]
    d = dec.set_index(["call_set", "draw", "sample"])
    np.testing.assert_allclose(d["factor_a_specific"] * d["factor_b_pattern"], csr.loc[d.index], rtol=1e-12)
    assert d["frac_a_clumped"].between(0, 1).all() and d["nn_self_ratio_b"].gt(0).all()


def test_figure(tiny_run, tmp_path):
    entries = UF.draw(tiny_run["results"], tiny_run["cfg"], tmp_path)
    assert [e["name"] for e in entries] == [n for n, _ in UF.FIGURES]
    e = entries[0]
    assert (tmp_path / e["png"]).is_file() and (tmp_path / e["pdf"]).is_file()
    assert e["group"] == common.GROUP and e["error"] is None and e["title"] == UF.FIGURES[0][1]
    nm = tiny_run["outputs"]["sensitivity_Figure4i_null_models.csv"].df
    man = nm[(nm["call_set"] == "pooled_dm") & (nm["draw"] == 0) & (nm["null_model"] == "csr_both")]
    for v in man["ratio"]:
        assert f"{v:.2f}" in e["caption"]
    cap = e["caption"]
    assert "fixed_b" in cap and "independent" in cap
    assert "only the clumping" not in cap and "confinement to part of it" in cap
    assert "one depth-matching draw (seed" in cap                 # the manuscript calls are pooled_dm here
    assert "(19 random configurations)" in cap and "(19 draws)" not in cap
    assert "habitat bins outside the calls" in cap and "not called in all" in cap
    assert set(UF.INPUTS) == set(tiny_run["outputs"])


def test_not_estimable_sections_do_not_break_run_or_figure(tiny_joint, tmp_path):
    a = tiny_joint["a"].copy()
    cfg = copy.deepcopy(tiny_joint["cfg"])
    cfg["uncoupling"].update(min_pos=10 ** 6, calls="native")
    cfg["sensitivity_analyses"]["uncoupling"].update(call_sets=["native"], n_perm=9)
    outs = U.run(a, cfg)
    nm = outs[0].df
    assert nm["status"].str.startswith("not estimable").all() and nm["ratio"].isna().all()
    assert "pooled_dm_pos_frac_osteoblast" not in a.obs
    for o in outs:
        o.df.to_csv(tmp_path / o.file, index=False)
    entries = UF.draw(tmp_path, cfg, tmp_path / "figs")
    assert len(entries) == 1 and (tmp_path / "figs" / entries[0]["png"]).is_file()
    assert "manuscript setup was not estimable" in entries[0]["caption"]
    assert "depth-matching draw (seed" not in entries[0]["caption"]      # native calls: no downsampling draw
    cfg["uncoupling"]["calls"] = "pooled_dm"
    assert "(pooled_dm) was not computed" in UF.draw(tmp_path, cfg, tmp_path / "figs")[0]["caption"]
    assert UF.draw(tmp_path / "missing", cfg, tmp_path / "figs2") == []
