"""Regression tests for the defects fixed in bone_analysis.py and bone_pending_masks.py (CHANGES.md)."""
from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp

import bone_analysis as ba
import bone_pending_masks as bpm
import pipeline_config as pc

GENES = [f"g{i}" for i in range(120)]


def _toy(n=400, seed=0, samples=("S1", "S2"), depth=(2.0, 4.0)):
    """Counts with two sections of different depth, normalized like the pipeline (layers counts and lognorm)."""
    rng = np.random.default_rng(seed)
    per = n // len(samples)
    X = np.vstack([rng.poisson(d, size=(per, len(GENES))) for d in depth[:len(samples)]]).astype(np.float32)
    a = ad.AnnData(sp.csr_matrix(X), obs=pd.DataFrame({"sample": np.repeat(samples, per)},
                                                      index=[f"b{i}" for i in range(per * len(samples))]),
                   var=pd.DataFrame(index=GENES))
    a.obs["sample"] = pd.Categorical(a.obs["sample"], categories=list(samples))
    a.obs["condition"] = np.where(a.obs["sample"] == samples[0], "control", "tumor")
    a.obs["total_counts"] = np.asarray(a.X.sum(axis=1)).ravel()
    side = int(np.ceil(np.sqrt(per)))
    grid = np.array([(i % side, i // side) for i in range(per)], dtype=float)
    a.obsm["spatial"] = np.vstack([grid * 27.0 for _ in samples])
    a.obs["array_col"] = np.tile(grid[:, 0].astype(int), len(samples))
    a.obs["array_row"] = np.tile(grid[:, 1].astype(int), len(samples))
    a.uns["microns_per_pixel"] = {s: 0.59 for s in samples}
    return ba.normalize_and_cluster(a, cluster=False)


SIGS = {"alpha": {"core": GENES[:4], "extended": GENES[:6], "competes": ["beta"], "excluded": {}},
        "beta": {"core": GENES[10:13], "extended": GENES[10:15], "competes": [], "excluded": {}},
        "ghost": {"core": ["absent1", "absent2"], "extended": ["absent1", "absent2"], "competes": [], "excluded": {}}}


# ---------------------------------------------------------------------------------------------- input
def test_norm_bc_keeps_visium_hd_coordinates():
    assert ba._norm_bc("s_016um_00012_00034-1") == "s_016um_00012_00034"
    assert ba._norm_bc("s_016um_00012_00034") == "s_016um_00012_00034"      # used to become "s_016um"
    assert ba._norm_bc(" AAACGT-2 ") == "AAACGT"


def test_microns_per_pixel_reads_scalefactors_not_sample_id():
    a = _toy()
    assert ba.microns_per_pixel(a, "S1") == pytest.approx(0.59)
    a.uns["microns_per_pixel"] = {}
    a.uns["spatial"] = {"S1": {"scalefactors": {"microns_per_pixel": 0.5}}}
    assert ba.microns_per_pixel(a, "S1") == 0.5
    assert ba.microns_per_pixel(a, "LL-2-1") == ba.UM_PER_PX["LL-2-1"]
    with pytest.raises(ValueError):
        ba.microns_per_pixel(a, "S2")


# ---------------------------------------------------------------------------------------------- scoring
def test_score_all_reads_the_prepared_matrix_not_a_stale_raw():
    a = _toy()
    fresh = a.copy()
    ba.score_all(fresh, SIGS)
    stale = a.copy()
    rng = np.random.default_rng(5)
    stale.raw = ad.AnnData(sp.csr_matrix(rng.random((a.n_obs, a.n_vars)).astype(np.float32)), var=a.var.copy())
    ba.score_all(stale, SIGS)
    np.testing.assert_array_equal(stale.obs["score_alpha"], fresh.obs["score_alpha"])


def test_null_threshold_follows_the_matrix():
    a = _toy()
    t1 = ba.null_threshold(a, "alpha", SIGS, n_perm=5, n_sub=300)
    b = a.copy()
    b.layers["lognorm"] = b.layers["lognorm"] * 3
    t2 = ba.null_threshold(b, "alpha", SIGS, n_perm=5, n_sub=300)
    assert t2 != pytest.approx(t1)
    assert np.isnan(ba.null_threshold(a, "ghost", SIGS))


def test_competitor_that_cannot_be_scored_does_not_veto():
    a = ba.score_all(_toy(), SIGS)
    a.obs["score_beta"] = np.nan
    a = ba.call_positive(a, ["alpha"], SIGS, min_core=0, n_perm=3, n_sub=300)
    assert a.obs["pos_alpha"].any()


def test_competitor_vetoes_only_when_strictly_higher():
    a = ba.score_all(_toy(), SIGS)
    a.obs["score_beta"] = a.obs["score_alpha"]                                   # ties keep the call
    with_tie = ba.call_positive(a.copy(), ["alpha"], SIGS, min_core=0, n_perm=3, n_sub=300).obs["pos_alpha"]
    a.obs["score_beta"] = a.obs["score_alpha"] + 1
    vetoed = ba.call_positive(a.copy(), ["alpha"], SIGS, min_core=0, n_perm=3, n_sub=300).obs["pos_alpha"]
    assert with_tie.any() and not vetoed.any()


def test_unscorable_signature_gets_a_false_column_and_overwrites_stale_calls():
    a = ba.score_all(_toy(), SIGS)
    a.obs["pos_ghost"] = True                                                    # stale call from a previous run
    progress = pc.Progress()
    a = ba.call_positive(a, ["ghost"], SIGS, n_perm=3, n_sub=300, progress=progress)
    assert a.obs["pos_ghost"].dtype == bool and not a.obs["pos_ghost"].any()
    assert np.isnan(a.uns["positivity_thresholds"]["ghost"])
    assert any("ghost" in w for w in progress.warnings)


def test_positivity_table_works_for_any_sample_ids_and_empty_sections():
    a = ba.call_positive(ba.score_all(_toy(samples=("A", "B")), SIGS), ["alpha"], SIGS, min_core=0, n_perm=3,
                         n_sub=300)
    a.obs["roi"] = (a.obs["sample"] == "A").to_numpy()
    t = ba.positivity_table(a, ["alpha"], roi_key="roi", reference="A", samples=["A", "B"])
    row_b = t[t["sample"] == "B"].iloc[0]
    assert row_b["n_bins"] == 0 and np.isnan(row_b["pct_positive"])
    assert "fold_tumor_over_control" in t and "fold_vs_reference" in t
    three = ba.positivity_table(a, ["alpha"], reference="A", samples=["A", "B", "C"])
    assert "fold_tumor_over_control" not in three


def test_fold_tumor_over_control_follows_the_conditions_and_the_exact_counts():
    a = ba.call_positive(ba.score_all(_toy(samples=("C", "T")), SIGS), ["alpha"], SIGS, min_core=0, n_perm=3,
                         n_sub=300)
    t = ba.positivity_table(a, ["alpha"], reference="T", samples=["C", "T"], tumor_samples=["T"])
    w = t.set_index("sample")
    rate = w["n_positive"] / w["n_bins"]                                         # exact, not the rounded pct_positive
    expected = round(rate["T"] / rate["C"], 2) if rate["C"] else np.nan
    assert t["fold_tumor_over_control"].iloc[0] == pytest.approx(expected, nan_ok=True)   # tumor / control, not ref
    expected_ref = round(rate["C"] / rate["T"], 4) if rate["T"] else np.nan     # reference without positives: NaN
    assert w.loc["C", "fold_vs_reference"] == pytest.approx(expected_ref, nan_ok=True)
    assert "fold_tumor_over_control" not in ba.positivity_table(a, ["alpha"], samples=["C", "T"], tumor_samples=[])


def test_config_validation_catches_crashing_settings():
    cfg = pc.default_config()
    cfg["clustering"]["stability_seeds"] = [0, 1, 1]
    cfg["steps"]["project_fine_bins"] = True
    cfg["bin_size"], cfg["qc_bin_size"] = "008um", "016um"
    cfg["samples"][1]["id"] = "ll-2-1"
    errors = " | ".join(pc.validate_config(cfg)[0])
    assert "seed twice" in errors and "qc_bin_size must be finer" in errors and "letter case" in errors
    ok = pc.default_config()
    ok["samples"][0]["id"] = "abc\n"
    assert any("sample id" in e for e in pc.validate_config(ok)[0])


def test_signatures_without_optional_fields_run(tiny_config, tmp_path):
    tiny_config["signatures"]["bare"] = {"core": ["Acp5", "Ctsk"], "extended": ["Acp5", "Ctsk", "Calcr"]}
    tiny_config["steps"] = {k: False for k in tiny_config["steps"]}
    a = ba.run_analysis_from_config(tiny_config, tmp_path)
    assert "score_bare" in a.obs and (tmp_path / "summary.json").is_file()


# ---------------------------------------------------------------------------------------------- depth matching
def test_downsampling_hits_the_target_and_is_independent_between_bins():
    X = sp.csr_matrix(np.tile(np.arange(1, 21, dtype=np.float32), (300, 1)))   # 300 identical bins, 210 UMI each
    D = ba.downsample_counts_matrix(X, 50, seed=0)
    assert np.all(np.asarray(D.sum(axis=1)).ravel() == 50)
    assert (D.toarray() <= X.toarray()).all()
    assert len({tuple(r) for r in D.toarray()}) > 250                            # scanpy kept one draw for all bins
    keep = ba.downsample_counts_matrix(X, 1000, seed=0)
    np.testing.assert_array_equal(keep.toarray(), X.toarray())


def test_quantile_depth_matching_matches_the_shallow_distribution():
    a = _toy(n=2000, depth=(1.0, 3.0))
    a.obs["total_counts"] = np.asarray(a.layers["counts"].sum(axis=1)).ravel()
    targets, target, med = ba.depth_match_targets(a, "quantile")
    shallow = (a.obs["sample"] == "S1").to_numpy()
    np.testing.assert_array_equal(targets[shallow], a.obs["total_counts"].to_numpy()[shallow])
    assert abs(np.median(targets[~shallow]) - np.median(targets[shallow])) <= 2
    cap, _, _ = ba.depth_match_targets(a, "cap")
    assert (cap == target).all()


def test_depth_matched_scores_come_from_the_downsampled_matrix():
    a = _toy(depth=(1.0, 4.0))
    a.obs["roi"] = True
    t = ba.depth_matched(a, ["alpha"], n_reps=2, roi_key="roi", sigs=SIGS, min_core=0, n_perm=3, n_sub=300)
    native = ba.positivity_table(ba.call_positive(ba.score_all(a.copy(), SIGS), ["alpha"], SIGS, min_core=0,
                                                  n_perm=3, n_sub=300), ["alpha"], roi_key="roi")
    assert {"n_positive_mean", "n_bins", "pct_positive_sd"} <= set(t.columns)
    deep = t[t["sample"] == "S2"].iloc[0]
    assert deep["mean_score"] != pytest.approx(native[native["sample"] == "S2"].iloc[0]["mean_score"])


# ---------------------------------------------------------------------------------------------- ROI and pooling
def test_roi_ignores_unscorable_exclusion_signatures():
    a = ba.score_all(_toy(), SIGS)
    a.obs["score_ghost"] = np.nan
    a = ba.define_marrow_roi(a, exclude=("ghost",), progress=pc.Progress())
    assert a.obs["in_marrow_roi"].all()
    a.obs.loc[a.obs_names[:5], "score_alpha"] = np.nan
    a = ba.define_marrow_roi(a, exclude=("alpha",))
    assert a.obs.loc[a.obs_names[:5], "in_marrow_roi"].all()


def test_neighbor_matrix_matches_brute_force_and_small_sections():
    from scipy.spatial import cKDTree
    a = _toy(n=200)
    P = ba.spatial_neighbor_matrix(a, 8)
    s1 = np.flatnonzero(a.obs["sample"] == "S1")
    _, idx = cKDTree(a.obsm["spatial"][s1]).query(a.obsm["spatial"][s1], k=9)
    for i in (0, 17, 99):
        assert set(P[s1[i]].indices) == set(s1[idx[i]])
    s2 = np.flatnonzero(a.obs["sample"] == "S2")
    tiny = a[a.obs_names[s1[:3]].tolist() + a.obs_names[s2].tolist()].copy()
    Pt = ba.spatial_neighbor_matrix(tiny, 8)                                    # 3 bins in S1: they pool each other
    assert Pt[0].nnz == 3


def test_pooled_scores_differ_from_native_scores():
    a = ba.call_positive(ba.score_all(_toy(), SIGS), ["alpha"], SIGS, n_perm=3, n_sub=300)
    a = ba.spatially_pooled_scores(a, pool_k=8, sigs=SIGS, names=["alpha"], n_perm=3, n_sub=300)
    assert not np.allclose(a.obs["pooled_score_alpha"], a.obs["score_alpha"])


# ---------------------------------------------------------------------------------------------- compartments
def test_control_osteoclasts_are_not_tumor_associated():
    a = ba.score_all(_toy(), {**SIGS, "osteoclast": SIGS["alpha"], "chondrocyte": SIGS["beta"]})
    a.obs["pos_osteoclast"] = True
    a = ba.split_osteoclast_compartments(a, 90)
    comp = a.obs.groupby("condition")["oc_compartment"].agg(lambda c: set(c.astype(str)))
    assert "tumor_associated" not in comp["control"]
    assert "non_growth_plate_control" in comp["control"] and "tumor_associated" in comp["tumor"]


def test_compartment_test_needs_two_bins_per_group():
    a = _toy()
    a.obs["oc_compartment"] = pd.Categorical(["negative"] * a.n_obs, categories=ba.COMPARTMENTS)
    a.obs.iloc[0, a.obs.columns.get_loc("oc_compartment")] = "growth_plate"
    progress = pc.Progress()
    de = ba.compare_osteoclast_states(a, progress=progress)
    assert de.empty and list(de.columns) == ba.COMPARTMENT_DE_COLUMNS and progress.warnings


# ---------------------------------------------------------------------------------------------- receptor panel
def test_receptor_panel_handles_missing_genes_empty_sections_and_reference():
    a = _toy(samples=("A", "B"))
    a.obs["roi"] = (a.obs["sample"] == "B").to_numpy()
    df = ba.receptor_panel(a, {"g1": "x", "nope": "y"}, n_reps=2, roi_key="roi", reference="B", samples=["A", "B"])
    assert set(df["gene"]) == {"g1"}
    assert np.isnan(df.loc[df["sample"] == "A", "pct_detected_dm"]).all()
    assert "expr_g1" not in a.obs


# ---------------------------------------------------------------------------------------------- spatial
def test_uncoupling_summary_has_every_column_when_not_estimable():
    a = _toy()
    a.obs["roi"] = True
    a.obs["A"] = False
    a.obs["B"] = False
    summ, dist = ba.uncoupling_metric(a, "A", "B", roi_key="roi", reference_sample="S1", n_perm=9)
    assert list(summ.columns) == ba.UNCOUPLING_COLUMNS
    assert summ["status"].str.startswith("not estimable").all() and dist.empty


def test_proximity_reports_co_positive_bins_and_can_exclude_self():
    a = _toy(n=200)
    a.obs["in_marrow_roi"] = True
    a.obs["pos_mac"] = False
    a.obs["pos_oc"] = False
    s1 = np.flatnonzero(a.obs["sample"] == "S1")
    a.obs.iloc[s1[0], a.obs.columns.get_loc("pos_mac")] = True
    a.obs.iloc[s1[0], a.obs.columns.get_loc("pos_oc")] = True                    # co-positive, far from others
    a.obs.iloc[s1[-1], a.obs.columns.get_loc("pos_oc")] = True
    t = ba.proximity_table(a, "mac", "oc", "positivity", 50.0, exclude_self=False, samples=["S1", "S2"])
    assert t.loc[0, "mac_also_oc_positive_bins"] == 1 and t.loc[0, "percent_mac_within_50um"] == 100.0
    t2 = ba.proximity_table(a, "mac", "oc", "positivity", 50.0, exclude_self=True, samples=["S1", "S2"])
    assert t2.loc[0, "percent_mac_within_50um"] == 0.0
    assert t2.loc[1, "oc_positive_bins"] == 0 and t2.loc[1, "percent_all_bins_within_50um"] == 0.0


def test_default_proximity_column_names_are_the_manuscript_ones():
    a = _toy(n=200)
    a.obs["in_marrow_roi"] = True
    a.obs["program"] = np.where(np.arange(a.n_obs) % 3 == 0, "fibroblast", "osteoblast")
    for n in ("macrophage", "osteoclast"):
        a.obs[f"pos_{n}"] = np.arange(a.n_obs) % 5 == 0
    fib, mac = ba.proximity_tables(a)
    assert list(fib.columns) == ["section", "fibroblast_assigned_bins", "osteoblast_assigned_bins",
                                 "region_of_interest_bins", "percent_fibroblast_within_50um_of_osteoblast",
                                 "percent_all_bins_within_50um_of_osteoblast"]
    assert list(mac.columns)[:7] == ["section", "osteoclast_positive_bins", "macrophage_positive_bins",
                                     "macrophage_within_50um_of_osteoclast", "percent_macrophage_within_50um",
                                     "region_of_interest_bins", "percent_all_bins_within_50um"]


def test_builtin_morans_i_equals_squidpy():
    sq = pytest.importorskip("squidpy")
    a = _toy(n=400)
    sub = a[a.obs["sample"] == "S1"].copy()
    sub.X = sub.layers["lognorm"]
    genes = GENES[:40]
    mine = ba._morans_i_builtin(sub, genes)
    sq.gr.spatial_neighbors(sub, coord_type="generic", delaunay=True)
    sq.gr.spatial_autocorr(sub, mode="moran", genes=genes, n_perms=None, n_jobs=1)
    ref = sub.uns["moranI"]
    np.testing.assert_allclose(mine.loc[ref.index, "I"], ref["I"], rtol=1e-6, atol=1e-9)
    ok = mine.loc[ref.index, "I"].notna()
    np.testing.assert_allclose(mine.loc[ref.index, "pval_norm"][ok], ref["pval_norm"][ok], rtol=1e-6)


def test_clustering_stability_cleans_up_uns():
    sc.settings.verbosity = 0
    a = _toy(n=300)
    sc.pp.pca(a, n_comps=10)
    a.uns["neighbors_rep"] = "X_pca"
    df = ba.clustering_stability(a, seeds=(0, 1), n_pcs=10)
    assert len(df) == 1 and not any(k.startswith(("_s", "n0", "n1")) for k in a.uns)


def test_projection_marks_missing_parents_separately():
    coarse = ad.AnnData(np.zeros((2, 1)), obs=pd.DataFrame({"sample": ["S", "S"], "array_row": [0, 0],
                                                             "array_col": [0, 1], "leiden": ["unassigned", "3"]}))
    fine = ad.AnnData(np.zeros((3, 1)), obs=pd.DataFrame({"sample": ["S"] * 3, "array_row": [0, 1, 5],
                                                           "array_col": [1, 2, 9]}))
    out = ba.map_domains_from_coarser(fine, coarse, 8, 16, keys=("leiden",))
    assert list(out.obs["leiden"].astype(str)) == ["unassigned", "3", "no_parent_bin"]


# ---------------------------------------------------------------------------------------------- drawn ROI
def test_rasterise_drawn_roi_scales_frames_and_respects_holes(tmp_path):
    a = ad.AnnData(np.zeros((4, 1)), obs=pd.DataFrame({"sample": ["S", "S", "S", "T"]}, index=list("abcd")))
    a.obsm["spatial"] = np.array([[100, 100], [500, 500], [150, 150], [100, 100]], dtype=float)  # full resolution
    a.uns["spatial"] = {"S": {"scalefactors": {"tissue_hires_scalef": 0.5}}}
    square = [[0, 0], [100, 0], [100, 100], [0, 100], [0, 0]]                    # hires = 0..200 full resolution
    hole = [[70, 70], [80, 70], [80, 80], [70, 80], [70, 70]]                    # hires 70..80 = 140..160 full res
    path = tmp_path / "S.geojson"
    path.write_text(json.dumps({"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"frame": "hires"}, "geometry": {"type": "Polygon",
                                                                           "coordinates": [square, hole]}}]}))
    keep, unrestricted = bpm.rasterise_drawn_roi(a, {"S": path}, "hires")
    assert keep.tolist() == [True, False, False, True] and unrestricted == ["T"]


# ---------------------------------------------------------------------------------------------- regression vs original
def test_native_results_equal_the_original_pipeline(tiny_demo, tiny_config, original_ba):
    """Native scores, null thresholds, positivity calls and the ROI were not affected by the defects: the fixed code
    must reproduce the manuscript code exactly on them. Run on the tiny demo sections."""
    sc.settings.verbosity = 0
    for s in tiny_demo["samples"]:
        original_ba.SAMPLES.setdefault(s["id"], {"condition": s["condition"], "limb": "-", "label": s["label"]})
    old = {}
    for s in tiny_demo["samples"]:
        p = pc.resolve_bin_dir(tiny_config, s)
        old[s["id"]] = original_ba.qc_filter(original_ba.qc_metrics(original_ba.read_visium_hd(p, s["id"])))
    A = ad.concat(old, label="sample_key", index_unique="-", uns_merge="unique", join="outer")
    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    A.layers["lognorm"] = A.X.copy()
    A.raw = A                                                                   # as the manuscript code did
    A = original_ba.call_positive(original_ba.score_all(A), n_perm=5)
    A = original_ba.define_marrow_roi(A)

    new, _, _ = ba.load_sections(tiny_config, "016um")
    B = ad.concat(new, label="sample_key", index_unique="-", uns_merge="unique", join="outer")
    B = ba.normalize_and_cluster(B, cluster=False)
    B = ba.call_positive(ba.score_all(B), n_perm=5)
    B = ba.define_marrow_roi(B)

    assert list(A.obs_names) == list(B.obs_names)
    for n in ba.SIG:
        np.testing.assert_array_equal(A.obs[f"score_{n}"].to_numpy(), B.obs[f"score_{n}"].to_numpy(), err_msg=n)
        np.testing.assert_array_equal(A.obs[f"ncore_{n}"].to_numpy(), B.obs[f"ncore_{n}"].to_numpy(), err_msg=n)
    for n in ba.PRIMARY:
        assert A.uns["positivity_thresholds"][n] == B.uns["positivity_thresholds"][n], n
        np.testing.assert_array_equal(A.obs[f"pos_{n}"].to_numpy(bool), B.obs[f"pos_{n}"].to_numpy(bool), err_msg=n)
    np.testing.assert_array_equal(A.obs["in_marrow_roi"].to_numpy(), B.obs["in_marrow_roi"].to_numpy())

    # identical inputs to the distance metric give identical results (same random-placement draws)
    for c in ("pos_osteoblast", "pos_osteoclast"):
        A.obs[c] = B.obs[c].to_numpy()
    ids = [s["id"] for s in tiny_demo["samples"]]
    old_sum, old_dist = original_ba.uncoupling_metric(A, "pos_osteoblast", "pos_osteoclast", n_perm=19, min_pos=3,
                                                      reference_sample=ids[0])
    new_sum, new_dist = ba.uncoupling_metric(B, "pos_osteoblast", "pos_osteoclast", n_perm=19, min_pos=3,
                                             reference_sample=ids[0], samples=ids)
    cols = ["ratio_obs_null_mean", "ratio_obs_null_lo", "ratio_obs_null_hi", "p_greater_than_null_mean"]
    np.testing.assert_allclose(new_sum[cols].to_numpy(float), old_sum[cols].to_numpy(float))
    pd.testing.assert_frame_equal(new_dist.reset_index(drop=True), old_dist.reset_index(drop=True), check_dtype=False)
