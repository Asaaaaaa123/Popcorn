"""Figure 4h sensitivity (bone_sensitivity/compartments.py and compartments_figures.py): synthetic checks of the
threshold diagnostics, the anatomical proxy, the depth adjustment and the pseudobulk units; a run on the tiny demo that
reproduces the manuscript test; the figure."""
from __future__ import annotations

import copy
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import pipeline_config as pc
from bone_sensitivity import common
from bone_sensitivity import compartments as cm
from bone_sensitivity import compartments_figures as cmf

GP, TA = cm.GROUPS


# ---- synthetic -------------------------------------------------------------------------------------------------------

def test_threshold_collapse_on_zero_inflated_score():
    rng = np.random.default_rng(0)
    score = np.r_[np.zeros(445), -0.01 - rng.random(475), 0.05 + rng.random(80)]   # 44.5% zero, 8% positive
    thr, above, zero = cm.effective_threshold(score, 90)
    assert thr == 0.0 and zero
    assert np.array_equal(above, score > 0) and above.sum() == 80
    thr, above, zero = cm.effective_threshold(rng.normal(size=1000), 90)
    assert not zero and thr > 0 and above.sum() == 100
    thr, above, zero = cm.effective_threshold(np.full(4, np.nan), 90)
    assert np.isnan(thr) and not above.any() and not zero


def test_connected_regions_eight_connectivity_and_minimum_size():
    r, c = np.meshgrid(np.arange(12), np.arange(12), indexing="ij")
    r, c = r.ravel(), c.ravel()
    block = (r >= 2) & (r <= 4) & (c >= 2) & (c <= 5)                       # 12 bins
    diag = (r == 5) & (c == 6)                                              # touches the block at a corner only
    small = (r == 9) & (c >= 9) & (c <= 10)                                 # 2 bins
    mask = block | diag | small
    keep = mask | ((r + c) % 3 != 0)                                        # QC holes outside the regions
    label, sizes = cm.connected_regions(r[keep], c[keep], mask[keep], min_bins=5)
    assert sizes == [13]
    lab = np.full(len(r), -9)
    lab[keep] = label
    assert set(lab[block | diag]) == {0} and (lab[small] == -1).all() and (lab[keep & ~mask] == -1).all()
    _, sizes = cm.connected_regions(r, c, mask, min_bins=2)
    assert sizes == [13, 2]
    label, sizes = cm.connected_regions(r, c, np.zeros(len(r), bool), min_bins=1)
    assert sizes == [] and (label == -1).all()


def test_near_regions_distance_rule():
    rows = np.r_[np.arange(10), [5, 5, 5, 5]]
    cols = np.r_[np.zeros(10), [2, 6, 7, 1]]
    xy = np.column_stack([cols, rows]) * 16.0
    region = np.r_[np.ones(10, bool), np.zeros(4, bool)]
    query = np.r_[np.zeros(10, bool), [True, True, True, False]]
    near = cm.near_regions(xy, region, query, 100.0)                        # 32, 96, 112 um; the last is not queried
    assert near.tolist() == [False] * 10 + [True, True, False, False]
    assert cm.near_regions(xy, region, query, 96.0)[11]                     # inclusive
    assert not cm.near_regions(xy, np.zeros(14, bool), query, 1e9).any()


def test_anatomical_growth_plate_stays_within_a_section():
    r, c = np.meshgrid(np.arange(10), np.arange(10), indexing="ij")
    r, c = np.tile(r.ravel(), 2), np.tile(c.ravel(), 2)
    smp = np.repeat(["A", "B"], 100)
    dense = (smp == "A") & (r <= 2)                                          # 30-bin band in section A only
    oc = (r == 5) & (c == 4)                                                 # one bin per section, 48 um from it
    a = ad.AnnData(obs=pd.DataFrame({"sample": smp, "array_row": r, "array_col": c, "pooled_pos_chondrocyte": dense},
                                    index=[f"b{i}" for i in range(200)]))
    near, regions, calls = cm.anatomical_growth_plate(a, oc, np.zeros(200, bool), "chondrocyte", 50.0, 20, 16.0,
                                                      ["A", "B"])
    assert calls == "pooled_pos_chondrocyte"
    assert near.tolist() == ((smp == "A") & oc).tolist()
    assert regions["n_regions"].tolist() == [1, 0] and regions["region_sizes"].tolist() == ["30", ""]
    near, _, _ = cm.anatomical_growth_plate(a, oc, np.zeros(200, bool), "chondrocyte", 40.0, 20, 16.0, ["A", "B"])
    assert not near.any()
    _, regions, _ = cm.anatomical_growth_plate(a, oc, np.zeros(200, bool), "chondrocyte", 50.0, 31, 16.0, ["A", "B"])
    assert regions["n_regions"].tolist() == [0, 0]


def test_depth_adjusted_targets_follow_the_shallower_group():
    import bone_analysis as ba
    rng = np.random.default_rng(1)
    deep, shallow = rng.poisson(750, 400), rng.poisson(350, 300)
    for first, second in ((GP, TA), (TA, GP)):
        totals, labels = np.r_[deep, shallow], np.array([first] * 400 + [second] * 300)
        targets, deeper = cm.depth_adjusted_targets(totals, labels)
        assert deeper == first
        assert np.array_equal(targets[400:], shallow)                        # the shallower group keeps its counts
        t = targets[:400]
        assert (t <= deep).all()
        qs = np.linspace(.02, .98, 49)
        assert np.abs(np.quantile(t, qs) - np.quantile(shallow, qs)).max() <= 3
        assert abs(np.median(t) - np.median(shallow)) <= 2
        assert np.corrcoef(pd.Series(t).rank(), pd.Series(deep).rank())[0, 1] > .99      # ranks are kept
    counts = sp.csr_matrix(np.vstack([rng.multinomial(n, np.full(30, 1 / 30)) for n in np.r_[deep, shallow]]))
    targets, _ = cm.depth_adjusted_targets(np.r_[deep, shallow], np.array([GP] * 400 + [TA] * 300))
    X = ba.downsample_counts_matrix(counts, targets, 0)
    assert np.array_equal(np.asarray(X.sum(axis=1)).ravel().astype(int), targets)
    assert np.array_equal(cm.quantile_targets([5, 1], [100, 200]), [5, 1])     # never above the own total


def test_block_units_sum_counts_per_section_block_and_compartment():
    sections = ["A", "A", "A", "A", "B", "B"]
    xy = [[0, 0], [16, 0], [200, 0], [32, 16], [0, 0], [16, 16]]
    labels = [GP, TA, GP, GP, GP, TA]
    unit, ulab, usec = cm.block_units(sections, xy, labels, 128.0)
    assert len(ulab) == 5 and unit[0] == unit[3] and len(set(unit[[0, 1, 2, 4, 5]])) == 5
    assert ulab[unit[1]] == TA and usec[unit[4]] == "B"
    counts = sp.csr_matrix(np.arange(12, dtype=float).reshape(6, 2))
    U = cm.sum_units(counts, unit, len(ulab)).toarray()
    assert np.array_equal(U[unit[0]], [0 + 6, 1 + 7]) and np.array_equal(U[unit[5]], [10, 11])


def test_block_pseudobulk_balances_unit_depth_under_the_null():
    """Identical expression profile; the growth-plate units sum 2 deep bins each, the tumor-associated units 15
    shallow bins each, so the tumor-associated units are ~3.5x deeper. Unbalanced units call many genes 'higher in
    tumor-associated'; balanced units call none."""
    rng = np.random.default_rng(0)
    n_genes, n_gp, k_gp, n_ta, k_ta = 2000, 50, 2, 150, 15
    p = rng.gamma(0.3, size=n_genes)
    p /= p.sum()
    labels = np.array([GP] * (n_gp * k_gp) + [TA] * (n_ta * k_ta))
    unit = np.r_[np.repeat(np.arange(n_gp), k_gp), n_gp + np.repeat(np.arange(n_ta), k_ta)]
    within = np.r_[np.tile(np.arange(k_gp), n_gp), np.tile(np.arange(k_ta), n_ta)]
    xy = np.column_stack([unit * 128.0 + (within % 8) * 16.0, (within // 8) * 16.0])   # one 128 um block per unit
    depth = np.r_[rng.poisson(750, n_gp * k_gp), rng.poisson(350, n_ta * k_ta)]
    counts = sp.csr_matrix(rng.poisson(np.outer(depth, p)).astype(np.float32))
    var, secs = pd.Index([f"g{i}" for i in range(n_genes)]), np.full(len(labels), "S")
    u, ulab, _ = cm.block_units(secs, xy, labels, 128.0)
    U = cm.sum_units(counts, u, len(ulab))
    raw = cm.de_table(cm.log_normalize(U, 1e4), U, ulab, var, 1e4, [])
    assert (cm.significant(raw) & (raw["group"] == TA)).sum() >= 50                # the confound the balance removes
    de, info = cm.block_pseudobulk(counts, secs, xy, labels, 128.0, var, 1e4, [], 3, 0)
    assert info["n"] == {GP: n_gp, TA: n_ta} and info["bins_per_unit"] == {GP: k_gp, TA: k_ta}
    assert info["deeper"] == TA and info["umi"][TA] > 3 * info["umi"][GP]
    assert abs(info["target"][TA] - info["target"][GP]) <= 0.02 * info["target"][GP]
    assert info["target"][GP] == info["umi"][GP]                                    # the shallower units keep theirs
    assert cm.significant(de).sum() <= 2 and de["frac_reps_significant"].between(0, 1).all()
    assert len(de) == 2 * n_genes
    none, info = cm.block_pseudobulk(counts[:k_gp + 1], secs[:k_gp + 1], xy[:k_gp + 1], labels[:k_gp + 1], 128.0,
                                     var, 1e4, [], 3, 0)
    assert none is None and info["n"] == {GP: 2, TA: 0}


def test_de_table_statistics_on_known_counts():
    counts = sp.csr_matrix(np.array([[9, 1, 0], [8, 2, 0], [10, 0, 1], [1, 9, 1], [0, 10, 1], [2, 8, 1]], float))
    labels = np.array([GP] * 3 + [TA] * 3)
    var = pd.Index(["g0", "g1", "g2"])
    de = cm.de_table(cm.log_normalize(counts, 1e4), counts, labels, var, 1e4, ["g2"])
    cp10k = counts.toarray() / counts.toarray().sum(1, keepdims=True) * 1e4
    mg, mt = cp10k[:3].mean(0), cp10k[3:].mean(0)
    d = de.set_index(["group", "names"])
    for j, g in enumerate(var):
        assert d.loc[(GP, g), "log2fc_mean_cp10k"] == pytest.approx(np.log2((mg[j] + 1) / (mt[j] + 1)))
        assert d.loc[(TA, g), "log2fc_mean_cp10k"] == pytest.approx(np.log2((mt[j] + 1) / (mg[j] + 1)))
        assert d.loc[(GP, g), "scores"] == pytest.approx(-d.loc[(TA, g), "scores"])
    assert d.loc[(GP, "g2"), "pct_group"] == pytest.approx(100 / 3) and d.loc[(GP, "g2"), "pct_rest"] == 100
    assert d.loc[(GP, "g0"), "scores"] > 0 and d["selection_gene"].groupby(level="names").first().tolist() == \
        [False, False, True]
    reps = [de.assign(pvals_adj=p) for p in (0.01, 0.01, 0.5)]
    agg = cm.aggregate_replicates(reps).set_index(["group", "names"])
    assert agg.loc[(GP, "g0"), "frac_reps_significant"] == pytest.approx(2 / 3)
    assert agg.loc[(TA, "g0"), "frac_reps_significant"] == 0          # significant only towards the growth plate


# ---- run on the tiny demo --------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def tiny_run(tiny_joint):
    a = tiny_joint["a"].copy()
    outs = cm.run(a, tiny_joint["cfg"])
    return a, {o.file: o.df for o in outs}, outs


def test_run_outputs_and_columns(tiny_run):
    _, t, outs = tiny_run
    assert [o.file for o in outs] == [cm.DIAG_CSV, cm.REGIONS_CSV, cm.VARIANTS_CSV, cm.DE_CSV, cm.TOP_CSV]
    assert list(t[cm.DIAG_CSV].columns) == cm.DIAG_COLUMNS
    assert list(t[cm.REGIONS_CSV].columns) == cm.REGION_COLUMNS
    assert list(t[cm.VARIANTS_CSV].columns) == cm.VARIANT_COLUMNS
    assert list(t[cm.DE_CSV].columns) == cm.DE_COLUMNS
    assert list(t[cm.TOP_CSV].columns[:3]) == ["gene", "group", "selection_gene"]
    assert {f"padj_{v}" for v in cm.VARIANTS} | {f"significant_{v}" for v in cm.VARIANTS} <= set(t[cm.TOP_CSV].columns)
    assert t[cm.VARIANTS_CSV]["variant"].tolist() == list(cm.VARIANTS)
    status = t[cm.VARIANTS_CSV].set_index("variant")["status"]
    assert status["drawn"] == "not configured"
    assert (status.drop("drawn") == "estimable").all()
    assert all(o.title and o.description for o in outs)


def test_run_leaves_the_object_unchanged(tiny_joint, tiny_run):
    a, _, _ = tiny_run
    ref = tiny_joint["a"]
    assert list(a.obs.columns) == list(ref.obs.columns) and a.obs.equals(ref.obs)
    assert sorted(a.uns.keys()) == sorted(ref.uns.keys()) and sorted(a.layers.keys()) == sorted(ref.layers.keys())
    assert (a.layers["counts"] != ref.layers["counts"]).nnz == 0


def test_manuscript_variant_reproduces_figure4h(tiny_joint, tiny_run):
    a, t, _ = tiny_run
    cfg = tiny_joint["cfg"]
    man = pd.read_csv(Path(tiny_joint["out"]) / "Figure4h_osteoclast_compartment_DE.csv")
    de = t[cm.DE_CSV]
    got = de[de["variant"] == "manuscript"].merge(man, on=["group", "names"], suffixes=("", "_ms"), validate="1:1")
    assert len(got) == (de["variant"] == "manuscript").sum() > 0
    np.testing.assert_allclose(got["scores"], got["scores_ms"], rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(got["pvals_adj"], got["pvals_adj_ms"], rtol=1e-9, atol=1e-12)
    sig = man[man["pvals_adj"] < 0.05]
    assert set(zip(sig["group"], sig["names"])) <= set(zip(got["group"], got["names"]))
    comp = a.obs["oc_compartment"].astype(str)
    assert np.array_equal(cm.manuscript_compartments(a, cfg["compartments"]), comp.to_numpy())
    row = t[cm.VARIANTS_CSV].set_index("variant").loc["manuscript"]
    assert cfg["compartments"]["de_sections"] == "all"
    assert row["n_group"] == (comp == GP).sum() and row["n_rest"] == (comp == TA).sum()
    n_sig = ((man["group"] == GP) & (man["pvals_adj"] < .05) & (man["scores"] > 0)).sum()
    assert row["n_significant_growth_plate"] == n_sig
    top = t[cm.TOP_CSV]
    for g in (GP, TA):
        expected = man[man["group"] == g].nlargest(cm.FIGURE_TOP_N, "scores")["names"].tolist()
        assert top.loc[top["group"] == g, "gene"].tolist() == expected
    np.testing.assert_allclose(top["padj_manuscript"],
                               man.set_index(["group", "names"])["pvals_adj"].reindex(
                                   pd.MultiIndex.from_frame(top[["group", "gene"]])).to_numpy(), rtol=1e-9)


def test_diagnostics_regions_and_variants_are_consistent(tiny_joint, tiny_run):
    a, t, _ = tiny_run
    cfg = tiny_joint["cfg"]
    st = common.settings(cfg, "compartments")
    d = t[cm.DIAG_CSV].set_index("sample")
    score = a.obs[f"score_{cfg['compartments']['split_by']}"].to_numpy(float)
    assert d.loc["all", "threshold"] == pytest.approx(np.nanpercentile(score, cfg["compartments"]["percentile"]))
    assert bool(d.loc["all", "threshold_is_zero"]) == (d.loc["all", "threshold"] == 0)
    comp = a.obs["oc_compartment"].astype(str)
    assert d.loc["all", "n_growth_plate"] == (comp == GP).sum()
    assert d.loc["all", "n_tumor_associated"] == (comp == TA).sum()
    secs = d.drop("all")
    for col in ("n_growth_plate", "n_tumor_associated", "n_non_growth_plate_control", "n_bins",
                "n_bins_with_split_gene_umi", "n_growth_plate_outside_roi"):
        assert secs[col].sum() == d.loc["all", col]
    roi = a.obs["in_marrow_roi"].astype(bool).to_numpy()
    assert d.loc["all", "n_growth_plate_outside_roi"] == ((comp == GP).to_numpy() & ~roi).sum()
    umi = a.obs["total_counts"].to_numpy(float)
    assert d.loc["all", "median_umi_growth_plate"] == pytest.approx(np.median(umi[(comp == GP).to_numpy()]))
    reg = t[cm.REGIONS_CSV]
    assert reg["sample"].tolist() == common.sample_ids(a, cfg) and (reg["calls"] == "pooled_pos_chondrocyte").all()
    for sizes in reg["region_sizes"].astype(str):
        assert all(int(s) >= st["min_region_bins"] for s in sizes.split(";") if s)
    v = t[cm.VARIANTS_CSV].set_index("variant")
    tumor_secs, ref = common.tumor_and_reference(cfg)
    on_tumor = reg["sample"].isin(tumor_secs)
    oc = a.obs["pos_osteoclast"].astype(bool).to_numpy()
    tumor = a.obs["condition"].astype(str).isin(cfg["compartments"]["tumor_conditions"]).to_numpy()
    assert v.loc["anatomical", "n_group"] == reg.loc[on_tumor, "n_osteoclast_bins_near_regions"].sum()  # tumor only
    assert v.loc["anatomical", "n_rest"] == (oc & tumor).sum() - v.loc["anatomical", "n_group"]
    n_ref_near = int(reg.loc[reg["sample"] == ref, "n_osteoclast_bins_near_regions"].sum())
    assert n_ref_near > 0 and f"left out, not tumor sections: {ref} {n_ref_near}" in v.loc["anatomical", "note"]
    n_ref_gp = int(d.loc[ref, "n_growth_plate"])
    assert f"{n_ref_gp} of {int(d.loc['all', 'n_growth_plate'])} from non-tumor sections" in v.loc["manuscript", "note"]
    assert v.loc["roi_only", "n_group"] == d.loc["all", "n_growth_plate"] - d.loc["all", "n_growth_plate_outside_roi"]
    dm = v.loc["depth_adjusted"]
    assert dm["n_group"] == v.loc["manuscript", "n_group"]
    assert abs(dm["median_umi_group"] - dm["median_umi_rest"]) <= 0.05 * max(dm["median_umi_group"], 1) + 2
    assert dm["median_umi_group_unadjusted"] == v.loc["manuscript", "median_umi_group"]
    bp = v.loc["block_pseudobulk"]
    assert bp["units"].startswith("blocks of")
    m = (a.obs["oc_compartment"].astype(str).isin(cm.GROUPS)).to_numpy()
    u, ulab, _ = cm.block_units(a.obs["sample"].astype(str)[m], common.array_xy_um(a, pc.bin_um(cfg["bin_size"]))[m],
                                a.obs["oc_compartment"].astype(str)[m], st["pseudobulk_block_um"])
    uumi = np.asarray(cm.sum_units(a.layers["counts"][m], u, len(ulab)).sum(axis=1)).ravel()
    assert (bp["n_group"], bp["n_rest"]) == ((ulab == GP).sum(), (ulab == TA).sum())
    assert bp["median_umi_group_unadjusted"] == pytest.approx(np.median(uumi[ulab == GP]))
    assert bp["median_umi_rest_unadjusted"] == pytest.approx(np.median(uumi[ulab == TA]))
    assert abs(bp["median_umi_group"] - bp["median_umi_rest"]) <= 0.1 * max(bp["median_umi_group"], 1) + 2
    low = "group" if bp["median_umi_group_unadjusted"] <= bp["median_umi_rest_unadjusted"] else "rest"
    assert bp[f"median_umi_{low}"] == bp[f"median_umi_{low}_unadjusted"]          # the shallower units keep theirs
    assert max(bp["median_umi_group_unadjusted"], bp["median_umi_rest_unadjusted"]) > 2 * bp[f"median_umi_{low}"]
    de = t[cm.DE_CSV]
    reps = de["variant"].isin(["depth_adjusted", "block_pseudobulk"])
    assert de.loc[reps, "frac_reps_significant"].between(0, 1).all()
    assert de.loc[~reps, "frac_reps_significant"].isna().all()
    sel = set(cfg["signatures"][cfg["compartments"]["split_by"]]["extended"])
    assert (de["selection_gene"] == de["names"].isin(sel)).all()


def test_drawn_and_not_estimable_variants(tiny_joint):
    a = tiny_joint["a"].copy()
    cfg = copy.deepcopy(tiny_joint["cfg"])
    ctrl, tum = cfg["reference_sample"], common.tumor_and_reference(cfg)[0][0]
    oc = a.obs["pos_osteoclast"].astype(bool).to_numpy()
    tumor = a.obs["condition"].astype(str).isin(cfg["compartments"]["tumor_conditions"]).to_numpy()
    rings, left = {}, {}
    for sid, half in ((ctrl, False), (tum, True)):        # the whole control section; the left half of the tumor one
        m = (a.obs["sample"].astype(str) == sid).to_numpy()
        f = float(a.uns["spatial"][sid]["scalefactors"]["tissue_hires_scalef"])
        xy = a.obsm["spatial"][m] * f
        lo, hi = xy.min(0) - 5, xy.max(0) + 5
        if half:
            xs = np.unique(xy[oc[m], 0])
            hi[0] = (xs[len(xs) // 2 - 1] + xs[len(xs) // 2]) / 2
        left[sid] = np.zeros(a.n_obs, bool)
        left[sid][m] = xy[:, 0] < hi[0]
        rings[sid] = {"type": "Polygon", "coordinates": [[[lo[0], lo[1]], [hi[0], lo[1]], [hi[0], hi[1]],
                                                          [lo[0], hi[1]], [lo[0], lo[1]]]]}
    cfg["sensitivity_analyses"]["compartments"].update(growth_plate_masks=rings, min_region_bins=10 ** 6)
    s = a.obs["score_chondrocyte"].to_numpy(float)                       # zero-inflated: the 90th percentile is 0
    s = np.where(s > np.quantile(s, .93), s, np.where(s > np.median(s), 0.0, -np.abs(s) - .01))
    a.obs["score_chondrocyte"] = s
    progress = pc.Progress()
    t = {o.file: o.df for o in cm.run(a, cfg, progress)}
    d = t[cm.DIAG_CSV].set_index("sample").loc["all"]
    assert d["threshold"] == 0 and d["threshold_is_zero"] and d["n_bins_above_threshold"] == (s > 0).sum()
    assert d["n_growth_plate"] == (a.obs["pos_osteoclast"].astype(bool).to_numpy() & (s > 0)).sum()
    assert any("exactly 0" in w for w in progress.warnings)
    assert any("anatomical" in w and "not estimable" in w for w in progress.warnings)
    assert "exactly 0" in t[cm.VARIANTS_CSV].set_index("variant").loc["manuscript", "note"]
    v = t[cm.VARIANTS_CSV].set_index("variant")
    assert v.loc["drawn", "status"] == "estimable"                 # the control polygon is left out, not tested
    assert v.loc["drawn", "n_group"] == (oc & left[tum]).sum() >= 2
    assert v.loc["drawn", "n_rest"] == (oc & tumor & ~left[tum]).sum() >= 2
    assert f"left out, not tumor sections: {ctrl} {(oc & left[ctrl]).sum()}" in v.loc["drawn", "note"]
    assert v.loc["anatomical", "status"] == "not estimable" and v.loc["anatomical", "n_group"] == 0
    assert (t[cm.REGIONS_CSV]["n_regions"] == 0).all()
    top = t[cm.TOP_CSV]
    assert top["padj_anatomical"].isna().all() and top["padj_drawn"].notna().all()


# ---- figure ----------------------------------------------------------------------------------------------------------

def test_draw_on_run_tables(tiny_joint, tiny_run, tmp_path):
    _, _, outs = tiny_run
    res = tmp_path / "results"
    res.mkdir()
    for o in outs:
        o.df.to_csv(res / o.file, index=False)
    assert set(cmf.INPUTS) <= {o.file for o in outs}
    entries = cmf.draw(res, tiny_joint["cfg"], tmp_path / "figures")
    assert [e["name"] for e in entries] == [n for n, _ in cmf.FIGURES]
    e = entries[0]
    assert set(e) == {"name", "group", "title", "caption", "png", "pdf", "error"} and e["group"] == common.GROUP
    assert (tmp_path / "figures" / e["png"]).is_file() and (tmp_path / "figures" / e["pdf"]).is_file()
    diag = outs[0].df.set_index("sample").loc["all"]
    assert f"{int(diag['n_growth_plate']):,} growth-plate bins" in e["caption"]
    assert f"is {diag['threshold']:.4g}" in e["caption"]
    assert "Not run: drawn growth plate" in e["caption"]
    for v in ("Manuscript", "ROI only", "Depth-adjusted", "Block pseudobulk", "Anatomical growth plate"):
        assert f"{v}:" in e["caption"]
    bp = outs[2].df.set_index("variant").loc["block_pseudobulk"]
    assert (f"median UMI per unit, growth plate vs tumor-associated: {bp['median_umi_group_unadjusted']:g} vs "
            f"{bp['median_umi_rest_unadjusted']:g} before, {bp['median_umi_group']:g} vs {bp['median_umi_rest']:g} "
            f"after the adjustment") in e["caption"]
    n_ref = int(outs[0].df.set_index("sample").loc[tiny_joint["cfg"]["reference_sample"], "n_growth_plate"])
    n_all = int(diag["n_growth_plate"])
    assert f"growth-plate bins of every section ({n_ref} of {n_all} from non-tumor" in e["caption"]
    assert "anatomical and drawn variants take both compartments from tumor sections only" in e["caption"]
    assert "anatomical growth-plate bins of non-tumor sections left out" in e["caption"]


def test_draw_survives_not_estimable_tables(tiny_joint, tmp_path):
    res = tmp_path / "results"
    res.mkdir()
    pd.DataFrame([{"variant": v, "status": "not estimable", "units": "bins", "n_group": 1, "n_rest": 0,
                   "median_umi_group": np.nan, "median_umi_rest": np.nan, "n_significant_growth_plate": np.nan,
                   "n_significant_tumor_associated": np.nan, "note": "fewer than two bins"} for v in cm.VARIANTS]
                 ).to_csv(res / cm.VARIANTS_CSV, index=False)
    pd.DataFrame(columns=["gene", "group", "selection_gene"]).to_csv(res / cm.TOP_CSV, index=False)
    entries = cmf.draw(res, tiny_joint["cfg"], tmp_path / "figures")
    assert len(entries) == 1 and (tmp_path / "figures" / entries[0]["png"]).is_file()
    assert "not estimable" in entries[0]["caption"]
    assert cmf.draw(tmp_path / "empty", tiny_joint["cfg"], tmp_path / "figures2") == []


def test_figure_module_needs_no_scanpy():
    import re
    src = Path(cmf.__file__).read_text(encoding="utf-8")
    imported = set(re.findall(r"^\s*(?:import|from)\s+([\w.]+)", src, flags=re.M))
    assert imported <= {"__future__", "numpy", "pandas", "matplotlib.figure", "matplotlib.lines", "pipeline_config",
                        "bone_sensitivity"}
