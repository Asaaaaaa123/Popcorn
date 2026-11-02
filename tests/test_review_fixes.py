"""Tests of the reporting and code fixes of the second scientific review (2026-11-02): exact folds, Wilson columns of
the depth-matched table, the Figure 4h table and compartment summary, Moran's I FDR, the unresolved share per section,
gene-set sensitivity gene lists, and the captions of Figures 3a, 4g/6a, 4h, 4i and Supplementary S1/S5."""
from __future__ import annotations

import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

import bone_analysis as ba
import bone_figures as bf
import pipeline_config as pc
from tests.test_figures import make_results

GENES = [f"g{i}" for i in range(60)]


def _obs_only(obs):
    """AnnData holding only `obs` (positivity_table, compartment_summary and the unresolved table read obs only)."""
    obs = obs.copy()
    obs.index = [f"b{i}" for i in range(len(obs))]
    return ad.AnnData(np.zeros((len(obs), 1), dtype=np.float32), obs=obs)


def _toy(n_per=200, seed=0, samples=("S1", "S2")):
    """Two sections of Poisson counts on a grid, normalized like the pipeline (layers counts and lognorm)."""
    rng = np.random.default_rng(seed)
    X = np.vstack([rng.poisson(2.0 + i, size=(n_per, len(GENES))) for i in range(len(samples))]).astype(np.float32)
    a = ad.AnnData(sp.csr_matrix(X), obs=pd.DataFrame({"sample": np.repeat(samples, n_per)},
                                                      index=[f"b{i}" for i in range(X.shape[0])]),
                   var=pd.DataFrame(index=GENES))
    a.obs["sample"] = pd.Categorical(a.obs["sample"], categories=list(samples))
    a.obs["condition"] = np.where(a.obs["sample"] == samples[0], "control", "tumor")
    a.obs["total_counts"] = np.asarray(a.X.sum(axis=1)).ravel()
    side = int(np.ceil(np.sqrt(n_per)))
    grid = np.array([(i % side, i // side) for i in range(n_per)], dtype=float)
    a.obsm["spatial"] = np.vstack([grid * 27.0 for _ in samples])
    a.obs["array_col"] = np.tile(grid[:, 0].astype(int), len(samples))
    a.obs["array_row"] = np.tile(grid[:, 1].astype(int), len(samples))
    a.uns["microns_per_pixel"] = {s: 0.59 for s in samples}
    return ba.normalize_and_cluster(a, cluster=False)


# ---------------------------------------------------------------------------------------------- item 2: exact folds
def test_folds_are_ratios_of_exact_rates_not_of_rounded_percentages():
    n_c, n_t, k_c, k_t = 3000, 7000, 1, 400
    pos = np.r_[np.arange(n_c) < k_c, np.arange(n_t) < k_t]
    a = _obs_only(pd.DataFrame({"sample": ["C"] * n_c + ["T"] * n_t, "pos_x": pos, "score_x": 0.0}))
    t = ba.positivity_table(a, ["x"], reference="C", samples=["C", "T"], tumor_samples=["T"]).set_index("sample")
    exact = (k_t / n_t) / (k_c / n_c)                                             # 171.4286
    rounded = t.loc["T", "pct_positive"] / t.loc["C", "pct_positive"]              # 5.714 / 0.033 = 173.15
    assert abs(exact - rounded) > 1
    assert t.loc["T", "fold_tumor_over_control"] == pytest.approx(round(exact, 2), abs=1e-9)
    assert t.loc["T", "fold_vs_reference"] == pytest.approx(round(exact, 4), abs=1e-9)
    assert t.loc["C", "fold_vs_reference"] == 1.0
    assert t.loc["C", "pct_positive"] == round(100 * k_c / n_c, 3)                 # stored percentages still rounded


# ---------------------------------------------------------------------------------------------- item 3: Wilson columns
def test_depth_matched_table_carries_wilson_intervals_of_the_mean_count():
    a = _toy()
    a.obs["roi"] = True
    sigs = {"alpha": {"core": GENES[:4], "extended": GENES[:6], "competes": [], "excluded": {}}}
    t = ba.depth_matched(a, ["alpha"], n_reps=2, roi_key="roi", sigs=sigs, min_core=0, n_perm=3, n_sub=300)
    assert list(t.columns) == ["signature", "sample", "pct_positive_mean", "pct_positive_sd", "n_positive_mean",
                               "n_bins", "mean_score", "ci_lo", "ci_hi"]
    for r in t.itertuples():
        lo, hi = proportion_confint(int(round(r.n_positive_mean)), int(r.n_bins), method="wilson")
        assert r.ci_lo == pytest.approx(round(100 * lo, 4)) and r.ci_hi == pytest.approx(round(100 * hi, 4))


# ---------------------------------------------------------------------------------------------- items 4, 5: Figure 4h
def _compartment_toy():
    """20 bins: 10 growth-plate, 10 tumor-associated; gA only in growth-plate bins, gB only in tumor-associated ones."""
    rng = np.random.default_rng(1)
    X = rng.poisson(3.0, size=(20, len(GENES))).astype(np.float32)
    X[10:, 0], X[:10, 0] = 0, rng.integers(5, 15, 10)                             # gene g0: growth-plate only
    X[:10, 1], X[10:, 1] = 0, rng.integers(5, 15, 10)                             # gene g1: tumor-associated only
    a = ad.AnnData(sp.csr_matrix(X), obs=pd.DataFrame({"sample": "T", "condition": "tumor"},
                                                      index=[f"b{i}" for i in range(20)]),
                   var=pd.DataFrame(index=GENES))
    a = ba.normalize_and_cluster(a, cluster=False)
    a.obs["oc_compartment"] = pd.Categorical(["growth_plate"] * 10 + ["tumor_associated"] * 10,
                                             categories=ba.COMPARTMENTS)
    return a


def test_compartment_de_reports_mean_cp10k_fold_detection_and_selection_genes():
    a = _compartment_toy()
    de = ba.compare_osteoclast_states(a, selection_genes=("g0",))
    assert list(de.columns) == ba.COMPARTMENT_DE_COLUMNS and "logfoldchanges" not in de
    cp = np.expm1(a.layers["lognorm"].toarray().astype(np.float64))
    gp = np.arange(20) < 10
    row = de[(de["group"] == "growth_plate") & (de["names"] == "g0")].iloc[0]
    assert row["log2fc_mean_cp10k"] == pytest.approx(np.log2((cp[gp, 0].mean() + 1) / (cp[~gp, 0].mean() + 1)),
                                                     abs=1e-4)
    assert row["pct_group"] == 100.0 and row["pct_rest"] == 0.0 and bool(row["selection_gene"])
    other = de[(de["group"] == "tumor_associated") & (de["names"] == "g1")].iloc[0]
    assert other["pct_group"] == 100.0 and other["pct_rest"] == 0.0 and not other["selection_gene"]
    assert de["selection_gene"].sum() == 2                                         # g0 in each group's rows
    # scanpy's logfoldchanges of the gene undetected in the other compartment is the 1e-9-floor artefact
    sub = a.copy()
    sub.X = sub.layers["lognorm"]
    sub.obs["oc_compartment"] = sub.obs["oc_compartment"].cat.remove_unused_categories()
    sc.tl.rank_genes_groups(sub, groupby="oc_compartment", method="wilcoxon", use_raw=False)
    ref = sc.get.rank_genes_groups_df(sub, group=None)
    assert ref.loc[(ref["group"] == "growth_plate") & (ref["names"] == "g0"), "logfoldchanges"].iloc[0] > 19
    assert abs(row["log2fc_mean_cp10k"]) < 19
    merged = de.merge(ref, on=["group", "names"], suffixes=("", "_ref"))
    np.testing.assert_allclose(merged["scores"], merged["scores_ref"], rtol=1e-6)    # the test itself is unchanged
    np.testing.assert_allclose(merged["pvals_adj"], merged["pvals_adj_ref"], rtol=1e-6)


def test_split_stores_the_threshold_and_warns_when_it_is_zero():
    n = 200
    score = np.where(np.arange(n) % 20 == 0, 0.5, 0.0)                             # 95% exactly zero
    a = _obs_only(pd.DataFrame({"sample": "T", "condition": "tumor", "score_chondrocyte": score,
                                "pos_osteoclast": np.arange(n) % 2 == 0}))
    progress = pc.Progress()
    a = ba.split_osteoclast_compartments(a, 90, progress=progress)
    assert a.uns["compartment_threshold"] == 0.0
    np.testing.assert_array_equal(a.obs["growth_plate"].to_numpy(), score > 0)
    assert any("is 0" in w and "at least one UMI" in w for w in progress.warnings)
    b = _obs_only(pd.DataFrame({"sample": "T", "condition": "tumor", "score_chondrocyte": np.linspace(0, 1, n),
                                "pos_osteoclast": True}))
    quiet = pc.Progress()
    b = ba.split_osteoclast_compartments(b, 90, progress=quiet)
    assert b.uns["compartment_threshold"] == pytest.approx(np.percentile(np.linspace(0, 1, n), 90))
    assert quiet.warnings == [] and b.obs["growth_plate"].sum() == 20


def test_compartment_summary_counts_compartments_roi_and_depth():
    comp = ["growth_plate"] * 3 + ["tumor_associated"] * 4 + ["negative"] * 3 + ["non_growth_plate_control"] * 2 + \
        ["growth_plate"] + ["negative"] * 7
    a = _obs_only(pd.DataFrame({
        "sample": ["T"] * 10 + ["C"] * 10, "oc_compartment": comp,
        "growth_plate": [True] * 3 + [False] * 7 + [False] * 2 + [True] + [False] * 7,
        "in_marrow_roi": [False, True, True, False, False, True, True, True, True, True] + [True] * 10,
        "total_counts": [100, 200, 300, 10, 20, 30, 40, 1, 1, 1] + [5] * 2 + [700] + [5] * 7}))
    a.uns["compartment_threshold"] = 0.0
    s = ba.compartment_summary(a, samples=["C", "T"]).set_index("sample")
    assert list(s.index) == ["all", "C", "T"] and list(s.reset_index().columns) == ba.COMPARTMENT_SUMMARY_COLUMNS
    assert s.loc["T", "n_growth_plate"] == 3 and s.loc["T", "n_tumor_associated"] == 4
    assert s.loc["T", "n_growth_plate_outside_roi"] == 1 and s.loc["T", "n_tumor_associated_outside_roi"] == 2
    assert s.loc["C", "n_non_growth_plate_control"] == 2 and s.loc["C", "n_growth_plate"] == 1
    assert s.loc["T", "median_umi_growth_plate"] == 200 and s.loc["T", "median_umi_tumor_associated"] == 25
    assert s.loc["all", "n_growth_plate"] == 4 and s.loc["all", "median_umi_growth_plate"] == 250
    assert s.loc["all", "pct_bins_above_threshold"] == 20.0 and s.loc["T", "pct_bins_above_threshold"] == 30.0
    assert bool(s.loc["all", "threshold_is_zero"]) and np.isnan(s.loc["C", "median_umi_tumor_associated"])


# ---------------------------------------------------------------------------------------------- item 7: Moran FDR
@pytest.mark.parametrize("path", ["squidpy", "builtin"])
def test_morans_i_fdr_is_finite_next_to_a_zero_variance_gene(path, monkeypatch):
    if path == "squidpy":
        pytest.importorskip("squidpy")
    else:
        monkeypatch.setitem(sys.modules, "squidpy", None)                       # import squidpy -> ImportError
    sc.settings.verbosity = 0
    a = _toy(n_per=196)
    counts = a.layers["counts"].toarray()
    s1 = (a.obs["sample"] == "S1").to_numpy()
    counts[:, 0] = 0
    counts[s1, 0] = np.where(a.obs["array_col"].to_numpy()[s1] < 7, 40, 0)       # g0: spatial block in S1, absent in S2
    a.layers["counts"] = sp.csr_matrix(counts)
    a.X = a.layers["counts"].copy()
    a = ba.normalize_and_cluster(a, cluster=False)
    moran, _ = ba.morans_i(a, n_svg=40, n_consensus=10, samples=["S1", "S2"], progress=pc.Progress())
    s2 = moran[moran["sample"] == "S2"].set_index("gene")
    assert "g0" in s2.index and np.isnan(s2.loc["g0", "pval_norm"])
    others = s2.drop(index="g0")
    assert np.isfinite(others["pval_norm_fdr_bh"]).all() and np.isnan(s2.loc["g0", "pval_norm_fdr_bh"])
    expected = multipletests(others["pval_norm"].to_numpy(), method="fdr_bh")[1]
    np.testing.assert_allclose(others["pval_norm_fdr_bh"].to_numpy(), expected, rtol=1e-9)


# ---------------------------------------------------------------------------------------------- item 12: unresolved
def test_unresolved_table_has_a_pooled_row_then_one_row_per_section():
    a = _obs_only(pd.DataFrame({
        "sample": ["S1"] * 10 + ["S2"] * 10, "score_x": np.arange(20, dtype=float),
        "leiden": ["unresolved"] * 2 + ["0"] * 8 + ["unresolved"] * 6 + ["1"] * 4,
        "in_marrow_roi": ([True] * 5 + [False] * 5) * 2}))
    _, un = ba.leiden_signature_matrix(a, primary=["x"], sigs={"x": {}}, samples=["S1", "S2"])
    assert list(un.columns) == ba.UNRESOLVED_COLUMNS and list(un["sample"]) == ["all", "S1", "S2"]
    u = un.set_index("sample")
    assert (u.loc["all", "unresolved_bins"], u.loc["all", "all_bins"], u.loc["all", "percent"]) == (8, 20, 40.0)
    assert (u.loc["S1", "percent"], u.loc["S2", "percent"]) == (20.0, 60.0)
    assert (u.loc["S1", "unresolved_bins_roi"], u.loc["S1", "all_bins_roi"], u.loc["S1", "percent_roi"]) == (2, 5, 40.0)
    assert u.loc["S2", "percent_roi"] == 100.0


# ---------------------------------------------------------------------------------------------- item 13: gene sets
def test_geneset_sensitivity_names_used_and_missing_genes():
    a = _toy()
    a.obs["roi"] = True
    sigs = {"alpha": {"core": GENES[:4], "extended": GENES[:6], "competes": [], "excluded": {}},
            "beta": {"core": GENES[10:13], "extended": GENES[10:15], "competes": [], "excluded": {}}}
    progress = pc.Progress()
    t = ba.geneset_sensitivity(a, n_reps=2, variants={"v": {"alpha": ["g0", "g1", "absentX"]}, "table1": {}},
                               names=("alpha",), comp=("beta",), sigs=sigs, roi_key="roi", min_core=0, n_perm=3,
                               n_perm_resampled=3, n_sub=300, progress=progress)
    assert list(t.columns) == ba.GENESET_SENSITIVITY_COLUMNS
    v, t1 = t[t["gene_set"] == "v"].iloc[0], t[t["gene_set"] == "table1"].iloc[0]
    assert (v["genes_used"], v["genes_missing"]) == ("g0,g1", "absentX")
    assert (t1["genes_used"], t1["genes_missing"]) == (",".join(GENES[:6]), "")
    assert any("absentX" in w for w in progress.warnings)


# ---------------------------------------------------------------------------------------------- run on tiny_joint
def test_tiny_joint_tables_reproduce_from_the_object(tiny_joint):
    a, cfg, out = tiny_joint["a"].copy(), tiny_joint["cfg"], tiny_joint["out"]
    nat = pd.read_csv(out / "Figure4g_6a_positivity_native_depth_roi.csv")
    rate = 100 * nat["n_positive"] / nat["n_bins"]                                # as positivity_table computes it
    is_ref = nat["sample"] == cfg["reference_sample"]
    ref = rate[is_ref].set_axis(nat.loc[is_ref, "signature"]).replace(0, np.nan)
    expected = (rate / nat["signature"].map(ref)).round(4)
    np.testing.assert_allclose(nat["fold_vs_reference"], expected, rtol=1e-12, equal_nan=True)
    dm = pd.read_csv(out / "Figure4g_6a_positivity_depth_matched.csv")
    for r in dm.itertuples():
        lo, hi = proportion_confint(int(round(r.n_positive_mean)), int(r.n_bins), method="wilson")
        assert r.ci_lo == pytest.approx(round(100 * lo, 4)) and r.ci_hi == pytest.approx(round(100 * hi, 4))

    cp = cfg["compartments"]
    assert a.uns["compartment_threshold"] == pytest.approx(np.nanpercentile(a.obs[f"score_{cp['split_by']}"],
                                                                            cp["percentile"]))
    summ = pd.read_csv(out / "Figure4h_compartment_summary.csv", dtype={"sample": str}, keep_default_na=False,
                       na_values=["", "NaN", "nan"])
    again = ba.compartment_summary(a, samples=pc.sample_ids(cfg))
    pd.testing.assert_frame_equal(summ, again, check_dtype=False)
    sections = summ[summ["sample"] != "all"]
    for c in ("n_growth_plate", "n_tumor_associated", "n_growth_plate_outside_roi"):
        assert summ.loc[summ["sample"] == "all", c].iloc[0] == sections[c].sum()
    de = pd.read_csv(out / "Figure4h_osteoclast_compartment_DE.csv")
    assert list(de.columns) == ba.COMPARTMENT_DE_COLUMNS
    if len(de):
        chondro = set(cfg["signatures"][cp["split_by"]]["extended"])
        assert set(de.loc[de["selection_gene"].astype(bool), "names"]) == chondro & set(de["names"])
        assert de["log2fc_mean_cp10k"].abs().max() < 19 and de["pct_group"].between(0, 100).all()

    un = pd.read_csv(out / "FigureS5_unresolved_bins.csv", dtype={"sample": str}, keep_default_na=False)
    assert list(un["sample"]) == ["all"] + pc.sample_ids(cfg)
    ld = a.obs["leiden"].astype(str)
    for r in un.itertuples():
        m = np.ones(a.n_obs, bool) if r.sample == "all" else (a.obs["sample"].astype(str) == r.sample).to_numpy()
        assert r.unresolved_bins == int((ld[m] == "unresolved").sum()) and r.all_bins == int(m.sum())
        roi = m & a.obs["in_marrow_roi"].to_numpy(bool)
        assert r.unresolved_bins_roi == int((ld[roi] == "unresolved").sum())
    moran = pd.read_csv(out / "morans_I_per_section.csv")
    assert np.isfinite(moran.loc[moran["pval_norm"].notna(), "pval_norm_fdr_bh"]).all()
    manifest = {f["file"]: f for f in json.loads((out / "manifest.json").read_text())["files"]}
    assert "threshold" in manifest["Figure4h_compartment_summary.csv"]["description"]


def test_module_docstring_no_longer_claims_every_transcript_is_murine():
    assert "every transcript is murine" not in ba.__doc__.split("[FIX]")[0]
    assert "cross-capture" in ba.__doc__ and "Krt7" in ba.__doc__


# ---------------------------------------------------------------------------------------------- captions
@pytest.fixture
def drawn(tmp_path):
    """make_figures on the synthetic results of test_figures (default configuration), with summary.json carrying the
    analysis-bin retention; returns (results folder, captions by figure name)."""
    cfg = pc.default_config()
    folder = make_results(tmp_path / "results", cfg)
    summary = json.loads((folder / "summary.json").read_text())
    for s, pct in zip(summary["samples"], (78.1, 90.4)):
        s["pct_retained"] = pct
    (folder / "summary.json").write_text(json.dumps(summary))
    entries = bf.make_figures(folder, None, pc.Progress())
    assert all(e["error"] is None for e in entries)
    return folder, {e["name"]: e["caption"] for e in entries}


def test_s1_caption_gives_the_analysis_bin_retention(drawn):
    folder, cap = drawn
    s1 = cap["FigureS1_quality_control"]
    sweep = pd.read_csv(folder / "qc_threshold_sweep_16um.csv").query("min_genes == 10").set_index("sample")
    assert "applied to the 16 µm analysis bins" in s1 and "open squares, dashed: the 16 µm analysis bins" in s1
    assert (f"more than 10 genes alone control {sweep.loc['LL-2-1', 'pct_retained']:.1f}%, tumor-bearing "
            f"{sweep.loc['RL-1-5', 'pct_retained']:.1f}%") in s1
    assert "all three criteria control 78.1%, tumor-bearing 90.4%" in s1
    assert "The same gene criterion on the 8 µm bins would keep" in s1


def test_s1_draws_old_results_without_the_analysis_sweep(tmp_path):
    cfg = pc.default_config()
    folder = make_results(tmp_path / "results", cfg)
    (folder / "qc_threshold_sweep_16um.csv").unlink()
    entries = {e["name"]: e for e in bf.make_figures(folder, cfg, pc.Progress())}
    s1 = entries["FigureS1_quality_control"]
    assert s1["error"] is None and "open squares" not in s1["caption"]
    assert "applied to the 16 µm analysis bins" in s1["caption"]


def test_figure3a_caption_names_threshold_zero_signatures_and_removed_shares(drawn):
    folder, cap = drawn
    pb = pd.read_csv(folder / "per_bin_signature_calls_16um.csv.gz", usecols=["sample_id", "tissue_class"])
    share = pb.assign(r=pb["tissue_class"] == "removed").groupby("sample_id")["r"].mean() * 100
    c = cap["Figure3a_1_control_tissue_classes_denoise"]
    assert "Skeletal muscle, smooth muscle and nerve have no positivity call and compete with a threshold of 0" in c
    assert f"{share['LL-2-1']:.1f}% of the bins of this section" in c
    assert f"control {share['LL-2-1']:.1f}%; tumor-bearing {share['RL-1-5']:.1f}%" in c


def test_positivity_captions_explain_the_intervals(drawn):
    _, cap = drawn
    for name in ("Figure4g_bone_lineage_positivity", "Figure6a_myeloid_positivity"):
        c = cap[name]
        assert "treat bins as independent" in c and "understate the uncertainty" in c
        assert "One section per condition" in c and "SD over the downsampling replicates only" in c
        assert "ci_lo and ci_hi columns" in c


def test_figure4h_caption_and_daggers(drawn, tmp_path, monkeypatch):
    _, cap = drawn
    c = cap["Figure4h_osteoclast_compartment_DE"]
    for fact in ("is 0 on these data", "at least one UMI of a chondrocyte gene", "(2.1% of all bins)",
                 "Bins tested: 116 growth-plate against 2,252 tumor-associated", "39 of 116 growth-plate",
                 "595 of 2,252 tumor-associated", "growth-plate 751.5, tumor-associated 349 (ratio 2.15)",
                 "all tested bins come from one section", "descriptive", "†: genes of the chondrocyte signature "
                 "(Acan and Comp)", "Not tested: 31 non-growth-plate"):
        assert fact in c, fact
    assert "growth-plate bins from every section" not in c
    figs = []
    real_save = bf.save
    monkeypatch.setattr(bf, "save", lambda fig, out, name, dpi=200: (figs.append(fig), real_save(fig, out, name, dpi)))
    de = pd.DataFrame({"group": ["growth_plate"] * 3 + ["tumor_associated"] * 2,
                       "names": ["Acan", "Col2a1", "Ibsp", "Ctsk", "Acp5"], "scores": [9.0, 8.0, 7.0, 6.0, 5.0],
                       "logfoldchanges": [25.0, 3.0, 2.0, 1.0, 1.0], "pvals": [1e-9] * 5, "pvals_adj": [1e-8] * 5})
    old = bf.figure4h(de, pc.default_config(), tmp_path)                            # old table, no summary
    assert old[0]["error"] is None and "Effective rule" not in old[0]["caption"]
    labels = [t.get_text() for t in figs[-1].axes[0].get_yticklabels()]
    assert "Acan†" in labels and "Col2a1†" in labels and "Ibsp" in labels             # fallback: Table 1 list
    new = de.drop(columns="logfoldchanges").assign(selection_gene=[True, False, False, False, False])
    bf.figure4h(new, pc.default_config(), tmp_path)
    labels = [t.get_text() for t in figs[-1].axes[0].get_yticklabels()]
    assert "Acan†" in labels and "Col2a1" in labels                                  # the column decides


def _figure4i_caption(folder, cfg):
    """The Figure 4i caption drawn by bone_figures' own preparation step on results `folder`."""
    ctx = bf._Context(folder, pc.merge_with_defaults(cfg), pc.Progress(), None)
    ctx.out.mkdir(parents=True, exist_ok=True)
    [(_, draw)] = bf._prep_4i(ctx)
    [entry] = draw()
    return entry["caption"]


def test_figure4i_caption_states_the_null_model(drawn):
    folder, cap = drawn
    c = cap["Figure4i_osteoblast_osteoclast_distance"]
    assert "independently and uniformly at random" in c and "complete spatial randomness" in c
    assert "random placement only" in c and "single depth-matching draw" in c
    # the step flag is on by default, but the results hold no sensitivity table: the caption must not point to one
    assert pc.default_config()["steps"]["sensitivity_analyses"]
    assert not (folder / bf.SENSITIVITY_4I_CSV).exists()
    assert "Sensitivity_Figure4i" not in c and "sensitivity_Figure4i" not in c
    pd.DataFrame({"sample": ["LL-2-1", "RL-1-5"]}).to_csv(folder / bf.SENSITIVITY_4I_CSV, index=False)
    c = _figure4i_caption(folder, pc.default_config())
    assert "Sensitivity_Figure4i_null_models" in c and "sensitivity_Figure4i_summary.csv" in c
    cfg = pc.default_config()
    cfg["steps"]["sensitivity_analyses"] = False                     # a table left over from another run is not cited
    assert "Sensitivity_Figure4i" not in _figure4i_caption(folder, cfg)


def test_figure_s5_caption_gives_the_share_per_section(drawn, tmp_path):
    folder, cap = drawn
    un = pd.read_csv(folder / "FigureS5_unresolved_bins.csv", dtype={"sample": str}).set_index("sample")
    c = cap["FigureS5_leiden_domain_heatmap"]
    assert f"Unresolved: 1,234 of {int(un.loc['all', 'all_bins']):,} bins ({un.loc['all', 'percent']:.1f}%)." in c
    assert (f"control {un.loc['LL-2-1', 'percent']:.1f}% ({un.loc['LL-2-1', 'percent_roi']:.1f}% of its region of "
            f"interest)") in c
    M = pd.read_csv(folder / "FigureS5_domain_by_signature_zscores.csv")
    old = pd.DataFrame([{"unresolved_bins": 10, "all_bins": 40, "percent": 25.0}])        # single-row table
    e = bf.figureS5(M, pc.default_config(), tmp_path, old)[0]
    assert "Unresolved: 10 of 40 bins (25.0%)." in e["caption"] and "per section" not in e["caption"]
