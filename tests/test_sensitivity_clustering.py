"""bone_sensitivity.clustering (Supplementary Figure S5 sensitivity): synthetic checks of the statistics, a run on the
tiny demo analysis that reproduces the manuscript S5 outputs, and its figure."""
from __future__ import annotations

import copy

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from bone_sensitivity import clustering as cl
from bone_sensitivity import clustering_figures as clf


# ---- synthetic -------------------------------------------------------------------------------------------------

def _toy_obs():
    """Section A: 10 bins, 2 unresolved (1 in the ROI), 1 unassigned; section B: 20 bins, 10 unresolved (8 in the
    ROI), lesion = 6 of B's unresolved ROI bins + 2 resolved B bins."""
    leiden = ["unresolved"] * 2 + ["unassigned"] + ["0"] * 7 + ["unresolved"] * 10 + ["1"] * 10
    roi = [True, False] + [True] * 8 + [True] * 8 + [False] * 2 + [True] * 5 + [False] * 5
    sample = ["A"] * 10 + ["B"] * 20
    lesion = np.zeros(30, dtype=bool)
    lesion[10:16] = True
    lesion[20:22] = True
    return pd.DataFrame({"leiden": leiden, "sample": sample, "in_marrow_roi": roi}), lesion


def test_unresolved_by_section_known_shares():
    obs, lesion = _toy_obs()
    df = cl.unresolved_by_section(obs, ["A", "B"], lesion).set_index(["sample", "region"])
    assert list(df.columns[:4]) == ["n_bins", "n_unresolved", "pct_unresolved", "n_unassigned"]
    assert df.loc[("A", "all"), ["n_bins", "n_unresolved", "n_unassigned"]].tolist() == [10, 2, 1]
    assert df.loc[("A", "all"), "pct_unresolved"] == 20.0
    assert df.loc[("A", "roi"), ["n_bins", "n_unresolved"]].tolist() == [9, 1]
    assert df.loc[("A", "outside_roi"), "pct_unresolved"] == 100.0
    assert df.loc[("B", "all"), "pct_unresolved"] == 50.0
    assert df.loc[("B", "roi"), ["n_bins", "n_unresolved"]].tolist() == [13, 8]
    assert df.loc[("all", "all"), ["n_bins", "n_unresolved", "n_unassigned"]].tolist() == [30, 12, 1]
    assert df.loc[("all", "all"), "pct_unresolved"] == 40.0
    # lesion: B has 8 lesion bins, 6 of them unresolved; they are 6 of B's 10 unresolved bins
    assert df.loc[("B", "all"), "n_lesion"] == 8
    assert df.loc[("B", "all"), "pct_lesion_unresolved"] == 75.0
    assert df.loc[("B", "all"), "pct_unresolved_lesion"] == 60.0
    assert np.isnan(df.loc[("A", "all"), "pct_lesion_unresolved"])        # no lesion bin in A
    assert df.loc[("all", "all"), "pct_unresolved_lesion"] == 50.0         # 6 of 12
    no_les = cl.unresolved_by_section(obs.drop(columns="in_marrow_roi"), ["A", "B"])
    assert set(no_les["region"]) == {"all"} and "n_lesion" not in no_les.columns


def test_domain_order():
    assert cl.domain_order(["unresolved", "10", "2", "unassigned", "0", "x"]) == ["0", "2", "10", "x", "unresolved",
                                                                                  "unassigned"]


def test_marker_summary_and_flags_on_a_toy_marker_table():
    m = pd.DataFrame({
        "group": ["A"] * 4 + ["B"] * 3 + ["C"] * 3 + ["D"] * 2,
        "names": ["a1", "a2", "a3", "a4", "b1", "b2", "b3", "c1", "c2", "c3", "d1", "d2"],
        "scores": [9, 8, 7, 1, 1.5, 1.0, -3, 12, 0.5, -1, -5, -6],
        "pvals_adj": [1e-9, 1e-6, 1e-3, .5, .2, .6, 1e-5, 1e-12, .3, .9, 1e-8, 1e-7]})
    s = cl.summarize_markers(m, 0.05)
    assert s.loc["A", "n_significant_markers"] == 3 and s.loc["A", "top_marker"] == "a1"
    assert s.loc["B", "n_significant_markers"] == 0 and s.loc["B", "top_marker"] == ""     # no marker: none named
    assert s.loc["C", "n_significant_markers"] == 1 and s.loc["C", "top_marker"] == "c1"
    assert s.loc["D", "n_significant_markers"] == 0                     # significant but negative: not a marker
    assert s.loc["D", "top_marker"] == ""                               # d1 is depleted, not a marker
    assert cl.summarize_markers(m, 1e-4).loc["A", "n_significant_markers"] == 2
    flags = cl.flag_domains([3, 0, 1, 1, 2, 1], [95, 99, 95, 80, 90, 95], [2, 1, 2, 2, 5, 6])
    assert list(flags) == ["", "no significant marker", "defined by one gene", "", "defined by one gene", ""]
    assert list(cl.flag_domains([np.nan], [np.nan], [np.nan])) == [""]


def test_top_programs_text():
    prog = ["osteoblast"] * 23 + ["osteoclast"] * 20 + ["skeletal_muscle"] * 18 + ["nerve"] * 39 + ["ambiguous"] * 50
    assert cl.top_programs(prog) == "nerve 39%; osteoblast 23%; osteoclast 20%"
    assert cl.top_programs(["ambiguous"] * 3) == ""


def test_partition_section_stats_known_answer():
    # cluster 0: 9 of S1 + 1 of S2 (specific to S1); cluster 1: 5 + 5; cluster 2: 10 of S2 (specific to S2)
    labels = ["0"] * 10 + ["1"] * 10 + ["2"] * 10
    sections = ["S1"] * 9 + ["S2"] + ["S1"] * 5 + ["S2"] * 5 + ["S2"] * 10
    i = np.array([0, 1, 2, 9, 20])
    j = np.array([1, 2, 9, 20, 21])
    g = sp.coo_matrix((np.ones(5), (i, j)), shape=(30, 30)).tocsr()
    g = g + g.T + sp.identity(30, format="csr")                         # self-loops are not edges
    st = cl.partition_section_stats(labels, sections, g, ["S1", "S2"])
    assert st["n_clusters"] == 3
    assert st["n_section_specific_clusters"] == 2
    assert st["pct_bins_in_section_specific_clusters"] == round(100 * 20 / 30, 2)
    assert st["pct_S1_in_section_specific_clusters"] == round(100 * 9 / 14, 2)
    assert st["pct_S2_in_section_specific_clusters"] == round(100 * 10 / 16, 2)
    assert st["same_section_edge_fraction"] == 0.8                      # 0-1, 1-2, 9-20, 20-21 same; 2-9 differs
    assert st["same_section_edge_fraction_random"] == round((14 / 30) ** 2 + (16 / 30) ** 2, 4)


def _two_section_embedding(offset, n=150, seed=0):
    """Three Gaussian groups seen in two sections; `offset` shifts section S2 along the first component."""
    rng = np.random.default_rng(seed)
    centres = np.eye(3, 8) * 8
    grp = np.tile(np.arange(3), 2 * n // 3)[:2 * n]
    sec = np.repeat(["S1", "S2"], n)
    X = centres[grp] + rng.normal(size=(2 * n, 8))
    X[sec == "S2", 3] += offset
    return X, sec, grp


def test_leiden_partition_identical_embeddings_and_harmony_check():
    import anndata as ad
    from sklearn.metrics import adjusted_rand_score
    X, sec, grp = _two_section_embedding(0.0)
    cfg = {"seed": 0, "samples": [{"id": "S1"}, {"id": "S2"}],
           "clustering": {"n_neighbors": 15, "n_pcs": 8, "resolution": 0.3, "leiden_n_iterations": 2}}
    lab1, g1 = cl.leiden_partition(X, 15, 8, 0.3, 2, 0)
    lab2, _ = cl.leiden_partition(X.copy(), 15, 8, 0.3, 2, 0)
    assert (lab1 == lab2).all() and g1.shape == (len(X), len(X))
    assert adjusted_rand_score(lab1, grp) > 0.9                          # the three groups are recovered
    a = ad.AnnData(obs=pd.DataFrame({"sample": sec, "leiden": lab1}, index=[f"b{i}" for i in range(len(X))]),
                   obsm={"X_pca": X, "X_pca_harmony": X.copy()}, uns={"neighbors_rep": "X_pca_harmony"})
    df = cl.harmony_check(a, cfg).set_index("variant")
    assert list(df.index) == ["harmony", "no_harmony"]
    assert df.loc["harmony", "ari_vs_harmony"] == 1.0 and df.loc["no_harmony", "ari_vs_harmony"] == 1.0
    assert df.loc["harmony", "reproduces_manuscript_domains"] is True
    assert df.loc["harmony", "n_section_specific_clusters"] == 0
    # a section offset in the uncorrected embedding only: the PCA partition splits by section
    Xs, _, _ = _two_section_embedding(12.0)
    b = ad.AnnData(obs=pd.DataFrame({"sample": sec}, index=a.obs_names), obsm={"X_pca": Xs, "X_pca_harmony": X})
    df = cl.harmony_check(b, cfg).set_index("variant")
    assert df.loc["no_harmony", "same_section_edge_fraction"] > 0.99
    assert df.loc["harmony", "same_section_edge_fraction"] < 0.7
    assert df.loc["no_harmony", "pct_bins_in_section_specific_clusters"] == 100.0
    assert df.loc["no_harmony", "ari_vs_harmony"] < 0.9
    assert pd.isna(df.loc["no_harmony", "reproduces_manuscript_domains"])
    only_pca = ad.AnnData(obs=b.obs.copy(), obsm={"X_pca": Xs})
    df = cl.harmony_check(only_pca, cfg, progress=None)
    assert list(df["variant"]) == ["no_harmony"] and df["ari_vs_harmony"].isna().all()


# ---- tiny demo analysis -----------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def run_out(tiny_joint, tmp_path_factory):
    """clustering.run on a copy of the tiny joint object (with a synthetic lesion), tables written to a folder."""
    a = tiny_joint["a"].copy()
    tum = a.obs["condition"].astype(str).to_numpy() == "tumor"
    rng = np.random.default_rng(0)
    a.obs["in_lesion"] = tum & ((a.obs["leiden"].astype(str).to_numpy() == "unresolved") | (rng.random(a.n_obs) < .1))
    before = (list(a.obs.columns), sorted(a.uns.keys()), sorted(a.layers.keys()), a.obs.copy())
    outs = cl.run(a, tiny_joint["cfg"])
    res = tmp_path_factory.mktemp("clustering_results")
    for o in outs:
        o.df.to_csv(res / o.file, index=False)
    return {"a": a, "outs": {o.file: o for o in outs}, "res": res, "before": before}


def test_run_leaves_the_object_unchanged(run_out):
    a, (cols, uns, layers, obs) = run_out["a"], run_out["before"]
    assert list(a.obs.columns) == cols and sorted(a.uns.keys()) == uns and sorted(a.layers.keys()) == layers
    pd.testing.assert_frame_equal(a.obs, obs)


def test_run_outputs_and_columns(run_out, tiny_joint):
    outs, samples = run_out["outs"], [s["id"] for s in tiny_joint["cfg"]["samples"]]
    assert set(outs) == {cl.UNRESOLVED_CSV, cl.ANNOTATION_CSV, cl.HARMONY_CSV}
    for o in outs.values():
        assert o.title and o.description and len(o.df)
    assert list(outs[cl.UNRESOLVED_CSV].df.columns) == ["sample", "region", "n_bins", "n_unresolved", "pct_unresolved",
                                                        "n_unassigned", "n_lesion", "pct_lesion_unresolved",
                                                        "pct_unresolved_lesion"]
    assert list(outs[cl.ANNOTATION_CSV].df.columns) == (
        ["domain", "manuscript_name", "n_bins"] + [f"pct_{s}" for s in samples]
        + ["median_umi", "median_genes_per_umi", "top_programs", "ambiguous_pct", "n_significant_markers", "top_marker",
           "pct_domain_detecting_top_marker", "pct_rest_detecting_top_marker", "flag", "suggested_label"])
    hm = outs[cl.HARMONY_CSV].df
    for c in ["variant", "n_clusters", "ari_vs_harmony", "same_section_edge_fraction", "n_section_specific_clusters",
              "pct_bins_in_section_specific_clusters"] + [f"pct_{s}_in_section_specific_clusters" for s in samples]:
        assert c in hm.columns
    assert "raw partitions" in outs[cl.HARMONY_CSV].description.lower()


def test_unresolved_table_reproduces_the_manuscript_share(run_out, tiny_joint):
    un = run_out["outs"][cl.UNRESOLVED_CSV].df.set_index(["sample", "region"])
    ref = pd.read_csv(tiny_joint["out"] / "FigureS5_unresolved_bins.csv").iloc[0]
    assert un.loc[("all", "all"), "n_unresolved"] == ref["unresolved_bins"]
    assert un.loc[("all", "all"), "n_bins"] == ref["all_bins"]
    assert round(un.loc[("all", "all"), "pct_unresolved"], 1) == ref["percent"]
    a = run_out["a"]
    for s in a.obs["sample"].astype(str).unique():
        m = a.obs["sample"].astype(str) == s
        assert un.loc[(s, "all"), "n_unresolved"] == int((a.obs.loc[m, "leiden"].astype(str) == "unresolved").sum())
        assert un.loc[(s, "roi"), "n_bins"] + un.loc[(s, "outside_roi"), "n_bins"] == un.loc[(s, "all"), "n_bins"]
    les = a.obs["in_lesion"].to_numpy(bool)
    unres = a.obs["leiden"].astype(str).to_numpy() == "unresolved"
    assert un.loc[("all", "all"), "n_lesion"] == les.sum()
    assert un.loc[("all", "all"), "pct_lesion_unresolved"] == round(100 * (les & unres).sum() / les.sum(), 2)


def test_annotation_reproduces_the_manuscript_names_and_markers(run_out, tiny_joint):
    ann = run_out["outs"][cl.ANNOTATION_CSV].df.set_index("domain")
    a, samples = run_out["a"], [s["id"] for s in tiny_joint["cfg"]["samples"]]
    M = pd.read_csv(tiny_joint["out"] / "FigureS5_domain_by_signature_zscores.csv")
    assert sorted(ann["manuscript_name"]) == sorted(M["leiden_domain"].astype(str))
    for d, name in ann["manuscript_name"].items():                      # each name sits on its own domain
        assert name == d or name.endswith(f"(domain {d})")
    assert ann["n_bins"].sum() == a.n_obs
    assert np.allclose(ann[[f"pct_{s}" for s in samples]].sum(axis=1), 100, atol=0.02)
    med = a.obs.groupby(a.obs["leiden"].astype(str))["total_counts"].median()
    assert np.allclose(ann["median_umi"], med.reindex(ann.index))
    markers = pd.read_csv(tiny_joint["out"] / "leiden_domain_markers.csv")
    n_top = tiny_joint["cfg"]["clustering"]["n_top_markers"]
    for d, r in ann.iterrows():
        sub = markers[markers["group"].astype(str) == d]
        sig = sub[(sub["pvals_adj"] < 0.05) & (sub["scores"] > 0)]
        assert len(sig) == min(r["n_significant_markers"], n_top)        # significant markers: a prefix by score
        assert r["flag"] in ("", "no significant marker", "defined by one gene")
        if not r["n_significant_markers"]:                              # no marker: no top gene, no detection
            assert r["top_marker"] == "" and r["flag"] == "no significant marker"
            assert np.isnan(r["pct_domain_detecting_top_marker"]) and np.isnan(r["pct_rest_detecting_top_marker"])
            assert "domain" in r["suggested_label"] or d in ("unresolved", "unassigned")
            continue
        top = r["top_marker"]
        assert top == sub.sort_values("scores", ascending=False)["names"].iloc[0]
        det = np.asarray(a.layers["counts"][:, a.var_names.get_loc(top)].todense()).ravel() > 0
        inside = a.obs["leiden"].astype(str).to_numpy() == d
        assert abs(r["pct_domain_detecting_top_marker"] - 100 * det[inside].mean()) <= 0.006
        assert abs(r["pct_rest_detecting_top_marker"] - 100 * det[~inside].mean()) <= 0.006
        assert r["suggested_label"]


def test_harmony_check_reproduces_the_manuscript_partition(run_out, tiny_joint):
    hm = run_out["outs"][cl.HARMONY_CSV].df.set_index("variant")
    screen = pd.read_csv(tiny_joint["out"] / "leiden_domain_screen.csv")
    assert hm.loc["harmony", "reproduces_manuscript_domains"] is True
    assert hm.loc["harmony", "n_clusters"] == len(screen)                # raw domains before merging
    assert hm.loc["harmony", "ari_vs_harmony"] == 1.0
    assert 0 <= hm.loc["no_harmony", "ari_vs_harmony"] <= 1
    assert 0 <= hm.loc["no_harmony", "same_section_edge_fraction"] <= 1


def test_harmony_check_can_be_switched_off(tiny_joint):
    cfg = copy.deepcopy(tiny_joint["cfg"])
    cfg.setdefault("sensitivity_analyses", {}).setdefault("clustering", {})["harmony_check"] = False
    outs = cl.run(tiny_joint["a"].copy(), cfg)
    assert [o.file for o in outs] == [cl.UNRESOLVED_CSV, cl.ANNOTATION_CSV]
    assert "n_lesion" not in outs[0].df.columns                          # the demo has no Krt7/Krt18


def test_lesion_fallback_matches_a_recomputation(tiny_joint):
    """Without obs in_lesion the lesion comes from the shared definition (two demo genes renamed to the lesion genes)
    and is not stored."""
    import bone_analysis as ba
    from bone_sensitivity import common
    a = tiny_joint["a"].copy()
    assert "in_lesion" not in a.obs
    cfg = copy.deepcopy(tiny_joint["cfg"])
    cfg.setdefault("sensitivity_analyses", {}).setdefault("clustering", {})["harmony_check"] = False
    les = common.settings(cfg, None)["lesion"]
    assert not set(les["genes"]) & set(a.var_names)
    X = a.layers["counts"].tocsc()
    A = ba.spatial_neighbor_matrix(a, int(les["pool_k"]))
    tot = np.asarray(X.sum(axis=0)).ravel()
    order = [j for j in np.argsort(tot, kind="stable") if tot[j] > 0]
    flag = None
    for i, j in zip(order[:-1], order[1:]):                              # rarest pair flagging some, not all, bins
        pooled = np.asarray(A @ np.asarray(X[:, [i, j]].sum(axis=1)).ravel()).ravel()
        f = pooled >= float(les["min_pooled_umi"])
        if 0 < f.sum() < a.n_obs:
            flag, pair = f, (i, j)
            break
    assert flag is not None
    names = a.var_names.to_numpy().astype(object)
    names[list(pair)] = les["genes"][:2]
    a.var_names = pd.Index(names)
    outs = {o.file: o for o in cl.run(a, cfg)}
    assert "in_lesion" not in a.obs                                     # computed, not stored
    o = outs[cl.UNRESOLVED_CSV]
    assert "+".join(les["genes"][:2]) in o.description
    un = o.df.set_index(["sample", "region"])
    unres = a.obs["leiden"].astype(str).to_numpy() == "unresolved"
    smp = a.obs["sample"].astype(str).to_numpy()
    for s in list(pd.unique(smp)) + ["all"]:
        m = np.ones(a.n_obs, bool) if s == "all" else smp == s
        fl, both = flag & m, flag & m & unres
        assert un.loc[(s, "all"), "n_lesion"] == fl.sum()
        exp = round(100 * both.sum() / fl.sum(), 2) if fl.sum() else np.nan
        assert un.loc[(s, "all"), "pct_lesion_unresolved"] == pytest.approx(exp, nan_ok=True)
        exp = round(100 * both.sum() / (m & unres).sum(), 2) if (m & unres).sum() else np.nan
        assert un.loc[(s, "all"), "pct_unresolved_lesion"] == pytest.approx(exp, nan_ok=True)


# ---- figure -----------------------------------------------------------------------------------------------------

def _check_entries(entries, out):
    assert [e["name"] for e in entries] == [n for n, _ in clf.FIGURES]
    for e in entries:
        assert set(e) >= {"name", "group", "title", "caption", "png", "pdf", "error"}
        assert e["error"] is None and e["caption"]
        assert (out / e["png"]).stat().st_size > 0 and (out / e["pdf"]).stat().st_size > 0


def test_draw_from_run_outputs(run_out, tiny_joint, tmp_path):
    assert set(clf.INPUTS) == set(run_out["outs"])
    entries = clf.draw(run_out["res"], tiny_joint["cfg"], tmp_path / "fig")
    _check_entries(entries, tmp_path / "fig")
    cap = entries[0]["caption"]
    pooled = run_out["outs"][cl.UNRESOLVED_CSV].df.query("sample == 'all' and region == 'all'")["pct_unresolved"]
    assert f"pooled {pooled.iloc[0]:.1f}%" in cap
    ari = run_out["outs"][cl.HARMONY_CSV].df.set_index("variant").loc["no_harmony", "ari_vs_harmony"]
    assert f"adjusted Rand index {ari:.2f}" in cap and "lesion bins" in cap
    assert "differ only by seed" not in cap                              # no seed table in this folder
    # with the manuscript seed-to-seed table next to them, the caption quotes it as the reference of the ARI
    res = tmp_path / "res"
    res.mkdir()
    for f in clf.INPUTS + clf.OPTIONAL_INPUTS:
        src = (run_out["res"] if f in clf.INPUTS else tiny_joint["out"]) / f
        (res / f).write_bytes(src.read_bytes())
    seeds = pd.read_csv(res / cl.STABILITY_CSV)
    cap = clf.draw(res, tiny_joint["cfg"], tmp_path / "fig_seeds")[0]["caption"]
    assert f"median adjusted Rand index of {seeds['ARI'].median():.2f}" in cap and cl.STABILITY_CSV in cap


def test_seed_reference_caption():
    seeds = pd.DataFrame({"ARI": [0.70, 0.72, 0.80]})
    assert "more than changing the seed" in clf._caption_seed_reference(0.44, seeds)
    assert "within the seed-to-seed range" in clf._caption_seed_reference(0.75, seeds)
    assert "less than changing the seed" in clf._caption_seed_reference(0.95, seeds)
    assert "median adjusted Rand index of 0.72 (range 0.70-0.80, 3 pairs)" in clf._caption_seed_reference(0.44, seeds)
    assert clf._caption_seed_reference(0.44, None) == "" and clf._caption_seed_reference(0.44, seeds.iloc[:0]) == ""


def test_draw_with_missing_tables(run_out, tiny_joint, tmp_path):
    res = tmp_path / "res"
    res.mkdir()
    (res / cl.UNRESOLVED_CSV).write_bytes((run_out["res"] / cl.UNRESOLVED_CSV).read_bytes())
    _check_entries(clf.draw(res, tiny_joint["cfg"], tmp_path / "fig1"), tmp_path / "fig1")
    (res / cl.UNRESOLVED_CSV).unlink()
    ann = run_out["outs"][cl.ANNOTATION_CSV].df.copy()
    ann.loc[ann.index[0], ["flag", "n_significant_markers"]] = ["no significant marker", 0]
    ann.loc[ann.index[1], "top_programs"] = ""
    ann.to_csv(res / cl.ANNOTATION_CSV, index=False)
    entries = clf.draw(res, tiny_joint["cfg"], tmp_path / "fig2")
    _check_entries(entries, tmp_path / "fig2")
    assert "no significant marker" in entries[0]["caption"]
    (res / cl.ANNOTATION_CSV).unlink()
    with pytest.raises(FileNotFoundError):
        clf.draw(res, tiny_joint["cfg"], tmp_path / "fig3")
