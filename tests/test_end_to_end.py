"""The whole analysis and every figure on the tiny synthetic dataset, checked against docs/ARCHITECTURE.md section 3."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

import bone_analysis as ba
import bone_figures as bf
import pipeline_config as pc

pytestmark = pytest.mark.slow

COLUMNS = {
    "qc_threshold_sweep_16um.csv": ["sample", "min_genes", "n_retained", "pct_retained", "median_genes"],
    "qc_threshold_sweep_8um.csv": ["sample", "min_genes", "n_retained", "pct_retained", "median_genes", "bin_size"],
    "FigureS1_qc_per_bin_8um.csv.gz": ["sample", "bin_size", "n_genes_by_counts", "total_counts", "pct_counts_mt"],
    "Figure6b_receptor_panel_depth_matched.csv": ["gene", "role", "sample", "roi", "mean_raw", "pct_detected_raw",
                                                  "mean_dm", "pct_detected_dm", "n_positive_dm", "n_bins", "ci_lo_dm",
                                                  "ci_hi_dm", "median_umi_per_bin", "fold_dm"],
    "Figure4g_6a_positivity_native_depth_roi.csv": ["signature", "sample", "roi", "n_positive", "n_bins",
                                                    "pct_positive", "ci_lo", "ci_hi", "mean_score",
                                                    "fold_vs_reference", "fold_tumor_over_control"],
    "Figure4g_6a_positivity_depth_matched.csv": ["signature", "sample", "pct_positive_mean", "pct_positive_sd",
                                                 "n_positive_mean", "n_bins", "mean_score", "ci_lo", "ci_hi"],
    "leiden_domain_screen.csv": ["domain", "n_bins", "ubiquitous_fraction", "median_genes_per_umi", "top_markers",
                                 "de_uninterpretable"],
    "Figure4h_osteoclast_compartment_DE.csv": ["group", "names", "scores", "pvals", "pvals_adj", "log2fc_mean_cp10k",
                                               "pct_group", "pct_rest", "selection_gene"],
    "Figure4h_compartment_summary.csv": ["sample", "threshold", "threshold_is_zero", "pct_bins_above_threshold",
                                         "n_growth_plate", "n_tumor_associated", "n_non_growth_plate_control",
                                         "n_growth_plate_outside_roi", "n_tumor_associated_outside_roi",
                                         "median_umi_growth_plate", "median_umi_tumor_associated"],
    "Figure3a_program_composition_native.csv": ["sample", "program", "n_bins", "pct_of_all_bins",
                                                "pct_of_assigned_bins", "sample_ambiguous_pct"],
    "leiden_seed_stability_adjusted_rand_index.csv": ["seed_a", "seed_b", "n_clusters_a", "n_clusters_b", "ARI"],
    "FigureS5_unresolved_bins.csv": ["sample", "unresolved_bins", "all_bins", "percent", "unresolved_bins_roi",
                                     "all_bins_roi", "percent_roi"],
    "Figure4i_observed_to_null_summary.csv": ba.UNCOUPLING_COLUMNS,
    "Figure4i_distance_distributions_observed_and_null.csv": ["sample", "source", "distance_um", "count"],
    "fibroblast_to_osteoblast_proximity.csv": ["section", "fibroblast_assigned_bins", "osteoblast_assigned_bins",
                                               "region_of_interest_bins", "percent_fibroblast_within_50um_of_osteoblast",
                                               "percent_all_bins_within_50um_of_osteoblast"],
    "Table1_gene_signatures.csv": ["signature", "canonical_markers", "signature_genes",
                                   "genes_counted_toward_two_gene_rule", "excluded_genes"],
}


class RecordingProgress(pc.Progress):
    """Progress that remembers every label it reports."""

    def __init__(self):
        super().__init__()
        self.reported = []

    def step(self, label):
        self.reported.append(label)
        super().step(label)

    def skip(self, label, reason):
        self.reported.append(label)
        super().skip(label, reason)


@pytest.fixture(scope="module")
def run(tmp_path_factory, tiny_demo):
    from tests.conftest import CHEAP
    samples = tiny_demo["samples"]
    cfg = pc.merge_with_defaults({**CHEAP, "data_root": str(tiny_demo["root"]), "samples": samples,
                                  "reference_sample": samples[0]["id"]})
    out = tmp_path_factory.mktemp("results")
    progress = RecordingProgress()
    a = ba.run_analysis_from_config(cfg, out, progress)
    figs = bf.make_figures(out, cfg, progress)
    return {"cfg": cfg, "out": out, "a": a, "progress": progress, "figures": figs}


def test_every_documented_table_exists_with_its_columns(run):
    for name, cols in COLUMNS.items():
        df = pd.read_csv(run["out"] / name, keep_default_na=False) if "distance" in name else pd.read_csv(run["out"] / name)
        assert list(df.columns) == cols, name
        assert len(df), name
    for name in ("config.json", "summary.json", "manifest.json", "joint_16um.h5ad", "labels_8um_projected.csv",
                 "morans_I_per_section.csv", "FigureS2_spatially_variable_genes_morans_I.csv",
                 "FigureS4_adipocyte_genes_depth_matched.csv", "leiden_domain_markers.csv",
                 "macrophage_to_osteoclast_proximity.csv", "Figure3b_program_composition_depth_matched.csv",
                 "FigureS5_domain_by_signature_zscores.csv"):
        assert (run["out"] / name).is_file(), name


def test_manifest_lists_every_written_file(run):
    manifest = json.loads((run["out"] / "manifest.json").read_text())["files"]
    listed = {f["file"] for f in manifest}
    written = {p.name for p in run["out"].iterdir() if p.is_file() and p.name != "manifest.json"}
    assert written <= listed
    assert all(f["title"] and f["description"] and f["kind"] in ("table", "object", "data") for f in manifest)


def test_summary_shape(run):
    s = json.loads((run["out"] / "summary.json").read_text())
    assert s["version"] == 1 and s["reference_sample"] == run["cfg"]["reference_sample"]
    assert [x["id"] for x in s["samples"]] == pc.sample_ids(run["cfg"])
    assert {"Sections", "Positivity, native depth (ROI)", "Positivity, depth-matched (ROI)"} <= {m["group"] for m in s["metrics"]}
    assert all(set(m) == {"group", "label", "value", "unit", "detail"} for m in s["metrics"])
    assert s["environment"]["scanpy"] and s["timings"]


def test_reported_steps_follow_the_plan(run):
    assert run["progress"].reported[:len(ba.planned_steps(run["cfg"]))] == ba.planned_steps(run["cfg"])
    assert run["progress"].reported[len(ba.planned_steps(run["cfg"])):] == bf.planned_steps(run["cfg"])


def test_per_bin_table(run):
    pb = pd.read_csv(run["out"] / "per_bin_signature_calls_16um.csv.gz")
    assert set(pb["sample_id"]) == set(pc.sample_ids(run["cfg"]))
    assert {"positive_granulocyte", "positive_osteoclast_pooled_depthmatched"} <= set(pb.columns)
    assert set(pb["tissue_class"]) <= {"bone lineage", "chondrocyte", "hematopoietic marrow", "stroma", "removed",
                                       "unassigned"}
    assert pb["in_region_of_interest"].dtype == bool and len(pb) == run["a"].n_obs


def test_pooled_and_depth_matched_scores_are_recomputed(run):
    obs = run["a"].obs
    for n in ("osteoclast", "osteoblast", "macrophage"):
        assert not np.allclose(obs[f"pooled_score_{n}"], obs[f"score_{n}"])
        assert not np.allclose(obs[f"pooled_dm_score_{n}"], obs[f"pooled_score_{n}"])
    assert run["a"].raw is None


SENSITIVITY_FIGURE_MODULES = ("positivity_figures", "proximity_figures", "uncoupling_figures", "compartments_figures",
                              "clustering_figures")


def test_sensitivity_analyses_ran(run):
    import importlib
    names = {f["name"] for f in run["figures"]}
    for m in SENSITIVITY_FIGURE_MODULES:
        mod = importlib.import_module(f"bone_sensitivity.{m}")
        assert (run["out"] / mod.INPUTS[0]).is_file(), m
        for name, _ in mod.FIGURES:
            assert name in names, name
    s = json.loads((run["out"] / "summary.json").read_text())
    assert any(m["group"] == "Sensitivity analyses" for m in s["metrics"])
    assert not any("sensitivity analysis" in w and "failed" in w for w in s["warnings"])
    pb = pd.read_csv(run["out"] / "per_bin_signature_calls_16um.csv.gz", nrows=5)
    assert {"positive_fraction_osteoclast_pooled_depthmatched_draws", "in_lesion"} <= set(pb.columns)


def test_sensitivity_analyses_change_no_manuscript_output(run, tiny_joint):
    """The same analysis with the sensitivity step off (tiny_joint) writes identical manuscript tables."""
    skip = {"summary.json", "manifest.json", "config.json", "per_bin_signature_calls_16um.csv.gz"}
    compared = 0
    for path in sorted(tiny_joint["out"].glob("*.csv*")):
        if path.name in skip or path.name.startswith("sensitivity_"):
            continue
        other = run["out"] / path.name
        assert other.is_file(), path.name
        pd.testing.assert_frame_equal(pd.read_csv(path), pd.read_csv(other), check_exact=False, rtol=1e-9,
                                      obj=path.name)
        compared += 1
    assert compared > 15
    off = pd.read_csv(tiny_joint["out"] / "per_bin_signature_calls_16um.csv.gz")
    on = pd.read_csv(run["out"] / "per_bin_signature_calls_16um.csv.gz")
    pd.testing.assert_frame_equal(off, on[off.columns])


def test_figures_all_drawn(run):
    figs = run["figures"]
    assert figs and all(f["error"] is None for f in figs)
    fig_dir = run["out"] / "figures"
    for f in figs:
        assert (fig_dir / f["png"]).stat().st_size > 1000 and (fig_dir / f["pdf"]).is_file()
    names = {f["name"] for f in figs}
    for expected in ("Figure3b_composition_denoise", "Figure4a-f_bone_lineage_maps", "Figure4g_bone_lineage_positivity",
                     "Figure6a_myeloid_positivity", "Figure4h_osteoclast_compartment_DE",
                     "Figure4i_osteoblast_osteoclast_distance", "Figure6b_receptor_panel", "FigureS1_quality_control",
                     "FigureS2_spatially_variable_genes", "FigureS4_adipocyte_genes", "FigureS5_leiden_domain_heatmap"):
        assert expected in names, expected
    assert sum(n.startswith("Figure3a_") for n in names) == len(run["cfg"]["samples"])
