"""
tests/test_other_signatures.py
==============================
A lab that does not study osteoclasts and osteoblasts: the bone signatures deleted in the Signatures tab (S1), or
every signature replaced by other ones (S2), then the settings changed as the configuration check says (S3: the
bone-specific steps switched off, or pointed at the new signatures). Validation must name every setting that still
names a deleted signature and the ways out; a run must finish without a crash or a failed figure, skipping what has
nothing to compute with a reason. The slow tests run the whole analysis and the figures on the tiny synthetic data.

delete_signature() and add_signature() follow deleteSignature() and the Add button of
webapp/static/js/components/signature-editor.js (there is no JavaScript test runner).
"""
from __future__ import annotations

import copy
import json
import shutil
import time
import warnings

import anndata as ad
import pandas as pd
import pytest

import bone_analysis as ba
import bone_figures as bf
import pipeline_config as pc
from bone_sensitivity import uncoupling as U

with warnings.catch_warnings():  # Starlette 1.x prefers httpx2 for its TestClient; httpx still works
    warnings.filterwarnings("ignore", message=".*httpx.*")
    from fastapi.testclient import TestClient

from webapp.server import create_app  # noqa: E402

NEW_COLORS = ["#1F77B4", "#2CA02C", "#D62728", "#9467BD", "#8C564B"]
# name: (class, core genes, scored genes, competes); every gene is in the tiny synthetic probe set
NONBONE = {
    "myeloid": ("immune", ["Lyz2", "C1qa", "C1qb", "Cd68", "Adgre1"],
                ["Lyz2", "Apoe", "C1qa", "C1qb", "C1qc", "Cd68", "Adgre1", "Mrc1", "Trem2", "Ctss"], ["neutro"]),
    "neutro": ("immune", ["S100a8", "S100a9", "Ngp", "Camp", "Ltf"],
               ["S100a8", "S100a9", "Ngp", "Camp", "Ltf", "Lcn2", "Retnlg", "Mpo", "Elane", "Ly6g"], []),
    "endo": ("vascular", ["Pecam1", "Cdh5", "Emcn"], ["Pecam1", "Cdh5", "Emcn", "Flt1", "Kdr", "Cldn5"], []),
    "fibro": ("stroma", ["Dcn", "Lum", "Pdgfra"], ["Col3a1", "Dcn", "Lum", "Postn", "Fbn1", "Pdgfra", "Pdgfrb", "Thy1"],
              []),
    "rbc": ("blood", ["Hba-a1", "Hbb-bs", "Alas2"],
            ["Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2", "Slc4a1", "Gypa", "Klf1"], []),
}
CLASSES = [{"name": "immune", "color": "#FF7F00", "emphasis": True},
           {"name": "vascular", "color": "#17BECF", "emphasis": False},
           {"name": "stroma", "color": "#E6AB02", "emphasis": False},
           {"name": "blood", "color": "#A65628", "emphasis": False}]
# cheap but not degenerate (the demo project's scoring settings; few sensitivity resamples)
SETTINGS = {
    "scoring": {"n_perm": 20, "n_perm_resampled": 10, "n_depth_reps": 2, "null_n_sub": 5000},
    "clustering": {"stability_seeds": [0, 1]},
    "uncoupling": {"n_perm": 49},
    "steps": {"moran": False, "write_h5ad": False, "project_fine_bins": False},
    "sensitivity_analyses": {"n_boot": 50, "uncoupling": {"n_draws": 2, "n_perm": 19},
                             "compartments": {"n_depth_reps": 2}},
}


# ============================================================================
# THE SIGNATURES TAB, IN PYTHON
# ============================================================================

def delete_signature(cfg, name):
    """deleteSignature(): the signature and every list entry that names it; proximity pairs that use it are dropped.
    Single-signature settings (compartments, uncoupling, figures.lineage_maps) are left for the user."""
    if name not in cfg["signatures"]:
        return
    del cfg["signatures"][name]
    drop = lambda lst: [x for x in lst if x != name]  # noqa: E731
    for s in cfg["signatures"].values():
        s["competes"] = drop(s.get("competes") or [])
    cfg["primary"], cfg["programs"] = drop(cfg["primary"]), drop(cfg["programs"])
    cfg["roi"]["exclude"] = drop(cfg["roi"]["exclude"])
    cfg["proximity"]["pairs"] = [p for p in cfg["proximity"]["pairs"] if name not in (p["from"], p["to"])]
    sens = cfg["sensitivity"]
    sens["signatures"], sens["competitors"] = drop(sens["signatures"]), drop(sens["competitors"])
    for sets in sens["variants"].values():
        sets.pop(name, None)
    for p in cfg["figures"]["positivity_panels"]:
        p["signatures"] = drop(p["signatures"])
    cfg["figures"]["composition_signatures"] = drop(cfg["figures"]["composition_signatures"])


def add_signature(cfg, name, klass, core, extended, competes=()):
    """Add (class 'removed', empty lists), then the editor: genes, class, competitors and both switches on."""
    used = {str(s.get("color")).lower() for s in cfg["signatures"].values()}
    color = next((c for c in NEW_COLORS if c.lower() not in used), NEW_COLORS[0])
    cfg["signatures"][name] = {"label": name.replace("_", " "), "table_label": name.replace("_", " "), "color": color,
                               "class": klass, "canonical": list(core[:2]), "core": list(core),
                               "extended": list(extended), "competes": list(competes), "excluded": {}}
    for key in ("primary", "programs"):
        if name not in cfg[key]:
            cfg[key].append(name)


def s1(cfg):
    """S1: the osteoclast and osteoblast signatures deleted."""
    delete_signature(cfg, "osteoclast")
    delete_signature(cfg, "osteoblast")
    return cfg


def s2(cfg):
    """S2: every signature deleted, five others added (all with a positivity call and in the assignment), tissue
    classes replaced."""
    for name in list(cfg["signatures"]):
        delete_signature(cfg, name)
    for name, (klass, core, ext, comp) in NONBONE.items():
        add_signature(cfg, name, klass, core, ext, comp)
    cfg["classes"] = copy.deepcopy(CLASSES)
    return cfg


def s3_repoint(cfg):
    """S3, the first way out: every setting the check names pointed at the new signatures, every step kept."""
    s2(cfg)
    cfg["compartments"].update(signature="myeloid", split_by="endo")
    cfg["uncoupling"].update(a="fibro", b="endo")
    cfg["figures"]["lineage_maps"] = {"a": "myeloid", "b": "endo"}
    cfg["figures"]["positivity_panels"] = [{"name": "Figure_positivity_1", "signatures": list(NONBONE)}]
    cfg["figures"]["composition_signatures"] = list(NONBONE)
    cfg["proximity"]["pairs"] = [{"from": "fibro", "to": "endo", "mode": "assignment"},
                                 {"from": "myeloid", "to": "neutro", "mode": "positivity"}]
    cfg["sensitivity_analyses"]["uncoupling"]["habitat_classes"] = ["vascular"]
    cfg["sensitivity_analyses"]["lesion"]["genes"] = []          # no A549 lesion here
    return cfg


def s3_off(cfg):
    """S3, the second way out: the steps the errors name switched off, nothing else changed."""
    s2(cfg)
    cfg["steps"].update(compartments=False, uncoupling=False)
    return cfg


def starting_with(messages, prefix):
    found = [m for m in messages if m.startswith(prefix)]
    assert len(found) == 1, (prefix, messages)
    return found[0]


# ============================================================================
# VALIDATION (no data)
# ============================================================================

def test_the_default_configuration_reports_nothing():
    assert pc.validate_config(pc.default_config()) == ([], [])


def test_s1_names_each_setting_its_field_and_both_ways_out():
    errors, warns = pc.validate_config(s1(pc.default_config()))
    assert len(errors) == 3, errors
    for prefix, field, step in (
            ("compartments.signature = 'osteoclast'", "Osteoclast compartments (Figure 4h) > Signature split into "
                                                      "compartments", "Osteoclast compartments (Figure 4h)"),
            ("uncoupling.a = 'osteoblast'", "Distance to nearest partner (Figure 4i) > From (signature)",
             "Distance to nearest partner (Figure 4i)"),
            ("uncoupling.b = 'osteoclast'", "Distance to nearest partner (Figure 4i) > To (signature)",
             "Distance to nearest partner (Figure 4i)")):
        e = starting_with(errors, prefix)
        assert "is not a signature of this project" in e and f"Parameters > {field}" in e
        assert f'switch off the step "{step}" (Parameters > Analysis steps)' in e
    for k in ("a", "b"):
        assert "Figure 4a-f is not drawn" in starting_with(warns, f"figures.lineage_maps.{k} = ")
    assert "Parameters > Proximity pairs" in starting_with(warns, "proximity.pairs is empty")
    assert not any(w.startswith(("figures.composition", "positivity panel")) for w in warns)  # others are left

    # way out 1: the steps switched off; nothing about them is reported any more
    off = s1(pc.default_config())
    off["steps"].update(compartments=False, uncoupling=False)
    errors, warns = pc.validate_config(off)
    assert errors == [] and not any(w.startswith(("compartments.", "uncoupling.")) for w in warns)
    # way out 2: other signatures, which need a positivity call
    other = s1(pc.default_config())
    other["compartments"]["signature"] = "macrophage"
    other["uncoupling"].update(a="fibroblast", b="macrophage")
    other["figures"]["lineage_maps"] = {"a": "macrophage", "b": "fibroblast"}
    assert pc.validate_config(other)[0] == []
    other["uncoupling"]["a"] = "skeletal_muscle"
    e = starting_with(pc.validate_config(other)[0], "uncoupling.a = 'skeletal_muscle'")
    assert "has no positivity call" in e and "Signatures tab" in e and "switch off the step" in e


def test_s2_reports_every_setting_left_on_the_bone_signatures():
    errors, warns = pc.validate_config(s2(pc.default_config()))
    assert sorted(e.split(" = ")[0] for e in errors) == ["compartments.signature", "compartments.split_by",
                                                         "uncoupling.a", "uncoupling.b"]
    assert "Figure 3b is not drawn" in starting_with(warns, "figures.composition_signatures lists no signature")
    for panel in ("Figure4g_bone_lineage_positivity", "Figure6a_myeloid_positivity"):
        assert "Positivity panels" in starting_with(warns, f"positivity panel {panel!r} lists no signature")
    starting_with(warns, "figures.lineage_maps.a = 'osteoclast'")
    starting_with(warns, "proximity.pairs is empty")
    assert "habitat_fixed_b" in starting_with(warns, "sensitivity_analyses.uncoupling.habitat_classes (bone lineage)")


def test_s3_either_way_out_leaves_a_clean_configuration():
    assert pc.validate_config(s3_repoint(pc.default_config())) == ([], [])
    errors, warns = pc.validate_config(s3_off(pc.default_config()))
    assert errors == [] and not any(w.startswith(("compartments.", "uncoupling.", "sensitivity_analyses.uncoupling"))
                                    for w in warns)
    # without the bone signatures, Krt7 / Krt18 would flag ordinary epithelium as A549 lesion
    w = starting_with(warns, "sensitivity_analyses.lesion.genes (Krt7, Krt18)")
    assert "A549 lesion flag" in w and "empty the list" in w


def test_settings_that_lose_figures_or_tests_are_warnings():
    cfg = s3_repoint(pc.default_config())
    cfg["sensitivity_analyses"]["uncoupling"]["habitat_classes"] = ["bone lineage"]
    assert "not estimable" in starting_with(pc.validate_config(cfg)[1], "sensitivity_analyses.uncoupling.habitat")

    cfg = s3_repoint(pc.default_config())
    cfg["samples"][0]["condition"], cfg["samples"][1]["condition"] = "untreated", "treated"
    w = starting_with(pc.validate_config(cfg)[1], "compartments.tumor_conditions (tumor) lists the condition of no")
    assert "Tumor conditions" in w and "switch the step off" in w
    cfg["compartments"]["tumor_conditions"] = ["treated"]
    assert pc.validate_config(cfg) == ([], [])

    cfg = s3_repoint(pc.default_config())
    cfg["steps"]["sensitivity"] = True                         # the osteoclast / osteoblast three-gene-set check
    assert "switch it off" in starting_with(pc.validate_config(cfg)[1], 'the step "Three-gene-set sensitivity')
    cfg["figures"]["positivity_panels"][0]["signatures"] = ["myeloid", "ghost"]
    cfg["figures"]["composition_signatures"] = ["myeloid", "ghost"]
    warns = pc.validate_config(cfg)[1]
    assert not any(w.startswith("positivity panel") for w in warns)           # one usable signature is enough
    starting_with(warns, "figures.composition_signatures: 'ghost' not a signature")


def test_a_receptor_panel_of_other_genes_needs_the_figure_genes_too():
    cfg = s3_repoint(pc.default_config())
    cfg["receptor_panel"] = {"Il6ra": "IL-6R", "Socs3": "STAT3 target"}
    w = starting_with(pc.validate_config(cfg)[1], "figures.receptor_genes (Tnfsf11, ")
    assert "Figure 6b is not drawn" in w and "Parameters > Figures > Receptor genes shown in Figure 6b" in w
    cfg["figures"]["receptor_genes"] = ["Socs3", "Il6ra"]
    assert pc.validate_config(cfg) == ([], [])
    cfg["figures"]["receptor_genes"] = []
    cfg["steps"]["receptor_panel"] = False                     # no Figure 6b anyway
    assert pc.validate_config(cfg) == ([], [])
    assert "Distance to nearest partner (Figure 4i)" in pc.STEP_LABELS.values()   # one name for the step


def test_other_compartments_are_not_called_growth_plate():
    assert pc.compartment_names(pc.default_config()) == ("growth-plate", "tumor-associated")
    assert pc.compartment_names(s3_repoint(pc.default_config())) == ("high-endo", "low-endo")


def test_a_signature_that_cannot_be_scored_is_not_drawn_as_zero(tmp_path):
    """A signature with fewer than two of its genes in the data has no score (mean_score NaN): a positivity panel
    leaves it out and says so, instead of a 0% bar with a Wilson interval; a panel of such signatures is skipped."""
    cfg = s3_repoint(pc.default_config())
    rows = []
    for sig, score in (("myeloid", 0.4), ("endo", float("nan"))):
        for smp, n_pos in (("LL-2-1", 3), ("RL-1-5", 9)):
            rows.append({"signature": sig, "sample": smp, "roi": "in_marrow_roi", "n_positive": 0 if score != score
                         else n_pos, "n_bins": 500, "pct_positive": 0.0, "ci_lo": 0.0, "ci_hi": 0.7,
                         "mean_score": score})
    results = tmp_path / "results"
    results.mkdir()
    pd.DataFrame(rows).to_csv(results / bf.NATIVE_CSV, index=False)
    cfg["figures"]["positivity_panels"] = [{"name": "Figure_positivity_1", "signatures": ["myeloid", "endo"]},
                                           {"name": "Figure_positivity_2", "signatures": ["endo"]}]
    progress = Recorder()
    entries = {e["name"]: e for e in bf.make_figures(results, cfg, progress)}
    caption = entries["Figure_positivity_1"]["caption"]
    assert entries["Figure_positivity_1"]["error"] is None
    assert "fewer than two of its genes are in the data, so it cannot be scored: endo." in caption
    skipped = {x["step"]: x["reason"] for x in progress.skipped}
    reason = next(r for step, r in skipped.items() if "positivity 2" in step)
    assert "no signature of positivity panel Figure_positivity_2 can be scored" in reason


def test_the_analysis_refuses_what_validation_refuses(tmp_path):
    with pytest.raises(ValueError, match=r"compartments\.signature = 'osteoclast' is not a signature"):
        ba.run_analysis_from_config(s1(pc.default_config()), tmp_path / "out")


def test_step_switches_explain_what_they_are_for():
    steps = pc.PARAM_SCHEMA[-1]
    assert steps["key"] == "steps" and [f["key"] for f in steps["fields"]] == list(pc.STEP_LABELS)
    assert all(f["type"] == "bool" and len(f["help"]) > 40 and f["label"] == pc.STEP_LABELS[f["key"]]
               for f in steps["fields"])
    assert {f["key"] for f in steps["fields"] if f.get("study") == "bone"} == set(pc.BONE_STUDY_STEPS)
    fields = [(sec["key"], f) for sec in pc.PARAM_SCHEMA[:-1] for f in sec["fields"]]
    assert all(f["uses_step"] in pc.STEP_LABELS for _, f in fields if "uses_step" in f)
    assert ("compartments", "tumor_conditions") in {(s, f["key"]) for s, f in fields}
    titles = {sec["key"]: sec["title"] for sec in pc.PARAM_SCHEMA}
    for key in ("compartments", "uncoupling"):                 # the messages name the cards the form shows
        assert titles[key] in " ".join(pc.validate_config(s1(pc.default_config()))[0])


def test_figure4i_sensitivity_skips_a_deleted_signature():
    cfg = s1(pc.default_config())
    a = ad.AnnData(obs=pd.DataFrame({"sample": ["LL-2-1"] * 3, "pos_osteocyte": [True, False, True]},
                                    index=["b1", "b2", "b3"]))
    progress = pc.Progress()
    assert U.run(a, cfg, progress) == []
    assert any("Figure 4i sensitivity needs positivity calls" in w and "skipped" in w for w in progress.warnings)
    assert list(U.scoring_signatures(cfg, ["osteoblast", "osteocyte"])) == ["osteocyte"]


def test_unsaved_settings_and_gene_lists_are_checked(tmp_path, tiny_demo):
    """POST /validate with data false: the configuration check of unsaved settings (the project page's chips); the
    data check also lists the ubiquitous transcripts and the A549 lesion genes missing from the probe set."""
    with TestClient(create_app(tmp_path / "ws", worker_module="tests.fake_worker"),
                    base_url="http://127.0.0.1:8765") as client:
        samples = [dict(s) for s in tiny_demo["samples"]]
        p = client.post("/api/projects", json={"name": "draft", "config": {
            "data_root": str(tiny_demo["root"]), "samples": samples, "reference_sample": samples[0]["id"]}}).json()
        url = f"/api/projects/{p['id']}/validate"
        draft = s1(copy.deepcopy(p["config"]))
        r = client.post(url, json={"config": draft, "data": False}).json()
        assert set(r) == {"ok", "errors", "warnings"} and r["ok"] is False
        assert (r["errors"], r["warnings"]) == pc.validate_config(pc.merge_with_defaults(draft))
        assert client.post(url, json={"data": False}).json() == {"ok": True, "errors": [], "warnings": []}
        assert client.post(url, json={"data": "no"}).status_code == 400

        cfg = copy.deepcopy(p["config"])
        cfg["ubiquitous"] = ["ACTB", "GAPDH"]                  # human symbols on mouse data
        rep = client.post(url, json={"config": cfg}).json()
        assert rep["genes"]["ubiquitous"]["missing"] == ["ACTB", "GAPDH"]
        assert "the Leiden screen marks no domain unresolved" in starting_with(
            rep["warnings"], "ubiquitous transcripts: none of the 2 genes")
        assert rep["genes"]["lesion_genes"]["missing"] == ["Krt7", "Krt18"]   # not in the tiny probe set
        starting_with(rep["warnings"], "A549 lesion flag: none of its genes (Krt7, Krt18)")
        cfg["sensitivity_analyses"]["lesion"]["genes"] = []
        cfg["ubiquitous"] = ["Actb", "ACTB"]
        rep = client.post(url, json={"config": cfg}).json()
        assert not any(w.startswith("A549") for w in rep["warnings"])
        assert "1 of 2 genes missing from the probe set (ACTB)" in starting_with(rep["warnings"], "ubiquitous")


def test_the_app_starts_a_run_only_after_one_way_out(tmp_path):
    with TestClient(create_app(tmp_path / "ws", worker_module="tests.fake_worker"),
                    base_url="http://127.0.0.1:8765") as client:
        p = client.post("/api/projects", json={"name": "no bone", "config": {}}).json()
        cfg = s1(p["config"])
        put = client.put(f"/api/projects/{p['id']}", json={"config": cfg}).json()
        assert any(e.startswith("compartments.signature") for e in put["validation"]["errors"])
        r = client.post(f"/api/projects/{p['id']}/runs")
        assert r.status_code == 409 and "uncoupling.a = 'osteoblast'" in r.json()["detail"]
        cfg["steps"].update(compartments=False, uncoupling=False)
        assert client.put(f"/api/projects/{p['id']}", json={"config": cfg}).json()["validation"]["errors"] == []
        rid = client.post(f"/api/projects/{p['id']}/runs").json()["id"]
        deadline = time.monotonic() + 60
        while client.get(f"/api/runs/{rid}").json()["status"] in ("queued", "running"):
            assert time.monotonic() < deadline
            time.sleep(0.05)
        assert client.get(f"/api/runs/{rid}").json()["status"] == "succeeded"


# ============================================================================
# WHOLE RUNS ON THE TINY DATA
# ============================================================================

class Recorder(pc.Progress):
    """Progress that keeps the label of every reported step (run or skipped)."""

    def __init__(self):
        super().__init__()
        self.reported = []

    def step(self, label):
        self.reported.append(label)
        super().step(label)

    def skip(self, label, reason):
        self.reported.append(label)
        super().skip(label, reason)


def _run(tmp_path_factory, tiny_demo, edit, name):
    samples = tiny_demo["samples"]
    cfg = edit(pc.merge_with_defaults({**copy.deepcopy(SETTINGS), "data_root": str(tiny_demo["root"]),
                                       "samples": [dict(s) for s in samples], "reference_sample": samples[0]["id"]}))
    validation = pc.validate_config(cfg)
    out = tmp_path_factory.mktemp(name)
    progress = Recorder()
    ba.run_analysis_from_config(cfg, out, progress)
    figures = bf.make_figures(out, cfg, progress)
    progress.done()
    return {"cfg": cfg, "validation": validation, "out": out, "progress": progress, "figures": figures,
            "skipped": {s["step"]: s["reason"] for s in progress.skipped}}


def _no_failure(run):
    assert run["progress"].reported == ba.planned_steps(run["cfg"]) + bf.planned_steps(run["cfg"])
    assert all(f["error"] is None for f in run["figures"]), [f for f in run["figures"] if f["error"]]
    bad = [w for w in run["progress"].warnings if "failed" in w or "could not be drawn" in w or "Error" in w]
    assert bad == []
    summary = json.loads((run["out"] / "summary.json").read_text(encoding="utf-8"))
    assert summary["skipped"] == run["progress"].skipped[:len(summary["skipped"])]


def _skip_reason(run, start):
    found = [r for label, r in run["skipped"].items() if label.startswith(start)]
    assert len(found) == 1, (start, run["skipped"])
    return found[0]


@pytest.fixture(scope="module")
def s1_run(tmp_path_factory, tiny_demo):
    """S1 with the way out that switches the steps off; the three-gene-set check switched on (nothing to test)."""
    def edit(cfg):
        s1(cfg)
        cfg["steps"].update(compartments=False, uncoupling=False, sensitivity=True)
        return cfg
    return _run(tmp_path_factory, tiny_demo, edit, "s1")


@pytest.fixture(scope="module")
def s3_run(tmp_path_factory, tiny_demo):
    """S2 with every setting pointed at the new signatures (S3) and every step on, sensitivity analyses included."""
    return _run(tmp_path_factory, tiny_demo, s3_repoint, "s3")


@pytest.mark.slow
def test_s1_run_finishes_and_skips_with_reasons(s1_run):
    errors, warns = s1_run["validation"]
    assert errors == [] and any(w.startswith("figures.lineage_maps.a") for w in warns)
    _no_failure(s1_run)
    assert _skip_reason(s1_run, "Osteoclast compartments") == "disabled in the configuration"
    assert _skip_reason(s1_run, "Distance to nearest partner") == "disabled in the configuration"
    assert "no proximity pairs" in _skip_reason(s1_run, "Proximity tables")
    assert "sensitivity.signatures names no signature" in _skip_reason(s1_run, "Gene-set sensitivity")
    assert "figures.lineage_maps.a = 'osteoclast' is not a signature" in _skip_reason(s1_run, "Draw Figures 4a-f")
    assert _skip_reason(s1_run, "Draw Sensitivity — Figure 4i") == \
        "step disabled: Distance to nearest partner (Figure 4i)"
    assert _skip_reason(s1_run, "Draw Sensitivity — Figure 4h").startswith("step disabled: Osteoclast compartments")
    assert "no proximity pairs" in _skip_reason(s1_run, "Draw Sensitivity — proximity")
    names = {f["name"] for f in s1_run["figures"]}
    assert {"Figure3b_composition_denoise", "Figure4g_bone_lineage_positivity", "Figure6a_myeloid_positivity",
            "Figure6b_receptor_panel", "FigureS5_leiden_domain_heatmap",
            "Sensitivity_positivity_folds"} <= names
    out = s1_run["out"]
    assert not (out / "Figure4h_osteoclast_compartment_DE.csv").exists()
    assert not (out / "sensitivity_Figure4i_null_models.csv").exists()
    assert (out / "sensitivity_positivity_folds.csv").is_file()
    pb = pd.read_csv(out / "per_bin_signature_calls_16um.csv.gz", nrows=5)
    assert "positive_osteocyte" in pb.columns and "positive_osteoclast" not in pb.columns


@pytest.mark.slow
def test_s3_run_draws_every_figure_for_the_new_signatures(s3_run):
    assert s3_run["validation"] == ([], [])
    _no_failure(s3_run)
    assert not any("lesion" in w for w in s3_run["progress"].warnings)        # the gene list was emptied
    # Figure 4h of a myeloid split by the endo score: no "growth plate" or "tumor-associated" in its text
    words = [w for w in s3_run["progress"].warnings if "compartments" in w and "Figure 4h" not in w]
    assert all("high-endo" in w and "growth" not in w for w in words if w.startswith("myeloid compartments"))
    for f in s3_run["figures"]:
        if f["name"] in ("Figure4h_osteoclast_compartment_DE", "Sensitivity_Figure4h_variants"):
            caption = f["caption"].lower().replace("growth_plate_masks", "")       # the setting keeps its name
            assert "high-endo" in caption and "growth" not in caption, f["caption"]
            assert "tumor-associated" not in f["caption"]
    reason = s3_run["skipped"].get("Draw Sensitivity — Figure 4h under other definitions and tests")
    if reason is not None and "has no rows" in reason:                           # the test was not estimable
        assert "Figure 4h was not estimable, so there is no gene to draw" in reason
    out = s3_run["out"]
    for name in ("Figure4g_6a_positivity_native_depth_roi.csv", "Figure3b_program_composition_depth_matched.csv",
                 "Figure4h_compartment_summary.csv", "Figure4i_observed_to_null_summary.csv",
                 "fibro_to_endo_proximity.csv", "myeloid_to_neutro_proximity.csv", "sensitivity_proximity.csv",
                 "sensitivity_Figure4i_null_models.csv", "sensitivity_Figure4h_variants.csv"):
        assert (out / name).is_file(), name
    native = pd.read_csv(out / "Figure4g_6a_positivity_native_depth_roi.csv")
    assert set(native["signature"]) == set(NONBONE)
    pb = pd.read_csv(out / "per_bin_signature_calls_16um.csv.gz")
    assert {f"positive_{n}" for n in NONBONE} <= set(pb.columns)
    assert set(pb["tissue_class"]) <= {c["name"] for c in CLASSES} | {"unassigned"}
    nulls = pd.read_csv(out / "sensitivity_Figure4i_null_models.csv")
    hab = nulls[nulls["null_model"] == "habitat_fixed_b"]
    assert len(hab) and (hab["n_habitat_extra"].fillna(0) > 0).any()           # the 'vascular' habitat exists
    names = {f["name"] for f in s3_run["figures"]}
    assert {"Figure3b_composition_denoise", "Figure4a-f_bone_lineage_maps", "Figure_positivity_1",
            "Figure4i_osteoblast_osteoclast_distance", "Figure6b_receptor_panel",
            "Sensitivity_positivity_folds", "Sensitivity_proximity_baselines",
            "Sensitivity_Figure4i_null_models"} <= names
    captions = {f["name"]: f["caption"] for f in s3_run["figures"]}
    assert "myeloid" in captions["Figure_positivity_1"] and "fibro" in captions["Figure4i_osteoblast_osteoclast_distance"]
    # every planned figure step drew its figure or says why not (e.g. a compartment with too few bins)
    assert all(reason for reason in s3_run["skipped"].values())


@pytest.mark.slow
def test_figures_with_emptied_lists_are_skipped_not_failed(s3_run, tmp_path):
    """S3 results drawn with the figure lists S2 leaves behind: nothing fails, each skip names the setting."""
    results = tmp_path / "results"
    shutil.copytree(s3_run["out"], results)
    shutil.rmtree(results / "figures")
    cfg = copy.deepcopy(s3_run["cfg"])
    cfg["figures"].update(composition_signatures=[], lineage_maps={"a": "osteoclast", "b": "osteoblast"},
                          positivity_panels=[{"name": "Figure4g_bone_lineage_positivity", "signatures": []}])
    cfg["receptor_panel"] = {"Socs3": "STAT3 target", "Il6st": "gp130"}       # a panel of other genes
    progress = Recorder()
    entries = bf.make_figures(results, cfg, progress)
    assert all(e["error"] is None for e in entries) and progress.warnings == []
    skipped = {s["step"]: s["reason"] for s in progress.skipped}
    assert "composition bars in Parameters > Figures" in skipped["Draw Figure 3b — tissue composition"]
    assert "Positivity panels" in skipped["Draw Figure 4g — bone-lineage positivity"]
    assert "figures.lineage_maps.a = 'osteoclast'" in skipped["Draw Figures 4a-f — osteoclast and osteoblast maps"]
    assert "Receptor genes shown in Figure 6b" in skipped["Draw Figure 6b — host receptor panel"]
    assert progress.reported == bf.planned_steps(cfg)
