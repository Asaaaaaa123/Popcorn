"""Shared pytest fixtures: a tiny synthetic Visium HD dataset and a cheap configuration that points at it."""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pipeline_config as pc  # noqa: E402

CHEAP = {
    "scoring": {"n_perm": 5, "n_perm_resampled": 5, "n_depth_reps": 2, "null_n_sub": 2000},
    "clustering": {"stability_seeds": [0, 1]},
    "moran": {"n_genes": 200},
    "uncoupling": {"n_perm": 19},
    "steps": {"sensitivity": False, "project_fine_bins": True},
}


@pytest.fixture(scope="session")
def tiny_demo(tmp_path_factory):
    """{"root": Path, "samples": [...]}: two tiny synthetic sections (CTRL-1 control, TUM-1 tumor-bearing)."""
    from tools.make_demo_data import generate
    root = tmp_path_factory.mktemp("demo")
    return {"root": root, "samples": generate(root, size="tiny", seed=0)}


@pytest.fixture
def tiny_config(tiny_demo):
    """A full configuration (defaults + cheap settings) analysing the tiny demo sections."""
    samples = tiny_demo["samples"]
    return pc.merge_with_defaults({**CHEAP, "data_root": str(tiny_demo["root"]), "samples": samples,
                                   "reference_sample": samples[0]["id"]})


@pytest.fixture(scope="session")
def tiny_joint(tmp_path_factory, tiny_demo):
    """{"a", "cfg", "out"}: the joint AnnData of the whole analysis of the tiny demo sections (cheap settings, no
    sensitivity analyses, no figures), run once per session. Tests must not modify `a` in place: copy it first."""
    import bone_analysis as ba
    samples = tiny_demo["samples"]
    cfg = pc.merge_with_defaults({**CHEAP, "data_root": str(tiny_demo["root"]), "samples": samples,
                                  "reference_sample": samples[0]["id"],
                                  "steps": {**CHEAP["steps"], "project_fine_bins": False, "sensitivity_analyses": False}})
    out = tmp_path_factory.mktemp("tiny_joint")
    a = ba.run_analysis_from_config(cfg, out, pc.Progress())
    return {"a": a, "cfg": pc.load_config(out / "config.json"), "out": out}


@pytest.fixture(scope="session")
def original_ba():
    """The untouched manuscript version of bone_analysis.py (_original/), imported under another name."""
    spec = importlib.util.spec_from_file_location("original_bone_analysis", ROOT / "_original" / "bone_analysis.py")
    module = importlib.util.module_from_spec(spec)
    saved = dict(os.environ)
    spec.loader.exec_module(module)
    # its import-time block sets NUMBA_NUM_THREADS after numba has started (defect A15): undo it for the other tests
    os.environ.clear()
    os.environ.update(saved)
    return module
