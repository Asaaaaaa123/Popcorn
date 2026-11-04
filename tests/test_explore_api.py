"""
tests/test_explore_api.py
=========================
Data explorer (POST /api/projects/{pid}/explore/...) and the gene maps of the run viewer (/api/runs/{rid}/genes,
/gene, /bins/image) on the tiny synthetic dataset. The server reads the raw sections with h5py and scipy; these tests
prove with bone_analysis (scanpy) that it reads them exactly as the analysis does: the same bins in the same order, the
same coordinates, per-bin totals and quality-control decisions, and that a gene map of a run holds the raw counts of the
run's barcodes.
"""
from __future__ import annotations

import copy
import gzip
import json
import os
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image

import pipeline_config as pc
from tests.test_api import FAKE, LOCAL, TestClient, client, demo_config, new_project, start_run, wait_for  # noqa: F401
from webapp import server
from webapp.server import create_app

ROOT = Path(__file__).resolve().parents[1]


# ============================================================================
# HELPERS
# ============================================================================

def bin_dir(sample: dict, bin_size: str = "016um") -> Path:
    return Path(sample["path"]) / f"square_{bin_size}"


def scalefactors(sample: dict, bin_size: str = "016um") -> dict:
    return json.loads((bin_dir(sample, bin_size) / "spatial" / "scalefactors_json.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def raw(tiny_demo):
    """Sample id -> the section as the analysis loads it (bone_analysis.read_visium_hd + qc_metrics), 16 um bins."""
    import bone_analysis as ba
    return {s["id"]: ba.qc_metrics(ba.read_visium_hd(bin_dir(s), s["id"], s["condition"], s["microns_per_pixel"]))
            for s in tiny_demo["samples"]}


def explore(client, pid: str, sid: str, what: str, **body):
    return client.post(f"/api/projects/{pid}/explore/{sid}/{what}", json=body)


def gene_counts(a, genes: list[str]) -> np.ndarray:
    """Per-bin sum of the raw counts of `genes` in an AnnData (int64)."""
    return np.rint(np.asarray(a[:, genes].X.sum(axis=1)).ravel()).astype(np.int64)


def dense(vector: dict) -> np.ndarray:
    out = np.zeros(vector["n"], np.int64)
    out[vector["indices"]] = vector["values"]
    return out


def copy_sample(tiny_demo: dict, dest: Path, index: int = 0) -> Path:
    """A copy of a demo section (its binned_outputs folder) at dest; returns the copied binned_outputs folder."""
    target = dest / "outs" / "binned_outputs"
    shutil.copytree(tiny_demo["samples"][index]["path"], target)
    return target


# ============================================================================
# READER UNITS
# ============================================================================

def test_names_are_made_unique_like_anndata():
    from anndata.utils import make_index_unique
    for names in (["a", "b"], ["a", "a", "b", "a"], ["a", "a-1", "a", "a"], ["x-1", "x", "x", "x-2", "x"], []):
        assert server.make_unique(names) == list(make_index_unique(pd.Index(names, dtype=object))), names
    assert server.norm_barcodes([" s_016um_00001_00002-1 ", "s_016um_00001_00002", "AAAC-12", "a-b"]).tolist() == [
        "s_016um_00001_00002", "s_016um_00001_00002", "AAAC", "a-b"]


def test_cache_builds_a_key_once():
    import threading
    import time
    cache, calls = server.LRUCache(2), []

    def build():
        calls.append(1)
        time.sleep(0.2)
        return object()

    results = []
    threads = [threading.Thread(target=lambda: results.append(cache.get("k", build))) for _ in range(5)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert len(calls) == 1 and len(results) == 5 and all(r is results[0] for r in results)
    for key in ("a", "b"):
        cache.get(key, object)
    cache.get("k", build)  # evicted: built again
    assert len(calls) == 2 and cache._gates == {}


def test_cache_shares_a_failed_build_and_keeps_to_its_budget():
    import threading
    import time
    cache, calls = server.LRUCache(4), []

    def failing():
        calls.append(1)
        time.sleep(0.2)
        raise server.DataError("unreadable")

    errors = []

    def ask():
        try:
            cache.get("k", failing)
        except server.DataError as exc:
            errors.append(str(exc))

    threads = [threading.Thread(target=ask) for _ in range(5)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    # the requests that waited for the failing build get its error instead of building again
    assert len(calls) == 1 and errors == ["unreadable"] * 5 and cache._gates == {}
    assert cache.get("k", lambda: "ok") == "ok"  # the error is not kept: a later request tries again

    sized = server.LRUCache(10, budget=100, weigh=len)
    for key in ("a", "b", "c"):
        sized.get(key, lambda: "x" * 40)
    assert list(sized._items) == ["b", "c"] and sized._weight == 80  # a third value of 40 would exceed 100
    sized.get("big", lambda: "x" * 500)  # the newest value is kept whatever its weight
    assert list(sized._items) == ["big"] and sized._weight == 500
    sized.discard(lambda key: key == "big")
    assert not sized._items and sized._weight == 0


# ============================================================================
# DATA EXPLORER
# ============================================================================

def test_explore_samples(client, tiny_demo):
    p = new_project(client, **demo_config(tiny_demo))
    url = f"/api/projects/{p['id']}/explore/samples"
    r = client.post(url).json()
    assert r["bin_size"] == "016um" and [s["id"] for s in r["samples"]] == ["CTRL-1", "TUM-1"]
    for entry, sample in zip(r["samples"], tiny_demo["samples"]):
        positions = pd.read_parquet(bin_dir(sample) / "spatial" / "tissue_positions.parquet")
        n = int((positions["in_tissue"] == 1).sum())
        assert entry == {"id": sample["id"], "label": sample["label"], "condition": sample["condition"],
                         "color": sample["color"], "n_bins": n, "has_image": True, "ok": True,
                         "message": f"{n:,} bins at 016um"}

    # unsaved edits are explored instead of the saved configuration
    cfg = p["config"]
    cfg["bin_size"] = "008um"
    cfg["samples"][1]["path"] = "nowhere"
    r = client.post(url, json={"config": cfg}).json()
    ctrl, tum = r["samples"]
    positions = pd.read_parquet(bin_dir(tiny_demo["samples"][0], "008um") / "spatial" / "tissue_positions.parquet")
    assert r["bin_size"] == "008um" and ctrl["n_bins"] == int((positions["in_tissue"] == 1).sum())
    assert tum["n_bins"] is None and tum["ok"] is False and "square_008um" in tum["message"]
    assert client.get(f"/api/projects/{p['id']}").json()["config"]["bin_size"] == "016um"

    cfg["bin_size"] = "16um"
    assert client.post(url, json={"config": cfg}).status_code == 400
    assert client.post(url, json={"config": [1]}).status_code == 400
    assert client.post("/api/projects/p-00000000/explore/samples").status_code == 404


def test_explore_bins_are_the_bins_of_the_analysis(client, tiny_demo, raw):
    import bone_analysis as ba
    p = new_project(client, **demo_config(tiny_demo))
    n_failing = 0
    for sample in tiny_demo["samples"]:
        sid, a = sample["id"], raw[sample["id"]]
        r = explore(client, p["id"], sid, "bins")
        assert r.status_code == 200 and r.headers["content-type"].startswith("application/json")
        b = r.json()
        sf = scalefactors(sample)
        mpp = sf["microns_per_pixel"]
        assert b["sample"] == sid and b["n"] == a.n_obs == len(b["x"]) == len(b["y"]) == len(b["barcodes"])
        assert b["barcodes"] == list(a.obs_names)
        # x_um, y_um exactly as bone_analysis.per_bin_table writes them
        xy = np.asarray(a.obsm["spatial"], dtype=float) * mpp
        assert b["x"] == xy[:, 0].round(2).tolist() and b["y"] == xy[:, 1].round(2).tolist()
        f = b["fields"]
        assert set(f) == {"total_umi", "n_genes", "pct_mt", "passes_qc"}
        assert f["total_umi"] == {"type": "number", "values": a.obs["total_counts"].astype(int).tolist()}
        assert f["n_genes"] == {"type": "number", "values": a.obs["n_genes_by_counts"].astype(int).tolist()}
        assert f["pct_mt"]["type"] == "number"
        assert np.allclose(f["pct_mt"]["values"], a.obs["pct_counts_mt"].to_numpy(float), atol=0.006)
        kept = ba.qc_filter(a, 10, 15, 20.0).obs_names
        assert f["passes_qc"] == {"type": "bool", "values": a.obs_names.isin(kept).astype(int).tolist()}
        n_failing += a.n_obs - len(kept)
        s = b["summary"]
        totals = a.obs["total_counts"].to_numpy(float)
        assert s["n_bins"] == s["n_barcodes"] == a.n_obs and s["n_without_position"] == 0
        assert s["n_features"] == a.n_vars and s["total_umi"] == int(totals.sum())
        assert s["median_umi"] == float(np.median(totals)) and s["mean_umi"] == round(float(totals.mean()), 2)
        assert s["median_genes"] == float(np.median(a.obs["n_genes_by_counts"]))
        assert s["median_pct_mt"] == pytest.approx(float(np.median(a.obs["pct_counts_mt"])), abs=1e-3)
        assert s["n_pass_qc"] == len(kept) and s["pct_pass_qc"] == round(100 * len(kept) / a.n_obs, 2)
        assert s["bin_size"] == "016um" and s["microns_per_pixel"] == mpp
        assert s["qc"] == {"min_genes": 10, "min_counts": 15, "max_pct_mt": 20.0}
        assert Path(s["matrix"]) == bin_dir(sample) / "filtered_feature_bc_matrix.h5"
        assert Path(s["positions"]) == bin_dir(sample) / "spatial" / "tissue_positions.parquet"
        assert b["colors"] == {} and b["notes"] == []
        # the H&E image under the bins: hires pixels per micrometre of x, y
        image = b["image"]
        png = bin_dir(sample) / "spatial" / "tissue_hires_image.png"
        with Image.open(png) as im:
            assert (image["width"], image["height"]) == im.size
        assert image["px_per_um"] == pytest.approx(sf["tissue_hires_scalef"] / mpp, rel=1e-12)
        assert image["url"].startswith("/api/explore-images/")
        got = client.get(image["url"])
        assert got.status_code == 200 and got.headers["content-type"] == "image/png"
        assert got.content == png.read_bytes()
        assert explore(client, p["id"], sid, "bins").content == r.content  # cached
    assert n_failing > 0  # both outcomes occur
    assert explore(client, p["id"], "NOPE", "bins").status_code == 404
    assert client.post("/api/projects/p-00000000/explore/CTRL-1/bins").status_code == 404


def test_explore_gene_list(client, tiny_demo, raw):
    p = new_project(client, **demo_config(tiny_demo))
    a = raw["CTRL-1"]
    totals = np.rint(np.asarray(a.X.sum(axis=0)).ravel()).astype(np.int64)
    detected = np.asarray((a.X > 0).sum(axis=0)).ravel()
    order = np.lexsort((np.arange(a.n_vars), -totals))
    r = explore(client, p["id"], "CTRL-1", "genes").json()
    assert r["n_bins"] == a.n_obs and r["n_features"] == r["n_matches"] == a.n_vars and len(r["genes"]) == 50
    assert r["genes"] == [{"gene": a.var_names[i], "total_umi": int(totals[i]), "n_bins_detected": int(detected[i]),
                           "pct_bins_detected": round(100 * int(detected[i]) / a.n_obs, 2)} for i in order[:50]]
    # case-insensitive substring
    mt = explore(client, p["id"], "CTRL-1", "genes", query="MT-", limit=500).json()
    expected = [a.var_names[i] for i in order if "mt-" in a.var_names[i].lower()]
    assert [g["gene"] for g in mt["genes"]] == expected and mt["n_matches"] == len(expected) == 13
    assert len(explore(client, p["id"], "CTRL-1", "genes", limit=500).json()["genes"]) == 500
    assert explore(client, p["id"], "CTRL-1", "genes", query="no such gene").json()["genes"] == []
    for bad in ({"limit": 0}, {"limit": 501}, {"limit": "5"}, {"limit": True}, {"limit": 2.5}, {"query": 3},
                {"query": "x" * 201}):
        assert explore(client, p["id"], "CTRL-1", "genes", **bad).status_code == 400, bad
    assert explore(client, p["id"], "NOPE", "genes").status_code == 404


def test_explore_gene_map(client, tiny_demo, raw):
    p = new_project(client, **demo_config(tiny_demo))
    a = raw["TUM-1"]
    one = explore(client, p["id"], "TUM-1", "gene", genes=["Ctsk"]).json()
    expected = gene_counts(a, ["Ctsk"])
    assert one["genes_used"] == ["Ctsk"] and one["genes_missing"] == [] and one["n"] == a.n_obs
    assert one["indices"] == np.flatnonzero(expected).tolist() and np.array_equal(dense(one), expected)
    assert one["total"] == int(expected.sum()) and one["n_detected"] == int((expected > 0).sum()) > 0
    assert one["max"] == int(expected.max())
    # several genes are summed; any case; a gene named twice counts once; unknown names are reported
    two = explore(client, p["id"], "TUM-1", "gene", genes=["ctsk", "Acp5", "Nope", "Ctsk"]).json()
    assert two["genes_used"] == ["Ctsk", "Acp5"] and two["genes_missing"] == ["Nope"]
    assert np.array_equal(dense(two), gene_counts(a, ["Ctsk", "Acp5"]))
    none = explore(client, p["id"], "TUM-1", "gene", genes=["Nope"]).json()
    assert none == {"genes_used": [], "genes_missing": ["Nope"], "n": a.n_obs, "indices": [], "values": [],
                    "total": 0, "n_detected": 0, "max": 0}
    url = f"/api/projects/{p['id']}/explore/TUM-1/gene"
    assert client.post(url).status_code == 400
    for bad in ({}, {"genes": []}, {"genes": "Ctsk"}, {"genes": [1]}, {"genes": [" "]}, {"genes": ["x" * 201]},
                {"genes": [f"G{i}" for i in range(51)]}):
        assert client.post(url, json=bad).status_code == 400, bad
    assert explore(client, p["id"], "NOPE", "gene", genes=["Ctsk"]).status_code == 404


def test_explore_uses_the_unsaved_settings(client, tiny_demo, raw):
    import bone_analysis as ba
    p = new_project(client, **demo_config(tiny_demo))
    a = raw["CTRL-1"]
    cfg = copy.deepcopy(p["config"])
    cfg["qc"] = {**cfg["qc"], "min_genes": 30, "min_counts": 60, "max_pct_mt": 5}
    b = explore(client, p["id"], "CTRL-1", "bins", config=cfg).json()
    kept = ba.qc_filter(a, 30, 60, 5).obs_names
    assert b["fields"]["passes_qc"]["values"] == a.obs_names.isin(kept).astype(int).tolist()
    assert b["summary"]["qc"] == {"min_genes": 30, "min_counts": 60, "max_pct_mt": 5}
    saved = explore(client, p["id"], "CTRL-1", "bins").json()
    assert saved["summary"]["qc"] == {"min_genes": 10, "min_counts": 15, "max_pct_mt": 20.0}
    # a threshold that is not a number falls back to its default, with a note
    cfg["qc"]["min_genes"] = "ten"
    b = explore(client, p["id"], "CTRL-1", "bins", config=cfg).json()
    assert b["summary"]["qc"]["min_genes"] == 10 and b["notes"] == ["qc.min_genes is not a number; the default 10 "
                                                                    "is used"]
    cfg["qc"]["min_genes"] = 10 ** 30  # beyond int64
    assert explore(client, p["id"], "CTRL-1", "bins", config=cfg).json()["summary"]["n_pass_qc"] == 0
    # another bin size, another microns per pixel, a sample that is not saved yet
    cfg = copy.deepcopy(p["config"])
    cfg["bin_size"] = "008um"
    cfg["samples"].append({**cfg["samples"][0], "id": "EXTRA", "microns_per_pixel": 0.25})
    b = explore(client, p["id"], "EXTRA", "bins", config=cfg).json()
    fine = ba.read_visium_hd(bin_dir(tiny_demo["samples"][0], "008um"), "EXTRA")
    assert b["n"] == fine.n_obs and b["summary"]["bin_size"] == "008um" and b["summary"]["microns_per_pixel"] == 0.25
    assert b["x"] == (np.asarray(fine.obsm["spatial"], dtype=float)[:, 0] * 0.25).round(2).tolist()
    sf = scalefactors(tiny_demo["samples"][0], "008um")
    assert b["image"]["px_per_um"] == pytest.approx(sf["tissue_hires_scalef"] / 0.25)
    genes = explore(client, p["id"], "EXTRA", "genes", config=cfg).json()
    assert genes["n_bins"] == fine.n_obs
    m = explore(client, p["id"], "EXTRA", "gene", config=cfg, genes=["Ctsk"]).json()
    assert np.array_equal(dense(m), gene_counts(fine, ["Ctsk"]))
    assert explore(client, p["id"], "EXTRA", "bins").status_code == 404  # not in the saved configuration


def test_explore_problems(client, tiny_demo, tmp_path, monkeypatch):
    binned = copy_sample(tiny_demo, tmp_path / "C")
    p = new_project(client, samples=[{"id": "C", "label": "copy", "path": str(binned)}], reference_sample="C")
    pid = p["id"]
    first = explore(client, pid, "C", "bins").json()
    n_genes = explore(client, pid, "C", "genes").json()["n_bins"]
    assert first["n"] == n_genes

    # a changed file is read again (the caches are keyed by path, modification time and size)
    matrix = binned / "square_016um" / "filtered_feature_bc_matrix.h5"
    shutil.copyfile(bin_dir(tiny_demo["samples"][1]) / "filtered_feature_bc_matrix.h5", matrix)
    second = explore(client, pid, "C", "bins").json()
    assert second["n"] != first["n"] and explore(client, pid, "C", "genes").json()["n_bins"] == second["n"]

    # no H&E image at any bin size: the maps still work, without an image
    for png in binned.glob("square_*/spatial/tissue_hires_image.png"):
        png.unlink()
    b = explore(client, pid, "C", "bins").json()
    assert b["image"] is None and b["notes"] == ["no H&E image (spatial/tissue_hires_image.png)"]
    assert client.post(f"/api/projects/{pid}/explore/samples").json()["samples"][0]["has_image"] is False

    # microns per pixel neither configured nor recorded
    sf_path = binned / "square_016um" / "spatial" / "scalefactors_json.json"
    sf = json.loads(sf_path.read_text(encoding="utf-8"))
    sf_path.write_text(json.dumps({k: v for k, v in sf.items() if k != "microns_per_pixel"}), encoding="utf-8")
    r = explore(client, pid, "C", "bins")
    assert r.status_code == 422 and "microns_per_pixel" in r.json()["detail"]
    sf_path.write_text(json.dumps(sf), encoding="utf-8")

    # positions of another bin size: almost no barcode has a position
    positions = binned / "square_016um" / "spatial" / "tissue_positions.parquet"
    saved_positions = positions.read_bytes()
    shutil.copyfile(binned / "square_008um" / "spatial" / "tissue_positions.parquet", positions)
    r = explore(client, pid, "C", "bins")
    assert r.status_code == 422 and "have a position" in r.json()["detail"]
    positions.write_bytes(saved_positions)

    # an unreadable matrix
    matrix.write_bytes(b"not an HDF5 file")
    for what, body in (("bins", {}), ("genes", {}), ("gene", {"genes": ["Ctsk"]})):
        r = explore(client, pid, "C", what, **body)
        assert r.status_code == 422 and "not a readable feature-barcode matrix" in r.json()["detail"], what

    # a section larger than the explorer reads (002um sections have millions of bins): refused before its counts are
    # read, and reported by explore/samples
    matrix.write_bytes((bin_dir(tiny_demo["samples"][0]) / "filtered_feature_bc_matrix.h5").read_bytes())
    n = explore(client, pid, "C", "bins").json()["n"]
    monkeypatch.setattr(server, "RAW_SECTION_MAX_BINS", n - 1)
    st = matrix.stat()  # a new modification time: not the cached section
    os.utime(matrix, ns=(st.st_atime_ns, st.st_mtime_ns + 10 ** 9))
    r = explore(client, pid, "C", "genes")
    assert r.status_code == 422 and f"at most {n - 1:,} bins" in r.json()["detail"]
    entry = client.post(f"/api/projects/{pid}/explore/samples").json()["samples"][0]
    assert entry["ok"] is False and "choose a larger bin size" in entry["message"]
    monkeypatch.undo()

    # no positions, then no matrix
    positions.unlink()
    r = explore(client, pid, "C", "bins")
    assert r.status_code == 404 and "no tissue positions file" in r.json()["detail"]
    entry = client.post(f"/api/projects/{pid}/explore/samples").json()["samples"][0]
    assert entry["n_bins"] is None and entry["ok"] is False and "no tissue positions file" in entry["message"]
    matrix.unlink()
    r = explore(client, pid, "C", "genes")
    assert r.status_code == 404 and "square_016um/filtered_feature_bc_matrix.h5" in r.json()["detail"]


def test_explore_keeps_the_sections_in_use(client, tiny_demo, monkeypatch):
    """Gene maps that cycle over four sections (two samples at two bin sizes: the 4-panel viewer, a changed draft bin
    size) read each section once."""
    built, section = [], server.RawSection

    def counting(matrix, positions):
        built.append(matrix)
        return section(matrix, positions)

    monkeypatch.setattr(server, "RawSection", counting)
    p = new_project(client, **demo_config(tiny_demo))
    fine = copy.deepcopy(p["config"])
    fine["bin_size"] = "008um"
    for _ in range(3):
        for cfg in (p["config"], fine):
            for sid in ("CTRL-1", "TUM-1"):
                assert explore(client, p["id"], sid, "gene", config=cfg, genes=["Ctsk"]).status_code == 200
    assert len(built) == len(set(built)) == 4


def test_explore_requests_share_a_failed_read(client, tiny_demo, tmp_path, monkeypatch):
    """The requests a map sends together (bins, genes, gene) wait for one read of an unreadable matrix and all get its
    error; the error is not kept, so a later request reads the file again."""
    import threading
    import time
    binned = copy_sample(tiny_demo, tmp_path / "B")
    p = new_project(client, samples=[{"id": "B", "label": "broken", "path": str(binned)}], reference_sample="B")
    (binned / "square_016um" / "filtered_feature_bc_matrix.h5").write_bytes(b"not an HDF5 file")
    reads, release, read_counts = [], threading.Event(), server.read_counts
    caches, cache_get = set(), server.LRUCache.get

    def slow(path, max_bins=None):
        reads.append(path)
        release.wait(10)
        return read_counts(path, max_bins)

    def get(self, key, build):
        caches.add(self)
        return cache_get(self, key, build)

    def joined() -> int:
        """Requests waiting on (or doing) a build: each notes the failed builds so far as it joins, under the lock."""
        n = 0
        for cache in list(caches):
            with cache._lock:
                n += sum(gate[1] for gate in cache._gates.values())
        return n

    monkeypatch.setattr(server, "read_counts", slow)
    monkeypatch.setattr(server.LRUCache, "get", get)
    asks = [("bins", {}), ("genes", {}), ("gene", {"genes": ["Ctsk"]})] * 2
    answers = [None] * len(asks)

    def ask(i: int, what: str, body: dict):
        answers[i] = explore(client, p["id"], "B", what, **body)

    threads = [threading.Thread(target=ask, args=(i, *a)) for i, a in enumerate(asks)]
    for t in threads:
        t.start()
    # all of them joined before the read fails: the two bins requests on the payload, the payload's build and the four
    # others on the section
    try:
        deadline = time.monotonic() + 10
        while joined() < len(asks) + 1 and time.monotonic() < deadline:
            time.sleep(0.01)
        assert joined() == len(asks) + 1 and len(reads) == 1
    finally:
        release.set()
        for t in threads:
            t.join()
    assert len(reads) == 1
    for (what, _), r in zip(asks, answers):
        assert r.status_code == 422 and "not a readable feature-barcode matrix" in r.json()["detail"], what
    assert explore(client, p["id"], "B", "genes").status_code == 422 and len(reads) == 2


def test_explore_images_are_served_by_key_only(client, tiny_demo):
    p = new_project(client, **demo_config(tiny_demo))
    url = explore(client, p["id"], "CTRL-1", "bins").json()["image"]["url"]
    key = url.rsplit("/", 1)[1]
    assert len(key) == 32 and client.get(url).status_code == 200
    for other in ("0" * 32, key.upper(), key[:-1], "..%2F..%2Fpyproject.toml", "tissue_hires_image.png"):
        assert client.get(f"/api/explore-images/{other}").status_code == 404, other
    # a key stays the same while the file is unchanged
    assert explore(client, p["id"], "CTRL-1", "bins").json()["image"]["url"] == url


def test_explore_keeps_the_analysis_out_of_the_server(tmp_path, tiny_demo):
    code = textwrap.dedent("""
        import json, sys, warnings
        warnings.simplefilter("ignore")
        from fastapi.testclient import TestClient
        from webapp.server import create_app
        config = json.loads(sys.argv[2])
        with TestClient(create_app(sys.argv[1]), base_url="http://127.0.0.1") as c:
            pid = c.post("/api/projects", json={"name": "x", "config": config}).json()["id"]
            for what, body in (("bins", {}), ("genes", {}), ("gene", {"genes": ["Ctsk"]})):
                assert c.post(f"/api/projects/{pid}/explore/CTRL-1/{what}", json=body).status_code == 200
        print(sorted(m for m in ("scanpy", "anndata", "squidpy", "bone_analysis", "bone_figures") if m in sys.modules))
    """)
    out = subprocess.run([sys.executable, "-c", code, str(tmp_path / "ws"), json.dumps(demo_config(tiny_demo))],
                         cwd=ROOT, capture_output=True, text=True, timeout=180)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip().splitlines()[-1] == "[]"


# ============================================================================
# GENE MAPS OF A RUN
# ============================================================================

def joint_per_bin(cfg: dict):
    """(joint AnnData of raw counts, per-bin table) built as bone_analysis.run_analysis_from_config builds them:
    load_sections, anndata.concat(index_unique="-"), per_bin_table."""
    import anndata as ad
    import bone_analysis as ba
    adatas, _, _ = ba.load_sections(cfg, cfg["bin_size"])
    mpp = {sid: float(x.uns["spatial"][sid]["scalefactors"]["microns_per_pixel"]) for sid, x in adatas.items()}
    a = ad.concat(adatas, label="sample_key", index_unique="-", uns_merge="unique", join="outer")
    a.obs["sample"] = pd.Categorical(a.obs["sample"].astype(str), categories=list(adatas))
    a.uns["microns_per_pixel"] = mpp
    a.obs["in_marrow_roi"] = True
    return a, ba.per_bin_table(a, cfg)


def write_per_bin(results: Path, table: pd.DataFrame) -> None:
    with gzip.open(results / "per_bin_signature_calls_16um.csv.gz", "wt", newline="", encoding="utf-8") as fh:
        table.to_csv(fh, index=False)


@pytest.fixture(scope="module")
def gene_run(tmp_path_factory, tiny_demo):
    """A succeeded (fake) run of the tiny demo whose per-bin table is the analysis' own for the raw data, rows
    shuffled, a few left out, one barcode altered: (client, run id, joint AnnData, per-bin table written)."""
    app = create_app(tmp_path_factory.mktemp("ws"), worker_module=FAKE)
    with TestClient(app, base_url=LOCAL) as c:
        p = new_project(c, "Genes", **demo_config(tiny_demo))
        rid = start_run(c, p["id"])["id"]
        assert wait_for(c, rid, {"succeeded", "failed"})["status"] == "succeeded"
        a, table = joint_per_bin(pc.merge_with_defaults(demo_config(tiny_demo)))
        table = table.sample(frac=1.0, random_state=0).iloc[5:].reset_index(drop=True)
        table.loc[3, "barcode"] = "s_016um_99999_99999-" + table.loc[3, "sample_id"]
        write_per_bin(app.state.storage.results_dir(rid), table)
        yield c, rid, a, table


def test_run_gene_map_holds_the_raw_counts_of_the_run_bins(gene_run):
    c, rid, a, table = gene_run
    for sid in ("CTRL-1", "TUM-1"):
        rows = table[table["sample_id"] == sid]
        assert c.get(f"/api/runs/{rid}/bins", params={"sample": sid}).json()["n"] == len(rows)
        r = c.post(f"/api/runs/{rid}/gene", json={"sample": sid, "genes": ["Ctsk", "acp5", "Nope"]})
        assert r.status_code == 200, r.text
        m = r.json()
        unmatched = ~rows["barcode"].isin(a.obs_names).to_numpy()
        expected = np.zeros(len(rows), np.int64)
        expected[~unmatched] = gene_counts(a[rows["barcode"][~unmatched].tolist()], ["Ctsk", "Acp5"])
        assert m["genes_used"] == ["Ctsk", "Acp5"] and m["genes_missing"] == ["Nope"] and m["n"] == len(rows)
        assert np.array_equal(dense(m), expected) and m["total"] == int(expected.sum()) > 0
        assert m["n_unmatched"] == int(unmatched.sum())
    assert sum(c.post(f"/api/runs/{rid}/gene", json={"sample": s, "genes": ["Ctsk"]}).json()["n_unmatched"]
               for s in ("CTRL-1", "TUM-1")) == 1


def test_run_gene_list(gene_run):
    c, rid, a, table = gene_run
    rows = table[(table["sample_id"] == "TUM-1") & table["barcode"].isin(a.obs_names)]
    sub = a[rows["barcode"].tolist()]
    totals = np.rint(np.asarray(sub.X.sum(axis=0)).ravel()).astype(np.int64)
    detected = np.asarray((sub.X > 0).sum(axis=0)).ravel()
    n = int((table["sample_id"] == "TUM-1").sum())
    r = c.get(f"/api/runs/{rid}/genes", params={"sample": "TUM-1", "query": "ct", "limit": 5}).json()
    order = [i for i in np.lexsort((np.arange(a.n_vars), -totals)) if "ct" in a.var_names[i].lower()]
    assert r["n_bins"] == n and r["n_features"] == a.n_vars and r["n_matches"] == len(order)
    assert r["n_unmatched"] == int((~table.loc[table["sample_id"] == "TUM-1", "barcode"].isin(a.obs_names)).sum())
    assert r["genes"] == [{"gene": a.var_names[i], "total_umi": int(totals[i]), "n_bins_detected": int(detected[i]),
                           "pct_bins_detected": round(100 * int(detected[i]) / n, 2)} for i in order[:5]]
    assert len(c.get(f"/api/runs/{rid}/genes", params={"sample": "TUM-1"}).json()["genes"]) == 50


def test_run_bins_image(gene_run, tiny_demo):
    c, rid, _, _ = gene_run
    sample = tiny_demo["samples"][1]
    info = c.get(f"/api/runs/{rid}/bins/image", params={"sample": sample["id"]}).json()
    sf = scalefactors(sample)
    png = bin_dir(sample) / "spatial" / "tissue_hires_image.png"
    with Image.open(png) as im:
        assert (info["width"], info["height"]) == im.size
    assert info["px_per_um"] == pytest.approx(sf["tissue_hires_scalef"] / sf["microns_per_pixel"], rel=1e-12)
    assert c.get(info["url"]).content == png.read_bytes()
    assert c.get(f"/api/runs/{rid}/bins/image", params={"sample": "NOPE"}).status_code == 404
    assert c.get(f"/api/runs/{rid}/bins/image").status_code == 400


def test_run_gene_requests_are_checked(gene_run):
    c, rid, _, _ = gene_run
    url = f"/api/runs/{rid}/gene"
    for bad in ({"genes": ["Ctsk"]}, {"sample": 1, "genes": ["Ctsk"]}, {"sample": "TUM-1"},
                {"sample": "TUM-1", "genes": []}, {"sample": "TUM-1", "genes": [f"G{i}" for i in range(51)]}):
        assert c.post(url, json=bad).status_code == 400, bad
    assert c.post(url, json={"sample": "NOPE", "genes": ["Ctsk"]}).status_code == 404
    assert c.post("/api/runs/r-20000101-000000-0000/gene", json={"sample": "TUM-1", "genes": ["Ctsk"]}
                  ).status_code == 404
    genes = f"/api/runs/{rid}/genes"
    for params in ({"sample": "TUM-1", "limit": 0}, {"sample": "TUM-1", "limit": 501},
                   {"sample": "TUM-1", "limit": "x"}, {}):
        assert c.get(genes, params=params).status_code == 400, params
    assert c.get(genes, params={"sample": "NOPE"}).status_code == 404


def test_run_gene_maps_when_the_raw_data_moved(client, tiny_demo, tmp_path):
    binned = copy_sample(tiny_demo, tmp_path / "moving")
    p = new_project(client, samples=[{"id": "M", "label": "moving", "path": str(binned)}], reference_sample="M")
    rid = start_run(client, p["id"])["id"]
    assert wait_for(client, rid, {"succeeded", "failed"})["status"] == "succeeded"
    shutil.rmtree(tmp_path / "moving")
    r = client.post(f"/api/runs/{rid}/gene", json={"sample": "M", "genes": ["Ctsk"]})
    assert r.status_code == 422 and "can no longer be found" in r.json()["detail"]
    assert client.get(f"/api/runs/{rid}/genes", params={"sample": "M"}).status_code == 422
    assert client.get(f"/api/runs/{rid}/bins/image", params={"sample": "M"}).status_code == 404
    assert client.get(f"/api/runs/{rid}/bins", params={"sample": "M"}).status_code == 200  # the table still loads


@pytest.mark.slow
def test_run_gene_map_of_a_real_analysis(tiny_joint, tmp_path, tiny_demo, raw):
    """The per-bin table of a whole analysis run (tiny_joint): every barcode matches, and the gene map and the totals
    are the raw counts of the bins."""
    app = create_app(tmp_path / "ws", worker_module=FAKE)
    with TestClient(app, base_url=LOCAL) as c:
        p = new_project(c, **demo_config(tiny_demo))
        rid = start_run(c, p["id"])["id"]
        assert wait_for(c, rid, {"succeeded", "failed"})["status"] == "succeeded"
        shutil.copyfile(tiny_joint["out"] / "per_bin_signature_calls_16um.csv.gz",
                        app.state.storage.results_dir(rid) / "per_bin_signature_calls_16um.csv.gz")
        table = pd.read_csv(tiny_joint["out"] / "per_bin_signature_calls_16um.csv.gz")
        for sid, a in raw.items():
            rows = table[table["sample_id"] == sid]
            barcodes = [b[:-len(sid) - 1] for b in rows["barcode"]]  # anndata.concat appended -<sample id>
            sub = a[barcodes]
            m = c.post(f"/api/runs/{rid}/gene", json={"sample": sid, "genes": ["Ctsk"]}).json()
            assert m["n_unmatched"] == 0 and m["n"] == len(rows)
            assert np.array_equal(dense(m), gene_counts(sub, ["Ctsk"]))
            assert np.array_equal(rows["total_umi"].to_numpy(float), sub.obs["total_counts"].to_numpy(float))
