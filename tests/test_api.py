"""
tests/test_api.py
=================
HTTP API of Bone Spatial Studio (docs/ARCHITECTURE.md, sections 4-6) against temporary workspaces, with
tests/fake_worker.py standing in for the analysis worker (no scanpy needed). The data endpoints use the tiny synthetic
dataset of the tiny_demo fixture (tests/conftest.py).
"""
from __future__ import annotations

import gzip
import io
import json
import platform
import shutil
import subprocess
import sys
import textwrap
import time
import warnings
from pathlib import Path

import h5py
import pandas as pd
import pytest
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]  # on sys.path through tests/conftest.py

with warnings.catch_warnings():  # Starlette 1.x prefers httpx2 for its TestClient; httpx still works
    warnings.filterwarnings("ignore", message=".*httpx.*")
    from fastapi.testclient import TestClient

import pipeline_config as pc  # noqa: E402
from tests import fake_worker  # noqa: E402
from webapp import __version__, server  # noqa: E402
from webapp.server import DEMO_NAME, create_app  # noqa: E402
from webapp.storage import Storage  # noqa: E402

FAKE = "tests.fake_worker"
LOCAL = "http://127.0.0.1:8765"
TIMEOUT = 60.0
TOKEN = "s3cret-Token_1"
SQUARE = {"type": "FeatureCollection", "features": [
    {"type": "Feature", "properties": {"frame": "hires"},
     "geometry": {"type": "Polygon", "coordinates": [[[10, 10], [300, 10], [300, 300], [10, 300]]]}}]}


# ============================================================================
# HELPERS
# ============================================================================

def demo_config(tiny_demo: dict, **extra) -> dict:
    """Project configuration (before merging with the defaults) of the tiny demo sections."""
    samples = [dict(s) for s in tiny_demo["samples"]]
    return {"data_root": str(tiny_demo["root"]), "samples": samples, "reference_sample": samples[0]["id"], **extra}


@pytest.fixture
def client(tmp_path):
    with TestClient(create_app(tmp_path / "ws", worker_module=FAKE), base_url=LOCAL) as c:
        yield c


def new_project(client, name: str = "P", **config) -> dict:
    r = client.post("/api/projects", json={"name": name, "config": config})
    assert r.status_code == 200, r.text
    return r.json()


def start_run(client, pid: str, **body) -> dict:
    r = client.post(f"/api/projects/{pid}/runs", json=body)
    assert r.status_code == 200, r.text
    return r.json()


def wait_for(client, rid: str, statuses: set[str], label: str | None = None) -> dict:
    """Poll a run until its status is one of `statuses` (and, with `label`, its current step is `label`)."""
    deadline = time.monotonic() + TIMEOUT
    while True:
        rec = client.get(f"/api/runs/{rid}").json()
        if rec["status"] in statuses and (label is None or rec["progress"]["label"] == label):
            return rec
        assert time.monotonic() < deadline, f"run {rid} stuck: {rec['status']} {rec['progress']} {rec['error']}"
        time.sleep(0.05)


def full_log(client, rid: str) -> tuple[str, list[dict]]:
    """The whole log of a finished run read through GET /log, and the responses."""
    offset, parts, responses = 0, [], []
    while True:
        r = client.get(f"/api/runs/{rid}/log", params={"offset": offset}).json()
        responses.append(r)
        parts.append(r["text"])
        if r["complete"]:
            return "".join(parts), responses
        assert r["offset"] > offset, "the log did not advance although the run has finished"
        offset = r["offset"]


# ============================================================================
# ENVIRONMENT, FRONT END, DOCUMENTS
# ============================================================================

def test_health(client):
    h = client.get("/api/health").json()
    assert h["ok"] is True and h["app_version"] == __version__ == "1.0.0"
    assert h["python"] == platform.python_version() and h["platform"]
    assert Path(h["workspace"]) == client.app.state.storage.root
    assert list(h["packages"]) == list(server.PACKAGES)
    assert h["packages"]["pandas"] == pd.__version__ and h["packages"]["h5py"] == h5py.__version__
    assert h["missing_required"] == [n for n in server.REQUIRED_PACKAGES if h["packages"][n] is None]
    assert h["missing_optional"] == [n for n in server.OPTIONAL_PACKAGES if h["packages"][n] is None]
    assert h["active_run"] is None and h["queued"] == 0


def test_server_never_imports_the_analysis(tmp_path):
    code = textwrap.dedent("""
        import sys, warnings
        warnings.simplefilter("ignore")
        from fastapi.testclient import TestClient
        from webapp.server import create_app
        with TestClient(create_app(sys.argv[1]), base_url="http://127.0.0.1") as c:
            assert c.get("/api/health").status_code == 200
            pid = c.post("/api/projects", json={"name": "x"}).json()["id"]
            assert c.post(f"/api/projects/{pid}/validate").status_code == 200
        print(sorted(m for m in ("scanpy", "anndata", "squidpy", "bone_analysis", "bone_figures") if m in sys.modules))
    """)
    out = subprocess.run([sys.executable, "-c", code, str(tmp_path / "ws")], cwd=ROOT, capture_output=True,
                         text=True, timeout=180)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip().splitlines()[-1] == "[]"


def test_defaults_and_docs(client):
    d = client.get("/api/defaults").json()
    assert d["config"] == json.loads(json.dumps(pc.default_config()))
    assert d["param_schema"] == json.loads(json.dumps(pc.PARAM_SCHEMA))
    assert d["step_labels"] == pc.STEP_LABELS
    for name, rel in server.DOCS.items():
        r = client.get(f"/api/docs/{name}")
        if (ROOT / rel).is_file():
            assert r.status_code == 200 and r.json() == {"name": name,
                                                          "markdown": (ROOT / rel).read_text(encoding="utf-8")}
        else:
            assert r.status_code == 404
    assert client.get("/api/docs/architecture").json()["markdown"].startswith("# Bone Spatial Studio")
    r = client.get("/api/docs/secrets")
    assert r.status_code == 404 and "readme" in r.json()["detail"]


def test_front_end(client):
    r = client.get("/")
    assert r.status_code == 200 and r.headers["content-type"].startswith("text/html")
    assert r.headers["cache-control"] == "no-cache" and '<div id="app">' in r.text
    for path, media_type in (("/static/js/core.js", "text/javascript"), ("/static/css/app.css", "text/css")):
        r = client.get(path)
        assert r.status_code == 200 and r.headers["content-type"].startswith(media_type)
        assert r.headers["cache-control"] == "no-cache"
    assert client.get("/static/nothing.js").status_code == 404
    r = client.get("/api/nothing")
    assert r.status_code == 404 and "detail" in r.json()


# ============================================================================
# PROJECTS AND REGIONS OF INTEREST
# ============================================================================

def test_projects_crud(client):
    assert client.get("/api/projects").json() == []
    a = new_project(client, "Alpha")
    assert a["id"].startswith("p-") and len(a["id"]) == 10
    assert a["config"] == json.loads(json.dumps({**pc.default_config(), "name": "Alpha"}))
    assert a["rois"] == [] and a["validation"] == {"errors": [], "warnings": []}
    one = {"id": "S1", "label": "only", "condition": "control", "color": "#123456", "path": "S1",
           "microns_per_pixel": 0.5}
    b = new_project(client, "Beta", qc={"min_genes": 5}, samples=[one], reference_sample="S1")
    assert b["config"]["qc"] == {**pc.default_config()["qc"], "min_genes": 5}
    assert b["config"]["samples"] == [one]

    listing = client.get("/api/projects").json()
    assert [p["id"] for p in listing] == [b["id"], a["id"]]  # newest first
    assert listing[0] == {"id": b["id"], "name": "Beta", "created": b["created"], "updated": b["updated"],
                          "n_samples": 1, "data_root": ".", "last_run": None}
    assert listing[1]["n_samples"] == 2

    # drafts with errors are saved; the response lists the errors
    cfg = a["config"]
    cfg["qc"]["min_genes"] = -1
    r = client.put(f"/api/projects/{a['id']}", json={"name": "Alpha 2", "config": cfg})
    assert r.status_code == 200
    put = r.json()
    assert put["name"] == "Alpha 2" and put["updated"] >= a["updated"]
    assert any("qc.min_genes" in e for e in put["validation"]["errors"])
    got = client.get(f"/api/projects/{a['id']}").json()
    assert got["config"]["qc"]["min_genes"] == -1 and got["validation"] == put["validation"]
    # a partial config is merged with the defaults
    put = client.put(f"/api/projects/{a['id']}", json={"config": {"moran": {"n_genes": 300}}}).json()
    assert put["name"] == "Alpha 2" and put["config"]["moran"] == {"n_genes": 300, "n_consensus": 200}
    assert put["config"]["qc"]["min_genes"] == 10 and put["validation"]["errors"] == []

    client.put(f"/api/projects/{b['id']}/roi/S1", json=SQUARE)
    dup = client.post(f"/api/projects/{b['id']}/duplicate").json()
    assert dup["id"] not in (a["id"], b["id"]) and dup["name"] == "Beta (copy)"
    assert dup["config"] == b["config"] and dup["rois"] == ["S1"]
    roi = client.get(f"/api/projects/{b['id']}/roi/S1").json()
    assert client.get(f"/api/projects/{dup['id']}/roi/S1").json() == roi

    assert client.delete(f"/api/projects/{a['id']}").json() == {"ok": True}
    assert client.get(f"/api/projects/{a['id']}").status_code == 404
    assert client.delete(f"/api/projects/{a['id']}").status_code == 404
    assert [p["id"] for p in client.get("/api/projects").json()] == [dup["id"], b["id"]]


def test_project_request_errors(client):
    pid = new_project(client)["id"]
    url = f"/api/projects/{pid}"
    for body in ([1, 2], "text", None):
        assert client.put(url, json=body).status_code == 400
    r = client.put(url, content=b"{not json", headers={"content-type": "application/json"})
    assert r.status_code == 400 and "not valid JSON" in r.json()["detail"]
    assert client.put(url, json={"config": [1]}).status_code == 400
    assert client.put(url, json={"name": " "}).status_code == 400
    assert client.post("/api/projects", json={}).status_code == 400
    assert client.post("/api/projects", json={"name": "x" * 201}).status_code == 400
    assert client.post("/api/projects", json={"config": {"name": "From config"}}).json()["name"] == "From config"
    for missing in ("p-00000000", "nonsense", "..%2F..%2Fruns"):
        assert client.get(f"/api/projects/{missing}").status_code == 404
    assert client.put("/api/projects/p-00000000", json={"name": "y"}).status_code == 404
    # a malformed draft is saved and reported, never a server error
    r = client.put(url, json={"config": {"samples": "not a list", "roi": 5}})
    assert r.status_code == 200 and r.json()["validation"]["errors"]
    assert client.post(f"{url}/validate").json()["ok"] is False


def test_roi(client):
    p = new_project(client)
    sid = p["config"]["samples"][0]["id"]
    url = f"/api/projects/{p['id']}/roi/{sid}"
    assert client.get(url).status_code == 404
    fc = {"type": "FeatureCollection", "features": [
        SQUARE["features"][0],
        {"type": "Feature", "properties": {"name": "lesion"}, "geometry": {"type": "MultiPolygon", "coordinates": [
            [[[400, 400], [500, 400], [450, 480], [400, 400]]],
            [[[0, 0], [90, 0], [90, 90], [0, 90]], [[10, 10], [20, 10], [20, 20]]]]}}]}
    assert client.put(url, json=fc).json() == {"ok": True, "n_polygons": 3}
    got = client.get(url).json()
    assert got["type"] == "FeatureCollection" and len(got["features"]) == 2
    ring = got["features"][0]["geometry"]["coordinates"][0]
    assert ring[0] == ring[-1] == [10.0, 10.0] and len(ring) == 5  # closed
    assert got["features"][1]["properties"] == {"name": "lesion", "frame": "hires"}
    assert got["features"][1]["geometry"]["coordinates"][1][1][-1] == [10.0, 10.0]  # the hole, closed
    assert client.get(f"/api/projects/{p['id']}").json()["rois"] == [sid]

    bad = [
        {"type": "Feature", "geometry": SQUARE["features"][0]["geometry"]},
        {"type": "FeatureCollection", "features": [{"type": "Feature", "geometry": {"type": "Point",
                                                                                     "coordinates": [1, 2]}}]},
        {"type": "FeatureCollection", "features": [{"type": "Feature", "geometry": {
            "type": "Polygon", "coordinates": [[[0, 0], [1, 1], [0, 0]]]}}]},
        {"type": "FeatureCollection", "features": [{"type": "Feature", "geometry": {
            "type": "Polygon", "coordinates": [[[0, 0], [1, float("nan")], [2, 2]]]}}]},
        {"type": "FeatureCollection", "features": [{"type": "Feature", "geometry": {
            "type": "Polygon", "coordinates": [[[0, 0], ["1", 1], [2, 2]]]}}]},
        {"type": "FeatureCollection", "features": [{"type": "Feature", "properties": {"frame": "screen"},
                                                    "geometry": SQUARE["features"][0]["geometry"]}]},
        {"type": "FeatureCollection", "features": ["polygon"]},
    ]
    for body in bad:  # raw JSON: NaN is not valid JSON, but Python's json module writes and reads it
        r = client.put(url, content=json.dumps(body), headers={"content-type": "application/json"})
        assert r.status_code == 400 and r.json()["detail"], body
    assert client.put(url, json=[SQUARE]).status_code == 400
    assert client.get(url).json() == got  # unchanged by the rejected uploads
    assert client.put(f"/api/projects/{p['id']}/roi/bad%20id", json=SQUARE).status_code == 400
    assert client.put(f"/api/projects/p-00000000/roi/{sid}", json=SQUARE).status_code == 404

    # an empty collection removes the region
    assert client.put(url, json={"type": "FeatureCollection", "features": []}).json() == {"ok": True,
                                                                                         "n_polygons": 0}
    assert client.get(url).status_code == 404
    client.put(url, json=SQUARE)
    assert client.delete(url).json() == {"ok": True}
    assert client.get(url).status_code == 404 and client.delete(url).status_code == 200
    assert client.get(f"/api/projects/{p['id']}").json()["rois"] == []


# ============================================================================
# DATA CHECKS
# ============================================================================

def test_fs_list(client, tmp_path):
    base = tmp_path / "browse"
    for d in ("SR1_A/outs", "plain", "bins/square_016um", ".hidden", "Zeta"):
        (base / d).mkdir(parents=True)
    (base / "file.txt").write_text("x", encoding="utf-8")
    r = client.get("/api/fs/list", params={"path": str(base)}).json()
    assert r["path"] == str(base) and r["parent"] == str(base.parent) and r["roots"]
    assert [e["name"] for e in r["entries"]] == ["bins", "plain", "SR1_A", "Zeta"]
    assert {e["name"]: e["looks_like_spaceranger"] for e in r["entries"]} == {"bins": True, "plain": False,
                                                                              "SR1_A": True, "Zeta": False}
    assert all(e["is_dir"] and e["path"] == str(base / e["name"]) for e in r["entries"])
    assert client.get("/api/fs/list").json()["path"] == str(Path.home())
    root = client.get("/api/fs/list", params={"path": r["roots"][0]}).json()
    assert root["parent"] is None
    assert client.get("/api/fs/list", params={"path": str(base / "missing")}).status_code == 404
    assert client.get("/api/fs/list", params={"path": str(base / "file.txt")}).status_code == 400


def test_detect_samples(client, tiny_demo, tmp_path):
    p = new_project(client, **demo_config(tiny_demo))
    r = client.post(f"/api/projects/{p['id']}/detect-samples").json()
    assert Path(r["data_root"]) == Path(tiny_demo["root"])
    assert r["candidates"] == [
        {"folder": "DEMO_CTRL-1", "path": "DEMO_CTRL-1/outs/binned_outputs", "suggested_id": "DEMO_CTRL-1",
         "bin_sizes": ["008um", "016um"]},
        {"folder": "DEMO_TUM-1", "path": "DEMO_TUM-1/outs/binned_outputs", "suggested_id": "DEMO_TUM-1",
         "bin_sizes": ["008um", "016um"]}]

    runs = tmp_path / "sr"
    for d in ("SR007144_LL-2-1/outs/binned_outputs/square_002um", "SR007144_LL-2-1/outs/binned_outputs/square_016um",
              "flat/square_008um", "a/b/c/d/e/square_016um"):  # the last one is too deep
        (runs / d).mkdir(parents=True)
        (runs / d / "filtered_feature_bc_matrix.h5").write_bytes(b"")
    (runs / "empty/square_016um").mkdir(parents=True)  # no matrix
    r = client.post(f"/api/projects/{p['id']}/detect-samples", json={"data_root": str(runs)}).json()
    assert r["candidates"] == [
        {"folder": "flat", "path": "flat", "suggested_id": "flat", "bin_sizes": ["008um"]},
        {"folder": "SR007144_LL-2-1", "path": "SR007144_LL-2-1/outs/binned_outputs", "suggested_id": "LL-2-1",
         "bin_sizes": ["002um", "016um"]}]
    r = client.post(f"/api/projects/{p['id']}/detect-samples", json={"data_root": str(runs / "nowhere")})
    assert r.status_code == 400
    assert client.post("/api/projects/p-00000000/detect-samples").status_code == 404


def test_validate(client, tiny_demo):
    p = new_project(client, **demo_config(tiny_demo))
    r = client.post(f"/api/projects/{p['id']}/validate").json()
    assert r["ok"] is True and r["errors"] == []
    ctrl = r["samples"][0]
    bin_dir = Path(tiny_demo["samples"][0]["path"]) / "square_016um"
    with h5py.File(bin_dir / "filtered_feature_bc_matrix.h5", "r") as f:
        n_bins, n_genes = int(f["matrix/shape"][1]), int(f["matrix/shape"][0])
    scalefactors = json.loads((bin_dir / "spatial" / "scalefactors_json.json").read_text(encoding="utf-8"))
    assert ctrl["id"] == "CTRL-1" and Path(ctrl["bin_dir"]) == bin_dir
    assert ctrl["bin_sizes"] == ["008um", "016um"] and ctrl["n_bins"] == n_bins
    assert ctrl["has_positions"] and ctrl["has_scalefactors"] and ctrl["has_hires_image"]
    assert ctrl["microns_per_pixel"] == scalefactors["microns_per_pixel"]
    assert f"{n_bins:,} bins at 016um" in ctrl["message"]
    genes = r["genes"]
    assert genes["n_features"] == n_genes
    assert genes["signatures"]["osteocyte"]["missing"] == ["Fgf23"]
    assert genes["signatures"]["osteocyte"]["core_present"] == 3
    assert genes["signatures"]["osteoclast"]["missing"] == [] and genes["receptor_panel"]["missing"] == []
    assert any(w.startswith("signature osteocyte: 1 of 5 genes missing") for w in r["warnings"])

    # unsaved edits are checked instead of the saved configuration
    cfg = p["config"]
    cfg["samples"][1]["path"] = "nowhere"
    cfg["samples"][0]["microns_per_pixel"] = 0.25
    r = client.post(f"/api/projects/{p['id']}/validate", json={"config": cfg}).json()
    assert r["ok"] is False and any("TUM-1" in e for e in r["errors"])
    assert r["samples"][0]["microns_per_pixel"] == 0.25
    tum = r["samples"][1]
    assert tum["bin_dir"] is None and tum["n_bins"] is None and "square_016um" in tum["message"]
    assert client.get(f"/api/projects/{p['id']}").json()["config"]["samples"][1]["path"] != "nowhere"
    # no readable sample: no gene report
    cfg["samples"] = cfg["samples"][1:]
    cfg["reference_sample"] = "TUM-1"
    genes = client.post(f"/api/projects/{p['id']}/validate", json={"config": cfg}).json()["genes"]
    assert genes == {"n_features": None, "signatures": {}, "receptor_panel": {"missing": []},
                     "adipocyte_genes": {"missing": []}, "ubiquitous": {"missing": []}, "lesion_genes": {"missing": []}}


def test_sample_images(client, tiny_demo, tmp_path, monkeypatch):
    p = new_project(client, **demo_config(tiny_demo))
    base = f"/api/projects/{p['id']}/samples/CTRL-1"
    image = client.get(f"{base}/image")
    assert image.status_code == 200 and image.headers["content-type"] == "image/png"
    with Image.open(io.BytesIO(image.content)) as im:
        size = im.size
    meta = client.get(f"{base}/image-meta").json()
    assert (meta["width"], meta["height"]) == size and meta["frame"] == "hires"
    scalef = meta["scalefactors"]["tissue_hires_scalef"]
    assert scalef > 0 and meta["microns_per_pixel"] == meta["scalefactors"]["microns_per_pixel"]

    spatial = Path(tiny_demo["samples"][0]["path"]) / "square_016um" / "spatial"
    positions = pd.read_parquet(spatial / "tissue_positions.parquet")
    inside = positions[positions["in_tissue"] == 1]
    preview = client.get(f"{base}/bins-preview").json()
    assert preview["n_total"] == len(inside) == len(preview["x"]) == len(preview["y"])
    assert preview["x"][0] == pytest.approx(inside["pxl_col_in_fullres"].iloc[0] * scalef, abs=0.051)
    assert preview["y"][0] == pytest.approx(inside["pxl_row_in_fullres"].iloc[0] * scalef, abs=0.051)
    assert 0 <= min(preview["x"]) and max(preview["x"]) <= size[0]
    monkeypatch.setattr(server, "PREVIEW_MAX_BINS", 100)
    small = client.get(f"{base}/bins-preview").json()
    assert len(small["x"]) == len(small["y"]) == 100 and small["n_total"] == len(inside)
    assert small["x"][0] == preview["x"][0] and small["x"][-1] == preview["x"][-1]

    for suffix in ("image", "image-meta", "bins-preview"):
        assert client.get(f"/api/projects/{p['id']}/samples/NOPE/{suffix}").status_code == 404

    # a copy of the sample without the 16 um image and positions: the image of 8 um is used, no preview
    copy = tmp_path / "copy"
    shutil.copytree(Path(tiny_demo["samples"][0]["path"]).parent, copy / "outs")  # .../outs/binned_outputs -> outs
    (copy / "outs" / "binned_outputs" / "square_016um" / "spatial" / "tissue_hires_image.png").unlink()
    (copy / "outs" / "binned_outputs" / "square_016um" / "spatial" / "tissue_positions.parquet").unlink()
    q = new_project(client, samples=[{"id": "C", "label": "copy", "path": str(copy)}], reference_sample="C")
    assert client.get(f"/api/projects/{q['id']}/samples/C/image").content == image.content
    assert client.get(f"/api/projects/{q['id']}/samples/C/bins-preview").status_code == 404
    report = client.post(f"/api/projects/{q['id']}/validate").json()
    assert report["ok"] is False and report["samples"][0]["has_positions"] is False
    assert report["samples"][0]["has_hires_image"] is True


# ============================================================================
# A FINISHED RUN AND ITS RESULTS
# ============================================================================

@pytest.fixture(scope="module")
def finished(tmp_path_factory, tiny_demo):
    """One succeeded run of a project with a drawn region: (client, project, run record)."""
    app = create_app(tmp_path_factory.mktemp("ws"), worker_module=FAKE)
    with TestClient(app, base_url=LOCAL) as c:
        p = new_project(c, "Lifecycle", **demo_config(tiny_demo, roi={"use_drawn": True}))
        assert c.put(f"/api/projects/{p['id']}/roi/CTRL-1", json=SQUARE).status_code == 200
        run = start_run(c, p["id"], label="first", deterministic=False, steps={"moran": False})
        assert run["status"] == "queued" and run["label"] == "first" and run["config"]["deterministic"] is False
        yield c, p, wait_for(c, run["id"], {"succeeded", "failed"})


def test_run_record(finished):
    c, p, rec = finished
    assert rec["status"] == "succeeded", rec["error"]
    assert rec["project_id"] == p["id"] and rec["project_name"] == "Lifecycle" and rec["label"] == "first"
    assert rec["started"] and rec["finished"] and rec["duration_s"] >= 0 and isinstance(rec["pid"], int)
    assert rec["error"] is None and rec["traceback"] is None and rec["warnings"] == [fake_worker.WARNING]
    assert [(s["label"], s["status"]) for s in rec["steps"]] == [
        (label, "skipped" if label == fake_worker.STEPS[3] else "done") for label in fake_worker.STEPS]
    assert rec["steps"][3]["reason"] == "disabled" and all(s["seconds"] is not None for s in rec["steps"][:3])
    assert rec["progress"] == {"index": 5, "total": 5, "label": "", "fraction": 1.0}

    cfg = rec["config"]
    assert cfg["deterministic"] is False and cfg["steps"]["moran"] is False and cfg["steps"]["leiden"] is True
    storage = c.app.state.storage
    copy = storage.run_dir(rec["id"]) / "rois" / "CTRL-1.geojson"
    assert cfg["roi"]["masks"] == {"CTRL-1": str(copy)}
    assert json.loads(copy.read_text(encoding="utf-8")) == c.get(f"/api/projects/{p['id']}/roi/CTRL-1").json()
    assert json.loads((storage.run_dir(rec["id"]) / "config.json").read_text(encoding="utf-8")) == cfg

    listed = c.get("/api/runs", params={"project": p["id"]}).json()
    assert [r["id"] for r in listed] == [rec["id"]] and "config" not in listed[0]
    assert listed[0]["status"] == "succeeded" and c.get("/api/runs", params={"project": "p-00000000"}).json() == []
    last = c.get("/api/projects").json()[0]["last_run"]
    assert last == {"id": rec["id"], "status": "succeeded", "created": rec["created"]}
    health = c.get("/api/health").json()
    assert health["active_run"] is None and health["queued"] == 0
    assert c.get("/api/runs/r-20000101-000000-abcd").status_code == 404
    assert c.get("/api/runs/whatever").status_code == 404


def test_run_log(finished):
    c, _, rec = finished
    text, responses = full_log(c, rec["id"])
    assert len(responses) == 1
    raw = (c.app.state.storage.log_path(rec["id"])).read_bytes()
    assert text == raw.decode("utf-8") and responses[0]["offset"] == len(raw)
    assert "@@EVENT" not in text and "16 µm bins" in text and "a library writes to stderr" in text
    assert text.startswith("--- started") and text.rstrip().endswith("--- succeeded (worker exit code 0)")
    assert f"WARNING: {fake_worker.WARNING}" in text
    tail = c.get(f"/api/runs/{rec['id']}/log", params={"offset": len(raw)}).json()
    assert tail == {"offset": len(raw), "text": "", "complete": True}
    part = c.get(f"/api/runs/{rec['id']}/log", params={"offset": 10}).json()
    assert part["text"] == raw[10:].decode("utf-8") and part["complete"] is True
    beyond = c.get(f"/api/runs/{rec['id']}/log", params={"offset": len(raw) + 999}).json()
    assert beyond == {"offset": len(raw), "text": "", "complete": True}
    assert c.get(f"/api/runs/{rec['id']}/log", params={"offset": "x"}).status_code == 400


def test_summary_and_manifest(finished):
    c, _, rec = finished
    rid = rec["id"]
    summary = c.get(f"/api/runs/{rid}/summary").json()
    assert summary["version"] == 1 and summary["bin_size"] == "016um"
    assert summary["metrics"][0]["value"] is None  # NaN in summary.json
    m = c.get(f"/api/runs/{rid}/manifest").json()
    results = c.app.state.storage.results_dir(rid)
    tables = {t["file"]: t for t in m["tables"]}
    assert list(tables) == [fake_worker.POSITIVITY_TABLE, "per_bin_signature_calls_16um.csv.gz",
                            "unlisted_table.csv"]
    assert tables[fake_worker.POSITIVITY_TABLE] == {
        "file": fake_worker.POSITIVITY_TABLE, "title": "Positivity, native depth (ROI)",
        "description": "Figures 4g and 6a, left panels", "rows": 6,
        "size": (results / fake_worker.POSITIVITY_TABLE).stat().st_size}
    unlisted = tables["unlisted_table.csv"]
    assert unlisted["title"] == "unlisted_table.csv" and unlisted["rows"] is None and unlisted["description"] == ""
    assert m["objects"] == [{"file": "joint_16um.h5ad", "title": "Joint AnnData", "description": "all sections",
                             "size": (results / "joint_16um.h5ad").stat().st_size}]
    assert {o["file"] for o in m["other"]} == {"summary.json", "config.json", "notes.txt"}
    assert all(set(o) == {"file", "size"} for o in m["other"])
    fig, failed = m["figures"]
    assert fig["name"] == fake_worker.FIGURE and fig["pdf_url"] is None
    assert fig["png_url"] == f"/api/runs/{rid}/files/figures/{fake_worker.FIGURE}.png"
    assert failed["png_url"] is None and failed["error"]
    png = c.get(fig["png_url"])
    assert png.status_code == 200 and png.headers["content-type"] == "image/png"
    assert png.headers["content-disposition"].startswith("inline")


def test_result_files(finished):
    c, _, rec = finished
    base = f"/api/runs/{rec['id']}/files"
    r = c.get(f"{base}/{fake_worker.POSITIVITY_TABLE}", params={"download": 1})
    assert r.status_code == 200 and r.headers["content-type"].startswith("text/csv")
    assert r.headers["content-disposition"] == f'attachment; filename="{fake_worker.POSITIVITY_TABLE}"'
    gz = c.get(f"{base}/per_bin_signature_calls_16um.csv.gz")
    assert gz.headers["content-type"] == "application/gzip"
    assert gzip.decompress(gz.content).decode("utf-8").startswith("barcode,sample_id")
    for bad in ("..%2Frun.json", "figures/..%2F..%2Frun.json", "%2E%2E%2Fconfig.json", "%5C..%5Crun.json",
                "C:%2FWindows%2Fwin.ini", "figures%2F%2Fx.png"):
        r = c.get(f"{base}/{bad}")
        assert r.status_code == 400, bad
    assert c.get(f"{base}/missing.csv").status_code == 404
    assert c.get(f"{base}/figures").status_code == 404  # a folder
    assert c.get("/api/runs/r-20000101-000000-abcd/files/notes.txt").status_code == 404


def test_table(finished):
    c, _, rec = finished
    url = f"/api/runs/{rec['id']}/table"
    t = c.get(url, params={"file": fake_worker.POSITIVITY_TABLE}).json()
    assert t["file"] == fake_worker.POSITIVITY_TABLE and t["total"] == 6 and t["offset"] == 0 and t["limit"] == 100
    assert t["columns"][:3] == ["signature", "sample", "roi"]
    assert t["dtypes"]["signature"] == "text" and t["dtypes"]["pct_positive"] == "number"
    pct = t["columns"].index("pct_positive")
    assert [row[pct] for row in t["rows"]].count(None) == 1  # NaN -> null

    def column(params, name="pct_positive"):
        data = c.get(url, params={"file": fake_worker.POSITIVITY_TABLE, **params}).json()
        return [row[data["columns"].index(name)] for row in data["rows"]], data

    up, _ = column({"sort": "pct_positive"})
    assert up[:-1] == sorted(up[:-1]) and up[-1] is None
    down, _ = column({"sort": "pct_positive", "desc": 1})
    assert down[:-1] == sorted(down[:-1], reverse=True) and down[-1] is None
    names, data = column({"q": "OSTEOC"}, "signature")
    assert data["total"] == 4 and set(names) == {"osteoclast", "osteocyte"}
    names, data = column({"q": "tum-1", "sort": "n_positive", "desc": "1"}, "n_positive")
    assert data["total"] == 3 and names == [61, 51, 41]
    _, data = column({"offset": 4, "limit": 10})
    assert data["total"] == 6 and len(data["rows"]) == 2 and data["offset"] == 4
    _, data = column({"limit": 100000})
    assert data["limit"] == server.TABLE_MAX_LIMIT
    assert column({"q": "no such text"})[1]["total"] == 0

    per_bin = c.get(url, params={"file": "per_bin_signature_calls_16um.csv.gz", "sort": "leiden_domain",
                                 "limit": 1000}).json()
    assert per_bin["total"] == 80 and per_bin["dtypes"]["in_region_of_interest"] == "bool"
    domains = [row[per_bin["columns"].index("leiden_domain")] for row in per_bin["rows"]]
    order = [d for d in dict.fromkeys(domains)]
    assert order == [str(i) for i in range(11)] + ["unassigned"]  # numbers in numeric order, then text

    assert c.get(url, params={"file": fake_worker.POSITIVITY_TABLE, "sort": "nope"}).status_code == 400
    assert c.get(url, params={"file": "summary.json"}).status_code == 400
    assert c.get(url, params={"file": "missing.csv"}).status_code == 404
    assert c.get(url, params={"file": "../run.json"}).status_code == 400
    assert c.get(url).status_code == 400


def test_bins(finished):
    c, _, rec = finished
    rid = rec["id"]
    samples = c.get(f"/api/runs/{rid}/bins/samples").json()
    assert samples == [{"id": "CTRL-1", "label": "control", "n": 40}, {"id": "TUM-1", "label": "tumor-bearing",
                                                                        "n": 40}]
    b = c.get(f"/api/runs/{rid}/bins", params={"sample": "TUM-1"}).json()
    assert b["sample"] == "TUM-1" and b["n"] == 40 == len(b["x"]) == len(b["y"])
    assert b["x"][:3] == [0.0, 16.0, 32.0] and b["y"][8] == 16.1  # rounded to 0.1 um
    f = b["fields"]
    assert f["tissue_class"]["type"] == "category"
    assert f["tissue_class"]["categories"] == ["bone lineage", "hematopoietic marrow", "removed", "unassigned"]
    assert f["assigned_signature"]["categories"] == ["osteoclast", "osteoblast", "granulocyte", "skeletal_muscle",
                                                     "ambiguous"]
    assert f["leiden_domain"]["categories"] == [str(i) for i in range(11)] + ["unassigned"]
    codes = f["assigned_signature"]["codes"]
    assert codes[:5] == [0, 1, 2, 3, 4] and len(codes) == 40
    assert f["leiden_domain"]["codes"][:3] == [11, 1, 2]  # bin 0 is "unassigned"
    assert f["in_region_of_interest"] == {"type": "bool", "values": [0 if k % 3 == 0 else 1 for k in range(40)]}
    assert f["positive_osteoclast"]["values"][:5] == [1, 0, 0, 0, 0]
    assert f["positive_granulocyte"]["type"] == "bool"
    assert f["total_umi"] == {"type": "number", "values": [10 + k for k in range(40)]}
    assert set(f) >= {"positive_osteoblast", "positive_osteoblast_pooled_depthmatched",
                      "positive_osteoclast_pooled_depthmatched"}
    assert b["colors"]["tissue_class"] == {"bone lineage": "#984EA3", "chondrocyte": "#33A02C",
                                           "hematopoietic marrow": "#FF7F00", "stroma": "#E6AB02",
                                           "removed": "#b8b8b8", "unassigned": "#dedede"}
    sig_colors = b["colors"]["assigned_signature"]
    assert sig_colors["granulocyte"] == "#FF7F00" and sig_colors["skeletal_muscle"] == "#8C8C8C"
    assert sig_colors["ambiguous"] == "#dedede" and "neutrophil" not in sig_colors
    assert c.get(f"/api/runs/{rid}/bins", params={"sample": "NOPE"}).status_code == 404
    assert c.get(f"/api/runs/{rid}/bins").status_code == 400


def test_rerun_and_delete(finished):
    c, p, rec = finished
    r = c.post(f"/api/runs/{rec['id']}/rerun")
    assert r.status_code == 200, r.text
    again = r.json()
    assert again["id"] != rec["id"] and again["status"] == "queued" and again["label"] == "first"
    storage = c.app.state.storage
    copy = storage.run_dir(again["id"]) / "rois" / "CTRL-1.geojson"
    assert again["config"]["roi"]["masks"] == {"CTRL-1": str(copy)}
    strip = lambda cfg: {**cfg, "roi": {**cfg["roi"], "masks": {}}}  # noqa: E731
    assert strip(again["config"]) == strip(rec["config"])
    assert copy.read_bytes() == (storage.run_dir(rec["id"]) / "rois" / "CTRL-1.geojson").read_bytes()
    assert wait_for(c, again["id"], {"succeeded", "failed"})["status"] == "succeeded"
    assert [x["id"] for x in c.get("/api/runs").json()] == [again["id"], rec["id"]]

    assert c.delete(f"/api/runs/{rec['id']}").json() == {"ok": True}
    assert c.get(f"/api/runs/{rec['id']}").status_code == 404
    assert c.get(f"/api/runs/{rec['id']}/summary").status_code == 404
    assert not storage.run_dir(rec["id"]).exists()
    assert c.delete(f"/api/runs/{rec['id']}").status_code == 404
    # deleting the project keeps its runs
    assert c.delete(f"/api/projects/{p['id']}").json() == {"ok": True}
    kept = c.get(f"/api/runs/{again['id']}").json()
    assert kept["status"] == "succeeded" and kept["project_name"] == "Lifecycle"
    assert c.post(f"/api/runs/{again['id']}/rerun").json()["project_name"] == "Lifecycle"


# ============================================================================
# RUN LIFECYCLE
# ============================================================================

def test_failed_run(client):
    p = new_project(client, fake_fail="synthetic failure")
    rec = wait_for(client, start_run(client, p["id"])["id"], {"succeeded", "failed"})
    assert rec["status"] == "failed" and rec["error"] == "synthetic failure"
    assert "RuntimeError: synthetic failure" in rec["traceback"]
    assert [s["status"] for s in rec["steps"]] == ["done", "done", "failed", "pending", "pending"]
    assert rec["steps"][2]["reason"] == "synthetic failure" and rec["progress"]["fraction"] == 0.6
    text, _ = full_log(client, rec["id"])
    assert "Traceback (most recent call last)" in text and text.rstrip().endswith("--- failed (worker exit code 1)")
    assert client.get(f"/api/runs/{rec['id']}/summary").status_code == 404
    assert client.get(f"/api/runs/{rec['id']}/manifest").json() == {"tables": [], "figures": [], "objects": [],
                                                                     "other": []}
    assert client.get(f"/api/runs/{rec['id']}/bins/samples").status_code == 404


def test_cancel_queue_and_delete_rules(client):
    p = new_project(client, fake_sleep=60)
    first = start_run(client, p["id"], label="A")
    second = start_run(client, p["id"], label="B")
    running = wait_for(client, first["id"], {"running"}, label=fake_worker.STEPS[1])
    assert running["started"] and running["pid"] and running["progress"]["index"] == 2
    assert running["steps"][0]["status"] == "done" and running["steps"][1]["status"] == "running"
    assert client.get(f"/api/runs/{second['id']}").json()["status"] == "queued"
    health = client.get("/api/health").json()
    assert health["active_run"] == first["id"] and health["queued"] == 1
    log = client.get(f"/api/runs/{first['id']}/log").json()
    assert log["complete"] is False and "=== [2/5]" in log["text"]
    for rid in (first["id"], second["id"]):
        r = client.delete(f"/api/runs/{rid}")
        assert r.status_code == 409 and "cancel" in r.json()["detail"]

    queued = client.post(f"/api/runs/{second['id']}/cancel").json()
    assert queued["status"] == "cancelled" and queued["started"] is None and queued["finished"]
    t0 = time.monotonic()
    cancelled = client.post(f"/api/runs/{first['id']}/cancel").json()
    assert time.monotonic() - t0 < 15
    assert cancelled["status"] == "cancelled" and cancelled["config"]["fake_sleep"] == 60
    assert [s["status"] for s in cancelled["steps"]] == ["done", "skipped", "pending", "pending", "pending"]
    assert cancelled["steps"][1]["reason"] == "cancelled"
    text, _ = full_log(client, first["id"])
    assert text.rstrip().splitlines()[-1].startswith("--- cancelled")
    assert full_log(client, second["id"])[0] == "--- cancelled before it started\n"
    assert client.post(f"/api/runs/{first['id']}/cancel").json()["status"] == "cancelled"  # idempotent
    health = client.get("/api/health").json()
    assert health["active_run"] is None and health["queued"] == 0
    for rid in (first["id"], second["id"]):
        assert client.delete(f"/api/runs/{rid}").json() == {"ok": True}
    assert client.get("/api/runs").json() == []


def test_log_paging(client):
    p = new_project(client, fake_log_kb=700)
    rec = wait_for(client, start_run(client, p["id"])["id"], {"succeeded", "failed"})
    assert rec["status"] == "succeeded", rec["error"]
    text, responses = full_log(client, rec["id"])
    raw = client.app.state.storage.log_path(rec["id"]).read_bytes()
    assert len(responses) >= 3 and text == raw.decode("utf-8")
    sizes = [len(r["text"].encode("utf-8")) for r in responses]
    assert max(sizes) <= server.LOG_CHUNK and [r["complete"] for r in responses][-2:] == [False, True]


def test_run_refused(client, monkeypatch):
    bad = new_project(client, qc={"min_genes": -1})
    r = client.post(f"/api/projects/{bad['id']}/runs")
    assert r.status_code == 409 and "qc.min_genes" in r.json()["detail"]
    good = new_project(client)
    for body in ({"steps": {"no_such_step": True}}, {"steps": {"moran": "yes"}}, {"deterministic": "yes"},
                 {"label": 5}, {"steps": [1]}):
        assert client.post(f"/api/projects/{good['id']}/runs", json=body).status_code == 400, body
    versions = server.package_versions()
    monkeypatch.setattr(server, "package_versions", lambda: {**versions, "scanpy": None, "leidenalg": None})
    r = client.post(f"/api/projects/{good['id']}/runs")
    assert r.status_code == 409 and "scanpy, leidenalg" in r.json()["detail"]
    assert client.get("/api/health").json()["missing_required"] == ["scanpy", "leidenalg"]
    monkeypatch.undo()
    # a region drawn on the hires image while roi.mask_frame says lowres
    frame = new_project(client, roi={"use_drawn": True, "mask_frame": "lowres"})
    client.put(f"/api/projects/{frame['id']}/roi/LL-2-1", json=SQUARE)
    r = client.post(f"/api/projects/{frame['id']}/runs")
    assert r.status_code == 409 and "mask_frame" in r.json()["detail"]
    assert any("mask_frame" in e for e in client.post(f"/api/projects/{frame['id']}/validate").json()["errors"])
    assert client.post("/api/projects/p-00000000/runs").status_code == 404
    assert client.get("/api/runs").json() == []


def test_restart_recovery(tmp_path):
    ws = tmp_path / "ws"
    storage = Storage(ws)
    project = storage.create_project("Restart")
    cfg = pc.merge_with_defaults(project["config"])
    left_running = storage.create_run(project_id=project["id"], project_name="Restart", label="was running",
                                      config=cfg, roi_files={})
    left_running.update(status="running", started=left_running["created"], pid=999999, steps=[
        {"label": fake_worker.STEPS[0], "status": "done", "seconds": 1.0, "reason": None},
        {"label": fake_worker.STEPS[1], "status": "running", "seconds": None, "reason": None}])
    storage.save_run(left_running)
    left_queued = storage.create_run(project_id=project["id"], project_name="Restart", label="was queued",
                                     config=cfg, roi_files={})
    with TestClient(create_app(ws, worker_module=FAKE), base_url=LOCAL) as c:
        rec = c.get(f"/api/runs/{left_running['id']}").json()
        assert rec["status"] == "interrupted" and rec["pid"] is None and "rerun" in rec["error"]
        assert rec["steps"][1] == {"label": fake_worker.STEPS[1], "status": "failed", "seconds": None,
                                   "reason": "interrupted"}
        assert "interrupted" in full_log(c, left_running["id"])[0]
        assert wait_for(c, left_queued["id"], {"succeeded", "failed"})["status"] == "succeeded"
        # stopping the server interrupts the active run
        c.put(f"/api/projects/{project['id']}", json={"config": {"fake_sleep": 60}})
        sleeper = start_run(c, project["id"])
        wait_for(c, sleeper["id"], {"running"}, label=fake_worker.STEPS[1])
    record = Storage(ws).get_run(sleeper["id"])
    assert record["status"] == "interrupted" and record["finished"]
    with TestClient(create_app(ws, worker_module=FAKE), base_url=LOCAL) as c:
        assert c.get(f"/api/runs/{sleeper['id']}").json()["status"] == "interrupted"
        assert c.post(f"/api/runs/{sleeper['id']}/rerun").status_code == 200


def test_other_websites_cannot_use_the_api(tmp_path):
    with TestClient(create_app(tmp_path / "ws", worker_module=FAKE), base_url=LOCAL) as c:
        assert c.get("/api/health").status_code == 200
        assert c.get("/api/health", headers={"Origin": LOCAL}).status_code == 200
        assert c.get("/api/health", headers={"Sec-Fetch-Site": "same-origin"}).status_code == 200
        # DNS rebinding: a foreign host name that resolves to 127.0.0.1
        assert c.get("/api/fs/list", headers={"Host": "attacker.example:8765"}).status_code == 403
        assert c.get("/", headers={"Host": "attacker.example:8765"}).status_code == 403
        # cross-site requests
        assert c.get("/api/projects", headers={"Origin": "http://attacker.example"}).status_code == 403
        assert c.get("/api/projects", headers={"Sec-Fetch-Site": "cross-site"}).status_code == 403
        r = c.post("/api/projects", content='{"name": "x"}', headers={"Content-Type": "text/plain;charset=UTF-8"})
        assert r.status_code == 403 and c.get("/api/projects").json() == []
        assert c.post("/api/projects", json={"name": "ok"}).status_code == 200
        assert c.get("/api/health", headers={"Host": "localhost:8765"}).status_code == 200
        assert c.get("/api/health", headers={"Host": "[::1]:8765"}).status_code == 200


def test_case_only_rename_keeps_the_region(client):
    pid = new_project(client, "Case")["id"]
    fc = {"type": "FeatureCollection", "features": [{"type": "Feature", "properties": {"frame": "hires"}, "geometry": {
        "type": "Polygon", "coordinates": [[[0, 0], [10, 0], [10, 10], [0, 0]]]}}]}
    assert client.put(f"/api/projects/{pid}/roi/ctrl-1", json=fc).status_code == 200
    assert client.put(f"/api/projects/{pid}/roi/CTRL-1", json=fc).status_code == 200
    assert client.delete(f"/api/projects/{pid}/roi/ctrl-1").status_code == 200
    assert client.get(f"/api/projects/{pid}/roi/CTRL-1").status_code == 200
    assert client.get(f"/api/projects/{pid}").json()["rois"] == ["CTRL-1"]


def test_ids_with_a_trailing_newline_are_rejected(client):
    pid = new_project(client, "Newline")["id"]
    assert client.get(f"/api/projects/{pid}%0A").status_code == 404
    assert client.get("/api/runs/r-20261028-000000-abcd%0A").status_code == 404


def test_token(tmp_path):
    with pytest.raises(ValueError):
        create_app(tmp_path / "bad", token="has space")
    with TestClient(create_app(tmp_path / "ws", token=TOKEN, worker_module=FAKE)) as c:
        for path in ("/api/health", "/", "/static/js/core.js", "/api/projects"):
            assert c.get(path).status_code == 401, path
        r = c.get("/api/health")
        assert "token" in r.json()["detail"]
        assert "Access token required" in c.get("/").text
        assert c.get("/api/health", headers={"X-Token": "wrong"}).status_code == 401
        assert c.get("/api/health", headers={"X-Token": TOKEN}).status_code == 200
        assert c.get("/?token=wrong", follow_redirects=False).status_code == 401
        r = c.get(f"/?token={TOKEN}", follow_redirects=False)
        assert r.status_code == 303 and r.headers["location"] == "/"
        cookie = r.headers["set-cookie"].lower()
        assert cookie.startswith(f"bss_token={TOKEN.lower()}") and "httponly" in cookie and "samesite=lax" in cookie
        assert c.get("/api/health").status_code == 200 and c.get("/").status_code == 200
        assert c.post("/api/projects", json={"name": "with cookie"}).status_code == 200
        c.cookies.clear()
        assert c.get("/api/projects").status_code == 401


def test_demo(client, tiny_demo, monkeypatch):
    import tools.make_demo_data as generator

    calls = []

    def fake_generate(out_dir, size="small", seed=0, force=False):
        calls.append((Path(out_dir), size))
        return [dict(s) for s in tiny_demo["samples"]]

    monkeypatch.setattr(generator, "generate", fake_generate)
    r = client.post("/api/demo")
    assert r.status_code == 200, r.text
    project = r.json()["project"]
    storage = client.app.state.storage
    assert calls == [(storage.demo_dir, "small")]
    cfg = project["config"]
    assert project["name"] == DEMO_NAME and cfg["name"] == DEMO_NAME
    assert cfg["data_root"] == str(storage.demo_dir) and cfg["samples"] == tiny_demo["samples"]
    assert cfg["reference_sample"] == "CTRL-1" and project["validation"]["errors"] == []
    assert cfg["scoring"] == {**pc.default_config()["scoring"], "n_perm": 20, "n_perm_resampled": 10,
                              "n_depth_reps": 2, "null_n_sub": 5000}
    assert cfg["clustering"]["stability_seeds"] == [0, 1] and cfg["moran"]["n_genes"] == 500
    assert cfg["uncoupling"]["n_perm"] == 99
    again = client.post("/api/demo").json()["project"]
    assert again["id"] == project["id"] and len(calls) == 2
    assert [p["name"] for p in client.get("/api/projects").json()] == [DEMO_NAME]

    def broken(out_dir, size="small", seed=0, force=False):
        raise OSError("disk full")

    monkeypatch.setattr(generator, "generate", broken)
    r = client.post("/api/demo")
    assert r.status_code == 500 and "disk full" in r.json()["detail"]
