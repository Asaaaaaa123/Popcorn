"""
tests/test_results_export.py
============================
Download of a run's results as one zip file (GET /api/runs/{rid}/export/results?part=all|figures|tables), with the
fake worker of tests/test_api.py: which files each part holds, that AnnData objects are left out, that the run's
configuration and log come with "all", and the refusals.
"""
from __future__ import annotations

import io
import tempfile
import zipfile

from tests.test_api import FAKE, LOCAL, TestClient, client, new_project, start_run, wait_for  # noqa: F401
from webapp import server


def finished_run(client) -> str:
    rid = start_run(client, new_project(client)["id"])["id"]
    assert wait_for(client, rid, {"succeeded", "failed"})["status"] == "succeeded"
    return rid


def get_zip(client, rid: str, part: str | None = None) -> tuple[object, zipfile.ZipFile]:
    r = client.get(f"/api/runs/{rid}/export/results", params={"part": part} if part else None)
    assert r.status_code == 200, r.text
    return r, zipfile.ZipFile(io.BytesIO(r.content))


def test_each_part_holds_its_files(client):
    rid = finished_run(client)
    on_disk = server.scan_results(client.app.state.storage.results_dir(rid))

    r, zf = get_zip(client, rid)
    assert r.headers["content-type"] == "application/zip" and r.headers["cache-control"] == "no-store"
    assert r.headers["content-disposition"] == f'attachment; filename="{rid}_results.zip"'
    names = zf.namelist()
    assert sorted(n for n in names if not n.startswith("run_")) == sorted(n for n in on_disk if not n.endswith(".h5ad"))
    assert not any(n.endswith(".h5ad") for n in names)
    assert {"run_config.json", "run_log.txt"} <= set(names)
    assert zf.read("summary.json") == (client.app.state.storage.results_dir(rid) / "summary.json").read_bytes()

    r, zf = get_zip(client, rid, "figures")
    assert r.headers["content-disposition"] == f'attachment; filename="{rid}_figures.zip"'
    figures = [n for n in on_disk if n.startswith("figures/") and n.endswith((".png", ".pdf"))]
    assert figures and sorted(zf.namelist()) == sorted(figures)
    assert all(i.compress_type == zipfile.ZIP_STORED for i in zf.infolist() if i.filename.endswith(".png"))

    r, zf = get_zip(client, rid, "tables")
    assert r.headers["content-disposition"] == f'attachment; filename="{rid}_tables.zip"'
    tables = [n for n in on_disk if server.is_table_file(n)]
    assert tables and sorted(zf.namelist()) == sorted(tables)
    assert any(n.endswith(".csv.gz") for n in zf.namelist())


def test_refusals_and_no_temporary_file_left(client, tmp_path, monkeypatch):
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path / "tmp"))
    (tmp_path / "tmp").mkdir()
    r = client.get("/api/runs/r-20000101-000000-0000/export/results")
    assert r.status_code == 404 and "not found" in r.json()["detail"]

    rid = finished_run(client)
    r = client.get(f"/api/runs/{rid}/export/results", params={"part": "objects"})
    assert r.status_code == 400 and "all, figures, tables" in r.json()["detail"]

    for path in (client.app.state.storage.results_dir(rid) / "figures").iterdir():
        path.unlink()
    r = client.get(f"/api/runs/{rid}/export/results", params={"part": "figures"})
    assert r.status_code == 404 and "has no figures" in r.json()["detail"]

    failed = wait_for(client, start_run(client, new_project(client, fake_fail="boom")["id"])["id"], {"failed"})
    r = client.get(f"/api/runs/{failed['id']}/export/results", params={"part": "tables"})
    assert r.status_code == 404 and "(failed) before writing any" in r.json()["detail"]

    assert client.get(f"/api/runs/{rid}/export/results").status_code == 200
    assert list((tmp_path / "tmp").iterdir()) == []  # the temporary zip goes once sent, and after a refusal

    running = start_run(client, new_project(client, fake_sleep=60)["id"])
    wait_for(client, running["id"], {"running"})
    r = client.get(f"/api/runs/{running['id']}/export/results")
    assert r.status_code == 409 and "when it has finished" in r.json()["detail"]
    client.post(f"/api/runs/{running['id']}/cancel")
