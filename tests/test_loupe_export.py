"""
tests/test_loupe_export.py
==========================
Export of a run for 10x Genomics Loupe Browser (GET /api/runs/{rid}/export/loupe) on the tiny synthetic dataset: one
category CSV per section and bin size whose barcode column is exactly the barcodes of that bin size's
filtered_feature_bc_matrix.h5 (the bins of the section's .cloupe file), a text label in every cell ("not analyzed"
for a bin the run has no calls for), 8 um bins carrying the calls of the 16 um bin that holds them, and README.txt.
The per-bin table of the run is written by the tests from the raw barcodes (no scanpy needed).
"""
from __future__ import annotations

import gzip
import io
import json
import shutil
import tempfile
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy.spatial import cKDTree

import pipeline_config as pc
from tests.test_api import FAKE, LOCAL, TestClient, client, demo_config, new_project, start_run, wait_for  # noqa: F401
from webapp import server
from webapp.server import create_app

HEADER = ["barcode", "Tissue class", "Assigned signature", "Leiden domain", "Region of interest",
          "Positive osteoclast", "Positive osteoblast", "Positive osteoblast pooled depthmatched", "In lesion",
          "Custom call"]
QC = "8 um QC"  # the 8 um files' last column: the 8 um bin's own quality control
CLASSES = [("osteoclast", "bone lineage"), ("granulocyte", "hematopoietic marrow"), ("skeletal_muscle", "removed"),
           ("ambiguous", "unassigned")]


# ============================================================================
# HELPERS
# ============================================================================

def matrix_barcodes(sample: dict, bin_size: str) -> list[str]:
    with h5py.File(Path(sample["path"]) / f"square_{bin_size}" / "filtered_feature_bc_matrix.h5", "r") as f:
        return [b.decode() for b in f["matrix/barcodes"][()]]


def per_bin_table(samples: list[dict]) -> pd.DataFrame:
    """A per-bin table as the analysis writes it for the 16 um bins of the sections (barcode = matrix barcode
    without -1, plus -<sample id>), every seventh bin left out (as if removed by quality control), with every kind of
    column: calls, a yes/no column with a missing value, numbers, a text call with a comma and a quote in its labels,
    a text column of too many values and an empty one."""
    parts = []
    for s in samples:
        raw = matrix_barcodes(s, "016um")
        kept = [b for i, b in enumerate(raw) if i % 7 != 3]
        k = np.arange(len(kept))
        positive = (k % 5 == 0).astype(object)
        positive[1] = np.nan
        parts.append(pd.DataFrame({
            "barcode": [f"{b[:-2]}-{s['id']}" for b in kept], "sample_id": s["id"], "section": s.get("label") or s["id"],
            "x_um": k * 16.0, "y_um": k * 8.0, "in_region_of_interest": k % 3 != 0, "total_umi": 10 + k,
            "assigned_signature": [CLASSES[i % 4][0] for i in k], "tissue_class": [CLASSES[i % 4][1] for i in k],
            "leiden_domain": [("unassigned" if i % 13 == 0 else str(i % 11)) for i in k],
            "positive_osteoclast": k % 4 == 0, "positive_osteoblast": positive,
            "positive_osteoblast_pooled_depthmatched": k % 6 == 0,
            "positive_fraction_osteoblast_pooled_depthmatched_draws": (k % 3) / 2,
            "in_lesion": k % 2 == 0, "custom_call": [("a, b" if i % 2 else 'say "x"') for i in k],
            "free_text": [f"note {s['id']} {i}" for i in k], "empty_col": np.nan,
        }))
    return pd.concat(parts, ignore_index=True)


def write_per_bin(results: Path, table: pd.DataFrame) -> None:
    with gzip.open(results / "per_bin_signature_calls_16um.csv.gz", "wt", newline="", encoding="utf-8") as fh:
        table.to_csv(fh, index=False)


def read_csv(zf: zipfile.ZipFile, name: str) -> pd.DataFrame:
    return pd.read_csv(zf.open(name), dtype=str, keep_default_na=False)


def expected_labels(table: pd.DataFrame) -> pd.DataFrame:
    """The labels the export writes for the rows of a per-bin table, by matrix barcode."""
    yes_no = lambda s, yes, no: s.map(lambda v: "NA" if not isinstance(v, (bool, np.bool_)) else yes if v else no)
    out = pd.DataFrame({
        "Tissue class": table["tissue_class"], "Assigned signature": table["assigned_signature"],
        "Leiden domain": [v if v == "unassigned" else f"domain {v}" for v in table["leiden_domain"].astype(str)],
        "Region of interest": yes_no(table["in_region_of_interest"], "inside", "outside"),
        "Positive osteoclast": yes_no(table["positive_osteoclast"], "positive", "negative"),
        "Positive osteoblast": yes_no(table["positive_osteoblast"], "positive", "negative"),
        "Positive osteoblast pooled depthmatched": yes_no(table["positive_osteoblast_pooled_depthmatched"],
                                                          "positive", "negative"),
        "In lesion": yes_no(table["in_lesion"], "inside", "outside"),
        "Custom call": table["custom_call"].map({"a, b": "a; b", 'say "x"': "say 'x'"}),
    })
    out.index = [b[: -len(sid) - 1] + "-1" for b, sid in zip(table["barcode"], table["sample_id"])]
    return out


def qc_passes(sample: dict, qc: dict) -> list[str]:
    """"pass" / "fail" of every 8 um bin (matrix order) by the thresholds qc, computed here from the dense counts."""
    with h5py.File(Path(sample["path"]) / "square_008um" / "filtered_feature_bc_matrix.h5", "r") as f:
        m = f["matrix"]
        n_genes, n = m["shape"][()]
        X = np.zeros((n, n_genes), np.int64)
        indptr, indices, data = m["indptr"][()], m["indices"][()], m["data"][()]
        for i in range(n):
            X[i, indices[indptr[i]:indptr[i + 1]]] = data[indptr[i]:indptr[i + 1]]
        mt = np.array([g.decode().lower().startswith("mt-") for g in m["features/name"][()]])
    total, genes, mito = X.sum(axis=1), (X > 0).sum(axis=1), X[:, mt].sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        ok = (genes > qc["min_genes"]) & (total > qc["min_counts"]) & (100 * mito / total < qc["max_pct_mt"])
    return ["pass" if v else "fail" for v in ok]


def parent_barcode(barcode: str) -> str:
    _, size, row, col = barcode[:-2].split("_")
    return f"s_016um_{int(row) // 2:05d}_{int(col) // 2:05d}-1"


@pytest.fixture(scope="module")
def loupe_run(tmp_path_factory, tiny_demo):
    """A succeeded (fake) run of the tiny demo with the per-bin table of per_bin_table(), and its export."""
    app = create_app(tmp_path_factory.mktemp("ws"), worker_module=FAKE)
    with TestClient(app, base_url=LOCAL) as c:
        p = new_project(c, "Loupe", **demo_config(tiny_demo))
        rid = start_run(c, p["id"], label="for Loupe")["id"]
        assert wait_for(c, rid, {"succeeded", "failed"})["status"] == "succeeded"
        table = per_bin_table(tiny_demo["samples"])
        write_per_bin(app.state.storage.results_dir(rid), table)
        r = c.get(f"/api/runs/{rid}/export/loupe")
        assert r.status_code == 200, r.text
        yield {"client": c, "rid": rid, "table": table, "response": r, "zip": zipfile.ZipFile(io.BytesIO(r.content))}


# ============================================================================
# THE EXPORT OF A RUN
# ============================================================================

def test_zip_files_and_format(loupe_run):
    r, zf, rid = loupe_run["response"], loupe_run["zip"], loupe_run["rid"]
    assert r.headers["content-type"] == "application/zip" and r.headers["cache-control"] == "no-store"
    assert r.headers["content-disposition"] == f'attachment; filename="{rid}_loupe.zip"'
    assert zf.namelist() == ["CTRL-1_016um_categories.csv", "CTRL-1_008um_categories.csv",
                             "TUM-1_016um_categories.csv", "TUM-1_008um_categories.csv", "README.txt"]
    for name in zf.namelist()[:-1]:
        raw = zf.read(name)
        assert not raw.startswith(b"\xef\xbb\xbf") and b"\r" not in raw and b'"' not in raw and raw.endswith(b"\n")
        lines = raw.decode("utf-8").split("\n")[:-1]
        header = HEADER + ([QC] if "_008um_" in name else [])
        assert lines[0].split(",") == header
        for line in lines:
            cells = line.split(",")
            assert len(cells) == len(header) and all(cells), line  # never an empty cell


def test_barcodes_are_those_of_the_matrix(loupe_run, tiny_demo):
    zf = loupe_run["zip"]
    for s in tiny_demo["samples"]:
        for bin_size in ("016um", "008um"):
            assert read_csv(zf, f"{s['id']}_{bin_size}_categories.csv")["barcode"].tolist() == \
                matrix_barcodes(s, bin_size)


def test_labels_of_the_run_bins(loupe_run, tiny_demo):
    zf, table = loupe_run["zip"], loupe_run["table"]
    for s in tiny_demo["samples"]:
        got = read_csv(zf, f"{s['id']}_016um_categories.csv").set_index("barcode")
        expected = expected_labels(table[table["sample_id"] == s["id"]])
        analysed = got.index.isin(expected.index)
        assert analysed.sum() == len(expected) and (~analysed).sum() == len(got) // 7 + (len(got) % 7 > 3)
        assert got.columns.tolist() == expected.columns.tolist()
        assert got[analysed].to_numpy().tolist() == expected.loc[got.index[analysed]].to_numpy().tolist()
        assert (got[~analysed] == "not analyzed").all().all()
        assert (got["Positive osteoblast"] == "NA").sum() == 1  # the missing value of an analysed bin
    assert set(got["Custom call"]) == {"a; b", "say 'x'", "not analyzed"}


def test_8um_bins_take_the_calls_of_their_16um_bin(loupe_run, tiny_demo):
    zf = loupe_run["zip"]
    for s in tiny_demo["samples"]:
        coarse = read_csv(zf, f"{s['id']}_016um_categories.csv").set_index("barcode")
        fine = read_csv(zf, f"{s['id']}_008um_categories.csv")
        parents = [parent_barcode(b) for b in fine["barcode"]]
        expected = coarse.reindex(parents).fillna("not analyzed")
        assert fine.drop(columns=["barcode", QC]).to_numpy().tolist() == expected.to_numpy().tolist()
        assert (fine["Tissue class"] != "not analyzed").sum() > 0.7 * len(fine)
        # the 8 um bin's own quality control, by the run's thresholds, whatever its calls
        qc = pc.merge_with_defaults(demo_config(tiny_demo))["qc"]
        assert fine[QC].tolist() == qc_passes(s, qc)
        assert {"pass", "fail"} <= set(fine.loc[fine["Tissue class"] != "not analyzed", QC])
        # the parent by row // 2, column // 2 is the 16 um bin whose centre is nearest (the grids nest)
        base = Path(s["path"])
        p16 = pd.read_parquet(base / "square_016um" / "spatial" / "tissue_positions.parquet")
        p8 = pd.read_parquet(base / "square_008um" / "spatial" / "tissue_positions.parquet").set_index("barcode")
        xy = p8.loc[fine["barcode"], ["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
        _, nearest = cKDTree(p16[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()).query(xy)
        assert p16["barcode"].to_numpy()[nearest].tolist() == parents
    readme = zf.read("README.txt").decode("utf-8")
    n8 = len(matrix_barcodes(tiny_demo["samples"][0], "008um"))
    assert f"CTRL-1: the centre of each of the {n8:,} 8 um bins lies inside its 16 um bin" in readme
    fine = read_csv(zf, "CTRL-1_008um_categories.csv")
    labelled = fine["Tissue class"] != "not analyzed"
    assert (f'column "8 um QC": {(fine[QC] == "pass").sum():,} bins pass the run\'s quality control at 8 um; '
            f'{(labelled & (fine[QC] == "fail")).sum():,} of the {labelled.sum():,} with calls fail it') in readme
    assert 'The last column, "8 um QC", says "pass" or "fail"' in readme and "In every column but 8 um QC:" in readme


def test_readme(loupe_run, tiny_demo):
    readme = loupe_run["zip"].read("README.txt").decode("utf-8")
    assert readme.isascii() and "\r" not in readme and loupe_run["rid"] in readme and "for Loupe" in readme
    for s in tiny_demo["samples"]:
        assert f"Section {s['id']} ({s['label']})" in readme
        for bin_size in ("016um", "008um"):
            assert f"{s['id']}_{bin_size}_categories.csv" in readme
    for text in ("outs/cloupe_008um.cloupe", "--custom-bin-size 16", "+ Create a new group", "Import Categories",
                 "Empty barcode list supplied", '"not analyzed"', '"NA"', "Positive osteoblast pooled depthmatched",
                 "positivity call after spatial pooling and depth matching", "UTF-8 without a byte-order mark"):
        assert text in readme, text
    left_out = readme.split("Not exported:")[1].split("The columns follow")[0]
    for col in ("sample_id", "section", "x_um", "y_um", "total_umi",
                "positive_fraction_osteoblast_pooled_depthmatched_draws", "free_text", "empty_col"):
        assert f"  {col}: " in left_out, col
    assert "different values, too many for a category" in left_out.split("free_text: ")[1].split("\n")[0]
    assert "empty_col: no values in this run" in left_out
    # the size below which a Leiden domain is "unassigned" is the run's own (clustering.min_domain_bins)
    min_bins = pc.merge_with_defaults(demo_config(tiny_demo))["clustering"]["min_domain_bins"]
    assert f'"unassigned": domains under {min_bins:,} bins' in readme


def test_temporary_zip_is_removed(loupe_run, tmp_path, monkeypatch):
    c, rid = loupe_run["client"], loupe_run["rid"]
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    assert c.get(f"/api/runs/{rid}/export/loupe").status_code == 200
    assert list(tmp_path.iterdir()) == []


# ============================================================================
# THE EXPORT FUNCTION
# ============================================================================

def test_a_run_at_8um_has_one_file_per_section(tiny_demo):
    samples = tiny_demo["samples"]
    cfg = server.anchored(pc.merge_with_defaults({**demo_config(tiny_demo), "bin_size": "008um"}))
    raw = matrix_barcodes(samples[1], "008um")
    table = pd.DataFrame({"barcode": [f"{b[:-2]}-{samples[1]['id']}" for b in raw[::2]], "sample_id": samples[1]["id"],
                          "tissue_class": "bone lineage", "positive_osteoclast": True})
    buf = io.BytesIO()
    report = server.loupe_export(table, cfg, "008um", buf)
    zf = zipfile.ZipFile(buf)
    assert zf.namelist() == ["TUM-1_008um_categories.csv", "README.txt"]
    got = read_csv(zf, "TUM-1_008um_categories.csv")
    assert got.columns.tolist() == ["barcode", "Tissue class", "Positive osteoclast"] and got["barcode"].tolist() == raw
    assert got["Tissue class"].tolist() == ["bone lineage" if i % 2 == 0 else "not analyzed" for i in range(len(raw))]
    readme = zf.read("README.txt").decode("utf-8")
    assert "own Loupe file at 8 um" in readme and "8 UM FILES" not in readme and "custom-bin-size" not in readme
    assert report["sections"][0]["files"][0]["analyzed"] == len(raw[::2])
    assert "8 um QC" not in readme  # the run's own bins: "not analyzed" already marks those that fail


def test_grids_that_do_not_nest_use_the_bin_centres(tiny_demo, tmp_path):
    """When an 8 um centre lies outside the 16 um bin of its row // 2 and column // 2, every 8 um bin takes the calls
    of the 16 um bin whose centre is nearest."""
    sample = tiny_demo["samples"][0]
    binned = tmp_path / "M" / "outs" / "binned_outputs"
    shutil.copytree(sample["path"], binned)
    positions = binned / "square_008um" / "spatial" / "tissue_positions.parquet"
    mpp = json.loads((binned / "square_008um" / "spatial" / "scalefactors_json.json").read_text())["microns_per_pixel"]
    p8 = pd.read_parquet(positions)
    p8["pxl_col_in_fullres"] += 7.0 / mpp  # 7 um: about half of the 8 um centres move into another 16 um bin
    p8.to_parquet(positions, index=False)
    cfg = server.anchored(pc.merge_with_defaults({"data_root": str(tmp_path), "reference_sample": "M",
                                                  "samples": [{"id": "M", "path": str(binned)}]}))
    table = per_bin_table([{"id": "M", "path": str(binned)}])
    buf = io.BytesIO()
    report = server.loupe_export(table, cfg, "016um", buf)
    check = report["sections"][0]["files"][1]["parents"]
    assert check["method"] == "position" and 0 < check["inside"] < check["checked"] == len(p8[p8["barcode"].isin(
        matrix_barcodes({"path": str(binned)}, "008um"))])
    zf = zipfile.ZipFile(buf)
    coarse = read_csv(zf, "M_016um_categories.csv").set_index("barcode")
    fine = read_csv(zf, "M_008um_categories.csv")
    p16 = pd.read_parquet(binned / "square_016um" / "spatial" / "tissue_positions.parquet")
    xy = p8.set_index("barcode").loc[fine["barcode"], ["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    _, nearest = cKDTree(p16[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()).query(xy)
    expected = coarse.reindex(p16["barcode"].to_numpy()[nearest]).fillna("not analyzed")
    assert fine.drop(columns=["barcode", QC]).to_numpy().tolist() == expected.to_numpy().tolist()
    assert "took the calls of the 16 um bin whose centre is nearest" in zf.read("README.txt").decode("utf-8")


def test_columns_follow_the_table():
    """Whatever the signatures and steps of a run: every yes/no and text column is a category, numbers are not, a
    name taken twice gets a number, and the Leiden domains are words."""
    df = pd.DataFrame({"barcode": ["a", "b", "c"], "sample_id": "S", "leiden_domain": [0, 1, np.nan],
                       "positive_T_cell": [True, False, True], "Positive T cell": ["x", "y", "x"],
                       "in_tumour": [True, np.nan, False], "score": [0.1, 0.2, 0.3], "label": ["p,q", " ", "r\nz"]})
    columns, left_out = server.loupe_columns(df)
    assert [c.header for c in columns] == ["Leiden domain", "Positive T cell", "Positive T cell 2", "In tumour",
                                           "Label"]
    labels = {c.header: [c.categories[k] if k >= 0 else None for k in c.codes] for c in columns}
    assert labels == {"Leiden domain": ["domain 0", "domain 1", None], "Positive T cell": ["positive", "negative",
                                                                                          "positive"],
                      "Positive T cell 2": ["x", "y", "x"], "In tumour": ["inside", None, "outside"],
                      "Label": ["p;q", None, "r z"]}
    assert left_out == [("sample_id", "the section (one file per section)"), ("score", "a number, not a category")]


def test_export_refusals(tiny_demo, tmp_path):
    cfg = server.anchored(pc.merge_with_defaults(demo_config(tiny_demo)))
    moved = {**cfg, "samples": [dict(s, path=str(tmp_path / s["id"])) for s in cfg["samples"]]}
    table = per_bin_table(tiny_demo["samples"])
    for bad_cfg, df, bin_size, message in (
            (cfg, table, "002um", "8 um and larger"),
            (cfg, table.drop(columns="barcode"), "016um", "no barcode column"),
            (cfg, table.iloc[:0], "016um", "no bins"),
            (cfg, table.assign(sample_id="OTHER"), "016um", "not in the run's configuration"),
            (cfg, table.assign(barcode="s_016um_99999_99999-CTRL-1"), "016um", "is this the data the run analysed"),
            (moved, table, "016um", "can no longer be found")):
        with pytest.raises(server.DataError, match=message):
            server.loupe_export(df, bad_cfg, bin_size, io.BytesIO())


@pytest.mark.slow
def test_export_of_a_real_analysis(tiny_joint):
    """The per-bin table of a whole analysis (tiny_joint): every bin is placed, and the 8 um files hold the labels
    the analysis itself projects onto 8 um bins (bone_analysis.map_domains_from_coarser, step project_fine_bins)."""
    import anndata as ad
    import bone_analysis as ba
    cfg = server.anchored(tiny_joint["cfg"])
    table = pd.read_csv(tiny_joint["out"] / "per_bin_signature_calls_16um.csv.gz", low_memory=False)
    buf = io.BytesIO()
    report = server.loupe_export(table, cfg, "016um", buf)
    zf = zipfile.ZipFile(buf)
    for section in report["sections"]:
        first = section["files"][0]
        assert first["analyzed"] == (table["sample_id"] == section["id"]).sum() and first["unplaced"] == 0
    fine, _, _ = ba.load_sections(cfg, "008um")
    f = ad.concat(fine, label="sample_key", index_unique="-", join="outer")
    f.obs["sample"] = f.obs["sample"].astype(str)
    f = ba.map_domains_from_coarser(f, tiny_joint["a"], 8, 16)
    for sid in pc.sample_ids(cfg):
        got = read_csv(zf, f"{sid}_008um_categories.csv").set_index("barcode")
        obs = f.obs[f.obs["sample"] == sid]
        rows = got.loc[[b[: -len(sid) - 1] + "-1" for b in obs.index]]
        leiden = obs["leiden"].astype(str)
        expected = np.where(leiden == "no_parent_bin", "not analyzed",
                            np.where(leiden.str.isdigit(), "domain " + leiden, leiden))
        assert rows["Leiden domain"].tolist() == expected.tolist()
        program = obs["program"].astype(str)
        expected = [p if p == "ambiguous" else "not analyzed" if p == "no_parent_bin" else
                    pc.signature_label(cfg, p).replace(" ", "_") for p in program]
        assert rows["Assigned signature"].tolist() == expected
        # the analysis keeps exactly the 8 um bins whose "8 um QC" is "pass"
        assert (rows[QC] == "pass").all() and (got[QC] == "pass").sum() == len(obs)


def test_endpoint_errors(client, tiny_demo, tmp_path, monkeypatch):
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path / "tmp"))
    (tmp_path / "tmp").mkdir()
    r = client.get("/api/runs/r-20000101-000000-0000/export/loupe")
    assert r.status_code == 404 and "not found" in r.json()["detail"]

    failed = wait_for(client, start_run(client, new_project(client, fake_fail="boom")["id"])["id"], {"failed"})
    r = client.get(f"/api/runs/{failed['id']}/export/loupe")
    assert r.status_code == 404 and "no per-bin table" in r.json()["detail"] and "(failed)" in r.json()["detail"]

    # the fake worker's per-bin table names bins that are not in the demo data
    done = wait_for(client, start_run(client, new_project(client, **demo_config(tiny_demo))["id"])["id"],
                    {"succeeded"})
    r = client.get(f"/api/runs/{done['id']}/export/loupe")
    assert r.status_code == 422 and "is this the data the run analysed" in r.json()["detail"]

    binned = tmp_path / "moving" / "outs" / "binned_outputs"
    shutil.copytree(tiny_demo["samples"][0]["path"], binned)
    p = new_project(client, samples=[{"id": "M", "label": "moving", "path": str(binned)}], reference_sample="M")
    rid = wait_for(client, start_run(client, p["id"])["id"], {"succeeded"})["id"]
    write_per_bin(client.app.state.storage.results_dir(rid), per_bin_table([{"id": "M", "path": str(binned)}]))
    assert client.get(f"/api/runs/{rid}/export/loupe").status_code == 200
    shutil.rmtree(tmp_path / "moving")
    r = client.get(f"/api/runs/{rid}/export/loupe")
    assert r.status_code == 422 and "can no longer be found" in r.json()["detail"]
    assert list((tmp_path / "tmp").iterdir()) == []  # no temporary zip left behind, also after an error

    running = start_run(client, new_project(client, fake_sleep=60)["id"])
    wait_for(client, running["id"], {"running"})
    r = client.get(f"/api/runs/{running['id']}/export/loupe")
    assert r.status_code == 409 and "when it has finished" in r.json()["detail"]
    client.post(f"/api/runs/{running['id']}/cancel")
