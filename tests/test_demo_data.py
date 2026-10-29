"""The synthetic demo data looks like Space Ranger 3 Visium HD binned outputs and loads with every reader."""
from __future__ import annotations

import json

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image

import bone_analysis as ba
import pipeline_config as pc
from tools.make_demo_data import generate


def _bin_dir(sample, size):
    return pc.resolve_bin_dir({"data_root": ".", "samples": [sample]}, sample, size)


def test_folders_and_files(tiny_demo):
    for s in tiny_demo["samples"]:
        for size in ("008um", "016um"):
            d = _bin_dir(s, size)
            assert d is not None, (s["id"], size)
            for f in ("filtered_feature_bc_matrix.h5", "spatial/tissue_positions.parquet",
                      "spatial/scalefactors_json.json", "spatial/tissue_hires_image.png",
                      "spatial/tissue_lowres_image.png"):
                assert (d / f).is_file(), f


def test_reads_with_scanpy_and_both_readers(tiny_demo, original_ba):
    for s in tiny_demo["samples"]:
        d = _bin_dir(s, "016um")
        raw = sc.read_10x_h5(d / "filtered_feature_bc_matrix.h5")
        new = ba.read_visium_hd(d, s["id"], s["condition"])
        original_ba.SAMPLES.setdefault(s["id"], {"condition": s["condition"], "limb": "-", "label": s["label"]})
        old = original_ba.read_visium_hd(d, s["id"])
        assert new.n_obs == raw.n_obs == old.n_obs                               # every barcode has a position
        assert list(new.obs_names) == list(old.obs_names)
        np.testing.assert_array_equal(new.obsm["spatial"], old.obsm["spatial"])
        assert "Fgf23" not in new.var_names and "Bglap" not in new.var_names    # absent from the probe set
        assert new.var_names.str.lower().str.startswith("mt-").sum() >= 10


def test_16um_bins_are_sums_of_their_8um_children(tiny_demo):
    s = tiny_demo["samples"][0]
    fine = ba.read_visium_hd(_bin_dir(s, "008um"), s["id"])
    coarse = ba.read_visium_hd(_bin_dir(s, "016um"), s["id"])
    parent = list(zip(fine.obs["array_row"] // 2, fine.obs["array_col"] // 2))
    sums = pd.DataFrame(fine.X.toarray(), index=pd.MultiIndex.from_tuples(parent)).groupby(level=[0, 1]).sum()
    key = pd.MultiIndex.from_arrays([coarse.obs["array_row"], coarse.obs["array_col"]])
    np.testing.assert_array_equal(sums.reindex(key).fillna(0).to_numpy(), coarse.X.toarray())


def test_scalefactors_and_images_are_consistent(tiny_demo):
    for s in tiny_demo["samples"]:
        d = _bin_dir(s, "016um")
        sf = json.loads((d / "spatial" / "scalefactors_json.json").read_text())
        assert 0.4 < sf["microns_per_pixel"] < 0.8 and sf["bin_size_um"] == 16
        tp = pd.read_parquet(d / "spatial" / "tissue_positions.parquet")
        step = np.diff(np.sort(tp.loc[tp["array_row"] == tp["array_row"].min(), "pxl_col_in_fullres"].to_numpy()))
        assert np.allclose(step * sf["microns_per_pixel"], 16, atol=0.5)        # adjacent bins 16 um apart
        hires = Image.open(d / "spatial" / "tissue_hires_image.png")
        assert max(hires.size) <= 2000
        x_max = tp["pxl_col_in_fullres"].max() * sf["tissue_hires_scalef"]
        y_max = tp["pxl_row_in_fullres"].max() * sf["tissue_hires_scalef"]
        assert x_max <= hires.size[0] and y_max <= hires.size[1]               # bins fall on the image
        with h5py.File(d / "filtered_feature_bc_matrix.h5") as f:
            assert list(f["matrix/shape"][()])[1] == int(tp["in_tissue"].sum())


def test_marker_and_determinism(tmp_path, tiny_demo):
    first = generate(tmp_path, size="tiny", seed=3)
    h5 = _bin_dir(first[0], "016um") / "filtered_feature_bc_matrix.h5"
    stamp = h5.stat().st_mtime_ns
    assert generate(tmp_path, size="tiny", seed=3) == first                      # .complete short-circuits
    assert h5.stat().st_mtime_ns == stamp
    with h5py.File(h5) as f:
        data3 = f["matrix/data"][()]
    other = generate(tmp_path / "again", size="tiny", seed=3)
    with h5py.File(_bin_dir(other[0], "016um") / "filtered_feature_bc_matrix.h5") as f:
        np.testing.assert_array_equal(f["matrix/data"][()], data3)
