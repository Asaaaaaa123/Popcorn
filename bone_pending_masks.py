#!/usr/bin/env python3
"""
bone_pending_masks.py
=====================
Hand-drawn region of interest (Methods): polygons drawn on a section's H&E image, stored as GeoJSON, turned into a
per-bin mask that bone_analysis.define_marrow_roi intersects with the transcriptional region of interest.

[APP] Added 2026-10-28. bone_analysis.py imported rasterise_drawn_roi from this module, which did not exist, so the
--roi-mask option always stopped the run with an ImportError ([FIX]). The polygons can be drawn in the web application
(project page, "Region of interest" tab) or in any tool that writes GeoJSON in image pixel coordinates (x = column,
y = row), such as QuPath or napari.

Frames: hires = spatial/tissue_hires_image.png (the web application's frame), lowres = tissue_lowres_image.png,
fullres = the full-resolution image that the Space Ranger pixel coordinates refer to.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath

FRAME_SCALEFACTOR = {"hires": "tissue_hires_scalef", "lowres": "tissue_lowres_scalef", "fullres": None}


def _polygons_of(geometry: dict) -> list:
    kind = (geometry or {}).get("type")
    if kind == "Polygon":
        return [geometry.get("coordinates") or []]
    if kind == "MultiPolygon":
        return list(geometry.get("coordinates") or [])
    if kind == "GeometryCollection":
        return [p for g in geometry.get("geometries", []) for p in _polygons_of(g)]
    return []


def load_polygons(source) -> list[tuple[np.ndarray, list[np.ndarray]]]:
    """
    GeoJSON (a path, JSON text, or the parsed object: FeatureCollection, Feature, Polygon, MultiPolygon or
    GeometryCollection) -> [(outer ring as an N x 2 array, [hole rings])] in (x, y) image coordinates.
    Rings with fewer than three points are ignored.
    """
    if isinstance(source, dict):
        obj = source
    elif isinstance(source, str) and source.lstrip().startswith("{"):
        obj = json.loads(source)
    else:
        obj = json.loads(Path(source).read_text(encoding="utf-8"))
    kind = obj.get("type")
    if kind == "FeatureCollection":
        geometries = [f.get("geometry") or {} for f in obj.get("features", [])]
    elif kind == "Feature":
        geometries = [obj.get("geometry") or {}]
    else:
        geometries = [obj]
    polygons = []
    for geometry in geometries:
        for polygon in _polygons_of(geometry):
            rings = [np.asarray(r, dtype=float)[:, :2] for r in polygon if len(r) >= 3]
            if rings:
                polygons.append((rings[0], rings[1:]))
    return polygons


def points_in_polygons(xy: np.ndarray, polygons) -> np.ndarray:
    """True for points inside an outer ring and outside that polygon's holes."""
    xy = np.asarray(xy, dtype=float)
    inside = np.zeros(len(xy), dtype=bool)
    for outer, holes in polygons:
        m = MplPath(outer).contains_points(xy)
        for hole in holes:
            m &= ~MplPath(hole).contains_points(xy)
        inside |= m
    return inside


def rasterise_drawn_roi(a, masks: dict, frame: str = "hires", scalefactors: dict | None = None):
    """
    Per-bin mask of the hand-drawn region of interest.
    a.obsm["spatial"] holds bin centres in full-resolution pixels (pxl_col, pxl_row = x, y); polygon coordinates are
    divided by the section's tissue_<frame>_scalef to reach that frame. masks maps sample id -> GeoJSON (path, text or
    object). scalefactors defaults to a.uns["spatial"][sample]["scalefactors"].
    Returns (mask, sections_without_region): a section without a region, or whose region holds no polygon, is not
    restricted (all True); the caller reports it.
    """
    if frame not in FRAME_SCALEFACTOR:
        raise ValueError(f"frame must be one of {sorted(FRAME_SCALEFACTOR)}, got {frame!r}")
    smp = a.obs["sample"].astype(str).to_numpy()
    xy = np.asarray(a.obsm["spatial"], dtype=float)
    keep = np.ones(a.n_obs, dtype=bool)
    unrestricted = []
    key = FRAME_SCALEFACTOR[frame]
    for s in pd.unique(smp):
        source = (masks or {}).get(s)
        polygons = load_polygons(source) if source is not None else []
        if not polygons:
            unrestricted.append(s)
            continue
        if key is None:
            factor = 1.0
        else:
            sf = (scalefactors or {}).get(s) or a.uns.get("spatial", {}).get(s, {}).get("scalefactors", {})
            factor = float(sf.get(key) or 0)
            if factor <= 0:
                raise ValueError(f"{s}: scalefactors_json.json has no {key}, so a region drawn on the {frame} image "
                                 f"cannot be placed on the bins")
        polygons = [(outer / factor, [h / factor for h in holes]) for outer, holes in polygons]
        m = smp == s
        keep[m] = points_in_polygons(xy[m], polygons)
    return keep, unrestricted
