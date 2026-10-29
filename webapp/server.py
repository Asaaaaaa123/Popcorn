"""
webapp/server.py
================
FastAPI application of Bone Spatial Studio (docs/ARCHITECTURE.md, section 5): the JSON API under /api/, the front end
under /static/ and its index.html at /.

The server process stays light and responsive: it never imports scanpy, anndata or bone_analysis (an analysis runs in
a worker subprocess started by webapp.jobs). Package versions come from the installed package metadata, the data
checks read the Space Ranger files directly (h5py for the feature-barcode matrix, pandas/pyarrow for the tissue
positions, Pillow for the image size), and every JSON response is strict JSON: numpy values become plain numbers,
NaN and infinities become null.
"""
from __future__ import annotations

import collections
import contextlib
import functools
import importlib
import importlib.metadata
import json
import logging
import math
import mimetypes
import os
import platform
import re
import secrets
import stat
import string
import threading
from pathlib import Path
from urllib.parse import quote

import h5py
import numpy as np
import pandas as pd
from fastapi import Depends, FastAPI, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from PIL import Image
from starlette.exceptions import HTTPException as StarletteHTTPException

import pipeline_config as pc
from webapp import __version__
from webapp.jobs import JobManager
from webapp.storage import (ACTIVE_STATUSES, FINAL_STATUSES, BadRequest, NotFound, Storage, StorageError,
                            check_sample_id, json_safe)

log = logging.getLogger("bss.server")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
STATIC_DIR = Path(__file__).resolve().parent / "static"

# Import name -> installed version, in the order of GET /api/health. Required packages must be present for a run to
# start; the analysis has fallbacks for the optional ones.
PACKAGES = ("scanpy", "anndata", "numpy", "pandas", "scipy", "sklearn", "statsmodels", "harmonypy", "igraph",
            "leidenalg", "squidpy", "pyarrow", "h5py", "matplotlib", "PIL")
REQUIRED_PACKAGES = ("scanpy", "anndata", "numpy", "pandas", "scipy", "sklearn", "statsmodels", "igraph", "leidenalg",
                     "pyarrow", "h5py", "matplotlib", "PIL")
OPTIONAL_PACKAGES = ("harmonypy", "squidpy")
_DISTRIBUTIONS = {"sklearn": ("scikit-learn",), "PIL": ("pillow",), "igraph": ("igraph", "python-igraph")}

DOCS = {"readme": "README.md", "changes": "CHANGES.md", "architecture": "docs/ARCHITECTURE.md"}

DEMO_NAME = "Demo — synthetic tibiae"
# Fewer permutations and replicates than the manuscript, so that the demo runs in a few minutes.
DEMO_SETTINGS = {
    "scoring": {"n_perm": 20, "n_perm_resampled": 10, "n_depth_reps": 2, "null_n_sub": 5000},
    "clustering": {"stability_seeds": [0, 1]},
    "moran": {"n_genes": 500},
    "uncoupling": {"n_perm": 99},
}

TOKEN_COOKIE = "bss_token"
TOKEN_HEADER = "x-token"
TOKEN_RE = re.compile(r"^[A-Za-z0-9._~-]{1,256}\Z")
TOKEN_COOKIE_AGE = 30 * 24 * 3600

MAX_BODY = 32 * 1024 * 1024        # bytes of a JSON request body
LOG_CHUNK = 256 * 1024             # bytes of log.txt per GET /log
PREVIEW_MAX_BINS = 60_000          # bin centres per GET /bins-preview
TABLE_MAX_LIMIT = 1000             # rows per GET /table
MAX_ROI_VERTICES = 500_000
MAX_ERRORS_SHOWN = 8               # validation errors quoted in a 409
SCAN_MAX_DEPTH = 4                 # detect-samples: folder levels below data_root
SCAN_MAX_DIRS = 20_000             # detect-samples: folders visited at most
LISTING_PEEK = 500                 # fs/list: entries of a child folder inspected for Space Ranger outputs

POSITION_FILES = ("tissue_positions.parquet", "tissue_positions.csv", "tissue_positions_list.csv")
POSITION_COLUMNS = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
SQUARE_RE = re.compile(r"^square_(\d{3}um)$")
RUN_PREFIX_RE = re.compile(r"^SR\d+_")
TABLE_SUFFIXES = (".csv", ".csv.gz", ".tsv", ".tsv.gz")
FRAMES = ("hires", "lowres", "fullres")

# per-bin table (docs/ARCHITECTURE.md, section 3)
BIN_SKIP_COLUMNS = {"barcode", "sample_id", "section", "x_um", "y_um"}
BIN_CATEGORY_COLUMNS = ("tissue_class", "assigned_signature", "leiden_domain")
BIN_MAX_CATEGORIES = 500
REMOVED_COLOR = "#b8b8b8"
UNASSIGNED_COLOR = "#dedede"

# The Windows registry can map .js to text/plain (browsers then refuse to load the ES modules of the front end) and
# .csv to Excel; serve the standard types whatever the registry says.
for _media_type, _ext in (("text/javascript", ".js"), ("text/javascript", ".mjs"), ("text/css", ".css"),
                          ("text/html", ".html"), ("application/json", ".json"), ("application/geo+json", ".geojson"),
                          ("image/svg+xml", ".svg"), ("image/png", ".png"), ("application/pdf", ".pdf"),
                          ("text/csv", ".csv"), ("text/tab-separated-values", ".tsv"), ("text/markdown", ".md")):
    mimetypes.add_type(_media_type, _ext)
_COMPRESSED_TYPES = {"gzip": "application/gzip", "bzip2": "application/x-bzip2", "xz": "application/x-xz"}


class DataError(Exception):
    """A data file exists but cannot be read (HTTP 422)."""


JSON_MEDIA_TYPE = "application/json; charset=utf-8"  # the charset helps clients that do not assume UTF-8 for JSON


class JSONOut(JSONResponse):
    """JSON response whose content is made strictly JSON-safe first (numpy values, NaN and infinities -> null)."""

    media_type = JSON_MEDIA_TYPE

    def render(self, content) -> bytes:
        return json.dumps(json_safe(content), ensure_ascii=False, allow_nan=False,
                          separators=(",", ":")).encode("utf-8")


class NoCacheStaticFiles(StaticFiles):
    """Static files; the browser revalidates the front-end code on every load, so an update is never half-applied."""

    def file_response(self, full_path, stat_result, scope, status_code: int = 200) -> Response:
        response = super().file_response(full_path, stat_result, scope, status_code)
        if str(full_path).endswith((".html", ".js", ".mjs", ".css")):
            response.headers["Cache-Control"] = "no-cache"
        return response


# ============================================================================
# ENVIRONMENT
# ============================================================================

@functools.lru_cache(maxsize=1)
def package_versions() -> dict[str, str | None]:
    """Installed version of every package of the analysis (import name -> version or None), read from the package
    metadata: nothing is imported, so the server never loads scanpy."""
    out: dict[str, str | None] = {}
    for name in PACKAGES:
        out[name] = None
        for dist in _DISTRIBUTIONS.get(name, (name,)):
            try:
                out[name] = importlib.metadata.version(dist)
                break
            except importlib.metadata.PackageNotFoundError:
                continue
    return out


def missing_required() -> list[str]:
    versions = package_versions()
    return [name for name in REQUIRED_PACKAGES if versions.get(name) is None]


@functools.lru_cache(maxsize=1)
def platform_text() -> str:
    return platform.platform()


# ============================================================================
# CONFIGURATION AND SAMPLE FILES
# ============================================================================

def check_config(cfg) -> tuple[list[str], list[str]]:
    """pipeline_config.validate_config, also for drafts too malformed for it to inspect."""
    if not isinstance(cfg, dict):
        return ["the configuration must be a JSON object"], []
    try:
        errors, warnings = pc.validate_config(cfg)
    except Exception as exc:  # validate_config expects the shapes of the defaults; a saved draft may hold anything
        return [f"the configuration is malformed ({type(exc).__name__}: {exc})"], []
    return list(errors), list(warnings)


def summarize(errors: list[str]) -> str:
    shown = "; ".join(errors[:MAX_ERRORS_SHOWN])
    return shown + (f"; and {len(errors) - MAX_ERRORS_SHOWN} more" if len(errors) > MAX_ERRORS_SHOWN else "")


def anchored(cfg: dict) -> dict:
    """`cfg` with an absolute data_root. A relative data_root is relative to the project root, the working directory
    of the worker, so the server finds the same files as the analysis whatever its own working directory."""
    root = Path(str(cfg.get("data_root") or ".")).expanduser()
    if not root.is_absolute():
        root = PROJECT_ROOT / root
    return {**cfg, "data_root": os.path.abspath(root)}


def config_samples(cfg: dict) -> list[dict]:
    """The sample entries that can locate data (a dict with a string id)."""
    samples = cfg.get("samples")
    if not isinstance(samples, list):
        return []
    return [s for s in samples if isinstance(s, dict) and isinstance(s.get("id"), str) and s["id"]]


def sample_base(cfg: dict, sample: dict) -> Path:
    """The folder a sample entry points at (pipeline_config resolves square_XXXum/ below it)."""
    p = Path(str(sample.get("path") or sample["id"])).expanduser()
    return p if p.is_absolute() else Path(cfg["data_root"]) / p


def spatial_dirs(cfg: dict, sample: dict) -> list[Path]:
    """spatial/ folders of a sample: the analysis bin size first, then the QC bin size, then the others."""
    sizes = []
    for b in (cfg.get("bin_size"), cfg.get("qc_bin_size"), *pc.available_bin_sizes(cfg, sample)):
        if isinstance(b, str) and b not in sizes:
            sizes.append(b)
    out = []
    for b in sizes:
        d = pc.resolve_bin_dir(cfg, sample, b)
        if d is not None and (d / "spatial").is_dir():
            out.append(d / "spatial")
    return out


def hires_image(cfg: dict, sample: dict) -> Path | None:
    """tissue_hires_image.png of the analysis bin size, else of any other bin size of the sample."""
    for d in spatial_dirs(cfg, sample):
        if (d / "tissue_hires_image.png").is_file():
            return d / "tissue_hires_image.png"
    return None


def positions_file(spatial: Path) -> Path | None:
    return next((spatial / n for n in POSITION_FILES if (spatial / n).is_file()), None)


def read_scalefactors(spatial: Path) -> dict | None:
    path = spatial / "scalefactors_json.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None
    except (OSError, ValueError) as exc:
        raise DataError(f"{path} cannot be read ({exc})") from exc
    if not isinstance(data, dict):
        raise DataError(f"{path} is not a JSON object")
    return data


def positive_number(value) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        return None
    return float(value)


def microns_per_pixel(sample: dict, scalefactors: dict | None) -> float | None:
    """The configured microns per full-resolution pixel, else the one Space Ranger recorded."""
    return positive_number(sample.get("microns_per_pixel")) or positive_number((scalefactors or {}).get(
        "microns_per_pixel"))


@functools.lru_cache(maxsize=32)
def _matrix_info(path: str, mtime_ns: int, size: int) -> tuple[int, tuple[str, ...]]:
    with h5py.File(path, "r") as f:
        matrix = f["matrix"]
        shape = matrix["shape"][()]
        names = matrix["features"]["name"][()] if "features" in matrix else matrix["gene_names"][()]
    return int(shape[1]), tuple(n.decode("utf-8", "replace") if isinstance(n, bytes) else str(n) for n in names)


def matrix_info(bin_dir: Path) -> tuple[int, tuple[str, ...]]:
    """(number of bins, gene names) of square_XXXum/filtered_feature_bc_matrix.h5, read from the HDF5 header."""
    path = bin_dir / "filtered_feature_bc_matrix.h5"
    try:
        st = path.stat()
        return _matrix_info(str(path), st.st_mtime_ns, st.st_size)
    except (OSError, KeyError, ValueError, IndexError, TypeError) as exc:
        raise DataError(f"{path} is not a readable feature-barcode matrix ({exc})") from exc


@functools.lru_cache(maxsize=6)
def _in_tissue_xy(path: str, mtime_ns: int, size: int) -> tuple[np.ndarray, np.ndarray]:
    cols = ["in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    if path.endswith(".parquet"):
        df = pd.read_parquet(path, columns=cols)
    elif Path(path).name == "tissue_positions_list.csv":  # Space Ranger 1.x: no header
        df = pd.read_csv(path, header=None, names=POSITION_COLUMNS, usecols=cols)
    else:
        df = pd.read_csv(path, usecols=cols)
    keep = pd.to_numeric(df["in_tissue"], errors="coerce").fillna(0).to_numpy() == 1
    x = pd.to_numeric(df["pxl_col_in_fullres"], errors="coerce").to_numpy(float)[keep]
    y = pd.to_numeric(df["pxl_row_in_fullres"], errors="coerce").to_numpy(float)[keep]
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    x.setflags(write=False)
    y.setflags(write=False)
    return x, y


def in_tissue_xy(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Full-resolution pixel coordinates (x = pxl_col, y = pxl_row) of the in-tissue bins of a positions file."""
    try:
        st = path.stat()
        return _in_tissue_xy(str(path), st.st_mtime_ns, st.st_size)
    except (OSError, KeyError, ValueError, ImportError) as exc:
        raise DataError(f"{path} is not a readable tissue positions table ({exc})") from exc


def mask_frames(path: Path) -> set[str]:
    """Image frames the polygons of a stored region declare (properties.frame, hires when absent)."""
    data = json.loads(path.read_text(encoding="utf-8"))
    frames = set()
    for feature in data.get("features") or []:
        props = feature.get("properties") if isinstance(feature, dict) else None
        frames.add(str((props or {}).get("frame") or "hires"))
    return frames


# ============================================================================
# DATA CHECKS (POST /validate, POST /detect-samples, GET /fs/list)
# ============================================================================

def sample_report(cfg: dict, sample: dict) -> tuple[dict, list[str], list[str], tuple[str, ...] | None]:
    """Report of one sample of POST /validate: (entry, errors, warnings, gene names or None)."""
    sid = sample["id"]
    entry = {"id": sid, "bin_dir": None, "bin_sizes": [], "has_positions": False, "has_scalefactors": False,
             "has_hires_image": False, "microns_per_pixel": None, "n_bins": None, "message": ""}
    errors, warnings = [], []
    bin_size = cfg.get("bin_size")
    entry["bin_sizes"] = pc.available_bin_sizes(cfg, sample)
    bin_dir = pc.resolve_bin_dir(cfg, sample)
    if bin_dir is None:
        found = f" (bin sizes found: {', '.join(entry['bin_sizes'])})" if entry["bin_sizes"] else ""
        entry["message"] = f"no square_{bin_size}/filtered_feature_bc_matrix.h5 in {sample_base(cfg, sample)}{found}"
        errors.append(f"sample {sid}: {entry['message']}")
        return entry, errors, warnings, None
    entry["bin_dir"] = str(bin_dir)
    spatial = bin_dir / "spatial"
    problems, notes = [], []
    entry["has_positions"] = positions_file(spatial) is not None
    if not entry["has_positions"]:
        problems.append("no tissue positions (spatial/tissue_positions.parquet)")
        errors.append(f"sample {sid}: no tissue positions file in {spatial}")
    scalefactors = None
    try:
        scalefactors = read_scalefactors(spatial)
    except DataError as exc:
        problems.append(str(exc))
        errors.append(f"sample {sid}: {exc}")
    entry["has_scalefactors"] = scalefactors is not None
    entry["has_hires_image"] = hires_image(cfg, sample) is not None
    if not entry["has_hires_image"]:
        notes.append("no H&E image, so no region can be drawn")
        warnings.append(f"sample {sid}: no spatial/tissue_hires_image.png; a region of interest cannot be drawn")
    entry["microns_per_pixel"] = microns_per_pixel(sample, scalefactors)
    if entry["microns_per_pixel"] is None:
        problems.append("microns per pixel unknown")
        errors.append(f"sample {sid}: microns_per_pixel is not set and spatial/scalefactors_json.json does not "
                      f"record it")
    qc = cfg.get("qc_bin_size")
    steps = cfg.get("steps") if isinstance(cfg.get("steps"), dict) else {}
    if steps.get("qc_figure_bins") and qc != bin_size and qc not in entry["bin_sizes"]:
        notes.append(f"no square_{qc}")
        warnings.append(f"sample {sid}: no square_{qc}/ for the QC figure; the analysis bin size is used instead")
    genes = None
    try:
        entry["n_bins"], genes = matrix_info(bin_dir)
    except DataError as exc:
        problems.append("the feature-barcode matrix cannot be read")
        errors.append(f"sample {sid}: {exc}")
    if problems:
        entry["message"] = "; ".join(problems)
    else:
        mpp = entry["microns_per_pixel"]
        entry["message"] = "; ".join([f"{entry['n_bins']:,} bins at {bin_size}, {mpp:.4g} µm per pixel"] + notes)
    return entry, errors, warnings, genes


def gene_list(value) -> list[str]:
    return [g for g in value if isinstance(g, str)] if isinstance(value, list) else []


def genes_report(cfg: dict, panels: list[tuple[str, ...]]) -> tuple[dict, list[str]]:
    """Gene part of POST /validate over the gene panels of the readable samples: (report, warnings)."""
    report = {"n_features": None, "signatures": {}, "receptor_panel": {"missing": []},
              "adipocyte_genes": {"missing": []}}
    if not panels:
        return report, []
    warnings = []
    present = set(panels[0]).intersection(*panels[1:])
    report["n_features"] = len(panels[0])
    if any(set(p) != set(panels[0]) for p in panels[1:]):
        warnings.append(f"the sections have different gene panels; {len(set().union(*panels) - present)} genes are "
                        f"not measured in every section")
    scoring = cfg.get("scoring") if isinstance(cfg.get("scoring"), dict) else {}
    min_core = scoring.get("min_core", 0)
    primary = cfg.get("primary") if isinstance(cfg.get("primary"), list) else []
    signatures = cfg.get("signatures") if isinstance(cfg.get("signatures"), dict) else {}
    for name, sig in signatures.items():
        if not isinstance(sig, dict):
            continue
        core = gene_list(sig.get("core"))
        scored = list(dict.fromkeys(gene_list(sig.get("extended")) + core))
        hit = [g for g in scored if g in present]
        miss = [g for g in scored if g not in present]
        n_core = sum(g in present for g in core)
        report["signatures"][name] = {"present": hit, "missing": miss, "core_present": n_core}
        if miss:
            warnings.append(f"signature {name}: {len(miss)} of {len(scored)} genes missing from the probe set "
                            f"({', '.join(miss)})")
        if len(hit) < 2:
            warnings.append(f"signature {name}: fewer than two of its genes are measured; it cannot be scored")
        elif (name in primary and isinstance(min_core, (int, float)) and not isinstance(min_core, bool)
              and n_core < min_core):
            warnings.append(f"signature {name}: only {n_core} core genes are measured, fewer than scoring.min_core "
                            f"({min_core}); no bin can be called positive")
    panel = cfg.get("receptor_panel") if isinstance(cfg.get("receptor_panel"), dict) else {}
    report["receptor_panel"]["missing"] = [g for g in panel if g not in present]
    report["adipocyte_genes"]["missing"] = [g for g in gene_list(cfg.get("adipocyte_genes")) if g not in present]
    for key, label in (("receptor_panel", "host receptor panel"), ("adipocyte_genes", "adipocyte genes")):
        miss = report[key]["missing"]
        if miss:
            warnings.append(f"{label}: {len(miss)} genes missing from the probe set ({', '.join(miss)})")
    return report, warnings


def suggested_sample_id(folder: str) -> str:
    """A valid sample id from a Space Ranger output folder name, without its leading SR<digits>_ run prefix."""
    text = re.sub(r"[^A-Za-z0-9_.\-]+", "-", RUN_PREFIX_RE.sub("", folder))
    text = re.sub(r"^[^A-Za-z0-9]+", "", text)[:64].rstrip(".-_")
    return text or "sample"


def sample_folder_name(path: Path) -> str:
    """Name of the sample folder of a binned-outputs folder (.../<sample>/outs/binned_outputs -> <sample>)."""
    parts = list(path.parts)
    while len(parts) > 1 and parts[-1].lower() in ("binned_outputs", "outs"):
        parts.pop()
    return parts[-1] if parts else str(path)


def hidden(entry: os.DirEntry) -> bool:
    if entry.name.startswith((".", "$")):
        return True
    try:
        attrs = getattr(entry.stat(follow_symlinks=False), "st_file_attributes", 0)
    except OSError:
        return True
    return bool(attrs & (getattr(stat, "FILE_ATTRIBUTE_HIDDEN", 0) | getattr(stat, "FILE_ATTRIBUTE_SYSTEM", 0)))


def subfolders(path: Path, follow_links: bool = True) -> list[os.DirEntry]:
    """Visible subfolders of `path`. follow_links=False leaves out symbolic links and junctions, so that a recursive
    scan cannot loop."""
    out = []
    with os.scandir(path) as it:
        for e in it:
            try:
                if not e.is_dir(follow_symlinks=follow_links) or hidden(e):
                    continue
                if not follow_links and (e.is_symlink() or getattr(e, "is_junction", lambda: False)()):
                    continue
            except OSError:
                continue
            out.append(e)
    return out


def detect_samples(root: Path) -> list[dict]:
    """Folders at most SCAN_MAX_DEPTH levels below `root` (root included) that contain
    square_XXXum/filtered_feature_bc_matrix.h5."""
    found = []
    queue = collections.deque([(root, 0)])
    visited = 0
    while queue and visited < SCAN_MAX_DIRS:
        folder, depth = queue.popleft()
        visited += 1
        try:
            children = subfolders(folder, follow_links=False)
        except OSError:
            continue
        bins = sorted(m.group(1) for e in children if (m := SQUARE_RE.match(e.name))
                      and os.path.isfile(os.path.join(e.path, "filtered_feature_bc_matrix.h5")))
        if bins:
            found.append((folder, bins))
            continue
        if depth < SCAN_MAX_DEPTH:
            queue.extend((Path(e.path), depth + 1) for e in sorted(children, key=lambda e: e.name.casefold())
                         if e.name != "spatial" and not SQUARE_RE.match(e.name))
    candidates, used = [], collections.Counter()
    for folder, bins in found:
        try:
            rel = Path(os.path.relpath(folder, root)).as_posix()
        except ValueError:  # another drive (Windows)
            rel = str(folder)
        name = sample_folder_name(folder)
        sid = suggested_sample_id(name)
        used[sid] += 1
        if used[sid] > 1:
            sid = f"{sid[:60]}-{used[sid]}"
        candidates.append({"folder": name, "path": rel, "suggested_id": sid, "bin_sizes": bins})
    return sorted(candidates, key=lambda c: c["path"].casefold())


def looks_like_spaceranger(path: str) -> bool:
    """The folder contains outs/, binned_outputs/ or square_XXXum/."""
    try:
        with os.scandir(path) as it:
            for i, e in enumerate(it):
                if i >= LISTING_PEEK:
                    break
                if (e.name in ("outs", "binned_outputs") or SQUARE_RE.match(e.name)) and e.is_dir():
                    return True
    except OSError:
        return False
    return False


def filesystem_roots() -> list[str]:
    if os.name != "nt":
        return ["/"]
    try:
        return sorted(os.listdrives())
    except (AttributeError, OSError):  # Python < 3.12
        return [f"{c}:\\" for c in string.ascii_uppercase if os.path.exists(f"{c}:\\")]


# ============================================================================
# REGIONS OF INTEREST
# ============================================================================

def _ring(ring, where: str) -> list[list[float]]:
    if not isinstance(ring, list):
        raise HTTPException(400, f"{where}: a ring must be a list of [x, y] positions")
    points = []
    for p in ring:
        if (not isinstance(p, list) or len(p) < 2
                or any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) for v in p[:2])):
            raise HTTPException(400, f"{where}: every position must be [x, y] with finite numbers")
        points.append([float(p[0]), float(p[1])])
    if len(points) > 1 and points[0] == points[-1]:
        points.pop()
    if len({tuple(p) for p in points}) < 3:
        raise HTTPException(400, f"{where}: a ring needs at least three distinct points")
    return points + [list(points[0])]


def clean_feature_collection(body: dict) -> tuple[dict, int]:
    """A drawn region checked and normalised (rings closed, frame recorded): (FeatureCollection, number of
    polygons). HTTP 400 for anything that is not a FeatureCollection of Polygon or MultiPolygon features."""
    if body.get("type") != "FeatureCollection" or not isinstance(body.get("features"), list):
        raise HTTPException(400, "expected a GeoJSON FeatureCollection with a list of features")
    features, n_polygons, n_vertices = [], 0, 0
    for i, feature in enumerate(body["features"], 1):
        where = f"feature {i}"
        if not isinstance(feature, dict) or feature.get("type") != "Feature":
            raise HTTPException(400, f"{where} is not a GeoJSON Feature")
        geometry = feature.get("geometry")
        kind = geometry.get("type") if isinstance(geometry, dict) else None
        if kind not in ("Polygon", "MultiPolygon"):
            raise HTTPException(400, f"{where}: the geometry must be a Polygon or a MultiPolygon, got {kind!r}")
        coords = geometry.get("coordinates")
        polygons = [coords] if kind == "Polygon" else coords
        if not isinstance(polygons, list) or not polygons:
            raise HTTPException(400, f"{where}: no coordinates")
        cleaned = []
        for polygon in polygons:
            if not isinstance(polygon, list) or not polygon:
                raise HTTPException(400, f"{where}: a polygon must be a list of rings (outer ring first)")
            cleaned.append([_ring(r, where) for r in polygon])
        n_polygons += len(cleaned)
        n_vertices += sum(len(r) for poly in cleaned for r in poly)
        if n_vertices > MAX_ROI_VERTICES:
            raise HTTPException(400, f"the region has more than {MAX_ROI_VERTICES:,} vertices; simplify it")
        props = feature.get("properties") if isinstance(feature.get("properties"), dict) else {}
        props = {**props, "frame": props.get("frame") or "hires"}
        if props["frame"] not in FRAMES:
            raise HTTPException(400, f"{where}: properties.frame must be one of {', '.join(FRAMES)}")
        features.append({"type": "Feature", "properties": props,
                         "geometry": {"type": kind, "coordinates": cleaned[0] if kind == "Polygon" else cleaned}})
    return {"type": "FeatureCollection", "features": features}, n_polygons


# ============================================================================
# RESULT TABLES
# ============================================================================
_DIGITS = re.compile(r"([0-9]+)")


def natural_key(text: str) -> tuple:
    """Sort key that orders the digit runs of a text by value ("domain 2" before "domain 10")."""
    return tuple((0, int(part), "") if part[0] in "0123456789" else (1, 0, part.casefold())
                 for part in _DIGITS.split(text) if part)


def column_kind(s: pd.Series) -> str:
    if pd.api.types.is_bool_dtype(s):
        return "bool"
    if pd.api.types.is_numeric_dtype(s):
        return "number"
    present = s.dropna()
    if len(present) and all(isinstance(v, (bool, np.bool_)) for v in present.tolist()):
        return "bool"
    return "text"


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.name.lower().endswith((".tsv", ".tsv.gz")) else ","
    try:
        return pd.read_csv(path, sep=sep, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


class Table:
    """A result table loaded once, with its sort orders and search text computed on first use."""

    def __init__(self, df: pd.DataFrame):
        self.df = df
        self.columns = [str(c) for c in df.columns]
        self.df.columns = self.columns
        self.dtypes = {c: column_kind(df[c]) for c in self.columns}
        self._lock = threading.Lock()
        self._orders: dict[tuple[str, bool], np.ndarray] = {}
        self._search: pd.Series | None = None
        self._last_query: tuple[str, np.ndarray] | None = None

    def sort_values(self, col: str) -> np.ndarray | None:
        """Float values to sort a column by numerically (NaN = missing), or None to sort it as text."""
        s = self.df[col]
        kind = self.dtypes[col]
        if kind == "bool":
            return s.map({True: 1.0, False: 0.0}).to_numpy(float, na_value=np.nan)
        num = pd.to_numeric(s, errors="coerce")
        if kind == "number" or num.notna().sum() == s.notna().sum():
            return num.to_numpy(float, na_value=np.nan)
        return None

    def order(self, col: str, desc: bool) -> np.ndarray:
        """Row positions sorted by `col` (stable; missing values last in both directions)."""
        key = (col, desc)
        with self._lock:
            cached = self._orders.get(key)
        if cached is not None:
            return cached
        values = self.sort_values(col)
        if values is not None:
            missing = np.isnan(values)
            present = np.flatnonzero(~missing)
            ranked = present[np.argsort(-values[present] if desc else values[present], kind="stable")]
            result = np.concatenate([ranked, np.flatnonzero(missing)])
        else:
            s = self.df[col]
            texts = s.tolist()
            na = s.isna().to_numpy()
            present = [i for i in range(len(texts)) if not na[i]]
            present.sort(key=lambda i: natural_key(str(texts[i])), reverse=desc)
            result = np.array(present + np.flatnonzero(na).tolist(), dtype=np.int64)
        with self._lock:
            self._orders[key] = result
        return result

    def matches(self, query: str) -> np.ndarray:
        """Boolean mask of the rows where any column contains `query` (case-insensitive)."""
        q = query.casefold()
        with self._lock:
            if self._last_query is not None and self._last_query[0] == q:
                return self._last_query[1]
            search = self._search
        if search is None:
            joined = None
            for col in self.columns:
                s = self.df[col]
                text = s.astype(str).str.casefold().where(s.notna(), "")
                joined = text if joined is None else joined + "\x1f" + text
            search = joined if joined is not None else pd.Series([], dtype=object)
        mask = search.str.contains(q, regex=False).to_numpy(bool, na_value=False)
        with self._lock:
            self._search = search
            self._last_query = (q, mask)
        return mask

    def rows(self, positions: np.ndarray) -> list[list]:
        page = self.df.iloc[positions]
        values = page.to_numpy(dtype=object)
        values[page.isna().to_numpy()] = None
        return json_safe(values.tolist())


class LRUCache:
    """Small thread-safe LRU mapping; values are built outside the lock."""

    def __init__(self, size: int):
        self.size = size
        self._items: collections.OrderedDict = collections.OrderedDict()
        self._lock = threading.Lock()

    def get(self, key, build):
        with self._lock:
            if key in self._items:
                self._items.move_to_end(key)
                return self._items[key]
        value = build()
        with self._lock:
            self._items[key] = value
            self._items.move_to_end(key)
            while len(self._items) > self.size:
                self._items.popitem(last=False)
        return value

    def discard(self, predicate) -> None:
        with self._lock:
            for key in [k for k in self._items if predicate(k)]:
                del self._items[key]


def file_key(path: Path) -> tuple[str, int, int]:
    st = path.stat()
    return str(path), st.st_mtime_ns, st.st_size


# ============================================================================
# PER-BIN MAPS (GET /bins)
# ============================================================================

def signature_slug(cfg: dict, name: str) -> str:
    return pc.signature_label(cfg, name).replace(" ", "_")


def bin_colors(cfg: dict) -> dict:
    classes = {c["name"]: c["color"] for c in cfg.get("classes") or []
               if isinstance(c, dict) and c.get("name") and c.get("color")}
    signatures = {signature_slug(cfg, n): d["color"] for n, d in (cfg.get("signatures") or {}).items()
                  if isinstance(d, dict) and d.get("color")}
    return {"tissue_class": {**classes, "removed": REMOVED_COLOR, "unassigned": UNASSIGNED_COLOR},
            "assigned_signature": {**signatures, "ambiguous": UNASSIGNED_COLOR}}


def category_label(value) -> str | None:
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return None
    if isinstance(value, float) and value.is_integer():
        return str(int(value))  # a Leiden domain read as 3.0 is "3"
    return str(value)


def bool01(s: pd.Series) -> list[int]:
    if pd.api.types.is_bool_dtype(s):
        return s.astype(int).tolist()
    if pd.api.types.is_numeric_dtype(s):
        return (s.fillna(0) != 0).astype(int).tolist()
    return s.map(lambda v: int(v is True or str(v).strip().lower() in ("true", "1", "yes"))).tolist()


def finite_list(values: np.ndarray, digits: int | None = None) -> list:
    if digits is not None:
        values = np.round(values, digits)
    return [v if math.isfinite(v) else None for v in values.tolist()]


class BinsTable:
    """The per-bin table of one run, with the JSON of GET /bins encoded once per sample."""

    def __init__(self, df: pd.DataFrame, cfg: dict):
        if "sample_id" not in df.columns:
            raise DataError("the per-bin table has no sample_id column")
        self.df = df
        self.cfg = cfg
        sample_ids = df["sample_id"].map(category_label)
        self.groups = dict(sample_ids.groupby(sample_ids, sort=False).indices)
        self.kinds: dict[str, str] = {}
        self.labels: dict[str, list] = {}
        self.categories: dict[str, list[str]] = {}
        for col in df.columns:
            if col in BIN_SKIP_COLUMNS:
                continue
            s = df[col]
            if col == "in_region_of_interest" or col.startswith("positive_") or column_kind(s) == "bool":
                self.kinds[col] = "bool"
            elif col not in BIN_CATEGORY_COLUMNS and pd.api.types.is_numeric_dtype(s):
                self.kinds[col] = "number"
            else:
                labels = [category_label(v) for v in s.tolist()]
                values = {v for v in labels if v is not None}
                if col not in BIN_CATEGORY_COLUMNS and len(values) > BIN_MAX_CATEGORIES:
                    continue
                self.kinds[col] = "category"
                self.labels[col] = labels
                self.categories[col] = self._category_order(col, values)
        self.colors = bin_colors(cfg)
        self._payloads: dict[str, bytes] = {}
        self._lock = threading.Lock()

    def _category_order(self, col: str, values: set[str]) -> list[str]:
        preferred, last = [], ["unassigned"]
        if col == "tissue_class":
            preferred = [c.get("name") for c in self.cfg.get("classes") or [] if isinstance(c, dict)]
            last = ["removed", "unassigned"]
        elif col == "assigned_signature":
            names = [*(self.cfg.get("programs") or []), *(self.cfg.get("signatures") or {})]
            preferred = [signature_slug(self.cfg, n) for n in names if isinstance(n, str)]
            last = ["ambiguous"]
        head = [v for v in dict.fromkeys(preferred) if v in values and v not in last]
        tail = [v for v in last if v in values]
        rest = sorted(values - set(head) - set(tail), key=natural_key)
        return head + rest + tail

    def samples(self) -> list[dict]:
        out, seen = [], set()
        sections = self.df["section"] if "section" in self.df.columns else None
        for s in config_samples(self.cfg):
            if s["id"] in self.groups:
                out.append({"id": s["id"], "label": s.get("label") or s["id"], "n": len(self.groups[s["id"]])})
                seen.add(s["id"])
        for sid, idx in self.groups.items():
            if sid is not None and sid not in seen:
                label = category_label(sections.iloc[idx[0]]) if sections is not None else None
                out.append({"id": sid, "label": label or sid, "n": len(idx)})
        return out

    def payload(self, sid: str) -> bytes:
        with self._lock:
            cached = self._payloads.get(sid)
        if cached is not None:
            return cached
        if sid not in self.groups:
            raise NotFound(f"sample {sid!r} is not in the per-bin table")
        idx = self.groups[sid]
        sub = self.df.iloc[idx]
        fields = {}
        for col, kind in self.kinds.items():
            if kind == "category":
                lookup = {v: i for i, v in enumerate(self.categories[col])}
                labels = self.labels[col]
                fields[col] = {"type": "category", "categories": self.categories[col],
                               "codes": [lookup.get(labels[i], -1) for i in idx.tolist()]}
            elif kind == "bool":
                fields[col] = {"type": "bool", "values": bool01(sub[col])}
            else:
                values = sub[col]
                fields[col] = {"type": "number", "values": (values.tolist() if pd.api.types.is_integer_dtype(values)
                                                            else finite_list(values.to_numpy(float)))}
        coords = {axis: finite_list(pd.to_numeric(sub[f"{axis}_um"], errors="coerce").to_numpy(float), 1)
                  if f"{axis}_um" in sub.columns else [None] * len(idx) for axis in ("x", "y")}
        body = {"sample": sid, "n": len(idx), "x": coords["x"], "y": coords["y"], "fields": fields,
                "colors": self.colors}
        encoded = json.dumps(body, ensure_ascii=False, allow_nan=False, separators=(",", ":")).encode("utf-8")
        with self._lock:
            self._payloads[sid] = encoded
        return encoded


def per_bin_file(results: Path, cfg: dict) -> Path:
    try:
        path = results / f"per_bin_signature_calls_{int(pc.bin_um(cfg.get('bin_size')))}um.csv.gz"
    except ValueError:
        path = None
    if path is not None and path.is_file():
        return path
    found = sorted(results.glob("per_bin_signature_calls_*um.csv.gz")) if results.is_dir() else []
    if found:
        return found[0]
    raise NotFound("this run has no per-bin table (per_bin_signature_calls_*.csv.gz) yet")


# ============================================================================
# RUN RESULTS
# ============================================================================

def scan_results(root: Path) -> dict[str, int]:
    """Every file below results/ (relative POSIX path -> size), except hidden and temporary files."""
    out = {}
    if not root.is_dir():
        return out
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        for name in sorted(filenames):
            if name.startswith(".") or name.endswith(".tmp"):
                continue
            path = os.path.join(dirpath, name)
            try:
                out[Path(os.path.relpath(path, root)).as_posix()] = os.path.getsize(path)
            except OSError:
                continue
    return out


def is_table_file(name: str) -> bool:
    return name.lower().endswith(TABLE_SUFFIXES)


def read_manifest_entries(path: Path, key: str) -> list[dict]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    entries = data.get(key) if isinstance(data, dict) else None
    return [e for e in entries if isinstance(e, dict)] if isinstance(entries, list) else []


def media_type_of(name: str) -> str:
    media_type, encoding = mimetypes.guess_type(name)
    if encoding:
        return _COMPRESSED_TYPES.get(encoding, "application/octet-stream")
    return media_type or "application/octet-stream"


def utf8_prefix(data: bytes) -> int:
    """Length of the longest prefix of `data` that does not end inside a multi-byte UTF-8 character."""
    for back in range(1, min(4, len(data)) + 1):
        byte = data[-back]
        if byte < 0x80:
            return len(data)
        if byte >= 0xC0:
            need = 2 if byte < 0xE0 else 3 if byte < 0xF0 else 4
            return len(data) if back >= need else len(data) - back
    return len(data)


# ============================================================================
# ACCESS TOKEN
# ============================================================================
_DENIED_PAGE = """<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Bone Spatial Studio</title></head>
<body style="font-family:system-ui,sans-serif;max-width:40em;margin:4em auto;line-height:1.5">
<h1>Access token required</h1>
<p>This Bone Spatial Studio server accepts other computers on the network, so it needs its access token. Open the
address printed in the terminal where the server was started; it ends with <code>?token=...</code>.</p>
</body></html>"""


class TokenGate:
    """ASGI middleware: every request must carry the access token, in the X-Token header or the bss_token cookie.
    Visiting /?token=<token> sets the cookie and redirects to /, so the token leaves the address bar."""

    def __init__(self, app, token: str):
        self.app = app
        self.token = token.encode()

    def _valid(self, value: str | None) -> bool:
        return value is not None and secrets.compare_digest(value.encode(), self.token)

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        request = Request(scope)
        query_token = request.query_params.get("token")
        if scope["path"] == "/" and query_token is not None and self._valid(query_token):
            response = RedirectResponse("/", status_code=303)
            response.set_cookie(TOKEN_COOKIE, query_token, max_age=TOKEN_COOKIE_AGE, path="/", httponly=True,
                                samesite="lax")
        elif self._valid(request.headers.get(TOKEN_HEADER)) or self._valid(request.cookies.get(TOKEN_COOKIE)):
            await self.app(scope, receive, send)
            return
        elif scope["path"].startswith("/api/"):
            response = JSONOut({"detail": "access token required: open the address printed by the server (it ends "
                                          "with ?token=...) or send the X-Token header"}, status_code=401)
        else:
            response = HTMLResponse(_DENIED_PAGE, status_code=401)
        await response(scope, receive, send)


LOOPBACK_HOSTS = frozenset({"127.0.0.1", "localhost", "::1"})


def _host_name(value: str) -> str:
    """'127.0.0.1:8765' -> '127.0.0.1', '[::1]:8765' -> '::1' (lower case)."""
    value = value.strip().lower()
    if value.startswith("["):
        return value[1:value.find("]")] if "]" in value else value
    return value.rsplit(":", 1)[0] if value.count(":") == 1 else value


class LocalOriginGate:
    """ASGI middleware of the default, token-less server on 127.0.0.1: only pages served by this server may use the
    API. Without it any website open in the same browser could call it: DNS rebinding (a foreign name resolving to
    127.0.0.1, so the Host header is foreign) reads every response, and cross-site requests with a text/plain body
    create projects or start runs. Rejected: a Host that is not a loopback name; on /api/, an Origin that is not this
    server or a Sec-Fetch-Site other than same-origin/none; a non-empty request body that is not JSON (forces a CORS
    preflight, which this server never grants)."""

    def __init__(self, app, hosts=LOOPBACK_HOSTS):
        self.app = app
        self.hosts = frozenset(hosts)

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        headers = {k.decode("latin-1").lower(): v.decode("latin-1") for k, v in scope.get("headers", [])}
        host = headers.get("host", "")
        reason = None
        if _host_name(host) not in self.hosts:
            reason = "this server only answers requests addressed to 127.0.0.1 or localhost"
        elif scope["path"].startswith("/api/"):
            origin = headers.get("origin")
            site = headers.get("sec-fetch-site")
            if origin is not None and origin.lower() not in (f"http://{host.lower()}", f"https://{host.lower()}"):
                reason = "requests from other websites are refused"
            elif site is not None and site not in ("same-origin", "none"):
                reason = "requests from other websites are refused"
            elif scope["method"] in ("POST", "PUT", "PATCH", "DELETE") and (
                    headers.get("content-length", "0") != "0" or "transfer-encoding" in headers):
                media = headers.get("content-type", "").split(";")[0].strip().lower()
                if media != "application/json":
                    reason = "request bodies must be sent as application/json"
        if reason is None:
            await self.app(scope, receive, send)
            return
        await JSONOut({"detail": reason}, status_code=403)(scope, receive, send)


# ============================================================================
# REQUEST BODIES
# ============================================================================

async def _read_object(request: Request, required: bool) -> dict:
    declared = request.headers.get("content-length", "")
    if declared.isdigit() and int(declared) > MAX_BODY:
        raise HTTPException(413, f"the request body is larger than {MAX_BODY // (1024 * 1024)} MB")
    chunks, size = [], 0
    async for chunk in request.stream():
        size += len(chunk)
        if size > MAX_BODY:
            raise HTTPException(413, f"the request body is larger than {MAX_BODY // (1024 * 1024)} MB")
        chunks.append(chunk)
    raw = b"".join(chunks)
    if not raw.strip():
        if required:
            raise HTTPException(400, "the request body must be a JSON object")
        return {}
    try:
        data = json.loads(raw)
    except ValueError as exc:
        raise HTTPException(400, f"the request body is not valid JSON ({exc})") from None
    if not isinstance(data, dict):
        raise HTTPException(400, "the request body must be a JSON object")
    return data


async def json_body(request: Request) -> dict:
    return await _read_object(request, required=True)


async def optional_json_body(request: Request) -> dict:
    return await _read_object(request, required=False)


def optional_object(body: dict, key: str) -> dict | None:
    value = body.get(key)
    if value is not None and not isinstance(value, dict):
        raise HTTPException(400, f"{key} must be a JSON object")
    return value


def optional_text(body: dict, key: str, max_len: int = 200, allow_empty: bool = True) -> str | None:
    value = body.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or (not allow_empty and not value.strip()):
        raise HTTPException(400, f"{key} must be a{'' if allow_empty else ' non-empty'} string")
    value = value.strip()
    if len(value) > max_len:
        raise HTTPException(400, f"{key} must be at most {max_len} characters")
    return value


# ============================================================================
# APPLICATION
# ============================================================================

def create_app(workspace: Path, token: str | None = None, worker_module: str = "webapp.worker") -> FastAPI:
    """The Bone Spatial Studio application for one workspace folder. With `token`, every request must carry it (see
    TokenGate). worker_module is the module the job runner starts for a run (the tests use a fake one). The job runner
    starts and stops with the application (lifespan)."""
    if token is not None and not TOKEN_RE.match(token):
        raise ValueError("the access token must be 1-256 letters, digits, '.', '_', '~' or '-'")
    storage = Storage(workspace)
    jobs = JobManager(storage.root, worker_module=worker_module, storage=storage)
    tables = LRUCache(6)
    bins_cache = LRUCache(4)
    demo_lock = threading.Lock()

    @contextlib.asynccontextmanager
    async def lifespan(_app: FastAPI):
        jobs.start()
        try:
            yield
        finally:
            jobs.stop()

    app = FastAPI(title="Bone Spatial Studio", version=__version__, lifespan=lifespan, default_response_class=JSONOut,
                  docs_url=None, redoc_url=None, openapi_url=None)
    app.state.storage = storage
    app.state.jobs = jobs
    if token:
        app.add_middleware(TokenGate, token=token)
    else:
        app.add_middleware(LocalOriginGate)

    # ---- errors ------------------------------------------------------------------------------------------------
    @app.exception_handler(StarletteHTTPException)
    async def http_error(_request: Request, exc: StarletteHTTPException):
        return JSONOut({"detail": exc.detail}, status_code=exc.status_code, headers=getattr(exc, "headers", None))

    @app.exception_handler(StorageError)
    async def storage_error(_request: Request, exc: StorageError):
        status = 404 if isinstance(exc, NotFound) else 400 if isinstance(exc, BadRequest) else 500
        if status == 500:
            log.error("workspace error: %s", exc)
        return JSONOut({"detail": str(exc)}, status_code=status)

    @app.exception_handler(DataError)
    async def data_error(_request: Request, exc: DataError):
        return JSONOut({"detail": str(exc)}, status_code=422)

    @app.exception_handler(RequestValidationError)
    async def invalid_request(_request: Request, exc: RequestValidationError):
        parts = []
        for e in exc.errors():
            where = ".".join(str(x) for x in e.get("loc", ()) if x not in ("query", "path", "body"))
            parts.append(f"{where}: {e.get('msg')}" if where else str(e.get("msg")))
        return JSONOut({"detail": "invalid request: " + "; ".join(parts)}, status_code=400)

    @app.exception_handler(Exception)
    async def internal_error(_request: Request, exc: Exception):
        # Starlette re-raises the exception after this response, so the server logs the traceback
        return JSONOut({"detail": f"internal server error: {type(exc).__name__}: {exc}"}, status_code=500)

    # ---- shared helpers --------------------------------------------------------------------------------------------
    def project_view(project: dict) -> dict:
        stored = project.get("config")
        cfg = pc.merge_with_defaults(stored if isinstance(stored, dict) else {})
        errors, warnings = check_config(cfg)
        return {"id": project["id"], "name": project.get("name", ""), "created": project.get("created"),
                "updated": project.get("updated"), "config": cfg, "rois": storage.list_rois(project["id"]),
                "validation": {"errors": errors, "warnings": warnings}}

    def project_sample(pid: str, sid: str) -> tuple[dict, dict]:
        project = storage.get_project(pid)
        cfg = anchored(pc.merge_with_defaults(project.get("config") or {}))
        for s in config_samples(cfg):
            if s["id"] == sid:
                return cfg, s
        raise NotFound(f"sample {sid!r} is not in project {pid}")

    def run_view(rid: str) -> dict:
        record = jobs.get_run(rid)
        try:
            record["config"] = storage.read_run_config(rid)
        except NotFound:
            record["config"] = None
        return record

    def require_run(rid: str) -> Path:
        run_dir = storage.run_dir(rid)
        if not (run_dir / "run.json").is_file():
            raise NotFound(f"run {rid} not found")
        return run_dir

    def run_config(rid: str) -> dict:
        try:
            return pc.merge_with_defaults(storage.read_run_config(rid))
        except NotFound:
            return pc.default_config()

    def results_file(rid: str, rel: str) -> Path:
        """A file inside the run's results/ folder; anything that could leave it is rejected."""
        root = (require_run(rid) / "results").resolve()
        parts = rel.split("/")
        if (not rel or "\\" in rel or ":" in rel or "\0" in rel or rel.startswith("/")
                or any(p in ("", ".", "..") for p in parts)):
            raise HTTPException(400, f"invalid file path {rel!r}")
        path = root.joinpath(*parts).resolve()
        if root not in path.parents:
            raise HTTPException(400, f"invalid file path {rel!r}")
        if not path.is_file():
            raise NotFound(f"{rel} not found in the results of run {rid}")
        return path

    def resolve_mask(cfg: dict, given) -> Path | None:
        p = Path(str(given)).expanduser()
        candidates = [p] if p.is_absolute() else [Path(cfg["data_root"]) / p, PROJECT_ROOT / p]
        return next((c for c in candidates if c.is_file()), None)

    def roi_sources(pid: str, cfg: dict) -> tuple[dict[str, Path], list[str]]:
        """Sample id -> GeoJSON file to copy into a new run (the region drawn in the app, else a file named in
        roi.masks), and the problems that must stop the run."""
        roi = cfg.get("roi") if isinstance(cfg.get("roi"), dict) else {}
        use_drawn = roi.get("use_drawn") is True
        frame = roi.get("mask_frame") or "hires"
        masks = roi.get("masks") if isinstance(roi.get("masks"), dict) else {}
        acfg = anchored(cfg)
        files, problems = {}, []
        for s in config_samples(cfg):
            sid = s["id"]
            try:
                drawn = storage.roi_path(pid, sid)
            except BadRequest:
                continue
            if drawn.is_file():
                files[sid] = drawn
                if not use_drawn:
                    continue
                try:
                    other = sorted(mask_frames(drawn) - {frame})
                except (OSError, ValueError, AttributeError) as exc:
                    problems.append(f"the region drawn for sample {sid} cannot be read ({exc})")
                    continue
                if other:
                    problems.append(f"the region drawn for sample {sid} is in {other[0]} image pixels but "
                                    f"roi.mask_frame is {frame!r}; set the image the regions were drawn on to "
                                    f"{other[0]}")
            elif masks.get(sid):
                path = resolve_mask(acfg, masks[sid])
                if path is not None:
                    files[sid] = path
                elif use_drawn:
                    problems.append(f"roi.masks: the region file of sample {sid} does not exist: {masks[sid]}")
        return files, problems

    def require_runnable(cfg: dict, problems: list[str]) -> None:
        errors = check_config(cfg)[0] + problems
        if errors:
            raise HTTPException(409, "the configuration has errors; fix them before running: " + summarize(errors))
        missing = missing_required()
        if missing:
            raise HTTPException(409, f"required packages are missing: {', '.join(missing)}; install them (see "
                                     f"requirements.txt) and restart the server")

    def submit(project_id: str, project_name: str, label: str, cfg: dict, roi_files: dict[str, Path]) -> dict:
        record = storage.create_run(project_id=project_id, project_name=project_name, label=label, config=cfg,
                                    roi_files=roi_files)
        snapshot = storage.read_run_config(record["id"])
        jobs.submit(record["id"])
        return {**record, "config": snapshot}

    def forget_run(rid: str) -> None:
        prefix = str(storage.runs_dir / rid) + os.sep
        tables.discard(lambda key: key[0].startswith(prefix))
        bins_cache.discard(lambda key: key[0] == rid)

    # ---- front end ---------------------------------------------------------------------------------------------
    app.mount("/static", NoCacheStaticFiles(directory=STATIC_DIR), name="static")

    @app.get("/", include_in_schema=False)
    def index():
        return FileResponse(STATIC_DIR / "index.html", media_type="text/html",
                            headers={"Cache-Control": "no-cache"})

    # ---- environment -------------------------------------------------------------------------------------------
    @app.get("/api/health")
    def health():
        versions = package_versions()
        return JSONOut({
            "ok": True, "app_version": __version__, "python": platform.python_version(), "platform": platform_text(),
            "workspace": str(storage.root), "packages": versions,
            "missing_required": [n for n in REQUIRED_PACKAGES if versions.get(n) is None],
            "missing_optional": [n for n in OPTIONAL_PACKAGES if versions.get(n) is None],
            "active_run": jobs.active_run_id, "queued": jobs.queued_count,
        })

    @app.get("/api/defaults")
    def defaults():
        return JSONOut({"config": pc.default_config(), "param_schema": pc.PARAM_SCHEMA, "step_labels": pc.STEP_LABELS})

    @app.get("/api/docs/{name}")
    def docs(name: str):
        if name not in DOCS:
            raise NotFound(f"no document {name!r}; available: {', '.join(DOCS)}")
        path = PROJECT_ROOT / DOCS[name]
        try:
            text = path.read_text(encoding="utf-8")
        except FileNotFoundError:
            raise NotFound(f"{DOCS[name]} is missing") from None
        return JSONOut({"name": name, "markdown": text})

    # ---- projects ----------------------------------------------------------------------------------------------
    @app.get("/api/projects")
    def list_projects():
        last_runs: dict[str, dict] = {}
        for run in jobs.list_runs():  # newest first
            last_runs.setdefault(run.get("project_id"), run)
        out = []
        for p in storage.list_projects():
            cfg = p.get("config") if isinstance(p.get("config"), dict) else {}
            last = last_runs.get(p["id"])
            out.append({"id": p["id"], "name": p.get("name", ""), "created": p.get("created"),
                        "updated": p.get("updated"),
                        "n_samples": len(cfg["samples"]) if isinstance(cfg.get("samples"), list) else 0,
                        "data_root": cfg.get("data_root"),
                        "last_run": ({"id": last["id"], "status": last.get("status"), "created": last.get("created")}
                                     if last else None)})
        return JSONOut(out)

    @app.post("/api/projects")
    def create_project(body: dict = Depends(json_body)):
        config = optional_object(body, "config")
        name = optional_text(body, "name")
        if not name and config and isinstance(config.get("name"), str):
            name = config["name"].strip()[:200]
        if not name:
            raise HTTPException(400, "a project name is required")
        return JSONOut(project_view(storage.create_project(name, config)))

    @app.get("/api/projects/{pid}")
    def get_project(pid: str):
        return JSONOut(project_view(storage.get_project(pid)))

    @app.put("/api/projects/{pid}")
    def update_project(pid: str, body: dict = Depends(json_body)):
        name = optional_text(body, "name", allow_empty=False)
        config = optional_object(body, "config")
        return JSONOut(project_view(storage.update_project(pid, name=name, config=config)))

    @app.delete("/api/projects/{pid}")
    def delete_project(pid: str):
        storage.delete_project(pid)
        return JSONOut({"ok": True})

    @app.post("/api/projects/{pid}/duplicate")
    def duplicate_project(pid: str):
        return JSONOut(project_view(storage.duplicate_project(pid)))

    @app.post("/api/projects/{pid}/validate")
    def validate_project(pid: str, body: dict = Depends(optional_json_body)):
        project = storage.get_project(pid)
        given = optional_object(body, "config")
        cfg = pc.merge_with_defaults(given if given is not None else project.get("config") or {})
        errors, warnings = check_config(cfg)
        acfg = anchored(cfg)
        reports, panels = [], []
        for sample in config_samples(acfg):
            entry, errs, warns, genes = sample_report(acfg, sample)
            reports.append(entry)
            errors += errs
            warnings += warns
            if genes is not None:
                panels.append(genes)
        genes, gene_warnings = genes_report(cfg, panels)
        warnings += gene_warnings
        errors += roi_sources(pid, cfg)[1]
        return JSONOut({"ok": not errors, "errors": errors, "warnings": warnings, "samples": reports,
                        "genes": genes})

    @app.post("/api/projects/{pid}/detect-samples")
    def detect(pid: str, body: dict = Depends(optional_json_body)):
        project = storage.get_project(pid)
        given = optional_text(body, "data_root", max_len=4096)
        cfg = anchored({"data_root": given} if given else project.get("config") or {})
        root = Path(cfg["data_root"])
        if not root.is_dir():
            raise HTTPException(400, f"not a folder: {root}")
        return JSONOut({"data_root": str(root), "candidates": detect_samples(root)})

    @app.get("/api/fs/list")
    def fs_list(path: str = ""):
        text = path.strip()
        if re.fullmatch(r"[A-Za-z]:", text):
            text += "\\"
        target = Path(text).expanduser() if text else Path.home()
        if not target.is_absolute():
            target = PROJECT_ROOT / target
        target = Path(os.path.abspath(target))
        if not target.exists():
            raise NotFound(f"folder not found: {target}")
        if not target.is_dir():
            raise HTTPException(400, f"not a folder: {target}")
        try:
            children = subfolders(target)
        except PermissionError:
            raise HTTPException(403, f"permission denied: {target}") from None
        except OSError as exc:
            raise HTTPException(400, f"cannot list {target}: {exc}") from None
        children.sort(key=lambda e: e.name.casefold())
        entries = [{"name": e.name, "path": os.path.join(str(target), e.name), "is_dir": True,
                    "looks_like_spaceranger": looks_like_spaceranger(e.path)} for e in children]
        parent = target.parent
        return JSONOut({"path": str(target), "parent": str(parent) if parent != target else None,
                        "roots": filesystem_roots(), "entries": entries})

    # ---- sample images -----------------------------------------------------------------------------------------
    @app.get("/api/projects/{pid}/samples/{sid}/image")
    def sample_image(pid: str, sid: str):
        cfg, sample = project_sample(pid, sid)
        image = hires_image(cfg, sample)
        if image is None:
            raise NotFound(f"sample {sid} has no spatial/tissue_hires_image.png")
        return FileResponse(image, media_type="image/png", headers={"Cache-Control": "no-cache"})

    @app.get("/api/projects/{pid}/samples/{sid}/image-meta")
    def sample_image_meta(pid: str, sid: str):
        cfg, sample = project_sample(pid, sid)
        image = hires_image(cfg, sample)
        if image is None:
            raise NotFound(f"sample {sid} has no spatial/tissue_hires_image.png")
        try:
            with Image.open(image) as im:
                width, height = im.size
        except (OSError, ValueError, Image.DecompressionBombError) as exc:
            raise DataError(f"{image} cannot be read as an image ({exc})") from exc
        scalefactors = read_scalefactors(image.parent) or {}
        return JSONOut({"width": width, "height": height, "frame": "hires", "scalefactors": scalefactors,
                        "microns_per_pixel": microns_per_pixel(sample, scalefactors)})

    @app.get("/api/projects/{pid}/samples/{sid}/bins-preview")
    def bins_preview(pid: str, sid: str):
        cfg, sample = project_sample(pid, sid)
        bin_dir = pc.resolve_bin_dir(cfg, sample)
        positions = positions_file(bin_dir / "spatial") if bin_dir is not None else None
        if positions is None:
            raise NotFound(f"sample {sid} has no tissue positions for bin size {cfg.get('bin_size')}")
        scalefactors = read_scalefactors(positions.parent)
        if scalefactors is None:
            image = hires_image(cfg, sample)
            scalefactors = read_scalefactors(image.parent) if image is not None else None
        scalef = positive_number((scalefactors or {}).get("tissue_hires_scalef"))
        if scalef is None:
            raise NotFound(f"sample {sid} has no tissue_hires_scalef in spatial/scalefactors_json.json")
        x, y = in_tissue_xy(positions)
        n = len(x)
        if n > PREVIEW_MAX_BINS:
            keep = np.unique(np.linspace(0, n - 1, PREVIEW_MAX_BINS).round().astype(np.int64))
            x, y = x[keep], y[keep]
        return JSONOut({"x": finite_list(x * scalef, 1), "y": finite_list(y * scalef, 1), "n_total": n})

    # ---- regions of interest -----------------------------------------------------------------------------------
    @app.get("/api/projects/{pid}/roi/{sid}")
    def get_roi(pid: str, sid: str):
        return JSONOut(storage.read_roi(pid, sid))

    @app.put("/api/projects/{pid}/roi/{sid}")
    def put_roi(pid: str, sid: str, body: dict = Depends(json_body)):
        storage.get_project(pid)
        check_sample_id(sid)
        collection, n_polygons = clean_feature_collection(body)
        if n_polygons:
            storage.write_roi(pid, sid, collection)
        else:
            storage.delete_roi(pid, sid)
        return JSONOut({"ok": True, "n_polygons": n_polygons})

    @app.delete("/api/projects/{pid}/roi/{sid}")
    def delete_roi(pid: str, sid: str):
        storage.delete_roi(pid, sid)
        return JSONOut({"ok": True})

    # ---- runs --------------------------------------------------------------------------------------------------
    @app.post("/api/projects/{pid}/runs")
    def start_run(pid: str, body: dict = Depends(optional_json_body)):
        project = storage.get_project(pid)
        cfg = pc.merge_with_defaults(project.get("config") or {})
        label = optional_text(body, "label") or ""
        steps = optional_object(body, "steps")
        if steps:
            unknown = sorted(set(steps) - set(pc.STEP_LABELS))
            if unknown:
                raise HTTPException(400, f"unknown steps: {', '.join(unknown)}")
            if not all(isinstance(v, bool) for v in steps.values()):
                raise HTTPException(400, "steps values must be true or false")
            if not isinstance(cfg.get("steps"), dict):
                cfg["steps"] = {}
            cfg["steps"].update(steps)
        deterministic = body.get("deterministic")
        if deterministic is not None:
            if not isinstance(deterministic, bool):
                raise HTTPException(400, "deterministic must be true or false")
            cfg["deterministic"] = deterministic
        roi_files, problems = roi_sources(pid, cfg)
        require_runnable(cfg, problems)
        return JSONOut(submit(pid, project.get("name", ""), label, cfg, roi_files))

    @app.get("/api/runs")
    def list_runs(project: str | None = None):
        return JSONOut([{k: v for k, v in r.items() if k != "config"} for r in jobs.list_runs(project or None)])

    @app.get("/api/runs/{rid}")
    def get_run(rid: str):
        return JSONOut(run_view(rid))

    @app.get("/api/runs/{rid}/log")
    def run_log(rid: str, offset: int = 0):
        finished = jobs.get_run(rid).get("status") in FINAL_STATUSES  # before reading: the log is then complete
        try:
            with open(storage.log_path(rid), "rb") as fh:
                size = os.fstat(fh.fileno()).st_size
                start = min(max(0, offset), size)
                fh.seek(start)
                data = fh.read(LOG_CHUNK)
        except FileNotFoundError:
            size, start, data = 0, 0, b""
        end = start + len(data)
        cut = len(data) if finished and end >= size else utf8_prefix(data)
        return JSONOut({"offset": start + cut, "text": data[:cut].decode("utf-8", errors="replace"),
                        "complete": finished and start + cut >= size})

    @app.post("/api/runs/{rid}/cancel")
    def cancel_run(rid: str):
        jobs.cancel(rid)
        return JSONOut(run_view(rid))

    @app.post("/api/runs/{rid}/rerun")
    def rerun(rid: str):
        old = jobs.get_run(rid)
        cfg = pc.merge_with_defaults(storage.read_run_config(rid))
        run_dir = storage.run_dir(rid)
        roi = cfg.get("roi") if isinstance(cfg.get("roi"), dict) else {}
        masks = roi.get("masks") if isinstance(roi.get("masks"), dict) else {}
        roi_files, problems = {}, []
        for sid, given in masks.items():
            kept = run_dir / "rois" / f"{sid}.geojson"  # the copy of this run, wherever the workspace now is
            path = kept if kept.is_file() else Path(str(given))
            if path.is_file():
                roi_files[sid] = path
            else:
                problems.append(f"the region of sample {sid} of run {rid} is missing ({given})")
        require_runnable(cfg, problems)
        project_id = old.get("project_id") or ""
        try:
            project_name = storage.get_project(project_id).get("name", "")
        except NotFound:
            project_name = old.get("project_name") or ""
        return JSONOut(submit(project_id, project_name, old.get("label") or "", cfg, roi_files))

    @app.delete("/api/runs/{rid}")
    def delete_run(rid: str):
        record = jobs.get_run(rid)
        if jobs.is_busy(rid) or record.get("status") in ACTIVE_STATUSES:
            raise HTTPException(409, f"run {rid} is {record.get('status')}; cancel it before deleting it")
        forget_run(rid)
        try:
            storage.delete_run(rid)
        except OSError as exc:
            raise HTTPException(500, f"cannot delete run {rid}; a file may be open in another program ({exc})") \
                from None
        return JSONOut({"ok": True})

    # ---- results -----------------------------------------------------------------------------------------------
    @app.get("/api/runs/{rid}/summary")
    def run_summary(rid: str):
        require_run(rid)
        return JSONOut(storage.read_json(storage.results_dir(rid) / "summary.json", f"summary of run {rid}"))

    @app.get("/api/runs/{rid}/manifest")
    def run_manifest(rid: str):
        results = require_run(rid) / "results"
        files = scan_results(results)
        listed = {"manifest.json", "figures/manifest.json"}
        tables_out, objects, other, figures = [], [], [], []

        def url(rel: str) -> str:
            return f"/api/runs/{rid}/files/{quote(rel, safe='/')}"

        for entry in read_manifest_entries(results / "manifest.json", "files"):
            rel = entry.get("file")
            if not isinstance(rel, str) or rel not in files or rel in listed:
                continue
            listed.add(rel)
            item = {"file": rel, "title": entry.get("title") or rel, "description": entry.get("description") or "",
                    "size": files[rel]}
            rows = entry.get("rows")
            if entry.get("kind") == "object":
                objects.append(item)
            elif entry.get("kind") in ("table", "data") and is_table_file(rel):
                tables_out.append({**item, "rows": rows if isinstance(rows, int) and not isinstance(rows, bool)
                                   else None})
            else:
                other.append({"file": rel, "size": files[rel]})
        for entry in read_manifest_entries(results / "figures" / "manifest.json", "figures"):
            item = dict(entry)
            for key in ("png", "pdf"):
                name = entry.get(key)
                rel = None
                if isinstance(name, str) and name:
                    rel = name if "/" in name else f"figures/{name}"
                item[f"{key}_url"] = url(rel) if rel in files else None
                if rel in files:
                    listed.add(rel)
            figures.append(item)
        for rel, size in files.items():
            if rel in listed:
                continue
            title = rel.rsplit("/", 1)[-1]
            if is_table_file(rel):
                tables_out.append({"file": rel, "title": title, "description": "", "rows": None, "size": size})
            elif rel.lower().endswith(".h5ad"):
                objects.append({"file": rel, "title": title, "description": "", "size": size})
            else:
                other.append({"file": rel, "size": size})
        return JSONOut({"tables": tables_out, "figures": figures, "objects": objects, "other": other})

    @app.get("/api/runs/{rid}/files/{path:path}")
    def run_file(rid: str, path: str, download: str = ""):
        file = results_file(rid, path)
        attach = download.strip().lower() in ("1", "true", "yes")
        return FileResponse(file, media_type=media_type_of(file.name), filename=file.name,
                            content_disposition_type="attachment" if attach else "inline",
                            headers={"Cache-Control": "no-cache"})

    @app.get("/api/runs/{rid}/table")
    def run_table(rid: str, file: str, offset: int = 0, limit: int = 100, sort: str | None = None, desc: str = "0",
                  q: str | None = None):
        path = results_file(rid, file)
        if not is_table_file(path.name):
            raise HTTPException(400, f"{file} is not a CSV or TSV table")
        try:
            table = tables.get(file_key(path), lambda: Table(read_table(path)))
        except (OSError, ValueError) as exc:  # ValueError includes pandas' ParserError and UnicodeDecodeError
            raise DataError(f"{file} cannot be read as a table ({exc})") from exc
        offset = max(0, offset)
        limit = min(max(1, limit), TABLE_MAX_LIMIT)
        n = len(table.df)
        if sort:
            if sort not in table.dtypes:
                raise HTTPException(400, f"unknown column {sort!r}")
            positions = table.order(sort, desc.strip().lower() in ("1", "true", "yes"))
        else:
            positions = np.arange(n)
        if q:
            positions = positions[table.matches(q)[positions]]
        return JSONOut({"file": file, "columns": table.columns, "dtypes": table.dtypes,
                        "rows": table.rows(positions[offset:offset + limit]), "total": int(len(positions)),
                        "offset": offset, "limit": limit})

    def bins_table(rid: str) -> BinsTable:
        results = require_run(rid) / "results"
        cfg = run_config(rid)
        path = per_bin_file(results, cfg)

        def build() -> BinsTable:
            try:
                df = pd.read_csv(path, low_memory=False)
            except (OSError, ValueError) as exc:
                raise DataError(f"{path.name} cannot be read ({exc})") from exc
            return BinsTable(df, cfg)

        return bins_cache.get((rid, *file_key(path)), build)

    @app.get("/api/runs/{rid}/bins/samples")
    def bin_samples(rid: str):
        return JSONOut(bins_table(rid).samples())

    @app.get("/api/runs/{rid}/bins")
    def bins(rid: str, sample: str):
        return Response(bins_table(rid).payload(sample), media_type=JSON_MEDIA_TYPE)

    # ---- demo --------------------------------------------------------------------------------------------------
    @app.post("/api/demo")
    def demo():
        with demo_lock:
            try:
                generator = importlib.import_module("tools.make_demo_data")
                samples = generator.generate(out_dir=storage.demo_dir, size="small")
            except (ImportError, OSError, ValueError) as exc:
                raise HTTPException(500, f"cannot generate the demo data: {exc}") from exc
            for project in storage.list_projects():  # newest first
                if project.get("name") == DEMO_NAME:
                    return JSONOut({"project": project_view(project)})
            config = {"name": DEMO_NAME, "data_root": str(storage.demo_dir), "samples": samples,
                      "reference_sample": samples[0]["id"], **DEMO_SETTINGS}
            return JSONOut({"project": project_view(storage.create_project(DEMO_NAME, config))})

    return app
