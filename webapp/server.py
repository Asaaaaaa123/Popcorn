"""
webapp/server.py
================
FastAPI application of Bone Spatial Studio (docs/ARCHITECTURE.md, section 5): the JSON API under /api/, the front end
under /static/ and its index.html at /.

The server process stays light and responsive: it never imports scanpy, anndata or bone_analysis (an analysis runs in
a worker subprocess started by webapp.jobs). Package versions come from the installed package metadata, the data
checks read the Space Ranger files directly (h5py for the feature-barcode matrix, pandas/pyarrow for the tissue
positions, Pillow for the image size), and every JSON response is strict JSON: numpy values become plain numbers,
NaN and infinities become null. The data explorer and the gene maps of the viewer read the raw counts the same way
(h5py + scipy.sparse), exactly as bone_analysis.read_visium_hd loads them.
"""
from __future__ import annotations

import collections
import contextlib
import functools
import hashlib
import hmac
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
import tempfile
import threading
import time
import zipfile
from pathlib import Path
from typing import NamedTuple
from urllib.parse import quote

import h5py
import numpy as np
import pandas as pd
from fastapi import Depends, FastAPI, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from PIL import Image
from scipy import sparse
from scipy.spatial import cKDTree
from starlette.background import BackgroundTask
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

# raw sections (data explorer, gene maps of the viewer)
RAW_SECTIONS_CACHED = 8            # sections whose counts are held in memory at most, and their estimated size at
RAW_SECTIONS_BUDGET = 1 << 30      # most (19,000 genes: about 100 MB for 105,000 16 um bins, 200 MB for 412,000 8 um)
RAW_SECTION_MAX_BINS = 1_000_000   # bins of one section at most (002um sections have millions: the server stays light)
EXPLORE_PAYLOADS_CACHED = 8        # encoded bins payloads of the explorer held at most, and their size at most
EXPLORE_PAYLOADS_BUDGET = 128 << 20  # (about 5 MB at 16 um, 20 MB at 8 um)
GENE_LIST_LIMIT = 50               # genes per gene list by default
GENE_LIST_MAX = 500                # genes per gene list at most
MAX_GENES_PER_MAP = 50             # genes summed into one gene map at most
MAX_GENE_NAME = 200
GENE_EXPRESSION = "Gene Expression"
MT_PREFIX = "mt-"                  # mitochondrial genes, case-insensitive (bone_analysis.qc_metrics)
BARCODE_SUFFIX_RE = re.compile(r"-\d+$")
IMAGE_KEY_RE = re.compile(r"^[0-9a-f]{32}\Z")
IMAGE_KEYS_MAX = 256

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
              "adipocyte_genes": {"missing": []}, "ubiquitous": {"missing": []}, "lesion_genes": {"missing": []}}
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
    # the Leiden screen's ubiquitous transcripts (Signatures tab) and the A549 lesion genes of the sensitivity analyses
    ubiquitous = gene_list(cfg.get("ubiquitous"))
    report["ubiquitous"]["missing"] = miss = [g for g in ubiquitous if g not in present]
    if ubiquitous and len(miss) == len(ubiquitous):
        warnings.append(f"ubiquitous transcripts: none of the {len(ubiquitous)} genes is in the probe set, so the "
                        f"Leiden screen marks no domain unresolved: list abundant housekeeping transcripts of your data "
                        f"(Signatures tab, Ubiquitous transcripts)")
    elif miss:
        warnings.append(f"ubiquitous transcripts: {len(miss)} of {len(ubiquitous)} genes missing from the probe set "
                        f"({', '.join(miss)})")
    sa = cfg.get("sensitivity_analyses") if isinstance(cfg.get("sensitivity_analyses"), dict) else {}
    lesion = gene_list((sa.get("lesion") if isinstance(sa.get("lesion"), dict) else {}).get("genes"))
    report["lesion_genes"]["missing"] = miss = [g for g in lesion if g not in present]
    steps = cfg.get("steps") if isinstance(cfg.get("steps"), dict) else {}
    if steps.get("sensitivity_analyses") is True and lesion and len(miss) == len(lesion):
        warnings.append(f"A549 lesion flag: none of its genes ({', '.join(lesion)}) is in the probe set, so the "
                        f"sensitivity analyses flag no bin as lesion")
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
    """Small thread-safe LRU mapping; values are built outside the lock, one build per key at a time (a request for a
    key being built waits for that build and gets its value, or its error, instead of building it again; an error is
    not kept, so a later request tries again). At most `size` values are kept and, with `weigh` (value -> bytes), at
    most `budget` bytes of them; the newest value is kept whatever its weight."""

    def __init__(self, size: int, budget: int | None = None, weigh=None):
        self.size, self.budget, self.weigh = size, budget, weigh
        self._items: collections.OrderedDict = collections.OrderedDict()  # key -> (value, weight)
        self._weight = 0
        self._lock = threading.Lock()
        # key -> [lock held while building, number of requests using it, failed builds, the last build's error]
        self._gates: dict = {}

    def _cached(self, key):
        with self._lock:
            if key in self._items:
                self._items.move_to_end(key)
                return True, self._items[key][0]
        return False, None

    def _store(self, key, value) -> None:
        weight = int(self.weigh(value)) if self.weigh is not None else 0
        with self._lock:
            old = self._items.pop(key, None)
            if old is not None:
                self._weight -= old[1]
            self._items[key] = (value, weight)
            self._weight += weight
            while len(self._items) > 1 and (len(self._items) > self.size
                                            or (self.budget is not None and self._weight > self.budget)):
                self._weight -= self._items.popitem(last=False)[1][1]

    def get(self, key, build):
        found, value = self._cached(key)
        if found:
            return value
        with self._lock:
            gate = self._gates.setdefault(key, [threading.Lock(), 0, 0, None])
            gate[1] += 1
            failures = gate[2]
        try:
            with gate[0]:
                found, value = self._cached(key)
                if found:
                    return value
                with self._lock:
                    error = gate[3] if gate[2] != failures else None
                if error is not None:  # the build this request waited for failed: the same answer, not another try
                    raise error
                try:
                    value = build()
                except Exception as exc:
                    with self._lock:
                        gate[2] += 1
                        gate[3] = exc
                    raise
                self._store(key, value)
                return value
        finally:
            with self._lock:
                gate[1] -= 1
                if gate[1] == 0:
                    self._gates.pop(key, None)

    def discard(self, predicate) -> None:
        with self._lock:
            for key in [k for k in self._items if predicate(k)]:
                self._weight -= self._items.pop(key)[1]


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
    """The per-bin table of one run, with the JSON of GET /bins encoded once per sample. `source` identifies the file
    it was read from (cache keys of what is derived from it)."""

    def __init__(self, df: pd.DataFrame, cfg: dict, source: tuple = ()):
        if "sample_id" not in df.columns:
            raise DataError("the per-bin table has no sample_id column")
        self.df = df
        self.cfg = cfg
        self.source = source
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

    def barcodes(self, sid: str) -> list[str]:
        """Barcodes of the bins of a sample, in the order of GET /bins."""
        if sid not in self.groups:
            raise NotFound(f"sample {sid!r} is not in the per-bin table")
        if "barcode" not in self.df.columns:
            raise DataError("the per-bin table has no barcode column")
        values = self.df["barcode"].iloc[self.groups[sid]].tolist()
        return ["" if v is None or (not isinstance(v, str) and pd.isna(v)) else str(v) for v in values]


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
# RAW SECTIONS (data explorer: POST /explore/...; gene maps of the viewer: /runs/{rid}/genes, /gene, /bins/image)
# ============================================================================

def decode_all(values) -> list[str]:
    return [v.decode("utf-8", "replace") if isinstance(v, bytes) else str(v) for v in values]


def norm_barcodes(values) -> pd.Series:
    """bone_analysis._norm_bc of every barcode: surrounding blanks and the -<n> suffix removed."""
    return pd.Series(list(values), dtype=object).astype(str).str.strip().str.replace(BARCODE_SUFFIX_RE, "",
                                                                                      regex=True)


def make_unique(values: list[str]) -> list[str]:
    """anndata's make_index_unique (var_names_make_unique, obs_names_make_unique): the second and later copies of a
    name become <name>-1, <name>-2, ..., skipping names already in use."""
    seen, dups = set(), []
    for i, v in enumerate(values):
        if v in seen:
            dups.append(i)
        seen.add(v)
    if not dups:
        return list(values)
    out, counter = list(values), collections.Counter()
    for i in dups:
        v = values[i]
        while True:
            counter[v] += 1
            name = f"{v}-{counter[v]}"
            if name not in seen:
                seen.add(name)
                out[i] = name
                break
    return out


def too_many_bins(path: Path, n: int) -> DataError:
    return DataError(f"{path} holds {n:,} bins; the explorer and the gene maps read at most "
                     f"{RAW_SECTION_MAX_BINS:,} bins per section: choose a larger bin size, such as 008um or 016um")


def read_counts(path: Path, max_bins: int | None = None) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    """(bins x genes counts, barcodes, gene names) of a feature-barcode matrix as scanpy.read_10x_h5 reads it: the
    'Gene Expression' features only, explicit zeros removed. Names are not made unique yet. A matrix of more than
    max_bins barcodes is refused (DataError) before its counts are read."""
    try:
        with h5py.File(path, "r") as f:
            m = f["matrix"]
            n_features, n_barcodes = (int(v) for v in m["shape"][()][:2])
            if max_bins is not None and n_barcodes > max_bins:
                raise too_many_bins(path, n_barcodes)
            barcodes = decode_all(m["barcodes"][()])
            features = m["features"]
            names = decode_all(features["name"][()])
            kinds = decode_all(features["feature_type"][()]) if "feature_type" in features else None
            data, indices, indptr = m["data"][()], m["indices"][()], m["indptr"][()]
        if data.dtype.kind not in "iu":
            if data.dtype.kind != "f" or not np.array_equal(data, np.round(data)):
                raise ValueError(f"the counts are not integers ({data.dtype})")
            data = data.astype(np.int64)
        X = sparse.csr_matrix((data, indices, indptr), shape=(n_barcodes, n_features))
    except (OSError, KeyError, ValueError, IndexError, TypeError) as exc:
        raise DataError(f"{path} is not a readable feature-barcode matrix ({exc})") from exc
    if len(barcodes) != n_barcodes or len(names) != n_features:
        raise DataError(f"{path}: the barcodes and features do not match the matrix shape")
    if kinds is not None and any(k != GENE_EXPRESSION for k in kinds):
        keep = np.array([k == GENE_EXPRESSION for k in kinds])
        X = X[:, keep]
        names = [n for n, k in zip(names, kinds) if k == GENE_EXPRESSION]
    X.eliminate_zeros()
    return X, barcodes, names


def read_positions(path: Path, normalise: bool = True) -> pd.DataFrame:
    """Every row of a tissue positions file (in tissue or not): barcode (normalised like the matrix barcodes, or as
    written with normalise=False) and the full-resolution pixel coordinates x = pxl_col_in_fullres,
    y = pxl_row_in_fullres (NaN when not a number)."""
    try:
        if path.suffix == ".parquet":
            df = pd.read_parquet(path)
        elif path.name == "tissue_positions_list.csv":  # Space Ranger 1.x: no header
            df = pd.read_csv(path, header=None, names=POSITION_COLUMNS)
        else:
            df = pd.read_csv(path)
        bc = "barcode" if "barcode" in df.columns else df.columns[0]
        barcodes = norm_barcodes(df[bc]) if normalise else df[bc].astype(str)
        return pd.DataFrame({"barcode": barcodes.to_numpy(object),
                             "x": pd.to_numeric(df["pxl_col_in_fullres"], errors="coerce").to_numpy(float),
                             "y": pd.to_numeric(df["pxl_row_in_fullres"], errors="coerce").to_numpy(float)})
    except (OSError, KeyError, ValueError, IndexError, ImportError) as exc:
        raise DataError(f"{path} is not a readable tissue positions table ({exc})") from exc


def section_key(matrix: Path, positions: Path) -> tuple:
    try:
        return (*file_key(matrix), *file_key(positions))
    except OSError as exc:
        raise DataError(f"the section files cannot be read ({exc})") from exc


class RawSection:
    """The raw counts of one section as bone_analysis.read_visium_hd loads them: the barcodes of
    filtered_feature_bc_matrix.h5 without their -<n> suffix (made unique), joined to the tissue positions by barcode,
    bins without coordinates dropped, matrix order kept; the 'Gene Expression' features with their names made unique.
    Per-bin totals as bone_analysis.qc_metrics computes them (mitochondrial genes: names starting with mt-, any case).
    The counts are held by gene (CSC) for the gene maps. Sections of more than RAW_SECTION_MAX_BINS bins are refused."""

    def __init__(self, matrix: Path, positions: Path):
        self.key = section_key(matrix, positions)
        self.matrix_path, self.positions_path = matrix, positions
        X, barcodes, names = read_counts(matrix, RAW_SECTION_MAX_BINS)
        barcodes = make_unique(norm_barcodes(barcodes).tolist())
        tp = read_positions(positions)
        tp = tp[tp["barcode"].isin(barcodes)].drop_duplicates("barcode").set_index("barcode").reindex(barcodes)
        ok = (tp["x"].notna() & tp["y"].notna()).to_numpy()
        self.n_barcodes = len(barcodes)
        if self.n_barcodes and ok.mean() < 0.5:  # the same check as bone_analysis.read_visium_hd
            raise DataError(f"only {int(ok.sum()):,} of {self.n_barcodes:,} barcodes of {matrix} have a position in "
                            f"{positions}; are both files from the same run and bin size?")
        rows = np.flatnonzero(ok)
        self.barcodes = [barcodes[i] for i in rows.tolist()]
        self.index = pd.Index(self.barcodes)
        self.x_px = tp["x"].to_numpy(float)[rows]
        self.y_px = tp["y"].to_numpy(float)[rows]
        X = X[rows] if len(rows) < X.shape[0] else X
        self.n_bins = X.shape[0]
        self.genes = make_unique(names)
        self.total = np.asarray(X.sum(axis=1, dtype=np.int64)).ravel()
        self.n_genes = np.diff(X.indptr).astype(np.int64)
        mt = np.array([g.lower().startswith(MT_PREFIX) for g in self.genes], dtype=bool)
        self.mt = (np.asarray(X[:, mt].sum(axis=1, dtype=np.int64)).ravel() if mt.any()
                   else np.zeros(self.n_bins, np.int64))
        X = X.tocsc()
        if X.nnz and X.data.min() >= 0 and X.data.max() <= np.iinfo(np.uint16).max:
            X.data = X.data.astype(np.uint16)  # half the memory of the int32 counts; sums are taken in int64
        self.X = X
        self.gene_total = np.asarray(X.sum(axis=0, dtype=np.int64)).ravel()
        self.gene_bins = np.diff(X.indptr).astype(np.int64)
        self.folded = [g.casefold() for g in self.genes]
        self._exact = {g: i for i, g in enumerate(self.genes)}
        self._by_folded: dict[str, int] = {}
        for i, g in enumerate(self.folded):
            self._by_folded.setdefault(g, i)
        # estimated memory (the section cache's budget): the counts, the per-bin and per-gene arrays, and about 300
        # bytes per barcode (string and index) and per gene name (three strings and two dictionary entries)
        arrays = (X.data, X.indices, X.indptr, self.x_px, self.y_px, self.total, self.n_genes, self.mt,
                  self.gene_total, self.gene_bins)
        self.nbytes = sum(a.nbytes for a in arrays) + 300 * (self.n_bins + len(self.genes))

    def resolve(self, names: list[str]) -> tuple[list[str], list[str], list[int]]:
        """(genes used, names not found, their columns): an exact name, else the first gene of that name in any case;
        each gene counted once."""
        used, missing, cols = [], [], []
        for name in names:
            col = self._exact.get(name)
            if col is None:
                col = self._by_folded.get(name.casefold())
            if col is None:
                missing.append(name)
            elif col not in cols:
                cols.append(col)
                used.append(self.genes[col])
        return used, missing, cols

    def counts(self, cols: list[int]) -> np.ndarray:
        """Per-bin sum of the counts of the given genes (int64, section order)."""
        if not cols:
            return np.zeros(self.n_bins, np.int64)
        return np.asarray(self.X[:, cols].sum(axis=1, dtype=np.int64)).ravel()

    def gene_stats(self, rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(total UMI, number of bins detecting it) of every gene over the given rows (a row listed twice counts
        twice)."""
        weights = np.bincount(rows, minlength=self.n_bins).astype(np.int64)
        genes = np.repeat(np.arange(len(self.genes)), np.diff(self.X.indptr))
        w = weights[self.X.indices]
        total = np.bincount(genes, weights=self.X.data * w, minlength=len(self.genes))
        detected = np.bincount(genes, weights=w, minlength=len(self.genes))
        return np.rint(total).astype(np.int64), np.rint(detected).astype(np.int64)

    def rows_of(self, barcodes: list[str], sid: str) -> np.ndarray:
        """Row of each barcode of a per-bin table (-1 when it has none); see run_barcode_rows."""
        return run_barcode_rows(self.index, barcodes, sid)


def run_barcode_rows(index: pd.Index, barcodes: list[str], sid: str) -> np.ndarray:
    """Position in `index` (a section's matrix barcodes, normalised and made unique) of each barcode of a per-bin
    table, -1 when it has none. The analysis joins the sections with anndata.concat(index_unique="-"), so a per-bin
    barcode is <normalised matrix barcode>-<sample id>; a barcode without that suffix is looked up as it is, then
    normalised like the matrix barcodes."""
    suffix = f"-{sid}"
    rows = index.get_indexer([b[:-len(suffix)] if b.endswith(suffix) else b for b in barcodes])
    for attempt in (lambda bs: bs, lambda bs: norm_barcodes(bs).tolist()):
        missing = np.flatnonzero(rows < 0)
        if not len(missing):
            break
        rows[missing] = index.get_indexer(attempt([barcodes[i] for i in missing.tolist()]))
    return rows


class RunGeneLink:
    """The bins of one sample of a run (GET /bins order) matched by barcode to the rows of its raw section. It keeps
    no reference to the section, so that the section cache alone bounds the memory."""

    def __init__(self, section: RawSection, barcodes: list[str], sid: str):
        self.rows = section.rows_of(barcodes, sid)
        self.n = len(self.rows)
        self.n_unmatched = int((self.rows < 0).sum())
        self._stats: tuple[np.ndarray, np.ndarray] | None = None
        self._lock = threading.Lock()

    def counts(self, section: RawSection, cols: list[int]) -> np.ndarray:
        values = np.zeros(self.n, np.int64)
        ok = self.rows >= 0
        values[ok] = section.counts(cols)[self.rows[ok]]
        return values

    def gene_stats(self, section: RawSection) -> tuple[np.ndarray, np.ndarray]:
        with self._lock:
            if self._stats is None:
                self._stats = section.gene_stats(self.rows[self.rows >= 0])
            return self._stats


class ImageKeys:
    """Opaque keys of the H&E images the maps may show (GET /api/explore-images/{key}). A key is issued only for an
    image the server resolved itself, never for a path a client sent; it stays the same while the file is unchanged
    (so a cached map keeps a valid address: every request that returns it issues it again)."""

    def __init__(self, size: int = IMAGE_KEYS_MAX):
        self.size = size
        self._secret = secrets.token_bytes(16)
        self._paths: collections.OrderedDict[str, Path] = collections.OrderedDict()
        self._lock = threading.Lock()

    def issue(self, path: Path) -> str:
        key = hmac.new(self._secret, repr(file_key(path)).encode("utf-8"), hashlib.sha256).hexdigest()[:32]
        with self._lock:
            self._paths[key] = path
            self._paths.move_to_end(key)
            while len(self._paths) > self.size:
                self._paths.popitem(last=False)
        return key

    def path(self, key: str) -> Path | None:
        with self._lock:
            return self._paths.get(key)


def section_files(cfg: dict, sample: dict) -> tuple[Path, Path]:
    """(filtered_feature_bc_matrix.h5, tissue positions file) of a sample at the analysis bin size; NotFound when
    either is missing."""
    sid = sample["id"]
    bin_dir = pc.resolve_bin_dir(cfg, sample)
    if bin_dir is None:
        raise NotFound(f"sample {sid}: no square_{cfg.get('bin_size')}/filtered_feature_bc_matrix.h5 in "
                       f"{sample_base(cfg, sample)}")
    positions = positions_file(bin_dir / "spatial")
    if positions is None:
        raise NotFound(f"sample {sid}: no tissue positions file in {bin_dir / 'spatial'}")
    return bin_dir / "filtered_feature_bc_matrix.h5", positions


def section_mpp(sample: dict, spatial: Path) -> float:
    """Microns per full-resolution pixel of a section, as the analysis takes it (the configured value, else
    scalefactors_json.json of the bin size)."""
    mpp = microns_per_pixel(sample, read_scalefactors(spatial))
    if mpp is None:
        raise DataError(f"sample {sample['id']}: microns_per_pixel is not set and "
                        f"{spatial / 'scalefactors_json.json'} does not record it; set it in the sample's settings")
    return mpp


def section_image(cfg: dict, sample: dict, mpp: float, keys: ImageKeys) -> tuple[dict | None, str | None]:
    """({"url", "px_per_um", "width", "height"} of the hires H&E image, None) or (None, why there is none).
    px_per_um = tissue_hires_scalef / microns per pixel: image pixels per micrometre of the bin coordinates."""
    image = hires_image(cfg, sample)
    if image is None:
        return None, "no H&E image (spatial/tissue_hires_image.png)"
    try:
        scalefactors = read_scalefactors(image.parent) or {}
    except DataError as exc:
        return None, str(exc)
    scalef = positive_number(scalefactors.get("tissue_hires_scalef"))
    if scalef is None:
        return None, (f"no tissue_hires_scalef in {image.parent / 'scalefactors_json.json'}; the H&E image cannot be "
                      f"placed")
    try:
        with Image.open(image) as im:
            width, height = im.size
        key = keys.issue(image)
    except (OSError, ValueError, Image.DecompressionBombError) as exc:
        return None, f"{image} cannot be read as an image ({exc})"
    return {"url": f"/api/explore-images/{key}", "px_per_um": scalef / mpp, "width": width, "height": height}, None


def qc_thresholds(cfg: dict) -> tuple[dict, list[str]]:
    """The qc thresholds of a (draft) configuration, a default for any that is not a number: (thresholds, notes)."""
    qc = cfg.get("qc") if isinstance(cfg.get("qc"), dict) else {}
    defaults = pc.default_config()["qc"]
    out, notes = {}, []
    for key in ("min_genes", "min_counts", "max_pct_mt"):
        value = qc.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            value = defaults[key]
            notes.append(f"qc.{key} is not a number; the default {value} is used")
        out[key] = value
    return out, notes


def qc_pct_mt(mt: np.ndarray, total: np.ndarray) -> np.ndarray:
    """Mitochondrial share per bin in percent, float32 like scanpy's calculate_qc_metrics (NaN without counts)."""
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return (mt.astype(np.float32) / total.astype(np.float32)) * np.float32(100)


def passes_qc(n_genes: np.ndarray, total: np.ndarray, pct_mt: np.ndarray, qc: dict) -> np.ndarray:
    """bone_analysis.qc_filter's decision per bin, to the last bin: genes > qc.min_genes, UMIs > qc.min_counts and
    mitochondrial share < qc.max_pct_mt, in float32 (Python-float thresholds: compared in float32 as pandas does there,
    and no overflow for a huge draft value)."""
    with np.errstate(invalid="ignore", over="ignore"):
        return ((n_genes > float(qc["min_genes"])) & (total.astype(np.float32) > float(qc["min_counts"]))
                & (pct_mt < float(qc["max_pct_mt"])))


def explore_payload(section: RawSection, sid: str, bin_size: str, mpp: float, qc: dict, image: dict | None,
                    notes: list[str]) -> bytes:
    """JSON of POST /explore/{sid}/bins: the bin payload of GET /runs/{rid}/bins (same shape) for the raw section."""
    total = section.total
    pct_mt = qc_pct_mt(section.mt, total)
    passes = passes_qc(section.n_genes, total, pct_mt, qc)
    n = section.n_bins
    finite_pct = pct_mt[np.isfinite(pct_mt)]
    n_pass = int(passes.sum())
    summary = {
        "n_bins": n, "n_barcodes": section.n_barcodes, "n_without_position": section.n_barcodes - n,
        "n_features": len(section.genes), "total_umi": int(total.sum()),
        "median_umi": float(np.median(total)) if n else None, "mean_umi": round(float(total.mean()), 2) if n else None,
        "median_genes": float(np.median(section.n_genes)) if n else None,
        "median_pct_mt": round(float(np.median(finite_pct)), 3) if len(finite_pct) else None,
        "n_pass_qc": n_pass, "pct_pass_qc": round(100.0 * n_pass / n, 2) if n else None, "bin_size": bin_size,
        "microns_per_pixel": mpp, "qc": qc, "matrix": str(section.matrix_path),
        "positions": str(section.positions_path),
    }
    body = {
        "sample": sid, "n": n, "x": finite_list(section.x_px * mpp, 2), "y": finite_list(section.y_px * mpp, 2),
        "fields": {"total_umi": {"type": "number", "values": total.tolist()},
                   "n_genes": {"type": "number", "values": section.n_genes.tolist()},
                   "pct_mt": {"type": "number", "values": finite_list(pct_mt.astype(np.float64), 2)},
                   "passes_qc": {"type": "bool", "values": passes.astype(np.int8).tolist()}},
        "colors": {}, "barcodes": section.barcodes, "summary": json_safe(summary), "image": json_safe(image),
        "notes": notes,
    }
    return json.dumps(body, ensure_ascii=False, allow_nan=False, separators=(",", ":")).encode("utf-8")


def gene_table(section: RawSection, total: np.ndarray, detected: np.ndarray, n_bins: int, query: str,
               limit: int) -> dict:
    """Genes by total UMI (descending, then panel order) whose name contains `query` (any case), at most `limit`."""
    order = np.lexsort((np.arange(len(total)), -total)).tolist()
    if query:
        q = query.casefold()
        order = [i for i in order if q in section.folded[i]]
    genes = [{"gene": section.genes[i], "total_umi": int(total[i]), "n_bins_detected": int(detected[i]),
              "pct_bins_detected": round(100.0 * int(detected[i]) / n_bins, 2) if n_bins else None}
             for i in order[:limit]]
    return {"n_bins": n_bins, "n_features": len(section.genes), "n_matches": len(order), "genes": genes}


def sparse_vector(values: np.ndarray, used: list[str], missing: list[str]) -> dict:
    """A gene map: per-bin counts as the positions and values of the bins with a count above zero."""
    idx = np.flatnonzero(values > 0)
    return {"genes_used": used, "genes_missing": missing, "n": int(len(values)), "indices": idx.tolist(),
            "values": values[idx].tolist(), "total": int(values.sum()), "n_detected": int(len(idx)),
            "max": int(values.max()) if len(values) else 0}


def requested_genes(value) -> list[str]:
    if not isinstance(value, list) or not value:
        raise HTTPException(400, "genes must be a non-empty list of gene names")
    if len(value) > MAX_GENES_PER_MAP:
        raise HTTPException(400, f"at most {MAX_GENES_PER_MAP} genes can be summed into one map")
    names = []
    for v in value:
        if not isinstance(v, str) or not v.strip() or len(v.strip()) > MAX_GENE_NAME:
            raise HTTPException(400, f"every gene must be a name of 1 to {MAX_GENE_NAME} characters")
        names.append(v.strip())
    return list(dict.fromkeys(names))


def gene_query(query, limit) -> tuple[str, int]:
    if query is None:
        query = ""
    if not isinstance(query, str) or len(query) > MAX_GENE_NAME:
        raise HTTPException(400, f"query must be a text of at most {MAX_GENE_NAME} characters")
    if limit is None:
        limit = GENE_LIST_LIMIT
    if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= GENE_LIST_MAX:
        raise HTTPException(400, f"limit must be an integer from 1 to {GENE_LIST_MAX}")
    return query.strip(), limit


# ============================================================================
# LOUPE BROWSER EXPORT (GET /runs/{rid}/export/loupe)
# ============================================================================
# The calls of a run as category CSV files for 10x Genomics Loupe Browser, which imports one into the .cloupe file of
# one section at one bin size. The first column, "barcode", lists every barcode of that section's
# filtered_feature_bc_matrix.h5 (the bins of its .cloupe file) exactly as Space Ranger wrote it, -1 included; every
# other column is a category (its header is the name Loupe shows) of text labels. No cell is empty: an analysed bin
# without a value is "NA", a bin the run did not analyse "not analyzed". Comma-separated, LF line endings, UTF-8
# without a byte-order mark, and no comma, quote, tab or line break inside a label (10x documents no quoting).
LOUPE_BIN_UM = 8                    # Space Ranger writes the .cloupe file of a Visium HD section at 8 um
LOUPE_BIN_SIZE = "008um"            # (outs/cloupe_008um.cloupe); other sizes only with --custom-bin-size
LOUPE_NOT_ANALYZED = "not analyzed"
LOUPE_MISSING = "NA"
LOUPE_MAX_VALUES = 32_768           # distinct values of one category Loupe Browser supports (10x)
LOUPE_CHUNK_ROWS = 50_000           # CSV rows formatted at a time
LOUPE_QC_CHUNK = 1_000_000          # stored counts of the 8 um matrix read at a time for its quality control
LOUPE_QC_HEADER = "8 um QC"         # the 8 um bin's own quality control (pass / fail) in the 8 um files
LOUPE_HEADERS = {"tissue_class": "Tissue class", "assigned_signature": "Assigned signature",
                 "leiden_domain": "Leiden domain", "in_region_of_interest": "Region of interest"}
LOUPE_LEFT_OUT = {"sample_id": "the section (one file per section)", "section": "the section (one file per section)",
                  "x_um": "a coordinate (Loupe Browser has the bin positions)",
                  "y_um": "a coordinate (Loupe Browser has the bin positions)"}
HD_BARCODE_RE = re.compile(r"^s_(\d{3})um_(\d{1,6})_(\d{1,6})(?:-\d+)?$")  # Visium HD bin: size, row, column
GRID_KEY = 1_000_000                # key of a grid position: row * GRID_KEY + column
PER_BIN_RE = re.compile(r"^per_bin_signature_calls_(\d+)um\.csv\.gz$")
_LOUPE_UNSAFE = str.maketrans({",": ";", '"': "'", "\r": " ", "\n": " ", "\t": " "})


def loupe_text(value) -> str:
    """A label or category name as Loupe Browser can read it: commas, quotes, tabs and line breaks replaced, no
    surrounding blanks."""
    return str(value).translate(_LOUPE_UNSAFE).strip()


def loupe_header(col: str) -> str:
    """Category name of a per-bin column: tissue_class -> "Tissue class", positive_osteoblast -> "Positive
    osteoblast"."""
    text = LOUPE_HEADERS.get(col) or loupe_text(col.replace("_", " "))
    return text[:1].upper() + text[1:]


def yes_no_labels(col: str) -> tuple[str, str]:
    """The labels of True and False in a yes/no column of the per-bin table."""
    if col.startswith("positive_"):
        return "positive", "negative"
    if col.startswith("in_"):
        return "inside", "outside"
    return "yes", "no"


class LoupeColumn(NamedTuple):
    header: str                     # category name in Loupe Browser
    source: str                     # column of the per-bin table
    kind: str                       # "bool" | "text"
    categories: list[str]           # its labels
    codes: np.ndarray               # per row of the per-bin table: index into categories, -1 = no value


def loupe_columns(df: pd.DataFrame) -> tuple[list[LoupeColumn], list[tuple[str, str]]]:
    """The categories of a per-bin table, whatever its signatures and steps: tissue class, assigned signature, Leiden
    domain and region of interest, then every other yes/no or text column in table order (one positivity call per
    signature, ...), and the columns left out with the reason. Numbers are left out: each value would be a cluster."""
    columns, left_out, used = [], [], {"barcode"}
    for col in [c for c in LOUPE_HEADERS if c in df.columns] + [c for c in df.columns if c not in LOUPE_HEADERS]:
        if col == "barcode":
            continue
        if col in LOUPE_LEFT_OUT:
            left_out.append((col, LOUPE_LEFT_OUT[col]))
            continue
        s = df[col]
        if not s.notna().any():
            left_out.append((col, "no values in this run"))
            continue
        kind = "text" if col in BIN_CATEGORY_COLUMNS else column_kind(s)
        if kind == "number":
            left_out.append((col, "a number, not a category"))
            continue
        if kind == "bool":
            yes, no = yes_no_labels(col)
            labels = [(yes if v else no) if isinstance(v, (bool, np.bool_)) else None for v in s.tolist()]
        else:
            labels = []
            for v in s.tolist():
                text = category_label(v)
                text = loupe_text(text) if text is not None else ""
                if col == "leiden_domain" and text.isdigit():
                    text = f"domain {text}"  # a word, like Loupe's own "Cluster 1": never read as a number
                labels.append(text or None)
        codes, uniques = pd.factorize(pd.Series(labels, dtype=object))
        limit = LOUPE_MAX_VALUES - 2 if col in BIN_CATEGORY_COLUMNS else BIN_MAX_CATEGORIES
        if not len(uniques):
            left_out.append((col, "no values in this run"))
            continue
        if len(uniques) > limit:
            left_out.append((col, f"{len(uniques):,} different values, too many for a category"))
            continue
        base = loupe_header(col) or "Column"
        header, n = base, 2
        while header.casefold() in used:
            header, n = f"{base} {n}", n + 1
        used.add(header.casefold())
        columns.append(LoupeColumn(header, col, kind, [str(u) for u in uniques], np.asarray(codes, np.int32)))
    return columns, left_out


def loupe_texts(columns: list[LoupeColumn], rows: np.ndarray) -> tuple[list[str], np.ndarray]:
    """The CSV text that follows the barcode on the row of a bin, for the per-bin table rows `rows`: (texts, the index
    of each row's text). Rows with the same labels share one text (a section has a few thousand combinations, so a
    file of 400,000 bins holds little besides its barcodes); the last text is that of a bin the run did not analyse."""
    codes = np.zeros((len(rows), len(columns)), np.int32)
    for j, col in enumerate(columns):
        codes[:, j] = col.codes[rows]
    combo, combos = pd.factorize(pd.Series([r.tobytes() for r in codes], dtype=object))
    texts = []
    for combination in combos:
        ks = np.frombuffer(combination, np.int32).tolist()
        texts.append("".join("," + (col.categories[k] if k >= 0 else LOUPE_MISSING) for col, k in zip(columns, ks)))
    texts.append("".join("," + LOUPE_NOT_ANALYZED for _ in columns))
    return texts, np.asarray(combo, np.int64)


def write_loupe_csv(zf: zipfile.ZipFile, name: str, headers: list[str], barcodes: list[str], src: np.ndarray,
                    texts: list[str], combo: np.ndarray, extra: tuple[str, list[str]] | None = None) -> None:
    """A CSV entry of the zip file: the header, then a row per barcode with the labels of the section row src[i] of
    the per-bin table (-1: not analysed; texts and combo from loupe_texts), LOUPE_CHUNK_ROWS rows at a time. extra:
    (header, label per barcode) of a last column that describes the bin itself (the 8 um QC). LF line endings whatever
    the platform, UTF-8 without a byte-order mark."""
    order = np.full(len(src), len(texts) - 1, np.int64)
    ok = src >= 0
    order[ok] = combo[src[ok]]
    order = order.tolist()
    if extra:
        shared = {v: f",{v}" for v in set(extra[1])}  # one text per label, not per row
        tails = [shared[v] for v in extra[1]]
    else:
        tails = [""] * len(barcodes)
    with zf.open(name, "w") as fh:
        fh.write(("barcode" + "".join("," + h for h in headers) + (f",{extra[0]}" if extra else "") + "\n")
                 .encode("utf-8"))
        for start in range(0, len(barcodes), LOUPE_CHUNK_ROWS):
            stop = start + LOUPE_CHUNK_ROWS
            fh.write("".join(f"{b}{texts[t]}{x}\n" for b, t, x in
                             zip(barcodes[start:stop], order[start:stop], tails[start:stop])).encode("utf-8"))


def matrix_barcodes(path: Path) -> list[str]:
    """The barcodes of a feature-barcode matrix as Space Ranger wrote them (the counts are not read)."""
    try:
        with h5py.File(path, "r") as f:
            m = f["matrix"]
            n = int(m["shape"][()][1])
            if n > RAW_SECTION_MAX_BINS:
                raise DataError(f"{path} holds {n:,} bins; the Loupe Browser export reads at most "
                                f"{RAW_SECTION_MAX_BINS:,} bins per section")
            barcodes = decode_all(m["barcodes"][()])
    except (OSError, KeyError, ValueError, IndexError, TypeError) as exc:
        raise DataError(f"{path} is not a readable feature-barcode matrix ({exc})") from exc
    if len(barcodes) != n:
        raise DataError(f"{path}: the barcodes do not match the matrix shape")
    joined = "\0".join(barcodes)
    if any(ch in joined for ch in ',"\r\n'):
        raise DataError(f"{path}: a barcode holds a comma, a quote or a line break, which a category file for Loupe "
                        f"Browser cannot carry")
    return barcodes


def hd_grid_keys(barcodes) -> tuple[np.ndarray, np.ndarray]:
    """(bin size in um, grid key row * GRID_KEY + column) of Visium HD bin barcodes (s_008um_00269_00526, with or
    without the -1 suffix); -1 for a barcode that is not one. The layout Space Ranger writes (five-digit row and
    column) is read as bytes, a few hundred thousand barcodes in a fraction of a second; any other with the regex."""
    values = [str(b) for b in barcodes]
    size = np.full(len(values), -1, np.int64)
    key = np.full(len(values), -1, np.int64)
    rest = np.arange(len(values))
    try:
        b = np.array(values, dtype="S") if values else None
    except UnicodeEncodeError:
        b = None
    if b is not None and b.dtype.itemsize >= 19:
        m = b.view(np.uint8).reshape(len(values), b.dtype.itemsize)  # one byte per character, 0 after the end
        digit = (m >= ord("0")) & (m <= ord("9"))
        fixed = digit[:, [2, 3, 4, *range(8, 13), *range(14, 19)]].all(axis=1)
        for pos, char in ((0, "s"), (1, "_"), (5, "u"), (6, "m"), (7, "_"), (13, "_")):
            fixed &= m[:, pos] == ord(char)
        if m.shape[1] > 19:  # after the column: nothing, or -<digits>
            ends = m[:, 19] == 0
            if m.shape[1] > 20:
                ends |= (m[:, 19] == ord("-")) & digit[:, 20] & (digit[:, 21:] | (m[:, 21:] == 0)).all(axis=1)
            fixed &= ends

        def number(sel: np.ndarray, cols: range) -> np.ndarray:
            return (m[sel][:, cols].astype(np.int64) - ord("0")) @ (10 ** np.arange(len(cols) - 1, -1, -1))

        size[fixed] = number(fixed, range(2, 5))
        key[fixed] = number(fixed, range(8, 13)) * GRID_KEY + number(fixed, range(14, 19))
        rest = np.flatnonzero(~fixed)
    for i in rest.tolist():
        match = HD_BARCODE_RE.match(values[i])
        if match:
            size[i] = int(match[1])
            key[i] = int(match[2]) * GRID_KEY + int(match[3])
    return size, key


def parent_keys(keys: np.ndarray, fine_um: int, coarse_um: int) -> np.ndarray:
    """Grid key of the coarse bin that holds each fine bin, from the grid indices alone: the grids of all bin sizes
    start at the same corner, so fine index i lies in coarse bin i * fine_um // coarse_um (8 into 16 um: i // 2, as
    bone_analysis.map_domains_from_coarser projects labels). -1 stays -1."""
    row, col = np.divmod(keys, GRID_KEY)
    out = (row * fine_um // coarse_um) * GRID_KEY + col * fine_um // coarse_um
    return np.where(keys >= 0, out, -1)


def position_keys(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """(grid key of each row of a tissue positions file, in tissue or not, -1 when its barcode is not a Visium HD bin
    barcode; the full-resolution centres x, y as an n x 2 array, NaN when not a number). The barcodes are dropped once
    parsed: an 8 um positions file has about 700,000 rows."""
    pos = read_positions(path, normalise=False)
    _, keys = hd_grid_keys(pos["barcode"])
    return keys, pos[["x", "y"]].to_numpy(float)


def nearest_parent_keys(fine_positions: Path, coarse_positions: Path, fine_keys: np.ndarray) -> np.ndarray:
    """Grid key of the coarse bin whose centre is nearest to the centre of each fine bin (fine_keys: the grid keys of
    the fine matrix barcodes, -1 for none), i.e. the coarse square the fine centre lies in; every coarse bin of the
    positions file counts, in tissue or not. -1 for a fine bin without a position. Bins are found in the positions
    files by grid key, not by barcode text (a few hundred thousand barcodes would cost far more memory)."""
    keys, xy = position_keys(coarse_positions)
    ok = (keys >= 0) & np.isfinite(xy).all(axis=1)
    out = np.full(len(fine_keys), -1, np.int64)
    if not ok.any():
        return out
    tree, keys = cKDTree(xy[ok]), keys[ok]
    del xy
    pos_keys, pos_xy = position_keys(fine_positions)
    rows = np.flatnonzero(pos_keys >= 0)
    index = pd.Index(pos_keys[rows])
    if not index.is_unique:
        first = ~index.duplicated()
        index, rows = index[first], rows[first]
    at = np.where(fine_keys >= 0, index.get_indexer(fine_keys), -1)
    fxy = np.full((len(fine_keys), 2), np.nan)
    fxy[at >= 0] = pos_xy[rows[at[at >= 0]]]
    have = np.isfinite(fxy).all(axis=1)
    if have.any():
        _, nearest = tree.query(fxy[have])
        out[have] = keys[nearest]
    return out


def matrix_qc(path: Path, qc: dict) -> tuple[np.ndarray, np.ndarray]:
    """(whether each barcode of a feature-barcode matrix passes the quality control `qc` (passes_qc), its UMIs), in
    matrix order, from the Gene Expression counts as RawSection reads them: LOUPE_QC_CHUNK stored counts at a time,
    since an 8 um section holds tens of millions."""
    try:
        with h5py.File(path, "r") as f:
            m = f["matrix"]
            n_features, n = (int(v) for v in m["shape"][()][:2])
            indptr = m["indptr"][()].astype(np.int64)
            features = m["features"]
            names = decode_all(features["name"][()])
            kinds = decode_all(features["feature_type"][()]) if "feature_type" in features else None
            if len(indptr) != n + 1 or len(names) != n_features:
                raise ValueError("the barcodes and features do not match the matrix shape")
            gene = (np.array([k == GENE_EXPRESSION for k in kinds], dtype=bool) if kinds is not None
                    else np.ones(n_features, dtype=bool))
            mito = gene & np.array([str(g).lower().startswith(MT_PREFIX) for g in names], dtype=bool)
            total, n_genes, mt = (np.zeros(n, np.int64) for _ in range(3))
            start = 0
            while start < n:
                stop = int(np.searchsorted(indptr, indptr[start] + LOUPE_QC_CHUNK, side="right")) - 1
                stop = min(max(stop, start + 1), n)
                lo, hi = int(indptr[start]), int(indptr[stop])
                data, cols = m["data"][lo:hi], m["indices"][lo:hi].astype(np.int64)
                row = np.repeat(np.arange(stop - start), np.diff(indptr[start:stop + 1]))
                keep = gene[cols] & (data != 0)
                row, data, cols = row[keep], data[keep].astype(np.float64), cols[keep]
                total[start:stop] = np.rint(np.bincount(row, weights=data, minlength=stop - start))
                n_genes[start:stop] = np.bincount(row, minlength=stop - start)
                is_mt = mito[cols]
                mt[start:stop] = np.rint(np.bincount(row[is_mt], weights=data[is_mt], minlength=stop - start))
                start = stop
    except (OSError, KeyError, ValueError, IndexError, TypeError) as exc:
        raise DataError(f"{path} is not a readable feature-barcode matrix ({exc})") from exc
    return passes_qc(n_genes, total, qc_pct_mt(mt, total), qc), total


def fine_bin_rows(bin_dir: Path, fine_dir: Path, run_um: int, run_raw: list[str], run_src: np.ndarray
                  ) -> tuple[list[str], np.ndarray, dict]:
    """(8 um barcodes as in square_008um/filtered_feature_bc_matrix.h5, for each the section row of the per-bin
    table whose calls it takes or -1, how the parents were found). run_raw, run_src: the barcodes of the run's matrix
    and the section row of each. The parent of an 8 um bin comes from the grid indices (parent_keys), checked against
    the bin centres of the two tissue positions files; when any bin's centre lies outside the square of its parent by
    index, every bin takes the bin whose centre is nearest instead. DataError when the barcodes are not Visium HD bin
    barcodes."""
    raw = matrix_barcodes(fine_dir / "filtered_feature_bc_matrix.h5")
    fine_size, fine_keys = hd_grid_keys(raw)
    run_size, run_keys = hd_grid_keys(run_raw)
    fine_keys[fine_size != LOUPE_BIN_UM] = -1
    run_keys[run_size != run_um] = -1
    if not (fine_keys >= 0).any() or not (run_keys >= 0).any():
        raise DataError("the barcodes are not Visium HD bin barcodes (s_<size>um_<row>_<column>)")
    parents = parent_keys(fine_keys, LOUPE_BIN_UM, run_um)
    check = {"method": "grid", "checked": 0, "inside": 0, "note": None}
    fine_pos, run_pos = positions_file(fine_dir / "spatial"), positions_file(bin_dir / "spatial")
    if fine_pos is None or run_pos is None:
        check["note"] = "not checked against the bin centres: a tissue positions file is missing"
    else:
        try:
            by_position = nearest_parent_keys(fine_pos, run_pos, fine_keys)
        except DataError as exc:
            check["note"] = f"not checked against the bin centres ({exc})"
        else:
            both = (by_position >= 0) & (parents >= 0)
            check["checked"] = int(both.sum())
            check["inside"] = int((both & (by_position == parents)).sum())
            if not check["checked"]:
                check["note"] = "not checked against the bin centres: no bin has a position in both files"
            elif check["inside"] < check["checked"]:
                check["method"] = "position"
                parents = by_position
    valid = np.flatnonzero(run_keys >= 0)
    lookup = pd.Index(run_keys[valid])
    if not lookup.is_unique:
        first = ~lookup.duplicated()
        lookup, valid = lookup[first], valid[first]
    hit = lookup.get_indexer(parents)
    src = np.where(hit >= 0, run_src[valid[np.maximum(hit, 0)]], -1)
    return raw, src, check


def safe_file_part(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", str(text)).strip("-.") or "section"


def remove_quietly(path) -> None:
    with contextlib.suppress(OSError):
        os.unlink(path)


def loupe_export(df: pd.DataFrame, cfg: dict, bin_size: str, out, title: str = "", source: str = "") -> dict:
    """Write the Loupe Browser export of a run into the zip file `out` (a path or a binary file): for each section a
    category CSV at the run's bin size and, when that is a whole multiple of 8 um and the section has
    square_008um/filtered_feature_bc_matrix.h5, one at 8 um (the bin size of Space Ranger's .cloupe file), in which
    every 8 um bin takes the calls of the run bin that holds it and a last column, LOUPE_QC_HEADER, says whether the
    8 um bin itself passes the run's quality control; then README.txt. df: the per-bin table
    (per_bin_signature_calls_{B}.csv.gz, named `source`); cfg: the run's configuration with an absolute data_root.
    Returns what README.txt reports. DataError when a section's data cannot be found or read (the zip is unusable
    then)."""
    try:
        run_um = int(pc.bin_um(bin_size))
    except ValueError as exc:
        raise DataError(str(exc)) from None
    if run_um < LOUPE_BIN_UM:
        raise DataError(f"this run analysed {run_um} um bins; Loupe Browser files (.cloupe) of Visium HD data exist "
                        f"for bins of {LOUPE_BIN_UM} um and larger only")
    for col in ("barcode", "sample_id"):
        if col not in df.columns:
            raise DataError(f"the per-bin table has no {col} column")
    sample_ids = df["sample_id"].map(category_label)
    groups = dict(sample_ids.groupby(sample_ids, sort=False).indices)
    if not groups:
        raise DataError("the per-bin table has no bins")
    entries = {s["id"]: s for s in config_samples(cfg)}
    plan = []
    for sid in [s for s in entries if s in groups] + [s for s in groups if s not in entries]:
        sample = entries.get(sid)
        if sample is None:
            raise DataError(f"section {sid} of the per-bin table is not in the run's configuration, so its Space "
                            f"Ranger outputs cannot be found")
        bin_dir = pc.resolve_bin_dir(cfg, sample, bin_size)
        if bin_dir is None:
            raise DataError(f"the raw data of section {sid} can no longer be found (no square_{bin_size}/"
                            f"filtered_feature_bc_matrix.h5 in {sample_base(cfg, sample)}); the Loupe Browser export "
                            f"lists the barcodes of the Space Ranger outputs the run analysed")
        fine_dir = None
        if run_um > LOUPE_BIN_UM and run_um % LOUPE_BIN_UM == 0:
            fine_dir = pc.resolve_bin_dir(cfg, sample, LOUPE_BIN_SIZE)
        plan.append((sid, sample, bin_dir, fine_dir))

    columns, left_out = loupe_columns(df)
    qc, _ = qc_thresholds(cfg)
    clustering = cfg.get("clustering") if isinstance(cfg.get("clustering"), dict) else {}
    report = {"title": title, "source": source, "bin_size": bin_size, "bin_um": run_um,
              "columns": [{"header": c.header, "source": c.source, "kind": c.kind, "categories": c.categories}
                          for c in columns],
              "left_out": left_out, "sections": [], "qc": qc,
              "min_domain_bins": clustering.get("min_domain_bins")}
    headers = [c.header for c in columns]
    qc_header = LOUPE_QC_HEADER
    while qc_header.casefold() in {h.casefold() for h in headers}:
        qc_header += " 2"
    barcode_values = df["barcode"].to_numpy(object)
    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for sid, sample, bin_dir, fine_dir in plan:
            section = {"id": sid, "label": sample.get("label") or sid, "folder": str(bin_dir.parent), "files": [],
                       "notes": []}
            idx = groups[sid]
            matrix = bin_dir / "filtered_feature_bc_matrix.h5"
            raw = matrix_barcodes(matrix)
            table = ["" if v is None or (not isinstance(v, str) and pd.isna(v)) else str(v)
                     for v in barcode_values[idx].tolist()]
            rows = run_barcode_rows(pd.Index(make_unique(norm_barcodes(raw).tolist())), table, sid)
            matched = np.flatnonzero(rows >= 0)
            if len(idx) and not len(matched):
                raise DataError(f"none of the {len(idx):,} bins of section {sid} in the per-bin table is a barcode of "
                                f"{matrix}; is this the data the run analysed?")
            src = np.full(len(raw), -1, np.int64)  # per matrix barcode: its row among the section's rows of the table
            src[rows[matched]] = matched
            texts, combo = loupe_texts(columns, idx)
            name = f"{safe_file_part(sid)}_{bin_size}_categories.csv"
            write_loupe_csv(zf, name, headers, raw, src, texts, combo)
            analyzed = int((src >= 0).sum())
            section["files"].append({"file": name, "bin_size": bin_size, "bin_um": run_um, "rows": len(raw),
                                     "analyzed": analyzed, "unplaced": int(len(idx) - analyzed)})
            if fine_dir is not None:
                try:
                    fine_raw, fine_src, check = fine_bin_rows(bin_dir, fine_dir, run_um, raw, src)
                except DataError as exc:
                    section["notes"].append(f"no {LOUPE_BIN_UM} um file: {exc}")
                else:
                    # the bin's own quality control: its calls are those of the run bin that holds it, whatever its
                    # own depth, so that a low-quality 8 um bin can be told apart in Loupe Browser
                    entry = {"file": f"{safe_file_part(sid)}_{LOUPE_BIN_SIZE}_categories.csv",
                             "bin_size": LOUPE_BIN_SIZE, "bin_um": LOUPE_BIN_UM, "rows": len(fine_raw),
                             "analyzed": int((fine_src >= 0).sum()), "parents": check, "qc": None}
                    extra = None
                    try:
                        passes, umis = matrix_qc(fine_dir / "filtered_feature_bc_matrix.h5", qc)
                    except DataError as exc:
                        section["notes"].append(f"no {qc_header} column: {exc}")
                    else:
                        labelled = fine_src >= 0
                        entry["qc"] = {"header": qc_header, "pass": int(passes.sum()),
                                       "labelled_fail": int((labelled & ~passes).sum()),
                                       "labelled_no_umi": int((labelled & (umis == 0)).sum()),
                                       "unlabelled_pass": int((~labelled & passes).sum())}
                        extra = (qc_header, ["pass" if ok else "fail" for ok in passes.tolist()])
                    write_loupe_csv(zf, entry["file"], headers, fine_raw, fine_src, texts, combo, extra)
                    section["files"].append(entry)
            elif run_um > LOUPE_BIN_UM:
                section["notes"].append(
                    f"no {LOUPE_BIN_UM} um file: {run_um} um is not a whole multiple of {LOUPE_BIN_UM} um"
                    if run_um % LOUPE_BIN_UM else f"no {LOUPE_BIN_UM} um file: the section has no "
                    f"square_{LOUPE_BIN_SIZE}/filtered_feature_bc_matrix.h5")
            report["sections"].append(section)
        zf.writestr("README.txt", loupe_readme(report).encode("utf-8"))
    return report


def loupe_readme(report: dict) -> str:
    """README.txt of the export: which file goes with which .cloupe file, how to import it, what the columns mean."""
    run_um = report["bin_um"]
    lines = ["BONE SPATIAL STUDIO - CATEGORIES FOR LOUPE BROWSER", ""]
    if report.get("title"):
        lines.append(f"Run: {report['title']}")
    lines += [f"Calls from: {report.get('source') or 'the per-bin table'} ({run_um} um bins)",
              f"Written: {time.strftime('%Y-%m-%d %H:%M')}", "",
              "These CSV files put the calls of the run on the bins of 10x Genomics Loupe Browser, next to Loupe's",
              "own clusters. Loupe imports a file into the .cloupe file of one section at one bin size, and only",
              "into that one: the barcodes must be the bins of that .cloupe file.", "",
              "FILES", "-----"]
    for section in report["sections"]:
        lines += [f"Section {section['id']} ({section['label']})", f"  Space Ranger output: {section['folder']}"]
        for f in section["files"]:
            lines.append(f"  {f['file']}")
            if f["bin_um"] == LOUPE_BIN_UM:
                lines.append(f"    for the section's own Loupe file at {LOUPE_BIN_UM} um, outs/cloupe_{LOUPE_BIN_SIZE}"
                             f".cloupe (Space Ranger writes it for every Visium HD run)")
            else:
                lines += [f"    for a {f['bin_um']} um Loupe file, if you have one: Space Ranger writes it only when "
                          f"run with --custom-bin-size {f['bin_um']}",
                          f"    (spaceranger count 3.1 or later); not for cloupe_{LOUPE_BIN_SIZE}.cloupe"]
            if "parents" in f:
                lines.append(f"    {f['rows']:,} bins (every barcode of square_{f['bin_size']}/"
                             f"filtered_feature_bc_matrix.h5): {f['analyzed']:,} take the calls of their "
                             f"{run_um} um bin, {f['rows'] - f['analyzed']:,} are \"{LOUPE_NOT_ANALYZED}\"")
                q = f.get("qc")
                if q:
                    lines.append(f"    column \"{q['header']}\": {q['pass']:,} bins pass the run's quality control "
                                 f"at {LOUPE_BIN_UM} um; {q['labelled_fail']:,} of the {f['analyzed']:,} with calls "
                                 f"fail it ({q['labelled_no_umi']:,} have no UMI at all), and {q['unlabelled_pass']:,} "
                                 f"\"{LOUPE_NOT_ANALYZED}\" bins pass it")
            else:
                lines.append(f"    {f['rows']:,} bins (every barcode of square_{f['bin_size']}/"
                             f"filtered_feature_bc_matrix.h5): {f['analyzed']:,} analysed, "
                             f"{f['rows'] - f['analyzed']:,} \"{LOUPE_NOT_ANALYZED}\"")
                if f.get("unplaced"):
                    lines.append(f"    {f['unplaced']:,} bins of the per-bin table are not in that matrix and "
                                 f"are left out")
        lines += [f"  Note: {n}" for n in section["notes"]]
        lines.append("")

    fine = [(s, f) for s in report["sections"] for f in s["files"] if "parents" in f]
    qc_col = next((f["qc"]["header"] for _, f in fine if f.get("qc")), None)
    if fine:
        ratio = run_um // LOUPE_BIN_UM
        lines += [f"{LOUPE_BIN_UM} UM FILES", "-" * len(f"{LOUPE_BIN_UM} UM FILES"),
                  f"An {LOUPE_BIN_UM} um bin takes the calls of the {run_um} um bin that holds it. Space Ranger's grids",
                  f"of all bin sizes start at the same corner, so {LOUPE_BIN_UM} um bin (row r, column c) lies in "
                  f"{run_um} um bin",
                  f"(r // {ratio}, c // {ratio}). This was checked against the bin centres of the tissue positions "
                  f"files:"]
        for s, f in fine:
            p = f["parents"]
            if p["note"]:
                lines.append(f"  {s['id']}: {p['note']}")
            elif p["method"] == "grid":
                lines.append(f"  {s['id']}: the centre of each of the {p['checked']:,} {LOUPE_BIN_UM} um bins lies "
                             f"inside its {run_um} um bin (its nearest {run_um} um centre)")
            else:
                lines.append(f"  {s['id']}: {p['checked'] - p['inside']:,} of {p['checked']:,} {LOUPE_BIN_UM} um "
                             f"centres lie outside the {run_um} um bin of that row and column, so every "
                             f"{LOUPE_BIN_UM} um bin took the calls of the {run_um} um bin whose centre is nearest")
        lines += ["", f"The calls are those of the {run_um} um bin; the {LOUPE_BIN_UM} um bin's own counts play no part. "
                  f"So an {LOUPE_BIN_UM} um bin",
                  "that would fail quality control on its own takes calls too, even a bin without any UMI."]
        if qc_col:
            qc = report.get("qc") or {}
            lines += [f"The last column, \"{qc_col}\", says \"pass\" or \"fail\" for each {LOUPE_BIN_UM} um bin "
                      f"by the run's own thresholds",
                      f"(more than {qc.get('min_genes'):g} genes, more than {qc.get('min_counts'):g} UMIs and less "
                      f"than {qc.get('max_pct_mt'):g}% mitochondrial UMIs in the bin): leave out",
                      "the \"fail\" bins before you compare groups in Loupe Browser (a differential expression, for "
                      "example).",
                      f"The run's own projection onto {LOUPE_BIN_UM} um bins (step \"{pc.STEP_LABELS['project_fine_bins']}"
                      f"\") labels only the",
                      "bins that pass."]
        lines.append("")

    lines += ["HOW TO IMPORT", "-------------",
              "1. Open the section's .cloupe file in Loupe Browser 8 or later (9.1 or later for 11 mm slides).",
              "2. In the Clusters panel on the left, click the upload icon next to \"+ Create a new group\".",
              "   Older versions: in Categories mode, open the action menu next to the category name and choose",
              "   \"Import Categories\".",
              "3. Choose the CSV file of that section and bin size. Each column becomes a group (category) of its",
              "   own; choose one to color the bins by its labels.",
              "If Loupe Browser reports \"Empty barcode list supplied\", the file is not the one for this .cloupe file",
              "(another section or bin size), or a program changed it after the export (a spreadsheet program that",
              "saved it again, for example): use the file as it is in this zip.",
              "In Loupe Browser, the \"barcodes\" (and \"cells\") of Visium HD data are bins.", "",
              "COLUMNS", "-------"]
    described = [("barcode", "the bin's barcode in Space Ranger's filtered matrix, e.g. s_008um_00269_00526-1 "
                             "(bin size, row, column)")]
    for c in report["columns"]:
        described.append((c["header"], f"{c['source']}: {loupe_column_text(c, report.get('min_domain_bins'))}"))
    width = max(len(name) for name, _ in described) + 2
    lines += [f"{name.ljust(width)}{text}" for name, text in described]
    pad = " " * 18
    if qc_col:
        lines.append(f"{qc_col.ljust(width)}{LOUPE_BIN_SIZE} files only: \"pass\" or \"fail\", the "
                     f"{LOUPE_BIN_UM} um bin's own quality control (see {LOUPE_BIN_UM} UM FILES)")
    lines += ["", f"In every column{f' but {qc_col}' if qc_col else ''}:",
              f"  \"{LOUPE_NOT_ANALYZED}\"  the bin is in Space Ranger's filtered matrix (so in the .cloupe file), but "
              f"the run has no calls",
              f"{pad}for it: in a {run_um} um file, the bin did not pass quality control or has no position"
              + (";" if fine else ".")]
    if fine:
        lines.append(f"{pad}in an {LOUPE_BIN_UM} um file, the {run_um} um bin holding it was not analysed (the "
                     f"{LOUPE_BIN_UM} um bin itself may pass).")
    lines.append(f"  \"{LOUPE_MISSING}\"{' ' * 12}the bin was analysed but has no value in this column.")
    if report["left_out"]:
        lines += ["Not exported:"] + [f"  {col}: {why}" for col, why in report["left_out"]]
    lines += ["The columns follow the run's configuration: a run with other signatures, classes or steps has other",
              "columns, read the same way.", "",
              "FORMAT", "------",
              "Comma-separated text, UTF-8 without a byte-order mark, LF line endings. The first row holds the",
              "column names; the first column is \"barcode\"; one row per bin; every cell holds a label (none is",
              "empty); labels hold no commas, quotes or line breaks (a comma in a name became \";\").",
              "This follows 10x Genomics' description of category CSV files. Where it says nothing, the most",
              "compatible choice was made: the header \"barcode\" in lower case as documented (Loupe's own exports",
              "write \"Barcode\"), the -1 suffix kept on every barcode, a label in every cell. The files were not",
              "tried in Loupe Browser itself.", ""]
    return "\n".join(lines)


def loupe_column_text(column: dict, min_domain_bins=None) -> str:
    """What a category column of the export holds (README.txt); min_domain_bins: the run's
    clustering.min_domain_bins (domains smaller than it are "unassigned")."""
    source, cats = column["source"], column["categories"]
    small = (f"domains under {min_domain_bins:,} bins" if isinstance(min_domain_bins, int)
             and not isinstance(min_domain_bins, bool) else "domains smaller than clustering.min_domain_bins")
    known = {
        "tissue_class": "the drawn tissue class of the bin's assigned signature (\"removed\": a signature whose class "
                        "is not drawn; \"unassigned\": no single best signature)",
        "assigned_signature": "the signature the bin was assigned to (\"ambiguous\": no single best signature)",
        "leiden_domain": "the bin's Leiden domain (\"domain 0\", \"domain 1\", ...; \"unresolved\": domains of "
                         f"ubiquitous transcripts merged; \"unassigned\": {small})",
        "in_region_of_interest": "\"inside\" or \"outside\" the region of interest",
    }
    if source in known:
        return known[source]
    if column["kind"] == "bool":
        yes, no = yes_no_labels(source)
        what = ""
        if source.endswith("_pooled_depthmatched"):
            what = "positivity call after spatial pooling and depth matching: "
        elif source.endswith("_pooled"):
            what = "positivity call after spatial pooling: "
        elif source.startswith("positive_"):
            what = "positivity call: "
        return f"{what}\"{yes}\" or \"{no}\""
    shown = ", ".join(f"\"{c}\"" for c in cats[:6])
    return f"{len(cats)} labels ({shown}{', ...' if len(cats) > 6 else ''})"


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
    # RawSection by matrix + positions file (path, mtime, size), bounded by their estimated memory
    sections = LRUCache(RAW_SECTIONS_CACHED, RAW_SECTIONS_BUDGET, lambda s: s.nbytes)
    explore_cache = LRUCache(EXPLORE_PAYLOADS_CACHED, EXPLORE_PAYLOADS_BUDGET, len)  # JSON of POST /explore/../bins
    gene_links = LRUCache(8)                   # RunGeneLink by run, sample, per-bin table and section
    image_keys = ImageKeys()
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
        gene_links.discard(lambda key: key[0] == rid)

    def raw_section(matrix: Path, positions: Path) -> RawSection:
        return sections.get(section_key(matrix, positions), lambda: RawSection(matrix, positions))

    def explore_config(pid: str, body: dict) -> dict:
        """Configuration of the data explorer: the unsaved edits sent with the request (body.config, as for POST
        /validate) or else the saved one, merged with the defaults, data_root made absolute."""
        project = storage.get_project(pid)
        given = optional_object(body, "config")
        cfg = anchored(pc.merge_with_defaults(given if given is not None else project.get("config") or {}))
        try:
            pc.bin_um(cfg.get("bin_size"))
        except ValueError:
            raise HTTPException(400, f"bin_size must look like '016um', got {cfg.get('bin_size')!r}") from None
        return cfg

    def explore_sample(pid: str, sid: str, body: dict) -> tuple[dict, dict]:
        cfg = explore_config(pid, body)
        for s in config_samples(cfg):
            if s["id"] == sid:
                return cfg, s
        raise NotFound(f"sample {sid!r} is not in the configuration of project {pid}")

    def run_sample(rid: str, sid: str) -> tuple[dict, dict]:
        """(configuration of the run with an absolute data_root, sample entry)."""
        require_run(rid)
        cfg = anchored(run_config(rid))
        for s in config_samples(cfg):
            if s["id"] == sid:
                return cfg, s
        raise NotFound(f"sample {sid!r} is not in run {rid}")

    def run_gene_link(rid: str, sid: str) -> tuple[RunGeneLink, RawSection]:
        """The bins of a sample of a run matched to its raw section, read again from the run's data_root."""
        table = bins_table(rid)
        if sid not in table.groups:
            raise NotFound(f"sample {sid!r} is not in the per-bin table of run {rid}")
        cfg, sample = run_sample(rid, sid)
        try:
            matrix, positions = section_files(cfg, sample)
        except NotFound as exc:
            raise DataError(f"the raw data of run {rid} can no longer be found ({exc}); gene maps need the Space "
                            f"Ranger outputs the run analysed") from None
        section = raw_section(matrix, positions)
        link = gene_links.get((rid, sid, *table.source, *section.key),
                              lambda: RunGeneLink(section, table.barcodes(sid), sid))
        return link, section

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
        data = body.get("data", True)
        if not isinstance(data, bool):
            raise HTTPException(400, "data must be true or false")
        if not data:  # the configuration check alone (the project page checks the unsaved settings as they change)
            return JSONOut({"ok": not errors, "errors": errors, "warnings": warnings})
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

    # ---- data explorer (raw sections; unsaved settings like POST /validate) ----------------------------------------
    @app.post("/api/projects/{pid}/explore/samples")
    def explore_samples(pid: str, body: dict = Depends(optional_json_body)):
        cfg = explore_config(pid, body)
        out = []
        for sample in config_samples(cfg):
            entry = {"id": sample["id"], "label": sample.get("label") or sample["id"],
                     "condition": sample.get("condition"), "color": sample.get("color"), "n_bins": None,
                     "has_image": hires_image(cfg, sample) is not None, "ok": False, "message": ""}
            try:
                matrix, positions = section_files(cfg, sample)
                entry["n_bins"] = len(in_tissue_xy(positions)[0])  # positions only: the matrix is not read
                if entry["n_bins"] > RAW_SECTION_MAX_BINS:
                    raise too_many_bins(matrix, entry["n_bins"])
            except (NotFound, DataError) as exc:
                entry["message"] = str(exc)
            else:
                entry["ok"] = True
                entry["message"] = f"{entry['n_bins']:,} bins at {cfg['bin_size']}"
            out.append(entry)
        return JSONOut({"bin_size": cfg["bin_size"], "samples": out})

    @app.post("/api/projects/{pid}/explore/{sid}/bins")
    def explore_bins(pid: str, sid: str, body: dict = Depends(optional_json_body)):
        cfg, sample = explore_sample(pid, sid, body)
        matrix, positions = section_files(cfg, sample)
        mpp = section_mpp(sample, matrix.parent / "spatial")
        image, image_note = section_image(cfg, sample, mpp, image_keys)  # issues the image key on every request
        qc, notes = qc_thresholds(cfg)
        notes += [image_note] if image_note else []
        key = (sid, cfg["bin_size"], *section_key(matrix, positions), mpp, tuple(qc.items()), json.dumps(image),
               tuple(notes))

        def build() -> bytes:
            return explore_payload(raw_section(matrix, positions), sid, cfg["bin_size"], mpp, qc, image, notes)

        return Response(explore_cache.get(key, build), media_type=JSON_MEDIA_TYPE)

    @app.post("/api/projects/{pid}/explore/{sid}/genes")
    def explore_genes(pid: str, sid: str, body: dict = Depends(optional_json_body)):
        query, limit = gene_query(body.get("query"), body.get("limit"))
        cfg, sample = explore_sample(pid, sid, body)
        section = raw_section(*section_files(cfg, sample))
        return JSONOut(gene_table(section, section.gene_total, section.gene_bins, section.n_bins, query, limit))

    @app.post("/api/projects/{pid}/explore/{sid}/gene")
    def explore_gene(pid: str, sid: str, body: dict = Depends(json_body)):
        names = requested_genes(body.get("genes"))
        cfg, sample = explore_sample(pid, sid, body)
        section = raw_section(*section_files(cfg, sample))
        used, missing, cols = section.resolve(names)
        return JSONOut(sparse_vector(section.counts(cols), used, missing))

    @app.get("/api/explore-images/{key}")
    def explore_image(key: str):
        path = image_keys.path(key) if IMAGE_KEY_RE.match(key) else None
        if path is None or not path.is_file():
            raise NotFound("unknown image; load the map again")
        return FileResponse(path, media_type="image/png", headers={"Cache-Control": "no-cache"})

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
            return BinsTable(df, cfg, source=file_key(path))

        return bins_cache.get((rid, *file_key(path)), build)

    @app.get("/api/runs/{rid}/bins/samples")
    def bin_samples(rid: str):
        return JSONOut(bins_table(rid).samples())

    @app.get("/api/runs/{rid}/bins")
    def bins(rid: str, sample: str):
        return Response(bins_table(rid).payload(sample), media_type=JSON_MEDIA_TYPE)

    @app.get("/api/runs/{rid}/bins/image")
    def run_bins_image(rid: str, sample: str):
        cfg, entry = run_sample(rid, sample)
        image = hires_image(cfg, entry)
        if image is None:
            raise NotFound(f"sample {sample} of run {rid} has no H&E image (spatial/tissue_hires_image.png)")
        bin_dir = pc.resolve_bin_dir(cfg, entry)
        try:  # the microns per pixel of the per-bin table's x_um, y_um
            mpp = section_mpp(entry, bin_dir / "spatial" if bin_dir is not None else image.parent)
        except DataError as exc:
            raise NotFound(str(exc)) from None
        info, note = section_image(cfg, entry, mpp, image_keys)
        if info is None:
            raise NotFound(f"sample {sample} of run {rid}: {note}")
        return JSONOut(info)

    @app.get("/api/runs/{rid}/genes")
    def run_genes(rid: str, sample: str, query: str = "", limit: int = GENE_LIST_LIMIT):
        query, limit = gene_query(query, limit)
        link, section = run_gene_link(rid, sample)
        total, detected = link.gene_stats(section)
        return JSONOut({**gene_table(section, total, detected, link.n, query, limit),
                        "n_unmatched": link.n_unmatched})

    @app.post("/api/runs/{rid}/gene")
    def run_gene(rid: str, body: dict = Depends(json_body)):
        sid = body.get("sample")
        if not isinstance(sid, str) or not sid:
            raise HTTPException(400, "sample must be the id of a sample of the run")
        names = requested_genes(body.get("genes"))
        link, section = run_gene_link(rid, sid)
        used, missing, cols = section.resolve(names)
        return JSONOut({**sparse_vector(link.counts(section, cols), used, missing), "n_unmatched": link.n_unmatched})

    @app.get("/api/runs/{rid}/export/loupe")
    def export_loupe(rid: str):
        require_run(rid)
        record = jobs.get_run(rid)
        status = record.get("status")
        if jobs.is_busy(rid) or status in ACTIVE_STATUSES:
            raise HTTPException(409, f"run {rid} is {status}; export it for Loupe Browser when it has finished")
        cfg = run_config(rid)
        try:
            path = per_bin_file(storage.results_dir(rid), cfg)
        except NotFound:
            raise NotFound(f"run {rid} has no per-bin table (per_bin_signature_calls_*.csv.gz)"
                           + ("" if status == "succeeded" else f": it stopped ({status}) before writing one")) from None
        m = PER_BIN_RE.match(path.name)
        bin_size = f"{int(m.group(1)):03d}um" if m else cfg.get("bin_size")
        df = bins_table(rid).df
        title = " - ".join(str(v) for v in dict.fromkeys([rid, record.get("label"), record.get("project_name")]) if v)
        # built in a temporary file, so that an error is still a JSON answer and only one CSV is in memory at a time
        fd, tmp = tempfile.mkstemp(prefix="bss-loupe-", suffix=".zip")
        try:
            with os.fdopen(fd, "wb") as fh:
                loupe_export(df, anchored(cfg), bin_size, fh, title=title, source=path.name)
        except BaseException:
            remove_quietly(tmp)
            raise
        return FileResponse(tmp, media_type="application/zip", filename=f"{rid}_loupe.zip",
                            headers={"Cache-Control": "no-store"}, background=BackgroundTask(remove_quietly, tmp))

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
