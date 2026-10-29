"""
webapp/storage.py
=================
Projects and runs of Bone Spatial Studio on disk (docs/ARCHITECTURE.md, section 4):

    <workspace>/projects/<pid>/project.json      {"id", "name", "created", "updated", "config"}
    <workspace>/projects/<pid>/rois/<sid>.geojson
    <workspace>/runs/<rid>/run.json              run record (section 5) without the configuration
    <workspace>/runs/<rid>/config.json           effective configuration snapshot
    <workspace>/runs/<rid>/rois/<sid>.geojson    copies of the drawn regions at submission time
    <workspace>/runs/<rid>/log.txt               worker output (webapp.jobs)
    <workspace>/runs/<rid>/results/              analysis outputs
    <workspace>/demo_data/                       synthetic data (tools/make_demo_data.py)

Every JSON file is written atomically (temporary file + os.replace) while holding a lock shared by all Storage objects
of one workspace, and read under the same lock, so this process never sees a half-written file (and, on Windows, never
holds a file open while it is being replaced).
"""
from __future__ import annotations

import copy
import json
import logging
import math
import os
import re
import secrets
import shutil
import stat
import sys
import threading
import time
from datetime import date, datetime, timezone
from pathlib import Path

from pipeline_config import merge_with_defaults

log = logging.getLogger("bss.storage")

PROJECT_ID_RE = re.compile(r"^p-[0-9a-f]{8}\Z")
RUN_ID_RE = re.compile(r"^r-\d{8}-\d{6}-[0-9a-f]{4}\Z")
# The rule pipeline_config.validate_config applies to sample ids; it also keeps rois/<sid>.geojson inside rois/.
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]{0,63}\Z")
_WINDOWS_DEVICE_NAMES = {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
                         *(f"LPT{i}" for i in range(1, 10))}

FINAL_STATUSES = frozenset({"succeeded", "failed", "cancelled", "interrupted"})
ACTIVE_STATUSES = frozenset({"queued", "running"})
_EPOCH = datetime(1970, 1, 1, tzinfo=timezone.utc)


class StorageError(Exception):
    """A workspace file cannot be read or written."""


class NotFound(StorageError):
    """No such project, run or region of interest."""


class BadRequest(StorageError):
    """An id or a value that cannot be stored."""


# ============================================================================
# HELPERS
# ============================================================================

def now_iso() -> str:
    """Local time with its UTC offset, e.g. 2026-10-28T15:30:00.123-05:00 (JavaScript's Date parses it)."""
    return datetime.now().astimezone().isoformat(timespec="milliseconds")


def parse_time(text) -> datetime | None:
    """An ISO timestamp as an aware datetime (naive ones are taken as local time), or None."""
    try:
        dt = datetime.fromisoformat(str(text))
    except (TypeError, ValueError):
        return None
    return dt if dt.tzinfo else dt.astimezone()


def created_key(record: dict) -> tuple:
    """Sort key of a project or run record by creation time (then id)."""
    return parse_time(record.get("created")) or _EPOCH, str(record.get("id", ""))


def _newest_first(records: list[dict]) -> list[dict]:
    return sorted(records, key=created_key, reverse=True)


def json_safe(obj):
    """`obj` as strict JSON data: NaN and infinities become None, numpy scalars and arrays plain Python values, paths
    and datetimes strings, mapping keys strings."""
    if obj is None or isinstance(obj, (bool, int, str)):
        return obj
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {str(k): json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set, frozenset)):
        return [json_safe(v) for v in obj]
    if isinstance(obj, (datetime, date)):
        return obj.isoformat()
    if isinstance(obj, Path):
        return str(obj)
    tolist = getattr(obj, "tolist", None)  # numpy scalars and arrays, pandas objects
    if callable(tolist):
        return json_safe(tolist())
    return str(obj)


def check_sample_id(sid) -> str:
    """`sid` when it can name a file of rois/; BadRequest otherwise."""
    if not isinstance(sid, str) or not SAMPLE_ID_RE.match(sid):
        raise BadRequest(f"invalid sample id {sid!r}: use 1-64 letters, digits, '.', '_' or '-'")
    if os.name == "nt" and sid.split(".")[0].upper() in _WINDOWS_DEVICE_NAMES:
        raise BadRequest(f"sample id {sid!r} is a reserved device name on Windows; rename the sample")
    return sid


def _replace(src: Path, dst: Path, attempts: int = 20) -> None:
    for i in range(attempts):
        try:
            os.replace(src, dst)
            return
        except PermissionError:
            # Windows refuses to replace a file another process (virus scanner, indexer, editor) holds open
            if i == attempts - 1:
                raise
            time.sleep(0.05)


def rmtree(path: Path, attempts: int = 5) -> None:
    """Delete a folder; read-only files are made writable and files locked for a moment are retried (Windows)."""
    def on_error(func, target, _info):
        try:
            os.chmod(target, stat.S_IWRITE)
            func(target)
        except FileNotFoundError:
            pass

    for i in range(attempts):
        try:
            if sys.version_info >= (3, 12):
                shutil.rmtree(path, onexc=on_error)
            else:
                shutil.rmtree(path, onerror=on_error)
            return
        except FileNotFoundError:
            return
        except OSError:
            if i == attempts - 1:
                raise
            time.sleep(0.2)


_LOCKS: dict[str, threading.RLock] = {}
_LOCKS_GUARD = threading.Lock()


def _workspace_lock(root: Path) -> threading.RLock:
    key = os.path.normcase(str(root))
    with _LOCKS_GUARD:
        return _LOCKS.setdefault(key, threading.RLock())


# ============================================================================
# STORAGE
# ============================================================================

class Storage:
    """Thread-safe access to one workspace folder."""

    def __init__(self, workspace):
        self.root = Path(workspace).expanduser().resolve()
        self.projects_dir = self.root / "projects"
        self.runs_dir = self.root / "runs"
        self.demo_dir = self.root / "demo_data"
        self.lock = _workspace_lock(self.root)
        for d in (self.projects_dir, self.runs_dir):
            d.mkdir(parents=True, exist_ok=True)

    # ---- JSON ------------------------------------------------------------------------------------------------------
    def read_json(self, path: Path, what: str | None = None):
        path = Path(path)
        try:
            with self.lock:
                text = path.read_text(encoding="utf-8")
        except FileNotFoundError:
            raise NotFound(f"{what or path.name} not found") from None
        except OSError as exc:
            raise StorageError(f"cannot read {path}: {exc}") from exc
        try:
            return json.loads(text)
        except ValueError as exc:
            raise StorageError(f"{path} is not valid JSON ({exc})") from exc

    def write_json(self, path: Path, data) -> None:
        path = Path(path)
        text = json.dumps(json_safe(data), indent=2, ensure_ascii=False, allow_nan=False)
        tmp = path.with_name(f".{path.name}.{os.getpid()}-{threading.get_ident()}.tmp")
        with self.lock:
            try:
                path.parent.mkdir(parents=True, exist_ok=True)
                with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
                    fh.write(text)
                    fh.flush()
                    os.fsync(fh.fileno())
                _replace(tmp, path)
            except OSError as exc:
                try:
                    tmp.unlink(missing_ok=True)
                except OSError:
                    pass
                raise StorageError(f"cannot write {path}: {exc}") from exc

    # ---- projects --------------------------------------------------------------------------------------------------
    def project_dir(self, pid) -> Path:
        if not isinstance(pid, str) or not PROJECT_ID_RE.match(pid):
            raise NotFound(f"project {pid!r} not found")
        return self.projects_dir / pid

    def get_project(self, pid) -> dict:
        return self.read_json(self.project_dir(pid) / "project.json", f"project {pid}")

    def list_projects(self) -> list[dict]:
        out = []
        for d in self.projects_dir.iterdir():
            path = d / "project.json"
            if not PROJECT_ID_RE.match(d.name) or not path.is_file():
                continue
            try:
                out.append(self.read_json(path))
            except StorageError as exc:
                log.warning("skipping project %s: %s", d.name, exc)
        return _newest_first(out)

    def _new_project(self, name: str, config: dict) -> dict:
        with self.lock:
            while True:
                pid = f"p-{secrets.token_hex(4)}"
                try:
                    (self.projects_dir / pid).mkdir()
                    break
                except FileExistsError:
                    continue
            now = now_iso()
            project = {"id": pid, "name": name, "created": now, "updated": now, "config": config}
            self.write_json(self.projects_dir / pid / "project.json", project)
        return project

    def create_project(self, name: str, config: dict | None = None) -> dict:
        """A new project; `config` is merged with the defaults and named after the project unless it has a name."""
        cfg = merge_with_defaults(config or {})
        if not (config or {}).get("name"):
            cfg["name"] = name
        return self._new_project(name, cfg)

    def update_project(self, pid, name: str | None = None, config: dict | None = None) -> dict:
        with self.lock:
            project = self.get_project(pid)
            if name is not None:
                project["name"] = name
            if config is not None:
                project["config"] = merge_with_defaults(config)
            project["updated"] = now_iso()
            self.write_json(self.project_dir(pid) / "project.json", project)
        return project

    def delete_project(self, pid) -> None:
        d = self.project_dir(pid)
        with self.lock:
            if not (d / "project.json").is_file():
                raise NotFound(f"project {pid} not found")
            rmtree(d)

    def duplicate_project(self, pid) -> dict:
        with self.lock:
            src = self.get_project(pid)
            new = self._new_project(f"{src['name']} (copy)", copy.deepcopy(src["config"]))
            rois = self.project_dir(pid) / "rois"
            if rois.is_dir():
                shutil.copytree(rois, self.project_dir(new["id"]) / "rois")
        return new

    # ---- regions of interest ---------------------------------------------------------------------------------------
    def roi_path(self, pid, sid) -> Path:
        return self.project_dir(pid) / "rois" / f"{check_sample_id(sid)}.geojson"

    def list_rois(self, pid) -> list[str]:
        d = self.project_dir(pid) / "rois"
        if not d.is_dir():
            return []
        return sorted(p.stem for p in d.glob("*.geojson") if SAMPLE_ID_RE.match(p.stem))

    def read_roi(self, pid, sid) -> dict:
        self.get_project(pid)
        return self.read_json(self.roi_path(pid, sid), f"region of interest of sample {sid}")

    def _roi_names_ignoring_case(self, pid, sid) -> list[str]:
        """File names in rois/ that a case-insensitive file system (Windows, macOS) treats as <sid>.geojson."""
        d = self.project_dir(pid) / "rois"
        want = f"{check_sample_id(sid)}.geojson".casefold()
        return [n for n in os.listdir(d) if n.casefold() == want] if d.is_dir() else []

    def write_roi(self, pid, sid, feature_collection: dict) -> None:
        with self.lock:
            self.get_project(pid)
            # a case-only rename (ctrl-1 -> CTRL-1): drop the old spelling first, or the file keeps its old name
            for name in self._roi_names_ignoring_case(pid, sid):
                if name != f"{sid}.geojson":
                    (self.project_dir(pid) / "rois" / name).unlink(missing_ok=True)
            self.write_json(self.roi_path(pid, sid), feature_collection)

    def delete_roi(self, pid, sid) -> None:
        with self.lock:
            self.get_project(pid)
            # only the exact spelling: on a case-insensitive file system "ctrl-1" would also delete "CTRL-1"
            if f"{check_sample_id(sid)}.geojson" not in self._roi_names_ignoring_case(pid, sid):
                return
            try:
                self.roi_path(pid, sid).unlink(missing_ok=True)
            except OSError as exc:
                raise StorageError(f"cannot delete the region of interest of {sid}: {exc}") from exc

    # ---- runs ------------------------------------------------------------------------------------------------------
    def run_dir(self, rid) -> Path:
        if not isinstance(rid, str) or not RUN_ID_RE.match(rid):
            raise NotFound(f"run {rid!r} not found")
        return self.runs_dir / rid

    def results_dir(self, rid) -> Path:
        return self.run_dir(rid) / "results"

    def log_path(self, rid) -> Path:
        return self.run_dir(rid) / "log.txt"

    def get_run(self, rid) -> dict:
        return self.read_json(self.run_dir(rid) / "run.json", f"run {rid}")

    def save_run(self, record: dict) -> None:
        data = {k: v for k, v in record.items() if k != "config"}
        self.write_json(self.run_dir(record["id"]) / "run.json", data)

    def read_run_config(self, rid) -> dict:
        return self.read_json(self.run_dir(rid) / "config.json", f"configuration of run {rid}")

    def list_runs(self, project_id: str | None = None) -> list[dict]:
        out = []
        for d in self.runs_dir.iterdir():
            path = d / "run.json"
            if not RUN_ID_RE.match(d.name) or not path.is_file():
                continue
            try:
                record = self.read_json(path)
            except StorageError as exc:
                log.warning("skipping run %s: %s", d.name, exc)
                continue
            if project_id is None or record.get("project_id") == project_id:
                out.append(record)
        return _newest_first(out)

    def create_run(self, *, project_id: str, project_name: str, label: str, config: dict,
                   roi_files: dict[str, Path], extra: dict | None = None) -> dict:
        """A queued run: the configuration snapshot with roi.masks pointing at copies of `roi_files` (sample id ->
        GeoJSON file) inside the run folder, and its run.json."""
        with self.lock:
            while True:
                rid = datetime.now().strftime("r-%Y%m%d-%H%M%S-") + secrets.token_hex(2)
                try:
                    (self.runs_dir / rid).mkdir()
                    break
                except FileExistsError:
                    continue
        run_dir = self.runs_dir / rid
        try:
            masks = {}
            for sid, src in roi_files.items():
                dst = run_dir / "rois" / f"{check_sample_id(sid)}.geojson"
                dst.parent.mkdir(exist_ok=True)
                shutil.copyfile(src, dst)
                masks[sid] = str(dst)
            cfg = copy.deepcopy(config)
            if not isinstance(cfg.get("roi"), dict):
                cfg["roi"] = {}
            cfg["roi"]["masks"] = masks
            self.write_json(run_dir / "config.json", cfg)
            record = {
                "id": rid, "project_id": project_id, "project_name": project_name, "label": label,
                "status": "queued", "created": now_iso(), "started": None, "finished": None, "duration_s": None,
                "pid": None, "progress": {"index": 0, "total": 0, "label": "", "fraction": 0.0},
                "steps": [], "warnings": [], "error": None, "traceback": None, **(extra or {}),
            }
            self.save_run(record)
        except (OSError, StorageError) as exc:
            rmtree(run_dir)
            if isinstance(exc, StorageError):
                raise
            raise StorageError(f"cannot create the run folder: {exc}") from exc
        return record

    def delete_run(self, rid) -> None:
        d = self.run_dir(rid)
        with self.lock:
            if not (d / "run.json").is_file():
                raise NotFound(f"run {rid} not found")
            rmtree(d)
