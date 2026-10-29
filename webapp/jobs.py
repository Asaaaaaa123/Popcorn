"""
webapp/jobs.py
==============
Job runner of Bone Spatial Studio (docs/ARCHITECTURE.md, section 6).

One analysis at a time (they need a lot of memory), first in first out. A daemon thread takes the next queued run and
starts `python -u -m webapp.worker <run_dir>` in the project root with stdout and stderr merged. Every output line is
appended to the run's log.txt, except the `@@EVENT {json}` lines of pipeline_config.Progress, which update the run
record: in memory at once (what the API serves while the run is active) and in run.json at most every WRITE_INTERVAL
seconds. Cancelling kills the worker's whole process tree. When the server starts, runs left "running" by a previous
server become "interrupted" and runs left "queued" are queued again.
"""
from __future__ import annotations

import collections
import copy
import json
import logging
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

from pipeline_config import EVENT_PREFIX
from webapp.storage import NotFound, Storage, StorageError, created_key, now_iso

log = logging.getLogger("bss.jobs")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
THREAD_VARS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
               "NUMBA_NUM_THREADS")
WRITE_INTERVAL = 0.5      # seconds between two run.json writes of the active run
ERROR_TAIL = 20           # log lines reported as the error when the worker dies without an error event
MAX_WARNINGS = 1000       # warnings kept in the run record (a pathological loop must not bloat run.json)
_WINDOWS_FLAGS = (getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0) | getattr(subprocess, "CREATE_NO_WINDOW", 0))


def user_cache_dir() -> Path:
    """Per-user cache folder of the application (platform convention)."""
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
    elif sys.platform == "darwin":
        base = str(Path.home() / "Library" / "Caches")
    else:
        base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "BoneSpatialStudio"


def worker_env(deterministic: bool) -> dict:
    """Environment of a worker process."""
    env = dict(os.environ)
    # UTF-8 output whatever the console code page, so labels with non-ASCII characters cannot crash the worker
    env.update(PYTHONUNBUFFERED="1", MPLBACKEND="Agg", PYTHONIOENCODING="utf-8")
    # importable even when PYTHONSAFEPATH keeps `-m` from adding the working directory to sys.path
    env["PYTHONPATH"] = os.pathsep.join([str(PROJECT_ROOT)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    # numba (umap, pynndescent) caches compiled functions next to their source by default; in a deeply nested
    # environment those file names exceed the 260-character path limit of Windows and the analysis fails. A short
    # per-user folder also keeps site-packages untouched.
    env.setdefault("NUMBA_CACHE_DIR", str(user_cache_dir() / "numba"))
    if deterministic:
        env.update(PYTHONHASHSEED="0", VISIUM_DETERMINISTIC="1", **{v: "1" for v in THREAD_VARS})
    else:
        env["VISIUM_DETERMINISTIC"] = "0"
    return env


_JOB_HANDLE = None


def kill_workers_with_server() -> None:
    """Make workers die with the server, so a closed or crashed server does not leave an analysis running unseen (the
    run would be marked interrupted on the next start while its process kept computing).
    Windows: the server joins a Job Object with JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE; every process it starts later
    inherits the job, and the system kills them all when the server's handle closes. Linux: see _parent_death_signal.
    Elsewhere workers stop at their next output line (broken pipe)."""
    global _JOB_HANDLE
    if os.name != "nt" or _JOB_HANDLE is not None:
        return
    import ctypes
    from ctypes import wintypes

    class _Basic(ctypes.Structure):
        _fields_ = [("PerProcessUserTimeLimit", ctypes.c_int64), ("PerJobUserTimeLimit", ctypes.c_int64),
                    ("LimitFlags", wintypes.DWORD), ("MinimumWorkingSetSize", ctypes.c_size_t),
                    ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", wintypes.DWORD),
                    ("Affinity", ctypes.c_size_t), ("PriorityClass", wintypes.DWORD), ("SchedulingClass", wintypes.DWORD)]

    class _Extended(ctypes.Structure):
        _fields_ = [("BasicLimitInformation", _Basic), ("IoInfo", ctypes.c_uint64 * 6),
                    ("ProcessMemoryLimit", ctypes.c_size_t), ("JobMemoryLimit", ctypes.c_size_t),
                    ("PeakProcessMemoryUsed", ctypes.c_size_t), ("PeakJobMemoryUsed", ctypes.c_size_t)]

    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.CreateJobObjectW.restype = wintypes.HANDLE
    k32.GetCurrentProcess.restype = wintypes.HANDLE
    job = k32.CreateJobObjectW(None, None)
    info = _Extended()
    info.BasicLimitInformation.LimitFlags = 0x2000                               # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
    if job and k32.SetInformationJobObject(wintypes.HANDLE(job), 9, ctypes.byref(info), ctypes.sizeof(info)) \
            and k32.AssignProcessToJobObject(wintypes.HANDLE(job), wintypes.HANDLE(k32.GetCurrentProcess())):
        _JOB_HANDLE = job                                                        # kept open for the server's lifetime
    else:
        log.warning("workers will not be stopped automatically if the server dies (Job Object error %s)",
                    ctypes.get_last_error())


def _parent_death_signal() -> None:
    """Linux: the worker receives SIGTERM when the server dies (prctl PR_SET_PDEATHSIG); runs in the child."""
    import ctypes
    try:
        ctypes.CDLL("libc.so.6", use_errno=True).prctl(1, signal.SIGTERM)
    except OSError:
        pass


def kill_tree(proc: subprocess.Popen) -> None:
    """Terminate a worker and every process it started."""
    if os.name == "nt":
        if proc.poll() is None:
            try:
                subprocess.run(["taskkill", "/T", "/F", "/PID", str(proc.pid)], stdin=subprocess.DEVNULL,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30,
                               creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            except (OSError, subprocess.SubprocessError) as exc:
                log.warning("taskkill failed for pid %s: %s", proc.pid, exc)
        if proc.poll() is None:
            try:
                proc.kill()
            except OSError:
                pass
        return
    # POSIX: the worker leads its own session (start_new_session=True), so its process group id is its pid
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except (ProcessLookupError, PermissionError):
        return
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass


def _new_step(label: str) -> dict:
    return {"label": label, "status": "pending", "seconds": None, "reason": None}


class _Active:
    """State of the run whose worker is (about to be) running."""

    def __init__(self, rid: str, record: dict):
        self.rid = rid
        self.record = record
        self.proc: subprocess.Popen | None = None
        self.cancel = False            # cancel requested by the user
        self.interrupt = False         # the server is stopping
        self.done = False              # a done event arrived
        self.error: str | None = None  # message of an error event
        self.tail: collections.deque[str] = collections.deque(maxlen=ERROR_TAIL)
        self.dirty = False
        self.t0 = time.monotonic()
        self.finished = threading.Event()


class JobManager:
    """FIFO queue of runs executed one at a time in worker subprocesses."""

    def __init__(self, workspace, worker_module: str = "webapp.worker", storage: Storage | None = None):
        self.storage = storage or Storage(workspace)
        self.worker_module = worker_module
        self._lock = threading.RLock()
        self._wake = threading.Condition(self._lock)
        self._queue: collections.deque[str] = collections.deque()
        self._active: _Active | None = None
        self._stopping = False
        self._thread: threading.Thread | None = None

    # ---- lifecycle -------------------------------------------------------------------------------------------------
    def start(self) -> None:
        """Recover the runs of a previous server, then start dispatching."""
        with self._lock:
            if self._thread is not None:
                return
            self._stopping = False
            self._recover()
            self._thread = threading.Thread(target=self._dispatch, name="bss-jobs", daemon=True)
            self._thread.start()

    def stop(self, timeout: float = 20.0) -> None:
        """Stop dispatching; a running worker is killed and its run marked interrupted. Queued runs stay queued on disk
        and are queued again by the next start()."""
        with self._wake:
            self._stopping = True
            active = self._active
            if active is not None:
                active.interrupt = True
            self._wake.notify_all()
        if active is not None and active.proc is not None:
            kill_tree(active.proc)
        thread = self._thread
        if thread is not None:
            thread.join(timeout)
        with self._lock:
            self._thread = None
            self._queue.clear()

    def _recover(self) -> None:
        queued = []
        for record in self.storage.list_runs():
            if record.get("status") == "running":
                for step in record.get("steps") or []:
                    if step.get("status") == "running":
                        step["status"], step["reason"] = "failed", "interrupted"
                record.update(status="interrupted", pid=None,
                              error="The server stopped while this run was in progress; rerun it.")
                self._append_log(record["id"], "--- interrupted: the server stopped while the worker was running")
                try:
                    self.storage.save_run(record)
                except StorageError as exc:
                    log.warning("cannot mark run %s interrupted: %s", record.get("id"), exc)
            elif record.get("status") == "queued":
                queued.append(record)
        self._queue.extend(r["id"] for r in sorted(queued, key=created_key))

    # ---- public API ------------------------------------------------------------------------------------------------
    @property
    def active_run_id(self) -> str | None:
        with self._lock:
            return self._active.rid if self._active is not None else None

    @property
    def queued_count(self) -> int:
        with self._lock:
            return len(self._queue)

    def submit(self, rid: str) -> None:
        with self._wake:
            if rid not in self._queue:
                self._queue.append(rid)
            self._wake.notify_all()

    def is_busy(self, rid: str) -> bool:
        """True while the run is queued or running."""
        with self._lock:
            return rid in self._queue or (self._active is not None and self._active.rid == rid)

    def get_run(self, rid: str) -> dict:
        """The run record; the live one while the run is active."""
        with self._lock:
            if self._active is not None and self._active.rid == rid:
                return copy.deepcopy(self._active.record)
        return self.storage.get_run(rid)

    def list_runs(self, project_id: str | None = None) -> list[dict]:
        records = self.storage.list_runs(project_id)
        with self._lock:
            if self._active is not None:
                records = [copy.deepcopy(self._active.record) if r.get("id") == self._active.rid else r
                           for r in records]
        return records

    def cancel(self, rid: str, wait: float = 10.0) -> dict:
        """Cancel a run: a queued run leaves the queue, a running worker is killed with its process tree. Returns the
        record (final once the worker has exited, which normally takes well under `wait` seconds)."""
        with self._lock:
            if rid in self._queue:
                self._queue.remove(rid)
                record = self.storage.get_run(rid)
                record.update(status="cancelled", finished=now_iso())
                self.storage.save_run(record)
                self._append_log(rid, "--- cancelled before it started")
                return record
            active = self._active if self._active is not None and self._active.rid == rid else None
            if active is None:
                return self.storage.get_run(rid)
            active.cancel = True
            active.dirty = True
            proc = active.proc
        if proc is not None:
            kill_tree(proc)
        active.finished.wait(wait)
        return self.get_run(rid)

    # ---- dispatcher ------------------------------------------------------------------------------------------------
    def _dispatch(self) -> None:
        while True:
            with self._wake:
                while not self._queue and not self._stopping:
                    self._wake.wait()
                if self._stopping:
                    return
                rid = self._queue.popleft()
                try:
                    record = self.storage.get_run(rid)
                except StorageError as exc:
                    log.warning("queued run %s cannot be read and is dropped: %s", rid, exc)
                    continue
                if record.get("status") != "queued":
                    continue
                active = _Active(rid, record)
                self._active = active
            try:
                self._execute(active)
            except Exception as exc:  # the dispatcher must survive anything one run does
                log.exception("run %s: job runner error", rid)
                self._finalize(active, None, error=f"internal error of the job runner: {exc}")
            finally:
                with self._lock:
                    self._active = None
                active.finished.set()

    def command(self, run_dir: Path) -> list[str]:
        return [sys.executable, "-u", "-m", self.worker_module, str(run_dir)]

    def _execute(self, active: _Active) -> None:
        rid = active.rid
        run_dir = self.storage.run_dir(rid)
        try:
            cfg = self.storage.read_run_config(rid)
        except StorageError as exc:
            self._finalize(active, None, error=f"cannot read the configuration snapshot: {exc}")
            return
        deterministic = bool(cfg.get("deterministic", True)) if isinstance(cfg, dict) else True
        with self._lock:
            # a run interrupted before its worker started is still "queued" on disk: the next server start runs it
            if active.interrupt:
                return
            cancelled = active.cancel
            if not cancelled:
                active.record.update(status="running", started=now_iso(), error=None, traceback=None)
                active.t0 = time.monotonic()
        if cancelled:
            self._finalize(active, None)
            return
        cmd = self.command(run_dir)
        self._append_log(rid, f"--- started {active.record['started']}: {' '.join(cmd)} "
                              f"(deterministic: {'yes' if deterministic else 'no'})")
        popen_kwargs = {"creationflags": _WINDOWS_FLAGS} if os.name == "nt" else {"start_new_session": True}
        if sys.platform.startswith("linux"):
            popen_kwargs["preexec_fn"] = _parent_death_signal
        try:
            proc = subprocess.Popen(cmd, cwd=str(PROJECT_ROOT), env=worker_env(deterministic), stdin=subprocess.DEVNULL,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **popen_kwargs)
        except OSError as exc:
            self._finalize(active, None, error=f"cannot start the worker ({' '.join(cmd)}): {exc}")
            return
        with self._lock:
            active.proc = proc
            active.record["pid"] = proc.pid
            stop_now = active.cancel or active.interrupt
        self._write(active)
        if stop_now:
            kill_tree(proc)
        reader = threading.Thread(target=self._read_output, args=(active, proc.stdout, self.storage.log_path(rid)),
                                  name=f"bss-log-{rid}", daemon=True)
        reader.start()
        last_write = time.monotonic()
        while True:
            try:
                returncode = proc.wait(timeout=0.25)
                break
            except subprocess.TimeoutExpired:
                pass
            if active.dirty and time.monotonic() - last_write >= WRITE_INTERVAL:
                self._write(active)
                last_write = time.monotonic()
        # a grandchild that inherited the pipe keeps it open; do not wait for it forever
        reader.join(timeout=30)
        self._finalize(active, returncode)

    # ---- worker output ---------------------------------------------------------------------------------------------
    def _read_output(self, active: _Active, stream, log_path: Path) -> None:
        try:
            with open(log_path, "ab") as log_file:
                for raw in iter(stream.readline, b""):
                    line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
                    if "\r" in line:  # a progress bar redrawn with carriage returns: keep its last state
                        line = line.rsplit("\r", 1)[-1]
                    text = self._consume(active, line)
                    if text is not None:
                        log_file.write(text.encode("utf-8") + b"\n")
                        log_file.flush()
        except (OSError, ValueError) as exc:
            log.warning("run %s: reading the worker output failed: %s", active.rid, exc)
        finally:
            try:
                stream.close()
            except OSError:
                pass

    def _consume(self, active: _Active, line: str) -> str | None:
        """Apply the event on `line`, if any; returns the text to log (None when the line was only an event)."""
        i = line.find(EVENT_PREFIX)
        if i >= 0:
            rest = line[i + len(EVENT_PREFIX):]
            try:
                event, end = json.JSONDecoder().raw_decode(rest)
            except ValueError:
                event, end = None, 0
            if isinstance(event, dict):
                self._apply(active, event)
                # text another stream wrote onto the same line stays in the log
                line = (line[:i] + rest[end:]).strip()
                if not line:
                    return None
        if line.strip():
            active.tail.append(line)
        return line

    def _apply(self, active: _Active, event: dict) -> None:
        kind = event.get("type")
        with self._lock:
            rec = active.record
            steps = rec["steps"]

            def step_of(label) -> dict:
                label = str(label)
                for s in steps:
                    if s["label"] == label:
                        return s
                s = _new_step(label)
                steps.append(s)
                return s

            if kind == "plan":
                known = {s["label"]: s for s in steps}
                labels = list(dict.fromkeys(str(x) for x in event.get("steps") or []))
                steps[:] = [known.pop(lbl, None) or _new_step(lbl) for lbl in labels] + list(known.values())
            elif kind == "step":
                label = str(event.get("label", ""))
                for s in steps:  # a step starts when the previous one has finished
                    if s["status"] == "running" and s["label"] != label:
                        s["status"] = "done"
                s = step_of(label)
                s["status"], s["reason"] = "running", None
                # 1-based: "step 3 of 12"
                rec["progress"].update(index=_as_int(event.get("index"), 0) + 1, label=label,
                                       total=max(_as_int(event.get("total"), 0), len(steps)))
            elif kind == "step_done":
                s = step_of(event.get("label", ""))
                s["status"], s["seconds"] = "done", _as_float(event.get("seconds"))
            elif kind == "skip":
                s = step_of(event.get("label", ""))
                s["status"], s["reason"] = "skipped", str(event.get("reason") or "")
            elif kind == "warning":
                if len(rec["warnings"]) < MAX_WARNINGS:
                    rec["warnings"].append(str(event.get("message", "")))
            elif kind == "error":
                active.error = str(event.get("message") or "the analysis failed")
                rec["error"] = active.error
                rec["traceback"] = str(event.get("traceback") or "") or None
                for s in steps:
                    if s["status"] == "running":
                        s["status"], s["reason"] = "failed", active.error.splitlines()[0][:300]
            elif kind == "done":
                active.done = True
                for s in steps:
                    if s["status"] == "running":
                        s["status"] = "done"
            else:
                return
            _update_fraction(rec)
            active.dirty = True

    # ---- records ---------------------------------------------------------------------------------------------------
    def _write(self, active: _Active) -> None:
        with self._lock:
            snapshot = copy.deepcopy(active.record)
            active.dirty = False
        try:
            self.storage.save_run(snapshot)
        except StorageError as exc:
            log.warning("run %s: cannot write run.json: %s", active.rid, exc)

    def _finalize(self, active: _Active, returncode: int | None, error: str | None = None) -> None:
        with self._lock:
            rec = active.record
            if active.cancel:
                status = "cancelled"
            elif active.interrupt:
                status = "interrupted"
            elif error is None and returncode == 0 and active.done:
                status = "succeeded"
            else:
                status = "failed"
            if status == "failed":
                rec["error"] = error or active.error or self._tail_message(active, returncode)
            elif status == "interrupted":
                rec["error"] = "The server was stopped while this run was in progress; rerun it."
            running_to, reason = {"succeeded": ("done", None), "failed": ("failed", None),
                                  "cancelled": ("skipped", "cancelled"),
                                  "interrupted": ("failed", "interrupted")}[status]
            for s in rec["steps"]:
                if s["status"] == "running":
                    s["status"] = running_to
                    s["reason"] = reason or s["reason"] or (rec["error"] or "").splitlines()[0][:300] or None
                elif s["status"] == "pending" and status == "succeeded":
                    s["status"], s["reason"] = "skipped", "not reported by the pipeline"
            _update_fraction(rec)
            if status == "succeeded":
                total = len(rec["steps"])
                rec["progress"] = {"index": total, "total": total, "label": "", "fraction": 1.0}
        # the footer goes to the log before the final status is visible, so a reader that sees the run finished and
        # reads the log to its end has everything
        footer = f"--- {status}" + (f" (worker exit code {returncode})" if returncode is not None else "")
        if status == "failed" and error:
            footer += f": {error}"
        self._append_log(active.rid, footer)
        with self._lock:
            started = rec.get("started")
            rec.update(status=status, finished=now_iso(),
                       duration_s=round(time.monotonic() - active.t0, 1) if started else None)
            active.dirty = True
        self._write(active)

    @staticmethod
    def _tail_message(active: _Active, returncode: int | None) -> str:
        if active.tail:
            return "\n".join(active.tail)
        if returncode is None:
            return "The worker did not run."
        return f"The worker exited with code {returncode} without any output."

    def _append_log(self, rid: str, line: str) -> None:
        try:
            with open(self.storage.log_path(rid), "ab") as fh:
                fh.write(line.encode("utf-8") + b"\n")
        except (OSError, NotFound) as exc:
            log.warning("run %s: cannot append to log.txt: %s", rid, exc)


def _update_fraction(rec: dict) -> None:
    steps = rec["steps"]
    total = max(len(steps), _as_int(rec["progress"].get("total"), 0))
    finished = sum(s["status"] in ("done", "skipped", "failed") for s in steps)
    rec["progress"]["total"] = total
    rec["progress"]["fraction"] = round(finished / total, 4) if total else 0.0


def _as_int(value, default: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _as_float(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
