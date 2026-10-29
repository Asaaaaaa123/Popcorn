"""
`python -m webapp`: start Bone Spatial Studio and open it in the browser.

    python -m webapp [--host 127.0.0.1] [--port 8765] [--workspace DIR] [--token TOKEN] [--no-browser]

The workspace (projects, runs, demo data) defaults to $BSS_WORKSPACE, else <project root>/workspace. The server listens
on this computer only (127.0.0.1) unless --host says otherwise; then every request must carry an access token, which
is generated when --token is not given and printed as part of the address to open.
"""
from __future__ import annotations

import argparse
import ipaddress
import logging
import os
import secrets
import socket
import sys
import threading
import time
import webbrowser
from pathlib import Path

from webapp import __version__

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8765
BROWSER_WAIT = 60.0  # seconds to wait for the server to accept connections before giving up on the browser


def is_loopback(host: str) -> bool:
    if host.lower() == "localhost":
        return True
    try:
        return ipaddress.ip_address(host.strip("[]")).is_loopback
    except ValueError:
        return False


def url_host(host: str) -> str:
    """Host part of a URL for `host`; a wildcard address is reached through the loopback address."""
    host = host.strip("[]")
    if host in ("", "0.0.0.0"):
        return "127.0.0.1"
    if host == "::":
        return "[::1]"
    return f"[{host}]" if ":" in host else host


def token_arg(value: str) -> str:
    from webapp.server import TOKEN_RE

    if not TOKEN_RE.match(value):
        raise argparse.ArgumentTypeError("use 1-256 letters, digits, '.', '_', '~' or '-'")
    return value


def port_arg(value: str) -> int:
    try:
        port = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a port number: {value!r}") from None
    if not 1 <= port <= 65535:
        raise argparse.ArgumentTypeError("the port must be between 1 and 65535")
    return port


def bind_problem(host: str, port: int) -> str | None:
    """Why the server cannot listen on host:port, or None when it can."""
    family = socket.AF_INET6 if ":" in host.strip("[]") else socket.AF_INET
    with socket.socket(family, socket.SOCK_STREAM) as s:
        if os.name != "nt":  # as uvicorn does; on Windows SO_REUSEADDR would let two servers share the port
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            s.bind((host.strip("[]"), port))
        except socket.gaierror as exc:
            return f"unknown host {host!r} ({exc})"
        except OSError as exc:
            return (f"cannot listen on {host}:{port} ({exc.strerror or exc}). If Bone Spatial Studio is already "
                    f"running, open http://{url_host(host)}:{port}/; otherwise choose another port with --port.")
    return None


def open_browser_when_ready(url: str, host: str, port: int) -> None:
    def wait_and_open():
        deadline = time.monotonic() + BROWSER_WAIT
        while time.monotonic() < deadline:
            try:
                with socket.create_connection((host.strip("[]"), port), timeout=1):
                    break
            except OSError:
                time.sleep(0.25)
        else:
            return
        try:
            webbrowser.open(url)
        except webbrowser.Error as exc:
            print(f"Could not open a browser ({exc}); open {url} yourself.", flush=True)

    threading.Thread(target=wait_and_open, name="bss-browser", daemon=True).start()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m webapp", description="Bone Spatial Studio: a local web application "
                                 "for the Visium HD bone-lesion analysis.")
    ap.add_argument("--host", default=DEFAULT_HOST,
                    help=f"address to listen on (default {DEFAULT_HOST}: this computer only; 0.0.0.0 = every "
                         f"network interface, which requires an access token)")
    ap.add_argument("--port", type=port_arg, default=DEFAULT_PORT, help=f"port (default {DEFAULT_PORT})")
    ap.add_argument("--workspace", default=os.environ.get("BSS_WORKSPACE") or str(PROJECT_ROOT / "workspace"),
                    help="folder for projects, runs and demo data (default $BSS_WORKSPACE, else "
                         "<project root>/workspace)")
    ap.add_argument("--token", type=token_arg, default=None,
                    help="access token required on every request (generated when --host is not a loopback address)")
    ap.add_argument("--no-browser", action="store_true", help="do not open the browser")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    token = args.token or (None if is_loopback(args.host) else secrets.token_urlsafe(24))
    workspace = Path(args.workspace).expanduser()

    problem = bind_problem(args.host, args.port)
    if problem:
        print(f"error: {problem}", file=sys.stderr)
        return 1

    from webapp.server import create_app

    try:
        app = create_app(workspace, token=token)
    except OSError as exc:
        print(f"error: cannot use the workspace {workspace}: {exc}", file=sys.stderr)
        return 1

    query = f"?token={token}" if token else ""
    url = f"http://{url_host(args.host)}:{args.port}/{query}"
    print(f"Bone Spatial Studio {__version__}", flush=True)
    print(f"  workspace: {app.state.storage.root}", flush=True)
    print(f"  open:      {url}", flush=True)
    if not is_loopback(args.host):
        print(f"  from other computers: http://{socket.gethostname()}:{args.port}/{query}", flush=True)
        print("  Anyone with this address can use the server and read the files it can read.", flush=True)
    print("Press Ctrl+C to stop.", flush=True)
    if not args.no_browser:
        open_browser_when_ready(url, url_host(args.host), args.port)

    import uvicorn

    from webapp.jobs import kill_workers_with_server
    kill_workers_with_server()
    uvicorn.run(app, host=args.host.strip("[]"), port=args.port, log_level="info", access_log=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
