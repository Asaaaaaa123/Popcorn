#!/usr/bin/env bash
# Start Bone Spatial Studio on macOS or Linux.
#
#   ./start.sh                            (or: bash start.sh)
#   ./start.sh --port 8800 --no-browser   every argument is passed on to "python -m webapp"
#
# The first start creates a virtual environment in .venv (with uv when it is installed, else with python3.12, python3
# or python) and installs requirements.txt; requirements.txt is installed again whenever it changes.
set -euo pipefail

cd "$(dirname "$0")"
VENV=.venv
PY="$VENV/bin/python"
STAMP="$VENV/.requirements.sha256"

find_python() {
  local candidate
  for candidate in python3.12 python3 python; do
    if command -v "$candidate" >/dev/null 2>&1 &&
       "$candidate" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' >/dev/null 2>&1; then
      echo "$candidate"
      return 0
    fi
  done
  return 1
}

fail() {
  echo "$1" >&2
  exit 1
}

if [ ! -x "$PY" ]; then
  if command -v uv >/dev/null 2>&1; then
    echo "Creating the virtual environment .venv with uv (Python 3.12) ..."
    uv venv --seed --python 3.12 "$VENV" || fail "Creating the virtual environment failed (see above). Delete the .venv folder and try again."
  else
    BASE=$(find_python) || fail "Python 3.10 or newer (3.12 recommended) is required but was not found.
Install it from https://www.python.org/downloads/ or with your package manager
(macOS: brew install python@3.12; Debian/Ubuntu: sudo apt install python3 python3-venv), then run ./start.sh again."
    echo "Creating the virtual environment .venv with $BASE ..."
    "$BASE" -m venv "$VENV" || fail "Creating the virtual environment failed (see above; on Debian/Ubuntu install python3-venv). Delete the .venv folder and try again."
  fi
fi

HASH=$("$PY" -c 'import hashlib; print(hashlib.sha256(open("requirements.txt", "rb").read()).hexdigest())')
if [ ! -f "$STAMP" ] || [ "$(cat "$STAMP")" != "$HASH" ]; then
  echo "Installing the requirements (the first time this takes several minutes) ..."
  if command -v uv >/dev/null 2>&1; then
    uv pip install --python "$PY" -r requirements.txt || fail "Installing the requirements failed (see above). Check the internet connection and run ./start.sh again."
  else
    "$PY" -m pip install --disable-pip-version-check -r requirements.txt || fail "Installing the requirements failed (see above). Check the internet connection and run ./start.sh again."
  fi
  echo "$HASH" > "$STAMP"
fi

# ${1+"$@"}: the arguments, also with none under "set -u" in the bash 3.2 of macOS
exec "$PY" -m webapp ${1+"$@"}
