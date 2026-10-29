# Bone Spatial Studio

The Visium HD bone-lesion analysis of the manuscript (`bone_analysis.py`, `bone_figures.py`), fixed and turned into a
local web application. Configure the sections, parameters and gene signatures in the browser, draw a region of interest
on the H&E image, run the analysis in the background with live progress, and browse the summary, figures, tables and
per-bin spatial maps. The command line keeps working as before.

What was wrong in the original scripts, and which manuscript numbers change because of it: **[CHANGES.md](CHANGES.md)**.

## Start

**Windows** (PowerShell, in this folder):

```
powershell -ExecutionPolicy Bypass -File .\start.ps1
```

**macOS / Linux**:

```
./start.sh
```

The first start creates `.venv`, installs `requirements.txt` (a few minutes) and opens http://127.0.0.1:8765.
Python 3.10 or newer is needed (3.12 recommended): https://www.python.org/downloads/ or
`winget install Python.Python.3.12`. Arguments are passed on to the server, e.g. `.\start.ps1 --port 9000`.

Without the scripts: `python -m pip install -r requirements.txt`, then `python -m webapp`
(`--port`, `--workspace`, `--no-browser`, `--host`, `--token`). The exact versions the app was tested with are in
`requirements-lock.txt`.

## Use

1. **Projects → New project** (or **Try the demo**, which generates two small synthetic tibiae).
2. **Data**: choose the folder that holds the Space Ranger outputs (e.g. `data/`, containing `SR007144_LL-2-1/` and
   `SR007144_RL-1-5/`), **Detect samples**, pick the reference (control) section; the data check confirms the bin
   folders, images, microns per pixel and which signature genes are in the probe set.
3. **Parameters** and **Signatures**: the defaults are the manuscript setup; every value has a help text and a
   reset button. Gene signatures (Table 1) can be edited, added, renamed or removed.
4. **Region of interest** (optional): draw polygons on the H&E image of each section (D draw, Enter close, Esc cancel,
   Backspace remove a corner, Del remove the polygon, Ctrl+Z undo, Ctrl+S save) and switch on "Use the drawn region".
5. **Run → Start analysis**. The run page shows each step, the live log, then the summary, figures (PNG/PDF), tables
   (sortable, searchable, downloadable), a spatial viewer of every bin and the exact configuration used.

Runs, projects and drawn regions are stored in `workspace/` (or `--workspace`, or the `BSS_WORKSPACE` variable).
Each run keeps a frozen copy of its configuration, so it can be reproduced with **Re-run**.

On a remote server, keep the default `--host 127.0.0.1` and open an SSH tunnel
(`ssh -L 8765:127.0.0.1:8765 user@server`); with `--host 0.0.0.0` a token is required and printed with the URL.

## Command line

```
python bone_analysis.py --root data --out results/run1        # the manuscript analysis
python bone_figures.py --results results/run1                 # its figures
python pipeline_config.py > my_config.json                    # write the default configuration
python bone_analysis.py --config my_config.json --out results/run2
python tools/compare_with_original.py --data data --results results/run1   # fixed vs original code
```

Old flags still work: `--skip-moran`, `--geneset-sensitivity`, `--project-8um`, `--roi-mask SAMPLE=PATH.geojson`.
The outputs and their columns are listed in [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md), section 3.

## Files

| Path | What |
|---|---|
| `bone_analysis.py` | analysis pipeline (fixed; `[FIX]` / `[APP]` notes mark every change) |
| `bone_figures.py` | figures (fixed) |
| `bone_pending_masks.py` | hand-drawn region of interest (GeoJSON → per-bin mask) |
| `pipeline_config.py` | the configuration: defaults = manuscript setup, validation, parameter help |
| `webapp/` | the web application (FastAPI server, job runner, browser interface) |
| `tools/make_demo_data.py` | synthetic Visium HD data for the demo and the tests |
| `tools/compare_with_original.py` | reruns the original code on the same data and compares |
| `tests/` | `python -m pytest` (84 tests; the slow end-to-end ones run the whole analysis on tiny data) |
| `_original/` | the two scripts exactly as received |

## Notes

* Reproducibility: runs are single-threaded by default (`deterministic`), so a rerun reproduces every number.
* `harmonypy` and `squidpy` are optional (PCA without Harmony; a built-in Moran's I with the same numbers).
* `pynndescent` must stay below 0.6 with scanpy 1.11 (pinned in `requirements.txt`): newer versions make scanpy's
  neighbour graph fail on data of 4,096 bins or more.
