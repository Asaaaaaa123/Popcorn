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
3. **Explore** (or click the bin count of a section in the Data tab): read the raw data before running anything,
   with the settings as they are (saved or not). A table gives the bins, genes and UMIs of each section, the median
   UMIs, genes and mitochondrial share per bin and the bins that pass the quality control of the Parameters tab; the
   maps show these per-bin numbers or any gene on the H&E image, in one, two or four linked panels; the gene list of
   a section is ranked by UMIs and searchable: click a gene to put it on the map, or its layers button to add it as
   an overlay.
4. **Parameters** and **Signatures**: the defaults are the manuscript setup; every value has a help text and a
   reset button. Gene signatures (Table 1) can be edited, added, renamed or removed (see
   [Using it for other gene groups](#using-it-for-other-gene-groups)).
5. **Region of interest** (optional): draw polygons on the H&E image of each section (D draw, Enter close, Esc cancel,
   Backspace remove a corner, Del remove the polygon, Ctrl+Z undo, Ctrl+S save) and switch on "Use the drawn region".
6. **Run → Start analysis**. The run page shows each step, the live log, then the summary, figures (PNG/PDF), tables
   (sortable, searchable, downloadable; **Export** at the top downloads all results, the figures or the tables as
   one zip), a spatial viewer of every bin (tissue classes, calls, Leiden domains or any gene; one, two or four
   linked panels, overlays of calls, region and genes, H&E underneath) and the exact configuration used.

Runs, projects and drawn regions are stored in `workspace/` (or `--workspace`, or the `BSS_WORKSPACE` variable).
Each run keeps a frozen copy of its configuration, so it can be reproduced with **Re-run**.

On a remote server, keep the default `--host 127.0.0.1` and open an SSH tunnel
(`ssh -L 8765:127.0.0.1:8765 user@server`); with `--host 0.0.0.0` a token is required and printed with the URL.

## Using it for other gene groups

The defaults are the manuscript's: mouse bone and marrow signatures (Table 1) and a few steps built around osteoclasts
and osteoblasts in an A549 bone lesion. Quality control, Leiden domains, scores and positivity calls, depth matching,
the region of interest, the tissue assignment (Figure 3), Moran's I, the per-bin table, the summary and the maps work
with any signatures.

1. **Signatures** tab: delete the signatures you do not need and add your own (scored genes, core genes, competing
   signatures, and the switches *Positivity call* and *Tissue assignment*); replace the tissue classes. Deleting a
   signature removes it from every list; the dialog names the settings that still select it and where to change them.
2. The **configuration check** (top of the project page; it checks your unsaved settings as you edit them, while runs
   use the saved ones) lists every setting that still names a missing signature:
   * errors, the analysis cannot start: *Osteoclast compartments (Figure 4h)* and *Distance to nearest partner
     (Figure 4i)* need signatures with a positivity call. Choose yours under **Parameters**, or switch the step off
     under **Parameters → Analysis steps**.
   * warnings, the figure or table is skipped: empty signature lists of the Figure 3b composition bars and the
     positivity panels, the Figure 4a–f maps (**Parameters → Figures**); Figure 6b genes that are not in the host
     receptor panel; no proximity pairs; no section whose condition is a tumor condition (*Osteoclast compartments →
     Tumor conditions*); habitat classes of the Figure 4i sensitivity analysis that none of your signatures has; the
     A549 lesion genes still listed in a project without the osteoclast and osteoblast signatures.
3. Steps of the bone/A549 study are marked *bone study* under **Parameters → Analysis steps**: host receptor panel
   (Figure 6b), adipocyte genes (Supplementary Figure S4), osteoclast compartments (Figure 4h), distance to the nearest
   partner (Figure 4i; osteoblast to osteoclast in the manuscript) and the three-gene-set sensitivity (off by default;
   it tests only the osteoclast and osteoblast gene sets). The receptor panel and the adipocyte genes take any gene list
   (Signatures tab); with a panel of your own, also choose the genes Figure 6b draws (**Parameters → Figures →
   Receptor genes shown in Figure 6b**) and switch off the shading of the first two (RANKL and OPG). Figure 4h calls its
   compartments growth plate and tumor-associated only for the manuscript's osteoclast bins split by the chondrocyte
   score; any other split is named after the split signature (*high-endothelium* and *low-endothelium*, say). In the
   sensitivity analyses, the Figure 4h and 4i parts run only with their steps, and the A549 lesion flag counts the
   cross-captured human genes Krt7 and Krt18: empty that list (*Sensitivity analyses: A549 lesion flag*) in other
   studies, where they are ordinary epithelial genes (the run log says how many bins were flagged).
4. Gene names must match the probe set exactly. The defaults are mouse symbols (`Lyz2`, `Ptprc`); human data needs
   human symbols (`LYZ`, `PTPRC`), also in the ubiquitous transcripts of the Leiden screen (Signatures tab). The data
   check lists the missing genes of every signature, the receptor panel, the adipocyte genes, the ubiquitous
   transcripts and (with the sensitivity analyses on) the A549 lesion genes; a run warns when none of the ubiquitous
   transcripts is in the data. A positivity panel leaves out, and names in its caption, a signature with fewer than two
   of its genes in the data (it cannot be scored). Mitochondrial genes are those whose name starts with `mt-` in any
   case, so `MT-` works too. File names (e.g. `Figure4h_osteoclast_compartment_DE.csv`) stay those of the manuscript;
   figure titles and captions use your signature labels.

## Viewing results in Loupe Browser

**Export → Loupe Browser** (at the top of a finished run) downloads `<run>_loupe.zip`: the run's per-bin calls as
category files for 10x Genomics Loupe Browser, one CSV per section and bin size, and a README.txt that says which file
goes with which `.cloupe` file.

* `<section>_008um_categories.csv` is for the section's own Loupe file, `outs/cloupe_008um.cloupe`, which Space Ranger
  writes for every Visium HD run (at 8 µm only). Each 8 µm bin carries the calls of the 16 µm bin that holds it (row
  and column divided by 2, checked against the bin centres of the tissue positions files; on the manuscript data every
  8 µm centre lies inside its 16 µm bin), whatever its own depth: on the manuscript data 115,048 of the 223,186
  labelled 8 µm bins of LL-2-1 and 154,257 of the 373,080 of RL-1-5 would fail quality control on their own (257 and
  477 have no UMI). The last column, **8 um QC**, says *pass* or *fail* for each 8 µm bin by the run's quality-control
  thresholds: leave out the *fail* bins before comparing groups in Loupe Browser. The run's own projection onto 8 µm
  bins (step *Project labels onto the QC bin size*) keeps only the bins that pass.
* `<section>_016um_categories.csv` (the run's own bin size) fits a 16 µm `.cloupe` file, which exists only if Space
  Ranger 3.1 or later was run with `--custom-bin-size 16`.

In Loupe Browser 8 or later, open the section's `.cloupe` file, click the upload icon next to "+ Create a new group" in
the Clusters panel (older versions: Categories → action menu → Import Categories) and choose the CSV. Each column
becomes a group: tissue class, assigned signature, Leiden domain ("domain 0", ...), region of interest
(inside/outside) and one positivity call (positive/negative) per signature, plus any other yes/no or text call of the
run; numbers such as UMIs are left out. The columns follow the run's configuration, so a study with other signatures
exports the same way. Every bin of the `.cloupe` file has a label: bins the run did not analyse (removed by quality
control; in an 8 µm file, bins whose 16 µm bin was removed) are "not analyzed".

The files follow 10x Genomics' description of category CSVs, choosing the most compatible option where it is silent:
lower-case `barcode` header first, barcodes exactly as in `filtered_feature_bc_matrix.h5` (with `-1`), a label in every
cell (no empty cells), no commas or quotes inside labels, LF line endings, UTF-8 without BOM. They have not been tried
in Loupe Browser itself (no `.cloupe` file was available), and the menu names above may differ between versions.

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
| `tests/` | `python -m pytest` (the slow end-to-end ones run the whole analysis on tiny data; `-m "not slow"` skips them) |
| `_original/` | the two scripts exactly as received |

## Notes

* Reproducibility: runs are single-threaded by default (`deterministic`), so a rerun reproduces every number.
* `harmonypy` and `squidpy` are optional (PCA without Harmony; a built-in Moran's I with the same numbers).
* `pynndescent` must stay below 0.6 with scanpy 1.11 (pinned in `requirements.txt`): newer versions make scanpy's
  neighbour graph fail on data of 4,096 bins or more.
