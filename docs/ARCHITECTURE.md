# Bone Spatial Studio — architecture and contracts

Bone Spatial Studio turns the manuscript pipeline (`bone_analysis.py`, `bone_figures.py`) into a local web
application: configure samples, parameters and gene signatures in the browser, draw a region of interest on the H&E
image, run the analysis in the background with live progress, and browse the figures, tables and per-bin maps.

Everything below is a contract between components. Names, signatures, file names and JSON shapes must match exactly.

```
popcorn/
  pipeline_config.py        configuration: defaults, schema, validation, path resolution, Progress     (shared contract)
  bone_analysis.py          analysis library + CLI          (fixed and generalised from the manuscript version)
  bone_figures.py           figure library + CLI            (fixed and generalised)
  bone_pending_masks.py     hand-drawn region of interest: GeoJSON polygons -> per-bin mask
  webapp/
    __init__.py
    __main__.py             `python -m webapp` = start the server (argparse: --host --port --workspace --no-browser --token)
    server.py               FastAPI application (API + static files)
    storage.py              projects and runs on disk (workspace/)
    jobs.py                 job queue: one analysis at a time in a subprocess, log capture, event parsing, cancel
    worker.py               `python -m webapp.worker <run_dir>`: runs the analysis and the figures for one run
    static/                 single-page front end (no build step, no CDN; works offline)
      index.html, css/app.css (base design system), css/project.css, css/run.css
      js/core.js (helpers), js/main.js (shell + routes)
      js/views/{projects,project,runs,run,help}.js, js/components/*.js
  tools/make_demo_data.py   synthetic Space Ranger-like data (two small "tibiae") for the demo and the tests
  tests/                    pytest: unit, end-to-end (demo data), API
  _original/                untouched copies of the two scripts as received
  requirements.txt, README.md, CHANGES.md, start.ps1, start.sh
```

Python 3.10+; package versions as stated in `bone_analysis.py` (scanpy 1.11.5, anndata 0.12.9, pandas 2.3.3, numpy
2.3.5, scipy 1.16.3, ...). Code must also run on pandas 3 / anndata 0.13 where that costs nothing (avoid chained
assignment, `inplace=True` on column subsets, relying on object dtype for strings).

---------------------------------------------------------------------------------------------------------------------

## 1. Configuration (`pipeline_config.py`, already written — read it)

`default_config()` returns the manuscript setup. `merge_with_defaults(cfg)`, `validate_config(cfg) -> (errors,
warnings)`, `effective_config(cfg)` (merge + raise on errors), `load_config`, `save_config`, `resolve_bin_dir(cfg,
sample, bin_size=None) -> Path|None`, `available_bin_sizes(cfg, sample)`, `bin_um('016um') -> 16.0`, `slugify`,
`sample_by_id`, `sample_ids`, `signature_label`, `PARAM_SCHEMA`, `STEP_LABELS`, `Progress`.

Key facts used everywhere:
* `cfg["samples"]`: ordered list of `{id, label, condition, color, path, microns_per_pixel}`; N >= 1 sections
  (default two). `cfg["reference_sample"]` is the denominator of every fold change / ratio of ratios.
* Gene sets: `cfg["signatures"][name] = {label, table_label, color, class, canonical, core, extended, competes,
  excluded}`; `cfg["primary"]` (ordered: signatures with a positivity call), `cfg["programs"]` (ordered: signatures
  competing in the tissue assignment), `cfg["classes"]` (drawn tissue classes, ordered; a signature whose class is not
  listed, e.g. "removed", is assigned but not drawn), `cfg["receptor_panel"]` (gene -> role), `cfg["adipocyte_genes"]`,
  `cfg["ubiquitous"]`.
* `cfg["steps"]` booleans decide what runs (see STEP_LABELS); `cfg["figures"]` drives bone_figures.
* Method options added with the fixes (defaults = manuscript): `scoring.depth_matching` ("cap" | "quantile"),
  `proximity.exclude_self` (false), `compartments.tumor_conditions` (["tumor"]) and `compartments.de_sections`
  ("all" | "tumor").
* `Progress(events=True)` prints `@@EVENT {json}` lines (plan, step, step_done, skip, warning, done, error) next to
  human-readable lines; `events=False` prints only the readable lines (CLI).

## 2. Python entry points

### bone_analysis.py
```python
planned_steps(cfg: dict) -> list[str]
    # Labels of every step run_analysis_from_config will report via progress.step()/progress.skip(), in order,
    # for this cfg (respecting cfg["steps"]). Must match the labels actually used.
run_analysis_from_config(cfg: dict, out, progress: Progress | None = None) -> anndata.AnnData
    # cfg: any partial or full config; it calls pipeline_config.effective_config(cfg) itself.
    # Writes every output of section 3 into `out` (created if needed). Raises on fatal errors.
run_analysis(root, out, roi_masks=None, mask_frame="hires", do_moran=True, do_sensitivity=False, project_8um=False)
    # Backward-compatible wrapper (old signature): default config with data_root=root.
```
CLI (old flags keep working): `python bone_analysis.py --root . --out results [--config cfg.json] [--skip-moran]
[--geneset-sensitivity] [--project-8um] [--roi-mask SAMPLE=PATH ...] [--mask-frame hires] [--events]`.
Module-level constants (SAMPLES, SIG, PRIMARY, PROGRAMS, RECEPTOR_PANEL, ADIPOCYTE_GENES, UBIQUITOUS, PALETTE,
UM_PER_PX, MIN_GENES, ...) stay importable; they are derived from / consistent with `pipeline_config.DEFAULT_CONFIG`.

### bone_figures.py
```python
planned_steps(cfg: dict) -> list[str]
make_figures(results_dir, cfg: dict | None = None, progress: Progress | None = None, root=None,
             manuscript_extras: bool = False) -> list[dict]
    # cfg None -> results_dir/config.json if present, else defaults. Draws every figure whose inputs exist into
    # results_dir/figures/ (PNG + PDF), each in its own try/except (a failure becomes progress.warn + a manifest
    # entry with "error"; the others still run). Writes figures/manifest.json and returns its "figures" list.
    # manuscript_extras=True also draws Figure 1a and the Figure 2c H&E crops (needs `root`); never in the web app.
```
CLI: `python bone_figures.py --results results [--root .] [--config cfg.json] [--manuscript-extras]`.
`planned_steps(cfg, manuscript_extras=False)` also plans the extras' labels; figure progress labels read
"Draw Figure …" / "Draw Supplementary Figure …".
Figures never read raw Space Ranger data and never need the .h5ad: they read the CSV outputs (the per-bin table has
coordinates and calls) and `config.json`.

### bone_pending_masks.py
```python
load_polygons(path_or_geojson) -> list[tuple[np.ndarray, list[np.ndarray]]]
    # GeoJSON FeatureCollection / Feature / Polygon / MultiPolygon -> [(outer ring Nx2, [hole rings])], (x, y) order
rasterise_drawn_roi(a, masks: dict[str, str | Path | dict], frame: str = "hires",
                    scalefactors: dict[str, dict] | None = None) -> tuple[np.ndarray, list[str]]
    # (bool mask of len a.n_obs, ids of sections without a region). Bin centres a.obsm["spatial"] (full-resolution
    # pixels, columns = pxl_col, pxl_row = x, y) are tested against the polygons of that bin's section after converting
    # polygon coordinates to full resolution: hires: / tissue_hires_scalef, lowres: / tissue_lowres_scalef, fullres:
    # unchanged. scalefactors default: a.uns["spatial"][sample]["scalefactors"]. Inside = inside any outer ring and
    # not inside one of its holes. A section without a region (or with an empty one) is all True; the caller warns.
```

### webapp/worker.py
`python -m webapp.worker <run_dir>` (cwd = project root): reads `<run_dir>/config.json`, makes
`Progress(events=True)`, calls `progress.plan(ba.planned_steps(cfg) + (bf.planned_steps(cfg) if
cfg["steps"]["figures"] else []))`, then `ba.run_analysis_from_config(cfg, <run_dir>/results, progress)`, then
`bf.make_figures(<run_dir>/results, cfg, progress)` if enabled, then `progress.done()`. On any exception:
`progress.error(str(exc), traceback.format_exc())` and exit code 1.

## 3. Outputs (`<out>/`, the run's `results/`)

`{B}` = analysis bin size as in the original names ("16um" for bin_size "016um": `f"{int(bin_um)}um"`);
`{Q}` = qc_bin_size likewise ("8um").

| file | content (columns) |
|---|---|
| `config.json` | the effective configuration of this run |
| `summary.json` | key numbers, see below |
| `manifest.json` | `{"files": [{"file", "kind": "table"|"object"|"data", "title", "description", "rows": int|null}]}` — one entry per file written (including itself is not needed) |
| `qc_threshold_sweep_{B}.csv` | sample, min_genes, n_retained, pct_retained, median_genes (analysis bins, before filtering) |
| `qc_threshold_sweep_{Q}.csv` | same at the QC bin size (Supplementary Figure S1, right panel), plus bin_size |
| `FigureS1_qc_per_bin_{Q}.csv.gz` | sample, bin_size, n_genes_by_counts, total_counts, pct_counts_mt — every bin of the QC bin size before filtering (S1 histograms). A section without square_{qc_bin_size} falls back to the analysis bin size (bin_size says which), with a warning |
| `Figure6b_receptor_panel_depth_matched.csv` | gene, role, sample, roi, mean_raw, pct_detected_raw, mean_dm, pct_detected_dm, n_positive_dm, n_bins, ci_lo_dm, ci_hi_dm, median_umi_per_bin, fold_dm |
| `FigureS4_adipocyte_genes_depth_matched.csv` | same columns |
| `Figure4g_6a_positivity_native_depth_roi.csv` | signature, sample, roi, n_positive, n_bins, pct_positive, ci_lo, ci_hi, mean_score, fold_vs_reference; plus fold_tumor_over_control (same value on every row of a signature) when there are exactly two samples. Both folds are ratios of the exact rates n_positive / n_bins (fold_vs_reference rounded to 4 decimals, fold_tumor_over_control to 2), not of the rounded pct_positive |
| `Figure4g_6a_positivity_depth_matched.csv` | signature, sample, pct_positive_mean, pct_positive_sd, n_positive_mean, n_bins, mean_score, ci_lo, ci_hi. pct_positive_sd = SD over the downsampling replicates only; ci_lo, ci_hi = Wilson 95% interval (percent) of round(n_positive_mean) of n_bins |
| `leiden_domain_screen.csv` | domain, n_bins, ubiquitous_fraction, median_genes_per_umi, top_markers, de_uninterpretable |
| `leiden_domain_markers.csv` | scanpy rank_genes_groups_df columns (group, names, scores, logfoldchanges, pvals, pvals_adj, pct_nz_group, pct_nz_reference) |
| `Figure4h_osteoclast_compartment_DE.csv` | group, names, scores, pvals, pvals_adj, log2fc_mean_cp10k, pct_group, pct_rest, selection_gene. log2fc_mean_cp10k = log2((mean CP10k in the group + 1) / (mean CP10k in the other compartment + 1)), CP10k = expm1(lognorm); pct_group / pct_rest = % of bins detecting the gene; selection_gene = gene of the split_by signature's scored list (defines the growth-plate compartment, higher there by construction). scanpy's logfoldchanges is not written (1e-9 floor: \|values\| of 19-32 for genes undetected in one compartment) |
| `Figure4h_compartment_summary.csv` | sample ("all" first, then one row per section), threshold (the split_by score threshold used: the compartments.percentile of all bins; also `uns["compartment_threshold"]` of the h5ad), threshold_is_zero, pct_bins_above_threshold, n_growth_plate, n_tumor_associated, n_non_growth_plate_control, n_growth_plate_outside_roi, n_tumor_associated_outside_roi, median_umi_growth_plate, median_umi_tumor_associated |
| `Figure3a_program_composition_native.csv`, `Figure3b_program_composition_depth_matched.csv` | sample, program, n_bins, pct_of_all_bins, pct_of_assigned_bins, sample_ambiguous_pct |
| `leiden_seed_stability_adjusted_rand_index.csv` | seed_a, seed_b, n_clusters_a, n_clusters_b, ARI |
| `FigureS5_domain_by_signature_zscores.csv` | leiden_domain, then one column per primary signature (label) |
| `FigureS5_unresolved_bins.csv` | sample, unresolved_bins, all_bins, percent, unresolved_bins_roi, all_bins_roi, percent_roi; first row sample "all" (all sections pooled, the number quoted before), then one row per section; the _roi columns count bins of the region of interest (in_marrow_roi). Older results have one row and only unresolved_bins, all_bins, percent |
| `Figure4i_observed_to_null_summary.csv` | sample, n_a, n_b, n_units_roi, ratio_obs_null_mean, ratio_obs_null_lo, ratio_obs_null_hi, p_greater_than_null_mean, median_nn_um, mean_nn_um, null_mean_um, status, ratio_of_ratios_vs_reference, note |
| `Figure4i_distance_distributions_observed_and_null.csv` | sample, source ("observed" / "null"), distance_um, count |
| `{from}_to_{to}_proximity.csv` | one per proximity pair (defaults `fibroblast_to_osteoblast_proximity.csv`, `macrophage_to_osteoclast_proximity.csv`, original columns with the signature names substituted) |
| `morans_I_per_section.csv` | sample, gene, I, pval_norm, var_norm, pval_norm_fdr_bh (when available; Benjamini-Hochberg over the genes with a defined pval_norm of that section, NaN for genes with zero variance there), rank_in_sample |
| `FigureS2_spatially_variable_genes_morans_I.csv` | gene, `moranI_<sample id>` per sample, mean_rank |
| `Table1_gene_signatures.csv` | signature, canonical_markers, signature_genes, genes_counted_toward_two_gene_rule, excluded_genes |
| `per_bin_signature_calls_{B}.csv.gz` | barcode, sample_id, section, x_um, y_um, in_region_of_interest, total_umi, assigned_signature, tissue_class, leiden_domain, `positive_<slug>` per primary signature, `positive_<b>_pooled_depthmatched`, `positive_<a>_pooled_depthmatched` (uncoupling pair; suffix `_pooled` for calls "pooled", no columns for "native") |
| `gene_set_sensitivity_three_gene_sets.csv` | optional (steps.sensitivity): gene_set, signature, sample, n_units_roi, n_positive_native, pct_native, ci_lo, ci_hi, pct_depth_matched_mean, pct_depth_matched_sd, genes_used, genes_missing (comma-separated gene lists of the variant; missing genes are also reported as warnings) |
| `joint_{B}.h5ad` | optional (steps.write_h5ad) |
| `labels_{Q}_projected.csv` | optional (steps.project_fine_bins) |

Per-bin conventions: `<slug>` of a signature = its label with spaces as underscores (neutrophil -> granulocyte,
skeletal_muscle -> skeletal_muscle); `assigned_signature` uses the same slug or "ambiguous"; `tissue_class` = the
assigned signature's class ("removed" when that class is not drawn), "unassigned" for ambiguous bins; `x_um`, `y_um`
= full-resolution pixel coordinates (pxl_col, pxl_row) times the section's microns per pixel; `section` = sample label.
A step that is disabled or cannot run leaves its columns out (or empty) and the reason goes to summary.json.

### summary.json
```json
{
  "version": 1,
  "bin_size": "016um",
  "reference_sample": "LL-2-1",
  "depth_match_target_umi": 23,
  "samples": [{"id": "LL-2-1", "label": "control", "condition": "control", "color": "#4C72B0",
               "n_bins_loaded": 0, "n_bins_qc": 0, "pct_retained": 0.0, "median_umi": 0.0, "median_genes": 0.0,
               "n_bins_roi": 0, "microns_per_pixel": 0.59}],
  "metrics": [{"group": "Positivity, native depth (ROI)", "label": "osteoclast — control", "value": 0.42,
               "unit": "%", "detail": "42 of 10,000 bins; Wilson 95% CI 0.30–0.57"}],
  "warnings": ["signature osteocyte: 2 of 5 genes missing from the probe set (Fgf23, Pdpn)"],
  "skipped": [{"step": "Moran's I (Supplementary Figure S2, slow)", "reason": "disabled"}],
  "timings": [{"step": "1. Load sections and quality control", "seconds": 12.3}],
  "environment": {"python": "3.12.16", "scanpy": "1.11.5", "anndata": "0.12.9", "numpy": "...", "pandas": "...",
                  "squidpy": "1.8.1" or null, "harmonypy": "0.0.10" or null},
  "finished": "2026-10-28T15:30:00"
}
```
`metrics` groups (as available): "Sections", "Positivity, native depth (ROI)", "Positivity, depth-matched (ROI)",
"Tissue assignment", "Distance to nearest partner", "Proximity", "Leiden domains" (including the unresolved share of
each section), "Host receptor panel".

### figures/manifest.json (written by bone_figures.make_figures)
```json
{"figures": [{"name": "Figure4g_bone_lineage_positivity", "group": "Figure 4", "title": "Figure 4g — bone-lineage positivity",
              "caption": "...", "png": "Figure4g_bone_lineage_positivity.png", "pdf": "Figure4g_bone_lineage_positivity.pdf",
              "error": null}]}
```
Figure names: `Figure3a_<i>_<slug(label)>_tissue_classes_denoise` (one per sample, i = 1-based), `Figure3b_composition_denoise`
(+ `Figure3b_composition_denoise.csv`), `Figure4a-f_bone_lineage_maps`, one per `cfg.figures.positivity_panels`
entry (defaults `Figure4g_bone_lineage_positivity`, `Figure6a_myeloid_positivity`), `Figure4h_osteoclast_compartment_DE`,
`Figure4i_osteoblast_osteoclast_distance`, `Figure6b_receptor_panel`, `FigureS1_quality_control`,
`FigureS2_spatially_variable_genes`, `FigureS4_adipocyte_genes`, `FigureS5_leiden_domain_heatmap`.

## 4. Workspace (default `<project root>/workspace`, env `BSS_WORKSPACE` or `--workspace`)
```
workspace/projects/<pid>/project.json      {"id", "name", "created", "updated", "config"}
workspace/projects/<pid>/rois/<sid>.geojson
workspace/runs/<rid>/run.json              run record (section 5, GET /api/runs/{rid})
workspace/runs/<rid>/config.json           effective config snapshot; roi.masks rewritten to absolute paths of rois/ below
workspace/runs/<rid>/rois/<sid>.geojson    copies of the project's drawn regions at submission time
workspace/runs/<rid>/log.txt               full stdout+stderr of the worker (event lines removed)
workspace/runs/<rid>/results/              outputs (section 3), results/figures/
workspace/demo_data/                       synthetic data from tools/make_demo_data.py
```
Ids: project `p-<8 hex>`; run `r-YYYYMMDD-HHMMSS-<4 hex>` (sortable). All JSON written atomically (tmp + os.replace).

## 5. HTTP API (FastAPI, JSON; errors are `{"detail": "message"}` with 4xx/5xx)

```
GET  /                      -> static/index.html            GET /static/...  -> static files
GET  /api/health            -> {"ok": true, "app_version": "1.0.0", "python": "3.12.16", "platform": "...",
                                "workspace": "...", "packages": {"scanpy": "1.11.5"|null, "anndata", "numpy", "pandas",
                                "scipy", "sklearn", "statsmodels", "harmonypy", "igraph", "leidenalg", "squidpy",
                                "pyarrow", "h5py", "matplotlib", "PIL"}, "missing_required": [..], "missing_optional": [..],
                                "active_run": "<rid>"|null, "queued": 0}
                               required: scanpy anndata numpy pandas scipy sklearn statsmodels igraph leidenalg pyarrow
                               h5py matplotlib PIL; optional: harmonypy squidpy (fallbacks exist). Versions are computed
                               once (importlib.metadata, never importing scanpy in the server process).
GET  /api/defaults          -> {"config": default_config(), "param_schema": PARAM_SCHEMA, "step_labels": STEP_LABELS}
GET  /api/docs/{name}       -> {"name", "markdown"} for name in readme | changes | architecture

GET  /api/projects          -> [{"id", "name", "created", "updated", "n_samples", "data_root",
                                 "last_run": {"id", "status", "created"} | null}]   newest first
POST /api/projects          body {"name": str, "config"?: dict}  -> project (config merged with defaults)
GET  /api/projects/{pid}    -> {"id", "name", "created", "updated", "config", "rois": ["<sid>", ...],
                                "validation": {"errors": [..], "warnings": [..]}}   (validate_config, no data access)
PUT  /api/projects/{pid}    body {"name"?: str, "config"?: dict} -> project + "validation": {"errors": [..], "warnings": [..]}
                               The config is merged with the defaults and saved even when validate_config reports
                               errors (drafts are allowed); runs refuse to start until the errors are fixed.
                               400 only when the body is not a JSON object.
DELETE /api/projects/{pid}  -> {"ok": true}            (runs keep their own snapshot and stay)
POST /api/projects/{pid}/duplicate -> project          ("<name> (copy)", same config and rois)
POST /api/projects/{pid}/validate  body {"config"?: dict} (unsaved edits; default = saved config) -> report:
     {"ok": bool, "errors": [str], "warnings": [str],
      "samples": [{"id", "bin_dir": str|null, "bin_sizes": ["008um", "016um"], "has_positions": bool,
                   "has_scalefactors": bool, "has_hires_image": bool, "microns_per_pixel": float|null,
                   "n_bins": int|null, "message": str}],
      "genes": {"n_features": int|null,
                "signatures": {"<name>": {"present": [..], "missing": [..], "core_present": int}},
                "receptor_panel": {"missing": [..]}, "adipocyte_genes": {"missing": [..]}}}
     n_bins / gene names come from filtered_feature_bc_matrix.h5 read with h5py (matrix/shape, matrix/features/name);
     fast, no scanpy.
POST /api/projects/{pid}/detect-samples  body {"data_root"?: str}
     -> {"data_root": str, "candidates": [{"folder": "SR007144_LL-2-1", "path": "SR007144_LL-2-1/outs/binned_outputs",
         "suggested_id": "LL-2-1", "bin_sizes": ["002um", "008um", "016um"]}]}
     scans data_root and its children (depth <= 4) for folders containing square_*um/filtered_feature_bc_matrix.h5;
     path is relative to data_root when possible; suggested_id strips a leading "SR<digits>_" run prefix.
GET  /api/fs/list?path=<dir>  -> {"path": abs, "parent": abs|null, "roots": ["C:\\", ...] | ["/"],
     "entries": [{"name", "path", "is_dir": true, "looks_like_spaceranger": bool}]}   directories only, sorted;
     empty path -> the user's home folder. looks_like_spaceranger: contains outs/, binned_outputs/ or square_*um/.
GET  /api/projects/{pid}/samples/{sid}/image        -> image/png  tissue_hires_image.png of the analysis bin
                                                       (falls back to any square_*um/spatial/); 404 if none
GET  /api/projects/{pid}/samples/{sid}/image-meta   -> {"width", "height", "frame": "hires", "scalefactors": {...},
                                                         "microns_per_pixel": float|null}
GET  /api/projects/{pid}/samples/{sid}/bins-preview -> {"x": [..], "y": [..], "n_total": int}
     in-tissue bin centres of the analysis bin size in hires pixels (pxl * tissue_hires_scalef), at most 60,000
     (evenly subsampled); 404 when positions are missing
GET    /api/projects/{pid}/roi/{sid}  -> GeoJSON FeatureCollection (404 if none)
PUT    /api/projects/{pid}/roi/{sid}  body FeatureCollection of Polygon features, coordinates [x, y] in hires image
                                      pixels, properties {"frame": "hires"} -> {"ok": true, "n_polygons": int}
DELETE /api/projects/{pid}/roi/{sid}  -> {"ok": true}

POST /api/projects/{pid}/runs   body {"label"?: str, "steps"?: {key: bool}, "deterministic"?: bool}
     -> run record (status "queued"); 409 {"detail"} if the config has validation errors or required packages are
        missing. Snapshot = saved project config with the overrides applied; drawn ROIs copied into the run.
GET  /api/runs?project=<pid>    -> [run record without "config"] newest first
GET  /api/runs/{rid}            -> run record:
     {"id", "project_id", "project_name", "label", "status": "queued"|"running"|"succeeded"|"failed"|"cancelled"|"interrupted",
      "created", "started", "finished", "duration_s", "pid",
      "progress": {"index": int, "total": int, "label": str, "fraction": float},
      "steps": [{"label", "status": "pending"|"running"|"done"|"skipped"|"failed", "seconds": float|null, "reason": str|null}],
      "warnings": [str], "error": str|null, "traceback": str|null, "config": {...}}
GET  /api/runs/{rid}/log?offset=<bytes>  -> {"offset": int, "text": str, "complete": bool}  (complete = run finished
                                             and offset reached the end); at most 256 KB per call
POST /api/runs/{rid}/cancel     -> run record   (queued: removed from the queue; running: process tree killed)
POST /api/runs/{rid}/rerun      -> new run record (same config snapshot, same ROI copies)
DELETE /api/runs/{rid}          -> {"ok": true}  409 while queued or running

GET  /api/runs/{rid}/summary    -> summary.json (404 until written)
GET  /api/runs/{rid}/manifest   -> {"tables": [{"file", "title", "description", "rows", "size"}],
                                    "figures": [figures/manifest.json entries + {"png_url", "pdf_url"}],
                                    "objects": [{"file", "title", "description", "size"}],
                                    "other": [{"file", "size"}]}      built from the two manifests + a directory scan
GET  /api/runs/{rid}/files/{path:path}[?download=1]  -> the file (path inside results/; traversal rejected)
GET  /api/runs/{rid}/table?file=<csv or csv.gz>&offset=0&limit=100&sort=<col>&desc=0&q=<text>
     -> {"file", "columns": [..], "dtypes": {col: "number"|"text"|"bool"}, "rows": [[..]], "total": int,
         "offset", "limit"}   q = case-insensitive substring over all columns; NaN -> null
GET  /api/runs/{rid}/bins/samples -> [{"id", "label", "n"}]
GET  /api/runs/{rid}/bins?sample=<sid>
     -> {"sample", "n", "x": [..], "y": [..],                                    x_um, y_um
         "fields": {"tissue_class": {"type": "category", "categories": [..], "codes": [..]},
                    "assigned_signature": {...category...}, "leiden_domain": {...category...},
                    "in_region_of_interest": {"type": "bool", "values": [0|1, ..]},
                    "positive_<slug>": {"type": "bool", "values": [..]}, ...,
                    "total_umi": {"type": "number", "values": [..]}},
         "colors": {"tissue_class": {name: "#hex"}, "assigned_signature": {slug: "#hex"}}}
     from per_bin_signature_calls_{B}.csv.gz (cached per run); colors from the run's config.json (classes, signatures;
     "removed" "#b8b8b8", "unassigned"/"ambiguous" "#dedede").

POST /api/demo                  -> {"project": project}   generates workspace/demo_data (once; ~1 min) with
                                   tools.make_demo_data.generate and creates a project "Demo — synthetic tibiae"
                                   whose config points at it (cheap settings: n_perm 20, n_depth_reps 2, stability seeds
                                   [0, 1], moran.n_genes 500)
```
Details fixed during implementation: `progress.index` is the 1-based number of the current step; project create and
duplicate return the GET shape; PUT of an empty ROI FeatureCollection deletes the region; POST runs also answers 409
when `roi.use_drawn` is on and a drawn region's frame differs from `roi.mask_frame` or a mask file is missing; 422 for
an unreadable data file; `/table` returns at most 1,000 rows per call; in `/bins` a category code of -1 is a missing
value; `png_url`/`pdf_url` are null when the file is missing; every error `detail` is a string.

Security: binds 127.0.0.1 by default. With `--host` other than loopback a token is required (generated if not given and
printed with the URL); requests carry it as cookie `bss_token` (set by visiting `/?token=...`) or header `X-Token`.

## 6. Job runner (webapp/jobs.py)
* One analysis at a time (memory); FIFO queue; a daemon thread starts the next run when the previous ends.
* Subprocess: `[sys.executable, "-u", "-m", "webapp.worker", run_dir]`, cwd = project root, stdout+stderr merged,
  env adds PYTHONUNBUFFERED=1, MPLBACKEND=Agg, and when deterministic: PYTHONHASHSEED=0, VISIUM_DETERMINISTIC=1 and
  OMP/MKL/OPENBLAS/NUMEXPR/NUMBA_NUM_THREADS=1 (otherwise VISIUM_DETERMINISTIC=0). POSIX: start_new_session=True,
  cancel = killpg(SIGTERM, then SIGKILL after 5 s). Windows: CREATE_NEW_PROCESS_GROUP, cancel = `taskkill /T /F /PID`.
* Each stdout line: `@@EVENT ` lines update run.json (plan -> steps pending; step -> running + progress; step_done ->
  done + seconds; skip -> skipped + reason; warning -> warnings; error -> error/traceback; done); every other line is
  appended to log.txt.
* Final status: cancelled if cancel was requested; succeeded if exit code 0 and a done event arrived; else failed
  (error = error event message, or the last 20 log lines).
* On server start: runs left "running" become "interrupted"; runs left "queued" are queued again.

## 7. Front end (webapp/static, ES modules, no framework, no network dependencies)
`core.js` (already written — read it) provides: api(), apiUrl(), ApiError, h(), svg(), mount(), clear(), icon(),
button(), card(), field(), switchToggle(), tabs(), statusBadge(), loading(), emptyState(), errorPanel(), banner(),
toast(), modal(), confirmDialog(), promptDialog(), fmt.{int,num,pct,bytes,duration,date,ago}, debounce(), sleep(),
poll(), copyText(), clone(), state.health, refreshHealth(), route(), navigate(), setPathSilently(), currentPath().
`css/app.css` (already written — read it) defines the design system: .page, .page-header, .card, .btn (primary,
danger, ghost, small, icon-only), .input/.select/.textarea, .field, .form-grid, .switch, .tabs, .badge(-status),
.chip(.ok/.warn/.bad/.core), .table-wrap/.table, .kv, .metric, .progress, .log, .banner, .modal, .toast, .empty,
.grid-2/.grid-3/.grid-auto, .row, .stack, .legend, .prose. Use them; put view-specific rules in css/project.css
(projects, project editor, ROI editor) or css/run.css (runs, results, viewers). The interface has one dark theme; its
colors come from the CSS variables on `:root` only — never hard-code colors except data colors (signature/sample
palettes).

Routes (main.js, already written) and the view modules they call — each export renders into `container` and may
return a cleanup function (stop timers, remove listeners):
```
#/projects                 views/projects.js  export async function renderProjects(container, {})
#/projects/:id/:tab?       views/project.js   export async function renderProject(container, {id, tab})
                           tabs: data | parameters | signatures | roi | run
#/runs                     views/runs.js      export async function renderRuns(container, {})
#/runs/:id/:tab?           views/run.js       export async function renderRun(container, {id, tab})
                           tabs: progress | summary | figures | tables | viewer | config  (tab null -> progress
                           while queued/running, summary when finished)
#/help/:topic?             views/help.js      export async function renderHelp(container, {topic})
                           topics: guide (default) | changes | architecture | readme
```
Components live in js/components/ and are imported by the views.
