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
`sample_by_id`, `sample_ids`, `signature_label`, `compartment_names(cfg) -> (high, other)` (the Figure 4h compartments
in figure text and messages: "growth-plate", "tumor-associated" for the manuscript's osteoclast / chondrocyte split,
else "high-<split_by label>", "low-<split_by label>"; file and column names keep growth_plate / tumor_associated),
`PARAM_SCHEMA`, `STEP_LABELS`, `Progress`.

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
* Settings that name single signatures are checked against `cfg["signatures"]` (other gene groups, README "Using it
  for other gene groups"), only while the step that uses them is on. Errors (the run is refused), each naming the
  setting, its Parameters field and both ways out (another signature, or the step switched off):
  `compartments.signature` (needs a positivity call) and `compartments.split_by` (steps.compartments),
  `uncoupling.a` / `.b` (positivity calls; steps.uncoupling). Warnings (the output is skipped): `figures.lineage_maps`
  unknown or without a call, `figures.composition_signatures` with no known signature (or unknown names), a positivity
  panel without a called signature, `proximity.pairs` empty (or a pair naming an unknown or, for mode positivity, an
  uncalled signature), no section whose condition is in `compartments.tumor_conditions`, habitat classes of
  `sensitivity_analyses.uncoupling` that no tissue-assignment signature has, `sensitivity.signatures` with no known
  signature while steps.sensitivity is on, `figures.receptor_genes` naming no gene of `receptor_panel` (steps
  receptor_panel and figures on: Figure 6b is skipped), `sensitivity_analyses.lesion.genes` not empty in a configuration
  without the osteoclast and osteoblast signatures (steps.sensitivity_analyses on: Krt7 / Krt18 would flag ordinary
  epithelium as A549 lesion).
* `STEP_HELP` (what each step is for) and `BONE_STUDY_STEPS` (receptor_panel, adipocyte_panel, compartments,
  uncoupling, sensitivity: specific to the bone-lesion study). `PARAM_SCHEMA` ends with a section `{"key": "steps"}`
  that is not a configuration section: one bool field per STEP_LABELS key with its help and `"study": "bone"`; the
  Analysis steps card and the Run tab read it. Signature fields may carry `"uses_step"` (a steps key) and
  `"needs_call": true`; `compartments.tumor_conditions` is a str_list field.

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
An enabled step that has nothing to compute is reported through `progress.skip(label, reason)`, not as done: Osteoclast
compartments without `pos_<signature>` / `score_<split_by>`, Distance to nearest partner without the calls of `a` and
`b`, Proximity tables without pairs, Gene-set sensitivity when `sensitivity.signatures` names no signature (unknown
names are left out with a warning). Each second-review sensitivity module runs only when the step it re-examines ran
(positivity; proximity; uncoupling; compartments; leiden and domain_matrix), otherwise a log line says why; the
uncoupling module also returns nothing, with a warning, when `uncoupling.a` / `.b` have no calls. An empty
`sensitivity_analyses.lesion.genes` flags no bin and gives no warning; when its genes are found, the log says how many
bins are flagged, and a configuration without the osteoclast and osteoblast signatures also gets a warning. The Leiden
step warns when none of `cfg["ubiquitous"]` is in the data (the screen can mark no domain unresolved). The Figure 4h
warnings and the compartment DE table's description use `pc.compartment_names` and the signature's label
("<label> compartments: ..."), unchanged for the manuscript configuration.
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
Skipped, with the setting to change in the reason: Figure 3b when `figures.composition_signatures` names no signature
of the configuration (unknown names are left out), Figure 4a–f when a `figures.lineage_maps` signature is unknown or
has no call, a positivity panel that names no signature, Figure 6b when `figures.receptor_genes` names no gene of
`receptor_panel`. A sensitivity figure whose analysis step was switched off says so; the proximity, uncoupling and
positivity sensitivity figures are skipped when their main table has no rows, the Figure 4h one when
`sensitivity_Figure4h_top_genes.csv` has none (the Figure 4h test was not estimable: no gene to draw). A positivity
panel leaves out a signature that could not be scored (mean_score NaN in every section with bins: fewer than two of its
genes in the data) and names it in the caption; a panel of such signatures only is skipped. Figure 4h and its
sensitivity figure name the compartments with `pc.compartment_names` (growth plate / tumor-associated only for the
manuscript's split).
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
POST /api/projects/{pid}/validate  body {"config"?: dict (unsaved edits; default = saved config), "data"?: bool}
     -> report:
     {"ok": bool, "errors": [str], "warnings": [str],
      "samples": [{"id", "bin_dir": str|null, "bin_sizes": ["008um", "016um"], "has_positions": bool,
                   "has_scalefactors": bool, "has_hires_image": bool, "microns_per_pixel": float|null,
                   "n_bins": int|null, "message": str}],
      "genes": {"n_features": int|null,
                "signatures": {"<name>": {"present": [..], "missing": [..], "core_present": int}},
                "receptor_panel": {"missing": [..]}, "adipocyte_genes": {"missing": [..]},
                "ubiquitous": {"missing": [..]}, "lesion_genes": {"missing": [..]}}}
     n_bins / gene names come from filtered_feature_bc_matrix.h5 read with h5py (matrix/shape, matrix/features/name);
     fast, no scanpy. Warnings for the missing genes of each list; for the ubiquitous transcripts a stronger one when
     none is present (the Leiden screen marks no domain unresolved); for the A549 lesion genes
     (sensitivity_analyses.lesion.genes) only when none is present and steps.sensitivity_analyses is on.
     "data": false -> {"ok", "errors", "warnings"} of the configuration check alone (validate_config, as the project's
     "validation"; no data access): the project page checks its unsaved settings with it as they change. 400 when
     "data" is not a boolean.
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

Data explorer: the raw sections before any run. Every body may carry {"config": dict} (unsaved edits, exactly as
/validate; default = saved config); the analysis bin size of that config is read. A section is read as
bone_analysis.read_visium_hd reads it, with h5py + scipy.sparse in the server (no scanpy): the "Gene Expression"
features of filtered_feature_bc_matrix.h5 (names made unique like var_names_make_unique), barcodes without their -<n>
suffix (made unique), joined to the positions file by barcode, bins without coordinates dropped, matrix order kept.
POST /api/projects/{pid}/explore/samples
     -> {"bin_size", "samples": [{"id", "label", "condition", "color", "n_bins": int|null, "has_image": bool,
         "ok": bool, "message": str}]}   cheap: n_bins = in-tissue rows of the positions file (the matrix is not read);
        message = "1,234 bins at 016um" or why the section cannot be read (ok false also for more than 1,000,000 bins,
        the most a section the explorer and the gene maps read: 002um sections have millions)
POST /api/projects/{pid}/explore/{sid}/bins
     -> the bin payload of GET /api/runs/{rid}/bins (sample, n, x, y in µm = full-res pixel * microns_per_pixel rounded
        to 0.01 like the per-bin table's x_um/y_um, fields, "colors": {}) with fields total_umi, n_genes, pct_mt
        (number; pct_mt rounded to 0.01, null for a bin without counts) and passes_qc (bool: n_genes > qc.min_genes,
        total_umi > qc.min_counts, pct_mt < qc.max_pct_mt in float32 like bone_analysis.qc_metrics/qc_filter;
        mitochondrial = gene name starting with "mt-", any case), plus
        "barcodes": [..] (normalised: no -1 suffix),
        "summary": {"n_bins", "n_barcodes" (matrix), "n_without_position", "n_features", "total_umi", "median_umi",
                    "mean_umi", "median_genes", "median_pct_mt", "n_pass_qc", "pct_pass_qc", "bin_size",
                    "microns_per_pixel", "qc": {"min_genes", "min_counts", "max_pct_mt"}, "matrix": path,
                    "positions": path},
        "image": {"url": "/api/explore-images/<key>", "px_per_um", "width", "height"} | null   (hires H&E;
                 px_per_um = tissue_hires_scalef / microns_per_pixel: image pixel = x_um * px_per_um),
        "notes": [str]   (no image and why; a qc threshold that is not a number, replaced by its default)
     404 sample not in the config / no matrix / no positions; 422 unreadable files, microns_per_pixel unknown,
     fewer than half of the barcodes with a position (the analysis' own check), or a matrix of more than 1,000,000
     bins (refused from matrix/shape before its counts are read; also for /genes, /gene and the run gene maps)
POST /api/projects/{pid}/explore/{sid}/genes  body {"config"?, "query"?: str, "limit"?: 1..500 (default 50)}
     -> {"n_bins", "n_features", "n_matches", "genes": [{"gene", "total_umi", "n_bins_detected",
         "pct_bins_detected"}]}   by total_umi descending (ties: panel order); query = case-insensitive substring
POST /api/projects/{pid}/explore/{sid}/gene   body {"config"?, "genes": [1..50 names]}
     -> sparse gene vector {"genes_used", "genes_missing", "n", "indices", "values", "total", "n_detected", "max"}:
        per-bin sum of the raw counts of the genes, in the order of the bins payload; indices = bins with a count > 0.
        A name matches exactly, else the first gene of that name in any case; a gene named twice counts once.
GET  /api/explore-images/{key}  -> image/png   only keys the server issued while resolving an image (an HMAC of path,
                                   modification time and size; a client never sends a path); 404 otherwise
Caches (LRUCache: each key is built once, requests that arrive meanwhile wait for that build and get its value or
its error; an error is not kept): the counts of the sections last used (CSC, keyed by matrix + positions path, mtime,
size), at most 8 and 1 GB of their estimated memory (counts, per-bin arrays, about 300 bytes per barcode and gene; a
16 µm section of 105,000 bins is about 100 MB, an 8 µm one of 412,000 bins about 200 MB; the newest is kept whatever
its size); the encoded bins payloads, at most 8 and 128 MB (about 5 MB at 16 µm, 20 MB at 8 µm); 8 run-to-raw barcode
matches (one per run, sample, per-bin table and section).

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
     "removed" "#b8b8b8", "unassigned"/"ambiguous" "#dedede"). A positive_* column is yes/no only when its values
     are (booleans, text such as True/False, or numbers all 0 or 1: is_yes_no());
     positive_fraction_<slug>_pooled_depthmatched_draws (the share of depth-matching draws in which a bin is positive,
     0 to 1) is always a number, also with one draw (every share 0 or 1).
GET  /api/runs/{rid}/bins/image?sample=<sid>  -> {"url", "px_per_um", "width", "height"} of the hires H&E image under
     the bins of /bins (as in the explorer, with the run's config.json); 404 when there is none
GET  /api/runs/{rid}/genes?sample=<sid>&query=&limit=50  -> as POST /explore/{sid}/genes plus "n_unmatched", over the
     bins of the run's per-bin table (n_bins = its bins of that sample)
POST /api/runs/{rid}/gene  body {"sample", "genes": [1..50]}  -> sparse gene vector + "n_unmatched", aligned to the
     /bins order. The raw section is read again from the run's data_root; each per-bin barcode
     (<normalised matrix barcode>-<sample id>, as anndata.concat(index_unique="-") wrote it) is matched to its raw bin
     (fallbacks: as is, then normalised); unmatched bins count 0. 422 when the raw data can no longer be found.
GET  /api/runs/{rid}/export/loupe  -> application/zip, attachment "<rid>_loupe.zip": the per-bin calls as category CSV
     files for 10x Genomics Loupe Browser and README.txt (see "Loupe Browser export" below). 404 unknown run or no
     per-bin table (detail says when the run stopped before writing one); 409 while the run is queued or running; 422
     when a section's raw data can no longer be found or read, no bin of the table is a barcode of its matrix, a
     section of the table is not in the run's config, or the run's bins are smaller than 8 µm (no .cloupe exists).
     The UI: "Loupe Browser" in the run page's Export menu
GET  /api/runs/{rid}/export/results?part=all|figures|tables  -> application/zip, attachment "<rid>_results.zip",
     "<rid>_figures.zip" or "<rid>_tables.zip" (server.results_zip): files at their paths under results/. figures =
     figures/*.png and *.pdf; tables = every CSV/TSV (also gzipped); all = every result file except .h5ad objects (often
     gigabytes; downloaded on their own from the Tables tab) plus run_config.json and run_log.txt. PNG, PDF and .gz
     are stored, the rest deflated; built in a temporary file removed once sent. 400 unknown part; 404 unknown run or
     nothing of that part (detail says when the run stopped before writing any); 409 while queued or running. The UI:
     "All results", "Figures" and "Tables" in the run page's Export menu; "Download all" of the Figures tab (figures)

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

Loupe Browser export (`server.loupe_export(df, cfg, bin_size, out, title, source) -> report`, which GET /export/loupe
calls with the per-bin table, the run's config with an absolute data_root and the bin size of the table's file name):
* Zip entries, sections in config order: `<sid>_<B>_categories.csv` (the run's bin size, e.g. `016um`); when that is a
  whole multiple of 8 µm above 8 and the section has square_008um/filtered_feature_bc_matrix.h5, also
  `<sid>_008um_categories.csv` (Space Ranger writes a Visium HD section's .cloupe file, outs/cloupe_008um.cloupe, at 8 µm
  only; other sizes need `--custom-bin-size`); then README.txt (plain text: which file goes with which .cloupe file, rows,
  the parent check, the 8 µm quality control, import steps, the columns, the columns left out, the format choices; the
  Leiden "unassigned" size is the run's clustering.min_domain_bins). `<sid>` keeps letters, digits, `.`, `_`, `-`.
* CSV: comma-separated, LF line endings, UTF-8 without BOM, no quoting. Header `barcode`, then one column per category
  (its header is the name Loupe shows). One row per barcode of that bin size's filtered_feature_bc_matrix.h5, in matrix
  order and exactly as written there (-1 kept): every bin of the .cloupe file. Every cell holds a label: "not analyzed"
  in every column (but "8 um QC") when the run has no calls for the bin (per-bin barcodes are matched to the matrix as
  /gene matches them), "NA" for an analysed bin without a value in that column.
* Columns (loupe_columns), whatever the configuration: tissue_class "Tissue class", assigned_signature "Assigned
  signature", leiden_domain "Leiden domain" (a label of digits becomes "domain N"), in_region_of_interest "Region of
  interest" (inside/outside), then in table order every other yes/no column (positive_* positive/negative, in_*
  inside/outside, else yes/no) and text column of at most 500 values, headed by its name with spaces and a capital
  ("Positive osteoblast pooled depthmatched"; a name taken twice gets " 2"). Left out (listed in README.txt): numbers
  (total_umi, positive_fraction_*, lesion_umi_pooled), sample_id, section, x_um, y_um, empty columns, text columns of
  more than 500 values, and tissue_class, assigned_signature or leiden_domain beyond 32,766 values (Loupe Browser takes
  at most 32,768 per category, "NA" and "not analyzed" included). Labels: "," becomes ";", '"' becomes "'", tabs and
  line breaks become spaces.
* 8 µm rows take the labels of the run bin that holds them. Parent from the grid indices of the barcodes
  (s_008um_RRRRR_CCCCC): index * 8 // bin µm, as bone_analysis.map_domains_from_coarser. Checked against the bin centres
  of both tissue positions files (bins found there by grid key, not by barcode text): the nearest run-bin centre (all
  bins of the file, in tissue or not) must be the parent's. If any 8 µm centre is nearer another run bin, every 8 µm
  bin takes the run bin whose centre is nearest (README.txt says so). An 8 µm bin whose parent is not analysed is
  "not analyzed".
* The labels of an 8 µm bin are its parent's whatever its own counts, so 8 µm bins that would fail quality control on
  their own (even bins without a UMI) are labelled too, unlike the analysis' own projection (step project_fine_bins,
  which keeps the bins that pass). The 8 µm files therefore end with a column "8 um QC" (`matrix_qc`): "pass" or "fail"
  for every 8 µm bin, from its Gene Expression counts and the run's qc thresholds as `passes_qc` decides them (the
  explorer's passes_qc; float32 like bone_analysis), read 1,000,000 stored counts at a time; README.txt gives per file
  how many pass, how many labelled bins fail (and have no UMI) and how many "not analyzed" bins pass. If the counts
  cannot be read the file is written without the column, with a note. Manuscript data: 115,048 of 223,186 labelled
  LL-2-1 bins and 154,257 of 373,080 RL-1-5 bins fail (257 and 477 without a UMI), the same as a recount from the
  matrices.
* The zip is built in a temporary file (deleted after the response or on error), one CSV at a time; rows with the same
  labels share one text. Manuscript result (2 sections, 16 µm): 4.4 s, 7.3 MB zip (CSV 14 and 19 MB at 16 µm, 57 and
  78 MB at 8 µm); the process grows by about 0.45 GB at the peak (working set; measured), mostly while the parent check
  reads the 8 µm tissue positions file (702,000 rows with their barcodes for RL-1-5) and the 16 µm one, and by up to
  about 60 MB for the 8 µm quality control. Repeated exports do not accumulate memory.

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
button(), menuButton(), popover(), card(), field(), switchToggle(), tabs(), statusBadge(), loading(), emptyState(), errorPanel(),
banner(), toast(), modal(), confirmDialog(), promptDialog(), fmt.{int,num,pct,bytes,duration,date,ago}, debounce(),
sleep(), poll(), copyText(), clone(), state.health, refreshHealth(), route(), navigate(), setPathSilently(),
currentPath().
`menuButton({items, label, icon, kind, small, title, align, cls}) -> button`: a button (aria-haspopup="menu") that
opens a role="menu" popover (.menu-pop) of items {label, icon, hint, onClick, href, download, target, danger, disabled,
title} | 'separator' | {heading}; `items` may be a function, called at each opening. onClick may be async (a rejection
becomes a toast); an item with href is a link (download: true or a file name saves the file). Without a label it is
the ghost icon-only "More actions" button (aria-label = title), with one a button with a caret. One menu is open at a
time; it follows its button while the page scrolls and closes on a choice, Escape, Tab, a click elsewhere, a resize or
navigation, giving the focus back to its button; ↑/↓ on the button open it, ↑/↓/Home/End move in it; modal() leaves
the keys to a menu open in a dialog.
`popover(trigger, content, {label, align, focus, onClose, cls}) -> {el, place(), close({focus}), isOpen()}`: a small
panel that is not modal (role="dialog", .menu-pop.popover) under its button, for a short form (the spatial viewer's
"+ Add layer"). It shares the menus' slot (one menu or popover open at a time) and their placing (kept inside the
window when it fits on neither side of its button) and closing, closes on Escape (focus back on its button; a field
that uses Escape first, closing its own list, keeps it open), on Tab past its last field (focus on the page's next
stop after its button) or Shift+Tab before its first (focus on its button), and when the focus leaves it otherwise
(a click in another field); it focuses its first field on opening. The caller
sets aria-haspopup="dialog" on the button and toggles it with isOpen()/close().
Convention for secondary and destructive actions (move, remove, duplicate, delete, reset): a menuButton, or
`.row-actions`, a span in a table row or list item shown only while the row (tr, li, .has-row-actions) is hovered or
holds the focus (always on touch screens, and while its menu is open; `.keep` keeps it shown, e.g. while a toggle in
it is on). A row or card that opens something on click keeps a real link or button for the keyboard.
`css/app.css` (already written — read it) defines the design system: .page, .page-header, .card, .btn (primary,
danger, ghost, ghost.danger, small, icon-only, toggle[aria-pressed]), .menu-pop/.menu-item (menuButton), .menu-pop.popover (popover), .row-actions,
.input/.select/.textarea, .field, .form-grid, .switch, .tabs, .badge(-status), .chip(.ok/.warn/.bad/.core),
.table-wrap/.table, .kv, .metric, .progress, .log, .banner, .modal, .toast, .empty, .grid-2/.grid-3/.grid-auto, .row,
.stack, .legend, .prose. Use them; put view-specific rules in css/project.css
(projects, project editor, ROI editor, data explorer) or css/run.css (runs, results, viewers). The interface has one dark theme; its
colors come from the CSS variables on `:root` only — never hard-code colors except data colors (signature/sample
palettes).

Routes (main.js, already written) and the view modules they call — each export renders into `container` and may
return a cleanup function (stop timers, remove listeners):
```
#/projects                 views/projects.js  export async function renderProjects(container, {})
                           one card per project: a click on it opens the project's Data tab (its name is the link);
                           Duplicate and Delete project are in its "Project actions" menu (.row-actions). Header: Try
                           the demo (ghost) and New project, hidden while there is no project (the intro has them)
#/projects/:id/:tab?       views/project.js   export async function renderProject(container, {id, tab})
                           tabs: data | explore | parameters | signatures | roi | run
                           explore = dataExplorer() on the draft config (section "Data explorer" below)
                           header: the configuration-check chip only (see "Signatures and Parameters tabs"); Save
                           and Discard are in the bar at the bottom shown while there are unsaved changes (and
                           Ctrl+S). data: the sample table (a row's "Sample actions" menu: Move up, Move down, Remove
                           sample) and the Data check card ("Check again"). run: the project's runs (a click on a
                           row opens the run) and a "See all runs" link to #/runs
#/runs                     views/runs.js      export async function renderRuns(container, {})
                           no Refresh button: polls GET /api/runs every 4 s while a run is queued or running (not
                           while the page is hidden), refreshes when the window or tab is back and on the shell's
                           health events. A click on a row opens the run; the last column has Cancel (queued or
                           running) and the row's "Run actions" menu (Open project, Delete run), which a poll leaves
                           open
#/runs/:id/:tab?           views/run.js       export async function renderRun(container, {id, tab})
                           tabs: progress | summary | figures | tables | viewer | config  (tab null -> progress
                           while queued/running, summary when finished); viewer = spatialViewer(runSource(id)).
                           Header actions, built once (the menus read the run when they open): Cancel (only while
                           queued or running), Export (menu: All results, Figures, Tables = GET /export/results;
                           Loupe Browser = GET /export/loupe), Re-run, "Run actions" menu (Open project, Delete run;
                           Delete disabled while active). An export item is disabled while the run is active, while
                           that zip is being prepared, or when the manifest already fetched has nothing of its kind;
                           the zip is fetched (not linked) and saved, and a refusal shows the server's detail in a
                           toast. figures = figureGallery: a card has no buttons, a click opens the viewer (fit or
                           1:1, previous / next, Download PNG / PDF); "Download all" in the gallery's top line = the
                           figures zip. tables = tableBrowser; each data object and other file is one download link
                           (.h5ad objects are in no zip). config: Copy and Download are icon-only buttons beside
                           "Create project from this configuration"
#/help/:topic?             views/help.js      export async function renderHelp(container, {topic})
                           topics: guide (default) | changes | architecture | readme
```
Components live in js/components/ and are imported by the views.

### Signatures and Parameters tabs (components/signature-editor.js, components/param-form.js)
* `deleteSignature(cfg, name)` removes the signature from every list (competitors, primary, programs, roi.exclude,
  sensitivity sets, positivity panels, composition bars) and drops the proximity pairs that use it; single-signature
  settings stay. The delete dialog lists `signatureUsage(cfg, name)` -> `[{where, step}]` (the Parameters field that
  still selects it; step = the steps key that needs it, null for the Figure 4a–f maps), says that the analysis cannot
  start until another signature is chosen there or the step is switched off, and `deletionEffects(cfg, name)`: lists
  it empties (composition bars, a positivity panel, roi.exclude, sensitivity.signatures with the step on) and pairs it
  drops.
* A signature field shows, under it, why its value cannot be used (unknown, or without a positivity call when
  `needs_call`), with the same ways out, unless its `uses_step` is off; it is checked again when a step is switched.
  Empty composition and panel lists say that the figure is not drawn; the Tumor conditions chips are marked when no
  section has that condition. The Analysis steps card shows each step's help and a "bone study" badge from the
  schema's "steps" section (the Run tab shows the help as the step's tooltip, and the badge).
* Signatures tab buttons: Delete signature is in the "Signature actions" menu at the head of the signature, "Reset gene
  sets to manuscript defaults" in the "Gene set actions" menu of the Signatures card, beside Add (and in the empty
  state, beside Add signature); both ask for confirmation. The move and remove buttons of the calling-order lists and
  the tissue classes, and the remove buttons of excluded genes, receptor-panel rows, proximity pairs and positivity
  panels, are .row-actions.
* Parameters navigation: the section links on the left (hidden below 960 px) are grouped and the cards follow their
  order: Analysis (qc, clustering, scoring, roi), Figures (compartments, uncoupling, proximity, moran, the proximity
  pairs, figures, steps), Sensitivity analyses (key sensitivity_analyses[.*], or a title starting "Sensitivity
  analyses:"; linked without that prefix), Other (Reproducibility and any section added to the schema later). The
  link of the first card in the band from 16 px to a third of the window's height is current (aria-current=
  "location"); a click scrolls to the card, focuses its first visible control and keeps that link current until the
  next scroll, click or key.
* The configuration-check chip, the one control in the header of the project page, shows the check of the unsaved
  settings while there are any: POST /validate with "data": false, 300 ms after the last change (aborted when another
  follows), kept with the JSON of the config it checked; until it answers, and when it fails, it shows the saved
  settings' check and says so. The Run tab and starting a run use the saved settings' check, since runs use those.
  The chip ("Settings valid", or the counts of errors and warnings; a spinner while a data check runs) opens and
  closes the check panel below the header: the errors and warnings, and "Check the data again" (POST /validate with
  the data; aria-disabled while a check runs) with one line on the last data check.

### Spatial viewer (components/spatial-viewer.js, styles in css/run.css)
```
spatialViewer({source, store, binUm, emptyText}) -> {el, destroy(), show({panel?, sample?, field?, genes?, overlays?})}
runSource(runId) -> source        the Viewer tab of a run (GET /bins/samples, /bins, /genes, /bins/image; POST /gene)
geneInput({id, label, placeholder, submitLabel, search, onSubmit, onCancel, suggestOnFocus}) -> {el, input, set, focus, destroy}
                                  onSubmit(names, text typed); commas, spaces or + separate the names
fieldLabel(key) -> display name of a per-bin field
```
* source = {name?, emptyText?, samples(signal) -> [{id, label, n}], bins(sampleId, signal) -> bin payload (the shape
  of GET /api/runs/{rid}/bins, optional "barcodes"), genes?: {search(sampleId, query, signal) -> [{gene, total_umi,
  n_bins_detected, pct_bins_detected}], values(sampleId, genes, signal) -> sparse gene vector}, image?(sampleId,
  signal, {retry}) -> {url, px_per_um, width, height} | null (retry: asked again by the Retry of an image that failed
  to load)}. Every method gets an AbortSignal; a 404 from samples() shows the empty state. The run viewer uses
  runSource; the data explorer a source over the /explore endpoints.
* 1, 2 or 4 panels (2 x 2; one column below 760 px), each with its own section (a new panel shows the first section
  no panel before it shows) and two parts: "Color by", the base map, and layers on top.
  "Color by": "None (gray)" (every bin in the neutral color, under layers; its one legend entry "Every bin (gray)"
  hides them like any legend entry, which leaves the layers alone, e.g. only the bins passing QC over the H&E image;
  the legend's "x of y bins shown" counts the bins of the coloring and of the layers), the labels (category fields), the
  numbers (number fields), or "Gene…": the summed counts of one or more genes, 0 in the neutral color, detected bins
  on viridis with a log(1 + x) scale clamped at the 99th percentile of the nonzero counts; numbers likewise from their
  smallest value to their 99th percentile, log(1 + x) for UMIs and genes per bin; the color scale starts at viridis'
  quarter point, so the lowest step is brighter than the neutral color on the dark background. The yes/no fields are
  not in it: they are layers. A panel colored by one (a view kept from before, or viewer.show({field})) shows that
  field as a layer over no coloring, the same picture (its negatives hidden: the gray bins hidden).
  Layers (the panel's "overlays" in the store and the code): genes, each marking the bins detecting it, so that
  several genes show on one map, and boolean fields: calls, region of interest, quality control; each in the
  signature's configured color when free, else a palette color apart from the other layers and, while one is left,
  from the colors of the categories the panel is colored by; bins in two or more layers in the theme's text color, in
  the dark text color of the exported image over the H&E image and in that image; the coloring under them can be
  dimmed. A legend with counts and click-to-hide entries, and a title (section + what is shown: 'Tissue class + Ctsk,
  Acp5', or the layers alone over no coloring). Layers are added with "+ Add layer" (the panel's last control, under
  the label "Layers"), a popover() with a gene box (geneInput, suggestions only once typing; Enter or a suggestion
  adds the genes typed, one layer each, 'Ctsk, Acp5', or one for genes joined by +, 'Ctsk + Acp5' = the bins
  detecting any of them) and the boolean fields not shown yet (a click adds one). A layer asked for again (any case)
  flashes its chip instead; at most 11 per panel (one palette color each; viewer.show() adds no more either, with the
  same warning). Each layer is a chip under the controls (its color, its name, a remove button; marked when its genes
  are not in the data or failed to load), with "Dim the coloring under the layers" beside them. A gene layer typed in
  another case takes the data's names once its counts are there ('ctsk' -> 'Ctsk'). Adding, removing or renaming
  layers dispatches a bubbling "sv-change"
  event from the panel (the popover lives outside the viewer's element). "Color scales: Shared" (default on): panels colored by
  the same genes or the same number use one range (the smallest low end, the largest 99th percentile), the legend
  says whose range it shows. Pan and zoom are linked by default (one view in µm for all panels, with a crosshair at
  the hovered position in the others); unlinked, each panel keeps its own view per section. The H&E image
  (source.image) can be drawn under the bins, the bins then see-through by the chosen amount; it is asked for again
  with {retry: true} only by the Retry of an image that failed to load. Save PNG (a ghost button at the end of the
  toolbar) writes one composite image of all panels (title, map, scale bar, legend, with the notes of a shared or
  clamped color scale) on white. Zoom in, zoom out and fit are one compact group of icon-only ghost buttons
  (role="group", "Zoom") in the top-right corner of each map.
* Layout: the maps take the window's height less what must be seen with them (the viewer measures it: its place on
  the page when one row of panels starts in the top half of the first screen, else its toolbar and controls; measured
  again when the window is resized or the page above the viewer grows or shrinks; with four panels the second row is
  measured on its own, and a map is at least 200 px high); with
  two or four panels the maps of a row line up (the panels share the rows of a CSS subgrid), the controls of a panel
  sit on one line, and a wide panel (500 px with four panels, 700 px with two) has its legend in a 180 px column
  beside the map (an entry's count under its label), else under it. The tooltip of a bin (barcode, genes shown, the field colored by and the layers in
  bold, the other fields, the remaining positivity calls on one line) is placed in the window beside the pointer, or
  above or below it, never cut by the panel.
* store (kept by the caller between visits): layout, link, sameScale, view (linked), pointSize, image {on, opacity},
  panels [{sample, field, genes, overlays, dim, hidden, views}], and the caches samples, data (prepared bins per
  section, shared by the panels), geneData (the 40 gene vectors used last), images.
* Bins are drawn into one offscreen layer per panel only when the data, the zoom or the style change; panning copies
  the layers. Requests are shared by key across panels and aborted when no panel needs them any more.

### Data explorer (components/data-explorer.js, styles in css/project.css)
```
exploreState(cfg, key, previous) -> state        kept in the project's draft between visits
dataExplorer({projectId, name, state, getConfig, focus, onOpenData}) -> {el, destroy()}
```
* The Explore tab of a project: the raw sections before any run, read through the explorer endpoints (section 5) with
  the draft config, so unsaved data settings are explored as /validate checks them. project.js keys the state by
  dataKey(cfg) plus the qc thresholds (min_genes, min_counts, max_pct_mt); when that key changes, exploreState() starts
  a new state (nothing read; the viewer choices and the gene filter are kept). The state holds the config snapshot the
  requests send, the explore/samples answer, one entry per section ({status idle | loading | ok | error, summary, notes,
  image}), the spatialViewer store, the highlighted panel and the gene filter. "Read again" (a ghost button in the
  Sections card's header) re-reads everything.
* Sections table: one row per configured section, filled as each one is read (sections one at a time, the ones on the
  map first; a section that cannot be read shows why, a failed read has Retry): bins, genes in the matrix, UMIs,
  median UMIs, genes and % mitochondrial per bin, bins passing QC, microns per pixel; the qc thresholds used, barcodes
  without a position and the server's notes underneath. A click in the row of a readable section shows it in the
  highlighted panel (rowAction(): the section name is the row's button for the keyboard; a click that ends a text
  selection in the row does nothing, and the row's other controls, such as Retry, keep their own clicks).
* Map: spatialViewer on a projectSource (samples = the readable sections, bins = POST /explore/{sid}/bins, genes =
  /genes and /gene, image(sid, signal, {retry}) = the "image" of the section's bins payload, kept with the section;
  only retry (the viewer's Retry of an image that failed to load) posts /bins again, for a new image address). One
  bins request per section serves the table, the viewer and the image; payloads read ahead for the table are kept
  (at most 2) until the viewer takes them. Default: one panel, UMIs per bin, H&E available.
* Highlighted panel: the panel clicked or focused last, outlined when there are two or four. markTarget() keeps the
  outline, the section of the gene table, the active row of the sections table and the gene rows in step with it; it
  runs once the viewer is mounted (on every visit: with sections known from an earlier visit the viewer makes its
  panels at once), when panels are made or removed (a MutationObserver), after a press, key, focus or change in the
  viewer, on its "sv-change" event (layers added from "+ Add layer" or renamed; the highlighted panel stays) and after each viewer.show() from the
  explorer (section rows, gene rows, a section opened from the Data tab). Gene table: the section's genes by UMIs
  (POST /genes, case-insensitive filter, 50 / 200 / 500 rows) with UMIs, bins detecting the gene and their share. A
  click in a gene row (rowAction(), as in the sections table; the gene name is the row's button) =
  viewer.show({panel, sample, genes: [gene]}); the row is marked and its name pressed (aria-pressed) while the
  highlighted panel is colored by that gene alone. The + toggle at the end of the row (.btn.toggle, icon-only, in
  .row-actions) adds 'gene:<gene>' to the panel's layers (overlays), on top of its coloring, or removes it: pressed
  (aria-pressed, its icon a check mark) and kept visible (.keep) while the highlighted panel has it; its aria-label "Show <gene> as a
  layer" is the same in both states, its title says what a click does. The rows follow the coloring and the
  overlays of the highlighted panel however they change, updated in place (a focused button keeps the focus). Below
  760 px the count of bins detecting the gene is left out (the share stays, its title gives the count). A "Gene
  table" jump above the map (with four panels, and below 760 px where the panels stack) scrolls to the gene table and
  moves the focus to it.
* Data tab: the bin count of every section the data check found is a button (sampleEditor onExplore; aria-label
  "<bins>: explore section <id>") that opens the Explore tab with that section in the highlighted panel; the Data
  check card's pointer has an "Explore tab" text link (reportPanel onExplore) that opens the tab.
