# What was fixed (2026-10-28)

The two scripts were audited line by line (five independent reviews, each finding checked by a second, adversarial
review against the source code of the exact library versions you use), then fixed and turned into a configurable
pipeline with a web interface. The untouched originals are in `_original/`. In the code, every change carries a
`[FIX]` note (a defect, with what it changed) or an `[APP]` note (new, for the configurable pipeline and the web app);
your `[SHAMIM]`, `[SHAMIM, MODIFIED]` and `[RUN]` provenance notes are kept as they were.

## What this means for the manuscript

**Numbers that change** (the original code did not compute what the Methods describe):

| Result | Why it changes |
|---|---|
| Figure 3a/3b, tissue assignment (`program`, `program_dm`) | pooled scores were the native per-bin scores (defect A1) |
| Figures 4g and 6a, right panels (depth-matched positivity) | depth-matched scores and thresholds were native (A1, A2); per-bin subsampling (A4) |
| Figure 4i (osteoblast-to-osteoclast distance) | uses the pooled, depth-matched calls (A1, A2, A4) |
| Fibroblast-to-osteoblast proximity ("41.8/41.3/67.3/5.6%") | uses the pooled tissue assignment (A1) |
| Supplementary Figure S5 domain names | named from the pooled tissue assignment (A1) |
| Figure 6b and Supplementary S4, depth-matched columns | per-bin subsampling (A4) |
| Figure 4h (osteoclast compartment DE) | control-limb osteoclasts are no longer in the "tumor-associated" group (A5) |
| Gene-set sensitivity, depth-matched columns (Methods) | scored from native data (A1) |

**Numbers that do not change**, verified bit for bit against the original code (tests/test_analysis_units.py,
`test_native_results_equal_the_original_pipeline`): QC and retention, native gene-set scores, null thresholds and
positivity calls (Figures 4g/6a left panels), the region of interest, Leiden domains and the unresolved fraction,
the seed-stability ARI, Moran's I (Supplementary S2), Table 1, the raw columns of Figure 6b/S4, and the
macrophage-to-osteoclast proximity ("115 of 405"), which uses native calls.

### Checked on the SR007144 data (2026-10-28)

The fixed pipeline was run on both sections with the manuscript configuration (`results/SR007144_manuscript_config/`,
21 min), and the original code was rerun on the same data in the same Python environment
(`tools/compare_with_original.py`).

| Number quoted in the code comments | Manuscript | Original code, this computer | Fixed code, this computer |
|---|---|---|---|
| 8 µm bins kept at > 100 genes (control / tumor-bearing) | 0.99% / 3.85% | — | 0.99% / 3.85% |
| Moran's I consensus genes (Suppl. S2) | 89 | — | 89 |
| Macrophage-positive bins in the tumor-bearing ROI | 405 | 405 | 405 |
| …within 50 µm of an osteoclast-positive bin | 115 (28.4% vs 17.1%) | 114 (28.1% vs 15.9%) | 114 (28.1% vs 15.9%) |
| Unresolved bins (Suppl. S5) | 50,273 (33.1%) | 50,527 (33.3%) | 50,527 (33.3%) |
| Median seed-stability ARI | 0.72 | 0.614 | 0.614 |

On the same computer the original and the fixed code give identical Leiden domains, ARI, native positivity calls,
region of interest and macrophage table, so the small differences from the manuscript values do not come from these
fixes. For the clustering numbers (unresolved bins, ARI) the likely cause is the computing environment (library
versions behind the neighbour graph and Leiden). The proximity numbers do not depend on clustering, so the quoted
"115 … 17.1%" was probably computed with an earlier version of the code (for example before the hypoxia signature was
removed on 2026-10-27). Rerunning in the manuscript's own environment and code version would settle it.

### How much the changed results move on the SR007144 data

Original code vs fixed code, same data and environment
(`results/SR007144_manuscript_config/comparison_with_original/`):

| Result | Original code | Fixed code |
|---|---|---|
| Figure 3a, ambiguous bins (control / tumor-bearing) | 23.6% / 53.3% | 10.9% / 43.6% |
| Figure 3a, tumor-bearing bins assigned osteoclast / osteoblast | 6.0% / 5.7% | 9.1% / 9.0% |
| Figure 3a, bins keeping the same assignment | | 83.5% |
| Figure 4g, depth-matched osteoclast-positive (control / tumor-bearing) | 0.0056% / 0.595% | 0.0067% / 0.582% |
| Figure 6a, depth-matched macrophage-positive, tumor-bearing | 0.081% | 0.072% |
| Figure 4h, compartment sizes | 7 control bins counted as tumor-associated | they are left out; top 8 genes of each compartment identical |

The tissue assignment changes most: pooling now really adds depth, so far fewer bins stay ambiguous. The depth-matched
positivity changes only slightly, and Figure 4h not materially.

## Defects in `bone_analysis.py`

| # | Severity | What was wrong | Fix |
|---|---|---|---|
| A1 | critical | `scanpy.tl.score_genes` reads `adata.raw` whenever it exists, and `normalize_and_cluster` set `a.raw = a` (native log-normalized). `AnnData.copy()` carries `.raw`, so in `depth_matched`, `spatially_pooled_scores` and `geneset_sensitivity` every "depth-matched" or "pooled" score was bit-for-bit the native score; only the core-gene count used the new matrix. | Scores read `layers["lognorm"]` explicitly (`layer=`, `use_raw=False`); `.raw` is no longer stored (the same values are in `layers["lognorm"]`). |
| A2 | major | `null_threshold` had the same problem, so every pooled/depth-matched null threshold was a native-depth threshold (and pooled = pooled depth-matched thresholds). | Same fix; the null is scored on a light object instead of a full copy. Random draws unchanged, so native thresholds are identical. |
| A3 | major | Spatial pooling assigned `A @ X[m]` into a `lil_matrix` by fancy row index, which densifies each section (tens of GB on 16 µm data); a section with fewer than 9 bins crashed the neighbor query. | One sparse neighbor matrix product (identical neighbors); sections smaller than k pool all their bins. |
| A4 | minor | `sc.pp.downsample_counts` re-seeds with the same seed for every bin, so all bins with the same UMI total kept the same UMI positions: one draw shared by a section, not an independent subsample per bin. | Own downsampler: an independent multivariate-hypergeometric draw per bin from one generator per replicate (same replicate seeds). |
| A5 | major | Every osteoclast-positive bin below the growth-plate threshold was "tumor_associated", including bins of the tumor-free contralateral limb, which then entered the Figure 4h test. | Only sections whose condition is a tumor condition hold tumor-associated bins; others become `non_growth_plate_control` and are not tested. |
| A6 | major | Competitor rule `score >= competitor`: when a competitor cannot be scored (NaN) every comparison is False, so the signature silently has no positive bin. | A competitor vetoes a bin only where it scores strictly higher. |
| A7 | major | A signature that could not be scored got no `pos_` column: copies kept a stale column from the previous run (the sensitivity analysis then reported the main calls for a variant whose genes were absent) and later steps crashed. | Always written (False) with a NaN threshold and a warning. |
| A8 | major | `_norm_bc` stripped `_<row>_<col>` from Visium HD barcodes lacking the `-1` suffix, collapsing every barcode to `s_016um`. | Only the `-<n>` suffix is removed (identical for your barcodes); a join that places < 50% of barcodes stops with an explanation. |
| A9 | major | `--roi-mask` imported `bone_pending_masks`, which did not exist (crash after the slow steps). | `bone_pending_masks.py` written: GeoJSON polygons (with holes) on the hires/lowres/full-res image → per-bin mask. Regions can be drawn in the web app. |
| A10 | major | Sample ids, folder names and microns per pixel hard-coded (`positivity_table` used `w["RL-1-5"] / w["LL-2-1"]`, `proximity_tables` looped over `SAMPLES` with `UM_PER_PX`, the receptor panel fold, `load_sections` used `SR007144_*`): any other dataset crashed, and a section with no ROI bins crashed too. | Everything from the configuration; microns per pixel read from each section's `scalefactors_json.json`. |
| A11 | minor | `define_marrow_roi`: an exclusion signature that cannot be scored made the percentile NaN and the ROI empty; bins with a NaN score were excluded. | Unscorable signatures are skipped with a warning; NaN never excludes. |
| A12 | minor | Groups with < 2 bins crash `rank_genes_groups` (domain markers; a compartment with 0–1 bins crashed or gave all-zero scores). | Untestable groups are left out with a warning; the compartment test needs ≥ 2 bins per side. |
| A13 | minor | `positivity_table` / `receptor_panel` / `geneset_sensitivity`: division by zero for a section without ROI bins. | NaN with a warning. |
| A14 | minor | Gene-set sensitivity used 20 null gene sets for its native calls (main analysis: 50), so its "Table 1" rows did not reproduce Figures 4g/6a. | Native calls use the main analysis' null size; the Table 1 rows now equal the main results exactly. |
| A15 | minor | Determinism: `PYTHONHASHSEED` set inside Python never applies; setting `NUMBA_NUM_THREADS` after numba started its threads makes numba raise on its next config reload (a notebook importing scanpy first). | Thread limits enforced at run time (`numba.set_num_threads`, `threadpoolctl`); the web app sets the variables before the worker starts. |
| A16 | minor | `clustering_stability` left `uns["_s0".."_s4"]` in the h5ad; `receptor_panel` added `expr_<gene>` columns to `obs`; `morans_i` copied the whole joint object. | Cleaned up / light objects. |
| A17 | minor | `map_domains_from_coarser` labelled bins without a 16 µm parent "unassigned", the label of fragment domains. | `no_parent_bin`. |
| A18 | minor | Docstring: score_genes' background is not "50 genes" but up to 50 per expression bin hit. | Reworded. |

## Defects in `bone_figures.py`

| # | Severity | What was wrong | Fix |
|---|---|---|---|
| F1 | major | `OUT` was a module global set only under `__main__`: importing the module and calling any figure function failed. | Every function takes the output folder; `make_figures()` is the entry point. |
| F2 | major | Supplementary S1 imported `pymupdf` (not a dependency) after writing the PDF, stopping the script before S2, S4 and S5; it also re-read the raw 8 µm matrices. | PNG written by matplotlib; S1 reads the per-bin QC table the analysis writes. |
| F3 | major | Figure 4i: a section that is not estimable has no observed distances → `xlim(0, NaN)` raised; an all-zero section gave identical limits; an empty file raised. | Handled; the panel says "not estimable". |
| F4 | major | Figures 4g/6a/6b/S4 raised `KeyError` for a gene or signature absent from the table (e.g. `Il11ra1`, not in the probe set). | Absent entries are shown as "not in data" instead of crashing. |
| F5 | minor | Figure 4h used fixed x-limits (-13.5, 10) that clip longer bars and their p-value labels; an empty or one-sided table crashed. | Limits follow the data. |
| F6 | minor | Figure 3b attached legend labels only for the first section, so a signature absent from the control lost its legend entry; a zero total divided by zero. | Legend built from the configured signatures. |
| F7 | minor | Error bars: Wilson bounds against rounded pseudo-counts gave tiny negative whiskers (matplotlib raises); depth-matched SD whiskers were centred on the clipped bar, not on the mean; count labels were rebuilt from rounded percentages. | Whiskers clipped at 0 and centred on the mean; labels from `n_positive_mean`. |
| F8 | minor | Sample ids, section names, microns per pixel and the Moran's I column names were hard-coded; the 1 mm scale bar did not suit small sections; the backed h5ad handle was never closed. | From the configuration; figures read only the result tables (no h5ad, no raw data). |

Each figure is drawn independently (one failure no longer stops the rest), and every figure is listed with a caption in
`figures/manifest.json`. Deliberate visual differences: Figure 4h's x range follows the data, Figure 4a-f uses
`figures.map_dpi` (250) instead of a fixed 300 dpi.

## Found by the final review of the new code (fixed)

A second, independent review of the finished code (four reviewers, every finding reproduced) found 19 problems,
all fixed and covered by tests:

* Web server: with the default local setup any website open in the same browser could use the API (DNS rebinding,
  cross-site requests): the server now only answers requests addressed to 127.0.0.1/localhost from its own pages and
  JSON bodies. On Windows a case-only rename of a sample (`ctrl-1` → `CTRL-1`) deleted its drawn region; ids that
  differ only in case are now refused. Ids with a trailing newline gave server errors. A worker could outlive a
  closed server (now killed with it on Windows and Linux).
* Analysis: settings that passed validation but crashed late (duplicate stability seeds, a QC bin size coarser than
  the analysis bin with label projection on, signatures without `competes`/`excluded`) are now caught or handled;
  `fold_tumor_over_control` is again tumor over control computed as in the original (it was the non-reference section
  over the reference, rounded twice).
* Figures: captions now describe the compartments and the depth matching as computed; bars use exact counts
  (a 1-of-73,140 bin was drawn at the axis floor); genes outside the analysed panel are no longer called "not in the
  probe set"; long gene lists no longer overlap.
* Web pages: unsaved drawings of a removed or renamed section could never be cleared; dialogs did not take keyboard
  focus; the spatial viewer's canvas was 2 px too large (slight stretch, hover offset); pooled (not depth-matched)
  calls were coloured and labelled as native ones.

## Methodological notes (defaults unchanged, options added)

* **Depth matching caps, it does not match.** Bins above the shallower section's median are downsampled to it and bins
  below are kept, so the deeper (tumor-bearing) section keeps more bins at the target depth: a residual depth advantage
  in every depth-matched comparison. `scoring.depth_matching = "quantile"` downsamples each bin of the deeper section to
  the shallower section's depth at the same quantile, so the per-bin depth distributions match. Worth a sensitivity check.
* **Proximity counts self-matches.** A bin positive for both signatures is its own nearest target (0 µm), and the
  "all ROI bins" baseline includes the target bins. `proximity.exclude_self = true` measures such bins to the nearest
  other target bin; the number of co-positive bins is now reported in every table.
* **Figure 4h mixes sections.** Growth-plate bins come from both limbs (two-fold capture difference);
  `compartments.de_sections = "tumor"` restricts the contrast to the tumor-bearing section.
* **Null size.** Native calls use 50 random gene sets, depth-matched and pooled calls 20 (`scoring.n_perm`,
  `scoring.n_perm_resampled`).

## New

* `pipeline_config.py`: one configuration (samples, parameters, gene signatures, steps, figures); its defaults are the
  manuscript setup, value for value.
* `python -m webapp`: the web application (projects, data validation, parameter and signature editors, region-of-interest
  drawing, background runs with live progress and log, summary, figures, tables, spatial viewer).
* Optional dependencies: without `harmonypy` the graph uses PCA; without `squidpy` Moran's I is computed by an identical
  built-in implementation (tested against squidpy).
* `tools/make_demo_data.py`: synthetic Visium HD tibiae for the demo and the tests; `tests/`: unit, regression,
  end-to-end and API tests.
