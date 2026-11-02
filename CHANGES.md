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

Three of these notes were wrong or misleading as first written; the second review (below) corrected them.

* **Depth matching caps both sections.** Bins above the shallower section's median UMI per bin are downsampled to it,
  bins below are kept. On SR007144 this lowers the control too (its ROI mean falls from 81.9 to 47.5 UMI); the residual
  advantage of the deeper section is small (+4.8%). The main effect is that depth-matched magnitudes depend on the scheme.
  `scoring.depth_matching = "quantile"` matches whole-section distributions, which leaves the tumor ROI about 15%
  shallower than the control ROI (conservative for tumor-enriched programs); it is not a correction of `cap`.
* **Proximity counts self-matches, and the baseline includes the targets.** A bin positive for both signatures is its own
  nearest target (0 µm), and the "all ROI bins" baseline includes the target bins. `proximity.exclude_self = true` only
  changes the target bins' own distances; it does not correct the baseline. Non-target and depth-matched baselines are in
  the sensitivity analyses.
* **Figure 4h does not mix sections on SR007144.** Every growth-plate bin is in the tumor-bearing section, so
  `compartments.de_sections = "tumor"` gives the same test. The confound that matters is within that section: the
  growth-plate bins are about 2.15 times deeper than the tumor-associated bins.
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

---

# Second scientific review (2026-11-02)

The fixed code was reviewed again, this time for the analysis itself rather than for coding defects. Nine reviewers each
took one area. Every finding was then checked twice, once against the code and once for its scientific impact, on the
SR007144 results. Of 34 findings, 33 were confirmed (most with a lower severity than first claimed) and 1 was not a
problem. 148 calculations were confirmed correct, including the data loading, units, QC, depth-matching subsampler,
scores and calls, Wilson intervals, distances and Moran's I; every reported number reproduces from the raw data.

**The headline directions hold under every check**: osteoclast, osteoblast and macrophage positivity far higher, and
granulocyte and erythroid positivity lower, in the tumor-bearing tibia; RANKL up. They hold within depth strata too. The
problems are in how some results are tested and described.

## What the review found

| # | Result | Finding (SR007144) | Severity | Done |
|---|---|---|---|---|
| F1 | Figure 4i "uncoupling" | The null places both osteoblast and osteoclast calls at random over the ROI. Pooled calls are clumped and both lineages sit on bone, so the ratio exceeds 1 whatever the osteoblasts do. The 1.73 / 2.88 split into an osteoblast-specific part of 0.40 / 0.44 (closer, not farther) and an osteoclast-pattern part of about 4.4 / 6.5. Under other nulls, the 1.66 between-limb ratio becomes 0.61–1.60. | major | sensitivity, caption |
| F7, F8, F26 | Figure 4h compartments | "Above the 90th percentile" collapses to "at least one chondrocyte-gene UMI": 44.5% of bins score exactly 0, so the percentile is 0. Only about 2 of the 116 "growth-plate" bins are near the growth plate; the rest sit at the lesion and periosteal margin. Acan and Comp define the split, so they are higher there by construction. | major | sensitivity, summary table, caption, flags |
| F9, F10, G1.2 | Figure 4h test | Growth-plate bins are 2.15 times deeper. P values treat neighbouring bins of one section as independent. 27% of the tested bins lie outside the ROI. Acp5 (and probably Ctsk, Mmp9) is robust; mt-Co2, mt-Cytb and Itm2b are not. | minor–major | sensitivity, caption |
| G1.1 | "macrophages near osteoclasts" (115/405, 28.4% vs 17.1%) | Calls sit in deep bins (median 271 vs 65 UMI). Against ROI bins of the same depth there is no enrichment: −3.3 pp, 95% block CI −8.9 to +3.1. | major | sensitivity |
| F11 | Fibroblast–osteoblast proximity | The baseline includes the target bins. Against non-target bins the tumor row is 1.14 times, not 0.92 times. The control's strong co-location holds. | minor | sensitivity, help text |
| F4, F5 | Region of interest | The control ROI keeps 25–45% skeletal muscle (tumor 8–17%), which inflates folds 1.2–1.5 times (osteoclast 171 → 113–139). For nerve and smooth muscle the 90th percentile is 0, so a single UMI of Acta2, Tagln, Myl9 or Pmp22 excludes a bin. This removes 27% of tumor osteoclast-positive bins, which is conservative. | minor | sensitivity |
| F3, F23 | Positivity rule | The 95% null threshold almost never decides; a call means "≥ 2 core genes detected, not out-competed". The competitor veto acts almost only in the control (native folds 100 → 171, 19 → 70). Osteocyte and fibroblast contrasts depend on min_core = 2. | minor | sensitivity |
| F6, F24, F25 | Depth matching | Depth-matched magnitudes depend on the scheme (RANKL 31 times with cap vs 23 times with distribution matching; the adipocyte-gene depletion roughly doubles under cap). The depth-matched SD covers the downsampling only. | minor | sensitivity, Wilson columns, help text |
| G2.1, G2.2, G4.2 | Figure 4i calls | The calls come from one depth-matching draw. The ratio of ratios spans 1.42–1.76 over draws, and about 27–30% of tumor per-bin calls change between draws. | minor | sensitivity |
| F12 | All bin-level intervals | One section per condition; Wilson intervals are 2–6 times too narrow for the tumor section. | limitation | captions, block-bootstrap intervals |
| F2 | Methods sentence | "Human A549 transcripts are not captured" is false. Mouse probes cross-capture the lesion (Krt7 in 2,523 vs 17 bins). It lies inside the ROI (about 17% of the tumor ROI) and makes up much of the tumor-section "unresolved" domain. Effect on rates: conservative. Itm2b in Figure 4h is probably tumor RNA. | major (text) | fixed, lesion flag, sensitivity |
| F13–F18, F20–F22, F19, G3 | Reporting | See the next section. Also: S5 names favour rare programs; domain 7 has no marker; Harmony on the section label confounds batch and condition; the three-gene-set check was never run in the delivered results. | minor | fixed or sensitivity |
| G4.1 | Pooling across the ROI edge | Not a problem: the borrowed neighbours are the same tissue. | — | — |

## Fixed (outputs that change)

* `Figure4g_6a_positivity_native_depth_roi.csv`: `fold_vs_reference` and `fold_tumor_over_control` come from the exact
  counts, not from rounded percentages (osteoclast 169.36 → 171.30; megakaryocyte `fold_vs_reference` was off by 26%).
* `Figure4g_6a_positivity_depth_matched.csv`: new `ci_lo` and `ci_hi` columns, a Wilson interval on the mean count.
* `Figure4h_osteoclast_compartment_DE.csv`: the scanpy `logfoldchanges` column is replaced. In 27% of rows it held
  ±19–31 values, because a gene undetected in one group divides by a 1e-9 floor. It is now `log2fc_mean_cp10k`, with
  `pct_group`, `pct_rest` and `selection_gene`, which marks the split signature's genes.
* New `Figure4h_compartment_summary.csv`, with the threshold actually used, the bins per compartment and section, the
  bins outside the ROI and the median UMI per compartment. A warning is raised when the threshold is 0.
* `morans_I_per_section.csv`: `pval_norm_fdr_bh` was NaN for every tumor-section gene on the squidpy path. It is now BH
  over the finite P values.
* `FigureS5_unresolved_bins.csv`: one row per section after the pooled row (control 18.3%, tumor-bearing 42.3%).
* `gene_set_sensitivity_three_gene_sets.csv`: new `genes_used` and `genes_missing` columns, and a warning for absent
  variant genes (Bglap).
* Captions:
  * S1: the 10-gene criterion was applied to the 16 µm bins. Retention is now given at 16 µm (78.1% / 90.4%), not read
    off the 8 µm curve.
  * 3a: the threshold-0 classes and the share of bins not drawn (48.6% / 29.0%).
  * 4g and 6a: what the Wilson intervals and the SD cover.
  * 4h: the effective rule, the bins per section, the depth difference, descriptive P values and the selection genes (†).
  * 4i: the null model, the random-placement-only interval and the single draw.
  * S5: the per-section unresolved shares.
* The module docstring, the `collapse_unresolved` note and the help glossary on A549 cross-capture. The parameter help
  texts for depth matching, `exclude_self`, the compartment percentile, `de_sections` and the ROI percentile.

## Sensitivity analyses added (`bone_sensitivity/`)

A new step, "Sensitivity analyses (second review)" (`steps.sensitivity_analyses`, on by default; settings in
`sensitivity_analyses`, editable on the Parameters page). It writes extra tables and figures and changes no manuscript
output; a test runs the analysis with the step on and off and compares every manuscript table. A module that fails
becomes a warning. The step first flags the cross-captured A549 lesion: per-bin `lesion_umi_pooled` and `in_lesion`
(Krt7 + Krt18 UMIs over a bin and its 8 neighbours ≥ 2).

| Module | Tables (`sensitivity_*.csv`) | Figure | Answers |
|---|---|---|---|
| uncoupling | `Figure4i_null_models`, `_summary`, `_decomposition`, `_call_stability` | `Sensitivity_Figure4i_null_models` | F1, G2: five null models (csr_both = manuscript, fixed_b, fixed_a, habitat_fixed_b, random_labels) × three call sets × five depth-matching draws; the decomposition of the ratio; per-bin call stability (also as per-bin `positive_fraction_*_pooled_depthmatched_draws`) |
| proximity | `proximity` | `Sensitivity_proximity_baselines` | G1.1, F11: non-target and depth-matched baselines (strata at the ROI and the source quantiles), 200 µm block-bootstrap intervals, the `program_dm` assignment |
| positivity | `positivity`, `positivity_folds`, `roi_composition`, `roi_exclusion_rules`, `positivity_rule_diagnostics`, `positivity_depth_targets` | `Sensitivity_positivity_folds` | F3–F6, F12, F23, F24: four ROI variants (manuscript, call-rule exclusion, assignment exclusion, without the lesion) × four rules (manuscript, min_core 1 and 3, no veto) × native / cap / quantile / ROI-quantile depth; Katz and spatial block-bootstrap fold intervals |
| compartments | `Figure4h_diagnostics`, `_regions`, `_variants`, `_de`, `_top_genes` | `Sensitivity_Figure4h_variants` | F7–F10, F26, G1.2: effective threshold; anatomical growth plate (osteoclast bins within 100 µm of a chondrocyte-dense region); ROI-only, depth-adjusted and depth-balanced spatial pseudobulk tests; drawn growth-plate polygons when configured |
| clustering | `FigureS5_unresolved_by_section`, `leiden_domain_annotation`, `leiden_harmony` | `Sensitivity_FigureS5_domains` | F2, F19–F22: unresolved share per section and ROI, overlap with the lesion, domain composition and flags, and the clustering without Harmony against the seed-to-seed ARI |

## Status

* Tested on the synthetic demo data: every module, its figure and the integration. 196 tests pass. One end-to-end
  comparison fails only when run after the API tests: for one zero-variance gene, a Moran's I P value is NaN in one run and
  0 in the other, which suggests the squidpy and built-in paths differ. It does not involve the sensitivity step. Not yet
  resolved.
* **Not yet run on SR007144.** The numbers quoted above come from the review's own checks on the existing results. The
  sensitivity tables for the real data will come from the next run (about 25 minutes).

## Recommended next changes

These are decisions for the authors. Each changes manuscript numbers or text, so none is applied by default. Run the
sensitivity analyses on SR007144 first and decide from their tables.

1. **Figure 4i: replace the null, or the claim.**
   * Report the distance against a null that keeps the osteoclast pattern (`fixed_b`, or `habitat_fixed_b` restricted to
     bone-lineage bins), averaged over depth-matching draws (mean and range), next to the manuscript's complete
     spatial randomness.
   * Withdraw "osteoblasts farther than chance" and "uncoupling in the tumor-bearing limb" unless the pattern-preserving
     nulls support them. On the review's checks they do not.
   * For a formal test, use a cross-type K or pair-correlation function with spatial-block inference.
2. **Figure 4h: define the growth plate by anatomy, or rename the contrast.**
   * Draw the physis on the H&E and set `sensitivity_analyses.compartments.growth_plate_masks`. The ROI editor would
     need a second region type for this; it is not built yet.
   * If only about 2 bins qualify, drop the growth-plate interpretation and describe the panel as osteoclast bins with
     vs without a chondrocyte-gene transcript at the lesion margin.
   * Test with the depth-adjusted or pseudobulk variant, remove the bin-level P values (or label them descriptive),
     drop the selection genes (Acan, Comp) from the bars, and annotate Itm2b as likely tumor RNA.
3. **Macrophage–osteoclast proximity.**
   * Report the depth-matched comparison.
   * Remove "enriched near osteoclasts" unless it survives: on the review's checks it does not.
   * Make the non-target baseline the default column of both proximity tables.
4. **Region of interest.**
   * Exclude skeletal muscle, smooth muscle and nerve by the pipeline's call rule (or a drawn ROI) instead of a 90th
     percentile, which is 0 for sparse scores.
   * Report the muscle share of each section's ROI and recompute the quoted folds and counts (405 macrophage-positive
     bins, 1,734 osteoclast-positive bins).
5. **Positivity: describe the rule as computed.**
   * Write it as "≥ 2 core genes co-detected, not out-competed; the score threshold rarely binds".
   * Add the min_core 1/2/3 table.
   * Quote folds whose control count is single-digit as lower bounds or with block-bootstrap intervals (osteoclast
     171× rests on 7 control bins, macrophage 70× on 4).
6. **Depth matching.** Keep `cap`, but report the ROI-quantile scheme next to it and do not quote precise depth-matched
   folds (for example "RANKL > 20-fold" rather than 31.07).
7. **Lesion.**
   * Delineate the A549 lesion (the `in_lesion` flag or a drawn region) and report the tumor-section rates with and
     without it.
   * Correct the Methods sentence on species specificity.
8. **Supplementary S5.**
   * Report the unresolved share per section.
   * Name domains by their top programs (`leiden_domain_annotation.csv`) and flag the domain without markers.
   * Disclose that Harmony on the section label cannot separate batch from condition.
9. **Statistics throughout.**
   * One animal and one section per condition: avoid inferential language about the tumor effect, and present bin-level
     intervals as within-section only.
   * More animals are the only remedy for population-level claims.
10. **Software.**
    * Run the analysis on SR007144 with the sensitivity step.
    * Resolve the Moran's I test-order issue.
    * Add growth-plate regions to the ROI editor.
    * If the authors adopt items 3–5, move the non-target baseline and the block-bootstrap intervals into the main tables.
