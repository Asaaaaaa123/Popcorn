"""
bone_sensitivity
================
Sensitivity analyses added after the second scientific review of the pipeline (2026-11-02; CHANGES.md). Each module
re-examines one manuscript result under the alternatives the review asked for, on the joint AnnData object that
bone_analysis.run_analysis_from_config builds, and writes its own tables and figures. None of them changes a
manuscript output: the defaults of pipeline_config.py remain the manuscript analysis.

  uncoupling     Figure 4i: null models that keep the spatial pattern of the calls, call sets, depth-matching draws
  proximity      "within 50 um" tables: non-target and depth-matched baselines, spatial block intervals
  positivity     Figures 4g/6a: regions of interest, calling rules, depth-matching schemes, block intervals, lesion
  compartments   Figure 4h: effective threshold, anatomical definition, depth-adjusted and spatially aggregated tests
  clustering     Supplementary S5: domains per section, domain annotation, clustering without Harmony

Every compute module exposes run(a, cfg, progress) -> list[Output]; every <module>_figures.py exposes FIGURES and
draw(results, cfg, out) and needs no scanpy (bone_figures imports it). Shared helpers are in common.py.
"""
