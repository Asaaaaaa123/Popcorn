#!/usr/bin/env python3
"""
Compare a finished run of the fixed pipeline with the original manuscript code (_original/bone_analysis.py) on the
same data and in the same Python environment.

  1. Leiden domains, the unresolved fraction and the seed-stability ARI recomputed with the original code: identical
     to the run if (as intended) the fixes did not touch clustering. Differences from the numbers quoted in the
     manuscript then come from the computing environment (library versions), not from the fixes.
  2. Native positivity and the macrophage -> osteoclast proximity table with the original code: identical.
  3. Results the fixes change, before (original code) and after (fixed code):
       Figures 4g/6a right panels (depth-matched positivity), Figure 4h (osteoclast compartments),
       Figure 3a (pooled tissue assignment; reconstructed exactly from native scores, native thresholds with 20 null
       gene sets and pooled core-gene counts, which is what the original computed - its own pooling needs tens of GB).

Usage: python tools/compare_with_original.py --data data --results results/SR007144_manuscript_config
Writes <results>/comparison_with_original/ (CSV tables + summary.json + report.md).
"""
from __future__ import annotations

import os

for _v in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import importlib.util
import json
import sys
import time
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")


def _md(df: pd.DataFrame) -> str:
    """Markdown table without the optional tabulate dependency."""
    cells = [[str(c) for c in df.columns]] + [["" if pd.isna(v) else str(v) for v in row] for row in df.to_numpy()]
    lines = ["| " + " | ".join(cells[0]) + " |", "|" + "---|" * len(cells[0])]
    return "\n".join(lines + ["| " + " | ".join(r) + " |" for r in cells[1:]])


def load_original():
    spec = importlib.util.spec_from_file_location("original_bone_analysis", ROOT / "_original" / "bone_analysis.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="folder with SR007144_LL-2-1/ and SR007144_RL-1-5/")
    ap.add_argument("--results", required=True, help="results folder of the fixed pipeline (with joint_16um.h5ad)")
    args = ap.parse_args()
    orig = load_original()
    import bone_analysis as ba
    R = Path(args.results)
    out = R / "comparison_with_original"
    out.mkdir(parents=True, exist_ok=True)
    report, summary = [], {}
    t0 = time.time()

    def log(msg):
        print(f"[{time.time() - t0:7.1f} s] {msg}", flush=True)

    new = ad.read_h5ad(R / "joint_16um.h5ad")
    log(f"fixed run loaded: {new.n_obs:,} bins")

    # ---- 1. clustering with the original code ----------------------------------------------------------------
    adatas, _ = orig.load_sections(args.data, "016um")
    A = ad.concat(adatas, label="sample_key", index_unique="-", uns_merge="unique", join="outer")
    A.obs["sample"] = A.obs["sample"].astype("category")
    A = orig.normalize_and_cluster(A, resolution=0.3)
    log("original normalize_and_cluster done")
    de = orig.domain_markers(A)
    A = orig.collapse_fragments(orig.collapse_unresolved(A, orig.domain_screen(A, de)))
    same_bins = list(A.obs_names) == list(new.obs_names)
    same_leiden = same_bins and bool((A.obs["leiden"].astype(str).to_numpy() == new.obs["leiden"].astype(str).to_numpy()).all())
    unresolved = int((A.obs["leiden"].astype(str) == "unresolved").sum())
    ari = orig.clustering_stability(A)
    new_ari = pd.read_csv(R / "leiden_seed_stability_adjusted_rand_index.csv")
    summary["clustering"] = {"same_bins": same_bins, "identical_leiden_domains": same_leiden,
                             "unresolved_original_code": unresolved,
                             "unresolved_fixed_code": int((new.obs["leiden"].astype(str) == "unresolved").sum()),
                             "median_ari_original_code": round(float(ari["ARI"].median()), 4),
                             "median_ari_fixed_code": round(float(new_ari["ARI"].median()), 4)}
    log(f"clustering: {summary['clustering']}")

    # ---- 2. native calls and the macrophage -> osteoclast table with the original code ------------------------------
    A = orig.define_marrow_roi(orig.call_positive(orig.score_all(A)))
    native_same = {n: bool((A.obs[f"pos_{n}"].to_numpy(bool) == new.obs[f"pos_{n}"].to_numpy(bool)).all())
                   for n in orig.PRIMARY}
    roi_same = bool((A.obs["in_marrow_roi"].to_numpy() == new.obs["in_marrow_roi"].to_numpy()).all())
    A.obs["program"] = new.obs["program"].to_numpy()                            # only the macrophage table is compared
    _, mac_old = orig.proximity_tables(A)
    mac_new = pd.read_csv(R / "macrophage_to_osteoclast_proximity.csv")
    cols = list(mac_old.columns)
    summary["native"] = {"identical_positivity_calls": native_same, "identical_region_of_interest": roi_same,
                         "identical_macrophage_to_osteoclast_table": bool(
                             np.allclose(mac_old[cols[1:]].to_numpy(float), mac_new[cols[1:]].to_numpy(float)))}
    log(f"native: {summary['native']}")

    # ---- 3a. depth-matched positivity, before and after ------------------------------------------------------------
    before = orig.depth_matched(A, orig.PRIMARY, roi_key="in_marrow_roi")
    after = pd.read_csv(R / "Figure4g_6a_positivity_depth_matched.csv")
    dm = before.merge(after, on=["signature", "sample"], suffixes=("_original", "_fixed"))
    dm = dm[["signature", "sample", "pct_positive_mean_original", "pct_positive_mean_fixed", "pct_positive_sd_original",
             "pct_positive_sd_fixed", "mean_score_original", "mean_score_fixed"]]
    dm.to_csv(out / "Figure4g_6a_depth_matched_before_after.csv", index=False)
    log("depth-matched positivity compared")

    # ---- 3b. osteoclast compartments, before and after --------------------------------------------------------------
    A = orig.split_osteoclast_compartments(A)
    comp_before = A.obs.groupby(["sample", "oc_compartment"], observed=True).size().rename("n_bins").reset_index()
    comp_after = new.obs.groupby(["sample", "oc_compartment"], observed=True).size().rename("n_bins").reset_index()
    comp = comp_before.merge(comp_after, on=["sample", "oc_compartment"], how="outer", suffixes=("_original", "_fixed"))
    comp.to_csv(out / "Figure4h_compartment_sizes_before_after.csv", index=False)
    de_before = orig.compare_osteoclast_states(A)
    de_after = pd.read_csv(R / "Figure4h_osteoclast_compartment_DE.csv")
    rows = []
    for g in ("tumor_associated", "growth_plate"):
        b = de_before[de_before.group == g].nlargest(8, "scores")
        a_ = de_after[de_after.group == g].nlargest(8, "scores")
        rows.append({"group": g, "top8_original": ", ".join(b["names"]), "top8_fixed": ", ".join(a_["names"]),
                     "shared": len(set(b["names"]) & set(a_["names"]))})
    pd.DataFrame(rows).to_csv(out / "Figure4h_top_genes_before_after.csv", index=False)
    log("compartments compared")

    # ---- 3c. tissue assignment (Figure 3a) as the original code computed it -----------------------------------------
    for n in ba.PRIMARY:
        A.obs[f"orig_pooled_score_{n}"] = A.obs[f"score_{n}"]                   # the bug: native scores
    for n in ba.PROGRAMS:
        A.obs[f"orig_pooled_score_{n}"] = A.obs[f"score_{n}"]
        A.obs[f"orig_pooled_ncore_{n}"] = new.obs[f"pooled_ncore_{n}"].to_numpy()  # pooled core-gene counts
    A.uns["orig_pooled_thresholds"] = {n: orig.null_threshold(A, n, n_perm=20) for n in ba.PRIMARY}
    A = ba.assign_dominant_program(A, "orig_pooled_", "program_original_code")
    ct = (pd.crosstab(A.obs["sample"].astype(str), A.obs["program_original_code"].astype(str), normalize="index") * 100).round(2)
    ct_new = (pd.crosstab(new.obs["sample"].astype(str), new.obs["program"].astype(str), normalize="index") * 100).round(2)
    before_pct, after_pct = ct.stack(), ct_new.stack()
    before_pct.index = before_pct.index.set_names(["sample", "program"])
    after_pct.index = after_pct.index.set_names(["sample", "program"])
    fig3 = pd.concat({"pct_original": before_pct, "pct_fixed": after_pct}, axis=1).fillna(0)
    fig3.reset_index().to_csv(out / "Figure3a_assignment_before_after.csv", index=False)
    agree = float((A.obs["program_original_code"].astype(str).to_numpy() == new.obs["program"].astype(str).to_numpy()).mean())
    summary["figure3a_bins_with_same_assignment_pct"] = round(100 * agree, 2)
    log(f"figure 3a: {summary['figure3a_bins_with_same_assignment_pct']}% of bins keep their assignment")

    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    report.append("# Fixed pipeline vs original code (same data, same environment)\n")
    report.append(f"```json\n{json.dumps(summary, indent=2)}\n```\n")
    report.append("## Figures 4g/6a, depth-matched positivity (% of ROI bins)\n")
    report.append(_md(dm.round(4)))
    report.append("\n\n## Figure 3a, tissue assignment (% of bins)\n")
    report.append(_md(fig3.reset_index()))
    report.append("\n\n## Figure 4h, compartment sizes\n")
    report.append(_md(comp))
    report.append("\n\n## Figure 4h, top genes\n")
    report.append(_md(pd.DataFrame(rows)))
    (out / "report.md").write_text("\n".join(report), encoding="utf-8")
    log(f"done -> {out}")


if __name__ == "__main__":
    main()
