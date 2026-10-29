#!/usr/bin/env python3
"""
bone_figures.py
===============
Draws every code-generated panel of manuscript09102026.docx from the outputs of bone_analysis.py.
Files are named by the manuscript's figure numbering and written as PNG and PDF to <results>/figures/.

  Figure 1a        experimental design schematic
  Figure 2c        crops of the Space Ranger high-resolution H&E images (no drawing on the photographs)
  Figure 3a, 3b    deterministic tissue assignment ("denoise" version: skeletal muscle, nerve, smooth muscle and
                   endothelium not drawn; unassigned bins as a gray background; composition renormalized over the ten
                   retained signatures)
  Figure 4a-f      osteoclast- and osteoblast-positive bins within the region of interest
  Figure 4g, 6a    positivity, native depth and depth-matched, logarithmic axes
  Figure 4h        osteoclast compartment differential expression
  Figure 4i        osteoblast-to-osteoclast distance, observed against the random-placement null
  Figure 6b        host receptor panel, depth-matched
  Supplementary S1 quality control (8 um bins), S2 Moran's I, S4 adipocyte genes, S5 Leiden domain heatmap

Not drawn by code: Figures 1b, 1d, 2a-b. Figure 1c is a chart inside the Word file. Figures 2d, 2e, 5 and S3 are placeholders (see bone_pending_masks.py).

PROVENANCE: every function here is [RUN]. Shamim Mollah's original pipeline (00_shamim_original/) has plotting steps
(QC violin, PCA variance, UMAP, spatial cluster plot, marker plot), none of which is used in the manuscript.

Usage:  python Manuscript/code/bone_figures.py --root . --results Manuscript/code/results
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bone_analysis as ba                                                     

# [RUN] house style of every panel (8 pt type, embedded TrueType fonts in PDF)
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8, "axes.linewidth": 0.6,
                            "figure.dpi": 150, "savefig.bbox": "tight"})

S2 = ["LL-2-1", "RL-1-5"]
COL ={"LL-2-1": ba.PALETTE["control"], "RL-1-5": ba.PALETTE["tumor"]}
NAME = {"LL-2-1": "control", "RL-1-5": "tumor-bearing"}
UMPX = ba.UM_PER_PX
OUT = None


def lbl(c):
    return c.replace("_", " ").replace("neutrophil", "granulocyte")


def save(fig, name, dpi=200):
    fig.savefig(OUT / f"{name}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def section_xy_mm(obs_sample, xy, s):
    m = obs_sample == s
    return m, xy[m, 0] * UMPX[s] / 1000, xy[m, 1] * UMPX[s] / 1000


def scalebar(ax, xs, ys, fs=8, lw=3):
    x0 = xs.min() + .03 * (xs.max() - xs.min())
    y0 = ys.max() - .02 * (ys.max() - ys.min())
    ax.plot([x0, x0 + 1], [y0, y0], "k-", lw=lw, solid_capstyle="butt")
    ax.text(x0 + .5, y0 - .03 * (ys.max() - ys.min()), "1 mm", ha="center", va="bottom", fontsize=fs)


# ---------------------------------------------------------------------------- Figure 1a DO NOT NEED FOR PIPELINE SETUP
def figure1a():
    """[RUN] Timeline of the experiment."""
    fig, ax = plt.subplots(figsize=(8.5, 2.2))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.annotate("", xy=(.965, .62), xytext=(.035, .62), arrowprops=dict(arrowstyle="-|>", lw=2.5, color="k", mutation_scale=14))
    for x, day, lab in [(.06, "day −1", "µCT\n(baseline)"), (.26, "day 0", "A549-luc cells,\nintratibial, right tibia"),
                        (.455, "day 7", "BLI"), (.67, "day 20", "µCT"), (.885, "day 21", "BLI; endpoint")]:
        ax.plot([x, x], [.57, .67], "k-", lw=2.5)
        ax.text(x, .72, day, ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.text(x, .52, lab, ha="center", va="top", fontsize=8.5)
    ax.text(.5, .06, "Both tibiae of one animal profiled by Visium HD on a single capture slide", ha="center", va="center",
            fontsize=9, style="italic")
    save(fig, "Figure1a_experimental_design", dpi=250)


# ---------------------------------------------------------------------------- Figure 2c DO NOT NEED FOR PIPELINE SETUP
def figure2c(root):
    """[RUN] Plain crops of the Space Ranger high-resolution H&E images (pixel boxes in Manuscript/figures/README.md);
    arrows are to be placed by the pathologist, not drawn here."""
    from PIL import Image
    for s, box, name in [("RL-1-5", (1850, 1350, 2650, 1950), "Figure2c_top_tumor-bearing_tumor_bone_interface"),
                         ("LL-2-1", (700, 650, 1300, 1150), "Figure2c_bottom_control_chondro_osseous_junction")]:
        Image.open(Path(root) / f"SR007144_{s}/outs/binned_outputs/square_008um/spatial/tissue_hires_image.png") \
            .convert("RGB").crop(box).save(OUT / f"{name}.png")


# ---------------------------------------------------------------------------- Figure 3 CAN BE ADDED TO PIPELINE SETUP BUT THE GENE MARKERS SHOULD BE CHANGED BASED ON THE SITUATION
CLASS = {"osteoclast": "bone lineage", "osteoblast": "bone lineage", "osteocyte": "bone lineage", "chondrocyte": "chondrocyte",
         "neutrophil": "hematopoietic marrow", "erythroid": "hematopoietic marrow", "megakaryocyte": "hematopoietic marrow",
         "macrophage": "hematopoietic marrow", "fibroblast": "stroma", "adipocyte": "stroma"}
CCOL = {"bone lineage": ba.PALETTE["osteocyte"], "chondrocyte": ba.PALETTE["chondrocyte"],
        "hematopoietic marrow": ba.PALETTE["neutrophil"], "stroma": "#E6AB02"}
KEEP10 = ["osteoclast", "osteoblast", "osteocyte", "chondrocyte", "fibroblast", "macrophage", "neutrophil", "erythroid",
          "megakaryocyte", "adipocyte"]


def figure3a(a):
    """[RUN] Tissue class per bin (assignment 'program', native depth, pooled). Skeletal muscle, nerve, smooth muscle
    and endothelium are not drawn; ambiguous bins are a gray background without a legend entry."""
    smp = a.obs["sample"].astype(str).values
    xy = np.asarray(a.obsm["spatial"])
    prog = a.obs["program"].astype(str).values
    for s, name in zip(S2, ["Figure3a_top_control_tissue_classes_denoise", "Figure3a_bottom_tumor-bearing_tissue_classes_denoise"]):
        m, xs, ys = section_xy_mm(smp, xy, s)
        p = prog[m]
        w = 6.5
        fig, ax = plt.subplots(figsize=(w, w * (ys.max() - ys.min()) / (xs.max() - xs.min()) + .6))
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
        bg = p == "ambiguous"
        ax.scatter(xs[bg], ys[bg], s=.6, color="#e3e3e3", lw=0, rasterized=True)
        for c, col in CCOL.items():
            k = np.isin(p, [g for g, cc in CLASS.items() if cc == c])
            ax.scatter(xs[k], ys[k], s=1.4 if c in ("bone lineage", "chondrocyte") else .6, color=col, lw=0,
                       rasterized=True, label=c)
        x0 = xs.min() + .03 * (xs.max() - xs.min()); y0 = ys.max() - .02 * (ys.max() - ys.min())
        ax.plot([x0, x0 + 1], [y0, y0], "k-", lw=3)
        ax.text(x0 + .5, y0 - .025 * (ys.max() - ys.min()), "1 mm", ha="center", va="bottom", fontsize=8)
        ax.set_xlim(xs.min() - .05, xs.max() + .05); ax.set_ylim(ys.max() + .05, ys.min() - .05)
        ax.legend(markerscale=8, fontsize=7.5, ncol=4, loc="upper center", bbox_to_anchor=(.5, -.01), frameon=False)
        save(fig, name, dpi=250)


def figure3b(comp_dm):
    """[RUN] Composition after depth matching (assignment 'program_dm'), renormalized over the ten retained
    signatures; ambiguous, skeletal muscle, nerve, smooth muscle and endothelium excluded (counts written alongside)."""
    rows = []
    for s in S2:
        g = comp_dm[comp_dm["sample"] == s].set_index("program")["n_bins"]
        tot = int(g.reindex(KEEP10).fillna(0).sum())
        for p in KEEP10:
            n = int(g.get(p, 0))
            rows.append({"sample_id": s, "section": NAME[s], "signature": lbl(p), "bins_depth_matched": n,
                         "retained_bins": tot, "percent_of_retained": round(100 * n / tot, 2)})
    comp = pd.DataFrame(rows)
    comp.to_csv(OUT / "Figure3b_composition_denoise.csv", index=False)
    fig, ax = plt.subplots(figsize=(3.8, 3.6))
    for i, s in enumerate(S2):
        c = comp[comp["sample_id"] == s].set_index("signature")
        b0 = 0
        for p in KEEP10:
            v = c.loc[lbl(p), "percent_of_retained"]
            if v > 0:
                ax.bar(i, v, bottom=b0, color=ba.PALETTE[p], width=.6, label=lbl(p) if i == 0 else None)
                b0 += v
    ax.set_xticks([0, 1]); ax.set_xticklabels([NAME[s] for s in S2], fontsize=8)
    ax.set_ylabel("% of retained bins", fontsize=8); ax.set_ylim(0, 100)
    h, l = ax.get_legend_handles_labels()
    fig.legend(h[::-1], l[::-1], fontsize=6.5, loc="center left", bbox_to_anchor=(1.0, .5), frameon=False)
    fig.tight_layout()
    save(fig, "Figure3b_composition_denoise", dpi=250)


# ---------------------------------------------------------------------------- Figure 4a-f
def figure4af(a):
    """[RUN] Osteoclast-positive (blue), osteoblast-positive (red) and both (black) bins within the ROI; one spatial
    scale for all six panels."""
    obs = a.obs
    xy = np.asarray(a.obsm["spatial"])
    smp = obs["sample"].astype(str).values
    roi = obs["in_marrow_roi"].values.astype(bool)
    oc = obs["pos_osteoclast"].fillna(False).values.astype(bool) & roi
    ob = obs["pos_osteoblast"].fillna(False).values.astype(bool) & roi
    OUTC, IN, BLUE, RED = "#f2f2f2", "#cfcfcf", ba.PALETTE["osteoclast"], ba.PALETTE["osteoblast"]
    ext = {s: section_xy_mm(smp, xy, s)[1:] for s in S2}
    rng = {s: (ext[s][0].max() - ext[s][0].min(), ext[s][1].max() - ext[s][1].min()) for s in S2}
    labw, gap = .32, .06
    colw = (6.5 - labw - 2 * gap) / 3
    sc_ = colw / max(r[0] for r in rng.values())
    rowh = {s: rng[s][1] * sc_ for s in S2}
    head, vgap, foot = .25, .12, .45
    H = head + rowh["LL-2-1"] + vgap + rowh["RL-1-5"] + foot
    fig = plt.figure(figsize=(6.5, H))
    letters = iter("abcdef")
    y_top = H - head
    for s in S2:
        m, xs, ys = section_xy_mm(smp, xy, s)
        o_, b_, r_ = oc[m], ob[m], roi[m]
        y0 = y_top - rowh[s]
        for jx, col in enumerate(["osteoclast", "osteoblast", "overlay"]):
            ax = fig.add_axes([(labw + jx * (colw + gap)) / 6.5, y0 / H, rng[s][0] * sc_ / 6.5, rowh[s] / H])
            ax.set_xlim(xs.min(), xs.max()); ax.set_ylim(ys.max(), ys.min()); ax.axis("off")
            ax.scatter(xs[~r_], ys[~r_], s=.08, color=OUTC, lw=0, rasterized=True)
            ax.scatter(xs[r_], ys[r_], s=.08, color=IN, lw=0, rasterized=True)
            if col in ("osteoclast", "overlay"):
                k = o_ & ~(b_ & (col == "overlay"))
                ax.scatter(xs[k], ys[k], s=1.2, color=BLUE, lw=0, rasterized=True)
            if col in ("osteoblast", "overlay"):
                k = b_ & ~(o_ & (col == "overlay"))
                ax.scatter(xs[k], ys[k], s=1.2, color=RED, lw=0, rasterized=True)
            if col == "overlay":
                ax.scatter(xs[o_ & b_], ys[o_ & b_], s=1.6, color="black", lw=0, rasterized=True)
            scalebar(ax, xs, ys, fs=6, lw=2)
            ax.text(-.02, 1.0, next(letters), transform=ax.transAxes, fontsize=10, fontweight="bold", va="top", ha="right")
            if s == "LL-2-1":
                fig.text((labw + jx * (colw + gap) + colw / 2) / 6.5, (H - head / 2) / H, col, ha="center", va="center", fontsize=9)
        fig.text(.12 / 6.5, (y0 + rowh[s] / 2) / H, NAME[s], rotation=90, ha="center", va="center", fontsize=9)
        y_top = y0 - vgap
    hs = [Line2D([], [], marker="o", ls="", color=c, ms=5) for c in (BLUE, RED, "black", IN, OUTC)]
    fig.legend(hs, ["osteoclast-positive", "osteoblast-positive", "positive for both", "other unit in region of interest",
                    "outside region of interest"], loc="lower center", ncol=3, fontsize=7, frameon=False, bbox_to_anchor=(.5, 0))
    fig.savefig(OUT / "Figure4a-f_bone_lineage_maps.png", dpi=300)
    fig.savefig(OUT / "Figure4a-f_bone_lineage_maps.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------- Figures 4g, 6a CAN BE ADDED TO PIPELINE SETUP
def positivity_panel(t1, t2, sigs, name, xt):
    """[RUN] Left: positive bins at native depth with Wilson intervals. Right: depth-matched mean +/- SD over five
    replicates. Values in percent; anything below 0.001% (including zero) sits at the axis floor and shows no bar."""
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.2), sharey=True)
    for j, s in enumerate(S2):
        a1 = t1[t1["sample"] == s].set_index("signature").loc[sigs]
        a2 = t2[t2["sample"] == s].set_index("signature").loc[sigs]
        x = np.arange(len(sigs)) + (j - .5) * .38
        y1 = a1.pct_positive.clip(lower=1e-3)
        axes[0].bar(x, y1, .38, color=COL[s], label=NAME[s], yerr=[y1 - a1.ci_lo.clip(lower=5e-4), a1.ci_hi - y1], capsize=2)
        for xi, k_, hi_ in zip(x, a1.n_positive, a1.ci_hi):
            axes[0].text(xi, hi_ * 1.25, f"{k_:,}", ha="center", fontsize=6)
        mdm = a2.pct_positive_mean / 100 * a1.n_bins
        y2 = a2.pct_positive_mean.clip(lower=1e-3)
        axes[1].bar(x, y2, .38, color=COL[s], label=NAME[s], yerr=a2.pct_positive_sd, capsize=2)
        for xi, k_, v_ in zip(x, mdm, y2 + a2.pct_positive_sd):
            axes[1].text(xi, v_ * 1.25, f"{k_:,.1f}" if k_ < 100 else f"{k_:,.0f}", ha="center", fontsize=6)
    for ax, yl in zip(axes, ["% positive units, native depth (log)", "% positive units, depth-matched (log)"]):
        ax.set_yscale("log"); ax.set_xticks(range(len(sigs)))
        ax.set_xticklabels(xt, fontsize=7.5, rotation=30 if len(sigs) > 3 else 0, ha="right" if len(sigs) > 3 else "center")
        ax.set_ylabel(yl, fontsize=7.5); ax.set_ylim(1e-3, 150)
    h_, l_ = axes[0].get_legend_handles_labels()
    fig.legend(h_, l_, fontsize=7, frameon=False, loc="lower center", ncol=2, bbox_to_anchor=(.5, 1.0))
    fig.tight_layout()
    save(fig, name)


# ---------------------------------------------------------------------------- Figure 4h CAN BE ADDED TO PIPELINE SETUP 
def figure4h(t4):
    """[RUN] Top eight genes per osteoclast compartment by Wilcoxon score; gray bars did not reach adjusted P < 0.05."""
    gp = t4[t4.group == "growth_plate"].nlargest(8, "scores")
    ta = t4[t4.group == "tumor_associated"].nlargest(8, "scores")
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    yl = []
    for k, (_, r) in enumerate(pd.concat([ta, gp.iloc[::-1]]).iterrows()):
        v = r.scores if r.group == "tumor_associated" else -r.scores
        col = (ba.PALETTE["osteoclast"] if r.group == "tumor_associated" else ba.PALETTE["chondrocyte"]) if r.pvals_adj < .05 else "#bdbdbd"
        ax.barh(k, v, color=col)
        yl.append(r.names)
        ax.text(v + (.2 if v > 0 else -.2), k, f"{r.pvals_adj:.1e}", va="center", ha="left" if v > 0 else "right", fontsize=6)
    ax.set_yticks(range(len(yl))); ax.set_yticklabels(yl, fontsize=7.5, style="italic")
    ax.axvline(0, color="k", lw=.6); ax.set_xlim(-13.5, 10)
    ax.set_xlabel("Wilcoxon score: growth-plate compartment (left) / tumor-associated (right)", fontsize=7.5)
    save(fig, "Figure4h_osteoclast_compartment_DE")


# ---------------------------------------------------------------------------- Figure 4i CAN BE ADDED TO PIPELINE SETUP
def figure4i(dist):
    """[RUN] Cumulative distributions of osteoblast-to-nearest-osteoclast distance, observed (solid) against the
    random-placement null (dashed), per limb."""
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.0), sharey=True)
    for ax, s in zip(axes, S2):
        ds = dist[dist["sample"] == s]
        for src, sty in [("observed", "-"), ("null", "--")]:
            dd = ds[ds.source == src].sort_values("distance_um")
            ax.step(dd.distance_um, dd["count"].cumsum() / dd["count"].sum(), sty, where="post", color=COL[s],
                    label=f"{NAME[s]}, {src}")
        w = ds[ds.source == "observed"].sort_values("distance_um")
        cw = w["count"].cumsum() / w["count"].sum()
        ax.set_xlim(0, float(w.distance_um[cw >= .99].min()) * 1.1)
        ax.set_xlabel("distance to nearest osteoclast-positive unit (µm)", fontsize=7.5)
        ax.legend(fontsize=7, frameon=False, loc="lower right")
    axes[0].set_ylabel("cumulative fraction of osteoblast-positive units", fontsize=7.5)
    fig.tight_layout()
    save(fig, "Figure4i_osteoblast_osteoclast_distance")


# ---------------------------------------------------------------------------- Figure 6b, Supplementary S4 CAN BE ADDED TO PIPELINE SETUP
def detection_panel(t, genes, name, shade_first_two=False, offset=.02, comma=False):
    """[RUN] Percentage of ROI bins detecting each gene after depth matching, Wilson 95% intervals."""
    fig, ax = plt.subplots(figsize=(6.5, 3.2 if shade_first_two else 3.0))
    for j, s in enumerate(S2):
        sub = t[t["sample"] == s].set_index("gene").loc[genes]
        x = np.arange(len(genes)) + (j - .5) * .38
        ax.bar(x, sub.pct_detected_dm, .38, yerr=[sub.pct_detected_dm - sub.ci_lo_dm, sub.ci_hi_dm - sub.pct_detected_dm],
               color=COL[s], label=NAME[s], capsize=2)
        for xi, k_, v_ in zip(x, sub.n_positive_dm, sub.ci_hi_dm):
            ax.text(xi, v_ + offset, f"{int(k_):,}" if comma else str(k_), fontsize=6, ha="center")
    if shade_first_two:
        ax.axvspan(-.5, 1.5, color="0.93", zorder=-1)                          # RANKL/OPG pair
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, fontsize=8, style="italic")
    ax.set_ylabel("% units detecting (depth-matched)", fontsize=8)
    ax.legend(fontsize=7, frameon=False)
    save(fig, name)


# ---------------------------------------------------------------------------- Supplementary S1 CAN BE ADDED TO PIPELINE SETUP
def figureS1(root):
    """[RUN] Per-bin QC distributions at 8 um and retention against the minimum gene count (solid line: the 10-gene
    criterion used; dashed: 100 genes). Reads the Space Ranger 8 um matrices through bone_analysis.read_visium_hd."""
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.2))
    unit = " per 8 µm bin"
    for s in S2:
        a = ba.qc_metrics(ba.read_visium_hd(Path(root) / f"SR007144_{s}/outs/binned_outputs/square_008um", s))
        sw = ba.qc_sweep(a, s)
        c = COL[s]
        ax[0].hist(np.log10(a.obs["n_genes_by_counts"] + 1), bins=60, alpha=.55, label=NAME[s], color=c, density=True)
        ax[1].hist(np.log10(a.obs["total_counts"] + 1), bins=60, alpha=.55, label=NAME[s], color=c, density=True)
        ax[2].hist(a.obs["pct_counts_mt"], bins=60, alpha=.55, label=NAME[s], color=c, density=True)
        sw = sw[sw["min_genes"] > 0]
        ax[3].plot(sw["min_genes"], sw["pct_retained"], "o-", color=c, label=NAME[s])
    ax[0].set_xlabel("log10 genes" + unit); ax[1].set_xlabel("log10 UMIs" + unit); ax[2].set_xlabel("% mitochondrial counts")
    for x in ax[:3]:
        x.set_ylabel("density")
    ax[3].axvline(100, ls="--", c="k", lw=.8); ax[3].axvline(ba.MIN_GENES, ls="-", c="0.4", lw=.8)
    ax[3].set_xscale("log"); ax[3].set_xlabel("minimum genes" + unit + " (log)"); ax[3].set_ylabel("% bins retained")
    for x in ax:
        x.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "FigureS1_quality_control.pdf")
    plt.close(fig)
    import pymupdf                                   # the PNG is a 220 dpi rendering of the PDF
    pymupdf.open(OUT / "FigureS1_quality_control.pdf")[0].get_pixmap(dpi=220).save(OUT / "FigureS1_quality_control.png")


# ---------------------------------------------------------------------------- Supplementary S2 CAN BE ADDED TO PIPELINE SETUP
def figureS2(cons):
    """[RUN] Moran's I per section for the top 30 (by mean rank) of the genes in the top 200 of both sections."""
    c6 = cons.sort_values("mean_rank").head(30)
    fig, ax = plt.subplots(figsize=(4.6, 6.2))
    y = np.arange(len(c6))[::-1]
    ax.barh(y + .2, c6["moranI_LL-2-1"], .4, color=COL["LL-2-1"], label="control")
    ax.barh(y - .2, c6["moranI_RL-1-5"], .4, color=COL["RL-1-5"], label="tumor-bearing")
    ax.set_yticks(y); ax.set_yticklabels(c6.gene, fontsize=6.5, style="italic")
    ax.set_xlabel("Moran's I (per section)", fontsize=8); ax.legend(fontsize=7, frameon=False)
    save(fig, "FigureS2_spatially_variable_genes")


# ---------------------------------------------------------------------------- Supplementary S5 CAN BE ADDED TO PIPELINE SETUP
def figureS5(M):
    """[RUN] Mean signature score per Leiden domain, z-scored across domains (color limits -2 to 2)."""
    M = M.set_index("leiden_domain")
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    im_ = ax.imshow(M.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index, fontsize=7)
    ax.set_xlabel("signature score", fontsize=8); ax.set_ylabel("Leiden domain (named by most enriched signature)", fontsize=8)
    cb = fig.colorbar(im_, ax=ax, shrink=.7)
    cb.set_label("z-scored mean across domains", fontsize=7)
    fig.tight_layout()
    save(fig, "FigureS5_leiden_domain_heatmap")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--results", default="Manuscript/code/results", help="--out of bone_analysis.py")
    args = ap.parse_args()
    R = Path(args.results)
    OUT = R / "figures"
    OUT.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(R / "joint_16um.h5ad", backed="r")
    nat = pd.read_csv(R / "Figure4g_6a_positivity_native_depth_roi.csv")
    dm = pd.read_csv(R / "Figure4g_6a_positivity_depth_matched.csv")

    figure1a()
    figure2c(args.root)
    figure3a(a)
    figure3b(pd.read_csv(R / "Figure3b_program_composition_depth_matched.csv"))
    figure4af(a)
    positivity_panel(nat, dm, ["osteoclast", "osteoblast", "osteocyte"], "Figure4g_bone_lineage_positivity",
                     ["osteoclast", "osteoblast", "osteocyte"])
    figure4h(pd.read_csv(R / "Figure4h_osteoclast_compartment_DE.csv"))
    figure4i(pd.read_csv(R / "Figure4i_distance_distributions_observed_and_null.csv", keep_default_na=False))   # "null" is a label, not NA
    positivity_panel(nat, dm, ["macrophage", "fibroblast", "neutrophil", "erythroid", "megakaryocyte"],
                     "Figure6a_myeloid_positivity", ["macrophage", "fibroblast", "granulocyte", "erythroid", "megakaryocyte"])
    detection_panel(pd.read_csv(R / "Figure6b_receptor_panel_depth_matched.csv"),
                    ["Tnfsf11", "Tnfrsf11b", "Tnfrsf11a", "Pth1r", "Csf1", "Csf1r", "Serpine1"], "Figure6b_receptor_panel",
                    shade_first_two=True)
    figureS1(args.root)
    if (R / "FigureS2_spatially_variable_genes_morans_I.csv").exists():
        figureS2(pd.read_csv(R / "FigureS2_spatially_variable_genes_morans_I.csv"))
    detection_panel(pd.read_csv(R / "FigureS4_adipocyte_genes_depth_matched.csv"), ba.ADIPOCYTE_GENES,
                    "FigureS4_adipocyte_genes", offset=.04, comma=True)
    figureS5(pd.read_csv(R / "FigureS5_domain_by_signature_zscores.csv"))
    print(f"done -> {OUT}")
