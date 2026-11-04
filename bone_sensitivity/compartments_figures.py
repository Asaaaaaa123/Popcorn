"""
Figure of the Figure 4h sensitivity analysis (compartments.py): the genes drawn in Figure 4h (top eight per
compartment by manuscript score) against every variant of the compartment test, as a dot matrix. Filled = BH-adjusted
P < 0.05 in that variant, colour = the compartment the gene is higher in, open = not significant; genes of the split_by
extended list (which select the growth-plate bins) are marked. Reads only the result tables (no scanpy).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import pipeline_config as pc
from bone_sensitivity import common

NAME = "Sensitivity_Figure4h_variants"
TITLE = "Sensitivity — Figure 4h under other definitions and tests"
FIGURES = [(NAME, TITLE)]
DIAG_CSV = "sensitivity_Figure4h_diagnostics.csv"
REGIONS_CSV = "sensitivity_Figure4h_regions.csv"
VARIANTS_CSV = "sensitivity_Figure4h_variants.csv"
TOP_CSV = "sensitivity_Figure4h_top_genes.csv"
INPUTS = [VARIANTS_CSV, TOP_CSV, DIAG_CSV, REGIONS_CSV]
GROUPS = ("growth_plate", "tumor_associated")
ALPHA = 0.05
LABELS = {"manuscript": "Manuscript", "roi_only": "ROI only", "depth_adjusted": "Depth-\nadjusted",
          "block_pseudobulk": "Block\npseudobulk", "anatomical": "Anatomical\ngrowth plate",
          "drawn": "Drawn\ngrowth plate"}
GRAY, LIGHT = "#9a9a9a", "#ececec"


def _names(cfg):
    """[APP] (adjective, noun, other compartment) as the text calls the two compartments: the manuscript's
    growth-plate / growth plate / tumor-associated for its osteoclast / chondrocyte split, else the names of
    pipeline_config.compartment_names (high-<split_by> / low-<split_by>)."""
    gp, ta = pc.compartment_names(cfg)
    return gp, ("growth plate" if gp == "growth-plate" else gp), ta


def _labels(cfg):
    """Column labels of the variants (the anatomical and drawn variants named after the compartment)."""
    noun = _names(cfg)[1]
    return {**LABELS, "anatomical": f"Anatomical\n{noun}", "drawn": f"Drawn\n{noun}"}


def _color(cfg, name, fallback):
    return (cfg.get("signatures", {}).get(name) or {}).get("color") or fallback


def _num(x, fmt="{:,.0f}", missing="n/a"):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return missing
    return fmt.format(v) if np.isfinite(v) else missing


def _label(variant, labels=LABELS):
    """Column label of a variant on one line."""
    return labels.get(variant, variant).replace("-\n", "-").replace("\n", " ")


def _units_short(units):
    return "blocks" if str(units).startswith("block") else "bins"


def _what(variant, S, P, st, n_reps, GPn="growth plate"):
    return {"manuscript": "the bins and test of Figure 4h",
            "roi_only": "both groups inside the marrow region of interest",
            "depth_adjusted": f"the deeper group downsampled bin by bin to the other's UMI distribution, median of "
                              f"{n_reps} replicates",
            "block_pseudobulk": f"units = {float(st['pseudobulk_block_um']):g} um blocks x compartment per section "
                                f"(summed counts), against pseudo-replication of contiguous bins; the deeper "
                                f"compartment's units downsampled to the other's unit-UMI distribution (a dense "
                                f"compartment sums more bins per block), median of {n_reps} replicates",
            "anatomical": f"{GPn} = {S}-positive bins within {float(st['anatomical_distance_um']):g} um of a "
                          f"{P}-dense region, both compartments from tumor sections",
            "drawn": f"{GPn} = {S}-positive bins inside hand-drawn polygons, both compartments from tumor "
                     f"sections"}.get(variant, variant)


def _sections(cfg, diag, reg):
    """Which sections each variant draws its compartments from, with the growth-plate bins of non-tumor sections."""
    cp = cfg["compartments"]
    GP, _, TA = _names(cfg)
    tumor_secs = set(common.tumor_and_reference(cfg)[0])
    tumor_txt = ", ".join(cp["tumor_conditions"])

    def off(df, col):
        if df is None or col not in df.columns:
            return np.nan
        rows = df[(df["sample"].astype(str) != "all") & ~df["sample"].astype(str).isin(tumor_secs)]
        return pd.to_numeric(rows[col], errors="coerce").fillna(0).sum()

    if cp["de_sections"] == "tumor":
        first = f"take both compartments from {tumor_txt} sections (de_sections = tumor)"
    else:
        allrow = diag[diag["sample"].astype(str) == "all"] if diag is not None else None
        n_all = allrow["n_growth_plate"].iloc[0] if allrow is not None and len(allrow) else np.nan
        first = (f"follow de_sections = all: {GP} bins of every section ({_num(off(diag, 'n_growth_plate'))} "
                 f"of {_num(n_all)} from non-{tumor_txt} sections), {TA} bins of {tumor_txt} sections")
    k = off(reg, "n_osteoclast_bins_near_regions")
    return (f"Sections: the manuscript, ROI-only, depth-adjusted and block-pseudobulk variants {first}; the anatomical "
            f"and drawn variants take both compartments from {tumor_txt} sections only, a contrast within the "
            f"tumor-bearing sections"
            + (f" ({_num(k)} anatomical {GP} bins of non-{tumor_txt} sections left out)"
               if np.isfinite(k) and k > 0 else "") + ".")


def _caption(cfg, var, top, diag, reg):
    cp, st = cfg["compartments"], common.settings(cfg, "compartments")
    S, P = pc.signature_label(cfg, cp["signature"]), pc.signature_label(cfg, cp["split_by"])
    GP, GPn, TA = _names(cfg)
    labels = _labels(cfg)
    parts = [f"{S[:1].upper() + S[1:]} compartments of Figure 4h under other definitions and tests."]
    allrow = diag[diag["sample"].astype(str) == "all"].iloc[0] if diag is not None and \
        (diag["sample"].astype(str) == "all").any() else None
    if allrow is not None:
        zero = common.as_bool([allrow.get("threshold_is_zero")])[0]
        parts.append(
            f"Effective rule: the {float(cp['percentile']):g}th percentile of the {P} score over all "
            f"{_num(allrow.get('n_bins'))} bins is {_num(allrow.get('threshold'), '{:.4g}')}"
            + (f", exactly 0, so a bin is '{GPn}' whenever its {P} score is above 0, in practice when a "
               f"{P} gene is detected" if zero else "")
            + f": {_num(allrow.get('pct_bins_above_threshold'), '{:.1f}')}% of bins are above it, and "
              f"{_num(allrow.get('n_bins_with_split_gene_umi'))} bins have at least one UMI of a {P} gene. "
              f"Manuscript compartments: {_num(allrow.get('n_growth_plate'))} {GP} bins "
              f"({_num(allrow.get('n_growth_plate_outside_roi'))} outside the marrow region of interest, "
              f"{_num(allrow.get('n_growth_plate_one_split_umi'))} with a single {P} UMI, "
              f"{_num(allrow.get('n_growth_plate_without_core_gene'))} without any {P} core gene) and "
              f"{_num(allrow.get('n_tumor_associated'))} {TA} bins "
              f"({_num(allrow.get('n_tumor_associated_outside_roi'))} outside the region of interest).")
        g, t = allrow.get("median_umi_growth_plate"), allrow.get("median_umi_tumor_associated")
        ratio = common.fold(g, t)
        parts.append(f"Median UMI per bin: {GPn} {_num(g, '{:g}')}, {TA} {_num(t, '{:g}')} "
                     f"({_num(ratio, '{:.2f}')}-fold).")
    if reg is not None and len(reg):
        sec = "; ".join(f"{r['sample']}: {_num(r.get('n_regions'))} region(s), "
                        f"{_num(r.get('n_osteoclast_bins_near_regions'))} {S}-positive bins within reach "
                        f"({_num(r.get('n_osteoclast_bins_near_regions_manuscript_growth_plate'))} of them manuscript "
                        f"{GP} bins, of {_num(r.get('n_manuscript_growth_plate'))})"
                        for _, r in reg.iterrows())
        calls = ", ".join(sorted({str(c) for c in reg.get("calls", pd.Series(dtype=str)).dropna() if str(c)})) \
            or "no calls"
        parts.append(f"Anatomical proxy: {P}-dense regions = 8-connected {P}-positive bins ({calls}) of at least "
                     f"{int(st['min_region_bins'])} bins on the array grid, {GPn} = {S}-positive bins within "
                     f"{float(st['anatomical_distance_um']):g} um of one ({sec}).")
    an = var[var["variant"] == "anatomical"]
    if len(an):
        r = an.iloc[0]
        parts.append(f"The anatomical contrast is {r['status']} ({_num(r['n_group'])} {GP} vs "
                     f"{_num(r['n_rest'])} {TA} bins of tumor sections).")
    parts.append(_sections(cfg, diag, reg))
    n_reps = int(st["n_depth_reps"])
    desc = []
    for _, r in var.iterrows():
        v = str(r["variant"])
        if r["status"] == "not configured":
            continue
        txt = (f"{_label(v, labels)}: {_what(v, S, P, st, n_reps, GPn)}; "
               f"n = {_num(r['n_group'])} / {_num(r['n_rest'])} {_units_short(r['units'])}")
        if r["status"] == "estimable":
            txt += (f", genes significantly higher in the {GPn} / {TA} compartment: "
                    f"{_num(r['n_significant_growth_plate'])} / {_num(r['n_significant_tumor_associated'])}")
            if v in ("depth_adjusted", "block_pseudobulk") and \
                    np.isfinite(pd.to_numeric(r["median_umi_group"], errors="coerce")):
                txt += (f" (median UMI per {'bin' if v == 'depth_adjusted' else 'unit'}, {GPn} vs "
                        f"{TA}: {_num(r.get('median_umi_group_unadjusted'), '{:g}')} vs "
                        f"{_num(r.get('median_umi_rest_unadjusted'), '{:g}')} before, "
                        f"{_num(r['median_umi_group'], '{:g}')} vs {_num(r['median_umi_rest'], '{:g}')} after the "
                        f"adjustment)")
        else:
            txt += f", {r['status']}"
        desc.append(txt)
    parts.append(f"Columns (n = {GP} / {TA} units): " + "; ".join(desc) + ".")
    for _, r in var[var["status"] == "not configured"].iterrows():
        parts.append(f"Not run: {_label(r['variant'], labels).lower()}, {r['note']}.")
    sig_rows = []
    if top is not None and len(top):
        for v in var.loc[var["status"] == "estimable", "variant"].astype(str):
            col = f"significant_{v}"
            if col in top.columns:
                sig_rows.append(f"{_label(v, labels)} {int(common.as_bool(top[col]).sum())}")
        sel = top.loc[common.as_bool(top["selection_gene"]), "gene"].astype(str).tolist()
        n_gp = int((top["group"].astype(str) == GROUPS[0]).sum())
        parts.append(f"Rows: the {len(top)} genes of Figure 4h ({n_gp} {GP}, {len(top) - n_gp} "
                     f"{TA}: the top genes per compartment by manuscript Wilcoxon score); filled = "
                     f"BH-adjusted P < {ALPHA:g} in that variant (depth-adjusted and block pseudobulk: median adjusted "
                     f"P over replicates), "
                     f"colour = the compartment the gene is higher in ({P} colour = {GPn}, {S} colour = "
                     f"{TA}), open = not significant, gray band = variant not estimable. Of the "
                     f"{len(top)} genes, significant in their Figure 4h direction: " + ", ".join(sig_rows) + ". "
                     + (f"† = gene of the {P} extended list ({', '.join(sel)}), which selects the {GP} bins "
                        f"and is higher there by construction." if sel else
                        f"No drawn gene belongs to the {P} extended list."))
    else:
        parts.append("The manuscript test was not estimable, so no gene is drawn.")
    return " ".join(parts)


def draw(results, cfg, out):
    """Sensitivity_Figure4h_variants.{png,pdf}: dot matrix of the Figure 4h genes x variants. Returns the manifest
    entries ([] when the variants table is absent)."""
    var = common.read_csv(results, VARIANTS_CSV)
    if var is None or "variant" not in var.columns:
        return []
    top, diag, reg = (common.read_csv(results, f) for f in (TOP_CSV, DIAG_CSV, REGIONS_CSV))
    var = var.assign(variant=var["variant"].astype(str), status=var["status"].astype(str))
    shown = var[var["status"] != "not configured"].reset_index(drop=True)
    cp = cfg["compartments"]
    S, P = pc.signature_label(cfg, cp["signature"]), pc.signature_label(cfg, cp["split_by"])
    _, GPn, TA = _names(cfg)
    labels = _labels(cfg)
    c_gp, c_ta = _color(cfg, cp["split_by"], "#33A02C"), _color(cfg, cp["signature"], "#C44E52")
    genes = top if top is not None else pd.DataFrame(columns=["gene", "group", "selection_gene"])
    genes = pd.concat([genes[genes["group"].astype(str) == g] for g in GROUPS], ignore_index=True)
    n_g, n_v = len(genes), max(len(shown), 1)
    with common.house_style():
        fig = Figure(figsize=(1.6 + 0.78 * n_v, 1.5 + 0.2 * max(n_g, 3)))
        ax = fig.subplots()
        for j, r in shown.iterrows():
            v = r["variant"]
            if r["status"] != "estimable":
                ax.axvspan(j - .45, j + .45, color=LIGHT, lw=0, zorder=0)
                ax.text(j, (max(n_g, 1) - 1) / 2, "not estimable", rotation=90, ha="center", va="center", fontsize=7,
                        color=GRAY)
                continue
            if not n_g or f"padj_{v}" not in genes.columns:
                continue
            padj = pd.to_numeric(genes[f"padj_{v}"], errors="coerce").to_numpy()
            score = pd.to_numeric(genes.get(f"score_{v}", pd.Series(np.nan, index=genes.index)),
                                  errors="coerce").to_numpy()
            in_gp = np.where(genes["group"].astype(str).to_numpy() == GROUPS[0], score > 0, score < 0)
            for i in range(n_g):
                if not np.isfinite(padj[i]) or not np.isfinite(score[i]):
                    ax.plot(j, i, marker="x", ms=4, color=GRAY, mew=.8)
                    continue
                col = c_gp if in_gp[i] else c_ta
                ax.plot(j, i, marker="o", ms=7, mec=col, mew=1.1, mfc=col if padj[i] < ALPHA else "white")
        if n_g:
            sel = common.as_bool(genes["selection_gene"])
            ax.set_yticks(range(n_g))
            ax.set_yticklabels([f"{g} †" if s else str(g) for g, s in zip(genes["gene"], sel)], style="italic",
                               fontsize=7.5)
            n_first = int((genes["group"].astype(str) == GROUPS[0]).sum())
            if 0 < n_first < n_g:
                ax.axhline(n_first - .5, color="k", lw=.5)
            for lo, hi, txt in ((0, n_first, f"higher in {GPn}\n(Figure 4h)"),
                                (n_first, n_g, f"higher in {TA}\n(Figure 4h)")):
                if hi > lo:
                    ax.text(n_v - .35, (lo + hi - 1) / 2, txt, rotation=270, ha="left", va="center", fontsize=6.5)
        else:
            ax.set_yticks([])
            ax.text((n_v - 1) / 2, 0, "manuscript test not estimable: no gene to show", ha="center", va="center",
                    fontsize=7.5, color=GRAY)
        ax.set_xticks(range(len(shown)))
        ax.set_xticklabels([f"{labels.get(r['variant'], r['variant'])}\n{_num(r['n_group'])} / {_num(r['n_rest'])} "
                            f"{_units_short(r['units'])}" for _, r in shown.iterrows()], fontsize=6.5)
        ax.xaxis.tick_top()
        ax.tick_params(length=0)
        ax.set_xlim(-.6, n_v - .4)
        ax.set_ylim(max(n_g, 1) - .5, -.5)
        for side in ("top", "right", "bottom", "left"):
            ax.spines[side].set_visible(False)
        handles = [Line2D([], [], ls="", marker="o", ms=6, mfc=c_gp, mec=c_gp,
                          label=f"adj. P < {ALPHA:g}, higher in {GPn} ({P})"),
                   Line2D([], [], ls="", marker="o", ms=6, mfc=c_ta, mec=c_ta,
                          label=f"adj. P < {ALPHA:g}, higher in {TA} ({S})"),
                   Line2D([], [], ls="", marker="o", ms=6, mfc="white", mec=GRAY,
                          label=f"adj. P ≥ {ALPHA:g} (edge colour = direction)"),
                   Line2D([], [], ls="", marker="x", ms=4, color=GRAY, label="gene not tested")]
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=2, frameon=False, fontsize=6.5)
        common.save_figure(fig, out, NAME, dpi=cfg["figures"]["dpi"])
    return [common.figure_entry(NAME, TITLE, _caption(cfg, var, top, diag, reg))]
