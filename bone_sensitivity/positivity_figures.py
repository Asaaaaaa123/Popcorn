"""
Figure of the positivity sensitivity analysis (bone_sensitivity.positivity): tumor/control positivity folds of the
Figure 4g/6a signatures under every region of interest, calling rule and depth-matching scheme, as a forest plot.
Reads sensitivity_positivity_folds.csv (and sensitivity_roi_composition.csv for the caption when present).
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

import pipeline_config as pc
from bone_sensitivity import common

NAME, TITLE = "Sensitivity_positivity_folds", "Sensitivity — tumor/control positivity folds"
FIGURES = [(NAME, TITLE)]
FOLDS_CSV, COMPOSITION_CSV = "sensitivity_positivity_folds.csv", "sensitivity_roi_composition.csv"
INPUTS = [FOLDS_CSV]
DEPTH_LABELS = {"native": "native", "cap": "cap", "quantile": "quantile", "quantile_roi": "quantile (ROI)"}
ROI_LABELS = {"call_rule_exclusion": "ROI: call-rule exclusion", "assignment_exclusion": "ROI: assignment exclusion",
              "no_lesion": "ROI: lesion removed"}
ROI_TEXT = {"call_rule_exclusion": "muscle and nerve bins removed by the positivity rule, not the percentile rule",
            "assignment_exclusion": "bins assigned to muscle or nerve in Figure 3 removed",
            "no_lesion": "A549 lesion bins removed"}
DEPTH_TEXT = {"cap": "cap, bins above the shallower section's median UMI downsampled to it",
              "quantile": "quantile, whole sections matched to the shallower section's UMI distribution",
              "quantile_roi": "quantile (ROI), region-of-interest bins matched to the shallower section's ROI "
                              "distribution"}
_FALLBACK_COLORS = ("#C44E52", "#55A868", "#8172B3", "#CCB974", "#64B5CD", "#DD8452")


def _rule_label(rule):
    m = re.fullmatch(r"min_core_(\d+)", str(rule))
    if m:
        k = int(m.group(1))
        return f"rule: ≥{k} core gene{'s' if k != 1 else ''}"
    return {"no_competitor_veto": "rule: no competitor veto"}.get(rule, f"rule: {rule}")


def _row_label(depth, rule, roi):
    d = DEPTH_LABELS.get(depth, depth)
    if rule == "manuscript" and roi == "manuscript":
        return d
    return f"{d} · {_rule_label(rule) if rule != 'manuscript' else ROI_LABELS.get(roi, f'ROI: {roi}')}"


def variant_rows(t):
    """(depth, rule, roi) rows of the forest plot in drawing order: native manuscript, its region and rule variants;
    every depth-matched scheme's manuscript analysis; then the region and rule variants of each depth-matched scheme
    that has them. Combinations of a non-manuscript region with a non-manuscript rule are not drawn."""
    present = set(zip(t["depth"], t["rule"], t["roi"]))
    depths = list(dict.fromkeys(t["depth"]))
    depths = (["native"] if "native" in depths else []) + [d for d in depths if d != "native"]
    rois = [r for r in dict.fromkeys(t["roi"]) if r != "manuscript"]
    rules = [r for r in dict.fromkeys(t["rule"]) if r != "manuscript"]

    def variants(d):
        return [(d, "manuscript", r) for r in rois] + [(d, r, "manuscript") for r in rules]

    order = []
    if "native" in depths:
        order += [("native", "manuscript", "manuscript")] + variants("native")
    order += [(d, "manuscript", "manuscript") for d in depths if d != "native"]
    for d in depths:
        if d != "native":
            order += variants(d)
    return [r for r in order if r in present]


def _log_ticks(lo, hi):
    """Three to about six readable ticks between lo and hi on a log axis (1-3-10, 1-2-5 or finer steps)."""
    ticks = []
    for steps in ((1,), (1, 3), (1, 2, 5), (1, 1.5, 2, 3, 5, 7)):
        ticks = [m * 10.0 ** e for e in range(int(np.floor(np.log10(lo))) - 1, int(np.ceil(np.log10(hi))) + 1)
                 for m in steps]
        ticks = [v for v in ticks if lo <= v <= hi]
        if len(ticks) >= 3:
            break
    return ticks


def _signatures(cfg, t):
    names = list(dict.fromkeys(s for p in (cfg.get("figures") or {}).get("positivity_panels") or []
                               for s in p.get("signatures") or []))
    have = set(t["signature"])
    names = [s for s in names if s in have]
    return names or list(dict.fromkeys(t["signature"]))


def _label(cfg, name):
    try:
        return pc.signature_label(cfg, name)
    except (KeyError, TypeError):
        return str(name).replace("_", " ")


def _fmt(v):
    return f"{v:.2f}" if np.isfinite(v) else "n/a"


def _who(cfg, sections):
    """How the caption names the plotted (non-reference) sections."""
    try:
        tumor = set(common.tumor_and_reference(cfg)[0])
    except (KeyError, TypeError):
        tumor = set()
    if sections and set(sections) <= tumor:
        return "the tumor-bearing section" if len(sections) == 1 else "each tumor-bearing section"
    return "the non-reference section" if len(sections) == 1 else "each non-reference section"


def draw(results_dir, cfg, out_dir):
    """Forest plot of the folds (log axis), one panel per signature of figures.positivity_panels: every variant row
    (variant_rows), each non-reference section in its color, the manuscript analyses shaded, Katz intervals (thin bars)
    and, on the native manuscript-rule rows, spatial block-bootstrap intervals (thick pale bars) where they are finite.
    A fold of 0 sits at the left edge (open triangle); a fold that is not estimable is marked n.e."""
    t = common.read_csv(results_dir, FOLDS_CSV)
    if t is None or t.empty or not {"depth", "rule", "roi", "signature", "sample", "fold"} <= set(t.columns):
        return []
    for c in ("depth", "rule", "roi", "signature", "sample"):
        t[c] = t[c].astype(str)
    for c in ("fold", "fold_ci_lo", "fold_ci_hi", "block_fold_ci_lo", "block_fold_ci_hi", "n_positive",
              "n_positive_reference"):
        t[c] = pd.to_numeric(t[c], errors="coerce") if c in t else np.nan
    t["manuscript"] = common.as_bool(t["manuscript"]) if "manuscript" in t else False
    sigs = _signatures(cfg, t)
    rows = variant_rows(t)
    pos = {r: i for i, r in enumerate(rows)}
    samples = {s["id"]: s for s in cfg.get("samples") or []}
    sections = list(dict.fromkeys(t["sample"]))
    color = {s: samples.get(s, {}).get("color") or _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)]
             for i, s in enumerate(sections)}
    name = {s: samples.get(s, {}).get("label") or s for s in sections}
    ref_id = cfg.get("reference_sample")
    ref_name = (samples.get(ref_id) or {}).get("label") or ref_id or "the reference section"
    offset = {s: (i - (len(sections) - 1) / 2) * min(.5 / max(len(sections), 1), .22) for i, s in enumerate(sections)}

    ncol = min(4, len(sigs))
    nrow = int(np.ceil(len(sigs) / ncol))
    height = .55 + .14 * len(rows)
    dpi = int(((cfg.get("figures") or {}).get("dpi")) or 200)
    with common.house_style():
        fig = Figure(figsize=(7.2, nrow * height + .6))
        axes = np.atleast_1d(fig.subplots(nrow, ncol, sharey=True, squeeze=False)).ravel()
        manuscript_rows = [i for i, r in enumerate(rows) if t[(t["depth"] == r[0]) & (t["rule"] == r[1]) &
                                                               (t["roi"] == r[2])]["manuscript"].any()]
        n_native = sum(1 for r in rows if r[0] == "native")
        for ax, sig in zip(axes, sigs):
            sub = t[t["signature"] == sig]
            vals = np.concatenate([sub[c].to_numpy(float) for c in ("fold", "fold_ci_lo", "fold_ci_hi",
                                                                     "block_fold_ci_lo", "block_fold_ci_hi")]
                                  + [[1.0]])
            vals = vals[np.isfinite(vals) & (vals > 0)]
            lo_lim, hi_lim = (vals.min() / 1.35, vals.max() * 1.35) if vals.size else (.5, 2.0)
            for i in manuscript_rows:
                ax.axhspan(i - .45, i + .45, color="#ececec", zorder=0, lw=0)
            if 0 < n_native < len(rows):
                ax.axhline(n_native - .5, color="#9a9a9a", lw=.5, ls=":")
            ax.axvline(1.0, color="#555555", lw=.6, ls="--", zorder=1)
            for r in sub.itertuples(index=False):
                key = (r.depth, r.rule, r.roi)
                if key not in pos:
                    continue
                y = pos[key] + offset[r.sample]
                c, filled = color[r.sample], bool(r.manuscript)
                if np.isfinite(r.fold) and r.fold > 0:
                    if np.isfinite(r.block_fold_ci_lo):                  # upper end NaN: unbounded, to the edge
                        bhi = min(r.block_fold_ci_hi, hi_lim) if np.isfinite(r.block_fold_ci_hi) else hi_lim
                        ax.plot([max(r.block_fold_ci_lo, lo_lim), bhi], [y, y], color=c, lw=2.8, alpha=.3,
                                solid_capstyle="butt", zorder=1.5)
                    if np.isfinite(r.fold_ci_lo) and np.isfinite(r.fold_ci_hi):
                        ax.plot([max(r.fold_ci_lo, lo_lim), min(r.fold_ci_hi, hi_lim)], [y, y], color=c, lw=.8,
                                zorder=2)
                    ax.plot([r.fold], [y], marker="s" if filled else "o", ms=4.2 if filled else 3.4,
                            mfc=c if filled else "white", mec="black" if filled else c, mew=.7, zorder=3, ls="none")
                elif np.isfinite(r.fold):
                    ax.plot([lo_lim * 1.04], [y], marker="<", ms=3.6, mfc="white", mec=c, mew=.7, ls="none", zorder=3)
                else:
                    ax.text(hi_lim / 1.04, y, "n.e.", fontsize=5.5, color="#777777", ha="right", va="center")
            ax.set_xscale("log")
            ax.set_xlim(lo_lim, hi_lim)
            ax.xaxis.set_major_locator(FixedLocator(_log_ticks(lo_lim, hi_lim)))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_title(_label(cfg, sig), fontsize=8)
            ax.tick_params(labelsize=6.5)
            ax.set_xlabel("fold vs " + ref_name + " (log)", fontsize=6.5)
        for ax in axes[len(sigs):]:
            ax.set_visible(False)
        axes[0].set_yticks(range(len(rows)))
        axes[0].set_yticklabels([_row_label(*r) for r in rows])
        axes[0].set_ylim(len(rows) - .5, -.5)
        has_block = bool(np.isfinite(t["block_fold_ci_lo"].to_numpy(float)).any())
        handles = [Line2D([], [], marker="s", ls="none", mfc="#888888", mec="black", label="manuscript analysis"),
                   Line2D([], [], marker="o", ls="none", mfc="white", mec="#888888", label="sensitivity variant"),
                   Line2D([], [], color="#888888", lw=.8, label="Katz 95%")]
        if has_block:
            handles.append(Line2D([], [], color="#888888", lw=2.8, alpha=.3, solid_capstyle="butt",
                                  label="block bootstrap 95%"))
        if len(sections) > 1:
            handles += [Line2D([], [], marker="o", ls="none", mfc=color[s], mec=color[s], label=name[s])
                        for s in sections]
        fig.legend(handles=handles, fontsize=6.5, frameon=False, loc="lower center",
                   ncol=len(handles) if len(handles) <= 5 else int(np.ceil(len(handles) / 2)), bbox_to_anchor=(.5, 1.0))
        fig.tight_layout()
        common.save_figure(fig, out_dir, NAME, dpi=dpi)
    caption = _caption(results_dir, cfg, t, rows, sigs, sections, name, ref_name, has_block)
    return [common.figure_entry(NAME, TITLE, caption)]


def _caption(results_dir, cfg, t, rows, sigs, sections, name, ref_name, has_block=False):
    sco = cfg.get("scoring") or {}
    reps = int(sco.get("n_depth_reps", 0) or 0)
    depths = list(dict.fromkeys(r[0] for r in rows if r[0] != "native"))
    rois = list(dict.fromkeys(r[2] for r in rows if r[2] != "manuscript"))
    rules = list(dict.fromkeys(r[1] for r in rows if r[1] != "manuscript"))
    full = list(dict.fromkeys(r[0] for r in rows if r[0] != "native" and (r[1], r[2]) != ("manuscript", "manuscript")))
    text = (f"Fold of the percentage of region-of-interest bins called positive in {_who(cfg, sections)} over "
            f"{ref_name}, per signature, under alternative analysis choices. Native depth: the manuscript analysis "
            f"(transcriptional region of interest, score above the {sco.get('null_percentile', 95):g}th-percentile "
            f"null threshold, at least {sco.get('min_core', 2)} core genes, competitor veto)")
    if rois:
        text += "; regions of interest " + "; ".join(ROI_TEXT.get(r, r) for r in rois)
    if rules:
        text += "; calling rules " + "; ".join(_rule_label(r).replace("rule: ", "") for r in rules)
    text += ". "
    if depths:
        spread = f"mean over {reps} replicates" if reps > 1 else "one replicate"
        text += (f"Depth-matched ({spread}; scores and null thresholds recomputed per replicate): "
                 + "; ".join(DEPTH_TEXT.get(d, d) for d in depths))
        if full:
            text += f"; the region and rule variants under {', '.join(DEPTH_LABELS.get(d, d) for d in full)}"
        text += ". "
    text += ("Squares on gray, the manuscript analyses (native and the configured depth matching); circles, "
             "sensitivity variants. Thin bars, Katz 95% intervals of a ratio of two independent proportions (bins "
             "treated as independent, so too narrow where positive bins cluster; on the depth-matched rows computed "
             "from the replicate-mean counts, so the spread over the downsampling replicates is not included; not "
             "shown when a count is 0). ")
    if has_block:
        top = common.settings(cfg, None)
        text += (f"Thick pale bars (native depth, manuscript calling rule), spatial block-bootstrap 95% intervals: "
                 f"{float(top['block_um']):g} µm blocks of each section resampled independently "
                 f"({int(top['n_boot'])} resamples), percentile interval of the ratio, which allows for the clustering "
                 f"of positive bins; a thick bar running to the right edge has no upper bound (no positive reference "
                 f"bin in more than 2.5% of resamples). ")
    text += ("Logarithmic axis; dashed line, fold 1; open triangles, fold 0; n.e., not estimable (no positive "
             "reference bin). ")
    drawn = t[[(d, r, o) in set(rows) for d, r, o in zip(t["depth"], t["rule"], t["roi"])]]
    parts = []
    for s in sections:
        items = []
        for sig in sigs:
            sub = drawn[(drawn["signature"] == sig) & (drawn["sample"] == s)]
            f = sub["fold"].to_numpy(float)
            f = f[np.isfinite(f)]
            nat = sub[(sub["depth"] == "native") & sub["manuscript"]]
            dm = sub[(sub["depth"] != "native") & sub["manuscript"]]["fold"].to_numpy(float)
            if not f.size:
                items.append(f"{_label(cfg, sig)} not estimable")
                continue
            item = f"{_label(cfg, sig)} {_fmt(float(nat['fold'].iloc[0])) if len(nat) else 'n/a'}"
            inner = []
            if len(nat):
                r = nat.iloc[0]
                if np.isfinite(r["fold_ci_lo"]) and np.isfinite(r["fold_ci_hi"]):
                    inner.append(f"Katz {r['fold_ci_lo']:.2f}–{r['fold_ci_hi']:.2f}")
                if np.isfinite(r["block_fold_ci_lo"]):
                    hi = f"{r['block_fold_ci_hi']:.2f}" if np.isfinite(r["block_fold_ci_hi"]) else "unbounded"
                    inner.append(f"block bootstrap {r['block_fold_ci_lo']:.2f}–{hi}")
            if dm.size:
                inner.append(f"depth-matched {_fmt(dm[0])}")
            if inner:
                item += " (" + "; ".join(inner) + ")"
            items.append(item + f", range {f.min():.2f}–{f.max():.2f}")
        prefix = f"{name[s]}: " if len(sections) > 1 else ""
        parts.append(prefix + "; ".join(items))
    text += ("Manuscript native fold (its intervals; the configured depth-matched fold) and range over the variants "
             "drawn: " + " | ".join(parts) + ".")
    comp = common.read_csv(results_dir, COMPOSITION_CSV)
    if comp is not None and {"roi", "sample", "pct_skeletal_muscle"} <= set(comp.columns):
        labels = {s["id"]: s.get("label") or s["id"] for s in cfg.get("samples") or []}
        m = comp[comp["roi"].astype(str) == "manuscript"]
        pct = pd.to_numeric(m["pct_skeletal_muscle"], errors="coerce")
        shares = [f"{labels.get(sid, sid)} {v:.1f}%" for sid, v in zip(m["sample"].astype(str), pct) if np.isfinite(v)]
        if shares:
            text += (" Share of the manuscript region of interest assigned to skeletal muscle (Figure 3): "
                     + ", ".join(shares) + ".")
    return text
