"""
Figure of the proximity sensitivity table (bone_sensitivity/proximity.py, sensitivity_proximity.csv): for every
proximity pair, call set and section, the observed "within R um" percentage next to the all-ROI, non-target and
depth-matched baselines, and the depth-matched difference with its spatial block-bootstrap interval (also with the
ROI-quantile strata alone, which leave residual depth confounding). The wording about target bins follows the
table's exclude_self column (cfg proximity.exclude_self for tables without it).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import pipeline_config as pc
from bone_sensitivity import common

FILE = "sensitivity_proximity.csv"
NAME = "Sensitivity_proximity_baselines"
FIGURES = [(NAME, "Sensitivity — proximity against other baselines")]
INPUTS = [FILE]

BASELINES = [("pct_all_within", "all ROI bins (manuscript)", "0.82"),
             ("pct_nontarget_within", "non-target ROI bins", "0.62"),
             ("expected_depth_matched_pct", "depth-matched expectation", "0.38")]
CALLS = {"program": "assignment", "program_dm": "depth-matched assignment", "pos": "positivity calls"}
_FALLBACK = ("#4C72B0", "#C44E52", "#55A868", "#8172B3", "#CCB974", "#64B5CD")


def _num(df, col):
    return pd.to_numeric(df[col], errors="coerce").to_numpy(float) if col in df.columns else np.full(len(df), np.nan)


def _f(x, fmt=".1f", suffix=""):
    return f"{x:{fmt}}{suffix}" if np.isfinite(x) else "n/a"


def _v(r, col):
    """Value of a table row as a float (NaN when missing or not numeric)."""
    x = pd.to_numeric(pd.Series([r.get(col, np.nan)]), errors="coerce").iloc[0]
    return float(x) if pd.notna(x) else np.nan


def _i(r, col):
    x = _v(r, col)
    return int(x) if np.isfinite(x) else 0


def _sections(cfg, ids):
    """(id, label, color) of the sections in the table, in configuration order."""
    conf = {s["id"]: s for s in cfg.get("samples", [])}
    order = [s for s in conf if s in ids] + [s for s in ids if s not in conf]
    return [(s, conf.get(s, {}).get("label") or s, conf.get(s, {}).get("color") or _FALLBACK[i % len(_FALLBACK)])
            for i, s in enumerate(order)]


def _set_title(cfg, r):
    src, dst = pc.signature_label(cfg, str(r["source"])), pc.signature_label(cfg, str(r["target"]))
    return f"{src} → {dst} ({CALLS.get(str(r['assignment_column']), str(r['assignment_column']))})"


def _row_text(cfg, r, label):
    head = f"{_set_title(cfg, r)}, {label}"
    obs, n_s, n_t = _v(r, "pct_source_within"), _i(r, "n_source"), _i(r, "n_target")
    if str(r["status"]) != "estimable":
        return (f"{head}: not estimable ({n_s} source, {n_t} target bins; observed {_f(obs)}% vs "
                f"{_f(_v(r, 'pct_all_within'))}% of all ROI bins)")
    text = (f"{head}: {_f(obs)}% of {n_s} source bins vs {_f(_v(r, 'pct_all_within'))}% all ROI, "
            f"{_f(_v(r, 'pct_nontarget_within'))}% non-target, {_f(_v(r, 'expected_depth_matched_pct'))}% "
            f"depth-matched (difference {_f(_v(r, 'diff_pp'), '+.1f')} pp, 95% CI "
            f"{_f(_v(r, 'diff_ci_lo'), '+.1f')} to {_f(_v(r, 'diff_ci_hi'), '+.1f')}")
    text += (f"; median UMI {_f(_v(r, 'median_umi_source'), '.0f')} source vs "
             f"{_f(_v(r, 'median_umi_matched_reference'), '.0f')} matched reference, "
             f"{_f(_v(r, 'median_umi_roi'), '.0f')} ROI")
    if _i(r, "n_copositive") > 0:
        text += (f"; excluding the {_i(r, 'n_copositive')} co-positive source bins "
                 f"{_f(_v(r, 'diff_excl_copositive_pp'), '+.1f')} pp, 95% CI "
                 f"{_f(_v(r, 'diff_excl_copositive_ci_lo'), '+.1f')} to "
                 f"{_f(_v(r, 'diff_excl_copositive_ci_hi'), '+.1f')}")
    if np.isfinite(_v(r, "diff_roi_strata_pp")):
        text += (f"; ROI-quantile strata alone {_f(_v(r, 'expected_roi_strata_pct'))}% expected, "
                 f"{_f(_v(r, 'diff_roi_strata_pp'), '+.1f')} pp (95% CI {_f(_v(r, 'diff_roi_strata_ci_lo'), '+.1f')} "
                 f"to {_f(_v(r, 'diff_roi_strata_ci_hi'), '+.1f')}; matched reference median "
                 f"{_f(_v(r, 'median_umi_matched_reference_roi_strata'), '.0f')} UMI)")
    return text + f"; largest ratio the all-ROI baseline allows {_f(_v(r, 'max_possible_ratio_vs_all'), '.2f')})"


def _exclude_self(df, cfg):
    """exclude_self of the table (its column), else of the configuration."""
    if "exclude_self" in df.columns and df["exclude_self"].notna().any():
        return bool(common.as_bool(df.loc[df["exclude_self"].notna(), "exclude_self"]).any())
    return bool(common.as_bool([cfg.get("proximity", {}).get("exclude_self", False)])[0])


def _bars(ax, sub, sections):
    """Observed (section color) and the three baselines per section; 'n.e.' where the expectation is missing."""
    w = .19
    top = 0.0
    for j, (sid, label, color) in enumerate(sections):
        r = sub[sub["sample"] == sid]
        if r.empty:
            continue
        vals = [_num(r, "pct_source_within")[0]] + [_num(r, c)[0] for c, _, _ in BASELINES]
        colors = [color] + [g for _, _, g in BASELINES]
        for i, (v, c) in enumerate(zip(vals, colors)):
            x = j + (i - 1.5) * w
            if np.isfinite(v):
                ax.bar(x, v, w * .92, color=c, linewidth=0)
                top = max(top, v)
            elif i == 3:
                ax.text(x, 1.5, "n.e.", ha="center", va="bottom", fontsize=6, color="0.3", rotation=90)
    ax.set_xticks(range(len(sections)))
    ax.set_xticklabels([f"{label}\nn = {_i(sub[sub['sample'] == sid].iloc[0], 'n_source')}"
                        if (sub["sample"] == sid).any() else label for sid, label, _ in sections], fontsize=7)
    ax.set_xlim(-.6, len(sections) - .4)
    ax.set_ylim(0, min(100.0, max(10.0, top * 1.15)))
    ax.set_ylabel("% within radius", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)


def _forest(ax, sub, sections):
    """Depth-matched difference (filled circle), the same without co-positive source bins (open circle) and with the
    ROI-quantile strata alone (grey square), 95% block intervals."""
    ys = np.arange(len(sections))[::-1].astype(float)
    ax.axvline(0, color="0.5", linewidth=.6, linestyle="--", zorder=0)
    for y, (sid, label, color) in zip(ys, sections):
        r = sub[sub["sample"] == sid]
        d = _num(r, "diff_pp")[0] if len(r) else np.nan
        if not np.isfinite(d):
            ax.text(0, y, "not estimable", ha="center", va="center", fontsize=6, color="0.3",
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
            continue
        series = [("diff_pp", "diff_ci_lo", "diff_ci_hi", color, color, "o", 0.0),
                  ("diff_roi_strata_pp", "diff_roi_strata_ci_lo", "diff_roi_strata_ci_hi", "0.6", "0.6", "s", .22)]
        if _num(r, "n_copositive")[0] > 0:
            series.append(("diff_excl_copositive_pp", "diff_excl_copositive_ci_lo", "diff_excl_copositive_ci_hi",
                           "white", color, "o", -.22))
        for c, lo_c, hi_c, face, edge, marker, dy in series:
            v, lo, hi = _num(r, c)[0], _num(r, lo_c)[0], _num(r, hi_c)[0]
            if not np.isfinite(v):
                continue
            if np.isfinite(lo) and np.isfinite(hi):
                ax.plot([lo, hi], [y + dy, y + dy], color=edge, linewidth=1.2, solid_capstyle="butt")
            ax.plot([v], [y + dy], marker, markersize=4.5 if marker == "o" else 3.8, markerfacecolor=face,
                    markeredgecolor=edge, markeredgewidth=1)
    ax.set_yticks(ys)
    ax.set_yticklabels([label for _, label, _ in sections], fontsize=7)
    ax.set_ylim(-.7, len(sections) - .3)
    ax.set_xlabel("observed − depth-matched (pp)", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    lim = max(5.0, np.nanmax(np.abs(np.r_[ax.get_xlim(), 0.0])))
    ax.set_xlim(-lim, lim)


def draw(results, cfg, out):
    """Sensitivity_proximity_baselines from sensitivity_proximity.csv; [] when the table is absent or empty."""
    df = common.read_csv(results, FILE)
    if df is None or df.empty:
        return []
    for c in ("pair", "assignment_column", "sample", "status", "source", "target"):
        df[c] = df[c].astype(str)
    sets = list(dict.fromkeys(zip(df["pair"], df["assignment_column"])))
    sections = _sections(cfg, list(pd.unique(df["sample"])))
    labels = {sid: label for sid, label, _ in sections}
    st, top = common.settings(cfg, "proximity"), common.settings(cfg, None)
    radius = float(np.nanmax(_num(df, "radius_um"))) if np.isfinite(_num(df, "radius_um")).any() else np.nan
    exclude_self = _exclude_self(df, cfg)
    texts = []
    with common.house_style():
        height = 1.1 + 1.85 * len(sets)
        fig = Figure(figsize=(6.5, height))
        axes = fig.subplots(len(sets), 2, squeeze=False, gridspec_kw={"width_ratios": [2.1, 1], "wspace": .35,
                                                                       "hspace": .75})
        fig.subplots_adjust(top=1 - .95 / height)
        for k, (pair, col) in enumerate(sets):
            sub = df[(df["pair"] == pair) & (df["assignment_column"] == col)]
            _bars(axes[k, 0], sub, sections)
            _forest(axes[k, 1], sub, sections)
            axes[k, 0].set_title(_set_title(cfg, sub.iloc[0]), fontsize=8, loc="left")
            texts += [_row_text(cfg, r, labels.get(r["sample"], r["sample"])) for _, r in sub.iterrows()]
        handles = [Patch(color=sections[0][2] if sections else "0.2", label="observed: source bins (section color)")]
        handles += [Patch(color=g, label=lab) for _, lab, g in BASELINES]
        handles += [Line2D([], [], marker="o", color="0.3", markerfacecolor="0.3", linewidth=1,
                           label="difference, 95% block-bootstrap CI"),
                    Line2D([], [], marker="o", color="0.3", markerfacecolor="white", linewidth=1,
                           label="same, co-positive source bins excluded"),
                    Line2D([], [], marker="s", color="0.6", markerfacecolor="0.6", markersize=4, linewidth=1,
                           label="ROI-quantile strata alone")]
        fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, 1.0), ncol=3, fontsize=6.5,
                   frameon=False)
        common.save_figure(fig, out, NAME, dpi=cfg.get("figures", {}).get("dpi", 200))
    if exclude_self:
        all_roi = ("the manuscript baseline; with exclude_self each target bin is measured to the nearest other "
                   "target bin")
        copos = "co-positive bins, measured to the nearest other target bin"
    else:
        all_roi = "the manuscript baseline, which counts the target bins themselves, at 0 um, as within"
        copos = "co-positive bins, which sit at 0 um"
    caption = (f"Proximity tables against other baselines. Left: percentage of source bins in the marrow region of "
               f"interest whose nearest target bin lies within {_f(radius, 'g')} um (center to center; section color; "
               f"n = source bins) next to three baselines: all ROI bins ({all_roi}), ROI bins that are not targets, "
               f"and the depth-matched expectation (each source bin is given the within-radius rate of the ROI bins "
               f"that are neither source nor target in its total-UMI stratum; strata cut at the {st['n_strata']} "
               f"quantiles of total UMIs over the section's ROI bins and at the same quantiles over its source bins; "
               f"n.e., not estimable). Right: observed minus depth-matched expectation in percentage points, with 95% "
               f"spatial block-bootstrap intervals ({top['n_boot']} resamples of {st['block_um']:g} um square blocks); "
               f"open circles: the same for source bins that are not themselves targets (positivity calls with "
               f"{copos}); grey squares: the same with strata at the ROI quantiles only. Source bins are much deeper "
               f"than most ROI bins, so with ROI-quantile strata alone most of them share the top stratum, sit at its "
               f"deep end and are given the rate of shallower reference bins; that difference is biased upward "
               f"wherever depth raises the within-radius rate, as the median UMI of the reference bins weighted by "
               f"each match (matched reference, against the source median) shows. Because the all-ROI baseline "
               f"differs between sections, the largest ratio it allows (100 / baseline) differs too, and ratios are "
               f"not comparable between sections. " + "; ".join(texts) + ".")
    return [common.figure_entry(NAME, FIGURES[0][1], caption)]
