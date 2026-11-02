"""
Figure of the Figure 4i sensitivity analysis (bone_sensitivity.uncoupling): the observed/null distance ratio of every
section under each null model and call set, with its range over the depth-matching draws, and the ratio of ratios
against the reference section.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatter

import pipeline_config as pc
from bone_sensitivity import common

NULL_CSV = "sensitivity_Figure4i_null_models.csv"
SUMMARY_CSV = "sensitivity_Figure4i_summary.csv"
DECOMPOSITION_CSV = "sensitivity_Figure4i_decomposition.csv"
STABILITY_CSV = "sensitivity_Figure4i_call_stability.csv"
FIG = "Sensitivity_Figure4i_null_models"
FIGURES = [(FIG, "Sensitivity — Figure 4i under other null models")]
INPUTS = [NULL_CSV, SUMMARY_CSV, DECOMPOSITION_CSV, STABILITY_CSV]
CALL_NAMES = {"pooled_dm": "pooled, depth-matched", "pooled": "pooled, native depth", "native": "native (unpooled)"}
_FALLBACK_COLORS = ("#4C72B0", "#C44E52", "#55A868", "#8172B3", "#CCB974", "#64B5CD", "#DD8452", "#937860")


def _samples(cfg):
    """Configured sections as {id, label, color} (bone_figures._samples: a missing color takes an unused fallback)."""
    used = {str(s["color"]).lower() for s in cfg["samples"] if s.get("color")}
    spare = (c for c in _FALLBACK_COLORS if c.lower() not in used)
    return [{"id": str(s["id"]), "label": str(s.get("label") or s["id"]),
             "color": s.get("color") or next(spare, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)])}
            for i, s in enumerate(cfg["samples"])]


class _PlainLog(LogFormatter):
    """Log-axis labels as plain numbers (0.6, 1, 2) instead of 6e-01; which minor ticks are labelled is
    LogFormatter's choice."""

    def _num_to_string(self, x, vmin, vmax):
        return f"{x:.3g}"


def _model_labels(cfg):
    u = cfg["uncoupling"]
    A, B = pc.signature_label(cfg, u["a"]), pc.signature_label(cfg, u["b"])
    return {"csr_both": "both random\n(manuscript)", "fixed_b": f"{B} kept", "fixed_a": f"{A} kept",
            "habitat_fixed_b": f"{B} kept,\n{A} in habitat", "random_labels": "labels\npermuted"}


def _model_meanings(cfg, classes):
    u = cfg["uncoupling"]
    A, B = pc.signature_label(cfg, u["a"]), pc.signature_label(cfg, u["b"])
    return {"csr_both": f"both call sets placed at random over the region of interest (complete spatial randomness, "
                        f"the Figure 4i test)",
            "fixed_b": f"{B} calls kept, {A} calls placed at random",
            "fixed_a": f"{A} calls kept, {B} calls placed at random",
            "habitat_fixed_b": f"{B} calls kept, {A} calls placed at random among the bins assigned to "
                               f"{classes} signatures (plus every called bin)",
            "random_labels": f"the called positions kept, the {A}-only and {B}-only labels permuted over them "
                             f"(co-positive bins kept)"}


def _num(x):
    return pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float)


def _fmt_range(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if not len(v):
        return None
    lo, hi = f"{v.min():.2f}", f"{v.max():.2f}"
    return lo if lo == hi else f"{lo}–{hi}"


def _panel(ax, df, models, sections, value, lo, hi, man_value, empty_text):
    """Points (mean over draws) with min-max whiskers per null model and section; the manuscript value as a star."""
    width = 0.6 / max(len(sections), 1)
    drawn = False
    for j, s in enumerate(sections):
        g = df[df["sample"] == s["id"]].set_index("null_model")
        off = (j - (len(sections) - 1) / 2) * width
        for i, m in enumerate(models):
            if m not in g.index:
                continue
            r = g.loc[m]
            v, a, b = (float(pd.to_numeric(r.get(c), errors="coerce")) for c in (value, lo, hi))
            if not (np.isfinite(v) and v > 0):
                continue
            drawn = True
            if np.isfinite(a) and np.isfinite(b) and b > a > 0:
                ax.vlines(i + off, a, b, color=s["color"], lw=1.2)
            ax.plot(i + off, v, "o", ms=4.5, color=s["color"], mec="white", mew=0.5, zorder=3)
            mv = man_value.get((s["id"], m))
            if mv is not None and np.isfinite(mv) and mv > 0:
                ax.plot(i + off, mv, "*", ms=9, mfc="none", mec="black", mew=0.8, zorder=4)
    ax.axhline(1.0, ls="--", lw=0.7, color="0.4", zorder=1)
    ax.set_xticks(range(len(models)))
    ax.set_xlim(-0.6, len(models) - 0.4)
    if not drawn:
        ax.text(0.5, 0.5, empty_text, transform=ax.transAxes, ha="center", va="center", fontsize=7, color="0.35")
    else:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(_PlainLog(labelOnlyBase=False))
        ax.yaxis.set_minor_formatter(_PlainLog(labelOnlyBase=False, minor_thresholds=(2, 0.5)))
    ax.grid(axis="y", which="major", lw=0.3, color="0.88")
    ax.spines[["top", "right"]].set_visible(False)
    return drawn


def draw(results_dir, cfg, out_dir):
    """Sensitivity_Figure4i_null_models: one column per call set; top, observed/null ratio per null model and section;
    bottom, the ratio of ratios vs the reference section. Returns [figure_entry] (none without the tables)."""
    cfg = pc.merge_with_defaults(cfg)
    nulls = common.read_csv(results_dir, NULL_CSV)
    if nulls is None or nulls.empty:
        return []
    summ = common.read_csv(results_dir, SUMMARY_CSV)
    dec = common.read_csv(results_dir, DECOMPOSITION_CSV)
    stab = common.read_csv(results_dir, STABILITY_CSV)
    st = common.settings(cfg, "uncoupling")
    u, ref = cfg["uncoupling"], str(cfg["reference_sample"])
    nulls = nulls.assign(sample=nulls["sample"].astype(str), draw=_num(nulls["draw"]),
                         ratio=_num(nulls["ratio"]), rr=_num(nulls["ratio_of_ratios_vs_reference"]),
                         ok=nulls["status"].astype(str).eq("estimable"))
    present_cs = list(pd.unique(nulls["call_set"].astype(str)))
    call_sets = [c for c in st["call_sets"] if c in present_cs] + [c for c in present_cs if c not in st["call_sets"]]
    present_m = list(pd.unique(nulls["null_model"].astype(str)))
    models = [m for m in st["null_models"] if m in present_m] + [m for m in present_m if m not in st["null_models"]]
    sections = [s for s in _samples(cfg) if s["id"] in set(nulls["sample"])]
    others = [s for s in sections if s["id"] != ref]
    if summ is None or summ.empty:              # rebuild the per-draw summary the figure needs
        g = nulls[nulls["ok"]].groupby(["call_set", "null_model", "sample"])
        summ = g.agg(ratio_mean=("ratio", "mean"), ratio_min=("ratio", "min"), ratio_max=("ratio", "max"),
                     ratio_of_ratios_mean=("rr", "mean"), ratio_of_ratios_min=("rr", "min"),
                     ratio_of_ratios_max=("rr", "max"), n_draws=("ratio", "size")).reset_index()
    summ = summ.assign(sample=summ["sample"].astype(str), call_set=summ["call_set"].astype(str),
                       null_model=summ["null_model"].astype(str))
    man = nulls[(nulls["call_set"] == u["calls"]) & (nulls["draw"] == 0) & (nulls["null_model"] == "csr_both")]
    man_ratio = {(r["sample"], "csr_both"): r["ratio"] for _, r in man.iterrows() if r["ok"]}
    man_rr = {(r["sample"], "csr_both"): r["rr"] for _, r in man.iterrows() if r["ok"]}
    labels = _model_labels(cfg)

    with common.house_style():
        n = len(call_sets)
        fig = Figure(figsize=(max(2.7 * n, 3.6), 5.4))
        axes = fig.subplots(2, n, squeeze=False, sharey="row")
        for k, cs in enumerate(call_sets):
            sm = summ[summ["call_set"] == cs]
            nd = _num(sm["n_draws"]) if "n_draws" in sm else np.array([])
            nd = int(np.nanmax(nd)) if np.isfinite(nd).any() else 1
            top, bottom = axes[0, k], axes[1, k]
            _panel(top, sm, models, sections, "ratio_mean", "ratio_min", "ratio_max",
                   man_ratio if cs == u["calls"] else {}, "not estimable")
            top.set_title(CALL_NAMES.get(cs, cs) + (f"\n(mean, range over {nd} draws)" if nd > 1 else ""),
                          fontsize=7.5)
            _panel(bottom, sm[sm["sample"] != ref], models, others, "ratio_of_ratios_mean", "ratio_of_ratios_min",
                   "ratio_of_ratios_max", man_rr if cs == u["calls"] else {}, "no estimable section pair")
            top.set_xticklabels([])
            bottom.set_xticklabels([labels.get(m, m) for m in models], rotation=45, ha="right", fontsize=6.5)
        axes[0, 0].set_ylabel("observed / null mean distance", fontsize=7.5)
        ref_label = next((s["label"] for s in sections if s["id"] == ref), ref)
        axes[1, 0].set_ylabel(f"ratio of ratios vs {ref_label}", fontsize=7.5)
        handles = [Line2D([], [], marker="o", ls="", color=s["color"], mec="white", label=s["label"]) for s in sections]
        handles.append(Line2D([], [], marker="*", ls="", mfc="none", mec="black", ms=9, label="manuscript (Figure 4i)"))
        fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, fontsize=7,
                   bbox_to_anchor=(0.5, 1.02))
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        common.save_figure(fig, out_dir, FIG, dpi=cfg["figures"]["dpi"])
    return [common.figure_entry(FIG, FIGURES[0][1], _caption(cfg, st, nulls, dec, stab, sections, call_sets, models))]


def _caption(cfg, st, nulls, dec, stab, sections, call_sets, models):
    u, ref = cfg["uncoupling"], str(cfg["reference_sample"])
    A, B = pc.signature_label(cfg, u["a"]), pc.signature_label(cfg, u["b"])
    label = {s["id"]: s["label"] for s in sections}
    classes = " / ".join(st["habitat_classes"])
    meaning = _model_meanings(cfg, classes)
    perms = {}
    if "n_perm" in nulls:
        perms = {m: int(v) for m, v in nulls.groupby("null_model")["n_perm"].max().items() if np.isfinite(v)}
    man_calls = CALL_NAMES.get(u["calls"], u["calls"])
    if u["calls"] == "pooled_dm":
        man_calls += f" calls of one depth-matching draw (seed {cfg['seed']}; star)"
    else:
        man_calls += " calls (star)"
    parts = [f"Mean distance from each {A}-called bin to the nearest {B}-called bin inside the region of interest, "
             f"over its mean under random configurations (log scale; 1 = as expected under the null). Null models: "
             + "; ".join(f"{m} = {meaning.get(m, m)}" + (f" ({perms[m]} random configurations)" if m in perms else "")
                         for m in models)
             + f". The manuscript used csr_both with the {man_calls}. Call sets: "
             + ", ".join(f"{c} = {CALL_NAMES.get(c, c)}" for c in call_sets)
             + f"; pooled_dm points are means over the depth-matching draws, whiskers their range. Lower row: ratio "
               f"over the {label.get(ref, ref)} ratio."]
    man = nulls[(nulls["call_set"] == u["calls"]) & (nulls["draw"] == 0)]
    m_csr = man[(man["null_model"] == "csr_both") & man["ok"]]
    if len(m_csr):
        txt = ", ".join(f"{label.get(r['sample'], r['sample'])} {r['ratio']:.2f}" for _, r in m_csr.iterrows())
        rr = [f"{label.get(r['sample'], r['sample'])} {r['rr']:.2f}" for _, r in m_csr.iterrows()
              if r["sample"] != ref and np.isfinite(r["rr"])]
        parts.append(f"Manuscript ratios: {txt}" + (f" (ratio of ratios {', '.join(rr)})" if rr else "") + ".")
    elif len(man):
        parts.append(f"The manuscript setup was not estimable ({'; '.join(pd.unique(man['status'].astype(str)))}).")
    else:
        parts.append(f"The manuscript call set ({u['calls']}) was not computed.")
    if dec is not None and not dec.empty:
        d = dec[(dec["call_set"].astype(str) == u["calls"]) & (_num(dec["draw"]) == 0)]
        fa, fb = _num(d["factor_a_specific"]), _num(d["factor_b_pattern"])
        txt = [f"{label.get(str(s), s)} {x:.2f} x {y:.2f}" for s, x, y in zip(d["sample"], fa, fb)
               if np.isfinite(x) and np.isfinite(y)]
        if txt:
            parts.append(f"It factorises exactly into (observed / mean distance from every ROI bin to the real "
                         f"{B} calls) x (that baseline / csr null mean): {'; '.join(txt)}; the second factor measures "
                         f"how far the {B} pattern departs from random placement in the region of interest, by local "
                         f"clumping or by confinement to part of it.")
        fc, nn = _num(d["frac_b_clumped"]), _num(d["nn_self_ratio_b"])
        txt = [f"{label.get(str(s), s)} {100 * x:.0f}% / {y:.2f}" for s, x, y in zip(d["sample"], fc, nn)
               if np.isfinite(x) and np.isfinite(y)]
        if txt:
            parts.append(f"Local clumping alone (share of the {B} calls with another among their 8 neighbouring bins "
                         f"/ their nearest-neighbour spacing over that of random sets of the same size): "
                         f"{'; '.join(txt)}.")
    for model in ("fixed_b", "habitat_fixed_b"):
        g = nulls[(nulls["call_set"] == u["calls"]) & (nulls["null_model"] == model) & nulls["ok"]]
        txt = [f"{label.get(s, s)} {_fmt_range(gg['ratio'])}" for s, gg in g.groupby("sample", sort=False)
               if _fmt_range(gg["ratio"])]
        if not txt:
            continue
        pool = ""
        if model == "habitat_fixed_b" and "n_habitat_extra" in g:
            g0 = g[g["draw"] == 0]
            ex, na = _num(g0["n_habitat_extra"]), _num(g0["n_a"])
            pool = [f"{label.get(s, s)} {x:.0f} for {y:.0f} {A} calls" for s, x, y in zip(g0["sample"], ex, na)
                    if np.isfinite(x) and np.isfinite(y)]
            pool = f"; habitat bins outside the calls: {', '.join(pool)}" if pool else ""
        parts.append(f"With the {B} calls kept ({model}): {', '.join(txt)}"
                     + (" (range over draws)" if g["draw"].nunique() > 1 else "") + pool + ".")
    g = nulls[(nulls["call_set"] == u["calls"]) & (nulls["null_model"] == "csr_both") & nulls["ok"]]
    if g["draw"].nunique() > 1:
        txt = [f"{label.get(s, s)} {_fmt_range(gg['ratio'])}" for s, gg in g.groupby("sample", sort=False)]
        rr = _fmt_range(g.loc[g["sample"] != ref, "rr"])
        parts.append(f"Across {g['draw'].nunique()} depth-matching draws the csr_both ratio ranges "
                     f"{', '.join(t for t in txt if t)}" + (f", the ratio of ratios {rr}" if rr else "") + ".")
    rr_all = []
    for cs in call_sets:
        r = _fmt_range(nulls.loc[(nulls["call_set"] == cs) & nulls["ok"] & (nulls["sample"] != ref), "rr"])
        if r:
            rr_all.append(f"{cs} {r}")
    if rr_all:
        parts.append("Ratio of ratios across null models and draws: " + "; ".join(rr_all) + ".")
    if stab is not None and not stab.empty and {"region", "signature", "n_draws"} <= set(stab.columns):
        r = stab[(stab["region"].astype(str) == "roi") & (_num(stab["n_draws"]) > 1)]
        txt = []
        for sig in dict.fromkeys([u["a"], u["b"]]):
            g = r[r["signature"].astype(str) == sig]
            items = [f"{label.get(str(s), s)} {x:.0f}%" + (f" ({y:.0f}%)" if np.isfinite(y) else "")
                     for s, x, y in zip(g["sample"], _num(g["pct_unstable_of_ever_called"]),
                                        _num(g["pct_changed_vs_draw0_mean"])) if np.isfinite(x)]
            if items:
                txt.append(f"{pc.signature_label(cfg, sig)} {', '.join(items)}")
        if txt:
            parts.append(f"Pooled depth-matched calls between the {int(np.nanmax(_num(r['n_draws'])))} draws (ROI): "
                         f"share of the bins called in any draw that are not called in all (in brackets, the mean % "
                         f"of calls that differ from draw 0, 1 - Jaccard overlap): " + "; ".join(txt) + ".")
    parts.append("P values (tables) treat bins as independent units; there is one section per condition, so no "
                 "test here generalises beyond these two sections.")
    return " ".join(parts)
