"""
Figure of the Supplementary Figure S5 sensitivity tables (bone_sensitivity.clustering): the 'unresolved' share per
section and region with its lesion overlap, and the section composition and annotation of every Leiden domain; the
Harmony check goes to the caption, next to the manuscript seed-to-seed adjusted Rand index (the Leiden call's own
variability) when that table exists. Needs no scanpy.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.patches import Patch

import pipeline_config as pc
from bone_sensitivity import common

UNRESOLVED_CSV = "sensitivity_FigureS5_unresolved_by_section.csv"
ANNOTATION_CSV = "sensitivity_leiden_domain_annotation.csv"
HARMONY_CSV = "sensitivity_leiden_harmony.csv"
STABILITY_CSV = "leiden_seed_stability_adjusted_rand_index.csv"
NAME = "Sensitivity_FigureS5_domains"
FIGURES = [(NAME, "Sensitivity — Leiden domains per section and annotation")]
INPUTS = [UNRESOLVED_CSV, ANNOTATION_CSV, HARMONY_CSV]
OPTIONAL_INPUTS = [STABILITY_CSV]       # manuscript table, read when present: seed-to-seed reference of the ARI

_FALLBACK_COLORS = ("#4C72B0", "#C44E52", "#55A868", "#8172B3", "#CCB974", "#64B5CD", "#DD8452", "#937860", "#8C8C8C")
_REGIONS = (("all", "all bins", "#4d4d4d"), ("roi", "region of interest", "#7a9cc6"),
            ("outside_roi", "outside the region", "#c9c9c9"))
_LESION_COLOR, _FLAG_COLOR = "#8c510a", "#b2182b"


def _samples(cfg, ids):
    """{id: (label, color)} of the configured sections (bone_figures colors), then any other id in `ids`."""
    used = {str(s["color"]).lower() for s in cfg["samples"] if s.get("color")}
    spare = iter(c for c in _FALLBACK_COLORS if c.lower() not in used)
    out = {}
    for i, s in enumerate(cfg["samples"]):
        out[str(s["id"])] = (str(s.get("label") or s["id"]),
                             s.get("color") or next(spare, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)]))
    for i, s in enumerate(x for x in ids if x not in out):
        out[s] = (s, next(spare, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)]))
    return out


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def _p(v):
    v = _num(v)
    return f"{v:.1f}%" if np.isfinite(v) else "n/a"


def _not_available(ax, text):
    ax.text(.5, .5, text, ha="center", va="center", transform=ax.transAxes, fontsize=7, color="0.4")
    ax.set_xticks([]); ax.set_yticks([])


def _panel_unresolved(ax, un, smp):
    """Bars: % unresolved per section and region, and the % of lesion bins that are unresolved."""
    order = [s for s in smp if s in set(un["sample"])] + ["all"]
    lesion = "pct_lesion_unresolved" in un.columns
    keys = [r for r in _REGIONS if r[0] in set(un["region"])] + ([("lesion", "lesion bins", _LESION_COLOR)]
                                                                   if lesion else [])
    w = .8 / len(keys)
    for j, (key, label, color) in enumerate(keys):
        vals = []
        for s in order:
            r = un[(un["sample"] == s) & (un["region"] == ("all" if key == "lesion" else key))]
            if not len(r):
                vals.append(np.nan)
            elif key == "lesion":
                vals.append(_num(r["pct_lesion_unresolved"].iloc[0]) if _num(r["n_lesion"].iloc[0]) > 0 else np.nan)
            else:
                vals.append(_num(r["pct_unresolved"].iloc[0]))
        x = np.arange(len(order)) + (j - (len(keys) - 1) / 2) * w
        v = np.asarray(vals, dtype=float)
        ok = np.isfinite(v)
        ax.bar(x[ok], v[ok], w, color=color, label=label, hatch="//" if key == "lesion" else None,
               edgecolor="white" if key == "lesion" else None, lw=0)
        for xi, vi in zip(x[ok], v[ok]):
            ax.text(xi, vi + 1, f"{vi:.0f}", ha="center", va="bottom", fontsize=5.5)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([smp[s][0] if s in smp else "pooled" for s in order], fontsize=7)
    ax.set_ylabel("% of bins in the 'unresolved' domain", fontsize=8)
    top = np.nanmax(pd.to_numeric(un["pct_unresolved"], errors="coerce").to_numpy(float), initial=0)
    if lesion:
        top = max(top, np.nanmax(pd.to_numeric(un["pct_lesion_unresolved"], errors="coerce").to_numpy(float),
                                 initial=0))
    ax.set_ylim(0, min(100, top * 1.25 + 5) if top > 0 else 10)
    ax.legend(fontsize=6, frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("a  unresolved share per section", loc="left", fontsize=8, fontweight="bold")


def _domain_label(r):
    n = _num(r.get("n_bins"))
    return f"{r['manuscript_name']}  (n={n:,.0f})" if np.isfinite(n) else str(r["manuscript_name"])


def _panel_domains(ax, ann, smp, shares):
    """Horizontal 100% bars of each domain's section composition, labelled with its name, top programs and flags."""
    cols = [s for s in smp if f"pct_{s}" in ann.columns]
    n = len(ann)
    y = np.arange(n)[::-1]
    left = np.zeros(n)
    for s in cols:
        v = pd.to_numeric(ann[f"pct_{s}"], errors="coerce").fillna(0).to_numpy(float)
        ax.barh(y, v, left=left, color=smp[s][1], height=.7, lw=0)
        left += v
    cum = 0.0
    for s in cols[:-1]:
        cum += shares.get(s, np.nan)
        if np.isfinite(cum):
            ax.axvline(cum, color="k", lw=.6, ls="--")
    for yi, (_, r) in zip(y, ann.iterrows()):
        ax.text(102, yi + .17, str(r.get("top_programs") or "no assigned bins"), va="center", ha="left", fontsize=6,
                clip_on=False)
        nsig, umi = _num(r.get("n_significant_markers")), _num(r.get("median_umi"))
        info = f"median {umi:,.0f} UMI" if np.isfinite(umi) else ""
        if np.isfinite(nsig):
            info += (f"; {nsig:.0f} significant marker{'' if nsig == 1 else 's'}"
                     + (f", top {r['top_marker']}" if nsig > 0 and r.get("top_marker") else ""))
        ax.text(102, yi - .2, info, va="center", ha="left", fontsize=5.5, color="0.4", clip_on=False)
        if r.get("flag"):
            ax.text(-1, yi - .36, r["flag"], va="center", ha="right", fontsize=5.5, color=_FLAG_COLOR,
                    style="italic", clip_on=False)
    ax.set_yticks(y)
    ax.set_yticklabels([_domain_label(r) for _, r in ann.iterrows()], fontsize=6.5)
    ax.set_xlim(0, 100); ax.set_ylim(-.7, n - .3)
    ax.set_xlabel("% of the domain's bins from each section", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend([Patch(color=smp[s][1]) for s in cols], [smp[s][0] for s in cols], fontsize=6, frameon=False,
              loc="lower left", bbox_to_anchor=(0, 1.0), ncol=len(cols))
    ax.text(102, n - .3, "most frequent assigned programs (non-ambiguous bins)", fontsize=6, color="0.35",
            va="bottom", ha="left", clip_on=False)
    ax.set_title("b  Leiden domains: section composition and annotation", loc="left", fontsize=8, fontweight="bold",
                 pad=14)


def _caption_unresolved(un, smp):
    parts = []
    for s in [s for s in smp if s in set(un["sample"])] + ["all"]:
        r = un[un["sample"] == s].set_index("region")
        if "all" not in r.index:
            continue
        name = smp[s][0] if s in smp else "pooled"
        txt = f"{name} {_p(r.at['all', 'pct_unresolved'])}"
        if "roi" in r.index and "outside_roi" in r.index:
            txt += (f" (region of interest {_p(r.at['roi', 'pct_unresolved'])}, "
                    f"outside {_p(r.at['outside_roi', 'pct_unresolved'])})")
        parts.append(txt)
    cap = ("(a) Share of bins in the 'unresolved' Leiden domain per section and pooled over sections, over all bins, "
           "inside the region of interest and outside it (the manuscript quotes only the pooled share): "
           + "; ".join(parts) + ".")
    if "pct_lesion_unresolved" in un.columns:
        les = []
        for s in [s for s in smp if s in set(un["sample"])]:
            r = un[(un["sample"] == s) & (un["region"] == "all")]
            if len(r) and _num(r["n_lesion"].iloc[0]) > 0:
                r = r.iloc[0]
                les.append(f"{smp[s][0]}: {_p(r['pct_lesion_unresolved'])} of its {_num(r['n_lesion']):,.0f} lesion "
                           f"bins are unresolved and they make up {_p(r['pct_unresolved_lesion'])} of its unresolved "
                           f"bins")
        cap += (" Hatched: share of the cross-captured A549 lesion bins that is unresolved; "
                + ("; ".join(les) if les else "no lesion bins") + ".")
    return cap


def _caption_domains(ann, cfg, smp):
    cap = (f" (b) Section composition of each of the {len(ann)} final Leiden domains (dashed: the sections' shares of "
           f"all bins), named as in Supplementary Figure S5 (the assigned program most enriched relative to the whole "
           f"dataset), with its most frequent assigned programs among its non-ambiguous bins.")
    differ = []
    for _, r in ann.iterrows():
        name = str(r["manuscript_name"])
        if "-enriched (domain" not in name or not str(r.get("top_programs") or ""):
            continue
        first = str(r["top_programs"]).split(";")[0].strip().rsplit(" ", 1)
        top = pc.signature_label(cfg, first[0])
        if top != name.split("-enriched (domain")[0]:
            differ.append(f"{name}: {top} {first[-1]}")
    if differ:
        cap += " Named after a program that is not the most frequent one: " + "; ".join(differ) + "."
    flags = ann[ann["flag"].astype(str).str.len() > 0] if "flag" in ann.columns else ann.iloc[:0]
    if len(flags):
        txt = []
        for _, r in flags.iterrows():
            t = f"{r['manuscript_name']} ({r['flag']}"
            if "one gene" in str(r["flag"]):
                t += (f": {r['top_marker']} in {_p(r['pct_domain_detecting_top_marker'])} of its bins, "
                      f"{_p(r['pct_rest_detecting_top_marker'])} elsewhere")
            t += f"; {_num(r['n_bins']):,.0f} bins, median {_num(r['median_umi']):,.0f} UMI)"
            txt.append(t)
        cap += " Flagged domains (final Wilcoxon markers): " + "; ".join(txt) + "."
    else:
        cap += " No domain lacks a significant marker or is defined by one gene."
    return cap


def _caption_seed_reference(ari, seeds):
    """The Leiden call's own variability (manuscript seed-to-seed ARI on the same graph and call) next to the ARI of
    the partition without Harmony; '' without the seed table."""
    v = pd.to_numeric(seeds["ARI"], errors="coerce").dropna() if seeds is not None and "ARI" in seeds else []
    if not len(v):
        return ""
    lo, med, hi = float(v.min()), float(v.median()), float(v.max())
    cap = (f" For reference, Leiden runs on the manuscript embedding that differ only by seed (graph and partition "
           f"re-seeded together; {STABILITY_CSV}) agree with a median adjusted Rand index of {med:.2f} (range "
           f"{lo:.2f}-{hi:.2f}, {len(v)} pair{'' if len(v) == 1 else 's'})")
    if not np.isfinite(ari):
        return cap + "."
    if ari < lo:
        return cap + ": removing Harmony changes the partition more than changing the seed."
    if ari > hi:
        return cap + ": removing Harmony changes the partition less than changing the seed."
    return cap + (": the value without Harmony lies within the seed-to-seed range, so it cannot be told apart from "
                  "the Leiden call's own variability.")


def _caption_harmony(hm, seeds=None):
    if hm is None or not len(hm):
        return " The clustering without Harmony was not computed."
    r = {v: row for v, row in hm.set_index("variant").iterrows()}
    if "harmony" not in r or "no_harmony" not in r:
        only = next(iter(r))
        x = r[only]
        return (f" Leiden rerun on {x['embedding']} only ({_num(x['n_clusters']):.0f} clusters): no Harmony "
                f"embedding to compare with.")
    h, n = r["harmony"], r["no_harmony"]
    cap = (f" Leiden rerun with the manuscript settings on the Harmony-corrected and on the uncorrected PCA embedding "
           f"(raw partitions, before the unresolved and fragment merging): {_num(h['n_clusters']):.0f} and "
           f"{_num(n['n_clusters']):.0f} clusters, adjusted Rand index {_num(n['ari_vs_harmony']):.2f}.")
    cap += _caption_seed_reference(_num(n["ari_vs_harmony"]), seeds)
    cap += (f" Graph edges joining bins of the same section: {100 * _num(h['same_section_edge_fraction']):.1f}% with "
            f"and {100 * _num(n['same_section_edge_fraction']):.1f}% without Harmony "
            f"({100 * _num(h['same_section_edge_fraction_random']):.1f}% if edges ignored the section); clusters with "
            f">= 90% of their bins from one section: {_num(h['n_section_specific_clusters']):.0f} "
            f"({_p(h['pct_bins_in_section_specific_clusters'])} of bins) with and "
            f"{_num(n['n_section_specific_clusters']):.0f} ({_p(n['pct_bins_in_section_specific_clusters'])}) without "
            f"Harmony. Harmony corrects on the section label, and each section is its own condition.")
    rep = [row for row in (h, n) if str(row.get("reproduces_manuscript_domains")).lower() == "false"]
    if rep:
        cap += " The rerun on the manuscript embedding did not reproduce the manuscript domains exactly."
    return cap


def draw(results_dir, cfg, out_dir):
    """Sensitivity_FigureS5_domains from the clustering sensitivity tables; a missing table leaves its panel (or
    caption sentence) out. Raises FileNotFoundError when neither main table exists."""
    un = common.read_csv(results_dir, UNRESOLVED_CSV)
    ann = common.read_csv(results_dir, ANNOTATION_CSV)
    hm = common.read_csv(results_dir, HARMONY_CSV)
    seeds = common.read_csv(results_dir, STABILITY_CSV)
    if un is None and ann is None:
        raise FileNotFoundError(f"{UNRESOLVED_CSV} and {ANNOTATION_CSV} are missing")
    for df in (un, ann, hm):
        if df is not None:
            for c in ("sample", "region", "domain", "manuscript_name", "top_programs", "top_marker", "flag",
                      "suggested_label", "variant", "embedding"):
                if c in df.columns:
                    df[c] = df[c].fillna("").astype(str)
    ids = list(un["sample"]) if un is not None else []
    if ann is not None:
        ids += [c[4:] for c in ann.columns if c.startswith("pct_") and c not in
                ("pct_domain_detecting_top_marker", "pct_rest_detecting_top_marker")]
    smp = _samples(cfg, [s for s in dict.fromkeys(ids) if s != "all"])
    shares = {}
    if un is not None:
        tot = un[(un["sample"] != "all") & (un["region"] == "all")].set_index("sample")["n_bins"]
        tot = pd.to_numeric(tot, errors="coerce")
        shares = (100 * tot / tot.sum()).to_dict() if tot.sum() > 0 else {}
    n_dom = len(ann) if ann is not None else 0
    with common.house_style():
        fig = Figure(figsize=(12.5, max(3.4, 1.4 + .36 * n_dom)))
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 2.1], wspace=.55)
        ax_a, ax_b = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        if un is not None and len(un):
            _panel_unresolved(ax_a, un, smp)
        else:
            _not_available(ax_a, "unresolved share per section not available")
        if ann is not None and len(ann):
            _panel_domains(ax_b, ann, smp, shares)
        else:
            _not_available(ax_b, "domain annotation not available")
        common.save_figure(fig, out_dir, NAME, dpi=cfg["figures"]["dpi"])
    caption = _caption_unresolved(un, smp) if un is not None and len(un) else "(a) Not available."
    caption += _caption_domains(ann, cfg, smp) if ann is not None and len(ann) else " (b) Not available."
    caption += _caption_harmony(hm, seeds)
    return [common.figure_entry(NAME, FIGURES[0][1], caption)]
