#!/usr/bin/env python3
"""
bone_figures.py
===============
Draws every code-generated panel of manuscript09102026.docx from the outputs of bone_analysis.py.
Files are named by the manuscript's figure numbering and written as PNG and PDF to <results>/figures/, together with
figures/manifest.json (title, caption and files of every figure, or the error that stopped it), which the web
application reads.

  Figure 1a        experimental design schematic                                         (manuscript extra)
  Figure 2c        crops of the Space Ranger high-resolution H&E images (no drawing on the photographs)
                                                                                         (manuscript extra)
  Figure 3a, 3b    deterministic tissue assignment ("denoise" version: skeletal muscle, nerve, smooth muscle and
                   endothelium not drawn; unassigned bins as a gray background; composition renormalized over the ten
                   retained signatures)
  Figure 4a-f      osteoclast- and osteoblast-positive bins within the region of interest
  Figure 4g, 6a    positivity, native depth and depth-matched, logarithmic axes
  Figure 4h        osteoclast compartment differential expression
  Figure 4i        osteoblast-to-osteoclast distance, observed against the random-placement null
  Figure 6b        host receptor panel, depth-matched
  Supplementary S1 quality control (8 um bins), S2 Moran's I, S4 adipocyte genes, S5 Leiden domain heatmap

Not drawn by code: Figures 1b, 1d, 2a-b. Figure 1c is a chart inside the Word file. Figures 2d, 2e, 5 and S3 are
placeholders (see bone_pending_masks.py).

Inputs are the CSV outputs of bone_analysis.py (the per-bin table carries coordinates and calls) and its config.json:
no Space Ranger data and no .h5ad, except the H&E images cropped for Figure 2c. Samples, signatures, tissue classes,
panels and colors come from the configuration (pipeline_config.py); the list above names what the default
(manuscript) configuration draws, and with it the figures are those of the manuscript.

PROVENANCE
----------
  [RUN]  written by Run Zhang for the manuscript; every function of the manuscript version is [RUN]. Shamim Mollah's
         original pipeline (00_shamim_original/) has plotting steps (QC violin, PCA variance, UMAP, spatial cluster
         plot, marker plot), none of which is used in the manuscript.
  [FIX]  defect corrected on 2026-10-28 (the comment says what was wrong)
  [APP]  added for the configurable web application

Usage:  python bone_figures.py --results results [--root .] [--config cfg.json] [--manuscript-extras]
"""
from __future__ import annotations

import argparse
import functools
import json
import os
import re
import sys
import textwrap
import traceback
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.textpath import TextPath

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import pipeline_config as pc  # noqa: E402

# [RUN] house style of every panel (8 pt type, embedded TrueType fonts in PDF)
HOUSE_STYLE = {"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8, "axes.linewidth": 0.6, "figure.dpi": 150,
               "savefig.bbox": "tight"}

FLOOR = 1e-3                  # [RUN] floor of the logarithmic positivity axes (%); smaller values show no bar
DE_TOP_N = 8                  # [RUN] genes per compartment in Figure 4h
SINGLE_CELL_MIN_GENES = 100   # [RUN] dashed line of Supplementary Figure S1: a gene threshold from single-cell work
UNASSIGNED_GRAY, OUTSIDE_GRAY, INSIDE_GRAY, NOT_SIGNIFICANT_GRAY = "#e3e3e3", "#f2f2f2", "#cfcfcf", "#bdbdbd"
# [APP] colors of samples and tissue classes configured without one (signatures without one are gray)
_FALLBACK_COLORS = ("#4C72B0", "#C44E52", "#55A868", "#8172B3", "#CCB974", "#64B5CD", "#DD8452", "#937860", "#8C8C8C")

# [APP] inputs (bone_analysis.py outputs, ARCHITECTURE.md section 3) and fixed figure names
NATIVE_CSV = "Figure4g_6a_positivity_native_depth_roi.csv"
DEPTH_MATCHED_CSV = "Figure4g_6a_positivity_depth_matched.csv"
COMPOSITION_CSV = "Figure3b_program_composition_depth_matched.csv"
DE_CSV = "Figure4h_osteoclast_compartment_DE.csv"
DISTANCE_CSV = "Figure4i_distance_distributions_observed_and_null.csv"
DISTANCE_SUMMARY_CSV = "Figure4i_observed_to_null_summary.csv"
RECEPTOR_CSV = "Figure6b_receptor_panel_depth_matched.csv"
ADIPOCYTE_CSV = "FigureS4_adipocyte_genes_depth_matched.csv"
MORAN_CSV = "FigureS2_spatially_variable_genes_morans_I.csv"
DOMAIN_CSV = "FigureS5_domain_by_signature_zscores.csv"
UNRESOLVED_CSV = "FigureS5_unresolved_bins.csv"
FIG3B, FIG4AF, FIG4H, FIG4I = ("Figure3b_composition_denoise", "Figure4a-f_bone_lineage_maps",
                               "Figure4h_osteoclast_compartment_DE", "Figure4i_osteoblast_osteoclast_distance")
FIG6B, FIGS1, FIGS2, FIGS4, FIGS5 = ("Figure6b_receptor_panel", "FigureS1_quality_control",
                                     "FigureS2_spatially_variable_genes", "FigureS4_adipocyte_genes",
                                     "FigureS5_leiden_domain_heatmap")
FIG1A = "Figure1a_experimental_design"
FIG2C = (("RL-1-5", (1850, 1350, 2650, 1950), "Figure2c_top_tumor-bearing_tumor_bone_interface",
          "Figure 2c (top) — tumor-bone interface, tumor-bearing tibia"),
         ("LL-2-1", (700, 650, 1300, 1150), "Figure2c_bottom_control_chondro_osseous_junction",
          "Figure 2c (bottom) — chondro-osseous junction, control tibia"))
_TITLES = {"Figure4g_bone_lineage_positivity": "Figure 4g — bone-lineage positivity",
           "Figure6a_myeloid_positivity": "Figure 6a — myeloid positivity"}

# [APP] columns read verbatim: ids, gene names and labels ("NA" or "null" are names here, e.g. the Figure 4i source
# "null"); every other column treats the usual missing-value spellings as missing
_TEXT_COLUMNS = {"sample", "sample_id", "section", "barcode", "signature", "program", "gene", "role", "roi", "group",
                 "names", "source", "status", "note", "leiden_domain", "tissue_class", "assigned_signature",
                 "top_markers", "bin_size"}
_NA_VALUES = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>",
              "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null"]


def _house_style(func):
    """[APP] Runs a figure function under HOUSE_STYLE. The manuscript script changed matplotlib's global settings on
    import; as a library it leaves the importing program's settings alone."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with matplotlib.rc_context(HOUSE_STYLE):
            return func(*args, **kwargs)
    return wrapper


# ============================================================================
# HELPERS
# ============================================================================

def lbl(cfg, name):
    """[RUN] Signature name as drawn ('_' -> ' ', neutrophil -> granulocyte).
    [APP] the signature's label in the configuration."""
    return pc.signature_label(cfg, name)


def _slug(cfg, name):
    """[APP] Suffix of a signature's columns in the per-bin table: its label with spaces as underscores."""
    return lbl(cfg, name).replace(" ", "_")


def _sig_color(cfg, name):
    return (cfg["signatures"].get(name) or {}).get("color") or "#7F7F7F"


def _samples(cfg):
    """[APP] The configured samples in order, as {id, label, color}; a missing label is the id, a missing color the
    next fallback color that no other sample uses."""
    used = {str(s["color"]).lower() for s in cfg["samples"] if s.get("color")}
    spare = (c for c in _FALLBACK_COLORS if c.lower() not in used)
    out = []
    for i, s in enumerate(cfg["samples"]):
        sid = str(s["id"])
        color = s.get("color") or next(spare, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)])
        out.append({"id": sid, "label": str(s.get("label") or sid), "color": color})
    return out


def _reference(cfg):
    """[APP] The reference sample of fold changes and ratios, or None when it is not configured."""
    return next((s for s in _samples(cfg) if s["id"] == str(cfg.get("reference_sample"))), None)


def _signature_names(cfg):
    """[APP] Maps the ways a table may name a signature (name, label, label with underscores) to its name."""
    names = {n: n for n in cfg["signatures"]}
    for n in cfg["signatures"]:
        label = lbl(cfg, n)
        names.setdefault(label, n)
        names.setdefault(label.replace(" ", "_"), n)
    return names


def _bin_tag(bin_size):
    """[APP] '016um' -> '16um', the bin size as it appears in output file names."""
    return f"{int(pc.bin_um(bin_size))}um"


def _bin_text(bin_size):
    """[APP] '016um' -> '16 µm', the bin size in captions."""
    return f"{pc.bin_um(bin_size):g} µm"


def _per_bin_csv(cfg):
    return f"per_bin_signature_calls_{_bin_tag(cfg['bin_size'])}.csv.gz"


def _dpi(cfg, key, default):
    try:
        v = int(cfg["figures"][key])
    except (KeyError, TypeError, ValueError):
        return default
    return v if 20 <= v <= 1200 else default


def _safe_name(name):
    """[APP] A configured figure name as a file name (letters, digits, '_', '-', '.')."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name or "")).strip("._") or "positivity_panel"


def _title_from_name(name):
    """[APP] 'Figure4g_bone_lineage_positivity' -> 'Figure 4g — bone-lineage positivity' (configured panels)."""
    if name in _TITLES:
        return _TITLES[name]
    m = re.match(r"^Figure(S?)([0-9]+[A-Za-z]*(?:-[A-Za-z]+)?)_(.+)$", name)
    if not m:
        return name.replace("_", " ")
    return f"{'Supplementary Figure S' if m.group(1) else 'Figure '}{m.group(2)} — {m.group(3).replace('_', ' ')}"


def _figure3a_name(index, label):
    return f"Figure3a_{index}_{pc.slugify(label)}_tissue_classes_denoise"


def _group(name):
    m = re.match(r"^Figure(S?)(\d+)", name)
    if not m:
        return "Other figures"
    return "Supplementary figures" if m.group(1) else f"Figure {int(m.group(2))}"


def _entry(name, title, caption):
    """[APP] One figures/manifest.json entry of a figure that was drawn."""
    return {"name": name, "group": _group(name), "title": title, "caption": caption, "png": f"{name}.png",
            "pdf": f"{name}.pdf", "error": None}


def _error_entry(name, title, message):
    return {"name": name, "group": _group(name), "title": title, "caption": "", "png": None, "pdf": None,
            "error": message}


def _join(items):
    items = [str(i) for i in items]
    if len(items) < 2:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def _count(v):
    """[APP] A bin count as text (integers without decimals, thousands separated); '' when missing."""
    v = float(v)
    return f"{int(round(v)):,}" if np.isfinite(v) else ""


def _mean_count(v):
    """[RUN] Mean positive bins over the depth-matching replicates: one decimal below 100.
    [FIX] zero printed as '0.0' and a mean below 0.05 bins also as '0.0'; they print '0' and '<0.1'."""
    v = float(v)
    if not np.isfinite(v):
        return ""
    if v == 0:
        return "0"
    if v < .05:
        return "<0.1"
    return f"{v:,.1f}" if v < 100 else f"{v:,.0f}"


def _col(df, name):
    """[APP] A column as floats (NaN where missing or not numeric; all NaN when the column is absent)."""
    if name in df.columns:
        return pd.to_numeric(df[name], errors="coerce").to_numpy(float)
    return np.full(len(df), np.nan)


def _as_bool(s):
    """[APP] A True/False column of a CSV as booleans (missing = False)."""
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).to_numpy(bool)
    if pd.api.types.is_numeric_dtype(s):
        return s.fillna(0).to_numpy() != 0
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "1.0", "yes", "t"]).to_numpy()


def _wilson(k, n, z=1.959963984540054):
    """[APP] Wilson 95% interval in percent, for tables that lack theirs (NaN where n is 0 or missing)."""
    k, n = np.asarray(k, float), np.asarray(n, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = k / n
        den = 1 + z ** 2 / n
        centre = (p + z ** 2 / (2 * n)) / den
        half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return 100 * np.clip(centre - half, 0, 1), 100 * np.clip(centre + half, 0, 1)


def _per_sample(df, col):
    """[APP] The first finite value of `col` per sample (e.g. n_bins, the same on every row of a sample)."""
    out = {}
    if col in df.columns and "sample" in df.columns:
        for sid, v in zip(df["sample"], pd.to_numeric(df[col], errors="coerce")):
            if sid not in out and np.isfinite(v):
                out[sid] = float(v)
    return out


def _sizes_text(cfg, values, unit="bins"):
    """[APP] 'control 12,345; tumor-bearing 23,456' in sample order (samples without a value left out)."""
    parts = [f"{s['label']} {_count(values[s['id']])}" for s in _samples(cfg) if s["id"] in values]
    return f"{'; '.join(parts)} {unit}" if parts else ""


def _fold_text(cfg, values, items, name_of, what="Fold"):
    """[APP] Fold of every other sample over the reference sample, per item (values: sample -> item -> value)."""
    ref = _reference(cfg)
    others = [s for s in _samples(cfg) if ref and s["id"] != ref["id"] and s["id"] in values]
    if not ref or ref["id"] not in values or not others:
        return ""
    parts = []
    for s in others:
        folds = []
        for it in items:
            r, v = values[ref["id"]].get(it, np.nan), values[s["id"]].get(it, np.nan)
            folds.append(f"{name_of(it)} {v / r:.2f}" if np.isfinite(r) and np.isfinite(v) and r > 0
                         else f"{name_of(it)} n/a")
        parts.append(f"{s['label']}: {', '.join(folds)}")
    return f" {what} vs {ref['label']}: {'; '.join(parts)}."


def _read_table(path, usecols=None, nrows=None):
    """[APP] A CSV output of bone_analysis.py (plain or .gz). Text columns are read verbatim (_TEXT_COLUMNS); an empty
    file gives an empty table instead of raising EmptyDataError."""
    try:
        header = list(pd.read_csv(path, nrows=0).columns)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    cols = header if usecols is None else [c for c in header if c in set(usecols)]
    if not cols:
        return pd.DataFrame()
    return pd.read_csv(path, usecols=cols, nrows=nrows, dtype={c: str for c in cols if c in _TEXT_COLUMNS},
                       keep_default_na=False, na_values={c: _NA_VALUES for c in cols if c not in _TEXT_COLUMNS})


def _read_json(path):
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def save(fig, out, name, dpi=200):
    """[RUN] Writes one figure as PNG and PDF.
    [FIX] the folder is an argument: it was the module global OUT, set only when the file ran as a script, so
    importing the module and calling a figure function failed."""
    out = Path(out)
    fig.savefig(out / f"{name}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")


def section_xy_mm(pb, sample_id):
    """[RUN] Bin centres of one section in mm.
    [FIX] read from the per-bin table (x_um, y_um, already scaled by the section's own microns per pixel) instead of
    the .h5ad with microns per pixel hard-coded for the two manuscript sections. Returns (rows, x, y)."""
    sub = pb[pb["sample_id"] == sample_id]
    x, y = _col(sub, "x_um"), _col(sub, "y_um")
    ok = np.isfinite(x) & np.isfinite(y)
    return sub[ok], x[ok] / 1000, y[ok] / 1000


def _span(v, minimum):
    return max(float(np.max(v) - np.min(v)), minimum) if len(v) else minimum


def _scalebar_mm(width_mm):
    """[APP] Scale-bar length: 1 mm, or for a section narrower than about 2.5 mm the longest 1-2-5 length that fits
    in 40% of its width."""
    if not np.isfinite(width_mm) or width_mm >= 2.5:
        return 1.0
    return next((v for v in (0.5, 0.2, 0.1, 0.05, 0.02) if v <= 0.4 * width_mm), 0.01)


def _length_text(mm):
    return f"{mm:g} mm" if mm >= 1 else f"{mm * 1000:g} µm"


def scalebar(ax, xs, ys, fs=8, lw=3, length_mm=None):
    """[RUN] Scale bar at the lower left of a section (y axis pointing down). Returns its label.
    [FIX] it was always 1 mm, longer than a section narrower than about 2.5 mm allows: the length is now 1 mm or a
    shorter 1-2-5 length (or `length_mm`, shared by the panels of one figure). Figure 3a drew its own bar with
    projecting caps, which made it longer than the length written under it; every bar now has butt caps."""
    xr, yr = _span(xs, 1e-3), _span(ys, 1e-3)
    v = min(xr, yr)                        # vertical offsets: the section height, or its width when it is tall
    length = length_mm or _scalebar_mm(xr)
    x0 = xs.min() + .03 * xr
    y0 = ys.max() - .02 * v
    ax.plot([x0, x0 + length], [y0, y0], "k-", lw=lw, solid_capstyle="butt")
    ax.text(x0 + length / 2, y0 - .03 * v, _length_text(length), ha="center", va="bottom", fontsize=fs)
    return _length_text(length)


def _marker_area(s0, pitch_pt, rel):
    """[APP] Scatter marker area (points^2) of the maps: the manuscript's s0, or, where bins lie farther apart on the
    page than in the manuscript (small sections, coarse bins; pitch_pt = bin spacing in points), a marker `rel` bin
    spacings across. The rel values are chosen so that s0 wins at the manuscript's scale (16 um bins of sections about
    5 mm or more across)."""
    return max(s0, (rel * pitch_pt) ** 2)


def _letter(i):
    """[APP] Panel letter: 0 -> a, 25 -> z, 26 -> aa."""
    s, i = "", i + 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(97 + r) + s
    return s


def _chars(room_in, fontsize):
    """[APP] Characters of text at `fontsize` that fit `room_in` inches (average glyph width 0.6 em)."""
    return max(1, int(room_in * 72 / (0.6 * fontsize)))


def _text_in(text, size):
    """[APP] Width of one line of text in inches at `size` points (current font)."""
    return TextPath((0, 0), text, size=size).get_extents().width / 72 if text.strip() else 0.0


def _fit_text(text, room_in, size, smallest=6.0):
    """[APP] (text, font size) of a label that must fit `room_in` inches: one line at the largest size from `size`
    down to `smallest` that fits, else wrapped at spaces and hyphens at `smallest`."""
    fs = size
    while fs > smallest and _text_in(text, fs) > room_in:
        fs -= .5
    if _text_in(text, fs) <= room_in:
        return text, fs
    lines = textwrap.wrap(text, _chars(room_in, fs), break_long_words=False, break_on_hyphens=True)
    return "\n".join(lines or [text]), fs


def _label_bars(ax, labels, n_samples, spacing_in, log, axes_h_in, ymin, ymax):
    """[APP] Writes count labels (x, y, text) above bars at 6 pt: horizontal as in the manuscript, or rotated by 90
    degrees when there are three or more samples and a label is wider than the spacing of neighboring bars. Returns
    the upper y limit (at least ymax) at which every label is inside the axes (axes_h_in: approximate axes height)."""
    if not labels:
        return ymax
    widths = [_text_in(t, 6) for _, _, t in labels]
    rotate = n_samples >= 3 and max(widths) > spacing_in
    extent = widths if rotate else [6 * 1.3 / 72] * len(labels)          # inches above the anchor
    top = ymax
    for _ in range(5):
        if log:
            per_in = (np.log10(top) - np.log10(ymin)) / axes_h_in
            need = 10 ** (max(np.log10(y) + h * per_in for (_, y, _), h in zip(labels, extent)) + .03)
        else:
            per_in = (top - ymin) / axes_h_in
            need = 1.02 * max(y + h * per_in for (_, y, _), h in zip(labels, extent))
        if need <= top:
            break
        top = need
    for x, y, t in labels:
        ax.text(x, y, t, ha="center", va="bottom" if rotate else "baseline", fontsize=6, rotation=90 if rotate else 0)
    return top


def _depth_text(cfg, summary):
    """[FIX] How the depth matching worked, for captions: "to N UMI per bin" was wrong for bins below the target (kept
    as they are) and for scoring.depth_matching = "quantile" (per-bin targets)."""
    target = (summary or {}).get("depth_match_target_umi")
    if cfg["scoring"].get("depth_matching") == "quantile":
        return (" (each bin of the deeper sections downsampled to the UMI count at the same quantile of the shallowest "
                "section's per-bin distribution)")
    return f" (bins above {target} UMI, the median of the shallowest section, downsampled to it)" if target else ""


def _tick_labels(ax, labels, fontsize, slot_in):
    """[APP] x tick labels wrapped to the width of their slot (slot_in: inches per tick); slanted when a single word
    is too long for it."""
    n = _chars(slot_in, fontsize)
    wrapped = [textwrap.wrap(str(t), n, break_long_words=False, break_on_hyphens=True) or [""] for t in labels]
    ax.set_xticks(range(len(labels)))
    if any(len(line) > n for w in wrapped for line in w) or any(len(w) > 3 for w in wrapped):
        ax.set_xticklabels(labels, fontsize=fontsize, rotation=30, ha="right")
    else:
        ax.set_xticklabels(["\n".join(w) for w in wrapped], fontsize=fontsize)


# ============================================================================
# FIGURE 1a, 2c (manuscript extras; never drawn by the web application)
# ============================================================================

@_house_style
def figure1a(out, dpi=250):
    """[RUN] Timeline of the experiment.
    [APP] manuscript extra (make_figures(..., manuscript_extras=True)): it describes the manuscript's animal experiment
    and is not drawn for other configurations."""
    fig = Figure(figsize=(8.5, 2.2))
    ax = fig.subplots()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.annotate("", xy=(.965, .62), xytext=(.035, .62),
                arrowprops=dict(arrowstyle="-|>", lw=2.5, color="k", mutation_scale=14))
    for x, day, lab in [(.06, "day −1", "µCT\n(baseline)"), (.26, "day 0", "A549-luc cells,\nintratibial, right tibia"),
                        (.455, "day 7", "BLI"), (.67, "day 20", "µCT"), (.885, "day 21", "BLI; endpoint")]:
        ax.plot([x, x], [.57, .67], "k-", lw=2.5)
        ax.text(x, .72, day, ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.text(x, .52, lab, ha="center", va="top", fontsize=8.5)
    ax.text(.5, .06, "Both tibiae of one animal profiled by Visium HD on a single capture slide", ha="center",
            va="center", fontsize=9, style="italic")
    save(fig, out, FIG1A, dpi=dpi)
    return [_entry(FIG1A, "Figure 1a — experimental design",
                   "Experimental design: baseline µCT on day −1; A549-luc cells injected into the right tibia on day 0; "
                   "bioluminescence imaging (BLI) on days 7 and 21; µCT on day 20; endpoint on day 21. Both tibiae of "
                   "the animal were profiled by Visium HD on a single capture slide.")]


def _hires_image(root, cfg, sample_id):
    """[APP] tissue_hires_image.png of a sample under `root`, through the configured sample path (any bin size: the
    high-resolution image is the same), else the manuscript's folder layout."""
    rcfg = {**cfg, "data_root": str(root)}
    sample = pc.sample_by_id(cfg, sample_id)
    candidates = []
    for b in dict.fromkeys(["008um", cfg["bin_size"], *pc.available_bin_sizes(rcfg, sample)]):
        d = pc.resolve_bin_dir(rcfg, sample, b)
        if d is not None:
            candidates.append(d / "spatial" / "tissue_hires_image.png")
    candidates.append(Path(root) / f"SR007144_{sample_id}/outs/binned_outputs/square_008um/spatial/"
                                    "tissue_hires_image.png")
    for c in candidates:
        if c.is_file():
            return c
    raise FileNotFoundError(f"no tissue_hires_image.png of sample {sample_id} under {root}")


def figure2c(root, cfg, out):
    """[RUN] Plain crops of the Space Ranger high-resolution H&E images (pixel boxes in Manuscript/figures/README.md);
    arrows are to be placed by the pathologist, not drawn here.
    [APP] manuscript extra: needs `root` and the manuscript sample ids; the images are found through the configured
    sample paths; a PDF is written next to each PNG; a box outside the image raises instead of being padded with
    black."""
    from PIL import Image
    out = Path(out)
    entries = []
    for sid, box, name, title in FIG2C:
        with Image.open(_hires_image(root, cfg, sid)) as im:
            if box[2] > im.width or box[3] > im.height:
                raise ValueError(f"crop box {box} lies outside the {im.width} x {im.height} image of {sid}")
            crop = im.convert("RGB").crop(box)
        crop.save(out / f"{name}.png")
        crop.save(out / f"{name}.pdf", "PDF", resolution=300)
        entries.append(_entry(name, title, f"Crop {box} (left, top, right, bottom; high-resolution image pixels) of "
                                           f"the Space Ranger H&E image of section {sid}, unmodified."))
    return entries


# ============================================================================
# FIGURE 3 (deterministic tissue assignment)
# ============================================================================

def _class_members(cfg):
    """[APP] Drawn tissue classes with the competing signatures (cfg['programs']) of each, and the signatures whose
    class is not drawn."""
    drawn = [c["name"] for c in cfg["classes"] if c.get("name")]
    members = {c: [p for p in cfg["programs"] if (cfg["signatures"].get(p) or {}).get("class") == c] for c in drawn}
    removed = [p for p in cfg["programs"] if (cfg["signatures"].get(p) or {}).get("class") not in drawn]
    return members, removed


@_house_style
def figure3a(pb, cfg, out, sample_id, index):
    """[RUN] Tissue class per bin (assignment 'program', native depth, pooled). Skeletal muscle, nerve, smooth muscle
    and endothelium are not drawn; ambiguous bins are a gray background without a legend entry.
    [FIX] read from the per-bin table (tissue_class) instead of the .h5ad, whose backed handle was never closed; the
    two sections, their names and microns per pixel were hard-coded.
    [APP] one figure per configured sample (index = 1-based position); classes, colors, legend order and the larger
    markers of emphasized classes from cfg['classes']. Sections taller than wide are drawn at most 8 in high."""
    smp = next((s for s in _samples(cfg) if s["id"] == sample_id), None)
    if smp is None:
        raise ValueError(f"sample {sample_id} is not configured")
    sub, xs, ys = section_xy_mm(pb, sample_id)
    if not len(sub):
        raise ValueError(f"the per-bin table has no bins of sample {sample_id}")
    cls = sub["tissue_class"].fillna("").astype(str).to_numpy()
    classes = [c for c in cfg["classes"] if c.get("name")]
    if not np.isin(cls, ["unassigned"] + [c["name"] for c in classes]).any():
        raise ValueError(f"no bin of sample {sample_id} has a tissue class (tissue_class is empty)")
    bin_mm = pc.bin_um(cfg["bin_size"]) / 1000
    xr, yr = _span(xs, bin_mm), _span(ys, bin_mm)
    w, h = 6.5, 6.5 * yr / xr
    if h > 8.0:
        w, h = max(8.0 * xr / yr, 2.5), 8.0
    pitch = bin_mm * min(w / xr, h / yr) * 72
    fig = Figure(figsize=(w, h + .6))
    ax = fig.subplots()
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
    bg = cls == "unassigned"
    ax.scatter(xs[bg], ys[bg], s=_marker_area(.6, pitch, .51), color=UNASSIGNED_GRAY, lw=0, rasterized=True)
    shares, handles = [], []
    for i, c in enumerate(classes):
        m = cls == c["name"]
        s0, rel = (1.4, .78) if c.get("emphasis") else (.6, .51)
        color = c.get("color") or _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)]
        ax.scatter(xs[m], ys[m], s=_marker_area(s0, pitch, rel), color=color, lw=0, rasterized=True)
        handles.append(Line2D([], [], marker="o", ls="", color=color, mew=0, ms=8 * np.sqrt(s0)))   # markerscale 8
        shares.append(100 * m.mean())
    bar = scalebar(ax, xs, ys, fs=8, lw=3)
    ax.set_xlim(xs.min() - .05, xs.min() + xr + .05); ax.set_ylim(ys.min() + yr + .05, ys.min() - .05)
    if classes:
        ax.legend(handles, [c["name"] for c in classes], fontsize=7.5, ncol=min(4 if w >= 5 else 2, len(classes)),
                  loc="upper center", bbox_to_anchor=(.5, -.01), frameon=False)
    name = _figure3a_name(index, smp["label"])
    save(fig, out, name, dpi=_dpi(cfg, "map_dpi", 250))
    members, removed = _class_members(cfg)
    parts = [f"{c['name']} ({', '.join(lbl(cfg, p) for p in members[c['name']]) or 'no signature'}; {pct:.1f}%)"
             for c, pct in zip(classes, shares)]
    caption = (f"Deterministic tissue assignment of the {smp['label']} section ({len(sub):,} "
               f"{_bin_text(cfg['bin_size'])} bins after quality control). Each bin takes the signature with the "
               f"highest spatially pooled score (the bin and its {cfg['scoring']['pool_k']} nearest bins, native depth) "
               f"among those that clear their null threshold with at least {cfg['scoring']['min_core']} core genes "
               f"detected, and is colored by that signature's tissue class: {'; '.join(parts)}. ")
    if removed:
        caption += f"Bins assigned to {_join(lbl(cfg, p) for p in removed)} are not drawn. "
    caption += (f"Unassigned bins (no signature clears its threshold; {100 * bg.mean():.1f}%) form the gray "
                f"background. Scale bar {bar}.")
    return [_entry(name, f"Figure 3a — tissue classes, {smp['label']}", caption)]


@_house_style
def figure3b(comp_dm, cfg, out, summary=None):
    """[RUN] Composition after depth matching (assignment 'program_dm'), renormalized over the ten retained
    signatures; ambiguous, skeletal muscle, nerve, smooth muscle and endothelium excluded (counts written alongside).
    [FIX] legend entries were attached while drawing the first sample only, and only for signatures present in it, so
    a signature absent from the control had no legend entry; a section without retained bins divided by zero.
    [APP] samples, retained signatures (cfg['figures']['composition_signatures']) and colors from the
    configuration."""
    out = Path(out)
    keep = list(dict.fromkeys(cfg["figures"]["composition_signatures"]))
    if not keep:
        raise ValueError("figures.composition_signatures is empty")
    names = _signature_names(cfg)
    prog = comp_dm["program"].astype(str).map(names)
    smp = _samples(cfg)
    rows, ambiguous = [], _per_sample(comp_dm, "sample_ambiguous_pct")
    for s in smp:
        mine = comp_dm["sample"] == s["id"]
        n_bins = pd.to_numeric(comp_dm.loc[mine, "n_bins"], errors="coerce").fillna(0)
        counts = n_bins.groupby(prog[mine]).sum()
        tot = int(sum(counts.get(p, 0) for p in keep))
        for p in keep:
            n = int(counts.get(p, 0))
            rows.append({"sample_id": s["id"], "section": s["label"], "signature": lbl(cfg, p),
                         "bins_depth_matched": n, "retained_bins": tot,
                         "percent_of_retained": round(100 * n / tot, 2) if tot else np.nan})
    comp = pd.DataFrame(rows)
    comp.to_csv(out / f"{FIG3B}.csv", index=False)
    fig = Figure(figsize=(2.7 + .55 * len(smp), 3.6))
    ax = fig.subplots()
    for i, s in enumerate(smp):
        c = comp[comp["sample_id"] == s["id"]]
        if not int(c["retained_bins"].iloc[0]):
            ax.text(i, 2, "no retained bins", rotation=90, ha="center", va="bottom", fontsize=6.5, color="0.4")
            continue
        b0 = 0.0
        for p, v in zip(keep, c["percent_of_retained"]):
            if v > 0:
                ax.bar(i, v, bottom=b0, color=_sig_color(cfg, p), width=.6)
                b0 += v
    _tick_labels(ax, [s["label"] for s in smp], 8, (fig.get_figwidth() - 1.2) / len(smp))
    ax.set_xlim(-.6, len(smp) - .4)
    ax.set_ylabel("% of retained bins", fontsize=8); ax.set_ylim(0, 100)
    handles = [Patch(color=_sig_color(cfg, p), label=lbl(cfg, p)) for p in keep]
    fig.legend(handles[::-1], [h.get_label() for h in handles[::-1]], fontsize=6.5, loc="center left",
               bbox_to_anchor=(1.0, .5), frameon=False)
    fig.tight_layout()
    save(fig, out, FIG3B, dpi=_dpi(cfg, "map_dpi", 250))
    excluded = [p for p in cfg["programs"] if p not in keep]
    totals = {r["sample_id"]: r["retained_bins"] for r in rows}
    caption = (f"Tissue composition after depth matching{_depth_text(cfg, summary)} "
               f"(assignment on spatially pooled, depth-matched scores), as a percentage of the bins assigned to one "
               f"of the {len(keep)} retained signatures")
    caption += f"; unassigned bins and bins assigned to {_join(lbl(cfg, p) for p in excluded)} are excluded. " \
        if excluded else "; unassigned bins are excluded. "
    caption += f"Retained bins: {_sizes_text(cfg, totals)}."
    if ambiguous:
        caption += " Unassigned: " + "; ".join(f"{s['label']} {ambiguous[s['id']]:.1f}%" for s in smp
                                               if s["id"] in ambiguous) + " of all bins."
    caption += f" Counts and percentages: {FIG3B}.csv."
    return [_entry(FIG3B, "Figure 3b — tissue composition after depth matching", caption)]


# ============================================================================
# FIGURE 4a-f (positive bins within the region of interest)
# ============================================================================

@_house_style
def figure4af(pb, cfg, out):
    """[RUN] Osteoclast-positive (blue), osteoblast-positive (red) and both (black) bins within the ROI; one spatial
    scale for all six panels.
    [FIX] read from the per-bin table instead of the .h5ad, whose backed handle was never closed; the two sections,
    their names, sizes and microns per pixel were hard-coded, and the 1 mm scale bar did not fit narrow sections.
    [APP] one row per configured sample; the two signatures and their colors from cfg['figures']['lineage_maps'];
    without an in_region_of_interest column every bin counts as inside."""
    a, b = cfg["figures"]["lineage_maps"]["a"], cfg["figures"]["lineage_maps"]["b"]
    ca, cb = f"positive_{_slug(cfg, a)}", f"positive_{_slug(cfg, b)}"
    for c in (ca, cb):
        if c not in pb.columns:
            raise ValueError(f"the per-bin table has no column {c}")
    col_a, col_b = _sig_color(cfg, a), _sig_color(cfg, b)
    bin_mm = pc.bin_um(cfg["bin_size"]) / 1000
    rows = []
    for s in _samples(cfg):
        sub, xs, ys = section_xy_mm(pb, s["id"])
        if len(sub):
            rows.append((s, sub, xs, ys, _span(xs, bin_mm), _span(ys, bin_mm)))
    if not rows:
        raise ValueError("the per-bin table has no bins of any configured sample")
    has_roi = "in_region_of_interest" in pb.columns
    labw, gap, head, vgap, foot = .32, .06, .25, .12, .45
    colw = (6.5 - labw - 2 * gap) / 3
    scale = min(colw / max(r[4] for r in rows), 9.0 / sum(r[5] for r in rows))     # inches per mm, every panel
    H = head + sum(r[5] * scale for r in rows) + vgap * (len(rows) - 1) + foot
    pitch = bin_mm * scale * 72
    s_bin, s_pos, s_both = (_marker_area(s0, pitch, rel) for s0, rel in ((.08, .6), (1.2, 1.0), (1.6, 1.15)))
    bar_mm = _scalebar_mm(min(r[4] for r in rows))
    fig = Figure(figsize=(6.5, H))
    titles = [lbl(cfg, a), lbl(cfg, b), "overlay"]
    y_top, n_panel, counts = H - head, 0, []
    for r, (s, sub, xs, ys, xr, yr) in enumerate(rows):
        roi = _as_bool(sub["in_region_of_interest"]) if has_roi else np.ones(len(sub), bool)
        pa, pb_ = _as_bool(sub[ca]) & roi, _as_bool(sub[cb]) & roi
        counts.append((s, int(roi.sum()), int(pa.sum()), int(pb_.sum()), int((pa & pb_).sum())))
        rowh = yr * scale
        y0 = y_top - rowh
        for jx, col in enumerate(titles):
            ax = fig.add_axes([(labw + jx * (colw + gap)) / 6.5, y0 / H, xr * scale / 6.5, rowh / H])
            ax.set_xlim(xs.min(), xs.min() + xr); ax.set_ylim(ys.min() + yr, ys.min()); ax.axis("off")
            ax.scatter(xs[~roi], ys[~roi], s=s_bin, color=OUTSIDE_GRAY, lw=0, rasterized=True)
            ax.scatter(xs[roi], ys[roi], s=s_bin, color=INSIDE_GRAY, lw=0, rasterized=True)
            overlay = jx == 2
            if jx in (0, 2):
                m = pa & ~(pb_ & overlay)
                ax.scatter(xs[m], ys[m], s=s_pos, color=col_a, lw=0, rasterized=True)
            if jx in (1, 2):
                m = pb_ & ~(pa & overlay)
                ax.scatter(xs[m], ys[m], s=s_pos, color=col_b, lw=0, rasterized=True)
            if overlay:
                ax.scatter(xs[pa & pb_], ys[pa & pb_], s=s_both, color="black", lw=0, rasterized=True)
            scalebar(ax, xs, ys, fs=6, lw=2, length_mm=bar_mm)
            ax.text(-.02, 1.0, _letter(n_panel), transform=ax.transAxes, fontsize=10, fontweight="bold", va="top",
                    ha="right")
            n_panel += 1
            if r == 0:
                text, fs = _fit_text(col, colw, 9)
                fig.text((labw + jx * (colw + gap) + colw / 2) / 6.5, (H - head / 2) / H, text, ha="center",
                         va="center", fontsize=fs)
        text, fs = _fit_text(s["label"], rowh + vgap, 9)       # may reach into the gaps, not into the next row
        fig.text(.12 / 6.5, (y0 + rowh / 2) / H, text, rotation=90, ha="center", va="center", fontsize=fs,
                 linespacing=1.0)
        y_top = y0 - vgap
    hs = [Line2D([], [], marker="o", ls="", color=c, ms=5) for c in (col_a, col_b, "black", INSIDE_GRAY, OUTSIDE_GRAY)]
    fig.legend(hs, [f"{titles[0]}-positive", f"{titles[1]}-positive", "positive for both",
                    "other unit in region of interest", "outside region of interest"], loc="lower center", ncol=3,
               fontsize=7, frameon=False, bbox_to_anchor=(.5, 0))
    save(fig, out, FIG4AF, dpi=_dpi(cfg, "map_dpi", 250))
    A, B = titles[0], titles[1]
    caption = (f"Bins positive for {A}, for {B} and for both (black) at native depth, within the region of interest "
               f"of each section (rows); columns: {A}, {B}, overlay. Gray: other bins of the region of interest; "
               f"light gray: bins outside it. One spatial scale for all panels; scale bar {_length_text(bar_mm)}. ")
    caption += " ".join(f"{s['label']}: {n_roi:,} bins in the region of interest, {na:,} {A}-positive, {nb:,} "
                        f"{B}-positive, {nab:,} both." for s, n_roi, na, nb, nab in counts)
    if not has_roi:
        caption += " The per-bin table has no region-of-interest flag; every bin is counted as inside."
    absent = [s["label"] for s in _samples(cfg) if s["id"] not in {c[0]["id"] for c in counts}]
    if absent:
        caption += f" No bins of {_join(absent)} in the per-bin table."
    return [_entry(FIG4AF, f"Figures 4a-f — {A}- and {B}-positive bins", caption)]


# ============================================================================
# FIGURES 4g, 6a (positivity, native depth and depth-matched)
# ============================================================================

def _by_signature(t, cfg, sid, sigs):
    """[APP] Rows of one sample indexed by signature name, in the order `sigs` (missing ones as NaN rows)."""
    names = _signature_names(cfg)
    sub = t[t["sample"] == sid].copy()
    sub["_sig"] = sub["signature"].astype(str).map(names)
    return sub.dropna(subset=["_sig"]).drop_duplicates("_sig").set_index("_sig").reindex(sigs)


@_house_style
def positivity_panel(t1, t2, cfg, out, name, sigs, summary=None):
    """[RUN] Left: positive bins at native depth with Wilson intervals. Right: depth-matched mean +/- SD over five
    replicates. Values in percent; anything below 0.001% (including zero) sits at the axis floor and shows no bar.
    [FIX] .loc[sigs] raised KeyError for a signature absent from a table (left out now, and named in the caption);
    whiskers could be slightly negative after rounding, which matplotlib refuses (every whisker is clipped at 0); the
    depth-matched SD whiskers were centred on the floor-clipped bar instead of the mean (now from max(mean - SD,
    floor / 2) to max(mean + SD, floor)), and a single replicate (SD undefined) broke them; the depth-matched labels
    were pct_positive_mean * n_bins of the native table instead of the replicates' mean count n_positive_mean.
    [APP] any number of samples (bars side by side, colors and labels from cfg['samples']); t2 None draws the native
    panel alone; Wilson intervals computed when the table has none."""
    in_table = set(t1["signature"].astype(str).map(_signature_names(cfg)).dropna())
    present = [s for s in dict.fromkeys(sigs) if s in in_table]
    if not present:
        raise ValueError(f"none of {', '.join(sigs)} is in the native-depth positivity table")
    smp = _samples(cfg)
    w = .76 / len(smp)
    fig = Figure(figsize=(6.5, 3.2) if t2 is not None else (3.6, 3.2))
    axes = fig.subplots(1, 2 if t2 is not None else 1, sharey=True, squeeze=False)[0]
    x0 = np.arange(len(present))
    labels, pct_native, n_roi = [[], []], {}, {}
    for j, s in enumerate(smp):
        x = x0 + (j - (len(smp) - 1) / 2) * w
        a1 = _by_signature(t1, cfg, s["id"], present)
        k, n, pct = _col(a1, "n_positive"), _col(a1, "n_bins"), _col(a1, "pct_positive")
        with np.errstate(divide="ignore", invalid="ignore"):                  # exact counts beat the rounded column
            pct = np.where(np.isfinite(k) & (n > 0), 100 * k / n, pct)
        wlo, whi = _wilson(k, n)
        lo = np.where(np.isfinite(_col(a1, "ci_lo")), _col(a1, "ci_lo"), wlo)
        hi = np.where(np.isfinite(_col(a1, "ci_hi")), _col(a1, "ci_hi"), whi)
        ok = np.isfinite(pct)
        y = np.maximum(pct, FLOOR)
        e_lo = np.nan_to_num(np.clip(y - np.maximum(lo, FLOOR / 2), 0, None))
        e_hi = np.nan_to_num(np.clip(np.maximum(hi, FLOOR) - y, 0, None))
        if ok.any():
            axes[0].bar(x[ok], y[ok], w, color=s["color"], yerr=[e_lo[ok], e_hi[ok]], capsize=2)
        labels[0] += [(xi, top * 1.25, _count(ki)) for xi, ki, top in zip(x[ok], k[ok], (y + e_hi)[ok])
                      if np.isfinite(ki)]
        pct_native[s["id"]] = dict(zip(present, pct))
        if np.isfinite(n).any():
            n_roi[s["id"]] = float(n[np.isfinite(n)][0])
        if t2 is None:
            continue
        a2 = _by_signature(t2, cfg, s["id"], present)
        mean = _col(a2, "pct_positive_mean")
        sd = np.nan_to_num(_col(a2, "pct_positive_sd"))              # one replicate: SD undefined, no whisker
        kk = _col(a2, "n_positive_mean")
        n2 = np.where(np.isfinite(_col(a2, "n_bins")), _col(a2, "n_bins"), n)
        kk = np.where(np.isfinite(kk), kk, mean / 100 * n2)
        ok2 = np.isfinite(mean)
        y2 = np.maximum(mean, FLOOR)
        top2 = np.maximum(mean + sd, FLOOR)
        e2 = [np.nan_to_num(np.clip(y2 - np.maximum(mean - sd, FLOOR / 2), 0, None)),
              np.nan_to_num(np.clip(top2 - y2, 0, None))]
        if ok2.any():
            axes[1].bar(x[ok2], y2[ok2], w, color=s["color"], yerr=[e2[0][ok2], e2[1][ok2]], capsize=2)
        labels[1] += [(xi, top * 1.25, _mean_count(ki)) for xi, ki, top in zip(x[ok2], kk[ok2], top2[ok2])
                      if np.isfinite(ki)]
    spacing = w * ((fig.get_figwidth() - .8) / len(axes)) / len(present)
    ymax = max(_label_bars(ax, lab, len(smp), spacing, True, 2.2, FLOOR, 150.0) for ax, lab in zip(axes, labels))
    xt = [lbl(cfg, p) for p in present]
    many = len(present) > 3
    for ax, yl in zip(axes, ["% positive units, native depth (log)", "% positive units, depth-matched (log)"]):
        ax.set_yscale("log"); ax.set_xticks(x0)
        ax.set_xticklabels(xt, fontsize=7.5, rotation=30 if many else 0, ha="right" if many else "center")
        ax.set_ylabel(yl, fontsize=7.5); ax.set_ylim(FLOOR, ymax)
    handles = [Patch(color=s["color"], label=s["label"]) for s in smp]
    fig.legend(handles, [s["label"] for s in smp], fontsize=7, frameon=False, loc="lower center",
               ncol=min(len(smp), 4), bbox_to_anchor=(.5, 1.0))
    fig.tight_layout()
    save(fig, out, name, dpi=_dpi(cfg, "dpi", 200))
    sc_ = cfg["scoring"]
    reps = int(sc_["n_depth_reps"])
    caption = (f"Percentage of region-of-interest bins positive for each signature (score above the "
               f"{sc_['null_percentile']:g}th percentile of a random-gene-set null, at least {sc_['min_core']} core "
               f"genes detected, no competing signature scoring higher). "
               f"{'Left, native' if t2 is not None else 'Native'} depth: error bars, Wilson 95% confidence intervals; "
               f"numbers above bars, positive bins. ")
    if t2 is not None:
        spread = f"mean ± SD over {reps} replicates" if reps > 1 else "a single replicate (no SD)"
        caption += (f"Right, after depth matching{_depth_text(cfg, summary)}: {spread}; numbers, "
                    f"mean positive bins. ")
    else:
        caption += "Depth-matched positivity was not computed for this run. "
    caption += ("Logarithmic axes; values below 0.001% (including zero) sit at the axis floor without a bar. "
                f"Region of interest: {_sizes_text(cfg, n_roi)}.")
    caption += _fold_text(cfg, pct_native, present, lambda p: lbl(cfg, p), "Native-depth fold")
    missing = [lbl(cfg, p) for p in sigs if p not in present]
    if missing:
        caption += f" Not in the positivity table (no positivity call): {_join(missing)}."
    return [_entry(name, _title_from_name(name), caption)]


# ============================================================================
# FIGURE 4h (osteoclast compartments)
# ============================================================================

@_house_style
def figure4h(t4, cfg, out, top_n=DE_TOP_N):
    """[RUN] Top eight genes per osteoclast compartment by Wilcoxon score; gray bars did not reach adjusted P < 0.05.
    [FIX] the x axis was fixed at (-13.5, 10), which clipped longer bars and their P-value labels; it now follows the
    data. A table with one compartment only is drawn with that side; rows without a score are ignored.
    [APP] bar colors from the configured compartment signature (tumor-associated, right) and split_by signature
    (growth-plate, left); adjusted P values fall back to raw P values when the table has none."""
    p_col = next((c for c in ("pvals_adj", "pvals") if c in t4.columns), None)
    t = t4.assign(scores=_col(t4, "scores"), _p=_col(t4, p_col) if p_col else np.nan)
    t = t[np.isfinite(t["scores"])]
    ta = t[t["group"] == "tumor_associated"].nlargest(top_n, "scores")
    gp = t[t["group"] == "growth_plate"].nlargest(top_n, "scores")
    rows = pd.concat([ta, gp.iloc[::-1]])
    if rows.empty:
        raise ValueError("no gene of either compartment has a score")
    sig, split = cfg["compartments"]["signature"], cfg["compartments"]["split_by"]
    is_ta = (rows["group"] == "tumor_associated").to_numpy()
    v = np.where(is_ta, rows["scores"], -rows["scores"])
    lo, hi = min(0.0, float(v.min())), max(0.0, float(v.max()))
    span = (hi - lo) or 1.0
    fig = Figure(figsize=(5.2, max(1.6, 1.2 + .2 * len(rows))))
    ax = fig.subplots()
    for k, (vk, right, p) in enumerate(zip(v, is_ta, rows["_p"])):
        col = (_sig_color(cfg, sig) if right else _sig_color(cfg, split)) if p < .05 else NOT_SIGNIFICANT_GRAY
        ax.barh(k, vk, color=col)
        pos = vk > 0 or (vk == 0 and right)
        ax.text(vk + (.01 if pos else -.01) * span, k, f"{p:.1e}" if np.isfinite(p) else "n/a", va="center",
                ha="left" if pos else "right", fontsize=6)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows["names"].astype(str), fontsize=7.5, style="italic")
    ax.axvline(0, color="k", lw=.6)
    ax.set_xlim(lo - (.15 if lo < 0 else .03) * span, hi + (.15 if hi > 0 else .03) * span)
    ax.set_xlabel("Wilcoxon score: growth-plate compartment (left) / tumor-associated (right)", fontsize=7.5)
    save(fig, out, FIG4H, dpi=_dpi(cfg, "dpi", 200))
    S, P, comp = lbl(cfg, sig), lbl(cfg, split), cfg["compartments"]
    tumor = _join(comp.get("tumor_conditions") or ["tumor"])
    where = (f"both compartments from sections whose condition is {tumor}" if comp.get("de_sections") == "tumor"
             else "growth-plate bins from every section")
    caption = (f"Genes higher in each {S} compartment: Wilcoxon rank-sum test between the two compartments on the "
               f"log-normalized expression of {S}-positive bins. Growth-plate compartment (left): {S}-positive bins "
               f"whose {P} score is above the {comp['percentile']:g}th percentile of all bins; tumor-associated (right): "
               f"the other {S}-positive bins of sections whose condition is {tumor} (non-growth-plate bins of other "
               f"sections are not tested); {where}. Top {top_n} genes per compartment by score "
               f"(shown: {len(gp)} growth-plate, {len(ta)} tumor-associated); ")
    caption += {"pvals_adj": "Benjamini-Hochberg adjusted P next to each bar; gray bars: adjusted P ≥ 0.05.",
                "pvals": "unadjusted P next to each bar (the table has no adjusted P); gray bars: P ≥ 0.05.",
                None: "the table has no P values (gray bars, 'n/a')."}[p_col]
    return [_entry(FIG4H, f"Figure 4h — {S} compartment differential expression", caption)]


# ============================================================================
# FIGURE 4i (distance to the nearest partner)
# ============================================================================

@_house_style
def figure4i(dist, summ, cfg, out):
    """[RUN] Cumulative distributions of osteoblast-to-nearest-osteoclast distance, observed (solid) against the
    random-placement null (dashed), per limb.
    [FIX] a section reported 'not estimable' has no observed distances and the x limit became NaN (ValueError); all
    distances 0 gave identical limits. Such a section now gets a panel that states the reason.
    [APP] one panel per configured sample; the two signatures from cfg['uncoupling']; observed/null ratios and the
    ratio of ratios vs the reference sample in the caption (summ: Figure4i_observed_to_null_summary.csv, may be
    empty)."""
    a, b = cfg["uncoupling"]["a"], cfg["uncoupling"]["b"]
    A, B = lbl(cfg, a), lbl(cfg, b)
    smp = _samples(cfg)
    fig = Figure(figsize=(min(3.25 * len(smp), 13.0) if len(smp) > 1 else 3.6, 3.0))
    axes = fig.subplots(1, len(smp), sharey=True, squeeze=False)[0]
    status = dict(zip(summ["sample"], summ["status"].fillna(""))) if {"sample", "status"} <= set(summ.columns) else {}
    d = dist.assign(distance_um=_col(dist, "distance_um"), count=_col(dist, "count"))
    d = d[np.isfinite(d["distance_um"]) & (d["count"] > 0)]
    bin_um = pc.bin_um(cfg["bin_size"])

    def q99(dd):                                   # smallest distance reaching 99% of the units
        cw = dd["count"].cumsum() / dd["count"].sum()
        return float(dd["distance_um"][cw >= .99 - 1e-9].min())

    for ax, s in zip(axes, smp):
        ds = d[d["sample"] == s["id"]]
        obs = ds[ds["source"] == "observed"].sort_values("distance_um")
        null = ds[ds["source"] == "null"].sort_values("distance_um")
        if not len(obs):
            ax.set_xticks([])
            reason = "\n".join(textwrap.wrap(status.get(s["id"]) or "no observed distances", 40))
            ax.text(.5, .5, f"{s['label']}:\n{reason}", transform=ax.transAxes, ha="center", va="center", fontsize=7,
                    color="0.35")
            ax.set_xlabel(f"distance to nearest {B}-positive unit (µm)", fontsize=7.5)
            continue
        upper = 1.1 * q99(obs)
        if upper <= 0:                             # every distance 0: show the null's range
            upper = 1.1 * q99(null) if len(null) else 0.0
        upper = upper if upper > 0 else bin_um
        for dd, src, sty in ((obs, "observed", "-"), (null, "null", "--")):
            if not len(dd):
                continue
            x, y = dd["distance_um"].to_numpy(), (dd["count"].cumsum() / dd["count"].sum()).to_numpy()
            if len(x) == 1:                        # one distance: draw its jump and plateau (a step from one point is empty)
                x, y = np.r_[x, x, max(upper, x[0])], np.r_[0.0, y, y]
            ax.step(x, y, sty, where="post", color=s["color"], label=f"{s['label']}, {src}")
        ax.set_xlim(0, upper)
        ax.legend(fontsize=7, frameon=False, loc="lower right")
        ax.set_xlabel(f"distance to nearest {B}-positive unit (µm)", fontsize=7.5)
    axes[0].set_ylabel(f"cumulative fraction of {A}-positive units", fontsize=7.5)
    fig.tight_layout()
    save(fig, out, FIG4I, dpi=_dpi(cfg, "dpi", 200))
    u = cfg["uncoupling"]
    calls = {"pooled_dm": "spatially pooled, depth-matched", "pooled": "spatially pooled", "native": "native-depth"}
    caption = (f"Distance from each {A}-positive bin to the nearest {B}-positive bin within the region of interest "
               f"({calls.get(u['calls'], u['calls'])} positivity calls), observed (solid) against a random-placement "
               f"null (dashed: {u['n_perm']} placements of the same numbers of positive bins among the bins of the same "
               f"region, pooled). ")
    parts, rr = [], []
    by_sample = {r["sample"]: r for r in summ.to_dict("records")} if "sample" in summ.columns else {}
    for s in smp:
        r = by_sample.get(s["id"])
        if r is None:
            continue
        ratio, rlo, rhi = (float(pd.to_numeric(r.get(c), errors="coerce"))
                           for c in ("ratio_obs_null_mean", "ratio_obs_null_lo", "ratio_obs_null_hi"))
        if str(r.get("status", "")).startswith("estimable") and np.isfinite(ratio):
            text = f"{s['label']} {ratio:.2f}"
            if np.isfinite(rlo) and np.isfinite(rhi):
                text += f" (95% null interval {rlo:.2f}–{rhi:.2f})"
            na, nb = r.get("n_a"), r.get("n_b")
            if pd.notna(na) and pd.notna(nb):
                text += f", {_count(na)} {A}- and {_count(nb)} {B}-positive bins"
            parts.append(text)
            rr_v = float(pd.to_numeric(r.get("ratio_of_ratios_vs_reference"), errors="coerce"))
            if s["id"] != str(cfg.get("reference_sample")) and np.isfinite(rr_v):
                rr.append(f"{s['label']} {rr_v:.2f}")
        else:
            parts.append(f"{s['label']} {r.get('status') or 'not estimable'}")
    if parts:
        caption += f"Observed / null mean distance: {'; '.join(parts)}."
    ref = _reference(cfg)
    if rr and ref:
        caption += f" Ratio of ratios vs {ref['label']}: {', '.join(rr)}."
    caption += " Raw distances depend on the density of positive bins and are not compared between sections."
    return [_entry(FIG4I, f"Figure 4i — {A}-to-{B} distance", caption)]


# ============================================================================
# FIGURE 6b, SUPPLEMENTARY S4 (detection of single genes, depth-matched)
# ============================================================================

@_house_style
def detection_panel(t, genes, cfg, out, name, title, intro, shade_first_two=False, comma=False, summary=None,
                    panel=None):
    """[RUN] Percentage of ROI bins detecting each gene after depth matching, Wilson 95% intervals.
    [FIX] .loc[genes] raised KeyError for a gene absent from the table (genes outside the probe set are left out and
    named in the caption); an undetected gene gave a slightly negative whisker (Wilson bounds of the rounded count
    against the unrounded rate), which matplotlib refuses: every whisker is clipped at 0.
    [APP] any number of samples; the shaded pair is the first two configured genes that are drawn; the count labels
    sit a fixed fraction of the axis above the whiskers (they were 0.02 or 0.04 percentage points above)."""
    tg = t.copy()
    tg["gene"] = tg["gene"].astype(str)
    present = [g for g in dict.fromkeys(genes) if g in set(tg["gene"])]
    if not present:
        raise ValueError(f"none of {', '.join(genes)} is in the table")
    smp = _samples(cfg)
    w = .76 / len(smp)
    fig = Figure(figsize=(6.5, 3.2 if shade_first_two else 3.0))
    ax = fig.subplots()
    x0 = np.arange(len(present))
    drawn, pct_by_sample, n_roi = [], {}, {}
    for j, s in enumerate(smp):
        sub = tg[tg["sample"] == s["id"]].drop_duplicates("gene").set_index("gene").reindex(present)
        p, k, n = _col(sub, "pct_detected_dm"), _col(sub, "n_positive_dm"), _col(sub, "n_bins")
        wlo, whi = _wilson(k, n)
        lo = np.where(np.isfinite(_col(sub, "ci_lo_dm")), _col(sub, "ci_lo_dm"), wlo)
        hi = np.where(np.isfinite(_col(sub, "ci_hi_dm")), _col(sub, "ci_hi_dm"), whi)
        e_lo, e_hi = np.nan_to_num(np.clip(p - lo, 0, None)), np.nan_to_num(np.clip(hi - p, 0, None))
        x = x0 + (j - (len(smp) - 1) / 2) * w
        ok = np.isfinite(p)
        if ok.any():
            ax.bar(x[ok], p[ok], w, yerr=[e_lo[ok], e_hi[ok]], color=s["color"], capsize=2)
        drawn.append((x[ok], k[ok], (p + e_hi)[ok]))
        pct_by_sample[s["id"]] = dict(zip(present, p))
        if np.isfinite(n).any():
            n_roi[s["id"]] = float(n[np.isfinite(n)][0])
    top = max([float(np.max(tops)) for _, _, tops in drawn if len(tops)] + [0.0])
    offset = .015 * (top or 1.0)
    labels = [(xi, ti + offset, f"{int(round(ki)):,}" if comma else str(int(round(ki))))
              for xs_, ks, tops in drawn for xi, ki, ti in zip(xs_, ks, tops) if np.isfinite(ki)]
    ymax = _label_bars(ax, labels, len(smp), w * 5.7 / len(present), False, 2.3, 0.0,
                       ax.get_ylim()[1] if top else 1.0)
    shaded = [g for g in genes[:2] if g in present] if shade_first_two else []
    if shaded:
        ax.axvspan(-.5, len(shaded) - .5, color="0.93", zorder=-1)            # RANKL/OPG pair
    _tick_labels(ax, present, 8, 5.5 / len(present))                           # slanted when names do not fit
    for tick in ax.get_xticklabels():
        tick.set_style("italic")
    ax.set_xlim(-.5, len(present) - .5)
    ax.set_ylabel("% units detecting (depth-matched)", fontsize=8)
    ax.set_ylim(0, ymax)
    ax.legend([Patch(color=s["color"]) for s in smp], [s["label"] for s in smp], fontsize=7, frameon=False)
    save(fig, out, name, dpi=_dpi(cfg, "dpi", 200))
    reps = int(cfg["scoring"]["n_depth_reps"])
    caption = (f"{intro} Percentage of region-of-interest bins detecting each gene after depth matching"
               f"{_depth_text(cfg, summary)} (mean over {reps} "
               f"replicate{'s' if reps != 1 else ''}); error bars, Wilson 95% confidence intervals; numbers above "
               f"bars, detecting bins. Region of interest: {_sizes_text(cfg, n_roi)}.")
    if shaded:
        roles = dict(zip(tg["gene"], tg["role"].astype(str))) if "role" in tg.columns else {}
        caption += " Shaded: " + _join(f"{g} ({roles[g]})" if roles.get(g) else g for g in shaded) + "."
    caption += _fold_text(cfg, pct_by_sample, present, str)
    missing = [g for g in genes if g not in present]
    unanalysed = [g for g in missing if panel is not None and g not in panel]
    absent = [g for g in missing if g not in unanalysed]
    if absent:
        caption += f" Not in the probe set (not drawn): {_join(absent)}."
    if unanalysed:
        caption += f" Not in the analysed gene panel (not drawn): {_join(unanalysed)}."
    return [_entry(name, title, caption)]


# ============================================================================
# SUPPLEMENTARY FIGURES S1, S2, S5
# ============================================================================

def _retention(qc, thresholds):
    """[APP] Retention sweep (sample, min_genes, pct_retained) computed from the per-bin QC table, for runs without
    qc_threshold_sweep_<bin>.csv (bins kept: more than min_genes genes, as in bone_analysis.qc_sweep)."""
    rows = []
    genes = _col(qc, "n_genes_by_counts")
    for sid in pd.unique(qc["sample"]):
        g = genes[(qc["sample"] == sid).to_numpy()]
        for t in thresholds:
            rows.append({"sample": sid, "min_genes": t, "pct_retained": 100 * float(np.mean(g > t)) if len(g) else np.nan})
    return pd.DataFrame(rows)


@_house_style
def figureS1(qc, sweep, cfg, out, bin_size=None):
    """[RUN] Per-bin QC distributions at 8 um and retention against the minimum gene count (solid line: the 10-gene
    criterion used; dashed: 100 genes).
    [FIX] it re-read the raw Space Ranger 8 um matrices (bone_analysis.read_visium_hd) and, after writing the PDF,
    imported pymupdf (not a declared dependency) to render the PNG; the ImportError also stopped Supplementary Figures
    S2, S4 and S5. It now draws from the per-bin QC table of bone_analysis.py (FigureS1_qc_per_bin_<bin>.csv.gz) and
    matplotlib writes both files.
    [APP] bin size of the table (default cfg['qc_bin_size']; an optional bin_size column gives it per sample, for
    sections without bins of the QC size), samples and the criterion line from the configuration; bins other than
    the QC size are named in the caption; sweep None -> the retention curve is computed from the per-bin table."""
    qc_size = cfg["qc_bin_size"]
    sizes = {k: v for k, v in zip(qc["sample"], qc["bin_size"]) if v} if "bin_size" in qc.columns else {}
    mixed = len(set(sizes.values())) > 1
    if len(set(sizes.values())) == 1:
        bin_size = next(iter(sizes.values()))
    bin_size = bin_size or qc_size
    unit = " per bin" if mixed else f" per {_bin_text(bin_size)} bin"
    if sweep is None or sweep.empty:
        sweep = _retention(qc, sorted({int(t) for t in cfg["qc"]["sweep_thresholds"]}))
    smp = _samples(cfg)
    fig = Figure(figsize=(13, 3.2))
    ax = fig.subplots(1, 4)
    n_bins = {}
    genes, umis, mt = _col(qc, "n_genes_by_counts"), _col(qc, "total_counts"), _col(qc, "pct_counts_mt")
    for s in smp:
        m = (qc["sample"] == s["id"]).to_numpy()
        n_bins[s["id"]] = int(m.sum())
        c = s["color"]
        with np.errstate(invalid="ignore"):
            values = (np.log10(genes[m] + 1), np.log10(umis[m] + 1), mt[m])
        for axis, v in zip(ax[:3], values):
            v = v[np.isfinite(v)]
            if len(v):
                axis.hist(v, bins=60, alpha=.55, label=s["label"], color=c, density=True)
        sw = sweep[sweep["sample"] == s["id"]]
        sw = sw.assign(min_genes=_col(sw, "min_genes"), pct_retained=_col(sw, "pct_retained"))
        sw = sw[(sw["min_genes"] > 0) & np.isfinite(sw["pct_retained"])].sort_values("min_genes")
        if len(sw):
            ax[3].plot(sw["min_genes"], sw["pct_retained"], "o-", color=c, label=s["label"])
    ax[0].set_xlabel("log10 genes" + unit); ax[1].set_xlabel("log10 UMIs" + unit)
    ax[2].set_xlabel("% mitochondrial counts")
    for x in ax[:3]:
        x.set_ylabel("density")
    min_genes = float(cfg["qc"]["min_genes"])
    ax[3].axvline(SINGLE_CELL_MIN_GENES, ls="--", c="k", lw=.8)
    if min_genes > 0:
        ax[3].axvline(min_genes, ls="-", c="0.4", lw=.8)
    ax[3].set_xscale("log"); ax[3].set_xlabel("minimum genes" + unit + " (log)"); ax[3].set_ylabel("% bins retained")
    for x in ax:
        if x.get_legend_handles_labels()[0]:
            x.legend(frameon=False, fontsize=7)
        else:                                      # column absent from the per-bin QC table
            x.text(.5, .5, "not available", transform=x.transAxes, ha="center", va="center", fontsize=7,
                   color="0.35")
    fig.tight_layout()
    save(fig, out, FIGS1, dpi=_dpi(cfg, "dpi", 200))
    q = cfg["qc"]
    counted = _sizes_text(cfg, {k: v for k, v in n_bins.items() if v})
    if mixed:
        caption = ("Quality control of the bins before filtering (" + "; ".join(
            f"{s['label']}: {_bin_text(sizes[s['id']])}" for s in smp if s["id"] in sizes) + f"; {counted}).")
    else:
        caption = f"Quality control of the {_bin_text(bin_size)} bins before filtering ({counted})."
    caption += (f" Left to right: genes per bin and UMIs per bin (log10(x + 1)), % mitochondrial counts (densities), "
                f"and % of bins retained against the minimum gene count (log axis; dashed line: "
                f"{SINGLE_CELL_MIN_GENES} genes, a threshold carried over from single-cell work")
    caption += f"; solid line: the {q['min_genes']:g}-gene criterion used)." if min_genes > 0 else ")."
    caption += (f" Bins kept for analysis: more than {q['min_genes']:g} genes, more than {q['min_counts']:g} UMIs and "
                f"less than {q['max_pct_mt']:g}% mitochondrial counts.")
    if mixed:
        other = [s["label"] for s in smp if sizes.get(s["id"], qc_size) != qc_size]
        caption += f" No {_bin_text(qc_size)} bins for {_join(other)}: the analysis bins are shown."
    elif bin_size != qc_size:
        caption += (f" The {_bin_text(qc_size)} bins (qc_bin_size) were not available; the figure shows the "
                    f"{_bin_text(bin_size)} analysis bins instead.")
    return [_entry(FIGS1, "Supplementary Figure S1 — quality control", caption)]


@_house_style
def figureS2(cons, cfg, out, top_n=None):
    """[RUN] Moran's I per section for the top 30 (by mean rank) of the genes in the top 200 of both sections.
    [FIX] the columns moranI_LL-2-1 and moranI_RL-1-5 and the two sections' names were hard-coded.
    [APP] one bar per configured sample with a moranI_<sample id> column; the number of genes from
    cfg['figures']['moran_top_n']; the height grows with it (6.2 in for 30 genes)."""
    top_n = int(top_n or cfg["figures"]["moran_top_n"])
    smp = [s for s in _samples(cfg) if f"moranI_{s['id']}" in cons.columns]
    if not smp:
        raise ValueError("no moranI_<sample id> column of a configured sample")
    c6 = cons.assign(mean_rank=_col(cons, "mean_rank")).sort_values("mean_rank", kind="stable").head(top_n)
    n = len(c6)
    hb = .8 / len(smp)
    fig = Figure(figsize=(4.6, max(1.6, .5 + .19 * n)))
    ax = fig.subplots()
    y = np.arange(n)[::-1]
    for j, s in enumerate(smp):
        v = _col(c6, f"moranI_{s['id']}")
        ok = np.isfinite(v)
        ax.barh(y[ok] + ((len(smp) - 1) / 2 - j) * hb, v[ok], hb, color=s["color"], label=s["label"])
    ax.set_yticks(y); ax.set_yticklabels(c6["gene"].astype(str), fontsize=6.5, style="italic")
    ax.set_ylim(-.6, n - .4)
    ax.set_xlabel("Moran's I (per section)", fontsize=8); ax.legend(fontsize=7, frameon=False)
    save(fig, out, FIGS2, dpi=_dpi(cfg, "dpi", 200))
    m = cfg["moran"]
    every = "every section" if len(cfg["samples"]) > 1 else "the section"
    caption = (f"Moran's I of each section (computed separately per section, Delaunay graph of bin centres) for the "
               f"top {n} genes by mean rank among the {len(cons)} genes in the top {m['n_consensus']} of {every} "
               f"(of the {m['n_genes']} highly variable genes tested).")
    missing = [s["label"] for s in _samples(cfg) if f"moranI_{s['id']}" not in cons.columns]
    if missing:
        caption += f" No Moran's I for {_join(missing)}."
    return [_entry(FIGS2, "Supplementary Figure S2 — spatially variable genes (Moran's I)", caption)]


@_house_style
def figureS5(M, cfg, out, unresolved=None):
    """[RUN] Mean signature score per Leiden domain, z-scored across domains (color limits -2 to 2).
    [APP] cells without a value (a signature constant across domains) are light gray; the height grows with the
    number of domains; the unresolved share (FigureS5_unresolved_bins.csv) goes to the caption."""
    M = M.set_index("leiden_domain")
    vals = M.apply(pd.to_numeric, errors="coerce").to_numpy(float)
    fig = Figure(figsize=(6.5, max(5.2, 1.6 + .26 * len(M))))
    ax = fig.subplots()
    im_ = ax.imshow(vals, cmap=matplotlib.colormaps["RdBu_r"].with_extremes(bad="0.85"), vmin=-2, vmax=2,
                    aspect="auto")
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index, fontsize=7)
    ax.set_xlabel("signature score", fontsize=8)
    ax.set_ylabel("Leiden domain (named by most enriched signature)", fontsize=8)
    cb = fig.colorbar(im_, ax=ax, shrink=.7)
    cb.set_label("z-scored mean across domains", fontsize=7)
    fig.tight_layout()
    save(fig, out, FIGS5, dpi=_dpi(cfg, "dpi", 200))
    caption = (f"Mean signature score of the bins of each Leiden domain ({len(M)} domains), z-scored across domains "
               f"(color limits −2 to 2; gray: no variation across domains). Domains are named by the signature most "
               f"enriched among their assigned bins; 'unresolved' merges domains whose top markers are dominated by "
               f"ubiquitous transcripts, 'unassigned' collects domains of fewer than "
               f"{cfg['clustering']['min_domain_bins']} bins.")
    if unresolved is not None and len(unresolved):
        r = unresolved.iloc[0]
        un, tot, pct = (float(pd.to_numeric(r.get(c), errors="coerce")) for c in ("unresolved_bins", "all_bins",
                                                                                   "percent"))
        if np.isfinite(un) and np.isfinite(tot):
            caption += f" Unresolved: {_count(un)} of {_count(tot)} bins"
            caption += f" ({pct:.1f}%)." if np.isfinite(pct) else "."
    return [_entry(FIGS5, "Supplementary Figure S5 — Leiden domain heatmap", caption)]


# ============================================================================
# DRIVER  [APP]
# ============================================================================

class _Skip(Exception):
    """[APP] A figure's input is absent or has nothing to draw; the message is the reason given to progress.skip."""


def _absent(cfg, name, step):
    if step and not cfg["steps"].get(step, True):
        return f"step disabled: {pc.STEP_LABELS.get(step, step)}"
    return f"input missing: {name}"


class _Context:
    """[APP] What the figure steps share: results and figures folders, configuration, summary.json, root, progress,
    and the tables read so far."""

    def __init__(self, results, cfg, progress, root):
        self.results, self.cfg, self.progress, self.root = Path(results), cfg, progress, root
        self.out = self.results / "figures"
        self.summary = _read_json(self.results / "summary.json")
        self._cache = {}

    def load(self, name, columns=(), step=None, nrows=None, usecols=None):
        """results/<name> as a table; raises _Skip when the file is absent or empty, lacks one of `columns` or has
        no rows. Full reads are cached (figure functions do not modify their inputs)."""
        path = self.results / name
        if not path.is_file():
            raise _Skip(_absent(self.cfg, name, step))
        if nrows is not None:
            return self._check(_read_table(path, usecols=usecols, nrows=nrows), name, columns, step)
        key = (name, tuple(usecols) if usecols else None)
        if key not in self._cache:
            self._cache[key] = _read_table(path, usecols=usecols)
        return self._check(self._cache[key], name, columns, step)

    def _check(self, df, name, columns, step):
        if not len(df.columns):
            raise _Skip(f"{name} is empty")
        missing = [c for c in columns if c not in df.columns]
        if missing:
            disabled = step and not self.cfg["steps"].get(step, True)
            raise _Skip(f"{name} has no column {', '.join(missing)}"
                        + (f" (step disabled: {pc.STEP_LABELS.get(step, step)})" if disabled else ""))
        if df.empty:
            raise _Skip(f"{name} has no rows")
        return df

    def optional(self, name, columns=()):
        """A table that only adds to a figure: None when absent, empty or incomplete."""
        try:
            return self.load(name, columns)
        except _Skip:
            return None

    def per_bin(self, columns, peek=False):
        """The per-bin table restricted to the columns the maps use (read once; peek=True reads one row, enough to
        check `columns` and that there are rows)."""
        cfg = self.cfg
        lm = cfg["figures"]["lineage_maps"]
        cols = ["sample_id", "x_um", "y_um", "in_region_of_interest", "tissue_class",
                f"positive_{_slug(cfg, lm['a'])}", f"positive_{_slug(cfg, lm['b'])}"]
        return self.load(_per_bin_csv(cfg), columns, nrows=1 if peek else None, usecols=cols)


class _Step:
    """[APP] One progress step of make_figures: its label, the figures it draws as (name, title), and prepare(ctx),
    which reads and checks the inputs (raising _Skip) and returns the drawing tasks [((name, title) list, draw)]."""

    def __init__(self, label, figures, prepare):
        self.label, self.figures, self.prepare = label, figures, prepare


def _prep_1a(ctx):
    return [([(FIG1A, "Figure 1a — experimental design")],
             lambda: figure1a(ctx.out, dpi=_dpi(ctx.cfg, "map_dpi", 250)))]


def _prep_2c(ctx):
    ids = {s["id"] for s in ctx.cfg["samples"]}
    reason = None
    if ctx.root is None:
        reason = "the H&E crops need the Space Ranger folder (root)"
    elif not all(sid in ids for sid, *_ in FIG2C):
        reason = f"the H&E crops are defined for the manuscript sections {_join(sid for sid, *_ in FIG2C)} only"
    if reason:
        ctx.progress.warn(f"Figure 2c not drawn: {reason}")
        raise _Skip(reason)
    return [([(name, title) for _, _, name, title in FIG2C], lambda: figure2c(ctx.root, ctx.cfg, ctx.out))]


def _prep_3a(ctx):
    cfg = ctx.cfg
    if not cfg["steps"].get("assignment", True):
        raise _Skip(_absent(cfg, "", "assignment"))
    cols = ["sample_id", "x_um", "y_um", "tissue_class"]
    ctx.per_bin(cols, peek=True)
    tasks = []
    for i, s in enumerate(_samples(cfg), start=1):
        def draw(sid=s["id"], index=i):
            return figure3a(ctx.per_bin(cols), cfg, ctx.out, sid, index)
        tasks.append(([(_figure3a_name(i, s["label"]), f"Figure 3a — tissue classes, {s['label']}")], draw))
    return tasks


def _prep_3b(ctx):
    comp = ctx.load(COMPOSITION_CSV, ["sample", "program", "n_bins"], step="assignment")
    return [([(FIG3B, "Figure 3b — tissue composition after depth matching")],
             lambda: figure3b(comp, ctx.cfg, ctx.out, ctx.summary))]


def _prep_4af(ctx):
    cfg = ctx.cfg
    lm = cfg["figures"]["lineage_maps"]
    for sig in (lm["a"], lm["b"]):
        if sig not in cfg["primary"]:
            raise _Skip(f"{lbl(cfg, sig)} has no positivity call (not a primary signature)")
    cols = ["sample_id", "x_um", "y_um", f"positive_{_slug(cfg, lm['a'])}", f"positive_{_slug(cfg, lm['b'])}"]
    ctx.per_bin(cols, peek=True)
    return [([(FIG4AF, f"Figures 4a-f — {lbl(cfg, lm['a'])}- and {lbl(cfg, lm['b'])}-positive bins")],
             lambda: figure4af(ctx.per_bin(cols), cfg, ctx.out))]


def _prep_positivity(ctx, name, sigs):
    cfg = ctx.cfg
    native = ctx.load(NATIVE_CSV, ["signature", "sample"])
    if "pct_positive" not in native.columns and not {"n_positive", "n_bins"} <= set(native.columns):
        raise _Skip(f"{NATIVE_CSV} has neither pct_positive nor n_positive and n_bins")
    names = set(native["signature"].astype(str).map(_signature_names(cfg)).dropna())
    if not any(s in names for s in sigs):
        raise _Skip(f"no positivity calls for {_join(lbl(cfg, s) for s in sigs) or 'any signature'} in {NATIVE_CSV}")
    dm = ctx.optional(DEPTH_MATCHED_CSV, ["signature", "sample", "pct_positive_mean"])
    return [([(name, _title_from_name(name))],
             lambda: positivity_panel(native, dm, cfg, ctx.out, name, sigs, ctx.summary))]


def _prep_4h(ctx):
    de = ctx.load(DE_CSV, ["group", "names", "scores"], step="compartments")
    if not de["group"].isin(["growth_plate", "tumor_associated"]).any():
        raise _Skip(f"{DE_CSV} has no gene of the growth_plate or tumor_associated compartment")
    sig = ctx.cfg["compartments"]["signature"]
    return [([(FIG4H, f"Figure 4h — {lbl(ctx.cfg, sig)} compartment differential expression")],
             lambda: figure4h(de, ctx.cfg, ctx.out))]


def _prep_4i(ctx):
    cfg = ctx.cfg
    summ = ctx.optional(DISTANCE_SUMMARY_CSV, ["sample", "status"])
    if summ is None:
        summ = pd.DataFrame(columns=["sample", "status"])
    dist = ctx.optional(DISTANCE_CSV, ["sample", "source", "distance_um", "count"])
    if dist is None and summ.empty:
        if not (ctx.results / DISTANCE_CSV).is_file():
            raise _Skip(_absent(cfg, DISTANCE_CSV, "uncoupling"))
        raise _Skip(f"{DISTANCE_CSV} has no distances")
    ids = {s["id"] for s in _samples(cfg)}
    if dist is None or not ((dist["source"] == "observed") & dist["sample"].isin(ids)).any():
        statuses = [f"{s['label']}: {st}" for s in _samples(cfg)
                    for st in summ.loc[summ["sample"] == s["id"], "status"].head(1)]
        raise _Skip("no section has observed distances" + (f" ({'; '.join(statuses)})" if statuses else ""))
    A, B = lbl(cfg, cfg["uncoupling"]["a"]), lbl(cfg, cfg["uncoupling"]["b"])
    return [([(FIG4I, f"Figure 4i — {A}-to-{B} distance")], lambda: figure4i(dist, summ, cfg, ctx.out))]


def _prep_detection(ctx, csv, step, genes, name, title, intro, shade, comma, panel):
    t = ctx.load(csv, ["gene", "sample", "pct_detected_dm"], step=step)
    if not any(g in set(t["gene"].astype(str)) for g in genes):
        raise _Skip(f"none of {_join(genes) or 'the configured genes'} is in {csv} (outside the probe set or the "
                    f"analysed panel)")
    return [([(name, title)], lambda: detection_panel(t, genes, ctx.cfg, ctx.out, name, title, intro,
                                                      shade_first_two=shade, comma=comma, summary=ctx.summary,
                                                      panel=set(panel)))]


def _prep_6b(ctx):
    cfg = ctx.cfg
    return _prep_detection(ctx, RECEPTOR_CSV, "receptor_panel", list(cfg["figures"]["receptor_genes"]), FIG6B,
                           "Figure 6b — host receptor panel", "Host receptor panel.",
                           bool(cfg["figures"].get("receptor_shade_first_two")), False, cfg["receptor_panel"])


def _prep_s4(ctx):
    return _prep_detection(ctx, ADIPOCYTE_CSV, "adipocyte_panel", list(ctx.cfg["adipocyte_genes"]), FIGS4,
                           "Supplementary Figure S4 — adipocyte genes", "Adipocyte genes.", False, True,
                           ctx.cfg["adipocyte_genes"])


def _prep_s1(ctx):
    cfg = ctx.cfg
    q_size, b_size = cfg["qc_bin_size"], cfg["bin_size"]
    cols = ["sample", "n_genes_by_counts", "total_counts", "pct_counts_mt", "bin_size"]
    size = q_size
    name = f"FigureS1_qc_per_bin_{_bin_tag(q_size)}.csv.gz"
    fallback = f"FigureS1_qc_per_bin_{_bin_tag(b_size)}.csv.gz"
    if not (ctx.results / name).is_file() and (ctx.results / fallback).is_file():
        name, size = fallback, b_size                                   # analysis-bin fallback
    ctx.load(name, cols[:2], step="qc_figure_bins", nrows=1)
    sweep = ctx.optional(f"qc_threshold_sweep_{_bin_tag(size)}.csv", ["sample", "min_genes", "pct_retained"])

    def draw():
        qc = _read_table(ctx.results / name, usecols=cols)
        return figureS1(qc, sweep, cfg, ctx.out, bin_size=size)
    return [([(FIGS1, "Supplementary Figure S1 — quality control")], draw)]


def _prep_s2(ctx):
    cons = ctx.load(MORAN_CSV, ["gene"], step="moran")
    if not any(f"moranI_{s['id']}" in cons.columns for s in _samples(ctx.cfg)):
        raise _Skip(f"{MORAN_CSV} has no moranI_<sample id> column of a configured sample")
    return [([(FIGS2, "Supplementary Figure S2 — spatially variable genes (Moran's I)")],
             lambda: figureS2(cons, ctx.cfg, ctx.out))]


def _prep_s5(ctx):
    M = ctx.load(DOMAIN_CSV, ["leiden_domain"], step="leiden")
    if M.shape[1] < 2:
        raise _Skip(f"{DOMAIN_CSV} has no signature column")
    unresolved = ctx.optional(UNRESOLVED_CSV, ["unresolved_bins", "all_bins"])
    return [([(FIGS5, "Supplementary Figure S5 — Leiden domain heatmap")],
             lambda: figureS5(M, ctx.cfg, ctx.out, unresolved))]


def _positivity_panels(cfg):
    """[APP] (name, signatures) of the configured positivity panels, and the reasons for those left out: entries
    that are not {name, signatures}, and names that repeat or that another figure uses (they would overwrite its
    files)."""
    reserved = {FIG3B, FIG4AF, FIG4H, FIG4I, FIG6B, FIGS1, FIGS2, FIGS4, FIGS5, FIG1A} | {n for _, _, n, _ in FIG2C}
    panels, ignored = [], []
    for p in cfg["figures"].get("positivity_panels") or []:
        sigs = p.get("signatures") if isinstance(p, dict) else None
        if isinstance(sigs, str):
            sigs = [sigs]
        if not isinstance(sigs, list):
            ignored.append(f"positivity panel {p!r} ignored: expected {{name, signatures: [...]}}")
            continue
        name = _safe_name(p.get("name"))
        if name in reserved or name.startswith("Figure3a_") or name in {n for n, _ in panels}:
            ignored.append(f"positivity panel {name!r} ignored: another figure has that name")
            continue
        panels.append((name, [str(s) for s in sigs]))
    return panels, ignored


def _plan(cfg, manuscript_extras=False):
    """[APP] The figure steps, in order, for a configuration (merged with the defaults); no data access."""
    smp = _samples(cfg)
    steps = []
    if manuscript_extras:
        steps.append(_Step("Draw Figure 1a — experimental design", [(FIG1A, "Figure 1a — experimental design")],
                           _prep_1a))
        steps.append(_Step("Draw Figure 2c — H&E crops", [(n, t) for _, _, n, t in FIG2C], _prep_2c))
    steps.append(_Step("Draw Figure 3a — tissue classes per section",
                       [(_figure3a_name(i, s["label"]), f"Figure 3a — tissue classes, {s['label']}")
                        for i, s in enumerate(smp, start=1)], _prep_3a))
    steps.append(_Step("Draw Figure 3b — tissue composition",
                       [(FIG3B, "Figure 3b — tissue composition after depth matching")], _prep_3b))
    lm = cfg["figures"]["lineage_maps"]
    A, B = lbl(cfg, lm["a"]), lbl(cfg, lm["b"])
    steps.append(_Step(f"Draw Figures 4a-f — {A} and {B} maps", [(FIG4AF, f"Figures 4a-f — {A}- and {B}-positive bins")],
                       _prep_4af))
    for name, sigs in _positivity_panels(cfg)[0]:
        title = _title_from_name(name)
        steps.append(_Step(f"Draw {title}", [(name, title)], functools.partial(_prep_positivity, name=name, sigs=sigs)))
    S = lbl(cfg, cfg["compartments"]["signature"])
    steps.append(_Step(f"Draw Figure 4h — {S} compartment DE",
                       [(FIG4H, f"Figure 4h — {S} compartment differential expression")], _prep_4h))
    UA, UB = lbl(cfg, cfg["uncoupling"]["a"]), lbl(cfg, cfg["uncoupling"]["b"])
    steps.append(_Step(f"Draw Figure 4i — {UA}-to-{UB} distance", [(FIG4I, f"Figure 4i — {UA}-to-{UB} distance")],
                       _prep_4i))
    steps.append(_Step("Draw Figure 6b — host receptor panel", [(FIG6B, "Figure 6b — host receptor panel")], _prep_6b))
    steps.append(_Step("Draw Supplementary Figure S1 — quality control",
                       [(FIGS1, "Supplementary Figure S1 — quality control")], _prep_s1))
    steps.append(_Step("Draw Supplementary Figure S2 — Moran's I",
                       [(FIGS2, "Supplementary Figure S2 — spatially variable genes (Moran's I)")], _prep_s2))
    steps.append(_Step("Draw Supplementary Figure S4 — adipocyte genes",
                       [(FIGS4, "Supplementary Figure S4 — adipocyte genes")], _prep_s4))
    steps.append(_Step("Draw Supplementary Figure S5 — Leiden domain heatmap",
                       [(FIGS5, "Supplementary Figure S5 — Leiden domain heatmap")], _prep_s5))
    seen = set()
    for s in steps:                                    # labels are the progress keys: keep them unique
        if s.label in seen:
            s.label = f"{s.label} ({s.figures[0][0]})"
        seen.add(s.label)
    return steps


def planned_steps(cfg: dict, manuscript_extras: bool = False) -> list[str]:
    """[APP] Labels of every step make_figures reports through progress.step() or progress.skip(), in order, for
    this configuration (partial configurations are merged with the defaults)."""
    return [s.label for s in _plan(pc.merge_with_defaults(cfg), manuscript_extras)]


def _failed(progress, figures, exc):
    """[APP] A figure that raised: a warning (traceback in the log) and manifest entries carrying the error."""
    message = f"{type(exc).__name__}: {exc}" if str(exc) else type(exc).__name__
    progress.warn(f"{_join(name for name, _ in figures)} could not be drawn: {message}")
    progress.info("".join(traceback.format_exception(type(exc), exc, exc.__traceback__)).rstrip())
    return [_error_entry(name, title, message) for name, title in figures]


def _write_manifest(out, entries):
    tmp = out / "manifest.json.tmp"
    tmp.write_text(json.dumps({"figures": entries}, indent=2, ensure_ascii=False), encoding="utf-8")
    os.replace(tmp, out / "manifest.json")


def make_figures(results_dir, cfg: dict | None = None, progress: pc.Progress | None = None, root=None,
                 manuscript_extras: bool = False) -> list[dict]:
    """[APP] Draws every figure whose inputs exist into results_dir/figures/ (PNG and PDF) and writes
    figures/manifest.json; returns its "figures" list.
    cfg None -> results_dir/config.json if present, else the defaults; a given cfg is merged with the defaults.
    Each step reports progress.step(label), or progress.skip(label, reason) when an input is absent or empty (its
    analysis step was disabled); the labels are those of planned_steps(cfg, manuscript_extras). A figure that fails
    becomes progress.warn and a manifest entry with "error"; the other figures still run.
    manuscript_extras=True also draws Figure 1a and the Figure 2c H&E crops (needs `root`, the folder of the Space
    Ranger outputs, and the manuscript sample ids); the web application never does."""
    results = Path(results_dir)
    if not results.is_dir():
        raise FileNotFoundError(f"results folder not found: {results}")
    if cfg is None:
        cfg = pc.load_config(results / "config.json") if (results / "config.json").is_file() else pc.default_config()
    else:
        cfg = pc.merge_with_defaults(cfg)
    progress = progress or pc.Progress()
    ctx = _Context(results, cfg, progress, root)
    ctx.out.mkdir(parents=True, exist_ok=True)
    for message in _positivity_panels(cfg)[1]:
        progress.warn(message)
    entries = []
    for step in _plan(cfg, manuscript_extras):
        try:
            tasks = step.prepare(ctx)
        except _Skip as exc:
            progress.skip(step.label, str(exc))
            continue
        except Exception as exc:                       # unreadable input: the step fails, the others run
            progress.step(step.label)
            entries += _failed(progress, step.figures, exc)
            continue
        progress.step(step.label)
        for figures, draw in tasks:
            try:
                entries += draw()
            except Exception as exc:
                entries += _failed(progress, figures, exc)
    _write_manifest(ctx.out, entries)
    return entries


def main(argv=None) -> int:
    """[APP] Command line; exit code 1 when a figure could not be drawn."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="results", help="output folder of bone_analysis.py (its --out); the "
                                                         "figures go to <results>/figures")
    ap.add_argument("--root", default=".", help="folder with the Space Ranger outputs (Figure 2c crops only)")
    ap.add_argument("--config", help="configuration JSON (default: <results>/config.json, else the manuscript "
                                     "defaults)")
    ap.add_argument("--manuscript-extras", action="store_true", help="also draw Figure 1a and the Figure 2c crops")
    args = ap.parse_args(argv)
    results = Path(args.results)
    if not results.is_dir():
        ap.error(f"results folder not found: {results}")
    if args.config:
        cfg = pc.load_config(args.config)
    elif (results / "config.json").is_file():
        cfg = pc.load_config(results / "config.json")
    else:
        cfg = pc.default_config()
    progress = pc.Progress()
    progress.plan(planned_steps(cfg, manuscript_extras=args.manuscript_extras))
    entries = make_figures(results, cfg, progress, root=args.root, manuscript_extras=args.manuscript_extras)
    progress.done()
    failed = [e["name"] for e in entries if e["error"]]
    print(f"{len(entries) - len(failed)} figures written, {len(failed)} failed -> {results / 'figures'}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
