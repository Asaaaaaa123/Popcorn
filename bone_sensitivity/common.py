"""
Helpers shared by the sensitivity analyses. This module must not import scanpy or bone_analysis at import time: the
figure modules use it, and bone_figures runs without scanpy. Functions that need bone_analysis import it inside.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import pipeline_config as pc

GROUP = "Sensitivity analyses"          # figures/manifest.json group of every sensitivity figure
Z95 = 1.959963984540054


@dataclass
class Output:
    """One table a sensitivity module writes: results/<file>, its title and description (manifest.json)."""
    file: str
    title: str
    description: str
    df: pd.DataFrame


def warn(progress, message):
    """Report a warning through progress (pipeline_config.Progress) or print it."""
    if progress is not None:
        progress.warn(message)
    else:
        print(f"WARNING: {message}")


def settings(cfg, section):
    """cfg["sensitivity_analyses"][section] merged over its defaults (configurations written before 2026-11-02 lack
    the section)."""
    full = pc.merge_with_defaults({"sensitivity_analyses": (cfg or {}).get("sensitivity_analyses") or {}})
    return full["sensitivity_analyses"][section] if section else full["sensitivity_analyses"]


def as_bool(values) -> np.ndarray:
    """A boolean, 0/1 or True/False-string column as a numpy bool array (missing = False)."""
    s = pd.Series(values)
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).to_numpy(bool)
    if pd.api.types.is_numeric_dtype(s):
        return s.fillna(0).to_numpy() != 0
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "1.0", "yes", "t"]).to_numpy()


def sample_ids(a, cfg=None):
    """Section ids in configuration order (or order of appearance)."""
    present = list(pd.unique(a.obs["sample"].astype(str)))
    if cfg:
        return [s for s in pc.sample_ids(cfg) if s in present]
    return present


def tumor_and_reference(cfg):
    """(tumor section ids, reference id): tumor sections are those whose condition is in
    compartments.tumor_conditions."""
    tc = set(cfg["compartments"]["tumor_conditions"])
    return [s["id"] for s in cfg["samples"] if s.get("condition") in tc], cfg["reference_sample"]


def pixel_xy_um(a):
    """Bin centres in um: full-resolution pixels times each section's microns per pixel (proximity tables)."""
    xy = np.asarray(a.obsm["spatial"], dtype=float).copy()
    mpp = a.uns.get("microns_per_pixel", {})
    smp = a.obs["sample"].astype(str).to_numpy()
    for s in pd.unique(smp):
        xy[smp == s] *= float(mpp[s])
    return xy


def array_xy_um(a, bin_um):
    """Bin centres in um from the array grid (Figure 4i): array_col, array_row times the bin size."""
    return np.column_stack([a.obs["array_col"].to_numpy(float), a.obs["array_row"].to_numpy(float)]) * float(bin_um)


def wilson(k, n, z=Z95):
    """Wilson 95% interval in percent (NaN where n is 0)."""
    k, n = np.asarray(k, float), np.asarray(n, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = k / n
        den = 1 + z ** 2 / n
        centre = (p + z ** 2 / (2 * n)) / den
        half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
        lo, hi = 100 * np.clip(centre - half, 0, 1), 100 * np.clip(centre + half, 0, 1)
    bad = ~(n > 0)
    return np.where(bad, np.nan, lo), np.where(bad, np.nan, hi)


def spatial_blocks(xy_um, block_um):
    """Integer label of the square spatial block (side block_um) of every point."""
    xy = np.asarray(xy_um, float)
    ij = np.floor((xy - np.nanmin(xy, axis=0)) / float(block_um)).astype(np.int64)
    return pd.MultiIndex.from_arrays([ij[:, 0], ij[:, 1]]).factorize()[0]


def block_bootstrap(stat, blocks, n_boot=500, seed=0, alpha=0.05):
    """Spatial block bootstrap: resample whole blocks with replacement (as many as there are), concatenate their
    points and evaluate stat(indices) -> float or 1-D array. Returns (lo, hi, draws): the percentile interval and the
    bootstrap values (NaN interval with fewer than two blocks)."""
    blocks = np.asarray(blocks)
    uniq, inv = np.unique(blocks, return_inverse=True)
    members = [np.flatnonzero(inv == i) for i in range(len(uniq))]
    est = np.asarray(stat(np.arange(len(blocks))), float)
    if len(members) < 2:
        nan = np.full(est.shape, np.nan)
        return nan, nan.copy(), np.empty((0,) + est.shape)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(int(n_boot)):
        pick = rng.integers(0, len(members), len(members))
        draws.append(np.asarray(stat(np.concatenate([members[i] for i in pick])), float))
    draws = np.asarray(draws)
    lo = np.nanpercentile(draws, 100 * alpha / 2, axis=0)
    hi = np.nanpercentile(draws, 100 * (1 - alpha / 2), axis=0)
    return lo, hi, draws


def lesion_signal(a, genes=("Krt7", "Krt18"), pool_k=8, min_pooled_umi=2):
    """Cross-captured A549 lesion: UMIs of `genes` summed over each bin and its pool_k nearest bins of the same
    section (bone_analysis.spatial_neighbor_matrix), and the bins where that sum reaches min_pooled_umi.
    Returns (pooled UMI per bin, in-lesion flag per bin, genes found in the data)."""
    import bone_analysis as ba
    found = [g for g in genes if g in a.var_names]
    if not found:
        return np.zeros(a.n_obs), np.zeros(a.n_obs, dtype=bool), found
    X = a.layers["counts"][:, a.var_names.get_indexer(found)]
    per_bin = np.asarray(X.sum(axis=1), dtype=float).ravel()
    pooled = np.asarray(ba.spatial_neighbor_matrix(a, pool_k) @ per_bin, dtype=float).ravel()
    return pooled, pooled >= float(min_pooled_umi), found


def fold(numerator, denominator):
    """numerator / denominator, NaN where the denominator is not positive."""
    n, d = np.asarray(numerator, float), np.asarray(denominator, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(d > 0, n / d, np.nan)


def figure_entry(name, title, caption):
    """One figures/manifest.json entry of a drawn sensitivity figure (bone_figures format)."""
    return {"name": name, "group": GROUP, "title": title, "caption": caption, "png": f"{name}.png",
            "pdf": f"{name}.pdf", "error": None}


def read_csv(results, name):
    """results/<name> as a table, or None when it is absent or empty. Text columns stay text."""
    path = Path(results) / name
    if not path.is_file():
        return None
    try:
        df = pd.read_csv(path, keep_default_na=False, na_values=["", "NA", "NaN", "nan", "None"])
    except pd.errors.EmptyDataError:
        return None
    return df if len(df.columns) else None


def save_figure(fig, out, name, dpi=200):
    """Write a matplotlib Figure as PNG and PDF into `out` (bone_figures.save)."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")


def house_style():
    """matplotlib rc_context of the manuscript figures (bone_figures.HOUSE_STYLE)."""
    import matplotlib
    from bone_figures import HOUSE_STYLE
    return matplotlib.rc_context(HOUSE_STYLE)
