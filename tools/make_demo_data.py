#!/usr/bin/env python3
"""
tools/make_demo_data.py
=======================
Synthetic Visium HD data for Bone Spatial Studio: Space Ranger 3 binned outputs (square_008um/ and square_016um/) of
two mouse tibia sections, a control and a tumor-bearing one, used by the demo project of the web application and by
the tests. The data are simulated; nothing in them is measured.

Each section is drawn as a 1 um label raster of a long bone in longitudinal section, proximal end at the top:
articular cartilage, epiphysis, growth plate, chondro-osseous junction, metaphyseal trabeculae, diaphyseal marrow with
adipocytes distally, cortical walls with periosteum, and skeletal muscle with a neurovascular bundle on one side. The
tumor-bearing section (the contralateral limb, hence mirrored) adds an osteolytic lesion in the proximal-mid marrow,
with osteoclasts along its bone interface and a fibroblast/macrophage rim. Every 8 um bin mixes cell-type expression
profiles according to what it covers; counts are Poisson draws around lognormal library sizes (median about 25 UMI per
8 um bin in the control and 50 in the tumor-bearing section), and every 16 um bin is the exact sum of its four 8 um
children. The H&E-like images are rendered from the same raster, so a region drawn on them matches the expression.

Usage:  python tools/make_demo_data.py --out DIR [--size tiny|small] [--seed 0] [--force]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import shutil
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage, sparse, special

try:
    import pipeline_config
except ModuleNotFoundError:  # run as a script: the project root is not on sys.path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import pipeline_config

log = logging.getLogger(__name__)

GENERATOR_VERSION = 1
COMPLETE_MARKER = ".complete"
BIN_SIZES = ("008um", "016um")
LOWRES_SIDE = 600

# Physical layout of a section per size (um): bone length in the section, outer bone width at the diaphysis and at the
# proximal end, widest part of the muscle band, diaphyseal cortical thickness, and the extent along the bone of the
# epiphysis, growth plate, chondro-osseous junction and metaphysis. hires: longest side of tissue_hires_image.png.
SIZES = {
    "tiny": {"length": 760.0, "width_dia": 210.0, "width_prox": 280.0, "muscle": 84.0, "cortex": 28.0,
             "epiphysis": 120.0, "growth_plate": 56.0, "junction": 24.0, "metaphysis": 190.0, "margin": 64.0,
             "hires": 800},
    "small": {"length": 1800.0, "width_dia": 380.0, "width_prox": 520.0, "muscle": 150.0, "cortex": 40.0,
              "epiphysis": 200.0, "growth_plate": 72.0, "junction": 28.0, "metaphysis": 420.0, "margin": 96.0,
              "hires": 2000},
}
_PERIOSTEUM = 14.0
_ARTICULAR = 18.0

# The two sections. depth: median UMI per 8 um bin; microns_per_pixel and rotation: registration of the array onto the
# (virtual) full-resolution image, close to the two sections of the manuscript; mirror: the contralateral limb.
_SECTIONS = (
    {"id": "CTRL-1", "label": "control", "condition": "control", "color": "#4C72B0", "lesion": False,
     "depth": 25.0, "microns_per_pixel": 0.5936, "rotation_deg": 0.35, "mirror": False},
    {"id": "TUM-1", "label": "tumor-bearing", "condition": "tumor", "color": "#C44E52", "lesion": True,
     "depth": 50.0, "microns_per_pixel": 0.5555, "rotation_deg": -0.6, "mirror": True},
)

# ---- gene panel ---------------------------------------------------------------------------------------------------
# Named in the configuration but left out of the simulated probe set, so that missing-gene handling is exercised.
ABSENT_GENES = ("Bglap", "Bglap2", "Fgf23")
# mitochondrial transcripts and their relative abundance
_MT_GENES = {"mt-Co1": 25, "mt-Co3": 15, "mt-Co2": 12, "mt-Atp6": 12, "mt-Cytb": 10, "mt-Nd1": 6, "mt-Nd4": 6,
             "mt-Nd2": 5, "mt-Nd5": 3, "mt-Nd3": 2, "mt-Atp8": 2, "mt-Nd4l": 1, "mt-Nd6": 1}
_FILLER_PREFIX = "Demo"
_N_FILLER = 1500
_N_PATTERNED = 40       # filler genes with a spatial pattern of their own (signal for Moran's I)
_AMBIENT = 1.5e-5       # level of a named gene in a cell type that does not express it

# ---- cell types ---------------------------------------------------------------------------------------------------
CELL_TYPES = ("osteoclast", "osteoblast", "osteocyte", "chondrocyte", "hypertrophic_chondrocyte", "fibroblast",
              "macrophage", "granulocyte", "erythroid", "megakaryocyte", "adipocyte", "endothelium",
              "skeletal_muscle", "smooth_muscle", "nerve", "lesion")
_T = {name: i for i, name in enumerate(CELL_TYPES)}
# Marker transcripts per cell type, per mille of the nuclear transcriptome.
_MARKERS = {
    "osteoclast": {"Ctsk": 50, "Acp5": 40, "Mmp9": 30, "Atp6v0d2": 14, "Dcstamp": 10, "Ocstamp": 7, "Calcr": 7,
                   "Nfatc1": 5, "Car2": 8, "Tnfrsf11a": 4, "Csf1r": 3, "Ctss": 4, "Cd68": 3, "Spp1": 2},
    "osteoblast": {"Col1a1": 70, "Col1a2": 55, "Ibsp": 30, "Alpl": 18, "Sp7": 12, "Runx2": 12, "Spp1": 10,
                   "Pth1r": 6, "Tnfrsf11b": 3, "Tnfsf11": 1.5, "Csf1": 2, "Dmp1": 2, "Phex": 0.3, "Il6st": 1,
                   "Tgfbr1": 1, "Tgfbr2": 1},
    "osteocyte": {"Dmp1": 35, "Sost": 18, "Mepe": 14, "Pdpn": 6, "Phex": 0.4, "Col1a1": 25, "Col1a2": 18,
                  "Ibsp": 4, "Pth1r": 3, "Tnfsf11": 3, "Tnfrsf11b": 2, "Csf1": 1},
    "chondrocyte": {"Col2a1": 80, "Acan": 30, "Comp": 14, "Matn3": 10, "Sox9": 6, "Col10a1": 1, "Pth1r": 3,
                    "Mmp13": 0.5, "Il6st": 1},
    "hypertrophic_chondrocyte": {"Col10a1": 70, "Col2a1": 40, "Acan": 12, "Mmp13": 8, "Comp": 6, "Matn3": 4,
                                 "Sox9": 3, "Spp1": 4, "Alpl": 3, "Ibsp": 2, "Pth1r": 2},
    "fibroblast": {"Col3a1": 40, "Col1a1": 35, "Col1a2": 30, "Dcn": 30, "Lum": 20, "Postn": 10, "Fbn1": 6,
                   "Pdgfra": 5, "Pdgfrb": 4, "Thy1": 4, "Fap": 2, "Acta2": 3, "Serpine1": 2, "Tgfbr2": 3,
                   "Tgfbr1": 1, "Il6st": 2, "Il11ra1": 2, "Csf1": 2, "Pth1r": 1, "Tnfsf11": 0.5, "Tnfrsf11b": 1,
                   "Socs3": 1, "Tagln": 2, "Myl9": 1},
    "macrophage": {"Lyz2": 45, "Apoe": 30, "C1qa": 15, "C1qb": 15, "C1qc": 12, "Ctss": 10, "Cd68": 8, "Csf1r": 7,
                   "Adgre1": 5, "Mrc1": 6, "Trem2": 4, "Spp1": 4, "Fabp4": 2, "Tnfrsf11a": 2, "Il6ra": 1,
                   "Socs3": 1, "Tgfbr1": 1, "Mmp9": 1, "Acp5": 1},
    "granulocyte": {"S100a8": 60, "S100a9": 60, "Ngp": 30, "Camp": 20, "Ltf": 20, "Lcn2": 20, "Retnlg": 10,
                    "Mpo": 8, "Elane": 6, "Ly6g": 5, "Lyz2": 10, "Cxcr2": 3, "Il6ra": 2, "Socs3": 1},
    "erythroid": {"Hbb-bs": 100, "Hba-a1": 80, "Hba-a2": 60, "Hbb-bt": 40, "Alas2": 10, "Car2": 10, "Slc4a1": 8,
                  "Blvrb": 6, "Gypa": 4, "Dmtn": 4, "Klf1": 3},
    "megakaryocyte": {"Pf4": 50, "Ppbp": 50, "Gp1bb": 8, "Itga2b": 8, "Vwf": 6, "Gp9": 5, "Thbs1": 4},
    "adipocyte": {"Cfd": 50, "Fabp4": 40, "Adipoq": 30, "Car3": 20, "Plin1": 10, "Cidec": 8, "Retn": 6,
                  "Plin4": 6, "Lipe": 4, "Lep": 2, "Pparg": 2},
    "endothelium": {"Emcn": 12, "Pecam1": 10, "Cdh5": 10, "Cldn5": 10, "Flt1": 8, "Kdr": 6, "Vwf": 4, "Fabp4": 4,
                    "Tgfbr2": 2, "Il6st": 2, "Thbs1": 1},
    "skeletal_muscle": {"Acta1": 60, "Ckm": 60, "Mylpf": 30, "Tnnt3": 20, "Myh1": 20, "Tnnc2": 20, "Pvalb": 15,
                        "Atp2a1": 15, "Tpm2": 10, "Ttn": 10, "Car3": 10},
    "smooth_muscle": {"Acta2": 50, "Tagln": 30, "Myl9": 15, "Myh11": 15, "Cnn1": 12},
    "nerve": {"Mpz": 30, "Mbp": 15, "Pmp22": 15, "Plp1": 10, "Sox10": 3},
    # Human tumor cells are not captured by the mouse probes; what the lesion yields is dominated by conserved,
    # ubiquitous transcripts (_UBIQUITOUS_MASS).
    "lesion": {},
}
# Share of the nuclear transcriptome taken by the configuration's ubiquitous transcripts (default 0.22).
_UBIQUITOUS_MASS = {"lesion": 0.72, "skeletal_muscle": 0.12, "erythroid": 0.14, "granulocyte": 0.18}
_UBIQUITOUS_LEVEL = {"Actb": 9, "Tmsb4x": 9, "Eef1a1": 8, "Ftl1": 6, "Fth1": 6, "Tpt1": 6, "Gapdh": 4, "Tmsb10": 4,
                     "Myl6": 3, "Rack1": 3, "Uba52": 3, "Eef2": 3, "Hspa8": 3}
# RNA captured per unit area (relative), and mitochondrial share of the counts (default 0.04).
_RNA_DENSITY = {"osteoclast": 1.1, "osteoblast": 1.0, "osteocyte": 0.35, "chondrocyte": 0.6,
                "hypertrophic_chondrocyte": 0.45, "fibroblast": 0.75, "macrophage": 1.0, "granulocyte": 1.0,
                "erythroid": 0.85, "megakaryocyte": 1.2, "adipocyte": 0.2, "endothelium": 0.7,
                "skeletal_muscle": 0.85, "smooth_muscle": 0.75, "nerve": 0.55, "lesion": 0.5}
_MITO_SHARE = {"skeletal_muscle": 0.09, "osteoclast": 0.06, "lesion": 0.05}
# Host genes induced around the lesion: fold change at its edge, decaying over about 100 um.
_LESION_RESPONSE = {"Tnfsf11": 5.0, "Socs3": 4.0, "Serpine1": 4.0, "Spp1": 3.0, "Csf1": 2.0, "Il11ra1": 2.0,
                    "Pth1r": 1.6, "Il6st": 1.5, "Tnfrsf11b": 0.6}

# ---- tissue classes of the label raster ---------------------------------------------------------------------------
# name: (cell-type composition, H&E colour, share of the area covered by nuclei in the image). "marrow" is mixed per
# bin from _MARROW and smooth random fields (erythroid islands, granulocyte-rich areas).
_CLASSES = {
    "background": ({}, (243, 240, 244), 0.0),
    "muscle": ({"skeletal_muscle": .88, "fibroblast": .04, "endothelium": .05, "smooth_muscle": .02,
                "macrophage": .01}, (219, 104, 138), .03),
    "perimysium": ({"fibroblast": .75, "endothelium": .15, "skeletal_muscle": .10}, (238, 178, 200), .15),
    "nerve": ({"nerve": .85, "fibroblast": .10, "endothelium": .05}, (232, 170, 194), .08),
    "vessel_wall": ({"smooth_muscle": .75, "endothelium": .15, "fibroblast": .10}, (206, 106, 146), .08),
    "vessel_lumen": ({"endothelium": .45, "erythroid": .45, "granulocyte": .10}, (196, 66, 86), .03),
    "periosteum": ({"fibroblast": .85, "osteoblast": .05, "endothelium": .05, "macrophage": .05},
                   (226, 150, 182), .12),
    "cortex": ({"osteocyte": .95, "osteoblast": .03, "endothelium": .02}, (234, 128, 166), .012),
    "trabecula": ({"osteocyte": .85, "osteoblast": .10, "endothelium": .05}, (237, 136, 172), .012),
    "osteoblast": ({"osteoblast": .85, "fibroblast": .05, "endothelium": .05, "macrophage": .05},
                   (172, 100, 172), .30),
    "osteoclast": ({"osteoclast": .92, "macrophage": .08}, (198, 116, 182), .22),
    "cartilage": ({"chondrocyte": .95, "fibroblast": .05}, (198, 172, 216), .07),
    "hypertrophic": ({"hypertrophic_chondrocyte": .85, "chondrocyte": .15}, (214, 192, 224), .05),
    "junction": ({"hypertrophic_chondrocyte": .30, "osteoblast": .25, "endothelium": .15, "granulocyte": .10,
                  "erythroid": .10, "macrophage": .05, "fibroblast": .05}, (186, 118, 170), .20),
    "marrow": ({}, (126, 70, 154), .32),
    "megakaryocyte": ({"megakaryocyte": .90, "granulocyte": .05, "erythroid": .05}, (214, 140, 188), .18),
    "adipocyte": ({"adipocyte": .95, "endothelium": .05}, (251, 249, 251), 0.0),
    "sinusoid": ({"endothelium": .60, "erythroid": .40}, (192, 70, 92), .04),
    "lesion": ({"lesion": .84, "macrophage": .07, "fibroblast": .05, "endothelium": .04}, (102, 62, 150), .42),
    "lesion_rim": ({"fibroblast": .45, "macrophage": .35, "endothelium": .10, "granulocyte": .10},
                   (186, 120, 170), .18),
}
_CODE = {name: i for i, name in enumerate(_CLASSES)}
_MARROW = {"granulocyte": .40, "erythroid": .36, "macrophage": .08, "endothelium": .07, "fibroblast": .03,
           "megakaryocyte": .01}
_NUCLEUS_RGB = np.array([58, 32, 96], np.float32)


def _type_vector(shares: dict) -> np.ndarray:
    v = np.zeros(len(CELL_TYPES))
    for name, share in shares.items():
        v[_T[name]] = share
    return v


_CLASS_MATRIX = np.array([_type_vector(c[0]) for c in _CLASSES.values()])          # classes x cell types
_PALETTE = np.array([c[1] for c in _CLASSES.values()], np.float32)
_NUCLEUS_THRESHOLD = np.array([special.ndtri(1 - c[2]) if c[2] > 0 else np.inf for c in _CLASSES.values()],
                              np.float32)
_RNA = np.array([_RNA_DENSITY[t] for t in CELL_TYPES])
_MITO = np.array([_MITO_SHARE.get(t, 0.04) for t in CELL_TYPES])


# ============================================================================
# GENES AND PROFILES
# ============================================================================

def gene_panel() -> list[str]:
    """Feature names of the simulated probe set in file order: every gene named in pipeline_config.DEFAULT_CONFIG
    except ABSENT_GENES, the mitochondrial genes and the filler genes. The same for every seed and section."""
    cfg = pipeline_config.DEFAULT_CONFIG
    named = []
    for sig in cfg["signatures"].values():
        named += [*sig.get("canonical", []), *sig.get("core", []), *sig.get("extended", []), *sig.get("excluded", {})]
    named += [*cfg["receptor_panel"], *cfg["adipocyte_genes"], *cfg["ubiquitous"], *cfg["figures"]["receptor_genes"]]
    for variant in cfg["sensitivity"]["variants"].values():
        for genes in variant.values():
            named += genes
    genes = [g for g in dict.fromkeys(named) if g not in ABSENT_GENES]
    genes += list(_MT_GENES) + [f"{_FILLER_PREFIX}{i:04d}" for i in range(1, _N_FILLER + 1)]
    return sorted(genes, key=lambda g: hashlib.md5(g.encode()).hexdigest())     # a fixed pseudo-genomic order


def _profiles(genes: list[str], rng: np.random.Generator) -> np.ndarray:
    """Expression profile (each row sums to 1) of every cell type over the genes, and a last row for the
    mitochondrial transcripts."""
    idx = {g: i for i, g in enumerate(genes)}
    ubi = [idx[g] for g in pipeline_config.DEFAULT_CONFIG["ubiquitous"] if g in idx]
    ubi_level = np.array([_UBIQUITOUS_LEVEL.get(genes[i], 1.0) for i in ubi])
    mito = [idx[g] for g in _MT_GENES]
    filler = np.array([i for g, i in idx.items() if g.startswith(_FILLER_PREFIX)])
    named = np.setdiff1d(np.arange(len(genes)), np.concatenate([np.asarray(ubi, int), np.asarray(mito, int), filler]))
    base = rng.lognormal(0.0, 1.3, filler.size)
    P = np.zeros((len(CELL_TYPES) + 1, len(genes)))
    for t, cell_type in enumerate(CELL_TYPES):
        P[t, named] = _AMBIENT
        for gene, permille in _MARKERS[cell_type].items():
            if gene in idx:
                P[t, idx[gene]] = permille / 1000.0
        u = ubi_level * rng.lognormal(0.0, 0.25, len(ubi))
        P[t, ubi] = u / u.sum() * _UBIQUITOUS_MASS.get(cell_type, 0.22)
        f = base * rng.lognormal(0.0, 0.8, filler.size)
        f[rng.choice(filler.size, 12, replace=False)] *= 8.0                    # filler genes specific to this type
        P[t, filler] = f / f.sum() * max(1.0 - P[t].sum(), 0.05)
        P[t] /= P[t].sum()
    P[-1, mito] = np.array(list(_MT_GENES.values()), float) / sum(_MT_GENES.values())
    return P


# ============================================================================
# TISSUE LAYOUT
# ============================================================================

def _smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _smooth_noise(rng: np.random.Generator, shape, sigma) -> np.ndarray:
    """Gaussian random field with unit standard deviation (sigma in pixels, per axis or one for all)."""
    f = ndimage.gaussian_filter(rng.standard_normal(shape, dtype=np.float32), sigma, truncate=3.0)
    return f / max(float(f.std()), 1e-6)


def _ellipse(shape, cy, cx, ry, rx):
    """Window (row slice, column slice) of a raster and the boolean mask of an axis-aligned ellipse inside it."""
    y0, y1 = max(int(cy - ry) - 1, 0), min(int(cy + ry) + 2, shape[0])
    x0, x1 = max(int(cx - rx) - 1, 0), min(int(cx + rx) + 2, shape[1])
    yy, xx = np.ogrid[y0:y1, x0:x1]
    return (slice(y0, y1), slice(x0, x1)), ((yy + 0.5 - cy) / ry) ** 2 + ((xx + 0.5 - cx) / rx) ** 2 <= 1.0


def _paint_cells(lab, rate, name, radius, allowed, gap, rng) -> int:
    """Paint round cells of class `name`: as many as the sum of `rate` (expected cells per pixel), centred on pixels
    drawn with probability proportional to `rate`, at least gap * (r1 + r2) apart, over pixels of the `allowed`
    classes only. Returns the number painted."""
    ys, xs = np.nonzero(rate > 0)
    n = int(round(float(rate[ys, xs].sum())))
    if n == 0:
        return 0
    m = min(ys.size, 8 * n)
    keys = np.log(rate[ys, xs].astype(float)) + rng.gumbel(size=ys.size)   # weighted draw without replacement
    order = np.argpartition(-keys, m - 1)[:m]
    order = order[np.argsort(-keys[order], kind="stable")]
    radii = rng.uniform(radius[0], radius[1], m)
    ok = np.isin(lab, [_CODE[a] for a in allowed])
    cy, cx, cr = np.empty(n), np.empty(n), np.empty(n)
    k = 0
    for i, r in zip(order, radii):
        y, x = ys[i] + 0.5, xs[i] + 0.5
        if k and np.any((cy[:k] - y) ** 2 + (cx[:k] - x) ** 2 < (gap * (cr[:k] + r)) ** 2):
            continue
        win, inside = _ellipse(lab.shape, y, x, r, r)
        view = lab[win]
        view[inside & ok[win]] = _CODE[name]
        cy[k], cx[k], cr[k] = y, x, r
        k += 1
        if k == n:
            break
    return k


def _layout(dims: dict, lesion: bool, rng: np.random.Generator):
    """Label raster of one section at 1 um per pixel (rows run along the bone, proximal end at the top) and the
    distance of every pixel to the lesion (None without one)."""
    m, L = dims["margin"], dims["length"]
    wp, wd, mus = dims["width_prox"] / 2, dims["width_dia"] / 2, dims["muscle"]
    ep, gpt, jt, met = dims["epiphysis"], dims["growth_plate"], dims["junction"], dims["metaphysis"]
    bend, ext = 0.02 * L, 40.0              # bend of the shaft; extraosseous extent of the lesion
    rows16 = int(np.ceil((L + 2 * m) / 16))
    cols16 = int(np.ceil((mus + 2 * wp + 2 * _PERIOSTEUM + bend + ext + 2 * m) / 16))
    H, W = 16 * rows16, 16 * cols16
    yy = np.arange(H, dtype=np.float32) + 0.5 - m           # distance from the proximal end of the bone, per row
    Y = yy[:, None]
    cap, met_end = 0.6 * ep, ep + gpt + jt + met

    # outline: rounded proximal cap, proximal width through the metaphysis, smooth taper, slight bend of the shaft;
    # the shaft continues beyond the distal cut, so the cut is not a bone surface
    hw = wp + (wd - wp) * _smoothstep((yy - met_end) / (0.25 * L))
    hw = np.where(yy < cap, wp * np.sqrt(np.clip(1 - ((cap - yy) / cap) ** 2, 0, 1)), hw)
    hw = np.where(yy < 0, 0.0, hw).astype(np.float32)
    xc = m + mus + _PERIOSTEUM + wp + bend * np.sin(np.pi * np.clip(yy / L, 0, 1))
    xs = np.arange(W, dtype=np.float32) + 0.5
    dx = xs[None, :] - xc[:, None]
    shaft = np.abs(dx) < hw[:, None]
    bone = shaft & (Y <= L)
    d_in = ndimage.distance_transform_edt(shaft)
    d_out = ndimage.distance_transform_edt(~bone)

    # wavy growth plate across the whole width; articular cartilage on the cap; cortex thinner proximally
    rel = (Y - ep - 0.12 * gpt * np.sin(2 * np.pi * dx / (1.1 * wp) + rng.uniform(0, 2 * np.pi))) / gpt
    plate = bone & (rel >= 0) & (rel < 1)
    ct = dims["cortex"] * (0.6 + 0.4 * _smoothstep((yy - met_end + 0.1 * L) / (0.2 * L)))
    articular = bone & (d_in < _ARTICULAR) & (Y < cap)
    shell = bone & ~plate & ~articular & (d_in < (ct + np.where(yy < cap, _ARTICULAR, 0.0))[:, None])
    cavity = bone & ~plate & ~articular & ~shell
    junction = cavity & (rel >= 1) & (rel < 1 + jt / gpt)
    depth = np.clip((Y - ep - gpt - jt) / met, 0, None)     # 0 at the junction, 1 at the end of the metaphysis

    # trabeculae: struts along the bone below the junction, thinning out distally; spongy bone in the epiphysis
    trab = cavity & ~junction & (rel >= 1) & (_smooth_noise(rng, (H, W), (14.0, 3.5)) > 0.5 + 1.7 * depth)
    trab |= cavity & (rel < 0) & (_smooth_noise(rng, (H, W), 5.0) > 0.8)

    # muscle band on one side, with perimysium, neurovascular bundles and small vessels
    mw = mus * _smoothstep((yy - 0.1 * L) / (0.2 * L)) * (1 - 0.3 * _smoothstep((yy - 0.75 * L) / (0.25 * L)))
    mw = np.where(yy <= L, mw * (1 + 0.08 * _smooth_noise(rng, (H,), 40.0)), 0.0)
    periosteum = ~bone & (d_out <= _PERIOSTEUM) & (Y > 0.5 * cap) & (Y <= L)
    muscle = (~bone & ~periosteum & (dx < 0) & (-dx < (hw + _PERIOSTEUM + mw)[:, None]) & (mw > 1.0)[:, None]
              & (Y >= 0) & (Y <= L))
    inner = xc - hw - _PERIOSTEUM
    perimysium = np.zeros_like(muscle)
    n_lines = int(mus // 34) + 1
    for k in range(n_lines):
        frac = (k + 0.5) / n_lines + rng.uniform(-0.08, 0.08)
        wave = 4.0 * np.sin(2 * np.pi * yy / rng.uniform(150, 320) + rng.uniform(0, 2 * np.pi))
        perimysium |= muscle & (np.abs(xs[None, :] - (inner - frac * mw + wave)[:, None]) < 1.2)
    nerve, wall, lumen = np.zeros_like(muscle), np.zeros_like(muscle), np.zeros_like(muscle)
    thick = np.nonzero(mw > 0.6 * mus)[0]
    for k in range(max(1, round(L / 900)) + max(1, round(L / 700))):
        r = int(rng.choice(thick))
        x = inner[r] - rng.uniform(0.3, 0.7) * mw[r]
        ry, rx = rng.uniform(35, 75), rng.uniform(9, 13)
        if k < max(1, round(L / 900)):                       # neurovascular bundle: a nerve beside an artery
            win, inside = _ellipse(muscle.shape, r, x, ry, rng.uniform(10, 14))
            nerve[win] |= inside
            x -= 2 * rx + 6
        r += rng.uniform(-10, 10)
        win, inside = _ellipse(muscle.shape, r, x, ry, rx)
        wall[win] |= inside
        win, inside = _ellipse(muscle.shape, r, x, ry - 5, rx - 5)
        lumen[win] |= inside
    nerve &= muscle
    wall &= muscle
    lumen &= muscle

    les = np.zeros_like(bone)
    d_les = None
    if lesion:
        cy, ry, rx = met_end + 0.04 * L, 0.17 * L, 0.95 * (wd - dims["cortex"])
        cx = xc[int(cy + m)] + 0.35 * wd                     # off-centre: breaches the cortex on that side
        ty, tx = (Y - cy) / ry, (xs[None, :] - cx) / rx
        theta = np.arctan2(ty, tx)
        ph = rng.uniform(0, 2 * np.pi, 3)
        edge = 1 + 0.10 * np.sin(3 * theta + ph[0]) + 0.06 * np.sin(5 * theta + ph[1]) + 0.04 * np.sin(9 * theta + ph[2])
        les = (np.hypot(ty, tx) < edge) & (bone | ((dx > 0) & (d_out < ext))) & ~plate & ~articular
        shell &= ~les
        trab &= ~les
        junction &= ~les
        periosteum &= ~les
        d_les = ndimage.distance_transform_edt(~les).astype(np.float32)

    # bone surfaces: osteoblast lining (active proximally, sparse on the endosteum, reactive around the lesion)
    mineral = shell | trab
    d_bone = ndimage.distance_transform_edt(~mineral)
    marrow = cavity & ~trab & ~junction & ~les
    thr = np.where(rel < 0, 0.3, np.where(depth < 1, -0.4 + 0.9 * depth, 1.1))
    if lesion:
        thr = thr - 1.4 * np.exp(-d_les / 60.0)
    lining = marrow & (d_bone <= 4.5) & (_smooth_noise(rng, (H, W), 18.0) > thr)
    rim = (marrow | junction | periosteum) & (d_les < 18) if lesion else np.zeros_like(bone)

    lab = np.zeros((H, W), np.uint8)
    for name, mask in (("muscle", muscle), ("perimysium", perimysium), ("nerve", nerve), ("vessel_wall", wall),
                       ("vessel_lumen", lumen), ("periosteum", periosteum), ("cortex", shell),
                       ("cartilage", articular | (plate & (rel < 0.55))), ("hypertrophic", plate & (rel >= 0.55)),
                       ("marrow", marrow), ("junction", junction), ("trabecula", trab), ("lesion", les),
                       ("lesion_rim", rim), ("osteoblast", lining)):
        lab[mask] = _CODE[name]
    lab[(lab == _CODE["marrow"]) & (d_bone > 6) & (_smooth_noise(rng, (H, W), 3.0) > 2.3)] = _CODE["sinusoid"]

    # discrete cells: adipocytes (increasingly distal), megakaryocytes, osteoclasts (dense at the junction, on the
    # metaphyseal trabeculae and along the lesion; sparse on epiphyseal and endosteal surfaces)
    is_marrow = lab == _CODE["marrow"]
    adipocytes = 0.35 / (np.pi * 11.0 ** 2) * _smoothstep((Y - met_end - 0.05 * L) / (0.45 * L))
    _paint_cells(lab, is_marrow * adipocytes, "adipocyte", (8.0, 14.0), ("marrow",), 0.9, rng)
    _paint_cells(lab, (is_marrow & (d_bone > 10)) / 9000.0, "megakaryocyte", (6.0, 9.0), ("marrow",), 1.2, rng)
    host = ("marrow", "osteoblast", "junction", "sinusoid", "lesion", "lesion_rim")
    surface = np.isin(lab, [_CODE[c] for c in host]) & (d_bone >= 1) & (d_bone <= 7)
    oc = surface * np.where(rel < 0, 1 / 2500, np.where(depth < 1, (1 - depth) / 1200 + 1 / 6000, 1 / 6000))
    oc = np.where(lab == _CODE["junction"], 1 / 700, oc)
    if lesion:
        oc = np.where(surface & (d_les < 20), 1 / 300, oc)
    _paint_cells(lab, oc, "osteoclast", (8.0, 12.0), host, 0.8, rng)
    return lab, d_les


# ============================================================================
# EXPRESSION
# ============================================================================

def _class_fractions(lab: np.ndarray) -> np.ndarray:
    """Share of each class in every 8 um bin (row-major over the 8 um grid)."""
    H, W = lab.shape
    rows8, cols8, K = H // 8, W // 8, len(_CLASSES)
    b = (np.arange(H)[:, None] // 8) * cols8 + np.arange(W)[None, :] // 8
    counts = np.bincount((b * K + lab).ravel(), minlength=rows8 * cols8 * K)
    return counts.reshape(rows8 * cols8, K) / 64.0


def _simulate(frac, rows8, cols8, section, d_les, P, gene_folds, rng):
    """Counts of the in-tissue 8 um bins (CSR, bins x genes, row-major order) and the in-tissue mask of the grid.
    gene_folds(in_tissue_index) -> (gene columns, per-bin fold changes) for the spatially modulated genes."""
    n, T = frac.shape[0], len(CELL_TYPES)
    cover = 1.0 - frac[:, _CODE["background"]]
    in_tissue = cover >= 0.5

    # cell-type composition: fixed per class, marrow mixed per bin; then a Dirichlet jitter per bin
    mix = np.tile(_type_vector(_MARROW), (n, 1))
    islands = _smooth_noise(rng, (rows8, cols8), 4.0).ravel()
    mix[:, _T["granulocyte"]] *= np.exp(0.6 * islands)
    mix[:, _T["erythroid"]] *= np.exp(-0.6 * islands)
    if d_les is not None:
        near = np.exp(-d_les / 150.0)
        mix[:, _T["macrophage"]] *= 1 + 1.5 * near
        mix[:, _T["granulocyte"]] *= 1 + 0.5 * near
    comp = frac @ _CLASS_MATRIX + frac[:, [_CODE["marrow"]]] * (mix / mix.sum(1, keepdims=True))
    comp /= np.maximum(cover, 1e-9)[:, None]
    jitter = np.where(comp > 0, rng.gamma(np.maximum(30.0 * comp, 1e-6)), 0.0)
    comp = jitter / np.maximum(jitter.sum(1, keepdims=True), 1e-12)

    # library size and mitochondrial share; a few damaged bins, more of them at the tissue edge
    lib = cover * (comp @ _RNA) * rng.lognormal(0.0, 0.45, n)
    mt = (comp @ _MITO) * rng.lognormal(0.0, 0.3, n)
    poor = in_tissue & (rng.random(n) < np.where(cover < 0.9, 0.05, 0.012))
    lib[poor] *= rng.uniform(0.05, 0.25, int(poor.sum()))
    mt[poor] = rng.uniform(0.22, 0.55, int(poor.sum()))
    mt = np.clip(mt, 0.005, 0.6)
    lib[~in_tissue] = 0.0
    mass = np.column_stack([comp * (lib * (1 - mt))[:, None], lib * mt]).reshape(rows8, cols8, T + 1)

    # lateral diffusion of transcripts into the four neighbouring bins; what lands outside the tissue is lost
    d = 0.05
    spread = (1 - 4 * d) * mass
    spread[1:] += d * mass[:-1]
    spread[:-1] += d * mass[1:]
    spread[:, 1:] += d * mass[:, :-1]
    spread[:, :-1] += d * mass[:, 1:]
    mass = spread.reshape(n, T + 1)
    mass[~in_tissue] = 0.0
    mass *= section["depth"] / np.median(mass[in_tissue].sum(1))

    idx = np.nonzero(in_tissue)[0]
    cols, folds = gene_folds(in_tissue)
    blocks = []
    for start in range(0, idx.size, 2048):
        b = idx[start:start + 2048]
        lam = mass[b] @ P
        lam[:, cols] *= folds[b]
        blocks.append(sparse.csr_matrix(rng.poisson(lam).astype(np.int32)))
    return sparse.vstack(blocks, format="csr"), in_tissue


def _aggregate(X8: sparse.csr_matrix, in8: np.ndarray, rows8: int, cols8: int):
    """16 um counts, each the exact sum of its four 8 um children, and the in-tissue mask of the 16 um grid (a 16 um
    bin is in tissue when any child is)."""
    cols16 = cols8 // 2
    r, c = np.divmod(np.nonzero(in8)[0], cols8)
    parent = (r // 2) * cols16 + c // 2
    in16 = np.zeros((rows8 // 2) * cols16, bool)
    in16[parent] = True
    row_of = np.cumsum(in16) - 1
    A = sparse.csr_matrix((np.ones(parent.size, np.int32), (row_of[parent], np.arange(parent.size))),
                          shape=(int(in16.sum()), parent.size))
    X16 = (A @ X8).tocsr()
    X16.sort_indices()
    return X16, in16


# ============================================================================
# IMAGE AND SPACE RANGER FILES
# ============================================================================

def _frame(section: dict, rows16: int, cols16: int, hires_side: int) -> dict:
    """Placement of the array on the (virtual) full-resolution image: microns per pixel, a small rotation, an
    offset, the image size and the hires/lowres scale factors."""
    mpp = section["microns_per_pixel"]
    th = np.deg2rad(section["rotation_deg"])
    cos, sin = float(np.cos(th)), float(np.sin(th))
    cx = np.array([0.0, 16.0 * cols16, 0.0, 16.0 * cols16])
    cy = np.array([0.0, 0.0, 16.0 * rows16, 16.0 * rows16])
    px, py = (cos * cx - sin * cy) / mpp, (sin * cx + cos * cy) / mpp
    pad = 160.0
    f = {"mpp": mpp, "cos": cos, "sin": sin, "ox": pad - px.min(), "oy": pad - py.min(),
         "full_w": float(px.max() - px.min() + 2 * pad), "full_h": float(py.max() - py.min() + 2 * pad)}
    longest = max(f["full_w"], f["full_h"])
    for key, side in (("hires", hires_side), ("lowres", LOWRES_SIDE)):
        f[f"{key}_scalef"] = side / longest
        f[f"{key}_size"] = (int(round(f["full_w"] * side / longest)), int(round(f["full_h"] * side / longest)))
    return f


def _to_fullres(f: dict, x_um, y_um):
    return (f["ox"] + (f["cos"] * x_um - f["sin"] * y_um) / f["mpp"],
            f["oy"] + (f["sin"] * x_um + f["cos"] * y_um) / f["mpp"])


def _to_array(f: dict, px, py):
    dx, dy = (px - f["ox"]) * f["mpp"], (py - f["oy"]) * f["mpp"]
    return f["cos"] * dx + f["sin"] * dy, -f["sin"] * dx + f["cos"] * dy


def _render_he(lab: np.ndarray, f: dict, rng: np.random.Generator) -> Image.Image:
    """H&E-like hires image of the label raster: eosin-pink bone and muscle, hematoxylin-purple marrow, pale
    cartilage, white adipocytes, a dense purple lesion, light background, nuclei and fibre texture."""
    w, h = f["hires_size"]
    s = f["hires_scalef"]
    um_px = f["mpp"] / s
    py, px = np.mgrid[0:h, 0:w].astype(np.float32)
    x_um, y_um = _to_array(f, (px + 0.5) / s, (py + 0.5) / s)
    ix, iy = np.floor(x_um).astype(np.int64), np.floor(y_um).astype(np.int64)
    inside = (ix >= 0) & (ix < lab.shape[1]) & (iy >= 0) & (iy < lab.shape[0])
    code = np.zeros((h, w), np.uint8)
    code[inside] = lab[iy[inside], ix[inside]]

    rgb = _PALETTE[code] * (1 + 0.035 * _smooth_noise(rng, (h, w), 0.7))[..., None]
    fibrous = np.isin(code, [_CODE[c] for c in ("muscle", "perimysium", "nerve", "periosteum", "cortex",
                                                  "trabecula", "vessel_wall")])
    rgb *= (1 + 0.07 * fibrous * _smooth_noise(rng, (h, w), (9.0 / um_px, 1.1 / um_px)))[..., None]
    lacunae = (np.isin(code, [_CODE["cartilage"], _CODE["hypertrophic"]])
               & (_smooth_noise(rng, (h, w), 2.4 / um_px) > 0.9))
    rgb[lacunae] = 0.55 * rgb[lacunae] + 0.45 * 250.0
    nuclei = _smooth_noise(rng, (h, w), 1.6 / um_px) > _NUCLEUS_THRESHOLD[code]
    rgb[nuclei] = 0.35 * rgb[nuclei] + 0.65 * _NUCLEUS_RGB
    fat = code == _CODE["adipocyte"]
    rgb[fat & ~ndimage.binary_erosion(fat)] = (222.0, 150.0, 184.0)
    rgb = ndimage.gaussian_filter(rgb, (0.6, 0.6, 0.0))
    return Image.fromarray(np.clip(rgb + 0.5, 0, 255).astype(np.uint8), "RGB")


def _barcodes(bin_um: float, rows: int, cols: int) -> list[str]:
    r, c = np.divmod(np.arange(rows * cols), cols)
    return [f"s_{int(bin_um):03d}um_{a:05d}_{b:05d}-1" for a, b in zip(r.tolist(), c.tolist())]


def _write_h5(path: Path, X: sparse.csr_matrix, barcodes: list[str], genes: list[str], library_id: str) -> None:
    """Feature-barcode matrix in the 10x HDF5 v3 layout (CSC by barcode: indptr per barcode, indices = feature)."""
    n = len(genes)
    gz = {"compression": "gzip", "compression_opts": 4, "shuffle": True}
    with h5py.File(path, "w") as f:
        f.attrs["filetype"] = "matrix"
        f.attrs["version"] = 2
        f.attrs["library_ids"] = np.array([library_id.encode()])
        f.attrs["original_gem_groups"] = np.array([1], np.int64)
        f.attrs["chemistry_description"] = "Visium HD v1 (synthetic)"
        f.attrs["software_version"] = f"tools/make_demo_data.py {GENERATOR_VERSION}"
        g = f.create_group("matrix")
        g.create_dataset("barcodes", data=np.array([b.encode() for b in barcodes]), **gz)
        g.create_dataset("data", data=X.data.astype(np.int32), **gz)
        g.create_dataset("indices", data=X.indices.astype(np.int64), **gz)
        g.create_dataset("indptr", data=X.indptr.astype(np.int64), **gz)
        g.create_dataset("shape", data=np.array([n, X.shape[0]], np.int32))
        ft = g.create_group("features")
        ft.create_dataset("_all_tag_keys", data=np.array([b"genome"]))
        ft.create_dataset("feature_type", data=np.array([b"Gene Expression"] * n), **gz)
        ft.create_dataset("genome", data=np.array([b"mm10"] * n), **gz)
        ft.create_dataset("id", data=np.array([f"ENSMUSG9{i:010d}".encode() for i in range(1, n + 1)]), **gz)
        ft.create_dataset("name", data=np.array([g_.encode() for g_ in genes]), **gz)


def _write_positions(path: Path, barcodes: list[str], bin_um: float, rows: int, cols: int, in_tissue: np.ndarray,
                     f: dict) -> None:
    r, c = np.divmod(np.arange(rows * cols), cols)
    px, py = _to_fullres(f, (c + 0.5) * bin_um, (r + 0.5) * bin_um)
    pd.DataFrame({"barcode": barcodes, "in_tissue": in_tissue.astype(np.int32), "array_row": r.astype(np.int32),
                  "array_col": c.astype(np.int32), "pxl_row_in_fullres": py, "pxl_col_in_fullres": px}
                 ).to_parquet(path, index=False)


def _scalefactors(f: dict, bin_um: float) -> dict:
    return {"spot_diameter_fullres": bin_um / f["mpp"], "bin_size_um": bin_um, "microns_per_pixel": f["mpp"],
            "regist_target_img_scalef": f["hires_scalef"], "tissue_lowres_scalef": f["lowres_scalef"],
            "fiducial_diameter_fullres": 330.0 / f["mpp"], "tissue_hires_scalef": f["hires_scalef"]}


def _write_sample(folder: Path, section: dict, index: int, size: str, seed: int, genes: list[str], P: np.ndarray,
                  gene_folds) -> None:
    dims = SIZES[size]
    geo, expr, img = (np.random.default_rng([seed, index, k]) for k in range(3))
    lab, d_les = _layout(dims, section["lesion"], geo)
    if section["mirror"]:
        lab = np.ascontiguousarray(lab[:, ::-1])
        d_les = None if d_les is None else np.ascontiguousarray(d_les[:, ::-1])
    rows8, cols8 = lab.shape[0] // 8, lab.shape[1] // 8
    d_bins = None if d_les is None else d_les[4::8, 4::8].ravel()
    X8, in8 = _simulate(_class_fractions(lab), rows8, cols8, section, d_bins, P,
                        lambda in_tissue: gene_folds(in_tissue, rows8, cols8, d_bins, expr), expr)
    X16, in16 = _aggregate(X8, in8, rows8, cols8)
    frame = _frame(section, rows8 // 2, cols8 // 2, dims["hires"])
    image = _render_he(lab, frame, img)

    first = None
    for bin_size, X, in_t, rows, cols in (("008um", X8, in8, rows8, cols8), ("016um", X16, in16, rows8 // 2,
                                                                               cols8 // 2)):
        um = pipeline_config.bin_um(bin_size)
        spatial = folder / "outs" / "binned_outputs" / f"square_{bin_size}" / "spatial"
        spatial.mkdir(parents=True)
        bcs = _barcodes(um, rows, cols)
        _write_h5(spatial.parent / "filtered_feature_bc_matrix.h5", X, [b for b, t in zip(bcs, in_t) if t], genes,
                  folder.name)
        _write_positions(spatial / "tissue_positions.parquet", bcs, um, rows, cols, in_t, frame)
        (spatial / "scalefactors_json.json").write_text(json.dumps(_scalefactors(frame, um), indent=2),
                                                        encoding="utf-8")
        if first is None:
            image.save(spatial / "tissue_hires_image.png")
            image.resize(frame["lowres_size"], Image.Resampling.LANCZOS).save(spatial / "tissue_lowres_image.png")
            first = spatial
        else:
            for name in ("tissue_hires_image.png", "tissue_lowres_image.png"):
                shutil.copyfile(first / name, spatial / name)
    umi8, umi16 = np.asarray(X8.sum(1)).ravel(), np.asarray(X16.sum(1)).ravel()
    log.info("%s: %d in-tissue 8 um bins (median %.0f UMI), %d 16 um bins (median %.0f UMI)", section["id"],
             X8.shape[0], np.median(umi8), X16.shape[0], np.median(umi16))


# ============================================================================
# PUBLIC ENTRY POINTS
# ============================================================================

def _entry(out: Path, section: dict) -> dict:
    return {"id": section["id"], "label": section["label"], "condition": section["condition"],
            "color": section["color"], "path": str(out / f"DEMO_{section['id']}" / "outs" / "binned_outputs"),
            "microns_per_pixel": None}


def _up_to_date(marker: Path, stamp: dict, samples: list[dict]) -> bool:
    try:
        done = json.loads(marker.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    if not isinstance(done, dict) or any(done.get(k) != v for k, v in stamp.items()):
        return False
    return all((Path(s["path"]) / f"square_{b}" / "filtered_feature_bc_matrix.h5").is_file()
               for s in samples for b in BIN_SIZES)


def generate(out_dir: Path | str, size: str = "small", seed: int = 0, force: bool = False) -> list[dict]:
    """
    Write the two synthetic sections into out_dir (DEMO_<id>/outs/binned_outputs/square_008um/ and square_016um/,
    each with filtered_feature_bc_matrix.h5 and spatial/) and return their configuration sample entries
    [{id, label, condition, color, path, microns_per_pixel: None}] (control first). Deterministic for a seed.
    size "tiny" (about 1,000 16 um bins per section, for tests) or "small" (about 4,000, for the demo).
    Nothing is written when out_dir/.complete records the same size, seed and generator version and the matrices
    exist, unless force; .complete is written last.
    """
    if size not in SIZES:
        raise ValueError(f"size must be one of {', '.join(SIZES)}, got {size!r}")
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError(f"seed must be a non-negative integer, got {seed!r}")
    seed = int(seed)
    out = Path(out_dir).expanduser().resolve()
    samples = [_entry(out, s) for s in _SECTIONS]
    genes = gene_panel()
    stamp = {"generator_version": GENERATOR_VERSION, "size": size, "seed": seed,
             "panel": hashlib.sha1("\n".join(genes).encode()).hexdigest()[:16]}
    marker = out / COMPLETE_MARKER
    if not force and _up_to_date(marker, stamp, samples):
        log.info("demo data in %s is up to date (size %s, seed %d)", out, size, seed)
        return samples

    out.mkdir(parents=True, exist_ok=True)
    marker.unlink(missing_ok=True)
    rng = np.random.default_rng([seed, 1000])
    P = _profiles(genes, rng)
    patterned = np.sort(rng.choice([i for i, g in enumerate(genes) if g.startswith(_FILLER_PREFIX)],
                                   _N_PATTERNED, replace=False))
    amplitude = rng.uniform(0.6, 1.4, _N_PATTERNED)
    response = np.array([genes.index(g) for g in _LESION_RESPONSE if g in genes], int)
    response_fold = np.array([fold for g, fold in _LESION_RESPONSE.items() if g in genes])

    def gene_folds(in_tissue, rows8, cols8, d_les, r):
        """Patterned filler genes (smooth fields of their own) and, around the lesion, the induced host genes."""
        fields = np.column_stack([_smooth_noise(r, (rows8, cols8), s).ravel()
                                  for s in r.uniform(3.0, 12.0, _N_PATTERNED)])
        fold = np.exp(amplitude * fields)
        fold /= fold[in_tissue].mean(0)
        if d_les is None:
            return patterned, fold
        induced = 1 + (response_fold - 1) * np.exp(-d_les / 100.0)[:, None]
        return np.concatenate([patterned, response]), np.column_stack([fold, induced])

    t0 = time.perf_counter()
    for i, section in enumerate(_SECTIONS):
        folder = out / f"DEMO_{section['id']}"
        if folder.exists():
            shutil.rmtree(folder)
        _write_sample(folder, section, i, size, seed, genes, P, gene_folds)
    tmp = marker.with_name(marker.name + ".tmp")
    tmp.write_text(json.dumps({**stamp, "samples": [s["id"] for s in samples],
                               "created": time.strftime("%Y-%m-%dT%H:%M:%S")}, indent=2), encoding="utf-8")
    os.replace(tmp, marker)
    log.info("demo data (size %s, seed %d) written to %s in %.1f s", size, seed, out, time.perf_counter() - t0)
    return samples


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="output folder (created if needed)")
    ap.add_argument("--size", choices=list(SIZES), default="small", help="tiny: tests; small: demo (default)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--force", action="store_true", help="regenerate even when the folder is up to date")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    try:
        samples = generate(args.out, size=args.size, seed=args.seed, force=args.force)
    except (ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(samples, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
