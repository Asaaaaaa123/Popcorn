"""
tests/test_figures.py
=====================
bone_figures.make_figures on synthetic results folders that follow the output contract of bone_analysis.py
(docs/ARCHITECTURE.md, section 3), for the default two-sample configuration, one sample and three samples.

make_results() fabricates every CSV the figures read, with plausible values and the edge cases the figures must
survive: a section in which the Figure 4i distance is not estimable (and which has no positive bin of the second
mapped signature), a panel gene outside the probe set, a panel signature without positivity calls, zero positives,
a single depth-matching replicate (SD undefined), a composition signature absent from the reference section, a
Wilson bound that rounding put on the wrong side of its estimate, and a sample id that pandas would read as missing.
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import bone_figures as bf  # noqa: E402
import pipeline_config as pc  # noqa: E402

SHAPES = ((3.2, 1.6), (1.2, 2.0), (2.0, 1.0), (2.6, 1.4))        # section width, height in mm
DE_CSV = "Figure4h_osteoclast_compartment_DE.csv"
DIST_CSV = "Figure4i_distance_distributions_observed_and_null.csv"
SUMMARY_4I_CSV = "Figure4i_observed_to_null_summary.csv"
GENES_TA = ["Ctsk", "Acp5", "Mmp9", "Atp6v0d2", "Spp1", "Ccl3", "Il1b", "Cxcl2", "Lyz2", "Ptprc", "Tyrobp", "Fcer1g"]
GENES_GP = ["Col2a1", "Col10a1", "Ibsp", "Mmp13", "Acan", "Sox9", "Comp", "Matn3", "Ihh", "Pth1r", "Runx2", "Sp7"]
MORAN_GENES = ["Hbb-bs", "Hba-a1", "Col1a1", "Col2a1", "S100a8", "S100a9", "Camp", "Ltf", "Ngp", "Car2", "Ctsk",
               "Acp5", "Mmp9", "Ibsp", "Bglap", "Sparc", "Col1a2", "Acan", "Col10a1", "Lyz2", "Cd74", "H2-Aa",
               "Apoe", "C1qa", "C1qb", "Ttn", "Acta1", "Ckm", "Mb", "Myl1", "Tnnt3", "Pf4", "Ppbp", "Adipoq",
               "Cfd", "Plin1", "Lpl", "Cxcl12", "Kitl", "Lepr"]


def _wilson_pct(k, n, z=1.959963984540054):
    p = k / n
    den = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return 100 * (c - h), 100 * (c + h)


def make_results(out: Path, cfg: dict, seed: int = 0, shapes=SHAPES) -> Path:
    """Writes a complete results folder for `cfg` (every input of bone_figures, contract columns and names)."""
    cfg = pc.merge_with_defaults(cfg)
    rng = np.random.default_rng(seed)
    out.mkdir(parents=True, exist_ok=True)
    pc.save_config(cfg, out / "config.json")
    smp = cfg["samples"]
    ids = [s["id"] for s in smp]
    ref = cfg["reference_sample"]
    last = ids[-1] if len(ids) > 1 else None            # not estimable in Figure 4i, no lineage_maps.b positives
    B, Q = bf._bin_tag(cfg["bin_size"]), bf._bin_tag(cfg["qc_bin_size"])
    step = pc.bin_um(cfg["bin_size"])
    sig = cfg["signatures"]
    drawn = {c["name"] for c in cfg["classes"]}

    def slug(n):
        return pc.signature_label(cfg, n).replace(" ", "_")

    # ---- per-bin table: elliptical sections on the bin grid, ROI = inner ellipse --------------------------------
    frames = []
    progs = list(cfg["programs"]) + ["ambiguous"]
    weights = np.r_[np.ones(len(progs) - 1), 6.0]
    for i, s in enumerate(smp):
        w, h = shapes[i % len(shapes)]
        gx, gy = (g.ravel() for g in np.meshgrid(np.arange(0, w * 1000, step), np.arange(0, h * 1000, step)))
        r = ((gx - w * 500) / (w * 500)) ** 2 + ((gy - h * 500) / (h * 500)) ** 2
        gx, gy, r = gx[r <= 1], gy[r <= 1], r[r <= 1]
        n = len(gx)
        assigned = rng.choice(progs, size=n, p=weights / weights.sum())
        df = pd.DataFrame({
            "barcode": [f"{s['id']}_{j:06d}" for j in range(n)], "sample_id": s["id"], "section": s["label"],
            "x_um": (gx + 250 + 37 * i).round(2), "y_um": (gy + 400).round(2), "in_region_of_interest": r <= 0.6,
            "total_umi": rng.poisson(40, n),
            "assigned_signature": [a if a == "ambiguous" else slug(a) for a in assigned],
            "tissue_class": ["unassigned" if a == "ambiguous" else
                             (sig[a]["class"] if sig[a]["class"] in drawn else "removed") for a in assigned],
            "leiden_domain": rng.choice(["0", "1", "2", "unresolved"], n)})
        for p in cfg["primary"]:
            pos = rng.random(n) < (0.02 + 0.015 * i) * (2.0 if p == "osteoclast" else 1.0)
            if s["id"] == last and p == cfg["figures"]["lineage_maps"]["b"]:
                pos[:] = False
            df[f"positive_{slug(p)}"] = pos
        for p in (cfg["uncoupling"]["a"], cfg["uncoupling"]["b"]):
            df[f"positive_{slug(p)}_pooled_depthmatched"] = rng.random(n) < 0.05
        frames.append(df)
    pb = pd.concat(frames, ignore_index=True)
    pb.to_csv(out / f"per_bin_signature_calls_{B}.csv.gz", index=False)
    n_roi = pb[pb["in_region_of_interest"]].groupby("sample_id").size().to_dict()

    # ---- Figure 3 composition (depth-matched table names neutrophil by its label) ---------------------------------
    by_slug = {slug(n): n for n in sig}
    rows = []
    for s in smp:
        counts = pb.loc[pb["sample_id"] == s["id"], "assigned_signature"].value_counts()
        tot, amb = int(counts.sum()), int(counts.get("ambiguous", 0))
        for key, nb in counts.items():
            prog = key if key == "ambiguous" else by_slug[key]
            if s["id"] == ref and prog == "adipocyte":
                continue                                   # absent from the reference section
            rows.append({"sample": s["id"], "program": prog, "n_bins": int(nb),
                         "pct_of_all_bins": round(100 * nb / tot, 3),
                         "pct_of_assigned_bins": round(100 * nb / (tot - amb), 3) if prog != "ambiguous" else np.nan,
                         "sample_ambiguous_pct": round(100 * amb / tot, 2)})
    comp = pd.DataFrame(rows)
    comp.to_csv(out / "Figure3a_program_composition_native.csv", index=False)
    comp.assign(program=comp["program"].replace({"neutrophil": "granulocyte"})).to_csv(
        out / "Figure3b_program_composition_depth_matched.csv", index=False)

    # ---- Figures 4g, 6a positivity ---------------------------------------------------------------------------------
    rows = []
    for p in cfg["primary"]:
        if p == "megakaryocyte":
            continue                                       # a panel signature without positivity calls
        for s in smp:
            sub = pb[(pb["sample_id"] == s["id"]) & pb["in_region_of_interest"]]
            n, k = len(sub), int(sub[f"positive_{slug(p)}"].sum())
            if s["id"] == ref and p == "osteocyte":
                k = 0
            lo, hi = _wilson_pct(k, n)
            pct = round(100 * k / n, 3)
            if s["id"] == ref and p == "erythroid":
                lo = pct + 0.0004                          # rounding put the lower bound above the estimate
            rows.append({"signature": p, "sample": s["id"], "roi": "in_marrow_roi", "n_positive": k, "n_bins": n,
                         "pct_positive": pct, "ci_lo": round(lo, 4), "ci_hi": round(hi, 3),
                         "mean_score": float(rng.normal(0.1, 0.05))})
    nat = pd.DataFrame(rows)
    ref_pct = nat[nat["sample"] == ref].set_index("signature")["pct_positive"].replace(0, np.nan)
    nat["fold_vs_reference"] = (nat["pct_positive"] / nat["signature"].map(ref_pct)).round(2)
    if len(smp) == 2:
        other = nat[nat["sample"] != ref].set_index("signature")["pct_positive"]
        nat["fold_tumor_over_control"] = (nat["signature"].map(other) / nat["signature"].map(ref_pct)).round(2)
    nat.to_csv(out / "Figure4g_6a_positivity_native_depth_roi.csv", index=False)
    rows = []
    for r in nat.itertuples():
        mean = 0.6 * r.pct_positive
        if r.sample == ids[-1] and r.signature == "osteoblast":
            mean = 0.0004                                  # below the axis floor, with an SD reaching above it
        single = r.sample == ids[0]                        # a single depth-matching replicate: SD undefined
        rows.append({"signature": r.signature, "sample": r.sample, "pct_positive_mean": mean,
                     "pct_positive_sd": np.nan if single else 0.2 * mean + 0.002,
                     "n_positive_mean": mean / 100 * r.n_bins, "n_bins": r.n_bins, "mean_score": r.mean_score})
    pd.DataFrame(rows).to_csv(out / "Figure4g_6a_positivity_depth_matched.csv", index=False)

    # ---- Figure 4h ------------------------------------------------------------------------------------------------
    p_ta = [1e-15, 1e-12, 3e-9, 2e-7, 4e-6, 1e-4, 2e-3, .01, .03, .2, np.nan, .6]
    p_gp = [1e-11, 5e-9, 1e-6, 3e-5, 2e-3, .02, .04, .07, .3, .5, .8, .9]
    pd.DataFrame({"group": ["tumor_associated"] * 12 + ["growth_plate"] * 12, "names": GENES_TA + GENES_GP,
                  "scores": np.r_[np.linspace(26, 1, 12), np.linspace(19, .5, 12)],
                  "logfoldchanges": rng.normal(1.5, .5, 24), "pvals": np.r_[p_ta, p_gp],
                  "pvals_adj": np.r_[p_ta, p_gp]}).to_csv(out / DE_CSV, index=False)

    # ---- Figure 4i ------------------------------------------------------------------------------------------------
    summ, dist = [], []
    for i, s in enumerate(smp):
        if s["id"] == last:
            summ.append({"sample": s["id"], "n_a": 0, "n_b": 41, "status": "not estimable (n_positive < 10)"})
            continue
        if len(smp) > 2 and i == 1:
            obs = np.zeros(60)                             # every distance 0: identical x limits before the fix
        else:
            near = np.r_[np.full(7, 6.0), np.ones(7)]
            obs = rng.choice(np.arange(0, 14) * step, size=300, p=near / near.sum())
        null = rng.choice(np.arange(0, 20) * step, size=4000)
        for src, d, scale in (("observed", obs, 1.0), ("null", null, 1 / 199)):
            v, c = np.unique(d, return_counts=True)
            dist += [{"sample": s["id"], "source": src, "distance_um": float(a), "count": float(b) * scale}
                     for a, b in zip(v, c)]
        ratio = float(obs.mean() / null.mean())
        summ.append({"sample": s["id"], "n_a": len(obs), "n_b": 120, "n_units_roi": n_roi[s["id"]],
                     "ratio_obs_null_mean": ratio, "ratio_obs_null_lo": .9 * ratio, "ratio_obs_null_hi": 1.1 * ratio,
                     "p_greater_than_null_mean": .005, "median_nn_um": float(np.median(obs)),
                     "mean_nn_um": float(obs.mean()), "null_mean_um": float(null.mean()), "status": "estimable"})
    summ = pd.DataFrame(summ)
    est = summ[summ["status"] == "estimable"].set_index("sample")["ratio_obs_null_mean"]
    if ref in est.index and est[ref] > 0:
        summ["ratio_of_ratios_vs_reference"] = summ["ratio_obs_null_mean"] / est[ref]
    summ["note"] = "raw distances are density-dominated; do not compare between sections"
    summ.to_csv(out / SUMMARY_4I_CSV, index=False)
    pd.DataFrame(dist).to_csv(out / DIST_CSV, index=False)

    # ---- Figure 6b, Supplementary Figure S4 (gene missing from the probe set, an undetected gene) ------------------
    def detection(genes, roles, missing, undetected):
        rows = []
        for g in genes:
            if g == missing:
                continue
            for s in smp:
                n = n_roi[s["id"]]
                p = 0.0 if (g == undetected and s["id"] == ref) else float(rng.uniform(.2, 4.0))
                k = int(round(p / 100 * n))
                lo, hi = _wilson_pct(k, n)
                if k == 0:
                    lo = 1e-17                             # floating-point Wilson bound of a zero count
                rows.append({"gene": g, "role": roles[g], "sample": s["id"], "roi": "in_marrow_roi",
                             "mean_raw": p / 50, "pct_detected_raw": 1.3 * p, "mean_dm": p / 60,
                             "pct_detected_dm": p, "n_positive_dm": k, "n_bins": n, "ci_lo_dm": lo,
                             "ci_hi_dm": hi, "median_umi_per_bin": 23.0, "fold_dm": np.nan})
        return pd.DataFrame(rows)

    detection(list(cfg["receptor_panel"]), cfg["receptor_panel"], "Csf1", "Tnfsf11").to_csv(
        out / "Figure6b_receptor_panel_depth_matched.csv", index=False)
    detection(cfg["adipocyte_genes"], {g: "adipocyte" for g in cfg["adipocyte_genes"]}, "Lep", "Retn").to_csv(
        out / "FigureS4_adipocyte_genes_depth_matched.csv", index=False)

    # ---- Supplementary Figure S1 ------------------------------------------------------------------------------------
    qc = []
    for i, s in enumerate(smp):
        n = 3000 + 500 * i
        umis = np.round(rng.lognormal(3.2 + .4 * i, .8, n))
        mt = rng.uniform(0, 12, n)
        mt[:5] = np.nan
        qc.append(pd.DataFrame({"sample": s["id"], "n_genes_by_counts": np.round(umis * rng.uniform(.6, .95, n)),
                                "total_counts": umis, "pct_counts_mt": mt}))
    qc = pd.concat(qc, ignore_index=True)
    qc.to_csv(out / f"FigureS1_qc_per_bin_{Q}.csv.gz", index=False)
    for tag, factor in ((Q, 1.0), (B, 2.5)):
        rows = []
        for s in smp:
            g = qc.loc[qc["sample"] == s["id"], "n_genes_by_counts"].to_numpy() * factor
            for t in cfg["qc"]["sweep_thresholds"]:
                keep = g > t
                rows.append({"sample": s["id"], "min_genes": t, "n_retained": int(keep.sum()),
                             "pct_retained": round(100 * keep.mean(), 2),
                             "median_genes": float(np.median(g[keep])) if keep.any() else np.nan})
        pd.DataFrame(rows).to_csv(out / f"qc_threshold_sweep_{tag}.csv", index=False)

    # ---- Supplementary Figures S2, S5 -------------------------------------------------------------------------------
    moran = {"gene": MORAN_GENES}
    for i, s in enumerate(smp):
        moran[f"moranI_{s['id']}"] = np.clip(np.linspace(.8, .15, len(MORAN_GENES)) + rng.normal(0, .04, 40), 0, 1)
    moran["mean_rank"] = np.arange(1, len(MORAN_GENES) + 1, dtype=float)
    pd.DataFrame(moran).to_csv(out / "FigureS2_spatially_variable_genes_morans_I.csv", index=False)
    labels = [pc.signature_label(cfg, p) for p in cfg["primary"]]
    domains = ["osteoclast-enriched (domain 0)", "erythroid-enriched (domain 1)", "granulocyte-enriched (domain 2)",
               "chondrocyte-enriched (domain 3)", "4", "unresolved", "unassigned"]
    z = rng.normal(0, 1, (len(domains), len(labels)))
    z[:, -1] = np.nan                                      # a signature constant across domains
    S5 = pd.DataFrame(z, columns=labels)
    S5.insert(0, "leiden_domain", domains)
    S5.to_csv(out / "FigureS5_domain_by_signature_zscores.csv", index=False)
    pd.DataFrame([{"unresolved_bins": 1234, "all_bins": len(pb), "percent": round(100 * 1234 / len(pb), 1)}]).to_csv(
        out / "FigureS5_unresolved_bins.csv", index=False)

    (out / "summary.json").write_text(json.dumps({
        "version": 1, "bin_size": cfg["bin_size"], "reference_sample": ref, "depth_match_target_umi": 23,
        "samples": [{"id": s["id"], "label": s["label"], "n_bins_qc": int((pb["sample_id"] == s["id"]).sum()),
                     "n_bins_roi": int(n_roi[s["id"]])} for s in smp]}), encoding="utf-8")
    return out


class Recorder(pc.Progress):
    """Progress that remembers the order of reported step and skip labels."""

    def __init__(self):
        super().__init__(events=False, stream=io.StringIO())
        self.reported = []

    def step(self, label):
        self.reported.append(label)
        super().step(label)

    def skip(self, label, reason):
        self.reported.append(label)
        super().skip(label, reason)


def expected_names(cfg):
    cfg = pc.merge_with_defaults(cfg)
    names = [f"Figure3a_{i}_{pc.slugify(s['label'])}_tissue_classes_denoise" for i, s in enumerate(cfg["samples"], 1)]
    names += ["Figure3b_composition_denoise", "Figure4a-f_bone_lineage_maps"]
    names += [p["name"] for p in cfg["figures"]["positivity_panels"]]
    names += ["Figure4h_osteoclast_compartment_DE", "Figure4i_osteoblast_osteoclast_distance",
              "Figure6b_receptor_panel", "FigureS1_quality_control", "FigureS2_spatially_variable_genes",
              "FigureS4_adipocyte_genes", "FigureS5_leiden_domain_heatmap"]
    return names


def one_sample_config():
    cfg = pc.default_config()
    cfg["samples"] = [dict(cfg["samples"][1], id="TUM-1")]
    cfg["reference_sample"] = "TUM-1"
    return cfg


def three_sample_config():
    cfg = pc.default_config()
    base = {"condition": "control", "path": "x", "microns_per_pixel": None}
    cfg["samples"] = [dict(base, id="CTRL-1", label="control", color="#4C72B0"),
                      dict(base, id="NA", label="contralateral control", color="#55A868"),
                      dict(base, id="TUM-2", label="tumor-bearing", condition="tumor", color=None)]
    cfg["reference_sample"] = "CTRL-1"
    cfg["figures"]["lineage_maps"] = {"a": "macrophage", "b": "osteoclast"}
    cfg["figures"]["positivity_panels"].append({"name": "FigureS9_stroma_positivity",
                                                "signatures": ["fibroblast", "adipocyte", "chondrocyte"]})
    return cfg


def _check_run(folder, cfg, progress, entries):
    names = expected_names(cfg)
    figs = folder / "figures"
    for n in names:
        assert (figs / f"{n}.png").stat().st_size > 5000, n
        assert (figs / f"{n}.pdf").stat().st_size > 1000, n
    manifest = json.loads((figs / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["figures"] == entries
    assert [e["name"] for e in entries] == names
    for e in entries:
        assert e["error"] is None, (e["name"], e["error"])
        assert e["png"] == f"{e['name']}.png" and e["pdf"] == f"{e['name']}.pdf"
        assert e["title"] and len(e["caption"]) > 80 and e["group"]
    assert progress.reported == bf.planned_steps(cfg)
    assert progress.warnings == []
    assert progress.skipped == []


@pytest.mark.parametrize("make_cfg", [pc.default_config, one_sample_config, three_sample_config],
                         ids=["default-2-samples", "1-sample", "3-samples"])
def test_make_figures_draws_every_figure(tmp_path, make_cfg):
    cfg = make_cfg()
    folder = make_results(tmp_path / "results", cfg, shapes=((1.0, 2.6),) if len(cfg["samples"]) == 1 else SHAPES)
    if len(cfg["samples"]) == 3:                # the last section has no 8 um bins (optional bin_size column)
        qc_path = folder / "FigureS1_qc_per_bin_8um.csv.gz"
        qc = pd.read_csv(qc_path, dtype={"sample": str}, keep_default_na=False)
        qc.assign(bin_size=np.where(qc["sample"] == cfg["samples"][-1]["id"], "016um", "008um")).to_csv(
            qc_path, index=False)
    progress = Recorder()
    entries = bf.make_figures(folder, None, progress)          # cfg from results/config.json
    _check_run(folder, cfg, progress, entries)
    comp = pd.read_csv(folder / "figures" / "Figure3b_composition_denoise.csv")
    n_keep = len(cfg["figures"]["composition_signatures"])
    assert len(comp) == n_keep * len(cfg["samples"])
    for _, sub in comp.groupby("sample_id"):
        assert sub["percent_of_retained"].sum() == pytest.approx(100, abs=0.1)
    captions = {e["name"]: e["caption"] for e in entries}
    assert "Csf1" in captions["Figure6b_receptor_panel"]               # panel gene outside the probe set
    assert "megakaryocyte" in captions["Figure6a_myeloid_positivity"]   # panel signature without calls
    if len(cfg["samples"]) > 1:
        assert "not estimable" in captions["Figure4i_osteoblast_osteoclast_distance"]
    if len(cfg["samples"]) == 1:
        assert "200 µm" in captions["Figure3a_1_tumor-bearing_tissue_classes_denoise"]   # narrow section: short bar
    if len(cfg["samples"]) == 3:
        assert "No 8 µm bins for tumor-bearing" in captions["FigureS1_quality_control"]


def test_absent_and_empty_inputs_are_skipped(tmp_path):
    cfg = pc.default_config()
    cfg["steps"]["moran"] = False
    folder = make_results(tmp_path / "results", cfg)
    (folder / "FigureS2_spatially_variable_genes_morans_I.csv").unlink()        # step disabled
    (folder / "FigureS4_adipocyte_genes_depth_matched.csv").unlink()            # step enabled, file absent
    pd.read_csv(folder / DE_CSV).head(0).to_csv(folder / DE_CSV, index=False)     # header only
    (folder / DIST_CSV).write_text("")                                            # zero bytes
    summ = pd.read_csv(folder / SUMMARY_4I_CSV)
    summ.assign(status="not estimable (n_positive < 10)").to_csv(folder / SUMMARY_4I_CSV, index=False)
    (folder / "Figure4g_6a_positivity_depth_matched.csv").unlink()              # native panel only
    per_bin = folder / "per_bin_signature_calls_16um.csv.gz"
    pd.read_csv(per_bin).drop(columns="tissue_class").to_csv(per_bin, index=False)
    qc = pd.read_csv(folder / "FigureS1_qc_per_bin_8um.csv.gz")                 # analysis-bin fallback, no mt column
    qc.drop(columns="pct_counts_mt").to_csv(folder / "FigureS1_qc_per_bin_16um.csv.gz", index=False)
    (folder / "FigureS1_qc_per_bin_8um.csv.gz").unlink()
    progress = Recorder()
    entries = bf.make_figures(folder, cfg, progress)
    assert progress.reported == bf.planned_steps(cfg)
    assert progress.warnings == []
    reasons = {s["step"]: s["reason"] for s in progress.skipped}

    def reason(fragment):
        return next(r for label, r in reasons.items() if fragment in label)
    assert "step disabled" in reason("Figure S2")
    assert "input missing" in reason("Figure S4")
    assert "no rows" in reason("Figure 4h")
    assert "not estimable" in reason("Figure 4i")
    assert "tissue_class" in reason("Figure 3a")
    assert len(reasons) == 5
    names = [e["name"] for e in entries]
    assert not any(n.startswith(("Figure3a_", "Figure4h", "Figure4i", "FigureS2", "FigureS4")) for n in names)
    assert all(e["error"] is None for e in entries)
    caps = {e["name"]: e["caption"] for e in entries}
    assert "not computed" in caps["Figure4g_bone_lineage_positivity"]
    assert "16 µm analysis bins" in caps["FigureS1_quality_control"]
    assert {"Figure4a-f_bone_lineage_maps", "Figure6b_receptor_panel", "FigureS5_leiden_domain_heatmap"} <= set(names)


def test_failing_figure_is_reported_and_the_others_still_run(tmp_path, monkeypatch):
    cfg = pc.default_config()
    folder = make_results(tmp_path / "results", cfg)
    for f in folder.glob("*"):
        if f.name not in ("config.json", DE_CSV, "FigureS5_domain_by_signature_zscores.csv"):
            f.unlink()

    def broken(*args, **kwargs):
        raise RuntimeError("synthetic failure")
    monkeypatch.setattr(bf, "figure4h", broken)
    progress = Recorder()
    entries = bf.make_figures(folder, cfg, progress)
    assert progress.reported == bf.planned_steps(cfg)
    by_name = {e["name"]: e for e in entries}
    assert set(by_name) == {"Figure4h_osteoclast_compartment_DE", "FigureS5_leiden_domain_heatmap"}
    bad = by_name["Figure4h_osteoclast_compartment_DE"]
    assert "synthetic failure" in bad["error"] and bad["png"] is None and bad["pdf"] is None
    assert by_name["FigureS5_leiden_domain_heatmap"]["error"] is None
    assert len(progress.warnings) == 1 and "Figure4h" in progress.warnings[0]


def test_figure_functions_work_on_import_without_module_state(tmp_path):
    font_size = matplotlib.rcParams["font.size"]
    cfg = pc.default_config()
    de = pd.DataFrame({"group": ["tumor_associated"] * 3, "names": ["Ctsk", "Acp5", "Mmp9"],
                       "scores": [40.0, 12.0, 3.0], "logfoldchanges": [3.0, 2.0, 1.0], "pvals": [1e-9, 1e-4, .2],
                       "pvals_adj": [1e-8, 1e-3, .4]})
    entries = bf.figure4h(de, cfg, tmp_path)
    assert (tmp_path / "Figure4h_osteoclast_compartment_DE.png").is_file()
    assert (tmp_path / "Figure4h_osteoclast_compartment_DE.pdf").is_file()
    assert entries[0]["error"] is None and "shown: 0 growth-plate, 3 tumor-associated" in entries[0]["caption"]
    assert matplotlib.rcParams["font.size"] == font_size     # the house style does not leak into the caller


def test_planned_steps_and_cli(tmp_path, capsys):
    default = bf.planned_steps(pc.default_config())
    assert bf.planned_steps({}) == default
    assert len(set(default)) == len(default)
    extras = bf.planned_steps({}, manuscript_extras=True)
    assert extras[2:] == default and "Figure 1a" in extras[0] and "Figure 2c" in extras[1]
    folder = make_results(tmp_path / "results", pc.default_config())
    for f in folder.glob("*"):
        if f.name not in ("config.json", "FigureS5_domain_by_signature_zscores.csv"):
            f.unlink()
    assert bf.main(["--results", str(folder)]) == 0
    assert "1 figures written, 0 failed" in capsys.readouterr().out
    manifest = json.loads((folder / "figures" / "manifest.json").read_text(encoding="utf-8"))
    assert [e["name"] for e in manifest["figures"]] == ["FigureS5_leiden_domain_heatmap"]
    with pytest.raises(SystemExit):
        bf.main(["--results", str(tmp_path / "missing")])


def test_bad_positivity_panels_are_ignored_with_a_warning(tmp_path):
    cfg = pc.default_config()
    cfg["figures"]["positivity_panels"] += [{"name": "Figure6b_receptor_panel", "signatures": ["osteoclast"]},
                                            "not a panel", {"name": "FigureS8 odd/name", "signatures": "osteoclast"}]
    labels = bf.planned_steps(cfg)
    assert "Draw Supplementary Figure S8 — odd name" in labels and len(labels) == len(set(labels))
    folder = tmp_path / "results"
    folder.mkdir()
    progress = Recorder()
    assert bf.make_figures(folder, cfg, progress) == []
    assert progress.reported == labels
    assert sum("ignored" in w for w in progress.warnings) == 2


def test_manuscript_extras_need_root(tmp_path):
    folder = tmp_path / "results"
    folder.mkdir()
    progress = Recorder()
    entries = bf.make_figures(folder, None, progress, root=None, manuscript_extras=True)
    assert progress.reported == bf.planned_steps({}, manuscript_extras=True)
    assert [e["name"] for e in entries] == ["Figure1a_experimental_design"]
    assert any("Figure 2c" in w for w in progress.warnings)


def test_manuscript_extras_crop_the_hires_images(tmp_path):
    from PIL import Image
    cfg = pc.default_config()
    for s, folder in zip(cfg["samples"], ("L", "R")):                 # images found through the configured paths
        s["path"] = folder
        square = tmp_path / folder / "square_016um"
        (square / "spatial").mkdir(parents=True)
        (square / "filtered_feature_bc_matrix.h5").write_bytes(b"")
        Image.new("RGB", (3000, 2400), (200, 120, 160)).save(square / "spatial" / "tissue_hires_image.png")
    folder = tmp_path / "results"
    folder.mkdir()
    progress = Recorder()
    entries = bf.make_figures(folder, cfg, progress, root=tmp_path, manuscript_extras=True)
    names = [e["name"] for e in entries]
    assert names == ["Figure1a_experimental_design", "Figure2c_top_tumor-bearing_tumor_bone_interface",
                     "Figure2c_bottom_control_chondro_osseous_junction"]
    assert all(e["error"] is None for e in entries) and progress.warnings == []
    with Image.open(folder / "figures" / "Figure2c_top_tumor-bearing_tumor_bone_interface.png") as im:
        assert im.size == (800, 600)
    assert all((folder / "figures" / f"{n}.pdf").is_file() for n in names)


@pytest.mark.parametrize("width, expected", [(6.0, 1.0), (2.5, 1.0), (2.4, 0.5), (1.2, 0.2), (0.4, 0.1), (0.2, 0.05)])
def test_scalebar_length_suits_the_section(width, expected):
    assert bf._scalebar_mm(width) == expected
