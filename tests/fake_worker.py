"""
tests/fake_worker.py
====================
Stand-in for webapp.worker in the API tests: create_app(..., worker_module="tests.fake_worker") makes the job runner
start `python -u -m tests.fake_worker <run_dir>`. It reads no data and imports no analysis package. Like the real
worker it reports pipeline_config.Progress events (plan, steps, a warning, a skipped step, done or error) next to
ordinary log lines, and it writes a small results/ folder with every kind of file the API serves: config.json,
summary.json (with a NaN), manifest.json, tables listed in the manifest and one that is not, the per-bin table, a
stand-in .h5ad, figures/manifest.json and a PNG.

Switches read from the run configuration (merge_with_defaults keeps unknown keys):
    fake_sleep    seconds to sleep inside the second step (cancel and queue tests)
    fake_fail     message of an exception raised in the third step (failure tests)
    fake_log_kb   extra log output in kilobytes, with multi-byte characters (log paging tests)
"""
from __future__ import annotations

import csv
import gzip
import json
import struct
import sys
import time
import traceback
import zlib
from pathlib import Path

from pipeline_config import STEP_LABELS, Progress, bin_um, merge_with_defaults, signature_label

STEPS = ["1. Load sections and quality control", "2. Gene-set scoring and positivity", "3. Tissue assignment",
         STEP_LABELS["moran"], STEP_LABELS["figures"]]
WARNING = "signature osteocyte: 1 of 5 genes missing from the probe set (Fgf23)"
BINS_PER_SAMPLE = 40
# assigned signature -> tissue class, as the analysis writes them (a "removed" class and an ambiguous bin included)
ASSIGNMENTS = [("osteoclast", "bone lineage"), ("osteoblast", "bone lineage"), ("granulocyte", "hematopoietic marrow"),
               ("skeletal_muscle", "removed"), ("ambiguous", "unassigned")]
POSITIVITY_TABLE = "Figure4g_6a_positivity_native_depth_roi.csv"
FIGURE = "Figure4g_bone_lineage_positivity"


def write_csv(path: Path, header: list[str], rows: list[list], compress: bool = False) -> None:
    opener = gzip.open if compress else open
    with opener(path, "wt", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


def write_png(path: Path, width: int = 4, height: int = 3) -> None:
    def chunk(kind: bytes, data: bytes) -> bytes:
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data))

    pixels = b"".join(b"\x00" + bytes([200, 60, 60]) * width for _ in range(height))
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
                     + chunk(b"IDAT", zlib.compress(pixels)) + chunk(b"IEND", b""))


def per_bin_rows(cfg: dict) -> tuple[list[str], list[list]]:
    slugs = [signature_label(cfg, n).replace(" ", "_") for n in ("osteoclast", "osteoblast", "neutrophil")]
    header = ["barcode", "sample_id", "section", "x_um", "y_um", "in_region_of_interest", "total_umi",
              "assigned_signature", "tissue_class", "leiden_domain", *[f"positive_{s}" for s in slugs],
              "positive_osteoblast_pooled_depthmatched", "positive_osteoclast_pooled_depthmatched"]
    rows = []
    for s in cfg["samples"]:
        for k in range(BINS_PER_SAMPLE):
            signature, tissue = ASSIGNMENTS[k % len(ASSIGNMENTS)]
            domain = "unassigned" if k % 13 == 0 else k % 11
            rows.append([f"s_016um_{k:05d}-1", s["id"], s.get("label") or s["id"], (k % 8) * 16.0 + 0.04,
                         (k // 8) * 16.0 + 0.06, k % 3 != 0, 10 + k, signature, tissue, domain,
                         *[signature == slug for slug in slugs], k % 4 == 0, k % 5 == 0])
    return header, rows


def write_tables(results: Path, cfg: dict) -> None:
    b = f"{int(bin_um(cfg['bin_size']))}um"
    (results / "config.json").write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    header = ["signature", "sample", "roi", "n_positive", "n_bins", "pct_positive", "ci_lo", "ci_hi", "mean_score",
              "fold_vs_reference"]
    rows = []
    for i, name in enumerate(("osteoclast", "osteoblast", "osteocyte")):
        for j, s in enumerate(cfg["samples"]):
            pct = float("nan") if (i, j) == (2, 0) else round(0.4 + 1.7 * i - 0.3 * j, 3)
            rows.append([name, s["id"], "marrow_bone", 40 + 10 * i + j, 1000, pct, 0.1, 2.5, 0.01 * i, 1.0 + j])
    write_csv(results / POSITIVITY_TABLE, header, rows)
    write_csv(results / "unlisted_table.csv", ["key", "value"], [["a", 1], ["b", 2]])
    per_bin_header, per_bin = per_bin_rows(cfg)
    write_csv(results / f"per_bin_signature_calls_{b}.csv.gz", per_bin_header, per_bin, compress=True)
    (results / f"joint_{b}.h5ad").write_bytes(b"stand-in for an AnnData file")
    (results / "notes.txt").write_text("written by tests/fake_worker.py\n", encoding="utf-8")
    summary = {
        "version": 1, "bin_size": cfg["bin_size"], "reference_sample": cfg["reference_sample"],
        "depth_match_target_umi": 23,
        "samples": [{"id": s["id"], "label": s.get("label"), "condition": s.get("condition"), "color": s.get("color"),
                     "n_bins_loaded": BINS_PER_SAMPLE, "n_bins_qc": BINS_PER_SAMPLE, "pct_retained": 100.0,
                     "median_umi": 29.5, "median_genes": 21.0, "n_bins_roi": 26, "microns_per_pixel": None}
                    for s in cfg["samples"]],
        "metrics": [{"group": "Sections", "label": "median UMI per bin", "value": float("nan"), "unit": "",
                     "detail": "not computed by the fake worker"}],
        "warnings": [WARNING], "skipped": [{"step": STEPS[3], "reason": "disabled"}], "timings": [],
        "environment": {"python": sys.version.split()[0]}, "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    (results / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    manifest = {"files": [
        {"file": POSITIVITY_TABLE, "kind": "table", "title": "Positivity, native depth (ROI)",
         "description": "Figures 4g and 6a, left panels", "rows": len(rows)},
        {"file": f"per_bin_signature_calls_{b}.csv.gz", "kind": "data", "title": "Per-bin calls",
         "description": "one row per bin", "rows": len(per_bin)},
        {"file": f"joint_{b}.h5ad", "kind": "object", "title": "Joint AnnData", "description": "all sections",
         "rows": None},
        {"file": "summary.json", "kind": "data", "title": "Summary", "description": "key numbers", "rows": None},
    ]}
    (results / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def write_figures(results: Path) -> None:
    figures = results / "figures"
    figures.mkdir(exist_ok=True)
    write_png(figures / f"{FIGURE}.png")
    manifest = {"figures": [
        {"name": FIGURE, "group": "Figure 4", "title": "Figure 4g — bone-lineage positivity", "caption": "Fake.",
         "png": f"{FIGURE}.png", "pdf": f"{FIGURE}.pdf", "error": None},
        {"name": "FigureS2_spatially_variable_genes", "group": "Supplementary", "title": "Figure S2",
         "caption": "", "png": None, "pdf": None, "error": "morans_I_per_section.csv is missing"},
    ]}
    (figures / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def run(run_dir: Path, progress: Progress) -> None:
    cfg = merge_with_defaults(json.loads((run_dir / "config.json").read_text(encoding="utf-8")))
    results = run_dir / "results"
    results.mkdir(exist_ok=True)
    progress.plan(STEPS)

    progress.step(STEPS[0])
    progress.info(f"loading {len(cfg['samples'])} sections at {cfg['bin_size']} (16 µm bins)")
    print("a library writes to stderr", file=sys.stderr, flush=True)
    progress.warn(WARNING)

    progress.step(STEPS[1])
    time.sleep(float(cfg.get("fake_sleep") or 0))
    kb = int(cfg.get("fake_log_kb") or 0)
    written, i = 0, 0
    while written < kb * 1024:
        line = f"filler {i:06d} µ{'é' * 40}"
        progress.info(line)
        written += len(line.encode("utf-8")) + 1
        i += 1

    progress.step(STEPS[2])
    if cfg.get("fake_fail"):
        raise RuntimeError(str(cfg["fake_fail"]))
    write_tables(results, cfg)

    progress.skip(STEPS[3], "disabled")
    progress.step(STEPS[4])
    write_figures(results)


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 1:
        print("usage: python -m tests.fake_worker <run_dir>", file=sys.stderr)
        return 2
    progress = Progress(events=True)
    try:
        run(Path(args[0]), progress)
    except Exception as exc:
        tb = traceback.format_exc()
        print(tb, file=sys.stderr, flush=True)
        progress.error(str(exc) or type(exc).__name__, tb)
        return 1
    progress.done()
    return 0


if __name__ == "__main__":
    sys.exit(main())
