// spatial-viewer.js — the spatial viewer: every bin of a section drawn at its position (x, y in µm, y down, equal
// aspect), colored by a label, a call, a number or the summed UMI counts of genes, with yes/no layers (positivity calls,
// region of interest, "gene detected") drawn on top and the H&E image underneath. One, two or four panels side by side,
// each with its own section and coloring; pan and zoom are linked across the panels (toggle), and the position under
// the mouse shows as a crosshair in the other panels. Used by the Viewer tab of a run (runSource below) and by the data
// explorer of a project (a projectSource on the same interface).
//
// Data source: every method receives an AbortSignal and may reject with an ApiError (a 404 from samples() shows the
// empty state).
//   source = {
//     name?                        file name stem of Save PNG (run id, project name)
//     emptyText?                   text of the empty state when there is nothing to map (the emptyText option wins)
//     samples(signal)              -> [{id, label, n}]  sections in display order; n = number of bins or null
//                                     (n_bins is accepted for n)
//     bins(sampleId, signal)       -> bin payload {sample, n, x[], y[] (µm, y down), fields: {name: {type: 'number',
//                                     values[]} | {type: 'bool', values[0|1|null]} | {type: 'category', categories[],
//                                     codes[] (-1 = no value)}}, colors: {field: {category: '#hex'}}, barcodes?: []}
//     genes?: {                    optional: coloring by genes and "detected" overlays
//       search(sampleId, query, signal) -> [{gene, total_umi, n_bins_detected, pct_bins_detected}]
//       values(sampleId, genes[], signal) -> {genes_used[], genes_missing[], n, indices[], values[], total,
//                                     n_detected, max, n_unmatched?}: summed UMIs per bin in the bin payload's order
//     }
//     image?(sampleId, signal, {retry}) -> {url, px_per_um, width, height} | null: H&E under the bins, image pixel
//                                     (u, v) at (u / px_per_um, v / px_per_um) µm; retry: asked again after the image
//                                     failed to load (a source may then get a new address from its server)
//   }
// Store (owned by the caller and kept between visits, so leaving the tab and coming back restores everything):
//   layout 1 | 2 | 4; link (pan and zoom shared by the panels); view {cx, cy, scale, auto?} = the linked view (centre in
//   µm, css px per µm; auto = a fit nobody has moved since, fitted again when the sections or the sizes change);
//   sameScale (panels colored by the same genes or number share one color range: the smallest low end, the largest
//   99th percentile); pointSize; image {on, opacity (how much of the image shows through the bins)};
//   panels [{sample, field ('' | field name | '__gene__'), genes [names, summed], overlays [bool field | 'gene:a,b'],
//            dim (dim the coloring under overlays), hidden {field: [legend keys]}, views {sample: {cx, cy, scale,
//            auto?}} (views when not linked)}];
//   caches shared by the panels: samples, data (Map sample -> prepared bins), geneData (Map 'sample\ngenes' ->
//   prepared counts, the last 40 used), images (Map sample -> loaded H&E | {none: true}).
// Bins are drawn into one offscreen layer per panel (one fillRect batch per color) only when the data, the zoom or the
// style change; panning and hovering only copy the layers, so four panels of 100,000 bins or more stay smooth.
import { h, api, mount, icon, emptyState, errorPanel, loading, fmt, toast, debounce, switchToggle } from '../core.js';
import { downloadBlob } from './log-view.js';
import { formatNumber } from './table-viewer.js';
import { parseTokens } from './chips.js';

const FIELD_LABELS = {
  tissue_class: 'Tissue class',
  assigned_signature: 'Assigned signature',
  leiden_domain: 'Leiden domain',
  in_region_of_interest: 'In region of interest',
  passes_qc: 'Passes quality control',
  total_umi: 'UMIs per bin',
  n_genes: 'Genes per bin',
  pct_mt: 'Mitochondrial UMIs (%)',
};
const CATEGORY_ORDER = ['tissue_class', 'assigned_signature', 'leiden_domain'];
const NUMBER_ORDER = ['total_umi', 'n_genes', 'pct_mt'];
const DEFAULT_FIELDS = ['tissue_class', 'assigned_signature', 'leiden_domain', 'in_region_of_interest', 'total_umi'];
const LOG_FIELDS = new Set(['total_umi', 'n_genes']);
const BOOL_TEXT = { in_region_of_interest: ['Inside', 'Outside'], passes_qc: ['Passes', 'Fails'] };
// suffix of pooled calls in the per-bin table: _pooled_depthmatched (uncoupling.calls "pooled_dm") or _pooled ("pooled")
const POOLED_RE = /_pooled(_depthmatched)?$/;
// categories without a configured color (Leiden domains, ...): Tableau 20 without its grays
const PALETTE = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf',
  '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5'];
const GRAY_NAMES = { unassigned: '#dedede', ambiguous: '#dedede', removed: '#b8b8b8', unresolved: '#c4c4c4' };
const BOOL_COLORS = { in_region_of_interest: '#1f78b4', passes_qc: '#2ca02c' };
const POSITIVE_COLOR = '#d62728';
// overlays without a configured color: well separated colors (from Sasha Trubetskoy's list of 20; the pale ones, too
// close to the color of bins in several overlays, left out)
const OVERLAY_PALETTE = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#bfef45',
  '#911eb4', '#469990', '#9a6324'];
const VIRIDIS = ['#440154', '#472c7a', '#3b518b', '#2c718e', '#21908d', '#27ad81', '#5cc863', '#aadc32', '#fde725'];
const LEVELS = 64; // color steps of a numeric field (one fillRect batch each)
const NEUTRAL = 'neutral'; // group color token: the theme's neutral (negatives, missing values, genes not detected)
// group color token: bins in two or more overlays, in the color with the most contrast to what lies under the bins: the
// theme's text color on the dark plot background, the dark text color of the exported image over the H&E image (mostly
// white) and in the exported image itself
const SEVERAL = 'several';
const MISSING = '__missing__'; // legend key of bins without a value
const ZERO = '__zero__'; // legend key of bins where none of the genes is detected
const GENE_FIELD = '__gene__'; // "Color by" value of the gene coloring
const GENE_OVERLAY = 'gene:'; // overlay key prefix: 'gene:Ctsk,Acp5' = bins where the summed count is above 0
const MAX_GENES = 50;
const GENE_CACHE = 40;
const DIM = 0.72; // share of the background mixed into the coloring under overlays
const LAYOUTS = [1, 2, 4];
const LAYOUT_TITLES = { 1: 'One panel', 2: 'Two panels side by side', 4: 'Four panels (2 × 2)' };
// css px: a panel at least this wide has its legend beside the map; four panels need the height (two rows), two keep
// the width of their maps unless there is plenty
const WIDE_PANEL = { 2: 700, 4: 500 };
const MIN_SIDE = 1.5; // css px: bins never vanish when zoomed out
const HOVER_PX = 6;
const MARGIN = 160; // css px drawn around the viewport, so short pans need no redraw
const MAX_LAYER_PX = 16e6;
const EXPORT_PALETTE = { bg: '#ffffff', neutral: '#d0d5dc', text: '#18212c', muted: '#5b6675', several: '#18212c' };
const FONT = 'system-ui, -apple-system, "Segoe UI", Roboto, Arial, sans-serif';

let uidSeq = 0;
const uid = (prefix) => `${prefix}-${(uidSeq += 1)}`;
const isColor = (c) => typeof c === 'string' && /^#(?:[0-9a-f]{3}|[0-9a-f]{6})$/i.test(c.trim());
const goodView = (v) => !!v && Number.isFinite(v.cx) && Number.isFinite(v.cy) && Number.isFinite(v.scale) && v.scale > 0;
const geneText = (genes) => genes.join(' + ');
const genesKey = (genes) => genes.join(',');
const clampedMark = (scale) => (scale && scale.clamped ? '≥ ' : ''); // a scale cut below the largest value
const rangeText = (r) => (r ? `${r.lo}\n${r.hi}` : ''); // a shared color range, '' = none

/** [r, g, b] of '#rgb', '#rrggbb' or 'rgb(r, g, b)'; null for anything else. */
function parseColor(c) {
  const s = String(c || '').trim();
  let m = /^#([0-9a-f]{3})$/i.exec(s);
  if (m) return [...m[1]].map((ch) => parseInt(ch + ch, 16));
  m = /^#([0-9a-f]{6})$/i.exec(s);
  if (m) {
    const v = parseInt(m[1], 16);
    return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
  }
  m = /^rgba?\(\s*([\d.]+)[\s,]+([\d.]+)[\s,]+([\d.]+)/i.exec(s);
  return m ? [Number(m[1]), Number(m[2]), Number(m[3])] : null;
}

/** Color c mixed with bg (t = share of bg) as an opaque color; c itself when either cannot be read. */
function mixColor(c, bg, t) {
  const a = parseColor(c);
  const b = parseColor(bg);
  if (!a || !b) return c;
  const m = (k) => Math.round(a[k] + (b[k] - a[k]) * t);
  return `rgb(${m(0)}, ${m(1)}, ${m(2)})`;
}

/** True when color c is hard to tell from one of `colors` (RGB distance below 90 of 441). */
function tooClose(c, colors) {
  const a = parseColor(c);
  if (!a) return false;
  return colors.some((o) => {
    const b = parseColor(o);
    return !!b && Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]) < 90;
  });
}

/** Color at t in [0, 1] along a list of hex stops. */
function ramp(stops, t) {
  const x = Math.min(1, Math.max(0, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(x));
  const f = x - i;
  const a = parseColor(stops[i]);
  const b = parseColor(stops[i + 1]);
  return `rgb(${Math.round(a[0] + (b[0] - a[0]) * f)}, ${Math.round(a[1] + (b[1] - a[1]) * f)}, ${Math.round(a[2] + (b[2] - a[2]) * f)})`;
}
// numbers and gene counts: viridis from its quarter point, so that the lowest step (often 1 UMI, the commonest count of
// a sparse gene) is brighter than the neutral color of the bins without a count on the dark plot background, and still
// dark on the white of the exported image
const LUT_START = 0.25;
const LUT = Array.from({ length: LEVELS }, (_, k) => ramp(VIRIDIS, LUT_START + ((1 - LUT_START) * k) / (LEVELS - 1)));

export function fieldLabel(key) {
  if (FIELD_LABELS[key]) return FIELD_LABELS[key];
  if (key.startsWith('positive_')) {
    const raw = key.slice('positive_'.length);
    const m = raw.match(POOLED_RE);
    const name = m ? raw.slice(0, -m[0].length) : raw;
    const suffix = m ? (m[1] ? ' (pooled, depth-matched)' : ' (pooled)') : '';
    return `Positive: ${name.replace(/_/g, ' ')}${suffix}`;
  }
  const text = key.replace(/_/g, ' ');
  return text.charAt(0).toUpperCase() + text.slice(1);
}

function boolText(key) {
  if (BOOL_TEXT[key]) return BOOL_TEXT[key];
  return key.startsWith('positive_') ? ['Positive', 'Negative'] : ['Yes', 'No'];
}

function overlayLabel(key) {
  return key.startsWith(GENE_OVERLAY) ? `Detected: ${geneText(key.slice(GENE_OVERLAY.length).split(','))}` : fieldLabel(key);
}

/** Gene names of typed text (commas, semicolons or spaces separate), case-insensitive duplicates dropped. */
function parseGenes(text) {
  const seen = new Set();
  const out = [];
  for (const t of parseTokens(text)) {
    const k = t.toLowerCase();
    if (!seen.has(k)) { seen.add(k); out.push(t); }
  }
  return out;
}

/** Changes of a panel state from viewer.show(): {sample?, field?, genes? (array or text), overlays?}. */
function patchState(st, patch) {
  if (patch.sample !== undefined && patch.sample !== null) st.sample = String(patch.sample);
  const genes = typeof patch.genes === 'string' ? parseGenes(patch.genes)
    : Array.isArray(patch.genes) ? parseGenes(patch.genes.join(',')) : null;
  if (genes) {
    st.genes = genes.slice(0, MAX_GENES);
    st.field = GENE_FIELD;
  } else if (typeof patch.field === 'string') {
    st.field = patch.field;
  }
  if (Array.isArray(patch.overlays)) st.overlays = [...new Set(patch.overlays.map(String))];
  return st;
}

function normalizeSamples(res) {
  return (Array.isArray(res) ? res : []).filter((s) => s && s.id !== null && s.id !== undefined).map((s) => {
    const n = Number.isFinite(s.n) ? s.n : (Number.isFinite(s.n_bins) ? s.n_bins : null);
    return { ...s, id: String(s.id), label: s.label === null || s.label === undefined ? '' : String(s.label), n };
  });
}

// ------------------------------------------------------------------------------------------------ data
function toFloat32(values, n) {
  const out = new Float32Array(n);
  for (let i = 0; i < n; i += 1) {
    const v = values ? values[i] : null;
    out[i] = v === null || v === undefined || v === '' ? NaN : Number(v);
  }
  return out;
}

function orderFields(fields) {
  const rank = (k) => {
    const f = fields.get(k);
    if (f.type === 'category') {
      const i = CATEGORY_ORDER.indexOf(k);
      return [0, i < 0 ? CATEGORY_ORDER.length : i];
    }
    if (f.type === 'bool') return [1, BOOL_TEXT[k] ? 0 : (POOLED_RE.test(k) ? 2 : 1)];
    const i = NUMBER_ORDER.indexOf(k);
    return [2, i < 0 ? NUMBER_ORDER.length : i];
  };
  return [...fields.keys()].map((k, i) => ({ k, i, r: rank(k) }))
    .sort((a, b) => a.r[0] - b.r[0] || a.r[1] - b.r[1] || a.i - b.i)
    .map((e) => e.k);
}

/** Bin payload -> typed arrays, bounding box and per-field summaries. */
function prepare(raw) {
  const xs = raw && Array.isArray(raw.x) ? raw.x : [];
  const ys = raw && Array.isArray(raw.y) ? raw.y : [];
  const n = Math.min(xs.length, ys.length);
  const x = toFloat32(xs, n);
  const y = toFloat32(ys, n);
  const ok = new Uint8Array(n);
  const bbox = { x0: Infinity, y0: Infinity, x1: -Infinity, y1: -Infinity };
  let valid = 0;
  for (let i = 0; i < n; i += 1) {
    if (!Number.isFinite(x[i]) || !Number.isFinite(y[i])) continue;
    ok[i] = 1;
    valid += 1;
    if (x[i] < bbox.x0) bbox.x0 = x[i];
    if (x[i] > bbox.x1) bbox.x1 = x[i];
    if (y[i] < bbox.y0) bbox.y0 = y[i];
    if (y[i] > bbox.y1) bbox.y1 = y[i];
  }
  const fields = new Map();
  for (const [key, f] of Object.entries((raw && raw.fields) || {})) {
    if (!f || typeof f !== 'object') continue;
    if (f.type === 'category') {
      const categories = Array.isArray(f.categories) ? f.categories.map((c) => String(c)) : [];
      const src = Array.isArray(f.codes) ? f.codes : [];
      const codes = new Int32Array(n);
      for (let i = 0; i < n; i += 1) {
        const c = src[i];
        codes[i] = Number.isInteger(c) && c >= 0 && c < categories.length ? c : -1;
      }
      fields.set(key, { key, type: 'category', categories, codes });
    } else if (f.type === 'bool') {
      const src = Array.isArray(f.values) ? f.values : [];
      const values = new Uint8Array(n); // 0 no, 1 yes, 2 no value
      for (let i = 0; i < n; i += 1) {
        const v = src[i];
        values[i] = v === null || v === undefined ? 2 : (v === true || Number(v) === 1 ? 1 : 0);
      }
      fields.set(key, { key, type: 'bool', values });
    } else if (f.type === 'number') {
      const values = toFloat32(f.values, n);
      const finite = new Float32Array(n);
      let m = 0;
      let min = Infinity;
      let max = -Infinity;
      for (let i = 0; i < n; i += 1) {
        const v = values[i];
        if (!ok[i] || !Number.isFinite(v)) continue;
        finite[m] = v;
        m += 1;
        if (v < min) min = v;
        if (v > max) max = v;
      }
      // the color scale stops at the 99th percentile, as for genes: a few extreme bins do not flatten the map
      const sorted = finite.subarray(0, m).sort();
      fields.set(key, { key, type: 'number', values, min, max, q99: m ? sorted[Math.floor(0.99 * (m - 1))] : NaN });
    }
  }
  return {
    sample: raw && raw.sample, n, x, y, ok, valid, bbox, fields, order: orderFields(fields),
    colors: raw && raw.colors && typeof raw.colors === 'object' ? raw.colors : {},
    barcodes: raw && Array.isArray(raw.barcodes) && raw.barcodes.length >= n ? raw.barcodes : null,
  };
}

/** Uniform grid over the bins (CSR layout): the bins of cell c are items[start[c] .. start[c + 1]). */
function buildIndex(d, cellSize) {
  if (!d.valid) return null;
  const { x, y, ok, n, bbox } = d;
  let cell = cellSize > 0 && Number.isFinite(cellSize) ? cellSize : 1;
  let nx;
  let ny;
  for (;;) {
    nx = Math.floor((bbox.x1 - bbox.x0) / cell) + 1;
    ny = Math.floor((bbox.y1 - bbox.y0) / cell) + 1;
    if (nx * ny <= 4e6) break;
    cell *= 2;
  }
  const cells = nx * ny;
  const start = new Int32Array(cells + 1);
  const cellOf = new Int32Array(n);
  for (let i = 0; i < n; i += 1) {
    if (!ok[i]) { cellOf[i] = -1; continue; }
    const c = Math.floor((y[i] - bbox.y0) / cell) * nx + Math.floor((x[i] - bbox.x0) / cell);
    cellOf[i] = c;
    start[c + 1] += 1;
  }
  for (let c = 0; c < cells; c += 1) start[c + 1] += start[c];
  const items = new Int32Array(start[cells]);
  const fill = start.slice(0, cells);
  for (let i = 0; i < n; i += 1) {
    const c = cellOf[i];
    if (c >= 0) { items[fill[c]] = i; fill[c] += 1; }
  }
  return { x0: bbox.x0, y0: bbox.y0, cell, nx, ny, start, items };
}

/** Index of the bin nearest to (px, py) within radius r (data units), or -1. */
function nearest(g, d, px, py, r, visible, exclude = -1) {
  if (!g) return -1;
  const c0 = Math.max(0, Math.floor((px - r - g.x0) / g.cell));
  const c1 = Math.min(g.nx - 1, Math.floor((px + r - g.x0) / g.cell));
  const r0 = Math.max(0, Math.floor((py - r - g.y0) / g.cell));
  const r1 = Math.min(g.ny - 1, Math.floor((py + r - g.y0) / g.cell));
  if (c1 < c0 || r1 < r0 || (c1 - c0 + 1) * (r1 - r0 + 1) > 250000) return -1;
  const { x, y } = d;
  let best = -1;
  let bestD = r * r;
  for (let row = r0; row <= r1; row += 1) {
    for (let col = c0; col <= c1; col += 1) {
      const c = row * g.nx + col;
      for (let k = g.start[c]; k < g.start[c + 1]; k += 1) {
        const i = g.items[k];
        if (i === exclude || (visible && !visible[i])) continue;
        const dx = x[i] - px;
        const dy = y[i] - py;
        const dd = dx * dx + dy * dy;
        if (dd <= bestD) { bestD = dd; best = i; }
      }
    }
  }
  return best;
}

/** Bin pitch (µm) estimated as the median nearest-neighbor distance of up to ~300 bins. */
function estimatePitch(d) {
  const { bbox, n, ok, x, y } = d;
  const area = Math.max(1e-9, (bbox.x1 - bbox.x0) * (bbox.y1 - bbox.y0));
  const guess = Math.sqrt(area / Math.max(1, d.valid)) || 1;
  const g = buildIndex(d, guess * 2);
  const step = Math.max(1, Math.floor(n / 300));
  const dist = [];
  for (let i = 0; i < n; i += step) {
    if (!ok[i]) continue;
    const j = nearest(g, d, x[i], y[i], guess * 4, null, i);
    if (j >= 0) dist.push(Math.hypot(x[j] - x[i], y[j] - y[i]));
  }
  dist.sort((a, b) => a - b);
  return dist.length && dist[dist.length >> 1] > 0 ? dist[dist.length >> 1] : guess;
}

/**
 * Sparse gene vector -> dense counts in the bin order of the map, with the numbers of the legend: bins detecting the
 * genes, smallest and largest nonzero count, and the color clamp (99th percentile of the nonzero counts).
 */
function prepareGene(raw, genes, n) {
  if (!raw || typeof raw !== 'object') throw new Error('The server returned no gene counts.');
  if (Number.isFinite(raw.n) && raw.n !== n) {
    throw new Error(`The gene counts cover ${fmt.int(raw.n)} bins but the map has ${fmt.int(n)}; reload the page.`);
  }
  const counts = new Float32Array(n);
  const idx = Array.isArray(raw.indices) ? raw.indices : [];
  const vals = Array.isArray(raw.values) ? raw.values : [];
  for (let k = 0; k < idx.length; k += 1) {
    const i = idx[k];
    const v = Number(vals[k]);
    if (Number.isInteger(i) && i >= 0 && i < n && v > 0) counts[i] = v;
  }
  let detected = 0;
  let max = 0;
  let lo = Infinity;
  let total = 0;
  for (let i = 0; i < n; i += 1) {
    const v = counts[i];
    if (!(v > 0)) continue;
    detected += 1;
    total += v;
    if (v > max) max = v;
    if (v < lo) lo = v;
  }
  const nz = new Float32Array(detected);
  let k = 0;
  for (let i = 0; i < n; i += 1) if (counts[i] > 0) { nz[k] = counts[i]; k += 1; }
  nz.sort();
  return {
    requested: genes.slice(),
    used: Array.isArray(raw.genes_used) ? raw.genes_used.map(String) : genes.slice(),
    missing: Array.isArray(raw.genes_missing) ? raw.genes_missing.map(String) : [],
    unmatched: Number.isFinite(raw.n_unmatched) ? raw.n_unmatched : 0,
    counts, detected, total, max,
    lo: detected ? lo : 0,
    clamp: detected ? nz[Math.floor(0.99 * (detected - 1))] : 0,
  };
}

// ------------------------------------------------------------------------------------------------ styles
function categoryColors(f, colors) {
  const configured = colors[f.key] && typeof colors[f.key] === 'object' ? colors[f.key] : {};
  let next = 0;
  return f.categories.map((name) => {
    if (isColor(configured[name])) return configured[name].trim();
    const gray = GRAY_NAMES[name.toLowerCase()];
    if (gray) return gray;
    const c = PALETTE[next % PALETTE.length];
    next += 1;
    return c;
  });
}

/** Configured signature color of a positive_<slug> field (from the payload's colors), else null. */
function signatureColor(key, colors) {
  if (!key.startsWith('positive_')) return null;
  const slug = key.slice('positive_'.length).replace(POOLED_RE, '');
  const sig = colors.assigned_signature && colors.assigned_signature[slug];
  return isColor(sig) ? sig.trim() : null;
}

function boolColor(key, colors) {
  return BOOL_COLORS[key] || signatureColor(key, colors) || POSITIVE_COLOR;
}

/**
 * Bins of `d` sorted into nG groups by groupOf(i) in [0, nG): counts per group, the bin indices of every group not
 * hidden by off(g) (null for hidden ones), visible = 1 per drawn bin.
 */
function sortIntoGroups(d, nG, groupOf, off) {
  const { n, ok } = d;
  const counts = new Int32Array(nG);
  const groupIdx = new Int32Array(n);
  for (let i = 0; i < n; i += 1) {
    const g = ok[i] ? groupOf(i) : -1;
    groupIdx[i] = g;
    if (g >= 0) counts[g] += 1;
  }
  const idx = Array.from({ length: nG }, (_, g) => (off(g) ? null : new Int32Array(counts[g])));
  const fill = new Int32Array(nG);
  const visible = new Uint8Array(n);
  let shown = 0;
  for (let i = 0; i < n; i += 1) {
    const g = groupIdx[i];
    if (g < 0 || !idx[g]) continue;
    idx[g][fill[g]] = i;
    fill[g] += 1;
    visible[i] = 1;
    shown += 1;
  }
  return { counts, idx, visible, shown };
}

/** Color range {lo, hi} of a number field on its own: its smallest value to its 99th percentile; null without values. */
function numberRange(f) {
  if (!Number.isFinite(f.min)) return null;
  return { lo: f.min, hi: Number.isFinite(f.q99) && f.q99 > f.min ? f.q99 : f.max };
}

/**
 * Draw groups and legend of one field: groups = [{color, idx}] in drawing order (color may be NEUTRAL),
 * legend = [{key, label, color, count, off}], visible = 1 per drawn bin, scale = color bar of a number.
 * range {lo, hi}: the color range of a number shared with other panels (default: numberRange(f)).
 */
function buildStyle(d, f, hidden, range = null) {
  if (!f) {
    const s = sortIntoGroups(d, 1, () => 0, () => false);
    return { kind: 'none', groups: [{ color: NEUTRAL, idx: s.idx[0] }], legend: [], visible: s.visible, shown: s.shown, scale: null };
  }
  // group of every bin, its color and legend entry; the last group holds bins without a value
  let groupOf;
  let entries;
  let scale = null;
  if (f.type === 'category') {
    const colors = categoryColors(f, d.colors);
    entries = f.categories.map((name, c) => ({ key: name, label: name, color: colors[c], gray: !!GRAY_NAMES[name.toLowerCase()] }));
    groupOf = (i) => (f.codes[i] < 0 ? entries.length : f.codes[i]);
  } else if (f.type === 'bool') {
    const [yes, no] = boolText(f.key);
    entries = [
      { key: 'yes', label: yes, color: boolColor(f.key, d.colors) },
      { key: 'no', label: no, color: NEUTRAL, gray: true },
    ];
    groupOf = (i) => (f.values[i] === 2 ? 2 : 1 - f.values[i]);
  } else {
    const r = range || numberRange(f);
    const log = LOG_FIELDS.has(f.key) && !!r && r.lo >= 0;
    const tf = log ? Math.log1p : (v) => v;
    const lo = r ? tf(r.lo) : 0;
    const span = r && r.hi > r.lo ? tf(r.hi) - lo : 0;
    entries = LUT.map((color, k) => ({ key: `level-${k}`, color }));
    groupOf = (i) => {
      const v = f.values[i];
      if (!Number.isFinite(v)) return LEVELS;
      return span > 0 ? Math.min(LEVELS - 1, Math.max(0, Math.floor(((tf(v) - lo) / span) * LEVELS))) : 0;
    };
    if (r) scale = { min: r.lo, max: r.hi, log, clamped: f.max > r.hi, shared: !!range };
  }
  const nG = entries.length + 1;
  const off = (g) => (g === nG - 1 ? hidden.has(MISSING) : f.type !== 'number' && hidden.has(entries[g].key));
  const s = sortIntoGroups(d, nG, groupOf, off);
  // drawing order: missing values, then gray / negative groups, then the rest (categories: larger first, so the rarer
  // ones stay visible on top; numbers: low to high)
  const order = [nG - 1];
  const rest = entries.map((e, g) => g);
  if (f.type === 'category') {
    rest.sort((a, b) => (entries[b].gray ? 1 : 0) - (entries[a].gray ? 1 : 0) || s.counts[b] - s.counts[a]);
  } else if (f.type === 'bool') {
    rest.reverse();
  }
  order.push(...rest);
  const groups = order.filter((g) => s.idx[g] && s.idx[g].length)
    .map((g) => ({ color: g === nG - 1 ? NEUTRAL : entries[g].color, idx: s.idx[g] }));
  const legend = f.type === 'number' ? [] : entries.map((e, g) => ({ key: e.key, label: e.label, color: e.color, count: s.counts[g], off: off(g) }));
  if (s.counts[nG - 1]) legend.push({ key: MISSING, label: 'no value', color: NEUTRAL, count: s.counts[nG - 1], off: off(nG - 1) });
  return { kind: f.type, groups, legend, visible: s.visible, shown: s.shown, scale };
}

/**
 * Coloring by the summed counts of genes: 0 in the neutral color, detected bins on viridis with a log(1 + x) scale
 * from the smallest nonzero count to the 99th percentile of the nonzero counts (higher counts take the top color), or
 * over `range` {lo, hi} shared with other panels.
 */
function buildGeneStyle(d, gv, hidden, range = null) {
  const { counts } = gv;
  const r = range || { lo: gv.lo, hi: gv.clamp };
  const lo = Math.log1p(r.lo);
  const span = Math.log1p(r.hi) - lo;
  // group 0: not detected; groups 1 .. LEVELS: detected, low to high
  const groupOf = (i) => {
    const v = counts[i];
    if (!(v > 0)) return 0;
    if (!(span > 0)) return LEVELS;
    return 1 + Math.min(LEVELS - 1, Math.max(0, Math.floor(((Math.log1p(v) - lo) / span) * LEVELS)));
  };
  const off = (g) => g === 0 && hidden.has(ZERO);
  const s = sortIntoGroups(d, LEVELS + 1, groupOf, off);
  const groups = [];
  for (let g = 0; g <= LEVELS; g += 1) {
    if (s.idx[g] && s.idx[g].length) groups.push({ color: g === 0 ? NEUTRAL : LUT[g - 1], idx: s.idx[g] });
  }
  return {
    kind: 'gene', gene: gv, groups, visible: s.visible, shown: s.shown,
    scale: gv.detected ? { min: r.lo, max: r.hi, log: true, clamped: gv.max > r.hi, shared: !!range } : null,
    legend: [{ key: ZERO, label: 'Not detected (0 UMIs)', color: NEUTRAL, count: s.counts[0], off: off(0) }],
  };
}

/**
 * Overlays drawn on top of the coloring. layers = [{key, label, color, mask: (i) => bool | null (not available),
 * status}]; a bin in exactly one layer takes that layer's color, a bin in two or more the SEVERAL color.
 * -> {layers (+ count: bins in the layer, overlap included), groups, several (bins in two or more), active,
 *     member (Int16Array: first layer of every bin, -1 = none)}
 */
function buildOverlays(d, layers) {
  const { n, ok } = d;
  const member = new Int16Array(n).fill(-1);
  const many = new Uint8Array(n);
  const counts = layers.map(() => 0);
  layers.forEach((l, k) => {
    if (!l.mask) return;
    const m = l.mask;
    let c = 0;
    for (let i = 0; i < n; i += 1) {
      if (!ok[i] || !m(i)) continue;
      c += 1;
      if (member[i] < 0) member[i] = k;
      else many[i] = 1;
    }
    counts[k] = c;
  });
  const only = new Int32Array(layers.length);
  let several = 0;
  for (let i = 0; i < n; i += 1) {
    if (member[i] < 0) continue;
    if (many[i]) several += 1;
    else only[member[i]] += 1;
  }
  const idx = layers.map((l, k) => new Int32Array(only[k]));
  const both = new Int32Array(several);
  const fill = new Int32Array(layers.length);
  let fb = 0;
  for (let i = 0; i < n; i += 1) {
    const k = member[i];
    if (k < 0) continue;
    if (many[i]) { both[fb] = i; fb += 1; } else { idx[k][fill[k]] = i; fill[k] += 1; }
  }
  // larger layers first, so the rarer ones stay visible on top; bins in several layers last
  const groups = layers.map((l, k) => k).filter((k) => only[k] > 0).sort((a, b) => only[b] - only[a])
    .map((k) => ({ color: layers[k].color, idx: idx[k] }));
  if (several) groups.push({ color: SEVERAL, idx: both });
  return {
    layers: layers.map((l, k) => ({ ...l, count: l.mask ? counts[k] : null })),
    groups, several, member, active: layers.some((l) => l.mask),
  };
}

// ------------------------------------------------------------------------------------------------ drawing
function groupColor(c, palette) {
  if (c === NEUTRAL) return palette.neutral;
  if (c === SEVERAL) return palette.several || palette.text;
  return c;
}

/**
 * Bins of `groups` as squares of `side` css px; view {scale, tx, ty} in css px; ratio = device px per css px;
 * dim = share of the background mixed into every color. Square edges snap outward to whole pixels: squares stay crisp
 * and neighboring bins always meet (no seams).
 */
function drawBins(g, groups, d, view, ratio, side, palette, dim = 0) {
  const W = g.canvas.width;
  const H = g.canvas.height;
  const s = view.scale * ratio;
  const tx = view.tx * ratio;
  const ty = view.ty * ratio;
  const half = Math.max(0.5, (side * ratio) / 2);
  const { x, y } = d;
  for (const grp of groups) {
    const c = groupColor(grp.color, palette);
    g.fillStyle = dim ? mixColor(c, palette.bg, dim) : c;
    const idx = grp.idx;
    for (let k = 0; k < idx.length; k += 1) {
      const i = idx[k];
      const cx = x[i] * s + tx;
      const left = Math.floor(cx - half);
      const right = Math.ceil(cx + half);
      if (right < 0 || left > W) continue;
      const cy = y[i] * s + ty;
      const top = Math.floor(cy - half);
      const bottom = Math.ceil(cy + half);
      if (bottom < 0 || top > H) continue;
      g.fillRect(left, top, right - left, bottom - top);
    }
  }
}

/**
 * The part of the H&E image inside the css px rectangle (x0, y0, w, h) for the view t {scale, tx, ty}: image pixel u
 * is at u / pxPerUm µm, so at u * scale / pxPerUm + tx css px. Only the visible pixels are read.
 */
function drawImageUnder(g, im, t, x0, y0, w, h) {
  const k = t.scale / im.pxPerUm; // css px per image pixel
  if (!(k > 0)) return;
  const sx0 = Math.max(0, Math.floor((x0 - t.tx) / k));
  const sy0 = Math.max(0, Math.floor((y0 - t.ty) / k));
  const sx1 = Math.min(im.width, Math.ceil((x0 + w - t.tx) / k));
  const sy1 = Math.min(im.height, Math.ceil((y0 + h - t.ty) / k));
  if (sx1 <= sx0 || sy1 <= sy0) return;
  g.drawImage(im.img, sx0, sy0, sx1 - sx0, sy1 - sy0, t.tx + sx0 * k, t.ty + sy0 * k, (sx1 - sx0) * k, (sy1 - sy0) * k);
}

/** Scale bar of a "nice" length (1, 2 or 5 x 10^k µm, about 120 px long) with its bottom-left corner at (x, y). */
function drawScaleBar(g, scale, x, y, palette) {
  const raw = 120 / scale;
  const p = 10 ** Math.floor(Math.log10(raw));
  const m = raw / p;
  const um = Number(((m >= 5 ? 5 : m >= 2 ? 2 : 1) * p).toPrecision(3));
  const len = um * scale;
  const label = um >= 1000 ? `${Number((um / 1000).toPrecision(3))} mm` : `${um} µm`;
  g.save();
  g.font = `600 11px ${FONT}`;
  const boxW = Math.max(len, g.measureText(label).width) + 16;
  g.globalAlpha = 0.85;
  g.fillStyle = palette.bg;
  g.fillRect(x, y - 34, boxW, 34);
  g.globalAlpha = 1;
  g.strokeStyle = palette.text;
  g.fillStyle = palette.text;
  g.lineWidth = 2;
  g.beginPath();
  g.moveTo(x + 8, y - 13);
  g.lineTo(x + 8, y - 8);
  g.lineTo(x + 8 + len, y - 8);
  g.lineTo(x + 8 + len, y - 13);
  g.stroke();
  g.textBaseline = 'alphabetic';
  g.fillText(label, x + 8, y - 17);
  g.restore();
}

/** Dashed crosshair through (x, y) across a w x h css px canvas, with a ring at the point. */
function drawCrosshair(g, x, y, w, hgt, palette) {
  const px = Math.round(x) + 0.5;
  const py = Math.round(y) + 0.5;
  g.save();
  g.globalAlpha = 0.7;
  g.lineWidth = 1;
  g.strokeStyle = palette.text;
  g.setLineDash([4, 4]);
  g.beginPath();
  g.moveTo(0, py);
  g.lineTo(w, py);
  g.moveTo(px, 0);
  g.lineTo(px, hgt);
  g.stroke();
  g.setLineDash([]);
  g.globalAlpha = 1;
  g.beginPath();
  g.arc(x, y, 6, 0, Math.PI * 2);
  g.lineWidth = 3;
  g.strokeStyle = palette.bg;
  g.stroke();
  g.lineWidth = 1.5;
  g.strokeStyle = palette.text;
  g.stroke();
  g.restore();
}

function readTheme(el) {
  const cs = getComputedStyle(el);
  const v = (name, fallback) => cs.getPropertyValue(name).trim() || fallback;
  const text = v('--text', '#18212c');
  return {
    bg: v('--plot-bg', '#ffffff'), neutral: v('--plot-neutral', '#c9d0da'),
    text, muted: v('--text-muted', '#566273'), several: text,
  };
}

/** css px per µm that fits the box b (plus one bin) into a z.w x z.h area with a small padding. */
function fitScaleFor(b, pitch, z) {
  const pad = Math.min(24, z.w / 8, z.h / 8);
  const s = Math.min((z.w - 2 * pad) / (b.x1 - b.x0 + pitch), (z.h - 2 * pad) / (b.y1 - b.y0 + pitch));
  return Number.isFinite(s) && s > 0 ? s : 1;
}

/** The image at url, loaded (rejects with an AbortError when signal aborts). */
function fetchImage(url, signal) {
  return new Promise((resolve, reject) => {
    const img = new Image();
    const done = () => {
      img.onload = null;
      img.onerror = null;
      signal.removeEventListener('abort', onAbort);
    };
    const onAbort = () => {
      done();
      img.removeAttribute('src');
      reject(new DOMException('The request was aborted.', 'AbortError'));
    };
    img.onload = () => { done(); resolve(img); };
    img.onerror = () => { done(); reject(new Error('The H&E image could not be loaded.')); };
    if (signal.aborted) { onAbort(); return; }
    signal.addEventListener('abort', onAbort);
    img.decoding = 'async';
    img.src = url;
  });
}

// ------------------------------------------------------------------------------------------------ gene box
/**
 * Text box for gene names with suggestions: names separated by commas (or spaces); the last one is completed from
 * search(query, signal) -> [{gene, total_umi, pct_bins_detected}] (debounced, stale requests aborted; an empty query
 * lists the most abundant genes). ↑ / ↓ pick a suggestion, Enter or a click takes it and submits, Enter alone (or the
 * button) submits onSubmit(names), Esc closes the list (then calls onCancel when given).
 * Returns {el, input, set(text), focus(), destroy()}.
 */
export function geneInput({ id, label = 'Genes', placeholder = 'e.g. Ctsk, Acp5', submitLabel = 'Show', search, onSubmit, onCancel } = {}) {
  const inputId = id || uid('sv-genes');
  const listId = `${inputId}-list`;
  const input = h('input', {
    class: 'input mono', type: 'text', id: inputId, placeholder, spellcheck: false, autocomplete: 'off',
    role: 'combobox', 'aria-autocomplete': 'list', 'aria-expanded': 'false', 'aria-controls': listId,
  });
  const list = h('ul', { class: 'sv-suggest', id: listId, role: 'listbox', 'aria-label': `${label}: suggestions`, hidden: true });
  const el = h('div', { class: 'sv-gene-box' },
    h('div', { class: 'sv-gene-field' }, input, list),
    h('button', { class: 'btn small', type: 'button', onClick: () => submit() }, submitLabel),
    onCancel ? h('button', {
      class: 'btn small ghost icon-only', type: 'button', title: 'Cancel', 'aria-label': 'Cancel', onClick: () => onCancel(),
    }, icon('x', { size: 14 })) : null);
  let items = [];
  let active = -1;
  let ctrl = null;
  let seq = 0;

  const lastToken = () => {
    const m = /[^\s,;]*$/.exec(input.value);
    return m ? m[0] : '';
  };

  const lookup = debounce(async () => {
    if (!search || document.activeElement !== input) return;
    const my = ++seq;
    if (ctrl) ctrl.abort();
    const c = new AbortController();
    ctrl = c;
    try {
      const res = await search(lastToken(), c.signal);
      if (my !== seq || document.activeElement !== input) return;
      items = (Array.isArray(res) ? res : []).filter((r) => r && r.gene).slice(0, 12);
      active = -1;
      paintList();
    } catch (err) {
      // suggestions are a help only: names can still be typed
      if (my === seq && !(err && err.name === 'AbortError')) close();
    } finally {
      if (ctrl === c) ctrl = null;
    }
  }, 180);

  function meta(r) {
    const parts = [];
    if (Number.isFinite(r.total_umi)) parts.push(`${fmt.int(r.total_umi)} UMIs`);
    if (Number.isFinite(r.pct_bins_detected)) parts.push(`${fmt.num(r.pct_bins_detected, 1)}% of bins`);
    return parts.join(' · ');
  }

  function paintList() {
    if (!items.length) { close(); return; }
    mount(list, items.map((r, k) => h('li', {
      id: `${inputId}-opt-${k}`, role: 'option', class: k === active ? 'active' : null,
      'aria-selected': k === active ? 'true' : 'false',
      onMousedown: (e) => { e.preventDefault(); take(k); }, // before the input's blur
    }, h('span', { class: 'sv-suggest-gene' }, String(r.gene)), h('span', { class: 'sv-suggest-meta' }, meta(r)))));
    list.hidden = false;
    input.setAttribute('aria-expanded', 'true');
    if (active >= 0) {
      input.setAttribute('aria-activedescendant', `${inputId}-opt-${active}`);
      const li = list.children[active];
      if (li) li.scrollIntoView({ block: 'nearest' });
    } else {
      input.removeAttribute('aria-activedescendant');
    }
  }

  function close() {
    items = [];
    active = -1;
    list.hidden = true;
    mount(list);
    input.setAttribute('aria-expanded', 'false');
    input.removeAttribute('aria-activedescendant');
  }

  function stop() {
    lookup.cancel();
    seq += 1;
    if (ctrl) ctrl.abort();
    ctrl = null;
  }

  function take(k) {
    const r = items[k];
    if (!r) return;
    const tok = lastToken();
    input.value = `${input.value.slice(0, input.value.length - tok.length)}${r.gene}`;
    submit();
  }

  function submit() {
    stop();
    close();
    if (onSubmit) onSubmit(parseGenes(input.value));
  }

  input.addEventListener('input', () => lookup());
  input.addEventListener('focus', () => lookup());
  input.addEventListener('blur', () => { stop(); close(); });
  input.addEventListener('keydown', (e) => {
    const open = !list.hidden && items.length > 0;
    if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
      e.preventDefault();
      if (!open) { lookup(); return; }
      const next = active + (e.key === 'ArrowDown' ? 1 : -1);
      active = next < -1 ? items.length - 1 : (next >= items.length ? -1 : next);
      paintList();
    } else if (e.key === 'Enter') {
      e.preventDefault();
      if (open && active >= 0) take(active);
      else submit();
    } else if (e.key === 'Escape') {
      if (open) {
        e.preventDefault();
        e.stopPropagation();
        stop();
        close();
      } else if (onCancel) {
        e.preventDefault();
        onCancel();
      }
    }
  });

  return {
    el,
    input,
    set(text) { input.value = text; },
    focus() { input.focus(); },
    destroy() { stop(); },
  };
}

// ------------------------------------------------------------------------------------------------ sources
/** Data source of the Viewer tab of a run: the per-bin table of its results and the raw data of its config.json. */
export function runSource(runId) {
  const enc = encodeURIComponent(runId);
  return {
    name: runId,
    emptyText: 'This run has no per-bin table (per_bin_signature_calls), so there is nothing to map.',
    samples: (signal) => api(`/api/runs/${enc}/bins/samples`, { signal }),
    bins: (sampleId, signal) => api(`/api/runs/${enc}/bins`, { query: { sample: sampleId }, signal }),
    genes: {
      async search(sampleId, query, signal) {
        const res = await api(`/api/runs/${enc}/genes`, { query: { sample: sampleId, query, limit: 50 }, signal });
        return res && Array.isArray(res.genes) ? res.genes : [];
      },
      values: (sampleId, genes, signal) => api(`/api/runs/${enc}/gene`, { method: 'POST', body: { sample: sampleId, genes }, signal }),
    },
    async image(sampleId, signal) {
      try {
        return await api(`/api/runs/${enc}/bins/image`, { query: { sample: sampleId }, signal });
      } catch (err) {
        if (err && err.status === 404) return null;
        throw err;
      }
    },
  };
}

// ------------------------------------------------------------------------------------------------ component
function initStore(store) {
  if (!LAYOUTS.includes(store.layout)) store.layout = 1;
  if (typeof store.link !== 'boolean') store.link = true;
  if (typeof store.sameScale !== 'boolean') store.sameScale = true;
  if (!Number.isFinite(store.pointSize)) store.pointSize = 1;
  if (!store.image || typeof store.image !== 'object') store.image = { on: false, opacity: 0.5 };
  if (!Number.isFinite(store.image.opacity)) store.image.opacity = 0.5;
  if (!(store.data instanceof Map)) store.data = new Map();
  if (!(store.geneData instanceof Map)) store.geneData = new Map();
  if (!(store.images instanceof Map)) store.images = new Map();
  if (!Array.isArray(store.panels)) store.panels = [];
  if (!goodView(store.view)) store.view = null;
}

/**
 * spatialViewer({source, store, binUm, emptyText}) -> {el, destroy(), show(patch)}
 * source and store: see the top of this file. binUm: bin size in µm (estimated from the positions when missing).
 * emptyText: text of the empty state (no sections, or samples() answered 404).
 * show({panel = 0, sample?, field?, genes? (names or text; colors the panel by their summed counts), overlays?})
 * changes what a panel shows, e.g. after a click in a gene list.
 */
export function spatialViewer({ source, store, binUm: binUmOption, emptyText }) {
  initStore(store);
  const hasGenes = !!(source.genes && typeof source.genes.search === 'function' && typeof source.genes.values === 'function');
  const hasImage = typeof source.image === 'function';

  let destroyed = false;
  let samples = Array.isArray(store.samples) && store.samples.length ? store.samples : null;
  let theme = null;
  let frame = 0;
  let pointer = null; // {panel, x, y}: the position (µm) under the mouse, a crosshair in the other panels
  let lastPanel = 0; // the panel used last: its view becomes the linked one when linking is turned on
  let startCtrl = null;
  const panels = [];
  const inflight = new Map(); // request key -> {ctrl, promise}; one request per key, shared by the panels

  // ---- toolbar
  const layoutBtns = LAYOUTS.map((k) => h('button', {
    class: 'btn small', type: 'button', title: LAYOUT_TITLES[k], 'aria-label': LAYOUT_TITLES[k], 'aria-pressed': 'false',
    onClick: () => setLayout(k),
  }, String(k)));
  const linkWrap = h('div', { class: 'sv-control sv-control-auto' },
    h('span', { class: 'sv-control-label' }, 'Pan and zoom'),
    h('div', { class: 'sv-switch-row' }, switchToggle({
      checked: store.link, label: 'Linked', title: 'Pan and zoom every panel together', onChange: (on) => setLink(on),
    })));
  const scaleWrap = h('div', { class: 'sv-control sv-control-auto' },
    h('span', { class: 'sv-control-label' }, 'Color scales'),
    h('div', { class: 'sv-switch-row' }, switchToggle({
      checked: store.sameScale, label: 'Shared',
      title: 'Panels colored by the same genes or the same number use one color scale, so that their colors compare',
      onChange: (on) => {
        store.sameScale = on;
        for (const p of panels) p.restyle();
      },
    })));
  const sizeId = uid('sv-size');
  const sizeOut = h('output', { class: 'sv-size-out', for: sizeId });
  const sizeInput = h('input', {
    type: 'range', id: sizeId, min: 0.5, max: 3, step: 0.25, value: store.pointSize, class: 'sv-range',
    onInput: () => {
      store.pointSize = Number(sizeInput.value) || 1;
      sizeOut.textContent = `${store.pointSize}×`;
      for (const p of panels) p.requestRender();
    },
  });
  sizeOut.textContent = `${store.pointSize}×`;
  let imageCtl = null;
  if (hasImage) {
    const opacity = h('input', {
      type: 'range', min: 0, max: 1, step: 0.05, value: store.image.opacity, class: 'sv-range', disabled: !store.image.on,
      'aria-label': 'H&E image: how much of it shows through the bins', title: 'How much of the image shows through the bins',
      onInput: () => {
        const v = Number(opacity.value);
        store.image.opacity = Number.isFinite(v) ? v : 0.5;
        for (const p of panels) p.queuePaint();
      },
    });
    imageCtl = h('div', { class: 'sv-control sv-control-auto' },
      h('span', { class: 'sv-control-label' }, 'H&E image'),
      h('div', { class: 'sv-range-row' }, switchToggle({
        checked: store.image.on, label: 'Show', title: 'Draw the H&E image under the bins',
        onChange: (on) => {
          store.image.on = on;
          opacity.disabled = !on;
          for (const p of panels) p.imageChanged();
          prune();
        },
      }), opacity));
  }
  // a quiet button at the end of the toolbar: the viewer's controls come first
  const saveBtn = h('button', {
    class: 'btn small ghost sv-save', type: 'button', onClick: savePng, disabled: true,
    title: 'Save the panels, with their legends, as a PNG image',
  }, icon('download', { size: 14 }), 'Save PNG');
  const toolbar = h('div', { class: 'sv-toolbar' },
    h('div', { class: 'sv-control sv-control-auto' },
      h('span', { class: 'sv-control-label' }, 'Panels'),
      h('div', { class: 'btn-group sv-layout', role: 'group', 'aria-label': 'Number of panels' }, layoutBtns)),
    linkWrap,
    scaleWrap,
    h('div', { class: 'sv-control' }, h('label', { for: sizeId }, 'Point size'), h('div', { class: 'sv-range-row' }, sizeInput, sizeOut)),
    imageCtl,
    h('span', { class: 'spacer' }),
    saveBtn);
  const panelsEl = h('div', { class: 'sv-panels' });
  const el = h('div', { class: 'sv', dataset: { layout: String(store.layout) } }, toolbar, panelsEl,
    h('p', { class: 'sv-note' }, 'Drag to pan, scroll or double-click to zoom (Shift + double-click zooms out); with a map focused, the arrow keys pan, + and − zoom and 0 fits. Hover a bin for its values; click a legend entry to hide or show it, Alt + click to show only that entry.'));

  // ---- shared observers
  const resizeObserver = new ResizeObserver((entries) => {
    for (const e of entries) {
      const p = panels.find((q) => q.stage === e.target);
      if (p) p.resize();
      const q = panels.find((x) => x.el === e.target);
      if (q) q.markWidth();
    }
    if (store.link && store.view && store.view.auto) fitLinked();
  });
  const darkQuery = window.matchMedia('(prefers-color-scheme: dark)');
  const onTheme = () => {
    theme = null;
    for (const p of panels) p.requestRender();
  };
  darkQuery.addEventListener('change', onTheme);
  const themeObserver = new MutationObserver(onTheme);
  themeObserver.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

  // ---- frames
  function schedule() {
    if (!frame && !destroyed) frame = requestAnimationFrame(tick);
  }

  function tick() {
    frame = 0;
    if (destroyed) return;
    if (!theme && el.isConnected) theme = readTheme(el);
    for (const p of panels) p.step();
  }

  // ---- requests (shared by the panels; aborted once no panel waits for them)
  function request(key, fn) {
    let entry = inflight.get(key);
    if (!entry) {
      const ctrl = new AbortController();
      const mine = { ctrl, promise: null };
      inflight.set(key, mine);
      mine.promise = (async () => {
        try {
          return await fn(ctrl.signal);
        } finally {
          if (inflight.get(key) === mine) inflight.delete(key);
        }
      })();
      entry = mine;
    }
    return entry.promise;
  }

  function prune() {
    const wanted = new Set();
    for (const p of panels) for (const k of p.wants()) wanted.add(k);
    for (const [key, r] of inflight) {
      if (wanted.has(key)) continue;
      r.ctrl.abort();
      inflight.delete(key);
    }
  }

  function loadBins(sid) {
    const cached = store.data.get(sid);
    if (cached) return Promise.resolve(cached);
    return request(`bins\n${sid}`, async (signal) => {
      const d = prepare(await source.bins(sid, signal));
      d.pitch = Number.isFinite(binUmOption) && binUmOption > 0 ? binUmOption : (d.valid > 1 ? estimatePitch(d) : 16);
      d.grid = d.valid ? buildIndex(d, d.pitch * 2) : null;
      store.data.set(sid, d);
      return d;
    });
  }

  /** A gene vector of the store's cache, made the most recently used one; undefined when not cached. */
  function cachedGene(ck) {
    const gv = store.geneData.get(ck);
    if (gv) {
      store.geneData.delete(ck);
      store.geneData.set(ck, gv);
    }
    return gv;
  }

  function loadGene(sid, d, genes) {
    const ck = `${sid}\n${genesKey(genes)}`;
    const cached = cachedGene(ck);
    if (cached) return Promise.resolve(cached);
    return request(`gene\n${ck}`, async (signal) => {
      const gv = prepareGene(await source.genes.values(sid, genes, signal), genes, d.n);
      store.geneData.set(ck, gv);
      while (store.geneData.size > GENE_CACHE) store.geneData.delete(store.geneData.keys().next().value);
      return gv;
    });
  }

  /** The H&E image of a section; retry: after a failure (a source may then ask its server for a new address). */
  function loadImage(sid, retry = false) {
    const cached = store.images.get(sid);
    if (cached) return Promise.resolve(cached);
    return request(`image\n${sid}`, async (signal) => {
      const meta = await source.image(sid, signal, { retry });
      const ppu = meta ? Number(meta.px_per_um) : NaN;
      if (!meta || !meta.url || !(ppu > 0)) {
        const none = { none: true };
        store.images.set(sid, none);
        return none;
      }
      const img = await fetchImage(String(meta.url), signal);
      // a server that sent a smaller copy of the image than its stated size: scale px_per_um with it
      const k = Number(meta.width) > 0 && img.naturalWidth ? img.naturalWidth / Number(meta.width) : 1;
      const entry = { img, pxPerUm: ppu * k, width: img.naturalWidth, height: img.naturalHeight };
      store.images.set(sid, entry);
      return entry;
    });
  }

  // ---- linked views
  /** View that fits the sections of every panel into the smallest panel, or null before any section is loaded. */
  function linkedFit() {
    let b = null;
    let pitch = 0;
    for (const p of panels) {
      const d = p.data();
      if (!d || !d.valid) continue;
      b = b ? {
        x0: Math.min(b.x0, d.bbox.x0), y0: Math.min(b.y0, d.bbox.y0),
        x1: Math.max(b.x1, d.bbox.x1), y1: Math.max(b.y1, d.bbox.y1),
      } : { ...d.bbox };
      pitch = Math.max(pitch, d.pitch);
    }
    if (!b) return null;
    let s = Infinity;
    for (const p of panels) {
      const z = p.size();
      if (z.w && z.h) s = Math.min(s, fitScaleFor(b, pitch, z));
    }
    return Number.isFinite(s) ? { cx: (b.x0 + b.x1) / 2, cy: (b.y0 + b.y1) / 2, scale: s } : null;
  }

  function fitLinked() {
    const v = linkedFit();
    if (!v) return null;
    store.view = { ...v, auto: true };
    for (const p of panels) p.requestRender();
    return store.view;
  }

  /** The view of `from` moved: copy the layers now (all panels when linked); a zoom also redraws once it settles. */
  function viewMoved(from, zooming) {
    for (const p of store.link ? panels : [from]) {
      p.queuePaint();
      if (zooming) p.crisp();
    }
  }

  function viewSettled(from) {
    for (const p of store.link ? panels : [from]) p.requestRender();
  }

  function setPointer(next) {
    if (!next && !pointer) return;
    const from = next ? next.panel : pointer.panel;
    pointer = next;
    if (store.link) for (const p of panels) if (p.index !== from) p.queuePaint();
  }

  function setLink(on) {
    if (on === store.link) return;
    setPointer(null);
    if (on) {
      const src = panels[lastPanel] || panels[0];
      const v = src ? src.ownView() : null;
      store.link = true;
      store.view = null;
      if (goodView(v) && !v.auto) store.view = { cx: v.cx, cy: v.cy, scale: v.scale };
      else fitLinked();
    } else {
      store.link = false;
      if (goodView(store.view)) for (const p of panels) p.adoptView(store.view);
    }
    for (const p of panels) p.requestRender();
  }

  // ---- shared color scales: panels colored by the same genes or the same number use the widest of their ranges
  let syncingScales = false;

  /** {lo, hi, n} over the panels whose coloring has the scale key `key` (see panel.ownScale), or null. */
  function sharedRange(key) {
    let r = null;
    for (const p of panels) {
      const s = p.ownScale();
      if (!s || s.key !== key) continue;
      r = r ? { lo: Math.min(r.lo, s.lo), hi: Math.max(r.hi, s.hi), n: r.n + 1 } : { lo: s.lo, hi: s.hi, n: 1 };
    }
    return r;
  }

  /** The range a panel colors with: shared when another panel shows the same thing, else null (its own). */
  function scaleRangeFor(own) {
    if (!own || !store.sameScale || panels.length < 2) return null;
    const r = sharedRange(own.key);
    return r && r.n > 1 ? { lo: r.lo, hi: r.hi } : null;
  }

  /** The coloring of `from` changed (or a panel came or went): restyle the panels whose shared range moved. */
  function syncScales(from) {
    if (syncingScales) return;
    syncingScales = true;
    try {
      for (const p of panels) if (p !== from) p.restyleIfScaleMoved();
    } finally {
      syncingScales = false;
    }
  }

  // ---- room on screen: the maps take the height of the window, less what has to be seen with them
  /**
   * --sv-lead: the space above the first maps once the viewer is in view: its place on the page when one row of
   * panels starts in the top half of the first screen (all of it then fits without scrolling), else its toolbar and
   * the panels' controls (the viewer scrolled to the top); --sv-row: what a second row of panels adds above its maps
   * (panel head, controls, overlays, gap), measured on that row (with four panels its rows differ when one shows a gene
   * box or overlays and the other does not).
   */
  function measureRoom() {
    const first = panels[0];
    if (destroyed || !first || !el.isConnected) return;
    if (!pageWatched) {
      const page = el.closest('.page') || el.parentElement;
      if (page) { pageObserver.observe(page); pageWatched = true; }
    }
    const top = el.getBoundingClientRect().top + window.scrollY;
    measuredTop = Math.round(top);
    const stageTop = first.stage.getBoundingClientRect().top + window.scrollY;
    if (!(stageTop > top)) return; // not laid out (hidden)
    const lead = store.layout !== 4 && top < window.innerHeight / 2 ? stageTop : stageTop - top + 12;
    const second = store.layout === 4 && panels[2] ? panels[2] : first;
    const row = second.stage.getBoundingClientRect().top - second.el.getBoundingClientRect().top + 16;
    el.style.setProperty('--sv-lead', `${Math.round(lead)}px`);
    el.style.setProperty('--sv-row', `${Math.round(row)}px`);
  }
  let roomFrame = 0;
  const measureSoon = () => {
    if (!roomFrame && !destroyed) roomFrame = requestAnimationFrame(() => { roomFrame = 0; measureRoom(); });
  };
  window.addEventListener('resize', measureSoon);
  // the page above the viewer grows or shrinks (a table filling in, a banner): measure again once the viewer moved
  let pageWatched = false;
  let measuredTop = NaN;
  const pageObserver = new ResizeObserver(() => {
    if (el.isConnected && Math.round(el.getBoundingClientRect().top + window.scrollY) !== measuredTop) measureSoon();
  });
  // the tooltip is placed in the window (so that no panel edge cuts it): hide it when the page scrolls under it
  const onScroll = () => { for (const p of panels) p.hideTip(); };
  document.addEventListener('scroll', onScroll, { capture: true, passive: true });

  // ---- panels
  function sampleText(sid) {
    const s = samples && samples.find((x) => x.id === sid);
    if (!s) return String(sid ?? '');
    return s.label && s.label !== s.id ? `${s.label} (${s.id})` : s.id;
  }

  /** Section of a new panel i: the first one that none of the panels before it shows (else they cycle). */
  function freshSample(i) {
    const ids = samples.map((s) => s.id);
    const shown = store.panels.slice(0, i).map((p) => (p && typeof p === 'object' ? p.sample : null));
    return ids.find((id) => !shown.includes(id)) ?? ids[i % ids.length];
  }

  function ensurePanelState(i) {
    const ids = samples.map((s) => s.id);
    let st = store.panels[i];
    if (!st || typeof st !== 'object') {
      // a new panel shows a section not shown yet, with the coloring of the first panel
      const base = i > 0 && store.panels[0] && typeof store.panels[0] === 'object' ? store.panels[0] : null;
      st = {
        sample: freshSample(i), field: base ? base.field : '', genes: base && Array.isArray(base.genes) ? base.genes.slice() : [],
        overlays: base && Array.isArray(base.overlays) ? base.overlays.slice() : [], dim: base ? base.dim !== false : true,
        hidden: {}, views: {},
      };
      store.panels[i] = st;
    }
    if (!ids.includes(st.sample)) st.sample = freshSample(i);
    if (typeof st.field !== 'string') st.field = '';
    st.genes = Array.isArray(st.genes) ? st.genes.map(String).slice(0, MAX_GENES) : [];
    st.overlays = Array.isArray(st.overlays) ? [...new Set(st.overlays.map(String))] : [];
    if (typeof st.dim !== 'boolean') st.dim = true;
    if (!st.hidden || typeof st.hidden !== 'object') st.hidden = {};
    if (!st.views || typeof st.views !== 'object') st.views = {};
    return st;
  }

  function markLayout(k) {
    el.dataset.layout = String(k);
    layoutBtns.forEach((b, j) => {
      const on = LAYOUTS[j] === k;
      b.classList.toggle('active', on);
      b.setAttribute('aria-pressed', on ? 'true' : 'false');
    });
    linkWrap.hidden = k === 1;
    scaleWrap.hidden = k === 1;
  }

  function setLayout(k) {
    if (!LAYOUTS.includes(k)) return;
    store.layout = k;
    markLayout(k);
    if (!samples || !samples.length) return;
    for (let i = 0; i < k; i += 1) ensurePanelState(i);
    while (panels.length > k) panels.pop().destroy();
    const added = [];
    while (panels.length < k) {
      const p = createPanel(panels.length);
      panels.push(p);
      added.push(p);
    }
    if (pointer && pointer.panel >= k) pointer = null;
    if (lastPanel >= k) lastPanel = 0;
    mount(panelsEl, panels.map((p) => p.el));
    for (const p of added) p.start();
    prune();
    updateSave();
    // a panel that went or came may change a shared color scale (and the legends say whose scale they show)
    for (const p of panels) if (!added.includes(p)) p.restyle();
    for (const p of panels) p.markWidth();
    for (const p of panels) p.requestRender();
    measureRoom();
  }

  function updateSave() {
    saveBtn.disabled = !panels.some((p) => p.hasData());
  }

  function showEmpty() {
    mount(el, h('div', { class: 'card' }, emptyState({
      icon: 'map', title: 'No per-bin data', text: emptyText || source.emptyText || 'There are no bins to map.',
    })));
  }

  async function start() {
    if (!samples) {
      mount(panelsEl, h('div', { class: 'card' }, loading('Loading sections…')));
      const my = new AbortController();
      startCtrl = my;
      try {
        const res = await source.samples(my.signal);
        if (destroyed) return;
        samples = normalizeSamples(res);
        if (samples.length) store.samples = samples;
      } catch (err) {
        if (destroyed || (err && err.name === 'AbortError')) return;
        if (err && err.status === 404) { showEmpty(); return; }
        mount(panelsEl, h('div', { class: 'card' }, h('div', { class: 'card-body' }, errorPanel(err, { retry: start }))));
        return;
      } finally {
        if (startCtrl === my) startCtrl = null;
      }
    }
    if (!samples.length) { showEmpty(); return; }
    setLayout(store.layout);
  }

  // ---- export
  function savePng() {
    if (!panels.some((p) => p.canExport())) return;
    const multi = panels.length > 1;
    const ratio = multi ? 2 : Math.min(3, Math.max(2, window.devicePixelRatio || 1));
    const legendW = 240;
    const titleH = 40;
    const gap = 24;
    const cols = multi ? 2 : 1;
    const rows = Math.ceil(panels.length / cols);
    const mapW = Math.max(200, ...panels.map((p) => p.size().w));
    const mapH = Math.max(200, ...panels.map((p) => Math.max(p.size().h, p.legendHeight())));
    const cellW = mapW + legendW;
    const cellH = titleH + mapH;
    const W = cols * cellW + (cols - 1) * gap;
    const H = rows * cellH + (rows - 1) * gap;
    const out = document.createElement('canvas');
    out.width = Math.round(W * ratio);
    out.height = Math.round(H * ratio);
    const g = out.getContext('2d');
    if (!g) { toast('This browser cannot create the image', 'error'); return; }
    g.fillStyle = EXPORT_PALETTE.bg;
    g.fillRect(0, 0, out.width, out.height);
    panels.forEach((p, i) => p.exportTo(g, (i % cols) * (cellW + gap), Math.floor(i / cols) * (cellH + gap), ratio, legendW, titleH));
    const stem = source.name || 'spatial';
    const name = `${stem}_${multi ? `${panels.length}_panels` : panels[0].fileTag()}.png`.replace(/[^\w.-]+/g, '_');
    out.toBlob((blob) => {
      if (blob) downloadBlob(name, blob);
      else toast('The image could not be created', 'error');
    }, 'image/png');
  }

  // ------------------------------------------------------------------------------------------------ one panel
  function createPanel(index) {
    const st = store.panels[index];
    const self = { index };
    let alive = true;
    let loadToken = 0;
    let data = null;
    let style = null;
    let scaleUsed = ''; // rangeText of the shared color range the style was built with ('' = its own range)
    let ovl = null;
    let vis = null; // 1 per bin drawn by the coloring or by an overlay (hover targets)
    let geneVecs = new Map(); // genes key -> prepared counts of the current section
    let geneErrors = new Map(); // genes key -> error message
    let genePending = new Set();
    let imageEntry = null;
    let imageError = null;
    let imagePending = false;
    let size = { w: 0, h: 0 };
    let dpr = window.devicePixelRatio || 1;
    let hover = -1;
    let dirty = true;
    let needPaint = true;
    let crispTimer = null;
    let layerState = null; // the view the layer was drawn with
    const layer = document.createElement('canvas');
    const pointers = new Map();
    let drag = null;
    let pinch = null;

    // ---- DOM
    const ids = { sample: uid('sv-sample'), field: uid('sv-field'), genes: uid('sv-genes'), overlay: uid('sv-overlay'), ovlGenes: uid('sv-ovl-genes') };
    const titleEl = h('h3', { class: 'sv-panel-title' });
    const sampleSel = h('select', { class: 'select', id: ids.sample, onChange: () => setSample(sampleSel.value) });
    const fieldSel = h('select', {
      class: 'select', id: ids.field,
      onChange: () => {
        setField(fieldSel.value);
        if (st.field === GENE_FIELD && !st.genes.length && geneBox) geneBox.focus();
      },
    });
    const searchGenes = (q, signal) => source.genes.search(st.sample, q, signal);
    const geneBox = hasGenes ? geneInput({ id: ids.genes, label: 'Genes', submitLabel: 'Show', search: searchGenes, onSubmit: setGenes }) : null;
    const geneCtl = geneBox ? h('div', { class: 'sv-control sv-control-wide', hidden: true },
      h('label', { for: ids.genes }, 'Genes (several are summed)'), geneBox.el) : null;
    const overlaySel = h('select', { class: 'select', id: ids.overlay, onChange: onOverlaySelect });
    const ovlGeneBox = hasGenes ? geneInput({
      id: ids.ovlGenes, label: 'Overlay genes', submitLabel: 'Add', search: searchGenes, onSubmit: addGeneOverlay, onCancel: closeGeneOverlay,
    }) : null;
    const ovlGeneCtl = ovlGeneBox ? h('div', { class: 'sv-control sv-control-wide', hidden: true },
      h('label', { for: ids.ovlGenes }, 'Overlay: bins detecting these genes'), ovlGeneBox.el) : null;
    const dimInput = h('input', { type: 'checkbox', checked: st.dim, onChange: () => { st.dim = dimInput.checked; requestRender(); } });
    const ovlChips = h('ul', { class: 'chip-list sv-ovl-chips', 'aria-label': 'Overlays' });
    const ovlRow = h('div', { class: 'sv-ovl-row', hidden: true }, ovlChips,
      h('label', { class: 'check sv-dim' }, dimInput, 'Dim the coloring under the overlays'));
    const controls = h('div', { class: 'sv-panel-controls' },
      h('div', { class: 'sv-control' }, h('label', { for: ids.sample }, 'Section'), sampleSel),
      h('div', { class: 'sv-control' }, h('label', { for: ids.field }, 'Color by'), fieldSel),
      geneCtl,
      h('div', { class: 'sv-control' }, h('label', { for: ids.overlay }, 'Overlays'), overlaySel),
      ovlGeneCtl);
    const canvas = h('canvas', {
      class: 'sv-canvas', tabIndex: 0, role: 'img',
      'aria-label': 'Map of the bins of the section. Drag or use the arrow keys to pan; scroll, + and − to zoom; 0 fits the section.',
    });
    const ctx = canvas.getContext('2d');
    // placed in the window (position: fixed) next to the pointer, so that the edges of a small panel do not cut it
    const tip = h('div', { class: 'sv-tip', hidden: true, role: 'tooltip' });
    const msg = h('div', { class: 'sv-msg', hidden: true });
    // compact icon buttons in one floating group on the map (run.css), named by aria-label and explained by title
    const zoomBtn = (iconName, title, onClick) => h('button', { class: 'btn small ghost icon-only', type: 'button', title, 'aria-label': title, onClick }, icon(iconName, { size: 14 }));
    const stage = h('div', { class: 'sv-stage' }, canvas,
      h('div', { class: 'sv-zoom', role: 'group', 'aria-label': 'Zoom' },
        zoomBtn('zoom-in', 'Zoom in (+)', () => zoomBy(1.6)),
        zoomBtn('zoom-out', 'Zoom out (−)', () => zoomBy(1 / 1.6)),
        zoomBtn('maximize', 'Fit (0); linked views fit every panel\'s section', fit)),
      tip, msg);
    const legendEl = h('div', { class: 'sv-legend-body' });
    const el = h('section', { class: 'card sv-panel', 'aria-label': `Panel ${index + 1}` },
      h('div', { class: 'sv-panel-head' }, titleEl), controls, ovlRow, h('div', { class: 'sv-panel-body' }, stage, legendEl));

    mount(sampleSel, samples.map((s) => h('option', { value: s.id },
      `${sampleText(s.id)}${Number.isFinite(s.n) ? ` — ${fmt.int(s.n)} bins` : ''}`)));
    sampleSel.disabled = samples.length < 2;
    sampleSel.value = st.sample;
    paintTitle();

    // ---- events
    canvas.addEventListener('pointerdown', onDown);
    canvas.addEventListener('pointermove', onMove);
    canvas.addEventListener('pointerup', onUp);
    canvas.addEventListener('pointercancel', onUp);
    canvas.addEventListener('pointerleave', () => {
      if (pointers.size) return;
      setHover(-1);
      if (pointer && pointer.panel === index) setPointer(null);
    });
    canvas.addEventListener('wheel', onWheel, { passive: false });
    canvas.addEventListener('dblclick', (e) => {
      if (!data) return;
      const p = local(e);
      zoomAt(e.shiftKey ? 1 / 2 : 2, p.x, p.y);
      viewMoved(self, true);
    });
    canvas.addEventListener('keydown', onKey);
    canvas.addEventListener('blur', () => setHover(-1));
    resizeObserver.observe(stage);
    resizeObserver.observe(el);

    /** A panel wide enough has its legend in a narrow column beside the map (two or four panels; see run.css). */
    function markWidth() {
      const w = el.clientWidth;
      if (w) el.classList.toggle('sv-wide', w >= (WIDE_PANEL[store.layout] || Infinity));
    }

    // ---- view geometry: a view is {cx, cy (µm at the centre of the canvas), scale (css px per µm)}
    function getView() {
      if (store.link) return goodView(store.view) ? store.view : fitLinked();
      const v = st.views[st.sample];
      return goodView(v) ? v : fitOwn();
    }

    function xf(v) {
      return { scale: v.scale, tx: size.w / 2 - v.cx * v.scale, ty: size.h / 2 - v.cy * v.scale };
    }

    function setView(v) {
      if (store.link) store.view = v;
      else st.views[st.sample] = v;
    }

    function fitOwn() {
      if (!data || !data.valid || !size.w || !size.h) return null;
      const b = data.bbox;
      const v = { cx: (b.x0 + b.x1) / 2, cy: (b.y0 + b.y1) / 2, scale: fitScaleFor(b, data.pitch, size), auto: true };
      st.views[st.sample] = v;
      return v;
    }

    const sidePx = (scale) => Math.max(MIN_SIDE, (data ? data.pitch : 16) * scale * store.pointSize);

    function local(e) {
      const r = canvas.getBoundingClientRect();
      return { x: e.clientX - r.left, y: e.clientY - r.top };
    }

    function clampScale(s) {
      let fitS = 1;
      if (store.link) {
        const f = linkedFit();
        if (f) fitS = f.scale;
      } else if (data && data.valid && size.w && size.h) {
        fitS = fitScaleFor(data.bbox, data.pitch, size);
      }
      return Math.min(Math.max(s, fitS / 4), Math.max(fitS, 240 / (data ? data.pitch : 16)));
    }

    function zoomAt(factor, px, py) {
      const v = data ? getView() : null;
      if (!v) return;
      const t = xf(v);
      const s1 = clampScale(v.scale * factor);
      const mx = (px - t.tx) / t.scale;
      const my = (py - t.ty) / t.scale;
      // the point under (px, py) stays there
      setView({ scale: s1, cx: mx + (size.w / 2 - px) / s1, cy: my + (size.h / 2 - py) / s1 });
    }

    function panBy(dx, dy) {
      const v = data ? getView() : null;
      if (!v) return;
      setView({ scale: v.scale, cx: v.cx - dx / v.scale, cy: v.cy - dy / v.scale });
    }

    function zoomBy(factor) {
      lastPanel = index;
      zoomAt(factor, size.w / 2, size.h / 2);
      viewMoved(self, true);
    }

    function fit() {
      if (!data || !data.valid) return;
      lastPanel = index;
      setHover(-1);
      if (store.link) fitLinked();
      else { fitOwn(); requestRender(); }
    }

    function resize() {
      // the drawable area inside the stage's border (the canvas fills it)
      const w = Math.max(0, Math.floor(stage.clientWidth));
      const hh = Math.max(0, Math.floor(stage.clientHeight));
      const ratio = window.devicePixelRatio || 1;
      if (!w || !hh) return; // hidden or detached: keep the current geometry
      if (w === size.w && hh === size.h && ratio === dpr) return;
      size = { w, h: hh };
      dpr = ratio;
      canvas.width = Math.max(1, Math.round(w * dpr));
      canvas.height = Math.max(1, Math.round(hh * dpr));
      if (!store.link) {
        const v = st.views[st.sample];
        if (v && v.auto) fitOwn();
      }
      requestRender();
    }

    // ---- rendering
    function requestRender() {
      dirty = true;
      needPaint = true;
      schedule();
    }

    function queuePaint() {
      needPaint = true;
      schedule();
    }

    function crisp() {
      clearTimeout(crispTimer);
      crispTimer = setTimeout(requestRender, 140);
    }

    function step() {
      if (!alive) return;
      if ((window.devicePixelRatio || 1) !== dpr) resize(); // moved to a screen with another pixel density
      if (!size.w || !size.h) return;
      if (dirty) renderLayer();
      if (needPaint) {
        needPaint = false;
        paint();
      }
    }

    const dimOf = () => (ovl && ovl.active && st.dim ? DIM : 0);

    function showImage() {
      return hasImage && store.image.on && imageEntry && imageEntry.img ? imageEntry : null;
    }

    /** Draw every visible bin into the offscreen layer for the current view (plus a margin around it). */
    function renderLayer() {
      dirty = false;
      clearTimeout(crispTimer);
      const v = data && style ? getView() : null;
      if (!v) { layerState = null; return; }
      const t = xf(v);
      let margin = MARGIN;
      while (margin > 0 && (size.w + 2 * margin) * (size.h + 2 * margin) * dpr * dpr > MAX_LAYER_PX) margin = Math.floor(margin / 2);
      const lw = Math.round((size.w + 2 * margin) * dpr);
      const lh = Math.round((size.h + 2 * margin) * dpr);
      if (layer.width !== lw || layer.height !== lh) { layer.width = lw; layer.height = lh; }
      const g = layer.getContext('2d');
      g.setTransform(1, 0, 0, 1, 0, 0);
      g.clearRect(0, 0, lw, lh);
      const base = theme || EXPORT_PALETTE;
      // over the H&E image (mostly white) the bins in several overlays take the dark color of the exported image
      const palette = showImage() ? { ...base, several: EXPORT_PALETTE.several } : base;
      const lt = { scale: t.scale, tx: t.tx + margin, ty: t.ty + margin };
      const side = sidePx(t.scale);
      drawBins(g, style.groups, data, lt, dpr, side, palette, dimOf());
      if (ovl) drawBins(g, ovl.groups, data, lt, dpr, side, palette, 0);
      layerState = { scale: t.scale, tx: t.tx, ty: t.ty, margin, dpr, w: size.w, h: size.h };
    }

    /**
     * Visible canvas = background + H&E image + the layer moved / scaled to the current view (see-through when the
     * image is shown) + hover marker + the crosshair of another panel + scale bar.
     */
    function paint() {
      const palette = theme || EXPORT_PALETTE;
      ctx.setTransform(1, 0, 0, 1, 0, 0);
      ctx.globalAlpha = 1;
      ctx.fillStyle = palette.bg;
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      const v = data ? getView() : null;
      if (!v) return;
      const t = xf(v);
      const im = showImage();
      if (im) {
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.imageSmoothingEnabled = true;
        drawImageUnder(ctx, im, t, 0, 0, size.w, size.h);
        ctx.setTransform(1, 0, 0, 1, 0, 0);
      }
      if (layerState && layerState.dpr === dpr) {
        const k = t.scale / layerState.scale;
        const left = t.tx - (layerState.tx + layerState.margin) * k; // css px
        const top = t.ty - (layerState.ty + layerState.margin) * k;
        const lw = (layerState.w + 2 * layerState.margin) * k;
        const lh = (layerState.h + 2 * layerState.margin) * k;
        ctx.imageSmoothingEnabled = k !== 1;
        ctx.globalAlpha = im ? Math.max(0, 1 - store.image.opacity) : 1;
        ctx.drawImage(layer, left * dpr, top * dpr, lw * dpr, lh * dpr);
        ctx.globalAlpha = 1;
        // the layer no longer covers the viewport (long pan, zooming out): redraw on the next frame
        if (left > 0.5 || top > 0.5 || left + lw < size.w - 0.5 || top + lh < size.h - 0.5) requestRender();
      } else if (style) {
        requestRender();
      }
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      if (hover >= 0 && hover < data.n) {
        const half = Math.max(sidePx(t.scale) / 2, 3) + 2;
        const cx = data.x[hover] * t.scale + t.tx;
        const cy = data.y[hover] * t.scale + t.ty;
        ctx.lineWidth = 3;
        ctx.strokeStyle = palette.bg;
        ctx.strokeRect(cx - half, cy - half, half * 2, half * 2);
        ctx.lineWidth = 1.5;
        ctx.strokeStyle = palette.text;
        ctx.strokeRect(cx - half, cy - half, half * 2, half * 2);
      }
      if (store.link && pointer && pointer.panel !== index) {
        const cx = pointer.x * t.scale + t.tx;
        const cy = pointer.y * t.scale + t.ty;
        if (cx >= 0 && cx <= size.w && cy >= 0 && cy <= size.h) drawCrosshair(ctx, cx, cy, size.w, size.h, palette);
      }
      drawScaleBar(ctx, t.scale, 10, size.h - 10, palette);
    }

    // ---- pointer, wheel and keyboard
    function onDown(e) {
      if (!data) return;
      if (e.pointerType === 'mouse' && e.button !== 0) return;
      lastPanel = index;
      try { canvas.setPointerCapture(e.pointerId); } catch { /* pointer already released */ }
      const p = local(e);
      pointers.set(e.pointerId, p);
      if (pointers.size === 1) {
        drag = { x: p.x, y: p.y, moved: false };
      } else if (pointers.size === 2) {
        drag = null;
        pinch = pinchState();
      }
    }

    function pinchState() {
      const [a, b] = [...pointers.values()];
      return { d: Math.max(1, Math.hypot(a.x - b.x, a.y - b.y)), x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 };
    }

    function onMove(e) {
      const p = local(e);
      if (pointers.has(e.pointerId)) {
        pointers.set(e.pointerId, p);
        if (pinch && pointers.size === 2) {
          const now = pinchState();
          zoomAt(now.d / pinch.d, now.x, now.y);
          panBy(now.x - pinch.x, now.y - pinch.y);
          pinch = now;
          viewMoved(self, true);
        } else if (drag) {
          const dx = p.x - drag.x;
          const dy = p.y - drag.y;
          if (!drag.moved && Math.abs(dx) + Math.abs(dy) < 3) return;
          if (!drag.moved) {
            drag.moved = true;
            stage.classList.add('dragging');
            setHover(-1);
            setPointer(null);
          }
          panBy(dx, dy);
          drag.x = p.x;
          drag.y = p.y;
          viewMoved(self, false);
        }
        return;
      }
      if (e.pointerType !== 'touch') hoverAt(p.x, p.y);
    }

    function onUp(e) {
      if (!pointers.has(e.pointerId)) return;
      pointers.delete(e.pointerId);
      if (pointers.size < 2) pinch = null;
      if (pointers.size === 1) {
        const [p] = pointers.values();
        drag = { x: p.x, y: p.y, moved: true };
        return;
      }
      if (pointers.size) return;
      const moved = drag ? drag.moved : true;
      drag = null;
      stage.classList.remove('dragging');
      if (moved) {
        viewSettled(self);
      } else if (e.type === 'pointerup') {
        const p = local(e); // a click or tap shows the bin under it
        hoverAt(p.x, p.y);
      }
    }

    function onWheel(e) {
      if (!data) return;
      e.preventDefault();
      lastPanel = index;
      const unit = e.deltaMode === 1 ? 16 : e.deltaMode === 2 ? 400 : 1;
      const p = local(e);
      // trackpad pinches arrive as ctrl + wheel with small deltas
      zoomAt(Math.exp(-e.deltaY * unit * (e.ctrlKey ? 0.01 : 0.0015)), p.x, p.y);
      viewMoved(self, true);
      hoverAt(p.x, p.y);
    }

    function onKey(e) {
      if (!data || e.altKey || e.ctrlKey || e.metaKey) return;
      const stepPx = 60;
      lastPanel = index;
      switch (e.key) {
        case '+': case '=': zoomBy(1.4); break;
        case '-': case '_': zoomBy(1 / 1.4); break;
        case '0': fit(); break;
        case 'ArrowLeft': panBy(stepPx, 0); viewMoved(self, true); break;
        case 'ArrowRight': panBy(-stepPx, 0); viewMoved(self, true); break;
        case 'ArrowUp': panBy(0, stepPx); viewMoved(self, true); break;
        case 'ArrowDown': panBy(0, -stepPx); viewMoved(self, true); break;
        case 'Escape': setHover(-1); setPointer(null); return;
        default: return;
      }
      e.preventDefault();
    }

    // ---- hover
    function hoverAt(px, py) {
      const v = data ? getView() : null;
      if (!v) { setHover(-1); return; }
      const t = xf(v);
      const mx = (px - t.tx) / t.scale;
      const my = (py - t.ty) / t.scale;
      setPointer({ panel: index, x: mx, y: my });
      if (!style || !vis || !data.grid) { setHover(-1); return; }
      const r = Math.max(HOVER_PX, sidePx(t.scale) * 0.71 + 1) / t.scale;
      setHover(nearest(data.grid, data, mx, my, r, vis), px, py);
    }

    function valueText(f, i) {
      if (f.type === 'category') return f.codes[i] < 0 ? '—' : f.categories[f.codes[i]];
      if (f.type === 'bool') return f.values[i] === 2 ? '—' : (f.values[i] ? 'yes' : 'no');
      return Number.isFinite(f.values[i]) ? formatNumber(f.values[i]) : '—';
    }

    /**
     * Tooltip: position, barcode, the counts of the genes shown, the field colored by and the overlays (in bold), the
     * other labels and numbers, and the positivity calls that are not shown, on one line (the ones the bin has).
     */
    function tipContent(i) {
      const rows = [];
      const row = (label, value, on, cls) => {
        rows.push(h('dt', { class: on ? 'on' : null }, label), h('dd', { class: [on ? 'on' : null, cls] }, value));
      };
      if (data.barcodes) row('Barcode', String(data.barcodes[i] ?? '—'), false, 'mono');
      const colorKey = st.field === GENE_FIELD && st.genes.length ? genesKey(st.genes) : null;
      for (const gk of neededGeneKeys()) {
        const gv = geneVecs.get(gk);
        if (gv) row(`${geneText(gv.used.length ? gv.used : gk.split(','))} (UMIs)`, formatNumber(gv.counts[i]), gk === colorKey);
      }
      const calls = [];
      let anyCall = false;
      for (const key of data.order) {
        const f = data.fields.get(key);
        const on = key === st.field || st.overlays.includes(key);
        if (!on && f.type === 'bool' && key.startsWith('positive_')) {
          anyCall = true;
          if (f.values[i] === 1) calls.push(fieldLabel(key).replace(/^Positive: /, ''));
        } else {
          row(fieldLabel(key), valueText(f, i), on);
        }
      }
      if (anyCall) row('Positive for', calls.length ? calls.join(', ') : 'none of the other calls', false, 'sv-tip-calls');
      return [
        h('div', { class: 'sv-tip-head' }, `x ${fmt.int(data.x[i])} µm · y ${fmt.int(data.y[i])} µm`),
        rows.length ? h('dl', { class: 'sv-tip-grid' }, rows) : null,
      ];
    }

    function hideTip() {
      if (hover >= 0) setHover(-1);
    }

    /**
     * Tooltip next to the pointer (px, py: css px in the canvas), in window coordinates: beside it when it fits on the
     * right or the left, else above or below it, never over it when the window has room.
     */
    function placeTip(px, py) {
      const r = canvas.getBoundingClientRect();
      const ax = r.left + px;
      const ay = r.top + py;
      const vw = document.documentElement.clientWidth;
      const vh = document.documentElement.clientHeight;
      const pad = 14;
      const edge = 8;
      const tw = tip.offsetWidth;
      const th = tip.offsetHeight;
      const clampTo = (v, lo, hi) => Math.max(lo, Math.min(v, hi));
      let left;
      let top;
      if (ax + pad + tw <= vw - edge || ax - pad - tw >= edge) {
        left = ax + pad + tw <= vw - edge ? ax + pad : ax - pad - tw;
        top = clampTo(ay + pad, edge, vh - th - edge);
      } else {
        left = clampTo(ax - tw / 2, edge, vw - tw - edge);
        if (ay + pad + th <= vh - edge) top = ay + pad;
        else if (ay - pad - th >= edge) top = ay - pad - th;
        else top = ay > vh / 2 ? edge : vh - th - edge;
      }
      tip.style.transform = `translate(${Math.round(left)}px, ${Math.round(top)}px)`;
    }

    function setHover(i, px, py) {
      if (i !== hover) {
        hover = i;
        queuePaint();
        if (i < 0) {
          tip.hidden = true;
          return;
        }
        mount(tip, tipContent(i));
        tip.hidden = false;
      }
      if (i >= 0 && px !== undefined) placeTip(px, py);
    }

    // ---- title and legend
    function geneNames() {
      const gv = st.genes.length ? geneVecs.get(genesKey(st.genes)) : null;
      return gv && gv.used.length ? gv.used : st.genes;
    }

    function shownText() {
      let what;
      if (st.field === GENE_FIELD) what = st.genes.length ? `${geneText(geneNames())} (UMIs)` : 'gene expression';
      else what = st.field && data.fields.has(st.field) ? fieldLabel(st.field) : 'bins';
      const n = st.overlays.length;
      if (n === 1 || n === 2) what += ` + ${st.overlays.map(overlayLabel).join(', ')}`;
      else if (n > 2) what += ` + ${n} overlays`;
      return what;
    }

    function titleText() {
      return data ? `${sampleText(st.sample)} — ${shownText()}` : sampleText(st.sample);
    }

    function paintTitle() {
      const text = titleText();
      titleEl.textContent = text;
      titleEl.title = text;
    }

    function legendTitle() {
      if (st.field === GENE_FIELD) return st.genes.length ? `${geneText(geneNames())} (UMIs)` : 'Gene expression';
      return st.field && data.fields.has(st.field) ? fieldLabel(st.field) : 'Bins';
    }

    /** Color bar of a scale; whose range it is when there are several panels (shared or this panel's own). */
    function colorbar(scale) {
      const { min, max, log, shared } = scale;
      const bar = h('canvas', { class: 'sv-colorbar', width: 256, height: 1, 'aria-hidden': 'true' });
      const g = bar.getContext('2d');
      for (let k = 0; k < 256; k += 1) {
        g.fillStyle = LUT[Math.min(LEVELS - 1, Math.floor((k / 256) * LEVELS))];
        g.fillRect(k, 0, 1, 1);
      }
      // the value at the middle of the bar; log scales are counts (UMIs, genes): a whole number, left out when it
      // would repeat an end
      let mid = log ? Math.round(Math.expm1((Math.log1p(min) + Math.log1p(max)) / 2))
        : Math.round(((min + max) / 2) * 100) / 100;
      if (!(mid > min && mid < max)) mid = null;
      const whose = panels.length < 2 ? '' : (shared ? '; the same in the panels showing it' : '; this panel’s own');
      return [bar,
        h('div', { class: 'sv-colorbar-labels' },
          h('span', null, formatNumber(min)), h('span', null, mid === null ? '' : formatNumber(mid)),
          h('span', null, `${clampedMark(scale)}${formatNumber(max)}`)),
        h('div', { class: 'sv-note' }, `${log ? 'Logarithmic color scale (log 1 + x)' : 'Linear color scale'}${whose}`)];
    }

    /** "Colors stop at x, ..." when the scale is cut below the largest value. */
    function clampFact(scale, whose) {
      if (!scale || !scale.clamped) return null;
      return h('li', null, scale.shared
        ? `Colors stop at ${formatNumber(scale.max)}, the largest 99th percentile of the panels showing it`
        : `Colors stop at ${formatNumber(scale.max)}, the 99th percentile of ${whose}`);
    }

    function retryLine(text, retry) {
      return h('div', { class: 'sv-error' }, icon('alert', { size: 14 }), h('span', null, text),
        h('button', { class: 'btn ghost small', type: 'button', onClick: retry }, 'Retry'));
    }

    function busyLine(text) {
      return h('div', { class: 'sv-busy' }, h('span', { class: 'spinner' }), text);
    }

    function retryGenes(gk) {
      geneErrors.delete(gk);
      ensureGenes();
      restyle();
    }

    function geneLegend() {
      if (!st.genes.length) return h('p', { class: 'sv-hint' }, 'Type one or more genes above; the counts of several genes are summed per bin.');
      const gk = genesKey(st.genes);
      const err = geneErrors.get(gk);
      if (err) return retryLine(err, () => retryGenes(gk));
      const gv = geneVecs.get(gk);
      if (!gv) return busyLine('Loading the counts…');
      const facts = [
        h('li', null, `Detected in ${fmt.int(gv.detected)} of ${fmt.int(data.n)} bins (${fmt.pct(data.n ? (100 * gv.detected) / data.n : 0, 1)})`),
      ];
      if (gv.detected) facts.push(h('li', null, `Up to ${formatNumber(gv.max)} UMIs in a bin, ${formatNumber(gv.total)} in all`));
      const scale = style && style.kind === 'gene' ? style.scale : null;
      const cut = clampFact(scale, 'the bins detecting it');
      if (cut) facts.push(cut);
      if (gv.missing.length) facts.push(h('li', { class: 'sv-warn' }, `Not in the data: ${gv.missing.join(', ')}`));
      if (!gv.used.length) facts.push(h('li', { class: 'sv-warn' }, 'None of these genes is in the data.'));
      if (gv.unmatched) facts.push(h('li', { class: 'sv-warn' }, `${fmt.int(gv.unmatched)} bins of the map were not found in the raw data (counted as 0)`));
      return [scale ? colorbar(scale) : null, h('ul', { class: 'sv-facts' }, facts)];
    }

    function numberLegend() {
      const cut = clampFact(style.scale, 'the bins');
      return [colorbar(style.scale), cut ? h('ul', { class: 'sv-facts' }, cut) : null];
    }

    /** Color of the bins in two or more overlays as drawn now (see SEVERAL; its swatch has a visible rim when dark). */
    const severalColor = () => (showImage() ? EXPORT_PALETTE.several : 'var(--text)');

    function overlayLegend() {
      const layers = ovl ? ovl.layers : [];
      const items = layers.map((l) => h('li', { class: 'sv-ovl' },
        h('span', { class: 'swatch', style: { background: l.color } }),
        h('span', { class: 'sv-item-label' }, l.label),
        h('span', { class: 'sv-item-count' }, l.count === null ? l.status : fmt.int(l.count))));
      if (ovl && ovl.several) {
        items.push(h('li', { class: 'sv-ovl' },
          h('span', { class: 'swatch sv-swatch-ink', style: { background: severalColor() } }),
          h('span', { class: 'sv-item-label' }, 'In two or more overlays'),
          h('span', { class: 'sv-item-count' }, fmt.int(ovl.several))));
      }
      const errors = st.overlays.filter((k) => k.startsWith(GENE_OVERLAY)).map((k) => k.slice(GENE_OVERLAY.length))
        .filter((gk) => geneErrors.has(gk) && !(st.field === GENE_FIELD && genesKey(st.genes) === gk))
        .map((gk) => retryLine(`${geneText(gk.split(','))}: ${geneErrors.get(gk)}`, () => retryGenes(gk)));
      return h('div', { class: 'sv-legend-section' },
        h('h4', null, 'Overlays'),
        h('ul', { class: 'sv-ovl-list', title: 'Bins in each overlay (overlap included)' }, items),
        errors);
    }

    function imageLegend() {
      if (!hasImage || !store.image.on) return null;
      if (imageError) {
        return retryLine(`H&E image: ${imageError}`, () => {
          imageError = null;
          ensureImage(true);
          paintLegend();
        });
      }
      if (imagePending) return busyLine('Loading the H&E image…');
      if (imageEntry && imageEntry.none) return h('p', { class: 'sv-hint' }, 'No H&E image was found for this section.');
      return null;
    }

    function paintLegend() {
      if (!data) { mount(legendEl); return; }
      const parts = [h('h3', { class: 'sv-legend-title' }, legendTitle())];
      if (st.field === GENE_FIELD) parts.push(geneLegend());
      else if (style && style.kind === 'number' && style.scale) parts.push(numberLegend());
      if (style && style.legend.length) {
        parts.push(h('div', { class: 'sv-legend', role: 'group', 'aria-label': `${legendTitle()}: click an entry to hide or show it` },
          style.legend.map((e) => h('button', {
            class: e.off ? 'sv-item off' : 'sv-item', type: 'button', 'aria-pressed': e.off ? 'false' : 'true',
            dataset: { key: e.key }, title: 'Click to hide or show; Alt+click to show only this entry',
            onClick: (ev) => toggle(e.key, ev.altKey),
          },
          h('span', { class: 'swatch', style: { background: e.color === NEUTRAL ? 'var(--plot-neutral)' : e.color } }),
          h('span', { class: 'sv-item-label' }, e.label),
          h('span', { class: 'sv-item-count' }, fmt.int(e.count))))));
      }
      if (style) {
        const anyOff = style.legend.some((e) => e.off);
        parts.push(h('div', { class: 'sv-legend-foot' },
          h('span', null, `${fmt.int(style.shown)} of ${fmt.int(data.valid)} bins shown`),
          anyOff ? h('button', { class: 'btn ghost small', type: 'button', onClick: showAll }, 'Show all') : null));
      }
      if (st.overlays.length) parts.push(overlayLegend());
      parts.push(imageLegend());
      mount(legendEl, parts);
    }

    const fieldKey = () => (st.field === GENE_FIELD ? GENE_FIELD : st.field);

    function toggle(key, solo) {
      if (!style) return;
      const fk = fieldKey();
      let hidden = new Set(st.hidden[fk] || []);
      if (solo) {
        const keys = style.legend.map((e) => e.key);
        const alone = keys.every((k) => (k === key ? !hidden.has(k) : hidden.has(k)));
        hidden = alone ? new Set() : new Set(keys.filter((k) => k !== key));
      } else if (hidden.has(key)) {
        hidden.delete(key);
      } else {
        hidden.add(key);
      }
      st.hidden[fk] = [...hidden];
      restyle();
      const again = [...legendEl.querySelectorAll('.sv-item')].find((b) => b.dataset.key === key);
      if (again) again.focus();
    }

    function showAll() {
      st.hidden[fieldKey()] = [];
      restyle();
      canvas.focus();
    }

    function paintOverlayRow() {
      ovlRow.hidden = !st.overlays.length;
      const layers = ovl ? ovl.layers : [];
      mount(ovlChips, st.overlays.map((key) => {
        const l = layers.find((x) => x.key === key);
        const label = overlayLabel(key);
        return h('li', { class: 'chip sv-chip', title: l && l.count === null ? `${label}: ${l.status}` : label },
          h('span', { class: 'swatch', style: { background: l ? l.color : 'transparent' } }),
          h('span', { class: 'chip-text' }, label),
          h('button', { type: 'button', 'aria-label': `Remove the overlay ${label}`, title: 'Remove', onClick: () => removeOverlay(key) },
            icon('x', { size: 12 })));
      }));
    }

    // ---- coloring and overlays
    /** Overlay layers of this panel with their colors: the configured signature color when free, else the palette. */
    function overlayLayers() {
      const layers = st.overlays.map((key) => {
        if (key.startsWith(GENE_OVERLAY)) {
          const gk = key.slice(GENE_OVERLAY.length);
          const gv = geneVecs.get(gk);
          return {
            key, label: overlayLabel(key), color: null, mask: gv ? (i) => gv.counts[i] > 0 : null,
            status: gv ? '' : (geneErrors.has(gk) ? 'failed' : 'loading…'),
          };
        }
        const f = data.fields.get(key);
        const isBool = !!f && f.type === 'bool';
        return {
          key, label: overlayLabel(key), color: BOOL_COLORS[key] || signatureColor(key, data.colors),
          mask: isBool ? (i) => f.values[i] === 1 : null, status: isBool ? '' : 'not in this section',
        };
      });
      const used = [];
      for (const l of layers) {
        if (l.color && !tooClose(l.color, used)) used.push(l.color);
        else l.color = null; // e.g. the pooled and the native call of one signature: the second takes a palette color
      }
      // palette colors apart from the other overlays, first the ones also apart from the colors of the coloring under
      // them (categories, calls); configured colors stay as they are (one signature, one color)
      const under = style ? style.legend.map((e) => e.color).filter((c) => c !== NEUTRAL) : [];
      const free = OVERLAY_PALETTE.filter((c) => !tooClose(c, used));
      const apart = free.filter((c) => !tooClose(c, under));
      const order = [...apart, ...free.filter((c) => !apart.includes(c))];
      let next = 0;
      for (const l of layers) {
        if (l.color) continue;
        l.color = order.length ? order[next % order.length] : OVERLAY_PALETTE[next % OVERLAY_PALETTE.length];
        next += 1;
      }
      return layers;
    }

    /**
     * What the color scale of this panel's coloring is fitted to, on its own: {key (genes summed or number field),
     * lo, hi} or null (labels, calls, nothing loaded yet).
     */
    function ownScale() {
      if (!data) return null;
      if (st.field === GENE_FIELD) {
        const gv = st.genes.length ? geneVecs.get(genesKey(st.genes)) : null;
        if (!gv || !gv.detected) return null;
        const names = (gv.used.length ? gv.used : st.genes).map((g) => g.toLowerCase()).sort();
        return { key: `gene\n${genesKey(names)}`, lo: gv.lo, hi: gv.clamp };
      }
      const f = data.fields.get(st.field);
      const r = f && f.type === 'number' ? numberRange(f) : null;
      return r ? { key: `field\n${f.key}`, lo: r.lo, hi: r.hi } : null;
    }

    function restyleIfScaleMoved() {
      if (data && rangeText(scaleRangeFor(ownScale())) !== scaleUsed) restyle();
    }

    function restyle() {
      if (!data) {
        style = null;
        ovl = null;
        vis = null;
        scaleUsed = '';
        paintTitle();
        paintLegend();
        requestRender();
        syncScales(self);
        return;
      }
      const range = scaleRangeFor(ownScale());
      scaleUsed = rangeText(range);
      if (st.field === GENE_FIELD) {
        const gv = st.genes.length ? geneVecs.get(genesKey(st.genes)) : null;
        style = gv ? buildGeneStyle(data, gv, new Set(st.hidden[GENE_FIELD] || []), range) : buildStyle(data, null, new Set());
      } else {
        const f = data.fields.get(st.field) || null;
        style = buildStyle(data, f, new Set(f ? st.hidden[f.key] || [] : []), range);
      }
      ovl = st.overlays.length ? buildOverlays(data, overlayLayers()) : null;
      if (ovl && ovl.groups.length) {
        vis = style.visible.slice();
        for (let i = 0; i < data.n; i += 1) if (ovl.member[i] >= 0) vis[i] = 1;
      } else {
        vis = style.visible;
      }
      if (hover >= 0 && !vis[hover]) setHover(-1);
      paintTitle();
      paintLegend();
      paintOverlayRow();
      requestRender();
      syncScales(self);
      measureSoon(); // the overlay row may have come or gone
    }

    function normalizeField() {
      if (st.field === GENE_FIELD && hasGenes) return;
      if (!data.fields.has(st.field)) st.field = DEFAULT_FIELDS.find((k) => data.fields.has(k)) || data.order[0] || '';
    }

    function fillFieldSelect() {
      const keysOf = (type) => data.order.filter((k) => data.fields.get(k).type === type);
      const groups = [
        { label: 'Labels', keys: keysOf('category') },
        { label: 'Region, quality control and calls', keys: keysOf('bool') },
        { label: 'Numbers', keys: keysOf('number') },
      ].filter((g) => g.keys.length);
      const opts = groups.map((g) => h('optgroup', { label: g.label }, g.keys.map((k) => h('option', { value: k }, fieldLabel(k)))));
      if (hasGenes) opts.push(h('optgroup', { label: 'Expression' }, h('option', { value: GENE_FIELD }, 'Gene…')));
      if (!opts.length) opts.push(h('option', { value: '' }, 'Position only'));
      mount(fieldSel, opts);
      fieldSel.disabled = !groups.length && !hasGenes;
    }

    function fillOverlaySelect() {
      const opts = [h('option', { value: '' }, st.overlays.length ? 'Add another…' : 'Add an overlay…')];
      const bools = data ? data.order.filter((k) => data.fields.get(k).type === 'bool' && !st.overlays.includes(k)) : [];
      if (bools.length) {
        opts.push(h('optgroup', { label: 'Region, quality control and calls' }, bools.map((k) => h('option', { value: k }, fieldLabel(k)))));
      }
      if (hasGenes) {
        const geneOpts = [];
        const current = st.field === GENE_FIELD && st.genes.length ? GENE_OVERLAY + genesKey(st.genes) : null;
        if (current && !st.overlays.includes(current)) geneOpts.push(h('option', { value: current }, overlayLabel(current)));
        geneOpts.push(h('option', { value: GENE_OVERLAY }, 'Gene detected…'));
        opts.push(h('optgroup', { label: 'Expression' }, geneOpts));
      }
      mount(overlaySel, opts);
      overlaySel.value = '';
      overlaySel.disabled = opts.length < 2;
    }

    function syncControls() {
      fillFieldSelect();
      fieldSel.value = st.field;
      if (geneCtl) geneCtl.hidden = st.field !== GENE_FIELD;
      if (geneBox) geneBox.set(st.genes.join(', '));
      fillOverlaySelect();
      dimInput.checked = st.dim;
    }

    function setField(key) {
      st.field = key; // chosen while a section loads: kept, and checked against its fields once it is there
      if (!data) {
        if (geneCtl) geneCtl.hidden = st.field !== GENE_FIELD; // genes can be typed meanwhile
        return;
      }
      normalizeField();
      syncControls();
      ensureGenes();
      prune();
      restyle();
    }

    function setGenes(genes) {
      let list = genes;
      if (list.length > MAX_GENES) {
        toast(`At most ${MAX_GENES} genes are summed; the first ${MAX_GENES} are shown.`, 'warning');
        list = list.slice(0, MAX_GENES);
      }
      st.genes = list;
      if (hasGenes) st.field = GENE_FIELD;
      geneErrors.delete(genesKey(list)); // a gene list that failed before is tried again
      if (!data) return;
      syncControls();
      ensureGenes();
      prune();
      restyle();
    }

    function onOverlaySelect() {
      const v = overlaySel.value;
      overlaySel.value = '';
      if (!v) return;
      if (v === GENE_OVERLAY) {
        if (ovlGeneCtl) {
          ovlGeneCtl.hidden = false;
          ovlGeneBox.focus();
        }
        return;
      }
      addOverlay(v);
    }

    function addOverlay(key) {
      if (!st.overlays.includes(key)) st.overlays.push(key);
      if (key.startsWith(GENE_OVERLAY)) geneErrors.delete(key.slice(GENE_OVERLAY.length));
      if (!data) { paintOverlayRow(); return; } // a section is loading: the chips follow now, the map once it is there
      fillOverlaySelect();
      ensureGenes();
      prune();
      restyle();
    }

    function addGeneOverlay(genes) {
      if (!genes.length) return;
      if (genes.length > MAX_GENES) toast(`At most ${MAX_GENES} genes are summed; the first ${MAX_GENES} are used.`, 'warning');
      ovlGeneBox.set('');
      ovlGeneCtl.hidden = true;
      addOverlay(GENE_OVERLAY + genesKey(genes.slice(0, MAX_GENES)));
      overlaySel.focus();
    }

    function closeGeneOverlay() {
      ovlGeneBox.set('');
      ovlGeneCtl.hidden = true;
      overlaySel.focus();
    }

    function removeOverlay(key) {
      st.overlays = st.overlays.filter((k) => k !== key);
      if (data) {
        fillOverlaySelect();
        prune();
        restyle();
      } else {
        paintOverlayRow();
      }
      overlaySel.focus();
    }

    // ---- data
    function neededGeneKeys() {
      const keys = [];
      if (st.field === GENE_FIELD && st.genes.length) keys.push(genesKey(st.genes));
      for (const k of st.overlays) if (k.startsWith(GENE_OVERLAY)) keys.push(k.slice(GENE_OVERLAY.length));
      return [...new Set(keys)];
    }

    /** Request keys this panel waits for (requests nobody waits for are aborted). */
    function wants() {
      const keys = [`bins\n${st.sample}`];
      for (const gk of neededGeneKeys()) keys.push(`gene\n${st.sample}\n${gk}`);
      if (hasImage && store.image.on) keys.push(`image\n${st.sample}`);
      return keys;
    }

    /** Load the counts of the genes shown (color and overlays) that are not loaded yet. */
    function ensureGenes() {
      if (!data || !hasGenes) return;
      const sid = st.sample;
      const d = data;
      const token = loadToken;
      for (const gk of neededGeneKeys()) {
        if (geneVecs.has(gk) || geneErrors.has(gk) || genePending.has(gk)) continue;
        const cached = cachedGene(`${sid}\n${gk}`);
        if (cached) { geneVecs.set(gk, cached); continue; }
        genePending.add(gk);
        loadGene(sid, d, gk.split(',')).then((gv) => {
          if (!alive || token !== loadToken) return;
          genePending.delete(gk);
          geneVecs.set(gk, gv);
          restyle();
        }, (err) => {
          if (!alive || token !== loadToken) return;
          genePending.delete(gk);
          if (err && err.name === 'AbortError') return;
          geneErrors.set(gk, (err && err.message) || String(err));
          restyle();
        });
      }
    }

    /** Load the H&E image of the section when it is shown; retry: after the Retry of a failed load. */
    function ensureImage(retry = false) {
      if (!hasImage || !store.image.on || !data) return;
      const sid = st.sample;
      const cached = store.images.get(sid);
      if (cached) { imageEntry = cached; return; }
      if (imagePending || imageError) return;
      imagePending = true;
      const token = loadToken;
      loadImage(sid, retry).then((e) => {
        if (!alive || token !== loadToken) return;
        imagePending = false;
        imageEntry = e;
        paintLegend();
        requestRender(); // the bins in several overlays change color over the image
      }, (err) => {
        if (!alive || token !== loadToken) return;
        imagePending = false;
        if (err && err.name === 'AbortError') return;
        imageError = (err && err.message) || String(err);
        paintLegend();
      });
    }

    function imageChanged() {
      ensureImage();
      paintLegend();
      requestRender();
    }

    function showMsg(content) {
      mount(msg, h('div', { class: 'sv-msg-body' }, content));
      msg.hidden = false;
    }

    function hideMsg() {
      mount(msg);
      msg.hidden = true;
    }

    function setSample(sid) {
      st.sample = sid;
      if (sampleSel.value !== sid) sampleSel.value = sid;
      return load();
    }

    async function load() {
      const token = ++loadToken;
      const sid = st.sample;
      setHover(-1);
      geneVecs = new Map();
      geneErrors = new Map();
      genePending = new Set();
      imageEntry = null;
      imageError = null;
      imagePending = false;
      prune();
      let d = store.data.get(sid);
      if (!d) {
        data = null;
        layerState = null;
        restyle();
        updateSave();
        const meta = samples.find((s) => s.id === sid);
        showMsg(loading(meta && Number.isFinite(meta.n) ? `Loading ${fmt.int(meta.n)} bins…` : 'Loading bins…'));
        try {
          d = await loadBins(sid);
        } catch (err) {
          if (!alive || token !== loadToken || (err && err.name === 'AbortError')) return;
          showMsg(errorPanel(err, { retry: load }));
          return;
        }
        if (!alive || token !== loadToken) return;
      }
      data = d;
      layerState = null;
      normalizeField();
      syncControls();
      if (d.valid) hideMsg();
      else showMsg(emptyState({ icon: 'map', title: 'No bins to draw', text: 'This section has no bins with coordinates.' }));
      // the linked view fits the sections shown until someone moves it
      if (store.link && (!goodView(store.view) || store.view.auto)) fitLinked();
      ensureGenes();
      ensureImage();
      prune();
      restyle();
      updateSave();
    }

    /** viewer.show(): change the section, the coloring or the overlays of this panel. */
    function apply(patch) {
      const before = st.sample;
      patchState(st, patch);
      if (!samples.some((s) => s.id === st.sample)) st.sample = before;
      if (Array.isArray(patch.genes) || typeof patch.genes === 'string') geneErrors.delete(genesKey(st.genes));
      if (st.sample !== before || !data) {
        sampleSel.value = st.sample;
        load();
        return;
      }
      normalizeField();
      syncControls();
      ensureGenes();
      prune();
      restyle();
    }

    // ---- export
    function fileTag() {
      const what = st.field === GENE_FIELD ? genesKey(st.genes) || 'genes' : st.field || 'bins';
      return `${st.sample}_${what}`;
    }

    function legendHeight() {
      if (!data) return 0;
      let hgt = 8 + 24;
      if (style && style.scale) hgt += 36 + 18 + (style.scale.clamped ? 18 : 0);
      if (style && style.kind === 'gene') hgt += 18 * (1 + (style.gene.missing.length ? 1 : 0));
      hgt += 19 * (style ? style.legend.filter((e) => !e.off).length : 0);
      if (ovl && ovl.layers.length) hgt += 26 + 19 * (ovl.layers.length + (ovl.several ? 1 : 0));
      return hgt + 8;
    }

    function drawLegendTo(g, x, y, width, P) {
      g.save();
      g.textBaseline = 'top';
      let top = y;
      const swatchRow = (color, label, count) => {
        g.fillStyle = groupColor(color, P);
        g.fillRect(x, top + 1, 12, 12);
        g.font = `12px ${FONT}`;
        const text = count === null || count === undefined ? '—' : fmt.int(count);
        const cw = g.measureText(text).width;
        g.fillStyle = P.text;
        g.fillText(label, x + 18, top, Math.max(10, width - 26 - cw));
        g.fillStyle = P.muted;
        g.fillText(text, x + width - cw, top);
        top += 19;
      };
      g.fillStyle = P.text;
      g.font = `600 13px ${FONT}`;
      g.fillText(legendTitle(), x, top, width);
      top += 24;
      g.font = `12px ${FONT}`;
      if (style && style.scale) {
        for (let k = 0; k < width; k += 1) {
          g.fillStyle = LUT[Math.min(LEVELS - 1, Math.floor((k / width) * LEVELS))];
          g.fillRect(x + k, top, 1, 12);
        }
        g.fillStyle = P.muted;
        g.fillText(formatNumber(style.scale.min), x, top + 16);
        const maxText = `${clampedMark(style.scale)}${formatNumber(style.scale.max)}`;
        g.fillText(maxText, x + width - g.measureText(maxText).width, top + 16);
        top += 36;
        // what the screen's legend says under its color bar, in short
        const { log, shared, clamped } = style.scale;
        g.fillText(`${log ? 'log' : 'linear'} color scale${shared ? ', shared by the panels' : ''}`, x, top, width);
        top += 18;
        if (clamped) { g.fillText('colors stop at the 99th percentile', x, top, width); top += 18; }
      }
      if (style && style.kind === 'gene') {
        const gv = style.gene;
        g.fillStyle = P.muted;
        g.fillText(`detected in ${fmt.int(gv.detected)} of ${fmt.int(data.n)} bins`, x, top, width);
        top += 18;
        if (gv.missing.length) {
          g.fillText(`not in the data: ${gv.missing.join(', ')}`, x, top, width);
          top += 18;
        }
      }
      for (const e of style ? style.legend : []) if (!e.off) swatchRow(e.color, e.label, e.count);
      if (ovl && ovl.layers.length) {
        top += 6;
        g.fillStyle = P.text;
        g.font = `600 12px ${FONT}`;
        g.fillText('Overlays', x, top, width);
        top += 20;
        for (const l of ovl.layers) swatchRow(l.color, l.label, l.count);
        if (ovl.several) swatchRow(SEVERAL, 'In two or more overlays', ovl.several);
      }
      g.restore();
    }

    /** This panel into the composite PNG: title, map (H&E, bins, frame, scale bar) and legend, at (ox, oy) css px. */
    function exportTo(g, ox, oy, ratio, legendW, titleH) {
      const P = EXPORT_PALETTE;
      g.setTransform(ratio, 0, 0, ratio, 0, 0);
      g.fillStyle = P.text;
      g.font = `600 15px ${FONT}`;
      g.textBaseline = 'middle';
      g.fillText(titleText(), ox + 12, oy + titleH / 2, Math.max(40, size.w + legendW - 24));
      const v = data && style ? getView() : null;
      if (!v || !size.w || !size.h) return;
      const t = xf(v);
      const mx = ox;
      const my = oy + titleH;
      g.save();
      g.beginPath();
      g.rect(mx, my, size.w, size.h);
      g.clip();
      const im = showImage();
      if (im) {
        g.imageSmoothingEnabled = true;
        drawImageUnder(g, im, { scale: t.scale, tx: t.tx + mx, ty: t.ty + my }, mx, my, size.w, size.h);
      }
      // the bins are drawn apart and copied in one go, so that a see-through layer has no darker seams
      const tmp = document.createElement('canvas');
      tmp.width = Math.round(size.w * ratio);
      tmp.height = Math.round(size.h * ratio);
      const tg = tmp.getContext('2d');
      const side = sidePx(t.scale);
      drawBins(tg, style.groups, data, t, ratio, side, P, dimOf());
      if (ovl) drawBins(tg, ovl.groups, data, t, ratio, side, P, 0);
      g.setTransform(1, 0, 0, 1, 0, 0);
      g.globalAlpha = im ? Math.max(0, 1 - store.image.opacity) : 1;
      g.drawImage(tmp, Math.round(mx * ratio), Math.round(my * ratio));
      g.restore();
      tmp.width = 0;
      tmp.height = 0;
      g.strokeStyle = P.neutral;
      g.lineWidth = 1;
      g.strokeRect(mx + 0.5, my + 0.5, size.w - 1, size.h - 1);
      drawScaleBar(g, t.scale, mx + 10, my + size.h - 10, P);
      drawLegendTo(g, mx + size.w + 16, my + 8, legendW - 32, P);
    }

    function destroy() {
      alive = false;
      loadToken += 1;
      clearTimeout(crispTimer);
      resizeObserver.unobserve(stage);
      resizeObserver.unobserve(el);
      tip.hidden = true;
      if (geneBox) geneBox.destroy();
      if (ovlGeneBox) ovlGeneBox.destroy();
      layer.width = 0;
      layer.height = 0;
    }

    return Object.assign(self, {
      el,
      stage,
      start: load,
      data: () => data,
      size: () => size,
      hasData: () => !!(data && data.valid),
      canExport: () => !!(data && data.valid && style && size.w && size.h),
      ownView: () => st.views[st.sample],
      adoptView(v) { st.views[st.sample] = { cx: v.cx, cy: v.cy, scale: v.scale }; },
      wants,
      resize,
      requestRender,
      queuePaint,
      crisp,
      step,
      imageChanged,
      ownScale,
      restyle,
      restyleIfScaleMoved,
      hideTip,
      markWidth,
      apply,
      fileTag,
      legendHeight,
      exportTo,
      destroy,
    });
  }

  markLayout(store.layout);
  start();

  return {
    el,
    show(patch = {}) {
      const i = Number.isInteger(patch.panel) && patch.panel >= 0 && patch.panel < LAYOUTS[LAYOUTS.length - 1] ? patch.panel : 0;
      if (panels[i]) {
        panels[i].apply(patch);
        return;
      }
      // a panel not shown (yet): remember the change for when it is
      store.panels[i] = patchState(store.panels[i] && typeof store.panels[i] === 'object' ? store.panels[i] : {}, patch);
    },
    destroy() {
      destroyed = true;
      if (startCtrl) startCtrl.abort();
      for (const r of inflight.values()) r.ctrl.abort();
      inflight.clear();
      cancelAnimationFrame(frame);
      cancelAnimationFrame(roomFrame);
      for (const p of panels) p.destroy();
      window.removeEventListener('resize', measureSoon);
      document.removeEventListener('scroll', onScroll, { capture: true });
      resizeObserver.disconnect();
      pageObserver.disconnect();
      themeObserver.disconnect();
      darkQuery.removeEventListener('change', onTheme);
    },
  };
}
