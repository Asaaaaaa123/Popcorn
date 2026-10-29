// spatial-viewer.js — the Viewer tab of a run: every bin of one section drawn at its position (x, y in µm, y down,
// equal aspect) and colored by a label, a positivity call or the UMI count (GET /api/runs/<id>/bins).
// Bins are drawn into an offscreen layer (one fillRect batch per color) only when the data, the zoom or the style
// change; panning and hovering only copy that layer, so sections of 80,000 bins and more stay smooth.
import { h, api, mount, icon, emptyState, errorPanel, loading, fmt, toast } from '../core.js';
import { downloadBlob } from './log-view.js';
import { formatNumber } from './table-viewer.js';

const FIELD_LABELS = {
  tissue_class: 'Tissue class',
  assigned_signature: 'Assigned signature',
  leiden_domain: 'Leiden domain',
  in_region_of_interest: 'In region of interest',
  total_umi: 'UMIs per bin',
};
const CATEGORY_ORDER = ['tissue_class', 'assigned_signature', 'leiden_domain'];
const DEFAULT_FIELDS = ['tissue_class', 'assigned_signature', 'leiden_domain', 'in_region_of_interest', 'total_umi'];
const LOG_FIELDS = new Set(['total_umi']);
// suffix of pooled calls in the per-bin table: _pooled_depthmatched (uncoupling.calls "pooled_dm") or _pooled ("pooled")
const POOLED_RE = /_pooled(_depthmatched)?$/;
// categories without a configured color (Leiden domains, ...): Tableau 20 without its grays
const PALETTE = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf',
  '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5'];
const GRAY_NAMES = { unassigned: '#dedede', ambiguous: '#dedede', removed: '#b8b8b8', unresolved: '#c4c4c4' };
const ROI_COLOR = '#1f78b4';
const POSITIVE_COLOR = '#d62728';
const VIRIDIS = ['#440154', '#472c7a', '#3b518b', '#2c718e', '#21908d', '#27ad81', '#5cc863', '#aadc32', '#fde725'];
const LEVELS = 64; // color steps of a numeric field (one fillRect batch each)
const NEUTRAL = 'neutral'; // group color token: the theme's neutral (negatives, missing values)
const MISSING = '__missing__'; // legend key of bins without a value
const MIN_SIDE = 1.5; // css px: bins never vanish when zoomed out
const HOVER_PX = 6;
const MARGIN = 160; // css px drawn around the viewport, so short pans need no redraw
const MAX_LAYER_PX = 16e6;
const EXPORT_PALETTE = { bg: '#ffffff', neutral: '#d0d5dc', text: '#18212c', muted: '#5b6675' };
const FONT = 'system-ui, -apple-system, "Segoe UI", Roboto, Arial, sans-serif';

let uidSeq = 0;
const uid = (prefix) => `${prefix}-${(uidSeq += 1)}`;
const isColor = (c) => typeof c === 'string' && /^#(?:[0-9a-f]{3}|[0-9a-f]{6})$/i.test(c.trim());

function rgb(hex) {
  const v = parseInt(hex.slice(1), 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}

/** Color at t in [0, 1] along a list of hex stops. */
function ramp(stops, t) {
  const x = Math.min(1, Math.max(0, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(x));
  const f = x - i;
  const a = rgb(stops[i]);
  const b = rgb(stops[i + 1]);
  return `rgb(${Math.round(a[0] + (b[0] - a[0]) * f)}, ${Math.round(a[1] + (b[1] - a[1]) * f)}, ${Math.round(a[2] + (b[2] - a[2]) * f)})`;
}
const LUT = Array.from({ length: LEVELS }, (_, k) => ramp(VIRIDIS, k / (LEVELS - 1)));

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
    if (f.type === 'bool') return [1, k === 'in_region_of_interest' ? 0 : (POOLED_RE.test(k) ? 2 : 1)];
    return [2, 0];
  };
  return [...fields.keys()].map((k, i) => ({ k, i, r: rank(k) }))
    .sort((a, b) => a.r[0] - b.r[0] || a.r[1] - b.r[1] || a.i - b.i)
    .map((e) => e.k);
}

/** GET /bins response -> typed arrays, bounding box and per-field summaries. */
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
      let min = Infinity;
      let max = -Infinity;
      for (let i = 0; i < n; i += 1) {
        const v = values[i];
        if (!ok[i] || !Number.isFinite(v)) continue;
        if (v < min) min = v;
        if (v > max) max = v;
      }
      fields.set(key, { key, type: 'number', values, min, max });
    }
  }
  return {
    sample: raw && raw.sample, n, x, y, ok, valid, bbox, fields, order: orderFields(fields),
    colors: raw && raw.colors && typeof raw.colors === 'object' ? raw.colors : {},
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

function boolColor(key, colors) {
  if (key === 'in_region_of_interest') return ROI_COLOR;
  if (key.startsWith('positive_')) {
    const slug = key.slice('positive_'.length).replace(POOLED_RE, '');
    const sig = colors.assigned_signature && colors.assigned_signature[slug];
    if (isColor(sig)) return sig.trim();
  }
  return POSITIVE_COLOR;
}

/**
 * Draw groups and legend of one field: groups = [{color, idx}] in drawing order (color may be NEUTRAL),
 * legend = [{key, label, color, count, off}], visible = 1 per drawn bin.
 */
function buildStyle(d, f, hidden) {
  const { n, ok } = d;
  const visible = new Uint8Array(n);
  if (!f) {
    const idx = new Int32Array(d.valid);
    let k = 0;
    for (let i = 0; i < n; i += 1) if (ok[i]) { idx[k] = i; k += 1; visible[i] = 1; }
    return { kind: 'none', groups: [{ color: NEUTRAL, idx }], legend: [], visible, shown: d.valid };
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
    entries = [
      { key: 'yes', label: f.key === 'in_region_of_interest' ? 'Inside' : 'Positive', color: boolColor(f.key, d.colors) },
      { key: 'no', label: f.key === 'in_region_of_interest' ? 'Outside' : 'Negative', color: NEUTRAL, gray: true },
    ];
    groupOf = (i) => (f.values[i] === 2 ? 2 : 1 - f.values[i]);
  } else {
    const log = LOG_FIELDS.has(f.key) && f.min >= 0;
    const tf = log ? Math.log1p : (v) => v;
    const lo = Number.isFinite(f.min) ? tf(f.min) : 0;
    const span = Number.isFinite(f.max) && f.max > f.min ? tf(f.max) - lo : 0;
    entries = LUT.map((color, k) => ({ key: `level-${k}`, color }));
    groupOf = (i) => {
      const v = f.values[i];
      if (!Number.isFinite(v)) return LEVELS;
      return span > 0 ? Math.min(LEVELS - 1, Math.max(0, Math.floor(((tf(v) - lo) / span) * LEVELS))) : 0;
    };
    if (Number.isFinite(f.min)) scale = { min: f.min, max: f.max, log };
  }
  const nG = entries.length + 1;
  const counts = new Int32Array(nG);
  const groupIdx = new Int32Array(n);
  for (let i = 0; i < n; i += 1) {
    if (!ok[i]) { groupIdx[i] = -1; continue; }
    const g = groupOf(i);
    groupIdx[i] = g;
    counts[g] += 1;
  }
  const off = (g) => (g === nG - 1 ? hidden.has(MISSING) : f.type !== 'number' && hidden.has(entries[g].key));
  const idx = Array.from({ length: nG }, (_, g) => (off(g) ? null : new Int32Array(counts[g])));
  const fill = new Int32Array(nG);
  let shown = 0;
  for (let i = 0; i < n; i += 1) {
    const g = groupIdx[i];
    if (g < 0 || !idx[g]) continue;
    idx[g][fill[g]] = i;
    fill[g] += 1;
    visible[i] = 1;
    shown += 1;
  }
  // drawing order: missing values, then gray / negative groups, then the rest (categories: larger first, so the rarer
  // ones stay visible on top; numbers: low to high)
  const order = [nG - 1];
  const rest = entries.map((e, g) => g);
  if (f.type === 'category') {
    rest.sort((a, b) => (entries[b].gray ? 1 : 0) - (entries[a].gray ? 1 : 0) || counts[b] - counts[a]);
  } else if (f.type === 'bool') {
    rest.reverse();
  }
  order.push(...rest);
  const groups = order.filter((g) => idx[g] && idx[g].length)
    .map((g) => ({ color: g === nG - 1 ? NEUTRAL : entries[g].color, idx: idx[g] }));
  const legend = f.type === 'number' ? [] : entries.map((e, g) => ({ key: e.key, label: e.label, color: e.color, count: counts[g], off: off(g) }));
  if (counts[nG - 1]) legend.push({ key: MISSING, label: 'no value', color: NEUTRAL, count: counts[nG - 1], off: off(nG - 1) });
  return { kind: f.type, groups, legend, visible, shown, scale };
}

// ------------------------------------------------------------------------------------------------ drawing
/**
 * Bins of `groups` as squares of `side` css px; view {scale, tx, ty} in css px; ratio = device px per css px.
 * Square edges snap outward to whole pixels: squares stay crisp and neighboring bins always meet (no seams).
 */
function drawBins(g, groups, d, view, ratio, side, palette) {
  const W = g.canvas.width;
  const H = g.canvas.height;
  const s = view.scale * ratio;
  const tx = view.tx * ratio;
  const ty = view.ty * ratio;
  const half = Math.max(0.5, (side * ratio) / 2);
  const { x, y } = d;
  for (const grp of groups) {
    g.fillStyle = grp.color === NEUTRAL ? palette.neutral : grp.color;
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

function readTheme(el) {
  const cs = getComputedStyle(el);
  const v = (name, fallback) => cs.getPropertyValue(name).trim() || fallback;
  return {
    bg: v('--plot-bg', '#ffffff'), neutral: v('--plot-neutral', '#c9d0da'),
    text: v('--text', '#18212c'), muted: v('--text-muted', '#566273'),
  };
}

// ------------------------------------------------------------------------------------------------ component
/**
 * spatialViewer({runId, store, binUm, emptyText}) -> {el, destroy()}
 * store (kept by the caller between visits): samples, data (Map sample -> prepared bins), sample, field,
 * hidden {field: [legend keys]}, pointSize, views {sample: {scale, tx, ty, w, h}}.
 * binUm: bin size of the run in µm (estimated from the positions when missing).
 */
export function spatialViewer({ runId, store, binUm: binUmOption, emptyText }) {
  if (!(store.data instanceof Map)) store.data = new Map();
  if (!store.views) store.views = {};
  if (!store.hidden) store.hidden = {};
  if (!Number.isFinite(store.pointSize)) store.pointSize = 1;
  const enc = encodeURIComponent(runId);

  let destroyed = false;
  let ctrl = null;
  let loadToken = 0;
  let data = null;
  let style = null;
  let grid = null;
  let pitch = 16;
  let view = null;
  let size = { w: 0, h: 0 };
  let dpr = window.devicePixelRatio || 1;
  let theme = null;
  let hover = -1;
  let frame = 0;
  let dirty = true;
  let crispTimer = null;
  let layerState = null; // the view the layer was drawn with
  const layer = document.createElement('canvas');
  const pointers = new Map();
  let drag = null;
  let pinch = null;

  // ---- DOM
  const sampleId = uid('sv-sample');
  const fieldId = uid('sv-field');
  const sizeId = uid('sv-size');
  const sampleSel = h('select', { class: 'select', id: sampleId, onChange: () => selectSample(sampleSel.value) });
  const fieldSel = h('select', { class: 'select', id: fieldId, onChange: () => selectField(fieldSel.value) });
  const sizeOut = h('output', { class: 'sv-size-out', for: sizeId });
  const sizeInput = h('input', {
    type: 'range', id: sizeId, min: 0.5, max: 3, step: 0.25, value: store.pointSize, class: 'sv-range',
    onInput: () => {
      store.pointSize = Number(sizeInput.value) || 1;
      sizeOut.textContent = `${store.pointSize}×`;
      requestRender();
    },
  });
  sizeOut.textContent = `${store.pointSize}×`;
  const saveBtn = h('button', { class: 'btn small', type: 'button', onClick: savePng, disabled: true },
    icon('download', { size: 14 }), 'Save PNG');
  const toolbar = h('div', { class: 'sv-toolbar' },
    h('div', { class: 'sv-control' }, h('label', { for: sampleId }, 'Section'), sampleSel),
    h('div', { class: 'sv-control sv-control-wide' }, h('label', { for: fieldId }, 'Color by'), fieldSel),
    h('div', { class: 'sv-control' }, h('label', { for: sizeId }, 'Point size'), h('div', { class: 'sv-range-row' }, sizeInput, sizeOut)),
    h('span', { class: 'spacer' }),
    saveBtn);

  const canvas = h('canvas', {
    class: 'sv-canvas', tabIndex: 0, role: 'img',
    'aria-label': 'Map of the bins of the section. Drag or use the arrow keys to pan; scroll, + and − to zoom; 0 fits the section.',
  });
  const ctx = canvas.getContext('2d');
  const tip = h('div', { class: 'sv-tip', hidden: true, role: 'tooltip' });
  const overlay = h('div', { class: 'sv-overlay', hidden: true });
  const zoomBtn = (iconName, title, onClick) => h('button', { class: 'btn small icon-only', type: 'button', title, 'aria-label': title, onClick }, icon(iconName, { size: 15 }));
  const stage = h('div', { class: 'sv-stage' }, canvas,
    h('div', { class: 'sv-zoom' },
      zoomBtn('zoom-in', 'Zoom in (+)', () => zoomBy(1.6)),
      zoomBtn('zoom-out', 'Zoom out (−)', () => zoomBy(1 / 1.6)),
      zoomBtn('maximize', 'Fit the section (0)', resetView)),
    tip, overlay);
  const legendEl = h('div', { class: 'sv-legend-body' });
  const side = h('aside', { class: 'card sv-side' }, legendEl);
  const el = h('div', { class: 'sv' }, toolbar, h('div', { class: 'sv-body' }, stage, side));

  // ---- events
  canvas.addEventListener('pointerdown', onDown);
  canvas.addEventListener('pointermove', onMove);
  canvas.addEventListener('pointerup', onUp);
  canvas.addEventListener('pointercancel', onUp);
  canvas.addEventListener('pointerleave', () => { if (!pointers.size) setHover(-1); });
  canvas.addEventListener('wheel', onWheel, { passive: false });
  canvas.addEventListener('dblclick', (e) => {
    if (!view) return;
    const p = local(e);
    zoomAt(e.shiftKey ? 1 / 2 : 2, p.x, p.y);
    afterViewChange(true);
  });
  canvas.addEventListener('keydown', onKey);
  canvas.addEventListener('blur', () => setHover(-1));

  const resizeObserver = new ResizeObserver(() => resize());
  resizeObserver.observe(stage);
  const darkQuery = window.matchMedia('(prefers-color-scheme: dark)');
  const onTheme = () => { theme = null; requestRender(); };
  darkQuery.addEventListener('change', onTheme);
  const themeObserver = new MutationObserver(onTheme);
  themeObserver.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

  // ---- view geometry
  const sidePx = () => (view ? Math.max(MIN_SIDE, pitch * view.scale * store.pointSize) : MIN_SIDE);

  function local(e) {
    const r = canvas.getBoundingClientRect();
    return { x: e.clientX - r.left, y: e.clientY - r.top };
  }

  function fitScale() {
    const b = data.bbox;
    const pad = 24;
    const bw = b.x1 - b.x0 + pitch;
    const bh = b.y1 - b.y0 + pitch;
    const s = Math.min((size.w - 2 * pad) / bw, (size.h - 2 * pad) / bh);
    return Number.isFinite(s) && s > 0 ? s : 1;
  }

  function fitView() {
    const s = fitScale();
    const b = data.bbox;
    view = { scale: s, tx: size.w / 2 - ((b.x0 + b.x1) / 2) * s, ty: size.h / 2 - ((b.y0 + b.y1) / 2) * s };
    saveView();
  }

  function initView() {
    const saved = store.views[store.sample];
    if (saved && Number.isFinite(saved.scale) && saved.scale > 0) {
      view = { scale: saved.scale, tx: saved.tx + (size.w - saved.w) / 2, ty: saved.ty + (size.h - saved.h) / 2 };
      saveView();
    } else {
      fitView();
    }
  }

  function saveView() {
    if (view && store.sample !== undefined && size.w && size.h) store.views[store.sample] = { ...view, w: size.w, h: size.h };
  }

  function zoomAt(factor, cx, cy) {
    if (!view || !data) return;
    const fit = fitScale();
    const s1 = Math.min(Math.max(view.scale * factor, fit / 4), Math.max(fit, 240 / pitch));
    const f = s1 / view.scale;
    view.tx = cx - (cx - view.tx) * f;
    view.ty = cy - (cy - view.ty) * f;
    view.scale = s1;
  }

  function zoomBy(factor) {
    zoomAt(factor, size.w / 2, size.h / 2);
    afterViewChange(true);
  }

  function resetView() {
    if (!data || !data.valid) return;
    fitView();
    setHover(-1);
    requestRender();
  }

  /** The view moved: copy the layer now; a zoom also redraws crisply once it settles. */
  function afterViewChange(zooming) {
    saveView();
    queuePaint();
    if (zooming) {
      clearTimeout(crispTimer);
      crispTimer = setTimeout(requestRender, 140);
    }
  }

  function resize() {
    // the drawable area inside the stage's border (the canvas fills it): the border box made the bitmap 2 px too large,
    // stretching the map unevenly and offsetting bins from hover and the scale bar
    const w = Math.max(0, Math.floor(stage.clientWidth));
    const hh = Math.max(0, Math.floor(stage.clientHeight));
    const ratio = window.devicePixelRatio || 1;
    if (!w || !hh) return; // hidden or detached: keep the current geometry
    if (w === size.w && hh === size.h && ratio === dpr) return;
    if (view && size.w && size.h) {
      view.tx += (w - size.w) / 2;
      view.ty += (hh - size.h) / 2;
    }
    size = { w, h: hh };
    dpr = ratio;
    canvas.width = Math.max(1, Math.round(w * dpr));
    canvas.height = Math.max(1, Math.round(hh * dpr));
    if (data && data.valid && !view && w && hh) initView();
    saveView();
    requestRender();
  }

  // ---- rendering
  function requestRender() {
    dirty = true;
    queuePaint();
  }

  function queuePaint() {
    if (!frame && !destroyed) frame = requestAnimationFrame(tick);
  }

  function tick() {
    frame = 0;
    if (destroyed) return;
    if ((window.devicePixelRatio || 1) !== dpr) resize(); // moved to a screen with another pixel density
    if (!size.w || !size.h) return;
    if (!theme && stage.isConnected) theme = readTheme(stage);
    if (dirty) renderLayer();
    paint();
  }

  /** Draw every visible bin into the offscreen layer for the current view (plus a margin around it). */
  function renderLayer() {
    dirty = false;
    clearTimeout(crispTimer);
    if (!data || !style || !view) { layerState = null; return; }
    let margin = MARGIN;
    while (margin > 0 && (size.w + 2 * margin) * (size.h + 2 * margin) * dpr * dpr > MAX_LAYER_PX) margin = Math.floor(margin / 2);
    const lw = Math.round((size.w + 2 * margin) * dpr);
    const lh = Math.round((size.h + 2 * margin) * dpr);
    if (layer.width !== lw || layer.height !== lh) { layer.width = lw; layer.height = lh; }
    const g = layer.getContext('2d');
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.clearRect(0, 0, lw, lh);
    const palette = theme || EXPORT_PALETTE;
    drawBins(g, style.groups, data, { scale: view.scale, tx: view.tx + margin, ty: view.ty + margin }, dpr, sidePx(), palette);
    layerState = { scale: view.scale, tx: view.tx, ty: view.ty, margin, dpr, w: size.w, h: size.h };
  }

  /** Visible canvas = background + the layer moved / scaled to the current view + hover marker + scale bar. */
  function paint() {
    const palette = theme || EXPORT_PALETTE;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.fillStyle = palette.bg;
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    if (!data || !view) return;
    if (layerState && layerState.dpr === dpr) {
      const k = view.scale / layerState.scale;
      const left = view.tx - (layerState.tx + layerState.margin) * k; // css px
      const top = view.ty - (layerState.ty + layerState.margin) * k;
      const lw = (layerState.w + 2 * layerState.margin) * k;
      const lh = (layerState.h + 2 * layerState.margin) * k;
      ctx.imageSmoothingEnabled = k !== 1;
      ctx.drawImage(layer, left * dpr, top * dpr, lw * dpr, lh * dpr);
      // the layer no longer covers the viewport (long pan, zooming out): redraw on the next frame
      if (left > 0.5 || top > 0.5 || left + lw < size.w - 0.5 || top + lh < size.h - 0.5) requestRender();
    } else if (style) {
      requestRender();
    }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    if (hover >= 0 && hover < data.n) {
      const half = Math.max(sidePx() / 2, 3) + 2;
      const cx = data.x[hover] * view.scale + view.tx;
      const cy = data.y[hover] * view.scale + view.ty;
      ctx.lineWidth = 3;
      ctx.strokeStyle = palette.bg;
      ctx.strokeRect(cx - half, cy - half, half * 2, half * 2);
      ctx.lineWidth = 1.5;
      ctx.strokeStyle = palette.text;
      ctx.strokeRect(cx - half, cy - half, half * 2, half * 2);
    }
    drawScaleBar(ctx, view.scale, 10, size.h - 10, palette);
  }

  // ---- pointer, wheel and keyboard
  function onDown(e) {
    if (!view || !data) return;
    if (e.pointerType === 'mouse' && e.button !== 0) return;
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
        view.tx += now.x - pinch.x;
        view.ty += now.y - pinch.y;
        pinch = now;
        afterViewChange(true);
      } else if (drag) {
        const dx = p.x - drag.x;
        const dy = p.y - drag.y;
        if (!drag.moved && Math.abs(dx) + Math.abs(dy) < 3) return;
        if (!drag.moved) { drag.moved = true; stage.classList.add('dragging'); setHover(-1); }
        view.tx += dx;
        view.ty += dy;
        drag.x = p.x;
        drag.y = p.y;
        afterViewChange(false);
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
      requestRender();
    } else if (e.type === 'pointerup') {
      const p = local(e); // a click or tap shows the bin under it
      hoverAt(p.x, p.y);
    }
  }

  function onWheel(e) {
    if (!view || !data) return;
    e.preventDefault();
    const unit = e.deltaMode === 1 ? 16 : e.deltaMode === 2 ? 400 : 1;
    const p = local(e);
    // trackpad pinches arrive as ctrl + wheel with small deltas
    zoomAt(Math.exp(-e.deltaY * unit * (e.ctrlKey ? 0.01 : 0.0015)), p.x, p.y);
    afterViewChange(true);
    hoverAt(p.x, p.y);
  }

  function onKey(e) {
    if (!view || !data || e.altKey || e.ctrlKey || e.metaKey) return;
    const step = 60;
    switch (e.key) {
      case '+': case '=': zoomBy(1.4); break;
      case '-': case '_': zoomBy(1 / 1.4); break;
      case '0': resetView(); break;
      case 'ArrowLeft': view.tx += step; afterViewChange(true); break;
      case 'ArrowRight': view.tx -= step; afterViewChange(true); break;
      case 'ArrowUp': view.ty += step; afterViewChange(true); break;
      case 'ArrowDown': view.ty -= step; afterViewChange(true); break;
      case 'Escape': setHover(-1); return;
      default: return;
    }
    e.preventDefault();
  }

  // ---- hover
  function hoverAt(px, py) {
    if (!data || !style || !view || !grid) { setHover(-1); return; }
    const r = Math.max(HOVER_PX, sidePx() * 0.71 + 1) / view.scale;
    const i = nearest(grid, data, (px - view.tx) / view.scale, (py - view.ty) / view.scale, r, style.visible);
    setHover(i, px, py);
  }

  function valueText(f, i) {
    if (f.type === 'category') return f.codes[i] < 0 ? '—' : f.categories[f.codes[i]];
    if (f.type === 'bool') return f.values[i] === 2 ? '—' : (f.values[i] ? 'yes' : 'no');
    return Number.isFinite(f.values[i]) ? formatNumber(f.values[i]) : '—';
  }

  function setHover(i, px, py) {
    if (i !== hover) {
      hover = i;
      queuePaint();
      if (i < 0) {
        tip.hidden = true;
        return;
      }
      const rows = data.order.map((key) => {
        const f = data.fields.get(key);
        const on = key === store.field;
        return [h('dt', { class: on ? 'on' : null }, fieldLabel(key)), h('dd', { class: on ? 'on' : null }, valueText(f, i))];
      });
      mount(tip,
        h('div', { class: 'sv-tip-head' }, `x ${fmt.int(data.x[i])} µm · y ${fmt.int(data.y[i])} µm`),
        rows.length ? h('dl', { class: 'sv-tip-grid' }, rows) : null);
      tip.hidden = false;
    }
    if (i >= 0 && px !== undefined) {
      const pad = 14;
      const tw = tip.offsetWidth;
      const th = tip.offsetHeight;
      let left = px + pad;
      let top = py + pad;
      if (left + tw > size.w - 4) left = px - pad - tw;
      if (top + th > size.h - 4) top = size.h - th - 4;
      tip.style.transform = `translate(${Math.max(4, left)}px, ${Math.max(4, top)}px)`;
    }
  }

  // ---- legend
  function paintLegend() {
    if (!data) { mount(legendEl); return; }
    const f = data.fields.get(store.field);
    const title = f ? fieldLabel(f.key) : 'Bins';
    const parts = [h('h3', { class: 'sv-legend-title' }, title)];
    if (!style) {
      mount(legendEl, parts);
      return;
    }
    if (style.kind === 'number' && style.scale) {
      const bar = h('canvas', { class: 'sv-colorbar', width: 256, height: 1, 'aria-hidden': 'true' });
      const g = bar.getContext('2d');
      for (let k = 0; k < 256; k += 1) {
        g.fillStyle = LUT[Math.min(LEVELS - 1, Math.floor((k / 256) * LEVELS))];
        g.fillRect(k, 0, 1, 1);
      }
      const { min, max, log } = style.scale;
      const mid = log ? Math.expm1((Math.log1p(min) + Math.log1p(max)) / 2) : (min + max) / 2;
      parts.push(bar, h('div', { class: 'sv-colorbar-labels' },
        h('span', null, formatNumber(min)), h('span', null, formatNumber(Math.round(mid * 100) / 100)), h('span', null, formatNumber(max))),
      h('div', { class: 'sv-note' }, log ? 'Logarithmic color scale (log 1 + x)' : 'Linear color scale'));
    }
    if (style.legend.length) {
      parts.push(h('div', { class: 'sv-legend', role: 'group', 'aria-label': `${title}: click an entry to hide or show it` },
        style.legend.map((e) => h('button', {
          class: e.off ? 'sv-item off' : 'sv-item', type: 'button', 'aria-pressed': e.off ? 'false' : 'true',
          dataset: { key: e.key }, title: 'Click to hide or show; Alt+click to show only this entry',
          onClick: (ev) => toggle(e.key, ev.altKey),
        },
        h('span', { class: 'swatch', style: { background: e.color === NEUTRAL ? 'var(--plot-neutral)' : e.color } }),
        h('span', { class: 'sv-item-label' }, e.label),
        h('span', { class: 'sv-item-count' }, fmt.int(e.count))))));
    }
    const anyOff = style.legend.some((e) => e.off);
    parts.push(h('div', { class: 'sv-legend-foot' },
      h('span', null, `${fmt.int(style.shown)} of ${fmt.int(data.valid)} bins shown`),
      anyOff ? h('button', { class: 'btn ghost small', type: 'button', onClick: showAll }, 'Show all') : null));
    parts.push(h('p', { class: 'sv-note' }, 'Drag to pan, scroll or double-click to zoom; hover a bin for its values.'));
    mount(legendEl, parts);
  }

  function toggle(key, solo) {
    const f = data && data.fields.get(store.field);
    if (!f || !style) return;
    let hidden = new Set(store.hidden[f.key] || []);
    if (solo) {
      const keys = style.legend.map((e) => e.key);
      const alone = keys.every((k) => (k === key ? !hidden.has(k) : hidden.has(k)));
      hidden = alone ? new Set() : new Set(keys.filter((k) => k !== key));
    } else if (hidden.has(key)) {
      hidden.delete(key);
    } else {
      hidden.add(key);
    }
    store.hidden[f.key] = [...hidden];
    restyle();
    const again = [...legendEl.querySelectorAll('.sv-item')].find((b) => b.dataset.key === key);
    if (again) again.focus();
  }

  function showAll() {
    const f = data && data.fields.get(store.field);
    if (!f) return;
    store.hidden[f.key] = [];
    restyle();
    canvas.focus();
  }

  function restyle() {
    if (!data) { style = null; paintLegend(); requestRender(); return; }
    const f = data.fields.get(store.field) || null;
    style = buildStyle(data, f, new Set(f ? store.hidden[f.key] || [] : []));
    if (hover >= 0 && !style.visible[hover]) setHover(-1);
    paintLegend();
    requestRender();
  }

  // ---- data
  function showOverlay(content) {
    mount(overlay, content);
    overlay.hidden = false;
  }

  function hideOverlay() {
    mount(overlay);
    overlay.hidden = true;
  }

  function showEmpty() {
    mount(el, h('div', { class: 'card' }, emptyState({
      icon: 'map', title: 'No per-bin data',
      text: emptyText || 'This run has no per-bin table (per_bin_signature_calls), so there is nothing to map.',
    })));
  }

  function fillFieldSelect() {
    const groups = [
      { label: 'Labels', keys: data.order.filter((k) => data.fields.get(k).type === 'category') },
      { label: 'Region and positivity calls', keys: data.order.filter((k) => data.fields.get(k).type === 'bool') },
      { label: 'Numbers', keys: data.order.filter((k) => data.fields.get(k).type === 'number') },
    ].filter((g) => g.keys.length);
    mount(fieldSel, groups.map((g) => h('optgroup', { label: g.label },
      g.keys.map((k) => h('option', { value: k }, fieldLabel(k))))));
    fieldSel.disabled = !groups.length;
    if (!groups.length) mount(fieldSel, h('option', { value: '' }, 'Position only'));
  }

  function selectField(key) {
    if (!data) return;
    store.field = data.fields.has(key) ? key : (DEFAULT_FIELDS.find((k) => data.fields.has(k)) || data.order[0] || '');
    if (fieldSel.value !== store.field) fieldSel.value = store.field;
    restyle();
  }

  async function selectSample(sid) {
    const token = ++loadToken;
    if (ctrl) ctrl.abort();
    ctrl = null;
    store.sample = sid;
    if (sampleSel.value !== sid) sampleSel.value = sid;
    setHover(-1);
    let d = store.data.get(sid);
    if (!d) {
      const meta = (store.samples || []).find((s) => String(s.id) === sid);
      showOverlay(loading(meta && Number.isFinite(meta.n) ? `Loading ${fmt.int(meta.n)} bins…` : 'Loading bins…'));
      saveBtn.disabled = true;
      const my = new AbortController();
      ctrl = my;
      try {
        const raw = await api(`/api/runs/${enc}/bins`, { query: { sample: sid }, signal: my.signal });
        if (destroyed || token !== loadToken) return;
        d = prepare(raw);
        store.data.set(sid, d);
      } catch (err) {
        if (destroyed || token !== loadToken || (err && err.name === 'AbortError')) return;
        data = null;
        style = null;
        paintLegend();
        requestRender();
        showOverlay(h('div', { class: 'sv-overlay-msg' }, errorPanel(err, { retry: () => selectSample(sid) })));
        return;
      } finally {
        if (ctrl === my) ctrl = null;
      }
    }
    data = d;
    if (!d.pitch) d.pitch = Number.isFinite(binUmOption) && binUmOption > 0 ? binUmOption : (d.valid > 1 ? estimatePitch(d) : 16);
    pitch = d.pitch;
    if (d.valid && !d.grid) d.grid = buildIndex(d, pitch * 2);
    grid = d.grid || null;
    fillFieldSelect();
    view = null;
    layerState = null;
    if (d.valid && size.w && size.h) initView();
    selectField(store.field);
    saveBtn.disabled = !d.valid;
    if (d.valid) hideOverlay();
    else showOverlay(h('div', { class: 'sv-overlay-msg' }, emptyState({ icon: 'map', title: 'No bins to draw', text: 'This section has no bins with coordinates.' })));
  }

  async function start() {
    let samples = store.samples;
    if (!samples) {
      showOverlay(loading('Loading sections…'));
      const my = new AbortController();
      ctrl = my;
      try {
        const res = await api(`/api/runs/${enc}/bins/samples`, { signal: my.signal });
        if (destroyed) return;
        samples = (Array.isArray(res) ? res : []).filter((s) => s && s.id !== null && s.id !== undefined);
        if (samples.length) store.samples = samples;
      } catch (err) {
        if (destroyed || (err && err.name === 'AbortError')) return;
        if (err.status === 404) { showEmpty(); return; }
        showOverlay(h('div', { class: 'sv-overlay-msg' }, errorPanel(err, { retry: start })));
        return;
      } finally {
        if (ctrl === my) ctrl = null;
      }
    }
    if (!samples.length) { showEmpty(); return; }
    mount(sampleSel, samples.map((s) => h('option', { value: String(s.id) },
      `${s.label && s.label !== s.id ? `${s.label} (${s.id})` : s.id}${Number.isFinite(s.n) ? ` — ${fmt.int(s.n)} bins` : ''}`)));
    sampleSel.disabled = samples.length < 2;
    const wanted = samples.some((s) => String(s.id) === store.sample) ? store.sample : String(samples[0].id);
    await selectSample(wanted);
  }

  // ---- export
  function drawLegendTo(g, x, y, width, palette) {
    const f = data.fields.get(store.field);
    g.save();
    g.fillStyle = palette.text;
    g.font = `600 13px ${FONT}`;
    g.textBaseline = 'top';
    g.fillText(f ? fieldLabel(f.key) : 'Bins', x, y, width);
    let top = y + 24;
    g.font = `12px ${FONT}`;
    if (style.kind === 'number' && style.scale) {
      for (let k = 0; k < width; k += 1) {
        g.fillStyle = LUT[Math.min(LEVELS - 1, Math.floor((k / width) * LEVELS))];
        g.fillRect(x + k, top, 1, 12);
      }
      g.fillStyle = palette.muted;
      g.fillText(formatNumber(style.scale.min), x, top + 16);
      const maxText = formatNumber(style.scale.max);
      g.fillText(maxText, x + width - g.measureText(maxText).width, top + 16);
      top += 36;
      if (style.scale.log) { g.fillText('log color scale', x, top); top += 20; }
    }
    for (const e of style.legend) {
      if (e.off) continue;
      g.fillStyle = e.color === NEUTRAL ? palette.neutral : e.color;
      g.fillRect(x, top + 1, 12, 12);
      g.fillStyle = palette.text;
      const count = fmt.int(e.count);
      const cw = g.measureText(count).width;
      g.fillText(e.label, x + 18, top, Math.max(10, width - 26 - cw));
      g.fillStyle = palette.muted;
      g.fillText(count, x + width - cw, top);
      top += 19;
    }
    g.restore();
  }

  function savePng() {
    if (!data || !style || !view || !size.w || !size.h) return;
    const ratio = Math.min(3, Math.max(2, dpr));
    const legendW = 240;
    const titleH = 40;
    const legendH = 80 + style.legend.length * 19;
    const W = size.w + legendW;
    const H = Math.max(size.h, legendH) + titleH;
    const out = document.createElement('canvas');
    out.width = Math.round(W * ratio);
    out.height = Math.round(H * ratio);
    const g = out.getContext('2d');
    if (!g) { toast('This browser cannot create the image', 'error'); return; }
    g.fillStyle = EXPORT_PALETTE.bg;
    g.fillRect(0, 0, out.width, out.height);
    g.save();
    g.beginPath();
    g.rect(0, titleH * ratio, size.w * ratio, size.h * ratio);
    g.clip();
    drawBins(g, style.groups, data, { scale: view.scale, tx: view.tx, ty: view.ty + titleH }, ratio, sidePx(), EXPORT_PALETTE);
    g.restore();
    g.setTransform(ratio, 0, 0, ratio, 0, 0);
    const meta = (store.samples || []).find((s) => String(s.id) === store.sample);
    const f = data.fields.get(store.field);
    g.fillStyle = EXPORT_PALETTE.text;
    g.font = `600 15px ${FONT}`;
    g.textBaseline = 'middle';
    g.fillText(`${meta && meta.label ? meta.label : store.sample} — ${f ? fieldLabel(f.key) : 'bins'}`, 12, titleH / 2, W - 24);
    drawScaleBar(g, view.scale, 10, titleH + size.h - 10, EXPORT_PALETTE);
    drawLegendTo(g, size.w + 16, titleH + 8, legendW - 32, EXPORT_PALETTE);
    const name = `${runId}_${store.sample}_${store.field || 'bins'}.png`.replace(/[^\w.-]+/g, '_');
    out.toBlob((blob) => {
      if (blob) downloadBlob(name, blob);
      else toast('The image could not be created', 'error');
    }, 'image/png');
  }

  start();

  return {
    el,
    destroy() {
      destroyed = true;
      if (ctrl) ctrl.abort();
      cancelAnimationFrame(frame);
      clearTimeout(crispTimer);
      resizeObserver.disconnect();
      themeObserver.disconnect();
      darkQuery.removeEventListener('change', onTheme);
    },
  };
}
