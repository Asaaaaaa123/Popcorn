// roi-editor.js — the Region of interest tab: polygons drawn on a section's H&E image (spatial/tissue_hires_image.png,
// coordinates in hires pixels), saved as GeoJSON with the project and, when cfg.roi.use_drawn is on, intersected with
// the transcriptional region of interest by the analysis (bone_pending_masks.py).
import { h, svg, api, apiUrl, ApiError, mount, icon, button, switchToggle, emptyState, errorPanel, banner, toast, confirmDialog, loading, fmt } from '../core.js';
import { iconButton } from './chips.js';

const HANDLE_PX = 5;      // vertex handle radius on screen
const SNAP_PX = 10;       // a click this close to the first corner closes the polygon
const CLICK_SLOP_PX = 4;  // a press that moves further than this pans instead of clicking
const DOUBLE_MS = 450;
const MAX_UNDO = 200;
const MAX_ZOOM = 40;      // screen pixels per image pixel
const BIN_CANVAS_MAX = 2400;

// ------------------------------------------------------------------------------------------------ geometry
export function ringArea(ring) {
  let a = 0;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i, i += 1) a += ring[j][0] * ring[i][1] - ring[i][0] * ring[j][1];
  return Math.abs(a) / 2;
}

export function pointInRing([x, y], ring) {
  let inside = false;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i, i += 1) {
    const [xi, yi] = ring[i];
    const [xj, yj] = ring[j];
    if ((yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

/** Inside the outer ring and not inside one of its holes. */
export function pointInPolygon(p, poly) {
  const [outer, ...holes] = poly.rings;
  return !!outer && pointInRing(p, outer) && !holes.some((r) => pointInRing(p, r));
}

export function polygonArea(poly) {
  const [outer, ...holes] = poly.rings;
  return outer ? Math.max(0, ringArea(outer) - holes.reduce((s, r) => s + ringArea(r), 0)) : 0;
}

function ringFrom(coords) {
  if (!Array.isArray(coords)) return null;
  const pts = coords.filter((c) => Array.isArray(c) && c.length >= 2 && Number.isFinite(Number(c[0])) && Number.isFinite(Number(c[1])))
    .map((c) => [Number(c[0]), Number(c[1])]);
  if (pts.length > 1) {
    const a = pts[0];
    const b = pts[pts.length - 1];
    if (a[0] === b[0] && a[1] === b[1]) pts.pop();
  }
  return pts.length >= 3 ? pts : null;
}

/** GeoJSON (FeatureCollection, Feature, GeometryCollection, Polygon, MultiPolygon) -> [{rings: [outer, ...holes]}], open rings. */
export function polygonsFromGeoJSON(gj) {
  const out = [];
  const addPolygon = (coords) => {
    if (!Array.isArray(coords)) return;
    const outer = ringFrom(coords[0]);
    if (outer) out.push({ rings: [outer, ...coords.slice(1).map(ringFrom).filter(Boolean)] });
  };
  const visit = (g) => {
    if (!g || typeof g !== 'object') return;
    if (g.type === 'FeatureCollection') (g.features || []).forEach(visit);
    else if (g.type === 'Feature') visit(g.geometry);
    else if (g.type === 'GeometryCollection') (g.geometries || []).forEach(visit);
    else if (g.type === 'Polygon') addPolygon(g.coordinates);
    else if (g.type === 'MultiPolygon') (g.coordinates || []).forEach(addPolygon);
  };
  visit(gj);
  return out;
}

const round2 = (v) => Math.round(v * 100) / 100;

/** [{rings}] -> FeatureCollection of Polygon features with closed rings and properties {frame}. */
export function polygonsToGeoJSON(polys, frame = 'hires') {
  return {
    type: 'FeatureCollection',
    features: polys.map((p) => ({
      type: 'Feature',
      properties: { frame },
      geometry: {
        type: 'Polygon',
        coordinates: p.rings.map((r) => {
          const ring = r.map(([x, y]) => [round2(x), round2(y)]);
          ring.push([...ring[0]]);
          return ring;
        }),
      },
    })),
  };
}

const pathD = (poly) => poly.rings.map((r) => `M${r.map(([x, y]) => `${x} ${y}`).join('L')}Z`).join('');
const plural = (n, one, many = `${one}s`) => `${n} ${n === 1 ? one : many}`;

// ------------------------------------------------------------------------------------------------ editor
/**
 * store: per-project editor state kept by the caller across renders ({sample, opacity, showBins, entries, meta, bins});
 * savedRois: Set of sample ids with a saved region; getVersion(): cache key of the images (changes when the project is
 * saved); sampleStatus(sid) -> {ok, title, text, note, canSave}; onConfigChange() after cfg.roi edits; onRoisChange()
 * after any region edit; onSaveProject() -> Promise<boolean>. Returns {el, destroy()}.
 */
export function roiEditor({ projectId, cfg, store, savedRois, getVersion, sampleStatus, onConfigChange, onRoisChange, onSaveProject }) {
  store.entries = store.entries || {};
  store.meta = store.meta || {};
  store.bins = store.bins || {};
  if (typeof store.opacity !== 'number') store.opacity = 1;
  const samples = () => (Array.isArray(cfg.samples) ? cfg.samples : []);
  if (!samples().some((s) => s.id === store.sample)) store.sample = samples().length ? samples()[0].id : null;

  const enc = encodeURIComponent;
  const sampleBase = (sid) => `/api/projects/${enc(projectId)}/samples/${enc(sid)}`;
  const roiUrl = (sid) => `/api/projects/${enc(projectId)}/roi/${enc(sid)}`;
  const entryOf = (sid) => {
    if (!store.entries[sid]) store.entries[sid] = { polys: [], saved: '[]', undo: [], redo: [], loaded: false, view: null };
    return store.entries[sid];
  };
  const entry = () => entryOf(store.sample);
  const isUnsaved = (sid) => {
    const en = store.entries[sid];
    return !!en && en.loaded && JSON.stringify(en.polys) !== en.saved;
  };

  let destroyed = false;
  let loadSeq = 0;
  let ready = false;
  let W = 0;
  let H = 0;
  let meta = null;
  let mode = 'select';
  let drawing = null;
  let selected = -1;
  let drag = null;
  let hover = null;
  let lastClick = null;
  let saving = false;
  let blocked = false;  // sampleStatus() kept the section from loading
  let shownNote = '';

  // ---------------------------------------------------------------- DOM
  const sampleSel = h('select', { class: 'select roi-sample', 'aria-label': 'Section' });
  const selectBtn = button({ label: 'Select', small: true, onClick: () => setMode('select') });
  const drawBtn = button({ label: 'Draw', icon: 'polygon', small: true, onClick: () => setMode('draw') });
  selectBtn.title = 'Select and adjust polygons; drag to pan (V)';
  drawBtn.title = 'Draw a polygon (D)';
  const zoomOut = iconButton('zoom-out', 'Zoom out (-)', () => zoomCenter(1 / 1.4));
  const zoomIn = iconButton('zoom-in', 'Zoom in (+)', () => zoomCenter(1.4));
  const fitBtn = iconButton('maximize', 'Fit the image (0)', () => { if (ready) fit(); });
  const opacity = h('input', { type: 'range', min: 0, max: 100, step: 5, value: Math.round(store.opacity * 100), class: 'roi-opacity', 'aria-label': 'Image opacity' });
  const binsSwitch = switchToggle({ checked: !!store.showBins, label: 'Bins', title: 'Show the in-tissue bins of the analysis bin size', onChange: (v) => { store.showBins = v; loadBins(); } });
  const binsInput = binsSwitch.querySelector('input');

  const img = h('img', { class: 'roi-img', alt: '', draggable: 'false' });
  const canvas = h('canvas', { class: 'roi-bins hidden' });
  const gPolys = svg('g');
  const gHandles = svg('g');
  const gDraw = svg('g');
  const overlay = svg('svg', { class: 'roi-svg' }, gPolys, gHandles, gDraw);
  const stage = h('div', { class: 'roi-stage' }, img, canvas, overlay);
  const msg = h('div', { class: 'roi-msg' });
  const note = h('div', { class: 'roi-note hidden' });
  const viewport = h('div', {
    class: 'roi-viewport', tabIndex: 0, role: 'application', dataset: { mode },
    'aria-label': 'H&E image. D draws, V selects, Enter closes the polygon, Escape cancels, Backspace removes the last corner, Delete removes the selected polygon, Ctrl+Z undoes, arrow keys pan, + and - zoom, 0 fits.',
  }, stage, note, msg);
  const hint = h('span', { class: 'roi-hint' });
  const coords = h('span', { class: 'roi-coords mono' });
  const sideRegion = h('div', { class: 'stack' });
  const sideAnalysis = h('div', { class: 'stack' });

  opacity.addEventListener('input', () => {
    store.opacity = Number(opacity.value) / 100;
    img.style.opacity = String(store.opacity);
  });
  img.style.opacity = String(store.opacity);
  sampleSel.addEventListener('change', () => {
    if (drawing && drawing.length >= 3) closeDrawing();
    store.sample = sampleSel.value;
    loadSample();
  });

  const el = h('div', { class: 'roi-layout' },
    h('section', { class: 'card roi-main' },
      h('div', { class: 'toolbar roi-toolbar' },
        h('label', { class: 'roi-tool-label' }, 'Section', sampleSel),
        h('div', { class: 'btn-group', role: 'group', 'aria-label': 'Tool' }, selectBtn, drawBtn),
        h('div', { class: 'row roi-zoom' }, zoomOut, zoomIn, fitBtn),
        h('label', { class: 'roi-tool-label' }, 'Image', opacity),
        binsSwitch),
      viewport,
      h('div', { class: 'roi-status' }, hint, coords)),
    h('div', { class: 'roi-side stack' },
      card('Region', sideRegion),
      card('Analysis', sideAnalysis),
      card('Keys', h('dl', { class: 'kv roi-keys' }, [
        ['Click', 'add a corner (Draw) or select a polygon'], ['Drag', 'pan, or move a corner of the selected polygon'],
        ['Wheel', 'zoom at the cursor'], ['Enter', 'close the polygon (or double-click)'], ['Esc', 'cancel the polygon'],
        ['Bksp', 'remove the last corner'], ['Del', 'remove the selected polygon'], ['Ctrl+Z', 'undo (Ctrl+Y: redo)'],
        ['Ctrl+S', 'save the region'],
      ].map(([k, v]) => [h('dt', null, h('kbd', null, k)), h('dd', null, v)])))));

  function card(title, body) {
    return h('section', { class: 'card' }, h('div', { class: 'card-header' }, h('h2', null, title)), h('div', { class: 'card-body' }, body));
  }

  // ---------------------------------------------------------------- view and painting
  const view = () => {
    const en = entry();
    if (!en.view) en.view = { s: 1, tx: 0, ty: 0, auto: true };
    return en.view;
  };

  function toImage(e) {
    const r = viewport.getBoundingClientRect();
    const v = view();
    return [(e.clientX - r.left - v.tx) / v.s, (e.clientY - r.top - v.ty) / v.s];
  }

  const clampPt = ([x, y]) => [Math.min(Math.max(x, 0), W), Math.min(Math.max(y, 0), H)];
  const screenDist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1]) * view().s;

  function fitScale() {
    const r = viewport.getBoundingClientRect();
    if (!W || !H || !r.width || !r.height) return 1;
    return Math.min(r.width / W, r.height / H);
  }

  function fit() {
    const r = viewport.getBoundingClientRect();
    const s = fitScale() * 0.96;
    entry().view = { s, tx: (r.width - W * s) / 2, ty: (r.height - H * s) / 2, auto: true };
    applyView();
  }

  function zoomAt(cx, cy, factor) {
    if (!ready) return;
    const v = view();
    const f = fitScale();
    const s2 = Math.min(Math.max(v.s * factor, f * 0.5), Math.max(MAX_ZOOM, f));
    v.tx = cx - (cx - v.tx) * (s2 / v.s);
    v.ty = cy - (cy - v.ty) * (s2 / v.s);
    v.s = s2;
    v.auto = false;
    applyView();
  }

  function zoomCenter(factor) {
    const r = viewport.getBoundingClientRect();
    zoomAt(r.width / 2, r.height / 2, factor);
  }

  function applyView() {
    const v = view();
    stage.style.transform = `translate(${v.tx}px, ${v.ty}px) scale(${v.s})`;
    paintPolys();
    paintHandles();
    paintDrawing();
  }

  function paintPolys() {
    mount(gPolys, entry().polys.map((p, i) => svg('path', { d: pathD(p), class: i === selected ? 'roi-poly selected' : 'roi-poly' })));
  }

  function paintHandles() {
    const poly = entry().polys[selected];
    const r = HANDLE_PX / view().s;
    mount(gHandles, poly ? poly.rings.flatMap((ring) => ring.map(([x, y]) => svg('circle', { cx: x, cy: y, r, class: 'roi-handle' }))) : []);
  }

  function paintDrawing() {
    if (!drawing || !drawing.length) { mount(gDraw); return; }
    const s = view().s;
    const pts = hover ? [...drawing, clampPt(hover)] : drawing;
    const snap = hover && drawing.length >= 3 && screenDist(clampPt(hover), drawing[0]) <= SNAP_PX;
    mount(gDraw,
      svg('polyline', { points: pts.map(([x, y]) => `${x},${y}`).join(' '), class: 'roi-draft' }),
      drawing.map(([x, y], i) => svg('circle', { cx: x, cy: y, r: (i === 0 ? HANDLE_PX + 2 : HANDLE_PX - 1) / s, class: i === 0 ? `roi-first${snap ? ' snap' : ''}` : 'roi-vertex' })));
  }

  function paintHint() {
    let text = '';
    if (ready) {
      if (drawing && drawing.length) {
        text = `${plural(drawing.length, 'corner')}. Click the first corner, double-click or press Enter to close; Backspace removes the last corner; Escape cancels.`;
      } else if (mode === 'draw') {
        text = 'Click to place the first corner. Drag to pan, scroll to zoom.';
      } else if (selected >= 0) {
        text = 'Drag a corner to move it. Delete removes the polygon; Escape deselects.';
      } else {
        text = entry().polys.length ? 'Click a polygon to select it. Drag to pan, scroll to zoom.' : 'No region yet: choose Draw (D) and outline the region of interest.';
      }
    }
    hint.textContent = text;
    selectBtn.classList.toggle('active', mode === 'select');
    drawBtn.classList.toggle('active', mode === 'draw');
    selectBtn.setAttribute('aria-pressed', String(mode === 'select'));
    drawBtn.setAttribute('aria-pressed', String(mode === 'draw'));
    viewport.dataset.mode = mode;
  }

  function paintSamples() {
    mount(sampleSel, samples().length ? samples().map((s) => h('option', { value: s.id, selected: s.id === store.sample },
      `${s.label || s.id}${s.label && s.label !== s.id ? ` (${s.id})` : ''}${isUnsaved(s.id) ? ' • unsaved' : savedRois.has(s.id) ? ' • region' : ''}`))
      : h('option', { value: '' }, 'No sections'));
    sampleSel.disabled = !samples().length;
  }

  function areaText(poly) {
    const a = polygonArea(poly);
    const sf = (meta && meta.scalefactors) || {};
    const scalef = Number(sf.tissue_hires_scalef);
    const mpp = Number((meta && meta.microns_per_pixel) || sf.microns_per_pixel);
    if (scalef > 0 && mpp > 0) return `${fmt.num((a / (scalef * scalef)) * mpp * mpp / 1e6, 3)} mm²`;
    return `${fmt.int(a)} px²`;
  }

  function paintSide() {
    paintSamples();
    const sid = store.sample;
    const en = sid ? entry() : null;
    const unsaved = sid ? isUnsaved(sid) : false;
    const savedHere = sid ? savedRois.has(sid) : false;
    const n = en ? en.polys.length : 0;
    const saveBtn = button({ label: saving ? 'Saving…' : 'Save region', icon: 'save', kind: 'primary', small: true, disabled: !unsaved || saving, onClick: saveRegion });
    const revertBtn = button({ label: 'Revert', icon: 'undo', small: true, disabled: !unsaved || saving, onClick: revert });
    revertBtn.title = 'Back to the saved region (Ctrl+Z restores your edits)';
    const delBtn = button({ label: 'Delete region', icon: 'trash', kind: 'ghost danger', small: true, disabled: !savedHere || saving, onClick: deleteRegion });
    mount(sideRegion,
      sid ? h('div', { class: 'row between' },
        h('span', null, savedHere ? 'Saved region' : 'No saved region'),
        unsaved ? h('span', { class: 'badge badge-unsaved', title: 'Edits not saved yet' }, 'Unsaved') : null) : h('p', { class: 'muted' }, 'No section.'),
      sid ? h('div', { class: 'row' }, saveBtn, revertBtn, delBtn) : null,
      n ? h('ol', { class: 'roi-polys' }, en.polys.map((p, i) => h('li', { class: i === selected ? 'active' : '' },
        h('button', {
          type: 'button', class: 'linklike', onClick: () => { selected = i; setMode('select'); paintPolys(); paintHandles(); paintSide(); viewport.focus({ preventScroll: true }); },
        }, `Polygon ${i + 1}`),
        h('span', { class: 'muted small' }, `${plural(p.rings[0].length, 'corner')}${p.rings.length > 1 ? `, ${plural(p.rings.length - 1, 'hole')}` : ''} · ${areaText(p)}`),
        iconButton('trash', `Remove polygon ${i + 1}`, () => removePolygon(i))))) : (sid ? h('p', { class: 'muted small' }, 'No polygon drawn for this section.') : null));
    paintAnalysis();
  }

  const roiCfg = () => {
    if (!cfg.roi || typeof cfg.roi !== 'object') cfg.roi = {};
    return cfg.roi;
  };
  const useSwitch = switchToggle({
    checked: !!roiCfg().use_drawn, label: 'Use the drawn region in the analysis',
    onChange: (v) => { roiCfg().use_drawn = v; onConfigChange(); paintAnalysis(); },
  });
  const analysisNotes = h('div', { class: 'stack' });
  mount(sideAnalysis, useSwitch,
    h('p', { class: 'muted small' }, 'When on, each section is restricted to the bins inside its drawn region, intersected with the '
      + 'transcriptional region of interest (Parameters → Region of interest). Sections without a drawn region are not restricted. '
      + 'This switch is a project setting (Save in the bar at the bottom of the page, or Ctrl+S); regions are saved with “Save region”.'),
    analysisNotes);

  function paintAnalysis() {
    const roi = roiCfg();
    useSwitch.querySelector('input').checked = !!roi.use_drawn;
    const withRegion = samples().filter((s) => savedRois.has(s.id));
    const pending = samples().filter((s) => isUnsaved(s.id));
    mount(analysisNotes,
      roi.use_drawn && roi.mask_frame && roi.mask_frame !== 'hires' ? banner('warning', h('span', null,
        `Regions drawn here are in hires image pixels, but “Image the regions were drawn on” is ${roi.mask_frame}. `,
        h('button', { type: 'button', class: 'linklike', onClick: () => { roi.mask_frame = 'hires'; onConfigChange(); paintAnalysis(); } }, 'Set it to hires'))) : null,
      roi.use_drawn && samples().length && !withRegion.length ? banner('info', 'No section has a saved region yet: the analysis is not restricted.') : null,
      pending.length ? banner('warning', `Unsaved region${pending.length === 1 ? '' : 's'}: ${pending.map((s) => s.id).join(', ')}. Runs use the saved regions only.`) : null,
      samples().length ? h('ul', { class: 'roi-summary' }, samples().map((s) => h('li', null,
        h('span', { class: 'swatch', style: { background: /^#[0-9A-Fa-f]{6}$/.test(String(s.color)) ? s.color : 'transparent' } }),
        h('span', { class: 'mono' }, s.id),
        h('span', { class: 'muted' }, savedRois.has(s.id) ? 'region saved' : 'no region (not restricted)')))) : null);
  }

  // ---------------------------------------------------------------- edits
  function pushUndo(snapshot) {
    const en = entry();
    en.undo.push(snapshot);
    if (en.undo.length > MAX_UNDO) en.undo.shift();
    en.redo = [];
  }

  function changed() {
    paintPolys();
    paintHandles();
    paintDrawing();
    paintSide();
    paintHint();
    onRoisChange();
  }

  function commit(mutate) {
    pushUndo(JSON.stringify(entry().polys));
    mutate(entry());
    changed();
  }

  function undo() {
    const en = entry();
    if (!en.undo.length) return;
    en.redo.push(JSON.stringify(en.polys));
    en.polys = JSON.parse(en.undo.pop());
    selected = -1;
    changed();
  }

  function redo() {
    const en = entry();
    if (!en.redo.length) return;
    en.undo.push(JSON.stringify(en.polys));
    en.polys = JSON.parse(en.redo.pop());
    selected = -1;
    changed();
  }

  function closeDrawing() {
    if (!drawing) return;
    if (drawing.length < 3) {
      toast('A region needs at least three corners.', 'warning', { timeout: 2500 });
      return;
    }
    const ring = drawing;
    drawing = null;
    commit((en) => {
      en.polys.push({ rings: [ring] });
      selected = en.polys.length - 1;
    });
  }

  function removePolygon(i) {
    if (i < 0 || i >= entry().polys.length) return;
    selected = -1;
    commit((en) => { en.polys.splice(i, 1); });
  }

  function popVertex() {
    if (!drawing) return;
    drawing.pop();
    if (!drawing.length) drawing = null;
    paintDrawing();
    paintHint();
  }

  function setMode(m) {
    if (m === mode) { paintHint(); return; }
    if (mode === 'draw' && drawing) {
      if (drawing.length >= 3) closeDrawing();
      else drawing = null;
    }
    mode = m;
    if (m === 'draw') selected = -1;
    paintPolys();
    paintHandles();
    paintDrawing();
    paintSide();
    paintHint();
  }

  async function saveRegion() {
    const sid = store.sample;
    if (!sid || saving || !isUnsaved(sid)) return;
    const en = entry();
    const sent = JSON.stringify(en.polys);
    const polys = JSON.parse(sent);
    saving = true;
    paintSide();
    try {
      if (polys.length) {
        const res = await api(roiUrl(sid), { method: 'PUT', body: polygonsToGeoJSON(polys, (meta && meta.frame) || 'hires') });
        savedRois.add(sid);
        toast(`Region of ${sid} saved (${plural(res && Number.isFinite(res.n_polygons) ? res.n_polygons : polys.length, 'polygon')})`, 'success');
      } else {
        try { await api(roiUrl(sid), { method: 'DELETE' }); } catch (e) { if (e.status !== 404) throw e; }
        savedRois.delete(sid);
        toast(`Region of ${sid} removed`, 'success');
      }
      en.saved = sent;
    } catch (e) {
      toast(`The region could not be saved: ${e.message}`, 'error');
    } finally {
      saving = false;
      if (!destroyed) paintSide();
      onRoisChange();
    }
  }

  function revert() {
    const en = entry();
    if (!isUnsaved(store.sample)) return;
    drawing = null;
    selected = -1;
    commit((e) => { e.polys = JSON.parse(en.saved); });
  }

  async function deleteRegion() {
    const sid = store.sample;
    if (!sid || !savedRois.has(sid)) return;
    const ok = await confirmDialog({
      title: `Delete the region of ${sid}?`,
      message: 'The saved region is removed from the project; this section is then not restricted by a drawn region. Runs already started keep their copy.',
      confirmLabel: 'Delete region', danger: true,
    });
    if (!ok || destroyed) return;
    try {
      await api(roiUrl(sid), { method: 'DELETE' });
    } catch (e) {
      if (e.status !== 404) { toast(`The region could not be deleted: ${e.message}`, 'error'); return; }
    }
    savedRois.delete(sid);
    const en = entryOf(sid);
    en.saved = '[]';
    if (store.sample === sid && !destroyed) {
      selected = -1;
      drawing = null;
      commit((e) => { e.polys = []; });
    } else {
      en.polys = [];
      onRoisChange();
    }
    toast(`Region of ${sid} deleted`, 'success');
  }

  // ---------------------------------------------------------------- pointer and keyboard
  function nearVertex(p) {
    const poly = entry().polys[selected];
    if (!poly) return null;
    const tol = (HANDLE_PX + 3) / view().s;
    let best = null;
    poly.rings.forEach((ring, ri) => ring.forEach(([x, y], vi) => {
      const d = Math.hypot(x - p[0], y - p[1]);
      if (d <= tol && (!best || d < best.d)) best = { ring: ri, index: vi, d };
    }));
    return best;
  }

  function polygonAt(p) {
    const polys = entry().polys;
    for (let i = polys.length - 1; i >= 0; i -= 1) if (pointInPolygon(p, polys[i])) return i;
    return -1;
  }

  function handleClick(p, e) {
    const now = performance.now();
    const isDouble = !!lastClick && now - lastClick.t < DOUBLE_MS && Math.hypot(e.clientX - lastClick.x, e.clientY - lastClick.y) <= CLICK_SLOP_PX * 2;
    lastClick = isDouble ? null : { t: now, x: e.clientX, y: e.clientY };
    if (mode === 'draw') {
      const pt = clampPt(p);
      if (isDouble) { if (drawing && drawing.length >= 3) closeDrawing(); return; }
      if (drawing && drawing.length >= 3 && screenDist(pt, drawing[0]) <= SNAP_PX) { closeDrawing(); return; }
      if (!drawing) drawing = [];
      drawing.push(pt);
      paintDrawing();
      paintHint();
      return;
    }
    const i = polygonAt(p);
    if (i !== selected) {
      selected = i;
      paintPolys();
      paintHandles();
      paintSide();
      paintHint();
    }
  }

  viewport.addEventListener('pointerdown', (e) => {
    if (!ready || (e.button !== 0 && e.button !== 1)) return;
    e.preventDefault();
    viewport.focus({ preventScroll: true });
    const p = toImage(e);
    if (e.button === 0 && mode === 'select') {
      const hit = nearVertex(p);
      if (hit) {
        drag = { kind: 'vertex', ring: hit.ring, index: hit.index, before: JSON.stringify(entry().polys), moved: false };
        viewport.setPointerCapture(e.pointerId);
        return;
      }
    }
    const v = view();
    drag = { kind: 'pan', x0: e.clientX, y0: e.clientY, tx0: v.tx, ty0: v.ty, moved: false, button: e.button };
    viewport.setPointerCapture(e.pointerId);
  });

  viewport.addEventListener('pointermove', (e) => {
    if (!ready) return;
    const p = toImage(e);
    hover = p;
    coords.textContent = p[0] >= 0 && p[1] >= 0 && p[0] <= W && p[1] <= H ? `x ${fmt.int(p[0])} · y ${fmt.int(p[1])} px` : '';
    if (drag && drag.kind === 'pan') {
      const dx = e.clientX - drag.x0;
      const dy = e.clientY - drag.y0;
      if (!drag.moved && Math.hypot(dx, dy) < CLICK_SLOP_PX) return;
      drag.moved = true;
      viewport.classList.add('panning');
      const v = view();
      v.tx = drag.tx0 + dx;
      v.ty = drag.ty0 + dy;
      v.auto = false;
      applyView();
      return;
    }
    if (drag && drag.kind === 'vertex') {
      const poly = entry().polys[selected];
      if (!poly) return;
      poly.rings[drag.ring][drag.index] = clampPt(p);
      drag.moved = true;
      paintPolys();
      paintHandles();
      return;
    }
    if (mode === 'select') viewport.classList.toggle('over-handle', !!nearVertex(p));
    if (drawing) paintDrawing();
  });

  const endDrag = (e) => {
    if (!drag) return;
    const d = drag;
    drag = null;
    viewport.classList.remove('panning');
    if (viewport.hasPointerCapture(e.pointerId)) viewport.releasePointerCapture(e.pointerId);
    if (d.kind === 'vertex') {
      if (d.moved) { pushUndo(d.before); changed(); }
      return;
    }
    if (e.type === 'pointerup' && !d.moved && d.button === 0) handleClick(toImage(e), e);
  };
  viewport.addEventListener('pointerup', endDrag);
  viewport.addEventListener('pointercancel', endDrag);
  viewport.addEventListener('pointerleave', () => {
    if (drag) return;
    hover = null;
    coords.textContent = '';
    paintDrawing();
  });
  viewport.addEventListener('wheel', (e) => {
    if (!ready) return;
    e.preventDefault();
    const r = viewport.getBoundingClientRect();
    const unit = e.deltaMode === 1 ? 16 : e.deltaMode === 2 ? r.height : 1;
    zoomAt(e.clientX - r.left, e.clientY - r.top, Math.exp(-e.deltaY * unit * 0.0015));
  }, { passive: false });
  viewport.addEventListener('contextmenu', (e) => { if (drawing) e.preventDefault(); });

  viewport.addEventListener('keydown', (e) => {
    const mod = e.ctrlKey || e.metaKey;
    const key = e.key;
    if (mod && (key === 's' || key === 'S')) { e.preventDefault(); e.stopPropagation(); saveRegion(); return; }
    if (!ready) return;
    if (mod && (key === 'z' || key === 'Z')) {
      e.preventDefault();
      if (e.shiftKey) redo(); else if (drawing) popVertex(); else undo();
      return;
    }
    if (mod && (key === 'y' || key === 'Y')) { e.preventDefault(); redo(); return; }
    if (mod || e.altKey) return;
    const v = view();
    let handled = true;
    switch (key) {
      case 'Escape':
        if (drawing) { drawing = null; paintDrawing(); paintHint(); } else if (mode === 'draw') setMode('select');
        else if (selected >= 0) { selected = -1; paintPolys(); paintHandles(); paintSide(); paintHint(); } else handled = false;
        break;
      case 'Enter': if (drawing) closeDrawing(); else handled = false; break;
      case 'Backspace': if (drawing) popVertex(); else if (selected >= 0) removePolygon(selected); else handled = false; break;
      case 'Delete': if (selected >= 0) removePolygon(selected); else handled = false; break;
      case 'd': case 'D': setMode('draw'); break;
      case 'v': case 'V': setMode('select'); break;
      case '+': case '=': zoomCenter(1.4); break;
      case '-': case '_': zoomCenter(1 / 1.4); break;
      case '0': fit(); break;
      case 'ArrowLeft': v.tx += 60; v.auto = false; applyView(); break;
      case 'ArrowRight': v.tx -= 60; v.auto = false; applyView(); break;
      case 'ArrowUp': v.ty += 60; v.auto = false; applyView(); break;
      case 'ArrowDown': v.ty -= 60; v.auto = false; applyView(); break;
      default: handled = false;
    }
    if (handled) { e.preventDefault(); e.stopPropagation(); }
  });

  const resizer = new ResizeObserver(() => { if (ready && view().auto) fit(); });
  resizer.observe(viewport);

  // ---------------------------------------------------------------- loading
  function showMessage(node) {
    ready = false;
    mount(msg, node);
    msg.classList.remove('hidden');
    stage.classList.add('hidden');
    coords.textContent = '';
    paintHint();
  }

  function hideMessage() {
    mount(msg);
    msg.classList.add('hidden');
    stage.classList.remove('hidden');
  }

  // savedRois comes from the project ("rois"), so sections without a region need no request.
  async function loadRoi(sid, en) {
    en.polys = [];
    if (savedRois.has(sid)) {
      try {
        en.polys = polygonsFromGeoJSON(await api(roiUrl(sid)));
      } catch (e) {
        if (e.status !== 404) throw e;
        savedRois.delete(sid);
      }
    }
    en.saved = JSON.stringify(en.polys);
    en.loaded = true;
  }

  function loadImage(sid) {
    const src = apiUrl(`${sampleBase(sid)}/image`, { v: getVersion() });
    return new Promise((resolve, reject) => {
      if (img.getAttribute('src') === src && img.complete && img.naturalWidth) { resolve(); return; }
      img.onload = () => resolve();
      img.onerror = () => reject(new ApiError('The H&E image could not be loaded', 404));
      img.setAttribute('src', src);
    });
  }

  function binPx() {
    const sf = (meta && meta.scalefactors) || {};
    const scalef = Number(sf.tissue_hires_scalef);
    const mpp = Number((meta && meta.microns_per_pixel) || sf.microns_per_pixel);
    const m = /^(\d{3})um$/.exec(String(cfg.bin_size));
    if (m && scalef > 0 && mpp > 0) return (Number(m[1]) / mpp) * scalef;
    if (scalef > 0 && Number(sf.spot_diameter_fullres) > 0) return Number(sf.spot_diameter_fullres) * scalef;
    return 2;
  }

  function drawBins(data) {
    const xs = Array.isArray(data.x) ? data.x : [];
    const ys = Array.isArray(data.y) ? data.y : [];
    const k = Math.min(1, BIN_CANVAS_MAX / Math.max(W, H, 1));
    canvas.width = Math.max(1, Math.round(W * k));
    canvas.height = Math.max(1, Math.round(H * k));
    canvas.style.width = `${W}px`;
    canvas.style.height = `${H}px`;
    const ctx = canvas.getContext('2d');
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = getComputedStyle(viewport).getPropertyValue('--roi-bins').trim() || 'rgba(0, 150, 210, 0.55)';
    const d = Math.max(1, binPx() * k * 0.85);
    const n = Math.min(xs.length, ys.length);
    for (let i = 0; i < n; i += 1) ctx.fillRect(xs[i] * k - d / 2, ys[i] * k - d / 2, d, d);
    canvas.classList.remove('hidden');
    canvas.title = data.n_total > n ? `${fmt.int(n)} of ${fmt.int(data.n_total)} bins shown` : '';
  }

  async function loadBins() {
    const sid = store.sample;
    binsInput.checked = !!store.showBins;
    if (!store.showBins || !ready || !sid) { canvas.classList.add('hidden'); return; }
    const key = `${sid}@${getVersion()}`;
    try {
      if (!store.bins[key]) {
        binsSwitch.classList.add('busy');
        store.bins[key] = await api(`${sampleBase(sid)}/bins-preview`);
      }
      if (destroyed || store.sample !== sid || !store.showBins || !ready) return;
      drawBins(store.bins[key]);
    } catch (e) {
      if (destroyed) return;
      toast(e.status === 404 ? `No bin positions for ${sid}.` : `The bins could not be loaded: ${e.message}`, 'warning');
      store.showBins = false;
      binsInput.checked = false;
      canvas.classList.add('hidden');
    } finally {
      binsSwitch.classList.remove('busy');
    }
  }

  async function loadSample() {
    const sid = store.sample;
    const my = ++loadSeq;
    drawing = null;
    selected = -1;
    drag = null;
    meta = null;
    blocked = false;
    shownNote = '';
    canvas.classList.add('hidden');
    mount(note);
    note.classList.add('hidden');
    paintSide();
    if (!sid) {
      showMessage(emptyState({ icon: 'image', title: 'No sections', text: 'Add samples in the Data tab first.' }));
      return;
    }
    const st = sampleStatus(sid);
    if (!st.ok) {
      blocked = true;
      showMessage(emptyState({
        icon: 'image', title: st.title, text: st.text,
        action: st.canSave ? button({ label: 'Save the project', icon: 'save', kind: 'primary', onClick: async () => { if (await onSaveProject() && !destroyed) loadSample(); } }) : null,
      }));
      return;
    }
    showMessage(loading('Loading the H&E image…'));
    const en = entry();
    try {
      const metaKey = `${sid}@${getVersion()}`;
      const [m] = await Promise.all([
        store.meta[metaKey] || api(`${sampleBase(sid)}/image-meta`),
        en.loaded ? null : loadRoi(sid, en),
      ]);
      if (my !== loadSeq || destroyed) return;
      store.meta[metaKey] = m;
      meta = m;
      await loadImage(sid);
      if (my !== loadSeq || destroyed) return;
      W = Number(m.width) || img.naturalWidth;
      H = Number(m.height) || img.naturalHeight;
      stage.style.width = `${W}px`;
      stage.style.height = `${H}px`;
      img.width = W;
      img.height = H;
      overlay.setAttribute('viewBox', `0 0 ${W} ${H}`);
      overlay.setAttribute('width', W);
      overlay.setAttribute('height', H);
      hideMessage();
      ready = true;
      shownNote = st.note || '';
      if (shownNote) {
        mount(note, icon('alert', { size: 14 }), h('span', null, shownNote));
        note.classList.remove('hidden');
      }
      if (view().auto) fit(); else applyView();
      paintSide();
      paintHint();
      if (store.showBins) loadBins();
    } catch (err) {
      if (my !== loadSeq || destroyed) return;
      if (err.status === 404) {
        showMessage(emptyState({
          icon: 'image', title: 'No H&E image for this section',
          text: `${err.message}. Regions are drawn on spatial/tissue_hires_image.png of the section's Space Ranger output; check its folder in the Data tab.`,
        }));
      } else {
        showMessage(errorPanel(err, { retry: () => loadSample() }));
      }
      paintSide();
    }
  }

  paintHint();
  loadSample();

  return {
    el,
    /** Reload when the section's status changed (e.g. the project was saved): blocked -> loadable, or the note changed. */
    refresh() {
      if (destroyed || !store.sample) return;
      const st = sampleStatus(store.sample);
      if ((blocked && st.ok) || (ready && (st.note || '') !== shownNote)) loadSample();
    },
    destroy() {
      destroyed = true;
      loadSeq += 1;
      resizer.disconnect();
      img.onload = null;
      img.onerror = null;
    },
  };
}
