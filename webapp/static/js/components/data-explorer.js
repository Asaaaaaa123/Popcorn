// data-explorer.js — the Explore tab of a project: the raw Space Ranger data of every section, read with the current
// settings (saved or not, like POST /validate) before any run. A table of the sections (bins, genes, UMIs, medians per
// bin, mitochondrial share, bins passing the configured quality control), the spatial viewer (one, two or four panels)
// on a projectSource over the explore endpoints, and the genes of the section shown in the highlighted panel, ranked by
// UMIs, each one a click away from the map.
// Endpoints: POST /api/projects/{pid}/explore/samples | /{sid}/bins | /{sid}/genes | /{sid}/gene, each with the draft
// config in the body.
//
// State (owned by the caller and kept between visits; exploreState() makes it, a new one when the data settings change):
//   {key, config (the settings the requests send), samples (POST explore/samples) | null, binSize,
//    sections: Map sid -> {status 'idle' | 'loading' | 'ok' | 'error', summary, notes, image, error},
//    viewer (the spatialViewer store), target (the panel the summary's and the gene table's actions apply to),
//    genes {query, limit}}
import { h, api, mount, icon, button, card, loading, emptyState, errorPanel, fmt, debounce, clone } from '../core.js';
import { spatialViewer } from './spatial-viewer.js';
import { binLabel } from './sample-editor.js';

const KEEP_AHEAD = 2; // bins payloads read for the table and kept for the viewer until it takes them
const GENE_LIMITS = [50, 200, 500];
const GENE_OVERLAY = 'gene:'; // overlay key of "detected" layers in the viewer store: 'gene:Ctsk' (spatial-viewer.js)
const HEX_RE = /^#[0-9A-Fa-f]{6}$/;

const enc = encodeURIComponent;
const aborted = () => new DOMException('The request was aborted.', 'AbortError');

/** '016um' -> 16; NaN for anything else (the viewer then estimates the bin size from the positions). */
function binUmOf(bin) {
  const m = /^(\d{3})um$/.exec(String(bin || ''));
  return m ? Number(m[1]) : NaN;
}

/** p, rejected with an AbortError as soon as `signal` aborts; p itself goes on (someone else may wait for it). */
function untilAborted(p, signal) {
  if (!signal) return p;
  return new Promise((resolve, reject) => {
    const onAbort = () => reject(aborted());
    if (signal.aborted) { onAbort(); return; }
    signal.addEventListener('abort', onAbort, { once: true });
    p.then((v) => { signal.removeEventListener('abort', onAbort); resolve(v); },
      (e) => { signal.removeEventListener('abort', onAbort); reject(e); });
  });
}

/**
 * A viewer store with the choices of `old` (layout, linking, shared color scales, point size, H&E, what each panel
 * shows) and no data.
 */
function carryViewer(old) {
  if (!old || typeof old !== 'object') return {};
  const keep = { layout: old.layout, link: old.link, sameScale: old.sameScale, pointSize: old.pointSize };
  if (old.image && typeof old.image === 'object') keep.image = { on: !!old.image.on, opacity: old.image.opacity };
  if (Array.isArray(old.panels)) {
    keep.panels = old.panels.map((p) => (p && typeof p === 'object' ? {
      sample: p.sample, field: p.field, genes: Array.isArray(p.genes) ? p.genes.slice() : [],
      overlays: Array.isArray(p.overlays) ? p.overlays.slice() : [], dim: p.dim, hidden: {}, views: {},
    } : null));
  }
  return keep;
}

/**
 * The explorer state for the settings `cfg` whose data settings are summed up by `key`: `previous` itself when its
 * key is the same, else a new state (nothing read yet) that keeps the viewer choices and the gene filter of `previous`.
 */
export function exploreState(cfg, key, previous = null) {
  if (previous && previous.key === key) return previous;
  return {
    key,
    config: clone(cfg),
    samples: null,
    binSize: null,
    sections: new Map(),
    viewer: carryViewer(previous && previous.viewer),
    target: previous && Number.isInteger(previous.target) ? previous.target : 0,
    genes: { query: previous && previous.genes ? previous.genes.query || '' : '', limit: GENE_LIMITS[0] },
  };
}

/**
 * dataExplorer({projectId, name, state, getConfig, focus, onOpenData}) -> {el, destroy()}
 * state: exploreState(); getConfig() -> the current draft (labels and colors are read from it); focus: a section to
 * show first (the Data tab's Explore action); onOpenData(): switch to the Data tab.
 */
export function dataExplorer({ projectId, name = '', state: st, getConfig, focus = null, onOpenData }) {
  let destroyed = false;
  let gen = 0; // bumped by reload() and destroy(): answers of an earlier generation are dropped
  let ctrl = new AbortController();
  let viewer = null;
  let observer = null;
  let payloads = new Map(); // sid -> promise of its bins payload (this generation), until the viewer takes it
  let ahead = []; // sids whose resolved payload is kept for the viewer, oldest first
  let rows = new Map(); // sid -> <tr> of the section table
  let geneSid = null;
  let geneRes = null;
  let geneErr = null;
  let geneBusy = false;
  let geneCtrl = null;
  let geneSeq = 0;
  let geneRows = []; // [{gene, tr, show, toggle}] of the gene table painted
  let genesShownFor = ''; // panel count, highlighted panel and its overlays when the gene rows were painted
  let titlesFor = 1; // panel count the section names' titles were written for

  const url = (path) => `/api/projects/${enc(projectId)}/explore${path}`;
  const body = (extra) => ({ config: st.config, ...extra });
  const draftSamples = () => {
    const cfg = getConfig() || {};
    return Array.isArray(cfg.samples) ? cfg.samples.filter((s) => s && typeof s === 'object') : [];
  };
  const draftSample = (sid) => draftSamples().find((s) => s.id === sid) || null;
  const labelOf = (s) => {
    const ds = draftSample(s.id);
    return String((ds && ds.label) || s.label || s.id);
  };
  const colorOf = (s) => {
    const ds = draftSample(s.id);
    const c = ds && HEX_RE.test(String(ds.color)) ? ds.color : s.color;
    return HEX_RE.test(String(c)) ? c : 'transparent';
  };
  const readable = () => (st.samples || []).filter((s) => s.ok);
  const sampleById = (sid) => (st.samples || []).find((s) => s.id === sid) || null;
  const sectionText = (sid) => {
    const s = sampleById(sid);
    const label = s ? labelOf(s) : sid;
    return label !== sid ? `${label} (${sid})` : sid;
  };

  function sectionState(sid) {
    let sec = st.sections.get(sid);
    if (!sec) {
      sec = { status: 'idle', summary: null, notes: [], image: null, error: null };
      st.sections.set(sid, sec);
    }
    return sec;
  }

  // ---------------------------------------------------------------- DOM
  const tableHost = h('div', { class: 'dx-sections' });
  const qcLine = h('p', { class: 'dx-qc' });
  const notesHost = h('div');
  const reloadBtn = button({
    label: 'Read again', icon: 'refresh', small: true, onClick: () => reload(),
    title: 'Read the sections again (after the files on disk changed)',
  });
  const summaryCard = card({
    title: 'Sections',
    subtitle: 'The Space Ranger outputs as the analysis reads them, with the current settings (saved or not), before quality control and any run.',
    actions: reloadBtn,
    body: h('div', { class: 'stack' }, tableHost, qcLine, notesHost),
  });
  const mapHost = h('div', { class: 'dx-map' });
  const geneTitle = h('span', null, 'Genes');
  const geneSub = h('span');
  const geneFilterId = `dx-genes-${projectId}`;
  const geneFilter = h('input', {
    class: 'input', type: 'search', id: geneFilterId, value: st.genes.query, spellcheck: false, autocomplete: 'off',
    placeholder: 'Filter by name, e.g. Ctsk or mt-',
  });
  const geneCount = h('span', { class: 'muted small' });
  const geneBody = h('div', { class: 'dx-gene-body' });
  const geneCard = card({
    title: geneTitle, subtitle: geneSub, cls: 'dx-genes hidden',
    body: h('div', { class: 'stack' },
      h('div', { class: 'dx-gene-bar' }, h('label', { class: 'sr-only', for: geneFilterId }, 'Filter the genes by name'), geneFilter, geneCount),
      geneBody),
  });
  geneCard.tabIndex = -1; // the target of the jump below
  // stacked or several maps push the gene table far down: a way to it from the top of the map (css/project.css)
  const geneJump = h('div', { class: 'dx-jump hidden' }, button({
    label: 'Gene table', icon: 'chevron-down', small: true, kind: 'ghost', onClick: () => jumpToGenes(),
    title: 'Go to the genes of the section on the map, below it',
  }));
  const el = h('div', { class: 'dx' }, summaryCard, geneJump, mapHost, geneCard);

  // ---------------------------------------------------------------- reading the sections
  function forget(sid, p) {
    if (payloads.get(sid) === p) payloads.delete(sid);
    ahead = ahead.filter((x) => x !== sid);
  }

  /** The bins payload of a section: one request per section, shared by the table, the viewer and the image. */
  function binsOf(sid, fresh = false) {
    if (!fresh && payloads.has(sid)) return payloads.get(sid);
    const g = gen;
    const sec = sectionState(sid);
    if (!sec.summary) {
      sec.status = 'loading';
      sec.error = null;
      paintRow(sid);
    }
    const p = api(url(`/${enc(sid)}/bins`), { method: 'POST', body: body(), signal: ctrl.signal });
    payloads.set(sid, p);
    p.then((payload) => {
      if (g !== gen) return;
      sec.status = 'ok';
      sec.error = null;
      sec.summary = payload && payload.summary && typeof payload.summary === 'object' ? payload.summary : {};
      sec.notes = payload && Array.isArray(payload.notes) ? payload.notes.map(String) : [];
      sec.image = payload && payload.image && typeof payload.image === 'object' ? payload.image : null;
      // read ahead for the table: kept for the viewer, the oldest dropped (the server keeps them too)
      if (payloads.get(sid) === p) {
        ahead = [...ahead.filter((x) => x !== sid), sid];
        while (ahead.length > KEEP_AHEAD) payloads.delete(ahead.shift());
      }
      paintRow(sid);
      paintFoot();
    }, (err) => {
      if (payloads.get(sid) === p) payloads.delete(sid);
      if (g !== gen) return;
      if (err && err.name === 'AbortError') {
        if (sec.status === 'loading') sec.status = 'idle';
        return;
      }
      if (!sec.summary) {
        sec.status = 'error';
        sec.error = (err && err.message) || String(err);
      }
      paintRow(sid);
    });
    return p;
  }

  /** Read every readable section in turn (the ones on the map first), so the table fills up. */
  async function readAll(first) {
    const g = gen;
    const ids = readable().map((s) => s.id);
    const order = [...new Set([...first.filter((sid) => ids.includes(sid)), ...ids])];
    for (const sid of order) {
      if (g !== gen) return;
      const sec = sectionState(sid);
      if (sec.summary || sec.status === 'error') continue;
      try { await binsOf(sid); } catch { /* the row shows it */ }
    }
  }

  const source = {
    name: `${name || projectId} raw data`,
    emptyText: 'None of the sections can be read with the current data settings; the table above says why.',
    samples: async () => readable().map((s) => ({ id: s.id, label: labelOf(s), n: s.n_bins })),
    async bins(sid, signal) {
      const p = binsOf(sid);
      const payload = await untilAborted(p, signal);
      forget(sid, p); // the viewer keeps its own copy
      return payload;
    },
    genes: {
      async search(sid, query, signal) {
        const res = await api(url(`/${enc(sid)}/genes`), { method: 'POST', body: body({ query, limit: GENE_LIMITS[0] }), signal });
        return res && Array.isArray(res.genes) ? res.genes : [];
      },
      values: (sid, genes, signal) => api(url(`/${enc(sid)}/gene`), { method: 'POST', body: body({ genes }), signal }),
    },
    async image(sid, signal, { retry = false } = {}) {
      const sec = sectionState(sid);
      if (sec.summary && !retry) return sec.image;
      // the Retry of an image that failed to load: ask the server again, since an image link only lasts while the
      // server that issued it runs
      const p = binsOf(sid, !!sec.summary);
      const payload = await untilAborted(p, signal);
      forget(sid, p);
      return payload && payload.image && typeof payload.image === 'object' ? payload.image : null;
    },
  };

  // ---------------------------------------------------------------- section table
  /** Section name (a button that shows it on the map when it can be read), id and the matrix it comes from. */
  function sectionCell(s) {
    const sec = st.sections.get(s.id);
    const matrix = sec && sec.summary ? sec.summary.matrix : '';
    const label = labelOf(s);
    const nameEl = s.ok
      ? h('button', {
        type: 'button', class: 'linklike dx-section-label', title: `Show ${sectionText(s.id)} on ${where()}`,
        onClick: () => showSection(s.id),
      }, label)
      : h('span', { class: 'dx-section-label' }, label);
    return h('td', { title: matrix ? `Matrix: ${matrix}` : null },
      h('div', { class: 'dx-section' },
        h('span', { class: 'swatch', style: { background: colorOf(s) } }),
        h('span', { class: 'dx-section-text' }, nameEl, label !== s.id ? h('span', { class: 'mono small muted' }, s.id) : null)));
  }

  function sectionRow(s) {
    const head = sectionCell(s);
    if (!s.ok) {
      return h('tr', { class: 'row-bad' }, head,
        h('td', { colSpan: 8 }, h('div', { class: 'dx-row-msg bad' }, icon('alert', { size: 14 }), h('span', null, s.message || 'This section cannot be read.'))));
    }
    const sec = sectionState(s.id);
    const sum = sec.summary;
    if (!sum) {
      const msg = sec.status === 'error'
        ? h('div', { class: 'dx-row-msg bad' }, icon('alert', { size: 14 }), h('span', null, sec.error || 'The section could not be read.'),
          button({ label: 'Retry', icon: 'refresh', small: true, kind: 'ghost', onClick: () => { binsOf(s.id, true).catch(() => {}); } }))
        : h('div', { class: 'dx-row-msg' }, sec.status === 'loading' ? h('span', { class: 'spinner' }) : icon('clock', { size: 14 }),
          h('span', null, sec.status === 'loading' ? 'Reading the counts…' : 'Waiting for its turn…'));
      return h('tr', null, head,
        h('td', { class: 'num muted', title: 'Bins under tissue in the positions file' }, fmt.int(s.n_bins)),
        h('td', { colSpan: 7 }, msg));
    }
    const missing = Number(sum.n_without_position) || 0;
    return h('tr', null, head,
      h('td', { class: 'num', title: `${fmt.int(sum.n_barcodes)} barcodes in the matrix${missing ? `, ${fmt.int(missing)} of them without a position under tissue (left out)` : ''}` }, fmt.int(sum.n_bins)),
      h('td', { class: 'num' }, fmt.int(sum.n_features)),
      h('td', { class: 'num', title: `Mean ${fmt.num(sum.mean_umi, 1)} UMIs per bin` }, fmt.int(sum.total_umi)),
      h('td', { class: 'num' }, fmt.num(sum.median_umi, 1)),
      h('td', { class: 'num' }, fmt.num(sum.median_genes, 1)),
      h('td', { class: 'num' }, fmt.pct(sum.median_pct_mt, 2)),
      h('td', { class: 'num' }, h('div', null, fmt.int(sum.n_pass_qc)), h('div', { class: 'small muted' }, fmt.pct(sum.pct_pass_qc, 1))),
      h('td', { class: 'num' }, fmt.num(sum.microns_per_pixel, 4)));
  }

  function paintRow(sid) {
    if (destroyed) return;
    const old = rows.get(sid);
    const s = sampleById(sid);
    if (!old || !s) return;
    const tr = sectionRow(s);
    tr.classList.toggle('active', sid === geneSid);
    rows.set(sid, tr);
    old.replaceWith(tr);
  }

  function paintTable() {
    const samples = st.samples || [];
    rows = new Map();
    if (!samples.length) {
      mount(tableHost, emptyState({
        icon: 'folder', title: 'No sections yet', text: 'Add the Space Ranger outputs in the Data tab, then come back here to look at them.',
        action: onOpenData ? button({ label: 'Go to the Data tab', icon: 'folder', kind: 'primary', onClick: () => onOpenData() }) : null,
      }));
      return;
    }
    const th = (text, title, num = true) => h('th', { class: num ? 'num' : null, title }, text);
    const trs = samples.map((s) => {
      const tr = sectionRow(s);
      tr.classList.toggle('active', s.id === geneSid);
      rows.set(s.id, tr);
      return tr;
    });
    mount(tableHost, h('div', { class: 'table-wrap' }, h('table', { class: 'table compact dx-table' },
      h('thead', null, h('tr', null,
        th('Section', 'Click a name to show the section on the map', false),
        th('Bins', 'Bins under tissue with a position: the bins the analysis loads'),
        th('Genes', 'Genes (Gene Expression features) in the matrix'),
        th('UMIs', 'All UMIs of the section (transcripts captured)'),
        th('Median UMIs per bin'),
        th('Median genes per bin', 'Median number of genes detected in a bin'),
        th('Median % mitochondrial', 'Median share of the UMIs of a bin on mitochondrial genes (names starting with mt-)'),
        th('Bins passing QC', 'Bins that pass the quality control of the current settings (see below)'),
        th('µm per pixel', 'Microns per full-resolution pixel'))),
      h('tbody', null, trs))));
  }

  function paintFoot() {
    if (destroyed) return;
    const ok = readable().map((s) => st.sections.get(s.id)).find((sec) => sec && sec.summary && sec.summary.qc);
    const cfgQc = st.config && st.config.qc && typeof st.config.qc === 'object' ? st.config.qc : {};
    const qc = ok ? ok.summary.qc : cfgQc;
    const bin = st.binSize || (st.config && st.config.bin_size);
    mount(qcLine,
      bin ? `Analysis bin size ${binLabel(bin)} (square_${bin}). ` : null,
      'Quality control (Parameters → Quality control): a bin passes with more than ', fmt.num(qc.min_genes), ' genes, more than ',
      fmt.num(qc.min_counts), ' UMIs and less than ', fmt.num(qc.max_pct_mt), '% of its UMIs on mitochondrial genes.');
    const notes = [];
    for (const s of st.samples || []) {
      const sec = st.sections.get(s.id);
      if (!sec || !sec.summary) continue;
      const missing = Number(sec.summary.n_without_position) || 0;
      if (missing) notes.push(`${sectionText(s.id)}: ${fmt.int(missing)} barcodes of the matrix have no position under tissue and are left out, as in the analysis.`);
      for (const n of sec.notes || []) notes.push(`${sectionText(s.id)}: ${n}`);
    }
    mount(notesHost, notes.length ? h('ul', { class: 'dx-notes' }, notes.map((t) => h('li', null, icon('info', { size: 14 }), h('span', null, t)))) : null);
  }

  function paintActiveRow() {
    for (const [sid, tr] of rows) tr.classList.toggle('active', sid === geneSid);
  }

  // ---------------------------------------------------------------- the viewer and its highlighted panel
  const panelEls = () => (viewer ? [...viewer.el.querySelectorAll('.sv-panels > .sv-panel')] : []);
  const panelCount = () => {
    const k = st.viewer && st.viewer.layout;
    return k === 2 || k === 4 ? k : 1;
  };
  const where = () => (panelCount() > 1 ? 'the highlighted panel' : 'the map'); // what the actions apply to

  /** The section of the highlighted panel: the one it shows, or the one a new panel would show (spatial-viewer.js). */
  function panelSample() {
    const ids = readable().map((s) => s.id);
    if (!ids.length) return null;
    const panels = Array.isArray(st.viewer.panels) ? st.viewer.panels : [];
    const p = panels[st.target];
    if (p && ids.includes(p.sample)) return p.sample;
    const shown = panels.slice(0, st.target).map((q) => (q && typeof q === 'object' ? q.sample : null));
    return ids.find((id) => !shown.includes(id)) ?? ids[st.target % ids.length];
  }

  /** Overlays of the highlighted panel (the gene table marks the genes shown as overlays). */
  function targetOverlays() {
    const p = Array.isArray(st.viewer.panels) ? st.viewer.panels[st.target] : null;
    return p && Array.isArray(p.overlays) ? p.overlays : [];
  }

  /** What the gene rows show depends on: the panel count (their titles), the highlighted panel and its overlays. */
  const rowsKey = () => `${panelCount()}\n${st.target}\n${targetOverlays().join('\n')}`;

  /**
   * Everything that follows the highlighted panel: its outline (two or four panels), the section of the gene table and
   * the active row of the section table, and the pressed overlay toggles of the gene rows. Called after every change
   * of the viewer: panels made or removed (the observer), a click, focus or change in it, each viewer.show() from
   * here, and once the viewer is mounted (its panels may come before the observer).
   */
  function markTarget() {
    if (destroyed) return;
    if (st.target >= panelCount()) st.target = 0;
    const list = panelEls();
    list.forEach((p, i) => p.classList.toggle('dx-target', i === st.target && list.length > 1));
    geneJump.dataset.panels = String(panelCount());
    if (titlesFor !== panelCount()) paintSectionTitles();
    syncGeneSection();
    // the overlays of the highlighted panel changed (added here or in the viewer), or another panel is highlighted
    if (rowsKey() !== genesShownFor) paintGeneRows();
  }

  /** The section names say where they show their section: the map, or the highlighted one of several panels. */
  function paintSectionTitles() {
    titlesFor = panelCount();
    for (const [sid, tr] of rows) {
      const b = tr.querySelector('button.dx-section-label');
      if (b) b.title = `Show ${sectionText(sid)} on ${where()}`;
    }
  }

  /** The gene table follows the section of the highlighted panel. */
  function syncGeneSection() {
    const sid = panelSample();
    if (sid !== geneSid) {
      geneSid = sid;
      geneRes = null;
      geneErr = null;
      st.genes.limit = GENE_LIMITS[0];
      if (sid) searchGenes();
    }
    paintActiveRow();
    paintGeneHead();
  }

  // a press, a focus or a change inside a panel makes it the highlighted one; layout buttons, section menus and the
  // gene boxes (Enter, or the mousedown of a suggestion, adds an overlay) are handled by the viewer first (bubbling),
  // so the store is current here
  function onViewerEvent(e) {
    const panel = e.target instanceof Element ? e.target.closest('.sv-panel') : null;
    const i = panel ? panelEls().indexOf(panel) : -1;
    if (i >= 0) st.target = i;
    markTarget();
  }

  /** Bring the highlighted panel into view unless most of its map is visible already. */
  function revealPanel() {
    const p = panelEls()[st.target];
    const stage = p && p.querySelector('.sv-stage');
    if (!stage) return;
    const r = stage.getBoundingClientRect();
    const visible = Math.min(r.bottom, window.innerHeight) - Math.max(r.top, 0);
    if (visible < Math.min(r.height, window.innerHeight) * 0.4) p.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  /** viewer.show() on the highlighted panel, then the table, the outline and the gene rows follow it. */
  function showOnTarget(patch, reveal = true) {
    if (!viewer) return;
    viewer.show({ panel: st.target, ...patch });
    markTarget();
    if (reveal) revealPanel();
  }

  /** The jump above the map: scroll to the gene table and move the focus to it. */
  function jumpToGenes() {
    geneCard.scrollIntoView({ behavior: 'smooth', block: 'start' });
    geneCard.focus({ preventScroll: true });
  }

  const showSection = (sid) => showOnTarget({ sample: sid });

  function showGene(gene) {
    if (geneSid) showOnTarget({ sample: geneSid, genes: [gene] });
  }

  /** Add the "detected" overlay of a gene to the highlighted panel, or remove it when it is there. */
  function toggleGeneOverlay(gene) {
    if (!geneSid) return;
    const key = GENE_OVERLAY + gene;
    const overlays = targetOverlays();
    if (overlays.includes(key)) showOnTarget({ overlays: overlays.filter((k) => k !== key) }, false);
    else showOnTarget({ sample: geneSid, overlays: [...overlays, key] });
  }

  function mountViewer() {
    if (!readable().length) {
      mount(mapHost, card({ body: emptyState({ icon: 'map', title: 'Nothing to map', text: source.emptyText }) }));
      return;
    }
    // labels edited in the Data tab since the sections were listed
    if (Array.isArray(st.viewer.samples)) {
      for (const s of st.viewer.samples) {
        const known = sampleById(s.id);
        if (known) s.label = labelOf(known);
      }
    }
    viewer = spatialViewer({ source, store: st.viewer, binUm: binUmOf(st.binSize || st.config.bin_size) });
    mount(mapHost, viewer.el);
    if (st.target >= panelCount()) st.target = 0;
    if (focus && readable().some((s) => s.id === focus)) viewer.show({ panel: st.target, sample: focus });
    const panelsEl = viewer.el.querySelector('.sv-panels');
    if (panelsEl) {
      observer = new MutationObserver(markTarget);
      observer.observe(panelsEl, { childList: true });
    }
    for (const type of ['pointerdown', 'mousedown', 'keydown', 'focusin', 'click', 'change']) viewer.el.addEventListener(type, onViewerEvent);
    geneCard.classList.remove('hidden');
    geneJump.classList.remove('hidden');
    markTarget(); // panels the viewer made at once (sections known from an earlier visit) came before the observer
  }

  // ---------------------------------------------------------------- gene table
  function paintGeneHead() {
    if (destroyed) return;
    geneTitle.textContent = geneSid ? `Genes of ${sectionText(geneSid)}` : 'Genes';
    geneSub.textContent = panelCount() > 1
      ? 'Ranked by UMIs in the section of the highlighted panel (click a panel to choose it). “Show on map” colors that panel by the gene; “Overlay” marks the bins detecting it (checked while they are marked; click it again to remove the overlay).'
      : 'Ranked by UMIs in the section on the map. “Show on map” colors the map by the gene; “Overlay” marks the bins detecting it (checked while they are marked; click it again to remove the overlay).';
  }

  const searchSoon = debounce(() => searchGenes(), 250);

  async function searchGenes() {
    searchSoon.cancel();
    if (!geneSid || destroyed) return;
    if (geneCtrl) geneCtrl.abort();
    const c = new AbortController();
    geneCtrl = c;
    const my = ++geneSeq;
    const sid = geneSid;
    geneBusy = true;
    geneErr = null;
    paintGenes();
    try {
      const res = await api(url(`/${enc(sid)}/genes`), {
        method: 'POST', body: body({ query: st.genes.query, limit: st.genes.limit }), signal: c.signal,
      });
      if (destroyed || my !== geneSeq) return;
      geneRes = { ...(res && typeof res === 'object' ? res : {}), sid };
    } catch (err) {
      if (destroyed || my !== geneSeq || (err && err.name === 'AbortError')) return;
      geneErr = err;
    } finally {
      if (geneCtrl === c) geneCtrl = null;
    }
    geneBusy = false;
    paintGenes();
  }

  function pctBar(p) {
    const v = Number(p);
    return h('div', { class: 'dx-pct' },
      h('span', { class: 'dx-pct-track', 'aria-hidden': 'true' }, h('span', { class: 'dx-pct-fill', style: { width: `${Number.isFinite(v) ? Math.min(100, Math.max(0, v)) : 0}%` } })),
      h('span', { class: 'dx-pct-text' }, fmt.pct(v, 1)));
  }

  /** Icon and text of a row button; the text gives way to the icon in narrow windows (css/project.css). */
  const buttonFace = (iconName, text) => [icon(iconName, { size: 14 }), h('span', { class: 'dx-btn-text' }, text)];

  /**
   * A gene row: "Show on map" and the "Overlay" toggle, each named by its text and the gene (aria-label; the text
   * comes first, so the name holds what is seen, and the toggle keeps its text and name in both states: its state is
   * aria-pressed and the icon). paintGeneRow() writes what depends on the highlighted panel.
   */
  function geneRow(g) {
    const gene = String(g.gene);
    const show = h('button', {
      class: 'btn small', type: 'button', 'aria-label': `Show on map: ${gene}`, onClick: () => showGene(gene),
    }, buttonFace('map', 'Show on map'));
    const toggle = h('button', {
      class: 'btn small ghost', type: 'button', 'aria-label': `Overlay: ${gene}`, 'aria-pressed': 'false',
      onClick: () => toggleGeneOverlay(gene),
    });
    const tr = h('tr', null,
      h('td', { class: 'mono dx-gene' }, gene),
      h('td', { class: 'num' }, fmt.int(g.total_umi)),
      h('td', { class: 'num dx-col-bins' }, fmt.int(g.n_bins_detected)),
      h('td', { title: `${fmt.int(g.n_bins_detected)} bins detect ${gene}` }, pctBar(g.pct_bins_detected)),
      h('td', { class: 'dx-act' }, h('div', { class: 'dx-gene-actions' }, show, toggle)));
    geneRows.push({ gene, tr, show, toggle });
    return tr;
  }

  /** A gene row as the highlighted panel has it: the toggle pressed, with a check, while the overlay is there. */
  function paintGeneRow({ gene, tr, show, toggle }, overlays, at) {
    const on = overlays.includes(GENE_OVERLAY + gene);
    tr.classList.toggle('dx-on-map', on);
    show.title = `Color ${at} by the UMIs of ${gene}`;
    toggle.title = on ? `The bins of ${at} that detect ${gene} are marked; click to remove the overlay` : `Mark the bins of ${at} that detect ${gene}`;
    if (toggle.firstChild && toggle.getAttribute('aria-pressed') === String(on)) return;
    toggle.setAttribute('aria-pressed', String(on));
    toggle.classList.toggle('active', on);
    mount(toggle, buttonFace(on ? 'check' : 'layers', 'Overlay'));
  }

  /** Bring the rows painted up to date with the highlighted panel, in place (a row button keeps the focus). */
  function paintGeneRows() {
    genesShownFor = rowsKey();
    const overlays = targetOverlays();
    const at = where();
    for (const row of geneRows) paintGeneRow(row, overlays, at);
  }

  function paintGenes() {
    if (destroyed) return;
    geneRows = [];
    genesShownFor = rowsKey();
    const res = geneRes && geneRes.sid === geneSid ? geneRes : null;
    const q = st.genes.query;
    geneBody.classList.toggle('busy', geneBusy && !!res);
    if (geneErr) {
      mount(geneCount);
      mount(geneBody, errorPanel(geneErr, { retry: () => searchGenes() }));
      return;
    }
    if (!res) {
      mount(geneCount);
      mount(geneBody, loading(geneSid ? 'Reading the genes…' : 'Choose a section on the map.'));
      return;
    }
    const genes = Array.isArray(res.genes) ? res.genes.filter((g) => g && g.gene !== undefined) : [];
    const matches = Number.isFinite(res.n_matches) ? res.n_matches : genes.length;
    mount(geneCount, geneBusy ? h('span', { class: 'spinner' }) : null,
      q ? `${fmt.int(matches)} of ${fmt.int(res.n_features)} genes contain “${q}”` : `${fmt.int(res.n_features)} genes`,
      genes.length < matches ? `; the ${fmt.int(genes.length)} with the most UMIs are listed` : '');
    if (!genes.length) {
      mount(geneBody, h('p', { class: 'muted' }, q ? `No gene of this section contains “${q}”.` : 'This section has no genes.'));
      return;
    }
    const next = GENE_LIMITS.find((k) => k > st.genes.limit);
    const more = genes.length < matches && next
      ? button({
        label: `Show up to ${fmt.int(next)}`, icon: 'chevron-down', small: true,
        onClick: () => { st.genes.limit = next; return searchGenes(); },
      })
      : null;
    mount(geneBody,
      h('div', { class: 'table-wrap dx-gene-wrap' }, h('table', { class: 'table compact dx-gene-table' },
        h('thead', null, h('tr', null,
          h('th', null, 'Gene'),
          h('th', { class: 'num', title: 'UMIs of the gene in the whole section' }, 'UMIs'),
          h('th', { class: 'num dx-col-bins', title: 'Bins with at least one UMI of the gene' }, 'Bins detecting it'),
          h('th', { title: 'Share of the bins of the section that detect the gene' }, '% of bins'),
          h('th', null, h('span', { class: 'sr-only' }, 'Actions')))),
        h('tbody', null, genes.map(geneRow)))),
      more ? h('div', { class: 'row' }, more) : null);
    paintGeneRows();
  }

  geneFilter.addEventListener('input', () => {
    st.genes.query = geneFilter.value.trim();
    st.genes.limit = GENE_LIMITS[0];
    searchSoon();
  });
  geneFilter.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') { e.preventDefault(); searchGenes(); }
  });

  // ---------------------------------------------------------------- start, reload, destroy
  async function loadSamples() {
    const g = gen;
    mount(tableHost, loading('Looking for the sections…'));
    try {
      const res = await api(url('/samples'), { method: 'POST', body: body(), signal: ctrl.signal });
      if (g !== gen) return false;
      st.samples = (res && Array.isArray(res.samples) ? res.samples : [])
        .filter((s) => s && s.id !== null && s.id !== undefined)
        .map((s) => ({ ...s, id: String(s.id), ok: s.ok !== false, n_bins: Number.isFinite(s.n_bins) ? s.n_bins : null }));
      st.binSize = res && res.bin_size ? String(res.bin_size) : null;
      return true;
    } catch (err) {
      if (g !== gen || (err && err.name === 'AbortError')) return false;
      mount(tableHost, errorPanel(err, { retry: () => start() }));
      return false;
    }
  }

  async function start() {
    const g = gen;
    if (!st.samples && !(await loadSamples())) return;
    if (g !== gen) return;
    // a visit that ended while a section was being read: read it again
    for (const sec of st.sections.values()) if (sec.status === 'loading') sec.status = 'idle';
    paintTable();
    if (!st.samples.length) {
      mount(qcLine);
      mount(notesHost);
      return;
    }
    paintFoot();
    mountViewer();
    const shown = Array.isArray(st.viewer.panels) ? st.viewer.panels.slice(0, panelCount()).map((p) => p && p.sample) : [];
    readAll([focus, panelSample(), ...shown].filter(Boolean));
  }

  function stop() {
    gen += 1;
    ctrl.abort();
    ctrl = new AbortController();
    if (geneCtrl) geneCtrl.abort();
    geneCtrl = null;
    searchSoon.cancel();
    if (observer) observer.disconnect();
    observer = null;
    if (viewer) viewer.destroy();
    viewer = null;
    payloads = new Map();
    ahead = [];
  }

  function reload() {
    stop();
    st.config = clone(getConfig());
    st.samples = null;
    st.binSize = null;
    st.sections = new Map();
    st.viewer = carryViewer(st.viewer);
    focus = null;
    geneSid = null;
    geneRes = null;
    geneErr = null;
    geneCard.classList.add('hidden');
    geneJump.classList.add('hidden');
    mount(mapHost);
    start();
  }

  start();

  return {
    el,
    destroy() {
      destroyed = true;
      stop();
    },
  };
}
