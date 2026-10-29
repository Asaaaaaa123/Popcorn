// run.js — one run: header (status, times, actions) and the tabs progress | summary | figures | tables | viewer |
// config. The run record is polled every 1.5 s while the run is queued or running. Results are fetched when their tab
// is first opened and kept in memory while the view lives; they are fetched again once after the run finishes.
import {
  h, api, mount, button, icon, card, tabs, statusBadge, loading, emptyState, errorPanel, banner, toast,
  confirmDialog, promptDialog, fmt, copyText, navigate, setPathSilently, clone, refreshHealth, state,
} from '../core.js';
import { logView, downloadText } from '../components/log-view.js';
import { tableBrowser, fileUrl, formatNumber } from '../components/table-viewer.js';
import { figureGallery } from '../components/figure-gallery.js';
import { spatialViewer } from '../components/spatial-viewer.js';

const ACTIVE = new Set(['queued', 'running']);
const FAILED = new Set(['failed', 'interrupted']);
const POLL_MS = 1500;
const TABS = [
  { id: 'progress', label: 'Progress', icon: 'activity' },
  { id: 'summary', label: 'Summary', icon: 'chart' },
  { id: 'figures', label: 'Figures', icon: 'image' },
  { id: 'tables', label: 'Tables', icon: 'table' },
  { id: 'viewer', label: 'Spatial viewer', icon: 'map' },
  { id: 'config', label: 'Configuration', icon: 'sliders' },
];
const RESULT_TABS = new Set(['summary', 'figures', 'tables', 'viewer']);
const STEP_TEXT = { pending: 'pending', running: 'running', done: 'done', skipped: 'skipped', failed: 'failed' };
const WARNINGS_SHOWN = 8;

const isColor = (c) => typeof c === 'string' && /^#(?:[0-9a-f]{3}|[0-9a-f]{6})$/i.test(c.trim());
const numOrNull = (v) => (v === null || v === undefined || v === '' || !Number.isFinite(Number(v)) ? null : Number(v));

function sinceSeconds(iso) {
  const t = iso ? new Date(iso).getTime() : NaN;
  return Number.isFinite(t) ? Math.max(0, (Date.now() - t) / 1000) : null;
}

function secondsText(s) {
  const v = numOrNull(s);
  if (v === null) return '';
  return v < 10 ? `${v.toFixed(1)} s` : fmt.duration(v);
}

function binUmOf(cfg) {
  const m = cfg && /^(\d{3})um$/.exec(String(cfg.bin_size || ''));
  return m ? Number(m[1]) : NaN;
}

/** Element that shows the time since `iso` and is kept current by the view's one-second ticker. */
function liveDuration(iso) {
  return h('span', { dataset: { since: iso } }, fmt.duration(sinceSeconds(iso)));
}

/** Pretty-printed JSON with keys, strings, numbers and literals marked for coloring (text nodes only). */
function jsonView(value) {
  const text = JSON.stringify(value, null, 2);
  const re = /("(?:\\.|[^"\\])*")(\s*:)?|\b(true|false|null)\b|-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?/g;
  const pre = h('pre', { class: 'run-json', tabIndex: 0, 'aria-label': 'Configuration (JSON)' });
  let last = 0;
  let m;
  while ((m = re.exec(text))) {
    if (m.index > last) pre.append(text.slice(last, m.index));
    if (m[1]) {
      pre.append(h('span', { class: m[2] ? 'j-key' : 'j-str' }, m[1]));
      if (m[2]) pre.append(m[2]);
    } else {
      pre.append(h('span', { class: m[3] ? 'j-lit' : 'j-num' }, m[0]));
    }
    last = re.lastIndex;
  }
  pre.append(text.slice(last));
  return pre;
}

// ------------------------------------------------------------------------------------------------ progress tab
function stepIcon(status) {
  if (status === 'done') return icon('check', { size: 15 });
  if (status === 'failed') return icon('x', { size: 15 });
  if (status === 'running') return h('span', { class: 'spinner' });
  if (status === 'skipped') return h('span', { class: 'run-dash' });
  return h('span', { class: 'run-ring' });
}

function progressHeadline(r) {
  const p = r.progress || {};
  const steps = Array.isArray(r.steps) ? r.steps : [];
  const total = Number(p.total) || steps.length;
  const done = steps.filter((s) => s.status === 'done').length;
  const skipped = steps.filter((s) => s.status === 'skipped').length;
  const failedAt = steps.findIndex((s) => s.status === 'failed');
  switch (r.status) {
    case 'queued': return 'Waiting in the queue';
    case 'running':
      return total && Number(p.index) ? `Step ${Math.min(Number(p.index), total)} of ${total} — ${p.label || 'working'}` : 'Starting the analysis…';
    case 'succeeded':
      return `Finished: ${done} step${done === 1 ? '' : 's'} done${skipped ? `, ${skipped} skipped` : ''}`;
    case 'cancelled': return 'Cancelled';
    case 'interrupted': return failedAt >= 0 ? `Interrupted at step ${failedAt + 1} of ${total} — ${steps[failedAt].label}` : 'Interrupted';
    case 'failed': return failedAt >= 0 ? `Failed at step ${failedAt + 1} of ${total} — ${steps[failedAt].label}` : 'Failed';
    default: return String(r.status || 'Unknown status');
  }
}

// ------------------------------------------------------------------------------------------------ summary tab
function metricTile(m) {
  const v = m.value;
  let text;
  if (v === null || v === undefined || v === '') text = '—';
  else if (typeof v === 'number') text = formatNumber(v);
  else if (typeof v === 'boolean') text = v ? 'yes' : 'no';
  else text = String(v);
  return h('div', { class: 'metric' },
    h('div', { class: 'metric-label' }, m.label || ''),
    h('div', { class: 'metric-value' }, text, m.unit && text !== '—' ? h('span', { class: 'metric-unit' }, m.unit) : null),
    m.detail ? h('div', { class: 'metric-detail' }, String(m.detail)) : null);
}

function sectionTile(sm, isReference) {
  const color = isColor(sm.color) ? sm.color.trim() : 'var(--accent)';
  const loaded = numOrNull(sm.n_bins_loaded);
  const qc = numOrNull(sm.n_bins_qc);
  const roi = numOrNull(sm.n_bins_roi);
  const pctQc = numOrNull(sm.pct_retained) ?? (loaded && qc !== null ? (100 * qc) / loaded : null);
  const pctRoi = qc && roi !== null ? (100 * roi) / qc : null;
  const width = (n) => (loaded && n !== null ? Math.max(0, Math.min(100, (100 * n) / loaded)) : 0);
  const stage = (label, n, note) => h('li', null,
    h('div', { class: 'sum-funnel-row' }, h('span', { class: 'muted' }, label), h('strong', null, n === null ? '—' : fmt.int(n))),
    h('div', { class: 'sum-funnel-bar', 'aria-hidden': 'true' }, h('span', { style: { width: `${width(n)}%`, background: color } })),
    note ? h('div', { class: 'faint small' }, note) : null);
  const mpp = numOrNull(sm.microns_per_pixel);
  return h('section', { class: 'sum-section' },
    h('div', { class: 'sum-section-head' },
      h('span', { class: 'swatch', style: { background: color } }),
      h('h3', null, sm.label || sm.id || 'section'),
      isReference ? h('span', { class: 'badge', title: 'Fold changes and ratios use this section as the denominator' }, 'reference') : null),
    h('div', { class: 'faint small' }, [sm.id, sm.condition].filter(Boolean).join(' · ')),
    h('ol', { class: 'sum-funnel' },
      stage('Bins loaded', loaded),
      stage('After quality control', qc, pctQc !== null ? `${fmt.num(pctQc, 1)}% retained` : null),
      stage('In region of interest', roi, pctRoi !== null ? `${fmt.num(pctRoi, 1)}% of the bins after QC` : null)),
    h('dl', { class: 'kv sum-section-kv' },
      h('dt', null, 'Median UMIs per bin'), h('dd', null, formatNumber(numOrNull(sm.median_umi))),
      h('dt', null, 'Median genes per bin'), h('dd', null, formatNumber(numOrNull(sm.median_genes))),
      h('dt', null, 'Microns per pixel'), h('dd', null, mpp === null ? '—' : fmt.num(mpp, 4))));
}

function summaryView(s) {
  const samples = (Array.isArray(s.samples) ? s.samples : []).filter((x) => x && typeof x === 'object');
  const metrics = (Array.isArray(s.metrics) ? s.metrics : []).filter((m) => m && typeof m === 'object');
  const warnings = (Array.isArray(s.warnings) ? s.warnings : []).filter(Boolean).map(String);
  const skipped = (Array.isArray(s.skipped) ? s.skipped : []).filter((x) => x && typeof x === 'object');
  const timings = (Array.isArray(s.timings) ? s.timings : []).filter((x) => x && typeof x === 'object');
  const env = s.environment && typeof s.environment === 'object' ? Object.entries(s.environment) : [];
  const ref = samples.find((x) => x.id === s.reference_sample);
  const bin = /^(\d{3})um$/.exec(String(s.bin_size || ''));
  const target = numOrNull(s.depth_match_target_umi);
  const parts = [];

  let warningsCard = null;
  if (warnings.length) {
    warningsCard = card({
      cls: 'sum-warnings', title: `Warnings (${warnings.length})`, subtitle: 'Reported by the pipeline; the numbers above are complete unless a warning says otherwise',
      body: h('ul', { class: 'run-warning-list' }, warnings.map((w) => h('li', null, w))),
    });
  }
  parts.push(card({
    title: 'Overview',
    actions: warnings.length ? h('button', {
      class: 'btn small', type: 'button', onClick: () => warningsCard.scrollIntoView({ behavior: 'smooth', block: 'start' }),
    }, icon('alert', { size: 14 }), `${warnings.length} warning${warnings.length === 1 ? '' : 's'}`) : null,
    body: h('dl', { class: 'kv sum-overview' },
      h('dt', null, 'Finished'), h('dd', null, fmt.date(s.finished)),
      h('dt', null, 'Bin size'), h('dd', null, bin ? `${Number(bin[1])} µm (square_${s.bin_size})` : (s.bin_size || '—')),
      h('dt', null, 'Reference section'), h('dd', null, ref ? `${ref.label || ref.id} (${ref.id})` : (s.reference_sample || '—')),
      h('dt', null, 'Depth-matched to'), h('dd', null, target === null ? '—' : [
        `${fmt.int(target)} UMIs per bin`,
        h('span', { class: 'faint' }, ' — the median of the shallowest section; depth matching downsamples the bins above it (see scoring.depth_matching)'),
      ])),
  }));

  const groups = new Map();
  for (const m of metrics) {
    const g = m.group || 'Other';
    if (!groups.has(g)) groups.set(g, []);
    groups.get(g).push(m);
  }
  if (samples.length) {
    // the pipeline's "Sections" metrics belong with the section tiles
    const sectionMetrics = groups.get('Sections');
    groups.delete('Sections');
    parts.push(card({
      title: 'Sections', subtitle: 'Bins loaded, kept by quality control and inside the region of interest',
      body: [
        h('div', { class: 'sum-sections' }, samples.map((sm) => sectionTile(sm, sm.id === s.reference_sample))),
        sectionMetrics ? h('div', { class: 'sum-metrics sum-section-metrics' }, sectionMetrics.map(metricTile)) : null,
      ],
    }));
  }
  if (groups.size) {
    // small groups share a row; groups with many numbers take the full width
    parts.push(h('div', { class: 'sum-groups' }, [...groups].map(([name, list]) => card({
      cls: list.length > 3 ? 'sum-group wide' : 'sum-group', title: name,
      body: h('div', { class: 'sum-metrics' }, list.map(metricTile)),
    }))));
  } else if (!metrics.length) {
    parts.push(card({ body: emptyState({ icon: 'chart', title: 'No metrics', text: 'The summary lists no key numbers; the Tables tab has every output.' }) }));
  }
  if (warningsCard) parts.push(warningsCard);

  const extra = [];
  if (skipped.length) {
    extra.push(card({
      title: `Skipped steps (${skipped.length})`, flush: true,
      body: h('div', { class: 'table-wrap borderless' }, h('table', { class: 'table compact' },
        h('thead', null, h('tr', null, h('th', { scope: 'col' }, 'Step'), h('th', { scope: 'col' }, 'Reason'))),
        h('tbody', null, skipped.map((k) => h('tr', null, h('td', null, k.step || '—'), h('td', { class: 'muted' }, k.reason || '—')))))),
    }));
  }
  if (timings.length) {
    const secs = timings.map((t) => Math.max(0, numOrNull(t.seconds) || 0));
    const max = Math.max(0, ...secs);
    const total = secs.reduce((a, b) => a + b, 0);
    extra.push(card({
      title: 'Timings', subtitle: `${fmt.duration(total)} in total`, flush: true,
      body: h('div', { class: 'table-wrap borderless' }, h('table', { class: 'table compact sum-timings' },
        h('thead', null, h('tr', null, h('th', { scope: 'col' }, 'Step'), h('th', { scope: 'col' }, h('span', { class: 'sr-only' }, 'Share')), h('th', { scope: 'col', class: 'num' }, 'Time'))),
        h('tbody', null, timings.map((t, i) => h('tr', null,
          h('td', null, t.step || '—'),
          h('td', { class: 'sum-bar-cell' }, h('div', { class: 'sum-bar', 'aria-hidden': 'true' },
            h('span', { style: { width: `${max ? (100 * secs[i]) / max : 0}%` } }))),
          h('td', { class: 'num nowrap' }, secondsText(secs[i]))))))),
    }));
  }
  if (extra.length) parts.push(h('div', { class: extra.length > 1 ? 'grid-2 sum-extra' : 'stack' }, extra));
  if (env.length) {
    parts.push(card({
      title: 'Environment', subtitle: 'Package versions this run used',
      body: h('dl', { class: 'kv sum-env' }, env.map(([k, v]) => [
        h('dt', null, k), h('dd', { class: v ? 'mono' : 'faint' }, v ? String(v) : 'not installed')])),
    }));
  }
  return h('div', { class: 'stack-lg run-summary' }, parts);
}

// ------------------------------------------------------------------------------------------------ view
export async function renderRun(container, { id, tab }) {
  const enc = encodeURIComponent(id);
  let first;
  try {
    first = await api(`/api/runs/${enc}`);
  } catch (err) {
    if (err.status === 404) {
      mount(container, h('div', { class: 'page' }, card({
        body: emptyState({
          icon: 'alert', title: 'Run not found', text: `There is no run ${id}; it may have been deleted.`,
          action: h('a', { class: 'btn primary', href: '#/runs' }, 'All runs'),
        }),
      })));
      return undefined;
    }
    throw err;
  }

  const st = {
    run: first, destroyed: false, gone: false, cache: {},
    logStore: {}, tableStore: {}, viewerStore: {}, tracebackOpen: false, warningsOpen: false,
  };
  const known = TABS.some((t) => t.id === tab);
  st.tab = known ? tab : (ACTIVE.has(first.status) ? 'progress' : 'summary');
  if (tab && !known) setPathSilently(`/runs/${enc}/${st.tab}`);

  const headerEl = h('div');
  const noticeEl = h('div');
  const tabsEl = h('div');
  const bodyEl = h('div', { class: 'run-body' });
  mount(container, h('div', { class: 'page run-page' }, headerEl, noticeEl, tabsEl, bodyEl));

  let pane = null;
  let paneToken = 0;
  let headerSig = null;
  let tabsSig = null;
  let pollTimer = null;
  let pollCtrl = null;
  let pollFailures = 0;
  let ticker = null;

  const isActive = () => ACTIVE.has(st.run.status);

  // ---- results, fetched lazily; partial results of an active run are never kept
  function fetchResult(key, path) {
    const slot = st.cache[key];
    if (slot && 'value' in slot) return Promise.resolve(slot.value);
    if (slot && slot.promise) return slot.promise;
    const entry = {};
    const keep = !isActive();
    st.cache[key] = entry;
    entry.promise = api(path)
      .catch((err) => { if (err.status === 404) return null; throw err; })
      .then((v) => { if (keep && st.cache[key] === entry) entry.value = v; return v; })
      .finally(() => {
        delete entry.promise;
        if (!('value' in entry) && st.cache[key] === entry) delete st.cache[key];
      });
    return entry.promise;
  }
  const getSummary = () => fetchResult('summary', `/api/runs/${enc}/summary`);
  const getManifest = () => fetchResult('manifest', `/api/runs/${enc}/manifest`);
  const cachedManifest = () => (st.cache.manifest && st.cache.manifest.value) || null;

  // ---- header
  function paintHeader() {
    const r = st.run;
    const active = isActive();
    const title = r.label || r.project_name || r.id;
    const sig = JSON.stringify([r.status, title, r.project_id, r.project_name, fmt.ago(r.created), r.started, r.finished, r.duration_s]);
    if (sig === headerSig) return;
    headerSig = sig;
    const meta = [];
    if (r.project_id) {
      meta.push(h('span', null, 'Project ', h('a', { href: `#/projects/${encodeURIComponent(r.project_id)}` }, r.project_name || r.project_id)));
    }
    meta.push(h('span', { title: fmt.date(r.created) }, `Created ${fmt.ago(r.created)}`));
    if (r.started) meta.push(h('span', null, `Started ${fmt.date(r.started)}`));
    if (r.finished) meta.push(h('span', null, `Finished ${fmt.date(r.finished)}`));
    if (r.status === 'running' && r.started) meta.push(h('span', null, 'Running for ', liveDuration(r.started)));
    else if (numOrNull(r.duration_s) !== null) meta.push(h('span', null, `Duration ${fmt.duration(r.duration_s)}`));
    const actions = [];
    if (active) actions.push(button({ label: 'Cancel', icon: 'stop', kind: 'danger', onClick: cancelRun }));
    actions.push(button({ label: 'Re-run', icon: 'refresh', title: 'Start a new run with the same configuration snapshot and regions', onClick: rerun }));
    if (r.project_id) {
      actions.push(h('a', { class: 'btn', href: `#/projects/${encodeURIComponent(r.project_id)}` }, icon('folder', { size: 16 }), 'Open project'));
    }
    if (!active) actions.push(button({ label: 'Delete', icon: 'trash', kind: 'danger', onClick: deleteRun }));
    mount(headerEl, h('div', { class: 'page-header run-header' },
      h('div', { class: 'titles' },
        h('nav', { class: 'breadcrumb', 'aria-label': 'Breadcrumb' },
          h('a', { href: '#/runs' }, 'Runs'), h('span', { 'aria-hidden': 'true' }, '/'), h('span', { class: 'mono' }, r.id)),
        h('div', { class: 'run-title' }, h('h1', null, title), statusBadge(r.status)),
        h('div', { class: 'subtitle run-meta' }, meta)),
      h('div', { class: 'actions' }, actions)));
  }

  function setNotice(node) {
    mount(noticeEl, node);
  }

  // ---- tabs
  function badgeFor(tid) {
    const r = st.run;
    if (tid === 'progress') {
      if (FAILED.has(r.status)) return { badge: '!', badgeKind: 'bad' };
      const n = Array.isArray(r.warnings) ? r.warnings.length : 0;
      return n ? { badge: n, badgeKind: 'warn' } : {};
    }
    const m = cachedManifest();
    if (!m) return {};
    if (tid === 'figures' && Array.isArray(m.figures) && m.figures.length) {
      return { badge: m.figures.length, badgeKind: m.figures.some((f) => f && f.error) ? 'warn' : '' };
    }
    if (tid === 'tables' && Array.isArray(m.tables) && m.tables.length) return { badge: m.tables.length };
    return {};
  }

  function paintTabs(force = false) {
    const items = TABS.map((t) => ({ ...t, ...badgeFor(t.id) }));
    const sig = JSON.stringify([st.tab, items.map((i) => [i.badge, i.badgeKind])]);
    if (!force && sig === tabsSig) return;
    tabsSig = sig;
    const hadFocus = tabsEl.contains(document.activeElement);
    const strip = tabs({ items, active: st.tab, onSelect: selectTab });
    strip.setAttribute('aria-label', 'Run sections');
    strip.addEventListener('keydown', (e) => {
      const i = TABS.findIndex((t) => t.id === st.tab);
      let j;
      if (e.key === 'ArrowRight') j = (i + 1) % TABS.length;
      else if (e.key === 'ArrowLeft') j = (i - 1 + TABS.length) % TABS.length;
      else if (e.key === 'Home') j = 0;
      else if (e.key === 'End') j = TABS.length - 1;
      else return;
      e.preventDefault();
      selectTab(TABS[j].id);
      const activeBtn = tabsEl.querySelector('.tab.active');
      if (activeBtn) activeBtn.focus();
    });
    mount(tabsEl, strip);
    if (hadFocus) {
      const activeBtn = tabsEl.querySelector('.tab.active');
      if (activeBtn) activeBtn.focus();
    }
  }

  function selectTab(t) {
    if (st.destroyed || !TABS.some((x) => x.id === t)) return;
    if (t === st.tab && pane) return;
    st.tab = t;
    setPathSilently(`/runs/${enc}/${t}`);
    paintTabs(true);
    renderPane();
  }

  // ---- panes: each returns {el, update?, mounted?, destroy?} (or a promise of one)
  function renderPane() {
    if (pane && pane.destroy) pane.destroy();
    pane = null;
    const token = ++paneToken;
    const builders = { progress: progressPane, summary: summaryPane, figures: figuresPane, tables: tablesPane, viewer: viewerPane, config: configPane };
    const show = (p) => {
      pane = p;
      mount(bodyEl, p.el);
      if (p.mounted) p.mounted();
    };
    let result;
    try {
      result = builders[st.tab]();
    } catch (err) {
      mount(bodyEl, errorPanel(err, { retry: renderPane }));
      return;
    }
    if (!result || typeof result.then !== 'function') { show(result); return; }
    mount(bodyEl, loading());
    result.then((p) => {
      if (token !== paneToken || st.destroyed) { if (p && p.destroy) p.destroy(); return; }
      show(p);
    }, (err) => {
      if (token !== paneToken || st.destroyed) return;
      mount(bodyEl, errorPanel(err, { retry: renderPane }));
    });
  }

  function progressPane() {
    const top = h('div', { class: 'stack-lg' });
    const lv = logView({ runId: id, store: st.logStore, isFinished: () => !isActive() });
    let sig = null;
    const paint = () => {
      const r = st.run;
      const next = JSON.stringify([r.status, r.progress, r.steps, r.warnings, r.error, r.traceback, r.started, r.finished, r.duration_s]);
      if (next === sig) return;
      sig = next;
      mount(top, progressContent());
    };
    paint();
    return {
      el: h('div', { class: 'stack-lg' }, top, lv.el),
      update: paint,
      mounted: () => lv.mounted(),
      destroy: () => lv.destroy(),
    };
  }

  function progressContent() {
    const r = st.run;
    const p = r.progress || {};
    const steps = (Array.isArray(r.steps) ? r.steps : []).filter((s) => s && typeof s === 'object');
    const warnings = (Array.isArray(r.warnings) ? r.warnings : []).filter(Boolean).map(String);
    const total = Number(p.total) || steps.length;
    const frac = r.status === 'succeeded' ? 1 : Math.min(1, Math.max(0, Number(p.fraction) || 0));
    const indeterminate = r.status === 'queued' || (r.status === 'running' && !total);
    const times = [];
    if (r.status === 'queued') {
      times.push(h('span', null, 'Waiting for ', liveDuration(r.created)));
    } else if (r.status === 'running' && r.started) {
      times.push(h('span', null, 'Elapsed ', liveDuration(r.started)));
    } else if (numOrNull(r.duration_s) !== null) {
      times.push(h('span', null, `Duration ${fmt.duration(r.duration_s)}`));
    }
    if (r.started) times.push(h('span', null, `Started ${fmt.date(r.started)}`));
    if (r.finished) times.push(h('span', null, `Finished ${fmt.date(r.finished)}`));

    const out = [];
    if (r.status === 'queued') {
      const activeRun = state.health && state.health.active_run;
      out.push(banner('info', ['One analysis runs at a time; this run starts when the ones ahead of it have finished.',
        activeRun && activeRun !== r.id ? [' ', h('a', { href: `#/runs/${encodeURIComponent(activeRun)}` }, 'Open the running analysis')] : null]));
    } else if (FAILED.has(r.status)) {
      const tb = r.traceback ? h('details', {
        class: 'run-traceback', open: st.tracebackOpen,
        onToggle: (e) => { st.tracebackOpen = e.currentTarget.open; },
      }, h('summary', null, 'Traceback'),
      h('div', { class: 'run-traceback-body' },
        button({ label: 'Copy', icon: 'copy', small: true, onClick: () => copyText(r.traceback) }),
        h('pre', { class: 'run-pre' }, r.traceback))) : null;
      out.push(h('div', { class: 'banner danger run-failure', role: 'alert' }, icon('alert', { size: 18 }),
        h('div', { class: 'banner-body' },
          h('strong', null, r.status === 'interrupted' ? 'The run was interrupted.' : 'The analysis failed.'),
          r.error ? h('pre', { class: 'run-pre run-error' }, String(r.error)) : null,
          tb,
          h('div', { class: 'faint small' }, 'Fix the cause (often a path, a missing file or a parameter), then use Re-run, or change the project and start a new run.'))));
    } else if (r.status === 'cancelled') {
      out.push(banner('warning', 'This run was cancelled before it finished; outputs written until then are kept. Re-run starts it again from the beginning.'));
    }

    out.push(card({
      cls: 'run-progress-card',
      body: h('div', { class: 'stack' },
        h('div', { class: 'run-headline' }, h('strong', null, progressHeadline(r)),
          h('span', { class: 'muted' }, total && !indeterminate ? `${Math.round(frac * 100)}%` : '')),
        h('div', {
          class: ['progress', 'run-progress', indeterminate ? 'indeterminate' : '', `is-${r.status}`], role: 'progressbar',
          'aria-label': 'Run progress', 'aria-valuemin': '0', 'aria-valuemax': '100',
          'aria-valuenow': indeterminate ? null : String(Math.round(frac * 100)),
        }, h('div', { class: 'bar', style: { width: `${(frac * 100).toFixed(1)}%` } })),
        times.length ? h('div', { class: 'run-times muted small' }, times) : null),
    }));

    if (steps.length) {
      out.push(card({
        title: 'Steps', subtitle: `${total} planned`, flush: true,
        body: h('ol', { class: 'run-steps' }, steps.map((s) => {
          const status = STEP_TEXT[s.status] ? s.status : 'pending';
          return h('li', { class: `run-step is-${status}` },
            h('span', { class: 'run-step-icon', 'aria-hidden': 'true' }, stepIcon(status)),
            h('div', { class: 'run-step-main' },
              h('div', { class: 'run-step-label' }, s.label || '—', h('span', { class: 'sr-only' }, ` (${STEP_TEXT[status]})`)),
              s.reason ? h('div', { class: 'run-step-reason' }, s.reason) : null),
            h('span', { class: 'run-step-time' }, status === 'running' ? 'running…' : secondsText(s.seconds)));
        })),
      }));
    } else if (r.status === 'running') {
      out.push(h('p', { class: 'muted' }, 'The worker is starting; the list of steps appears in a moment.'));
    }

    if (warnings.length) {
      const list = (items) => h('ul', { class: 'run-warning-list' }, items.map((w) => h('li', null, w)));
      out.push(card({
        cls: 'run-warnings', title: `Warnings (${warnings.length})`,
        body: [list(warnings.slice(0, WARNINGS_SHOWN)),
          warnings.length > WARNINGS_SHOWN ? h('details', {
            open: st.warningsOpen, onToggle: (e) => { st.warningsOpen = e.currentTarget.open; },
          }, h('summary', null, `Show ${warnings.length - WARNINGS_SHOWN} more`), list(warnings.slice(WARNINGS_SHOWN))) : null],
      }));
    }
    return out;
  }

  async function summaryPane() {
    const s = await getSummary();
    if (s && typeof s === 'object') return { el: summaryView(s) };
    const r = st.run;
    let text;
    let iconName = 'alert';
    if (isActive()) {
      iconName = 'clock';
      text = 'Results appear when the run finishes. The Progress tab shows how far it has got.';
    } else if (r.status === 'succeeded') {
      text = 'The run finished but wrote no summary (summary.json). The Tables tab lists the outputs it wrote.';
    } else {
      text = 'This run stopped before writing its summary. The Progress tab shows what happened; the Tables tab lists any outputs written before it stopped.';
    }
    return {
      el: card({
        body: emptyState({
          icon: iconName, title: isActive() ? 'No results yet' : 'No summary', text,
          action: button({ label: 'Show progress', icon: 'activity', small: true, onClick: () => selectTab('progress') }),
        }),
      }),
    };
  }

  async function figuresPane() {
    const m = await getManifest();
    paintTabs();
    const figures = m && Array.isArray(m.figures) ? m.figures : [];
    const cfg = st.run.config || {};
    if (!figures.length) {
      let text;
      if (cfg.steps && cfg.steps.figures === false) {
        text = 'Figures were not drawn because the “Draw figures” step was turned off for this run. Turn it on and run again, or draw them from the command line with bone_figures.py.';
      } else if (isActive()) {
        text = 'Figures are drawn at the end of the run, from its tables.';
      } else if (st.run.status === 'succeeded') {
        text = 'This run has no figures.';
      } else {
        text = 'The run stopped before any figure was drawn.';
      }
      return { el: card({ body: emptyState({ icon: 'image', title: 'No figures', text }) }) };
    }
    const g = figureGallery({ figures });
    return {
      el: h('div', { class: 'stack' },
        isActive() ? banner('info', 'The run is still in progress; more figures may appear when it finishes.') : null, g.el),
      destroy: g.destroy,
    };
  }

  async function tablesPane() {
    const m = await getManifest();
    paintTabs();
    const has = m && ['tables', 'objects', 'other'].some((k) => Array.isArray(m[k]) && m[k].length);
    if (!has) {
      return {
        el: card({
          body: emptyState({
            icon: 'table', title: 'No tables',
            text: isActive() ? 'Tables appear as the analysis writes them.' : 'This run wrote no result files.',
          }),
        }),
      };
    }
    const tb = tableBrowser({ runId: id, manifest: m, store: st.tableStore });
    return {
      el: h('div', { class: 'stack' },
        isActive() ? banner('info', 'The run is still in progress; this list shows the files written so far.') : null, tb.el),
      destroy: tb.destroy,
    };
  }

  function viewerPane() {
    let emptyText;
    if (isActive()) emptyText = 'The per-bin map is available when the run finishes (it reads the per-bin signature calls).';
    else if (st.run.status !== 'succeeded') emptyText = 'The run stopped before writing the per-bin table (per_bin_signature_calls), so there is nothing to map.';
    const v = spatialViewer({ runId: id, store: st.viewerStore, binUm: binUmOf(st.run.config), emptyText });
    return { el: v.el, destroy: v.destroy };
  }

  function configPane() {
    const cfg = st.run.config;
    if (!cfg || typeof cfg !== 'object') {
      return { el: card({ body: emptyState({ icon: 'sliders', title: 'No configuration', text: 'The configuration snapshot of this run could not be read.' }) }) };
    }
    const text = JSON.stringify(cfg, null, 2);
    const fileName = `${id}-config.json`;
    // results/config.json exists once the analysis has started writing; before that the snapshot is downloaded as shown
    const download = st.run.status === 'succeeded'
      ? h('a', { class: 'btn small', href: fileUrl(id, 'config.json', true), download: fileName }, icon('download', { size: 14 }), 'Download')
      : button({ label: 'Download', icon: 'download', small: true, onClick: () => downloadText(fileName, text, 'application/json') });
    const samples = (Array.isArray(cfg.samples) ? cfg.samples : []).filter((x) => x && typeof x === 'object');
    const stepsOff = Object.entries(cfg.steps && typeof cfg.steps === 'object' ? cfg.steps : {}).filter(([, on]) => on === false).map(([k]) => k);
    const masks = cfg.roi && cfg.roi.masks && typeof cfg.roi.masks === 'object' ? Object.keys(cfg.roi.masks) : [];
    const overview = card({
      title: 'At a glance',
      body: h('div', { class: 'stack' },
        h('dl', { class: 'kv' },
          h('dt', null, 'Data folder'), h('dd', { class: 'mono' }, String(cfg.data_root ?? '—')),
          h('dt', null, 'Bin size'), h('dd', null, `${cfg.bin_size || '—'} (QC figure: ${cfg.qc_bin_size || '—'})`),
          h('dt', null, 'Reference section'), h('dd', null, String(cfg.reference_sample ?? '—')),
          h('dt', null, 'Hand-drawn region'), h('dd', null, cfg.roi && cfg.roi.use_drawn
            ? (masks.length ? `used for ${masks.join(', ')}` : 'on, but no region was drawn') : 'not used'),
          h('dt', null, 'Steps turned off'), h('dd', null, stepsOff.length ? stepsOff.join(', ') : 'none'),
          h('dt', null, 'Deterministic'), h('dd', null, cfg.deterministic === false ? 'no' : 'yes')),
        samples.length ? h('div', { class: 'table-wrap' }, h('table', { class: 'table compact' },
          h('thead', null, h('tr', null, ['Section', 'Label', 'Condition', 'Path'].map((c) => h('th', { scope: 'col' }, c)))),
          h('tbody', null, samples.map((sm) => h('tr', null,
            h('td', null, h('span', { class: 'swatch', style: { background: isColor(sm.color) ? sm.color : 'transparent' } }), ' ', h('span', { class: 'mono' }, sm.id || '—')),
            h('td', null, sm.label || '—'), h('td', null, sm.condition || '—'),
            h('td', { class: 'mono small' }, sm.path || '—')))))) : null),
    });
    return {
      el: h('div', { class: 'stack-lg' }, overview, card({
        title: 'config.json', subtitle: 'The configuration snapshot this run used, defaults included',
        actions: [
          button({ label: 'Copy', icon: 'copy', small: true, onClick: () => copyText(text) }),
          download,
          button({ label: 'Create project from this configuration', icon: 'plus', small: true, kind: 'primary', onClick: createProject }),
        ],
        flush: true, body: jsonView(cfg),
      })),
    };
  }

  // ---- actions
  async function cancelRun() {
    const running = st.run.status === 'running';
    const ok = await confirmDialog({
      title: 'Cancel this run?',
      message: running
        ? 'The analysis process is stopped. Outputs written so far stay in the run folder, but the run is incomplete.'
        : 'The run is taken out of the queue and will not start.',
      confirmLabel: 'Cancel run', danger: true,
    });
    if (!ok || st.destroyed) return;
    const rec = await api(`/api/runs/${enc}/cancel`, { method: 'POST' });
    refreshHealth();
    if (!st.destroyed && rec && rec.id === id) applyRun(rec);
  }

  async function rerun() {
    const rec = await api(`/api/runs/${enc}/rerun`, { method: 'POST' });
    refreshHealth();
    if (!rec || !rec.id) throw new Error('The server did not return the new run.');
    toast(`Run ${rec.id} queued`, 'success');
    navigate(`/runs/${encodeURIComponent(rec.id)}`);
  }

  async function deleteRun() {
    const ok = await confirmDialog({
      title: 'Delete this run?',
      message: `Run ${id} and all of its results, figures and log are removed from the workspace. This cannot be undone.`,
      confirmLabel: 'Delete run', danger: true,
    });
    if (!ok) return;
    await api(`/api/runs/${enc}`, { method: 'DELETE' });
    toast('Run deleted', 'success');
    navigate('/runs');
  }

  async function createProject() {
    const cfg = clone(st.run.config);
    const hadMasks = !!(cfg.roi && cfg.roi.masks && Object.keys(cfg.roi.masks).length);
    // the masks point at this run's copies of the drawn regions; a new project draws its own
    if (cfg.roi && typeof cfg.roi === 'object') cfg.roi.masks = {};
    const name = await promptDialog({
      title: 'Create a project from this configuration',
      label: 'Name of the new project (samples, parameters and gene signatures are copied)',
      value: `${st.run.project_name || 'Project'} (from run ${id})`,
      confirmLabel: 'Create project',
    });
    if (!name) return;
    const project = await api('/api/projects', { method: 'POST', body: { name, config: cfg } });
    if (!project || !project.id) throw new Error('The server did not return the new project.');
    toast(hadMasks
      ? `Project “${project.name || name}” created. Hand-drawn regions are not copied; draw them again in its ROI tab.`
      : `Project “${project.name || name}” created`, 'success', { timeout: hadMasks ? 9000 : 4500 });
    navigate(`/projects/${encodeURIComponent(project.id)}`);
  }

  // ---- live updates
  function onFinished() {
    st.cache = {};
    st.viewerStore = {};
    refreshHealth();
    const r = st.run;
    if (r.status === 'cancelled') {
      toast('The run was cancelled.', 'warning');
      return;
    }
    const ok = r.status === 'succeeded';
    const target = ok ? 'summary' : 'progress';
    let close = () => {};
    const go = () => {
      close();
      if (st.destroyed) navigate(`/runs/${enc}/${target}`);
      else selectTab(target);
    };
    const message = ok
      ? `Analysis finished${numOrNull(r.duration_s) !== null ? ` in ${fmt.duration(r.duration_s)}` : ''}.`
      : `The analysis ${r.status === 'interrupted' ? 'was interrupted' : 'failed'}.`;
    close = toast(h('div', { class: 'run-toast' }, h('span', null, message),
      h('button', { class: 'btn small', type: 'button', onClick: go }, ok ? 'View summary' : 'Show details')),
    ok ? 'success' : 'error', { timeout: 12000 });
  }

  function applyRun(next) {
    if (!next || typeof next !== 'object') return;
    const prev = st.run;
    st.run = { ...next, config: next.config !== undefined ? next.config : prev.config };
    const finished = ACTIVE.has(prev.status) && !isActive();
    if (finished) onFinished();
    paintHeader();
    paintTabs();
    if (pane && pane.update) pane.update();
    else if (finished && RESULT_TABS.has(st.tab)) renderPane();
    syncTicker();
    if (finished) getManifest().then(() => { if (!st.destroyed) paintTabs(); }, () => {});
  }

  function schedulePoll(ms = POLL_MS) {
    clearTimeout(pollTimer);
    pollTimer = null;
    if (st.destroyed || st.gone || !isActive()) return;
    pollTimer = setTimeout(pollOnce, ms);
  }

  async function pollOnce() {
    pollTimer = null;
    const my = new AbortController();
    pollCtrl = my;
    try {
      const next = await api(`/api/runs/${enc}`, { signal: my.signal });
      if (st.destroyed) return;
      if (pollFailures) setNotice(null);
      pollFailures = 0;
      applyRun(next);
    } catch (err) {
      if (st.destroyed || (err && err.name === 'AbortError')) return;
      if (err.status === 404) {
        st.gone = true;
        setNotice(banner('danger', ['This run no longer exists; it was deleted. ', h('a', { href: '#/runs' }, 'Back to all runs')]));
        return;
      }
      pollFailures += 1;
      setNotice(banner('warning', `Lost contact with the server (${err.message}). Retrying…`));
    } finally {
      if (pollCtrl === my) pollCtrl = null;
    }
    schedulePoll(pollFailures ? Math.min(POLL_MS * 2 ** pollFailures, 15000) : POLL_MS);
  }

  function syncTicker() {
    if (isActive() && !ticker) {
      ticker = setInterval(() => {
        for (const el of container.querySelectorAll('[data-since]')) el.textContent = fmt.duration(sinceSeconds(el.dataset.since));
      }, 1000);
    } else if (!isActive() && ticker) {
      clearInterval(ticker);
      ticker = null;
    }
  }

  paintHeader();
  paintTabs(true);
  renderPane();
  syncTicker();
  schedulePoll();
  if (!isActive()) getManifest().then(() => { if (!st.destroyed) paintTabs(); }, () => {});

  return () => {
    st.destroyed = true;
    paneToken += 1;
    clearTimeout(pollTimer);
    if (pollCtrl) pollCtrl.abort();
    if (ticker) clearInterval(ticker);
    if (pane && pane.destroy) pane.destroy();
    pane = null;
  };
}
