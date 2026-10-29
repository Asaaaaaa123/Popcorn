// table-viewer.js — the Tables tab of a run: the list of result tables and data objects, and a server-paged table
// viewer (GET /api/runs/<id>/table) with search, sorting, paging and number formatting.
import { h, api, apiUrl, mount, icon, card, emptyState, errorPanel, debounce, fmt } from '../core.js';

const PAGE = 100;
const nfInt = new Intl.NumberFormat(undefined, { maximumFractionDigits: 0 });
const nfSig = new Intl.NumberFormat(undefined, { maximumSignificantDigits: 4 });

/** Integers with separators, other numbers to 4 significant digits, very small or very large ones in scientific notation. */
export function formatNumber(value) {
  if (value === null || value === undefined || value === '') return '—';
  const x = typeof value === 'number' ? value : Number(value);
  if (Number.isNaN(x)) return String(value);
  if (!Number.isFinite(x)) return x > 0 ? '∞' : '−∞';
  const a = Math.abs(x);
  if (Number.isInteger(x)) return a >= 1e15 ? x.toExponential(3) : nfInt.format(x);
  if (a < 1e-3 || a >= 1e15) return x.toExponential(2);
  return nfSig.format(x);
}

function asBool(v) {
  if (typeof v === 'boolean') return v;
  if (typeof v === 'number') return v !== 0;
  const s = String(v).trim().toLowerCase();
  if (s === 'true' || s === '1' || s === 'yes') return true;
  if (s === 'false' || s === '0' || s === 'no' || s === '') return false;
  return null;
}

/** <td> for one value of a column of the given dtype ("number" | "text" | "bool"). */
function cell(value, dtype) {
  if (value === null || value === undefined) return h('td', { class: dtype === 'number' ? 'num tv-null' : 'tv-null' }, '—');
  if (dtype === 'bool' || typeof value === 'boolean') {
    const b = asBool(value);
    if (b !== null) {
      return h('td', { class: b ? 'tv-bool tv-true' : 'tv-bool tv-false' },
        h('span', { 'aria-hidden': 'true' }, b ? '✓' : '–'), h('span', { class: 'sr-only' }, b ? 'yes' : 'no'));
    }
  }
  if (dtype === 'number' || typeof value === 'number') {
    const text = formatNumber(value);
    const raw = String(value);
    return h('td', { class: 'num', title: raw !== text ? raw : null }, text);
  }
  const text = String(value);
  return h('td', { class: text.length > 60 ? 'tv-long' : null }, text); // ids stay on one line, long lists wrap
}

/** URL of a file inside a run's results/ folder; download = ask the browser to save it. */
export function fileUrl(runId, file, download) {
  const path = String(file).split('/').map(encodeURIComponent).join('/');
  return apiUrl(`/api/runs/${encodeURIComponent(runId)}/files/${path}`, download ? { download: 1 } : null);
}

function downloadLink(runId, file, label = 'Download') {
  return h('a', { class: 'btn small', href: fileUrl(runId, file, true), download: String(file).split('/').pop() },
    icon('download', { size: 14 }), label);
}

/** One table, fetched a page at a time. state = {q, sort, desc, offset}, kept by the caller between visits. */
function tableView({ runId, table, state }) {
  let destroyed = false;
  let ctrl = null;
  let reqId = 0;

  const search = h('input', {
    class: 'input tv-search', type: 'search', placeholder: 'Search all columns', value: state.q || '',
    'aria-label': `Search ${table.title || table.file}`,
  });
  const spinner = h('span', { class: 'spinner tv-spinner', hidden: true, 'aria-hidden': 'true' });
  const counter = h('div', { class: 'tv-count', role: 'status' });
  let total = 0;
  const pageBtn = (iconName, double, title, target) => h('button', {
    class: double ? 'btn small icon-only tv-double' : 'btn small icon-only', type: 'button', title, 'aria-label': title,
    onClick: () => go(target()),
  }, icon(iconName, { size: 14 }), double ? icon(iconName, { size: 14 }) : null);
  const firstBtn = pageBtn('chevron-left', true, 'First page', () => 0);
  const prevBtn = pageBtn('chevron-left', false, 'Previous page', () => state.offset - PAGE);
  const nextBtn = pageBtn('chevron-right', false, 'Next page', () => state.offset + PAGE);
  const lastBtn = pageBtn('chevron-right', true, 'Last page', () => Math.floor(Math.max(0, total - 1) / PAGE) * PAGE);
  const wrap = h('div', { class: 'table-wrap tv-wrap' }, h('div', { class: 'loading' }, h('span', { class: 'spinner' }), 'Loading table…'));

  const el = card({
    cls: 'tv-card',
    title: table.title || table.file,
    subtitle: h('span', null, h('span', { class: 'mono' }, table.file),
      table.rows !== null && table.rows !== undefined ? ` · ${fmt.int(table.rows)} rows` : '',
      table.size !== null && table.size !== undefined ? ` · ${fmt.bytes(table.size)}` : ''),
    actions: downloadLink(runId, table.file),
    body: h('div', { class: 'stack' },
      table.description ? h('p', { class: 'muted tv-desc' }, table.description) : null,
      h('div', { class: 'tv-toolbar' },
        h('div', { class: 'tv-search-box' }, icon('search', { size: 15 }), search),
        spinner, h('span', { class: 'spacer' }), counter,
        h('div', { class: 'tv-pager', role: 'group', 'aria-label': 'Pages' }, firstBtn, prevBtn, nextBtn, lastBtn)),
      wrap),
  });

  const onSearch = debounce(() => {
    const q = search.value.trim();
    if (q === (state.q || '')) return;
    state.q = q;
    state.offset = 0;
    load();
  }, 300);
  search.addEventListener('input', onSearch);
  search.addEventListener('keydown', (e) => {
    if (e.key === 'Escape' && search.value) {
      e.preventDefault();
      search.value = '';
      onSearch.cancel();
      state.q = '';
      state.offset = 0;
      load();
    }
  });

  function go(offset) {
    const next = Math.max(0, Math.min(offset, Math.floor(Math.max(0, total - 1) / PAGE) * PAGE));
    if (next === state.offset) return;
    state.offset = next;
    load();
  }

  function sortBy(col) {
    if (state.sort === col) state.desc = !state.desc;
    else { state.sort = col; state.desc = false; }
    state.offset = 0;
    load();
  }

  function setBusy(busy) {
    spinner.hidden = !busy;
    wrap.classList.toggle('is-loading', busy);
  }

  function paintPager() {
    const atStart = state.offset <= 0;
    const atEnd = state.offset + PAGE >= total;
    firstBtn.disabled = atStart;
    prevBtn.disabled = atStart;
    nextBtn.disabled = atEnd;
    lastBtn.disabled = atEnd;
  }

  function paint(data) {
    const cols = Array.isArray(data.columns) ? data.columns : [];
    const dtypes = data.dtypes && typeof data.dtypes === 'object' ? data.dtypes : {};
    const rows = Array.isArray(data.rows) ? data.rows : [];
    total = Number.isFinite(data.total) ? data.total : rows.length;
    const offset = Number.isFinite(data.offset) ? data.offset : state.offset;
    const head = h('tr', null, h('th', { class: 'num tv-rownum', scope: 'col' }, h('span', { class: 'sr-only' }, 'Row')),
      cols.map((c) => {
        const active = state.sort === c;
        return h('th', {
          class: [dtypes[c] === 'number' ? 'num' : '', 'sortable', active ? 'tv-sorted' : ''], scope: 'col',
          'aria-sort': active ? (state.desc ? 'descending' : 'ascending') : 'none',
        }, h('button', { class: 'tv-sort', type: 'button', title: `Sort by ${c}`, onClick: () => sortBy(c) },
          h('span', null, c), icon(active && state.desc ? 'chevron-down' : 'chevron-up', { size: 13, cls: active ? 'tv-sort-icon on' : 'tv-sort-icon' })));
      }));
    const body = rows.length
      ? rows.map((r, i) => h('tr', null, h('td', { class: 'num tv-rownum' }, fmt.int(offset + i + 1)),
        cols.map((c, k) => cell(Array.isArray(r) ? r[k] : null, dtypes[c]))))
      : h('tr', null, h('td', { class: 'tv-norows', colSpan: cols.length + 1 },
        state.q ? `No rows match “${state.q}”.` : 'This table has no rows.'));
    const scrollLeft = wrap.scrollLeft;
    mount(wrap, h('table', { class: 'table compact tv-table' }, h('thead', null, head), h('tbody', null, body)));
    wrap.scrollLeft = scrollLeft;
    wrap.scrollTop = 0;
    counter.textContent = total
      ? `Rows ${fmt.int(offset + 1)}–${fmt.int(offset + rows.length)} of ${fmt.int(total)}${state.q ? ' matching' : ''}`
      : '';
    paintPager();
  }

  async function load() {
    const my = ++reqId;
    if (ctrl) ctrl.abort();
    ctrl = new AbortController();
    setBusy(true);
    try {
      const data = await api(`/api/runs/${encodeURIComponent(runId)}/table`, {
        signal: ctrl.signal,
        query: {
          file: table.file, offset: state.offset, limit: PAGE, q: state.q || null,
          sort: state.sort || null, desc: state.sort ? (state.desc ? 1 : 0) : null,
        },
      });
      if (destroyed || my !== reqId) return;
      const n = data && Number.isFinite(data.total) ? data.total : 0;
      if (state.offset > 0 && state.offset >= n) {
        // the remembered page is past the end (fewer matches now): show the last page instead
        state.offset = Math.floor(Math.max(0, n - 1) / PAGE) * PAGE;
        load();
        return;
      }
      paint(data || {});
      setBusy(false);
    } catch (err) {
      if (destroyed || my !== reqId || (err && err.name === 'AbortError')) return;
      setBusy(false);
      total = 0;
      counter.textContent = '';
      paintPager();
      mount(wrap, h('div', { class: 'tv-error' }, errorPanel(err, { retry: load })));
    }
  }

  load();
  return {
    el,
    destroy() {
      destroyed = true;
      onSearch.cancel();
      if (ctrl) ctrl.abort();
    },
  };
}

/**
 * tableBrowser({runId, manifest, store}) -> {el, destroy()}
 * manifest: GET /api/runs/<id>/manifest. store: {selected, filter, states: {file: {q, sort, desc, offset}}}.
 */
export function tableBrowser({ runId, manifest, store }) {
  if (!store.states) store.states = {};
  const tables = Array.isArray(manifest.tables) ? manifest.tables.filter((t) => t && t.file) : [];
  const objects = Array.isArray(manifest.objects) ? manifest.objects.filter((o) => o && o.file) : [];
  const other = Array.isArray(manifest.other) ? manifest.other.filter((o) => o && o.file) : [];
  let viewer = null;

  const listEl = h('div', { class: 'tv-list' });
  const main = h('div', { class: 'tv-main' });
  const filter = tables.length > 8 ? h('input', {
    class: 'input tv-filter', type: 'search', placeholder: `Filter ${tables.length} tables`, value: store.filter || '',
    'aria-label': 'Filter tables',
  }) : null;

  function paintList() {
    const q = (store.filter || '').toLowerCase();
    const shown = tables.filter((t) => !q || `${t.title || ''} ${t.file} ${t.description || ''}`.toLowerCase().includes(q));
    mount(listEl, shown.length ? shown.map((t) => {
      const active = t.file === store.selected;
      return h('button', {
        class: active ? 'tv-item active' : 'tv-item', type: 'button', 'aria-pressed': active ? 'true' : 'false',
        title: t.description || t.file,
        onClick: () => select(t.file),
      },
      h('span', { class: 'tv-item-title' }, t.title || t.file),
      h('span', { class: 'tv-item-file' }, t.file),
      h('span', { class: 'tv-item-meta' },
        t.rows !== null && t.rows !== undefined ? `${fmt.int(t.rows)} rows` : 'rows unknown',
        t.size !== null && t.size !== undefined ? ` · ${fmt.bytes(t.size)}` : ''),
      t.description ? h('span', { class: 'tv-item-desc' }, t.description) : null);
    }) : h('div', { class: 'tv-none muted' }, 'No table matches the filter.'));
  }

  function select(file) {
    if (file === store.selected && viewer) return;
    store.selected = file;
    paintList();
    if (viewer) viewer.destroy();
    viewer = null;
    const t = tables.find((x) => x.file === file);
    if (!t) {
      mount(main, card({ body: emptyState({ icon: 'table', title: 'No table selected', text: 'Choose a table on the left.' }) }));
      return;
    }
    if (!store.states[file]) store.states[file] = { q: '', sort: null, desc: false, offset: 0 };
    viewer = tableView({ runId, table: t, state: store.states[file] });
    mount(main, viewer.el);
  }

  if (filter) {
    filter.addEventListener('input', () => { store.filter = filter.value.trim(); paintList(); });
  }

  const objectsCard = objects.length ? card({
    title: 'Data objects', subtitle: 'AnnData and other files for further analysis',
    body: h('ul', { class: 'tv-files' }, objects.map((o) => h('li', null,
      h('div', { class: 'tv-file-main' },
        h('div', { class: 'tv-item-title' }, o.title || o.file),
        h('div', { class: 'tv-item-file' }, o.file, o.size !== null && o.size !== undefined ? ` · ${fmt.bytes(o.size)}` : ''),
        o.description ? h('div', { class: 'tv-item-desc' }, o.description) : null),
      downloadLink(runId, o.file)))),
  }) : null;

  const otherBlock = other.length ? h('details', { class: 'card tv-other' },
    h('summary', null, `Other files (${other.length})`),
    h('ul', { class: 'tv-files' }, other.map((o) => h('li', null,
      h('div', { class: 'tv-file-main' }, h('div', { class: 'tv-item-file' }, o.file,
        o.size !== null && o.size !== undefined ? ` · ${fmt.bytes(o.size)}` : '')),
      downloadLink(runId, o.file))))) : null;

  const side = h('div', { class: 'tv-side' },
    tables.length ? card({ title: 'Tables', subtitle: `${tables.length} CSV file${tables.length === 1 ? '' : 's'}`, flush: true, body: h('div', null, filter ? h('div', { class: 'tv-filter-box' }, filter) : null, listEl) }) : null,
    objectsCard, otherBlock);

  const el = h('div', { class: tables.length ? 'tv-layout' : 'tv-layout no-tables' }, side, main);

  if (tables.length) {
    paintList();
    select(tables.some((t) => t.file === store.selected) ? store.selected : tables[0].file);
  } else {
    mount(main, card({ body: emptyState({ icon: 'table', title: 'No tables', text: 'This run wrote no CSV tables.' }) }));
  }

  return {
    el,
    destroy() {
      if (viewer) viewer.destroy();
      viewer = null;
    },
  };
}
