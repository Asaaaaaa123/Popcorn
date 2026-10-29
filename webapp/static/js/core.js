// core.js — shared helpers for every view of Bone Spatial Studio:
// API client, DOM builder, icons, toasts, dialogs, small components, formatting, polling and the hash router.

// ------------------------------------------------------------------------------------------------ API client
export class ApiError extends Error {
  constructor(message, status = 0, detail = null) {
    super(message);
    this.name = 'ApiError';
    this.status = status;
    this.detail = detail;
  }
}

/** URL of an API path with query parameters (undefined, null and '' values are dropped). */
export function apiUrl(path, query) {
  const u = new URL(path, window.location.origin);
  if (query) {
    for (const [k, v] of Object.entries(query)) {
      if (v === undefined || v === null || v === '') continue;
      u.searchParams.set(k, String(v));
    }
  }
  return u.pathname + u.search;
}

function errorFrom(status, data) {
  let message = `Request failed (HTTP ${status})`;
  let detail = data;
  if (data && typeof data === 'object' && 'detail' in data) {
    detail = data.detail;
    if (typeof data.detail === 'string') message = data.detail;
    else if (Array.isArray(data.detail)) {
      message = data.detail
        .map((d) => (d && d.msg ? `${(d.loc || []).filter((x) => x !== 'body').join('.')}: ${d.msg}` : String(d)))
        .join('; ');
    } else if (data.detail && typeof data.detail === 'object') {
      message = data.detail.message || message;
    }
  } else if (typeof data === 'string' && data.trim()) {
    message = data.trim().slice(0, 400);
  }
  return new ApiError(message, status, detail);
}

/**
 * JSON request to the backend. Paths are absolute ('/api/projects').
 * body: object (sent as JSON) or FormData. Resolves to the parsed JSON (or text); rejects with ApiError.
 */
export async function api(path, { method = 'GET', body, query, signal } = {}) {
  const opts = { method, headers: { Accept: 'application/json' }, signal, credentials: 'same-origin' };
  if (body !== undefined) {
    if (body instanceof FormData) opts.body = body;
    else {
      opts.headers['Content-Type'] = 'application/json';
      opts.body = JSON.stringify(body);
    }
  }
  let res;
  try {
    res = await fetch(apiUrl(path, query), opts);
  } catch (e) {
    if (e && e.name === 'AbortError') throw e;
    throw new ApiError('Cannot reach the server. Is it still running?', 0);
  }
  const text = await res.text();
  let data = null;
  if (text) {
    try { data = JSON.parse(text); } catch { data = text; }
  }
  if (!res.ok) throw errorFrom(res.status, data);
  return data;
}

// ------------------------------------------------------------------------------------------------ DOM
const SVG_NS = 'http://www.w3.org/2000/svg';
const PROPS = new Set(['value', 'checked', 'disabled', 'selected', 'multiple', 'readOnly', 'indeterminate',
  'hidden', 'title', 'placeholder', 'type', 'min', 'max', 'step', 'rows', 'tabIndex', 'href', 'src', 'alt',
  'id', 'name', 'colSpan', 'rowSpan', 'spellcheck', 'autocomplete', 'download', 'target', 'rel', 'width', 'height']);

function setAttrs(el, attrs) {
  if (!attrs) return;
  for (const [k, v] of Object.entries(attrs)) {
    if (v === undefined || v === null || v === false) continue;
    if (k === 'class' || k === 'className') {
      el.setAttribute('class', Array.isArray(v) ? v.filter(Boolean).join(' ') : String(v));
    } else if (k === 'style') {
      if (typeof v === 'string') el.style.cssText = v;
      else Object.assign(el.style, v);
    } else if (k === 'dataset') {
      Object.assign(el.dataset, v);
    } else if (k.startsWith('on') && typeof v === 'function') {
      el.addEventListener(k.slice(2).toLowerCase(), v);
    } else if (k === 'html') {
      el.innerHTML = v; // trusted strings only
    } else if (PROPS.has(k) && !(el instanceof SVGElement)) {
      el[k] = v;
    } else {
      el.setAttribute(k, v === true ? '' : String(v));
    }
  }
}

function appendChildren(el, children) {
  for (const c of children.flat(Infinity)) {
    if (c === null || c === undefined || c === false || c === true) continue;
    el.append(c instanceof Node ? c : document.createTextNode(String(c)));
  }
}

/** h('div', {class: 'card', onClick: fn}, 'text', childNode, [more, nodes]) */
export function h(tag, attrs, ...children) {
  const el = document.createElement(tag);
  setAttrs(el, attrs);
  appendChildren(el, children);
  return el;
}

/** SVG element: svg('circle', {cx: 4, cy: 4, r: 2}) */
export function svg(tag, attrs, ...children) {
  const el = document.createElementNS(SVG_NS, tag);
  setAttrs(el, attrs);
  appendChildren(el, children);
  return el;
}

export function clear(el) {
  while (el.firstChild) el.removeChild(el.firstChild);
  return el;
}

/** Replace the children of `el`. */
export function mount(el, ...children) {
  clear(el);
  appendChildren(el, children);
  return el;
}

// ------------------------------------------------------------------------------------------------ icons
const ICONS = {
  folder: '<path d="M3 7a2 2 0 0 1 2-2h4l2 2h8a2 2 0 0 1 2 2v8a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2z"/>',
  play: '<path d="M7 4.5v15l12-7.5z"/>',
  stop: '<rect x="6" y="6" width="12" height="12" rx="1.5"/>',
  refresh: '<path d="M21 12a9 9 0 1 1-2.6-6.4L21 8"/><path d="M21 3v5h-5"/>',
  download: '<path d="M12 3v12"/><path d="m7 10 5 5 5-5"/><path d="M5 21h14"/>',
  upload: '<path d="M12 21V9"/><path d="m7 14 5-5 5 5"/><path d="M5 3h14"/>',
  trash: '<path d="M3 6h18"/><path d="M8 6V4h8v2"/><path d="M6 6l1 14h10l1-14"/><path d="M10 11v6M14 11v6"/>',
  plus: '<path d="M12 5v14M5 12h14"/>',
  edit: '<path d="M12 20h9"/><path d="M16.5 3.5a2.1 2.1 0 0 1 3 3L7 19l-4 1 1-4z"/>',
  copy: '<rect x="9" y="9" width="12" height="12" rx="2"/><path d="M5 15H4a1 1 0 0 1-1-1V4a1 1 0 0 1 1-1h10a1 1 0 0 1 1 1v1"/>',
  'chevron-right': '<path d="m9 18 6-6-6-6"/>',
  'chevron-left': '<path d="m15 18-6-6 6-6"/>',
  'chevron-down': '<path d="m6 9 6 6 6-6"/>',
  'chevron-up': '<path d="m18 15-6-6-6 6"/>',
  check: '<path d="M20 6 9 17l-5-5"/>',
  x: '<path d="M18 6 6 18M6 6l12 12"/>',
  alert: '<path d="M10.3 3.9 1.8 18a2 2 0 0 0 1.7 3h17a2 2 0 0 0 1.7-3L13.7 3.9a2 2 0 0 0-3.4 0z"/><path d="M12 9v4M12 17h.01"/>',
  info: '<circle cx="12" cy="12" r="10"/><path d="M12 16v-4M12 8h.01"/>',
  image: '<rect x="3" y="3" width="18" height="18" rx="2"/><circle cx="8.5" cy="8.5" r="1.5"/><path d="m21 15-5-5L5 21"/>',
  table: '<rect x="3" y="3" width="18" height="18" rx="2"/><path d="M3 9h18M3 15h18M9 3v18"/>',
  map: '<path d="M1 6v16l7-4 8 4 7-4V2l-7 4-8-4-7 4z"/><path d="M8 2v16M16 6v16"/>',
  sliders: '<path d="M4 21v-7M4 10V3M12 21v-9M12 8V3M20 21v-5M20 12V3M1 14h6M9 8h6M17 16h6"/>',
  home: '<path d="m3 10 9-7 9 7v10a1 1 0 0 1-1 1h-5v-7H9v7H4a1 1 0 0 1-1-1z"/>',
  list: '<path d="M8 6h13M8 12h13M8 18h13M3 6h.01M3 12h.01M3 18h.01"/>',
  help: '<circle cx="12" cy="12" r="10"/><path d="M9.1 9a3 3 0 0 1 5.8 1c0 2-3 3-3 3M12 17h.01"/>',
  external: '<path d="M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6"/><path d="M15 3h6v6M10 14 21 3"/>',
  sun: '<circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/>',
  moon: '<path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/>',
  layers: '<path d="m12 2 10 5-10 5L2 7z"/><path d="m2 17 10 5 10-5M2 12l10 5 10-5"/>',
  undo: '<path d="M3 7v6h6"/><path d="M21 17a9 9 0 0 0-15-6.7L3 13"/>',
  'zoom-in': '<circle cx="11" cy="11" r="8"/><path d="m21 21-4.3-4.3M11 8v6M8 11h6"/>',
  'zoom-out': '<circle cx="11" cy="11" r="8"/><path d="m21 21-4.3-4.3M8 11h6"/>',
  maximize: '<path d="M8 3H5a2 2 0 0 0-2 2v3M21 8V5a2 2 0 0 0-2-2h-3M3 16v3a2 2 0 0 0 2 2h3M16 21h3a2 2 0 0 0 2-2v-3"/>',
  file: '<path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/><path d="M14 2v6h6"/>',
  activity: '<path d="M22 12h-4l-3 9L9 3l-3 9H2"/>',
  search: '<circle cx="11" cy="11" r="8"/><path d="m21 21-4.3-4.3"/>',
  clock: '<circle cx="12" cy="12" r="10"/><path d="M12 6v6l4 2"/>',
  flask: '<path d="M9 3h6M10 3v6L4.5 18.5A2 2 0 0 0 6.2 21h11.6a2 2 0 0 0 1.7-2.5L14 9V3"/><path d="M7 15h10"/>',
  polygon: '<path d="M12 2.5 21 9l-3.5 11.5h-11L3 9z"/>',
  eye: '<path d="M1 12s4-8 11-8 11 8 11 8-4 8-11 8S1 12 1 12z"/><circle cx="12" cy="12" r="3"/>',
  save: '<path d="M19 21H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h11l5 5v11a2 2 0 0 1-2 2z"/><path d="M17 21v-8H7v8M7 3v5h8"/>',
  more: '<circle cx="12" cy="12" r="1.2"/><circle cx="19" cy="12" r="1.2"/><circle cx="5" cy="12" r="1.2"/>',
  terminal: '<path d="m4 17 6-6-6-6M12 19h8"/>',
  chart: '<path d="M12 20V10M18 20V4M6 20v-4"/>',
  dna: '<path d="M4 2c0 5 16 5 16 10S4 17 4 22"/><path d="M20 2c0 5-16 5-16 10s16 5 16 10"/><path d="M7 5h10M7 19h10"/>',
  grid: '<rect x="3" y="3" width="7" height="7" rx="1"/><rect x="14" y="3" width="7" height="7" rx="1"/><rect x="3" y="14" width="7" height="7" rx="1"/><rect x="14" y="14" width="7" height="7" rx="1"/>',
};

/** Inline SVG icon (stroke style). icon('play', {size: 16}) */
export function icon(name, { size = 16, cls = '', title } = {}) {
  const el = document.createElementNS(SVG_NS, 'svg');
  el.setAttribute('viewBox', '0 0 24 24');
  el.setAttribute('width', size);
  el.setAttribute('height', size);
  el.setAttribute('fill', 'none');
  el.setAttribute('stroke', 'currentColor');
  el.setAttribute('stroke-width', '2');
  el.setAttribute('stroke-linecap', 'round');
  el.setAttribute('stroke-linejoin', 'round');
  el.setAttribute('class', `icon ${cls}`.trim());
  el.setAttribute('aria-hidden', title ? 'false' : 'true');
  el.innerHTML = (title ? `<title>${escapeHtml(title)}</title>` : '') + (ICONS[name] || ICONS.info);
  return el;
}

export function escapeHtml(s) {
  return String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
}

// ------------------------------------------------------------------------------------------------ small components
/**
 * Button. onClick may be async: the button shows a spinner and is disabled until it settles;
 * a rejected promise is reported with toast().
 */
export function button({ label = '', icon: iconName, kind = '', small = false, onClick, title, disabled = false, type = 'button', cls = '' } = {}) {
  const classes = ['btn', kind, small ? 'small' : '', !label && iconName ? 'icon-only' : '', cls].filter(Boolean).join(' ');
  const b = h('button', { class: classes, type, title: title || (label ? undefined : iconName), disabled });
  const paint = (busy) => {
    mount(b, busy ? h('span', { class: 'spinner' }) : (iconName ? icon(iconName, { size: small ? 14 : 16 }) : null), label || null);
  };
  paint(false);
  if (onClick) {
    b.addEventListener('click', async (ev) => {
      if (b.disabled) return;
      const ret = onClick(ev);
      if (ret && typeof ret.then === 'function') {
        b.disabled = true;
        paint(true);
        try { await ret; } catch (err) { toast(err.message || String(err), 'error'); console.error(err); } finally {
          if (b.isConnected) { b.disabled = disabled; paint(false); }
        }
      }
    });
  }
  return b;
}

export function card({ title, subtitle, actions, body, footer, cls = '', flush = false } = {}) {
  return h('section', { class: `card ${cls}`.trim() },
    (title || actions) ? h('div', { class: 'card-header' },
      h('div', { class: 'card-titles' }, title ? h('h2', null, title) : null, subtitle ? h('div', { class: 'card-sub' }, subtitle) : null),
      actions ? h('div', { class: 'row' }, actions) : null) : null,
    body !== undefined ? h('div', { class: flush ? 'card-body flush' : 'card-body' }, body) : null,
    footer ? h('div', { class: 'card-footer' }, footer) : null);
}

/** Labeled form field around any control. */
export function field({ label, help, control, error, inline = false, forId } = {}) {
  return h('div', { class: inline ? 'field inline' : 'field' },
    label ? h('label', { for: forId }, label) : null,
    control,
    error ? h('div', { class: 'error-text' }, error) : null,
    help ? h('div', { class: 'help' }, help) : null);
}

export function switchToggle({ checked = false, onChange, label, disabled = false, title } = {}) {
  const input = h('input', { type: 'checkbox', checked, disabled });
  if (onChange) input.addEventListener('change', () => onChange(input.checked));
  return h('label', { class: 'switch', title }, input, h('span', { class: 'track' }), label ? h('span', null, label) : null);
}

/** Tab strip. items: [{id, label, icon, badge, badgeKind: 'warn'|'bad'}]; onSelect(id). */
export function tabs({ items, active, onSelect }) {
  const nav = h('nav', { class: 'tabs', role: 'tablist' });
  for (const it of items) {
    nav.append(h('button', {
      class: it.id === active ? 'tab active' : 'tab', role: 'tab', 'aria-selected': it.id === active ? 'true' : 'false',
      onClick: () => onSelect && onSelect(it.id),
    }, it.icon ? icon(it.icon, { size: 15 }) : null, it.label,
    it.badge !== undefined && it.badge !== null && it.badge !== '' ? h('span', { class: `tab-badge ${it.badgeKind || ''}`.trim() }, String(it.badge)) : null));
  }
  return nav;
}

const STATUS_LABEL = {
  queued: 'Queued', running: 'Running', succeeded: 'Succeeded', failed: 'Failed', cancelled: 'Cancelled',
  interrupted: 'Interrupted', pending: 'Pending', skipped: 'Skipped', done: 'Done',
};
export function statusBadge(status) {
  const s = String(status || 'pending');
  return h('span', { class: `badge badge-${s === 'done' ? 'succeeded' : s}` }, h('span', { class: 'dot' }), STATUS_LABEL[s] || s);
}

export function loading(text = 'Loading…') {
  return h('div', { class: 'loading' }, h('span', { class: 'spinner' }), text);
}

export function emptyState({ icon: iconName = 'info', title, text, action } = {}) {
  return h('div', { class: 'empty' },
    h('div', { class: 'empty-icon' }, icon(iconName, { size: 34 })),
    title ? h('h3', null, title) : null,
    text ? h('p', null, text) : null,
    action || null);
}

export function errorPanel(err, { retry } = {}) {
  const msg = (err && err.message) || String(err);
  return h('div', { class: 'banner danger' }, icon('alert', { size: 18 }),
    h('div', { class: 'banner-body' }, h('strong', null, 'Something went wrong. '), msg),
    retry ? button({ label: 'Retry', icon: 'refresh', small: true, onClick: retry }) : null);
}

/** Banner: kind = info | warning | danger | success; body = string or nodes. */
export function banner(kind, body, { iconName } = {}) {
  const ic = iconName || ({ info: 'info', warning: 'alert', danger: 'alert', success: 'check' }[kind] || 'info');
  return h('div', { class: `banner ${kind}` }, icon(ic, { size: 18 }), h('div', { class: 'banner-body' }, body));
}

// ------------------------------------------------------------------------------------------------ toasts and dialogs
export function toast(message, kind = 'info', { timeout = 4500 } = {}) {
  let root = document.getElementById('toasts');
  if (!root) { root = h('div', { id: 'toasts' }); document.body.append(root); }
  const close = () => { el.remove(); };
  const el = h('div', { class: `toast ${kind}`, role: kind === 'error' ? 'alert' : 'status' },
    h('div', { class: 'toast-body' }, message),
    h('button', { 'aria-label': 'Dismiss', onClick: close }, icon('x', { size: 14 })));
  root.append(el);
  if (timeout > 0) setTimeout(close, kind === 'error' ? Math.max(timeout, 8000) : timeout);
  return close;
}

/**
 * Modal dialog. actions: [{label, kind, icon, onClick}] — onClick may be async; return false to keep it open.
 * size: '' | 'wide' | 'full'. Resolves nothing; use the returned close().
 */
let modalCount = 0;

export function modal({ title, body, actions = [], size = '', onClose, dismissable = true } = {}) {
  const root = document.getElementById('modal-root') || document.body;
  const opener = document.activeElement;
  let closed = false;
  const close = () => {
    if (closed) return;
    closed = true;
    backdrop.remove();
    document.removeEventListener('keydown', onKey, true);
    modalCount -= 1;
    const app = document.getElementById('app');
    if (app && modalCount === 0) app.inert = false;
    if (onClose) onClose();
    // after the opener's own async click handler has re-enabled it
    setTimeout(() => { if (opener && opener.isConnected && typeof opener.focus === 'function') opener.focus(); }, 0);
  };
  const focusables = () => [...dialog.querySelectorAll('button, [href], input, select, textarea, [tabindex]:not([tabindex="-1"])')]
    .filter((el) => !el.disabled && el.offsetParent !== null);
  const onKey = (e) => {
    if (root.lastElementChild !== backdrop) return;                             // a dialog opened on top of this one
    if (e.key === 'Escape' && dismissable) { e.stopPropagation(); close(); return; }
    if (e.key !== 'Tab') return;
    const items = focusables();                                                 // keep Tab inside the dialog
    if (!items.length) { e.preventDefault(); dialog.focus(); return; }
    const first = items[0];
    const last = items[items.length - 1];
    if (e.shiftKey && (document.activeElement === first || !dialog.contains(document.activeElement))) {
      e.preventDefault(); last.focus();
    } else if (!e.shiftKey && (document.activeElement === last || !dialog.contains(document.activeElement))) {
      e.preventDefault(); first.focus();
    }
  };
  const footer = actions.length ? h('div', { class: 'modal-footer' }, actions.map((a) => button({
    label: a.label, icon: a.icon, kind: a.kind || '',
    onClick: async () => {
      const r = a.onClick ? await a.onClick({ close }) : undefined;
      if (r !== false) close();
    },
  }))) : null;
  const titleId = `modal-title-${Math.random().toString(36).slice(2, 9)}`;
  const dialog = h('div', { class: `modal ${size}`.trim(), role: 'dialog', 'aria-modal': 'true', 'aria-labelledby': titleId, tabindex: '-1' },
    h('div', { class: 'modal-header' }, h('h2', { id: titleId }, title || ''),
      dismissable ? button({ icon: 'x', kind: 'ghost', small: true, title: 'Close', onClick: close }) : null),
    h('div', { class: 'modal-body' }, body),
    footer);
  const backdrop = h('div', { class: 'modal-backdrop', onMousedown: (e) => { if (e.target === backdrop && dismissable) close(); } }, dialog);
  root.append(backdrop);
  modalCount += 1;
  const app = document.getElementById('app');
  if (app) app.inert = true;                                                    // the page behind is not reachable
  document.addEventListener('keydown', onKey, true);
  // first field, else the confirming button (primary or danger), else the dialog itself
  const target = dialog.querySelector('.modal-body input, .modal-body select, .modal-body textarea')
    || dialog.querySelector('.modal-footer .btn:last-child') || dialog;
  setTimeout(() => { if (!closed) target.focus(); }, 30);
  return { close, el: dialog };
}

export function confirmDialog({ title = 'Are you sure?', message = '', confirmLabel = 'Confirm', danger = false } = {}) {
  return new Promise((resolve) => {
    let answered = false;
    modal({
      title, body: typeof message === 'string' ? h('p', null, message) : message,
      onClose: () => { if (!answered) resolve(false); },
      actions: [
        { label: 'Cancel', onClick: () => { answered = true; resolve(false); } },
        { label: confirmLabel, kind: danger ? 'danger solid' : 'primary', onClick: () => { answered = true; resolve(true); } },
      ],
    });
  });
}

export function promptDialog({ title = '', label = '', value = '', placeholder = '', confirmLabel = 'OK', validate } = {}) {
  return new Promise((resolve) => {
    let answered = false;
    const input = h('input', { class: 'input', value, placeholder });
    const err = h('div', { class: 'error-text' });
    const submit = () => {
      const v = input.value.trim();
      const msg = validate ? validate(v) : (v ? '' : 'Required');
      if (msg) { err.textContent = msg; input.classList.add('invalid'); return false; }
      answered = true;
      resolve(v);
      return true;
    };
    const m = modal({
      title, body: h('div', { class: 'field' }, label ? h('label', null, label) : null, input, err),
      onClose: () => { if (!answered) resolve(null); },
      actions: [
        { label: 'Cancel', onClick: () => { answered = true; resolve(null); } },
        { label: confirmLabel, kind: 'primary', onClick: () => submit() },
      ],
    });
    input.addEventListener('keydown', (e) => { if (e.key === 'Enter' && submit()) m.close(); });
  });
}

// ------------------------------------------------------------------------------------------------ formatting
const nf0 = new Intl.NumberFormat(undefined, { maximumFractionDigits: 0 });
export const fmt = {
  int(x) { return x === null || x === undefined || Number.isNaN(Number(x)) ? '—' : nf0.format(Math.round(Number(x))); },
  num(x, digits = 2) {
    if (x === null || x === undefined || x === '' || Number.isNaN(Number(x))) return '—';
    const v = Number(x);
    if (!Number.isFinite(v)) return v > 0 ? '∞' : '−∞';
    if (v !== 0 && Math.abs(v) < 10 ** -digits) return v.toExponential(1);
    return v.toLocaleString(undefined, { maximumFractionDigits: digits, minimumFractionDigits: 0 });
  },
  pct(x, digits = 2) { return x === null || x === undefined || Number.isNaN(Number(x)) ? '—' : `${fmt.num(x, digits)}%`; },
  bytes(n) {
    if (n === null || n === undefined) return '—';
    const u = ['B', 'KB', 'MB', 'GB', 'TB'];
    let i = 0; let v = Number(n);
    while (v >= 1024 && i < u.length - 1) { v /= 1024; i += 1; }
    return `${v.toFixed(i && v < 10 ? 1 : 0)} ${u[i]}`;
  },
  duration(seconds) {
    if (seconds === null || seconds === undefined || Number.isNaN(Number(seconds))) return '—';
    let s = Math.max(0, Math.round(Number(seconds)));
    const hh = Math.floor(s / 3600); s -= hh * 3600;
    const mm = Math.floor(s / 60); s -= mm * 60;
    if (hh) return `${hh} h ${mm} min`;
    if (mm) return `${mm} min ${s} s`;
    return `${s} s`;
  },
  date(iso) {
    if (!iso) return '—';
    const d = new Date(iso);
    return Number.isNaN(d.getTime()) ? String(iso) : d.toLocaleString(undefined, { dateStyle: 'medium', timeStyle: 'short' });
  },
  ago(iso) {
    if (!iso) return '—';
    const d = new Date(iso);
    if (Number.isNaN(d.getTime())) return String(iso);
    const s = (Date.now() - d.getTime()) / 1000;
    if (s < 45) return 'just now';
    if (s < 3600) return `${Math.round(s / 60)} min ago`;
    if (s < 86400) return `${Math.round(s / 3600)} h ago`;
    if (s < 86400 * 14) return `${Math.round(s / 86400)} d ago`;
    return fmt.date(iso);
  },
};

// ------------------------------------------------------------------------------------------------ utilities
export function debounce(fn, ms = 250) {
  let t = null;
  const d = (...args) => { clearTimeout(t); t = setTimeout(() => fn(...args), ms); };
  d.cancel = () => clearTimeout(t);
  return d;
}

export const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

/**
 * Call fn() now and then every `ms` after each call settles, until fn returns false or stop() is called.
 * Errors are logged and polling continues. Returns {stop}.
 */
export function poll(fn, ms = 1500) {
  let stopped = false;
  let timer = null;
  const tick = async () => {
    if (stopped) return;
    let again = true;
    try { again = (await fn()) !== false; } catch (e) { console.error(e); }
    if (!stopped && again) timer = setTimeout(tick, ms);
  };
  tick();
  return { stop() { stopped = true; clearTimeout(timer); } };
}

export async function copyText(text) {
  try {
    await navigator.clipboard.writeText(text);
  } catch {
    const ta = h('textarea', { style: 'position:fixed;left:-9999px' }, text);
    document.body.append(ta);
    ta.select();
    document.execCommand('copy');
    ta.remove();
  }
  toast('Copied to clipboard', 'success', { timeout: 1800 });
}

/** Deep copy of plain JSON data. */
export const clone = (x) => (x === undefined ? undefined : JSON.parse(JSON.stringify(x)));

// ------------------------------------------------------------------------------------------------ shared state
export const state = { health: null };

export async function refreshHealth() {
  try { state.health = await api('/api/health'); } catch (e) { state.health = { ok: false, error: e.message }; }
  window.dispatchEvent(new CustomEvent('health', { detail: state.health }));
  return state.health;
}

// ------------------------------------------------------------------------------------------------ router
// Hash routes: '#/projects/abc/data'. route('/projects/:id/:tab?', handler). A handler receives (container, params);
// it renders into `container` (a fresh element per navigation, so a slow stale render can never overwrite a newer
// view) and may return a cleanup function that runs before the next navigation.
const routes = [];
const routeListeners = [];
let viewRoot = null;
let navToken = 0;
let currentCleanup = null;

export function route(pattern, handler) {
  const keys = [];
  const source = pattern.replace(/\/:([A-Za-z_]+)(\?)?/g, (_, key, optional) => {
    keys.push(key);
    return optional ? '(?:/([^/]+))?' : '/([^/]+)';
  });
  routes.push({ re: new RegExp(`^${source}/?$`), keys, handler });
}

export function currentPath() {
  const p = window.location.hash.replace(/^#/, '');
  return p || '/';
}

/** Navigate to a path ('/runs/abc'). Navigating to the current path re-renders it. */
export function navigate(path, { replace = false } = {}) {
  const target = `#${path}`;
  if (window.location.hash === target) { resolveRoute(); return; }
  if (replace) { history.replaceState(null, '', target); resolveRoute(); } else { window.location.hash = target; }
}

/** Change the address bar without re-rendering (e.g. switching tabs inside a view that keeps its own state). */
export function setPathSilently(path) {
  history.replaceState(null, '', `#${path}`);
  routeListeners.forEach((fn) => fn(path));
}

export function onRouteChange(fn) { routeListeners.push(fn); }

export function startRouter(root) {
  viewRoot = root;
  window.addEventListener('hashchange', resolveRoute);
  resolveRoute();
}

export async function resolveRoute() {
  const path = currentPath();
  const token = ++navToken;
  if (typeof currentCleanup === 'function') {
    try { currentCleanup(); } catch (e) { console.error(e); }
  }
  currentCleanup = null;
  routeListeners.forEach((fn) => fn(path));
  const container = h('div', { class: 'view' }, loading());
  mount(viewRoot, container);
  let match = null;
  for (const r of routes) {
    const m = r.re.exec(path);
    if (m) {
      const params = {};
      r.keys.forEach((k, i) => { params[k] = m[i + 1] === undefined ? undefined : decodeURIComponent(m[i + 1]); });
      match = { handler: r.handler, params };
      break;
    }
  }
  if (!match) {
    mount(container, h('div', { class: 'page' }, emptyState({
      icon: 'alert', title: 'Page not found', text: `Nothing lives at ${path}.`,
      action: h('a', { class: 'btn primary', href: '#/projects' }, 'Go to projects'),
    })));
    return;
  }
  try {
    const cleanup = await match.handler(container, match.params);
    if (token !== navToken) { if (typeof cleanup === 'function') cleanup(); return; }
    currentCleanup = cleanup;
  } catch (err) {
    console.error(err);
    if (token === navToken) mount(container, h('div', { class: 'page' }, errorPanel(err, { retry: resolveRoute })));
  }
}
