// runs.js — the Runs page: every analysis run, newest first, with its status, live progress and actions.
// The list keeps itself current, so it has no Refresh button: it polls GET /api/runs every 4 s while a run is queued
// or running (paused while the page is hidden), refreshes when the window gets the focus back or the page becomes
// visible, and the shell's health events wake it up when a run is started elsewhere. A row opens its run; its other
// actions are in the row's "Run actions" menu (shown on hover or focus), and Cancel stays visible on an active run.
import {
  h, api, mount, button, menuButton, icon, card, statusBadge, emptyState, errorPanel, banner, loading, fmt, navigate,
  confirmDialog, toast, refreshHealth,
} from '../core.js';

const ACTIVE = new Set(['queued', 'running']);
const POLL_MS = 4000;
const WAKE_GAP_MS = 1500; // focus and visibilitychange often come together: one refresh answers both
const STATUS_FILTERS = [
  { id: '', label: 'All statuses', test: () => true },
  { id: 'active', label: 'Queued or running', test: (r) => ACTIVE.has(r.status) },
  { id: 'succeeded', label: 'Succeeded', test: (r) => r.status === 'succeeded' },
  { id: 'failed', label: 'Failed or interrupted', test: (r) => r.status === 'failed' || r.status === 'interrupted' },
  { id: 'cancelled', label: 'Cancelled', test: (r) => r.status === 'cancelled' },
];

const runHref = (id) => `#/runs/${encodeURIComponent(id)}`;
const projectHref = (id) => `#/projects/${encodeURIComponent(id)}`;

function elapsedSince(iso) {
  const t = iso ? new Date(iso).getTime() : NaN;
  return Number.isFinite(t) ? Math.max(0, (Date.now() - t) / 1000) : null;
}

function durationCell(r) {
  if (r.status === 'running') {
    const s = elapsedSince(r.started);
    return s === null ? '—' : h('span', { title: 'Running for' }, fmt.duration(s));
  }
  if (r.status === 'queued') return h('span', { class: 'faint' }, '—');
  return fmt.duration(r.duration_s);
}

function progressCell(r, queuePos) {
  const p = r.progress || {};
  if (r.status === 'running') {
    const total = Number(p.total) || 0;
    const frac = Math.min(1, Math.max(0, Number(p.fraction) || 0));
    const text = total
      ? `Step ${Math.min(Number(p.index) || 0, total)} of ${total}${p.label ? ` — ${p.label}` : ''}`
      : 'Starting…';
    return h('div', { class: 'runs-progress' },
      h('div', {
        class: total ? 'progress' : 'progress indeterminate', role: 'progressbar', 'aria-label': 'Progress',
        'aria-valuemin': '0', 'aria-valuemax': '100', 'aria-valuenow': total ? String(Math.round(frac * 100)) : null,
      }, h('div', { class: 'bar', style: { width: `${(frac * 100).toFixed(1)}%` } })),
      h('div', { class: 'runs-step truncate', title: text }, text));
  }
  if (r.status === 'queued') {
    return h('span', { class: 'muted' }, queuePos === 1 ? 'Next in the queue' : `Waiting (${queuePos} in the queue)`);
  }
  if ((r.status === 'failed' || r.status === 'interrupted') && r.error) {
    const first = String(r.error).split('\n').find((l) => l.trim()) || '';
    return h('span', { class: 'runs-error truncate', title: String(r.error) }, first);
  }
  const steps = Array.isArray(r.steps) ? r.steps : [];
  if (!steps.length) return h('span', { class: 'faint' }, '—');
  const done = steps.filter((s) => s.status === 'done').length;
  const skipped = steps.filter((s) => s.status === 'skipped').length;
  return h('span', { class: 'muted' }, `${done} step${done === 1 ? '' : 's'} done${skipped ? `, ${skipped} skipped` : ''}`);
}

export async function renderRuns(container) {
  let destroyed = false;
  let runs = null;
  let timer = null;
  let ctrl = null;
  let failures = 0;
  let refreshError = null;
  let lastLoad = 0;
  let repaintTimer = null;
  const filters = { project: '', status: '' };

  const subtitle = h('div', { class: 'subtitle' }, 'Every analysis started from a project, newest first.');
  const notice = h('div');
  const body = h('div', null, loading('Loading runs…'));
  mount(container, h('div', { class: 'page runs-page' },
    h('div', { class: 'page-header' },
      h('div', { class: 'titles' }, h('h1', null, 'Runs'), subtitle),
      h('div', { class: 'actions' },
        h('a', { class: 'btn primary', href: '#/projects' }, icon('play', { size: 16 }), 'Start a run'))),
    notice, body));

  async function cancelRun(r) {
    const running = r.status === 'running';
    const ok = await confirmDialog({
      title: 'Cancel this run?',
      message: running
        ? 'The analysis process is stopped. Outputs written so far stay in the run folder, but the run is incomplete.'
        : 'The run is taken out of the queue and will not start.',
      confirmLabel: 'Cancel run', danger: true,
    });
    if (!ok) return;
    await api(`/api/runs/${encodeURIComponent(r.id)}/cancel`, { method: 'POST' });
    toast(running ? 'Run cancelled' : 'Run removed from the queue', 'success');
    refreshHealth();
    await refresh();
  }

  async function deleteRun(r) {
    const ok = await confirmDialog({
      title: 'Delete this run?',
      message: `Run ${r.id}${r.label ? ` (“${r.label}”)` : ''} and all of its results, figures and log are removed from the workspace. This cannot be undone.`,
      confirmLabel: 'Delete run', danger: true,
    });
    if (!ok) return;
    await api(`/api/runs/${encodeURIComponent(r.id)}`, { method: 'DELETE' });
    toast('Run deleted', 'success');
    await refresh();
  }

  // the row's secondary actions; read when the menu opens, from the run as this row shows it
  function rowMenuItems(r) {
    const active = ACTIVE.has(r.status);
    return [
      r.project_id ? { label: 'Open project', icon: 'folder', hint: r.project_name || r.project_id, href: projectHref(r.project_id) } : null,
      'separator',
      {
        label: 'Delete run', icon: 'trash', danger: true, disabled: active,
        hint: active ? 'Cancel the run before deleting it' : null, onClick: () => deleteRun(r),
      },
    ];
  }

  function row(r, queuePos) {
    const active = ACTIVE.has(r.status);
    // an active run keeps a visible Cancel; everything else is in the row's menu, revealed on hover or focus
    const cancel = active
      ? button({ label: 'Cancel', icon: 'stop', small: true, kind: 'ghost danger', title: 'Cancel this run', onClick: () => cancelRun(r) })
      : null;
    if (cancel) cancel.dataset.action = 'cancel';
    const more = menuButton({ small: true, title: 'Run actions', items: () => rowMenuItems(r) });
    more.setAttribute('aria-label', `Actions for run ${r.label || r.id}`);
    more.dataset.action = 'more';
    return h('tr', {
      class: 'clickable', dataset: { id: r.id },
      onClick: (e) => {
        if (e.target instanceof Element && e.target.closest('a, button')) return;
        navigate(`/runs/${encodeURIComponent(r.id)}`);
      },
    },
    h('td', null, statusBadge(r.status)),
    h('td', { class: 'runs-name' },
      h('a', { href: runHref(r.id), dataset: { action: 'open' } }, r.label || 'Untitled run'),
      h('div', { class: 'mono faint small' }, r.id)),
    h('td', null, r.project_id
      ? h('a', { href: projectHref(r.project_id), dataset: { action: 'project' } }, r.project_name || r.project_id)
      : h('span', { class: 'muted' }, r.project_name || '—')),
    h('td', { class: 'nowrap', title: fmt.date(r.created) }, fmt.ago(r.created)),
    h('td', { class: 'nowrap num' }, durationCell(r)),
    h('td', { class: 'runs-progress-cell' }, progressCell(r, queuePos)),
    h('td', { class: 'runs-actions-cell' }, h('div', { class: 'runs-actions' }, cancel, h('span', { class: 'row-actions' }, more))));
  }

  // filter bar and table are built once and updated in place, so a poll never steals focus or closes a select
  const projectSel = h('select', { class: 'select', 'aria-label': 'Filter by project', onChange: () => { filters.project = projectSel.value; paint(); } });
  const statusSel = h('select', { class: 'select', 'aria-label': 'Filter by status', onChange: () => { filters.status = statusSel.value; paint(); } },
    STATUS_FILTERS.map((f) => h('option', { value: f.id }, f.label)));
  const tbody = h('tbody');
  const table = h('div', { class: 'table-wrap borderless' }, h('table', { class: 'table runs-table' },
    h('thead', null, h('tr', null,
      h('th', { scope: 'col' }, 'Status'), h('th', { scope: 'col' }, 'Run'), h('th', { scope: 'col' }, 'Project'),
      h('th', { scope: 'col' }, 'Created'), h('th', { scope: 'col', class: 'num' }, 'Duration'),
      h('th', { scope: 'col' }, 'Progress'), h('th', { scope: 'col' }, h('span', { class: 'sr-only' }, 'Actions')))),
    tbody));
  const noMatch = emptyState({
    icon: 'search', title: 'No runs match the filters',
    action: button({ label: 'Clear filters', small: true, onClick: () => { filters.project = ''; filters.status = ''; statusSel.value = ''; paint(); } }),
  });
  const tableHost = h('div');
  const listCard = card({ flush: true, body: h('div', null, h('div', { class: 'toolbar runs-filters' }, projectSel, statusSel), tableHost) });
  let projectKey = null;

  function paintFilters(projects) {
    const key = JSON.stringify(projects);
    if (key !== projectKey) {
      projectKey = key;
      mount(projectSel, h('option', { value: '' }, 'All projects'), projects.map(([pid, name]) => h('option', { value: pid }, name)));
    }
    projectSel.hidden = projects.length < 2;
    projectSel.value = filters.project;
    statusSel.value = filters.status;
  }

  /** Replace only the rows whose content changed; a changed set or order of runs rebuilds the body. */
  function paintRows(rows) {
    const focused = document.activeElement;
    const focusTr = focused && tbody.contains(focused) ? focused.closest('tr') : null;
    const focusKey = focusTr ? [focusTr.dataset.id, focused.dataset.action] : null;
    const old = [...tbody.children];
    let held = false;
    if (old.length === rows.length && old.every((tr, i) => tr.dataset.id === rows[i].dataset.id)) {
      old.forEach((tr, i) => {
        // a row whose menu is open keeps its button (the open menu belongs to it) and is updated a moment later
        if (tr.querySelector('.menu-trigger[aria-expanded="true"]')) { held = true; return; }
        if (tr.outerHTML !== rows[i].outerHTML) tr.replaceWith(rows[i]);
      });
    } else if (tbody.querySelector('.menu-trigger[aria-expanded="true"]')) {
      held = true; // a changed set of runs is rebuilt once the open menu has closed, not under it
    } else {
      mount(tbody, rows);
    }
    clearTimeout(repaintTimer);
    repaintTimer = held ? setTimeout(() => { repaintTimer = null; paint(); }, 1000) : null;
    if (focusKey && !tbody.contains(document.activeElement)) {
      const tr = [...tbody.children].find((t) => t.dataset.id === focusKey[0]);
      const target = tr && focusKey[1] ? tr.querySelector(`[data-action="${focusKey[1]}"]`) : null;
      if (target) target.focus();
    }
  }

  function paint() {
    if (destroyed || !runs) return;
    const active = runs.filter((r) => ACTIVE.has(r.status)).length;
    subtitle.textContent = runs.length
      ? `${runs.length} run${runs.length === 1 ? '' : 's'}${active ? ` · ${active} queued or running` : ''}. Every analysis started from a project, newest first.`
      : 'Every analysis started from a project, newest first.';
    mount(notice, refreshError ? banner('warning', `Could not refresh the list (${refreshError.message}). Retrying…`) : null);

    if (!runs.length) {
      mount(body, card({
        body: emptyState({
          icon: 'activity', title: 'No runs yet',
          text: 'A run is one analysis of a project: its configuration is frozen, the pipeline runs in the background and the results stay here. Open a project, check its data and start a run from its Run tab.',
          // secondary: "Start a run" in the header is the page's one primary action and leads to the same place
          action: h('a', { class: 'btn', href: '#/projects' }, icon('folder', { size: 16 }), 'Go to projects'),
        }),
      }));
      return;
    }
    if (!body.contains(listCard)) mount(body, listCard);

    const projects = [...new Map(runs.filter((r) => r.project_id).map((r) => [r.project_id, r.project_name || r.project_id])).entries()]
      .sort((a, b) => String(a[1]).localeCompare(String(b[1])));
    if (filters.project && !projects.some(([pid]) => pid === filters.project)) filters.project = '';
    paintFilters(projects);
    const statusTest = (STATUS_FILTERS.find((f) => f.id === filters.status) || STATUS_FILTERS[0]).test;
    const shown = runs.filter((r) => (!filters.project || r.project_id === filters.project) && statusTest(r));
    // queue position: queued runs start oldest first
    const queued = runs.filter((r) => r.status === 'queued').sort((a, b) => String(a.created).localeCompare(String(b.created)));
    const queuePos = new Map(queued.map((r, i) => [r.id, i + 1]));

    if (!shown.length) {
      mount(tableHost, noMatch);
      return;
    }
    if (!tableHost.contains(table)) mount(tableHost, table);
    paintRows(shown.map((r) => row(r, queuePos.get(r.id))));
  }

  // poll while a run is queued or running, and retry after a failed refresh (the notice says "Retrying…"); a hidden
  // page waits for onWake instead
  function schedule() {
    clearTimeout(timer);
    timer = null;
    if (destroyed || !runs || document.hidden) return;
    if (!failures && !runs.some((r) => ACTIVE.has(r.status))) return;
    timer = setTimeout(refresh, failures ? Math.min(POLL_MS * 2 ** failures, 30000) : POLL_MS);
  }

  async function refresh() {
    if (destroyed) return;
    clearTimeout(timer);
    timer = null;
    if (ctrl) ctrl.abort();
    const my = new AbortController();
    ctrl = my;
    lastLoad = Date.now();
    try {
      const list = await api('/api/runs', { signal: my.signal });
      if (destroyed || ctrl !== my) return;
      runs = Array.isArray(list) ? list.filter((r) => r && r.id) : [];
      failures = 0;
      refreshError = null;
      paint();
    } catch (err) {
      if (destroyed || ctrl !== my || (err && err.name === 'AbortError')) return;
      failures += 1;
      if (!runs) {
        mount(body, errorPanel(err, { retry: () => { mount(body, loading('Loading runs…')); return refresh(); } }));
        return;
      }
      refreshError = err;
      paint();
    } finally {
      if (ctrl === my) ctrl = null;
    }
    schedule();
  }

  // a run started or finished elsewhere (another tab, the project page): the shell reports it every few seconds
  const onHealth = (ev) => {
    const hl = ev.detail || {};
    if (!runs || timer || ctrl) return;
    const serverActive = !!hl.active_run || Number(hl.queued) > 0;
    if (serverActive !== runs.some((r) => ACTIVE.has(r.status))) refresh();
  };
  window.addEventListener('health', onHealth);

  // back in the window or tab (runs may have been started, finished or deleted meanwhile): refresh at once
  const onWake = () => {
    if (destroyed || document.hidden || Date.now() - lastLoad < WAKE_GAP_MS) return;
    refresh();
  };
  window.addEventListener('focus', onWake);
  document.addEventListener('visibilitychange', onWake);

  await refresh();

  return () => {
    destroyed = true;
    clearTimeout(timer);
    clearTimeout(repaintTimer);
    if (ctrl) ctrl.abort();
    window.removeEventListener('health', onHealth);
    window.removeEventListener('focus', onWake);
    document.removeEventListener('visibilitychange', onWake);
  };
}
