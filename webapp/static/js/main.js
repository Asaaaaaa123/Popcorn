// main.js — application shell: sidebar navigation, environment status, theme switch and the route table.
import { h, icon, mount, route, startRouter, onRouteChange, refreshHealth, state, poll, svg, button } from './core.js';
import { renderProjects } from './views/projects.js';
import { renderProject } from './views/project.js';
import { renderRuns } from './views/runs.js';
import { renderRun } from './views/run.js';
import { renderHelp } from './views/help.js';

const THEME_KEY = 'bss-theme';
const THEMES = ['auto', 'light', 'dark'];

function storedTheme() {
  try { return localStorage.getItem(THEME_KEY) || 'auto'; } catch { return 'auto'; }
}

function applyTheme(theme) {
  if (theme === 'light' || theme === 'dark') document.documentElement.dataset.theme = theme;
  else delete document.documentElement.dataset.theme;
}

function brandMark() {
  // a 4 x 4 grid of capture bins, a few of them "positive"
  const colors = ['#9fd3ea', '#ffffff', '#ffffff', '#f4a6a6', '#ffffff', '#f6c667', '#9fd3ea', '#ffffff',
    '#ffffff', '#9fd3ea', '#f4a6a6', '#ffffff', '#f6c667', '#ffffff', '#ffffff', '#9fd3ea'];
  const s = svg('svg', { viewBox: '0 0 32 32', class: 'brand-mark', 'aria-hidden': 'true' },
    svg('rect', { width: 32, height: 32, rx: 8, fill: '#2b6c8c' }));
  colors.forEach((c, i) => {
    s.append(svg('rect', { x: 5 + (i % 4) * 6, y: 5 + Math.floor(i / 4) * 6, width: 4, height: 4, rx: 1.2, fill: c, opacity: c === '#ffffff' ? 0.35 : 1 }));
  });
  return s;
}

const NAV = [
  { href: '#/projects', label: 'Projects', icon: 'folder', match: /^\/(projects)?(\/|$)/ },
  { href: '#/runs', label: 'Runs', icon: 'activity', match: /^\/runs/, counter: true },
  { href: '#/help', label: 'Guide', icon: 'help', match: /^\/help/ },
];

function buildShell() {
  const app = document.getElementById('app');
  const links = NAV.map((n) => {
    const a = h('a', { href: n.href }, icon(n.icon, { size: 17 }), h('span', null, n.label));
    if (n.counter) a.append(h('span', { class: 'nav-count hidden', id: 'nav-run-count' }));
    return { ...n, el: a };
  });
  const envDot = h('span', { class: 'dot' });
  const envText = h('span', null, 'Checking environment…');
  const envSub = h('div', { class: 'faint', style: 'font-size:11px;margin-bottom:8px' });
  const themeBtn = h('button', { class: 'btn ghost small', title: 'Switch between automatic, light and dark theme' });
  const paintTheme = () => {
    const t = storedTheme();
    mount(themeBtn, icon(t === 'dark' ? 'moon' : t === 'light' ? 'sun' : 'layers', { size: 14 }),
      h('span', null, `Theme: ${t}`));
  };
  themeBtn.addEventListener('click', () => {
    const next = THEMES[(THEMES.indexOf(storedTheme()) + 1) % THEMES.length];
    try { localStorage.setItem(THEME_KEY, next); } catch { /* storage unavailable */ }
    applyTheme(next);
    paintTheme();
  });
  paintTheme();

  const banner = h('div', { id: 'global-banner' });
  const view = h('main', { id: 'view' });
  mount(app,
    h('aside', { class: 'sidebar' },
      h('div', { class: 'brand' }, brandMark(),
        h('div', null, h('div', { class: 'brand-name' }, 'Bone Spatial Studio'), h('div', { class: 'brand-sub' }, 'Visium HD bone lesion analysis'))),
      h('nav', { class: 'nav' }, links.map((l) => l.el)),
      h('div', { class: 'sidebar-foot' }, h('div', { class: 'env-line' }, envDot, envText), envSub, themeBtn)),
    h('div', { style: 'flex:1 1 auto;min-width:0;display:flex;flex-direction:column' }, banner, view));

  onRouteChange((path) => {
    for (const l of links) l.el.classList.toggle('active', l.match.test(path));
  });

  window.addEventListener('health', (ev) => {
    const hl = ev.detail || {};
    const missing = hl.missing_required || [];
    if (!hl.ok && hl.error) {
      envDot.className = 'dot bad';
      envText.textContent = 'Server unreachable';
      envSub.textContent = '';
    } else if (missing.length) {
      envDot.className = 'dot bad';
      envText.textContent = 'Packages missing';
      envSub.textContent = `Python ${hl.python || '?'}`;
    } else {
      envDot.className = hl.active_run ? 'dot ok pulse' : 'dot ok';
      envText.textContent = hl.active_run ? 'Analysis running' : 'Ready';
      envSub.textContent = `Python ${hl.python || '?'} · scanpy ${(hl.packages || {}).scanpy || '?'}`;
    }
    const count = (hl.active_run ? 1 : 0) + (hl.queued || 0);
    const badge = document.getElementById('nav-run-count');
    if (badge) {
      badge.textContent = String(count);
      badge.classList.toggle('hidden', count === 0);
    }
    mount(banner, missing.length ? h('div', { class: 'page', style: 'padding-bottom:0' },
      h('div', { class: 'banner danger' }, icon('alert', { size: 18 }), h('div', { class: 'banner-body' },
        h('strong', null, 'Missing Python packages: '), missing.join(', '),
        '. Analyses cannot start until they are installed: ', h('code', null, 'pip install -r requirements.txt')))) : null);
  });
  return view;
}

route('/', async (container) => renderProjects(container, {}));
route('/projects', async (container) => renderProjects(container, {}));
route('/projects/:id/:tab?', async (container, p) => renderProject(container, { id: p.id, tab: p.tab || 'data' }));
route('/runs', async (container) => renderRuns(container, {}));
route('/runs/:id/:tab?', async (container, p) => renderRun(container, { id: p.id, tab: p.tab || null }));
route('/help/:topic?', async (container, p) => renderHelp(container, { topic: p.topic || null }));

applyTheme(storedTheme());
const viewRoot = buildShell();
poll(async () => { await refreshHealth(); return true; }, 5000);
startRouter(viewRoot);

// keep `state` reachable from the browser console for debugging
window.bss = { state, button };
