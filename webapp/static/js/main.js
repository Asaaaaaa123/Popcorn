// main.js — application shell: sidebar navigation, environment status and the route table.
import { h, icon, mount, route, startRouter, onRouteChange, refreshHealth, state, poll, svg, button } from './core.js';
import { renderProjects } from './views/projects.js';
import { renderProject } from './views/project.js';
import { renderRuns } from './views/runs.js';
import { renderRun } from './views/run.js';
import { renderHelp } from './views/help.js';

function brandMark() {
  // a 3 x 3 grid of capture bins on the accent tile, the "positive" ones solid
  const on = [1, 0, 1, 0, 1, 1, 1, 1, 0];
  const s = svg('svg', { viewBox: '0 0 28 28', class: 'brand-mark', 'aria-hidden': 'true' },
    svg('defs', null, svg('linearGradient', { id: 'brand-grad', x1: 0, y1: 0, x2: 1, y2: 1 },
      svg('stop', { offset: '0', 'stop-color': '#4d8bff' }), svg('stop', { offset: '1', 'stop-color': '#2353c9' }))),
    svg('rect', { width: 28, height: 28, rx: 7, fill: 'url(#brand-grad)' }));
  on.forEach((v, i) => {
    s.append(svg('rect', { x: 6 + (i % 3) * 6, y: 6 + Math.floor(i / 3) * 6, width: 4.2, height: 4.2, rx: 1, fill: '#ffffff', opacity: v ? 0.95 : 0.3 }));
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
  const envSub = h('div', { class: 'env-sub' });

  const banner = h('div', { id: 'global-banner' });
  const view = h('main', { id: 'view' });
  mount(app,
    h('aside', { class: 'sidebar' },
      h('div', { class: 'brand' }, brandMark(),
        h('div', null, h('div', { class: 'brand-name' }, 'Bone Spatial Studio'), h('div', { class: 'brand-sub' }, 'Visium HD bone lesion analysis'))),
      h('div', { class: 'nav-label' }, 'Workspace'),
      h('nav', { class: 'nav' }, links.map((l) => l.el)),
      h('div', { class: 'sidebar-foot' },
        h('div', { class: 'env-card' }, h('div', { class: 'env-line' }, envDot, envText), envSub))),
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

const viewRoot = buildShell();
poll(async () => { await refreshHealth(); return true; }, 5000);
startRouter(viewRoot);

// keep `state` reachable from the browser console for debugging
window.bss = { state, button };
