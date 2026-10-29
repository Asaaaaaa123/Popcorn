// projects.js — the project list: create a project or the demo, and one card per project (data root, samples, last
// run) with Open, Duplicate and Delete.
import { h, api, mount, icon, button, statusBadge, loading, emptyState, errorPanel, toast, modal, confirmDialog, promptDialog, fmt, navigate } from '../core.js';
import { hasUnsavedChanges, forgetDraft } from './project.js';

const projectPath = (id) => `/projects/${encodeURIComponent(id)}/data`;

const STEPS = [
  { icon: 'sliders', title: 'Configure', text: 'Point the project at the Space Ranger outputs of your Visium HD sections, adjust the gene signatures and parameters, and draw a region of interest on the H&E image.' },
  { icon: 'play', title: 'Run', text: 'The bone-lesion analysis runs in the background, one at a time, with live progress and a full log.' },
  { icon: 'chart', title: 'Explore', text: 'Browse the figures, the tables and the per-bin spatial maps of every run, and download them.' },
];

export async function renderProjects(container) {
  let alive = true;
  let projects = [];
  const body = h('div', null, loading('Loading projects…'));
  const filter = h('input', { class: 'input projects-filter', type: 'search', placeholder: 'Filter projects', 'aria-label': 'Filter projects', autocomplete: 'off' });
  filter.addEventListener('input', () => paint());

  const actions = () => [
    button({ label: 'Try the demo', icon: 'flask', onClick: () => { tryDemo(); } }),
    button({ label: 'New project', icon: 'plus', kind: 'primary', onClick: () => { newProject().catch((e) => toast(e.message, 'error')); } }),
  ];

  mount(container, h('div', { class: 'page' },
    h('div', { class: 'page-header' },
      h('div', { class: 'titles' }, h('h1', null, 'Projects'),
        h('div', { class: 'subtitle' }, 'A project holds the sections, gene signatures, parameters and drawn regions of one analysis.')),
      h('div', { class: 'actions' }, actions())),
    body));

  async function newProject() {
    const name = await promptDialog({
      title: 'New project', label: 'Name', placeholder: 'e.g. Tibia lesions, batch 2', confirmLabel: 'Create project',
      validate: (v) => (v ? '' : 'The project needs a name'),
    });
    if (!name) return;
    const p = await api('/api/projects', { method: 'POST', body: { name } });
    navigate(projectPath(p.id));
  }

  async function tryDemo() {
    const started = Date.now();
    const elapsed = h('span', { class: 'muted' }, '0 s');
    let background = false;
    const tick = setInterval(() => { elapsed.textContent = fmt.duration((Date.now() - started) / 1000); }, 1000);
    const wait = modal({
      title: 'Preparing the demo', dismissable: false,
      body: h('div', { class: 'demo-wait' }, h('span', { class: 'spinner' }), h('div', null,
        h('p', null, 'Generating two small synthetic tibiae in the Space Ranger layout and creating a demo project with fast settings.'),
        h('p', { class: 'muted' }, 'The first time this takes about a minute; afterwards the data is reused. Elapsed: ', elapsed))),
      actions: [{ label: 'Continue in the background', onClick: () => { background = true; } }],
      onClose: () => clearInterval(tick),
    });
    try {
      const res = await api('/api/demo', { method: 'POST' });
      const p = res && res.project;
      if (!p || !p.id) throw new Error('The server did not return the demo project.');
      wait.close();
      if (background) {
        toast(h('span', null, 'The demo project is ready. ', h('a', { href: `#${projectPath(p.id)}` }, 'Open it')), 'success', { timeout: 15000 });
        if (alive) load();
      } else {
        navigate(projectPath(p.id));
      }
    } catch (e) {
      wait.close();
      toast(`The demo could not be prepared: ${e.message}`, 'error');
    }
  }

  async function duplicate(p) {
    const copy = await api(`/api/projects/${encodeURIComponent(p.id)}/duplicate`, { method: 'POST' });
    toast(`Created “${copy.name}”`, 'success');
    if (alive) await load();
  }

  async function remove(p) {
    const ok = await confirmDialog({
      title: `Delete “${p.name || p.id}”?`,
      message: h('div', null,
        h('p', null, 'The project settings and its drawn regions are deleted.'),
        h('p', { class: 'muted' }, 'Its runs are kept: each run holds its own copy of the configuration and stays under Runs.')),
      confirmLabel: 'Delete project', danger: true,
    });
    if (!ok) return;
    await api(`/api/projects/${encodeURIComponent(p.id)}`, { method: 'DELETE' });
    forgetDraft(p.id);
    toast(`Deleted “${p.name || p.id}”`, 'success');
    if (alive) await load();
  }

  function projectCard(p) {
    const href = `#${projectPath(p.id)}`;
    const run = p.last_run;
    return h('article', { class: 'card project-card' },
      h('div', { class: 'pc-head' },
        h('a', { class: 'pc-name', href }, p.name || p.id),
        hasUnsavedChanges(p.id) ? h('span', { class: 'badge badge-unsaved', title: 'Edits made in this browser tab are not saved yet' }, 'Unsaved edits') : null),
      h('div', { class: 'pc-root', title: p.data_root || '' }, icon('folder', { size: 13 }), h('span', { class: 'mono truncate' }, p.data_root || '—')),
      h('dl', { class: 'pc-meta' },
        h('div', null, h('dt', null, 'Samples'), h('dd', null, fmt.int(p.n_samples))),
        h('div', null, h('dt', null, 'Updated'), h('dd', { title: fmt.date(p.updated) }, fmt.ago(p.updated))),
        h('div', null, h('dt', null, 'Last run'), h('dd', null, run
          ? h('a', { class: 'pc-run', href: `#/runs/${encodeURIComponent(run.id)}`, title: `${run.id} · ${fmt.date(run.created)}` },
            statusBadge(run.status), h('span', { class: 'muted small' }, fmt.ago(run.created)))
          : h('span', { class: 'muted' }, 'none')))),
      h('div', { class: 'pc-actions' },
        h('a', { class: 'btn small primary', href }, icon('chevron-right', { size: 14 }), 'Open'),
        button({ label: 'Duplicate', icon: 'copy', small: true, onClick: () => duplicate(p) }),
        h('span', { class: 'spacer' }),
        button({ label: 'Delete', icon: 'trash', small: true, kind: 'danger', onClick: () => remove(p) })));
  }

  function intro() {
    return h('div', { class: 'stack-lg' },
      emptyState({
        icon: 'folder', title: 'No projects yet',
        text: 'Bone Spatial Studio runs the bone-lesion analysis of Visium HD sections. Create a project for your own data, or try the demo with two small synthetic sections.',
        action: h('div', { class: 'row', style: 'justify-content:center' }, actions()),
      }),
      h('div', { class: 'grid-3 intro-steps' }, STEPS.map((s, i) => h('div', { class: 'card intro-step' },
        h('div', { class: 'intro-icon' }, icon(s.icon, { size: 20 })),
        h('h3', null, `${i + 1}. ${s.title}`),
        h('p', { class: 'muted' }, s.text)))));
  }

  function paint() {
    if (!projects.length) { mount(body, intro()); return; }
    const q = filter.value.trim().toLowerCase();
    const shown = projects.filter((p) => !q || `${p.name} ${p.data_root} ${p.id}`.toLowerCase().includes(q));
    mount(body,
      projects.length > 6 ? h('div', { class: 'row projects-tools' }, filter, h('span', { class: 'muted small' }, `${shown.length} of ${projects.length}`)) : null,
      shown.length ? h('div', { class: 'project-grid' }, shown.map(projectCard)) : h('p', { class: 'muted' }, 'No project matches the filter.'));
  }

  async function load() {
    try {
      const list = await api('/api/projects');
      if (!alive) return;
      projects = Array.isArray(list) ? list : [];
      paint();
    } catch (e) {
      if (!alive) return;
      mount(body, errorPanel(e, { retry: () => { mount(body, loading('Loading projects…')); return load(); } }));
    }
  }

  load();
  return () => { alive = false; };
}
