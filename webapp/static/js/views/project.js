// project.js — the project editor. A working copy (draft) per project id lives in this module, so switching tabs or
// leaving and coming back never loses edits; Save writes it (PUT), Discard drops it. Tabs: data | parameters |
// signatures | roi | run.
import {
  h, api, mount, icon, button, card, field, tabs, switchToggle, loading, emptyState, errorPanel, banner, toast,
  confirmDialog, fmt, clone, navigate, setPathSilently, resolveRoute, statusBadge, poll, state,
} from '../core.js';
import { openFolderPicker } from '../components/folder-picker.js';
import { sampleEditor, chooseDetectedSamples, applyDetectedSamples, reportPanel, binLabel } from '../components/sample-editor.js';
import { paramForm } from '../components/param-form.js';
import { signatureEditor, geneCoverage } from '../components/signature-editor.js';
import { roiEditor } from '../components/roi-editor.js';

const TABS = [
  { id: 'data', label: 'Data', icon: 'folder' },
  { id: 'parameters', label: 'Parameters', icon: 'sliders' },
  { id: 'signatures', label: 'Signatures', icon: 'dna' },
  { id: 'roi', label: 'Region of interest', icon: 'polygon' },
  { id: 'run', label: 'Run', icon: 'play' },
];
const SLOW_STEPS = ['moran', 'sensitivity', 'stability'];
const DEFAULT_BINS = ['002um', '008um', '016um'];
const VALIDATE_DELAY_MS = 800;
const RUNS_POLL_MS = 4000;

const drafts = new Map();
let defaultsPromise = null;
let guardInstalled = false;

const enc = encodeURIComponent;
const has = (o, k) => !!o && Object.prototype.hasOwnProperty.call(o, k);
const plural = (n, one, many = `${one}s`) => `${n} ${n === 1 ? one : many}`;
const uniq = (list) => [...new Set(list)];

// ------------------------------------------------------------------------------------------------ drafts
function normValidation(v) {
  return {
    errors: v && Array.isArray(v.errors) ? v.errors : [],
    warnings: v && Array.isArray(v.warnings) ? v.warnings : [],
  };
}

const snapshot = (name, config) => JSON.stringify({ name, config });

function setRois(d, list) {
  d.rois.clear();
  (Array.isArray(list) ? list : []).forEach((sid) => d.rois.add(sid));
}

/** The server's project becomes the saved state and the working copy. Report, UI state and regions are kept. */
function adopt(d, project) {
  d.name = project.name || '';
  d.config = clone(project.config || {});
  d.saved = { name: d.name, config: clone(d.config) };
  d.savedJson = snapshot(d.name, d.config);
  d.dirty = false;
  d.validation = normValidation(project.validation);
  d.created = project.created;
  d.updated = project.updated;
  d.remoteUpdated = project.updated;
  if (Array.isArray(project.rois)) setRois(d, project.rois);
  d.originalIds = new WeakMap((Array.isArray(d.config.samples) ? d.config.samples : [])
    .filter((s) => s && typeof s === 'object').map((s) => [s, s.id]));
}

function newDraft(project) {
  const d = {
    id: project.id, rois: new Set(), report: null, reportKey: null, reportAt: null, reportError: null, checking: false,
    saving: false, ui: { signature: null, runLabel: '', runSteps: {}, runDeterministic: null }, roi: {},
  };
  adopt(d, project);
  return d;
}

function roiUnsaved(d) {
  // only sections that still exist: drawings of a removed section can be neither saved nor discarded
  const ids = new Set((Array.isArray(d.config.samples) ? d.config.samples : []).map((s) => s && s.id));
  return Object.entries((d.roi && d.roi.entries) || {})
    .filter(([sid, en]) => ids.has(sid) && en && en.loaded && JSON.stringify(en.polys) !== en.saved).map(([sid]) => sid);
}

/** True while a project has unsaved settings or regions in this browser tab. */
export function hasUnsavedChanges(id) {
  const d = drafts.get(id);
  return !!d && (d.dirty || roiUnsaved(d).length > 0);
}

/** Drop the working copy of a deleted project. */
export function forgetDraft(id) {
  drafts.delete(id);
}

function installGuard() {
  if (guardInstalled) return;
  guardInstalled = true;
  window.addEventListener('beforeunload', (e) => {
    for (const d of drafts.values()) {
      if (d.dirty || roiUnsaved(d).length) {
        e.preventDefault();
        e.returnValue = '';
        return;
      }
    }
  });
}

function getDefaults() {
  if (!defaultsPromise) defaultsPromise = api('/api/defaults').catch((e) => { defaultsPromise = null; throw e; });
  return defaultsPromise;
}

/** What a data check depends on; a report made for another key is stale. */
function dataKey(cfg) {
  return JSON.stringify([cfg.data_root, (Array.isArray(cfg.samples) ? cfg.samples : []).map((s) => [s.id, s.path, s.microns_per_pixel]),
    cfg.bin_size, cfg.qc_bin_size]);
}

// ------------------------------------------------------------------------------------------------ view
export async function renderProject(container, { id, tab }) {
  installGuard();
  const active = TABS.some((t) => t.id === tab) ? tab : 'data';
  if (active !== tab) setPathSilently(`/projects/${enc(id)}/${active}`);
  mount(container, h('div', { class: 'page' }, loading('Loading the project…')));

  let defaults;
  let project;
  try {
    [defaults, project] = await Promise.all([getDefaults(), api(`/api/projects/${enc(id)}`)]);
  } catch (e) {
    if (e.status === 404) {
      drafts.delete(id);
      mount(container, h('div', { class: 'page' }, emptyState({
        icon: 'folder', title: 'Project not found', text: `There is no project “${id}”. It may have been deleted.`,
        action: h('a', { class: 'btn primary', href: '#/projects' }, 'All projects'),
      })));
    } else {
      mount(container, h('div', { class: 'page' }, errorPanel(e, { retry: resolveRoute })));
    }
    return null;
  }

  let d = drafts.get(id);
  if (!d) {
    d = newDraft(project);
    drafts.set(id, d);
  } else if (!d.dirty) {
    adopt(d, project);
  } else {
    d.remoteUpdated = project.updated;
    d.validation = normValidation(project.validation);
    if (Array.isArray(project.rois)) setRois(d, project.rois);
  }
  return projectView(container, d, defaults, active);
}

function projectView(container, d, defaults, initialTab) {
  let alive = true;
  let current = initialTab;
  let tabApi = {};
  let vOpen = false;
  let validateTimer = null;
  let validateAbort = null;
  const defCfg = (defaults && defaults.config) || {};

  // ---------------------------------------------------------------- header
  const nameInput = h('input', { class: 'title-input', value: d.name, 'aria-label': 'Project name', spellcheck: false, autocomplete: 'off', maxlength: 200 });
  const crumbName = h('span', { class: 'truncate' });
  const meta = h('div', { class: 'subtitle project-meta' });
  const vchips = h('div', { class: 'vchips' });
  const vpanelId = `vpanel-${d.id}`;
  const vpanel = h('div', { class: 'card vpanel hidden', id: vpanelId });
  const conflict = h('div');
  const tabsHost = h('div', { class: 'project-tabs' });
  const body = h('div', { class: 'project-body' });
  const saveBtn = button({ label: 'Save', icon: 'save', kind: 'primary', onClick: () => { save(); } });
  const discardBtn = button({ label: 'Discard', icon: 'undo', onClick: () => { discard(); } });
  const validateBtn = button({ label: 'Validate data', icon: 'check', onClick: () => { validateNow(true); } });
  saveBtn.title = 'Save the project settings (Ctrl+S)';
  discardBtn.title = 'Drop the unsaved changes';
  validateBtn.title = 'Check the sample folders and the genes against the data, with your current (unsaved) settings';
  const barSave = button({ label: 'Save', icon: 'save', kind: 'primary', small: true, onClick: () => { save(); } });
  const barDiscard = button({ label: 'Discard', small: true, onClick: () => { discard(); } });
  const unsavedBar = h('div', { class: 'unsaved-bar hidden', role: 'status' },
    icon('alert', { size: 16 }), h('strong', null, 'Unsaved changes'),
    h('span', { class: 'muted small' }, 'Saved settings are what runs use. Ctrl+S saves.'),
    h('span', { class: 'spacer' }), barDiscard, barSave);

  nameInput.addEventListener('input', () => { d.name = nameInput.value; refreshChrome(); });
  nameInput.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') { e.preventDefault(); nameInput.blur(); } else if (e.key === 'Escape') {
      d.name = d.saved.name;
      nameInput.value = d.name;
      refreshChrome();
      nameInput.blur();
    }
  });

  mount(container, h('div', { class: 'page project-page' },
    h('div', { class: 'breadcrumb' }, h('a', { href: '#/projects' }, 'Projects'), icon('chevron-right', { size: 12 }), crumbName),
    h('div', { class: 'page-header project-head' },
      h('div', { class: 'titles' }, nameInput, meta),
      h('div', { class: 'actions' }, vchips, validateBtn, discardBtn, saveBtn)),
    vpanel, conflict, tabsHost, body, unsavedBar));

  function setButton(btn, busy, iconName, label) {
    mount(btn, busy ? h('span', { class: 'spinner' }) : icon(iconName, { size: btn.classList.contains('small') ? 14 : 16 }), label);
  }

  function computeDirty() {
    d.dirty = snapshot(d.name, d.config) !== d.savedJson;
    return d.dirty;
  }

  let chipsSig = '';
  let tabsSig = '';

  function paintChips() {
    const { errors, warnings } = d.validation;
    const sig = JSON.stringify([errors, warnings, vOpen, d.dirty]);
    if (sig === chipsSig) return;
    chipsSig = sig;
    const scope = d.dirty ? 'Check of the saved settings; it is updated when you save.' : 'Check of the saved settings.';
    const chip = (kind, text) => h('button', {
      type: 'button', class: `vchip ${kind}`, 'aria-expanded': String(vOpen), 'aria-controls': vpanelId, title: scope,
      onClick: () => { vOpen = !vOpen; paintChips(); },
    }, icon('alert', { size: 13 }), text);
    mount(vchips,
      errors.length ? chip('bad', plural(errors.length, 'error')) : null,
      warnings.length ? chip('warn', plural(warnings.length, 'warning')) : null,
      !errors.length && !warnings.length ? h('span', { class: 'vchip ok', title: scope }, icon('check', { size: 13 }), 'Settings valid') : null);
    if (vOpen && (errors.length || warnings.length)) {
      mount(vpanel, h('div', { class: 'card-body stack' },
        h('div', { class: 'row between' }, h('h3', null, 'Configuration check'),
          h('span', { class: 'muted small' }, d.dirty ? 'Saved settings; your unsaved changes are checked when you save.' : 'Saved settings.'),
          button({ icon: 'x', kind: 'ghost', small: true, title: 'Close', onClick: () => { vOpen = false; paintChips(); } })),
        errors.length ? banner('danger', [h('strong', null, `${plural(errors.length, 'error')}: runs cannot start until they are fixed.`), h('ul', null, errors.map((t) => h('li', null, t)))]) : null,
        warnings.length ? banner('warning', [h('strong', null, plural(warnings.length, 'warning')), h('ul', null, warnings.map((t) => h('li', null, t)))]) : null));
      vpanel.classList.remove('hidden');
    } else {
      vOpen = false;
      mount(vpanel);
      vpanel.classList.add('hidden');
    }
  }

  function signaturesWithMissingGenes() {
    const cov = geneCoverage(d.report);
    if (!cov.known) return 0;
    return Object.values(d.config.signatures || {}).filter((s) => s && [...(s.extended || []), ...(s.core || [])].some((g) => cov.missing.has(g))).length;
  }

  function paintTabs() {
    const repErrors = d.report ? (d.report.errors || []).length : 0;
    const sigMissing = signaturesWithMissingGenes();
    const pendingRois = roiUnsaved(d).length;
    const cfgErrors = d.validation.errors.length;
    const sig = JSON.stringify([current, repErrors, sigMissing, pendingRois, cfgErrors]);
    if (sig === tabsSig) return;
    tabsSig = sig;
    const hadFocus = tabsHost.contains(document.activeElement);
    const badge = {
      data: repErrors ? { badge: repErrors, badgeKind: 'bad' } : {},
      signatures: sigMissing ? { badge: sigMissing, badgeKind: 'warn' } : {},
      roi: pendingRois ? { badge: pendingRois, badgeKind: 'warn' } : {},
      run: cfgErrors ? { badge: cfgErrors, badgeKind: 'bad' } : {},
    };
    mount(tabsHost, tabs({ items: TABS.map((t) => ({ ...t, ...(badge[t.id] || {}) })), active: current, onSelect: switchTab }));
    if (hadFocus) {
      const b = tabsHost.querySelector('.tab.active');
      if (b) b.focus();
    }
  }

  function paintConflict() {
    mount(conflict, d.dirty && d.remoteUpdated && d.updated && d.remoteUpdated !== d.updated
      ? banner('warning', h('span', null, `This project was saved elsewhere (${fmt.ago(d.remoteUpdated)}) after you started editing; saving overwrites that version. `,
        h('button', { type: 'button', class: 'linklike', onClick: () => { discard(); } }, 'Load the saved version')))
      : null);
  }

  function refreshChrome() {
    if (!alive) return;
    computeDirty();
    const named = d.name.trim().length > 0;
    nameInput.classList.toggle('invalid', !named);
    nameInput.title = named ? 'Click to rename' : 'The project needs a name';
    crumbName.textContent = d.name.trim() || d.id;
    mount(meta, h('span', { class: 'mono' }, d.id), ` · created ${fmt.date(d.created)} · saved ${fmt.ago(d.updated)}`);
    for (const b of [saveBtn, barSave]) b.disabled = !d.dirty || d.saving;
    for (const b of [discardBtn, barDiscard]) b.disabled = !d.dirty || d.saving;
    setButton(saveBtn, d.saving, 'save', d.saving ? 'Saving…' : 'Save');
    setButton(barSave, d.saving, 'save', d.saving ? 'Saving…' : 'Save');
    validateBtn.disabled = d.checking;
    setButton(validateBtn, d.checking, 'check', d.checking ? 'Checking…' : 'Validate data');
    unsavedBar.classList.toggle('hidden', !d.dirty);
    paintChips();
    paintTabs();
    paintConflict();
    if (tabApi.onChrome) tabApi.onChrome();
  }

  function changed({ data = false } = {}) {
    refreshChrome();
    if (data) {
      clearTimeout(validateTimer);
      validateTimer = setTimeout(() => validateNow(false), VALIDATE_DELAY_MS);
    }
  }

  function notifyReport() {
    if (!alive) return;
    refreshChrome();
    if (tabApi.onReport) tabApi.onReport();
  }

  // ---------------------------------------------------------------- data check, save, discard
  async function validateNow(manual) {
    clearTimeout(validateTimer);
    if (validateAbort) validateAbort.abort();
    const ctrl = new AbortController();
    validateAbort = ctrl;
    const key = dataKey(d.config);
    d.checking = true;
    d.reportError = null;
    notifyReport();
    try {
      const report = await api(`/api/projects/${enc(d.id)}/validate`, { method: 'POST', body: { config: d.config }, signal: ctrl.signal });
      d.report = report && typeof report === 'object' ? report : null;
      d.reportKey = key;
      d.reportAt = new Date().toISOString();
      if (manual && alive && d.report) {
        const samples = d.report.samples || [];
        const found = samples.filter((s) => s.bin_dir).length;
        const errs = (d.report.errors || []).length;
        const warns = (d.report.warnings || []).length;
        toast(`Data check: ${found} of ${plural(samples.length, 'section')} found${errs ? `, ${plural(errs, 'error')}` : ''}${warns ? `, ${plural(warns, 'warning')}` : ''}.`,
          d.report.ok ? 'success' : 'warning');
      }
    } catch (e) {
      if (e && e.name === 'AbortError') return;
      d.reportError = e;
      if (manual && alive) toast(`The data check failed: ${e.message}`, 'error');
    } finally {
      if (validateAbort === ctrl) {
        validateAbort = null;
        d.checking = false;
        notifyReport();
      }
    }
  }

  const reportStale = () => !!d.report && d.reportKey !== dataKey(d.config);

  /** [from, to] for every sample whose id changed since the last save. */
  function idRenames() {
    const out = [];
    for (const s of Array.isArray(d.config.samples) ? d.config.samples : []) {
      const orig = s && typeof s === 'object' ? d.originalIds.get(s) : undefined;
      if (orig !== undefined && s.id && orig !== s.id) out.push([orig, s.id]);
    }
    return out;
  }

  /** Renames whose saved region on the server has to move. */
  function sampleRenames() {
    return idRenames().filter(([from, to]) => d.rois.has(from) && !d.rois.has(to));
  }

  /** In-browser drawings (saved or not) follow their section's new id. */
  function moveRoiEntries(renames) {
    const entries = (d.roi && d.roi.entries) || {};
    for (const [from, to] of renames) {
      if (entries[from] && !(entries[to] && entries[to].loaded)) { entries[to] = entries[from]; delete entries[from]; }
      if (d.roi && d.roi.sample === from) d.roi.sample = to;
    }
  }

  // A drawn region is stored under its sample id; when the id is renamed, move the region with it. On a
  // case-insensitive file system a case-only rename (ctrl-1 -> CTRL-1) is one file: the PUT renames it (the server
  // drops the old spelling) and no DELETE is sent. The list of saved regions is then reloaded from the server.
  async function moveRenamedRois(renames) {
    const roiUrl = (sid) => `/api/projects/${enc(d.id)}/roi/${enc(sid)}`;
    for (const [from, to] of renames) {
      try {
        const gj = await api(roiUrl(from));
        await api(roiUrl(to), { method: 'PUT', body: gj });
        if (from.toLowerCase() !== to.toLowerCase()) await api(roiUrl(from), { method: 'DELETE' });
        toast(`The drawn region of ${from} now belongs to ${to}.`, 'info');
      } catch (e) {
        toast(`The drawn region of ${from} could not be moved to ${to}: ${e.message}`, 'warning');
      }
    }
    try {
      const fresh = await api(`/api/projects/${enc(d.id)}`);
      if (Array.isArray(fresh.rois)) setRois(d, fresh.rois);
    } catch (e) {
      toast(`Could not reload the list of drawn regions: ${e.message}`, 'warning');
    }
  }

  /** Save the draft (PUT). Resolves to true when the saved settings now match the draft; a call made while a save is
   *  in flight shares its result. */
  function save() {
    if (d.saving && d.savePromise) return d.savePromise;
    d.savePromise = writeProject();
    return d.savePromise;
  }

  async function writeProject() {
    computeDirty();
    if (!d.dirty) return true;
    const name = d.name.trim();
    if (!name) {
      toast('The project needs a name.', 'warning');
      if (alive) nameInput.focus();
      return false;
    }
    d.saving = true;
    refreshChrome();
    const sent = snapshot(d.name, d.config);
    const renames = sampleRenames();
    const entryRenames = idRenames();
    let replaced = false;
    try {
      const res = await api(`/api/projects/${enc(d.id)}`, { method: 'PUT', body: { name, config: d.config } });
      const serverName = res.name || name;
      const serverConfig = res.config || d.config;
      const editedMeanwhile = snapshot(d.name, d.config) !== sent;
      d.saved = { name: serverName, config: clone(serverConfig) };
      d.savedJson = snapshot(serverName, serverConfig);
      d.validation = normValidation(res.validation);
      d.updated = res.updated || d.updated;
      d.remoteUpdated = d.updated;
      if (Array.isArray(res.rois)) setRois(d, res.rois);
      if (!editedMeanwhile) {
        if (snapshot(d.name, d.config) !== d.savedJson) {
          d.name = serverName;
          d.config = clone(serverConfig);
          replaced = true;
        }
        d.originalIds = new WeakMap((Array.isArray(d.config.samples) ? d.config.samples : [])
          .filter((s) => s && typeof s === 'object').map((s) => [s, s.id]));
      }
      moveRoiEntries(entryRenames);
      if (renames.length) await moveRenamedRois(renames);
      const errs = d.validation.errors.length;
      if (errs) toast(`Saved with ${plural(errs, 'configuration error')}: runs cannot start until ${errs === 1 ? 'it is' : 'they are'} fixed.`, 'warning');
      else toast('Project saved', 'success', { timeout: 2000 });
      return true;
    } catch (e) {
      toast(`The project could not be saved: ${e.message}`, 'error');
      return false;
    } finally {
      d.saving = false;
      if (alive) {
        if (replaced) {
          if (document.activeElement !== nameInput) nameInput.value = d.name;
          renderBody();
        }
        refreshChrome();
      }
    }
  }

  async function discard() {
    computeDirty();
    if (!d.dirty || d.saving) return;
    const ok = await confirmDialog({
      title: 'Discard your changes?',
      message: 'The project returns to its saved settings. Drawn regions are not affected (they are saved separately).',
      confirmLabel: 'Discard changes', danger: true,
    });
    if (!ok) return;
    try {
      adopt(d, await api(`/api/projects/${enc(d.id)}`));
    } catch (e) {
      toast(`The saved version could not be loaded (${e.message}); restored the version loaded earlier.`, 'warning');
      const saved = d.saved;
      d.name = saved.name;
      d.config = clone(saved.config);
    }
    if (!alive) return;
    nameInput.value = d.name;
    renderBody();
    refreshChrome();
    toast('Changes discarded', 'info', { timeout: 2000 });
  }

  // ---------------------------------------------------------------- tabs
  function switchTab(tabId) {
    if (tabId === current) return;
    current = tabId;
    setPathSilently(`/projects/${enc(d.id)}/${tabId}`);
    renderBody();
    paintTabs();
  }

  function renderBody() {
    if (tabApi.destroy) {
      try { tabApi.destroy(); } catch (e) { console.error(e); }
    }
    tabApi = {};
    const renderers = { data: dataTab, parameters: parametersTab, signatures: signaturesTab, roi: roiTab, run: runTab };
    try {
      tabApi = renderers[current](body) || {};
      if (tabApi.start) tabApi.start();
    } catch (e) {
      console.error(e);
      mount(body, errorPanel(e, { retry: renderBody }));
    }
  }

  // ---- Data
  function dataTab(host) {
    const cfg = d.config;
    const rootInput = h('input', { class: 'input mono', id: `data-root-${d.id}`, value: cfg.data_root ?? '', spellcheck: false, autocomplete: 'off', placeholder: '.' });
    rootInput.addEventListener('input', () => { cfg.data_root = rootInput.value.trim(); changed({ data: true }); });
    const browse = button({
      label: 'Browse…', icon: 'folder',
      onClick: async () => {
        const picked = await openFolderPicker({ title: 'Data root', start: cfg.data_root || '.', note: 'The folder that contains the sample folders (Space Ranger outputs).' });
        if (!picked) return;
        cfg.data_root = picked;
        if (alive) rootInput.value = picked;
        changed({ data: true });
      },
    });
    const detect = button({
      label: 'Detect samples', icon: 'search',
      onClick: async () => {
        const res = await api(`/api/projects/${enc(d.id)}/detect-samples`, { method: 'POST', body: { data_root: cfg.data_root || '.' } });
        if (alive) chooseSamples(res);
      },
    });
    const binId = `bin-size-${d.id}`;
    const qcId = `qc-bin-size-${d.id}`;
    const binSel = h('select', { class: 'select', id: binId });
    const qcSel = h('select', { class: 'select', id: qcId });
    binSel.addEventListener('change', () => { cfg.bin_size = binSel.value; changed({ data: true }); });
    qcSel.addEventListener('change', () => { cfg.qc_bin_size = qcSel.value; changed({ data: true }); });
    const editor = sampleEditor({ cfg, report: d.report, onChange: () => changed({ data: true }) });
    const reportBox = h('div');
    const recheck = button({ label: 'Check again', icon: 'refresh', small: true, onClick: () => { validateNow(true); } });

    function fillBins() {
      const found = ((d.report && d.report.samples) || []).filter((s) => s.bin_dir && Array.isArray(s.bin_sizes));
      const union = uniq(found.flatMap((s) => s.bin_sizes)).sort();
      const options = union.length ? union : DEFAULT_BINS;
      for (const [sel, key] of [[binSel, 'bin_size'], [qcSel, 'qc_bin_size']]) {
        const cur = cfg[key];
        const items = uniq([...options, ...(cur ? [cur] : [])]).sort().map((b) => {
          let note = '';
          if (union.length && !union.includes(b)) note = ' — not found';
          else if (found.length > 1 && !found.every((s) => s.bin_sizes.includes(b))) note = ' — missing in some sections';
          return { b, text: `${binLabel(b)} (${b})${note}` };
        });
        const sig = `${items.map((i) => i.text).join('|')}#${cur}`;
        if (sel.dataset.sig === sig) continue;
        sel.dataset.sig = sig;
        mount(sel, items.map((i) => h('option', { value: i.b, selected: i.b === cur }, i.text)));
      }
    }

    function paintReport() {
      mount(reportBox, reportPanel({ report: d.report, cfg, checking: d.checking, stale: reportStale(), error: d.reportError, checkedAt: d.reportAt }));
      recheck.disabled = d.checking;
      fillBins();
      editor.setReport(d.report);
    }

    function suggestReplace() {
      const samples = Array.isArray(d.config.samples) ? d.config.samples : [];
      const pairs = (list) => JSON.stringify((list || []).map((s) => [s.id, s.path]));
      const reported = d.report && !reportStale() ? d.report.samples || [] : [];
      return !samples.length || pairs(samples) === pairs(defCfg.samples) || (reported.length > 0 && reported.every((s) => !s.bin_dir));
    }

    // The detection request may outlive this tab (tab switch, or a save that re-rendered it): apply the choice to the
    // current draft and repaint whichever sample table is on screen.
    async function chooseSamples(res) {
      const choice = await chooseDetectedSamples({ result: res, cfg: d.config, suggestReplace: suggestReplace() });
      if (!choice) return;
      const added = applyDetectedSamples(d.config, choice.candidates, { replace: choice.replace });
      if (alive && tabApi === tabHandle) editor.render();
      else if (alive && current === 'data') renderBody();
      changed({ data: true });
      toast(`${plural(added.length, 'sample')} added. Check their labels and conditions.`, 'success');
    }

    mount(host, h('div', { class: 'stack-lg' },
      card({
        title: 'Data location',
        body: h('div', { class: 'stack' },
          field({
            label: 'Data root', forId: rootInput.id, control: h('div', { class: 'input-group' }, rootInput, browse),
            help: 'The folder that holds the sample folders; relative sample folders are resolved against it. A relative data root is resolved against the folder the server was started in.',
          }),
          h('div', { class: 'row' }, detect, h('span', { class: 'muted small' }, 'Finds the Space Ranger outputs under the data root (up to four levels deep) and adds them as samples.')),
          h('div', { class: 'form-grid' },
            field({ label: 'Analysis bin size', forId: binId, control: binSel, help: 'The square_<size> folder that is analysed.' }),
            field({ label: 'Quality-control bin size', forId: qcId, control: qcSel, help: 'Per-bin QC of Supplementary Figure S1.' }))),
      }),
      card({ title: 'Samples', subtitle: 'One row per Visium HD section, in figure order.', body: editor.el, flush: true }),
      card({ title: 'Data check', subtitle: 'Sample folders, bin sizes, images and gene coverage, checked with the settings above (saved or not).', actions: recheck, body: reportBox })));
    paintReport();
    const tabHandle = { onReport: paintReport, start: () => { if (!d.checking) validateNow(false); } };
    return tabHandle;
  }

  // ---- Parameters
  function parametersTab(host) {
    const form = paramForm({
      cfg: d.config, defaults: defCfg, schema: defaults.param_schema || [], stepLabels: defaults.step_labels || {},
      onChange: () => changed(),
    });
    mount(host, form.el);
    return {};
  }

  // ---- Signatures
  function signaturesTab(host) {
    const hint = h('div');
    const paintHint = () => mount(hint, !d.report && !d.checking ? banner('info', h('span', null,
      'Gene chips are colored once the data has been checked (green: in the probe set, red: missing). ',
      h('button', { type: 'button', class: 'linklike', onClick: () => { validateNow(true); } }, 'Check the data now'))) : null);
    const editor = signatureEditor({ cfg: d.config, defaults: defCfg, report: d.report, ui: d.ui, onChange: () => changed() });
    mount(host, h('div', null, hint, editor.el));
    paintHint();
    return { onReport: () => { paintHint(); editor.setReport(d.report); } };
  }

  // ---- Region of interest
  function roiSampleStatus(sid) {
    const saved = d.saved.config || {};
    const savedSample = (Array.isArray(saved.samples) ? saved.samples : []).find((s) => s.id === sid);
    if (!savedSample) {
      return {
        ok: false, canSave: true, title: 'Save the project first',
        text: `Section ${sid} is not in the saved settings yet, and its H&E image is found through the saved data settings.`,
      };
    }
    const draftSample = (Array.isArray(d.config.samples) ? d.config.samples : []).find((s) => s.id === sid) || {};
    const differs = saved.data_root !== d.config.data_root || savedSample.path !== draftSample.path || saved.bin_size !== d.config.bin_size;
    return { ok: true, note: differs ? 'The data settings of this section have unsaved changes; the image shown comes from the saved settings.' : '' };
  }

  function roiTab(host) {
    const editor = roiEditor({
      projectId: d.id, cfg: d.config, store: d.roi, savedRois: d.rois,
      getVersion: () => d.updated || '',
      sampleStatus: roiSampleStatus,
      onConfigChange: () => changed(),
      onRoisChange: () => refreshChrome(),
      onSaveProject: () => save(),
    });
    mount(host, editor.el);
    return { destroy: editor.destroy, onChrome: editor.refresh };
  }

  // ---- Run
  function runTab(host) {
    let destroyed = false;
    const cfg = d.config;
    const stepLabels = defaults.step_labels || {};
    const ui = d.ui;
    const steps = () => (cfg.steps && typeof cfg.steps === 'object' ? cfg.steps : {});
    const stepOn = (k) => (has(ui.runSteps, k) ? ui.runSteps[k] : steps()[k] === true);
    const deterministic = () => (typeof ui.runDeterministic === 'boolean' ? ui.runDeterministic : cfg.deterministic === true);
    const statusBox = h('div', { class: 'stack' });
    const errorBox = h('div');
    const runsBox = h('div', null, loading('Loading runs…'));
    const labelId = `run-label-${d.id}`;
    const labelInput = h('input', { class: 'input', id: labelId, value: ui.runLabel || '', placeholder: 'e.g. first pass with the default settings', maxlength: 200 });
    labelInput.addEventListener('input', () => { ui.runLabel = labelInput.value; });
    const startBtn = h('button', { type: 'button', class: 'btn primary start-btn' });
    const paintStart = (busy) => mount(startBtn, busy ? h('span', { class: 'spinner' }) : icon('play', { size: 16 }), busy ? 'Starting…' : 'Start analysis');
    paintStart(false);
    let starting = false;
    startBtn.addEventListener('click', async () => {
      if (startBtn.disabled || starting) return;
      starting = true;
      startBtn.disabled = true;
      paintStart(true);
      try {
        await start();
      } catch (e) {
        showStartError(e.message || String(e));
      } finally {
        starting = false;
        if (!destroyed) { paintStart(false); paintStatus(); }
      }
    });

    function showStartError(message) {
      if (destroyed) { toast(message, 'error'); return; }
      mount(errorBox, message ? banner('danger', [h('strong', null, 'The analysis did not start. '), message]) : null);
    }

    const missingPackages = () => (state.health && Array.isArray(state.health.missing_required) ? state.health.missing_required : []);

    function paintStatus() {
      if (destroyed) return;
      const items = [];
      const missing = missingPackages();
      const errs = d.validation.errors;
      if (missing.length) items.push(banner('danger', [h('strong', null, 'Missing Python packages: '), missing.join(', '), '. The analysis cannot start until they are installed.']));
      if (errs.length) {
        items.push(banner(d.dirty ? 'warning' : 'danger', [
          h('strong', null, d.dirty ? `The saved settings have ${plural(errs.length, 'error')}; they are checked again when you start: `
            : `${plural(errs.length, 'configuration error')} must be fixed before the analysis can start:`),
          h('ul', null, errs.map((t) => h('li', null, t)))]));
      }
      if (d.dirty) items.push(banner('info', 'You have unsaved changes: starting saves them first.'));
      const rep = d.report;
      if (d.checking) items.push(banner('info', 'Checking the data…'));
      else if (!rep) {
        items.push(banner('info', h('span', null, 'The data has not been checked. ',
          h('button', { type: 'button', class: 'linklike', onClick: () => { validateNow(true); } }, 'Check it now'))));
      } else if ((rep.errors || []).length) {
        items.push(banner('warning', [h('strong', null, `The data check found ${plural(rep.errors.length, 'problem')}`, reportStale() ? ' (settings changed since)' : '', ':'),
          h('ul', null, rep.errors.map((t) => h('li', null, t)))]));
      } else if (reportStale()) {
        items.push(banner('info', h('span', null, 'The data settings changed since the last check. ',
          h('button', { type: 'button', class: 'linklike', onClick: () => { validateNow(true); } }, 'Check again'))));
      }
      const roi = cfg.roi || {};
      const pending = roiUnsaved(d);
      if (roi.use_drawn && pending.length) items.push(banner('warning', `Unsaved drawn region${pending.length === 1 ? '' : 's'} (${pending.join(', ')}): runs use the saved regions only. Save them in the Region of interest tab.`));
      if (roi.use_drawn && !d.rois.size) items.push(banner('info', 'Drawn regions are switched on, but no section has a saved region: the analysis is not restricted.'));
      if (!items.length) items.push(banner('success', 'Ready to start.'));
      mount(statusBox, items);
      if (!starting) startBtn.disabled = missing.length > 0 || (errs.length > 0 && !d.dirty);
    }

    function stepRow(k) {
      const label = stepLabels[k] || k;
      if (!SLOW_STEPS.includes(k)) {
        const on = steps()[k] === true;
        return h('li', { class: on ? 'plan-step' : 'plan-step off' }, icon(on ? 'check' : 'x', { size: 14 }), h('span', null, label));
      }
      const note = h('span', { class: 'muted small' });
      const paintNote = () => {
        note.textContent = has(ui.runSteps, k) ? ` (this run only; project setting: ${steps()[k] === true ? 'on' : 'off'})` : '';
      };
      paintNote();
      return h('li', { class: 'plan-step slow' }, switchToggle({
        checked: stepOn(k), label,
        onChange: (v) => {
          if (v === (steps()[k] === true)) delete ui.runSteps[k]; else ui.runSteps[k] = v;
          paintNote();
        },
      }), note);
    }

    function planCard() {
      const samples = Array.isArray(cfg.samples) ? cfg.samples : [];
      const keys = uniq([...Object.keys(stepLabels), ...Object.keys(steps())]);
      const detNote = h('span', { class: 'muted small' });
      const paintDet = () => {
        detNote.textContent = typeof ui.runDeterministic === 'boolean' ? ` (this run only; project setting: ${cfg.deterministic === true ? 'on' : 'off'})` : '';
      };
      paintDet();
      const roi = cfg.roi || {};
      return card({
        title: 'What will run',
        subtitle: 'The saved settings, with the switches below applied to this run only.',
        body: h('div', { class: 'stack' },
          h('dl', { class: 'kv' },
            h('dt', null, 'Sections'), h('dd', null, samples.length ? h('ul', { class: 'plan-samples' }, samples.map((s) => h('li', null,
              h('span', { class: 'swatch', style: { background: /^#[0-9A-Fa-f]{6}$/.test(String(s.color)) ? s.color : 'transparent' } }),
              h('span', null, s.label || s.id), h('span', { class: 'mono small muted' }, s.id), h('span', { class: 'muted small' }, s.condition || ''),
              s.id === cfg.reference_sample ? h('span', { class: 'badge' }, 'reference') : null))) : h('span', { class: 'muted' }, 'none')),
            h('dt', null, 'Bin size'), h('dd', null, `${binLabel(cfg.bin_size)} (QC figure: ${binLabel(cfg.qc_bin_size)})`),
            h('dt', null, 'Drawn regions'), h('dd', null, roi.use_drawn
              ? (samples.length ? samples.map((s) => `${s.id}: ${d.rois.has(s.id) ? 'saved region' : 'none (not restricted)'}`).join(' · ') : '—')
              : 'not used')),
          h('h4', null, 'Steps'),
          h('ul', { class: 'plan-steps' }, keys.map(stepRow)),
          h('div', { class: 'plan-det' }, switchToggle({
            checked: deterministic(), label: 'Deterministic (single-threaded, reproducible; slower)',
            onChange: (v) => {
              ui.runDeterministic = v === (cfg.deterministic === true) ? null : v;
              paintDet();
            },
          }), detNote)),
      });
    }

    async function start() {
      showStartError(null);
      if (missingPackages().length) { showStartError(`Missing Python packages: ${missingPackages().join(', ')}.`); return; }
      if (computeDirty()) {
        const ok = await confirmDialog({
          title: 'Save changes and start?',
          message: 'The analysis runs on the saved settings. Your unsaved changes are saved first.',
          confirmLabel: 'Save and start',
        });
        if (!ok) return;
        if (!(await save())) return;
      }
      if (d.validation.errors.length) {
        showStartError(`Fix the ${plural(d.validation.errors.length, 'configuration error')} listed above first.`);
        return;
      }
      if (d.report && !reportStale() && (d.report.errors || []).length) {
        const ok = await confirmDialog({
          title: 'Start despite data problems?',
          message: `The data check reported ${plural(d.report.errors.length, 'problem')}; the analysis will probably fail.`,
          confirmLabel: 'Start anyway', danger: true,
        });
        if (!ok) return;
      }
      const pending = roiUnsaved(d);
      if (pending.length && cfg.roi && cfg.roi.use_drawn) {
        const ok = await confirmDialog({
          title: 'Start without the unsaved regions?',
          message: `The drawn region${pending.length === 1 ? '' : 's'} of ${pending.join(', ')} ${pending.length === 1 ? 'is' : 'are'} not saved; the run uses the saved regions only.`,
          confirmLabel: 'Start anyway',
        });
        if (!ok) return;
      }
      const runSteps = Object.fromEntries(SLOW_STEPS.map((k) => [k, stepOn(k)]));
      try {
        const run = await api(`/api/projects/${enc(d.id)}/runs`, {
          method: 'POST', body: { label: (ui.runLabel || '').trim(), steps: runSteps, deterministic: deterministic() },
        });
        ui.runLabel = '';
        ui.runSteps = {};
        ui.runDeterministic = null;
        toast('Analysis queued', 'success', { timeout: 2500 });
        navigate(`/runs/${enc(run.id)}`);
      } catch (e) {
        if (e.status === 409) { showStartError(e.message); return; }
        throw e;
      }
    }

    function paintRuns(runs) {
      if (!runs.length) {
        mount(runsBox, h('p', { class: 'muted' }, 'No runs yet.'));
        return;
      }
      mount(runsBox, h('div', { class: 'table-wrap' }, h('table', { class: 'table compact' },
        h('thead', null, h('tr', null, h('th', null, 'Run'), h('th', null, 'Status'), h('th', null, 'Created'), h('th', { class: 'num' }, 'Duration'), h('th', null, 'Progress'))),
        h('tbody', null, runs.map((r) => {
          const href = `#/runs/${enc(r.id)}`;
          const prog = r.progress || {};
          const active = r.status === 'running' || r.status === 'queued';
          return h('tr', null,
            h('td', null, h('a', { href }, r.label || r.id), r.label ? h('div', { class: 'mono small muted' }, r.id) : null),
            h('td', null, h('a', { href, class: 'plain-link' }, statusBadge(r.status))),
            h('td', { title: fmt.date(r.created) }, fmt.ago(r.created)),
            h('td', { class: 'num' }, fmt.duration(r.duration_s)),
            h('td', null, active && r.status === 'running' && Number.isFinite(prog.fraction)
              ? h('div', { class: 'run-progress', title: prog.label || '' }, h('div', { class: 'progress' }, h('div', { class: 'bar', style: { width: `${Math.round(prog.fraction * 100)}%` } })),
                h('span', { class: 'muted small' }, `${Math.round(prog.fraction * 100)}%`))
              : h('span', { class: 'muted small' }, r.status === 'failed' && r.error ? String(r.error).split('\n')[0].slice(0, 90) : '')));
        })))));
    }

    let poller = null;
    function startPolling() {
      if (poller) poller.stop();
      poller = poll(async () => {
        try {
          const runs = await api('/api/runs', { query: { project: d.id } });
          if (destroyed) return false;
          const list = Array.isArray(runs) ? runs : [];
          paintRuns(list);
          return list.some((r) => r.status === 'queued' || r.status === 'running');
        } catch (e) {
          if (!destroyed) mount(runsBox, errorPanel(e, { retry: () => { mount(runsBox, loading('Loading runs…')); startPolling(); } }));
          return false;
        }
      }, RUNS_POLL_MS);
    }

    const onHealth = () => paintStatus();
    window.addEventListener('health', onHealth);

    mount(host, h('div', { class: 'stack-lg' },
      h('div', { class: 'grid-2 run-grid' },
        planCard(),
        card({
          title: 'Start',
          body: h('div', { class: 'stack' }, statusBox,
            field({ label: 'Run label (optional)', forId: labelId, control: labelInput, help: 'Shown in the list of runs.' }),
            h('div', { class: 'row' }, startBtn), errorBox),
        })),
      card({ title: 'Runs of this project', actions: h('a', { class: 'btn small ghost', href: '#/runs' }, 'All runs'), body: runsBox })));
    paintStatus();
    startPolling();
    return {
      onChrome: paintStatus,
      onReport: paintStatus,
      destroy() {
        destroyed = true;
        if (poller) poller.stop();
        window.removeEventListener('health', onHealth);
      },
    };
  }

  // ---------------------------------------------------------------- keyboard, first paint, cleanup
  const onKey = (e) => {
    if (!(e.ctrlKey || e.metaKey) || e.altKey || (e.key !== 's' && e.key !== 'S')) return;
    if (document.querySelector('.modal-backdrop')) return;
    e.preventDefault();
    save();
  };
  document.addEventListener('keydown', onKey);

  renderBody();
  refreshChrome();

  return () => {
    alive = false;
    clearTimeout(validateTimer);
    if (validateAbort) {
      validateAbort.abort();
      validateAbort = null;
      d.checking = false;
    }
    if (tabApi.destroy) {
      try { tabApi.destroy(); } catch (e) { console.error(e); }
    }
    tabApi = {};
    document.removeEventListener('keydown', onKey);
  };
}
