// sample-editor.js — the Data tab's sample table, the dialog that turns detected Space Ranger folders into samples,
// and the panel that shows the data-check report (POST /api/projects/{pid}/validate).
import { h, mount, icon, button, modal, banner, fmt, api } from '../core.js';
import { iconButton } from './chips.js';
import { openFolderPicker } from './folder-picker.js';

export const SAMPLE_COLORS = ['#4C72B0', '#C44E52', '#55A868', '#8172B2', '#CCB974', '#64B5CD', '#DD8452', '#937860',
  '#DA8BC3', '#8C8C8C'];
export const CONDITIONS = ['control', 'tumor', 'other'];
const ID_RE = /^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$/;
const HEX_RE = /^#[0-9A-Fa-f]{6}$/;
const CONTROL_WORDS = new Set(['ctrl', 'ctl', 'control', 'contra', 'contralateral', 'sham', 'naive', 'healthy', 'normal',
  'vehicle', 'veh', 'wt', 'wildtype']);
const TUMOR_WORDS = new Set(['tum', 'tumor', 'tumour', 'tumorbearing', 'lesion', 'met', 'mets', 'metastasis', 'cancer',
  'carcinoma', 'inj', 'injected', 'a549']);
let editorSeq = 0;

// ------------------------------------------------------------------------------------------------ pure helpers
export function isValidSampleId(id) { return ID_RE.test(String(id ?? '')); }

/** A valid sample id from free text: other characters become '-', at most 64 characters. */
export function sanitizeId(text) {
  const s = String(text ?? '').replace(/[^A-Za-z0-9_.-]+/g, '-').replace(/^[^A-Za-z0-9]+/, '').slice(0, 64);
  return s || 'sample';
}

/** `base`, or `base-2`, `base-3`, ... (`sep` between) — the first value not in `taken`. */
export function uniqueValue(base, taken, sep = '-') {
  if (!taken.has(base)) return base;
  for (let i = 2; ; i += 1) {
    const suffix = `${sep}${i}`;
    const v = `${base.slice(0, 64 - suffix.length)}${suffix}`;
    if (!taken.has(v)) return v;
  }
}

/** First palette color not used yet; cycles through the palette once every color is taken. */
export function nextColor(used) {
  const taken = new Set([...(used || [])].filter(Boolean).map((c) => String(c).toLowerCase()));
  return SAMPLE_COLORS.find((c) => !taken.has(c.toLowerCase())) || SAMPLE_COLORS[taken.size % SAMPLE_COLORS.length];
}

/** 'control', 'tumor' or null from the words of a folder name or id ("CTRL-1", "tumour_B", "SR01_sham2"). */
export function guessCondition(text) {
  const words = String(text ?? '').toLowerCase().split(/[^a-z0-9]+/).map((w) => w.replace(/\d+$/, '')).filter(Boolean);
  if (words.some((w) => CONTROL_WORDS.has(w))) return 'control';
  if (words.some((w) => TUMOR_WORDS.has(w))) return 'tumor';
  return null;
}

/** '016um' -> '16 µm' */
export function binLabel(bin) {
  const m = /^(\d{3})um$/.exec(String(bin));
  return m ? `${Number(m[1])} µm` : String(bin);
}

export function isAbsolutePath(p) { return /^([A-Za-z]:[\\/]|[\\/]|~)/.test(String(p ?? '')); }

export function joinPath(a, b) {
  const x = String(a ?? '').trim();
  const y = String(b ?? '').trim();
  if (!y || y === '.') return x || '.';
  if (!x || x === '.') return y;
  return `${x.replace(/[\\/]+$/, '')}/${y.replace(/^[\\/]+/, '')}`;
}

/** `target` relative to `root` ('/' separators) when it lies inside it, '.' when equal, null otherwise. */
export function relativePath(root, target) {
  const norm = (p) => String(p ?? '').replace(/\\/g, '/').replace(/\/+$/, '');
  const r = norm(root);
  const t = norm(target);
  if (!t) return null;
  const caseless = /^[A-Za-z]:/.test(r) || r.startsWith('//');
  const k = (x) => (caseless ? x.toLowerCase() : x);
  if (k(t) === k(r)) return '.';
  if (k(t).startsWith(`${k(r)}/`)) return t.slice(r.length + 1);
  return null;
}

const normPath = (p) => String(p ?? '').replace(/\\/g, '/').replace(/\/+$/, '').toLowerCase();

/**
 * Add detected candidates ({folder, path, suggested_id}) to cfg.samples (or replace them): unique ids, conditions
 * guessed from the folder names (the first sample without a hint is the control), distinct colors. Keeps
 * cfg.reference_sample valid. Returns the samples added.
 */
export function applyDetectedSamples(cfg, candidates, { replace = false } = {}) {
  const base = replace ? [] : (Array.isArray(cfg.samples) ? cfg.samples.slice() : []);
  const ids = new Set(base.map((s) => s.id));
  const labels = new Set(base.map((s) => s.label).filter(Boolean));
  const colors = base.map((s) => s.color);
  let hasControl = base.some((s) => s.condition === 'control');
  const added = [];
  for (const c of candidates || []) {
    const id = uniqueValue(sanitizeId(c.suggested_id || c.folder), ids);
    ids.add(id);
    const condition = guessCondition(`${c.folder || ''} ${c.suggested_id || ''}`) || (hasControl ? 'tumor' : 'control');
    if (condition === 'control') hasControl = true;
    const label = uniqueValue(condition === 'control' ? 'control' : 'tumor-bearing', labels, ' ');
    labels.add(label);
    const color = nextColor(colors);
    colors.push(color);
    added.push({ id, label, condition, color, path: c.path || c.folder || '', microns_per_pixel: null });
  }
  cfg.samples = [...base, ...added];
  if (!cfg.samples.some((s) => s.id === cfg.reference_sample)) {
    const ref = cfg.samples.find((s) => s.condition === 'control') || cfg.samples[0];
    cfg.reference_sample = ref ? ref.id : '';
  }
  return added;
}

function moveMaskKey(cfg, from, to) {
  const masks = cfg.roi && cfg.roi.masks;
  if (!masks || typeof masks !== 'object' || from === to || !(from in masks) || to in masks) return;
  masks[to] = masks[from];
  delete masks[from];
}

// ------------------------------------------------------------------------------------------------ sample table
/**
 * Editable table of cfg.samples. Edits change `cfg` in place and call onChange(); add, remove and reorder rebuild
 * the rows. setReport(report) colors the rows with a data-check report; render() rebuilds after outside changes.
 * onExplore(sampleId), when given, adds an Explore button to every section the report found.
 */
export function sampleEditor({ cfg, report = null, onChange, onExplore } = {}) {
  const uid = `smp-${++editorSeq}`;
  let currentReport = report;
  let refs = [];
  const emit = () => { if (onChange) onChange(); };
  const tbody = h('tbody');
  const addBtn = button({ label: 'Add sample', icon: 'plus', small: true, onClick: () => addSample() });
  const el = h('div', { class: 'sample-editor' },
    h('div', { class: 'table-wrap borderless' }, h('table', { class: 'table compact sample-table' },
      h('thead', null, h('tr', null,
        h('th', { title: 'Reference sample: fold changes and ratios are every other sample over this one' }, 'Ref.'),
        h('th', null, 'Id'), h('th', null, 'Label'), h('th', null, 'Condition'), h('th', null, 'Color'),
        h('th', null, 'Folder (relative to the data root)'),
        h('th', { title: 'Microns per full-resolution pixel. Empty: read from spatial/scalefactors_json.json' }, 'µm / px'),
        h('th', null, 'Data'), h('th', null, h('span', { class: 'sr-only' }, 'Order and removal')))),
      tbody)),
    h('div', { class: 'sample-foot' }, addBtn,
      h('span', { class: 'muted small' }, 'The reference sample (Ref.) is the denominator of every fold change and ratio of ratios.')));

  function checkIds() {
    // ids that differ only in letter case share one region file on Windows and macOS
    const key = (id) => String(id || '').toLowerCase();
    const counts = new Map();
    refs.forEach((r) => counts.set(key(r.s.id), (counts.get(key(r.s.id)) || 0) + 1));
    for (const r of refs) {
      const msg = !isValidSampleId(r.s.id)
        ? 'Use 1–64 letters, digits, “.”, “_” or “-”, starting with a letter or a digit'
        : counts.get(key(r.s.id)) > 1 ? 'Another sample has this id (letter case does not count)' : '';
      r.idInput.classList.toggle('invalid', !!msg);
      r.idInput.title = msg;
      r.idInput.setAttribute('aria-invalid', msg ? 'true' : 'false');
    }
  }

  function paintStatus(r) {
    const rep = currentReport && Array.isArray(currentReport.samples)
      ? currentReport.samples.find((x) => x.id === r.s.id) : null;
    r.mpp.placeholder = rep && rep.microns_per_pixel ? `auto (${fmt.num(rep.microns_per_pixel, 4)})` : 'auto';
    r.tr.classList.toggle('row-bad', !!rep && !rep.bin_dir);
    if (!rep) {
      mount(r.status, h('span', { class: 'faint', title: 'Not checked yet' }, '—'));
      return;
    }
    const found = !!rep.bin_dir;
    let explore = null;
    if (found && onExplore) {
      explore = button({ label: 'Explore', icon: 'eye', small: true, kind: 'ghost', onClick: () => onExplore(r.s.id) });
      explore.title = 'Look at the data of this section: summary, maps of UMIs, quality control and any gene, gene list';
      explore.setAttribute('aria-label', `Explore the data of section ${r.s.id}`);
    }
    mount(r.status, h('div', { class: 'sample-status-row' },
      h('span', { class: `state ${found ? 'ok' : 'bad'}`, title: rep.message || (found ? rep.bin_dir : 'Not found') },
        icon(found ? 'check' : 'alert', { size: 13 }),
        found ? (rep.n_bins !== null && rep.n_bins !== undefined ? `${fmt.int(rep.n_bins)} bins` : 'Found') : 'Not found'),
      explore));
  }

  async function browse(s, pathInput) {
    const root = String(cfg.data_root ?? '').trim() || '.';
    const cur = String(s.path ?? '').trim();
    const picked = await openFolderPicker({
      title: `Folder of sample ${s.id || ''}`.trim(),
      start: isAbsolutePath(cur) ? cur : joinPath(root, cur),
      note: 'Choose the sample folder, its outs/ or outs/binned_outputs/ folder, or one square_XXXum/ folder.',
    });
    if (!picked) return;
    let value = picked;
    try {
      const rootAbs = (await api('/api/fs/list', { query: { path: root } })).path;
      value = relativePath(rootAbs, picked) ?? picked;
    } catch {
      // the data root cannot be listed: keep the absolute path
    }
    s.path = value;
    pathInput.value = value;
    emit();
  }

  function buildRow(s, i, n) {
    const name = () => s.id || `sample ${i + 1}`;
    const radio = h('input', { type: 'radio', name: `${uid}-ref`, checked: cfg.reference_sample === s.id, 'aria-label': `Reference sample: ${name()}` });
    const idInput = h('input', { class: 'input mono', value: s.id ?? '', spellcheck: false, autocomplete: 'off', 'aria-label': `Sample ${i + 1}: id` });
    const label = h('input', { class: 'input', value: s.label ?? '', 'aria-label': `Sample ${i + 1}: label`, placeholder: 'shown in figures' });
    const conditions = CONDITIONS.includes(s.condition) || !s.condition ? CONDITIONS : [...CONDITIONS, s.condition];
    const cond = h('select', { class: 'select', 'aria-label': `Sample ${i + 1}: condition` },
      conditions.map((c) => h('option', { value: c, selected: c === (s.condition || 'other') }, c)));
    const color = h('input', { type: 'color', class: 'color-input', value: HEX_RE.test(String(s.color)) ? s.color : '#888888', 'aria-label': `Sample ${i + 1}: color` });
    const pathInput = h('input', { class: 'input mono', value: s.path ?? '', spellcheck: false, autocomplete: 'off', 'aria-label': `Sample ${i + 1}: folder`, placeholder: s.id ? `${s.id} (default)` : '' });
    const browseBtn = button({ icon: 'folder', small: true, title: 'Browse…', onClick: () => browse(s, pathInput) });
    browseBtn.setAttribute('aria-label', `Browse for the folder of sample ${i + 1}`);
    const mpp = h('input', { class: 'input num', type: 'number', min: 0, step: 'any', value: s.microns_per_pixel ?? '', placeholder: 'auto', 'aria-label': `Sample ${i + 1}: microns per pixel` });
    const status = h('td', { class: 'sample-status' });
    const up = iconButton('chevron-up', `Move sample ${i + 1} up`, () => move(i, -1), { disabled: i === 0 });
    const down = iconButton('chevron-down', `Move sample ${i + 1} down`, () => move(i, 1), { disabled: i === n - 1 });
    const remove = iconButton('trash', `Remove sample ${i + 1}`, () => removeAt(i), { cls: 'danger-hover' });

    radio.addEventListener('change', () => { if (radio.checked) { cfg.reference_sample = s.id; emit(); } });
    idInput.addEventListener('input', () => {
      const old = s.id;
      s.id = idInput.value.trim();
      if (radio.checked) cfg.reference_sample = s.id;
      moveMaskKey(cfg, old, s.id);
      radio.setAttribute('aria-label', `Reference sample: ${name()}`);
      pathInput.placeholder = s.id ? `${s.id} (default)` : '';
      checkIds();
      paintStatus(ref);
      emit();
    });
    label.addEventListener('input', () => { s.label = label.value; emit(); });
    cond.addEventListener('change', () => { s.condition = cond.value; emit(); });
    color.addEventListener('input', () => { s.color = color.value; emit(); });
    pathInput.addEventListener('input', () => { s.path = pathInput.value.trim(); emit(); });
    mpp.addEventListener('input', () => {
      const raw = mpp.value.trim();
      const v = Number(raw);
      const bad = raw !== '' && !(Number.isFinite(v) && v > 0);
      mpp.classList.toggle('invalid', bad || mpp.validity.badInput);
      mpp.title = bad || mpp.validity.badInput ? 'A positive number, or empty to read it from the data' : '';
      if (raw === '' && !mpp.validity.badInput) { s.microns_per_pixel = null; emit(); } else if (!bad && raw !== '') { s.microns_per_pixel = v; emit(); }
    });
    mpp.addEventListener('change', () => {
      if (mpp.classList.contains('invalid')) {
        mpp.value = s.microns_per_pixel ?? '';
        mpp.classList.remove('invalid');
        mpp.title = '';
      }
    });

    const tr = h('tr', null,
      h('td', { class: 'center' }, radio),
      h('td', { class: 'col-id' }, idInput),
      h('td', { class: 'col-label' }, label),
      h('td', { class: 'col-cond' }, cond),
      h('td', null, color),
      h('td', { class: 'col-path' }, h('div', { class: 'input-group' }, pathInput, browseBtn)),
      h('td', { class: 'col-mpp' }, mpp),
      status,
      h('td', { class: 'col-actions' }, h('div', { class: 'row nowrap' }, up, down, remove)));
    const ref = { s, tr, radio, idInput, mpp, status, up, down, remove };
    return ref;
  }

  function render(focus) {
    if (!Array.isArray(cfg.samples)) cfg.samples = [];
    const n = cfg.samples.length;
    refs = cfg.samples.map((s, i) => buildRow(s, i, n));
    mount(tbody, refs.length ? refs.map((r) => r.tr)
      : h('tr', null, h('td', { colSpan: 9, class: 'empty-row' }, 'No samples yet. Detect them in the data root, or add one by hand.')));
    checkIds();
    refs.forEach(paintStatus);
    if (!focus) return;
    const r = refs[focus.index];
    if (!r) { addBtn.focus(); return; }
    const target = r[focus.control];
    if (target && !target.disabled) target.focus();
    else if (focus.control === 'up' || focus.control === 'down') (focus.control === 'up' ? r.down : r.up).focus();
  }

  function move(i, delta) {
    const list = cfg.samples;
    const j = i + delta;
    if (j < 0 || j >= list.length) return;
    [list[i], list[j]] = [list[j], list[i]];
    render({ index: j, control: delta < 0 ? 'up' : 'down' });
    emit();
  }

  function removeAt(i) {
    const [gone] = cfg.samples.splice(i, 1);
    if (gone && cfg.reference_sample === gone.id && !cfg.samples.some((s) => s.id === gone.id)) {
      const next = cfg.samples.find((s) => s.condition === 'control') || cfg.samples[0];
      cfg.reference_sample = next ? next.id : '';
    }
    render({ index: Math.min(i, cfg.samples.length - 1), control: 'remove' });
    emit();
  }

  function addSample() {
    const list = cfg.samples;
    const id = uniqueValue(`sample-${list.length + 1}`, new Set(list.map((s) => s.id)));
    list.push({
      id, label: id, condition: list.some((s) => s.condition === 'control') ? 'tumor' : 'control',
      color: nextColor(list.map((s) => s.color)), path: '', microns_per_pixel: null,
    });
    if (!list.some((s) => s.id === cfg.reference_sample)) cfg.reference_sample = id;
    render({ index: list.length - 1, control: 'idInput' });
    emit();
  }

  render();
  return {
    el,
    render: () => render(),
    setReport(r) { currentReport = r; refs.forEach(paintStatus); },
  };
}

// ------------------------------------------------------------------------------------------------ detection dialog
/**
 * Dialog over a POST detect-samples result. Resolves to {candidates, replace} (the folders to add, and whether they
 * replace the current samples) or null.
 */
export function chooseDetectedSamples({ result, cfg, suggestReplace = false }) {
  return new Promise((resolve) => {
    let answered = false;
    const finish = (value) => { if (!answered) { answered = true; resolve(value); } };
    const candidates = (result && Array.isArray(result.candidates)) ? result.candidates : [];
    const root = (result && result.data_root) || cfg.data_root || '.';
    if (!candidates.length) {
      modal({
        title: 'No samples found',
        body: h('div', null,
          h('p', null, 'No Space Ranger output was found under ', h('code', null, root), '.'),
          h('p', { class: 'muted' }, 'The search looks in this folder and up to four levels of subfolders for '
            + 'square_XXXum/filtered_feature_bc_matrix.h5. Check the data root, or add the samples by hand.')),
        actions: [{ label: 'Close', kind: 'primary' }],
        onClose: () => finish(null),
      });
      return;
    }
    const current = Array.isArray(cfg.samples) ? cfg.samples : [];
    const paths = new Set(current.map((s) => normPath(s.path)));
    const ids = new Set(current.map((s) => s.id));
    const rows = candidates.map((c, i) => {
      const known = paths.has(normPath(c.path)) || ids.has(sanitizeId(c.suggested_id || c.folder));
      const id = `detect-${editorSeq}-${i}`;
      const box = h('input', { type: 'checkbox', id, checked: !known });
      const tr = h('tr', { class: known ? 'known' : '' },
        h('td', { class: 'center' }, box),
        h('td', null, h('label', { for: id, class: 'detect-folder' }, c.folder || c.path),
          known ? h('span', { class: 'badge', title: 'A sample with this folder or id is configured already' }, 'already added') : null),
        h('td', { class: 'mono small' }, c.path),
        h('td', { class: 'mono' }, sanitizeId(c.suggested_id || c.folder)),
        h('td', null, h('div', { class: 'pill-list' }, (c.bin_sizes || []).map((b) => h('span', { class: 'chip no-x' }, b)))));
      return { c, box, tr };
    });
    const all = h('input', { type: 'checkbox', 'aria-label': 'Select every folder' });
    const syncAll = () => {
      const n = rows.filter((r) => r.box.checked).length;
      all.checked = n === rows.length;
      all.indeterminate = n > 0 && n < rows.length;
    };
    all.addEventListener('change', () => { rows.forEach((r) => { r.box.checked = all.checked; }); });
    rows.forEach((r) => r.box.addEventListener('change', syncAll));
    syncAll();
    const replaceBox = h('input', { type: 'checkbox', checked: suggestReplace && current.length > 0 });
    const err = h('div', { class: 'error-text', role: 'alert' });
    modal({
      title: 'Detected samples', size: 'wide',
      onClose: () => finish(null),
      body: h('div', { class: 'stack detect' },
        h('p', null, `${candidates.length} Space Ranger output${candidates.length === 1 ? '' : 's'} found under `, h('code', null, root), '.'),
        h('div', { class: 'table-wrap' }, h('table', { class: 'table compact' },
          h('thead', null, h('tr', null, h('th', { class: 'center' }, all), h('th', null, 'Folder'), h('th', null, 'Path'),
            h('th', null, 'Sample id'), h('th', null, 'Bin sizes'))),
          h('tbody', null, rows.map((r) => r.tr)))),
        current.length ? h('label', { class: 'check' }, replaceBox,
          `Replace the current sample${current.length === 1 ? '' : 's'} (${current.map((s) => s.id).join(', ')})`) : null,
        h('p', { class: 'muted small' }, 'Labels, conditions and colors are guessed from the folder names; check them in the sample table afterwards.'),
        err),
      actions: [
        { label: 'Cancel', onClick: () => finish(null) },
        {
          label: 'Add selected', kind: 'primary', icon: 'plus',
          onClick: () => {
            const chosen = rows.filter((r) => r.box.checked).map((r) => r.c);
            if (!chosen.length) { err.textContent = 'Select at least one folder.'; return false; }
            finish({ candidates: chosen, replace: replaceBox.checked });
            return true;
          },
        },
      ],
    });
  });
}

// ------------------------------------------------------------------------------------------------ report panel
const yesNo = (ok, yes = 'yes', no = 'no') => h('span', { class: `state ${ok ? 'ok' : 'bad'}` }, icon(ok ? 'check' : 'x', { size: 12 }), ok ? yes : no);

function listBanner(kind, title, items) {
  return banner(kind, [h('strong', null, title), h('ul', null, items.map((t) => h('li', null, t)))]);
}

/**
 * The data-check report: errors and warnings, one row per sample (found, folder, bin sizes, bins, microns per pixel,
 * H&E image) and gene coverage (signatures with genes missing from the probe set, receptor panel, adipocyte genes).
 * onExplore(), when given, adds a pointer to the Explore tab once a section was found.
 */
export function reportPanel({ report, cfg, checking = false, stale = false, error = null, checkedAt = null, onExplore } = {}) {
  const status = checking
    ? h('div', { class: 'report-status' }, h('span', { class: 'spinner' }), 'Checking the data…')
    : h('div', { class: 'report-status' },
      report ? icon(report.ok ? 'check' : 'alert', { size: 16, cls: report.ok ? 'tone-ok' : 'tone-bad' }) : icon('info', { size: 16 }),
      h('span', null, report ? (report.ok ? 'Data found and readable.' : 'The data check found problems.') : 'The data has not been checked yet.'),
      report && checkedAt ? h('span', { class: 'muted' }, `Checked ${fmt.ago(checkedAt)}${stale ? ', settings changed since' : ''}.`) : null);
  const parts = [status];
  if (error) parts.push(banner('danger', `The data check failed: ${error.message || error}`));
  if (!report) return h('div', { class: 'report' }, parts);

  if (onExplore && (report.samples || []).some((s) => s.bin_dir)) {
    parts.push(h('div', { class: 'report-explore' }, icon('eye', { size: 16 }),
      h('span', null, h('strong', null, 'Read the data before any run. '),
        'The Explore tab shows each section: bins, genes, UMIs and quality control, maps of them and of any gene on the H&E image, and the genes ranked by UMIs.'),
      button({ label: 'Explore the data', icon: 'eye', small: true, onClick: () => onExplore() })));
  }

  const errors = report.errors || [];
  const warnings = report.warnings || [];
  if (errors.length) parts.push(listBanner('danger', `${errors.length} error${errors.length === 1 ? '' : 's'}`, errors));
  if (warnings.length) parts.push(listBanner('warning', `${warnings.length} warning${warnings.length === 1 ? '' : 's'}`, warnings));

  const samples = report.samples || [];
  parts.push(h('h4', { class: 'report-h' }, 'Sections'));
  parts.push(samples.length ? h('div', { class: 'table-wrap' }, h('table', { class: 'table compact report-table' },
    h('thead', null, h('tr', null, ['Sample', 'Found', 'Bin folder', 'Bin sizes', 'Bins', 'µm / px', 'Positions', 'Scale factors', 'H&E image']
      .map((t) => h('th', { class: t === 'Bins' || t === 'µm / px' ? 'num' : '' }, t)))),
    h('tbody', null, samples.map((s) => h('tr', { class: s.bin_dir ? '' : 'row-bad' },
      h('td', { class: 'mono' }, s.id),
      h('td', { title: s.message || '' }, yesNo(!!s.bin_dir, 'found', 'not found'), s.message && !s.bin_dir ? h('div', { class: 'small muted' }, s.message) : null),
      h('td', { class: 'mono small path-cell', title: s.bin_dir || '' }, s.bin_dir || '—'),
      h('td', null, h('div', { class: 'pill-list' }, (s.bin_sizes || []).map((b) => h('span', {
        class: `chip no-x ${cfg && b === cfg.bin_size ? 'ok' : ''}`.trim(), title: cfg && b === cfg.bin_size ? 'analysis bin size' : '',
      }, b)))),
      h('td', { class: 'num' }, fmt.int(s.n_bins)),
      h('td', { class: 'num' }, fmt.num(s.microns_per_pixel, 4)),
      h('td', null, yesNo(!!s.has_positions)),
      h('td', null, yesNo(!!s.has_scalefactors)),
      h('td', null, yesNo(!!s.has_hires_image))))))) : h('p', { class: 'muted' }, 'No sections in the report.'));

  const genes = report.genes || {};
  parts.push(h('h4', { class: 'report-h' }, 'Gene coverage'));
  if (genes.n_features === null || genes.n_features === undefined) {
    parts.push(h('p', { class: 'muted' }, 'Gene names are not available: no readable filtered_feature_bc_matrix.h5 was found.'));
  } else {
    parts.push(h('p', null, `${fmt.int(genes.n_features)} genes in the probe set.`));
    const minCore = cfg && cfg.scoring ? Number(cfg.scoring.min_core) : NaN;
    const sigRows = Object.entries(genes.signatures || {}).filter(([name, v]) => v && ((v.missing || []).length
      || (Number.isFinite(minCore) && cfg.signatures && cfg.signatures[name] && (cfg.signatures[name].core || []).length && v.core_present < minCore)));
    if (sigRows.length) {
      parts.push(h('div', { class: 'table-wrap' }, h('table', { class: 'table compact' },
        h('thead', null, h('tr', null, h('th', null, 'Signature'), h('th', null, 'Missing from the probe set'), h('th', { class: 'num' }, 'Core genes detected'))),
        h('tbody', null, sigRows.map(([name, v]) => {
          const sig = (cfg && cfg.signatures && cfg.signatures[name]) || {};
          const nCore = (sig.core || []).length;
          const few = Number.isFinite(minCore) && nCore && v.core_present < minCore;
          return h('tr', null,
            h('td', null, sig.color ? h('span', { class: 'swatch', style: { background: sig.color } }) : null, ' ', sig.label || name),
            h('td', null, h('div', { class: 'pill-list' }, (v.missing || []).map((g) => h('span', { class: 'chip no-x bad' }, g)))),
            h('td', { class: 'num', title: few ? `Fewer than scoring.min_core = ${minCore}: no bin can be called positive` : '' },
              few ? h('span', { class: 'state bad' }, icon('alert', { size: 12 }), `${v.core_present} of ${nCore}`) : `${v.core_present} of ${nCore}`));
        })))));
    } else {
      parts.push(h('p', { class: 'state ok' }, icon('check', { size: 13 }), 'Every signature gene is in the probe set.'));
    }
    const rec = (genes.receptor_panel && genes.receptor_panel.missing) || [];
    const adi = (genes.adipocyte_genes && genes.adipocyte_genes.missing) || [];
    for (const [title, miss] of [['Host receptor panel', rec], ['Adipocyte genes', adi]]) {
      parts.push(h('div', { class: 'report-genes' }, h('span', { class: 'report-genes-title' }, `${title}:`),
        miss.length ? h('div', { class: 'pill-list' }, miss.map((g) => h('span', { class: 'chip no-x bad', title: 'missing from the probe set' }, g)))
          : h('span', { class: 'state ok' }, icon('check', { size: 12 }), 'all present')));
    }
  }
  return h('div', { class: 'report' }, parts);
}
