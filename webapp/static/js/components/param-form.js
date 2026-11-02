// param-form.js — the Parameters tab: one card per PARAM_SCHEMA section (a control per field type, inline invalid
// state, reset to the default), then proximity pairs, figure options, analysis steps and reproducibility.
// Every edit changes the draft configuration in place and calls onChange().
import { h, mount, icon, card, switchToggle, button, clone } from '../core.js';
import { chipInput, chipSelect, signatureOptions, optionLabel, iconButton, parseTokens, uniqueList } from './chips.js';

let seq = 0;
const nextId = (prefix) => `${prefix}-${++seq}`;
const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
const PANEL_NAME_RE = /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;

const FIGURE_FIELDS = [
  { key: 'dpi', label: 'Figure resolution (dpi)', type: 'int', min: 50, max: 1200, help: 'Resolution of the PNG files.' },
  { key: 'map_dpi', label: 'Spatial map resolution (dpi)', type: 'int', min: 50, max: 1200, help: 'Resolution of the per-bin maps (Figures 3a and 4a–f).' },
  { key: 'moran_top_n', label: 'Genes shown in Supplementary Figure S2', type: 'int', min: 1, max: 500 },
  { key: 'receptor_shade_first_two', label: 'Shade the first two receptor genes (Figure 6b)', type: 'bool' },
];
const REPRO_FIELDS = [
  { key: 'seed', label: 'Random seed', type: 'int', min: 0, max: 2147483647, help: 'Seeds every random step: subsampling, null gene sets, Leiden and random placement.' },
  { key: 'deterministic', label: 'Deterministic run', type: 'bool', help: 'Single-threaded numerics and PYTHONHASHSEED=0, so that a rerun reproduces every number (slower).' },
];

/** "0, 5, 10" -> [0, 5, 10]; null when an item is not a whole number. */
export function parseIntList(text) {
  const out = [];
  for (const t of parseTokens(text)) {
    if (!/^[-+]?\d+$/.test(t)) return null;
    out.push(Number(t));
  }
  return out;
}

/** {ok, value, msg} for a number input: ok = the value can be stored; msg = what is wrong (range errors keep ok). */
export function checkNumber(raw, spec, badInput = false) {
  const text = String(raw ?? '').trim();
  if (text === '') return { ok: false, msg: badInput ? 'Not a number' : 'Required' };
  const v = Number(text);
  if (!Number.isFinite(v)) return { ok: false, msg: 'Not a number' };
  if (spec.type === 'int' && !Number.isInteger(v)) return { ok: false, msg: 'A whole number is required' };
  let msg = '';
  if (typeof spec.min === 'number' && v < spec.min) msg = `Must be at least ${spec.min}`;
  else if (typeof spec.max === 'number' && v > spec.max) msg = `Must be at most ${spec.max}`;
  return { ok: true, value: v, msg };
}

const numText = (v) => (v === null || v === undefined ? '' : String(v));
const listText = (v) => (Array.isArray(v) ? v.join(', ') : '');

function describeValue(v) {
  let text;
  if (typeof v === 'boolean') text = v ? 'on' : 'off';
  else if (Array.isArray(v)) text = v.length ? v.join(', ') : 'empty';
  else text = String(v);
  return text.length > 70 ? `${text.slice(0, 67)}…` : text;
}

function obj(parent, key) {
  if (!parent[key] || typeof parent[key] !== 'object' || Array.isArray(parent[key])) parent[key] = {};
  return parent[key];
}

/** Accessor of cfg[section][key] (section null = top level) with its default. */
function accessor(ctx, section, key) {
  // a section may be nested: "sensitivity_analyses.uncoupling" is cfg.sensitivity_analyses.uncoupling
  const parts = section ? section.split('.') : [];
  const parent = () => parts.reduce((o, p) => obj(o, p), ctx.cfg);
  const defParent = parts.reduce((o, p) => (o && o[p]) || {}, ctx.defaults);
  return { get: () => parent()[key], set: (v) => { parent()[key] = v; }, def: defParent[key], path: section ? `${section}.${key}` : key };
}

/** One labelled control for a schema field; returns {el}. */
function makeField(spec, acc, ctx) {
  const id = nextId('pf');
  const helpId = spec.help ? `${id}-help` : null;
  const err = h('div', { class: 'error-text hidden', id: `${id}-err`, role: 'alert' });
  const hasDefault = acc.def !== undefined;
  const reset = h('button', {
    type: 'button', class: 'pf-reset hidden', title: hasDefault ? `Reset to the default (${describeValue(acc.def)})` : '',
    'aria-label': `Reset “${spec.label}” to the default`,
  }, icon('undo', { size: 13 }), h('span', null, 'Default'));
  const sync = () => reset.classList.toggle('hidden', !hasDefault || same(acc.get(), acc.def));
  const write = (v) => { acc.set(v); ctx.onChange(); sync(); };
  const setError = (msg, control) => {
    err.textContent = msg || '';
    err.classList.toggle('hidden', !msg);
    if (control) {
      control.classList.toggle('invalid', !!msg);
      control.setAttribute('aria-invalid', msg ? 'true' : 'false');
    }
  };
  let control;
  let show;
  let labelled = true;

  switch (spec.type) {
    case 'int':
    case 'float': {
      const input = h('input', {
        class: 'input pf-num', type: 'number', id, min: spec.min, max: spec.max,
        step: spec.step ?? (spec.type === 'int' ? 1 : 'any'), inputmode: spec.type === 'int' ? 'numeric' : 'decimal',
      });
      input.addEventListener('input', () => {
        const r = checkNumber(input.value, spec, input.validity.badInput);
        setError(r.msg, input);
        if (r.ok) write(r.value);
      });
      input.addEventListener('change', () => {
        if (!checkNumber(input.value, spec, input.validity.badInput).ok) show();
      });
      show = () => {
        input.value = numText(acc.get());
        const r = checkNumber(input.value, spec);
        setError(r.ok ? r.msg : '', input);
      };
      control = input;
      break;
    }
    case 'bool': {
      labelled = false;
      control = switchToggle({ checked: !!acc.get(), label: spec.label, onChange: (v) => write(v) });
      const box = control.querySelector('input');
      box.id = id;
      show = () => { box.checked = !!acc.get(); };
      break;
    }
    case 'select': {
      const sel = h('select', { class: 'select', id });
      show = () => {
        const cur = acc.get();
        const options = spec.options || [];
        mount(sel, [...options, ...(options.includes(cur) || cur === undefined || cur === null ? [] : [cur])]
          .map((o) => h('option', { value: o, selected: o === cur }, options.includes(o) ? o : `${o} (not supported)`)));
        sel.classList.toggle('invalid', !options.includes(cur));
      };
      sel.addEventListener('change', () => { write(sel.value); show(); });
      control = sel;
      break;
    }
    case 'int_list': {
      const input = h('input', { class: 'input mono', id, placeholder: 'e.g. 0, 5, 10', spellcheck: false, autocomplete: 'off' });
      input.addEventListener('input', () => {
        const list = parseIntList(input.value);
        if (list === null) setError('Whole numbers separated by commas', input);
        else if (!list.length) setError('At least one value', input);
        else { setError('', input); write(list); }
      });
      input.addEventListener('change', () => show());
      show = () => { input.value = listText(acc.get()); setError('', input); };
      control = input;
      break;
    }
    case 'str_list': {
      const ci = chipInput({
        values: acc.get() || [], id, label: spec.label, placeholder: 'Type and press Enter, or paste a list',
        classify: spec.classify, describe: spec.describe, onChange: (v) => write(v),
      });
      show = () => ci.set(acc.get() || []);
      control = ci.el;
      break;
    }
    case 'signature': {
      const sel = h('select', { class: 'select', id });
      show = () => {
        const cur = acc.get();
        const options = signatureOptions(ctx.cfg);
        const known = options.some((o) => o.value === cur);
        mount(sel, known ? null : h('option', { value: cur ?? '', selected: true }, cur ? `${cur} (unknown signature)` : '(none)'),
          options.map((o) => h('option', { value: o.value, selected: o.value === cur }, optionLabel(o))));
        sel.classList.toggle('invalid', !known);
        sel.title = known ? '' : 'Not a signature of this project: the step is skipped';
      };
      sel.addEventListener('change', () => { write(sel.value); show(); });
      control = sel;
      break;
    }
    case 'signature_list': {
      const cs = chipSelect({
        options: signatureOptions(ctx.cfg), values: acc.get() || [], id, label: spec.label, placeholder: 'Add a signature…',
        unknownText: 'not a signature of this project', onChange: (v) => write(v),
      });
      show = () => cs.set(acc.get() || []);
      control = cs.el;
      break;
    }
    default: {
      const input = h('input', { class: 'input', id, spellcheck: false });
      input.addEventListener('input', () => write(input.value));
      show = () => { input.value = acc.get() ?? ''; };
      control = input;
    }
  }

  const describedBy = [helpId, `${id}-err`].filter(Boolean).join(' ');
  const focusable = control.matches('input, select') ? control : control.querySelector(`#${id}`);
  if (focusable) focusable.setAttribute('aria-describedby', describedBy);
  reset.addEventListener('click', () => {
    acc.set(clone(acc.def));
    ctx.onChange();
    show();
    sync();
  });
  show();
  sync();
  const wide = spec.type === 'signature_list' || spec.type === 'str_list';
  const head = labelled
    ? h('div', { class: 'pf-head' }, h('label', { for: id, title: acc.path }, spec.label), reset)
    : h('div', { class: 'pf-head', title: acc.path }, control, reset);
  return {
    el: h('div', { class: `field pf-field${wide ? ' pf-wide' : ''}` }, head, labelled ? control : null, err,
      spec.help ? h('div', { class: 'help', id: helpId }, spec.help) : null),
  };
}

function section(title, subtitle, body, actions) {
  return { title, el: card({ title, subtitle, body, actions, cls: 'pf-card' }) };
}

// ------------------------------------------------------------------------------------------------ extra cards
function pairsSection(ctx) {
  const body = h('div', { class: 'stack' });
  const pairs = () => {
    const prox = obj(ctx.cfg, 'proximity');
    if (!Array.isArray(prox.pairs)) prox.pairs = [];
    return prox.pairs;
  };
  const sigSelect = (p, key, label, onPick) => {
    const options = signatureOptions(ctx.cfg);
    const known = options.some((o) => o.value === p[key]);
    const sel = h('select', { class: `select${known ? '' : ' invalid'}`, 'aria-label': label },
      known ? null : h('option', { value: p[key] ?? '', selected: true }, p[key] ? `${p[key]} (unknown)` : '(choose)'),
      options.map((o) => h('option', { value: o.value, selected: o.value === p[key] }, optionLabel(o))));
    sel.addEventListener('change', () => { p[key] = sel.value; sel.classList.remove('invalid'); ctx.onChange(); onPick(); });
    return sel;
  };
  const paint = (focusAdd) => {
    const list = pairs();
    const rows = list.map((p, i) => {
      const file = h('td', { class: 'mono small muted' });
      const paintFile = () => { file.textContent = `${p.from || '?'}_to_${p.to || '?'}_proximity.csv`; };
      paintFile();
      const mode = h('select', { class: 'select', 'aria-label': `Pair ${i + 1}: counted by` },
        [['assignment', 'tissue assignment'], ['positivity', 'positivity calls']].map(([v, t]) => h('option', { value: v, selected: p.mode === v }, t)),
        ['assignment', 'positivity'].includes(p.mode) ? null : h('option', { value: p.mode ?? '', selected: true }, `${p.mode} (not supported)`));
      mode.addEventListener('change', () => { p.mode = mode.value; ctx.onChange(); });
      return h('tr', null,
        h('td', null, sigSelect(p, 'from', `Pair ${i + 1}: from`, paintFile)),
        h('td', { class: 'center muted' }, icon('chevron-right', { size: 14 })),
        h('td', null, sigSelect(p, 'to', `Pair ${i + 1}: to`, paintFile)),
        h('td', null, mode),
        file,
        h('td', { class: 'right' }, iconButton('trash', `Remove pair ${i + 1}`, () => { list.splice(i, 1); ctx.onChange(); paint(true); })));
    });
    const add = button({
      label: 'Add pair', icon: 'plus', small: true,
      onClick: () => {
        const names = Object.keys(ctx.cfg.signatures || {});
        list.push({ from: names[0] || '', to: names[1] || names[0] || '', mode: 'assignment' });
        ctx.onChange();
        paint();
      },
    });
    mount(body, list.length ? h('div', { class: 'table-wrap' }, h('table', { class: 'table compact pf-pairs' },
      h('thead', null, h('tr', null, h('th', null, 'From'), h('th', null, h('span', { class: 'sr-only' }, 'to')), h('th', null, 'To'),
        h('th', null, 'Counted by'), h('th', null, 'Output table'), h('th', null, h('span', { class: 'sr-only' }, 'Remove')))),
      h('tbody', null, rows))) : h('p', { class: 'muted' }, 'No proximity pairs: no proximity table is written.'),
    h('div', { class: 'row' }, add));
    if (focusAdd) add.focus();
  };
  paint();
  return section('Proximity pairs', 'One table per pair: how many “from” bins have a “to” bin within the radius set under '
    + 'Proximity tables, counted by the tissue assignment or by the positivity calls.', body);
}

function panelsEditor(ctx, onPanelsChange) {
  const box = h('div', { class: 'pf-panels' });
  const panels = () => {
    const fig = obj(ctx.cfg, 'figures');
    if (!Array.isArray(fig.positivity_panels)) fig.positivity_panels = [];
    return fig.positivity_panels;
  };
  const changed = () => { ctx.onChange(); onPanelsChange(); };
  const paint = (focusName) => {
    const list = panels();
    const inputs = [];
    const check = () => {
      const counts = new Map();
      list.forEach((p) => counts.set(p.name, (counts.get(p.name) || 0) + 1));
      inputs.forEach(({ input, p }) => {
        const msg = !PANEL_NAME_RE.test(String(p.name ?? '')) ? 'Letters, digits, “_”, “.” or “-” (a file name)'
          : counts.get(p.name) > 1 ? 'Another panel has this name' : '';
        input.classList.toggle('invalid', !!msg);
        input.title = msg;
        input.setAttribute('aria-invalid', msg ? 'true' : 'false');
      });
    };
    const cards = list.map((p, i) => {
      const input = h('input', { class: 'input mono', value: p.name ?? '', spellcheck: false, autocomplete: 'off', 'aria-label': `Positivity panel ${i + 1}: figure name` });
      input.addEventListener('input', () => { p.name = input.value.trim(); changed(); check(); });
      inputs.push({ input, p });
      const sigs = chipSelect({
        options: signatureOptions(ctx.cfg), values: p.signatures || [], label: `Signatures of panel ${i + 1}`,
        placeholder: 'Add a signature…', unknownText: 'not a signature of this project',
        onChange: (v) => { p.signatures = v; changed(); },
      });
      return h('div', { class: 'pf-panel' },
        h('div', { class: 'pf-panel-head' }, input,
          iconButton('trash', `Remove positivity panel ${p.name || i + 1}`, () => { list.splice(i, 1); changed(); paint(); })),
        sigs.el);
    });
    const add = button({
      label: 'Add panel', icon: 'plus', small: true,
      onClick: () => {
        const names = new Set(list.map((p) => p.name));
        let n = list.length + 1;
        while (names.has(`Figure_positivity_${n}`)) n += 1;
        list.push({ name: `Figure_positivity_${n}`, signatures: [] });
        changed();
        paint(true);
      },
    });
    mount(box, cards.length ? cards : h('p', { class: 'muted' }, 'No positivity panels.'), h('div', { class: 'row' }, add));
    check();
    if (focusName && inputs.length) inputs[inputs.length - 1].input.focus();
  };
  paint();
  return { el: box, paint };
}

function figuresSection(ctx) {
  const fig = (key) => accessor(ctx, 'figures', key);
  const lineage = (key) => ({
    get: () => obj(obj(ctx.cfg, 'figures'), 'lineage_maps')[key],
    set: (v) => { obj(obj(ctx.cfg, 'figures'), 'lineage_maps')[key] = v; },
    def: ((ctx.defaults.figures || {}).lineage_maps || {})[key],
    path: `figures.lineage_maps.${key}`,
  });
  const panelsDefault = (ctx.defaults.figures || {}).positivity_panels;
  const panelsReset = button({
    label: 'Default panels', icon: 'undo', small: true, kind: 'ghost',
    onClick: () => { obj(ctx.cfg, 'figures').positivity_panels = clone(panelsDefault || []); ctx.onChange(); panels.paint(); syncPanelsReset(); },
  });
  const syncPanelsReset = () => panelsReset.classList.toggle('hidden', same(obj(ctx.cfg, 'figures').positivity_panels, panelsDefault));
  const panels = panelsEditor(ctx, syncPanelsReset);
  syncPanelsReset();
  const receptorSpec = {
    key: 'receptor_genes', label: 'Receptor genes shown in Figure 6b', type: 'str_list',
    help: 'Genes of the host receptor panel (Signatures tab), in drawing order.',
    classify: (g) => (ctx.cfg.receptor_panel && Object.prototype.hasOwnProperty.call(ctx.cfg.receptor_panel, g) ? '' : 'warn'),
    describe: (g) => (ctx.cfg.receptor_panel && Object.prototype.hasOwnProperty.call(ctx.cfg.receptor_panel, g)
      ? ctx.cfg.receptor_panel[g] : 'Not in the host receptor panel: the figure has no data for this gene'),
  };
  const body = h('div', { class: 'stack-lg' },
    h('div', { class: 'stack' },
      h('div', { class: 'row between' }, h('h3', null, 'Positivity panels'), panelsReset),
      h('p', { class: 'muted small' }, 'One figure per panel (defaults: Figures 4g and 6a): the percentage of positive bins of each signature, per section.'),
      panels.el),
    h('div', { class: 'form-grid pf-grid' },
      makeField(receptorSpec, fig('receptor_genes'), ctx).el,
      makeField({ key: 'composition_signatures', label: 'Signatures in the composition bars (Figure 3b)', type: 'signature_list' }, fig('composition_signatures'), ctx).el,
      makeField({ key: 'a', label: 'Figure 4a–f: first signature (blue)', type: 'signature' }, lineage('a'), ctx).el,
      makeField({ key: 'b', label: 'Figure 4a–f: second signature (red)', type: 'signature' }, lineage('b'), ctx).el,
      FIGURE_FIELDS.map((f) => makeField(f, fig(f.key), ctx).el)));
  return section('Figures', 'What bone_figures.py draws after the analysis.', body);
}

function stepsSection(ctx, stepLabels) {
  const grid = h('div', { class: 'pf-steps' });
  const defaults = ctx.defaults.steps || {};
  const restore = button({
    label: 'Restore defaults', icon: 'undo', small: true, kind: 'ghost',
    onClick: () => { ctx.cfg.steps = clone(defaults); ctx.onChange(); paint(); },
  });
  const paint = () => {
    const steps = obj(ctx.cfg, 'steps');
    const keys = uniqueList([...Object.keys(defaults), ...Object.keys(steps)]);
    mount(grid, keys.map((k) => h('div', { class: 'pf-step', title: `steps.${k}` }, switchToggle({
      checked: steps[k] === true, label: stepLabels[k] || k,
      onChange: (v) => { obj(ctx.cfg, 'steps')[k] = v; ctx.onChange(); syncRestore(); },
    }))));
    syncRestore();
  };
  const syncRestore = () => restore.classList.toggle('hidden', same(obj(ctx.cfg, 'steps'), defaults));
  paint();
  return section('Analysis steps', 'Steps that are switched off are skipped; their tables and figures are not written. '
    + 'The Run tab can switch the slow steps for a single run.', grid, restore);
}

function reproSection(ctx) {
  return section('Reproducibility', null, h('div', { class: 'form-grid pf-grid' },
    REPRO_FIELDS.map((f) => makeField(f, accessor(ctx, null, f.key), ctx).el)));
}

/**
 * The whole Parameters tab. schema = PARAM_SCHEMA, defaults = the default configuration, stepLabels = STEP_LABELS.
 * Returns {el}.
 */
export function paramForm({ cfg, defaults, schema, stepLabels, onChange }) {
  const ctx = { cfg, defaults: defaults || {}, onChange: onChange || (() => {}) };
  const sections = (schema || []).map((sec) => section(sec.title, sec.description,
    h('div', { class: 'form-grid pf-grid' }, (sec.fields || []).map((f) => makeField(f, accessor(ctx, sec.key, f.key), ctx).el))));
  sections.push(pairsSection(ctx), figuresSection(ctx), stepsSection(ctx, stepLabels || {}), reproSection(ctx));
  const nav = h('nav', { class: 'pf-nav', 'aria-label': 'Parameter sections' },
    sections.map((s) => h('button', {
      type: 'button', class: 'pf-nav-item',
      onClick: () => {
        s.el.scrollIntoView({ behavior: 'smooth', block: 'start' });
        const first = s.el.querySelector('input, select, button.btn');
        if (first) first.focus({ preventScroll: true });
      },
    }, s.title)));
  return { el: h('div', { class: 'pf-layout' }, nav, h('div', { class: 'pf-cards' }, sections.map((s) => s.el))) };
}
