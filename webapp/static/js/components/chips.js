// chips.js — small input widgets of the project editor: chip inputs for gene lists, ordered chip multi-selects for
// signature names, and an icon-only button with an accessible name.
import { h, icon, mount, button } from '../core.js';

const SEPARATORS = /[\s,;]+/;
const QUOTES = /^["'“”‘’]+|["'“”‘’]+$/g;

/** Tokens of typed or pasted text: whitespace, commas and semicolons separate; surrounding quotes are dropped. */
export function parseTokens(text) {
  return String(text ?? '')
    .split(SEPARATORS)
    .map((t) => t.replace(QUOTES, '').trim())
    .filter(Boolean);
}

/** The distinct values of a list as strings, first occurrence kept. */
export function uniqueList(list) {
  return [...new Set((Array.isArray(list) ? list : []).map(String))];
}

/** chipSelect options for the signatures of a configuration, in configuration order: [{value, label, color}]. */
export function signatureOptions(cfg, only) {
  const keep = only ? new Set(only) : null;
  return Object.entries((cfg && cfg.signatures) || {})
    .filter(([name]) => !keep || keep.has(name))
    .map(([name, s]) => ({ value: name, label: (s && s.label) || name.replace(/_/g, ' '), color: s && s.color }));
}

/** "granulocyte (neutrophil)" when the label differs from the key, the label otherwise. */
export function optionLabel(o) {
  const plain = String(o.value).replace(/_/g, ' ');
  return o.label && o.label !== o.value && o.label !== plain ? `${o.label} (${o.value})` : o.label || o.value;
}

/** Icon-only button with an accessible name (button() of core.js only sets a title). */
export function iconButton(iconName, label, onClick, { kind = 'ghost', disabled = false, cls = '' } = {}) {
  const b = button({ icon: iconName, kind, small: true, title: label, disabled, cls, onClick });
  b.setAttribute('aria-label', label);
  return b;
}

function restartFlash(el) {
  el.classList.remove('flash');
  void el.offsetWidth; // reflow, so that the animation runs again
  el.classList.add('flash');
}

/**
 * Token input for gene names. Enter, Tab, comma, semicolon or space adds the typed text; a pasted list adds every
 * item; duplicates are skipped (the existing chip flashes); Backspace in the empty field removes the last chip;
 * Backspace or Delete on a chip's remove button removes that chip.
 * classify(value) -> '' | 'ok' | 'warn' | 'bad'; isCore(value) -> bold chip; describe(value) -> tooltip.
 * Returns {el, set(values), refresh(), values(), focus()}.
 */
export function chipInput({ values = [], onChange, classify, isCore, describe, label = 'Genes', placeholder = 'Add…', id } = {}) {
  let items = uniqueList(values);
  const nodes = new Map();
  const list = h('ul', { class: 'chip-list', 'aria-label': label });
  const entry = h('input', {
    class: 'chip-entry', type: 'text', id, placeholder, spellcheck: false, autocomplete: 'off',
    'aria-label': `Add to ${label}`,
  });
  const el = h('div', { class: 'chip-input' }, list, entry);
  el.addEventListener('mousedown', (e) => {
    if (e.target === el || e.target === list) { e.preventDefault(); entry.focus(); }
  });

  const emit = () => { if (onChange) onChange(items.slice()); };

  function node(v) {
    let li = nodes.get(v);
    if (!li) {
      const remove = h('button', { type: 'button', 'aria-label': `Remove ${v}` }, icon('x', { size: 12 }));
      remove.addEventListener('click', (e) => removeValue(v, e.detail === 0));
      remove.addEventListener('keydown', (e) => {
        if (e.key === 'Backspace' || e.key === 'Delete') { e.preventDefault(); removeValue(v, true); }
      });
      li = h('li', { dataset: { value: v } }, h('span', { class: 'chip-text' }, v), remove);
      li.addEventListener('animationend', () => li.classList.remove('flash'));
      nodes.set(v, li);
    }
    const kind = classify ? classify(v) || '' : '';
    li.className = ['chip', kind, isCore && isCore(v) ? 'core' : ''].filter(Boolean).join(' ');
    const tip = describe ? describe(v) : '';
    if (tip) li.title = tip; else li.removeAttribute('title');
    return li;
  }

  function paint() {
    for (const v of [...nodes.keys()]) if (!items.includes(v)) nodes.delete(v);
    list.replaceChildren(...items.map(node));
    el.classList.toggle('is-empty', items.length === 0);
  }

  function removeValue(v, keyboard) {
    const i = items.indexOf(v);
    if (i < 0) return;
    items.splice(i, 1);
    paint();
    emit();
    const next = items[i] ?? items[i - 1];
    const target = keyboard && next !== undefined ? nodes.get(next).querySelector('button') : entry;
    target.focus();
  }

  function add(tokens) {
    let added = false;
    for (const t of tokens) {
      if (items.includes(t)) {
        const li = nodes.get(t);
        if (li) restartFlash(li);
      } else {
        items.push(t);
        added = true;
      }
    }
    if (added) { paint(); emit(); }
  }

  function commit() {
    const tokens = parseTokens(entry.value);
    entry.value = '';
    if (tokens.length) add(tokens);
  }

  entry.addEventListener('keydown', (e) => {
    if (e.key === 'Enter' || e.key === ',' || e.key === ';' || e.key === ' ') {
      if (e.key !== 'Enter' || entry.value.trim()) e.preventDefault();
      commit();
    } else if (e.key === 'Tab') {
      if (entry.value.trim()) commit();
    } else if (e.key === 'Backspace' && entry.value === '' && items.length) {
      e.preventDefault();
      items.pop();
      paint();
      emit();
    } else if (e.key === 'Escape' && entry.value) {
      e.stopPropagation();
      entry.value = '';
    }
  });
  entry.addEventListener('paste', (e) => {
    const text = e.clipboardData ? e.clipboardData.getData('text') : '';
    if (!SEPARATORS.test(text)) return;
    e.preventDefault();
    const { value } = entry;
    const start = entry.selectionStart ?? value.length;
    const end = entry.selectionEnd ?? value.length;
    entry.value = '';
    add(parseTokens(value.slice(0, start) + text + value.slice(end)));
  });
  entry.addEventListener('blur', commit);

  paint();
  return {
    el,
    set(next) { items = uniqueList(next); paint(); },
    refresh: paint,
    values: () => items.slice(),
    focus: () => entry.focus(),
  };
}

/**
 * Ordered multi-select of known values shown as chips (with a color swatch), plus a select that adds one more.
 * options: [{value, label, color}]. Selected values missing from the options are kept and shown as bad chips with
 * `unknownText`. Returns {el, set(values)}.
 */
export function chipSelect({ options = [], values = [], onChange, label = 'Selection', placeholder = 'Add…', emptyText = 'None selected', unknownText = 'not available here', id } = {}) {
  let items = uniqueList(values);
  const el = h('div', { class: 'chip-select' });
  const emit = () => { if (onChange) onChange(items.slice()); };

  function paint(focusAdd) {
    const byValue = new Map(options.map((o) => [o.value, o]));
    const chips = items.map((v) => {
      const o = byValue.get(v);
      const text = o ? o.label || v : v;
      return h('li', { class: o ? 'chip' : 'chip bad', title: o ? optionLabel(o) : `${v}: ${unknownText}` },
        o && o.color ? h('span', { class: 'swatch', style: { background: o.color } }) : null,
        h('span', { class: 'chip-text' }, text),
        h('button', {
          type: 'button', 'aria-label': `Remove ${text}`,
          onClick: () => { items = items.filter((x) => x !== v); paint(true); emit(); },
        }, icon('x', { size: 12 })));
    });
    const rest = options.filter((o) => !items.includes(o.value));
    const add = h('select', { class: 'select chip-add', id, 'aria-label': `Add to ${label}`, disabled: rest.length === 0 },
      h('option', { value: '' }, rest.length ? placeholder : 'All selected'),
      rest.map((o) => h('option', { value: o.value }, optionLabel(o))));
    add.addEventListener('change', () => {
      if (!add.value) return;
      items.push(add.value);
      paint(true);
      emit();
    });
    mount(el, h('ul', { class: 'chip-list', 'aria-label': label },
      chips.length ? chips : h('li', { class: 'chip-empty' }, emptyText)), add);
    if (focusAdd) add.focus();
  }

  paint(false);
  return { el, set(next) { items = uniqueList(next); paint(false); } };
}
