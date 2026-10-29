// folder-picker.js — modal folder browser over GET /api/fs/list (directories only). Used for the data root and for
// sample folders; resolves to the chosen absolute path, or null when the dialog is cancelled.
import { h, api, modal, icon, mount, loading, banner, button } from '../core.js';

let pickerSeq = 0;

/** Breadcrumb segments of an absolute path, root first: [{label, path}]. Handles Windows, UNC and POSIX paths. */
export function pathCrumbs(path) {
  const p = String(path || '');
  if (!p) return [];
  const win = /^[A-Za-z]:([\\/]|$)/.test(p) || p.startsWith('\\\\');
  if (!win) {
    const out = [{ label: '/', path: '/' }];
    let acc = '';
    for (const part of p.split('/').filter(Boolean)) {
      acc += `/${part}`;
      out.push({ label: part, path: acc });
    }
    return out;
  }
  const norm = p.replace(/\//g, '\\');
  let root;
  let rest;
  if (norm.startsWith('\\\\')) {
    const parts = norm.slice(2).split('\\').filter(Boolean);
    root = `\\\\${parts.slice(0, 2).join('\\')}\\`;
    rest = parts.slice(2);
  } else {
    root = `${norm.slice(0, 2)}\\`;
    rest = norm.slice(2).split('\\').filter(Boolean);
  }
  const out = [{ label: root.replace(/\\$/, ''), path: root }];
  let acc = root;
  for (const part of rest) {
    acc = acc.endsWith('\\') ? acc + part : `${acc}\\${part}`;
    out.push({ label: part, path: acc });
  }
  return out;
}

/**
 * Folder chooser. start: folder shown first (absolute, relative to the server's working folder, or '' for the home
 * folder); when it cannot be listed the home folder is shown instead. Resolves to an absolute path or null.
 */
export function openFolderPicker({ title = 'Choose a folder', start = '', note = '', selectLabel = 'Select this folder' } = {}) {
  return new Promise((resolve) => {
    const uid = `fp-${++pickerSeq}`;
    let current = null;
    let seq = 0;
    let answered = false;

    const pathInput = h('input', {
      class: 'input mono', id: `${uid}-path`, spellcheck: false, autocomplete: 'off', placeholder: 'Type a folder path and press Enter',
    });
    const rootSelect = h('select', { class: 'select fp-roots', 'aria-label': 'Drive or root folder' });
    const upBtn = button({ icon: 'chevron-up', small: true, title: 'Parent folder (Alt+Up)', onClick: () => { if (current && current.parent) go(current.parent); } });
    upBtn.setAttribute('aria-label', 'Parent folder');
    const goBtn = button({ label: 'Go', small: true, onClick: () => { const v = pathInput.value.trim(); if (v) go(v); } });
    const crumbs = h('nav', { class: 'fp-crumbs', 'aria-label': 'Current folder' });
    const filter = h('input', { class: 'input fp-filter', type: 'search', placeholder: 'Filter folders', 'aria-label': 'Filter folders', autocomplete: 'off' });
    const list = h('ul', { class: 'fp-list', 'aria-label': 'Subfolders' });
    const status = h('div', { class: 'fp-status' });
    const message = h('div', { class: 'fp-message' });

    rootSelect.addEventListener('change', () => { if (rootSelect.value) go(rootSelect.value); });
    pathInput.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') { e.preventDefault(); const v = pathInput.value.trim(); if (v) go(v); }
    });
    filter.addEventListener('input', paintList);

    const body = h('div', { class: 'fp', onKeydown: (e) => {
      if (e.altKey && e.key === 'ArrowUp' && current && current.parent) { e.preventDefault(); go(current.parent); }
    } },
    note ? h('p', { class: 'muted small fp-note' }, note) : null,
    h('div', { class: 'fp-bar' }, rootSelect, upBtn, h('div', { class: 'input-group fp-path' }, pathInput, goBtn)),
    crumbs, message,
    h('div', { class: 'fp-listwrap' }, filter, list),
    status);

    const m = modal({
      title, body, size: 'wide',
      onClose: () => { if (!answered) { answered = true; resolve(null); } },
      actions: [
        { label: 'Cancel', onClick: () => { answered = true; resolve(null); } },
        {
          label: selectLabel, kind: 'primary', icon: 'check',
          onClick: () => {
            if (!current) return false;
            answered = true;
            resolve(current.path);
            return true;
          },
        },
      ],
    });
    const selectBtn = m.el.querySelector('.modal-footer .btn.primary');

    function paintList() {
      const q = filter.value.trim().toLowerCase();
      const all = (current && current.entries) || [];
      const shown = all.filter((e) => !q || String(e.name).toLowerCase().includes(q));
      if (!shown.length) {
        mount(list, h('li', { class: 'fp-empty' }, all.length ? 'No folder matches the filter.' : 'This folder has no subfolders.'));
        return;
      }
      mount(list, shown.map((e) => h('li', null, h('button', {
        type: 'button', class: 'fp-entry', title: e.path, onClick: () => go(e.path),
      }, icon('folder', { size: 15 }), h('span', { class: 'fp-name' }, e.name),
      e.looks_like_spaceranger ? h('span', { class: 'badge fp-sr', title: 'Contains outs/, binned_outputs/ or square_*um/ folders' }, 'Space Ranger') : null,
      icon('chevron-right', { size: 14, cls: 'fp-go' })))));
    }

    function paint() {
      const roots = (current && current.roots) || [];
      const crumbList = pathCrumbs(current ? current.path : '');
      const rootPath = crumbList.length ? crumbList[0].path : '';
      mount(rootSelect, roots.map((r) => h('option', { value: r, selected: r.toLowerCase() === rootPath.toLowerCase() }, r)));
      rootSelect.classList.toggle('hidden', roots.length < 2);
      upBtn.disabled = !(current && current.parent);
      pathInput.value = current ? current.path : '';
      mount(crumbs, crumbList.map((c, i) => [
        i ? icon('chevron-right', { size: 12, cls: 'fp-sep' }) : null,
        i === crumbList.length - 1
          ? h('span', { class: 'fp-crumb current', 'aria-current': 'location' }, c.label)
          : h('button', { type: 'button', class: 'fp-crumb', onClick: () => go(c.path) }, c.label),
      ]));
      const entries = (current && current.entries) || [];
      const sr = entries.filter((e) => e.looks_like_spaceranger).length;
      mount(status, `${entries.length} folder${entries.length === 1 ? '' : 's'}`,
        sr ? h('span', { class: 'muted' }, ` · ${sr} look${sr === 1 ? 's' : ''} like Space Ranger output`) : null);
      if (selectBtn) selectBtn.disabled = !current;
      paintList();
    }

    async function go(path, fallbackNote) {
      const mySeq = ++seq;
      mount(message);
      if (!current) mount(list, h('li', null, loading('Reading folders…')));
      else list.classList.add('busy');
      try {
        const res = await api('/api/fs/list', { query: { path } });
        if (mySeq !== seq || answered) return;
        const first = !current;
        current = res;
        filter.value = '';
        if (fallbackNote) mount(message, banner('warning', fallbackNote));
        paint();
        list.scrollTop = 0;
        if (first && (!body.contains(document.activeElement) || document.activeElement === rootSelect)) filter.focus();
      } catch (err) {
        if (mySeq !== seq || answered) return;
        if (!current && path) {
          go('', `Could not open “${path}” (${err.message}). Showing your home folder.`);
          return;
        }
        mount(message, banner('danger', `Could not open “${path || 'home folder'}”: ${err.message}`));
        if (!current) mount(list, h('li', { class: 'fp-empty' }, 'No folder to show.'));
      } finally {
        if (mySeq === seq) list.classList.remove('busy');
      }
    }

    if (selectBtn) selectBtn.disabled = true;
    go(String(start || '').trim());
  });
}
