// figure-gallery.js — the Figures tab of a run: the figures of figures/manifest.json grouped as in the manuscript,
// with lazily loaded thumbnails and a full-size viewer (fit or 1:1, PNG / PDF download, previous / next with the
// arrow keys).
import { h, mount, icon, modal, banner } from '../core.js';

function downloadUrl(url) {
  const u = new URL(url, window.location.origin);
  u.searchParams.set('download', '1');
  return u.pathname + u.search;
}

function downloadButton(url, label) {
  if (!url) return null;
  return h('a', { class: 'btn small', href: downloadUrl(url), download: '', title: `Download ${label}` },
    icon('download', { size: 13 }), label);
}

function brokenImage(text) {
  return h('div', { class: 'fg-broken' }, icon('image', { size: 22 }), h('span', null, text));
}

/** figureGallery({figures}) -> {el, destroy()}; figures = the "figures" list of GET /api/runs/<id>/manifest. */
export function figureGallery({ figures }) {
  const groups = new Map();
  for (const f of Array.isArray(figures) ? figures : []) {
    if (!f || !(f.name || f.title)) continue;
    const g = f.group || 'Other figures';
    if (!groups.has(g)) groups.set(g, []);
    groups.get(g).push(f);
  }
  const ordered = [...groups.values()].flat(); // the viewer steps through the figures in gallery order
  let viewer = null;

  function card(f, index) {
    const title = f.title || f.name;
    const thumb = h('button', { class: 'fg-thumb', type: 'button', 'aria-label': `Open ${title}`, onClick: () => open(index, thumb) });
    if (f.error || !f.png_url) {
      thumb.classList.add('is-broken');
      mount(thumb, h('div', { class: 'fg-broken bad' }, icon('alert', { size: 22 }),
        h('span', null, f.error ? 'This figure could not be drawn' : 'No image was written')));
    } else {
      thumb.classList.add('is-loading');
      const img = h('img', { loading: 'lazy', decoding: 'async', alt: title, src: f.png_url });
      img.addEventListener('load', () => { thumb.classList.remove('is-loading'); thumb.classList.add('loaded'); });
      img.addEventListener('error', () => {
        thumb.classList.remove('is-loading');
        thumb.classList.add('is-broken');
        mount(thumb, brokenImage('Image not available'));
      });
      thumb.append(img);
    }
    return h('article', { class: f.error ? 'fg-card has-error' : 'fg-card' },
      thumb,
      h('div', { class: 'fg-info' },
        h('div', { class: 'fg-title' }, title),
        f.error ? h('div', { class: 'fg-error' }, String(f.error)) : null,
        !f.error && f.caption ? h('div', { class: 'fg-caption' }, f.caption) : null,
        h('div', { class: 'fg-links' }, downloadButton(f.png_url, 'PNG'), downloadButton(f.pdf_url, 'PDF'))));
  }

  function open(start, returnFocus) {
    if (viewer) viewer.close();
    let i = start;
    const stage = h('div', { class: 'fg-stage fit' });
    const caption = h('div', { class: 'fg-modal-caption' });
    const counter = h('span', { class: 'fg-counter' });
    const links = h('div', { class: 'row' });
    const navButton = (iconName, title, step) => h('button', {
      class: 'btn small icon-only', type: 'button', title, 'aria-label': title, onClick: () => show(i + step),
    }, icon(iconName, { size: 15 }));
    const prev = navButton('chevron-left', 'Previous figure (left arrow)', -1);
    const next = navButton('chevron-right', 'Next figure (right arrow)', 1);
    const zoom = h('button', { class: 'btn small', type: 'button', 'aria-pressed': 'false', onClick: () => toggleZoom() }, 'Actual size');
    const body = h('div', { class: 'fg-modal' },
      h('div', { class: 'fg-modal-bar' }, prev, counter, next, h('span', { class: 'spacer' }), zoom, links),
      stage, caption);

    const onKey = (e) => {
      if (e.defaultPrevented || e.altKey || e.ctrlKey || e.metaKey) return;
      if (e.target instanceof Element && e.target.closest('input, textarea, select')) return;
      if (e.key === 'ArrowLeft') { e.preventDefault(); show(i - 1); }
      if (e.key === 'ArrowRight') { e.preventDefault(); show(i + 1); }
    };
    const m = modal({
      title: '', body, size: 'full',
      onClose: () => {
        document.removeEventListener('keydown', onKey);
        viewer = null;
        if (returnFocus && returnFocus.isConnected) returnFocus.focus();
      },
    });
    viewer = m;
    document.addEventListener('keydown', onKey);
    const heading = m.el.querySelector('.modal-header h2');

    function setZoom(actual) {
      stage.classList.toggle('actual', actual);
      stage.classList.toggle('fit', !actual);
      zoom.textContent = actual ? 'Fit to window' : 'Actual size';
      zoom.setAttribute('aria-pressed', actual ? 'true' : 'false');
    }

    function toggleZoom() {
      setZoom(!stage.classList.contains('actual'));
    }

    function show(k) {
      if (k < 0 || k >= ordered.length) return;
      i = k;
      const f = ordered[i];
      const title = f.title || f.name;
      if (heading) heading.textContent = title;
      counter.textContent = `${i + 1} of ${ordered.length}`;
      prev.disabled = i === 0;
      next.disabled = i === ordered.length - 1;
      setZoom(false);
      if (f.error || !f.png_url) {
        zoom.disabled = true;
        mount(stage, h('div', { class: 'fg-stage-msg' }, banner('danger', [h('strong', null, 'This figure could not be drawn. '),
          f.error ? String(f.error) : 'No image was written.'])));
      } else {
        zoom.disabled = false;
        const img = h('img', { alt: title, src: f.png_url });
        img.addEventListener('error', () => { zoom.disabled = true; mount(stage, brokenImage('Image not available')); });
        mount(stage, h('button', { class: 'fg-zoom', type: 'button', title: 'Click to switch between fit and actual size', onClick: toggleZoom }, img));
      }
      mount(caption, f.caption ? String(f.caption) : null);
      mount(links, downloadButton(f.png_url, 'PNG'), downloadButton(f.pdf_url, 'PDF'));
      const ahead = ordered[i + 1];
      if (ahead && ahead.png_url && !ahead.error) new Image().src = ahead.png_url; // warm the cache for "next"
    }

    show(i);
    m.el.tabIndex = -1;
    setTimeout(() => { if (m.el.isConnected) m.el.focus(); }, 40);
  }

  const failed = ordered.filter((f) => f.error).length;
  let index = 0;
  const el = h('div', { class: 'fg' },
    failed ? banner('warning', `${failed} figure${failed === 1 ? '' : 's'} could not be drawn; the card${failed === 1 ? ' shows' : 's show'} the reason. The other figures are complete.`) : null,
    [...groups.entries()].map(([name, figs]) => h('section', { class: 'fg-group' },
      h('h3', { class: 'fg-group-title' }, name, h('span', { class: 'fg-group-count' }, String(figs.length))),
      h('div', { class: 'fg-grid' }, figs.map((f) => card(f, index++))))));

  return {
    el,
    destroy() {
      if (viewer) viewer.close();
    },
  };
}
