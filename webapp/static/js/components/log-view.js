// log-view.js — live log of a run: incremental GET /api/runs/<id>/log?offset=..., line highlighting, auto-scroll that
// pauses while the reader is scrolled up, copy and download. The fetched text lives in a caller-owned `store`, so
// leaving the Progress tab and coming back neither refetches nor loses the log.
import { h, api, button, icon, copyText, toast, fmt } from '../core.js';

const MAX_DOM_CHARS = 1500000; // older lines leave the page beyond this (Download always has everything)
const BIG_CHUNK = 64 * 1024; // a response this large means more is waiting: ask again at once
const ERR_RE = /Traceback|[A-Za-z]*Error\b|[A-Za-z]*Exception\b/;
const WARN_RE = /WARNING|Warning:/;
const FOOTER_RE = /^--- (succeeded|failed|cancelled|interrupted)\b/;
const FOOTER_CLASS = { succeeded: 'l-step', failed: 'l-err', cancelled: 'l-warn', interrupted: 'l-err' };

/** Offer a Blob as a file download (no server round trip). */
export function downloadBlob(filename, blob) {
  const url = URL.createObjectURL(blob);
  const a = h('a', { href: url, download: filename, style: 'display:none' });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

/** Offer `text` as a file download. */
export function downloadText(filename, text, type = 'text/plain;charset=utf-8') {
  downloadBlob(filename, new Blob([text], { type }));
}

function utf8Length(text) {
  return new TextEncoder().encode(text).length;
}

/**
 * logView({runId, store, isFinished}) -> {el, mounted(), destroy()}
 * store: {text, offset, complete} — filled by the view; pass the same object again to resume.
 * isFinished(): whether the run has ended (a missing log is then final instead of "not written yet").
 */
export function logView({ runId, store, isFinished = () => false, interval = 1500 }) {
  if (typeof store.text !== 'string') store.text = '';
  if (!Number.isFinite(store.offset)) store.offset = 0;
  store.complete = !!store.complete;

  let destroyed = false;
  let timer = null;
  let ctrl = null;
  let failures = 0;
  let doneLen = 0; // characters of store.text already rendered as complete lines
  let inTraceback = false;
  let domChars = 0;

  const pre = h('pre', { class: 'log lv-log', tabIndex: 0, role: 'log', 'aria-label': 'Run log' });
  const tail = h('span', { class: 'lv-tail' }); // the last, still incomplete line
  pre.append(tail);
  const placeholder = h('div', { class: 'lv-empty' });
  const trimNote = h('div', { class: 'lv-trim', hidden: true },
    icon('info', { size: 14 }), 'Earlier lines are hidden to keep the page fast. Download has the complete log.');
  const jump = h('button', { class: 'btn small lv-jump', type: 'button', hidden: true, onClick: () => scrollToBottom() },
    icon('chevron-down', { size: 14 }), 'Jump to latest');
  const meta = h('span', { class: 'lv-meta' });
  const errorLine = h('span', { class: 'lv-error', role: 'status' });

  const el = h('section', { class: 'card lv' },
    h('div', { class: 'card-header' },
      h('div', { class: 'card-titles' }, h('h2', null, 'Log'), h('div', { class: 'card-sub' }, meta)),
      h('div', { class: 'row' }, errorLine,
        button({ label: 'Copy', icon: 'copy', small: true, onClick: copyLog }),
        button({ label: 'Download', icon: 'download', small: true, onClick: downloadLog }))),
    h('div', { class: 'lv-body' }, trimNote, pre, placeholder, jump));

  function lineClass(line) {
    if (line.startsWith('=== [')) { inTraceback = false; return 'l-step'; }
    const footer = FOOTER_RE.exec(line); // the job runner's last line: "--- succeeded (worker exit code 0)"
    if (footer) { inTraceback = false; return FOOTER_CLASS[footer[1]]; }
    if (line.startsWith('Traceback')) { inTraceback = true; return 'l-err'; }
    if (inTraceback) {
      if (/^\s/.test(line)) return 'l-err';
      inTraceback = false; // the exception line closes the traceback
      return 'l-err';
    }
    if (WARN_RE.test(line)) return 'l-warn';
    if (ERR_RE.test(line)) return 'l-err';
    return '';
  }

  /** Complete lines -> fragment; plain lines are merged into text nodes, highlighted lines get a span each. */
  function build(chunk) {
    const frag = document.createDocumentFragment();
    let plain = '';
    const lines = chunk.split('\n');
    lines.pop(); // the chunk ends with "\n"
    for (const line of lines) {
      const cls = lineClass(line);
      if (cls) {
        if (plain) { frag.append(plain); plain = ''; }
        frag.append(h('span', { class: cls }, `${line}\n`));
      } else {
        plain += `${line}\n`;
      }
    }
    if (plain) frag.append(plain);
    return frag;
  }

  function render() {
    const text = store.text;
    const lastNl = text.lastIndexOf('\n');
    if (lastNl >= doneLen) {
      let from = doneLen;
      if (lastNl + 1 - from > MAX_DOM_CHARS) {
        // a long log arriving at once (or shown again): only its end goes into the page
        from = text.indexOf('\n', lastNl + 1 - MAX_DOM_CHARS) + 1;
        while (pre.firstChild && pre.firstChild !== tail) pre.firstChild.remove();
        domChars = 0;
        trimNote.hidden = false;
      }
      const chunk = text.slice(from, lastNl + 1);
      pre.insertBefore(build(chunk), tail);
      domChars += chunk.length;
      doneLen = lastNl + 1;
    }
    const rest = text.slice(doneLen).replace(/\r$/, '');
    tail.textContent = rest;
    tail.className = `lv-tail ${rest ? lineClassPreview(rest) : ''}`.trim();
    while (domChars > MAX_DOM_CHARS && pre.firstChild && pre.firstChild !== tail) {
      domChars -= pre.firstChild.textContent.length;
      pre.firstChild.remove();
      trimNote.hidden = false;
    }
    paintMeta();
  }

  function lineClassPreview(line) {
    if (line.startsWith('=== [')) return 'l-step';
    if (WARN_RE.test(line)) return 'l-warn';
    return ERR_RE.test(line) || inTraceback ? 'l-err' : '';
  }

  /** Append a chunk to the store, normalizing CRLF / CR line ends (a CR at a chunk end may be half of a CRLF). */
  function ingest(chunk) {
    let t = chunk;
    if (store.text.endsWith('\r')) {
      store.text = store.text.slice(0, -1);
      if (!t.startsWith('\n')) t = `\n${t}`;
    }
    const pendingCR = t.endsWith('\r');
    if (pendingCR) t = t.slice(0, -1);
    store.text += t.replace(/\r\n?/g, '\n') + (pendingCR ? '\r' : '');
  }

  function paintMeta() {
    const finished = store.complete;
    meta.textContent = `${fmt.bytes(store.offset)}${finished ? ' · complete' : ' · live'}`;
    el.classList.toggle('lv-live', !finished);
    const empty = !store.text.replace(/\r$/, '');
    placeholder.hidden = !empty;
    if (empty) {
      placeholder.textContent = finished ? 'The log is empty.' : (isFinished() ? 'Reading the log…' : 'Waiting for output…');
    }
  }

  function atBottom() {
    return pre.scrollHeight - pre.scrollTop - pre.clientHeight < 24;
  }

  function scrollToBottom() {
    pre.scrollTop = pre.scrollHeight;
    jump.hidden = true;
  }

  pre.addEventListener('scroll', () => { jump.hidden = atBottom(); });

  function schedule(ms) {
    clearTimeout(timer);
    timer = null;
    if (!destroyed && !store.complete) timer = setTimeout(fetchMore, ms);
  }

  async function fetchMore() {
    if (destroyed || store.complete) return;
    ctrl = new AbortController();
    try {
      const res = await api(`/api/runs/${encodeURIComponent(runId)}/log`, { query: { offset: store.offset }, signal: ctrl.signal });
      if (destroyed) return;
      failures = 0;
      errorLine.textContent = '';
      const text = res && typeof res.text === 'string' ? res.text : '';
      // the response offset is the new read position; fall back to counting bytes if it did not advance
      store.offset = res && Number.isFinite(res.offset) && res.offset > store.offset ? res.offset : store.offset + utf8Length(text);
      store.complete = !!(res && res.complete);
      if (text) {
        const follow = atBottom();
        ingest(text);
        render();
        if (follow) scrollToBottom();
        else jump.hidden = false;
      } else {
        paintMeta();
      }
      schedule(text.length >= BIG_CHUNK ? 0 : interval);
    } catch (err) {
      if (destroyed || (err && err.name === 'AbortError')) return;
      if (err.status === 404 && isFinished()) {
        store.complete = true;
        paintMeta();
        placeholder.textContent = 'No log was written for this run.';
        return;
      }
      if (err.status === 404) {
        schedule(interval); // the log file appears when the worker starts
        return;
      }
      failures += 1;
      errorLine.textContent = `Log unavailable (${err.message}); retrying`;
      schedule(Math.min(interval * 2 ** failures, 15000));
    } finally {
      ctrl = null;
    }
  }

  async function copyLog() {
    const text = store.text.replace(/\r$/, '');
    if (!text) { toast('The log is empty', 'info', { timeout: 2000 }); return; }
    await copyText(text);
  }

  function downloadLog() {
    downloadText(`${runId}-log.txt`, store.text.replace(/\r$/, ''));
  }

  render();
  if (!store.complete) fetchMore();

  return {
    el,
    mounted() { scrollToBottom(); },
    destroy() {
      destroyed = true;
      clearTimeout(timer);
      if (ctrl) ctrl.abort();
    },
  };
}
