// markdown.js — a small, safe Markdown renderer for the documentation pages: CommonMark blocks and inlines plus the
// GitHub extensions the project documents use (tables, strikethrough, task lists, bare URLs). It only ever creates
// elements and text nodes, so HTML in a document is shown as text and never interpreted. Images are shown as links
// (the app makes no network requests on behalf of a document).
import { h } from '../core.js';

const FENCE_RE = /^( {0,3})(`{3,}|~{3,})(.*)$/;
const ATX_RE = /^ {0,3}(#{1,6})(?:[ \t]+(.*?))?(?:[ \t]+#+)?[ \t]*$/;
const HR_RE = /^ {0,3}([-*_])(?:[ \t]*\1){2,}[ \t]*$/;
const QUOTE_RE = /^ {0,3}> ?(.*)$/;
const SETEXT_RE = /^ {0,3}(=+|-+)[ \t]*$/;
const LIST_RE = /^( {0,3})([-*+]|(\d{1,9})([.)]))(?:([ \t]+)(.*))?$/;
const REF_DEF_RE = /^ {0,3}\[([^\]]+)\]:[ \t]*<?([^\s>]+)>?(?:[ \t]+(?:"([^"]*)"|'([^']*)'|\(([^)]*)\)))?[ \t]*$/;
const COMMENT_RE = /^ {0,3}<!--/;
const TASK_RE = /^\[([ xX])\][ \t]+/;
const PUNCT_RE = /[!-/:-@[-`{-~]/;
const ENTITY_RE = /&(#[0-9]{1,7}|#[xX][0-9a-fA-F]{1,6}|[A-Za-z][A-Za-z0-9]{1,31});/y;
const AUTOLINK_RE = /<((?:https?:\/\/|mailto:)[^\s<>]+)>/iy;
const EMAIL_RE = /<([^\s<>@]+@[^\s<>@]+\.[^\s<>@]+)>/y;
const BARE_URL_RE = /https?:\/\/[^\s<>]+/iy;

const NAMED_ENTITIES = {
  amp: '&', lt: '<', gt: '>', quot: '"', apos: "'", nbsp: ' ', ensp: ' ', emsp: ' ', thinsp: ' ',
  mdash: '—', ndash: '–', hellip: '…', copy: '©', reg: '®', trade: '™', times: '×',
  divide: '÷', minus: '−', plusmn: '±', le: '≤', ge: '≥', ne: '≠', asymp: '≈',
  larr: '←', rarr: '→', uarr: '↑', darr: '↓', harr: '↔', rArr: '⇒', micro: 'µ',
  deg: '°', middot: '·', bull: '•', sect: '§', para: '¶', laquo: '«', raquo: '»',
  lsquo: '‘', rsquo: '’', ldquo: '“', rdquo: '”', alpha: 'α', beta: 'β',
  gamma: 'γ', delta: 'δ', mu: 'μ', sigma: 'σ', check: '✓',
};

// ------------------------------------------------------------------------------------------------ helpers
const isBlank = (line) => /^\s*$/.test(line);
const isSpace = (ch) => ch === undefined || /\s/.test(ch);
const isAlnum = (ch) => ch !== undefined && /[\p{L}\p{N}]/u.test(ch);
const indentOf = (line) => /^ */.exec(line)[0].length;
const normLabel = (label) => label.trim().toLowerCase().replace(/\s+/g, ' ');

/** Leading tabs become spaces (tab stops of 4), so indentation can be measured in columns. */
function expandTabs(line) {
  if (!line.includes('\t')) return line;
  let out = '';
  let i = 0;
  for (; i < line.length; i += 1) {
    if (line[i] === ' ') out += ' ';
    else if (line[i] === '\t') out += ' '.repeat(4 - (out.length % 4));
    else break;
  }
  return out + line.slice(i);
}

function runLength(s, i, ch) {
  let n = 0;
  while (s[i + n] === ch) n += 1;
  return n;
}

/** Index of the next run of exactly `n` backticks at or after `from`, or -1. */
function findBacktickRun(s, from, n) {
  let j = from;
  while (j < s.length) {
    const k = s.indexOf('`', j);
    if (k < 0) return -1;
    const r = runLength(s, k, '`');
    if (r === n) return k;
    j = k + r;
  }
  return -1;
}

/** Index just past the code span that starts at s[i] (a backtick), or i + run length when it is not closed. */
function skipCodeSpan(s, i) {
  const n = runLength(s, i, '`');
  const close = findBacktickRun(s, i + n, n);
  return close < 0 ? i + n : close + n;
}

export function slugify(text) {
  return String(text).trim().toLowerCase().replace(/[^\p{L}\p{N}\s_-]/gu, '').replace(/\s/g, '-');
}

function decodeEntity(body) {
  if (body[0] === '#') {
    const hex = body[1] === 'x' || body[1] === 'X';
    const code = hex ? parseInt(body.slice(2), 16) : parseInt(body.slice(1), 10);
    if (!Number.isFinite(code) || code <= 0 || code > 0x10ffff || (code >= 0xd800 && code <= 0xdfff)) return '�';
    return String.fromCodePoint(code);
  }
  return Object.prototype.hasOwnProperty.call(NAMED_ENTITIES, body) ? NAMED_ENTITIES[body] : null;
}

function isFence(m) {
  return !!m && !(m[2][0] === '`' && m[3].includes('`'));
}

function isClosingFence(line, marker) {
  const m = FENCE_RE.exec(line);
  return !!m && m[2][0] === marker[0] && m[2].length >= marker.length && !m[3].trim();
}

function listItemStart(line) {
  if (HR_RE.test(line)) return null;
  const m = LIST_RE.exec(line);
  if (!m) return null;
  const ordered = m[3] !== undefined;
  const spaces = m[5] || '';
  const rest = m[6] ?? '';
  const empty = !rest.trim();
  const pad = empty || spaces.length > 4 ? 1 : spaces.length;
  return {
    ordered,
    kind: ordered ? m[4] : m[2],
    start: ordered ? Number(m[3]) : 1,
    width: m[1].length + m[2].length + pad,
    content: empty ? '' : (spaces.length > 4 ? spaces.slice(1) + rest : rest),
  };
}

/** Lines that end a paragraph (CommonMark: an ordered list interrupts one only when it starts at 1). */
function interrupts(line) {
  if (ATX_RE.test(line) || HR_RE.test(line) || QUOTE_RE.test(line) || COMMENT_RE.test(line)) return true;
  if (isFence(FENCE_RE.exec(line))) return true;
  const li = listItemStart(line);
  return !!li && li.content.trim() !== '' && (!li.ordered || li.start === 1);
}

// ------------------------------------------------------------------------------------------------ tables
/** Cells of a table row: split on pipes outside code spans; "\|" is a literal pipe. */
function splitRow(line) {
  let s = line.trim();
  if (s.startsWith('|')) s = s.slice(1);
  if (s.endsWith('|') && !s.endsWith('\\|')) s = s.slice(0, -1);
  const cells = [];
  let cur = '';
  let i = 0;
  while (i < s.length) {
    const c = s[i];
    if (c === '\\' && s[i + 1] === '|') { cur += '|'; i += 2; continue; }
    if (c === '`') { const end = skipCodeSpan(s, i); cur += s.slice(i, end); i = end; continue; }
    if (c === '|') { cells.push(cur.trim()); cur = ''; i += 1; continue; }
    cur += c;
    i += 1;
  }
  cells.push(cur.trim());
  return cells;
}

function isDelimiterRow(line) {
  if (!line.includes('|') || !line.includes('-') || !/^[\s|:-]+$/.test(line)) return false;
  return splitRow(line).every((c) => /^:?-+:?$/.test(c));
}

function isTableStart(lines, i) {
  if (i + 1 >= lines.length || !lines[i].includes('|') || !isDelimiterRow(lines[i + 1])) return false;
  return splitRow(lines[i]).length === splitRow(lines[i + 1]).length;
}

// ------------------------------------------------------------------------------------------------ blocks
function parseBlocks(lines, ctx) {
  const out = [];
  let i = 0;
  while (i < lines.length) {
    const line = lines[i];
    if (isBlank(line)) { i += 1; continue; }
    let r;
    const fence = FENCE_RE.exec(line);
    const atx = ATX_RE.exec(line);
    if (COMMENT_RE.test(line)) {
      i = skipComment(lines, i);
      continue;
    } else if (isFence(fence)) {
      r = fencedCode(lines, i, fence);
    } else if (atx) {
      r = { node: heading(atx[1].length, atx[2] || '', ctx), next: i + 1 };
    } else if (HR_RE.test(line)) {
      r = { node: h('hr'), next: i + 1 };
    } else if (QUOTE_RE.test(line)) {
      r = blockquote(lines, i, ctx);
    } else if (listItemStart(line)) {
      r = list(lines, i, ctx);
    } else if (indentOf(line) >= 4) {
      r = indentedCode(lines, i);
    } else if (isTableStart(lines, i)) {
      r = table(lines, i, ctx);
    } else {
      r = paragraph(lines, i, ctx);
    }
    out.push(r.node);
    i = r.next;
  }
  return out;
}

function skipComment(lines, i) {
  for (let j = i; j < lines.length; j += 1) {
    if (lines[j].includes('-->') && (j > i || lines[j].indexOf('-->') > lines[j].indexOf('<!--'))) return j + 1;
  }
  return lines.length;
}

function heading(level, text, ctx) {
  const el = h(`h${level}`, null, parseInline(text.trim(), ctx));
  const base = slugify(el.textContent) || 'section';
  const seen = ctx.ids.get(base) || 0;
  ctx.ids.set(base, seen + 1);
  el.id = `${ctx.idPrefix}${seen ? `${base}-${seen}` : base}`;
  return el;
}

function fencedCode(lines, i, m) {
  const indent = m[1].length;
  const marker = m[2];
  const lang = (m[3].trim().split(/\s+/)[0] || '').replace(/[^\w+#.-]/g, '');
  const strip = new RegExp(`^ {0,${indent}}`);
  const body = [];
  let j = i + 1;
  while (j < lines.length && !isClosingFence(lines[j], marker)) {
    body.push(lines[j].replace(strip, ''));
    j += 1;
  }
  const node = h('pre', null, h('code', { class: lang ? `language-${lang}` : null }, body.join('\n')));
  return { node, next: Math.min(j + 1, lines.length) };
}

function indentedCode(lines, i) {
  const body = [];
  let j = i;
  while (j < lines.length && (isBlank(lines[j]) || indentOf(lines[j]) >= 4)) {
    body.push(lines[j].slice(Math.min(4, indentOf(lines[j]))));
    j += 1;
  }
  while (body.length && isBlank(body[body.length - 1])) body.pop();
  return { node: h('pre', null, h('code', null, body.join('\n'))), next: j };
}

function blockquote(lines, i, ctx) {
  const inner = [];
  let j = i;
  while (j < lines.length) {
    const line = lines[j];
    const m = QUOTE_RE.exec(line);
    if (m) {
      inner.push(m[1]);
    } else if (!isBlank(line) && inner.length && !isBlank(inner[inner.length - 1]) && !interrupts(line)) {
      inner.push(line); // lazy continuation of a quoted paragraph
    } else {
      break;
    }
    j += 1;
  }
  return { node: h('blockquote', null, parseBlocks(inner, ctx)), next: j };
}

/** True when a blank line separates two blocks at the top level of a list item (which makes the list loose). */
function hasBlockGap(body) {
  let fence = null;
  for (let k = 0; k < body.length; k += 1) {
    const line = body[k];
    if (fence) {
      if (isClosingFence(line, fence)) fence = null;
      continue;
    }
    const f = FENCE_RE.exec(line);
    if (isFence(f)) { fence = f[2]; continue; }
    if (isBlank(line) && k > 0 && !isBlank(body[k - 1]) && k + 1 < body.length && !isBlank(body[k + 1])
      && !listItemStart(body[k + 1]) && indentOf(body[k + 1]) < 4) return true;
  }
  return false;
}

function list(lines, i, ctx) {
  const first = listItemStart(lines[i]);
  const items = [];
  let loose = false;
  let j = i;
  while (j < lines.length) {
    const st = listItemStart(lines[j]);
    if (!st || st.ordered !== first.ordered || st.kind !== first.kind) break;
    const body = [st.content];
    let k = j + 1;
    while (k < lines.length) {
      const line = lines[k];
      if (isBlank(line)) {
        body.push('');
      } else if (indentOf(line) >= st.width) {
        body.push(line.slice(st.width));
      } else if (isBlank(body[body.length - 1]) || listItemStart(line) || interrupts(line)) {
        break;
      } else {
        body.push(line.trimStart()); // lazy continuation line
      }
      k += 1;
    }
    let trailingBlank = false;
    while (body.length > 1 && isBlank(body[body.length - 1])) { body.pop(); trailingBlank = true; }
    if (hasBlockGap(body)) loose = true;
    items.push(body);
    j = k;
    if (trailingBlank && j < lines.length) {
      const next = listItemStart(lines[j]);
      if (next && next.ordered === first.ordered && next.kind === first.kind) loose = true;
    }
  }
  const el = h(first.ordered ? 'ol' : 'ul', { start: first.ordered && first.start !== 1 ? String(first.start) : null });
  for (const body of items) {
    const task = TASK_RE.exec(body[0]);
    const lines2 = task ? [body[0].slice(task[0].length), ...body.slice(1)] : body;
    const li = h('li', { class: task ? 'md-task' : null },
      task ? h('input', { type: 'checkbox', checked: task[1] !== ' ', disabled: true, 'aria-label': task[1] !== ' ' ? 'done' : 'open' }) : null);
    for (const block of parseBlocks(lines2, ctx)) {
      if (!loose && block.tagName === 'P') li.append(...block.childNodes);
      else li.append(block);
    }
    el.append(li);
  }
  return { node: el, next: j };
}

function table(lines, i, ctx) {
  const head = splitRow(lines[i]);
  const aligns = splitRow(lines[i + 1]).map((c) => {
    const left = c.startsWith(':');
    const right = c.endsWith(':');
    if (left && right) return 'center';
    if (right) return 'right';
    return left ? 'left' : null;
  });
  const rows = [];
  let j = i + 2;
  while (j < lines.length && !isBlank(lines[j]) && lines[j].includes('|') && !interrupts(lines[j])) {
    rows.push(splitRow(lines[j]));
    j += 1;
  }
  const cell = (tag, text, k) => h(tag, { style: aligns[k] ? `text-align:${aligns[k]}` : null }, parseInline(text, ctx));
  const node = h('div', { class: 'md-table-wrap' }, h('table', null,
    h('thead', null, h('tr', null, head.map((c, k) => cell('th', c, k)))),
    rows.length ? h('tbody', null, rows.map((r) => h('tr', null, head.map((_, k) => cell('td', r[k] ?? '', k))))) : null));
  return { node, next: j };
}

function paragraph(lines, i, ctx) {
  const buf = [lines[i].trimStart()];
  let j = i + 1;
  while (j < lines.length) {
    const line = lines[j];
    if (isBlank(line)) break;
    const setext = SETEXT_RE.exec(line);
    if (setext) return { node: heading(setext[1][0] === '=' ? 1 : 2, buf.join('\n'), ctx), next: j + 1 };
    if (interrupts(line) || isTableStart(lines, j)) break;
    buf.push(line.trimStart());
    j += 1;
  }
  return { node: h('p', null, parseInline(buf.join('\n').replace(/[ \t]+$/, ''), ctx)), next: j };
}

/** Reference definitions ("[label]: url") are collected and removed from the text, outside fenced code. */
function collectReferences(lines, ctx) {
  const out = [];
  let fence = null;
  for (const line of lines) {
    if (fence) {
      if (isClosingFence(line, fence)) fence = null;
      out.push(line);
      continue;
    }
    const f = FENCE_RE.exec(line);
    if (isFence(f)) {
      fence = f[2];
      out.push(line);
      continue;
    }
    const m = REF_DEF_RE.exec(line);
    if (m) {
      const label = normLabel(m[1]);
      if (!ctx.refs.has(label)) ctx.refs.set(label, { url: m[2], title: m[3] ?? m[4] ?? m[5] ?? null });
      continue;
    }
    out.push(line);
  }
  return out;
}

// ------------------------------------------------------------------------------------------------ inlines
function safeHref(url, ctx) {
  const u = String(url || '').trim();
  if (!u) return null;
  if (ctx.resolveLink) {
    const r = ctx.resolveLink(u);
    if (r !== undefined) return r;
  }
  if (/^(https?:\/\/|mailto:)/i.test(u)) return u;
  if (u.startsWith('#') && !u.startsWith('#/')) return u;
  return null;
}

function makeLink(url, content, title, ctx, cls) {
  const href = safeHref(url, ctx);
  if (!href) return h('span', { class: ['md-link-text', cls], title: url }, content);
  if (href.startsWith('#') && !href.startsWith('#/')) {
    let target = href.slice(1);
    try { target = decodeURIComponent(target); } catch { /* keep the raw fragment */ }
    const id = `${ctx.idPrefix}${slugify(target)}`;
    return h('a', {
      href: ctx.anchorHref, title, class: cls,
      onClick: (e) => {
        e.preventDefault();
        const el = document.getElementById(id);
        if (el) el.scrollIntoView({ behavior: 'smooth', block: 'start' });
      },
    }, content);
  }
  const external = /^https?:/i.test(href);
  return h('a', { href, title, class: cls, target: external ? '_blank' : null, rel: external ? 'noopener noreferrer' : null }, content);
}

function matchBracket(s, i) {
  let depth = 0;
  for (let j = i; j < s.length; j += 1) {
    const ch = s[j];
    if (ch === '\\') { j += 1; continue; }
    if (ch === '`') { j = skipCodeSpan(s, j) - 1; continue; }
    if (ch === '[') depth += 1;
    else if (ch === ']') {
      depth -= 1;
      if (depth === 0) return j;
    }
  }
  return -1;
}

function skipWhitespace(s, j) {
  let k = j;
  while (k < s.length && (s[k] === ' ' || s[k] === '\t' || s[k] === '\n')) k += 1;
  return k;
}

/** Inline link destination and optional title starting after "(", up to and including ")". */
function parseDestination(s, start) {
  let j = skipWhitespace(s, start);
  let url;
  if (s[j] === '<') {
    const end = s.indexOf('>', j + 1);
    if (end < 0 || s.slice(j + 1, end).includes('\n')) return null;
    url = s.slice(j + 1, end);
    j = end + 1;
  } else {
    const from = j;
    let depth = 0;
    while (j < s.length) {
      const ch = s[j];
      if (ch === '\\' && j + 1 < s.length) { j += 2; continue; }
      if (ch === '(') depth += 1;
      else if (ch === ')') {
        if (depth === 0) break;
        depth -= 1;
      } else if (/\s/.test(ch)) break;
      j += 1;
    }
    url = s.slice(from, j);
  }
  j = skipWhitespace(s, j);
  let title = null;
  const q = s[j];
  if (q === '"' || q === "'" || q === '(') {
    const end = s.indexOf(q === '(' ? ')' : q, j + 1);
    if (end < 0) return null;
    title = s.slice(j + 1, end);
    j = skipWhitespace(s, end + 1);
  }
  if (s[j] !== ')') return null;
  return { url: url.replace(/\\([!-/:-@[-`{-~])/g, '$1'), title, end: j + 1 };
}

/** [text](url "title"), [text][ref], [text][] or [ref] at s[i]; image = the "!" form. */
function parseLink(s, i, ctx, image) {
  if (ctx.inLink) return null;
  const close = matchBracket(s, i);
  if (close < 0) return null;
  const text = s.slice(i + 1, close);
  let dest;
  let end;
  if (s[close + 1] === '(') {
    const d = parseDestination(s, close + 2);
    if (!d) return null;
    dest = d;
    end = d.end;
  } else {
    let label = text;
    end = close + 1;
    if (s[close + 1] === '[') {
      const c2 = s.indexOf(']', close + 2);
      if (c2 < 0) return null;
      if (s.slice(close + 2, c2).trim()) label = s.slice(close + 2, c2);
      end = c2 + 1;
    }
    dest = ctx.refs.get(normLabel(label));
    if (!dest) return null;
  }
  const content = image
    ? document.createTextNode(text.trim() ? `Image: ${text.trim()}` : 'Image')
    : parseInline(text, { ...ctx, inLink: true });
  return { node: makeLink(dest.url, content, dest.title || null, ctx, image ? 'md-image' : null), end };
}

/** Next closing delimiter run of `k` characters `c` at or after `from` (k = 3 accepts longer runs), or -1. */
function findCloser(s, from, c, k) {
  let j = from;
  while (j < s.length) {
    const ch = s[j];
    if (ch === '\\') { j += 2; continue; }
    if (ch === '`') { j = skipCodeSpan(s, j); continue; }
    if (ch === c) {
      const r = runLength(s, j, c);
      const fits = k === 3 ? r >= 3 : r === k;
      if (fits && j > from && !isSpace(s[j - 1]) && !(c === '_' && isAlnum(s[j + r]))) return k === 3 ? j + r - 3 : j;
      j += r;
      continue;
    }
    j += 1;
  }
  return -1;
}

function parseEmphasis(s, i, ctx) {
  const c = s[i];
  const n = runLength(s, i, c);
  if (c === '~') {
    if (n !== 2 || isSpace(s[i + 2])) return null;
    const close = findCloser(s, i + 2, '~', 2);
    return close < 0 ? null : { node: h('del', null, parseInline(s.slice(i + 2, close), ctx)), end: close + 2 };
  }
  if (isSpace(s[i + n]) || (c === '_' && isAlnum(s[i - 1]))) return null;
  const sizes = n >= 3 ? [3, 2, 1] : n === 2 ? [2, 1] : [1];
  for (const k of sizes) {
    const close = findCloser(s, i + k, c, k);
    if (close < 0) continue;
    const inner = parseInline(s.slice(i + k, close), ctx);
    let node;
    if (k === 1) node = h('em', null, inner);
    else if (k === 2) node = h('strong', null, inner);
    else node = h('strong', null, h('em', null, inner));
    return { node, end: close + k };
  }
  return null;
}

/** Trailing punctuation and unbalanced closing parentheses are not part of a bare URL. */
function trimUrl(url) {
  let u = url.replace(/[.,:;!?'"*_~]+$/, '');
  while (u.endsWith(')') && (u.match(/\)/g) || []).length > (u.match(/\(/g) || []).length) {
    u = u.slice(0, -1).replace(/[.,:;!?'"*_~]+$/, '');
  }
  return u;
}

function stickyMatch(re, s, i) {
  re.lastIndex = i;
  return re.exec(s);
}

function parseInline(src, ctx) {
  const s = String(src);
  const frag = document.createDocumentFragment();
  let buf = '';
  const push = (node) => {
    if (buf) { frag.append(buf); buf = ''; }
    frag.append(node);
  };
  let i = 0;
  while (i < s.length) {
    const c = s[i];
    let r = null;
    let m = null;
    if (c === '\\') {
      const next = s[i + 1];
      if (next === '\n') { push(h('br')); i += 2; continue; }
      if (next !== undefined && PUNCT_RE.test(next)) { buf += next; i += 2; continue; }
    } else if (c === '`') {
      const n = runLength(s, i, '`');
      const close = findBacktickRun(s, i + n, n);
      if (close < 0) { buf += '`'.repeat(n); i += n; continue; }
      let code = s.slice(i + n, close).replace(/\n/g, ' ');
      if (code.length > 2 && code.startsWith(' ') && code.endsWith(' ') && code.trim()) code = code.slice(1, -1);
      push(h('code', null, code));
      i = close + n;
      continue;
    } else if (c === '*' || c === '_' || (c === '~' && s[i + 1] === '~')) {
      r = parseEmphasis(s, i, ctx);
      if (!r) { const n = runLength(s, i, c); buf += c.repeat(n); i += n; continue; }
    } else if (c === '!' && s[i + 1] === '[') {
      r = parseLink(s, i + 1, ctx, true);
    } else if (c === '[') {
      r = parseLink(s, i, ctx, false);
    } else if (c === '<' && !ctx.inLink) {
      if ((m = stickyMatch(AUTOLINK_RE, s, i))) r = { node: makeLink(m[1], m[1], null, ctx), end: i + m[0].length };
      else if ((m = stickyMatch(EMAIL_RE, s, i))) r = { node: makeLink(`mailto:${m[1]}`, m[1], null, ctx), end: i + m[0].length };
    } else if ((c === 'h' || c === 'H') && !ctx.inLink && !isAlnum(s[i - 1]) && (m = stickyMatch(BARE_URL_RE, s, i))) {
      const url = trimUrl(m[0]);
      if (url.length > 8) r = { node: makeLink(url, url, null, ctx), end: i + url.length };
    } else if (c === '&' && (m = stickyMatch(ENTITY_RE, s, i))) {
      const decoded = decodeEntity(m[1]);
      if (decoded !== null) { buf += decoded; i += m[0].length; continue; }
    } else if (c === '\n') {
      const hard = / {2,}$/.test(buf);
      buf = buf.replace(/[ \t]+$/, '');
      if (hard) push(h('br'));
      else buf += '\n';
      i += 1;
      continue;
    }
    if (r) {
      push(r.node);
      i = r.end;
    } else {
      buf += c;
      i += 1;
    }
  }
  if (buf) frag.append(buf);
  return frag;
}

// ------------------------------------------------------------------------------------------------ entry point
/**
 * Markdown text -> DocumentFragment.
 * resolveLink(href) may map a link target to an app URL (return a string), forbid it (null) or leave the default
 * handling (undefined: http(s)/mailto links open in a new tab, "#heading" scrolls to that heading, anything else is
 * shown as plain text). Heading ids are `${idPrefix}${slug}`; in-document anchors use `anchorHref` as their href so
 * that the hash router is not triggered.
 */
export function renderMarkdown(source, { resolveLink, idPrefix = 'md-', anchorHref = '#' } = {}) {
  const lines = String(source ?? '').replace(/\r\n?/g, '\n').split('\n').map(expandTabs);
  const ctx = { refs: new Map(), ids: new Map(), resolveLink, idPrefix, anchorHref, inLink: false };
  const frag = document.createDocumentFragment();
  frag.append(...parseBlocks(collectReferences(lines, ctx), ctx));
  return frag;
}
