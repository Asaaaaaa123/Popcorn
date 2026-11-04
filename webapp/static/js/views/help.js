// help.js — the Guide: how to use Bone Spatial Studio and what the analysis terms mean, plus the project documents
// (CHANGES.md, docs/ARCHITECTURE.md, README.md from GET /api/docs/<name>) rendered by markdown.js, which never
// interprets HTML. Topics switch in place; the address bar follows.
import { h, api, mount, icon, card, loading, emptyState, errorPanel, setPathSilently } from '../core.js';
import { renderMarkdown } from '../components/markdown.js';

const TOPICS = [
  { id: 'guide', label: 'Guide', icon: 'help' },
  { id: 'changes', label: 'What was fixed', icon: 'check', doc: 'changes' },
  { id: 'architecture', label: 'Architecture', icon: 'layers', doc: 'architecture' },
  { id: 'readme', label: 'README', icon: 'file', doc: 'readme' },
];
const DOC_FILES = { 'readme.md': 'readme', 'changes.md': 'changes', 'architecture.md': 'architecture' };

/** Links between the project documents ("docs/ARCHITECTURE.md", "../README.md#usage") open the matching topic. */
function resolveLink(href) {
  const path = href.split(/[?#]/)[0].replace(/^(\.\.?\/)+/, '').toLowerCase();
  const name = path.split('/').pop();
  if (DOC_FILES[name] && (path === name || path === `docs/${name}`)) return `#/help/${DOC_FILES[name]}`;
  return undefined;
}

// ------------------------------------------------------------------------------------------------ the guide
const link = (href, text) => h('a', { href }, text);
const code = (text) => h('code', null, text);

const workflow = () => [
  {
    title: 'Create a project, or try the demo',
    body: [
      'A project holds one analysis setup: the sections, the parameters, the gene signatures and an optional hand-drawn region of interest. Create one on the ', link('#/projects', 'Projects'),
      ' page. The demo offered there generates two small synthetic tibiae (about a minute, once) and a project with fast settings, so the whole workflow can be explored without real data.',
    ],
  },
  {
    title: 'Point it at the Space Ranger binned outputs',
    body: [
      'In the ', h('em', null, 'Data'), ' tab, choose the folder that contains one folder per Visium HD section and detect the samples: the app finds every folder with ',
      code('square_*um/filtered_feature_bc_matrix.h5'), ' (the ', code('outs/binned_outputs'),
      ' of Space Ranger). Give each section an id, a label, a condition and a color, and pick the reference section: every fold change and ratio is another section over it. The analysis bin size is 16 µm by default; Supplementary Figure S1 uses the 8 µm bins.',
    ],
  },
  {
    title: 'Validate',
    body: [
      h('em', null, 'Validate'), ' checks every section without running the analysis: the bin folders, positions, scale factors and H&E image, the number of bins, and which genes of each signature are in the probe set. ',
      h('strong', null, 'Errors'), ' must be fixed before a run can start; ', h('strong', null, 'warnings'), ' point at things worth a look (a missing gene, an unusual parameter).',
    ],
  },
  {
    title: 'Look at the raw data',
    body: [
      'The ', h('em', null, 'Explore'), ' tab reads every section as Space Ranger wrote it, with the current settings (saved or not) and before any run: a table of the bins, genes and UMIs of each section, the median UMIs, genes and mitochondrial share per bin, and how many bins pass the quality control set under Parameters; maps of these per-bin numbers or of any gene (several are summed), in one, two or four linked panels with the H&E image underneath; and the genes of a section ranked by UMIs, each with ',
      h('em', null, 'Show on map'), ' and ', h('em', null, 'Overlay'), '. The ', h('em', null, 'Explore'), ' button next to a section in the Data tab opens it there.',
    ],
  },
  {
    title: 'Review the parameters and gene signatures',
    body: [
      'The ', h('em', null, 'Parameters'), ' tab covers quality control, normalization and Leiden domains, scoring and positivity, the region of interest, compartments, distances, proximity and Moran’s I; every field explains itself and the defaults reproduce the manuscript. The ',
      h('em', null, 'Signatures'), ' tab edits the gene sets of Table 1: scored genes, core genes, competing signatures and excluded genes with their reasons, which signatures get a positivity call and which compete in the tissue assignment.',
    ],
  },
  {
    title: 'Optionally draw a region of interest',
    body: [
      'In the ', h('em', null, 'ROI'), ' tab, draw polygons (holes allowed) on the H&E image of a section, then turn on ', h('em', null, 'Intersect with the hand-drawn region'),
      ' under Parameters → Region of interest. The drawn region is intersected with the transcriptional one; sections without a drawing are not restricted further.',
    ],
  },
  {
    title: 'Run',
    body: [
      'The ', h('em', null, 'Run'), ' tab starts the analysis with an optional label and a choice of steps (Moran’s I and the gene-set sensitivity are the slow ones). One analysis runs at a time; others wait in the queue on the ',
      link('#/runs', 'Runs'), ' page. The Progress tab of a run shows every step, warnings and the live log, and a run can be cancelled at any time.',
    ],
  },
  {
    title: 'Explore the results',
    body: [
      h('em', null, 'Summary'), ': bins per section through quality control and the region of interest, and the key numbers with their confidence intervals. ',
      h('em', null, 'Figures'), ': the manuscript figures as PNG and PDF. ',
      h('em', null, 'Tables'), ': every CSV output, searchable and sortable, plus the AnnData object for further work in Python; ',
      h('em', null, 'Export for Loupe Browser'), ' (above the tables) downloads the per-bin calls as category files for 10x Genomics Loupe Browser — the ', code('_008um'),
      ' file of a section goes into its ', code('outs/cloupe_008um.cloupe'), ', and README.txt in the zip says how to import it. ',
      h('em', null, 'Spatial viewer'), ': every bin at its position, colored by tissue class, assigned signature, Leiden domain, region of interest, a positivity call, the UMI count or the counts of any genes (numbers and counts are colored up to their 99th percentile, marked ≥ in the legend, so that a few extreme bins do not flatten the map); one, two or four panels side by side with linked pan and zoom and, for panels colored by the same gene or number, one shared color scale (Color scales → Shared); yes/no overlays (calls, region of interest, “gene detected”) on top and the H&E image underneath.',
    ],
  },
  {
    title: 'Keep it reproducible',
    body: [
      'Each run freezes its own copy of the configuration and of the drawn regions, so editing the project later never changes a finished run; ', h('em', null, 'Re-run'),
      ' repeats it from that snapshot. Deterministic mode (on by default) runs the numerical libraries single-threaded with fixed seeds, so a rerun reproduces every number. The Summary records the package versions, and the Configuration tab downloads ',
      code('config.json'), ' (which also drives the command-line scripts described in the README) or turns it into a new project.',
    ],
  },
];

// using the app for tissues and signatures other than the bone lesion of the manuscript
const otherTissues = () => [
  h('p', null, 'The defaults are the manuscript’s: mouse bone and marrow signatures (Table 1) and a few steps built around osteoclasts and osteoblasts in an A549 bone lesion. The rest — quality control, Leiden domains, signature scores and positivity calls, depth matching, the region of interest, the tissue assignment and composition, Moran’s I, the per-bin table, the summary and the maps — works with any gene signatures.'),
  h('ol', null,
    h('li', null, h('strong', null, 'Signatures tab. '), 'Delete the signatures you do not need and add your own: scored genes, core genes (at least scoring.min_core of them must be detected for a positive call), competing signatures, and the switches ', h('em', null, 'Positivity call'), ' and ', h('em', null, 'Tissue assignment'), '. Replace the tissue classes, the host receptor panel and the adipocyte genes, or leave them. After a data check the gene chips show which genes are in the probe set.'),
    h('li', null, h('strong', null, 'Configuration check. '), 'Settings that still name a deleted signature are listed with the field to change. The Figure 4h compartments and the Figure 4i distance need signatures with a positivity call: choose yours under Parameters, or switch the steps off under Parameters → Analysis steps (they are marked “bone study”); until then the analysis cannot start. Figures whose signature lists are empty (composition bars, positivity panels, the Figure 4a–f maps) are not drawn: choose their signatures under Parameters → Figures.'),
    h('li', null, h('strong', null, 'Steps of the bone study. '), 'The host receptor panel (Figure 6b), the adipocyte genes (Supplementary Figure S4), the osteoclast compartments (Figure 4h), the distance to the nearest partner (Figure 4i; osteoblast to osteoclast in the manuscript) and the three-gene-set sensitivity follow the manuscript; each has a switch under Analysis steps. With a receptor panel of your own, also choose the genes Figure 6b draws (Parameters → Figures). With other signatures, Figure 4h calls its two compartments high-<split signature> and low-<split signature> instead of growth plate and tumor-associated. Under the sensitivity analyses, empty the A549 lesion genes (Krt7 and Krt18 are ordinary epithelial genes elsewhere) and set the habitat classes of Figure 4i to your own tissue classes. Sections count as tumor sections when their condition is listed under Osteoclast compartments → Tumor conditions.'),
    h('li', null, h('strong', null, 'Gene names. '), 'They must match the probe set exactly: the defaults are mouse symbols (', code('Lyz2'), ', ', code('Ptprc'), '), so human data needs human symbols (', code('LYZ'), ', ', code('PTPRC'), '), also in the ubiquitous transcripts of the Leiden screen (Signatures tab). Mitochondrial genes are recognised by the prefix ', code('mt-'), ' in any case, so ', code('MT-'), ' works too.')),
];

const GLOSSARY = [
  ['Bin', 'Visium HD captures expression on a continuous grid of 2 × 2 µm squares; Space Ranger sums them into square bins of 8 and 16 µm (square_008um/, square_016um/). A bin is the unit of the analysis — a small patch of tissue that usually holds parts of several cells, not a cell.'],
  ['UMI', 'Unique molecular identifier: every captured transcript is counted once, so the UMIs of a bin are the transcripts detected in it, its capture depth. A gene is detected in a bin when it has at least one UMI there.'],
  ['Quality control', 'A bin enters the analysis when it detects more than qc.min_genes genes, has more than qc.min_counts UMIs and has less than qc.max_pct_mt percent of its UMIs on mitochondrial genes (names starting with mt-). The Explore tab shows how many bins pass with the current thresholds before anything runs.'],
  ['Region of interest', 'The bins the section-level numbers are computed on. Its transcriptional part drops bins scoring above the 90th percentile for skeletal muscle, nerve or smooth muscle (tissue around the bone); a region drawn on the H&E image can be intersected with it. Positivity rates, fold changes and distances are computed inside it.'],
  ['Depth matching', 'Sections differ in capture depth (UMIs per bin), and a deeper section detects more of every signature. Depth matching downsamples the raw counts of every bin above the median UMIs per bin of the shallowest section to that median (in every section; bins below it are kept), several times, and reports the mean (± SD) over these replicates, so rates are compared at equal depth.'],
  ['Null threshold', 'For each signature, random gene sets of the same size are scored the same way on a sample of bins. The 95th percentile of those null scores is the threshold: a bin above it scores higher than 95% of random gene sets would.'],
  ['Positivity call', 'A bin is positive for a signature when its score exceeds the null threshold, at least two of the signature’s core genes are detected in it, and no competing signature scores higher. Rates are reported with Wilson 95% confidence intervals.'],
  ['Core genes', 'The lineage-restricted genes of a signature (Table 1). Requiring at least two of them to be detected (the two-gene rule) keeps a single shared or ambient transcript from producing a positive call.'],
  ['Competing signature', 'A signature that shares transcripts with another, for example the erythroid and granulocyte signatures for osteoclasts. A bin where a competitor scores higher is not called positive.'],
  ['Tissue assignment', 'Figure 3: the counts of each bin are pooled with its 8 nearest neighbors and every program signature is scored; the bin takes the highest-scoring signature that clears its null threshold with at least two core genes detected, otherwise it stays unassigned (“ambiguous”). Its tissue class is that signature’s class. Deterministic: no graph, seed or resolution parameter.'],
  ['Observed/null distance ratio', 'Figure 4i: the mean distance from each osteoblast-positive bin to the nearest osteoclast-positive bin inside the region of interest, divided by the same mean when the same numbers of positive bins are placed at random in that region. Below 1 the partners are closer than chance, above 1 farther; sections are compared by the ratio of these ratios, because raw distances depend on how many positive bins there are.'],
  ['Leiden domain', 'An unsupervised cluster of bins: Leiden community detection on a nearest-neighbor graph of the principal components (Harmony-integrated across sections). Domains summarize the tissue (Supplementary Figure S5); scores, calls and tests do not depend on them.'],
  ['Unresolved domain', 'A Leiden domain whose top markers are mostly transcripts abundant in every cell (Actb, Tmsb4x, Eef1a1, …) has no interpretable identity; such domains are merged into “unresolved”. Domains of fewer than 100 bins become “unassigned”. The label is defined without reference to the tumor, but it is not tumor-free: the mouse probes of conserved genes also capture the human A549 tumor RNA (Krt7, Krt18), and Tmsb10 and Itm2b track the lesion, so in the tumor-bearing section the unresolved domain is strongly enriched for it. Its share differs between sections (FigureS5_unresolved_bins.csv gives it per section).'],
];

const SHORTCUTS = [
  ['Figures', 'Click a figure to enlarge it; ← and → step through the figures, a click on the image switches between fitting the window and actual size, Esc closes.'],
  ['Tables', 'Click a column header to sort (again to reverse); the search box filters every column; Esc clears it.'],
  ['Maps (Spatial viewer, Explore)', 'Drag to pan, scroll or double-click to zoom (Shift + double-click zooms out); with a map focused, the arrow keys pan, + and − zoom and 0 fits the section (every panel’s section while the views are linked). Click a legend entry to hide or show it, Alt + click to show only that entry. In a gene box, ↑ and ↓ pick a suggestion and Enter shows it; commas separate genes whose counts are summed.'],
  ['Log', 'The log follows new output until you scroll up; “Jump to latest” returns to the end.'],
];

function guideSection(id, title, ...content) {
  return h('section', { class: 'help-section' }, h('h2', { id: `guide-${id}` }, title), content);
}

function guide() {
  return h('article', { class: 'prose help-guide' },
    h('p', { class: 'help-lead' },
      'Bone Spatial Studio runs the Visium HD bone-lesion analysis of the manuscript on your own sections: configure the samples, parameters and gene signatures, optionally draw a region of interest on the H&E image, run the pipeline in the background, and browse the summary, figures, tables and per-bin maps.'),
    h('div', { class: 'row help-actions' },
      h('a', { class: 'btn primary', href: '#/projects' }, icon('folder', { size: 16 }), 'Go to projects'),
      h('a', { class: 'btn', href: '#/runs' }, icon('activity', { size: 16 }), 'See runs')),
    guideSection('workflow', 'Workflow',
      h('ol', { class: 'help-steps' }, workflow().map((s) => h('li', null, h('h3', null, s.title), h('p', null, s.body))))),
    guideSection('other-tissues', 'Other tissues and gene signatures', otherTissues()),
    guideSection('glossary', 'Glossary',
      h('dl', { class: 'help-glossary' }, GLOSSARY.map(([term, text]) => [h('dt', null, term), h('dd', null, text)]))),
    guideSection('shortcuts', 'Mouse and keyboard',
      h('dl', { class: 'help-glossary' }, SHORTCUTS.map(([term, text]) => [h('dt', null, term), h('dd', null, text)]))),
    guideSection('storage', 'Where things are stored',
      h('p', null, 'Everything lives in the workspace folder (', code('workspace/'), ' next to the application unless it was started with ', code('--workspace'),
        '): ', code('projects/<id>/'), ' holds a project and its drawn regions; ', code('runs/<id>/'), ' holds a run: ', code('run.json'), ' (status and steps), ',
        code('config.json'), ' (the frozen configuration), ', code('rois/'), ', ', code('log.txt'), ' and ', code('results/'), ' with the tables, ', code('summary.json'),
        ' and ', code('figures/'), '. Deleting a project keeps its runs.')));
}

// ------------------------------------------------------------------------------------------------ view
export async function renderHelp(container, { topic }) {
  let destroyed = false;
  let token = 0;
  let ctrl = null;
  let current = TOPICS.some((t) => t.id === topic) ? topic : 'guide';
  if (topic && topic !== current) setPathSilently('/help');
  const docs = new Map();

  const nav = h('nav', { class: 'help-nav', 'aria-label': 'Help topics' });
  const toc = h('div', { class: 'help-toc' });
  const content = h('div', { class: 'help-content', 'aria-live': 'polite' });
  mount(container, h('div', { class: 'page help-page' },
    h('div', { class: 'page-header' }, h('div', { class: 'titles' }, h('h1', null, 'Help'),
      h('div', { class: 'subtitle' }, 'How to use Bone Spatial Studio, what the analysis terms mean, and the project documents.'))),
    h('div', { class: 'help-layout' }, h('aside', { class: 'help-side' }, card({ flush: true, body: nav }), toc), content)));

  function paintNav() {
    mount(nav, TOPICS.map((t) => h('a', {
      class: t.id === current ? 'help-nav-item active' : 'help-nav-item', href: `#/help/${t.id}`,
      'aria-current': t.id === current ? 'page' : null,
      onClick: (e) => {
        if (e.button !== 0 || e.ctrlKey || e.metaKey || e.shiftKey || e.altKey) return; // let the browser open a new tab
        e.preventDefault();
        select(t.id);
      },
    }, icon(t.icon, { size: 16 }), h('span', null, t.label))));
  }

  /** "On this page": the sections of the current topic. */
  function paintToc() {
    let heads = [...content.querySelectorAll('h2[id]')];
    if (heads.length < 2) heads = [...content.querySelectorAll('h2[id], h3[id]')];
    if (heads.length < 2) { mount(toc); return; }
    mount(toc, h('div', { class: 'help-toc-title' }, 'On this page'),
      h('ul', null, heads.map((el) => h('li', { class: el.tagName === 'H3' ? 'sub' : null }, h('a', {
        href: `#/help/${current}`,
        onClick: (e) => {
          e.preventDefault();
          el.scrollIntoView({ behavior: 'smooth', block: 'start' });
        },
      }, el.textContent)))));
  }

  function select(id) {
    if (destroyed || id === current) return;
    setPathSilently(`/help/${id}`);
    show(id);
    window.scrollTo(0, 0);
  }

  async function show(id) {
    current = id;
    paintNav();
    const my = ++token;
    if (ctrl) ctrl.abort();
    ctrl = null;
    const t = TOPICS.find((x) => x.id === id);
    if (!t.doc) {
      mount(content, guide());
      paintToc();
      return;
    }
    let markdown = docs.get(t.doc);
    if (markdown === undefined) {
      mount(content, loading(`Loading ${t.label}…`));
      mount(toc);
      const c = new AbortController();
      ctrl = c;
      try {
        const res = await api(`/api/docs/${t.doc}`, { signal: c.signal });
        markdown = res && typeof res.markdown === 'string' ? res.markdown : '';
        docs.set(t.doc, markdown);
      } catch (err) {
        if (destroyed || my !== token || (err && err.name === 'AbortError')) return;
        mount(content, err.status === 404
          ? card({ body: emptyState({ icon: 'file', title: 'Document not found', text: `The ${t.label} document is not part of this installation.` }) })
          : errorPanel(err, { retry: () => show(id) }));
        return;
      } finally {
        if (ctrl === c) ctrl = null;
      }
      if (destroyed || my !== token) return;
    }
    mount(content, markdown.trim()
      ? h('article', { class: 'prose help-doc' }, renderMarkdown(markdown, { resolveLink, idPrefix: `doc-${t.id}-`, anchorHref: `#/help/${t.id}` }))
      : card({ body: emptyState({ icon: 'file', title: 'Empty document', text: `The ${t.label} document has no content yet.` }) }));
    paintToc();
  }

  // links from one document to another switch the topic in place
  content.addEventListener('click', (e) => {
    const a = e.target instanceof Element ? e.target.closest('a[href^="#/help/"]') : null;
    if (!a || e.defaultPrevented || e.button !== 0 || e.ctrlKey || e.metaKey || e.shiftKey || e.altKey) return;
    const id = a.getAttribute('href').slice('#/help/'.length);
    if (!TOPICS.some((t) => t.id === id)) return;
    e.preventDefault();
    select(id);
  });

  await show(current);

  return () => {
    destroyed = true;
    if (ctrl) ctrl.abort();
  };
}
