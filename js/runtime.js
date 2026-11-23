// Minimal stand-alone runtime that provides the `dash` API used by the page scripts.
// Loads CSV files from ./data, runs calculations, and redraws on data / filter / resize changes.
(function () {
  const DATASETS = {
    deg: 'data/deg_significant.csv',
    compare: 'data/ai_vs_classic.csv',
    samples: 'data/samples_qc.csv',
    studies: 'data/studies.csv',
    ml_metrics: 'data/ml_metrics.csv',
    ml_features: 'data/ml_features.csv',
    // classical analysis (meningitis), written by analysis/classical_meningitis.py
    cl_methods: 'data/classical_methods.csv',
    cl_summary: 'data/classical_summary.csv',
    cl_deg: 'data/classical_deg.csv',
    cl_genes: 'data/classical_genes.csv',
    cl_counts: 'data/classical_counts.csv',
    cl_top: 'data/classical_top.csv',
    cl_paper_top: 'data/classical_paper_top.csv',
    cl_pca: 'data/classical_pca.csv',
    cl_sample_dist: 'data/classical_sample_dist.csv',
    cl_dendrogram: 'data/classical_dendrogram.csv',
    cl_sample_heatmap: 'data/classical_sample_heatmap.csv',
    cl_pathways: 'data/classical_pathways.csv',
    cl_pathway_genes: 'data/classical_pathway_genes.csv',
  };
  const DEFAULTS = { tab: 'overview', disease: 'stroke', lfc: 1, gene: 'Selp', clsrc: 'pa', clrank: 'padj', clhl: 'paper', clheat: 'z', clpca: 'top500', cldb: 'KEGG', clpaper: 'all', clpanel: '' };

  const store = {};   // id -> {status, data, columns}
  const calcs = {};   // id -> {inputs, fn, cache}
  const draws = {};   // name -> fn
  const listeners = [];
  let params = { ...DEFAULTS };
  let scheduled = false;

  for (const [id, url] of Object.entries(DATASETS)) {
    store[id] = { status: 'loading', data: [], columns: [] };
    d3.csv(url).then((rows) => {
      store[id] = { status: 'ok', data: rows, columns: rows.columns || [] };
      for (const c of Object.values(calcs)) c.cache = undefined;
      schedule();
    }).catch((e) => {
      store[id] = { status: 'error', data: [], columns: [], message: String(e) };
      schedule();
    });
  }

  function data(id) {
    if (store[id]) return store[id];
    const c = calcs[id];
    if (!c) return { status: 'missing', data: [] };
    if (c.cache) return c.cache;
    const inputs = c.inputs.map(data);
    const bad = inputs.find((d) => d.status !== 'ok');
    if (bad) return { status: bad.status, data: [] };
    try {
      c.cache = { status: 'ok', data: c.fn(...inputs.map((d) => d.data)) };
    } catch (e) {
      console.error('calc ' + id, e);
      c.cache = { status: 'error', data: [], message: String(e) };
    }
    return c.cache;
  }

  function isShown(el) { return !el.closest('[hidden]'); }

  function runDraw(name) {
    const fn = draws[name];
    if (!fn) return;
    document.querySelectorAll(`[data-draw="${name}"]`).forEach((el) => {
      if (!isShown(el)) return;
      try { fn(el); } catch (e) { console.error('draw ' + name, e); }
    });
  }

  function render() {
    scheduled = false;
    for (const fn of listeners) { try { fn(); } catch (e) { console.error(e); } }
    for (const name of Object.keys(draws)) runDraw(name);
  }

  function schedule() {
    if (scheduled) return;
    scheduled = true;
    requestAnimationFrame(render);
  }

  const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

  window.dash = {
    data,
    calc(id, { inputs, fn }) { calcs[id] = { inputs, fn, cache: undefined }; },
    draw(name, fn) { if (fn) { draws[name] = fn; schedule(); } else runDraw(name); },
    onData(fn) { listeners.push(fn); schedule(); },
    params() { return { ...params }; },
    setParams(p) { params = { ...params, ...p }; schedule(); },
    resetParams() { params = { ...DEFAULTS }; schedule(); },
    setLink(name) { history.replaceState(null, '', '#' + name); },
    refresh() {},
    get colors() { return [1, 2, 3, 4, 5, 6, 7, 8].map((i) => css('--series-' + i)); },
  };

  let w = window.innerWidth;
  window.addEventListener('resize', () => { if (window.innerWidth !== w) { w = window.innerWidth; schedule(); } });
})();
