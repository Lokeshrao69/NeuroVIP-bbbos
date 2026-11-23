// ---------- shared helpers (window.BX) ----------
(function () {
  const BX = (window.BX = window.BX || {});
  BX.NAMES = { stroke: 'Stroke', meningitis: 'Meningitis', gb: 'Glioblastoma' };
  BX.ORDER = ['stroke', 'meningitis', 'gb'];
  BX.fmtInt = d3.format(',');
  BX.fmtFc = (v) => (v == null || v === '' || isNaN(+v) ? '' : d3.format('+.2f')(+v));
  BX.fmtP = (v) => {
    if (v === '' || v == null) return '';
    const n = +v;
    if (isNaN(n)) return String(v);
    if (n === 0) return '0';
    return n < 0.001 ? n.toExponential(1) : d3.format('.3f')(n);
  };
  BX.num = (v) => (v === '' || v == null ? null : +v);
  BX.css = (name, el) => getComputedStyle(el || document.querySelector('#dash-root') || document.body).getPropertyValue(name).trim();
  BX.upColor = () => BX.css('--color-bad') || '#c0392b';
  BX.downColor = () => (dash.colors && dash.colors[0]) || '#2b6cb0';
  BX.state = { deSearch: '', sharedSearch: '' };

  // A drawing that shows a failed / waiting source
  BX.status = (el, ids) => {
    for (const id of ids) {
      const s = dash.data(id).status;
      if (s === 'error' || s === 'connect' || s === 'declined' || s === 'missing') {
        el.classList.remove('dash-skeleton');
        const f = document.createElement('div');
        f.className = 'bx-fail';
        f.textContent = 'Failed to load';
        el.replaceChildren(f);
        return false;
      }
      if (s !== 'ok') {
        el.replaceChildren();
        el.classList.add('dash-skeleton');
        return false;
      }
    }
    el.classList.remove('dash-skeleton');
    return true;
  };

  // tooltip inside a chart box
  BX.tip = (box) => {
    const t = document.createElement('div');
    t.className = 'bx-tip';
    t.hidden = true;
    box.appendChild(t);
    return {
      show(x, y, head, rows) {
        t.replaceChildren();
        const h = document.createElement('div');
        h.className = 'bx-tip-head';
        h.textContent = head;
        t.appendChild(h);
        for (const r of rows) {
          const d = document.createElement('div');
          d.className = 'bx-tip-row';
          d.textContent = r;
          t.appendChild(d);
        }
        t.hidden = false;
        const bw = box.clientWidth, tw = t.offsetWidth, th = t.offsetHeight;
        let left = x + 12, top = y - th - 8;
        if (left + tw > bw) left = Math.max(0, x - tw - 12);
        if (top < 0) top = y + 12;
        t.style.left = left + 'px';
        t.style.top = top + 'px';
      },
      hide() { t.hidden = true; },
    };
  };

  BX.svg = (el, w, h) =>
    d3.select(el).append('svg').attr('width', '100%').attr('viewBox', `0 0 ${w} ${h}`).attr('preserveAspectRatio', 'xMinYMin meet');

  BX.cleanAxis = (g) => g.attr('font-size', null).attr('font-family', null).call((s) => s.select('.domain').attr('class', 'bx-axis-line'));

  BX.textWidth = (svg, strings) => {
    const t = svg.append('text').attr('visibility', 'hidden');
    let w = 0;
    for (const s of strings) { t.text(s); w = Math.max(w, t.node().getComputedTextLength()); }
    t.remove();
    return w;
  };
})();

// ---------- calculations ----------
dash.calc('deg_summary', {
  title: 'DEG counts per disease',
  description: 'Counts significant genes (padj < 0.05) per disease, split by sign of log2FC.',
  inputs: ['deg'],
  fn: (deg) => ['stroke', 'meningitis', 'gb'].map((disease) => {
    const rows = deg.filter((r) => r.disease === disease);
    return {
      disease,
      n_sig: rows.length,
      n_up: rows.filter((r) => +r.log2fc > 0).length,
      n_down: rows.filter((r) => +r.log2fc < 0).length,
    };
  }),
});

dash.calc('overlap', {
  title: 'DEG overlap between diseases',
  description: 'Counts each significant gene once, in the exact combination of diseases where it is significant.',
  inputs: ['deg'],
  fn: (deg) => {
    const sets = {};
    for (const r of deg) (sets[r.gene_id] = sets[r.gene_id] || new Set()).add(r.disease);
    const key = (s) => ['stroke', 'meningitis', 'gb'].filter((d) => s.has(d)).join('+');
    const names = {
      'stroke': 'stroke_only', 'meningitis': 'meningitis_only', 'gb': 'gb_only',
      'stroke+meningitis': 'stroke_meningitis', 'stroke+gb': 'stroke_gb',
      'meningitis+gb': 'meningitis_gb', 'stroke+meningitis+gb': 'all_three',
    };
    const n = {};
    for (const s of Object.values(sets)) { const k = names[key(s)]; n[k] = (n[k] || 0) + 1; }
    return Object.values(names).map((set) => ({ set, n: n[set] || 0 }));
  },
});

dash.calc('shared_all', {
  title: 'Genes significant in all three',
  description: 'Genes with padj < 0.05 in stroke, meningitis and GB, with each log2FC; sorted by the smallest |log2FC|.',
  inputs: ['deg'],
  fn: (deg) => {
    const by = {};
    for (const r of deg) {
      const g = (by[r.gene_id] = by[r.gene_id] || { gene_id: r.gene_id, symbol: r.symbol });
      g[r.disease] = +r.log2fc;
    }
    return Object.values(by)
      .filter((g) => g.stroke != null && g.meningitis != null && g.gb != null)
      .map((g) => {
        const v = [g.stroke, g.meningitis, g.gb];
        const direction = v.every((x) => x > 0) ? 'up in all' : v.every((x) => x < 0) ? 'down in all' : 'mixed';
        return { ...g, min_abs: Math.min(...v.map(Math.abs)), direction };
      })
      .sort((a, b) => b.min_abs - a.min_abs);
  },
});

dash.calc('shared_direction', {
  title: 'Direction of shared genes',
  description: 'Counts the genes significant in all three diseases by whether they go up everywhere, down everywhere, or differ.',
  inputs: ['shared_all'],
  fn: (shared_all) => ['up in all', 'down in all', 'mixed'].map((direction) => ({
    direction,
    n: shared_all.filter((g) => g.direction === direction).length,
  })),
});

dash.calc('concord', {
  title: 'This pipeline vs classic',
  description: 'Per disease: genes significant in both analyses or only one, and the Pearson correlation of log2FC for genes found by both.',
  inputs: ['compare'],
  fn: (compare) => ['stroke', 'meningitis', 'gb'].map((disease) => {
    const rows = compare.filter((r) => r.disease === disease);
    const both = rows.filter((r) => r.found_by === 'both');
    const x = both.map((r) => +r.ai_log2fc), y = both.map((r) => +r.mk_log2fc);
    const mx = x.reduce((a, b) => a + b, 0) / x.length, my = y.reduce((a, b) => a + b, 0) / y.length;
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < x.length; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
    return {
      disease,
      both: both.length,
      ai_only: rows.filter((r) => r.found_by === 'ai_only').length,
      classic_only: rows.filter((r) => r.found_by === 'classic_only').length,
      r: Math.round((sxy / Math.sqrt(sxx * syy)) * 100000) / 100000,
    };
  }),
});

// ---------- tabs and controls (listen + set attributes only) ----------
(function () {
  const TABS = ['overview', 'de', 'gene', 'cross', 'check', 'classical', 'ml', 'help'];
  const startTab = TABS.includes(location.hash.slice(1)) ? location.hash.slice(1) : null;
  if (startTab) dash.setParams({ tab: startTab });

  const goTab = (t) => { dash.setParams({ tab: t }); dash.setLink(t); window.scrollTo(0, 0); };
  document.querySelectorAll('.bx-tab').forEach((b) =>
    b.addEventListener('click', () => goTab(b.dataset.tab)));
  document.querySelectorAll('[data-goto]').forEach((a) =>
    a.addEventListener('click', (e) => { e.preventDefault(); goTab(a.dataset.goto); }));
  document.querySelectorAll('.nv-logo').forEach((a) =>
    a.addEventListener('click', (e) => { e.preventDefault(); goTab('overview'); }));
  // header gene search -> Gene search page
  document.querySelectorAll('#header-search').forEach((inp) => {
    const go = () => { const v = inp.value.trim(); if (v) { dash.setParams({ gene: v }); goTab('gene'); inp.value = ''; inp.blur(); } };
    inp.addEventListener('keydown', (e) => { if (e.key === 'Enter') go(); });
    document.querySelectorAll('#header-search-btn').forEach((b) => b.addEventListener('click', go));
  });
  document.querySelectorAll('.bx-seg-btn[data-disease]').forEach((b) =>
    b.addEventListener('click', () => dash.setParams({ disease: b.dataset.disease })));
  document.querySelectorAll('.bx-seg-btn[data-lfc]').forEach((b) =>
    b.addEventListener('click', () => dash.setParams({ lfc: +b.dataset.lfc })));
  document.querySelectorAll('.bx-chip[data-gene]').forEach((b) =>
    b.addEventListener('click', () => dash.setParams({ gene: b.dataset.gene })));

  document.querySelectorAll('#gene-input').forEach((inp) => {
    const go = () => { const v = inp.value.trim(); if (v) dash.setParams({ gene: v }); };
    inp.addEventListener('change', go);
    inp.addEventListener('keydown', (e) => { if (e.key === 'Enter') go(); });
  });
  document.querySelectorAll('#de-search').forEach((inp) =>
    inp.addEventListener('input', () => { BX.state.deSearch = inp.value.trim().toLowerCase(); dash.draw('de-body'); }));
  document.querySelectorAll('#shared-search').forEach((inp) =>
    inp.addEventListener('input', () => { BX.state.sharedSearch = inp.value.trim().toLowerCase(); dash.draw('shared-body'); }));

  // click a gene row -> gene lookup
  document.querySelectorAll('.bx-clickable').forEach((t) =>
    t.addEventListener('click', (e) => {
      const tr = e.target.closest('tr[data-gene]');
      if (tr) { dash.setParams({ gene: tr.dataset.gene }); goTab('gene'); }
    }));

  dash.onData(() => {
    const p = dash.params();
    const tab = TABS.includes(p.tab) ? p.tab : 'overview';
    document.querySelectorAll('.bx-page').forEach((s) => { s.hidden = s.dataset.page !== tab; });
    document.querySelectorAll('.bx-tab').forEach((b) => b.setAttribute('aria-selected', String(b.dataset.tab === tab)));
    document.querySelectorAll('.bx-seg-btn[data-disease]').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.disease === (p.disease || 'stroke'))));
    document.querySelectorAll('.bx-seg-btn[data-lfc]').forEach((b) => b.setAttribute('aria-pressed', String(+b.dataset.lfc === +(p.lfc ?? 1))));
    document.querySelectorAll('#gene-input').forEach((inp) => { if (inp !== document.activeElement) inp.value = p.gene || ''; });
  });
})();
