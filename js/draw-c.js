// ---------- Classical analysis (meningitis, Devraj et al. 2020) ----------
// Reads the classical_*.csv files written by analysis/classical_meningitis.py.
(function () {
  const BX = window.BX;
  const LFC = Math.log2(1.5);   // paper threshold, |FC| >= 1.5
  const PADJ = 0.05;
  const SRC = { pa: 'Paper counts (re-run)', nv: 'NeuroVIP pipeline', cl: 'Classic analysis' };
  const SRC_LONG = { pa: 'Paper count matrix (GEO), DESeq2 re-run here', nv: 'NeuroVIP pipeline', cl: 'Classic analysis' };
  const src = () => { const s = dash.params().clsrc; return SRC[s] ? s : 'pa'; };
  const P = (k, def) => dash.params()[k] || def;
  const empty = (el, msg) => { const e = document.createElement('div'); e.className = 'bx-empty'; e.textContent = msg; el.appendChild(e); };
  const td = (tr, v, cls) => { const c = document.createElement('td'); c.textContent = v; if (cls) c.className = cls; tr.appendChild(c); return c; };
  const pill = (txt, kind) => { const s = document.createElement('span'); s.className = 'bx-pill bx-pill-' + kind; s.textContent = txt; return s; };
  const muted = () => BX.css('--cds-chart-muted') || '#bbb';
  const callColor = (c) => (c === 'up' ? BX.upColor() : c === 'down' ? BX.downColor() : muted());
  const fcCls = (v) => (v === '' || v == null ? 'bx-num' : 'bx-num ' + (+v > 0 ? 'bx-up' : 'bx-down'));
  const ok = (id) => dash.data(id).status === 'ok';
  const unavailable = (el, msg) => {
    const box = document.createElement('div');
    box.className = 'bx-unavail';
    box.textContent = msg;
    el.appendChild(box);
  };
  // gene_id -> deg row, built once per loaded table
  let degIndex = null, degIndexSrc = null;
  const degById = () => {
    const d = dash.data('cl_deg').data;
    if (degIndexSrc !== d) { degIndex = new Map(d.map((r) => [r.gene_id, r])); degIndexSrc = d; }
    return degIndex;
  };
  const sampleOrder = (rows, key = 'sample') => Array.from(new Set(rows.map((r) => r[key])))
    .sort((a, b) => (a.startsWith('Healthy') === b.startsWith('Healthy') ? a.localeCompare(b) : a.startsWith('Healthy') ? 1 : -1));
  const sampleShort = (s) => (s.startsWith('Healthy') ? 'H' : 'M') + s.split('_')[1];

  // ----- controls (classical page only) -----
  const seg = (attr, key) => document.querySelectorAll(`.bx-seg-btn[data-${attr}]`).forEach((b) => b.addEventListener('click', () => dash.setParams({ [key]: b.dataset[attr] })));
  seg('clsrc', 'clsrc'); seg('clrank', 'clrank'); seg('clhl', 'clhl'); seg('clheat', 'clheat'); seg('clpca', 'clpca'); seg('cldb', 'cldb'); seg('clpaper', 'clpaper');
  document.querySelectorAll('#cl-panel-select').forEach((s) => s.addEventListener('change', () => dash.setParams({ clpanel: s.value })));
  dash.onData(() => {
    const cur = { clsrc: src(), clrank: P('clrank', 'padj'), clhl: P('clhl', 'paper'), clheat: P('clheat', 'z'), clpca: P('clpca', 'top500'), cldb: P('cldb', 'KEGG'), clpaper: P('clpaper', 'all') };
    for (const [k, v] of Object.entries(cur)) document.querySelectorAll(`.bx-seg-btn[data-${k}]`).forEach((b) => b.setAttribute('aria-pressed', String(b.dataset[k] === v)));
  });

  // ----- tiles -----
  dash.draw('cl-tile', (el) => {
    const d = dash.data('cl_summary');
    if (d.status !== 'ok') { el.textContent = '—'; return; }
    const tkey = el.dataset.tkey; // strict | padj | reported
    const s = tkey === 'reported' ? 'paper' : src();
    const row = d.data.find((r) => r.source === s && r.tkey === tkey);
    el.textContent = row ? BX.fmtInt(+row[el.dataset.field]) : '—';
  });
  dash.draw('cl-tile-sub', (el) => {
    const d = dash.data('cl_summary');
    if (d.status !== 'ok') return;
    const tkey = el.dataset.tkey;
    const row = d.data.find((r) => r.source === (tkey === 'reported' ? 'paper' : src()) && r.tkey === tkey);
    if (row) el.textContent = `up / down · ${row.threshold} · of ${BX.fmtInt(+row.genes_tested)} genes`;
  });
  dash.draw('cl-src-name', (el) => { el.textContent = SRC[src()]; });

  // ----- methods table -----
  dash.draw('cl-methods-body', (tb) => {
    tb.replaceChildren();
    const d = dash.data('cl_methods');
    if (d.status !== 'ok') return;
    const kind = { applied: 'ok', compared: 'ok', reused: 'info', 'partly applied': 'warn', adapted: 'warn', 'not possible': 'bad' };
    for (const r of d.data) {
      const tr = document.createElement('tr');
      tr.dataset.key = r.step;
      td(tr, r.step);
      td(tr, r.paper_method, 'bx-muted');
      td(tr, '').appendChild(pill(r.status, kind[r.status] || 'info'));
      td(tr, r.detail);
      tb.appendChild(tr);
    }
  });

  // ----- DEG counts vs paper: diverging bars -----
  dash.draw('cl-counts', (el) => {
    if (!BX.status(el, ['cl_summary'])) return;
    el.replaceChildren();
    const rows = dash.data('cl_summary').data.map((r) => ({
      label: r.source === 'paper' ? r.analysis : `${SRC[r.source]} · ${r.threshold}`,
      up: +r.n_up, down: +r.n_down, paper: r.source === 'paper', sel: r.source === src(), pu: r.pct_up, pd: r.pct_down, tested: +r.genes_tested,
    }));
    const W = Math.max(320, el.clientWidth), rowH = 30;
    const narrow = W < 640;
    const top = narrow ? 14 : 0;
    const H = 24 + rows.length * (rowH + top) + 4;
    const svg = BX.svg(el, W, H);
    const lw = narrow ? 0 : BX.textWidth(svg, rows.map((r) => r.label)) + 12;
    const mx = d3.max(rows, (r) => Math.max(r.up, r.down));
    const half = (W - lw) / 2;
    const x0 = lw + half;
    const x = d3.scaleLinear([0, mx], [0, half - 46]);
    svg.append('text').attr('x', x0 - 6).attr('y', 12).attr('text-anchor', 'end').text('◀ down-regulated');
    svg.append('text').attr('x', x0 + 6).attr('y', 12).text('up-regulated ▶');
    const tip = BX.tip(el);
    rows.forEach((r, i) => {
      const yy = 24 + i * (rowH + top);
      const fw = r.paper || r.sel ? 600 : null;
      if (narrow) svg.append('text').attr('x', 0).attr('y', yy + 2).attr('font-size', 10.5).attr('font-weight', fw).text(r.label);
      else svg.append('text').attr('class', r.paper || r.sel ? 'bx-lab' : null).attr('x', 0).attr('y', yy + 14).attr('font-weight', fw).text(r.label);
      const by = yy + top + 2;
      const show = (ev) => { const [px, py] = d3.pointer(ev, el); tip.show(px, py, r.label, [`${BX.fmtInt(r.up)} up (${r.pu} %)`, `${BX.fmtInt(r.down)} down (${r.pd} %)`, `of ${BX.fmtInt(r.tested)} genes`]); };
      const op = r.paper || r.sel ? 1 : 0.55;
      svg.append('rect').attr('x', x0 - x(r.down)).attr('y', by).attr('width', x(r.down)).attr('height', 18).attr('rx', 2)
        .attr('fill', BX.downColor()).attr('fill-opacity', op).on('pointermove', show).on('pointerleave', () => tip.hide());
      svg.append('rect').attr('x', x0).attr('y', by).attr('width', x(r.up)).attr('height', 18).attr('rx', 2)
        .attr('fill', BX.upColor()).attr('fill-opacity', op).on('pointermove', show).on('pointerleave', () => tip.hide());
      svg.append('text').attr('class', 'bx-lab').attr('x', x0 - x(r.down) - 4).attr('y', by + 13).attr('text-anchor', 'end').text(BX.fmtInt(r.down));
      svg.append('text').attr('class', 'bx-lab').attr('x', x0 + x(r.up) + 4).attr('y', by + 13).text(BX.fmtInt(r.up));
    });
    svg.append('line').attr('class', 'bx-axis-line').attr('x1', x0).attr('x2', x0).attr('y1', 18).attr('y2', H - 4);
  });
  // generated from the numbers: which threshold matches the paper's reported counts best
  dash.draw('cl-counts-note', (el) => {
    const d = dash.data('cl_summary');
    if (d.status !== 'ok') { el.textContent = ''; return; }
    const paper = d.data.find((r) => r.tkey === 'reported');
    const runs = d.data.filter((r) => r.source !== 'paper');
    const err = (r) => (Math.abs(+r.n_up - +paper.n_up) + Math.abs(+r.n_down - +paper.n_down)) / (+paper.n_up + +paper.n_down);
    const best = runs.slice().sort((a, b) => err(a) - err(b))[0];
    const byT = (t) => d3.mean(runs.filter((r) => r.tkey === t), err);
    const pct = d3.format('.1%');
    el.textContent = `Closest to the paper's ${BX.fmtInt(+paper.n_up)} up / ${BX.fmtInt(+paper.n_down)} down: ${SRC[best.source]} at “${best.threshold}” `
      + `(${BX.fmtInt(+best.n_up)} / ${BX.fmtInt(+best.n_down)}, ${pct(err(best))} off). Averaged over the three runs, padj-only counts are ${pct(byT('padj'))} off `
      + `and counts with the stated 1.5-fold cut-off are ${pct(byT('strict'))} off${byT('padj') < byT('strict') ? ', so the reported numbers most likely did not apply the fold-change cut-off' : ''}. `
      + 'The NeuroVIP and classic runs use a different genome build (GRCm39 vs mm10), aligner version and annotation, so small differences are expected.';
  });

  // ----- volcano with paper thresholds -----
  const highlightSet = () => {
    const mode = P('clhl', 'paper');
    if (mode === 'top' && ok('cl_paper_top')) return new Set(dash.data('cl_paper_top').data.map((r) => r.symbol));
    if (mode === 'panel' && ok('cl_pathway_genes')) {
      const p = dash.params().clpanel;
      const ids = new Set(dash.data('cl_pathway_genes').data.filter((r) => r.pathway_id === p).map((r) => r.gene_id));
      const idx = degById();
      return new Set([...ids].map((g) => idx.get(g)?.symbol).filter(Boolean));
    }
    return new Set(ok('cl_genes') ? dash.data('cl_genes').data.map((r) => r.symbol) : []);
  };
  dash.draw('cl-volcano', (el) => {
    if (!BX.status(el, ['cl_deg', 'cl_genes', 'cl_pathway_genes', 'cl_paper_top'])) return;
    el.replaceChildren();
    const s = src();
    const rows = dash.data('cl_deg').data.filter((r) => r[s + '_padj'] !== '').map((r) => {
      const p = +r[s + '_padj'];
      return { symbol: r.symbol, fc: +r[s + '_log2fc'], padj: p, y: -Math.log10(Math.max(p, 1e-300)), call: r[s + '_call'] };
    });
    if (!rows.length) { empty(el, 'No genes'); return; }
    const hlSet = highlightSet();
    const W = Math.max(320, el.clientWidth), H = Math.max(380, Math.min(540, W * 0.55));
    const m = { t: 12, r: 16, b: 40, l: 48 };
    const svg = BX.svg(el, W, H);
    const xMax = d3.max(rows, (d) => Math.abs(d.fc)) * 1.05;
    const x = d3.scaleLinear([-xMax, xMax], [m.l, W - m.r]);
    const y = d3.scaleLinear([0, d3.max(rows, (d) => d.y) * 1.06], [H - m.b, m.t]);
    svg.append('g').selectAll('line').data(y.ticks(6)).join('line').attr('class', 'bx-grid-line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y).attr('y2', y);
    svg.append('g').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(Math.max(3, Math.floor(W / 80))).tickSizeOuter(0)).call(BX.cleanAxis);
    svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6).tickSizeOuter(0)).call(BX.cleanAxis);
    svg.append('text').attr('x', (m.l + W - m.r) / 2).attr('y', H - 6).attr('text-anchor', 'middle').text('log2 fold change (infected vs healthy)');
    svg.append('text').attr('transform', `translate(12,${(m.t + H - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text('−log10 padj');
    for (const v of [-LFC, LFC]) svg.append('line').attr('class', 'bx-ref').attr('x1', x(v)).attr('x2', x(v)).attr('y1', m.t).attr('y2', H - m.b);
    const yp = y(-Math.log10(PADJ));
    if (yp > m.t && yp < H - m.b) svg.append('line').attr('class', 'bx-ref').attr('x1', m.l).attr('x2', W - m.r).attr('y1', yp).attr('y2', yp);
    const base = rows.filter((d) => !hlSet.has(d.symbol));
    const hi = rows.filter((d) => hlSet.has(d.symbol));
    svg.append('g').selectAll('circle').data(base).join('circle')
      .attr('cx', (d) => x(d.fc)).attr('cy', (d) => y(d.y)).attr('r', 2)
      .attr('fill', (d) => callColor(d.call)).attr('fill-opacity', hi.length ? 0.25 : 0.55);
    const fg = BX.css('--color-fg') || '#222';
    svg.append('g').selectAll('circle').data(hi).join('circle')
      .attr('cx', (d) => x(d.fc)).attr('cy', (d) => y(d.y)).attr('r', 4.2)
      .attr('fill', (d) => callColor(d.call)).attr('stroke', fg).attr('stroke-width', 1);
    const placed = [];
    for (const d of hi.slice().sort((a, b) => b.y - a.y)) {
      const px = x(d.fc), py = y(d.y) - 7;
      const t = svg.append('text').attr('class', 'bx-lab').attr('x', px).attr('y', py).attr('font-weight', 600)
        .attr('text-anchor', px > W - 80 ? 'end' : px < m.l + 60 ? 'start' : 'middle').text(d.symbol);
      const b = t.node().getBBox();
      if (placed.some((p) => !(b.x + b.width < p.x || p.x + p.width < b.x || b.y + b.height < p.y || p.y + p.height < b.y))) t.remove();
      else placed.push(b);
    }
    const ring = svg.append('circle').attr('r', 6).attr('fill', 'none').attr('stroke', fg).attr('visibility', 'hidden');
    const tip = BX.tip(el);
    const dl = d3.Delaunay.from(rows, (d) => x(d.fc), (d) => y(d.y));
    const lbl = { up: 'up (passes paper thresholds)', down: 'down (passes paper thresholds)', below_fc: 'padj ≤ 0.05 but |FC| < 1.5', '': 'not significant' };
    svg.on('pointermove', (ev) => {
      const [px, py] = d3.pointer(ev);
      const d = rows[dl.find(px, py)];
      if (!d || Math.hypot(x(d.fc) - px, y(d.y) - py) > 24) { ring.attr('visibility', 'hidden'); tip.hide(); return; }
      ring.attr('cx', x(d.fc)).attr('cy', y(d.y)).attr('visibility', 'visible');
      const k = el.clientWidth / W;
      tip.show(x(d.fc) * k, y(d.y) * k, d.symbol, ['log2FC ' + BX.fmtFc(d.fc) + '  (FC ' + d3.format('.2f')(2 ** d.fc) + ')', 'padj ' + BX.fmtP(d.padj), lbl[d.call] || d.call]);
    }).on('pointerleave', () => { ring.attr('visibility', 'hidden'); tip.hide(); });
  });
  dash.draw('cl-volcano-note', (el) => {
    const s = src();
    const n = ok('cl_deg') ? dash.data('cl_deg').data.filter((r) => r[s + '_padj'] !== '').length : 0;
    const what = s === 'pa'
      ? `All ${BX.fmtInt(n)} genes tested in the re-run are shown, including non-significant ones.`
      : `This run stores only significant genes (${BX.fmtInt(n)}), so the cloud of unchanged genes near zero is absent; switch to “${SRC.pa}” for the full plot.`;
    const mode = P('clhl', 'paper');
    const hl = mode === 'panel' ? 'genes of the selected pathway (see below)' : mode === 'top' ? 'the paper’s Fig. 6c top-25 genes' : 'genes named in the paper';
    el.textContent = `${what} Highlighted: ${hl}.`;
  });

  // ----- genes named in the paper -----
  dash.draw('cl-genes-chart', (el) => {
    if (!BX.status(el, ['cl_genes'])) return;
    el.replaceChildren();
    const s = src();
    const rows = dash.data('cl_genes').data.map((r) => ({ ...r, v: BX.num(r[s + '_log2fc']), p: r[s + '_padj'], obs: r[s + '_observed'], ag: r[s + '_agreement'] }));
    const W = Math.max(320, el.clientWidth), rowH = 22;
    const cats = Array.from(new Set(rows.map((r) => r.category)));
    const H = rows.length * rowH + cats.length * 18 + 40;
    const svg = BX.svg(el, W, H);
    const lw = BX.textWidth(svg, rows.map((r) => r.paper_label)) + 12;
    const rw = Math.min(150, W * 0.28);
    const ext = Math.max(2, d3.max(rows, (r) => Math.abs(r.v || 0))) * 1.08;
    const x = d3.scaleLinear([-ext, ext], [lw, W - rw]).nice();
    svg.append('g').attr('transform', `translate(0,${H - 22})`).call(d3.axisBottom(x).ticks(Math.max(3, Math.floor((W - lw - rw) / 70))).tickSizeOuter(0)).call(BX.cleanAxis);
    svg.append('text').attr('x', (lw + W - rw) / 2).attr('y', H - 0).attr('text-anchor', 'middle').attr('font-size', 10).text('log2 fold change');
    svg.append('line').attr('class', 'bx-axis-line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 4).attr('y2', H - 22);
    for (const v of [-LFC, LFC]) svg.append('line').attr('class', 'bx-ref').attr('x1', x(v)).attr('x2', x(v)).attr('y1', 4).attr('y2', H - 22);
    svg.append('text').attr('x', W - rw + 10).attr('y', 12).attr('font-weight', 600).text('Paper → here');
    const arrow = { up: '▲', down: '▼', unchanged: '≈', housekeeping: 'ref', marker: 'marker' };
    const okc = BX.css('--color-ok') || '#2e7d32', bad = BX.css('--color-warn') || '#b25e00';
    const tip = BX.tip(el);
    let yy = 6;
    for (const c of cats) {
      yy += 16;
      svg.append('text').attr('x', 0).attr('y', yy).attr('font-size', 10.5).attr('font-weight', 600).text(c);
      yy += 2;
      for (const r of rows.filter((q) => q.category === c)) {
        const cy = yy + rowH / 2;
        svg.append('text').attr('class', 'bx-lab').attr('x', lw - 8).attr('y', cy + 4).attr('text-anchor', 'end').text(r.paper_label);
        const sig = r.p !== '' && +r.p <= PADJ;
        if (r.v == null || !sig) {
          svg.append('text').attr('x', x(0) + 6).attr('y', cy + 4).attr('font-size', 10).text(r.v == null ? 'not significant (not in results table)' : `n.s. (log2FC ${BX.fmtFc(r.v)})`);
        } else {
          svg.append('rect').attr('x', Math.min(x(0), x(r.v))).attr('y', cy - 6).attr('width', Math.max(1, Math.abs(x(r.v) - x(0)))).attr('height', 12).attr('rx', 2)
            .attr('fill', r.v > 0 ? BX.upColor() : BX.downColor()).attr('fill-opacity', Math.abs(r.v) >= LFC ? 0.9 : 0.45)
            .on('pointermove', (ev) => { const [px, py] = d3.pointer(ev, el); tip.show(px, py, `${r.paper_label} (${r.symbol})`, [SRC_LONG[s], 'log2FC ' + BX.fmtFc(r.v) + '  (FC ' + d3.format('.2f')(2 ** r.v) + ')', 'padj ' + BX.fmtP(r.p), 'Paper: ' + r.paper_claim + ' · ' + r.paper_source]); })
            .on('pointerleave', () => tip.hide());
          svg.append('text').attr('x', r.v > 0 ? x(r.v) + 4 : x(r.v) - 4).attr('y', cy + 4).attr('text-anchor', r.v > 0 ? 'start' : 'end').attr('font-size', 10).text(BX.fmtFc(r.v));
        }
        svg.append('text').attr('x', W - rw + 10).attr('y', cy + 4).text(`${arrow[r.paper_claim] || r.paper_claim} → ${r.obs}`);
        if (r.ag && r.ag !== '—') svg.append('text').attr('x', W - 4).attr('y', cy + 4).attr('text-anchor', 'end').attr('font-weight', 600)
          .style('fill', r.ag === 'agrees' || r.ag === 'stable' ? okc : bad).text(r.ag === 'agrees' || r.ag === 'stable' ? '✓' : '!');
        yy += rowH;
      }
    }
  });

  dash.draw('cl-genes-body', (tb) => {
    tb.replaceChildren();
    const d = dash.data('cl_genes');
    if (d.status !== 'ok') return;
    const s = src();
    const f = (v) => (v === '' ? '—' : BX.fmtFc(v));
    for (const r of d.data) {
      const tr = document.createElement('tr');
      tr.dataset.key = r.symbol;
      tr.dataset.gene = r.symbol;
      td(tr, r.symbol);
      td(tr, r.paper_label);
      td(tr, r.paper_claim + ' · ' + r.paper_source, 'bx-muted');
      for (const k of ['pa', 'nv', 'cl']) {
        td(tr, f(r[k + '_log2fc']), fcCls(r[k + '_log2fc']) + (k === s ? ' bx-strong' : ''));
        td(tr, r[k + '_padj'] === '' ? '—' : BX.fmtP(r[k + '_padj']), 'bx-num' + (k === s ? ' bx-strong' : ''));
      }
      const ag = r[s + '_agreement'];
      const kind = { agrees: 'ok', stable: 'ok', differs: 'bad', changed: 'warn' }[ag];
      td(tr, '').appendChild(kind ? pill(ag, kind) : document.createTextNode('—'));
      tb.appendChild(tr);
    }
  });

  // ----- read counts per sample (paper Fig. 6d) -----
  dash.draw('cl-reads', (el) => {
    if (!BX.status(el, ['cl_counts', 'cl_genes'])) return;
    el.replaceChildren();
    const data = dash.data('cl_counts').data;
    if (!data.length) { unavailable(el, 'No per-sample counts in this build.'); return; }
    const s = src();
    const genes = dash.data('cl_genes').data;
    const byGene = d3.group(data, (r) => r.symbol);
    const order = genes.map((g) => g.symbol).filter((g) => byGene.has(g));
    const W = Math.max(300, el.clientWidth);
    const cols = Math.max(2, Math.min(6, Math.floor(W / 150)));
    const cw = W / cols, ch = 164;
    const H = Math.ceil(order.length / cols) * ch + 6;
    const svg = BX.svg(el, W, H);
    const tip = BX.tip(el);
    order.forEach((sym, i) => {
      const g = genes.find((q) => q.symbol === sym);
      const rows = byGene.get(sym);
      const gx = (i % cols) * cw, gy = Math.floor(i / cols) * ch;
      const box = svg.append('g').attr('transform', `translate(${gx},${gy})`);
      const m = { t: 38, b: 22, l: 42, r: 8 };
      const y = d3.scaleLinear([0, d3.max(rows, (r) => +r.count) * 1.1 || 1], [ch - m.b, m.t]).nice();
      const x = d3.scalePoint(['control', 'case'], [m.l + 18, cw - m.r - 18]);
      box.append('text').attr('class', 'bx-lab').attr('x', m.l).attr('y', 14).attr('font-weight', 600).text(g.paper_label);
      const p = g[s + '_padj'];
      const sigTxt = p === '' ? 'not in results' : +p <= PADJ ? `padj ${BX.fmtP(p)} · ${BX.fmtFc(g[s + '_log2fc'])}` : `n.s. (padj ${BX.fmtP(p)})`;
      box.append('text').attr('x', m.l).attr('y', 26).attr('font-size', 9.5).text(sigTxt);
      box.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(3).tickFormat(d3.format('~s')).tickSizeOuter(0)).call(BX.cleanAxis);
      box.append('g').attr('transform', `translate(0,${ch - m.b})`).call(d3.axisBottom(x).tickFormat((c) => (c === 'case' ? 'Infected' : 'Healthy')).tickSizeOuter(0)).call(BX.cleanAxis);
      for (const cond of ['control', 'case']) {
        const v = rows.filter((r) => r.condition === cond);
        const mean = d3.mean(v, (r) => +r.count);
        const col = cond === 'case' ? BX.upColor() : BX.downColor();
        box.append('rect').attr('x', x(cond) - 14).attr('width', 28).attr('y', y(mean)).attr('height', y(0) - y(mean)).attr('fill', col).attr('fill-opacity', 0.25).attr('rx', 2);
        v.forEach((r, j) => box.append('circle').attr('cx', x(cond) - 9 + j * 6).attr('cy', y(+r.count)).attr('r', 3).attr('fill', col)
          .on('pointermove', (ev) => { const [px, py] = d3.pointer(ev, el); tip.show(px, py, `${g.paper_label} · ${r.label}`, [`${BX.fmtInt(+r.count)} normalised reads`, r.run]); })
          .on('pointerleave', () => tip.hide()));
      }
    });
  });

  // ----- top 25 up / down heatmap -----
  dash.draw('cl-top-heat', (el) => {
    if (!BX.status(el, ['cl_top', 'cl_sample_heatmap', 'cl_deg'])) return;
    el.replaceChildren();
    const s = src();
    const rankBy = P('clrank', 'padj');
    const perSample = P('clheat', 'z') === 'z';
    const sampleData = dash.data('cl_sample_heatmap').data;
    const rows = dash.data('cl_top').data.filter((r) => r.source === s && r.rank_by === rankBy)
      .sort((a, b) => (a.direction === b.direction ? +a.rank - +b.rank : a.direction === 'up' ? -1 : 1));
    if (!rows.length) { empty(el, 'No genes pass the paper thresholds.'); return; }
    const zmap = new Map(sampleData.map((r) => [r.gene_id + '|' + r.sample, +r.z]));
    const idx = degById();
    const cols = perSample ? sampleOrder(sampleData) : ['pa', 'nv', 'cl'];
    const colName = (c) => (perSample ? c.replace('_', ' ') : SRC[c]);
    const val = perSample ? (r, c) => zmap.get(r.gene_id + '|' + c) : (r, c) => BX.num(idx.get(r.gene_id)?.[c + '_log2fc'] ?? '');
    const W = Math.max(300, el.clientWidth);
    const cellH = 15;
    const svg = BX.svg(el, W, 10);
    const gw = BX.textWidth(svg, rows.map((r) => r.symbol + ' ●')) + 12;
    const cw = Math.min(perSample ? 60 : 150, (W - gw - 8) / cols.length);
    const vals = rows.flatMap((r) => cols.map((c) => val(r, c))).filter((v) => v != null);
    const lim = perSample ? 2 : Math.min(10, d3.max(vals, Math.abs) || 1);
    const basec = BX.css('--color-panel') || '#fff';
    const pos = d3.interpolateRgb(basec, BX.upColor()), neg = d3.interpolateRgb(basec, BX.downColor());
    const color = (v) => (v == null ? (BX.css('--cds-border') || '#eee') : v >= 0 ? pos(Math.min(1, v / lim)) : neg(Math.min(1, -v / lim)));
    cols.forEach((c, j) => svg.append('text').attr('x', gw + j * cw + cw / 2).attr('y', 14).attr('text-anchor', 'middle')
      .attr('font-size', perSample ? 9.5 : 11).text(perSample ? sampleShort(c) : colName(c)));
    if (perSample) cols.forEach((c, j) => svg.append('rect').attr('x', gw + j * cw + 1).attr('y', 20).attr('width', cw - 2).attr('height', 6)
      .attr('fill', c.startsWith('Healthy') ? muted() : BX.upColor()));
    const tip = BX.tip(el);
    let yy = 32;
    let lastDir = null;
    rows.forEach((r) => {
      if (r.direction !== lastDir) {
        const n = rows.filter((q) => q.direction === r.direction);
        const inPaper = n.filter((q) => q.in_paper_top25 === 'True').length;
        svg.append('text').attr('x', 0).attr('y', yy + 11).attr('font-size', 10.5).attr('font-weight', 600)
          .text(`Top ${n.length} ${r.direction} · ${inPaper} in paper Fig. 6c`);
        yy += 16;
        lastDir = r.direction;
      }
      const inPaper = r.in_paper_top25 === 'True';
      svg.append('text').attr('class', 'bx-lab').attr('x', gw - 6).attr('y', yy + cellH - 4).attr('text-anchor', 'end').attr('font-size', 10.5)
        .attr('font-weight', inPaper ? 600 : null).text(r.symbol + (inPaper ? ' ●' : ''));
      cols.forEach((c, j) => {
        const v = val(r, c);
        svg.append('rect').attr('x', gw + j * cw + 1).attr('y', yy + 1).attr('width', cw - 2).attr('height', cellH - 2).attr('rx', 2).attr('fill', color(v))
          .on('pointermove', (ev) => {
            const [px, py] = d3.pointer(ev, el);
            const head = `${r.symbol} · ${colName(c)}`;
            tip.show(px, py, head, perSample
              ? [v == null ? 'not in the paper count matrix' : 'Z-score ' + d3.format('+.2f')(v), `${SRC[s]}: log2FC ${BX.fmtFc(r.log2fc)}, padj ${BX.fmtP(r.padj)}`]
              : ['log2FC ' + (v == null ? 'n.s. / absent' : BX.fmtFc(v)), 'padj ' + BX.fmtP(idx.get(r.gene_id)?.[c + '_padj'] ?? '')]);
          })
          .on('pointerleave', () => tip.hide());
        if (!perSample && cw > 46 && v != null) svg.append('text').attr('x', gw + j * cw + cw / 2).attr('y', yy + cellH - 4).attr('text-anchor', 'middle').attr('font-size', 9.5).attr('pointer-events', 'none')
          .style('fill', Math.abs(v) / lim > 0.6 ? '#fff' : null).text(BX.fmtFc(v));
      });
      yy += cellH;
    });
    svg.attr('viewBox', `0 0 ${W} ${yy + 6}`);
  });
  dash.draw('cl-top-note', (el) => {
    const s = src();
    const per = P('clheat', 'z') === 'z';
    el.textContent = (per
      ? 'Row Z-scores of VST expression across the 8 samples (as in paper Fig. 6c); M = infected, H = healthy; bar above columns: red = infected, grey = healthy.'
      : 'log2 fold change of the same genes in each DESeq2 run; empty = not significant or absent in that run.')
      + ` Genes are the top 25 of “${SRC[s]}” passing padj ≤ 0.05 and |FC| ≥ 1.5; unannotated Ensembl IDs are skipped. ● = also in the paper's Fig. 6c list.`;
  });

  // ----- paper Fig. 6c vs here -----
  dash.draw('cl-ptop-body', (tb) => {
    tb.replaceChildren();
    const d = dash.data('cl_paper_top');
    if (d.status !== 'ok') return;
    const s = src();
    for (const r of d.data) {
      const tr = document.createElement('tr');
      tr.dataset.key = r.direction + ',' + r.paper_label;
      tr.dataset.gene = r.symbol;
      td(tr, r.direction, r.direction === 'up' ? 'bx-up' : 'bx-down');
      td(tr, r.paper_label === r.symbol ? r.symbol : `${r.paper_label} → ${r.symbol}`);
      for (const k of ['pa', 'nv', 'cl']) {
        const rk = r[k + '_rank'];
        td(tr, rk === '' ? (r[k + '_call'] === 'absent' ? 'absent' : 'not ' + r.direction) : '#' + rk, 'bx-num' + (k === s ? ' bx-strong' : '') + (rk !== '' && +rk <= 25 ? '' : ' bx-muted'));
      }
      td(tr, BX.fmtFc(r[s + '_log2fc']), fcCls(r[s + '_log2fc']));
      td(tr, r[s + '_padj'] === '' ? '—' : BX.fmtP(r[s + '_padj']), 'bx-num');
      tb.appendChild(tr);
    }
  });
  dash.draw('cl-ptop-note', (el) => {
    const d = dash.data('cl_paper_top');
    if (d.status !== 'ok') return;
    const s = src();
    const cnt = (dir) => d.data.filter((r) => r.direction === dir && r[s + '_rank'] !== '' && +r[s + '_rank'] <= 25).length;
    const tot = (dir) => d.data.filter((r) => r.direction === dir).length;
    el.textContent = `${SRC[s]}: ${cnt('up')} of the paper's ${tot('up')} up genes and ${cnt('down')} of ${tot('down')} down genes are in the top 25 here (rank by significance among genes passing the paper thresholds). `
      + 'The paper prints “40057” for one gene: that is Sept9 after Excel date conversion.';
  });

  // ----- PCA -----
  dash.draw('cl-pca', (el) => {
    if (!BX.status(el, ['cl_pca'])) return;
    el.replaceChildren();
    const method = P('clpca', 'top500');
    const rows = dash.data('cl_pca').data.filter((r) => r.method === method).map((r) => ({ ...r, x: +r.PC1, y: +r.PC2 }));
    if (!rows.length) { unavailable(el, 'No PCA in this build.'); return; }
    const W = Math.max(300, el.clientWidth), H = Math.min(420, Math.max(300, W * 0.55));
    const m = { t: 30, r: 16, b: 40, l: 52 };
    const svg = BX.svg(el, W, H);
    const pad = (a) => { const e = d3.extent(a); const d = (e[1] - e[0]) * 0.15 || 1; return [e[0] - d, e[1] + d]; };
    const x = d3.scaleLinear(pad(rows.map((r) => r.x)), [m.l, W - m.r]), y = d3.scaleLinear(pad(rows.map((r) => r.y)), [H - m.b, m.t]);
    svg.append('g').selectAll('line').data(y.ticks(5)).join('line').attr('class', 'bx-grid-line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y).attr('y2', y);
    svg.append('g').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(5).tickSizeOuter(0)).call(BX.cleanAxis);
    svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickSizeOuter(0)).call(BX.cleanAxis);
    svg.append('text').attr('x', (m.l + W - m.r) / 2).attr('y', H - 6).attr('text-anchor', 'middle').text(`PC1 (${rows[0].pct_PC1} % variance)`);
    svg.append('text').attr('transform', `translate(12,${(m.t + H - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text(`PC2 (${rows[0].pct_PC2} % variance)`);
    const tip = BX.tip(el);
    svg.append('g').selectAll('circle').data(rows).join('circle').attr('cx', (d) => x(d.x)).attr('cy', (d) => y(d.y)).attr('r', 7)
      .attr('fill', (d) => (d.condition === 'case' ? BX.upColor() : BX.downColor())).attr('fill-opacity', 0.85)
      .on('pointermove', (ev, d) => { const [px, py] = d3.pointer(ev, el); tip.show(px, py, d.label, [d.run, 'PC1 ' + d3.format('.2f')(d.x), 'PC2 ' + d3.format('.2f')(d.y)]); })
      .on('pointerleave', () => tip.hide());
    svg.append('g').selectAll('text').data(rows).join('text').attr('class', 'bx-lab').attr('x', (d) => x(d.x) + 10).attr('y', (d) => y(d.y) + 4)
      .attr('font-size', 10).text((d) => sampleShort(d.sample));
    let lx = m.l + 8;
    for (const [lab, c] of [['S. pneumoniae infected', BX.upColor()], ['Healthy', BX.downColor()]]) {
      svg.append('circle').attr('cx', lx + 5).attr('cy', 10).attr('r', 5).attr('fill', c);
      const t = svg.append('text').attr('x', lx + 14).attr('y', 14).text(lab);
      lx += 26 + t.node().getComputedTextLength();
    }
  });
  dash.draw('cl-pca-note', (el) => {
    const d = dash.data('cl_pca');
    if (d.status !== 'ok' || !d.data.length) { el.textContent = ''; return; }
    const r = d.data.find((q) => q.method === P('clpca', 'top500'));
    const how = r.method === 'top500'
      ? `DESeq2 plotPCA convention: VST values, ${BX.fmtInt(+r.n_genes)} most variable genes, centred.`
      : `All ${BX.fmtInt(+r.n_genes)} tested genes, VST values, each gene scaled to unit variance (FactoMineR default).`;
    el.textContent = `${how} Here PC1 ${r.pct_PC1} % / PC2 ${r.pct_PC2} %; the paper reports Dim 1 ${r.paper_pct_PC1} % / Dim 2 ${r.paper_pct_PC2} % without stating its settings. `
      + 'Both separate infected from healthy samples on PC1, as in the paper.';
  });

  // ----- sample clustering (Supp. Fig. 8) -----
  dash.draw('cl-cluster', (el) => {
    if (!BX.status(el, ['cl_sample_dist', 'cl_dendrogram'])) return;
    el.replaceChildren();
    const dist = dash.data('cl_sample_dist').data;
    const seg = dash.data('cl_dendrogram').data;
    if (!dist.length) { unavailable(el, 'No clustering in this build.'); return; }
    const leaves = Array.from(d3.group(dist, (r) => r.sample_a), ([k, v]) => ({ sample: k, order: +v[0].order_a, label: v[0].label_a, cond: v[0].condition_a }))
      .sort((a, b) => a.order - b.order);
    const n = leaves.length;
    const W = Math.max(300, el.clientWidth);
    const lw = 90, dh = 90;
    const cell = Math.min(46, (W - lw - 8) / n);
    const H = dh + 20 + n * cell + 8;
    const svg = BX.svg(el, W, H);
    const x0 = lw;
    const ymax = d3.max(seg, (r) => Math.max(+r.y1, +r.y2)) || 1;
    const dy = d3.scaleLinear([0, ymax], [dh, 6]);
    const lx = (v) => x0 + (+v + 0.5) * cell;
    for (const r of seg) svg.append('line').attr('class', 'bx-axis-line').attr('x1', lx(r.x1)).attr('x2', lx(r.x2)).attr('y1', dy(+r.y1)).attr('y2', dy(+r.y2));
    leaves.forEach((l, j) => svg.append('rect').attr('x', x0 + j * cell + 1).attr('y', dh + 4).attr('width', cell - 2).attr('height', 8)
      .attr('fill', l.cond === 'case' ? BX.upColor() : muted()));
    const dmax = d3.max(dist, (r) => +r.dist) || 1;
    const col = d3.interpolateRgb(BX.css('--brand') || '#1f4e79', BX.css('--color-panel') || '#fff');
    const tip = BX.tip(el);
    const top = dh + 18;
    leaves.forEach((a, i) => {
      svg.append('text').attr('class', 'bx-lab').attr('x', lw - 6).attr('y', top + i * cell + cell / 2 + 4).attr('text-anchor', 'end').attr('font-size', 10.5).text(a.label);
      leaves.forEach((b, j) => {
        const r = dist.find((q) => q.sample_a === a.sample && q.sample_b === b.sample);
        svg.append('rect').attr('x', x0 + j * cell + 1).attr('y', top + i * cell + 1).attr('width', cell - 2).attr('height', cell - 2).attr('rx', 2)
          .attr('fill', col(+r.dist / dmax))
          .on('pointermove', (ev) => { const [px, py] = d3.pointer(ev, el); tip.show(px, py, `${a.label} ↔ ${b.label}`, [`Euclidean distance ${d3.format('.1f')(+r.dist)}`]); })
          .on('pointerleave', () => tip.hide());
      });
    });
  });

  // ----- pathway enrichment -----
  const pathRows = () => {
    const s = src(), db = P('cldb', 'KEGG'), only = P('clpaper', 'all') === 'paper';
    return dash.data('cl_pathways').data.filter((r) => r.source === s && r.db === db && (!only || r.paper !== ''))
      .sort((a, b) => +a.p_enrich - +b.p_enrich);
  };
  dash.draw('cl-path-chart', (el) => {
    if (!BX.status(el, ['cl_pathways'])) return;
    el.replaceChildren();
    const all = pathRows();
    if (!all.length) { empty(el, 'No pathways for this selection.'); return; }
    const rows = all.slice(0, 25).map((r) => ({ ...r, sc: -Math.log10(Math.max(+r.q_enrich, 1e-300)) }));
    const sel = dash.params().clpanel;
    const W = Math.max(320, el.clientWidth), rowH = 24, H = rows.length * rowH + 48;
    const svg = BX.svg(el, W, H);
    const short = (t) => (t.length > 52 ? t.slice(0, 50) + '…' : t);
    const lw = Math.min(W * 0.5, BX.textWidth(svg, rows.map((r) => short(r.pathway) + ' ★')) + 12);
    const nw = 70;
    const x = d3.scaleLinear([0, Math.max(1.5, d3.max(rows, (r) => r.sc)) * 1.05], [lw, W - nw]);
    svg.append('g').attr('transform', `translate(0,${H - 34})`).call(d3.axisBottom(x).ticks(5).tickSizeOuter(0)).call(BX.cleanAxis);
    svg.append('text').attr('x', (lw + W - nw) / 2).attr('y', H - 2).attr('text-anchor', 'middle').attr('font-size', 10).text('−log10 q (over-representation, BH within database)');
    const q05 = x(-Math.log10(0.05));
    svg.append('line').attr('class', 'bx-ref').attr('x1', q05).attr('x2', q05).attr('y1', 0).attr('y2', H - 34);
    const col = (r) => (r.direction === 'mostly up' ? BX.upColor() : r.direction === 'mostly down' ? BX.downColor() : muted());
    const tip = BX.tip(el);
    rows.forEach((r, i) => {
      const yy = 4 + i * rowH;
      const g = svg.append('g').style('cursor', 'pointer').on('click', () => dash.setParams({ clpanel: r.pathway_id }))
        .on('pointermove', (ev) => { const [px, py] = d3.pointer(ev, el); tip.show(px, py, `${r.pathway} (${r.pathway_id})`, [`${r.n_sig} of ${r.n_genes} genes selected (${r.n_up} up, ${r.n_down} down)`, `fold enrichment ${r.fold_enrichment}`, `q = ${BX.fmtP(r.q_enrich)}; direction q = ${BX.fmtP(r.q_direction)}`, r.paper ? `In the paper (${r.paper === 'text' ? 'named in text / arrow' : 'shown in figure'})` : 'Not in the paper figure']); })
        .on('pointerleave', () => tip.hide());
      g.append('rect').attr('x', 0).attr('y', yy).attr('width', W).attr('height', rowH - 2).attr('fill', r.pathway_id === sel ? (BX.css('--brand-light') || '#eef') : 'transparent');
      g.append('text').attr('class', 'bx-lab').attr('x', lw - 8).attr('y', yy + 15).attr('text-anchor', 'end').attr('font-weight', r.pathway_id === sel || r.paper === 'text' ? 600 : null)
        .text(short(r.pathway) + (r.paper ? ' ★' : ''));
      g.append('rect').attr('x', lw).attr('y', yy + 4).attr('height', 14).attr('width', Math.max(1, x(r.sc) - lw)).attr('rx', 2).attr('fill', col(r));
      g.append('text').attr('x', x(r.sc) + 4).attr('y', yy + 15).attr('font-size', 10).text(`${r.n_sig}/${r.n_genes}`);
    });
  });
  dash.draw('cl-path-note', (el) => {
    const d = dash.data('cl_pathways');
    if (d.status !== 'ok') { el.textContent = ''; return; }
    const s = src(), db = P('cldb', 'KEGG');
    const rows = d.data.filter((r) => r.source === s && r.db === db);
    const paper = rows.filter((r) => r.paper !== '');
    const sel = rows[0] ? BX.fmtInt(+rows[0].n_selected) : '—';
    const how = s === 'pa'
      ? `${sel} genes regulated ≥ 1.5-fold in all four infected samples with p < 0.05 (the paper's KOBAS selection)`
      : `${sel} genes passing padj ≤ 0.05 and |FC| ≥ 1.5 in this run (no per-sample data for the “all four samples” filter)`;
    const tail = db === 'Literature' ? ' Literature panels are BBB gene groups curated for this dashboard, not database pathways.'
      : paper.length ? ` ${paper.filter((r) => +r.q_enrich <= 0.05).length} of the ${paper.length} ${db} pathways in the paper's figure reach q ≤ 0.05 here (★).` : '';
    el.textContent = `Selected: ${how}; universe: ${rows[0] ? BX.fmtInt(+rows[0].universe) : '—'} genes tested. Showing the 25 most significant of ${rows.length} kept pathways.${tail}`;
  });

  dash.draw('cl-path-body', (tb) => {
    tb.replaceChildren();
    if (!ok('cl_pathways')) return;
    const sel = dash.params().clpanel;
    for (const r of pathRows()) {
      const tr = document.createElement('tr');
      tr.dataset.key = r.pathway_id;
      tr.className = r.pathway_id === sel ? 'bx-selected' : '';
      tr.addEventListener('click', () => dash.setParams({ clpanel: r.pathway_id }));
      td(tr, r.pathway);
      td(tr, r.paper === 'text' ? 'named / arrow' : r.paper === 'figure' ? 'in figure' : '—', r.paper ? '' : 'bx-muted');
      td(tr, `${r.n_sig} / ${r.n_genes}`, 'bx-num');
      td(tr, `${r.n_up} / ${r.n_down}`, 'bx-num');
      td(tr, d3.format('.2f')(+r.fold_enrichment), 'bx-num');
      td(tr, BX.fmtP(r.q_enrich), 'bx-num' + (+r.q_enrich <= 0.05 ? ' bx-strong' : ''));
      td(tr, r.direction, r.direction === 'mostly up' ? 'bx-up' : r.direction === 'mostly down' ? 'bx-down' : 'bx-muted');
      td(tr, BX.fmtP(r.q_direction), 'bx-num');
      tb.appendChild(tr);
    }
  });

  dash.draw('cl-panel-select', (s) => {
    if (!ok('cl_pathways')) return;
    const rows = pathRows();
    const key = src() + '|' + P('cldb', 'KEGG') + '|' + P('clpaper', 'all');
    if (s.dataset.built !== key) {
      s.replaceChildren();
      for (const r of rows) { const o = document.createElement('option'); o.value = r.pathway_id; o.textContent = r.pathway; s.appendChild(o); }
      s.dataset.built = key;
    }
    const p = dash.params().clpanel;
    if (!rows.length) return;
    if (!p || !rows.some((r) => r.pathway_id === p)) { dash.setParams({ clpanel: rows[0].pathway_id }); return; }
    s.value = p;
  });

  dash.draw('cl-panel-info', (el) => {
    if (!ok('cl_pathways')) { el.textContent = ''; return; }
    const p = dash.params().clpanel;
    const r = dash.data('cl_pathways').data.find((x) => x.pathway_id === p && x.source === src());
    if (!r) { el.textContent = ''; return; }
    const link = { KEGG: `https://www.kegg.jp/pathway/${r.pathway_id}`, Reactome: `https://reactome.org/content/detail/${r.pathway_id}`, PANTHER: `https://www.pantherdb.org/pathway/pathDetail.do?clsAccession=${r.pathway_id}` }[r.db];
    el.replaceChildren(document.createTextNode(`${r.db} ${r.pathway_id}: ${r.n_sig} of ${r.n_genes} tested genes selected; q = ${BX.fmtP(r.q_enrich)}. `));
    if (link) { const a = document.createElement('a'); a.href = link; a.target = '_blank'; a.rel = 'noopener'; a.textContent = 'Open in ' + r.db; el.appendChild(a); }
  });

  dash.draw('cl-members-body', (tb) => {
    tb.replaceChildren();
    if (!ok('cl_pathway_genes') || !ok('cl_deg')) return;
    const s = src();
    const p = dash.params().clpanel;
    const idx = degById();
    const order = { up: 0, down: 1, below_fc: 2, '': 3 };
    const rows = dash.data('cl_pathway_genes').data.filter((r) => r.pathway_id === p).map((r) => idx.get(r.gene_id) || { gene_id: r.gene_id, symbol: r.gene_id })
      .sort((a, b) => (order[a[s + '_call'] || ''] - order[b[s + '_call'] || '']) || Math.abs(+b[s + '_log2fc'] || 0) - Math.abs(+a[s + '_log2fc'] || 0));
    const lbl = { up: 'up', down: 'down', below_fc: '|FC| < 1.5', '': 'n.s.' };
    for (const r of rows) {
      const tr = document.createElement('tr');
      tr.dataset.key = p + ',' + r.gene_id;
      tr.dataset.gene = r.symbol;
      td(tr, r.symbol);
      const v = r[s + '_log2fc'] ?? '';
      td(tr, v === '' ? '—' : BX.fmtFc(v), fcCls(v));
      td(tr, (r[s + '_padj'] ?? '') === '' ? '—' : BX.fmtP(r[s + '_padj']), 'bx-num');
      const c = r[s + '_call'] || '';
      td(tr, lbl[c] ?? c, c === 'up' ? 'bx-up' : c === 'down' ? 'bx-down' : 'bx-muted');
      td(tr, r.pa_all4 === '1' ? 'up in all 4' : r.pa_all4 === '-1' ? 'down in all 4' : '—', 'bx-muted');
      tb.appendChild(tr);
    }
  });
})();
