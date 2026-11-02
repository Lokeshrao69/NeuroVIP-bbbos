// ---------- Gene lookup, Cross-disease, Pipeline check, ML drawings ----------
(function () {
  const BX = window.BX;
  const geneRows = () => {
    const g = (dash.params().gene || 'Selp').toLowerCase();
    return dash.data('compare').data.filter((r) => r.symbol.toLowerCase() === g || r.gene_id.toLowerCase() === g);
  };
  const empty = (el, msg) => { const e = document.createElement('div'); e.className = 'bx-empty'; e.textContent = msg; el.appendChild(e); };

  // grouped bar chart: per disease, this pipeline vs classic
  dash.draw('gene-chart', (el) => {
    if (!BX.status(el, ['compare'])) return;
    el.replaceChildren();
    const gene = dash.params().gene || 'Selp';
    const rows = geneRows();
    if (!rows.length) { empty(el, `"${gene}" is not significant in any disease (or the symbol was not found).`); return; }
    const series = [{ key: 'ai', name: 'This pipeline', color: dash.colors[0] }, { key: 'mk', name: 'Classic analysis', color: dash.colors[2] || dash.colors[1] }];
    const data = [];
    for (const dz of BX.ORDER) {
      const r = rows.find((x) => x.disease === dz);
      for (const s of series) {
        const v = r ? BX.num(r[s.key + '_log2fc']) : null;
        const p = r ? r[s.key + '_padj'] : '';
        data.push({ dz, s, v, p, symbol: r ? r.symbol : gene });
      }
    }
    const W = Math.max(320, el.clientWidth), H = 260;
    const m = { t: 28, r: 12, b: 34, l: 44 };
    const svg = BX.svg(el, W, H);
    const vals = data.map((d) => d.v).filter((v) => v != null);
    const ext = Math.max(1, d3.max(vals, Math.abs) || 1) * 1.15;
    const x0 = d3.scaleBand(BX.ORDER, [m.l, W - m.r]).padding(0.25);
    const x1 = d3.scaleBand(series.map((s) => s.key), [0, x0.bandwidth()]).padding(0.12);
    const y = d3.scaleLinear([d3.min(vals) < 0 ? -ext : 0, d3.max(vals) > 0 ? ext : 0.5], [H - m.b, m.t]).nice();
    svg.append('g').selectAll('line').data(y.ticks(5)).join('line').attr('class', 'bx-grid-line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y).attr('y2', y);
    svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickSizeOuter(0)).call(BX.cleanAxis);
    svg.append('line').attr('class', 'bx-axis-line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(0)).attr('y2', y(0));
    svg.append('g').selectAll('text').data(BX.ORDER).join('text').attr('x', (d) => x0(d) + x0.bandwidth() / 2).attr('y', H - 10).attr('text-anchor', 'middle').text((d) => BX.NAMES[d]);
    const bw = Math.min(32, x1.bandwidth());
    const tip = BX.tip(el);
    const g = svg.append('g');
    for (const d of data) {
      const cx = x0(d.dz) + x1(d.s.key) + (x1.bandwidth() - bw) / 2;
      if (d.v == null) {
        g.append('text').attr('x', cx + bw / 2).attr('y', y(0) - 6).attr('text-anchor', 'middle').attr('font-size', 10).text('n.s.');
        continue;
      }
      g.append('rect').attr('x', cx).attr('width', bw).attr('y', Math.min(y(0), y(d.v))).attr('height', Math.max(1, Math.abs(y(d.v) - y(0)))).attr('rx', 3).attr('fill', d.s.color)
        .on('pointermove', (ev) => { const [mx, my] = d3.pointer(ev, el); tip.show(mx, my, `${d.symbol} · ${BX.NAMES[d.dz]}`, [`${d.s.name}`, `log2FC ${BX.fmtFc(d.v)}`, `padj ${BX.fmtP(d.p)}`]); })
        .on('pointerleave', () => tip.hide());
      g.append('text').attr('class', 'bx-lab').attr('x', cx + bw / 2).attr('y', d.v >= 0 ? y(d.v) - 4 : y(d.v) + 12).attr('text-anchor', 'middle').attr('font-size', 10).text(BX.fmtFc(d.v));
    }
    // legend
    let lx = m.l;
    for (const s of series) {
      svg.append('rect').attr('x', lx).attr('y', 6).attr('width', 10).attr('height', 10).attr('rx', 2).attr('fill', s.color);
      const t = svg.append('text').attr('x', lx + 14).attr('y', 15).text(s.name);
      lx += 24 + t.node().getComputedTextLength();
    }
  });

  dash.draw('gene-body', (tb) => {
    tb.replaceChildren();
    if (dash.data('compare').status !== 'ok') return;
    const rows = geneRows().sort((a, b) => BX.ORDER.indexOf(a.disease) - BX.ORDER.indexOf(b.disease));
    for (const r of rows) {
      const tr = document.createElement('tr');
      tr.dataset.key = encodeURIComponent(r.disease) + ',' + encodeURIComponent(r.gene_id);
      [[BX.NAMES[r.disease]], [r.symbol], [r.gene_id, 'bx-muted'], [BX.fmtFc(r.ai_log2fc), 'bx-num'], [BX.fmtP(r.ai_padj), 'bx-num'], [BX.fmtFc(r.mk_log2fc), 'bx-num'], [BX.fmtP(r.mk_padj), 'bx-num']].forEach(([v, cls]) => {
        const td = document.createElement('td'); td.textContent = v; if (cls) td.className = cls; tr.appendChild(td);
      });
      tb.appendChild(tr);
    }
  });

  // UpSet-style exclusive intersections
  dash.draw('upset', (el) => {
    if (!BX.status(el, ['overlap'])) return;
    el.replaceChildren();
    const LBL = { stroke_only: ['stroke'], meningitis_only: ['meningitis'], gb_only: ['gb'], stroke_meningitis: ['stroke', 'meningitis'], stroke_gb: ['stroke', 'gb'], meningitis_gb: ['meningitis', 'gb'], all_three: ['stroke', 'meningitis', 'gb'] };
    const rows = dash.data('overlap').data.slice().sort((a, b) => b.n - a.n);
    const W = Math.max(300, el.clientWidth), dotH = 18, barH = 170, H = barH + 3 * dotH + 30;
    const svg = BX.svg(el, W, H);
    const lw = BX.textWidth(svg, BX.ORDER.map((d) => BX.NAMES[d])) + 10;
    const m = { l: lw, r: 8, t: 16 };
    const x = d3.scaleBand(rows.map((r) => r.set), [m.l, W - m.r]).padding(0.3);
    const y = d3.scaleLinear([0, d3.max(rows, (r) => r.n)], [barH, m.t]);
    const accent = dash.colors[0], muted = BX.css('--cds-chart-muted') || '#bbb', fg = BX.css('--color-fg') || '#222';
    const bw = Math.min(32, x.bandwidth());
    const tip = BX.tip(el);
    for (const r of rows) {
      const cx = x(r.set) + (x.bandwidth() - bw) / 2;
      svg.append('rect').attr('x', cx).attr('width', bw).attr('y', y(r.n)).attr('height', barH - y(r.n)).attr('rx', 3).attr('fill', r.set === 'all_three' ? accent : muted)
        .on('pointermove', (ev) => { const [mx, my] = d3.pointer(ev, el); tip.show(mx, my, LBL[r.set].map((d) => BX.NAMES[d]).join(' + ') + (LBL[r.set].length === 1 ? ' only' : ''), [BX.fmtInt(r.n) + ' genes']); })
        .on('pointerleave', () => tip.hide());
      svg.append('text').attr('class', 'bx-lab').attr('x', cx + bw / 2).attr('y', y(r.n) - 4).attr('text-anchor', 'middle').attr('font-size', 10).text(BX.fmtInt(r.n));
      const cxm = cx + bw / 2;
      const ys = BX.ORDER.map((d, i) => barH + 16 + i * dotH);
      BX.ORDER.forEach((d, i) => svg.append('circle').attr('cx', cxm).attr('cy', ys[i]).attr('r', 5).attr('fill', LBL[r.set].includes(d) ? fg : muted).attr('fill-opacity', LBL[r.set].includes(d) ? 1 : 0.35));
      const on = BX.ORDER.map((d, i) => (LBL[r.set].includes(d) ? ys[i] : null)).filter((v) => v != null);
      if (on.length > 1) svg.append('line').attr('x1', cxm).attr('x2', cxm).attr('y1', d3.min(on)).attr('y2', d3.max(on)).attr('stroke', fg).attr('stroke-width', 2);
    }
    BX.ORDER.forEach((d, i) => svg.append('text').attr('x', 0).attr('y', barH + 20 + i * dotH).text(BX.NAMES[d]));
  });

  // direction of shared genes: horizontal bars
  dash.draw('direction', (el) => {
    if (!BX.status(el, ['shared_direction'])) return;
    el.replaceChildren();
    const rows = dash.data('shared_direction').data;
    const W = Math.max(280, el.clientWidth), rowH = 34, H = rows.length * rowH + 8;
    const svg = BX.svg(el, W, H);
    const lw = BX.textWidth(svg, rows.map((r) => r.direction)) + 10;
    const nw = BX.textWidth(svg, rows.map((r) => BX.fmtInt(r.n))) + 8;
    const x = d3.scaleLinear([0, d3.max(rows, (r) => r.n) || 1], [lw, W - nw]);
    const col = { 'up in all': BX.upColor(), 'down in all': BX.downColor(), mixed: BX.css('--cds-chart-muted') || '#bbb' };
    rows.forEach((r, i) => {
      const yy = i * rowH + 6;
      svg.append('text').attr('x', 0).attr('y', yy + 15).text(r.direction);
      svg.append('rect').attr('x', lw).attr('y', yy + 2).attr('height', 18).attr('width', Math.max(1, x(r.n) - lw)).attr('rx', 3).attr('fill', col[r.direction]);
      svg.append('text').attr('class', 'bx-lab').attr('x', x(r.n) + 4).attr('y', yy + 15).text(BX.fmtInt(r.n));
    });
  });

  // heatmap of the strongest shared genes
  dash.draw('heatmap', (el) => {
    if (!BX.status(el, ['shared_all'])) return;
    el.replaceChildren();
    const rows = dash.data('shared_all').data.slice(0, 40);
    if (!rows.length) { empty(el, 'No shared genes'); return; }
    const W = Math.max(300, el.clientWidth);
    const narrow = W < 640;
    const cellH = narrow ? 14 : 16;
    const svg = BX.svg(el, W, 10);
    const gw = BX.textWidth(svg, rows.map((r) => r.symbol)) + 10;
    const H = rows.length * cellH + 40;
    svg.attr('viewBox', `0 0 ${W} ${H}`);
    const cw = Math.min(120, (W - gw - 8) / 3);
    const lim = Math.min(8, d3.max(rows, (r) => Math.max(Math.abs(r.stroke), Math.abs(r.meningitis), Math.abs(r.gb))));
    const base = BX.css('--color-panel') || '#fff';
    const pos = d3.interpolateRgb(base, BX.upColor()), neg = d3.interpolateRgb(base, BX.downColor());
    const color = (v) => (v >= 0 ? pos(Math.min(1, v / lim)) : neg(Math.min(1, -v / lim)));
    BX.ORDER.forEach((d, j) => svg.append('text').attr('x', gw + j * cw + cw / 2).attr('y', 14).attr('text-anchor', 'middle').text(BX.NAMES[d]));
    const tip = BX.tip(el);
    rows.forEach((r, i) => {
      const yy = 24 + i * cellH;
      svg.append('text').attr('class', 'bx-lab').attr('x', gw - 6).attr('y', yy + cellH - 4).attr('text-anchor', 'end').attr('font-size', 10.5).text(r.symbol);
      BX.ORDER.forEach((d, j) => {
        svg.append('rect').attr('x', gw + j * cw + 1).attr('y', yy + 1).attr('width', cw - 2).attr('height', cellH - 2).attr('rx', 2).attr('fill', color(r[d]))
          .on('pointermove', (ev) => { const [mx, my] = d3.pointer(ev, el); tip.show(mx, my, `${r.symbol} · ${BX.NAMES[d]}`, ['log2FC ' + BX.fmtFc(r[d])]); })
          .on('pointerleave', () => tip.hide());
        if (cw > 46) svg.append('text').attr('x', gw + j * cw + cw / 2).attr('y', yy + cellH - 4).attr('text-anchor', 'middle').attr('font-size', 9.5).attr('pointer-events', 'none')
          .style('fill', Math.abs(r[d]) / lim > 0.6 ? '#fff' : null).text(BX.fmtFc(r[d]));
      });
    });
  });

  dash.draw('shared-body', (tb) => {
    tb.replaceChildren();
    const d = dash.data('shared_all');
    if (d.status !== 'ok') return;
    const q = BX.state.sharedSearch;
    const rows = d.data.filter((r) => !q || r.symbol.toLowerCase().includes(q) || r.gene_id.toLowerCase().includes(q)).slice(0, 400);
    for (const r of rows) {
      const tr = document.createElement('tr');
      tr.dataset.key = r.gene_id;
      tr.dataset.gene = r.symbol;
      const cells = [[r.symbol]];
      for (const dz of BX.ORDER) cells.push([BX.fmtFc(r[dz]), 'bx-num ' + (r[dz] > 0 ? 'bx-up' : 'bx-down')]);
      cells.push([r.direction, 'bx-muted']);
      cells.forEach(([v, cls]) => { const td = document.createElement('td'); td.textContent = v; if (cls) td.className = cls; tr.appendChild(td); });
      tb.appendChild(tr);
    }
  });

  // pipeline check: stacked bars per disease
  dash.draw('concord-bars', (el) => {
    if (!BX.status(el, ['concord'])) return;
    el.replaceChildren();
    const rows = dash.data('concord').data;
    const parts = [{ k: 'both', name: 'Both', color: dash.colors[0] }, { k: 'ai_only', name: 'This pipeline only', color: dash.colors[1] }, { k: 'classic_only', name: 'Classic only', color: dash.colors[3] || dash.colors[2] }];
    const W = Math.max(300, el.clientWidth), rowH = 36, H = rows.length * rowH + 34;
    const svg = BX.svg(el, W, H);
    const lw = BX.textWidth(svg, rows.map((r) => BX.NAMES[r.disease])) + 10;
    const tot = (r) => parts.reduce((s, p) => s + r[p.k], 0);
    const nw = BX.textWidth(svg, rows.map((r) => BX.fmtInt(tot(r)))) + 8;
    const x = d3.scaleLinear([0, d3.max(rows, tot)], [lw, W - nw]);
    const tip = BX.tip(el);
    rows.forEach((r, i) => {
      const yy = 30 + i * rowH;
      svg.append('text').attr('x', 0).attr('y', yy + 14).text(BX.NAMES[r.disease]);
      let acc = 0;
      for (const p of parts) {
        const v = r[p.k];
        if (v > 0) svg.append('rect').attr('x', x(acc)).attr('y', yy).attr('height', 20).attr('width', Math.max(1, x(acc + v) - x(acc))).attr('fill', p.color)
          .on('pointermove', (ev) => { const [mx, my] = d3.pointer(ev, el); tip.show(mx, my, BX.NAMES[r.disease], [`${p.name}: ${BX.fmtInt(v)}`]); })
          .on('pointerleave', () => tip.hide());
        acc += v;
      }
      svg.append('text').attr('class', 'bx-lab').attr('x', x(acc) + 4).attr('y', yy + 14).text(BX.fmtInt(acc));
    });
    let lx = lw;
    for (const p of parts) {
      svg.append('rect').attr('x', lx).attr('y', 4).attr('width', 10).attr('height', 10).attr('rx', 2).attr('fill', p.color);
      const t = svg.append('text').attr('x', lx + 14).attr('y', 13).text(p.name);
      lx += 26 + t.node().getComputedTextLength();
    }
  });

  dash.draw('ck', (el) => {
    const disease = dash.params().disease || 'stroke';
    el.setAttribute('data-where', 'disease=' + disease);
    const d = dash.data('concord');
    if (d.status !== 'ok') { el.textContent = '—'; return; }
    const row = d.data.find((r) => r.disease === disease);
    const f = el.dataset.field;
    el.textContent = !row ? '—' : f === 'r' ? d3.format('.5f')(row.r) : BX.fmtInt(row[f]);
  });

  dash.draw('ck-scatter', (el) => {
    if (!BX.status(el, ['compare'])) return;
    el.replaceChildren();
    const disease = dash.params().disease || 'stroke';
    const rows = dash.data('compare').data.filter((r) => r.disease === disease && r.found_by === 'both').map((r) => ({ symbol: r.symbol, a: +r.ai_log2fc, b: +r.mk_log2fc }));
    const W = Math.max(300, el.clientWidth), H = Math.min(460, Math.max(320, W * 0.6));
    const m = { t: 10, r: 14, b: 40, l: 48 };
    const svg = BX.svg(el, W, H);
    const ext = d3.max(rows, (d) => Math.max(Math.abs(d.a), Math.abs(d.b))) * 1.05 || 1;
    const x = d3.scaleLinear([-ext, ext], [m.l, W - m.r]), y = d3.scaleLinear([-ext, ext], [H - m.b, m.t]);
    svg.append('g').selectAll('line').data(y.ticks(6)).join('line').attr('class', 'bx-grid-line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y).attr('y2', y);
    svg.append('g').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(Math.max(3, Math.floor(W / 80))).tickSizeOuter(0)).call(BX.cleanAxis);
    svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6).tickSizeOuter(0)).call(BX.cleanAxis);
    svg.append('line').attr('class', 'bx-ref').attr('x1', x(-ext)).attr('y1', y(-ext)).attr('x2', x(ext)).attr('y2', y(ext));
    svg.append('text').attr('x', (m.l + W - m.r) / 2).attr('y', H - 6).attr('text-anchor', 'middle').text('log2FC, this pipeline');
    svg.append('text').attr('transform', `translate(12,${(m.t + H - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text('log2FC, classic analysis');
    svg.append('g').selectAll('circle').data(rows).join('circle').attr('cx', (d) => x(d.a)).attr('cy', (d) => y(d.b)).attr('r', 2.2).attr('fill', dash.colors[0]).attr('fill-opacity', 0.5);
    const tip = BX.tip(el);
    const hl = svg.append('circle').attr('r', 5).attr('fill', 'none').attr('stroke', BX.css('--color-fg') || '#000').attr('visibility', 'hidden');
    const dl = d3.Delaunay.from(rows, (d) => x(d.a), (d) => y(d.b));
    svg.on('pointermove', (ev) => {
      const [mx, my] = d3.pointer(ev);
      const d = rows[dl.find(mx, my)];
      if (!d || Math.hypot(x(d.a) - mx, y(d.b) - my) > 24) { hl.attr('visibility', 'hidden'); tip.hide(); return; }
      hl.attr('cx', x(d.a)).attr('cy', y(d.b)).attr('visibility', 'visible');
      const k = el.clientWidth / W;
      tip.show(x(d.a) * k, y(d.b) * k, d.symbol, ['this pipeline ' + BX.fmtFc(d.a), 'classic ' + BX.fmtFc(d.b)]);
    }).on('pointerleave', () => { hl.attr('visibility', 'hidden'); tip.hide(); });
  });

  // ML
  dash.draw('ml-body', (tb) => {
    tb.replaceChildren();
    const d = dash.data('ml_metrics');
    if (d.status !== 'ok') return;
    const f = (v) => (v === '' || v == null ? '—' : d3.format('.2f')(+v));
    for (const r of d.data) {
      const tr = document.createElement('tr');
      tr.dataset.key = r.id;
      [[r.model], [r.evaluation], [r.n_test, 'bx-num'], [f(r.accuracy), 'bx-num'], [f(r.balanced_accuracy), 'bx-num'], [f(r.roc_auc), 'bx-num']].forEach(([v, cls]) => {
        const td = document.createElement('td'); td.textContent = v; if (cls) td.className = cls; tr.appendChild(td);
      });
      tb.appendChild(tr);
    }
  });

  dash.draw('feat', (el) => {
    if (!BX.status(el, ['ml_features'])) return;
    el.replaceChildren();
    const which = new URLSearchParams(el.dataset.where || '').get('analysis');
    const rows = dash.data('ml_features').data.filter((r) => r.analysis === which).map((r) => ({ s: r.symbol, v: +r.coefficient })).sort((a, b) => b.v - a.v);
    const W = Math.max(260, el.clientWidth), rowH = 20, H = rows.length * rowH + 8;
    const svg = BX.svg(el, W, H);
    const lw = BX.textWidth(svg, rows.map((r) => r.s)) + 10;
    const nw = BX.textWidth(svg, rows.map((r) => d3.format('.2f')(r.v))) + 8;
    const x = d3.scaleLinear([0, d3.max(rows, (r) => r.v) || 1], [lw, W - nw]);
    rows.forEach((r, i) => {
      const yy = i * rowH + 4;
      svg.append('text').attr('class', 'bx-lab').attr('x', lw - 6).attr('y', yy + 12).attr('text-anchor', 'end').text(r.s);
      svg.append('rect').attr('x', lw).attr('y', yy + 2).attr('height', 13).attr('width', Math.max(1, x(r.v) - lw)).attr('rx', 2).attr('fill', dash.colors[0]);
      svg.append('text').attr('x', x(r.v) + 4).attr('y', yy + 12).text(d3.format('.2f')(r.v));
    });
  });
})();
