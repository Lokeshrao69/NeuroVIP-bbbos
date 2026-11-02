// ---------- Overview + Differential expression drawings ----------
(function () {
  const BX = window.BX;
  const whereOf = (el) => Object.fromEntries(new URLSearchParams(el.dataset.where || ''));

  // headline numbers: tiles read their own data-where
  dash.draw('tile', (el) => {
    const id = el.dataset.source, field = el.dataset.field, w = whereOf(el);
    const d = dash.data(id);
    if (d.status !== 'ok') { el.textContent = '—'; return; }
    const row = d.data.find((r) => Object.entries(w).every(([k, v]) => String(r[k]) === v));
    el.textContent = row ? BX.fmtInt(row[field]) : '—';
  });

  dash.draw('take-shared', (el) => {
    const d = dash.data('overlap');
    if (d.status !== 'ok') { el.textContent = '—'; return; }
    const row = d.data.find((r) => r.set === 'all_three');
    el.textContent = row ? BX.fmtInt(row.n) : '—';
  });

  // datasets table
  dash.draw('studies-body', (tb) => {
    const d = dash.data('studies');
    tb.replaceChildren();
    if (d.status !== 'ok') return;
    for (const r of d.data) {
      const tr = document.createElement('tr');
      tr.dataset.key = r.disease;
      const cells = [r.disease_name, r.model, r.tissue, r.study, r.geo || '—', r.n_case, r.n_control, BX.fmtInt(+r.genes_tested)];
      cells.forEach((v, i) => {
        const td = document.createElement('td');
        td.textContent = v;
        if (i >= 5) td.className = 'bx-num';
        tr.appendChild(td);
      });
      tb.appendChild(tr);
    }
  });

  // QC table with an in-cell bar for assigned %
  dash.draw('qc-body', (tb) => {
    const d = dash.data('samples');
    tb.replaceChildren();
    if (d.status !== 'ok') return;
    const LOW = 70; // assigned-to-genes % below this is flagged
    for (const r of d.data) {
      const tr = document.createElement('tr');
      const td = (txt, cls) => { const c = document.createElement('td'); c.textContent = txt; if (cls) c.className = cls; tr.appendChild(c); return c; };
      td(r.run);
      td(BX.NAMES[r.disease] || r.disease);
      td(r.group);
      td(BX.fmtInt(+r.raw_reads), 'bx-num');
      td(d3.format('.1f')(+r.adapter_pct), 'bx-num');
      td(d3.format('.1f')(+r.unique_pct), 'bx-num');
      const c = document.createElement('td');
      c.className = 'bx-num';
      const wrap = document.createElement('span');
      wrap.className = 'bx-bar-cell';
      const bar = document.createElement('span');
      const v = +r.assigned_pct;
      bar.className = 'bx-bar' + (v < LOW ? ' bx-warn' : '');
      bar.style.width = Math.round(v * 0.6) + 'px';
      const num = document.createElement('span');
      num.textContent = d3.format('.1f')(v);
      if (v < LOW) num.className = 'bx-warn-num';
      wrap.append(bar, num);
      c.appendChild(wrap);
      tr.appendChild(c);
      tb.appendChild(tr);
    }
  });

  // DE summary sentence for the selected disease
  dash.draw('de-sum', (el) => {
    const disease = dash.params().disease || 'stroke';
    el.setAttribute('data-where', 'disease=' + disease);
    const d = dash.data('deg_summary');
    if (d.status !== 'ok') { el.textContent = '—'; return; }
    const row = d.data.find((r) => r.disease === disease);
    el.textContent = row ? BX.fmtInt(row[el.dataset.field]) : '—';
  });

  // volcano plot of significant genes
  dash.draw('volcano', (el) => {
    if (!BX.status(el, ['deg'])) return;
    el.replaceChildren();
    const disease = dash.params().disease || 'stroke';
    const lfc = +(dash.params().lfc ?? 1);
    const rows = dash.data('deg').data.filter((r) => r.disease === disease).map((r) => {
      const p = +r.padj;
      return { symbol: r.symbol, gene_id: r.gene_id, fc: +r.log2fc, padj: p, y: -Math.log10(Math.max(p, 1e-320)) };
    });
    if (!rows.length) { const e = document.createElement('div'); e.className = 'bx-empty'; e.textContent = 'No genes'; el.appendChild(e); return; }
    const W = Math.max(320, el.clientWidth), H = Math.max(360, Math.min(520, W * 0.55));
    const m = { t: 12, r: 16, b: 40, l: 48 };
    const svg = BX.svg(el, W, H);
    const xMax = d3.max(rows, (d) => Math.abs(d.fc)) * 1.05;
    const x = d3.scaleLinear([-xMax, xMax], [m.l, W - m.r]);
    const y = d3.scaleLinear([0, d3.max(rows, (d) => d.y) * 1.06], [H - m.b, m.t]);
    svg.append('g').selectAll('line').data(y.ticks(6)).join('line').attr('class', 'bx-grid-line')
      .attr('x1', m.l).attr('x2', W - m.r).attr('y1', y).attr('y2', y);
    svg.append('g').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(Math.max(3, Math.floor(W / 80))).tickSizeOuter(0)).call(BX.cleanAxis);
    svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6).tickSizeOuter(0)).call(BX.cleanAxis);
    svg.append('text').attr('x', (m.l + W - m.r) / 2).attr('y', H - 6).attr('text-anchor', 'middle').text('log2 fold change (case vs control)');
    svg.append('text').attr('transform', `translate(12,${(m.t + H - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text('−log10 padj');
    if (lfc > 0) for (const v of [-lfc, lfc]) if (Math.abs(v) < xMax)
      svg.append('line').attr('class', 'bx-ref').attr('x1', x(v)).attr('x2', x(v)).attr('y1', m.t).attr('y2', H - m.b);
    const up = BX.upColor(), down = BX.downColor(), muted = BX.css('--cds-chart-muted') || '#bbb';
    const col = (d) => (Math.abs(d.fc) < lfc ? muted : d.fc > 0 ? up : down);
    const sorted = rows.slice().sort((a, b) => (Math.abs(a.fc) >= lfc) - (Math.abs(b.fc) >= lfc));
    svg.append('g').selectAll('circle').data(sorted).join('circle')
      .attr('cx', (d) => x(d.fc)).attr('cy', (d) => y(d.y)).attr('r', 2.4)
      .attr('fill', col).attr('fill-opacity', 0.7);
    // labels for strongest genes
    const lab = rows.filter((d) => Math.abs(d.fc) >= lfc).sort((a, b) => b.y * Math.abs(b.fc) - a.y * Math.abs(a.fc)).slice(0, 12);
    const placed = [];
    for (const d of lab) {
      const px = x(d.fc), py = y(d.y) - 5;
      const t = svg.append('text').attr('class', 'bx-lab').attr('x', px).attr('y', py)
        .attr('text-anchor', px > W - 80 ? 'end' : px < m.l + 60 ? 'start' : 'middle').text(d.symbol);
      const b = t.node().getBBox();
      if (placed.some((p) => !(b.x + b.width < p.x || p.x + p.width < b.x || b.y + b.height < p.y || p.y + p.height < b.y))) t.remove();
      else placed.push(b);
    }
    const hl = svg.append('circle').attr('r', 5).attr('fill', 'none').attr('stroke', BX.css('--color-fg') || '#000').attr('visibility', 'hidden');
    const tip = BX.tip(el);
    const delaunay = d3.Delaunay.from(rows, (d) => x(d.fc), (d) => y(d.y));
    svg.on('pointermove', (ev) => {
      const [mx, my] = d3.pointer(ev);
      const i = delaunay.find(mx, my);
      const d = rows[i];
      if (!d || Math.hypot(x(d.fc) - mx, y(d.y) - my) > 24) { hl.attr('visibility', 'hidden'); tip.hide(); return; }
      hl.attr('cx', x(d.fc)).attr('cy', y(d.y)).attr('visibility', 'visible');
      const k = el.clientWidth / W;
      tip.show(x(d.fc) * k, y(d.y) * k, d.symbol, [d.gene_id, 'log2FC ' + BX.fmtFc(d.fc), 'padj ' + BX.fmtP(d.padj)]);
    }).on('pointerleave', () => { hl.attr('visibility', 'hidden'); tip.hide(); });
  });

  // DEG table for the selected disease
  dash.draw('de-body', (tb) => {
    const d = dash.data('deg');
    tb.replaceChildren();
    if (d.status !== 'ok') return;
    const disease = dash.params().disease || 'stroke';
    const lfc = +(dash.params().lfc ?? 1);
    const q = BX.state.deSearch;
    const rows = d.data
      .filter((r) => r.disease === disease && Math.abs(+r.log2fc) >= lfc)
      .filter((r) => !q || r.symbol.toLowerCase().includes(q) || r.gene_id.toLowerCase().includes(q))
      .sort((a, b) => +a.padj - +b.padj || Math.abs(+b.log2fc) - Math.abs(+a.log2fc))
      .slice(0, 250);
    for (const r of rows) {
      const tr = document.createElement('tr');
      tr.dataset.key = encodeURIComponent(r.disease) + ',' + encodeURIComponent(r.gene_id);
      tr.dataset.gene = r.symbol;
      const fc = +r.log2fc;
      [[r.symbol], [r.gene_id, 'bx-muted'], [BX.fmtFc(fc), 'bx-num ' + (fc > 0 ? 'bx-up' : 'bx-down')], [BX.fmtP(r.padj), 'bx-num']].forEach(([v, cls]) => {
        const td = document.createElement('td');
        td.textContent = v;
        if (cls) td.className = cls;
        tr.appendChild(td);
      });
      tb.appendChild(tr);
    }
  });
})();
