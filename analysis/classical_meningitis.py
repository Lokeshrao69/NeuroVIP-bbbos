#!/usr/bin/env python3
"""
Classical analysis of the meningitis brain-microvessel RNA-seq (GEO GSE122952),
following the computational methods of:

    Devraj G, et al. HIF-1α is involved in blood–brain barrier dysfunction and
    paracellular migration of bacteria in pneumococcal meningitis.
    Acta Neuropathol. 2020;140:183–208. doi:10.1007/s00401-020-02174-2

Three DESeq2 runs of the same 8 samples are compared (column prefixes):
  pa_  the paper's own count matrix from GEO (GSE122952_counts.matrix.norm_anno.txt.gz,
       mm10 / STAR 2.4 / featureCounts), re-analysed here with DESeq2 (pydeseq2) and the
       paper's "combined mean >= 5 reads" filter. This is the closest reproduction.
  nv_  the NeuroVIP pipeline (GRCm39 / Ensembl 115), from data/ai_vs_classic.csv
  cl_  the classic analysis (same alignment), from data/ai_vs_classic.csv
The nv_ and cl_ runs only store significant genes; everything per-sample (PCA, clustering,
read counts, per-sample heatmap, the "regulated in all four samples" filter) comes from the
paper's matrix.

Paper method                                   -> here
  DESeq2; BH padj <= 0.05; |FC| >= 1.5; mean>=5 -> pydeseq2 on the GEO matrix, same thresholds
  PCA (Fig. 6b)                                 -> VST, top 500 variable genes (plotPCA) and all genes
  Hierarchical clustering (Supp. Fig. 8)        -> Euclidean distance of VST, complete linkage
  Top 25 up / down Z-score heatmap (Fig. 6c)    -> same, plus overlap with the paper's gene list
  Read counts of named genes (Fig. 6d)          -> normalised counts per sample from the GEO matrix
  KOBAS: KEGG, PANTHER, Reactome (Fig. 6e–f)    -> hypergeometric test + BH on official gene sets,
                                                   genes regulated in all four infected samples

Downloads (cached in data/raw/, which is git-ignored; use --offline to forbid downloads):
  GEO count matrix, KEGG REST (pathways, genes), Reactome Ensembl2Reactome (mouse rows),
  PANTHER SequenceAssociationPathway (mouse rows), MGI MRK_ENSEMBL.rpt (ID mapping).

Outputs: data/classical_*.csv. No result values are typed into this file; the only numbers
from the paper are what it reports (PAPER_REPORTED, PAPER_TOP25, PAPER_PATHWAYS), used for
comparison only.
"""
from __future__ import annotations

import argparse
import io
import math
import os
import sys
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform

DISEASE = "meningitis"
SOURCES = {"pa": "Paper counts (GEO), DESeq2 here", "nv": "NeuroVIP pipeline", "cl": "Classic analysis"}

# ---------------------------------------------------------------- paper parameters
PADJ_MAX = 0.05                     # Benjamini–Hochberg adjusted p <= 0.05
FC_MIN = 1.5                        # fold change >= 1.5 in either direction
LFC_MIN = math.log2(FC_MIN)         # = 0.585
MEAN_MIN = 5                        # "minimum combined mean of five reads"
KOBAS_P_MAX = 0.05                  # enrichment gene selection: p < 0.05 (unadjusted), |log2FC| > 0.585
PCA_TOP_GENES = 500                 # DESeq2 plotPCA default
HEATMAP_TOP = 25                    # Fig. 6c: top 25 up and 25 down, ranked by significance
GENESET_MIN = 5                     # ignore pathways with fewer annotated genes in the universe

# ---------------------------------------------------------------- what the paper reports
# Results text and the Fig. 6b volcano legend.
PAPER_REPORTED = {"n_up": 4382, "n_down": 4354, "n_unchanged": 12670}
# Fig. 6b PCA axis labels.
PAPER_PCA = {"pct_PC1": 39.2, "pct_PC2": 14.1}
# Fig. 6c gene labels, top to bottom. "40057" is printed in the figure: it is Sept9 after
# Excel date conversion (SEPT9 -> Sep-09 -> 1 Sep 2009 -> serial 40057).
PAPER_TOP25 = {
    "up": "Wars Gm8730 Slc41a2 Batf2 Tmem252 Gbp2 Zbp1 Pik3r5 Tgtp1 Igtp Samhd1 Pfkfb3 Gbp5 Ifit2 "
          "Trim30d Apol10b Cd274 Slc15a3 Gbp4 Slc10a6 Nlrc5 Rnd1 Rgs16 Ido1 Saa3",
    "down": "Palm Ptprb 40057 Fzd6 Cxcl12 Abcc6 Slc16a4 Notum Gtf2i Afap1l1 Gm684 Pcp4l1 Slc6a6 Golm1 "
            "Bche Hey1 Crispld1 Gstm7 Hmcn1 Rhpn2 Tcf7 Lfng Prickle1 Adgre5 Jup",
}
PAPER_TOP25_ALIAS = {"40057": "Sept9"}

# Pathways on the x-axis of Fig. 6f (KEGG) and Fig. 6e (PANTHER), read from the published figure.
# 'text' = named in the Results text or marked with an arrow; 'figure' = only plotted.
PAPER_PATHWAYS = {
    "KEGG": {
        "mmu05200": "text", "mmu04668": "figure", "mmu05168": "figure", "mmu05166": "figure",
        "mmu05164": "figure", "mmu04380": "figure", "mmu04210": "figure", "mmu05152": "figure",
        "mmu04060": "figure", "mmu05205": "figure", "mmu04612": "figure", "mmu04145": "figure",
        "mmu01100": "figure", "mmu05169": "figure", "mmu04620": "figure", "mmu05145": "figure",
        "mmu05140": "figure", "mmu05162": "figure", "mmu05416": "figure", "mmu04933": "figure",
        "mmu05142": "figure", "mmu04064": "text", "mmu04010": "figure", "mmu04151": "text",
        "mmu04630": "text", "mmu05134": "figure", "mmu04015": "figure", "mmu05133": "figure",
        "mmu04062": "figure", "mmu04621": "figure", "mmu05161": "figure", "mmu05160": "figure",
        "mmu04066": "text", "mmu04014": "text", "mmu04662": "figure", "mmu01524": "figure",
        "mmu05323": "figure", "mmu05332": "figure", "mmu04510": "figure", "mmu04514": "figure",
        "mmu05150": "figure", "mmu05144": "figure", "mmu05203": "figure", "mmu05330": "figure",
        "mmu05321": "figure", "mmu04940": "figure", "mmu05202": "figure", "mmu04390": "figure",
        "mmu05215": "figure", "mmu04670": "figure",
    },
    # Fig. 6e pathways left of the significance cut-off; points near the line are hard to read.
    "PANTHER": {
        "P00006": "figure", "P00005": "text", "P00031": "figure", "P00036": "figure",
        "P00054": "figure", "P00047": "text", "P00004": "figure", "P00019": "figure",
        "P00018": "figure", "P06664": "figure", "P00034": "figure", "P00056": "text",
        "P00046": "figure", "P00021": "figure", "P00059": "text",
    },
}
# Pathways central to the paper's hypothesis, always kept in the output.
FOCUS = {"mmu04066", "mmu04370", "mmu04064", "mmu04630", "mmu04115", "mmu04151", "mmu04014",
         "mmu04530", "mmu04670", "P00056", "P00047", "P00059", "P00030", "P00005"}

# ---------------------------------------------------------------- genes named in the paper
# claim: what the paper reports for infected vs healthy microvessels.
#   marker       : used as an endothelial marker, no expression claim
#   housekeeping : reference gene in the paper (qRT-PCR in cultured cells, or Western blot)
PAPER_GENES = [
    # symbol, label in paper, category, claim, where
    ("Hif1a", "HIF-1α", "HIF-1α / VEGF axis", "up", "Fig. 6d, Results"),
    ("Vegfa", "VEGF", "HIF-1α / VEGF axis", "up", "Fig. 6d, Results"),
    ("Angpt2", "ANGPT2", "Angiopoietin–Tie2", "up", "Fig. 6d, Results"),
    ("Tek", "Tie2", "Angiopoietin–Tie2", "down", "Fig. 6d"),
    ("Ptprb", "PTPRB (VE-PTP)", "Angiopoietin–Tie2", "down", "Fig. 6c ('dysregulated' in Results)"),
    ("Il1b", "IL-1β", "Inflammation", "up", "Fig. 6d, Results"),
    ("Tnf", "TNF-α", "Inflammation", "up", "Results, Discussion"),
    ("Snai1", "Snail", "EndMT", "up", "Discussion ('twofold upregulation')"),
    ("Cldn5", "Claudin-5", "BBB junctions", "down", "Fig. 6d, Results"),
    ("Ocln", "Occludin", "BBB junctions", "down", "Fig. 6d, Results"),
    ("Tjp1", "ZO-1", "BBB junctions", "down", "Fig. 6d, Results"),
    ("Cav1", "Caveolin-1", "BBB transcytosis", "down", "Discussion ('potentially compensatory decrease')"),
    ("Plvap", "PLVAP", "BBB transcytosis", "unchanged", "Discussion"),
    ("Pecam1", "PECAM-1 / CD31", "Endothelial markers", "down", "Fig. 6d"),
    ("Vwf", "vWF", "Endothelial markers", "unchanged", "Fig. 6d (n.s.)"),
    ("Cdh5", "VE-cadherin", "Endothelial markers", "marker", "Methods (IF marker)"),
    ("Rplp0", "RPLP0", "Reference genes", "housekeeping", "Methods (qRT-PCR, cultured cells)"),
    ("G6pdx", "G6PDX", "Reference genes", "housekeeping", "Methods (qRT-PCR, cultured cells)"),
    ("Actb", "β-actin", "Reference genes", "housekeeping", "Methods (Western blot loading control)"),
]

# ---------------------------------------------------------------- literature panels
# BBB-specific gene groups that are not KEGG/PANTHER/Reactome pathways. Curated from the
# literature; listed here so they can be checked.
LIT_PANELS = {
    "LIT_JUNC": ("BBB junctions", "Tight and adherens junction genes",
                 "Cldn5 Ocln Tjp1 Tjp2 Tjp3 Cldn1 Cldn3 Cldn12 Marveld2 Jam2 Jam3 F11r Esam Cdh5 "
                 "Ctnnb1 Ctnna1 Ctnnd1 Pecam1 Lsr Cgn Cgnl1"),
    "LIT_TRANS": ("BBB transporters & identity", "BBB endothelial transporters, Wnt/β-catenin identity genes",
                  "Slc2a1 Mfsd2a Abcb1a Abcg2 Slc7a5 Slc16a1 Slc3a2 Tfrc Lrp1 Slco1a4 Slco1c1 "
                  "Slc22a8 Slc38a5 Lef1 Foxf2 Zic3 Foxq1 Apcdd1 Axin2"),
    "LIT_CYTO": ("Transcytosis", "Caveolar transcytosis and its regulators", "Cav1 Cav2 Cavin1 Mfsd2a Plvap"),
    "LIT_ANGPT": ("Angiopoietin–Tie2 axis", "Angiopoietins and their receptors", "Angpt1 Angpt2 Angpt4 Tek Tie1 Ptprb"),
    "LIT_GLYC": ("HIF-1 glycolytic targets", "Classic HIF-1 target genes in glycolysis",
                 "Slc2a1 Slc2a3 Hk1 Hk2 Gpi1 Pfkl Pfkm Pfkp Pfkfb3 Aldoa Tpi1 Gapdh Pgk1 Pgam1 Eno1 "
                 "Pkm Ldha Pdk1 Slc16a3 Bnip3 Car9"),
    "LIT_ADH": ("Leukocyte adhesion", "Endothelial adhesion molecules and leukocyte ligands",
                "Sele Selp Icam1 Vcam1 Ackr1 Madcam1 Itgb2 Itgal Itgam Itga4 Sell Cd99 Pecam1"),
}

# ---------------------------------------------------------------- downloads
GEO_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122952/suppl/"
           "GSE122952_counts.matrix.norm_anno.txt.gz")
COUNTS_FILE = "GSE122952_counts.matrix.norm_anno.txt.gz"
GENESETS_FILE = "genesets.tsv"
URLS = {
    "kegg_link": "https://rest.kegg.jp/link/pathway/mmu",
    "kegg_list": "https://rest.kegg.jp/list/pathway/mmu",
    "kegg_genes": "https://rest.kegg.jp/list/mmu",
    "reactome": "https://reactome.org/download/current/Ensembl2Reactome.txt",
    "panther": "https://data.pantherdb.org/ftp/pathway/current_release/SequenceAssociationPathway3.6.8.txt",
    "mgi": "https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt",
}
# GEO column -> SRA run (GSM3489074–81 / SRR8244647–54, from ENA sample titles)
GEO_TO_RUN = {
    "Meningitis_1": "SRR8244647", "Meningitis_2": "SRR8244648", "Meningitis_3": "SRR8244649",
    "Meningitis_4": "SRR8244650", "Healthy_1": "SRR8244651", "Healthy_2": "SRR8244652",
    "Healthy_3": "SRR8244653", "Healthy_4": "SRR8244654",
}


def fetch(url: str, dest: str, offline: bool, keep=None) -> str:
    """Download url to dest once. keep(line) filters lines of large text files."""
    if os.path.exists(dest):
        return dest
    if offline:
        sys.exit(f"missing {dest} and --offline was given (source: {url})")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    print(f"  downloading {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "classical-analysis-script"})
    tmp = dest + ".part"
    with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as f:
        if keep is None:
            while chunk := r.read(1 << 20):
                f.write(chunk)
        else:
            for line in io.TextIOWrapper(r, encoding="utf-8", errors="replace"):
                if keep(line):
                    f.write(line.encode())
    os.replace(tmp, dest)
    return dest


# ================================================================ helpers
def bh(p) -> np.ndarray:
    """Benjamini–Hochberg adjusted p-values (same as R p.adjust(method='BH'))."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def call(lfc, padj, lfc_min=LFC_MIN):
    """up / down = passes paper thresholds; below_fc = padj ok but |FC| < 1.5; '' = not significant."""
    if pd.isna(padj) or pd.isna(lfc) or padj > PADJ_MAX:
        return ""
    if abs(lfc) < lfc_min:
        return "below_fc"
    return "up" if lfc > 0 else "down"


def read_csv(path):
    return pd.read_csv(path, keep_default_na=True)


# ================================================================ inputs
def load_counts(raw: str, offline: bool):
    """Paper count matrix (genes x 8 samples), gene annotation, sample table."""
    path = fetch(GEO_URL, os.path.join(raw, COUNTS_FILE), offline)
    d = pd.read_csv(path, sep="\t", low_memory=False, dtype={"Ensembl gene id": str})
    d = d.set_index("Ensembl gene id")
    cols = [c for c in d.columns if c in GEO_TO_RUN]
    if len(cols) != 8:
        sys.exit(f"{path}: expected 8 sample columns, found {cols}")
    counts = d[cols].apply(pd.to_numeric, errors="coerce").fillna(0).round().astype(int)
    annot = pd.DataFrame({"symbol": d["Ensembl gene"].fillna(pd.Series(d.index, index=d.index))})
    samples = pd.DataFrame({
        "sample": cols,
        "run": [GEO_TO_RUN[c] for c in cols],
        "condition": ["case" if c.startswith("Meningitis") else "control" for c in cols],
    })
    samples["label"] = [("Infected " if c.startswith("Men") else "Healthy ") + c.split("_")[1] for c in cols]
    return counts, annot, samples


def run_deseq2(counts: pd.DataFrame, samples: pd.DataFrame):
    """DESeq2 (pydeseq2) infected vs healthy on genes with mean >= 5 reads; also VST."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    keep = counts.mean(axis=1) >= MEAN_MIN
    m = counts[keep]
    meta = pd.DataFrame({"condition": np.where(samples.condition == "case", "infected", "healthy")},
                        index=samples["sample"].values)
    dds = DeseqDataSet(counts=m.T.loc[meta.index], metadata=meta, design="~condition", quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "infected", "healthy"], alpha=0.1, quiet=True)
    st.summary()
    res = st.results_df[["baseMean", "log2FoldChange", "pvalue", "padj"]].copy()
    dds.vst(use_design=False)
    vst = pd.DataFrame(np.asarray(dds.layers["vst_counts"]).T, index=m.index, columns=meta.index)
    return res, vst, int(keep.sum())


def regulated_in_all_four(counts: pd.DataFrame, samples: pd.DataFrame) -> pd.Series:
    """+1 / -1 if every infected sample is >= 1.5-fold above / below the healthy mean, else 0.

    This is how the paper's KOBAS gene selection ("commonly regulated in all four meningitis
    samples compared to all the healthy ones") is interpreted here.
    """
    inf = counts[samples[samples.condition == "case"]["sample"]]
    hea = counts[samples[samples.condition == "control"]["sample"]].mean(axis=1)
    lfc = np.log2((inf + 0.5).div(hea + 0.5, axis=0))
    up = (lfc >= LFC_MIN).all(axis=1)
    down = (lfc <= -LFC_MIN).all(axis=1)
    return pd.Series(np.where(up, 1, np.where(down, -1, 0)), index=counts.index)


def build_deg(data_dir: str, res: pd.DataFrame, annot: pd.DataFrame, all4: pd.Series) -> pd.DataFrame:
    """One row per gene: paper-count re-run (pa_), NeuroVIP (nv_) and classic (cl_) DESeq2."""
    pa = pd.DataFrame({
        "gene_id": res.index, "symbol": annot.symbol.reindex(res.index).values,
        "base_mean": res.baseMean.values, "pa_log2fc": res.log2FoldChange.values,
        "pa_pvalue": res.pvalue.values, "pa_padj": res.padj.values,
        "pa_all4": all4.reindex(res.index).fillna(0).astype(int).values,
    })
    cmp = read_csv(os.path.join(data_dir, "ai_vs_classic.csv"))
    cmp = cmp[cmp.disease == DISEASE]
    pipe = pd.DataFrame({"gene_id": cmp.gene_id, "symbol_pipe": cmp.symbol,
                         "nv_log2fc": cmp.ai_log2fc, "nv_padj": cmp.ai_padj,
                         "cl_log2fc": cmp.mk_log2fc, "cl_padj": cmp.mk_padj})
    deg = pa.merge(pipe, on="gene_id", how="outer")
    deg["in_paper_matrix"] = deg.pa_log2fc.notna() | deg.gene_id.isin(annot.index)
    deg["symbol"] = deg.symbol.fillna(deg.pop("symbol_pipe")).fillna(deg.gene_id)
    for s in SOURCES:
        deg[f"{s}_call"] = [call(a, b) for a, b in zip(deg[f"{s}_log2fc"], deg[f"{s}_padj"])]
    return deg.reset_index(drop=True)


# ================================================================ analyses
def summary(deg: pd.DataFrame, counts: pd.DataFrame, samples, n_tested: dict) -> pd.DataFrame:
    rows = []
    for s, name in SOURCES.items():
        fc, p = deg[f"{s}_log2fc"], deg[f"{s}_padj"]
        sig = p.notna() & (p <= PADJ_MAX)
        big = sig & (fc.abs() >= LFC_MIN)
        rows.append(dict(analysis=name, source=s, tkey="padj", threshold="padj ≤ 0.05",
                         n_up=int((sig & (fc > 0)).sum()), n_down=int((sig & (fc < 0)).sum()),
                         genes_tested=n_tested[s]))
        rows.append(dict(analysis=name, source=s, tkey="strict", threshold="padj ≤ 0.05 and |FC| ≥ 1.5",
                         n_up=int((big & (fc > 0)).sum()), n_down=int((big & (fc < 0)).sum()),
                         genes_tested=n_tested[s]))
    rows.append(dict(analysis="Devraj et al. 2020 (reported)", source="paper", tkey="reported",
                     threshold="as reported (Results, Fig. 6b)", n_up=PAPER_REPORTED["n_up"],
                     n_down=PAPER_REPORTED["n_down"], genes_tested=sum(PAPER_REPORTED.values())))
    out = pd.DataFrame(rows)
    out["n_total"] = out.n_up + out.n_down
    out["pct_up"] = (100 * out.n_up / out.genes_tested).round(1)
    out["pct_down"] = (100 * out.n_down / out.genes_tested).round(1)
    paper = out[out.source == "paper"].iloc[0]
    out["diff_up_pct"] = (100 * (out.n_up - paper.n_up) / paper.n_up).round(1)
    out["diff_down_pct"] = (100 * (out.n_down - paper.n_down) / paper.n_down).round(1)
    return out


def sensitivity(counts, samples) -> pd.DataFrame:
    """Same DESeq2 run without the >= 5 reads filter, to show what the filter changes."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    meta = pd.DataFrame({"condition": np.where(samples.condition == "case", "infected", "healthy")},
                        index=samples["sample"].values)
    m = counts[counts.sum(axis=1) > 0]
    dds = DeseqDataSet(counts=m.T.loc[meta.index], metadata=meta, design="~condition", quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "infected", "healthy"], alpha=0.1, quiet=True)
    st.summary()
    r = st.results_df
    sig = r.padj <= PADJ_MAX
    return dict(n_up=int((sig & (r.log2FoldChange > 0)).sum()), n_down=int((sig & (r.log2FoldChange < 0)).sum()),
                genes_tested=len(m))


def paper_genes(deg: pd.DataFrame) -> pd.DataFrame:
    by_sym = deg.drop_duplicates("symbol").set_index("symbol")
    rows = []
    for sym, label, cat, claim, where in PAPER_GENES:
        r = by_sym.loc[sym] if sym in by_sym.index else None
        row = dict(symbol=sym, paper_label=label, category=cat, paper_claim=claim, paper_source=where,
                   gene_id=r.gene_id if r is not None else "")
        for s in SOURCES:
            fc = r[f"{s}_log2fc"] if r is not None else np.nan
            p = r[f"{s}_padj"] if r is not None else np.nan
            sig = pd.notna(p) and p <= PADJ_MAX
            obs = ("up" if fc > 0 else "down") if sig else "n.s."
            if claim in ("up", "down"):
                agree = "agrees" if obs == claim else "differs"
            elif claim == "unchanged":
                agree = "agrees" if obs == "n.s." else "differs"
            elif claim == "housekeeping":
                agree = "stable" if obs == "n.s." else "changed"
            else:
                agree = "—"
            row.update({f"{s}_log2fc": fc, f"{s}_padj": p, f"{s}_observed": obs, f"{s}_agreement": agree,
                        f"{s}_passes_fc": bool(sig and abs(fc) >= LFC_MIN)})
        rows.append(row)
    return pd.DataFrame(rows)


def read_counts(counts, annot, samples) -> pd.DataFrame:
    """Normalised read counts of the named genes per sample (paper Fig. 6d)."""
    sym2id = {s: g for g, s in annot.symbol.items()}
    rows = []
    for sym, label, cat, *_ in PAPER_GENES:
        g = sym2id.get(sym)
        if g is None or g not in counts.index:
            continue
        for s in samples.itertuples():
            rows.append(dict(symbol=sym, paper_label=label, category=cat, gene_id=g, sample=s.sample,
                             run=s.run, label=s.label, condition=s.condition, count=int(counts.at[g, s.sample])))
    return pd.DataFrame(rows)


def top_genes(deg: pd.DataFrame) -> pd.DataFrame:
    """Top 25 up and top 25 down passing paper thresholds, per run, ranked two ways.

    Ranking by significance uses padj, then the Wald p-value / |log2FC| to break ties
    (padj underflows to 0 for the strongest genes).
    """
    out = []
    for s in SOURCES:
        ok = deg[deg[f"{s}_call"].isin(["up", "down"])].copy()
        ok = ok[~ok.symbol.astype(str).str.startswith("ENSMUSG")]
        for rank_by in ("padj", "log2fc"):
            for direction in ("up", "down"):
                d = ok[ok[f"{s}_call"] == direction].copy()
                d["_abs"] = d[f"{s}_log2fc"].abs()
                if rank_by == "padj":
                    keys = [f"{s}_padj"] + (["pa_pvalue"] if s == "pa" else []) + ["_abs"]
                    d = d.sort_values(keys, ascending=[True] * (len(keys) - 1) + [False])
                else:
                    d = d.sort_values("_abs", ascending=False)
                d = d.head(HEATMAP_TOP)
                for i, r in enumerate(d.itertuples(), 1):
                    out.append(dict(source=s, rank_by=rank_by, direction=direction, rank=i, symbol=r.symbol,
                                    gene_id=r.gene_id, log2fc=getattr(r, f"{s}_log2fc"),
                                    padj=getattr(r, f"{s}_padj")))
    top = pd.DataFrame(out)
    paper = {g: d for d, lst in PAPER_TOP25.items() for g in (PAPER_TOP25_ALIAS.get(x, x) for x in lst.split())}
    top["in_paper_top25"] = [paper.get(sym) == d for sym, d in zip(top.symbol, top.direction)]
    return top


def paper_top(deg: pd.DataFrame) -> pd.DataFrame:
    """Each gene of paper Fig. 6c and its rank by significance in each run here."""
    rows = []
    for direction, lst in PAPER_TOP25.items():
        for order, label in enumerate(lst.split(), 1):
            sym = PAPER_TOP25_ALIAS.get(label, label)
            row = dict(direction=direction, paper_order=order, paper_label=label, symbol=sym)
            hit = deg[deg.symbol == sym]
            for s in SOURCES:
                d = deg[deg[f"{s}_call"] == direction].copy()
                d = d[~d.symbol.astype(str).str.startswith("ENSMUSG")]
                keys = [f"{s}_padj"] + (["pa_pvalue"] if s == "pa" else [])
                d["_abs"] = -d[f"{s}_log2fc"].abs()
                d = d.sort_values(keys + ["_abs"]).reset_index(drop=True)
                pos = d.index[d.symbol == sym]
                row[f"{s}_rank"] = int(pos[0]) + 1 if len(pos) else np.nan
                row[f"{s}_log2fc"] = hit[f"{s}_log2fc"].iloc[0] if len(hit) else np.nan
                row[f"{s}_padj"] = hit[f"{s}_padj"].iloc[0] if len(hit) else np.nan
                row[f"{s}_call"] = hit[f"{s}_call"].iloc[0] if len(hit) else "absent"
            rows.append(row)
    return pd.DataFrame(rows)


def sample_level(vst: pd.DataFrame, samples: pd.DataFrame, top: pd.DataFrame):
    """PCA, sample clustering and per-sample Z-scores, all from the VST matrix."""
    sm = samples.set_index("sample")
    pca_rows = []
    # top500: DESeq2 plotPCA (centred, top 500 variable genes).
    # scaled : all genes, each gene scaled to unit variance (FactoMineR PCA default; the paper's
    #          axis labels "Dim 1 / Dim 2" suggest FactoMineR, but its settings are not stated).
    for method, n, scale in (("top500", PCA_TOP_GENES, False), ("scaled", len(vst), True)):
        sel = vst.loc[vst.var(axis=1).sort_values(ascending=False).index[:n]]
        sel = sel[sel.var(axis=1) > 0]
        x = sel.T - sel.T.mean(axis=0)
        if scale:
            x = x / sel.T.std(axis=0, ddof=1)
        x = x.values
        u, s, _ = np.linalg.svd(x, full_matrices=False)
        pcs = u * s
        pct = 100 * s ** 2 / np.sum(s ** 2)
        # orient PC1 so that infected samples have positive scores (sign of a PC is arbitrary)
        case = (sm.loc[sel.columns, "condition"] == "case").values
        if pcs[case, 0].mean() < 0:
            pcs[:, 0] *= -1
        for i, smp in enumerate(sel.columns):
            pca_rows.append(dict(method=method, n_genes=len(sel), sample=smp, run=sm.at[smp, "run"],
                                 label=sm.at[smp, "label"], condition=sm.at[smp, "condition"],
                                 PC1=round(pcs[i, 0], 4), PC2=round(pcs[i, 1], 4),
                                 pct_PC1=round(pct[0], 1), pct_PC2=round(pct[1], 1),
                                 paper_pct_PC1=PAPER_PCA["pct_PC1"], paper_pct_PC2=PAPER_PCA["pct_PC2"]))
    pca = pd.DataFrame(pca_rows)

    # hierarchical clustering of samples: Euclidean distance on VST, complete linkage (R hclust default)
    dvec = pdist(vst.T.values, metric="euclidean")
    link = hierarchy.linkage(dvec, method="complete")
    dendro = hierarchy.dendrogram(link, no_plot=True, labels=list(vst.columns))
    order = dendro["ivl"]
    dist = squareform(dvec)
    cols = list(vst.columns)
    dist_rows = [dict(sample_a=a, sample_b=b, dist=round(dist[i, j], 3),
                      order_a=order.index(a), order_b=order.index(b),
                      label_a=sm.at[a, "label"], label_b=sm.at[b, "label"],
                      condition_a=sm.at[a, "condition"], condition_b=sm.at[b, "condition"])
                 for i, a in enumerate(cols) for j, b in enumerate(cols)]
    seg = []
    for k, (xs, ys) in enumerate(zip(dendro["icoord"], dendro["dcoord"])):
        # scipy places leaves at x = 5, 15, 25, ...; convert to leaf index 0..n-1
        xs = [(v - 5) / 10 for v in xs]
        for i in range(3):
            seg.append(dict(link=k, x1=xs[i], y1=round(ys[i], 3), x2=xs[i + 1], y2=round(ys[i + 1], 3)))
    seg = pd.DataFrame(seg)

    # per-sample Z-scores for every gene in any top-25 list
    ids = top.gene_id.unique()
    hm = []
    for g in ids:
        if g not in vst.index:
            continue
        v = vst.loc[g].values
        z = (v - v.mean()) / (v.std(ddof=1) or 1)
        for smp, vv, zz in zip(vst.columns, v, z):
            hm.append(dict(gene_id=g, sample=smp, z=round(zz, 3), vst=round(vv, 3)))
    hm = pd.DataFrame(hm, columns=["gene_id", "sample", "z", "vst"])
    return pca, pd.DataFrame(dist_rows), seg, hm


# ================================================================ gene sets
def build_genesets(raw: str, offline: bool) -> pd.DataFrame:
    """Official KEGG, PANTHER and Reactome mouse pathways as (db, pathway_id, pathway, gene_id)."""
    out_path = os.path.join(raw, GENESETS_FILE)
    if os.path.exists(out_path):
        return pd.read_csv(out_path, sep="\t", dtype=str)
    src = os.path.join(raw, "genesets_src")
    f = {k: fetch(u, os.path.join(src, k + ".txt"), offline,
                  keep=(lambda l: "ENSMUSG" in l) if k == "reactome" else
                  (lambda l: "\tMOUSE|" in l) if k == "panther" else None)
         for k, u in URLS.items()}

    mgi = pd.read_csv(f["mgi"], sep="\t", header=None, usecols=[0, 1, 5], names=["mgi", "symbol", "gene_id"], dtype=str)
    sym2ens = dict(zip(mgi.symbol, mgi.gene_id))
    mgi2ens = dict(zip(mgi.mgi, mgi.gene_id))

    # KEGG: entrez -> symbol (first name in 'list/mmu') -> Ensembl via MGI
    kg = pd.read_csv(f["kegg_genes"], sep="\t", header=None, dtype=str, names=["kid", "type", "pos", "desc"])
    kg["symbol"] = kg.desc.str.split(";").str[0].str.split(",").str[0].str.strip()
    k2ens = {k: sym2ens.get(s) for k, s in zip(kg.kid, kg.symbol)}
    link = pd.read_csv(f["kegg_link"], sep="\t", header=None, names=["kid", "pid"], dtype=str)
    names = pd.read_csv(f["kegg_list"], sep="\t", header=None, names=["pid", "name"], dtype=str)
    names["name"] = names.name.str.replace(r" - Mus musculus \(house mouse\)$", "", regex=True)
    pname = dict(zip("path:" + names.pid, names.name))
    kegg = pd.DataFrame({"db": "KEGG", "pathway_id": link.pid.str.replace("path:", ""),
                         "pathway": link.pid.map(pname), "gene_id": link.kid.map(k2ens)})

    # PANTHER: MGI accession -> Ensembl
    pth = pd.read_csv(f["panther"], sep="\t", header=None, dtype=str, usecols=[0, 1, 4],
                      names=["pathway_id", "pathway", "ids"])
    pth["mgi"] = "MGI:" + pth.ids.str.extract(r"MGI=MGI[=:](\d+)")[0]
    panther = pd.DataFrame({"db": "PANTHER", "pathway_id": pth.pathway_id, "pathway": pth.pathway,
                            "gene_id": pth.mgi.map(mgi2ens)})

    # Reactome: Ensembl gene IDs directly (lowest-level pathways)
    rc = pd.read_csv(f["reactome"], sep="\t", header=None, dtype=str, usecols=[0, 1, 3, 5],
                     names=["gene_id", "pathway_id", "pathway", "species"])
    rc = rc[rc.species == "Mus musculus"]
    reactome = pd.DataFrame({"db": "Reactome", "pathway_id": rc.pathway_id, "pathway": rc.pathway,
                             "gene_id": rc.gene_id})

    gs = pd.concat([kegg, panther, reactome]).dropna().drop_duplicates()
    gs.to_csv(out_path, sep="\t", index=False)
    return gs


def add_literature(gs: pd.DataFrame, sym2id: dict) -> pd.DataFrame:
    rows = [dict(db="Literature", pathway_id=pid, pathway=name, gene_id=sym2id[g])
            for pid, (name, _desc, genes) in LIT_PANELS.items() for g in genes.split() if g in sym2id]
    return pd.concat([gs, pd.DataFrame(rows)], ignore_index=True)


def enrichment(deg: pd.DataFrame, gs: pd.DataFrame, universe: set) -> tuple[pd.DataFrame, pd.DataFrame]:
    """KOBAS-style over-representation (hypergeometric, BH per database) + direction (binomial).

    Universe: genes tested in the paper-count DESeq2 run (mean >= 5 reads). Selected genes:
      pa: p < 0.05, |log2FC| >= 0.585 and regulated in all four infected samples (paper's selection)
      nv, cl: padj <= 0.05 and |FC| >= 1.5 (no per-sample data for these runs), restricted to the universe
    """
    gs = gs[gs.gene_id.isin(universe)]
    sizes = gs.groupby(["db", "pathway_id"]).gene_id.nunique()
    keep = sizes[sizes >= GENESET_MIN].index
    gs = gs.set_index(["db", "pathway_id"]).loc[keep].reset_index()
    names = gs.drop_duplicates(["db", "pathway_id"]).set_index(["db", "pathway_id"]).pathway
    members = gs.groupby(["db", "pathway_id"]).gene_id.apply(set)
    N = len(universe)
    paper_ref = {pid: v for db in PAPER_PATHWAYS.values() for pid, v in db.items()}
    rows = []
    for s in SOURCES:
        d = deg[deg.gene_id.isin(universe)]
        if s == "pa":
            sel = d[(d.pa_pvalue < KOBAS_P_MAX) & (d.pa_log2fc.abs() >= LFC_MIN) &
                    (np.sign(d.pa_log2fc) == d.pa_all4)]
        else:
            sel = d[d[f"{s}_call"].isin(["up", "down"])]
        up_ids = set(sel.gene_id[sel[f"{s}_log2fc"] > 0])
        hit_ids = set(sel.gene_id)
        n = len(hit_ids)
        p_up = len(up_ids) / n if n else 0.5
        for (db, pid), genes in members.items():
            K = len(genes)
            hits = genes & hit_ids
            k = len(hits)
            up = len(hits & up_ids)
            rows.append(dict(source=s, db=db, pathway_id=pid, pathway=names[(db, pid)],
                             n_genes=K, n_sig=k, n_up=up, n_down=k - up, n_selected=n, universe=N,
                             expected=round(K * n / N, 2),
                             fold_enrichment=round((k / K) / (n / N), 3) if n else np.nan,
                             p_enrich=stats.hypergeom.sf(k - 1, N, K, n) if k else 1.0,
                             p_direction=stats.binomtest(up, k, p_up).pvalue if k else 1.0))
    out = pd.DataFrame(rows)
    out["q_enrich"] = out.groupby(["source", "db"]).p_enrich.transform(bh)
    out["q_direction"] = out.groupby(["source", "db"]).p_direction.transform(bh)
    out["direction"] = np.where(out.n_sig == 0, "—", np.where(out.n_up > out.n_down, "mostly up",
                                np.where(out.n_down > out.n_up, "mostly down", "balanced")))
    out["rank"] = out.groupby(["source", "db"]).p_enrich.rank(method="min").astype(int)
    out["paper"] = out.pathway_id.map(paper_ref).fillna("")
    out["focus"] = out.pathway_id.isin(FOCUS) | (out.db == "Literature")
    # keep the file small: significant somewhere, in the paper, or a focus pathway
    sig_any = out[out.q_enrich <= 0.05].pathway_id.unique()
    out = out[out.pathway_id.isin(sig_any) | (out.paper != "") | out.focus]
    out = out.sort_values(["source", "db", "p_enrich"]).reset_index(drop=True)
    mem = gs[gs.pathway_id.isin(out.pathway_id.unique())][["db", "pathway_id", "gene_id"]].drop_duplicates()
    return out, mem


# ================================================================ methods table
def methods_table(summ, sens, n_tested, pca, top, ptop, paths) -> pd.DataFrame:
    pa = summ[(summ.source == "pa") & (summ.tkey == "padj")].iloc[0]
    pas = summ[(summ.source == "pa") & (summ.tkey == "strict")].iloc[0]
    p5 = pca[pca.method == "top500"].iloc[0]
    pall = pca[pca.method == "scaled"].iloc[0]
    ov = {d: int(ptop[(ptop.direction == d) & (ptop.pa_rank <= HEATMAP_TOP)].shape[0]) for d in ("up", "down")}
    pk = paths[(paths.source == "pa") & (paths.db == "KEGG") & (paths.paper != "")]
    pp = paths[(paths.source == "pa") & (paths.db == "PANTHER") & (paths.paper != "")]
    n_sel = int(paths[paths.source == "pa"].n_selected.iloc[0])
    rows = [
        ("Read QC, trimming, alignment, counting",
         "FastQC; Reaper (Q20, 30–150 nt); STAR 2.4.0a to mm10; featureCounts 1.4.5",
         "reused",
         "The paper's own count matrix from GEO (GSE122952) is used, so its alignment and counting are kept. "
         "NeuroVIP's runs re-aligned the same 8 SRA runs to GRCm39 / Ensembl 115 with STAR 2.7.10a."),
        ("Differential expression",
         "DESeq2; BH padj ≤ 0.05; |FC| ≥ 1.5; combined mean ≥ 5 reads",
         "applied",
         f"DESeq2 (pydeseq2) on the GEO matrix after the ≥ 5-read filter ({n_tested['pa']:,} genes): "
         f"{int(pa.n_up):,} up / {int(pa.n_down):,} down at padj ≤ 0.05, {int(pas.n_up):,} / {int(pas.n_down):,} "
         f"with the 1.5-fold cut-off. Without the read filter: {sens['n_up']:,} / {sens['n_down']:,} of "
         f"{sens['genes_tested']:,}. The GEO values are already DESeq2-normalised and rounded, so size factors are ≈ 1."),
        ("PCA of samples (Fig. 6b)",
         "PCA of meningitis vs healthy microvessels",
         "applied",
         f"VST, top {PCA_TOP_GENES} most variable genes (DESeq2 plotPCA): PC1 {p5.pct_PC1} %, PC2 {p5.pct_PC2} %. "
         f"All {int(pall.n_genes):,} genes scaled to unit variance: PC1 {pall.pct_PC1} %, PC2 {pall.pct_PC2} %. "
         f"Paper: Dim 1 {PAPER_PCA['pct_PC1']} %, Dim 2 {PAPER_PCA['pct_PC2']} % (exact PCA settings not stated)."),
        ("Hierarchical clustering (Supp. Fig. 8)",
         "Complete hierarchical clustering of the samples",
         "applied",
         "Euclidean distances between samples on VST values, complete linkage (R hclust default). "
         "The supplementary figure was not available for comparison."),
        ("Top 25 up / down genes (Fig. 6c)",
         "Row Z-scores of the 25 most significant up- and down-regulated genes",
         "applied",
         f"Per-sample VST Z-scores. Of the paper's 25 up and 25 down genes, {ov['up']} and {ov['down']} are in the "
         "top 25 here (paper-count run). The paper's label '40057' is Sept9 after Excel date conversion."),
        ("Read counts of named genes (Fig. 6d)",
         "Normalised reads of HIF-1α, VEGF, IL-1β, ANGPT2, Tie2, vWF, PECAM, ZO-1, claudin-5, occludin",
         "applied",
         "Per-sample normalised counts from the GEO matrix, with DESeq2 results from all three runs."),
        ("Pathway enrichment (Fig. 6e–f, Supp. Fig. 9)",
         "KOBAS 2.0 (KEGG, PANTHER, Reactome) on genes regulated in all four infected samples (p < 0.05, |log2FC| > 0.585)",
         "applied",
         f"Hypergeometric test with BH correction per database (as KOBAS) on current official KEGG, PANTHER 3.6.8 and "
         f"Reactome mouse gene sets; universe = {n_tested['pa']:,} tested genes; {n_sel:,} selected genes. "
         f"'Regulated in all four' is read as: every infected sample ≥ 1.5-fold from the healthy mean. "
         f"{int((pk.q_enrich <= 0.05).sum())} of {len(pk)} KEGG and {int((pp.q_enrich <= 0.05).sum())} of {len(pp)} "
         "PANTHER pathways from the paper's figure reach q ≤ 0.05. Database versions differ from the paper's (2019), "
         "so membership and rankings shift."),
        ("DEG counts vs paper",
         f"{PAPER_REPORTED['n_up']:,} up / {PAPER_REPORTED['n_down']:,} down reported",
         "compared",
         f"Paper-count run at padj ≤ 0.05: {int(pa.n_up):,} up ({pa.diff_up_pct:+.1f} %) / {int(pa.n_down):,} down "
         f"({pa.diff_down_pct:+.1f} %). With the stated 1.5-fold cut-off: {int(pas.n_up):,} / {int(pas.n_down):,}, so the "
         "reported counts appear to use padj only."),
    ]
    return pd.DataFrame(rows, columns=["step", "paper_method", "status", "detail"])


# ================================================================ main
def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", default=os.path.join(here, "..", "data"))
    ap.add_argument("--raw", default=None, help="cache for downloads (default: <data>/raw)")
    ap.add_argument("--offline", action="store_true", help="never download; fail if a cached file is missing")
    a = ap.parse_args(argv)
    d = os.path.abspath(a.data)
    raw = os.path.abspath(a.raw or os.path.join(d, "raw"))

    studies = read_csv(os.path.join(d, "studies.csv")).set_index("disease")
    n_pipe = int(studies.loc[DISEASE, "genes_tested"])

    counts, annot, samples = load_counts(raw, a.offline)
    print(f"Paper count matrix: {counts.shape[0]:,} genes x {counts.shape[1]} samples")
    res, vst, n_pa = run_deseq2(counts, samples)
    sens = sensitivity(counts, samples)
    all4 = regulated_in_all_four(counts, samples)
    deg = build_deg(d, res, annot, all4)
    n_tested = {"pa": n_pa, "nv": n_pipe, "cl": n_pipe}

    summ = summary(deg, counts, samples, n_tested)
    genes = paper_genes(deg)
    rc = read_counts(counts, annot, samples)
    top = top_genes(deg)
    ptop = paper_top(deg)
    pca, dist, dendro, hm = sample_level(vst, samples, top)

    sym2id = {s: g for g, s in annot.symbol.reindex(res.index).items()}
    gs = add_literature(build_genesets(raw, a.offline), sym2id)
    paths, members = enrichment(deg, gs, set(res.index))
    meth = methods_table(summ, sens, n_tested, pca, top, ptop, paths)

    def save(df, name):
        df.to_csv(os.path.join(d, name), index=False, float_format="%.6g")
        print(f"  wrote {name:34s} {len(df):6d} rows")

    deg_out = deg[deg.pa_padj.notna() | deg.nv_padj.notna() | deg.cl_padj.notna()]
    save(meth, "classical_methods.csv")
    save(summ, "classical_summary.csv")
    save(deg_out[["gene_id", "symbol", "base_mean", "pa_log2fc", "pa_padj", "pa_all4", "nv_log2fc", "nv_padj",
                  "cl_log2fc", "cl_padj", "pa_call", "nv_call", "cl_call"]], "classical_deg.csv")
    save(genes, "classical_genes.csv")
    save(rc, "classical_counts.csv")
    save(top, "classical_top.csv")
    save(ptop, "classical_paper_top.csv")
    save(pca, "classical_pca.csv")
    save(dist, "classical_sample_dist.csv")
    save(dendro, "classical_dendrogram.csv")
    save(hm, "classical_sample_heatmap.csv")
    save(paths, "classical_pathways.csv")
    save(members, "classical_pathway_genes.csv")


if __name__ == "__main__":
    main()
