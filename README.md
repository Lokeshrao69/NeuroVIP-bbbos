# NeuroVIP

**NeuroVIP is an interactive portal for blood–brain barrier (BBB) microvessel transcriptomes in mouse models of ischemic stroke, bacterial meningitis and glioblastoma.**

🔗 **Live site:** https://neurovip-bbbos.vercel.app

---

## 1. What this is

NeuroVIP presents the results of an RNA-seq analysis pipeline built for this project. It lets you browse them in a web browser without running any code. It is designed around one question:

> **Do BBB-disrupting diseases of very different origin (ischemic, infectious and neoplastic) share a common transcriptional response in the brain microvasculature?**

All three datasets are RNA-seq of **brain microvessels** from mouse disease models, each compared with its own controls. The raw reads were processed with one uniform pipeline, so the results are directly comparable across diseases.

### Headline findings

| | Stroke | Meningitis | Glioblastoma |
|---|---:|---:|---:|
| Significant genes (padj < 0.05) | 5,227 | 9,291 | 8,687 |
| Up-regulated | 2,804 | 4,725 | 4,468 |
| Down-regulated | 2,423 | 4,566 | 4,219 |

- **1,647 genes** are significantly changed in **all three** diseases.
- Of these, **1,322 (80 %)** change in the **same direction** everywhere: 727 up and 595 down.
- The strongest shared genes form an **endothelial activation and inflammation programme**:
  - adhesion molecules: *Selp*, *Sele*, *Icam1*
  - *Ackr1*
  - acute-phase and antimicrobial genes: *Lcn2*, *Saa3*
  - inflammatory signalling: *Il6*, *Il1b*
  - matrix remodelling: *Mmp3*
- This pipeline agrees almost perfectly with the lab's independent classic analysis. The log2 fold-change correlation is r ≈ 0.9999 for genes found by both.

---

## 2. Data

| Disease | Model | Tissue | Cases / controls | Source |
|---|---|---|---|---|
| Ischemic stroke | tMCAO, 24 h, ipsilateral hemisphere, vs sham | Brain microvessels | 3 / 3 | Kestner et al., *Front. Neurosci.* 2020 — GEO **GSE131193** |
| Bacterial meningitis | Haematogenous *S. pneumoniae* infection vs healthy | Brain microvessels | 4 / 4 | Devraj G. et al., *Acta Neuropathol.* 2020 — GEO **GSE122952** |
| Glioblastoma | Tumour-bearing region vs sham | Brain microvessels | 3 / 3 | *Source study to be added* |

20 single-end Illumina libraries in total (NCBI SRA run accessions are listed on the Overview tab).

## 3. How the results were produced

```
Raw FASTQ (SRA)
  → FastQC 0.11.9                    read quality control
  → Trim Galore 0.6.7 / Cutadapt 5.2 adapter + quality trimming (Q20, min length 30 bp)
  → STAR 2.7.10a                     splice-aware alignment to mouse GRCm39 (Ensembl 115), unique reads only
  → featureCounts (Subread 2.0.3)    reads per gene (exons, MAPQ ≥ 60)
  → DESeq2 1.50.2                    case vs control within each study; Benjamini–Hochberg padj < 0.05
  → cross-disease comparison, machine learning (scikit-learn)
```

The pipeline ran in Google Colab with data on Google Drive. The same datasets were also analysed with the lab's existing ("classic") workflow, and NeuroVIP compares the two.

---

## 4. How to use NeuroVIP

The site has seven pages (tabs) and a gene search box in the top bar. Type a gene symbol there from any page to jump to Gene search; suggestions appear as you type. The pages along the top. Charts show details when you hover over them, and every table can be searched.

### Home
- **Tiles:** the number of significant genes (DEGs) per disease, split into up and down, plus the number shared by all three.
- **Datasets:** the model, tissue, source paper and GEO accession for each disease.
- **Sample quality:** reads, adapter content, alignment rate and the percentage of reads assigned to genes, for every sample.
  - Values highlighted in **amber** are low. The glioblastoma libraries assign only 47–67 % of reads to genes, against about 89 % for the others, which suggests a different library preparation.

### Differential expression
1. Choose a disease: **Stroke / Meningitis / Glioblastoma**.
2. Choose a fold-change threshold: **|log2FC| ≥ 0, 1 or 2**.
3. Read the **volcano plot**:
   - x-axis: log2 fold change (right = higher in disease)
   - y-axis: −log10 adjusted p-value (higher = more significant)
   - colours: red = up, blue = down; grey points are significant but below the fold-change threshold
   - the 12 strongest genes are labelled
4. The **gene table** below lists genes by significance. Type in the search box to find a gene, and click a row to open it in *Gene search*.

> The volcano plots show only significant genes, because the exported tables contain only those. The usual grey cloud of unchanged genes near zero is therefore absent.

### Gene search
Type a gene symbol (e.g. `Selp`, `Icam1`, `Cldn5`, `Spp1`) and press Enter, or click one of the example buttons.
- **Bars:** the log2 fold change in each disease, from this pipeline (blue) and the classic analysis (purple).
- **"n.s.":** not significant in that disease.
- **Table:** the exact values and adjusted p-values.

### Shared genes
- **Overlap chart (UpSet plot):** each bar is one exact combination of diseases, and the dots underneath show which ones. For example, the highlighted bar is genes significant in all three (1,647).
- **Direction of the shared genes:** how many go up everywhere, down everywhere, or in mixed directions.
- **Heatmap:** the 40 shared genes that are strongly changed in *every* disease (ranked by their smallest |log2FC| across the three). Red = up, blue = down.
- **All shared genes:** a searchable table of all 1,647. Click a row to open it in *Gene search*.

### Pipeline validation
Compares this pipeline with the lab's classic analysis of the same data.
- **Stacked bars:** genes found significant by both, by this pipeline only, or by the classic analysis only.
- **Tiles and scatter plot (pick a disease):** points lying on the diagonal mean both pipelines give the same fold change.

| | Both | This pipeline only | Classic only |
|---|---:|---:|---:|
| Stroke | 4,383 | 844 | 14 |
| Meningitis | 8,553 | 738 | 4 |
| Glioblastoma | 7,637 | 1,050 | 0 |

Fold changes are essentially identical. This pipeline calls somewhat more genes significant, probably because of differences in filtering or counting settings.

### Help
Definitions of all terms (log2 fold change, adjusted p-value, up/down-regulated), a summary of the processing pipeline, limitations and how to cite NeuroVIP.

### Classical analysis (meningitis)
Reproduces the RNA-seq analysis of Devraj et al. 2020 (*Acta Neuropathol.* 140:183–208) on the same experiment (GSE122952, 4 infected vs 4 healthy mice). The study's own count matrix from GEO is re-analysed with DESeq2 and the paper's thresholds (padj ≤ 0.05, FC ≥ 1.5, mean ≥ 5 reads) and compared with NeuroVIP's two DESeq2 runs.
- **DESeq2 run** switch (paper counts re-run, NeuroVIP pipeline, classic analysis): every chart and table on the tab follows it.
- **Methods table:** each paper method and how it was reproduced, with text generated from the results.
- **DEG counts vs paper:** counts with and without the 1.5-fold cut-off next to the paper's 4,382 up / 4,354 down.
- **Volcano plot** (all genes for the paper-count run), highlighting the paper's named genes, its top-25 genes or a selected pathway.
- **Genes named in the paper** and **read counts per sample** (Fig. 6d): HIF-1α, VEGF, ANGPT2/Tie2/PTPRB, IL-1β, TNF-α, Snail, claudin-5, occludin, ZO-1, caveolin-1, PLVAP, PECAM-1, vWF, VE-cadherin and the reference genes.
- **Top 25 up / 25 down** (Fig. 6c) as per-sample Z-scores, with a table ranking the paper's own top-25 genes in each run.
- **PCA** (Fig. 6b) and **hierarchical clustering of samples** (Supp. Fig. 8).
- **Pathway enrichment** (Fig. 6e–f): KOBAS-style hypergeometric test on official KEGG, PANTHER and Reactome mouse gene sets, plus curated BBB panels. ★ marks the pathways in the paper's figure. Click a pathway to list its genes.

### Machine learning
- **Model performance:** results for classifiers trained to tell disease samples from controls.
  - **Leave-one-disease-out** validation trains on two diseases and tests on the third.
  - **Leave-one-out (LOOCV)** is used for the smaller analyses.
- **Top genes:** the genes with the largest logistic-regression weights.

**Please read these results with care.**
- The leave-one-disease-out models **ranked** held-out samples perfectly (ROC-AUC = 1.0). However, their decision threshold did not transfer between studies, so accuracy on held-out stroke and glioblastoma was only 0.5.
- The perfect scores in the glioblastoma-only and multiclass analyses mostly show that these groups are very easy to separate. They are not evidence of a validated biomarker.

---

## 5. Limitations

- **Disease is confounded with study.** Each disease comes from a different lab, protocol and sequencing run. Within-study comparisons (case vs control) are reliable. Differences *between* diseases may partly reflect batch effects.
- **Small sample size.** There are 3–4 replicates per group.
- **Only significant genes are included** in this version, so expression of individual samples, PCA and pathway enrichment are not shown yet (see §7).
- **Microvessel preparations** contain some non-endothelial cells (pericytes, adherent leukocytes), so not every change is endothelial-specific.

---

## 6. Running it yourself

The site is plain HTML, CSS, JavaScript and CSV, so there is nothing to install or build.

### View it locally
Browsers block a page from reading CSV files when it is opened by double-clicking. Start a small local server in this folder instead:

```bash
python3 -m http.server 8000      # Windows: python -m http.server 8000
```

Then open **http://localhost:8000**. Press `Ctrl + C` in the terminal to stop the server.

### Publish it on GitHub Pages
1. Push this folder to the root of a public GitHub repository. Include the hidden `.nojekyll` file.
2. Go to **Settings → Pages**: set Source to *Deploy from a branch*, Branch to `main`, and the folder to `/ (root)`.
3. The site appears at `https://<username>.github.io/<repo-name>/` within a couple of minutes.

### Custom domain (neurovip.bbbos.app)
The repository contains a `CNAME` file with `neurovip.bbbos.app`, which tells GitHub Pages to serve the site at that address. Two more steps are needed:
1. **DNS:** whoever manages the `bbbos.app` domain adds a **CNAME record** with name `neurovip` pointing to `<username>.github.io`.
2. **GitHub:** under **Settings → Pages → Custom domain**, enter `neurovip.bbbos.app`, save, and tick **Enforce HTTPS** once the certificate is issued (this can take up to an hour). `.app` domains only work over HTTPS.

---

## 7. Updating or extending the data

All numbers on the site are read from the CSV files in `data/`. To update the site, replace a file with a new export that has the **same column names**, then commit.

| File | Contents | Columns |
|---|---|---|
| `deg_significant.csv` | Significant genes from this pipeline | `disease, gene_id, symbol, log2fc, padj` |
| `ai_vs_classic.csv` | This pipeline vs classic analysis | `disease, gene_id, symbol, ai_log2fc, ai_padj, mk_log2fc, mk_padj, found_by` |
| `samples_qc.csv` | Per-sample read statistics | `run, disease, condition, group, tissue, raw_reads, adapter_pct, short_removed_pct, unique_aligned, unique_pct, assigned_pct` |
| `studies.csv` | Dataset descriptions | `disease, disease_name, model, tissue, study, geo, n_case, n_control, genes_tested` |
| `ml_metrics.csv` | Model performance | `model, evaluation, held_out, n_test, accuracy, balanced_accuracy, roc_auc, id` |
| `ml_features.csv` | Top model genes | `analysis, symbol, gene_id, coefficient` |

`disease` is always one of `stroke`, `meningitis` or `gb`.

### Classical analysis files

`data/classical_*.csv` are generated, not edited by hand. To regenerate them:

```bash
python3 -m venv .venv                             # once; .venv/ is git-ignored
.venv/bin/pip install -r analysis/requirements.txt   # numpy, pandas, scipy, pydeseq2, pytest
.venv/bin/python analysis/classical_meningitis.py    # about 30 s after the first download
.venv/bin/python -m pytest analysis/tests            # synthetic-data tests, no network
```

A virtual environment is needed on systems that block system-wide `pip` (PEP 668, e.g. Arch, Debian 12+, Homebrew Python).

On its first run the script downloads, and caches in `data/raw/` (git-ignored, about 40 MB):

| Source | Used for |
|---|---|
| GEO `GSE122952_counts.matrix.norm_anno.txt.gz` | the study's normalised counts for DESeq2, PCA, clustering, heatmaps and read counts |
| KEGG REST (`list/pathway/mmu`, `link/pathway/mmu`, `list/mmu`) | KEGG pathways |
| PANTHER `SequenceAssociationPathway3.6.8.txt` (mouse rows) | PANTHER pathways |
| Reactome `Ensembl2Reactome.txt` (mouse rows) | Reactome pathways |
| MGI `MRK_ENSEMBL.rpt` | mapping KEGG symbols and PANTHER MGI IDs to Ensembl |

Use `--offline` to forbid downloads. Delete `data/raw/genesets.tsv` to rebuild the gene sets from the cached sources, or delete `data/raw/` to download fresh releases. The NeuroVIP and classic runs are read from `ai_vs_classic.csv`.

**Planned additions:** pathway enrichment, PCA and sample-level heatmaps, plus volcano plots that include all genes, for all three diseases (meningitis already has them on the Classical analysis tab). These need the following files from the pipeline's `final_output/` folder:

- `*_all_genes.csv`
- `ml/matrices/vst_counts_gene_x_sample.csv`
- `ml/matrices/combined_metadata.csv`
- `extended_analysis/pathway_analysis/`

### Code layout

```
index.html        page layout and all text
CNAME             custom domain for GitHub Pages (neurovip.bbbos.app)
style.css         layout and component styles
theme.css         colours, fonts, spacing (light and dark mode)
js/d3.min.js      D3 v7 charting library
js/runtime.js     loads the CSV files and redraws charts when a filter changes
js/core.js        shared helpers, calculations (DEG counts, overlaps, shared genes, pipeline agreement) and tab controls
js/draw-a.js      Home and Differential expression
js/draw-b.js      Gene search, Shared genes, Pipeline validation, Machine learning and the gene autocomplete list
js/draw-c.js      Classical analysis (meningitis)
analysis/classical_meningitis.py   computes data/classical_*.csv (downloads its inputs to data/raw/)
analysis/tests/   tests for the script on a small synthetic dataset
data/             the CSV files above
```

---

## 8. References

- Kestner R-I, et al. Gene expression dynamics at the neurovascular unit during early regeneration after cerebral ischemia/reperfusion injury in mice. *Front. Neurosci.* 2020;14:280.
- Devraj G, et al. HIF-1α is involved in blood–brain barrier dysfunction and paracellular migration of bacteria in pneumococcal meningitis. *Acta Neuropathol.* 2020;140:183–208.
- Munji RN, et al. Profiling the mouse brain endothelial transcriptome in health and disease models reveals a core blood–brain barrier dysfunction module. *Nat. Neurosci.* 2019;22:1892–1902.
- Dobin A, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 2013;29:15–21.
- Liao Y, Smyth GK, Shi W. featureCounts. *Bioinformatics* 2014;30:923–930.
- Love MI, Huber W, Anders S. DESeq2. *Genome Biol.* 2014;15:550.
