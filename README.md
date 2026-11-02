# BBB RNA-seq Explorer

Interactive dashboard for the BIO F376 Design Project **"Integrated RNA-seq and Machine Learning Pipeline for Profiling Blood–Brain Barrier-Associated Transcriptomic Changes in Ischemic Stroke, Bacterial Meningitis and Glioblastoma"**.

J Lokesh Rao (2023B1A10829H) · Supervisor: Dr. Kavi Devraj, Department of Biological Sciences, BITS Pilani Hyderabad Campus.

## What it shows

| Tab | Contents |
|---|---|
| Overview | DEG counts per disease, source studies, per-sample read QC |
| Differential expression | Volcano plot and searchable gene table per disease, with a fold-change filter |
| Gene lookup | Any gene's log2FC across the three diseases, from this pipeline and the classic analysis |
| Cross-disease | Overlap (UpSet) of significant genes, direction of shared genes, heatmap and table of the 1,647 shared genes |
| Pipeline check | Agreement between this pipeline and the lab's classic workflow |
| Machine learning | Cross-validated model scores and top genes |

It is a plain static site (HTML, CSS, JavaScript and CSV), so it needs no server or build step.

## Host it on GitHub Pages

1. Create a new **public** repository on GitHub (for example `bbb-rnaseq-explorer`).
2. Upload **everything in this folder** to the root of the repository, keeping the folder structure (`index.html`, `style.css`, `theme.css`, `.nojekyll`, `js/`, `data/`). On the GitHub website: **Add file → Upload files**, drag the contents in, then **Commit changes**.
3. Go to **Settings → Pages**. Under **Build and deployment**, set **Source** to *Deploy from a branch*, **Branch** to `main` and the folder to `/ (root)`, then **Save**.
4. After a minute or two the site is live at `https://<your-username>.github.io/<repository-name>/`.

Or from the command line:

```bash
git init
git add .
git commit -m "BBB RNA-seq Explorer"
git branch -M main
git remote add origin https://github.com/<your-username>/<repository-name>.git
git push -u origin main
```

`.nojekyll` is a hidden file; make sure it is uploaded (it tells GitHub Pages to serve the files as they are).

## Run it locally

Browsers block reading CSV files from `file://`, so start a small local server in this folder:

```bash
python3 -m http.server 8000
```

then open <http://localhost:8000>.

## Files

```
index.html          page layout and text
style.css           page styles
theme.css           colours, fonts and spacing (light and dark mode)
js/d3.min.js        D3 v7 (charts)
js/runtime.js       loads the CSV files and redraws the page when filters change
js/core.js          shared helpers, calculations (DEG summary, overlap, shared genes, pipeline agreement) and tab controls
js/draw-a.js        Overview and Differential expression charts/tables
js/draw-b.js        Gene lookup, Cross-disease, Pipeline check and ML charts/tables
data/deg_significant.csv   significant genes (padj < 0.05) per disease from this pipeline
data/ai_vs_classic.csv     this pipeline vs the classic analysis, per disease
data/samples_qc.csv        per-sample read statistics
data/studies.csv           datasets and source studies
data/ml_metrics.csv        model performance
data/ml_features.csv       top model genes
```

## Updating the data

Replace a CSV in `data/` with a new export that keeps the same column names, commit, and the site updates. To add pathway enrichment or PCA tabs, the pipeline's `final_output/*/*_all_genes.csv`, `ml/matrices/vst_counts_gene_x_sample.csv`, `ml/matrices/combined_metadata.csv` and `extended_analysis/pathway_analysis/` files are needed.

## Data sources

- Ischemic stroke: brain microvessels, tMCAO 24 h vs sham — Kestner et al., *Front. Neurosci.* 2020 (GEO GSE131193)
- Bacterial meningitis: brain microvessels, *S. pneumoniae* vs healthy — Devraj G. et al., *Acta Neuropathol.* 2020 (GEO GSE122952)
- Glioblastoma: brain microvessels, tumour-bearing region vs sham
- Processing: FastQC → Trim Galore/Cutadapt → STAR (GRCm39, Ensembl 115) → featureCounts → DESeq2
