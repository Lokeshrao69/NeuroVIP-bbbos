BBB RNA-seq Explorer

An interactive dashboard of blood–brain barrier (BBB) microvessel transcriptomes in mouse models of ischemic stroke, bacterial meningitis and glioblastoma.

🔗 Live site: https://bbb-rnaseq.vercel.app/

BIO F376 Design Project, BITS Pilani Hyderabad Campus Student: J Lokesh Rao (2023B1A10829H) Supervisor: Dr. Kavi Devraj, Department of Biological Sciences

1. What this is

This dashboard presents the results of an RNA-seq analysis pipeline built for this project. It lets you browse them in a web browser without running any code. It is designed around one question:

Do BBB-disrupting diseases of very different origin (ischemic, infectious and neoplastic) share a common transcriptional response in the brain microvasculature?

All three datasets are RNA-seq of brain microvessels from mouse disease models, each compared with its own controls. The raw reads were processed with one uniform pipeline, so the results are directly comparable across diseases.

Headline findings
	Stroke	Meningitis	Glioblastoma
Significant genes (padj < 0.05)	5,227	9,291	8,687
Up-regulated	2,804	4,725	4,468
Down-regulated	2,423	4,566	4,219
1,647 genes are significantly changed in all three diseases.
Of these, 1,322 (80 %) change in the same direction everywhere: 727 up and 595 down.
The strongest shared genes form an endothelial activation and inflammation programme:
adhesion molecules: Selp, Sele, Icam1
Ackr1
acute-phase and antimicrobial genes: Lcn2, Saa3
inflammatory signalling: Il6, Il1b
matrix remodelling: Mmp3
This pipeline agrees almost perfectly with the lab's independent classic analysis. The log2 fold-change correlation is r ≈ 0.9999 for genes found by both.
2. Data
Disease	Model	Tissue	Cases / controls	Source
Ischemic stroke	tMCAO, 24 h, ipsilateral hemisphere, vs sham	Brain microvessels	3 / 3	Kestner et al., Front. Neurosci. 2020 — GEO GSE131193
Bacterial meningitis	Haematogenous S. pneumoniae infection vs healthy	Brain microvessels	4 / 4	Devraj G. et al., Acta Neuropathol. 2020 — GEO GSE122952
Glioblastoma	Tumour-bearing region vs sham	Brain microvessels	3 / 3	(source study to be added)

20 single-end Illumina libraries in total (NCBI SRA run accessions are listed on the Overview tab).

3. How the results were produced
Raw FASTQ (SRA)
  → FastQC 0.11.9                    read quality control
  → Trim Galore 0.6.7 / Cutadapt 5.2 adapter + quality trimming (Q20, min length 30 bp)
  → STAR 2.7.10a                     splice-aware alignment to mouse GRCm39 (Ensembl 115), unique reads only
  → featureCounts (Subread 2.0.3)    reads per gene (exons, MAPQ ≥ 60)
  → DESeq2 1.50.2                    case vs control within each study; Benjamini–Hochberg padj < 0.05
  → cross-disease comparison, machine learning (scikit-learn)

The pipeline ran in Google Colab with data on Google Drive. The same datasets were also analysed with the lab's existing ("classic") workflow, and the dashboard compares the two.

4. How to use the dashboard

The dashboard has six tabs along the top. Charts show details when you hover over them, and every table can be searched.

Overview
Tiles: the number of significant genes (DEGs) per disease, split into up and down, plus the number shared by all three.
Datasets: the model, tissue, source paper and GEO accession for each disease.
Sample quality: reads, adapter content, alignment rate and the percentage of reads assigned to genes, for every sample.
Values highlighted in amber are low. The glioblastoma libraries assign only 47–67 % of reads to genes, against about 89 % for the others, which suggests a different library preparation.
Differential expression
Choose a disease: Stroke / Meningitis / Glioblastoma.
Choose a fold-change threshold: |log2FC| ≥ 0, 1 or 2.
Read the volcano plot:
x-axis: log2 fold change (right = higher in disease)
y-axis: −log10 adjusted p-value (higher = more significant)
colours: red = up, blue = down; grey points are significant but below the fold-change threshold
the 12 strongest genes are labelled
The gene table below lists genes by significance. Type in the search box to find a gene, and click a row to open it in Gene lookup.

The volcano plots show only significant genes, because the exported tables contain only those. The usual grey cloud of unchanged genes near zero is therefore absent.

Gene lookup

Type a gene symbol (e.g. Selp, Icam1, Cldn5, Spp1) and press Enter, or click one of the example buttons.

Bars: the log2 fold change in each disease, from this pipeline (blue) and the classic analysis (purple).
"n.s.": not significant in that disease.
Table: the exact values and adjusted p-values.
Cross-disease
Overlap chart (UpSet plot): each bar is one exact combination of diseases, and the dots underneath show which ones. For example, the highlighted bar is genes significant in all three (1,647).
Direction of the shared genes: how many go up everywhere, down everywhere, or in mixed directions.
Heatmap: the 40 shared genes that are strongly changed in every disease (ranked by their smallest |log2FC| across the three). Red = up, blue = down.
All shared genes: a searchable table of all 1,647. Click a row to open it in Gene lookup.
Pipeline check

Compares this pipeline with the lab's classic analysis of the same data.

Stacked bars: genes found significant by both, by this pipeline only, or by the classic analysis only.
Tiles and scatter plot (pick a disease): points lying on the diagonal mean both pipelines give the same fold change.
	Both	This pipeline only	Classic only
Stroke	4,383	844	14
Meningitis	8,553	738	4
Glioblastoma	7,637	1,050	0

Fold changes are essentially identical. This pipeline calls somewhat more genes significant, probably because of differences in filtering or counting settings.

Machine learning
Model performance: results for classifiers trained to tell disease samples from controls.
Leave-one-disease-out validation trains on two diseases and tests on the third.
Leave-one-out (LOOCV) is used for the smaller analyses.
Top genes: the genes with the largest logistic-regression weights.

Please read these results with care.

The leave-one-disease-out models ranked held-out samples perfectly (ROC-AUC = 1.0). However, their decision threshold did not transfer between studies, so accuracy on held-out stroke and glioblastoma was only 0.5.
The perfect scores in the glioblastoma-only and multiclass analyses mostly show that these groups are very easy to separate. They are not evidence of a validated biomarker.
5. Limitations
Disease is confounded with study. Each disease comes from a different lab, protocol and sequencing run. Within-study comparisons (case vs control) are reliable. Differences between diseases may partly reflect batch effects.
Small sample size. There are 3–4 replicates per group.
Only significant genes are included in this version, so expression of individual samples, PCA and pathway enrichment are not shown yet (see §7).
Microvessel preparations contain some non-endothelial cells (pericytes, adherent leukocytes), so not every change is endothelial-specific.
6. Running it yourself

The site is plain HTML, CSS, JavaScript and CSV, so there is nothing to install or build.

View it locally

Browsers block a page from reading CSV files when it is opened by double-clicking. Start a small local server in this folder instead:

bash
python3 -m http.server 8000      # Windows: python -m http.server 8000

Then open http://localhost:8000. Press Ctrl + C in the terminal to stop the server.

Publish it on GitHub Pages
Push this folder to the root of a public GitHub repository. Include the hidden .nojekyll file.
Go to Settings → Pages: set Source to Deploy from a branch, Branch to main, and the folder to / (root).
The site appears at https://<username>.github.io/<repo-name>/ within a couple of minutes.
7. Updating or extending the data

All numbers on the site are read from the CSV files in data/. To update the site, replace a file with a new export that has the same column names, then commit.

File	Contents	Columns
deg_significant.csv	Significant genes from this pipeline	disease, gene_id, symbol, log2fc, padj
ai_vs_classic.csv	This pipeline vs classic analysis	disease, gene_id, symbol, ai_log2fc, ai_padj, mk_log2fc, mk_padj, found_by
samples_qc.csv	Per-sample read statistics	run, disease, condition, group, tissue, raw_reads, adapter_pct, short_removed_pct, unique_aligned, unique_pct, assigned_pct
studies.csv	Dataset descriptions	disease, disease_name, model, tissue, study, geo, n_case, n_control, genes_tested
ml_metrics.csv	Model performance	model, evaluation, held_out, n_test, accuracy, balanced_accuracy, roc_auc, id
ml_features.csv	Top model genes	analysis, symbol, gene_id, coefficient

disease is always one of stroke, meningitis or gb.

Planned additions: pathway enrichment, PCA and sample-level heatmaps, plus volcano plots that include all genes. These need the following files from the pipeline's final_output/ folder:

*_all_genes.csv
ml/matrices/vst_counts_gene_x_sample.csv
ml/matrices/combined_metadata.csv
extended_analysis/pathway_analysis/
Code layout
index.html        page layout and all text
style.css         layout and component styles
theme.css         colours, fonts, spacing (light and dark mode)
js/d3.min.js      D3 v7 charting library
js/runtime.js     loads the CSV files and redraws charts when a filter changes
js/core.js        shared helpers, calculations (DEG counts, overlaps, shared genes, pipeline agreement) and tab controls
js/draw-a.js      Overview and Differential expression
js/draw-b.js      Gene lookup, Cross-disease, Pipeline check and Machine learning
data/             the CSV files above