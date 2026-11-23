"""Tests for analysis/classical_meningitis.py on a small synthetic dataset (no network).

Run:  python -m pytest analysis/tests
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import classical_meningitis as cm  # noqa: E402

N_GENES = 400
N_UP = N_DOWN = 40
NAMED = [g[0] for g in cm.PAPER_GENES]


def gene_ids():
    return [f"ENSMUSG{i:011d}" for i in range(N_GENES)]


def symbols():
    s = [f"Gene{i}" for i in range(N_GENES)]
    # named paper genes: Hif1a strongly up (index 0), Cldn5 strongly down (index N_UP)
    s[0], s[N_UP] = "Hif1a", "Cldn5"
    for k, g in enumerate(n for n in NAMED if n not in ("Hif1a", "Cldn5")):
        s[200 + k] = g          # unchanged region
    return s


def make_counts(seed=1):
    rng = np.random.default_rng(seed)
    base = rng.uniform(50, 2000, N_GENES)
    fc = np.ones(N_GENES)
    fc[:N_UP] = 8.0
    fc[N_UP:N_UP + N_DOWN] = 1 / 8.0
    cols = {}
    for name in cm.GEO_TO_RUN:
        mu = base * (fc if name.startswith("Meningitis") else 1.0)
        cols[name] = rng.negative_binomial(20, 20 / (20 + mu))
    return pd.DataFrame(cols, index=gene_ids())


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    d = tmp_path_factory.mktemp("data")
    raw = d / "raw"
    raw.mkdir()
    counts = make_counts()
    ids, syms = gene_ids(), symbols()
    m = counts.copy()
    m.insert(0, "Ensembl gene", syms)
    m.insert(0, "Ensembl gene id", ids)
    with gzip.open(raw / cm.COUNTS_FILE, "wt") as f:
        m.to_csv(f, sep="\t", index=False)

    pd.DataFrame([dict(disease="meningitis", disease_name="Bacterial meningitis", genes_tested=N_GENES)]) \
        .to_csv(d / "studies.csv", index=False)
    # pipeline results: only significant genes are stored, like the real file
    rows = [dict(disease="meningitis", gene_id=ids[i], symbol=syms[i],
                 ai_log2fc=3.0 if i < N_UP else -3.0, ai_padj=1e-10,
                 mk_log2fc=2.9 if i < N_UP else -2.9, mk_padj=1e-9, found_by="both")
            for i in range(N_UP + N_DOWN)]
    rows.append(dict(disease="stroke", gene_id=ids[0], symbol=syms[0], ai_log2fc=9, ai_padj=1e-9,
                     mk_log2fc=9, mk_padj=1e-9, found_by="both"))
    pd.DataFrame(rows).to_csv(d / "ai_vs_classic.csv", index=False)

    # gene sets: one pathway made of the up genes, one of the down genes, one random
    gs = [dict(db="KEGG", pathway_id="mmu04066", pathway="HIF-1 signaling pathway", gene_id=ids[i]) for i in range(20)]
    gs += [dict(db="PANTHER", pathway_id="P00056", pathway="VEGF signaling pathway", gene_id=ids[N_UP + i]) for i in range(20)]
    gs += [dict(db="Reactome", pathway_id="R-MMU-1", pathway="Background", gene_id=ids[300 + i]) for i in range(20)]
    pd.DataFrame(gs).to_csv(raw / cm.GENESETS_FILE, sep="\t", index=False)

    cm.main(["--data", str(d), "--offline"])
    return d


def out(d, name):
    return pd.read_csv(d / f"classical_{name}.csv", keep_default_na=False, na_values=[""])


# ---------------------------------------------------------------- helpers
def test_bh_matches_r():
    # R: p.adjust(c(0.01, 0.04, 0.03, 0.2), "BH") -> 0.04 0.05333 0.05333 0.20
    np.testing.assert_allclose(cm.bh([0.01, 0.04, 0.03, 0.2]), [0.04, 0.16 / 3, 0.16 / 3, 0.2])


@pytest.mark.parametrize("lfc,padj,expected", [
    (1.0, 0.01, "up"), (-1.0, 0.01, "down"), (0.3, 0.01, "below_fc"),
    (2.0, 0.2, ""), (np.nan, 0.01, ""), (1.0, np.nan, ""),
])
def test_call(lfc, padj, expected):
    assert cm.call(lfc, padj) == expected


def test_regulated_in_all_four():
    c = pd.DataFrame({"Meningitis_1": [100, 10, 100], "Meningitis_2": [100, 10, 15],
                      "Meningitis_3": [100, 10, 100], "Meningitis_4": [100, 10, 100],
                      "Healthy_1": [10, 100, 10], "Healthy_2": [10, 100, 10],
                      "Healthy_3": [10, 100, 10], "Healthy_4": [10, 100, 10]}, index=["a", "b", "c"])
    s = pd.DataFrame({"sample": c.columns, "condition": ["case"] * 4 + ["control"] * 4})
    # c: one infected sample is below 1.5-fold (log2(15.5 / 10.5) = 0.56) -> not regulated in all four
    assert list(cm.regulated_in_all_four(c, s)) == [1, -1, 0]


# ---------------------------------------------------------------- end to end
def test_all_outputs_written(dataset):
    for name in ["methods", "summary", "deg", "genes", "counts", "top", "paper_top", "pca", "sample_dist",
                 "dendrogram", "sample_heatmap", "pathways", "pathway_genes"]:
        assert (dataset / f"classical_{name}.csv").exists(), name


def test_summary_has_stable_keys(dataset):
    s = out(dataset, "summary")
    keys = set(zip(s.source, s.tkey))
    for src in cm.SOURCES:
        assert (src, "padj") in keys and (src, "strict") in keys
    assert ("paper", "reported") in keys


def test_deseq2_recovers_planted_genes(dataset):
    deg = out(dataset, "deg").set_index("gene_id")
    ids = gene_ids()
    up = deg.loc[ids[:N_UP], "pa_call"]
    down = deg.loc[ids[N_UP:N_UP + N_DOWN], "pa_call"]
    assert (up == "up").mean() > 0.9
    assert (down == "down").mean() > 0.9
    assert (deg.loc[ids[N_UP + N_DOWN:], "pa_call"].isin(["up", "down"])).mean() < 0.05
    # planted genes are also regulated in all four infected samples
    assert (deg.loc[ids[:N_UP], "pa_all4"] == 1).mean() > 0.9


def test_named_genes_compared_per_run(dataset):
    g = out(dataset, "genes").set_index("symbol")
    assert g.at["Hif1a", "pa_observed"] == "up" and g.at["Hif1a", "pa_agreement"] == "agrees"
    assert g.at["Cldn5", "nv_observed"] == "down" and g.at["Cldn5", "nv_agreement"] == "agrees"
    # Tnf is unchanged in the fake data but the paper says up -> must be reported as differing
    assert g.at["Tnf", "pa_agreement"] == "differs"


def test_read_counts_are_raw_matrix_values(dataset):
    rc = out(dataset, "counts")
    c = make_counts()
    r = rc[(rc.symbol == "Hif1a") & (rc["sample"] == "Meningitis_1")].iloc[0]
    assert r["count"] == c.at[gene_ids()[0], "Meningitis_1"]
    assert r["run"] == "SRR8244647"


def test_top_genes_per_source(dataset):
    t = out(dataset, "top")
    assert set(t.source) == set(cm.SOURCES)
    pa = t[(t.source == "pa") & (t.rank_by == "padj") & (t.direction == "up")]
    assert len(pa) == cm.HEATMAP_TOP and (pa.log2fc > 0).all()


def test_pca_separates_groups(dataset):
    p = out(dataset, "pca")
    for method in ("top500", "scaled"):
        x = p[p.method == method]
        assert x[x.condition == "case"].PC1.min() > x[x.condition == "control"].PC1.max()


def test_clustering_splits_groups(dataset):
    d = out(dataset, "sample_dist")
    order = d.drop_duplicates("sample_a").sort_values("order_a")
    cond = list(order.condition_a)
    # leaves of each condition are contiguous in the dendrogram order
    assert cond == sorted(cond) or cond == sorted(cond, reverse=True)
    assert len(out(dataset, "dendrogram")) == 3 * 7    # 7 merges for 8 samples


def test_enrichment_finds_planted_pathways(dataset):
    p = out(dataset, "pathways")
    pa = p[p.source == "pa"].set_index("pathway_id")
    assert pa.at["mmu04066", "q_enrich"] < 0.01 and pa.at["mmu04066", "direction"] == "mostly up"
    assert pa.at["P00056", "q_enrich"] < 0.01 and pa.at["P00056", "direction"] == "mostly down"
    assert pa.at["mmu04066", "paper"] == "text"
    members = out(dataset, "pathway_genes")
    assert set(members[members.pathway_id == "mmu04066"].gene_id) == set(gene_ids()[:20])


def test_methods_text_generated_from_numbers(dataset):
    m = out(dataset, "methods").set_index("step")
    s = out(dataset, "summary")
    pa = s[(s.source == "pa") & (s.tkey == "padj")].iloc[0]
    assert f"{int(pa.n_up):,} up" in m.at["Differential expression", "detail"]
