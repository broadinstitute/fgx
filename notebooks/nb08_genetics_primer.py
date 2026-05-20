# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo<0.23.4",
#     "jedi<0.20.0",
#     "polars",
#     "httpx",
#     "altair",
#     "python-dotenv",
# ]
# ///

import marimo

__generated_with = "0.23.3"
app = marimo.App(width="medium")

with app.setup:
    import os
    import sys
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import polars as pl
    from dotenv import load_dotenv

    load_dotenv(Path(__file__).resolve().parent.parent / ".env")
    FINNGENIE_TOKEN = os.environ.get("FINNGENIE_TOKEN")
    BASE = "https://finngenie.fi/api/v1"

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    from nb01_pcsk9_walkthrough import client, fetch_json, fetch_tsv


@app.cell
def _():
    mo.md(r"""
    # nb08: How human genetics works

    This notebook walks through the core ideas behind modern human genetics using one gene
    as a running example. Every concept is grounded in real data from FinnGen and collaborating
    biobanks, queried live through the FinnGenie API.

    The approach is bottom-up: start with a gene, ask what the genome tells us about it, and
    build up the layers of evidence that connect a stretch of DNA to a disease and a drug.
    By the end, you'll have seen:

    - **GWAS** - scanning genomes for statistical associations between common variants and traits
    - **Fine-mapping** - resolving which variant at a locus is likely causal
    - **Credible sets** - the probabilistic framework behind fine-mapping
    - **Colocalization** - testing whether two traits share the same causal variant
    - **Rare-variant evidence** - large-effect coding mutations that complement GWAS
    - **QTL layers** - connecting variants to gene expression and protein levels
    - **Gene-disease curation** - expert clinical classifications

    Each section is one or two API calls plus a chart. The concepts generalize to any gene;
    LDLR is chosen because its biology is clear, its data is rich, and its story ends with
    one of the most prescribed drug classes in medicine.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ## 1. The gene: LDLR

    LDLR encodes the low-density lipoprotein receptor, a cell-surface protein that pulls
    LDL cholesterol out of the bloodstream and into cells for degradation. When LDLR works
    well, circulating LDL stays low. When it doesn't - through inherited mutations or reduced
    expression - LDL accumulates in the blood, deposits in artery walls, and drives
    atherosclerosis.

    This is not a hypothesis. Familial hypercholesterolemia (FH), caused by loss-of-function
    mutations in LDLR, was one of the first genetic diseases understood at the molecular level
    (Goldstein and Brown, Nobel Prize 1985). Statins - the most prescribed drug class for
    cardiovascular prevention - work by upregulating LDLR expression. The gene is a textbook
    case where genetics led directly to therapy.

    Let's see what population-scale genomic data adds to that picture.
    """)
    return


@app.cell
def _():
    GENE = "LDLR"
    mo.md(f"**Target gene:** `{GENE}` (chromosome 19, ~44 kb)")
    return (GENE,)


@app.cell
def _(GENE):
    gd = pl.DataFrame(fetch_json(f"/gene_disease/{GENE}"))
    mo.vstack(
        [
            mo.md(r"""
    ## 2. What do experts already know?

    Before looking at statistical associations, check what clinical expert panels have
    curated about this gene. The `gene_disease` endpoint pulls classifications from
    GenCC (Gene Curation Coalition) and Monarch - organizations that formally evaluate
    whether a gene causes a disease based on published evidence, family studies, and
    functional data. Classifications range from "Definitive" (overwhelming evidence)
    to "Limited" (suggestive but incomplete).
    """),
            gd.select(
                "resource",
                "disease_title",
                "classification",
                "mode_of_inheritance",
                "submitter",
            ),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    ClinGen rates LDLR as **"Definitive"** for familial hypercholesterolemia with
    autosomal dominant inheritance. Multiple independent panels agree. The genetics
    community considers this settled.

    The rest of this notebook shows what population-scale data looks like when the
    answer is already known. That's the point: by seeing what a clean signal looks like
    for a well-understood gene, you'll know what to look for in genes where the answer
    isn't settled yet.
    """)
    return


@app.cell
def _(GENE):
    cs = fetch_tsv(f"/credible_sets_by_gene/{GENE}")
    by_type = cs.group_by("data_type").agg(pl.len().alias("n_rows")).sort("n_rows", descending=True)
    mo.vstack(
        [
            mo.md(r"""
    ## 3. GWAS: scanning for common-variant associations

    A genome-wide association study (GWAS) tests hundreds of thousands of genetic variants
    across hundreds of thousands of people, asking: does carrying a particular allele
    correlate with a measurable trait? Each test produces a p-value. With ~10 million
    variants tested, only signals surviving a stringent threshold (p < 5e-8, or
    -log10(p) > 7.3) count as genome-wide significant.

    FinnGenie stores the results of fine-mapped GWAS as **credible sets** - groups of
    variants that, together, have high probability of containing the true causal variant.
    """),
            mo.md(f"### {len(cs):,} credible-set rows for **{GENE}** across all data types"),
            by_type,
        ]
    )
    return (cs,)


@app.cell
def _(cs):
    finngen_gwas = (
        cs.filter(pl.col("data_type") == "GWAS")
        .filter(pl.col("resource").str.contains("finngen"))
        .filter(pl.col("mlog10p").is_not_null())
    )
    trait_summary = (
        finngen_gwas.group_by("trait")
        .agg(
            pl.col("mlog10p").max().alias("max_mlog10p"),
            pl.col("pip").max().alias("max_pip"),
            pl.len().alias("n_variants"),
        )
        .sort("max_mlog10p", descending=True)
        .head(15)
    )
    trait_chart = (
        alt.Chart(trait_summary.select("trait", "max_mlog10p"))
        .mark_bar()
        .encode(
            x=alt.X("max_mlog10p:Q", title="-log10(p)"),
            y=alt.Y("trait:N", sort="-x", title=None),
        )
        .properties(height=380, title="Top FinnGen traits associated near LDLR")
    )
    mo.vstack(
        [
            mo.md(r"""
    Most rows are GWAS - common-variant associations with disease and lab-measured traits.
    But notice the other data types: **eQTL** (variants affecting gene expression),
    **pQTL** (variants affecting protein levels), **caQTL** (chromatin accessibility),
    **sQTL** (splicing). These are the molecular layers that explain *how* a variant
    affects biology, not just *whether* it correlates with a phenotype. We will come back
    to them.

    ### Which traits light up?

    Filter to GWAS rows from FinnGen resources and rank traits by their strongest signal.
    """),
            mo.ui.altair_chart(trait_chart),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    Read that chart as a causal chain made visible by genetics:

    1. **Lipid lab values** (3019900, 3001308, nonHDL) - the direct biochemical consequence
       of LDLR function. These are the strongest signals because they're closest to the gene's
       molecular role.
    2. **Diagnostic codes** (Disorders of lipoprotein metabolism, Pure hypercholesterolaemia,
       E4_LIPOPROT) - clinical recognition of the same biology, one step removed.
    3. **Drug prescriptions** (RX_STATIN, ATC_C10AA) - statins prescribed *because* of high
       cholesterol. The drug signal follows the diagnosis signal.
    4. **Disease endpoints** (Ischaemic heart disease, Coronary atherosclerosis, Myocardial
       infarction, Angina) - the downstream consequences of decades of elevated LDL.

    This is the power of biobank-scale genetics: the entire pathway from molecule to disease
    to treatment is encoded in the association statistics of one gene.
    """)
    return


@app.cell
def _(cs):
    top_pip = (
        cs.filter(pl.col("data_type") == "GWAS")
        .filter(pl.col("pip").is_not_null())
        .sort("pip", descending=True)
        .head(15)
        .select(
            "trait",
            "resource",
            "chr",
            "pos",
            "ref",
            "alt",
            "pip",
            "cs_size",
            "mlog10p",
            "most_severe",
        )
    )
    mo.vstack(
        [
            mo.md(r"""
    ## 4. Fine-mapping: which variant is causal?

    GWAS identifies a *region* where a signal lives, but neighboring variants are correlated
    (in linkage disequilibrium, or LD) and all show association. Fine-mapping uses statistical
    models to assign each variant a **posterior inclusion probability (PIP)** - the probability
    that *this* variant, not its neighbors, is the true causal one.

    A **credible set** is the smallest group of variants whose PIPs sum to a target (typically
    95%). A credible set of size 1 with PIP near 1.0 means fine-mapping has pinpointed the
    causal variant with high confidence. A credible set of size 50 means the signal is real
    but the causal variant is ambiguous.

    ### Highest-PIP GWAS variants near LDLR
    """),
            top_pip,
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    Variants with PIP close to 1.0 in small credible sets are the ones fine-mapping is most
    confident about. Notice the `most_severe` column: it tells you the predicted consequence
    of each variant on the protein (missense, intron_variant, etc.). An intronic variant with
    high PIP might affect gene regulation rather than protein structure - that's where the
    QTL layers (eQTL, caQTL) become essential for interpretation.
    """)
    return


@app.cell
def _(cs):
    lead = (
        cs.filter(pl.col("data_type") == "GWAS")
        .filter(pl.col("mlog10p").is_not_null())
        .sort("mlog10p", descending=True)
        .head(1)
    )
    lead_variant = f"{lead['chr'][0]}:{lead['pos'][0]}:{lead['ref'][0]}:{lead['alt'][0]}"
    near = pl.DataFrame(fetch_json(f"/nearest_genes/{lead_variant}", n=3))
    mo.vstack(
        [
            mo.md(r"""
    ## 5. From variant to gene: closing the loop

    A GWAS hit is a genomic coordinate, not a gene name. The variant with the strongest
    signal might sit between two genes, inside an intron, or in a regulatory element far
    from its target. Assigning the right gene is one of the hardest problems in genetics.

    The simplest heuristic: which protein-coding gene is nearest?
    """),
            mo.md(f"### Nearest genes to lead variant `{lead_variant}`"),
            near.select("gene_name", "gene_type", "distance", "hgnc_name"),
        ]
    )
    return (lead_variant,)


@app.cell
def _():
    mo.md(r"""
    Distance zero means the variant falls *inside* the gene body. For LDLR, the answer is
    straightforward: the top variant is in the gene itself. But for many GWAS loci, the
    nearest gene isn't the causal one - the signal may act through a distal enhancer. That's
    why nearest-gene assignment is a starting point, not a conclusion, and why the QTL and
    colocalization layers below add real information.
    """)
    return


@app.cell
def _(lead_variant):
    try:
        coloc = pl.DataFrame(fetch_json(f"/colocalization_by_variant/{lead_variant}"))
        coloc_header = mo.md(f"### {len(coloc):,} colocalization pairs at `{lead_variant}`")
    except Exception as e:
        coloc = pl.DataFrame()
        coloc_header = mo.md(f"### Colocalization lookup failed for `{lead_variant}` ({e})")

    mo.vstack(
        [
            mo.md(r"""
    ## 6. Colocalization: do two traits share a causal variant?

    When two traits both show a GWAS signal at the same locus, are they driven by the
    *same* causal variant, or by two different variants that happen to be nearby?
    Colocalization tests answer this by comparing the fine-mapped credible sets.

    High colocalization (PP.H4 > 0.8, or high CLPP) means the two traits likely share a
    mechanism at this locus. This is how genetics connects risk factors to disease outcomes:
    if LDL cholesterol and coronary heart disease colocalize at LDLR, the shared variant
    provides Mendelian randomization-style evidence that LDL *causes* heart disease, not
    just that the two are correlated.
    """),
            coloc_header,
        ]
    )
    return (coloc,)


@app.cell
def _(coloc, lead_variant):
    mo.stop(
        coloc.is_empty(),
        mo.md("_No colocalization data available for this variant._"),
    )
    top_coloc = (
        coloc.filter(pl.col("hit1_mlog10p").is_not_null())
        .sort("hit1_mlog10p", descending=True)
        .head(20)
        .with_columns(pl.concat_str([pl.col("trait1"), pl.lit(" / "), pl.col("trait2")]).alias("pair"))
    )
    coloc_chart = (
        alt.Chart(top_coloc.select("pair", "hit1_mlog10p", "data_type1"))
        .mark_bar()
        .encode(
            x=alt.X("hit1_mlog10p:Q", title="-log10(p) for trait 1"),
            y=alt.Y("pair:N", sort="-x", title=None),
            color=alt.Color("data_type1:N", title="Data type"),
        )
        .properties(height=450, title=f"Top colocalizing trait pairs at {lead_variant}")
    )
    mo.vstack(
        [
            mo.ui.altair_chart(coloc_chart),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    The colocalization results read like the causal chain from section 3, but now with
    statistical backing: lipid lab values colocalize with statin prescriptions, which
    colocalize with cardiovascular endpoints. The shared causal variant at LDLR ties
    them together mechanistically, not just statistically.

    When you see eQTL traits in the coloc list, that's especially informative: it means
    the same variant that drives disease risk *also* changes LDLR gene expression, pointing
    to a regulatory mechanism rather than a protein-coding change.
    """)
    return


@app.cell
def _(GENE):
    exome = fetch_tsv(f"/exome_results_by_gene/{GENE}")
    mo.vstack(
        [
            mo.md(r"""
    ## 7. Rare variants: the other end of the frequency spectrum

    GWAS finds common variants (carried by >1% of people) with small effects. But the same
    gene can also harbor rare variants (carried by <0.1%) with much larger effects.

    Exome sequencing catalogs coding variants - the ones that change the protein sequence.
    These come in flavors:
    - **pLoF** (predicted loss of function) - stop-gains, frameshifts, splice donors -
      variants expected to break the protein entirely
    - **missense** - single amino acid changes that may or may not impair function
    - **synonymous** - no protein change; useful as negative controls

    For LDLR, rare pLoF variants cause familial hypercholesterolemia with large effect sizes.
    """),
            mo.md(f"### {len(exome):,} exome variant-trait rows for **{GENE}**"),
        ]
    )
    return (exome,)


@app.cell
def _(exome):
    from nb04_gene_exome_burden import prepare_deleterious

    deleterious = prepare_deleterious(exome)
    top_exome = deleterious.sort("mlog10p", descending=True).head(20)

    points = (
        alt.Chart(top_exome.select("trait_variant", "beta", "ci_lo", "ci_hi", "annotation"))
        .mark_circle(size=80)
        .encode(
            x=alt.X("beta:Q", title="Per-allele effect (beta)"),
            y=alt.Y("trait_variant:N", sort="-x", title=None),
            color=alt.Color(
                "annotation:N",
                scale=alt.Scale(domain=["pLoF", "missense"], range=["#d62728", "#1f77b4"]),
                title="VEP class",
            ),
        )
    )
    bars = (
        alt.Chart(top_exome.select("trait_variant", "ci_lo", "ci_hi"))
        .mark_rule()
        .encode(x="ci_lo:Q", x2="ci_hi:Q", y=alt.Y("trait_variant:N", sort="-x"))
    )
    zero = alt.Chart(pl.DataFrame({"x": [0.0]})).mark_rule(strokeDash=[4, 3], color="black").encode(x="x:Q")
    forest = (bars + points + zero).properties(
        height=500,
        title="Forest plot: top rare coding variants in LDLR (beta +/- 1.96 SE)",
    )
    mo.vstack(
        [
            mo.ui.altair_chart(forest),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    Compare this forest plot to the GWAS bar chart above. The GWAS common variants have
    tiny effect sizes (betas in the hundredths). The rare coding variants have effects
    10-100x larger - a single pLoF allele can shift cholesterol by a full standard deviation.

    This is the **allelic series**: the same gene, a spectrum of variant frequencies and
    effect sizes, all pointing at the same biology. Common variants nudge LDLR expression;
    rare variants break it. Both end up at cholesterol. This complementarity is why modern
    genetics uses both GWAS and exome analysis together - each fills gaps the other can't.
    """)
    return


@app.cell
def _(cs):
    qtl_types = (
        cs.filter(pl.col("data_type") != "GWAS")
        .group_by("data_type")
        .agg(
            pl.len().alias("n_rows"),
            pl.col("resource").n_unique().alias("n_resources"),
            pl.col("trait").n_unique().alias("n_distinct_traits"),
        )
        .sort("n_rows", descending=True)
    )
    mo.vstack(
        [
            mo.md(r"""
    ## 8. Molecular layers: eQTL, pQTL, and chromatin

    The data types beyond GWAS in the credible-set table are **quantitative trait loci**
    (QTLs) - variants associated not with disease, but with molecular measurements:

    - **eQTL**: gene expression levels (mRNA) in specific tissues or cell types
    - **pQTL**: protein levels in plasma
    - **sQTL**: splicing patterns
    - **caQTL**: chromatin accessibility (whether DNA is open for transcription)

    When a GWAS variant also shows up as an eQTL for the same gene, it suggests a regulatory
    mechanism: the variant changes how much protein the cell makes, which changes the trait.

    ### QTL layers available near LDLR
    """),
            qtl_types,
        ]
    )
    return


@app.cell
def _(cs):
    eqtl_cells = (
        cs.filter(pl.col("data_type") == "eQTL")
        .filter(pl.col("pip").is_not_null())
        .group_by("cell_type")
        .agg(
            pl.col("pip").max().alias("max_pip"),
            pl.len().alias("n_variants"),
        )
        .sort("max_pip", descending=True)
        .head(15)
    )
    eqtl_chart = (
        alt.Chart(eqtl_cells.select("cell_type", "max_pip"))
        .mark_bar()
        .encode(
            x=alt.X("max_pip:Q", title="Max PIP (eQTL)"),
            y=alt.Y("cell_type:N", sort="-x", title=None),
        )
        .properties(height=350, title="Cell types with eQTL signal near LDLR")
    )
    mo.vstack(
        [
            mo.md(
                "### In which cell types does this locus regulate LDLR expression?\n\n"
                "Each bar is the strongest eQTL PIP in that cell type. High PIP in liver "
                "would be the expected biology (hepatocytes are where LDLR clears LDL), but "
                "context-specific eQTL data is patchy - the absence of liver here reflects "
                "data availability, not biology."
            ),
            mo.ui.altair_chart(eqtl_chart),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    ## 9. Synthesis: convergent evidence

    Step back and see what one gene query produced:

    | Layer | What it shows | LDLR result |
    |---|---|---|
    | Gene-disease curation | Expert clinical classification | "Definitive" for familial hypercholesterolemia |
    | GWAS | Common-variant associations | Lipids, statin Rx, heart disease - the full causal chain |
    | Fine-mapping | Which variant is causal | High-PIP variants inside the gene body |
    | Colocalization | Shared mechanism across traits | Lipid labs and cardiovascular endpoints share a signal |
    | Exome | Rare coding variants | pLoF with 10-100x larger effect on cholesterol |
    | eQTL/pQTL | Molecular mechanism | Variants affect gene expression in relevant tissues |

    No single layer is conclusive on its own. GWAS gives association, not causation.
    Fine-mapping gives a variant, not a gene. Exome gives a gene, not a mechanism. eQTL
    gives a mechanism, not a disease. But when all layers converge on the same story -
    LDLR variants reduce receptor function, raise LDL, cause atherosclerosis - the
    combined evidence is far stronger than any single line.

    This convergence is exactly what drug developers look for when evaluating a new target.
    PCSK9 inhibitors (nb01's story) went through the same logic: GWAS found the locus,
    rare-variant studies confirmed the direction of effect, and the drug was designed to
    mimic what loss-of-function carriers get for free.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ## To extend

    - **Swap LDLR for PCSK9** and re-run. PCSK9 loss-of-function *lowers* LDL (opposite
      direction from LDLR loss-of-function) because PCSK9 normally degrades LDLR. Same
      pathway, opposite allelic-series direction - visible in the forest plot.
    - **Try a gene with weaker evidence** - e.g. `SORT1`, `HMGCR`, or `NPC1L1` - and see
      how the layers look when the signal is noisier or the biology less clear.
    - **Start from a phenotype instead of a gene.** Use nb03's approach: pick a FinnGen
      phenocode like `I9_CHD` (coronary heart disease) and pull all credible sets for that
      trait. You'll find LDLR alongside dozens of other loci - that's the genome-wide view
      this notebook deliberately avoided to keep the focus on one gene.
    - **Explore the full FinnGenie catalog.** nb07 lists all 29 datasets the API exposes.
      This notebook touched GWAS, exome, colocalization, gene-disease, eQTL, and
      nearest-genes - there are more (summary stats, peak-to-genes, chromatin peaks).
    """)
    return


if __name__ == "__main__":
    app.run()
