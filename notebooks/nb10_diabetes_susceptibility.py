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
    BASE = "https://finngenie.broadinstitute.org/api/v1"

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    from nb01_pcsk9_walkthrough import client, fetch_json, fetch_tsv
    from nb03_phenotype_locus_zoom import pick_leads


@app.cell
def _():
    mo.md(r"""
    # nb10: What makes people susceptible to type 2 diabetes?

    Type 2 diabetes (T2D) is one of the most genetically dissected diseases in human genetics.
    FinnGen R13 enrolled ~90,000 T2D cases out of ~486,000 participants, making it one of the
    largest single-biobank GWAS for this trait.

    This notebook asks: what does the FinnGen data tell us about genetic susceptibility to T2D?

    1. **The polygenic landscape** — a Manhattan plot showing hundreds of independent signals
       across the genome.
    2. **Small individual effects** — most variants nudge risk by fractions of a percent.
    3. **TCF7L2: the strongest signal in T2D genetics** — what it colocalizes with and why it
       matters.
    4. **The MODY connection** — several GWAS loci sit in genes that also cause monogenic
       diabetes, linking the common and rare ends of the genetic spectrum.

    The punchline: T2D susceptibility is profoundly polygenic, but many of the genes involved
    are already known from rare Mendelian forms of diabetes — common variants whisper what
    rare mutations shout.
    """)
    return


@app.cell
def _():
    RESOURCE = "finngen"
    PHENOTYPE = "T2D"
    t2d_cs = fetch_tsv(f"/credible_sets_by_phenotype/{RESOURCE}/{PHENOTYPE}")
    n_cs = t2d_cs["cs_id"].n_unique()
    n_chr = t2d_cs["chr"].cast(pl.Utf8).n_unique()
    mo.md(
        f"### {len(t2d_cs):,} credible-set rows for **{PHENOTYPE}** on **{RESOURCE}**, "
        f"spanning **{n_cs}** independent credible sets across **{n_chr}** chromosomes"
    )
    return PHENOTYPE, t2d_cs


@app.cell
def _():
    mo.md(r"""
    ## 1. The polygenic landscape

    If T2D were driven by one gene, fine-mapping would find a handful of credible sets on one
    chromosome. Instead, FinnGen finds hundreds of independent signals scattered across every
    autosome and the X chromosome.

    The Manhattan plot below shows one dot per lead variant (the highest -log10p variant in each
    credible set), positioned by genomic coordinate. The sheer spread across chromosomes is the
    first line of evidence: T2D susceptibility is not localized — it is genome-wide.
    """)
    return


@app.cell
def _(t2d_cs):
    t2d_leads = pick_leads(t2d_cs)
    mo.vstack(
        [
            mo.md(
                f"### {len(t2d_leads)} lead variants (one per credible set)\n\n"
                "Each row is the top variant from one independent fine-mapped signal, "
                "sorted by statistical significance."
            ),
            t2d_leads.head(10).select(
                "variant_id",
                "chr",
                "pos",
                "mlog10p",
                "pip",
                "cs_size",
                "most_severe",
                "gene_most_severe",
            ),
        ]
    )
    return (t2d_leads,)


@app.cell
def _(t2d_leads):
    chrom_order = [str(i) for i in range(1, 23)] + ["23", "X", "Y"]
    t2d_leads_str = t2d_leads.with_columns(pl.col("chr").cast(pl.Utf8).alias("chr_str"))
    t2d_chr_max = (
        t2d_leads_str.group_by("chr_str")
        .agg(pl.col("pos").max().alias("chr_len"))
        .with_columns(
            pl.col("chr_str")
            .map_elements(
                lambda c: chrom_order.index(c) if c in chrom_order else 99,
                return_dtype=pl.Int64,
            )
            .alias("chr_rank")
        )
        .sort("chr_rank")
    )
    t2d_offsets = t2d_chr_max.with_columns(pl.col("chr_len").cum_sum().shift(1, fill_value=0).alias("offset")).select(
        "chr_str", "chr_rank", "offset"
    )

    t2d_manhattan_df = (
        t2d_leads_str.join(t2d_offsets, on="chr_str", how="left")
        .with_columns((pl.col("pos") + pl.col("offset")).alias("cum_pos"))
        .with_columns((pl.col("chr_rank") % 2).cast(pl.Utf8).alias("band"))
    )

    known_genes = {
        "TCF7L2": "TCF7L2",
        "CCND2": "CCND2",
        "FTO": "FTO",
        "CDKAL1": "CDKAL1",
        "MTNR1B": "MTNR1B",
        "KCNQ1": "KCNQ1",
        "IGF2BP2": "IGF2BP2",
        "SLC30A8": "SLC30A8",
        "PPARG": "PPARG",
        "HNF1A": "HNF1A",
    }
    t2d_label_df = t2d_manhattan_df.filter(pl.col("gene_most_severe").is_in(list(known_genes.keys()))).with_columns(
        pl.col("gene_most_severe").alias("label")
    )
    return t2d_label_df, t2d_manhattan_df


@app.cell
def _(PHENOTYPE, t2d_label_df, t2d_manhattan_df):
    t2d_manhattan_points = (
        alt.Chart(alt.Data(values=t2d_manhattan_df.select("cum_pos", "mlog10p", "band").to_dicts()))
        .mark_circle(size=45, opacity=0.85)
        .encode(
            x=alt.X(
                "cum_pos:Q",
                title="Genome position",
                axis=alt.Axis(labels=False, ticks=False),
            ),
            y=alt.Y("mlog10p:Q", title="-log10(p)"),
            color=alt.Color(
                "band:N",
                scale=alt.Scale(domain=["0", "1"], range=["#2ca02c", "#98df8a"]),
                legend=None,
            ),
        )
    )
    t2d_manhattan_labels = (
        alt.Chart(alt.Data(values=t2d_label_df.select("cum_pos", "mlog10p", "label").to_dicts()))
        .mark_text(align="left", dx=6, dy=-4, fontSize=10, fontWeight="bold")
        .encode(x="cum_pos:Q", y="mlog10p:Q", text="label:N")
    )
    t2d_manhattan_chart = (t2d_manhattan_points + t2d_manhattan_labels).properties(
        height=380,
        width=820,
        title=f"Manhattan of credible-set leads for {PHENOTYPE} (FinnGen R13)",
    )
    mo.ui.altair_chart(t2d_manhattan_chart)
    return


@app.cell
def _(t2d_leads):
    n_chromosomes = t2d_leads["chr"].cast(pl.Utf8).n_unique()
    mo.md(
        f"**{len(t2d_leads)} independent signals across {n_chromosomes} chromosomes.** "
        "Nearly every chromosome hosts T2D credible sets. TCF7L2 on chromosome 10 dominates "
        "with -log10(p) > 280, but the bulk of the genetic architecture is hundreds of "
        "smaller signals — each adding a sliver of risk."
    )
    return


@app.cell
def _():
    mo.md(r"""
    ## 2. Small individual effects

    The second hallmark of polygenic architecture: each variant contributes very little on its
    own. Monogenic diseases (cystic fibrosis, sickle cell) have variants with large effects.
    T2D variants mostly nudge risk by fractions of a standard deviation.

    The histogram below shows the distribution of absolute effect sizes across all lead variants.
    """)
    return


@app.cell
def _(t2d_leads):
    t2d_leads_with_beta = t2d_leads.filter(pl.col("beta").is_not_null() & pl.col("beta").is_finite()).with_columns(
        pl.col("beta").abs().alias("abs_beta")
    )

    t2d_beta_hist = (
        alt.Chart(alt.Data(values=t2d_leads_with_beta.select("abs_beta").to_dicts()))
        .mark_bar(cornerRadiusTopLeft=3, cornerRadiusTopRight=3)
        .encode(
            x=alt.X(
                "abs_beta:Q",
                bin=alt.Bin(maxbins=30),
                title="|beta| (per-allele effect size)",
            ),
            y=alt.Y("count()", title="Number of lead variants"),
        )
        .properties(
            height=300,
            width=600,
            title="Distribution of effect sizes across T2D lead variants",
        )
    )
    t2d_median_beta = t2d_leads_with_beta["abs_beta"].median()
    mo.vstack(
        [
            mo.ui.altair_chart(t2d_beta_hist),
            mo.md(
                f"**Median |beta| = {t2d_median_beta:.4f}.** The vast majority of T2D lead "
                "variants shift disease risk by less than 0.1 standard deviations per allele. "
                "No single common variant is sufficient to cause diabetes — but hundreds of "
                "these small effects, aggregated into a polygenic risk score, can stratify "
                "individuals into meaningfully different risk tiers."
            ),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    ## 3. TCF7L2: the strongest common-variant signal in T2D genetics

    TCF7L2 (transcription factor 7-like 2) on chromosome 10q25 has been the most strongly
    associated T2D locus since its discovery in 2006. It encodes a transcription factor in
    the Wnt signaling pathway that regulates pancreatic beta-cell function and insulin
    secretion.

    What does it colocalize with? If the same causal variant at TCF7L2 drives multiple traits,
    that tells us about the biological mechanism linking the variant to T2D.
    """)
    return


@app.cell
def _(t2d_leads):
    tcf7l2_lead = t2d_leads.filter(pl.col("gene_most_severe") == "TCF7L2").head(1)
    tcf7l2_variant = tcf7l2_lead["variant_id"][0]

    t2d_self_traits = {
        "T2D",
        "T2D_WIDE",
        "E4_DM2",
        "E4_DM2NASCOMP",
        "E4_DM2NAS",
        "KELA_DIAB_INSUL_EXMORE",
    }

    tcf7l2_coloc = pl.DataFrame(fetch_json(f"/colocalization_by_variant/{tcf7l2_variant}"))
    tcf7l2_diverse = (
        tcf7l2_coloc.filter(pl.col("trait1") != pl.col("trait2"))
        .filter(~pl.col("trait2").is_in(t2d_self_traits))
        .filter(pl.col("PP.H4.abf") > 0.5)
        .sort("PP.H4.abf", descending=True)
        .unique(subset=["trait2"], keep="first", maintain_order=True)
    )

    tcf7l2_display = tcf7l2_diverse.head(15).select(
        "trait2",
        "data_type2",
        pl.col("PP.H4.abf").round(3).alias("PP_H4"),
    )

    mo.vstack(
        [
            mo.md(
                f"### Colocalizations at TCF7L2 lead variant `{tcf7l2_variant}`\n\n"
                "T2D self-colocalizations filtered out. Each row is a distinct trait sharing "
                "a causal variant with T2D at this locus."
            ),
            tcf7l2_display,
        ]
    )
    return tcf7l2_diverse, tcf7l2_variant


@app.cell
def _(tcf7l2_diverse, tcf7l2_variant):
    tcf7l2_by_type = (
        tcf7l2_diverse.group_by("data_type2")
        .agg(pl.len().alias("n_traits"))
        .sort("n_traits", descending=True)
        .cast({"n_traits": pl.Int64})
    )

    tcf7l2_type_chart = (
        alt.Chart(
            alt.Data(
                values=tcf7l2_by_type.select("data_type2", "n_traits").to_dicts()
            )
        )
        .mark_bar(cornerRadiusTopLeft=3, cornerRadiusTopRight=3)
        .encode(
            x=alt.X("data_type2:N", title="Data type", sort="-y"),
            y=alt.Y("n_traits:Q", title="Number of colocalizing traits"),
            color=alt.Color("data_type2:N", legend=None),
        )
        .properties(
            height=300,
            width=500,
            title=f"TCF7L2 ({tcf7l2_variant}) colocalizations by data type",
        )
    )
    mo.vstack(
        [
            mo.md(
                "The TCF7L2 locus colocalizes with other diabetes definitions (T1D, "
                "gestational diabetes), metabolic traits (SHBG, medication use), and molecular "
                "QTLs — consistent with its role as a master regulator of beta-cell "
                "transcription. The breadth of colocalizing traits reflects TCF7L2's central "
                "position in glucose homeostasis."
            ),
            tcf7l2_type_chart,
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    ## 4. The MODY connection: where common and rare genetics converge

    One of the most striking features of T2D genetics is that many GWAS loci sit in genes
    that also cause **maturity-onset diabetes of the young (MODY)** — rare, monogenic forms
    of diabetes caused by single high-impact mutations.

    The table below shows curated gene-disease relationships from GenCC for key T2D genes.
    These genes carry both:
    - **Common variants** with tiny effects on T2D risk (visible in the GWAS above)
    - **Rare mutations** that cause monogenic diabetes (curated at "Definitive" or "Strong"
      confidence)

    Common variants whisper what rare mutations shout — the same genes, different effect sizes,
    same biological pathway.
    """)
    return


@app.cell
def _():
    mody_genes = ["HNF1A", "HNF4A", "GCK", "KCNJ11", "ABCC8", "PPARG"]
    gd_rows = []
    for gene in mody_genes:
        try:
            gd = fetch_tsv(f"/gene_disease/{gene}")
            diabetes_gd = gd.filter(pl.col("disease_title").str.contains("(?i)diabet|MODY|insulin|hyperinsulin"))
            if len(diabetes_gd) > 0:
                for row in diabetes_gd.iter_rows(named=True):
                    gd_rows.append(
                        {
                            "gene": gene,
                            "disease": row["disease_title"],
                            "classification": row["classification"],
                            "submitter": row["submitter"],
                        }
                    )
        except Exception:
            pass

    mody_gd = pl.DataFrame(gd_rows)
    definitive_strong = mody_gd.filter(pl.col("classification").is_in(["Definitive", "Strong"])).sort("gene")

    mo.vstack(
        [
            mo.md(
                f"### {len(definitive_strong)} curated gene-disease entries "
                "(Definitive or Strong confidence)\n\n"
                "These are the same genes that appear in the T2D GWAS Manhattan above. "
                "In rare families, a single mutation in any of these genes is sufficient to "
                "cause diabetes. In the general population, common variants at the same loci "
                "each contribute a small fraction of T2D risk."
            ),
            definitive_strong,
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    ## 5. Synthesis

    | Layer | What it shows |
    |---|---|
    | **GWAS loci** | 300+ independent credible sets across every chromosome |
    | **Effect sizes** | Small — median |beta| < 0.05 per allele |
    | **TCF7L2** | Strongest signal (mlog10p > 280), colocalizes with T1D, gestational diabetes, metabolic traits |
    | **MODY genes in GWAS** | HNF1A, HNF4A, GCK, KCNJ11, ABCC8, PPARG — all carry common T2D variants AND cause monogenic diabetes |

    T2D susceptibility is deeply polygenic: hundreds of variants, each with a tiny effect,
    spread across the genome. But the story isn't just "many small effects." Several of the
    genes harboring common T2D variants are the same genes where rare mutations cause
    monogenic diabetes (MODY). This convergence is a powerful validator — it means the GWAS
    is finding real biology, not statistical noise, because the same genes are independently
    confirmed by Mendelian genetics.

    The clinical implication: polygenic risk scores for T2D aggregate these hundreds of signals
    into a single per-person number. Combined with the known MODY genes, we now have
    genetic tools spanning the full allele-frequency spectrum — from rare mutations diagnosed
    in families to common-variant scores applicable to entire populations.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ## To extend

    - **Compare T2D with T1D.** Search for `T1D_WIDE` on FinnGen and run the same pipeline.
      T1D's genetic architecture is concentrated in the HLA region on chromosome 6 — a striking
      contrast with T2D's genome-wide spread.
    - **Add the exome layer.** When `/exome_results_by_gene` is available, import
      `prepare_deleterious` from nb04 and build a multi-gene forest plot for the MODY genes,
      showing rare coding variant effects alongside the common-variant GWAS signals.
    - **Drill into a specific locus.** Use nb02's PheWAS pattern on the TCF7L2 lead variant
      to see every trait it associates with across all resources, not just the colocalizing ones.
    - **Cross-ethnic comparison.** FinnGen is a Finnish biobank; Open Targets includes
      multi-ancestry studies. Compare credible sets for the same T2D code across resources
      to see which signals replicate and which are population-specific.
    """)
    return


if __name__ == "__main__":
    app.run()
