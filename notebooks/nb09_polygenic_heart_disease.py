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
    from nb03_phenotype_locus_zoom import annotate_with_nearest_gene, pick_leads


@app.cell
def _():
    mo.md(r"""
    # nb09: What is the evidence that heart disease is polygenic?

    A disease is **monogenic** when one gene breaks and the disease follows — cystic fibrosis,
    sickle cell, Huntington's. A disease is **polygenic** when hundreds of genetic variants,
    each with a tiny effect, collectively shift risk across the population. Heart disease is
    the textbook polygenic case.

    This notebook builds that argument from data, using coronary heart disease (`I9_CHD`) on
    FinnGen R13:

    1. **Many independent loci** — a Manhattan plot showing dozens of genome-wide significant
       signals spread across chromosomes, not clustered at one gene.
    2. **Small individual effects** — a histogram showing that most lead variants nudge risk
       by fractions of a percent.
    3. **Diverse biological pathways** — colocalization revealing that CHD loci overlap with
       lipid metabolism, blood pressure, inflammation, and coagulation — not one mechanism.
    4. **Multiple genes in the rare-variant layer** — exome data showing that rare coding
       mutations in several different genes (LDLR, PCSK9, APOB, LPA) independently affect
       cardiovascular risk.

    If heart disease were monogenic, you'd see one locus, one pathway, one gene. Instead,
    you see the opposite at every layer.
    """)
    return


@app.cell
def _():
    RESOURCE = "finngen"
    PHENOTYPE = "I9_CHD"
    chd_cs = fetch_tsv(f"/credible_sets_by_phenotype/{RESOURCE}/{PHENOTYPE}")
    n_credible_sets = chd_cs["cs_id"].n_unique()
    mo.md(
        f"### {len(chd_cs):,} credible-set rows for **{PHENOTYPE}** on **{RESOURCE}**, "
        f"spanning **{n_credible_sets}** independent credible sets"
    )
    return PHENOTYPE, chd_cs


@app.cell
def _():
    mo.md(r"""
    ## 1. Many independent loci across the genome

    If a disease is driven by one gene, fine-mapping finds one or two credible sets on one
    chromosome. Coronary heart disease has dozens of independent credible sets scattered
    across the genome. Each credible set represents a statistically independent signal —
    a distinct region of the genome contributing to disease risk.

    The Manhattan plot below shows one dot per lead variant (the highest -log10p variant in
    each credible set), positioned by genomic coordinate. The horizontal spread IS the
    evidence: heart disease lives everywhere in the genome, not in one place.
    """)
    return


@app.cell
def _(chd_cs):
    chd_leads = pick_leads(chd_cs)
    mo.vstack(
        [
            mo.md(
                f"### {len(chd_leads)} lead variants (one per credible set)\n\n"
                "Each row is the top variant from one independent fine-mapped signal."
            ),
            chd_leads.head(10).select(
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
    return (chd_leads,)


@app.cell
def _(chd_leads):
    chd_top10 = chd_leads.head(10)
    chd_nearest = annotate_with_nearest_gene(chd_top10["variant_id"].to_list())
    chd_top10_annotated = chd_top10.join(chd_nearest, on="variant_id", how="left").select(
        "variant_id",
        "mlog10p",
        "pip",
        "most_severe",
        "gene_most_severe",
        "nearest_gene",
        "distance",
    )
    return chd_top10, chd_top10_annotated


@app.cell
def _(chd_leads, chd_top10_annotated):
    chrom_order = [str(i) for i in range(1, 23)] + ["23", "X", "Y"]
    leads_str = chd_leads.with_columns(pl.col("chr").cast(pl.Utf8).alias("chr_str"))
    chr_max = (
        leads_str.group_by("chr_str")
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
    offsets = chr_max.with_columns(pl.col("chr_len").cum_sum().shift(1, fill_value=0).alias("offset")).select(
        "chr_str", "chr_rank", "offset"
    )
    chd_manhattan_df = (
        leads_str.join(offsets, on="chr_str", how="left")
        .with_columns((pl.col("pos") + pl.col("offset")).alias("cum_pos"))
        .with_columns((pl.col("chr_rank") % 2).cast(pl.Utf8).alias("band"))
    )
    chd_label_df = (
        chd_manhattan_df.join(
            chd_top10_annotated.select("variant_id", "nearest_gene"),
            on="variant_id",
            how="inner",
        )
        .with_columns(pl.coalesce([pl.col("nearest_gene"), pl.col("gene_most_severe")]).alias("label"))
        .filter(pl.col("label").is_not_null())
    )
    return chd_label_df, chd_manhattan_df


@app.cell
def _(PHENOTYPE, chd_label_df, chd_manhattan_df):
    manhattan_points = (
        alt.Chart(chd_manhattan_df.select("cum_pos", "mlog10p", "band"))
        .mark_circle(size=55, opacity=0.85)
        .encode(
            x=alt.X("cum_pos:Q", title="Genome position", axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("mlog10p:Q", title="-log10(p)"),
            color=alt.Color(
                "band:N",
                scale=alt.Scale(domain=["0", "1"], range=["#1f77b4", "#9ecae1"]),
                legend=None,
            ),
        )
    )
    manhattan_labels = (
        alt.Chart(chd_label_df.select("cum_pos", "mlog10p", "label"))
        .mark_text(align="left", dx=6, dy=-4, fontSize=11, fontWeight="bold")
        .encode(x="cum_pos:Q", y="mlog10p:Q", text="label:N")
    )
    manhattan_chart = (manhattan_points + manhattan_labels).properties(
        height=380,
        width=820,
        title=f"Manhattan of credible-set leads for {PHENOTYPE} (FinnGen R13)",
    )
    mo.ui.altair_chart(manhattan_chart)
    return


@app.cell
def _(chd_leads):
    n_chromosomes = chd_leads["chr"].cast(pl.Utf8).n_unique()
    mo.md(
        f"**{len(chd_leads)} independent signals across {n_chromosomes} chromosomes.** "
        "A monogenic disease would show one cluster on one chromosome. Instead, nearly every "
        "chromosome hosts at least one fine-mapped credible set for coronary heart disease. "
        "This genome-wide spread is the first line of evidence for polygenicity."
    )
    return


@app.cell
def _():
    mo.md(r"""
    ## 2. Small individual effects

    The second hallmark of a polygenic disease: each variant contributes very little on its
    own. Monogenic diseases have variants with large effects (odds ratios of 5-50x). Polygenic
    diseases have variants with tiny effects (odds ratios of 1.01-1.15x).

    The histogram below shows the distribution of effect sizes (absolute beta) across all lead
    variants. If heart disease were driven by a few large-effect loci, you'd see a long right
    tail. Instead, the effects cluster near zero — each variant nudges risk by a fraction of
    a standard deviation.
    """)
    return


@app.cell
def _(chd_leads):
    leads_with_beta = chd_leads.filter(pl.col("beta").is_not_null() & pl.col("beta").is_finite()).with_columns(
        pl.col("beta").abs().alias("abs_beta")
    )
    beta_hist = (
        alt.Chart(leads_with_beta.select("abs_beta"))
        .mark_bar(cornerRadiusTopLeft=3, cornerRadiusTopRight=3)
        .encode(
            x=alt.X("abs_beta:Q", bin=alt.Bin(maxbins=30), title="|beta| (per-allele effect size)"),
            y=alt.Y("count()", title="Number of lead variants"),
        )
        .properties(height=300, width=600, title="Distribution of effect sizes across CHD lead variants")
    )
    median_beta = leads_with_beta["abs_beta"].median()
    mo.vstack(
        [
            mo.ui.altair_chart(beta_hist),
            mo.md(
                f"**Median |beta| = {median_beta:.4f}.** Most lead variants shift CHD risk "
                "by less than a tenth of a standard deviation. No single variant comes close "
                "to explaining the disease — but collectively, hundreds of these small nudges "
                "add up to substantial heritable risk."
            ),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    ## 3. Diverse biological pathways

    If heart disease were polygenic but mechanistically simple — say, all loci acted through
    one pathway — it would be "polygenic but mono-mechanistic." The stronger claim is that
    CHD is polygenic AND pleiotropic: the contributing loci span fundamentally different
    biological pathways.

    We can test this with **colocalization**. For each of the top CHD loci, ask: what *other*
    traits share the same causal variant? If the answers span lipids, blood pressure,
    inflammation, and coagulation, the disease is driven by multiple independent biological
    processes.
    """)
    return


@app.cell
def _(chd_top10):
    chd_self_traits = {
        "I9_CHD",
        "I9_IHD",
        "I9_CORATHER",
        "I9_CVD_HARD",
        "I9_REVASC",
        "I9_POSTAMI",
        "I9_MI_STRICT",
        "I9_ANGINA",
        "I9_ANGIO",
        "I9_CABG",
        "I9_UAP",
        "I9_ATHSCLE",
        "I9_HEARTFAIL_AND_CHD",
        "CAD",
    }
    coloc_rows = []
    for v in chd_top10["variant_id"].to_list()[:5]:
        try:
            pairs = pl.DataFrame(fetch_json(f"/colocalization_by_variant/{v}"))
            if len(pairs) > 0:
                diverse_pairs = (
                    pairs.filter(pl.col("trait1") != pl.col("trait2"))
                    .filter(~pl.col("trait2").is_in(chd_self_traits))
                    .filter(pl.col("PP.H4.abf") > 0.5)
                    .sort("PP.H4.abf", descending=True)
                    .unique(subset=["trait2"], keep="first", maintain_order=True)
                    .head(4)
                    .with_columns(pl.lit(v).alias("chd_variant"))
                )
                coloc_rows.append(diverse_pairs)
        except Exception:
            pass

    if coloc_rows:
        coloc_all = pl.concat(coloc_rows, how="diagonal")
    else:
        coloc_all = pl.DataFrame()
    return (coloc_all,)


@app.cell
def _(coloc_all):
    mo.stop(
        coloc_all.is_empty(),
        mo.md("_No colocalization data available for the top CHD loci._"),
    )
    coloc_display = (
        coloc_all.sort("PP.H4.abf", descending=True)
        .head(25)
        .select(
            "chd_variant",
            "trait2",
            "data_type2",
            pl.col("PP.H4.abf").round(3).alias("PP_H4"),
        )
    )
    mo.vstack(
        [
            mo.md(
                "### Cross-pathway colocalizations at the top 5 CHD loci\n\n"
                "CHD-on-CHD self-colocalizations filtered out to surface the diverse biology. "
                "Each row is a non-CHD trait sharing a causal variant with coronary heart disease. "
                "Notice the spread: lipid labs, blood pressure, migraine, stroke, drug prescriptions, "
                "and molecular QTL layers (eQTL, caQTL, metaboQTL) — each pointing at a different "
                "mechanism feeding into the same disease."
            ),
            coloc_display,
        ]
    )
    return


@app.cell
def _(coloc_all):
    mo.stop(
        coloc_all.is_empty(),
        mo.md(""),
    )
    coloc_by_locus = (
        coloc_all.with_columns(
            pl.concat_str(
                [
                    pl.col("trait2"),
                    pl.lit("  ("),
                    pl.col("data_type2"),
                    pl.lit(")"),
                ]
            ).alias("trait_label")
        )
        .sort("PP.H4.abf", descending=True)
        .unique(subset=["trait_label"], keep="first", maintain_order=True)
        .head(20)
    )

    gene_labels = {
        "9:22103814:A:G": "CDKN2B-AS1",
        "6:12903725:A:G": "PHACTR1",
        "6:160564494:G:A": "LPA",
        "1:55039974:G:T": "PCSK9",
        "15:88887075:G:T": "chr15 locus",
    }
    coloc_labeled = coloc_by_locus.with_columns(
        pl.col("chd_variant").replace_strict(gene_labels, default=pl.col("chd_variant")).alias("locus_gene")
    )

    pathway_chart = (
        alt.Chart(coloc_labeled.select("trait_label", "data_type2", "locus_gene"))
        .mark_circle(size=120, opacity=0.9)
        .encode(
            x=alt.X(
                "data_type2:N",
                title="Data type",
                sort=["GWAS", "metaboQTL", "eQTL", "pQTL", "sQTL", "caQTL"],
            ),
            y=alt.Y(
                "trait_label:N",
                title=None,
                sort=alt.EncodingSortField(field="data_type2", order="ascending"),
            ),
            color=alt.Color("locus_gene:N", title="CHD locus"),
        )
        .properties(
            height=480,
            width=550,
            title="Each CHD locus connects to different traits through different mechanisms",
        )
    )
    mo.vstack(
        [
            mo.md(
                "Each dot is a trait that shares a causal variant with a CHD locus (PP_H4 > 0.5, "
                "CHD self-colocs removed). The x-axis groups traits by data type — GWAS phenotypes, "
                "metabolomic QTLs, expression/splicing/chromatin QTLs. The colors spread across columns: "
                "no single locus explains every mechanism, and no single mechanism captures every locus."
            ),
            mo.ui.altair_chart(pathway_chart),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    ## 4. Multiple genes in the rare-variant layer

    Common-variant GWAS is one half of the genetic architecture. The other half is rare
    coding variants — mutations that change or break a protein, carried by <0.1% of people
    but with much larger individual effects.

    If heart disease were monogenic, rare-variant evidence would point at one gene. Instead,
    loss-of-function and missense variants in **multiple genes** independently shift
    cardiovascular risk: LDLR (LDL receptor), PCSK9 (LDL receptor recycling), APOB
    (LDL particle protein), and LPA (lipoprotein(a)). Each gene tells a different mechanistic
    story, and all four now have FDA-approved or late-stage therapies.
    """)
    return


@app.cell
def _():
    from nb04_gene_exome_burden import prepare_deleterious

    cardio_genes = ["LDLR", "PCSK9", "APOB", "LPA"]
    exome_parts = []
    for gene in cardio_genes:
        try:
            raw = fetch_tsv(f"/exome_results_by_gene/{gene}")
            prepped = prepare_deleterious(raw)
            if len(prepped) > 0:
                top_per_gene = prepped.sort("mlog10p", descending=True).head(5).with_columns(pl.lit(gene).alias("gene"))
                exome_parts.append(top_per_gene)
        except Exception:
            pass

    if exome_parts:
        exome_multi = pl.concat(exome_parts, how="diagonal")
    else:
        exome_multi = pl.DataFrame()
    return (exome_multi,)


@app.cell
def _(exome_multi):
    mo.stop(
        exome_multi.is_empty(),
        mo.md("_No exome data available for the selected genes._"),
    )
    forest_df = (
        exome_multi.sort("mlog10p", descending=True)
        .head(20)
        .with_columns(
            pl.concat_str([pl.col("gene"), pl.lit("  "), pl.col("trait_variant")]).alias("gene_trait_variant")
        )
    )
    forest_points = (
        alt.Chart(forest_df.select("gene_trait_variant", "beta", "annotation", "gene"))
        .mark_circle(size=80)
        .encode(
            x=alt.X("beta:Q", title="Per-allele effect (beta)"),
            y=alt.Y("gene_trait_variant:N", sort="-x", title=None),
            color=alt.Color(
                "gene:N",
                title="Gene",
            ),
        )
    )
    forest_bars = (
        alt.Chart(forest_df.select("gene_trait_variant", "ci_lo", "ci_hi"))
        .mark_rule()
        .encode(x="ci_lo:Q", x2="ci_hi:Q", y=alt.Y("gene_trait_variant:N", sort="-x"))
    )
    forest_zero = alt.Chart(pl.DataFrame({"x": [0.0]})).mark_rule(strokeDash=[4, 3], color="black").encode(x="x:Q")
    multi_gene_forest = (forest_bars + forest_points + forest_zero).properties(
        height=500,
        width=650,
        title="Rare coding variants across 4 cardiovascular genes (beta ± 1.96 SE)",
    )
    mo.vstack(
        [
            mo.md(
                "### Forest plot: top rare coding variants across LDLR, PCSK9, APOB, LPA\n\n"
                "Each color is a different gene. The fact that multiple genes carry rare variants "
                "with measurable cardiovascular effects — and that those genes operate in different "
                "parts of lipid biology — is rare-variant evidence for polygenicity."
            ),
            mo.ui.altair_chart(multi_gene_forest),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    ## 5. Synthesis: four layers of polygenic evidence

    | Layer | What a monogenic disease looks like | What CHD looks like |
    |---|---|---|
    | **GWAS loci** | 1-2 signals on one chromosome | Dozens of signals across nearly every chromosome |
    | **Effect sizes** | Large (OR > 5) | Small (median |beta| < 0.1) |
    | **Colocalizing traits** | One trait, one pathway | Lipids, blood pressure, inflammation, coagulation |
    | **Rare-variant genes** | One gene | LDLR, PCSK9, APOB, LPA, and more |

    This is what polygenic means in practice: no single variant, gene, or pathway explains
    the disease. Instead, genetic risk is distributed across hundreds of common variants with
    tiny effects spanning diverse biology, plus rarer coding variants in multiple genes that
    each contribute larger but still incomplete effects.

    The clinical implication is direct: polygenic risk scores (PRS) aggregate all these small
    signals into a single per-person number. People in the top few percent of a CHD PRS have
    heart attack risk comparable to carriers of monogenic familial hypercholesterolemia
    mutations — but the mechanism is fundamentally different. FH is one broken gene; high
    polygenic risk is the unlucky end of a distribution shaped by hundreds of small pushes.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ## To extend

    - **Compare with a truly monogenic disease.** Run nb03 for a Mendelian phenotype (if one
      is available in FinnGen) and see how the Manhattan changes: fewer loci, larger effects,
      concentrated on one chromosome.
    - **Add more genes to the exome panel.** Try `HMGCR` (statin target), `NPC1L1`
      (ezetimibe target), or `ANGPTL3` (emerging target) to see how deep the multi-gene
      architecture goes.
    - **Break out colocalizations by data type.** Filter the coloc table to eQTL or pQTL
      pairs only — these tell you not just *that* a locus is shared, but *how* the variant
      acts (through gene expression or protein levels).
    - **Build a polygenic risk score.** Collect all lead variant betas from the Manhattan
      and weight them into a single score — that's the conceptual foundation of a PRS,
      and the reason polygenic architecture matters clinically.
    """)
    return


if __name__ == "__main__":
    app.run()
