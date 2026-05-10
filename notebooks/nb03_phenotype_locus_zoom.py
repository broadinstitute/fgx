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
    import httpx
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


@app.function
def pick_leads(cs: pl.DataFrame) -> pl.DataFrame:
    """One row per credible set: the variant with the highest -log10p inside each `cs_id`.

    Adds a `variant_id` column in `chr:pos:ref:alt` form. Sorts the result by
    `mlog10p` descending so `.head(N)` gives the top loci. Generic for any
    credible-set DataFrame returned by `credible_sets_by_*` (same schema across
    resources).
    """
    return (
        cs.filter(pl.col("mlog10p").is_not_null())
        .sort("mlog10p", descending=True)
        .group_by("cs_id", maintain_order=True)
        .head(1)
        .with_columns(
            pl.concat_str(
                [pl.col("chr"), pl.col("pos"), pl.col("ref"), pl.col("alt")],
                separator=":",
            ).alias("variant_id")
        )
        .sort("mlog10p", descending=True)
    )


@app.function
def annotate_with_nearest_gene(variants: list[str]) -> pl.DataFrame:
    """For each `chr:pos:ref:alt` variant, fetch its nearest protein-coding gene.

    Calls `/nearest_genes/{variant}?n=1` per variant; returns a DataFrame with
    columns `variant_id`, `nearest_gene`, `distance`. Swallows
    `httpx.HTTPStatusError` per variant and emits null annotations for those rows
    so the caller can still join cleanly.
    """
    rows = []
    for v in variants:
        try:
            g = fetch_json(f"/nearest_genes/{v}", n=1)
            rows.append(
                {
                    "variant_id": v,
                    "nearest_gene": g[0]["gene_name"] if g else None,
                    "distance": g[0]["distance"] if g else None,
                }
            )
        except httpx.HTTPStatusError:
            rows.append(
                {"variant_id": v, "nearest_gene": None, "distance": None}
            )
    return pl.DataFrame(rows)


@app.cell
def _():
    mo.md(r"""
    # nb03: phenotype-first locus discovery

    Pick a phenotype, ask FinnGenie for every fine-mapped credible set genome-wide, then draw a
    Manhattan from the lead variants and annotate the top loci with their nearest gene. The
    entry point flips relative to nb01: instead of "what does this gene do?", we're asking
    "where in the genome does this disease live?".

    Default phenotype is `I9_CHD` on FinnGen R13 -- coronary heart disease, well-powered, plenty
    of independent loci to make a Manhattan worth drawing.
    """)
    return


@app.cell
def _():
    RESOURCE = "finngen"
    PHENOTYPE = "I9_CHD"
    cs = fetch_tsv(f"/credible_sets_by_phenotype/{RESOURCE}/{PHENOTYPE}")
    n_cs = cs["cs_id"].n_unique()
    mo.md(
        f"### Pulled {len(cs):,} credible-set rows for **{PHENOTYPE}** "
        f"on **{RESOURCE}**, spanning **{n_cs}** independent credible sets"
    )
    return PHENOTYPE, cs


@app.cell
def _(cs):
    # One row per credible set: the variant with the highest -log10p inside each cs_id.
    # That variant is the locus's "lead" -- what gets a Manhattan dot and a gene label.
    leads = pick_leads(cs)
    mo.vstack(
        [
            mo.md(
                f"### {len(leads)} lead variants (one per credible set)\n\n"
                "Top 10 by `-log10p`. `gene_most_severe` is the gene the lead variant's "
                "consequence falls in -- not necessarily the causal gene, but a useful prior."
            ),
            leads.head(10).select(
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
    return (leads,)


@app.cell
def _(leads):
    # Annotate the top 10 loci with their nearest gene from the gene-model endpoint.
    # gene_most_severe is variant-consequence-based; nearest_genes is purely positional and
    # works even when the lead falls in an intergenic region.
    top10 = leads.head(10)
    nearest = annotate_with_nearest_gene(top10["variant_id"].to_list())
    top10_annotated = top10.join(nearest, on="variant_id", how="left").select(
        "variant_id",
        "mlog10p",
        "pip",
        "most_severe",
        "gene_most_severe",
        "nearest_gene",
        "distance",
    )
    mo.vstack(
        [
            mo.md(
                "### Top 10 loci with nearest-gene annotation\n\n"
                "`gene_most_severe` and `nearest_gene` usually agree -- when they don't, the lead "
                "is intergenic and the nearest gene is the better label."
            ),
            top10_annotated,
        ]
    )
    return top10, top10_annotated


@app.cell
def _(leads, top10_annotated):
    # Manhattan plot. x = cumulative position across chromosomes, y = -log10p.
    # Build chromosome offsets so positions stack left-to-right without gaps.
    chrom_order = [str(i) for i in range(1, 23)] + ["23", "X", "Y"]
    leads_str = leads.with_columns(pl.col("chr").cast(pl.Utf8).alias("chr_str"))
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
    # Cumulative offset = sum of all prior chromosome lengths.
    offsets = chr_max.with_columns(
        pl.col("chr_len").cum_sum().shift(1, fill_value=0).alias("offset")
    ).select("chr_str", "chr_rank", "offset")

    manhattan_df = (
        leads_str.join(offsets, on="chr_str", how="left")
        .with_columns((pl.col("pos") + pl.col("offset")).alias("cum_pos"))
        .with_columns((pl.col("chr_rank") % 2).cast(pl.Utf8).alias("band"))
    )

    # Label = top-10 leads, joined with their nearest_gene.
    label_df = (
        manhattan_df.join(
            top10_annotated.select("variant_id", "nearest_gene"),
            on="variant_id",
            how="inner",
        )
        .with_columns(
            pl.coalesce([pl.col("nearest_gene"), pl.col("gene_most_severe")]).alias(
                "label"
            )
        )
        .filter(pl.col("label").is_not_null())
    )
    return label_df, manhattan_df


@app.cell
def _(PHENOTYPE, label_df, manhattan_df):
    points = (
        alt.Chart(manhattan_df)
        .mark_circle(size=55, opacity=0.85)
        .encode(
            x=alt.X("cum_pos:Q", title="Genome position", axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("mlog10p:Q", title="-log10(p)"),
            color=alt.Color(
                "band:N",
                scale=alt.Scale(domain=["0", "1"], range=["#1f77b4", "#9ecae1"]),
                legend=None,
            ),
            tooltip=[
                "variant_id",
                "chr_str",
                "pos",
                "mlog10p",
                "pip",
                "gene_most_severe",
                "most_severe",
            ],
        )
    )
    labels = (
        alt.Chart(label_df)
        .mark_text(align="left", dx=6, dy=-4, fontSize=11, fontWeight="bold")
        .encode(x="cum_pos:Q", y="mlog10p:Q", text="label:N")
    )
    chart = (points + labels).properties(
        height=380,
        width=820,
        title=f"Manhattan of credible-set leads for {PHENOTYPE} (FinnGen R13)",
    )
    chart
    return


@app.cell
def _(top10):
    # Drill into the strongest locus: pull every coloc pair sharing its lead variant
    # and show the ten with the highest H4 (probability of one shared causal signal).
    top_variant = top10["variant_id"][0]
    coloc = pl.DataFrame(fetch_json(f"/colocalization_by_variant/{top_variant}"))
    if len(coloc) == 0:
        mo.md(f"### No colocalization pairs found at `{top_variant}`")
    else:
        top_coloc = (
            coloc.filter(pl.col("trait1") != pl.col("trait2"))
            .sort("PP.H4.abf", descending=True)
            .head(10)
            .select(
                "trait1",
                "dataset1",
                "trait2",
                "dataset2",
                pl.col("PP.H4.abf").round(3).alias("PP_H4"),
                pl.col("clpp").round(4).alias("clpp"),
                "hit1_mlog10p",
                "hit2_mlog10p",
            )
        )
        mo.vstack(
            [
                mo.md(
                    f"### Strongest coloc pairs at the top locus `{top_variant}`\n\n"
                    "Sorted by H4 (posterior that the two traits share a single causal variant). "
                    "This tells you which other phenotypes ride the same signal."
                ),
                top_coloc,
            ]
        )
    return


@app.cell
def _():
    mo.md(r"""
    ## To extend

    - swap `I9_CHD` for `E4_DM2` (type 2 diabetes) or `J10_ASTHMA` and re-run -- the Manhattan
      shape will tell you immediately whether the disease is locus-rich or polygenic-flat.
    - flip from `finngen` to a different resource (`ukbb`, `ebi`) on the same phenocode where
      available, and overlay both Manhattans to see which loci replicate.
    - swap the Manhattan for a single-locus zoom: filter `cs` to the top `cs_id`, plot
      `pos` vs `mlog10p` for every variant in that credible set, and color by `pip`.
    """)
    return


if __name__ == "__main__":
    app.run()
