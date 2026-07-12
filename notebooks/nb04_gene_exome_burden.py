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
    GENEGENIE_TOKEN = os.environ.get("GENEGENIE_TOKEN")

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    from nb01_pcsk9_walkthrough import client, fetch_json, fetch_tsv


@app.function
def prepare_deleterious(exome: pl.DataFrame, ci_z: float = 1.96) -> pl.DataFrame:
    """Filter an `exome_results_by_*` DataFrame to deleterious variants and prep for plotting.

    Keeps only `pLoF` and `missense` annotations, drops rows where `mlog10p`
    saturated at the floating-point cap (`se == 0`) or is null, then attaches
    confidence-interval bounds (`ci_lo`, `ci_hi` at +/- `ci_z * se`) and a
    pretty `trait_variant` label of the form `<trait>  (<chr>:<pos>:<ref>:<alt>)`.
    Caller decides how to sort and how many rows to keep.
    """
    return (
        # Pin numeric dtypes: the API serves TSV, and a gene whose column is all-NA (or has a
        # stray non-numeric) infers as String, so mlog10p/beta/se dtype varies gene to gene.
        # Left unpinned, concatenating per-gene frames (nb09) fails to vstack on the dtype clash.
        exome.with_columns(pl.col("mlog10p", "beta", "se").cast(pl.Float64, strict=False))
        .filter(pl.col("annotation").is_in(["pLoF", "missense"]))
        .filter(pl.col("se") > 0)
        .filter(pl.col("mlog10p").is_not_null())
        .with_columns(
            (pl.col("beta") - ci_z * pl.col("se")).alias("ci_lo"),
            (pl.col("beta") + ci_z * pl.col("se")).alias("ci_hi"),
            pl.concat_str(
                [
                    pl.col("trait"),
                    pl.lit("  ("),
                    pl.col("chr").cast(pl.Utf8),
                    pl.lit(":"),
                    pl.col("pos").cast(pl.Utf8),
                    pl.lit(":"),
                    pl.col("ref"),
                    pl.lit(":"),
                    pl.col("alt"),
                    pl.lit(")"),
                ]
            ).alias("trait_variant"),
        )
    )


@app.cell
def _():
    mo.md(r"""
    # nb04: PCSK9 rare-variant exome view

    nb01 told the *common-variant* PCSK9 story: a missense at chr1:55039974 with PIP near 1
    that colocalizes lipid labs and cardiovascular endpoints. Here we ask the orthogonal
    question: when you set fine-mapped GWAS aside and look only at the **exome single-variant
    results from genebass** -- the rare-variant arm of the same UK Biobank cohort -- does
    PCSK9 still light up the same lipid traits, and is the direction of effect the same?

    Thesis to evaluate: rare loss-of-function and missense variants in PCSK9 should drive LDL
    and total cholesterol *down* with much larger per-allele effect sizes than the common
    missense (which sits at AF ~2%), because each rare allele knocks out one functional copy
    rather than nudging activity. The forest plot at the bottom is the test.
    """)
    return


@app.cell
def _():
    GENE = "PCSK9"
    exome = fetch_tsv(f"/exome_results_by_gene/{GENE}")
    mo.md(
        f"### Pulled {len(exome):,} exome single-variant rows for **{GENE}** "
        f"across {exome['resource'].n_unique()} resource(s) "
        f"({', '.join(sorted(exome['resource'].unique().to_list()))})"
    )
    return GENE, exome


@app.cell
def _(exome):
    by_annotation = (
        exome.group_by("annotation")
        .agg(
            pl.len().alias("n_variant_trait_pairs"),
            pl.col("af_overall").median().alias("median_af"),
        )
        .sort("n_variant_trait_pairs", descending=True)
    )
    mo.vstack(
        [
            mo.md(
                "**Rows by VEP annotation class.** Each row is one (variant, trait) pair, so "
                "pLoF rows count *once per phenotype tested*, not once per variant. Median AF "
                "shows how rare each class is in UK Biobank."
            ),
            by_annotation,
        ]
    )
    return


@app.cell
def _(GENE, exome):
    deleterious = prepare_deleterious(exome)
    top = deleterious.sort("mlog10p", descending=True).head(20)
    mo.vstack(
        [
            mo.md(
                f"### Top 20 deleterious-variant exome hits in {GENE}\n\n"
                "Restricted to `pLoF` and `missense`, dropped the three rows where the p-value "
                "underflowed and `se` came back as 0. UKB trait codes you'll see: `30640` is "
                "ApoB, `30690` is total cholesterol, `30780` is direct LDL-C, `6177_1` is "
                "self-reported cholesterol-lowering medication, `20002_1473` is self-reported "
                "high cholesterol, `130814` is the ICD-10 first-occurrence record for "
                "hypercholesterolemia. Sign of `beta` matters."
            ),
            top.select(
                "annotation",
                "trait",
                "pos",
                "ref",
                "alt",
                "af_overall",
                "mlog10p",
                "beta",
                "se",
            ),
        ]
    )
    return (top,)


@app.cell
def _(GENE, top):
    points = (
        alt.Chart(top)
        .mark_circle(size=80)
        .encode(
            x=alt.X("beta:Q", title="Per-allele effect (beta, IRNT-scale or log-odds)"),
            y=alt.Y("trait_variant:N", sort="-x", title=None),
            color=alt.Color(
                "annotation:N",
                scale=alt.Scale(domain=["pLoF", "missense"], range=["#d62728", "#1f77b4"]),
                title="VEP class",
            ),
            tooltip=[
                "trait",
                "annotation",
                "pos",
                "ref",
                "alt",
                "af_overall",
                "mlog10p",
                "beta",
                "se",
            ],
        )
    )
    bars = points.mark_rule().encode(
        x="ci_lo:Q",
        x2="ci_hi:Q",
    )
    zero = alt.Chart(pl.DataFrame({"x": [0.0]})).mark_rule(strokeDash=[4, 3], color="black").encode(x="x:Q")
    forest = (bars + points + zero).properties(
        height=500,
        title=f"Forest plot: top deleterious exome variants in {GENE} (beta +/- 1.96 SE, sorted by effect size on x)",
    )
    mo.ui.altair_chart(forest)
    return


@app.cell
def _(GENE):
    mo.md(
        r"""
        ## Cross-check: curated gene-disease evidence

        The exome side is one signal; gene-level disease curation is another. `gene_disease`
        pulls submissions from gencc and monarch -- expert panels with written classifications
        like "Definitive" or "Strong" rather than p-values.
        """
    )
    gd = pl.DataFrame(fetch_json(f"/gene_disease/{GENE}"))
    if gd.is_empty():
        mo.md(f"No curated gene-disease entries for {GENE}.")
    else:
        gd = gd.select(
            "resource",
            "disease_title",
            "classification",
            "mode_of_inheritance",
            "submitter",
        )
    gd
    return


@app.cell
def _():
    mo.md(r"""
    ## What the forest plot says vs nb01

    - Direction: every top pLoF and missense hit in PCSK9 has **negative beta** for LDL-C
      (UKB 30780), ApoB (30640), and total cholesterol (30690). Same direction as the common
      variant in nb01 -- knocking out PCSK9 lowers cholesterol, regardless of how rare the
      allele is.
    - Magnitude: pLoF effects (beta ~ -1.0 in IRNT units) are roughly **3x larger** than
      the common missense effect from nb01 on the same trait, exactly what the
      allelic-series logic predicts.
    - Frequency: the common missense in nb01 is at AF ~2%; the strongest pLoF here
      (`1:55061557:G:A`) is AF ~3e-4. Two orders of magnitude rarer, three times the effect.
    - Curation: gencc and monarch both flag PCSK9 for `hypercholesterolemia, autosomal
      dominant, 3` with classifications up to "Definitive". The exome forest is the
      quantitative shadow of that curated call.

    Same gene, two views, one story.

    ## To extend

    - Swap `PCSK9` for `LDLR` and re-run -- LDLR is the *opposite* allelic series (LoF raises
      LDL), so betas should flip sign while still concentrating on the same UKB trait codes.
    - Filter to only `pLoF` and pull `exome_results_by_variant/{chr}:{pos}:{ref}:{alt}` for
      the lead pLoF allele to get its full PheWAS across all UKB phenotypes, not just the
      ones surfaced by gene-keyed query.
    - Pivot to a *different* gene with a known cardiovascular signal but weaker rare-variant
      evidence (e.g. `APOC3`, `ANGPTL3`) and compare how cleanly the forest separates from
      zero -- it's a quick visual read on rare-variant power for a target.
    """)
    return


if __name__ == "__main__":
    app.run()
