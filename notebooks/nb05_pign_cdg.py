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

    from nb01_pcsk9_walkthrough import client, fetch_json, fetch_tsv  # noqa: F401
    from nb04_gene_exome_burden import prepare_deleterious  # noqa: F401


@app.cell
def _():
    mo.md(r"""
    # nb05: PIGN-CDG -- two windows on a recessive Mendelian gene

    PIGN encodes a GPI-anchor biosynthesis enzyme. Biallelic loss-of-function causes
    **MCAHS1** (multiple congenital anomalies-hypotonia-seizures syndrome 1, MONDO:0013563),
    a severe pediatric autosomal-recessive disorder; this is what clinicians call
    "PIGN-CDG" (PIGN-related congenital disorder of glycosylation). The curated arm of
    FinnGenie carries that story cleanly: gencc and monarch return eight independent
    submissions, with ClinGen and G2P both at "Definitive".

    But UK Biobank doesn't really see MCAHS1. It's a recessive pediatric syndrome whose
    affected probands typically don't reach adulthood biobank enrollment, and adult
    heterozygous carriers don't manifest the recessive phenotype. So what does the
    *quantitative* arm (genebass single-variant exome scores in adult UKB heterozygotes)
    say about PIGN, and how does it relate -- if at all -- to the curated MCAHS1 picture?

    Thesis to evaluate: in a recessive Mendelian gene, exome single-variant scores in
    adult heterozygotes will not reproduce the recessive disease phenotype. They will
    instead surface whatever quantitative trait happens to be most sensitive to a
    one-allele dose change in PIGN's biology. The forest plot at the bottom is the
    test, and the answer turns out to be unrelated to MCAHS1.
    """)
    return


@app.cell
def _():
    GENE = "PIGN"
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
                "shows how rare each class is in UK Biobank -- PIGN's pLoF alleles are "
                "extremely rare, consistent with strong purifying selection on a recessive "
                "Mendelian disease gene."
            ),
            by_annotation,
        ]
    )
    return


@app.cell
def _(GENE, exome):
    # Same prep pipeline as nb04 -- filter to deleterious classes, drop mlog10p
    # underflow rows, attach CI bounds and pretty trait_variant labels. Imported
    # rather than copy-pasted so the next time this prep changes, both notebooks
    # pick it up.
    deleterious = prepare_deleterious(exome)
    top = deleterious.sort("mlog10p", descending=True).head(20)
    mo.vstack(
        [
            mo.md(
                f"### Top 20 deleterious-variant exome hits in {GENE}\n\n"
                "Restricted to `pLoF` and `missense`, dropped rows where the p-value "
                "underflowed to the floating-point cap. The trait names look cryptic "
                "because UKB encodes ophthalmology phenotypes as continuous codes: "
                "`Corneal_hysteresis_left/right`, `Corneal_resistance_factor_left/right`, "
                "`IOP_Corneal_compensated_left/right`, `Corneal_curvature` (codes `5085`-"
                "`5263`), and a few visual-acuity fields (`20018`, `20020`). What you do "
                "**not** see in this top-20: hypotonia, seizures, intellectual disability, "
                "or any MCAHS1 phenotype. That's the answer to the thesis -- adult "
                "heterozygotes don't reproduce the recessive disease."
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
    zero = (
        alt.Chart(pl.DataFrame({"x": [0.0]}))
        .mark_rule(strokeDash=[4, 3], color="black")
        .encode(x="x:Q")
    )
    forest = (bars + points + zero).properties(
        height=500,
        title=f"Forest plot: top deleterious exome variants in {GENE} "
        f"(beta +/- 1.96 SE, sorted by effect size on x)",
    )
    forest
    return


@app.cell
def _(GENE):
    mo.md(
        rf"""
        ## Cross-check: curated gene-disease evidence

        The exome forest told one story (corneal biomechanics in heterozygotes); the
        curated arm tells the orthogonal one. `gene_disease` pulls submissions from
        gencc and monarch -- expert panels with written classifications like
        "Definitive" or "Strong" rather than p-values.
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
    ## What the two arms say together

    - **Curated arm:** seven gencc submissions and one monarch entry, all converging on
      MCAHS1 (MONDO:0013563). ClinGen and G2P both at "Definitive"; PanelApp Australia
      and Labcorp at "Strong"; mode of inheritance autosomal recessive across every
      submitter. One Orphanet "Supportive" call for Fryns syndrome -- almost certainly
      noise from overlapping clinical features in the dysmorphology literature.
    - **Quantitative arm:** zero MCAHS1 phenotypes in the top exome hits. Instead, a
      tight cluster of corneal biomechanics traits -- corneal hysteresis, corneal
      resistance factor, IOP-corneal-compensated, corneal curvature -- with rare
      missense variants (AF ~1e-5 to 1e-4) at very large effect sizes. These are the
      adult-heterozygote-carrier signals you can see in UKB.
    - **Why this is the expected pattern.** PIGN is a recessive Mendelian gene. UKB
      doesn't enroll affected MCAHS1 probands. Heterozygous carriers don't manifest
      the recessive disease, but a single missense allele in a GPI-anchor enzyme is
      apparently enough to perturb the GPI-anchored proteins that scaffold corneal
      stroma -- giving you a quantitative phenotype that nobody curates as a "PIGN
      disease" but that genebass picks up cleanly.
    - **Contrast with nb04 (PCSK9).** PCSK9's two arms *converge* on the same
      LDL/cholesterol traits because PCSK9 has both common-variant common-disease
      action (lipid biology in everyone) and rare-variant Mendelian action
      (familial hypercholesterolemia). PIGN's two arms *diverge* because the
      Mendelian disease (MCAHS1) is recessive and pediatric, while the heterozygous
      quantitative trait (corneal biomechanics) is a different layer of biology.
      Both genes have signal in both arms; the *relationship* between the arms is
      what changes.

    Same gene, two views, two stories that don't overlap -- the orthogonal-signal
    pattern is itself the finding.

    ## To extend

    - Pull `exome_results_by_variant/18:62161260:T:C` (or whichever lead missense the
      forest highlights) for the top corneal-biomechanics variant and check its full
      PheWAS across all UKB phenotypes. If the corneal cluster is the *only* signal
      for that variant, it strengthens the "PIGN-heterozygote -> corneal stroma"
      reading.
    - Repeat this notebook for another GPI-anchor biosynthesis gene with a CDG
      diagnosis -- `PIGA`, `PIGN`'s upstream partner `PIGV`, or `PIGT` (also
      MCAHS-spectrum). If they all light up corneal biomechanics in heterozygotes,
      the GPI-anchor pathway has a quantitative ophthalmologic phenotype that the
      curated arm has missed.
    - Pull `credible_sets_by_gene/PIGN` (common-variant arm) and check whether the
      PIGN locus has cis-eQTLs whose tissues match the corneal phenotype -- e.g. is
      PIGN expression QTL'd in fibroblast / connective tissue? That would make the
      heterozygous quantitative signal mechanistically plausible.
    """)
    return


if __name__ == "__main__":
    app.run()
