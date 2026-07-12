# /// script
# requires-python = ">=3.12,<3.14"
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
    import httpx  # noqa: F401
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
def pqtl_credible_sets(variant: str) -> pl.DataFrame:
    """All pQTL credible sets that contain `variant`, with this variant's per-cs stats.

    Hits `/credible_sets_by_variant/{variant}?format=json` and filters to
    `data_type == "pQTL"`. Each row is one credible set (in some pQTL resource,
    for some protein trait) that includes `variant`; the row's `pip`, `beta`,
    and `mlog10p` are for *this* variant's contribution to that credible set,
    not the credible-set lead. Use `pip` to gauge whether the variant is a
    plausible driver of the protein signal vs incidental membership.
    """
    rows = fetch_json(f"/credible_sets_by_variant/{variant}")
    if not rows:
        return pl.DataFrame()
    df = pl.DataFrame(rows)
    if "data_type" not in df.columns:
        return pl.DataFrame()
    return df.filter(pl.col("data_type") == "pQTL")


@app.function
def direction_consensus(df: pl.DataFrame, beta_col: str = "beta") -> dict:
    """Summarize sign agreement across a set of beta estimates.

    Returns a dict with `n`, `n_negative`, `n_positive`, `n_zero`, and a
    `verdict` string -- one of `"all_negative"`, `"all_positive"`, `"mixed"`,
    or `"empty"`. Used to ask whether a variant pushes a panel of proteins
    in the same direction (consistent with a single mechanism) or splits
    them (not a single mechanism, or noise).
    """
    if df.is_empty():
        return {"n": 0, "n_negative": 0, "n_positive": 0, "n_zero": 0, "verdict": "empty"}
    betas = df[beta_col].drop_nulls()
    n_neg = int((betas < 0).sum())
    n_pos = int((betas > 0).sum())
    n_zero = int((betas == 0).sum())
    if n_pos == 0 and n_neg > 0:
        verdict = "all_negative"
    elif n_neg == 0 and n_pos > 0:
        verdict = "all_positive"
    else:
        verdict = "mixed"
    return {
        "n": len(betas),
        "n_negative": n_neg,
        "n_positive": n_pos,
        "n_zero": n_zero,
        "verdict": verdict,
    }


@app.cell
def _():
    mo.md(r"""
    # nb06: Variant -> pQTL credible sets -> mechanism

    The hero example from the GeneGenie demo (Karjalainen, VOA 2026-05-05): take a GWAS-significant
    variant in a gene of interest, ask which **plasma proteins** have a credible set that includes
    this variant, then read the **direction of effect** across those proteins. When a variant sits
    in a sheddase or a regulator and pushes multiple downstream proteins in the same direction,
    that's a far cleaner functional read than any single colocalization test.

    Companion to nb01 (gene-first, GWAS->coloc) and nb02 (variant PheWAS across all data types).
    Where nb02 surveys the full PheWAS dot plot, nb06 zooms into the pQTL bucket and asks:
    *which proteins, in what direction, by how much?* The pQTL filter is what lets you go from
    "this variant is associated with disease" to a candidate molecular mechanism.

    Note: this is **not** formal colocalization. We're using credible-set co-membership plus
    sign-of-effect as a lighter, faster proxy. Use `colocalization_by_variant` (see nb01) when
    you need PP.H4 / CLPP for a regulatory claim.
    """)
    return


@app.cell
def _():
    VARIANT = "2:9521321:A:G"
    GENE_CONTEXT = "ADAM17"
    DISEASE_CONTEXT = "inflammatory bowel disease (IBD)"
    mo.md(
        f"### Entry point: `{VARIANT}` (in {GENE_CONTEXT})\n\n"
        f"A splice-region variant in {GENE_CONTEXT}, the canonical sheddase that cleaves "
        "extracellular domains off membrane-anchored proteins (TNF, EGFR ligands, several "
        f"receptors). Carriers have higher risk of {DISEASE_CONTEXT}. The biological question: "
        f"if {GENE_CONTEXT} activity is altered by this variant, which substrate proteins should "
        "we expect to move, and in what direction?"
    )
    return GENE_CONTEXT, VARIANT


@app.cell
def _(VARIANT):
    pqtl = pqtl_credible_sets(VARIANT)
    if pqtl.is_empty():
        header = mo.md(
            f"No pQTL credible sets contain `{VARIANT}`. Either the variant doesn't index any "
            "protein signal at the resources GeneGenie holds, or the allele encoding differs from "
            f"the credible-set store -- try `nb02.alt_alleles({VARIANT!r})` and re-run with each "
            "candidate."
        )
    else:
        header = mo.md(
            f"### `{VARIANT}` sits in **{len(pqtl):,} pQTL credible sets** "
            f"across {pqtl['resource'].n_unique()} resource(s) "
            f"({', '.join(sorted(pqtl['resource'].unique().to_list()))})."
        )
    header
    return (pqtl,)


@app.cell
def _(pqtl):
    if pqtl.is_empty():
        per_protein_view = mo.md("No pQTL rows -- skipping table.")
    else:
        cols = [
            c
            for c in (
                "trait",
                "resource",
                "cs_id",
                "cs_size",
                "pip",
                "beta",
                "se",
                "mlog10p",
            )
            if c in pqtl.columns
        ]
        ranked = pqtl.sort("pip", descending=True).select(cols)
        per_protein_view = mo.vstack(
            [
                mo.md(
                    "**Per-protein view.** One row per credible set containing this variant. "
                    "`pip` is the variant's posterior inclusion probability *within* the credible "
                    "set -- high pip means the variant is plausibly the driver, low pip means it's "
                    "one of several variants tagging the signal. `beta` is the per-allele effect on "
                    "the protein's plasma level (alt allele relative to ref)."
                ),
                ranked,
            ]
        )
    per_protein_view
    return


@app.cell
def _(pqtl):
    if pqtl.is_empty():
        summary = {"verdict": "empty", "n": 0, "n_negative": 0, "n_positive": 0, "n_zero": 0}
        consensus_view = mo.md("_No pQTL rows; direction consensus not computed._")
    else:
        summary = direction_consensus(pqtl)
        readout = {
            "all_negative": "**all proteins move down** with the alt allele -- consistent with "
            "loss of a positive regulator (or gain of an inhibitor).",
            "all_positive": "**all proteins move up** with the alt allele -- consistent with "
            "gain of a positive regulator (or loss of an inhibitor).",
            "mixed": "proteins move in **mixed directions** -- not a single shared mechanism, "
            "or the panel is heterogeneous.",
            "empty": "no rows to summarize.",
        }[summary["verdict"]]
        consensus_view = mo.md(
            f"### Direction consensus across {summary['n']} pQTL hits\n\n"
            f"- {summary['n_negative']} negative beta(s), "
            f"{summary['n_positive']} positive, {summary['n_zero']} zero\n"
            f"- Verdict: {readout}"
        )
    consensus_view
    return


@app.cell
def _(VARIANT, pqtl):
    if pqtl.is_empty():
        chart = mo.md("No data to plot.")
    else:
        plot_df = (
            pqtl.sort("pip", descending=True)
            .head(20)
            .with_columns(
                pl.when(pl.col("beta") > 0)
                .then(pl.lit("higher in alt"))
                .otherwise(pl.lit("lower in alt"))
                .alias("direction"),
            )
        )
        if "se" in plot_df.columns:
            plot_df = plot_df.with_columns(
                (pl.col("beta") - 1.96 * pl.col("se")).alias("ci_lo"),
                (pl.col("beta") + 1.96 * pl.col("se")).alias("ci_hi"),
            )
        points = (
            alt.Chart(plot_df)
            .mark_circle(size=140, opacity=0.85)
            .encode(
                x=alt.X("beta:Q", title="Per-allele effect on plasma protein (beta)"),
                y=alt.Y("trait:N", sort="-x", title=None),
                color=alt.Color(
                    "direction:N",
                    scale=alt.Scale(
                        domain=["lower in alt", "higher in alt"],
                        range=["#2b83ba", "#d7191c"],
                    ),
                    title="Effect direction",
                ),
                size=alt.Size("pip:Q", scale=alt.Scale(range=[60, 400]), title="PIP"),
                tooltip=[
                    c
                    for c in (
                        "trait",
                        "resource",
                        "cs_id",
                        "cs_size",
                        "pip",
                        "beta",
                        "se",
                        "mlog10p",
                    )
                    if c in plot_df.columns
                ],
            )
        )
        zero = alt.Chart(pl.DataFrame({"x": [0.0]})).mark_rule(strokeDash=[4, 3], color="black").encode(x="x:Q")
        layers = [points, zero]
        if "ci_lo" in plot_df.columns:
            bars = (
                alt.Chart(plot_df)
                .mark_rule(opacity=0.5)
                .encode(x="ci_lo:Q", x2="ci_hi:Q", y=alt.Y("trait:N", sort="-x"))
            )
            layers.insert(0, bars)
        chart = alt.layer(*layers).properties(
            height=420,
            title=f"pQTL effects at {VARIANT} (top 20 by PIP)",
        )
    mo.ui.altair_chart(chart)
    return


@app.cell
def _(GENE_CONTEXT):
    mo.md(rf"""
    ## Reading the panel

    A variant in **{GENE_CONTEXT}** that pushes a set of plasma proteins in one shared
    direction is the signature of a single upstream mechanism. The interpretation step --
    which the GeneGenie chatbot does via literature search -- is to ask whether the affected
    proteins share a known relationship to {GENE_CONTEXT}'s function. For ADAM17 specifically,
    the demo's reading was: ERBB4 and TNFRSF1A are documented ADAM17 substrates whose
    ectodomains are cleaved off the cell surface; if the variant **reduces** ADAM17 sheddase
    activity, fewer ectodomains are released into plasma, so plasma levels of the cleaved
    forms drop -- consistent with the all-negative beta panel above.

    That last step is not in this notebook on purpose. Substrate-level interpretation is a
    literature task, not an API call -- if you want it agentic, point an LLM at the protein
    list with `{GENE_CONTEXT}` as context and ask which are known substrates / regulated
    partners. The notebook's job is to surface the panel and the direction; the biology call
    is the reader's.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ## To extend

    - Swap `2:9521321:A:G` for a variant in **PCSK9** (e.g. `1:55039974:G:T`, the LoF missense
      from nb01) and re-run -- you should see plasma PCSK9 itself drop and downstream LDL-pathway
      proteins shift, a different mechanism flavor than ADAM17's substrate panel.
    - Take the top-PIP pQTL credible set above, pull its `cs_id`, and feed it into
      `/colocalization_by_credible_set_id/{cs_id}` to get the formal coloc verdict
      (PP.H4 / CLPP) for the protein vs the GWAS trait that motivated the lookup. That bridges
      this view back to nb01's coloc step.
    - Pivot to the **eQTL** bucket (filter `data_type == "eQTL"` instead of `"pQTL"`) and ask
      whether the same variant moves the gene's *transcript* in the same direction as it moves
      the protein. Disagreement is informative: post-translational regulation, splicing, or
      protein stability rather than transcription.
    """)
    return


if __name__ == "__main__":
    app.run()
