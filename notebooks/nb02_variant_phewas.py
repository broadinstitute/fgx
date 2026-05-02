# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo<0.23.4",
#     "jedi<0.20.0",
#     "polars",
#     "httpx",
#     "altair",
#     "python-dotenv",
#     "pyarrow==24.0.0",
# ]
# ///

import marimo

__generated_with = "0.23.3"
app = marimo.App(width="medium")

with app.setup:
    import io
    import os
    from pathlib import Path

    import altair as alt
    import httpx
    import marimo as mo
    import polars as pl
    from dotenv import load_dotenv

    load_dotenv(Path(__file__).resolve().parent.parent / ".env")
    FINNGENIE_TOKEN = os.environ.get("FINNGENIE_TOKEN")
    BASE = "https://finngenie.fi/api/v1"


@app.function
def client() -> httpx.Client:
    """Authenticated httpx client. Bearer token comes from .env at the repo root."""
    if not FINNGENIE_TOKEN:
        raise RuntimeError(
            "FINNGENIE_TOKEN not set. Copy .env.example to .env and paste your key."
        )
    return httpx.Client(
        headers={"Authorization": f"Bearer {FINNGENIE_TOKEN}"}, timeout=60
    )


@app.function
def fetch_tsv(path: str, **params) -> pl.DataFrame:
    """GET a FinnGenie endpoint as TSV and return a polars DataFrame."""
    with client() as c:
        r = c.get(f"{BASE}{path}", params=params)
        r.raise_for_status()
    return pl.read_csv(io.BytesIO(r.content), separator="\t", null_values="NA")


@app.function
def fetch_json(path: str, **params) -> list[dict]:
    """GET a FinnGenie endpoint as JSON. Most endpoints accept ?format=json."""
    with client() as c:
        r = c.get(f"{BASE}{path}", params={**params, "format": "json"})
        r.raise_for_status()
    return r.json()


@app.cell
def _():
    mo.md(r"""
    # nb02: PheWAS for an rsID

    Start from an rsID, resolve it to a chr:pos:ref:alt variant, then ask: which traits across
    which resources have a credible set that contains this variant? That's a PheWAS view, but
    powered by fine-mapping rather than marginal association — every dot below is a variant
    sitting inside a 95% credible set, so the signal has already survived LD shrinkage.

    Companion to nb01 (gene-first). Same two helpers (`fetch_tsv`, `fetch_json`), different
    entry point and a different chart pattern.
    """)
    return


@app.cell
def _():
    RSID = "rs11591147"
    mo.md(
        f"### Entry point: `{RSID}`\n\n"
        "PCSK9 R46L — a missense loss-of-function allele, one of the cleanest "
        "natural-experiment variants in cardiometabolic genetics. Carriers have lower LDL "
        "cholesterol and lower CHD risk. We expect the PheWAS to light up lipid traits and "
        "atherosclerotic disease."
    )
    return (RSID,)


@app.cell
def _(RSID):
    # Step 1: rsID -> variant. The endpoint returns a list of variants per rsID
    # (multi-allelic sites can map to several).
    resolved = fetch_json("/rsid/variants", rsids=RSID)
    candidate_variants = resolved[0]["variants"] if resolved else []

    # Step 2: not every variant string the resolver returns has credible-set data
    # (allele-strand / ref-allele conventions differ between dbSNP and FinnGen).
    # Try each candidate plus the strand-flipped alt; keep the first that returns rows.
    def alt_alleles(v: str) -> list[str]:
        chrom, pos, ref, alt = v.split(":")
        # FinnGen sometimes stores the same site under a different alt. Try the
        # resolver's answer first, then a small fallback set.
        tries = [v]
        for a in ("A", "C", "G", "T"):
            if a not in (ref, alt):
                tries.append(f"{chrom}:{pos}:{ref}:{a}")
        return tries

    tries: list[str] = []
    for v in candidate_variants:
        tries.extend(alt_alleles(v))

    chosen, cs = None, []
    for v in tries:
        cs = fetch_json(f"/credible_sets_by_variant/{v}")
        if cs:
            chosen = v
            break

    variant_id = chosen
    cs_df = pl.DataFrame(cs) if cs else pl.DataFrame()
    mo.md(
        f"**Resolved {RSID}** -> resolver returned `{candidate_variants}`; "
        f"credible-set data is keyed under `{variant_id}` "
        f"({len(cs_df):,} rows across resources/data types)."
    )
    return cs_df, variant_id


@app.cell
def _(cs_df):
    by_resource = (
        cs_df.group_by(["resource", "data_type"])
        .agg(pl.len().alias("n_rows"))
        .sort("n_rows", descending=True)
    )
    mo.vstack(
        [
            mo.md(
                "**Rows by resource and data type.** GWAS rows dominate; QTL hits "
                "(pQTL/eQTL/caQTL) are the molecular layer that lets us reason about "
                "mechanism if we want to chain into colocalization later."
            ),
            by_resource,
        ]
    )
    return


@app.cell
def _(cs_df, variant_id):
    # Filter to GWAS rows with a finite -log10p, dedupe to one row per (trait, resource)
    # taking the strongest hit, then take the top N for the dot plot.
    phewas = (
        cs_df.filter(pl.col("data_type") == "GWAS")
        .filter(pl.col("mlog10p").is_not_null())
        .sort("mlog10p", descending=True)
        .unique(subset=["trait", "resource"], keep="first")
        .head(25)
        .with_columns(
            pl.when(pl.col("beta") > 0)
            .then(pl.lit("higher in alt"))
            .otherwise(pl.lit("lower in alt"))
            .alias("direction"),
            pl.concat_str(
                [pl.col("trait"), pl.lit(" ("), pl.col("resource"), pl.lit(")")]
            ).alias("trait_label"),
        )
    )
    mo.md(
        f"### Top 25 GWAS traits whose credible set contains `{variant_id}`\n\n"
        "Trait codes are FinnGen / UKBB phenocodes — opaque short strings rather than "
        "human labels. That's authentic to the API; resolving them to plain English would "
        "be its own vignette (see `/phenotype/{resource}/{phenocode}/markdown`)."
    )
    return (phewas,)


@app.cell
def _(phewas, variant_id):
    chart = (
        alt.Chart(phewas)
        .mark_circle(size=140, opacity=0.85)
        .encode(
            x=alt.X("mlog10p:Q", title="-log10(p)"),
            y=alt.Y("trait_label:N", sort="-x", title=None),
            color=alt.Color(
                "direction:N",
                scale=alt.Scale(
                    domain=["lower in alt", "higher in alt"],
                    range=["#2b83ba", "#d7191c"],
                ),
                title="Effect direction",
            ),
            size=alt.Size("pip:Q", scale=alt.Scale(range=[40, 300]), title="PIP"),
            tooltip=[
                "trait",
                "resource",
                "dataset",
                "mlog10p",
                "beta",
                "pip",
                "most_severe",
            ],
        )
        .properties(
            height=520,
            title=f"PheWAS at {variant_id} (credible-set membership)",
        )
    )
    chart
    return


@app.cell
def _(variant_id):
    # Light annotation: which gene is this variant nearest to? Useful when an rsID
    # is unfamiliar and you want a one-line "what region is this".
    near = pl.DataFrame(fetch_json(f"/nearest_genes/{variant_id}", n=3))
    mo.vstack(
        [
            mo.md(f"### Nearest protein-coding genes to `{variant_id}`"),
            near.select(
                "gene_name", "distance", "hgnc_name", "gene_strand", "gene_type"
            ),
        ]
    )
    return


@app.cell
def _():
    mo.md(r"""
    ## To extend

    - Swap `rs11591147` for `rs7412` (APOE, drives LDL via a different mechanism) and re-run —
      the dot plot should still be lipid-anchored but the QTL co-signals will look different.
    - Take the top GWAS trait above and feed its `cs_id` into `/colocalization_by_credible_set_id/{id}`
      to ask which other credible sets share this signal — that bridges PheWAS into the
      colocalization view from nb01.
    - Resolve the FinnGen phenocodes to human labels by calling
      `/phenotype/finngen/{phenocode}/markdown` for each top trait, and re-render the y-axis
      with readable names.
    """)
    return


if __name__ == "__main__":
    app.run()
