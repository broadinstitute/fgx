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

__generated_with = "0.10.0"
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
    mo.md(
        r"""
        # nb01: PCSK9 walkthrough

        Pick a gene, find its lead missense variant in fine-mapped GWAS, then ask which traits share
        its causal signal. Demonstrates the full pattern fgx is built on: an `httpx.get` against
        `https://finngenie.fi/api/v1/*`, parse the response, plot.
        """
    )
    return


@app.cell
def _():
    mo.md(
        r"""
        ## How fgx talks to FinnGenie

        The whole "library" is `httpx.get` against `https://finngenie.fi/api/v1/*` with a bearer token.
        Every Python call below is the equivalent of one of these shell one-liners; if you'd rather
        skip Python, the same data comes out either way:

        ```bash
        # list all 29 datasets
        curl -H "Authorization: Bearer $FINNGENIE_TOKEN" https://finngenie.fi/api/v1/datasets

        # credible sets near a gene (default Content-Type is text/tab-separated-values; pipe to duckdb)
        curl -H "Authorization: Bearer $FINNGENIE_TOKEN" \
             "https://finngenie.fi/api/v1/credible_sets_by_gene/PCSK9"

        # colocalization at a variant (JSON, opt-in via ?format=json)
        curl -H "Authorization: Bearer $FINNGENIE_TOKEN" \
             "https://finngenie.fi/api/v1/colocalization_by_variant/1:55039974:G:T?format=json"
        ```

        That's the entire data-access layer. The two helpers below (`fetch_tsv`, `fetch_json`) are
        thin wrappers that fold the bearer header in, parse the response, and hand back a polars
        DataFrame or a list of dicts. There is no SDK, no MCP server, no schema cache.

        For the **complete surface** -- 26 paths across 13 tags (credible sets, colocalization,
        exome, summary stats, search, rsid, etc.) -- the API ships its own OpenAPI 3.1 spec:

        - <https://finngenie.fi/api/v1/openapi.json> (machine-readable, ~74 KB)
        - <https://finngenie.fi/api/v1/docs> (Swagger UI -- browseable; live + version-correct)

        When composing a new notebook against an endpoint nb01 doesn't already touch, read the
        spec rather than guessing the path; both links above are the source of truth.
        """
    )
    return


@app.cell
def _():
    GENE = "PCSK9"
    cs = fetch_tsv(f"/credible_sets_by_gene/{GENE}")
    mo.md(f"### Pulled {len(cs):,} credible-set rows for **{GENE}** across all data types")
    return GENE, cs


@app.cell
def _(cs):
    by_type = (
        cs.group_by("data_type")
        .agg(pl.len().alias("n_rows"))
        .sort("n_rows", descending=True)
    )
    mo.vstack([mo.md("**Rows by data type:**"), by_type])
    return (by_type,)


@app.cell
def _(GENE, cs):
    top_gwas = (
        cs.filter(pl.col("data_type") == "GWAS")
        .filter(pl.col("mlog10p").is_not_null())
        .sort("mlog10p", descending=True)
        .head(20)
        .select(
            "trait",
            "resource",
            "dataset",
            "chr",
            "pos",
            "ref",
            "alt",
            "mlog10p",
            "beta",
            "pip",
            "most_severe",
        )
    )
    mo.vstack([mo.md(f"### Top 20 GWAS hits near {GENE}"), top_gwas])
    return (top_gwas,)


@app.cell
def _(GENE, top_gwas):
    chart = (
        alt.Chart(top_gwas)
        .mark_bar()
        .encode(
            x=alt.X("mlog10p:Q", title="-log10(p)"),
            y=alt.Y("trait:N", sort="-x", title=None),
            color=alt.Color("pip:Q", scale=alt.Scale(scheme="viridis"), title="PIP"),
            tooltip=["trait", "resource", "dataset", "mlog10p", "beta", "pip"],
        )
        .properties(height=400, title=f"Top traits associated near {GENE}")
    )
    chart
    return (chart,)


@app.cell
def _(cs):
    mo.md(
        r"""
        ## Lead missense variant

        PCSK9's canonical hit is the loss-of-function missense at chr1:55039974. Pull the highest-PIP
        missense variant from the GWAS bucket -- that's the variant we'll feed into colocalization.
        """
    )
    lead = (
        cs.filter(pl.col("data_type") == "GWAS")
        .filter(pl.col("most_severe").str.contains("missense"))
        .filter(pl.col("mlog10p").is_not_null())
        .sort("mlog10p", descending=True)
        .head(1)
    )
    variant_id = (
        f"{lead['chr'][0]}:{lead['pos'][0]}:{lead['ref'][0]}:{lead['alt'][0]}"
    )
    mo.md(
        f"**Lead missense variant:** `{variant_id}` "
        f"(PIP={lead['pip'][0]:.3f}, -log10p={lead['mlog10p'][0]:.0f}, "
        f"trait={lead['trait'][0]}, resource={lead['resource'][0]})"
    )
    return lead, variant_id


@app.cell
def _(variant_id):
    coloc = pl.DataFrame(fetch_json(f"/colocalization_by_variant/{variant_id}"))
    mo.md(
        f"### {len(coloc):,} colocalization pairs at `{variant_id}`"
        " (across all dataset combinations)"
    )
    return (coloc,)


@app.cell
def _(coloc):
    top_coloc = (
        coloc.filter(pl.col("hit1_mlog10p").is_not_null())
        .sort("hit1_mlog10p", descending=True)
        .head(25)
        .select(
            "trait1",
            "dataset1",
            "trait2",
            "dataset2",
            "hit1_mlog10p",
            "hit2_mlog10p",
            "hit1_beta",
            "hit2_beta",
        )
    )
    mo.vstack(
        [
            mo.md(
                "### Top colocalizing trait pairs (by `hit1_mlog10p`)\n\n"
                "Each row is one pair of traits whose credible sets share this variant. "
                "Expected biology: lipid labs (LDL, nonHDL) colocalizing with hyperlipidemia "
                "and cardiovascular endpoints, all driven by the PCSK9 LoF missense."
            ),
            top_coloc,
        ]
    )
    return (top_coloc,)


@app.cell
def _(coloc, variant_id):
    pairs = (
        coloc.filter(pl.col("hit1_mlog10p").is_not_null())
        .with_columns(
            pl.concat_str([pl.col("trait1"), pl.lit(" / "), pl.col("trait2")]).alias(
                "pair"
            )
        )
        .sort("hit1_mlog10p", descending=True)
        .head(25)
    )
    coloc_chart = (
        alt.Chart(pairs)
        .mark_bar()
        .encode(
            x=alt.X("hit1_mlog10p:Q", title="-log10(p) for trait1 at this variant"),
            y=alt.Y("pair:N", sort="-x", title=None),
            color=alt.Color("dataset1:N", title="Dataset 1"),
            tooltip=[
                "trait1",
                "trait2",
                "dataset1",
                "dataset2",
                "hit1_mlog10p",
                "hit1_beta",
            ],
        )
        .properties(
            height=500,
            title=f"Top colocalizing trait pairs at {variant_id}",
        )
    )
    coloc_chart
    return (coloc_chart,)


if __name__ == "__main__":
    app.run()
