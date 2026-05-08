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

    import altair as alt  # noqa: F401
    import httpx  # noqa: F401
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


@app.function
def list_datasets() -> pl.DataFrame:
    """All datasets FinnGenie indexes, flattened to one row per dataset.

    Hits `/datasets` and pulls the top-level fields out of the nested
    `products` (`credible_sets`, `summary_stats`, `colocalization.partners`)
    and `stats` (`n_phenotypes`, `n_samples_median`) substructures into flat
    columns suitable for sorting, filtering, and grouping. Unknown extras
    are dropped on purpose -- print the raw response with
    `fetch_json("/datasets")` if you need them.

    Returned columns: `dataset_id`, `resource`, `version`, `data_type`,
    `trait_type`, `n_phenotypes`, `n_samples`, `has_credible_sets`,
    `has_summary_stats`, `n_coloc_partners`, `description`, `author`,
    `publication_date`.
    """
    rows = fetch_json("/datasets")
    if not rows:
        return pl.DataFrame()

    flat = []
    for r in rows:
        products = r.get("products") or {}
        stats = r.get("stats") or {}
        coloc = products.get("colocalization") or {}
        flat.append({
            "dataset_id": r.get("dataset_id"),
            "resource": r.get("resource"),
            "version": r.get("version"),
            "data_type": r.get("data_type"),
            "trait_type": r.get("trait_type"),
            "n_phenotypes": stats.get("n_phenotypes"),
            "n_samples": r.get("n_samples") or stats.get("n_samples_median"),
            "has_credible_sets": bool(products.get("credible_sets")),
            "has_summary_stats": bool(products.get("summary_stats")),
            "n_coloc_partners": len(coloc.get("partners") or []) if isinstance(coloc, dict) else 0,
            "description": r.get("description"),
            "author": r.get("author"),
            "publication_date": r.get("publication_date"),
        })
    return pl.DataFrame(flat)


@app.function
def list_resources() -> pl.DataFrame:
    """All resource entries grouped by product family, flattened.

    Hits `/resources`, which returns a dict keyed by family
    (`credible_sets`, `colocalization`, `expression`, `chromatin_peaks`,
    `exome_results`, `gene_based`, `gene_disease`). Each family maps to a
    list of resource entries with `id`, `resource`, `gencode_version`, and a
    nested `metadata` block. This helper unrolls the dict into one row per
    (family, resource entry) with `family`, `id`, `resource`,
    `gencode_version`, `version_label`, `author`, `publication_date`.
    """
    payload = fetch_json("/resources")
    if not payload:
        return pl.DataFrame()
    if isinstance(payload, list):
        return pl.DataFrame(payload)

    flat = []
    for family, entries in payload.items():
        for e in entries or []:
            if isinstance(e, str):
                flat.append({
                    "family": family, "id": e, "resource": e,
                    "gencode_version": None, "version_label": None,
                    "author": None, "publication_date": None,
                })
                continue
            meta = e.get("metadata") or {}
            flat.append({
                "family": family,
                "id": e.get("id") or e.get("name"),
                "resource": e.get("resource") or e.get("name"),
                "gencode_version": e.get("gencode_version"),
                "version_label": meta.get("version_label") or e.get("version"),
                "author": meta.get("author"),
                "publication_date": meta.get("publication_date"),
            })
    return pl.DataFrame(flat) if flat else pl.DataFrame()


@app.function
def resource_metadata(resource: str) -> pl.DataFrame:
    """Per-phenotype metadata for one resource.

    Hits `/resource_metadata/{resource}`. Response is large for high-phenotype
    resources -- finngen R13 returns ~2,754 rows, open_targets ~19,691.
    Filter or `.head(N)` before printing; never let the raw response land
    in the chat at full size.

    Typical columns (finngen): `phenotype_code`, `phenotype_string`,
    `n_samples`, `n_cases`, `n_controls`, `trait_type`, `resource`, `version`.
    """
    rows = fetch_json(f"/resource_metadata/{resource}")
    if not rows:
        return pl.DataFrame()
    return pl.DataFrame(rows, infer_schema_length=None)


@app.cell
def _():
    mo.md(r"""
    # nb07: What's in the catalog?

    The FinnGenie demo (Karjalainen, VOA 2026-05-05) opened with a slide cataloguing all 29
    datasets the API exposes -- GWAS, QTL families, exome / gene-based, expression, gene-disease,
    chromatin peaks, colocalization-only -- with phenotype counts and sample sizes. This notebook
    is the agentic equivalent: hit `/datasets` and `/resources`, group the response, and let an
    agent (or a human) introspect what's available before composing an analysis.

    Use this as the **first** notebook to read when a new question lands. "Can FinnGenie answer
    X for trait Y?" reduces to "is there a dataset whose `resource`/`type` covers Y?" -- and the
    answer is one cell away. The dataset names you find here are the strings you pass to the
    `resource=` parameter on most other endpoints (nb01-nb04 helpers).
    """)
    return


@app.cell
def _():
    datasets = list_datasets()
    if datasets.is_empty():
        header = mo.md("`/datasets` returned no rows -- check `FINNGENIE_TOKEN` and `BASE`.")
    else:
        header = mo.md(
            f"### `/datasets` returned {len(datasets):,} rows "
            f"with columns: `{', '.join(datasets.columns)}`"
        )
    header
    return (datasets,)


@app.cell
def _(datasets):
    if datasets.is_empty():
        category_view = mo.md("No datasets to summarize.")
    else:
        group_col = next(
            (c for c in ("type", "data_type", "category") if c in datasets.columns),
            None,
        )
        if group_col is None:
            summary = datasets
            note = (
                "No obvious category column (`type` / `data_type` / `category`) -- "
                "showing the raw catalog. Adjust `group_col` once the schema is known."
            )
        else:
            agg_exprs = [pl.len().alias("n_datasets")]
            if "resource" in datasets.columns:
                agg_exprs.append(pl.col("resource").n_unique().alias("n_resources"))
            summary = (
                datasets.group_by(group_col)
                .agg(agg_exprs)
                .sort("n_datasets", descending=True)
            )
            note = (
                "Maps the slide's category bands (GWAS / QTL / Exome / Expression / "
                "Gene-disease / Chromatin peaks / Colocalization-only) onto whatever "
                "category column the API actually exposes."
            )
        category_view = mo.vstack([
            mo.md(f"### Datasets grouped by type\n\n{note}"),
            summary,
        ])
    category_view
    return


@app.cell
def _(datasets):
    if datasets.is_empty():
        catalog_table = None
    else:
        cols = [c for c in (
            "dataset_id", "resource", "version", "data_type", "trait_type",
            "n_phenotypes", "n_samples", "has_credible_sets", "has_summary_stats",
            "n_coloc_partners",
        ) if c in datasets.columns]
        catalog_view = datasets.select(cols) if cols else datasets
        catalog_table = mo.ui.table(catalog_view, selection="single", page_size=20)
    catalog_table
    return (catalog_table,)


@app.cell
def _(catalog_table, datasets):
    if catalog_table is None or datasets.is_empty():
        selected_resource = None
        selection_note = mo.md(
            "_Select a row above to see per-phenotype metadata for that resource._"
        )
    else:
        sel = catalog_table.value
        if sel is None or sel.is_empty():
            selected_resource = None
            selection_note = mo.md(
                "_Select a row above to see per-phenotype metadata for that resource._"
            )
        else:
            row = sel.row(0, named=True)
            selected_resource = row.get("resource") or row.get("dataset_id")
            selection_note = mo.md(f"### Drilling into `{selected_resource}`")
    selection_note
    return (selected_resource,)


@app.cell
def _(selected_resource):
    if selected_resource is None:
        meta = pl.DataFrame()
        meta_header = mo.md("")
    else:
        meta = resource_metadata(selected_resource)
        if meta.is_empty():
            meta_header = mo.md(f"No metadata returned for `{selected_resource}`.")
        else:
            meta_header = mo.md(
                f"`/resource_metadata/{selected_resource}` returned {len(meta):,} rows. "
                "Showing the first 20 -- use `meta.filter(...)` to narrow."
            )
    meta_header
    return (meta,)


@app.cell
def _(meta):
    if meta.is_empty():
        meta_view = mo.md("_No rows to display._")
    else:
        meta_view = meta.head(20)
    meta_view
    return


@app.cell
def _():
    mo.md(r"""
    ## Phenotype counts at a glance

    The other introspection question the FinnGenie demo answered live was *"how many phenotypes
    are in FinnGen R13?"* -- 2,754 core GWAS, 388 Kanta lab, 126 drug-purchase. That number lives
    inside `stats.n_phenotypes` on each `/datasets` row, which `list_datasets()` has already
    pulled out into the flat `n_phenotypes` column. The cell below filters to the FinnGen family;
    swap the filter to ask the same question for ukbb / open_targets / genebass.
    """)
    return


@app.cell
def _(datasets):
    if datasets.is_empty() or "n_phenotypes" not in datasets.columns:
        finngen_view = mo.md(
            "_No `n_phenotypes` column on `/datasets` -- skip this cell or update the schema._"
        )
    else:
        if "resource" in datasets.columns:
            finngen = datasets.filter(
                pl.col("resource").str.contains("finngen", literal=False)
            )
        else:
            finngen = datasets
        finngen = finngen.select(
            [c for c in ("dataset_id", "resource", "version", "n_phenotypes", "n_samples")
             if c in datasets.columns]
        ).sort("n_phenotypes", descending=True, nulls_last=True)
        finngen_view = mo.vstack([
            mo.md("### FinnGen-family phenotype counts"),
            finngen,
        ])
    finngen_view
    return


@app.cell
def _():
    mo.md(r"""
    ## To extend

    - Add a `mo.ui.dropdown` over the category column, wire its `.value` into the catalog table's
      filter, and you've built a one-cell catalog browser an agent can drive in a single step.
    - Pull `/trait_name_mapping` once and join it to the metadata table so phenotype codes
      (`I9_CHD`, `T2D`, `K11_CROHN`, ...) resolve to human labels in place.
    - Cache `/datasets` to a local parquet (it's small and changes only on FinnGenie releases).
      An agent can then introspect the catalog without hitting the API on every cold start.
    - When the OpenAPI spec adds a new endpoint family, the surface this notebook surveys still
      stays one cell -- the spec at `/api/v1/openapi.json` is the source of truth, and adding a
      `list_<thing>()` helper keeps the introspection pattern uniform.
    """)
    return


if __name__ == "__main__":
    app.run()
