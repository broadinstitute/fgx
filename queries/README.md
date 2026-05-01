# fgx queries

Self-contained [ggsql](https://ggsql.org) charts against a local DuckDB cached from FinnGenie. Mirrors jx's `queries/` pattern: each `q*.gsql` is one chart with no Python in the loop.

## How it composes with the rest of fgx

```
recipes/q*.sh        ->  one-shot curl + jq/duckdb on stdout (humans/agents discovering the API)
queries/q*.gsql      ->  declarative SQL + Vega-Lite spec (this dir)
notebooks/nb*.py     ->  marimo, composes multiple endpoints with reactive plots
```

The auth boundary lives in `just fetch`. The ggsql files run against a local cached DuckDB and don't see the bearer token.

## Prerequisites

- [`ggsql`](https://ggsql.org/get_started/installation.html) - runs each `.gsql` and emits a Vega-Lite JSON.
- [`duckdb`](https://duckdb.org/docs/installation/) - used by `just build-db` to cache TSVs into `data/cache.duckdb`.
- A populated `.env` with `FINNGENIE_TOKEN` (used only by `just fetch`).

## Running

From the fgx repo root:

```bash
just fetch       # curl TSVs from finngenie.fi/api/v1/* into queries/data/
just build-db    # load TSVs into queries/data/cache.duckdb
just render      # run every q*.gsql -> queries/rendered/q*.json (Vega-Lite spec)
```

Paste a rendered JSON into [vega.github.io/editor](https://vega.github.io/editor) for the SVG. (jx's `render.py` regenerates a README catalog with embedded SVGs via `vl-convert`; fgx skips that for now -- one less moving part.)

## Catalog

| # | File | Question |
|---|---|---|
| q01 | [`q01_top_traits_near_pcsk9.gsql`](q01_top_traits_near_pcsk9.gsql) | Top 20 GWAS associations within the PCSK9 credible-set window, ranked by `mlog10p`. |
| q02 | [`q02_top_chd_loci.gsql`](q02_top_chd_loci.gsql) | Top 20 fine-mapped loci for `I9_CHD` (FinnGen R13), annotated with nearest gene + consequence. |

## Adding a query

1. If your question needs a TSV that isn't cached yet, add a `curl` line to the `fetch:` recipe in the root `Justfile`.
2. Add a matching `CREATE TABLE` to `build-db:`.
3. Drop a new `q??_*.gsql` in this dir referencing the new table.
4. Append a row to the catalog table above.

The data layer is stupidly simple by design: TSVs in, DuckDB tables out, declarative SQL on top.
