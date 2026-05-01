# fgx recipes

Each `q*.sh` is one self-contained question against the FinnGenie REST API. Run any of them with `bash recipes/qNN_*.sh` from the repo root, after `cp .env.example .env` and pasting your token.

## Catalog

| # | Recipe | Endpoint | Question it answers |
|---|---|---|---|
| q01 | [`q01_list_datasets.sh`](q01_list_datasets.sh) | `GET /v1/datasets` | What datasets are available, and which products (credible sets / summary stats / colocalization) does each support? |
| q02 | [`q02_credible_sets_by_gene.sh`](q02_credible_sets_by_gene.sh) | `GET /v1/credible_sets_by_gene/{gene}` | What credible-set fine-mapped associations exist near a gene, across all data types? Walked through with PCSK9 as the worked example. |
| q03 | [`q03_credible_sets_by_phenotype.sh`](q03_credible_sets_by_phenotype.sh) | `GET /v1/credible_sets_by_phenotype` | What genes / loci drive a given phenotype? Worked example: coronary heart disease (`I9_CHD`). |
| q04 | [`q04_colocalization.sh`](q04_colocalization.sh) | `GET /v1/colocalization` | Which traits share a causal signal at a given variant -- the standard mechanism question. |

## Conventions

- Every recipe sources `.env` and reads `FINNGENIE_TOKEN`. None of them embed credentials.
- Every recipe defaults to `?format=tsv` because the API returns TSV by default and TSV pipes into `duckdb` directly. JSON-shaped recipes pass `?format=json` explicitly and use `jq`.
- Every recipe sets sensible parameters at the top of the script (gene name, phenotype code, variant ID) so editing it for your own question is one line.
- Every recipe prints the URL it's about to hit. This makes it self-documenting: copy the URL into a browser with the bearer header and you get the same response.

## Adding a recipe

1. Pick a tool you used in the FinnGenie MCP or chat that you'd want to repeat.
2. Find the matching `/v1/*` endpoint in the [FinnGenie API docs](https://finngenie.fi/docs).
3. Write the curl version, parameterized at the top.
4. Append a row to the table above.

The recipe is the artifact. There is no Python wrapper to update.
