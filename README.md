# fgx - FinnGen eXplore

An experiment in agent-driven scientific data exploration, built around [FinnGen](https://www.finngen.fi/) and partner human-genetics datasets exposed through the [FinnGenie](https://finngenie.broadinstitute.org/) results API (FinnGen R13 + UK Biobank + MVP meta-analyses, eQTL Catalogue R7, Open Targets 25.12, GTEx, Genebass, GenCC, Monarch -- 29 datasets in total at last count).

## The hypothesis

fgx is a catalog of **REST recipes**, not a tool layer. The FinnGenie team operates a clean REST API at `https://finngenie.fi/api/v1/*`: bearer auth, predictable paths, TSV by default (which streams straight into `duckdb` / `polars`). Given that, every wrapper above the API -- MCP server, Python SDK, skill file -- is overhead. The shorter path is to publish worked queries that humans copy and agents read.

Three properties make this work:

- **The recipes are the docs.** Each `recipes/q*.sh` is a self-contained shell script that sources a token, hits one endpoint, and pipes the result through `jq` or `duckdb` to answer one specific question. The script is simultaneously a tutorial, a regression test, and a thing you can pipe into a notebook.
- **TSV-by-default beats JSON.** The API's default `Content-Type: text/tab-separated-values` means a one-line `curl ... | duckdb -c "..."` does what a Python SDK would do in twenty. The agent reads the recipe, learns the pattern, and composes new queries against the same surface.
- **No local data layer.** Unlike sibling pilots that build a local DuckDB from a snapshot, fgx is internet-first by design. The FinnGenie REST API is the source of truth; recipes hit it live, with the bearer token as the only piece of state.

A thin layer of structure (recipes named `q01_*`, `q02_*`, ...; one paragraph of context per recipe; a shared `.env` for the token) is the entire library. There is no Python package, no MCP server, no schema cache.

If this works for FinnGenie, the pattern transfers: any dataset with a clean REST surface = new recipe catalog, same machinery.

## Why not skills, an MCP, or a Python SDK?

We measured. The FinnGenie MCP server (`fulltiltgenomics/genetics-mcp-server`) registers 29 tools; every one is a thin `httpx` wrapper around `https://finngenie.fi/api/v1/*` or an internal BigQuery proxy. None use the two MCP-specific primitives (sampling, elicitation). A direct `curl` reproduces every result. See [the evaluator pass](../jx-dev/reference/2026-05-01__voa__finngenie-mcp-evaluator.md) for the full per-tool verdict. (The DepMap MCP scored the same way in April; two-for-two on Broad-built MCP servers wrapping internal REST APIs.)

A skill file would be the next-lightest layer above the API. fgx goes one step further: the recipes themselves are the skill. An agent reads `recipes/README.md`, sees the available queries, picks the closest one, edits parameters. No separate `SKILL.md` to keep in sync.

## Getting started

1. Create an API key at [finngenie.broadinstitute.org](https://finngenie.broadinstitute.org/) (`MCP/API KEYS` -> `Create key`). The same key works for both the MCP and the REST API.
2. `cp .env.example .env` and paste your key.
3. `bash recipes/q01_list_datasets.sh` to confirm the token works.
4. Browse `recipes/` and edit any script in place to ask your own variant of the question.

Prerequisites: `curl`, `jq`, and (for the recipes that aggregate) [`duckdb`](https://duckdb.org/docs/installation/). All three are one `brew install` away.

## Layout

Three layers of the same surface, climbing in capability and dependencies:

```
fgx/
  README.md           # this file
  Justfile            # `just fetch | build-db | render | notebook`
  .env.example        # FINNGENIE_TOKEN= placeholder
  recipes/            # one-shot bash + curl + jq/duckdb
    q01_list_datasets.sh
    q02_credible_sets_by_gene.sh         # PCSK9 walkthrough
    q03_credible_sets_by_phenotype.sh    # I9_CHD walkthrough
    q04_colocalization.sh                # variant -> shared causal signals
  queries/            # ggsql charts against a local DuckDB cached from API
    q01_top_traits_near_pcsk9.gsql
    q02_top_chd_loci.gsql
  notebooks/          # marimo, composes endpoints with reactive plots
    nb01_pcsk9_walkthrough.py            # composes recipes q02 + q04
```

Each layer has its own audience:

- **`recipes/`** -- the smallest demonstration that the API needs no tool layer. Bash + curl + jq/duckdb. Read-once-and-go for humans, copy-and-edit for agents.
- **`queries/`** -- declarative charts. ggsql files run against a local DuckDB cached by `just fetch`. The auth boundary lives entirely in `fetch`; the .gsql files are auth-free SQL.
- **`notebooks/`** -- composition + plots. Marimo notebooks chain multiple endpoints (e.g. gene -> credible sets -> lead variant -> colocalization). Reactive kernel feedback while exploring; reusable functions across notebooks.

## License

BSD 3-Clause. See [LICENSE](LICENSE).

## Citation

If fgx is useful in your work, please cite both fgx and the underlying FinnGenie service. See `CITATION.cff` (TODO).
