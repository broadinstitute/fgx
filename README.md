# fgx - FinnGen eXplore

An experiment in agent-driven scientific data exploration, built around [FinnGen](https://www.finngen.fi/) and partner human-genetics datasets exposed through the [FinnGenie](https://finngenie.broadinstitute.org/) results API (FinnGen R13 + UK Biobank + MVP meta-analyses, eQTL Catalogue R7, Open Targets 25.12, GTEx, Genebass, GenCC, Monarch -- 29 datasets at last count).

## The hypothesis

fgx is a catalog of marimo notebooks against the FinnGenie REST API, and nothing else. No Python SDK, no MCP server, no schema cache. The substrate is the [`marimo`](https://marimo.io) notebook itself; the API is the contract; `httpx.get` is the entire data-access layer. The reusable functions inside each notebook are the only thing that gets shared.

The argument behind that minimalism is the FinnGenie REST API itself: bearer auth, predictable paths, TSV by default (which streams straight into [`duckdb`](https://duckdb.org)). Given that, every wrapper above the API -- MCP server, Python SDK, custom CLI -- is overhead. The shorter path is to write the notebook directly against `https://finngenie.fi/api/v1/*`. Each notebook becomes simultaneously a tutorial, a regression test, and a reusable function library for the next analysis.

## Why notebooks instead of MCP, skills, or a Python SDK?

We measured. The FinnGenie MCP server (`fulltiltgenomics/genetics-mcp-server`) registers 29 tools; every one is a thin `httpx` wrapper around `https://finngenie.fi/api/v1/*` or an internal BigQuery proxy. None use the two MCP-specific primitives (sampling, elicitation). A direct `curl` reproduces every result. See [the evaluator pass](../jx-dev/reference/2026-05-01__voa__finngenie-mcp-evaluator.md) for the full per-tool verdict. (The DepMap MCP scored the same way in April.)

A skill file would be the next-lightest layer above the API. fgx goes one step further: the marimo notebook is the skill -- its first cell shows the bare-`curl` equivalent of every API call below, so an agent reading the notebook learns the API surface from the same artifact it composes against. There's no separate `SKILL.md` to drift out of sync with the notebook code.

A Python SDK would be lighter still in syntax (`finngenie.credible_sets_by_gene("PCSK9")` vs three lines of `httpx`), but it adds a packaging layer that has to keep up with the API. The TSV-by-default response means `duckdb`'s `read_csv_auto` over HTTPS turns "API access" into "SQL access" without an SDK. We accept the few extra `httpx.get` lines per notebook in exchange for not maintaining a dependency.

## Getting started

1. Create an API key at [finngenie.broadinstitute.org](https://finngenie.broadinstitute.org/) (`MCP/API KEYS` -> `Create key`). The same key works for both the MCP and the REST API.
2. `cp .env.example .env` and paste your key.
3. From a fresh Claude Code session in this clone, ask: *help me get started*. The [`getting-started`](.claude/skills/getting-started/SKILL.md) skill installs `uv` and `marimo-pair`, launches nb01 in a live marimo kernel, and hands off to interactive composition.

If you want to skip the agent-assist and just open the notebook: `just notebook` (uses `uvx marimo edit --sandbox`, provisions deps from the PEP 723 header on first launch).

Prerequisites: [`uv`](https://docs.astral.sh/uv/) is the only thing you need to install; the notebook's PEP 723 header pulls in marimo, polars, httpx, altair, and python-dotenv. Optional: [`duckdb`](https://duckdb.org/docs/installation/) on `$PATH` if you want to run ad-hoc SQL against the same TSV endpoints from your shell.

## Layout

```
fgx/
  README.md                          # this file
  Justfile                           # `just notebook`
  .env.example                       # FINNGENIE_TOKEN= placeholder
  .claude/skills/
    getting-started/SKILL.md         # cold-clone bootstrap
  notebooks/
    nb01_pcsk9_walkthrough.py        # gene -> credible sets -> lead variant -> colocalization
```

That's the whole repo. New analyses become `notebooks/nbNN_*.py`. When the catalog grows past the point where the README's index can carry it, fgx will likely add a [`compose-notebook`](https://github.com/broadinstitute/jx/blob/main/.claude/skills/compose-notebook/SKILL.md) skill (the same pattern jx uses); that's a scale-driven addition, not a starting point.

## License

BSD 3-Clause. See [LICENSE](LICENSE).

## Citation

If fgx is useful in your work, please cite both fgx and the underlying FinnGenie service. CITATION.cff TBD.
