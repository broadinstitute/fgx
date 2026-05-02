# fgx - FinnGen eXplore

An experiment in agent-driven scientific data exploration, built around [FinnGen](https://www.finngen.fi/) and partner human-genetics datasets exposed through the [FinnGenie](https://finngenie.broadinstitute.org/) results API (FinnGen R13 + UK Biobank + MVP meta-analyses, eQTL Catalogue R7, Open Targets 25.12, GTEx, Genebass, GenCC, Monarch -- 29 datasets at last count).

## The hypothesis

Same shape as [jx](https://github.com/broadinstitute/jx) -- a catalog of marimo notebooks plus thin operational skills -- with one difference: jx queries local DuckDB files, fgx hits a REST API. FinnGenie ships [`https://finngenie.fi/api/v1/*`](https://finngenie.broadinstitute.org/) with bearer auth, predictable paths, and TSV by default. That collapses the data-access layer to `httpx.get`, so there is no Python SDK, no MCP server, no schema cache -- nb01's first cell shows the bare-`curl` equivalent of every API call below it. (The official [`genetics-mcp-server`](https://github.com/fulltiltgenomics/genetics-mcp-server) registers 29 tools that are all thin `httpx` wrappers around the same REST endpoints; none use MCP-specific primitives like sampling or elicitation, and a direct `curl` reproduces every result.)

A Python SDK would be lighter in syntax (`finngenie.credible_sets_by_gene("PCSK9")` vs three lines of `httpx`), but it adds a packaging layer that has to keep up with the API. The TSV-by-default response means [`duckdb`](https://duckdb.org)'s `read_csv_auto` over HTTPS turns "API access" into "SQL access" without an SDK -- we accept the extra `httpx.get` lines per notebook in exchange for not maintaining a dependency.

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
    compose-notebook/SKILL.md        # edit nb01 in place vs copy-fork to nb02_*
  notebooks/
    nb01_pcsk9_walkthrough.py        # gene -> credible sets -> lead variant -> colocalization
```

That's the whole repo. `getting-started` walks a fresh clone to a running marimo kernel; `compose-notebook` covers picking the right path (edit-in-place vs copy-fork) and the right endpoint when the user asks for a new genetics analysis. Both skills are intentionally thin -- when the catalog grows past one notebook with cross-file function reuse, `compose-notebook` will need to grow into a per-module catalog table the way [jx's compose-notebook](https://github.com/broadinstitute/jx/blob/main/.claude/skills/compose-notebook/SKILL.md) does.

## License

BSD 3-Clause. See [LICENSE](LICENSE).

## Citation

If fgx is useful in your work, please cite both fgx and the underlying FinnGenie service. CITATION.cff TBD.
