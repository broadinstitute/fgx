# fgx — FinnGen eXplore

An experiment in agent-driven scientific data exploration, built around [FinnGen](https://www.finngen.fi/) and partner human-genetics datasets exposed through the [FinnGenie](https://finngenie.broadinstitute.org/) results API (FinnGen R13 + UK Biobank + MVP meta-analyses, eQTL Catalogue R7, Open Targets 25.12, GTEx, Genebass, GenCC, Monarch -- 29 datasets at last count).

Same shape as [jx](https://github.com/broadinstitute/jx) -- a catalog of marimo notebooks plus thin operational skills -- with one difference: jx queries local DuckDB files, fgx hits a REST API.
FinnGenie ships [`https://finngenie.fi/api/v1/*`](https://finngenie.broadinstitute.org/) with bearer auth, predictable paths, and TSV by default.
That collapses the data-access layer to `httpx.get`, so there is no Python SDK, no MCP server, no schema cache -- nb01's first cell shows the bare-`curl` equivalent of every API call below it.
(The official [`genetics-mcp-server`](https://github.com/fulltiltgenomics/genetics-mcp-server) registers 29 tools that are all thin `httpx` wrappers around the same REST endpoints; none use MCP-specific primitives like sampling or elicitation, and a direct `curl` reproduces every result.)

A Python SDK would be lighter in syntax (`finngenie.credible_sets_by_gene("PCSK9")` vs three lines of `httpx`), but it adds a packaging layer that has to keep up with the API.
The TSV-by-default response means [`duckdb`](https://duckdb.org)'s `read_csv_auto` over HTTPS turns "API access" into "SQL access" without an SDK -- we accept the extra `httpx.get` lines per notebook in exchange for not maintaining a dependency.

## The catalog

Each notebook ships with a committed session snapshot under [`notebooks/__marimo__/session/`](notebooks/__marimo__/session/) so the molab preview renders cell outputs (plots, tables) without re-executing.
Click a badge to view the rendered notebook in [molab](https://docs.marimo.io/guides/molab/), or fork it from there.

| Notebook | Role | Preview |
|---|---|---|
| [`nb01_pcsk9_walkthrough.py`](notebooks/nb01_pcsk9_walkthrough.py) | Gene -> credible sets -> lead variant -> colocalization | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb01_pcsk9_walkthrough.py) |
| [`nb02_variant_phewas.py`](notebooks/nb02_variant_phewas.py) | Variant PheWAS across endpoints | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb02_variant_phewas.py) |
| [`nb03_phenotype_locus_zoom.py`](notebooks/nb03_phenotype_locus_zoom.py) | Phenotype-driven locus zoom | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb03_phenotype_locus_zoom.py) |
| [`nb04_gene_exome_burden.py`](notebooks/nb04_gene_exome_burden.py) | Gene-based exome burden results | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb04_gene_exome_burden.py) |
| [`nb05_pign_cdg.py`](notebooks/nb05_pign_cdg.py) | Gene -> exome + curated `gene_disease` -> recessive-Mendelian companion to nb04 (PIGN-CDG / MCAHS1); imports `prepare_deleterious` from nb04 | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb05_pign_cdg.py) |
| [`nb06_variant_pqtl_function.py`](notebooks/nb06_variant_pqtl_function.py) | Variant -> pQTL credible sets -> direction-of-effect across proteins (ADAM17 / IBD demo replay) | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb06_variant_pqtl_function.py) |
| [`nb07_data_catalog.py`](notebooks/nb07_data_catalog.py) | Catalog introspection (`/datasets`, `/resources`, `/resource_metadata`) for "what's available?" | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb07_data_catalog.py) |

Related public catalogs of the same pattern: [jx](https://github.com/broadinstitute/jx) for JUMP Cell Painting, [prx](https://github.com/broadinstitute/prx) for PROSPECT chemical genetics, and [dmx](https://github.com/broadinstitute/dmx) for DepMap Breadbox.

## Getting started

1. Create an API key at [finngenie.broadinstitute.org](https://finngenie.broadinstitute.org/) (`MCP/API KEYS` -> `Create key`).
   The same key works for both the MCP and the REST API.
2. `cp .env.example .env` and paste your key.
3. From a fresh Claude Code session in this clone, ask: *help me get started*.
   The [`getting-started`](.claude/skills/getting-started/SKILL.md) skill installs `uv` and `marimo-pair`, launches nb01 in a live marimo kernel, and hands off to interactive composition.

If you want to skip the agent-assist and just open the notebook: `just notebook` (uses `uvx marimo edit --sandbox`, provisions deps from the PEP 723 header on first launch).

Prerequisites: [`uv`](https://docs.astral.sh/uv/) is the only thing you need to install; the notebook's PEP 723 header pulls in marimo, polars, httpx, altair, and python-dotenv.
Optional: [`duckdb`](https://duckdb.org/docs/installation/) on `$PATH` if you want to run ad-hoc SQL against the same TSV endpoints from your shell.

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
    nb02_variant_phewas.py           # variant PheWAS across endpoints
    nb03_phenotype_locus_zoom.py     # phenotype-driven locus zoom
    nb04_gene_exome_burden.py        # gene-based exome burden results
    nb05_pign_cdg.py                 # recessive-Mendelian companion to nb04 (PIGN/MCAHS1)
    nb06_variant_pqtl_function.py    # variant -> pQTL credible sets -> direction-of-effect (ADAM17/IBD)
    nb07_data_catalog.py             # /datasets + /resources introspection ("what's available?")
    __marimo__/session/              # committed session snapshots for molab previews
```

That's the whole repo.
`getting-started` walks a fresh clone to a running marimo kernel; `compose-notebook` covers picking the right path (edit-in-place vs copy-fork) and the right endpoint when the user asks for a new genetics analysis.
Both skills are intentionally thin -- as the catalog grows and cross-file function reuse appears, `compose-notebook` will need to grow into a per-module catalog table the way [jx's compose-notebook](https://github.com/broadinstitute/jx/blob/main/.claude/skills/compose-notebook/SKILL.md) does.

## Citation

If fgx is useful in your work, please cite both fgx and the underlying FinnGenie service.
CITATION.cff TBD.

## License

BSD 3-Clause — see [LICENSE](LICENSE).
