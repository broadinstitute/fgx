# fgx — FinnGen eXplore

> [!NOTE]
> Part of the [jx](https://github.com/broadinstitute/jx) experiment — not an official FinnGen or GeneGenie project.

An experiment in agent-driven scientific data exploration, built around [FinnGen](https://www.finngen.fi/) and partner human-genetics datasets exposed through the [GeneGenie](https://genegenie.broadinstitute.org/) results API — FinnGen R13 + UK Biobank + MVP meta-analyses, eQTL Catalogue R7, Open Targets 25.12, GTEx, Genebass, GenCC, Monarch (29 datasets at last count).

fgx is a curated catalog of [marimo](https://marimo.io) notebooks for human-genetics analysis, plus a thin skill that lets an agent compose new analyses from them.
Each notebook is both a runnable demonstration and a source of pure functions other notebooks can [import and reuse](https://docs.marimo.io/guides/reusing_functions/) directly.
Given a new human-genetics question, the agent picks relevant notebooks, composes their functions into a new notebook, executes it in a live kernel, and hands back a self-contained, re-runnable result.

GeneGenie ships [`https://genegenie.broadinstitute.org/api/v1/*`](https://genegenie.broadinstitute.org/) with bearer auth, predictable paths, and TSV by default — that collapses the data-access layer to `httpx.get`, so there is no Python SDK, no MCP server, no schema cache.

## The catalog

Each notebook ships with a committed session snapshot under [`notebooks/__marimo__/session/`](notebooks/__marimo__/session/) so the molab preview renders cell outputs without re-executing.

| Notebook | Role | Preview |
|---|---|---|
| [`nb01_pcsk9_walkthrough.py`](notebooks/nb01_pcsk9_walkthrough.py) | Gene -> credible sets -> lead variant -> colocalization | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb01_pcsk9_walkthrough.py) |
| [`nb02_variant_phewas.py`](notebooks/nb02_variant_phewas.py) | Variant PheWAS across endpoints | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb02_variant_phewas.py) |
| [`nb03_phenotype_locus_zoom.py`](notebooks/nb03_phenotype_locus_zoom.py) | Phenotype-driven locus zoom | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb03_phenotype_locus_zoom.py) |
| [`nb04_gene_exome_burden.py`](notebooks/nb04_gene_exome_burden.py) | Gene-based exome burden results | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb04_gene_exome_burden.py) |
| [`nb05_pign_cdg.py`](notebooks/nb05_pign_cdg.py) | Gene -> exome + curated `gene_disease` -> recessive-Mendelian companion to nb04 (PIGN-CDG / MCAHS1); imports `prepare_deleterious` from nb04 | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb05_pign_cdg.py) |
| [`nb06_variant_pqtl_function.py`](notebooks/nb06_variant_pqtl_function.py) | Variant -> pQTL credible sets -> direction-of-effect across proteins (ADAM17 / IBD demo replay) | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb06_variant_pqtl_function.py) |
| [`nb07_data_catalog.py`](notebooks/nb07_data_catalog.py) | Catalog introspection (`/datasets`, `/resources`, `/resource_metadata`) for "what's available?" | [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/broadinstitute/fgx/blob/main/notebooks/nb07_data_catalog.py) |

The agent-facing catalog table in `.claude/skills/compose-notebook/SKILL.md` is the detailed contract: it lists reusable helpers, import patterns, and current gotchas.

Related public catalogs of the same pattern: [jx](https://github.com/broadinstitute/jx) for JUMP Cell Painting, [prx](https://github.com/broadinstitute/prx) for PROSPECT chemical genetics, and [dmx](https://github.com/broadinstitute/dmx) for DepMap Breadbox.

## Getting started

Create an API key at [genegenie.broadinstitute.org](https://genegenie.broadinstitute.org/) (`MCP/API KEYS` -> `Create key`), then `cp .env.example .env` and paste your key in.

Clone this repo, open Claude Code inside it, and ask: *help me get started*.
The `getting-started` skill installs prereqs ([uv](https://docs.astral.sh/uv/) and the [marimo-pair](https://github.com/marimo-team/marimo-pair) skill), launches `nb01_pcsk9_walkthrough` in a live marimo kernel, and hands off to the `compose-notebook` skill for the actual analysis.

If you prefer to run setup by hand:

```bash
uv --version  # or: curl -LsSf https://astral.sh/uv/install.sh | sh
AGENT=claude-code  # or: codex
npx skills add marimo-team/marimo-pair -g --agent "$AGENT" -y
uvx marimo edit --sandbox notebooks/nb01_pcsk9_walkthrough.py
```

The skills reference in-repo notebooks and assets, so they only work in the cloned repo — there's no `npx skills add broadinstitute/fgx` flow.

## License

BSD 3-Clause — see [LICENSE](LICENSE).
