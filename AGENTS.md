# AGENTS.md - fgx

Project-specific guidance for agents working in this repository. This is the
public, runnable catalog of marimo notebooks for FinnGenie human-genetics
analysis. Planning and cross-instance coordination live in the primary
[`jx`](https://github.com/broadinstitute/jx) repo.

`README.md` is the human entry point. The skills under `.claude/skills/` are
the operational entry points: `getting-started` for first-run setup and
`compose-notebook` for adding or composing analyses.

## Validation Rule

After composing or editing any notebook in `notebooks/`, launch it in a
marimo sandbox kernel and run all cells before reporting the task complete.
Static checks do not catch wrong outputs, empty tables, stale endpoint
assumptions, auth mistakes, or plots that render but encode the wrong thing.

Minimal launch:

```bash
PORT=$(python -c "import socket; s=socket.socket(); s.bind(('127.0.0.1',0)); print(s.getsockname()[1])")
env -u PYTHONPATH uvx marimo edit --sandbox --headless --no-token --port $PORT notebooks/nbNN_*.py
```

Then run static checks:

```bash
uvx ruff check notebooks/
uvx ruff format notebooks/
uvx marimo check notebooks/*.py
```

**Then, last, refresh the molab session snapshot** for any notebook whose source
changed in this task:

```bash
env -u PYTHONPATH uvx marimo export session --sandbox notebooks/nbNN_*.py
```

Order matters. Session snapshots store a `code_hash` per cell, and molab
attaches the stored output only when the snapshot hash matches the source
cell. Any later edit to the notebook source - including a `ruff format`
whitespace pass - shifts every `code_hash` and silently strips outputs in
the public molab preview. Always regenerate snapshots **after** the final
formatter / source edit, and commit the regenerated `.json` files in the
same change that touched the `.py` files.

## Architecture

- Catalog over library. Helpers live as `@app.function` cells in numbered
  notebooks. Later notebooks import from earlier notebooks by adding
  `notebooks/` to `sys.path`.
- FinnGenie access is direct REST via `httpx`; there is no SDK, no schema
  cache, and no MCP dependency in the notebook path.
- `FINNGENIE_TOKEN` lives only in local `.env`. Never commit or paste it.
- API reads are live by design; do not add a committed cache unless the data
  surface changes.
- Do not add a Python package until repeated cross-notebook imports make the
  notebook-as-library pattern painful.

## When the Question Fits the Catalog

Almost every FinnGenie question should compose existing helpers:

- gene -> credible sets / colocalization -> `nb01_pcsk9_walkthrough`
- variant / rsID PheWAS -> `nb02_variant_phewas`
- phenotype locus zoom -> `nb03_phenotype_locus_zoom`
- gene exome results -> `nb04_gene_exome_burden`
- recessive Mendelian companion story -> `nb05_pign_cdg`
- variant pQTL direction of effect -> `nb06_variant_pqtl_function`
- available datasets/resources -> `nb07_data_catalog`

Read `.claude/skills/compose-notebook/SKILL.md` before writing new analysis
code.
