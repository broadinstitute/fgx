# AGENTS.md - fgx

Project-specific guidance for agents working in this repository.
This is the public, runnable catalog of marimo notebooks for GeneGenie human-genetics analysis.
Planning and cross-instance coordination live in the primary [`jx`](https://github.com/broadinstitute/jx) repo.

`README.md` is the human entry point.
This catalog uses the shared [vignette-catalog-skills](https://github.com/carpenter-singh-lab/vignette-catalog-skills) (`vignette-catalog-setup` for first-run setup, `vignette-catalog-compose-notebook` for adding or composing analyses); its specifics live in `catalog.toml`.
The skills are installed via `npx skills add carpenter-singh-lab/vignette-catalog-skills --agent claude-code -y`, recorded in the tracked `skills-lock.json`, but **not vendored** - `.claude/skills/*` is gitignored, so restore them on a fresh clone before use.

## Launching notebooks

Always use `--sandbox` so the PEP 723 inline metadata is provisioned:

```bash
uvx marimo edit --sandbox notebooks/nbNN_*.py
```

Do not improvise alternative launch commands.
`--sandbox` is what makes `uvx marimo` read each notebook's `/// script` dependency block; without it every notebook fails with `ModuleNotFoundError`.

## Validation Rule

After composing or editing any notebook in `notebooks/`, launch it in a marimo sandbox kernel and run all cells before reporting the task complete.
Static checks do not catch wrong outputs, empty tables, stale endpoint assumptions, auth mistakes, or plots that render but encode the wrong thing.

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

**Then, last, refresh the molab session snapshot** for any notebook whose source changed in this task:

```bash
env -u PYTHONPATH uvx marimo export session --sandbox notebooks/nbNN_*.py
```

Order matters.
Session snapshots store a `code_hash` per cell, and molab attaches the stored output only when the snapshot hash matches the source cell.
Any later edit to the notebook source - including a `ruff format` whitespace pass - shifts every `code_hash` and silently strips outputs in the public molab preview.
Always regenerate snapshots **after** the final formatter / source edit, and commit the regenerated `.json` files in the same change that touched the `.py` files.

## Architecture

- Catalog over library.
  Helpers live as `@app.function` cells in numbered notebooks.
  Later notebooks import from earlier notebooks by adding `notebooks/` to `sys.path`.
- GeneGenie access is direct REST via `httpx`; there is no SDK, no schema cache, and no MCP dependency in the notebook path.
- `GENEGENIE_TOKEN` lives only in local `.env`.
  Never commit or paste it.
- API reads are live by design; do not add a committed cache unless the data surface changes.
- Do not add a Python package until repeated cross-notebook imports make the notebook-as-library pattern painful.

## Conventions

- Prose in `.md` files uses semantic line breaks: one sentence per line, no hard wrapping at a column count.
  Markdown collapses single newlines inside a paragraph, so the rendered output is unchanged, but diffs stay local to the edited sentence instead of re-flowing every line below it.
  Applies to `AGENTS.md`, `.claude/skills/**/SKILL.md`, and any other prose-heavy markdown we revise often.
  `ruff`'s `line-length = 120` setting is for Python only; there is no column rule on markdown.

## When the Question Fits the Catalog

Almost every GeneGenie question should compose existing helpers:

- gene -> credible sets / colocalization -> `nb01_pcsk9_walkthrough`
- variant / rsID PheWAS -> `nb02_variant_phewas`
- phenotype locus zoom -> `nb03_phenotype_locus_zoom`
- gene exome results -> `nb04_gene_exome_burden`
- recessive Mendelian companion story -> `nb05_pign_cdg`
- variant pQTL direction of effect -> `nb06_variant_pqtl_function`
- available datasets/resources -> `nb07_data_catalog`

Read the installed `vignette-catalog-compose-notebook` skill (and `catalog.toml`'s `[[vignette]]` table) before writing new analysis code.

## GeneGenie API gotchas

These are fgx-specific endpoint-semantics papercuts (they do not live in the shared skill; the generic marimo/molab gotchas do).

- **TSV is the API default** (`Content-Type: text/tab-separated-values`). `fetch_json` *unconditionally* injects `format=json`, so pick the helper by the response shape you want; you cannot flip format via a kwarg.
- **Filter before printing.** `credible_sets_by_gene/PCSK9` is ~3,200 rows / ~1.3 MB and `resource_metadata/finngen` ~3,300 rows. Summarize in-cell (`.filter`, `.head(20)`, `group_by+agg`) and never let a raw full-size response land in the transcript.
- **Do not infer semantics from a path name.** `exome_results_by_gene` returns per-variant single-variant rows, not gene-level burden. `mlog10p` saturates at floating-point limits (`mlog10p == 324` with `se == 0`), `pip` can be null, and `most_severe` is a VEP prediction, not a clinical classification. Read the actual response (or the OpenAPI spec) first.
- **`/search` uses `q`, not `query`.** It returns JSON; call `client().get(f"{BASE}/search", params={"q": term})` directly and parse `r.json()` - `fetch_json` (which forces `format=json`) can 422 here.
- **`/credible_sets_by_phenotype` needs `{resource}/{phenotype}`** in the path (e.g. `finngen/T2D`, `open_targets/GCST...`). Omitting the resource is a 404.
- **Identifiers are system-specific.** Phenotype codes (`I9_CHD`, `T2D`) are FinnGen definitions - look them up via `fetch_json("/trait_name_mapping")`. An rsID can resolve to a different alt allele than the credible-set store indexed; when a known-good identifier returns 0 rows, suspect the mapping, not the data.
