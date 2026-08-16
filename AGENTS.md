# AGENTS.md - fgx

Project-specific guidance for agents working in this repository.
This is the public, runnable catalog of marimo notebooks for GeneGenie human-genetics analysis.
Planning and cross-instance coordination live in the primary [`jx`](https://github.com/broadinstitute/jx) repo.

`README.md` is the human entry point.
This catalog uses the shared [vignette-catalog-skills](https://github.com/carpenter-singh-lab/vignette-catalog-skills), with `vignette-catalog-compose-notebook` handling setup, execution, and composition; its specifics live in `catalog.toml`.
The skills are recorded in the tracked `skills-lock.json` but not vendored; restore them with the exact commands in `README.md` after cloning.

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

Then run static checks, with ruff pinned:

```bash
uvx ruff@0.16.2 check notebooks/
uvx ruff@0.16.2 format notebooks/
uvx marimo check notebooks/*.py
```

The pin is load-bearing, not fussiness.
`ruff format` output shifts between ruff releases, and any reformat changes every `code_hash` below, which silently strips outputs from the molab preview.
An unpinned `uvx ruff` therefore turns a routine lint into a snapshot-stripping edit on somebody else's machine.
The lint rule set is pinned separately by an explicit `select` in `pyproject.toml`, because ruff's implicit defaults have widened across releases.
When you do bump the pin, expect to reformat and regenerate every snapshot in the same change.

**Then, last, refresh the molab session snapshot** for any notebook whose source changed in this task:

```bash
env -u PYTHONPATH uvx marimo export session --sandbox notebooks/nbNN_*.py
```

Order matters.
Session snapshots store a `code_hash` per cell, and molab attaches the stored output only when the snapshot hash matches the source cell.
Any later edit to the notebook source - including a `ruff format` whitespace pass - shifts every `code_hash` and silently strips outputs in the public molab preview.
The PEP 723 `/// script` header counts too: it feeds the `setup` cell's hash, so bumping a dependency pin is a snapshot-invalidating edit even though no cell body changed.
Always regenerate snapshots **after** the final formatter / source edit, and commit the regenerated `.json` files in the same change that touched the `.py` files.

This failure is silent by construction - a stale snapshot looks fine in the diff and simply renders blank on molab - so verify rather than trust.
This check reads hashes offline and needs no kernel or API token:

```bash
env -u PYTHONPATH uv run --no-project --with 'marimo<0.23.4' python3 - <<'EOF'
import json
from pathlib import Path
from marimo._ast.load import load_app
from marimo._utils.code import hash_code

for py in sorted(Path("notebooks").glob("nb*.py")):
    snap = Path("notebooks/__marimo__/session") / (py.name + ".json")
    stored = {c["id"]: c["code_hash"] for c in json.loads(snap.read_text())["cells"]}
    live = {cid: hash_code(cd.code) for cid, cd in load_app(str(py))._cell_manager._cell_data.items()}
    stale = [k for k in stored if stored[k] != live.get(k)]
    print(f"{py.name}: {'STALE ' + ','.join(stale) if stale else 'ok'}")
EOF
```

## Architecture

- Catalog over library.
  Helpers live as `@app.function` cells in numbered notebooks.
  Later notebooks import from earlier notebooks by adding `notebooks/` to `sys.path`.
- GeneGenie access is direct REST via `httpx`; there is no SDK, no schema cache, and no MCP dependency in the notebook path.
- `GENEGENIE_TOKEN` lives only in local `.env`.
  Never commit or paste it.
- API reads are live by design; do not add a committed cache unless the data surface changes.
- Every PEP 723 dependency carries an upper bound.
  The notebooks are the deployable unit and there is no lockfile, so an unbounded requirement means a future release picks itself; that is how a polars deprecation and an altair/Python-3.14 incompatibility both landed as runtime breaks.
  Bound the next major (`polars<2`), or the next minor for pre-1.0 packages whose minors break (`httpx<0.29`).
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

- **TSV is the API default** (`Content-Type: text/tab-separated-values`).
  `fetch_json` *unconditionally* injects `format=json`, so pick the helper by the response shape you want; you cannot flip format via a kwarg.
- **Filter before printing.** `credible_sets_by_gene/PCSK9` is ~3,200 rows / ~1.3 MB and `resource_metadata/finngen` ~3,300 rows.
  Summarize in-cell (`.filter`, `.head(20)`, `group_by+agg`) and never let a raw full-size response land in the transcript.
- **Do not infer semantics from a path name.** `exome_results_by_gene` returns per-variant single-variant rows, not gene-level burden.
  `mlog10p` saturates at floating-point limits (`mlog10p == 324` with `se == 0`), `pip` can be null, and `most_severe` is a VEP prediction, not a clinical classification.
  Read the actual response (or the OpenAPI spec) first.
- **`/search` uses `q`, not `query`.** Passing `query=` 422s with `Unknown query parameter(s): query`.
  `/search` does accept `format`, so plain `fetch_json("/search", q=term)` works - the 422 you will hit here is the parameter name, not the format.
- **`/credible_sets_by_phenotype` needs `{resource}/{phenotype}`** in the path (e.g.
  `finngen/T2D`, `open_targets/GCST...`).
  Omitting the resource is a 404.
- **Identifiers are system-specific.** Phenotype codes (`I9_CHD`, `T2D`) are FinnGen definitions - look them up via `fetch_json_raw("/trait_name_mapping")` (it is one of the no-`format` endpoints).
  An rsID can resolve to a different alt allele than the credible-set store indexed; when a known-good identifier returns 0 rows, suspect the mapping, not the data.
- **`/credible_sets_by_phenotype_leads` answers JSON, not TSV**, unlike most of the API.
  One row per credible set (the highest-PIP variant in the 95% set), which makes it the right entry point for counting *independent signals* rather than tag variants.
  It also returns sets regardless of significance - filter on `mlog10p` yourself.
- **The batch `POST /credible_sets_by_variant` separator is a NEWLINE.** The body is `{"variants": "<v1>\n<v2>\n..."}`; commas, spaces and semicolons all 422 with "variant needs to contain four fields".
  Chunk at ~40 variants.
- **Eleven endpoints reject `format` but answer JSON anyway.** The ones this catalog touches are `/datasets`, `/resources`, `/rsid/variants`, `/trait_name_mapping` and `/dataset_display_names`.
  They 422 with `Unknown query parameter(s): format`, so `fetch_json` cannot reach them (it injects `format=json` unconditionally) and `fetch_tsv` cannot parse the JSON that comes back.
  Neither helper works: use nb01's **`fetch_json_raw`**.
  Before writing against an endpoint, check whether `/openapi.json` gives it a `format` parameter; the 422 body also names everything it does accept.
- **`/search` caps `limit` at 100** and its `types` values are plural (`phenotypes`, `genes`).
  Use it for discovery, never for exhaustive enumeration - to enumerate, probe the endpoint and record the status codes.
- **`/resource_metadata` dtypes vary by resource.** Count columns come back Int64 for some resources and String for others, so `pl.concat` across resources raises `SchemaError`.
  Cast, or pass `how="vertical_relaxed"`.
- **The phenotype list is not the fine-mapping store.** A phenotype in `/resource_metadata` frequently 404s on the credible-set endpoints (in FinnGen R14, 371 of 602 cancer endpoints do).
  Never infer coverage from the metadata list; probe and keep the failures.
