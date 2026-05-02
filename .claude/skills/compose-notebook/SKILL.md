---
name: compose-notebook
description: Compose a new genetics analysis in the fgx repo by reusing the marimo notebook(s) under notebooks/ as the substrate, hitting FinnGenie's REST API directly via httpx.get. Trigger when the user asks for a notebook, analysis, figure, or vignette against human-genetics data exposed by FinnGenie -- "credible sets near gene X", "loci driving phenotype Y", "what colocalizes at variant Z", "pull the GWAS signal for...", "compare PIPs across resources", "PheWAS for rsid", "look up an rsID", "search for a phenotype", or any composition question against /api/v1/* -- credible_sets_by_{gene,phenotype,variant,qtl_gene,region,id} and credible_sets/{id}/stats; colocalization_by_{variant,credible_set_id}; exome_results_by_{gene,phenotype,region,variant}; summary_stats/{resource}/{data_type}; nearest_genes; peak_to_genes; gene_disease; rsid/variants; search; phenotype/{resource}/{phenocode}/markdown; resource_metadata; resources; datasets. Authoritative path list at /api/v1/openapi.json (browseable at /api/v1/docs); confirm against the spec rather than guessing if the user asks about an endpoint not listed here. Use this instead of writing standalone httpx code from scratch when the existing helpers already cover the call shape, and instead of asking the user to wire auth themselves.
---

# Compose a new marimo notebook in fgx

## What this skill is for

fgx is a single-substrate repo: marimo notebooks against FinnGenie's REST API at `https://finngenie.fi/api/v1/*` via `httpx.get`, and nothing else. There is no Python SDK, no MCP server, no schema cache. Every notebook carries its own copy of the same three `@app.function` helpers (`client`, `fetch_tsv`, `fetch_json`) -- the helpers are tiny, self-contained, and easier to fork than to import across files.

When the user asks a genetics question that the FinnGenie API can answer, your job is to (a) decide whether to edit the closest existing notebook in place or copy-fork to a new `nbNN_*.py`, and (b) compose against the existing helpers rather than re-implementing auth and TSV/JSON parsing from scratch.

## The catalog

Each notebook is a self-contained vignette covering one entry-point + visualization combination. Pick the one whose shape matches the user's question; if none does, copy-fork the closest.

| File | Entry point | Endpoint chain | Visualization |
|---|---|---|---|
| `notebooks/nb01_pcsk9_walkthrough.py` | gene | `credible_sets_by_gene` -> `colocalization_by_variant` | bar (top traits + colocalizing pairs) |
| `notebooks/nb02_variant_phewas.py` | rsID / variant | `rsid/variants` -> `credible_sets_by_variant` -> `nearest_genes` | PheWAS dot plot |
| `notebooks/nb03_phenotype_locus_zoom.py` | phenocode + resource | `credible_sets_by_phenotype` -> `nearest_genes` (per top locus) | Manhattan with gene labels |
| `notebooks/nb04_gene_exome_burden.py` | gene (rare-variant arm) | `exome_results_by_gene` -> `gene_disease` | forest plot of betas with CIs |

Reusable helpers in every notebook: `fetch_tsv(path, **params) -> pl.DataFrame` and `fetch_json(path, **params) -> list[dict] | dict`. `path` is the part after `/api/v1`, kwargs become query parameters. The underlying `client()` builder is private to the notebook -- don't call it from cells, don't rename it `_client` (see the underscore-mangling Gotcha).

## Portfolio principle

The library is *minimal-spanning by design*: each notebook earns its place by adding an axis -- a new entry point (gene / variant / phenotype / region) or a new visualization shape (bar / PheWAS / Manhattan / forest / locus-zoom / network) -- not just a new question. When the user asks for nb05+, ask "what axis does this add?" first. If the answer is "same shape as an existing notebook, just different parameters", that's Path A: edit the existing notebook. If it's a genuinely new entry point or viz, that's Path B: copy-fork. Resist the urge to spawn a near-duplicate -- it dilutes the catalog and gives the next composer two notebooks to choose between when one would do.

## Two paths: edit-in-place vs copy-fork

### Path A: edit the closest-matching notebook in place

Use this when the new question has the **same shape** as an existing notebook, just a different parameter:

- "do nb01 for APOE" -> change `GENE = "PCSK9"` to `GENE = "APOE"` in nb01, the rest cascades reactively.
- "do nb02 for `rs7412`" -> change `RSID` in nb02.
- "do nb03 for T2D on UKBB" -> change `PHENOTYPE` and `RESOURCE` in nb03.
- "do nb04 for LDLR" -> change `GENE` in nb04.

Marimo's reactivity reruns dependent cells; don't manually kick anything. This is the cheapest path and preserves the per-notebook canonical example.

### Path B: copy-fork to `notebooks/nbNN_<topic>.py`

Use this when the analysis shape **differs** from every existing notebook. Symptoms:

- Different endpoint chain (e.g. region-first instead of gene-first, or chaining `summary_stats` after `credible_sets`).
- Different output (locus-zoom around a single region, multi-resource comparison, network of trait-trait colocalization).
- A reusable helper would need to live in two notebooks, and copy-pasting it twice is uglier than letting one notebook import it from a sibling.

When forking: copy the closest-matching notebook -> `nbNN_<topic>.py`, edit the PEP 723 header if you need new deps, replace the analysis cells, keep the setup block (loads `.env` and re-defines the helpers verbatim). Keep the `marimo<0.23.4` and `jedi<0.20.0` pins in the header verbatim -- both work around `uvx marimo edit --sandbox` resolver bugs. If you trim them, the next launch will fail at venv provisioning, not at runtime.

## Process for a new composition

1. **Read the closest existing notebook.** The first cell explains what the notebook does; the helpers (`client`, `fetch_tsv`, `fetch_json`) are at the top as `@app.function`. Skim before composing.
2. **Identify the endpoint.** Map the user's question to a `/api/v1/*` path. The full surface is in the OpenAPI spec at <https://finngenie.fi/api/v1/openapi.json> (browseable Swagger at <https://finngenie.fi/api/v1/docs>) -- 28 operations across 13 tags. Read the spec before guessing; the description's path list is a reminder, not a contract. If you're unsure of the response shape, run `curl -fsS -H "Authorization: Bearer $FINNGENIE_TOKEN" "https://finngenie.fi/api/v1/<path>?format=json"` once to confirm, then write the cell.
3. **Pick the path (A or B).** When in doubt, copy-fork.
4. **Validate against the live kernel via marimo-pair.** This is non-negotiable, and the rest of this step exists because there are three failure modes that look like shortcuts but aren't.
    - **Browser session is a prerequisite, not optional.** `uvx marimo edit --headless` launches a kernel, but `execute-code.sh` won't attach until the user opens the URL in a browser tab and a websocket session is alive. If you call `execute-code.sh` and see `No active sessions on the server`, the right move is to tell the user "open `http://127.0.0.1:<port>/` in your browser, then ping me" -- not to fall back to running the same logic in a standalone `uv run python3` script. A parallel script imports nothing from the notebook and runs with its own module scope, so it can't catch `@app.function` vs `@app.cell` scoping bugs, marimo auto-formatter rewrites, or stale compiled cells. False green there is worse than no signal.
    - **Hot-reload reloads the file but does not re-run cells.** When you edit a `.py` file, marimo notices and re-imports module-level definitions, but the kernel still holds the old compiled cell objects in memory until something triggers a rerun. So after every code edit, ask the user to hit "Run all" (Cmd+Shift+R, or the play-all button in the toolbar) before you re-validate. Otherwise `execute-code.sh` will report errors against the *previous* version of the code, which is maddening to debug.
    - **One cell at a time.** Once the session is attached and cells have been re-run, drive validation through `execute-code.sh --port <port> -c "..."` cell by cell -- not by re-pasting whole notebook bodies. The setup block (`with app.setup:`) is already provisioned by the PEP 723 header, and `@app.function` helpers like `fetch_tsv` and `fetch_json` are visible in the scratchpad's namespace; a single cell returns a `polars.DataFrame` you can summarize inline.
5. **End with extension prompts.** Every vignette closes with a `## To extend` markdown cell listing 2-3 concrete next prompts (e.g. "swap PCSK9 for LDLR", "add a second resource and stack the dot plots"). This turns each notebook into a launchpad rather than a dead end.
6. **Plot when it adds signal.** Altair is in the PEP 723 deps; a one-cell chart is often more informative than a 20-row table head.
7. **Don't fabricate.** If a helper doesn't exist or an endpoint returns an empty result, say so and ask the user how to proceed -- don't invent a fallback path.

## Gotchas

The first three are tooling traps that cost real session time when missed; the middle three are API papercuts; the last three are biology / endpoint-semantics nuances.

- **Validate via marimo-pair against a live kernel, never via a standalone `uv run python3` script.** A parallel script can't see `@app.function` mangling, can't catch auto-formatter rewrites, and gives you a false green that costs more than no signal. If `execute-code.sh` reports `No active sessions`, the fix is to ask the user to open the browser URL -- not to bypass.
- **Polars-only; never call `.to_pandas()`.** Altair 6+ accepts polars DataFrames directly via narwhals (`alt.Chart(df_polars)` works -- no conversion). The PEP 723 deps don't include `pandas` or `pyarrow`, and modern polars routes `.to_pandas()` through Arrow, so any `.to_pandas()` call raises `ModuleNotFoundError: pa.Table requires 'pyarrow' module`. The fix is the cleaner one: pass polars DataFrames straight to Altair. Don't add pandas / pyarrow to the deps to "fix" it -- the dependency is the bug, not the missing package.
- **Marimo mangles any name starting with `_`.** Any name beginning with an underscore is treated as cell-local and rewritten at compile time to `_cell_<id>_<name>`, so it cannot be referenced from another cell -- including from `app.setup` to a regular cell. Two flavors hit fgx specifically: (1) an `@app.function` helper like `_client` called from another `@app.function` like `fetch_tsv` fails with `NameError: name '_cell_<id>_client' is not defined`; (2) a setup-block constant like `_DEFAULT_TOKEN` referenced from a cell's `mo.ui.text(value=_DEFAULT_TOKEN)` raises `NameError` and cascades into "ancestor raised an exception" failures down the cell graph. The fix is uniform: drop the leading underscore for anything visible across cells (`client`, `DEFAULT_TOKEN`). Reserve `_x` for names you genuinely want cell-private -- scratch in a single cell. `_env_path` defined and used only inside `app.setup` is fine; `_DEFAULT_TOKEN` defined in setup and read elsewhere is not.
- **Auto-save races when the user is typing in the browser.** A live marimo browser session writes the `.py` file every time the user edits a cell. If you `Write`/`Edit` the same file from CLI at the same moment, marimo's saver may interleave the two and produce file states like `app._unparsable_cell(...)` wrapping half-typed code, or move `with app.setup:` into a cell, or split `@app.function` decorators across cells. When this happens, don't try to surgical-edit the broken state -- re-`Write` the file fresh from a clean version, then ask the user to "Run all" to reload. If you're doing more than a one-line tweak while the browser is open, ask the user to switch tabs / pause typing for a beat first.
- **TSV is the API default.** `Content-Type: text/tab-separated-values`; pass `format=json` only when the response is a nested object (e.g. `colocalization_by_variant`). `fetch_tsv` adds nothing; `fetch_json` *unconditionally* injects `format=json` into the query string -- if you pass `format="tsv"` to `fetch_json` it gets silently overridden. Pick the right helper for the response shape you want, don't try to flip format via kwargs.
- **Response sizes can blow up context.** `credible_sets_by_gene/PCSK9` returns ~3,200 rows, ~1.3 MB. `resource_metadata/finngen` returns ~3,300 rows. Filter or summarize in-cell (`pl.DataFrame.filter`, `.head(20)`, `group_by + agg`) before printing the result -- and never let raw API output land in the chat transcript at full size.
- **No data caching.** Every notebook run hits the API live. If a parameter is fixed (e.g. you're iterating on the plotting code, not the gene), bind the dataframe to a variable in an upstream cell so reactive cells downstream don't re-fetch.
- **`exome_results_by_*` returns per-variant exome scores, not collapsed burden tests.** The path name reads like "gene-level burden" but each row is one (variant, trait) pair from genebass single-variant association testing -- there's no SKAT/SKAT-O/burden statistic in the response. If you need a per-gene rare-variant headline (one beta per gene-trait), aggregate the variant rows yourself: filter to `pLoF` and/or `missense` annotation classes, drop rows where `mlog10p` saturates at the floating-point cap (`mlog10p == 324` with `se == 0`), and either pick the top variant or do a meta-analytic combine. Don't promise the user "burden test results" in the notebook intro -- frame it as "per-variant exome scores" and let the reader see the allelic series. nb04 is the worked example.
- **Lead-variant logic is per-context.** nb01 picks the highest-PIP missense in the GWAS bucket for PCSK9; nb03 picks the highest-`mlog10p` variant per `cs_id`. Neither is universally right. For a new gene or phenotype, pick the strategy that matches the question: "what's the canonical coding hit?" -> missense filter; "what's the lead variant in each independent locus?" -> per-`cs_id` top-mlog10p. Don't reuse one notebook's heuristic in another without thinking.
- **Phenotype codes are FinnGen-specific.** `I9_CHD`, `T2D`, `K11_CROHN` etc. come from FinnGen's endpoint definitions, not free text. Hit `fetch_json("/trait_name_mapping")` once to get the full code -> name mapping, then pick the code you actually want before passing strings to phenotype-keyed paths. UKBB-style trait codes (`30640` ApoB, `30780` LDL-C, `6177_1` cholesterol-lowering medication) come from genebass / Open Targets and have their own conventions; don't mix them with FinnGen codes in the same prompt.
- **rsID -> variant resolution can return the wrong allele.** `/rsid/variants?rsid=...` returns chr:pos:ref:alt, but FinnGen's credible-set store may have indexed the variant under a different alt allele than the resolver returns. nb02 falls back across alt alleles to handle this; if you see 0 rows for a variant that should clearly have hits, try the other alt allele before giving up.

## When to revise this skill

Two thresholds, with different responses:

- **A new notebook lands.** Add a row to the catalog table above (file, entry point, endpoint chain, viz). The Portfolio principle should already have flagged this as a deliberate axis-add, not a duplicate; if it's a duplicate, push back before merging.
- **A function defined in nbN is imported by nbM.** Replace the four-column catalog table (entry / chain / viz) with a per-module table that lists each notebook's `@app.function` helpers and what they do, the way [jx's compose-notebook](https://github.com/broadinstitute/jx/blob/main/.claude/skills/compose-notebook/SKILL.md) does it. Until then, the helpers stay duplicated by design -- they're tiny and self-contained -- and the catalog stays compact.
