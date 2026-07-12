---
name: compose-notebook
description: Compose a new genetics analysis in the fgx repo by reusing the marimo notebook(s) under notebooks/ as the substrate, hitting GeneGenie's REST API directly via httpx.get. Trigger when the user asks for a notebook, analysis, figure, or vignette against human-genetics data exposed by GeneGenie -- "credible sets near gene X", "loci driving phenotype Y", "what colocalizes at variant Z", "pull the GWAS signal for...", "compare PIPs across resources", "PheWAS for rsid", "look up an rsID", "search for a phenotype", or any composition question against /api/v1/* -- credible_sets_by_{gene,phenotype,variant,qtl_gene,region,id} and credible_sets/{id}/stats; colocalization_by_{variant,credible_set_id}; exome_results_by_{gene,phenotype,region,variant}; summary_stats/{resource}/{data_type}; nearest_genes; peak_to_genes; gene_disease; rsid/variants; search; phenotype/{resource}/{phenocode}/markdown; resource_metadata; resources; datasets. Authoritative path list at /api/v1/openapi.json (browseable at /api/v1/docs); confirm against the spec rather than guessing if the user asks about an endpoint not listed here. Use this instead of writing standalone httpx code from scratch when the existing helpers already cover the call shape, and instead of asking the user to wire auth themselves.
---

# Compose a new marimo notebook in fgx

## What this skill is for

fgx is a single-substrate repo: marimo notebooks against GeneGenie's REST API at `https://genegenie.broadinstitute.org/api/v1/*` via `httpx.get`, and nothing else.
There is no Python SDK, no MCP server, no schema cache, and no shared library file -- every helper lives in the notebook that introduced it, and other notebooks reach across via plain Python imports.
The `nbNN_<topic>.py` filename convention exists *specifically* so a sibling notebook can do `sys.path.insert(0, "notebooks") + from nb01_pcsk9_walkthrough import fetch_tsv`.
Python forbids module names starting with a digit; `nbNN_*` sidesteps that.
This is the same layering the [jx repo](https://github.com/broadinstitute/jx) uses for JUMP Cell Painting.

When the user asks a genetics question that the GeneGenie API can answer, your job is to (a) decide whether to edit a canonical notebook in place (parameter swap), or compose a fresh exploration notebook that imports from the canonical homes, and (b) reuse the existing helpers rather than re-implementing auth, TSV/JSON parsing, or any other plumbing.

## The catalog

nb01-04 are the canonical core: each is a self-contained vignette that demonstrates one entry-point + visualization combination AND owns the helpers it introduced.
Future explorations are nb05+ -- they import from this core and are not expected to back-port their own helpers unless a third notebook ends up needing them.

| Module | Reusable functions | What the vignette demonstrates |
|---|---|---|
| `nb01_pcsk9_walkthrough` | `client()`, `fetch_tsv(path, **params)`, `fetch_json(path, **params)` | Gene -> `credible_sets_by_gene` -> `colocalization_by_variant` -> bar chart of top traits + colocalizing pairs. **Every other notebook imports the three helpers from here.** |
| `nb02_variant_phewas` | `alt_alleles(variant)` | rsID / variant -> `rsid/variants` -> `credible_sets_by_variant` -> `nearest_genes` -> PheWAS dot plot. `alt_alleles` returns strand/allele-fallback candidates when an rsID resolution misses the credible-set index under one alt allele. |
| `nb03_phenotype_locus_zoom` | `pick_leads(cs)`, `annotate_with_nearest_gene(variants)` | Phenocode + resource -> `credible_sets_by_phenotype` -> Manhattan with gene labels. `pick_leads` reduces a credible-set DataFrame to one lead row per `cs_id`; `annotate_with_nearest_gene` joins each variant with its nearest protein-coding gene. |
| `nb04_gene_exome_burden` | `prepare_deleterious(exome, ci_z=1.96)` | Gene (rare-variant arm) -> `exome_results_by_gene` -> `gene_disease` -> forest plot of pLoF/missense betas with CIs. `prepare_deleterious` is the polars prep pipeline (filter to pLoF/missense, drop `mlog10p` underflow rows, attach CI bounds + `trait_variant` label). nb05 imports it to avoid duplicating the pipeline. |
| `nb05_pign_cdg` | -- | Gene (rare-variant arm, recessive Mendelian) -> `exome_results_by_gene` -> `gene_disease` -> forest plot + curated cross-check. Same shape as nb04, different story: where PCSK9's two arms *converge* on lipids, PIGN's *diverge* -- gencc/monarch converge on MCAHS1 ("Definitive" by ClinGen and G2P, autosomal recessive) while genebass surfaces corneal biomechanics in adult heterozygotes. UKB doesn't enroll affected MCAHS1 probands, so the rare-variant arm shows what one PIGN dose perturbs in carriers, not the recessive disease itself. First concrete demo of the cross-notebook composition pattern -- imports `prepare_deleterious` from nb04 rather than duplicating the polars chain. |
| `nb06_variant_pqtl_function` | `pqtl_credible_sets(variant)`, `direction_consensus(df, beta_col)` | Variant -> `credible_sets_by_variant` filtered to `data_type == "pQTL"` -> per-protein PIP/beta panel + direction-of-effect summary. Hero replay of the GeneGenie demo (Karjalainen, 2026-05-05): `chr2:9521321:A:G` in ADAM17 -> 4 plasma proteins all with negative beta -> consistent loss-of-sheddase mechanism. Distinct from nb01 (gene-first) and nb02 (full PheWAS) by being variant-first AND pQTL-only AND by surfacing sign-of-effect across the panel rather than top-mlog10p across all data types. |
| `nb07_data_catalog` | `list_datasets()`, `list_resources()`, `resource_metadata(resource)` | `/datasets` + `/resources` + `/resource_metadata/{resource}` -> grouped catalog table with `mo.ui.table` selection driving a per-resource metadata drill-down. The "what's available?" introspection notebook -- read it first when a new question lands to confirm a dataset covers it. Replays the slide-05/06 catalog walk from the GeneGenie demo. |
| `nb08_genetics_primer` | -- | Educational walkthrough: LDLR as a case study through GWAS, fine-mapping, colocalization, exome, eQTL/pQTL, and gene-disease curation. Imports `prepare_deleterious` from nb04. Teaches genetics concepts by showing what a clean signal looks like for a well-understood gene. |
| `nb09_polygenic_heart_disease` | -- | Phenotype (CHD) -> `credible_sets_by_phenotype` -> Manhattan + effect-size histogram + cross-pathway colocalization + multi-gene exome forest plot. Builds the case that heart disease is polygenic from four layers of evidence: many independent loci, small individual effects, diverse colocalizing traits (with CHD self-coloc filtering), and rare coding variants in LDLR/PCSK9/APOB/LPA. Imports `pick_leads` and `annotate_with_nearest_gene` from nb03, `prepare_deleterious` from nb04. |
| `nb10_diabetes_susceptibility` | -- | Phenotype (T2D) -> `credible_sets_by_phenotype` -> Manhattan + effect-size histogram + TCF7L2 colocalization deep dive + MODY gene-disease curation. Same polygenic shape as nb09 but with a distinct angle: the GWAS loci overlap heavily with monogenic MODY genes (HNF1A, HNF4A, GCK, KCNJ11, ABCC8, PPARG), showing that common variants whisper what rare mutations shout. Imports `pick_leads` from nb03. First notebook to use the `alt.Data(values=df.to_dicts())` pattern for altair chart data — see the polars/altair gotcha. |

When the user's question matches a row's "demonstrates" column with only a parameter changed, see Path A below.
When it matches the *shape* of a row but tells a different story, or when it doesn't match any row, compose -- see Path C.

## Cross-notebook import recipe

Every non-nb01 notebook's `with app.setup:` block ends with this snippet:

```python
NOTEBOOK_DIR = Path(__file__).resolve().parent
if str(NOTEBOOK_DIR) not in sys.path:
    sys.path.insert(0, str(NOTEBOOK_DIR))

from nb01_pcsk9_walkthrough import client, fetch_json, fetch_tsv
```

Ruff's static analysis can't see all marimo references across `with app.setup`, `@app.function`, `@app.cell`, and cross-notebook imports.
Rather than peppering per-line `# noqa` comments, fgx keeps notebook-specific ignores centralized in `[tool.ruff.lint.per-file-ignores]` in `pyproject.toml`.
Per-line comments are fragile -- a manual `ruff check --fix --unsafe-fixes` run, or an aggressive IDE save action, can strip them and re-trigger warnings on the next lint run.
The pyproject-level rule survives both.

Importing nb01 also runs nb01's own `with app.setup:` block (idempotent: same `load_dotenv`, same `import polars`); cell *bodies* don't execute on import because `if __name__ == "__main__": app.run()` guards the kernel start.
Local `BASE` and `GENEGENIE_TOKEN` constants stay in each notebook's setup so .env diagnostics surface locally rather than chasing imports.

## Two paths: parameter swap vs compose

### Path A: edit a canonical notebook in place

Use this when the new question has the **same shape** as nb01-04 with only a parameter changed:

- "do nb01 for APOE" -> change `GENE = "PCSK9"` to `GENE = "APOE"` in nb01.
- "do nb02 for `rs7412`" -> change `RSID` in nb02.
- "do nb03 for T2D on UKBB" -> change `PHENOTYPE` and `RESOURCE` in nb03.
- "do nb04 for LDLR" -> change `GENE` in nb04.

Marimo's reactivity reruns dependent cells; don't manually kick anything.
This is the cheapest path.
Use it sparingly though -- editing a canonical in place overwrites its narrative and the next reader will see the new gene's story instead of the original.
If the swap *also* needs different markdown to make sense, consider Path C instead.

### Path C: compose a new exploration notebook

This is the default for everything that isn't a parameter-swap.
Create `notebooks/nbNN_<topic>.py`.
The starting template:

```python
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo<0.23.4",
#     "jedi<0.20.0",
#     "polars",
#     "httpx",
#     "altair",
#     "python-dotenv",
# ]
# ///

import marimo

__generated_with = "0.23.3"
app = marimo.App(width="medium")

with app.setup:
    import os
    import sys
    from pathlib import Path

    import altair as alt
    import httpx  # only if you actually call httpx.X directly
    import marimo as mo
    import polars as pl
    from dotenv import load_dotenv

    load_dotenv(Path(__file__).resolve().parent.parent / ".env")
    GENEGENIE_TOKEN = os.environ.get("GENEGENIE_TOKEN")
    BASE = "https://genegenie.broadinstitute.org/api/v1"

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    from nb01_pcsk9_walkthrough import client, fetch_json, fetch_tsv


@app.cell
def _():
    mo.md(r"""
    # nbNN: <topic>
    ...
    """)
    return


# Pattern for cells that combine narrative with data output:
# WRONG - intro markdown is silently dropped:
#   mo.md("## Section title\n\nExplanation...")
#   result = fetch_tsv(...)
#   mo.vstack([mo.md("### Subtitle"), result])  # only this renders
#
# RIGHT - everything in one mo.vstack:
#   result = fetch_tsv(...)
#   mo.vstack([
#       mo.md("## Section title\n\nExplanation..."),
#       mo.md("### Subtitle"),
#       result,
#   ])
#
# Also RIGHT - split into separate cells:
#   Cell 1: mo.md("## Section title\n\nExplanation...")
#   Cell 2: result = fetch_tsv(...); mo.vstack([mo.md("### Subtitle"), result])

# Use unique variable names across cells (trait_chart, coloc_chart, not chart).


if __name__ == "__main__":
    app.run()
```

Keep the `marimo<0.23.4` and `jedi<0.20.0` pins verbatim -- both work around `uvx marimo edit --sandbox` resolver bugs.
If you trim them, the next launch will fail at venv provisioning, not at runtime.

If the new notebook genuinely needs a fresh helper -- e.g. a region-around-gene utility, a coloc-row deduper -- define it as a top-level `@app.function` *in nbNN itself*.
It then becomes importable by future composers via the same `from nbNN_<topic> import <helper>` pattern.
**No central library, ever.** Helpers stay where they were born; importers reach across notebooks.

## Process for a new composition

1. **Read the closest existing notebook.** The first cell explains what it does; the helpers and import recipe sit at the top.
   Skim before composing.
2. **Identify the endpoint.** Map the user's question to a `/api/v1/*` path.
   The full surface is in the OpenAPI spec at <https://genegenie.broadinstitute.org/api/v1/openapi.json> (browseable Swagger at <https://genegenie.broadinstitute.org/api/v1/docs>) -- 28 operations across 13 tags.
   Read the spec before guessing; the description's path list is a reminder, not a contract.
   If you're unsure of the response shape, run `curl -fsS -H "Authorization: Bearer $GENEGENIE_TOKEN" "https://genegenie.broadinstitute.org/api/v1/<path>?format=json"` once to confirm, then write the cell.
3. **Pick the path.** Parameter-swap on a canonical -> Path A. Anything else -> Path C.
4. **See the data before writing the story.** Run the key API calls and inspect the actual values *before* drafting narrative cells.
   Characterize what the data contains: unique value counts, group-by distributions, range of the column you plan to sort or color by.
   If the notebook's thesis is "X is diverse / spread / multi-pathway," confirm the data actually shows diversity — not the same category repeated under different labels.
   If the data doesn't support the planned narrative, either adjust the query (filter, deduplicate, sort by a different column, widen the search) or rewrite the story to match what the data actually says.
   Narratives written before seeing data are hypotheses; narratives written after are observations.
   The notebook should contain observations.
5. **Validate against the live kernel via marimo-pair.** This is non-negotiable, and the rest of this step exists because there are three failure modes that look like shortcuts but aren't. - **Browser session is a prerequisite, not optional.** `uvx marimo edit --headless` launches a kernel, but `execute-code.sh` won't attach until the user opens the URL in a browser tab and a websocket session is alive.
   If you call `execute-code.sh` and see `No active sessions on the server`, the right move is to tell the user "open `http://127.0.0.1:<port>/` in your browser, then ping me" -- not to fall back to running the same logic in a standalone `uv run python3` script.
   A parallel script imports nothing from the notebook and runs with its own module scope, so it can't catch `@app.function` vs `@app.cell` scoping bugs, marimo auto-formatter rewrites, or stale compiled cells.
   False green there is worse than no signal. - **Hot-reload reloads the file but does not re-run cells.** When you edit a `.py` file, marimo notices and re-imports module-level definitions, but the kernel still holds the old compiled cell objects in memory until something triggers a rerun.
   So after every code edit, ask the user to hit "Run all" (Cmd+Shift+R, or the play-all button in the toolbar) before you re-validate.
   Otherwise `execute-code.sh` will report errors against the *previous* version of the code, which is maddening to debug. - **"Run all" does NOT re-execute `with app.setup:` or re-register `@app.function` defs.** Both run once at kernel-start and stay registered for the kernel's lifetime. "Run all" only re-runs `@app.cell` bodies.
   So if you edit imports, the setup block, or any `@app.function` helper, the kernel still holds the *old* registrations and `fetch_tsv.__module__` (or whatever) reports stale results.
   The fix is a full kernel restart -- kill the marimo server (`kill $(lsof -ti :<port>)`), relaunch with `uvx marimo edit --sandbox`, and have the user refresh the browser tab.
   Symptom that you've hit this: a function imported from another notebook reports `__module__ == "__main__"` instead of the expected sibling-notebook name -- the local definition is still in memory from the pre-edit kernel. - **One cell at a time.** Once the session is attached and cells have been re-run, drive validation through `execute-code.sh --port <port> -c "..."` cell by cell -- not by re-pasting whole notebook bodies.
   The setup block (`with app.setup:`) is already provisioned by the PEP 723 header, and `fetch_tsv` / `fetch_json` are visible in the scratchpad's namespace via the import recipe; a single cell returns a `polars.DataFrame` you can summarize inline.
5. **End with extension prompts.** Every vignette closes with a `## To extend` markdown cell listing 2-3 concrete next prompts (e.g. "swap PCSK9 for LDLR", "add a second resource and stack the dot plots").
   This turns each notebook into a launchpad rather than a dead end.
6. **Match chart encodings to the narrative claim.** If the prose says "diverse" or "spread across," the chart must encode the *diversity dimension* visually — not a metric that happens to be uniform.
   A bar chart where every bar is the same length carries no information.
   Before choosing a chart type, ask: what is the visual variable the reader should compare?
   If it's categorical spread (which loci, which data types, which pathways), use color, position, or faceting on those categories — not a quantitative axis that's uniformly high.
   If every value in a quantitative column falls within a narrow band (e.g. all PP_H4 between 0.95 and 1.0), that column is not informative as an axis — it's a filter criterion, not an encoding.
7. **Plot when it adds signal.** Altair is in the PEP 723 deps; a one-cell chart is often more informative than a 20-row table head.
8. **Don't fabricate.** If a helper doesn't exist or an endpoint returns an empty result, say so and ask the user how to proceed -- don't invent a fallback path.
9. **Visually verify every chart, not just its presence.** The accessibility tree confirms a chart *renders* ("Vega visualization" node) but not that it's *readable*.
   After the notebook runs end-to-end, save each chart standalone from the kernel (`chart.save("/tmp/chart_name.html")`), open it in agent-browser, and screenshot it.
   Check: (a) bars/dots/lines are visually distinguishable, (b) colors in the legend actually appear in the chart, (c) axis ranges make the data variation visible, (d) the chart tells the story the caption claims.
   A chart that renders but shows uniform bars, a single color, or an axis that hides all variation is worse than no chart — it looks like it was never looked at.
10. **Verify structural completeness via HTML export.** After the notebook runs end-to-end, export it as static HTML and inspect the rendered accessibility tree to confirm every section's narrative text, tables, and charts actually appear. The command is:
   ```bash
   env -u PYTHONPATH uvx marimo export html --sandbox notebooks/nbNN_<topic>.py -o /tmp/nbNN_<topic>.html
   ```
   Then open the HTML in agent-browser and run `agent-browser snapshot` to get the full accessibility tree. Walk the tree and confirm: (a) every section heading appears, (b) every `mo.md()` intro is visible (not silently dropped by a subsequent expression), (c) every table has data rows, (d) every chart renders as a "Vega visualization" node. This catches the silent-drop bug from the cell-output gotcha above, which is invisible in the live editor (where the user scrolls past) but glaring in the exported HTML (where sections are simply missing). Do this before declaring the notebook done.

## Generate a session snapshot for the molab preview

Once the notebook runs end-to-end, commit a session snapshot so the molab badge in the README serves cached outputs instead of triggering live re-execution.
One-liner:

```bash
env -u PYTHONPATH uvx marimo export session --sandbox notebooks/nbNN_<topic>.py
```

This executes the notebook in an isolated venv (deps from the PEP 723 header) and writes `notebooks/__marimo__/session/nbNN_<topic>.py.json` -- ~40k of cell outputs that the marimo viewer can replay statelessly.
Don't reach for the alternative (`marimo edit` + browser + Run all + autosave) -- it requires a live websocket session and is much heavier than this one-shot.
Run it after every meaningful change to the notebook so the committed preview tracks the source.
Skip if the notebook hits an endpoint that's slow or rate-limited; the cell will execute again on each export and a stale snapshot is worse than no snapshot for a flaky path.

**Always export the snapshot AFTER the final source edit / `ruff format` pass.** Each cell in the snapshot carries a `code_hash` of the cell source it was generated from, and molab attaches the stored output to a source cell only if the hashes match -- otherwise the cell renders empty in the public preview.
A whitespace-only reformat is enough to invalidate every snapshot in the file, so if you change *anything* about the `.py` (including running `ruff format`) you have to regenerate the `.json` too.
Commit the refreshed snapshot in the same change that touched the source.

**Wrap altair charts in `mo.ui.altair_chart(...)` or molab will not render them.** A bare `chart` expression emits the raw vega-lite spec as `application/vnd.vegalite.v6+json`, which expects the viewer to ship its own vega-lite renderer.
The marimo live editor does; molab's static viewer does not bundle vega-lite v6, so the cell paints blank even though the spec is fully present in the snapshot.
Wrapping in `mo.ui.altair_chart(chart)` flips the output to `text/html` containing a `<marimo-vega>` custom element, and the marimo runtime that molab *does* load brings its own vega renderer with it.
The general rule: raw third-party mimetypes depend on the viewer's renderer registry, while marimo widgets carry their renderer with them.
If a chart shows up in the live editor but not in molab, the first thing to try is the widget wrap.

**Project the DataFrame to the columns the chart actually encodes before passing it to altair.** vega-lite embeds the entire input frame inline in the chart spec, including columns you never reference in `encode(...)` or in a `transform_density` / `transform_filter`.
A many-thousand-row frame with several long string columns easily blows past molab's `output_max_bytes` ceiling (default ~10 MB) and the cell renders as a "Your output is too large" callout instead of a chart.
The live editor's limit is higher, so the same notebook works locally and fails in molab.
Fix: `df.select(["x_col", "color_col", ...])` before `alt.Chart(df)`.
The chart is mathematically identical and the embedded spec shrinks by an order of magnitude.
Same principle for tooltips: a `tooltip=[...]` list on a many-row chart is dead weight in a static preview (no hover) and can be the difference between fitting and not.

## Gotchas

The first one is the most important and applies to every notebook. The next five are marimo cell-authoring traps that silently produce wrong output; then tooling traps; then API papercuts; then biology / endpoint-semantics nuances.

- **Narrative must follow from data, not precede it.** The most damaging composition bug isn't a code error — it's prose that claims something the data doesn't show.
  The failure pattern: decide on a story, write markdown that asserts the conclusion, write code that fetches data, and never check whether the output matches.
  Any time a cell's narrative makes a qualitative claim — "diverse," "spread across," "multiple pathways," "large effect," "concentrated" — the data in that cell or its chart must visibly support it.
  The fix is process, not code: after the fetch-and-filter pipeline runs, inspect the actual values (unique counts, group-by distributions, the column you plan to sort or color by) and confirm they match the story.
  Common traps: a sort-and-head(N) selection surfaces N copies of the dominant category, not N diverse categories; overlapping ontology codes inflate apparent variety while carrying the same signal; a single group's entries crowd out the rest when per-group counts are uneven.
  When the data contradicts the planned narrative, either fix the query (filter, deduplicate, re-sort, widen the search) or rewrite the prose. A chart that contradicts its own caption is worse than no chart.

- **Only the last expression before `return` renders as cell output.** If a cell has `mo.md(intro_text)` on line 1, then `mo.vstack([header, table])` on line 10, the intro markdown is evaluated but its result is silently discarded -- the reader never sees it. The fix: wrap everything in a single `mo.vstack([mo.md(intro_text), header, table])` so both the explanation and the data appear together, or split the intro into its own cell. This is the single most common composition bug. Every cell that mixes narrative markdown with a data output must use `mo.vstack()` to combine them into one rendered element.
- **Each top-level variable name must be unique across all cells.** Marimo tracks every name assigned in a cell; if two cells both define `chart = alt.Chart(...)`, the second cell errors with `MultipleDefinitionError` even though neither cell exports `chart` via `return`. The fix: give each chart a descriptive name scoped to its section -- `trait_chart`, `coloc_chart`, `eqtl_chart`, `forest`. This applies to any commonly reused name: `df`, `top`, `summary`, `chart`, `table`. Underscore-prefixed names (`_chart`) are cell-local and don't collide, but can't be referenced from other cells.
- **Validate via marimo-pair against a live kernel, never via a standalone `uv run python3` script.** A parallel script can't see `@app.function` mangling, can't catch auto-formatter rewrites, and gives you a false green that costs more than no signal.
  If `execute-code.sh` reports `No active sessions`, the fix is to ask the user to open the browser URL -- not to bypass.
- **Polars DataFrames must be converted to dicts before passing to altair.** Although altair 6+ nominally accepts polars DataFrames via narwhals, the resulting vega-lite spec uses a named-dataset reference (`"data": {"name": "data-abc123"}`) instead of inlining the values.
  The `mo.ui.altair_chart()` wrapper and marimo's static HTML export cannot resolve these named references — the chart renders axes and labels but **no marks** (no bars, no dots, no lines).
  The standalone `chart.save("file.html")` works because it embeds the full vega-lite renderer with dataset resolution, but marimo's `marimo-vega` custom element does not.
  The fix: always convert to dicts before passing to altair — `alt.Chart(alt.Data(values=df.select("col_a", "col_b").to_dicts()))`.
  This forces inline `"data": {"values": [...]}` in the spec, which every renderer can handle.
  The `.select()` projection before `.to_dicts()` is not optional — without it, vega-lite embeds the entire frame (all columns), which bloats the spec and can exceed molab's `output_max_bytes` limit.
  Never add pandas / pyarrow to the deps — the PEP 723 header doesn't include them, and `.to_pandas()` will raise `ModuleNotFoundError`.
- **`mo.ui.altair_chart()` can silently break bar charts with few data points.** When a chart has very few rows (e.g. 3-row bar chart from a `group_by().agg()`), wrapping in `mo.ui.altair_chart()` renders the chart frame (axes, title, labels) but no marks — the bars are invisible.
  The same chart renders correctly as a bare expression (without the wrapper) or via `chart.save()`.
  The fix: for charts with small datasets (under ~10 rows), output the bare chart object instead of wrapping in `mo.ui.altair_chart()`.
  Reserve `mo.ui.altair_chart()` for charts with enough data points that the interactive selection adds value (Manhattan plots, PheWAS dots, histograms).
  This trades molab static-preview rendering (which needs the wrapper) for actually showing data in the live editor (which is what the user sees).
  If you need both, test the wrapped version visually before committing — don't assume it works because the data is in the spec.
- **Marimo mangles any name starting with `_`.** Any name beginning with an underscore is treated as cell-local and rewritten at compile time to `_cell_<id>_<name>`, so it cannot be referenced from another cell -- including from `app.setup` to a regular cell.
  Two flavors hit fgx specifically: (1) an `@app.function` helper like `_client` called from another `@app.function` like `fetch_tsv` fails with `NameError: name '_cell_<id>_client' is not defined`; (2) a setup-block constant like `_DEFAULT_TOKEN` referenced from a cell's `mo.ui.text(value=_DEFAULT_TOKEN)` raises `NameError` and cascades into "ancestor raised an exception" failures down the cell graph.
  The fix is uniform: drop the leading underscore for anything visible across cells (`client`, `DEFAULT_TOKEN`).
  Reserve `_x` for names you genuinely want cell-private -- scratch in a single cell.
  `_env_path` defined and used only inside `app.setup` is fine; `_DEFAULT_TOKEN` defined in setup and read elsewhere is not.
- **Auto-save races when the user is typing in the browser.** A live marimo browser session writes the `.py` file every time the user edits a cell.
  If you `Write`/`Edit` the same file from CLI at the same moment, marimo's saver may interleave the two and produce file states like `app._unparsable_cell(...)` wrapping half-typed code, or move `with app.setup:` into a cell, or split `@app.function` decorators across cells.
  When this happens, don't try to surgical-edit the broken state -- re-`Write` the file fresh from a clean version, then ask the user to "Run all" to reload.
  If you're doing more than a one-line tweak while the browser is open, ask the user to switch tabs / pause typing for a beat first.
- **TSV is the API default.** `Content-Type: text/tab-separated-values`; pass `format=json` only when the response is a nested object (e.g.
  `colocalization_by_variant`).
  `fetch_tsv` adds nothing; `fetch_json` *unconditionally* injects `format=json` into the query string -- if you pass `format="tsv"` to `fetch_json` it gets silently overridden.
  Pick the right helper for the response shape you want, don't try to flip format via kwargs.
- **Response sizes can blow up context.** `credible_sets_by_gene/PCSK9` returns ~3,200 rows, ~1.3 MB.
  `resource_metadata/finngen` returns ~3,300 rows.
  Filter or summarize in-cell (`pl.DataFrame.filter`, `.head(20)`, `group_by + agg`) before printing the result -- and never let raw API output land in the chat transcript at full size.
- **No data caching.** Every notebook run hits the API live.
  If a parameter is fixed (e.g. you're iterating on the plotting code, not the gene), bind the dataframe to a variable in an upstream cell so reactive cells downstream don't re-fetch.
- **Don't infer response semantics from path names.** Endpoint paths like `exome_results_by_gene` sound like they return gene-level burden statistics, but the actual response is per-variant single-variant association rows.
  Read the actual response shape (hit the endpoint once, or check the OpenAPI spec) before writing code or narrative that assumes a particular granularity or aggregation level.
  The same principle applies to column names: `mlog10p` can saturate at floating-point limits (`mlog10p == 324` with `se == 0`), `pip` can be null, `most_severe` is a VEP prediction not a clinical classification.
  Frame your notebook's descriptions around what the data actually contains, not what you infer from the label.
- **Heuristics from one notebook don't transfer automatically.** Each notebook picks a filtering and sorting strategy matched to its question — nb01 picks the highest-PIP missense for PCSK9, nb03 picks the highest-`mlog10p` per credible set for a phenotype-wide Manhattan.
  When composing a new notebook, choose the heuristic that matches *this* question rather than copy-pasting from the nearest existing notebook.
  "What's the canonical coding hit?" needs a different selection than "what's the lead variant in each independent locus?" which is different from "which loci show cross-trait sharing?"
- **The `/search` endpoint uses `q`, not `query`.** The parameter name is `q` (required), with optional `types` (comma-separated: `phenotypes`, `genes`) and `limit`.
  The response is JSON (not TSV) — don't use `fetch_tsv` or `fetch_json` (which injects `format=json` and can cause 422s); call `client().get(f"{BASE}/search", params={"q": term})` directly and parse `r.json()`.
- **`/credible_sets_by_phenotype` requires `resource/phenotype` in the path.** The path is `/credible_sets_by_phenotype/{resource}/{phenotype_or_study}`, not `/credible_sets_by_phenotype/{phenotype}`.
  Omitting the resource gives a 404.
  Use the FinnGen resource name (`finngen`) for FinnGen phenotype codes (`T2D`, `I9_CHD`), or `open_targets` for GCST codes.
- **Identifiers are system-specific and may not resolve as you expect.** Phenotype codes (`I9_CHD`, `T2D`) are FinnGen endpoint definitions, not free text — use `fetch_json("/trait_name_mapping")` to look up the right code before passing strings to phenotype-keyed paths.
  UKBB-style trait codes come from a different namespace; don't mix them with FinnGen codes.
  rsIDs can resolve to the wrong alt allele: `/rsid/variants` returns chr:pos:ref:alt, but the credible-set store may have indexed a different alt allele.
  The general principle: any time you bridge between identifier systems (rsID → variant, gene name → phenotype, trait code → human-readable label), verify the resolution matched before building downstream logic on it.
  When a query returns 0 rows for an identifier you know should have hits, the identifier mapping — not the data — is usually the problem.

## When to revise this skill

- **A new notebook lands.** Add a row to the catalog table.
  If the notebook introduces a reusable `@app.function` helper, fill the "Reusable functions" column.
  If it just composes existing helpers and adds a vignette, leave that column with `--`.
  The Portfolio principle stays implicit: if a *third* notebook reaches for a helper that lives in nbN, that helper has earned a row entry.
  If `from nb05_... import ...` ever appears in nb06+, update nb05's row to list the imported helpers.
- **A canonical helper migrates between notebooks.** If `client / fetch_tsv / fetch_json` ever moves out of nb01 (e.g. nb01 gets renamed or split), update the import recipe snippet in this skill to point at the new home -- and update the import lines in every importer.
  There's no central library to grep, but there are only a handful of importers; `grep -rn "from nb01_" notebooks/` finds them all.
- **A central library actually becomes warranted.** This skill currently *commits* to "no library."
  If five notebooks start importing the same three helpers from nb01 and nb01 itself becomes incidental scaffolding for those helpers, that's the moment to revisit -- and at that point, port to a `notebooks/_fgx.py` plus a single import line, not to a packaged `fgx/` directory.
  Keep the bar high; the cross-notebook import works fine for the foreseeable future.
