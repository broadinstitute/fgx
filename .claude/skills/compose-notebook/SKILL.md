---
name: compose-notebook
description: Compose a new genetics analysis in the fgx repo by reusing the marimo notebook(s) under notebooks/ as the substrate, hitting FinnGenie's REST API directly via httpx.get. Trigger when the user asks for a notebook, analysis, figure, or vignette against human-genetics data exposed by FinnGenie -- "credible sets near gene X", "loci driving phenotype Y", "what colocalizes at variant Z", "pull the GWAS signal for...", "compare PIPs across resources", "PheWAS for rsid", "look up an rsID", "search for a phenotype", or any composition question against /api/v1/* -- credible_sets_by_{gene,phenotype,variant,qtl_gene,region,id} and credible_sets/{id}/stats; colocalization_by_{variant,credible_set_id}; exome_results_by_{gene,phenotype,region,variant}; summary_stats/{resource}/{data_type}; nearest_genes; peak_to_genes; gene_disease; rsid/variants; search; phenotype/{resource}/{phenocode}/markdown; resource_metadata; resources; datasets. Authoritative path list at /api/v1/openapi.json (browseable at /api/v1/docs); confirm against the spec rather than guessing if the user asks about an endpoint not listed here. Use this instead of writing standalone httpx code from scratch when the existing helpers already cover the call shape, and instead of asking the user to wire auth themselves.
---

# Compose a new marimo notebook in fgx

## What this skill is for

fgx is a single-substrate repo: marimo notebooks against FinnGenie's REST API at `https://finngenie.fi/api/v1/*` via `httpx.get`, and nothing else. There is no Python SDK, no MCP server, no schema cache, and no shared library file -- every helper lives in the notebook that introduced it, and other notebooks reach across via plain Python imports. The `nbNN_<topic>.py` filename convention exists *specifically* so a sibling notebook can do `sys.path.insert(0, "notebooks") + from nb01_pcsk9_walkthrough import fetch_tsv`. Python forbids module names starting with a digit; `nbNN_*` sidesteps that. This is the same layering the [jx repo](https://github.com/broadinstitute/jx) uses for JUMP Cell Painting.

When the user asks a genetics question that the FinnGenie API can answer, your job is to (a) decide whether to edit a canonical notebook in place (parameter swap), or compose a fresh exploration notebook that imports from the canonical homes, and (b) reuse the existing helpers rather than re-implementing auth, TSV/JSON parsing, or any other plumbing.

## The catalog

nb01-04 are the canonical core: each is a self-contained vignette that demonstrates one entry-point + visualization combination AND owns the helpers it introduced. Future explorations are nb05+ -- they import from this core and are not expected to back-port their own helpers unless a third notebook ends up needing them.

| Module | Reusable functions | What the vignette demonstrates |
|---|---|---|
| `nb01_pcsk9_walkthrough` | `client()`, `fetch_tsv(path, **params)`, `fetch_json(path, **params)` | Gene -> `credible_sets_by_gene` -> `colocalization_by_variant` -> bar chart of top traits + colocalizing pairs. **Every other notebook imports the three helpers from here.** |
| `nb02_variant_phewas` | `alt_alleles(variant)` | rsID / variant -> `rsid/variants` -> `credible_sets_by_variant` -> `nearest_genes` -> PheWAS dot plot. `alt_alleles` returns strand/allele-fallback candidates when an rsID resolution misses the credible-set index under one alt allele. |
| `nb03_phenotype_locus_zoom` | `pick_leads(cs)`, `annotate_with_nearest_gene(variants)` | Phenocode + resource -> `credible_sets_by_phenotype` -> Manhattan with gene labels. `pick_leads` reduces a credible-set DataFrame to one lead row per `cs_id`; `annotate_with_nearest_gene` joins each variant with its nearest protein-coding gene. |
| `nb04_gene_exome_burden` | `prepare_deleterious(exome, ci_z=1.96)` | Gene (rare-variant arm) -> `exome_results_by_gene` -> `gene_disease` -> forest plot of pLoF/missense betas with CIs. `prepare_deleterious` is the polars prep pipeline (filter to pLoF/missense, drop `mlog10p` underflow rows, attach CI bounds + `trait_variant` label). nb05 imports it to avoid duplicating the pipeline. |
| `nb05_pign_cdg` | -- | Gene (rare-variant arm, recessive Mendelian) -> `exome_results_by_gene` -> `gene_disease` -> forest plot + curated cross-check. Same shape as nb04, different story: where PCSK9's two arms *converge* on lipids, PIGN's *diverge* -- gencc/monarch converge on MCAHS1 ("Definitive" by ClinGen and G2P, autosomal recessive) while genebass surfaces corneal biomechanics in adult heterozygotes. UKB doesn't enroll affected MCAHS1 probands, so the rare-variant arm shows what one PIGN dose perturbs in carriers, not the recessive disease itself. First concrete demo of the cross-notebook composition pattern -- imports `prepare_deleterious` from nb04 rather than duplicating the polars chain. |
| `nb06_variant_pqtl_function` | `pqtl_credible_sets(variant)`, `direction_consensus(df, beta_col)` | Variant -> `credible_sets_by_variant` filtered to `data_type == "pQTL"` -> per-protein PIP/beta panel + direction-of-effect summary. Hero replay of the FinnGenie demo (Karjalainen, 2026-05-05): `chr2:9521321:A:G` in ADAM17 -> 4 plasma proteins all with negative beta -> consistent loss-of-sheddase mechanism. Distinct from nb01 (gene-first) and nb02 (full PheWAS) by being variant-first AND pQTL-only AND by surfacing sign-of-effect across the panel rather than top-mlog10p across all data types. |
| `nb07_data_catalog` | `list_datasets()`, `list_resources()`, `resource_metadata(resource)` | `/datasets` + `/resources` + `/resource_metadata/{resource}` -> grouped catalog table with `mo.ui.table` selection driving a per-resource metadata drill-down. The "what's available?" introspection notebook -- read it first when a new question lands to confirm a dataset covers it. Replays the slide-05/06 catalog walk from the FinnGenie demo. |

When the user's question matches a row's "demonstrates" column with only a parameter changed, see Path A below. When it matches the *shape* of a row but tells a different story, or when it doesn't match any row, compose -- see Path C.

## Cross-notebook import recipe

Every non-nb01 notebook's `with app.setup:` block ends with this snippet:

```python
NOTEBOOK_DIR = Path(__file__).resolve().parent
if str(NOTEBOOK_DIR) not in sys.path:
    sys.path.insert(0, str(NOTEBOOK_DIR))

from nb01_pcsk9_walkthrough import client, fetch_json, fetch_tsv  # noqa: F401
```

The `noqa: F401` matters because ruff sees the imported names as unused -- `@app.function` and `@app.cell` cells reference them statically-invisibly. Importing nb01 also runs nb01's own `with app.setup:` block (idempotent: same `load_dotenv`, same `import polars`); cell *bodies* don't execute on import because `if __name__ == "__main__": app.run()` guards the kernel start. Local `BASE` and `FINNGENIE_TOKEN` constants stay in each notebook's setup so .env diagnostics surface locally rather than chasing imports.

## Two paths: parameter swap vs compose

### Path A: edit a canonical notebook in place

Use this when the new question has the **same shape** as nb01-04 with only a parameter changed:

- "do nb01 for APOE" -> change `GENE = "PCSK9"` to `GENE = "APOE"` in nb01.
- "do nb02 for `rs7412`" -> change `RSID` in nb02.
- "do nb03 for T2D on UKBB" -> change `PHENOTYPE` and `RESOURCE` in nb03.
- "do nb04 for LDLR" -> change `GENE` in nb04.

Marimo's reactivity reruns dependent cells; don't manually kick anything. This is the cheapest path. Use it sparingly though -- editing a canonical in place overwrites its narrative and the next reader will see the new gene's story instead of the original. If the swap *also* needs different markdown to make sense, consider Path C instead.

### Path C: compose a new exploration notebook

This is the default for everything that isn't a parameter-swap. Create `notebooks/nbNN_<topic>.py`. The starting template:

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
    FINNGENIE_TOKEN = os.environ.get("FINNGENIE_TOKEN")
    BASE = "https://finngenie.fi/api/v1"

    NOTEBOOK_DIR = Path(__file__).resolve().parent
    if str(NOTEBOOK_DIR) not in sys.path:
        sys.path.insert(0, str(NOTEBOOK_DIR))

    from nb01_pcsk9_walkthrough import client, fetch_json, fetch_tsv  # noqa: F401


@app.cell
def _():
    mo.md(r"""
    # nbNN: <topic>
    ...
    """)
    return


# analysis cells: call fetch_tsv / fetch_json directly, no boilerplate.


if __name__ == "__main__":
    app.run()
```

Keep the `marimo<0.23.4` and `jedi<0.20.0` pins verbatim -- both work around `uvx marimo edit --sandbox` resolver bugs. If you trim them, the next launch will fail at venv provisioning, not at runtime.

If the new notebook genuinely needs a fresh helper -- e.g. a region-around-gene utility, a coloc-row deduper -- define it as a top-level `@app.function` *in nbNN itself*. It then becomes importable by future composers via the same `from nbNN_<topic> import <helper>` pattern. **No central library, ever.** Helpers stay where they were born; importers reach across notebooks.

## Process for a new composition

1. **Read the closest existing notebook.** The first cell explains what it does; the helpers and import recipe sit at the top. Skim before composing.
2. **Identify the endpoint.** Map the user's question to a `/api/v1/*` path. The full surface is in the OpenAPI spec at <https://finngenie.fi/api/v1/openapi.json> (browseable Swagger at <https://finngenie.fi/api/v1/docs>) -- 28 operations across 13 tags. Read the spec before guessing; the description's path list is a reminder, not a contract. If you're unsure of the response shape, run `curl -fsS -H "Authorization: Bearer $FINNGENIE_TOKEN" "https://finngenie.fi/api/v1/<path>?format=json"` once to confirm, then write the cell.
3. **Pick the path.** Parameter-swap on a canonical -> Path A. Anything else -> Path C.
4. **Validate against the live kernel via marimo-pair.** This is non-negotiable, and the rest of this step exists because there are three failure modes that look like shortcuts but aren't.
    - **Browser session is a prerequisite, not optional.** `uvx marimo edit --headless` launches a kernel, but `execute-code.sh` won't attach until the user opens the URL in a browser tab and a websocket session is alive. If you call `execute-code.sh` and see `No active sessions on the server`, the right move is to tell the user "open `http://127.0.0.1:<port>/` in your browser, then ping me" -- not to fall back to running the same logic in a standalone `uv run python3` script. A parallel script imports nothing from the notebook and runs with its own module scope, so it can't catch `@app.function` vs `@app.cell` scoping bugs, marimo auto-formatter rewrites, or stale compiled cells. False green there is worse than no signal.
    - **Hot-reload reloads the file but does not re-run cells.** When you edit a `.py` file, marimo notices and re-imports module-level definitions, but the kernel still holds the old compiled cell objects in memory until something triggers a rerun. So after every code edit, ask the user to hit "Run all" (Cmd+Shift+R, or the play-all button in the toolbar) before you re-validate. Otherwise `execute-code.sh` will report errors against the *previous* version of the code, which is maddening to debug.
    - **"Run all" does NOT re-execute `with app.setup:` or re-register `@app.function` defs.** Both run once at kernel-start and stay registered for the kernel's lifetime. "Run all" only re-runs `@app.cell` bodies. So if you edit imports, the setup block, or any `@app.function` helper, the kernel still holds the *old* registrations and `fetch_tsv.__module__` (or whatever) reports stale results. The fix is a full kernel restart -- kill the marimo server (`kill $(lsof -ti :<port>)`), relaunch with `uvx marimo edit --sandbox`, and have the user refresh the browser tab. Symptom that you've hit this: a function imported from another notebook reports `__module__ == "__main__"` instead of the expected sibling-notebook name -- the local definition is still in memory from the pre-edit kernel.
    - **One cell at a time.** Once the session is attached and cells have been re-run, drive validation through `execute-code.sh --port <port> -c "..."` cell by cell -- not by re-pasting whole notebook bodies. The setup block (`with app.setup:`) is already provisioned by the PEP 723 header, and `fetch_tsv` / `fetch_json` are visible in the scratchpad's namespace via the import recipe; a single cell returns a `polars.DataFrame` you can summarize inline.
5. **End with extension prompts.** Every vignette closes with a `## To extend` markdown cell listing 2-3 concrete next prompts (e.g. "swap PCSK9 for LDLR", "add a second resource and stack the dot plots"). This turns each notebook into a launchpad rather than a dead end.
6. **Plot when it adds signal.** Altair is in the PEP 723 deps; a one-cell chart is often more informative than a 20-row table head.
7. **Don't fabricate.** If a helper doesn't exist or an endpoint returns an empty result, say so and ask the user how to proceed -- don't invent a fallback path.

## Generate a session snapshot for the molab preview

Once the notebook runs end-to-end, commit a session snapshot so the molab badge in the README serves cached outputs instead of triggering live re-execution. One-liner:

```bash
env -u PYTHONPATH uvx marimo export session --sandbox notebooks/nbNN_<topic>.py
```

This executes the notebook in an isolated venv (deps from the PEP 723 header) and writes `notebooks/__marimo__/session/nbNN_<topic>.py.json` -- ~40k of cell outputs that the marimo viewer can replay statelessly. Don't reach for the alternative (`marimo edit` + browser + Run all + autosave) -- it requires a live websocket session and is much heavier than this one-shot. Run it after every meaningful change to the notebook so the committed preview tracks the source. Skip if the notebook hits an endpoint that's slow or rate-limited; the cell will execute again on each export and a stale snapshot is worse than no snapshot for a flaky path.

## Gotchas

The first three are tooling traps that cost real session time when missed; the next three are API papercuts; the rest are biology / endpoint-semantics nuances.

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

- **A new notebook lands.** Add a row to the catalog table. If the notebook introduces a reusable `@app.function` helper, fill the "Reusable functions" column. If it just composes existing helpers and adds a vignette, leave that column with `--`. The Portfolio principle stays implicit: if a *third* notebook reaches for a helper that lives in nbN, that helper has earned a row entry. If `from nb05_... import ...` ever appears in nb06+, update nb05's row to list the imported helpers.
- **A canonical helper migrates between notebooks.** If `client / fetch_tsv / fetch_json` ever moves out of nb01 (e.g. nb01 gets renamed or split), update the import recipe snippet in this skill to point at the new home -- and update the import lines in every importer. There's no central library to grep, but there are only a handful of importers; `grep -rn "from nb01_" notebooks/` finds them all.
- **A central library actually becomes warranted.** This skill currently *commits* to "no library." If five notebooks start importing the same three helpers from nb01 and nb01 itself becomes incidental scaffolding for those helpers, that's the moment to revisit -- and at that point, port to a `notebooks/_fgx.py` plus a single import line, not to a packaged `fgx/` directory. Keep the bar high; the cross-notebook import works fine for the foreseeable future.
