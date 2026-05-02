---
name: compose-notebook
description: Compose a new genetics analysis in the fgx repo by reusing the marimo notebook(s) under notebooks/ as the substrate, hitting FinnGenie's REST API directly via httpx.get. Trigger when the user asks for a notebook, analysis, figure, or vignette against human-genetics data exposed by FinnGenie -- "credible sets near gene X", "loci driving phenotype Y", "what colocalizes at variant Z", "pull the GWAS signal for...", "compare PIPs across resources", "PheWAS for rsid", "look up an rsID", "search for a phenotype", or any composition question against /api/v1/* -- credible_sets_by_{gene,phenotype,variant,qtl_gene,region,id} and credible_sets/{id}/stats; colocalization_by_{variant,credible_set_id}; exome_results_by_{gene,phenotype,region,variant}; summary_stats/{resource}/{data_type}; nearest_genes; peak_to_genes; gene_disease; rsid/variants; search; phenotype/{resource}/{phenocode}/markdown; resource_metadata; resources; datasets. Authoritative path list at /api/v1/openapi.json (browseable at /api/v1/docs); confirm against the spec rather than guessing if the user asks about an endpoint not listed here. Use this instead of writing standalone httpx code from scratch when nb01's helpers already cover the call shape, and instead of asking the user to wire auth themselves.
---

# Compose a new marimo notebook in fgx

## What this skill is for

fgx is a single-substrate repo: marimo notebooks against FinnGenie's REST API at `https://finngenie.fi/api/v1/*` via `httpx.get`, and nothing else. There is no Python SDK, no MCP server, no schema cache. The reusable functions you compose with live as `@app.function` helpers inside the notebooks themselves; the API surface is documented inline in nb01's first cell as bare `curl` one-liners.

When the user asks a genetics question that the FinnGenie API can answer, your job is to (a) decide whether to edit nb01 in place or copy-fork to `nb02_*.py`, then (b) compose against the existing helpers rather than re-implementing auth and TSV/JSON parsing from scratch.

## The catalog at this scale

Today fgx has one notebook. The reusable helpers live there:

| File | Reusable functions | What they do |
|---|---|---|
| `notebooks/nb01_pcsk9_walkthrough.py` | `fetch_tsv(path, **params)`, `fetch_json(path, **params)`, `_client()` | Authenticated `httpx.Client` against `https://finngenie.fi/api/v1/*`. `fetch_tsv` returns a `polars.DataFrame`; `fetch_json` returns a list of dicts. Bearer token loaded from `.env` at the repo root. |

Both helpers are pure on the inputs they're given: `path` is the part after `/api/v1`, kwargs become query parameters, response shape follows what the endpoint returns. New notebooks should `from nb01_pcsk9_walkthrough import fetch_tsv, fetch_json` rather than re-deriving the auth + parsing logic.

## Two paths: edit-in-place vs copy-fork

### Path A: edit nb01 in place

Use this when the new question has the **same analysis shape** as nb01's PCSK9 walkthrough, just a different parameter:

- "do the same thing for APOE" -> change `GENE = "PCSK9"` to `GENE = "APOE"`, the rest cascades reactively.
- "use a 200 kb window instead" -> add `window=200000` to the `fetch_tsv("/credible_sets_by_gene/{GENE}", ...)` call.
- "plot the colocalization for a different lead variant" -> change the missense filter or the `head(1)` selector.

Marimo's reactivity reruns dependent cells; don't manually kick anything. This is the cheapest path and preserves a single canonical notebook.

### Path B: copy-fork to `notebooks/nb02_*.py`

Use this when the analysis shape **differs** from nb01:

- Different endpoint chain (e.g. start from a phenotype rather than a gene): `fetch_tsv("/credible_sets_by_phenotype/{resource}/{phenotype}")` -> filter loci -> `fetch_json("/nearest_genes/{variant}", n=3)` to annotate -> chart.
- Different output (manhattan plot, PheWAS layout, multi-resource comparison) -> the cell DAG changes shape, not just parameters.
- Reusable helpers needed by both notebooks (e.g. a small wrapper that fetches `/trait_name_mapping` once and indexes phenotype codes -> names) -> add it to nb01 first as `@app.function`, then `from nb01_pcsk9_walkthrough import that_helper` in nb02.

When forking: copy `nb01_pcsk9_walkthrough.py` -> `nb02_<topic>.py`, edit the PEP 723 header if you need new deps, replace the analysis cells, keep the setup block (loads `.env` and re-exposes the helpers).

## Process for a new composition

1. **Read nb01.** The first cell is the bare-`curl` reference for every endpoint the notebook hits; the second cell defines `_client()`, `fetch_tsv`, `fetch_json`. Skim before composing.
2. **Identify the endpoint.** Map the user's question to a `/api/v1/*` path. The full surface is in the OpenAPI spec at <https://finngenie.fi/api/v1/openapi.json> (browseable Swagger at <https://finngenie.fi/api/v1/docs>) -- 26 paths across 13 tags as of May 2026. Read the spec before guessing; the description's path list is a reminder, not a contract. If you're unsure of the response shape, run `curl -fsS -H "Authorization: Bearer $FINNGENIE_TOKEN" "https://finngenie.fi/api/v1/<path>?format=json"` once to confirm, then write the cell.
3. **Pick the path (A or B).** When in doubt, copy-fork.
4. **Run the kernel.** Use `marimo-pair`'s `execute-code.sh` against the running kernel rather than asking the user to paste output. The setup block is already provisioned from the PEP 723 header; running a single cell returns a `polars.DataFrame` you can summarize inline.
5. **Plot when it adds signal.** Altair is in the PEP 723 deps; a one-cell chart is often more informative than a 20-row table head.
6. **Don't fabricate.** If a helper doesn't exist or an endpoint returns an empty result, say so and ask the user how to proceed -- don't invent a fallback path.

## Gotchas

- **TSV is the API default.** `Content-Type: text/tab-separated-values`; pass `format=json` only when the response is a nested object (e.g. `colocalization_by_variant`). `fetch_tsv` adds nothing; `fetch_json` *unconditionally* injects `format=json` into the query string -- if you pass `format="tsv"` to `fetch_json` it gets silently overridden. Pick the right helper for the response shape you want, don't try to flip format via kwargs.
- **Response sizes can blow up context.** `credible_sets_by_gene/PCSK9` returns 3,224 rows, ~1.3 MB. Filter or summarize in-cell (`pl.DataFrame.filter`, `.head(20)`, `group_by + agg`) before printing the result -- and never let raw API output land in the chat transcript at full size.
- **No data caching.** Every notebook run hits the API live. If a parameter is fixed (e.g. you're iterating on the plotting code, not the gene), bind the dataframe to a variable in an upstream cell so reactive cells downstream don't re-fetch.
- **Lead-variant logic is fragile.** nb01 picks the highest-PIP missense in the GWAS bucket as a stand-in for "lead variant"; that's correct for PCSK9 but not generally. For other genes, look at the credible set explicitly (`get_credible_set_by_id` or filter by `cs_id`) rather than re-using the missense heuristic.
- **Phenotype codes are FinnGen-specific.** `I9_CHD`, `T2D`, `K11_CROHN` etc. come from FinnGen's endpoint definitions, not free text. Hit `fetch_json("/trait_name_mapping")` once to get the full code -> name mapping, then pick the code you actually want before passing strings to phenotype-keyed paths.

## When to revise this skill

If a second notebook lands and starts importing real analysis functions (not just `fetch_tsv`/`fetch_json`) from a sibling, fgx has crossed the threshold where the catalog needs a real composition table -- the way [jx's compose-notebook](https://github.com/broadinstitute/jx/blob/main/.claude/skills/compose-notebook/SKILL.md) carries one. Symptoms that the threshold has been crossed:

- A function defined in nbN is used by nbM via `from nbN_topic import helper`.
- A user asks "is there a notebook that does X?" and the answer requires checking three files.
- Reusable helpers start to outnumber analysis cells in any single notebook.

When any of those land, replace the single-row catalog table above with a per-module table (one row per notebook, listing its `@app.function` helpers and what they do), and add a "Process for a new composition" section that orients around picking-and-gluing across notebooks rather than the edit-vs-fork binary.
