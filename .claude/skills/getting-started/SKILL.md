---
name: getting-started
description: >-
  Walk a first-time fgx user from a fresh clone to a running marimo kernel
  with agent composition enabled. Trigger when the user says "help me get
  started", "onboard me", "set me up", "I'm new to fgx", "first time using
  fgx", "what do I do next", or asks any FinnGenie / human-genetics
  composition question before a marimo kernel is running and marimo-pair
  is connected. Verifies the FINNGENIE_TOKEN is set, sets up uv, prompts
  the user to install the marimo-pair skill, launches the marimo server
  on nb01.
---

# Getting started with fgx

Your job: get this user from a cold clone to a live marimo kernel pointed at
nb01, then hand off to interactive composition. fgx has only one notebook
today (`nb01_pcsk9_walkthrough.py`), so there is no `compose-notebook` skill
yet -- the user composes new analyses by editing nb01 or copy-forking it.

## Setup flow

### 1. Verify FINNGENIE_TOKEN is set

Check `.env` exists and contains a non-empty `FINNGENIE_TOKEN`:

    grep -q "^FINNGENIE_TOKEN=." .env 2>/dev/null && echo OK || echo MISSING

If `MISSING`, walk the user through:

1. Open <https://finngenie.broadinstitute.org/> -> `MCP/API KEYS` -> `Create key`.
2. `cp .env.example .env`.
3. Paste the key after `FINNGENIE_TOKEN=` in `.env`.

Confirm the same token works for both the FinnGenie MCP and the REST API
fgx uses; one key, both paths.

### 2. Verify uv is installed

Run `uv --version`. If it fails, tell the user to run:

    curl -LsSf https://astral.sh/uv/install.sh | sh

Then have them source their shell profile (`. ~/.zshrc`) or open a new
terminal. Re-check `uv --version`.

### 3. Install the marimo-pair skill

marimo-pair is a self-contained skill (SKILL.md + bundled `scripts/`)
distributed via [skills.sh](https://skills.sh). Install it globally for
the user:

    npx skills add marimo-team/marimo-pair -g --agent claude-code -y

After install, the user should restart their Claude Code session so the
skill loads. On the next session, marimo-pair's tools (e.g. running
`scripts/execute-code.sh`) are available via its `allowed-tools`
frontmatter. If you don't see the skill, it isn't loaded yet.

### 4. Launch the marimo server

From the fgx repo root, pick a free port and start nb01 in `--sandbox`
mode so the PEP 723 header provisions a venv automatically:

    PORT=$(python -c "import socket; s=socket.socket(); s.bind(('127.0.0.1',0)); print(s.getsockname()[1])")
    env -u PYTHONPATH uvx marimo edit --sandbox --headless --no-token \
        --port $PORT notebooks/nb01_pcsk9_walkthrough.py

Use `run_in_background=true` on the Bash call so you can poll while
deps install. First launch installs marimo + polars + httpx + altair +
python-dotenv (~30 seconds). Verify the server is up with:

    curl -s -o /dev/null -w "%{http_code}" http://127.0.0.1:$PORT/

Expect HTTP 200. Tell the user the URL so they can open the browser UI
alongside your session if they want to watch cells render.

### 5. Hand off

Once the kernel is live and marimo-pair is connected:

1. Ask the user what genetics question they want to explore (a gene, a
   phenotype code like `I9_CHD`, a variant in `chr:pos:ref:alt` form).
2. Read nb01's first cell with the user -- it shows the bare-`curl`
   equivalent of every API call the notebook makes. The reusable
   functions (`fetch_tsv`, `fetch_json`) are defined in the setup block
   and can be called from any new cell.
3. For their question, decide whether to edit nb01 in place (small
   tweak, same gene/phenotype shape) or copy it to `notebooks/nb02_*.py`
   (different analysis pattern). When in doubt, copy-fork.
4. Use marimo-pair's `execute-code.sh` to run cells against the live
   kernel rather than asking the user to paste output.

## Gotchas

- **Nix shells poison PYTHONPATH** with a bad websockets shim that
  crashes `marimo edit` on startup. The `env -u PYTHONPATH` prefix in
  step 4 avoids this. Apply it to any marimo invocation on Nix.
- **`--sandbox` is required.** Without it, `uvx marimo edit` opens the
  file picker but does not provision a venv from the PEP 723 header --
  opening any notebook then fails with `ModuleNotFoundError: polars`.
- **Don't share the `FINNGENIE_TOKEN`.** It's personal and gates access
  to embargoed FinnGen data. Never commit `.env`; never paste the token
  into chat. The portal's "Created" timestamp lets the user audit which
  keys are active.
- **No data caching.** Unlike jx, fgx hits the API live every cell.
  Every notebook run re-fetches; that's intentional (no stale cache to
  manage). If the API is slow today, suggest the user pin a
  parameter (gene/phenotype) so reactive cells re-fetch less.

## When the catalog grows

If a second notebook appears that imports functions from nb01 (e.g.
`from nb01_pcsk9_walkthrough import fetch_tsv`), fgx has crossed the
threshold where a `compose-notebook` skill earns its keep. The pattern
to copy is at
<https://github.com/broadinstitute/jx/blob/main/.claude/skills/compose-notebook/SKILL.md>.
Until then, this skill is the whole onboarding path.
