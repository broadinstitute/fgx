ggsql := env_var_or_default("GGSQL", "ggsql")
data  := "queries/data"
db    := "queries/data/cache.duckdb"

default:
    @just --list

# Fetch all TSVs needed by the queries/ layer into queries/data/. Auth boundary lives here.
fetch:
    @[ -f .env ] || (echo "Missing .env. Copy .env.example and paste your token." && exit 1)
    @set -a; . ./.env; set +a; \
        if [ -z "$FINNGENIE_TOKEN" ]; then echo "FINNGENIE_TOKEN not set in .env" >&2; exit 1; fi; \
        mkdir -p {{data}}; \
        echo ">>> fetching credible sets near PCSK9"; \
        curl -fsS -H "Authorization: Bearer $FINNGENIE_TOKEN" \
             "https://finngenie.fi/api/v1/credible_sets_by_gene/PCSK9?format=tsv" \
             -o {{data}}/pcsk9_cs.tsv; \
        echo ">>> fetching credible sets for I9_CHD (FinnGen R13)"; \
        curl -fsS -H "Authorization: Bearer $FINNGENIE_TOKEN" \
             "https://finngenie.fi/api/v1/credible_sets_by_phenotype/finngen/I9_CHD?format=tsv" \
             -o {{data}}/chd_cs.tsv; \
        wc -l {{data}}/*.tsv

# Load TSVs into a local DuckDB so ggsql files can run without re-fetching.
build-db: fetch
    @rm -f {{db}}
    duckdb {{db}} -c "\
        CREATE TABLE pcsk9_cs AS SELECT * FROM read_csv_auto('{{data}}/pcsk9_cs.tsv', header=true, sep='\t'); \
        CREATE TABLE chd_cs   AS SELECT * FROM read_csv_auto('{{data}}/chd_cs.tsv',   header=true, sep='\t'); \
    "
    @duckdb {{db}} -c "SELECT 'pcsk9_cs' AS t, COUNT(*) FROM pcsk9_cs UNION ALL SELECT 'chd_cs', COUNT(*) FROM chd_cs;"

# Render every q*.gsql in queries/ to a Vega-Lite JSON (paste into vega editor for SVG).
render: build-db
    @mkdir -p queries/rendered
    @for f in queries/q*.gsql; do \
        name=$(basename "$f" .gsql); \
        echo ">>> rendering $name"; \
        {{ggsql}} run "$f" --reader "duckdb://{{db}}" --output "queries/rendered/$name.json" || exit 1; \
    done
    @ls queries/rendered/*.json

# Open nb01 in a live marimo kernel.
notebook:
    uv run --with marimo --with httpx --with polars --with altair --with python-dotenv \
        marimo edit notebooks/nb01_pcsk9_walkthrough.py

# Trash fetched data, cached DB, and rendered specs.
clean:
    trash {{data}}/*.tsv {{db}} queries/rendered/*.json 2>/dev/null || true
