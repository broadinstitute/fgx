#!/usr/bin/env bash
# q03: every gene/locus driving a phenotype (genome-wide fine-mapping summary).
# Endpoint: GET /api/v1/credible_sets_by_phenotype/{resource}/{phenotype}
#
# Worked example: I9_CHD (coronary heart disease) in FinnGen R13.
# Use list_datasets (q01) to see available resources; phenotype codes
# come from FinnGen endpoint definitions.

set -euo pipefail

RESOURCE="${RESOURCE:-finngen}"
PHENOTYPE="${1:-I9_CHD}"

cd "$(dirname "$0")/.."
[ -f .env ] || { echo "Missing .env. Copy .env.example and paste your token." >&2; exit 1; }
set -a; . ./.env; set +a
: "${FINNGENIE_TOKEN:?FINNGENIE_TOKEN not set in .env}"

URL="https://finngenie.fi/api/v1/credible_sets_by_phenotype/${RESOURCE}/${PHENOTYPE}"
echo ">>> GET ${URL} (TSV)" >&2

TMP=$(mktemp -t "fgx-q03.${PHENOTYPE}.XXXXXX.tsv")
trap 'rm -f "$TMP"' EXIT
curl -fsS -H "Authorization: Bearer $FINNGENIE_TOKEN" "$URL?format=tsv" -o "$TMP"

echo ">>> rows fetched: $(($(wc -l < "$TMP") - 1))" >&2

# Top 20 loci by significance, with the nearest-gene annotation that ships
# in the response. PIP is the per-variant causal probability within its credible set.
duckdb -c "
  SELECT chr, pos, ref, alt,
         mlog10p, beta, pip,
         most_severe, gene_most_severe,
         cs_id, cs_size
  FROM read_csv_auto('$TMP', header=true, sep='\t')
  ORDER BY mlog10p DESC NULLS LAST
  LIMIT 20;
"
