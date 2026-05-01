#!/usr/bin/env bash
# q02: credible-set fine-mapped associations near a gene (GWAS, eQTL, pQTL, sQTL, caQTL).
# Endpoint: GET /api/v1/credible_sets_by_gene/{gene}
#
# Worked example: PCSK9. Returns the canonical loss-of-function story --
# the missense variant 1:55039974:G:T (PIP=1.0, mlog10p~512) underlying
# evolocumab/alirocumab. Edit GENE below for your own query.

set -euo pipefail

GENE="${1:-PCSK9}"
WINDOW="${WINDOW:-100000}"   # bp around the gene; default matches API default

cd "$(dirname "$0")/.."
[ -f .env ] || { echo "Missing .env. Copy .env.example and paste your token." >&2; exit 1; }
set -a; . ./.env; set +a
: "${FINNGENIE_TOKEN:?FINNGENIE_TOKEN not set in .env}"

URL="https://finngenie.fi/api/v1/credible_sets_by_gene/${GENE}?window=${WINDOW}"
echo ">>> GET ${URL} (TSV)" >&2

# TSV streams straight into duckdb. Show the top 20 GWAS associations by mlog10p.
TMP=$(mktemp -t "fgx-q02.${GENE}.XXXXXX.tsv")
trap 'rm -f "$TMP"' EXIT
curl -fsS -H "Authorization: Bearer $FINNGENIE_TOKEN" "$URL" -o "$TMP"

echo ">>> rows fetched: $(($(wc -l < "$TMP") - 1))" >&2

duckdb -c "
  SELECT data_type, trait, resource, dataset, chr, pos, ref, alt,
         mlog10p, beta, pip, most_severe, gene_most_severe
  FROM read_csv_auto('$TMP', header=true, sep='\t')
  WHERE data_type = 'GWAS'
  ORDER BY mlog10p DESC NULLS LAST
  LIMIT 20;
"
