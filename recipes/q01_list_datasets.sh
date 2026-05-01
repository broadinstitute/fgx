#!/usr/bin/env bash
# q01: list every dataset FinnGenie exposes, plus its products and stats.
# Endpoint: GET /api/v1/datasets
#
# Use this first to confirm the token works and to discover dataset_ids
# (resource + version pairs) you can pass to other recipes.

set -euo pipefail

cd "$(dirname "$0")/.."
[ -f .env ] || { echo "Missing .env. Copy .env.example and paste your token." >&2; exit 1; }
set -a; . ./.env; set +a
: "${FINNGENIE_TOKEN:?FINNGENIE_TOKEN not set in .env}"

URL="https://finngenie.fi/api/v1/datasets"
echo ">>> GET $URL" >&2

curl -fsS -H "Authorization: Bearer $FINNGENIE_TOKEN" "$URL" \
  | jq -r '
      "dataset_id\tresource\tversion\tdata_type\tn_phenotypes\tproducts",
      (.[] | [
        .dataset_id,
        .resource,
        .version,
        .data_type,
        (.stats.n_phenotypes // .n_phenotypes // "-"),
        ((.products // {}) | keys | join(","))
      ] | @tsv)
  ' \
  | column -t -s $'\t'
