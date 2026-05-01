#!/usr/bin/env bash
# q04: traits that share a causal signal at a variant -- the standard mechanism question.
# Endpoint: GET /api/v1/colocalization_by_variant/{variant}
#
# Worked example: 1:55039974:G:T (PCSK9 missense, the LDL-lowering variant).
# Expected hit: PCSK9-mediated colocalization between LDL/nonHDL labs and
# familial hypercholesterolemia / cardiovascular disease endpoints.

set -euo pipefail

VARIANT="${1:-1:55039974:G:T}"

cd "$(dirname "$0")/.."
[ -f .env ] || { echo "Missing .env. Copy .env.example and paste your token." >&2; exit 1; }
set -a; . ./.env; set +a
: "${FINNGENIE_TOKEN:?FINNGENIE_TOKEN not set in .env}"

URL="https://finngenie.fi/api/v1/colocalization_by_variant/${VARIANT}?format=json"
echo ">>> GET ${URL}" >&2

# Group colocalized trait pairs by trait1 and show the strongest partners.
# This recipe uses jq because the colocalization endpoint defaults to JSON
# (the response is already nested as trait pairs).
curl -fsS -H "Authorization: Bearer $FINNGENIE_TOKEN" "$URL" \
  | jq -r '
      "trait1\ttrait2\tdataset1\tdataset2\thit1_mlog10p\thit2_mlog10p\thit1_beta\thit2_beta",
      (.[] | [
        .trait1, .trait2,
        .dataset1, .dataset2,
        (.hit1_mlog10p // "-"), (.hit2_mlog10p // "-"),
        (.hit1_beta // "-"),    (.hit2_beta // "-")
      ] | @tsv)
  ' \
  | (head -1; tail -n +2 | sort -t $'\t' -k5,5gr) \
  | head -25 \
  | column -t -s $'\t'
