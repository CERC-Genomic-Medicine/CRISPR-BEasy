#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 <input.json>" >&2
  exit 1
fi

INPUT="$1"

# iterate each top‐level key
jq -r 'keys[]' "$INPUT" | while read -r ID; do
  # extract values
  PAM=$(jq -r --arg id "$ID" '.[$id].pam'       "$INPUT")
  LENGTH=$(jq -r --arg id "$ID" '.[$id].length'    "$INPUT")
  CFD=$(jq -r --arg id "$ID" '.[$id].CFD_access' "$INPUT")

  # make folder
  mkdir -p "$ID"

  # write the JSON file
  cat > "$ID/$ID.json" <<EOF
{
  "$ID": {
    "pam":       "$PAM",
    "length":    $LENGTH,
    "CFD_access":"$CFD",
    "Pam_side":  "3prime",
    "rs3_access":"$CFD"
  }
}
EOF

  # write cloudgene.yaml
  cat > "$ID/cloudgene.yaml" <<EOF
id: $ID
name: ${LENGTH} bp - ${PAM} - ${ID}
version: 0.0
category: cas_database
properties:
  database_url: \${CLOUDGENE_APP_LOCATION}/$ID.json
EOF

done
