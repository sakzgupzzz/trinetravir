#!/usr/bin/env bash
# Stream-extract infected_monocytes_cluster_singlets.rds from Zenodo 10.5281/zenodo.4273999.
# 36.6GB tar.gz; selective extraction via bsdtar --include.
set -euo pipefail

URL="https://zenodo.org/api/records/4273999/files/inputs.tar.gz/content"
DEST="/Users/sakshamgupta/Documents/coding_projects/trinetravir/data/raw/randolph_2021/zenodo_inputs"
LOG="/Users/sakshamgupta/Documents/coding_projects/trinetravir/data/raw/randolph_2021/zenodo_extract.log"

mkdir -p "$DEST"
cd "$DEST"

echo "$(date) start streaming extract from $URL" | tee -a "$LOG"
echo "$(date) dest=$DEST" | tee -a "$LOG"

# Stream-extract: curl pipes to bsdtar, filter matching files only.
# Patterns:
#   - infected_monocytes_cluster_singlets.rds (primary target)
#   - any other infected_*.rds (informative for v1.5)
#   - monocyte demux / metadata files (cross-check)
curl --location --silent --show-error --retry 3 --retry-delay 30 \
    --connect-timeout 60 --max-time 21600 \
    "$URL" \
  | tar -xzvf - \
      --include '*infected_monocytes_cluster_singlets*' \
      --include '*infected_*.rds' \
      2>&1 | tee -a "$LOG"

echo "$(date) extract complete; listing:" | tee -a "$LOG"
find "$DEST" -type f -name '*.rds' -exec ls -lh {} \; | tee -a "$LOG"
