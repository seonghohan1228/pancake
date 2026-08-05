#!/usr/bin/env bash
# Mirror cases/ into the remote repo clone: afml-r6:~/pancake/cases.
# Manifest-based: ships exactly the files whose (path, size) differ or are
# absent remotely, so it is resumable and idempotent. Append-only archive
# (results are never deleted), so no deletion propagation by design.
# Usage: bash cases/sync_to_r6.sh
set -euo pipefail
HOST=user@143.248.99.6
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "[sync] building manifests"
ssh "$HOST" "mkdir -p ~/pancake/cases && cd ~/pancake && find cases -type f -printf '%s %p\n' | sort -k2" > /tmp/r6_manifest.txt
find cases -type f -printf '%s %p\n' | sort -k2 > /tmp/local_manifest.txt

# files present locally whose size+path pair is not matched remotely
join -v1 -j2 -o 1.2 /tmp/local_manifest.txt /tmp/r6_manifest.txt > /tmp/ship_pathonly.txt || true
awk 'NR==FNR{r[$2]=$1;next} {if(!($2 in r) || r[$2]!=$1) print $2}' \
    /tmp/r6_manifest.txt /tmp/local_manifest.txt > /tmp/ship_list.txt

N=$(wc -l < /tmp/ship_list.txt)
echo "[sync] $N files to ship"
if [ "$N" -gt 0 ]; then
    tar -cf - -T /tmp/ship_list.txt | ssh "$HOST" "tar -C ~/pancake -xf -"
fi

# verify + report remote-only files (candidates to pull back)
ssh "$HOST" "cd ~/pancake && find cases -type f -printf '%s %p\n' | sort -k2" > /tmp/r6_manifest2.txt
awk 'NR==FNR{l[$2]=$1;next} {if(!($2 in l)) print $2}' \
    /tmp/local_manifest.txt /tmp/r6_manifest2.txt > /tmp/remote_only.txt
MISSING=$(awk 'NR==FNR{r[$2]=$1;next} {if(!($2 in r) || r[$2]!=$1) print $2}' \
    /tmp/r6_manifest2.txt /tmp/local_manifest.txt | wc -l)
echo "[sync] after ship: $MISSING local files still unmatched remotely (want 0)"
echo "[sync] remote-only files (exist on afml-r6, not local): $(wc -l < /tmp/remote_only.txt) -> /tmp/remote_only.txt"
echo "[sync] done $(date '+%F %T')"
