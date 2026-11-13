#!/usr/bin/env bash
# Ship case inputs to afml-r6:~/pancake/cases.
#
# Sends only files git tracks under cases/ (configs, scripts, reference data,
# READMEs) whose (path, size) differ remotely. Run output is never shipped and
# nothing is ever pulled: the remote is a place to compute, not to store.
# Retrieve finished runs with cases/retrieve_from_r6.py, which verifies them
# into D:\Seongho\Sim\pancake and then removes them from the remote.
#
# Usage: bash cases/sync_to_r6.sh [--dry-run]
set -euo pipefail
HOST=user@143.248.99.6
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

git ls-files -z cases | xargs -0 stat -c '%s %n' | sort -k2 > "$WORK/local.txt"
ssh "$HOST" "mkdir -p ~/pancake/cases && cd ~/pancake && find cases -type f -printf '%s %p\n'" \
    | sort -k2 > "$WORK/remote.txt"
awk 'NR==FNR{r[$2]=$1;next} {if(!($2 in r) || r[$2]!=$1) print $2}' \
    "$WORK/remote.txt" "$WORK/local.txt" > "$WORK/ship.txt"

N=$(wc -l < "$WORK/ship.txt")
echo "[sync] $N tracked case input files to ship"
if [ "${1:-}" = "--dry-run" ] || [ "$N" -eq 0 ]; then
    cat "$WORK/ship.txt"; exit 0
fi
tar -cf - -T "$WORK/ship.txt" | ssh "$HOST" "tar -C ~/pancake -xf -"
echo "[sync] done $(date '+%F %T')"
