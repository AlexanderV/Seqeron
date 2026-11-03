#!/usr/bin/env bash
# Sourced by every run_*.sh.  Variables you may override from the environment:
#   BENCH_DIR  directory that contains data/ queries/ MANIFEST.tsv lib/   (default: parent of this folder)
#   WORK_DIR   scratch for clones/builds                                  (default: competitors/_work)
#   RAW_FILE   raw per-repeat rows are APPENDED here                      (default: $BENCH_DIR/runs/raw_competitors.tsv)
#   REPEATS    repeats per measurement                                    (default: 3)
#   INPUTS_DNA inputs run for DNA tools                                   (default: "dna_1m dna_10m rep_10m dna_20m")
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH="${BENCH_DIR:-$(cd "$HERE/.." && pwd)}"
WORK="${WORK_DIR:-$HERE/_work}"
RAW="${RAW_FILE:-$BENCH/runs/raw_competitors.tsv}"
REPEATS="${REPEATS:-3}"
INPUTS_DNA="${INPUTS_DNA:-dna_1m dna_10m rep_10m dna_20m}"
BL="python3 $BENCH/lib/benchlib.py"
mkdir -p "$WORK" "$(dirname "$RAW")"
cd "$BENCH"

meas() { $BL measure --raw "$RAW" --data-dir "$BENCH" --repeats "$REPEATS" "$@" || echo "[BLOCKED/FAILED] meas $* -- continuing" >&2; }
need() { command -v "$1" >/dev/null 2>&1 || { echo "missing required tool: $1" >&2; exit 1; }; }
# clone_at <url> <sha> <dir>   (shallow-fetch of one commit would need server support; plain clone + checkout is the portable way)
clone_at() {
  [ -d "$3/.git" ] || git clone --quiet "$1" "$3"
  git -C "$3" checkout --quiet "$2"
  echo "$1 @ $(git -C "$3" rev-parse HEAD)" >> "$WORK/SOURCES.txt"
}
# record the licence file text head + sha for provenance
note_licence() { # dir
  for f in LICENSE LICENSE.md LICENSE.txt COPYING COPYING.md; do
    [ -f "$1/$f" ] && { echo "licence file $1/$f: $(head -c 120 "$1/$f" | tr '\n' ' ')" >> "$WORK/SOURCES.txt"; return; }
  done
}
cpu_model() { grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //'; }
