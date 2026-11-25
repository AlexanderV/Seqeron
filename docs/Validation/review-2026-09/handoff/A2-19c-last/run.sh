#!/usr/bin/env bash
# One command: builds the image (Debian 12 mafft 7.505-1 + LAST candidates), runs the sets, zips ./out into ./returned/. Apple Silicon: linux/amd64 (emulated; the build compiles LAST, the run takes a while).
#   LAST_QUICK=1 ./run.sh        a check of the package with 6 sets
#   LAST_ONLY="set200 set201" ./run.sh   only these sets (the 2-sequence sets are the cheapest: set184 ... set203)
#   RUN_TIMEOUT=3600 (seconds per mafft run), TWICE_MAX=900 (a set whose first run took longer is run once only), SMOKE_FAST=90 (seconds a LAST candidate may need for the real lastal call of mafft)
# Whatever happens (an error, Ctrl-C), out/ is zipped into returned/A2-19c-last-out.zip.
set -uo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-last:7.505
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
pack() {   # the next session reads whatever is in the zip (out/index.tsv and the sets are written set by set)
  cp build.log out/ 2>/dev/null
  rm -f returned/A2-19c-last-out.zip returned/A2-19c-last-out.tar.gz
  if command -v zip > /dev/null; then zip -qr returned/A2-19c-last-out.zip out; else tar czf returned/A2-19c-last-out.tar.gz out; fi
  ls -l returned
  echo "Send back: docs/Validation/review-2026-09/handoff/A2-19c-last/returned/ (the archive above)."
}
python3 make_inputs.py 2> /dev/null || true
rm -rf out; mkdir -p out returned
trap pack EXIT
docker build --platform linux/amd64 --build-arg "SMOKE_FAST=${SMOKE_FAST:-90}" -t "$IMAGE" . 2>&1 | tee build.log
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -e LAST_QUICK="${LAST_QUICK:-0}" -e RUN_TIMEOUT="${RUN_TIMEOUT:-3600}" -e TWICE_MAX="${TWICE_MAX:-900}" -e LAST_ONLY="${LAST_ONLY:-}" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
