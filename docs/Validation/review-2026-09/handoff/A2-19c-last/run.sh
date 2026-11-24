#!/usr/bin/env bash
# One command: builds the image (Debian 12 mafft 7.505-1 + LAST candidates), runs the sets (each twice), zips ./out into ./returned/. Apple Silicon: linux/amd64 (emulated; the build compiles LAST, the run takes about 30 - 90 minutes).
#   LAST_QUICK=1 ./run.sh      a check of the package with 6 sets (about 10 minutes)
# The docker run is wrapped in `set +e`: whatever happens inside, out/ is zipped into returned/A2-19c-last-out.zip.
set -uo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-last:7.505
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
python3 make_inputs.py 2> /dev/null || true
set +e
docker build --platform linux/amd64 -t "$IMAGE" . 2>&1 | tee build.log
rm -rf out; mkdir -p out returned
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -e LAST_QUICK="${LAST_QUICK:-0}" -e RUN_TIMEOUT="${RUN_TIMEOUT:-300}" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
cp build.log out/ 2>/dev/null
set -e
rm -f returned/A2-19c-last-out.zip returned/A2-19c-last-out.tar.gz
if command -v zip > /dev/null; then zip -qr returned/A2-19c-last-out.zip out; else tar czf returned/A2-19c-last-out.tar.gz out; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/A2-19c-last/returned/ (the archive above)."
