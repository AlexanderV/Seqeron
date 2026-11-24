#!/usr/bin/env bash
# One command: builds the image (MAFFT 7.505 with extensions + ViennaRNA), runs the sets (each twice), zips ./out into ./returned/A2-19c-rna-out.zip. Apple Silicon: linux/amd64 (emulated).
#   RNA_QUICK=1 ./run.sh      a check of the package with 6 sets
# The docker run is wrapped in `set +e`: whatever happens inside, out/ is zipped into returned/A2-19c-rna-out.zip.
set -uo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-rna:7.505
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
python3 make_inputs.py 2> /dev/null || true
set +e
docker build --platform linux/amd64 -t "$IMAGE" . 2>&1 | tee build.log
rm -rf out; mkdir -p out returned
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -e RNA_QUICK="${RNA_QUICK:-0}" -e RUN_TIMEOUT="${RUN_TIMEOUT:-300}" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
cp build.log out/ 2>/dev/null
set -e
rm -f returned/A2-19c-rna-out.zip returned/A2-19c-rna-out.tar.gz
if command -v zip > /dev/null; then zip -qr returned/A2-19c-rna-out.zip out; else tar czf returned/A2-19c-rna-out.tar.gz out; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/A2-19c-rna/returned/ (the archive above)."
