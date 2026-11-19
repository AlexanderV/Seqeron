#!/usr/bin/env bash
# One command: build the image, run the structure sets through the REAL SAP / DaliLite and t_coffee, zip ./out into ./returned/.
# Apple Silicon: the image is linux/amd64 (Docker Desktop emulates it).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-a2-12:sap-dali
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
# the generator of the committed oracle is shared with the repository
cp ../../../../../scripts/alignment/generate_tcoffee_external_oracle.py ./generate_tcoffee_external_oracle.py
docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out; mkdir -p out returned
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
rm -f returned/A2-12-out.zip returned/A2-12-out.tar.gz
if command -v zip > /dev/null; then zip -qr returned/A2-12-out.zip out; else tar czf returned/A2-12-out.tar.gz out; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/A2-12-sap-dali/returned/ (the archive above)."
