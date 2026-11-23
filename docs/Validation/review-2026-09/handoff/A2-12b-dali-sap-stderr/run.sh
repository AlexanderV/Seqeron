#!/usr/bin/env bash
# One command: build the image, run the sets, zip ./out into ./returned/. Apple Silicon: the image is linux/amd64 (Docker Desktop emulates it; SAP / DaliLite builds take a few minutes).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-a2-12b:dali-sap-stderr
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
cp ../../../../../scripts/alignment/generate_tcoffee_external_oracle.py ./generate_tcoffee_external_oracle.py
docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out; mkdir -p out returned
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
rm -f returned/A2-12b-out.run7.zip returned/A2-12b-out.run7.tar.gz
if command -v zip > /dev/null; then zip -qr returned/A2-12b-out.run7.zip out; else tar czf returned/A2-12b-out.run7.tar.gz out; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/A2-12b-dali-sap-stderr/returned/ (the archive above)."
