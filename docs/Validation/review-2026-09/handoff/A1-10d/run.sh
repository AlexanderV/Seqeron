#!/usr/bin/env bash
# One command: build the image (compiles ClustalW 2.1 with the 3-line -ENDGAPS patch), run every set in four modes, zip ./out into ./returned/.
# Apple Silicon: the image is linux/amd64 (Docker Desktop emulates it).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-clustalw-endgaps:2.1
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out; mkdir -p out returned
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
rm -f returned/A1-10d-out.zip
if command -v zip > /dev/null; then zip -qr returned/A1-10d-out.zip out; else tar czf returned/A1-10d-out.tar.gz out; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/A1-10d/returned/ (the archive above)."
