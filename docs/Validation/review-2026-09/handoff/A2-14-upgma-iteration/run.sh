#!/usr/bin/env bash
# One command: build the image (compiles ClustalW 2.1 from the Debian source at -O0), run every set of inputs/sets.tsv, zip ./out into ./returned/. Apple Silicon: the image is linux/amd64 (emulated).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-clustalw-upgma-iteration:2.1
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out; mkdir -p out returned
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" bash /work/in_container.sh 2>&1 | tee out/run_console.log
rm -f returned/A2-14-upgma-iteration-out.zip returned/A2-14-upgma-iteration-out.tar.gz
if command -v zip > /dev/null; then zip -qr returned/A2-14-upgma-iteration-out.zip out; else tar czf returned/A2-14-upgma-iteration-out.tar.gz out; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/A2-14-upgma-iteration/returned/ (the archive above)."
