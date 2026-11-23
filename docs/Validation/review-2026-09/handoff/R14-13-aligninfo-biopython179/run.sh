#!/usr/bin/env bash
# One command: builds the image (Python 3.9 + Biopython 1.79), writes cases.jsonl, runs every call of every case, zips ./out into ./returned/.
# Apple Silicon: the image is linux/amd64 (Docker Desktop emulates it; a few minutes at most).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-r14-13-aligninfo:bp179
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out; mkdir -p out returned
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" bash /work/in_container.sh 2>&1 | tee out/run_console.log
rm -f returned/R14-13-aligninfo-biopython179-out.zip returned/R14-13-aligninfo-biopython179-out.tar.gz
if command -v zip > /dev/null; then zip -qr returned/R14-13-aligninfo-biopython179-out.zip out; else tar czf returned/R14-13-aligninfo-biopython179-out.tar.gz out; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/R14-13-aligninfo-biopython179/returned/ (the archive above)."
