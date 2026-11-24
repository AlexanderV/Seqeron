#!/usr/bin/env bash
# One command: builds the image (Debian 12 mafft 7.505-1 + fasta34), runs the 180 sets, zips ./out into ./returned/. Apple Silicon: linux/amd64 (emulated; the run takes a while, about 20 - 40 minutes).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-fastapair2:7.505
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
cp ../A2-19c-fastaparttree/inputs/fasta-34.26.5.shar.Z inputs/
docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out; mkdir -p out returned
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
rm -f inputs/fasta-34.26.5.shar.Z returned/A2-19c-fastapair2-out.zip returned/A2-19c-fastapair2-out.tar.gz
if command -v zip > /dev/null; then zip -qr returned/A2-19c-fastapair2-out.zip out; else tar czf returned/A2-19c-fastapair2-out.tar.gz out; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/A2-19c-fastapair2/returned/ (the archive above)."
