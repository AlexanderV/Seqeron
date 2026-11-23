#!/usr/bin/env bash
# One command: builds the image (Debian 12 mafft 7.505-1 + fasta34), runs the battery and the sets, zips ./out into ./returned/. Apple Silicon: linux/amd64 (emulated; the run takes a while).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-fasta34-battery:7.505
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
cp ../A2-19c-fastaparttree/inputs/fasta-34.26.5.shar.Z inputs/
docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out; mkdir -p out returned
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
rm -f inputs/fasta-34.26.5.shar.Z returned/A2-19c-fasta34-battery-out.zip returned/A2-19c-fasta34-battery-out.tar.gz
if command -v zip > /dev/null; then zip -qr returned/A2-19c-fasta34-battery-out.zip out; else tar czf returned/A2-19c-fasta34-battery-out.tar.gz out; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/A2-19c-fasta34-battery/returned/ (the archive above)."
