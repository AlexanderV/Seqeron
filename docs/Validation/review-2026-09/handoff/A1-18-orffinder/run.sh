#!/usr/bin/env bash
# One command: build the image, run NCBI ORFfinder 0.4.3 on all test inputs, collect ./out and zip it into ./returned/.
#   ./run.sh                 full run
#   QUICK=1 ./run.sh         reduced smoke test (~700 runs)
#   JOBS=8 ./run.sh          parallel processes inside the container (default: CPU count of the Docker VM)
# Apple Silicon: the NCBI binary is x86-64 Linux; this script always uses --platform linux/amd64 (Docker Desktop emulates it).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-orffinder-oracle:0.4.3
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }

docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out
mkdir -p out returned
docker run --rm --platform linux/amd64 \
  -e JOBS="${JOBS:-}" -e QUICK="${QUICK:-}" -e SINGLE="${SINGLE:-}" -e FORCE="${FORCE:-}" \
  -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" \
  -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log

if command -v zip > /dev/null; then
  rm -f returned/A1-18-orffinder-out.zip
  zip -qr returned/A1-18-orffinder-out.zip out
  ls -l returned/A1-18-orffinder-out.zip
else
  tar czf returned/A1-18-orffinder-out.tar.gz out
  ls -l returned/A1-18-orffinder-out.tar.gz
fi
echo "Send back: docs/Validation/review-2026-09/handoff/A1-18-orffinder/returned/ (the archive above)."
