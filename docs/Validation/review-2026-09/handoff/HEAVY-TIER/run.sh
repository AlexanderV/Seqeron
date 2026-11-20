#!/usr/bin/env bash
# One command. Usage:  ./run.sh            (heavy tier: Fuzzing/Properties/Metamorphic of Seqeron.Genomics.Tests)
#                      MODE=ALL ./run.sh   (heavy + full fast tier + all MCP/Integration/SuffixTree.Mcp test projects)
#                      MODE=MCP ./run.sh   (only the MCP projects)   MODE=FAST ./run.sh
#                      WORKERS=4 ./run.sh  (limit parallelism if Docker Desktop has < 12 GB RAM)
# Apple Silicon: the image is native arm64 (no emulation). Docker Desktop -> Settings -> Resources: give it all CPUs and >= 12 GB memory.
set -uo pipefail
cd "$(dirname "$0")"
REPO="$(cd ../../../../.. && pwd)"
IMAGE=seqeron-heavy-tier:1
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
docker build -t "$IMAGE" . 2>&1 | tee build_image.log
rm -rf out; mkdir -p out returned
docker run --rm -e MODE="${MODE:-HEAVY}" -e WORKERS="${WORKERS:-}" -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" \
  -v "$REPO":/src:ro -v "$PWD/out":/out "$IMAGE" in_container.sh 2>&1 | tee out/run_console.log
rm -f returned/heavy-tier-out.zip
if command -v zip > /dev/null; then zip -qr returned/heavy-tier-out.zip out build_image.log; else tar czf returned/heavy-tier-out.tar.gz out build_image.log; fi
ls -l returned; echo; cat out/summary.txt 2>/dev/null
echo "Send back: docs/Validation/review-2026-09/handoff/HEAVY-TIER/returned/ (the archive above)."
