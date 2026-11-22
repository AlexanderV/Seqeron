#!/usr/bin/env bash
# One command: builds the image (EMBOSS 6.6.0 + the real PROSITE database), runs every row of inputs/rows.txt, zips ./out into ./returned/.
# Apple Silicon: the image is linux/amd64 (Docker Desktop emulates it; slower).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-wp62-patmatmotifs:prosite
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out; mkdir -p out returned
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
rm -f returned/WP62-patmatmotifs-prosite-out.zip returned/WP62-patmatmotifs-prosite-out.tar.gz
if command -v zip > /dev/null; then zip -qr returned/WP62-patmatmotifs-prosite-out.zip out inputs; else tar czf returned/WP62-patmatmotifs-prosite-out.tar.gz out inputs; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/WP62-patmatmotifs-prosite/returned/ (the archive above)."
