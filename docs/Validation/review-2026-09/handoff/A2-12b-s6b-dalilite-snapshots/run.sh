#!/usr/bin/env bash
# One command: build the small image, probe the snapshots, zip ./out into ./returned/ (the zip is ALWAYS written, also after a failure). Nothing that is downloaded is compiled or run.
set -uo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-a2-12b:s6b-dalilite-snapshots
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
rm -rf out; mkdir -p out returned
docker build --platform linux/amd64 -t "$IMAGE" . 2>&1 | tee out/docker_build.log
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
rm -f returned/A2-12b-s6b-out.zip returned/A2-12b-s6b-out.tar.gz
if command -v zip > /dev/null; then zip -qr returned/A2-12b-s6b-out.zip out; else tar czf returned/A2-12b-s6b-out.tar.gz out; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/A2-12b-s6b-dalilite-snapshots/returned/ (the archive above)."
