#!/usr/bin/env bash
# One command: builds the image (the MAFFT 7.505 source built WITH dash_client), runs the battery against the real DASH server through a logging proxy, zips ./out into ./returned/.
# Apple Silicon: linux/amd64 (emulated; about 30 - 90 minutes, mostly waiting for the server).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-dash:7.505
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
docker build --platform linux/amd64 -t "$IMAGE" .
[ "${DASH_RESUME:-}" = 1 ] || rm -rf out
mkdir -p out returned
set +e
docker run --rm --platform linux/amd64 -e DASH_RESUME="${DASH_RESUME:-}" -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" bash /work/in_container.sh 2>&1 | tee -a out/run_console.log
set -e   # the archive is made even when the container reports a build failure
rm -f returned/A2-19c-dash-out.zip returned/A2-19c-dash-out.tar.gz
if command -v zip > /dev/null; then zip -qr returned/A2-19c-dash-out.zip out; else tar czf returned/A2-19c-dash-out.tar.gz out; fi
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/A2-19c-dash/returned/ (the archive above)."
