#!/usr/bin/env bash
# One command: build the small image, probe the snapshots (about 2 hours, progress on the console), zip ./out into ./returned/ (the zip is ALWAYS written, also after a failure or Ctrl-C).
# Nothing that is downloaded is compiled or run; the downloaded archives stay in ./out/found on your disk and are NOT in the zip (third-party licences; do not commit them).
set -uo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-a2-12b:s6b-dalilite-snapshots
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
rm -rf out; mkdir -p out returned
pack() {   # also after Ctrl-C (requests.log is flushed request by request; the script writes report.txt on SIGTERM / SIGINT)
  ./pack_out.sh out returned/A2-12b-s6b-out.run3.zip
  ls -l returned
  echo "Send back: docs/Validation/review-2026-09/handoff/A2-12b-s6b-dalilite-snapshots/returned/ (the archive above). Do NOT commit out/found/ (third-party archives)."
}
trap pack EXIT
docker build --platform linux/amd64 -t "$IMAGE" . 2>&1 | tee out/docker_build.log
# --init: tini forwards SIGTERM / SIGINT to in_container.sh, which hands them to the script
docker run --rm --init --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
