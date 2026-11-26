#!/usr/bin/env bash
# One command: builds the image (Debian 12 mafft 7.505-1 + fasta34), runs the sets 254 - 288 (RUN 4; sets 1 - 253 are integrated and are NOT repeated), zips ./out into ./returned/. Apple Silicon: linux/amd64 (emulated; the run takes a while, about 20 - 60 minutes).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-fastapair2:7.505
NAME=seqeron-fastapair2-run4
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
cp ../A2-19c-fastaparttree/inputs/fasta-34.26.5.shar.Z inputs/
docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out; mkdir -p out returned
# RUN 4: the complete `mafft` and the direct dndfast7 run for the sets SETS_FROM .. SETS_TO (default 254 .. 288, the ones that runs 2 and 3 did not reach: run 3 stopped at set 254 with an error of the script, now fixed:
# non-FASTA output such as --clustalout is stored raw, and every set runs in its own try/except). Every program run has a time limit (SET_TIMEOUT, 1200 s).
# The zip is written also when the run is interrupted (Ctrl-C): the container is stopped, then ./out is zipped - the next step of the session reads whatever is in it.
pack() {
  trap - EXIT INT TERM
  docker stop -t 30 "$NAME" > /dev/null 2>&1 || true
  rm -f inputs/fasta-34.26.5.shar.Z
  mkdir -p returned
  if command -v zip > /dev/null; then rm -f returned/A2-19c-fastapair2-run4-out.zip; zip -qr returned/A2-19c-fastapair2-run4-out.zip out; else rm -f returned/A2-19c-fastapair2-run4-out.tar.gz; tar czf returned/A2-19c-fastapair2-run4-out.tar.gz out; fi
  ls -l returned
}
trap pack EXIT
trap 'exit 130' INT TERM
docker rm -f "$NAME" > /dev/null 2>&1 || true
docker run --rm --init --name "$NAME" --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -e SETS_FROM="${SETS_FROM:-254}" -e SETS_TO="${SETS_TO:-288}" -e DIRECT_FROM="${DIRECT_FROM:-${SETS_FROM:-254}}" -e SET_TIMEOUT="${SET_TIMEOUT:-1200}" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
echo "Send back: docs/Validation/review-2026-09/handoff/A2-19c-fastapair2/returned/A2-19c-fastapair2-run4-out.zip (next to the archives of runs 2 and 3)."
