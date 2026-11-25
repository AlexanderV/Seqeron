#!/usr/bin/env bash
# One command: builds the image (Debian 12 mafft 7.505-1 + fasta34), runs the 180 sets, zips ./out into ./returned/. Apple Silicon: linux/amd64 (emulated; the run takes a while, about 20 - 40 minutes).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-fastapair2:7.505
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
cp ../A2-19c-fastaparttree/inputs/fasta-34.26.5.shar.Z inputs/
docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out; mkdir -p out returned
# RUN 3 (run 2 was stopped after set 253): all 288 sets run `mafft` again (results.tsv is now written row by row, every program run has a time limit); the direct dndfast7 runs start at set 254 (those of 1 - 253 are in run 2).
# The zip is written also when the run is interrupted (Ctrl-C) - the next step of the session reads whatever is in it.
pack() {
  rm -f inputs/fasta-34.26.5.shar.Z
  mkdir -p returned
  if command -v zip > /dev/null; then rm -f returned/A2-19c-fastapair2-run3-out.zip; zip -qr returned/A2-19c-fastapair2-run3-out.zip out; else rm -f returned/A2-19c-fastapair2-run3-out.tar.gz; tar czf returned/A2-19c-fastapair2-run3-out.tar.gz out; fi
  ls -l returned
}
trap pack EXIT
docker run --rm --platform linux/amd64 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -e DIRECT_FROM="${DIRECT_FROM:-254}" -e SET_TIMEOUT="${SET_TIMEOUT:-1200}" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
ls -l returned
echo "Send back: docs/Validation/review-2026-09/handoff/A2-19c-fastapair2/returned/A2-19c-fastapair2-run3-out.zip (next to the archive of run 2)."
