#!/usr/bin/env bash
# ONE-MINUTE DIAGNOSTIC (after the image is built): what the system says about the shipped Fortran programs of DaliLite 3.3 (`puu`, `puutos`, `dalicon`, `dssp`). Returns out/diag_puu.txt (+ the fetch log).
set -euo pipefail
cd "$(dirname "$0")"
IMAGE=seqeron-a2-12b:dali-sap-stderr
command -v docker > /dev/null || { echo "docker not found: install Docker Desktop and start it" >&2; exit 1; }
cp ../../../../../scripts/alignment/generate_tcoffee_external_oracle.py ./generate_tcoffee_external_oracle.py
docker build --platform linux/amd64 -t "$IMAGE" .
rm -rf out; mkdir -p out returned
docker run --rm --platform linux/amd64 -e DIAG_ONLY=1 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$PWD":/work "$IMAGE" /work/in_container.sh 2>&1 | tee out/run_console.log
rm -f returned/A2-12b-diag.zip returned/A2-12b-diag.tar.gz
if command -v zip > /dev/null; then zip -qr returned/A2-12b-diag.zip out; else tar czf returned/A2-12b-diag.tar.gz out; fi
ls -l returned
