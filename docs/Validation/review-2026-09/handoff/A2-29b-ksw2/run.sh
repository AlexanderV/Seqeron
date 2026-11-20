#!/bin/bash
# One command: builds the image, runs the ksw2 functions on cases.txt and writes ./out/ (ksw2_oracle.txt = the result file to send back, plus versions / checksums / logs).
# Apple Silicon: x86_64 emulation is required (the code uses SSE intrinsics); the script passes --platform linux/amd64.
set -eu
cd "$(dirname "$0")"
mkdir -p out
docker build --platform linux/amd64 -t seqeron-ksw2-review .
docker run --rm --platform linux/amd64 -v "$PWD/out:/out" seqeron-ksw2-review
echo "done: $(wc -l < out/ksw2_oracle.txt) result lines in out/ksw2_oracle.txt"
