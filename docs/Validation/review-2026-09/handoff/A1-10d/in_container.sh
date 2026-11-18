#!/usr/bin/env bash
set -euo pipefail
cd /work
{ echo "date: $(date -u)"; echo "uname -m: $(uname -m)"; echo "binary sha256:"; cat /opt/clustalw/clustalw.sha256; echo "base: $(grep PRETTY /etc/os-release)"; } > out/ENV.txt
python3 make_inputs.py > /dev/null 2>&1 || true   # inputs/sets.tsv is committed; regenerated only if the oracle TSVs are mounted
python3 run_sets.py | tee out/run.log
cp /opt/clustalw/CommandLineParser.patched.cpp out/ 2>/dev/null || true
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null || true
