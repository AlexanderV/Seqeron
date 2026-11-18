#!/usr/bin/env bash
set -euo pipefail
cd /work
{ echo "date: $(date -u)"; echo "uname -m: $(uname -m)"; echo "binary sha256:"; cat /opt/muscle3/muscle3.sha256; echo "base: $(grep PRETTY /etc/os-release)"; \
  echo "version:"; /opt/muscle3/muscle3 -version 2>&1 || true; } > out/ENV.txt
cp /opt/muscle3/fetch_src.sha256 out/
python3 run_sets.py | tee out/run.log
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null || true
