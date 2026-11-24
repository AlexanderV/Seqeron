#!/usr/bin/env bash
cd /work
mkdir -p out
{ echo "date: $(date -u +%FT%TZ)"; echo "uname -m: $(uname -m)"; } > out/ENV.txt
python3 -I /opt/fetch_old_dalilite.py /work/out > out/fetch.log 2>&1 || echo "fetch failed"
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null
ls -l out; tail -20 out/fetch.log
