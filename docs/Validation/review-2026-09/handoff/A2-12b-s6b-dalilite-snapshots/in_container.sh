#!/usr/bin/env bash
# Runs inside the image; never aborts (every step is logged, failures included).
cd /work
mkdir -p out
{ echo "date: $(date -u +%FT%TZ)"; echo "uname -m: $(uname -m)"; python3 --version; curl --version | head -1; } > out/ENV.txt 2>&1
# connectivity facts first (so that a blocked network is visible in the archive)
for u in https://web.archive.org/ https://archive.org/wayback/available?url=example.com https://www.ebi.ac.uk/dali/DaliLite/ http://ekhidna2.biocenter.helsinki.fi/dali/; do
  echo "== $u" >> out/connectivity.log; curl -sS -I -L -A Seqeron-review/1.0 --max-time 60 "$u" >> out/connectivity.log 2>&1 || echo "rc=$?" >> out/connectivity.log
done
python3 -I /opt/fetch_dalilite_snapshots.py /work/out > out/fetch.log 2>&1 || echo "fetch script ended with rc=$?" >> out/fetch.log
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null
ls -l out; tail -20 out/fetch.log
exit 0
