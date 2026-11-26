#!/usr/bin/env bash
# Runs inside the image (PID 1 is tini: `docker run --init`). Never aborts (every step is logged, failures included). Everything the script prints goes to the console AND out/fetch.log.
# SIGTERM / SIGINT (docker stop, Ctrl-C): forwarded to the python script, which writes report.txt / summary.json and exits; then this script ends too.
cd /work
mkdir -p out
{ echo "date: $(date -u +%FT%TZ)"; echo "uname -m: $(uname -m)"; python3 --version; curl --version | head -1; } > out/ENV.txt 2>&1
# connectivity facts first (so that a blocked network is visible in the archive)
for u in https://web.archive.org/ https://archive.org/wayback/available?url=example.com https://www.ebi.ac.uk/dali/DaliLite/ http://ekhidna2.biocenter.helsinki.fi/dali/; do
  echo "== $u" >> out/connectivity.log; curl -sS -I -L -A Seqeron-review/1.0 --max-time 60 "$u" >> out/connectivity.log 2>&1 || echo "rc=$?" >> out/connectivity.log
done
python3 -I -u /opt/fetch_dalilite_snapshots.py /work/out > >(tee -a out/fetch.log) 2>&1 &
PY=$!
trap 'echo "[in_container] signal received: asking the script to write its summary" >&2; kill -TERM $PY 2>/dev/null' TERM INT
wait $PY; rc=$?
while kill -0 $PY 2>/dev/null; do wait $PY; rc=$?; done      # a trapped signal interrupts the first wait; the script is still writing its summary
[ "$rc" -ne 0 ] && echo "fetch script ended with rc=$rc" | tee -a out/fetch.log
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null
ls -l out | head -30; tail -5 out/report.txt 2>/dev/null | cut -c1-300
exit 0
