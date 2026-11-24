#!/usr/bin/env bash
# Runs inside the image (never aborts on a failing step).
cd /work
mkdir -p out out/s4
{
  echo "date: $(date -u +%FT%TZ)"; echo "uname -m: $(uname -m)"; echo "t_coffee: $(t_coffee -version 2>&1 | head -1)"
  echo "DaliLite v3: $( [ -f /opt/dalilite3/BUILD_OK ] && cat /opt/dalilite3/BUILD_OK || echo 'not available (see dalilite3_fetch.log)' )"
  echo "DaliLite tarball: $(sha256sum /opt/dalilite3/src/DaliLite.tar.gz 2>/dev/null)"
} > out/ENV.txt
cp /opt/dalilite3/fetch.log out/dalilite3_fetch.log 2>/dev/null
if [ -f /opt/dalilite3/BUILD_OK ]; then
  mkdir -p /tmp/dali_log
  # the 15 sets of run 10 that stopped (positions in the dali_pair job list); same seeds, same structures (FULL_BACKBONE), one worker as in run 10; t_coffee's stderr / stdout are kept in extra.stderr / extra.stdout
  ONLY=5,8,10,15,16,17,18,20,21,23,24,25,29,32,38 TC_STDERR=1 DALI_REAL=1 FULL_BACKBONE=1 WORKERS=1 python3 -I /work/generate_tcoffee_external_oracle.py /work/out/s4 dali_pair > out/gen_s4.log 2>&1 || echo "s4 failed"
fi
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null
ls -l out
