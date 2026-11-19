#!/usr/bin/env bash
# Diagnostic by design: never aborts; every failure leaves its trace in ./out (build.log, preflight.txt, run.log).
cd /work
{ echo "date: $(date -u)"; echo "uname -m: $(uname -m)"; echo "binary sha256:"; cat /opt/muscle3/muscle3.sha256; echo "base: $(grep PRETTY /etc/os-release)"; \
  echo "version:"; /opt/muscle3/muscle3 -version 2>&1; } > out/ENV.txt
cp /opt/muscle3/fetch_src.sha256 /opt/muscle3/build.log out/ 2>/dev/null
if [ -e /opt/muscle3/BUILD_FAILED ] || [ ! -x /opt/muscle3/muscle3 ]; then
  echo "MUSCLE BUILD FAILED - see out/build.log" | tee out/run.log
else
  { printf '>a\nMKVLAAGIVGLL\n>b\nMKVLAGIVGLL\n>c\nMRVLAAGIIGLL\n>d\nMKVLAAGIVGLI\n' > /tmp/pre.fa; /opt/muscle3/muscle3 -in /tmp/pre.fa -out /tmp/pre.out -quiet 2>&1; echo "exit $?"; cat /tmp/pre.out 2>&1; \
    echo "--- -maxiters 8 -tree1"; /opt/muscle3/muscle3 -in /tmp/pre.fa -out /tmp/pre2.out -quiet -maxiters 8 -tree1 /tmp/t1.txt 2>&1; echo "exit $?"; cat /tmp/pre2.out /tmp/t1.txt 2>&1; } > out/preflight.txt
  python3 run_sets.py 2>&1 | tee out/run.log
fi
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null || true
exit 0
