#!/usr/bin/env bash
# Diagnostic by design: never aborts; every failure leaves its trace in ./out.
cd /work
{ echo "date: $(date -u)"; echo "uname -m: $(uname -m)"; echo "base: $(grep PRETTY /etc/os-release)"; echo "binaries:"; cat /opt/cw/binaries.sha256; ls /opt/cw; } > out/ENV.txt
cp /opt/cw/build_*.log /opt/cw/fetch_src.sha256 out/ 2>/dev/null
if ! ls /opt/cw/clustalw_O0 > /dev/null 2>&1; then
  echo "CLUSTALW BUILD FAILED - see out/build_O0.log" | tee out/run.log
else
  { printf '>a\nACGTACGTAAGCT\n>b\nACGTACGTCAGCT\n>c\nACGTTCGTCAGGT\n>d\nTCGTTCGTCAGGT\n>e\nACGTACGTCAGCTAA\n' > /tmp/pre.fa
    echo "--- clustalw_O0 UPGMA + iteration preflight"; /opt/cw/clustalw_O0 -INFILE=/tmp/pre.fa -TYPE=DNA -ALIGN -CLUSTERING=UPGMA -ITERATION=alignment -OUTPUT=FASTA -OUTFILE=/tmp/pre.out -OUTORDER=INPUT > /tmp/pre.log 2>&1; echo "exit $?"; tail -5 /tmp/pre.log; cat /tmp/pre.out; } > out/preflight.txt 2>&1
  python3 -I run_sets.py 2>&1 | tee out/run.log
fi
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null || true
exit 0
