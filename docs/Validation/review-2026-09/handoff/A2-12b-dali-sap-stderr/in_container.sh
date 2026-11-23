#!/usr/bin/env bash
# Runs inside the image (never aborts on a failing step).
cd /work
mkdir -p out
{
  echo "date: $(date -u +%FT%TZ)"; echo "uname -m: $(uname -m)"; echo "t_coffee: $(t_coffee -version 2>&1 | head -1)"
  echo "SAP build: $( [ -f /opt/sap/BUILD_FAILED ] && echo FAILED || cat /opt/sap/sap.sha256 )"
  echo "DaliLite v5 tarball: $(cat /opt/dali5/tarball.sha256 2>/dev/null) build: $( [ -f /opt/dali5/BUILD_FAILED ] && echo FAILED || echo ok )"
  echo "DaliLite v5 bin/ has a DaliLite command: $( [ -e /opt/dali5/DaliLite.v5/bin/DaliLite ] && echo YES || echo NO )"
  echo "DaliLite v3: $( [ -f /opt/dalilite3/BUILD_OK ] && cat /opt/dalilite3/BUILD_OK || echo 'not available (see dalilite3_fetch.log)' )"
} > out/ENV.txt
cp /opt/dalilite3/fetch.log out/dalilite3_fetch.log 2>/dev/null; cp /opt/dali5/bin.list out/dali5_bin.list 2>/dev/null; cp /opt/dali5/bin_after_build.list out/dali5_bin_after_build.list 2>/dev/null; cp /opt/dali5/build.log out/dali5_build.log 2>/dev/null
cp /opt/sap/build.log out/sap_build.log 2>/dev/null
mkdir -p out/s1 out/s2 out/s3 out/s4 out/s5
# (1) the SAME 130 sets as the integrated A2-12 run (feature sap_real, same seeds), with t_coffee's stderr / stdout recorded for every set: why the 5 sets with a missing SAP result stop (exit 1).
export SAP_REAL=1 TC_STDERR=1 PATH=/opt/sap/bin:$PATH
python3 -I /work/generate_tcoffee_external_oracle.py /work/out/s1 sap_real > out/gen_s1.log 2>&1 || echo "s1 failed"
# (2) -mode 3dcoffee with the real SAP (8 sets; with a network the name look-ups may work)
python3 -I /work/generate_tcoffee_external_oracle.py /work/out/s2 sap_3dcoffee > out/gen_s2.log 2>&1 || echo "s2 failed"
unset SAP_REAL
# (3) dali_pair (the 60 sets of the committed oracle, same seeds) with DaliLite v5's bin/ on the PATH: no `DaliLite` command there, so t_coffee must stop as it does without DaliLite
PATH=/opt/dali5/DaliLite.v5/bin:$PATH python3 -I /work/generate_tcoffee_external_oracle.py /work/out/s3 dali_pair > out/gen_s3.log 2>&1 || echo "s3 failed"
# (4) dali_pair with the real DaliLite v3 behind the logging wrapper (only when fetch_dali3.py built it); extra.dali = every call of `DaliLite` (arguments, return code, stdout, stderr)
if [ -f /opt/dalilite3/BUILD_OK ]; then
  mkdir -p /tmp/dali_log
  DALI_REAL=1 python3 -I /work/generate_tcoffee_external_oracle.py /work/out/s4 dali_pair > out/gen_s4.log 2>&1 || echo "s4 failed"
fi
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null
ls -l out
