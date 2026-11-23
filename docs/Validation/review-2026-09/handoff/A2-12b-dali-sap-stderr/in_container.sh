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
# RUN 10 diagnostics (1 minute): runs 8 and 9 ended with `sh: puu: not found` / `puutos: not found` although `ldd -r` finds the programs and only `G77_fdate_0` was missing (the shim now has it).
# What the kernel / loader / shell say about the shipped Fortran programs, directly and through sh (the way DaliLite calls them).
D=/data/backup/zope_data/dali_server/DaliLite_3.3/Bin
{
  echo "== install path"; ls -ld /data/backup/zope_data/dali_server/DaliLite_3.3; ls -l $D/puu $D/puutos $D/dssp $D/dalicon $D/DaliLite 2>&1
  for f in puu puutos dalicon dssp; do
    echo "== $f"; file -L $D/$f; readelf -l $D/$f 2>&1 | grep -i -A1 interp; readelf -d $D/$f 2>&1 | grep -i "NEEDED\|RPATH\|RUNPATH"; ldd -r $D/$f 2>&1
    echo "-- direct"; (cd /tmp && $D/$f < /dev/null > /tmp/diag_out.txt 2> /tmp/diag_err.txt; echo "rc=$?"; head -c 400 /tmp/diag_out.txt; head -c 400 /tmp/diag_err.txt)
    echo "-- via sh -c"; (cd /tmp && sh -c "$D/$f < /dev/null" > /tmp/diag_out.txt 2> /tmp/diag_err.txt; echo "rc=$?"; head -c 400 /tmp/diag_out.txt; head -c 400 /tmp/diag_err.txt)
    echo "-- strace execve"; (cd /tmp && strace -f -e trace=execve,openat sh -c "$D/$f < /dev/null" 2>&1 | grep -v ENOENT | tail -12)
    echo "-- strace failures"; (cd /tmp && strace -f -e trace=execve,openat sh -c "$D/$f < /dev/null" 2>&1 | grep ENOENT | tail -8)
  done
  echo "== loader"; ls -l /lib64/ld-linux-x86-64.so.2 /usr/lib/x86_64-linux-gnu/libg2c.so.0 /usr/lib/x86_64-linux-gnu/libf2c.so.2 2>&1
  echo "== g2c shim"; cat /opt/dalilite3/g2cshim.c 2>&1
} > out/diag_puu.txt 2>&1
[ "${DIAG_ONLY:-0}" = 1 ] && { chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null; echo "diag only: out/diag_puu.txt"; exit 0; }
# Steps (1) - (3) are FINAL (runs 1 - 6 returned and integrated by B09 WP68); they run only with RUN_S1_S3=1 (docker run -e RUN_S1_S3=1 ...).
if [ "${RUN_S1_S3:-0}" = 1 ]; then
# (1) the SAME 130 sets as the integrated A2-12 run (feature sap_real, same seeds), with t_coffee's stderr / stdout recorded for every set: why the 5 sets with a missing SAP result stop (exit 1).
export SAP_REAL=1 TC_STDERR=1 PATH=/opt/sap/bin:$PATH
python3 -I /work/generate_tcoffee_external_oracle.py /work/out/s1 sap_real > out/gen_s1.log 2>&1 || echo "s1 failed"
# (2) -mode 3dcoffee with the real SAP (8 sets; with a network the name look-ups may work)
python3 -I /work/generate_tcoffee_external_oracle.py /work/out/s2 sap_3dcoffee > out/gen_s2.log 2>&1 || echo "s2 failed"
unset SAP_REAL
# (3) dali_pair (the 60 sets of the committed oracle, same seeds) with DaliLite v5's bin/ on the PATH: no `DaliLite` command there, so t_coffee must stop as it does without DaliLite
PATH=/opt/dali5/DaliLite.v5/bin:$PATH python3 -I /work/generate_tcoffee_external_oracle.py /work/out/s3 dali_pair > out/gen_s3.log 2>&1 || echo "s3 failed"
fi
# (4) dali_pair with the real DaliLite v3.3 behind the logging wrapper (only when fetch_dali3.py built it); extra.dali = every call of `DaliLite` (arguments, return code, stdout, stderr; call_*.txt.saved = the files that DaliLite removes,
#     the directory listing and the files left in the working directory). RUN 7: FULL_BACKBONE=1 - the structures have N, CA, C, O of every residue (the CA-only files of runs 5 and 6 gave "suspect pdb1.pdb is CA-only", no alignment).
if [ -f /opt/dalilite3/BUILD_OK ]; then
  mkdir -p /tmp/dali_log
  DALI_REAL=1 FULL_BACKBONE=1 WORKERS=1 python3 -I /work/generate_tcoffee_external_oracle.py /work/out/s4 dali_pair > out/gen_s4.log 2>&1 || echo "s4 failed"
fi
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null
ls -l out
