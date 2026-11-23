#!/bin/bash
# runs INSIDE the container: downloads the pinned sources, verifies their checksums, compiles the harness (SSE4.1, SSE2-only and scalar builds) and runs the cases
set -u
mkdir -p /work/mm2 /work/ksw2 /out
cd /work
UA="Seqeron-review/1.0"
MM2=https://raw.githubusercontent.com/lh3/minimap2/v2.26
KSW=https://raw.githubusercontent.com/lh3/ksw2/master
for f in ksw2.h ksw2_extz2_sse.c ksw2_extd2_sse.c; do curl -fsSL -A "$UA" -o mm2/$f $MM2/$f || echo "DOWNLOAD FAILED $f"; done
for f in ksw2.h ksw2_extz.c ksw2_extd.c; do curl -fsSL -A "$UA" -o ksw2/$f $KSW/$f || echo "DOWNLOAD FAILED $f"; done
cp /in/harness.c mm2/; cp /in/harness.c ksw2/
(sha256sum mm2/* ksw2/*; uname -m; gcc --version | head -1; grep -m1 -o -w -E 'sse4_1|sse2' /proc/cpuinfo | sort -u; date -u +%FT%TZ) > /out/versions.txt 2>&1
# expected checksums of the files the Seqeron port was read from (sha256 of what was fetched on 2026-11-19; a mismatch is reported, not fatal)
cat > /work/expected.sha256 <<'EOS'
e6499bda435ebd55940d2fa5f4c70c108b603c5077587fb5528ed1fbe9adfaf0  mm2/ksw2.h
74c491c834396b574df1246beb27186118e4e7a2d3373b66134b713998409f7d  mm2/ksw2_extz2_sse.c
46911ad50e3b03e02354371bb756f6c0ce30d0f22feb87c30822874abe7452f1  mm2/ksw2_extd2_sse.c
23265ae6e7e62ffb3a8c92e14b3f65ea9539f3f46eb8f47f9c393c0f48b33a20  ksw2/ksw2_extz.c
006d64ffe710a3e207c0686ebb33b9f2596011d54ee9e0e6d5c1a4e1a4ef3bfa  ksw2/ksw2_extd.c
EOS
sha256sum -c /work/expected.sha256 > /out/checksums.txt 2>&1 || echo "CHECKSUM MISMATCH (see checksums.txt)" >> /out/versions.txt
(cd mm2 && gcc -O2 -msse4.1 -o ../h_sse41 harness.c ksw2_extz2_sse.c ksw2_extd2_sse.c > /out/build_sse41.log 2>&1)
(cd mm2 && gcc -O2 -msse2 -DKSW_SSE2_ONLY -o ../h_sse2 harness.c ksw2_extz2_sse.c ksw2_extd2_sse.c > /out/build_sse2.log 2>&1)
(cd ksw2 && gcc -O2 -DSCALAR -o ../h_scalar harness.c ksw2_extz.c ksw2_extd.c > /out/build_scalar.log 2>&1)
# run_all BINARY OUT ERR: runs the whole case file; when the program under test aborts (the last "CASE id" line of the error file names the case), records the id in OUT.crashes and resumes after it
run_all() {
  local bin=$1 outf=$2 errf=$3 start=0 rc
  : > "$outf"; : > "$errf"; : > "$outf.crashes"
  while :; do
    "$bin" "$start" < /in/cases.txt >> "$outf" 2>> "$errf"; rc=$?
    [ $rc -eq 0 ] && break
    local last; last=$(grep '^CASE ' "$errf" | tail -1 | cut -d' ' -f2)
    echo "$last $rc" >> "$outf.crashes"
    [ -z "$last" ] && break
    start=$((last + 1))
  done
}
run_all ./h_sse41 /out/ksw2_sse41.txt /out/run_sse41.err
run_all ./h_sse2  /out/ksw2_sse2.txt  /out/run_sse2.err
run_all ./h_scalar /out/ksw2_scalar.txt /out/run_scalar.err
# diagnosis of any crash: the same scalar harness under AddressSanitizer / UBSan, first report only (the case id is the last CASE line)
(cd ksw2 && gcc -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -DSCALAR -o ../h_scalar_asan harness.c ksw2_extz.c ksw2_extd.c > /out/build_scalar_asan.log 2>&1)
./h_scalar_asan 0 < /in/cases.txt > /dev/null 2> /out/asan_scalar.txt || true
(cd mm2 && gcc -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -msse4.1 -o ../h_sse41_asan harness.c ksw2_extz2_sse.c ksw2_extd2_sse.c > /out/build_sse41_asan.log 2>&1)
./h_sse41_asan 0 < /in/cases.txt > /dev/null 2> /out/asan_sse41.txt || true
(command -v valgrind > /dev/null && valgrind -q --error-limit=no ./h_scalar 0 < /in/cases.txt 2>&1 > /dev/null | head -200 > /out/valgrind_scalar.txt) || true
cat /out/ksw2_sse41.txt /out/ksw2_scalar.txt | sort -n -s -k1,1 > /out/ksw2_oracle.txt
wc -l /out/*.txt >> /out/versions.txt
