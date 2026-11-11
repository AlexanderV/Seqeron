#!/bin/bash
# E3b protocol (outside the shipped code): the external-memory PH1 (difference-cover rank file) on the inputs of the prototype work package.
# usage: run_e3b_protocol.sh <work dir with the packed texts *.pk and the in-RAM oracle rank files *.oracle> <rows file> <group>
# groups: g100 (100 Mbp: 3 runs of the baseline, threads 1 / 2, three bucket capacities), small (10 M inputs and adversarial families), slope (10 / 50 / 100 / 200 M, plain and with --touch-buffers),
#         fam (the small committed adversarial families at 10^5 and 10^6 with their libsais names: threads 1 / 2 / 4, three bucket capacities), big (unscaled parameters at 520 Mbp a^n, 1 / 2 Gbp: see the RESULTS-E3B file)
# Every run: dropped page cache (cold), fresh process, DOTNET_GCHeapHardLimit = 96 MiB, scale-down emulation by n / 3e9, output compared byte for byte with the oracle (--verify).
set -u
W=$1; ROWS=$2; GROUP=$3
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e3b_run/SuffixTree.EmProto.dll}
cd "$W" || exit 1
mkdir -p scratch
run() { # text label run extra-args...
  T=$1; L=$2; R=$3; shift 3
  sync; echo 3 > /proc/sys/vm/drop_caches
  dotnet "$DLL" e3b --text "$T.pk" --out "$T.rank" --scratch scratch --scale-n 3000000000 --verify "$T.oracle" --case "$L" --run "$R" --rows "$ROWS" "$@" > /dev/null
  echo "$L run $R rc=$?"
}
prep() { # family label [data dir]: pack and in-RAM oracle (no heap limit for these two commands)
  ( unset DOTNET_GCHeapHardLimit
    [ -f "$2.pk" ] || dotnet "$DLL" pack --family "$1" --data /tmp/stbench/data --out "$2.pk" > /dev/null
    [ -f "$2.oracle" ] || dotnet "$DLL" mkranks --text "$2.pk" --out "$2.oracle" > /dev/null )
}
case $GROUP in
g100)
  for T in dna_100m genome_100m; do
    for i in 1 2 3; do run $T "$T.t4" $i --threads 4; done
    run $T "$T.t1" 1 --threads 1
    run $T "$T.t2" 1 --threads 2
    run $T "$T.t4.Khalf" 1 --threads 4 --k-s 17476
    run $T "$T.t4.Kdouble" 1 --threads 4 --k-s 69905
  done ;;
small)
  for T in dna_10m rep10 a10 fib10 nrun10 genome10 tandem; do
    for i in 1 2 3; do run $T "$T.t4" $i --threads 4; done
    run $T "$T.t1" 1 --threads 1
    run $T "$T.t2" 1 --threads 2
  done ;;
slope)
  for T in dna10 dna50 dna_100m dna200 genome10 genome50 genome_100m; do
    for i in 1 2 3; do run $T "$T.slope" $i --threads 4; done
  done ;;
slopetouch)
  for T in dna10 dna50 dna_100m dna200 genome10 genome50 genome_100m; do
    for i in 1 2 3; do run $T "$T.slopetouch" $i --threads 4 --touch-buffers; done
  done ;;
fam)
  while read -r FAM LABEL; do
    prep "$FAM" "$LABEL"
    run $LABEL "$LABEL.t4" 1 --threads 4
    run $LABEL "$LABEL.t1" 1 --threads 1
    run $LABEL "$LABEL.t2" 1 --threads 2
    N=$(grep '^n=' "$LABEL.pk.meta" | cut -d= -f2)
    run $LABEL "$LABEL.t4.Kx2" 1 --threads 4 --k-s $(python3 -I -c "print(max(3, round(1048576*$N/3000000000*2)))")
    run $LABEL "$LABEL.t4.Khalf" 1 --threads 4 --k-s $(python3 -I -c "print(max(3, round(1048576*$N/3000000000/2)))")
  done <<'LIST'
a:100000 fam_a_1e5
a:1000000 fam_a_1e6
fib:100000 fam_fib_1e5
fib:1000000 fam_fib_1e6
tm:100000 fam_tm_1e5
tm:1000000 fam_tm_1e6
dna_1m fam_dna_1m
prot_1m fam_prot_1m
tandem:171x600 fam_tandem_1e5
nrun:1000000:100000 fam_nrun_1e6
LIST
  ;;
*) echo "unknown group $GROUP"; exit 2 ;;
esac
echo protocol-done
