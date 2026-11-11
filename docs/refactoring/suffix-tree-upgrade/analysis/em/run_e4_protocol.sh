#!/bin/bash
# E4 protocol (outside the shipped code): PH4 (the LCP array from a suffix array file with bounded RAM) on the inputs prepared by run_e4_setup.sh.
# usage: run_e4_protocol.sh <work dir with <name>.pk, <name>.sa, <name>.oracle.lcp/.ovf> <rows file> <group>
# groups: g100 (100 Mbp random and genome-like: 3 runs of each mode, threads 1 / 2, two more chunk sizes, rounds 1 / 2 / 4 of the bounded mode),
#         small (10 M inputs, a^n 50 M, Fibonacci, tandem, N run: the same set), fam (the small committed adversarial families at 10^5 / 10^6),
#         slope (10 / 50 / 100 / 200 M random DNA and 10 / 50 / 100 M genome-like, plain and with every buffer touched, both modes)
# Every run: dropped page cache (cold), fresh process, DOTNET_GCHeapHardLimit = 96 MiB, scale-down emulation by n / 3e9 (R = ceil(8,000,000 n / 3e9) => P = 375 chunks), output compared byte for byte with
# the LcpKernel oracle (--verify). The label of a case names mode, rounds, threads and chunk size; the run column is the repetition.
set -u
W=$1; ROWS=$2; GROUP=$3
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e4_run/SuffixTree.EmProto.dll}
cd "$W" || exit 1
mkdir -p scratch
printf 'R\tE4\tprotocol\thead\t%s\t\tgit rev-parse HEAD of the binary under test; group %s\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" "$GROUP" >> "$ROWS"
run() { # input label run extra-args...
  T=$1; L=$2; R=$3; shift 3
  sync; echo 3 > /proc/sys/vm/drop_caches
  dotnet "$DLL" e4 --text "$T.pk" --sa "$T.sa" --out "$T.out" --scratch scratch --scale-n 3000000000 --verify "$T.oracle" --case "$L" --run "$R" --rows "$ROWS" "$@" > /dev/null
  echo "$L run $R rc=$?"
}
chunk() { # input factor -> R scaled by n/3e9 (rounded up) times factor
  N=$(grep '^n=' "$1.pk.meta" | cut -d= -f2)
  python3 -I -c "import math; print(max(3, math.ceil(8000000*$N/3000000000*$2)))"
}
full_set() { # input
  T=$1
  for i in 1 2 3; do run $T "$T.dist.t4" $i --threads 4; done
  for i in 1 2 3; do run $T "$T.b2.t4" $i --threads 4 --mode bounded --rounds 2; done
  for i in 1 2 3; do run $T "$T.b1.t4" $i --threads 4 --mode bounded --rounds 1; done
  run $T "$T.b4.t4" 1 --threads 4 --mode bounded --rounds 4
  run $T "$T.dist.t1" 1 --threads 1
  run $T "$T.dist.t2" 1 --threads 2
  run $T "$T.b2.t1" 1 --threads 1 --mode bounded --rounds 2
  run $T "$T.b2.t2" 1 --threads 2 --mode bounded --rounds 2
  run $T "$T.dist.t4.Rhalf" 1 --threads 4 --chunk "$(chunk $T 0.5)"
  run $T "$T.dist.t4.Rdouble" 1 --threads 4 --chunk "$(chunk $T 2)"
  run $T "$T.b2.t4.Rhalf" 1 --threads 4 --mode bounded --rounds 2 --chunk "$(chunk $T 0.5)"
  run $T "$T.b2.t4.Rdouble" 1 --threads 4 --mode bounded --rounds 2 --chunk "$(chunk $T 2)"
}
case $GROUP in
g100)
  for T in dna_100m genome_100m; do full_set $T; done ;;
small)
  for T in dna10 rep10 a50 fib10 tandem nrun10 genome10; do full_set $T; done ;;
fam)
  for T in fam_a_1e5 fam_a_1e6 fam_fib_1e5 fam_fib_1e6 fam_tm_1e5 fam_tm_1e6 fam_tandem_1e5 fam_nrun_1e6 dna1m prot1m; do
    run $T "$T.dist.t4" 1 --threads 4
    run $T "$T.dist.t1" 1 --threads 1
    run $T "$T.dist.t2" 1 --threads 2
    run $T "$T.b2.t4" 1 --threads 4 --mode bounded --rounds 2
    run $T "$T.b1.t1" 1 --threads 1 --mode bounded --rounds 1
    run $T "$T.dist.t4.Rhalf" 1 --threads 4 --chunk "$(chunk $T 0.5)"
    run $T "$T.dist.t4.Rdouble" 1 --threads 4 --chunk "$(chunk $T 2)"
    run $T "$T.b2.t4.Rdouble" 1 --threads 4 --mode bounded --rounds 2 --chunk "$(chunk $T 2)"
  done ;;
slope)
  for T in dna10r dna50r dna_100m dna200r genome10 genome50 genome_100m; do
    for i in 1 2 3; do run $T "$T.dist.slope" $i --threads 4; done
    for i in 1 2 3; do run $T "$T.dist.slopetouch" $i --threads 4 --touch-buffers; done
    for i in 1 2 3; do run $T "$T.b2.slope" $i --threads 4 --mode bounded --rounds 2; done
    for i in 1 2 3; do run $T "$T.b2.slopetouch" $i --threads 4 --mode bounded --rounds 2 --touch-buffers; done
  done ;;
*) echo "unknown group $GROUP"; exit 2 ;;
esac
echo protocol-done
