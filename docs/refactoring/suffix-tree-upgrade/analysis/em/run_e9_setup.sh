#!/bin/bash
# E9 setup (outside the shipped code): for every input name the files the Psi experiments start from, and the in-RAM references of gate G1, with the SHA-256 of every reference.
#   <name>.sa        suffix array file (n + 1 little-endian int32): the `oracle` command (SaisKernel in RAM) for the small inputs; for dna_100m the SA of the E6 pipeline (e6 --scale-n 3000000000), checked against the
#                    hash of the E4-E7 references (/tmp/stbench/e7/dna_100m.refs.sha)
#   <name>.pk (.meta) packed text (the `pack` command; for dna_100m the E6 pipeline's, whose .meta carries PH0's symbol counts)
#   <name>.txt       one byte per symbol (the input of the pipeline and of the smoke test)
#   <name>.bwt (.meta) the BWT of the pass under test (made by e9-bwt, so it is an output of the code under test; its correctness is gate G1 through the Psi it yields)
#   <name>.ref.isa / <name>.ref.psi   ISA and Psi of the P1 prototype (LinkTables in RAM), PSI[0] = 0
#   <name>.refs.sha  "<kind> <sha256>"  (sa pk isa psi)
# Names: dna_10m prot_1m rep_10m a_1m fib_1m tm_1m dna_100m rnd256_10m (10^7 bytes of random.Random(1).randbytes: sigma = 256, the multi-pass case). Every step is marked by <name>.step.<kind> when complete, so an interrupted step is redone from scratch.
# usage: run_e9_setup.sh <work dir> <rows file> <names...>
set -u
W=$1; ROWS=$2; shift 2
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e9_run/SuffixTree.EmProto.dll}
DATA=${E9_DATA:-/tmp/stbench/data}
E6TXT=${E9_E6_TEXT:-/tmp/stbench/e6/dna_100m.txt}
cd "$W" || exit 1
printf 'R\tE9\tsetup\thead\t%s\t\tgit rev-parse HEAD of the checkout when the inputs were made; binary %s\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" "$DLL" >> "$ROWS"
hash1() { sha256sum "$1" | cut -d' ' -f1; }

for X in "$@"; do
  unset DOTNET_GCHeapHardLimit
  case $X in
    a_1m)   F=a:1000000 ;;
    fib_1m) F=fib:1000000 ;;
    tm_1m)  F=tm:1000000 ;;
    *)      F=$X ;;
  esac
  DX=$DATA
  if [ "$X" = rnd256_10m ]; then
    [ -f $X.step.gen ] || { rm -f $X.txt; python3 -I -c 'import random, sys; sys.stdout.buffer.write(random.Random(1).randbytes(10**7))' > $X.txt && touch $X.step.gen; }
    DX=$W
  fi
  if [ "$X" = dna_100m ]; then
    if [ ! -f "$X.step.sa" ]; then
      rm -f $X.sa $X.pk $X.pk.meta $X.lcp $X.ovf $X.cld $X.stats; rm -rf scr_$X
      ln -sf "$E6TXT" $X.txt
      dotnet "$DLL" e6 --input $X.txt --out $X --scratch scr_$X --threads 4 --scale-n 3000000000 --case $X.setup --rows "$ROWS" > /dev/null || { echo "$X: the E6 pipeline failed" >&2; exit 3; }
      rm -f $X.lcp $X.ovf $X.cld $X.stats; rm -rf scr_$X
      touch $X.step.sa $X.step.pk
    fi
    WANT_SA=$(awk '$1=="sa" {print $2}' /tmp/stbench/e7/dna_100m.refs.sha); WANT_PK=$(awk '$1=="pk" {print $2}' /tmp/stbench/e7/dna_100m.refs.sha)
    [ "$(hash1 $X.sa)" = "$WANT_SA" ] || { echo "$X: the SA differs from the E4-E7 reference hash" >&2; exit 3; }
    [ "$(hash1 $X.pk)" = "$WANT_PK" ] || { echo "$X: the packed text differs from the E4-E7 reference hash" >&2; exit 3; }
    printf 'R\tE9\t%s\tsetup\tsa_pk_equal_e7_reference_hashes\tyes\t\tSHA-256 of the pipeline SA and packed text equal the hashes of /tmp/stbench/e7/dna_100m.refs.sha\n' "$X" >> "$ROWS"
  else
    [ -f "$X.step.txt" ] || { [ "$X" = rnd256_10m ] || rm -f "$X.txt"; case $X in dna_10m|prot_1m|rep_10m) ln -sf "$DATA/$X.txt" "$X.txt" ;; rnd256_10m) : ;; *) dotnet "$DLL" e6-mktext --family "$F" --data "$DX" --out "$X.txt" --case "$X" --rows "$ROWS" > /dev/null ;; esac && touch "$X.step.txt"; }
    [ -f "$X.step.sa" ] || { rm -f "$X.sa"; dotnet "$DLL" oracle --family "$F" --data "$DX" --out "$X.sa" --case "$X" --rows "$ROWS" > /dev/null && touch "$X.step.sa"; }
    [ -f "$X.step.pk" ] || { rm -f "$X.pk" "$X.pk.meta"; dotnet "$DLL" pack --family "$F" --data "$DX" --out "$X.pk" --case "$X" --rows "$ROWS" > /dev/null && touch "$X.step.pk"; }
  fi
  if [ ! -f "$X.step.ref" ]; then
    rm -f $X.ref.isa $X.ref.psi
    dotnet "$DLL" e9-oracle --sa $X.sa --out $X.ref --case $X --rows "$ROWS" > /dev/null && touch $X.step.ref
  fi
  if [ ! -f "$X.step.bwt" ]; then
    rm -f $X.bwt $X.bwt.meta
    dotnet "$DLL" e9-bwt --text $X.pk --sa $X.sa --out $X.bwt --case $X.setup --rows "$ROWS" > /dev/null && touch $X.step.bwt
  fi
  { echo "sa $(hash1 $X.sa)"; echo "pk $(hash1 $X.pk)"; echo "isa $(hash1 $X.ref.isa)"; echo "psi $(hash1 $X.ref.psi)"; echo "bwt $(hash1 $X.bwt)"; } > $X.refs.sha
  echo "$X ready"
done
echo setup-done
