#!/bin/bash
# E9 gates G1 (identity) and G3 (smoke test) outside the shipped code, on the inputs prepared by run_e9_setup.sh. Per name, in one directory, no timing claims (warm cache, nothing sampled):
#   e9-psi (default parameters)                        cmp with <name>.ref.psi (LinkTables in RAM, PSI[0] = 0)
#   e9-psi --window-kb 1 --symbols-per-pass 2 --read-kb 64   (several passes and small windows on any alphabet)   cmp with <name>.ref.psi
#   e9-psi --cyclic                                    entry 0 differs from the reference only in the first four bytes, which hold the row of the suffix 0 (checked)
#   e9-isa --scale-n 3000000000 (DIST)                 cmp with <name>.ref.isa
#   e9-psi-isa                                         cmp with <name>.ref.psi
#   e9-smoke (dna_10m rep_10m prot_1m, and the 10^6 families with --sample 5000 --link-budget 20000)   every kernel x source of P1 with the in-RAM tables and with the streamed ones, sample A of the P1 benchmark
# usage: run_e9_identity.sh <work dir> <rows file> <names...>
set -u
W=$1; ROWS=$2; shift 2
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e9_run/SuffixTree.EmProto.dll}
cd "$W" || exit 1
printf 'R\tE9\tidentity\thead\t%s\t\tgit rev-parse HEAD of the checkout when the group started; binary %s\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" "$DLL" >> "$ROWS"
purge() { awk -F'\t' -v l="$1" '!($3==l)' "$ROWS" > "$ROWS.tmp" && mv "$ROWS.tmp" "$ROWS"; }
verdict() { # label what file ref
  if cmp -s "$3" "$4"; then V=yes; else V=NO; fi
  printf 'R\tE9\t%s\tcmp_%s\t%s\t\tcmp %s %s (%s bytes)\n' "$1" "$2" "$V" "$3" "$4" "$(stat -L -c %s "$3")" >> "$ROWS"
}
for X in "$@"; do
  L=$X.id
  grep -q -P "^R\tE9\t$L\trun_done\t" "$ROWS" && { echo "$L already done"; continue; }
  purge "$L"
  rm -rf scr_$X.id out.$X.*
  export DOTNET_GCHeapHardLimit=0x6000000
  dotnet "$DLL" e9-psi --bwt $X.bwt --out out.$X.psi --pk-meta $X.pk.meta --case "$L" --rows "$ROWS" > /dev/null; verdict "$L" psi out.$X.psi $X.ref.psi
  dotnet "$DLL" e9-psi --bwt $X.bwt --out out.$X.psi_small --window-kb 1 --symbols-per-pass 2 --read-kb 64 --case "$L.small" --rows "$ROWS" > /dev/null; verdict "$L" psi_small_windows_2_symbols_per_pass out.$X.psi_small $X.ref.psi
  dotnet "$DLL" e9-psi --bwt $X.bwt --out out.$X.psi_cyc --cyclic --case "$L.cyclic" --rows "$ROWS" > /dev/null
  if [ "$(tail -c +5 out.$X.psi_cyc | sha256sum | cut -d' ' -f1)" = "$(tail -c +5 $X.ref.psi | sha256sum | cut -d' ' -f1)" ] && [ "$(stat -L -c %s out.$X.psi_cyc)" = "$(stat -L -c %s $X.ref.psi)" ]; then V=yes; else V=NO; fi
  P=$(awk -F= '$1=="primary" {print $2}' $X.bwt.meta); E=$(od -An -tu4 -N4 out.$X.psi_cyc | tr -d ' ')
  if [ "$E" = "$P" ]; then V2=yes; else V2=NO; fi
  printf 'R\tE9\t%s\tcmp_psi_cyclic_tail\t%s\t\tbytes 5.. of the cyclic Psi equal the reference; entry 0 = %s equals the primary row %s: %s\n' "$L" "$V" "$E" "$P" "$V2" >> "$ROWS"
  N=$(awk -F= '$1=="n" {print $2}' $X.pk.meta); KB=$(( 64 * N / 100000000 )); [ $KB -lt 1 ] && KB=1; [ $KB -gt 64 ] && KB=64
  dotnet "$DLL" e9-isa --sa $X.sa --out out.$X.isa --scratch scr_$X.id --scale-n 3000000000 --io-kb $KB --case "$L" --rows "$ROWS" > /dev/null; verdict "$L" isa_by_dist out.$X.isa $X.ref.isa
  dotnet "$DLL" e9-psi-isa --sa $X.sa --isa out.$X.isa --out out.$X.psi2 --case "$L" --rows "$ROWS" > /dev/null; verdict "$L" psi_from_isa out.$X.psi2 $X.ref.psi
  unset DOTNET_GCHeapHardLimit   # the smoke test holds the in-RAM index and the in-RAM tables of P1: it is a test of the kernels, not a bounded-memory run
  case $X in
    dna_10m|rep_10m|prot_1m|rnd256_10m) dotnet "$DLL" e9-smoke --input $X.txt --sa $X.sa --isa out.$X.isa --psi out.$X.psi --psi2 out.$X.psi2 --case "$L" --rows "$ROWS" > /dev/null ;;
    a_1m|fib_1m|tm_1m)       dotnet "$DLL" e9-smoke --input $X.txt --sa $X.sa --isa out.$X.isa --psi out.$X.psi --psi2 out.$X.psi2 --sample 5000 --link-budget 20000 --case "$L" --rows "$ROWS" > /dev/null ;;
  esac
  printf 'R\tE9\t%s\trun_done\tyes\t\tall rows of this group are written\n' "$L" >> "$ROWS"
  echo "$L: $(awk -F'\t' -v l="$L" '$3==l && $4 ~ /^cmp_/ {printf "%s=%s ", $4, $5}' "$ROWS")"
  rm -rf scr_$X.id out.$X.*
done
echo identity-done
