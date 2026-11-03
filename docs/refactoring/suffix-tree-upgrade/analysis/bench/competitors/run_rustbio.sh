#!/usr/bin/env bash
# rust-bio (crate `bio`, MIT): suffix_array (SA-IS) + BWT + FM-index (Occ k=3), locate with the full SA.
# Needs: cargo (rustc >= 1.7x) and access to crates.io (static.crates.io / index.crates.io)
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
need cargo
( cd "$HERE/rustbio" && cargo build --release --quiet )
BIN="$HERE/rustbio/target/release/rustbio-bench"
BIOV=$(cd "$HERE/rustbio" && cargo metadata --format-version 1 2>/dev/null | python3 -c "import json,sys;m=json.load(sys.stdin);print([p['version'] for p in m['packages'] if p['name']=='bio'][0])")
RV=$(rustc --version | cut -d' ' -f2)
echo "rust-bio crate bio $BIOV (crates.io; licence MIT per crate metadata), rustc $RV" >> "$WORK/SOURCES.txt"
INPUTS_RB="${INPUTS_RB:-dna_1m dna_10m rep_10m prot_1m}"
for x in $INPUTS_RB; do
  meas --tool rust-bio --version "bio $BIOV rustc-$RV --release lto" --input $x --phase all --timeout "${TMO:-1800}" -- "$BIN" data/$x.txt queries/$x
done
echo "done: rows appended to $RAW"
