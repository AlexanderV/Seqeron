#!/usr/bin/env bash
# Python: pydivsufsort (libdivsufsort bindings, MIT) and `suffix-trees` (pure-Python Ukkonen suffix tree).
# Needs: python3 (venv, pip) and access to PyPI.
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
VENV="$WORK/pyvenv"
[ -d "$VENV" ] || python3 -m venv "$VENV"
"$VENV/bin/pip" install --quiet pydivsufsort suffix-trees
PDV=$("$VENV/bin/pip" show pydivsufsort | awk '/^Version:/{print $2}')
STV=$("$VENV/bin/pip" show suffix-trees | awk '/^Version:/{print $2}')
{ "$VENV/bin/pip" show pydivsufsort suffix-trees numpy | grep -E '^(Name|Version|License|Home-page):'; } >> "$WORK/SOURCES.txt"
PYV=$("$VENV/bin/python" -c 'import sys;print("%d.%d.%d"%sys.version_info[:3])')
ONLY="${ONLY:-all}"   # all | pydiv | stree
if [ "$ONLY" = all ] || [ "$ONLY" = pydiv ]; then
for x in ${INPUTS_PYSA:-dna_1m dna_10m rep_10m prot_1m}; do
  meas --tool pydivsufsort --version "pydivsufsort $PDV py$PYV" --input $x --phase all --timeout "${TMO:-1800}" -- "$VENV/bin/python" "$HERE/python/pydivsufsort_bench.py" data/$x.txt queries/$x
done
fi
# suffix-trees (pure-Python Ukkonen): prefixes of the DNA texts, growing until a size does not complete (BLOCKED row) or hits the memory cap.
# "file:n" pairs; the address-space cap (ulimit -v, STREE_VMEM_KB) turns an out-of-memory situation into a clean MemoryError instead of an OOM kill.
if [ "$ONLY" = all ] || [ "$ONLY" = stree ]; then
for fn in ${STREE_SIZES:-dna_1m:10000 dna_1m:30000 dna_1m:100000 dna_1m:300000 dna_1m:1000000 dna_10m:3000000 dna_10m:10000000 dna_20m:20000000}; do
  f=${fn%%:*}; n=${fn##*:}
  $BL measure --raw "$RAW" --data-dir "$BENCH" --repeats "${STREE_REPEATS:-3}" --timeout "${STREE_TMO:-900}" --tool suffix-trees-py --version "suffix-trees $STV py$PYV" \
      --input "${f}_prefix_$n" --n "$n" --phase all -- bash -c 'ulimit -v "$0"; exec "$@"' "${STREE_VMEM_KB:-12000000}" \
      "$VENV/bin/python" "$HERE/python/stree_bench.py" data/$f.txt queries/$f "$n" || { echo "suffix-trees: $f n=$n did not complete, larger sizes skipped"; break; }
done
fi
echo "done: rows appended to $RAW"
