#!/bin/bash
# Work package 1b.9d: a list of cells of run_1b9d_cache.sh, one after the other (the machine does nothing else).
# usage: run_1b9d_cache_sweep.sh <rows file> <work dir> <input name> <input file> <text offset> <text length> <variant D|C|U> <runs> <MB limits...> -- [extra args of the sesa build command]
# cells are named <input>.<variant><limit MB>[.<tag>] ; the tag (env B9D_TAG) distinguishes configurations (e.g. lean)
set -u
ROWS=$1; W=$2; NAME=$3; INPUT=$4; TOFF=$5; TLEN=$6; VAR=$7; RUNS=$8; shift 8
LIMITS=()
while [ "$1" != "--" ]; do LIMITS+=("$1"); shift; done
shift
HERE=$(cd "$(dirname "$0")" && pwd)
TAG=${B9D_TAG:+.$B9D_TAG}
for mb in "${LIMITS[@]}"; do
  for r in $(seq 1 "$RUNS"); do
    if [ "$mb" = max ]; then LIM=max; CELL=$NAME.U$TAG; else LIM=$((mb * 1000000)); CELL=$NAME.$VAR$mb$TAG; fi
    if [ "$VAR" = D ] && [ "$mb" != max ]; then export B9D_RAISE=1073741824; else unset B9D_RAISE; fi
    "$HERE/run_1b9d_cache.sh" "$W" "$ROWS" "$CELL" "$r" "$LIM" "$INPUT" "$TOFF" "$TLEN" "$@"
  done
done
