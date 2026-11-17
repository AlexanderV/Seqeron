#!/bin/bash
# E8 entry point of the image (ENTRYPOINT). UNBUILT AND UNTESTED IN A CONTAINER, ON macOS AND ON arm64.
#   e8.sh preflight [--sizes "1g 2.2g 3g"] [--quick]   what this machine / container is, whether the data volume is fit, what to change
#   e8.sh plan                                          the protocol with expected durations and the state of every step
#   e8.sh gen <size>...                                 generate the genome-like texts (1g 2.2g 3g)
#   e8.sh run a|b|c|e [--runs N] [--threads T] ...      the generous steps (1, 2.2, 3 Gbp) with the full verification, and the thread-count identity
#   e8.sh run d --size 1g --memory-mb M                 one cell of the limited-cache sweep (the container must have been started with that --memory)
#   e8.sh verify <prefix> [<input text>] [options]      the verifier on any set of e6 outputs
#   e8.sh report                                        the tables of RESULTS-E8, the gates, and the return archive (rows, sums, logs)
#   e8.sh clean                                         remove what interrupted runs left behind
#   e8.sh shell                                         a shell in the container
# The data directory is E8_DATA (/data: the named Docker volume). The rows file is $E8_DATA/rows/e8_rows.tsv.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
export E8_DATA=${E8_DATA:-/data}
export E8_DLL=${E8_DLL:-/opt/e8/bin/SuffixTree.EmProto.dll}
cmd=${1:-help}; [ $# -gt 0 ] && shift
case $cmd in
  preflight) exec "$HERE/preflight.sh" "$@" ;;
  plan|gen|clean|cells) exec "$HERE/run_e8.sh" "$cmd" "$@" ;;
  run)
    if [ "${E8_SKIP_PREFLIGHT:-0}" != 1 ]; then
      "$HERE/preflight.sh" --quick --no-rows > "$E8_DATA/preflight.last.txt" 2>&1; rc=$?
      if [ $rc -ne 0 ]; then cat "$E8_DATA/preflight.last.txt"; echo "run refused: fix the [REFUSE] lines above (e8.sh preflight shows the details)" >&2; exit $rc; fi
    fi
    exec "$HERE/run_e8.sh" run "$@" ;;
  verify)
    prefix=${1:?usage: e8.sh verify <prefix> [<input text>] [options]}; shift
    input=; if [ $# -gt 0 ] && [ "${1#--}" = "$1" ]; then input=$1; shift; fi
    exec dotnet "$E8_DLL" e8-verify --prefix "$prefix" ${input:+--input "$input"} --threads "${E8_THREADS:-4}" "$@" ;;
  report)
    ROWS=$E8_DATA/rows/e8_rows.tsv
    [ -f "$ROWS" ] || { echo "no rows file ($ROWS): nothing to report" >&2; exit 1; }
    python3 -I "$HERE/tables_e8.py" "$ROWS" --out "$E8_DATA/rows/RESULTS-E8.tables.md" --json "$E8_DATA/rows/gates.json"
    mkdir -p "$E8_DATA/return"
    ARCHIVE=$E8_DATA/return/e8-return-$(date -u +%Y%m%dT%H%M%SZ).tar.gz
    ( cd "$E8_DATA" && tar -czf "$ARCHIVE" rows sums logs preflight.last.txt 2>/dev/null )
    echo
    echo "return archive: $ARCHIVE ($(du -h "$ARCHIVE" | cut -f1))  = the rows file, the hashes, the logs and the tables. Copy it out with: e8-docker.sh export <folder on the Mac>"
    ;;
  shell) exec /bin/bash ;;
  help|-h|--help|*) sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//' ;;
esac
