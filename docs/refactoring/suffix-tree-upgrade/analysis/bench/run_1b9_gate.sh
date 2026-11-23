#!/usr/bin/env bash
# The S1 gate measurement (ST-S1 work package 1b.9): the in-memory build, the query throughput, the matching operations, the container and the in-process A/B runs.
# Every repeat is a FRESH process (lib/benchlib.py: wall time, CPU, peak RSS from wait4 ru_maxrss). Run ONE group at a time on a quiet machine (uptime, ps); nothing else may run.
#   usage: run_1b9_gate.sh <group> [args]
#     drift <label>                          the libsais reference harness (competitors/sa_c/sa_bench.c flow), 5 fresh processes on dna_10m; <label> = start | end
#     build <input> [repeats]                gate1b9 build: prepare / SA-IS / LCP / child scan / total, peak RSS, index bytes, the libsais vector check (default 5)
#     lean <input> [repeats]                 gate1b9 lean: the same stages with nothing read or counted before the build (rows with the operation prefix lean_)
#     steady <input> [repeats]               gate1b9 steady: 5 builds in one process, the last three reported (default 3 processes)
#     entry <input> [repeats]                gate1b9 entry: EsaSuffixTreeFactory.CreateEsa, one call, cold (default 5)
#     query <input> [repeats]                Bench query: 20-mer count (present, absent) and 8-mer locate of the in-memory ESA, one thread (default 5)
#     nav <input> <op> [repeats] [ukkonen]   Bench nav: MEM / MUM / anchors of memq100k on the production navigator, ESA only (default 5); a 4th argument "ukkonen" adds the in-memory Ukkonen tree in the same process
#     container <input> [work dir]           builds the container of a benchmark text with the production builder (sesa build, 4 threads) and prints its sections
#     open <container> <input> [repeats] [cold] [patternsfirst]   gate1b9 open: first (cold JIT) and warm Open; with "cold" the page cache is dropped before every process; "patternsfirst" = the command openp (op prefix p...)
#     patterns <input> <file>                the abapp pattern file of a benchmark input: P = present20, L = locate8 (locate4 on protein), A = absent20
#     ab <rows file> <processes> <pattern file> <passes> <threads> <seconds> <specA> <specB> [label]    abapp in several processes (run_ab.sh)
#   Environment: GATE_BIN (the gate1b9 build), BENCH_BIN (the Esa.Bench build), ESA_BIN (SuffixTree.Esa.dll and its dependencies), ABAPP_BIN, LIBSAIS_BIN, D (the bench dir, default /tmp/stbench),
#   OUT (the scratch directory, default /tmp/stbench/b9). Raw rows -> runs_esa/raw_1b9.tsv, aggregated by `aggregate` into results_1b9_rows.tsv.
#   Build the pieces first:
#     dotnet build apps/SuffixTree.Esa.Bench/SuffixTree.Esa.Bench.csproj -c Release -o $OUT/bin/bench
#     dotnet build src/SuffixTree/Algorithms/SuffixTree.Esa/SuffixTree.Esa.csproj -c Release -o $OUT/bin/esa
#     dotnet build docs/refactoring/suffix-tree-upgrade/analysis/bench/gate1b9/Gate1b9.csproj -c Release -p:EsaDir=$OUT/bin/esa -o $OUT/bin/gate
#     dotnet build docs/refactoring/suffix-tree-upgrade/analysis/bench/abapp/abapp.csproj -c Release -o $OUT/bin/abapp
#     libsais: clone IlyaGrebnov/libsais at 1aae6344b8a58f00626f001b1fe33bdbecd5e469, cmake -DLIBSAIS_USE_OPENMP=ON, then
#       gcc -O3 -fopenmp -DUSE_LIBSAIS -DLIBSAIS_OPENMP -DLIBSAIS_OMP -I<libsais>/include analysis/bench/competitors/sa_c/sa_bench.c <liblibsais.a> -lpthread -o $OUT/sa_bench_libsais
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$B/../../../../.." && pwd)"
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 DOTNET_gcServer=0
unset DOTNET_GCHeapHardLimit
D=${D:-/tmp/stbench}; OUT=${OUT:-/tmp/stbench/b9}
GATE_BIN=${GATE_BIN:-$OUT/bin/gate}; BENCH_BIN=${BENCH_BIN:-$OUT/bin/bench}; ESA_BIN=${ESA_BIN:-$OUT/bin/esa}; ABAPP_BIN=${ABAPP_BIN:-$OUT/bin/abapp}; LIBSAIS_BIN=${LIBSAIS_BIN:-$OUT/sa_bench_libsais}
GATE="dotnet $GATE_BIN/SuffixTree.Esa.Bench.dll"; BENCH="dotnet $BENCH_BIN/SuffixTree.Esa.Bench.dll"
BL="python3 $B/lib/benchlib.py"
RAW=${RAW:-$B/runs_esa/raw_1b9.tsv}
VEC=$B/../sa_vectors.tsv
HEAD=$(git -C "$ROOT" rev-parse --short=12 HEAD)
SV="Seqeron@$HEAD net10.0 Release 1b9 production (SA-IS pointer + direct LCP + child scan, defaults, 1 thread)"
mkdir -p "$B/runs_esa" "$OUT/rows" "$OUT/log"
m() { $BL measure --raw "$RAW" --data-dir "$B" --tool seqeron-esa --version "$SV" "$@"; }

group=${1:-}; shift || true
case $group in
  drift)
    label=${1:?label}
    $BL measure --raw "$B/runs_esa/raw_1b9_libsais.tsv" --data-dir "$B" --tool libsais --version "libsais v2.9.1 (1aae634) gcc-13.3.0 -O3 (OpenMP build; serial rows use libsais()) session $label" \
      --input dna_10m --repeats 5 --phase all -- "$LIBSAIS_BIN" "$D/data/dna_10m.txt" "$D/queries/dna_10m" ;;
  build)   x=$1; m --input "$x" --repeats "${2:-5}" --phase build  -- $GATE build  "$D/data/$x.txt" "$x" "$VEC" ;;
  lean)    x=$1; m --input "$x" --repeats "${2:-5}" --phase lean --op-prefix lean_ -- $GATE lean "$D/data/$x.txt" "$x" "$VEC" ;;
  steady)  x=$1; m --input "$x" --repeats "${2:-3}" --phase steady -- $GATE steady "$D/data/$x.txt" "$x" 5 ;;
  entry)   x=$1; m --input "$x" --repeats "${2:-5}" --phase entry  -- $GATE entry  "$D/data/$x.txt" "$x" ;;
  query)   x=$1; m --input "$x" --repeats "${2:-5}" --phase query  -- $BENCH query --input "$D/data/$x.txt" --name "$x" --queries "$D" --min-seconds 2 ;;
  nav)     x=$1; op=$2; uk=false; ph="nav_$op"; [ "${4:-}" = ukkonen ] && { uk=true; ph="nav_${op}_ukk"; }
           m --input "$x" --repeats "${3:-5}" --phase "$ph" $( [ $uk = true ] && echo --op-prefix ukk_ ) -- $BENCH nav --input "$D/data/$x.txt" --name "$x" --queries "$D" --op "$op" --min-length 20 --calls 10 --ukkonen $uk ;;
  container)
    x=$1; W=${2:-$OUT/work_$x}; mkdir -p "$W"
    $BENCH sesa build --input "$D/data/$x.txt" --work "$W" --label "$x" --threads 4 --keep true --rss true | cut -c1-200
    python3 "$B/../em/sesa_sections.py" "$W/out.sesa" hash --no-crc | cut -c1-120 ;;
  open)
    c=$1; x=$2; r=${3:-5}; cold=${4:-}
    cmd=open; pre=; [ "${5:-}" = patternsfirst ] && { cmd=openp; pre=p; }
    if [ "$cold" = cold ]; then m --input "$x" --repeats "$r" --phase ${pre}open_cold --op-prefix ${pre}cold_ --drop-caches -- $GATE $cmd "$c" "$D/queries/$x.present20.txt"
    else m --input "$x" --repeats "$r" --phase ${pre}open_warmcache --op-prefix ${pre}warm_ -- $GATE $cmd "$c" "$D/queries/$x.present20.txt"; fi ;;
  patterns)
    x=$1; out=$2; loc=locate8; case $x in prot*) loc=locate4 ;; esac
    { sed 's/^/P /' "$D/queries/$x.present20.txt"; sed 's/^/L /' "$D/queries/$x.$loc.txt"; sed 's/^/A /' "$D/queries/$x.absent20.txt"; } > "$out" ;;
  ab)
    ABAPP="$ABAPP_BIN" "$B/abapp/run_ab.sh" "$@" ;;
  aggregate)
    $BL aggregate --raw "$RAW" --out "$B/results_1b9_rows.tsv"
    $BL aggregate --raw "$B/runs_esa/raw_1b9_libsais.tsv" --out "$B/results_1b9_libsais.tsv"
    [ -f "$B/runs_esa/raw_1b9_instr.tsv" ] && $BL aggregate --raw "$B/runs_esa/raw_1b9_instr.tsv" --out "$B/results_1b9_instr.tsv" ;;
  *) echo "usage: see the header of $0" >&2; exit 2 ;;
esac
