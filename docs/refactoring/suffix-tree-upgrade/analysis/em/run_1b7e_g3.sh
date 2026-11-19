#!/bin/bash
# Gate G3 / G5 of work package 1b.7e (outside the shipped code). One run = one process (dotnet test of the explicit class SesaBuildCostTests, which times the phases and samples RssAnon every 10 ms with the production
# sampler; an independent sampler, sample_rss_1b7a.py, polls the test host from outside), cold cache (sync + drop_caches), DOTNET_GCHeapHardLimit = 96 MiB and DOTNET_gcConcurrent = 0 as in the E-series, 4 threads, the
# variants alternate within a run number so that neither always goes first.
#   usage: run_1b7e_g3.sh <rows file> <inputs: dna genome> [runs per cell, default 3] [cells]
#   cells: C-A-cold  the whole container build (container, PSI on, fsync off) against the sum of the separate phases + the checksum pass (separate), geometry A (ScaleN 3x10^9: the E6 parameters)
#          C-B-cold  the same at the real geometry (ScaleN 0)
#          P6-cold   PH6 alone: the prototype of E9 (ph6proto, windows of 1 MiB), the same with its own 64 KB windows (ph6proto64) and the production (ph6prod), on the E6 files of the input
#          G5-psi    container with PSI against without (A cold)         G5-sparse  Auto (fallocate) against sparse creation of the container (A cold)
#          G5-fsync  fsync of the publication on against off (A cold)      G5-copy    the separate phases + the copy of their outputs into one file (A cold)      G5-crc  the checksum pass, hardware and table (warm, one build)
set -u
ROWS=$(readlink -f "$1")
INPUTS=${2:-"dna genome"}
RUNS=${3:-3}
CELLS=${4:-"C-A-cold C-B-cold P6-cold"}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
PROJ=$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj
D=${G3_DIR:-/tmp/stbench/c7}
W=${G3_WORK:-/tmp/stbench/g3e}
mkdir -p "$W"
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" > "$ROWS"
printf '#container: label variant scan create pack ph1 ph23 ph45 bwt psi crc seal total read write mapped scratch_peak arena_peak rss_ph0 rss_ph1 rss_ph23 rss_ph45 rss_ph6 rss_seal rss_whole verdict file_size grown reservation overflow lcpkind rss_before_kb rss_max_kb samples\n' >> "$ROWS"
printf '#separate: label variant ph0 ph1 ph23 ph45 bwt psi crc total read write mapped scratch_peak bytes copy_s copied rss_before_kb rss_max_kb samples\n' >> "$ROWS"
printf '#ph6: label variant bwt_s psi_s read1 write1 mapped1 read2 write2 mapped2 windows primary rss_before_kb rss_max_kb samples   (production: bwt_read bwt_write bwt_mapped psi_read psi_write psi_mapped windows primary)\n' >> "$ROWS"

meta() { grep "^$2=" "$1.pk.meta" | cut -d= -f2; }

run() { # name label variant scale cache run [extra env...]
  NAME=$1; L=$2; V=$3; SC=$4; CACHE=$5; K=$6
  case $NAME in dna) TXT=/tmp/stbench/e6/dna_100m.txt; BASE=$D/dna_100m ;; *) TXT=/tmp/stbench/e7/${NAME}_100m.txt; BASE=$D/${NAME}_100m ;; esac
  rm -rf "$W"/*
  sync
  if [ "$CACHE" = cold ]; then echo 3 > /proc/sys/vm/drop_caches; else cat "$BASE.pk" > /dev/null; fi
  python3 -I "$HERE/sample_rss_1b7a.py" "$ROWS" "$L.$CACHE.$K" "$V" &
  SP=$!
  env G3_INPUT=$TXT G3_TARGET=$W/out.sesa G3_WORK=$W G3_ROWS=$ROWS G3_LABEL="$L.$CACHE.$K" G3_SCALE=$SC G3_THREADS=4 \
      G3_PK=$BASE.pk G3_SA=$BASE.sa G3_N=$(meta "$BASE" n) G3_BITS=$(meta "$BASE" bits) G3_SIGMA=$(meta "$BASE" sigma) G3_COUNTS=$(meta "$BASE" counts) \
      ${EXTRA:-} G3_VARIANT=$V \
    dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaBuildCostTests" > "$W/last.log" 2>&1 || { echo "run failed: $L $V $CACHE $K"; tail -20 "$W/last.log"; }
  wait $SP
  for F in q.bwt q.psi p.bwt p.psi; do [ -f "$W/$F" ] && printf '#sha\t%s\t%s\t%s\t%s\n' "$L.$CACHE.$K" "$V" "$F" "$(sha256sum "$W/$F" | cut -d' ' -f1)" >> "$ROWS"; done
  rm -rf "$W"/*
}

for NAME in $INPUTS; do
  for CELL in $CELLS; do
    case $CELL in
      C-A-cold)   SC=3000000000; VARS="container separate"; EX0="G3_FSYNC=0" ;;
      C-B-cold)   SC=0;          VARS="container separate"; EX0="G3_FSYNC=0" ;;
      P6-cold)    SC=0;          VARS="ph6proto ph6prod ph6proto64"; EX0="" ;;
      G5-psi)     SC=3000000000; VARS="container containernopsi"; EX0="G3_FSYNC=0" ;;
      G5-sparse)  SC=3000000000; VARS="container containersparse"; EX0="G3_FSYNC=0" ;;
      G5-fsync)   SC=3000000000; VARS="container containernosync"; EX0="" ;;
      G5-copy)    SC=3000000000; VARS="separate separatecopy"; EX0="G3_FSYNC=0" ;;
      *) echo "unknown cell $CELL"; continue ;;
    esac
    for K in $(seq 1 "$RUNS"); do
      set -- $VARS
      ORDER=$(for i in $(seq 0 $(($#-1))); do eval echo \${$(( (i + K - 1) % $# + 1 ))}; done)
      for V in $ORDER; do
        EXTRA="$EX0"; VV=$V
        case $V in
          containernopsi)  VV=container; EXTRA="$EX0 G3_PSI=0" ;;
          containersparse) VV=container; EXTRA="$EX0 G3_RESERVATION=sparse" ;;
          containernosync) VV=container; EXTRA="G3_FSYNC=0" ;;
          separatecopy)    VV=separate;  EXTRA="$EX0 G3_COPY=1" ;;
          ph6proto64)      VV=ph6proto;  EXTRA="G3_WINDOW=65536" ;;
          container)       [ "$CELL" = G5-fsync ] && EXTRA="G3_FSYNC=1" ;;
        esac
        EXTRA=$EXTRA run "$NAME" "$NAME.$CELL.$V" "$VV" "$SC" cold "$K"
      done
    done
  done
done
echo g3-done
