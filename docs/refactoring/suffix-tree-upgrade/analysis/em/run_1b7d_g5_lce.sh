#!/bin/bash
# G5 of work package 1b.7d (not a gate): the cost of the clamped MappedTextReader.Lce in the comparator of PH3 against an unchecked variant. The unchecked variant exists only in a throwaway worktree (LceUnchecked: the same
# loop without the check of the arguments and without the clamp to the shorter suffix, called by SuffixOrder, which proves its arguments in range); the production tree is not changed. The production variant of the cost test
# (SaBuilder.Build, real geometry, cold cache, 4 threads) runs alternately in the two builds, 3 times each per input; both write the suffix array, which is cmp'ed with the reference.
#   usage: run_1b7d_g5_lce.sh <rows file> <worktree for the unchecked build> [inputs, default "dna genome"] [runs, default 3]
set -u
ROWS=$(readlink -f "$1"); U=$2; INPUTS=${3:-"dna genome"}; RUNS=${4:-3}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
D=/tmp/stbench/c7
W=/tmp/stbench/g3c
mkdir -p "$W"
[ -d "$U" ] || git -C "$ROOT" worktree add --detach "$U" HEAD > /dev/null
git -C "$U" checkout --detach -q "$(git -C "$ROOT" rev-parse HEAD)"
python3 -I - "$U" <<'PYEOF'
import sys
u = sys.argv[1] + '/src/SuffixTree/Algorithms/SuffixTree.Esa/External/'
p = u + 'MappedTextReader.cs'
s = open(p).read()
if 'LceUnchecked' not in s:
    a = s.index('    internal long Lce(long i, long j, long max, ref long steps)')
    b = s.index('        long k = 0;', a)
    e = s.index('    /// <summary>\n    /// Lexicographic comparison', a)
    body = s[b:e]
    s = s[:e] + '    internal long LceUnchecked(long i, long j, long max, ref long steps)\n    {\n' + body + s[e:]
    open(p, 'w').write(s)
p = u + 'Sort/SuffixOrder.cs'
s = open(p).read()
s = s.replace('t.Lce(i + s, j + s, lim - s, ref c.LceSteps)', 't.LceUnchecked(i + s, j + s, lim - s, ref c.LceSteps)')
open(p, 'w').write(s)
PYEOF
for T in "$ROOT" "$U"; do dotnet build "$T/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj" -c Release --nologo -v q -p:RunAnalyzers=false -p:EnforceCodeStyleInBuild=false 2>&1 | grep -E " error " | head -3; done
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" > "$ROWS"
printf '#label\tvariant\tph1_s\tph23_s\ttotal_s\tph1_read\tph1_write\tph1_mapped\tph23_read\tph23_write\tph23_mapped\tscratch_peak_ph1\tscratch_peak\tarena_ph1\tarena_ph23\tlevels\tverdict\trssanon_before_kb\trssanon_max_kb\trss_samples\n' >> "$ROWS"
meta() { grep "^$2=" "$1.pk.meta" | cut -d= -f2; }
run() { # name tree variantlabel run
  NAME=$1; T=$2; V=$3; K=$4
  N=$(meta "$D/${NAME}_100m" n); BITS=$(meta "$D/${NAME}_100m" bits)
  rm -rf "$W/scratch" "$W/out".*; sync; echo 3 > /proc/sys/vm/drop_caches
  G3_VARIANT=prod G3_PK=$D/${NAME}_100m.pk G3_N=$N G3_BITS=$BITS G3_SA=$W/out.sa G3_SCRATCH=$W/scratch G3_SCALE=0 G3_ROWS=$ROWS G3_LABEL="$NAME.lce.$V.$K" \
    dotnet test "$T/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj" -c Release --no-build --nologo --filter "FullyQualifiedName~SaBuildCostTests" > "$W/last.log" 2>&1 || { echo "run failed: $NAME $V $K"; tail -20 "$W/last.log"; }
  cmp -s "$W/out.sa" "$D/${NAME}_100m.sa" && R=equal-to-reference || R="DIFFERENT"
  echo "#cmp	$NAME.lce.$V.$K	$V	$R" >> "$ROWS"
  rm -rf "$W/scratch" "$W/out".*
}
for NAME in $INPUTS; do
  for K in $(seq 1 "$RUNS"); do
    if [ $((K % 2)) = 1 ]; then A="clamped $ROOT"; B="unchecked $U"; else A="unchecked $U"; B="clamped $ROOT"; fi
    set -- $A; run "$NAME" "$2" "$1" "$K"
    set -- $B; run "$NAME" "$2" "$1" "$K"
  done
done
echo g5-lce-done
