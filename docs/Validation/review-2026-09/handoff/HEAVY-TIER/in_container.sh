#!/usr/bin/env bash
# Runs inside the container. /src = read-only checkout of the repository, /out = result directory (bind-mounted).
# Never aborts: every step writes its log and the final summary is always produced.
set -uo pipefail
MODE="${MODE:-HEAVY}"          # HEAVY | FAST | MCP | ALL
WORKERS="${WORKERS:-}"         # optional: limit NUnit workers (Docker memory < 12 GB: use 2-4)
mkdir -p /out /work
{ echo "date: $(date -u)"; echo "uname -m: $(uname -m)"; echo "nproc: $(nproc)"; free -m | head -2; dotnet --version; echo "MODE=$MODE WORKERS=$WORKERS"; } > /out/ENV.txt
echo "== copying the checkout (no .git, no bin/obj)"
rsync -a --delete --exclude .git --exclude bin --exclude obj /src/ /work/repo/ > /out/rsync.log 2>&1
cd /work/repo || exit 1
{ cat /src/.git/HEAD 2>/dev/null; git -C /src rev-parse HEAD 2>/dev/null; } >> /out/ENV.txt 2>&1 || true
echo "== build"
dotnet build Seqeron.sln -c Debug > /out/build.log 2>&1; echo "build rc=$?" | tee -a /out/summary.txt
EXTRA=(-- "NUnit.AssemblySelectLimit=1000000")
[ -n "$WORKERS" ] && EXTRA+=("NUnit.NumberOfTestWorkers=$WORKERS")
run() {  # name, project, filter
  local name="$1" proj="$2" filt="${3:-}"
  echo "== $name"
  local t0=$(date +%s)
  if [ -n "$filt" ]; then
    dotnet test "$proj" -c Debug --no-build --filter "$filt" --logger "trx;LogFileName=$name.trx" --results-directory /out/trx "${EXTRA[@]}" > "/out/$name.log" 2>&1
  else
    dotnet test "$proj" -c Debug --no-build --logger "trx;LogFileName=$name.trx" --results-directory /out/trx "${EXTRA[@]}" > "/out/$name.log" 2>&1
  fi
  local rc=$? t1=$(date +%s)
  echo "$name rc=$rc seconds=$((t1-t0)) :: $(grep -E '^(Passed!|Failed!)' "/out/$name.log" | tail -1)" | tee -a /out/summary.txt
  grep -E "^\s+Failed " "/out/$name.log" | sed 's/^\s*//' > "/out/$name.failed.txt" 2>/dev/null || true
}
G=tests/Seqeron/Seqeron.Genomics.Tests/Seqeron.Genomics.Tests.csproj
HEAVYF='FullyQualifiedName~.Fuzzing.|FullyQualifiedName~.Properties.|FullyQualifiedName~.Metamorphic.'
FASTF='FullyQualifiedName!~.Fuzzing.&FullyQualifiedName!~.Properties.&FullyQualifiedName!~.Metamorphic.'
case "$MODE" in
  HEAVY|ALL) run genomics_heavy "$G" "$HEAVYF" ;;
esac
case "$MODE" in
  FAST|ALL) run genomics_fast "$G" "$FASTF" ;;
esac
case "$MODE" in
  MCP|ALL)
    for p in tests/Seqeron/Seqeron.Mcp.*.Tests/*.csproj tests/SuffixTree/SuffixTree.Mcp.Core.Tests/*.csproj tests/Seqeron/Seqeron.Genomics.Integration/*.csproj; do
      [ -f "$p" ] && run "$(basename "$(dirname "$p")")" "$p"
    done ;;
esac
echo "== done" | tee -a /out/summary.txt
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" /out 2>/dev/null || true
