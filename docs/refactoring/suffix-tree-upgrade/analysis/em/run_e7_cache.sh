#!/bin/bash
# E7 part C (the whole-pipeline page-cache test of ADR-004 section 6b): one cold run of the end-to-end pipeline inside a cgroup-v1 memory sub-group with a limit, sampled by sample_e7_cache.py.
# usage: run_e7_cache.sh <work dir> <rows file> <name> <label> <run> <limit bytes | max> [extra args of the e6 command...]
#   env: E7_GC          DOTNET_GCHeapHardLimit of the run (default 0x6000000 = 96 MiB, as in E6)
#        E7_WATCHDOG    seconds after which the run is killed and recorded as "> limit" (default 900 = 15 min)
#        E7_ENV         extra "NAME=value" settings for the process (space separated), e.g. MALLOC_ARENA_MAX=1
# The limit is the TOTAL of the sub-group (anonymous memory + page cache incl. mapped file pages + kernel memory charged to it); the process joins the group before it starts, so everything it reads and writes is
# charged to it. The sampler is outside the group. A run whose rows end in run_done is skipped; the sub-group is removed at the end. Outputs are hashed (outside the group, after the run) against <name>.refs.sha.
set -u
W=$1; ROWS=$2; T=$3; L=$4; R=$5; LIMIT=$6; shift 6
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e7_run/SuffixTree.EmProto.dll}
HERE=$(cd "$(dirname "$0")" && pwd)
WD=${E7_WATCHDOG:-900}
cd "$W" || exit 1
if grep -q -P "^R\tE7\t$L\t$R\trun_done\t" "$ROWS"; then echo "$L run $R already done"; exit 0; fi
awk -F'\t' -v l="$L" -v r="$R" '!($3==l && $4==r)' "$ROWS" > "$ROWS.tmp" && mv "$ROWS.tmp" "$ROWS"
CG=/sys/fs/cgroup/memory/emproto_e7_$L.$R
O=out.$L.$R; SCR=scratch.$L.$R; PIDF=pid.$L.$R; SERIES=series.$L.$R.tsv
rm -rf "$SCR" $O.* $PIDF 2>/dev/null
rmdir "$CG" 2>/dev/null
mkdir "$CG" || { echo "cannot create $CG" >&2; exit 2; }
if [ "$LIMIT" != "max" ]; then echo "$LIMIT" > "$CG/memory.limit_in_bytes" || { echo "cannot set the limit" >&2; rmdir "$CG"; exit 2; }; fi
[ -w "$CG/memory.swappiness" ] && echo 0 > "$CG/memory.swappiness"
printf 'R\tE7\t%s\t%s\tsettings\tlimit=%s gc=%s env=%s args=%s\t\tcgroup limit (bytes), DOTNET_GCHeapHardLimit, extra environment, extra program arguments\n' "$L" "$R" "$LIMIT" "${E7_GC:-0x6000000}" "${E7_ENV:-}" "$*" >> "$ROWS"
sync; echo 3 > /proc/sys/vm/drop_caches
python3 -I "$HERE/sample_e7_cache.py" $PIDF "$ROWS" "$L" "$R" "$SCR" "$O.pk" "$CG" "$SERIES" &
SP=$!
sleep 0.2
(
  echo $BASHPID > "$CG/cgroup.procs"
  echo $BASHPID > $PIDF
  export DOTNET_GCHeapHardLimit=${E7_GC:-0x6000000}
  for kv in ${E7_ENV:-}; do export "$kv"; done
  exec dotnet "$DLL" e6 --input "$T.txt" --out "$O" --scratch "$SCR" --case "$L" --run "$R" --rows "$ROWS.prog.$L.$R" "$@" > /dev/null
) &
PID=$!
( sleep "$WD"; touch "$PIDF.watchdog"; kill -9 $PID 2>/dev/null ) &
WDG=$!
wait $PID
RC=$?
kill $WDG 2>/dev/null; wait $WDG 2>/dev/null
wait $SP
FAIL=$(cat "$CG/memory.failcnt"); MAXU=$(cat "$CG/memory.max_usage_in_bytes")
STAT=$(grep -E "^(total_)?(cache|rss|mapped_file|pgpgin|pgpgout|pgmajfault) " "$CG/memory.stat" | tr '\n' ' ')
KILLED=no; [ -e "$PIDF.watchdog" ] && KILLED=watchdog
[ "$RC" = 137 ] && [ "$KILLED" = no ] && KILLED=oom-or-signal
[ -f "$ROWS.prog.$L.$R" ] && { sed 's/^R\tE6\t/R\tE7\t/' "$ROWS.prog.$L.$R" >> "$ROWS"; rm -f "$ROWS.prog.$L.$R"; }
printf 'R\tE7\t%s\t%s\texit_code\t%s\t\tprogram exit code (0 = healthy)\n' "$L" "$R" "$RC" >> "$ROWS"
printf 'R\tE7\t%s\t%s\tkilled\t%s\t\tno | watchdog (%s s) | oom-or-signal (exit 137 without the watchdog: the cgroup OOM killer)\n' "$L" "$R" "$KILLED" "$WD" >> "$ROWS"
printf 'R\tE7\t%s\t%s\tcg_failcnt_final\t%s\tn\tmemory.failcnt after the run\n' "$L" "$R" "$FAIL" >> "$ROWS"
printf 'R\tE7\t%s\t%s\tcg_max_usage_final\t%s\tbytes\tmemory.max_usage_in_bytes after the run; memory.stat: %s\n' "$L" "$R" "$MAXU" "$STAT" >> "$ROWS"
rmdir "$CG" 2>/dev/null
# identity (outside the group, after the run)
if [ "$RC" = 0 ]; then
  for K in pk sa lcp ovf cld; do ( sha256sum "$O.$K" 2>/dev/null | cut -d' ' -f1 > "$O.$K.sha" ) & done; wait
  for K in pk sa lcp ovf cld; do
    GOT=$(cat "$O.$K.sha"); WANT=$(awk -v k=$K '$1==k {print $2}' "$T.refs.sha")
    if [ -n "$GOT" ] && [ "$GOT" = "$WANT" ]; then V=yes; else V=NO; fi
    printf 'R\tE7\t%s\t%s\tcmp_%s\t%s\t\tSHA-256 %s against the reference hash in %s.refs.sha\n' "$L" "$R" "$K" "$V" "${GOT:0:16}" "$T" >> "$ROWS"
  done
  if [ "$(grep -E '^(nodecount|internal|lrsdepth|lrsrank)=' "$O.stats" 2>/dev/null | sort | sed 's/^/stat /')" = "$(grep '^stat ' "$T.refs.sha")" ]; then V=yes; else V=NO; fi
  printf 'R\tE7\t%s\t%s\tcmp_stats\t%s\t\tstatistics against %s.refs.sha\n' "$L" "$R" "$V" "$T" >> "$ROWS"
fi
printf 'R\tE7\t%s\t%s\tscratch_left\t%s\t\tentries of %s after the run (killed runs leave their scratch; it is removed here)\n' "$L" "$R" "$( [ -e "$SCR" ] && ls -A "$SCR" | wc -l || echo 0)" "$SCR" >> "$ROWS"
printf 'R\tE7\t%s\t%s\trun_done\tyes\t\tall rows of this run are written\n' "$L" "$R" >> "$ROWS"
echo "$L run $R rc=$RC killed=$KILLED failcnt=$FAIL $(grep -P "^R\tE7\t$L\t$R\t(total_s|eff_resident_ph3|rssanon_peak)\t" "$ROWS" | awk -F'\t' '{printf "%s=%.4g ", $5, $6}')"
rm -rf "$SCR" $O.* $PIDF $PIDF.watchdog
