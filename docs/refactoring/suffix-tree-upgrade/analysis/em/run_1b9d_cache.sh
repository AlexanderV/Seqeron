#!/bin/bash
# Work package 1b.9d, the page-cache-limited run of the PRODUCTION chain (the E7 method of run_e7_cache.sh: one cold run inside a cgroup-v1 memory sub-group with a limit, sampled from outside by sample_1b9d_cache.py),
# the process being `SuffixTree.Esa.Bench sesa build` = SesaBuilder.Build with its default parameters and four threads (the options after the limit are added to the command, e.g. --scale-n 3000000000).
# usage: run_1b9d_cache.sh <work dir> <rows file> <case> <run> <limit bytes | max> <input text file> <text offset> <text length> [extra args of the sesa build command...]
#   env: B9D_DLL       the Bench assembly (default /tmp/stbench/b9d/bin/SuffixTree.Esa.Bench.dll)
#        B9D_GC        DOTNET_GCHeapHardLimit of the run (default: unset, as for every production build of the package)
#        B9D_WATCHDOG  seconds after which the run is killed and recorded as "> limit" (default 600)
#        B9D_RAISE     limit in bytes to set when PH3 has finished (the rank file is deleted): "limited PH0-PH3, PH4 and later not limited" (variant D of E7); default: the limit stays for the whole run (variant C)
#        B9D_PREWARM   1 (default): read the .NET runtime image into the page cache from outside the sub-group first (its pages are then charged to the parent group), as in E7
#        B9D_REF       a file that holds the SHA-256 of the container of the first (unlimited) run of this input; created by the run when it does not exist, compared with when it does
# The limit is the TOTAL of the sub-group (anonymous memory + page cache incl. mapped pages + kernel memory charged to it), as in E7. The text is the TEXT section of the container <work>/w.<case>.<run>/out.sesa.tmp
# (offset and length are the arguments), the rank file is <...>/out.sesa.tmp.d/ranks.bin. A run whose rows end in run_done is skipped; the sub-group is removed at the end.
set -u
W=$1; ROWS=$2; L=$3; R=$4; LIMIT=$5; T=$6; TOFF=$7; TLEN=$8; shift 8
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
unset DOTNET_GCHeapHardLimit
DLL=${B9D_DLL:-/tmp/stbench/b9d/bin/SuffixTree.Esa.Bench.dll}
HERE=$(cd "$(dirname "$0")" && pwd)
WD=${B9D_WATCHDOG:-600}
cd "$W" || exit 1
touch "$ROWS"
if grep -q -P "^R\t1b9d\t$L\t$R\trun_done\t" "$ROWS"; then echo "$L run $R already done"; exit 0; fi
awk -F'\t' -v l="$L" -v r="$R" '!($3==l && $4==r)' "$ROWS" > "$ROWS.tmp" && mv "$ROWS.tmp" "$ROWS"
CG=/sys/fs/cgroup/memory/b9d_$L.$R
D=$W/w.$L.$R; SCR=$D/out.sesa.tmp.d; PIDF=$W/pid.$L.$R; SERIES=$W/series.$L.$R.tsv; PROG=$W/prog.$L.$R
rm -rf "$D" $PIDF $PROG 2>/dev/null
mkdir -p "$D"
rmdir "$CG" 2>/dev/null
mkdir "$CG" || { echo "cannot create $CG" >&2; exit 2; }
if [ "$LIMIT" != "max" ]; then echo "$LIMIT" > "$CG/memory.limit_in_bytes" || { echo "cannot set the limit" >&2; rmdir "$CG"; exit 2; }; fi
[ -w "$CG/memory.swappiness" ] && echo 0 > "$CG/memory.swappiness"
printf 'R\t1b9d\t%s\t%s\tsettings\tlimit=%s raise=%s prewarm=%s gc=%s args=%s\t\tcgroup limit (bytes), limit after PH3, runtime pre-read outside the group, DOTNET_GCHeapHardLimit, extra program arguments\n' "$L" "$R" "$LIMIT" "${B9D_RAISE:-none}" "${B9D_PREWARM:-1}" "${B9D_GC:-unset}" "$*" >> "$ROWS"
sync; echo 3 > /proc/sys/vm/drop_caches
if [ "${B9D_PREWARM:-1}" = 1 ]; then
  FW=$(dirname "$(dotnet --list-runtimes | awk '/Microsoft.NETCore.App/ {gsub(/[\[\]]/,"",$3); print $3 "/" $2 "/x"; exit}')")
  dotnet "$DLL" > /dev/null 2>&1
  cat "$FW"/*.so "$FW"/System.Private.CoreLib.dll "$FW"/System.Runtime.dll "$FW"/System.Collections.dll "$FW"/System.Linq.dll "$FW"/System.Threading.dll "$FW"/System.IO*.dll "$FW"/System.Memory.dll "$FW"/System.Runtime.InteropServices.dll "$FW"/System.Diagnostics*.dll "$FW"/System.Text*.dll "$FW"/System.Security.Cryptography.dll "$FW"/System.Numerics*.dll "$FW"/System.Console.dll "$(dirname "$DLL")"/*.dll > /dev/null 2>&1
fi
python3 -I "$HERE/sample_1b9d_cache.py" $PIDF "$ROWS" "$L" "$R" "$SCR" "$D/out.sesa.tmp" "$TOFF" "$TLEN" "$CG" "$SERIES" &
SP=$!
sleep 0.2
(
  echo $BASHPID > "$CG/cgroup.procs"
  echo $BASHPID > $PIDF
  [ -n "${B9D_GC:-}" ] && export DOTNET_GCHeapHardLimit=$B9D_GC
  exec dotnet "$DLL" sesa build --input "$T" --work "$D" --label "$L.$R" --keep true --threads 4 "$@" > "$PROG" 2>"$PROG.err"
) &
PID=$!
( exec > /dev/null 2>&1; for ((w = 0; w < WD; w++)); do kill -0 $PID 2>/dev/null || exit 0; sleep 1; done; touch "$PIDF.watchdog"; kill -9 $PID 2>/dev/null ) &   # polls once a second, ends with the run
WDG=$!
if [ -n "${B9D_RAISE:-}" ]; then
  ( exec > /dev/null 2>&1
    while kill -0 $PID 2>/dev/null && [ ! -e "$SCR/ranks.bin" ]; do sleep 0.05; done   # PH1 has made the rank file (the directories ph1 and ph3 exist from the start of the production build, unlike the prototype's)
    while kill -0 $PID 2>/dev/null && [ -e "$SCR/ranks.bin" ]; do sleep 0.05; done     # the rank file is deleted as soon as PH3 has used it
    kill -0 $PID 2>/dev/null && echo "$B9D_RAISE" > "$CG/memory.limit_in_bytes" && touch "$PIDF.raised" ) &
  RZ=$!
fi
wait $PID
RC=$?
wait $WDG 2>/dev/null
[ -n "${RZ:-}" ] && wait $RZ 2>/dev/null
wait $SP
FAIL=$(cat "$CG/memory.failcnt"); MAXU=$(cat "$CG/memory.max_usage_in_bytes")
STAT=$(grep -E "^(total_)?(cache|rss|mapped_file|pgpgin|pgpgout|pgmajfault) " "$CG/memory.stat" | tr '\n' ' ')
KILLED=no; [ -e "$PIDF.watchdog" ] && KILLED=watchdog
[ "$RC" = 137 ] && [ "$KILLED" = no ] && KILLED=oom-or-signal
[ -f "$PROG" ] && awk -F'\t' -v l="$L" -v r="$R" 'BEGIN{OFS="\t"} $1=="R" && $2=="sesa" {print "R","1b9d",l,r,$3,$4,$5,$6}' "$PROG" >> "$ROWS"
printf 'R\t1b9d\t%s\t%s\texit_code\t%s\t\tprogram exit code (0 = healthy)\n' "$L" "$R" "$RC" >> "$ROWS"
[ -e "$PIDF.raised" ] && printf 'R\t1b9d\t%s\t%s\tlimit_raised\tyes\t\tthe limit was raised to %s bytes when PH3 had finished\n' "$L" "$R" "${B9D_RAISE:-}" >> "$ROWS"
printf 'R\t1b9d\t%s\t%s\tkilled\t%s\t\tno | watchdog (%s s) | oom-or-signal (exit 137 without the watchdog: the cgroup OOM killer)\n' "$L" "$R" "$KILLED" "$WD" >> "$ROWS"
printf 'R\t1b9d\t%s\t%s\tcg_failcnt_final\t%s\tn\tmemory.failcnt after the run\n' "$L" "$R" "$FAIL" >> "$ROWS"
printf 'R\t1b9d\t%s\t%s\tcg_max_usage_final\t%s\tbytes\tmemory.max_usage_in_bytes after the run; memory.stat: %s\n' "$L" "$R" "$MAXU" "$STAT" >> "$ROWS"
rmdir "$CG" 2>/dev/null
# identity (outside the group, after the run): the container of an unlimited run is the reference of the others
if [ "$RC" = 0 ] && [ -f "$D/out.sesa" ]; then
  GOT=$(sha256sum "$D/out.sesa" | cut -d' ' -f1)
  REF=${B9D_REF:-$W/$(basename "$T").sesa.sha}
  if [ ! -f "$REF" ]; then echo "$GOT" > "$REF"; V=reference; elif [ "$GOT" = "$(cat "$REF")" ]; then V=yes; else V=NO; fi
  printf 'R\t1b9d\t%s\t%s\tcmp_container\t%s\t\tSHA-256 %s of the container against the reference in %s (reference = this run made it)\n' "$L" "$R" "$V" "${GOT:0:16}" "$REF" >> "$ROWS"
fi
printf 'R\t1b9d\t%s\t%s\tleft_behind\t%s\t\tentries of the work directory of the run besides the container (a build that succeeded leaves the container only)\n' "$L" "$R" "$( [ -d "$D" ] && ls -A "$D" | grep -v '^out.sesa$' | wc -l || echo 0)" >> "$ROWS"
printf 'R\t1b9d\t%s\t%s\trun_done\tyes\t\tall rows of this run are written\n' "$L" "$R" >> "$ROWS"
echo "$L run $R rc=$RC killed=$KILLED failcnt=$FAIL $(grep -P "^R\t1b9d\t$L\t$R\t(total_s|eff_resident_ph3|ext_rssanon_peak)\t" "$ROWS" | awk -F'\t' '{printf "%s=%.4g ", $5, $6}')"
rm -rf "$D" $PIDF $PIDF.watchdog $PIDF.raised $PROG $PROG.err
