#!/bin/bash
# E8 runner (INSIDE the container; also works on a plain Linux host). UNBUILT AND UNTESTED IN A CONTAINER, ON macOS AND ON arm64.
# One process per run, rows appended to the rows file, a run counts only after its run_done row; an interrupted run (container restart, Mac asleep) is purged and repeated; the scratch is always cleaned.
# usage:
#   run_e8.sh plan                          the protocol, the expected durations, the abandon rules and what is done / pending in the rows file
#   run_e8.sh gen <size>...                 generate the texts (sizes: 1g 2.2g 3g 500m 123456; decimal units)
#   run_e8.sh run a|b|c [options]           generous runs of the three sizes (step a: 1 Gbp, b: 2.2 Gbp, c: 3 Gbp), each followed by the full verification
#   run_e8.sh run e [options]               thread-count identity at the size of step a (1 and 2 threads, generous memory), hashes only
#   run_e8.sh run d --size 1g --memory-mb M [options]
#                                           one cell of the limited-cache sweep: the CONTAINER must have been started with --memory M MiB (--memory-swap equal); hashes against the generous run of the same text
#   run_e8.sh clean                         remove the leftovers of interrupted runs (work directory)
# options: --runs N (replicates of the cell, default 1)  --threads T (default min(4, cores))  --keep (keep the outputs of every run)  --no-verify (generous steps: skip the verifier; the gate G1 then fails)
#          --lce-cap C (verifier, default 65536)
# environment: E8_DATA (/data), E8_DLL, E8_SIZE_A / _B / _C (1000000000 / 2200000000 / 3000000000), E8_THREADS, E8_GC (0x6000000 = 96 MiB, as E6 / E7), E8_SCALE_N (3000000000), E8_WATCHDOG_S (hard limit of a generous run, 14400),
#              E8_ABANDON_FACTOR (6), E8_EXTRA (extra arguments of the e6 command), E8_KEEP (1 = keep outputs), E8_IMAGE_COMMIT
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=e8lib.sh
. "$HERE/e8lib.sh"
DATA=${E8_DATA:-/data}
DLL=${E8_DLL:-/opt/e8/bin/SuffixTree.EmProto.dll}
ROWS=$DATA/rows/e8_rows.tsv
TEXTDIR=$DATA/text; WORK=$DATA/work; LOGS=$DATA/logs; SUMS=$DATA/sums
SIZE_A=${E8_SIZE_A:-1000000000}; SIZE_B=${E8_SIZE_B:-2200000000}; SIZE_C=${E8_SIZE_C:-3000000000}
CORES=$(nproc 2>/dev/null || echo 1); THREADS=${E8_THREADS:-$(( CORES > 4 ? 4 : CORES ))}
GC=${E8_GC:-0x6000000}; SCALE_N=${E8_SCALE_N:-3000000000}
HARD_WATCHDOG=${E8_WATCHDOG_S:-14400}; ABANDON=${E8_ABANDON_FACTOR:-6}
KEEP=${E8_KEEP:-0}; VERIFY=1; RUNS=1; LCE_CAP=65536; MEMORY_MB=; SIZE_ARG=
mkdir -p "$TEXTDIR" "$WORK" "$LOGS" "$SUMS" "$(dirname "$ROWS")"
touch "$ROWS"
CHILD_PIDS=""
cleanup_children() { for p in $CHILD_PIDS; do kill -9 "$p" 2>/dev/null; done; }
trap 'cleanup_children; echo "interrupted: the run in progress is purged and repeated by the next call" >&2; exit 130' INT TERM

est_minutes() { awk -v n="$1" 'BEGIN { printf "%.0f", n * 250e-9 / 60 }'; }          # ESTIMATE: 250 ns per symbol on the prototype machine (4 vCPU x86, virtio disk); E7: 205-300 ns
est_verify_minutes() { awk -v n="$1" 'BEGIN { m = n * 60e-9 / 60; printf "%.0f", (m < 1 ? 1 : m) }'; }   # measured 58 ns per symbol at 100 Mbp with everything cached (4 threads); more when the files are not cached

grep_done() { grep -q -P "^R\tE8\t$1\t$2\trun_done\t" "$ROWS"; }
purge() { awk -F'\t' -v l="$1" -v r="$2" '!($3==l && $4==r)' "$ROWS" > "$ROWS.tmp" && mv "$ROWS.tmp" "$ROWS"; }
rowval() { awk -F'\t' -v l="$1" -v r="$2" -v m="$3" '$3==l && $4==r && $5==m {print $6; exit}' "$ROWS"; }
addrow() { printf 'R\tE8\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "${5:-}" "${6:-}" >> "$ROWS"; }   # case run metric value unit notes

drop_caches() {
  sync
  if [ -w "$E8_PROC/sys/vm/drop_caches" ] && echo 3 > "$E8_PROC/sys/vm/drop_caches" 2>/dev/null; then echo yes; else echo no; fi
}

do_gen() { # n
  local n=$1 lbl file
  lbl=$(lbl_n "$n"); file=$TEXTDIR/genome_$lbl.txt
  if [ -f "$file" ] && [ -f "$file.sums" ] && [ "$(stat -c %s "$file")" = "$n" ] && grep -q "^n=$n\$" "$file.sums"; then
    echo "text $lbl exists ($(human "$n"), sums file present): skipped"; return 0
  fi
  local avail; avail=$(df_avail_bytes "$DATA")
  if [ "$avail" -lt $(( n + 2000000000 )) ]; then echo "REFUSED: $(human $((n + 2000000000))) needed for the text of $n symbols, $(human "$avail") free" >&2; return 4; fi
  echo "generating $lbl ($n symbols): ESTIMATE about $(awk -v n="$n" 'BEGIN { s = n * 10e-9; printf "%.0f", (s < 5 ? 5 : s) }') s of CPU plus the disk write (measured here: 7 ns per symbol, 200 Mbp in 1.4 s)"
  dotnet "$DLL" e8-gen --n "$n" --out "$file" --force --case "gen.$lbl" --run 0 --rows "$ROWS" > "$LOGS/gen.$lbl.log" 2>&1 || { echo "generation failed, see $LOGS/gen.$lbl.log" >&2; return 1; }
  grep -E "gen_(sha256|seconds)" "$LOGS/gen.$lbl.log" | cut -f5-6 | sed 's/^/  /'
}

# the abandon limit of a limited run: ABANDON x the process wall time of the generous run of the same text (run 1)
generous_wall() { rowval "$1.gen" 1 process_wall_s; }

# hash_outputs LABEL RUN PREFIX: sha256 of the files, rows sha256_<ext>; sets HASH_EQUAL (yes / NO / none) against the rows of the generous run of the same text
hash_outputs() {
  local L=$1 R=$2 O=$3 ext g="${1%%.*}.gen"
  for ext in pk pk.meta sa lcp ovf cld stats; do ( sha256sum "$O.$ext" 2>/dev/null | cut -d' ' -f1 > "$O.$ext.sha" ) & done; wait
  HASH_EQUAL=yes; local any=0
  for ext in pk pk.meta sa lcp ovf cld stats; do
    local got key ref; got=$(cat "$O.$ext.sha"); key=sha256_${ext//./_}; rm -f "$O.$ext.sha"
    addrow "$L" "$R" "$key" "$got" "" "SHA-256 of ${O##*/}.$ext"
    ref=$(rowval "$g" 1 "$key")
    if [ -n "$ref" ]; then any=1; [ "$ref" = "$got" ] || HASH_EQUAL=NO; fi
  done
  [ $any = 1 ] || HASH_EQUAL=none
}

run_cell() { # n cell-kind(gen|mNNN|tK) run
  local n=$1 kind=$2 R=$3 lbl L file O SCR PIDF SERIES
  lbl=$(lbl_n "$n"); file=$TEXTDIR/genome_$lbl.txt
  case $kind in gen) L=$lbl.gen ;; t*) L=$lbl.gen.$kind ;; m*) L=$lbl.$kind ;; esac
  if grep_done "$L" "$R"; then echo "$L run $R already done"; return 0; fi
  [ -f "$file" ] || do_gen "$n" || return $?
  local threads=$THREADS limit_check=
  case $kind in t*) threads=${kind#t} ;; esac
  # a limited cell must run in a container that has exactly that limit
  local watchdog=$HARD_WATCHDOG
  if [ "${kind#m}" != "$kind" ]; then
    local want=$(( ${kind#m} * 1048576 )) got; got=$(cg_limit_bytes)
    if [ "$got" = max ] || [ "$got" = none ] || ! awk -v g="$got" -v w="$want" 'BEGIN { d = g - w; if (d < 0) d = -d; exit !(d <= 0.05 * w) }'; then
      echo "REFUSED: this cell needs a container started with --memory ${kind#m}m (limit seen here: $got). Use e8-docker.sh sweep, or docker run --memory ${kind#m}m --memory-swap ${kind#m}m" >&2; return 5
    fi
    local gw; gw=$(generous_wall "$lbl"); if [ -z "$gw" ]; then echo "REFUSED: no generous run of $lbl in the rows file: run the generous step first (the abandon limit is $ABANDON x its time)" >&2; return 6; fi
    watchdog=$(awk -v g="$gw" -v f="$ABANDON" 'BEGIN { printf "%.0f", g * f }')
  fi
  # restart safety: whatever an interrupted attempt left behind goes
  purge "$L" "$R"
  O=$WORK/o.$L.$R; SCR=$WORK/scratch.$L.$R; PIDF=$WORK/pid.$L.$R; SERIES=$LOGS/series.$L.$R.tsv
  rm -rf "$SCR" "$O".* "$PIDF" "$PIDF".watchdog 2>/dev/null
  # disk
  local avail need; avail=$(df_avail_bytes "$DATA"); need=$(need_run_bytes "$n")
  addrow "$L" "$R" disk_free_at_start "$avail" bytes "free space of the data volume before the run; need $need (n x 21.5 + margin)"
  if [ "$avail" -lt "$need" ]; then addrow "$L" "$R" refused_disk 1 "" "free $avail < need $need"; echo "REFUSED: $(human "$need") needed for $L, $(human "$avail") free (remove other texts: rm $TEXTDIR/genome_n*m.txt, they are regenerated by gen)" >&2; return 4; fi
  local cold; cold=$(drop_caches)
  local img=${E8_IMAGE_COMMIT:-unknown} cgl swl; cgl=$(cg_limit_bytes); swl=$(cg_swap_limit_bytes); cg_detect
  addrow "$L" "$R" settings "n=$n kind=$kind threads=$threads gc=$GC scale_n=$SCALE_N cgroup=$CG_VERSION memory_limit=$cgl swap_allowance=$swl cold_start=$cold watchdog_s=$watchdog arch=$(uname -m) kernel=$(uname -r) image_commit=$img extra=${E8_EXTRA:-}" "" "run settings; cold_start=no means the page cache was not dropped before the start (warm start)"
  addrow "$L" "$R" cold_start "$cold" "" "page cache dropped before the start (needs a privileged container)"
  echo "== $L run $R: n=$n threads=$threads; ESTIMATE $(est_minutes "$n") min on the prototype machine (4 vCPU x86, virtio disk; E7), abandon after $((watchdog / 60)) min"
  python3 -I "$HERE/sample_e8.py" "$PIDF" "$ROWS" "$L" "$R" "$SCR" "$O" "$file" "$DATA" "$SERIES" &
  local SP=$!; CHILD_PIDS="$SP"
  sleep 0.2
  local S0 S1; S0=$(date +%s.%N)
  ( export DOTNET_GCHeapHardLimit=$GC; exec dotnet "$DLL" e6 --input "$file" --out "$O" --scratch "$SCR" --threads "$threads" --scale-n "$SCALE_N" --case "$L" --run "$R" --rows "$ROWS.prog.$L.$R" ${E8_EXTRA:-} > "$LOGS/e6.$L.$R.log" 2>&1 ) &
  local PID=$!; echo "$PID" > "$PIDF"; CHILD_PIDS="$CHILD_PIDS $PID"
  ( exec > /dev/null 2>&1; for ((w = 0; w < watchdog; w++)); do kill -0 $PID 2>/dev/null || exit 0; sleep 1; done; touch "$PIDF.watchdog"; kill -9 $PID 2>/dev/null ) &
  local WD=$!; CHILD_PIDS="$CHILD_PIDS $WD"
  wait $PID; local RC=$?
  S1=$(date +%s.%N)
  wait $WD 2>/dev/null
  wait $SP
  CHILD_PIDS=""
  [ -f "$ROWS.prog.$L.$R" ] && { sed 's/^R\tE6\t/R\tE8\t/' "$ROWS.prog.$L.$R" >> "$ROWS"; rm -f "$ROWS.prog.$L.$R"; }
  local wall; wall=$(echo "$S1 - $S0" | bc -l)
  addrow "$L" "$R" process_wall_s "$(printf '%.2f' "$wall")" s "from the launch of dotnet to its exit, including the runtime start (the program rows give the stages alone)"
  addrow "$L" "$R" exit_code "$RC" "" "program exit code (137 = killed by the watchdog or by the OOM killer)"
  local KILLED=no; [ -e "$PIDF.watchdog" ] && KILLED=watchdog
  [ "$RC" = 137 ] && [ "$KILLED" = no ] && KILLED=oom-or-signal
  addrow "$L" "$R" killed "$KILLED" "" "no | watchdog (more than the abandon limit: reported as 'does not finish in practice at this memory') | oom-or-signal"
  addrow "$L" "$R" threads "$threads" n "threads of the e6 command"
  HASH_EQUAL=none
  if [ "$RC" = 0 ]; then
    if [ "$kind" = gen ] && [ "$VERIFY" = 1 ]; then
      echo "   verifying the outputs of $L (ESTIMATE $(est_verify_minutes "$n") min with the files cached, more when they are not)"
      local V0; V0=$(date +%s.%N)
      dotnet "$DLL" e8-verify --prefix "$O" --input "$file" --threads "$threads" --lce-cap "$LCE_CAP" --scratch "$WORK" --sums "$SUMS/$L.$R.sums" --case "$L" --run "$R" --rows "$ROWS" > "$LOGS/verify.$L.$R.out" 2> "$LOGS/verify.$L.$R.log"
      local VRC=$?
      addrow "$L" "$R" verify_exit_code "$VRC" "" "e8-verify exit code (0 = every check passed)"
      addrow "$L" "$R" verify_process_wall_s "$(printf '%.2f' "$(echo "$(date +%s.%N) - $V0" | bc -l)")" s "wall time of the verifier process"
      HASH_EQUAL=none
    else
      hash_outputs "$L" "$R" "$O"
      addrow "$L" "$R" identity_vs_generous "$HASH_EQUAL" "" "the seven output files have the same SHA-256 as those of the generous run of the same text (none = no generous run to compare with, or this is the generous run)"
      if [ "$HASH_EQUAL" = NO ] && [ "$VERIFY" = 1 ]; then
        echo "   outputs DIFFER from the generous run: running the verifier on them"
        dotnet "$DLL" e8-verify --prefix "$O" --input "$file" --threads "$threads" --lce-cap "$LCE_CAP" --scratch "$WORK" --case "$L" --run "$R" --rows "$ROWS" > "$LOGS/verify.$L.$R.out" 2> "$LOGS/verify.$L.$R.log"
      fi
    fi
  fi
  [ "$KEEP" = 1 ] || rm -rf "$O".* 2>/dev/null
  addrow "$L" "$R" scratch_left "$( [ -e "$SCR" ] && ls -A "$SCR" | wc -l || echo 0)" "" "entries of the scratch directory after the run (a killed run leaves its scratch; it is removed here)"
  rm -rf "$SCR" "$PIDF" "$PIDF.watchdog"
  addrow "$L" "$R" run_done yes "" "all rows of this run are written"
  echo "   $L run $R: rc=$RC killed=$KILLED $(grep -P "^R\tE8\t$L\t$R\t(total_s|ext_rssanon_peak|eff_resident_ph3|identity_vs_generous|verify_all)\t" "$ROWS" | awk -F'\t' '{printf "%s=%s ", $5, $6}')"
  [ "$RC" = 0 ] || return 7        # the run did not finish (abandoned by the watchdog, killed, failed): the caller of a sweep stops after two such cells in a row
  return 0
}

parse_options() {
  while [ $# -gt 0 ]; do
    case $1 in
      --runs) RUNS=$2; shift 2 ;;
      --threads) THREADS=$2; shift 2 ;;
      --keep) KEEP=1; shift ;;
      --no-verify) VERIFY=0; shift ;;
      --lce-cap) LCE_CAP=$2; shift 2 ;;
      --size) SIZE_ARG=$(size_n "$2"); shift 2 ;;
      --memory-mb) MEMORY_MB=$2; shift 2 ;;
      *) echo "unknown option $1" >&2; exit 2 ;;
    esac
  done
}

status_of() { if grep_done "$1" "$2"; then echo done; else echo pending; fi; }

plan() {
  say "E8 protocol: cheap first. Times are ESTIMATES (250 ns per symbol on the prototype machine: 4 vCPU x86, virtio disk; E7 measured 205-300 ns up to 500 Mbp); on this machine they can be 1-3x longer, or more with a small page cache."
  say "Every step is restart safe: run it again after an interruption (the run in progress is purged and repeated; finished runs are skipped)."
  say
  local A B C la lb lc
  A=$SIZE_A; B=$SIZE_B; C=$SIZE_C; la=$(lbl_n "$A"); lb=$(lbl_n "$B"); lc=$(lbl_n "$C")
  printf '%-4s %-46s %-12s %-10s %s\n' step what estimate state "abandon rule / gate"
  printf '%-4s %-46s %-12s %-10s %s\n' a "$la generous memory, then the full verification" "$(est_minutes "$A") + $(est_verify_minutes "$A") min" "$(status_of "$la.gen" 1)" "none (the baseline of d)"
  printf '%-4s %-46s %-12s %-10s %s\n' b "$lb generous: crosses 2^31 = 2147483648" "$(est_minutes "$B") + $(est_verify_minutes "$B") min" "$(status_of "$lb.gen" 1)" "none; positions above 2^31 are unverified until this run"
  printf '%-4s %-46s %-12s %-10s %s\n' c "$lc generous" "$(est_minutes "$C") + $(est_verify_minutes "$C") min" "$(status_of "$lc.gen" 1)" "none; G2, G3 are judged here"
  printf '%-4s %-46s %-12s %-10s %s\n' e "$la with 1 and 2 threads (identity of the outputs)" "$(est_minutes "$A") x 1.8 + x 1.2" "$(status_of "$la.gen.t1" 1) / $(status_of "$la.gen.t2" 1)" "outputs must be SHA-256 identical"
  printf '%-4s %-46s %-12s %-10s %s\n' d "$la at decreasing --memory (see below)" "1-6 x step a per cell" "see below" "stopped after $ABANDON x the generous time: does not finish in practice at this memory"
  say
  say "step d: memory limits of the container in MiB, generous first (the cache floor of ADR-004 6b is about 1x (text + rank file) = $(awk -v n="$A" 'BEGIN { printf "%.0f", n * 0.5 / 1048576 }') MiB at this size for the text and the rank file alone, and PH4's scratch comes on top):"
  local m
  for m in ${E8_CELLS:-$(default_cells "$A")}; do printf '   %6s MiB  %s\n' "$m" "$(status_of "$la.m$m" 1)"; done
  say
  say "the gates (verbatim in E8-RUNBOOK.md): G1 verifier + identical SHA-256; G2 RssAnon <= 256 MB; G3 scratch <= 12.5 B/char and total disk <= 70 GB at 3 Gbp; G4 reported."
}

# default cells for a text of n symbols (MiB of container memory, largest first): f x (text + rank file) + 160 MiB for the process, the runtime and its pages, f = 16 (a limit that still holds everything), 6, 4, 3, 2.2, 1.6, 1.25, 1.0
# (the text and the rank file together are 0.5 B/char for 3-bit symbols: 477 MiB at 1 Gbp; the last cell is about the "1x (text + rank) + the process" of ADR-004 6b, with no room for PH4's scratch pages)
default_cells() {
  awk -v n="$1" 'BEGIN { base = n * 0.5 / 1048576; split("16 6 4 3 2.2 1.6 1.25 1.0", f, " ");
    for (i = 1; i <= 8; i++) { m = base * f[i] + 160; printf "%d ", int((m + 15) / 16) * 16 } }'
}

cmd=${1:-help}; [ $# -gt 0 ] && shift
case $cmd in
  plan) plan ;;
  gen) [ $# -gt 0 ] || { echo "usage: run_e8.sh gen <size>..." >&2; exit 2; }; for s in "$@"; do do_gen "$(size_n "$s")" || exit $?; done ;;
  clean) rm -rf "$WORK"/scratch.* "$WORK"/o.* "$WORK"/pid.* "$WORK"/*.verify.bits 2>/dev/null; echo "work directory cleaned" ;;
  run)
    step=${1:-}; [ $# -gt 0 ] && shift; parse_options "$@"
    case $step in
      a|b|c)
        n=$SIZE_A; [ "$step" = b ] && n=$SIZE_B; [ "$step" = c ] && n=$SIZE_C
        for r in $(seq 1 "$RUNS"); do run_cell "$n" gen "$r" || exit $?; done ;;
      e) for t in 1 2; do for r in $(seq 1 "$RUNS"); do run_cell "$SIZE_A" "t$t" "$r" || exit $?; done; done ;;
      d)
        n=${SIZE_ARG:-$SIZE_A}; [ -n "$MEMORY_MB" ] || { echo "usage: run_e8.sh run d [--size 1g] --memory-mb M" >&2; exit 2; }
        for r in $(seq 1 "$RUNS"); do run_cell "$n" "m$MEMORY_MB" "$r" || exit $?; done ;;
      *) echo "usage: run_e8.sh run a|b|c|d|e ..." >&2; exit 2 ;;
    esac ;;
  cells) default_cells "${1:-$SIZE_A}"; echo ;;
  *) sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//' ;;
esac
