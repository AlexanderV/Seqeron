#!/bin/bash
# Entry point of the hand-off package s1b-handoff (work package 1b.9f): the production chain of the SESA builder at up to 3x10^9 symbols and the test sets, on the owner's machine. It runs INSIDE the Docker image
# (the wrapper s1b-docker.sh starts it there) and natively on Linux. Nothing is sent anywhere: the results are files in the data directory; the owner returns form.json by hand.
#
#   s1b.sh preflight [--n N]       what this machine is, whether the run at N symbols fits, what to change (REFUSES with numbers)
#   s1b.sh plan [--n N]            the disk, RAM and time tables for 10^8, 10^9, 3x10^9 and N (estimates are marked)
#   s1b.sh gen    --n N [--kind genome|dna|rep|prot]       the text, streamed, deterministic, SHA-256 in <text>.sums (the text is a FILE of N bytes: the build reads it twice)
#   s1b.sh build  --n N [--kind K] [--resume auto|on|off]   one attempt of the production build (sesa resume build) under the independent sampler; run it again after a kill: with Resume on it continues
#   s1b.sh hash | verify | structure | queries --n N [--kind K]    the steps after the build, each on its own
#   s1b.sh chain  --n N [--kind K]   gen, build, hash, e8-verify --container, the Structure level (auto), the query set, the form, the checker
#   s1b.sh tests <set>...|all      the test sets of the owner's CI policy (esa-default esa-exhaustive esa-fuzz esa-checked emproto persistent); needs the SDK image
#   s1b.sh form | check            (re)write form.json / form.md from what the steps left; run the checker on it (PASS / FAIL per gate)
#   s1b.sh status | clean [--n N --kind K | --all]    what is done; remove a run's text, container and scratch
# options (anywhere): --data DIR  --kind K (genome)  --n N (decimal suffixes k m g)  --tag T (a second run of the same text)  --resume auto|on|off (auto: on from 10^9)  --threads T (4)  --fsync true|false (true)
#   --psi true|false (true)  --structure auto|yes|no (auto)  --queries yes|no (yes)  --rebuild  --kill-after SECONDS (rehearse a kill)  --dry-run  --force-disk  --interval-ms MS (sampler, 100)  --filter EXPR (tests)
# exit codes: 0 ok, 2 usage, 3 refused (pre-flight, disk), 4 a step failed. The chain always ends with the form and the checker, also after a failure.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=s1blib.sh
. "$HERE/s1blib.sh"
DATA=${S1B_DATA:-/data}
BENCH=${S1B_BENCH:-/opt/s1b/bin/bench/SuffixTree.Esa.Bench.dll}
EMPROTO=${S1B_EMPROTO:-/opt/s1b/bin/emproto/SuffixTree.EmProto.dll}
SRC=${S1B_SRC:-/src}
KIND=genome; N=; TAG=; RESUME=auto; THREADS=4; FSYNC=true; PSI=true; STRUCT=auto; QUERIES=yes; REBUILD=0; KILL_AFTER=; DRY=0; FORCE_DISK=0; INTERVAL=${S1B_INTERVAL_MS:-100}; FILTER=
HEARTBEAT=${S1B_HEARTBEAT:-30}; ALL=0; POS=()
CMD=${1:-help}; [ $# -gt 0 ] && shift
while [ $# -gt 0 ]; do
  case $1 in
    --data) DATA=$2; shift 2 ;;
    --kind) KIND=$2; shift 2 ;;
    --n) N=$(size_n "$2"); shift 2 ;;
    --tag) TAG=$2; shift 2 ;;
    --resume) RESUME=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --fsync) FSYNC=$2; shift 2 ;;
    --psi) PSI=$2; shift 2 ;;
    --structure) STRUCT=$2; shift 2 ;;
    --queries) QUERIES=$2; shift 2 ;;
    --interval-ms) INTERVAL=$2; shift 2 ;;
    --kill-after) KILL_AFTER=$2; shift 2 ;;
    --filter) FILTER=$2; shift 2 ;;
    --rebuild) REBUILD=1; shift ;;
    --dry-run) DRY=1; shift ;;
    --force-disk) FORCE_DISK=1; shift ;;
    --all) ALL=1; shift ;;
    -h|--help) CMD=help; shift ;;
    --*) die "unknown option $1 (s1b.sh help)" ;;
    *) POS+=("$1"); shift ;;
  esac
done
case $KIND in genome|dna|rep|prot) ;; *) die "--kind must be genome, dna, rep or prot" ;; esac
export S1B_DATA=$DATA S1B_BENCH=$BENCH S1B_EMPROTO=$EMPROTO S1B_SRC=$SRC

setup_run() { # sets LBL RUN WORK TEXT after N and KIND are known
  [ -n "$N" ] || die "--n is required for this command"
  LBL=$(run_label "$KIND" "$N" "$TAG"); RUN=$DATA/runs/$LBL; WORK=$DATA/work/$LBL; TEXT=$DATA/text/${KIND}_n$N.txt
  case $RESUME in
    auto) if [ "$N" -ge 1000000000 ]; then RESUME_BOOL=true; else RESUME_BOOL=false; fi ;;
    on|true) RESUME_BOOL=true ;;
    off|false) RESUME_BOOL=false ;;
    *) die "--resume must be auto, on or off" ;;
  esac
}
run_or_echo() { if [ "$DRY" = 1 ]; then echo "+ $*"; return 0; fi; "$@"; }
step() { echo; echo "==== $* ($(date -u +%H:%M:%SZ))"; }

write_run_json() {
  mkdir -p "$RUN"
  [ -f "$RUN/run.json" ] && return 0
  printf '{"kind": "%s", "n": %s, "threads": %s, "psi": %s, "resume": %s, "fsync": %s, "tag": "%s"}\n' "$KIND" "$N" "$THREADS" "$PSI" "$RESUME_BOOL" "$FSYNC" "$TAG" > "$RUN/run.json"
}

# ------------------------------------------------------------------------------------------------- plan
cmd_plan() {
  local n kinds
  echo "S1B plan. MEASURED = taken on the 4-vCPU x86-64 Linux box of the orchestrator at 10^7 / 10^8 (the ADR-004 subsection 1b.9f); ESTIMATE = scaled from it, not measured."
  echo
  echo "DISK at the peak of one run (the text + the container, reserved from the start + the scratch at its peak; the scratch falls to about 9 bytes per symbol once the suffix array is complete)"
  printf '  %-16s %-10s %-12s %-12s %-12s %-12s\n' "n" "kind" "text" "container" "scratch peak" "total peak"
  for n in 100000000 1000000000 3000000000 ${N:+$N}; do
    for kinds in genome dna; do
      printf '  %-16s %-10s %-12s %-12s %-12s %-12s\n' "$n" "$kinds" "$(human "$n")" "$(human "$(awk -v n="$n" -v c="$(container_bpc "$kinds")" 'BEGIN { printf "%.0f", n * c }')")" "$(human $((n * S1B_SCRATCH_BPC)))" "$(human "$(need_bytes "$n" "$kinds")")"
    done
  done
  echo "  (container bytes per symbol, MEASURED: dna 13.25, genome-like 13.49; protein-like and repetitive: see the ADR table. Peak total = n x (1 + container + 13) x 1.02 + 2 GB, ESTIMATE for n above 10^8)"
  echo
  echo "RAM: the gate is RssAnon <= $S1B_GATE_MB MB at every n (default parameters, four threads). MEASURED peaks: 262 MB (dna 10^8), 296 MB (genome 10^8). The page cache is elastic; 8 GB for the VM is plenty."
  echo "     Structure level of the reader (n <= $S1B_READER_LIMIT only): n / 8 bytes of bitmap + 34 MB (MEASURED 42 MB at 10^7, 109 MB at 10^8; ESTIMATE 0.45 GB at 2x10^9)."
  echo
  echo "TIME of the chain on this box (MEASURED at 10^7 and 10^8, ESTIMATE above; the Mac is 1.5-2x faster per the E8 prototype run: 123-137 ns per symbol there against 205-300 on the box):"
  echo "  see the table of ADR-004 1b.9f and the RUNBOOK section 6; the step durations are printed by 'chain' as it goes."
  echo
  echo "STEP RULES: --resume auto = on from 10^9 symbols (a kill after 40 minutes costs the stages that are done, not the run); --structure auto = on up to 3x10^8, off up to the reader's limit (a timed optional step, about one hour at 3x10^9 if it could run), on above it (it records the refusal: the reader takes at most $S1B_READER_LIMIT symbols)."
}

# ------------------------------------------------------------------------------------------------- gen
cmd_gen() {
  setup_run
  local sums=$TEXT.sums
  if [ "$DRY" = 1 ]; then
    if [ "$KIND" = genome ]; then echo "+ dotnet $EMPROTO e8-gen --n $N --out $TEXT --force"; else echo "+ python3 s1b_gen.py --kind $KIND --n $N --out $TEXT --force"; fi
    return 0
  fi
  mkdir -p "$DATA/text"
  if [ -f "$TEXT" ] && [ -f "$sums" ] && [ "$(fsize "$TEXT")" = "$N" ] && grep -q "^n=$N\$" "$sums"; then
    local want have; want=$(sed -n 's/^sha256=//p' "$sums"); have=$(sha256_of "$TEXT")
    if [ "$want" = "$have" ]; then echo "text ${KIND}_n$N exists, SHA-256 $have matches its sums file: reused"; return 0; fi
    echo "text ${KIND}_n$N exists but its SHA-256 differs from the sums file: generated again" >&2
  fi
  local avail; avail=$(free_bytes_of "$DATA")
  if [ "$avail" -lt $((N + 1000000000)) ]; then echo "REFUSED: $(human $((N + 1000000000))) needed for the text, $(human "$avail") free" >&2; return 3; fi
  local t0; t0=$(now_s)
  echo "generating $KIND, $N symbols, into $TEXT (a file of $N bytes: the build reads it twice, so it cannot be streamed into the builder)"
  if [ "$KIND" = genome ]; then
    dotnet "$EMPROTO" e8-gen --n "$N" --out "$TEXT" --force --case "gen.$KIND" > "$DATA/gen.$KIND.n$N.log" 2>&1 || { echo "generation failed:" >&2; tail -5 "$DATA/gen.$KIND.n$N.log" >&2; return 4; }
  else
    python3 -I "$HERE/s1b_gen.py" --kind "$KIND" --n "$N" --out "$TEXT" --force > "$DATA/gen.$KIND.n$N.log" 2>&1 || { echo "generation failed:" >&2; tail -5 "$DATA/gen.$KIND.n$N.log" >&2; return 4; }
  fi
  echo "  done in $(( $(now_s) - t0 )) s; $(grep -E '^sha256=' "$sums")"
}

# ------------------------------------------------------------------------------------------------- build
heartbeat() { # pid start
  local pid=$1 start=$2 last=0 now el m rss disk
  while kill -0 "$pid" 2> /dev/null; do
    sleep 1
    now=$(now_s); el=$((now - start))
    if [ $((el - last)) -ge "$HEARTBEAT" ]; then
      last=$el
      m=$(last_marker "$RUN/build.$K.out"); rss=$(tail -1 "$RUN/build.$K.series.tsv" 2> /dev/null | cut -f2); disk=$(tail -1 "$RUN/build.$K.series.tsv" 2> /dev/null | awk -F'\t' '{ printf "%.1f", ($6 + $7 + $8) / 1e9 }')
      printf '  [%5d s] last marker: %-14s RssAnon %s MB, work directory %s GB\n' "$el" "${m:-none yet}" "$(awk -v k="${rss:-0}" 'BEGIN { printf "%.0f", k * 1024 / 1e6 }')" "${disk:-?}"
    fi
  done
}

cmd_build() {
  setup_run
  if [ "$DRY" = 1 ]; then
    step "build attempt of $LBL: SesaBuilder.Build, $THREADS threads, PSI $PSI, fsync $FSYNC, Resume $RESUME_BOOL"
    echo "+ dotnet $BENCH sesa resume build --input $TEXT --work $WORK --threads $THREADS --psi $PSI --fsync $FSYNC --resume $RESUME_BOOL --keep true --hash false --markers true --label $LBL   (under s1b_sampler.py)"
    return 0
  fi
  write_run_json
  mkdir -p "$WORK" "$RUN"
  if [ "$REBUILD" = 1 ]; then rm -rf "$WORK" "$RUN/build."* "$RUN/container.sha256" "$RUN/verify."* "$RUN/structure."* "$RUN/queries."* "$RUN/form."*; mkdir -p "$WORK"; fi
  if [ -f "$RUN/container.sha256" ] && [ -f "$WORK/out.sesa" ]; then echo "the container of $LBL is built ($(human "$(fsize "$WORK/out.sesa")")): nothing to do (--rebuild starts again)"; return 0; fi
  [ -f "$TEXT" ] || die "the text $TEXT is missing: s1b.sh gen --n $N --kind $KIND first"
  # the disk: what this attempt still needs (a resumed attempt needs less: the stages that are done keep their files)
  local avail need existing; avail=$(free_bytes_of "$DATA"); need=$(need_bytes "$N" "$KIND"); need=$((need - N))
  existing=$(du -sk "$WORK" 2> /dev/null | cut -f1); existing=$(( ${existing:-0} * 1024 ))
  if [ "$FORCE_DISK" = 0 ] && [ "$avail" -lt $((need - existing)) ]; then
    echo "REFUSED: this attempt needs about $(human $((need - existing))) more and $(human "$avail") are free (s1b.sh preflight --n $N shows the table; --force-disk goes on at your risk: a full disk is an IOException, nothing is corrupted, and with Resume on the next call continues)" >&2; return 3
  fi
  K=$(( $(ls "$RUN"/build.*.out 2> /dev/null | wc -l) + 1 ))
  local state="no earlier attempt"
  [ -d "$WORK/out.sesa.tmp.d" ] && state="an earlier attempt left $(du -sh "$WORK/out.sesa.tmp.d" | cut -f1) of scratch"
  step "build attempt $K of $LBL: SesaBuilder.Build, $THREADS threads, PSI $PSI, fsync $FSYNC, Resume $RESUME_BOOL ($state)"
  local cmd=(dotnet "$BENCH" sesa resume build --input "$TEXT" --work "$WORK" --threads "$THREADS" --psi "$PSI" --fsync "$FSYNC" --resume "$RESUME_BOOL" --keep true --hash false --markers true --label "$LBL")
  echo "  ${cmd[*]##*/}" | sed "s#$DATA#<data>#g"
  printf '%s\n' "${cmd[*]}" | sed "s#$DATA#<data>#g" > "$RUN/build.$K.cmd"
  local start; start=$(now_s)
  "${cmd[@]}" > "$RUN/build.$K.out" 2> "$RUN/build.$K.err" &
  local pid=$! killer=
  python3 -I "$HERE/s1b_sampler.py" --pid "$pid" --out "$RUN/build.$K.series.tsv" --summary "$RUN/build.$K.sampler.json" --interval-ms "$INTERVAL" --disk-interval-ms 1000 --dir "$WORK" &
  local spid=$!
  if [ -n "$KILL_AFTER" ]; then ( sleep "$KILL_AFTER"; kill -9 "$pid" 2> /dev/null; echo "  (--kill-after: SIGKILL sent to the build after $KILL_AFTER s)" ) & killer=$!; fi
  trap 'kill -9 "$pid" 2>/dev/null; kill "$spid" 2>/dev/null; echo "interrupted: the attempt is recorded; run the same command again" >&2; exit 130' INT TERM
  [ "$HEARTBEAT" -gt 0 ] && heartbeat "$pid" "$start"
  wait "$pid"; local rc=$?
  wait "$spid" 2> /dev/null
  [ -n "$killer" ] && { kill "$killer" 2> /dev/null; wait "$killer" 2> /dev/null; }
  trap - INT TERM
  echo "$rc" > "$RUN/build.$K.rc"; echo $(( $(now_s) - start )) > "$RUN/build.$K.wall_s"
  if [ "$rc" = 0 ]; then
    grep -a '^M	' "$RUN/build.$K.out" | cut -f2,3 | tr '\n' ' ' | fold -w 150 | sed 's/^/  markers: /'
    awk -F'\t' '$1 == "R" && $3 ~ /^(total_s|ph1_s|ph23_s|ph45_s|bwt_s|psi_s|file_size|scratch_peak|rssanon_peak_mb|skipped)$/ { printf "  %-18s %s %s\n", $3, $4, $5 }' "$RUN/build.$K.out"
    echo "  independent sampler: RssAnon peak $(python3 -I -c "import json;print('%.1f MB' % json.load(open('$RUN/build.$K.sampler.json'))['rssanon_peak_mb'])") (gate $S1B_GATE_MB MB)"
    return 0
  fi
  echo "BUILD ATTEMPT $K ENDED WITH EXIT CODE $rc" >&2
  case $rc in
    137) echo "  exit code 137 = SIGKILL: the process was killed (an out-of-memory kill by the Docker VM, docker stop, or --kill-after). Docker Desktop: more memory (Settings > Resources) if it was the OOM killer." >&2 ;;
    3) echo "  exit code 3: an IOException, usually a full disk ($(head -c 300 "$RUN/build.$K.err"; grep -a 'io_exception' "$RUN/build.$K.out" | cut -f4 | head -c 300))" >&2 ;;
  esac
  tail -3 "$RUN/build.$K.err" 2> /dev/null | sed 's/^/  /' >&2
  if [ "$RESUME_BOOL" = true ]; then echo "  Resume is ON: run the same command again; it skips the stages that are done and verified (it checks the CRC-32C of what it keeps)." >&2
  else echo "  Resume is OFF: the next call starts the build again from the beginning (the scratch is removed first). Use --resume on for a long run." >&2; fi
  return 4
}

cmd_hash() {
  setup_run
  if [ "$DRY" = 1 ]; then step "SHA-256 of the container"; echo "+ sha256sum $WORK/out.sesa > $RUN/container.sha256"; return 0; fi
  [ -f "$WORK/out.sesa" ] || die "no container for $LBL"
  step "SHA-256 of the container"
  local t0 sha; t0=$(now_s); sha=$(sha256_of "$WORK/out.sesa")
  printf '%s  %s\n' "$sha" "$(fsize "$WORK/out.sesa")" > "$RUN/container.sha256"
  echo "  $sha  ($(human "$(fsize "$WORK/out.sesa")") in $(( $(now_s) - t0 )) s)"
}

cmd_verify() {
  setup_run
  step "e8-verify --container (eleven checks; no in-RAM reference; reads the sections in place)"
  if [ "$DRY" = 1 ]; then echo "+ dotnet $EMPROTO e8-verify --container $WORK/out.sesa --input $TEXT --threads $THREADS --scratch $WORK/vscratch --case $LBL"; return 0; fi
  [ -f "$WORK/out.sesa" ] || die "no container for $LBL"
  mkdir -p "$WORK/vscratch" "$RUN"
  local t0; t0=$(now_s)
  dotnet "$EMPROTO" e8-verify --container "$WORK/out.sesa" --input "$TEXT" --threads "$THREADS" --scratch "$WORK/vscratch" --case "$LBL" > "$RUN/verify.out" 2> "$RUN/verify.err"
  local rc=$?; echo "$rc" > "$RUN/verify.rc"; rm -rf "$WORK/vscratch"
  grep -E "verify_all|verify_wall_s" "$RUN/verify.out" | cut -f4,5,7 | sed 's/^/  /'
  echo "  exit code $rc in $(( $(now_s) - t0 )) s"
  [ "$rc" = 0 ] || { tail -5 "$RUN/verify.err" >&2; return 4; }
}

cmd_structure() {
  setup_run
  [ "$DRY" = 1 ] || [ -f "$WORK/out.sesa" ] || die "no container for $LBL"
  local note="at the reader's limit"
  is_reader_n "$N" && note="about $(awk -v n="$N" 'BEGIN { printf "%.0f", (n / 1e8) * 27 / 60 + 0.1 }') minutes at this size by the 10^8 figure (27 s), more when the SA does not fit the page cache: 2x10^9 about 20-60 min, ESTIMATE"
  is_reader_n "$N" || note="n is above $S1B_READER_LIMIT: the reader refuses (S2); this records the refusal after the checksum pass over the sections"
  step "Structure level of the reader ($note)"
  if [ "$DRY" = 1 ]; then echo "+ dotnet $BENCH sesa handoff open --container $WORK/out.sesa --level structure --sha true --label $LBL"; return 0; fi
  local t0; t0=$(now_s)
  dotnet "$BENCH" sesa handoff open --container "$WORK/out.sesa" --level structure --sha true --label "$LBL" > "$RUN/structure.out" 2> "$RUN/structure.err"
  local rc=$?; echo "$rc" > "$RUN/structure.rc"
  awk -F'\t' '$1 == "R" && $3 ~ /^open_(result|seconds|structure_s|sections_s|rssanon_peak_mb)$/ { printf "  %-22s %s %s\n", $3, $4, $6 }' "$RUN/structure.out" | cut -c1-200
  echo "  exit code $rc in $(( $(now_s) - t0 )) s (0 pass, 3 refused as corrupt, 4 not applicable above the reader's limit)"
  [ "$rc" = 0 ] || [ "$rc" = 4 ] || return 4
}

cmd_queries() {
  setup_run
  [ "$DRY" = 1 ] || [ -f "$WORK/out.sesa" ] || die "no container for $LBL"
  step "query set (count, exists, locate) against the oracle: queries cut from the text and random absent ones, one rolling-hash pass over the text"
  if [ "$DRY" = 1 ]; then echo "+ dotnet $BENCH sesa handoff queries --container $WORK/out.sesa --input $TEXT --label $LBL"; return 0; fi
  local t0; t0=$(now_s)
  dotnet "$BENCH" sesa handoff queries --container "$WORK/out.sesa" --input "$TEXT" --label "$LBL" > "$RUN/queries.out" 2> "$RUN/queries.err"
  local rc=$?; echo "$rc" > "$RUN/queries.rc"
  awk -F'\t' '$1 == "R" && $3 ~ /^q_(total|oracle_scan_s|oracle_hits_total|sa_s|sa_mismatches|engine_status|engine_s|engine_mismatches|result|mismatch)$/ { printf "  %-22s %s %s\n", $3, $4, substr($6, 1, 110) }' "$RUN/queries.out"
  echo "  exit code $rc in $(( $(now_s) - t0 )) s"
  [ "$rc" = 0 ] || { tail -3 "$RUN/queries.err" >&2; return 4; }
}

# ------------------------------------------------------------------------------------------------- tests
cmd_tests() {
  local sets=("${POS[@]}")
  [ "${#sets[@]}" -gt 0 ] || die "tests: name the sets (esa-default esa-exhaustive esa-fuzz esa-checked emproto persistent) or all"
  [ "${sets[0]}" = all ] && sets=(esa-default esa-exhaustive esa-fuzz esa-checked emproto persistent)
  command -v dotnet > /dev/null && dotnet --list-sdks 2> /dev/null | grep -q . || die "no .NET SDK here: the test sets need the SDK image (s1b-docker.sh image builds it)"
  [ -f "$SRC/Seqeron.sln" ] || die "no repository checkout at $SRC (S1B_SRC): the test sets run from the sources"
  mkdir -p "$DATA/tests"
  local flags=${S1B_MSBUILD_FLAGS:--p:NuGetAudit=false} T=$SRC/tests/SuffixTree rc=0 set dir args envs
  for set in "${sets[@]}"; do
    envs=(); args=()
    case $set in
      esa-default) dir=$T/SuffixTree.Esa.Tests; args=(test --nologo) ;;
      esa-exhaustive) dir=$T/SuffixTree.Esa.Tests; args=(test -c Release --nologo --filter "Category=Exhaustive") ;;
      esa-fuzz) dir=$T/SuffixTree.Esa.Tests; args=(test -c Release --nologo --filter "Category=Fuzz") ;;
      esa-checked) dir=$T/SuffixTree.Esa.Tests; args=(test --nologo -p:EsaBounds=true); envs=(ESA_EXPECT_BOUNDS=true) ;;
      emproto) dir=$T/SuffixTree.EmProto.Tests; args=(test -c Release --nologo) ;;
      persistent) dir=$T/SuffixTree.Persistent.Tests; args=(test --nologo) ;;
      *) die "unknown test set '$set'" ;;
    esac
    if [ -n "$FILTER" ]; then args=("${args[@]/--filter=*/}"); args+=(--filter "$FILTER"); fi
    # shellcheck disable=SC2206
    args+=($flags)
    step "test set $set"
    echo "  (cd tests/SuffixTree/$(basename "$dir") && ${envs[*]:-} dotnet ${args[*]})"
    [ "$DRY" = 1 ] && continue
    local t0; t0=$(now_s)
    ( cd "$dir" && env "${envs[@]+"${envs[@]}"}" dotnet "${args[@]}" ) > "$DATA/tests/$set.log" 2>&1
    local trc=$? secs=$(( $(now_s) - t0 ))
    local line; line=$(grep -aE '(Passed|Failed)!  - Failed:' "$DATA/tests/$set.log" | tail -1)
    local passed failed skipped
    passed=$(echo "$line" | sed -n 's/.*Passed: *\([0-9]*\).*/\1/p'); failed=$(echo "$line" | sed -n 's/.*Failed: *\([0-9]*\),.*/\1/p'); skipped=$(echo "$line" | sed -n 's/.*Skipped: *\([0-9]*\).*/\1/p')
    local status=pass; { [ "$trc" != 0 ] || [ -z "$line" ] || [ "${failed:-1}" != 0 ]; } && status=fail
    printf '{"status": "%s", "passed": %s, "failed": %s, "skipped": %s, "seconds": %s, "exit_code": %s, "configuration": "%s", "filter": "%s"}\n' "$status" "${passed:-0}" "${failed:-0}" "${skipped:-0}" "$secs" "$trc" "$([ "${args[*]}" != "${args[*]/Release/}" ] && echo Release || echo Debug)" "$(json_escape "$FILTER")" > "$DATA/tests/$set.json"
    echo "  $set: $status, ${line:-no summary line (see $DATA/tests/$set.log)} in $secs s"
    [ "$status" = pass ] || { tail -15 "$DATA/tests/$set.log" | sed 's/^/    /'; rc=4; }
  done
  return $rc
}

# ------------------------------------------------------------------------------------------------- form, check, status, clean
cmd_form() { python3 -I "$HERE/s1b_form.py" "$DATA" --commit "${S1B_COMMIT:-unknown}"; }
cmd_check() { python3 -I "$HERE/s1b_check.py" "$DATA/form.json" "$@"; }

cmd_status() {
  echo "data directory: $DATA"
  local d
  for d in "$DATA"/runs/*/; do
    [ -f "$d/run.json" ] || continue
    printf '  %-34s attempts %s, container %s, verify %s, structure %s, queries %s\n' "$(basename "$d")" "$(ls "$d"/build.*.out 2> /dev/null | wc -l)" "$([ -f "$d/container.sha256" ] && echo done || echo -)" "$([ -f "$d/verify.rc" ] && echo "rc $(cat "$d/verify.rc")" || echo -)" "$([ -f "$d/structure.rc" ] && echo "rc $(cat "$d/structure.rc")" || echo -)" "$([ -f "$d/queries.rc" ] && echo "rc $(cat "$d/queries.rc")" || echo -)"
  done
  ls "$DATA/tests/"*.json 2> /dev/null | sed 's#.*/#  tests: #; s#\.json##'
  echo "free: $(human "$(free_bytes_of "$DATA")")"
}

cmd_clean() {
  if [ "$ALL" = 1 ]; then rm -rf "$DATA/work" "$DATA/text"; echo "removed all texts and containers (runs/ and form.json stay)"; return 0; fi
  setup_run
  rm -rf "$WORK"; rm -f "$TEXT" "$TEXT.sums"
  echo "removed the container, scratch and text of $LBL (the rows, logs and hashes in runs/$LBL stay)"
}

# ------------------------------------------------------------------------------------------------- chain
cmd_chain() {
  setup_run
  step "chain $LBL: n = $N, kind $KIND, Resume $RESUME_BOOL"
  if [ "$DRY" = 1 ]; then
    echo "DRY RUN: nothing heavy runs; the pre-flight below is real, the steps are listed."
    "$HERE/preflight.sh" --n "$N" --kind "$KIND" --quick --no-form || true
    local s; for s in gen build hash verify structure queries; do "cmd_$s" || true; done
    return 0
  fi
  mkdir -p "$DATA"; write_run_json
  "$HERE/preflight.sh" --n "$N" --kind "$KIND" --quick > "$DATA/preflight.last.txt" 2>&1; local prc=$?
  sed -n '/^RESULT/p;/\[REFUSE\]/p;/\[warn\]/p' "$DATA/preflight.last.txt"
  if [ "$prc" != 0 ] && [ "$FORCE_DISK" = 0 ]; then cat "$DATA/preflight.last.txt"; echo "chain refused: fix the [REFUSE] lines (s1b.sh preflight --n $N shows everything)" >&2; return 3; fi
  local rc=0
  cmd_gen || rc=$?
  [ $rc = 0 ] && { cmd_build || rc=$?; }
  [ $rc = 0 ] && { cmd_hash || rc=$?; }
  [ $rc = 0 ] && { cmd_verify || rc=$?; }
  local do_struct=0
  case $STRUCT in
    yes) do_struct=1 ;;
    no) do_struct=0 ;;
    *) if [ "$N" -le 300000000 ] || ! is_reader_n "$N"; then do_struct=1; fi ;;
  esac
  if [ $rc = 0 ] && [ $do_struct = 1 ]; then cmd_structure || rc=$?
  elif [ $rc = 0 ]; then echo; echo "Structure level not run here (--structure auto: only up to 3x10^8 symbols, or above the reader's limit). Run it on its own when you want it: s1b.sh structure --n $N --kind $KIND (about $(awk -v n="$N" 'BEGIN { printf "%.0f", (n / 1e8) * 27 / 60 + 0.1 }') min by the 10^8 figure, ESTIMATE)."; fi
  if [ $rc = 0 ] && [ "$QUERIES" = yes ]; then cmd_queries || rc=$?; fi
  step "form and checker"
  cmd_form
  cmd_check || { local crc=$?; [ $rc = 0 ] && rc=$crc; }
  return $rc
}

case $CMD in
  preflight) exec "$HERE/preflight.sh" ${N:+--n "$N"} --kind "$KIND" ;;
  plan) cmd_plan ;;
  gen) cmd_gen ;;
  build) cmd_build ;;
  hash) cmd_hash ;;
  verify) cmd_verify ;;
  structure) cmd_structure ;;
  queries) cmd_queries ;;
  chain) cmd_chain ;;
  tests) cmd_tests ;;
  form) cmd_form ;;
  check) cmd_check "${POS[@]+"${POS[@]}"}" ;;
  status) cmd_status ;;
  clean) cmd_clean ;;
  shell) exec /bin/bash ;;
  help|-h|--help) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//' ;;
  *) echo "unknown command '$CMD'" >&2; sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 2 ;;
esac
