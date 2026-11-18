#!/bin/bash
# E8 preflight (runs INSIDE the container; also runs on a plain Linux host). UNBUILT AND UNTESTED IN A CONTAINER, ON macOS AND ON arm64.
# usage: preflight.sh [--sizes "1g 2.2g 3g"] [--quick] [--no-rows]
#   --sizes   the text sizes the protocol will use (default: E8_SIZE_A / _B / _C = 1g 2.2g 3g); the disk need is computed for the largest
#   --quick   skip the 64 MB write / mmap / msync / fallocate test (the runner uses this before every run)
#   --cell    a cell of the limited-cache sweep (step d): the small memory limit of the container is the experiment, not a problem
# Prints [ ok ] / [warn] / [REFUSE] lines and, at the end, what to change in Docker Desktop. Exit code 0 = no refusal (warnings allowed), 3 = at least one refusal.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=e8lib.sh
. "$HERE/e8lib.sh"
DATA=${E8_DATA:-/data}
DLL=${E8_DLL:-/opt/e8/bin/SuffixTree.EmProto.dll}
ROWS=$DATA/rows/e8_rows.tsv
SIZES="${E8_SIZE_A:-1000000000} ${E8_SIZE_B:-2200000000} ${E8_SIZE_C:-3000000000}"
QUICK=0; WRITE_ROWS=1; CELL=0
while [ $# -gt 0 ]; do
  case $1 in
    --sizes) SIZES=""; for s in $2; do SIZES="$SIZES $(size_n "$s")"; done; shift 2 ;;
    --quick) QUICK=1; shift ;;
    --no-rows) WRITE_ROWS=0; shift ;;
    --cell) CELL=1; shift ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
REFUSALS=0; WARNINGS=0
ok()     { say "  [ ok ] $*"; }
warn()   { say "  [warn] $*"; WARNINGS=$((WARNINGS + 1)); }
refuse() { say "  [REFUSE] $*"; REFUSALS=$((REFUSALS + 1)); }
row()    { # metric value unit notes
  [ "$WRITE_ROWS" = 1 ] || return 0
  mkdir -p "$(dirname "$ROWS")" 2>/dev/null && printf 'R\tE8\tpreflight\t0\t%s\t%s\t%s\t%s\n' "$1" "$2" "${3:-}" "${4:-}" >> "$ROWS"
}
NEED_LINES=""
MAXN=0; for n in $SIZES; do [ "$n" -gt "$MAXN" ] && MAXN=$n; done

say "E8 preflight  ($(date -u +%Y-%m-%dT%H:%M:%SZ))   sizes: $SIZES   data directory: $DATA"
say
say "1. Processor architecture and emulation"
UNAME=$(uname -m); row uname_m "$UNAME" "" "uname -m"
IMGARCH=${E8_IMAGE_ARCH:-unknown}; row image_arch "$IMGARCH" "" "architecture the image was built for (TARGETARCH)"
CPUFMT=unknown
if grep -qE '^(CPU implementer|Features)' "$E8_PROC/cpuinfo" 2>/dev/null; then CPUFMT=arm; elif grep -q '^flags' "$E8_PROC/cpuinfo" 2>/dev/null; then CPUFMT=x86; fi
CPUMODEL=$(grep -m1 -E '^(model name|Model|Hardware|cpu model)' "$E8_PROC/cpuinfo" 2>/dev/null | sed 's/^[^:]*: *//')
say "     uname -m: $UNAME, image built for: $IMGARCH, /proc/cpuinfo format: $CPUFMT${CPUMODEL:+, model: $CPUMODEL}"
row cpuinfo_format "$CPUFMT" "" "x86 = a flags line, arm = CPU implementer / Features lines"
EMU=no
case "$UNAME:$CPUFMT" in x86_64:arm) EMU="yes: this is an x86-64 process on an ARM CPU (Rosetta or QEMU)" ;; aarch64:x86) EMU="yes: an ARM process on an x86 CPU" ;; esac
grep -qaE 'rosetta|qemu' "$E8_PROC/self/maps" 2>/dev/null && EMU="yes: an emulator is mapped into this process (/proc/self/maps)"
case "$IMGARCH:$UNAME" in amd64:aarch64|arm64:x86_64) EMU="yes: the image is for $IMGARCH but the kernel reports $UNAME" ;; esac
row emulation "$EMU" "" "refusal when the process runs under emulation"
if [ "$EMU" = no ]; then ok "native: the process runs on the architecture it was built for ($UNAME)"; else refuse "emulation detected ($EMU). Timings under emulation are meaningless. Build and run the image for the CPU of this Mac (Apple Silicon: linux/arm64; Intel: linux/amd64): e8-docker.sh does this by itself; never run the amd64 image on Apple Silicon."; fi
case $UNAME in aarch64|arm64) say "     note: arm64 timings are not comparable with the x86-64 timings of the prototype machine; correctness (identical outputs) is what carries over";; esac
say
say "2. Processors"
CORES=$(nproc 2>/dev/null || echo 1); CPULIM=$(cg_cpu_limit); row cores "$CORES" n "nproc"; row cpu_limit "$CPULIM" cpus "cgroup CPU quota (max = none)"
say "     cores visible: $CORES, CPU quota of the container: $CPULIM"
if [ "$CORES" -ge 4 ]; then ok "$CORES cores (the protocol uses 4 threads; E7: 1 / 2 / 4 threads took 39 / 29 / 25 s at 100 Mbp on the prototype machine)"; else warn "only $CORES core(s): runs use $CORES thread(s) and take longer; the thread-count identity step needs at least 2"; fi
if [ "$CPULIM" != max ] && [ "$CPULIM" != none ] && awk -v c="$CPULIM" 'BEGIN { exit !(c < 3.5) }'; then warn "the container is limited to $CPULIM CPUs (--cpus); the protocol assumes 4"; fi
say
say "3. Memory"
MEMT=$(mem_total_bytes); CGL=$(cg_limit_bytes); cg_detect; SWAPL=$(cg_swap_limit_bytes); SWAPT=$(swap_total_bytes)
row mem_total "$MEMT" bytes "MemTotal visible in the container (on Docker Desktop: the memory of the Linux VM)"; row cgroup_version "$CG_VERSION" "" "cgroup hierarchy seen by the container (v1 and v2 are both handled)"; row cgroup_dir "$CG_DIR" "" "directory read for the memory counters"
row cg_memory_limit "$CGL" bytes "memory limit of the container (max = none)"; row cg_swap_limit "$SWAPL" bytes "swap the container may add to its memory limit (0 = none)"; row swap_total "$SWAPT" bytes "SwapTotal of the VM"
say "     memory visible (VM): $(human "$MEMT"), cgroup $CG_VERSION limit of this container: $( [ "$CGL" = max ] && echo none || human "$CGL"), swap allowance: $SWAPL, VM swap: $(human "${SWAPT:-0}")"
EFF=$MEMT; if [ "$CGL" != max ] && [ "$CGL" != none ] && [ "$CGL" -lt "$MEMT" ]; then EFF=$CGL; fi
GENEROUS=$((6 * 1000 * 1000 * 1000))
if [ "$CG_VERSION" = none ]; then warn "no cgroup memory files are visible: the limited-memory cells (step d) cannot be checked or sampled"; fi
if [ "$CELL" = 1 ]; then ok "step d cell: the container is limited to $(human "$EFF") on purpose (the VM has $(human "$MEMT"))";
elif [ "$EFF" -ge "$GENEROUS" ]; then ok "$(human "$EFF") available to the container: enough for the 'generous' steps a, b, c (the text and the rank file together are about 1.5 GB at 3 Gbp; they and the scratch pages stay cached)";
else warn "only $(human "$EFF") available: the 'generous' runs of steps a-c would be cache-limited themselves, and every slowdown in step d is then relative to a slow baseline. Give Docker at least 8 GB."; fi
if [ "$CELL" = 0 ] && [ "$EFF" -lt $((3 * 1000 * 1000 * 1000)) ]; then refuse "less than 3 GB of memory: not even the generous 1 Gbp run is meaningful"; fi
case $SWAPL in
  0|none) : ;;
  *) if [ "${SWAPT:-0}" -gt 0 ]; then warn "swap is available to the container ($SWAPL): under a memory limit the kernel may swap the program's anonymous memory instead of dropping cache; start the container with --memory-swap equal to --memory (e8-docker.sh does)"; fi ;;
esac
say
say "4. Disk of the data volume"
mkdir -p "$DATA" 2>/dev/null
if [ ! -d "$DATA" ] || [ ! -w "$DATA" ]; then refuse "$DATA does not exist or is not writable"; else
  AVAIL=$(df_avail_bytes "$DATA"); SIZE=$(df_size_bytes "$DATA")
  PRESENT=0; [ -f "$DATA/text/genome_$(lbl_n "$MAXN").txt" ] && PRESENT=$(stat -c %s "$DATA/text/genome_$(lbl_n "$MAXN").txt" 2>/dev/null || echo 0)
  NEED=$(need_total_bytes "$MAXN"); NEED=$((NEED - PRESENT))
  row disk_avail "$AVAIL" bytes "free space of the data volume"; row disk_size "$SIZE" bytes "size of the file system holding the data volume (on Docker Desktop: the virtual disk)"; row disk_need "$NEED" bytes "need of the largest run: n x (1 + 12 + 9.5) + 2 GB + 3 % margin (ESTIMATES: scratch 12 B/char peak, outputs 9.25-9.50 B/char, ASCII text 1 B/char)"
  say "     free: $(human "$AVAIL") of $(human "$SIZE"); need for the largest run ($MAXN symbols): $(human "$NEED") (ASCII text 1 + scratch 12 + outputs 9.5 B/char, plus a margin; plus $(human $(( $(awk -v t="$SIZES" 'BEGIN{n=split(t,a," "); s=0; for(i=1;i<=n;i++) s+=a[i]; print s}') - MAXN ))) if the texts of the smaller sizes stay on the volume)"
  if [ "$AVAIL" -ge "$NEED" ]; then ok "the largest run fits ($(human "$AVAIL") free, $(human "$NEED") needed)"; else refuse "the largest run ($MAXN symbols) needs about $(human "$NEED") and only $(human "$AVAIL") are free. Enlarge the disk image of Docker Desktop (see the end of this report), or run smaller sizes (--sizes)."; fi
  if [ "$SIZE" -lt 100000000000 ]; then warn "the file system is $(human "$SIZE") big; for the 3 Gbp run the runbook asks for a disk image of at least 100 GB (the default of Docker Desktop is smaller)"; fi
fi
say
say "5. File system of the data volume"
FST=$(fs_type "$DATA"); row fs_type "$FST" "" "stat -f / mount type of the data directory"
say "     type: $FST"
if bad_fs "$FST" && [ "${E8_ALLOW_BADFS:-0}" != 1 ]; then refuse "$FST is a network or user-space file system: mmap with a large dirty set is slow or unreliable on it (a macOS bind mount shows up as virtiofs, grpcfuse or osxfs). Use a named Docker volume (e8-docker.sh does)."; else
  case "$FST" in
    tmpfs*|ramfs*) if [ "${E8_ALLOW_TMPFS:-0}" = 1 ]; then warn "tmpfs: the data lives in RAM (allowed by E8_ALLOW_TMPFS, for tests only)"; else refuse "tmpfs: the data volume is in RAM and the sizes of this protocol do not fit; use a named Docker volume"; fi ;;
    overlay*) warn "overlay: the data directory is in the container's writable layer, not in a volume; it works but is slower and is lost with the container. Use the named volume of e8-docker.sh." ;;
    *) ok "$FST" ;;
  esac
  if command -v mountpoint > /dev/null && ! mountpoint -q "$DATA" 2>/dev/null; then warn "$DATA is not a mount point of its own (a named volume is expected there)"; fi
fi
say
say "6. Write, mmap, msync and fallocate test on the data volume (64 MB)"
if [ "$QUICK" = 1 ]; then say "     skipped (--quick)"; else
  IOOUT=$(python3 -I - "$DATA/.e8-io-test" <<'PY'
import mmap, os, sys, time
path = sys.argv[1]
size = 64 << 20
res = {}
try:
    fd = os.open(path, os.O_RDWR | os.O_CREAT | os.O_TRUNC, 0o644)
    t = time.time()
    try:
        os.posix_fallocate(fd, 0, size)
        res['fallocate'] = 'ok'
    except OSError as e:
        res['fallocate'] = 'unsupported (%s)' % e.strerror
        os.ftruncate(fd, size)
    mm = mmap.mmap(fd, size)
    pattern = bytes(range(256)) * 4096
    for off in range(0, size, len(pattern)):
        mm[off:off + len(pattern)] = pattern
    mm.flush()
    res['mmap_write_msync_s'] = '%.3f' % (time.time() - t)
    mm.close()
    os.fsync(fd)
    t = time.time()
    with open(path, 'rb') as f:
        data = f.read()
    res['readback'] = 'ok' if data == pattern * (size // len(pattern)) else 'MISMATCH'
    res['write_mb_s'] = '%.0f' % (size / 1e6 / max(1e-6, float(res['mmap_write_msync_s'])))
    os.close(fd)
except Exception as e:  # noqa
    res['error'] = repr(e)
finally:
    try:
        os.unlink(path)
    except OSError:
        pass
print(' '.join('%s=%s' % kv for kv in res.items()))
PY
)
  say "     $IOOUT"
  row io_test "$IOOUT" "" "64 MB write + mmap + msync + fallocate + read-back on the data volume"
  case "$IOOUT" in
    *error=*|*MISMATCH*) refuse "the write / mmap / msync test failed: $IOOUT" ;;
    *) ok "mmap write, msync and read-back work"
       case "$IOOUT" in *"fallocate=unsupported"*) warn "fallocate is not supported here (the pipeline sets file lengths with SetLength, so this is informative)" ;; esac
       MBS=$(echo "$IOOUT" | sed -n 's/.*write_mb_s=\([0-9]*\).*/\1/p'); if [ -n "$MBS" ] && [ "$MBS" -lt 100 ]; then warn "sequential write through a mapping was $MBS MB/s (including msync): below 100 MB/s the disk, not the algorithm, will dominate the times"; fi ;;
  esac
fi
say
say "7. Dropping the page cache (cold-cache runs)"
if [ -w "$E8_PROC/sys/vm/drop_caches" ] && sync && echo 1 > "$E8_PROC/sys/vm/drop_caches" 2>/dev/null; then COLD=yes; ok "drop_caches works: every run can start cold"; else COLD=no; warn "drop_caches is not writable (the container is not privileged): runs start 'warm', the first reads may be served from the page cache of the VM; rows say cold_start=no. Start the container with --privileged for cold starts (e8-docker.sh --cold)."; fi
row cold_start_possible "$COLD" "" "drop_caches writable"
say
say "8. .NET and the program"
if [ ! -f "$DLL" ]; then refuse "$DLL not found (the program is not in this image)"; else
  command -v dotnet > /dev/null && { say "     $(dotnet --list-runtimes 2>/dev/null | grep Microsoft.NETCore.App | head -1)"; row dotnet_runtime "$(dotnet --list-runtimes 2>/dev/null | grep Microsoft.NETCore.App | head -1 | tr '\t' ' ')" "" "dotnet --list-runtimes"; }
  if [ "$WRITE_ROWS" = 1 ]; then ENVOUT=$(dotnet "$DLL" e8-env --case preflight --run 0 --rows "$ROWS" 2>&1) || ENVOUT="e8-env failed: $ENVOUT"; else ENVOUT=$(dotnet "$DLL" e8-env --case preflight --run 0 2>&1) || ENVOUT="e8-env failed: $ENVOUT"; fi
  echo "$ENVOUT" | awk -F'\t' '$1=="R" && NF >= 7 {printf "     %-22s %s\n", $5 ":", $6}' | grep -E "env_framework|env_process_arch|env_os_arch|env_cores|env_gc_available|env_page|env_lce" | sed 's/^/ /'
  if echo "$ENVOUT" | grep -q "env_lce_path	scalar"; then warn "no AVX2: the suffix comparison uses the scalar 64-bit path (same results, slower on tie-rich text; normal on arm64)"; fi
  if echo "$ENVOUT" | grep -q "env_process_arch" && [ "$(echo "$ENVOUT" | awk -F'\t' '$5=="env_process_arch" {print $6}')" != "$(echo "$ENVOUT" | awk -F'\t' '$5=="env_os_arch" {print $6}')" ]; then refuse ".NET reports a process architecture different from the OS architecture: emulation"; fi
fi
say
say "What to set in Docker Desktop (Settings > Resources) and why"
say "  Memory:      at least 8 GB for the generous steps (a, b, c); 16 GB if the Mac has it. Step d limits the memory per container below that, so the VM must be bigger than the biggest container limit."
say "  CPUs:        4 or more (the protocol uses 4 threads)."
say "  Swap:        0 or the smallest value (the containers are started without swap)."
say "  Disk image size: at least 100 GB (the 3 Gbp run needs about $(human "$(need_total_bytes 3000000000)") at its peak; the default of Docker Desktop is smaller). After enlarging it Docker restarts."
say "  File sharing: nothing to do. The data lives in a named volume inside the VM (not in a macOS folder): a bind mount through virtiofs / gRPC FUSE is slow and unreliable for mmap."
say "  Apple Silicon: the image is built for linux/arm64 and runs natively. Never run the amd64 image there (emulation); the wrapper refuses it."
say
if [ "$REFUSALS" -gt 0 ]; then say "RESULT: REFUSED ($REFUSALS problem(s), $WARNINGS warning(s)). Fix the lines marked [REFUSE] first."; row preflight_result REFUSED "" "$REFUSALS refusals, $WARNINGS warnings"; exit 3; fi
say "RESULT: OK ($WARNINGS warning(s))"; row preflight_result OK "" "$WARNINGS warnings"; exit 0
