#!/bin/bash
# Pre-flight of the hand-off package (work package 1b.9f). Runs INSIDE the container (or natively on Linux) and says, with numbers, whether this machine can do the run at --n symbols, and what to change.
# usage: preflight.sh [--n N] [--kind genome|dna|rep|prot] [--quick] [--no-form]
#   --n      symbols of the run to plan for (default 3000000000); decimal suffixes k m g allowed
#   --quick  skip the write / mmap / msync test of the data directory
# Prints [ ok ] / [warn] / [REFUSE] lines, writes $S1B_DATA/env/env.tsv (what the form takes) and ends with RESULT: OK or REFUSED. Exit code 0 = no refusal (warnings allowed), 3 = at least one refusal.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=s1blib.sh
. "$HERE/s1blib.sh"
DATA=${S1B_DATA:-/data}
BENCH=${S1B_BENCH:-/opt/s1b/bin/bench/SuffixTree.Esa.Bench.dll}
N=3000000000; KIND=genome; QUICK=0; WRITE=1
while [ $# -gt 0 ]; do
  case $1 in
    --n) N=$(size_n "$2"); shift 2 ;;
    --kind) KIND=$2; shift 2 ;;
    --quick) QUICK=1; shift ;;
    --no-form) WRITE=0; shift ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
REFUSALS=0; WARNINGS=0
ok()     { say "  [ ok ] $*"; }
warn()   { say "  [warn] $*"; WARNINGS=$((WARNINGS + 1)); }
refuse() { say "  [REFUSE] $*"; REFUSALS=$((REFUSALS + 1)); }
if [ "$WRITE" = 1 ]; then ENVDIR=$DATA/env; mkdir -p "$ENVDIR" 2> /dev/null; else ENVDIR=$(mktemp -d); fi
ENVF=$ENVDIR/env.tsv
: > "$ENVF" 2> /dev/null
kv() { printf '%s\t%s\n' "$1" "$2" >> "$ENVF"; }

SYS=$(uname -s); UNAME=$(uname -m)
say "S1B pre-flight  ($(date -u +%Y-%m-%dT%H:%M:%SZ))   plan for n = $N, text kind $KIND, data directory $DATA"
say
say "1. Where this runs, processor architecture, emulation"
INC=no; { [ -f /.dockerenv ] || grep -qaE 'docker|containerd|kubepods' /proc/1/cgroup 2> /dev/null; } && INC=yes
kv in_container "$INC"; kv arch "$UNAME"; kv kernel "$(uname -r)"
IMGARCH=${S1B_IMAGE_ARCH:-unknown}; kv image_arch "$IMGARCH"
say "     system $SYS, uname -m $UNAME, in a container: $INC, image built for: $IMGARCH"
if [ "$SYS" != Linux ]; then
  refuse "this is $SYS, not Linux. The 320 MB gate is a Linux RssAnon gate (/proc/<pid>/status has no counterpart on macOS: ps shows RSS with the mapped file pages), and the disk reservation of the build (fallocate) is Linux code. Run the Docker image (the wrapper s1b-docker.sh): inside it everything is Linux, also on a Mac"
  kv rssanon_source ps-rss
else
  kv rssanon_source proc
fi
CPUFMT=unknown
grep -qE '^(CPU implementer|Features)' /proc/cpuinfo 2> /dev/null && CPUFMT=arm
grep -q '^flags' /proc/cpuinfo 2> /dev/null && CPUFMT=x86
EMU=no
case "$UNAME:$CPUFMT" in x86_64:arm) EMU="yes: an x86-64 process on an ARM CPU (Rosetta or QEMU)" ;; aarch64:x86) EMU="yes: an ARM process on an x86 CPU" ;; esac
grep -qaE 'rosetta|qemu' /proc/self/maps 2> /dev/null && EMU="yes: an emulator is mapped into this process (/proc/self/maps)"
case "$IMGARCH:$UNAME" in amd64:aarch64|arm64:x86_64) EMU="yes: the image is for $IMGARCH but the kernel reports $UNAME" ;; esac
kv emulation "$EMU"
if [ "$EMU" = no ]; then ok "native: the process runs on the architecture it was built for ($UNAME)"; else refuse "emulation detected ($EMU). Timings under emulation mean nothing and the arm64 instruction paths are not run. Build the image for the Docker daemon's architecture (s1b-docker.sh image) and remove the other one"; fi
case $UNAME in aarch64|arm64) say "     arm64: the CRC-32C instruction path, the scalar common-prefix path and the Vector128 window decode are the ones that run here; the container must be byte-identical to the x86-64 reference";; esac
say
say "2. Processors"
CORES=$(nproc 2> /dev/null || sysctl -n hw.ncpu 2> /dev/null || echo 1); CPULIM=$(cg_cpu_limit 2> /dev/null || echo none)
kv cores "$CORES"; kv cpu_limit "$CPULIM"
say "     cores visible: $CORES, CPU quota of the container: $CPULIM"
if [ "$CORES" -ge 4 ]; then ok "$CORES cores (the build runs four threads; the gate is the figure at four threads)"; else warn "only $CORES core(s): the build still runs four threads (the gate is at four), only slower"; fi
say
say "3. Memory"
MEMT=$(mem_total_bytes 2> /dev/null || echo 0); CGL=$(cg_limit_bytes 2> /dev/null || echo none); cg_detect 2> /dev/null; SWAPT=$(swap_total_bytes 2> /dev/null || echo 0)
kv mem_total_bytes "$MEMT"; kv cgroup_version "${CG_VERSION:-none}"; kv cgroup_memory_limit "$CGL"; kv swap_total_bytes "${SWAPT:-0}"
say "     memory visible: $(human "$MEMT") (inside Docker Desktop: the memory of its Linux VM), limit of this container: $([ "$CGL" = max ] || [ "$CGL" = none ] && echo none || human "$CGL"), swap of the VM: $(human "${SWAPT:-0}")"
EFF=$MEMT; if [ "$CGL" != max ] && [ "$CGL" != none ] && [ "$CGL" -lt "$MEMT" ]; then EFF=$CGL; fi
if [ "$EFF" -ge 4000000000 ]; then ok "$(human "$EFF") available: the build holds about 0.3 GB of anonymous memory (the gate), the rest is page cache for the text, the rank file and the container being written"
elif [ "$EFF" -ge 1000000000 ]; then warn "only $(human "$EFF") available: the build needs about 0.35 GB of anonymous memory and runs, but the page cache it relies on is small: expect it slower (the E8 sweep: below text + rank file the time grows). Docker Desktop: Settings > Resources > Memory, 8 GB or more"
else refuse "only $(human "$EFF") of memory: less than the anonymous memory of the build (0.3 GB) plus a working page cache"; fi
if [ "${SWAPT:-0}" -gt 0 ]; then warn "swap is on ($(human "$SWAPT")): under memory pressure the kernel may swap anonymous pages out, and RssAnon (resident) then reads LOWER than the memory the build asked for; a run with swap in use proves nothing about the gate. The sampler records it only as RssAnon. Prefer Docker Desktop's swap at its minimum"; fi
if [ "$N" -le "$S1B_READER_LIMIT" ]; then STRUCT_RSS=$(awk -v n="$N" 'BEGIN { printf "%.0f", n / 8 + 40e6 }'); say "     the Structure level of the reader needs about $(human "$STRUCT_RSS") of anonymous memory at this size (a bitmap of n / 8 bytes and 34 MB of buffers; measured 42 MB at 10^7, 109 MB at 10^8): it is not the build's gate"; fi
say
say "4. Disk of the data directory"
mkdir -p "$DATA" 2> /dev/null
if [ ! -d "$DATA" ] || [ ! -w "$DATA" ]; then refuse "$DATA does not exist or is not writable"; else
  AVAIL=$(free_bytes_of "$DATA"); SIZE=$(size_bytes_of "$DATA")
  kv disk_free_bytes "$AVAIL"; kv disk_size_bytes "$SIZE"
  say "     free: $(human "$AVAIL") of $(human "$SIZE") on the file system of the data directory"
  say "     need at the peak of one run (the text n x 1, the container n x $(container_bpc "$KIND") reserved from the start, the scratch n x $S1B_SCRATCH_BPC at the peak, plus 2 GB and 2 %):"
  for n in $(printf '%s\n' 100000000 1000000000 3000000000 "$N" | sort -un); do
    nd=$(need_bytes "$n" "$KIND"); say "       n = $(printf '%13s' "$n"): $(human "$nd") (text $(human "$n"), container $(human "$(awk -v n="$n" -v c="$(container_bpc "$KIND")" 'BEGIN { printf "%.0f", n * c }')"), scratch $(human $((n * S1B_SCRATCH_BPC))))"
  done
  NEED=$(need_bytes "$N" "$KIND")
  if [ "$AVAIL" -ge $((NEED + NEED / 4)) ]; then ok "the run at n = $N fits with room ($(human "$AVAIL") free, $(human "$NEED") needed)"
  elif [ "$AVAIL" -ge "$NEED" ]; then warn "the run at n = $N fits barely ($(human "$AVAIL") free, $(human "$NEED") needed): a second run beside it will not"
  else refuse "the run at n = $N needs about $(human "$NEED") and only $(human "$AVAIL") are free. Free disk (Docker Desktop: the disk image size, Settings > Resources) or plan a smaller --n"; fi
fi
say
say "5. File system of the data directory (a macOS bind mount is slow and unreliable for the mapped files of the build)"
FST=$(fs_type "$DATA" 2> /dev/null || echo unknown); kv fs_type "$FST"
say "     type: $FST"
if bad_fs "$FST" && [ "${S1B_ALLOW_BADFS:-0}" != 1 ]; then refuse "$FST is a network or user-space file system. A folder of the Mac mounted into the container (a bind mount) looks like this (virtiofs, fuse, grpcfuse, osxfs): mmap with a large dirty set is slow and not reliable on it, and the numbers would be those of the mount. Use the named volume of s1b-docker.sh (the default), not -v <folder>:/data"
else
  case "$FST" in
    tmpfs*|ramfs*) if [ "${S1B_ALLOW_TMPFS:-0}" = 1 ]; then warn "tmpfs: the data is in RAM (allowed by S1B_ALLOW_TMPFS, tests only)"; else refuse "tmpfs: the data directory is in RAM"; fi ;;
    overlay*) warn "overlay: the data directory is in the container's writable layer, not in a volume: slower, lost with the container, and limited by the VM disk. Use the named volume of s1b-docker.sh" ;;
    *) ok "$FST" ;;
  esac
  if [ "$INC" = yes ] && command -v mountpoint > /dev/null && ! mountpoint -q "$DATA" 2> /dev/null; then warn "$DATA is not a mount point of its own (a named volume is expected there)"; fi
fi
say
say "6. Write, mmap, msync and fallocate on the data directory (64 MB)"
if [ "$QUICK" = 1 ]; then say "     skipped (--quick)"; else
  IOOUT=$(python3 -I - "$DATA/.s1b-io-test" <<'PY'
import mmap, os, sys, time
path = sys.argv[1]; size = 64 << 20; res = {}
try:
    fd = os.open(path, os.O_RDWR | os.O_CREAT | os.O_TRUNC, 0o644)
    t = time.time()
    try:
        os.posix_fallocate(fd, 0, size); res['fallocate'] = 'ok'
    except (OSError, AttributeError) as e:
        res['fallocate'] = 'unsupported'; os.ftruncate(fd, size)
    mm = mmap.mmap(fd, size); pattern = bytes(range(256)) * 4096
    for off in range(0, size, len(pattern)):
        mm[off:off + len(pattern)] = pattern
    mm.flush(); res['mmap_write_msync_s'] = '%.3f' % (time.time() - t); mm.close(); os.fsync(fd)
    with open(path, 'rb') as f:
        res['readback'] = 'ok' if f.read() == pattern * (size // len(pattern)) else 'MISMATCH'
    res['write_mb_s'] = '%.0f' % (size / 1e6 / max(1e-6, float(res['mmap_write_msync_s']))); os.close(fd)
except Exception as e:
    res['error'] = repr(e)
finally:
    try: os.unlink(path)
    except OSError: pass
print(' '.join('%s=%s' % kv for kv in res.items()))
PY
)
  say "     $IOOUT"; kv io_test "$IOOUT"
  case "$IOOUT" in
    *error=*|*MISMATCH*|'') refuse "the write / mmap / msync test failed: $IOOUT" ;;
    *) ok "mmap write, msync and read-back work"
       case "$IOOUT" in *fallocate=unsupported*) warn "fallocate is not supported here: the container is then not reserved up front (a full disk would show as SIGBUS instead of an IOException)" ;; esac ;;
  esac
fi
say
say "7. .NET, the program, the CRC-32C path"
if ! command -v dotnet > /dev/null; then refuse "dotnet is not on the PATH"; else
  DOTNETV=$(dotnet --list-runtimes 2> /dev/null | grep Microsoft.NETCore.App | tail -1 | cut -d' ' -f2); kv dotnet "$DOTNETV"
  say "     runtime: Microsoft.NETCore.App $DOTNETV"
  if [ ! -f "$BENCH" ]; then refuse "$BENCH not found (the program is not here)"; else
    dotnet "$BENCH" sesa handoff env --label env > "$ENVDIR/handoff.out" 2> "$ENVDIR/handoff.err" || refuse "the program's environment check failed: $(head -c 300 "$ENVDIR/handoff.err")"
    awk -F'\t' '$1 == "R" { printf "     %-34s %s\n", $3 ":", $4 }' "$ENVDIR/handoff.out" | grep -E "process_arch|crc32c|vector" | sed 's/^/ /'
    grep -q "env_crc32c_check_value	ok" "$ENVDIR/handoff.out" || refuse "CRC-32C does not give the catalogue's check value E3069283 on this machine"
    grep -q "env_crc32c_table_equals_hardware	yes" "$ENVDIR/handoff.out" || refuse "the CRC-32C table path and the hardware path give different values on this machine"
    grep -q "env_crc32c_hardware	no" "$ENVDIR/handoff.out" && warn "no CRC-32C instruction: the managed table path runs (same values, slower); the instruction path is not exercised"
  fi
fi
if [ -r /proc/self/status ] && grep -q '^RssAnon:' /proc/self/status; then ok "RssAnon is readable from /proc/<pid>/status (the figure of the 320 MB gate)"; else warn "no RssAnon in /proc/self/status: the independent sampler falls back to ps (RSS with mapped file pages): the gate cannot be judged here"; kv rssanon_source ps-rss; fi
command -v python3 > /dev/null && ok "python3 $(python3 -c 'import sys; print(sys.version.split()[0])')" || refuse "python3 is not available (the sampler, the generators and the checker are Python)"
dotnet --list-sdks 2> /dev/null | grep -q . && ok ".NET SDK present: the test sets (s1b.sh tests) can run" || say "     no .NET SDK: only the chain can run here, not the test sets"
say
if [ "$REFUSALS" -gt 0 ]; then say "RESULT: REFUSED ($REFUSALS problem(s), $WARNINGS warning(s)). Fix the lines marked [REFUSE] first."; kv preflight_result REFUSED; exit 3; fi
say "RESULT: OK ($WARNINGS warning(s))"; kv preflight_result OK; exit 0
