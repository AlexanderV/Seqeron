#!/bin/bash
# Self-test of the cgroup readers of E8 (shell: e8lib.sh, Python: sample_e8.py) against SYNTHETIC cgroup v1 and v2 trees, and against the real hierarchy of the machine it runs on (whichever version that is).
# usage: selftest_e8.sh        prints PASS / FAIL per check, exit code 1 on any failure.
# What this proves and what it does not: the parsing of the files in the layout the kernel documents; NOT that Docker Desktop's VM presents exactly these files (that is what the runbook asks the owner to look at in the preflight).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
FAILS=0
check() { if [ "$2" = "$3" ]; then echo "PASS $1"; else echo "FAIL $1: got '$2', expected '$3'"; FAILS=$((FAILS + 1)); fi; }

# ---- synthetic cgroup v2 (a container with its own cgroup namespace: the group is the mount root)
V2=$T/v2; mkdir -p "$V2/sys" "$V2/proc/self"
echo "cpuset cpu io memory pids" > "$V2/sys/cgroup.controllers"
echo 1073741824 > "$V2/sys/memory.max"; echo 500000000 > "$V2/sys/memory.current"; echo 777777777 > "$V2/sys/memory.peak"; echo 0 > "$V2/sys/memory.swap.max"
printf 'anon 100000000\nfile 300000000\nfile_mapped 50000000\nfile_dirty 1000000\nfile_writeback 2000\n' > "$V2/sys/memory.stat"
printf 'low 0\nhigh 0\nmax 42\noom 1\noom_kill 1\n' > "$V2/sys/memory.events"
echo "200000 100000" > "$V2/sys/cpu.max"
echo "0::/" > "$V2/proc/self/cgroup"
printf 'MemTotal:       16000000 kB\nSwapTotal:             0 kB\n' > "$V2/proc/meminfo"
( export E8_CGROUP_ROOT=$V2/sys E8_PROC=$V2/proc; . "$HERE/e8lib.sh"
  cg_detect; check "v2 version" "$CG_VERSION" 2
  check "v2 limit" "$(cg_limit_bytes)" 1073741824
  check "v2 swap" "$(cg_swap_limit_bytes)" 0
  check "v2 cpu" "$(cg_cpu_limit)" 2.00
  check "v2 usage" "$(cg_usage_bytes)" 500000000
  check "v2 failures (events max)" "$(cg_failures)" 42
  check "v2 mem total" "$(mem_total_bytes)" 16384000000
  echo max > "$E8_CGROUP_ROOT/memory.max"; check "v2 no limit" "$(cg_limit_bytes)" max )
# the same through the Python reader
V2PY=$(E8_CGROUP_ROOT=$V2/sys E8_PROC=$V2/proc python3 -I -c "
import importlib.util, os
spec = importlib.util.spec_from_file_location('s', '$HERE/sample_e8.py'); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
v, d = m.cg_find(); s = m.cg_snapshot(v, d)
print(v, m.cg_limit(v, d), s['usage'], s['anon'], s['file'], s['mapped'], s['dirty'], s['writeback'], s['fail'], s['oom_kill'], m.cg_peak(v, d))")
check "v2 python reader" "$V2PY" "2 None 500000000 100000000 300000000 50000000 1000000 2000 42 1 777777777"   # memory.max is 'max' after the shell test above
# ---- synthetic cgroup v1 (a container: the group is the 'memory' mount, paths below it are empty)
V1=$T/v1; mkdir -p "$V1/sys/memory" "$V1/sys/cpu" "$V1/proc/self"
echo 2147483648 > "$V1/sys/memory/memory.limit_in_bytes"; echo 2147483648 > "$V1/sys/memory/memory.memsw.limit_in_bytes"
echo 900000000 > "$V1/sys/memory/memory.usage_in_bytes"; echo 950000000 > "$V1/sys/memory/memory.max_usage_in_bytes"; echo 5 > "$V1/sys/memory/memory.failcnt"
printf 'cache 400000000\nrss 123000000\nmapped_file 60000000\ndirty 3000\nwriteback 0\ntotal_cache 400000000\n' > "$V1/sys/memory/memory.stat"
printf 'oom_kill_disable 0\nunder_oom 0\noom_kill 2\n' > "$V1/sys/memory/memory.oom_control"
echo 150000 > "$V1/sys/cpu/cpu.cfs_quota_us"; echo 100000 > "$V1/sys/cpu/cpu.cfs_period_us"
printf '4:memory:/docker/abc\n1:cpu:/docker/abc\n' > "$V1/proc/self/cgroup"
( export E8_CGROUP_ROOT=$V1/sys E8_PROC=$V1/proc; . "$HERE/e8lib.sh"
  cg_detect; check "v1 version" "$CG_VERSION" 1
  check "v1 limit" "$(cg_limit_bytes)" 2147483648
  check "v1 swap allowance" "$(cg_swap_limit_bytes)" 0
  check "v1 cpu" "$(cg_cpu_limit)" 1.50
  check "v1 usage" "$(cg_usage_bytes)" 900000000
  check "v1 failures (failcnt)" "$(cg_failures)" 5
  echo 9223372036854771712 > "$V1/sys/memory/memory.limit_in_bytes"; check "v1 no limit (2^63)" "$(cg_limit_bytes)" max )
V1PY=$(E8_CGROUP_ROOT=$V1/sys E8_PROC=$V1/proc python3 -I -c "
import importlib.util
spec = importlib.util.spec_from_file_location('s', '$HERE/sample_e8.py'); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
v, d = m.cg_find(); s = m.cg_snapshot(v, d)
print(v, m.cg_limit(v, d), s['usage'], s['anon'], s['file'], s['mapped'], s['dirty'], s['writeback'], s['fail'], s['oom_kill'], m.cg_peak(v, d))")
check "v1 python reader" "$V1PY" "1 None 900000000 123000000 400000000 60000000 3000 0 5 2 950000000"   # the limit was set to 2^63 above
# ---- the real hierarchy of this machine
( . "$HERE/e8lib.sh"; cg_detect; echo "INFO real hierarchy: version $CG_VERSION, dir $CG_DIR, limit $(cg_limit_bytes), swap allowance $(cg_swap_limit_bytes), cpu $(cg_cpu_limit)" )
REAL=$(python3 -I -c "
import importlib.util
spec = importlib.util.spec_from_file_location('s', '$HERE/sample_e8.py'); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
v, d = m.cg_find(); print(v, d, m.cg_limit(v, d) if d else None, m.cg_snapshot(v, d) if d else None)")
echo "INFO real hierarchy through python: $REAL"
[ "$FAILS" = 0 ] && echo "ALL PASS" || { echo "$FAILS FAILED"; exit 1; }
