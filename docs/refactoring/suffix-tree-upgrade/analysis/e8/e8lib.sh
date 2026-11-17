#!/bin/bash
# E8 shared shell helpers (sourced by e8.sh, preflight.sh, run_e8.sh). UNBUILT AND UNTESTED IN A CONTAINER, ON macOS AND ON arm64 (see E8-RUNBOOK.md).
# Everything that reads the kernel's view (cgroup, meminfo, mounts) takes its root from an environment variable so that the cgroup v2 reader can be tested against a synthetic tree:
#   E8_CGROUP_ROOT (default /sys/fs/cgroup), E8_PROC (default /proc).

E8_CGROUP_ROOT=${E8_CGROUP_ROOT:-/sys/fs/cgroup}
E8_PROC=${E8_PROC:-/proc}

# cg_detect: sets CG_VERSION (1, 2 or none) and CG_DIR (the directory that holds the memory files of THIS process's group; in a container with its own cgroup namespace that is the mount root)
cg_detect() {
  CG_VERSION=none; CG_DIR=
  local sub
  if [ -f "$E8_CGROUP_ROOT/cgroup.controllers" ]; then
    CG_VERSION=2
    sub=$(awk -F: '$1=="0" {print $3; exit}' "$E8_PROC/self/cgroup" 2>/dev/null)
    CG_DIR="$E8_CGROUP_ROOT$sub"
    [ -f "$CG_DIR/memory.max" ] || CG_DIR="$E8_CGROUP_ROOT"
    [ -f "$CG_DIR/memory.max" ] || [ -f "$CG_DIR/memory.current" ] || CG_DIR=
  elif [ -d "$E8_CGROUP_ROOT/memory" ] || [ -f "$E8_CGROUP_ROOT/memory.limit_in_bytes" ]; then
    CG_VERSION=1
    sub=$(awk -F: '$2 ~ /(^|,)memory(,|$)/ {print $3; exit}' "$E8_PROC/self/cgroup" 2>/dev/null)
    if [ -d "$E8_CGROUP_ROOT/memory" ]; then CG_DIR="$E8_CGROUP_ROOT/memory$sub"; [ -f "$CG_DIR/memory.limit_in_bytes" ] || CG_DIR="$E8_CGROUP_ROOT/memory"; else CG_DIR="$E8_CGROUP_ROOT"; fi
    [ -f "$CG_DIR/memory.limit_in_bytes" ] || CG_DIR=
  fi
  [ -n "$CG_DIR" ] || CG_VERSION=none
}

# cg_limit_bytes: the memory limit of the group in bytes, "max" when there is none (cgroup v1 reports a huge number: 2^60 and above counts as none)
cg_limit_bytes() {
  cg_detect
  local v
  case $CG_VERSION in
    2) v=$(cat "$CG_DIR/memory.max" 2>/dev/null) ;;
    1) v=$(cat "$CG_DIR/memory.limit_in_bytes" 2>/dev/null) ;;
    *) echo none; return ;;
  esac
  if [ -z "$v" ] || [ "$v" = max ] || [ "$v" -ge 1152921504606846976 ] 2>/dev/null; then echo max; else echo "$v"; fi
}

# cg_swap_limit_bytes: the swap the group may use IN ADDITION to its memory limit (v2: memory.swap.max; v1: memsw limit minus the limit); "max" = unlimited, "none" = not reported
cg_swap_limit_bytes() {
  cg_detect
  local v lim
  case $CG_VERSION in
    2) v=$(cat "$CG_DIR/memory.swap.max" 2>/dev/null); [ -n "$v" ] || { echo none; return; }; echo "$v" ;;
    1) v=$(cat "$CG_DIR/memory.memsw.limit_in_bytes" 2>/dev/null); lim=$(cat "$CG_DIR/memory.limit_in_bytes" 2>/dev/null)
       if [ -z "$v" ] || [ -z "$lim" ]; then echo none; elif [ "$v" -ge 1152921504606846976 ] 2>/dev/null; then echo max; else echo $((v - lim)); fi ;;
    *) echo none ;;
  esac
}

# cg_cpu_limit: CPUs the group may use ("max" = no quota), from cpu.max (v2) or cfs quota / period (v1)
cg_cpu_limit() {
  cg_detect
  local q p base
  case $CG_VERSION in
    2) read -r q p < "$CG_DIR/cpu.max" 2>/dev/null || { echo none; return; }
       if [ "$q" = max ]; then echo max; else awk -v q="$q" -v p="$p" 'BEGIN {printf "%.2f", q/p}'; fi ;;
    1) base=$E8_CGROUP_ROOT/cpu; [ -d "$base" ] || base=$E8_CGROUP_ROOT
       q=$(cat "$base/cpu.cfs_quota_us" 2>/dev/null); p=$(cat "$base/cpu.cfs_period_us" 2>/dev/null)
       if [ -z "$q" ] || [ -z "$p" ]; then echo none; elif [ "$q" -lt 0 ]; then echo max; else awk -v q="$q" -v p="$p" 'BEGIN {printf "%.2f", q/p}'; fi ;;
    *) echo none ;;
  esac
}

# cg_usage_bytes / cg_failures: current usage and the "hit the limit" counter (v1 memory.failcnt, v2 memory.events max)
cg_usage_bytes() {
  cg_detect
  case $CG_VERSION in
    2) cat "$CG_DIR/memory.current" 2>/dev/null ;;
    1) cat "$CG_DIR/memory.usage_in_bytes" 2>/dev/null ;;
  esac
}
cg_failures() {
  cg_detect
  case $CG_VERSION in
    2) awk '$1=="max" {print $2}' "$CG_DIR/memory.events" 2>/dev/null ;;
    1) cat "$CG_DIR/memory.failcnt" 2>/dev/null ;;
  esac
}

mem_total_bytes() { awk '/^MemTotal:/ {print $2 * 1024; exit}' "$E8_PROC/meminfo"; }
swap_total_bytes() { awk '/^SwapTotal:/ {print $2 * 1024; exit}' "$E8_PROC/meminfo"; }

# df_avail_bytes DIR / df_size_bytes DIR: free and total bytes of the file system of DIR
df_avail_bytes() { df -B1 --output=avail "$1" 2>/dev/null | tail -1 | tr -d ' '; }
df_size_bytes() { df -B1 --output=size "$1" 2>/dev/null | tail -1 | tr -d ' '; }

# fs_type DIR: the file system type as the kernel reports it (stat -f), and the mount entry's type when the two differ (fuse variants show up there)
fs_type() {
  local t m
  t=$(stat -f -c %T "$1" 2>/dev/null)
  m=$(awk -v d="$(readlink -f "$1")" '{ if (index(d, $2) == 1 && length($2) >= best) { best = length($2); typ = $3 } } END { print typ }' "$E8_PROC/mounts" 2>/dev/null)
  if [ -n "$m" ] && [ "$m" != "$t" ]; then echo "$t/$m"; else echo "$t"; fi
}

# human BYTES -> "12.3 GB" (decimal)
human() { awk -v b="$1" 'BEGIN { if (b >= 1e12) printf "%.2f TB", b/1e12; else if (b >= 1e9) printf "%.1f GB", b/1e9; else if (b >= 1e6) printf "%.0f MB", b/1e6; else printf "%d B", b }'; }

# size_n TEXT -> symbols: 1g, 2.2g, 3000m, 20m, 123456 ; decimal units
size_n() {
  awk -v s="$1" 'BEGIN { u = 1; l = tolower(substr(s, length(s), 1)); if (l == "k") u = 1e3; else if (l == "m") u = 1e6; else if (l == "g") u = 1e9; if (u > 1) s = substr(s, 1, length(s) - 1); printf "%.0f", s * u }'
}

# lbl_n N -> the label of a text size: n1000m, n2200m, n20m, n123456
lbl_n() { if [ $(( $1 % 1000000 )) -eq 0 ]; then echo "n$(( $1 / 1000000 ))m"; else echo "n$1"; fi; }

# disk need of one run at n symbols, in bytes (ESTIMATES of ADR-004: scratch 12 B/char at its peak, outputs 9.5 B/char (9.25-9.50 measured), the ASCII text is already on disk): peak total = n x (1 + 12 + 9.5)
need_run_bytes() { awk -v n="$1" 'BEGIN { printf "%.0f", n * 21.5 + 2e9 + n * 0.03 }'; }       # without the text (it exists) but with a 2 GB + 3 % margin
need_total_bytes() { awk -v n="$1" 'BEGIN { printf "%.0f", n * 22.5 + 2e9 + n * 0.03 }'; }     # with the text

# is_fuse_like TYPE: file systems that are slow or unreliable for mmap with a large dirty set
bad_fs() {
  case "$1" in
    *fuse*|*9p*|*virtiofs*|*nfs*|*cifs*|*smb*|*vboxsf*|*osxfs*|*grpcfuse*|*vmhgfs*|*afs*) return 0 ;;
  esac
  return 1
}

say() { printf '%s\n' "$*"; }
