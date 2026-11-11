#!/bin/bash
# E4 helper (outside the shipped code): run the EmProto E4 (PH4) build inside a cgroup-v1 memory sub-group with a memory limit, cold page cache.
# usage: run_e4_cgroup.sh <limit-bytes|max> <label> <emproto args for the e4 command...>
# The limit is the TOTAL of the sub-group (anonymous memory + page cache incl. mapped file pages + kernel memory charged to it).
# Needs root and a mounted cgroup-v1 memory controller (mkdir in /sys/fs/cgroup/memory is allowed in the prototype VM). Nothing outside the sub-group is limited.
set -u
LIMIT=$1; LABEL=$2; shift 2
CG=/sys/fs/cgroup/memory/emproto_e4_$LABEL
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
EMPROTO_DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e4_run/SuffixTree.EmProto.dll}
mkdir "$CG" || { echo "cannot create $CG" >&2; exit 2; }
if [ "$LIMIT" != "max" ]; then echo "$LIMIT" > "$CG/memory.limit_in_bytes"; fi
[ -w "$CG/memory.swappiness" ] && echo 0 > "$CG/memory.swappiness"
sync; echo 3 > /proc/sys/vm/drop_caches   # cold start: the text, SA and PLCP pages are charged to the sub-group when first read
(
  echo $BASHPID > "$CG/cgroup.procs"
  exec dotnet "$EMPROTO_DLL" e4 --cgroup-usage "$CG/memory.usage_in_bytes" "$@"
)
RC=$?
echo "# cgroup $LABEL limit=$LIMIT rc=$RC failcnt=$(cat $CG/memory.failcnt) max_usage=$(cat $CG/memory.max_usage_in_bytes)" >&2
grep -E "^(total_)?(cache|rss|mapped_file|pgpgin|pgpgout|pgfault|pgmajfault|writeback|dirty) " "$CG/memory.stat" | tr '\n' ' ' >&2; echo >&2
rmdir "$CG" 2>/dev/null
exit $RC
