#!/bin/bash
# Work package 1b.9d: the window of the package at a real size, in a child process (a signal shows as an exit code, not as the death of a test host): a tmpfs that holds the build with a little room is mounted, the build
# of a text of 10^8 symbols whose LCP is four bytes per entry (the container grows by 300 MB behind its reserved length) runs on it, and at the point lcp-placed (after the late pre-flight and the growth) the rest of the
# volume is taken by a ballast file. Two builds: `--sparse true` (the growth is a hole, the behaviour before 1b.9d: expected to die of SIGBUS, exit code 135) and the default (the growth is reserved: the build ends, exit code 0
# without PSI, or with an IOException, exit code 3, when the next file cannot get its blocks, with PSI).
# usage: run_1b9d_fullvolume.sh <rows file> <input text> <volume bytes>     (root: mounts a tmpfs on /tmp/stbench/b9d/vol and removes it at the end)
set -u
ROWS=$(readlink -f "$1"); INPUT=$2; BYTES=$3
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
unset DOTNET_GCHeapHardLimit
DLL=${B9D_DLL:-/tmp/stbench/b9d/bin/SuffixTree.Esa.Bench.dll}
V=/tmp/stbench/b9d/vol
mkdir -p $V /tmp/stbench/b9d/volref
# the reference: the same build without PSI on the ordinary disk, nothing filled; its SHA-256 is what the build on the full volume must give
dotnet "$DLL" sesa build --input "$INPUT" --work /tmp/stbench/b9d/volref --label ref --psi false --threads 4 --fsync false --keep true > /dev/null 2>&1
REF=$(sha256sum /tmp/stbench/b9d/volref/out.sesa | cut -d' ' -f1); rm -rf /tmp/stbench/b9d/volref
for cell in "sparse psi=false" "reserved psi=false" "sparse psi=true" "reserved psi=true"; do
  set -- $cell; MODE=$1; PSI=${2#psi=}
  mount -t tmpfs -o size=$BYTES tmpfs $V || { echo "cannot mount" >&2; exit 2; }
  ARGS=""; [ "$MODE" = sparse ] && ARGS="--sparse true"
  dotnet "$DLL" sesa build --input "$INPUT" --work $V --label "full.$MODE.psi_$PSI" --psi "$PSI" --threads 4 --fsync false --keep true --fill-at lcp-placed $ARGS > $V.out 2> $V.err
  RC=$?
  if [ -f $V/out.sesa ] && [ "$PSI" = false ]; then
    GOT=$(sha256sum $V/out.sesa | cut -d' ' -f1)
    printf 'R\tfullvolume\tfull.%s.psi_%s\tcontainer_equals_reference\t%s\t\tSHA-256 of the container built on the full volume against the build on the ordinary disk (%s)\n' "$MODE" "$PSI" "$([ "$GOT" = "$REF" ] && echo yes || echo NO)" "${REF:0:16}" >> "$ROWS"
  fi
  if [ $RC = 0 ]; then LEFT=$(ls -A $V | grep -vc '^ballast$\|^out.sesa$'); else LEFT=$(ls -A $V | grep -vc '^ballast$'); fi   # a build that succeeded leaves its container
  printf 'R\tfullvolume\tfull.%s.psi_%s\texit_code\t%s\t\t0 = built, 3 = IOException, 135 = SIGBUS\n' "$MODE" "$PSI" "$RC" >> "$ROWS"
  printf 'R\tfullvolume\tfull.%s.psi_%s\tleft_behind\t%s\t\tentries of the volume besides the ballast\n' "$MODE" "$PSI" "$LEFT" >> "$ROWS"
  grep -E "io_exception|total_s|grown_bytes" $V.out | sed "s/^R\tsesa/R\tfullvolume/" >> "$ROWS"
  tail -2 $V.err | sed 's/^/#stderr: /' >> "$ROWS"
  umount $V
done
rm -f $V.out $V.err
