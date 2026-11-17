#!/bin/bash
# E8 host-side wrapper for macOS (Docker Desktop) and Linux shells. Written for bash 3.2 (the bash of macOS). UNBUILT AND UNTESTED ON macOS AND ON arm64; --dry-run prints the docker commands without running them.
#
# usage: e8-docker.sh [global options] <command> [arguments]
#   build                      build the image natively for this machine's Docker architecture (never the other one: an amd64 image on Apple Silicon would run under emulation)
#   preflight [--sizes ".."]   run the preflight in a container (no memory limit)
#   plan | status              the protocol and its state
#   gen <size>...              generate texts (1g 2.2g 3g) into the data volume
#   run a|b|c|e [..]           generous steps (no memory limit), with the full verification
#   run d --size 1g --memory-mb M   ONE cache cell (starts the container with --memory M --memory-swap M)
#   sweep [--cells "M1 M2 .."] [--size 1g]   step a if needed, then every cell of step d, largest first; stops a cell after 6 x the generous time
#   report                     tables, gates and the return archive (inside the volume)
#   export <folder>            copy rows, sums, logs and the return archives out of the volume into <folder>
#   shell                      a shell in the container
#   clean-volume               delete the data volume (asks first)
# global options (before the command):
#   --image NAME (seqeron-e8)  --volume NAME (seqeron-e8-data)  --memory MiB (container limit; default none = generous)  --cpus N  --cold (start the container --privileged so that every run can drop the page cache)
#   --detach (run in the background: docker logs -f seqeron-e8-run follows it; survives a closed terminal)  --dry-run
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../../../.." && pwd)
DOCKERFILE=docs/refactoring/suffix-tree-upgrade/analysis/e8/Dockerfile
IMAGE=seqeron-e8; VOLUME=seqeron-e8-data; MEMORY=; CPUS=; COLD=0; DETACH=0; DRY=0
while [ $# -gt 0 ]; do
  case $1 in
    --image) IMAGE=$2; shift 2 ;;
    --volume) VOLUME=$2; shift 2 ;;
    --memory) MEMORY=$2; shift 2 ;;
    --cpus) CPUS=$2; shift 2 ;;
    --cold) COLD=1; shift ;;
    --detach) DETACH=1; shift ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) break ;;
  esac
done
CMD=${1:-help}; [ $# -gt 0 ] && shift

die() { echo "e8-docker: $*" >&2; exit 2; }
run_cmd() { if [ "$DRY" = 1 ]; then echo "+ $*"; else "$@"; fi; }
command -v docker > /dev/null || die "docker is not installed or not on the PATH (install Docker Desktop and start it)"

# the architecture of the Docker daemon (on Docker Desktop for Apple Silicon: aarch64 = arm64; on an Intel Mac: x86_64)
DARCH=$(docker info --format '{{.Architecture}}' 2>/dev/null)
if [ -z "$DARCH" ]; then
  [ "$DRY" = 1 ] || die "cannot reach the Docker daemon: start Docker Desktop (or the docker service) and try again"
  DARCH=$(uname -m)
fi
case $DARCH in
  aarch64|arm64) ARCH=arm64 ;;
  x86_64|amd64) ARCH=amd64 ;;
  *) die "unknown Docker architecture '$DARCH'" ;;
esac
PLATFORM=linux/$ARCH
if [ -n "${E8_PLATFORM:-}" ] && [ "$E8_PLATFORM" != "$PLATFORM" ]; then
  die "E8_PLATFORM=$E8_PLATFORM differs from the architecture of the Docker daemon ($PLATFORM): that image would run under emulation and its timings would be meaningless. Refused."
fi

image_arch() { docker image inspect -f '{{.Architecture}}' "$IMAGE" 2>/dev/null; }
need_image() {
  [ "$DRY" = 1 ] && return 0
  local a; a=$(image_arch) || true
  [ -n "$a" ] || die "image '$IMAGE' not found: run  e8-docker.sh build  first"
  [ "$a" = "$ARCH" ] || die "image '$IMAGE' is built for $a but the Docker daemon is $ARCH: rebuild with  e8-docker.sh build  (an image for the other architecture would run under emulation)"
}

# docker run arguments common to every container
base_args() {
  BASE="--rm --init --platform $PLATFORM -v $VOLUME:/data"
  [ "$COLD" = 1 ] && BASE="$BASE --privileged"
  [ -n "$CPUS" ] && BASE="$BASE --cpus $CPUS"
  if [ -n "$MEMORY" ]; then BASE="$BASE --memory ${MEMORY}m --memory-swap ${MEMORY}m"; fi   # no swap: --memory-swap equal to --memory
  if [ -t 0 ] && [ -t 1 ] && [ "$DETACH" = 0 ]; then BASE="$BASE -it"; fi
}

# in_container [extra docker args...] -- <e8.sh arguments>
in_container() {
  local extra=""
  while [ $# -gt 0 ] && [ "$1" != "--" ]; do extra="$extra $1"; shift; done
  [ "$1" = "--" ] && shift
  base_args
  local name=seqeron-e8-run
  if [ "$DETACH" = 1 ]; then
    if [ "$DRY" = 0 ] && docker ps -a --format '{{.Names}}' | grep -qx "$name"; then
      if docker ps --format '{{.Names}}' | grep -qx "$name"; then die "$name is still running (docker logs -f $name follows it; docker stop $name stops it, a restart of the same command continues where it stopped)"; fi
      docker rm "$name" > /dev/null   # a finished earlier container of the same name
    fi
    # shellcheck disable=SC2086
    run_cmd docker run -d --name $name ${BASE#--rm } $extra "$IMAGE" "$@"
    echo "started in the background: docker logs -f $name   (it stops by itself; run the same command again after an interruption: finished runs are skipped)"
  else
    # shellcheck disable=SC2086
    run_cmd docker run $BASE $extra "$IMAGE" "$@"
  fi
}

case $CMD in
  build)
    COMMIT=$(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)
    if [ "$COMMIT" != unknown ] && [ -n "$(git -C "$ROOT" status --porcelain --untracked-files=no 2>/dev/null)" ]; then echo "note: the checkout has uncommitted changes; the image is built from the files as they are (label revision: $COMMIT)"; fi
    echo "building $IMAGE for $PLATFORM from $(basename "$ROOT") (the build context is the repository root, filtered by Dockerfile.dockerignore)"
    run_cmd docker build --platform "$PLATFORM" -f "$ROOT/$DOCKERFILE" --build-arg "COMMIT=$COMMIT" -t "$IMAGE" "$ROOT" ;;
  preflight) need_image; in_container -- preflight "$@" ;;
  plan|status) need_image; in_container -- plan ;;
  gen) need_image; in_container -- gen "$@" ;;
  run) need_image; in_container -- run "$@" ;;
  sweep)
    need_image
    CELLS=; SIZE=1g
    while [ $# -gt 0 ]; do case $1 in --cells) CELLS=$2; shift 2 ;; --size) SIZE=$2; shift 2 ;; *) die "unknown option $1" ;; esac; done
    if [ -z "$CELLS" ]; then
      if [ "$DRY" = 1 ]; then CELLS="8448 2432 1952 1568 1344 1184 1088 992"; else CELLS=$(MEMORY= in_container -- cells 2>/dev/null | tail -1 | tr -d '\r'); fi
    fi
    VMMEM=$(docker info --format '{{.MemTotal}}' 2>/dev/null || echo 0)
    echo "step a first (generous, skipped when done), then the cells (MiB): $CELLS"
    MEMORY= in_container -- run a || die "step a did not finish: the cells need its time as the baseline"
    for M in $CELLS; do
      if [ "$DRY" = 0 ] && [ "$VMMEM" -gt 0 ] && [ $((M * 1048576)) -gt "$VMMEM" ]; then echo "cell $M MiB skipped: more than the memory of the Docker VM ($((VMMEM / 1048576)) MiB)"; continue; fi
      echo "---- cell $M MiB"
      MEMORY=$M in_container -- run d --size "$SIZE" --memory-mb "$M" || echo "cell $M MiB: the container returned an error (see the rows; the next cell follows)"
    done ;;
  report) need_image; in_container -- report ;;
  export)
    DEST=${1:?usage: e8-docker.sh export <folder>}
    need_image
    mkdir -p "$DEST" || die "cannot create $DEST"
    run_cmd docker run --rm --platform "$PLATFORM" -v "$VOLUME:/data:ro" -v "$DEST:/out" --entrypoint /bin/bash "$IMAGE" -c 'cp -r /data/rows /data/sums /data/logs /data/return /out/ 2>/dev/null; ls -la /out' ;;
  shell) need_image; MEMORY=$MEMORY in_container --entrypoint /bin/bash -- ;;
  clean-volume)
    printf 'delete the Docker volume %s with everything in it (texts, outputs, rows)? type the volume name to confirm: ' "$VOLUME"
    read -r ANSWER
    [ "$ANSWER" = "$VOLUME" ] && run_cmd docker volume rm "$VOLUME" || echo "not deleted" ;;
  help|*) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//' ;;
esac
