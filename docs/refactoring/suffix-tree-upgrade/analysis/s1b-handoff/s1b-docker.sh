#!/bin/bash
# Host-side wrapper of the hand-off package s1b-handoff (work package 1b.9f) for macOS (Docker Desktop) and Linux shells. Written for bash 3.2 (the bash of macOS). The image is built and run for the architecture of
# the Docker daemon only (on an Apple-silicon Mac: linux/arm64, natively; an Intel Mac would run linux/arm64 under emulation: the wrapper refuses that, because timings under emulation mean nothing).
# STATUS: written and run on a linux/amd64 daemon (--dry-run also without one); NOT run on arm64 or on Docker Desktop. --dry-run prints the docker commands without running them.
#
# usage: s1b-docker.sh [global options] <command> [arguments]
#   image [--chain-only]       build the image natively for this machine's Docker architecture (the full one with the .NET SDK and the test sets; --chain-only: the runtime and the chain, 0.4 GB)
#   preflight [--n N]          the pre-flight in a container: what this machine is and whether the run at N symbols fits (refuses with numbers)
#   plan [--n N]               the disk, RAM and time tables
#   chain --n N [--kind K]     gen, build (Resume on from 10^9), hash, e8-verify --container, Structure level (auto), query set, form, checker
#   gen | build | hash | verify | structure | queries | form | check | status | clean   the single steps (arguments as s1b.sh)
#   tests <set>...|all         the test sets: esa-default esa-exhaustive esa-fuzz esa-checked emproto persistent (--network none; needs the full image)
#   export <folder>            copy form.json, form.md and the small result files out of the volume into <folder> and make <folder>/s1b-return.tar.gz (nothing is sent anywhere)
#   shell                      a shell in the container
#   clean-volume               delete the data volume (asks first)
# global options (before the command):
#   --image NAME (seqeron-s1b)  --volume NAME (seqeron-s1b-data)  --memory MiB (container limit; default none)  --cpus N  --detach (run in the background: docker logs -f seqeron-s1b-run)  --dry-run
# The data lives in a NAMED Docker volume inside the VM, never in a macOS folder (a bind mount goes through virtiofs / gRPC FUSE, which is slow and unreliable for the mapped files of the build). Every container
# runs with --network none: the scripts do not need a network and send nothing anywhere.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../../../.." && pwd)
DOCKERFILE=docs/refactoring/suffix-tree-upgrade/analysis/s1b-handoff/Dockerfile
IMAGE=seqeron-s1b; VOLUME=seqeron-s1b-data; MEMORY=; CPUS=; DETACH=0; DRY=0
while [ $# -gt 0 ]; do
  case $1 in
    --image) IMAGE=$2; shift 2 ;;
    --volume) VOLUME=$2; shift 2 ;;
    --memory) MEMORY=$2; shift 2 ;;
    --cpus) CPUS=$2; shift 2 ;;
    --detach) DETACH=1; shift ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) break ;;
  esac
done
CMD=${1:-help}; [ $# -gt 0 ] && shift

die() { echo "s1b-docker: $*" >&2; exit 2; }
run_cmd() { if [ "$DRY" = 1 ]; then echo "+ $*"; else "$@"; fi; }
if [ "$CMD" = help ] || [ "$CMD" = -h ]; then sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 0; fi
if ! command -v docker > /dev/null; then
  [ "$DRY" = 1 ] || die "docker is not installed or not on the PATH (install Docker Desktop and start it)"
fi

# the architecture of the Docker daemon (Docker Desktop on Apple silicon: aarch64 = arm64; on an Intel Mac: x86_64)
DARCH=
command -v docker > /dev/null && DARCH=$(docker info --format '{{.Architecture}}' 2> /dev/null)
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
if [ -n "${S1B_PLATFORM:-}" ] && [ "$S1B_PLATFORM" != "$PLATFORM" ]; then
  die "S1B_PLATFORM=$S1B_PLATFORM differs from the architecture of the Docker daemon ($PLATFORM): that image would run under emulation and its timings would be meaningless. Refused."
fi

need_image() {
  [ "$DRY" = 1 ] && return 0
  local a; a=$(docker image inspect -f '{{.Architecture}}' "$IMAGE" 2> /dev/null) || true
  [ -n "$a" ] || die "image '$IMAGE' not found: run  s1b-docker.sh image  first"
  [ "$a" = "$ARCH" ] || die "image '$IMAGE' is built for $a but the Docker daemon is $ARCH: rebuild with  s1b-docker.sh image  (an image for the other architecture would run under emulation)"
}

# in_container [extra docker args...] -- <s1b.sh arguments>
in_container() {
  local extra=""
  while [ $# -gt 0 ] && [ "$1" != "--" ]; do extra="$extra $1"; shift; done
  [ $# -gt 0 ] && shift
  local base="--init --platform $PLATFORM --network none -v $VOLUME:/data"
  [ -n "$CPUS" ] && base="$base --cpus $CPUS"
  if [ -n "$MEMORY" ]; then base="$base --memory ${MEMORY}m --memory-swap ${MEMORY}m"; fi   # no swap: --memory-swap equal to --memory
  local name=seqeron-s1b-run
  if [ "$DETACH" = 1 ]; then
    if [ "$DRY" = 0 ] && docker ps -a --format '{{.Names}}' | grep -qx "$name"; then
      if docker ps --format '{{.Names}}' | grep -qx "$name"; then die "$name is still running (docker logs -f $name follows it; docker stop $name stops it; the same command again continues where it stopped)"; fi
      docker rm "$name" > /dev/null
    fi
    # shellcheck disable=SC2086
    run_cmd docker run -d --name $name $base $extra "$IMAGE" "$@"
    echo "started in the background: docker logs -f $name   (docker stop $name interrupts it; the same command again resumes, with Resume on)"
  else
    local tty=""; if [ -t 0 ] && [ -t 1 ]; then tty="-it"; fi
    # shellcheck disable=SC2086
    run_cmd docker run --rm $tty $base $extra "$IMAGE" "$@"
  fi
}

case $CMD in
  image)
    TARGET=full
    while [ $# -gt 0 ]; do case $1 in --chain-only) TARGET=chain; shift ;; *) die "unknown option $1" ;; esac; done
    COMMIT=$(git -C "$ROOT" rev-parse HEAD 2> /dev/null || echo unknown)
    if [ "$COMMIT" != unknown ] && [ -n "$(git -C "$ROOT" status --porcelain --untracked-files=no 2> /dev/null)" ]; then echo "note: the checkout has uncommitted changes; the image is built from the files as they are (label: $COMMIT)"; fi
    echo "building $IMAGE (target $TARGET) for $PLATFORM from the repository checkout (the build context is the repository root, filtered by Dockerfile.dockerignore). The first build downloads the base images and, for the full target, the NuGet packages of the test projects; ESTIMATE 10-20 minutes."
    # S1B_BUILD_ARGS: extra arguments for docker build (for example --build-arg HTTPS_PROXY=... behind a proxy); empty by default
    # shellcheck disable=SC2086
    run_cmd docker build --platform "$PLATFORM" --target "$TARGET" -f "$ROOT/$DOCKERFILE" --build-arg "COMMIT=$COMMIT" ${S1B_BUILD_ARGS:-} -t "$IMAGE" "$ROOT" ;;
  preflight|plan|gen|build|hash|verify|structure|queries|chain|form|check|status|clean|tests) need_image; in_container -- "$CMD" "$@" ;;
  export)
    DEST=${1:?usage: s1b-docker.sh export <folder>}
    need_image
    run_cmd mkdir -p "$DEST" || die "cannot create $DEST"
    run_cmd docker run --rm --platform "$PLATFORM" --network none -v "$VOLUME:/data:ro" -v "$DEST:/out" --entrypoint /bin/bash "$IMAGE" -c \
      'cd /data && mkdir -p /tmp/ret && cp -r form.json form.md env preflight.last.txt /tmp/ret/ 2>/dev/null; for d in runs/*/; do mkdir -p "/tmp/ret/$d"; cp "$d"* "/tmp/ret/$d" 2>/dev/null; done; mkdir -p /tmp/ret/tests; cp tests/*.json /tmp/ret/tests/ 2>/dev/null; cp -r /tmp/ret/. /out/ && tar -czf /out/s1b-return.tar.gz -C /tmp/ret . && ls -la /out' ;;
  shell) need_image; in_container --entrypoint /bin/bash -- ;;
  clean-volume)
    printf 'delete the Docker volume %s with everything in it (texts, containers, results)? type the volume name to confirm: ' "$VOLUME"
    read -r ANSWER
    if [ "$ANSWER" = "$VOLUME" ]; then run_cmd docker volume rm "$VOLUME"; else echo "not deleted"; fi ;;
  *) echo "unknown command '$CMD'" >&2; sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 2 ;;
esac
