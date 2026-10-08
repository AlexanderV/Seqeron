#!/usr/bin/env bash
# Stop (pause) a running local review batch. Its worktree and uncommitted files are kept;
# run-batch.sh with a resume prompt continues from there.
#   tools/review/stop-batch.sh <BATCH>
set -euo pipefail
BATCH="${1:?usage: stop-batch.sh <BATCH>}"
ROOT="$(git rev-parse --show-toplevel)"
LOGDIR="${REVIEW_WT_BASE:-$ROOT/../seqeron-review-wt}/_logs"
PIDFILE="$LOGDIR/$BATCH.pid"
[ -f "$PIDFILE" ] || { echo "$BATCH: no pid file" >&2; exit 1; }
PID="$(cat "$PIDFILE")"
if kill -0 "$PID" 2>/dev/null; then
  pkill -INT -P "$PID" 2>/dev/null || true; kill -INT "$PID" 2>/dev/null || true
  for _ in $(seq 1 20); do kill -0 "$PID" 2>/dev/null || break; sleep 1; done
  kill -0 "$PID" 2>/dev/null && { pkill -TERM -P "$PID" 2>/dev/null || true; kill -TERM "$PID" 2>/dev/null || true; }
  echo "$BATCH stopped (pid $PID)"
else
  echo "$BATCH was not running"
fi
