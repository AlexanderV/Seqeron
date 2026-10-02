#!/usr/bin/env bash
# Launch one review batch as a local headless Claude Code session in its own git worktree.
#   tools/review/run-batch.sh <BATCH> <prompt-file>
# Example: tools/review/run-batch.sh B07 tools/review/prompts/B07.md
# Env overrides: REVIEW_WT_BASE (worktree root), REVIEW_PERMISSION_MODE (default: auto),
#                REVIEW_VENV (python venv with reference tools, activated if present).
set -euo pipefail

BATCH="${1:?usage: run-batch.sh <BATCH> <prompt-file>}"
PROMPT="${2:?usage: run-batch.sh <BATCH> <prompt-file>}"
BRANCH="claude/stoic-maxwell-0olr2z"
ROOT="$(git rev-parse --show-toplevel)"
WT_BASE="${REVIEW_WT_BASE:-$ROOT/../seqeron-review-wt}"
WT="$WT_BASE/$BATCH"
LOGDIR="$WT_BASE/_logs"
PIDFILE="$LOGDIR/$BATCH.pid"
mkdir -p "$LOGDIR"

[ -f "$PROMPT" ] || { echo "prompt file not found: $PROMPT" >&2; exit 1; }
PROMPT="$(cd "$(dirname "$PROMPT")" && pwd)/$(basename "$PROMPT")"

if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
  echo "$BATCH is already running (pid $(cat "$PIDFILE"))" >&2; exit 1
fi

git -C "$ROOT" fetch -q origin "$BRANCH"
if [ ! -d "$WT" ]; then
  # Detached worktree: several batches can work on the same branch at once.
  git -C "$ROOT" worktree add -q --detach "$WT" "origin/$BRANCH"
  echo "created worktree $WT at origin/$BRANCH"
else
  # Keep whatever an interrupted run left behind (uncommitted work survives locally).
  echo "reusing worktree $WT ($(git -C "$WT" status --porcelain | wc -l | tr -d ' ') uncommitted paths kept)"
fi

if [ -n "${REVIEW_VENV:-}" ] && [ -f "$REVIEW_VENV/bin/activate" ]; then
  # shellcheck disable=SC1091
  . "$REVIEW_VENV/bin/activate"
fi

STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="$LOGDIR/$BATCH-$STAMP.jsonl"
FULL_PROMPT="$LOGDIR/$BATCH-$STAMP.prompt.md"
cat "$ROOT/tools/review/prompts/LOCAL_PREAMBLE.md" "$PROMPT" | sed "s/<BATCH>/$BATCH/g" > "$FULL_PROMPT"

KEEPAWAKE=()
command -v caffeinate >/dev/null 2>&1 && KEEPAWAKE=(caffeinate -is)   # macOS: no idle/system sleep while running

cd "$WT"
nohup ${KEEPAWAKE[@]+"${KEEPAWAKE[@]}"} claude -p \
  --permission-mode "${REVIEW_PERMISSION_MODE:-auto}" \
  --output-format stream-json --verbose \
  < "$FULL_PROMPT" > "$LOG" 2>&1 &
echo $! > "$PIDFILE"
ln -sf "$LOG" "$LOGDIR/$BATCH.latest.jsonl"
echo "started $BATCH pid $(cat "$PIDFILE")"
echo "  worktree: $WT"
echo "  log:      $LOG"
