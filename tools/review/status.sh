#!/usr/bin/env bash
# Show the state of every local review batch: running?, last activity, last result, cost, branch tips.
#   tools/review/status.sh            # all batches
#   tools/review/status.sh B07        # one batch, with the last assistant messages
set -euo pipefail
BRANCH="claude/stoic-maxwell-0olr2z"
ROOT="$(git rev-parse --show-toplevel)"
WT_BASE="${REVIEW_WT_BASE:-$ROOT/../seqeron-review-wt}"
LOGDIR="$WT_BASE/_logs"

git -C "$ROOT" fetch -q origin 2>/dev/null || true
echo "== branch tips =="
git -C "$ROOT" log -1 --format='main  %h %cd %s' --date=format:'%m-%d %H:%M' "origin/$BRANCH" | cut -c1-150
for wip in $(git -C "$ROOT" for-each-ref --format='%(refname:short)' "refs/remotes/origin/$BRANCH-wip-*"); do
  git -C "$ROOT" log -1 --format="${wip##*-wip-}   %h %cd %s" --date=format:'%m-%d %H:%M' "$wip" | cut -c1-150
done

[ -d "$LOGDIR" ] || { echo "no local batches started yet ($LOGDIR)"; exit 0; }
echo "== batches =="
for pidf in "$LOGDIR"/*.pid; do
  [ -e "$pidf" ] || continue
  b="$(basename "$pidf" .pid)"
  [ $# -gt 0 ] && [ "$1" != "$b" ] && continue
  pid="$(cat "$pidf")"
  if kill -0 "$pid" 2>/dev/null; then state="RUNNING"; else state="EXITED"; fi
  log="$LOGDIR/$b.latest.jsonl"
  python3 - "$b" "$state" "$log" "$([ $# -gt 0 ] && echo detail || echo brief)" <<'PY'
import json, os, sys, time
b, state, log, mode = sys.argv[1:5]
last_ts = time.strftime('%m-%d %H:%M', time.localtime(os.path.getmtime(log))) if os.path.exists(log) else '-'
texts, result = [], None
if os.path.exists(log):
    for line in open(log, errors='replace'):
        try: o = json.loads(line)
        except Exception:
            if line.strip(): texts.append(line.strip()[:300])
            continue
        if o.get('type') == 'assistant':
            for c in o.get('message', {}).get('content', []):
                if c.get('type') == 'text' and c['text'].strip(): texts.append(c['text'].strip())
        elif o.get('type') == 'result':
            result = o
line = f"{b:6} {state:8} last-log {last_ts}"
if result:
    line += f"  result={result.get('subtype')}  is_error={result.get('is_error')}  cost=${result.get('total_cost_usd', 0):.2f}"
print(line)
if result and result.get('result'):
    print('   final:', str(result['result']).strip().replace('\n', ' ')[:400])
n = 5 if mode == 'detail' else 1
for t in texts[-n:]:
    print('   >', t.replace('\n', ' ')[:400])
PY
done
