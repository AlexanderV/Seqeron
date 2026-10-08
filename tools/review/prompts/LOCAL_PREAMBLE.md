# Local run notes (prepended to every batch prompt by tools/review/run-batch.sh)

You run as a LOCAL headless Claude Code session on the user's Mac (Apple Silicon), not in a cloud container.
The campaign rules in docs/Validation/review-2026-09/README.md apply unchanged, with these local differences:

- **Worktree:** your working directory is a dedicated git worktree for batch <BATCH> in DETACHED HEAD mode
  (several batches share the branch). Never `git checkout` the branch by name. Use:
  - update:  `git fetch origin claude/stoic-maxwell-0olr2z && git rebase origin/claude/stoic-maxwell-0olr2z`
  - push:    `git push origin HEAD:refs/heads/claude/stoic-maxwell-0olr2z` (on rejection: fetch, rebase, push again)
  - WIP checkpoint (README rule, unchanged):
    `git add -A && git commit -q --no-verify -m "WIP(<BATCH>): <WP> <what> $(date -u +%H:%M)" && git push -q -f origin HEAD:refs/heads/claude/stoic-maxwell-0olr2z-wip-<BATCH> && git reset -q --soft HEAD~1`
- **Leftover local state:** if the worktree already has uncommitted changes when you start, they are from an
  interrupted run of this batch — keep them. First make a WIP checkpoint, then `git stash`, rebase, `git stash pop`,
  and continue that work package. Also check the WIP branch as README says.
- **Environment:** macOS arm64. .NET SDK and Python reference tools are pre-installed (do NOT use apt-get).
  If a Python reference package is missing: `python3 -m pip install <pkg>` (inside the active venv if any).
  Builds/tests are fast here; still follow the test tiers in README (targeted per WP, full fast tier per round, heavy once).
- **No human in the loop:** nobody answers questions. Never end your turn with a question. If you hit a usage/session
  limit, the orchestrator relaunches you later with a resume prompt — your committed work, the WIP branch and this
  worktree's files survive.
- **Concurrency:** another batch may build in its own worktree at the same time; that is fine (separate bin/obj).
  Inside your batch, run subagents STRICTLY SEQUENTIALLY (one build directory).

---

