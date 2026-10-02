# Running the review campaign on your own Mac

Same campaign, same rules (README.md), same prompts, same queue (QUEUE.md). Only *where* sessions run changes:
batches run as local headless `claude -p` sessions, each in its own git worktree, and the orchestrator is a local
Claude Code session you can watch from the phone (Remote Control).

| Cloud (before) | Local (now) |
|---|---|
| orchestrator = cloud session | orchestrator = local `claude` session in the repo (Remote Control on) |
| `create_session` | `tools/review/run-batch.sh <BATCH> <prompt-file>` |
| `get_session` | `tools/review/status.sh [BATCH]` (process alive? last messages, final result, cost, branch/WIP tips) |
| `interrupt_session` | `tools/review/stop-batch.sh <BATCH>` (worktree + uncommitted files kept) |
| cloud trigger every 45 min | `/loop 45m …` in the orchestrator session |
| container lost at idle | worktree stays on disk; WIP branch still used as backup |

## 1. One-time setup (macOS, Apple Silicon)

```bash
# Toolchains
brew install --cask dotnet-sdk          # must provide .NET 10:  dotnet --list-sdks
brew install git gh python@3.12
python3.12 -m venv ~/.venvs/seqeron-review
source ~/.venvs/seqeron-review/bin/activate
pip install biopython primer3-py ViennaRNA scikit-bio scipy numpy
export DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1

# Claude Code CLI (same account as claude.ai)
curl -fsSL https://claude.ai/install.sh | bash
claude auth login

# GitHub push access for git (sessions push to the campaign branch)
gh auth login && gh auth setup-git

# Repo
mkdir -p ~/src && cd ~/src
git clone https://github.com/AlexanderV/Seqeron && cd Seqeron
git checkout claude/stoic-maxwell-0olr2z
dotnet build Seqeron.sln -c Debug      # 0 errors expected
```

Add to `~/.zshrc` so every batch gets the reference tools:

```bash
export REVIEW_VENV=~/.venvs/seqeron-review
export DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
```

Worktrees and logs go to `~/src/seqeron-review-wt/` (override with `REVIEW_WT_BASE`).

**Permissions.** Batches run with `--permission-mode auto` (same as the cloud sessions). If your CLI rejects `auto`
in headless mode, run with `REVIEW_PERMISSION_MODE=bypassPermissions` instead — that skips all permission prompts,
so only do it in this dedicated checkout.

**Keep the Mac awake.** `run-batch.sh` wraps each session in `caffeinate -is`, which prevents idle sleep while the
batch runs. Keep the Mac on power; closing the lid still sleeps a MacBook unless it is in clamshell mode with an
external display.

## 2. Manual control

```bash
cd ~/src/Seqeron
tools/review/run-batch.sh B07 tools/review/prompts/B07.md   # start / resume a batch
tools/review/status.sh                                     # all batches
tools/review/status.sh B07                                 # one batch, last 5 messages
tools/review/stop-batch.sh B08                             # pause a batch
tail -f ../seqeron-review-wt/_logs/B07.latest.jsonl        # raw stream
```

A session that hit a usage limit exits; `status.sh` shows `EXITED` with an error result mentioning the limit.
Relaunch it with the same command after the reset — it continues from the worktree, the WIP branch and Bxx.md.

## 3. The orchestrator session

```bash
cd ~/src/Seqeron && claude remote-control      # or plain `claude`; Remote Control lets you follow from the phone
```

First message to the orchestrator:

> You are the orchestrator of the Seqeron review campaign. Read docs/Validation/review-2026-09/README.md,
> QUEUE.md and LOCAL_RUN.md. You do not edit library code yourself; you start, watch, pause and finish batches with
> tools/review/{run-batch,status,stop-batch}.sh, keep QUEUE.md current (commit + push it after each change), and
> brief me in Russian. Never start anything the queue's standing user rules forbid. Then run the check-in below once
> and start `/loop 45m` with it.

Check-in text for `/loop 45m`:

> Check-in (review campaign). Read QUEUE.md (standing rules + state). Run tools/review/status.sh.
> For every batch: RUNNING → note progress (new commits, WIP tip). EXITED with a usage/session-limit error → if the
> reset time has passed, relaunch it with its prompt (a resume prompt: read Bxx.md incl. "Audit worklist" + git log,
> cherry-pick a pending WIP tip first); otherwise wait. EXITED otherwise → read its Bxx.md "## Leftovers" yourself:
> any DOABLE item (LIMITATIONS proposals count) → write a FINISHER prompt from tools/review/prompts/TEMPLATE.md and
> launch it; zero DOABLE → mark it Done in QUEUE.md. If a slot is free (max 2) and the standing rules allow launches,
> take the next pending batch: build its prompt from TEMPLATE.md + the README batch table (units, owned files),
> add a RESUME block if Bxx.md already has progress, list the known cross-batch requests to it (grep other B*.md),
> save it as tools/review/prompts/<BATCH>.md, commit, and run-batch.sh it. Update QUEUE.md, commit + push,
> and brief me in Russian in a few lines. Never launch more than 2 batches; never launch while I asked to pause.

## 4. Moving B07 from the cloud

1. Tell the cloud orchestrator "переезжаем": it stops the cloud B07 right after its next WIP checkpoint and stops its
   cloud check-ins.
2. Locally: `tools/review/run-batch.sh B07 tools/review/prompts/B07.md` — the local session picks up the WIP tip and
   the round-3 worklist and continues.
3. Start the local orchestrator (section 3).
