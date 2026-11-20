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

## 5. Memory and parallelism of the fast tier (audit A3-15, B09 F144)

`tests/Seqeron/Seqeron.Genomics.Tests/Parallelization.cs` runs fixtures with `LevelOfParallelism(8)` (the container has 4 cores and 16 GB). Fixtures share no mutable state but they share the machine: the full fast tier once ended in an OOM kill and
once in a flaky 20 s timing test (`CancelledWhileRunning_GenomicAnalyzerScans_StopPromptly`, which waits for the uninterruptible suffix-tree build of a 3 M-base sequence while seven other workers compete). Measured per fixture
(`dotnet test --no-build --filter FullyQualifiedName~.<Class>`, one fixture per process; the peak working set of the whole process tree, the .NET runtime + test host alone is ~270 MB; build 2026-11-20, Debug):

| Fixture | Tests | Wall | Peak working set |
|---|---:|---:|---:|
| `SequenceAligner_Cancellation_Tests` | 23 | 11-17 s | **1.9-2.5 GB** (3 M-base suffix trees, 7 000 x 7 000 matrices) |
| `SequenceAligner_LinearSpace_Tests` | 20 | 38-73 s | **1.8-2.2 GB** (5 600 x 5 600 affine cases) |
| `GenomicAnalyzer_FindOpenReadingFrames_Oracle_Tests` | 18 | 4 s | 0.9 GB (mega-character oracle inputs) |
| `SequenceAligner_TopNLinearSpace_Tests` | 10 | 14-17 s | 0.6-0.7 GB |
| `SequenceAligner_TCoffeeOptions_Tests` | 42 | 47-51 s | 0.9 GB before / **0.54 GB after** the fix below |
| `GenomicAnalyzer_FindCommonRegion_Tests` | 25 | 8 s | 0.5 GB |
| `SequenceAligner_TCoffeeExternal_Tests` | 9 | 6 s | 0.48 GB before / 0.40 GB after |
| 31 other Alignment / Analysis fixtures measured (Mafft, Kalign, MUSCLE, POA, ProbCons, ClustalW, needle / water / stretcher, ksw2, anchors, tandems, similarity sketches ...) | 1 - 1 103 | 3-52 s | 0.27-0.46 GB (= the runtime baseline plus the oracle) |

What was changed: (1) the four fixtures above 0.5 GB plus `SequenceAligner_Cancellation_Tests` are `[NonParallelizable]` (they run alone, so eight 1-2 GB fixtures can no longer coincide), the rule being written into `Parallelization.cs`;
(2) one parsed copy of the oracle per fixture: `SequenceAligner_TCoffeeOptions_Tests.Oracle()` and `SequenceAligner_TCoffeeExternal_Tests.Oracle()` re-read and re-parsed the whole embedded TSV (6.0 / 3.2 MB as UTF-16 strings, about 20 000 rows) in
EVERY test, and twelve `GenomicAnalyzer_*Oracle*` fixtures called `ReadOracle()` two to six times: all now parse once (`Lazy<>`); (3) the timing test no longer has a latency bound that the uninterruptible build can exceed
(fixture runs alone; the limit of the two tree-building calls is 120 s). **Rule for new fixtures:** more than ~0.5 GB working set -> `[NonParallelizable]`; the measurement command above is the way to find out.

**Embedded test data (decision: keep plain TSV).** `TestData/` is 76 MB in 133 files (the 12 largest: `tcoffee_program_options_oracle` 6.0 MB, `emboss_aformat_oracle` 5.4, `orffinder_cli_oracle` 4.7, `mhcflurry_single_net.bin` 4.4 (B26), `ani_ci_oracle` 3.2,
`tcoffee_external_oracle` 3.2, `biopython_msa_oracle` 2.8, `clustalw_commands_oracle` 2.4, `poa_oracle` 2.2, `clustalw_output_oracle` 2.0, `mash_dist_oracle` 1.8, `emboss_multi_oracle` 1.8); the hand-off packages are 21 MB; `.git` is 210 MB (pack 91 MiB).
Compressing the TSVs as `.gz` resources was weighed and rejected: git already stores every blob deflated (the pack is 91 MiB; `gzip -9` of the three largest files gives 4x / 16x / 4x, the same ratio git achieves), so the repository would not shrink, and a binary `.gz`
re-generated by a script defeats delta compression (every regeneration is a full new blob); the 94 fixtures read their file through their own `GetManifestResourceStream` call (no common helper), so the change would touch all of them; and a compressed
oracle cannot be reviewed with `git diff` / `grep` row by row, which is the point of committing the oracles. The cost that mattered - one string copy per read - is removed by the per-fixture caches in (2).

