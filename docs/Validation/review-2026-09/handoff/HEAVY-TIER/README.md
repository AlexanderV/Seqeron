# STATUS: RUN 1 DONE (f419048e: 7993/0); RUN 2 REQUIRED on the final commit (code added since) - TEMPORARY hand-off package HEAVY-TIER

**Run 1 (user, 2026-11-20, aarch64, 18 cores, .NET 10.0.404, commit f419048e):** `genomics_heavy rc=0 seconds=240 :: Passed! - Failed: 0, Passed: 7993, Skipped: 0, Total: 7993, Duration: 3 m 57 s`; returned archive: `returned/heavy-tier-out.zip`. This was an early check; commits after f419048e (megamerger, nucmer clusters, minimap2 chaining, dotmatcher and every audit round-4 item) are NOT covered. Run 2 = this same package on the final commit, exactly once at the very end (campaign README rule); the package tests the commit that is checked out, so the next run needs no change to it.

Runs the .NET test tiers of this repository inside a Docker container on **your** machine (normally faster than the cloud sandbox, which has a shared
memory limit and has already OOM-killed test runs). Default = the **heavy tier** (`Fuzzing/`, `Properties/`, `Metamorphic/` of
`Seqeron.Genomics.Tests`) that has not been run yet in this batch and is the longest step; `MODE=ALL` also runs the full fast tier and every MCP /
Integration / SuffixTree.Mcp test project. Uses the official `mcr.microsoft.com/dotnet/sdk:10.0` image (native arm64 on Apple Silicon: no emulation).
Nothing is downloaded except NuGet packages and the apt packages python3/rsync/zip; no account data, no secrets.

## Commands (run from your clone; the package tests exactly the commit you have checked out)

```bash
git fetch origin claude/busy-sagan-8ztcyb
git checkout claude/busy-sagan-8ztcyb
git pull origin claude/busy-sagan-8ztcyb
cd docs/Validation/review-2026-09/handoff/HEAVY-TIER
./run.sh                 # heavy tier only (default)
# optional variants:
#   MODE=ALL ./run.sh    heavy + fast + MCP + integration (everything)
#   MODE=MCP ./run.sh    only the MCP projects
#   WORKERS=4 ./run.sh   fewer parallel test workers (use if Docker Desktop has < 12 GB RAM; fixtures peak at 2-2.5 GB each)
```

Docker Desktop -> Settings -> Resources: give it all CPUs and at least 12 GB of memory (8 GB works with `WORKERS=2`).

## What to send back

`returned/heavy-tier-out.zip` (contains `out/summary.txt` with pass/fail counts and durations per project, `out/*.log`, `out/*.failed.txt` with every
failing test name, `out/trx/*.trx`, `out/ENV.txt`, `out/build.log`). Put it in `returned/` (`git add -f`, commit, push to `claude/busy-sagan-8ztcyb`)
or attach it to the session. A run with failures is just as useful: the failing test names tell the next session what to fix.

## What the next session does with it

Fixes every failing heavy-tier test (updates a test only with a sourced reason - never weakens an invariant), adds heavy-tier tests for the new
behaviour of the batch, and asks for ONE final run of this same package on the final commit (README rule: heavy tier exactly once at the end; this
early run only finds problems in advance).
