# ST-UPGRADE execution process (S0 + S1)

Authorization (owner, 2026-11-03): "Робимо S0, S1. Зі своєю імплементацією без нативних бінарників. Мегажорстко оптимізуємо.
Не зупиняємось, якщо є хоч щось, що можна пофіксити і покращити. Робимо ітераціями: план, реалізація, тестування,
пошук покращень і фіксів (якщо знайшли - знову з першого пункта). Робимо в окремих контекстах."
Scope of this authorization: **S0 and S1 only** (S1-0 contract/ADR, S1a SA-IS + LCP, S1b child table + navigator + public API).
S2-S8 need a separate yes. Plan and gates: [PLAN.md](PLAN.md). Numbers: [BENCHMARK.md](BENCHMARK.md).

## Hard rules
1. Managed C# only. No native binaries, no P/Invoke, no copied libsais code in the repository. libsais, Go, sdsl may be used only as
   **test oracles / speed references outside the shipped code** (generated vectors, hashes, timings).
2. Existing Ukkonen / persistent code changes only through S0 work packages. Ukkonen is the oracle for S1; do not modify its behaviour.
3. Correctness is never traded for speed. An oracle mismatch stops everything until fixed.
4. No test is skipped, disabled, weakened or deleted to get green. No `#pragma` warning suppression to hide new warnings.
5. Never send the owner's identity (name, e-mail, account) anywhere. No network access is needed for S0/S1.
6. Large artifacts (built binaries, generated 50 Mbp inputs, sdsl/libsais builds) stay out of the repository.

## The loop (each step is a separate agent context)
PLAN -> IMPLEMENT -> TEST -> REVIEW (independent agent, tries to break it, profiles, lists improvements) -> FIX/IMPROVE -> back to PLAN.
- Memory between contexts is **only** the repository: `PROGRESS-S0.md`, `PROGRESS-S1.md`, `git log`. An agent reads its PROGRESS file first and
  updates it last (what was done, numbers, open issues, next step). Orchestrator: the main session.
- The reviewer is never the implementer of the same iteration.
- Final report of every agent: at most 300 words: changes, evidence (test counts, measured numbers), open issues.

## Optimization stopping rule ("do not stop while something can be improved", made measurable)
- An optimization iteration is **accepted** only if it improves the 10 Mbp random-DNA SA-build time by >= 3% (median of 5 runs, quiet machine) and all correctness
  tests/fuzz pass; otherwise it is reverted (kept in the log with its number).
- Stop optimizing after 3 consecutive iterations below 3%, or when build time <= 1.1x libsais (<= 0.36 s for SA alone, <= 0.62 s for SA+LCP at 10 Mbp on this machine). Record why.
  Amendment 2026-11-03 (after 1a.9): the optimized quantity is SA+LCP (the metric the harness reports), so the reference is libsais SA+PLCP+LCP = 0.564 s; the safe baseline is 1.41 s (2.5x). Accept/reject is measured on SA+LCP; the SA-only time is tracked alongside.
  Amendment 2026-11-03 (S1a review MAJ-2): identical binaries differed by 5.9% between sessions, so "median of 5 across sessions" is inside the noise. Accept/reject uses **paired, interleaved A/B runs in one session** (>= 7 pairs, alternate A and B, fresh process each), preceded by an **A/A calibration** with the same binary on both sides; the improvement counts only if the median paired ratio is >= 3% **and** outside the A/A noise band of that session. The libsais reference is re-run in the same session.
  Amendment 2026-11-03 (S1a review MAJ-1): pointer kernels lose the .NET bounds checks, which are currently the only memory-safety net. The first commit of the optimization phase is the safety infrastructure (checked-pointer build, canary buffers, kernel-vs-reference registry with trace equality, paired A/B runner); no `unsafe` kernel is accepted without it.
  Note (S1a review MAJ-3): the size gates (16 B/char index, 24 B/char peak) are defined on random DNA; on highly repetitive texts (a^n, abab, Fibonacci, Thue-Morse at 10 Mbp) the LCP becomes U32 and Full is about 17.4 B/char. That is reported as a documented worst case, not hidden; a denser LCP encoding is a candidate for S1b.
- Machine: 4 vCPU, 15 GB. Timing runs only while no other build/test runs (check `uptime`/`ps`). Report median and min/max.

## Environment
```
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1
dotnet test tests/SuffixTree/SuffixTree.Tests/SuffixTree.Tests.csproj -c Release --nologo
dotnet test tests/SuffixTree/SuffixTree.Persistent.Tests/SuffixTree.Persistent.Tests.csproj -c Release --nologo
dotnet test tests/SuffixTree/SuffixTree.Mcp.Core.Tests/SuffixTree.Mcp.Core.Tests.csproj -c Release --nologo
# format / perf sensitive: add  --filter Category=Parity   and   --filter Category=Performance
```
Benchmark inputs/harness: `analysis/bench` (`HANDOFF.md` section 0 regenerates inputs deterministically; hashes in `MANIFEST.sha256`).

## Git
Branch `claude/stoic-maxwell-0olr2z`. Small commits, message prefix `ST-S0:` or `ST-S1:`, push after every work package.
If the push is rejected: `git fetch origin claude/stoic-maxwell-0olr2z` and `git rebase origin/claude/stoic-maxwell-0olr2z` (own commits only), never force-push.
Commit trailers (required):
```
Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01RvCz2EUn2rHYPaGoqfozS2
```
