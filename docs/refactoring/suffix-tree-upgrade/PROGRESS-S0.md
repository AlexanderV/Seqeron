# S0 progress (Safety, honest baseline, perf harness)

Rules: [PROCESS.md](PROCESS.md). Work packages run in order; each = failing test first, fix, 3 suites green (+Parity/Performance where format/perf is touched), commit, push.
Status legend: [ ] todo, [~] in progress, [x] done (with commit hash), [!] blocked (with reason).

| WP | content | status |
|---|---|---|
| 0.0 | Baseline: run the 3 suites + Parity + Performance filters; record counts/durations here. If red: stop, report, do not fix unrelated failures | [x] green, see log |
| 0.1 | Documentation: wrong "approximately 20 bytes per character" (`SuffixTree.cs:76` and any copy in docs/README) -> measured 169 B/char in memory, ~56 B/char persistent file (10 Mbp random DNA, BENCHMARK.md) | [x] `b0c9479b` |
| 0.2 | int32 guards with clear exceptions: unchecked `_nodeCount++` (`PersistentSuffixTreeBuilder.Core.cs:180`), NodeCount overflow with custom `ITextSource` above ~1.07 Gchar, HeapStorageProvider 2 GB cap, `LoadStringInternal` above `int.MaxValue`. Tests must not allocate gigabytes (use a fake `ITextSource` reporting a huge `Length`) | [ ] |
| 0.3 | `PrintTree` capacity `estimatedNodes * 50` int overflow (`SuffixTree.Diagnostics.cs:82`) | [ ] |
| 0.4 | `GenomicAnalyzer.FindRepeats` uses O(n^2) `GetAllSuffixes` (`GenomicAnalyzer.cs:56`): replace with a tree/LCP-based method, identical output; the old implementation stays as a test-local reference oracle (random + adversarial inputs) | [ ] |
| 0.5 | `PersistentSuffixTree.Load`: header-only validation. Add opt-in strict structural validation (offsets/links in range, counts consistent); corrupt-file tests (flipped bytes -> exception, never AV/crash/hang). Default path and cold-load speed (~5 ms) unchanged; Parity green | [ ] |
| 0.6 | MCP `suffix_tree_*` tools rebuild per call with no size limit: documented, configurable limit chosen from measured memory (169 B/char in-memory), clear error, tests | [ ] |
| 0.7 | Perf harness in the repo (existing `apps/SuffixTree.Benchmarks` or a new project; decide by conventions): reproduces the 10 Mbp row of BENCHMARK.md within 15%; not part of default test runs. Sources of the measured harness: `analysis/bench/seqeron/SeqeronBench`, `analysis/bench/go` | [ ] |

## Log
(append newest last: date, WP, commit, test counts, numbers, open issues)
- 2026-11-03, WP 0.0 (baseline, commit: this entry's own commit, docs only). Branch tip before S0 work: `b0922125`. .NET SDK 10.0.403, 4 vCPU. 0 build warnings in all three builds.
  Default `dotnet test -c Release`: SuffixTree.Tests 469 passed (2 s), SuffixTree.Persistent.Tests 560 passed (19 s), SuffixTree.Mcp.Core.Tests 103 passed (0.1 s); 0 failed, 0 skipped.
  `--filter Category=Parity` (Persistent.Tests): 152 passed (10 s). `--filter Category=Performance` (SuffixTree.Tests): 12 passed (3 s).
  Note: the 12 Performance tests are `[Explicit]` (ComplexityGuardTests, PerformanceGuardTests), so they run only with the filter; the default run prints them as "Skipped" lines but reports Skipped: 0.
  Total baseline: 1132 default tests + 152 Parity (subset of the 560) + 12 Performance (explicit, in addition to the 469).
- 2026-11-03, WP 0.1 (commit `b0c9479b`): docs/comment only, no test (not a code defect). Fixed XML doc remark in `SuffixTree.cs` (was "approximately 20 bytes per character"; now 169 B/char random DNA, 190 protein, 193 repetitive, file 56 B/char). Only one copy of the claim existed in the repo; added a measured "2.1 Memory footprint" table to `docs/algorithms/Pattern_Matching/Suffix_Tree.md` and a size line to `SuffixTree.Persistent/README.md`.
  Suites after change: 469 / 560 / 103 passed, Parity 152, Performance 12; 0 warnings.
  Open: `wiki/concepts/suffix-tree.md` has `sources: docs/algorithms/Pattern_Matching/Suffix_Tree.md` and is now flagged stale by `/wiki:lint`; a `/wiki:ingest` refresh was not part of S0 and was not done.
