# S0 progress (Safety, honest baseline, perf harness)

Rules: [PROCESS.md](PROCESS.md). Work packages run in order; each = failing test first, fix, 3 suites green (+Parity/Performance where format/perf is touched), commit, push.
Status legend: [ ] todo, [~] in progress, [x] done (with commit hash), [!] blocked (with reason).

| WP | content | status |
|---|---|---|
| 0.0 | Baseline: run the 3 suites + Parity + Performance filters; record counts/durations here. If red: stop, report, do not fix unrelated failures | [ ] |
| 0.1 | Documentation: wrong "approximately 20 bytes per character" (`SuffixTree.cs:76` and any copy in docs/README) -> measured 169 B/char in memory, ~56 B/char persistent file (10 Mbp random DNA, BENCHMARK.md) | [ ] |
| 0.2 | int32 guards with clear exceptions: unchecked `_nodeCount++` (`PersistentSuffixTreeBuilder.Core.cs:180`), NodeCount overflow with custom `ITextSource` above ~1.07 Gchar, HeapStorageProvider 2 GB cap, `LoadStringInternal` above `int.MaxValue`. Tests must not allocate gigabytes (use a fake `ITextSource` reporting a huge `Length`) | [ ] |
| 0.3 | `PrintTree` capacity `estimatedNodes * 50` int overflow (`SuffixTree.Diagnostics.cs:82`) | [ ] |
| 0.4 | `GenomicAnalyzer.FindRepeats` uses O(n^2) `GetAllSuffixes` (`GenomicAnalyzer.cs:56`): replace with a tree/LCP-based method, identical output; the old implementation stays as a test-local reference oracle (random + adversarial inputs) | [ ] |
| 0.5 | `PersistentSuffixTree.Load`: header-only validation. Add opt-in strict structural validation (offsets/links in range, counts consistent); corrupt-file tests (flipped bytes -> exception, never AV/crash/hang). Default path and cold-load speed (~5 ms) unchanged; Parity green | [ ] |
| 0.6 | MCP `suffix_tree_*` tools rebuild per call with no size limit: documented, configurable limit chosen from measured memory (169 B/char in-memory), clear error, tests | [ ] |
| 0.7 | Perf harness in the repo (existing `apps/SuffixTree.Benchmarks` or a new project; decide by conventions): reproduces the 10 Mbp row of BENCHMARK.md within 15%; not part of default test runs. Sources of the measured harness: `analysis/bench/seqeron/SeqeronBench`, `analysis/bench/go` | [ ] |

## Log
(append newest last: date, WP, commit, test counts, numbers, open issues)
