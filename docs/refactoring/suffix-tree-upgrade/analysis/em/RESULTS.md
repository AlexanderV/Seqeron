# ADR-004 prototypes E0-E3: measured results

Prototype code: `apps/SuffixTree.EmProto` (managed C#, no native code), tests `tests/SuffixTree/SuffixTree.EmProto.Tests`. Raw rows: the TSV files in this directory
(`R <TAB> experiment <TAB> case <TAB> metric <TAB> value <TAB> unit <TAB> notes`). Analysis in the ADR: [ADR-004 section 5a](../../adr/ADR-004-external-memory-construction.md).
Every number below is **measured here** (4 vCPU Intel Xeon, 15 GB, kernel 6.18.44, .NET 10.0.14, Release, workstation non-concurrent GC, nothing else running; build servers shut down)
unless it is marked estimated or literature. Inputs: `/tmp/stbench/data` (outside the repository).

## E0 harness (2026-11-09)

Command: `dotnet apps/SuffixTree.EmProto/bin/Release/net10.0/SuffixTree.EmProto.dll e0 --dir <scratch> --out e0.tsv` (raw rows: `e0.tsv`).

| gate | measured | result |
|---|---|---|
| RssAnon sampler (`/proc/self/status`, 10 ms) overhead < 2 % | 0.64 % (median sampled/unsampled ratio 1.0064 of 15 interleaved pairs of `SaisKernel` on 10 M random DNA, 95 % interval 0.994-1.052); A/A median 0.9972, A/A noise band [0.9917, 1.0083], so the overhead is inside the noise of this session | PASS |
| sampler reads the right quantity | 200 MiB (209.7 MB) of touched native memory raised the sampled RssAnon maximum by 212.1 MB | PASS |
| explicit I/O counters within 1 % of `/proc/self/io` | writes: 301,989,888 bytes counted, `wchar` delta 301,989,888 (difference 0 %); reads: 335,544,320 counted, `rchar` delta 335,702,685 (0.047 %, the difference is the sampler's own reads of `/proc/self/status`) | PASS |
| `ArenaBudget` on `NativeArena` | cap enforced (a block that would exceed it is refused and reserves nothing), `ReleaseTo` and `Dispose` return the bytes, shared by arenas of several threads (unit tests) | PASS |

Notes. The counters cover explicit `pread`/`pwrite` calls (class `CountedFile`); accesses through memory mappings are not explicit calls and are reported separately in E3 (`RssFile`, block-layer
`read_bytes` of `/proc/self/io`, `Dirty` of `/proc/meminfo`). The `ArenaBudget` is a separate class that wraps a `NativeArena` (`BudgetedArena`): `NativeArena` itself is not modified.
