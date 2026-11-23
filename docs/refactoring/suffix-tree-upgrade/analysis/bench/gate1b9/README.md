# gate1b9: the in-memory build, entry and Open measurements of the S1 gate (work package 1b.9)

One file, outside the shipped code and not part of `Seqeron.sln` (as `readpath1b9g` and `abapp`). It takes the build under test from `-p:EsaDir=<dir>` and carries the assembly name `SuffixTree.Esa.Bench`,
which `SuffixTree.Esa` lets see its internals. It adds nothing to `src/` or `apps/`: the Bench command `build` measures SA + LCP only (with `SaisReference` by default) and no command of the app times the child-table scan
inside the whole build, which the gate of PLAN section S1 asks for (`build_child_s`, PROGRESS-S1 "Measurement protocol").

| command | what |
|---|---|
| `build <text> <name> [sa_vectors.tsv]` | the three calls of `EsaIndex.Build` one by one (prepare, `SaisKernel.Build`, `LcpKernel.Build(Direct)`) and the child-table scan of `EsaIndex.FromArrays` with the default layout of the length: `build_prepare_s`, `build_sais_s`, `build_sa_s` (prepare + SA-IS), `build_lcp_s`, `build_child_s`, `build_total_s` (the gate metric), VmHWM, GC counters, bytes of every in-memory array, SHA-256 of SA and LCP against the libsais vectors |
| `lean <text> <name> [sa_vectors.tsv]` | the same stages with nothing read or counted before the first one (no `/proc` read, no GC counter): the check of whether the instrumentation before the build changes the cold build (1b.9: it does not, beyond the spread of the machine; the interleaved check is `../run_1b9_instr.sh`) |
| `steady <text> <name> [builds]` | the same stages `builds` times in one process, the last three reported |
| `entry <text> <name>` | one call of the public `EsaSuffixTreeFactory.CreateEsa(text)`, cold; the first count after it |
| `open <container> <pattern file>` | the first `SuffixIndexFactory.Open` of the process (cold JIT; nothing before it but the arguments), the first query, 1000 queries, the median of 200 warm opens (index and facade), `Sections` and `Structure` opens; `openp` reads the pattern file before the first Open (the first session of 1b.9 used that order: 9.5 ms against 15 ms for the first Open of a fresh process) |

```
dotnet build src/SuffixTree/Algorithms/SuffixTree.Esa/SuffixTree.Esa.csproj -c Release -o <esa dir>
dotnet build Gate1b9.csproj -c Release -p:EsaDir=<esa dir> -o <out dir>
dotnet <out dir>/SuffixTree.Esa.Bench.dll build /tmp/stbench/data/dna_10m.txt dna_10m docs/refactoring/suffix-tree-upgrade/analysis/sa_vectors.tsv
```

Driven by `../run_1b9_gate.sh` (fresh processes through `lib/benchlib.py`).
