# readpath1b9g: three small tools of work package 1b.9g

Outside the shipped code and not part of `Seqeron.sln`. Each is one file, takes the build under test from `-p:EsaDir=<dir>` (a directory with `SuffixTree.Esa.dll`, `SuffixTree.Core.dll` and `SuffixTree.Persistent.dll`, for example the
output of `dotnet build src/SuffixTree/Algorithms/SuffixTree.Esa/SuffixTree.Esa.csproj -c Release -o <dir>`), and carries the assembly name `SuffixTree.Esa.Bench`, which `SuffixTree.Esa` lets see its internals. Build into a fresh output
directory for every build under test, or copy the three dlls of the build into the output directory afterwards (the build keeps the first references it resolved).

| tool | what it measures | used for |
|---|---|---|
| `FacadeCost` | `facade`: the cost of one call of the facade (`IsEmpty`: the registration of the call, one comparison, the exit) on a tree over a container and on an in-memory tree, 1 and 4 threads on one tree; `gate` (`-p:Gate=true`, a build of 1b.9g or later): the pieces of the new protocol and the striped count of 1b.8 emulated; `barrier`: the time of `Interlocked.MemoryBarrierProcessWide()` with 0, 1, 2 and 4 other threads spinning | the cost of the registration of a call before and after the protocol (ADR-004 "1b.9g") |
| `Decompose` | one thread, index level (no facade), the in-memory index and the index over the container in one process, in blocks of 0.25 s in alternating order: the count of present patterns, of absent patterns, locate, the descent alone (`TryMatch`) for the locate patterns, the materialisation of the list (`Sa.ToArray`, `CopyTo`, `new int[n]`) | the profile of item 2: where the cache-resident cost of the mapped path is (the descent, not the list) |
| `ColdStart` | a fresh process: the time of the first blocks of queries after the open, per block of N, for count (`P`) or locate (`L`) | the warm-up cost of leaving the mapped descent to tiered compilation |

```
dotnet build FacadeCost/FacadeCost.csproj -c Release -p:EsaDir=<dir> -o <out>             # add -p:Gate=true for the gate mode
dotnet <out>/SuffixTree.Esa.Bench.dll facade <container> <text file>      # or: gate | barrier
dotnet build Decompose/Decompose.csproj -c Release -p:EsaDir=<dir> -o <out>
dotnet <out>/SuffixTree.Esa.Bench.dll <pattern file> <container> <text file> [C|T]         # C = Checked (default), T = Trusted
dotnet build ColdStart/ColdStart.csproj -c Release -p:EsaDir=<dir> -o <out>
dotnet <out>/SuffixTree.Esa.Bench.dll <container> <pattern file> <P|L> <blocks> <block size>
```

The pattern files are those of `../abapp/genpat.py` (lines `P <pattern>`, `L <pattern>`, `A <pattern>`). `Decompose` and `ColdStart` use the internal types of the build, so they build against a `SuffixTree.Esa.dll` that names
`SuffixTree.Esa.Bench` in its `InternalsVisibleTo` (every build of the S1b line does).
