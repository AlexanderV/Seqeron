# abapp: the in-process A/B client

The instrument that resolves 1-2 % of throughput, used by the orchestrator's verifications of work packages 1b.9a, 1b.9b and by 1b.9g (ADR-004 section 5a). It is outside the shipped code and not part of `Seqeron.sln`
(`abapp.csproj` references nothing; it loads the builds it compares by reflection).

## Why it exists

Two processes of the same program on this machine differ by 1.4x in throughput in two groups, whatever the variant and the round, so a ratio of 5-10 % cannot be resolved across processes (ADR-004 "1b.9b": the null experiment
of the cross-process instrument has a 95 % interval of 0.74-1.39). Putting the two things to compare **in one process**, alternating in A B B A blocks, cancels the drift and the group of the process: the same build against itself
gives a ratio of 1.00 +- 1-2 % at one thread and +- 1-5 % at four.

## What it compares

`abapp <pattern file> <passes> <threads> <seconds> <specA> <specB>` loads two *specs* in one process, each in its own `AssemblyLoadContext` (so that two builds of `SuffixTree.Esa.dll` can sit side by side), and calls
`CountOccurrences(string)` and `FindAllOccurrences(string)` through delegates:

| spec | what |
|---|---|
| `<dir>\|mem:<text file>` | the in-memory tree (`EsaSuffixTreeFactory.Create`) of the build in `<dir>`, from the text file |
| `<dir>\|file:<container>` | the facade over a SESA container (the file-backed tree, `Checked`); `fileT:` is `Trusted` |
| `<dir>\|idx:<container>` | the index behind that facade: the same search engine without the facade (no reader registration, no disposal check); `idxT:` is `Trusted` |

`<dir>` holds a build of `SuffixTree.Esa.dll` and its dependencies (`dotnet build src/SuffixTree/Algorithms/SuffixTree.Esa/SuffixTree.Esa.csproj -c Release -o <dir>`). Two specs may name the same directory (the null experiment).

A pass is `<seconds>` long (with one thread and `<seconds>` = 0, one sweep of the pattern list); with `<threads>` > 1 the threads are released together and each walks the list from its own offset. A block is A B B A (or B A A B on odd
blocks); the figure of a block is the ratio of the means of its two B passes and its two A passes; the result is the median over `<passes>` blocks, with the quartiles. One line per operation:

```
AB count_present  A 4264873  B 4249459  B/A median-of-paired 0.9582  p25 0.9450 p75 0.9806  n 9
```

`count_present` counts the `P` patterns (present), `locate` lists the `L` patterns, `count_absent` counts the `A` patterns (absent); throughputs are calls per second summed over the threads.

## Pattern files

`genpat.py <text> <alphabet or @text> <len P> <len L> <len A> <n P> <n L> <n A> <seed> <out>` draws the lines `P <substring of the text>`, `L <substring>`, `A <random string over the alphabet>` with SplitMix64 (the benchmark
definition of ADR-004 "1b.9b": 20 / 8 / 24 symbols on DNA, 12 / 5 / 12 on the protein texts for the cache-resident sets of 2,000 / 500 / 2,000 patterns; 200,000 / 20,000 / 100,000 for the large sets).

## How to run

```
dotnet build abapp.csproj -c Release -o <abapp dir>
ABAPP=<abapp dir> ./run_ab.sh <rows file> <processes> <pattern file> <passes> <threads> <seconds> <specA> <specB> [label]
python3 summarize_ab.py <rows file>...
```

**Always run the null experiment first** (specA = specB, same build, same kind of tree, the same thread count) and report its spread: a figure is a pass against a gate only when it clears the gate by more than the spread of the
processes of the null experiment. The gate protocol of 1b.9g: four threads, 9 blocks of passes of 1 s, 3 processes per figure.
