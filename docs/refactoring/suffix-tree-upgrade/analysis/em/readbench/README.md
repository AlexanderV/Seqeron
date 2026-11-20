# ReadBench: the read-path benchmark of the suffix-tree upgrade (S1b)

This folder holds the measurement harness of work packages 1b.9a (the diagnosis of the file-backed read path) and 1b.9b (the reconciliation of its two instruments). It is **outside the shipped code**: `ReadBench.csproj` is not part of
`Seqeron.sln` (checked in 1b.9b: the solution file does not name it), it is built only by the scripts of `analysis/em/` (`dotnet build ReadBench.csproj -c Release -o <dir>`), and it must stay that way: it carries the assembly name
`SuffixTree.Esa.Bench` because `SuffixTree.Esa` lets that name see its internals (the experiments of the diagnosis need them), and the files `Experiments*.cs` are compiled only with `-p:Experiments=true`.

## Benchmark definition (work package 1b.9b)

One definition, used by every number that 1b.9b reports. The two trees compared are the **in-memory tree** (`EsaSuffixTreeFactory.Create(text)`) and the **file-backed tree** in `Checked` mode over the container that the builder of
the head writes for the same text (`EsaSuffixTreeFactory.Open(path)` at the `Header` level; `Trusted` is reported for the record in the in-process instrument).

### Texts (10 Mbp each, `n = 10,000,000` symbols, one byte per symbol; `ReadBench gen <kind> 10000000 20261120 <file>`)

| name | symbols | generator |
|---|---|---|
| `dna` | 4 (ACGT) | the S0 input `dna_10m` (`/tmp/stbench/data/dna_10m.txt`), as in 1b.9a |
| `sigma5` | 5 (ACGTN) | uniform random over the five symbols |
| `protu` | 20 | uniform random over the 20 amino-acid letters |
| `protr` | 20 | a proteome-like text: Swiss-Prot composition (L 9.7 %, A 8.3 %, G 7.1 %, V 6.9 %, E 6.8 % ... W 1.1 %, C 1.4 %), proteins of length 60 + exponential (mean 360) each starting with M, 35 % of the segments copies of an earlier stretch (40-300 residues) with a substitution rate of 2, 10, 30 or 50 % and indels at 1 % (paralogous domains), 4 % low-complexity runs (a residue 5-30 times or a unit of 2-6 residues repeated 3-15 times), otherwise 10-60 fresh residues |
| `text90` | 90 | a natural-language-like text: 20,000 random words of English letter frequencies drawn with Zipf probability 1 / rank, sentences of 5-20 words, capitals, punctuation, numbers (5 % of the tokens) and rare marks (1 %) |
| `bytes256` | 256 | uniform random bytes |

The container of each text is built by `SesaBuildCostTests` (variant `container`, the real geometry, with PSI), codec 2 (tagged child table), as the builder's default.

### Patterns (drawn with SplitMix64 from the text, seed 1234, the same patterns for every variant and every process of a text)

| kind | what | length | number | operation |
|---|---|---|---|---|
| present long | substrings of the text at uniformly random positions | 20 | 200,000 | `CountOccurrences` |
| present short | substrings of the text at uniformly random positions | 8 | 20,000 | `FindAllOccurrences` (locate: the list of positions is materialised) |
| absent | uniformly random strings over the alphabet of the text | 24 | 100,000 | `CountOccurrences` |

`RB_L20`, `RB_L8`, `RB_L24` and `RB_N20`, `RB_N8`, `RB_N24` change the lengths and the numbers (a sensitivity run on the protein texts uses 12 / 5 / 12 and 2,000 / 500 / 2,000, the protocol of the orchestrator's check of 1b.9a).

### Instruments (always both)

**(a) Cross-process.** One tree per process: the process builds or opens its tree, warms the JIT by time (1.5 s of count, 0.7 s of locate, 1.0 s of absent patterns), runs 7 passes of the three pattern sets and reports the median of each.
The variants run in a **rotated order** from round to round (`run_1b9a_matrix.sh`: round r starts with the (r-1) mod k-th variant), the page cache is warm, `DOTNET_GCHeapHardLimit` is **unset** (the in-memory tree of 10 Mbp needs more
than 96 MiB), at least **21 processes per variant** (`ROUNDS`). Reported: the ratio of the medians over all processes (file-backed / in-memory), the median over the rounds of the per-round ratio (paired by round), and a 95 % bootstrap
interval of the first (the processes of both variants resampled with replacement).

**(b) In process.** Both trees (and the `Trusted` tree) live in one process; every pass runs all of them, the order rotated from pass to pass; 11 passes after the warm-up, 6 processes; reported is the median over the processes of the
median over the passes of the per-pass throughput ratio to the in-memory tree (and the range over the processes).

Why two: the spread between processes of one variant on this machine is up to 2x and bimodal, so (a) resolves about 5 % only with many processes; (b) resolves 1-2 % but puts the trees in one process, where they share the caches, the
heap and the page tables. A disagreement between them is a result, and is explained by an isolated experiment (1b.9b, ADR-004 section 5a).

### Files

`run_1b9b_bench.sh` (the driver: `harness`, `texts`, `containers`, `cross <text>...`, `inproc <text>...`, `report`), `summarize_1b9b.py` (the table), `run_1b9a_matrix.sh` (the cross-process rounds), `Program.cs` (the harness), `Generate.cs` (the
texts), `Experiments*.cs` (the isolated experiments of 1b.9a and 1b.9b, `-p:Experiments=true`).
