# abcover: the in-process A/B client of the difference cover

Work package 1b.9h, item 4 (debt G(5) of ADR-004 section 5a "1b.9e"): `DifferenceCover.Create` costs 0.27 s per build at X = 993 (11 % of a build of 10^7 symbols, 1.6 % of 10^8). The instrument is outside the shipped
code and not part of `Seqeron.sln` (`abcover.csproj` references nothing; it loads the builds it compares by reflection, as `abapp` does for queries).

## What it measures

| command | what |
|---|---|
| `abcover split <dir> <nominal period> <calls>` | where the time of one `CoverFor` goes: the whole call, the search for the Singer set (`SingerSet(q)`, a private method of the build that has it), the rest (the tables `IndexOf` and `Witness`, the nearest Singer period), and the cost of `Verify` + `IsPerfect` (what a precomputed table that is verified when used costs) |
| `abcover ab <dirA> <dirB> <nominal period> <calls> <blocks>` | the A/B of the two builds in one process, each in its own `AssemblyLoadContext`: A B B A blocks (A and B swap on odd blocks) of `<calls>` calls of `SaBuilder.CoverFor(nominal, false)`; the figure is the median over the blocks of the ratio of the means, with the quartiles and the extremes; the same build on both sides is the null experiment |
| `abcover same <dirA> <dirB> <nominal period>` | the two covers residue by residue (exit code 1 when they differ) |

`<dir>` holds a build of `SuffixTree.Esa.dll` and its dependencies (`dotnet build src/SuffixTree/Algorithms/SuffixTree.Esa/SuffixTree.Esa.csproj -c Release -o <dir>`).

## What it is for and what it is not for

The cover is made once per build; in a fresh process per build (the production chain) a cache that lives in the process gains nothing, which is why the figure that matters is the per-call cost and the fresh-process figure of the
whole build (`sesa build`, or `s1b.sh build`, of the two builds). For the whole build in one process the existing instrument is `sesa resume abhead --old <dir>` (the build before the package against the build of the program, alternating,
Resume off; work package 1b.9e). Always run the null experiment first and report its spread.

## How to run

```
dotnet build abcover.csproj -c Release -o <abcover dir>
dotnet <abcover dir>/abcover.dll split <dir> 1024 20
dotnet <abcover dir>/abcover.dll ab <dirA> <dirA> 1024 20 9        # the null experiment
dotnet <abcover dir>/abcover.dll ab <dirA> <dirB> 1024 20 9
```
