# Large: tests of texts above 2^31 symbols (work package S2.1)

Everything that makes a text of 2^31, 2^32, 2^33 or 2^40 symbols testable on a 15 GB RAM / 16 GB disk machine without native code (ADR-005 section 3.7). Nothing here is a product class: the generators, the sparse container writer and the 5-byte re-encoder are test code.

## What is where

| file | what |
|---|---|
| `ClosedFormText.cs` | the contract of a closed-form text: `SuffixAt`, `RankOf`, `Lcp`, `CodeAt`, `RangeOf` / `RangeOfRun` / `RangeOfFactor` (count, first and last rank of a pattern), `Down` / `Up` / `Next` (the child table at the semantic level, the values of the reference layout), `PackedCell` (the cell of CHILD codec 1), `InternalCount`, `NodeCount`, `MaxDepth`, `LrsDepth`, `LrsRank`, `LrsPosition`, `LcpSum`; every argument and result is `long` |
| `PowerOfA.cs`, `PeriodicAkB.cs`, `DeBruijnRepeat.cs` | the families `a^n`; `(a^k b)*` truncated to any n (`(ab)^m` is k = 1); `u^m` of a de Bruijn word (no child table: it is the engine's scan over the generated LCP) |
| `ClosedFormChecks.cs` | the oracles the closed forms are compared with: the definitions (n up to 2,000), the in-memory index of S1 (n up to 10^6), the periodicity of the text (`PeriodicOracle`: n = 2^31 and above) |
| `Crc32CCombine.cs` | CRC-32C of N zero bytes and of a section with a few patches by the algebra of the polynomial (square and multiply), not by a pass |
| `SparseContainer.cs` | writes a SESA container of any n whose big sections are holes: valid header page and table at POS_WIDTH 4 and 5 (ADR-005 section 6.1), section CRCs by the combine, the patches a test chooses (`EndPatches` of a closed form) |
| `LargeDisk.cs` | the disk rules: `$SEQERON_LARGE_DIR` or a subdirectory of the temporary directory, the sparse probe, 64 MiB of dense bytes at most, 200 GB nominal at most |
| `Sesa5Encoder.cs`, `Sesa5Image.cs` | the test-only re-encoder of a container of POS_WIDTH 4 into POS_WIDTH 5 and the reader that knows only ADR-005 section 6.1 |
| `LargeFixtureBase.cs` | the base of the `Large` fixtures: the disk guard, the clean-up, the `LARGE-REPORT` line |
| `*Tests.cs` | default set (small n, fast, no disk beyond the temporary directory): `ClosedFormSelfTests`, `Crc32CCombineTests`, `SparseContainerTests`, `Sesa5EncoderTests`; explicit (`Large`): `LargeCrcTests`, `LargeSparseContainerTests`, `LargeRedTests` |

## Categories and how to run

* **Default run** (`dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release`): only the self-tests of the closed forms (3 oracles), of the CRC combine, of the sparse writer at small n and of the encoder (round trips at n <= 10^6). Nothing explicit runs.
* **`Large`** (explicit, never in the default run): `dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter Category=Large` (about 1 minute: sparse files of up to 78 GB nominal, the `Sections` pass over 26 GB of holes, 8.6 GB of zeros for the CRC). One class: `--filter FullyQualifiedName~LargeRedTests`. Add `--logger "console;verbosity=detailed"` to see the report lines.
  * `LargeRedTests` are the red suites of S2.3 - S2.8: they open a sparse container of 2^31 + 10 / 2^32 + 10 symbols through the entry points that exist at the S1 head and assert the stage-S2 behaviour; the property `TurnsGreenAt` names the package that makes each pass. The tests marked `Kind=guard` pass today and must keep passing.
* **`Width5`** (explicit): the container conformance suites (`IndexConformanceTests`, `IndexAliasingTests`, `IndexAliasingPropertyTests`) again over `Sesa5Encoder(container)`: `--filter Category=Width5`. They are red until S2.7 (the S1 reader refuses POS_WIDTH 5); S2.7 turns them into the default set by moving them into `IndexEngineFactories.All`.

## Disk safety (binding)

The volume has about 17 GB free and 252 GB nominal. A sparse file is made only after a probe (`SetLength` of 1 GiB, then `stat` blocks, which must be below 1 MiB); on a filesystem that does not make sparse files the whole `Large` fixture is inconclusive with the message, and nothing is written. No dense region above 64 MiB is ever written, the nominal size of one container is capped at 200 GB, at least 2 GB must stay free, every file a test made is deleted in the tear-down (a failed write deletes its own file), and a sparse file is never mapped and written to as a whole. Every `Large` test prints one line:

`LARGE-REPORT test=<name> outcome=<status> seconds=<s> rss_anon_peak_mb=<MB> disk_mb=<MB that the files occupy> nominal_gb=<GB>`

(`rss_anon_peak_mb` is the process total, sampled every 20 ms; the host itself is about 100 MB.)
