# A1-18 / A1-19 - integration plan for the returned NCBI ORFfinder 0.4.3 output

> **DONE (B09 WP10, F57)** - integrated as described in README.md; deviations from this plan: the converter is `scripts/analysis/generate_orffinder_oracle.py` (compact result-table TSV `orffinder_cli_oracle.tsv`), `OrfSearchOptions` members are `IgnoreNestedOrfs`, `OrfFinderMinimumLength`, `DualCodingStopsTerminate`, `Strands`, `RegionStart/RegionEnd`, profile `OrfFinderCompatible` (no `Circular`: the binary ignores `-c`).

> **Archive retention (B09 audit A2-45 c, 2026-11-19):** `returned/A1-18-orffinder-out.zip` (17.7 MB, sha256 `0a5cf9ce8759c44f1a6f9d23b04713e05500c04cadc6c016ca112fcd08fab09b`) is KEPT in the tree as provenance: `scripts/analysis/generate_orffinder_oracle.py` reads it by default and it is the only record of the raw output of the real NCBI `ORFfinder` 0.4.3 (not runnable in the review sandbox); the committed `orffinder_cli_oracle.tsv` is derived from it. It is reproducible from this package (`./run.sh` on a machine with Docker, `--platform linux/amd64` on Apple Silicon). The package is not "TEMPORARY" any more; do not delete the archive without first changing the generator to read a different source.

Written before the output exists. Every "decided by the oracle" item below is a question the 1750 runs answer; the plan fixes the
code shape, so the integrating session only fills in the measured rule. No rule is implemented before it is observed in `out/`.

## 0. Published facts already in hand (USAGE.txt, sha256 368d0fe3...b409, read in full)

`ORFfinder [-in file | -id acc] [-b begin] [-e end] [-c circular] [-g 1..31] [-s 0|1|2] [-ml N] [-n bool] [-strand both|plus|minus]
[-out file] [-outfmt 0|1|2|3]`; defaults `-b 1 -e 0 (to the end) -c false -g 1 -s 1 -ml 75 -n false -strand both -outfmt 0`;
`-s` 0 = ATG only, 1 = ATG + alternative initiation codons, 2 = any sense codon; `-ml` "Value less than 30 is automatically changed
by 30"; `-n` "Ignore nested ORFs (completely placed within another)"; `-c` "Is the sequence circular? *** Under development";
`-outfmt` 0 ORF list FASTA, 1 CDS FASTA, 2 text ASN.1, 3 feature table. CHANGELOG: v0.4.0 added `-s`, v0.3.0 `-strand`, v0.2.x multi-FASTA.

## 1. Step 0 - integrate and sanity-check the output (no library change)

1. Unzip `returned/A1-18-orffinder-out.zip`; read `out/ENV.txt` (`ORFfinder -version` must say 0.4.3; sha256 list), `out/rc_histogram.txt`,
   `out/runs.tsv` (timeouts rc 124 on the `circedge` runs are an answer for `-c`).
2. `python3 parse_orffinder_out.py --probe out/` shows real `-outfmt 0..3` text; fix `HEADER_PATTERNS` (the feature-table parser
   is already per INSDC spec). The converter asserts that fmt 0, fmt 1 and fmt 3 agree for every sequence before it emits a row.
3. Promote the converter to `scripts/analysis/generate_orffinder_cli_oracle.py` (reads this package's `inputs/` + the returned
   `out/`, writes `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/orffinder_cli_oracle.tsv`, an EmbeddedResource like
   `orf_oracle.tsv`). Row: `config  seqId  sequence  orfs` with `config = g,s,ml,n,strand,c,b,e` and orfs =
   `strand,start,end,partial5,partial3` (1-based inclusive forward coordinates, ORF 5'->3' partial flags). It additionally
   **asserts the binary against the three existing engines** of `generate_orf_oracle.py` (Biopython definition finder, `COrf` port,
   getorf) on every row whose definition coincides (`-n false -c false -b 1 -e 0`, any `-g/-s/-strand`, `-ml >= 30`):
   that equality is the acceptance test for the already-shipped options (and shows how ORFfinder treats 5'-partial ORFs, the stop in `-ml`, `-g 27/28/31`).
4. Everything the binary does that the library does not (rows differing from the shipped finder) is listed in the F-entry; it is the work list of step 2.

## 2. Planned C# changes (all additive, defaults unchanged)

Files: `src/Seqeron/Algorithms/Seqeron.Genomics.Core/OrfScanner.cs`,
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs` (`FindOpenReadingFramesCore`, `OrfInfo`, `OrfSearchOptions`),
MCP `src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs` (`find_open_reading_frames`, ~line 745).

| ORFfinder option | Today | Planned member | Decided by the oracle |
|---|---|---|---|
| `-n true` (A1-18) | not available | `OrfSearchOptions.IgnoreNestedOrfs` (bool, default false). Implemented as a post-filter on the *forward-coordinate spans* of the ORFs of both strands after the per-strand scan: sort by (start asc, end desc), sweep with the running maximum end, O(m log m). Lives in a new `OrfScanner`-level helper `RemoveNested(IReadOnlyList<(Strand, Start, End, ...)>)` so `Translator.FindOrfs`-style callers could reuse it. | strict vs non-strict containment (equal spans across strands only occur with `-s 2`: `strand_equal_span_s2_*`, `strand_share_boundary_s2_*`); whether the stop codon counts in the span (the `xnest_edge_*`, `xoverlap_*` designs have 1-4 nt boundary offsets); same-frame nesting; whether the filter sees the other strand when `-strand plus/minus` hides it; interaction with `-ml` (geometrically order-independent: a container shorter than `-ml` contains only shorter ORFs) and with 3'/5'-partial ORFs; interaction with `-s 1` alternative starts (`xnest_altstart_*`); whether an ORF removed as nested still suppresses others (transitivity). |
| `-ml` clamp | any value | `OrfSearchOptions.MinimumLengthFloor` (int, default 0) = 30 for ORFfinder; effective `minLength = max(minLength, floor)`. Convenience preset `OrfSearchOptions.NcbiOrfFinder` = {GeneticCodeStarts, nested false (longest per stop), IncludeOpenEnded as observed, MinLengthExcludesStop true, MinimumLengthFloor 30, both strands} with C# `minLength` default 75. | the `mllen_*` rows (6..90 nt with `-ml` 0,1,10,29,30,31,45,74,75,76,90): inclusive bound, stop counted or not (the existing `MinLengthExcludesStop` claim is re-checked), negative `-ml`, `-ml 29` = 30. |
| `-b / -e` (A1-19) | `GenomicAnalyzer` scans the whole input | `OrfSearchOptions.RegionStart` (1-based, default 1) and `RegionEnd` (1-based inclusive, 0 = to the end), exactly the CLI convention, validated as the CLI does (rc/stderr of `b0`, `b>e`, `e>len`, `-5`). `FindOpenReadingFramesCore` slices the region, scans it, maps coordinates back; `OrfInfo` gets additive `RegionStart`/`RegionLength` (defaults 0 / -1 = whole input) so the existing `GetForwardStrandRange(L)` keeps its meaning and a new parameterless overload returns the input-coordinate span. | cut-then-scan (ORFs open at the region ends, reverse frames anchored at the region end) vs scan-then-clip; coordinates in input or region space; what a region end does to a partial ORF; behaviour for the invalid pairs. |
| `-c true` (A1-19) | not available | **if the binary wraps ORFs around the origin:** `OrfSearchOptions.Circular` + `OrfScanner.FindCircularRegions` (scan the strand unrolled three times - 3L covers every frame class when L is not a multiple of 3 - report ORFs whose start lies in the first copy, dedupe by (start mod L, length), limit length < 3L; forward coordinates wrap: `OrfInfo.WrapsOrigin`, `GetForwardStrandRange` returns Start > End for a wrapped ORF or a new `WrappedRange`). **If `-c` has no observable effect** (output identical to `-c false` on `circular.fa`; USAGE says "under development"): record exactly that as the concrete proof in B09/GENOMIC-ORF-001/the algorithm doc, implement nothing, tick A1-19 for `-b/-e` only. `circedge` timeouts/exit codes decide the no-stop genomes. | wrap-around rows (`circ_*`: cut at start 1/2/3, mid, stop 1/2, before/after), nested/cross-frame ORFs across the origin, genomes of length not divisible by 3, partial flags, `-n` interplay. |
| `-strand plus/minus` | `SearchBothStrands` (both/plus) | `OrfSearchOptions.Strands` (`OrfStrands.Both/Plus/Minus`, default Both); `SearchBothStrands` kept as a compatible view (`get => Strands == Both`, `init => Strands = value ? Both : Plus`). | factorial rows (`-strand` x `-n` x `-s`). |
| `-s 0/1/2`, `-g` | `OrfSearchOptions.StartCodons`, `GeneticCode` | none expected; any deviation is a defect fixed in `OrfScanner`/adapter and noted as F-entry. `-g 27/28/31` (rejected today as dual-coding) and `-g 0/32/33`: take the binary's behaviour into the docs (+ MCP validation messages). | `gcode__*`, factorial rows vs the shipped finder. |
| partial ORFs | `IncludeOpenEnded` = 3'-partial only | if ORFfinder also reports 5'-partial ORFs (stop-to-stop from the sequence start without a start codon): `OrfRegion.IsStartPartial` + scanner option `includeOpenStart`, `OrfInfo.IsStartPartial`; else documented as measured. | `partial*`, `len_*`, region and `-c` rows; feature-table `<` / `>` markers. |
| N / IUPAC / lower case / U | n/a here | observation only, feeds A1-21 (string overloads): which codons containing N the binary treats as stop/start/neither (`e10..e15`, `e23`). | `edge__*`. |

MCP `find_open_reading_frames`: additive optional parameters `ignoreNested` (bool), `regionStart`, `regionEnd` (1-based, as
`-b/-e`; the existing `start`/`end` result fields stay), `circular` (only if implemented), `strand` ("both"|"plus"|"minus"),
`ncbiMinLength` (the floor of 30). Result items gain `isStartPartial` / `wrapsOrigin` only when those features exist. Then
`python3 scripts/skills/gen-catalog.py` and `--check`, link check as usual.

## 3. Tests that must pass

1. `GenomicAnalyzer_OrfTests.FindOpenReadingFrames_MatchesNcbiOrfFinderCliOracle` over **every row** of `orffinder_cli_oracle.tsv`
   (all configs of the factorial / ml / gcode / region / circular groups), options mapped exactly as in section 2; coordinates,
   strands, partial flags and order-insensitive set equality (ORFfinder's own ordering is recorded but not a contract).
2. One focused unit test per option with a worked example taken from the oracle (nested: an `xnest_*` sequence with and without
   `-n`; region: a `region.fa` row; clamp: `mllen_*` at `-ml 0/29/30/31`; circular: a wrap row or the "no effect" proof; strand).
3. Unchanged: the 1746-row `orf_oracle.tsv` test, `GENOMIC-ORF-001` unit tests, `Translator.FindOrfs` getorf tests,
   `Differential_*` tests. MCP: `FindOpenReadingFramesTests` + new tests for the added parameters.
4. Targeted tier: `FullyQualifiedName~GenomicAnalyzer|~Translator|~OrfScanner|~GenomeAnnotator` + `Seqeron.Mcp.Analysis.Tests`, then the
   lead's fast tier. Docs: `Open_Reading_Frame_Detection.md` (delete the BLOCKED/"not part of the method" text of §5.3/§6.2),
   `GENOMIC-ORF-001` TestSpec + Evidence (record the binary version, sha256, the oracle row counts), B09 F-entry + worklist ticks, `LIMITATIONS` text in B09.
5. Oracle provenance in the Evidence: `ORFfinder -version` output from `out/ENV.txt`, sha256 of ORFfinder.gz, run count, the package path.

## 4. If the binary cannot be run on the user's side

Fall back (only then) to the secondary sources listed in A1-18: the ProteinIQ ORFfinder wrapper description and orfipy; neither defines
tie rules, so a rule not measured on the binary stays unimplemented and A1-18 stays open with that proof.
