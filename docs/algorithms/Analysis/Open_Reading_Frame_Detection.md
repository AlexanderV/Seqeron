# Open Reading Frame (ORF) Detection — GenomicAnalyzer

| Field | Value |
|-------|-------|
| Algorithm Group | Analysis |
| Test Unit ID | GENOMIC-ORF-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Core (`OrfScanner`) |
| Implementation Status | Production |
| Last Reviewed | 2026-11-19 (B09 review: F21–F26; A1-18/A1-19: NCBI ORFfinder 0.4.3 measured, F57; audit round 1 line-by-line pass WP26) |

## 1. Overview

`GenomicAnalyzer.FindOpenReadingFrames` enumerates potential Open Reading Frames (ORFs) in a DNA
sequence across all six reading frames (three forward, three on the reverse complement). An ORF is
a span that begins at a start codon (ATG) and ends at the first in-frame stop codon (TAA, TAG, TGA)
with no internal in-frame stop [1][2]. Following the canonical definition used by Rosalind and NCBI
ORFfinder ("ATG only"), every in-frame ATG terminated by a downstream in-frame stop is reported,
including nested ORFs that share a stop codon [1][3][5]. The method is exact (deterministic, no
heuristics) and runs in linear time (O(n + output)). The genetic code (stop codons, alternative start
codons), the start-codon policy, nested vs. longest-ORF-per-stop reporting, open-ended (partial) ORFs and
the minimum-length convention are options (`OrfSearchOptions`) that reproduce NCBI ORFfinder's engine
`COrf::FindOrfs` [5] and EMBOSS getorf `-find 3` [6]. The options of the NCBI standalone program that `COrf` does not define —
cross-frame "ignore nested ORFs" (`-n`), search region (`-b/-e`), the `-ml` clamp, `-strand`, the 5'-partial flag and the dual-coding
tables — were **measured on the real `ORFfinder` 0.4.3 binary** (A1-18/A1-19, 1750 runs, see §5.5) and are implemented
as options; `OrfSearchOptions.OrfFinderCompatible` bundles the convention and the library equals the binary on all 54 820 oracle rows.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

During translation the ribosome reads mRNA in non-overlapping triplets (codons). Translation
initiates at a start codon and terminates at a stop codon. Because a double-stranded DNA molecule
can be read in three frames on each strand, there are six possible frame translations [2]. An ORF
is a candidate protein-coding region delimited by a start and an in-frame stop codon [1][2].

### 2.2 Core Model

For the standard genetic code (NCBI transl_table=1) the start codon is ATG and the stop codons are
TAA, TAG, TGA [4]. Given a strand `S`, for each frame offset `f ∈ {0,1,2}` the codons are
`S[f..f+3), S[f+3..f+6), …`. An ORF is a maximal-prefix span `S[a..b+3)` where `S[a..a+3) = ATG`,
`S[b..b+3)` is a stop codon, `b ≥ a`, `(b−a) mod 3 = 0`, and no codon strictly between `a` and `b`
is a stop [1][2]. The translated protein candidate is the translation of `S[a..b)` (start through
the codon before the stop) [1]. The reverse complement is searched identically [1][2]. The stop set is
the canonical `GeneticCode` stop set of the chosen table (never a private copy); tables 27/28/31, whose
stop codons also code an amino acid, are rejected as in `Translator.FindOrfs`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported ORF `Sequence` begins with `ATG` | scan only opens an ORF at ATG [1][4] |
| INV-02 | Every reported ORF `Sequence` ends with TAA/TAG/TGA | span terminates at the first in-frame stop [1][4] |
| INV-03 | `Length` is divisible by 3 | span runs whole codons from start to stop inclusive [2] |
| INV-04 | `Length ≥ minLength` | length filter applied before yielding [3] |
| INV-05 | `Frame ∈ {1,2,3}`; `IsReverseComplement` selects the strand | six-frame search [1][2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA to scan | non-null; `DnaSequence`: case normalized, A/C/G/T only; `string` overloads (A1-21): any IUPAC nucleotide string (DNA or RNA, `U` read as `T`), case-insensitive |
| `minLength` | `int` | 100 | minimum ORF length in nucleotides, **stop codon included** (inclusive lower bound) | any int (≤ 0 admits every ORF) |
| `options` | `OrfSearchOptions` | `Default` | see below (overload `FindOpenReadingFrames(sequence, minLength, options)`) | non-null |

`OrfSearchOptions` (additive record; defaults = the two-argument method):

| Option | Default | Meaning | Reference |
|--------|---------|---------|-----------|
| `GeneticCode` | `GeneticCode.Standard` | stop codons (and alternative starts) from NCBI table 1–33 except 27/28/31 | NCBI ORFfinder `-g`; [4][5] |
| `StartCodons` | `AtgOnly` | `AtgOnly`, `GeneticCodeStarts` (ATG + the table's alternative initiation codons), `AnySenseCodon` (stop-to-stop) | ORFfinder `-s`; [3][5] |
| `NestedOrfs` | `true` | every in-frame start (Rosalind; `COrf longest_orfs = false`); `false`: only the most upstream start per stop (ORFfinder default, getorf) | [1][5][6] |
| `IncludeOpenEnded` | `false` | also report ORFs that run to the sequence end without a stop (`IsTerminated = false`, end = last complete codon) | `COrf` partial ORFs; [5][6] |
| `MinLengthExcludesStop` | `false` | compare `minLength` with the length without the stop codon | `COrf` "not including the stop"; getorf `-minsize`; [5][6] |
| `SearchBothStrands` | `true` | scan the reverse complement; a view of `Strands` (`true` = `Both`, `false` = `Plus`) | ORFfinder `-strand` |
| `Strands` | `Both` | `Both`, `Plus` or `Minus` | ORFfinder `-strand both\|plus\|minus` (measured) |
| `IgnoreNestedOrfs` | `false` | drop an ORF whose forward span lies completely inside the span of another reported ORF of any frame/strand (non-strict; identical spans: the plus-strand ORF is dropped) | ORFfinder `-n true` (measured, §5.5) |
| `RegionStart` / `RegionEnd` | `1` / `0` | 1-based inclusive search region (`-b` / `-e`; 0 = to the end); the region is cut out and scanned alone | ORFfinder `-b/-e` (measured, §5.5) |
| `OrfFinderMinimumLength` | `false` | `minLength` below 30 becomes 30, then rounded down to whole codons | ORFfinder `-ml` (measured, §5.5) |
| `DualCodingStopsTerminate` | `false` | tables 27/28/31: the codons the table lists as stops end ORFs instead of the code being rejected | ORFfinder `-g 27/28/31` (measured) |
| `OrfSearchOptions.OrfFinderCompatible` | (profile) | `GeneticCodeStarts`, `NestedOrfs = false`, `IncludeOpenEnded`, `MinLengthExcludesStop`, `OrfFinderMinimumLength`, `DualCodingStopsTerminate`; add `GeneticCode`, `Strands`, `IgnoreNestedOrfs`, `RegionStart/End` with `with { }` | ORFfinder defaults `-s 1 -ml 75 -n false` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `OrfInfo.Sequence` | `string` | ORF span, start codon through stop codon inclusive (through the last complete codon for an open-ended ORF) |
| `OrfInfo.Position` | `int` | 0-based start offset within the **scanned strand**: the input for forward ORFs, its reverse complement for `IsReverseComplement` ORFs |
| `OrfInfo.Frame` | `int` | reading-frame number 1–3 within the scanned strand (`Position % 3 + 1`) |
| `OrfInfo.IsReverseComplement` | `bool` | true if the ORF was found on the reverse complement |
| `OrfInfo.IsTerminated` | `bool` | false only for an open-ended ORF (option `IncludeOpenEnded`) = ORFfinder 3'-partial |
| `OrfInfo.IsStartPartial` | `bool` | the ORF begins within the first 3 nt of the scanned strand / region (`Position < 3`) = ORFfinder 5'-partial (also true for an ORF that begins with ATG there) |
| `OrfInfo.RegionStart` / `RegionLength` | `int` | offset / length of the searched region (`Position` is relative to it; defaults 0 / whole input) |
| `OrfInfo.Length` / `CodonCount` | `int` | nucleotide length / codon count (includes stop) |
| `OrfInfo.GetForwardStrandRange(L)` / `GetForwardStrandRange()` | `(int, int)` | 0-based inclusive input-strand span, stop included: forward `(RegionStart + Position, RegionStart + Position + Length − 1)`; reverse `(RegionStart + R − Position − Length, RegionStart + R − 1 − Position)` with R = region length (the input length without a region); the parameterless overload needs no length |

Coordinate convention (reverse strand): `Position` is a coordinate of the reverse complement, as in
`Translator.FindOrfs`; use `GetForwardStrandRange` (MCP: `start`/`end`) for input coordinates. The reverse
complement of `input[Start..End]` is `Sequence` (checked on 200 random sequences and by the oracle).

### 3.3 Preconditions and Validation

Null `sequence` throws `ArgumentNullException`. Sequences shorter than one codon, or with no ATG, or
with an ATG that has no downstream in-frame stop (unless `IncludeOpenEnded`), yield no ORFs. A genetic code
with dual-coding stop codons (27, 28, 31) throws `ArgumentException` eagerly; a null `options` throws
`ArgumentNullException`. Indexing is 0-based; the `Sequence` span is inclusive of the stop codon. The DNA
alphabet/case is normalized by `DnaSequence` (A/C/G/T only, N and other IUPAC symbols are rejected there). The
string overloads `FindOpenReadingFrames(string, minLength[, options])` (B09 A1-21; MCP accepts IUPAC input) take any IUPAC
nucleotide string (DNA or RNA, case-insensitive, `U` read as `T`; anything else → `ArgumentException` naming the symbol and
position) through the same `OrfScanner` core: an ambiguous codon is a stop only when every concrete codon is a stop (`TAR`, `TRA`
in table 1), is never a start (also in `AnySenseCodon` mode) and is otherwise read through, so `ATGNNNTAA` is one ORF and `NNN` neither
opens nor ends one; the reverse strand is the canonical IUPAC reverse complement and `GetForwardStrandRange` takes the string length.

## 4. Algorithm

### 4.1 High-Level Steps

1. Build the 64-entry stop and start lookup tables once from the canonical `GeneticCode` (`IsStopCodon` /
   `IsStartCodon`), `ATG` only, or every sense codon.
2. For each frame offset 0,1,2 of a strand (`OrfScanner.FindRegions`): one pass over the codons. A start codon is
   appended to the list of pending starts (only the first one when `NestedOrfs = false`); a stop codon emits an ORF
   for every pending start (ascending, so lengths decrease) and clears the list. At the end of the frame the pending
   starts are emitted as open-ended ORFs when `IncludeOpenEnded`.
3. Apply the length filter (`Length ≥ minLength`, or `Length − 3 ≥ minLength` with `MinLengthExcludesStop`; with
   `OrfFinderMinimumLength` the threshold first becomes `⌊max(minLength, 30) / 3⌋ · 3`) and build
   `OrfInfo`; run steps 2–3 on the reverse complement with `IsReverseComplement = true`. A region (`RegionStart/End`) is cut out
   before step 2 and coordinates are mapped back afterwards; `Strands` decides which strands run.
4. `IgnoreNestedOrfs`: collect the surviving ORFs of both strands as forward spans, sort by (start ↑, end ↓, minus before plus) and sweep with the
   running maximum end (`OrfScanner.FindNestedSpans`, O(m log m)); an ORF with `end ≤ maximum` is nested and dropped.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindOpenReadingFrames` | O(n + output) | O(pending starts of one stop-free stretch) + output | Single pass per frame; pending starts are released at the next stop. Before the review each ATG rescanned to its stop (O(n²) on an ATG-rich stop-free stretch; 300 000 ATG + stop took minutes, now milliseconds). The output itself is quadratic for nested ORFs (each ORF carries its sequence). 5 Mbp random DNA, six frames: ≈ 1 s (Release). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int)` and `(DnaSequence, int, OrfSearchOptions)`: six-frame ORF enumeration; `FindOpenReadingFrames(string, int[, OrfSearchOptions])`: the same over an IUPAC nucleotide string (§3.3).
- `OrfScanner.FindRegions` ([OrfScanner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/OrfScanner.cs)): the shared single-strand scanning core
  (start predicate, nested/open-ended policies); `Translator.FindOrfs` (getorf `-find 1` model) is built on it too.
- MCP `find_open_reading_frames` delegates to the overload and adds the input-strand `start`/`end` columns.

### 5.2 Current Behavior

Forward frames are yielded first, then reverse-complement frames; within a frame, ORFs are yielded
in order of start position. Every in-frame ATG that reaches a stop is reported, so overlapping/nested
ORFs sharing a stop are all returned (matches the Rosalind worked example) [1].

**Relation to NCBI ORFfinder and EMBOSS getorf (verified, see Evidence):** the NCBI engine `COrf::FindOrfs`
[5] reports (a) each ORF from the most upstream start by default (`NestedOrfs = false`), (b) 3'-partial ORFs,
(c) a minimum length that does **not** count the stop codon. The default of this method keeps the Rosalind
convention (all starts, complete ORFs only) and the self-consistent length convention (the reported `Length`,
stop included, is compared with `minLength`); `OrfSearchOptions` switches each choice to the NCBI/getorf
behaviour. The earlier documentation claimed the default "matches NCBI ORFfinder's nucleotide length filter":
it does not (differs by the 3-nt stop); the claim is corrected.

**Suffix tree decision (search-reuse evaluation):** Not used. ORF detection is a fixed-stride
six-frame codon scan, not occurrence enumeration of a query pattern against a text. There is no
"find all positions of pattern P" sub-problem; the suffix tree (`Contains`/`FindAllOccurrences`)
offers no advantage and would not change the required output. A direct linear codon scan is the
correct structure here.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Six-frame search (forward + reverse complement) [1][2].
- Start codon ATG (default) or the genetic code's alternative initiation codons or any sense codon (ORFfinder start modes); stop codons from the NCBI genetic code table (default Standard TAA/TAG/TGA) [4][5].
- Every ATG opening an ORF terminated by an in-frame stop is reported, including nested ORFs sharing
  a stop (canonical Rosalind semantics; reproduces the Rosalind sample exactly) [1]; optional most-upstream-start-per-stop reporting (NCBI ORFfinder, getorf) [5][6].
- Open-ended (3'-partial) ORFs on request [5][6].
- Minimum-length filtering in nucleotides (inclusive), with or without the stop codon [3][5][6]; the `-ml` clamp / whole-codon rule of the standalone program [7].
- NCBI ORFfinder 0.4.3 `-n` (cross-frame nested filter), `-b/-e` (region), `-strand`, 5'-partial flag and dual-coding tables 27/28/31, as measured on the binary [7] (§5.5).

**Intentionally simplified:**

- (none)

**Not implemented (scope, not defects):**

- ORFfinder's circular option `-c` (wrap-around ORFs): the real binary ignores it (§5.5, measured); there is nothing to implement.
- (Implemented, kept here as a pointer:) ambiguity codes in the input are handled by the string overloads (§3.3 contract paragraph; oracle `orf_ambiguous_oracle.tsv`, 2844 rows: Biopython definition for all, EMBOSS getorf 6.6.0 for the 184 A/C/G/T/N rows it resolves exactly — getorf mis-resolves OR-ed IUPAC codons such as `TRA`, a documented difference of `Translator.FindOrfs(string)`).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | ORF `Sequence` includes the terminating stop codon | Assumption | reported span/length include the stop (Length % 3 == 0); protein candidate excludes it | accepted | Wikipedia "length divisible by three, bounded by stop codons" [2]; INV-02/03; NCBI ORFfinder also reports the stop inside the ORF range [5] |
| 2 | Default `minLength` counts the stop codon | Deviation from NCBI (documented) | an ORF whose coding part is `minLength − 3` is kept; NCBI `-ml` and getorf `-minsize` would drop it | accepted, option `MinLengthExcludesStop` | NCBI `COrf` "length in base pairs not including the stop" [5]; keeps INV-04 `Length ≥ minLength` self-consistent for existing callers |
| 3 | Pre-existing greedy scan corrected | Deviation (fixed) | old code missed nested ORFs sharing a stop | fixed | now matches Rosalind sample [1] |
| 4 | Hard-coded TAA/TAG/TGA | Defect (fixed, F21) | wrong for tables 2–33 | fixed | stop set from `GeneticCode` |
| 5 | Per-ATG rescan | Defect (fixed, F22) | O(n²) on ATG-rich stop-free stretches | fixed | single pass |

### 5.5 NCBI ORFfinder 0.4.3, measured on the binary (A1-18 / A1-19)

`COrf` (§5 source [5]) defines only the per-frame ORF search; the standalone program's cross-frame `-n` filter, the `-b/-e` region, the `-c`
option, the `-ml` clamp and the output flags are defined nowhere in the open-source library. They were decided by running the **unmodified NCBI
binary** (`https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz`, sha256 `f30a0806…09f2`, `ORFfinder -version` = `0.4.3`,
Build-Date May 4 2021, Debian 12 x86_64) in Docker over 629 deterministic inputs (designed nesting / boundary / partial / start-stop / length
cases, 170 random sequences, 122 origin-spanning genomes, 60 region sequences, 33 edge files) — 1750 invocations, 1738 rc 0 — and comparing the output
(`-outfmt 0` ORF FASTA, 1 CDS FASTA, 2 text ASN.1, 3 feature table: the converter asserts that the four formats agree) with this implementation.
Package: `docs/Validation/review-2026-09/handoff/A1-18-orffinder/` (README, `returned/A1-18-orffinder-out.zip`); converter
`scripts/analysis/generate_orffinder_oracle.py` → `tests/…/TestData/Analysis/orffinder_cli_oracle.tsv` (54 820 results, 232 162 ORFs, 260 option sets).
`GenomicAnalyzer_FindOpenReadingFrames_OrfFinderCli_Tests` asserts **equality on every row**: strand, forward coordinates, order, 5'/3' partial flags and
the protein sequences. Measured rules:

| Item | Rule (all rows of the oracle) |
|------|-------------------------------|
| Coordinates | ORFfinder prints **0-based inclusive** `from:to` including the stop (`-outfmt 1` headers are 1-based); `from > to` = minus strand. ORFs are ordered plus strand first, each strand by frame of the strand-coordinate start then start (= the library order). |
| `-n true` | an ORF is dropped when its forward span lies completely inside another ORF of the result (any frame, either strand): 67 373 ORFs dropped in 15 306 (-n false, -n true) row pairs — 58 732 strictly inside, 3 373 sharing only the start, 3 539 sharing only the end, **1 729 with an identical span** (possible only across strands, mostly `-s 2`): the **plus-strand ORF is dropped, the minus-strand ORF kept**. 23 231 dropped ORFs are contained only by ORFs of the other strand. The containers are the ORFs left by `-ml`/`-strand` (with `-strand plus` a minus ORF cannot nest anything); containment is transitive. |
| `-ml` | values below 30 become 30 (USAGE.txt), then the length is used **in whole codons**: ORFs with ≥ `3·⌊ml/3⌋` coding nt are reported (`-ml 31` → 30, `-ml 74` → 72, `-ml 76` → 75); the stop is not counted for a terminated ORF, an open ORF counts all its nt. `-ml −5 … 29` = 30. |
| `-b/-e` | 1-based, inclusive; `-b 0` = 1; `-e 0` or beyond the sequence = to the end; `-b < 0` ("Requested interval start -5 is out of range") and `-b > -e` ("Interval Start value is expected to be smaller than Stop") are errors for the sequence (360 error rows). The region is **cut out and scanned as a sequence of its own** (24 pairs × strands × `-n` × `-s`): ORFs crossing a region end are open (partial) there, the reverse frames are anchored at the region end, coordinates are in input space. |
| `-c true` | **no effect**: output, stderr and exit code are byte-identical to `-c false` on 36 pairs of 122 origin-spanning genomes, 20 region pairs and 12 no-stop genomes (no wrap-around ORF is ever reported; USAGE.txt: "Under development"). Nothing to implement; `OrfSearchOptions` has no circular option. |
| partial | `partial` in the `-outfmt 0` header and `fuzz-from/to`, `no-left/no-right` in `-outfmt 2`: **3'-partial** = no stop codon (ORF ends at the last complete codon); **5'-partial** = the ORF starts at strand offset 0, 1 or 2 (22 000 ORFs; also when it begins with ATG). Feature tables (`-outfmt 3`) carry no `<`/`>`. |
| Protein | translation with the table, stop dropped; a 5'-complete ORF begins with `M` whatever its first codon (alternative start, any sense codon in `-s 2`), a 5'-partial ORF keeps its first codon's amino acid. |
| `-s` / `-g` | `-s 0/1/2` and every table 1–6, 9–16, 21–33 equal the shipped scanner (`-s 1` = ATG + the table's start codons; 30 000 ORF rows); tables **27, 28, 31 are accepted** and the codons the table lists as stops end ORFs; `-g 0, 7, 8, 17–20` fail with "Unable to find genetic code number N"; `-s` outside 0–2 and `-outfmt` outside 0–3 are errors (rc 254 / 255). |
| `-strand` | `both`/`plus`/`minus` (case-insensitive); other values: warning, both strands. Boolean options accept `true/false/t/f/yes/no/TRUE/1`; `maybe` is a usage error (rc 1). |
| Input | upper/lower/mixed case accepted; `N` and IUPAC codons follow the `OrfScanner` rules (a codon with an ambiguity code is a stop only if every expansion is a stop — `TAR`; never a start); `U` (RNA), gaps and protein are rejected by the FASTA reader ("invalid residue(s)"), as are some very short or C-rich records (27 rows, no result recorded). The string overloads accept `U` (read as `T`) — a deliberate superset. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null sequence | `ArgumentNullException` | contract |
| sequence shorter than a codon | empty | no codon to scan |
| ATG with no in-frame stop | not reported | incomplete ORF [1] |
| nested ATGs sharing one stop | all reported | canonical ORF definition [1] |
| ORF only on reverse strand | reported with `IsReverseComplement=true`; `Position` in reverse-complement coordinates, `GetForwardStrandRange` gives input coordinates | six-frame search [1][2] |
| ATG at the last complete codon, no stop | not reported (with `IncludeOpenEnded`: one-codon open ORF, NCBI `COrf` reports it; getorf does not) | [1][5][6] |
| genetic code 27 / 28 / 31 | `ArgumentException`; with `DualCodingStopsTerminate` the table's stop codons end ORFs | dual-coding stops; ORFfinder accepts them [7] |
| `RegionStart` / `RegionEnd` | 0-based `Position` relative to the region; `RegionStart < 0`, `RegionEnd < 0` or `RegionStart > RegionEnd` → `ArgumentOutOfRangeException` (ORFfinder reports the error per sequence); `RegionStart > length` with `RegionEnd = 0` is the same error (not measured) | [7] |
| ORF begins at strand offset 0–2 | `IsStartPartial` = true (also with ATG) | [7] |

### 6.2 Limitations

Does not model codon bias, ribosome binding sites, or splicing; ORF presence is not evidence of a real
gene [2]. Input is a `DnaSequence` (A/C/G/T) or an IUPAC nucleotide string. `-c` (circular genomes) is not modelled because the reference program ignores it:
an ORF crossing the origin of a circular genome is not found (not an ORFfinder behaviour, so also not implemented as an extension here).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var dna = new DnaSequence("ATGAAAAAATAA");
var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 1).ToList();
// orfs[0].Sequence == "ATGAAAAAATAA"; Position == 0; Frame == 1; protein candidate == "MKK"
```

**Numerical / biological walk-through:**

Rosalind sample `AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG`,
searched in all six frames, yields the four distinct protein candidates
`MLLGSFRLIPKETLIQVAGSSPCNLS`, `M`, `MGMTPRLGLESLLE`, `MTPRLGLESLLE` (the last two share a stop) [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomicAnalyzer_FindOpenReadingFrames_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_FindOpenReadingFrames_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [GENOMIC-ORF-001-Evidence.md](../../../docs/Evidence/GENOMIC-ORF-001-Evidence.md)
- ORFfinder binary oracle: [GenomicAnalyzer_FindOpenReadingFrames_OrfFinderCli_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_FindOpenReadingFrames_OrfFinderCli_Tests.cs) (54 820 rows, 232 162 ORFs, all equal), generator `scripts/analysis/generate_orffinder_oracle.py`
- Oracle tests: [GenomicAnalyzer_FindOpenReadingFrames_Oracle_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_FindOpenReadingFrames_Oracle_Tests.cs) (1746 rows: Biopython finder = NCBI `COrf` port = EMBOSS getorf; differential tests against `Translator.FindOrfs` and `GenomeAnnotator.FindOrfs`), generator `scripts/analysis/generate_orf_oracle.py`
- Related ORF finders (same model, different policies — see the Evidence "three ORF engines" table): ANNOT-ORF-001 (`GenomeAnnotator.FindOrfs`), TRANS-SIXFRAME-001 (`Translator.FindOrfs`)

## 8. References

1. Rosalind. 2026. ORF — Open Reading Frames. https://rosalind.info/problems/orf/
2. Wikipedia contributors. 2026. Open reading frame. https://en.wikipedia.org/wiki/Open_reading_frame
3. NCBI. 2026. ORFfinder (Open Reading Frame Finder). https://www.ncbi.nlm.nih.gov/orffinder/
4. NCBI. 2026. The Genetic Codes (transl_table=1, Standard). https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
5. NCBI C++ Toolkit, `src/algo/sequence/orf.cpp` and `include/algo/sequence/orf.hpp` (`COrf::FindOrfs`, `FindForwardOrfs`, `s_FindOrfs`) — the ORF engine of NCBI ORFfinder. https://github.com/ncbi/ncbi-cxx-toolkit-public (opened 2026-11-18)
6. EMBOSS 6.6.0 `getorf` (`-find 3`, `-minsize`, `-table`). http://emboss.open-bio.org/rel/rel6/apps/getorf.html
7. NCBI ORFfinder 0.4.3 standalone program (`ORFfinder.gz`, `USAGE.txt`, `CHANGELOG.txt`), https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/ — run in Docker (linux/amd64, Debian 12) over the hand-off package `docs/Validation/review-2026-09/handoff/A1-18-orffinder/`; output archive `returned/A1-18-orffinder-out.zip` (1750 runs, 2026-11-18).

### Lazy sequence (B09 WP78, F361)

`OrfInfo.Sequence` is cut out of the scanned strand on first access (`SequenceMemory` is the view without a copy; `Length`, `CodonCount`, `Position` need no string). The nested-ORF input ATG x 8000 + TAA with `minLength` 3 (8 000 ORFs, 96 036 000 characters) is enumerated in 32 ms / 0.7 MB instead of 120 ms / 183.7 MB; `RepeatInfo.Sequence` (`FindRepeats`) and `CommonRegion.Sequence` (`FindCommonRegions`) are lazy in the same way. Values and equality are those of the eager records.
