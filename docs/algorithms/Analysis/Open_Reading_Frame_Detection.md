# Open Reading Frame (ORF) Detection — GenomicAnalyzer

| Field | Value |
|-------|-------|
| Algorithm Group | Analysis |
| Test Unit ID | GENOMIC-ORF-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Core (`OrfScanner`) |
| Implementation Status | Production |
| Last Reviewed | 2026-11-18 (B09 review: F21–F26) |

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
`COrf::FindOrfs` [5] and EMBOSS getorf `-find 3` [6].

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
| `sequence` | `DnaSequence` | required | DNA to scan | non-null; case normalized by `DnaSequence`, A/C/G/T only |
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
| `SearchBothStrands` | `true` | scan the reverse complement | ORFfinder `-strand` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `OrfInfo.Sequence` | `string` | ORF span, start codon through stop codon inclusive (through the last complete codon for an open-ended ORF) |
| `OrfInfo.Position` | `int` | 0-based start offset within the **scanned strand**: the input for forward ORFs, its reverse complement for `IsReverseComplement` ORFs |
| `OrfInfo.Frame` | `int` | reading-frame number 1–3 within the scanned strand (`Position % 3 + 1`) |
| `OrfInfo.IsReverseComplement` | `bool` | true if the ORF was found on the reverse complement |
| `OrfInfo.IsTerminated` | `bool` | false only for an open-ended ORF (option `IncludeOpenEnded`) |
| `OrfInfo.Length` / `CodonCount` | `int` | nucleotide length / codon count (includes stop) |
| `OrfInfo.GetForwardStrandRange(L)` | `(int, int)` | 0-based inclusive input-strand span, stop included: forward `(Position, Position + Length − 1)`; reverse `(L − Position − Length, L − 1 − Position)` |

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
3. Apply the length filter (`Length ≥ minLength`, or `Length − 3 ≥ minLength` with `MinLengthExcludesStop`) and build
   `OrfInfo`; run steps 2–3 on the reverse complement with `IsReverseComplement = true`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindOpenReadingFrames` | O(n + output) | O(pending starts of one stop-free stretch) + output | Single pass per frame; pending starts are released at the next stop. Before the review each ATG rescanned to its stop (O(n²) on an ATG-rich stop-free stretch; 300 000 ATG + stop took minutes, now milliseconds). The output itself is quadratic for nested ORFs (each ORF carries its sequence). 5 Mbp random DNA, six frames: ≈ 1 s (Release). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int)` and `(DnaSequence, int, OrfSearchOptions)`: six-frame ORF enumeration.
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
- Minimum-length filtering in nucleotides (inclusive), with or without the stop codon [3][5][6].

**Intentionally simplified:**

- (none)

**Not implemented (scope, not defects):**

- ORFfinder's *cross-frame* "ignore nested ORFs" filter (drop an ORF lying completely inside an ORF of another frame/strand) and its search-region / circular options: they are post-filters on this method's output (coordinates are provided via `GetForwardStrandRange`); NCBI itself marks circular handling "under development".
- Ambiguity codes in the input are handled by the string overloads (see §6.1 contract paragraph; oracle `orf_ambiguous_oracle.tsv`, 2844 rows: Biopython definition for all, EMBOSS getorf 6.6.0 for the 184 A/C/G/T/N rows it resolves exactly — getorf mis-resolves OR-ed IUPAC codons such as `TRA`, a documented difference of `Translator.FindOrfs(string)`).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | ORF `Sequence` includes the terminating stop codon | Assumption | reported span/length include the stop (Length % 3 == 0); protein candidate excludes it | accepted | Wikipedia "length divisible by three, bounded by stop codons" [2]; INV-02/03; NCBI ORFfinder also reports the stop inside the ORF range [5] |
| 2 | Default `minLength` counts the stop codon | Deviation from NCBI (documented) | an ORF whose coding part is `minLength − 3` is kept; NCBI `-ml` and getorf `-minsize` would drop it | accepted, option `MinLengthExcludesStop` | NCBI `COrf` "length in base pairs not including the stop" [5]; keeps INV-04 `Length ≥ minLength` self-consistent for existing callers |
| 3 | Pre-existing greedy scan corrected | Deviation (fixed) | old code missed nested ORFs sharing a stop | fixed | now matches Rosalind sample [1] |
| 4 | Hard-coded TAA/TAG/TGA | Defect (fixed, F21) | wrong for tables 2–33 | fixed | stop set from `GeneticCode` |
| 5 | Per-ATG rescan | Defect (fixed, F22) | O(n²) on ATG-rich stop-free stretches | fixed | single pass |

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
| genetic code 27 / 28 / 31 | `ArgumentException` | dual-coding stops |

### 6.2 Limitations

Does not model codon bias, ribosome binding sites, or splicing; ORF presence is not evidence of a real
gene [2]. Input is a `DnaSequence` (A/C/G/T) or an IUPAC nucleotide string; the ORFfinder cross-frame nested filter and region/circular options are
not part of the method (see §5.3).

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
- Oracle tests: [GenomicAnalyzer_FindOpenReadingFrames_Oracle_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_FindOpenReadingFrames_Oracle_Tests.cs) (1746 rows: Biopython finder = NCBI `COrf` port = EMBOSS getorf; differential tests against `Translator.FindOrfs` and `GenomeAnnotator.FindOrfs`), generator `scripts/analysis/generate_orf_oracle.py`
- Related ORF finders (same model, different policies — see the Evidence "three ORF engines" table): ANNOT-ORF-001 (`GenomeAnnotator.FindOrfs`), TRANS-SIXFRAME-001 (`Translator.FindOrfs`)

## 8. References

1. Rosalind. 2026. ORF — Open Reading Frames. https://rosalind.info/problems/orf/
2. Wikipedia contributors. 2026. Open reading frame. https://en.wikipedia.org/wiki/Open_reading_frame
3. NCBI. 2026. ORFfinder (Open Reading Frame Finder). https://www.ncbi.nlm.nih.gov/orffinder/
4. NCBI. 2026. The Genetic Codes (transl_table=1, Standard). https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
5. NCBI C++ Toolkit, `src/algo/sequence/orf.cpp` and `include/algo/sequence/orf.hpp` (`COrf::FindOrfs`, `FindForwardOrfs`, `s_FindOrfs`) — the ORF engine of NCBI ORFfinder. https://github.com/ncbi/ncbi-cxx-toolkit-public (opened 2026-11-18)
6. EMBOSS 6.6.0 `getorf` (`-find 3`, `-minsize`, `-table`). http://emboss.open-bio.org/rel/rel6/apps/getorf.html
