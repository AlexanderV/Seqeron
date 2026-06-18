# Open Reading Frame (ORF) Detection — GenomicAnalyzer

| Field | Value |
|-------|-------|
| Algorithm Group | Analysis |
| Test Unit ID | GENOMIC-ORF-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

`GenomicAnalyzer.FindOpenReadingFrames` enumerates potential Open Reading Frames (ORFs) in a DNA
sequence across all six reading frames (three forward, three on the reverse complement). An ORF is
a span that begins at a start codon (ATG) and ends at the first in-frame stop codon (TAA, TAG, TGA)
with no internal in-frame stop [1][2]. Following the canonical definition used by Rosalind and NCBI
ORFfinder ("ATG only"), every in-frame ATG terminated by a downstream in-frame stop is reported,
including nested ORFs that share a stop codon [1][3]. The method is exact (deterministic, no
heuristics) and runs in linear time.

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
the codon before the stop) [1]. The reverse complement is searched identically [1][2].

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
| `sequence` | `DnaSequence` | required | DNA to scan | non-null; case normalized by `DnaSequence` |
| `minLength` | `int` | 100 | minimum ORF length in nucleotides (inclusive lower bound) [3] | any int; ORFs shorter are excluded |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `OrfInfo.Sequence` | `string` | ORF span, start codon through stop codon inclusive |
| `OrfInfo.Position` | `int` | 0-based start offset within the scanned strand |
| `OrfInfo.Frame` | `int` | reading-frame number 1–3 within the scanned strand |
| `OrfInfo.IsReverseComplement` | `bool` | true if the ORF was found on the reverse complement |
| `OrfInfo.Length` / `CodonCount` | `int` | nucleotide length / codon count (includes stop) |

### 3.3 Preconditions and Validation

Null `sequence` throws `ArgumentNullException`. Sequences shorter than one codon, or with no ATG, or
with an ATG that has no downstream in-frame stop, yield no ORFs. Indexing is 0-based; the `Sequence`
span is inclusive of the stop codon. The DNA alphabet/case is normalized by `DnaSequence`; codons
containing characters that do not match ATG or a stop codon are simply not start/stop sites.

## 4. Algorithm

### 4.1 High-Level Steps

1. For each forward frame offset 0,1,2: scan codons; at each ATG, find the first in-frame stop
   downstream and report the ATG→stop span (if length ≥ minLength). Continue scanning from the next
   codon so every ATG is considered (nested ORFs are reported).
2. Compute the reverse complement and repeat for its three frames, marking results
   `IsReverseComplement = true`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindOpenReadingFrames` | O(n²) worst case | O(1) excl. output | Per ATG the inner loop scans to the next in-frame stop; in the pathological case of many ATGs before a stop this is quadratic, but is O(n) in the typical case (bounded inter-stop distance). n = sequence length. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int)`: six-frame ORF enumeration.
- `GenomicAnalyzer.FindOrfsInFrame(...)` (private): single-frame ATG→first-in-frame-stop scan over codons.

### 5.2 Current Behavior

Forward frames are yielded first, then reverse-complement frames; within a frame, ORFs are yielded
in order of start position. Every in-frame ATG that reaches a stop is reported, so overlapping/nested
ORFs sharing a stop are all returned (matches the Rosalind worked example) [1].

**Suffix tree decision (search-reuse evaluation):** Not used. ORF detection is a fixed-stride
six-frame codon scan, not occurrence enumeration of a query pattern against a text. There is no
"find all positions of pattern P" sub-problem; the suffix tree (`Contains`/`FindAllOccurrences`)
offers no advantage and would not change the required output. A direct linear codon scan is the
correct structure here.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Six-frame search (forward + reverse complement) [1][2].
- Start codon ATG; stop codons TAA/TAG/TGA (standard genetic code) [4].
- Every ATG opening an ORF terminated by an in-frame stop is reported, including nested ORFs sharing
  a stop (canonical Rosalind semantics; reproduces the Rosalind sample exactly) [1].
- Minimum-length filtering in nucleotides (inclusive) [3].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Alternative initiation codons (GTG/TTG) and non-standard genetic codes: feature limited to
  ATG/standard code (NCBI ORFfinder default "ATG only"); **users should rely on:**
  `GenomeAnnotator.FindOrfs` / `Translator.FindOrfs` for genetic-code-parameterized ORF finding [3].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | ORF `Sequence` includes the terminating stop codon | Assumption | reported span/length include the stop (Length % 3 == 0); protein candidate excludes it | accepted | Wikipedia "length divisible by three, bounded by stop codons" [2]; INV-02/03 |
| 2 | `minLength` is in nucleotides, inclusive | Assumption | filter boundary | accepted | NCBI ORFfinder nucleotide length options [3]; INV-04 |
| 3 | Pre-existing greedy scan corrected | Deviation (fixed) | old code missed nested ORFs sharing a stop | fixed | now matches Rosalind sample [1] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null sequence | `ArgumentNullException` | contract |
| sequence shorter than a codon | empty | no codon to scan |
| ATG with no in-frame stop | not reported | incomplete ORF [1] |
| nested ATGs sharing one stop | all reported | canonical ORF definition [1] |
| ORF only on reverse strand | reported with `IsReverseComplement=true` | six-frame search [1][2] |

### 6.2 Limitations

ATG-only / standard-code only (no alternative start codons or non-standard codes); does not model
codon bias, ribosome binding sites, or splicing. ORF presence is not evidence of a real gene [2].

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

- Tests: [GenomicAnalyzer_FindOpenReadingFrames_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindOpenReadingFrames_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [GENOMIC-ORF-001-Evidence.md](../../../docs/Evidence/GENOMIC-ORF-001-Evidence.md)
- Related algorithm (genetic-code-parameterized ORF finding): ANNOT-ORF-001 (`GenomeAnnotator.FindOrfs`)

## 8. References

1. Rosalind. 2026. ORF — Open Reading Frames. https://rosalind.info/problems/orf/
2. Wikipedia contributors. 2026. Open reading frame. https://en.wikipedia.org/wiki/Open_reading_frame
3. NCBI. 2026. ORFfinder (Open Reading Frame Finder). https://www.ncbi.nlm.nih.gov/orffinder/
4. NCBI. 2026. The Genetic Codes (transl_table=1, Standard). https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
