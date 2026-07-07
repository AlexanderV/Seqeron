# Regulatory Element Scan

| Field | Value |
|-------|-------|
| Algorithm Group | Motif Discovery (Matching) |
| Test Unit ID | MOTIF-REGULATORY-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

`FindRegulatoryElements` scans a DNA sequence for a fixed library of well-characterised regulatory consensus motifs (eukaryotic and prokaryotic promoter elements, translation-initiation signals, a polyadenylation signal, and several transcription-factor binding sites). Each library entry is matched as an exact or IUPAC-degenerate consensus string and every occurrence is reported with its 0-based start position. The algorithm is specification-driven: detection is exact pattern matching against published consensus sequences, not a probabilistic score, so a hit means the input literally contains the cited consensus (or, for the E-box, a member of its IUPAC family) [1][2][3][4][5][6][7][8][9].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Regulatory elements are short, conserved DNA sequences recognised by the transcription/translation machinery. The "consensus" of an element is the most-frequent base at each aligned position across many functional examples; individual instances may vary, but the consensus is the canonical signature used for recognition [2]. The library mixes eukaryotic core-promoter elements (TATA box, CCAAT box, GC box) [1][3], prokaryotic σ70 promoter hexamers (-10 Pribnow box, -35 box) [2], a bacterial ribosome-binding site (Shine-Dalgarno), the eukaryotic translation-initiation Kozak context [4], the poly(A) signal [5], and four transcription-factor binding sites (E-box, AP-1, NF-κB, CREB) [6][7][8][9].

### 2.2 Core Model

For a sequence `S` of length `n` and a consensus pattern `P` of length `m`, the scan reports every start index `i` with `0 <= i <= n - m` such that for all `j`, `S[i+j]` is in the IUPAC base set of `P[j]`. Plain bases match themselves; `N` matches any of A/C/G/T; the remaining IUPAC ambiguity codes match their defined subsets. The library consensus strings are (5'→3'):

| Element | Pattern | Source |
|---------|---------|--------|
| TATA Box | `TATAAA` | Bucher (1990) [1] |
| CAAT Box | `CCAAT` | Bucher (1990) [1] |
| GC Box | `GGGCGG` | Lundin, Nehlin & Ronne (1994) [3] |
| -10 Box (Pribnow) | `TATAAT` | Harley & Reynolds (1987) [2] |
| -35 Box | `TTGACA` | Harley & Reynolds (1987) [2] |
| Kozak | `GCCGCCACCATGG` | Kozak (1987) [4] |
| Shine-Dalgarno | `AGGAGG` | [10] |
| Poly(A) Signal | `AATAAA` | Proudfoot & Brownlee (1976) [5] |
| E-box | `CANNTG` | Massari & Murre (2000) [6] |
| AP-1 | `TGACTCA` | Lee, Mitchell & Tjian (1987) [7] |
| NF-κB | `GGGACTTTCC` | Sen & Baltimore (1986) [8] |
| CREB | `TGACGTCA` | Montminy et al. (1986) [9] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported element's matched `Sequence` has length equal to its `Pattern` and occurs at `Position` (0-based) in the input. | scan emits `S[i..i+m)` only at matching `i` |
| INV-02 | Each matched `Sequence` IUPAC-matches its `Pattern`. | per-position IUPAC membership test |
| INV-03 | The scan is exhaustive: all matching start indices `0 <= i <= n-m` are reported. | linear window over every offset |
| INV-04 | Consensus strings equal the published values; no fabricated constants. | each constant cites its primary source [1]–[9] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` | required | DNA sequence to scan | non-null; A/C/G/T (case-insensitive, normalised to upper) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Name` | `string` | Human-readable element name (e.g. "TATA Box", "AP-1"). |
| `Position` | `int` | 0-based start index of the occurrence in the sequence. |
| `Sequence` | `string` | The matched substring (length = pattern length). |
| `Pattern` | `string` | The consensus pattern that matched (IUPAC). |
| `Description` | `string` | Short biological role label. |

Returns an `IEnumerable<RegulatoryElement>`; elements are yielded library-entry by library-entry, in increasing position within each entry.

### 3.3 Preconditions and Validation

Null `sequence` → `ArgumentNullException`. Empty sequence → empty result. Matching is case-insensitive (`DnaSequence` normalises to uppercase). Coordinates are 0-based, inclusive start. Alphabet is DNA; the input must already be a valid `DnaSequence` (validated by its constructor).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `sequence` is non-null.
2. For each library entry `(Name, Pattern, Description)`:
3. Run the IUPAC degenerate scan of `Pattern` over the sequence.
4. For each match, yield a `RegulatoryElement` with the name, 0-based position, matched substring, pattern, and description.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The element library (§2.2 table) is the reference table; each consensus string is a named constant in `MotifFinder.KnownMotifs` carrying an inline source citation. IUPAC ambiguity codes follow the standard nucleotide code (only `N` is used by the current library, in the E-box `CANNTG`).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Full scan | O(n × r × m̄) | O(1) extra | n = sequence length, r = number of library entries (12), m̄ = mean pattern length; matches `O(n × r)` in the registry for bounded m |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs)

- `MotifFinder.FindRegulatoryElements(DnaSequence)`: scans the consensus library and yields `RegulatoryElement` records.
- `MotifFinder.KnownMotifs`: nested static class of source-cited consensus constants.
- `MotifFinder.FindDegenerateMotif(DnaSequence, string)`: the underlying IUPAC per-position scan reused for each entry.

### 5.2 Current Behavior

Each library entry is scanned independently via `FindDegenerateMotif`, so results are grouped by entry rather than globally sorted by position. The suffix tree (`SuffixTree.FindAllOccurrences`) was **not** used: the E-box pattern `CANNTG` is IUPAC-degenerate (not a single exact substring), and the library is a small, fixed set of short patterns scanned once over typically short promoter regions, so a single linear IUPAC scan is the correct and simplest algorithm. The exact-substring suffix tree would require enumerating all 16 E-box expansions and offers no asymptotic benefit at this scale.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- All 12 consensus strings copied from their primary sources (§2.2): TATAAA [1], CCAAT [1], GGGCGG [3], TATAAT/TTGACA [2], GCCGCCACCATGG [4], AGGAGG [10], AATAAA [5], CANNTG [6], TGACTCA [7], GGGACTTTCC [8], TGACGTCA [9].
- Exhaustive 0-based exact/IUPAC occurrence reporting (INV-01..INV-03).

**Intentionally simplified:**

- NF-κB: scanned as the single strong reference κB site `GGGACTTTCC` rather than the full degenerate consensus `GGGRNWYYCC`; **consequence:** weaker/variant κB sites are not reported [8].
- Kozak: scanned as the single most-preferred-base string `GCCGCCACCATGG` rather than expanding the -3 purine / +4 G degeneracy; **consequence:** Kozak contexts that differ from the optimal string are not reported [4].

**Not implemented:**

- Position-weight-matrix scoring of partial/weak matches; **users should rely on:** `MotifFinder.CreatePwm` / `MotifFinder.ScanWithPwm` for score-based motif detection.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | AP-1 consensus corrected `TGAGTCA` → `TGACTCA` | Deviation (fix) | prior value reported wrong AP-1 sites and missed real ones | fixed | Lee, Mitchell & Tjian (1987) [7] |
| 2 | Added -10 (`TATAAT`) and -35 (`TTGACA`) prokaryotic hexamers | Deviation (addition) | prokaryotic promoters now detected | fixed | Harley & Reynolds (1987) [2] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null sequence | `ArgumentNullException` | contract |
| Empty sequence | empty result | no offset satisfies `0 <= i <= n-m` |
| Multiple occurrences of one element | all reported with their positions | INV-03 exhaustiveness |
| Degenerate E-box (`CACGTG`, `CAGCTG`, …) | matched as `CANNTG` | IUPAC `N` membership |
| Sequence with no consensus | empty result | no match |

### 6.2 Limitations

Detects only the fixed library of consensus strings; it is not a general motif discovery method and does not score partial matches, account for strand (only the given strand is scanned), or use spacing constraints between the -35 and -10 hexamers. NF-κB and Kozak use single representative strings (see §5.3). For weak/variant sites use PWM scanning.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var seq = new DnaSequence("GGGTATAAAGGG");
var hits = MotifFinder.FindRegulatoryElements(seq).ToList();
// hits[0]: Name="TATA Box", Position=3, Sequence="TATAAA", Pattern="TATAAA"
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MotifFinder_FindRegulatoryElements_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/MotifFinder_FindRegulatoryElements_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [MOTIF-REGULATORY-001-Evidence.md](../../../docs/Evidence/MOTIF-REGULATORY-001-Evidence.md)
- Related algorithms: [Shared_Motifs](Shared_Motifs.md), [Known_Motif_Search](../Motif_Analysis/Known_Motif_Search.md)

## 8. References

1. Bucher P. 1990. Weight matrix descriptions of four eukaryotic RNA polymerase II promoter elements derived from 502 unrelated promoter sequences. J Mol Biol 212(4):563-578. https://doi.org/10.1016/0022-2836(90)90223-9
2. Harley C.B., Reynolds R.P. 1987. Analysis of E. coli promoter sequences. Nucleic Acids Res 15(5):2343-2361. https://doi.org/10.1093/nar/15.5.2343
3. Lundin M., Nehlin J.O., Ronne H. 1994. Importance of a flanking AT-rich region in target site recognition by the GC box-binding zinc finger protein MIG1. Mol Cell Biol 14(3):1979-1985. https://doi.org/10.1128/mcb.14.3.1979
4. Kozak M. 1987. An analysis of 5'-noncoding sequences from 699 vertebrate messenger RNAs. Nucleic Acids Res 15(20):8125-8148. https://doi.org/10.1093/nar/15.20.8125
5. Proudfoot N.J., Brownlee G.G. 1976. 3' non-coding region sequences in eukaryotic messenger RNA. Nature 263:211-214. https://doi.org/10.1038/263211a0
6. Massari M.E., Murre C. 2000. Helix-loop-helix proteins: regulators of transcription in eucaryotic organisms. Mol Cell Biol 20(2):429-440. https://doi.org/10.1128/MCB.20.2.429-440.2000
7. Lee W., Mitchell P., Tjian R. 1987. Purified transcription factor AP-1 interacts with TPA-inducible enhancer elements. Cell 49(6):741-752. https://doi.org/10.1016/0092-8674(87)90612-X
8. Sen R., Baltimore D. 1986. Multiple nuclear factors interact with the immunoglobulin enhancer sequences. Cell 46(5):705-716. https://doi.org/10.1016/0092-8674(86)90346-6
9. Montminy M.R., Sevarino K.A., Wagner J.A., Mandel G., Goodman R.H. 1986. Identification of a cyclic-AMP-responsive element within the rat somatostatin gene. PNAS 83(18):6682-6686. https://doi.org/10.1073/pnas.83.18.6682
10. Shine–Dalgarno sequence. Wikipedia (citing primaries). https://en.wikipedia.org/wiki/Shine%E2%80%93Dalgarno_sequence
