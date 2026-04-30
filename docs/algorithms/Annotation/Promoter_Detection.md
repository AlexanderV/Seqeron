# Promoter Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Annotation |
| Test Unit ID | ANNOT-PROM-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Promoter detection in this repository is a bacterial motif-search helper, not a full transcription-start-site predictor. `GenomeAnnotator.FindPromoterMotifs(...)` scans DNA for exact `-35` and `-10` consensus substrings and a small set of derived prefix/suffix variants, then assigns a literature-based score to each matched motif. The method is deterministic and inexpensive, but it does not pair `-35` and `-10` elements into complete promoter models and does not tolerate mismatches outside the hard-coded motif library. It should therefore be read as motif annotation rather than full promoter prediction.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Bacterial promoters commonly contain two short sequence elements upstream of the transcription start site: a `-35` element and a `-10` element (Pribnow box). The current document identifies their consensus sequences as `TTGACA` for the `-35` box and `TATAAT` for the `-10` box, and it notes the classical spacing guideline of about `17` base pairs between them. The same source material also notes that natural promoters often match only `3` to `4` of the `6` consensus positions and that complete conservation is not the only determinant of promoter strength. The `-10` element is AT-rich, which is consistent with its role in local strand separation during transcription initiation.

### 2.2 Core Model

The score basis used by the repository comes from the E. coli position-specific nucleotide occurrence frequencies documented in the current source and references:

- `-35` box `TTGACA`: `T(0.69)`, `T(0.79)`, `G(0.61)`, `A(0.56)`, `C(0.54)`, `A(0.54)`
- `-10` box `TATAAT`: `T(0.77)`, `A(0.76)`, `T(0.60)`, `A(0.61)`, `A(0.56)`, `T(0.82)`

For a matched consecutive substring of one consensus element, the score is computed as:

$$
score = \frac{\sum p_i\text{ for matched consensus positions}}{\sum p_i\text{ for all six consensus positions}}
$$

where the denominator is the total weight of the full six-position consensus for the corresponding box.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported score lies in `[0, 1]` | Each score is a normalized fraction of a full-consensus weight |
| INV-02 | A full `TTGACA` or `TATAAT` match has score `1.0` | All six consensus positions are present in the numerator |
| INV-03 | Every reported motif belongs to the repository's fixed `-35` or `-10` variant library | The scan compares the sequence only against those hard-coded motif strings |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `dnaSequence` | `string` | required | DNA sequence scanned for promoter motifs | Empty input yields no hits; the implementation does not add an explicit null guard before uppercasing |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `position` | `int` | 0-based start index of the matched motif |
| `type` | `string` | Either `-35 box` or `-10 box` |
| `sequence` | `string` | Exact motif string that matched the input |
| `score` | `double` | Hard-coded probability-derived score for the matched motif |

### 3.3 Preconditions and Validation

`FindPromoterMotifs(...)` uppercases the input sequence before scanning, so lowercase and mixed-case DNA are accepted. An empty string yields no matches because the substring loops never execute. The method does not validate the alphabet beyond exact substring comparison, and it does not add an explicit null guard before calling `ToUpperInvariant()`. Matching is exact: bases outside the hard-coded motif library, spacing relationships between motifs, and mismatch-tolerant alternatives are not considered.

## 4. Algorithm

### 4.1 High-Level Steps

1. Convert the input DNA sequence to uppercase.
2. Scan the sequence for each supported `-35` motif variant and emit every exact hit with the corresponding hard-coded score.
3. Scan the sequence for each supported `-10` motif variant and emit every exact hit with the corresponding hard-coded score.
4. Return the union of all emitted motif hits, allowing overlapping and adjacent matches.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Motif Class | Variant | Score |
|-------------|---------|-------|
| `-35 box` | `TTGACA` | `1.000` |
| `-35 box` | `TTGAC` | `0.855` |
| `-35 box` | `TGACA` | `0.815` |
| `-35 box` | `TTGA` | `0.710` |
| `-10 box` | `TATAAT` | `1.000` |
| `-10 box` | `TATAA` | `0.801` |
| `-10 box` | `ATAAT` | `0.813` |
| `-10 box` | `TATA` | `0.665` |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindPromoterMotifs(...)` | `O(n × v)` | `O(1)` plus output | `n` = sequence length, `v` = 8 hard-coded motif variants; motif lengths are bounded by `4` to `6` |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomeAnnotator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs)

- `GenomeAnnotator.FindPromoterMotifs(string)`: Scans for exact `-35` and `-10` motif variants and emits scored hits.

### 5.2 Current Behavior

Repository-specific behavior confirmed by source and tests:

- The method scans `-35` motifs and `-10` motifs independently across the entire sequence.
- A single full consensus occurrence can produce multiple overlapping hits because the full motif, prefix-5, suffix-5, and prefix-4 variants are all reported separately when present.
- Returned positions are 0-based indexes into the original sequence.
- No spacing check is applied between `-35` and `-10` elements, so the method does not attempt to assemble complete promoter pairs.
- All score constants are hard-coded in source and match the values exercised by `ANNOT-PROM-001`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Recognition of the bacterial `-35` consensus `TTGACA` and the `-10` consensus `TATAAT`.
- Probability-derived weighting of supported motif variants using the E. coli occurrence frequencies documented in the current source and references.

**Intentionally simplified:**

- Only exact substring matches to the supported motif library are considered; **consequence:** naturally occurring motifs with other mismatch patterns are not reported.
- `-35` and `-10` motifs are scanned independently; **consequence:** the method can report individual boxes without establishing a biologically plausible promoter pair.
- Short 4- and 5-base variants are emitted directly; **consequence:** sensitivity increases, but the result set can contain more partial-motif hits than a stricter full-consensus search.

**Not implemented:**

- Validation of the canonical `17`-bp spacing between `-35` and `-10` elements; **users should rely on:** downstream filtering or custom post-processing.
- PWM-, HMM-, or sigma-factor-specific promoter models; **users should rely on:** no current alternative in this repository.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | The hard-coded `-10` score table follows the values cited in the repository's `Promoter (genetics)` source/comments rather than alternate summary values from the `Pribnow box` page | Assumption | Score constants are internally consistent across code and tests but may differ from another published summary table | accepted | `ANNOT-PROM-001` explicitly standardizes on the repository's chosen table |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns no motifs | The scan loops have no valid start positions |
| Null input | Not explicitly handled before uppercasing | The method dereferences `dnaSequence` immediately |
| Lowercase or mixed-case DNA | Handled identically to uppercase DNA | The input is uppercased before scanning |
| No supported motifs present | Returns an empty collection | Matching is exact against the fixed motif library |
| Overlapping or adjacent motifs | All supported hits are reported | Each motif variant is scanned independently |
| Full consensus occurrence | Also yields matching shorter variants at the corresponding offsets | Prefix/suffix substrings are part of the emitted motif library |

### 6.2 Limitations

This helper identifies motif occurrences, not complete promoters. It does not validate `-35` / `-10` spacing, does not tolerate arbitrary mismatches, and does not estimate transcriptional strength beyond the fixed variant scores. The current design is also specific to the bacterial motifs documented here and does not cover broader promoter architectures.

## 8. References

1. Wikipedia contributors. Promoter (genetics). https://en.wikipedia.org/wiki/Promoter_(genetics)
2. Wikipedia contributors. Pribnow box. https://en.wikipedia.org/wiki/Pribnow_box
3. Wikipedia contributors. TATA box. https://en.wikipedia.org/wiki/TATA_box
4. Pribnow D. Nucleotide sequence of an RNA polymerase binding site at an early T7 promoter. Proceedings of the National Academy of Sciences. 1975;72(3):784-788.
5. Harley CB, Reynolds RP. Analysis of E. coli promoter sequences. Nucleic Acids Research. 1987;15(5):2343-2361.
