# Promoter Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Annotation |
| Test Unit ID | ANNOT-PROM-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Complete (consensus definition; delegates to MotifFinder) |
| Last Reviewed | 2026-11-25 (review-2026-09 B11, F7) |

## 1. Overview

`GenomeAnnotator.FindPromoterMotifs(...)` reports the −35 and −10 boxes of bacterial σ70 promoters. It does not
re-implement promoter search: it takes every forward-strand −35/−10 pair from the canonical
`MotifFinder.FindSigma70Promoters` (MOTIF-REGULATORY-001) and returns each box of those pairs once, as
`(position, type, sequence, score)` tuples. A box with no partner at a valid spacer is not a promoter and is not
reported. For the pairs themselves (strand, spacer, both mismatch counts) call `MotifFinder.FindSigma70Promoters`;
for a published quantitative promoter score (free energy, transcription rate) call
`MotifFinder.PredictSigma70Promoters` (Promoter Calculator v1.0, La Fleur, Hossain & Salis 2022).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

E. coli σ70 promoters carry two hexamers upstream of the transcription start site: the −35 element (consensus
`TTGACA`) and the −10 element / Pribnow box (consensus `TATAAT`). Harley & Reynolds (1987, 263 promoters) found
all bases of both hexamers highly conserved and an inter-region spacer of 17 ± 1 bp in 92 % of promoters
(Hawley & McClure 1983 compiled the earlier 112-promoter reference set). Natural promoters rarely match the full
consensus, so a promoter is recognised as a −35/−10 pair, each box within a few mismatches of its consensus,
separated by a spacer of 15–21 bp.

### 2.2 Core Model

For a sequence `S` (0-based), a σ70 pair is (i, j) with

- `d35 = Hamming(S[i, i+6), TTGACA) ≤ maxMismatches35`,
- `j = i + 6 + s`, `minSpacer ≤ s ≤ maxSpacer`,
- `d10 = Hamming(S[j, j+6), TATAAT) ≤ maxMismatches10`.

`FindPromoterMotifs` reports `(i, "-35 box", S[i,i+6), 1 − d35/6)` and `(j, "-10 box", S[j,j+6), 1 − d10/6)` for
every pair, deduplicated, ordered by position then −35 before −10. Defaults (single-argument overload) are those of
`FindSigma70Promoters`: ≤ 2 mismatches per box, spacer 15–21, forward strand.

The score is the fraction of the six consensus positions the box matches. It measures consensus identity only;
it is not a promoter-strength estimate (use `PredictSigma70Promoters` for that).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `score = 1 − d/6`, `d ≤` the box's mismatch limit | Computed from the Hamming distance `FindSigma70Promoters` reports |
| INV-02 | A full `TTGACA` / `TATAAT` box scores `1.0` | `d = 0` |
| INV-03 | Every −35 box has a −10 box 15–21 bp downstream (and vice versa) | Only boxes of pairs are reported |
| INV-04 | `0 ≤ position ≤ length − 6` and `S[position, position+6) = sequence` | Boxes are hexamers of the input |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `dnaSequence` | `string` | required | DNA sequence | A/C/G/T, case-insensitive; not null |
| `maxMismatches35` | `int` | 2 | Max mismatches of the −35 box to TTGACA | 0 … 6 |
| `maxMismatches10` | `int` | 2 | Max mismatches of the −10 box to TATAAT | 0 … 6 |
| `minSpacer` | `int` | 15 | Minimum spacer (bp) | ≥ 0 |
| `maxSpacer` | `int` | 21 | Maximum spacer (bp) | ≥ `minSpacer` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `position` | `int` | 0-based start of the hexamer |
| `type` | `string` | `-35 box` or `-10 box` |
| `sequence` | `string` | The hexamer (uppercase) |
| `score` | `double` | `1 − mismatches/6` |

### 3.3 Preconditions and Validation

Null → `ArgumentNullException`. A non-ACGT character → `ArgumentException` (input is wrapped in `DnaSequence`).
Mismatch limit outside 0 … 6, negative `minSpacer` or `minSpacer > maxSpacer` → `ArgumentOutOfRangeException`
(from `FindSigma70Promoters`). Empty or shorter than 6 + `minSpacer` + 6 → empty result. Validation is eager.

## 4. Algorithm

1. Wrap the input in `DnaSequence` (validates, uppercases).
2. Call `MotifFinder.FindSigma70Promoters(sequence, maxMismatches35, maxMismatches10, minSpacer, maxSpacer)`.
3. For each pair add its −35 and −10 boxes, with score `1 − d/6`, to a sorted set (dedup, order).
4. Return the set as a list.

| Operation | Time | Space |
|-----------|------|-------|
| `FindPromoterMotifs` | `O(n · (maxSpacer − minSpacer + 1))` | `O(boxes)` |

## 5. Implementation Notes

### 5.1 Location and Entry Points

[GenomeAnnotator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs):
`FindPromoterMotifs(string)` and `FindPromoterMotifs(string, int, int, int, int)`.
Canonical pairing: [MotifFinder.Sigma70Promoters.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.Sigma70Promoters.cs).
MCP: `find_promoter_motifs` (Annotation server) delegates to the single-argument overload.

### 5.2 Conformance to Theory / Spec

Implemented: Harley & Reynolds (1987) consensus pairing — TTGACA / TATAAT, spacer 15–21 bp (17 optimal),
mismatch-tolerant boxes. Not part of this method: strand (forward only; the tuple has no strand field — use
`FindSigma70Promoters(bothStrands: true)`), UP element / extended −10 / discriminator terms and promoter strength
(use `PredictSigma70Promoters`).

### 5.3 History

Before review-2026-09 (B11 F7) the method reported unpaired exact −35/−10 hexamer prefixes and suffixes
(`TTGAC`, `TGACA`, `TTGA`, `TATAA`, `ATAAT`, `TATA`) scored by a sum of per-position consensus percentages
(T69 T79 G61 A56 C54 A54 / T77 A76 T60 A61 A56 T82, from an uncited Wikipedia table). Neither the variant library
nor the score formula had a source, and the method did not apply the −35/−10 pairing that defines a σ70 promoter.

## 6. Edge Cases and Limitations

| Case | Behaviour |
|------|-----------|
| Empty / shorter than 27 nt (defaults) | Empty |
| Isolated −35 or −10 consensus | Empty (no pair) |
| Box paired with several partners | Reported once |
| Promoter on the reverse strand | Not reported (forward strand only) |
| Lowercase input | Same as uppercase |
| `N` or other non-ACGT | `ArgumentException` |

## 8. References

1. Harley CB, Reynolds RP. Analysis of E. coli promoter sequences. Nucleic Acids Res. 1987;15(5):2343-2361. doi:10.1093/nar/15.5.2343
2. Hawley DK, McClure WR. Compilation and analysis of Escherichia coli promoter DNA sequences. Nucleic Acids Res. 1983;11(8):2237-2255. doi:10.1093/nar/11.8.2237
3. Pribnow D. Nucleotide sequence of an RNA polymerase binding site at an early T7 promoter. PNAS. 1975;72(3):784-788.
4. La Fleur TL, Hossain A, Salis HM. Automated model-predictive design of synthetic promoters to control transcriptional profiles in bacteria. Nat Commun. 2022;13:5159 (via `PredictSigma70Promoters`).
