# Protein Secondary Structure Prediction (Chou-Fasman propensity profile)

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics / Protein sequence analysis |
| Test Unit ID | SEQ-SECSTRUCT-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

This algorithm scores a protein sequence for its local tendency to adopt α-helix,
β-sheet, or β-turn conformation using the Chou-Fasman conformational propensities
Pα, Pβ and Pt [1][2]. For each sliding window it reports the mean propensity of the
residues in the window, producing three parallel profiles (helix, sheet, turn) along
the chain. It is a heuristic, statistics-based predictor: a window mean above 1.0
indicates that the segment, on average, favours the corresponding conformation. It is
useful for quick, interpretable secondary-structure profiling; it is not a full
structure predictor and is known to be of modest accuracy (~50-60% Q3) [3].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Chou and Fasman analysed proteins of known structure and counted how often each of the
20 amino acids occurs in α-helix, β-sheet, and β-turn versus its overall occurrence,
yielding a conformational *propensity* (observed/expected) per residue per
conformation [1][4]. A propensity > 1 means the residue is over-represented in that
conformation ("former"); < 1 means under-represented ("breaker") [4][5].

### 2.2 Core Model

For residue *r*, let Pα(*r*), Pβ(*r*), Pt(*r*) be its helix, sheet, and turn
propensities [1]. For a window of residues w = (r₁ … r_k) the profile value for a
conformation *c* is the arithmetic mean of the member propensities:

> mean_c(w) = ( Σ_{r ∈ w, r known} P_c(r) ) / (number of known residues in w)

The full Chou-Fasman method additionally applies nucleation/extension rules
(helix: 4 of 6 contiguous formers, extend until 4 contiguous breakers; sheet: 3 of 5;
turn: product of position-specific bend frequencies above a cutoff) to assign discrete
secondary-structure segments [3][4][8]. The implementation here computes the windowed
mean-propensity profile only (see §5.3).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A window of one residue returns that residue's (Pα, Pβ, Pt) tuple. | Mean of one value is that value [1]. |
| INV-02 | Each emitted tuple is the per-component arithmetic mean of its window's residues. | Definition in §2.2. |
| INV-03 | For a length-n all-known sequence, the number of emitted windows is max(0, n − w + 1). | Step-1 sliding scan [4]. |
| INV-04 | Output is case-insensitive. | Input is upper-cased before lookup. |
| INV-05 | Unknown residues are excluded from a window's count and mean; an all-unknown window emits nothing. | Table defines only the 20 standard residues [5][7]. |
| INV-06 | Null/empty input, w > n, or w < 1 → empty result, no exception. | Validated precondition (§3.3). |

### 2.5 Comparison with Related Methods

| Aspect | Chou-Fasman propensity profile (this) | GOR method |
|--------|----------------------------------------|------------|
| Information used | Single-residue propensities, windowed mean | 17-residue window, pairwise conditional probabilities [6] |
| Output | Continuous per-window propensity profile | Per-residue conformation class |
| Accuracy (Q3) | ~50-60% [3] | Higher (information-theory based) [6] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| proteinSequence | string | required | Amino-acid sequence, one-letter code | Case-insensitive; non-standard residues skipped |
| windowSize | int | 7 | Sliding-window length in residues | ≥ 1 and ≤ sequence length, else empty result |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (Helix, Sheet, Turn) | IEnumerable of (double, double, double) | One mean-propensity tuple per window position, N-terminus → C-terminus. |

### 3.3 Preconditions and Validation

Input is upper-cased (case-insensitive). The accepted alphabet is the 20 standard
amino acids; any other character (X, B, Z, `*`, gaps, digits) is silently skipped and
excluded from the window average. Null or empty input, a window larger than the
sequence, or a window size below 1 yields an empty enumerable rather than an exception
(INV-06). Indexing is 0-based; windows step by exactly one residue.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input; return empty on null/empty, w > n, or w < 1.
2. Upper-case the sequence.
3. For each start index i from 0 to n − w, sum Pα/Pβ/Pt over the w residues, skipping
   unknown residues and counting only known ones.
4. If the window has ≥ 1 known residue, emit (helixSum, sheetSum, turnSum) ÷ count.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Chou-Fasman conformational parameters (propensity = published integer ÷ 100) [1][6][7]:

| AA | Pα | Pβ | Pt | AA | Pα | Pβ | Pt |
|----|----|----|----|----|----|----|----|
| A | 1.42 | 0.83 | 0.66 | M | 1.45 | 1.05 | 0.60 |
| R | 0.98 | 0.93 | 0.95 | F | 1.13 | 1.38 | 0.60 |
| N | 0.67 | 0.89 | 1.56 | P | 0.57 | 0.55 | 1.52 |
| D | 1.01 | 0.54 | 1.46 | S | 0.77 | 0.75 | 1.43 |
| C | 0.70 | 1.19 | 1.19 | T | 0.83 | 1.19 | 0.96 |
| E | 1.51 | 0.37 | 0.74 | W | 1.08 | 1.37 | 0.96 |
| Q | 1.11 | 1.10 | 0.98 | Y | 0.69 | 1.47 | 1.14 |
| G | 0.57 | 0.75 | 1.56 | V | 1.06 | 1.70 | 0.50 |
| H | 1.00 | 0.87 | 0.95 | I | 1.08 | 1.60 | 0.47 |
| L | 1.21 | 1.30 | 0.59 | K | 1.14 | 0.74 | 1.01 |

Lysine Pα is **1.14**: two retrieved sources [6][7] give 1.14 while one [5] gives 1.16;
the 1.14 majority is adopted (see Evidence Assumption 1).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Profile of length-n sequence, window w | O(n·w) | O(1) extra (streamed) | Lazily yields one tuple per position; table lookup is O(1). |

The unit is O(n·w) ≤ O(n²) and not a "find occurrences of X in Y" search, so the
repository suffix tree is not applicable (no substring search is performed).

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.PredictSecondaryStructure(string, int)`: yields per-window mean
  (Helix, Sheet, Turn) Chou-Fasman propensity tuples.

### 5.2 Current Behavior

The method is a lazy iterator (`yield return`) that streams one tuple per window
position. Unknown residues are skipped per window; a window containing only unknown
residues produces no output for that position. No substring search is performed, so the
repository suffix tree was evaluated and is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The Chou-Fasman Pα/Pβ/Pt parameter table for all 20 residues [1][5][6][7].
- Windowed mean-propensity scoring (the averaging step used during nucleation-region
  evaluation, §2.2) [4].

**Intentionally simplified:**

- Discrete secondary-structure assignment: only the continuous windowed mean profile is
  produced; **consequence:** the user gets propensity profiles, not labelled helix/sheet/turn
  segments.

**Not implemented:**

- Helix/sheet nucleation-and-extension state machine (4-of-6, 3-of-5, extend-until-breaker)
  and the β-turn position-frequency product rule [3][4][8]; **users should rely on:** no
  current in-repo alternative — interpret the profile directly or use a dedicated tool (GOR,
  PSIPRED) for discrete assignment.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Lysine Pα = 1.14 vs 1.16 | Deviation | Shifts any window mean containing K | fixed | Adopted 1.14 (2 sources vs 1); corrected from prior 1.16 |
| 2 | Default windowSize = 7 | Assumption | Default window is not a Chou-Fasman constant | accepted | Caller parameter; nucleation windows are 6/5 |
| 3 | Unknown-residue handling | Assumption | All-unknown window emits nothing | accepted | Skip-and-exclude; no source specifies otherwise |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null / empty sequence | Empty result | Precondition (INV-06) |
| windowSize > length | Empty result | No window fits [4] |
| windowSize < 1 | Empty result | Invalid window (INV-06) |
| Unknown residue inside window | Excluded from mean | Table covers 20 residues only (INV-05) |
| Window of only unknown residues | No tuple emitted | count = 0 (INV-05) |
| Lower-case input | Same as upper-case | Input upper-cased (INV-04) |

### 6.2 Limitations

Modest accuracy (~50-60% Q3); parameters derive from a small 29-protein 1970s sample [3].
Does not assign discrete secondary-structure elements, does not model long-range or
tertiary context, and does not implement the turn position-frequency product rule.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// "AE" with window 2 -> one tuple of the residue-pair means.
var profile = SequenceStatistics.PredictSecondaryStructure("AE", windowSize: 2).Single();
// profile.Helix = (1.42 + 1.51)/2 = 1.465
// profile.Sheet = (0.83 + 0.37)/2 = 0.60
// profile.Turn  = (0.66 + 0.74)/2 = 0.70
```

**Numerical walk-through:** For "AEV" with window 2, two windows are emitted:
window [A,E] → ((1.42+1.51)/2, (0.83+0.37)/2, (0.66+0.74)/2) = (1.465, 0.60, 0.70);
window [E,V] → ((1.51+1.06)/2, (0.37+1.70)/2, (0.74+0.50)/2) = (1.285, 1.035, 0.62).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_PredictSecondaryStructure_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_PredictSecondaryStructure_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [SEQ-SECSTRUCT-001-Evidence.md](../../../docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md)
- Related algorithms: [Hydrophobicity_Analysis](./Hydrophobicity_Analysis.md)

## 8. References

1. Chou PY, Fasman GD. 1978. Empirical predictions of protein conformation. Annual Review of Biochemistry 47:251-276. https://pubmed.ncbi.nlm.nih.gov/354496/
2. Chou PY, Fasman GD. 1974. Prediction of protein conformation. Biochemistry 13(2):222-245. https://pubmed.ncbi.nlm.nih.gov/4358940/
3. Wikipedia. Chou–Fasman method. https://en.wikipedia.org/wiki/Chou%E2%80%93Fasman_method (accessed 2026-06-13)
4. Kelley bioinfo. Protein 2° Structure: Chou-Fasman Algorithm. https://www.kelleybioinfo.org/algorithms/background/BCho.pdf (accessed 2026-06-13)
5. Jakubowski H. Chou-Fasman propensities (CSB|SJU CH331). https://employees.csbsju.edu/hjakubowski/classes/ch331/protstructure/tablechoufas.htm (accessed 2026-06-13)
6. Przytycka T. Protein secondary structure prediction (NCBI/NLM lecture). https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/download/lectures/CAMS_02_Prot_Sec_Str.pdf (accessed 2026-06-13)
7. ravihansa3000. ChouFasman reference implementation. https://raw.githubusercontent.com/ravihansa3000/ChouFasman/master/ChouFasman.py (accessed 2026-06-13)
8. Chen H, Gu F, Huang Z. 2006. Improved Chou-Fasman method for protein secondary structure prediction. BMC Bioinformatics 7(Suppl 4):S14. https://pmc.ncbi.nlm.nih.gov/articles/PMC1780123/
