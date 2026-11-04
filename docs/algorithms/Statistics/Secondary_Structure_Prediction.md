# Protein Secondary Structure Prediction (Chou-Fasman propensity profile and assignment)

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics / Protein sequence analysis |
| Test Unit ID | SEQ-SECSTRUCT-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Implemented (profile + full Chou-Fasman 1978 assignment; opt-in graded 1978 rule set) |
| Last Reviewed | 2026-10-14 |

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
secondary-structure segments [3][4][8]. `PredictSecondaryStructure` computes the windowed
mean-propensity profile; `PredictSecondaryStructureChouFasman` implements the discrete
assignment rules of §2.3.

### 2.3 Chou-Fasman assignment rules (`PredictSecondaryStructureChouFasman`)

Rules 1-3 as stated by Chen et al. (2006) "Methods" [8], plus the 1978 β-turn rule [1][9]:

1. **Nucleation.** Helix: any 6-residue window with ≥ 4 helix formers (Pα > 1.00).
   Sheet: any 5-residue window with ≥ 3 sheet formers (Pβ > 1.00).
2. **Extension.** Each nucleus is extended one residue at a time in both directions while the
   tetrapeptide formed by the new residue and the three adjacent segment residues has mean
   propensity ≥ 1.00 ("until the average 4-peptide propensity drops below 1" [8]).
3. **Acceptance.** Extended helix kept if ⟨Pα⟩ > 1.03 and ⟨Pα⟩ > ⟨Pβ⟩; extended strand kept
   if ⟨Pβ⟩ > 1.05 and ⟨Pβ⟩ > ⟨Pα⟩ [8]. Accepted segments of one type are unioned.
4. **Overlap.** Each maximal run of residues covered by both a helix and a strand is assigned
   helix if ⟨Pα⟩ > ⟨Pβ⟩ over the run, otherwise strand ("the conformation with higher average
   propensity" [8]).
5. **β-turn.** Tetrapeptide i..i+3 is a turn if p(t) = f(i)·f(i+1)·f(i+2)·f(i+3) > 7.5×10⁻⁵,
   ⟨Pt⟩ > 1.00 and ⟨Pα⟩ < ⟨Pt⟩ > ⟨Pβ⟩ [1][9]. All four residues are marked T; turns take
   precedence over helix/strand. Other residues are coil (C).

Bend frequencies f(i), f(i+1), f(i+2), f(i+3) (Chou & Fasman 1978 [1] via [9][10]):

| AA | f(i) | f(i+1) | f(i+2) | f(i+3) | AA | f(i) | f(i+1) | f(i+2) | f(i+3) |
|----|------|--------|--------|--------|----|------|--------|--------|--------|
| A | 0.060 | 0.076 | 0.035 | 0.058 | M | 0.068 | 0.082 | 0.014 | 0.055 |
| R | 0.070 | 0.106 | 0.099 | 0.085 | F | 0.059 | 0.041 | 0.065 | 0.065 |
| N | 0.161 | 0.083 | 0.191 | 0.091 | P | 0.102 | 0.301 | 0.034 | 0.068 |
| D | 0.147 | 0.110 | 0.179 | 0.081 | S | 0.120 | 0.139 | 0.125 | 0.106 |
| C | 0.149 | 0.050 | 0.117 | 0.128 | T | 0.086 | 0.108 | 0.065 | 0.079 |
| E | 0.056 | 0.060 | 0.077 | 0.064 | W | 0.077 | 0.013 | 0.064 | 0.167 |
| Q | 0.074 | 0.098 | 0.037 | 0.098 | Y | 0.082 | 0.065 | 0.114 | 0.125 |
| G | 0.102 | 0.085 | 0.190 | 0.152 | V | 0.062 | 0.048 | 0.028 | 0.053 |
| H | 0.140 | 0.047 | 0.093 | 0.054 | I | 0.043 | 0.034 | 0.013 | 0.056 |
| L | 0.061 | 0.025 | 0.036 | 0.070 | K | 0.055 | 0.115 | 0.072 | 0.095 |

All comparisons are done on the integer parameters (×100), so thresholds such as ⟨Pα⟩ > 1.03
are exact (Σ Pα > 103·L).

### 2.3a Graded 1978 rule set (`ChouFasmanRuleSet.Graded1978`, opt-in)

`PredictSecondaryStructureChouFasman(sequence, ChouFasmanRuleSet.Graded1978)` replaces step 1 (and
constrains step 2 for helices) with the graded-class rules of Chou & Fasman [1][2] as summarised by
Prevelige & Fasman [11] (abbreviated rules 1-4, Table I footnotes) and CFSSP / MIT 20.320 [12]; steps 3-5
are unchanged. Class assignments for the 1978 parameters (Prevelige & Fasman Table II [11], agrees with [5]):

| Class | Helix (Pα) | Sheet (Pβ) |
|-------|------------|------------|
| H strong former | E M A L | V I Y |
| h former | K F Q W I V | F W L C T Q M |
| I weak former | D H | — (no residue) |
| i indifferent | R T S C | R N H A |
| b breaker | Y N | S G K |
| B strong breaker | P G | P D E |

1. **Graded nucleation.** Helix: 6-residue window with H/h counted 1 and Iα counted ½ (≥ 4 formers;
   "Iα count as 0.5 hα" [12]) and at most one bα/Bα breaker ("unfavourable if ⅓ or more breakers" [12];
   "< 2 breakers" [13]). Sheet: 5-residue window, ≥ 3 Hβ/hβ formers (Iβ ½ — empty class for the 1978
   table), at most one bβ/Bβ breaker.
2. **Boundary residues.** "Iα assignments are also given to Pro and Asp (near the N-terminal helix) as well
   as Arg (near the C-terminal helix); bβ assignment is also given to Trp (near the C-terminal β region)"
   [11, Table I footnotes b, c] — the class form of "Pro, Asp⁻, Glu⁻ prefer the N-terminal helical end;
   His⁺, Lys⁺, Arg⁺ the C-terminal end" [12] (Asp/His are already Iα, Glu Hα, Lys hα in 1978).
   Convention: "near the N-/C-terminal end" = the first/last three residues of the helix window (the first
   helical turn, where Pro is permitted) and the last two residues of the sheet window.
3. **Proline rule.** "Proline cannot occur in the inner helix or at the C-terminal helical end but can occur
   within the last three residues at the N-terminal end" [11]: a helix window with Pro at offset ≥ 3 cannot
   nucleate, C-terminal extension never adds Pro, and N-terminal extension stops when it would move a Pro
   beyond the third helix residue.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A window of one residue returns that residue's (Pα, Pβ, Pt) tuple. | Mean of one value is that value [1]. |
| INV-02 | Each emitted tuple is the per-component arithmetic mean of its window's residues. | Definition in §2.2. |
| INV-03 | For a length-n all-known sequence, the number of emitted windows is max(0, n − w + 1). | Step-1 sliding scan [4]. |
| INV-04 | Output is case-insensitive. | Input is upper-cased before lookup (ASCII-only fold: U+017F 'ſ' and U+0131 'ı' are non-residues like 'X', EMBOSS convention — 2026-10-14, B03 F45). |
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

`PredictSecondaryStructureChouFasman(string proteinSequence)` — same alphabet/case rules;
returns a string of the input's length with one state per residue: `H` helix, `E` strand,
`T` β-turn, `C` coil; empty for null/empty input. Non-standard residues have no parameters:
they are never formers, any window/tetrapeptide containing one is ineligible, and they are `C`.

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
- `SequenceStatistics.PredictSecondaryStructureChouFasman(string)`: discrete per-residue
  Chou-Fasman assignment (H/E/T/C), §2.3. Public constants `ChouFasmanHelix/Sheet/Turn/Coil`.
- `SequenceStatistics.PredictSecondaryStructureChouFasman(string, ChouFasmanRuleSet)`:
  `Standard` (= one-argument overload) or `Graded1978` (§2.3a). Undefined enum value →
  `ArgumentOutOfRangeException`.

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

- Discrete Chou-Fasman assignment (`PredictSecondaryStructureChouFasman`): nucleation 4/6
  and 3/5, tetrapeptide extension, 1.03/1.05 acceptance with ⟨Pα⟩ vs ⟨Pβ⟩, overlap
  resolution, β-turn p(t) > 7.5e-5 rule with 1978 bend frequencies [1][8][9][10].

- Graded 1978 classes (Hα/hα/Iα/iα/bα/Bα, Hβ/hβ/iβ/bβ/Bβ) with fractional counting, the
  breaker limit, the Pro rule and the boundary-residue reassignments, opt-in via
  `ChouFasmanRuleSet.Graded1978` (§2.3a) [1][11][12][13]. The default stays the [8] formulation.

**Not implemented:** the hand "data reduction" judgements of [11] §V (including terminal class-i
residues by eye, preferring helix lengths that are multiples of four, resolving overlapping turns by
the higher local p(t)); they are explicitly left to the user by the authors and have no algorithmic form.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Lysine Pα = 1.14 vs 1.16 | Deviation | Shifts any window mean containing K | fixed | Adopted 1.14 (2 sources vs 1); corrected from prior 1.16 |
| 2 | Default windowSize = 7 | Assumption | Default window is not a Chou-Fasman constant | accepted | Caller parameter; nucleation windows are 6/5 |
| 3 | Unknown-residue handling | Assumption | All-unknown window emits nothing | accepted | Skip-and-exclude; no source specifies otherwise |
| 4 | Extension continues while tetrapeptide mean ≥ 1.00 | Convention | Boundary at exactly 1.00 | accepted | [8] "until … drops below 1"; ravihansa3000 uses > 1.00 |
| 5 | Overlap tie ⟨Pα⟩ = ⟨Pβ⟩ → strand | Convention | Ties only | accepted | Literal rule "helix if ⟨Pα⟩ > ⟨Pβ⟩" as in ravihansa3000 and kalliapap implementations |
| 6 | Turn marks all four residues; turns override H/E | Convention | Turn span | accepted | A β-turn is a tetrapeptide [1]; kalliapap marks 4 residues and applies turns last |
| 7 | Unknown residues in assignment are coil and block windows | Assumption | Sequences with X/B/Z | accepted | No parameters exist for them [1] |
| 8 | Graded1978 boundary zones: first/last 3 residues of the helix window, last 2 of the sheet window | Convention | Pro/Arg/Trp reassignments | accepted | Sources say only "near the N-/C-terminal end" [11]; 3 = first helical turn where Pro is allowed [11] |
| 9 | Graded1978 Pro rule applied per extended segment | Convention | Adjacent helices whose union contains an inner Pro | accepted | [11] rule 1 is stated for a single helix |

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
| Non-ASCII letters that `ToUpperInvariant` folds to a residue ('ſ' → S, 'ı' → I) | Treated as 'X' (no parameters; coil in the assignment) | ASCII-only fold, as every other protein path of `SequenceStatistics` (F22) and EMBOSS |

### 6.2 Limitations

Modest accuracy (~50-60% Q3); parameters derive from a small 29-protein 1970s sample [3].
The profile does not assign discrete elements (use `PredictSecondaryStructureChouFasman`).
Neither method models long-range or tertiary context.

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

- Tests: [SequenceStatistics_PredictSecondaryStructure_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceStatistics_PredictSecondaryStructure_Tests.cs) — covers `INV-01`..`INV-06`
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
9. Grinnell ExBioPy Project 7.5 turn rule and prowl.rockefeller.edu Chou-Fasman table, as reproduced in ravihansa3000/ChouFasman (`ChouFasman.py`, by residue name) — https://raw.githubusercontent.com/ravihansa3000/ChouFasman/master/ChouFasman.py (accessed 2026-09-28)
10. hassan11196/Chou-Fasman (`DS PROJECT Chou Fashman.cpp`, same table by one-letter code) — https://github.com/hassan11196/Chou-Fasman (git clone, 2026-09-28)
11. Prevelige P Jr, Fasman GD. 1989. Chou-Fasman prediction of the secondary structure of proteins: the Chou-Fasman-Prevelige algorithm. In: Fasman GD (ed.) Prediction of Protein Structure and the Principles of Protein Conformation, Plenum, pp. 391-416 (Table I, Table II, rules 1-4). PDF: https://cfpred.sourceforge.net/preveligeandfasman1989.pdf (accessed 2026-10-14)
12. Kumar TA. 2013. CFSSP: Chou and Fasman Secondary Structure Prediction server (rule statements "Iα count as 0.5 hα", "⅓ or more breakers"); MIT OCW 20.320 Problem Set 2 (2009), Criteria 1-3 — via search-engine extracts (sites not reachable from the sandbox), 2026-10-14
13. Maclin R, Shavlik J. UCI "protein-secondary-structure.theory" (Chou-Fasman domain theory: ≥ 4 formers with two high-indifferent = one former, < 2 breakers; sheet ≥ 3 of 5, < 2 breakers) — https://raw.githubusercontent.com/graingert/PyPred/master/protein-secondary-structure/protein-secondary-structure.theory (accessed 2026-10-14)
