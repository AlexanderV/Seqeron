# Validation Report: META-TETRA-001 — TETRA Tetranucleotide Z-score Signature

- **Validated:** 2026-06-25   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.CalculateTetranucleotideZScores`, `MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation`
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Canonical method(s)
`CalculateTetranucleotideZScores(string)`, `TetranucleotideZScoreCorrelation(string, string)`

- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs:674,706` (helpers `TetranucleotideZScore` :724, `ExtendWithReverseComplement` :758, `CountOligonucleotides` :791, `EnumerateTetranucleotides` :815, `PearsonCorrelation` :825).
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_TetranucleotideZScore_Tests.cs`

## Stage A — Description

### Sources opened this session
- **Teeling H, Waldmann J, Lombardot T, Bauer M, Glöckner FO (2004).** *TETRA: a web-service and a stand-alone program for the analysis and comparison of tetranucleotide usage patterns in DNA sequences.* BMC Bioinformatics 5:163. Open-access full text (PMC529438). Confirms in prose: 256 tetranucleotides counted; expected frequencies from a **maximal-order Markov model** from di-/trinucleotide composition; divergence converted to **z-scores using an approximation published by Schbath**; sequences **"extended by their reverse-complements"**; sequences compared **pairwise by the Pearson correlation coefficient of their z-scores**. (The BMC paper itself gives only prose; the explicit equations are in the companion Teeling et al. 2004 Environ Microbiol 6(9):938–947 and Schbath refs.)
- **Companion / corroborating sources** (BMC Genomics 2019 12864-019-6119-x; PLoS One 2009 0008113 genomic-signature literature) quote the explicit equations used by TETRA verbatim:
  - Expected: **E(n1n2n3n4) = N(n1n2n3)·N(n2n3n4) / N(n2n3)**
  - Variance (Schbath approximation): **var(n1n2n3n4) = E · [N(n2n3) − N(n1n2n3)]·[N(n2n3) − N(n2n3n4)] / N(n2n3)²**
  - Z-score: **Z(n1n2n3n4) = (N(n1n2n3n4) − E) / √var**
  - Comparison: Pearson r of the two 256-component z-vectors.

### Formula check
All four equations above match the cited sources symbol-for-symbol. The expected count is the standard maximal-order (2nd-order) Markov estimate; the variance is the Schbath (1995/1997) approximation for the count of an overlapping word under that Markov model. ✅

### Definitions & conventions
- Words counted on **overlapping** positions (standard). ✅
- Sequence is **extended by its own reverse complement** before counting → strand-symmetric signature. ✅
- Case-insensitive; non-ACGT characters dropped before counting. ✅
- Comparison is **Pearson** r ∈ [−1, 1]. ✅

### Edge-case semantics (sourced / mathematically defined)
- `N(n2n3)=0` → expected count undefined → z=0 (no signal). Mathematically well-defined.
- variance ≤ 0 → z=0. This occurs whenever a prefix/suffix trinucleotide count equals N(n2n3) (e.g. all counts = 1 on a very short strand), giving a degenerate (zero-variance) case → no over/under-representation signal. Correct.
- Empty / single-base / null → no 4-nt extended strand or all-degenerate → all-zero 256-map.
- Pearson of a zero-variance (all-zero) vector → defined as 0, not NaN.

### Independent cross-check (exact numbers, re-derived this session — NOT copied from repo)

**Invariant z(ACGT)=√5.** Input `"ACGTACGTGGCC"`; extended by its reverse complement →
`ACGTACGTGGCCGGCCACGTACGT` (24 nt). Overlapping counts (re-derived in pure Python):
N(ACGT)=4, N(ACG)=4, N(CGT)=4, N(CG)=5. Then
E = 4·4/5 = **3.2**; var = 3.2·(5−4)·(5−4)/5² = 3.2/25 = **0.128**;
z = (4 − 3.2)/√0.128 = 0.8/√0.128 = **2.2360679774997894 = √5**.
**Why √5 (algebraic, not numeric coincidence):** with N(ACGT)=N(ACG)=N(CGT)=k and N(CG)=k+1,
(N−E) = k − k²/(k+1) = k/(k+1); var = E·1·1/(k+1)² = k²/(k+1)³; so
z = [k/(k+1)] / √[k²/(k+1)³] = √(k+1). For k=4, z = √5. The √5 traces to the variance formula and the count structure, independent of the code.

**Pearson on a small worked example.** pearson([1,2,3,4],[2,4,6,8]) = **+1.0** (perfectly linear); pearson([1,2,3,4],[4,3,2,1]) = **−1.0**. Identical signatures → r = **1.0** (verified on `"ACGTACGTGGCCATGCATGCTTAA"`). r ∈ [−1,1] confirmed.

**Strand-symmetry invariant (re-derived).** On the extended strand s+rc(s), each occurrence of word w on one strand appears as rc(w) on the other; therefore the observed/prefix/suffix/middle counts feeding z(w) equal those feeding z(rc(w)). Hence **z(w) = z(rc(w)) for all 256 words** within a signature. Numerically verified: max|z(w)−z(rc(w))| = 0.0 (e.g. GGCC↔GGCC = 2.4495, ATGC↔GCAT = 1.3693). This is the precise meaning of "reverse-complement-merged counts identical".

### Findings / divergences (Stage A)
- **Minor doc imprecision (fixed):** the source XML doc claimed an input of ≥2 usable ACGT bases "yields … a non-trivial signature". Re-derivation shows `"AC"` → extended `"ACGT"` has every count = 1, so var = 0 and the whole signature is all-zero (trivial). The wording was corrected to describe the zero-variance fallback rather than asserting non-triviality. No algorithmic impact.

**Stage A verdict: ✅ PASS** — formulas, conventions, edge semantics, and invariants all match authoritative sources and independent derivation.

## Stage B — Implementation

### Code path reviewed
- `CalculateTetranucleotideZScores` (:674) → `ExtendWithReverseComplement` (:758, filters to upper ACGT, appends reverse complement) → `CountOligonucleotides` (:791, single-pass overlapping di/tri/tetra counts, skips non-ACGT) → per-word `TetranucleotideZScore` (:724).
- `TetranucleotideZScore` computes `expected = prefix3*suffix3/middle2`, `variance = expected*(middle2−prefix3)*(middle2−suffix3)/(middle2²)`, returns `(observed−expected)/√variance`; returns 0 when `middle2==0` or `variance<=0`. Matches the validated formulas exactly.
- `TetranucleotideZScoreCorrelation` (:706) builds both 256-vectors over the same fixed key order and calls `PearsonCorrelation` (:825), which returns 0 (not NaN) for a constant/zero vector.

### Formula realised correctly?
Yes — line-by-line correspondence with the Teeling/Schbath equations; integer counts promoted to `double` before division (no integer-truncation). ✅

### Cross-verification table recomputed vs code (independent Python oracle vs C# tests)

| Case | Independent value | C# result | Match |
|------|-------------------|-----------|-------|
| z(ACGT) for `ACGTACGTGGCC` | √5 = 2.2360679775 | 2.2360679775 (M-Z1) | ✅ |
| 256 keys returned | 256 | 256 (M-Z2) | ✅ |
| z(ACGT) for `AAAAAAAA` (N(CG)=0) | 0 | 0 (M-Z3) | ✅ |
| null/empty/`A`/`AC` → all-zero | all-zero | all-zero (M-Z4) | ✅ |
| noisy == clean | equal | equal (M-Z5) | ✅ |
| self-correlation | 1.0 | 1.0 (M-Z6) | ✅ |
| r_similar vs r_dissimilar | 0.6365 > ~0 | greater (M-Z7) | ✅ |
| symmetry corr(a,b)=corr(b,a) | equal | equal (M-Z8) | ✅ |
| z(w)=z(rc(w)) ∀w | max diff 0.0 | equal (M-Z9, added) | ✅ |
| corr vs empty | 0 (not NaN) | 0 (S-Z1) | ✅ |

### Variant/delegate consistency
`TetranucleotideZScoreCorrelation` reuses `CalculateTetranucleotideZScores`; both share `EnumerateTetranucleotides` key ordering → vectors align component-wise. Consistent. ✅

### Numerical robustness
Division guarded by `middle2==0` and `variance<=0`; counts promoted to double; Pearson guards zero denominator. No overflow/NaN on the stated ranges. ✅

### Test quality audit
- Expected values trace to the Teeling/Schbath formula and hand-derivation (√5 derived algebraically), not to code echoes. No green-washing; assertions check exact sourced values within tight tolerances.
- Coverage of every Stage-A path: z(ACGT)=√5 (M-Z1), 256-dim (M-Z2), absent middle dinucleotide / all-same-base AAAA (M-Z3), seq<4 bp & null/empty/single (M-Z4), non-ACGT + case (M-Z5), self-corr=1 / identical→1 (M-Z6), discriminative ordering (M-Z7), symmetry (M-Z8), zero-variance→0-not-NaN (S-Z1).
- **Gap closed this session:** the required "reverse-complement-merged counts identical" invariant was not asserted. Added **M-Z9** (`CalculateTetranucleotideZScores_IsReverseComplementSymmetric`) asserting z(w)=z(rc(w)) for all 256 words — a source-derived, code-independent invariant.

### Findings / defects (Stage B)
- No correctness defect. One test-coverage gap (RC-symmetry invariant) — fixed by adding M-Z9. One source XML-doc imprecision — fixed.

**Stage B verdict: ✅ PASS** — code faithfully realises the validated description; suite is real and now covers all required paths.

## Verdict & follow-ups
- **State: ✅ CLEAN.** Description and implementation independently validated against Teeling (2004) and the Schbath variance approximation; z(ACGT)=√5 re-derived from first principles; Pearson behaviour confirmed on hand examples; strand-symmetry invariant proven and now tested.
- Changes this session: added test **M-Z9**; corrected a minor source doc string. Full unfiltered `dotnet test Seqeron.sln -c Debug` = **Failed: 0** (Seqeron.Genomics.Tests: 18761 passed), 0 build warnings on changed files.
- No findings requiring follow-up beyond the (already-applied) doc fix and added test.
