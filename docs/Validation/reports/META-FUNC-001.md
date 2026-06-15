# Validation Report: META-FUNC-001 — Functional Prediction (homology transfer + pathway over-representation)

- **Validated:** 2026-06-15   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.PredictFunctions`, `FindPathwayEnrichment` (+ helpers `Blosum62SelfScore`, `FunctionalBitScore`, `ExpectedValue`, `HypergeometricUpperTail`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (no code defect; three Stage-B test-quality defects fixed in-session)

## Context

This unit is a *simplified, specification-driven* realization of two standard metagenomic
functional-profiling steps (the HUMAnN / PICRUSt / eggNOG-mapper concept of gene → function/KO
transfer plus pathway abundance/enrichment): (1) homology-based annotation transfer scored with
Karlin-Altschul BLAST statistics, and (2) pathway over-representation analysis (ORA) via the
hypergeometric test. The matching step is an intentional exact-substring model (ASM-01), not a full
Smith-Waterman/BLAST search — the cited bit-score / E-value / hypergeometric *formulas* are used
exactly. Stage A validates those formulas + constants; Stage B validates the code and its tests.

## Stage A — Description

### Sources opened this session (URLs + extracted numbers)

1. **NCBI BLAST tutorial — "The Statistics of Sequence Similarity Scores"**
   (https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html, WebFetched 2026-06-15).
   Confirmed verbatim, three formulas:
   - E-value: `E = K·m·n·e^(−λS)` (Formula 1)
   - Bit score: `S' = (λS − ln K) / ln 2` (Formula 2)
   - E from bit score: `E = m·n·2^(−S')` (Formula 3)
   All three match the spec/Evidence/doc exactly.
2. **Ungapped BLOSUM62 Karlin-Altschul parameters.** WebSearch located `blast_stat.c`; values
   independently corroborated via the O'Reilly *BLAST* book Table 4-2 "Effect of gaps on BLOSUM62"
   (https://doctorlib.org/medical/blast/5.html, WebFetched): **no-gaps row → λ = 0.318, K = 0.134,
   H = 0.40 nats** — matches the full-precision `λ = 0.3176, K = 0.134, H = 0.4012` in the code.
3. **NCBI BLOSUM62 matrix file** (https://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62, WebFetched).
   Diagonal (self-match) scores confirmed: **W = 11, C = 9, A = 4, M = 5, R = 5, N = 6** — match the
   code's `Blosum62Diagonal` table.
4. **PNNL Proteomics Data Analysis §8.2 ORA** (https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html,
   WebFetched). Formula `P(X≥x) = 1 − Σ_{i=0}^{x−1} C(M,i)C(N−M,n−i)/C(N,n)` confirmed verbatim.
   Worked example: the **R call `phyper(q = 20−1, m = 400, n = 8000−400, k = 100, lower.tail = FALSE)`
   → 7.88 × 10⁻⁸**. (The page's prose labels for M/n are transposed relative to the R call; the R
   call is authoritative and maps to N=8000, M=400, n=100, x=20 — see cross-check below.)

### Formula check

Every formula in `Functional_Prediction.md` §2 and the TestSpec matches sources 1 and 4 verbatim,
including the `ln 2` divisor (bit units) and the right-tail-only ORA convention. The ungapped λ/K and
the BLOSUM62 diagonals match sources 2 and 3. INV-01..INV-06 are genuine properties (S' linear in S
with slope λ/ln2 > 0; E ∝ e^(−λS); algebraic identity of the two E forms; probability ∈ [0,1];
empty-sum/degenerate margins ⇒ p = 1; E-value ranks best hit).

### Edge-case semantics

- ORA x = 0, M = 0, or n = 0 ⇒ p = 1 (empty sum / degenerate margins) — sourced (PNNL + standard
  hypergeometric). Right-tail only (no lower-tail enrichment) — sourced.
- BLAST: a non-empty exact self-match always has S > 0 (all diagonals ≥ 4), so the EVD applies.

### Independent cross-check (numbers I recomputed this session)

| Quantity | Independent computation | Source value | Match |
|----------|------------------------|--------------|-------|
| Bit score, "WWW" S=33 | `(0.3176·33 − ln0.134)/ln2 = 18.020293278753364` | 18.0202932787533 | ✅ |
| E-value, S=33, m=n=3 (form 1) | `0.134·9·e^(−0.3176·33) = 3.3852730346545964e−5` | 3.3852730346546e−5 | ✅ |
| E-value (form 3) | `9·2^(−18.0202932787533) = 3.3852730346545937e−5` | identical | ✅ (INV-02) |
| ORA P(X≥20), N=8000,M=400,n=100 | log-Γ tail sum → `7.884747217109681e−8`; `scipy.stats.hypergeom(8000,400,100).sf(19) = 7.884747217146755e−8` | 7.88×10⁻⁸ | ✅ |
| ORA small exact, N=10,M=4,n=5,x=3 | `[C(4,3)C(6,2)+C(4,4)C(6,1)]/C(10,5) = 66/252 = 11/42 = 0.2619047619047619` | hand-computed | ✅ |
| Best-hit E ("WW" vs "AAAA" in "AAAAWW") | WW: m=6,n=2,S=22 → `1.4851955539388528e−3`; AAAA: m=6,n=4,S=16 → `1.997e−2` | WW lower ⇒ best | ✅ |

The `scipy` value 7.8847e-8 and the from-scratch log-Gamma tail sum agree, and both round to the
published 7.88e-8. The PNNL prose M/n label transposition is benign — the R `phyper` call (m = 400
successes, k = 100 draws, lower.tail = FALSE, q = 19) unambiguously fixes the mapping.

### Stage A findings

PASS. All formulas, constants (ungapped λ/K, BLOSUM62 diagonals) and edge-case semantics trace to
authoritative sources retrieved this session, with two fully independent cross-checks of the ORA
worth (scipy + hand-fraction) and a closed-form check of the BLAST numbers. No biological or
mathematical defect in the description.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs`:
- Constants `Blosum62UngappedLambda=0.3176`, `Blosum62UngappedK=0.134`, `Ln2`, `Blosum62Diagonal` (lines 718–733).
- `PredictFunctions` (754–802): per-gene scan, skips null/whitespace sequences and empty signatures,
  ordinal `Contains`, raw score via `Blosum62SelfScore`, `S'` via `FunctionalBitScore`,
  `E = K·m·n·e^(−λS)` via `ExpectedValue(rawScore, sequence.Length, signature.Length)`, keeps the
  lowest-E candidate (best hit). Null args → `ArgumentNullException`.
- `Blosum62SelfScore` (809–822): sum of diagonals, upper-cases residues, unknown ⇒ 0.
- `FunctionalBitScore` (828–829) and `ExpectedValue` (839–841): exactly Formulas 2 and 1.
- `FindPathwayEnrichment` (873–916): builds query/background sets (null background ⇒ union of
  members; query always added), per pathway M = |members ∩ background|, x = overlap with query,
  p = `HypergeometricUpperTail`, sorted ascending by p.
- `HypergeometricUpperTail` (929–947): upper tail summed directly in log-space (log-Γ) i = x..min(n,M),
  degenerate cases (x≤0 ∨ M≤0 ∨ n≤0 ∨ N≤0) ⇒ 1.0, clamped to [0,1]. Summing the upper tail directly
  (rather than 1 − lower-tail) correctly avoids catastrophic cancellation for tiny tails (7.88e-8).

### Formula realised correctly?

Yes — verbatim. The code reproduces Formulas 1/2/3 and the hypergeometric tail exactly; my from-scratch
Python reimplementation of `HypergeometricUpperTail` (same log-Γ scheme) reproduced 7.884747…e-8, and
`scipy` independently confirmed it. Bit/E literals (18.0202932787533, 3.3852730346546e-5) reproduced.

### Cross-verification vs code

Every value in the table above was reproduced by the actual code via the test suite (M1, M2, M3, M4,
plus the new small-exact ORA and `Blosum62SelfScore` tests). Variants/helpers agree with the canonical
methods (E from `ExpectedValue` equals the `PredictFunctions` field; INV-02 holds to literal-rounding).

### Test quality audit (HARD gate)

Found **three Stage-B test-quality defects (test-only; no code defect)** and fixed them in-session:

1. **M1/M2/M3 code-echo (green-washing risk).** The expected bit score and E-value were *recomputed*
   in the test from locally re-declared `Lambda`/`K` constants using the **same formula as the code**
   (`(Lambda*33 − ln K)/ln2`, `K*3*3*exp(−Lambda*33)`, `K*6*2*exp(−Lambda*22)`). A test that recomputes
   the implementation's own formula can pass against a same-shaped wrong implementation. **Fix:**
   asserted against the **sourced literals** instead — `18.0202932787533`, `3.3852730346546e-5`,
   and the hand-computed best-hit E `0.0014851955539388528` — and removed the now-unused mirrored
   constants. (M2 retains the INV-02 cross-form check, now at `Within(1e-17)` because the second form
   is built from the 13-sig-fig bit-score literal and the two forms coincide only to that rounding
   ~1.5e-18; the primary assert against the published E literal stays at `Within(1e-18)`.)

2. **No direct `Blosum62SelfScore` coverage.** Added `Blosum62SelfScore_SumsDiagonalScores`: W=11,
   "WACMRN"=40 (11+4+9+5+5+6, NCBI diagonals), ""=0, unknown residues ⇒ 0.

3. **ORA tested only by the single published example + degenerate cases; null-background union path
   untested.** Added (a) `HypergeometricUpperTail_SmallExactCase_MatchesHandComputation` — N=10,M=4,
   n=5,x=3 ⇒ exact **11/42 = 0.2619047619047619** (independent of the published example), and (b)
   `FindPathwayEnrichment_NullBackground_UsesUnionOfPathwayMembers` — exercises the 2-arg overload's
   union-of-members fallback (N=M=4 ⇒ p=1).

No assertion was weakened, no tolerance widened beyond what literal-rounding requires, nothing skipped.
M4 (`7.88e-8 Within 1e-9`) is correctly sourced: the published value has 3 sig figs and the tolerance
accepts the full-precision result while rejecting any wrong implementation by orders of magnitude.

Coverage now exercises: both `PredictFunctions` outcomes (match/no-match/empty/null), best-hit
selection, field transfer, both `FindPathwayEnrichment` overloads, all ORA degenerate margins
(x=0, M=0, n=0), the helpers `Blosum62SelfScore`/`FunctionalBitScore`/`ExpectedValue`/
`HypergeometricUpperTail`, and INV-01..INV-06.

### Findings / defects

- **No implementation defect.** Code realizes the validated formulas exactly.
- Three test-quality defects (code-echo expectations; missing `Blosum62SelfScore` test; thin ORA /
  untested union-background path) — **all fixed in-session** with sourced/hand-computed literals.
- Fixture grew 14 → 17 tests. Build 0 errors (4 pre-existing warnings in unrelated files, none in
  this fixture). **Full unfiltered suite: 6555 passed, 0 failed.**

## Verdict & follow-ups

- **Stage A: PASS.** **Stage B: PASS-WITH-NOTES** (no code defect; three test-quality defects fixed).
- **End-state: ✅ CLEAN.** Algorithm fully functional within its documented exact-substring scope
  (ASM-01); production-scale divergent-homolog detection and FDR correction remain out of scope by
  design, as documented in `Functional_Prediction.md` §5.3.
- **Test-quality gate: PASS** after the in-session fixes (honest green; sourced, non-echoing
  expectations; all public methods/overloads and Stage-A branches exercised).
