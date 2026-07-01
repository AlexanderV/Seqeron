# Test Specification: META-TETRA-001

**Test Unit ID:** META-TETRA-001
**Area:** Metagenomics
**Algorithm:** TETRA Tetranucleotide Z-score Signature
**Status:** ☑ Validated — independent Stage A/B re-validation complete (2026-06-25), State CLEAN
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Teeling H, Waldmann J, Lombardot T, Bauer M, Glöckner FO (2004). *TETRA: a web-service and a stand-alone program …* BMC Bioinformatics 5:163 (open-access PMC529438) — method, reverse-complement extension, Pearson comparison. |
| 2 | Teeling H et al. (2004). *Application of tetranucleotide frequencies for the assignment of genomic fragments.* Environ Microbiol 6(9):938–947 — the explicit z-score equations. |
| 3 | Schbath S (1995/1997) — variance approximation for overlapping-word counts under a maximal-order Markov model. |
| 4 | Corroborating verbatim equations: BMC Genomics 2019 12864-019-6119-x; PLoS One 0008113 (genomic-signature literature). |

## 2. Canonical Method(s)

`CalculateTetranucleotideZScores(string)`, `TetranucleotideZScoreCorrelation(string, string)`

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs` (:674, :706; helpers :724/:758/:791/:815/:825)
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_TetranucleotideZScore_Tests.cs`

## 3. Formulas (validated against sources)

- Expected count: **E(n1n2n3n4) = N(n1n2n3)·N(n2n3n4) / N(n2n3)**
- Variance (Schbath): **var = E · [N(n2n3) − N(n1n2n3)]·[N(n2n3) − N(n2n3n4)] / N(n2n3)²**
- Z-score: **Z = (N(n1n2n3n4) − E) / √var**
- Sequence extended by its reverse complement before counting (strand-symmetry).
- Comparison: **Pearson r** of the two aligned 256-component z-vectors.

## 4. Contract / Invariants

- R: correlation ∈ [−1, 1]; identical signatures → r = 1; zero-variance vector → r = 0 (not NaN).
- INV (z(ACGT)=√5): for `ACGTACGTGGCC`, z(ACGT) = √5 = 2.2360679775 — derived from E=3.2, var=0.128.
- INV (strand symmetry): z(w) = z(reverse-complement(w)) for all 256 words.
- Degenerate: N(n2n3)=0 or var≤0 → z=0. Null/empty/single-base → all-zero 256-map.
- D: deterministic.

## 5. Cross-check / Differential Oracle

- **Reference:** independent Python reimplementation of the Teeling/Schbath equations (this session).
- **Numbers:** z(ACGT)=√5; Pearson ±1 on linear vectors; r_similar 0.6365 > r_dissimilar ≈0; strand-symmetry max diff 0.0. Comparison tolerance z-vector ±1e-10.

## 6. Test Map (fixture → spec)

| Test | Asserts |
|------|---------|
| M-Z1 | z(ACGT)=√5 (hand-derived formula) |
| M-Z2 | 256-component signature, keyed by 4-mers |
| M-Z3 | absent middle dinucleotide / all-same base → z=0 |
| M-Z4 | null/empty/single-base → all-zero |
| M-Z5 | non-ACGT + case ignored |
| M-Z6 | self-correlation = 1.0 |
| M-Z7 | similar > dissimilar (discriminative) |
| M-Z8 | correlation symmetric |
| M-Z9 | strand symmetry z(w)=z(rc(w)) ∀w |
| S-Z1 | correlation vs empty = 0 (not NaN) |

## 7. Validation Checklist (restored ☑)

- [x] Stage A: every source retrieved this session; formulas confirmed against publications + corroborating literature.
- [x] Stage B: implementation reviewed line-by-line against the formulas; cross-checked vs independent oracle.
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0 (Genomics 18761 passed).
- [x] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
