# Test Specification: PROBE-EVALUE-001

**Test Unit ID:** PROBE-EVALUE-001
**Area:** MolTools
**Algorithm:** Karlin–Altschul Off-Target E-value / Bit-Score
**Status:** ☑ Validated (B07 re-review 2026-10-01: Stage A ⚠ → fixed / Stage B ❌ → fixed; see docs/Validation/review-2026-09/B07.md F29–F31)
**Last Updated:** 2026-10-01

---

## 1. Evidence Summary

| # | Source | What it establishes |
|---|--------|---------------------|
| 1 | Karlin & Altschul (1990), PNAS 87:2264 | E = K·m·n·e^(−λS); λ = unique positive root of Σ p_i p_j e^{λ s_ij}=1 |
| 2 | Altschul et al. (1990), JMB 215:403 (BLAST) | local-alignment statistics; bit-score normalization |
| 3 | NCBI "Statistics of Sequence Similarity Scores" (Altschul), BLAST/tutorial/Altschul-1.html | verbatim S' = (λS − ln K)/ln 2; E = m·n·2^(−S'); negative-expected-score requirement |
| 4 | NCBI BLAST+ `blast_stat.c` / `blast_setup.c` / `blast_hits.c` (ncbi-cxx-toolkit-public, raw.githubusercontent.com) + NCBI blastn 2.12.0+ | ungapped λ/K/H computed (`BlastKarlinLHtoK`): 1/−3 → 1.374/0.711/1.31; 2/−3 → 0.634/0.408/0.912. Gapped `blastn_values_*` tables (2/−3 gap 5/2 → 0.625/0.41/0.78, α 0.8, β −2, even round-down). The 2/−3 {0,0} row 0.55/0.21 is the non-affine megablast entry, NOT the ungapped value (the earlier row here was wrong). Length adjustment `BLAST_ComputeLengthAdjustment`; effective search space `(m−ℓ)(n−Nℓ)` |

## 2. Canonical Method(s)

- `ProbeDesigner.ComputeLambdaNucleotide(int match, int mismatch, double baseFrequency = 0.25)`
- `ProbeDesigner.ComputeKarlinAltschul(double rawScore, int queryLength, long databaseLength, ScoringMatrix? scoring = null, double? k = null, double baseFrequency = 0.25)` — k null → K computed for the scheme.
- `KarlinAltschulStatistics` record struct.
- `ProbeDesigner.ComputeUngappedKarlinParameters(int, int, double = 0.25)` / `(int, int, IReadOnlyList<double>)` → `KarlinAltschulParameters` (λ, K, H, α, β, RoundDown, Gapped).
- `ProbeDesigner.GetBlastnGappedKarlinParameters(int reward, int penalty, int gapOpen, int gapExtend)`.
- `ProbeDesigner.ComputeLengthAdjustment(double k, double alphaOverLambda, double beta, int m, long n, int N = 1)`.
- `ProbeDesigner.ComputeBlastnStatistics(int rawScore, int m, long n, int N = 1, ScoringMatrix? = BlastDna, bool gapped = true)` and `(string probe, string subject, ScoringMatrix?)` (aligns with `SequenceAligner.LocalAlignAffine`) → `BlastnStatistics`.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs:1136-1286`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_ProbeValidation_Tests.cs` (KA1–KA21)

## 3. Contract / Invariants

- **R:** E-value ≥ 0; λ > 0 when defined.
- **Preconditions (λ undefined → throws):** match > 0 (≥1 positive score), mismatch < 0, expected per-pair score < 0.
- **M:** higher bit score / higher raw score → lower E (E monotonically decreasing in S).
- **M:** larger search space → higher E; E is **linear** in m, in n, and in K.
- **Identity:** E = K·m·n·e^(−λS) = m·n·2^(−S'), with S' = (λS − ln K)/ln 2.
- **Boundary:** S = 0 → E = K·m·n, S' = −ln K/ln 2.
- **Model note:** λ uses the **uniform-0.25** background; for a simple match/mismatch matrix this is the exact root of 0.25·e^{λ·match}+0.75·e^{λ·mismatch}=1. This reproduces the ungapped λ NCBI blastn prints for every scheme (1/−3 → 1.374, 2/−3 → 0.634). K is computed by the BLAST+ `BlastKarlinLHtoK` lattice formula on the gcd-reduced lattice (scale-invariant; default `KarlinKMethod.ReducedLattice`); `KarlinKMethod.NcbiBlast` reproduces the non-invariant K blastn 2.12 prints for gcd > 1 schemes (series indexed from the unreduced array: 4/−6 → 1.17, 4/−10 → 1.06, 6/−4 → 1.63) and is the default of `ComputeBlastnStatistics` / `GetBlastnGappedKarlinParameters` (ungapped and infinite-gap-domain blocks) — tests KA22–KA24 (audit round 3, A3-17, B07 F54); a caller K overrides it. Gapped scores use the BLAST+ tables via `ComputeBlastnStatistics`.

## 4. Cross-check / Differential Oracle

2026-10-01: NCBI blastn 2.12.0+ installed (apt `ncbi-blast+`) and a line-by-line Python port of blast_stat.c (λ NR, LtoH, LHtoK, length adjustment) that reproduces blastn's Lambda/K/H, effective search spaces and E-values exactly — see docs/Evidence/PROBE-VALID-001-Evidence.md "NCBI blastn 2.12.0+ statistics". Earlier (2026-06) numbers:

| Quantity | Oracle value |
|---|---|
| λ(1,−3, p=0.25) | 1.3740631224599755 (= NCBI 1.374) |
| λ(2,−3, p=0.25) | 0.6337314430979077 (default BlastDna; = blastn ungapped 0.634) |
| K(1,−3) / K(2,−3) computed | 0.7106027952162398 / 0.4081456625463167 (blastn 0.711 / 0.408) |
| blastn m 40, n 3079, 2/−3 5/2, S 80 | ℓ 11, space 88972, E 7.035793990394873e-18, 73.42105622960482 bits |
| Bit S' (S=30, K=0.711) | 59.962700114285006 |
| E (S=30, m=20, n=1000, K=0.711) | 1.7801583686083893e-14 |
| S=0 → E / S' | 14220 / 0.4920785350426718 |

## 5. Validation Checklist

- [x] Stage A: every source retrieved this session; formula & constants confirmed against NCBI/Karlin–Altschul.
- [x] Stage A independent cross-check: numerical re-solve of λ + worked example match code exactly.
- [x] Stage B: implementation reviewed against source; bisection well-posed; all values reproduced by code.
- [x] Test quality: KA1–KA12 trace to external oracle; coverage gaps (score 0, K param, m-monotonicity, root convergence, default scheme) closed; no green-washing.
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0; 0 warnings on changed file.
- [x] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the `docs/checklists/*.md`.
