# Test Specification: PROBE-EVALUE-001

**Test Unit ID:** PROBE-EVALUE-001
**Area:** MolTools
**Algorithm:** Karlin–Altschul Off-Target E-value / Bit-Score
**Status:** ☑ Validated (Stage A PASS / Stage B PASS / CLEAN) — 2026-06-25
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

| # | Source | What it establishes |
|---|--------|---------------------|
| 1 | Karlin & Altschul (1990), PNAS 87:2264 | E = K·m·n·e^(−λS); λ = unique positive root of Σ p_i p_j e^{λ s_ij}=1 |
| 2 | Altschul et al. (1990), JMB 215:403 (BLAST) | local-alignment statistics; bit-score normalization |
| 3 | NCBI "Statistics of Sequence Similarity Scores" (Altschul), BLAST/tutorial/Altschul-1.html | verbatim S' = (λS − ln K)/ln 2; E = m·n·2^(−S'); negative-expected-score requirement |
| 4 | NCBI BLAST+ manual appendices (NBK279684) + swipe `blastkar_partial.c` | published ungapped nucleotide λ/K: 1/−3 → λ=1.374, K=0.711; 2/−3 → λ=0.55, K=0.21 |

## 2. Canonical Method(s)

- `ProbeDesigner.ComputeLambdaNucleotide(int match, int mismatch, double baseFrequency = 0.25)`
- `ProbeDesigner.ComputeKarlinAltschul(double rawScore, int queryLength, long databaseLength, ScoringMatrix? scoring = null, double k = 0.711, double baseFrequency = 0.25)`
- `KarlinAltschulStatistics` record struct.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs:1136-1286`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_ProbeValidation_Tests.cs` (KA1–KA12)

## 3. Contract / Invariants

- **R:** E-value ≥ 0; λ > 0 when defined.
- **Preconditions (λ undefined → throws):** match > 0 (≥1 positive score), mismatch < 0, expected per-pair score < 0.
- **M:** higher bit score / higher raw score → lower E (E monotonically decreasing in S).
- **M:** larger search space → higher E; E is **linear** in m, in n, and in K.
- **Identity:** E = K·m·n·e^(−λS) = m·n·2^(−S'), with S' = (λS − ln K)/ln 2.
- **Boundary:** S = 0 → E = K·m·n, S' = −ln K/ln 2.
- **Model note:** λ uses the **uniform-0.25** background; for a simple match/mismatch matrix this is the exact root of 0.25·e^{λ·match}+0.75·e^{λ·mismatch}=1. This reproduces NCBI's published 1/−3 λ=1.374, but NOT NCBI's published 2/−3 λ=0.55 (that value comes from the full score lattice). K has no citable closed form → caller parameter (default 0.711, the published 1/−3 value).

## 4. Cross-check / Differential Oracle

`blastn`/`conda` not installable in this environment → oracle = an independent Python bisection re-solve (this session) + NCBI's published constants. Exact numbers:

| Quantity | Oracle value |
|---|---|
| λ(1,−3, p=0.25) | 1.3740631224599755 (= NCBI 1.374) |
| λ(2,−3, p=0.25) | 0.6337314430979077 (default BlastDna; ≠ NCBI published 0.55) |
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
