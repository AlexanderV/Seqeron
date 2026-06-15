# Validation Report: TRANS-EXPR-001 — Expression Quantification (TPM, FPKM/RPKM, Quantile Normalization)

- **Validated:** 2026-06-15   **Area:** Transcriptome
- **Canonical method(s):** `TranscriptomeAnalyzer.CalculateTPM`, `CalculateFPKM`, `QuantileNormalize`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** ✅ CLEAN

## Stage A — Description

### Sources opened this session (independent of the repo)

| # | Source | Retrieved | What it confirmed |
|---|--------|-----------|-------------------|
| 1 | Zhao, Ye & Stanton (2020), RNA 26(8) — PMC7373998 | WebFetch https://pmc.ncbi.nlm.nih.gov/articles/PMC7373998/ | Verbatim `RPKM = 10⁹ × reads / (Total reads × Transcript length)`; verbatim `TPM = 10⁶ × (reads/length) / Σ(reads/length)`; `TPM = 10⁶ × RPKM / Σ RPKM`; "The average TPM is equal to 10⁶ divided by the number of annotated transcripts … and thus is a constant" (⇒ within-sample sum = 10⁶). |
| 2 | Wikipedia "Quantile normalization" (cites Bolstad et al. 2003) | WebFetch https://en.wikipedia.org/wiki/Quantile_normalization | Full worked example: input C1=(5,2,3,4), C2=(4,1,4,2), C3=(3,4,6,8); rank means 2.00/3.00/4.67/5.67; tie rule verbatim; final matrix with both tied 4s in C2 → 5.17. |

Both fetched pages reproduce the formulas and the worked example exactly as the repo's Evidence doc claims, so the Evidence/TestSpec are corroborated by primary retrieval (not merely echoed).

### Formula check
- **TPM** `= (X_i/l_i)/Σ_j(X_j/l_j)·10⁶` — matches Source 1 verbatim.
- **FPKM/RPKM** `= X_i·10⁹/(l_i·N)` — matches Source 1 verbatim (`10⁹·reads/(total·length)`).
- **Quantile normalization** — sort each column, replace rank r by across-sample mean at rank r, place rank means back at original positions, with the Bolstad tie-average rule for ties — matches Source 2 verbatim.

### Edge-case semantics
- All-zero counts → Σ(X/l)=0 → TPM 0/0 undefined; convention emits 0 (ASSUMPTION-01, declared, not from literature — acceptable degenerate convention).
- length ≤ 0 or N ≤ 0 → FPKM undefined → 0 (degenerate convention).
- Empty input → empty output. All defined and reasonable.

### Independent cross-check (recomputed by hand / Python this session)
- TPM M1: RPK=(0.005,0.005,0.030), Σ=0.04 → TPM=(125000,125000,750000), Σ=10⁶. ✓
- TPM M3 (2 genes, equal RPK 0.005): Σ=0.01 → each TPM = **500000** (the TestSpec had erroneously written 125000 — corrected this session). ✓
- FPKM M4: 1000·10⁹/(2000·10⁶)=500. ✓
- QN rank means r0..r3 = 2, 3, 14/3, 17/3; tie mean (r2,r3)=31/6=5.166… → both tied 4s in C2 = 5.17. ✓ (matches fetched Wikipedia final matrix).

### Findings
Description is biologically and mathematically correct and matches primary sources verbatim. **PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs`
- `CalculateTPM` L112–145 — rate `X/l`, sum, divide·10⁶; all-zero guard; also fills FPKM using total raw count as N.
- `CalculateFPKM` L155–163 — `X·10⁹/(l·N)`; returns 0 for length ≤ 0 or N ≤ 0.
- `QuantileNormalize` L175–232 — rank means via sorted columns; tied-run averaging (Bolstad rule); original positions preserved.

### Formula realised correctly?
Yes. Each method computes the validated formula exactly (named constants `TpmScalingFactor=1e6`, `FpkmScalingFactor=1e9`). Rate uses `Math.Max(length,1)` to avoid div-by-zero on degenerate zero-length input; division is double (no integer truncation since `RawCount` is `double`). Tie-run loop groups equal values and assigns the average of their spanned rank means — exactly the Bolstad rule.

### Cross-verification table recomputed vs code
All M/S/C expected values (M1=125000/125000/750000; M4=500; M6/M7 = 14/3,17/3,31/6 etc.) match the actual code output (full suite green) and the externally fetched sources.

### Variant/delegate consistency
`CalculateTPM` reuses `CalculateFPKM` internally; FPKM-in-TPM values (A=B=83333.33…, C=500000 for N=60) re-derived from the RPKM formula and now locked by a new test.

### Test quality audit (HARD gate)
Defects found in the **tests** (not the algorithm), all fixed this session:
1. **M3 was a code-echoing/weak assertion** — it only asserted `result[0].TPM == result[1].TPM` (equality to each other), which a deliberately-wrong-but-consistent impl would pass. Rewritten to assert the exact sourced value **500000** for both. The TestSpec's stated expected value (125000) was also wrong and was corrected.
2. **FPKM field populated inside `CalculateTPM` was untested** (documented behaviour in §3.2/§5.2). Added `CalculateTPM_PopulatesFpkmField_UsingTotalCountAsDepth` with exact sourced values (83333.33…, 83333.33…, 500000; N=60).
3. **INV-05 (untied columns are permutations of the same rank-mean multiset) was listed but not directly tested.** Added `QuantileNormalize_UntiedColumns_ArePermutationsOfSameRankMeanMultiset` asserting columns 1 and 3 sorted = {2, 3, 14/3, 17/3}.

No assertions weakened, no tolerances widened, no tests skipped. Every expected value traces to a source fetched this session. Coverage now spans all three public methods, all formula paths, the all-zero/empty/non-positive edge cases, the tie rule, and INV-01..INV-05.

### Findings / defects
No implementation defect. Test-quality gaps fixed. **PASS-WITH-NOTES** (notes = the test/spec corrections above; the algorithm itself was already correct).

## Verdict & follow-ups
- **Stage A: PASS.** Formulas and worked example confirmed verbatim against PMC7373998 and Wikipedia/Bolstad fetched this session.
- **Stage B: PASS-WITH-NOTES.** Code realises the validated formulas exactly; tests strengthened (M3 exact value, FPKM-in-TPM branch, INV-05) and TestSpec M3 value corrected.
- **End-state: ✅ CLEAN.** `dotnet build` 0 errors; full unfiltered suite **6499 passed, 0 failed**.
- **Test-quality gate: PASS** (after fixes) — sourced expectations, no green-washing, all logic covered, honest green on the full suite.
