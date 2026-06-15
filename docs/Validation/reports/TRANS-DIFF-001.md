# Validation Report: TRANS-DIFF-001 — Differential Expression (log2 fold change, Welch t-test, BH FDR)

- **Validated:** 2026-06-15   **Area:** Transcriptome
- **Canonical method(s):** `TranscriptomeAnalyzer.CalculateFoldChange(expression1, expression2)`, `TranscriptomeAnalyzer.FindDifferentiallyExpressed(genes, alpha, log2FoldChangeThreshold)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (test-quality defect found and fixed in-session)

## Stage A — Description

### Sources opened & what they confirm
All sources were retrieved/recomputed this session:

- **log2 fold change** — DESeq2 (Love, Huber & Anders 2014, Genome Biology 15:550) and the Science Park RNA-seq lesson define LFC as `log2(mean_treatment / mean_control)`; positive = up in the numerator (treatment / condition 2). The doc/spec convention (condition 2 = numerator = treatment) matches.
- **Welch's t-test** (Wikipedia citing Welch 1947): `t = (X̄₂ − X̄₁)/√(s²₁/N₁ + s²₂/N₂)` with unbiased (N−1) variances; Welch–Satterthwaite `ν = (s²₁/N₁+s²₂/N₂)² / [s⁴₁/(N²₁(N₁−1)) + s⁴₂/(N²₂(N₂−1))]`. Matches the doc formula verbatim.
- **Two-sided p-value** (Student's t-distribution CDF): `P(|T|≥t) = I_{ν/(ν+t²)}(ν/2, ½)`. Matches the doc.
- **Benjamini–Hochberg** (R `p.adjust(method="BH")`): `pmin(1, cummin(n/i · p[o]))[ro]` — sort p descending, multiply by n/rank (rank m..1), running minimum, clamp to 1, restore order. Matches the doc.
- **Two-criterion DE gate** (Science Park / DESeq2): significant ⇔ `|log2FC| ≥ threshold` AND `adjusted p < alpha`. Matches the doc.

### Formula check
Every formula in `docs/algorithms/Transcriptome/Differential_Expression.md` §2.2 reproduces the cited source exactly: the Welch t-statistic, the Welch–Satterthwaite df, the regularized-incomplete-beta two-sided tail, and the BH step-up. Sign convention (positive = up in condition 2) is consistent across the doc, spec, and code.

### Edge-case semantics check
The three degenerate conventions are documented and self-consistent:
- `<2` replicates per group → variance undefined → p = 1 (Assumption 2).
- se = 0 with equal means → p = 1; with unequal means → p = 0 (Assumption 3; limit of the t-statistic).
- zero mean → pseudocount c = 1 keeps `log2((m2+c)/(m1+c))` finite (Assumption 1).
These are reasonable degenerate-input conventions; none alters any non-degenerate (N≥2, finite-variance) output.

### Independent cross-check (numbers — recomputed this session with SciPy 1.13.1)
| Quantity | Source value | This session |
|---|---|---|
| log2(41/11) (M1) | 1.8981203859807863 | 1.8981203859807863 ✓ |
| Welch t, {1,2,3}↔{7,8,9} (M4) | 7.3484692283495345 | scipy.stats.ttest_ind = 7.3484692283495345 ✓ |
| Welch df (M4) | 4.0 | 4.0 ✓ |
| two-sided p (M4) | 0.0018262606682599833 | scipy ttest_ind & special.betainc = 0.0018262606682599833 ✓ |
| BH(0.001,0.4,0.5,0.9) (M5b) | (0.004, 0.6666667, 0.6666667, 0.9) | R BH algorithm = (0.004, 0.66666667, 0.66666667, 0.9) ✓ |

### Findings / divergences
- **Minor doc/Evidence cosmetic note:** Evidence §"Benjamini-Hochberg adjusted p-values" includes a first table keyed on (0.01,0.02,0.03,0.04) whose prose mixes the two examples; the operative sourced example used by the spec is (0.001,0.4,0.5,0.9)→(0.004,0.66667,0.66667,0.9), which is correct. Not a correctness defect — the formula and the canonical example are right.

Stage A: the biology/maths is correct in the abstract. **PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs`:
- `CalculateFoldChange` (L331–346): `log2((mean2+c)/(mean1+c))`, c=1; null/empty → 0.
- `FindDifferentiallyExpressed` (L368–411): per-gene LFC + Welch p + BH + two-criterion gate + regulation label.
- `WelchTTestTwoSidedPValue` (L420–450): N<2 → p=1; se=0 → equal means p=1 / unequal p=0; t, Welch–Satterthwaite df, `I_{df/(df+t²)}(df/2,½)`.
- `BenjaminiHochbergAdjust` (L459–480): R `p.adjust` BH step-up (cummin from largest p down, m/rank, clamp to 1).
- `RegularizedIncompleteBeta` / `BetaContinuedFraction` / `LogGamma` (L490–577): Numerical-Recipes betai (Lentz CF) + Lanczos lgamma.

### Formula realised correctly? (evidence)
Yes. Production output matches SciPy/R to ≥1e-9 for every checked quantity:
- M1 fold change = `Math.Log2(41/11)` exactly.
- M4 raw p = 0.0018262606682599833 = scipy.
- Helper genes A–D raw p match scipy ttest_ind(equal_var=False) to 1e-9 (now asserted), and their BH-adjusted p match the R BH algorithm to 1e-9 (now asserted).
- Note the older `CalculateTTestPValue`/`BenjaminiHochberg`/`AnalyzeDifferentialExpression` methods (a different, normal-approx-style overload) are not the registry-canonical methods for this unit and are out of scope here; the canonical `FindDifferentiallyExpressed`/`CalculateFoldChange` use the exact t-tail path.

### Cross-verification table recomputed vs code
All MUST/SHOULD/COULD values reproduce against the actual code (15→16 DE tests green; full suite 6501 passed).

### Variant/delegate consistency
Single canonical path; no `*Fast` variant. The DE result's `AdjustedPValue` with one gene equals raw p (m/rank=1) — confirmed (M4).

### Test quality audit (TEST-QUALITY GATE)
**Defect found (fixed in-session):** the original M5 test `FindDifferentiallyExpressed_BenjaminiHochberg_ReturnsExactAdjustedPValues` was a **code-echo test** — it read the production raw p-values back out (`results.Select(r => r.PValue)`) and compared the production adjusted p to an in-test reimplementation fed those same raw values. It would have passed against a wrong BH implementation as long as the in-test helper shared the bug, and its expected values traced to code output, not to a source. This violates the "sourced expectations, not code echoes" rule.

**Fix applied:**
1. Rewrote M5 to hard-code **externally-sourced** expectations: raw p from SciPy `ttest_ind(equal_var=False)` and BH-adjusted p from the R `p.adjust(method="BH")` algorithm, asserted to 1e-9. The test now fails if EITHER the t-test OR the BH step were wrong.
2. Added M5b `BenjaminiHochbergReference_CanonicalSourcedVector_MatchesRPAdjust` locking the canonical sourced vector (0.001,0.4,0.5,0.9) → (0.004, 0.6666667, 0.6666667, 0.9), the textbook R p.adjust example.
3. Added `FindDifferentiallyExpressed_ZeroVarianceSeparatedMeans_PValueIsZero` to cover the previously-untested se=0 / unequal-means branch (Assumption 3): p = 0, significant.

After the fix all Stage-A branches are exercised: fold change (UP/DOWN/FLAT/zero-mean/null), Welch p (exact value, N<2, se=0 equal, se=0 separated), BH (sourced exact values + canonical vector + monotone/≥raw invariants), two-criterion gate (strong/flat/below-threshold), regulation label (up/down/unchanged), empty/null input.

**Gate result: PASS** — exact sourced values, no green-washing, all branches covered, honest green (full unfiltered suite Failed: 0, build warning-free for changed file).

### Findings / defects
- One test-quality defect (code-echo BH test) — fixed. No implementation defect.

## Verdict & follow-ups
- **Stage A: PASS.** **Stage B: PASS-WITH-NOTES** (test-quality defect found and completely fixed in-session).
- **End-state: CLEAN.** Implementation is correct against SciPy/R; tests rewritten to lock externally-sourced values and to cover the previously-untested se=0 branch.
- Full suite: **6501 passed, 0 failed.** Build: 0 errors (pre-existing warnings unrelated).
- Follow-up (cosmetic, non-blocking): tidy the Evidence BH section's first illustrative table so it does not mix the (0.01,0.02,…) and (0.001,0.4,…) examples.
