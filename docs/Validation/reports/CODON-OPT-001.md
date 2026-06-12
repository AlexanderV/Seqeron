# Validation Report: CODON-OPT-001 — Sequence (Codon) Optimization

- **Validated:** 2026-06-12   **Area:** Codon Optimization
- **Canonical method(s):** `CodonOptimizer.OptimizeSequence(string, CodonUsageTable, OptimizationStrategy, double, double, double)`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Codon usage bias"** (https://en.wikipedia.org/wiki/Codon_usage_bias): codon
  optimization adjusts a gene's nucleotides to match the host's tRNA/codon preferences "without
  altering the amino acid sequence" — i.e. **synonymous substitution that preserves the encoded
  protein**. Modern variants also balance mRNA folding, codon-pair bias, ramp, harmonization.
- **Wikipedia — "Codon Adaptation Index"** / **Sharp & Li (1987)**: relative adaptiveness
  `w_i = f_i / max(f_j)` over the synonymous codons of an amino acid (range 0–1);
  `CAI = (∏_{i=1}^{L} w_i)^(1/L)` — geometric mean of weights, range (0, 1].
- **Plotkin & Kudla (2011), "Synonymous but not the same"** (Nat Rev Genet 12:32–42): codon
  changes are synonymous (protein unchanged) but affect expression/translation — confirms the
  protein-preservation invariant as the defining property.

### Formula / definition check
- Canonical optimization = replace each codon with the **highest-usage synonymous codon** of the
  same amino acid in the target table (MaximizeCAI). Confirmed against Sharp & Li / CAI definition.
- The critical invariant — `translate(optimize(seq)) == translate(seq)` — is exactly the
  definition of synonymous substitution; sources confirm it must hold.
- CAI formula in the implementation matches Sharp & Li exactly (`w_i = f_i/max f_j`,
  geometric mean); zero-frequency codons clamped to 1e-6 to avoid ln(0) on partial custom tables.

### Edge-case semantics
Empty → empty/CAI 0; single-codon AAs (Met AUG, Trp UGG) unchanged; stop codons preserved;
length not multiple of 3 trimmed to complete codons; T→U RNA conversion; case-insensitive.
All defined and sourced (standard genetic code, RNA notation).

### Independent cross-check (numbers)
- E. coli K12 (Kazusa 316407): Leu max = CUG (0.50), Arg max = CGC (0.40), Pro max = CCG (0.53),
  Ala max = GCG (0.36). Confirmed against the table embedded in the source.
- CAI of `CUAAGACGA` (E. coli): w = 0.04/0.50, 0.04/0.40, 0.06/0.40 →
  CAI = exp((ln0.08+ln0.10+ln0.15)/3) ≈ **0.1063**. Matches spec and test M_CAI.

**Stage A findings:** Description is biologically and mathematically correct. PASS.

## Stage B — Implementation

### Code path reviewed
`OptimizeSequence` (CodonOptimizer.cs:239–319), `SelectOptimalCodon` (321–359),
`BalanceGcContent` (378–414), `CalculateCAI` (423–450),
`CalculateRelativeAdaptiveness` (452–468), `SplitIntoCodons` (687–695), `TranslateCodon` (697–700).

### KEY correctness check — protein-preservation invariant
**Guaranteed by construction.** `SelectOptimalCodon` only ever returns a member of
`AminoAcidToCodons[aminoAcid]` (synonyms of the same AA) or the current codon. Stop codons
(`aminoAcid == "*"`) bypass selection entirely and are copied verbatim (lines 273–278).
The BalancedOptimization GC-balancing phase (`BalanceGcContent`) likewise only substitutes within
`AminoAcidToCodons[aminoAcid]` and skips stops. Unknown/invalid codons translate to "X", have no
synonym set, so `SelectOptimalCodon` returns them unchanged. Therefore every replacement is
synonymous and translation is invariant. Independently confirmed by the test
`OptimizeSequence_Invariant_ProteinPreserved` (re-translates original vs optimized with a separate
genetic-code map) across 5 inputs and `OptimizeSequence_PreservesProtein_AllStrategies` (all 5
strategies).

### Max-weight synonym rule
`MaximizeCAI`: `synonymousCodons.OrderByDescending(freq).First()` → highest-usage synonym. ✓
Confirmed: `CUAAGACGA` → `CUGCGCCGC` (CUA→CUG, AGA→CGC, CGA→CGC), CAI 1.0
(`OptimizeSequence_MaximizeCAI_ProducesOptimalCAI`).

### Worked example (default BalancedOptimization), recomputed vs code
Input `AUGGCUUAA` (M-A-*):
- AUG → Met, single synonym → AUG (unchanged).
- GCU → Ala; goodCodons (freq ≥ 0.15) = GCG(0.36), GCC(0.27), GCA(0.21) desc → **GCG**.
- UAA → stop → preserved.
Result `AUGGCGUAA`, protein `MA*`, 1 change (3, GCU, GCG), OriginalCAI exp(ln(0.16/0.36)/2) ≈ 0.667,
OptimizedCAI 1.0, GC 3/9 → 4/9. Exactly matches `OptimizeSequence_ReturnsValidOptimizationResult`.

### Start/stop handling
Start (AUG/Met) is single-codon → never changed. All three stops (UAA/UAG/UGA) preserved verbatim
and never substituted; confirmed by `OptimizeSequence_StopCodons_Preserved` and the
`AllStopCodons_CorrectProtein` cases.

### Edge cases in code
Empty → early return (line 247). Non-multiple-of-3 trimmed (254–258). T→U + uppercase (252).
Met/Trp single-synonym short-circuit (326–327). Invalid codon → unchanged (graceful).
`SplitIntoCodons` loop `i+2 < len` is equivalent to `i+3 <= len`, so it keeps all complete codons
(verified for len 3, 9). All covered by tests.

### Variant / strategy consistency
All 5 strategies preserve protein. `HarmonizeExpression` uses weighted-random selection
(`new Random()` per call — non-deterministic, but tests only assert protein preservation + CAI
range, which hold for any synonym). `AvoidRareCodons` replaces only sub-threshold codons; `Balanced`
adds GC balancing with a Changes-list rebuild (prior bug fixed 2026-03-10, now covered).

### Test quality audit
35 OptimizeSequence tests. They assert exact sourced values (specific optimized sequences per
Kazusa, hand-computed CAI/GC, exact change tuples), independently re-translate to lock the
protein-preservation invariant, and cover every Stage-A edge case. Real assertions, not tautologies.

**Stage B findings:** No defects. Code faithfully realises the validated description; the
protein-preservation invariant holds by construction and is independently tested. PASS.

## Verdict & follow-ups
- Stage A PASS, Stage B PASS. **State: CLEAN.**
- No code changes required. No logged defects.
- Minor (non-defect) note: `HarmonizeExpression` is stochastic (`new Random()` each call); this is
  intentional (matches host distribution) and does not affect any invariant or deterministic test.
