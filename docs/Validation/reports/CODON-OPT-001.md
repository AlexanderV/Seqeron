# Validation Report: CODON-OPT-001 — Sequence (Codon) Optimization

- **Validated:** 2026-06-24   **Area:** Codon Optimization
- **Canonical method(s):** `CodonOptimizer.OptimizeSequence(string, CodonUsageTable, OptimizationStrategy, double gcTargetMin, double gcTargetMax, double rareCodonThreshold)`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs:239`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Codon Adaptation Index"** (fetched 2026-06-24): relative adaptiveness
  `w_i = f_i / max(f_j)` where `max(f_j)` is the frequency of the most frequent synonymous codon
  of that amino acid; `CAI = (∏_{i=1}^{L} w_i)^(1/L)` (geometric mean of weights). Cites
  Sharp & Li (1987), *Nucleic Acids Research* 15(3):1281–1295. Matches the spec/Evidence exactly.
- **Kazusa Codon Usage Database, species=316407** (E. coli W3110 K-12), fetched 2026-06-24:
  per-thousand and per-AA fraction columns used to re-derive the built-in E. coli table (below).
- Codon optimization = synonymous substitution replacing each codon with the highest-usage
  synonymous codon of the same amino acid for the target organism, *without altering the encoded
  protein* (Wikipedia "Codon usage bias"; Plotkin & Kudla 2011). The defining invariant is
  `translate(optimize(seq)) == translate(seq)`.

### Formula check
- CAI formula (`CalculateRelativeAdaptiveness` + geometric mean in `CalculateCAI`) reproduces
  Sharp & Li (1987) `w_i = f_i/max f_j`, `CAI = exp((1/L)·Σ ln w_i)` exactly.
- Canonical optimization (MaximizeCAI) = max-frequency synonymous codon per residue — confirmed.

### Edge-case semantics
Empty → empty result, CAI 0; single-codon AAs Met (AUG) and Trp (UGG) unchanged; stop codons
(UAA/UAG/UGA) preserved verbatim; length not a multiple of 3 trimmed to complete codons; T→U RNA
conversion; case-insensitive (uppercased). All defined and sourced (standard genetic code, RNA
notation).

### Independent cross-check (numbers) — Kazusa 316407 re-derived by hand
Per-AA fraction = (codon per-thousand) / (Σ synonymous per-thousand):
- **Leu** Σ=106.4: CUG 53.1/106.4 = 0.499 → 0.50 ✓; CUA 3.8/106.4 = 0.036 → 0.04 ✓;
  UUA 0.130 → 0.13 ✓; UUG 0.128 → 0.13 ✓; CUU 0.103 → 0.10 ✓; CUC 0.104 → 0.10 ✓.
- **Arg** Σ=55.3: CGC 22.3/55.3 = 0.403 → 0.40 ✓; CGU 0.380 → 0.38 ✓; AGA 0.036 → 0.04 ✓;
  AGG 0.020 → 0.02 ✓; CGA 0.063 → 0.06 ✓; CGG 0.098 → 0.10 ✓.
- **Ala** Σ=94.9: GCG 33.9/94.9 = 0.357 → 0.36 ✓; GCU 0.160 → 0.16 ✓; GCC 0.271 → 0.27 ✓;
  GCA 0.212 → 0.21 ✓.
- **Pro** Σ=44.3: CCG 23.4/44.3 = 0.528 → 0.53 ✓.
- **Thr** Σ=53.6: ACC 23.5/53.6 = 0.438 → 0.44 ✓; ACU 0.16 ✓; ACG 0.27 ✓; ACA 0.13 ✓.
Every built-in E. coli value tested matches Kazusa to the stored 2-decimal precision.

Hand-computed CAI (Python), matching the test assertions exactly:
- `CUAAGACGA` (E. coli) = exp((ln0.08+ln0.10+ln0.15)/3) = **0.1063**.
- `CUGCCGACC` (Human) = exp((ln1+ln(0.11/0.32)+ln1)/3) = **0.7005**.
- `AUGGCUUAA` original (E. coli) = exp(ln(0.16/0.36)/2) = **0.6667**.

### Findings / divergences
None affecting the canonical model. The CAI clamp to `1e-6` for codons *absent from the table*
diverges from Sharp & Li's 0.5/N convention, but it triggers only on incomplete custom tables;
all three built-in tables are complete, so it never fires for the canonical organisms. Documented.

## Stage B — Implementation

### Code path reviewed
`OptimizeSequence` (CodonOptimizer.cs:239–319), `SelectOptimalCodon` (321–359),
`BalanceGcContent` (378–414), `CalculateCAI` (423–450), `CalculateRelativeAdaptiveness` (452–468),
`SplitIntoCodons` (687–695), `TranslateCodon` (697–700).

### Protein-preservation invariant — holds by construction
`SelectOptimalCodon` only ever returns a member of `AminoAcidToCodons[aminoAcid]` (synonyms of the
same AA) or the current codon. Stop codons (`aminoAcid == "*"`) bypass selection and are copied
verbatim (273–278). The BalancedOptimization GC-balancing phase (`BalanceGcContent`) also only
substitutes within `AminoAcidToCodons[aminoAcid]` and skips stops. Unknown codons translate to "X"
(no synonym set) and are returned unchanged. Therefore every substitution is synonymous and
translation is invariant. Independently locked by `OptimizeSequence_Invariant_ProteinPreserved`
(re-translates original vs optimized via a separate genetic-code map, 5 inputs) and
`OptimizeSequence_PreservesProtein_AllStrategies` (all 5 strategies).

### Max-weight synonym rule
`MaximizeCAI`: `synonymousCodons.OrderByDescending(freq).First()` → highest-usage synonym. ✓
`CUAAGACGA` → `CUGCGCCGC` (CUA→CUG, AGA→CGC, CGA→CGC), CAI 1.0 — asserted by
`OptimizeSequence_MaximizeCAI_ProducesOptimalCAI`.

### Cross-verification table recomputed vs code (tests executed)
`CodonOptimizer_OptimizeSequence_Tests` + `CodonOptimizer_CAI_Tests`: **59 passed, 0 failed.**
- `AUGGCUUAA` BalancedOpt → `AUGGCGUAA`, 1 change (3,GCU,GCG), OrigCAI 0.667, OptCAI 1.0,
  GC 3/9→4/9 — matches my hand computation.
- E. coli vs Yeast `CUGAGA`: E. coli → `CUGCGC`, Yeast → `UUGAGA` (verified against both Kazusa
  tables: Leu best CUG 0.50 / UUG 0.29; Arg best CGC 0.40 / AGA 0.48).
- Human vs E. coli CAI of `CUGCCGACC`: 1.0 vs 0.700 — matches.

### Variant/delegate consistency
All 5 strategies preserve protein. `HarmonizeExpression` uses weighted-random selection
(`new Random()` per call → non-deterministic) but only over synonyms, so the invariant holds and
the test asserts only protein + CAI range. `AvoidRareCodeons` replaces only sub-threshold codons.
`MinimizeSecondary` falls through to BalancedOptimization in `SelectOptimalCodon`; the dedicated
`ReduceSecondaryStructure` method is separate (documented in TestSpec deviations).

### Numerical robustness
Geometric mean computed in log space (no underflow); `maxFreq <= 0` → NaN guard skips AAs with no
table data; empty/zero-codon inputs return 0. `SplitIntoCodons` loop `i+2 < len` keeps all complete
codons (verified len 3 and 9), drops only the incomplete tail.

### Test quality audit
Assertions check exact sourced values (specific optimized sequences per Kazusa, hand-computed
CAI/GC, exact change tuples), independently re-translate to lock protein preservation, and cover
every Stage-A edge case. Real assertions, not "no-throw" tautologies. The one stochastic strategy
asserts only order-independent invariants.

### Findings / defects
None.

## Verdict & follow-ups
- Stage A PASS, Stage B PASS. **State: CLEAN.**
- No code changes. No logged defects.
- Non-defect notes: (1) `HarmonizeExpression` is intentionally stochastic; (2) `MinimizeSecondary`
  selects codons via BalancedOptimization (secondary-structure reduction lives in the separate
  `ReduceSecondaryStructure`); (3) CAI `1e-6` clamp applies only to incomplete custom tables, never
  to the built-in organisms. All three are documented and none affects the canonical model or any
  invariant.
