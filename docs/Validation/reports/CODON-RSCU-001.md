# Validation Report: CODON-RSCU-001 — Relative Synonymous Codon Usage (RSCU) and codon counting

- **Validated:** 2026-06-15   **Area:** Codon
- **Canonical method(s):** `CodonUsageAnalyzer.CalculateRscu(DnaSequence)` / `(string)`; `CodonUsageAnalyzer.CountCodons(DnaSequence)` / `(string)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm (retrieved this session, 2026-06-15)

| Source | Retrieved | Confirms |
|--------|-----------|----------|
| LIRMM "RSCU RS" — https://www.lirmm.fr/~rivals/rscu/ (WebFetch) | verbatim formula | `RSCU_{i,j} = (n_i × x_{i,j}) / Σ_{j=1..n_i} x_{i,j}`, with `n_i` = number of codons coding amino acid i (degeneracy), `x_{i,j}` = occurrences of codon j. No-bias value = 1. |
| GenomicSig (CRAN) RSCU — https://rdrr.io/cran/GenomicSig/man/RSCU.html (WebFetch) + WebSearch synthesis | definition + bounds | "ratio of observed to the number observed under uniform synonymous usage"; RSCU is "a real value comprised between 0 and the number of synonymous codons for that amino acid"; no-bias = 1.00; >1 overused, <1 underused. |
| seqinr `uco` (via WebSearch over CRAN refman) | definition + primary ref | "observed relative to the number of times the codon would be observed for a uniform synonymous codon usage"; "In the absence of any codon usage bias, the RSCU values would be 1.00"; primary ref Sharp, Tuohy & Mosurski (1986). |

The formula obtained verbatim from LIRMM this session is **algebraically identical** to the form in
the implementation's docstring (`RSCU = (n·x_j)/Σx`) and to the TestSpec/Evidence (Source 2/3). The
degeneracy normalisation, the "observed/expected-under-uniform" interpretation, the bounds `[0, n_i]`,
and the no-bias value 1 are all confirmed by independently-retrieved sources.

### Formula check
`RSCU_{i,j} = (n_i · x_{i,j}) / Σ_{k=1}^{n_i} x_{i,k}` — matches LIRMM verbatim and GenomicSig.
Equivalent to `x_j / ((1/n)·Σx)`, i.e. observed / expected-under-uniform. Correct.

### Edge-case semantics check
- **Single-codon families (Met, Trp):** n=1 ⇒ RSCU=1 whenever present (degenerate case of the formula). Sourced (GenomicSig bounds; seqinr). Correct.
- **Absent family (0/0):** canonical RSCU undefined. cubar uses a pseudocount (default 1); repository convention is to return 0. This is documented (TestSpec Assumption #1) and only affects families that never occur — never the RSCU of an observed codon. Acceptable, documented convention (not a divergence for present families).
- **Stop codons:** repository groups TAA/TAG/TGA as one 3-fold synonymous family and applies the same family-ratio formula. Reference tools commonly *exclude* stops, but this choice cannot change the RSCU of any amino-acid codon. Documented (Assumption #2). Acceptable convention.

### Independent cross-check (hand computation against the sourced formula)
- **Leu (n=6), `CTGCTGCTGCTA`**, CTG×3 CTA×1 total=4: RSCU(CTG)=6·3/4=**4.5**, RSCU(CTA)=6·1/4=**1.5**, others=0; Σ=6=n_i. ✓
- **Phe (n=2), `TTTTTTTTC`**, TTT×2 TTC×1 total=3: RSCU(TTT)=2·2/3=**4/3**, RSCU(TTC)=2·1/3=**2/3**. ✓
- **Phe equal, `TTTTTC`**: RSCU=2·1/2=**1.0** each. ✓
- **Met (n=1), `ATGATG`**: RSCU(ATG)=1·2/2=**1.0**. ✓
- **Stop family `TAATAGTGA`** (each ×1, n=3): RSCU=3·1/3=**1.0** each. ✓
- **Stop family biased `TAATAATGA`** (TAA×2, TGA×1): RSCU(TAA)=3·2/3=**2.0**, RSCU(TGA)=**1.0**, RSCU(TAG)=**0.0**; Σ=3=n_i. ✓

All values trace to the sourced formula (LIRMM/GenomicSig), not to code output.

### Findings / divergences
None material. Two documented conventions (absent-family→0; stop-codon family) are sound and do
not affect amino-acid RSCU. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs`:
- `CalculateRscuCore` (lines 88–110): groups `CodonToAminoAcid` by amino acid, computes
  `expected = totalCount/numSynonymous`, then `rscu[codon] = expected>0 ? observed/expected : 0`.
  This is exactly `n_i·x_j/Σx` (since `observed/(Σx/n) = n·observed/Σx`) with the absent-family
  guard returning 0. Realises the validated formula faithfully.
- `CountCodonsCore` (lines 37–52): non-overlapping triplets from offset 0 (`i += 3`), trailing
  1–2 bases ignored (`i+3 <= length`), `IsValidCodon` excludes any non-ACGT triplet. Matches spec.
- String overloads (lines 29–35, 80–86): null/empty → empty dict; `ToUpperInvariant()` before
  counting (case-insensitive). DnaSequence overloads `ArgumentNullException.ThrowIfNull`.

### Formula realised correctly?
Yes. The group-by-amino-acid + per-family ratio computes the sourced formula exactly. The
single-codon degenerate case falls out naturally (n=1 ⇒ observed/observed=1). Stop codons form a
3-fold family per the table. The absent-family branch returns 0 (documented convention).

### Cross-verification table recomputed vs code
All six hand-computed cases above were verified by adding tests and running them; the code returns
the externally-derived values exactly (Within 1e-10): 4.5/1.5/0; 4/3, 2/3; 1.0/1.0; 1.0; absent→0;
stop 1.0/1.0/1.0 and biased 2.0/1.0/0.0.

### Variant/delegate consistency
String and DnaSequence overloads delegate to the same `*Core` (S5/S6 confirm exact agreement;
case-insensitivity confirmed).

### Test quality audit (test-quality gate)
- **Sourced expectations, not code echoes:** every MUST/COULD value (4.5, 1.5, 4/3, 2/3, 1.0,
  family sum=6, bounds [0,6]) traces to the LIRMM/GenomicSig formula independently retrieved this
  session, not to code output. A deliberately-wrong implementation (e.g. per-thousand normalisation,
  or treating each codon as its own family) would fail M1/M2/M5.
- **No green-washing:** assertions use exact values with `Within(1e-10)`; no weakened
  Greater/Contains/range substitutes for known exacts; no skips; no widened tolerances.
- **Cover all the logic — GAP FOUND AND FIXED:** the original 16-test fixture covered all formula
  paths, both invariants, all CountCodons paths and all input guards, BUT did **not** exercise two
  documented Stage-A branches: (1) the **absent-family `expected>0 ? … : 0`** code branch
  (Assumption #1), and (2) the **stop-codon 3-fold family** handling (Assumption #2). These are real
  code branches with defined, sourced expected values that no test locked — a coverage gap per the
  gate. Added three tests (sourced values, hand-computed independently):
  - `CalculateRscu_AbsentFamily_ReturnsZeroForEveryCodon` (`ATGATG` ⇒ all Leu/Phe = 0)
  - `CalculateRscu_StopCodonFamily_TreatedAsThreeFoldFamily` (`TAATAGTGA` ⇒ 1.0 each)
  - `CalculateRscu_StopFamilyBiased_ComputesFamilyRatioAndSumsToDegeneracy` (`TAATAATGA` ⇒ 2.0/1.0/0.0, Σ=3)
- **Honest green:** full unfiltered suite **6526 passed, 0 failed** (1 pre-existing skipped
  benchmark, unrelated); `dotnet build` 0 errors; the one changed file builds warning-free (the 4
  build warnings are in unrelated `ApproximateMatcher_EditDistance_Tests.cs`).

### Findings / defects
No algorithm defect. One **test-coverage gap** (two documented branches untested) — **fixed in
session** by adding three exact-value tests. Fixture 16 → 19 tests.

## Verdict & follow-ups
- **Stage A: PASS** — formula, bounds, no-bias value and conventions confirmed against
  independently-retrieved LIRMM + GenomicSig + seqinr.
- **Stage B: PASS-WITH-NOTES** — code faithfully realises the formula; minor test-coverage gap
  (absent-family + stop-codon-family branches) fixed this session.
- **Test-quality gate: PASS** (after fix).
- **End state: ✅ CLEAN** — coverage gap completely closed; full suite green.
