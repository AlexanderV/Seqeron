# Test Specification: RNA-HAIRPIN-001

**Test Unit ID:** RNA-HAIRPIN-001
**Area:** RnaStructure
**Algorithm:** Hairpin Loop and Stem Free-Energy Calculation (Turner 2004)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Mathews et al. (2004), PNAS 101:7287-7292 (primary parameters) | 1 | https://doi.org/10.1073/pnas.0401799101 | 2026-06-14 |
| 2 | Turner & Mathews (2010), NAR 38:D280-D282 (NNDB) | 1 | https://doi.org/10.1093/nar/gkp892 | 2026-06-14 |
| 3 | NNDB Turner 2004 — hairpin rules/formula | 2 | https://rna.urmc.rochester.edu/NNDB/turner04/hairpin.html | 2026-06-14 |
| 4 | NNDB loop.txt (initiation), hairpin-mismatch-parameters.html (bonuses), tstack.txt (terminal mismatch), wc-parameters.html (stacking + AU end), triloop/tloop/hexaloop.txt (special), hairpin-example-1/-2.html (worked examples) | 2 | https://rna.urmc.rochester.edu/NNDB/turner04/ | 2026-06-14 |

<!-- NNDB pages retrieved from Wayback Machine snapshot 20240709061712 (wc-parameters via 20240709061712); live server was in maintenance on the access date. -->

### 1.2 Key Evidence Points

1. Hairpin loop ΔG°37 (n>3) = initiation(n) + terminal mismatch + UU/GA bonus + GG bonus + special-GU closure + all-C penalty — NNDB hairpin.html.
2. Hairpin loop ΔG°37 (n=3) = initiation(3) + all-C penalty (NO first-mismatch term) — NNDB hairpin.html.
3. Initiation(kcal/mol): 3→5.4, 4→5.6, 5→5.7, 6→5.4, 7→6.0, 8→5.5, 9→6.4 — NNDB loop.txt.
4. Bonuses/penalties: UU/GA −0.9, GG −0.8, special-GU −2.2, C3 +1.5, all-C linear A=0.3 B=1.6 — NNDB hairpin-mismatch-parameters.html.
5. Terminal mismatch (closing A-U, A·A) = −0.8; (closing A-U, G·G) = −0.8 — NNDB tstack.txt.
6. WC stacking CA/GU=−2.11, AC/UG=−2.24; per-AU-end penalty = +0.45 — NNDB wc-parameters.html.
7. Worked Example 1: loop(A-U,6,A..A)=+4.6, helix=−6.01, total −1.4; Example 2: loop(A-U,5,G..G)=+4.1 — NNDB hairpin-example-1/-2.html.
8. Special loops (total override): CAACG=6.8, CCUCGG=2.5, CAACGG=5.5, ACAGUGUU=1.8 — NNDB triloop/tloop/hexaloop.txt.

### 1.3 Documented Corner Cases

- Loops < 3 nt are prohibited by the nearest-neighbor rules (no defined energy) — hairpin.html.
- 3-nt loops receive no sequence-dependent first-mismatch term — hairpin.html.
- Special tri/tetra/hexaloops replace the model with an experimental total — hairpin.html, special tables.
- special-GU closure applies to a `G-U` closing pair only, not `U-G` — hairpin.html.
- Stem of P pairs has P−1 stacks; ≤1 pair ⇒ 0 stacking energy — wc-parameters.html.

### 1.4 Known Failure Modes / Pitfalls

1. Treating a `U-G` closing pair as eligible for the special-GU −2.2 bonus (it is not) — hairpin.html.
2. Applying a first-mismatch term to a 3-nt loop (it must not be applied) — hairpin.html.
3. Forgetting the per-AU-end penalty (+0.45) at a helix terminus that ends in A-U/U-A or G-U/U-G — wc-parameters.html.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateHairpinLoopEnergy(string loopSequence, char closingBase5, char closingBase3, bool specialGUClosure=false)` | RnaSecondaryStructure | **Canonical** | Turner 2004 additive hairpin model + special-loop override |
| `CalculateStemEnergy(string sequence, IReadOnlyList<BasePair> basePairs)` | RnaSecondaryStructure | **Canonical** | Turner 2004 nearest-neighbor stacking + AU/GU end penalty |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Hairpin loop energy is deterministic and depends only on (loop sequence, closing pair, special-GU flag). | Yes | NNDB hairpin.html (formula) |
| INV-2 | Loops < 3 nt return a prohibitive (non-finite-purpose) energy; never a normal low value. | Yes | NNDB hairpin.html (loops <3 nt prohibited) |
| INV-3 | A special tri/tetra/hexaloop returns exactly its tabulated total, independent of the additive terms. | Yes | NNDB triloop/tloop/hexaloop.txt |
| INV-4 | Stem energy of an empty base-pair list is 0; a stem of P pairs sums P−1 stacking terms. | Yes | NNDB wc-parameters.html |
| INV-5 | The special-GU −2.2 bonus is applied only when the closing pair is G(5')-U(3'), not U-G. | Yes | NNDB hairpin.html |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Hairpin Example 1 | `CalculateHairpinLoopEnergy("AAAAAA",'A','U')` (6-nt, closing A-U, first/last A) | 4.6 (= init 5.4 + tm −0.8) | hairpin-example-1.html, loop.txt, tstack.txt |
| M2 | Hairpin Example 2 (GG) | `CalculateHairpinLoopEnergy("GAAAG",'A','U')` (5-nt, first/last G) | 4.1 (= 5.7 + tm −0.8 + GG −0.8) | hairpin-example-2.html, hairpin-mismatch-parameters.html |
| M3 | Special triloop | `CalculateHairpinLoopEnergy("AAC",'C','G')` (key CAACG) | 6.8 (override) | triloop.txt |
| M4 | Special tetraloop | `CalculateHairpinLoopEnergy("CUCG",'C','G')` (key CCUCGG) | 2.5 (override) | tloop.txt |
| M5 | 3-nt loop, no mismatch term | `CalculateHairpinLoopEnergy("AAA",'G','C')` | 5.4 (= init(3), no first-mismatch) | hairpin.html, loop.txt |
| M6 | all-C 3-nt penalty | `CalculateHairpinLoopEnergy("CCC",'G','C')` | 6.9 (= init(3) 5.4 + C3 1.5) | hairpin-mismatch-parameters.html |
| M7 | Loop < 3 nt prohibited | `CalculateHairpinLoopEnergy("AA",'G','C')` | ≥ 100 (prohibitive) | hairpin.html (loops <3 prohibited) |
| M8 | Stem energy Example 1 helix | `CalculateStemEnergy` on pairs C-G,A-U,C-G,A-U | −6.01 (3 stacks + 1 AU end) | wc-parameters.html, hairpin-example-1.html |
| M9 | Empty stem | `CalculateStemEnergy(seq, [])` | 0 | wc-parameters.html (P−1 stacks) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | special-GU closure applied | `CalculateHairpinLoopEnergy("AAAA",'G','U', specialGUClosure:true)` minus same with flag false = −2.2 | difference −2.2 | bonus value −2.2 (hairpin-mismatch-parameters.html) |
| S2 | special-GU NOT for U-G | flag true but closing U-G ⇒ no −2.2 applied | same as flag false | INV-5, hairpin.html |
| S3 | UU/GA first mismatch bonus | `CalculateHairpinLoopEnergy("UAAU",'C','G')` includes −0.9 extra vs a non-UU/GA mismatch | −0.9 component present | hairpin-mismatch-parameters.html |
| S4 | all-C >3 nt linear penalty | `CalculateHairpinLoopEnergy("CCCC",'G','C')` = 5.6 + tm(GCCC −0.7) + (0.3·4+1.6) | 7.7 | hairpin-mismatch-parameters.html, tstack.txt |
| S5 | Stem AU/GU end penalty count | stem with both ends A-U gets +0.45 twice | two +0.45 | wc-parameters.html |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Hexaloop override | `CalculateHairpinLoopEnergy("CAGUGU",'A','U')` (key ACAGUGUU) | 1.8 | hexaloop.txt |
| C2 | Determinism (INV-1) | same inputs twice ⇒ identical output | equal | property |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file for `CalculateHairpinLoopEnergy` / `CalculateStemEnergy`. Sibling unit RNA-PAIR-001 has `RnaSecondaryStructure_CanPair_Tests.cs` but does not exercise the energy methods. This is a brand-new canonical fixture.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new fixture |
| M2 | ❌ Missing | new fixture |
| M3 | ❌ Missing | new fixture |
| M4 | ❌ Missing | new fixture |
| M5 | ❌ Missing | new fixture |
| M6 | ❌ Missing | new fixture |
| M7 | ❌ Missing | new fixture |
| M8 | ❌ Missing | new fixture |
| M9 | ❌ Missing | new fixture |
| S1 | ❌ Missing | new fixture |
| S2 | ❌ Missing | new fixture |
| S3 | ❌ Missing | new fixture |
| S4 | ❌ Missing | new fixture |
| S5 | ❌ Missing | new fixture |
| C1 | ❌ Missing | new fixture |
| C2 | ❌ Missing | new fixture |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_HairpinEnergy_Tests.cs` — all cases for both in-scope methods.
- **Remove:** nothing (no prior tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `RnaSecondaryStructure_HairpinEnergy_Tests.cs` | Canonical (RNA-HAIRPIN-001) | 16 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | S1 | ❌ Missing | Implemented | ✅ Done |
| 11 | S2 | ❌ Missing | Implemented | ✅ Done |
| 12 | S3 | ❌ Missing | Implemented | ✅ Done |
| 13 | S4 | ❌ Missing | Implemented | ✅ Done |
| 14 | S5 | ❌ Missing | Implemented | ✅ Done |
| 15 | C1 | ❌ Missing | Implemented | ✅ Done |
| 16 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M9 | ✅ | Evidence-based exact values from NNDB worked examples and tables |
| S1–S5 | ✅ | Bonus/penalty and end-penalty terms isolated |
| C1–C2 | ✅ | Hexaloop override + determinism property |

Total in-scope cases: 16; ✅ = 16.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Rounding to 2 decimals is a display choice; tests assert exact param sums within 1e-9 | all numeric MUST/SHOULD |

---

## 7. Open Questions / Decisions

1. Loops < 3 nt: the source states they are "prohibited" but assigns no number. The implementation returns a large prohibitive energy (100.0) so an optimizer never selects them; the test asserts ≥ 100 rather than an exact value because the source gives no value. Decision: accepted (INV-2).
