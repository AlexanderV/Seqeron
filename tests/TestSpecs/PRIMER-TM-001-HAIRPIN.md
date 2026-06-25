# Test Specification: PRIMER-TM-001-HAIRPIN

**Test Unit ID:** PRIMER-TM-001 (hairpin / secondary-structure Tm extension)
**Area:** MolTools
**Algorithm:** DNA hairpin (stem + loop) MFE folding + unimolecular hairpin Tm
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | SantaLucia & Hicks (2004) Annu Rev Biophys 33:415 — Table 1 NN, Table 4 hairpin loops, Eqs 7–11 | 1 | https://doi.org/10.1146/annurev.biophys.32.110601.141800 | 2026-06-25 |
| 2 | SantaLucia (1998) PNAS 95:1460 — unified NN set (reused stem stacks) | 1 | https://doi.org/10.1073/pnas.95.4.1460 | 2026-06-25 |
| 3 | Vallone & Benight (1999) — DNA hairpin Tm concentration-independence | 1 | https://pubmed.ncbi.nlm.nih.gov/10423551/ | 2026-06-25 |

### 1.2 Key Evidence Points

1. Hairpin ΔG°37 = Σ stem NN stacks (Table 1) + hairpin-loop initiation ΔG°37 by size (Table 4) — Source 1, p.428 Eq.10.
2. Hairpin loop ΔG°37 by size (kcal/mol): 3→3.5, 4→3.5, 5→3.3, 6→4.0, 7→4.2, 8→4.3, 9→4.5, 10→4.6, 12→5.0, 14→5.1, 16→5.3, 18→5.5, 20→5.7, 25→6.1, 30→6.3 — Source 1, Table 4.
3. Loop ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15 (destabilising loop, ΔG°37 > 0) — Source 1, Table 4 footnote a.
4. Unimolecular hairpin Tm = ΔH°·1000/ΔS° − 273.15, NO concentration term — Source 1, Eq.11; Source 3.
5. Loop < 3 nt sterically prohibited; minimum stem 2 bp (≥1 NN stack) — Source 1, "Hairpin Loops".
6. Jacobson-Stockmayer extrapolation for non-tabulated sizes: ΔG°37(n) = ΔG°37(x) + 2.44·R·310.15·ln(n/x) — Source 1, Eq.7.

### 1.3 Documented Corner Cases

- Loop < 3 nt → prohibited. Homopolymer / no complementary stem → no hairpin. Length 3/4 special bonuses are
  supplementary (not bundled; opt-in increment). Non-ACGT / empty / null → not computable.

### 1.4 Known Failure Modes / Pitfalls

1. Including the bimolecular duplex-initiation term in a unimolecular hairpin (wrong) — Source 1, Eq.10/11.
2. Adding a concentration term to the hairpin Tm (wrong; the transition is unimolecular) — Source 1 Eq.11; Source 3.
3. Loop ΔS° sign error (loop is destabilising; ΔS° must be negative) — Source 1, Table 4 footnote a.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindMostStableHairpin(string, int, double)` | `PrimerDesigner` | Canonical | MFE hairpin: stem span, length, loop size, ΔH°/ΔS°/ΔG°37 |
| `CalculateHairpinMeltingTemperature(string, int, double)` | `PrimerDesigner` | Canonical | Unimolecular Tm = ΔH°·1000/ΔS° − 273.15 |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | The returned hairpin minimises ΔG°37 over all stem/loop placements (MFE). | Yes | Source 1 (folding = minimum free energy) |
| INV-2 | A homopolymer / oligo with no WC stem closing a ≥3-nt loop returns null (no hairpin). | Yes | Source 1 (no stem possible) |
| INV-3 | Hairpin Tm is concentration-independent: Tm = ΔH°·1000/ΔS° − 273.15, no C_T term. | Yes | Source 1 Eq.11; Source 3 |
| INV-4 | Loop ΔH° contribution is 0; loop ΔS° = −ΔG°37·1000/310.15. | Yes | Source 1 Table 4 footnote a |
| INV-5 | Loop size < 3 is never returned (sterically prohibited). | Yes | Source 1 |
| INV-6 | Stem ΔH°/ΔS° equal the SantaLucia Table 1 NN stacks summed over the stem (no bimolecular init). | Yes | Source 1 Table 1 + Eq.10 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Canonical hairpin ΔH°/ΔS° | `FindMostStableHairpin("GGGCTTTTGCCC")` stem ΔH°/ΔS° + loop | ΔH° = −25.8; ΔS° = −75.48486216346927 | Table 1 + Table 4 (size 4) |
| M2 | Canonical hairpin ΔG°37 | same input | ΔG°37 = −2.3883700000000054 | derived |
| M3 | Canonical hairpin Tm | `CalculateHairpinMeltingTemperature("GGGCTTTTGCCC")` | 68.6403836682880 °C | Eq.11 |
| M4 | Folding correctness | same input: structure found | StemLength = 4, LoopSize = 4, StemStart = 0, StemEnd = 11 | MFE + Table 4 |
| M5 | Poly-A non-hairpin | `FindMostStableHairpin("AAAAAAAAAAAA")` | null | no WC stem |
| M6 | Poly-A Tm | `CalculateHairpinMeltingTemperature("AAAAAAAAAAAA")` | NaN | no hairpin |
| M7 | Loop ΔS° rule | loop-4 ΔS° increment in the result | total ΔS° − stem(−64.2) = −11.28486216346929 | Table 4 footnote |
| M8 | Concentration independence | Tm equals ΔH°·1000/ΔS° − 273.15 (no C_T) | matches M3 exactly | Eq.11; Source 3 |
| M9 | 5-nt loop uses size-5 increment | `FindMostStableHairpin("GGGCAAAAAGCCC")` | LoopSize 5; loop ΔG°37 = 3.3 → ΔG°37 = −2.5883700000000054 | Table 4 (size 5) |
| M10 | Loop < 3 never returned | `FindMostStableHairpin("GCGC")` (a 2-bp would close a 0-nt loop) | LoopSize ≥ 3 in any returned result, or null | Source 1 (loop≥3) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null / empty input | `FindMostStableHairpin(null/"")` | null; Tm NaN | invalid input |
| S2 | Non-ACGT input | `FindMostStableHairpin("GGGCNNNNGCCC")` | null | non-ACGT base |
| S3 | minStemLength < 2 guard | `FindMostStableHairpin("GGGCTTTTGCCC", 1)` | null | needs ≥1 stack |
| S4 | Jacobson-Stockmayer monotonicity | private loop ΔG°37 via a long-loop hairpin > tabulated size 30 | ΔG°37 increases with loop size | Eq.7 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Loop bonus increment | non-zero `loopBonusDeltaG37` shifts ΔG°37 by that amount | ΔG°37 increases by the bonus | opt-in supplementary term |
| C2 | Self-folding picks the most stable of competing stems | a sequence with two possible stems returns the lower-ΔG°37 one | MFE stem chosen | INV-1 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- New capability (`FindMostStableHairpin`, `CalculateHairpinMeltingTemperature`) — no prior tests existed.
  Existing `PrimerDesigner_NearestNeighborTm_Tests.cs` covers the duplex NN Tm only (unchanged).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M10, S1–S4, C1–C2 | ❌ Missing | brand-new methods |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_HairpinTm_Tests.cs` — all hairpin tests.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|-----------|
| `PrimerDesigner_HairpinTm_Tests.cs` | canonical hairpin folding + Tm | 16 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | M7 | ❌ Missing | implemented | ✅ Done |
| 8 | M8 | ❌ Missing | implemented | ✅ Done |
| 9 | M9 | ❌ Missing | implemented | ✅ Done |
| 10 | M10 | ❌ Missing | implemented | ✅ Done |
| 11 | S1 | ❌ Missing | implemented | ✅ Done |
| 12 | S2 | ❌ Missing | implemented | ✅ Done |
| 13 | S3 | ❌ Missing | implemented | ✅ Done |
| 14 | S4 | ❌ Missing | implemented | ✅ Done |
| 15 | C1 | ❌ Missing | implemented | ✅ Done |
| 16 | C2 | ❌ Missing | implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact ΔH°/ΔS° asserted |
| M2 | ✅ | exact ΔG°37 asserted |
| M3 | ✅ | exact Tm asserted |
| M4 | ✅ | stem/loop span asserted |
| M5 | ✅ | null asserted |
| M6 | ✅ | NaN asserted |
| M7 | ✅ | loop ΔS° increment asserted |
| M8 | ✅ | Tm = ΔH°·1000/ΔS° − 273.15 asserted |
| M9 | ✅ | size-5 loop ΔG°37 asserted |
| M10 | ✅ | loop ≥ 3 / null asserted |
| S1 | ✅ | null/empty/NaN asserted |
| S2 | ✅ | non-ACGT null asserted |
| S3 | ✅ | minStemLength<2 null asserted |
| S4 | ✅ | J-S monotonicity asserted |
| C1 | ✅ | loop bonus shift asserted |
| C2 | ✅ | MFE stem selection asserted |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Bimolecular duplex-initiation term excluded from the unimolecular hairpin (loop init is the nucleation cost). | INV-6, M1 |
| 2 | Terminal-AT penalty not applied at the open stem end of the hairpin core (article body Eqs 8–10 add only stem stacks + loop). | M1, M2 |

---

## 7. Open Questions / Decisions

1. The supplementary triloop/tetraloop bonus and terminal-mismatch tables (length-3/4 loops) are not bundled;
   exposed as an opt-in `loopBonusDeltaG37` increment (default 0). This is the honest residual recorded in
   LIMITATIONS / the report — the stem-stack + loop-initiation core is exact and fully sourced.
