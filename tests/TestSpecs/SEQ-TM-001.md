# Test Specification: SEQ-TM-001

**Test Unit ID:** SEQ-TM-001
**Area:** Statistics
**Algorithm:** Melting Temperature (Wallace / Marmur-Doty GC formula, and nearest-neighbor Tm)
**Status:** ☑ Complete (resolved by consolidation — duplicate of SEQ-THERMO-001)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | SantaLucia, J. (1998). Unified NN thermodynamics. *PNAS* 95(4):1460–1465 | 1 | https://doi.org/10.1073/pnas.95.4.1460 | 2026-06-14 |
| 2 | Allawi & SantaLucia (1997). *Biochemistry* 36(34):10581–10594 (DNA_NN3, Table 1) | 1 | https://doi.org/10.1021/bi962590c | 2026-06-14 |
| 3 | Marmur & Doty (1962). *J Mol Biol* 5:109–118 (GC formula, via Biopython valueset 1) | 1 | (cited in source 5) | 2026-06-14 |
| 4 | Thein & Wallace (1986). Wallace rule (cited in Biopython `Tm_Wallace`) | 2 | (cited in source 5) | 2026-06-14 |
| 5 | Biopython `Bio.SeqUtils.MeltingTemp` (`Tm_Wallace`, `Tm_GC`, `Tm_NN`) | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py | 2026-06-14 |

### 1.2 Key Evidence Points

1. Wallace rule Tm = 4·(G+C) + 2·(A+T); `Tm_Wallace('ACGTTGCAATGCCGTA')` = 48.0 — source 5.
2. Marmur-Doty GC form Tm = 64.9 + 41·(GC − 16.4)/N (classic form; Biopython valueset 1 uses 69.3/650) — sources 3, 5.
3. NN Tm = (1000·ΔH)/(ΔS + R·ln(C_T/x)) − 273.15, R = 1.987, x = 4 (non-self-complementary) — sources 1, 5.
4. NN worked example `Tm_NN('CGTTCCAAAGATGTGGGCATGAGCTTAC')` = 60.32 at Na = 50 mM, C_T = 50 nM, method 5 — source 5.
5. Salt correction (method 5): ΔS(salt) = 0.368·(N−1)·ln[Na+] — sources 1, 5.

### 1.3 Documented Corner Cases

- Wallace rule documented only for 14–20 nt primers (source 5).
- NN model is undefined for length < 2 (no dinucleotide step) (source 5).
- NN Tm depends on C_T and self-complementarity factor x (sources 1, 5).

### 1.4 Known Failure Modes / Pitfalls

1. Using the wrong concentration factor (x = 1 self-complementary vs x = 4 non-self-complementary) shifts Tm — source 1.
2. Forgetting the +273.15 / −273.15 Kelvin conversion — source 5.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateThermodynamics(dnaSequence, naConc, primerConc)` | SequenceStatistics | **Canonical** | Nearest-neighbor ΔH/ΔS/ΔG/Tm — same method as SEQ-THERMO-001 |
| `CalculateMeltingTemperature(dnaSequence, useWallaceRule)` | SequenceStatistics | **Delegate** | Wallace / Marmur-Doty via `ThermoConstants` — same method as SEQ-THERMO-001 |

> These are the identical two methods already delivered, documented and tested under
> **SEQ-THERMO-001**. SEQ-TM-001 is a duplicate Registry entry — see §7.

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Wallace Tm = 2·(A+T) + 4·(G+C) for short oligos | Yes | Source 5 |
| INV-2 | Marmur-Doty Tm = 64.9 + 41·(GC − 16.4)/N | Yes | Sources 3, 5 |
| INV-3 | NN Tm = (1000·ΔH)/(ΔS + R·ln(C_T/4)) − 273.15, R = 1.987 | Yes | Sources 1, 5 |
| INV-4 | ΔG°37 = ΔH° − 310.15·ΔS°/1000 | Yes | Source 1 |
| INV-5 | Higher [Na+] raises Tm (salt term monotonic) | Yes | Source 1 |
| INV-6 | Empty / length-1 input → all-zero (NN undefined) | Yes | Source 5 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | NN reference value | `CGTTCCAAAGATGTGGGCATGAGCTTAC`, Na = 50 mM, C_T = 50 nM | Tm = 60.3 (60.32) | Source 5 |
| M2 | NN exact tuple | `GCGC` at repo defaults | ΔH −30.0, ΔS −84.91, ΔG −3.67, Tm −18.6 | Sources 1, 2 |
| M3 | Wallace short oligo | `ATGC`, useWallaceRule | 2·2 + 4·2 = 12.0 | Source 5 |
| M4 | Marmur-Doty GC formula | 20-mer, 10 GC | 64.9 + 41·(10 − 16.4)/20 = 51.78 | Sources 3, 5 |
| M5 | Empty / length-1 | `""`, `"A"` | all-zero | Source 5 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Salt monotonicity | `GCGC` at 0.05 M vs 1.0 M Na+ | Tm(1 M) > Tm(0.05 M) | INV-5 |
| S2 | Case-insensitivity | `gcgc` vs `GCGC` | equal | INV-* |
| S3 | Watson-Crick symmetry | `AAAA` vs `TTTT` | equal ΔH/ΔS | INV-* |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| (none) | | | | |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- The two canonical methods (`CalculateThermodynamics`, `CalculateMeltingTemperature`) are already implemented in `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs` and `src/.../Seqeron.Genomics.Infrastructure/ThermoConstants.cs`.
- The canonical test fixture `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateThermodynamics_Tests.cs` (delivered under SEQ-THERMO-001) already covers every MUST/SHOULD case above with exact, evidence-derived values (Wallace 48-family `ATGC`=12, Marmur-Doty 51.78, NN 60.3, GCGC tuple, salt monotonicity, symmetry, empty/length-1).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (NN reference 60.3) | ✅ Covered | Existing test `..._BiopythonWorkedExample_ReturnsTm60Point3` |
| M2 (GCGC tuple) | ✅ Covered | Existing `..._GcgcDefaults_ReturnsExactTuple` |
| M3 (Wallace `ATGC`=12) | ✅ Covered | Existing `..._ShortOligoWallace_...` |
| M4 (Marmur-Doty 51.78) | ✅ Covered | Existing `..._MarmurDoty_...` |
| M5 (empty / length-1) | ✅ Covered | Existing `..._EmptyInput_...`, `..._SingleBase_...` |
| S1 (salt monotonicity) | ✅ Covered | Existing `..._HigherSalt_RaisesMeltingTemperature` |
| S2 (case-insensitivity) | ✅ Covered | Existing `..._LowercaseInput_EqualsUppercase` |
| S3 (WC symmetry) | ✅ Covered | Existing `..._WatsonCrickSymmetry_...` |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateThermodynamics_Tests.cs` — already the single canonical fixture for both methods.
- **Remove:** nothing. Creating a second SEQ-TM-001 fixture would duplicate the existing tests, violating the duplicate-elimination rule. SEQ-TM-001 reuses the existing fixture.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceStatistics_CalculateThermodynamics_Tests.cs` | Canonical fixture (shared with SEQ-THERMO-001) | 13 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | (all) | ✅ Covered | None — already covered by canonical fixture | ✅ Done |

**Total items:** 0 missing/weak
**✅ Done:** 0 new | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M5, S1–S3 | ✅ | Covered by existing canonical fixture (verified passing this session) |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Repository NN default C_T = 250 nM differs from Biopython's 50 nM; the formula is identical, only the documented default parameter differs (passing `primerConcentration: 5e-8` reproduces 60.32). | Evidence §Assumptions; M1 uses the explicit 50 nM parameter |

---

## 7. Open Questions / Decisions

1. **DECISION — duplicate Registry entry.** The Processing Registry contains two entries for the identical pair of melting-temperature methods on the same class: **SEQ-THERMO-001** ("Thermodynamic Properties", ☑ Complete) and **SEQ-TM-001** ("Melting Temperature"). Both list `SequenceStatistics.CalculateThermodynamics(...)` as canonical and add `CalculateMeltingTemperature(...)`. SEQ-THERMO-001 already ships the implementation, the canonical test fixture, Evidence, and an algorithm doc. Per the prompt's duplicate-elimination rule ("one canonical test file per unit; no duplicate tests remain") and the workflow-control rule ("note the conflict in the TestSpec and update the checklist entry"), SEQ-TM-001 is **resolved by consolidation**: no new production code and no duplicate test file are created; this unit reuses the existing implementation and canonical fixture. Evidence was independently re-retrieved this session (Biopython `Tm_Wallace`/`Tm_GC`/`Tm_NN`, SantaLucia 1998, Allawi & SantaLucia 1997) to confirm the behavior is correct and source-backed. This mirrors the prior SEQ-COMPOSITION-001 ↔ SEQ-STATS-001 consolidation.
