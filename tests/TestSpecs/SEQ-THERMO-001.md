# Test Specification: SEQ-THERMO-001

**Test Unit ID:** SEQ-THERMO-001
**Area:** Statistics
**Algorithm:** DNA Duplex Thermodynamics (Nearest-Neighbor ΔH°/ΔS°/ΔG°/Tm)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Allawi & SantaLucia (1997), Biochemistry 36(34):10581-10594 | 1 | https://doi.org/10.1021/bi962590c | 2026-06-13 (corroborated via DNA_NN3) |
| 2 | SantaLucia (1998), PNAS 95(4):1460-1465 | 1 | https://doi.org/10.1073/pnas.95.4.1460 | 2026-06-13 (corroborated via MELTING/Wikipedia) |
| 3 | Biopython Bio.SeqUtils.MeltingTemp (DNA_NN3, Tm_NN) | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py | 2026-06-13 |
| 4 | MELTING 5 User Guide (Dumousseau et al. 2012) | 3 | https://www.ebi.ac.uk/biomodels/tools/melting/melting5-UserGuide.pdf | 2026-06-13 |
| 5 | Wikipedia — Nucleic acid thermodynamics | 4 | https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics | 2026-06-13 |

### 1.2 Key Evidence Points

1. NN table (kcal/mol ΔH, cal/(mol·K) ΔS): AA/TT −7.9/−22.2; AT −7.2/−20.4; TA −7.2/−21.3; CA/TG −8.5/−22.7; GT/AC −8.4/−22.4; CT/AG −7.8/−21.0; GA/TC −8.2/−22.2; CG −10.6/−27.2; GC −9.8/−24.4; GG/CC −8.0/−19.9. — Biopython DNA_NN3 (Allawi & SantaLucia 1997).
2. Initiation applied at BOTH termini: init_A/T (2.3, 4.1), init_G/C (0.1, −2.8). — Biopython Tm_NN (`ends = seq[0]+seq[-1]`); corroborated by Wikipedia terminal A/T and G/C entries.
3. Tm = (1000·ΔH)/(ΔS + R·ln(C_T/F)) − 273.15; R = 1.987; F = 4 default (non-self-complementary, equimolar). — Biopython Tm_NN; MELTING §4.2/§4.3.
4. Salt correction (method 5): ΔS += 0.368·(N−1)·ln[Na+], [Na+] in mol/L. — Biopython salt_correction; SantaLucia (1998).
5. Worked example Tm_NN('CGTTCCAAAGATGTGGGCATGAGCTTAC') = 60.32 °C at dnac1=dnac2=25 nM, Na=50 mM. — Biopython docstring.

### 1.3 Documented Corner Cases

- Length < 2 has no dinucleotide step ⇒ NN model undefined (Biopython). Repository contract: return `(0,0,0,0)`.
- Sequences processed case-insensitively (Biopython upper-cases input).
- F = 1 for self-complementary, F = 4 for non-self-complementary equimolar strands (MELTING §4.3); this unit uses F = 4 (default).

### 1.4 Known Failure Modes / Pitfalls

1. Applying initiation to only one terminus (first base) under-counts by one init term — defect corrected in this unit. — Biopython Tm_NN two-end logic.
2. Forgetting the kcal→cal conversion (×1000 on ΔH) in the Tm equation. — Biopython Tm_NN.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateThermodynamics(string, double, double)` | SequenceStatistics | Canonical | NN ΔH°/ΔS°/ΔG°/Tm; deep evidence-based testing. |
| `CalculateMeltingTemperature(string, bool)` | SequenceStatistics | Delegate | Simple Wallace / Marmur-Doty Tm via ThermoConstants; smoke verification. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | Initiation contributes once at each duplex terminus (2 init terms). | Yes | Biopython Tm_NN (`ends = seq[0]+seq[-1]`); Allawi & SantaLucia (1997) Table 1 |
| INV-02 | ΔG°₃₇ = ΔH° − 310.15·ΔS°/1000 (kcal/mol). | Yes | Gibbs relation; SantaLucia (1998) |
| INV-03 | Tm = (1000·ΔH°)/(ΔS° + R·ln(C_T/4)) − 273.15, R = 1.987. | Yes | Biopython Tm_NN; MELTING §4.2 |
| INV-04 | The NN table is Watson-Crick symmetric (AA=TT, CA=TG, GT=AC, CT=AG, GA=TC, GG=CC). | Yes | DNA_NN3 (Allawi & SantaLucia 1997) |
| INV-05 | Result is deterministic and case-insensitive. | Yes | Biopython upper-cases input |
| INV-06 | Empty or length-1 input returns (0,0,0,0). | Yes | NN model undefined for length < 2 (ASSUMPTION on return shape) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Biopython worked example | Tm of CGTTCCAAAGATGTGGGCATGAGCTTAC, Na=0.05, C_T=50 nM | Tm = 60.3 °C | Biopython Tm_NN docstring (60.32) |
| M2 | GCGC full result | Repo defaults Na=0.05, C_T=250 nM | ΔH=−30.0, ΔS=−84.91, ΔG=−3.67, Tm=−18.6 | Derivation from DNA_NN3 + SantaLucia (1998) |
| M3 | Two-end initiation (mixed ends) | ATCG (A…G) ΔH/ΔS reflect init_A/T + init_G/C | ΔH=−23.6, ΔS=−71.81 | Biopython two-end init |
| M4 | Both-AT ends | AATT init_A/T applied twice | ΔH=−18.4, ΔS=−59.91, ΔG=0.18, Tm=−75.0 | DNA_NN3 init_A/T; derivation |
| M5 | Empty input | "" returns all zero | (0,0,0,0) | NN undefined length<2 |
| M6 | Length-1 input | "A" returns all zero | (0,0,0,0) | NN undefined length<2 |
| M7 | ΔG relation (INV-02) | ΔG = ΔH − 310.15·ΔS/1000 holds for GCGC | matches rounded −3.67 | Gibbs relation |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitivity (INV-05) | lowercase "gcgc" equals "GCGC" | identical result | implementation upper-cases |
| S2 | Higher salt raises Tm | Na 1.0 M vs 0.05 M for same seq | Tm(1.0) > Tm(0.05) | salt entropy term monotonic |
| S3 | NN symmetry (INV-04) | AA-run vs TT-run give equal ΔH/ΔS | equal | DNA_NN3 symmetric |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Delegate Wallace Tm | CalculateMeltingTemperature short oligo | 2(A+T)+4(G+C) | smoke for delegate |
| C2 | Delegate Marmur-Doty | long sequence GC formula | 64.9 + 41(GC−16.4)/N | smoke for delegate |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior canonical test file for `CalculateThermodynamics` existed; searched `tests/Seqeron/Seqeron.Genomics.Tests/` (no `*Thermo*` file). `CalculateMeltingTemperature` had no dedicated unit either.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |
| C2 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateThermodynamics_Tests.cs` — all cases above.
- **Remove:** none (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| SequenceStatistics_CalculateThermodynamics_Tests.cs | Canonical (CalculateThermodynamics) + delegate smoke (CalculateMeltingTemperature) | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented worked-example Tm test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented GCGC full-result test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented mixed-end init test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented AATT test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented empty-input test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented length-1 test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented ΔG-relation test | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented case-insensitivity test | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented salt-monotonicity test | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented NN-symmetry test | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented Wallace delegate smoke | ✅ Done |
| 12 | C2 | ❌ Missing | Implemented Marmur-Doty delegate smoke | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact Tm 60.3 |
| M2 | ✅ Covered | exact 4-tuple |
| M3 | ✅ Covered | exact ΔH/ΔS |
| M4 | ✅ Covered | exact 4-tuple |
| M5 | ✅ Covered | (0,0,0,0) |
| M6 | ✅ Covered | (0,0,0,0) |
| M7 | ✅ Covered | Gibbs relation |
| S1 | ✅ Covered | case-insensitive |
| S2 | ✅ Covered | monotonic Tm |
| S3 | ✅ Covered | NN symmetry |
| C1 | ✅ Covered | Wallace smoke |
| C2 | ✅ Covered | Marmur-Doty smoke |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Empty/length-1 input returns (0,0,0,0) (API edge convention; no thermodynamic value affected). | M5, M6, INV-06 |

---

## 7. Open Questions / Decisions

1. Decision: this unit covers the default non-self-complementary case (F = 4). Self-complementary Tm (F = 1) and divalent-cation (Mg²⁺) corrections are out of scope and noted in the algorithm doc §5.3.
