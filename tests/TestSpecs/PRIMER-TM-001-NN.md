# Test Specification: PRIMER-TM-001 (Nearest-Neighbour Salt-Corrected Tm)

**Test Unit ID:** PRIMER-TM-001
**Area:** MolTools
**Algorithm:** Nearest-neighbour (SantaLucia 1998) salt-corrected melting temperature (opt-in)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-24

> Scope: the **opt-in NN salt-corrected design Tm** added under PRIMER-TM-001
> (`CalculateMeltingTemperatureNN`, `CalculateNearestNeighborThermodynamics`). The Wallace /
> Marmur-Doty default Tm and the Primer3 penalty objective are specified in
> `PRIMER-TM-001-Penalty.md` and are unchanged.

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | SantaLucia (1998) PNAS 95:1460 | 1 | https://doi.org/10.1073/pnas.95.4.1460 | 2026-06-24 |
| 2 | SantaLucia & Hicks (2004) Annu Rev Biophys 33:415 (Table 1, Eq. 3, Eq. 5) | 1 | https://doi.org/10.1146/annurev.biophys.32.110601.141800 | 2026-06-24 |
| 3 | Owczarzy et al. (2004) Biochemistry 43:3537 (monovalent) | 1 | https://doi.org/10.1021/bi034621r | 2026-06-24 |
| 4 | Owczarzy et al. (2008) Biochemistry 47:5336 (divalent) | 1 | https://doi.org/10.1021/bi702363u | 2026-06-24 |
| 5 | Biopython Bio.SeqUtils.MeltingTemp (DNA_NN4, salt_correction 6/7) | 3 | https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/MeltingTemp.py | 2026-06-24 |

### 1.2 Key Evidence Points

1. Unified NN ΔH°/ΔS° at 1 M NaCl (Table 1) — source 2 (cross-checked vs source 5 DNA_NN4).
2. `Tm = ΔH°·1000/(ΔS° + R·ln(C_T/x)) − 273.15`, R = 1.9872, x = 4 non-self-comp / x = 1 self-comp — source 2 Eq. 3.
3. Published worked example: ΔH°=−43.5, ΔS°=−122.5, 0.2 mM → 35.8 °C — source 2 p.419.
4. Owczarzy 2004 monovalent: `1/Tm[Na]=1/Tm[1M]+(4.29·fGC−3.95)·1e-5·ln[Na]+9.40e-6·ln[Na]²` — sources 3, 5.
5. SantaLucia Eq. 5 entropy: `ΔS°[Na]=ΔS°+0.368·(N/2)·ln[Na]`, N=2·(L−1) — source 2 Eq. 5 (6-bp → N=10).
6. Owczarzy 2008 divalent Mg²⁺ 1/Tm correction with R=√[Mg]/[Mon] regimes — sources 4, 5.

### 1.3 Documented Corner Cases

- Self-complementary duplex: x=1, add symmetry ΔS°=−1.4, no terminal-AT for G·C ends (source 2).
- Terminal A·T penalty applied per A·T-closed end (source 2).
- Non-ACGT base has no NN parameter → not computable (source 5).

### 1.4 Known Failure Modes / Pitfalls

1. Using x=4 for a self-complementary duplex (or omitting the symmetry term) — wrong Tm — source 2.
2. Forgetting the terminal-AT penalty on A·T ends — wrong ΔH°/ΔS° — source 2.
3. Applying the salt correction as a flat additive °C instead of the 1/Tm (Kelvin) form — source 3/5.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateMeltingTemperatureNN(string, double, double, double, double, SaltCorrectionMode)` | `PrimerDesigner` | **Canonical** | NN Tm + salt corrections |
| `CalculateNearestNeighborThermodynamics(string)` | `PrimerDesigner` | **Canonical** | ΔH°/ΔS° + self-comp flag |
| `CalculateMeltingTemperature(string)` | `PrimerDesigner` | **Internal** | unchanged default; guarded by E2 |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | ΔH°/ΔS° = init + Σ NN stacks + terminal-AT(per A·T end) + symmetry(self-comp) | Yes | Source 2 Eq. 1, Table 1 |
| INV-2 | Self-complementary ⇒ x=1; otherwise x=4 in the C_T term | Yes | Source 2 Eq. 3 |
| INV-3 | Lower [Na⁺] ⇒ lower (Owczarzy-2004-corrected) Tm | Yes | Source 3 |
| INV-4 | Adding Mg²⁺ ⇒ higher Tm than the Mg²⁺-free buffer | Yes | Source 4 |
| INV-5 | Divalent mode with [Mg²⁺]=0 ≡ monovalent 2004 mode | Yes | Source 5 method 7 fallback |
| INV-6 | Default Tm method (Wallace/Marmur-Doty) unchanged | Yes | Scope requirement |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | NN thermo non-self-comp | ATGCATGC ΔH°/ΔS° + self-comp flag | ΔH°=−57.1, ΔS°=−156.5, selfComp=false | Source 2 Table 1 |
| M2 | NN thermo self-comp | GCGCGC ΔH°/ΔS° w/ symmetry, no terminal-AT | ΔH°=−50.4, ΔS°=−134.7, selfComp=true | Source 2 Table 1 |
| M3 | Tm equation worked example | ΔH°=−43.5, ΔS°=−122.5, 0.2 mM, x=4 | 35.8 °C (within 0.05) | Source 2 p.419 |
| M4 | NN Tm no salt, self-comp | GCGCGC, C_T=0.5 µM, mode None | 35.0473059911 °C | Source 2 Eq. 3 |
| M5 | NN Tm no salt, non-self-comp | ATGCATGC, x=4, mode None | 30.4338060665 °C | Source 2 Eq. 3 |
| M6 | Owczarzy 2004, self-comp | GCGCGC, 50 mM Na | 28.1593085080 °C | Sources 3,5 |
| M7 | Owczarzy 2004, non-self-comp | ATGCATGC, 50 mM Na | 18.1899960529 °C | Sources 3,5 |
| M8 | Default salt mode | omit saltMode ≡ Owczarzy2004Monovalent | equal | Source 3 (default) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | non-ACGT thermo | "ATGN" | null | NN lookup fails |
| S2 | too short thermo | "A" | null | no NN stack |
| S3 | SantaLucia Eq. 5 | GCGCGC, 50 mM Na, entropy mode | 24.9976652723 °C | N=10 |
| S4 | EcoRI self-comp 12-mer | CGCGAATTCGCG, mode None | 61.1452300219 °C | x=1 |
| S5 | [Na⁺] monotonicity | 1.0 M vs 0.01 M | tmLow < tmHigh | INV-3 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | divalent no-Mg fallback | Mg=0 divalent ≡ monovalent | equal | INV-5 |
| C2 | add Mg²⁺ | 3 mM Mg raises Tm | withMg > noMg | INV-4 |
| E1 | invalid input | "", null, "ATGN", "A" | NaN | no crash |
| E2 | default Tm unchanged | Wallace ATATATAT | 16.0 | INV-6 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `PrimerDesigner_MeltingTemperature_Tests.cs` — legacy Wallace/Marmur-Doty + flat salt; no NN Tm.
- `PrimerDesigner_Primer3Penalty_Tests.cs` — Primer3 penalty objective (separate sub-unit).
- `SequenceStatistics_CalculateThermodynamics_Tests.cs` — SEQ-THERMO-001, uses the 1997 (Allawi)
  parameters and no salt correction; different unit, not reused here.
- No pre-existing test exercises the SantaLucia-1998 NN Tm or the Owczarzy salt corrections.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M8, S1–S5, C1–C2, E1–E2 | ❌ Missing | New method; no prior coverage |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_NearestNeighborTm_Tests.cs` — all NN Tm cases.
- **Remove:** nothing (legacy fixtures cover the unchanged default Tm and the penalty objective).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| PrimerDesigner_NearestNeighborTm_Tests.cs | Canonical NN salt-corrected Tm | 17 |

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
| 9 | S1 | ❌ Missing | Implemented | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented | ✅ Done |
| 12 | S4 | ❌ Missing | Implemented | ✅ Done |
| 13 | S5 | ❌ Missing | Implemented | ✅ Done |
| 14 | C1 | ❌ Missing | Implemented | ✅ Done |
| 15 | C2 | ❌ Missing | Implemented | ✅ Done |
| 16 | E1 | ❌ Missing | Implemented | ✅ Done |
| 17 | E2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 17
**✅ Done:** 17 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | CalculateNearestNeighborThermodynamics_NonSelfComp_MatchesTable1Sum |
| M2 | ✅ | CalculateNearestNeighborThermodynamics_SelfComp_AppliesSymmetryNoTerminalAt |
| M3 | ✅ | TmEquation_PublishedWorkedExample_Gives35Point8 |
| M4 | ✅ | CalculateMeltingTemperatureNN_NoSalt_SelfComp_MatchesEquation |
| M5 | ✅ | CalculateMeltingTemperatureNN_NoSalt_NonSelfComp_UsesX4 |
| M6 | ✅ | CalculateMeltingTemperatureNN_Owczarzy2004_50mM_LowersTm |
| M7 | ✅ | CalculateMeltingTemperatureNN_Owczarzy2004_NonSelfComp_MatchesDerivation |
| M8 | ✅ | CalculateMeltingTemperatureNN_DefaultSaltMode_IsOwczarzy2004 |
| S1 | ✅ | CalculateNearestNeighborThermodynamics_NonAcgt_ReturnsNull |
| S2 | ✅ | CalculateNearestNeighborThermodynamics_SingleBase_ReturnsNull |
| S3 | ✅ | CalculateMeltingTemperatureNN_SantaLuciaEntropy_50mM_MatchesEq5 |
| S4 | ✅ | CalculateMeltingTemperatureNN_EcoRiSelfComp_NoSalt_MatchesEquation |
| S5 | ✅ | CalculateMeltingTemperatureNN_LowerSodium_LowersTm |
| C1 | ✅ | CalculateMeltingTemperatureNN_DivalentNoMg_EqualsMonovalent |
| C2 | ✅ | CalculateMeltingTemperatureNN_AddMagnesium_RaisesTm |
| E1 | ✅ | CalculateMeltingTemperatureNN_InvalidInput_ReturnsNaN |
| E2 | ✅ | CalculateMeltingTemperature_DefaultMethod_Unchanged |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Default C_T = 0.5 µM (PCR working conc.) | method default |
| 2 | Eq. 5 N = 2·(length−1) phosphates (confirmed by 6-bp→N=10 example) | SantaLuciaEntropy mode |

---

## 7. Open Questions / Decisions

1. The Owczarzy 2004 monovalent coefficients are taken from the Biopython reference implementation
   (method 6) since the Biochemistry 43:3537 full text is paywalled; independently corroborated by web
   search and the OligoPool tutorial equation form. No fabricated values.
2. **Decision:** the NN Tm is placed in `PrimerDesigner` (PRIMER-TM-001's class) as an opt-in method;
   the default `CalculateMeltingTemperature` is unchanged. SEQ-THERMO-001's `CalculateThermodynamics`
   uses the older 1997 parameters and is a distinct unit, not modified here.
