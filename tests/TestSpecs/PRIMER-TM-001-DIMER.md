# Test Specification: PRIMER-TM-001-DIMER

**Test Unit ID:** PRIMER-TM-001 (self-/hetero-dimer Tm extension)
**Area:** MolTools
**Algorithm:** Self-dimer / hetero-dimer (intermolecular) Tm via thermodynamic alignment
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | SantaLucia J, Hicks D (2004), Annu Rev Biophys 33:415-440 | 1 | https://doi.org/10.1146/annurev.biophys.32.110601.141800 | 2026-06-25 |
| 2 | Untergasser A et al. (2012), Nucleic Acids Res 40:e115 (Primer3) | 1 | https://doi.org/10.1093/nar/gks596 | 2026-06-25 |
| 3 | Primer3 `thal.c` (ntthal), primer3-py vendored libprimer3 | 3 | https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/thal.c | 2026-06-25 |
| 4 | primer3-py 2.3.0 (`calc_homodimer`, `calc_heterodimer`) | 3 | https://pypi.org/project/primer3-py/ | 2026-06-25 |

### 1.2 Key Evidence Points

1. Dimer ΔH°/ΔS° = init(+0.2, −5.7) + Σ SantaLucia unified NN stacks + terminal-A·T penalty(+2.2, +6.9) per A·T-closed end — SantaLucia & Hicks (2004) Table 1; thal.c lines 128-129, 588-589.
2. Entropy salt correction baked into ΔS°: `ΔS° += 0.368·N_stacks·ln[Na⁺]` — SantaLucia & Hicks (2004) Eq. 5; thal.c lines 623-624, 1042.
3. Bimolecular Tm = ΔH°·1000/(ΔS° + R·ln(C_T/x)) − 273.15, R = 1.9872; x = 1 iff both oligos are reverse-complement palindromes, else x = 4 — SantaLucia & Hicks (2004) Eq. 3; thal.c lines 590-593, 2771.
4. ntthal default dimer concentration C_T = dna_conc = 50 nM — thal.c lines 829/844.
5. The engine returns the most stable (highest-Tm) structure — Untergasser et al. (2012).

### 1.3 Documented Corner Cases

- Poly-A self-dimer → `structure_found = False` (no stable dimer); method returns null/NaN.
- Self-dimer of a non-palindromic oligo uses x = 4 (it is not "symmetric" in ntthal's sense).

### 1.4 Known Failure Modes / Pitfalls

1. ntthal extends some duplexes into terminal mismatches/overhangs (`tstack2`); the plain contiguous-WC NN model deviates by a small terminal term for those (e.g. poly-A overhangs, ATCGTTAC/GTAACGAT) — documented model boundary. Parity is asserted only on cases whose optimal structure is a contiguous WC duplex.
2. Confusing the dimer C_T convention (50 nM) with the primer-Tm convention (0.5 µM) shifts Tm by ~6 °C.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindMostStableDimer(string, string, double, double)` | PrimerDesigner | **Canonical** | Thermodynamic alignment; returns `DimerResult?` (ΔH°/ΔS°/ΔG°37 + spans/bp). |
| `CalculateDimerMeltingTemperature(string, string, double, double)` | PrimerDesigner | **Canonical** | Bimolecular Tm of the most stable dimer. |
| `CalculateSelfDimerMeltingTemperature(string, double, double)` | PrimerDesigner | **Delegate** | Calls the two-argument method with the sequence twice. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | ΔH° = init + Σ stacks + A·T penalty (salt-independent); equals the hand-derived value | Yes | SantaLucia & Hicks (2004) Table 1 |
| INV-2 | Tm matches primer3/ntthal for any pair whose optimum is a contiguous WC duplex | Yes | primer3-py 2.3.0 |
| INV-3 | x = 1 iff both oligos are reverse-complement palindromes, else x = 4 | Yes | thal.c 590-593, 2771 |
| INV-4 | No Watson-Crick duplex of ≥ 2 bp ⇒ null / NaN | Yes | primer3 `structure_found=False` |
| INV-5 | Lower [Na⁺] strictly lowers Tm for a fixed duplex (Eq. 5 monotonicity) | Yes | SantaLucia & Hicks (2004) Eq. 5 |
| INV-6 | Self-dimer convenience method ≡ two-argument method with same sequence twice | Yes | delegation |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | GCGCGCGC self-dimer thermo | ΔH°/ΔS°/ΔG°37 of the full 8-bp self-dimer (x=1) | ΔH°=−70.8 kcal/mol; ΔS°=−192.61700633667505; ΔG°37=−11.059835484680235 | SantaLucia & Hicks (2004) Table 1 (hand-derived); primer3-py |
| M2 | GCGCGCGC self-dimer Tm | Bimolecular Tm at C_T=50 nM, x=1 | 40.09064476882935 °C | hand-derived; primer3-py 40.0906 |
| M3 | TGCATGCATG/CATGCATGCA hetero-dimer | Non-palindromic pair (x=4) ΔH°/ΔS°/Tm | ΔH°=−74.1 kcal/mol; ΔS°=−211.8218652900108; Tm=25.659587124835923 °C | hand-derived; primer3-py −74100/−211.8219/25.6596 |
| M4 | ATCGATCGATCG/CGATCGATCGAT Tm | primer3 parity (full 12-bp duplex, palindromic pair, x=1) | 32.6107 °C (±1e-3) | primer3-py 2.3.0 |
| M5 | CGATCGATCG self-dimer Tm | palindromic self-dimer (x=1) parity | 29.6600 °C (±1e-3) | primer3-py 2.3.0 |
| M6 | GGGGCCCC dimer Tm | parity | 29.0150 °C (±1e-3) | primer3-py 2.3.0 |
| M7 | GCATGC self-dimer Tm | parity near 0 °C | 0.6859 °C (±1e-3) | primer3-py 2.3.0 |
| M8 | ACGTACGTACGT self-dimer Tm | parity | 37.6251 °C (±1e-3) | primer3-py 2.3.0 |
| M9 | Poly-A self-dimer | no stable dimer | `FindMostStableDimer` null; Tm NaN | primer3 `structure_found=False` |
| M10 | Non-complementary pair (GGGGGGGG/AAAAAAAA) | no duplex | null / NaN | INV-4 |
| M11 | Invalid input | null, < 2 bases, non-ACGT | null / NaN | input contract |
| M12 | Base-pair count / spans | GCGCGCGC self-dimer reports 8 bp at start 0 | BasePairs=8, Strand1Start=0, Strand2Start=0 | alignment definition |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Self-dimer delegation | `CalculateSelfDimerMeltingTemperature(s)` ≡ `CalculateDimerMeltingTemperature(s,s)` | equal | INV-6 |
| S2 | Salt monotonicity | lower [Na⁺] ⇒ lower Tm for GCGCGCGC | Tm(10 mM) < Tm(1 M) | INV-5 |
| S3 | Concentration effect | higher C_T ⇒ higher bimolecular Tm | Tm(500 nM) > Tm(50 nM) | Eq. 3 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Most-stable selection | a longer/stronger duplex is chosen over a short one | longer-run result returned | ntthal max-Tm objective |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests cover self-/hetero-dimer Tm. Sibling NN tests live in `PrimerDesigner_NearestNeighborTm_Tests.cs` and `PrimerDesigner_HairpinTm_Tests.cs`; the legacy `HasPrimerDimer` (boolean complementarity heuristic) is unrelated to this thermodynamic dimer Tm.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M12, S1–S3, C1 | ❌ Missing | New unit; no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_DimerTm_Tests.cs` — all dimer-Tm tests.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `PrimerDesigner_DimerTm_Tests.cs` | Canonical dimer-Tm fixture | 19 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented (hand-derived ΔH/ΔS/ΔG37) | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented (hand-derived Tm) | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented (x=4 hetero-dimer) | ✅ Done |
| 4 | M4–M8 | ❌ Missing | Implemented (primer3 parity) | ✅ Done |
| 5 | M9 | ❌ Missing | Implemented (poly-A null/NaN) | ✅ Done |
| 6 | M10 | ❌ Missing | Implemented (non-complementary) | ✅ Done |
| 7 | M11 | ❌ Missing | Implemented (invalid input) | ✅ Done |
| 8 | M12 | ❌ Missing | Implemented (bp count / spans) | ✅ Done |
| 9 | S1–S3 | ❌ Missing | Implemented | ✅ Done |
| 10 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 10 groups (19 tests)
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M12 | ✅ Covered | Exact values, `.Within(1e-9)` (hand-derived) / `.Within(1e-3)` (primer3 parity) |
| S1–S3 | ✅ Covered | Delegation + monotonicity |
| C1 | ✅ Covered | Most-stable selection |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Gapless alignment only (no internal loops / terminal-overhang extension); parity asserted only on contiguous-WC-optimum cases | §1.4, §4.1 (M4–M8 tolerance) |

---

## 7. Open Questions / Decisions

1. None. Internal loops / terminal-overhang extension (the ntthal `tstack2` terms) are intentionally out of scope and documented as a model boundary in the algorithm doc §5.3.
