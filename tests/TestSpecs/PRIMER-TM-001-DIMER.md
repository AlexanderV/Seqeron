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

1. The full ntthal DP (`CalculateDimerThermodynamicsNtthal`, to which the Tm methods delegate) models internal mismatches/loops, bulges and terminal overhangs (`tstack2`/`dangle`), so parity holds for non-contiguous optima too (N1–N5). The legacy `FindMostStableDimer` `DimerResult` scorer remains contiguous-WC only (its `BasePairs`/spans/ΔH°/ΔS° fields); M12 asserts the contiguous case there.
2. Confusing the dimer C_T convention (50 nM) with the primer-Tm convention (0.5 µM) shifts Tm by ~6 °C.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindMostStableDimer(string, string, double, double)` | PrimerDesigner | **Canonical** | Contiguous-WC thermodynamic alignment; returns `DimerResult?` (ΔH°/ΔS°/ΔG°37 + spans/bp). |
| `CalculateDimerThermodynamicsNtthal(string, string, double, double)` | PrimerDesigner | **Canonical** | Full ntthal DP (internal mismatches/loops, bulges, terminal overhangs); returns `DimerThermodynamics?`. |
| `CalculateDimerMeltingTemperature(string, string, double, double)` | PrimerDesigner | **Canonical** | Bimolecular Tm of the most stable dimer (delegates to the full ntthal DP). |
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
| N1 | Full DP — 2×2 internal loop | GCGCATGCGC self-dimer (`CalculateDimerThermodynamicsNtthal`) | Tm=43.1572 °C; ΔH=−84.400 kcal/mol; ΔS=−233.42187; ΔG°37=−12.00421 kcal/mol (±1e-3) | primer3-py 2.3.0 `calc_homodimer` |
| N2 | Full DP — 3×3 internal loop | GCGCAAAGCGC/GCGCTTTGCGC | Tm=41.8816 °C; ΔH=−92.300; ΔS=−256.82429; ΔG°37=−12.64594 (±1e-3) | primer3-py 2.3.0 `calc_heterodimer` |
| N3 | Full DP — single-base bulge | GCGCGCGC/GCGCAGCGC | Tm=19.8125 °C; ΔH=−70.800; ΔS=−205.50701; ΔG°37=−7.06200 (±1e-3) | primer3-py 2.3.0 |
| N4 | Full DP — mixed 2×2 internal loop | GCGCACGCGC/GCGCTAGCGC | Tm=18.5604 °C; ΔH=−68.400; ΔG°37=−6.89198 (±1e-3) | primer3-py 2.3.0 |
| N5 | Full DP — terminal overhang | GCGCGCAAAA/AAAAGCGCGC | Tm=24.6547 °C; ΔH=−60.000; ΔG°37=−8.72844 (±1e-3) | primer3-py 2.3.0 (dangling-end/`tstack2`) |
| N6 | Full DP — contiguous regression | GCGCGCGC self-dimer through the full DP | Tm=40.0906 °C; ΔH=−70.8 kcal/mol; BasePairs=8 (±1e-3) | primer3-py; regression of contiguous case |
| N7 | Tm method delegates to full DP | `CalculateDimerMeltingTemperature` on GCGCAAAGCGC/GCGCTTTGCGC | equals full-DP Tm = 41.8816 °C | delegation |
| N8 | Full DP — no structure / invalid | poly-A/poly-A, non-ACGT, null | `CalculateDimerThermodynamicsNtthal` null | ntthal `no_structure`; input contract |

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
| M1–M12, S1–S3, C1 | ✅ Covered | Initial dimer round |
| N1–N8 (full ntthal DP) | ❌ Missing | Added with the full ntthal DP (2026-06-25) |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_DimerTm_Tests.cs` — all dimer-Tm tests.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `PrimerDesigner_DimerTm_Tests.cs` | Canonical dimer-Tm fixture | 27 |

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
| 11 | N1–N5 | ❌ Missing | Implemented (full ntthal DP non-contiguous parity: internal loops, bulge, overhang) | ✅ Done |
| 12 | N6 | ❌ Missing | Implemented (contiguous regression through the full DP) | ✅ Done |
| 13 | N7 | ❌ Missing | Implemented (Tm method delegates to full DP) | ✅ Done |
| 14 | N8 | ❌ Missing | Implemented (no-structure / invalid → null) | ✅ Done |

**Total items:** 14 groups (27 tests)
**✅ Done:** 14 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M12 | ✅ Covered | Exact values, `.Within(1e-9)` (hand-derived) / `.Within(1e-3)` (primer3 parity) |
| S1–S3 | ✅ Covered | Delegation + monotonicity |
| C1 | ✅ Covered | Most-stable selection |
| N1–N5 | ✅ Covered | Full ntthal DP parity (internal loops, bulge, overhang) vs primer3-py 2.3.0, `.Within(1e-3)` |
| N6–N8 | ✅ Covered | Contiguous regression; Tm delegation; no-structure/invalid → null |

---

## 6. Assumption Register

**Total assumptions:** 0

| # | Assumption | Used In |
|---|-----------|---------|
| — | None. The full ntthal dimer DP (internal mismatches/loops, bulges, terminal overhangs via `tstack2`/`dangle`/`stackmm`/`tstack`/loop tables, all verbatim from primer3 config) is implemented and matches primer3-py 2.3.0 to machine precision; every constant/table is source-backed. | — |

---

## 7. Open Questions / Decisions

1. None. The full ntthal dimer DP (internal loops, bulges, terminal-overhang `tstack2` extension) is now implemented and verified against primer3-py 2.3.0 for non-contiguous optima. The only ntthal capability not ported is the optional caller-supplied tri/tetraloop & terminal-mismatch hairpin bonus tables — a hairpin/monomer feature, not part of the dimer model.
