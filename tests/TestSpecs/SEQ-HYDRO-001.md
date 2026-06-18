# Test Specification: SEQ-HYDRO-001

**Test Unit ID:** SEQ-HYDRO-001
**Area:** Statistics
**Algorithm:** Hydrophobicity Analysis (Kyte-Doolittle GRAVY + sliding-window profile)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Kyte & Doolittle (1982), J Mol Biol 157:105–132 | 1 | https://doi.org/10.1016/0022-2836(82)90515-0 | 2026-06-13 |
| 2 | Biopython ProtParamData.py (kd scale) | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/ProtParamData.py | 2026-06-13 |
| 3 | Biopython ProtParam.py (gravy, protein_scale) | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/ProtParam.py | 2026-06-13 |
| 4 | Expasy ProtParam documentation (GRAVY) | 2 | https://web.expasy.org/protparam/protparam-doc.html | 2026-06-13 |
| 5 | GCAT Davidson — Kyte-Doolittle background | 4 | https://gcat.davidson.edu/DGPB/kd/kyte-doolittle-background.htm | 2026-06-13 |
| 6 | alakazam (CRAN) gravy doc | 3 | https://rdrr.io/cran/alakazam/man/gravy.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. Kyte-Doolittle scale values (A 1.8 … V 4.2, all 20 residues) — Biopython `kd` dict (src 2), cited to Kyte & Doolittle 1982 (src 1).
2. GRAVY = sum of hydropathy values / number of residues — Expasy doc (src 4) and Biopython `gravy()` = `total_gravy / self.length` (src 3).
3. Sliding-window profile returns N − W + 1 values, equal weight per position (edge=1.0) — Biopython `protein_scale` (src 3).
4. Window 9 for surface, window 19 with peaks > 1.6 for transmembrane segments — GCAT/Kyte-Doolittle (src 5).

### 1.3 Documented Corner Cases

- Window W > sequence length N → empty profile (`range(N−W+1)` yields 0 iterations) — Biopython (src 3).
- Unknown residues (B, Z, X, gaps) are undefined in the scale; Biopython raises `KeyError`. Kyte-Doolittle/Expasy define only the 20 canonical residues (src 1, 4).

### 1.4 Known Failure Modes / Pitfalls

1. Using a weighted window instead of an unweighted mean would change profile values — Biopython default is `edge=1.0` (unweighted) (src 3).
2. Dividing GRAVY by the full string length when unknown residues are present vs by the recognized-residue count is an undefined-input choice not fixed by sources (src 1, 4) → tracked as an assumption.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateHydrophobicity(string)` | SequenceStatistics | **Canonical** | GRAVY = sum(kd)/count |
| `CalculateHydrophobicityProfile(string, int windowSize=9)` | SequenceStatistics | **Canonical** | N−W+1 unweighted window means |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Empty/null sequence → GRAVY 0; empty profile | Yes | Biopython `range(N−W+1)`=0; contract |
| INV-2 | GRAVY of a single recognized residue equals that residue's kd value | Yes | GRAVY=sum/length, length 1 (src 3,4) |
| INV-3 | Profile length = N − W + 1 for W ≤ N, else 0 | Yes | Biopython `protein_scale` loop (src 3) |
| INV-4 | GRAVY is case-insensitive | Yes | scale on uppercase; impl uppercases input |
| INV-5 | Each profile value is the unweighted mean of its window's kd values | Yes | Biopython edge=1.0 (src 3) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | GRAVY single residue A | "A" → 1.8/1 | 1.8 | kd A=1.8 (src 2); GRAVY def (src 4) |
| M2 | GRAVY hydrophobic FLIV | (2.8+3.8+4.5+4.2)/4 | 3.825 | kd (src 2); GRAVY def (src 4) |
| M3 | GRAVY hydrophilic RKDE | (−4.5−3.9−3.5−3.5)/4 | −3.85 | kd (src 2); GRAVY def (src 4) |
| M4 | GRAVY case-insensitive | "fliv" == "FLIV" | 3.825 | INV-4 |
| M5 | GRAVY empty | "" → 0 | 0 | INV-1 |
| M6 | GRAVY null | null → 0 | 0 | INV-1 / contract |
| M7 | Profile FLIV W=3 | 2 windows: means | [3.7, 4.16666666667] | kd (src 2); protein_scale (src 3) |
| M8 | Profile length AG W=3 | W>N | empty | INV-3 (src 3) |
| M9 | Profile empty/null | "" / null → empty | empty | INV-1 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | GRAVY skips unknown residues | "AX" → divide by recognized count (1) | 1.8 | documents deviation from Biopython KeyError |
| S2 | Profile length general | "FLIVAG" W=3 → 4 windows | count = 4 | INV-3 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Transmembrane-style window | hydrophobic 19-mer "I"×19, W=19 → peak | 4.5 (> 1.6) | GCAT threshold (src 5) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatisticsTests.cs` — legacy aggregate fixture; grepped for `Hydrophobicity`: no dedicated, evidence-based hydrophobicity tests exist.
- No `SequenceStatistics_CalculateHydrophobicity_Tests.cs` present before this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M9, S1–S2, C1 | ❌ Missing | new unit; no prior evidence-based hydrophobicity tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateHydrophobicity_Tests.cs` — all GRAVY + profile cases.
- **Remove:** nothing (no pre-existing dedicated tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| SequenceStatistics_CalculateHydrophobicity_Tests.cs | canonical | 12 |

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
| 10 | S1 | ❌ Missing | implemented | ✅ Done |
| 11 | S2 | ❌ Missing | implemented | ✅ Done |
| 12 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact value 1.8 |
| M2 | ✅ Covered | exact value 3.825 |
| M3 | ✅ Covered | exact value −3.85 |
| M4 | ✅ Covered | case-insensitive |
| M5 | ✅ Covered | empty → 0 |
| M6 | ✅ Covered | null → 0 |
| M7 | ✅ Covered | [3.7, 4.16666666667] |
| M8 | ✅ Covered | W>N empty |
| M9 | ✅ Covered | empty/null empty |
| S1 | ✅ Covered | unknown skipped |
| S2 | ✅ Covered | profile length |
| C1 | ✅ Covered | W=19 peak 4.5 |

Total in-scope cases: 12; ✅ = 12.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Unknown residues are skipped (GRAVY divides by recognized count; profile treats as 0) rather than raising as in Biopython. Not correctness-affecting for the 20 canonical residues; sources define only those. | S1 |

---

## 7. Open Questions / Decisions

1. None. All GRAVY/profile values over the 20 standard residues are source-conformant; the only undefined-input behavior (non-standard residues) is documented as a deviation and does not affect canonical values.
