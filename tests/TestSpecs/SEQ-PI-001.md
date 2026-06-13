# Test Specification: SEQ-PI-001

**Test Unit ID:** SEQ-PI-001
**Area:** Statistics
**Algorithm:** Isoelectric Point (pI) Calculation
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | EMBOSS `iep` documentation (pKa table, pI definition) | 3 | https://emboss.sourceforge.net/emboss/apps/iep.html | 2026-06-13 |
| 2 | Peptides R pkg `charge_pI.cpp` (net charge formula) | 3 | https://raw.githubusercontent.com/cran/Peptides/master/src/charge_pI.cpp | 2026-06-13 |
| 3 | Peptides R pkg `charge()` doc (EMBOSS worked example) | 3 | https://rdrr.io/cran/Peptides/man/charge.html | 2026-06-13 |
| 4 | seqinr `computePI` doc (Bjellqvist worked example) | 3 | https://rdrr.io/cran/seqinr/man/computePI.html | 2026-06-13 |
| 5 | Bjellqvist et al. 1993, Electrophoresis 14:1023–1031 | 1 | https://doi.org/10.1002/elps.11501401163 | 2026-06-13 |

### 1.2 Key Evidence Points

1. pI is the pH where net charge = 0; found over [0,14] — EMBOSS iep [1].
2. EMBOSS pKa: Nterm 8.6, Cterm 3.6, C 8.5, D 3.9, E 4.1, H 6.5, K 10.8, R 12.5, Y 10.1 — EMBOSS iep [1].
3. Net charge: basic groups `+1/(1+10^(pH−pKa))`, acidic groups `−1/(1+10^(pKa−pH))` (Henderson–Hasselbalch, Moore 1985) — Peptides `charge_pI.cpp` [2].
4. Acidic = D,E,C,Y,Cterm; basic = R,K,H,Nterm — Peptides [2].
5. EMBOSS worked example: `FLPVLAGLTPSIVPKLVCLLTKKC` net charge 3.037398 / 2.914112 / 0.7184524 at pH 5/7/9 — Peptides `charge()` doc [3].

### 1.3 Documented Corner Cases

- "No electrostatic interactions": pI depends only on amino-acid composition, not order — EMBOSS iep [1].
- Small / highly basic proteins: predicted pI may be inaccurate (accuracy caveat, not a correctness rule) — ExPASy [doc].

### 1.4 Known Failure Modes / Pitfalls

1. Using a different pKa scale (Bjellqvist vs EMBOSS) gives a different pI for the same sequence; expected values are scale-specific — seqinr vs EMBOSS [1,4].
2. pI of an empty/zero-length protein is undefined in literature (a real protein always has termini) — no source defines it.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateIsoelectricPoint(string)` | SequenceStatistics | Canonical | pI via bisection on EMBOSS-scale net charge |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | 0 ≤ pI ≤ 14 for any input | Yes | Bisection interval [0,14] — EMBOSS iep [1] |
| INV-2 | pI is composition-only: permuting residues leaves pI unchanged | Yes | "no electrostatic interactions" — EMBOSS iep [1] |
| INV-3 | Net charge is monotonically non-increasing in pH (more acid → lower pI) | Yes | Henderson–Hasselbalch sign structure — Peptides [2] |
| INV-4 | Termini-only sequence → pI = midpoint of N/C-term pKa = (8.6+3.6)/2 = 6.10 | Yes | EMBOSS pKa [1] |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Charge of reference peptide @ pH 5/7/9 | Net charge of `FLPVLAGLTPSIVPKLVCLLTKKC` (exposed via internal/charge path or pI cross-check) | pI ≈ 9.67 (basic; charge +0.72 at pH 9) | Peptides EMBOSS [3] |
| M2 | pI bounds | Any sequence yields pI in [0,14] | 0 ≤ pI ≤ 14 | EMBOSS iep [1] |
| M3 | Termini-only `A` | No ionizable side chains | pI = 6.10 (±0.01) | EMBOSS pKa [1] |
| M4 | Acidic-only `DDDD` | Four Asp + termini | pI = 3.23 (±0.01) | Derived from EMBOSS pKa [1] |
| M5 | Basic-only `KKKK` | Four Lys + termini | pI = 11.27 (±0.01) | Derived from EMBOSS pKa [1] |
| M6 | All-20 `ACDEFGHIKLMNPQRSTVWY` | One of each residue, EMBOSS scale | pI = 7.36 (±0.01) | Derived from EMBOSS pKa [1] |
| M7 | Empty string | Input guard | 7.0 | ASSUMPTION (input-guard convention) |
| M8 | Null | Input guard | 7.0 | ASSUMPTION (input-guard convention) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Single Asp `D` | one acidic residue | pI = 3.75 (±0.01) | Derived EMBOSS pKa |
| S2 | Single Lys `K` | one basic residue | pI = 9.70 (±0.01) | Derived EMBOSS pKa |
| S3 | Lowercase input | Case-insensitivity | pI(`dddd`) = pI(`DDDD`) | Normalization convention |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Order-independence | Permutation of residues | pI(`DK`) = pI(`KD`) | INV-2 property test |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior test file for SEQ-PI-001. `CalculateIsoelectricPoint` existed in `SequenceStatistics.cs` with unsourced pKa values; no dedicated test fixture. (`SequenceStatistics_CalculateAminoAcidComposition*` indirectly touches the record field but does not assert pI values.)

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | New unit |
| M2 | ❌ Missing | New unit |
| M3 | ❌ Missing | New unit |
| M4 | ❌ Missing | New unit |
| M5 | ❌ Missing | New unit |
| M6 | ❌ Missing | New unit |
| M7 | ❌ Missing | New unit |
| M8 | ❌ Missing | New unit |
| S1 | ❌ Missing | New unit |
| S2 | ❌ Missing | New unit |
| S3 | ❌ Missing | New unit |
| C1 | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateIsoelectricPoint_Tests.cs` — all SEQ-PI-001 cases.
- **Remove:** none (no prior fixture).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| SequenceStatistics_CalculateIsoelectricPoint_Tests.cs | Canonical | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented (pI ≈ 9.67 for ref peptide) | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented (bounds 0..14) | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented (`A` → 6.10) | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented (`DDDD` → 3.23) | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented (`KKKK` → 11.27) | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented (all-20 → 7.36) | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented (empty → 7.0) | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented (null → 7.0) | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented (`D` → 3.75) | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented (`K` → 9.70) | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented (case-insensitive) | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented (order-independence) | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Covered |
| M2 | ✅ | Covered |
| M3 | ✅ | Covered |
| M4 | ✅ | Covered |
| M5 | ✅ | Covered |
| M6 | ✅ | Covered |
| M7 | ✅ | Covered |
| M8 | ✅ | Covered |
| S1 | ✅ | Covered |
| S2 | ✅ | Covered |
| S3 | ✅ | Covered |
| C1 | ✅ | Covered |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Empty/null → 7.0 (input-guard convention; pI undefined for zero-length protein) | M7, M8 |
| 2 | pKa scale = EMBOSS (single pKa per residue; Bjellqvist value not used as expected) | M1, M3–M6, S1, S2 |

---

## 7. Open Questions / Decisions

1. Decision: target the EMBOSS pKa scale (matches the repository's single-pKa-per-residue model), not Bjellqvist; recorded in Evidence Assumption 2. No open questions remain.
