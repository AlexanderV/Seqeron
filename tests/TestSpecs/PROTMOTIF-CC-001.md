# Test Specification: PROTMOTIF-CC-001

**Test Unit ID:** PROTMOTIF-CC-001
**Area:** ProteinMotif
**Algorithm:** Coiled-Coil Prediction (heptad-repeat a/d hydrophobic-core detection)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Mason JM, Arndt KM (2004), ChemBioChem 5(2):170–176 | 1 | https://doi.org/10.1002/cbic.200300781 | 2026-06-14 |
| 2 | Lupas A, Van Dyke M, Stock J (1991), Science 252:1162–1164 | 1 | https://doi.org/10.1126/science.252.5009.1162 | 2026-06-14 |
| 3 | Wikipedia "Coiled coil" (cites Mason & Arndt 2004) | 4 | https://en.wikipedia.org/wiki/Coiled_coil | 2026-06-14 |
| 4 | Wikipedia "Heptad repeat" (cites Chambers et al. 1990) | 4 | https://en.wikipedia.org/wiki/Heptad_repeat | 2026-06-14 |

### 1.2 Key Evidence Points

1. Heptad repeat denoted (abcdefg)n; positions **a** and **d** are the hydrophobic core. — Mason & Arndt (2004); Lupas (1991).
2. a/d positions "often being occupied by isoleucine, leucine, or valine" → hydrophobic-core set **{I, L, V}**. — Wikipedia/Coiled coil (Mason & Arndt 2004).
3. The register is unknown a priori; **seven heptad frames** must be tried and the best taken; canonical scoring **window = 28 residues** (4 heptads). — Lupas (1991).
4. Pattern hxxhcxc: h at a and d. — Wikipedia/Coiled coil; Wikipedia/Heptad repeat (HPPHCPC).

### 1.3 Documented Corner Cases

- Sequence shorter than the window cannot be scored → no prediction (Lupas window rule).
- Single hydrophobic position is not a coiled coil; multiple heptads are required (Mason & Arndt: (abcdefg)1-2-3…).
- Correct register unknown → evaluate all seven and take the max (Lupas).

### 1.4 Known Failure Modes / Pitfalls

1. Using a fixed register (frame 0) misses coiled coils offset from that frame — must max over 7 registers. — Lupas (1991).
2. Fabricated position-specific weights (a COILS-style PSSM) are not retrievable and must not be invented — this unit uses closed-form a/d occupancy instead. — unit directive.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictCoiledCoils(sequence, windowSize, threshold)` | ProteinMotifFinder | **Canonical** | Heptad a/d occupancy over 7 registers; returns (Start, End, Score) regions. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned `Score` ∈ [0, 1] (fraction of a/d positions occupied by {I,L,V}). | Yes | §1.2.2 + definition |
| INV-2 | Returned regions have `End − Start + 1 ≥ MinRegion (21)`. | Yes | §1.2 Mason & Arndt multi-heptad |
| INV-3 | `0 ≤ Start ≤ End ≤ length − 1`; regions are non-overlapping and in increasing Start order. | Yes | definition |
| INV-4 | Sequences shorter than `windowSize` return no regions. | Yes | §1.3 Lupas window |
| INV-5 | A returned region exists only if some covering window scores ≥ threshold (max over 7 registers). | Yes | §1.2.3 Lupas registers |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Perfect heptad | `("LAALAAA"×5)` (35 aa); a=L, d=L every heptad, register 0. W=28, thr=0.5. | Exactly one region `(0, 34, 1.0)`. | §1.2.1–.2; INV-1,3 |
| M2 | No core residues | `"G"×40` — no {I,L,V}. | No regions (score 0 < 0.5). | §1.2.2 |
| M3 | Below window length | `"LAALAAA"×3` (21 aa) with W=28. | Empty (no full window). | §1.3 Lupas; INV-4 |
| M4 | Off-frame register | Perfect heptad prefixed by 2 alanines: `"AA" + "LAALAAA"×5` (37 aa). Coiled core is in register r=2. | One region found via the 7-register max, score 1.0, covering the core. | §1.2.3 Lupas; INV-5 |
| M5 | Region below MinRegion | A single scoring window's residue span = exactly W (28) ≥ 21, so it IS kept; a construction whose contiguous scoring residue span < 21 is rejected. Use `"LAALAAA"×4` (28 aa) → one window, region length 28 (kept). Contrast: shortening below window already covered by M3; sub-MinRegion-but-≥window is impossible (W=28>21), so verify the boundary that any kept region length ≥ 21. | M5: `"LAALAAA"×4` → exactly one region `(0,27,1.0)`, length 28 ≥ 21. | INV-2,3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Half occupancy boundary | `"LAAAAAA"×5` (35 aa): L only at a, not d. Best register has 50% of a/d = L → score exactly 0.5 = threshold (≥) → region kept. | One region `(0,34,0.5)`. | Threshold boundary (≥). |
| S2 | Just below threshold | Same as S1 but threshold = 0.5001. | No region (0.5 < 0.5001). | Strict threshold check. |
| S5 | Core + hydrophilic tail | `"LAALAAA"×5 + "G"×35` (70 aa): core scores, tail does not — region must end inside the sequence. | One region `(0,48,1.0)` (ends at 48, before 69). | Mid-sequence drop-and-emit branch; INV-03. |
| S3 | Null input | `null`. | Empty. | Validation. |
| S4 | Empty input | `""`. | Empty. | Validation. |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Case-insensitive | lowercase `("laalaaa"×5)`. | Same as M1: `(0,34,1.0)`. | Input uppercased. |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `ProteinMotifFinder.PredictCoiledCoils` existed in `src/.../ProteinMotifFinder.cs` with a **fabricated position-specific weight table** (per-residue weights 0.9/0.8/0.7… for L/I/V/M/A and charged-position weights). No authoritative source provides those constants → invented values = defect (per unit directive, the COILS PSSM must not be fabricated). No existing canonical test file for this unit (`grep` of tests dir: none referencing `PredictCoiledCoils`).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new |
| M2 | ❌ Missing | new |
| M3 | ❌ Missing | new |
| M4 | ❌ Missing | new |
| M5 | ❌ Missing | new |
| S1 | ❌ Missing | new |
| S2 | ❌ Missing | new |
| S5 | ❌ Missing | new |
| S3 | ❌ Missing | new |
| S4 | ❌ Missing | new |
| C1 | ❌ Missing | new |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_PredictCoiledCoils_Tests.cs` — all cases for this unit.
- **Remove:** the fabricated weight table in the implementation (replaced by source-traceable a/d occupancy). No test files to remove.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ProteinMotifFinder_PredictCoiledCoils_Tests.cs` | Canonical | 11 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | S1 | ❌ Missing | Implemented | ✅ Done |
| 7 | S2 | ❌ Missing | Implemented | ✅ Done |
| 8 | S5 | ❌ Missing | Implemented | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented | ✅ Done |
| 10 | S4 | ❌ Missing | Implemented | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact region (0,34,1.0) |
| M2 | ✅ | empty |
| M3 | ✅ | empty (sub-window) |
| M4 | ✅ | off-frame found via 7 registers |
| M5 | ✅ | exact region (0,27,1.0) |
| S1 | ✅ | score exactly 0.5 |
| S2 | ✅ | empty above threshold |
| S5 | ✅ | region (0,48,1.0) ends inside sequence |
| S3 | ✅ | null → empty |
| S4 | ✅ | empty → empty |
| C1 | ✅ | lowercase → (0,34,1.0) |

**✅ count = 11 = total in-scope cases.**

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Hydrophobic-core set = {I, L, V} (verbatim from Wikipedia/Mason & Arndt 2004). | M1, M2, M4, S1, C1, impl |
| 2 | Defaults window=28, MinRegion=21, threshold=0.5 (Lupas 1991 window/registers; Mason & Arndt multi-heptad). | all, impl |

Both are source-grounded (the named residue set and the documented window/heptad counts); they are documented modeling choices with caller-overridable parameters, not invented constants.

---

## 7. Open Questions / Decisions

1. **COILS PSSM deliberately not implemented:** the 21×20 position-specific residue-frequency table from Lupas (1991) was not retrievable in this session; per the unit directive it must not be fabricated. This unit implements the fully-specified heptad-register a/d hydrophobic-core detection (registers + 28-window reused from Lupas) instead. Documented in algorithm doc §5.3.
