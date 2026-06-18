# Test Specification: RNA-MFE-001

**Test Unit ID:** RNA-MFE-001
**Area:** RnaStructure
**Algorithm:** Minimum Free Energy (Zuker–Stiegler DP with Turner 2004 nearest-neighbor parameters)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Zuker & Stiegler (1981), *Nucleic Acids Res.* 9(1):133–148 | 1 | https://doi.org/10.1093/nar/9.1.133 (full text PMC326673) | 2026-06-14 |
| 2 | Mathews et al. (2004), *PNAS* 101:7287–7292 (Turner 2004 params) | 1 | https://doi.org/10.1073/pnas.0401799101 | 2026-06-14 |
| 3 | Lorenz et al. (2011), ViennaRNA Package 2.0, *Algorithms Mol. Biol.* 6:26 | 3 | https://doi.org/10.1186/1748-7188-6-26 (PMC3319429) | 2026-06-14 |
| 4 | Ward et al. (2017), *Nucleic Acids Res.* 45(14):8541–8552 | 1 | https://doi.org/10.1093/nar/gkx512 (PMC5737859) | 2026-06-14 |
| 5 | NNDB Turner 2004 Hairpin Example 1 | 2 | https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-example-1.html | 2026-06-14 |
| 6 | NNDB Turner 2004 Hairpin Example 2 | 2 | https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-example-2.html | 2026-06-14 |
| 7 | NNDB Turner 2004 Watson-Crick Helix Example 2 | 2 | https://rna.urmc.rochester.edu/NNDB/turner04/wc-nsc-example.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. MFE folding is a dynamic program over loop decomposition (hairpin / stacked pair / bulge-interior / multibranch); total ΔG° is the sum of loop contributions — Zuker & Stiegler (1981).
2. Standard affine-multiloop Turner model runs in **O(n³) time, O(n²) space**; matrices C(i,j) (closed/paired), M/M1 (multiloop), F (exterior) — Ward et al. (2017); Lorenz et al. (2011).
3. A single hairpin's ΔG°37 = stem stacking + AU/GU end penalties + hairpin-loop energy; NNDB Example 1 = **−1.4** (sum −1.41), Example 2 = **−1.9** (sum −1.91) — NNDB hairpin examples.
4. For a unimolecular fold, the intermolecular helix-initiation term is **not** added — NNDB hairpin-example-1 note.
5. Pure helix nearest-neighbor stacking sum for 5'-GCACG/3'-CGUGC = **−10.13** (no AU end penalty; both ends G-C) — NNDB WC Example 2.

### 1.3 Documented Corner Cases

- Minimum hairpin loop = 3 nt; a pair (i,j) with `j−i−1 < 3` cannot close a hairpin; sequence shorter than `minLoopSize+2` forms no pair — Zuker & Stiegler / nearest-neighbor rules.
- Unfoldable sequence (no possible base pair) ⇒ MFE = 0 (open chain is always available, ΔG = 0).
- `null` / empty input ⇒ MFE = 0; `PredictStructure("")` ⇒ empty result.

### 1.4 Known Failure Modes / Pitfalls

1. Adding the bimolecular helix-initiation constant to a unimolecular hairpin would over-stabilize by a fixed offset — NNDB hairpin-example-1 (intermolecular initiation omitted).
2. Multiloop model choice changes both the optimum and complexity (affine O(n³) vs logarithmic O(n⁴)) — Ward et al. (2017). Seqeron uses the affine model.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateMinimumFreeEnergy(string, int)` | `RnaSecondaryStructure` | **Canonical** | Zuker-style O(n³) DP; deep evidence-based MFE tests. |
| `PredictStructure(string, int, int, int)` | `RnaSecondaryStructure` | **Canonical** | Greedy stem-loop structure assembly returning base pairs + dot-bracket; tested for pairing/dot-bracket on the worked example and empty/no-structure modes. |
| `CalculateMinimumFreeEnergyClassic(string, int)` | `RnaSecondaryStructure` | **Internal** | O(n³) baseline using a simplified per-pair energy model; exercised only as a timing comparator in the benchmark (INV-03 = determinism), not separately re-tested. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | MFE ≤ 0 for every input (open chain ΔG = 0 is always available; optimum cannot be positive). | Yes | Zuker & Stiegler (1981) |
| INV-02 | MFE is non-increasing under suffix extension: `MFE(prefix)` ≥ `MFE(prefix+suffix)`. | Yes | Zuker & Stiegler (1981) — extension only adds folding options |
| INV-03 | MFE score is deterministic (repeated evaluation yields identical value). | Yes | benchmark determinism assertion (RnaSecondaryStructure_MFE_Benchmark). The classic baseline uses a different simplified energy model and is NOT expected to match the Turner-model score. |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | `CalculateMinimumFreeEnergy_NndbHairpinExample1_ReturnsMinus1_41` | Sequence `CACAAAAAAAUGUG` folds to NNDB Example 1 stem-loop | −1.41 (NNDB −1.4) | NNDB hairpin-example-1.html |
| M2 | `CalculateMinimumFreeEnergy_NndbHairpinExample2_ReturnsMinus1_91` | Sequence `CACAGAAAGUGUG`, GG first mismatch | −1.91 (NNDB −1.9) | NNDB hairpin-example-2.html |
| M3 | `CalculateMinimumFreeEnergy_HomopolymerNoPairs_ReturnsZero` | `AAAAAAAA` has no complementary bases | 0 | Zuker–Stiegler (open chain) |
| M4 | `CalculateMinimumFreeEnergy_EmptyOrNull_ReturnsZero` | empty/null input | 0 | corner case |
| M5 | `CalculateMinimumFreeEnergy_ShorterThanMinLoop_ReturnsZero` | `GCGC` < `minLoopSize+2` | 0 | min-hairpin-loop rule |
| M6 | `CalculateMinimumFreeEnergy_Invariant_NeverPositive` | INV-01 on several deterministic sequences | every result ≤ 0 | Zuker–Stiegler |
| M7 | `CalculateMinimumFreeEnergy_Invariant_MonotonicUnderExtension` | INV-02: extend a hairpin sequence | `MFE(extended)` ≤ `MFE(prefix)` | Zuker–Stiegler |
| M8 | `PredictStructure_NndbHairpinExample1_PairsAndDotBracket` | structure of `CACAAAAAAAUGUG` | 4 pairs C-G/A-U/C-G/A-U, dot-bracket `((((......))))` | NNDB hairpin-example-1 structure |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | `PredictStructure_EmptySequence_ReturnsEmptyResult` | `""` | 0 pairs, empty dot-bracket, MFE 0 | documented empty mode |
| S2 | `PredictStructure_HomopolymerNoStructure_AllDots` | `AAAAAAAAAA` | 0 pairs, dot-bracket all `.` | no-structure mode |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | `CalculateMinimumFreeEnergy_GcStem_MoreStableThanAuStem` | a GC-rich stem is more stable (more negative) than an AU stem of equal geometry | `MFE(GC stem)` < `MFE(AU stem)` | stacking magnitudes (GC ≫ AU); performance baseline recorded via benchmark |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructureTests.cs` contains pre-existing MFE tests (`CalculateMinimumFreeEnergy_SimpleHairpin_ReturnsNegative` using `Is.LessThan(0)`; `_NoStructure_ReturnsZero`; `_EmptySequence_ReturnsZero`; `_LongerStem_MoreStable`) — all **permissive** (no exact evidence value).
- `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_MFE_Benchmark.cs` — `[Explicit]` performance baseline + optimized-engine determinism check (INV-03). The classic baseline uses a different simplified energy model, so the two scores are not expected to match.
- No canonical `{Class}_{Method}_Tests.cs` file existed for MFE.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (Example 1 exact MFE) | ❌ Missing | only `Is.LessThan(0)` existed |
| M2 (Example 2 exact MFE) | ❌ Missing | — |
| M3 (homopolymer → 0) | ⚠ Weak | `_NoStructure_ReturnsZero` exists but generic; rewrite with evidence note |
| M4 (empty/null → 0) | ⚠ Weak | `_EmptySequence_ReturnsZero` exists; rewrite in canonical file |
| M5 (too-short → 0) | ❌ Missing | — |
| M6 (INV-01 ≤ 0) | ❌ Missing | — |
| M7 (INV-02 monotonic) | ❌ Missing | — |
| M8 (PredictStructure pairs/dot-bracket) | ❌ Missing | no exact-structure MFE test |
| S1 (PredictStructure empty) | ❌ Missing | — |
| S2 (PredictStructure homopolymer) | ❌ Missing | — |
| C1 (GC vs AU stability) | ⚠ Weak | `_LongerStem_MoreStable` exists but tests length, not GC vs AU |
| INV-03 (determinism) | ✅ Covered | benchmark determinism assertion |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_MinimumFreeEnergy_Tests.cs` — all M/S/C cases for `CalculateMinimumFreeEnergy` and `PredictStructure`, evidence-based exact values.
- **Remove:** the four weak MFE tests in `RnaSecondaryStructureTests.cs` (`_SimpleHairpin_ReturnsNegative`, `_NoStructure_ReturnsZero`, `_EmptySequence_ReturnsZero`, `_LongerStem_MoreStable`) — superseded by the canonical file (no duplication left behind).
- **Keep:** `RnaSecondaryStructure_MFE_Benchmark.cs` (performance baseline + INV-03 determinism); corrected the broken classic-vs-optimized equality assertion (the two engines use different energy models) to a determinism check.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `RnaSecondaryStructure_MinimumFreeEnergy_Tests.cs` | canonical unit tests | 11 |
| `RnaSecondaryStructure_MFE_Benchmark.cs` | performance baseline + INV-03 determinism | 1 (`[Explicit]`) |
| `RnaSecondaryStructureTests.cs` | other methods (4 weak MFE tests removed) | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented exact −1.41 test | ✅ Done |
| 2 | M2 | ❌ Missing | implemented exact −1.91 test | ✅ Done |
| 3 | M3 | ⚠ Weak | rewritten in canonical file with evidence note | ✅ Done |
| 4 | M4 | ⚠ Weak | rewritten in canonical file | ✅ Done |
| 5 | M5 | ❌ Missing | implemented too-short test | ✅ Done |
| 6 | M6 | ❌ Missing | implemented INV-01 property test | ✅ Done |
| 7 | M7 | ❌ Missing | implemented INV-02 property test | ✅ Done |
| 8 | M8 | ❌ Missing | implemented PredictStructure pairs/dot-bracket | ✅ Done |
| 9 | S1 | ❌ Missing | implemented PredictStructure empty | ✅ Done |
| 10 | S2 | ❌ Missing | implemented PredictStructure homopolymer | ✅ Done |
| 11 | C1 | ⚠ Weak | rewritten as GC vs AU stability | ✅ Done |
| 12 | (old weak tests) | 🔁 Duplicate | removed 4 weak MFE tests from RnaSecondaryStructureTests.cs | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact −1.41 (NNDB Ex1) |
| M2 | ✅ | exact −1.91 (NNDB Ex2) |
| M3 | ✅ | homopolymer → 0 |
| M4 | ✅ | empty/null → 0 |
| M5 | ✅ | too-short → 0 |
| M6 | ✅ | INV-01 ≤ 0 |
| M7 | ✅ | INV-02 monotonic |
| M8 | ✅ | dot-bracket `((((......))))`, 4 pairs |
| S1 | ✅ | empty result |
| S2 | ✅ | all-dots, 0 pairs |
| C1 | ✅ | GC stem < AU stem |
| INV-03 | ✅ | benchmark determinism |

Total in-scope cases: 12 → ✅ = 12.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Multiloop per-unpaired cost `c = 0` (affine model simplification; offset/helix terms source-backed) | DP multiloop term; does not affect any cited hairpin worked example |
| 2 | Rounding to 2 decimals vs NNDB 1 decimal | M1/M2 assert exact sum `.Within(1e-9)` and `Math.Round(mfe,1)` == NNDB total |

---

## 7. Open Questions / Decisions

1. The MFE worked examples in scope are single hairpins (the only published per-term ΔG breakdowns retrievable). Multiloop/coaxial-stacking optima are not validated against a published numeric example (no full numeric MFE example with a single unambiguous optimum was retrievable for an intramolecular multiloop without coaxial stacking, which the DP does not model). This is captured by Assumption 1 and the multiloop simplification block in the algorithm doc; correctness of the single-hairpin path (the cited evidence) is fully validated.
