# Test Specification: RNA-PKRECURSIVE-001

**Test Unit ID:** RNA-PKRECURSIVE-001 (recursive-grammar extension of RNA-STRUCT-001 / RNA-PKPREDICT-001)
**Area:** RnaStructure
**Algorithm:** Recursive pknotsRG pseudoknot prediction (nested / multiple H-type pseudoknots)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-23

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Reeder & Giegerich 2004, BMC Bioinformatics 5:104 (pknotsRG; recursive class; penalties 9.0 / 0.3 / 0.0) | 1 | https://doi.org/10.1186/1471-2105-5-104 (https://pmc.ncbi.nlm.nih.gov/articles/PMC514697/) | 2026-06-23 |
| 2 | Reeder, Steffen & Giegerich 2007, NAR 35:W320 (per-interval competition; whole-sequence DP) | 1 | https://doi.org/10.1093/nar/gkm258 (https://pmc.ncbi.nlm.nih.gov/articles/PMC1933184/) | 2026-06-23 |
| 3 | pknotsRG source (Energy.lhs — pseudoknot penalties) | 3 | https://github.com/jensreeder/pknotsRG | 2026-06-23 |
| 4 | Antczak et al. 2018, Bioinformatics 34(8):1304 (crossing condition i<k<j<l) | 1 | https://doi.org/10.1093/bioinformatics/btx783 | 2026-06-23 |

### 1.2 Key Evidence Points

1. A simple recursive pseudoknot's three loops u, v, w "fold internally in an arbitrary way, **including simple recursive pseudoknots**" — loops may contain further knots — Source 1.
2. The pseudoknot value "competes with values of unknotted foldings for the interval (i, j)", folding the whole sequence so the optimum may contain several / nested knots — Source 2.
3. Energy: Turner stacking for both helices; initiation 9.0 kcal/mol; 0.3 kcal/mol per unpaired pseudoknot-loop nt; 0.0 per in-knot pair — Sources 1, 3.
4. Two base pairs (i,j),(k,l) cross iff i<k<j<l — Source 4.

### 1.3 Documented Corner Cases

- Spurious-knot suppression via the 9 kcal/mol penalty (a knot is taken only when it lowers ΔG for that interval) — Source 1.
- Two simultaneous strong G·C knots are the genuine MFE only when each region is isolated (e.g. A·U clamps), else a cross-region nested fold wins — energy-model property documented in the Evidence; two-knot tests are therefore engineered, not random.
- Excluded classes: kissing hairpins, triple-crossing / chained ("complex") helices, bulged/unequal-length pseudoknot helices — Source 1.

### 1.4 Known Failure Modes / Pitfalls

1. Asserting a universal "recursive beats single-knot" on random input would be false thermodynamics — assert on engineered cases only — Source 1.
2. Tertiary-stabilised knots (BWYV) are not the NN-thermodynamic MFE — an energy-model floor, not an algorithm gap — Source 1.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictStructurePseudoknotRecursive(string, int)` | `RnaSecondaryStructure` | **Canonical** | Recursive pknotsRG folding; nested / multiple / over-arching H-type knots. |
| `RecursivePkFolder.Fold` / `EvaluateHTypeRecursive` / `ScoreLoopRecursive` | `RnaSecondaryStructure` (private) | **Internal** | Tested indirectly via the public method. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Recursive ΔG ≤ `CalculateMfeStructure(seq)` ΔG for any sequence (MFE is the always-available fallback). | Yes | Evidence (fallback baseline) |
| INV-2 | If `HasPseudoknot`, the structure has ≥1 genuine crossing (i<k<j<l). | Yes | Source 4 |
| INV-3 | Every index is in range and each position is paired ≤ once. | Yes | rule 1, Source 1 |
| INV-4 | On a plain (non-pseudoknotted) sequence no spurious knot is introduced (ΔG ≤ MFE, no crossing). | Yes | Source 1 (9 kcal/mol penalty) |
| INV-5 | A reported pseudoknot strictly beats the plain MFE (else the MFE is returned). | Yes | Source 1 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Over-arching nested knot | `AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU` | structure `((((((((((((..[[[[..))))..]]]])))))))`, pairs as in Evidence dataset, `HasPseudoknot`, ΔG = −14.37 | Source 1 (loop recursion) |
| M2 | Nested knot beats single-knot method | same sequence | recursive ΔG (−14.37) < `PredictStructurePseudoknot` ΔG (−13.05); single method `HasPseudoknot == false` | Source 1 |
| M3 | Two separate knots recovered | `AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUUAAAAAAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU` | recursive recovers TWO crossing knots (crossing-count = 32), ΔG = −28.74 | Source 2 |
| M4 | Two-knot beats single-knot method | same 80-nt sequence | recursive ΔG (−28.74) < single/MFE (−27.14); single method `HasPseudoknot == false` (recovers none) | Source 2 |
| M5 | No spurious knot (hairpin) | `GGGGAAAACCCC` | recursive = MFE `((((....))))`, ΔG = −5.28, `HasPseudoknot == false`, crossing = 0 | Source 1 (INV-4) |
| M6 | INV-1 / INV-2 / INV-3 sweep | 150 random seqs (seed 20260623, len 12–38) | recursive ΔG ≤ MFE for all; any reported knot crosses and is valid; no spurious knots | INV-1/2/3 |
| M7 | Validity of recovered structure | M1 + M3 sequences | every index in range, each position paired ≤ once | INV-3 |
| M8 | Genuine crossing | M1 + M3 sequences | `DetectPseudoknots` finds ≥1 (M1) / ≥2 (M3) crossings | Source 4 (INV-2) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | null input | `null` | empty pseudoknot-free structure (no pairs, ΔG 0, `HasPseudoknot` false) | parity with `PredictStructurePseudoknot` |
| S2 | empty input | `""` | empty structure | parity |
| S3 | too short (< 11 nt) | `GGGCCC` | no knot; equals plain MFE | min knot length |
| S4 | single canonical H-type parity | `GGGGAACCCCAACCCCAAGGGG` | identical to `PredictStructurePseudoknot` (db, pairs, ΔG −8.76, `HasPseudoknot`) | recursion must not regress the single case |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | DNA spelling | M1 sequence with T for U | identical knot decision / pairs / ΔG | T read as U |
| C2 | minLoopSize < 3 clamped | M1 sequence, minLoopSize 0 | identical to default | NNDB minimum 3 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior tests for `PredictStructurePseudoknotRecursive` (new method). The sibling `RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs` covers the single-knot method only.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M8, S1–S4, C1–C2 | ❌ Missing | new method; no existing coverage |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_PredictStructurePseudoknotRecursive_Tests.cs` — all cases for this unit.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `RnaSecondaryStructure_PredictStructurePseudoknotRecursive_Tests.cs` | canonical | 14 |

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
| 9 | S1 | ❌ Missing | implemented | ✅ Done |
| 10 | S2 | ❌ Missing | implemented | ✅ Done |
| 11 | S3 | ❌ Missing | implemented | ✅ Done |
| 12 | S4 | ❌ Missing | implemented | ✅ Done |
| 13 | C1 | ❌ Missing | implemented | ✅ Done |
| 14 | C2 | ❌ Missing | implemented | ✅ Done |

**Total items:** 14
**✅ Done:** 14 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | nested-knot recovery + exact pairs |
| M2 | ✅ | beats single-knot method |
| M3 | ✅ | two knots recovered (crossing = 32) |
| M4 | ✅ | beats single-knot method (which recovers none) |
| M5 | ✅ | no spurious knot on hairpin |
| M6 | ✅ | random sweep INV-1/2/3 |
| M7 | ✅ | validity |
| M8 | ✅ | genuine crossings |
| S1 | ✅ | null |
| S2 | ✅ | empty |
| S3 | ✅ | too short |
| S4 | ✅ | single H-type parity |
| C1 | ✅ | DNA spelling |
| C2 | ✅ | minLoopSize clamp |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | PARTIAL coverage of the full pknotsRG O(n⁴) yield parser (faithful recursive *class* — nested/multiple/over-arching — via explicit helix scan + memoised interval recursion, not a bit-identical reproduction of every reference yield). | scope of M-cases |
| 2 | Two-simultaneous-knot cases are engineered (isolating A·U clamps), not random — asserting a universal random "beats single-knot" would be false thermodynamics. | M3, M4 |

---

## 7. Open Questions / Decisions

1. None. Kissing hairpins, triple-crossing / chained knots, and tertiary-stabilised knots remain out of scope and are documented in LIMITATIONS.md (RNA-STRUCT-001) and the algorithm doc §5.3 / §6.2.
