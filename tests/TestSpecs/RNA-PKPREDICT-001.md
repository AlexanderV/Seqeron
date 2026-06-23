# Test Specification: RNA-PKPREDICT-001

**Test Unit ID:** RNA-PKPREDICT-001
**Area:** RnaStructure
**Algorithm:** Pseudoknot Structure Prediction (canonical H-type, pknotsRG class)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-23

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Reeder & Giegerich (2004), BMC Bioinformatics 5:104 (pknotsRG) | 1 | https://doi.org/10.1186/1471-2105-5-104 (PMC514697) | 2026-06-23 |
| 2 | pknotsRG source (Energy.lhs penalties 9.0/0.3/0.0) | 3 | https://github.com/jensreeder/pknotsRG | 2026-06-23 |
| 3 | Pseudoknot — Wikipedia (H-type geometry; Rivas & Eddy 1999) | 4 | https://en.wikipedia.org/wiki/Pseudoknot | 2026-06-23 |
| 4 | Su et al. (1999), Nat Struct Biol 6(3):285–292; PDB 437D (BWYV) | 1 | https://www.rcsb.org/structure/437D | 2026-06-23 |
| 5 | Antczak et al. (2018), Bioinformatics 34(8):1304–1312 (crossing i<k<j<l) | 1 | https://doi.org/10.1093/bioinformatics/btx783 | 2026-06-23 |

### 1.2 Key Evidence Points

1. Canonical simple recursive pseudoknot = two crossing helices a·a', b·b' + three loops u,v,w — Reeder & Giegerich (2004).
2. Pseudoknot penalties: initiation 9.0, unpaired loop base 0.3, in-knot base pair 0.0 kcal/mol — pknotsRG Energy.lhs.
3. H-type 5'→3' order: stem1-5' · loop1 · stem2-5' · loop2 · stem1-3' · loop3 · stem2-3' — Wikipedia.
4. Crossing condition: pairs (i,j),(k,l) cross iff i<k<j<l — Antczak (2018).
5. Helices use Turner 2004 nearest-neighbour stacking (same as nested structures) — pknotsRG.

### 1.3 Documented Corner Cases

- A pseudoknot is reported only when the two crossing helices outweigh the 9 kcal/mol penalty plus the best nested alternative (no spurious knots).
- Tertiary-stabilised knots (BWYV / PDB 437D) are not the NN-thermodynamic MFE structure — documented limit of NN-only predictors.

### 1.4 Known Failure Modes / Pitfalls

1. Treating a strong nested hairpin as a pseudoknot — prevented by the strict-improvement acceptance rule.
2. Double-counting unpaired loop bases that the loop MFE pairs — `ScoreLoop` charges 0.3 only on the nucleotides the loop fold leaves unpaired.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictStructurePseudoknot(string, int)` | RnaSecondaryStructure | Canonical | Deep evidence-based testing |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-PK-01 | Returned ΔG ≤ plain-MFE ΔG for same input | Yes | MFE is the fallback baseline |
| INV-PK-02 | HasPseudoknot ⇒ base pairs contain a crossing pair | Yes | Antczak (2018); construction |
| INV-PK-03 | Each position paired ≤ once; indices in [0,n) | Yes | Canonization rule 1; disjoint spans |
| INV-PK-04 | No spurious pseudoknot (strict improvement only) | Yes | 9 kcal/mol penalty (Reeder & Giegerich 2004) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Designed H-type recovered | `GGGGAACCCCAACCCCAAGGGG` | HasPseudoknot=true; pairs {(0,15)(1,14)(2,13)(3,12)} ∪ {(6,21)(7,20)(8,19)(9,18)}; dot-bracket `((((..[[[[..))))..]]]]` | H-type geometry + pknotsRG penalties |
| M2 | Knot beats MFE | same sequence | FreeEnergy < CalculateMfeStructure(seq).FreeEnergy | INV-PK-01/04 |
| M3 | Crossing genuine | same sequence | DetectPseudoknots over returned pairs yields ≥1 crossing | INV-PK-02; Antczak (2018) |
| M4 | Valid structure | same sequence | every index in [0,n); each position paired ≤ once | INV-PK-03 |
| M5 | No spurious knot (hairpin) | `GGGGAAAACCCC` | HasPseudoknot=false; DotBracket & ΔG = plain MFE (`((((....))))`) | INV-PK-04 |
| M6 | Invariant ΔG ≤ MFE (property) | seeded random sweep (≥100 seqs) | for every seq, FreeEnergy ≤ MFE; any reported knot beats MFE and is valid | INV-PK-01/03/04 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | null input | `null` | empty seq, ".", no pairs, ΔG 0, HasPseudoknot=false | contract parity |
| S2 | empty input | `""` | empty result, HasPseudoknot=false | contract parity |
| S3 | too short | `GGGCCC` (<11 nt) | HasPseudoknot=false; = plain MFE | min knot length |
| S4 | BWYV real knot | `GGCGCGGCACCGUCCGCGGAACAAACGG` | HasPseudoknot=false (tertiary-stabilised, not NN-MFE); ΔG ≤ MFE | documents NN limit |
| S5 | DNA spelling parity | T-form of the M1 sequence | identical pairs/energy as the U-form | T read as U |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | minLoopSize clamp | minLoopSize=0 on M1 seq | clamped to 3; same result as default | NNDB minimum |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- New method; no prior tests. Canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M6 | ❌ Missing | new unit |
| S1–S5 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs` — all cases.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs | Canonical | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented (seeded sweep) | ✅ Done |
| 7 | S1 | ❌ Missing | Implemented | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented | ✅ Done |
| 10 | S4 | ❌ Missing | Implemented | ✅ Done |
| 11 | S5 | ❌ Missing | Implemented | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | DesignedHType_RecoversBothCrossingHelices |
| M2 | ✅ | DesignedHType_BeatsPlainMfe |
| M3 | ✅ | DesignedHType_ContainsGenuineCrossing |
| M4 | ✅ | DesignedHType_ValidStructure |
| M5 | ✅ | PlainHairpin_NoSpuriousPseudoknot |
| M6 | ✅ | RandomSweep_NeverWorseThanMfe_KnotsValid |
| S1 | ✅ | NullInput_EmptyStructure |
| S2 | ✅ | EmptyInput_EmptyStructure |
| S3 | ✅ | TooShort_NoPseudoknot |
| S4 | ✅ | Bwyv_NotRecoveredAsMfe_DocumentsTertiaryLimit |
| S5 | ✅ | DnaInput_FoldsIdenticallyToRna |
| C1 | ✅ | MinLoopSizeZero_ClampedToThree |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | PARTIAL: single canonical H-type knot only (no recursive/multiple knots; loops fold pseudoknot-free) | scope of the unit; §5.3 / Evidence |

---

## 7. Open Questions / Decisions

1. None. The class implemented (canonical single H-type) and its residual (recursive/multiple/non-canonical knots, tertiary stabilisation) are documented in LIMITATIONS.md §1.
