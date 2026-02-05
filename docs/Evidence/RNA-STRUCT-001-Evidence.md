# RNA-STRUCT-001 Evidence Document

**Test Unit:** RNA-STRUCT-001  
**Title:** Secondary Structure Prediction  
**Area:** RnaStructure  
**Created:** 2026-02-05  
**Status:** Complete

---

## 1. Overview

This document provides evidence for the RNA secondary structure prediction algorithm testing, including sources, reference data, test methodology, and corner cases.

---

## 2. Authoritative Sources

### 2.1 Primary References

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia - Nucleic acid structure prediction | Encyclopedia | Algorithm overview, dynamic programming, thermodynamic models |
| Wikipedia - Nussinov algorithm | Encyclopedia | O(n³) DP algorithm for maximizing base pairs |
| Wikipedia - Nucleic acid secondary structure | Encyclopedia | Base pairing rules, stem-loop motifs, dot-bracket notation |
| Nussinov & Jacobson (1980) | Peer-reviewed paper | Original Nussinov algorithm, PNAS 77(11):6309-6313 |
| Zuker & Stiegler (1981) | Peer-reviewed paper | MFE-based prediction, Nucleic Acids Research 9(1):133-148 |
| Turner (2004) | Thermodynamic parameters | Nearest-neighbor parameters for RNA folding |
| Mathews et al. (2004) | Peer-reviewed paper | PNAS 101(19):7287-7292, modern prediction constraints |

### 2.2 Secondary References

| Source | Usage |
|--------|-------|
| Rosetta Code | Algorithm implementation examples |
| ViennaRNA documentation | Reference implementation, energy parameters |
| MFOLD web server | Comparison tool for structure prediction |

---

## 3. Algorithm Specification

### 3.1 Core Concepts

**Secondary Structure Definition:**
- Base pairing interactions within an RNA molecule
- Watson-Crick pairs: A-U, G-C
- Wobble pairs: G-U (valid in RNA)
- Represented in dot-bracket notation

**Nussinov Algorithm:**
- Dynamic programming approach
- Time complexity: O(n³)
- Space complexity: O(n²)
- Maximizes number of base pairs
- Cannot detect pseudoknots

**MFE (Minimum Free Energy) Approach:**
- Uses thermodynamic parameters
- Turner 2004 nearest-neighbor model
- Returns structure with lowest free energy

### 3.2 Structural Motifs

| Motif | Description | Biological Significance |
|-------|-------------|------------------------|
| Stem | Double-stranded helix region | Structural stability |
| Hairpin Loop | Unpaired bases between stem | Protein binding sites |
| Internal Loop | Unpaired bases within helix | Flexibility |
| Bulge | One-sided unpaired bases | Structural kinks |
| Multi-loop | Junction of multiple stems | Complex structures like tRNA |
| Pseudoknot | Crossing base pairs | Catalytic activity |

### 3.3 Energy Parameters (Turner 2004)

**Stacking Energies (kcal/mol at 37°C):**
- GC/CG: -3.4
- CG/GC: -2.4
- GG/CC: -3.3
- AU/UA: -1.1
- UA/AU: -1.3
- GU/UG (wobble): -1.3

**Hairpin Loop Initiation:**
- 3 nt loop: +5.4 kcal/mol
- 4 nt loop: +5.6 kcal/mol
- 5 nt loop: +5.7 kcal/mol
- 6 nt loop: +5.4 kcal/mol

**Tetraloop Bonuses (GNRA, UNCG):**
- GAAA: -3.0 kcal/mol
- UUCG: -3.0 kcal/mol

---

## 4. Test Methodology

### 4.1 Method Coverage (RNA-STRUCT-001)

| Method | Type | Test Focus |
|--------|------|------------|
| `Predict(sequence)` | Canonical | Complete structure prediction |
| `PredictWithConstraints(seq, constraints)` | Constrained | Forced base pairs |
| `ToDotBracket(structure)` | Notation | Output format |
| `FromDotBracket(notation)` | Parse | Input format parsing |

### 4.2 Test Categories

**Must Tests (Evidence-based):**
1. Base pairing correctness (Watson-Crick, Wobble)
2. Stem-loop detection with known structures
3. MFE calculation validation
4. Dot-bracket notation correctness
5. Pseudoknot detection

**Should Tests:**
1. Complex structures (tRNA-like)
2. Energy parameter accuracy
3. Multiple structure handling
4. Non-overlapping structure selection

**Could Tests:**
1. Performance on long sequences
2. Comparison with external tools

---

## 5. Reference Test Data

### 5.1 Simple Hairpin

**Sequence:** `GGGGAAAACCCC`

**Expected Structure:**
- Dot-bracket: `((((....))))`
- Stem length: 4 bp
- Loop size: 4 nt
- MFE: negative (stabilizing)

**Source:** Basic complementarity, Wikipedia stem-loop

### 5.2 GNRA Tetraloop (Stable Hairpin)

**Sequence:** `GCGCGAAACGCGC` (G-GAAA-C closing)

**Expected:**
- Tetraloop bonus: -3.0 kcal/mol
- GAAA is a GNRA tetraloop (N=A, R=A)
- Extra stability vs non-GNRA loops

**Source:** Heus & Pardi (1991), Wikipedia tetraloop

### 5.3 tRNA-like Structure

**Sequence:** `GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA`

**Expected:**
- Multiple stem-loops
- Cloverleaf-like structure
- Valid dot-bracket notation

**Source:** Wikipedia transfer RNA structure

### 5.4 No Structure (Poly-A)

**Sequence:** `AAAAAAAAAAAA`

**Expected:**
- No base pairs
- MFE = 0
- Empty dot-bracket (all dots)

**Source:** A cannot pair with A

---

## 6. Edge Cases and Corner Cases

### 6.1 Input Validation

| Case | Input | Expected Behavior |
|------|-------|-------------------|
| Empty sequence | `""` | Return empty structure, MFE=0 |
| Null sequence | `null` | Return empty structure, MFE=0 |
| Too short | `"GC"` | No stem-loop possible |
| Case insensitivity | `"gggaaaaccc"` | Same result as uppercase |
| Invalid characters | `"GCXGC"` | Handle gracefully or reject |

### 6.2 Structural Edge Cases

| Case | Description | Expectation |
|------|-------------|-------------|
| Minimum hairpin | 3 bp stem + 3 nt loop | Should be detected |
| Maximum hairpin | Long sequence | Multiple non-overlapping |
| All-C loop penalty | Loop of only C's | Higher energy (penalty) |
| Wobble-only stem | G-U pairs only | Should work if wobble enabled |
| No complement | Poly-A, Poly-U | No structure |

### 6.3 Pseudoknot Detection

| Case | Pairs | Expected |
|------|-------|----------|
| Non-crossing | (0,5), (1,4) | No pseudoknot |
| Crossing | (0,6), (3,9) | Pseudoknot detected |
| Nested | (0,10), (2,8), (4,6) | No pseudoknot |

---

## 7. Invariants

### 7.1 Structural Invariants

1. **Dot-bracket balance:** Opening brackets = closing brackets
2. **MFE sign:** Structured RNA has MFE ≤ 0
3. **Stem energy:** Stacking interactions always negative (stabilizing)
4. **Loop energy:** Initiation always positive (destabilizing)
5. **No overlap:** Selected stem-loops must not overlap

### 7.2 Base Pairing Invariants

1. Watson-Crick pairs: A-U, U-A, G-C, C-G only
2. Wobble pairs: G-U, U-G only
3. Non-pairing: A-A, A-G, A-C, U-U, U-C, C-C, G-G

### 7.3 Probability Invariants

1. Structure probability: 0 ≤ P ≤ 1
2. MFE structure: highest probability
3. Boltzmann probability: exp(-E/RT) / Z

---

## 8. Known Limitations

1. **No pseudoknot prediction:** Standard DP cannot find pseudoknots efficiently
2. **Simplified energy model:** Implementation uses simplified Turner parameters
3. **Single sequence only:** No comparative structure prediction
4. **Minimum loop size:** Typically minimum 3 nt for biological relevance

---

## 9. Assumptions

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Turner 2004 parameters are sufficient | Widely used, well-validated |
| A2 | Minimum loop size of 3 | Biological constraint |
| A3 | No modified bases | Standard ACGU only |
| A4 | Greedy non-overlapping selection | Practical simplification |

---

## 10. Verification Checklist

- [x] Watson-Crick base pairing verified
- [x] Wobble G-U pairing verified
- [x] Stem-loop finding algorithm tested
- [x] MFE calculation verified
- [x] Dot-bracket notation validated
- [x] Pseudoknot detection tested
- [x] Edge cases covered
- [x] Invariants verified
- [x] Energy parameters from Turner 2004

---

## 11. References

1. Nussinov R, Jacobson AB (1980). "Fast algorithm for predicting the secondary structure of single-stranded RNA." PNAS 77(11):6309-6313.
2. Zuker M, Stiegler P (1981). "Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information." Nucleic Acids Res 9(1):133-148.
3. Mathews DH et al. (2004). "Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure." PNAS 101(19):7287-7292.
4. Wikipedia contributors. "Nucleic acid structure prediction." Wikipedia.
5. Wikipedia contributors. "Nucleic acid secondary structure." Wikipedia.
6. Wikipedia contributors. "Nussinov algorithm." Wikipedia.
