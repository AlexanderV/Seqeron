# RNA-ENERGY-001 Evidence Document

**Test Unit:** RNA-ENERGY-001  
**Title:** Free Energy Calculation  
**Area:** RnaStructure  
**Created:** 2026-02-05  
**Status:** Complete

---

## 1. Overview

This document provides evidence for the RNA free energy calculation algorithm testing, including thermodynamic parameters from the Turner 2004 nearest-neighbor model, test methodology, and reference data.

---

## 2. Authoritative Sources

### 2.1 Primary References

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia - Nucleic acid thermodynamics | Encyclopedia | Nearest-neighbor method, stacking interactions, free energy |
| Wikipedia - Nucleic acid secondary structure | Encyclopedia | Secondary structure prediction, thermodynamic models |
| Turner (2004) | Thermodynamic parameters | Definitive RNA nearest-neighbor parameters at 37°C |
| Mathews et al. (2004) | Peer-reviewed paper | PNAS 101(19):7287-7292, computational RNA folding |
| NNDB (Nearest Neighbor Database) | Reference database | https://rna.urmc.rochester.edu/NNDB/ |
| Xia et al. (1998) | Peer-reviewed paper | Biochemistry 37(42):14719-35, expanded NN parameters |
| SantaLucia (1998) | Peer-reviewed paper | PNAS 95(4):1460-5, unified view of NN thermodynamics |

### 2.2 Implementation References

| Source | Usage |
|--------|-------|
| ViennaRNA | Reference implementation, Turner parameters |
| MFOLD | Comparison tool for energy calculations |
| RNAstructure (Mathews Lab) | Reference implementation |

---

## 3. Algorithm Specification

### 3.1 Nearest-Neighbor Model

The nearest-neighbor (NN) model predicts RNA folding free energy by summing contributions from adjacent base pair stacks. The total free energy is:

**ΔG°total = ΔG°initiation + Σ ΔG°stacking + Σ ΔG°loops**

Where:
- **Initiation**: Cost of forming the first base pair
- **Stacking**: Stabilizing energy from adjacent base pairs
- **Loops**: Destabilizing energy from unpaired regions (hairpin, internal, bulge loops)

### 3.2 Stacking Energy Parameters (Turner 2004)

**Watson-Crick Stacking (ΔG°37 in kcal/mol):**

| Stack (5'→3'/3'→5') | ΔG°37 (kcal/mol) |
|---------------------|------------------|
| AA/UU | -0.93 |
| AU/UA | -1.10 |
| UA/AU | -1.33 |
| CU/GA (=AG/UC) | -2.08 |
| CA/GU (=UG/AC) | -2.11 |
| GU/CA (=AC/UG) | -2.24 |
| GA/CU (=UC/AG) | -2.35 |
| CG/GC | -2.36 |
| GG/CC | -3.26 |
| GC/CG | -3.42 |

**Source:** NNDB Turner 2004 (https://rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html)

**Key properties:**
- All Watson-Crick stacking energies are negative (stabilizing)
- GC-rich stacks are more stable (more negative)
- Stacking requires at least 2 adjacent base pairs
- Single base pair has no stacking contribution (ΔG° = 0)

### 3.3 G-U Wobble Pair Stacking

G-U wobble pairs have reduced stability compared to Watson-Crick:

| Stack | ΔG°37 (kcal/mol) |
|-------|------------------|
| GU/UG | -1.3 (approximate) |
| UG/GU | -1.3 (approximate) |

### 3.4 Hairpin Loop Initiation Energy

**Loop size-dependent initiation (ΔG°37 in kcal/mol):**

| Loop Size (nt) | ΔG°37 (kcal/mol) |
|----------------|------------------|
| 3 | +5.4 |
| 4 | +5.6 |
| 5 | +5.7 |
| 6 | +5.4 |
| 7 | +6.0 |
| 8 | +5.5 |
| 9 | +6.4 |
| ≥10 | Extrapolation formula |

**Important:** Loop initiation energies are **NOT monotonically increasing** with loop size. For example, size 6 (+5.4) is more stable than size 4 (+5.6), and size 8 (+5.5) is more stable than size 4 (+5.6). This reflects experimental thermodynamic measurements where optimal loop sizes (around 4-6 nt) can have varying energies based on conformational constraints.

**Extrapolation for n > 9:**
ΔG°37(n) = ΔG°37(9) + 1.75 × R × T × ln(n/9)

Where R = 1.987 cal/(mol·K), T = 310.15 K (37°C)

**Source:** NNDB Turner 2004 Hairpin Loops

### 3.5 Special Tetraloop Stability

Certain 4-nucleotide loops have exceptional stability due to tertiary interactions:

| Tetraloop | Bonus (kcal/mol) | Motif Type |
|-----------|------------------|------------|
| GAAA | -3.0 | GNRA |
| GCAA | -3.0 | GNRA |
| GGAA | -3.0 | GNRA |
| GUAA | -3.0 | GNRA |
| UUCG | -3.0 | UNCG |
| UACG | -3.0 | UNCG |
| UGCG | -3.0 | UNCG |
| UCCG | -3.0 | UNCG |
| CUUG | -2.0 | CUYG |
| CCUG | -2.0 | CUYG |

**Source:** Heus & Pardi (1991), Science 253(5016):191-194; Turner 2004 parameters

### 3.6 All-C Loop Penalty

Poly-C loops are destabilized due to electrostatic repulsion:

**ΔG°penalty = A × n + B**

Where n is loop size, and A ≈ 0.3 kcal/mol per C.

---

## 4. Test Methodology

### 4.1 Method Coverage (RNA-ENERGY-001)

| Method | Type | Test Focus |
|--------|------|------------|
| `CalculateStemEnergy` | Energy calculation | Stacking energy for stem |
| `CalculateHairpinLoopEnergy` | Energy calculation | Hairpin loop destabilization |
| `CalculateMinimumFreeEnergy` | Structure energy | Total MFE prediction |
| `CalculateStackingEnergy` | Helper | Individual base pair stack |

### 4.2 Reference Energy Values

**Test Case 1: GC Stem (3 bp)**
- Sequence: 5'-GCG...loop...CGC-3'
- Stacks: GC/CG + CG/GC
- Expected: ΔG° = -3.42 + -2.36 = -5.78 kcal/mol (approximately)

**Test Case 2: GNRA Tetraloop Hairpin**
- Full structure: GGGG-GAAA-CCCC (4bp stem + GAAA loop)
- Stem energy: GG/CC × 3 = -3.26 × 3 ≈ -9.78 kcal/mol
- Loop initiation: +5.6 kcal/mol
- GNRA bonus: -3.0 kcal/mol
- Total: ≈ -7.18 kcal/mol (approximately)

### 4.3 Edge Cases

| Case | Expected Behavior | Rationale |
|------|------------------|-----------|
| Empty sequence | ΔG° = 0 | No structure |
| Single base pair | ΔG° = 0 (no stacking) | Stacking requires 2 adjacent pairs |
| Poly-A (no structure) | ΔG° = 0 | No complementary pairing |
| All-C loop | Higher energy (penalty) | Electrostatic destabilization |
| Very long loop (>30) | Extrapolated energy | Jacobson-Stockmayer formula |

### 4.4 Invariants

1. **Negative stacking**: All Watson-Crick stacking energies must be negative
2. **Positive loop initiation**: Hairpin loop initiation is always positive (destabilizing)
3. **Monotonic loop penalty**: Longer loops have higher initiation penalties (to a point)
4. **GNRA stability**: GNRA tetraloops are more stable than random tetraloops
5. **All-C penalty**: Poly-C loops are less stable than equivalent poly-A loops

---

## 5. Corner Cases

### 5.1 Minimal Valid Cases

| Case | Description | Expected |
|------|-------------|----------|
| Minimum stem | 2 base pairs | Stem energy = one stacking term |
| Minimum loop | 3 nucleotides | Sterically minimum loop |
| Single pair | 1 base pair | Stacking energy = 0 |

### 5.2 Boundary Cases

| Case | Description | Expected |
|------|-------------|----------|
| Loop size = 3 | Minimum allowed | No terminal mismatch term |
| Loop size = 4 | Tetraloop region | Check for special bonus |
| Loop size = 30 | Extrapolation boundary | Use Jacobson-Stockmayer |

### 5.3 Error Cases

| Case | Expected |
|------|----------|
| Null input | Return 0 or handle gracefully |
| Empty string | Return 0 |
| Invalid bases | Undefined (implementation choice) |

---

## 6. ASSUMPTIONS

The following behaviors are not explicitly defined by primary sources and are implementation decisions:

1. **Default stacking energy**: When exact NN parameter not found, use -1.5 kcal/mol (ASSUMPTION)
2. **G-U wobble parameters**: Use approximate values when exact not available (ASSUMPTION)
3. **Rounding**: Energy values are rounded to 2 decimal places (ASSUMPTION)
4. **Temperature**: All calculations assume 37°C (310.15 K) standard conditions

---

## 7. References

1. Turner, D.H. (2004). Thermodynamics of RNA. In: RNA World, 3rd Ed.
2. Mathews, D.H., et al. (2004). PNAS 101(19):7287-7292.
3. Xia, T., et al. (1998). Biochemistry 37(42):14719-35.
4. SantaLucia, J. (1998). PNAS 95(4):1460-5.
5. NNDB: https://rna.urmc.rochester.edu/NNDB/
6. Heus, H.A. & Pardi, A. (1991). Science 253(5016):191-194.
7. Wikipedia - Nucleic acid thermodynamics
8. Wikipedia - Nucleic acid secondary structure
