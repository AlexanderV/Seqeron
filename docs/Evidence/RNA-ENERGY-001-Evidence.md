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

G-U wobble pairs have variable stability. All 20 entries from NNDB Turner 2004 are used:

| Stack | ΔG°37 (kcal/mol) | Notes |
|-------|-----------------|---------|
| AG/UU | -0.55 | |
| UG/AU | -1.00 | |
| GA/UU | -1.27 | |
| AU/UG | -1.36 | |
| CG/GU | -1.41 | |
| GG/CU | -1.53 | |
| CU/GG | -2.11 | |
| GU/CG | -2.51 | |
| GG/UU | -0.50 | Note a: set to -0.5 for prediction |
| UG/GU | +0.30 | Destabilizing |
| GU/UG | +1.29 | Destabilizing |

**Per GU end:** +0.45 kcal/mol (same as per AU end)

**Special GGUC/CUGG 3-stack context (NNDB note b):** When GU/UG is flanked by GC and CG
(5'GGUC/3'CUGG), the total for all 3 stacking interactions is **-4.12 kcal/mol**,
replacing the individual sum (-1.53 + 1.29 + (-1.53) = -1.77).

**Source:** NNDB Turner 2004 GU Parameters (rna.urmc.rochester.edu/NNDB/turner04/gu-parameters.html)

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

### 3.5 Special Hairpin Loops (NNDB Total Energies)

Special hairpin loops of 3, 4, and 6 nucleotides have stabilities poorly fit by the standard
model. NNDB provides **total experimental energies** that replace the model calculation entirely.
The key includes the closing base pair (e.g., CUUCGG = C-G closing + UUCG loop).

**Tetraloops (key = closing5' + loop + closing3'):**

| Key | ΔG°37 (kcal/mol) | Loop Motif |
|------|------------------|------------|
| CCUCGG | 2.5 | UNCG (most stable) |
| CUCCGG | 2.7 | UNCG |
| CUACGG | 2.8 | UNCG |
| CUGCGG | 2.8 | UNCG |
| CCAAGG | 3.3 | GNRA |
| CCCAGG | 3.4 | GNRA |
| CCGAGG | 3.5 | GNRA |
| CUUAGG | 3.5 | |
| CCGCGG | 3.6 | |
| CUAAGG | 3.6 | GNRA |
| CCUAGG | 3.7 | |
| CCACGG | 3.7 | |
| CUCAGG | 3.7 | |
| CUUCGG | 3.7 | UNCG |
| CUUUGG | 3.7 | |
| CAACGG | 5.5 | |

**Triloops:** CAACG (6.8), GUUAC (6.9)

**Hexaloops:** ACAGUGUU (1.8), ACAGUACU (2.8), ACAGUGCU (2.9), ACAGUGAU (3.6)

**Source:** NNDB Turner 2004 (rna.urmc.rochester.edu/NNDB/turner04/tloop.txt, triloop.txt, hexaloop.txt)

### 3.6 All-C Loop Penalty

Poly-C loops are destabilized due to electrostatic repulsion:

- **3nt all-C loops:** flat penalty = +1.5 kcal/mol
- **>3nt all-C loops:** penalty = 0.3n + 1.6 kcal/mol

**Source:** NNDB Turner 2004 (rna.urmc.rochester.edu/NNDB/turner04/hairpin-mismatch-parameters.html)

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

**Test Case 2: Inner AU Terminal Penalty (NNDB Hairpin Example 1)**
- Sequence: GGGAUAAAUCCC (4bp stem GC,GC,GC,AU + 4nt loop UAAA)
- Stacking: GG/CC(-3.26) + GG/CC(-3.26) + GA/CU(-2.35) = -8.87 kcal/mol
- Hairpin: init(4)=5.6 + tm(AUAU=-0.6) = 5.0 kcal/mol
- Inner AU end penalty: +0.45 kcal/mol (per NNDB)
- Total: -8.87 + 5.0 + 0.45 = **-3.42 kcal/mol**

**Test Case 3: GGUC/CUGG Special 3-Stack Context**
- 4 pairs: G-C, G-U, U-G, C-G (5'GGUC/3'CUGG)
- NNDB note b: total for 3 stacks = **-4.12 kcal/mol**
- Without special context: GG/CU(-1.53) + GU/UG(+1.29) + UC/GG(-1.53) = -1.77
- Source: NNDB GU parameters (rna.urmc.rochester.edu/NNDB/turner04/gu-parameters.html)

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

## 6. VERIFICATION STATUS

All parameters verified against NNDB Turner 2004 (rna.urmc.rochester.edu/NNDB/turner04/):

| Parameter Set | Entries | Status | Source URL |
|---------------|---------|--------|------------|
| WC stacking | 16 | ✅ Exact match | /wc-parameters.html |
| GU stacking | 20 | ✅ Exact match (incl. note a) | /gu-parameters.html |
| GGUC/CUGG 3-stack | 1 | ✅ Implemented (-4.12) | /gu-parameters.html (note b) |
| Hairpin initiation (3–30) | 28 | ✅ Exact match | /loop.txt |
| Terminal mismatch | 96 | ✅ Exact match | /tm-parameters.html |
| Special hairpins (tri/tetra/hexa) | 22 | ✅ Exact match | /tloop.txt, triloop.txt, hexaloop.txt |
| Hairpin bonuses/penalties | 6 | ✅ Exact match | /hairpin-mismatch-parameters.html |
| Terminal AU/GU penalty | 1 | ✅ +0.45 (both WC and GU) | /wc-parameters.html, /gu-parameters.html |
| Inner AU/GU penalty in MFE | — | ✅ Applied per NNDB example | /hairpin-example-1.html |
| Special GU closure in MFE | — | ✅ Detected from sequence | /hairpin-mismatch-parameters.html |

**No assumptions remain.** All parameters are sourced from NNDB Turner 2004.

**Standard conditions:** 37°C (310.15 K) — defined by Turner 2004 (not an assumption).
**Precision:** 2 decimal places — follows NNDB parameter precision.
**Unknown stacking pairs:** Contribute 0.0 kcal/mol — non-canonical interactions are outside the NN model scope.

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
