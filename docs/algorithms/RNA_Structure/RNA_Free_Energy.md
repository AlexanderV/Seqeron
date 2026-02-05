# RNA Free Energy Calculation

**Algorithm Group:** RNA_Structure  
**Test Unit:** RNA-ENERGY-001  
**Implementation:** `Seqeron.Genomics.Analysis.RnaSecondaryStructure`

---

## 1. Overview

RNA free energy calculation predicts the thermodynamic stability of RNA secondary structures using the nearest-neighbor model. The implementation uses Turner 2004 parameters, which are the standard for RNA folding predictions.

---

## 2. Theory

### 2.1 Nearest-Neighbor Model

The nearest-neighbor (NN) model treats RNA structure as a series of interactions between adjacent base pairs. The total free energy of a structure is the sum of:

1. **Stacking energies**: Stabilizing interactions between adjacent base pairs
2. **Loop penalties**: Destabilizing energies for unpaired regions
3. **Special motif contributions**: Bonuses/penalties for specific sequences

**Total Free Energy:**
```
ΔG°total = Σ ΔG°stacking + Σ ΔG°loops + Σ ΔG°special
```

### 2.2 Stacking Interactions

Stacking energies arise from:
- Van der Waals forces between adjacent bases
- Electrostatic interactions
- Hydrophobic effects

All Watson-Crick stacking energies are **negative** (stabilizing):
- GC-rich stacks: Most stable (most negative)
- AU-rich stacks: Less stable
- G-U wobble pairs: Reduced stability

### 2.3 Loop Energies

Hairpin loops have **positive** initiation energies (destabilizing):
- Minimum loop size: 3 nucleotides (steric constraint)
- Energy increases with loop size (to a limit)
- Special tetraloops (GNRA, UNCG) have stability bonuses

---

## 3. Implementation Details

### 3.1 Methods

| Method | Description |
|--------|-------------|
| `CalculateStemEnergy` | Sum of stacking energies for adjacent base pairs |
| `CalculateHairpinLoopEnergy` | Loop initiation + special bonuses/penalties |
| `CalculateMinimumFreeEnergy` | Dynamic programming for global MFE |

### 3.2 Parameters (Turner 2004)

**Stacking energies (kcal/mol at 37°C):**

| Stack | ΔG°37 |
|-------|-------|
| GC/CG | -3.42 |
| CG/GC | -2.36 |
| GG/CC | -3.26 |
| GA/CU | -2.35 |
| GU/CA | -2.24 |
| CA/GU | -2.11 |
| CU/GA | -2.08 |
| UA/AU | -1.33 |
| AU/UA | -1.10 |
| AA/UU | -0.93 |

**Hairpin loop initiation (kcal/mol):**

| Size | ΔG°37 |
|------|-------|
| 3 | +5.4 |
| 4 | +5.6 |
| 5 | +5.7 |
| 6 | +5.4 |

**Special tetraloops:**
- GNRA (GAAA, GCAA, etc.): -3.0 kcal/mol bonus
- UNCG: -3.0 kcal/mol bonus

### 3.3 MFE Algorithm

The Minimum Free Energy (MFE) calculation uses a simplified Nussinov-style dynamic programming approach:

1. Build DP table for all subsequences
2. For each position, consider:
   - Unpaired: Inherit energy from smaller subsequence
   - Paired: Add base pair energy + recursive subproblems
3. Return minimum energy from full sequence

**Time Complexity:** O(n³)  
**Space Complexity:** O(n²)

---

## 4. Biological Significance

### 4.1 Structure Stability

- Lower (more negative) free energy = more stable structure
- Structures with ΔG° ≈ 0 are unlikely to form
- Typical stable RNA hairpins: ΔG° < -5 kcal/mol

### 4.2 Special Motifs

**GNRA Tetraloops:**
- Extremely stable hairpin loops
- Found in ribosomal RNA, ribozymes
- Stabilized by unusual base stacking

**All-C Loops:**
- Destabilized by electrostatic repulsion
- Poly-C sequences avoid structure formation

---

## 5. Limitations

1. **Simplified model**: Does not include all Turner 2004 parameters
2. **No pseudoknots**: DP algorithm cannot handle crossing base pairs efficiently
3. **Temperature fixed**: Assumes 37°C (310.15 K)
4. **Ionic conditions**: Assumes standard buffer conditions

---

## 6. References

1. Turner, D.H. (2004). Thermodynamics of RNA. In: The RNA World, 3rd Ed.
2. Mathews, D.H., et al. (2004). Incorporating chemical modification constraints into a dynamic programming algorithm. PNAS 101(19):7287-7292.
3. Xia, T., et al. (1998). Thermodynamic parameters for an expanded nearest-neighbor model. Biochemistry 37(42):14719-35.
4. NNDB: https://rna.urmc.rochester.edu/NNDB/
5. Wikipedia - Nucleic acid thermodynamics
