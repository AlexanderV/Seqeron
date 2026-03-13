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

**MFE (Minimum Free Energy) DP:**
- Nussinov-style DP with weighted pair scores (WC: −2.0, Wobble: −1.0)
- Maximizes weighted pair count (not thermodynamic MFE)
- O(n³) time, O(n²) space
- Results indicate relative stability, not physical energy units

**Stem-Loop Energy Model:**
- Uses Turner 2004 nearest-neighbor stacking parameters from NNDB
- Terminal AU/GU penalty: +0.45 per helix end that terminates with AU/UA or GU/UG pair (NNDB)
- Hairpin loop initiation energies from NNDB (sizes 3–30; Jacobson–Stockmayer extrapolation for >30)
- Special hairpin loops (tri/tetra/hexaloops) replace model calculation with NNDB experimental data
- Mismatch bonuses: UU/GA first mismatch (−0.9), GG first mismatch (−0.8) per NNDB
- All-C loop penalty: +1.5 (3nt), 0.3n+1.6 (n>3) per NNDB

### 3.2 Structural Motifs

| Motif | Description | Biological Significance |
|-------|-------------|------------------------|
| Stem | Double-stranded helix region | Structural stability |
| Hairpin Loop | Unpaired bases between stem | Protein binding sites |
| Internal Loop | Unpaired bases within helix | Flexibility |
| Bulge | One-sided unpaired bases | Structural kinks |
| Multi-loop | Junction of multiple stems | Complex structures like tRNA |
| Pseudoknot | Crossing base pairs | Catalytic activity |

### 3.3 Energy Parameters (Turner 2004 — NNDB)

**Source:** rna.urmc.rochester.edu/NNDB/turner04/

**Watson-Crick Stacking Energies (kcal/mol at 37°C):**

| Stack (5’XY3’/3’X’Y’5’) | ΔG°37 |
|------|------|
| GC/CG | −3.42 |
| GG/CC | −3.26 |
| CG/GC | −2.36 |
| GA/CU | −2.35 |
| GU/CA | −2.24 |
| CA/GU | −2.11 |
| CU/GA | −2.08 |
| UA/AU | −1.33 |
| AU/UA | −1.10 |
| AA/UU | −0.93 |

**GU Wobble Stacking Energies (selected, kcal/mol at 37°C):**

| Stack | ΔG°37 |
|------|------|
| GU/CG | −2.51 |
| CU/GG | −2.11 |
| GG/CU | −1.53 |
| CG/GU | −1.41 |
| AU/UG | −1.36 |
| GA/UU | −1.27 |
| UG/AU | −1.00 |
| AG/UU | −0.55 |
| GG/UU | −0.50 |
| UG/GU | +0.30 |
| GU/UG | +1.29 |

**Terminal AU/GU Penalty:**

| Condition | Penalty | Source |
|-----------|---------|--------|
| Per AU/UA helix end | +0.45 | NNDB wc-parameters.html |
| Per GU/UG helix end | +0.45 | NNDB gu-parameters.html |

**Hairpin Loop Initiation (kcal/mol at 37°C):**

| Size | Energy | Size | Energy |
|------|--------|------|--------|
| 3 | +5.4 | 4 | +5.6 |
| 5 | +5.7 | 6 | +5.4 |
| 7 | +6.0 | 8 | +5.5 |
| 9 | +6.4 | 10–30 | NNDB extrapolated |

Extrapolation for n>9: ΔG°(n) = 6.4 + 1.75·R·T·ln(n/9) where R=1.987 cal/(mol·K), T=310.15K

**Special Hairpin Loops (total energies replacing model, selected):**

| Sequence | Loop | Closing | ΔG°37 |
|----------|------|---------|------|
| CCUCGG | CUCG (UNCG) | C-G | +2.5 |
| CUACGG | UACG (UNCG) | C-G | +2.8 |
| CCGAGG | CGAG (GNRA) | C-G | +3.5 |
| CUUCGG | UUCG (UNCG) | C-G | +3.7 |

**Mismatch Bonuses:**

| Type | ΔG°37 |
|------|------|
| UU or GA first mismatch | −0.9 |
| GG first mismatch | −0.8 |
| All-C loop (3nt) | +1.5 |
| All-C loop (n>3) | +0.3n + 1.6 |

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
- GA first mismatch bonus: −0.9 kcal/mol (NNDB)
- GNRA tetraloops with C-G closing are in the NNDB special table;
  with G-C closing, standard model + GA mismatch bonus applies

**Source:** NNDB Turner 2004, Heus & Pardi (1991)

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
2. **MFE sign:** Structured RNA has MFE ≤ 0 (Nussinov weighted pair score)
3. **Stem energy:** Watson-Crick stacking contributions are always negative per NNDB
4. **Loop energy:** Initiation always positive (destabilizing)
5. **No overlap:** Selected stem-loops must not overlap
6. **GU stacking caveat:** GU/UG and UG/GU tandem stackings are destabilizing per NNDB

### 7.2 Base Pairing Invariants

1. Watson-Crick pairs: A-U, U-A, G-C, C-G only
2. Wobble pairs: G-U, U-G only
3. Non-pairing: A-A, A-G, A-C, U-U, U-C, C-C, G-G

### 7.3 Probability Invariants

1. Structure probability: 0 ≤ P ≤ 1
2. MFE structure: highest probability
3. Boltzmann probability: exp(-E/RT) / Z

---

## 8. Deviations and Assumptions

### 8.1 Resolved

| # | Item | Status | Resolution |
|---|------|--------|------------|
| D1 | Single-nucleotide bulge missing −RT·ln(states) degeneracy | 🔧 FIXED | Added `numStates` parameter to `CalculateBulgeLoopEnergy`; DP computes states from adjacent identical bases. Verified against NNDB Example 1: 3 C's → −0.616·ln(3) = −0.68 kcal/mol. |
| D2 | Dangling ends (model d2) not in multiloop WM | 🔧 FIXED | Added Options A2/A3/A4 to WM recurrence: helix with 5' dangle, 3' dangle, or both. Matches exterior loop (W) treatment. |

### 8.2 Open

| # | Item | Status | Rationale |
|---|------|--------|-----------|
| D3 | 1×2 internal loop (int21) lookup table | ⛔ BLOCKED | NNDB int21 table has 2,304 entries. Currently uses generic initiation + asymmetry + mismatch model. Data too large for inline static table; would require external data file. |
| D4 | 2×2 internal loop (int22) lookup table | ⛔ BLOCKED | NNDB int22 table has 36,864 entries. Same issue as int21. |
| D5 | No Zuker traceback for optimal structure | ⚠ ASSUMPTION | `CalculateMinimumFreeEnergy` computes correct MFE value via Zuker DP, but `PredictStructure` uses greedy stem-loop selection instead of DP traceback. MFE value is correct; predicted structure may differ from global optimum. |

### 8.3 Design Decisions (not deviations)

- **No pseudoknot prediction** — O(n³) DP inherently cannot find pseudoknots; correct by design.
- **Single sequence only** — Comparative prediction (covariance) is a different algorithm class.
- **Minimum loop size 3** — Correct per NNDB steric constraint; not a limitation.

---

## 9. Design Decisions

| Decision | Rationale |
|----------|----------|
| Turner 2004 (NNDB) parameters for stacking and hairpin energies | Standard reference, well-validated, exact NNDB values used |
| Minimum loop size of 3 | Biological steric constraint (Wikipedia, NNDB) |
| Standard ACGU bases only | No modified bases |
| Greedy non-overlapping stem-loop selection | Standard approach for simple prediction |
| Zuker-style DP for MFE | Physical kcal/mol values; uses Turner stacking, hairpin, internal, bulge, multiloop energies |
| Special hairpin loops as total replacements | Matches NNDB approach: experimental data supersedes model calculation |
| Terminal mismatch table (96 entries) | Full NNDB table for all 6 closing pairs × 16 mismatch combinations |
| Special GU closure (−2.2) | NNDB: applied when G-U closing pair preceded by two Gs on 5' side |
| Energy methods as public statics | Allow direct testing and external use before full DP integration |

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
- [x] WC stacking energies match NNDB Turner 2004 exact values (16 entries)
- [x] GU wobble stacking energies from NNDB (20 entries, incl. destabilizing)
- [x] Hairpin loop initiation energies match NNDB for all sizes 3–30
- [x] Special hairpin loops use NNDB total energies (tri/tetra/hexaloops)
- [x] All-C penalty formula matches NNDB (+1.5 for 3nt, 0.3n+1.6 for n>3)
- [x] UU/GA and GG mismatch bonuses from NNDB (−0.9, −0.8)
- [x] Extrapolation formula for loops >30 uses NNDB Jacobson–Stockmayer formula
- [x] No default stacking energy fallback (unknown stackings contribute 0)
- [x] Terminal AU/GU penalty (+0.45 per end) applied per NNDB
- [x] Terminal mismatch stacking table (96 entries) from NNDB tm-parameters.html
- [x] Special GU closure (−2.2) for G-U closing preceded by two Gs
- [x] Internal loop energy model: initiation + asymmetry + AU/GU closure + mismatch
- [x] Bulge loop energy model: n=1 with stacking continuation + degeneracy −RT·ln(states), n>1 with terminal penalties
- [x] Multibranch loop energy model: offset + asymmetry + helix count + strain
- [x] Coaxial stacking: flush (WC stacking table) and mismatch-mediated (tm + base + bonus)
- [x] Dangling end energy tables: 24 3'-dangles + 24 5'-dangles from NNDB
- [x] Dangling ends (model d2) integrated into Zuker DP exterior loop (W) and multiloop (WM)
- [x] 1×1 internal loop int11 lookup table (576 entries from NNDB, includes AU/GU penalties)
- [x] Zuker DP MFE matches manual Turner 2004 calculations (verified: GGGAAACCC=-1.12, GGGGAAAACCCC=-5.28)
- [x] Zuker DP integrates internal loop, bulge loop, and multibranch loop energy models
- [x] NNDB-exact validation tests added for all parameter categories
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
