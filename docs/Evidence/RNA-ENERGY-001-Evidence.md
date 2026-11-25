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
3. **Loop initiation is not monotonic for small loops** (§3.4: n = 6 and 8 are below n = 4); only the n > 9 extrapolation grows with n
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
| Unknown base (N, IUPAC, gap) in a mismatch / dangle position | Maximum over A/C/G/U (ViennaRNA 2.7.2 'N' entries, §8.1) |
| DNA T, lower case | Read as U / upper case (ViennaRNA `vrna_nucleotide_encode`) |
| Non-canonical closing pair in a mismatch / dangle / coaxial term | 0 (outside the NN model) |

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
**Unknown unpaired bases** (mismatch / dangle positions): maximum over A/C/G/U, see §8.1.

---

## 7. References

1. Turner, D.H. (2004). Thermodynamics of RNA. In: RNA World, 3rd Ed.
2. Mathews, D.H., et al. (2004). PNAS 101(19):7287-7292.
3. Xia, T., et al. (1998). Biochemistry 37(42):14719-35.
4. SantaLucia, J. (1998). PNAS 95(4):1460-5.
5. NNDB: https://rna.urmc.rochester.edu/NNDB/
6. Heus, H.A. & Pardi, A. (1991). Science 253(5016):191-194.
7. Wikipedia - Nucleic acid thermodynamics (orientation only; not used for any value)
8. Wikipedia - Nucleic acid secondary structure (orientation only; not used for any value)
9. Mathews, D.H. & Turner, D.H. (2002). Biochemistry 41:869-880 (multibranch efn2 parameters).
10. Mathews, D.H., Sabina, J., Zuker, M. & Turner, D.H. (1999). J. Mol. Biol. 288:911-940 (coaxial stacking).
11. ViennaRNA 2.7.2 source (PyPI sdist): `params/default.c`, `eval/eval_exterior.c`, `eval/eval_internal.c`; Python `RNA` 2.7.2.

---

## 8. B12 review 2026-11 — remaining energy terms (Stage A + Stage B)

Sources opened: NNDB Turner 2004 pages (rna.urmc.rochester.edu/NNDB/rna_2004/): dangling ends, terminal mismatches,
bulge loops, coaxial stacking (§8.1–8.4) + `parameter_tables/turner_2004_rnastructure/turner_2004_coax.html`,
multibranch loops (§9.1–9.3) + `turner_2004_multi_branch` table (a 9.25, b 0.91, c −0.63, strain 3.14);
`turner_2004_terminal_mismatch_dg.txt`; ViennaRNA 2.7.2 sdist `params/default.c` (`mismatchExt37`, `mismatchM37`,
`dangle5_37`, `dangle3_37`, `stack37`, `bulge37`, `ML_closing37 = 930`, `ML_intern37 = −90`, `ML_BASE37 = 0`),
`misc/rna_turner2004.par` (`ML_params 0 0 930 3000 −90 −220`); Ward et al. 2017 NAR 45:8541 (linear model
9.3 / 0 / −0.9 of Mathews et al. 2004). Oracle: Python `RNA` 2.7.2, default file and the NNDB-precision file
(RNA-STRUCT-001 Evidence §12.2, file b). Scripts: `en/chk1.py` … `chk4.py` in the B12 scratchpad.

### 8.1 Unknown bases, DNA T, case in terminal mismatches / dangles (F13)

ViennaRNA's default file fills every 'N' row/column (index 0) of `mismatchExt37`, `mismatchM37`, `dangle5_37`,
`dangle3_37` with the maximum over A/C/G/U (e.g. `mismatchExt37[CG][N][A] = max(−110, −150, −130, −150) = −110`).
`GetTerminalMismatchEnergy(c5, c3, x, y)` = `mismatchExt37[pair(c3, c5)][y][x]`;
`GetDanglingEndEnergy(c5, c3, d, 3')` = `dangle3_37[pair(c3, c5)][d]`, 5' = `dangle5_37[pair(c3, c5)][d]`.

| Call | Before | After | ViennaRNA |
|------|--------|-------|-----------|
| `GetTerminalMismatchEnergy('C','G','N','N')` | 0 | −0.8 | `mismatchExt37[GC][N][N]` −80 |
| `GetDanglingEndEnergy('C','G','N', 3')` | 0 | −0.8 | `dangle3_37[GC][N]` −80 |
| `GetDanglingEndEnergy('A','T','A', 3')` | 0 | −0.8 | `dangle3_37[UA][A]` −80 |
| `CalculateMinimumFreeEnergy("NGGGAAACCCN")` | −1.12 | −1.92 | `RNA.fold` −1.92 (NNDB-precision file; default file −2.00) |

Systematic check: 6 canonical pairs × 5 × 5 (incl. N) for both mismatch tables, 6 × 5 dangles × 2, each with
lower-case / T / IUPAC (R, Y, X, '-') variants: 1,950 calls, 0 mismatches. Engine: 404 random RNAs (15–80 nt,
1/13 N): our MFE structure's energy = `eval_structure` (NNDB-precision file) + NNDB-only terms 404/404; our MFE =
`RNA.fold` MFE for every sequence whose optimum has no NNDB-only motif. A non-canonical closing pair keeps 0 (ViennaRNA
scores it only through its "non-standard" pair type 7, not part of the NNDB model).

### 8.2 Bulge special C without the strand; DNA T; size (F14)

NNDB bulge loops: "−0.9 kcal/mol for a bulged C adjacent to at least one paired C". The strand of a 1-nt bulge is not
derivable from the six arguments (the bulged base sits between pair5_base1/pair3_base1 or between
pair3_base2/pair5_base2). Before: `bulgeOn5PrimeStrand = null` accepted a C anywhere among the four closing bases
(`(1,'C','G','C','A','U')` = 0.55, bonus from the C on the opposite strand). After: null applies the bonus only when
it holds on both strands (least stabilising reading, the ViennaRNA unknown convention) → 1.45 = ViennaRNA
`eval_int_loop` (no special C); with the strand given the exact NNDB value (true → 1.45, false → 0.55). MCP
`bulge_loop_energy` gained the optional `bulgeOn5PrimeStrand` parameter. Also: T/lower case read as U
(`(1,'a','G','C','A','T',side 5')` 3.8 → 1.45; `(2,'A','A','T','G','C')` 2.8 → 3.25); `bulgeSize < 1` now throws
(was 3.8 / nonsense). Check vs `eval_int_loop` + NNDB special C: 6×6 pairs × n 1–30, 31, 35, 40 × A/C/G/U × both
strands × side given/null = 19,008 calls, 0 mismatches (n > 30: NNDB extrapolation from n = 6, ≤ 0.04 from
ViennaRNA's `lxc·ln(n/30)`).

### 8.3 Multibranch loop: two sourced models (F15)

Before: `CalculateMultibranchLoopEnergy(h, u, strain, stack)` = 9.25 + 0.91·(u/h) − 0.63·h — "u/h" is not the NNDB
average asymmetry (mean |5' − 3' unpaired| per helix, capped at 2.0), and the engines use a different model
(9.3 / 0 / −0.9), so the public value disagreed with both. After:
- `(int numHelices, int numUnpaired, …)` = linear Turner 2004 DP model 9.3 + 0·u − 0.9·h (+ stacking) — identical
  to the engine and to ViennaRNA: `eval_structure_verbose` "Multi loop" 660 (GC-closed 3-way, d0) = 6.6; 570
  (4-way, 0 or 11 unpaired) = 5.7; d2 250 = 6.6 + three terminal mismatches (−1.1 − 1.5 − 1.5). `hasStrain = true`
  throws (strain exists only in efn2); `h < 3` / `u < 0` throw.
- new `(IReadOnlyList<int> unpairedBetweenHelices, double stackingEnergy)` = NNDB efn2: 9.25 + 0.91·min(2, mean
  asymmetry) − 0.63·h + 3.14 (3-way, < 2 unpaired). NNDB §9.3: 4 helices, asymmetry capped at 2 → 8.55 (NNDB
  prints 8.6), total with stacking −5.9 → 2.65 (NNDB 2.7); gaps (2,2,2) → 7.36; (3,0,0) → 9.18; (0,0,0) → 10.50;
  (1,0,0) + stacking −2.0 → 9.11. MCP `multibranch_loop_energy` gained optional `unpairedBetweenHelices`.

### 8.4 Coaxial stacking (F16)

NNDB §8.3 / coax table: mismatch-mediated = terminal mismatch (continuous backbone) − 2.1, plus −0.4 if the
"mismatch" could form a Watson–Crick pair or −0.2 if a G·U pair. Before, the bonus was chosen by the CLOSING pair
(always −0.4, or −0.2 for a G·U closing pair): 76/96 pair × mismatch combinations wrong, e.g. NNDB §8.4 example
A·U + GG mismatch −3.3 → **−2.9** (NNDB −2.9); G·C + AA −3.6 → −3.2. Flush = helical stack (NNDB §8.4 AU→CG −2.24)
was correct for A/C/G/U; T/lower case now read as RNA (`('a','t','c','g')` 0 → −2.24). Neither function is used by the
engines (ViennaRNA's default d2 model has no coaxial stacking); they are public/MCP evaluators. Check: 6 pairs × 16
mismatches (+ case/T/N) and 36 flush combinations vs NNDB / `eval_int_loop` stacks: 256 calls, 0 mismatches.

### 8.5 Turner2004Parameters (F17)

`StackingEnergies` (36), `HairpinLoopInitiation` (3–30), `TerminalMismatchEnergies` (96 = NNDB txt, 96/96),
`TerminalAuGuPenalty` 0.45 re-verified. `IsTerminalPenaltyPair` was case/T-blind (`('a','t')` false); now
case-insensitive with T as U. The `TerminalMismatchEnergies` doc claimed use for internal loops (they use
`Turner2004InteriorLoopTables`) — corrected.

### 8.6 Structure evaluation `CalculateFreeEnergy` (F18)

The checklist canonical `CalculateFreeEnergy(structure)` did not exist (only the internal `EvaluateStructureEnergy`).
New public `CalculateFreeEnergy(rnaSequence, dotBracket)` (only `( ) .`; unbalanced / wrong length / non-canonical
pair throw) delegates to the shared evaluator; MCP `rna_structure_energy`. ViennaRNA `eval_structure` (d2,
NNDB-precision file): GGGAAACCC −1.12, GGGAUAAAUCCC −3.42 (NNDB hairpin example 1), GGGAAACCCAGGGAAACCC
`(((...))).(((...)))` −4.14, GAGGGAAAACCCAGGGAAAACCCAC 3-way −1.54, NGGGAAACCCN −1.92 — all identical.

### 8.7 Stage verdicts

Stage A: PASS-WITH-NOTES — §4.4 invariant 3 ("monotonic loop penalty") contradicted §3.4 (corrected); "invalid bases:
undefined" replaced by the ViennaRNA convention; Wikipedia rows kept for orientation only. Stage B: FAIL (fixed) —
F13–F17 above.
