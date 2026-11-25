# RNA Free Energy Calculation

| Field | Value |
|-------|-------|
| Algorithm Group | RNA Structure |
| Test Unit ID | RNA-ENERGY-001 |
| Related Projects | N/A |
| Implementation Status | Complete (Turner 2004 nearest-neighbour terms; NNDB values, ViennaRNA 2.7.2 conventions) |
| Last Reviewed | 2026-11-25 (B12 review, RNA-ENERGY-001 F13–F17) |

## 1. Overview

RNA free-energy calculation estimates the thermodynamic stability of RNA secondary-structure motifs and whole
sequences under the Turner 2004 nearest-neighbour model (NNDB). `RnaSecondaryStructure` exposes one public
evaluator per NNDB term — helix (stem) stacking, hairpin, interior, bulge and multibranch loops, dangling ends,
terminal mismatches, flush and mismatch-mediated coaxial stacking — and the MFE dynamic program that composes them
(RNA-MFE-001 / RNA-STRUCT-001). The folding engines call the same evaluators (one energy path), so a term's public
value is the value the MFE, traceback and partition function use. Shared tables live in
`Seqeron.Genomics.Core.Turner2004Parameters` (also used by `MiRnaAnalyzer`) and the generated
`Turner2004InteriorLoopTables`.

## 2. Scientific / Formal Basis

### 2.1 Core Model

ΔG°37(structure) = Σ stacks + Σ hairpin + Σ interior/bulge + Σ multibranch + exterior terms (NNDB Turner 2004;
Mathews et al. 2004). Every term is in kcal/mol at 37 °C; more negative = more stable.

### 2.2 Terms and sources

| Term | Public method | Rule | Source |
|------|---------------|------|--------|
| Helix | `CalculateStemEnergy` | Σ WC/GU stacks (GGUC/CUGG 3-stack −4.12) + 0.45 per AU/GU helix end | NNDB WC + GU pages |
| Hairpin | `CalculateHairpinLoopEnergy` | special tri/tetra/hexaloop totals, else initiation(n) (+1.75RT ln(n/9) for n > 30) + terminal mismatch + UU/GA −0.9, GG −0.8, special GU −2.2, all-C | NNDB hairpin page; ViennaRNA `eval_hairpin.c` |
| Interior | `CalculateInternalLoopEnergy` | int11/int21/int22 tables; 1×n, 2×3, generic rules; Ninio min(3.0, 0.6·\|n1−n2\|) | NNDB internal loops; ViennaRNA `vrna_E_internal` |
| Bulge | `CalculateBulgeLoopEnergy` | n = 1: init 3.8 + stack across the bulge, −0.9 if the bulged C is adjacent to a paired C on its strand, −RT ln(states); n > 1: init(n) + AU/GU closure penalties; n > 30: init(6) + 1.75RT ln(n/6) | NNDB bulge loops |
| Terminal mismatch | `GetTerminalMismatchEnergy` | NNDB terminal-mismatch table (= ViennaRNA `mismatchExt37` / `mismatchM37`) | NNDB terminal mismatches |
| Dangling end | `GetDanglingEndEnergy` | NNDB 3'/5' dangle tables (= ViennaRNA `dangle3_37` / `dangle5_37`) | NNDB dangling ends |
| Multibranch (engine) | `CalculateMultibranchLoopEnergy(numHelices, numUnpaired, …)` | linear DP model 9.3 + 0·unpaired − 0.9·helices (+ given stacking) | Mathews et al. 2004; ViennaRNA `default.c` `ML_closing37/ML_BASE37/ML_intern37` |
| Multibranch (efn2) | `CalculateMultibranchLoopEnergy(IReadOnlyList<int> gaps, stacking)` | 9.25 + 0.91·min(2, mean \|5'−3' unpaired\|) − 0.63·helices + 3.14 strain (3-way, < 2 unpaired) + stacking | NNDB multibranch §9.1, table turner_2004_multi_branch |
| Flush coaxial stack | `CalculateFlushCoaxialStacking` | WC/GU helical stack as if no backbone break | NNDB coaxial stacking §8.2 |
| Mismatch-mediated coaxial stack | `CalculateMismatchCoaxialStacking` | terminal mismatch − 2.1, − 0.4 / − 0.2 only if the mismatch could form a WC / GU pair | NNDB coaxial stacking §8.3, table turner_2004_coax |

The two multibranch models are both sourced: the folding engines use the linear model (an O(n³) recursion cannot
evaluate the efn2 average asymmetry), and the efn2 overload evaluates a given loop as NNDB's efn2 does.

### 2.3 Conventions

- Bases are case-insensitive; DNA T is read as U (ViennaRNA `vrna_nucleotide_encode`).
- An unpaired base other than A/C/G/U (N, IUPAC, gap) in a mismatch / dangle position is ViennaRNA's unknown base
  'N', whose Turner 2004 table entries are the maximum (least stabilising) over A/C/G/U; the evaluators take that
  maximum (exact for every ViennaRNA 'N' row/column).
- A non-canonical closing pair (not WC/GU) is outside the nearest-neighbour model: terminal-mismatch, dangle and
  coaxial terms return 0 (ViennaRNA only scores such pairs through its non-default "non-standard" type).
- Helix-end orientation for mismatches and dangles: `5'-closingBase5 x-3' / 3'-closingBase3 y-5'`; a 3' dangle is
  x (3' of closingBase5), a 5' dangle is y (5' of closingBase3).
- Bulge special C with the strand unknown (`bulgeOn5PrimeStrand = null`): the bonus is applied only when it holds
  on both strands (least stabilising reading, the same "maximum over the unknown" convention).

### 2.4 Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Hairpins shorter than 3 nt are prohibited | NNDB; generic loops < 3 return 100 |
| INV-02 | `CalculateStemEnergy([])` = 0 | no stack |
| INV-03 | Public helpers round to 2 decimals | `Math.Round(…, 2)` |
| INV-04 | A term with an unknown base never exceeds the stability of any concrete base choice | maximum over substitutions |

## 3. Contract

| Method | Validation |
|--------|-----------|
| `CalculateBulgeLoopEnergy` | `bulgeSize < 1` → `ArgumentOutOfRangeException` |
| `CalculateInternalLoopEnergy` | negative `n1`/`n2` → `ArgumentOutOfRangeException` |
| `CalculateMultibranchLoopEnergy(int,int,bool,double)` | `numHelices < 3`, `numUnpaired < 0` → `ArgumentOutOfRangeException`; `hasStrain = true` → `ArgumentException` (strain is an efn2 term) |
| `CalculateMultibranchLoopEnergy(IReadOnlyList<int>,double)` | null → `ArgumentNullException`; < 3 helices or a negative count → `ArgumentException` |
| `CalculateMinimumFreeEnergy` | empty/short → 0; `minLoopSize < 3` clamped to 3 |

## 4. Cross-verification (B12, 2026-11)

- `GetTerminalMismatchEnergy` / `GetDanglingEndEnergy` vs ViennaRNA 2.7.2 `default.c` `mismatchExt37`, `mismatchM37`,
  `dangle5_37`, `dangle3_37`: all 6 pairs × 5×5 bases incl. N, with lower case / T / IUPAC variants — 1,950 checks,
  0 mismatches. Terminal-mismatch table vs NNDB `turner_2004_terminal_mismatch_dg.txt`: 96/96 identical.
- `CalculateBulgeLoopEnergy` vs ViennaRNA `eval_int_loop` (NNDB-precision parameter file) + NNDB special C:
  6×6 closing pairs × sizes 1–30, 31, 35, 40 × A/C/G/U × both strands × side given / null — 19,008 checks, 0 mismatches.
- `CalculateMismatchCoaxialStacking` vs the NNDB rule on the NNDB tm file (6 pairs × 16 mismatches + case/T/N) and
  `CalculateFlushCoaxialStacking` vs the helical stack (36 pair combinations): 256 checks, 0 mismatches.
- Engine with unknown bases: 404 random RNAs (13 % N) — our MFE structure's energy = `RNA.eval_structure`
  (NNDB-precision file) + NNDB-only terms for all 404; MFE = `RNA.fold` for every sequence without NNDB-only motifs.
- Multibranch linear model vs `eval_structure_verbose` "Multi loop": 660 (3-way, d0), 570 (4-way, 0 or 11 unpaired),
  250 (3-way, d2) — identical.

Details and numbers: [RNA-ENERGY-001 Evidence](../../Evidence/RNA-ENERGY-001-Evidence.md) §8.

## 5. Implementation Notes

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs),
[Turner2004Parameters.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/Turner2004Parameters.cs).

- The unknown-base maximum is one private helper (`MaxOverUnknownBases`) shared by the hairpin first mismatch, the
  terminal mismatch, the dangles and the mismatch coaxial stack.
- The MFE / traceback / partition-function exterior and multiloop stem terms call `GetTerminalMismatchEnergy` /
  `GetDanglingEndEnergy`; interior and bulge terms call `CalculateInternalLoopEnergy` / `CalculateBulgeLoopEnergy`
  (the engine always passes the bulge strand).
- MCP tools (`stem_energy`, `hairpin_loop_energy`, `internal_loop_energy`, `bulge_loop_energy`,
  `multibranch_loop_energy`, `terminal_mismatch_energy`, `dangling_end_energy`, `flush_coaxial_stacking`,
  `mismatch_coaxial_stacking`, `minimum_free_energy`) delegate to these methods.

## 6. Edge Cases and Limitations

| Case | Behaviour | Rationale |
|------|-----------|-----------|
| Unknown base in a mismatch / dangle | maximum over A/C/G/U | ViennaRNA 'N' entries |
| Non-canonical closing pair | 0 for mismatch / dangle / coaxial terms | outside the NN model |
| Hairpin / bulge > 30 nt | NNDB logarithmic extrapolation (ViennaRNA truncates `lxc·ln(n/30)`; ≤ 0.04 kcal/mol apart) | NNDB formula kept |
| Multibranch strain without per-gap counts | rejected | strain is defined only in the efn2 model |

Temperature dependence (ΔH) and pseudoknot energies are not part of this unit.

## 7. References

1. NNDB, Turner 2004 RNA parameters: https://rna.urmc.rochester.edu/NNDB/ (pages 1–10 and parameter tables).
2. Mathews, D.H., Disney, M.D., Childs, J.L., Schroeder, S.J., Zuker, M., Turner, D.H. (2004). PNAS 101:7287–7292.
3. Mathews, D.H., Turner, D.H. (2002). Biochemistry 41:869–880 (multibranch efn2 parameters).
4. Mathews, D.H., Sabina, J., Zuker, M., Turner, D.H. (1999). J. Mol. Biol. 288:911–940 (coaxial stacking).
5. Lorenz, R. et al. (2011). ViennaRNA Package 2.0. Algorithms Mol. Biol. 6:26; ViennaRNA 2.7.2 source (`params/default.c`, `eval/*.c`).
6. Xia, T. et al. (1998). Biochemistry 37:14719–14735.
