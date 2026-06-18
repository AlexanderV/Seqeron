# Hairpin Loop and Stem Free-Energy Calculation (Turner 2004)

| Field | Value |
|-------|-------|
| Algorithm Group | RnaStructure |
| Test Unit ID | RNA-HAIRPIN-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

This unit computes the standard Gibbs free-energy change (ΔG°37, kcal/mol) of the two
elementary motifs of an RNA stem-loop: the **hairpin loop** (the terminal single-stranded
loop) and the **stem** (the double-stranded helix). It implements the Turner 2004
nearest-neighbor thermodynamic model as published by Mathews et al. (2004) [1] and tabulated
in the Nearest Neighbor Database (NNDB) [2][3]. The calculation is deterministic and
specification-driven: given a loop sequence with its closing base pair, or a list of stacked
base pairs, it returns the exact sum of the cited sequence-dependent parameters. It is used as
a building block for minimum-free-energy folding and stem-loop scoring.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

RNA folds back on itself to form helices (stems) closed by loops. The dominant model for the
stability of these motifs is the nearest-neighbor model, in which the free energy of a
structure is the sum of local sequence-dependent contributions: base-pair stacking inside
helices and length- and sequence-dependent initiation/mismatch terms for loops [1][3]. All
parameters are reported at 37 °C, 1 M NaCl [3].

### 2.2 Core Model

**Hairpin loop (n > 3 unpaired nt)** [3]:

```
ΔG°37(hairpin) = ΔG°37 initiation(n)
               + ΔG°37(terminal mismatch)
               + ΔG°37(UU or GA first mismatch)
               + ΔG°37(GG first mismatch)
               + ΔG°37(special GU closure)
               + ΔG°37 penalty(all-C loops)
```

**Hairpin loop (n = 3 unpaired nt)** [3] — no sequence-dependent first-mismatch term:

```
ΔG°37(hairpin, 3 nt) = ΔG°37 initiation(3) + ΔG°37 penalty(all-C loops)
```

The **terminal mismatch** is the sequence-dependent stacking of the first/last loop bases on
the closing pair [3]. The **special GU closure** term (−2.2) is applied only when the closing
pair is G(5')–U(3') (not U–G) preceded by two Gs [3]. The **all-C penalty** is +1.5 for a
3-nt all-C loop and `A·n + B` (A=+0.3, B=+1.6) for longer all-C loops [3]. Certain tri-,
tetra- and hexaloops are poorly fit by the model and instead take an experimentally measured
**total** that overrides the additive formula [3].

For n > 9 the initiation is extrapolated by `ΔG°37(n) = ΔG°37(9) + 1.75·R·T·ln(n/9)` [3].

**Stem (helix)** [4]: the helix free energy is the sum of nearest-neighbor stacking terms plus
a per-end penalty: for a helix of `P` uninterrupted base pairs there are `P−1` stacks; an
extra +0.45 ("per AU end") is added once for each helix terminus that closes with an A-U/U-A
or G-U/U-G pair [4]. (Intermolecular initiation and the self-complementary symmetry correction
from the full duplex equation are not part of a unimolecular stem and are excluded.)

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Nearest-neighbor independence: each motif's energy is a sum of local terms. | Long-range/tertiary effects (e.g. coaxial stacking, ions) are not captured; predicted ΔG drifts from experiment. |
| ASM-02 | Standard conditions 37 °C, 1 M NaCl. | At other temperatures/salts the tabulated parameters no longer apply. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Deterministic in (loop sequence, closing pair, special-GU flag). | Pure table lookups + arithmetic [3]. |
| INV-02 | Loops < 3 nt return a prohibitive energy, never a normal low value. | NN rules prohibit loops < 3 nt [3]; implementation returns 100.0 to exclude them. |
| INV-03 | A special tri/tetra/hexaloop returns exactly its tabulated total. | Special table overrides the additive model [3]. |
| INV-04 | Stem energy of 0 base pairs is 0; a stem of P pairs sums P−1 stacks. | P−1 stacks per helix [4]. |
| INV-05 | The −2.2 special-GU bonus applies only to a G(5')-U(3') closing pair. | "a GU closing pair (not UG)" [3]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `loopSequence` | string | required | Unpaired loop bases, 5'→3' | RNA alphabet; case-insensitive (uppercased) |
| `closingBase5` | char | required | 5' base of the closing pair | A/C/G/U |
| `closingBase3` | char | required | 3' base of the closing pair | A/C/G/U |
| `specialGUClosure` | bool | false | Whether the G-U closing pair is preceded by two Gs | applied only when closing pair is G-U |
| `sequence` | string | required (stem) | Sequence context (unused for energy; kept for API parity) | RNA |
| `basePairs` | IReadOnlyList\<BasePair\> | required (stem) | Stacked pairs, outer→inner | each pair has Base1 (5'), Base2 (3'), Type |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | double | ΔG°37 in kcal/mol, rounded to 2 decimals. |

### 3.3 Preconditions and Validation

Bases are uppercased; T is not auto-converted to U (RNA expected). For `CalculateHairpinLoopEnergy`,
loops shorter than 3 nt return a prohibitive 100.0 (INV-02); unknown terminal-mismatch keys
contribute 0. For `CalculateStemEnergy`, an empty `basePairs` list returns 0 (INV-04); unknown
stacking keys contribute 0. Indexing within `basePairs` is outer (helix end) → inner (loop end).

## 4. Algorithm

### 4.1 High-Level Steps

**Hairpin:**
1. Uppercase the loop; build the special-loop key `closing5' + loop + closing3'`. If found, return its total (INV-03).
2. Look up initiation(n); for n>30 use the log extrapolation; for n<3 return 100.0 (INV-02).
3. For n ≥ 4: add terminal mismatch, UU/GA bonus, GG bonus, and (if applicable) special-GU closure.
4. Add the all-C penalty (C3 for 3-nt, linear `A·n+B` for longer).
5. Round to 2 decimals.

**Stem:**
1. If no pairs, return 0.
2. Sum nearest-neighbor stacking over consecutive pairs (P−1 terms); apply the GGUC/CUGG 3-stack special context where it occurs.
3. Add +0.45 for each helix end (first/last pair) that closes with A-U/U-A or G-U/U-G.
4. Round to 2 decimals.

### 4.2 Decision Rules, Scoring, Reference Tables

All numeric parameters are the Turner 2004 NNDB tables [3][4]:
hairpin initiation (loop.txt); terminal mismatch (tstack.txt); first-mismatch bonuses and all-C
penalty (hairpin-mismatch-parameters.html); special tri/tetra/hexaloops (triloop/tloop/hexaloop.txt);
Watson-Crick/GU stacking and the +0.45 per-AU/GU-end penalty (wc-parameters.html, gu-parameters.html).
No search/matching is performed, so the repository suffix tree is **not applicable** to this unit.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateHairpinLoopEnergy` | O(n) | O(1) | n = loop length; table lookups + an all-C scan |
| `CalculateStemEnergy` | O(P) | O(1) | P = number of base pairs |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.CalculateHairpinLoopEnergy(string, char, char, bool)`: ΔG°37 of a hairpin loop.
- `RnaSecondaryStructure.CalculateStemEnergy(string, IReadOnlyList<BasePair>)`: ΔG°37 of a helix.

### 5.2 Current Behavior

Parameters are stored as static dictionaries / arrays initialized once. Loops < 3 nt yield a
deliberately prohibitive 100.0 so a downstream optimizer never selects them. The stem routine
recognises the NNDB GGUC/CUGG three-stack context (−4.12 total) in addition to the per-step
stacking. Results are rounded to 2 decimals. No suffix-tree search is used (§4.2).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Hairpin additive model for n ≥ 4 and the reduced n = 3 model (no first-mismatch term) [3].
- Hairpin initiation table and n>9/n>30 log extrapolation [3].
- Terminal mismatch, UU/GA bonus (−0.9), GG bonus (−0.8), special-GU closure (−2.2), all-C penalties (C3 +1.5; A=0.3,B=1.6) [3].
- Special tri/tetra/hexaloop totals overriding the model [3].
- Nearest-neighbor stacking and the per-AU/GU-end penalty (+0.45) for stems [4].

**Intentionally simplified:**

- Stem energy excludes intermolecular initiation and the self-complementary symmetry correction; **consequence:** the value is the unimolecular helix contribution only, which is the correct component inside a folded structure [4].

**Not implemented:**

- Coaxial stacking, dangling ends on the stem ends, and multibranch/exterior loop terms; **users should rely on:** `CalculateMinimumFreeEnergy` (RNA-MFE-001) for whole-structure energies.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Rounding to 2 decimals | Assumption | NNDB prints 1-decimal totals; tests assert exact parameter sums | accepted | display-precision only |
| 2 | Loops < 3 nt → 100.0 | Assumption | source gives no value; a sentinel is needed | accepted | INV-02 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Loop < 3 nt | Prohibitive energy (≥100) | NN rules prohibit such loops [3] |
| 3-nt loop | initiation(3) only (+all-C penalty if all-C) | no first-mismatch term for 3-nt loops [3] |
| Special tri/tetra/hexaloop | tabulated total overrides model | experimentally fit [3] |
| Empty base-pair list | 0 | P−1 stacks with P=0 [4] |
| U-G closing + specialGUClosure flag | no −2.2 applied | bonus is G-U only [3] |

### 6.2 Limitations

Restricted to the Turner 2004 parameter set at 37 °C / 1 M NaCl (ASM-02). Hairpin loop term
omits the salt/temperature corrections of full thermodynamic packages. Stem term is the helix
contribution only; it does not assemble a complete structure energy.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// NNDB Hairpin Example 1: closing A-U, 6-nt loop A..A
double loop = RnaSecondaryStructure.CalculateHairpinLoopEnergy("AAAAAA", 'A', 'U'); // +4.6
```

**Numerical walk-through (NNDB Example 1):** loop = initiation(6) +5.4 + terminal mismatch
(A·A on A-U) −0.8 = **+4.6**; helix (pairs C-G,A-U,C-G,A-U) = 3 stacks (−2.11,−2.24,−2.11) +
one AU-end penalty +0.45 = **−6.01**; total stem-loop = **−1.4 kcal/mol** [1][3][4].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_HairpinEnergy_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_HairpinEnergy_Tests.cs) — covers INV-01..INV-05.
- Evidence: [RNA-HAIRPIN-001-Evidence.md](../../../docs/Evidence/RNA-HAIRPIN-001-Evidence.md)
- Related algorithms: [RNA_Base_Pairing](../RnaStructure/RNA_Base_Pairing.md)

## 8. References

1. Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH. 2004. Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. Proc. Natl. Acad. Sci. USA 101:7287-7292. https://doi.org/10.1073/pnas.0401799101
2. Turner DH, Mathews DH. 2010. NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure. Nucleic Acids Res. 38:D280-D282. https://doi.org/10.1093/nar/gkp892
3. NNDB. Turner 2004 Hairpin Loop Parameters. https://rna.urmc.rochester.edu/NNDB/turner04/hairpin.html (accessed 2026-06-14)
4. NNDB. Turner 2004 Watson-Crick Helix Parameters. https://rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html (accessed 2026-06-14)
