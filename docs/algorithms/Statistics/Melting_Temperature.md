# Melting Temperature (Wallace / Marmur-Doty / Nearest-Neighbor)

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics |
| Test Unit ID | SEQ-TM-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Infrastructure |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Estimates the melting temperature (Tm) of a DNA oligonucleotide — the temperature at
which half the duplex is dissociated. Two complementary routines are provided:
`CalculateMeltingTemperature` gives a fast estimate via the **Wallace rule** for short
oligos (< 14 nt) or the **Marmur-Doty GC formula** for longer ones [1][2], and
`CalculateThermodynamics` gives the rigorous **nearest-neighbor (NN)** Tm of
Allawi & SantaLucia (1997) / SantaLucia (1998) [3][4]. The Wallace/GC formulas are
empirical rules of thumb; the NN model is specification-driven and exact for its
parameter table. SEQ-TM-001 documents the *melting-temperature* view of the same two
methods delivered under SEQ-THERMO-001 (see §5.4 and §7.3).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Tm depends on sequence length, base composition (G·C pairs are more stable than A·T),
and ionic conditions. Short-oligo rules approximate Tm from base counts; the NN model
sums experimentally measured stacking free energies over each adjacent base-pair step
and converts ΔH°/ΔS° to Tm at a given strand and salt concentration [3][4].

### 2.2 Core Model

**Wallace rule** (short oligos, rule of thumb for 14–20 nt) [1]:

> Tm = 2·(A + T) + 4·(G + C)

**Marmur-Doty GC formula** (longer oligos) [2]:

> Tm = 64.9 + 41·(G + C − 16.4) / N,  with N = sequence length

**Nearest-neighbor Tm** [3][4]:

> Tm = (1000·ΔH°) / (ΔS° + R·ln(C_T / x)) − 273.15

where R = 1.987 cal/(°C·mol), x = 4 for non-self-complementary duplexes (x = 1 if
self-complementary), C_T = total strand concentration, ΔH° and ΔS° are summed from the
unified NN parameter table plus helix-initiation terms, and ΔS° carries the salt
correction ΔS°(salt) = ΔS°(1 M) + 0.368·(N−1)·ln[Na+] [4]. (The full ΔH°/ΔS°/ΔG°
derivation is documented in [DNA_Thermodynamics.md](DNA_Thermodynamics.md).)

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Two-state (all-or-none) duplex melting at fixed buffer/temperature | NN Tm mis-estimates partially melted or structured sequences |
| ASM-02 | Non-self-complementary strands in equal amount (x = 4) | Self-complementary hairpins need x = 1; Tm shifts |
| ASM-03 | Wallace/GC rules valid only in their length regime (< 14 nt / longer) | Outside the regime the empirical estimate degrades [1] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Wallace Tm = 2·(A+T) + 4·(G+C) | Definition [1] |
| INV-02 | Marmur-Doty Tm = 64.9 + 41·(GC − 16.4)/N | Definition [2] |
| INV-03 | NN Tm = (1000·ΔH°)/(ΔS° + R·ln(C_T/4)) − 273.15, R = 1.987 | Definition [3][4] |
| INV-04 | Higher [Na+] raises NN Tm | Salt term 0.368·(N−1)·ln[Na+] increases ΔS magnitude [4] |
| INV-05 | Empty / length-1 input → 0 (NN undefined for < 2 nt) | NN sums over dinucleotides [3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| dnaSequence | string | required | DNA sequence (5'→3') | case-insensitive; A/C/G/T |
| useWallaceRule | bool | true | Wallace if length < 14, else Marmur-Doty | — (CalculateMeltingTemperature) |
| naConcentration | double | 0.05 | Na+ concentration, mol/L | > 0 (CalculateThermodynamics) |
| primerConcentration | double | 2.5e-7 | Total strand conc. C_T, mol/L (divided by F = 4) | > 0 (CalculateThermodynamics) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (CalculateMeltingTemperature) | double | Tm in °C (Wallace or Marmur-Doty) |
| MeltingTemperature | double | NN Tm in °C (rounded to 1 dp) |
| DeltaH / DeltaS / DeltaG | double | ΔH° (kcal/mol), ΔS° (cal/(mol·K)), ΔG°37 (kcal/mol) |

### 3.3 Preconditions and Validation

Input is upper-cased (case-insensitive). `CalculateMeltingTemperature` returns 0 for
null/empty. `CalculateThermodynamics` returns an all-zero tuple for input shorter than
2 nt (the NN model is undefined without a dinucleotide step). No exceptions are thrown
for these guarded inputs.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count A/T/G/C (and length N).
2. **CalculateMeltingTemperature:** if `useWallaceRule` and N < 14, apply the Wallace
   rule; otherwise apply Marmur-Doty.
3. **CalculateThermodynamics:** add helix-initiation terms at both termini, sum NN
   ΔH°/ΔS° over each dinucleotide, apply the salt correction to ΔS°, then evaluate the
   NN Tm equation and ΔG°37.

### 4.2 Decision Rules, Scoring, Reference Tables

Wallace contributions (2 per A·T, 4 per G·C) and the Marmur-Doty constants
(64.9, 41, 16.4) live in `ThermoConstants`; the NN ΔH°/ΔS° table and initiation terms
(Allawi & SantaLucia 1997, Table 1) live in `SequenceStatistics`. Each constant is
source-cited in code.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateMeltingTemperature | O(n) | O(1) | base counting |
| CalculateThermodynamics | O(n) | O(1) | single pass over dinucleotides |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs), [ThermoConstants.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/ThermoConstants.cs)

- `SequenceStatistics.CalculateMeltingTemperature(string, bool)`: Wallace / Marmur-Doty estimate.
- `SequenceStatistics.CalculateThermodynamics(string, double, double)`: NN ΔH°/ΔS°/ΔG°/Tm.
- `ThermoConstants.CalculateWallaceTm`, `CalculateMarmurDotyTm`: shared formula constants.

### 5.2 Current Behavior

`CalculateMeltingTemperature` delegates the formula constants to `ThermoConstants`
(the same constants used by `PrimerDesigner.CalculateMeltingTemperature` in the MolTools
project — the formula is consolidated in one place). The repository default total strand
concentration for the NN model is C_T = 250 nM (vs Biopython's 50 nM); the formula is
identical and the default is an explicit parameter (see §5.4). No substring search /
matching is involved, so the repository suffix tree is **not applicable** to this unit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Wallace rule 2·(A+T) + 4·(G+C) [1].
- Marmur-Doty GC formula 64.9 + 41·(GC − 16.4)/N [2].
- NN Tm equation with R = 1.987, x = 4, and the SantaLucia (1998) salt correction [3][4].

**Intentionally simplified:**

- Default C_T = 250 nM differs from Biopython's 50 nM; **consequence:** NN Tm differs by
  the fixed concentration offset unless `primerConcentration` is set (passing 5e-8
  reproduces Biopython's 60.32).

**Not implemented:**

- Mg²⁺ / dNTP / Tris salt corrections and mismatch/dangling-end NN terms;
  **users should rely on:** dedicated tools (Biopython `Tm_NN`, MELTING 5) for those.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Default C_T = 250 nM vs Biopython 50 nM | Assumption | NN Tm offset under defaults | accepted | Explicit `primerConcentration` parameter; formula identical (ASM-02) |
| 2 | SEQ-TM-001 ↔ SEQ-THERMO-001 same methods | Deviation (duplicate Registry entry) | Avoid duplicate code/tests | accepted | Consolidated; see §7.3 and TestSpec §7 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty / length-1 | NN returns all-zero; Wallace/GC returns 0 | NN undefined for < 2 nt [3] |
| Lowercase input | Same as uppercase | Case-insensitive (upper-cased internally) |
| All-A·T duplex | Low (possibly negative) Tm | A·T pairs least stable [4] |

### 6.2 Limitations

Wallace/GC formulas are rules of thumb; the NN model assumes two-state melting, fixed
buffer, and monovalent (Na⁺) salt only. Mismatches, dangling ends, divalent ions, and
secondary structure are not modeled.

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
double wallace = SequenceStatistics.CalculateMeltingTemperature("ACGTTGCAATGCCGTA"); // 48.0 °C [5]
var nn = SequenceStatistics.CalculateThermodynamics(
    "CGTTCCAAAGATGTGGGCATGAGCTTAC", naConcentration: 0.05, primerConcentration: 5e-8);
// nn.MeltingTemperature == 60.3 (Biopython Tm_NN reference 60.32) [5]
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateThermodynamics_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceStatistics_CalculateThermodynamics_Tests.cs) — covers `INV-01`–`INV-05` (shared canonical fixture)
- Evidence: [SEQ-TM-001-Evidence.md](../../../docs/Evidence/SEQ-TM-001-Evidence.md)
- Related algorithms: [DNA_Thermodynamics.md](DNA_Thermodynamics.md) (SEQ-THERMO-001, full ΔH°/ΔS°/ΔG° derivation)

## 8. References

1. Thein, S. L., & Wallace, R. B. 1986. The use of synthetic oligonucleotides as specific hybridization probes in the diagnosis of genetic disorders. In *Human Genetic Diseases: A Practical Approach*, 33–50. (Wallace rule, as cited in Biopython `Tm_Wallace` docstring)
2. Marmur, J., & Doty, P. 1962. Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature. *J Mol Biol* 5:109–118. (GC formula, as cited in Biopython `Tm_GC` valueset 1)
3. Allawi, H. T., & SantaLucia, J. 1997. Thermodynamics and NMR of internal G·T mismatches in DNA. *Biochemistry* 36(34):10581–10594. https://doi.org/10.1021/bi962590c
4. SantaLucia, J. 1998. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS* 95(4):1460–1465. https://doi.org/10.1073/pnas.95.4.1460
5. Cock, P. J. A. et al. Biopython, `Bio.SeqUtils.MeltingTemp` (`Tm_Wallace`, `Tm_GC`, `Tm_NN`). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py and https://biopython.org/docs/1.76/api/Bio.SeqUtils.MeltingTemp.html (accessed 2026-06-14)

- Related algorithm: [Melting_Temperature.md](../MolTools/Melting_Temperature.md) (PRIMER-TM-001 — the primer/nearest-neighbor Tm used by `PrimerDesigner`, distinct from this sequence-statistics implementation).
