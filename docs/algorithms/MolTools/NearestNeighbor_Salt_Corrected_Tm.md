# Nearest-Neighbour Salt-Corrected Melting Temperature (Primer Tm, opt-in)

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PRIMER-TM-001 |
| Related Projects | Seqeron.Genomics.MolTools |
| Implementation Status | Production |
| Last Reviewed | 2026-06-24 |

## 1. Overview

Computes the melting temperature (Tm) of a PCR primer / oligonucleotide using the
SantaLucia (1998) **unified nearest-neighbour (NN)** thermodynamics and the bimolecular Tm
equation, with a published salt correction. This is an **opt-in design Tm**: the default
`PrimerDesigner.CalculateMeltingTemperature` (Wallace rule / Marmur-Doty GC formula) is
unchanged. The NN model is specification-driven and exact for its parameter table: it sums
experimentally measured ΔH°/ΔS° stacking terms over each adjacent base-pair step, adds
helix-initiation, terminal-A·T and symmetry terms, then converts to Tm at a chosen strand
and salt concentration [1][2]. Monovalent ([Na⁺]) correction follows Owczarzy et al. (2004)
[3]; divalent ([Mg²⁺]/dNTP) correction follows Owczarzy et al. (2008) [4].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Tm depends on sequence, base composition (G·C pairs are more stable than A·T), strand
concentration, and ionic conditions. Empirical rules (Wallace, Marmur-Doty) approximate Tm
from base counts; the NN model is the accurate, literature-standard approach for short
oligos and is what primer-design tools (Primer3, IDT OligoAnalyzer) use [1][2].

### 2.2 Core Model

**Duplex thermodynamics (SantaLucia & Hicks 2004, Eq. 1, Table 1)** [2]:

> ΔH°(total) = ΔH°(init) + Σ ΔH°(NN stack) + ΔH°(terminal A·T per A·T-closed end)
> ΔS°(total) = ΔS°(init) + Σ ΔS°(NN stack) + ΔS°(terminal A·T) + ΔS°(symmetry, self-comp only)

**Melting temperature (Eq. 3)** [1][2]:

> Tm = ΔH° × 1000 / (ΔS° + R · ln(C_T / x)) − 273.15

with R = 1.9872 cal/(K·mol), x = 4 for non-self-complementary and x = 1 for
self-complementary duplexes, C_T = total molar strand concentration.

**Monovalent salt correction.** Two published forms are offered:
- SantaLucia & Hicks (2004) Eq. 5 entropy form [2]: ΔS°[Na] = ΔS°[1 M] + 0.368 · (N/2) · ln[Na⁺],
  N = total phosphates in the duplex = 2·(length − 1).
- Owczarzy et al. (2004) quadratic 1/Tm form [3]:
  1/Tm[Na] = 1/Tm[1 M] + (4.29·f(GC) − 3.95)·1e-5·ln[Na⁺] + 9.40e-6·(ln[Na⁺])².

**Divalent (Mg²⁺) correction (Owczarzy 2008)** [4]: a 1/Tm correction in ln[Mg²⁺] with a
GC-fraction and length term; the regime is selected by R = √[Mg²⁺]/[Mon] (R < 0.22 → monovalent
dominates; 0.22 ≤ R < 6 → mixed, with reparameterised coefficients; R ≥ 6 → divalent). Free
Mg²⁺ is reduced by dNTP chelation via Ka = 3×10⁴.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Two-state (all-or-none) duplex melting at fixed buffer | Tm mis-estimates partially melted/structured oligos |
| ASM-02 | Self-complementarity correctly detected (x, symmetry) | Wrong x or missing symmetry term shifts Tm |
| ASM-03 | Salt corrections valid in their published range ([Na⁺] 0.05–1.1 M, ≤16 bp for Eq. 5) | Extrapolation degrades accuracy [2][3] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | ΔH°/ΔS° = init + Σ NN + terminal-AT(per A·T end) + symmetry(self-comp) | Definition [2] Eq. 1, Table 1 |
| INV-02 | Self-complementary ⇒ x = 1, else x = 4 | Definition [2] Eq. 3 |
| INV-03 | Lower [Na⁺] ⇒ lower Owczarzy-2004 Tm | Quadratic correction sign [3] |
| INV-04 | Adding Mg²⁺ ⇒ higher Tm vs Mg²⁺-free buffer | Divalent stabilisation [4] |
| INV-05 | Divalent mode with [Mg²⁺]=0 ≡ monovalent 2004 mode | Method-7 fallback [4][5] |
| INV-06 | Empty / length-1 / non-ACGT input → not computable (null / NaN) | NN sums over ACGT dinucleotides [2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| primer | string | required | DNA primer 5'→3' | case-insensitive; ≥2 ACGT bases |
| strandConcentrationMolar | double | 0.5e-6 | Total strand conc. C_T (mol/L) | > 0 |
| sodiumMolar | double | 0.05 | Monovalent [Na⁺]+[K⁺]+[Tris]/2 (mol/L) | > 0 |
| magnesiumMolar | double | 0.0 | [Mg²⁺] (mol/L); divalent mode only | ≥ 0 |
| dntpMolar | double | 0.0 | Total dNTP (mol/L); chelates Mg²⁺ | ≥ 0 |
| saltMode | SaltCorrectionMode | Owczarzy2004Monovalent | None / SantaLuciaEntropy / Owczarzy2004Monovalent / Owczarzy2008Divalent | — |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `CalculateMeltingTemperatureNN` | double | Tm in °C, or NaN for invalid input |
| `CalculateNearestNeighborThermodynamics` | (double DeltaH, double DeltaS, bool IsSelfComplementary)? | ΔH° (kcal/mol), ΔS° (cal/(K·mol)), self-comp flag, or null |

### 3.3 Preconditions and Validation

Input is upper-cased (case-insensitive). Empty, null, length < 2, or any non-ACGT character
makes the NN lookup fail: `CalculateNearestNeighborThermodynamics` returns `null` and
`CalculateMeltingTemperatureNN` returns `double.NaN`. No exceptions are thrown for these
guarded inputs.

## 4. Algorithm

### 4.1 High-Level Steps

1. Sum NN ΔH°/ΔS° over each dinucleotide; add the duplex-initiation term.
2. Add the terminal-A·T penalty once per end closed by an A or T.
3. Detect self-complementarity (sequence equals its reverse complement); if so, add the
   symmetry ΔS° term and set x = 1, else x = 4.
4. (SantaLuciaEntropy mode) salt-correct ΔS° via Eq. 5 before the Tm equation.
5. Evaluate Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) − 273.15 (Kelvin internally).
6. Apply the Owczarzy 2004 (monovalent) or Owczarzy 2008 (divalent) 1/Tm correction, then
   convert Kelvin → °C.

### 4.2 Decision Rules, Scoring, Reference Tables

The unified NN ΔH°/ΔS° table (1 M NaCl), initiation (+0.2/−5.7), terminal-A·T (+2.2/+6.9)
and symmetry (0.0/−1.4) terms are SantaLucia & Hicks (2004) Table 1 [2], cross-checked
against Biopython `DNA_NN4` [5]. R = 1.9872, x ∈ {1,4}. Owczarzy 2004 coefficients
(4.29e-5, 3.95e-5, 9.40e-6) [3][5]; Owczarzy 2008 coefficients (a..g ×1e-5, Ka=3×10⁴, R
regime thresholds 0.22 / 6.0) [4][5]. Every constant is named and source-cited in code.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateNearestNeighborThermodynamics | O(n) | O(1) | single pass over dinucleotides |
| CalculateMeltingTemperatureNN | O(n) | O(1) | thermo pass + O(1) salt correction |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs)

- `PrimerDesigner.CalculateMeltingTemperatureNN(string, double, double, double, double, SaltCorrectionMode)`: opt-in NN salt-corrected Tm (perfectly complementary duplex).
- `PrimerDesigner.CalculateNearestNeighborThermodynamics(string)`: ΔH°/ΔS° + self-comp flag.
- `PrimerDesigner.CalculateMeltingTemperatureNNMismatch(string top, string bottom3to5, …, SaltCorrectionMode)`: **opt-in extension** — NN Tm for a probe–target duplex with a single internal mismatch and/or dangling end (mirrors Biopython `Tm_NN(imm_table=DNA_IMM, de_table=DNA_DE)`).
- `PrimerDesigner.CalculateNearestNeighborThermodynamicsMismatch(string top, string bottom3to5)`: ΔH°/ΔS° + self-comp flag for a mismatch/dangling-end duplex.
- `PrimerDesigner.SaltCorrectionMode`: None / SantaLuciaEntropy / Owczarzy2004Monovalent / Owczarzy2008Divalent.

### 5.2 Current Behavior

The default `CalculateMeltingTemperature` (Wallace/Marmur-Doty) and `Calculate3PrimeStability`
(NN ΔG°37 for 3'-end stability) are unchanged. The new NN Tm uses the **1998 unified**
parameters (DNA_NN4); note `SequenceStatistics.CalculateThermodynamics` (SEQ-THERMO-001) uses
the older **1997 (Allawi)** parameters and a different default concentration — a distinct unit,
not modified here. No substring search / matching is involved, so the repository suffix tree is
**not applicable** to this unit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Unified NN ΔH°/ΔS° table, initiation, terminal-A·T, symmetry terms [2].
- Tm = ΔH°·1000/(ΔS° + R·ln(C_T/x)) − 273.15 with R = 1.9872, x ∈ {1,4} [1][2].
- SantaLucia Eq. 5 entropy salt correction (N = 2·(L−1)) [2].
- Owczarzy 2004 monovalent quadratic 1/Tm correction [3][5].
- Owczarzy 2008 divalent Mg²⁺/dNTP correction with R-regime selection [4][5].
- Internal single-mismatch NN ΔH°/ΔS° (Allawi/SantaLucia 1997/1998; Peyret 1999) [6][7][9] and single
  dangling-end NN ΔH°/ΔS° (Bommarito 2000) [8], via `CalculateMeltingTemperatureNNMismatch` (opt-in).
  Convention mirrors Biopython `Tm_NN` (bottom strand 3'→5'; `top2/bottom2` keys tried forward then
  character-reversed; `.` marks the dangling base; terminal-AT from the un-dotted top termini). A
  perfectly paired duplex through this path equals the perfect-match `CalculateMeltingTemperatureNN`.

**Intentionally simplified:**

- Owczarzy 2004/2008 coefficients are taken from the Biopython reference implementation
  (the Biochemistry 43:3537 full text is paywalled); **consequence:** none — values are
  cross-corroborated and the published 35.8 °C worked example reproduces exactly.
- The internal-mismatch table covers a **single** internal mismatch (one mismatched column per NN step);
  two adjacent mismatches (a tandem mismatch) or a non-ACGT character yield no NN parameter → not
  computable (null/NaN). Terminal mismatches and coaxial stacking are out of scope.

**Not implemented:**

- Hairpin / secondary-structure Tm (folding-based melting); **users should rely on:** dedicated
  folding tools (UNAFold, ViennaRNA, MELTING 5) for those.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Default C_T = 0.5 µM | Assumption | Default operating point | accepted | Explicit parameter; formula identical (ASM-01) |
| 2 | Owczarzy coefficients via Biopython | Assumption | Source provenance | accepted | Cross-checked vs the 35.8 °C worked example and web search |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty / length-1 / non-ACGT | null (thermo) / NaN (Tm) | NN undefined without an ACGT dinucleotide step [2] |
| Self-complementary | x = 1 + symmetry ΔS° | Eq. 3 / Table 1 [2] |
| Both ends A·T | terminal-A·T penalty applied twice | Per-end penalty [2] |
| Lowercase input | same as uppercase | case-insensitive |
| Single internal mismatch | mismatch NN term applied | Allawi/SantaLucia/Peyret [6][7][9] |
| Single dangling end (`.` marker) | dangling-end NN term applied | Bommarito 2000 [8] |
| Perfect duplex via `*Mismatch` path | equals perfect-match path | strict superset |
| Tandem mismatch / unequal length / null | null (thermo) / NaN (Tm) | no NN parameter |

### 6.2 Limitations

The NN model assumes two-state melting and a fixed buffer. A **single** internal mismatch and a
**single** dangling end are now modelled (opt-in `*Mismatch` path); tandem/adjacent mismatches,
terminal mismatches, coaxial stacking, and hairpin/secondary-structure Tm are not — use a folding
tool (UNAFold, ViennaRNA, MELTING 5). The salt corrections are valid within their published ranges
([Na⁺] 0.05–1.1 M; Eq. 5 for ≤16 bp). Outside these ranges Tm is an extrapolation.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Opt-in NN Tm with the default Owczarzy-2004 monovalent correction at 50 mM Na+.
double tm = PrimerDesigner.CalculateMeltingTemperatureNN("ATGCATGC", sodiumMolar: 0.05);
// 18.19 °C  (x=4, fGC=0.5)

// Raw duplex thermodynamics.
var t = PrimerDesigner.CalculateNearestNeighborThermodynamics("GCGCGC");
// t.DeltaH == -50.4, t.DeltaS == -134.7, t.IsSelfComplementary == true

// Opt-in extension: a probe–target duplex with one internal T·G mismatch
// (bottom strand written 3'→5', aligned base-for-base under the top).
double tmMm = PrimerDesigner.CalculateMeltingTemperatureNNMismatch(
    "CGTGAC", "GCGCTG", saltMode: PrimerDesigner.SaltCorrectionMode.None);
// ΔH° = -35.5, ΔS° = -101.5  →  Tm = -6.41 °C (x=4)

// A 5'-dangling A on a GCGCGC core ('.' marks the unpaired base).
double tmDe = PrimerDesigner.CalculateMeltingTemperatureNNMismatch(
    "AGCGCGC", ".CGCGCG", saltMode: PrimerDesigner.SaltCorrectionMode.None);
// ΔH° = -51.9, ΔS° = -136.4  →  Tm = 35.80 °C
```

**Numerical walk-through:** ATGCATGC (non-self-comp, x=4): stacks AT+TG+GC+CA+AT+TG+GC give
ΔH° = 0.2 + (−7.2−8.5−9.8−8.5−7.2−8.5−9.8) + 2.2 + 2.2 = −57.1; ΔS° = −5.7 + (−20.4−22.7−24.4
−22.7−20.4−22.7−24.4) + 6.9 + 6.9 = −156.5. At C_T = 0.5 µM, 1 M NaCl: Tm = −57.1·1000/(−156.5
+ 1.9872·ln(0.5e-6/4)) − 273.15 = 30.43 °C.

### 7.2 Applications and Use Cases

- **PCR primer design:** accurate annealing-temperature estimation for primers and probes
  under realistic salt/Mg²⁺ buffers (the model Primer3/IDT use) [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PrimerDesigner_NearestNeighborTm_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_NearestNeighborTm_Tests.cs) — covers `INV-1`–`INV-9` (incl. mismatch/dangling-end MM1/DE1, perfect-match equivalence EQ1/EQ2)
- Evidence: [PRIMER-TM-001-NN-Evidence.md](../../../docs/Evidence/PRIMER-TM-001-NN-Evidence.md)
- Related algorithms: [Melting_Temperature.md](../Statistics/Melting_Temperature.md) (SEQ-TM-001), [DNA_Thermodynamics.md](../Statistics/DNA_Thermodynamics.md) (SEQ-THERMO-001)

### 7.4 Change History

| Date | Version | Changes |
|------|---------|---------|
| 2026-06-24 | 1.0 | Initial NN salt-corrected Tm (opt-in) under PRIMER-TM-001 |
| 2026-06-24 | 1.1 | Added internal single-mismatch (Allawi/SantaLucia/Peyret) + single dangling-end (Bommarito 2000) NN terms via `*Mismatch` path (opt-in extension) |

## 8. References

1. SantaLucia J. 1998. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS* 95(4):1460–1465. https://doi.org/10.1073/pnas.95.4.1460
2. SantaLucia J, Hicks D. 2004. The thermodynamics of DNA structural motifs. *Annu Rev Biophys Biomol Struct* 33:415–440. https://doi.org/10.1146/annurev.biophys.32.110601.141800
3. Owczarzy R, You Y, Moreira BG, et al. 2004. Effects of sodium ions on DNA duplex oligomers: improved predictions of melting temperatures. *Biochemistry* 43(12):3537–3554. https://doi.org/10.1021/bi034621r
4. Owczarzy R, Moreira BG, You Y, et al. 2008. Predicting stability of DNA duplexes in solutions containing magnesium and monovalent cations. *Biochemistry* 47(19):5336–5353. https://doi.org/10.1021/bi702363u
5. Cock PJA et al. Biopython, `Bio.SeqUtils.MeltingTemp` (`DNA_NN4`, `DNA_IMM`, `DNA_DE`, `Tm_NN`, `salt_correction`). https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/MeltingTemp.py (accessed 2026-06-24)
6. Allawi HT, SantaLucia J. 1997. Thermodynamics and NMR of internal G·T mismatches in DNA. *Biochemistry* 36(34):10581–10594. https://doi.org/10.1021/bi962590c
7. Allawi HT, SantaLucia J. 1998. Internal G·A (Biochemistry 37:2170), C·T (Nucleic Acids Res 26:2694), and A·C (Biochemistry 37:9435) mismatch parameters. https://doi.org/10.1021/bi9724873
8. Bommarito S, Peyret N, SantaLucia J. 2000. Thermodynamic parameters for DNA sequences with dangling ends. *Nucleic Acids Res* 28(9):1929–1934. https://doi.org/10.1093/nar/28.9.1929
9. Peyret N, Seneviratne PA, Allawi HT, SantaLucia J. 1999. Nearest-neighbor thermodynamics and NMR of DNA sequences with internal A·A, C·C, G·G, and T·T mismatches. *Biochemistry* 38(12):3468–3477. https://doi.org/10.1021/bi9825091
