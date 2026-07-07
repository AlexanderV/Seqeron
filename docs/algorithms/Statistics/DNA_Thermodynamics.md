# DNA Duplex Thermodynamics (Nearest-Neighbor)

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics |
| Test Unit ID | SEQ-THERMO-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Infrastructure |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Predicts the thermodynamic stability of a short DNA duplex from its sequence using the
unified nearest-neighbor (NN) model of Allawi & SantaLucia (1997) / SantaLucia (1998).
It returns the hybridization enthalpy ΔH° (kcal/mol), entropy ΔS° (cal/(mol·K)),
Gibbs free energy ΔG°₃₇ (kcal/mol) and melting temperature Tm (°C). The method is
specification-driven: outputs are determined exactly by the published NN parameter table,
the Na⁺ salt correction, and the Tm equation. A simpler `CalculateMeltingTemperature`
delegate provides Wallace-rule / Marmur-Doty estimates for quick screening.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The helix-coil transition of a nucleic acid duplex can be modeled as the sum of stacking
interactions between adjacent base pairs ("nearest neighbors"), plus end (initiation)
terms. Two duplexes with the same set of nearest-neighbor pairs have the same predicted
thermodynamics regardless of overall sequence [4].

### 2.2 Core Model

Total enthalpy and entropy are sums over the 10 unique Watson-Crick NN steps plus an
initiation contribution applied at **each** duplex terminus [1][2][3]:

- ΔH°(total) = Σ ΔH°(NN step) + ΔH°(init, 5′ end) + ΔH°(init, 3′ end)
- ΔS°(total) = Σ ΔS°(NN step) + ΔS°(init, 5′ end) + ΔS°(init, 3′ end)

with a per-terminus initiation chosen by whether that terminal base pair is G·C or A·T [1][3].
The Na⁺ salt entropy correction (SantaLucia 1998 "method 5") is [2][3]:

- ΔS° ← ΔS° + 0.368 · (N − 1) · ln[Na⁺],  [Na⁺] in mol/L.

Gibbs free energy at 37 °C (310.15 K) [2]:

- ΔG°₃₇ = ΔH° − 310.15 · ΔS° / 1000  (ΔS° converted cal→kcal).

Melting temperature [3][4]:

- Tm = (1000 · ΔH°) / (ΔS° + R · ln(C_T / F)) − 273.15,  R = 1.987 cal/(mol·K).

F = 4 for two non-self-complementary strands in equimolar amount (the default);
F = 1 for self-complementary duplexes [4].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Two-state (all-or-none) helix-coil transition. | Tm prediction degrades for sequences forming partial/hairpin structures. |
| ASM-02 | Duplex is non-self-complementary, equimolar strands (F = 4). | Self-complementary duplexes need F = 1; Tm would be off by R·ln(4) in the denominator. |
| ASM-03 | Fixed monovalent-cation buffer; only [Na⁺] entropy correction applied. | Mg²⁺ and other ions are not modeled; Tm under PCR-like Mg²⁺ conditions is approximate. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Initiation contributes once at each terminus (two init terms). | Model adds init for `seq[0]` and `seq[-1]` [1][3]. |
| INV-02 | ΔG°₃₇ = ΔH° − 310.15·ΔS°/1000. | Gibbs relation at 310.15 K [2]. |
| INV-03 | Tm = (1000·ΔH°)/(ΔS° + R·ln(C_T/4)) − 273.15, R = 1.987. | NN Tm equation [3][4]. |
| INV-04 | NN table is Watson-Crick symmetric (AA=TT, CA=TG, GT=AC, CT=AG, GA=TC, GG=CC). | DNA_NN3 parameter symmetry [1][3]. |
| INV-05 | Deterministic and case-insensitive. | Input upper-cased; no randomness [3]. |
| INV-06 | Empty or length-1 input returns (0,0,0,0). | NN model undefined for length < 2 (no dinucleotide). |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `dnaSequence` | `string` | required | DNA sequence 5′→3′. | length ≥ 2 for a non-zero result; A/C/G/T; case-insensitive; non-ACGT dinucleotides contribute 0. |
| `naConcentration` | `double` | 0.05 | Na⁺ concentration. | mol/L (0.05 = 50 mM); > 0. |
| `primerConcentration` | `double` | 2.5e-7 | Total strand concentration C_T. | mol/L (2.5e-7 = 250 nM); divided by F = 4. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `DeltaH` | `double` | ΔH° in kcal/mol, rounded to 2 decimals. |
| `DeltaS` | `double` | ΔS° in cal/(mol·K) including salt correction, rounded to 2 decimals. |
| `DeltaG` | `double` | ΔG°₃₇ in kcal/mol, rounded to 2 decimals. |
| `MeltingTemperature` | `double` | Tm in °C, rounded to 1 decimal. |

### 3.3 Preconditions and Validation

Input is upper-cased via `ToUpperInvariant` (case-insensitive). Empty or length-1 input
returns `(0,0,0,0)`. Indexing is 0-based; dinucleotides are overlapping windows of length 2.
Only A/C/G/T are recognized; an unrecognized dinucleotide adds 0 (no exception). No exceptions
are thrown for valid-typed input.

## 4. Algorithm

### 4.1 High-Level Steps

1. Return `(0,0,0,0)` if input is null/empty or length < 2.
2. Upper-case the sequence.
3. Add initiation (ΔH, ΔS) for the first base and for the last base by G·C vs A·T.
4. For each overlapping dinucleotide, add its NN (ΔH, ΔS) from the parameter table.
5. Apply the Na⁺ salt entropy correction: ΔS += 0.368·(N−1)·ln[Na⁺].
6. Compute ΔG°₃₇ = ΔH° − 310.15·ΔS°/1000.
7. Compute Tm = (1000·ΔH°)/(ΔS° + R·ln(C_T/4)) − 273.15.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

NN parameter table (1 M NaCl; ΔH kcal/mol, ΔS cal/(mol·K)) [1][3]:

| Step (and WC complement) | ΔH° | ΔS° |
|--------------------------|-----|-----|
| AA / TT | −7.9 | −22.2 |
| AT | −7.2 | −20.4 |
| TA | −7.2 | −21.3 |
| CA / TG | −8.5 | −22.7 |
| GT / AC | −8.4 | −22.4 |
| CT / AG | −7.8 | −21.0 |
| GA / TC | −8.2 | −22.2 |
| CG | −10.6 | −27.2 |
| GC | −9.8 | −24.4 |
| GG / CC | −8.0 | −19.9 |
| init. w/ term. G·C | +0.1 | −2.8 |
| init. w/ term. A·T | +2.3 | +4.1 |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateThermodynamics` | O(n) | O(1) | one pass over n−1 dinucleotides; fixed 16-entry table. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.CalculateThermodynamics(string, double, double)`: NN ΔH°/ΔS°/ΔG°/Tm.
- `SequenceStatistics.CalculateMeltingTemperature(string, bool)`: simple Wallace / Marmur-Doty Tm (delegates to [ThermoConstants](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/ThermoConstants.cs)).

### 5.2 Current Behavior

NN parameters live in a static dictionary; both Watson-Crick complements are listed so any
input dinucleotide resolves directly without computing the complement. Initiation is applied
to both termini via a private `AddTerminalInitiation` helper. Salt correction, R, the reference
temperature, the Kelvin offset, and F are named constants. No string search/matching is
performed (single linear pass over dinucleotides), so the repository suffix tree is **not
applicable** to this unit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Full 10-step unified NN ΔH°/ΔS° table (Allawi & SantaLucia 1997 / DNA_NN3) [1][3].
- Two-terminus initiation by G·C vs A·T [1][3].
- Na⁺ "method 5" salt entropy correction 0.368·(N−1)·ln[Na⁺] [2][3].
- Tm = (1000·ΔH°)/(ΔS° + R·ln(C_T/4)) − 273.15, R = 1.987, F = 4 [3][4].
- ΔG°₃₇ = ΔH° − 310.15·ΔS°/1000 [2].

**Intentionally simplified:**

- Self-complementarity: fixed F = 4 (non-self-complementary equimolar); **consequence:** Tm for a self-complementary duplex (true F = 1) is under-estimated by the R·ln(4) denominator difference.

**Not implemented:**

- Mismatch, dangling-end, and Mg²⁺/mixed-cation corrections; **users should rely on:** dedicated tools (e.g. MELTING 5, IDT OligoAnalyzer) for those conditions.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Empty/length-1 returns (0,0,0,0) | Assumption | API edge convention; no thermodynamic value affected | accepted | ASM/INV-06 |
| 2 | Single-terminus initiation in prior code | Deviation | Under-counted one init term (ΔH/ΔS/Tm wrong) | fixed | Corrected to two-end init in SEQ-THERMO-001 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty / length-1 | (0,0,0,0) | No dinucleotide; NN undefined [3]. |
| Lowercase / mixed case | Same as upper-case | Input upper-cased [3]. |
| Non-ACGT dinucleotide | Contributes 0 to ΔH/ΔS | Not in NN table; no exception. |

### 6.2 Limitations

Two-state model only; not valid for sequences dominated by hairpins or partial duplexes.
No Mg²⁺/mixed-cation correction. Self-complementary duplexes use the non-self-complementary
factor. Predictions are most accurate for short oligonucleotides under standard buffer.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var t = SequenceStatistics.CalculateThermodynamics("GCGC");
// t.DeltaH = -30.0, t.DeltaS = -84.91, t.DeltaG = -3.67, t.MeltingTemperature = -18.6
```

**Numerical walk-through (GCGC, Na = 0.05 M, C_T = 250 nM):**

- Init (G end + C end): ΔH = 2×0.1 = 0.2; ΔS = 2×(−2.8) = −5.6.
- NN steps GC, CG, GC: ΔH = −9.8−10.6−9.8 = −30.2; ΔS = −24.4−27.2−24.4 = −76.0.
- ΔH° = 0.2 − 30.2 = −30.0.
- Salt ΔS = 0.368×3×ln(0.05) = −3.307; ΔS° = −5.6 − 76.0 − 3.307 = −84.91.
- ΔG°₃₇ = −30.0 − 310.15×(−84.91)/1000 = −3.67.
- Tm = (−30000)/(−84.91 + 1.987·ln(2.5e-7/4)) − 273.15 = −18.6 °C.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateThermodynamics_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceStatistics_CalculateThermodynamics_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [SEQ-THERMO-001-Evidence.md](../../../docs/Evidence/SEQ-THERMO-001-Evidence.md)
- Related algorithms: [Molecular_Weight_Calculation](./Molecular_Weight_Calculation.md)

## 8. References

1. Allawi HT, SantaLucia J Jr. 1997. Thermodynamics and NMR of internal G·T mismatches in DNA. Biochemistry 36(34):10581–10594. https://doi.org/10.1021/bi962590c
2. SantaLucia J Jr. 1998. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. PNAS 95(4):1460–1465. https://doi.org/10.1073/pnas.95.4.1460
3. Biopython. Bio.SeqUtils.MeltingTemp (DNA_NN3 table, Tm_NN). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py
4. Dumousseau M, Rodriguez N, Juty N, Le Novère N. 2012. MELTING, a flexible platform to predict the melting temperatures of nucleic acids. BMC Bioinformatics 13:101. User guide: https://www.ebi.ac.uk/biomodels/tools/melting/melting5-UserGuide.pdf
5. Wikipedia. Nucleic acid thermodynamics. https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
