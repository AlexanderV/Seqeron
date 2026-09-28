# Melting Temperature Calculation

| Field | Value |
|-------|-------|
| Algorithm Group | Molecular Tools |
| Test Unit ID | PRIMER-TM-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Melting temperature (Tm) is the temperature at which 50% of the DNA duplex molecules are dissociated into single strands. In this repository, melting temperature estimation supports primer-oriented DNA calculations based on short-oligo and longer-oligo formulas, plus an optional sodium correction. The documented implementation is a closed-form estimator rather than a full nearest-neighbor thermodynamic model, so it is appropriate for fast screening and heuristic design workflows rather than detailed duplex thermodynamics.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

DNA duplex stability depends on hydrogen bonding, base stacking, ionic conditions, and sequence length. G-C pairs contribute more stability than A-T pairs because they form three hydrogen bonds rather than two, and cations stabilize the negatively charged phosphate backbone. These effects motivate the repository's use of base-composition formulas for quick Tm estimation. Source: Wikipedia (Nucleic acid thermodynamics), SantaLucia (1998).

### 2.2 Core Model

For short oligonucleotides, the documented model uses the Wallace rule:

$$
T_m = 2 \times (A + T) + 4 \times (G + C)
$$

where $(A + T)$ is the count of adenine and thymine bases and $(G + C)$ is the count of guanine and cytosine bases. Source: Thein & Wallace (1986), as cited in the original document.

For longer primers, the documented model uses the Marmur-Doty formula:

$$
T_m = 64.9 + \frac{41 \times (GC - 16.4)}{N}
$$

where $GC$ is the number of G and C bases and $N$ is the counted sequence length. Source: Marmur & Doty (1962).

Both basic formulas are OligoCalc's "Basic" Tm (Kibbe 2007, NAR 35:W43) and assume fixed conditions of
50 nM primer, 50 mM Na+ and pH 7.0 (the 64.9 constant already contains the salt term: 64.9 + 41(GC−16.4)/N =
64.9 + 0.41·%GC − 672.4/N).

The salt-adjusted variant is OligoCalc's "Salt Adjusted" Tm ([Na+] in mol/L):

$$
T_m^{salt} = 2(A+T) + 4(G+C) + 16.6\log_{10}\frac{[Na^+]}{0.050} \quad (N < 14)
$$

$$
T_m^{salt} = 100.5 + \frac{41(G+C)}{N} - \frac{820}{N} + 16.6\log_{10}[Na^+] \quad (N \ge 14)
$$

The 16.6·log10 slope is Schildkraut & Lifson (1965). An absolute 16.6·log10([Na+]) term is only valid for formulas
referenced to 1 M Na+ (e.g. 81.5 + 0.41·%GC − 600/N, Biopython `Tm_GC` valueset 7); adding it to the 50 mM-based
basic formulas double-counts salt (the pre-2026-09 implementation did this: 20-mer 50 %GC @ 50 mM → 30.2 °C).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The Wallace-rule estimate is linear in base counts with coefficients 2 for A/T and 4 for G/C | That is the documented closed-form rule |
| INV-02 | The Marmur-Doty estimate depends only on GC count and counted length | That is the documented formula |
| INV-03 | The salt-adjusted Tm rises by 16.6 °C per decade of [Na+]; for N<14 it equals the basic Tm at 50 mM | OligoCalc formula (relative correction to the 50 mM basis) |
| INV-04 | The implemented base estimator returns a value no lower than 0 for the longer-primer branch | `PrimerDesigner.CalculateMeltingTemperature` clamps the Marmur-Doty branch with `Math.Max(0, ...)` in source |

### 2.5 Comparison with Related Methods

| Aspect | This implementation | Nearest-neighbor thermodynamics |
|--------|---------------------|---------------------------------|
| Core inputs | Base counts, length, and optional sodium correction | Context-dependent dinucleotide stacking parameters |
| Computational model | Closed-form approximation | Detailed thermodynamic model |
| Intended use in current docs | Fast screening for primer-like oligos | Higher-fidelity duplex stability estimation |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `primer` | `string` | required | DNA sequence to score | Case-insensitive; only `A/C/G/T` contribute to counted length in `PrimerDesigner.CalculateMeltingTemperature(...)` |
| `naConcentration` | `double` | `50` | Sodium concentration in mM for the salt-adjusted estimate | Used only by `PrimerDesigner.CalculateMeltingTemperatureWithSalt(...)`; must be > 0 and finite (else `ArgumentOutOfRangeException`) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `tm` | `double` | Estimated melting temperature in degrees Celsius |

### 3.3 Preconditions and Validation

`PrimerDesigner.CalculateMeltingTemperature(...)` returns `0` for null or empty input and converts input to uppercase before counting bases. Only standard DNA bases `A/C/G/T` are counted; other characters are ignored when computing the valid length and GC count. `PrimerDesigner.CalculateMeltingTemperatureWithSalt(...)` also returns `0` for null or empty input (or no A/C/G/T), evaluates the OligoCalc salt-adjusted formula with [Na+] converted from mM to M, and rounds the result to one decimal place.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input sequence to uppercase.
2. Count A/T and G/C bases and derive the counted DNA length.
3. Return `0` if no counted DNA bases are present.
4. If the counted length is less than 14, apply the Wallace rule.
5. Otherwise apply the Marmur-Doty formula and clamp the result to `>= 0`.
6. For the salt-adjusted variant, evaluate the OligoCalc salt-adjusted formula for the same length branch and round to one decimal place.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The implementation centralizes the formula constants in `ThermoConstants`:

| Constant | Value | Description |
|----------|-------|-------------|
| `WallaceMaxLength` | 14 | Threshold between short- and longer-oligo formulas |
| `WallaceAtContribution` | 2 | Degrees Celsius per A/T base in the Wallace rule |
| `WallaceGcContribution` | 4 | Degrees Celsius per G/C base in the Wallace rule |
| `MarmurDotyBase` | 64.9 | Base temperature constant in the Marmur-Doty formula |
| `MarmurDotyGcCoefficient` | 41.0 | GC coefficient in the Marmur-Doty formula |
| `MarmurDotyGcOffset` | 16.4 | GC offset term in the Marmur-Doty formula |
| `SaltCoefficient` | 16.6 | Schildkraut–Lifson salt slope |
| `OligoCalcSaltAdjustedBase` | 100.5 | Constant of the OligoCalc salt-adjusted formula (N ≥ 14) |
| `OligoCalcSaltAdjustedLengthFactor` | 820 | Length term of the OligoCalc salt-adjusted formula |
| `WallaceReferenceNaMolar` | 0.050 | [Na+] (M) at which the Wallace rule is defined |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateMeltingTemperature` | `O(n)` | `O(1)` | Counts bases in a single pass over the input |
| `CalculateMeltingTemperatureWithSalt` | `O(n)` | `O(1)` | Reuses the base estimate and adds a closed-form correction |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs), [ThermoConstants.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/ThermoConstants.cs)

- `PrimerDesigner.CalculateMeltingTemperature(string)`: Estimates DNA primer Tm using the Wallace rule for counted lengths below 14 and the Marmur-Doty formula otherwise.
- `PrimerDesigner.CalculateMeltingTemperatureWithSalt(string, double)`: OligoCalc salt-adjusted Tm ([Na+] in mM), rounded to one decimal place.
- `ThermoConstants.CalculateWallaceTm(int, int)`: Applies the short-oligo closed form.
- `ThermoConstants.CalculateMarmurDotyTm(int, int)`: Applies the longer-primer closed form.
- `ThermoConstants.CalculateOligoCalcSaltAdjustedTm(int, int, double)`: OligoCalc salt-adjusted Tm from base counts and molar [Na+].
- `ThermoConstants.CalculateSaltCorrection(double)`: absolute 16.6·log10([Na+]) term (Biopython `salt_correction` method 1); not used by the primer Tm.

### 5.2 Current Behavior

The current implementation is DNA-oriented and case-insensitive. In `PrimerDesigner.CalculateMeltingTemperature(...)`, only `A/C/G/T` contribute to the counted length, so ambiguous or non-DNA characters are ignored rather than rejected. The short-sequence branch switches at fewer than 14 counted bases, and the longer-sequence branch clamps negative estimates to `0`. The salt-adjusted variant evaluates the OligoCalc salt-adjusted formula and rounds the final result to one decimal place.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Wallace-rule estimation with A/T and G/C contributions of 2 and 4, respectively.
- Marmur-Doty estimation using `64.9 + 41 * (GC - 16.4) / N`.
- OligoCalc salt-adjusted Tm (Kibbe 2007), cross-checked against the OligoCalc output 78 °C for a 39-mer at 50 mM Na+.

**Intentionally simplified:**

- The repository uses base-composition formulas instead of a full nearest-neighbor model; **consequence:** sequence-context effects from dinucleotide stacking are not reflected in the reported Tm.
- The branch point is fixed at fewer than 14 counted DNA bases; **consequence:** users may see different estimates from tools that switch formulas at a different threshold.
- Non-`ACGT` characters are ignored during counting; **consequence:** degenerate primers can yield estimates based only on the standard DNA subset.

**Not implemented:**

- Full nearest-neighbor thermodynamic melting-temperature estimation in this method; **users should rely on:** `PrimerDesigner.CalculateMeltingTemperatureNN` (PRIMER-NNTM-001).

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Fixed Wallace/Marmur-Doty switch at 14 counted bases | Assumption | Results may differ from tools that switch at another length threshold | accepted | The original document notes that some literature uses thresholds around 17-20 bp |
| 2 | Non-`ACGT` characters are excluded from counted length | Deviation | Degenerate or malformed symbols do not contribute to the estimate | accepted | Confirmed in `PrimerDesigner.CalculateMeltingTemperature(...)` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null or empty input | Returns `0` | Explicit guard in the implementation |
| No counted DNA bases after normalization | Returns `0` | `validLength == 0` short-circuits the calculation |
| Lowercase DNA input | Computes the same result as uppercase input | Input is uppercased before counting |
| Counted length below 14 | Uses the Wallace rule | `WallaceMaxLength` is 14 in source |
| Salt-corrected call | Returns a one-decimal-place value | `CalculateMeltingTemperatureWithSalt(...)` rounds the final result |

### 6.2 Limitations

The current implementation does not model nearest-neighbor stacking effects, mixed buffer chemistries, or ambiguity-code thermodynamics. The sodium-adjusted path accounts only for an additive sodium term, and the base estimator treats non-`ACGT` characters as non-contributing symbols rather than rejecting them.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical / biological walk-through:**

Wallace-rule example for `GCGCGCGC`:

- `A + T = 0`
- `G + C = 8`
- `T_m = 2 × 0 + 4 × 8 = 32°C`

Marmur-Doty example for `ACGTACGTACGTACGTACGT`:

- `GC = 10`
- `T_m = 64.9 + 41 × (10 - 16.4) / 20`
- `T_m = 64.9 - 13.12 = 51.78°C`

Salt-adjusted example for `ACGTACGTACGTACGTACGT` at `50 mM Na+` (N = 20 ≥ 14):

- `T_m = 100.5 + 41 × 10/20 − 820/20 + 16.6 × log10(0.050)`
- `T_m = 100.5 + 20.5 − 41 − 21.60 = 58.4°C`

OligoCalc cross-check: `GAGCAGGATCCCTATAGAGTGACAAAAGGATCTTGGTCC` (39 nt, 19 GC) at 50 mM — basic 67.63 °C, salt-adjusted
77.85 °C; OligoCalc reports 67.6 °C and 78 °C.

### 7.2 Applications and Use Cases (Optional)

Typical Tm ranges documented for PCR-related use cases:

| Application | Recommended Tm | Notes |
|-------------|----------------|-------|
| Standard PCR | 55-65°C | Optimal annealing |
| High-fidelity PCR | 60-72°C | Higher specificity |
| Colony PCR | 50-55°C | Lower stringency |
| Real-time PCR | 58-62°C | Narrow range preferred |

## 8. References

1. Marmur, J. & Doty, P. (1962). Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature. J Mol Biol 5:109-118.
2. SantaLucia, J. Jr. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proc Natl Acad Sci USA 95:1460-5.
3. Kibbe, W.A. (2007). OligoCalc: an online oligonucleotide properties calculator. Nucleic Acids Res 35:W43-W46.
4. Schildkraut, C. & Lifson, S. (1965). Dependence of the melting temperature of DNA on salt concentration. Biopolymers 3:195-208.
5. Biopython `Bio.SeqUtils.MeltingTemp` (1.88): `Tm_Wallace`, `Tm_GC`, `salt_correction`.
6. Wikipedia: Nucleic acid thermodynamics. https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
7. Wikipedia: DNA melting. https://en.wikipedia.org/wiki/DNA_melting

- Related algorithm: [Melting_Temperature.md](../Statistics/Melting_Temperature.md) (SEQ-TM-001 — the sequence-statistics Tm in `SequenceStatistics`, distinct from this simplified primer implementation).
