# Melting Temperature Calculation

| Field | Value |
|-------|-------|
| Algorithm Group | Molecular Tools |
| Test Unit ID | PRIMER-TM-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

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

The salt-corrected variant adds a sodium correction:

$$
T_m^{corrected} = T_m^{base} + 16.6 \times \log_{10}\left(\frac{[Na^+]}{1000}\right)
$$

where $[Na^+]$ is sodium concentration in mM. Source: Owczarzy et al. (2004), as cited in the original document.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The Wallace-rule estimate is linear in base counts with coefficients 2 for A/T and 4 for G/C | That is the documented closed-form rule |
| INV-02 | The Marmur-Doty estimate depends only on GC count and counted length | That is the documented formula |
| INV-03 | The sodium correction is additive in $\log_{10}([Na^+]/1000)$ | That is the documented salt-correction form |
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
| `naConcentration` | `double` | `50` | Sodium concentration in mM for the corrected estimate | Used only by `PrimerDesigner.CalculateMeltingTemperatureWithSalt(...)` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `tm` | `double` | Estimated melting temperature in degrees Celsius |

### 3.3 Preconditions and Validation

`PrimerDesigner.CalculateMeltingTemperature(...)` returns `0` for null or empty input and converts input to uppercase before counting bases. Only standard DNA bases `A/C/G/T` are counted; other characters are ignored when computing the valid length and GC count. `PrimerDesigner.CalculateMeltingTemperatureWithSalt(...)` also returns `0` for null or empty input, adds the sodium correction in mM units, and rounds the result to one decimal place.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input sequence to uppercase.
2. Count A/T and G/C bases and derive the counted DNA length.
3. Return `0` if no counted DNA bases are present.
4. If the counted length is less than 14, apply the Wallace rule.
5. Otherwise apply the Marmur-Doty formula and clamp the result to `>= 0`.
6. For the salt-corrected variant, add the sodium correction and round to one decimal place.

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
| `SaltCoefficient` | 16.6 | Sodium correction coefficient |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateMeltingTemperature` | `O(n)` | `O(1)` | Counts bases in a single pass over the input |
| `CalculateMeltingTemperatureWithSalt` | `O(n)` | `O(1)` | Reuses the base estimate and adds a closed-form correction |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PrimerDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs), [ThermoConstants.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/ThermoConstants.cs)

- `PrimerDesigner.CalculateMeltingTemperature(string)`: Estimates DNA primer Tm using the Wallace rule for counted lengths below 14 and the Marmur-Doty formula otherwise.
- `PrimerDesigner.CalculateMeltingTemperatureWithSalt(string, double)`: Adds a sodium correction in mM units and rounds the corrected value to one decimal place.
- `ThermoConstants.CalculateWallaceTm(int, int)`: Applies the short-oligo closed form.
- `ThermoConstants.CalculateMarmurDotyTm(int, int)`: Applies the longer-primer closed form.
- `ThermoConstants.CalculateSaltCorrection(double)`: Computes the additive sodium correction from mM concentration.

### 5.2 Current Behavior

The current implementation is DNA-oriented and case-insensitive. In `PrimerDesigner.CalculateMeltingTemperature(...)`, only `A/C/G/T` contribute to the counted length, so ambiguous or non-DNA characters are ignored rather than rejected. The short-sequence branch switches at fewer than 14 counted bases, and the longer-sequence branch clamps negative estimates to `0`. The corrected variant adds `16.6 * log10([Na+]/1000)` using the provided sodium concentration in mM and rounds the final result to one decimal place.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Wallace-rule estimation with A/T and G/C contributions of 2 and 4, respectively.
- Marmur-Doty estimation using `64.9 + 41 * (GC - 16.4) / N`.
- Additive sodium correction with coefficient `16.6` in log-space.

**Intentionally simplified:**

- The repository uses base-composition formulas instead of a full nearest-neighbor model; **consequence:** sequence-context effects from dinucleotide stacking are not reflected in the reported Tm.
- The branch point is fixed at fewer than 14 counted DNA bases; **consequence:** users may see different estimates from tools that switch formulas at a different threshold.
- Non-`ACGT` characters are ignored during counting; **consequence:** degenerate primers can yield estimates based only on the standard DNA subset.

**Not implemented:**

- Full nearest-neighbor thermodynamic melting-temperature estimation; **users should rely on:** no current alternative documented in this test unit.

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

Salt-correction example using `51.78°C` at `50 mM Na+`:

- `correction = 16.6 × log10(50 / 1000)`
- `correction = 16.6 × log10(0.05)`
- `correction = 16.6 × (-1.301) = -21.6°C`
- `T_m^{corrected} = 51.78 - 21.6 = 30.2°C`

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
3. Owczarzy, R. et al. (2004). Effects of sodium ions on DNA duplex oligomers: improved predictions of melting temperatures. Biochemistry 43:3537-3554.
4. Wikipedia: Nucleic acid thermodynamics. https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
5. Wikipedia: DNA melting. https://en.wikipedia.org/wiki/DNA_melting

- Related algorithm: [Melting_Temperature.md](../Statistics/Melting_Temperature.md) (SEQ-TM-001 — the sequence-statistics Tm in `SequenceStatistics`, distinct from this simplified primer implementation).
