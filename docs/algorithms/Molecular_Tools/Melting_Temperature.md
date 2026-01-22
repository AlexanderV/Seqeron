# Melting Temperature Calculation

## Overview

Melting temperature (Tm) is the temperature at which 50% of the DNA duplex molecules are dissociated into single strands. This is a fundamental parameter in molecular biology for designing PCR primers, hybridization probes, and other oligonucleotide-based applications.

## Theory

### Thermodynamic Basis

DNA double helix stability depends on:
1. **Hydrogen bonding**: G-C pairs form 3 hydrogen bonds (stronger), A-T pairs form 2 (weaker)
2. **Base stacking**: π-π interactions between adjacent bases
3. **Salt concentration**: Cations stabilize the negatively charged phosphate backbone
4. **Sequence length**: Longer duplexes have higher Tm

**Source:** Wikipedia (Nucleic acid thermodynamics), SantaLucia (1998) PNAS 95:1460-1465

### Formulas Implemented

#### 1. Wallace Rule (Short Oligonucleotides)

For primers shorter than 14 bp:

$$T_m = 2 \times (A + T) + 4 \times (G + C)$$

Where:
- (A + T) = count of adenine and thymine bases
- (G + C) = count of guanine and cytosine bases

**Rationale:** Simple approximation based on hydrogen bond count. Each A-T pair contributes ~2°C, each G-C pair contributes ~4°C.

**Source:** Thein & Wallace (1986), Oligonucleotide Synthesis: A Practical Approach

**Limitations:**
- Only accurate for very short oligonucleotides (< 14 bp)
- Does not account for base stacking or salt concentration
- Assumes ideal conditions

#### 2. Marmur-Doty Formula (Longer Primers)

For primers 14 bp or longer:

$$T_m = 64.9 + \frac{41 \times (GC - 16.4)}{N}$$

Where:
- GC = number of G and C bases
- N = total sequence length
- 64.9 = base temperature constant
- 41.0 = GC coefficient
- 16.4 = GC offset correction

**Source:** Marmur & Doty (1962) J Mol Biol 5:109-118, "Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature"

**Limitations:**
- Assumes standard salt conditions
- Does not account for nearest-neighbor effects
- Less accurate than nearest-neighbor thermodynamic methods

#### 3. Salt Correction

Adjusts Tm based on sodium ion concentration:

$$T_m^{corrected} = T_m^{base} + 16.6 \times \log_{10}\left(\frac{[Na^+]}{1000}\right)$$

Where:
- [Na+] = sodium concentration in mM
- Default: 50 mM (typical PCR conditions)

**Source:** Owczarzy et al. (2004) Biochemistry 43:3537-3554

## Implementation

### Constants (ThermoConstants.cs)

| Constant | Value | Description |
|----------|-------|-------------|
| `WallaceMaxLength` | 14 | Threshold for Wallace vs Marmur-Doty |
| `WallaceAtContribution` | 2 | °C per A/T base pair |
| `WallaceGcContribution` | 4 | °C per G/C base pair |
| `MarmurDotyBase` | 64.9 | Base temperature (°C) |
| `MarmurDotyGcCoefficient` | 41.0 | GC coefficient |
| `MarmurDotyGcOffset` | 16.4 | GC offset |
| `SaltCoefficient` | 16.6 | Salt adjustment coefficient |

### Methods (PrimerDesigner.cs)

#### CalculateMeltingTemperature(string primer)

```csharp
public static double CalculateMeltingTemperature(string primer)
```

**Behavior:**
1. Returns 0 for null or empty input
2. Converts input to uppercase for case-insensitive matching
3. If length < 14: uses Wallace rule
4. If length ≥ 14: uses Marmur-Doty formula
5. Result is clamped to ≥ 0

**Complexity:** O(n) where n = primer length

#### CalculateMeltingTemperatureWithSalt(string primer, double naConcentration)

```csharp
public static double CalculateMeltingTemperatureWithSalt(string primer, double naConcentration = 50)
```

**Behavior:**
1. Calculates base Tm using `CalculateMeltingTemperature`
2. Adds salt correction factor
3. Rounds to 1 decimal place

**Parameters:**
- `primer`: DNA sequence
- `naConcentration`: Na+ concentration in mM (default: 50)

## Example Calculations

### Wallace Rule Example

Primer: `GCGCGCGC` (8 bp)
- A + T = 0
- G + C = 8
- Tm = 2×0 + 4×8 = **32°C**

### Marmur-Doty Example

Primer: `ACGTACGTACGTACGTACGT` (20 bp)
- GC count = 10
- Tm = 64.9 + 41×(10 - 16.4)/20
- Tm = 64.9 + 41×(-6.4)/20
- Tm = 64.9 - 13.12 = **51.78°C**

### Salt Correction Example

Base Tm: 51.78°C at 50 mM Na+
- Correction = 16.6 × log10(50/1000)
- Correction = 16.6 × log10(0.05)
- Correction = 16.6 × (-1.301) = -21.6°C
- Final Tm = 51.78 - 21.6 = **30.2°C**

## Typical Tm Ranges for PCR Primers

| Application | Recommended Tm | Notes |
|-------------|----------------|-------|
| Standard PCR | 55-65°C | Optimal annealing |
| High-fidelity PCR | 60-72°C | Higher specificity |
| Colony PCR | 50-55°C | Lower stringency |
| Real-time PCR | 58-62°C | Narrow range preferred |

## Deviations from Literature

### ASSUMPTION: Threshold at 14 bp

The implementation uses < 14 bp as the threshold for Wallace rule. Some literature suggests 17-20 bp. This threshold is reasonable for simple Tm estimation but may differ from other tools.

### Implementation Note: Non-ACGT Characters

The implementation counts only A, C, G, T bases. Other characters (N, R, Y, etc.) are excluded from the count, which may affect Tm calculation for degenerate primers.

## References

1. Marmur, J. & Doty, P. (1962). Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature. J Mol Biol 5:109-118.

2. SantaLucia, J. Jr. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proc Natl Acad Sci USA 95:1460-1465.

3. Owczarzy, R. et al. (2004). Effects of sodium ions on DNA duplex oligomers: improved predictions of melting temperatures. Biochemistry 43:3537-3554.

4. Wikipedia: Nucleic acid thermodynamics. https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics

5. Wikipedia: DNA melting. https://en.wikipedia.org/wiki/DNA_melting
