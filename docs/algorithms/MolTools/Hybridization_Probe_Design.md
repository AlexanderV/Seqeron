# Hybridization Probe Design

## Overview

Hybridization probe design is a computational process for designing oligonucleotide probes that can detect specific nucleic acid sequences through complementary base pairing. Probes are used in various applications including FISH (Fluorescence In Situ Hybridization), DNA microarrays, Northern blots, and qPCR.

## Algorithm Description

### Core Algorithm: `DesignProbes`

The `DesignProbes` method generates hybridization probes for a target sequence:

1. **Validate Input**
   - Return empty if sequence is null, empty, or shorter than MinLength

2. **Generate Candidates**
   - Scan all possible probe positions and lengths within constraints
   - For each candidate: evaluate using `EvaluateProbe`
   - Collect valid probes (score > 0)

3. **Rank and Select**
   - Sort candidates by score (descending)
   - Return top `maxProbes` candidates

4. **Return Result**
   - Each probe includes: sequence, start/end positions, Tm, GC content, score, type, warnings

### Complexity

- **Time**: O(n × m) where n = sequence length, m = length range (optimized with prefix sums)
- **Space**: O(k) where k = number of valid candidates

### Optimization: GC Prefix Sums

The implementation precomputes GC content prefix sums for O(1) GC lookup at any position, reducing the inner loop complexity:

```csharp
// Precompute: O(n)
int[] gcPrefixSum = new int[n + 1];
for (int i = 0; i < n; i++)
    gcPrefixSum[i + 1] = gcPrefixSum[i] + (IsGC(seq[i]) ? 1 : 0);

// Query: O(1)
double gc = (double)(gcPrefixSum[end] - gcPrefixSum[start]) / length;
```

### Optimization: Suffix Tree for Specificity

For genome-wide probe design, the suffix tree enables O(m) uniqueness checking:

```csharp
// Build once: O(n) for entire genome
var genomeIndex = SuffixTree.Build(genome);

// Check specificity: O(m) per probe
var positions = genomeIndex.FindAllOccurrences(probeSequence);
bool isUnique = positions.Count == 1;
```

This is crucial for FISH probes where specificity is paramount.

## Quality Criteria

### Standard Parameters (Application-Specific)

| Application | Length (bp) | Tm (°C) | GC Content (%) | Source |
|-------------|-------------|---------|----------------|--------|
| Microarray | 50-70 | 75-85 | 40-60 | Wikipedia (DNA microarray), Implementation |
| FISH | 200-500 | 70-90 | 35-65 | Wikipedia (FISH), Implementation |
| Northern Blot | 100-300 | 65-80 | 40-60 | Implementation |
| qPCR | 20-30 | 68-72 | 40-60 | Standard practice |
| Southern Blot | 150-500 | 65-75 | 35-65 | Implementation |

### Probe Quality Factors

| Factor | Threshold | Penalty | Source |
|--------|-----------|---------|--------|
| GC Content | Outside min/max range | -0.3 | Wikipedia (Nucleic acid thermodynamics) |
| Tm | Outside min/max range | -0.3 | Wikipedia (Nucleic acid thermodynamics) |
| Homopolymer Run | > MaxHomopolymer | -0.2 | General practice |
| Self-Complementarity | > MaxSelfComplementarity | -0.2 | Implementation |
| Secondary Structure | Potential hairpin | -0.15 | Wikipedia (Hybridization probe) |
| Simple Repeats | Di/trinucleotide repeats | -0.1 | General practice |

## Thermodynamic Parameters

### Melting Temperature (Tm) Calculation

The implementation uses two methods based on probe length:

**Short Oligos (≤14 bp): Wallace Rule**
```
Tm = 2×(A+T) + 4×(G+C)
```
Source: Wallace et al. (1979)

**Longer Probes: Salt-Adjusted Formula**
```
Tm = 81.5 + 16.6×log10([Na+]) + 41×(GC fraction) - 675/length
```
Source: SantaLucia (1998), Wikipedia (Nucleic acid thermodynamics)

### GC Content

```
GC = (count(G) + count(C)) / length × 100%
```

Higher GC content increases Tm due to three hydrogen bonds (G≡C) vs two (A=T).
Source: Wikipedia (Nucleic acid thermodynamics)

## Scoring Algorithm

The `EvaluateProbe` method assigns a quality score (0.0 to 1.0):

```
score = 1.0
if (GC outside range): score -= 0.3
if (Tm outside range): score -= 0.3
if (homopolymer > max): score -= 0.2
if (selfComplementarity > max): score -= 0.2
if (hasSecondaryStructure): score -= 0.15
if (hasSimpleRepeats): score -= 0.1
if (starts/ends with G/C): score -= 0.02 each
if (score ≤ 0): reject probe
```

## Secondary Methods

### `DesignTilingProbes`

Designs overlapping probes for complete sequence coverage:
- Fixed probe length with configurable overlap
- Reports coverage statistics (covered bases, mean Tm, Tm range)
- Includes suboptimal probes for full coverage when needed

### `ScoreProbe` (via `EvaluateProbe`)

Evaluates a single probe sequence against parameters and returns:
- Tm, GC content, score
- List of quality warnings

## Self-Complementarity Detection

Measures potential for probe to form intramolecular structures:
1. Compute reverse complement
2. Align with original sequence
3. Count matching positions
4. Return fraction: matches / length

Values > 0.3-0.4 indicate significant self-complementarity risk.
Source: General molecular biology practice

## Secondary Structure Detection

Scans for potential hairpin formation:
1. Search for inverted repeat pairs with ≥4 bp stems
2. Require ≥3 bp loop between stems
3. Report if stem has ≥80% complementarity

Source: Wikipedia (Hybridization probe), General practice

## Implementation Notes

### Probe Types

| Type | Description |
|------|-------------|
| Standard | General hybridization probe |
| Tiling | Overlapping probes for coverage |
| Antisense | Reverse complement for mRNA detection |
| LNA | Locked Nucleic Acid (enhanced binding) |
| MolecularBeacon | Hairpin probe for real-time detection |

### Input Normalization

- All sequences are converted to uppercase before processing
- Supports DNA bases (A, C, G, T)

## References

1. Wikipedia. "Nucleic acid thermodynamics." https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
2. Wikipedia. "Hybridization probe." https://en.wikipedia.org/wiki/Hybridization_probe
3. Wikipedia. "Fluorescence in situ hybridization." https://en.wikipedia.org/wiki/Fluorescence_in_situ_hybridization
4. Wikipedia. "DNA microarray." https://en.wikipedia.org/wiki/DNA_microarray
5. SantaLucia, J. (1998). "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." PNAS 95(4):1460-5.
6. Breslauer, K.J. et al. (1986). "Predicting DNA Duplex Stability from the Base Sequence." PNAS 83:3746-3750.
