# GC Skew

**Test Unit ID:** SEQ-GCSKEW-001
**Algorithm Group:** Sequence Composition
**Implementation:** `GcSkewCalculator` class

---

## 1. Definition

GC skew is a measure of strand-specific nucleotide composition asymmetry, defined as:

$$
\text{GC skew} = \frac{G - C}{G + C}
$$

Where G and C are the counts of guanine and cytosine nucleotides in a sequence region.

### Value Range

- **Minimum:** -1 (when G = 0, sequence contains only C among G/C bases)
- **Maximum:** +1 (when C = 0, sequence contains only G among G/C bases)
- **Zero:** When G = C, or when neither G nor C is present

---

## 2. Biological Significance

### Replication Origin Detection

GC skew analysis is a standard bioinformatics method for identifying replication origins and termini in bacterial and archaeal genomes (Lobry 1996, Grigoriev 1998).

**Key observations:**
- The **leading strand** (synthesized continuously during DNA replication) typically shows **positive GC skew** (G > C)
- The **lagging strand** (synthesized discontinuously via Okazaki fragments) typically shows **negative GC skew** (C > G)
- GC skew **changes sign** at the boundaries between the two replichores (replication origin and terminus)

### Mechanism

The asymmetry arises from:
1. **Differential mutational pressure** between leading and lagging strands during replication
2. **Deamination of cytosine** to uracil (→ thymine) occurs more frequently on the single-stranded template during replication
3. The leading strand spends more time single-stranded and accumulates more C→T mutations

---

## 3. Cumulative GC Skew

The cumulative GC skew method (Grigoriev 1998) provides a more robust way to identify replication boundaries:

$$
\text{Cumulative GC skew}(n) = \sum_{i=1}^{n} \text{GC skew}(window_i)
$$

### Interpretation

For a typical bacterial circular chromosome:
- **Global minimum** of cumulative GC skew corresponds to the **origin of replication (oriC)**
- **Global maximum** corresponds to the **terminus of replication (ter)**
- The two extrema are typically separated by approximately half the chromosome length

---

## 4. Implementation Details

### GcSkewCalculator Class

The implementation provides several methods for GC skew analysis:

| Method | Description | Complexity |
|--------|-------------|------------|
| `CalculateGcSkew(sequence)` | Overall GC skew for entire sequence | O(n) |
| `CalculateWindowedGcSkew(sequence, windowSize, stepSize)` | Sliding window GC skew | O(n) |
| `CalculateCumulativeGcSkew(sequence, windowSize)` | Cumulative GC skew | O(n) |
| `PredictReplicationOrigin(sequence, windowSize)` | Predict ori/ter positions | O(n) |

### Edge Case Handling

| Condition | Behavior |
|-----------|----------|
| Empty sequence | Returns 0 |
| No G or C bases | Returns 0 (division by zero protection) |
| Null input | Throws `ArgumentNullException` |
| Invalid window/step size | Throws `ArgumentOutOfRangeException` |

### Window Position Reporting

Positions are reported at the **center** of each window:
```
Position = WindowStart + WindowSize / 2
```

---

## 5. Related Metrics

### AT Skew

AT skew follows the same formula for adenine and thymine:

$$
\text{AT skew} = \frac{A - T}{A + T}
$$

The `GcSkewCalculator` also provides `CalculateAtSkew()` methods.

### Relationship to GC Content

- **GC Content** = (G + C) / total × 100%
- **GC Skew** measures the relative abundance of G vs C
- High GC content with zero skew indicates equal G and C

---

## 6. References

1. Lobry, J.R. (1996). "Asymmetric substitution patterns in the two DNA strands of bacteria." *Molecular Biology and Evolution*, 13(5):660-665. [PMID: 8676740]

2. Grigoriev, A. (1998). "Analyzing genomes with cumulative skew diagrams." *Nucleic Acids Research*, 26(10):2286-2290. [DOI: 10.1093/nar/26.10.2286]

3. Tillier, E.R. & Collins, R.A. (2000). "The contributions of replication orientation, gene direction, and signal sequences to base-composition asymmetries in bacterial genomes." *Journal of Molecular Evolution*, 50:249-257.

4. Wikipedia contributors. "GC skew." *Wikipedia, The Free Encyclopedia*. (Accessed 2026-01-22)

---

## 7. Implementation Notes

### Current Implementation

The `GcSkewCalculator` in `Seqeron.Genomics` provides:
- Core skew calculation with O(n) complexity
- Sliding window analysis with configurable window and step sizes
- Cumulative skew for origin/terminus detection
- Combined GC content and skew analysis

### Deviations from Theory

None - the implementation follows the standard formula exactly.

### Known Limitations

1. Origin/terminus prediction assumes a single circular chromosome with bi-directional replication
2. Significance threshold for predictions is heuristic (amplitude > 1% of point count)
3. Does not account for horizontal gene transfer or recent inversions that may distort skew patterns
