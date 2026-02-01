# Microsatellite Detection (Short Tandem Repeats)

**Test Unit:** REP-STR-001
**Last Updated:** 2026-01-22
**Status:** Active

---

## 1. Definition

Microsatellites (also known as Short Tandem Repeats, STRs, or Simple Sequence Repeats, SSRs) are tracts of tandemly repeated DNA motifs ranging in length from 1 to 6 nucleotides, typically repeated 5–50 times.

### Terminology
- **Mononucleotide repeat:** 1 bp unit (e.g., AAAAAA)
- **Dinucleotide repeat:** 2 bp unit (e.g., CACACACA)
- **Trinucleotide repeat:** 3 bp unit (e.g., CAGCAGCAG)
- **Tetranucleotide repeat:** 4 bp unit (e.g., GATAGATAGATAGATA)
- **Pentanucleotide repeat:** 5 bp unit
- **Hexanucleotide repeat:** 6 bp unit

### Sources
- Wikipedia: [Microsatellite](https://en.wikipedia.org/wiki/Microsatellite)
- Richard GF et al. (2008) "Comparative genomics and molecular dynamics of DNA repeats in eukaryotes" MMBR
- Tóth G et al. (2000) "Microsatellites in different eukaryotic genomes: survey and analysis" Genome Res.

---

## 2. Biological Significance

### Medical Relevance
Trinucleotide repeat expansions are associated with over 30 genetic disorders:

| Disease | Gene | Repeat | Normal | Pathogenic |
|---------|------|--------|--------|------------|
| Huntington's disease | HTT | CAG | 6–35 | 36–250 |
| Fragile X syndrome | FMR1 | CGG | 6–53 | 230+ |
| Friedreich's ataxia | FXN | GAA | 7–34 | 100+ |
| Myotonic dystrophy 1 | DMPK | CTG | 5–34 | 50+ |

**Source:** Wikipedia: [Trinucleotide repeat disorder](https://en.wikipedia.org/wiki/Trinucleotide_repeat_disorder)

### Forensic Applications
- STRs are used for DNA profiling (genetic fingerprinting)
- Forensic markers use tetra- and pentanucleotide repeats (higher accuracy, less PCR stutter)
- FBI CODIS uses 13 core STR loci

---

## 3. Algorithm Description

### Input Parameters
- `sequence`: DNA sequence to search (DnaSequence or string)
- `minUnitLength`: Minimum repeat unit length (default: 1)
- `maxUnitLength`: Maximum repeat unit length (default: 6)
- `minRepeats`: Minimum number of consecutive repeats to report (default: 3)

### Output
Collection of `MicrosatelliteResult` containing:
- `Position`: 0-based start position
- `RepeatUnit`: The repeated motif
- `RepeatCount`: Number of consecutive repeats
- `TotalLength`: Total length in base pairs
- `RepeatType`: Classification (Mono-, Di-, Tri-, etc.)

### Algorithm (Implementation)

```
1. For each unit length from minUnitLength to maxUnitLength:
   a. For each position i in sequence where a valid repeat could start:
      i.   Extract candidate unit of length unitLen at position i
      ii.  Skip if unit is redundant (composed of smaller repeating pattern)
      iii. Count consecutive occurrences of the unit
      iv.  If count >= minRepeats and not overlapping with already reported:
           - Report the microsatellite
```

### Complexity
- **Time:** O(n × U × R) where n = sequence length, U = maxUnitLength, R = average repeat count
- **Space:** O(k) where k = number of results

---

## 4. Implementation Details

### Redundant Unit Detection
The implementation filters out redundant units where a larger pattern is composed of smaller repeating units:
- "ATAT" is redundant because it equals "AT" repeated twice
- "CAGCAG" is redundant because it equals "CAG" repeated twice

This prevents duplicate reporting of the same repeat with different unit sizes.

### Overlap Policy
Results are non-overlapping. When a shorter repeat would be contained within a longer one already reported, the shorter is skipped.

### Case Handling
- Input is normalized to uppercase before processing
- Both DnaSequence and raw string inputs are supported

---

## 5. Edge Cases

| Case | Input | Expected Behavior |
|------|-------|-------------------|
| Empty sequence | "" | Returns empty |
| No repeats found | "ACGT" with minRepeats=3 | Returns empty |
| Entire sequence is one repeat | "AAAAAA" | Single result with correct count |
| Boundary: minRepeats=2 | (any) | ArgumentOutOfRangeException |
| Null DnaSequence | null | ArgumentNullException |
| maxUnitLength < minUnitLength | - | ArgumentOutOfRangeException |

---

## 6. Validation Criteria

### Invariants
1. All results have `RepeatCount >= minRepeats`
2. All results have `minUnitLength <= RepeatUnit.Length <= maxUnitLength`
3. `TotalLength == RepeatUnit.Length × RepeatCount`
4. `FullSequence == RepeatUnit` repeated `RepeatCount` times
5. Results are non-overlapping (or overlap policy is documented)
6. `RepeatType` matches the unit length

### Property-Based Tests
- Repeat at any returned position should match the reported RepeatUnit
- Positions should be within valid range [0, sequence.Length - TotalLength]

---

## 7. References

1. Wikipedia: Microsatellite - https://en.wikipedia.org/wiki/Microsatellite
2. Wikipedia: Trinucleotide repeat disorder - https://en.wikipedia.org/wiki/Trinucleotide_repeat_disorder
3. Richard GF, Kerrest A, Dujon B (2008). "Comparative genomics and molecular dynamics of DNA repeats in eukaryotes". Microbiology and Molecular Biology Reviews. 72(4):686-727.
4. Tóth G, Gáspári Z, Jurka J (2000). "Microsatellites in different eukaryotic genomes: survey and analysis". Genome Research. 10(7):967-981.
5. Brinkmann B et al. (1998). "Mutation rate in human microsatellites". American Journal of Human Genetics.
