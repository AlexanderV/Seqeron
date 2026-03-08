# CHROM-TELO-001 Evidence Document

## Test Unit
**ID:** CHROM-TELO-001
**Title:** Telomere Analysis
**Area:** Chromosome
**Date:** 2026-02-01

---

## Sources

### Primary Sources

1. **Wikipedia - Telomere**
   - URL: https://en.wikipedia.org/wiki/Telomere
   - Last accessed: 2026-02-01
   - Authority: Encyclopedic, peer-reviewed references

2. **Meyne et al. (1989)** - Conservation of the human telomere sequence (TTAGGG)n among vertebrates
   - DOI: 10.1073/pnas.86.18.7049
   - PMID: 2780561
   - Key finding: TTAGGG is conserved across all vertebrates

3. **Cawthon (2002)** - Telomere measurement by quantitative PCR
   - DOI: 10.1093/nar/30.10.e47
   - PMID: 12000852
   - Key finding: T/S ratio is proportional to average telomere length

4. **Blackburn & Gall (1978)** - Nobel Prize-winning work on telomere structure
   - DOI: 10.1016/0022-2836(78)90294-2
   - PMID: 642006

---

## Key Biological Facts

### Telomere Structure
- **Vertebrate telomere repeat:** TTAGGG (5' to 3' toward chromosome end) — Wikipedia; Meyne et al. (1989)
- **Reverse complement (5' end):** CCCTAA
- **Repeat unit length:** 6 base pairs
- **Human telomere length:** Many kilobases; shortens with age at ~50–100 bp per cell division — Wikipedia
- **Critically short telomeres** trigger DNA damage response and cellular senescence — Wikipedia. Implementation default threshold: 3,000 bp (configurable parameter)

### Telomere Orientation
- **3' end (chromosome terminus):** Contains TTAGGG repeats extending toward the end
- **5' end (chromosome start):** Contains CCCTAA repeats (reverse complement of TTAGGG)
- **3' overhang:** 75–300 bases of single-stranded TTAGGG at very end

### Species Variation (from Wikipedia table)
| Organism | Telomeric Repeat |
|----------|------------------|
| Human, mouse, Xenopus | TTAGGG |
| Arabidopsis thaliana | TTTAGGG |
| Tetrahymena | TTGGGG |
| S. cerevisiae | TGTGGGTGTGGTG (from RNA template) |
| Bombyx mori | TTAGG |

---

## Measurement Methods

### T/S Ratio (qPCR Method)
- **Source:** Cawthon (2002)
- **Principle:** Ratio of telomere repeat copy number (T) to single-copy gene copy number (S)
- **Formula:** Telomere length ∝ T/S ratio
- **Linear relationship:** length = referenceLength × (tsRatio / referenceRatio)

### Reference Values
- Reference telomere length for T/S ratio calculations: configurable (default 7,000 bp)
- T/S ratio of 1.0 corresponds to reference sample length — Cawthon (2002)

---

## Edge Cases (Documented)

### Empty/Null Sequences
- Empty sequence → no telomere detected, marked as critically short
- Null handling → implementation-specific

### Minimum Telomere Length Thresholds
- Configurable detection threshold (default 500 bp)
- Lower thresholds increase sensitivity; higher thresholds reduce false positives

### Repeat Purity
- Biological telomeres show some divergence from perfect repeats
- Implementation uses 70% per-window similarity threshold: for 6 bp repeat, 5/6 bases must match (1 mismatch allowed); for 7 bp repeat (e.g. Arabidopsis TTTAGGG), 5/7 bases must match (2 mismatches allowed)
- Higher purity = younger/healthier telomere

### Critical Length Assessment
- Critically short telomeres trigger DNA damage response and cellular senescence — Wikipedia
- Implementation default critical threshold: 3,000 bp (configurable parameter, not a fixed biological constant)

---

## Implementation Notes

### Algorithm Design
The implementation should:
1. Search for CCCTAA repeats at 5' end (looking from start)
2. Search for TTAGGG repeats at 3' end (looking from end)
3. Allow configurable similarity threshold (default 70%)
4. Track both length and purity

### Test Datasets

| Test Case | Sequence Pattern | Expected Result |
|-----------|-----------------|-----------------|
| 3' telomere | [1000 A's] + [200× TTAGGG] | Has3PrimeTelomere=true, length=1200, purity=1.0 |
| 5' telomere | [200× CCCTAA] + [1000 A's] | Has5PrimeTelomere=true, length=1200, purity=1.0 |
| Both ends | [CCCTAA×200] + [2000 A's] + [TTAGGG×200] | Both detected, length=900 each |
| No telomere | [1000 A's] | Neither detected, lengths=0 |
| Empty | "" | Neither detected, critically short |
| Short telomere | [1000 A's] + [TTAGGG×50] | Detected if min threshold ≤ 300 |
| Divergent repeats | [1000 A's] + [TTAGGA×200] | Has3Prime=true, length=1200, purity=5/6≈0.833 |
| Long telomere | [1000 A's] + [TTAGGG×2000] | Has3Prime=true, length=12000, purity=1.0 |
| SearchLength limited | [1000 A's] + [TTAGGG×200], searchLen=600 | length=600 (truncated by search window) |

---

## T/S Ratio Test Cases

| T/S Ratio | Reference Ratio | Reference Length | Expected Length |
|-----------|-----------------|------------------|-----------------|
| 1.0 | 1.0 | 7000 | 7000 |
| 1.5 | 1.0 | 7000 | 10500 |
| 0.5 | 1.0 | 7000 | 3500 |
| 2.0 | 1.0 | 7000 | 14000 |
| 1.0 | 2.0 | 7000 | 3500 |
| 0.0 | 1.0 | 7000 | 0 |

---

## Invariants

1. **Length non-negative:** TelomereLength ≥ 0
2. **Purity range:** 0 ≤ RepeatPurity ≤ 1
3. **Threshold consistency:** Has5Prime/3Prime=true implies length ≥ minTelomereLength
4. **Critical flag logic:** IsCriticallyShort = (hasTelomere && length < criticalLength) OR empty
5. **T/S ratio linearity:** EstimatedLength = referenceLength × (tsRatio / referenceRatio)
6. **Repeat orientation:** 5' end expects reverse complement, 3' end expects forward repeat

---

## Deviations and Assumptions

None — all test rationale and implementation decisions are verified against cited sources:
- Telomere repeat sequences: Wikipedia telomere table; Meyne et al. (1989)
- T/S ratio proportionality: Cawthon (2002) — r² = 0.677 correlation with Southern blot TRF
- 5'/3' orientation: Wikipedia chromosome structure
- Configurable parameters (criticalLength, minTelomereLength, searchLength, referenceLength) are implementation defaults, not biological constants
