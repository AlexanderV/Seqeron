# CHROM-ANEU-001: Aneuploidy Detection - Evidence Document

**Test Unit ID:** CHROM-ANEU-001  
**Algorithm Group:** Chromosome Analysis  
**Date:** 2026-02-01  
**Status:** Complete

---

## 1. Authoritative Sources

### 1.1 Primary Sources

| Source | Type | URL/Reference |
|--------|------|---------------|
| Wikipedia (Aneuploidy) | Encyclopedia | https://en.wikipedia.org/wiki/Aneuploidy |
| Wikipedia (Copy Number Variation) | Encyclopedia | https://en.wikipedia.org/wiki/Copy_number_variation |
| Griffiths et al. (2000) | Textbook | "An Introduction to Genetic Analysis" 7th ed., W.H. Freeman |
| Santaguida & Amon (2015) | Review | doi:10.1038/nrm4025 |
| McCarroll & Altshuler (2007) | Review | doi:10.1038/ng2080 |

### 1.2 Key Definitions from Sources

**Aneuploidy** (Wikipedia):
> "The presence of an abnormal number of chromosomes in a cell, for example a human somatic cell having 45 or 47 chromosomes instead of the usual 46."

**Copy Number Terminology** (Wikipedia):
- **Nullisomy**: 0 copies (lethal for autosomes)
- **Monosomy**: 1 copy (e.g., Turner syndrome 45,X)
- **Disomy**: 2 copies (normal diploid state)
- **Trisomy**: 3 copies (e.g., Down syndrome, chr21)
- **Tetrasomy**: 4 copies
- **Pentasomy**: 5 copies

---

## 2. Algorithm Behavior

### 2.1 Copy Number Detection from Read Depth

The implementation uses depth ratio for copy number estimation:

```
logRatio = log2(observedDepth / medianDepth)
copyNumber = round(2^logRatio × 2)
```

**Depth to Copy Number Mapping:**
| Depth Ratio | Log2 Ratio | Estimated CN |
|-------------|------------|--------------|
| 0.0 | -∞ | 0 (Nullisomy) |
| 0.5 | -1.0 | 1 (Monosomy) |
| 1.0 | 0.0 | 2 (Disomy/Normal) |
| 1.5 | 0.58 | 3 (Trisomy) |
| 2.0 | 1.0 | 4 (Tetrasomy) |

### 2.2 Whole Chromosome Aneuploidy Classification

Classification requires dominant copy number across chromosome (default ≥80% of bins):
- CN=0 → Nullisomy
- CN=1 → Monosomy
- CN=2 → Normal (no aneuploidy)
- CN=3 → Trisomy
- CN=4 → Tetrasomy
- CN>4 → "Copy number = N"

---

## 3. Test Cases from Sources

### 3.1 Documented Clinical Examples

| Condition | Chromosome | Copy Number | Source |
|-----------|------------|-------------|--------|
| Down syndrome | chr21 | 3 (Trisomy) | Wikipedia |
| Edwards syndrome | chr18 | 3 (Trisomy) | Wikipedia |
| Patau syndrome | chr13 | 3 (Trisomy) | Wikipedia |
| Turner syndrome | chrX | 1 (Monosomy) | Wikipedia |
| Klinefelter syndrome | chrX | 3 (XXY) | Wikipedia |

### 3.2 Detection Thresholds

From implementation analysis:
- Copy number clamped to range [0, 10]
- Confidence = 1 - min(1, |expected - observed|)
  - where expected = copyNumber/2, observed = 2^logRatio

---

## 4. Edge Cases

### 4.1 From Sources

1. **Mosaicism**: Variable copy number across cells (Wikipedia)
   - Implementation: Handled by minFraction parameter
   
2. **Sex chromosomes**: Different baseline (Wikipedia)
   - Males: 1 copy of X, 1 copy of Y
   - Implementation: Does not special-case sex chromosomes

### 4.2 Implementation-Specific

1. **Empty input**: Returns empty enumerable
2. **Zero or negative median depth**: Returns empty (prevents division by zero)
3. **Multiple chromosomes**: Groups by chromosome name
4. **Binning**: Aggregates depth values by position/binSize

---

## 5. Invariants

### 5.1 DetectAneuploidy

| Invariant | Source |
|-----------|--------|
| CopyNumber ∈ [0, 10] | Implementation (clamped) |
| Confidence ∈ [0, 1] | Implementation |
| LogRatio = log2(depth/medianDepth) | Mathematical definition |
| Output bins are ordered | Implementation (OrderBy) |

### 5.2 IdentifyWholeChromosomeAneuploidy

| Invariant | Source |
|-----------|--------|
| Only returns chromosomes with CN ≠ 2 | Definition of aneuploidy |
| Requires minFraction of consistent CN | Implementation |
| Type names match standard terminology | Wikipedia |

---

## 6. Fallback Strategy

Not applicable - sufficient authoritative sources found.

---

## 7. Open Questions

1. **Sex chromosome handling**: Current implementation treats X/Y same as autosomes. 
   - For males, monosomic X is normal, not aneuploidy.
   - **Decision**: Document as limitation; do not change behavior for this test unit.

2. **Partial aneuploidy**: Wikipedia mentions partial monosomy/trisomy from translocations.
   - Current implementation detects regional CN changes.
   - Whole chromosome classification requires segment-level CN data.

---

## 8. References

1. Wikipedia contributors. "Aneuploidy." Wikipedia, The Free Encyclopedia. 25 Nov 2025.
2. Wikipedia contributors. "Copy number variation." Wikipedia, The Free Encyclopedia. 12 Dec 2025.
3. Griffiths AJF, et al. An Introduction to Genetic Analysis. 7th ed. New York: W.H. Freeman; 2000.
4. Santaguida S, Amon A. Short- and long-term effects of chromosome mis-segregation and aneuploidy. Nat Rev Mol Cell Biol. 2015;16(8):473-485.
5. McCarroll SA, Altshuler DM. Copy-number variation and association studies of human disease. Nat Genet. 2007;39(7 Suppl):S37-42.
