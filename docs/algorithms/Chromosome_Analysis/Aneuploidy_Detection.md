# Aneuploidy Detection

**Algorithm Group:** Chromosome Analysis  
**Test Unit ID:** CHROM-ANEU-001  
**Implementation:** `ChromosomeAnalyzer.DetectAneuploidy`, `ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy`

---

## 1. Overview

Aneuploidy is the presence of an abnormal number of chromosomes in a cell. In humans, the normal diploid complement is 46 chromosomes (2n). Aneuploidy involves gains or losses of individual chromosomes rather than complete sets.

**Sources:** Wikipedia (Aneuploidy), Griffiths et al. (2000)

---

## 2. Terminology

| Term | Definition | Copy Number |
|------|------------|-------------|
| Nullisomy | Complete absence of a chromosome pair | 0 |
| Monosomy | Single copy instead of pair | 1 |
| Disomy | Normal two copies (diploid) | 2 |
| Trisomy | Three copies | 3 |
| Tetrasomy | Four copies | 4 |
| Pentasomy | Five copies | 5 |

**Source:** Wikipedia (Aneuploidy)

---

## 3. Clinical Significance

Common human aneuploidies compatible with live birth:

| Condition | Karyotype | Description |
|-----------|-----------|-------------|
| Down syndrome | 47,XX,+21 or 47,XY,+21 | Trisomy 21 |
| Edwards syndrome | 47,XX,+18 or 47,XY,+18 | Trisomy 18 |
| Patau syndrome | 47,XX,+13 or 47,XY,+13 | Trisomy 13 |
| Turner syndrome | 45,X | Monosomy X (females) |
| Klinefelter syndrome | 47,XXY | Extra X (males) |

**Source:** Wikipedia (Aneuploidy)

---

## 4. Detection Algorithm

### 4.1 Copy Number from Read Depth

The implementation infers copy number from sequencing read depth:

```
logRatio = log₂(observedDepth / medianDepth)
copyNumber = round(2^logRatio × 2)
```

**Rationale:** In diploid cells, depth is proportional to copy number. A region with 3 copies (trisomy) has ~1.5× the depth of a 2-copy region.

### 4.2 Depth-to-Copy-Number Mapping

| Observed/Expected Ratio | Log₂ Ratio | Copy Number |
|------------------------|------------|-------------|
| 0.0 | -∞ | 0 |
| 0.5 | -1.0 | 1 |
| 1.0 | 0.0 | 2 |
| 1.5 | +0.58 | 3 |
| 2.0 | +1.0 | 4 |

### 4.3 Confidence Calculation

```
expected = copyNumber / 2.0
observed = 2^logRatio
confidence = 1.0 - min(1.0, |expected - observed|)
```

---

## 5. Implementation Details

### 5.1 DetectAneuploidy

**Signature:**
```csharp
public static IEnumerable<CopyNumberState> DetectAneuploidy(
    IEnumerable<(string Chromosome, int Position, double Depth)> depthData,
    double medianDepth,
    int binSize = 1000000)
```

**Parameters:**
- `depthData`: Read depth at genomic positions
- `medianDepth`: Expected depth for diploid regions
- `binSize`: Aggregation window size (default 1 Mb)

**Returns:** `CopyNumberState` records with chromosome, position range, copy number, log ratio, and confidence.

**Complexity:** O(n) where n = number of depth data points

### 5.2 IdentifyWholeChromosomeAneuploidy

**Signature:**
```csharp
public static IEnumerable<(string Chromosome, int CopyNumber, string Type)> 
    IdentifyWholeChromosomeAneuploidy(
        IEnumerable<CopyNumberState> copyNumberStates,
        double minFraction = 0.8)
```

**Parameters:**
- `copyNumberStates`: Per-bin copy number states
- `minFraction`: Minimum fraction of bins with consistent CN (default 80%)

**Returns:** Chromosomes with aneuploidy and their classification.

**Complexity:** O(n) where n = number of copy number states

---

## 6. Constraints and Limitations

1. **Copy number range:** Clamped to [0, 10]
2. **Sex chromosomes:** Not specially handled; males show monosomic X/Y
3. **Mosaicism:** Detected via minFraction < 1.0; mixed CN per chromosome
4. **Partial aneuploidy:** Regional CN changes detected but classified at bin level

---

## 7. Data Structures

```csharp
public readonly record struct CopyNumberState(
    string Chromosome,
    int Start,
    int End,
    int CopyNumber,
    double LogRatio,
    double Confidence);
```

---

## 8. References

1. Wikipedia. "Aneuploidy." https://en.wikipedia.org/wiki/Aneuploidy
2. Wikipedia. "Copy number variation." https://en.wikipedia.org/wiki/Copy_number_variation
3. Griffiths AJF, et al. An Introduction to Genetic Analysis. 7th ed. 2000.
