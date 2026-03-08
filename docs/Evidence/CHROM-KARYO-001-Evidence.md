# CHROM-KARYO-001: Karyotype Analysis Evidence

**Test Unit ID:** CHROM-KARYO-001  
**Area:** Chromosome Analysis  
**Date:** 2026-03-08  
**Status:** Complete  

---

## Authoritative Sources

### Primary Sources

1. **Wikipedia - Karyotype**
   - URL: https://en.wikipedia.org/wiki/Karyotype
   - Key concepts:
     - A karyotype is the general appearance of the complete set of chromosomes in a cell
     - Normal human diploid karyotype: 46 chromosomes (22 autosomal pairs + 2 sex chromosomes)
     - Karyotype notation: 46,XX (female) or 46,XY (male)
     - Autosomes numbered 1-22, from largest to smallest
     - Human chromosome groups A-G based on size and centromere position
   - Verified: 2026-03-08

2. **Wikipedia - Ploidy**
   - URL: https://en.wikipedia.org/wiki/Ploidy
   - Key concepts:
     - Ploidy is the number of complete sets of chromosomes in a cell
     - Monoploid (1 set), Diploid (2 sets), Triploid (3 sets), Tetraploid (4 sets), etc.
     - Humans: 2n = 46, n = 23
     - Polyploidy common in plants, rare in animals
     - Euploidy (normal set count) vs Aneuploidy (abnormal individual chromosome count)
   - Verified: 2026-03-08

3. **Wikipedia - Aneuploidy**
   - URL: https://en.wikipedia.org/wiki/Aneuploidy
   - Key concepts:
     - Aneuploidy = abnormal number of individual chromosomes (NOT whole-set changes)
     - Terminology is based on **absolute copy count**:
       - Nullisomy: 0 copies
       - Monosomy: 1 copy (e.g., Turner syndrome 45,X)
       - Disomy: 2 copies (normal for diploid)
       - Trisomy: 3 copies (e.g., Down syndrome: Trisomy 21)
       - Tetrasomy: 4 copies
       - Pentasomy: 5 copies
     - Sex chromosome tetrasomy and pentasomy (XXXX, XXXXY, etc.) documented in humans
     - Most autosomal trisomies are lethal; survivable: Trisomy 21, 18, 13
   - Verified: 2026-03-08

---

## Key Algorithms and Invariants

### 1. AnalyzeKaryotype

**Purpose:** Analyze karyotype from chromosome data to detect aneuploidy and classify chromosomes.

**Algorithm Steps (from sources):**
1. Separate sex chromosomes from autosomes
2. Group autosomes by base chromosome name (strip copy suffixes)
3. Count copies of each chromosome group
4. Compare counts against expected ploidy level
5. Label aneuploidy using **standard cytogenetic nomenclature** based on absolute copy count

**Invariants:**
- TotalChromosomes = AutosomeCount + SexChromosomeCount
- TotalGenomeSize = Σ(chromosome lengths)
- MeanChromosomeLength = TotalGenomeSize / TotalChromosomes
- HasAneuploidy = true IFF any chromosome group has count ≠ expectedPloidy

### 2. DetectPloidy

**Purpose:** Detect ploidy level from normalized read depth data.

**Algorithm Steps:**
1. If empty input → return (2, 0) — default diploid, zero confidence
2. Calculate **true median** from sorted depth values (average of two middle elements for even counts)
3. Compute ratio: medianDepth / expectedDiploidDepth
4. Estimate ploidy: round(ratio × 2)
5. Clamp to valid range [1, 8]
6. Calculate confidence: 1.0 − |ratio × 2 − ploidy| × 2

**Invariants:**
- PloidyLevel ∈ [1, 8]
- Confidence ∈ [0, 1]
- Diploid (ploidy=2): ratio ≈ 1.0
- Tetraploid (ploidy=4): ratio ≈ 2.0

---

## Test Datasets

### From Wikipedia / Standard Biology

| Scenario | Input | Expected Output | Source |
|----------|-------|-----------------|--------|
| Normal diploid human | 22 autosome pairs + XX/XY | 46 chromosomes, no aneuploidy | Wikipedia Karyotype |
| Trisomy 21 (Down syndrome) | 3 copies of chr21 | HasAneuploidy=true, "Trisomy" | Wikipedia Aneuploidy |
| Turner syndrome (45,X) | Single X, missing Y | HasAneuploidy=true, "Monosomy" | Wikipedia Aneuploidy |
| Disomy in tetraploid | 2 copies in tetraploid | HasAneuploidy=true, "Disomy" | Wikipedia Aneuploidy |
| Tetrasomy | 4 copies of a chr in diploid | HasAneuploidy=true, "Tetrasomy" | Wikipedia Aneuploidy |
| Pentasomy | 5 copies of a chr in diploid | HasAneuploidy=true, "Pentasomy" | Wikipedia Aneuploidy |
| Diploid depth | ratio ≈ 1.0 | ploidy=2 | Wikipedia Ploidy |
| Tetraploid depth | ratio ≈ 2.0 | ploidy=4 | Wikipedia Ploidy |
| Haploid depth | ratio ≈ 0.5 | ploidy=1 | Wikipedia Ploidy |

---

## Design Decisions

| ID | Decision | Rationale |
|----|----------|-----------|
| DD1 | Empty chromosome input → empty karyotype, no aneuploidy | Graceful degradation for edge case |
| DD2 | Empty depth input → (ploidy=2, confidence=0) | Diploid is most common default; zero confidence signals no data |
| DD3 | Ploidy clamped to [1, 8] | Practical limit; higher ploidy exists in nature (plants, polytene chromosomes up to 1024-ploid per Wikipedia) but is outside typical analysis scope |
| DD4 | Nullisomy (0 copies) is unreachable via `GroupBy` | Architecture detects only chromosomes present in input; absent chromosomes cannot form a group. Term is mapped for completeness |
| DD5 | Disomy (2 copies) is only aneuploidy in non-diploid contexts | In diploid organisms 2 copies is normal (never triggers); in polyploid contexts, 2 copies is correctly labeled Disomy per ISCN |

---

## Deviations and Assumptions

None.
