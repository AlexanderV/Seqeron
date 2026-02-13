# Evidence Artifact: EPIGEN-CPG-001

**Test Unit ID:** EPIGEN-CPG-001
**Algorithm:** CpG Site Detection
**Date Collected:** 2026-02-13

---

## Online Sources

### Wikipedia — CpG site

**URL:** https://en.wikipedia.org/wiki/CpG_site
**Accessed:** 2026-02-13
**Authority rank:** 4 (Wikipedia with citations to primary sources)

**Key Extracted Points:**

1. **CpG definition:** CpG is shorthand for 5'—C—phosphate—G—3', i.e., cytosine followed by guanine in the 5'→3' direction of a single-stranded DNA sequence. CpG should not be confused with GpC.
2. **CpG island definition:** Regions with a high frequency of CpG sites. The usual formal definition is a region with at least 200 bp, GC percentage > 50%, and observed-to-expected CpG ratio > 60% (citing Gardiner-Garden & Frommer 1987, ref [13]).
3. **O/E formula (Gardiner-Garden):** Observed = number of CpGs; Expected = (number of C × number of G) / length of sequence (ref [13]).
4. **Alternative O/E formula (Saxonov):** Expected = ((number of C + number of G) / 2)² / length of sequence (ref [14]).
5. **Takai & Jones (2002) stricter criteria:** Regions > 500 bp, GC content > 55%, O/E CpG ratio > 65% — more likely to be "true" CpG islands associated with 5' regions of genes (ref [18]).
6. **CpG island typical length:** 300–3,000 base pairs in mammalian genomes.

---

### Gardiner-Garden M, Frommer M (1987) — Primary source

**URL:** https://doi.org/10.1016/0022-2836(87)90689-9
**Accessed:** 2026-02-13
**Authority rank:** 1 (Peer-reviewed paper)

**Key Extracted Points:**

1. **CpG island criteria:** A region of DNA ≥200 bp with GC content >50% and CpG observed/expected ratio >0.6.
2. **O/E CpG formula:** O/E = CpG_count / ((C_count × G_count) / Length). This is the canonical formula used by UCSC Genome Browser and most bioinformatics tools.
3. **Scanning method:** Sliding window approach to identify CpG island boundaries.

---

### Takai D, Jones PA (2002) — Revised criteria

**URL:** https://doi.org/10.1073/pnas.052410099
**Accessed:** 2026-02-13
**Authority rank:** 1 (Peer-reviewed paper)

**Key Extracted Points:**

1. **Stricter criteria:** ≥500 bp, GC% >55%, O/E >0.65 — reduces false positives from Alu repeats.
2. **Validated on human chromosomes 21 and 22.**

---

### Saxonov S, Berg P, Brutlag DL (2006) — Alternative formula

**URL:** https://doi.org/10.1073/pnas.0510310103
**Accessed:** 2026-02-13
**Authority rank:** 1 (Peer-reviewed paper)

**Key Extracted Points:**

1. **Alternative expected formula:** Expected = ((C + G) / 2)² / L.
2. **Genome-wide analysis:** Distinguished two classes of promoters based on CpG content.

---

## Documented Corner Cases and Failure Modes

### From Wikipedia / Gardiner-Garden & Frommer (1987)

1. **No CpG in sequence:** When a sequence contains no CpG dinucleotides, the O/E ratio is 0.
2. **No C or G in sequence:** When either C or G count is zero, the expected value is 0; the O/E ratio should return 0 (division by zero guard).
3. **Sequence too short:** CpG island detection requires minimum window size (≥200 bp by definition). Sequences shorter than the minimum length cannot contain a CpG island.
4. **GpC vs CpG:** GpC (guanine then cytosine) is NOT a CpG site. Only the 5'→3' C-then-G pattern qualifies.
5. **Case sensitivity:** Biological sequence data may be in mixed case; implementation must handle both.
6. **Adjacent CpG sites:** The sequence "CGCG" contains two CpG sites at positions 0 and 2 — not overlapping, since each CpG is a distinct 2-nucleotide window.
7. **Single nucleotide input:** A sequence of length 1 cannot contain any dinucleotide, so it yields 0 CpG sites.

---

## Test Datasets

### Dataset 1: CGCG repeat (synthetic)

**Source:** Hand-derived from CpG definition

| Parameter | Value |
|-----------|-------|
| Sequence | `CGCGCGCGCGCGCGCGCGCG` (20 bp) |
| CpG count | 10 (at positions 0, 2, 4, 6, 8, 10, 12, 14, 16, 18) |
| C count | 10 |
| G count | 10 |
| Expected | (10 × 10) / 20 = 5.0 |
| O/E ratio | 10 / 5.0 = **2.0** |
| GC content | 100% |

### Dataset 2: AT-only (synthetic)

**Source:** Hand-derived — trivial case

| Parameter | Value |
|-----------|-------|
| Sequence | `AATTAATTAATTAATTAATT` (20 bp) |
| CpG count | 0 |
| O/E ratio | **0.0** |

### Dataset 3: Mixed sequence (synthetic)

**Source:** Hand-derived from CpG definition

| Parameter | Value |
|-----------|-------|
| Sequence | `ACGTCGACG` (9 bp) |
| CpG positions | 1, 4, 7 |
| CpG count | 3 |
| C count | 3 |
| G count | 3 |
| Expected | (3 × 3) / 9 = 1.0 |
| O/E ratio | 3 / 1.0 = **3.0** |

### Dataset 4: Minimal CpG (synthetic)

**Source:** Hand-derived

| Parameter | Value |
|-----------|-------|
| Sequence | `ACGT` (4 bp) |
| CpG positions | 1 |
| CpG count | 1 |
| C count | 1 |
| G count | 1 |
| Expected | (1 × 1) / 4 = 0.25 |
| O/E ratio | 1 / 0.25 = **4.0** |

### Dataset 5: Known CpG island behavior

**Source:** Gardiner-Garden & Frommer (1987) criteria

| Parameter | Value |
|-----------|-------|
| Input | 400 bp of "CGCG" repeats |
| Expected: is CpG island? | Yes — Length ≥200, GC = 100% > 50%, O/E = 2.0 > 0.6 |

---

## Assumptions

None. All algorithm behavior is formally defined by Gardiner-Garden & Frommer (1987) and confirmed by Wikipedia's cited primary sources.

---

## Test Coverage Traceability

All recommendations implemented in `EpigeneticsAnalyzer_CpGDetection_Tests.cs` (25 tests).

| Evidence Point | Test Coverage |
|---------------|---------------|
| CpG = C followed by G in 5'→3' | M1, M2, M3, M6, M7, M18 |
| Gardiner-Garden & Frommer O/E formula | M9, M11, M12 |
| CpG island: ≥200 bp, GC% > 50%, O/E > 0.6 | M15, M16, M17, S3, S4 |
| Null/empty input handling | M4a, M4b, M13, C2 |
| Case insensitivity | M5, S2 |
| Edge cases (single char, no G, minimal CG) | M8, M14, S1, C1 |
| GpC ≠ CpG | M18 |
| AT-only → O/E = 0 | M10 |

---

## References

1. Gardiner-Garden M, Frommer M (1987). "CpG islands in vertebrate genomes." J Mol Biol. 196(2):261–282. https://doi.org/10.1016/0022-2836(87)90689-9
2. Takai D, Jones PA (2002). "Comprehensive analysis of CpG islands in human chromosomes 21 and 22." Proc Natl Acad Sci USA. 99(6):3740–5. https://doi.org/10.1073/pnas.052410099
3. Saxonov S, Berg P, Brutlag DL (2006). "A genome-wide analysis of CpG dinucleotides in the human genome distinguishes two distinct classes of promoters." Proc Natl Acad Sci USA. 103(5):1412–1417. https://doi.org/10.1073/pnas.0510310103
4. Wikipedia. "CpG site." https://en.wikipedia.org/wiki/CpG_site (accessed 2026-02-13). Cites [13] Gardiner-Garden & Frommer (1987), [14] Saxonov et al. (2006), [18] Takai & Jones (2002).

---
