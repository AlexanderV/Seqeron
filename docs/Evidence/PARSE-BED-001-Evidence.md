# PARSE-BED-001: BED File Parsing — Evidence Document

**Test Unit ID:** PARSE-BED-001  
**Created:** 2026-02-05  
**Algorithm Group:** FileIO  
**Status:** Complete

---

## 1. Primary Sources

### 1.1 UCSC Genome Browser BED Format Specification (Authoritative)
**Source:** https://genome.ucsc.edu/FAQ/FAQformat.html#format1

#### Key Specifications:
1. **Columns**: BED has 3 required fields and 9 optional fields (12 total)
2. **Required Fields (BED3)**:
   - `chrom`: Chromosome/scaffold name (e.g., chr3, chrY)
   - `chromStart`: Start position (0-based)
   - `chromEnd`: End position (not included in display)
3. **Optional Fields (BED4-BED12)**:
   - `name`: Feature name
   - `score`: Score 0-1000 (used for grayscale display)
   - `strand`: "+" or "-" or "."
   - `thickStart`: Start of thick region (e.g., start codon)
   - `thickEnd`: End of thick region (e.g., stop codon)
   - `itemRgb`: RGB color (R,G,B format)
   - `blockCount`: Number of blocks (exons)
   - `blockSizes`: Comma-separated list of block sizes
   - `blockStarts`: Comma-separated list of block starts (relative to chromStart)

#### Coordinate System (Critical):
- **chromStart**: 0-based (first base is 0)
- **chromEnd**: 1-based, non-inclusive
- Example: "chr1:1-100" in browser = chromStart=0, chromEnd=100 spans bases 0-99
- **Zero-length features**: chromStart == chromEnd (used for insertions)

#### Block Constraints (BED12):
- First blockStart must be 0
- Final blockStart + final blockSize must equal chromEnd - chromStart
- Blocks may not overlap

#### Header Lines:
- Lines starting with "track " are track headers
- Lines starting with "browser " are browser control lines
- Lines starting with "#" are comments
- Headers only valid in custom tracks, not for bedToBigBed conversion

### 1.2 Wikipedia BED Format Article
**Source:** https://en.wikipedia.org/wiki/BED_(file_format)

#### Additional Specifications:
- Tab separation recommended for compatibility
- Space separation also allowed
- All rows in a file must have same number of columns
- Formal specification published 2021 by GA4GH (BEDv1.pdf)

### 1.3 BEDTools Manual (Usage Reference)
**Source:** Quinlan & Hall (2010), BEDTools Manual v4

#### Interval Operations:
- **Merge**: Combine overlapping intervals
- **Intersect**: Find overlapping regions between two BED files
- **Subtract**: Remove overlapping regions
- **Complement**: Find regions not covered

---

## 2. Test Data from Sources

### 2.1 UCSC Example - BED12 with blocks
```
chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
```

### 2.2 UCSC Example - BED9 with RGB colors
```
chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
chr7    127475864  127477031  Neg1  0  -  127475864  127477031  0,0,255
```

### 2.3 Edge Case: Zero-length insertion
```
chr1    0    0    insertion_before_first_base
```
Per UCSC: "chromStart=0, chromEnd=0 to represent an insertion before the first nucleotide"

---

## 3. Documented Corner Cases

### 3.1 Coordinate System Issues
| Case | chromStart | chromEnd | Interpretation |
|------|------------|----------|----------------|
| First 100 bases | 0 | 100 | Bases 0-99 (100 bases) |
| Single base at pos 1 | 0 | 1 | Base 0 only |
| Insertion point | 5 | 5 | Zero-length, insertion at position 5 |

### 3.2 Score Clamping
- UCSC: "A score between 0 and 1000"
- Scores outside this range should be clamped

### 3.3 Block Validation (BED12)
- blockCount must match number of items in blockSizes and blockStarts
- First blockStart must be 0
- Blocks must not overlap
- Sum of last blockStart + last blockSize must equal feature length

### 3.4 Strand Values
- Valid: "+", "-", "."
- "." indicates no strand or strand unknown

### 3.5 Delimiter Handling
- Tab-separated preferred
- Space-separated allowed but may cause issues with names containing spaces

---

## 4. Known Failure Modes

1. **Invalid coordinates**: chromStart > chromEnd (except zero-length)
2. **Non-numeric coordinates**: Characters in coordinate fields
3. **Insufficient columns**: Fewer than 3 columns
4. **Block mismatch**: blockCount doesn't match actual block arrays
5. **Block overflow**: blockStarts[n] + blockSizes[n] exceeds feature length
6. **Invalid strand**: Characters other than +, -, .

---

## 5. Established Testing Methodologies

### 5.1 From BEDTools Test Suite
- Round-trip parsing (parse → write → parse)
- Interval arithmetic correctness
- Edge case coverage for zero-length features
- Chromosome name normalization

### 5.2 From UCSC Browser
- Header line skipping
- Track/browser line handling
- Comment line filtering

---

## 6. Implementation-Specific Notes

### Current Implementation (BedParser.cs)
- Supports BED3-BED12 with auto-detection
- Score clamping to 0-1000 range
- Header line skipping (track, browser, #)
- Both tab and space delimiters
- GenomicInterval record for interval operations
- MergeOverlapping, Intersect, Subtract operations
- ExpandIntervals with strand-aware upstream/downstream
- BED12 block operations (ExpandBlocks, GetIntrons, GetTotalBlockLength)
- Statistics calculation
- Coverage depth calculation

---

## 7. References

1. Kent WJ, et al. (2002). "The Human Genome Browser at UCSC". Genome Research. 12(6): 996–1006. doi:10.1101/gr.229102
2. UCSC Genome Browser FAQ - Data File Formats. https://genome.ucsc.edu/FAQ/FAQformat.html
3. Quinlan AR, Hall IM (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26(6): 841-842.
4. GA4GH BED v1.0 Specification (2021). https://samtools.github.io/hts-specs/BEDv1.pdf
