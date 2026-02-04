# PARSE-FASTQ-001 Evidence

## Test Unit
**ID:** PARSE-FASTQ-001  
**Area:** FileIO  
**Algorithm:** FASTQ Parsing  

---

## Authoritative Sources

### Primary Sources

1. **Wikipedia - FASTQ format**
   - URL: https://en.wikipedia.org/wiki/FASTQ_format
   - Access Date: 2026-02-05
   - Key Information:
     - Format specification: 4 lines per record (@header, sequence, +, quality)
     - Quality encoding: Phred+33 (Sanger/Illumina 1.8+) and Phred+64 (Illumina 1.3-1.7)
     - Quality formula: Q = -10 × log₁₀(p), where p = error probability
     - ASCII ranges: Phred+33 uses 33-126, Phred+64 uses 64-126
     - Detection heuristic: chars < '@' → Phred+33; chars > 'I' → Phred+64

2. **Cock et al. (2009) - The Sanger FASTQ file format**
   - DOI: 10.1093/nar/gkp1137
   - Citation: Cock PJA, Fields CJ, Goto N, et al. Nucleic Acids Research. 2009;38(6):1767-1771
   - Key Information:
     - Authoritative definition of Sanger FASTQ format
     - Quality score encoding variations
     - Historical context of format evolution

3. **NCBI Sequence Read Archive - File Format Guide**
   - URL: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/
   - Key Information:
     - Official submission requirements
     - Paired-end FASTQ conventions (/1, /2 suffixes)
     - Quality encoding requirements for SRA

---

## Format Specification (from Wikipedia)

### Structure
```
@<identifier> [description]
<sequence>
+[optional identifier]
<quality>
```

### Quality Encodings

| Format | Offset | ASCII Range | Q Range | Era |
|--------|--------|-------------|---------|-----|
| Sanger | 33 | 33-126 | 0-93 | Current standard |
| Illumina 1.8+ | 33 | 33-126 | 0-41 | Current Illumina |
| Illumina 1.3-1.7 | 64 | 64-126 | 0-62 | Legacy |
| Solexa | 64 | 59-126 | -5-62 | Obsolete |

### Quality Score Mathematics
- **Phred score formula:** Q = -10 × log₁₀(p)
- **Error probability:** p = 10^(-Q/10)
- **Common values:**
  - Q10 → 10% error rate (1 in 10)
  - Q20 → 1% error rate (1 in 100)
  - Q30 → 0.1% error rate (1 in 1000)
  - Q40 → 0.01% error rate (1 in 10,000)

### Auto-Detection Heuristic
- Characters below '@' (ASCII 64) indicate Phred+33
- Characters above 'I' (ASCII 73) indicate Phred+64
- Ambiguous range '@'-'I' requires context or defaults to Phred+33

---

## Edge Cases from Sources

### Documented Edge Cases
1. **Multi-line sequences** (Wikipedia): Legacy Sanger files may split sequences across lines
2. **@ in quality string** (Wikipedia): Makes parsing ambiguous in multi-line files
3. **+ in sequence** (Wikipedia): Unusual but allowed in sequence data
4. **Empty records** (Implementation): Parser should skip blank lines gracefully
5. **Encoding detection failure** (Wikipedia): Default to Phred+33 when ambiguous

### Illumina-Specific Conventions
- Header format: `@INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y READ:FILTER:CONTROL:INDEX`
- Paired-end indicators: `/1`, `/2` or `1:`, `2:` in newer format
- Interleaved format: alternating R1/R2 records

---

## Test Dataset Patterns

### From Wikipedia Examples
```fastq
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

### Quality Score Boundary Values
- `!` (ASCII 33) = Q0 in Phred+33
- `I` (ASCII 73) = Q40 in Phred+33
- `@` (ASCII 64) = Q0 in Phred+64
- `h` (ASCII 104) = Q40 in Phred+64

---

## Implementation Notes

### Current Implementation (FastqParser.cs)
- Supports Phred+33 and Phred+64 encodings
- Auto-detection based on character range analysis
- Multi-line sequence/quality support via loop
- Paired-end support: interleaving and splitting
- Statistics: total reads, bases, GC content, Q20/Q30 percentages
- Filtering: by quality threshold, by length range
- Trimming: quality-based end trimming, adapter removal

### Invariants
1. Sequence length == Quality string length (per record)
2. Quality scores ∈ [0, 93] for Phred+33, [0, 62] for Phred+64
3. Encoding detection is deterministic for given input
4. Round-trip: Parse → Write → Parse yields equivalent data

---

## References

1. Wikipedia contributors. "FASTQ format." Wikipedia, The Free Encyclopedia. https://en.wikipedia.org/wiki/FASTQ_format
2. Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM. (2009). The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic Acids Research. 38(6):1767-1771.
3. NCBI. SRA File Format Guide. https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/
