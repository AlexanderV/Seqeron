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
  - Q0 → 100% error rate (p = 1.0)
  - Q10 → 10% error rate (1 in 10)
  - Q20 → 1% error rate (1 in 100)
  - Q30 → 0.1% error rate (1 in 1,000)
  - Q40 → 0.01% error rate (1 in 10,000)
  - Q50 → 0.001% error rate (1 in 100,000)
  - Q60 → 0.0001% error rate (1 in 1,000,000)

### Auto-Detection Heuristic (revised 2026-09)
- Decided per FILE (FastQC: lowest quality character over all reads), not per record
- Lowest character below '@' (ASCII 64) proves Phred+33 (outside the Phred+64 range)
- Otherwise a character above 'J' (ASCII 74 = Q41, Illumina 1.8+ ceiling) infers Phred+64
- All characters in '@'-'J' (ASCII 64-74): ambiguous, defaults to Phred+33 (LimitationPolicy)

---

## Edge Cases from Sources

### Documented Edge Cases
1. **Multi-line sequences** (Wikipedia): Legacy Sanger files may split sequences across lines
2. **@ in quality string** (Wikipedia): Makes parsing ambiguous in multi-line files
3. **+ in sequence** (Wikipedia): Unusual but allowed in sequence data
4. **Empty records** (Implementation): Parser should skip blank lines gracefully
5. **Encoding detection failure** (Wikipedia): Default to Phred+33 when ambiguous
6. **Malformed records** (Biopython 1.88 FastqGeneralIterator): missing '@', missing '+', '+' caption ≠ title, whitespace in sequence, seq/qual length mismatch, truncation at EOF → rejected
7. **Out-of-range quality symbol** (Biopython InvalidCharError): rejected, not clamped

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
- `~` (ASCII 126) = Q93 in Phred+33 (max representable)
- `@` (ASCII 64) = Q0 in Phred+64
- `h` (ASCII 104) = Q40 in Phred+64
- `~` (ASCII 126) = Q62 in Phred+64 (max representable)

---

## Implementation Notes

### Current Implementation (FastqParser.cs)
- Supports Phred+33 and Phred+64 encodings
- Phred+33: Q range 0-93 (ASCII 33-126), per Sanger/PacBio standard
- Phred+64: Q range 0-62 (ASCII 64-126), per Illumina 1.3-1.7 standard
- Auto-detection per file via canonical `QualityScoreAnalyzer.DetectEncoding(IEnumerable<string>)`
- Multi-line sequence/quality support via loop
- Paired-end support: interleaving and splitting
- Statistics: total reads, bases, GC content, Q20/Q30 percentages, position quality (1-based, population stddev)
- Filtering: by quality threshold, by length range
- Trimming: quality-based end trimming (Trimmomatic LEADING/TRAILING), adapter removal (cutadapt regular 3' adapter, exact match: leftmost full or 3'-partial occurrence)
- ErrorProbabilityToPhred: p=0 → Q93, every result capped at Q93 (max Sanger representable); p ∉ [0,1] or NaN throws
- Header parsing: first whitespace character separates Id and Description (Biopython `split(None, 1)`); special characters and unicode preserved

### Invariants
1. Sequence length == Quality string length (per record)
2. Quality scores ∈ [0, 93] for Phred+33, [0, 62] for Phred+64
3. Encoding detection is deterministic for given input
4. Round-trip: Parse → Write → Parse yields equivalent data

---

## Review 2026-09 — reference cross-checks (executed)

| Input | Reference output | Old code | Now |
|---|---|---|---|
| `@r1/hhhh`, `@r2/BBBB` (Phred+64 file) | Biopython `fastq-illumina`: r2 = [2,2,2,2] | r2 = [33,33,33,33] (per-record Phred+33) | [2,2,2,2] |
| `@r1/JJJJ`, `@r2/!!!!` | Biopython `fastq`: r1 = [41,41,41,41] | r1 = [10,10,10,10] (per-record Phred+64) | [41,…] |
| quality `J5` | Biopython `fastq`: [41,20] | Phred+64 → [10,0] (order-dependent scan) | [41,20] |
| `@rec/ACGTACGT/+/II` | Biopython ValueError "(8 and 2)" | record with 2 scores (TrimByQuality then IndexOutOfRange) | FormatException |
| `@r1/ACGT/+r2/IIII` | ValueError "captions differ" | accepted | FormatException |
| `@r1\tdesc here` | Biopython id `r1` | id `r1\tdesc` | `r1` |
| `DecodeQualityScores("J5hh", Phred64)` | Biopython InvalidCharError | [10,0,40,40] (clamped) | ArgumentOutOfRangeException |
| `ErrorProbabilityToPhred(1e-10)` | Q93 cap (Phred+33 max) | 100 (while p=0 → 93) | 93 |
| TrimAdapter `AGATCGGAAGAGCACACGTC`, adapter `AGATCGGAAGAGC` | cutadapt 5.2 `-e 0 -O 5`: empty | unchanged | empty |
| TrimAdapter `ACGTACGTAGATCGGAAGAGCTTTTTAGATC` | cutadapt: `ACGTACGT` | `ACGTACGTAGATCGGAAGAGCTTTTT` | `ACGTACGT` |
| WriteToFile first bytes | Biopython/cutadapt reject a BOM | EF BB BF 40 | 40 |

## References

1. Wikipedia contributors. "FASTQ format." Wikipedia, The Free Encyclopedia. https://en.wikipedia.org/wiki/FASTQ_format
2. Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM. (2009). The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic Acids Research. 38(6):1767-1771.
3. NCBI. SRA File Format Guide. https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/
4. Biopython 1.88 `Bio/SeqIO/QualityIO.py` (installed package source, read 2026-09-28).
5. FastQC `PhredEncoding.java`, `PerBaseQualityScores.java` — https://raw.githubusercontent.com/s-andrews/FastQC/master/uk/ac/babraham/FastQC/Sequence/QualityEncoding/PhredEncoding.java
6. cutadapt 5.2 (PyPI), run locally.
