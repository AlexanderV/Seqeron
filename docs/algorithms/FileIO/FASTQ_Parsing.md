# FASTQ Parsing

## Overview

FASTQ is a text-based format for storing nucleotide sequences with per-base quality scores. Originally developed at the Wellcome Trust Sanger Institute (~2000), it has become the de facto standard output format for high-throughput sequencing instruments.

---

## Format Specification

### Structure
Each record consists of exactly 4 lines:
1. **Header line:** Begins with `@`, followed by sequence identifier and optional description
2. **Sequence line:** Raw nucleotide sequence (ACGT/N)
3. **Separator line:** Begins with `+`, optionally followed by identifier
4. **Quality line:** ASCII-encoded quality scores (same length as sequence)

### Example
```
@SEQ_ID description text
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

---

## Quality Encoding

### Phred Quality Score
The Phred quality score Q represents the probability p of an incorrect base call:

$$Q = -10 \cdot \log_{10}(p)$$

Conversely:

$$p = 10^{-Q/10}$$

### Quality Score Reference

| Q Score | Error Probability | Accuracy |
|---------|-------------------|----------|
| 10 | 10% (1 in 10) | 90% |
| 20 | 1% (1 in 100) | 99% |
| 30 | 0.1% (1 in 1,000) | 99.9% |
| 40 | 0.01% (1 in 10,000) | 99.99% |

### ASCII Encoding Schemes

| Format | Offset | ASCII Range | Typical Q Range | Usage |
|--------|--------|-------------|-----------------|-------|
| **Phred+33** (Sanger) | 33 | '!' (33) to '~' (126) | 0-41 | Illumina 1.8+, PacBio, Nanopore |
| **Phred+64** (Illumina) | 64 | '@' (64) to 'h' (104) | 0-40 | Illumina 1.3-1.7 (legacy) |

### Encoding Detection Heuristic
```
For each character c in quality string:
  if c < '@' (ASCII 64): return Phred+33
  if c > 'I' (ASCII 73): return Phred+64
Default: return Phred+33 (most common modern format)
```

---

## Paired-End FASTQ

### Separate Files Format
Forward reads in `sample_R1.fastq`, reverse reads in `sample_R2.fastq`:
- Same read order in both files
- Matching identifiers (possibly with /1 and /2 suffixes)

### Interleaved Format
Alternating forward and reverse reads in single file:
```
@read1/1
<sequence1_forward>
+
<quality1_forward>
@read1/2
<sequence1_reverse>
+
<quality1_reverse>
@read2/1
...
```

---

## Implementation

### Class: `FastqParser`
**Namespace:** `Seqeron.Genomics.IO`

### Data Structures

```csharp
public readonly record struct FastqRecord(
    string Id,
    string Description,
    string Sequence,
    string QualityString,
    IReadOnlyList<int> QualityScores);

public enum QualityEncoding { Phred33, Phred64, Auto }

public readonly record struct FastqStatistics(
    int TotalReads,
    long TotalBases,
    double MeanReadLength,
    double MeanQuality,
    int MinReadLength,
    int MaxReadLength,
    double Q20Percentage,
    double Q30Percentage,
    double GcContent);
```

### Core Methods

| Method | Purpose | Complexity |
|--------|---------|------------|
| `Parse(content, encoding)` | Parse FASTQ from string | O(n) |
| `ParseFile(filePath, encoding)` | Parse FASTQ from file | O(n) |
| `DetectEncoding(qualityString)` | Auto-detect Phred offset | O(m) |
| `DecodeQualityScores(qualityString, encoding)` | Convert ASCII to Phred scores | O(m) |
| `EncodeQualityScores(scores, encoding)` | Convert Phred scores to ASCII | O(m) |
| `PhredToErrorProbability(phred)` | Q → p conversion | O(1) |
| `ErrorProbabilityToPhred(probability)` | p → Q conversion | O(1) |
| `FilterByQuality(records, minQ)` | Filter by average quality | O(n×m) |
| `FilterByLength(records, min, max)` | Filter by sequence length | O(n) |
| `TrimByQuality(record, minQ)` | Trim low-quality ends | O(m) |
| `TrimAdapter(record, adapter, minOverlap)` | Remove adapter sequences | O(m×a) |
| `CalculateStatistics(records)` | Compute summary statistics | O(n×m) |
| `InterleavePairedReads(r1, r2)` | Combine paired files | O(n) |
| `SplitInterleavedReads(interleaved)` | Separate interleaved file | O(n) |

Where: n = number of records, m = average sequence length, a = adapter length

---

## Edge Cases

### Parsing Edge Cases
1. **Empty content:** Returns empty enumerable
2. **Null content:** Returns empty enumerable (no exception)
3. **Multi-line sequences:** Supported (reads until '+' line)
4. **Multi-line quality:** Supported (reads until sequence length matched)
5. **'+' in sequence:** Potential ambiguity (handled by reading until standalone '+')
6. **Missing sequence after header:** Skipped
7. **Sequence/quality length mismatch:** Quality truncated/padded to match

### Encoding Edge Cases
1. **Empty quality string:** Defaults to Phred+33
2. **Ambiguous range ('@'-'I'):** Defaults to Phred+33
3. **Mixed encodings in file:** Detection uses first indicative character

### Filter/Trim Edge Cases
1. **All bases below threshold:** Returns empty sequence
2. **No adapter match:** Returns unchanged record
3. **Partial adapter at end:** Trimmed if overlap ≥ minOverlap

---

## Invariants

1. `|Sequence| == |QualityString|` for each valid record
2. `QualityScores[i] ∈ [0, 93]` for Phred+33
3. `QualityScores[i] ∈ [0, 62]` for Phred+64
4. `DetectEncoding(s) ∈ {Phred33, Phred64}` (deterministic)
5. `DecodeQualityScores(EncodeQualityScores(Q, E), E) == Q` (round-trip)

---

## Sources

1. Wikipedia. "FASTQ format." https://en.wikipedia.org/wiki/FASTQ_format
2. Cock PJA et al. (2009). Nucleic Acids Research. 38(6):1767-1771. DOI: 10.1093/nar/gkp1137
3. NCBI SRA File Format Guide. https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/
