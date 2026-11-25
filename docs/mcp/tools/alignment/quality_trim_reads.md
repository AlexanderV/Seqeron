# quality_trim_reads

Trim low-quality read ends with the BWA/cutadapt running-sum algorithm (cutadapt quality_trim_index), with optional separate 5'/3' cutoffs and Phred+33/Phred+64 decoding, dropping reads shorter than minLength after trimming.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `quality_trim_reads` |
| **Method ID** | `SequenceAssembler.QualityTrimReads` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Trim low-quality read ends with the BWA/cutadapt running-sum algorithm (cutadapt quality_trim_index), with optional separate 5'/3' cutoffs and Phred+33/Phred+64 decoding, dropping reads shorter than minLength after trimming.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1634](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1634)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reads` | array<object> | Yes | Reads, each an object with sequence and a quality string of equal length (Phred+phredOffset ASCII). |
| `minQuality` | integer | No | Quality cutoff subtracted from each Phred score at both ends (values < 1 disable trimming). (default: 20) |
| `minLength` | integer | No | Minimum length after trimming for a read to be kept. (default: 50) |
| `cutoffFront` | integer | No | Optional 5' cutoff overriding minQuality for the 5' end (< 1 leaves the 5' end untrimmed; cutadapt `-q 5',3'`). |
| `cutoffBack` | integer | No | Optional 3' cutoff overriding minQuality for the 3' end (< 1 leaves the 3' end untrimmed). A single cutadapt `-q C` / BWA `aln -q C` is 3'-only: `cutoffFront: 0, cutoffBack: C`. |
| `phredOffset` | integer | No | ASCII offset of the quality encoding: 33 (Sanger/Illumina 1.8+) or 64 (Illumina 1.3-1.7), as cutadapt `--quality-base`. (default: 33) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `trimmed` | array<string> | Surviving trimmed sequences, in input order. |

## Errors

| Code | Message |
|------|---------|
| 1011 | sequence and quality must be the same length |
| 1012 | read is null |
| 1013 | quality character outside [phredOffset, '~'] or phredOffset outside [33, 126] |

## Examples

### Example 1: High quality kept

**Tool Call:**
```json
{
  "tool": "quality_trim_reads",
  "arguments": {
    "reads": [
      {
        "sequence": "ACGTACGTAC",
        "quality": "IIIIIIIIII"
      }
    ],
    "minQuality": 20,
    "minLength": 5
  }
}
```

**Response:**
```json
{
  "trimmed": [
    "ACGTACGTAC"
  ]
}
```

### Example 2: Low quality dropped

**Tool Call:**
```json
{
  "tool": "quality_trim_reads",
  "arguments": {
    "reads": [
      {
        "sequence": "ACGTACGTAC",
        "quality": "!!!!!!!!!!"
      }
    ],
    "minQuality": 20,
    "minLength": 5
  }
}
```

**Response:**
```json
{
  "trimmed": []
}
```

### Example 3: cutadapt worked example in Phred+64

**Tool Call:**
```json
{
  "tool": "quality_trim_reads",
  "arguments": {
    "reads": [
      {
        "sequence": "ACGTACGTAC",
        "quality": "jhZ[GFJCAB"
      }
    ],
    "minQuality": 10,
    "minLength": 1,
    "phredOffset": 64
  }
}
```

**Response:**
```json
{
  "trimmed": [
    "ACGT"
  ]
}
```

## Worked Example

Phred+33 `I` = 40 (>= cutoff 20) so the read is kept whole; `!` = 0 (< cutoff) so the read is fully trimmed and dropped for being under minLength.

Cutadapt documentation example: qualities 42,40,26,27,8,7,11,4,2,3 (Phred+64 `jhZ[GFJCAB`), cutoff 10. 3' pass accumulates cutoff − q from the end: 7, 15, 21, 20, 23, 25, then 25 − 17 = 8 and 8 − 16 < 0 stops the scan; the last strict maximum (25) is at index 4, so the read keeps bases [0, 4) = `ACGT` (cutadapt 5.2 `quality_trim_index` returns (0, 4)).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1634](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs#L1634)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
