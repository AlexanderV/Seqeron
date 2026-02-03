# dna_reverse_complement

Get the reverse complement of a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `dna_reverse_complement` |
| **Method ID** | `DnaSequence.GetReverseComplementString` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the reverse complement of a DNA sequence. First, each nucleotide is complemented (A↔T, C↔G), then the sequence is reversed. This is essential for working with double-stranded DNA, primer design, and understanding gene orientation.

## Core Documentation Reference

- Source: [DnaSequence.cs#L149](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/DnaSequence.cs#L149)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `reverseComplement` | string | The reverse complement sequence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 2001 | Invalid DNA sequence |

## Examples

### Example 1: Simple reverse complement

**User Prompt:**
> What is the reverse complement of ATGC?

**Expected Tool Call:**
```json
{
  "tool": "dna_reverse_complement",
  "arguments": {
    "sequence": "ATGC"
  }
}
```

**Response:**
```json
{
  "reverseComplement": "GCAT"
}
```

### Example 2: Self-complementary sequence

**User Prompt:**
> Reverse complement of AATT

**Expected Tool Call:**
```json
{
  "tool": "dna_reverse_complement",
  "arguments": {
    "sequence": "AATT"
  }
}
```

**Response:**
```json
{
  "reverseComplement": "AATT"
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(n)

## See Also

- [dna_validate](dna_validate.md) - Validate DNA sequences
- [rna_from_dna](rna_from_dna.md) - Transcribe DNA to RNA
