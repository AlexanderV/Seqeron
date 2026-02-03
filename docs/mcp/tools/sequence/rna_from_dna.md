# rna_from_dna

Transcribe DNA to RNA.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `rna_from_dna` |
| **Method ID** | `RnaSequence.FromDna` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Transcribes a DNA sequence to RNA by replacing thymine (T) with uracil (U). This simulates the biological transcription process where DNA is used as a template to produce messenger RNA.

## Core Documentation Reference

- Source: [RnaSequence.cs#L147](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/RnaSequence.cs#L147)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to transcribe (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `rna` | string | The transcribed RNA sequence |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 2001 | Invalid DNA sequence |

## Examples

### Example 1: Simple transcription

**User Prompt:**
> Transcribe "ATGC" to RNA

**Expected Tool Call:**
```json
{
  "tool": "rna_from_dna",
  "arguments": {
    "sequence": "ATGC"
  }
}
```

**Response:**
```json
{
  "rna": "AUGC"
}
```

### Example 2: All thymine

**User Prompt:**
> Convert "TTTT" to RNA

**Expected Tool Call:**
```json
{
  "tool": "rna_from_dna",
  "arguments": {
    "sequence": "TTTT"
  }
}
```

**Response:**
```json
{
  "rna": "UUUU"
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(n)

## See Also

- [dna_validate](dna_validate.md) - Validate DNA sequences
- [rna_validate](rna_validate.md) - Validate RNA sequences
