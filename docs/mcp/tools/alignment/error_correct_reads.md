# error_correct_reads

Corrects substitution errors with the Musket k-mer spectrum corrector (Liu, Schmidt & Maskell 2013). k-mers with multiplicity ≥ `minKmerFrequency` are trusted. Default: Musket stage 1 (two-sided conservative correction — per iteration, the first untrusted position with exactly one alternative base making both its leftmost and rightmost covering k-mers trusted is substituted; ambiguous positions unchanged), iterated to convergence over a forward-strand spectrum. `musket: true` runs the full Musket 1.1 corrector (canonical both-strand k-mers, `-maxiter 2`, `-maxerr 4`: adds the aggressive one-sided stage and the voting stage); `countReverseComplement` / `maxIterations` / `maxErrorsPerKmer` override the chosen preset.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `error_correct_reads` |
| **Method ID** | `SequenceAssembler.ErrorCorrectReads` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Corrects substitution errors with the Musket k-mer spectrum corrector (Liu, Schmidt & Maskell 2013). k-mers with multiplicity ≥ `minKmerFrequency` are trusted. Default: Musket stage 1 (two-sided conservative correction — per iteration, the first untrusted position with exactly one alternative base making both its leftmost and rightmost covering k-mers trusted is substituted; ambiguous positions unchanged), iterated to convergence over a forward-strand spectrum. `musket: true` runs the full Musket 1.1 corrector (canonical both-strand k-mers, `-maxiter 2`, `-maxerr 4`: adds the aggressive one-sided stage and the voting stage); `countReverseComplement` / `maxIterations` / `maxErrorsPerKmer` override the chosen preset.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reads` | array<string> | Yes | Input reads. |
| `kmerSize` | integer | No | K-mer size used to build the frequency table (>= 1). (default: 21) |
| `minKmerFrequency` | integer | No | Trusted cut-off (>= 1), Musket `-minmulti`. (default: 3) |
| `musket` | boolean | No | Full Musket 1.1 corrector (both strands, `-maxiter 2`, `-maxerr 4`). (default: false) |
| `countReverseComplement` | boolean | No | Optional: canonical (both-strand) k-mer counting; overrides the preset. |
| `maxIterations` | integer | No | Optional: Musket `-maxiter` (>= 1); overrides the preset. |
| `maxErrorsPerKmer` | integer | No | Optional: Musket `-maxerr` (0 = two-sided stage only); overrides the preset. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `corrected` | array<string> | Corrected reads, upper-cased, in input order (length and count preserved). |

## Errors

| Code | Message |
|------|---------|
| 1004 | k-mer size must be at least 1 |
| 1014 | minKmerFrequency must be at least 1, maxIterations at least 1, maxErrorsPerKmer not negative |
| 1012 | read is null |

## Examples

### Example 1: Fix one error

**Tool Call:**
```json
{
  "tool": "error_correct_reads",
  "arguments": {
    "reads": [
      "ACGTACGT",
      "ACGTACGT",
      "ACGTACGT",
      "ACGTTCGT"
    ],
    "kmerSize": 4,
    "minKmerFrequency": 2
  }
}
```

**Response:**
```json
{
  "corrected": [
    "ACGTACGT",
    "ACGTACGT",
    "ACGTACGT",
    "ACGTACGT"
  ]
}
```

### Example 2: Empty input

**Tool Call:**
```json
{
  "tool": "error_correct_reads",
  "arguments": {
    "reads": [],
    "kmerSize": 4,
    "minKmerFrequency": 2
  }
}
```

**Response:**
```json
{
  "corrected": []
}
```

### Example 3: Reverse-strand reads need canonical k-mers (Musket)

**Tool Call:**
```json
{
  "tool": "error_correct_reads",
  "arguments": {
    "reads": ["GAAGCCTGTAATC", "GAAGCCTGTAATC", "GAAGCCTGTAATC", "GATTACGGGCTTC"],
    "kmerSize": 3,
    "minKmerFrequency": 2,
    "musket": true
  }
}
```

**Response:**
```json
{
  "corrected": ["GAAGCCTGTAATC", "GAAGCCTGTAATC", "GAAGCCTGTAATC", "GATTACAGGCTTC"]
}
```

The first three reads are the reverse complement of `GATTACAGGCTTC`. With canonical counting every genomic
3-mer is trusted and the C→G error at index 6 is fixed. Without it (`musket: false`, default) the forward-only
spectrum never sees the read's genomic k-mers, and the untrusted `GAT` at index 0 is changed to `AAT`
(`AATTACGGGCTTC`), so use `musket: true` or `countReverseComplement: true` for double-stranded data.

## Worked Example

Three copies of `ACGTACGT` make its 4-mers trusted (freq >= 2). The fourth read `ACGTTCGT` has one error at index 4; A is the only base making both its leftmost (`CGTT`→`CGTA`) and rightmost (`TCGT`→`ACGT`) covering 4-mers trusted, so it is corrected back to `ACGTACGT`.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAssembler.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
