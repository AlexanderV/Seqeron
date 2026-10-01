# complexity_linguistic

Calculate DNA linguistic complexity.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `complexity_linguistic` |
| **Method ID** | `SequenceComplexity.CalculateLinguisticComplexity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates linguistic complexity (LC) as the ratio of observed to possible subwords in a DNA sequence. LC = 1.0 indicates maximum complexity (all possible subwords present), while lower values indicate repetitive or low-complexity sequences. This metric is useful for detecting simple sequence repeats and low-complexity regions.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L121](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L121)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to analyze (min length: 1) |
| `maxWordLength` | integer | No | Maximum word length to consider (default: 10, minimum: 1) |
| `alphabetSize` | integer | No | Fixed alphabet size a (≥ 1 and ≥ distinct symbols; 4 = DNA as in Rosalind LING, 20 = protein). Default: inferred (A/C/G/T ∪ observed symbols) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `complexity` | number | Linguistic complexity (0 to 1) |
| `maxWordLength` | integer | Max word length used |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | Max word length must be at least 1 |
| 1001 | alphabetSize < 1 or smaller than the number of distinct symbols |

## Examples

### Example 1: Complex DNA sequence

**User Prompt:**
> Calculate linguistic complexity of "ATGCGATCGATCG"

**Expected Tool Call:**
```json
{
  "tool": "complexity_linguistic",
  "arguments": {
    "sequence": "ATGCGATCGATCG",
    "maxWordLength": 10
  }
}
```

**Response:**
```json
{
  "complexity": 0.7894736842105263,
  "maxWordLength": 10
}
```

### Example 2: Low complexity (repetitive) sequence

**User Prompt:**
> What's the linguistic complexity of "AAAAAAAAAA"?

**Expected Tool Call:**
```json
{
  "tool": "complexity_linguistic",
  "arguments": {
    "sequence": "AAAAAAAAAA"
  }
}
```

**Response:**
```json
{
  "complexity": 0.20408163265306123,
  "maxWordLength": 10
}
```

### Example 3: Fixed alphabet (Rosalind LING sample)

**Expected Tool Call:**
```json
{ "tool": "complexity_linguistic", "arguments": { "sequence": "ATTTGGATT", "maxWordLength": 100, "alphabetSize": 4 } }
```

**Response:**
```json
{ "complexity": 0.875, "maxWordLength": 100 }
```
35 observed / 40 possible substrings over a = 4 (Rosalind LING sample output 0.875).

## Performance

- **Time Complexity:** O(n × m) where n is sequence length and m is max word length
- **Space Complexity:** O(n × m) for unique substring storage

## See Also

- [linguistic_complexity](linguistic_complexity.md) - SequenceStatistics version
- [complexity_shannon](complexity_shannon.md) - Shannon entropy complexity
- [complexity_dust_score](complexity_dust_score.md) - DUST low-complexity score
