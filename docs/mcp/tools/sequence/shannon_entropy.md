# shannon_entropy

Calculate Shannon entropy of a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `shannon_entropy` |
| **Method ID** | `SequenceStatistics.CalculateShannonEntropy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the Shannon entropy `H = −Σ pᵢ log₂ pᵢ` (Shannon 1948; bits per symbol) of the per-letter composition — the same value as `scipy.stats.entropy(counts, base=2)`.

**Alphabet:** every letter (case-folded) is its own symbol, so DNA, RNA and protein input all work; `N` and IUPAC ambiguity codes are counted as distinct symbols, `T` and `U` are **not** merged, and non-letters (gaps `-`, digits, `*`, whitespace) are excluded from numerator and denominator. The maximum is `log₂ m` for m distinct letters (2 bits for A/C/G/T in equal proportion); a homopolymer (e.g. `AAAAAA`) scores 0. For the nucleotide-only {A, C, G, T/U} entropy use [`complexity_shannon`](complexity_shannon.md).

Example: `ATGCATGCAT` has A = 3, T = 3, G = 2, C = 2 → H = 1.9709505944546686 bits (not 2: the four bases are not equiprobable).

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L1884](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L1884)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to analyze (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entropy` | number | Shannon entropy value (bits per symbol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: DNA sequence entropy

**User Prompt:**
> What's the Shannon entropy of "ATGCATGCAT"?

**Expected Tool Call:**
```json
{
  "tool": "shannon_entropy",
  "arguments": {
    "sequence": "ATGCATGCAT"
  }
}
```

**Response:**
```json
{
  "entropy": 1.9709505944546686
}
```

### Example 2: Low complexity sequence

**User Prompt:**
> Calculate entropy of "AAAAAAAAAA"

**Expected Tool Call:**
```json
{
  "tool": "shannon_entropy",
  "arguments": {
    "sequence": "AAAAAAAAAA"
  }
}
```

**Response:**
```json
{
  "entropy": 0.0
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(k) where k is alphabet size

## See Also

- [linguistic_complexity](linguistic_complexity.md) - Linguistic complexity measure
- [nucleotide_composition](nucleotide_composition.md) - Base composition analysis
