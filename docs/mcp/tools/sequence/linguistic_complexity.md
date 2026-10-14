# linguistic_complexity

Calculate linguistic complexity of a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `linguistic_complexity` |
| **Method ID** | `SequenceStatistics.CalculateLinguisticComplexity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates linguistic complexity in the Orlov & Potapov (2004) summation form

`LC = Σ_{k=1..m} V_k / Σ_{k=1..m} V_max,k`, `V_max,k = min(a^k, N − k + 1)`, `m = min(maxK, N)`,

where `V_k` is the number of distinct k-words of the (upper-cased) sequence and `a` is the alphabet size — `alphabetSize` when given, otherwise `|{A, C, G, T (U for RNA)} ∪ observed symbols|` (canonical `SequenceComplexity`, B04 F20), so `LC ≤ 1` for any input (`ACGTN` → 1.0). With `maxK ≥ N` it equals Troyanskaya et al. (2002) / Rosalind LING. 1 means every possible word is present; repeats and low-complexity stretches lower it.

Relation to [`complexity_linguistic`](complexity_linguistic.md): the same canonical measure (both call `SequenceComplexity.CalculateLinguisticComplexity`); only the word-length parameter differs — `maxK` here (default 6) vs `maxWordLength` there (default 10).

Examples (exact fractions, Python reference): `ATGCATGCATGC` → V = 4+4+4+4+4+4 = 24 of V_max = 4+11+10+9+8+7 = 49 → 24/49; `AAAAAAAAAA` → 6 of 4+9+8+7+6+5 = 39 → 2/13.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L1993](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L1993)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to analyze (min length: 1) |
| `maxK` | integer | No | Maximum k-mer length to consider (default: 6) |
| `alphabetSize` | integer | No | Fixed alphabet size a (≥ 1 and ≥ distinct symbols; 4 = DNA as in Rosalind LING, 20 = protein). Default: inferred (A/C/G/T ∪ observed symbols) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `complexity` | number | Linguistic complexity (0-1) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | alphabetSize < 1 or smaller than the number of distinct symbols |

## Examples

### Example 1: Periodic sequence (ATGC×3)

**User Prompt:**
> What's the linguistic complexity of "ATGCATGCATGC"?

**Expected Tool Call:**
```json
{
  "tool": "linguistic_complexity",
  "arguments": {
    "sequence": "ATGCATGCATGC"
  }
}
```

**Response:**
```json
{
  "complexity": 0.4897959183673469
}
```

### Example 2: Low complexity sequence

**User Prompt:**
> Calculate complexity of "AAAAAAAAAA"

**Expected Tool Call:**
```json
{
  "tool": "linguistic_complexity",
  "arguments": {
    "sequence": "AAAAAAAAAA"
  }
}
```

**Response:**
```json
{
  "complexity": 0.15384615384615385
}
```

### Example 3: Protein with a fixed 20-letter alphabet

**Expected Tool Call:**
```json
{ "tool": "linguistic_complexity", "arguments": { "sequence": "MKVLAAGIVGLLLAA", "maxK": 15, "alphabetSize": 20 } }
```

**Response:**
```json
{ "complexity": 0.9 }
```

## Performance

- **Time Complexity:** O(n × maxK) where n is sequence length
- **Space Complexity:** O(n × maxK) for the k-mer sets (at most N − k + 1 distinct words per k)

## See Also

- [shannon_entropy](shannon_entropy.md) - Shannon entropy measure
- [nucleotide_composition](nucleotide_composition.md) - Base composition analysis
