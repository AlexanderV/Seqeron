# complexity_kmer_entropy

Calculate k-mer based entropy for DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `complexity_kmer_entropy` |
| **Method ID** | `SequenceComplexity.CalculateKmerEntropy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates Shannon entropy using k-mer frequencies in a DNA sequence. This provides a more nuanced complexity measure than single-base entropy by considering local sequence context. Higher entropy values indicate more diverse k-mer composition. Optional `correction` applies a finite-sample bias correction — Miller–Madow (1955: + (D − 1)/(2N) nats, D = observed k-mers; = R `entropy::entropy.MillerMadow`) or Grassberger (2003, arXiv:physics/0307138 eq. 35: ln N − (1/N) Σ n_i G(n_i)) — and `normalize` divides by log₂ N (BBTools `EntropyTracker` 0–1 scale; 0 when N ≤ 1).

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L283](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L283) — `CalculateKmerEntropy`
- Algorithm doc: [K-mer_Entropy.md](../../../algorithms/Complexity/K-mer_Entropy.md) (SEQ-COMPLEX-KMER-001)

H = −Σ p_i·log₂ p_i over the N = L − k + 1 overlapping k-mers, p_i = n_i/N (bits). Example 1: ATGCATGCAT, k=2 → AT=3, TG=2, GC=2, CA=2 of N=9 → 1.9749375 bits (scipy.stats.entropy cross-check).

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to analyze (must be valid DNA) |
| `k` | integer | No | K-mer size (default: 2 for dinucleotides, minimum: 1) |
| `correction` | string | No | `none` (plug-in, default), `millerMadow`, `grassberger` (case-insensitive) |
| `normalize` | boolean | No | Divide by log₂ N, N = L − k + 1 (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entropy` | number | K-mer entropy (bits) |
| `k` | integer | K-mer size used |
| `correction` | string | Correction applied (`none` / `millerMadow` / `grassberger`) |
| `normalized` | boolean | Whether the value was divided by log₂ N |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | K must be at least 1 |
| 2001 | Invalid DNA sequence |
| 1002 | Unknown k-mer entropy correction |

## Examples

### Example 1: Dinucleotide entropy

**User Prompt:**
> Calculate 2-mer entropy of "ATGCATGCAT"

**Expected Tool Call:**
```json
{
  "tool": "complexity_kmer_entropy",
  "arguments": {
    "sequence": "ATGCATGCAT",
    "k": 2
  }
}
```

**Response:**
```json
{
  "entropy": 1.974937501201927,
  "k": 2,
  "correction": "none",
  "normalized": false
}
```

### Example 2: Low complexity sequence

**User Prompt:**
> What's the k-mer entropy of "AAAAAAAAAA"?

**Expected Tool Call:**
```json
{
  "tool": "complexity_kmer_entropy",
  "arguments": {
    "sequence": "AAAAAAAAAA"
  }
}
```

**Response:**
```json
{
  "entropy": 0,
  "k": 2,
  "correction": "none",
  "normalized": false
}
```

### Example 3: Bias-corrected, normalised

**User Prompt:**
> Miller–Madow corrected 2-mer entropy of "ATATAT", and the normalised plug-in value

**Expected Tool Calls:**
```json
{ "tool": "complexity_kmer_entropy", "arguments": { "sequence": "ATATAT", "k": 2, "correction": "millerMadow" } }
{ "tool": "complexity_kmer_entropy", "arguments": { "sequence": "ATATAT", "k": 2, "normalize": true } }
```

**Responses:** `{"entropy": 1.115220098543565, "k": 2, "correction": "millerMadow", "normalized": false}` (R `entropy.MillerMadow(c(3,2), unit="log2")`); `{"entropy": 0.4181656600790516, "k": 2, "correction": "none", "normalized": true}` (BBTools `EntropyTracker.calcEntropy` 0.41816565). Grassberger: 1.2692841903863027.

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(min(4^k, L − k + 1)) distinct k-mers stored

## See Also

- [kmer_entropy](kmer_entropy.md) - KmerAnalyzer version
- [complexity_shannon](complexity_shannon.md) - Single-base entropy
- [complexity_linguistic](complexity_linguistic.md) - Linguistic complexity
