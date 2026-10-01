# complexity_mask_low

Mask low-complexity regions using DUST algorithm.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `complexity_mask_low` |
| **Method ID** | `SequenceComplexity.MaskLowComplexity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Masks low-complexity regions in a DNA sequence using the symmetric DUST (SDUST) algorithm of Morgulis et al. (2006) — output identical to lh3/sdust: every perfect interval of at most `windowSize` bases whose DUST score `Σ c(c−1)/2 / (ℓ−1)` exceeds `threshold` is masked. This is commonly used as a preprocessing step before BLAST searches or other sequence analyses to prevent spurious matches caused by simple/repetitive sequences. Low-complexity regions are replaced with a mask character (typically 'N'). N and other IUPAC codes are accepted (case-insensitive); each maximal A/C/G/T run is scanned as an independent sequence, as sdust specifies. `linker` reproduces NCBI dustmasker's `-linker` (intervals separated by fewer than `linker` unmasked bases are merged; default 1 = sdust) and `softMask` its `-outfmt fasta` lower-case masking.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L826](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L826)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to mask: A/C/G/T + IUPAC codes (N, R, Y, …), case-insensitive |
| `windowSize` | integer | No | Window size for analysis (default: 64; values 1–2 are rejected by the core, which requires ≥ 3) |
| `threshold` | number | No | DUST threshold above which to mask (default: 2.0) |
| `maskChar` | string | No | Character to use for masking (default: 'N'; ignored when `softMask`) |
| `linker` | integer | No | dustmasker linker, 1–32 (default: 1 = sdust) |
| `softMask` | boolean | No | Lower-case masking like dustmasker `-outfmt fasta` (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `maskedSequence` | string | Sequence with low-complexity regions masked |
| `originalLength` | integer | Original sequence length |
| `maskChar` | string | Character used for masking |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | Window size must be at least 1 |
| 2001 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1003 | windowSize < 3, threshold negative/NaN/infinite, or linker outside 1–32 (core validation) |

## Examples

### Example 1: Mask repetitive region

**User Prompt:**
> Mask the low-complexity regions in this sequence

**Expected Tool Call:**
```json
{
  "tool": "complexity_mask_low",
  "arguments": {
    "sequence": "ATGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGC",
    "windowSize": 64,
    "threshold": 2.0,
    "maskChar": "N"
  }
}
```

**Response:**
```json
{
  "maskedSequence": "ATGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATGC",
  "originalLength": 68,
  "maskChar": "N"
}
```

### Example 2: Custom threshold

**User Prompt:**
> Mask low-complexity regions with a stricter threshold of 1.0

**Expected Tool Call:**
```json
{
  "tool": "complexity_mask_low",
  "arguments": {
    "sequence": "ATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGC",
    "threshold": 1.0
  }
}
```

**Response:** lh3/sdust (`-w 64 -t 10`) reports [0, 69), so all 69 bases are masked:
```json
{
  "maskedSequence": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
  "originalLength": 69,
  "maskChar": "N"
}
```

### Example 3: N-containing input with 'X'

**Input:** `{ "sequence": "ACGTNNAAAAAAAAAAAANACGTACACACACACACACANNGGGCCCTAGGTCA", "maskChar": "X" }`
→ sdust on each ACGT run: [6,18) and [23,38) →
**Response:** `{ "maskedSequence": "ACGTNNXXXXXXXXXXXXNACGTXXXXXXXXXXXXXXXNNGGGCCCTAGGTCA", "originalLength": 53, "maskChar": "X" }`

## Performance

- **Time Complexity:** O(n × w) typical; O(n × w³) worst case (SDUST perfect-interval search), w = window size
- **Space Complexity:** O(n) for the masked sequence

## See Also

- [complexity_dust_score](complexity_dust_score.md) - Calculate DUST score
- [shannon_entropy](shannon_entropy.md) - Information-theoretic complexity
