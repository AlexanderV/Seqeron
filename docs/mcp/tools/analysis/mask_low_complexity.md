# mask_low_complexity

Mask low-complexity regions of a DNA sequence with symmetric DUST (SDUST).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mask_low_complexity` |
| **Method ID** | `SequenceComplexity.MaskLowComplexity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Masks low-complexity regions with the **symmetric DUST** algorithm (Morgulis et al. 2006;
identical output to lh3/sdust). Every *perfect interval* — a subsequence of
at most `windowSize` bases whose DUST score `Σ c(c−1)/2 / (ℓ−1)` exceeds `threshold` and is not
exceeded by any of its sub-intervals — is masked; overlapping/adjacent intervals
are merged. N and other IUPAC codes are accepted (case-insensitive): each maximal A/C/G/T run is scanned as an
independent sequence, as sdust specifies, and a non-ACGT base is never inside a perfect interval. `linker`
reproduces NCBI dustmasker's `-linker` (consecutive intervals separated by fewer than `linker` unmasked bases are
merged; default 1 = sdust/dustmasker default) and `softMask` its `-outfmt fasta` output (masked bases lower case,
the rest upper case; `maskChar` ignored). Otherwise masked bases become `maskChar` and the output is upper case.
The result is the same length as the input. `windowSize` must be ≥ 3, `threshold`
finite and ≥ 0 (default 2.0 = dustmasker level 20), `linker` 1–32.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L539](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L539)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence: A/C/G/T + IUPAC codes (N, R, Y, …), case-insensitive (min length 1) |
| `windowSize` | integer | No | Window size (default 64) |
| `threshold` | number | No | DUST threshold above which to mask (default 2.0) |
| `maskChar` | string | No | Mask character (default `N`; ignored when `softMask`) |
| `linker` | integer | No | dustmasker linker, 1–32 (default 1 = sdust) |
| `softMask` | boolean | No | Lower-case masking like dustmasker `-outfmt fasta` (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `masked` | string | Masked sequence (same length as input) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1003 | windowSize < 3, threshold negative/NaN/infinite, or linker outside 1–32 |

## Examples

### Example 1: Fully masked poly-A

**User Prompt:**
> Mask low-complexity regions of a 100-nt poly-A tract with 'X' (window 64, threshold 1.0).

**Expected Tool Call:**
```json
{
  "tool": "mask_low_complexity",
  "arguments": { "sequence": "AAAA…(100×)", "windowSize": 64, "threshold": 1.0, "maskChar": "X" }
}
```

**Response:**
```json
{ "masked": "XXXX…(100×)" }
```
lh3/sdust (`-w 64 -t 10`) reports the interval [0,100), so every position is masked.

### Example 2: High complexity preserved

**User Prompt:**
> Mask a varied 78-bp sequence at threshold 10.0.

**Expected Tool Call:**
```json
{
  "tool": "mask_low_complexity",
  "arguments": { "sequence": "ATGCTAGCATGCA…(78 bp)", "windowSize": 64, "threshold": 10.0 }
}
```

**Response:**
```json
{ "masked": "ATGCTAGCATGCA…(78 bp, unchanged)" }
```

### Example 3: Input with N, soft mask

**Input:** `{ "sequence": "ACGTNNAAAAAAAAAAAANACGTACACACACACACACANNGGGCCCTAGGTCA", "softMask": true }`
→ sdust on each ACGT run masks A×12 [6,18) and the (AC) run [23,38) →
**Response:** `{ "masked": "ACGTNNaaaaaaaaaaaaNACGTacacacacacacacaNNGGGCCCTAGGTCA" }`

### Example 4: dustmasker linker

**Input:** `{ "sequence": "ACGTGCATGCAAAAAAAAAAAAAAAAGCTAGCATCGACTGCAGCACACACACACACACAGATCGATCGTACGGTGCATGACAAAAAAAAAAAAACT", "linker": 18, "softMask": true }`
→ dustmasker 2.12.0 `-linker 18`: intervals 10–58 and 81–93 (17 unmasked bases between the first two runs are bridged) →
**Response:** `{ "masked": "ACGTGCATGCaaaaaaaaaaaaaaaagctagcatcgactgcagcacacacacacacacaGATCGATCGTACGGTGCATGACaaaaaaaaaaaaaCT" }`

## Performance

- **Time Complexity:** O(n · windowSize) typical; O(n · windowSize³) worst case (SDUST perfect-interval search).
- **Space Complexity:** O(n).

## See Also

- [find_low_complexity_regions](find_low_complexity_regions.md)
- [dust_score](dust_score.md)
