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

Masks low-complexity regions in a DNA sequence using the symmetric DUST (SDUST) algorithm of Morgulis et al. (2006) — output identical to lh3/sdust: every perfect interval of at most `windowSize` bases whose DUST score `Σ c(c−1)/2 / (ℓ−1)` exceeds `threshold` is masked. This is commonly used as a preprocessing step before BLAST searches or other sequence analyses to prevent spurious matches caused by simple/repetitive sequences. Low-complexity regions are replaced with a mask character (typically 'N'). N and other IUPAC codes are accepted (case-insensitive); each maximal A/C/G/T run is scanned as an independent sequence, as sdust specifies ("N effectively breaks input into pieces of independent sequences"; sdust's code itself carries its scoring window across N, which this tool does not reproduce). `linker` reproduces NCBI dustmasker's `-linker` (intervals separated by fewer than `linker` unmasked bases are merged; default 1 = sdust) and `softMask` its `-outfmt fasta` lower-case masking.

`engine: "dustmasker"` switches to a line-by-line port of NCBI `CSymDustMasker` (symdust.cpp) driven by dustmasker's `GetDustMasks_SkipNs` (dust_mask_app.cpp), reproducing dustmasker 2.12.0 output exactly (3 000+ randomized inputs with N/IUPAC, 0 mismatches; B04 F53): IUPAC codes are scanned as bases (C/G/T → 1/2/3, N → the toolkit's deterministic `CRandom` 2-bit code, every other code → A), only N runs longer than the window plus leading/trailing N runs cut the scan — and those N runs are themselves reported as masked —, and a window holding a single triplet value (homopolymer ≥ window) is masked without a score test.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L838](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L838)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to mask: A/C/G/T + IUPAC codes (N, R, Y, …), case-insensitive |
| `windowSize` | integer | No | Window size for analysis (default: 64; values 1–2 are rejected by the core, which requires ≥ 3) |
| `threshold` | number | No | DUST threshold above which to mask (default: 2.0) |
| `maskChar` | string | No | Character to use for masking (default: 'N'; ignored when `softMask`) |
| `linker` | integer | No | dustmasker linker, 1–32 (default: 1 = sdust) |
| `softMask` | boolean | No | Lower-case masking like dustmasker `-outfmt fasta` (default: false) |
| `engine` | string | No | `sdust` (default) or `dustmasker` — exact NCBI dustmasker 2.12.0 parity (IUPAC scanned as bases; N runs > window and leading/trailing N runs cut the scan and are masked; window 8–64, 10·threshold an integer 2–64) |

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
| 1003 | engine `dustmasker`: windowSize outside 8–64, 10·threshold not an integer 2–64, non-IUPAC input; unknown engine name |

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

### Example 4: dustmasker engine

**Input:** `{ "sequence": "NNACGTTGCAAAAAAAAAAAACGTRRRRRRRRRRRRTGCANNNNNNNNNNACGTTGCAANNNN", "windowSize": 8, "threshold": 2.0, "softMask": true, "engine": "dustmasker" }`
→ dustmasker 2.12.0 `-window 8 -level 20 -outfmt fasta` →
**Response:** `{ "maskedSequence": "nnACGTTGCaaaaaaaaaaaaCGTrrrrrrrrrrrrTGCAnnnnnnnnnnACGTTGCAAnnnn", "originalLength": 63, "maskChar": "N" }`

## Performance

- **Time Complexity:** O(n × w) typical; O(n × w³) worst case (SDUST perfect-interval search), w = window size
- **Space Complexity:** O(n) for the masked sequence

## See Also

- [complexity_dust_score](complexity_dust_score.md) - Calculate DUST score
- [shannon_entropy](shannon_entropy.md) - Information-theoretic complexity
