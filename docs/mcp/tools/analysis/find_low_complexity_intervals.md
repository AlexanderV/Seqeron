# find_low_complexity_intervals

SDUST low-complexity intervals (lh3/sdust, dustmasker -outfmt interval).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_low_complexity_intervals` |
| **Method ID** | `SequenceComplexity.FindLowComplexityIntervals` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the **SDUST** low-complexity intervals (Morgulis et al. 2006; symmetric DUST) as 0-based half-open
`[start, end)` pairs in ascending order — the native output of lh3/sdust (`name start end`) and, with `end − 1`,
of NCBI dustmasker `-outfmt interval`. These are exactly the positions `mask_low_complexity` masks. N and other
IUPAC codes split the input into independently scanned A/C/G/T runs (sdust's contract). `linker` reproduces
dustmasker `-linker` (intervals separated by fewer than `linker` unmasked bases are merged; default 1 = sdust).

`engine: "dustmasker"` switches to a line-by-line port of NCBI `CSymDustMasker` (symdust.cpp) driven by dustmasker's `GetDustMasks_SkipNs` (dust_mask_app.cpp), reproducing dustmasker 2.12.0 output exactly (3 000+ randomized inputs with N/IUPAC, 0 mismatches; B04 F53): IUPAC codes are scanned as bases (C/G/T → 1/2/3, N → the toolkit's deterministic `CRandom` 2-bit code, every other code → A), only N runs longer than the window plus leading/trailing N runs cut the scan — and those N runs are themselves reported as masked —, and a window holding a single triplet value (homopolymer ≥ window) is masked without a score test.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L872](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L872)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (A/C/G/T plus IUPAC codes such as N; case-insensitive) |
| `windowSize` | integer | No | SDUST window length in bases (default 64) |
| `threshold` | number | No | DUST score threshold; intervals scoring strictly above it are reported (2.0 = dustmasker level 20) (default 2.0) |
| `linker` | integer | No | dustmasker linker (default 1) |
| `engine` | string | No | `sdust` (default) or `dustmasker` — exact NCBI dustmasker 2.12.0 parity (IUPAC scanned as bases; N runs > window and leading/trailing N runs cut the scan and are masked; window 8–64, 10·threshold an integer 2–64) |

## Output Schema

`items`: `start, end` (0-based half-open), `length`

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1003 | windowSize >= 3; threshold finite >= 0; linker 1-32 |
| 1003 | engine `dustmasker`: windowSize outside 8–64, 10·threshold not an integer 2–64, non-IUPAC input; unknown engine name |

## Examples

### Example 1: sdust -w 64 -t 20

**Input:** `{"sequence": "ACGTGCATGCAAAAAAAAAAAAAAAAGCTAGCATCGACTGCAGCACACACACACACACAGATCGATCGTACGGTGCATGACAAAAAAAAAAAAACT"}`

**Output:**

```json
{"items": [{"start": 10, "end": 26, "length": 16}, {"start": 43, "end": 59, "length": 16}, {"start": 81, "end": 94, "length": 13}]}
```

### Example 2: dustmasker -linker 32

**Input:** `{"sequence": "ACGTGCATGCAAAAAAAAAAAAAAAAGCTAGCATCGACTGCAGCACACACACACACACAGATCGATCGTACGGTGCATGACAAAAAAAAAAAAACT", "linker": 32}`

**Output:**

```json
{"items": [{"start": 10, "end": 94, "length": 84}]}
```

### Example 3: N splits the scan

**Input:** `{"sequence": "ACGTNNAAAAAAAAAAAANACGTACACACACACACACANNGGGCCCTAGGTCA"}`

**Output:**

```json
{"items": [{"start": 6, "end": 18, "length": 12}, {"start": 23, "end": 38, "length": 15}]}
```

### Example 4: dustmasker engine, linker 5

**Input:** `{"sequence": "NNACGTTGCAAAAAAAAAAAACGTRRRRRRRRRRRRTGCANNNNNNNNNNACGTTGCAANNNN", "windowSize": 8, "linker": 5, "engine": "dustmasker"}`

**Output** (= dustmasker 2.12.0 `-window 8 -level 20 -linker 5 -outfmt interval` `0 - 1`, `9 - 49`, `59 - 62`; the 10-N run joins
because 35 + linker equals its start exactly — dustmasker's `s_InsertMerge`; linker 6 leaves `9 - 35`, `40 - 49` apart):

```json
{"items": [{"start": 0, "end": 2, "length": 2}, {"start": 9, "end": 50, "length": 41}, {"start": 59, "end": 63, "length": 4}]}
```

## Performance

- O(n · windowSize) worst case (SDUST perfect-interval list).

## See Also

- [mask_low_complexity](mask_low_complexity.md)
- [dust_score](dust_score.md)
- [find_longdust_regions](find_longdust_regions.md)
