# find_pre_mirna_hairpins

Identify pre-miRNA hairpin candidates in a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_pre_mirna_hairpins` |
| **Method ID** | `MiRnaAnalyzer.FindPreMiRnaHairpins` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Folds every window in `[minHairpinLength, maxHairpinLength]` with the canonical Zuker–Stiegler MFE engine
(Turner 2004) and accepts it as a pre-miRNA when the MFE structure is a single hairpin with ≥ 16 base pairs
(Ambros 2003), a 3–25 nt apical loop (Bartel 2004) and MFEI ≥ 0.85 (Zhang 2006). Bulges and internal loops are
tolerated (hsa-mir-21 and hsa-let-7a-1 are detected). Each candidate reports its span, the dot-bracket and ΔG°
of the MFE fold, and the mature (5p, `matureLength` nt from the outermost pair of the stem) and star
(opposite-arm product with the Dicer 2-nt 3′ overhang, miRDeep2 `find_star`) sequences. The tool is the
default of `MiRnaAnalyzer.FindPreMiRnaHairpins` (`minMfei` = Zhang 2006 MFEI threshold, default 0.85).

**Cost and the screen option.** The default folds every window (≈ 20 ms per 55–120-nt window, parallelised over
window starts), so scanning a long sequence is slow — keep the default for a locus of a few hundred nucleotides.
For long inputs set `useConsecutiveScreen = true`: the fast consecutive-pairing screen (window ends must pair
uninterrupted — Watson–Crick or G·U — for ≥ 18 bp with a 3–25 nt apical loop; no bulges, so e.g. hsa-mir-21 is
missed; blunt mature/star labels = first/last `matureLength` nt of the window; Structure/ΔG° from the canonical
stem-loop model; `minMfei` ignored). `useConsecutiveScreen`/`minMfei` are additive optional parameters: existing
calls are unchanged.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L2549](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L2549)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence to scan |
| `minHairpinLength` | integer | No | Minimum hairpin length in nt (default 55) |
| `maxHairpinLength` | integer | No | Maximum hairpin length in nt (default 120) |
| `matureLength` | integer | No | Mature miRNA length in nt (default 22) |
| `useConsecutiveScreen` | boolean | No | `false` (default) = MFE-fold detection; `true` = fast consecutive-pairing screen for long inputs |
| `minMfei` | number | No | MFEI threshold of the default mode (default 0.85, ≥ 0); ignored by the screen |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `hairpins[].start` | integer | Hairpin start |
| `hairpins[].end` | integer | Hairpin end |
| `hairpins[].sequence` | string | Full hairpin sequence |
| `hairpins[].matureSequence` | string | Extracted mature arm |
| `hairpins[].starSequence` | string | Extracted star arm |
| `hairpins[].structure` | string | Dot-bracket secondary structure |
| `hairpins[].freeEnergy` | number | Turner 2004 nearest-neighbour free energy |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | minMfei must be a finite value >= 0 |

## Examples

### Example 1: A 57 nt hairpin

A 25 bp-stem / 7 nt-loop / 25 bp-stem sequence with `minHairpinLength = 55` yields six accepted windows
(start 0–2, end 54–56); the full 57 nt window folds to `(((…(26 more)…...)))` with a 3-nt apical loop,
ΔG° −48.48 kcal/mol (ViennaRNA 2.8 −48.40).

**Response (abridged):**
```json
{ "hairpins": [ { "start": 0, "end": 54 }, { "start": 0, "end": 55 }, { "start": 0, "end": 56 }, { "start": 1, "end": 55 }, { "start": 1, "end": 56 }, { "start": 2, "end": 56 } ] }
```

### Example 2: Below minimum length

A 47 nt hairpin is rejected against a 55 nt minimum:

**Response:**
```json
{ "hairpins": [] }
```

## Performance

- **Time Complexity:** O(n · L · L³) for sequence length n and window length L (one MFE fold per window; with `useConsecutiveScreen = true` O(n · L²))
- **Space Complexity:** O(L²)

## See Also

- [create_mirna](create_mirna.md) — build a mature miRNA record
- [align_mirna_to_target](align_mirna_to_target.md) — miRNA-target duplex
