# dyad_analysis

RSAT dyad-analysis: over-represented spaced dyads (monad pairs with a fixed spacing).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `dyad_analysis` |
| **Method ID** | `OverrepresentedWordAnalyzer.AnalyzeDyads` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

RSAT dyad-analysis: over-represented spaced dyads M1 n{s} M2 (two monads of length m separated by s unspecified residues, s in [minSpacing, maxSpacing]) in DNA sequences — any dyad, direct repeats (dr), inverted repeats (ir) or both (rep); expected frequency from the input monad frequencies (or a monad / dyad frequency table), expected occurrences, z-score, and binomial significance occ_P / occ_E / occ_sig over the T_s possible positions. RSAT defaults: m = 3, spacing 0-20, both strands, non-overlapping counts.

## Core Documentation Reference

- Source: [Seqeron.Genomics.OverrepresentedWords/OverrepresentedWordAnalyzer.DyadAnalysis.cs#L56](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.OverrepresentedWords/OverrepresentedWordAnalyzer.DyadAnalysis.cs#L56)
- Algorithm: [Overrepresented_Kmer_Discovery.md](../../../algorithms/Motif_Discovery/Overrepresented_Kmer_Discovery.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | DNA sequences. |
| `monadLength` | integer | No | Monad length m (RSAT -l, >= 1, default 3). (default: 3) |
| `minSpacing` | integer | No | Minimal spacing (>= 0, default 0). (default: 0) |
| `maxSpacing` | integer | No | Maximal spacing (>= minSpacing, default 20). (default: 20) |
| `dyadType` | string | No | Dyad type: 'any' (default), 'dr' (direct repeats), 'ir' (inverted repeats) or 'rep' (both). |
| `strands` | string | No | Strands: 'both' (-2str, default; reverse-complement pairs) or 'single' (-1str). (default: both) |
| `countOverlapping` | boolean | No | Count overlapping occurrences (-ovlp); default false = RSAT -noov. (default: false) |
| `minCount` | integer | No | Occurrence threshold (RSAT -lth occ, default 1); <= 0 tests every combination of observed monads (RSAT without a threshold). (default: 1) |
| `background` | string | No | Background: 'monads' (default, input monad frequencies), 'monad_table' (needs monadFrequencies, RSAT -mncf) or 'dyad_table' (needs dyadFrequencies, RSAT -expfreq). |
| `monadFrequencies` | object | No | Monad -> expected frequency, for background 'monad_table'. |
| `dyadFrequencies` | object | No | Dyad (RSAT notation, e.g. 'acgn{4}tta') -> expected single-strand frequency, for background 'dyad_table'. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `dyads` | array<object> | Dyads with occ >= minCount: `pattern, firstMonad, spacing, secondMonad, reverseComplement, occurrences, overlaps, sequenceIndices, positions, observedFrequency, expectedFrequency, expectedOccurrences, expectedVariance, overlapCoefficient, zScore, ratio, occurrenceProbability, occurrenceEValue, occurrenceSignificance, isDirectRepeat, isReversePalindrome, expectedFromMonads` (null = RSAT NA / not tested) |
| `monadLength` | integer | m |
| `minSpacing` | integer | Minimal spacing |
| `maxSpacing` | integer | Maximal spacing |
| `dyadType` | string | Dyad type |
| `strands` | string | 'single' or 'both' |
| `countOverlapping` | boolean | Overlap mode |
| `sequenceCount` | integer | Number of sequences |
| `monadOccurrences` | integer | Monad windows (RSAT sum_oligo_count) |
| `spacings` | array<object> | Per spacing `{spacing, possiblePositions (T_s), occurrences (occ_sum), overlaps (ovl_sum)}` |
| `testedPatterns` | integer | Dyads with a P-value (occ_E multiplier) |
| `possibleDyads` | number|null | RSAT nb_possible_dyads |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required |
| 1001 | Invalid DNA sequence |
| 1002 | monadLength must be >= 1. |
| 1002 | Spacings must satisfy 0 <= minSpacing <= maxSpacing. |
| 1003 | dyadType must be 'any', 'dr', 'ir' or 'rep'. |
| 1004 | Background must be 'monads', 'monad_table' or 'dyad_table'. |
| 1005 | Background 'monad_table' requires monadFrequencies. |
| 1005 | Background 'dyad_table' requires dyadFrequencies. |
| 1006 | Strands must be 'single' or 'both'. |

## Examples

### Example 1: RSAT dyad-analysis -l 3 -sp 0-2 -lth occ 3 (defaults -2str -noov)

**User Prompt:**
> Which spaced trinucleotide pairs (spacing 0–2) are over-represented in this sequence, both strands?

**Tool Call:**
```json
{
  "tool": "dyad_analysis",
  "arguments": {
    "sequences": [
      "ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC"
    ],
    "monadLength": 3,
    "minSpacing": 0,
    "maxSpacing": 2,
    "minCount": 3
  }
}
```

**Response:**
```json
{
  "dyads": [
    {
      "pattern": "CATn{0}GCA",
      "firstMonad": "CAT",
      "spacing": 0,
      "secondMonad": "GCA",
      "reverseComplement": "TGCn{0}ATG",
      "occurrences": 3,
      "overlaps": 3,
      "sequenceIndices": [
        0,
        0,
        0
      ],
      "positions": [
        1,
        7,
        38
      ],
      "observedFrequency": 0.05172413793103448,
      "expectedFrequency": 0.015049717817790915,
      "expectedOccurrences": 0.7524858908895458,
      "expectedVariance": 0.7089265559961291,
      "overlapCoefficient": 1.00390625,
      "zScore": 2.6693271529711025,
      "ratio": 3.986785714285715,
      "occurrenceProbability": 0.057099595708850225,
      "occurrenceEValue": 0.057099595708850225,
      "occurrenceSignificance": 1.2433669667460436,
      "isDirectRepeat": false,
      "isReversePalindrome": false,
      "expectedFromMonads": true
    }
  ],
  "monadLength": 3,
  "minSpacing": 0,
  "maxSpacing": 2,
  "dyadType": "any",
  "strands": "both",
  "countOverlapping": false,
  "sequenceCount": 1,
  "monadOccurrences": 61,
  "spacings": [
    {
      "spacing": 0,
      "possiblePositions": 58,
      "occurrences": 50,
      "overlaps": 8
    },
    {
      "spacing": 1,
      "possiblePositions": 57,
      "occurrences": 48,
      "overlaps": 9
    },
    {
      "spacing": 2,
      "possiblePositions": 56,
      "occurrences": 53,
      "overlaps": 3
    }
  ],
  "testedPatterns": 1,
  "possibleDyads": 6240
}
```

## Worked Example

CATn{0}GCA|TGCn{0}ATG: 3 counted occurrences (3 overlapping ones discarded by -noov) among T_0 = 58 positions; exp_freq = f(CAT)·f(GCA) + f(ATG)·f(TGC) = (5·4 + 6·6)/61² = 0.015049717817790915; exp_occ = exp_freq · 50 = 0.7524858908895458; occ_P = P(Bin(58, 0.01505) ≥ 3) = 0.057099595708850225. The RSAT run of the same command (dyad-analysis 1.78) gives exp_freq 0.015049717817790915, occ_P 0.0570995957088501, 1 tested; on 200 single-strand runs (14,746 dyads) every statistic agrees with RSAT to ≤ 1.3e-14 (z-scores to RSAT's 2 printed decimals).

## See Also

- [oligo_analysis](oligo_analysis.md) — RSAT oligo-analysis (contiguous words)
- [shared_motifs_significance](shared_motifs_significance.md) — Matching-sequence significance

## References

- Algorithm source: [Seqeron.Genomics.OverrepresentedWords/OverrepresentedWordAnalyzer.DyadAnalysis.cs#L56](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.OverrepresentedWords/OverrepresentedWordAnalyzer.DyadAnalysis.cs#L56)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
