# oligo_analysis

RSAT oligo-analysis k-mer over-representation with binomial significance.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `oligo_analysis` |
| **Method ID** | `MotifFinder.DiscoverMotifs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

RSAT `oligo-analysis` (van Helden et al. 1998) over-representation of the length-k words of one DNA sequence (the RSAT overload of `MotifFinder.DiscoverMotifs`): for every word (or reverse-complement pair with `strands` = `both`, RSAT `-2str`) with occ ≥ `minCount`, the expected frequency under the background model, exp_occ, the ratio occ/exp_occ and the binomial right-tail significance occ_P = P(X ≥ occ), X ~ Bin(N−k+1, exp_freq), occ_E = occ_P × tested patterns, occ_sig = −log10 occ_E. Backgrounds: `equiprobable` (`-bg equi`), `input` (Bernoulli from input, RSAT default), `bernoulli` (given A,C,G,T), `markov` (order m from input, `-markov m`), `markov_table` (RSAT `-bgfile` oligo frequencies with pseudo-frequency). `countOverlapping` = false is RSAT `-noov`. Reproduces the RSAT perl code (≤ 1e-12). Values beyond the double range (`ratio`, `possibleOligos` for very long k) are returned as null.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L58](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L58)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence |
| `k` | integer | No | Oligonucleotide length (>= 1, default 6) |
| `minCount` | integer | No | Occurrence threshold (RSAT -lth occ, default 2) |
| `background` | string | No | Background model: 'input' (Bernoulli from input, default), 'equiprobable', 'bernoulli' (needs residueFrequencies), 'markov' (order markovOrder from input), 'markov_table' (needs oligoFrequencies) |
| `markovOrder` | integer | No | Markov order for background 'markov' (>= 0; <= k-2 when > 0; default 1) |
| `residueFrequencies` | array<number> | No | Residue probabilities A,C,G,T for background 'bernoulli' |
| `oligoFrequencies` | object | No | (m+1)-mer -> frequency table for background 'markov_table' (RSAT -bgfile oligos format) |
| `pseudoFrequency` | number | No | Pseudo-frequency in [0,1] for background 'markov_table' (default 0.01) |
| `strandInsensitive` | boolean | No | The 'markov_table' frequencies are strand-insensitive pair frequencies (default false) |
| `strands` | string | No | 'single' (-1str, default) or 'both' (-2str) |
| `countOverlapping` | boolean | No | Count overlapping occurrences (default true); false = -noov |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `motifs` | array<object> | Tested words in order of first occurrence: `sequence, reverseComplement (both strands), count, positions, expectedFrequency, expectedOccurrences, ratio, occurrenceProbability (occ_P), occurrenceEValue (occ_E), occurrenceSignificance (occ_sig)` |
| `oligoLength` | integer | k |
| `strands` | string | `single` or `both` |
| `countOverlapping` | boolean | Overlap mode |
| `totalOccurrences` | integer | Binomial trials N − k + 1 |
| `testedPatterns` | integer | Patterns tested (occ_E multiplier) |
| `possibleOligos` | number|null | RSAT nb_possible_oligos (null beyond the double range) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1002 | k must be >= 1 |
| 1003 | Background must be 'input', 'equiprobable', 'bernoulli', 'markov' or 'markov_table' |
| 1004 | Background 'bernoulli' requires 4 residueFrequencies (A,C,G,T) |
| 1005 | Background 'markov_table' requires oligoFrequencies |
| 1006 | Strands must be 'single' or 'both' |
| — | Markov order > k − 2 or a zero-probability observed word (library ArgumentException) |

## Examples

### Example 1: RSAT -l 4 -1str -bg equi -lth occ 5

**User Prompt:**
> Which 4-mers are over-represented in this sequence (equiprobable background)?

**Tool Call:**
```json
{
  "tool": "oligo_analysis",
  "arguments": {
    "sequence": "ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC",
    "k": 4,
    "minCount": 5,
    "background": "equiprobable"
  }
}
```

**Response:**
```json
{
  "motifs": [
    {
      "sequence": "ATGC",
      "reverseComplement": null,
      "count": 6,
      "positions": [
        0,
        4,
        8,
        24,
        39,
        43
      ],
      "expectedFrequency": 0.003906250000000001,
      "expectedOccurrences": 0.23437499999999997,
      "ratio": 25.6,
      "occurrenceProbability": 1.4844942103072874e-07,
      "occurrenceEValue": 2.968988420614575e-07,
      "occurrenceSignificance": 6.527391496197488
    },
    {
      "sequence": "CATG",
      "reverseComplement": null,
      "count": 5,
      "positions": [
        3,
        7,
        23,
        38,
        42
      ],
      "expectedFrequency": 0.003906250000000001,
      "expectedOccurrences": 0.23437499999999997,
      "ratio": 21.33333333333334,
      "occurrenceProbability": 4.153655783353449e-06,
      "occurrenceEValue": 8.307311566706898e-06,
      "occurrenceSignificance": 5.080539500963187
    }
  ],
  "oligoLength": 4,
  "strands": "single",
  "countOverlapping": true,
  "totalOccurrences": 60,
  "testedPatterns": 2,
  "possibleOligos": 256
}
```

### Example 2: RSAT -l 4 -2str -lth occ 6

**Tool Call:**
```json
{
  "tool": "oligo_analysis",
  "arguments": {
    "sequence": "ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC",
    "k": 4,
    "minCount": 6,
    "strands": "both"
  }
}
```

**Response:**
```json
{
  "motifs": [
    {
      "sequence": "ATGC",
      "reverseComplement": "GCAT",
      "count": 9,
      "positions": [
        0,
        2,
        4,
        6,
        8,
        24,
        39,
        41,
        43
      ],
      "expectedFrequency": 0.007777109332017005,
      "expectedOccurrences": 0.4666265599210202,
      "ratio": 19.287371900826464,
      "occurrenceProbability": 1.076064733019135e-09,
      "occurrenceEValue": 1.076064733019135e-09,
      "occurrenceSignificance": 8.968161601952728
    }
  ],
  "oligoLength": 4,
  "strands": "both",
  "countOverlapping": true,
  "totalOccurrences": 60,
  "testedPatterns": 1,
  "possibleOligos": 136
}
```

## Worked Example

occ_P for ATGC (6 of 60 windows, p = 1/256) equals scipy `binom.sf(5, 60, 1/256)` = 1.4845e-07 — the RSAT value; occ_E multiplies it by the number of tested patterns.

## See Also

- [discover_motifs](discover_motifs.md) — Observed/expected enrichment only
- [shared_motifs_significance](shared_motifs_significance.md) — Matching-sequence significance across sequences

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L58](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L58)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
