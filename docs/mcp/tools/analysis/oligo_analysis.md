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

RSAT `oligo-analysis` (van Helden et al. 1998) over-representation of the length-k words of one DNA sequence (the RSAT overload of `MotifFinder.DiscoverMotifs`): for every word (or reverse-complement pair with `strands` = `both`, RSAT `-2str`) with occ ≥ `minCount`, the expected frequency under the background model, exp_occ, the ratio occ/exp_occ and the binomial right-tail significance occ_P = P(X ≥ occ), X ~ Bin(N−k+1, exp_freq), occ_E = occ_P × tested patterns, occ_sig = −log10 occ_E. Backgrounds: `equiprobable` (`-bg equi`), `input` (Bernoulli from input, RSAT default), `bernoulli` (given A,C,G,T), `markov` (order m from input, `-markov m`), `markov_table` (RSAT `-bgfile` oligo frequencies with pseudo-frequency). `countOverlapping` = false is RSAT `-noov`. Reproduces the RSAT perl code (≤ 1e-12). Values beyond the double range (`ratio`, `possibleOligos` for very long k) are returned as null. RSAT options (routed to `MotifFinder.AnalyzeOligos`): `extraSequences` (multi-sequence input), `zscore` (RSAT `-return zscore`: exp_var = n·p·(2·ovlp − 1 − (2k+1)·p), or exp_occ with `-noov`; ovlp = Pevzner overlap coefficient), `expectedFrequencyPseudo` (RSAT `-pseudo`: exp_freq ← (1 − ψ)·exp_freq + ψ/NPO per strand), `degenerate` (`oneN` / `onedeg`: one IUPAC position, occ and exp_freq summed over matching words), background `lexicon` (RSAT `-lexicon` segmentation frequencies), `calibrationTable` + `calibrationMode` (RSAT `-calibN` / `-calib1`: P(occ ≤ X ≤ n) under a negative binomial when mean < variance, else Poisson). Patterns without a valid expected frequency (RSAT "NA") are omitted in this mode.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L58](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L58)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence |
| `k` | integer | No | Oligonucleotide length (>= 1, default 6) |
| `minCount` | integer | No | Occurrence threshold (RSAT -lth occ, default 2) |
| `background` | string | No | Background model: 'input' (Bernoulli from input, default), 'equiprobable', 'bernoulli' (needs residueFrequencies), 'markov' (order markovOrder from input), 'markov_table' (needs oligoFrequencies), 'lexicon' (RSAT -lexicon, k >= 2) |
| `markovOrder` | integer | No | Markov order for background 'markov' (>= 0; <= k-2 when > 0; default 1) |
| `residueFrequencies` | array<number> | No | Residue probabilities A,C,G,T for background 'bernoulli' |
| `oligoFrequencies` | object | No | (m+1)-mer -> frequency table for background 'markov_table' (RSAT -bgfile oligos format) |
| `pseudoFrequency` | number | No | Pseudo-frequency in [0,1] for background 'markov_table' (default 0.01) |
| `strandInsensitive` | boolean | No | The 'markov_table' frequencies are strand-insensitive pair frequencies (default false) |
| `strands` | string | No | 'single' (-1str, default) or 'both' (-2str) |
| `countOverlapping` | boolean | No | Count overlapping occurrences (default true); false = -noov |
| `extraSequences` | array<string> | No | Additional DNA sequences analysed together with 'sequence' (RSAT multi-sequence input) |
| `zscore` | boolean | No | Also report expectedVariance, overlapCoefficient and zScore (RSAT -return zscore; default false) |
| `expectedFrequencyPseudo` | number | No | RSAT -pseudo: pseudo-frequency psi in [0,1] on the expected frequencies (default 0 = off) |
| `degenerate` | string | No | 'none' (default), 'oneN' (RSAT -oneN) or 'onedeg' (RSAT -onedeg: one of R Y W S M K H B V D N per word) |
| `calibrationTable` | string | No | RSAT calibration file text (-calibN / -calib1; word [id] mean sd variance ...) |
| `calibrationMode` | string | No | 'set' (-calibN, default) or 'sequence' (-calib1: mean and variance × number of sequences) |

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
| `motifs[].sequenceIndices`, `overlaps`, `observedFrequency`, `fittedDistribution` | — | RSAT-option runs: sequence index per position, -noov overlaps (ovl_occ), obs_freq, Binomial / Poisson / NegativeBinomial |
| `motifs[].expectedVariance`, `overlapCoefficient`, `zScore` | number|null | With `zscore` (expectedVariance also with a calibration): RSAT exp_var, ovlp, z-score |
| `motifs[].lexiconSegmentation` | string|null | Best lexicon split `prefix|suffix` (background `lexicon`) |
| `sequenceCount`, `degenerate` | integer / string | RSAT-option runs: number of sequences, degenerate mode |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1002 | k must be >= 1 |
| 1003 | Background must be 'input', 'equiprobable', 'bernoulli', 'markov', 'markov_table' or 'lexicon' |
| 1004 | Background 'bernoulli' requires 4 residueFrequencies (A,C,G,T) |
| 1005 | Background 'markov_table' requires oligoFrequencies |
| 1006 | Strands must be 'single' or 'both' |
| 1007 | degenerate must be 'none', 'oneN' or 'onedeg'. |
| 1008 | expectedFrequencyPseudo must be in [0, 1]. |
| 1009 | calibrationMode must be 'set' or 'sequence'. / Calibration line N: expected <pattern> [id] <mean> <sd> <variance> |
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

### Example 3: RSAT -l 4 -1str -bg equi -lth occ 5 -return occ,proba,zscore

**Tool Call:**
```json
{
  "tool": "oligo_analysis",
  "arguments": {
    "sequence": "ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC",
    "k": 4,
    "minCount": 5,
    "background": "equiprobable",
    "zscore": true
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
      "occurrenceSignificance": 6.527391496197488,
      "sequenceIndices": [
        0,
        0,
        0,
        0,
        0,
        0
      ],
      "overlaps": 0,
      "observedFrequency": 0.1,
      "expectedVariance": 0.22613525390625006,
      "overlapCoefficient": 1,
      "zScore": 12.124455829017545,
      "fittedDistribution": "Binomial",
      "lexiconSegmentation": null
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
      "occurrenceSignificance": 5.080539500963187,
      "sequenceIndices": [
        0,
        0,
        0,
        0,
        0
      ],
      "overlaps": 0,
      "observedFrequency": 0.08333333333333333,
      "expectedVariance": 0.22613525390625006,
      "overlapCoefficient": 1,
      "zScore": 10.021569181166264,
      "fittedDistribution": "Binomial",
      "lexiconSegmentation": null
    }
  ],
  "oligoLength": 4,
  "strands": "single",
  "countOverlapping": true,
  "totalOccurrences": 60,
  "testedPatterns": 2,
  "possibleOligos": 256,
  "sequenceCount": 1,
  "degenerate": "none"
}
```

### Example 4: RSAT -l 3 -1str -calibN (negative binomial)

**Tool Call:**
```json
{
  "tool": "oligo_analysis",
  "arguments": {
    "sequence": "ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC",
    "k": 3,
    "minCount": 6,
    "calibrationTable": "; pattern\tavg\tstd\tvar\natg\t0.9\t1.2\t1.44\ntgc\t1.5\t1.1\t1.21\n"
  }
}
```

**Response:**
```json
{
  "motifs": [
    {
      "sequence": "ATG",
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
      "expectedFrequency": 0.014754098360655738,
      "expectedOccurrences": 0.9,
      "ratio": 6.666666666666666,
      "occurrenceProbability": 0.006714235225703919,
      "occurrenceEValue": 0.01342847045140784,
      "occurrenceSignificance": 1.8719734521372047,
      "sequenceIndices": [
        0,
        0,
        0,
        0,
        0,
        0
      ],
      "overlaps": 0,
      "observedFrequency": 0.09836065573770492,
      "expectedVariance": 1.44,
      "overlapCoefficient": null,
      "zScore": null,
      "fittedDistribution": "NegativeBinomial",
      "lexiconSegmentation": null
    },
    {
      "sequence": "TGC",
      "reverseComplement": null,
      "count": 6,
      "positions": [
        1,
        5,
        9,
        25,
        40,
        44
      ],
      "expectedFrequency": 0.02459016393442623,
      "expectedOccurrences": 1.5,
      "ratio": 4,
      "occurrenceProbability": 0.004455980775247846,
      "occurrenceEValue": 0.008911961550495694,
      "occurrenceSignificance": 2.0500266958736195,
      "sequenceIndices": [
        0,
        0,
        0,
        0,
        0,
        0
      ],
      "overlaps": 0,
      "observedFrequency": 0.09836065573770492,
      "expectedVariance": 1.21,
      "overlapCoefficient": null,
      "zScore": null,
      "fittedDistribution": "Poisson",
      "lexiconSegmentation": null
    }
  ],
  "oligoLength": 3,
  "strands": "single",
  "countOverlapping": true,
  "totalOccurrences": 61,
  "testedPatterns": 2,
  "possibleOligos": 64,
  "sequenceCount": 1,
  "degenerate": "none"
}
```
## Worked Example

occ_P for ATGC (6 of 60 windows, p = 1/256) equals scipy `binom.sf(5, 60, 1/256)` = 1.4845e-07 — the RSAT value; occ_E multiplies it by the number of tested patterns. With `zscore`, ATGC has exp_var = 60·(1/256)·(2·1 − 1 − 9/256) = 0.22613525390625 and z = (6 − 0.234375)/√0.22613525390625 = 12.124455829017547 (RSAT 12.124455829017547). In the calibration example ATG (mean 0.9 < variance 1.44) uses the negative binomial p = 0.6, size 1.5: P(6 ≤ X ≤ 61) = 0.0067142352257039 (mpmath; RSAT prints 0.0067141894 from 5-digit terms).

## See Also

- [discover_motifs](discover_motifs.md) — Observed/expected enrichment only
- [shared_motifs_significance](shared_motifs_significance.md) — Matching-sequence significance across sequences
- [dyad_analysis](dyad_analysis.md) — RSAT dyad-analysis (spaced dyads)

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L58](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L58)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
