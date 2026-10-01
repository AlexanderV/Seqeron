# shared_motifs_significance

RSAT matching-sequence significance of k-mers shared by several sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `shared_motifs_significance` |
| **Method ID** | `MotifFinder.FindSharedMotifs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

RSAT `oligo-analysis -return mseq,proba` (the RSAT overload of `MotifFinder.FindSharedMotifs`): for each word (or reverse-complement pair with `strands` = `both`) present in at least `minSequences` of the S input sequences, the expected frequency, the expected number of matching sequences exp_ms = S·(1 − (1 − p)^(nb_pos/S)), and the binomial significance ms_P = P(X ≥ mseq), X ~ Bin(S, exp_ms/S), ms_E = ms_P × possible oligos, ms_sig = −log10 ms_E. Same background models as `oligo_analysis`. Values beyond the double range are returned as null. `degenerate` = `oneN` / `onedeg` (RSAT `-oneN` / `-onedeg` with `-return mseq,proba`; `MotifFinder.FindSharedMotifs(…, pseudoFrequency, OligoDegeneracy)`): words with one N or one IUPAC code (R Y W S M K H B V D N); a sequence matches a degenerate word when one of its k-mers matches it (mseq = union of the matching sequences of its words), exp_freq = Σ over the matching words, NPO = k·|codes|·4^(k−1) (pairs with `-2str`).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L140](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L140)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | DNA sequences |
| `k` | integer | No | Word length (>= 1, default 6) |
| `minSequences` | integer | No | Matching-sequence quorum (>= 1, default 2) |
| `background` | string | No | Background model: 'input' (Bernoulli from input, default), 'equiprobable', 'bernoulli' (needs residueFrequencies), 'markov' (order markovOrder from input), 'markov_table' (needs oligoFrequencies) |
| `markovOrder` | integer | No | Markov order for background 'markov' (>= 0; <= k-2 when > 0; default 1) |
| `residueFrequencies` | array<number> | No | Residue probabilities A,C,G,T for background 'bernoulli' |
| `oligoFrequencies` | object | No | (m+1)-mer -> frequency table for background 'markov_table' (RSAT -bgfile oligos format) |
| `pseudoFrequency` | number | No | Pseudo-frequency in [0,1] for background 'markov_table' (default 0.01) |
| `strandInsensitive` | boolean | No | The 'markov_table' frequencies are strand-insensitive pair frequencies (default false) |
| `strands` | string | No | 'single' (default) or 'both' |
| `degenerate` | string | No | 'none' (default), 'oneN' (RSAT -oneN) or 'onedeg' (RSAT -onedeg: one of R Y W S M K H B V D N per word); a sequence matches a degenerate word when one of its k-mers matches it |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `motifs` | array<object> | Words meeting the quorum: `sequence, reverseComplement, sequenceIndices, prevalence, expectedFrequency, expectedMatchingSequences (exp_ms), matchingSequenceProbability (ms_P), matchingSequenceEValue (ms_E), matchingSequenceSignificance (ms_sig)` |
| `oligoLength` | integer | k |
| `strands` | string | `single` or `both` |
| `sequenceCount` | integer | S |
| `possiblePositions` | integer | nb_pos = Σ(Lᵢ − k + 1) |
| `possibleOligos` | number|null | ms_E multiplier (null beyond the double range) |
| `degenerate` | string|null | 'oneN' / 'onedeg' for degenerate runs, otherwise null |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required |
| 1001 | Invalid DNA sequence |
| 1002 | k must be >= 1 |
| 1003 | minSequences must be >= 1 |
| 1004 | Background must be 'input', 'equiprobable', 'bernoulli', 'markov' or 'markov_table' |
| 1005 | Strands must be 'single' or 'both' |
| 1006 | degenerate must be 'none', 'oneN' or 'onedeg'. |

## Examples

### Example 1: RSAT -l 4 -1str -lth mseq 3

**User Prompt:**
> Which 4-mers occur in at least 3 of these sequences, and how significant is that?

**Tool Call:**
```json
{
  "tool": "shared_motifs_significance",
  "arguments": {
    "sequences": [
      "ACGTACGTTAGC",
      "TTACGTAGCAAC",
      "GGTAGCACGTTT",
      "CATTTTACG"
    ],
    "k": 4,
    "minSequences": 3
  }
}
```

**Response:**
```json
{
  "motifs": [
    {
      "sequence": "ACGT",
      "reverseComplement": null,
      "sequenceIndices": [
        0,
        1,
        2
      ],
      "prevalence": 0.75,
      "expectedFrequency": 0.0037555250723974986,
      "expectedMatchingSequences": 0.12225827585979834,
      "matchingSequenceProbability": 0.00011159466135309416,
      "matchingSequenceEValue": 0.028568233306392098,
      "matchingSequenceSignificance": 1.5441166160753925
    },
    {
      "sequence": "TACG",
      "reverseComplement": null,
      "sequenceIndices": [
        0,
        1,
        3
      ],
      "prevalence": 0.75,
      "expectedFrequency": 0.0037555250723974986,
      "expectedMatchingSequences": 0.12225827585979834,
      "matchingSequenceProbability": 0.00011159466135309416,
      "matchingSequenceEValue": 0.028568233306392098,
      "matchingSequenceSignificance": 1.5441166160753925
    },
    {
      "sequence": "TAGC",
      "reverseComplement": null,
      "sequenceIndices": [
        0,
        1,
        2
      ],
      "prevalence": 0.75,
      "expectedFrequency": 0.0037555250723974986,
      "expectedMatchingSequences": 0.12225827585979834,
      "matchingSequenceProbability": 0.00011159466135309416,
      "matchingSequenceEValue": 0.028568233306392098,
      "matchingSequenceSignificance": 1.5441166160753925
    }
  ],
  "oligoLength": 4,
  "strands": "single",
  "sequenceCount": 4,
  "possiblePositions": 33,
  "possibleOligos": 256
}
```

### Example 2: RSAT -l 3 -2str -oneN -return mseq,proba -lth mseq 3

**Tool Call:**
```json
{
  "tool": "shared_motifs_significance",
  "arguments": {
    "sequences": ["ACGTACGGATCC", "ATGCATGAAC", "ACGATGTT"],
    "k": 3,
    "minSequences": 3,
    "strands": "both",
    "degenerate": "oneN"
  }
}
```

**Response:**
```json
{
  "motifs": [
    {
      "sequence": "ANG",
      "reverseComplement": "CNT",
      "sequenceIndices": [
        0,
        1,
        2
      ],
      "prevalence": 1,
      "expectedFrequency": 0.12444444444444448,
      "expectedMatchingSequences": 1.9639256569777699,
      "matchingSequenceProbability": 0.28055070694211726,
      "matchingSequenceEValue": 6.733216966610815,
      "matchingSequenceSignificance": -0.8282226091038667
    },
    {
      "sequence": "GTN",
      "reverseComplement": "NAC",
      "sequenceIndices": [
        0,
        1,
        2
      ],
      "prevalence": 1,
      "expectedFrequency": 0.12444444444444448,
      "expectedMatchingSequences": 1.9639256569777699,
      "matchingSequenceProbability": 0.28055070694211726,
      "matchingSequenceEValue": 6.733216966610815,
      "matchingSequenceSignificance": -0.8282226091038667
    },
    {
      "sequence": "NGA",
      "reverseComplement": "TCN",
      "sequenceIndices": [
        0,
        1,
        2
      ],
      "prevalence": 1,
      "expectedFrequency": 0.12444444444444448,
      "expectedMatchingSequences": 1.9639256569777699,
      "matchingSequenceProbability": 0.28055070694211726,
      "matchingSequenceEValue": 6.733216966610815,
      "matchingSequenceSignificance": -0.8282226091038667
    },
    {
      "sequence": "ATN",
      "reverseComplement": "NAT",
      "sequenceIndices": [
        0,
        1,
        2
      ],
      "prevalence": 1,
      "expectedFrequency": 0.14222222222222228,
      "expectedMatchingSequences": 2.120735413209536,
      "matchingSequenceProbability": 0.35326100605480804,
      "matchingSequenceEValue": 8.478264145315396,
      "matchingSequenceSignificance": -0.9283069431434835
    },
    {
      "sequence": "ANC",
      "reverseComplement": "GNT",
      "sequenceIndices": [
        0,
        1,
        2
      ],
      "prevalence": 1,
      "expectedFrequency": 0.12444444444444448,
      "expectedMatchingSequences": 1.9639256569777699,
      "matchingSequenceProbability": 0.28055070694211726,
      "matchingSequenceEValue": 6.733216966610815,
      "matchingSequenceSignificance": -0.8282226091038667
    },
    {
      "sequence": "GAN",
      "reverseComplement": "NTC",
      "sequenceIndices": [
        0,
        1,
        2
      ],
      "prevalence": 1,
      "expectedFrequency": 0.12444444444444448,
      "expectedMatchingSequences": 1.9639256569777699,
      "matchingSequenceProbability": 0.28055070694211726,
      "matchingSequenceEValue": 6.733216966610815,
      "matchingSequenceSignificance": -0.8282226091038667
    }
  ],
  "oligoLength": 3,
  "strands": "both",
  "sequenceCount": 3,
  "possiblePositions": 24,
  "possibleOligos": 24,
  "degenerate": "oneN"
}
```

Reference: RSAT 1.169 `-oneN`/`-onedeg` returns no rows (and its per-word mseq sum can exceed the number of sequences); the repaired oracle copy (`Degenerate` fixed, mseq = union of matching sequences) gives ANG|CNT exp_freq 0.12444444444444444, exp_ms 1.9639256569777699, ms_P 0.28055070694211726, NPO 24.

## Worked Example

ACGT is in sequences 0, 1, 2: RSAT reports ms_P 1.1159e-04 and ms_sig 1.5441 (ms_E = ms_P × 256).

## See Also

- [find_shared_motifs](find_shared_motifs.md) — Quorum only
- [oligo_analysis](oligo_analysis.md) — Single-sequence occurrence significance

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L140](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs#L140)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
