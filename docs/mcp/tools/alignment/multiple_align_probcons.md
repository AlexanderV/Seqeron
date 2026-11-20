# multiple_align_probcons

ProbCons 1.12 (port of the program)

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_probcons` |
| **Method ID** | `SequenceAligner.ProbConsRun` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

ProbCons 1.12's multiple alignment of protein sequences (Do, Mahabhashyam, Brudno & Batzoglou 2005, Genome Res 15:330; probcons in.fa), ported step by step from the program's source: pair-HMM posterior match probabilities of every pair (forward / backward in single precision with the program's log-sum approximation; a match state and two pairs of insertion states; emission tables of BLOCKS 5.0 / BLOSUM62, transitions trained on BAliBASE 2.0; upper and lower case alike, other characters at the table default), `consistency` rounds of the consistency transformation P'(x_i-y_j) = 1/|S| sum_z sum_k P(x_i-z_k) P(z_k-y_j) on the sparse matrices (cells >= 0.01), an expected-accuracy guide tree (UPGMA-like: merge the largest value, average weighted by the merged value), progressive maximum-expected-accuracy alignment of profiles, and `iterativeRefinement` rounds of randomized-partitioning refinement (the bipartitions come from glibc rand() at its default seed 1, as in the program). Optional: `preTraining` rounds of EM re-estimation of the transition probabilities on the input itself, `alignmentOrder` (rows in alignment order), the cosmetic gap penalties gapOpenPenalty / gapExtensionPenalty (<= 0; in the program they overflow a heap buffer whenever the two profiles differ in size, so only equal-size merges are verified), `cutoff` (subtracted from the posterior before the profile alignments), `parameters` (a parameter file as text, probcons -p) and `annotate` (the -annot column scores, 0 - 200). Letter case is kept. VERIFICATION: stdout, standard error (parameter set, guide tree, -v lines) and written files equal the installed Debian probcons 1.12-14 on 407 / 407 committed command lines (scripts/alignment/generate_probcons_oracle.py), see probcons_command. totalScore is the natural BLOSUM62 affine sum of pairs (the program reports no score).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.ProbCons.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.ProbCons.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section ProbCons 1.12)
- Oracle: `scripts/alignment/generate_probcons_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/probcons_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Protein sequences (letters, '*'; at least one). Letters outside the 20 amino acids have the program's default emissions. |
| `names` | array<string> | No | Sequence headers (the first word is the name in the guide tree and Clustal output); omit for S1, S2, ... |
| `consistency` | integer | No | probcons -c: passes of consistency transformation, 0 - 5 (default 2). |
| `iterativeRefinement` | integer | No | probcons -ir: passes of iterative refinement, 0 - 1000 (default 100). |
| `preTraining` | integer | No | probcons -pre: rounds of EM pre-training, 0 - 20 (default 0). |
| `alignmentOrder` | boolean | No | probcons -a: report the sequences in alignment order instead of input order. |
| `gapOpenPenalty` | number | No | probcons -go: cosmetic gap-open penalty (<= 0; 0 together with a zero extension penalty = none). |
| `gapExtensionPenalty` | number | No | probcons -ge: cosmetic gap-extension penalty (<= 0). |
| `cutoff` | number | No | probcons -co: value in [0, 1] subtracted from the posterior probabilities before the profile alignments. |
| `parameters` | string | No | probcons -p: the text of a parameter file (initial distribution, gap open, gap extension, alphabet, lower triangle of the match emissions, single emissions); omit for the built-in parameters. |
| `annotate` | boolean | No | probcons -annot: also return the per-column annotation score (0 - 200: twice the mean posterior of the aligned residue pairs of the column, in percent). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | Rows of the alignment in output order (input order unless alignmentOrder; letters keep their case, gaps are -). |
| `consensus` | string | Consensus string of the alignment. |
| `totalScore` | integer | Natural affine sum of pairs (BLOSUM62 10 / 0.2) rounded; 0 when a letter is outside the BLOSUM62 alphabet (the program reports no score). |
| `exactTotalScore` | number | The same score unrounded (NaN when it cannot be computed). |
| `order` | array<integer> | Input index of every row. |
| `guideTree` | string | The tree the program prints to standard error: ((a (c d)) b), leaves are the sequence headers. |
| `annotation` | array<integer> | With annotate: the -annot score of every column (0 - 200); omitted otherwise. |

## Errors

| Code | Message |
|------|---------|
| 1001 | sequences is null or empty, a sequence is empty, or a character is not a letter, * or - (the program: "Unknown character encountered") |
| 1002 | consistency outside 0 - 5, iterativeRefinement outside 0 - 1000, preTraining outside 0 - 20, a positive gap penalty, cutoff outside [0, 1] |
| 1003 | parameters is not a parameter file (fewer than four lines) |

## Examples

### Example 1: Four proteins with the program defaults (probcons in.fa)

**Tool Call:**
```json
{
  "tool": "multiple_align_probcons",
  "arguments": {
    "sequences": [
      "MKVLAAGIVGLLLAQPSTAEEKV",
      "MKVLAGIVGLLAQPSTAEEKVLL",
      "MKILAAGIAGLLLQPSAEEKV",
      "MRVLAAGLVGLLLAQSTAEEK"
    ],
    "names": [
      "a one",
      "b",
      "c",
      "d"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "MKVLAAGIVGLLLAQPSTAEEK--V",
    "MKVL-AGIVGLLA-QPSTAEEKVLL",
    "MKILAAGIAGLLL-QPS-AEEK--V",
    "MRVLAAGLVGLLL-AQSTAEEK---"
  ],
  "consensus": "MKVLAAGIVGLLL-QPSTAEEK--V",
  "totalScore": 339,
  "exactTotalScore": 339.2,
  "order": [
    0,
    1,
    2,
    3
  ],
  "guideTree": "((a one (c d)) b)",
  "annotation": null
}
```

`probcons in.fa` writes exactly these rows and prints `Alignment tree: ((a one (c d)) b)` to standard error.

### Example 2: No consistency, no refinement, column annotation (probcons -c 0 -ir 0 -annot ann.txt)

**Tool Call:**
```json
{
  "tool": "multiple_align_probcons",
  "arguments": {
    "sequences": [
      "MKVLAAGIVGLLLAQPSTAEEKV",
      "MKVLAGIVGLLAQPSTAEEKVLL",
      "MKILAAGIAGLLLQPSAEEKV",
      "MRVLAAGLVGLLLAQSTAEEK"
    ],
    "consistency": 0,
    "iterativeRefinement": 0,
    "annotate": true
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "MKVLAAGIVGLLLAQPSTAEEK--V",
    "MKVL-AGIVGLLAQP-STAEEKVLL",
    "MKILAAGIAGLLLQP-SA-EEK--V",
    "MRVLAAGLVGLLLAQ-ST-AEE--K"
  ],
  "consensus": "MKVLAAGIVGLLLAP-STAEEK--V",
  "totalScore": 276,
  "exactTotalScore": 276.4,
  "order": [
    0,
    1,
    2,
    3
  ],
  "guideTree": "((S1 (S3 S4)) S2)",
  "annotation": [
    99,
    99,
    98,
    95,
    99,
    75,
    96,
    99,
    99,
    99,
    92,
    86,
    79,
    68,
    65,
    0,
    83,
    70,
    99,
    62,
    64,
    65,
    0,
    0,
    59
  ]
}
```

The annotation values are the lines of the `-annot` file of the installed program.

## See Also

- [probcons_command](probcons_command.md)
- [multiple_align_tcoffee](multiple_align_tcoffee.md)
- [multiple_align_muscle](multiple_align_muscle.md)
