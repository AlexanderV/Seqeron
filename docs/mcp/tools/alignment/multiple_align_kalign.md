# multiple_align_kalign

Kalign 3.3.5's multiple sequence alignment, ported from the program.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_kalign` |
| **Method ID** | `SequenceAligner.KalignRun` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Kalign 3.3.5's multiple alignment (Lassmann 2019, Bioinformatics 36:1928; kalign -i in.fa -o out.fa), ported step by step from the program's source: the alphabet (DNA / protein) is detected from the letter counts (non-letters are gaps and are removed, sequences without residues are dropped, letter case is kept); distances = the bit-parallel edit distance of Myers (1999) of the shorter sequence (its first 255 residues) against the best window of the longer one, on the 13-letter reduced protein alphabet or ACGTN; guide tree = UPGMA (average linkage + 0.001) for fewer than 100 sequences, otherwise bisecting 2-means on the distances to 32 anchors; progressive alignment of profiles with the affine three-state dynamic program in Hirschberg's linear space (position-specific gap penalties, terminal gap extension), matrices RNA (the default for nucleotides: gap open 217, extension 39.4, terminal 292.6), dna (match 5 / mismatch -4, 8 / 6 / 0), CorBLOSUM66_13plus (default for proteins: 5.5 / 2 / 1) or GON250 (type divergent: 55 / 8 / 4). Quirks of 3.3.5 that the port reproduces and the oracle confirms: the terminal gap extension option is applied only together with a gap extension option (tgpe alone is ignored; gpe alone sets it to -1); the gap-open penalty is not charged for columns inserted into a profile; with threads = 1 the serial recursion is used (no shorter-first swap), with more threads the OpenMP recursion. VERIFICATION: rows equal the installed Debian kalign 1:3.3.5-1 (AVX2 build) on 526 / 526 committed runs (DNA / RNA / protein, every type and gap-penalty option, 1 / 2 / 4 threads, sequences of more than 500 columns, 100 - 1000 sequences, refusals). The bit-parallel distance of the portable (non-AVX2) build limits the pattern to 1 024 instead of 255 residues; sequences longer than 255 residues can therefore differ between builds. totalScore is the natural affine sum of pairs (BLOSUM62 10 / 0.2 for proteins, EDNAFULL 10 / 0.5 for nucleotides; 0 / NaN for thousands of sequences).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.Kalign.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Kalign.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section Kalign 3.3.5)
- Oracle: `scripts/alignment/generate_kalign_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/kalign_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences to align (at least two with residues): letters (the alphabet of the detected type), any other character is a gap and is removed. |
| `type` | string | No | Kalign --type: 'auto' (default), 'rna', 'dna', 'protein', 'divergent' (GON250, proteins only) or 'dna-internal' (library only). A type that does not fit the detected alphabet is refused like the program does. |
| `gapOpen` | number | No | Kalign --gpo: gap open penalty (positive = cost); omit or negative = the default of the type. |
| `gapExtend` | number | No | Kalign --gpe: gap extension penalty; omit or negative = the default of the type. |
| `terminalGapExtend` | number | No | Kalign --tgpe: terminal gap extension penalty; applied only when gapExtend is also given (a quirk of 3.3.5). |
| `threads` | integer | No | Kalign -n (at least 1, default 4): only 1 versus more than 1 matters (serial versus OpenMP recursion). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | Rows of the aligned sequences in input order (letters keep their case, gaps are -). |
| `consensus` | string | Consensus string of the alignment. |
| `totalScore` | integer | Natural affine sum of pairs (BLOSUM62 10 / 0.2 or EDNAFULL 10 / 0.5) rounded; 0 for thousands of sequences. |
| `exactTotalScore` | number | The same score unrounded (NaN when it cannot be computed). |
| `keptIndices` | array<integer> | Indices of the input sequences that were aligned (sequences without residues are dropped). |
| `alphabet` | string | 'dna' or 'protein': the alphabet detected from the letter counts. |
| `guideTree` | array<array<integer>> | Progressive alignment steps [a, b, c] in execution order: profiles a and b merged into c (0..n-1 sequences, n.. internal nodes). |

## Errors

| Code | Message |
|------|---------|
| 1001 | sequences is null, fewer than two sequences keep a residue, or a character is outside ASCII |
| 1002 | A letter that is not in the alphabet of the detected type (the program leaves its code undefined), or a type that does not fit the detected alphabet (the program stops with an error) |
| 1003 | threads < 1 or an unknown type word |

## Examples

### Example 1: Three DNA sequences, one thread (kalign -n 1)

**Tool Call:**
```json
{
  "tool": "multiple_align_kalign",
  "arguments": {
    "sequences": [
      "ACGTACGTTAGCTAGCTAGCATCGATCGATGCATGCATCG",
      "ACGTACGTAGCTAGCTAGCATCGATCGATGCATGCCATCG",
      "ACGTTCGTTAGCTAGCTTGCATCGATCGATGCATGCATCG"
    ],
    "threads": 1
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGTACGTTAGCTAGCTAGCATCGATCGATGCATG-CATCG",
    "ACGTACG-TAGCTAGCTAGCATCGATCGATGCATGCCATCG",
    "ACGTTCGTTAGCTAGCTTGCATCGATCGATGCATG-CATCG"
  ],
  "consensus": "ACGTACGTTAGCTAGCTAGCATCGATCGATGCATG-CATCG",
  "totalScore": 514,
  "exactTotalScore": 514,
  "keptIndices": [
    0,
    1,
    2
  ],
  "alphabet": "dna",
  "guideTree": [
    [
      0,
      1,
      3
    ],
    [
      3,
      2,
      4
    ]
  ]
}
```

`kalign -i in.fa -n 1 -q` writes exactly these three rows (the serial recursion; with the default 4 threads the first row ends `GCATGC-ATCG`).

### Example 2: Four proteins, GON250 matrix with a user gap-open penalty

**Tool Call:**
```json
{
  "tool": "multiple_align_kalign",
  "arguments": {
    "sequences": [
      "MKVLAAGIVGLLLAQWTR",
      "MKVLAGIVALLAQWTR",
      "MKVAAGIVGLLAQWTK",
      "MKLAAGIVGLLAQWR"
    ],
    "type": "divergent",
    "gapOpen": 40,
    "gapExtend": 6,
    "terminalGapExtend": 3
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "MKVLAAGIVGLLLAQWTR",
    "MK-VLAGIVALL-AQWTR",
    "MK-VAAGIVGLL-AQWTK",
    "MK-LAAGIVGLL-AQWR-"
  ],
  "consensus": "MK-LAAGIVGLL-AQWTR",
  "totalScore": 306,
  "exactTotalScore": 306,
  "keptIndices": [
    0,
    1,
    2,
    3
  ],
  "alphabet": "protein",
  "guideTree": [
    [
      2,
      3,
      4
    ],
    [
      0,
      4,
      5
    ],
    [
      5,
      1,
      6
    ]
  ]
}
```

`kalign -i in.fa --type divergent --gpo 40 --gpe 6 --tgpe 3 -q`.

## See Also

- [kalign_command](kalign_command.md)
- [multiple_align_muscle](multiple_align_muscle.md)
- [multiple_align_clustalw](multiple_align_clustalw.md)
- [multiple_align_tcoffee](multiple_align_tcoffee.md)
