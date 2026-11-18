# multiple_align_tcoffee

The multiple alignment of the T-Coffee 13.41 program (t_coffee; Notredame, Higgins & Heringa 2000) reproduced end to end: its own primary library built by the program's pair-wise methods - 'proba_pair' (the program's default: ProbCons pair-HMM posteriors > 0.01, at most 4 x min(length) pairs per pair of sequences), 'slow_pair' (Myers-Miller global alignment, dna_idmat 10 / 0 with gap open -10 / extend -1 for DNA, BLOSUM62 for proteins) and 'lalign_id_pair' (Huang-Miller SIM, the lalignTopAlignments best non-intersecting local alignments), several methods merged the way the program's forked library merge does (which is not their union) - then library relaxation, the program's guide tree (k-tuple distances, neighbour joining with the trifurcating root turned into a binary tree; fast neighbour joining from 100 sequences) and the sparse linked-pair-wise progressive stage. Rows equal t_coffee -seq s.fa [-method m1,m2] -output fasta_aln on 1 500+ committed DNA and protein sets (libraries pair by pair on 800+ more), for any -n_core. The type (DNA when more than 98 % of the letters are A G C T U N R Y M K S W H B V D, otherwise protein) selects the matrices, penalties and models. Not available: proba_pair of RNA (a U). totalScore is the natural affine sum of pairs (EDNAFULL / BLOSUM62), not the program's objective. Sequence names are irrelevant here; the program mishandles input names of which one is a prefix of another.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_tcoffee` |
| **Method ID** | `SequenceAligner.MultipleAlignTCoffeeProgram` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The multiple alignment of the T-Coffee 13.41 program (t_coffee; Notredame, Higgins & Heringa 2000) reproduced end to end: its own primary library built by the program's pair-wise methods - 'proba_pair' (the program's default: ProbCons pair-HMM posteriors > 0.01, at most 4 x min(length) pairs per pair of sequences), 'slow_pair' (Myers-Miller global alignment, dna_idmat 10 / 0 with gap open -10 / extend -1 for DNA, BLOSUM62 for proteins) and 'lalign_id_pair' (Huang-Miller SIM, the lalignTopAlignments best non-intersecting local alignments), several methods merged the way the program's forked library merge does (which is not their union) - then library relaxation, the program's guide tree (k-tuple distances, neighbour joining with the trifurcating root turned into a binary tree; fast neighbour joining from 100 sequences) and the sparse linked-pair-wise progressive stage. Rows equal t_coffee -seq s.fa [-method m1,m2] -output fasta_aln on 1 500+ committed DNA and protein sets (libraries pair by pair on 800+ more), for any -n_core. The type (DNA when more than 98 % of the letters are A G C T U N R Y M K S W H B V D, otherwise protein) selects the matrices, penalties and models. Not available: proba_pair of RNA (a U). totalScore is the natural affine sum of pairs (EDNAFULL / BLOSUM62), not the program's objective. Sequence names are irrelevant here; the program mishandles input names of which one is a prefix of another.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.TCoffeeLibrary.cs#L352](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.TCoffeeLibrary.cs#L352)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section 4.10)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences to align (at least one, none empty): letters only (the program removes gaps; here they are rejected). |
| `methods` | array<string> | No | Primary-library methods in command-line order (t_coffee -method m1,m2): 'proba_pair' (default; = plain t_coffee seqs.fa), 'slow_pair', 'lalign_id_pair'. |
| `lalignTopAlignments` | integer | No | t_coffee -lalign_n_top (default 10): alignments reported by lalign_id_pair for every pair of sequences. (default: 10) |
| `relaxLibrary` | boolean | No | t_coffee -relax_lib: true = 1 (default, library relaxation), false = 0. (default: true) |
| `type` | string | No | 'auto' (default; the program's decision), 'dna' or 'protein' (forces the matrices, penalties and models). (default: auto) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA (input order). |
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a residue beats a gap on ties). |
| `totalScore` | integer | Natural affine sum-of-pairs score of the alignment (EDNAFULL 10 / 0.5 for DNA, BLOSUM62 10 / 0.2 for proteins; the program's own objective is not reported). |

## Errors

| Code | Message |
|------|---------|
| 1010 | a sequence is empty or contains a character other than a letter |
| 1002 | methods must be proba_pair, slow_pair or lalign_id_pair; type must be auto, dna or protein; lalignTopAlignments must be at least 1; proba_pair of RNA (a U) is not available |

## Examples

### Example 1: Four DNA sequences, `-method slow_pair,lalign_id_pair`

**Tool Call:**
```json
{
  "tool": "multiple_align_tcoffee",
  "arguments": {
    "sequences": [
      "ACGTTGCATGCCGATAGC",
      "ACGTGCATGCGATAGGC",
      "ACTTGCATCCGATTGC",
      "CGTTGCATGGCCGATAGC"
    ],
    "methods": [
      "slow_pair",
      "lalign_id_pair"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGTTGCAT-GCCGATAG-C",
    "ACGT-GCAT-GC-GATAGGC",
    "AC-TTGCAT--CCGATTG-C",
    "-CGTTGCATGGCCGATAG-C"
  ],
  "consensus": "ACGTTGCAT-GCCGATAG-C",
  "totalScore": 238
}
```

### Example 2: Five related protein sequences, the program's default (`t_coffee seqs.fa` = proba_pair)

**Tool Call:**
```json
{
  "tool": "multiple_align_tcoffee",
  "arguments": {
    "sequences": [
      "MKVLAAGIVGLLLAQPWTTQAEAKP",
      "MKVLAGIVGLLAQPWTQAEAKPLL",
      "MRVLAAGIIGLLLQPWTTQAEAK",
      "KVLAAGIVALLLAQPWTTAEAKPS",
      "MKVAAGIVGLLLAQPWTTQEAKP"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "MKVLAAGIVGLLLAQPWTTQAEAK-P-",
    "MK-VLAGIVGLL-AQPWTQ-AEAKPLL",
    "MRVLAAGIIGLL-LQPWTTQAEA--K-",
    "KV-LAAGIVALLLAQPWTT-AEAKPS-",
    "MK-VAAGIVGLLLAQPWTT-QEAK-P-"
  ],
  "consensus": "MK-LAAGIVGLLLAQPWTT-AEAK-P-",
  "totalScore": 559
}
```

## Verification status

Rows equal the installed `t_coffee` 13.41 (Debian t-coffee 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1, `-n_core 1`) on 1 480 committed DNA and protein sets of 2-8 sequences (the plain default, `slow_pair,lalign_id_pair`, the single methods, the reverse order and a three-method merge: guide tree and every row), on two sets of 100 / 101 sequences (fast neighbour joining) and the primary libraries pair by pair, weight by weight on 820 further sets; `scripts/alignment/generate_tcoffee_program_oracle.py`. Not available: `proba_pair` of RNA, the other `-method`s of the program (open items A1-13c / A1-13d).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.TCoffeeLibrary.cs#L352](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.TCoffeeLibrary.cs#L352)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
- Notredame, C., Higgins, D. G., Heringa, J. (2000). T-Coffee: A novel method for fast and accurate multiple sequence alignment. J Mol Biol 302(1): 205-217.
