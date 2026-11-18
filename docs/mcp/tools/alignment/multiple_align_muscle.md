# multiple_align_muscle

MUSCLE 3.8.1551's multiple alignment (Edgar 2004, NAR 32:1792), ported step by step from the program's source for muscle -maxiters 1 (the progressive draft) and -maxiters 2 (plus the tree-dependent refinement): letters upper-cased, gaps / white space dropped, letters outside the alphabet replaced by the wildcard (X protein, N nucleotide), a leading M of a protein aligned as X and restored (methionine hack); k-mer distance (kmer6_6 over the six MAFFT residue groups for proteins, kmer4_6 for nucleotides), UPGMB guide tree (0.1 x average + 0.9 x minimum linkage), ClustalW tree weights, profile-profile dynamic programming with position-specific gap open / close penalties and half-open terminal gaps, the log-expectation score (VTML200 likelihood ratios, center -0.52, gap open -2.9) for proteins or MUSCLE's nucleotide sum-of-pairs score (BLASTZ matrix, gap open -400) for DNA / RNA (scoreModel sp = PAM200, sv = VTML240 for proteins); maxIters 2 rebuilds the tree from the Kimura distances of the draft and realigns only the subtrees whose leaf sets changed. Rows are returned in input order (treeOrder = the order muscle writes). totalScore is the natural affine sum of pairs (BLOSUM62 10 / 0.2 for proteins, EDNAFULL 10 / 0.5 for nucleotides), not MUSCLE's objective. NOT ported: iterations 3 and above (restricted-partition refinement; maxIters > 2 is an error), -diags / -anchors, -weight / -distance2 / -root / -usetree / -matrix variants. VERIFICATION STATUS: the installed muscle is 5.x (a different program); the row-for-row oracle is the MUSCLE 3.8.1551 binary of the Docker hand-off package docs/Validation/review-2026-09/handoff/A1-10c-muscle and is awaiting its run - until it is integrated the alignments are NOT verified against muscle 3.8 (the score formulas, weights, tree, k-mer distance, Kimura table and the DP optimality are tested against the published definitions and the source).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_muscle` |
| **Method ID** | `SequenceAligner.MultipleAlignMuscle` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUSCLE 3.8.1551's multiple alignment (Edgar 2004, NAR 32:1792), ported step by step from the program's source for muscle -maxiters 1 (the progressive draft) and -maxiters 2 (plus the tree-dependent refinement): letters upper-cased, gaps / white space dropped, letters outside the alphabet replaced by the wildcard (X protein, N nucleotide), a leading M of a protein aligned as X and restored (methionine hack); k-mer distance (kmer6_6 over the six MAFFT residue groups for proteins, kmer4_6 for nucleotides), UPGMB guide tree (0.1 x average + 0.9 x minimum linkage), ClustalW tree weights, profile-profile dynamic programming with position-specific gap open / close penalties and half-open terminal gaps, the log-expectation score (VTML200 likelihood ratios, center -0.52, gap open -2.9) for proteins or MUSCLE's nucleotide sum-of-pairs score (BLASTZ matrix, gap open -400) for DNA / RNA (scoreModel sp = PAM200, sv = VTML240 for proteins); maxIters 2 rebuilds the tree from the Kimura distances of the draft and realigns only the subtrees whose leaf sets changed. Rows are returned in input order (treeOrder = the order muscle writes). totalScore is the natural affine sum of pairs (BLOSUM62 10 / 0.2 for proteins, EDNAFULL 10 / 0.5 for nucleotides), not MUSCLE's objective. NOT ported: iterations 3 and above (restricted-partition refinement; maxIters > 2 is an error), -diags / -anchors, -weight / -distance2 / -root / -usetree / -matrix variants. VERIFICATION STATUS: the installed muscle is 5.x (a different program); the row-for-row oracle is the MUSCLE 3.8.1551 binary of the Docker hand-off package docs/Validation/review-2026-09/handoff/A1-10c-muscle and is awaiting its run - until it is integrated the alignments are NOT verified against muscle 3.8 (the score formulas, weights, tree, k-mer distance, Kimura table and the DP optimality are tested against the published definitions and the source).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.Muscle.cs#L67](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Muscle.cs#L67)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section 4.9)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences to align (at least one, none empty): letters, '-', '.' and white space only. |
| `type` | string | No | 'auto' (default; MUSCLE's guess: nucleotide when >= 95 % of the first 100 residues are A, C, G, T, U or N), 'protein', 'dna' or 'rna' (MUSCLE -seqtype). (default: auto) |
| `maxIters` | integer | No | MUSCLE -maxiters: 1 = progressive draft, 2 (default) = plus tree refinement; larger values are not implemented. (default: 2) |
| `scoreModel` | string | No | Profile-profile score (-le / -sp / -sv / -spn): 'default' (LE for proteins, SPN for nucleotides), 'le', 'sp' (PAM200), 'sv' (VTML240), 'spn'. (default: default) |
| `gapOpen` | number | No | MUSCLE -gapopen (a negative score; default -2.9 LE, -1439 SP, -300 SV, -400 SPN DNA, -420 RNA). |
| `gapExtend` | number | No | MUSCLE -gapextend (default 0). |
| `center` | number | No | MUSCLE -center (default -0.52 LE, 0 otherwise, -300 RNA). |
| `distance1` | string | No | MUSCLE -distance1: 'default' (kmer6_6 protein, kmer4_6 nucleotide), 'kmer6_6' or 'kmer4_6' (nucleotides only). (default: default) |
| `cluster1` | string | No | MUSCLE -cluster1 (first guide tree): 'upgmb' (default), 'upgma', 'upgmamin' or 'upgmamax'. (default: upgmb) |
| `cluster2` | string | No | MUSCLE -cluster2 (tree of iteration 2): 'upgmb' (default), 'upgma', 'upgmamin' or 'upgmamax'. (default: upgmb) |
| `sueff` | number | No | MUSCLE -sueff (default 0.1): weight of the average linkage in UPGMB. (default: 0.1) |
| `maxTrees` | integer | No | MUSCLE -maxtrees (default 1): 0 switches the tree refinement of iteration 2 off. (default: 1) |
| `treeOrder` | boolean | No | true: rows in the order muscle writes them (depth-first leaf order of the final guide tree, -group) instead of the input order. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA (input order; the order MUSCLE writes with treeOrder). |
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a residue beats a gap on ties). |
| `totalScore` | integer | Natural affine sum-of-pairs score of the alignment (BLOSUM62 10 / 0.2 for proteins, EDNAFULL 10 / 0.5 for nucleotides; MUSCLE itself reports no such score). |

## Errors

| Code | Message |
|------|---------|
| 1010 | sequence is empty or contains a character other than a letter, '-', '.' or white space |
| 1002 | type must be auto, protein, dna or rna; maxIters must be 1 or 2 (iterations >= 3 are not ported); scoreModel / distance1 / cluster1 / cluster2 must be a known name; scoreModel and distance1 must suit the alphabet |

## Examples

### Example 1: Five related protein sequences (maxIters 2)

**Tool Call:**
```json
{
  "tool": "multiple_align_muscle",
  "arguments": {
    "sequences": [
      "MKVLAAGIVGLLLAQPWTTQAEAKP",
      "MKVLAGIVGLLAQPWTQAEAKPLL",
      "MRVLAAGIIGLLLQPWTTQAEAK",
      "KVLAAGIVALLLAQPWTTAEAKPS",
      "MKVAAGIVGLLLAQPWTTQEAKP"
    ],
    "type": "protein"
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "MKVLAAGIVGLLLAQPWTTQAEAKP--",
    "MKV-LAGIVG-LLAQPWTQA-EAKPLL",
    "MRVLAAGIIG-LLLQPWTTQAEAK---",
    "-KVLAAGIVALLLAQPWTTA-EAKPS-",
    "MKV-AAGIVGLLLAQPWTTQ-EAKP--"
  ],
  "consensus": "MKVLAAGIVGLLLAQPWTTQ-EAKP--",
  "totalScore": 656
}
```

### Example 2: Four DNA sequences (progressive draft, maxIters 1)

**Tool Call:**
```json
{
  "tool": "multiple_align_muscle",
  "arguments": {
    "sequences": [
      "ACGTACGTTAGCTAGCTAGGATCGATCGGATC",
      "ACGTACGTAGCTAGCTAGGATCGATCGATC",
      "ACGTTCGTTAGCTAGCTAGGATCGGATCGGATC",
      "ACGTACGTTAGCTTAGCTAGGATCGATCGGATC"
    ],
    "type": "dna",
    "maxIters": 1
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGTACGTTAGC-TAGCTAGGATC-GATCGGATC",
    "ACGTACG-TAGC-TAGCTAGGATC-GATC-GATC",
    "ACGTTCGTTAGC-TAGCTAGGATCGGATCGGATC",
    "ACGTACGTTAGCTTAGCTAGGATC-GATCGGATC"
  ],
  "consensus": "ACGTACGTTAGC-TAGCTAGGATC-GATCGGATC",
  "totalScore": 783
}
```

## Verification status

The installed `muscle` of the audit machine is 5.x (the PPP / probcons-like pipeline), a different program, so it cannot be the oracle. The oracle is the MUSCLE 3.8.1551 binary built from the same Debian source by the Docker hand-off package `docs/Validation/review-2026-09/handoff/A1-10c-muscle` (420 sets): **awaiting its run**. Until it is integrated the alignments are not verified against muscle 3.8; the score formulas (log-expectation / SPN, VTML200 / BLASTZ tables), the ClustalW weights example of the source, the UPGMA test matrix of the source, the Kimura / Dayhoff distance, the k-mer distance definition, the wildcard rules and the optimality of the profile dynamic program (exhaustive enumeration, 480 random profile pairs) are tested.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.Muscle.cs#L67](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Muscle.cs#L67)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
- Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5): 1792-1797.
