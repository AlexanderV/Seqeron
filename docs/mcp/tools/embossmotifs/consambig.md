# consambig

EMBOSS consambig: ambiguity-code consensus of an alignment.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `consambig` |
| **Method ID** | `MotifFinder.GenerateEmbossAmbiguousConsensus` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Ambiguity-code consensus of an alignment, identical to EMBOSS 6.6.0 'consambig' (1 400 random alignments equal): per column the residue codes are OR-ed (nucleotides A, C, G, T as bits, so R = A|G, V = A|C|G, N = any; amino acids one bit each, B = D|N, Z = E|Q, J = I|L, X = all twenty) and written as the one letter that stands for exactly that set; there is no vote (unlike cons: one deviating base widens the code). A column that contains a gap ('-', '.', '~', or the padding of a shorter row) is written in lower case. An all-gap nucleotide column is 'n'; an amino-acid union that is no letter (A|D) is 'X'; a protein column without any scored residue (all gaps, or only '*') is the NUL character, as in the program. The set type is that of the first sequence (A/C/G/T/U + IUPAC letters = nucleotide) unless residueType 'nucleotide' forces it (residueType 'protein' cannot make a nucleotide-looking set protein, as in the program). Rows are upper-cased, '.'/'~' read as '-', '?' as X; shorter rows are padded with trailing gaps (the program always does). Rows are read as EMBOSS reads an alignment (A4-3): letters and '*' are kept, '.' '~' '-' are gaps, '?' '#' '+' are the unknown residue, every other character (digits, blanks, '_', punctuation) is dropped, and rows that become unequal are padded with trailing gaps.

## Core Documentation Reference

- Source: [MotifFinder.AmbiguousConsensus.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.Motifs/EmbossMotifPrograms.AmbiguousConsensus.cs)
- Algorithm doc: [Consensus_From_Alignment.md](../../../algorithms/Pattern_Matching/Consensus_From_Alignment.md)
- Oracle: [generate_consambig_oracle.py](../../../../scripts/analysis/generate_consambig_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | At least two aligned sequences (DNA or protein). |
| `residueType` | string | No | 'auto' (default: type of the first sequence, as the program), 'nucleotide' (forces nucleotide, -snucleotide) or 'protein'. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `consensus` | string | One symbol per alignment column: the IUPAC / amino-acid letter standing for the union of the column; lower case when the column has a gap. |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least two sequences are required |
| 1001 | Invalid alignment character |

## Examples

### Example 1: Nucleotide alignment with gaps and an ambiguity code

**User Prompt:**
> Make the ambiguous consensus of AC-GTN, AC-GTA and ARGGT-.

**Expected Tool Call:**
```json
{"tool": "consambig", "arguments": {"sequences": ["AC-GTN", "AC-GTA", "ARGGT-"]}}
```

**Response:**
```json
{"consensus": "AVgGTn"}
```

consambig prints AVgGTn: V = A|C|G in column 2, the gap columns are lower case, N = any base.

### Example 2: Protein alignment

**User Prompt:**
> Ambiguous consensus of the protein rows DI and NL.

**Expected Tool Call:**
```json
{"tool": "consambig", "arguments": {"sequences": ["DI", "NL"], "residueType": "protein"}}
```

**Response:**
```json
{"consensus": "BJ"}
```

D|N = B (aspartate or asparagine), I|L = J (leucine or isoleucine).

## See Also

- [generate_emboss_consensus](generate_emboss_consensus.md) - the plurality consensus of EMBOSS cons
- [generate_consensus](generate_consensus.md) - majority IUPAC consensus
- [multiple_alignment_statistics](multiple_alignment_statistics.md) - alignment statistics
