# multiple_align_mafft

MAFFT 7.505's alignment of two sequences, ported from the program (stage 1 of the MAFFT port).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_mafft` |
| **Method ID** | `SequenceAligner.MultipleAlignMafft` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MAFFT 7.505's alignment of TWO sequences (Katoh & Standley 2013, MBE 30:772; mafft --anysymbol in.fa), ported step by step from the program's source (stage 1 of the MAFFT port; three or more sequences need the 6-mer distances, the guide tree and the group-group dynamic program and are refused with the pointer B09.md A2-19b): the alphabet is guessed from the letters (nucleotide when more than 75 % of the letters are A, C, G, T, U or N; sequenceType 'nuc' / 'amino' forces it), with anySymbol (default, as T-Coffee's mafft_pair calls it) letters outside the alphabet - digits, '*', '#', B, Z, J, O, U in proteins - are aligned as X (n for nucleotides) and restored in the output, '-' characters of the input are removed; scoring = the BLOSUM matrix (bl 30, 45, 50, 62 = default, 80, 90, 100 or 0) rescaled to 600 / mean diagonal with the offset (ep) subtracted, for nucleotides the 200-PAM Kimura matrix (kappa 2) with IUPAC codes by averaging, the gap opening penalty op (default 1.53, scaled by 0.6, or 1.8 for nucleotides, to an integer), no gap extension penalty and free end gaps (the script passes -O); dynamic program = G__align11 (a gap-open-only affine recurrence over 'best earlier cell' trackers, tie rules >= for the trackers, > for the cells) with the end of the alignment chosen in the last row / column. retree / fft / reorder / thread are accepted (no effect for two sequences). The default FFT-NS-2 strategy gives the same rows as this dynamic program on every committed pair, up to 1 500 residues. VERIFICATION: rows equal the installed Debian mafft 7.505-1 on 884 / 884 committed pairs (protein, DNA / RNA, IUPAC and unusual symbols, op / ep / bl / nuc / amino, the plain command without --anysymbol where nucleotides come back in lower case and amino acids in upper case). T-Coffee's mafft_pair uses this port.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.Mafft.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Mafft.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section MAFFT 7.505, stage 1)
- Oracle: `scripts/alignment/generate_mafft_pair_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/mafft_pair_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Exactly two sequences (letters; with anySymbol any character except blanks and '-' is a residue). |
| `sequenceType` | string | No | 'auto' (default: guessed from the letters), 'nuc' (mafft --nuc) or 'amino' (mafft --amino). |
| `gapOpen` | number | No | mafft --op: gap opening penalty (default 1.53). |
| `offset` | number | No | mafft --ep: offset (default 0.0). |
| `blosum` | integer | No | mafft --bl: BLOSUM matrix: 30, 45, 50, 62 (default), 80, 90, 100 or 0. |
| `anySymbol` | boolean | No | true (default): mafft --anysymbol, the original letters are restored. false: the plain command (letters and '.' only; nucleotides lower case, amino acids upper case). |
| `retree` | integer | No | mafft --retree (accepted, no effect for two sequences). |
| `fft` | boolean | No | mafft --fft / --nofft (accepted, no effect on the committed pairs). |
| `threads` | integer | No | mafft --thread (accepted, no effect). |
| `reorder` | boolean | No | mafft --reorder (accepted, no effect for two sequences). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The two rows (original letters with --anysymbol, gaps are -). |
| `consensus` | string | Consensus string of the alignment. |
| `totalScore` | integer | Natural affine sum of pairs (BLOSUM62 10 / 0.2 or EDNAFULL 10 / 0.5), rounded; not the program's own score. |
| `exactTotalScore` | number | The same score unrounded. |

## Errors

| Code | Message |
|------|---------|
| 1001 | sequences is null, a sequence has no residues, only one sequence, a character above U+007F, or one of = < > with anySymbol |
| 1002 | sequenceType is not auto / nuc / amino, blosum is not one of 30 45 50 62 80 90 100 0, gapOpen / offset not finite |
| 1003 | NotSupportedException: three or more sequences (the progressive stage is registered as B09.md A2-19b) |

## Examples

### Example 1: A protein pair with a lower-case sequence (--anysymbol)

**Tool Call:**
```json
{
  "tool": "multiple_align_mafft",
  "arguments": {
    "sequences": [
      "MKVLAAgIVG-LLLAQ",
      "mkvlagivallaq"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "MKVLAAgIVGLLLAQ",
    "mkvl-agiva-llaq"
  ],
  "consensus": "MKVLAAgIVGLLLAQ",
  "totalScore": 33,
  "exactTotalScore": 33
}
```

`mafft --quiet --anysymbol --thread 1 in.fa` writes `MKVLAAgIVGLLLAQ` / `mkvl-agiva-llaq` (the gap of the input is removed).

### Example 2: Nucleotides with IUPAC codes, plain command

**Tool Call:**
```json
{
  "tool": "multiple_align_mafft",
  "arguments": {
    "sequences": [
      "ACGTACGTTAGCNNRY",
      "ACGTACGATAGCCC"
    ],
    "anySymbol": false
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "acgtacgttagcnnry",
    "acgtacgatagccc--"
  ],
  "consensus": "acgtacgatagcccry",
  "totalScore": 37,
  "exactTotalScore": 36.5
}
```

`mafft --quiet --thread 1 in.fa` writes the rows in lower case, the shorter sequence ends with a free end gap.

### Example 3: BLOSUM80 with a user gap opening penalty and offset

**Tool Call:**
```json
{
  "tool": "multiple_align_mafft",
  "arguments": {
    "sequences": [
      "MKVLAAGIVGLLLAQWTR",
      "MKVLAGIVALLAQWTR"
    ],
    "blosum": 80,
    "gapOpen": 2.0,
    "offset": 0.1
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "MKVLAAGIVGLLLAQWTR",
    "MKVL-AGIVA-LLAQWTR"
  ],
  "consensus": "MKVLAAGIVALLLAQWTR",
  "totalScore": 54,
  "exactTotalScore": 54
}
```

`mafft --anysymbol --bl 80 --op 2.0 --ep 0.1 in.fa`.

## See Also

- [multiple_align_kalign](multiple_align_kalign.md)
- [multiple_align_muscle](multiple_align_muscle.md)
- [multiple_align_tcoffee](multiple_align_tcoffee.md)
- [multiple_align_clustalw](multiple_align_clustalw.md)
