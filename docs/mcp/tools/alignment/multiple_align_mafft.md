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

MAFFT 7.505 (Katoh & Standley 2013, MBE 30:772; mafft [--anysymbol] [options] in.fa), ported step by step from the program's source, for ANY number of sequences: the progressive strategies FFT-NS-2 (default), FFT-NS-1 (retree 1), NW-NS-2 (fft false) and NW-NS-1 (fft false, retree 1; retree 3 builds the guide tree three times). The pipeline is the program's (script mafft -> disttbfast): the alphabet is guessed from the letters (nucleotide when more than 75 % of the letters are A, C, G, T, U or N; sequenceType 'nuc' / 'amino' forces it); with anySymbol letters outside the alphabet - digits, '*', '#', B, Z, J, O, U in proteins - are aligned as X (n for nucleotides) and restored in the output (the plain command stops on them), '-' characters of the input are removed; the distance of two sequences is the number of shared 6-mers (10-mers with tupleSize 10, nucleotides only; amino acids in the six classes of the program) scaled by the lengths, the guide tree is the UPGMA of MAFFT (mixed linkage: treeMixing 0.1 of average and 0.9 of minimum linkage; 1 = --averagelinkage, 0 = --minimumlinkage) with the weights of the sequences from the branch lengths (unweight switches them off), the groups are aligned by the fast Fourier transform of Falign (the polarity and volume of the amino acids, the letter 'a' of nucleotides: 20 candidate shifts, windows with a high score as anchors, a chain of anchors, the segments between them aligned by the dynamic program) when fft is true, and by the group-group dynamic program A__align (gap penalties scaled by the gap frequencies of the columns of the two groups, history of the profiles of the children for groups of more than 20 sequences; legacyGapPenalty = --leavegappyregion switches the scaling off) or the dynamic program G__align11 for two single sequences; for retree 2 the distances are then recomputed from the first alignment, a second tree is built and the sequences are aligned again. Scoring = the BLOSUM matrix (bl 30, 45, 50, 62 = default, 80, 90, 100 or 0), the JTT or the transmembrane PAM matrix (scoringModel 'jtt' / 'tm' with pam steps, default 200) rescaled to 600 / mean diagonal with the offset (ep) subtracted, for nucleotides the Kimura matrix (kappa 2, 200 PAM; pam sets the distance; frequencyModel uses the letter frequencies of the data; wildcardN scores n like a wildcard) with IUPAC codes by averaging; gap opening penalty op (default 1.53, scaled by 0.6, or 1.8 for nucleotides, to an integer; opdist = --opdist for the distance step), gap extension exp (default 0), free end gaps unless terminalGapPenalty. reorder returns the sequences in the order of the guide tree (--reorder). thread is accepted: the program's alignment does not depend on it; its progress text with threads is thread-order dependent and is written as of one thread. The result carries the program's standard error (the banner lines of disttbfast and the strategy text of the script) unless quiet. Not available (refused): --retree 0, the memory-saving mode for groups longer than 30 000 columns; the iterative (maxiterate), pair-library (localpair, globalpair, genafpair), parttree, add / merge / seed strategies are registered stages of the port. VERIFICATION: rows and standard error equal the installed Debian mafft 7.505-1 on every committed run of scripts/alignment/generate_mafft_progressive_oracle.py (6 530 sets in 28 strategies, at least 210 per strategy: protein and DNA, 3 - 45 sequences, tiny and unrelated sequences, unusual symbols, every option above), and on 884 / 884 pairs (generate_mafft_pair_oracle.py). T-Coffee's mafft_pair uses this port.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.Mafft.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Mafft.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section MAFFT 7.505, stage 1)
- Oracle: `scripts/alignment/generate_mafft_pair_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/mafft_pair_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Two or more sequences (letters; with anySymbol any character except blanks and '-' is a residue). Like the program: an empty list gives an empty result and ONE sequence is written back (gaps dropped; without anySymbol nucleotides in lower and amino acids in upper case). |
| `sequenceType` | string | No | 'auto' (default: guessed from the letters), 'nuc' (mafft --nuc) or 'amino' (mafft --amino). |
| `gapOpen` | number | No | mafft --op: gap opening penalty (default 1.53). |
| `offset` | number | No | mafft --ep: offset (default 0.0; a negative value adds to every substitution score, a positive one subtracts). |
| `blosum` | integer | No | mafft --bl: BLOSUM matrix of the amino-acid scoring: 30, 45, 50, 62 (default), 80, 90, 100 or 0. |
| `anySymbol` | boolean | No | true (default): mafft --anysymbol, the original letters are restored in the output. false: the plain command (letters and '.' only; nucleotides lower case, amino acids upper case in the rows). |
| `retree` | integer | No | mafft --retree: how many times the guide tree is built (1 = FFT-NS-1, 2 = FFT-NS-2 default, 3; 0 is refused). Two sequences always use 1. (default: 2) |
| `fft` | boolean | No | mafft --fft (true, default) / --nofft (false, NW-NS). (default: true) |
| `threads` | integer | No | mafft --thread (default 0; accepted, no effect on the alignment). |
| `reorder` | boolean | No | mafft --reorder: the rows in the order of the guide tree. (default: false) |
| `gapExtension` | number | No | mafft --exp: gap extension penalty (default 0). |
| `gapOpenDistance` | number | No | mafft --opdist: gap opening penalty of the distance step (default = gapOpen). |
| `scoringModel` | string | No | 'blosum' (default; mafft --bl), 'jtt' (--jtt / --kimura) or 'tm' (--tm): the scoring matrix of amino acids. |
| `pam` | integer | No | PAM steps of --jtt / --tm / --kimura (default 200). |
| `kappa` | integer | No | mafft --kappa: transition / transversion ratio of the nucleotide matrix (default: the program's 2). |
| `frequencyModel` | boolean | No | mafft --fmodel: letter frequencies of the data in the scoring matrix. (default: false) |
| `wildcardN` | boolean | No | mafft --nwildcard: n scores like a wildcard in nucleotide sequences. (default: false) |
| `forceFft` | boolean | No | mafft --fft given explicitly: the FFT is used for every pair of groups. (default: false) |
| `tupleSize` | integer | No | mafft --6merpair (6, default) / --10merpair (10, nucleotides only). (default: 6) |
| `treeMixing` | number | No | mafft --mixedlinkage x (0 = --minimumlinkage, 1 = --averagelinkage; default 0.1). |
| `legacyGapPenalty` | boolean | No | mafft --leavegappyregion / --legacygappenalty. (default: false) |
| `terminalGapPenalty` | boolean | No | mafft --termgappenalty: penalise gaps at the ends. (default: false) |
| `unweight` | boolean | No | mafft --unweight: no weighting of the sequences by the guide tree. (default: false) |
| `quiet` | boolean | No | mafft --quiet: no standard error text. (default: false) |
| `format` | string | No | Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text. |
| `names` | array of string | No | Optional names of the rows for the formatted output (default seq1, seq2, ...). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The two rows (original letters with --anysymbol, gaps are -). |
| `consensus` | string | Consensus string of the alignment. |
| `totalScore` | integer | Natural affine sum of pairs (BLOSUM62 10 / 0.2 or EDNAFULL 10 / 0.5), rounded; not the program's own score. |
| `exactTotalScore` | number | The same score unrounded. |
| `formatted` | string | No | The alignment written in the requested format (null when format is omitted). |

## Errors

| Code | Message |
|------|---------|
| 1001 | sequences is null, a sequence has no residues, only one sequence, a character above U+007F, or one of = < > with anySymbol |
| 1002 | sequenceType is not auto / nuc / amino, blosum is not one of 30 45 50 62 80 90 100 0, gapOpen / offset not finite |
| 1003 | NotSupportedException: three or more sequences (the progressive stage is not implemented) |

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
