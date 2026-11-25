# multiple_align_poa

POA 2.0 partial-order alignment (port of the program)

## Overview

| Property | Value |
|----------|-------|
| **Server** | MultipleAlignment |
| **Tool Name** | `multiple_align_poa` |
| **Method ID** | `SequenceAligner.PoaRun` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

POA 2.0's partial-order multiple alignment (Lee, Grasso & Sharlow 2002, Bioinformatics 18:452; poa -read_fasta in.fa -clustal out.aln matrix), ported step by step from the program's source: every sequence is aligned to the growing partial-order graph by dynamic programming over the graph (each node scored against each residue, predecessors of both graphs, the program's affine gap model with a truncation / decay length and a free initial gap in local mode, tie-breaking by the smallest end cell) and fused into it (aligned identical letters share a node, other aligned letters form an alignment ring, the ring is kept contiguous in the graph's linear order). Default is iterative local alignment in input order; `global` aligns end to end, `progressive` merges the pairs in descending order of their alignment scores (or of the `pairScores` text of name1 name2 score lines), `fuseAll` fuses identical letters on rings, `heaviestBundle` adds the heaviest-bundle consensus sequences CONSENS0, CONSENS1, ... (a sequence joins a bundle when at least `heaviestBundleMinimum` of its residues lie on the path), `preserveSequenceOrder` restores the input order after progressive merging. The scoring matrix is a POA matrix file as text (symbol line, integer rows, GAP-PENALTIES=open ext1 ext2, GAP-PENALTIES-X, GAP-TRUNCATION-LENGTH, GAP-DECAY-LENGTH); omitted = BLOSUM80 with gap penalties 12 / 6 / 6 (the blosum80.mat of the distribution). Letters that are not symbols of the matrix become its first symbol (so lower-case letters need `case` = 'upper' for the protein matrix). Returns the rows, the program's PO file (the graph), the FASTA-PIR text and the CLUSTAL text. VERIFICATION: rows and all three files equal the installed Debian poa 2.0+20060928 on 359 / 359 committed command lines (protein BLOSUM80, nucleotide matrices with truncation / decay gap models, every option; scripts/alignment/generate_poa_oracle.py), see poa_command. Not ported: MSA / PO file input (-read_msa, -read_msa2, -read_msa_list, -subset, -remove).

## Core Documentation Reference

- Source: [Seqeron.Genomics.PartialOrderAlignment/PoaAligner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.PartialOrderAlignment/PoaAligner.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section POA 2.0)
- Oracle: `scripts/alignment/generate_poa_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/poa_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences to align (at least one non-empty): residues; white space is removed. |
| `names` | array<string> | No | Sequence names (the first word of a FASTA header, at most 31 characters are kept); omit for S1, S2, ... |
| `matrix` | string | No | POA matrix file as text; omit for the BLOSUM80 matrix of the distribution (gap penalties 12 / 6 / 6). |
| `global` | boolean | No | poa -do_global: global alignment (default: local alignment with free end gaps). |
| `progressive` | boolean | No | poa -do_progressive: merge sequences / partial orders in descending order of their pairwise alignment scores (default: iterative in input order). |
| `pairScores` | string | No | poa -read_pairscores: the text of a pair-score file (lines 'name1 name2 score'); the merge order follows these scores. |
| `fuseAll` | boolean | No | poa -fuse_all: fuse identical letters on align rings. (default: false) |
| `heaviestBundle` | boolean | No | poa -hb: heaviest bundling: add the consensus sequences CONSENS0, CONSENS1, ... to the alignment. (default: false) |
| `heaviestBundleMinimum` | number | No | poa -hbmin: minimum fraction of a sequence's residues on a bundle's path for the sequence to join it (default 0.9). |
| `best` | boolean | No | poa -best: restrict the PIR text to the sequences of the heaviest bundle (bundle 0; every alignment column is still written). (default: false) |
| `preserveSequenceOrder` | boolean | No | poa -preserve_seqorder: write the sequences in their input order (progressive merging reorders them). (default: false) |
| `case` | string | No | Letter case of the input: 'keep' (default), 'lower' (poa -tolower) or 'upper' (poa -toupper). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `names` | array<string> | Name of every row: the sequences in the program's output order, then the CONSENS sequences of heaviestBundle. |
| `alignedSequences` | array<string> | The alignment, one row per name (gap -; letters are the symbols of the matrix). |
| `order` | array<integer> | Input index of every row (-1 for a consensus sequence). |
| `po` | string | The partial-order alignment in the program's PO file format (poa -po). |
| `pir` | string | The alignment in the program's FASTA-PIR format (poa -pir; gap ., 60 columns; restricted to bundle 0 with best). |
| `clustal` | string | The alignment in the program's CLUSTAL format (poa -clustal; names indented to column 36, 50 columns per block). |
| `stdout` | string | What the program writes to standard output: the permutation line of preserveSequenceOrder. |

## Errors

| Code | Message |
|------|---------|
| 1001 | sequences is null or empty, or a sequence has no residues |
| 1002 | matrix is not a POA matrix file (a row symbol that is not in the symbol line, a missing score), or case is not keep / lower / upper |

## Examples

### Example 1: Four proteins with the heaviest-bundle consensus (poa -hb, default BLOSUM80)

**Tool Call:**
```json
{
  "tool": "multiple_align_poa",
  "arguments": {
    "sequences": [
      "MKVLAAGIVGLLLAQPSTAEEKV",
      "MKVLAGIVGLLAQPSTAEEKVLL",
      "MKILAAGIAGLLLQPSAEEKV",
      "MRVLAAGLVGLLLAQSTAEEK"
    ],
    "names": [
      "a",
      "b",
      "c",
      "d"
    ],
    "heaviestBundle": true
  }
}
```

**Response:**
```json
{
  "names": [
    "a",
    "b",
    "c",
    "d",
    "CONSENS0",
    "CONSENS1"
  ],
  "alignedSequences": [
    "MKVLAAGIVGLLLAQPSTAEEKV--",
    "MKVLA-GIVGLL-AQPSTAEEKVLL",
    "MKILAAGIAGLL-LQPS-AEEKV--",
    "MRVLAAGLVGLLLAQ-STAEEK---",
    "MKVLAAGIVGLLLAQPSTAEEKVLL",
    "MKILAAGIAGLL-LQPS-AEEKV--"
  ],
  "order": [
    0,
    1,
    2,
    3,
    -1,
    -1
  ],
  "po": "VERSION=LPO.0.1\nNAME=a\nTITLE=untitled\nLENGTH=30\nSOURCECOUNT=6\nSOURCENAME=a\nSOURCEINFO=23 0 0 0 untitled\nSOURCENAME=b\nSOURCEINFO=23 0 0 0 untitled\nSOURCENAME=c\nSOURCEINFO=21 0 0 1 untitled\nSOURCENAME=d\nSOURCEINFO=21 0 0 0 untitled\nSOURCENAME=CONSENS0\nSOURCEINFO=25 0 0 0 consensus produced by heaviest_bundle, containing 3 seqs\nSOURCENAME=CONSENS1\nSOURCEINFO=21 0 0 1 consensus produced by heaviest_bundle, containing 1 seqs\nM:S0S1S2S3S4S5\nR:L0S3A2\nK:L0S0S1S2S4S5A1\nI:L2S2S5A4\nV:L2L1S0S1S3S4A3\nL:L4L3S0S1S2S3S4S5\nA:L5S0S1S2S3S4S5\nA:L6S0S2S3S4S5\nG:L7L6S0S1S2S3S4S5\nL:L8S3A10\nI:L8S0S1S2S4S5A9\nA:L10S2S5A12\nV:L10L9S0S1S3S4A11\nG:L12L11S0S1S2S3S4S5\nL:L13S0S1S2S3S4S5\nL:L14S0S1S2S3S4S5\nL:L15S0S3S4\nL:L15S2S5A18\nA:L16L15S0S1S3S4A17\nQ:L18L17S0S1S2S3S4S5\nP:L19S0S1S2S4S5\nS:L20L19S0S1S2S3S4S5\nT:L21S0S1S3S4\nA:L22L21S0S1S2S3S4S5\nE:L23S0S1S2S3S4S5\nE:L24S0S1S2S3S4S5\nK:L25S0S1S2S3S4S5\nV:L26S0S1S2S4S5\nL:L27S1S4\nL:L28S1S4\n",
  "pir": ">a untitled\nMKVLAAGIVGLLLAQPSTAEEKV..\n>b untitled\nMKVLA.GIVGLL.AQPSTAEEKVLL\n>c untitled\nMKILAAGIAGLL.LQPS.AEEKV..\n>d untitled\nMRVLAAGLVGLLLAQ.STAEEK...\n>CONSENS0 consensus produced by heaviest_bundle, containing 3 seqs\nMKVLAAGIVGLLLAQPSTAEEKVLL\n>CONSENS1 consensus produced by heaviest_bundle, containing 1 seqs\nMKILAAGIAGLL.LQPS.AEEKV..\n",
  "clustal": "CLUSTAL W (1.74) multiple sequence alignment\n\n\na                                   MKVLAAGIVGLLLAQPSTAEEKV--\nb                                   MKVLA-GIVGLL-AQPSTAEEKVLL\nc                                   MKILAAGIAGLL-LQPS-AEEKV--\nd                                   MRVLAAGLVGLLLAQ-STAEEK---\nCONSENS0                            MKVLAAGIVGLLLAQPSTAEEKVLL\nCONSENS1                            MKILAAGIAGLL-LQPS-AEEKV--\n\n\n",
  "stdout": ""
}
```

`poa -read_fasta in.fa -hb -clustal c.aln -pir p.pir -po p.po blosum80.mat` writes exactly these rows and files.

### Example 2: Global progressive alignment with the original sequence order restored

**Tool Call:**
```json
{
  "tool": "multiple_align_poa",
  "arguments": {
    "sequences": [
      "MKVLAAGIVGLLLAQPSTAEEKV",
      "MKVLAGIVGLLAQPSTAEEKVLL",
      "MKILAAGIAGLLLQPSAEEKV",
      "MRVLAAGLVGLLLAQSTAEEK"
    ],
    "names": [
      "a",
      "b",
      "c",
      "d"
    ],
    "global": true,
    "progressive": true,
    "preserveSequenceOrder": true
  }
}
```

**Response:**
```json
{
  "names": [
    "a",
    "b",
    "c",
    "d"
  ],
  "alignedSequences": [
    "MKVLAAGIVGLLLAQPSTAEEKV--",
    "MKVLA-GIVGLLA-QPSTAEEKVLL",
    "MKILAAGIAGLLL-QPS-AEEKV--",
    "MRVLAAGLVGLLLAQ-STAEEK---"
  ],
  "order": [
    0,
    1,
    2,
    3
  ],
  "po": "VERSION=LPO.0.1\nNAME=a\nTITLE=untitled\nLENGTH=30\nSOURCECOUNT=4\nSOURCENAME=a\nSOURCEINFO=23 0 1 -1 untitled\nSOURCENAME=b\nSOURCEINFO=23 0 1 -1 untitled\nSOURCENAME=c\nSOURCEINFO=21 0 1 -1 untitled\nSOURCENAME=d\nSOURCEINFO=21 0 1 -1 untitled\nM:S0S2S3S1\nR:L0S3A2\nK:L0S0S2S1A1\nI:L2S2A4\nV:L2L1S0S3S1A3\nL:L4L3S0S2S3S1\nA:L5S0S2S3S1\nA:L6S0S2S3\nG:L7L6S0S2S3S1\nL:L8S3A10\nI:L8S0S2S1A9\nA:L10S2A12\nV:L10L9S0S3S1A11\nG:L12L11S0S2S3S1\nL:L13S0S2S3S1\nL:L14S0S2S3S1\nA:L15S1A17\nL:L15S0S2S3A16\nA:L17S0S3\nQ:L18L17L16S0S2S3S1\nP:L19S0S2S1\nS:L20L19S0S2S3S1\nT:L21S0S3S1\nA:L22L21S0S2S3S1\nE:L23S0S2S3S1\nE:L24S0S2S3S1\nK:L25S0S2S3S1\nV:L26S0S2S1\nL:L27S1\nL:L28S1\n",
  "pir": ">a untitled\nMKVLAAGIVGLLLAQPSTAEEKV..\n>b untitled\nMKVLA.GIVGLLA.QPSTAEEKVLL\n>c untitled\nMKILAAGIAGLLL.QPS.AEEKV..\n>d untitled\nMRVLAAGLVGLLLAQ.STAEEK...\n",
  "clustal": "CLUSTAL W (1.74) multiple sequence alignment\n\n\na                                   MKVLAAGIVGLLLAQPSTAEEKV--\nb                                   MKVLA-GIVGLLA-QPSTAEEKVLL\nc                                   MKILAAGIAGLLL-QPS-AEEKV--\nd                                   MRVLAAGLVGLLLAQ-STAEEK---\n\n\n",
  "stdout": "0 2 3 1 \n"
}
```

The program prints the permutation `0 2 3 1` (merge order) to standard output.

## See Also

- [poa_command](poa_command.md)
- [multiple_align_probcons](multiple_align_probcons.md)
- [multiple_align_kalign](multiple_align_kalign.md)
