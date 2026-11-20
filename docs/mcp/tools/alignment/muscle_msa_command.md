# muscle_msa_command

The other command lines of MUSCLE 3.8.1551 that work on existing alignments, ported from the source: mode 'profile' (muscle -profile -in1 a -in2 b: profile-profile alignment of two alignments, rows of alignment followed by those of alignment2), 'refine' (-refine: maxIters rounds of restricted-partition refinement of an alignment, the program's default is 8), 'spscore' (-spscore: the sum-of-pairs objective score of an alignment, 'SP=%.4g'), 'ppscore' (-ppscore: the profile-profile score of two alignments of equal length; nucleotides only), 'maketree' (-maketree: the guide tree of an alignment as a Newick file text; the program never parses -distance2 in this mode, the distance is always pctidkimura), 'refinew' (-refinew: every window of refineWindow columns re-aligned from scratch with the full pipeline; nucleotides only), 'sw' (-sw: the local Smith-Waterman alignment of two profiles, residues outside the local alignment in lower case and padded with '.'; nucleotides only; needs distance1 pctidkimura when a profile has two or more sequences). For proteins the program never initialises the score matrix in ppscore / refinew / sw and crashes, so those modes are rejected for proteins; mode 'profdb' (-profdb) is always rejected: the program never initialises the alphabet there and crashes. The sequence weights are weight1 (default ClustalW weights of a guide tree built from the alignment with cluster2 / distance2 / root1 or root2). VERIFICATION STATUS: implemented from the program's source and verified against the MUSCLE 3.8.1551 binary on the second Docker round (committed oracle muscle3_r2_oracle.tsv: rows, trees, scores and refusals equal on every recorded set; protein ppscore / sw / refinew and profdb are refused because the program crashes).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `muscle_msa_command` |
| **Method ID** | `SequenceAligner.MuscleProfileAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The other command lines of MUSCLE 3.8.1551 that work on existing alignments, ported from the source: mode 'profile' (muscle -profile -in1 a -in2 b: profile-profile alignment of two alignments, rows of alignment followed by those of alignment2), 'refine' (-refine: maxIters rounds of restricted-partition refinement of an alignment, the program's default is 8), 'spscore' (-spscore: the sum-of-pairs objective score of an alignment, 'SP=%.4g'), 'ppscore' (-ppscore: the profile-profile score of two alignments of equal length; nucleotides only), 'maketree' (-maketree: the guide tree of an alignment as a Newick file text; the program never parses -distance2 in this mode, the distance is always pctidkimura), 'refinew' (-refinew: every window of refineWindow columns re-aligned from scratch with the full pipeline; nucleotides only), 'sw' (-sw: the local Smith-Waterman alignment of two profiles, residues outside the local alignment in lower case and padded with '.'; nucleotides only; needs distance1 pctidkimura when a profile has two or more sequences). For proteins the program never initialises the score matrix in ppscore / refinew / sw and crashes, so those modes are rejected for proteins; mode 'profdb' (-profdb) is always rejected: the program never initialises the alphabet there and crashes. The sequence weights are weight1 (default ClustalW weights of a guide tree built from the alignment with cluster2 / distance2 / root1 or root2). VERIFICATION STATUS: implemented from the program's source and verified against the MUSCLE 3.8.1551 binary on the second Docker round (committed oracle muscle3_r2_oracle.tsv: rows, trees, scores and refusals equal on every recorded set; protein ppscore / sw / refinew and profdb are refused because the program crashes).

## Core Documentation Reference

Modes and methods: `profile` = `SequenceAligner.MuscleProfileAlign`, `refine` = `MuscleRefine`, `spscore` = `MuscleScoreAlignment`, `ppscore` = `MuscleProfileScore`, `maketree` = `MuscleGuideTree`, `refinew` = `MuscleRefineWindows`, `sw` = `MuscleLocalAlign`, `profdb` = `MuscleProfileDatabase` (always rejected: the program crashes).

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.MuscleProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.MuscleProgram.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section 4.9)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `mode` | string | Yes | Mode: 'profile', 'refine', 'spscore', 'ppscore', 'maketree', 'refinew', 'sw' or 'profdb' (always rejected, the program crashes). |
| `alignment` | array<string> | Yes | The (first) alignment: equal-length rows, '-' = gap. |
| `alignment2` | array<string> | No | The second alignment (profile, ppscore, sw). |
| `type` | string | No | 'auto' (MUSCLE's guess from the alignment), 'protein', 'dna' or 'rna'. (default: auto) |
| `maxIters` | integer | No | MUSCLE -maxiters for refine / refinew (default 8, the program's). |
| `scoreModel` | string | No | 'default', 'le', 'sp', 'sv', 'spn' (profile-profile score). (default: default) |
| `gapOpen` | number | No | MUSCLE -gapopen. |
| `gapExtend` | number | No | MUSCLE -gapextend. |
| `center` | number | No | MUSCLE -center. |
| `distance1` | string | No | MUSCLE -distance1 (tree of the sw profiles; also the first tree of refinew). (default: default) |
| `distance2` | string | No | MUSCLE -distance2 (guide tree of an alignment): 'default' (pctidkimura), 'pctidlog', 'edit', 'scoredist'. (default: default) |
| `cluster1` | string | No | MUSCLE -cluster1. (default: upgmb) |
| `cluster2` | string | No | MUSCLE -cluster2 (guide tree of an alignment): 'upgmb', 'upgma', 'upgmamin', 'upgmamax', 'neighborjoining'. (default: upgmb) |
| `root1` | string | No | MUSCLE -root1 (profile / ppscore / sw trees). (default: pseudo) |
| `root2` | string | No | MUSCLE -root2 (refine / spscore / maketree tree). (default: pseudo) |
| `weight1` | string | No | MUSCLE -weight1 (sequence weights): 'clustalw' (default), 'none', 'henikoff', 'henikoffpb', 'gsc', 'threeway'. |
| `weight2` | string | No | MUSCLE -weight2. (default: clustalw) |
| `objScore` | string | No | MUSCLE -objscore (refinement objective). (default: spm) |
| `termGaps` | string | No | MUSCLE -termgaps. (default: half) |
| `anchors` | boolean | No | MUSCLE -anchors (true) / -noanchors (false). (default: true) |
| `refineWindow` | integer | No | MUSCLE -refinewindow (refinew, default 200). (default: 200) |
| `fromWindow` | integer | No | MUSCLE -fromwindow (refinew). (default: 0) |
| `toWindow` | integer | No | MUSCLE -towindow (refinew; 0 = the last window). (default: 0) |
| `windowOffset` | integer | No | MUSCLE -windowoffset (refinew). (default: 0) |
| `names` | array<string> | No | FASTA labels of the rows (default s0, s1, ...), the leaf names of the maketree tree. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | Rows of the result alignment (profile, refine, refinew, sw); empty for spscore / ppscore / maketree. |
| `score` | number | The score (spscore: sum-of-pairs objective; ppscore: profile-profile score), or null. |
| `tree` | string | The guide tree in the program's text format (maketree), or null. |

## Errors

| Code | Message |
|------|---------|
| 1002 | mode must be 'profile', 'refine', 'spscore', 'ppscore', 'maketree', 'refinew', 'sw' or 'profdb'; alignment2 is required for profile / ppscore / sw / profdb; option names must be known |
| 1010 | alignment must contain at least one row; rows of an alignment must have equal length |
| 1003 | the program refuses or crashes (ppscore of proteins: no protein score matrix in this mode; profdb: the alphabet is never initialised) - reported as an error |

## Examples

### Example 1: Sum-of-pairs score of a protein alignment (muscle -spscore)

**Tool Call:**
```json
{
  "tool": "muscle_msa_command",
  "arguments": {
    "mode": "spscore",
    "alignment": [
      "MKVLAAGIVGLLLAQPWTTQAEAKP--",
      "MKV-LAGIVG-LLAQPWTQA-EAKPLL",
      "MRVLAAGIIG-LLLQPWTTQAEAK---",
      "-KVLAAGIVALLLAQPWTTA-EAKPS-",
      "MKV-AAGIVGLLLAQPWTTQ-EAKP--"
    ],
    "type": "protein"
  }
}
```

**Response:**
```json
{
  "alignedSequences": [],
  "score": 33.679054260253906
}
```

### Example 2: Guide tree of an alignment (muscle -maketree)

**Tool Call:**
```json
{
  "tool": "muscle_msa_command",
  "arguments": {
    "mode": "maketree",
    "alignment": [
      "MKVLAAGIVGLLLAQPWTTQAEAKP--",
      "MKV-LAGIVG-LLAQPWTQA-EAKPLL",
      "MRVLAAGIIG-LLLQPWTTQAEAK---",
      "-KVLAAGIVALLLAQPWTTA-EAKPS-",
      "MKV-AAGIVGLLLAQPWTTQ-EAKP--"
    ],
    "type": "protein"
  }
}
```

**Response:**
```json
{
  "alignedSequences": [],
  "tree": "(\ns1:0.0820507\n,\n(\ns2:0.0757844\n,\n(\ns3:0.0464272\n,\n(\ns0:-0\n,\ns4:-0\n):0.0464272\n):0.0293572\n):0.00626627\n)\n;\n"
}
```

### Example 3: Profile-profile score of two nucleotide alignments (muscle -ppscore)

**Tool Call:**
```json
{
  "tool": "muscle_msa_command",
  "arguments": {
    "mode": "ppscore",
    "alignment": [
      "ACGTACGTAGCTAGCT",
      "ACGTACGAAGCTAGCT"
    ],
    "alignment2": [
      "ACGTACGTAGCTAGCT",
      "ACGAACGTAGCTAGCT"
    ],
    "type": "dna"
  }
}
```

**Response:**
```json
{
  "alignedSequences": [],
  "score": 2274
}
```

## Verification status

Ported from the source of MUSCLE 3.8.1551 (`-profile`, `-refine`, `-spscore`, `-ppscore`, `-maketree`, `-refinew`, `-sw`; `-profdb` crashes in the program and is always rejected) and verified against the program's binary in the user-run Docker rounds (committed oracle `muscle3_r2_oracle.tsv`: 630 / 630 recorded sets equal). The examples above are produced by this tool.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.MuscleProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.MuscleProgram.cs)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
- Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5): 1792-1797.
