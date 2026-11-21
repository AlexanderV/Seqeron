# show_snps

MUMmer 3.23 show-snps: the substitutions and indels of a nucmer delta file with buffer, distance, conflicts and context.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `show_snps` |
| **Method ID** | `MummerShowSnps.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 show-snps: the single-base substitutions and indels of every alignment of a nucmer delta file: [P1] reference position, [SUB] the two characters ('.' = indel), [P2] query position, [BUFF] the distance to the nearest other difference in the same alignment, [DIST] the distance to the nearest sequence end, [R] [Q] the number of other alignments over the position (showConflict = false is -C: only unambiguous SNPs and no [R] [Q] columns), optional sequence lengths (-l), context (-x characters each side), frames and ids; a difference between compatible IUPAC codes (N, R/A ...) is not a SNP (the program's asymmetric CompareIUPAC table). sortByReference (-r, default) / sortByQuery (-q), tabular (-T), header = false (-H), showIndels = false (-I), alignments (-S): only the alignments of these show-coords lines. Byte for byte equal to the installed show-snps on 450 + runs; a promer delta file is refused with NotSupportedException (registered A3-9d).

## Core Documentation Reference

- Source: [MummerShowSnps.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/MummerShowSnps.cs)
- Algorithm doc: [Nucmer_And_Delta_Programs.md](../../../algorithms/Alignment/Nucmer_And_Delta_Programs.md)
- Oracle: [generate_mummer_show_oracle.py](../../../../scripts/alignment/generate_mummer_show_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | What show-snps writes to standard output. |
| `messages` | string | What it writes to standard error (warnings). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Delta file text cannot be null or empty |
| 1001 | ERROR: SNP context must be a positive int |

## Examples

### Example 1: Three substitutions with two bases of context

**User Prompt:**
> List the SNPs of this alignment with 2 bases of context.

**Expected Tool Call:**
```json
{"tool": "show_snps", "arguments": {"delta": "ref.fa qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n", "referenceFasta": ">ref\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n", "queryFasta": ">qry\nAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTGGGGAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n", "context": 2}}
```

**Response:**
```json
{"text": "ref.fa qry.fa\nNUCMER\n\n    [P1]  [SUB]  [P2]      |   [BUFF]   [DIST]  |  [R]  [Q]  |  [CTX R]  [CTX Q]  | [FRM]  [TAGS]\n=====================================================================================================\n      71   A G   66        |        1       66  |    0    0  |    TTAAG    TTGGG  |  1  1  ref\tqry\n      72   A G   67        |        1       67  |    0    0  |    TAAGT    TGGGG  |  1  1  ref\tqry\n      74   T G   69        |        2       67  |    0    0  |    AGTAA    GGGAA  |  1  1  ref\tqry\n", "messages": ""}
```

show-snps -x 2 d1.delta.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
