# tcoffee_blast_msa_library

T-Coffee blastp_msa library built from the XML of an all-against-all BLAST search.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MultipleAlignment |
| **Tool Name** | `tcoffee_blast_msa_library` |
| **Method ID** | `SequenceAligner.TCoffeeBlastMsaLibraryEntries` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The library t_coffee -method blastp_msa builds (-lib_only -out_lib) from the XML of the all-against-all BLAST search of the sequences (blastall -p blastp -m7, or BLAST+ -outfmt 5): for every high-scoring pair (Hsp_qseq / Hsp_hseq) the residue pairs of its gap-free alignment columns with the HSP's percent identity (integer, truncated) as the weight, an entry and its mirror, in the printed format of tcoffee_primary_library (0-based positions, sequence index of the first smaller than that of the second, self hits not printed). The BLAST search itself is the external program and is NOT run here: pass its XML in blastXml (one Iteration per sequence in input order; the hits are matched to the sequences by their Hit_def name). VERIFICATION: the library equals the one the installed t_coffee writes on every committed set with a library (70 or more oracle sets, scripts/alignment/generate_tcoffee_external_oracle.py).

## Core Documentation Reference

- Source: [Seqeron.Genomics.TCoffee/TCoffeeAligner.Blast.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.TCoffee/TCoffeeAligner.Blast.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)
- Oracle: `scripts/alignment/generate_tcoffee_external_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/` (tcoffee_external_oracle.tsv)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `names` | array of string | Yes | Sequence names exactly as they stand in the BLAST database / query (one per sequence); the hits of the XML are matched by name. |
| `sequences` | array of string | Yes | The sequences (at least one, none empty). |
| `blastXml` | string | Yes | The BLAST XML output (-m7 / -outfmt 5) of the all-against-all search; one Iteration per sequence in input order. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entries` | array | Library entries: seqA, posA, seqB, posB (0-based), weight, count. |
| `count` | integer | Number of entries. |

## Errors

| Code | Message |
|------|---------|
| 1001 | a sequence argument is null or empty / names and matrix have the wrong shape |
| 1002 | a gap penalty is outside the range of the program / unknown matrix |

## Examples

### Example 1: One HSP of a two-protein BLAST search

**Tool Call:**
```json
{
  "tool": "tcoffee_blast_msa_library",
  "arguments": {
    "names": [
      "a",
      "b"
    ],
    "sequences": [
      "ACDEF",
      "ACDF"
    ],
    "blastXml": "<BlastOutput><BlastOutput_iterations><Iteration><Iteration_iter-num>1</Iteration_iter-num><Iteration_hits><Hit><Hit_num>1</Hit_num><Hit_def>b</Hit_def><Hit_hsps><Hsp><Hsp_query-from>1</Hsp_query-from><Hsp_hit-from>1</Hsp_hit-from><Hsp_identity>4</Hsp_identity><Hsp_align-len>5</Hsp_align-len><Hsp_qseq>ACDEF</Hsp_qseq><Hsp_hseq>ACD-F</Hsp_hseq></Hsp></Hit_hsps></Hit></Iteration_hits></Iteration><Iteration><Iteration_iter-num>2</Iteration_iter-num><Iteration_hits></Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>"
  }
}
```

**Response:**
```json
{
  "entries": [
    {
      "seqA": 0,
      "posA": 0,
      "seqB": 1,
      "posB": 0,
      "weight": 80,
      "count": 0
    },
    {
      "seqA": 0,
      "posA": 1,
      "seqB": 1,
      "posB": 1,
      "weight": 80,
      "count": 0
    },
    {
      "seqA": 0,
      "posA": 2,
      "seqB": 1,
      "posB": 2,
      "weight": 80,
      "count": 0
    },
    {
      "seqA": 0,
      "posA": 4,
      "seqB": 1,
      "posB": 3,
      "weight": 80,
      "count": 0
    }
  ],
  "count": 4
}
```

## See Also

- [tcoffee_primary_library](tcoffee_primary_library.md)
- [tcoffee_tblastx_library](tcoffee_tblastx_library.md)
- [multiple_align_tcoffee](multiple_align_tcoffee.md)
