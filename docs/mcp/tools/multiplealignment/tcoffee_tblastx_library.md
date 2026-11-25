# tcoffee_tblastx_library

T-Coffee tblastx_msa library built from the XML of an all-against-all tblastx search.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MultipleAlignment |
| **Tool Name** | `tcoffee_tblastx_library` |
| **Method ID** | `Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeTblastxLibraryEntries` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The library t_coffee -method tblastx_msa builds (tc_generic_method.pl -mode=tblastx_msa = seq2tblastx_lib) from the XML of the all-against-all tblastx search of nucleotide sequences: for every hit of a later sequence the nucleotide pairs of every high-scoring pair (the query / hit ranges of the translated frames, minus-strand frames mapped back to the plus strand, weight = the HSP's percent identity, integer truncated), an entry and its mirror, in the printed format of tcoffee_primary_library (0-based positions, first sequence index smaller than the second). The tblastx search is the external program and is NOT run here: pass its XML in blastXml (one Iteration per sequence in input order, sequences identified by their 1-based position as the program numbers them). VERIFICATION: the library equals the one the installed t_coffee writes on every committed tblastx_msa set with a library (scripts/alignment/generate_tcoffee_external_oracle.py).

## Core Documentation Reference

- Source: [Seqeron.Genomics.TCoffee/TCoffeeAligner.Blast.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.TCoffee/TCoffeeAligner.Blast.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)
- Oracle: `scripts/alignment/generate_tcoffee_external_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/` (tcoffee_external_oracle.tsv)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array of string | Yes | The nucleotide sequences (at least one, none empty). |
| `blastXml` | string | Yes | The tblastx XML output (-m7 / -outfmt 5) of the all-against-all search; one Iteration per sequence in input order. |

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

### Example 1: One HSP of a two-sequence tblastx search

**Tool Call:**
```json
{
  "tool": "tcoffee_tblastx_library",
  "arguments": {
    "sequences": [
      "ACGTACGTACGT",
      "ACGTACGTACGT"
    ],
    "blastXml": "<BlastOutput><BlastOutput_iterations><Iteration><Iteration_iter-num>1</Iteration_iter-num><Iteration_query-len>12</Iteration_query-len><Iteration_hits><Hit><Hit_num>1</Hit_num><Hit_accession>1</Hit_accession><Hit_len>12</Hit_len><Hit_hsps><Hsp><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>6</Hsp_query-to><Hsp_query-frame>1</Hsp_query-frame><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>6</Hsp_hit-to><Hsp_hit-frame>1</Hsp_hit-frame><Hsp_identity>2</Hsp_identity><Hsp_align-len>2</Hsp_align-len></Hsp></Hit_hsps></Hit></Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>"
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
      "weight": 100,
      "count": 0
    },
    {
      "seqA": 0,
      "posA": 1,
      "seqB": 1,
      "posB": 1,
      "weight": 100,
      "count": 0
    },
    {
      "seqA": 0,
      "posA": 2,
      "seqB": 1,
      "posB": 2,
      "weight": 100,
      "count": 0
    },
    {
      "seqA": 0,
      "posA": 3,
      "seqB": 1,
      "posB": 3,
      "weight": 100,
      "count": 0
    },
    {
      "seqA": 0,
      "posA": 4,
      "seqB": 1,
      "posB": 4,
      "weight": 100,
      "count": 0
    }
  ],
  "count": 5
}
```

## See Also

- [tcoffee_primary_library](tcoffee_primary_library.md)
- [tcoffee_blast_msa_library](tcoffee_blast_msa_library.md)
- [multiple_align_tcoffee](multiple_align_tcoffee.md)
