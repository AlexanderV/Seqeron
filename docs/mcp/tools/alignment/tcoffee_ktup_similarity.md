# tcoffee_ktup_similarity

The k-tuple similarity matrix of the T-Coffee guide-tree stage.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `tcoffee_ktup_similarity` |
| **Method ID** | `SequenceAligner.TCoffeeKtupSimilarity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The pairwise k-tuple similarity matrix T-Coffee 13.41 computes for its guide tree (-distance_matrix_mode ktup: shared k-tuples of the two sequences, k = 1 for protein / 2 for DNA reduced alphabet as in the program), as the integer matrix the program stores (symmetric; the diagonal is the self similarity). Verified together with the guide tree (see tcoffee_guide_tree).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.TCoffeeTree.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.TCoffeeTree.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences (at least one, none empty; letters only). |
| `protein` | boolean | No | true = protein, false (default) = DNA. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `similarity` | array | Symmetric integer matrix (the program stores 100 - min(pscore)) |

## Errors

| Code | Message |
|------|---------|
| 1001 | sequences must contain at least one sequence |

## Examples

### Example 1: Identical and disjoint

**User Prompt:**
> k-tuple similarity of AAAAAAA, AAAAAAA, CCCCCCC.

**Expected Tool Call:**
```json
{"tool": "tcoffee_ktup_similarity", "arguments": {"sequences": ["AAAAAAA", "AAAAAAA", "CCCCCCC"]}}
```

**Response:**
```json
{"similarity": [[100, 100, 70], [100, 100, 70], [70, 70, 100]]}
```

Identical sequences are 100, sequences without a common tuple 70.

## See Also

- [tcoffee_guide_tree](tcoffee_guide_tree.md)
