# tcoffee_guide_tree

The guide tree T-Coffee builds for a set of sequences, as Newick.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MultipleAlignment |
| **Tool Name** | `tcoffee_guide_tree` |
| **Method ID** | `Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeGuideTree` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The guide tree T-Coffee 13.41 builds for a set of sequences with its defaults (-tree_mode nj -distance_matrix_mode ktup): k-tuple similarity distances, the program's neighbour-joining routine and the rooting of its trifurcation into a binary tree. Returns the tree as Newick text (no branch lengths; leaf labels = names, default s0, s1, ...; left subtree = first profile of the merge) and the leaf indices left to right; the progressive merge order of multiple_align_tcoffee is the post-order of this tree. protein omitted = decided from the letters.

## Core Documentation Reference

- Source: [Seqeron.Genomics.TCoffee/TCoffeeAligner.Tree.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.TCoffee/TCoffeeAligner.Tree.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences (at least one, none empty; letters only). |
| `protein` | boolean | No | true = protein, false = DNA, omitted = decided from the letters. |
| `names` | array<string> | No | Leaf labels (default s0, s1, ...); one per sequence. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `newick` | string | Tree in Newick (no branch lengths) |
| `leafOrder` | array | Leaf indices left to right |

## Errors

| Code | Message |
|------|---------|
| 1001 | sequences must contain at least one sequence |
| 1001 | names must have one label per sequence |

## Examples

### Example 1: Four sequences

**User Prompt:**
> T-Coffee guide tree of four sequences.

**Expected Tool Call:**
```json
{"tool": "tcoffee_guide_tree", "arguments": {"sequences": ["ACGTACGTAC", "ACGTACGTAA", "TTTTGGGGCC", "TTTTGGGGCA"]}}
```

**Response:**
```json
{"newick": "((s3,s2),(s0,s1));", "leafOrder": [3, 2, 0, 1]}
```

## See Also

- [tcoffee_ktup_similarity](tcoffee_ktup_similarity.md)
- [multiple_align_tcoffee](multiple_align_tcoffee.md)
