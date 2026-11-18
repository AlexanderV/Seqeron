# multiple_align_exact

EXACT sum-of-pairs optimal multiple alignment: no other alignment of the sequences has a larger totalScore under the chosen gap model. gapModel='linear' (default): totalScore = column sum of pairs (residue pair match/mismatch, residue against gap gapExtend, gap against gap 0; gapOpen ignored). gapModel='affine': natural affine sum of pairs (Altschul 1989: every pair of rows is projected out of the alignment, gap-gap columns skipped, every gap run of the projection pays gapOpen + k*gapExtend) - the objective of the heuristic multiple_align* tools, which are all heuristics for three or more sequences (their totalScore never exceeds this tool's; the star / progressive / iterative aligners reach it on a large share of small related sets). Two sequences = global_align / global_align_affine; three or more = the k-dimensional dynamic program over the 2^k - 1 neighbours of every cell (Altschul & Lipman 1989; for k = 3 the classic 7-neighbour 3-D recurrence) with, for affine gaps, one value per gap-history state of every cell (the exact state-per-pair program, 13 / 75 / 541 states for k = 3 / 4 / 5), pruned by the Carrillo-Lipman (1988) bound: a cell is dropped when the sum over all pairs of the best pairwise alignment through it (forward + backward Gotoh tables) cannot reach the sum-of-pairs score of the best heuristic alignment (star, progressive or iterative); carrilloLipman=false evaluates the whole lattice (identical result, used as the oracle). The problem is NP-complete in the number of sequences: at most 8 sequences (linear) / 6 (affine), and the dynamic program is stopped with a clear error when it would hold more than maxStateCells values (cells x gap states; default 50 000 000, about 200 MB) - similar sequences are cheap (a thin corridor), unrelated ones explode. Verified equal to exhaustive enumeration of every alignment of tiny sets and to two independent Python exact programs on 804 random sets (k = 2..5), both objectives. Returns the alignment rows (input order), majority consensus, the optimal totalScore, lowerBound (score of the best heuristic seed), latticeCells (product of length+1), visitedCells and stateValues (size actually solved).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_exact` |
| **Method ID** | `SequenceAligner.MultipleAlignExact` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EXACT sum-of-pairs optimal multiple alignment: no other alignment of the sequences has a larger totalScore under the chosen gap model. gapModel='linear' (default): totalScore = column sum of pairs (residue pair match/mismatch, residue against gap gapExtend, gap against gap 0; gapOpen ignored). gapModel='affine': natural affine sum of pairs (Altschul 1989: every pair of rows is projected out of the alignment, gap-gap columns skipped, every gap run of the projection pays gapOpen + k*gapExtend) - the objective of the heuristic multiple_align* tools, which are all heuristics for three or more sequences (their totalScore never exceeds this tool's; the star / progressive / iterative aligners reach it on a large share of small related sets). Two sequences = global_align / global_align_affine; three or more = the k-dimensional dynamic program over the 2^k - 1 neighbours of every cell (Altschul & Lipman 1989; for k = 3 the classic 7-neighbour 3-D recurrence) with, for affine gaps, one value per gap-history state of every cell (the exact state-per-pair program, 13 / 75 / 541 states for k = 3 / 4 / 5), pruned by the Carrillo-Lipman (1988) bound: a cell is dropped when the sum over all pairs of the best pairwise alignment through it (forward + backward Gotoh tables) cannot reach the sum-of-pairs score of the best heuristic alignment (star, progressive or iterative); carrilloLipman=false evaluates the whole lattice (identical result, used as the oracle). The problem is NP-complete in the number of sequences: at most 8 sequences (linear) / 6 (affine), and the dynamic program is stopped with a clear error when it would hold more than maxStateCells values (cells x gap states; default 50 000 000, about 200 MB) - similar sequences are cheap (a thin corridor), unrelated ones explode. Verified equal to exhaustive enumeration of every alignment of tiny sets and to two independent Python exact programs on 804 random sets (k = 2..5), both objectives. Returns the alignment rows (input order), majority consensus, the optimal totalScore, lowerBound (score of the best heuristic seed), latticeCells (product of length+1), visitedCells and stateValues (size actually solved).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.MultipleExact.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.MultipleExact.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section 4.7)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | DNA sequences to align (at least one). |
| `match` | integer | No | Score for matching base. (default: 1) |
| `mismatch` | integer | No | Score for mismatching base. (default: -1) |
| `gapOpen` | integer | No | Gap-open penalty. Used only with gapModel='affine' (a gap of length k then costs gapOpen + k * gapExtend, Gotoh); ignored by the default linear model. (default: -2) |
| `gapExtend` | integer | No | Gap penalty per gap position (linear model: every gap position; affine model: the per-position part of gapOpen + k * gapExtend). (default: -1) |
| `gapModel` | string | No | 'linear' (default; column sum of pairs, gapOpen ignored) or 'affine' (natural affine sum of pairs, gap of length k = gapOpen + k*gapExtend). (default: linear) |
| `carrilloLipman` | boolean | No | true (default): prune the dynamic-programming lattice with the Carrillo-Lipman bound; false: evaluate every cell (same result, exponentially slower). (default: true) |
| `maxStateCells` | integer | No | Complexity guard: largest number of dynamic-programming values (surviving cells x gap states) the exact program may hold; above it the call fails with a clear error instead of exhausting memory. (default: 50000000) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The optimal alignment rows (input order). |
| `consensus` | string | Majority-vote consensus (gaps vote; a nucleotide beats a gap on ties). |
| `totalScore` | integer | The OPTIMAL sum-of-pairs score under the chosen gap model. |
| `lowerBound` | integer | Sum-of-pairs score of the best heuristic alignment (star / progressive / iterative) that seeds the Carrillo-Lipman bound (0 for fewer than three sequences or carrilloLipman = false). |
| `latticeCells` | integer | Number of cells of the full k-dimensional lattice, product of (length + 1). |
| `visitedCells` | integer | Number of lattice cells the dynamic program evaluated. |
| `stateValues` | integer | Dynamic-programming values held (cells x gap-history states); limited by maxStateCells. |

## Errors

| Code | Message |
|------|---------|
| 1010 | Invalid DNA at sequences[i] |
| 1011 | gapModel must be linear or affine |
| 1012 | Exact multiple alignment supports at most 8 sequences with linear gaps (6 with affine) |
| 1013 | Exact multiple alignment exceeds the complexity limit: more than maxStateCells dynamic-program values |

## Examples

### Example 1

**Tool Call:**
```json
{
  "tool": "multiple_align_exact",
  "arguments": {
    "sequences": [
      "ACGTTGCA",
      "ACGTGCA",
      "ACTTGCAA",
      "ACGTTCA"
    ],
    "match": 2,
    "mismatch": -3,
    "gapOpen": -5,
    "gapExtend": -2,
    "gapModel": "affine"
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGTTGC-A",
    "ACG-TGC-A",
    "ACT-TGCAA",
    "ACG-TTC-A"
  ],
  "consensus": "ACG-TGC-A",
  "totalScore": 12,
  "lowerBound": 12,
  "latticeCells": 5184,
  "visitedCells": 13,
  "stateValues": 975
}
```

## Worked Example

Four related 8-mers: the lattice has 9 x 8 x 9 x 8 = 5 184 cells and the Carrillo-Lipman bound (seeded by the best heuristic alignment, score 12) leaves 13 of them; the natural affine sum of pairs of the optimum is 12 (= the heuristic score here; exhaustive enumeration and two independent Python programs agree).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.MultipleExact.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.MultipleExact.cs)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
