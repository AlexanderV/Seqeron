# generate_decipher_consensus

Bioconductor DECIPHER `ConsensusSequence` consensus (DNA, RNA or protein).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `generate_decipher_consensus` |
| **Method ID** | `MotifFinder.GenerateDecipherConsensus` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Line-by-line port of DECIPHER 3.9.4 `ConsensusSequence` (`R/ConsensusSequence.R`, `src/ConsensusSequence.c`). Per column the characters are tallied as fractions (IUPAC codes split between their bases when `ambiguity` is true). With t = 1 − `threshold`, the first test that holds in source order chooses the symbol: a single residue strictly most frequent with fraction ≥ t; then (DNA/RNA) `Y K W S R M B D H V` — included bases all strictly above the excluded ones, sum ≥ t — then `N`; (protein) `B`/`Z`/`J`, then `X`. The chosen fraction must be ≥ the gap and mask fractions (else `-`/`+`). Positions whose consensus carries less than `minInformation` (default 1 − threshold) get `noConsensusChar`. Terminal gaps are ignored unless `includeTerminalGaps`; `.` is a gap, `+` a mask; unequal lengths are allowed. `includeNonLetters` behaves as in DECIPHER (passed to the C `ignoreNonLetters`: `true` leaves gaps/masks out). Cross-checked against DECIPHER's own R/C source built with R 4.3.3 + Biostrings 2.70.2: 10,000/10,000 identical.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.DecipherConsensus.cs#L110](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.DecipherConsensus.cs#L110)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Aligned sequences (unequal lengths allowed) |
| `sequenceType` | string | No | 'dna' (default), 'rna' or 'protein' |
| `threshold` | number | No | Fraction of information that may be lost per position, in [0, 1) (default 0.05) |
| `ambiguity` | boolean | No | Split IUPAC degeneracy codes between their residues (default true) |
| `noConsensusChar` | string | No | Single character of the alphabet for positions without consensus (default '+') |
| `minInformation` | number | No | Minimum fraction of information in (0, 1]; default 1 − threshold |
| `includeNonLetters` | boolean | No | DECIPHER includeNonLetters (default false: gaps/masks counted; true: left out) |
| `includeTerminalGaps` | boolean | No | Count leading/trailing gaps (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `consensus` | string | DECIPHER consensus (length of the longest sequence) |

## Errors

| Code | Message |
|------|---------|
| 1001 | At least one sequence is required |
| 1002 | noConsensusChar must be exactly one character |
| 1003 | sequenceType must be 'dna', 'rna' or 'protein' |
| — | threshold ∉ [0, 1), minInformation ∉ (0, 1] (ArgumentOutOfRangeException); character or noConsensusChar outside the alphabet, null element (ArgumentException) |

## Examples

### Example 1: Default threshold

**User Prompt:**
> DECIPHER consensus of A, A, A, T.

**Tool Call:**
```json
{
  "tool": "generate_decipher_consensus",
  "arguments": {
    "sequences": ["A", "A", "A", "T"]
  }
}
```

**Response:**
```json
{
  "consensus": "W"
}
```

### Example 2: Majority-based consensus (threshold 0.5)

**Tool Call:**
```json
{
  "tool": "generate_decipher_consensus",
  "arguments": {
    "sequences": ["GTT", "GAA", "CTG"],
    "threshold": 0.5
  }
}
```

**Response:**
```json
{
  "consensus": "GTD"
}
```

### Example 3: Protein

**Tool Call:**
```json
{
  "tool": "generate_decipher_consensus",
  "arguments": {
    "sequences": ["ANQIH-", "ADELW."],
    "sequenceType": "protein"
  }
}
```

**Response:**
```json
{
  "consensus": "ABZJX-"
}
```

## Worked Example

A, A, A, T: A = 0.75 < 0.95, T alone fails; W (A + T = 1 ≥ 0.95, both above C and G) passes → `W`. With `threshold: 0.3` (t = 0.7) A alone passes → `A`; adding `minInformation: 0.8` rejects it (0.75 < 0.8) → `+`.

## See Also

- [generate_consensus](generate_consensus.md) — per-base inclusion-threshold IUPAC consensus
- [generate_cavener_consensus](generate_cavener_consensus.md) — Cavener 1987 degenerate consensus
- [generate_emboss_consensus](generate_emboss_consensus.md) — EMBOSS cons consensus

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.DecipherConsensus.cs#L110](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.DecipherConsensus.cs#L110)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
- Wright E.S. DECIPHER `ConsensusSequence` (Bioconductor), source `src/ConsensusSequence.c`.
