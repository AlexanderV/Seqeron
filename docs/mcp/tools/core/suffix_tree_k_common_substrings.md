# suffix_tree_k_common_substrings

Find the longest substring(s) common to k strings, or present in at least minSupport of them.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_k_common_substrings` |
| **Method ID** | `SuffixTree.FindLongestCommonSubstrings` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds every distinct longest substring that occurs in at least `minSupport` of `texts` (default: all of them — the longest common substring of k strings, Rosalind LCSM), sorted ordinally. The texts are joined with k distinct separator characters absent from every input into a generalized suffix tree; each internal node's number of distinct source texts (Hui 1992 colour-set size) gives the k-common substring lengths of Gusfield 1997 §7.6. `lengthsBySupport[q-1]` is l(q), the length of the longest substring present in at least q texts (q = 1..k; l(1) is the longest text). Cross-checked against brute force on random sets and the Rosalind LCSM sample.

## Core Documentation Reference

- Source: [SuffixTree/SuffixTree.Repeats.cs#L27](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Repeats.cs#L27)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `texts` | array<string> | Yes | The input texts (at least one; none null) |
| `minSupport` | integer | No | Minimum number of texts that must contain the substring (1..k; default k = all) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `substrings` | array<string> | Every longest substring with the required support, sorted ordinally |
| `length` | integer | Their length (0 when none) |
| `minSupport` | integer | The support threshold used |
| `lengthsBySupport` | array<integer> | Element q-1 = longest length present in at least q texts |

## Errors

| Code | Message |
|------|---------|
| 1001 | Texts cannot be null or empty |
| 1002 | Texts cannot contain null |
| 1003 | minSupport must be between 1 and the number of texts |

## Examples

### Example 1: Rosalind LCSM sample

**User Prompt:**
> What is the longest common substring of GATTACA, TAGACCA and ATACA?

**Tool Call:**
```json
{
  "tool": "suffix_tree_k_common_substrings",
  "arguments": {
    "texts": [
      "GATTACA",
      "TAGACCA",
      "ATACA"
    ]
  }
}
```

**Response:**
```json
{
  "substrings": [
    "AC",
    "CA",
    "TA"
  ],
  "length": 2,
  "minSupport": 3,
  "lengthsBySupport": [
    7,
    4,
    2
  ]
}
```

### Example 2: Present in at least two

**User Prompt:**
> Longest substring shared by at least two of GATTACA, TAGACCA, ATACA?

**Tool Call:**
```json
{
  "tool": "suffix_tree_k_common_substrings",
  "arguments": {
    "texts": [
      "GATTACA",
      "TAGACCA",
      "ATACA"
    ],
    "minSupport": 2
  }
}
```

**Response:**
```json
{
  "substrings": [
    "TACA"
  ],
  "length": 4,
  "minSupport": 2,
  "lengthsBySupport": [
    7,
    4,
    2
  ]
}
```

## Worked Example

Rosalind LCSM sample: DATA GATTACA / TAGACCA / ATACA → sample output `AC` (any longest common substring is accepted). The tool returns all three of length 2: AC, CA, TA. `TACA` occurs in GATTACA and ATACA (support 2, l(2) = 4).

## See Also

- [suffix_tree_all_lcs](suffix_tree_all_lcs.md) — All longest common substrings of two texts
- [suffix_tree_lcs](suffix_tree_lcs.md) — One longest common substring

## References

- Algorithm source: [SuffixTree/SuffixTree.Repeats.cs#L27](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Repeats.cs#L27)
- Binding: [SuffixTreeCoreTools.cs](../../../../src/SuffixTree/Mcp/SuffixTree.Mcp.Core/Tools/SuffixTreeCoreTools.cs)
