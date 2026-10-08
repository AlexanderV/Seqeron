# suffix_tree_all_lcs

Find every distinct longest common substring of two texts (all ties) with all positions.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_all_lcs` |
| **Method ID** | `SuffixTree.FindAllDistinctLongestCommonSubstrings` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds every distinct longest common substring of `text1` (on which the suffix tree is built) and `text2` — all length ties, unlike `suffix_tree_lcs`, which returns one. Each substring carries all of its 0-based start positions in `text1` and in `text2` (ascending, overlapping occurrences included); substrings are ordered by first occurrence in `text2` (the first one is the canonical LCS of `suffix_tree_lcs`). Matching statistics (Chang & Lawler 1994; Gusfield 1997 §7.8) in O(|text1| + |text2| + occurrences): every occurrence in text2 of a longest common substring ends where the matching statistic equals the maximum. Empty when no character is shared.

## Core Documentation Reference

- Source: [SuffixTree/SuffixTree.Repeats.cs#L9](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Repeats.cs#L9)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `text1` | string | Yes | The first text (the suffix tree is built on it) |
| `text2` | string | Yes | The second text |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `substrings` | array<object> | `{substring, positionsInText1, positionsInText2}`, ordered by first occurrence in text2 |
| `length` | integer | Common length (0 when none) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Text1 cannot be null or empty |
| 1002 | Text2 cannot be null or empty |

## Examples

### Example 1: Two tied LCS

**User Prompt:**
> Which longest substrings do abcxyzabc and xyzqabcpxyz share, and where?

**Tool Call:**
```json
{
  "tool": "suffix_tree_all_lcs",
  "arguments": {
    "text1": "abcxyzabc",
    "text2": "xyzqabcpxyz"
  }
}
```

**Response:**
```json
{
  "substrings": [
    {
      "substring": "xyz",
      "positionsInText1": [
        3
      ],
      "positionsInText2": [
        0,
        8
      ]
    },
    {
      "substring": "abc",
      "positionsInText1": [
        0,
        6
      ],
      "positionsInText2": [
        4
      ]
    }
  ],
  "length": 3
}
```

## Worked Example

All substrings of `xyzqabcpxyz` that occur in `abcxyzabc` have length ≤ 3; the length-3 ones are `xyz` (text2 positions 0 and 8, text1 position 3) and `abc` (text2 position 4, text1 positions 0 and 6) — a brute-force enumeration gives the same lists. `suffix_tree_lcs` returns only `xyz`.

## See Also

- [suffix_tree_lcs](suffix_tree_lcs.md) — One longest common substring
- [suffix_tree_k_common_substrings](suffix_tree_k_common_substrings.md) — Longest common substring of k strings

## References

- Algorithm source: [SuffixTree/SuffixTree.Repeats.cs#L9](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Repeats.cs#L9)
- Binding: [SuffixTreeCoreTools.cs](../../../../src/SuffixTree/Mcp/SuffixTree.Mcp.Core/Tools/SuffixTreeCoreTools.cs)
