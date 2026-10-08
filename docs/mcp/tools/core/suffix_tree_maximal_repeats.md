# suffix_tree_maximal_repeats

Find every maximal repeated pair of a text (MUMmer repeat-match -f).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `suffix_tree_maximal_repeats` |
| **Method ID** | `SuffixTree.FindMaximalRepeatedPairs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds every maximal repeated pair of `text`: a triple (firstPosition i < secondPosition j, length L ≥ `minLength`) with text[i..i+L) = text[j..j+L) that can be extended neither to the left (i = 0 or the preceding characters differ) nor to the right (the copy at j ends the text or the following characters differ) — Gusfield 1997 §7.12. Copies may overlap (tandem repeats). Identical to MUMmer 3 `repeat-match -f -n minLength` (forward strand; cross-checked against the MUMmer 3.23 binary on 0.5 M pairs). Characters listed in `uniqueSymbols`, or every character other than upper-case A/C/G/T with `nonAcgtUnique`, never match — not even themselves — as separators / ambiguity codes in RepeatFinder.FindDirectRepeats (cross-checked) and Vmatch; by default every character matches itself like repeat-match (which also matches N with N). Gusfield §7.12.3 bottom-up left-character lists, O(n + z) enumeration; results 0-based, sorted by firstPosition then secondPosition.

## Core Documentation Reference

- Source: [SuffixTree/SuffixTree.Repeats.cs#L17](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Repeats.cs#L17)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `text` | string | Yes | The text to analyze |
| `minLength` | integer | No | Minimum repeat length (>= 1; repeat-match default 20) |
| `uniqueSymbols` | string | No | Characters treated as unique separators that never match, e.g. "N$" (default: none, like repeat-match) |
| `nonAcgtUnique` | boolean | No | Treat every character other than upper-case A, C, G, T as unique (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pairs` | array<object> | Maximal pairs `{firstPosition, secondPosition, length}`, sorted by first then second position |

## Errors

| Code | Message |
|------|---------|
| 1001 | Text cannot be null or empty |
| 1003 | minLength must be >= 1 |

## Examples

### Example 1: repeat-match -f -n 3 (Gusfield xabcyabcwabcyz)

**User Prompt:**
> What are the maximal repeats of length ≥ 3 in xabcyabcwabcyz?

**Tool Call:**
```json
{
  "tool": "suffix_tree_maximal_repeats",
  "arguments": {
    "text": "xabcyabcwabcyz",
    "minLength": 3
  }
}
```

**Response:**
```json
{
  "pairs": [
    {
      "firstPosition": 1,
      "secondPosition": 5,
      "length": 3
    },
    {
      "firstPosition": 1,
      "secondPosition": 9,
      "length": 4
    },
    {
      "firstPosition": 5,
      "secondPosition": 9,
      "length": 3
    }
  ]
}
```

### Example 2: N as a unique symbol

**User Prompt:**
> Find maximal repeats of length ≥ 3 in ACGTNACGTNACGT without matching across N.

**Tool Call:**
```json
{
  "tool": "suffix_tree_maximal_repeats",
  "arguments": {
    "text": "ACGTNACGTNACGT",
    "minLength": 3,
    "uniqueSymbols": "N"
  }
}
```

**Response:**
```json
{
  "pairs": [
    {
      "firstPosition": 0,
      "secondPosition": 5,
      "length": 4
    },
    {
      "firstPosition": 0,
      "secondPosition": 10,
      "length": 4
    },
    {
      "firstPosition": 5,
      "secondPosition": 10,
      "length": 4
    }
  ]
}
```

## Worked Example

MUMmer 3.23 `repeat-match -f -n 3` on `xabcyabcwabcyz` prints `2 10 4 / 6 10 3 / 2 6 3` (1-based Start1, Start2, Length, traversal order); converted to 0-based and sorted this is the first response. On `ACGTNACGTNACGT` repeat-match reports `1 6 9 / 1 11 4` (N matches N); with `uniqueSymbols: "N"` the 9-long N-spanning pair splits into the three ACGT pairs of the second example.

## See Also

- [suffix_tree_all_lrs](suffix_tree_all_lrs.md) — All longest repeated substrings
- [suffix_tree_find_mems](suffix_tree_find_mems.md) — Maximal exact matches between two texts

## References

- Algorithm source: [SuffixTree/SuffixTree.Repeats.cs#L17](../../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Repeats.cs#L17)
- Binding: [SuffixTreeCoreTools.cs](../../../../src/SuffixTree/Mcp/SuffixTree.Mcp.Core/Tools/SuffixTreeCoreTools.cs)
