# find_palindromes

DNA palindromes (restriction-site candidates).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_palindromes` |
| **Method ID** | `RepeatFinder.FindPalindromes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **DNA palindromes** — even-length subsequences that read the same 5'→3' on both
strands (i.e. identical to their reverse complement). These are the recognition sites
of many restriction enzymes (e.g. EcoRI `GAATTC`). `minLength` must be even and ≥ 4.
Results are enumerated by increasing length, then by position.

## Core Documentation Reference

- Source: [RepeatFinder.cs#L981](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L981)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minLength` | integer | No | Minimum palindrome length, even, ≥ 4 (default 4) |
| `maxLength` | integer | No | Maximum palindrome length (default 12) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Palindromes: `{ position, sequence, length }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | minLength must be even and ≥ 4 |

## Examples

### Example 1: EcoRI site (GAATTC)

**User Prompt:**
> Find restriction-site palindromes in "GAATTC".

**Expected Tool Call:**
```json
{
  "tool": "find_palindromes",
  "arguments": { "sequence": "GAATTC", "minLength": 4, "maxLength": 12 }
}
```

**Response:**
```json
{ "items": [ { "position": 1, "sequence": "AATT", "length": 4 }, { "position": 0, "sequence": "GAATTC", "length": 6 } ] }
```
AATT (length 4) and the full EcoRI site GAATTC (length 6) are both palindromic.

### Example 2: No palindrome

**User Prompt:**
> Palindromes in "AAAA"?

**Expected Tool Call:**
```json
{
  "tool": "find_palindromes",
  "arguments": { "sequence": "AAAA", "minLength": 4, "maxLength": 12 }
}
```

**Response:**
```json
{ "items": [] }
```
The reverse complement of AAAA is TTTT, so it is not palindromic.

## Performance

- **Time Complexity:** O(n · (maxLength − minLength)).
- **Space Complexity:** O(number of palindromes).

## See Also

- [find_inverted_repeats](find_inverted_repeats.md) — hairpin-forming inverted repeats
- [find_repeats](find_repeats.md) — any repeated substring
