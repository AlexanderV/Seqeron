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
of many restriction enzymes (e.g. EcoRI `GAATTC`). `minLength` must be even and ≥ 4;
an odd `maxLength` is rounded down (odd-length palindromes cannot exist).

Every palindromic window is reported, including windows nested inside longer ones
(Rosalind REVP convention). Positions are 0-based (REVP is 1-based). Results are ordered by
position, then length (the REVP sample-output order). Input is case-insensitive; only A/C/G/T
pair, so windows containing N, other IUPAC codes, U or gaps are never reported. For maximal
stems with a loop use `find_inverted_repeats`.

## Core Documentation Reference

- Source: [RepeatFinder.cs#L1266](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L1266)

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
{ "items": [ { "position": 0, "sequence": "GAATTC", "length": 6 }, { "position": 1, "sequence": "AATT", "length": 4 } ] }
```
The full EcoRI site GAATTC (length 6) and the nested AATT (length 4) are both palindromic.

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

### Example 3: Ambiguous bases never pair

**Expected Tool Call:**
```json
{
  "tool": "find_palindromes",
  "arguments": { "sequence": "GAATTCNNNNGAATTC", "minLength": 4, "maxLength": 12 }
}
```

**Response:**
```json
{ "items": [ { "position": 0, "sequence": "GAATTC", "length": 6 }, { "position": 1, "sequence": "AATT", "length": 4 }, { "position": 10, "sequence": "GAATTC", "length": 6 }, { "position": 11, "sequence": "AATT", "length": 4 } ] }
```
`NNNN` equals its IUPAC reverse complement symbolically but is not reported: an N may be any base.

## Performance

- **Time Complexity:** O(n · maxLength) comparisons (capped centre expansion) plus output.
- **Space Complexity:** O(number of palindromes).

## See Also

- [find_inverted_repeats](find_inverted_repeats.md) — hairpin-forming inverted repeats
- [find_repeats](find_repeats.md) — any repeated substring
