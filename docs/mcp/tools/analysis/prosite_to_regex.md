# prosite_to_regex

Translate a PROSITE pattern to a .NET regex.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `prosite_to_regex` |
| **Method ID** | `ProteinMotifFinder.ConvertPrositeToRegex` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Converts a **PROSITE PA-line pattern** into an equivalent **.NET regular expression**:

| PROSITE | Regex |
|---------|-------|
| `-` (separator) | *(removed)* |
| `x` | `.` |
| `x(n)` / `x(n,m)` | `.{n}` / `.{n,m}` |
| `[AC]` | `[AC]` |
| `{P}` (exclusion) | `[^P]` |
| `A(2)` | `A{2}` |
| `<` (N-term) | `^` |
| `>` (C-term) | `$` |

Unsupported metacharacters raise a format error rather than being silently dropped.

## Core Documentation Reference

- Source: [ProteinMotifFinder.cs#L277](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs#L277)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `prositePattern` | string | Yes | PROSITE-format pattern (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `regex` | string | .NET regex string |

## Errors

| Code | Message |
|------|---------|
| 1001 | PROSITE pattern cannot be null or empty |

## Examples

### Example 1: N-glycosylation site (PS00001)

**User Prompt:**
> Convert the PROSITE pattern "N-{P}-[ST]-{P}" to regex.

**Expected Tool Call:**
```json
{
  "tool": "prosite_to_regex",
  "arguments": { "prositePattern": "N-{P}-[ST]-{P}" }
}
```

**Response:**
```json
{ "regex": "N[^P][ST][^P]" }
```

### Example 2: Variable-gap pattern

**User Prompt:**
> Convert "[AC]-x(2)-V".

**Expected Tool Call:**
```json
{
  "tool": "prosite_to_regex",
  "arguments": { "prositePattern": "[AC]-x(2)-V" }
}
```

**Response:**
```json
{ "regex": "[AC].{2}V" }
```

## Performance

- **Time Complexity:** O(length of the pattern).
- **Space Complexity:** O(length of the pattern).

## See Also

- [find_motif_by_prosite](find_motif_by_prosite.md) — convert and scan
- [find_motif_by_pattern](find_motif_by_pattern.md) — scan with a raw regex
