# restriction_digest

Simulate a linear restriction digest.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `restriction_digest` |
| **Method ID** | `RestrictionAnalyzer.Digest` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Cuts a linear DNA molecule at the forward-strand cut positions of the named enzymes and returns the resulting fragments in 5′→3′ order. Each fragment carries its sequence, start position, length, flanking enzymes (`leftEnzyme`/`rightEnzyme`, null at the molecule ends), and fragment number. With `k` distinct forward-strand cut positions a linear molecule yields `k+1` fragments; with zero cuts it returns the whole molecule as one fragment.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L259](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L259)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to digest (non-empty). |
| `enzyme_names` | string[] | Yes | Enzyme names to use (≥1). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `fragments` | array | Fragments, each `{sequence, startPosition, length, leftEnzyme, rightEnzyme, fragmentNumber}`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | At least one enzyme name is required |

## Examples

### Example 1: EcoRI on `AAAGAATTCAAA`

Cut at position 4 → `AAAG` (0–4) and `AATTCAAA` (4–12).

### Example 2: No cut site → single whole-sequence fragment.

## See Also

- [digest_summary](digest_summary.md), [find_restriction_sites](find_restriction_sites.md)
