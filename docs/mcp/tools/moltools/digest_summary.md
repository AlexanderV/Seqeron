# digest_summary

Aggregate statistics over a simulated restriction digest.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `digest_summary` |
| **Method ID** | `RestrictionAnalyzer.GetDigestSummary` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Runs a linear digest and reports summary statistics instead of the full fragment sequences: total fragment count, fragment sizes (descending), largest and smallest fragment, average fragment size, and the enzymes used. Useful for a quick gel-like readout.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L442](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L442)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to digest (non-empty). |
| `enzyme_names` | string[] | Yes | Enzyme names to use (≥1). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalFragments` | integer | Number of fragments. |
| `fragmentSizes` | integer[] | Fragment sizes, descending. |
| `largestFragment` / `smallestFragment` | integer | Extremes. |
| `averageFragmentSize` | number | Mean fragment size. |
| `enzymesUsed` | string[] | Enzymes requested. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | At least one enzyme name is required |

## Examples

### Example 1: EcoRI on `AAAGAATTCAAA` → 2 fragments, sizes `[8, 4]`, average 6.

## See Also

- [restriction_digest](restriction_digest.md)
