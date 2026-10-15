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

Runs a digest and reports summary statistics instead of the full fragment sequences: total fragment count, fragment sizes (descending), largest and smallest fragment, average fragment size, the enzymes used, and the topology digested. Useful for a quick gel-like readout.

**Topology.** `topology=Linear` (the default) gives `k+1` fragments for `k` cut positions; `topology=Circular` searches for origin-spanning sites as well and gives exactly `k`. The fragment sizes always sum to the molecule length in both topologies (Addgene Restriction Digest Protocol: "the sum of the individual fragments should equal the size of the original"). See [restriction_digest](restriction_digest.md) for the full convention.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L950](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L950)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to digest (non-empty). |
| `enzyme_names` | string[] | Yes | Enzyme names to use (≥1). |
| `topology` | string | No | `Linear` (default) or `Circular`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalFragments` | integer | Number of fragments. |
| `fragmentSizes` | integer[] | Fragment sizes, descending. |
| `largestFragment` / `smallestFragment` | integer | Extremes. |
| `averageFragmentSize` | number | Mean fragment size. |
| `enzymesUsed` | string[] | Enzymes requested. |
| `topology` | string | Topology digested (`Linear` or `Circular`). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | At least one enzyme name is required |

## Examples

### Example 1: EcoRI on linear `AAAGAATTCAAA` → 2 fragments, sizes `[8, 4]`, average 6.

### Example 2: EcoRI on circular `TTCAAAAAAAGAA` → 1 fragment, sizes `[13]`, average 13 (the single `GAATTC` site spans the origin and linearizes the plasmid).

## See Also

- [restriction_digest](restriction_digest.md)
