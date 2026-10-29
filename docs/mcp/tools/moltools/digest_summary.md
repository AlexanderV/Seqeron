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

**Methylation (opt-in).** `methylation` (`none`|`dam`|`dcm`|`dam+dcm`|`emboss`) masks the bases the chosen methylase set methylates (both strands; the `N` spacer of EcoKI included) and drops every enzyme window containing one — the semantics of EMBOSS `restrict -methylation -mfile Emethylsites.dat`, confirmed against the binary (RESTR-FIND-001 Evidence §7). Default `none` leaves the output unchanged. Example: `GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA` with `MboI` and `methylation=dam` reports no site (without it, positions 10 and 36); EcoRI at 0 and 18 is unaffected.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L1411](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1411)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to digest (non-empty). |
| `enzyme_names` | string[] | Yes | Enzyme names to use (≥1). |
| `topology` | string | No | `Linear` (default) or `Circular`. |
| `methylation` | string | No | `none` (default), `dam` (`GmATC`), `dcm` (`mCmCAGG`, `mCmCTGG`), `dam+dcm` or `emboss` (all five `Emethylsites.dat` entries incl. EcoKI). Case-insensitive; an unknown value is an error. A recognition window containing a methylated base is not a site (EMBOSS `restrict -methylation`; a window rule, not an enzyme-sensitivity table — BamHI is blocked by Dam). |

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
