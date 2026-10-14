# find_restriction_sites

Find restriction sites for named enzymes on both strands.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_restriction_sites` |
| **Method ID** | `RestrictionAnalyzer.FindSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a DNA sequence on both strands for the recognition sites of the named built-in enzymes. IUPAC degenerate codes in recognition sequences are matched against ACGT input.

`position` and `cutPosition` are 0-based **forward-strand** coordinates for both strands. `cutPosition` is the number of forward-strand bases 5' of the cut the enzyme makes on that record's own strand, so a palindromic enzyme reports a forward record (top-strand cut) and a reverse record (bottom-strand cut) at the same `position`, and `recognizedSequence` is the matched window as read on that record's strand. Biopython `Bio.Restriction.<enzyme>.search()` reports the same top-strand cut 1-based, i.e. as `cutPosition + 1`. Overlapping sites are all reported.

With `topology = Linear` (the default) a site is reported only when **both** of its strand cuts fall on a phosphodiester bond that exists, i.e. both lie in `[1, length - 1]` — the same filter as Biopython's `AbstractCut._drop`. This excludes only an enzyme that cuts at the edge of its own recognition window (`Sau3AI` / `MboI`, `^GATC`) sitting at the very 5' or 3' end. With `topology = Circular` the recognition window may span the origin, cut positions are reduced modulo the sequence length, and nothing is dropped.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L237](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L237)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to scan (non-empty). |
| `enzyme_names` | string[] | Yes | Enzyme names to look for (≥1); resolved case-insensitively. |
| `topology` | string | No | `Linear` (default) or `Circular`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites` | array | Sites, each `{position, isForwardStrand, cutPosition, recognizedSequence}`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | At least one enzyme name is required |
| 4001 | Unknown enzyme name |

## Examples

### Example 1: EcoRI

`AAAGAATTCAAA` + `EcoRI` → one forward site at position 3 (cut 4) and one reverse site at position 3 (cut 8).

### Example 2: No site → `sites: []`.

### Example 3: Circular plasmid

`TTCAAAAAAAGAA` + `EcoRI` + `topology: Circular` → one forward site at position 10 (cut 11), spanning the origin; the same input with `topology: Linear` returns `sites: []`.

### Example 4: Linear cut feasibility

`GATCGATCGATC` + `Sau3AI` → one forward site at position 4 (cut 4). The GATC windows at 0 and 8 are not cuttable on a linear molecule (the top-strand cut of the first, and the bottom-strand cut of the last, would fall outside it) — matching `Bio.Restriction.Sau3AI.search(Seq("GATCGATCGATC")) == [5]`.

## See Also

- [find_all_restriction_sites](find_all_restriction_sites.md), [restriction_digest](restriction_digest.md), [restriction_map](restriction_map.md)
