# find_pattern_occurrences

Find every position where an IUPAC recognition pattern occurs, on one strand.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_pattern_occurrences` |
| **Method ID** | `RestrictionAnalyzer.FindPatternOccurrences` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The raw pattern scan underneath `find_restriction_sites`, exposed for a caller who holds a
**pattern** rather than an enzyme. It reports occurrences only: no enzyme, no cut arithmetic, no
reverse-strand scan, and no linear-molecule cut-feasibility filter. Positions are 0-based start
indices in ascending order, and overlapping occurrences are all reported.

**Matching rule (Biopython 1.88 `Bio.Restriction`).** A pattern `N` matches any of the 15 IUPAC
nucleotide codes in the sequence — an `N` included, so `find_pattern_occurrences("GGCCNNNNNGGCC",
"GGCCNNNNNGGCC")` yields `[0]`, exactly as `SfiI.search(Seq("GGCCNNNNNGGCC")) == [9]` (1-based).
Every other pattern code matches only a concrete A/C/G/T base it expands to, so an ambiguity code in
the **sequence** blocks a concrete pattern code rather than widening it: `"GAANTC"` does not match
`GAATTC`, and `"GTYRAC"` does not match `GTYRAC`. A character outside the 15 codes matches nothing.
Both arguments are case-insensitive.

**Why it exists over MCP (B08 H27/H28).** Every other restriction tool used to build a
`DnaSequence`, which accepts only A/C/G/T, so the IUPAC-tolerant scan the library documents was
unreachable: scanning an N-containing contig or a soft-masked assembly was an error, not a result.
Those tools now take the sequence as it is, and this tool exposes the pattern scan itself, which had
no MCP expression at all.

**One strand only.** To cover both, call again with the reverse complement of the pattern. For an
enzyme — both strands, cut positions, the feasibility filter and the circular wrap — use
`find_restriction_sites` instead.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L603](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L603)
- Algorithm: [Restriction_Site_Detection.md](../../../algorithms/MolTools/Restriction_Site_Detection.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string \| null | Yes | Sequence to scan; IUPAC codes accepted anywhere, null or empty yields no occurrences. |
| `recognition_pattern` | string | Yes | 1 or more of the 15 IUPAC nucleotide codes. |
| `start_index` | integer | No | First index to consider (default 0, must not be negative). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `positions` | integer[] | 0-based start indices on the given strand, ascending; overlaps included. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Recognition pattern cannot be null or empty. |
| 1002 | Recognition pattern '&lt;pattern&gt;' contains '&lt;char&gt;' at index &lt;i&gt;, which is not one of the 15 IUPAC nucleotide codes |
| 1003 | startIndex ('-1') must be a non-negative value. |

## Examples

### Example 1: `AAAGAATTCAAAGAATTC` + `GAATTC` → `{ positions: [3, 12] }`.

### Example 2: `GGCCNNNNNGGCC` + `GGCCNNNNNGGCC` → `{ positions: [0] }` — a pattern `N` matches a sequence `N`, which is what Biopython does (`SfiI.search` finds this site).

### Example 3: `GAANTC` + `GAATTC` → `{ positions: [] }` — the sequence's `N` blocks the concrete `T`.

### Example 4: `AAAA` + `AA` → `{ positions: [0, 1, 2] }` — overlapping occurrences are all reported.

### Example 5: `AAAA` + `AA` with `start_index: 2` → `{ positions: [2] }`.

### Example 6: `null` + `GAATTC` → `{ positions: [] }`; the pattern is still validated first, so `null` + `GAZTC` is an error.

## See Also

- [find_restriction_sites](find_restriction_sites.md) - The enzyme-based scan: both strands, cut positions, topology.
- [find_sites_custom_enzyme](find_sites_custom_enzyme.md) - A caller-supplied enzyme definition rather than a bare pattern.
