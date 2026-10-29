# restriction_digest_custom_enzyme

Digest a sequence with a caller-supplied enzyme definition.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `restriction_digest_custom_enzyme` |
| **Method ID** | `RestrictionAnalyzer.Digest` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

`restriction_digest` for an enzyme that cannot be resolved by name — not in the 754-enzyme REBASE
table, or one of the 334 whose cut position REBASE leaves undefined. The enzyme is given exactly as
in [find_sites_custom_enzyme](find_sites_custom_enzyme.md): the recognition sequence plus the two
cut offsets counted from the first base of the site, and optionally a second pair for a two-sided
(Type IIB) cutter.

Fragment conventions are `restriction_digest`'s: fragments in ascending start position, each with
its sequence, start, length, flanking enzymes and the blunt/5'/3' overhang of each cut end. On a
**circular** molecule the site search wraps past the origin and `k` cut positions give `k` fragments
(vs `k + 1` linear) — except `k = 0`, which leaves the circle whole as one fragment. IUPAC codes are
accepted in both the sequence and the recognition sequence (B08 H27).

**Why it exists (B08 H28).** The `IEnumerable<RestrictionEnzyme>` overloads of the C# API had no MCP
expression, so a caller-supplied or cut-undefined enzyme could not be digested over MCP at all.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.RawSequence.cs#L111](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.RawSequence.cs#L111)
- Algorithm: [Restriction_Digest_Simulation.md](../../../algorithms/MolTools/Restriction_Digest_Simulation.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to digest; IUPAC codes accepted. |
| `recognition_sequence` | string | Yes | The custom enzyme's site; IUPAC codes allowed, need not be palindromic. At most 64 bases (REBASE's longest is 15). |
| `cut_position_forward` | integer | Yes | Top-strand cut offset from the first site base. |
| `cut_position_reverse` | integer | Yes | Bottom-strand cut offset in the same frame. The two cuts may be at most 64 bases apart (the stagger, i.e. the overhang length); a larger one is refused. |
| `name` | string \| null | No | Name reported on the flanking-enzyme fields; defaults to `CustomEnzyme`. At most 64 characters. |
| `topology` | string | No | `Linear` (default) or `Circular`. |
| `second_cut_forward` | integer \| null | No | Second top-strand cut offset of a Type IIB cutter. |
| `second_cut_reverse` | integer \| null | No | Second bottom-strand cut offset; same stagger as the first pair. |
| `methylation` | string \| null | No | `none` (default), `dam`, `dcm`, `dam+dcm` or `emboss`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `fragments` | object[] | Fragments in the same shape as `restriction_digest`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Pass both second_cut_forward and second_cut_reverse, or neither. |
| 1002 | Enzyme '&lt;name&gt;' recognition sequence '&lt;site&gt;' contains '&lt;char&gt;' at index &lt;i&gt;, which is not one of the 15 IUPAC nucleotide codes |
| 1003 | Enzyme '&lt;name&gt;' second cut pair (&lt;f&gt;, &lt;r&gt;) leaves a different stagger from its first pair |
| 1004 | Sequence cannot be null or empty. |
| 1005 | Unknown methylation set '{name}'; expected none, dam, dcm, dam+dcm or emboss. |
| 1006 | The topology must be MoleculeTopology.Linear or MoleculeTopology.Circular. |
| 1007 | Enzyme '&lt;name&gt;' cuts its two strands &lt;n&gt; bases apart (&lt;f&gt;, &lt;r&gt;); the largest supported stagger is 64. |
| 1008 | Enzyme '&lt;name&gt;' has a recognition sequence of &lt;n&gt; bases; the longest supported is 64. |
| 1009 | The enzyme name is &lt;n&gt; characters long; the longest supported name is 64. |

## Examples

### Example 1: EcoRI by hand — `AAAGAATTCAAA`, site `GAATTC`, cuts `(1, 5)` → 2 fragments of 4 and 8 nt, the partition `restriction_digest` gives for `EcoRI`.

### Example 2: the same molecule with `topology: Circular` → 1 fragment of 12 nt — one cut on a circle linearises it rather than splitting it.

### Example 3: `second_cut_reverse` without `second_cut_forward` → **error** 1001.

## See Also

- [restriction_digest](restriction_digest.md) - For any of the 754 named REBASE enzymes.
- [find_sites_custom_enzyme](find_sites_custom_enzyme.md) - The sites rather than the fragments.
