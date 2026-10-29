# find_sites_custom_enzyme

Find the sites of a caller-supplied enzyme definition.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_sites_custom_enzyme` |
| **Method ID** | `RestrictionAnalyzer.FindSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans for an enzyme that is **not resolvable by name**: one that is not in the 754-enzyme REBASE
table at all, or one of the 334 REBASE enzymes whose cut position REBASE leaves undefined (e.g.
`Dpi3084I`) and which therefore cannot be digested by name.

**The cut frame.** Both offsets are counted from the first base of the recognition site:
`cut_position_forward` is the number of site bases 5' of the **top-strand** cut and
`cut_position_reverse` the same for the **bottom-strand** cut. So EcoRI `G^AATTC` / `GAATT^C` is
`(1, 5)`, AluI `AG^CT` is `(2, 2)` (blunt: both offsets equal), and the Type IIS BsaI
`GGTCTCN^NNNN` is `(7, 11)`. Offsets outside the site are how a **distal** cutter is expressed, and
a negative offset cuts upstream of the site. For a palindromic enzyme the two offsets sum to the
site length; a non-palindromic site yields distinct forward- and reverse-strand records, and the
reverse-strand record reports the cut the enzyme makes on the strand it read.

**Two-sided (Type IIB) cutters.** `second_cut_forward` / `second_cut_reverse` add the second cut
pair — BaeI is `(-10, -15)` and `(23, 18)` — and must leave the same stagger as the first pair, as
Biopython stores one `ovhg` per enzyme. Pass both or neither.

Coordinates, the two records per palindromic site, the linear cut-feasibility rule (both cuts must
fall on a bond that exists), the circular wrap and the methylation mode are exactly as in
`find_restriction_sites`. IUPAC codes are accepted in both the sequence and the recognition
sequence.

**Why it exists (B08 H28).** The `IEnumerable<RestrictionEnzyme>` overloads of the C# API had no MCP
expression, so a caller-supplied or cut-undefined enzyme could not be used over MCP at all.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.RawSequence.cs#L74](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.RawSequence.cs#L74)
- Algorithm: [Restriction_Site_Detection.md](../../../algorithms/MolTools/Restriction_Site_Detection.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to scan; IUPAC codes accepted. |
| `recognition_sequence` | string | Yes | The custom enzyme's site; IUPAC codes allowed, need not be palindromic. |
| `cut_position_forward` | integer | Yes | Top-strand cut offset from the first site base. |
| `cut_position_reverse` | integer | Yes | Bottom-strand cut offset in the same frame. |
| `name` | string \| null | No | Name reported on each record; defaults to `CustomEnzyme`. |
| `topology` | string | No | `Linear` (default) or `Circular`. |
| `second_cut_forward` | integer \| null | No | Second top-strand cut offset of a Type IIB cutter. |
| `second_cut_reverse` | integer \| null | No | Second bottom-strand cut offset; same stagger as the first pair. |
| `methylation` | string \| null | No | `none` (default), `dam`, `dcm`, `dam+dcm` or `emboss`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites` | object[] | Site records in the same shape as `find_restriction_sites` (`position`, `enzyme`, `isForwardStrand`, `cutPosition`, `recognizedSequence`, `isSecondCut`). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Pass both second_cut_forward and second_cut_reverse, or neither. |
| 1002 | Enzyme '&lt;name&gt;' recognition sequence '&lt;site&gt;' contains '&lt;char&gt;' at index &lt;i&gt;, which is not one of the 15 IUPAC nucleotide codes |
| 1003 | Enzyme '&lt;name&gt;' second cut pair (&lt;f&gt;, &lt;r&gt;) leaves a different stagger from its first pair |
| 1004 | Sequence cannot be null or empty. |
| 1005 | Unknown methylation set '{name}'; expected none, dam, dcm, dam+dcm or emboss. |
| 1006 | The topology must be MoleculeTopology.Linear or MoleculeTopology.Circular. |

## Examples

### Example 1: EcoRI by hand — `AAAGAATTCAAA`, site `GAATTC`, cuts `(1, 5)` → the same 2 records `find_restriction_sites` gives for `EcoRI`: position 3, forward cut 4, reverse cut 8.

### Example 2: a Type IIS distal cutter — `AAAAAAAAAAGGTCTCAAAAAAAAAA`, site `GGTCTC`, cuts `(7, 11)` → one forward record at position 10 with cut 17.

### Example 3: non-palindromic, both strands — `AAAGGTCTCAAAGAGACCAAA` with the same definition → one forward record and one reverse record, the latter from the `GAGACC` window.

### Example 4: IUPAC in the sequence — `NNNGAATTCNNN` with EcoRI's definition → 2 records at position 3. `DnaSequence` would have rejected this sequence, which is the whole point of B08 H27.

### Example 5: `second_cut_forward` without `second_cut_reverse` → **error** 1001.

## See Also

- [find_restriction_sites](find_restriction_sites.md) - For any of the 754 named REBASE enzymes.
- [restriction_digest_custom_enzyme](restriction_digest_custom_enzyme.md) - Fragments rather than sites.
- [get_enzyme](get_enzyme.md) - Check whether the enzyme is resolvable by name first.
