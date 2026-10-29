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

Scans a DNA sequence on both strands for the recognition sites of the named enzymes — any of the 754 REBASE enzymes with a defined cut (Biopython 1.88 REBASE 404 snapshot; the 39 common enzymes, Type IIS such as `BsaI` and two-sided Type IIB such as `BaeI` included; an unknown name is an error). IUPAC degenerate codes in recognition sequences are matched against ACGT input.

`position` and `cutPosition` are 0-based **forward-strand** coordinates for both strands. `cutPosition` is the number of forward-strand bases 5' of the cut the enzyme makes on that record's own strand, so a palindromic enzyme reports a forward record (top-strand cut) and a reverse record (bottom-strand cut) at the same `position`, and `recognizedSequence` is the matched window as read on that record's strand. Biopython `Bio.Restriction.<enzyme>.search()` reports the same top-strand cut 1-based, i.e. as `cutPosition + 1`. Overlapping sites are all reported.

With `topology = Linear` (the default) a site is reported only when **both** of its strand cuts fall on a phosphodiester bond that exists, i.e. both lie in `[1, length - 1]` — the same filter as Biopython's `AbstractCut._drop`. This excludes only an enzyme that cuts at the edge of its own recognition window (`Sau3AI` / `MboI`, `^GATC`) sitting at the very 5' or 3' end. With `topology = Circular` the recognition window may span the origin, cut positions are reduced modulo the sequence length, and nothing is dropped.

**Methylation (opt-in).** `methylation` (`none`|`dam`|`dcm`|`dam+dcm`|`emboss`) masks the bases the chosen methylase set methylates (both strands; the `N` spacer of EcoKI included) and drops every enzyme window containing one — the semantics of EMBOSS `restrict -methylation -mfile Emethylsites.dat`, confirmed against the binary (RESTR-FIND-001 Evidence §7). Default `none` leaves the output unchanged. Example: `GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA` with `MboI` and `methylation=dam` reports no site (without it, positions 10 and 36); EcoRI at 0 and 18 is unaffected.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L474](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L474)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to scan (non-empty). IUPAC ambiguity codes are accepted anywhere (a site `N` matches any of the 15 codes, every other site code only a concrete base it expands to), lower case too; the returned sequences are upper case. A character outside the 15 codes (`U`, a space, a FASTA header line) matches nothing: such text is not rejected, but it still occupies coordinates and cannot be part of a site. |
| `enzyme_names` | string[] | Yes | Enzyme names to look for (≥1); resolved case-insensitively. |
| `topology` | string | No | `Linear` (default) or `Circular`. |
| `methylation` | string | No | `none` (default), `dam` (`GmATC`), `dcm` (`mCmCAGG`, `mCmCTGG`), `dam+dcm` or `emboss` (all five `Emethylsites.dat` entries incl. EcoKI). Case-insensitive; an unknown value is an error. A recognition window containing a methylated base is not a site (EMBOSS `restrict -methylation`; a window rule, not an enzyme-sensitivity table — BamHI is blocked by Dam). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites` | array | Site records, each `{position, enzyme, isForwardStrand, cutPosition, recognizedSequence, isSecondCut}`: `position` is the 0-based forward-strand start of the window, `enzyme` the enzyme record (as in `enzymes_by_cut_length`), `isForwardStrand` whether the window reads as the site on the forward strand, `cutPosition` the top-strand cut of a forward record or the bottom-strand cut of a reverse record (reduced modulo the length on a circle), `recognizedSequence` the window on the strand that carries the site (upper-case). `isSecondCut` is `true` on the second record of a two-sided (Type IIB) enzyme such as BaeI or BcgI, whose one site is cleaved on both sides: it is the only way to tell the two records of one site apart (`false` for every other record). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | At least one enzyme name is required |
| 4001 | Unknown enzyme: {name} |
| 4002 | Unknown methylation set '{name}'; expected none, dam, dcm, dam+dcm or emboss. |
| 4003 | The topology must be MoleculeTopology.Linear or MoleculeTopology.Circular. |

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
