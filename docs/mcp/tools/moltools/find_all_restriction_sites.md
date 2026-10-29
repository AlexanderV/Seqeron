# find_all_restriction_sites

Find sites for every enzyme of a set (the 39 common enzymes by default, or all 754 REBASE enzymes).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `find_all_restriction_sites` |
| **Method ID** | `RestrictionAnalyzer.FindAllSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Performs a comprehensive restriction scan: every enzyme of `enzyme_set` — `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut) — is searched against both strands of the sequence. Use this when the caller does not know which enzymes to look for. Coordinates, the two records reported per palindromic site, overlapping-site reporting, the `Linear` cut-feasibility rule and the `Circular` origin-spanning search are all exactly as documented for [find_restriction_sites](find_restriction_sites.md).

**Methylation (opt-in).** `methylation` (`none`|`dam`|`dcm`|`dam+dcm`|`emboss`) masks the bases the chosen methylase set methylates (both strands; the `N` spacer of EcoKI included) once for the molecule and drops every enzyme window containing one, for every enzyme of `enzyme_set` — the semantics of EMBOSS `restrict -methylation -mfile Emethylsites.dat`, confirmed against the binary (RESTR-FIND-001 Evidence §7; B08 round 3 G4 / F84). Default `none` leaves the output unchanged. Example: `GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA` with `methylation=dam` keeps EcoRI (0, 18), HpaII/MspI (28), SmaI (27) and TaqI (38, 44) and loses DpnI, MboI and Sau3AI at 10 and 36 (EMBOSS over the 39 common enzymes gives exactly this).

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L547](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L547)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to scan (non-empty). |
| `topology` | string | No | `Linear` (default) or `Circular`. |
| `enzyme_set` | string | No | `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut). Case-insensitive. |
| `methylation` | string | No | `none` (default), `dam` (`GmATC`), `dcm` (`mCmCAGG`, `mCmCTGG`), `dam+dcm` or `emboss` (all five `Emethylsites.dat` entries incl. EcoKI). Case-insensitive; an unknown value is an error. A recognition window containing a methylated base is not a site (EMBOSS `restrict -methylation`; a window rule, not an enzyme-sensitivity table — BamHI is blocked by Dam). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites` | array | All restriction sites found across every enzyme of the set, each `{position, enzyme, isForwardStrand, cutPosition, recognizedSequence, isSecondCut}`: `position` is the 0-based forward-strand start of the window, `enzyme` the enzyme record (as in `enzymes_by_cut_length`), `isForwardStrand` whether the window reads as the site on the forward strand, `cutPosition` the top-strand cut of a forward record or the bottom-strand cut of a reverse record (reduced modulo the length on a circle), `recognizedSequence` the window on the strand that carries the site (upper-case). `isSecondCut` is `true` on the second record of a two-sided (Type IIB) enzyme such as BaeI or BcgI, whose one site is cleaved on both sides: it is the only way to tell the two records of one site apart (`false` for every other record). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1010 | Unknown enzyme set '&lt;value&gt;'; expected "common" or "rebase" |
| 1011 | Unknown methylation set '{name}'; expected none, dam, dcm, dam+dcm or emboss. |
| 1012 | The topology must be MoleculeTopology.Linear or MoleculeTopology.Circular. |

## Examples

### Example 1: `AAAGAATTCAAA` → includes the EcoRI site at position 3 (forward cut 4, plus reverse strand).

### Example 2: `AAAAAAAA` → `sites: []` (no recognition sequence matches).

### Example 3: `TTCAAAAAAAGAA` + `topology: Circular` → includes the origin-spanning EcoRI site at position 10 (forward cut 11), which the linear scan does not report.

### Example: the REBASE set

`find_all_restriction_sites(sequence: "AAGGTCTCAAAAAAAAAAAAAAAGAGACCAAAAAAAAAAA", enzyme_set: "rebase")`
→ 42 site records from 14 enzymes (`AbaSI`, `Alw26I`, `BcoDI`, `BmeDI`, `BsaI`, `BsmAI`, `Bso31I`,
`BspTNI`, `BstMAI`, `Eco31I`, `EcoO157SI`, `MspJI`, `SetI`, `YkrI`); BsaI's two records are
`{ position: 2, isForwardStrand: true, cutPosition: 9, recognizedSequence: "GGTCTC" }` and
`{ position: 23, isForwardStrand: false, cutPosition: 22, recognizedSequence: "GGTCTC" }`. The default
(common) set finds **no** site in this sequence. The set of enzymes with a site equals Biopython
1.88 `RestrictionBatch(754).search(...)` (locked in `EnzymeSetToolsTests` on a second sequence, linear
and circular).

### Example: unknown enzyme set

`find_all_restriction_sites(sequence: "GAATTC", enzyme_set: "all")` → error 1010.

## See Also

- [find_restriction_sites](find_restriction_sites.md), [restriction_map](restriction_map.md)
