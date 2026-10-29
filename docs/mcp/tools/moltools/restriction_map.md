# restriction_map

Build a restriction map of a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `restriction_map` |
| **Method ID** | `RestrictionAnalyzer.CreateMap` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Maps a DNA sequence: reports all forward+reverse restriction sites, the sites grouped by enzyme (positions), the total number of **cleavages**, the **unique cutters** (enzymes that cleave the molecule exactly once — ideal for cloning), and the **non-cutters**.

`totalSites` and `uniqueCutters` count cleavages, not recognition windows: every forward-strand
record counts, plus every reverse-strand record whose backwards reading cuts other bonds. So a
palindromic site with symmetric cuts (all 39 common enzymes) counts once, but a non-palindromic
enzyme's reverse-strand sites count on their own — a Type IIS enzyme such as `BsaI` is **not** a
unique cutter merely because it has one forward-strand site, since `GGTCTC` read on the other strand
is a second cleavage — and the five enzymes whose two strand readings of one window cut different
bonds (`HauII`, `SgeI`, `MspJI`, `LpnPI`, `SgrTI`) count twice for that one window. This paragraph
used to say "forward-strand site count", which read wrong for every Golden Gate enzyme
(review 2026-09 B08 H32; the cut rule itself is F101/H7).

The *considered set* is the enzyme set the map is built over: the names in `enzyme_names`, resolved case-insensitively and de-duplicated, or — when `enzyme_names` is empty — every enzyme of `enzyme_set`: `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut). Named enzymes may be any of the 754; `enzyme_set` is ignored when names are given. `sitesByEnzyme` and `nonCutters` always **partition** that set, and both are keyed by the table's own canonical enzyme name whatever casing the caller used. Site coordinates, the two records per palindromic site and the `Linear`/`Circular` rules are as documented for [find_restriction_sites](find_restriction_sites.md).

**Methylation (opt-in).** `methylation` (`none`|`dam`|`dcm`|`dam+dcm`|`emboss`) masks the bases the chosen methylase set methylates (both strands; the `N` spacer of EcoKI included) and drops every enzyme window containing one — the semantics of EMBOSS `restrict -methylation -mfile Emethylsites.dat`, confirmed against the binary (RESTR-FIND-001 Evidence §7). Default `none` leaves the output unchanged. Example: `GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA` with `MboI` and `methylation=dam` reports no site (without it, positions 10 and 36); EcoRI at 0 and 18 is unaffected.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L1558](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1558)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence to map (non-empty). |
| `enzyme_names` | string[] | Yes | Enzymes to consider (case-insensitive, de-duplicated; any of the 754 REBASE enzymes with a defined cut); empty = every enzyme of `enzyme_set`. |
| `topology` | string | No | `Linear` (default) or `Circular`. |
| `enzyme_set` | string | No | Table used when `enzyme_names` is empty: `"common"` (default, the 39) or `"rebase"` (all 754). Case-insensitive; ignored when names are given. |
| `methylation` | string | No | `none` (default), `dam` (`GmATC`), `dcm` (`mCmCAGG`, `mCmCTGG`), `dam+dcm` or `emboss` (all five `Emethylsites.dat` entries incl. EcoKI). Case-insensitive; an unknown value is an error. A recognition window containing a methylated base is not a site (EMBOSS `restrict -methylation`; a window rule, not an enzyme-sensitivity table — BamHI is blocked by Dam). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sequenceLength` | integer | Sequence length. |
| `sites` | array | All sites (both strands), in the shape of `find_restriction_sites` (including `isSecondCut` for two-sided Type IIB enzymes). |
| `sitesByEnzyme` | object | Enzyme → sorted positions. |
| `totalSites` | integer | Number of cleavages over the considered set (see above). |
| `uniqueCutters` | string[] | Enzymes that cleave the molecule exactly once. |
| `nonCutters` | string[] | Considered enzymes with no site, under their canonical name. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1010 | Unknown enzyme set '&lt;value&gt;'; expected "common" or "rebase" |
| 1011 | Unknown enzyme: {name} |
| 1012 | Unknown methylation set '{name}'; expected none, dam, dcm, dam+dcm or emboss. |
| 1013 | The topology must be MoleculeTopology.Linear or MoleculeTopology.Circular. |

## Examples

### Example 1: `AAAGAATTCAAA` with `[EcoRI, BamHI]`

EcoRI cuts once (unique cutter, positions `[3, 3]` for both strands); BamHI does not cut → non-cutter. `totalSites = 1`. The same holds for `[ecori, bamhi]`: `sitesByEnzyme` is keyed `EcoRI` and `nonCutters` is `["BamHI"]`.

### Example 2: `AAAGAATTCAAA` with `enzyme_names: []`

The considered set is the whole common table, so `sitesByEnzyme` is `{ "EcoRI": [3, 3] }` and `nonCutters` lists the other 38 common enzymes.

### Example: all REBASE enzymes

`restriction_map(sequence: "AAGGTCTCAAAAAAAAAAAAAAAGAGACCAAAAAAAAAAA", enzyme_names: [], enzyme_set: "rebase")`
→ `sitesByEnzyme` has 14 enzymes (e.g. `BsaI: [2, 23]` — its forward site and its reverse-strand
`GAGACC` site), `nonCutters` the other **740** (14 + 740 = 754), `totalSites: 41`,
`uniqueCutters: ["SetI"]`. The same call with the default set maps 0 sites and 39 non-cutters.

### Example: unknown enzyme set

`restriction_map(sequence: "GAATTC", enzyme_names: [], enzyme_set: "all")` → error 1010.

## See Also

- [find_restriction_sites](find_restriction_sites.md), [restriction_digest](restriction_digest.md)
