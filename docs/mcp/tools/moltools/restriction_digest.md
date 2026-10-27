# restriction_digest

Simulate a restriction digest of a linear or circular molecule.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `restriction_digest` |
| **Method ID** | `RestrictionAnalyzer.Digest` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Cuts the molecule at the top-strand (forward-strand) cut positions of the named enzymes and returns the resulting fragments in ascending `startPosition`. Each fragment carries its sequence, start position, length, flanking enzymes (`leftEnzyme`/`rightEnzyme`, and the complete lists `leftEnzymes`/`rightEnzymes`), fragment number, and the chemistry of each end created by a cut (`leftOverhangType`/`leftOverhang`, `rightOverhangType`/`rightOverhang`).

**Topology.** With `topology=Linear` (the default) `k` distinct cut positions yield `k+1` fragments; with zero cuts the whole molecule is returned as one fragment, and its two ends are the original termini (overhang fields `null`). With `topology=Circular` the site search wraps past the origin — so an origin-spanning recognition site is found and cut — and `k` cut positions yield exactly `k` fragments. The origin-spanning fragment is reported **last** (its `startPosition` is the largest cut), and every circular fragment end is a cut, so both flanks carry an enzyme and an overhang. A circular molecule with zero cuts stays whole (one full-length uncut fragment, no ends); a single cutter linearizes it into one full-length fragment.

**Fragment sequences are top-strand bases.** A digest partitions the top strand at the top-strand cuts, the same convention as Biopython `Bio.Restriction` `enzyme.catalyse(seq[, linear=False])`, so a sticky end's bases belong to the fragment downstream of the top-strand cut and the single-stranded extension is reported separately in the overhang fields (read 5′→3′ on the strand that carries it, taken from the molecule). Blunt cutters report `Blunt` with an empty overhang. Biopython emits the origin-spanning fragment first rather than last, so the circular fragment order differs from Biopython's by a rotation of one.

**Methylation (opt-in).** `methylation` (`none`|`dam`|`dcm`|`dam+dcm`|`emboss`) masks the bases the chosen methylase set methylates (both strands; the `N` spacer of EcoKI included) and drops every enzyme window containing one — the semantics of EMBOSS `restrict -methylation -mfile Emethylsites.dat`, confirmed against the binary (RESTR-FIND-001 Evidence §7). Default `none` leaves the output unchanged. Example: `GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA` with `MboI` and `methylation=dam` reports no site (without it, positions 10 and 36); EcoRI at 0 and 18 is unaffected.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L1124](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1124)

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
| `fragments` | array | Fragments, each `{sequence, startPosition, length, leftEnzyme, rightEnzyme, fragmentNumber, leftOverhangType, leftOverhang, rightOverhangType, rightOverhang, leftEnzymes, rightEnzymes}`. |

**Coinciding cuts.** Several enzymes can cut the top strand at the same position — isoschizomers (MspI/HpaII, MboI/Sau3AI) and nested sites (BamHI or BglII with MboI/Sau3AI, ClaI/TaqI, StuI/HaeIII, ScaI/RsaI, FseI/HaeIII; Biopython 1.88 `search` reports the same position for each). Such a cut is one boundary. `leftEnzymes`/`rightEnzymes` list every enzyme cutting at that end, once each, in ordinal order of the enzyme-table names (empty for an uncut end); `leftEnzyme`/`rightEnzyme` keep their meaning as the first name of that list, so the result no longer depends on the order of `enzyme_names`. The overhang fields give that first enzyme's end. All built-in collisions but one join identical ends; for FseI/HaeIII (FseI 3′ `CCGG`, HaeIII blunt) the reported end is FseI's alone, and the combined digest's physical end is not modelled.

`leftOverhangType`/`rightOverhangType` are `Blunt`, `FivePrime` or `ThreePrime`, or `null` when that end was not created by a cut. `leftOverhang`/`rightOverhang` give the extension's bases (empty string for a blunt end, `null` for an uncut end).

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | At least one enzyme name is required |

## Examples

### Example 1: EcoRI on linear `AAAGAATTCAAA`

Cut at position 4 → `AAAG` (0–4) and `AATTCAAA` (4–12); the cut ends are 5′ overhangs `AATT`.

### Example 2: EcoRI on circular `TTCAAAAAAAGAA`

The `GAATTC` site spans the origin (starting at position 10), so the plasmid is linearized into one 13 bp fragment `AATTCAAAAAAAG` starting at the cut, position 11 — matching `Bio.Restriction.EcoRI.catalyse(Seq("TTCAAAAAAAGAA"), linear=False) == (Seq("AATTCAAAAAAAG"),)`. With `topology=Linear` the same sequence is uncut.

### Example 3: No cut site → single whole-sequence fragment (linear or circular).

### Example 4: MspI + HpaII on linear `CCGGAAACCGG`

The isoschizomers cut at the same two positions (Biopython 1.88 `MspI.search == HpaII.search == [2, 9]`), so there are two boundaries and three fragments `C`, `CGGAAAC`, `CGG`. Each cut end lists `["HpaII", "MspI"]` in `leftEnzymes`/`rightEnzymes`, reports `HpaII` (ordinally first) as `leftEnzyme`/`rightEnzyme`, and carries the shared 5′ overhang `CG` — the same output whichever order the two names are given in.

## See Also

- [digest_summary](digest_summary.md), [find_restriction_sites](find_restriction_sites.md), [restriction_map](restriction_map.md)
