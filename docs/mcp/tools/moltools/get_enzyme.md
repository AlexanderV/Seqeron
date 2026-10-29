# get_enzyme

Look up a restriction enzyme by name (any of the 754 REBASE enzymes with a defined cut).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `get_enzyme` |
| **Method ID** | `RestrictionAnalyzer.GetEnzyme` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Case-insensitive lookup of a single enzyme among the **754 REBASE enzymes with a defined cut** (Biopython 1.88 REBASE 404 snapshot, vendored by B08 D4 / F59 — the 39 common laboratory enzymes included, as are Type IIS cutters such as `BsaI` and two-sided Type IIB cutters such as `BaeI`), or `enzyme = null` when the name is unknown **or** is one of the 334 REBASE enzymes whose cut position REBASE does not define (e.g. `Dpi3084I`) — such an enzyme cannot be digested.

The `enzyme` record carries the recognition sequence, the forward/reverse cut positions, the source
organism and three derived fields: `recognitionLength`, `isBluntEnd` and `overhangType`
(`Blunt` / `FivePrime` / `ThreePrime`). A two-sided (Type IIB) cutter also carries `secondCutForward` /
`secondCutReverse`, its second cut pair in the same convention (Biopython `scd5`, `size + scd3`;
`BaeI` → 23 / 18); the SDK omits both for a one-cut enzyme (they are `null`).

Alongside it the result reports:

- **`overhangSequence`** — the single-stranded overhang the enzyme is *defined* to leave: the
  top-strand bases between its two strand cuts, written 5′→3′. This is REBASE / Biopython
  `ovhgseq`: `EcoRI` → `"AATT"`, `PstI` → `"TGCA"`, `NdeI` → `"TA"`, `""` for a blunt cutter, and an
  IUPAC code where the enzyme's own site is degenerate inside the window (`SfiI` → `"NNN"` — the
  overhang of a divided palindrome genuinely is three unspecified bases; the bases a *particular*
  site leaves come from `restriction_digest`'s per-fragment overhang fields).
- **`equischizomers`** — other enzymes of `enzyme_set` (default `"common"` = the 39 common
  enzymes; `"rebase"` = all 754) with the same recognition sequence **and** the same cut
  positions, so they are interchangeable in a digest (`MspI` ↔ `HpaII`, `MboI` ↔ `Sau3AI`).
- **`neoschizomers`** — other enzymes of `enzyme_set` with the same recognition sequence but a different cut, so
  they leave different ends (`MboI` → `DpnI`).

- **`rebaseMethylationField`** / **`cognateMethylationSites`** — REBASE `allenz.104` field `<6>`
  ("METHYLATION SITE") verbatim, and its typed parse: one entry per modified base with `position`
  (1-based within the recognition sequence, counted from the 5′ end of the named strand; `null`
  for REBASE `?`), `isComplementaryStrand` (REBASE negative position) and `type`
  (`N4Methylcytosine` / `C5Methylcytosine` / `N6Methyladenine` = REBASE `(4)`/`(5)`/`(6)`).
  `EcoRI` → `"3(6)"` (m6A at the 3rd base of `GAATTC`), `MboI` → `"2(6)"`, `SfiI` → `"?(4)"`,
  `DpnI` → `""`. This is the base modified by the enzyme's **cognate** methyltransferase (the
  protection half of its R-M system); it is **not** a statement about whether Dam, Dcm or CpG
  methylation blocks the enzyme — a per-enzyme sensitivity table is not modelled (only REBASE/NEB carry one). The separate opt-in `methylation` parameter of `find_restriction_sites`, `find_all_restriction_sites`, `restriction_digest`, `partial_digest`, `digest_summary` and `restriction_map` applies EMBOSS's window rule driven by the *methylase* sites of the chosen set, not this field. All 39 values are transcribed from the
  file programmatically and locked by a test (B08 D7 / F52).

`recognitionLength` is the length of the full recognition string, i.e. Biopython `.size`, so the
interrupted palindrome `SfiI` (`GGCCNNNNNGGCC`) reports **13** although only 8 of its bases are
specified.

All 39 common entries were audited against Biopython 1.88 on 2026-10-14 — recognition string,
size, both cut offsets, blunt/sticky, 5′/3′ polarity and overhang sequence — with **0 wrong
entries**.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L98](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L98)
- Overhang (`GetOverhangSequence`): [RestrictionAnalyzer.cs#L1684](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1684)
- Isoschizomers (`FindIsoschizomers`): [RestrictionAnalyzer.cs#L1939](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1939)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `name` | string | Yes | Enzyme name (non-blank); resolved against all 754 REBASE enzymes with a defined cut. |
| `enzyme_set` | string | No | Table the `equischizomers` / `neoschizomers` lists are drawn from: `"common"` (default, the 39) or `"rebase"` (all 754). Case-insensitive. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `enzyme` | object \| null | The enzyme record (`name`, `recognitionSequence`, `cutPositionForward`, `cutPositionReverse`, `organism`, `rebaseMethylationField`, `cognateMethylationSites`, `secondCutForward` / `secondCutReverse` for a two-cut enzyme, `isBluntEnd`, `overhangType`, `recognitionLength`), or `null` if not found. |
| `overhangSequence` | string \| null | The enzyme-level overhang (REBASE `ovhgseq`); `""` for a blunt cutter, `null` when `enzyme` is `null`. |
| `equischizomers` | array \| null | Names of `enzyme_set` enzymes with the same site **and** cut positions; `null` when `enzyme` is `null`. |
| `neoschizomers` | array \| null | Names of `enzyme_set` enzymes with the same site but a different cut; `null` when `enzyme` is `null`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Enzyme name cannot be null or blank |
| 1010 | Unknown enzyme set '&lt;value&gt;'; expected "common" or "rebase" (`enzyme_set`) |

## Examples

### Example 1: EcoRI

`get_enzyme("EcoRI")` → `{ enzyme: { name: "EcoRI", recognitionSequence: "GAATTC", cutPositionForward: 1, cutPositionReverse: 5, organism: "Escherichia coli", isBluntEnd: false, overhangType: "FivePrime", recognitionLength: 6, rebaseMethylationField: "3(6)", cognateMethylationSites: [{ position: 3, isComplementaryStrand: false, type: "N6Methyladenine" }] }, overhangSequence: "AATT", equischizomers: [], neoschizomers: [] }`. Lookup is case-insensitive (`ecori` works too). Biopython: `EcoRI.ovhg == -4`, `EcoRI.ovhgseq == "AATT"`, `EcoRI.elucidate() == "G^AATT_C"`.

### Example 2: MboI — same site as two other common enzymes

`get_enzyme("MboI")` → `overhangSequence: "GATC"`, `equischizomers: ["Sau3AI"]` (both `GATC` 0/4),
`neoschizomers: ["DpnI"]` (`GATC` 2/2, blunt). Matches Biopython
`MboI.equischizomers(batch)` / `.neoschizomers(batch)` restricted to these 39 enzymes.

### Example 3: SfiI — a divided palindrome

`get_enzyme("SfiI")` → `recognitionLength: 13`, `overhangType: "ThreePrime"`,
`overhangSequence: "NNN"` (Biopython `SfiI.ovhg == 3`, `SfiI.ovhgseq == "NNN"`,
`SfiI.elucidate() == "GGCCN_NNN^NGGCC"`).

### Example 4: Unknown enzyme → `{ enzyme: null, overhangSequence: null, equischizomers: null, neoschizomers: null }`.

### Example 5: BsaI over the REBASE set

`get_enzyme("BsaI", enzyme_set: "rebase")` → `enzyme: { name: "BsaI", recognitionSequence: "GGTCTC",
cutPositionForward: 7, cutPositionReverse: 11, organism: "Bacillus stearothermophilus 6-55",
overhangType: "FivePrime", … }`, `overhangSequence: "NNNN"`,
`equischizomers: ["Bso31I", "BspTNI", "Eco31I"]`, `neoschizomers: []` — Biopython 1.88
`BsaI.equischizomers()` over the 754. With the default `enzyme_set` both lists are empty (no BsaI
isoschizomer among the 39).

### Example 6: BaeI — a two-sided (Type IIB) cutter

`get_enzyme("BaeI")` → `recognitionSequence: "ACNNNNGTAYC"`, `cutPositionForward: -10`,
`cutPositionReverse: -15`, `secondCutForward: 23`, `secondCutReverse: 18`,
`overhangType: "ThreePrime"`, `overhangSequence: "NNNNN"` (Biopython `fst5 = -10`,
`size + fst3 = 11 - 26`, `scd5 = 23`, `size + scd3 = 11 + 7`).

### Example 7: Unknown enzyme set

`get_enzyme("EcoRI", enzyme_set: "all")` → error 1010
`Unknown enzyme set 'all'; expected "common" or "rebase".`

## See Also

- [enzymes_by_cut_length](enzymes_by_cut_length.md), [find_restriction_sites](find_restriction_sites.md)
- [enzymes_compatible](enzymes_compatible.md) - Whether two enzymes' ends can be ligated.
- [compatible_enzymes](compatible_enzymes.md) - What a given enzyme's end can be ligated to.
