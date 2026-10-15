# get_enzyme

Look up a built-in restriction enzyme by name.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `get_enzyme` |
| **Method ID** | `RestrictionAnalyzer.GetEnzyme` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Case-insensitive lookup of a single enzyme in the built-in Type II restriction-enzyme database, or `enzyme = null` when the name is not present.

The `enzyme` record carries the recognition sequence, the forward/reverse cut positions, the source
organism and three derived fields: `recognitionLength`, `isBluntEnd` and `overhangType`
(`Blunt` / `FivePrime` / `ThreePrime`).

Alongside it the result reports:

- **`overhangSequence`** — the single-stranded overhang the enzyme is *defined* to leave: the
  top-strand bases between its two strand cuts, written 5′→3′. This is REBASE / Biopython
  `ovhgseq`: `EcoRI` → `"AATT"`, `PstI` → `"TGCA"`, `NdeI` → `"TA"`, `""` for a blunt cutter, and an
  IUPAC code where the enzyme's own site is degenerate inside the window (`SfiI` → `"NNN"` — the
  overhang of a divided palindrome genuinely is three unspecified bases; the bases a *particular*
  site leaves come from `restriction_digest`'s per-fragment overhang fields).
- **`equischizomers`** — other built-ins with the same recognition sequence **and** the same cut
  positions, so they are interchangeable in a digest (`MspI` ↔ `HpaII`, `MboI` ↔ `Sau3AI`).
- **`neoschizomers`** — other built-ins with the same recognition sequence but a different cut, so
  they leave different ends (`MboI` → `DpnI`).

- **`rebaseMethylationField`** / **`cognateMethylationSites`** — REBASE `allenz.104` field `<6>`
  ("METHYLATION SITE") verbatim, and its typed parse: one entry per modified base with `position`
  (1-based within the recognition sequence, counted from the 5′ end of the named strand; `null`
  for REBASE `?`), `isComplementaryStrand` (REBASE negative position) and `type`
  (`N4Methylcytosine` / `C5Methylcytosine` / `N6Methyladenine` = REBASE `(4)`/`(5)`/`(6)`).
  `EcoRI` → `"3(6)"` (m6A at the 3rd base of `GAATTC`), `MboI` → `"2(6)"`, `SfiI` → `"?(4)"`,
  `DpnI` → `""`. This is the base modified by the enzyme's **cognate** methyltransferase (the
  protection half of its R-M system); it is **not** a statement about whether Dam, Dcm or CpG
  methylation blocks the enzyme — that is not modelled. All 39 values are transcribed from the
  file programmatically and locked by a test (B08 D7 / F52).

`recognitionLength` is the length of the full recognition string, i.e. Biopython `.size`, so the
interrupted palindrome `SfiI` (`GGCCNNNNNGGCC`) reports **13** although only 8 of its bases are
specified.

All 39 built-in entries were audited against Biopython 1.88 on 2026-10-14 — recognition string,
size, both cut offsets, blunt/sticky, 5′/3′ polarity and overhang sequence — with **0 wrong
entries**.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L79](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L79)
- Overhang: [RestrictionAnalyzer.cs#L1153](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1153)
- Isoschizomers: [RestrictionAnalyzer.cs#L1413](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1413)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `name` | string | Yes | Enzyme name (non-blank). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `enzyme` | object \| null | The enzyme record (`name`, `recognitionSequence`, `cutPositionForward`, `cutPositionReverse`, `organism`, `rebaseMethylationField`, `cognateMethylationSites`, `isBluntEnd`, `overhangType`, `recognitionLength`), or `null` if not found. |
| `overhangSequence` | string \| null | The enzyme-level overhang (REBASE `ovhgseq`); `""` for a blunt cutter, `null` when `enzyme` is `null`. |
| `equischizomers` | array \| null | Names of built-ins with the same site **and** cut positions; `null` when `enzyme` is `null`. |
| `neoschizomers` | array \| null | Names of built-ins with the same site but a different cut; `null` when `enzyme` is `null`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Enzyme name cannot be null or blank |

## Examples

### Example 1: EcoRI

`get_enzyme("EcoRI")` → `{ enzyme: { name: "EcoRI", recognitionSequence: "GAATTC", cutPositionForward: 1, cutPositionReverse: 5, organism: "Escherichia coli", isBluntEnd: false, overhangType: "FivePrime", recognitionLength: 6, rebaseMethylationField: "3(6)", cognateMethylationSites: [{ position: 3, isComplementaryStrand: false, type: "N6Methyladenine" }] }, overhangSequence: "AATT", equischizomers: [], neoschizomers: [] }`. Lookup is case-insensitive (`ecori` works too). Biopython: `EcoRI.ovhg == -4`, `EcoRI.ovhgseq == "AATT"`, `EcoRI.elucidate() == "G^AATT_C"`.

### Example 2: MboI — same site as two other built-ins

`get_enzyme("MboI")` → `overhangSequence: "GATC"`, `equischizomers: ["Sau3AI"]` (both `GATC` 0/4),
`neoschizomers: ["DpnI"]` (`GATC` 2/2, blunt). Matches Biopython
`MboI.equischizomers(batch)` / `.neoschizomers(batch)` restricted to these 39 enzymes.

### Example 3: SfiI — a divided palindrome

`get_enzyme("SfiI")` → `recognitionLength: 13`, `overhangType: "ThreePrime"`,
`overhangSequence: "NNN"` (Biopython `SfiI.ovhg == 3`, `SfiI.ovhgseq == "NNN"`,
`SfiI.elucidate() == "GGCCN_NNN^NGGCC"`).

### Example 4: Unknown enzyme → `{ enzyme: null, overhangSequence: null, equischizomers: null, neoschizomers: null }`.

## See Also

- [enzymes_by_cut_length](enzymes_by_cut_length.md), [find_restriction_sites](find_restriction_sites.md)
- [enzymes_compatible](enzymes_compatible.md) - Whether two enzymes' ends can be ligated.
- [compatible_enzymes](compatible_enzymes.md) - What a given enzyme's end can be ligated to.
