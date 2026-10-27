# enzymes_by_cut_length

List restriction enzymes by recognition-sequence length (common set by default, or all 754 REBASE enzymes).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `enzymes_by_cut_length` |
| **Method ID** | `RestrictionAnalyzer.GetEnzymesByCutLength` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns every enzyme of `enzyme_set` — `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut) — whose recognition sequence is exactly `length` base pairs. Short recognition sites (4 bp) cut frequently, 6 bp are typical cloning enzymes, and 8 bp are rare cutters. Lengths with no matching enzyme return an empty list.

The length measured is the **full recognition string**, i.e. Biopython `.size` — so an interrupted
palindrome counts its degenerate spacer: `SfiI` (`GGCCNNNNNGGCC`) is returned by `length = 13` and
**not** by `length = 8`, although only 8 of its bases are specified and catalogues often call it an
8-base-specificity enzyme.

The common set contains **9** 4-cutters (MspI, HpaII, TaqI, AluI, RsaI, HaeIII, DpnI,
Sau3AI, MboI), **24** 6-cutters, **5** 8-cutters (NotI, PacI, AscI, FseI, SwaI) and the one 13-bp
interrupted palindrome SfiI — 39 in all. Every `recognitionLength` was checked against Biopython
`.size` on 2026-10-14: 39/39 agreement.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L127](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L127)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `length` | integer | Yes | Recognition-sequence length in bp (positive). |
| `enzyme_set` | string | No | `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut). Case-insensitive. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `enzymes` | array | Enzymes with recognition length == `length`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Recognition-sequence length must be positive |
| 1010 | Unknown enzyme set '&lt;value&gt;'; expected "common" or "rebase" (`enzyme_set`) |

## Examples

### Example 1: `length = 6` → 24 six-cutters (EcoRI, BamHI, HindIII, …).

### Example 2: `length = 8` → 5 rare cutters (NotI, PacI, AscI, FseI, SwaI). SfiI is **not** among them — ask for `length = 13`.

### Example 3: `length = 13` → the single interrupted palindrome SfiI (`GGCCNNNNNGGCC`, cut 8/5, 3′ `NNN` overhang).

### Example: the REBASE set

`enzymes_by_cut_length(length: 6, enzyme_set: "rebase")` → **385** enzymes (Biopython 1.88: `.size == 6`
over the 754). `length: 8` → 28 (`AbsI`, `AscI`, `AsiSI`, … `SwaI`); `length: 13` → 6: `AloI`,
`ArsI`, `BarI`, `PcsI`, `PsrI`, `SfiI`. With the default set, `length: 6` stays **24**.

### Example: unknown enzyme set

`enzymes_by_cut_length(length: 6, enzyme_set: "all")` → error 1010.

## See Also

- [get_enzyme](get_enzyme.md), [blunt_cutters](blunt_cutters.md), [sticky_cutters](sticky_cutters.md)
