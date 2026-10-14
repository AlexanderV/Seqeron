# compatible_enzymes

Enumerate pairs of built-in enzymes that produce ligatable ends, or one enzyme's compatible set.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `compatible_enzymes` |
| **Method ID** | `RestrictionAnalyzer.FindCompatibleEnzymes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Without arguments, returns every unordered pair of **distinct** built-in restriction enzymes whose
ends can be ligated together. Each result carries the shared end (`"blunt"` or the overhang string).

Two ends are compatible when **both are blunt**, or when they have the **same polarity** (both 5′ or
both 3′), the same overhang length, and one overhang is the reverse complement of the other (for an
overhang carrying IUPAC codes, matching is per-position code intersection). Identical overhang
letters with *opposite* polarity are **not** compatible — see `enzymes_compatible` for the full rule.

**With `enzyme_name` set**, returns what that one enzyme can be ligated to, **including itself**.
This is the query to use for an enzyme with a unique overhang: 15 of the 39 built-ins (`EcoRI`,
`PstI`, `KpnI`, `SacI`, `SphI`, `NcoI`, `NotI`, `ApaI`, `AscI`, `AgeI`, `FseI`, `NdeI`, `PacI`,
`SfiI`, `HindIII`) are compatible only with themselves and therefore appear in **no** pair of the
unfiltered listing. It is the per-enzyme form of Biopython's `compatible_end()`, without that
method's own bug of filtering candidates to 5′-overhang enzymes only (which makes it answer the
empty list for every blunt and every 3′ cutter).

Verified against Biopython 1.88 on 2026-10-14: 1521/1521 ordered built-in pairs and all 39
per-enzyme sets agree with the symmetrised `%` operator.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L1206](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1206)
- Per-enzyme overload: [RestrictionAnalyzer.cs#L1240](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1240)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `enzyme_name` | string \| null | No | Restrict the result to pairs involving this enzyme (case-insensitive, itself included). Omit or leave blank for every pair in the table. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pairs` | array | Compatible enzyme pairs, each `{enzyme1, enzyme2, compatibleEnd}`. With `enzyme_name` set, `enzyme1` is that enzyme (as given) in every row. |

## Errors

| Code | Message |
|------|---------|
| 1003 | Unknown enzyme: &lt;name&gt; (only when `enzyme_name` is given) |

## Examples

### Example 1: Sticky-end compatibility

BamHI (`GGATCC`, cut 1/5) and BglII (`AGATCT`, cut 1/5) both leave a 5′ `GATC` overhang, so they appear as a pair with `compatibleEnd = "GATC"`.

### Example 2: Blunt-end compatibility

All 10 blunt cutters are pairwise compatible, contributing C(10,2) = **45** pairs with `compatibleEnd = "blunt"`.

### Example 3: The complete listing

`compatible_enzymes()` returns **64** pairs: 45 `"blunt"`, 6 `"GATC"` (BamHI/BglII/MboI/Sau3AI),
6 `"CTAG"` (NheI/XbaI/SpeI/AvrII), 6 `"CG"` (ClaI/TaqI/MspI/HpaII) and 1 `"TCGA"` (SalI/XhoI). No
other shared end exists in this table.

### Example 4: Per-enzyme query

`compatible_enzymes(enzyme_name: "BamHI")` → 4 rows, `compatibleEnd = "GATC"`: `BamHI`, `BglII`,
`MboI`, `Sau3AI`. `compatible_enzymes(enzyme_name: "SfiI")` → 1 row, `SfiI` itself with
`compatibleEnd = "NNN"`. `compatible_enzymes(enzyme_name: "EcoRV")` → the 10 blunt cutters.

## Performance

- **Time Complexity:** O(n²) over the enzyme database (n = 39).
- **Space Complexity:** O(number of compatible pairs).

## See Also

- [enzymes_compatible](enzymes_compatible.md) - Test a single named pair.
- [blunt_cutters](blunt_cutters.md), [sticky_cutters](sticky_cutters.md)
