# compatible_enzymes

Enumerate pairs of enzymes that produce ligatable ends (the 39 common enzymes, or all 754 REBASE enzymes, paged), or one enzyme's compatible set.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `compatible_enzymes` |
| **Method ID** | `RestrictionAnalyzer.FindCompatibleEnzymes` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Without `enzyme_name`, returns every unordered pair of **distinct** enzymes of `enzyme_set` whose ends can
be ligated together (default `common`, the 39 enzymes: 64 pairs; `rebase`, the 754 enzymes with a defined
cut: **27 869** pairs, in table order). Each result carries the shared end (`"blunt"` or the overhang
string). The response is a bounded page: at most `max_pairs` pairs (default 1000, maximum 10 000) from
`offset`, plus `totalPairs` and `truncated`; the rebase listing takes 28 pages of 1000.

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

- Source: `RestrictionAnalyzer.FindCompatibleEnzymes()` / `FindCompatibleEnzymes(EnzymeSet)` ([RestrictionAnalyzer.cs#L1596](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1596), [RestrictionAnalyzer.Rebase.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.Rebase.cs))
- Per-enzyme overload: [RestrictionAnalyzer.cs#L1631](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1631)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `enzyme_name` | string \| null | No | Restrict the result to pairs involving this enzyme (case-insensitive, itself included when it is in `enzyme_set`; any of the 754 REBASE enzymes with a defined cut). Omit or leave blank for every pair of `enzyme_set`. |
| `enzyme_set` | string | No | Table to list: `"common"` (default, the 39) or `"rebase"` (all 754). Case-insensitive. |
| `max_pairs` | integer \| null | No | All-pairs listing only: maximum pairs returned, 1-10000 (default 1000, which holds the whole common set). Rejected together with `enzyme_name`. |
| `offset` | integer \| null | No | All-pairs listing only: pairs to skip before the page, >= 0 (default 0). Rejected together with `enzyme_name`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pairs` | array | Compatible enzyme pairs, each `{enzyme1, enzyme2, compatibleEnd}`. With `enzyme_name` set, `enzyme1` is that enzyme (as given) in every row. |
| `totalPairs` | integer | All-pairs listing only: number of compatible pairs in the whole `enzyme_set`, independent of the page (common 64, rebase 27 869). Omitted with `enzyme_name`. |
| `truncated` | boolean | All-pairs listing only: `true` when more pairs follow the returned page. Omitted with `enzyme_name`. |

## Errors

| Code | Message |
|------|---------|
| 1003 | Unknown enzyme: &lt;name&gt; (only when `enzyme_name` is given) |
| 1010 | Unknown enzyme set '&lt;value&gt;'; expected "common" or "rebase" (`enzyme_set`) |
| 1003 | `max_pairs` must be between 1 and 10000; `offset` cannot be negative (all-pairs listing) |
| 1011 | `max_pairs` and `offset` page the all-pairs listing; they cannot be combined with `enzyme_name` |

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

### Example 5: Partners over the REBASE set; the REBASE all-pairs listing

`compatible_enzymes(enzyme_name: "BamHI", enzyme_set: "rebase")` → **59** rows, exactly Biopython 1.88
`BamHI.compatible_end(RestrictionBatch(754))`: the `GATC` 5′ cutters (`BclI`, `BglII`, `BstYI`,
`MboI`, `Sau3AI`, `XhoII`, …) plus the distal Type IIS cutters whose overhang is unspecified (`BsaI`,
`BsmBI`, `FokI`, … — Biopython's `Ambiguous` class, F57), each with `compatibleEnd = "GATC"`.
`EcoRI` → 51 rows (`AATT`: `ApoI`, `MfeI`, `MunI`, … plus the distal 4-base 5′ cutters).
`compatible_enzymes(enzyme_set: "rebase")` → the first 1000 of **27 869** pairs (`totalPairs = 27869`,
`truncated = true`; first row `AanI`/`Acc16I` `blunt`); `offset: 27865, max_pairs: 10` returns the last 4
(`truncated = false`). `(EcoRI, MfeI, AATT)` is among them. On the 11 200 pairs among the 389 enzymes
with a fully defined overhang the set equals Biopython 1.88 `a % b` exactly (the remaining pairs involve
an ambiguous overhang, where `%` is unsound, and follow the derived rule).


## Performance

- **Time Complexity:** O(n²) over `enzyme_set` (n = 39 or 754; the 754-enzyme scan is 283 881 pair tests, well under a second).
- **Space Complexity:** O(`max_pairs`) — the listing is streamed and only the requested page is kept.

## See Also

- [enzymes_compatible](enzymes_compatible.md) - Test a single named pair.
- [blunt_cutters](blunt_cutters.md), [sticky_cutters](sticky_cutters.md)
