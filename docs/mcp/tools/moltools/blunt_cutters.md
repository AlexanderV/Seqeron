# blunt_cutters

List the restriction enzymes that produce blunt ends (common set by default, or all 754 REBASE enzymes).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `blunt_cutters` |
| **Method ID** | `RestrictionAnalyzer.GetBluntCutters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns every enzyme of `enzyme_set` — `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut) — that cuts both strands at the same position, leaving a blunt end (no single-stranded overhang). Blunt cutters are used for blunt-end cloning where any two blunt fragments can be ligated.

An enzyme is blunt when `CutPositionForward == CutPositionReverse` — equivalently, when REBASE /
Biopython `ovhg` is 0 and `is_blunt()` is true. Any blunt end ligates to any other blunt end. All
39 built-in entries were checked against `is_blunt()` on 2026-10-14: 39/39 agreement, **10** blunt
cutters.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L177](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L177)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `enzyme_set` | string | No | `"common"` (default, the 39 common laboratory enzymes) or `"rebase"` (all 754 REBASE enzymes with a defined cut). Case-insensitive. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `enzymes` | array | Blunt-cutting enzymes, each with `name`, `recognitionSequence`, `cutPositionForward`, `cutPositionReverse`, `organism`, `rebaseMethylationField` (REBASE allenz.104 field <6>, the cognate-methyltransferase site) with its typed parse `cognateMethylationSites` (see [get_enzyme](get_enzyme.md)), plus the derived `isBluntEnd`, `overhangType` and `recognitionLength`. |

## Errors

| Code | Message |
|------|---------|
| 1010 | Unknown enzyme set '&lt;value&gt;'; expected "common" or "rebase" (`enzyme_set`) |

## Examples

### Example 1: List blunt cutters

**User Prompt:**
> Which enzymes give me blunt ends?

**Expected Tool Call:**
```json
{ "tool": "blunt_cutters", "arguments": {} }
```

**Response (abridged):**
```json
{
  "enzymes": [
    { "name": "AluI", "recognitionSequence": "AGCT", "cutPositionForward": 2, "cutPositionReverse": 2, "organism": "Arthrobacter luteus" },
    { "name": "EcoRV", "recognitionSequence": "GATATC", "cutPositionForward": 3, "cutPositionReverse": 3, "organism": "Escherichia coli" }
  ]
}
```

The common set contains exactly **10** blunt cutters, listed in table order: SwaI, AluI, RsaI, HaeIII,
DpnI, EcoRV, SmaI, HincII, ScaI, StuI. Note `HincII` (`GTYRAC`, cut 3/3): its recognition site is degenerate but
the cut is central, so it is blunt and its overhang is empty (Biopython
`HincII.elucidate() == "GTY^_RAC"`, `ovhg == 0`). `DpnI` is blunt and methylation-dependent
(`GA^_TC`).

### Example 3: the REBASE set

`blunt_cutters(enzyme_set: "rebase")` → **137** enzymes (Biopython 1.88 `is_blunt()` over the 754),
first `{ "name": "AanI", "recognitionSequence": "TTATAA", "cutPositionForward": 3, "cutPositionReverse": 3,
"organism": "Arthrobacter aurescens RFL2", … }`, then `Acc16I`, `AccBSI`, `AccII`, `AcvI`, … The
default call still returns the 10 common blunt cutters.

### Example 4: unknown enzyme set

`blunt_cutters(enzyme_set: "all")` → error 1010.

## Performance

- **Time Complexity:** O(n) over the enzyme database (n = 39).
- **Space Complexity:** O(k) for the k blunt cutters returned.

## See Also

- [sticky_cutters](sticky_cutters.md) - Enzymes producing cohesive (overhang) ends.
- [get_enzyme](get_enzyme.md) - Look up a single enzyme by name, with its overhang sequence.
- [enzymes_compatible](enzymes_compatible.md) - Blunt-to-blunt and sticky-end compatibility.
