# blunt_cutters

List all built-in restriction enzymes that produce blunt ends.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `blunt_cutters` |
| **Method ID** | `RestrictionAnalyzer.GetBluntCutters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns every enzyme in the built-in Type II restriction-enzyme database that cuts both strands at the same position, leaving a blunt end (no single-stranded overhang). Blunt cutters are used for blunt-end cloning where any two blunt fragments can be ligated.

An enzyme is blunt when `CutPositionForward == CutPositionReverse`.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L107](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L107)

## Input Schema

_No parameters._

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `enzymes` | array | Blunt-cutting enzymes, each with `name`, `recognitionSequence`, `cutPositionForward`, `cutPositionReverse`, `organism`. |

## Errors

_None._

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

The built-in database contains exactly **10** blunt cutters: AluI, RsaI, HaeIII, DpnI, EcoRV, SmaI, HincII, ScaI, StuI, SwaI.

## Performance

- **Time Complexity:** O(n) over the enzyme database (n = 39).
- **Space Complexity:** O(k) for the k blunt cutters returned.

## See Also

- [sticky_cutters](sticky_cutters.md) - Enzymes producing cohesive (overhang) ends.
- [get_enzyme](get_enzyme.md) - Look up a single enzyme by name.
