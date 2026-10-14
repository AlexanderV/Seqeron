# isoelectric_point

Calculate isoelectric point (pI) of a protein sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `isoelectric_point` |
| **Method ID** | `SequenceStatistics.CalculateIsoelectricPoint` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the isoelectric point (pI) of a protein sequence. The pI is the pH at which the protein has no net electrical charge. Uses binary search to find the pH where net charge equals zero, considering ionizable side chains (D, E, C, Y, H, K, R) and terminal amino/carboxyl groups. Selenocysteine (U) and pyrrolysine (O) are accepted and treated as neutral (no pK), as EMBOSS 6.6 `iep` does (e.g. `MSKU` → 9.15, identical to `MSK`).

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L602](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L602)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The protein sequence (min length: 1) |
| `scale` | string | No | pK scale: `Emboss` (default; EMBOSS 6.6.0 `iep`) or `Bjellqvist` (ExPASy Compute pI/Mw, Biopython `IsoelectricPoint`); case-insensitive. E.g. INGAR: Bjellqvist 9.75 (Biopython 9.750021). |
| `window` | string | No | pH search convention (case-insensitive): `Default` (library [0,14] bisection, 2 dp), `Emboss` (`iep` `embIepPhConverge` over [1,14], unrounded; no pI → `found: false`), `Biopython` (`IsoelectricPoint.pi()` [4.05,12] starting at 7.775, unrounded). |
| `termini` | boolean | No | Count the N-/C-terminal charges (default true); `false` = `iep -notermini`. |
| `lysinesModified` | integer | No | Modified (uncharged) lysines, `iep -lysinemodified` (default 0; all K removed when fewer present). Must be ≥ 0. |
| `disulphides` | integer | No | Disulphide bridges, each removing 2 cysteine charges, `iep -disulphides` (default 0). Must be ≥ 0. |

Any non-default option delegates to `SequenceStatistics.TryCalculateIsoelectricPoint(sequence, IsoelectricPointOptions, out pI)` ([SequenceStatistics.cs#L760](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L760)); with all defaults the result is exactly `CalculateIsoelectricPoint(sequence, scale)`.

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pi` | number | Isoelectric point (pH). Omitted only when `window = "Emboss"` finds no pI. |
| `found` | boolean | Only for `window = "Emboss"`: `false` = `iep` "Isoelectric Point = None" (same strict charge sign at pH 1 and 14). Omitted otherwise, so the default output is unchanged. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 2002 | Invalid protein sequence |
| 1001 | Unknown pK scale '<scale>'. Expected one of: Emboss, Bjellqvist. |
| 1001 | Unknown pI window '<window>'. Expected one of: Default, Emboss, Biopython. |
| 1001 | lysinesModified must be non-negative / disulphides must be non-negative |

## Examples

### Example 1: Acidic protein

**User Prompt:**
> What's the isoelectric point of "MAEGEITTFT"?

**Expected Tool Call:**
```json
{
  "tool": "isoelectric_point",
  "arguments": {
    "sequence": "MAEGEITTFT"
  }
}
```

**Response:**
```json
{
  "pi": 3.61
}
```

### Example 2: Basic protein

**User Prompt:**
> Calculate pI of "MKRLLKR"

**Expected Tool Call:**
```json
{
  "tool": "isoelectric_point",
  "arguments": {
    "sequence": "MKRLLKR"
  }
}
```

**Response:**
```json
{
  "pi": 12.53
}
```

### Example 3: EMBOSS `iep` conventions

`iep -disulphides 2 -lysinemodified 2` on IFNA2_HUMAN (the `iep` documentation example, 5.7240):
`{"sequence": "<IFNA2_HUMAN>", "window": "Emboss", "lysinesModified": 2, "disulphides": 2}` →
`{"pi": 5.7240447998046875, "found": true}`. Poly-R (R×40) under `window: "Emboss"` → `{"found": false}`
(`iep`: "Isoelectric Point = None"). `{"sequence": "PETER", "window": "Emboss", "termini": false}` →
`{"pi": 4.099937438964844, "found": true}` (`iep -notermini` 4.0999). `{"sequence": "DDDD", "scale": "Bjellqvist",
"window": "Biopython"}` → `{"pi": 4.0500284194946286}` (Biopython 1.88 `IsoelectricPoint("DDDD").pi()`).

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [amino_acid_composition](amino_acid_composition.md) - Full protein composition analysis
- [protein_validate](protein_validate.md) - Validate protein sequences
