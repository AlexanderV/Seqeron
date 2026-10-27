# enzymes_compatible

Check whether two named enzymes produce ligatable ends.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `enzymes_compatible` |
| **Method ID** | `RestrictionAnalyzer.AreCompatible` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns `true` when the two named enzymes produce ligatable ends. Names are case-insensitive.

**The rule.** Two ends anneal only when their single-stranded extensions have the same length, the
**same polarity** (both 5′ or both 3′) and are reverse complements of each other; blunt ends ligate
to blunt ends. A cut leaves two ends — the downstream fragment's extension reads the overhang `W` on
the protruding strand and the upstream fragment's reads `revcomp(W)` on the other — so an enzyme
offers both, and two enzymes are compatible iff they share a polarity and `W1` matches `W2` or
`revcomp(W2)`. Where an overhang carries an IUPAC code the match is per-position code intersection,
so a degenerate verdict means "**possibly** compatible — some pair of real sites ligates".

**Same letters with opposite polarity are NOT compatible.** `SacI` leaves a 3′ `AGCT` and `HindIII`
a 5′ `AGCT`, yet the two cannot be joined: both extensions would lie on the same strand of the
joint. All four such pairs in the built-in table — `HindIII`/`SacI` (`AGCT`), `NcoI`/`SphI`
(`CATG`), `NotI`/`ApaI` (`GGCC`), `AgeI`/`FseI` (`CCGG`) — return `false`.

`compatibleEnd` is `"blunt"` for a blunt pair, otherwise the **first** enzyme's own overhang window,
and `null` when the two are incompatible; `end2` is the second enzyme's window. It is **not**
necessarily a string both enzymes leave: the two coincide whenever the windows intersect directly,
but 12 REBASE pairs anneal only through the reverse complement, and all 12 are `GN` against `CN` —
`Bse1I`/`BsmI`, `Bse1I`/`LmnI`, `Bse1I`/`Mva1269I`, `Bse1I`/`PctI`, `BseNI`/`BsmI`, `BseNI`/`LmnI`,
`BseNI`/`Mva1269I`, `BseNI`/`PctI`, `BsmI`/`BsrI`, `BsrI`/`LmnI`, `BsrI`/`Mva1269I`, `BsrI`/`PctI`.
`enzymes_compatible("BsmI", "BsrI")` answers `compatibleEnd: "CN", end2: "GN"` and the reversed call
answers `compatibleEnd: "GN", end2: "CN"`; BsrI never leaves `CN`. The paragraph here used to say the
value "echoes the shared end", and for those 12 pairs no single shared string exists
(review 2026-09 B08 H30). The other 18 541 sticky REBASE pairs and every blunt pair are unaffected.

**An unknown enzyme name is an ERROR** (changed 2026-10-14). It used to resolve to `false`, which
made a typo indistinguishable from a genuinely incompatible pair — and "not compatible" is the
answer that silently loses a real cloning option. `find_restriction_sites` already errors on an
unknown name; this tool now matches. Use `get_enzyme` first if the name may not exist.

Verified against Biopython 1.88 on 2026-10-14: all 1521 ordered pairs of the 39 built-ins agree with
`a % b or b % a`, and 75 466 pairs of the 389 REBASE enzymes with a fully-defined overhang agree
with `a % b` — 0 disagreements in both.

## Core Documentation Reference

- Source: [RestrictionAnalyzer.cs#L1840](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs#L1840)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `enzyme1_name` | string | Yes | First enzyme name (non-blank). |
| `enzyme2_name` | string | Yes | Second enzyme name (non-blank). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `compatible` | boolean | True if the two enzymes' ends can be ligated. |
| `compatibleEnd` | string \| null | The **first** enzyme's end — `"blunt"` or its own overhang window — or `null` when incompatible. Not necessarily a string the second enzyme also leaves. |
| `end2` | string \| null | The **second** enzyme's end, in the same convention; `null` when incompatible. Equal to `compatibleEnd` except for the 12 reverse-complement-only pairs. |

## Errors

| Code | Message |
|------|---------|
| 1001 | First enzyme name cannot be null or blank |
| 1002 | Second enzyme name cannot be null or blank |
| 1003 | Unknown enzyme: &lt;name&gt; |

## Examples

### Example 1: `EcoRV` + `SmaI` (both blunt) → `{ compatible: true, compatibleEnd: "blunt" }`.

### Example 2: `EcoRI` (5′ AATT) + `BamHI` (5′ GATC) → `{ compatible: false, compatibleEnd: null }` — same polarity and length, different sequences.

### Example 3: `BamHI` + `BglII` → `{ compatible: true, compatibleEnd: "GATC" }`. The whole 5′ `GATC` family — `BamHI`, `BglII`, `MboI`, `Sau3AI` — is mutually compatible, as are `SalI`/`XhoI` (`TCGA`), `NheI`/`XbaI`/`SpeI`/`AvrII` (`CTAG`) and `ClaI`/`TaqI`/`MspI`/`HpaII` (`CG`).

### Example 4: `HindIII` + `SacI` → `{ compatible: false, compatibleEnd: null }` — **both overhangs read `AGCT`**, but HindIII's is 5′ and SacI's is 3′.

### Example 5: `SfiI` + `PstI` → `false`. SfiI's 3′ overhang is 3 nt (`NNN`) and PstI's is 4 nt (`TGCA`); different lengths never anneal. Biopython agrees (`SfiI % PstI` is `False`).

### Example 7: `BsmI` + `BsrI` → `{ compatible: true, compatibleEnd: "CN", end2: "GN" }`, and `BsrI` + `BsmI` → `{ compatible: true, compatibleEnd: "GN", end2: "CN" }` — the pair anneals only through the reverse complement, so each call reports the first enzyme's own end (B08 H30).

### Example 6: `EcoRI` + `NotAnEnzyme` → **error** `Unknown enzyme: NotAnEnzyme` (it used to return `false`).

## See Also

- [compatible_enzymes](compatible_enzymes.md) - Enumerate all compatible pairs, or one enzyme's compatible set.
- [get_enzyme](get_enzyme.md) - The enzyme's overhang sequence and polarity.
