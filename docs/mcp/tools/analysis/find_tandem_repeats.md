# find_tandem_repeats

Consecutive repeating units (tandem repeats) in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_tandem_repeats` |
| **Method ID** | `GenomicAnalyzer.FindTandemRepeats` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **exact tandem repeats** — a unit repeated consecutively (e.g. `ATGATGATG`) — whose unit
length lies in `[minUnitLength, maxUnitLength]` and which repeat ≥ `minRepetitions` complete
times. For each unit length every **maximal run** (cannot be extended left or right with the same
period; Kolpakov & Kucherov 1999) with a **primitive** unit (not a power of a shorter word:
`ATATATAT` is `AT×4` at 0, not also `ATAT×2`; MISA "reject false type motifs") is reported once
at its left end with `⌊length / unit length⌋` complete copies. Each result carries the unit, the
0-based start, the number of copies and the total span (`unit.Length × repetitions`). Order: unit
length ascending, then start. Runs of different unit lengths may overlap and are reported
independently. Input is case-insensitive (A/C/G/T only).

Implementation: a facade over `RepeatFinder.FindMicrosatellites` (one linear pass per unit length)
that switches to a run-enumeration engine (Lyndon-array runs algorithm, Bannai et al. 2017,
O(n log² n)) when more than 256 unit lengths are admissible, so the default "no upper bound" is
fast on megabase sequences. **Approximate (interrupted) repeats** are not covered here: use
[`find_approximate_tandem_repeats`](find_approximate_tandem_repeats.md) (Tandem Repeats Finder
model, Benson 1999).

**Canonical rotation.** `canonicalRotation: true` reports each unit as its lexicographically least rotation (Booth 1980; the representative of MISA repeat-type classes), so `TGATGATGAT` gives `ATG` x3 at 0 with `originalUnit` `TGA`; `position`, `repetitions` and `totalLength` are unchanged. Verified against a brute-force least rotation on the MUMmer `exact-tandems` oracle (3630 rows).

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (min length 1): A/C/G/T plus IUPAC codes such as N, case-insensitive. A repeat unit must be spelled with A/C/G/T only, so a window containing N or another ambiguity code never forms a repeat (`ACACNACAC` → AC×2 at 0 and AC×2 at 5) |
| `minUnitLength` | integer | No | Minimum repeat-unit length (default 2, ≥ 1) |
| `minRepetitions` | integer | No | Minimum complete consecutive copies (default 2, ≥ 2) |
| `maxUnitLength` | integer | No | Maximum repeat-unit length (≥ `minUnitLength`); omit for no upper bound |
| `canonicalRotation` | boolean | No | Report each unit as its lexicographically least rotation (Booth 1980) instead of the spelling at the run's left end; originalUnit then holds the actual spelling. Default false. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ unit, position, repetitions, totalLength }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1002 | Minimum repeat-unit length must be at least 1 |
| 1002 | minRepetitions must be ≥ 2 |
| 1002 | Maximum repeat-unit length must be at least minUnitLength |

## Examples

### Example 1: Trinucleotide tandem (ATG ×3)

**User Prompt:**
> Find tandem 3-mers in "ATGATGATG".

**Expected Tool Call:**
```json
{
  "tool": "find_tandem_repeats",
  "arguments": { "sequence": "ATGATGATG", "minUnitLength": 3, "minRepetitions": 2 }
}
```

**Response:**
```json
{ "items": [ { "unit": "ATG", "position": 0, "repetitions": 3, "totalLength": 9 } ] }
```

### Example 2: Mononucleotide run (unit length 1)

**User Prompt:**
> Find tandem repeats in "AAAAA" with unit length ≥ 1.

**Expected Tool Call:**
```json
{
  "tool": "find_tandem_repeats",
  "arguments": { "sequence": "AAAAA", "minUnitLength": 1, "minRepetitions": 2 }
}
```

**Response:**
```json
{ "items": [ { "unit": "A", "position": 0, "repetitions": 5, "totalLength": 5 } ] }
```
The run is reported once with its primitive unit; `AA×2` is the same run tiled by a
non-primitive unit and is not reported.

## Performance

- **Time Complexity:** O(n) per admissible unit length (≤ 256 unit lengths), otherwise O(n log² n) once
  (1 Mbp random sequence, unbounded unit length: ≈ 3 s in a Debug build; 1 Mbp homopolymer / period-8 string ≈ 6 s).
- **Space Complexity:** O(1) plus output for the per-unit-length scan; O(n) for the run enumeration.

## See Also

- [find_microsatellites](find_microsatellites.md) — STRs with repeat-type classification (same engine)
- [find_approximate_tandem_repeats](find_approximate_tandem_repeats.md) — interrupted / imperfect repeats (TRF)
- [find_repeats](find_repeats.md) — any repeated substring
- [tandem_repeat_summary](tandem_repeat_summary.md) — aggregate statistics
