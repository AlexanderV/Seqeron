# find_microsatellites

Short Tandem Repeats (STRs / microsatellites) in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_microsatellites` |
| **Method ID** | `RepeatFinder.FindMicrosatellites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds **microsatellites** (Short Tandem Repeats): units of length
`minUnitLength..maxUnitLength` repeated at least `minRepeats` times consecutively.
Each hit is classified by repeat type (Mononucleotide, Dinucleotide, …). Redundant
units — those that are themselves a shorter unit repeated (e.g. `AA`) — are skipped.
Per unit length each maximal perfect run is reported once, at its left end, with the
number of complete copies (a trailing partial copy is not counted, and rotations of
the same run — e.g. `TA` inside `ATATATA` — are not re-reported); units containing
non-ACGT symbols (e.g. `N`) are never reported. Positions are 0-based.

Optional MISA behaviour (Thiel et al. 2003, `misa.pl` v1.0; verified against a real `perl misa.pl` run):

- `misaThresholds: true` — minimum copies per unit length from MISA's default `misa.ini`
  (`1-10 2-6 3-5 4-5 5-5 6-5`) for the unit lengths `minUnitLength..maxUnitLength` within 1–6;
  `minRepeats` is ignored. Items are ordered by unit length, then position.
- `misaScan = true` — reproduces misa.pl v1.0's SSR list exactly (`RepeatFinder.FindMicrosatellites(…, MicrosatelliteScanMode.MisaRegex)`):
  per unit length the leftmost greedy match of `([acgt]{p})\2{t-1,}` is taken and the scan resumes at its end;
  a non-primitive match (e.g. `CTCTCT`) is rejected after consuming its bases. Differs from the default only where
  a same-size run overlaps the previous match by < p bases or starts inside a rejected match (cross-checked:
  6 000 sequences × 6 misa.ini settings, 0 differences; SSRs sharing a start are chained in unit-length order —
  misa.pl orders them by Perl's randomised hash order).
- `maxCompoundInterruption ≥ 0` — also returns `compounds`: MISA compound microsatellites, i.e.
  chains of STRs (in start order) where each STR starts at most that many bases after the previous
  one ends (adjacent and overlapping STRs always join; MISA default 100). `type` is `c`, or `c*`
  when two components overlap; `notation` is MISA's SSR string (e.g. `(TA)6tccgt(GA)7ttttt(A)12`,
  interruptions lower-case); `end` is the end of the last component (MISA's `end` column).

## Core Documentation Reference

- Source: [RepeatFinder.cs#L86](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L86)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `minUnitLength` | integer | No | Minimum unit length (default 1, ≥ 1) |
| `maxUnitLength` | integer | No | Maximum unit length (default 6) |
| `minRepeats` | integer | No | Minimum repeats (default 3, ≥ 2); ignored when `misaThresholds` is true |
| `misaThresholds` | boolean | No | MISA default per-unit-size minimum copies `1-10 2-6 3-5 4-5 5-5 6-5` (default false) |
| `maxCompoundInterruption` | integer | No | ≥ 0: also return MISA compound microsatellites with at most this many interrupting bases (MISA default 100); default −1 = none |
| `misaScan` | boolean | No | misa.pl's regex scan instead of maximal primitive runs (default false); see below |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ position, repeatUnit, repeatCount, totalLength, repeatType }` |
| `compounds` | array/null | `{ start, end, length, type ("c"/"c*"), notation, components[] }`; null unless `maxCompoundInterruption ≥ 0` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Invalid unit-length or minRepeats bounds |

## Examples

### Example 1: Dinucleotide (CA)₄

**User Prompt:**
> Find microsatellites in "CACACACA" (unit 2–6, ≥ 3 repeats).

**Expected Tool Call:**
```json
{
  "tool": "find_microsatellites",
  "arguments": { "sequence": "CACACACA", "minUnitLength": 2, "maxUnitLength": 6, "minRepeats": 3 }
}
```

**Response:**
```json
{ "items": [ { "position": 0, "repeatUnit": "CA", "repeatCount": 4, "totalLength": 8, "repeatType": "Dinucleotide" } ] }
```

### Example 2: Trinucleotide (CAG)₃

**User Prompt:**
> Trinucleotide STRs in "CAGCAGCAG".

**Expected Tool Call:**
```json
{
  "tool": "find_microsatellites",
  "arguments": { "sequence": "CAGCAGCAG", "minUnitLength": 3, "maxUnitLength": 6, "minRepeats": 3 }
}
```

**Response:**
```json
{ "items": [ { "position": 0, "repeatUnit": "CAG", "repeatCount": 3, "totalLength": 9, "repeatType": "Trinucleotide" } ] }
```

### Example 3: MISA thresholds and compound microsatellite

**User Prompt:**
> MISA-style SSR search with compound SSRs in "ACGTATATATATATATccgtGAGAGAGAGAGAGAtttttAAAAAAAAAAAAT".

**Expected Tool Call:**
```json
{
  "tool": "find_microsatellites",
  "arguments": { "sequence": "ACGTATATATATATATccgtGAGAGAGAGAGAGAtttttAAAAAAAAAAAAT", "misaThresholds": true, "maxCompoundInterruption": 100 }
}
```

**Response** (`perl misa.pl`: `c (TA)6tccgt(GA)7ttttt(A)12 48 4 51`, 1-based):
```json
{ "items": [ { "position": 39, "repeatUnit": "A", "repeatCount": 12, "totalLength": 12, "repeatType": "Mononucleotide" },
             { "position": 3, "repeatUnit": "TA", "repeatCount": 6, "totalLength": 12, "repeatType": "Dinucleotide" },
             { "position": 20, "repeatUnit": "GA", "repeatCount": 7, "totalLength": 14, "repeatType": "Dinucleotide" } ],
  "compounds": [ { "start": 3, "end": 51, "length": 48, "type": "c", "notation": "(TA)6tccgt(GA)7ttttt(A)12", "components": [ "…the three items in start order…" ] } ] }
```

## Performance

- **Time Complexity:** O(n · (maxUnitLength − minUnitLength + 1)) character comparisons.
- **Space Complexity:** O(number of STRs).

## See Also

- [find_tandem_repeats](find_tandem_repeats.md) — general tandem repeats
- [tandem_repeat_summary](tandem_repeat_summary.md) — aggregate STR statistics
