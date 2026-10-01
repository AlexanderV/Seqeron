# tandem_repeat_summary

Aggregate statistics across all microsatellites in a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `tandem_repeat_summary` |
| **Method ID** | `RepeatFinder.GetTandemRepeatSummary` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Summarizes the perfect microsatellites (STRs, unit length 1–6, or the unit sizes of a custom `misaDefinition`) reported by
`find_microsatellites` for the same `minRepeats` (each maximal primitive run once per unit
length; no second scan). Conventions follow the MISA `.statistics` output (Thiel et al. 2003)
and Krait statistics (Du et al. 2018):

- `totalRepeats` — number of STRs; `countsByUnitLength` (one entry per searched unit size, e.g. `{"1": 0, …, "6": 1}`;
  every size of `misaDefinition`, incl. sizes > 6) sums to it — misa.pl `.statistics` "Distribution to different repeat
  type classes" (misa.pl omits sizes with 0 SSRs). The six per-class counts (mono … hexa) are its entries 1–6 and sum to
  `totalRepeats` unless sizes > 6 were searched.
- `totalRepeatBases` — sum of the STR lengths (Krait "Length (bp)"); runs of different unit
  lengths that overlap (e.g. `A×5` inside `AAAAATATATAT`) are each counted in full.
- `percentageOfSequence` — bases covered by at least one STR (union of spans) / length × 100,
  so always 0–100; equals `totalRepeatBases / length × 100` when no two STRs overlap.
- `longestRepeat` — largest `totalLength`; ties → shorter unit, then leftmost. `null` when none.
- `mostFrequentUnit` — the reported unit string (motif phase at the run start; not rotation- or
  reverse-complement-canonicalized, as in MISA's "Frequency of identified SSR motifs") occurring
  in the most STRs; ties → the unit whose first STR comes first in the unit-length/position order.
  `null` when none.
- `canonicalMotifCounts` — STR count per MISA repeat-type class "considering sequence
  complementary" (`misa.pl` `.statistics`): motif rotations and reverse-complement rotations are
  one class named `X/Y` (smallest rotation of each strand, smaller first), e.g. AC, CA, GT, TG →
  `AC/GT`; A, T → `A/T`. Verified against a real `perl misa.pl` run.
- `standardMotifCounts` — with `standardMotifLevel` 0–4, STR count per Krait standard motif
  (Du et al. 2018 `motif.py`; see `standardize_repeat_motif`), e.g. level 2 → `{ "A": 2, "AC": 4, "ATAC": 1 }`;
  `null` when `standardMotifLevel` is −1 (default).

With `misaThresholds: true` the STRs use MISA's default per-unit-size minimum copies
(`1-10 2-6 3-5 4-5 5-5 6-5`) instead of one `minRepeats`. With `misaScan: true` the STRs come from
misa.pl's regex scan (leftmost greedy match resumed after each match; non-primitive matches consumed
then rejected — `find_microsatellites` `misaScan`), so `totalRepeats` and the per-unit-size counts
equal misa.pl's `.statistics` ("Total number of identified SSRs", "Distribution to different repeat
type classes"; 6 000 N/IUPAC-containing sequences × 6 definitions, 0 mismatches).

With `misaDefinition` (misa.ini `def` syntax, e.g. `"1-10 2-6 3-5 4-5 5-5 6-5 7-5 8-5"`; any unit size, as misa.pl)
the STRs use exactly those unit sizes and thresholds (`minRepeats` ignored; cannot be combined with `misaThresholds`);
with `misaScan: true` the total and every `countsByUnitLength` entry equal misa.pl's `.statistics` (6 000 sequences ×
7 definitions with unit sizes 7–12, 0 mismatches; B04 F61). STRs longer than 6 bp have `repeatType` `Complex`.

N and other IUPAC codes are accepted (case-insensitive), as in `find_microsatellites`: only A/C/G/T
form units, so N/IUPAC symbols never belong to an STR and interrupt runs (MISA `[acgt]`). The
`percentageOfSequence` denominator is the full length, N included (misa.pl "Total size of examined
sequences" = `length $seq`).

## Core Documentation Reference

- Source: [RepeatFinder.cs#L6352](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs#L6352)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence: A/C/G/T plus IUPAC codes such as N, case-insensitive (min length 1) |
| `minRepeats` | integer | No | Minimum complete copies for every unit length (default 3, ≥ 2; smaller values throw `ArgumentOutOfRangeException`); ignored when `misaThresholds` is true |
| `misaThresholds` | boolean | No | MISA default per-unit-size minimum copies `1-10 2-6 3-5 4-5 5-5 6-5` (default false) |
| `standardMotifLevel` | integer | No | −1 (default, off) or Krait standardization level 0–4 for `standardMotifCounts` |
| `misaScan` | boolean | No | misa.pl regex scan instead of maximal primitive runs (default false) |
| `misaDefinition` | string/null | No | Custom MISA definition, misa.ini `def` syntax (`size-min` pairs separated by spaces/commas, any unit size ≥ 1, min ≥ 2, each size once; a leading `definition(...):` token is skipped), e.g. `"1-10 2-6 3-5 4-5 5-5 6-5 7-5"`; replaces `minRepeats`; not combinable with `misaThresholds` (default null) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalRepeats` | integer | Number of microsatellites found |
| `totalRepeatBases` | integer | Sum of STR lengths (overlapping STRs each counted in full) |
| `percentageOfSequence` | number | Percent of the sequence covered by the union of STR spans (0–100) |
| `mononucleotideRepeats` | integer | Count of mononucleotide STRs |
| `dinucleotideRepeats` | integer | Count of dinucleotide STRs |
| `trinucleotideRepeats` | integer | Count of trinucleotide STRs |
| `tetranucleotideRepeats` | integer | Count of tetranucleotide STRs |
| `pentanucleotideRepeats` | integer | Count of pentanucleotide STRs |
| `hexanucleotideRepeats` | integer | Count of hexanucleotide STRs |
| `countsByUnitLength` | object | STR count per searched unit size (keys = unit size in bp: 1–6, or every size of `misaDefinition`); sums to `totalRepeats` |
| `longestRepeat` | object/null | The longest microsatellite (or null) |
| `mostFrequentUnit` | string/null | The most frequent repeat unit (or null) |
| `canonicalMotifCounts` | object | STR count per MISA repeat-type class, e.g. `{ "AC/GT": 3, "A/T": 1 }` |
| `standardMotifCounts` | object/null | STR count per Krait standard motif (only when `standardMotifLevel` is 0–4) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1001 | `misaDefinition`: not a MISA `size-min` pair / size defined twice / size < 1 / minimum < 2 / combined with `misaThresholds` |

## Examples

### Example 1: Single trinucleotide STR

**User Prompt:**
> Summarize tandem repeats in "CAGCAGCAG".

**Expected Tool Call:**
```json
{
  "tool": "tandem_repeat_summary",
  "arguments": { "sequence": "CAGCAGCAG", "minRepeats": 3 }
}
```

**Response:**
```json
{ "totalRepeats": 1, "totalRepeatBases": 9, "percentageOfSequence": 100.0, "mononucleotideRepeats": 0, "dinucleotideRepeats": 0, "trinucleotideRepeats": 1, "tetranucleotideRepeats": 0, "pentanucleotideRepeats": 0, "hexanucleotideRepeats": 0, "countsByUnitLength": { "1": 0, "2": 0, "3": 1, "4": 0, "5": 0, "6": 0 }, "longestRepeat": { "position": 0, "repeatUnit": "CAG", "repeatCount": 3, "totalLength": 9, "repeatType": "Trinucleotide" }, "mostFrequentUnit": "CAG", "canonicalMotifCounts": { "AGC/CTG": 1 } }
```

### Example 2: No STRs

**User Prompt:**
> Tandem repeat summary of "ACGT".

**Expected Tool Call:**
```json
{
  "tool": "tandem_repeat_summary",
  "arguments": { "sequence": "ACGT", "minRepeats": 3 }
}
```

**Response:**
```json
{ "totalRepeats": 0, "totalRepeatBases": 0, "percentageOfSequence": 0.0, "mononucleotideRepeats": 0, "dinucleotideRepeats": 0, "trinucleotideRepeats": 0, "tetranucleotideRepeats": 0, "pentanucleotideRepeats": 0, "hexanucleotideRepeats": 0, "countsByUnitLength": { "1": 0, "2": 0, "3": 0, "4": 0, "5": 0, "6": 0 }, "longestRepeat": null, "mostFrequentUnit": null, "canonicalMotifCounts": {} }
```

### Example 3: Penta- and hexanucleotide STRs, overlapping homopolymers

**User Prompt:**
> Summarize the STRs in "AAAGAAAAGAAAAGAAAAGACCTTAGGGTTAGGGTTAGGGTTAGGG" (Penta E motif + telomere repeat).

**Expected Tool Call:**
```json
{
  "tool": "tandem_repeat_summary",
  "arguments": { "sequence": "AAAGAAAAGAAAAGAAAAGACCTTAGGGTTAGGGTTAGGGTTAGGG", "minRepeats": 3 }
}
```

**Response** (STRs: `A×3@0, A×4@4, A×4@9, A×4@14, G×3@25, G×3@31, G×3@37, G×3@43, (AAAGA)4@0, (TTAGGG)4@22`;
71 repeat bases, but only `[0,20) ∪ [22,46)` = 44 of 46 bases covered):
```json
{ "totalRepeats": 10, "totalRepeatBases": 71, "percentageOfSequence": 95.65217391304348, "mononucleotideRepeats": 8, "dinucleotideRepeats": 0, "trinucleotideRepeats": 0, "tetranucleotideRepeats": 0, "pentanucleotideRepeats": 1, "hexanucleotideRepeats": 1, "countsByUnitLength": { "1": 8, "2": 0, "3": 0, "4": 0, "5": 1, "6": 1 }, "longestRepeat": { "position": 22, "repeatUnit": "TTAGGG", "repeatCount": 4, "totalLength": 24, "repeatType": "Hexanucleotide" }, "mostFrequentUnit": "A", "canonicalMotifCounts": { "A/T": 4, "C/G": 4, "AAAAG/CTTTT": 1, "AACCCT/AGGGTT": 1 } }
```

### Example 4: N and IUPAC codes, misa.pl scan

**User Prompt:**
> MISA-style SSR statistics for "ACACACACANCACACACACACACACAGTNNNNNNNNNNNNNNNNNNNNATATATATATATATATATRYSWKMggggggggggggggg".

**Expected Tool Call:**
```json
{
  "tool": "tandem_repeat_summary",
  "arguments": { "sequence": "ACACACACANCACACACACACACACAGTNNNNNNNNNNNNNNNNNNNNATATATATATATATATATRYSWKMggggggggggggggg", "misaThresholds": true, "misaScan": true }
}
```

**Response** (misa.pl `.statistics`: size 87, 3 SSRs, unit sizes 1: 1, 2: 2; `.misa` row
`c (CA)8gtnnnnnnnnnnnnnnnnnnnn(AT)9ryswkm(G)15 77 11 87` → 16 + 18 + 15 = 49 bases, 49 / 87):
```json
{ "totalRepeats": 3, "totalRepeatBases": 49, "percentageOfSequence": 56.32183908045977, "mononucleotideRepeats": 1, "dinucleotideRepeats": 2, "trinucleotideRepeats": 0, "tetranucleotideRepeats": 0, "pentanucleotideRepeats": 0, "hexanucleotideRepeats": 0, "countsByUnitLength": { "1": 1, "2": 2, "3": 0, "4": 0, "5": 0, "6": 0 }, "longestRepeat": { "position": 48, "repeatUnit": "AT", "repeatCount": 9, "totalLength": 18, "repeatType": "Dinucleotide" }, "mostFrequentUnit": "G", "canonicalMotifCounts": { "C/G": 1, "AC/GT": 1, "AT/AT": 1 } }
```

### Example 5: Custom MISA definition with 7–10 bp units

**User Prompt:**
> misa.pl statistics with definition `1-10 2-6 3-5 4-5 5-5 6-5 7-3 8-3 9-2 10-2` for "ACGTTGCACGTTGCACGTTGCGAAAAAAAAAAAACCTTAGGCATTAGGCATTAGGCATTCGACGTACGTACGTACGTACGTACGTGGATCCATGCAATCCATGCAATTGCACCTTGAGACCTTGAGANNACCTTGAGAGG".

**Expected Tool Call:**
```json
{
  "tool": "tandem_repeat_summary",
  "arguments": { "sequence": "ACGTTGCACGTTGCACGTTGCGAAAAAAAAAAAACCTTAGGCATTAGGCATTAGGCATTCGACGTACGTACGTACGTACGTACGTGGATCCATGCAATCCATGCAATTGCACCTTGAGACCTTGAGANNACCTTGAGAGG", "misaScan": true, "misaDefinition": "1-10 2-6 3-5 4-5 5-5 6-5 7-3 8-3 9-2 10-2" }
}
```

**Response** (misa.pl `.statistics`: 5 SSRs; distribution 1: 1, 4: 1, 7: 2, 9: 1):
```json
{"totalRepeats":5,"totalRepeatBases":96,"percentageOfSequence":68.57142857142857,"mononucleotideRepeats":1,"dinucleotideRepeats":0,"trinucleotideRepeats":0,"tetranucleotideRepeats":1,"longestRepeat":{"position":61,"repeatUnit":"ACGT","repeatCount":6,"totalLength":24,"repeatType":"Tetranucleotide"},"mostFrequentUnit":"A","pentanucleotideRepeats":0,"hexanucleotideRepeats":0,"countsByUnitLength":{"1":1,"2":0,"3":0,"4":1,"5":0,"6":0,"7":2,"8":0,"9":1,"10":0},"canonicalMotifCounts":{"A/T":1,"ACGT/ACGT":1,"AACGTGC/ACGTTGC":1,"AATGCCT/AGGCATT":1,"AATCCATGC/ATGGATTGC":1},"standardMotifCounts":null}
```

## Performance

- **Time Complexity:** O(u·n) STR scan for u searched unit sizes (6 by default) (`FindMicrosatellites`) + O(k log k) grouping/sorting of the k STRs.
- **Space Complexity:** O(number of STRs).

## See Also

- [find_microsatellites](find_microsatellites.md) — the per-STR list
- [find_tandem_repeats](find_tandem_repeats.md) — general tandem repeats
- [standardize_repeat_motif](standardize_repeat_motif.md) — MISA class / Krait standard motif of one unit
