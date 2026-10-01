# find_longdust_regions

Low-complexity regions with longdust (Li & Li 2025), the k-mer generalisation of SDUST.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_longdust_regions` |
| **Method ID** | `SequenceComplexity.FindLongdustRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds low-complexity regions with **longdust** (Li & Li 2025; port of lh3/longdust 1.4-r97 `ld_dust1/ld_dust2`,
MIT) — the k-mer generalisation of SDUST for long windows (STRs, VNTRs, satellites). At each position a backward then
forward scan over the last `windowSize` k-mers finds a good interval with `S_L(x) − T·ℓ(x) > 0` (see
`longdust_score`), with X-drop and the reference's speed-ups; overlapping hits are merged; by default the union over
both strands is returned (`forwardOnly` = `-f`). Output = longdust BED intervals, 0-based half-open. Non-ACGT symbols
make the overlapping k-mers ambiguous. Defaults = `longdust -k7 -w5000 -t0.6 -e50 -b3`.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L1195](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L1195)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (A/C/G/T plus IUPAC codes such as N; case-insensitive) |
| `k` | integer | No | k-mer length (-k) (default 7) |
| `windowSize` | integer | No | Window size in k-mers (-w) (default 5000) |
| `threshold` | number | No | Score threshold T per k-mer (-t) (default 0.6) |
| `xdropLength` | integer | No | X-drop length (-e; 0 disables X-drop) (default 50) |
| `minStartCount` | integer | No | Minimum count of the current k-mer before a search starts (-b) (default 3) |
| `forwardOnly` | boolean | No | Forward strand only (-f) (default false) |
| `approximate` | boolean | No | Guaranteed O(Lw) one-pass mode (-a) (default false) |
| `gcContent` | number | No | Genome GC fraction in (0, 1) for GC correction (-g); omit = off |

## Output Schema

`items`: `start, end` (0-based half-open, longdust BED), `length`

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1003 | k 1-14; windowSize 1-65534; threshold finite > 0; xdropLength >= 0; minStartCount >= 2; gcContent in (0, 1) |

## Examples

### Example 1: VNTR, longdust defaults

**Input:** `{"sequence": "ATCAGTCATTAAACTATAAACCACTTGAACCACAACGATGTCGTTTATAGCGCGCGGGGACGGCAGCTGCGATACCCCCTCGAATCCCCGGCGGCTCTCACCTGCAGGGTGGACGTTTGGGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAACCGAGCCTCAACGGAAAGGCGGCATTGGGCGTAGATCATTGTAAGAATTGAGAGGACTGAGGGATAGGGAAAGGTACGGGCCCCGATTTCCCATGCAGGCATCTCCAAGTGTAAGCACG"}`

**Output:**

```json
{"items": [{"start": 120, "end": 288, "length": 168}]}
```

### Example 2: longdust -k5 -w100

**Input:** `{"sequence": "ATCAGTCATTAAACTATAAACCACTTGAACCACAACGATGTCGTTTATAGCGCGCGGGGACGGCAGCTGCGATACCCCCTCGAATCCCCGGCGGCTCTCACCTGCAGGGTGGACGTTTGGGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAACCGAGCCTCAACGGAAAGGCGGCATTGGGCGTAGATCATTGTAAGAATTGAGAGGACTGAGGGATAGGGAAAGGTACGGGCCCCGATTTCCCATGCAGGCATCTCCAAGTGTAAGCACG", "k": 5, "windowSize": 100}`

**Output:**

```json
{"items": [{"start": 117, "end": 291, "length": 174}]}
```

### Example 3: Only N

**Input:** `{"sequence": "NNNNNNNNNN"}`

**Output:**

```json
{"items": []}
```

## Performance

- O(n · windowSize) worst case (O(n · windowSize) guaranteed with approximate).

## See Also

- [longdust_score](longdust_score.md)
- [find_low_complexity_intervals](find_low_complexity_intervals.md)
- [mask_low_complexity](mask_low_complexity.md)
