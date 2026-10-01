# kmer_spectrum

Frequency-of-frequencies (k-mer spectrum) for a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `kmer_spectrum` |
| **Method ID** | `KmerAnalyzer.GetKmerSpectrum` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **k-mer spectrum** (frequency-of-frequencies): after counting every
overlapping k-mer, it reports, for each occurrence count `c`, how many *distinct*
k-mers occur exactly `c` times. This is the classic sequencing-coverage histogram
used to distinguish erroneous (low-count) from genomic (high-count) k-mers.
Counting is case-insensitive. When `k` exceeds the sequence length the spectrum is
empty.

Optional parameters (all additive; without them the result is unchanged):

- `canonical` / `acgtOnly` — Jellyfish counting modes as in `count_kmers` (`count -C`; skip non-ACGT windows).
- `low`, `high`, `increment`, `full` — the `jellyfish histo` options `-l` (default 1), `-h` (default 10000),
  `-i` (default 1) and `-f`. Giving any of them adds `histogram`: the rows `jellyfish histo` prints, computed
  exactly as `sub_commands/histo_main.cc` (base = increment ≥ low ? 0 : low − increment; bucket labels
  base + i·increment; counts above `high` pooled in the last bucket, labelled ≥ `high`; counts below `low`
  pooled in the first; empty buckets listed only with `full`). Cross-checked against Jellyfish 2.3.1 on 72 runs.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L423](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L423)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `canonical` | boolean | No | Canonical k-mers (`jellyfish count -C`, key = min(k-mer, reverse complement)); implies `acgtOnly`. Default `false` |
| `acgtOnly` | boolean | No | Skip windows containing a symbol other than A/C/G/T (Jellyfish convention). Default `false` |
| `low` | integer | No | `jellyfish histo -l` (≥ 0, default 1) |
| `high` | integer | No | `jellyfish histo -h` (≥ `low`, default 10000); higher counts go to the cap bin |
| `increment` | integer | No | `jellyfish histo -i` bucket width (≥ 1, default 1) |
| `full` | boolean | No | `jellyfish histo -f`: also list empty bins. Default `false` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `spectrum` | object | Map from occurrence count → number of distinct k-mers reaching that count |
| `histogram` | array \| null | `{bin, frequency}` rows of `jellyfish histo` (ascending); null unless a histo option is given |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |
| 1003 | High count value must be >= low count value; increment must be >= 1; low must be >= 0 |

## Examples

### Example 1: Monomer spectrum (GTAGAGCTGT, k=1)

**User Prompt:**
> Give me the k-mer spectrum of "GTAGAGCTGT" for k=1.

**Expected Tool Call:**
```json
{
  "tool": "kmer_spectrum",
  "arguments": { "sequence": "GTAGAGCTGT", "k": 1 }
}
```

**Response:**
```json
{ "spectrum": { "4": 1, "3": 1, "2": 1, "1": 1 } }
```
Monomer counts are G=4, T=3, A=2, C=1 — four distinct counts, one k-mer each.

### Example 2: Repeated trimer (ATGATG, k=3)

**User Prompt:**
> k-mer spectrum of "ATGATG" for k=3.

**Expected Tool Call:**
```json
{
  "tool": "kmer_spectrum",
  "arguments": { "sequence": "ATGATG", "k": 3 }
}
```

**Response:**
```json
{ "spectrum": { "2": 1, "1": 2 } }
```
Trimers: ATG×2, TGA×1, GAT×1 ⇒ one k-mer at count 2, two k-mers at count 1.

### Example 3: Canonical histo with empty bins (BA1B sample, k=4)

**Expected Tool Call:**
```json
{
  "tool": "kmer_spectrum",
  "arguments": { "sequence": "ACGTTGCATGTCGCATGATGCATGAGAGCT", "k": 4, "canonical": true, "low": 2, "high": 6, "full": true }
}
```

**Response** (equals `jellyfish count -C -m 4` + `jellyfish histo -f -l 2 -h 6`, Jellyfish 2.3.1):
```json
{
  "spectrum": { "1": 16, "2": 2, "3": 1, "4": 1 },
  "histogram": [ {"bin":1,"frequency":16}, {"bin":2,"frequency":2}, {"bin":3,"frequency":1}, {"bin":4,"frequency":1},
                 {"bin":5,"frequency":0}, {"bin":6,"frequency":0}, {"bin":7,"frequency":0} ]
}
```
low 2 > increment 1, so base = 1 and the first bucket is labelled 1; the last bucket (7 = high + increment) is the cap.

## Performance

- **Time Complexity:** O(n) to build the k-mer table.
- **Space Complexity:** O(distinct k-mers).

## See Also

- [count_kmers](count_kmers.md) — raw k-mer → count map
- [analyze_kmers](analyze_kmers.md) — aggregate statistics
- [kmer_frequencies](kmer_frequencies.md) — normalized frequencies
