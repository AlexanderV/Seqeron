# count_kmers_both_strands

Count k-mers over both strands of double-stranded DNA.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `count_kmers_both_strands` |
| **Method ID** | `KmerAnalyzer.CountKmersBothStrands` |
| **Version** | 1.1.0 |
| **Stability** | Stable |

## Description

Counts every k-mer over **both** strands of double-stranded DNA. The reported count of a
k-mer `w` is its overlapping occurrences on the forward strand plus its overlapping
occurrences on the reverse-complement strand, i.e. `count[w] = forward[w] + forward[RC(w)]`
(kPAL "balance" operation; generalized second Chargaff rule). Every observed k-mer keeps its
own key — `w` and `RC(w)` carry equal counts. The total over all k-mers is `2·(L − k + 1)`.
The input must be a valid DNA sequence.

With `acgtOnly: true` the tool follows kPAL (`kpal/klib.py`: `Profile.from_sequences` splits each sequence on
`[^AaCcGgTt]`, `balance()` adds each k-mer's count to its reverse complement's): sequences with N or IUPAC symbols are
accepted and every window containing a non-ACGT symbol is skipped on both strands
(`KmerAnalyzer.CountKmersBothStrands(sequence, k, KmerCountingOptions)`). Strand-collapsed canonical counts
(Jellyfish `count -C`, palindromes counted once) are `count_kmers` with `canonical: true`; the library rejects
`Canonical` for both-strand counting.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L3023](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L3023)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `acgtOnly` | boolean | No | kPAL rule: accept any symbols, skip windows containing a non-ACGT symbol. Default false (sequence must be A/C/G/T) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `counts` | object | Map of k-mer → summed forward + reverse-complement count |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: Homopolymer

**User Prompt:**
> Count 2-mers on both strands of "AAAA".

**Expected Tool Call:**
```json
{ "tool": "count_kmers_both_strands", "arguments": { "sequence": "AAAA", "k": 2 } }
```

**Response:**
```json
{ "counts": { "AA": 3, "TT": 3 } }
```
Forward AA=3; reverse complement "TTTT" contributes TT=3. Total 6 = 2·(L−k+1).

### Example 2: Mixed 2-mers

**Input:** `{ "sequence": "ATGC", "k": 2 }`
→ **Response:** `{ "counts": { "AT": 2, "TG": 1, "GC": 2, "CA": 1 } }`

Forward AT,TG,GC; reverse complement "GCAT" adds GC,CA,AT. AT and GC each reach 2.

### Example 3: kPAL ACGT-only counting (k=3)

**Expected Tool Call:**
```json
{
  "tool": "count_kmers_both_strands",
  "arguments": { "sequence": "ACGTNACGTAAcgtRTT", "k": 3, "acgtOnly": true }
}
```

**Response:**
```json
{ "counts": { "ACG": 6, "CGT": 6, "AAC": 1, "GTA": 1, "GTT": 1, "TAA": 1, "TAC": 1, "TTA": 1 } }
```
kPAL (run from source): `Profile.from_sequences(["ACGTNACGTAAcgtRTT"], 3)` then `balance()` gives the same profile (Σ 18 = 2 × 9 ACGT windows).

## Performance

- **Time Complexity:** O(n). **Space Complexity:** O(distinct k-mers).

## See Also

- [count_kmers](count_kmers.md)
- [kmer_frequencies](kmer_frequencies.md)
