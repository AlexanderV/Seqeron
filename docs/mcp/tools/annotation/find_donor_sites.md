# find_donor_sites

Find candidate 5' (donor) splice sites.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_donor_sites` |
| **Method ID** | `SpliceSitePredictor.FindDonorSites` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans a sequence for 5' splice donor sites at every GU/GT dinucleotide (GC too when `includeNonCanonical`).
`position` is the 0-based index of the G of GT, i.e. the first intron nucleotide (+1).

- **`model = "shapiro-senapathy"` (default):** the Shapiro & Senapathy (1987, NAR 15:7155) 5' splice-site score
  of the 8-nt subsequence −2..+6: `(t − mint)/(maxt − mint)`, `t` = sum of the Table 1 percentages of the
  observed nucleotides in the chosen `category` subtable (default primate). This is the published 0–100 score
  divided by 100 (primate `AAGTGAGT` → 0.843). A candidate whose window runs off the sequence or contains a
  non-A/C/G/T/U character is not scored. `confidence` = clamp((score − 0.5)/0.5, 0, 1).
- **`model = "maxent"`:** MaxEntScan score5 (Yeo & Burge 2004) of the 9-nt window −3..+6 in bits (same as
  [maxent_score](maxent_score.md)); `minScore` is then a threshold in bits; `confidence` = 2^s/(1+2^s).
- **U12-type donors** (`includeNonCanonical`, Shapiro–Senapathy model only): every AT/AU or GT/GU start is also
  scored as the fraction of the 8 positions +1..+8 matching `/RTATCCTY` (Sharp & Burge 1997; Burge et al. 1999)
  and reported as `U12Donor`; a GT can thus be reported both as `Donor` and `U12Donor`.

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L445](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L348)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA/DNA sequence to scan |
| `minScore` | number | No | Minimum score (default 0.5): [0,1] for shapiro-senapathy / U12, bits for maxent |
| `includeNonCanonical` | boolean | No | Include GC donors and U12-type (/RTATCCTY) donors (default false) |
| `model` | string | No | `shapiro-senapathy` (default) or `maxent` |
| `category` | string | No | S&S Table 1 subtable: `primate` (default), `total`, `rodent`, `othermammal`, `othervertebrate`, `invertebrate`, `plant`, `viral` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `sites[].position` | integer | 0-based index of the G of GT/GC (first intron nt) |
| `sites[].type` | string | `Donor` (U2-type GT/GC) or `U12Donor` |
| `sites[].motif` | string | Context −3..+6 (U12: −3..+8), RNA alphabet, truncated at sequence ends |
| `sites[].score` | number | S&S score in [0,1], U12 match fraction, or MaxEntScan bits |
| `sites[].confidence` | number | See Description |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Unknown model / Unknown category |

## Examples

### Example 1: Perfect consensus

`CAGGUAAGU`: every scored position −2..+6 carries its highest primate percentage:

```json
{ "sites": [ { "position": 3, "type": "Donor", "motif": "CAGGUAAGU", "score": 1.0, "confidence": 1.0 } ] }
```

### Example 2: Shapiro & Senapathy worked example

`AAGTGAGT` → score 462/548 = 0.843065693430657 (paper: 84.3).

### Example 3: MaxEntScan model

`CAGGTAAGT`, `model = "maxent"` → score 10.858313101356437 bits.

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(k) for k sites

## See Also

- [find_acceptor_sites](find_acceptor_sites.md) — 3' acceptor sites
- [predict_introns](predict_introns.md) — pair donors and acceptors into introns
- [maxent_score](maxent_score.md) — MaxEntScan score5 donor score (bits)
