# score_target_site_context_plus_plus

TargetScan 7 context++ score (Agarwal et al. 2015) of an 8mer, 7mer-m8, 7mer-A1 or 6mer site with every feature contribution.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `score_target_site_context_plus_plus` |
| **Method ID** | `MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Port of `targetscan_70_context_scores.pl` with the `Agarwal_2015_parameters.txt` coefficients (Agarwal et al. 2015, eLife 4:e05005). `mRnaSequence` is the 3′UTR; `siteStart`/`siteEnd` are 0-based inclusive and the span must match the site type (8 / 7 / 7 / 6 nt). The result lists the contribution of each feature (intercept, local AU, sRNA1/sRNA8 identity, Site8, SA, 3P pairing, Min_dist, Len_3UTR, Off6m, SPS, TA, Len_ORF, ORF8m, PCT), the unrounded `contextScorePartial` and `contextScore` as TargetScan reports it (3-decimal contributions, per-site-type ceiling −0.03 / −0.02 / −0.01 / 0, `null` when the site starts within 14 nt of the CDS). SPS and TA default to the bundled TargetScan 7.0 seed-region table. All optional inputs are sparse: `orfSequence` derives Len_ORF and ORF8m as `targetscan_count_8mers.pl` does (`orfLength` / `orf8mCount` override it); `plfoldUnpairedProbability` is the exact RNAplfold value (`-L 40 -W 80 -u 20`, row utrStart+7, column 14) that replaces the library's single-window SA estimate; `conservation` takes a Newick tree **or** a bundled TargetScan bin tree (`targetScanBin` 1–10) plus `speciesWithSite` (and optionally a PCT sigmoid row) and gives the Friedman 2009 branch length and PCT. A feature without data contributes 0 and is listed in `omittedFeatures`; because such a score is partial, the library treats it as a guarded limitation (MIRNA-TARGET-001) and throws `SeqeronLimitationException` outside Permissive mode unless every optional feature is supplied.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L1137](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L1137)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `mRnaSequence` | string | Yes | 3′UTR nucleotide sequence (T→U applied internally) |
| `miRnaSequence` | string | Yes | miRNA nucleotide sequence (≥ 8 nt) |
| `siteStart` | integer | Yes | 0-based inclusive site start |
| `siteEnd` | integer | Yes | 0-based inclusive site end (span 8 for Seed8mer, 7 for Seed7merM8/Seed7merA1, 6 for Seed6mer) |
| `siteType` | string | Yes | Seed8mer, Seed7merM8, Seed7merA1 or Seed6mer |
| `miRnaName` | string | No | Optional miRNA name (label only) |
| `sps` | number | No | Optional seed-pairing stability SPS (kcal/mol); overrides the bundled TargetScan value |
| `ta` | number | No | Optional target-site abundance TA (log10 site count; see compute_ta_3utr); overrides the bundled value |
| `orfLength` | number | No | Optional ORF length in nt (> 0); overrides the length derived from orfSequence |
| `orf8mCount` | integer | No | Optional count of cognate 8mer sites in the ORF (≥ 0); overrides the count derived from orfSequence |
| `orfSequence` | string | No | Optional ORF sequence (RNA/DNA, any case, '-'/'.' gaps accepted) |
| `plfoldUnpairedProbability` | number | No | Optional exact RNAplfold unpaired probability in [0, 1] of the site's 14-nt window; 0/omitted ⇒ library estimate |
| `conservation` | object | No | Optional multi-species conservation: `newick` OR `targetScanBin` (1–10), `speciesWithSite` (leaf names), optional `sigmoidParameters` {b0,b1,b2,b3} |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | object | Contributions: `siteType`, `intercept`, `localAuContribution`, `sRna1Contribution`, `sRna8Contribution`, `site8Contribution`, `saContribution`, `threePrimePairingContribution`, `minDistContribution`, `len3UtrContribution`, `off6mContribution`, `spsContribution`, `taContribution`, `lenOrfContribution`, `orf8mContribution`, `pctContribution`, `branchLengthScore`, `pct`, `pctIsNa`, `contextScorePartial`, `contextScore` (number | null), `omittedFeatures` (string[]) |

## Errors

| Code | Message |
|------|---------|
| 1001 | mRNA sequence cannot be null or empty |
| 1001 | miRNA sequence must have at least 8 nt (seed region nt 2-8) |
| 1001 | Site indices [s, e] are out of range for sequence of length L |
| 1001 | Unknown or unsupported siteType (valid values listed) |
| 1001 | A <type> site spans N nt; got [s, e] |
| 1001 | orfLength / orf8mCount / plfoldUnpairedProbability / sps / ta out of range |
| 1001 | conservation requires exactly one of newick or targetScanBin |
| 1001 | Invalid Newick tree |
| 1001 | SeqeronLimitationException MIRNA-TARGET-001 (partial score outside Permissive mode) |

## Examples

### Example 1: Worked example, 8mer (TargetScan 7.0 perl run)

181-nt 3′UTR of the TargetScan 7.0 reference run with let-7a: the 8mer at UTR 31–38 scores −0.742 (6mer 82–87: −0.118; 7mer-A1 125–131: −0.257).

**Request:**
```json
{"mRnaSequence": "GCUAAAGACAAUUACAUAACAUACACGUCACUACCUCAGCACGAAACUUGUUGGCCCAGUGUGAAUCGCUUAAGGGUUGCGUACCUCGAAGUAAGUGUGAUGCAUACGCCUUUACUUGCUGUGGUACCUCAUCCACCCCAUCGGACUGGCAUUUUUAUUACACUCAGAAACAGAACUCGGG", "miRnaSequence": "UGAGGUAGUAGGUUGUAUAGUU", "siteStart": 30, "siteEnd": 37, "siteType": "Seed8mer", "sps": -9.25, "ta": 3.393, "orfLength": 1500, "orf8mCount": 1, "plfoldUnpairedProbability": 0.009772092}
```

**Response (abridged):**
```json
{"score": {"siteType": "Seed8mer", "intercept": -0.589, "contextScore": -0.742}}
```

## Performance

- **Time Complexity:** O(n + m) plus one 80-nt McCaskill partition function when SA is estimated
- **Space Complexity:** O(n)

## See Also

- [score_target_site_grimson](score_target_site_grimson.md)
- [analyze_target_context](analyze_target_context.md)
- [compute_ta_3utr](compute_ta_3utr.md)
- [score_targetscan_conservation](score_targetscan_conservation.md)
- [find_mirna_target_sites](find_mirna_target_sites.md)
