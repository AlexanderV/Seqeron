# crispr_system_info

Return the PAM / guide-length metadata for a CRISPR nuclease system.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `crispr_system_info` |
| **Method ID** | `CrisprDesigner.GetSystem` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Looks up the built-in metadata record for a CRISPR system: its name, PAM recognition sequence (IUPAC), guide/spacer length, whether the PAM is 3′ of the target (Cas9-like) or 5′ of the target (Cas12a-like), the cut offset used by `design_guide_rnas` to place the cut, and a human description.

| System | PAM | Guide length | PAM after target | Cut offset |
|--------|-----|--------------|------------------|------------|
| SpCas9 | NGG | 20 | true | 3 |
| SpCas9_NAG | NAG | 20 | true | 3 |
| SpCas9_NGA | NGA | 20 | true | 3 |
| SaCas9 | NNGRRT | 21 | true | 3 |
| SaCas9_NNGRRN | NNGRRN | 21 | true | 3 |
| Cas12a / AsCas12a | TTTV | 23 | false | 18 |
| LbCas12a | TTTV | 23 | false | 18 |
| CasX | TTCN | 20 | false | 13 |

**Cut offset** (`cutOffset`) is the number of nucleotides between the PAM and the cut on the PAM-bearing strand: for PAM-after systems the cut lies that many bases 5′ of the PAM (SpCas9/SaCas9: 3 bp, the CRISPOR `cut = 3` convention), for PAM-before systems that many bases 3′ of the PAM (Cas12a: after the 18th protospacer base, CRISPOR; CasX: 13, the midpoint of the 12–14 nt non-target-strand cleavage window that Liu et al. 2019, *Nature* 566:218, report as "12–14 nucleotides after the PAM on the non-target strand and 22–25 nucleotides after the PAM on the target strand"; CasX makes a staggered cut, the tool reports the non-target-strand position).

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L26](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L26)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `system_type` | enum | Yes | `SpCas9` \| `SpCas9_NAG` \| `SpCas9_NGA` \| `SaCas9` \| `SaCas9_NNGRRN` \| `Cas12a` \| `AsCas12a` \| `LbCas12a` \| `CasX`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | System name. |
| `pamSequence` | string | PAM recognition sequence (IUPAC). |
| `guideLength` | integer | Guide/spacer length (nt). |
| `pamAfterTarget` | boolean | True if PAM is 3′ of the target. |
| `description` | string | Human description. |
| `cutOffset` | integer | Nucleotides between the PAM and the cut on the PAM-bearing strand (3 for the Cas9 systems, 18 for Cas12a, 13 for CasX). |

## Errors

| Code | Message |
|------|---------|
| 1001 | `Unknown CRISPR system: {system_type}` (`ArgumentException`; `system_type` is not one of the nine listed values, e.g. a numeric value outside the enum) |

## Examples

### Example 1: SpCas9 → NGG PAM, 20-nt guide, PAM after target.

### Example 2: Cas12a → TTTV PAM, 23-nt guide, PAM before target.

## See Also

- [find_pam_sites](find_pam_sites.md), [design_guide_rnas](design_guide_rnas.md)
