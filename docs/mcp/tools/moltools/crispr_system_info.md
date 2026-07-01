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

Looks up the built-in metadata record for a CRISPR system: its name, PAM recognition sequence (IUPAC), guide/spacer length, whether the PAM is 3′ of the target (Cas9-like) or 5′ of the target (Cas12a-like), and a human description.

| System | PAM | Guide length | PAM after target |
|--------|-----|--------------|------------------|
| SpCas9 | NGG | 20 | true |
| SpCas9_NAG | NAG | 20 | true |
| SaCas9 | NNGRRT | 21 | true |
| Cas12a / AsCas12a | TTTV | 23 | false |
| LbCas12a | TTTV | 24 | false |
| CasX | TTCN | 20 | false |

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L18](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L18)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `system_type` | enum | Yes | `SpCas9` \| `SpCas9_NAG` \| `SaCas9` \| `Cas12a` \| `AsCas12a` \| `LbCas12a` \| `CasX`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | System name. |
| `pamSequence` | string | PAM recognition sequence (IUPAC). |
| `guideLength` | integer | Guide/spacer length (nt). |
| `pamAfterTarget` | boolean | True if PAM is 3′ of the target. |
| `description` | string | Human description. |

## Errors

_None._ (Enum inputs are always valid.)

## Examples

### Example 1: SpCas9 → NGG PAM, 20-nt guide, PAM after target.

### Example 2: Cas12a → TTTV PAM, 23-nt guide, PAM before target.

## See Also

- [find_pam_sites](find_pam_sites.md), [design_guide_rnas](design_guide_rnas.md)
