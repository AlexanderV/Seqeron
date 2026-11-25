# fold_mirna_duplex

Minimum-free-energy miRNA·target duplex (RNAduplex / duplexfold model, Turner 2004): structure, base pairs and free energy.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `fold_mirna_duplex` |
| **Method ID** | `MiRnaAnalyzer.FoldMiRnaDuplex` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds the single contiguous miRNA·target duplex of minimum free energy (RNAduplex model of ViennaRNA), allowing internal loops and bulges of up to `maxLoopSize` unpaired nucleotides per side, scored with Turner 2004 nearest-neighbour parameters (kcal/mol, 37 °C). The structure uses RNAduplex notation: the miRNA span as `(`/`.`, then `&`, then the target span as `)`/`.`. Unlike [align_mirna_to_target](align_mirna_to_target.md) (fixed ungapped register) the fold recovers bulged/looped helices. `found = false` when no base pair is possible. DNA `T` is read as `U`.

## Core Documentation Reference

- Source: [MiRnaAnalyzer.cs#L690](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs#L690)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `miRnaSequence` | string | Yes | miRNA nucleotide sequence |
| `targetSequence` | string | Yes | Target nucleotide sequence |
| `maxLoopSize` | integer | No | Maximum unpaired nt per side of an internal loop / bulge (default 30) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `found` | boolean | False when no base pair is possible |
| `duplex.miRnaSequence / targetSequence` | string | Normalised sequences (upper case, T→U) |
| `duplex.miRnaStart / miRnaEnd` | integer | 0-based span of the duplex on the miRNA |
| `duplex.targetStart / targetEnd` | integer | 0-based span of the duplex on the target |
| `duplex.structure` | string | RNAduplex notation `((((&))))` |
| `duplex.basePairs[]` | object | `{miRnaIndex, targetIndex}` 0-based, outermost first |
| `duplex.freeEnergy` | number | Turner 2004 ΔG°37 (kcal/mol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | miRNA sequence cannot be null or empty |
| 1001 | Target sequence cannot be null or empty |
| 1001 | maxLoopSize must be >= 0 |

## Examples

### Example 1: Perfect let-7a duplex

`UGAGGUAGUAGGUUGUAUAGUU` against its reverse complement folds to a 22-bp helix, ΔG° −35.90 kcal/mol (ViennaRNA 2.8 duplexfold −35.6).

**Request:**
```json
{"miRnaSequence": "UGAGGUAGUAGGUUGUAUAGUU", "targetSequence": "AACUAUACAACCUACUACCUCA"}
```

**Response (abridged):**
```json
{"found": true, "duplex": {"structure": "((((((((((((((((((((((&))))))))))))))))))))))", "freeEnergy": -35.9}}
```

## Performance

- **Time Complexity:** O(n·m·L²) for miRNA length n, target length m and loop bound L
- **Space Complexity:** O(n·m)

## See Also

- [align_mirna_to_target](align_mirna_to_target.md)
- [find_mirna_target_sites](find_mirna_target_sites.md)
