# calculate_similarity

Calculate similarity between two DNA sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Core |
| **Tool Name** | `calculate_similarity` |
| **Method ID** | `GenomicAnalyzer.CalculateSimilarity` |
| **Version** | 1.0.0 |

## Description

Calculates the exact k-mer Jaccard index J = |A∩B| / |A∪B| of the two sequences' distinct k-mer sets (Jaccard 1901; Ondov et al. 2016), reported as a **percentage** between 0 (no shared k-mer) and 100 (identical k-mer sets). Both k-mer sets empty (both sequences empty or shorter than k) returns 0 (convention). Delegates to `GenomicAnalyzer.CalculateSimilarity`, which calls the canonical `KmerAnalyzer.JaccardSimilarity`.

By default the literal (strand-specific) k-mers are compared (Python `set`, `scipy.spatial.distance.jaccard` on presence vectors). With `canonical: true` each k-mer is replaced by min(k-mer, reverse complement) as in Mash and sourmash (exact, unsketched), so a sequence and its reverse complement score 100.

**Estimators, distances, ANI.** `sketchSize` estimates the Jaccard index from a Mash bottom-s MinHash sketch (as `mash dist`), `scaled` from a sourmash FracMinHash sketch; the result also carries `jaccardDistance` (1 - J), `mashDistance`, `jaccardAni` and - for the exact and FracMinHash estimators - the containment indices and containment ANI (`C^(1/k)`), verified against sourmash 4.9.4 on 4290 rows. Output fields: `estimator`, `jaccardDistance`, `mashDistance`, `jaccardAni`, `containmentAInB`, `containmentBInA`, `maxContainment`, `containmentAniAInB`, `containmentAniBInA`, `sharedHashes`, `denominatorHashes`, `mashPValue` (the last six are absent when not applicable).

**Confidence intervals.** With `scaled`, `aniEstimates` holds the sourmash ANI results `containmentAInB`, `containmentBInA`, `maxContainment`, `jaccard`, each `{ distance, ani, probabilityNothingInCommon, probabilityExceedsThreshold, sizeIsInaccurate, jaccardError? }` and, with `estimateAniInterval: true`, `{ distanceLow, distanceHigh, aniLow, aniHigh }` (sourmash `estimate_ci`, Hera et al. 2023; null when the sketch is too small). Verified against sourmash 4.9.4 (1 120 rows) and mash 2.3 (400 rows).

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs#L707](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs#L699)

## Input/Output

**Input:** `{ sequence1: string, sequence2: string, kmerSize?: integer (default 5, >= 1), canonical?: boolean (default false), sketchSize?: integer (Mash bottom-s sketch), scaled?: integer (sourmash FracMinHash), estimateAniInterval?: boolean (default false), confidence?: number (default 0.95) }`

**Output:** `{ similarity: number (0-100 percentage scale), estimator, jaccardDistance, mashDistance, jaccardAni, containmentAInB?, containmentBInA?, maxContainment?, containmentAniAInB?, containmentAniBInA?, sharedHashes?, denominatorHashes?, mashPValue? }`

## Example

**User Prompt:**
> How similar are ATGCATGC and ATGCATTT (k = 5)?

**Tool Call:**
```json
{ "tool": "calculate_similarity", "arguments": { "sequence1": "ATGCATGC", "sequence2": "ATGCATTT" } }
```

**Response:**
```json
{ "similarity": 33.333333333333336 }
```
