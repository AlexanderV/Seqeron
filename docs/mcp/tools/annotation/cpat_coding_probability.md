# cpat_coding_probability

CPAT coding probability of a transcript (Wang et al. 2013, CPAT 3.0.5).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `cpat_coding_probability` |
| **Method ID** | `GenomeAnnotator.PredictCodingProbability`, `GenomeAnnotator.PredictBestCodingOrf` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds CPAT ORF candidates (start ATG, stops TAG/TAA/TGA by default; ORF length in nt incl. stop must exceed
`minOrfLength`; unterminated ORFs run to the last complete codon; optional antisense; the `topOrfCount` longest are
kept), scores each with transcript length, ORF length, Fickett TESTCODE (CPAT table set) and hexamer score, and returns
the logistic probability `1/(1+e^-η)`, η = intercept + Σ coefficient·feature. `isCoding` = probability ≥ cutoff.
`bestIndex` points at the highest-probability candidate (or the longest with `bestByLength`), -1 when none.
The model (hexamer tables + coefficients) is supplied by the caller: CPAT's prebuilt models are GPL and not bundled.
CPAT Human: intercept −7.8284881840360567, mRNA −5.4831277888113823e-05, ORF 0.011187815405519514,
Fickett 3.4471695641643225, Hexamer 6.068275924110611, cutoff 0.364 (tables: CPAT `Human_Hexamer.tsv`).
Coordinates follow CPAT: 1-based; antisense ORFs have `start > end`.

## Core Documentation Reference

- Source: [GenomeAnnotator.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs) — `PredictCodingProbability`
- Algorithm: [Coding_Potential_Calculation.md §10](../../../algorithms/Extended_Annotation/Coding_Potential_Calculation.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Transcript DNA |
| `codingHexamerFrequencies` / `noncodingHexamerFrequencies` | object | Yes | 6-mer → frequency |
| `intercept`, `mrnaLengthCoefficient`, `orfLengthCoefficient`, `fickettCoefficient`, `hexamerCoefficient` | number | Yes | `coef(mylogit)` |
| `codingProbabilityCutoff` | number | Yes | in [0, 1] |
| `searchAntisense` | boolean | No | default false |
| `topOrfCount` | integer | No | default 5 |
| `minOrfLength` | integer | No | default 75 |
| `startCodons` / `stopCodons` | string[] | No | defaults ATG / TAG, TAA, TGA |
| `bestByLength` | boolean | No | default false (best by probability) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `candidates` | array | `{ strand, frame, start, end, orfLength, orfSequence, mrnaLength, fickettScore, hexamerScore, codingProbability, isCoding }` |
| `bestIndex` | integer | index into `candidates`, -1 if none |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | invalid coefficients / cutoff / options |

## See Also

- [coding_potential](coding_potential.md), [fickett_score](fickett_score.md)
