# predict_signal_peptide

von Heijne signal-peptide cleavage-site prediction.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `predict_signal_peptide` |
| **Method ID** | `ProteinMotifFinder.PredictSignalPeptide` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Predicts the **signal-peptide cleavage site** using the von Heijne (1986)
weight-matrix method, as implemented by EMBOSS `sigcleave`. Returns the best-scoring
cleavage site (1-based mature-peptide start), the weight-matrix score, the signal
sequence, the 15-residue scoring window, and whether the site is a likely signal
peptide. A separate prokaryotic weight matrix is used when `prokaryote` is true.

## Core Documentation Reference

- Source: [ProteinMotifFinder.cs#L606](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs#L606)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence (min length 1) |
| `prokaryote` | boolean | No | Use the prokaryotic matrix (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `found` | boolean | Whether a candidate cleavage site was found |
| `cleavagePosition` | integer | 1-based mature-peptide start |
| `score` | number | von Heijne weight-matrix score |
| `signalSequence` | string | Predicted signal peptide (residues before the mature start) |
| `windowSequence` | string | 15-residue scoring window |
| `isLikelySignalPeptide` | boolean | Whether the site passes the likelihood threshold |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: ACH2_DROME (EMBOSS sigcleave reference)

**User Prompt:**
> Predict the signal peptide of the ACH2_DROME precursor.

**Expected Tool Call:**
```json
{
  "tool": "predict_signal_peptide",
  "arguments": { "proteinSequence": "MAPGCCTTRPRPIALLAHIWRHCKPLCLLLVLLLLCETVQANP…" }
}
```

**Response:**
```json
{ "found": true, "cleavagePosition": 42, "score": 13.7390400704164, "windowSequence": "LLVLLLLCETVQANP", "isLikelySignalPeptide": true }
```
EMBOSS sigcleave reports maximum score 13.739 with the mature peptide starting at
residue 42.

### Example 2: Short hydrophobic N-terminus

**User Prompt:**
> Predict a signal peptide for "AAAAAAAAAAAAGAN".

**Expected Tool Call:**
```json
{
  "tool": "predict_signal_peptide",
  "arguments": { "proteinSequence": "AAAAAAAAAAAAGAN" }
}
```

**Response (fields):**
```json
{ "found": true }
```

## Performance

- **Time Complexity:** O(n) over cleavage-site windows.
- **Space Complexity:** O(1).

## See Also

- [predict_transmembrane_helices](predict_transmembrane_helices.md) — hydropathy-based TM prediction
- [hydrophobicity_profile](hydrophobicity_profile.md) — Kyte-Doolittle profile
