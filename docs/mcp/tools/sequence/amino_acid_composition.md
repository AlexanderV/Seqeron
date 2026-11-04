# amino_acid_composition

Calculate amino acid composition of a protein sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `amino_acid_composition` |
| **Method ID** | `SequenceStatistics.CalculateAminoAcidComposition` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the amino acid composition of a protein sequence, including counts of each amino acid, molecular weight, isoelectric point (pI), hydrophobicity, charged residue ratio, and aromatic residue ratio.

The charged (D+E+H+K+R) and aromatic (F+W+Y) ratios divide by the residue count by default
(`denominator = "Residues"`, EMBOSS pepstats Mole%: "MKW*" → 1/3). `denominator = "SequenceLength"`
divides by the full input length, `*` included (Biopython 1.88 `ProteinAnalysis.aromaticity()`:
"MKW*" → 0.25, "mkwfy**" → 0.4285714285714286). `includeFractions = true` adds `fractions`, the
per-residue count/denominator of every occurring residue letter (library
`CalculateAminoAcidFractions`; pepstats Mole%/100 or Biopython `amino_acids_percent`/100). `length`,
`counts`, MW, pI and GRAVY do not depend on the denominator.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L133](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L133)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The protein sequence to analyze (min length: 1) |
| `denominator` | string | No | `Residues` (default, EMBOSS pepstats) or `SequenceLength` (Biopython `len(seq)`), case-insensitive exact name |
| `includeFractions` | boolean | No | Also return `fractions` (default false: field omitted) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `length` | integer | Number of amino acids |
| `counts` | object | Count of each amino acid |
| `molecularWeight` | number | Molecular weight in Daltons |
| `isoelectricPoint` | number | Isoelectric point (pI), 0-14 |
| `hydrophobicity` | number | Average hydrophobicity |
| `chargedResidueRatio` | number | Ratio of charged residues (D, E, K, R, H) |
| `aromaticResidueRatio` | number | Ratio of aromatic residues (F, Y, W) |
| `fractions` | object | Only with `includeFractions`: residue → count/denominator |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 2002 | Invalid protein sequence |
| 1001 | Unknown composition denominator '…'. Expected one of: Residues, SequenceLength. |

## Examples

### Example 1: Protein composition

**User Prompt:**
> What's the amino acid composition of "MAEGEITTFT"?

**Expected Tool Call:**
```json
{
  "tool": "amino_acid_composition",
  "arguments": {
    "sequence": "MAEGEITTFT"
  }
}
```

**Response:**
```json
{
  "length": 10,
  "counts": {"M": 1, "A": 1, "E": 2, "G": 1, "I": 1, "T": 3, "F": 1},
  "molecularWeight": 1099.2116,
  "isoelectricPoint": 3.61,
  "hydrophobicity": 0.15000000000000002,
  "chargedResidueRatio": 0.2,
  "aromaticResidueRatio": 0.1
}
```

### Example 2: Biopython denominator with fractions

**Input:** `{ "sequence": "MKW*", "denominator": "SequenceLength", "includeFractions": true }`
→ **Response:** `{ "length": 3, "counts": {"M": 1, "K": 1, "W": 1}, "molecularWeight": 463.5935, "isoelectricPoint": 9.15, "hydrophobicity": -0.9666666666666667, "chargedResidueRatio": 0.25, "aromaticResidueRatio": 0.25, "fractions": {"M": 0.25, "K": 0.25, "W": 0.25} }`
Biopython 1.88 `ProteinAnalysis("MKW*")`: aromaticity 0.25, `amino_acids_percent` 25.0 for M/K/W.
With the default denominator the ratios are 1/3 (pepstats Mole% 33.333, Residues = 3).

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1) for standard amino acids

## See Also

- [protein_validate](protein_validate.md) - Validate protein sequences
- [nucleotide_composition](nucleotide_composition.md) - DNA/RNA composition
