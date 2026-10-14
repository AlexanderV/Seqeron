# codon_frequencies

Codon usage frequencies in a reading frame.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `codon_frequencies` |
| **Method ID** | `SequenceStatistics.CalculateCodonFrequencies` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes codon usage frequencies over the non-overlapping in-frame triplets of a DNA or RNA
sequence (case-insensitive; U is read as T and codons are reported in DNA spelling, as EMBOSS cusp), following the Kazusa CUTG definition `frequency = count(codon) / total
counted codons` (Nakamura et al. 2000). Reading starts at `readingFrame` (0, 1, or 2);
trailing bases that do not complete a triplet are ignored, and codons containing any
non-ACGT base are excluded from the count. Only observed codons appear as keys.

Optional `normalization` selects the EMBOSS 6.6.0 `cusp` column: `Total` (default, count/total),
`PerThousand` (count/total × 1000, cusp "Frequency") or `PerAminoAcid` (count / count of the codons
encoding the same amino acid, stop codons forming one `*` family, cusp "Fraction"). The synonymous
families follow the NCBI genetic code `genetic_code` (default 1 = Standard, as cusp; same contract as
the MolTools codon tools). Library overload: `CalculateCodonFrequencies(seq, frame, normalization, geneticCode)`.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L1830](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L1830)
- Evidence: `docs/Evidence/SEQ-CODON-FREQ-001-Evidence.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `dnaSequence` | string | Yes | DNA or RNA coding sequence, U read as T (min length 1) |
| `readingFrame` | integer | No | Reading frame 0, 1, or 2 (default 0) |
| `normalization` | string | No | `Total` (default), `PerThousand` or `PerAminoAcid` (case-insensitive exact name) |
| `genetic_code` | integer | No | NCBI translation table for `PerAminoAcid` families (default 1; one of 1-6, 9-16, 21-33; validated for every normalization) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `frequencies` | object | Map of codon → value: count/total (sums to 1), ×1000, or the per-amino-acid fraction (sums to 1 within each family) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | Reading frame must be 0, 1, or 2 |
| 1001 | Unknown normalization '…'. Expected one of: Total, PerThousand, PerAminoAcid. |
| 1001 | Unknown genetic code table: N. Supported NCBI tables: … |

## Examples

### Example 1: Frame 0

**User Prompt:**
> Codon frequencies of "ATGATGAAA" in frame 0.

**Expected Tool Call:**
```json
{
  "tool": "codon_frequencies",
  "arguments": { "dnaSequence": "ATGATGAAA", "readingFrame": 0 }
}
```

**Response:**
```json
{ "frequencies": { "ATG": 0.6666666666666666, "AAA": 0.3333333333333333 } }
```
Non-overlapping triplets ATG, ATG, AAA ⇒ ATG=2/3, AAA=1/3.

### Example 2: Frame 1

**Input:** `{ "dnaSequence": "ATGATGAAA", "readingFrame": 1 }`
→ **Response:** `{ "frequencies": { "TGA": 1.0 } }`
From index 1 the triplets are TGA, TGA; trailing `AA` is ignored.

### Example 3: cusp "Fraction" (per amino acid)

**Input:** `{ "dnaSequence": "ATGAAAAAGAAATAA", "normalization": "PerAminoAcid" }`
→ **Response:** `{ "frequencies": { "ATG": 1.0, "AAA": 0.6666666666666666, "AAG": 0.3333333333333333, "TAA": 1.0 } }`
Lys family AAA×2 + AAG×1. On the 26-codon cusp test sequence of SEQ-CODON-FREQ-001 the tool equals the
executed cusp table (GCT 0.400 / 76.923 per thousand, AAA 0.667, TAA/TAG/TGA 0.333); with
`"genetic_code": 2` AGA/AGG join the stop family (0.25 each) and TGA joins Trp.

## Performance

- **Time Complexity:** O(n). **Space Complexity:** O(distinct codons).

## See Also

- [dinucleotide_frequencies](dinucleotide_frequencies.md)
