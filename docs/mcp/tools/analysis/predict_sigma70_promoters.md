# predict_sigma70_promoters

σ70 promoter free energy and transcription rate (Promoter Calculator v1.0).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `predict_sigma70_promoters` |
| **Method ID** | `MotifFinder.PredictSigma70Promoters` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Port of the Promoter Calculator v1.0 (La Fleur, Hossain & Salis 2022, Nat Commun 13:5159; reference code hsalis/SalisLabCode `Promoter_Calculator_v1_0.py`, `util.py`, `free_energy_coeffs.npy`, `model_intercept.npy`). For every TSS it enumerates UP (24) · 1 · −35 · spacer (15–20) · −10 · discriminator (6–10) · ITR (20) configurations, evaluates ΔG_total = ΔG_−10 + ΔG_−35 + ΔG_disc + ΔG_ITR + ΔG_ext−10 + ΔG_spacer + ΔG_UP + intercept (trained 3-mer / dimer energies, quadratic spacer term, groove width + rigidity of the UP region, DNA:DNA − RNA:DNA hybrid of the ITR) and keeps the minimum; transcription rate = 42·exp(−β·ΔG_total), β = 1.636217004872062 (E. coli MG1655) or 0.81632623 (`inVitro`). Needs ≥ 78 nt. `tss` follows the reference (minus strand: n − t of the reverse-complement TSS t, first transcribed base at tss − 1); segment starts are forward-strand coordinates. Order: `+` by TSS, then `-` by TSS. `best` = lowest ΔG_total.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/MotifFinder.Sigma70Promoters.cs#L149](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.Sigma70Promoters.cs#L149)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (A/C/G/T) |
| `bothStrands` | boolean | No | Also predict on the minus strand (default true) |
| `inVitro` | boolean | No | Use the in-vitro β (default false: E. coli) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array<object> | Per-TSS minimum-ΔG configuration `{strand, tss, promoterSequence, up, minus35, spacer, minus10, discriminator, itr, upStart, minus35Start, spacerStart, minus10Start, discriminatorStart, deltaGTotal, deltaG10, deltaG35, deltaGDiscriminator, deltaGItr, deltaGExtended10, deltaGSpacer, deltaGUp, deltaGBind, transcriptionRate}` |
| `best` | object|null | Item with the lowest deltaGTotal (null when no TSS qualifies) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: lacUV5 core promoter, plus strand

**User Prompt:**
> Predict σ70 promoter strength along this lacUV5 fragment (plus strand).

**Tool Call:**
```json
{
  "tool": "predict_sigma70_promoters",
  "arguments": {
    "sequence": "AGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTCA",
    "bothStrands": false
  }
}
```

**Response:**
```json
{
  "items": [
    {
      "strand": "+",
      "tss": 58,
      "promoterSequence": "AGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAAC",
      "up": "AGCTCACTCATTAGGCACCCCAGG",
      "minus35": "TTTACA",
      "spacer": "CTTTATGCTTCCGGC",
      "minus10": "TCGTAT",
      "discriminator": "AATGTG",
      "itr": "TGGAATTGTGAGCGGATAAC",
      "upStart": 0,
      "minus35Start": 25,
      "spacerStart": 31,
      "minus10Start": 46,
      "discriminatorStart": 52,
      "deltaGTotal": -1.6393329933233844,
      "deltaG10": -0.5932690220174014,
      "deltaG35": -1.095566862007599,
      "deltaGDiscriminator": -0.004857517483925856,
      "deltaGItr": 0.07615303241093338,
      "deltaGExtended10": 0.027096935239089858,
      "deltaGSpacer": 0.36700000000000443,
      "deltaGUp": 0.5424688942938479,
      "deltaGBind": -0.7522700544920581,
      "transcriptionRate": 613.9872353549823
    },
    {
      "strand": "+",
      "tss": 59,
      "promoterSequence": "AGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACA",
      "up": "AGCTCACTCATTAGGCACCCCAGG",
      "minus35": "TTTACA",
      "spacer": "CTTTATGCTTCCGGC",
      "minus10": "TCGTAT",
      "discriminator": "AATGTGT",
      "itr": "GGAATTGTGAGCGGATAACA",
      "upStart": 0,
      "minus35Start": 25,
      "spacerStart": 31,
      "minus10Start": 46,
      "discriminatorStart": 52,
      "deltaGTotal": -1.7501695850501884,
      "deltaG10": -0.5932690220174014,
      "deltaG35": -1.095566862007599,
      "deltaGDiscriminator": -0.004857517483925856,
      "deltaGItr": -0.03468355931587059,
      "deltaGExtended10": 0.027096935239089858,
      "deltaGSpacer": 0.36700000000000443,
      "deltaGUp": 0.5424688942938479,
      "deltaGBind": -0.7522700544920581,
      "transcriptionRate": 736.0712011570919
    },
    {
      "strand": "+",
      "tss": 60,
      "promoterSequence": "AGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAA",
      "up": "AGCTCACTCATTAGGCACCCCAGG",
      "minus35": "TTTACA",
      "spacer": "CTTTATGCTTCCGGC",
      "minus10": "TCGTAT",
      "discriminator": "AATGTGTG",
      "itr": "GAATTGTGAGCGGATAACAA",
      "upStart": 0,
      "minus35Start": 25,
      "spacerStart": 31,
      "minus10Start": 46,
      "discriminatorStart": 52,
      "deltaGTotal": -1.6789174903686717,
      "deltaG10": -0.5932690220174014,
      "deltaG35": -1.095566862007599,
      "deltaGDiscriminator": -0.004857517483925856,
      "deltaGItr": 0.0365685353656462,
      "deltaGExtended10": 0.027096935239089858,
      "deltaGSpacer": 0.36700000000000443,
      "deltaGUp": 0.5424688942938479,
      "deltaGBind": -0.7522700544920581,
      "transcriptionRate": 655.0705670600657
    },
    {
      "strand": "+",
      "tss": 61,
      "promoterSequence": "AGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAAT",
      "up": "AGCTCACTCATTAGGCACCCCAGG",
      "minus35": "TTTACA",
      "spacer": "CTTTATGCTTCCGGCTCG",
      "minus10": "TATAAT",
      "discriminator": "GTGTGG",
      "itr": "AATTGTGAGCGGATAACAAT",
      "upStart": 0,
      "minus35Start": 25,
      "spacerStart": 31,
      "minus10Start": 49,
      "discriminatorStart": 55,
      "deltaGTotal": -3.031803752017029,
      "deltaG10": -1.781524450740899,
      "deltaG35": -1.095566862007599,
      "deltaGDiscriminator": 0.023900602351315935,
      "deltaGItr": -0.07087509947156181,
      "deltaGExtended10": 0.19135161731618874,
      "deltaGSpacer": 0.116800000000012,
      "deltaGUp": 0.5424688942938479,
      "deltaGBind": -2.026470801138449,
      "transcriptionRate": 5993.066356723385
    },
    {
      "strand": "+",
      "tss": 62,
      "promoterSequence": "AGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATT",
      "up": "AGCTCACTCATTAGGCACCCCAGG",
      "minus35": "TTTACA",
      "spacer": "CTTTATGCTTCCGGCTCG",
      "minus10": "TATAAT",
      "discriminator": "GTGTGGA",
      "itr": "ATTGTGAGCGGATAACAATT",
      "upStart": 0,
      "minus35Start": 25,
      "spacerStart": 31,
      "minus10Start": 49,
      "discriminatorStart": 55,
      "deltaGTotal": -2.897216462063053,
      "deltaG10": -1.781524450740899,
      "deltaG35": -1.095566862007599,
      "deltaGDiscriminator": 0.023900602351315935,
      "deltaGItr": 0.06371219048241457,
      "deltaGExtended10": 0.19135161731618874,
      "deltaGSpacer": 0.116800000000012,
      "deltaGUp": 0.5424688942938479,
      "deltaGBind": -2.026470801138449,
      "transcriptionRate": 4808.519215543328
    },
    {
      "strand": "+",
      "tss": 63,
      "promoterSequence": "AGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTT",
      "up": "AGCTCACTCATTAGGCACCCCAGG",
      "minus35": "TTTACA",
      "spacer": "CTTTATGCTTCCGGCTCG",
      "minus10": "TATAAT",
      "discriminator": "GTGTGGAA",
      "itr": "TTGTGAGCGGATAACAATTT",
      "upStart": 0,
      "minus35Start": 25,
      "spacerStart": 31,
      "minus10Start": 49,
      "discriminatorStart": 55,
      "deltaGTotal": -3.031803752017029,
      "deltaG10": -1.781524450740899,
      "deltaG35": -1.095566862007599,
      "deltaGDiscriminator": 0.023900602351315935,
      "deltaGItr": -0.07087509947156186,
      "deltaGExtended10": 0.19135161731618874,
      "deltaGSpacer": 0.116800000000012,
      "deltaGUp": 0.5424688942938479,
      "deltaGBind": -2.026470801138449,
      "transcriptionRate": 5993.066356723385
    },
    {
      "strand": "+",
      "tss": 64,
      "promoterSequence": "AGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTC",
      "up": "AGCTCACTCATTAGGCACCCCAGG",
      "minus35": "TTTACA",
      "spacer": "CTTTATGCTTCCGGCTCG",
      "minus10": "TATAAT",
      "discriminator": "GTGTGGAAT",
      "itr": "TGTGAGCGGATAACAATTTC",
      "upStart": 0,
      "minus35Start": 25,
      "spacerStart": 31,
      "minus10Start": 49,
      "discriminatorStart": 55,
      "deltaGTotal": -2.8806286728250274,
      "deltaG10": -1.781524450740899,
      "deltaG35": -1.095566862007599,
      "deltaGDiscriminator": 0.023900602351315935,
      "deltaGItr": 0.08029997972043977,
      "deltaGExtended10": 0.19135161731618874,
      "deltaGSpacer": 0.116800000000012,
      "deltaGUp": 0.5424688942938479,
      "deltaGBind": -2.026470801138449,
      "transcriptionRate": 4679.765297193107
    },
    {
      "strand": "+",
      "tss": 65,
      "promoterSequence": "AGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTCA",
      "up": "AGCTCACTCATTAGGCACCCCAGG",
      "minus35": "TTTACA",
      "spacer": "CTTTATGCTTCCGGCTCG",
      "minus10": "TATAAT",
      "discriminator": "GTGTGGAATT",
      "itr": "GTGAGCGGATAACAATTTCA",
      "upStart": 0,
      "minus35Start": 25,
      "spacerStart": 31,
      "minus10Start": 49,
      "discriminatorStart": 55,
      "deltaGTotal": -3.044998584365458,
      "deltaG10": -1.781524450740899,
      "deltaG35": -1.095566862007599,
      "deltaGDiscriminator": 0.023900602351315935,
      "deltaGItr": -0.08406993181999087,
      "deltaGExtended10": 0.19135161731618874,
      "deltaGSpacer": 0.116800000000012,
      "deltaGUp": 0.5424688942938479,
      "deltaGBind": -2.026470801138449,
      "transcriptionRate": 6123.861140216731
    }
  ],
  "best": {
    "strand": "+",
    "tss": 65,
    "promoterSequence": "AGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTCA",
    "up": "AGCTCACTCATTAGGCACCCCAGG",
    "minus35": "TTTACA",
    "spacer": "CTTTATGCTTCCGGCTCG",
    "minus10": "TATAAT",
    "discriminator": "GTGTGGAATT",
    "itr": "GTGAGCGGATAACAATTTCA",
    "upStart": 0,
    "minus35Start": 25,
    "spacerStart": 31,
    "minus10Start": 49,
    "discriminatorStart": 55,
    "deltaGTotal": -3.044998584365458,
    "deltaG10": -1.781524450740899,
    "deltaG35": -1.095566862007599,
    "deltaGDiscriminator": 0.023900602351315935,
    "deltaGItr": -0.08406993181999087,
    "deltaGExtended10": 0.19135161731618874,
    "deltaGSpacer": 0.116800000000012,
    "deltaGUp": 0.5424688942938479,
    "deltaGBind": -2.026470801138449,
    "transcriptionRate": 6123.861140216731
  }
}
```

## Worked Example

The best configuration (TSS 65 here = 108 in the full lac region) is −35 TTTACA, 18-nt spacer, −10 TATAAT, ΔG_total −3.044998584365458, rate 6123.861140216731 — identical to the reference Python. 60 random sequences (4634 per-TSS predictions, both strands, both β) agree bit for bit in every field.

## See Also

- [find_sigma70_promoters](find_sigma70_promoters.md) — Consensus −35 / −10 pairing
- [find_promoter_elements_by_matrix](find_promoter_elements_by_matrix.md) — Eukaryotic Pol II promoter matrices

## References

- Algorithm source: [Seqeron.Genomics.Analysis/MotifFinder.Sigma70Promoters.cs#L149](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.Sigma70Promoters.cs#L149)
- Binding: [AnalysisTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs)
