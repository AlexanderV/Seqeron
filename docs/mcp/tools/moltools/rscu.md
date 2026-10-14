# rscu

Relative Synonymous Codon Usage per codon.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `rscu` |
| **Method ID** | `CodonUsageAnalyzer.CalculateRscu` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes RSCU (Sharp, Tuohy & Mosurski 1986) for each codon: `RSCU = observed / (total_synonymous / n_synonymous)`, i.e. observed count divided by the count expected if all synonymous codons of that amino acid were used equally. `RSCU = 1` means no bias, `> 1` over-represented, `< 1` under-represented. Codons of amino acids not present in the sequence get RSCU 0 (CodonW).

Optional `pseudocount` `c` (default 0) is added to every codon count before the per-family normalisation: `RSCU = n·(x + c) / (Σx + n·c)` — the formula of codonbias 0.5.0 `RelativeSynonymousCodonUsage(pseudocount=c)` and cubar `est_rscu(pseudo_cnt = c, level = "amino_acid")` (both default to `c = 1`). With `c > 0` an absent amino acid is uniform (RSCU 1 for each codon).

Optional `family_level` (default `AminoAcid`) selects the synonymous grouping. `Subfamily` reproduces the R package cubar's default `est_rscu(level = "subfam")`: every codon is labelled `subfam = paste(amino_acid, substr(codon, 1, 2), sep = "_")`, so an amino acid whose codons differ in the first two bases splits (table 1: Leu → Leu_TT {TTA, TTG} + Leu_CT {CTN}; Ser → Ser_TC + Ser_AG; Arg → Arg_CG + Arg_AG; stops → \*_TA {TAA, TAG} + \*_TG {TGA}), and RSCU is normalised within each sub-family. With `pseudocount = 0` an absent sub-family reports 0 (cubar: NaN).

## Core Documentation Reference

- Source: [CodonUsageAnalyzer.cs#L222](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs#L222)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Coding DNA sequence (frame 0), non-empty. |
| `genetic_code` | integer | No | NCBI translation table number (default `1` = Standard; e.g. `2` vertebrate mitochondrial, `11` bacterial/plastid). Defines the synonymous families and stops (CodonW `-code`). One of 1–6, 9–16, 21–33. |
| `pseudocount` | number | No | Pseudocount added to every codon count, finite and ≥ 0 (default `0` = none, CodonW). |
| `family_level` | string | No | `AminoAcid` (default; Sharp & Li / CodonW families) or `Subfamily` (cubar `level = "subfam"`: amino acid + first two codon bases). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `rscu` | object | Codon → RSCU value. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Unknown genetic code table: N. Supported NCBI tables: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |
| 1003 | pseudocount must be a finite non-negative number (got X). |
| 1004 | Unknown RSCU family level: X. Supported: AminoAcid, Subfamily. |

## Examples

### Example 1: Biased usage

`TTTTTTTTC` → TTT×2, TTC×1; Phe expected = 3/2 = 1.5 → `RSCU[TTT] = 1.333`, `RSCU[TTC] = 0.667`.

### Example 2: Unbiased usage

`TTTTTC` → TTT×1, TTC×1 → `RSCU = 1.0` for each.

### Example 3: Non-standard genetic code (`genetic_code = 2`)

`ATAATGATGTGATGGTGGAGATAA` (ATA×1, ATG×2, TGA×1, TGG×2, AGA×1, TAA×1). Under NCBI table 2 (vertebrate mitochondrial) Met = {ATA, ATG}, Trp = {TGA, TGG}, stops = {TAA, TAG, AGA, AGG}: `RSCU[ATA] = 0.667`, `RSCU[ATG] = 1.333`, `RSCU[TGA] = 0.667`, `RSCU[TGG] = 1.333`, `RSCU[AGA] = RSCU[TAA] = 2.0`. With the default table 1 the same gene gives `RSCU[ATA] = 3.0` (Ile), `RSCU[ATG] = 1.0`, `RSCU[TGA] = 1.5` (stop), `RSCU[AGA] = 6.0` (Arg). Families from Biopython `CodonTable.unambiguous_dna_by_id[2]`.

### Example 4: Pseudocount (`pseudocount = 1`)

`TTTTTTTTC` (TTT×2, TTC×1): `RSCU[TTT] = 2·(2+1)/(3+2) = 1.2`, `RSCU[TTC] = 2·(1+1)/(3+2) = 0.8`; Met (absent) `= 1·1/1 = 1.0`. Identical to codonbias 0.5.0 `RelativeSynonymousCodonUsage(pseudocount=1, ignore_stop=False)`.

### Example 5: cubar sub-families (`family_level = "Subfamily"`)

`CTGCTGCTGCTATTA` (CTG×3, CTA×1, TTA×1): Leu_CT = {CTT, CTC, CTA, CTG} → `RSCU[CTG] = 4·3/4 = 3.0`, `RSCU[CTA] = 1.0`, `RSCU[CTT] = RSCU[CTC] = 0`; Leu_TT = {TTA, TTG} → `RSCU[TTA] = 2·1/1 = 2.0`, `RSCU[TTG] = 0`. At the default amino-acid level the six-fold Leu family gives `RSCU[CTG] = 6·3/5 = 3.6`, `RSCU[TTA] = 1.2`.

## See Also

- [count_codons](count_codons.md), [codon_usage_statistics](codon_usage_statistics.md), [codon_adaptation_index](codon_adaptation_index.md)
