# effective_number_of_codons

Wright's Effective Number of Codons (ENC / Nc).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `effective_number_of_codons` |
| **Method ID** | `CodonUsageAnalyzer.CalculateEnc` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes Wright's Effective Number of Codons (Nc), a single measure of how far a gene departs from uniform synonymous-codon usage. It uses per-amino-acid homozygosity `F̂` grouped by degeneracy class (Wright Eq. 1), class averaging (Eq. 4), the isoleucine Eq. 5a fallback, and `Nc = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆` (Eq. 3). The result lies in `[20, 61]` under the standard code (20 = extreme bias, 61 = no bias; overshoot re-adjusted to the code's sense-codon count). When Nc cannot be calculated (a synonymous class other than the single 3-fold Ile class has no estimable amino acid — CodonW prints `*****`) the tool returns `enc = 0` and `calculable = false` (`CodonUsageAnalyzer.TryCalculateEnc`).

Optional `method` selects another estimator (`TryCalculateEnc(string, GeneticCode, EncOptions, out double)`), both as implemented by codonbias 0.5.0 `EffectiveNumberOfCodons`: `Novembre2002` — Novembre's background-corrected Nc′ (per amino acid χ² of the codon frequencies against those expected from the gene's A/C/G/T composition, or that of `background`; `F = (N·χ² + N − k)/(k(N − 1))`; unweighted class means; no pseudocount) — and `SunYangXia2013` (pseudocount 1; `F = (χ² + 1)/k`, i.e. `Σp²` without background correction; class means weighted by amino-acid counts). For these, a class with no estimable amino acid is imputed with `F̄ = 1/z` (a missing 3-fold class with `(F̄₂ + F̄₄)/2`), so the value is always calculable when the gene has a sense codon; the result is re-adjusted down to the code's sense-codon count, with no lower clamp. `codonw_cap = true` (Wright only) re-adjusts to CodonW's fixed 61 (`enc_out` prints any larger value as `61.00`, whatever `-code`) instead of the code's sense-codon count (60 for table 2, 62 for tables 3/4/5/9/10, 63 for table 6).

`EncPrimeNc` / `EncPrimeNcPrime` reproduce the `Nc` / `Ncp` columns of ENCprime, Novembre's reference program (github.com/jnovembre/ENCprime, `src/ENCprime.c`): F (Wright Eq. 1) or F′ = `(χ² + n − k)/(k(n − 1))` only for amino acids observed **n > 5** times; expected codon frequencies from the A/C/G/T composition of the gene (as `SeqCount -n` on the same FASTA) or of `background` (Nc′ only), no pseudocount; unweighted class means; a missing class is imputed with `1/k` (3-fold: `(F̄₂ + F̄₄)/2` when both exist); each class mean is floored at `1/k`; no cap. Always calculable for a non-empty sequence (a gene with no amino acid seen more than 5 times reports the uninformative Σ K_z·z, e.g. 61). Agrees with the compiled program at its printed 4 decimals on 397 genes × 16 NCBI tables; for table 13 ENCprime's `get_code` drops V/A/D/E (n_aa bug) and the library keeps them.

## Core Documentation Reference

- Source: [CodonUsageAnalyzer.cs#L485](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs#L485) (`TryCalculateEnc(string, GeneticCode, out double)`; `CalculateEnc` delegates to it); [CodonUsageAnalyzer.cs#L580](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs#L580) (`TryCalculateEnc(string, GeneticCode, EncOptions, out double)`, used when `method`/`codonw_cap` are given)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Coding DNA sequence (frame 0), non-empty. |
| `genetic_code` | integer | No | NCBI translation table number (default `1` = Standard; e.g. `2` vertebrate mitochondrial, `11` bacterial/plastid). Defines the synonymous families and stops (CodonW `-code`). One of 1–6, 9–16, 21–33. |
| `method` | string | No | `Wright1990` (default; CodonW `enc_out`), `Novembre2002` or `SunYangXia2013` (codonbias 0.5.0), `EncPrimeNc` or `EncPrimeNcPrime` (ENCprime program `Nc` / `Ncp`). |
| `pseudocount` | number | No | Novembre2002/SunYangXia2013 only: pseudocount added to every sense-codon count, ≥ 0 (omitted = 0 for Novembre2002, 1 for SunYangXia2013). |
| `background_correction` | boolean | No | Novembre2002/SunYangXia2013 only: expected codon frequencies from the background base composition (omitted = true for Novembre2002, false for SunYangXia2013). |
| `background` | string | No | With background correction, or `EncPrimeNcPrime`: sequence whose A/C/G/T composition is the background (omitted = the gene itself; ENCprime: the `SeqCount -n` input). |
| `codonw_cap` | boolean | No | Wright1990 only: cap at CodonW's fixed 61 instead of the code's sense-codon count (default `false`). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `enc` | number | Effective Number of Codons, 20..61 (Wright, table 1); 0 when not calculable. |
| `calculable` | boolean | `false` when Nc cannot be calculated (CodonW `*****`), otherwise `true`. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Unknown genetic code table: N. Supported NCBI tables: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |
| 1003 | Unknown ENC method: X. Supported: Wright1990, Novembre2002, SunYangXia2013, EncPrimeNc, EncPrimeNcPrime. / pseudocount, background_correction and background apply to the Novembre2002 and SunYangXia2013 methods only. / pseudocount and background_correction do not apply to the ENCprime methods (fixed program rules). |
| 1004 | codonw_cap applies to the Wright1990 method only. |
| 1005 | pseudocount must be a finite non-negative number (got X). |
| 1006 | background must be a non-empty sequence and requires background correction. / background applies to EncPrimeNcPrime only (ENCprime's Nc ignores the base composition). |

## Examples

### Example 1: Fully-populated biased gene

With `F̂₂=0.6, F̂₃=0.2667, F̂₄=0.4333, F̂₆=0.3333`, `Nc = 2 + 9/0.6 + 1/0.2667 + 5/0.4333 + 3/0.3333 = 41.288461538461526` (value verified against an independent Wright/codonW reference).

### Example 2: Non-standard genetic code (`genetic_code = 2`)

The same gene under NCBI table 2 (vertebrate mitochondrial): Nc = 45.0 (CodonW 1.4.4 `-enc -code 1` → 45.00); table 3 → 43.185, table 9 → 43.833 (CodonW 43.19 / 43.83).

### Example 3: Nc not calculable

`ATGAAAGAGCTGTTCGCCAAA`: the 3-, 4- and 6-fold classes have no estimable amino acid; CodonW 1.4.4 `-enc` prints `*****`. Output: `{ "enc": 0, "calculable": false }`.

### Example 4: Novembre 2002 / Sun–Yang–Xia 2013

Gene of Example 1: `method = Novembre2002` → 45.49951134064658; `method = SunYangXia2013` → 50.957809330628805; `SunYangXia2013` + `background_correction = true` → 47.56338039024118 (codonbias 0.5.0 values, equal to ≤ 1e-14).

### Example 5: CodonW cap

Gene codon_i × ((7i mod 5) + 1) (all 64 codons, ACGT order) under `genetic_code = 2`: default `enc = 60` (60 sense codons); `codonw_cap = true` → 61 (CodonW 1.4.4 `-enc -code 1` prints 61.00).

### Example 6: ENCprime

Gene of Example 1: `method = EncPrimeNc` → 50.75 (ENCprime prints `50.7500`: Phe and Val, n = 5, are not estimated, so the empty 2-fold class is imputed with 1/2 and the 4-fold mean comes from Gly alone; Ile's F = 4/15 is floored at 1/3: Nc = 2 + 18 + 3 + 18.75 + 9), `method = EncPrimeNcPrime` → 48.95085190175721 (ENCprime `48.9509`).

## See Also

- [rscu](rscu.md), [codon_usage_statistics](codon_usage_statistics.md)
