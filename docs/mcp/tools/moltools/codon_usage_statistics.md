# codon_usage_statistics

Aggregate codon-usage statistics for a coding sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `codon_usage_statistics` |
| **Method ID** | `CodonUsageAnalyzer.GetStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Produces a one-shot codon-usage report: per-codon counts, RSCU (Relative Synonymous Codon Usage), ENC (Effective Number of Codons, Wright), total codons, GC% at codon positions 1/2/3, GC3s (GC% at the third position of synonymous codons only, excluding Met/Trp/stop), and overall GC = (GC1+GC2+GC3)/3. `gc_convention` selects how GC1/GC2/GC3 and overall GC treat stop codons: `EmbossCusp` (default; every valid codon, stops included — EMBOSS `cusp`, independent of the genetic code) or `CodonW` (stop codons of `genetic_code` excluded — CodonW 1.4.4 `gc_out`; overall GC then equals CodonW's `GC` index; CodonW reports fractions, the tool percent). `totalCodons`, GC3s, counts, RSCU and ENC do not depend on the convention. `isEncCalculable` is `false` when ENC cannot be calculated (CodonW `*****`; `enc` = 0).

With `include_fop_cbi = true` the tool also reports CodonW 1.4.4's **Fop** (frequency of optimal codons, Ikemura 1981: optimal / total) and **CBI** (codon bias index, Bennetzen & Hall 1982: (optimal − expected)/(total − expected), expected = Σ n_aa·k_aa/d_aa), both over the amino acids that have more than one codon in `genetic_code` (stops and Met/Trp excluded) and at least one optimal codon — exactly CodonW `-fop -cbi` (`codon_us.c fop_out`/`cbi_out`). Optimal codons come from CodonW's built-in set `fop_cbi_organism` (`-f_type`; default *E. coli*, CodonW's default) or the user list `optimal_codons` (`-fop_file`/`-cbi_file`). `fop`/`cbi` are `null` when not requested or not calculable (CodonW prints 0.000, e.g. a gene of only Met/Trp/stop codons). Verified on 7 280 Fop / 6 160 CBI values of the CodonW binary.

## Core Documentation Reference

- Source: [CodonUsageAnalyzer.cs#L662](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs#L662) (`GetStatistics(string, GeneticCode, CodonGcConvention)`)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Coding DNA sequence (frame 0), non-empty. |
| `genetic_code` | integer | No | NCBI translation table number (default `1` = Standard; e.g. `2` vertebrate mitochondrial, `11` bacterial/plastid). Defines the synonymous families and stops (CodonW `-code`). One of 1–6, 9–16, 21–33. |
| `gc_convention` | string | No | `EmbossCusp` (default) or `CodonW`. Convention for GC1/GC2/GC3 and overall GC: stops included (EMBOSS cusp) or stops of `genetic_code` excluded (CodonW `gc_out`). |
| `include_fop_cbi` | boolean | No | `false` (default): `fop`/`cbi` null, output otherwise unchanged. `true`: compute CodonW Fop and CBI. |
| `fop_cbi_organism` | string | No | CodonW built-in optimal-codon set (`-f_type`): `EscherichiaColi` (default), `BacillusSubtilis`, `DictyosteliumDiscoideum`, `AspergillusNidulans`, `SaccharomycesCerevisiae`, `DrosophilaMelanogaster`, `CaenorhabditisElegans`, `NeurosporaCrassa` (case/spaces ignored; `E. coli` style accepted). |
| `optimal_codons` | string[] | No | User optimal codons (DNA/RNA, any case), overriding `fop_cbi_organism` (CodonW `-fop_file`/`-cbi_file`); codons of stops or single-codon amino acids are ignored. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `codonCounts` | object | Codon → count. |
| `rscu` | object | Codon → RSCU. |
| `enc` | number | Effective Number of Codons (20..61; 0 when not calculable). |
| `isEncCalculable` | boolean | `false` when ENC cannot be calculated (CodonW `*****`). |
| `totalCodons` | integer | Number of complete codons. |
| `gc1` / `gc2` / `gc3` | number | GC% at codon position 1 / 2 / 3 (denominator: all codons for `EmbossCusp`, sense codons for `CodonW`). |
| `gc3s` | number | GC% at the third position of synonymous codons. |
| `overallGc` | number | (GC1+GC2+GC3)/3. |
| `fop` | number \| null | CodonW Fop in [0, 1]; null unless `include_fop_cbi` and calculable. |
| `cbi` | number \| null | CodonW CBI (1 = only optimal codons, ≈ 0 random, < 0 optimal codons avoided); null unless `include_fop_cbi` and calculable. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | Unknown GC convention: X. Supported: EmbossCusp, CodonW. |
| 1004 | Unknown Fop/CBI organism: X. Supported: EscherichiaColi, BacillusSubtilis, DictyosteliumDiscoideum, AspergillusNidulans, SaccharomycesCerevisiae, DrosophilaMelanogaster, CaenorhabditisElegans, NeurosporaCrassa. |
| 1005 | Invalid codon 'X': expected three of A, C, G, T/U. |
| 1002 | Unknown genetic code table: N. Supported NCBI tables: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |

## Examples

### Example 1: Three-codon gene `ATGTTTAAA`

Codons ATG, TTT, AAA. GC1=0, GC2=0, GC3=33.33% (only the G in ATG), overall GC=11.11%. GC3s=0 (synonymous codons TTT, AAA have T/A at position 3). RSCU[TTT]=2.0 (Phe family, only TTT used).

### Example 2: Non-standard genetic code (`genetic_code = 2`)

`ATAATG`: under the default table 1 ATA (Ile, T at position 3) is the only synonymous codon, so GC3s = 0; under NCBI table 2 Met = {ATA, ATG} is a two-codon family, so both codons are synonymous and GC3s = 50% (CodonW 1.4.4 `-gc3s`: 0.000 / 0.500). Counts and GC1/GC2/GC3 are unchanged.

### Example 3: GC convention (`gc_convention = "CodonW"`)

`ATGAAACTGTTCGCCAAATAA` (terminal TAA): default (EMBOSS cusp, 7 codons) overall GC = 33.33 %; `CodonW` (6 sense codons) GC1/GC2/GC3 = 33.33/16.67/66.67 %, overall GC = 38.89 % (CodonW 1.4.4 `-gc` 0.389, `-base` GC1-3 0.333/0.167/0.667). `ATGGCTTAAGCCNNNGCGTGAGGCCTGCTTTAGAAAAAGTGGTAA` with `genetic_code = 2` (TGA = Trp): CodonW overall GC 57.58 % (CodonW 0.576).

### Example 4: Fop and CBI (`include_fop_cbi = true`)

`ATGCTGCTGCTGCTGAAAAAGTTTTTCGGCGGTGGATAA`, default *E. coli* set: the scored codons are Leu CTG×4, Lys AAA/AAG, Phe TTT/TTC, Gly GGC/GGT/GGA (11; Met and the stop are excluded); optimal CTG×4, AAA, TTC, GGC, GGT = 8 ⇒ `fop` = 8/11 = 0.727. Expected optimal codons under random usage = 4·1/6 + 2·1/2 + 2·1/2 + 3·2/4 = 25/6 ⇒ `cbi` = (8 − 25/6)/(11 − 25/6) = 23/41 = 0.561 (CodonW 1.4.4 `-fop -cbi`: 0.727 / 0.561). With `fop_cbi_organism = "BacillusSubtilis"` and `genetic_code = 2`: `fop` 0.273, `cbi` −0.055 (CodonW `-f_type 1 -code 1`).

## Performance

- **Time Complexity:** O(n) in sequence length.
- **Space Complexity:** O(k) in distinct codons.

## See Also

- [rscu](rscu.md), [effective_number_of_codons](effective_number_of_codons.md), [count_codons](count_codons.md)
