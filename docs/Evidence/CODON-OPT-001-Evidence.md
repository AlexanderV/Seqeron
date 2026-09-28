# CODON-OPT-001: Codon Optimization Evidence

## Test Unit
- **ID:** CODON-OPT-001
- **Area:** Codon
- **Algorithm:** Sequence Optimization (OptimizeSequence)

## Sources

### Primary Sources

1. **Wikipedia - Codon Usage Bias**
   - URL: https://en.wikipedia.org/wiki/Codon_usage_bias
   - Key concepts:
     - Codon optimization adjusts codons to match host tRNA abundances for heterologous expression
     - Strategies include: local mRNA folding, codon pair bias, codon ramp, codon harmonization
     - Rare codons can lead to inefficient use and depletion of ribosomes
     - Protein folding is affected by translation rates (cotranslational folding)

2. **Wikipedia - Codon Adaptation Index (CAI)**
   - URL: https://en.wikipedia.org/wiki/Codon_Adaptation_Index
   - Reference: Sharp, P.M. & Li, W.H. (1987). "The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications." Nucleic Acids Research. 15(3): 1281–1295.
   - Key concepts:
     - CAI measures deviation from reference gene set
     - Relative adaptiveness (wi) = fi / max(fj) for synonymous codons
     - CAI = geometric mean of weights: CAI = (∏wi)^(1/L)
     - Range: 0 < CAI ≤ 1

3. **Kazusa Codon Usage Database**
   - URL: https://www.kazusa.or.jp/codon/
   - Provides organism-specific codon usage tables

### Key Academic References

4. **Sharp & Li (1987)**
   - "The codon adaptation index-a measure of directional synonymous codon usage bias"
   - Original CAI definition
   - PMC 340524, PMID 3547335

5. **Plotkin & Kudla (2011)**
   - "Synonymous but not the same: The causes and consequences of codon bias"
   - Nature Reviews Genetics 12(1): 32-42
   - Comprehensive review of codon optimization strategies

6. **Mignon et al. (2018)**
   - "Codon harmonization - going beyond the speed limit for protein expression"
   - FEBS Letters 592(9): 1554-1564
   - Describes HarmonizeExpression strategy

## Documented Behaviors

### Optimization Strategies (from Wikipedia & implementation)

1. **MaximizeCAI**: Use most frequent codons for each amino acid
   - Source: Sharp & Li (1987), CAI definition

2. **BalancedOptimization**: Balance CAI with GC content constraints
   - Source: Wikipedia - mRNA secondary structure affects translation

3. **HarmonizeExpression**: Match host codon usage distribution
   - Source: Mignon et al. (2018)

4. **AvoidRareCodons**: Replace only rare codons (frequency < threshold)
   - Source: Wikipedia - consecutive rare codons inhibit translation

5. **MinimizeSecondary**: Avoid mRNA secondary structures
   - Source: Wikipedia - 5' secondary structure inhibits translation

### Invariants (from theory)

1. **Protein preservation**: Optimization must not change encoded protein
   - Source: Definition of synonymous codon substitution

2. **CAI range**: 0 < CAI ≤ 1
   - Source: Sharp & Li (1987) - geometric mean of values 0 < wi ≤ 1

3. **Methionine/Tryptophan unchanged**: AUG and UGG are unique codons
   - Source: Standard genetic code

## Test Datasets

### Organism-Specific Codon Tables

1. **E. coli K12** preferred codons (Kazusa species=316407, W3110 K-12 substrain, 4332 CDS):
   - Leu: CUG (0.50)
   - Arg: CGC/CGU (0.40/0.38), AGA/AGG rare (0.04/0.02)
   - Pro: CCG (0.53)
   
2. **S. cerevisiae (Yeast)** preferred codons (Kazusa species=4932, 14411 CDS):
   - Leu: UUA/UUG (0.28/0.29)
   - Arg: AGA (0.48)
   - Pro: CCA (0.42)

3. **H. sapiens (Human)** preferred codons (Kazusa species=9606, 93487 CDS):
   - Leu: CUG (0.40)
   - Arg: AGA/AGG (0.21 each)
   - Pro: CCC (0.32)

### Edge Cases (from theory)

1. **Empty sequence**: Should return empty result
2. **Incomplete codons**: Trim to complete codons (length % 3 == 0)
3. **DNA input (T)**: Convert to RNA (U)
4. **Lowercase input**: Case-insensitive processing
5. **Stop codons**: Preserved, not optimized

## Testing Methodology

1. **Unit tests**: Verify individual method behaviors
2. **Invariant tests**: Verify protein preservation across all strategies
3. **CAI mathematical tests**: Verify CAI calculation formula
4. **Organism-specific tests**: Verify different organisms yield different optimizations
5. **Edge case coverage**: Empty, single codon, all-same codons

## Known Failure Modes

1. **Invalid codon**: Unknown codon should translate to 'X' or error
2. **Non-RNA characters**: Should handle gracefully
3. **Stop codon in middle**: May terminate protein prematurely

## Implementation Notes

- Implementation uses RNA notation (U not T)
- Automatically converts T → U
- Trims to complete codons
- GC content balancing in BalancedOptimization strategy (40-60% target)
- BalancedOptimization rebuilds Changes list after GC balancing to reflect all modifications
- CAI via the canonical CodonUsageAnalyzer core: Met/Trp and stops not scored, w < 0.0001 → 0.01 (CodonW `cai_out`) — review 2026-09, CODON-CAI-001 F12/F13
- MinimizeSecondary strategy delegates to BalancedOptimization in codon selection; dedicated `ReduceSecondaryStructure` method handles secondary structure reduction separately

## Review 2026-09 (campaign review-2026-09, batch B02)

### Sources actually opened

| Source | What was opened | What it establishes |
|---|---|---|
| DNA Chisel 3.2.16 (PyPI, installed and executed) | `dnachisel.builtin_specifications.codon_optimization.MaximizeCAI` / `MatchTargetCodonUsage` source, `CodonOptimize` docs page | `use_best_codon` = most frequent synonymous codon; `match_codon_usage` objective = −Σ_aa n_aa·Σ_codon |f_seq − f_table|; `harmonize_rca` needs a *source* host table |
| Biopython 1.88 (installed) | `Bio.SeqUtils.CodonAdaptationIndex.__init__` / `.calculate` / `.optimize` source | "we use a value of 0.5 for codons that do not appear in the reference sequences"; w = count/max(count in family); Met/Trp/stops excluded from CAI |
| python_codon_tables 0.1.18 | `e_coli_316407` table | Preset EColiK12 values (Leu CTG 0.50, Arg CGC 0.40 …) confirmed identical |
| Edinburgh Genome Foundry DnaChisel documentation (web) | `CodonOptimize` method descriptions | Wording of the three published methods |

### Numeric cross-checks (2026-09-28)

| Case | Reference | Reference value | Seqeron |
|---|---|---|---|
| `ATGAGCAAAGGTGAAGAACTGTTCACCGGTGTTGTTCCGATTCTGGTTGAACTGGATGGTGATGTTAAC`, E. coli 316407, `use_best_codon` | DNA Chisel 3.2.16 | `ATGAGCAAAGGCGAAGAACTGTTTACCGGCGTGGTGCCGATTCTGGTGGAACTGGATGGCGATGTGAAC` | identical (U spelling), CAI 1.0000 |
| same gene, `match_codon_usage` objective score | DNA Chisel optimizer output | −11.77 (original −23.25, best-codon −21.67) | HarmonizeExpression output scores −11.77 |
| reference set `ATGAAAGCGTTCAAGCGTACTGCGATGCCCAAAGGGTTTTAA` → relative adaptiveness | Biopython `CodonAdaptationIndex` | AAA 1.0, AAG 0.5, GCG 1.0, GCT 0.25, TTT 1.0, TTC 1.0, CGT 1.0, AGA 0.5, TAA 1.0, TAG 0.5 | identical (w = f / max f from the built table) |
| 43-nt window MFE cost (why the structure pass stays a heuristic) | `RnaSecondaryStructure.CalculateMinimumFreeEnergy` | ≈6.0 ms/window (Release) | — |

### Behaviour changes locked by tests

- `MaximizeCAI` ties broken by NCBI codon order (deterministic; Biopython only warns on ties).
- `HarmonizeExpression` deterministic (largest-remainder allocation) — was `new Random()` weighted sampling.
- `AvoidRareCodeons` / `BalancedOptimization` fall back to the best synonymous codon when no synonym reaches the threshold.
- GC pass: no neutral swaps, no overshoot-without-alternative, stops on entering the window, frequency floor = `rareCodonThreshold`.
- `RemoveRestrictionSites`: one codon per occurrence, highest-frequency substitution, IUPAC sites, both strands.
- `ReduceSecondaryStructure`: canonical `RnaSecondaryStructure.CanPair` (G·U wobble), re-evaluated baseline, normalised RNA output.
- `CreateCodonTableFromSequence`: Sharp & Li / Biopython 0.5 pseudo-count, all 64 codons, Standard-code `CodonToAminoAcid`.

## Date
2026-03-10 (reviewed 2026-09-28)
