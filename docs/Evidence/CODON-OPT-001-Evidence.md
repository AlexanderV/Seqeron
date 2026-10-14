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
- GC content balancing in BalancedOptimization strategy (40-60% target); global by default, per sliding window with `gcWindow` / MCP `gc_window` (review 2026-09 F43)
- BalancedOptimization rebuilds Changes list after GC balancing to reflect all modifications
- CAI via the canonical CodonUsageAnalyzer core: Met/Trp and stops not scored, w < 0.0001 → 0.01 (CodonW `cai_out`) — review 2026-09, CODON-CAI-001 F12/F13
- MinimizeSecondary strategy uses the BalancedOptimization codon selection and GC pass and then runs the `ReduceSecondaryStructure` structure pass (review 2026-09 F42; it previously stopped after the GC pass)

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

### Numeric cross-checks (2026-10-14, F41/F42 — Turner 2004 window scoring)

Reference implementation: ViennaRNA 2.7.2 (Python `RNA`), `md = RNA.md(); md.dangles = 2` (the RNAfold default), `RNA.fold_compound(window, md).mfe()`. Energy model: NNDB Turner 2004, evaluated here by the canonical engine `RnaSecondaryStructure.CalculateMinimumFreeEnergy` (owned by B12, locked against NNDB).

| 45-nt window | ViennaRNA 2.7.2 (dangles=2) | Seqeron | Δ |
|---|---|---|---|
| gene A [0..45) `AUGGGCGGC…GGC` | −25.80 | −25.82 | 0.02 |
| gene A [3..48) `GGCGGCGGC…GGC` | −29.00 | −29.34 | 0.34 |
| gene A [6..51) `GGCGGC…GGCAAG` | −26.50 | −25.82 | 0.68 |
| gene B [0..45) `AUGAAACUGGCA…CUG` | −9.20 | −9.22 | 0.02 |
| gene B [3..48) `AAACUGGCA…CGU` | −14.00 | −13.44 | 0.56 |
| gene B [6..51) `CUGGCAGAU…CAG` | −14.00 | −13.82 | 0.18 |

Over 40 random windows (20–60 nt, seed 42): mean |Δ| 0.50 kcal/mol, max 2.24, identical MFE dot-bracket in 30/40. The residual is the engine's documented open deviation (RNA-STRUCT-001 Evidence D3/D4: NNDB int21/int22 lookup tables not inlined), not a model difference; the tests lock our values exactly and assert agreement with ViennaRNA within 1.0 kcal/mol on the six windows above.

| Case | Reference | Reference value | Seqeron |
|---|---|---|---|
| DNA Chisel free-energy objective | dnachisel 3.2.16 (installed) | none — `AvoidHairpins` is the IDT reverse-complement stem rule; no module references ViennaRNA, `RNA.fold` or a free energy | energy model sourced from ViennaRNA RNAfold / NNDB instead |
| Single fold, 43-nt window | `RnaSecondaryStructure.CalculateMinimumFreeEnergy` (Release) | 4.42 ms | — |
| `ReduceSecondaryStructure` on a 999-nt CDS, `PairCount` | this implementation (Release) | 0.06 s | — |
| `ReduceSecondaryStructure` on a 999-nt CDS, `Turner2004Mfe` | this implementation (Release, memoised windows) | 9.02 s, 193 codons changed, summed 45-nt window MFE −2559.46 → −984.27 kcal/mol, worst window −22.79 → −11.56 | — |
| Turner pass, gene A (111 nt) | this implementation | summed window MFE −477.92 → −144.26, worst −37.34 → −16.04 kcal/mol, protein unchanged | — |
| Turner pass, gene B (99 nt) | this implementation | summed window MFE −252.88 → −96.45, worst −17.00 → −7.18 kcal/mol, protein unchanged | — |

Default: `PairCount` stays the default of `ReduceSecondaryStructure` and of the MinimizeSecondary strategy — it is ~150× cheaper (0.06 s vs 9.02 s per kilobase) and is the historical behaviour; `Turner2004Mfe` is the sourced energy model and is selected per call (`scoring` / `secondaryScoring`, MCP `scoring` / `secondary_scoring`).

### Behaviour changes locked by tests

- `MaximizeCAI` ties broken by NCBI codon order (deterministic; Biopython only warns on ties).
- `HarmonizeExpression` deterministic (largest-remainder allocation) — was `new Random()` weighted sampling.
- `AvoidRareCodeons` / `BalancedOptimization` fall back to the best synonymous codon when no synonym reaches the threshold.
- GC pass: no neutral swaps, no overshoot-without-alternative, stops on entering the window, frequency floor = `rareCodonThreshold`.
- `RemoveRestrictionSites`: one codon per occurrence, highest-frequency substitution, IUPAC sites, both strands.
- `ReduceSecondaryStructure`: canonical `RnaSecondaryStructure.CanPair` (G·U wobble), re-evaluated baseline, normalised RNA output.
- `ReduceSecondaryStructure`: opt-in `Turner2004Mfe` scoring (Turner 2004 window MFE, maximised), memoised window folds, `maxFoldEvaluations` budget; default `PairCount` behaviour unchanged (F41).
- `MinimizeSecondary`: runs the structure pass after the GC pass; output differs from BalancedOptimization on structured sequences (F42).
- Windowed GC constraint (F43): `EvaluateGcContent` equals dnachisel's `EnforceGCContent(...).evaluate(problem)` (score, breaching window starts, grouped breach segments) and `OptimizeSequence(..., gcWindow: w)` repairs every window where that is achievable; default (`gcWindow = null`) global behaviour unchanged.
- `CreateCodonTableFromSequence`: Sharp & Li / Biopython 0.5 pseudo-count, all 64 codons, Standard-code `CodonToAminoAcid`.

### Review 2026-09 — F43: windowed GC constraint (audit A1-7)

#### Sources actually opened

| Source | What was opened | What it establishes |
|---|---|---|
| DNA Chisel 3.2.16 (installed, executed) | `dnachisel/builtin_specifications/EnforceGCContent.py` | `evaluate` = "sum of breaches extent for all windowed breaches": per window `max(0, mini-gc) + max(0, gc-maxi)`, `score = -breaches.sum()`, breach starts = `(breaches > 0).nonzero()`, one breach → `[start, start+window]`, several → `group_nearby_segments(segments, max_start_spread=max(1, locations_span=50))` → `(group[0][0], group[-1][-1])`; windowless mode reports the whole location |
| DNA Chisel 3.2.16 | `dnachisel/biotools/gc_content.py` | Sliding GC via cumulative sums: window `i` = `sequence[i : i+window_size]`, `len(sequence) - window_size + 1` values, fraction = GC/window (so both bounds are inclusive; a window longer than the sequence yields no values) |
| DNA Chisel 3.2.16 | `dnachisel/biotools/indices_operations.py` (`group_nearby_segments`) | Grouping rule used for the breach locations (`seg[0] - current_group[0][0] < max_start_spread`) |
| DNA Chisel 3.2.16 (executed) | `DnaOptimizationProblem(constraints=[EnforceTranslation, EnforceGCContent(window=20)]).resolve_constraints()` | Feasibility reference: solvable for the mixed gene only by also mutating the stop codon (TAA→TAG); `NoSolutionError` for the AT-rich, GC-rich and two-segment genes |

#### Numeric cross-checks (2026-10-14) — bounds 0.40/0.60

| Sequence (DNA) | window | dnachisel score | Seqeron `EvaluateGcContent` | breach segments (both) |
|---|---|---|---|---|
| `ATGAAAGAAAAAATTCTTAAAGGCGGCGGCGCGGCGGCCATGAAAAAATTTTAA` (mixed, 54 nt) | 20 | −4.25 | −4.25 | 0-54 (starts 0–7, 13–27, 33, 34) |
| same | none | −0.0 | 0 | — ("Passed !") |
| same | 9 (bounds 0.30/0.70) | −7.655555555556 | identical | 0-54 |
| `ATGGCCGCCGCCGCGGCGGCCCTGCTGCTGCGCCGCGGCGGCGCCGCCTGA` (GC-rich, 51 nt) | 20 | −9.45 | −9.45 | 0-51 |
| same | none | −0.262745098039 | identical | 0-51 |
| `ATGAAATTTAAA…TTTTAA` (AT-rich, 54 nt) | 20 | −13.85 | −13.85 | 0-54 |
| same | none | −0.381481481481 | identical | 0-54 |
| AT-rich flanks + balanced middle (two-segment, 177 nt) | 20 | −20.1 | −20.1 | 0-51, 125-177 (two groups) |
| `ATG(AAAGAA)×7TAA` (M-(KE)₇-*, 48 nt) | 20 / 30 | −6.7 / −4.366667 | identical | 0-48 |
| `ATG(GGCAGC)×7TAA` (M-(GS)₇-*, 48 nt) | 20 / 30 | −6.4 / −4.166667 | identical | 0-48 |

#### Repair results (dnachisel `EnforceGCContent(0.40, 0.60, window=w).evaluate` applied to **our** output)

| Input | mode | output | dnachisel score of output |
|---|---|---|---|
| M-(KE)₇-* | `gcWindow: 20` | `AUGAAAGAGAAGGAGAAGGAGAAAGAGAAGGAGAAGGAGAAAGAGUAA` (11 codons changed) | **0** (input −6.7) |
| M-(KE)₇-* | `gcWindow: 30` | `AUGAAAGAAAAGGAGAAGGAGAAGGAGAAGGAGAAAGAGAAAGAGUAA` | **0** at w=30 (input −4.366667) |
| M-(GS)₇-* | `gcWindow: 20` | `AUGGGCUCUGGUUCUGGUUCUGGCUCUGGUUCUGGUUCUGGCAGCUAA` (10 codons changed) | **0** (input −6.4) |
| M-(GS)₇-* | default (global) | `AUGGGUUCUGGUUCUGGUUCUGGUUCUGGCAGCGGCAGCGGCAGCUAA` | −1.5 at w=20 (global GC pass satisfies the global bounds only) |
| AT-rich (Lys/Phe only) | `gcWindow: 20` | `AUGAAGUUCAAG…UUCUAA` | −2.4 (input −13.85); unreachable — dnachisel's own solver: `NoSolutionError` |
| GC-rich | `gcWindow: 20` | `AUGGCAGCAGCA…GCAUGA` | −2.0 (input −9.45); unreachable (`NoSolutionError`) |
| two-segment (177 nt) | `gcWindow: 20` | see tests | −2.6 (input −20.1; global-only default: −32.45); unreachable (`NoSolutionError`) |
| mixed | `gcWindow: 20` | `AUGAAGGAGAAGAUCCUGAAAGGUGGUGGUGCAGCAGCGAUGAAGAAGUUCUAA` | −0.2 (input −4.25); dnachisel reaches 0 only by mutating the stop codon TAA→TAG and using codons below our `rareCodonThreshold` floor, both of which `OptimizeSequence` never does |

Protein preserved in every row (`MKEKEKEKEKEKEKE*`, `MGSGSGSGSGSGSGS*`, `MKFKFKFKFKFKFKFKF*`, `MAAAAAALLLRRGGAA*`, `MKEKILKGGGAAAMKKF*`), sequence length unchanged, outputs deterministic (repeated calls identical).

#### Search

Deterministic greedy, not DNA Chisel's randomised constraint solver: rounds of a left-to-right pass (at most 8, each round only runs while the total breach strictly decreased) in which the left-most breaching window is repaired by the synonymous replacement — frequency ≥ `rareCodonThreshold` — of an overlapping codon that minimises, lexicographically, (breach of that window, total breach of the sequence, −codon frequency, codon position, NCBI codon order). A window that no synonymous replacement improves is left alone and the scan continues after it, so an unreachable band cannot loop or throw. Only the *evaluation* is claimed to equal DNA Chisel; the repaired sequences may differ from that library's (equally valid) random-search output.

## Date
2026-03-10 (reviewed 2026-09-28; F43 windowed GC 2026-10-14)
