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

Deterministic greedy, not DNA Chisel's randomised constraint solver: rounds of a left-to-right pass (at most 8, each round only runs while the total breach strictly decreased) in which the left-most breaching window is repaired by the synonymous replacement — frequency ≥ `rareCodonThreshold` — of an overlapping codon that minimises, lexicographically, (breach of that window, total breach of the sequence, −codon frequency, codon position, NCBI codon order). A window that no synonymous replacement improves is left alone and the scan continues after it, so an unreachable band cannot loop or throw. The pass returns the lowest-total-breach sequence it reached (review 2026-09 F62): a round that fixes the left-most window first can raise the total breach, and DNA Chisel's solver keeps a mutated sequence only when the constraint score improves (`ConstraintsSolverMixin.resolve_constraints_by_random_mutations`: `if new_score > score: score = new_score else: self.sequence = previous_sequence`), so the repaired output never breaches more than the codon selection it starts from (pre-fix repro `CTGGGTTGGCATTCCAATTGGGATATTGTATCAGCGTCGGTTGGCCGG`, E. coli K-12, window 12, [0.45, 0.55]: selection −3.48333, old output −3.63333; now −3.483333333333332 = dnachisel 3.2.16 `evaluate` of our output). Only the *evaluation* is claimed to equal DNA Chisel; the repaired sequences may differ from that library's (equally valid) random-search output.

## F48 — genetic-code-aware `CodonOptimizer` API (review 2026-09, audit A2-1)

### Sources opened

- dnachisel 3.2.16: `builtin_specifications/EnforceTranslation.py` (`genetic_table="default"` →
  `"Standard"`; `self.backtranslation_table = get_backtranslation_table(genetic_table)`),
  `biotools/biotables.py::get_backtranslation_table` (Biopython
  `CodonTable.unambiguous_dna_by_name[name].forward_table` grouped by amino acid, with
  `table.stop_codons` as the `*` family and `X`/`B`/`J`/`Z` unions), `biotools/sequences_operations.py::translate`
  (`table` is an NCBI table name/number resolved through Biopython), and
  `builtin_specifications/codon_optimization/BaseCodonOptimizationClass.py::get_codons_synonyms`
  (families are the entries of the codon-usage table, which `CodonOptimize` takes per species or
  verbatim — so a non-standard table is expressed through `EnforceTranslation(genetic_table=…)`
  plus a usage table built for that code).
- Biopython 1.88 `Bio.Data.CodonTable.unambiguous_dna_by_id` for tables 2, 4, 11, 27, 28, 31.

### Reference family sets (dnachisel `get_backtranslation_table`, verified by running it)

| Table | R | M | W | * | Q |
|-------|---|---|---|---|---|
| 2 Vertebrate Mitochondrial | CGT CGC CGA CGG | ATA ATG | TGA TGG | TAA TAG AGA AGG | CAA CAG |
| 11 Bacterial | CGT CGC CGA CGG AGA AGG | ATG | TGG | TAA TAG TGA | CAA CAG |
| 27 Karyorelict Nuclear (dual-coding) | CGT CGC CGA CGG AGA AGG | ATG | TGA TGG | TGA | TAA TAG CAA CAG |

Biopython table 11 has the same `forward_table` and `stop_codons` as table 1 (only `start_codons`
differ), so every optimizer output under table 11 equals the Standard one — asserted as a test.

### Fix

One `CodeFamilies` instance per NCBI table, built from the canonical `GeneticCode.CodonTable`
(codon → amino acid, and the amino-acid grouping in NCBI codon order) — the same grouping the
canonical `CodonUsageAnalyzer` core uses for RSCU/CAI/ENC, and the same one DNA Chisel builds from
Biopython. Every code-dependent public method gained an additive optional
`GeneticCode? geneticCode = null` (`null` = Standard): `OptimizeSequence`, `CalculateCAI`,
`RemoveRestrictionSites`, `ReduceSecondaryStructure`, `FindRareCodons` ×2,
`CalculateMinMaxProfile` ×2, `CreateCodonUsageTable`, `CreateCodonTableFromSequence`. No private
codon table remains in the class (the former `StandardCodonToAminoAcid` / `SynonymousCodons` static
copies are now views of the table-1 `CodeFamilies`).

### Reference numbers / locked outputs

| Input | Code | Output | Note |
|-------|------|--------|------|
| `AUGCGUAGACGCAGGUGAAUAUAA`, `MaximizeCAI`, E. coli K12 | 1 | `AUGCGCCGCCGCCGCUGAAUUUAA` | protein `MRRRR*I*` |
| same | 2 | `AUGCGCAGACGCAGGUGGAUGUAA` | protein `MR*R*WM*`: AGA/AGG kept (stops), UGA → UGG (Trp), AUA → AUG (Met) |
| `AUGCGUCGUCGUCGUUAA`, AGA-dominant table | 1 | `AUGAGAAGAAGAAGAUAA` | AGA is the best Arg codon of that table |
| same | 2 | `AUGCGCCGCCGCCGCUAA` | AGA/AGG are stops → best Arg codon is CGC |
| `CGU`×10, `HarmonizeExpression`, E. coli K12 | 2 | `CGUCGUCGUCGUCGCCGCCGCCGCCGACGG` | dnachisel `match_codon_usage` score −0.6 = its own optimizer's score (input −11.8) |
| `AUGCAACAGUGA`, `MaximizeCAI` | 27 | `AUGUAAUAAUGG` | protein `MQQW`: UAA/UAG are Gln codons of table 27 (Biopython / DNA Chisel family) |

DNA Chisel cross-check command (table 2): `DnaOptimizationProblem(sequence=…,
constraints=[EnforceTranslation(genetic_table="Vertebrate Mitochondrial")],
objectives=[CodonOptimize(codon_usage_table={'R': {'CGT':0.38,'CGC':0.40,'CGA':0.06,'CGG':0.10}},
method="match_codon_usage")])` — `objective_scores_sum()` −11.8 for `CGT`×10, −0.6 for both our
allocation and the sequence its own `optimize()` returns (`CGTCGTCGGCGCCGCCGCCGTCGCCGACGT`).

### Decisions

- **Dual-coding tables 27/28/31** are accepted, not rejected: following Biopython's `forward_table`
  and DNA Chisel's back-translation table, a context-dependent stop belongs to the family of the
  amino acid it encodes (table 27: `UAA`/`UAG` ∈ Gln, `UGA` ∈ Trp). Consequence, documented in the
  algorithm doc and the XML docs: no codon unambiguously terminates under those tables, so a
  terminator must be excluded by the caller (the same reason `Translator.FindOrfs` rejects them).
  Also note that the built-in presets' fractions are Standard-code fractions and are not meaningful
  under another code — `CreateCodonTableFromSequence(reference, name, code)` builds a table for it.
- **Code-independent methods** keep their signatures: `FindRareCodonClusters`,
  `FindConservedRareCodonClusters`, `CalculateCodonUsage`, `CompareCodonUsage` perform per-codon
  usage lookups / codon counting only.

## F54 — restriction-site avoidance chained in `OptimizeSequence` (review 2026-09, audit A2-13)

### Sources opened

- dnachisel 3.2.16 (installed): `builtin_specifications/AvoidPattern.py`, `biotools/enzymes_operations.py`
  (`EnzymeSitePattern`), `builtin_specifications/codon_optimization/*`, `DnaOptimizationProblem`
  (`resolve_constraints()` runs the **hard constraints**, `optimize()` then maximises the objectives
  without breaking them) — i.e. site avoidance is a constraint composed with the `CodonOptimize`
  objective in ONE problem, resolved last.
- Our existing `RemoveRestrictionSites` (F22: both strands, IUPAC matching, table-aware, one codon
  per removed occurrence) and `OptimizeSequence` (F42/F43 pass order).

### Defect

`OptimizeSequence` had no way to express the site constraint: the caller had to run
`RemoveRestrictionSites` afterwards by hand, and then the reported CAI/GC/`Changes` no longer
described the sequence actually used. Worse, codon selection itself **creates** sites: with the
E. coli K12 table, Leu `CUG` + Gln `CAG` = `CTGCAG` (PstI), His `CAU` + Met `AUG` = `CATATG` (NdeI),
Thr `ACC` + Met `AUG` = `CCATGG` (NcoI).

### Fix

Additive optional `IReadOnlyList<string>? avoidSites = null` (MCP: `avoid_sites`), applied after
selection → GC pass → structure pass, i.e. the hard constraint last, so no later pass can
re-introduce a site; the metrics and `Changes` are then rebuilt from the final sequence (as F42).

### Verification genes and reference numbers (dnachisel 3.2.16)

Problem per gene: `DnaOptimizationProblem(sequence, constraints=[AvoidPattern("CTGCAG"),
AvoidPattern("CATATG"), AvoidPattern("CCATGG"), EnforceTranslation()],
objectives=[CodonOptimize(species='e_coli', method='use_best_codon')])`.

| Gene | Input (len) | Our `MaximizeCAI` + `avoidSites` output (RNA) | our CAI | dnachisel `AvoidPattern.evaluate` | `CodonOptimize` score |
|------|-------------|-----------------------------------------------|---------|-----------------------------------|-----------------------|
| A `MKLQALHMTMAEF*` | `ATGAAGCTACAAGCATTACACATGACAATGGCTGAGTTTTAA` (42) | `AUGAAACUGCAAGCGCUGCACAUGACGAUGGCGGAAUUUUAA` | 0.870294 (plain 1.0) | 0 / 0 / 0, all pass; `EnforceTranslation` passes | −1.3892 |
| B GFP N-term + `LQHMTMA` | 150 nt | `AUGAGCAAAGGCGAAGAAUUUCUGUUUACCGGCGUGGUGCCGAUUCUGGUGGAACUGGAUGGCGAUGUGAACGGCCAUAAAUUUAGCGUGAGCGGCGAAGGCGAAGGCGAUGCGACCUAUGGCAAACUGCAACACAUGACGAUGGCGUAA` | 0.970251 (plain 1.0) | 0 / 0 / 0, all pass | −1.3892 |
| C `MLQHMTMAK*` | `ATGCTACAACACATGACAATGGCAAAGTAA` (30) | `AUGCUGCAACACAUGACGAUGGCGAAAUAA` | 0.793311 (plain 1.0) | 0 / 0 / 0, all pass | −1.3892 |

Without `avoidSites` the same three outputs score **−1 (fail) on each of the three
`AvoidPattern`s** (sites present), while `EnforceTranslation` passes — the defect reproduced
through the reference library.

For all three genes dnachisel's own `resolve_constraints() + optimize()` returns the
**byte-identical** sequence to ours (`ATGAAACTGCAAGCGCTGCACATGACGATGGCGGAATTTTAA`,
`…AAACTGCAACACATGACGATGGCGTAA`, `ATGCTGCAACACATGACGATGGCGAAATAA`) with the same objective score
−1.3892; exact equality with its stochastic solver is not required, but here the deterministic
greedy repair reaches the same optimum. Protein identical in all cases (`MKLQALHMTMAEF*`,
`MSKGEEFLFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLQHMTMA*`, `MLQHMTMAK*`).

Property (F22) re-checked through the new path: every occurrence removable by one synonymous change
is gone on **both strands**, for `MaximizeCAI`, `BalancedOptimization`, `MinimizeSecondary` and
`HarmonizeExpression`, and the pass is idempotent (re-running `OptimizeSequence(..., avoidSites)` on
its own output changes nothing). GC note: on these genes the replacements are GC-neutral
(0.452381 / 0.526667 / 0.433333 before and after). On other genes the one codon changed per
occurrence can move GC across a bound the GC pass had satisfied; since F60 the pass carries the GC
constraint and prefers GC-feasible candidates (see below).

## F55 — `AvoidRareCodons` spelling (review 2026-09, audit A2-14)

`OptimizationStrategy.AvoidRareCodeons` was misspelled. A correctly spelled member
`AvoidRareCodons = AvoidRareCodeons` (same underlying value 4) was added; the old name stays, kept
**un-attributed** by **policy: additive alias only — the misspelt member stays supported, no
deprecation (campaign backward-compatibility rule)**. `ToString()` still returns `AvoidRareCodeons`
(the field declared first), so no serialized value changes meaning. MCP parsing accepts both
spellings case-insensitively — `System.Text.Json`'s `JsonStringEnumConverter` (the converter the
MCP SDK uses for enum parameters) resolves `"AvoidRareCodons"`, `"avoidrarecodons"`,
`"AVOIDRARECODONS"`, `"AvoidRareCodeons"` and `"avoidrarecodeons"` all to that value (asserted in
the tests). `Enum.GetValues<OptimizationStrategy>()` now returns 6 entries with 5 distinct values;
the two consumers that iterate it (`Properties/CodonProperties`,
`Metamorphic/CodonMetamorphicTests`) run the rare-codon strategy twice and are unaffected — no
counting or indexing of strategies anywhere in the tests.

## F60 — GC-constraint-aware `avoidSites` pass (review 2026-09, audit A3-6)

**Defect.** `OptimizeSequence(..., avoidSites)` runs the site-removal pass *after* the GC pass, and
the pass ranked site-removing synonymous candidates by codon-usage frequency alone
(`if (frequency > bestFrequency)`), so the best-frequency candidate could move the GC fraction of a
window (or the global GC fraction) back outside the `gcTargetMin`/`gcTargetMax` bounds that the GC
pass had just satisfied. DNA Chisel 3.2.16 composes `AvoidPattern` **and** `EnforceGCContent` as
hard constraints of ONE `DnaOptimizationProblem` and resolves both, so breaking the GC constraint to
clear a site is not reference behaviour.

**Reference runs (dnachisel 3.2.16).** Four constructed E. coli genes. For each, "balanced" is our
`BalancedOptimization` output *without* `avoidSites` (GC constraint satisfied, one site present),
"old" the historical GC-blind chain (`RemoveRestrictionSites` on that sequence), "new" the
GC-constraint-aware output. Problem:
`DnaOptimizationProblem(sequence=<seq>, constraints=[AvoidPattern("CTGCAG"), AvoidPattern("CATATG"),
AvoidPattern("CCATGG"), EnforceGCContent(mini, maxi, window=w), EnforceTranslation()])`.

| Gene (protein) | bounds / window | old chain | new (ours) |
|---|---|---|---|
| `ATGATTAGACTTTGCACGCCAGCCTATGGACGAAGACGCTAA` (`MIRLCTPAYGRRR*`) | 0.40–0.60, w=24 | `AvoidPattern` ×3 pass, `EnforceGCContent` **−0.225 FAIL** ("Out of bound on segments 5-40") | all constraints **pass** |
| `ATGCGTATGGCCCTCCGTTCGGCTGTGGGCTGGCTAGTGTAA` (`MRMALRSAVGWLV*`) | 0.40–0.60, w=24 | `EnforceGCContent` **−0.225 FAIL** (segments 1-39) | all **pass** |
| `ATGGGTTCCGAAACGTCTGCTTACGGTGCAGCGTTGTGGTAA` (`MGSETSAYGAALW*`) | 0.40–0.60, w=24 | `EnforceGCContent` **−0.150 FAIL** (segments 2-40) | all **pass** |
| `ATGGTAACTGCGGTACCCATGAGGGGGATGGGCAGACGCTAA` (`MVTAVPMRGMGRR*`) | 0.40–0.60, global | `EnforceGCContent` **−0.019048 FAIL** (GC 0.619048) | all **pass** |

For all four genes our GC-aware output is **byte-identical** to the sequence dnachisel's own
`resolve_constraints() + optimize()` (objective `CodonOptimize(species='e_coli',
method='use_best_codon')`) returns from the same GC-balanced start:
`ATGATTCGCCTGTGCACTCCAGCTTATGGTCGTCGTCGCTAA` (objective −3.14),
`ATGCGTATGGCACTGCGTTCTGCTGTTGGTTGGCTGGTGTAA` (−2.62),
`ATGGGCAGCGAAACCTCTGCTTATGGTGCAGCGCTGTGGTAA` (−2.16),
`ATGGTTACTGCTGTTCCAATGCGCGGCATGGGCCGCCGCTAA`. Proteins unchanged in every case. Our
`EvaluateGcContent` breach numbers equal dnachisel's `EnforceGCContent` scores exactly
(0.225 / 0.225 / 0.150 / 0.019048, and 0.483333 / 1.066667 below).

**Fallback (decided, documented).** With unreachable bounds — gene 3 under
`EnforceGCContent(0.30, 0.45, window=24)` — dnachisel's `resolve_constraints()` raises
`NoSolutionError`. Our pass never throws: `AvoidPattern` stays the hard constraint, the site is
removed, and the least-breach candidate is chosen, i.e. the pre-existing breach is not increased
(dnachisel scores: GC-balanced start −0.483333, GC-blind chain −1.066667, ours −0.483333).

**Fix.** An internal constraint descriptor (`GcConstraint(mini, maxi, window)`, the
`EvaluateGcContent` semantics) is threaded `OptimizeSequence` → internal `RemoveRestrictionSites`
overload → `TryRemoveOccurrence`, for the strategies that declare a GC target
(`BalancedOptimization`, `MinimizeSecondary`) only. Candidate ranking is lexicographic:
(a) GC-constraint-satisfying candidates — or, when the input already breaches, candidates that do
not increase the breach score — then (b) codon frequency, then (c) the left-most codon, then
(d) NCBI codon order. The public `RemoveRestrictionSites` overloads are unchanged (GC-blind), and
`MaximizeCAI`/`AvoidRareCodons`/`HarmonizeExpression` thread no constraint, so the F54 locked
outputs are byte-identical. The structure pass (`MinimizeSecondary`) runs *before* the site pass and
is not re-run: the site pass is GC-aware only, documented as such.

## Date
2026-03-10 (reviewed 2026-09-28; F43 windowed GC, F54 avoidSites, F55 AvoidRareCodons, F60 GC-aware avoidSites 2026-10-14)
