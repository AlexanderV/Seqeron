# seqeron-rna-structure — tool map (all 17)

Server: **analysis**. One backing class: `RnaSecondaryStructure.*`.
This skill is **not** in `domain-map.json`, so it has **no** generated `_generated/tools.md` —
**this curated map is the index.** Verify schemas in `docs/mcp/tools/analysis/<tool>.md`.

> Coordinates: base-pair indices are **0-based** positions into the sequence; dot-bracket indices
> are `(openingIndex, closingIndex)`. Energies are **kcal/mol at 37 °C**, Turner 2004
> nearest-neighbour parameters. `minimum_free_energy` returns **0** when no structure forms.
> Always confirm exact I/O in the tool doc.

## Fold / predict structure

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `predict_rna_structure` | `RnaSecondaryStructure.PredictStructure` | Greedy non-overlapping stem-loops (most stable first) → sequence, dot-bracket, base pairs, stem-loops, pseudoknots, total MFE. Nested pairs only. | `predict_rna_structure.md` |
| `minimum_free_energy` | `RnaSecondaryStructure.CalculateMinimumFreeEnergy` | Zuker-style O(n³) DP MFE (Turner 2004), kcal/mol; 0 if no structure. `minLoopSize` default 3 (values < 3 clamped). | `minimum_free_energy.md` |

## Stem-loops / complementary arms

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_stem_loops` | `RnaSecondaryStructure.FindStemLoops` | Enumerate hairpin candidates: stem bounds/pairs/energy, loop type/bounds/size/seq, dot-bracket, Turner total energy. G-U allowed when `allowWobble`. | `find_stem_loops.md` |
| `find_rna_inverted_repeats` | `RnaSecondaryStructure.FindInvertedRepeats` | Antiparallel complementary arms (left arm reverse-complements right arm) across a `minSpacing..maxSpacing` loop → potential stems. **RNA** complement rules (distinct from DNA `find_inverted_repeats`). | `find_rna_inverted_repeats.md` |

## Pseudoknots

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `detect_pseudoknots` | `RnaSecondaryStructure.DetectPseudoknots` | Crossing base pairs: for ordered pairs `(i,j)`/`(k,l)`, a pseudoknot iff `i<k<j<l`; nested/disjoint excluded. Each crossing pair-of-pairs reported once. Input pairs carry `type` WatsonCrick/Wobble/NonCanonical. | `detect_pseudoknots.md` |

## Base pairing

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `can_pair` | `RnaSecondaryStructure.CanPair` | `true` for A-U, G-C, or G-U wobble (case-insensitive), else `false`. | `can_pair.md` |
| `base_pair_type` | `RnaSecondaryStructure.GetBasePairType` | Classify a pair → `WatsonCrick` (A-U, G-C) / `Wobble` (G-U) / `null` if it cannot pair. | `base_pair_type.md` |
| `rna_complement_base` | `RnaSecondaryStructure.GetComplement` | RNA complement of one base: A↔U, G↔C (IUPAC ambiguity mapped; T→A). | `rna_complement_base.md` |

## Dot-bracket notation

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `validate_dot_bracket` | `RnaSecondaryStructure.ValidateDotBracket` | Every bracket balanced + family-matched (`()[]{}<>` + letter pairs); `(]` rejected. | `validate_dot_bracket.md` |
| `parse_dot_bracket` | `RnaSecondaryStructure.ParseDotBracket` | Parse into `(openingIndex, closingIndex)` tuples, one stack per bracket family; unpaired `.,-:_~` skipped. | `parse_dot_bracket.md` |

## Free-energy terms (Turner 2004, kcal/mol at 37 °C)

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `stem_energy` | `RnaSecondaryStructure.CalculateStemEnergy` | Nearest-neighbour stacking over consecutive pairs + AU/GU terminal penalties. Single-pair stem = 0. | `stem_energy.md` |
| `hairpin_loop_energy` | `RnaSecondaryStructure.CalculateHairpinLoopEnergy` | Size initiation + terminal mismatch; tri/tetra/hexaloop bonuses, all-C penalty, G-U closure adj. 3-nt loop = initiation only. | `hairpin_loop_energy.md` |
| `internal_loop_energy` | `RnaSecondaryStructure.CalculateInternalLoopEnergy` | `int11` lookup when `n1=n2=1`; else size init + asymmetry + terminal mismatches. | `internal_loop_energy.md` |
| `bulge_loop_energy` | `RnaSecondaryStructure.CalculateBulgeLoopEnergy` | n=1: init + stacking + special-C bonus − RT·ln(states); n>1: size init. | `bulge_loop_energy.md` |
| `multibranch_loop_energy` | `RnaSecondaryStructure.CalculateMultibranchLoopEnergy` | Affine `a + b·unpaired + c·helices` (a=9.25, b=0.91, c=−0.63) + optional stacking + strain +3.14. | `multibranch_loop_energy.md` |
| `dangling_end_energy` | `RnaSecondaryStructure.GetDanglingEndEnergy` | 5'/3' dangling-end stacking (NNDB) keyed `closing5+dangling+closing3`; `is3Prime` selects table; 0 if absent. | `dangling_end_energy.md` |
| `terminal_mismatch_energy` | `RnaSecondaryStructure.GetTerminalMismatchEnergy` | Terminal-mismatch stacking for closing pair + first mismatch (NNDB). Used by hairpin/internal/coaxial. | `terminal_mismatch_energy.md` |

## Envelope

- **RNA-STRUCT-001** — documented-only limitation (no runtime guard). Tertiary-stabilised knots are
  not recoverable as the MFE by any NN model; non-csr-PK pseudoknots (kissing / loop-loop /
  triple-crossing) are not faithfully detectable. See
  [`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md) and the STOP rule in
  [`../SKILL.md`](../SKILL.md).

## Not a tool here (route elsewhere)

- **Base-pair probabilities / partition function (McCaskill):** documented C#-only algorithm
  (`docs/algorithms/RnaStructure/RNA_Partition_Function.md`,
  `Turner_McCaskill_Partition_Function.md`) — **no MCP tool**. Use `predict_rna_structure` /
  `can_pair` / `base_pair_type` for pairing questions.
- **Pre-miRNA hairpins, miRNA seed/target pairing, RNA-seq, splicing:** → [`bio-annotation`](../../bio-annotation/SKILL.md).
