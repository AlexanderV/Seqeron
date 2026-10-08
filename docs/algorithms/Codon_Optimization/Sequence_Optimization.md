# Codon Optimization (Sequence Optimization)

| Field | Value |
|-------|-------|
| Algorithm Group | Codon Optimization |
| Test Unit ID | CODON-OPT-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Codon optimization replaces codons with synonymous alternatives that are better suited to a target host while preserving the encoded protein sequence.[1][2] In this repository, `CodonOptimizer.OptimizeSequence` normalizes the input to RNA notation, trims incomplete codons, applies a strategy-specific synonymous-codon selection policy, optionally balances GC content, and returns the optimized sequence together with CAI, GC content, and explicit codon changes. The method is linear in sequence length under the repository's fixed-size codon tables and is designed for heuristic host adaptation rather than full mRNA design.[7]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Codon optimization is commonly used for heterologous expression because organisms prefer different synonymous codons, reflecting tRNA pools and broader codon-usage bias.[1][2] The original document preserved these example host preferences:

| Organism | Example preferred codon | Frequency | Example rare codon | Frequency |
|----------|--------------------------|-----------|---------------------|-----------|
| E. coli K12 | `CUG` (Leu) | `0.50` | `CUA` (Leu) | `0.04` |
| E. coli K12 | `CGC` (Arg) | `0.40` | `AGG` (Arg) | `0.02` |
| S. cerevisiae | `UUG` (Leu) | `0.29` | `CUC` (Leu) | `0.06` |
| S. cerevisiae | `AGA` (Arg) | `0.48` | `CGG` (Arg) | `0.04` |

### 2.2 Core Model

The optimization goal is to preserve the amino-acid sequence while changing codons according to one of the public strategy options:

| Strategy | Current API name | Repository behavior |
|----------|------------------|---------------------|
| Maximum CAI | `MaximizeCAI` | Choose the highest-frequency synonymous codon. |
| Balanced optimization | `BalancedOptimization` | Prefer non-rare codons, then rebalance GC toward the configured range. |
| Harmonized expression | `HarmonizeExpression` | Deterministic largest-remainder match of the target codon-usage profile (DNA Chisel `match_codon_usage`).[8] |
| Avoid rare codons | `AvoidRareCodeons` | Replace only codons whose current frequency is below the threshold. |
| Minimize secondary structure | `MinimizeSecondary` | Same codon selection and GC pass as `BalancedOptimization` inside `OptimizeSequence`; the dedicated pass is `ReduceSecondaryStructure`. |

CAI is used as one of the optimization metrics and is defined as in [CAI_Calculation.md](CAI_Calculation.md).[1]

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Synonymous substitutions preserve the relevant biological product. | Codon changes may be acceptable mathematically but biologically inappropriate for the use case. |
| ASM-02 | The supplied codon-usage table represents the desired host preferences. | Optimization can move the sequence toward the wrong codon profile. |
| ASM-03 | GC-content heuristics between `gcTargetMin` and `gcTargetMax` are sufficient for the use case. | The result may satisfy the heuristic but still be suboptimal for structure or expression. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The optimized sequence preserves the translated protein sequence. | Replacement codons are selected from the synonymous set of the original amino acid and stop codons are preserved. |
| INV-02 | `OptimizedSequence.Length % 3 == 0`. | The method trims to complete codons before optimization. |
| INV-03 | `OriginalSequence.Length == OptimizedSequence.Length`. | The result stores the trimmed normalized input and a codon-for-codon optimized sequence of the same length. |
| INV-04 | `MaximizeCAI` yields `OptimizedCAI >= OriginalCAI`. | The method replaces each codon with the highest-frequency synonymous codon in the supplied table. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `codingSequence` | `string` | required | DNA or RNA coding sequence to optimize. | Empty input returns an empty `OptimizationResult`. |
| `targetOrganism` | `CodonUsageTable` | required | Target codon-usage table. | Must contain the host frequencies used for selection. |
| `strategy` | `OptimizationStrategy` | `BalancedOptimization` | Strategy controlling synonymous replacement. | Public enum values are listed in Section 2.2. |
| `gcTargetMin` | `double` | `0.40` | Lower bound for GC balancing. | Used only during `BalancedOptimization`. |
| `gcTargetMax` | `double` | `0.60` | Upper bound for GC balancing. | Used only during `BalancedOptimization`. |
| `rareCodonThreshold` | `double` | `0.15` | Threshold used by rare-codon-aware branches. | Compared against table frequencies. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `OriginalSequence` | `string` | Normalized and trimmed RNA input. |
| `OptimizedSequence` | `string` | Final optimized RNA sequence. |
| `ProteinSequence` | `string` | Protein translation accumulated from the original codons. |
| `OriginalCAI` | `double` | CAI of the normalized trimmed input. |
| `OptimizedCAI` | `double` | CAI of the final optimized sequence. |
| `GcContentOriginal` | `double` | GC fraction of the normalized trimmed input. |
| `GcContentOptimized` | `double` | GC fraction of the optimized sequence. |
| `ChangedCodons` | `int` | Number of codon changes after any GC-balancing phase. |
| `Changes` | `IReadOnlyList<(int Position, string Original, string Optimized)>` | Nucleotide-position change list. |

### 3.3 Preconditions and Validation

The method uppercases the sequence, converts `T` to `U`, and trims any trailing incomplete codon. Empty input yields an all-empty result with zero-valued metrics. Stop codons are preserved, and single-codon amino acids such as Met and Trp remain unchanged because no synonymous alternatives exist.[7]

## 4. Algorithm

### 4.1 High-Level Steps

1. Return an empty result for `null` or empty input.
2. Normalize the sequence to uppercase RNA notation and trim incomplete trailing bases.
3. Split the sequence into codons and compute `OriginalCAI`.
4. Translate each codon to an amino acid and preserve stop codons.
5. Select a synonymous replacement according to the chosen strategy.
6. Build the optimized sequence and record codon changes.
7. If the strategy is `BalancedOptimization`, rebalance GC content and rebuild the change list from the final codons.
8. Recompute CAI and GC content and return the final `OptimizationResult`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository strategy behaviors are:

| Strategy | Selection rule |
|----------|----------------|
| `MaximizeCAI` | Choose the synonymous codon with the highest table frequency (ties: first in NCBI codon order). Identical to DNA Chisel `CodonOptimize(method="use_best_codon")` and Biopython `CodonAdaptationIndex.optimize`.[8][9] |
| `AvoidRareCodeons` | If the current codon frequency is below `rareCodonThreshold`, replace it with the highest-frequency synonymous codon; otherwise keep the current codon. |
| `HarmonizeExpression` | Per amino acid, allocate the synonymous codons by largest-remainder rounding of `frequency × residue count`; positions already holding an allotted codon keep it. This is the deterministic optimum of DNA Chisel's `MatchTargetCodonUsage` objective.[8] |
| `BalancedOptimization` | As `MaximizeCAI`, then the GC pass. |
| `MinimizeSecondary` | Same as `BalancedOptimization` in `SelectOptimalCodon`. |

GC pass (`BalanceGcContent`): while the overall GC fraction is outside `[gcTargetMin, gcTargetMax]`, walk the codons and apply only swaps that move GC in the required direction, preferring a swap that lands the sequence inside the window over one that overshoots and, among equals, the most frequent codon; codons below `rareCodonThreshold` are not used; the pass stops as soon as the sequence is inside the window (DNA Chisel `EnforceGCContent(mini, maxi)` resolved by synonymous mutations).[8]

Ambiguity: amino acids come from `GeneticCode.Standard`, so an IUPAC-ambiguous triplet that still encodes one amino acid (e.g. `GCN` → Ala) is optimised like any other codon, while `NNN`, `B`/`Z`/`J` cases and triplets containing non-IUPAC symbols are carried through unchanged and translate to `X`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `OptimizeSequence` | `O(n)` | `O(n)` | The repository test spec classifies the algorithm as linear in the sequence length under fixed codon tables.[7] |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CodonOptimizer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs)

- `CodonOptimizer.OptimizeSequence(...)`
- `CodonOptimizer.SelectOptimalCodon(...)` (private helper)
- `CodonOptimizer.BalanceGcContent(...)` (private helper)
- `CodonOptimizer.ReduceSecondaryStructure(...)` (separate public helper not called from `OptimizeSequence`)
- `CodonOptimizer.RemoveRestrictionSites(...)` (separate public helper: IUPAC site matching on both strands, one synonymous codon change per removed occurrence, replacement chosen by host codon usage)
- `CodonOptimizer.CreateCodonTableFromSequence(...)` / `CodonOptimizer.CreateCodonUsageTable(...)` (reference table from a gene set, Sharp & Li 0.5 pseudo-count; table factory carrying the Standard-code codon→amino-acid map)

### 5.2 Current Behavior

The result stores the normalized RNA input rather than the caller's original string. Stop codons are preserved unchanged. `Changes` and `ChangedCodons` are always rebuilt from the original vs. final codons, so they include the GC-balancing pass. `BalancedOptimization` and `MinimizeSecondary` share the same codon selection and the same GC pass; a separate `ReduceSecondaryStructure` helper exists for direct structure-reduction work and is not invoked automatically by `OptimizeSequence`.[7] Codon families and translation come from `GeneticCode.Standard` (no private genetic-code copy), in NCBI codon order, which is also the tie-break for equal-frequency codons.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Synonymous substitutions preserve the encoded protein while changing codon usage toward the target host.[1][2]
- `MaximizeCAI` uses codon frequencies to increase or maintain CAI relative to the original sequence.[1][7]

**Intentionally simplified:**

- GC balancing works on the global GC fraction rather than on sliding GC windows (DNA Chisel supports both); **consequence:** a sequence inside the global window may still contain GC-extreme stretches.
- `MinimizeSecondary` does not invoke a secondary-structure pass inside `OptimizeSequence`; **consequence:** callers must call `ReduceSecondaryStructure` for that.
- `ReduceSecondaryStructure` scores a window by counting positions that can form a canonical pair (Watson-Crick plus G·U wobble, via `RnaSecondaryStructure.CanPair`) rather than by folding free energy; **consequence:** it relaxes obvious self-complementarity only. A free-energy-guided search would need `RnaSecondaryStructure.CalculateMinimumFreeEnergy` per candidate codon — O(w³) each (≈6 ms per 43-nt window measured here, i.e. minutes per kilobase) — and is left to the caller.

**Not implemented:**

- Automatic invocation of `ReduceSecondaryStructure` as part of `OptimizeSequence`; **users should rely on:** `CodonOptimizer.ReduceSecondaryStructure(...)` when a separate structure pass is required.
- Constraint sets such as restriction-site removal or codon-pair scoring inside `OptimizeSequence`; **users should rely on:** separate repository APIs such as `RemoveRestrictionSites(...)` where applicable.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Public enum member is spelled `AvoidRareCodeons`. | Deviation | Documentation and callers must use the API spelling, not the natural-language spelling. | accepted | Confirmed in [CodonOptimizer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs). |
| 2 | `MinimizeSecondary` shares codon selection and the GC pass with `BalancedOptimization` inside `OptimizeSequence`. | Deviation | This strategy name does not currently provide a distinct optimization path in that method. | accepted | Documented in [CODON-OPT-001.md](../../../tests/TestSpecs/CODON-OPT-001.md). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns an empty `OptimizationResult` with zero metrics. | Explicit early return. |
| DNA input | Converted to RNA notation before optimization. | The method normalizes `T` to `U`. |
| Incomplete final codon | Trimmed away. | Optimization works only on complete codons. |
| Single-codon amino acids | Remain unchanged. | No synonymous alternatives exist for Met or Trp. |
| Stop codons | Preserved. | The loop keeps stop codons unchanged. |

### 6.2 Limitations

`OptimizeSequence` is a codon-selection routine over a single codon-usage table. It does not model mRNA folding, codon-pair bias, ramp effects or regulatory context, and it optimises one objective at a time rather than solving a multi-constraint design problem the way a full sequence optimiser (DNA Chisel) does. Restriction-site removal and structure reduction are separate passes (`RemoveRestrictionSites`, `ReduceSecondaryStructure`) and are not run automatically.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [CodonOptimizer_OptimizeSequence_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_OptimizeSequence_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [CODON-OPT-001.md](../../../tests/TestSpecs/CODON-OPT-001.md)
- Related algorithms: [CAI_Calculation.md](CAI_Calculation.md), [Rare_Codon_Detection.md](Rare_Codon_Detection.md), [Codon_Usage_Analysis.md](Codon_Usage_Analysis.md)

## 8. References

1. Sharp PM, Li WH. 1987. The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research. N/A
2. Plotkin JB, Kudla G. 2011. Synonymous but not the same: the causes and consequences of codon bias. Nature Reviews Genetics. N/A
3. Mignon C et al. 2018. Codon harmonization - going beyond the speed limit for protein expression. FEBS Letters. N/A
4. Athey J et al. 2017. A new and updated resource for codon usage tables. BMC Bioinformatics. N/A
5. Kazusa Codon Usage Database. E. coli K-12 substr. W3110, species 316407. https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=316407
6. Kazusa Codon Usage Database. Saccharomyces cerevisiae, species 4932. https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932
7. Test specification: [CODON-OPT-001.md](../../../tests/TestSpecs/CODON-OPT-001.md)
8. Zulkower V, Rosser S. 2020. DNA Chisel, a versatile sequence optimizer. Bioinformatics 36(16):4508-4509. Reference implementation 3.2.16: `CodonOptimize` (`use_best_codon` / `match_codon_usage` / `harmonize_rca`), `MatchTargetCodonUsage.codon_usage_matching_stats`, `EnforceGCContent`, `AvoidPattern`. https://edinburgh-genome-foundry.github.io/DnaChisel/
9. Biopython 1.88, `Bio.SeqUtils.CodonAdaptationIndex` — builds the reference index from a gene set with a count of 0.5 for codons absent from the reference set ("following the description in the original paper", Sharp & Li 1987) and `optimize()` for the most-preferred-codon rewrite.
10. Puigbò P, Guzmán E, Romeu A, Garcia-Vallvé S. 2007. OPTIMIZER: a web server for optimizing the codon usage of DNA sequences. Nucleic Acids Research 35:W126-W131 ("one amino acid – one codon" method).
