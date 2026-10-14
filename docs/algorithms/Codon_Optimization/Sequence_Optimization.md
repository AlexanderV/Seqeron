# Codon Optimization (Sequence Optimization)

| Field | Value |
|-------|-------|
| Algorithm Group | Codon Optimization |
| Test Unit ID | CODON-OPT-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-10-14 |

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
| Avoid rare codons | `AvoidRareCodons` (historical alias `AvoidRareCodeons`, same value) | Replace only codons whose current frequency is below the threshold. |
| Minimize secondary structure | `MinimizeSecondary` | `BalancedOptimization` codon selection and GC pass, then the `ReduceSecondaryStructure` sliding-window pass; scored by `secondaryScoring` (`PairCount` by default, `Turner2004Mfe` for the Turner 2004 energy model). |

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
| `geneticCode` | `GeneticCode?` | `null` | NCBI genetic code defining the encoded protein and the synonymous families (`null` = `GeneticCode.Standard`, bit-identical to the historical output). DNA Chisel `EnforceTranslation(genetic_table=…)` / `CodonOptimize(…)`, whose families are `get_backtranslation_table(genetic_table)` (Biopython `forward_table` grouped by amino acid, the table's `stop_codons` as `*`). Same parameter on `CalculateCAI`, `RemoveRestrictionSites`, `ReduceSecondaryStructure`, `FindRareCodons`, `CalculateMinMaxProfile`, `CreateCodonUsageTable`, `CreateCodonTableFromSequence` (review 2026-09, F48). | Any NCBI table of `GeneticCode.SupportedTableNumbers`. |
| `avoidSites` | `IReadOnlyList<string>?` | `null` | `null`/empty: no site avoidance (historical behaviour). A list of restriction recognition sequences (DNA or RNA, IUPAC codes allowed) runs `RemoveRestrictionSites(...)` with the same table and genetic code **after** codon selection, the GC pass and the structure pass — DNA Chisel composes `CodonOptimize(...)` with `AvoidPattern(EnzymeSitePattern(...))` + `EnforceTranslation()` in one problem and resolves the hard constraints last (review 2026-09, F54). | A site with a character that is not an IUPAC nucleotide code → `ArgumentException`. Applies to every strategy. |
| `gcWindow` | `int?` | `null` | `null`: the GC pass constrains the global GC fraction (historical behaviour). A positive value: every sliding window of that many nucleotides is constrained to `[gcTargetMin, gcTargetMax]`. | Must be positive; `ArgumentOutOfRangeException` otherwise. Used only during `BalancedOptimization` / `MinimizeSecondary`. |

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
7. If the strategy is `BalancedOptimization` or `MinimizeSecondary`, rebalance GC content (globally or per `gcWindow`).
8. If the strategy is `MinimizeSecondary`, run the `ReduceSecondaryStructure` pass.
9. If `avoidSites` is non-empty, run `RemoveRestrictionSites` on the result (hard constraint last).
10. Rebuild the change list from the original vs. final codons, recompute CAI and GC content, and return the final `OptimizationResult`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository strategy behaviors are:

| Strategy | Selection rule |
|----------|----------------|
| `MaximizeCAI` | Choose the synonymous codon with the highest table frequency (ties: first in NCBI codon order). Identical to DNA Chisel `CodonOptimize(method="use_best_codon")` and Biopython `CodonAdaptationIndex.optimize`.[8][9] |
| `AvoidRareCodeons` | If the current codon frequency is below `rareCodonThreshold`, replace it with the highest-frequency synonymous codon; otherwise keep the current codon. |
| `HarmonizeExpression` | Per amino acid, allocate the synonymous codons by largest-remainder rounding of `frequency × residue count`; positions already holding an allotted codon keep it. This is the deterministic optimum of DNA Chisel's `MatchTargetCodonUsage` objective.[8] |
| `BalancedOptimization` | As `MaximizeCAI`, then the GC pass. |
| `MinimizeSecondary` | Same as `BalancedOptimization` in `SelectOptimalCodon`, plus the `ReduceSecondaryStructure` pass on the GC-balanced sequence (default window 40 nt, scoring `secondaryScoring`). |

GC pass (`BalanceGcContent`): while the overall GC fraction is outside `[gcTargetMin, gcTargetMax]`, walk the codons and apply only swaps that move GC in the required direction, preferring a swap that lands the sequence inside the window over one that overshoots and, among equals, the most frequent codon; codons below `rareCodonThreshold` are not used; the pass stops as soon as the sequence is inside the window (DNA Chisel `EnforceGCContent(mini, maxi)` resolved by synonymous mutations).[8]

Genetic code: amino acids and synonymous families come from the `geneticCode` argument (`null` = `GeneticCode.Standard`), so a non-standard table changes the rewrite exactly as the code dictates — under NCBI table 2 (vertebrate mitochondrial) `AGA`/`AGG` are stops and are never used for Arg, `UGA` is Trp (family {UGA, UGG}) and `AUA` is Met (family {AUA, AUG}); under table 11 the codon table equals table 1's, so the output is identical to the Standard one. In the dual-coding tables 27/28/31 a context-dependent stop belongs to the family of the amino acid it encodes (`UAA`/`UAG` → Gln in table 27), exactly as in Biopython's `forward_table` and DNA Chisel's back-translation table; a caller that must keep a terminator in those tables has to exclude it itself, and a usage table built for the Standard code is not meaningful there (use `CreateCodonTableFromSequence(reference, name, code)`).

Ambiguity: amino acids come from that table, so an IUPAC-ambiguous triplet that still encodes one amino acid (e.g. `GCN` → Ala) is optimised like any other codon, while `NNN`, `B`/`Z`/`J` cases and triplets containing non-IUPAC symbols are carried through unchanged and translate to `X`.

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
- `CodonOptimizer.RemoveRestrictionSites(...)` (also chained inside `OptimizeSequence` via the optional `avoidSites`: IUPAC site matching on both strands, one synonymous codon change per removed occurrence, replacement chosen by host codon usage)
- `CodonOptimizer.CreateCodonTableFromSequence(...)` / `CodonOptimizer.CreateCodonUsageTable(...)` (reference table from a gene set, Sharp & Li 0.5 pseudo-count; table factory carrying the selected code's codon→amino-acid map, Standard by default)

### 5.2 Current Behavior

The result stores the normalized RNA input rather than the caller's original string. Stop codons are preserved unchanged. `Changes` and `ChangedCodons` are always rebuilt from the original vs. final codons, so they include the GC-balancing pass. `BalancedOptimization` and `MinimizeSecondary` share the same codon selection and the same GC pass; `MinimizeSecondary` then runs `ReduceSecondaryStructure` on the GC-balanced sequence, and the reported CAI, GC and `Changes` describe that final sequence (review 2026-09, F42).[7] With `avoidSites`, `RemoveRestrictionSites` runs last, after the GC and structure passes, and the reported CAI, GC and `Changes` describe that final sequence (review 2026-09, F54). Codon families and translation come from the canonical `GeneticCode` tables (no private genetic-code copy; `null` = Standard), in NCBI codon order, which is also the tie-break for equal-frequency codons.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Synonymous substitutions preserve the encoded protein while changing codon usage toward the target host.[1][2]
- `MaximizeCAI` uses codon frequencies to increase or maintain CAI relative to the original sequence.[1][7]

**Intentionally simplified:**

- GC balancing constrains the global GC fraction by default and, with `gcWindow = w`, every sliding window of `w` nucleotides — both modes of DNA Chisel's `EnforceGCContent(mini, maxi[, window])` (review 2026-09, F43). The windowed breach set, breach score and breach segments of `CodonOptimizer.EvaluateGcContent(...)` equal DNA Chisel 3.2.16's `EnforceGCContent(mini, maxi, window).evaluate(problem)`; the repair is a deterministic greedy synonymous-codon search rather than DNA Chisel's randomised constraint solver; **consequence:** the two libraries agree on which windows breach and by how much, but may land on different (equally valid) repaired sequences.
- `ReduceSecondaryStructure` defaults to the pair-count window score (Watson-Crick plus G·U wobble, via `RnaSecondaryStructure.CanPair`), a pairing-potential heuristic that ignores stacking, loop penalties and nesting; **consequence:** with the default scoring it relaxes obvious self-complementarity only. The thermodynamic alternative is opt-in: `scoring = Turner2004Mfe` scores every window with `RnaSecondaryStructure.CalculateMinimumFreeEnergy` (NNDB Turner 2004, the ViennaRNA RNAfold energy model) and maximises it. It is practical because window energies are memoised by window content and the fold count is capped by `maxFoldEvaluations` (default 20000), but it costs ~9.0 s per kilobase against ~0.06 s for the pair count (999-nt CDS, Release; 4.42 ms per 43-nt fold), which is why the heuristic remains the default of both the method and the `MinimizeSecondary` strategy.

- Restriction-site avoidance is composed into `OptimizeSequence` through the optional `avoidSites` list (review 2026-09, F54), mirroring DNA Chisel's single `DnaOptimizationProblem(constraints=[AvoidPattern(EnzymeSitePattern(...)), EnforceTranslation()], objectives=[CodonOptimize(...)])`: the hard site constraint is resolved **last**, after codon selection, the GC pass and the structure pass, so no later pass can re-introduce a site; **consequence:** the one codon changed per removed occurrence may perturb the GC fraction and the window structure that the preceding passes achieved, and an occurrence that no single synonymous substitution can remove is left in place (as for the standalone `RemoveRestrictionSites`). The repair is a deterministic greedy synonymous search rather than DNA Chisel's randomised solver, but on the three verification genes it lands on the byte-identical sequence (Evidence D-F54).

**Not implemented:**

- Codon-pair scoring / codon-pair bias inside `OptimizeSequence`; **users should rely on:** the single-codon strategies plus the composable `avoidSites`, `gcWindow` and `MinimizeSecondary` passes.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Public enum member was spelled `AvoidRareCodeons`. | Deviation | — | resolved (review 2026-09, F55) | The correctly spelled `AvoidRareCodons` member was added with the **same underlying value**, so both spellings are interchangeable (and both parse case-insensitively in the MCP layer); the old name is kept, un-attributed, because `TreatWarningsAsErrors` would turn an `[Obsolete]` attribute into build errors at the existing call sites. `ToString()` still returns `AvoidRareCodeons`, so no serialized value changes meaning. |
| 2 | `MinimizeSecondary` shares codon selection and the GC pass with `BalancedOptimization`, and adds the `ReduceSecondaryStructure` pass. | Design | The strategy is a distinct path: its output differs from `BalancedOptimization` whenever a window exceeds the structure threshold. | resolved (review 2026-09, F42) | Documented in [CODON-OPT-001.md](../../../tests/TestSpecs/CODON-OPT-001.md). |

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

`OptimizeSequence` is a codon-selection routine over a single codon-usage table. It does not model mRNA folding, codon-pair bias, ramp effects or regulatory context, and it optimises one objective at a time rather than solving a multi-constraint design problem the way a full sequence optimiser (DNA Chisel) does. Restriction-site removal runs inside `OptimizeSequence` only when `avoidSites` is supplied (it is still available standalone as `RemoveRestrictionSites`), and it is a last-pass hard constraint that can slightly perturb the GC fraction it follows; the structure pass (`ReduceSecondaryStructure`) runs only under the `MinimizeSecondary` strategy. Its thermodynamic scoring (`Turner2004Mfe`) is greedy and window-local: it maximises per-window Turner 2004 MFE, it does not minimise the folding energy of the whole transcript.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [CodonOptimizer_OptimizeSequence_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_OptimizeSequence_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Tests: [CodonOptimizer_GeneticCode_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_GeneticCode_Tests.cs) — optional `geneticCode` on every code-dependent API: default == explicit Standard, table 2 / 11 / 27 behaviour, DNA Chisel `match_codon_usage` cross-check under table 2
- Tests: [CodonOptimizer_SecondaryStructureScoring_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_SecondaryStructureScoring_Tests.cs) — Turner 2004 window energies vs ViennaRNA 2.7.2, window-MFE improvement, protein preservation, `MinimizeSecondary` structure pass
- Tests: [CodonOptimizer_AvoidSites_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_AvoidSites_Tests.cs) — the `avoidSites` pass (sites cleared on both strands, protein preserved, metrics of the final sequence, dnachisel-identical outputs) and the `AvoidRareCodons` spelling alias
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
9. Turner DH, Mathews DH. 2010. NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure. Nucleic Acids Research 38:D280-D282. Turner 2004 parameter set used by `RnaSecondaryStructure.CalculateMinimumFreeEnergy`. https://doi.org/10.1093/nar/gkp892
10. Lorenz R et al. 2011. ViennaRNA Package 2.0. Algorithms for Molecular Biology 6:26. Reference implementation 2.7.2 (`RNA.fold_compound(...).mfe()`, `md.dangles = 2`) used to cross-check the window free energies of the `Turner2004Mfe` scoring mode. https://doi.org/10.1186/1748-7188-6-26
