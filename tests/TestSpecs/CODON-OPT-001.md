# Test Specification: CODON-OPT-001

## Test Unit Information
- **ID:** CODON-OPT-001
- **Title:** Sequence Optimization
- **Canonical Method:** `CodonOptimizer.OptimizeSequence(...)`
- **Area:** Codon Optimization
- **Complexity:** O(n)
- **Status:** ☑ Complete

## Method Under Test

```csharp
public static OptimizationResult OptimizeSequence(
    string codingSequence,
    CodonUsageTable targetOrganism,
    OptimizationStrategy strategy = OptimizationStrategy.BalancedOptimization,
    double gcTargetMin = 0.40,
    double gcTargetMax = 0.60,
    double rareCodonThreshold = 0.15)
```

## Test Categories

### MUST Tests (Required for Completion)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| M1 | OptimizeSequence_PreservesProtein_AllStrategies | Optimized sequence must encode identical protein | Synonymous substitution definition |
| M2 | OptimizeSequence_EmptySequence_ReturnsEmptyResult | Empty input returns empty output with CAI=0 | Edge case definition |
| M3 | OptimizeSequence_ConvertsThymine_ToUracil | T is converted to U in RNA representation | RNA notation standard |
| M4 | OptimizeSequence_TrimsToCompleteCodons | Incomplete codons (length % 3 != 0) are trimmed | Codon definition |
| M5 | OptimizeSequence_MaximizeCAI_IncreasesOrMaintainsCAI | CAI after optimization >= CAI before | Sharp & Li (1987) |
| M6 | OptimizeSequence_SingleAminoAcidCodons_Unchanged | AUG (Met) and UGG (Trp) cannot be changed | Standard genetic code |
| M7 | OptimizeSequence_StopCodons_Preserved | Stop codons remain as stop codons | Translation termination |
| M8 | OptimizeSequence_DifferentOrganisms_DifferentResults | E. coli vs Yeast optimization differs | Organism-specific bias |
| M9 | OptimizeSequence_LowercaseInput_Handled | Case-insensitive processing | Robustness requirement |
| M10 | OptimizeSequence_ReturnsValidOptimizationResult | All result fields populated correctly | API contract |

### SHOULD Tests (Recommended)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| S1 | OptimizeSequence_AvoidRareCodons_OnlyReplacesRare | Only codons below threshold are changed | Strategy definition |
| S2 | OptimizeSequence_BalancedOptimization_AimsForTargetGc | GC content moves toward 40-60% range | Implementation spec |
| S3 | OptimizeSequence_TracksChanges_Correctly | Changes list accurately reflects modifications | API contract |
| S4 | OptimizeSequence_LongSequence_Completes | Performance acceptable for long sequences | Usability |
| S5 | OptimizeSequence_KnownSequence_VerifiedOutput | GFP or similar known sequence optimizes correctly | Integration test |

### COULD Tests (Optional)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| C1 | OptimizeSequence_HarmonizeExpression_MaintainsDistribution | Codon distribution matches host pattern | DNA Chisel `match_codon_usage` |

### Reference cross-checks (added review 2026-09)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| R1 | OptimizeSequence_MaximizeCAI_MatchesDnaChiselUseBestCodon | MaximizeCAI = DNA Chisel `CodonOptimize(method="use_best_codon")`, codon for codon | DNA Chisel 3.2.16 |
| R2 | OptimizeSequence_HarmonizeExpression_MatchesTargetUsageDeterministically | Deterministic largest-remainder allocation; same `match_codon_usage` score as DNA Chisel's optimizer (−11.77) | DNA Chisel 3.2.16 |
| R3 | OptimizeSequence_BalancedOptimization_MakesNoNeutralSwaps | GC pass only swaps codons that move GC toward the window and stops on entry | DNA Chisel `EnforceGCContent` |
| R4 | OptimizeSequence_AvoidRareCodons_NoSynonymAboveThreshold_UsesBestCodon | With no synonym above the threshold the best codon is still used | DNA Chisel `use_best_codon` |
| R5 | OptimizeSequence_AmbiguousCodons_HandledPerGeneticCode | GCN → best Ala codon; NNN and non-IUPAC triplets untouched, protein 'X' | GeneticCode / Biopython ambiguity |
| R6 | RemoveRestrictionSites_* (5 tests) | One codon per site, best-frequency substitution, both strands, IUPAC sites, unremovable site kept | REBASE/IUPAC, DNA Chisel `AvoidPattern` |
| R7 | CreateCodonTableFromSequence_MatchesBiopythonRelativeAdaptiveness | w from the built table = Biopython `CodonAdaptationIndex` (0.5 pseudo-count) | Biopython 1.88, Sharp & Li 1987 |

## Invariants to Verify

1. **Protein Preservation**: `original.ProteinSequence == optimized.ProteinSequence`
2. **CAI Improvement**: `MaximizeCAI strategy → OptimizedCAI >= OriginalCAI`
3. **Codon Length**: `OptimizedSequence.Length % 3 == 0`
4. **Same Length**: `OriginalSequence.Length == OptimizedSequence.Length`

## Edge Cases

| Case | Input | Expected Output |
|------|-------|-----------------|
| Empty | `""` | Empty result, CAI=0 |
| Single codon (Met) | `"AUG"` | Unchanged |
| Single codon (Trp) | `"UGG"` | Unchanged |
| DNA input | `"ATGGCT"` | Converts T→U, then optimizes |
| Incomplete codon | `"AUGGCUA"` (7 nt) | Trims to `"AUGGCU"` (6 nt) |
| All optimal | Already optimal E. coli | No changes, high CAI |
| All rare | All rare codons | Significant changes, CAI improves |

## Test Data

### Reference Sequences

```
# Short test sequence (M-A-Stop)
AUGGCUUAA → Protein: MA*

# E. coli rare codons for Arginine (Kazusa species=316407, W3110 K-12 substrain)
AGA, AGG → Frequency ~0.04, 0.02

# E. coli preferred codons for Leucine
CUG → Frequency 0.50

# Yeast preferred codons for Leucine (Kazusa species=4932)
UUA, UUG → Frequency 0.28, 0.29
```

### Known CAI Values (hand-verified against Sharp & Li formula)

- AUG alone: CAI = 1.0 (only codon for Met, wi = 1.0)
- CUGCCGACC (L-P-T, all optimal E. coli): CAI = 1.0
- CUAAGACGA (L-R-R, rare E. coli): CAI ≈ 0.106

## Audit Notes

### Existing Test Coverage Analysis

File: `CodonOptimizer_OptimizeSequence_Tests.cs`

| Region | Tests | Coverage | Assessment |
|--------|-------|----------|------------|
| Protein Preservation | 4+ | Good | All 5 strategies tested (incl. MinimizeSecondary) |
| CAI Behavior | 5 | Good | Covers increase/maintain, exact formula verification |
| Special Codons | 5 | Good | Met, Trp, Stop codons |
| Organism Specificity | 2 | Good | E. coli vs Yeast vs Human, exact optimized sequences per Kazusa |
| Input Handling | 2 | Good | Lowercase, exact result field values (hand-computed CAI, GC) |
| Strategy Specifics | 5 | Good | All strategies covered, AvoidRareCodons asserts replacement |
| Invariants | 6 | Good | Protein, CAI range, Sharp & Li formula, MaximizeCAI→1.0 |

### Coverage Classification Result (2026-03-10)

| Classification | Count | Details |
|---------------|-------|---------|
| ❌ Missing → Added | 1 | M1: MinimizeSecondary added to AllStrategies test |
| ⚠ Weak → Strengthened | 5 | M8: exact codon assertions per Kazusa; M10: exact hand-computed values; S2: GC enters target range; C1: CAI validity |
| 🔁 Duplicate → Removed | 19 | All OptimizeSequence duplicates in CodonOptimizerTests.cs |
| ✅ Covered | 24 | All remaining tests |

### Code Bug Fixed

- **BalancedOptimization Changes list**: `Changes` and `ChangedCodons` were not updated after GC content balancing phase. Fixed by rebuilding changes from original vs final codons.

### Data Source Traceability

| Organism | Kazusa Species ID | Kazusa Name | Dataset Size |
|----------|------------------|-------------|--------------|
| E. coli K12 | 316407 | E. coli W3110 (K-12 substrain) | 4332 CDS |
| S. cerevisiae | 4932 | Saccharomyces cerevisiae | 14411 CDS |
| H. sapiens | 9606 | Homo sapiens | 93487 CDS |

---

## Deviations and Assumptions

- **BalancedOptimization Changes rebuild (fixed 2026-03-10)**: Previously, `Changes` list only reflected the initial optimization pass, missing GC content balancing modifications. Fixed to rebuild changes by comparing original vs final codons.
- **Codon usage tables**: All three tables (E. coli, Yeast, Human) verified against Kazusa Codon Usage Database raw data (per-thousand frequencies → relative fractions per amino acid).
- **CAI formula**: Matches Sharp & Li (1987) definition: w_i = f_i / max(f_j), CAI = exp((1/L)·Σ ln(w_i)). Met/Trp and stops are not scored; w < 0.0001 → 0.01 (CodonW `cai_out`) — review 2026-09, CODON-CAI-001 F12/F13 (formerly Met/Trp scored with w = 1 and a 1e-6 clamp).
- **Standard genetic code**: All 64 codons verified correct.
- **Optimization strategies**: All thresholds exposed as configurable parameters (`rareCodonThreshold`, `gcTargetMin`, `gcTargetMax`); no hardcoded assumptions.
- **ReduceSecondaryStructure window width (review 2026-09, F64)**: the sliding window is the smallest whole-codon window at least `windowSize` nucleotides wide — ceil(`windowSize` / 3) codons. The default 40 is 14 codons = 42 nt (unchanged); before F64 the width was `windowSize / 3 + 1` codons, which over-rounded by one codon whenever `windowSize` was a multiple of 3 (42 was scored over 45 nt), contradicting the documented contract. A sequence shorter than `windowSize` is still returned normalised but otherwise unchanged (its trailing partial codon is kept — documented, deliberate).
- **MinimizeSecondary**: balanced codon selection + GC pass + the dedicated `ReduceSecondaryStructure` pass (review 2026-09, F42; the strategy used to stop after the GC pass, i.e. it was an alias of BalancedOptimization). The pass is scored by `secondaryScoring`: `PairCount` (default) or `Turner2004Mfe`. The structure pass carries the GC constraint of the GC pass and rejects any rewrite that would increase the GC breach (review 2026-09, F63; dnachisel resolves `EnforceGCContent` as a hard constraint and accepts objective mutations only while the constraints hold).

### Review 2026-09 (CODON-OPT-001, F21–F25)

- **Genetic code**: the private RNA-keyed `StandardGeneticCode` / `AminoAcidToCodons` copies were removed; amino acids and synonymous families now come from `GeneticCode.Standard` (NCBI table 1). Codon families are in NCBI order, which is the documented, deterministic tie-break for every "most frequent synonymous codon" choice. IUPAC-ambiguous triplets resolve through `GeneticCode.Translate` (GCN → Ala → best Ala codon); a triplet with non-IUPAC symbols is left untouched and contributes `X` to the protein.
- **Strategies**: `MaximizeCAI` = DNA Chisel `use_best_codon` (verified identical output). `AvoidRareCodons` (historical spelling `AvoidRareCodeons`, F55) / `BalancedOptimization` no longer keep a rare codon when *no* synonym reaches the threshold — the most frequent synonym is used. `HarmonizeExpression` is no longer weighted-random (`new Random()` per call): it is the deterministic largest-remainder match of the target codon-usage profile (DNA Chisel `match_codon_usage`), which reaches the same DNA Chisel objective score as that library's own optimizer.
- **GC balancing**: only swaps that move GC toward `[gcTargetMin, gcTargetMax]` are applied (no neutral swaps), a candidate that lands inside the window is preferred over one that overshoots, the direction is re-evaluated each step and the pass stops as soon as the sequence is inside the window; the replacement frequency floor is the caller's `rareCodonThreshold` instead of a hard-coded 0.1; the GC count is maintained incrementally (was an O(n²) full rescan).
- **RemoveRestrictionSites**: exactly one codon changes per removed occurrence (the old loop kept rewriting the following codons after the site was already gone), the substitution is the most frequent synonymous codon in the supplied table (the table was previously ignored), sites are matched with IUPAC semantics and on both strands, and an unremovable occurrence no longer aborts the remaining occurrences.
- **ReduceSecondaryStructure**: base pairing delegates to `RnaSecondaryStructure.CanPair` (ViennaRNA pair set, so the G·U wobble counts); the window baseline is re-evaluated after each accepted change; output is always normalised RNA.

### Review 2026-09 (CODON-OPT-001, F41–F42 — thermodynamic structure scoring)

- **`SecondaryStructureScoring` (F41, audit A1-6)**: `ReduceSecondaryStructure` gained an opt-in `scoring` mode. `PairCount` (default, unchanged) scores a window by the fraction of position pairs (i, j), j ≥ i + 4, that can form a canonical pair. `Turner2004Mfe` scores a window by its minimum free energy from the canonical engine `RnaSecondaryStructure.CalculateMinimumFreeEnergy` (NNDB Turner 2004 nearest-neighbour parameters, the energy model of ViennaRNA RNAfold) and **maximises** it (less negative = less structure); a window is rewritten when its MFE is below `mfeThreshold` (default 0 kcal/mol) and a synonymous swap is accepted only when it raises the window MFE. Window energies are memoised by window content and the pass stops after `maxFoldEvaluations` distinct folds (default 20000), which bounds the runtime and keeps the output deterministic.
- **Reference numbers**: our window MFE vs ViennaRNA 2.7.2 `RNA.fold_compound(window, md).mfe()` with `md.dangles = 2` (RNAfold default) on six 45-nt windows of the two test genes: −25.82/−25.80, −29.34/−29.00, −25.82/−26.50, −9.22/−9.20, −13.44/−14.00, −13.82/−14.00 kcal/mol — max deviation 0.68 kcal/mol, within the engine's documented residual (RNA-STRUCT-001 Evidence D3/D4, missing NNDB int21/int22 lookup tables). Over 40 random 20–60-nt windows: mean |Δ| 0.50, max 2.24 kcal/mol, identical MFE structure in 30/40.
- **Effect**: on the designed 111-nt gene A the Turner pass raises the summed 45-nt window MFE from −477.92 to −144.26 kcal/mol (worst window −37.34 → −16.04); on the 99-nt gene B from −252.88 to −96.45 (worst window −17.00 → −7.18). Protein preserved in both; exact outputs locked.
- **Runtime** (999-nt CDS, Release, this machine): `PairCount` 0.06 s, `Turner2004Mfe` 9.02 s (single-fold cost 4.42 ms per 43-nt window). `PairCount` therefore stays the default of both `ReduceSecondaryStructure` and the MinimizeSecondary strategy: ~150× cheaper and the historical behaviour; `Turner2004Mfe` is the sourced energy model and is opt-in per call.
- **DNA Chisel**: has no free-energy objective to mirror — `AvoidHairpins` is the IDT reverse-complement stem rule, and dnachisel 3.2.16 contains no ViennaRNA / `RNA.fold` / free-energy specification (verified by grep over the installed package). The energy model is therefore sourced from ViennaRNA RNAfold / NNDB Turner 2004 through this repository's own folding engine.
- **MinimizeSecondary (F42, audit A1-8)**: the strategy now runs the structure pass after the GC pass. On a CDS of two complementary blocks (8×Glu + 8×Phe) BalancedOptimization yields `GAGGAGGAGGAGGAGGAGGAGGAGUUCUUCUUCUUCUUUUUUUUUUUU` while MinimizeSecondary yields `GAGGAGGAGGAGGAGGAGGAGGAGUUCUUCUUCUUCUUCUUCUUCUUU` (protein `EEEEEEEEFFFFFFFF` preserved; before F63 it yielded `GAGGAAGAAGAAGAAGAAGAAGAAUUCUUCUUCUUCUUCUUCUUCUUU`, GC 0.3333, which breached the strategy's own GC band by 0.0667). The reported CAI/GC and `Changes` describe the sequence actually returned.
- **CreateCodonTableFromSequence**: codons absent from the reference set get the Sharp & Li (1987) count of 0.5 ("following the description in the original paper", Biopython `CodonAdaptationIndex` 1.88), so all 64 codons are present and the derived relative adaptiveness equals Biopython's index exactly. Previously an absent codon was missing from the table and was scored by `CalculateCAI` with the CodonW zero substitute 0.01.

### Review 2026-09 (CODON-OPT-001, F43 — windowed GC constraint, audit A1-7)

- **`gcWindow` (additive)**: `OptimizeSequence(..., int? gcWindow = null)` — `null` (default) keeps the historical global GC pass; a positive window size constrains **every** sliding window of that many nucleotides to `[gcTargetMin, gcTargetMax]`, i.e. DNA Chisel `EnforceGCContent(mini, maxi, window=…)`. MCP: optional `gc_window` on `optimize_codons`.
- **`CodonOptimizer.EvaluateGcContent(sequence, mini, maxi, window, locationsSpan = 50)`**: the DNA Chisel evaluation itself — window `i` spans `[i, i+w)` for `i = 0 … len-w` (`len-w+1` windows; a window longer than the sequence yields none, score 0), per-window breach `max(0, mini-gc) + max(0, gc-maxi)` with **inclusive** bounds, `Score` = −Σ breaches, and breach segments grouped as DNA Chisel's `group_nearby_segments(max_start_spread = max(1, locationsSpan))` → `[first start, last start + w)`. Global mode (`window = null`) reports the whole sequence as the breach segment.
- **Reference numbers** (dnachisel 3.2.16 `EnforceGCContent(...).evaluate(problem)`, bounds 0.40/0.60, window 20): AT-rich/GC-rich/mixed/two-segment test genes score −13.85 / −9.45 / −4.25 / −20.10; breach window starts `[0..7, 13..27, 33, 34]` (mixed) and `[0..31, 125..157]` (two-segment); segments `0-54` and `0-51, 125-177`; global mode −0.262745098039 (GC-rich) and −0.381481481481 (AT-rich). Our values are equal to all of these.
- **Repair**: deterministic greedy — the left-most breaching window is repaired by the synonymous replacement (frequency ≥ `rareCodonThreshold`) of an overlapping codon that lowers the whole-sequence breach total the most, ties broken by codon frequency, then position, then NCBI order; a window whose bounds are unreachable with synonymous codons is skipped (no exception, protein intact). After the pass dnachisel's own `evaluate` scores every locked output **0** on the requested window size, while the global-only default still scores below 0 on 20-nt windows.

### Review 2026-09 (CODON-OPT-001, F48 — genetic-code-aware API, audit A2-1)

- **Additive optional `GeneticCode? geneticCode = null`** on `OptimizeSequence`, `CalculateCAI`, `RemoveRestrictionSites`, `ReduceSecondaryStructure`, `FindRareCodons` (both overloads), `CalculateMinMaxProfile` (both overloads), `CreateCodonUsageTable` and `CreateCodonTableFromSequence`. `null` = `GeneticCode.Standard`, so every existing call site compiles and every existing expectation is bit-identical (locked by `CodonOptimizer_GeneticCode_Tests` default-vs-explicit-Standard assertions over four strategies).
- **Families**: one grouping of `GeneticCode.CodonTable` per NCBI table (no private tables) — the grouping DNA Chisel uses (`get_backtranslation_table(genetic_table)` = Biopython `forward_table` grouped by amino acid, `stop_codons` as `*`) and the one the canonical `CodonUsageAnalyzer` core uses for RSCU/CAI/ENC. Verified sets: table 2 R = CGU CGC CGA CGG, M = AUA AUG, W = UGA UGG, * = UAA UAG AGA AGG; table 11 is identical to table 1 (`forward_table` and `stop_codons` equal, only starts differ), so its output equals the Standard one; table 27 Q = UAA UAG CAA CAG, W = UGA UGG, * = UGA.
- **Locked behaviour**: `AUGCGUAGACGCAGGUGAAUAUAA` with `MaximizeCAI`/E. coli K12 → `AUGCGCCGCCGCCGCUGAAUUUAA` (table 1, protein `MRRRR*I*`) and `AUGCGCAGACGCAGGUGGAUGUAA` (table 2, protein `MR*R*WM*`: AGA/AGG kept as stops, UGA → UGG as Trp, AUA → AUG as Met). With an AGA-dominant usage table, `MaximizeCAI` yields `AUGAGAAGAAGAAGAUAA` under table 1 and `AUGCGCCGCCGCCGCUAA` under table 2 — no AGA/AGG is ever emitted for Arg under table 2. `HarmonizeExpression` on `CGU`×10 under table 2 → `CGUCGUCGUCGUCGCCGCCGCCGCCGACGG`, which DNA Chisel's `match_codon_usage` objective scores **−0.6**, the score its own randomised optimizer reaches (input −11.8) — the F23 cross-check repeated under a non-standard table.
- **CAI**: `CalculateCAI(seq, table, geneticCode: code)` equals `CodonUsageAnalyzer.CalculateCai(seq, table, code)` (already CodonW-locked) for every table — asserted directly in the tests.
- **Dual-coding tables 27/28/31**: treated per Biopython / DNA Chisel — a context-dependent stop sits in the family of the amino acid it encodes, so under table 27 `CAA`/`CAG`/`UAA`/`UAG` are one Gln family and a rewrite may move between them (locked: `AUGCAACAGUGA` → `AUGUAAUAAUGG`, protein `MQQW` under table 27). No codon unambiguously terminates in those tables, so a caller that must keep a terminator excludes it itself (`Translator.FindOrfs` rejects these tables for the same reason).
- **Code-independent by construction** (no parameter offered, documented): `FindRareCodonClusters`, `FindConservedRareCodonClusters`, `CalculateCodonUsage`, `CompareCodonUsage` — pause marking and codon counting are per-codon table lookups with no translation and no family.

### Review 2026-09 (CODON-OPT-001, F54/F55 — chained site avoidance, strategy spelling; audit A2-13 / A2-14)

- **`avoidSites` (additive, F54)**: `OptimizeSequence(..., IReadOnlyList<string>? avoidSites = null)` — `null`/empty keeps the historical output; a list of recognition sequences (DNA or RNA, IUPAC codes allowed) runs `RemoveRestrictionSites` with the same table and genetic code **after** codon selection, the GC pass and the structure pass, mirroring DNA Chisel's single `DnaOptimizationProblem(constraints=[AvoidPattern(EnzymeSitePattern(...)), EnforceTranslation()], objectives=[CodonOptimize(...)])`, which resolves the hard constraints last. MCP: optional `avoid_sites` on `optimize_codons` (recognition sequences, exactly like `remove_restriction_sites`). The reported CAI, GC and `Changes` describe the final sequence (as F42).
- **Reference numbers (F54)**: codon selection itself creates sites with the E. coli K12 table (Leu `CUG` + Gln `CAG` = PstI `CTGCAG`, His `CAU` + Met `AUG` = NdeI `CATATG`, Thr `ACC` + Met `AUG` = NcoI `CCATGG`). On three genes (A `MKLQALHMTMAEF*` 42 nt, B GFP N-terminus + `LQHMTMA` 150 nt, C `MLQHMTMAK*` 30 nt) the plain `MaximizeCAI` output fails all three dnachisel 3.2.16 `AvoidPattern(...).evaluate` checks (score −1 each) while `EnforceTranslation` passes; with `avoidSites` every check scores **0 (passes)** and the outputs are `AUGAAACUGCAAGCGCUGCACAUGACGAUGGCGGAAUUUUAA` (CAI 0.870294), `…CUGCAACACAUGACGAUGGCGUAA` (0.970251) and `AUGCUGCAACACAUGACGAUGGCGAAAUAA` (0.793311) — all three byte-identical to dnachisel's own `resolve_constraints() + optimize()` solution, same `CodonOptimize` objective score −1.3892. Protein preserved; the pass is idempotent; both strands are cleared (F22 property, re-checked through the new path for four strategies).
- **`AvoidRareCodons` (F55)**: the correctly spelled member was added with the **same underlying value** as the historical `AvoidRareCodeons`, which is kept un-attributed by policy: additive alias only — the misspelt member stays supported, no deprecation (campaign backward-compatibility rule). Both spellings behave identically and both parse case-insensitively through `System.Text.Json`'s enum converter, i.e. in the MCP layer. `ToString()` still returns `AvoidRareCodeons`; `Enum.GetValues<OptimizationStrategy>()` returns 6 entries / 5 distinct values, so the loops in `CodonProperties` and `CodonMetamorphicTests` simply run that strategy twice.
- Tests: `CodonOptimizer_AvoidSites_Tests` (8 cases), `OptimizeCodonsTests.OptimizeCodons_AvoidSites_ClearsSitesAfterTheOtherPasses` / `OptimizeCodons_AvoidRareCodons_BothSpellingsAccepted`.

### Review 2026-09 (CODON-OPT-001, F60 — GC-constraint-aware `avoidSites` pass; audit A3-6)

- **Behaviour**: for the strategies that declare a GC target (`BalancedOptimization`, `MinimizeSecondary`) the `avoidSites` pass carries the GC constraint the GC pass enforced — `gcTargetMin`/`gcTargetMax` globally, and every `gcWindow` window with `EvaluateGcContent` semantics — because DNA Chisel resolves `AvoidPattern` and `EnforceGCContent` as hard constraints of ONE problem. Site-removing synonymous candidates are ranked lexicographically: (a) GC-constraint-satisfying (or, when the input already breaches, not increasing the breach score), (b) codon frequency, (c) left-most codon, (d) NCBI codon order. Internal only: the public `RemoveRestrictionSites` overloads stay GC-blind, and `MaximizeCAI`/`AvoidRareCodons`/`HarmonizeExpression` thread no constraint, so the F54 locked outputs are byte-identical.
- **Documented fallback**: if no site-removing candidate satisfies the GC constraint, the least-breach candidate is applied and the site is still removed — `AvoidPattern` stays the hard constraint and nothing throws (dnachisel's `resolve_constraints()` raises `NoSolutionError` for such a problem). The structure pass runs *before* the site pass and is not re-run, so the site pass is GC-aware only.
- **Reference numbers (dnachisel 3.2.16)**: four constructed genes, `constraints=[AvoidPattern("CTGCAG"/"CATATG"/"CCATGG"), EnforceGCContent(0.40, 0.60[, window=24]), EnforceTranslation()]`. The GC-blind chain fails `EnforceGCContent` with score −0.225 / −0.225 / −0.150 (window 24) and −0.019048 (global) while clearing the sites; the GC-aware output passes **every** constraint and is byte-identical to dnachisel's own `resolve_constraints() + optimize()` result: `ATGATTCGCCTGTGCACTCCAGCTTATGGTCGTCGTCGCTAA` (objective −3.14), `ATGCGTATGGCACTGCGTTCTGCTGTTGGTTGGCTGGTGTAA` (−2.62), `ATGGGCAGCGAAACCTCTGCTTATGGTGCAGCGCTGTGGTAA` (−2.16), `ATGGTTACTGCTGTTCCAATGCGCGGCATGGGCCGCCGCTAA`. Unreachable bounds (gene 3, `EnforceGCContent(0.30, 0.45, window=24)`, dnachisel `NoSolutionError`): start −0.483333, GC-blind chain −1.066667, ours −0.483333.
- Tests: `CodonOptimizer_AvoidSites_Tests` GC region (5 cases: the defect through the public GC-blind overload with locked breach values, windowed and global GC-aware locked outputs, the unreachable-bounds fallback, non-GC strategies unchanged).

## Date
2026-03-10 (reviewed 2026-09-28, campaign review-2026-09 / B02; F54/F55/F60 2026-10-14)

