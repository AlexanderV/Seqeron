# Rare Codon Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Codon Optimization |
| Test Unit ID | CODON-RARE-001 |
| Related Projects | N/A |
| Implementation Status | Production |
| Last Reviewed | 2026-10-14 |

## 1. Overview

Rare codon detection identifies codons whose reference frequency in a target organism falls below a chosen threshold.[1][2][3] In this repository, `CodonOptimizer.FindRareCodons` streams the positions, codons, amino-acid translations, and table frequencies for codons whose frequency is strictly less than the threshold (per-codon detection, default behaviour; an overload with `inclusiveThreshold: true` uses ≤). Two opt-in companions detect rare-codon **clusters / runs**: `CodonOptimizer.CalculateMinMaxProfile` computes the Clarke & Clark (2008) %MinMax sliding-window profile,[6] and `CodonOptimizer.FindRareCodonClusters` reports rare-codon clusters using the Sherlocc rule of Chartier et al. (2012) — a 7-codon window containing at least 4 rare ("pause") codons.[8] The methods are designed for sequence inspection and optimization workflows, especially when expressing genes in a heterologous host. They are deterministic, linear in the number of complete codons, and rely entirely on the supplied `CodonUsageTable`.[5]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Rare codons can slow translation because they are associated with lower-abundance tRNAs, and they are often discussed in the context of expression optimization and co-translational folding.[2][3][4] The original document recorded these common rare codons for E. coli K12 from Kazusa-backed frequencies:

| Codon | Amino Acid | Frequency | tRNA Gene |
|-------|------------|-----------|-----------|
| `AGA` | Arginine | `0.04` | `argU` |
| `AGG` | Arginine | `0.02` | `argW` |
| `CGA` | Arginine | `0.06` | `argW` |
| `CUA` | Leucine | `0.04` | `leuV` |

### 2.2 Core Model

For each codon in the normalized coding sequence, the repository looks up its frequency in the supplied codon-usage table and emits a result when:

```text
frequency(codon) < threshold
```

Each reported item is a tuple of the nucleotide position, codon, translated amino acid, and frequency.

**Cluster / run detection (opt-in).** Two published sliding-window methods extend the per-codon view to consecutive rare-codon runs:

- **%MinMax (Clarke & Clark 2008).**[6][7] For amino acid `i` with `n` synonymous codons, let `Xij` be the usage frequency of the codon actually used, `Xmax,i` / `Xmin,i` the most / least common synonymous codon frequencies, and `Xavg,i = (1/n) Σ Xij` the per-family mean. Over a window of `w` codons:

  ```text
  if Σ Xij > Σ Xavg,i :  %Max = Σ(Xij − Xavg,i) / Σ(Xmax,i − Xavg,i) × 100   (positive)
  if Σ Xij < Σ Xavg,i :  %Min = Σ(Xavg,i − Xij) / Σ(Xavg,i − Xmin,i) × 100   (negative)
  ```

  A window of predominantly rare codons appears as a negative %Min value; −100 is encoded with only the rarest synonymous codons, +100 with only the most common.[6] The default window is 18 codons.[6] The formula is reproduced term-for-term from the Clark-lab reference code (`calculateMinMax` in CHARMING.py[9]), whose synonymous families are those of the Standard code including the stop family {UAA, UAG, UGA}; families here come from `geneticCode` (default `GeneticCode.Standard`, the reference behaviour; under NCBI table 2 the Arg family loses AGA/AGG — they are stops — while Met becomes {AUA, AUG} and Trp {UGA, UGG}, so the profile changes as the code dictates). The reference takes codon usage **per thousand codons** (Kazusa "/1000" column), so window sums weight each residue by its overall usage; `CalculateMinMaxProfile` therefore uses `CodonUsageTable.PerThousand` when the table carries it — the `Yeast` (Kazusa 4932) and `Human` (Kazusa 9606) presets do (review 2026-10), and `EColiK12` (Kazusa 316407) since 2026-10-14 (SRC RO1) — and otherwise `CodonFrequencies`; a table whose `CodonFrequencies` hold per-thousand values also reproduces it exactly, while per-amino-acid relative fractions (a table without `PerThousand`, or `MinMaxUsageBasis.CodonFrequencies`) weight residues equally (identical only for single-residue windows, and then only when the fractions are exact ratios of the per-thousand values).

- **Sherlocc rare-codon cluster (RCC) rule (Chartier et al. 2012).**[8] A "seven position-wide window … containing at least four pause positions out of seven" is a rare-codon cluster, where a "pause"/"slow" position has a codon usage frequency "equal or lower than the specified threshold" (Sherlocc README; `sherlocc.pl`: `if($avrg <= $avrg_lim){ #We have a pause`).[8][10] `FindRareCodonClusters` therefore marks a codon as a pause when its table frequency is ≤ the threshold (since review 2026-10; `inclusiveThreshold: false` restores the strict `<`). Defaults: window 7 codons, ≥ 4 rare codons, threshold 0.15 (the cutoff of `FindRareCodons`).

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The provided codon-usage table reflects the organism or host of interest. | Reported rare codons may not be biologically rare in the actual target system. |
| ASM-02 | The selected threshold matches the user's screening goal. | Results can be overly conservative or overly broad depending on the cutoff. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported position is a multiple of `3`. | Positions are emitted as `codonIndex * 3`. |
| INV-02 | Every frequency reported by the default `FindRareCodons` is strictly less than `threshold` (≤ with `inclusiveThreshold: true`). | The default comparison is `<`. |
| INV-03 | Every reported codon has length `3`. | Only complete triplets are split into codons. |
| INV-04 | Re-running the method with the same inputs yields the same output. | The method is a deterministic single pass over normalized codons. |
| INV-05 | Every `CalculateMinMaxProfile` value lies in `[-100, 100]`. | The %MinMax numerator never exceeds the denominator (actual deviation ≤ max/min deviation).[6] |
| INV-06 | `CalculateMinMaxProfile` produces `codonCount − w + 1` windows when `codonCount ≥ w`, else none. | The window slides one codon at a time.[6] |
| INV-07 | A single-codon amino acid (Met/Trp) or an ambiguous triplet contributes `0` to a %MinMax window (no NaN). | `Xmax = Xmin = Xavg = Xij`, so its numerator and denominator terms are both `0`.[6] |
| INV-08 | Every `FindRareCodonClusters` cluster contains `≥ minRareCodons` rare codons. | A cluster originates from a qualifying window; merged regions only add codons.[8] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `codingSequence` | `string` | required | DNA or RNA coding sequence to inspect. | Empty or `null` input yields no results. |
| `table` | `CodonUsageTable` | required | Reference codon-usage table used for frequency lookups. | A valid codon absent from the table has frequency `0`. |
| `threshold` | `double` | `0.15` | Frequency cutoff for reporting a codon as rare. | Default comparison is strict (`< threshold`); the 4-argument overload's `inclusiveThreshold: true` reports `≤ threshold` (Sherlocc). |
| `geneticCode` | `GeneticCode?` | `null` | NCBI genetic code used to translate each reported codon (`FindRareCodons`) and to form the synonymous families (`CalculateMinMaxProfile`); `null` = `GeneticCode.Standard`, identical to the historical output (review 2026-09, F48). | Any NCBI table of `GeneticCode.SupportedTableNumbers`. Rare-codon *selection* is code-independent (a per-codon table lookup), so only the reported amino acid changes: `UGA` is `*` under table 1 and `W` under table 2. `FindRareCodonClusters` / `FindConservedRareCodonClusters` are fully code-independent and take no such parameter. |

### 3.2 Output / Return Value

| Name | Type | Description |
|------|------|-------------|
| `Position` | `int` | Nucleotide index of the codon start. |
| `Codon` | `string` | Normalized RNA codon. |
| `AminoAcid` | `string` | Single-letter translation under `geneticCode` (default NCBI table 1); `*` for a stop codon of that code. |
| `Frequency` | `double` | Frequency retrieved from the supplied codon-usage table. |

### 3.3 Preconditions and Validation

The sequence is read in frame 0, case-insensitively, with DNA `T` normalized to RNA `U`. Codons are extracted in complete triplets only; trailing bases are ignored. Only the 64 unambiguous codons over {A,C,G,U} are screened: a triplet containing an IUPAC ambiguity code (N, R, Y, …) or any other symbol has no codon-usage frequency and is skipped without shifting the frame (the codon set of `CodonUsageAnalyzer.CountCodons`, per EMBOSS `ajCodSetTripletsS`: "Skips triplets with ambiguity codes and any incomplete triplet at the end"). A valid codon absent from the table (never observed in the reference genes) has frequency `0` and is reported whenever `threshold > 0`.[5]

## 4. Algorithm

### 4.1 High-Level Steps

1. Return an empty sequence for `null` or empty input.
2. Normalize the sequence to uppercase RNA notation.
3. Split the sequence into complete frame-0 codons (shared core `CodonUsageAnalyzer.SplitInFrameCodons`); skip ambiguous triplets without shifting the frame.
4. For each codon, read its table frequency.
5. If the frequency is strictly less than `threshold`, translate the codon with `geneticCode` (default `GeneticCode.Standard`) and yield `(position, codon, aminoAcid, frequency)`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The original document recorded the following threshold-selection guidance:

| Threshold | Description | Use Case |
|-----------|-------------|----------|
| `0.10` | Very rare | Critical optimization |
| `0.15` | Default | Standard analysis |
| `0.20` | Moderately rare | Broad screening |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindRareCodons` | `O(n)` | `O(1)` auxiliary | Streams results while scanning complete codons. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CodonOptimizer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs)

- `CodonOptimizer.FindRareCodons(string, CodonUsageTable, double)`: per-codon detection (default, strict `<`); `FindRareCodons(string, CodonUsageTable, double, bool inclusiveThreshold)` selects `≤`.
- `CodonOptimizer.CalculateMinMaxProfile(string, CodonUsageTable, int)`: %MinMax sliding-window profile (Clarke & Clark 2008).[6] Returns `IReadOnlyList<MinMaxWindow>` (`WindowStartCodon`, signed `PercentMinMax`). Uses `CodonUsageTable.PerThousand` when present; `CalculateMinMaxProfile(string, CodonUsageTable, int, MinMaxUsageBasis)` forces `CodonFrequencies`.
- `CodonUsageTable.PerThousand` (optional, `init`): usage per thousand codons of the same table entry; set on `Yeast` (Kazusa 4932, one-decimal column from CHARMING `ScerCUB.txt`) and `Human` (Kazusa 9606, 1000·count/40 662 582 from the Kazusa counts shipped in PyPI codon-optimizer 0.2.0); `EColiK12` (Kazusa 316407, published one-decimal column fetched from kazusa.or.jp 2026-10-14, SRC RO1; `null` before).
- `CodonOptimizer.FindRareCodonClusters(string, CodonUsageTable, double, int, int)`: Sherlocc rare-codon clusters (Chartier et al. 2012),[8] pause = frequency ≤ threshold; the 6-argument overload with `inclusiveThreshold: false` uses `<`. Returns `IReadOnlyList<RareCodonCluster>` (`StartCodon`, `EndCodon`, `RareCount`), maximal non-overlapping regions in codon-index order.
- `CodonOptimizer.FindConservedRareCodonClusters(IReadOnlyList<string>, CodonUsageTable | IReadOnlyList<CodonUsageTable>, double threshold, SherloccOptions?)`: the Sherlocc 1.0 stage-3 analysis of a codon alignment of homologous CDS (gaps `-`/`.`), ported from `sherlocc.pl`:[10] per-cell per-thousand usage (one table per sequence = Sherlocc's per-organism usage), centred 7-column average over ALL sequences, representativeness/conservation filter (≥ 0.75 × sequences × 7 valued cells), pause = average ≤ threshold, cluster = first…last pause of each 7-column window with ≥ 4 pauses, overlap-only merging. Returns `SherloccResult` (column averages / usage counts, pause columns, `ConservedRareCodonCluster(StartColumn, EndColumn, PauseCount, MeanUsage, MiddlePoint)`, total length, fraction, mean cluster usage). Options `SherloccOptions` (AveragingWindow 7, ClusterWindow 7, MinPauses 4, MinCoverage 0.75, UseCodonFrequencies false).

### 5.2 Current Behavior

The method converts `T` to `U`, ignores trailing incomplete codons, skips ambiguous triplets (frame preserved), and reports nucleotide positions as `i * 3`. A valid codon absent from the table is treated as frequency `0`. Stop codons are screened like any other table row (E. coli `UAG` 0.07 is rare at 0.10). The default threshold comparison of `FindRareCodons` is strict `<`, which means a codon exactly at the threshold is not reported;[5] `FindRareCodonClusters` uses Sherlocc's `≤` by default.[10]

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Codons are flagged as rare when their reference frequency is below the chosen threshold.[2][5]
- The result includes position, codon, amino acid, and frequency, matching the repository test specification.[5]
- The %MinMax window formula (per-amino-acid `Xij`, `Xmax,i`, `Xmin,i`, `Xavg,i`; `%Max`/`%Min` selection by summed-window comparison; 18-codon default) is reproduced exactly.[6][7]
- The Sherlocc rare-codon-cluster rule (7-codon window, ≥ 4 pause positions; pause = value ≤ threshold, `sherlocc.pl`) is reproduced, with tunable window/threshold parameters.[8][10]
- The Sherlocc 1.0 alignment analysis — family/window averaging of per-thousand usage, the 75 % representativeness (conservation/gap) filter, cluster spans, overlap merging, cluster mean usage, middle point and fraction — is reproduced by `FindConservedRareCodonClusters` and locked against `sherlocc.pl` output on three alignments (e.g. yeast alignment A, threshold 15: cluster 9–12, mean 9.83928571428571).[10]

**Intentionally simplified:**

- Per-codon detection is table-threshold-based and does not model ribosome dynamics or local codon context; **consequence:** the per-codon output is a screening list rather than a direct translation-efficiency prediction. (Cluster context is now available via the opt-in methods.)
- The %MinMax reference implementation takes codon usage per thousand codons. The `EColiK12`, `Yeast` and `Human` presets carry the per-thousand usage of their Kazusa entries (316407 / 4932 / 9606) and reproduce the reference by default (yeast CUG·AGA, window 2: 33.99 = CHARMING; E. coli: 58.10 = CHARMING 58.1); **consequence:** with `MinMaxUsageBasis.CodonFrequencies` (or a caller table of fractions without `PerThousand`) mixed-residue windows weight every amino acid equally and differ from the reference (yeast fractions: 58.62 vs 33.99; E. coli fractions: 36.47 vs 58.10). Supplying a table with per-thousand values reproduces the reference exactly.
- `FindRareCodonClusters` (single sequence) compares each codon's own table value (default per-amino-acid fraction, threshold 0.15) instead of Sherlocc's window-and-family average; **consequence:** its pause positions are per-codon, single-sequence. The full Sherlocc analysis (family averaging, 75 % filter, per-thousand threshold) is `FindConservedRareCodonClusters`, which reproduces `sherlocc.pl` (verified by running the script, Evidence "Review 2026-10 (FIN-B02 WP8)").[10]
- `FindConservedRareCodonClusters` keeps two script behaviours verbatim: a non-representative column inside a cluster span contributes 0 to the cluster's mean usage (Perl undefined value), and abutting clusters are not merged; the FoldIndex/Uversky unstructured-region output and the stage-1/2 UniProt/ENA retrieval of the script are not ported (the caller supplies the codon alignment and tables).
- Sherlocc (and `FindConservedRareCodonClusters`) reports a cluster from the first to the last pause inside each qualifying window and merges only overlapping clusters; `FindRareCodonClusters` reports the full window span and also merges adjacent windows; **consequence:** region boundaries can differ from Sherlocc's by up to window − 1 codons.[10]

**Not implemented:**

- Per-cluster P-values of Chartier et al. (2012): `sherlocc.pl` contains no P-value / randomisation code and the paper's statistical model could not be obtained (full text unreachable; Evidence "Review 2026-10 (FIN-B02 WP8)"), so no P-value is reported;[8] **users should rely on:** the published Sherlocc results for significance.
- Automatic integration with codon optimization inside these methods; **users should rely on:** [Sequence_Optimization.md](Sequence_Optimization.md).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns no results. | Explicit early exit. |
| Default-threshold call | Uses `0.15`. | Default parameter value in the public API. |
| Frequency exactly equal to threshold | `FindRareCodons`: not reported by default (`<`), reported with `inclusiveThreshold: true`; `FindRareCodonClusters`: a pause by default (Sherlocc `≤`). | Default contracts; Sherlocc `$avrg <= $avrg_lim`. |
| Ambiguous triplet (`NNN`, `RYU`, …) | Skipped (not reported, not a cluster pause position, no %MinMax contribution); later positions keep their frame. | No codon-usage frequency exists for it (EMBOSS `ajCodSetTripletsS`). |
| Valid codon absent from the table | Reported with frequency `0` when `threshold > 0`. | Never observed in the reference genes. |
| Stop codon | Screened like any table row; amino acid `*`. | Stops are rows of the Kazusa/cusp table. |
| Alignment column with too many gaps (`FindConservedRareCodonClusters`) | Column average `null` ("NA"), never a pause; inside a cluster span it contributes 0 to the mean. | `sherlocc.pl` l. 772/799 (`$freq_count >= $tresh`), l. 875–879. |
| Alignment rows of unequal length / length not a multiple of 3 / table without `PerThousand` | `ArgumentException` (`SherloccOptions.UseCodonFrequencies = true` scores fractions instead). | Codon alignment contract; Sherlocc uses CUTG per-thousand. |
| Incomplete trailing codon | Ignored. | Codon splitting only emits complete triplets. |

### 6.2 Limitations

Per-codon `FindRareCodons` is a simple table lookup over codon triplets and does not consider local sequence context or initiation effects. Rare-codon **clustering / runs** are now covered by the opt-in `CalculateMinMaxProfile` (%MinMax) and `FindRareCodonClusters` (Sherlocc) methods, but these still operate purely on the supplied frequency table: they do not model ribosome dynamics, do not weight 5'-vs-internal position, and no cluster P-value is computed (see §5.3 Not implemented). Family-level (alignment) conservation is available through `FindConservedRareCodonClusters`.

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
string sequence = "AUGAGAAGGCGA";
var rareList = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10);
// Per-codon: (3, "AGA", "R", 0.04), (6, "AGG", "R", 0.02), (9, "CGA", "R", 0.06)

// Cluster / run detection (opt-in):
string run = string.Concat(Enumerable.Repeat("AGA", 7)); // 7 rare Arg codons
var clusters = CodonOptimizer.FindRareCodonClusters(run, CodonOptimizer.EColiK12);
// Sherlocc 7/4 rule -> one cluster: RareCodonCluster(StartCodon: 0, EndCodon: 6, RareCount: 7)

var profile = CodonOptimizer.CalculateMinMaxProfile("AGAAGAAGA", CodonOptimizer.EColiK12, windowSize: 3);
// %MinMax: one window, PercentMinMax ≈ −88.9117 (Kazusa 316407 per thousand; a rare-codon %Min trough;
// MinMaxUsageBasis.CodonFrequencies gives −86.3636 from the fractions)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests (per-codon): [CodonOptimizer_FindRareCodons_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_FindRareCodons_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Tests (clusters / runs): [CodonOptimizer_RareCodonClusters_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_RareCodonClusters_Tests.cs) — covers `INV-05`, `INV-06`, `INV-07`, `INV-08`
- Tests (Sherlocc conserved clusters over an alignment): [CodonOptimizer_SherloccConservedClusters_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_SherloccConservedClusters_Tests.cs)
- Tests (Sherlocc threshold, per-thousand presets): [CodonOptimizer_SherloccThreshold_PerThousandPresets_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_SherloccThreshold_PerThousandPresets_Tests.cs)
- Test specification: [CODON-RARE-001.md](../../../tests/TestSpecs/CODON-RARE-001.md)
- Related algorithms: [Sequence_Optimization.md](Sequence_Optimization.md), [CAI_Calculation.md](CAI_Calculation.md)

## 8. References

1. Sharp PM, Li WH. 1987. The codon adaptation index-a measure of directional synonymous codon usage bias. Nucleic Acids Research. N/A
2. Shu P, Dai H, Gao W, Goldman E. 2006. Inhibition of translation by consecutive rare leucine codons in E. coli. Gene Expression. N/A
3. Plotkin JB, Kudla G. 2011. Synonymous but not the same: causes and consequences of codon bias. Nature Reviews Genetics. N/A
4. Kane JF. 1995. Effects of rare codon clusters on high-level expression of heterologous proteins in E. coli. Current Opinion in Biotechnology. N/A
5. Test specification: [CODON-RARE-001.md](../../../tests/TestSpecs/CODON-RARE-001.md)
6. Clarke TF, Clark PL. 2008. Rare Codons Cluster. PLoS ONE 3(10):e3412. https://doi.org/10.1371/journal.pone.0003412
7. Rodriguez A, Wright G, Emrich S, Clark PL. 2018. %MinMax: A versatile tool for calculating and comparing synonymous codon usage and its impact on protein folding. Protein Science. https://pmc.ncbi.nlm.nih.gov/articles/PMC5734269/
8. Chartier M, Gaudreault F, Najmanovich R. 2012. Large-scale analysis of conserved rare codon clusters suggests an involvement in co-translational molecular recognition events. Bioinformatics 28(11):1438–1445. https://doi.org/10.1093/bioinformatics/bts149
9. Wright G, Rodriguez A, Li J, Milenković T, Emrich SJ, Clark PL. 2022. CHARMING: Harmonizing synonymous codon usage to replicate a desired codon usage pattern. Protein Science 31(1), doi:10.1002/pro.4223 (PMID 34738275). Reference code https://github.com/wrightgs/CHARMING (`CHARMING.py` `calculateMinMax`, `ScerCUB.txt`).
10. Chartier M, Najmanovich R. SHERLOCC 1.0 reference implementation, `sherlocc.pl` and README. https://github.com/mtthchrtr/sherlocc (commit b3b4d15987c5f8205d1957db3a757dbd390d5049, 2020-06-25; opened 2026-10-14).
