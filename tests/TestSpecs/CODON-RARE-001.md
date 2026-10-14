# Test Specification: CODON-RARE-001

## Test Unit Information
- **ID**: CODON-RARE-001
- **Title**: Rare Codon Detection
- **Area**: Codon Optimization
- **Complexity**: O(n)
- **Status**: ☑ Complete

## Method Under Test

```csharp
public static IEnumerable<(int Position, string Codon, string AminoAcid, double Frequency)> 
    FindRareCodons(string codingSequence, CodonUsageTable table, double threshold = 0.15)
```

## Evidence Summary

Sources: Wikipedia (Codon usage bias), Kazusa Codon Usage Database, GenScript GenRCA, Shu et al. (2006), Sharp & Li (1987)

Key findings:
- Rare codons (AGA, AGG, CGA, CUA) have low tRNA abundance in E. coli
- Default threshold 0.15 is reasonable for most applications
- Position, codon, amino acid, and frequency are all valuable output fields

---

## Test Categories

### MUST Tests (Evidence-backed, Critical)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| M01 | Empty sequence returns empty | `FindRareCodons("", table)` → empty enumerable | Implementation contract |
| M02 | Single rare codon detection | Sequence with one AGA (0.04) detected below threshold 0.10 | Shu et al., Kazusa W3110 (316407) |
| M03 | Multiple rare codons detection | Sequence with AGA, AGG, CGA all detected | Kazusa frequencies |
| M04 | Position is nucleotide index | AGA at codon index 1 reports position 3 | Implementation spec |
| M05 | Threshold boundary - below | Codon at freq 0.07 detected with threshold 0.10 | Math invariant |
| M06 | Threshold boundary - at | Codon at freq 0.10 NOT detected with threshold 0.10 | Math invariant (< not ≤) |
| M07 | No rare codons found | Sequence with only common codons returns empty | Logic correctness |
| M08 | Amino acid translation correct | Reported AA matches genetic code | Sharp & Li |
| M09 | Frequency value matches table | Reported frequency equals table lookup | Implementation spec |
| M10 | DNA T converted to RNA U | Input with T handled correctly | Implementation spec |

### SHOULD Tests (Important behavior)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| S01 | Default threshold is 0.15 | Single-arg overload uses 0.15 | Implementation default |
| S02 | Incomplete codons ignored | 10-nucleotide sequence processes only 3 codons | Implementation spec |
| S03 | Case insensitivity | Lowercase input handled | Robustness |
| S04 | All codons rare scenario | Sequence of all rare codons returns all positions | Stress case |
| S05 | Mixed rare/common codons | CUG (common) + CUA (rare) - only CUA reported | Kazusa frequencies |

### COULD Tests (Edge cases, nice-to-have)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| C01 | Threshold zero | All codons with freq > 0 NOT reported as rare | Math edge |
| C02 | Threshold one | All codons reported as rare | Math edge |
| C03 | Different organism tables | Yeast table has different rare codons than E. coli | Kazusa data |
| C04 | Ambiguous / unknown triplet handling | Ambiguous triplet (NNN, RYU) skipped, frame preserved (review 2026-09; was: reported with freq 0 / `X`); valid codon absent from table → freq 0 | EMBOSS ajCodSetTripletsS; CodonUsageAnalyzer.CountCodons |

---

## Test Data

### E. coli K12 Rare Codons (freq < 0.10) — Kazusa K-12 W3110 (species=316407)
- AGA: 0.04 (Arginine)
- AGG: 0.02 (Arginine)  
- CGA: 0.06 (Arginine)
- CUA: 0.04 (Leucine)
- UAG: 0.07 (Stop)

### E. coli K12 Common Codons (freq > 0.40)
- CUG: 0.50 (Leucine)
- AUG: 1.00 (Methionine)
- UGG: 1.00 (Tryptophan)

### Test Sequences

| Name | Sequence | Description |
|------|----------|-------------|
| SingleRare | AUGAGA | M + R(rare) |
| MultiRare | AUGAGAAGGCGA | M + R(rare) + R(rare) + R(rare) |
| NoRare | AUGCUGUGG | M + L(common) + W(common) |
| AllRare | AGAAGGCGA | All rare arginine codons |
| MixedLeu | CUGCUA | L(common 0.50) + L(rare 0.04) |
| Borderline | CUGAUA | L(common) + I(0.07 near threshold) |

---

## Invariants

1. All reported positions are multiples of 3
2. All reported positions are in range [0, sequence.Length - 3]
3. All reported frequencies are < threshold
4. All reported codons have length 3
5. Same input → same output (deterministic)
6. Protein sequence is preserved (amino acid translation correct)

---

## Audit Notes

### Existing Tests (CodonOptimizerTests.cs)
- `FindRareCodons_IdentifiesRare` - basic detection ✓
- `FindRareCodons_ReturnsPosition` - position check ✓
- `FindRareCodons_EmptySequence_ReturnsEmpty` - empty input ✓

### Coverage Gaps (All Closed)
- ~~No threshold boundary tests~~ ✅ Covered
- ~~No threshold=0 or threshold=1 edge cases~~ ✅ Covered
- ~~No multi-rare codon sequence test~~ ✅ Covered
- ~~No case insensitivity test~~ ✅ Covered
- ~~No DNA (T) conversion test~~ ✅ Covered
- ~~No frequency value verification~~ ✅ Covered
- ~~No amino acid translation verification~~ ✅ Covered
- ~~No incomplete codon test~~ ✅ Covered

### Consolidation Plan
- ~~Create dedicated test file: `CodonOptimizer_FindRareCodons_Tests.cs`~~ ✅ Done
- ~~Move existing 3 tests to new file with consistent naming~~ ✅ Done
- ~~Add missing MUST tests~~ ✅ Done
- ~~Add selected SHOULD/COULD tests~~ ✅ Done

---

## Open Questions / Decisions

1. **Q**: Should threshold comparison be `<` or `<=`?  
   **A**: `FindRareCodons` uses `<` (strict less than) by default; since 2026-10 an overload with
   `inclusiveThreshold: true` uses `<=`. `FindRareCodonClusters` follows its source, Sherlocc
   (`sherlocc.pl` `if($avrg <= $avrg_lim)`), and uses `<=` by default (`inclusiveThreshold: false` → `<`).

2. **Q**: What happens with non-standard codons?  
   **A** (revised 2026-09): an ambiguous / non-ACGU triplet has no usage frequency and is skipped without shifting the frame; a valid codon absent from the table has frequency 0 and is flagged whenever threshold > 0.

---

# Addendum (2026-06-24): Rare-Codon Cluster / Run Detection

Closes the documented limitation that per-codon detection does not find consecutive
rare-codon **clusters / runs / ramps**. Two opt-in methods added; per-codon `FindRareCodons`
unchanged. Sources retrieved this session (see `docs/Evidence/CODON-RARE-001-Evidence.md`
Addendum): Clarke & Clark (2008) %MinMax; Chartier et al. (2012) Sherlocc RCC.

## Methods Under Test (clusters)

```csharp
public static IReadOnlyList<MinMaxWindow> CalculateMinMaxProfile(
    string codingSequence, CodonUsageTable table, int windowSize = 18);

public static IReadOnlyList<RareCodonCluster> FindRareCodonClusters(
    string codingSequence, CodonUsageTable table,
    double rareThreshold = 0.15, int windowSize = 7, int minRareCodons = 4);
```

| Method | Type | Source / parameters |
|--------|------|---------------------|
| `CalculateMinMaxProfile` | Canonical | Clarke & Clark (2008), PLoS ONE 3:e3412; window 18 codons |
| `FindRareCodonClusters` | Canonical | Chartier et al. (2012), Bioinformatics 28(11):1438; 7-codon window, ≥4 rare |

## Evidence-derived expected values (E. coli K12 Arg family; Xavg=0.16667, Xmax=0.40, Xmin=0.02)

| Case | Input | Method / params | Expected |
|------|-------|-----------------|----------|
| MM1 | AGA×3 | %MinMax, w=3 | −86.363636…% |
| MM2 | CGC×3 | %MinMax, w=3 | +100% |
| MM3 | CUG·AGA | %MinMax, w=2 | +36.470588…% |
| C1 | AGA×7 | RCC default | 1 cluster, codons 0–6, 7 rare |
| C2 | AGA×3+CGC×4 | RCC default | no cluster (3 < 4) |
| C3 | AGA×4+CGC×3 | RCC default | 1 cluster, 4 rare (boundary ≥4) |
| C8 | AGA×7 (yeast) | RCC default | no cluster (AGA 0.48 common) |

## Cluster Test Cases (file `CodonOptimizer_RareCodonClusters_Tests.cs`)

%MinMax: MM1 all-rare %Min (exact), MM2 all-common +100, MM3 mixed %Max (exact), MM4 sliding
windows count+order, MM5 short-sequence empty, MM6 empty/null, MM7 single-codon-AA 0 (no NaN),
MM8 bound `[-100,100]`, MM9 DNA→RNA, MM10 windowSize<1 throws.

RCC: C1 7-rare cluster, C2 below-threshold (3<4) no cluster, C3 exactly-4 boundary, C4 isolated
rare codons → per-codon flags 3 but zero clusters (the closed gap), C5 long run merges to one
cluster, C6 two separated runs → two clusters, C7 default-threshold = 0.15, C8 organism-specific,
C9 short-sequence empty, C10 empty/null, C11 DNA→RNA, C12 tunable window/threshold short run,
C13 invalid params throw, C14 determinism.

## Invariants (clusters)

- INV-05: every %MinMax value ∈ `[-100, 100]`.
- INV-06: `CalculateMinMaxProfile` yields `codonCount − w + 1` windows (or none).
- INV-07: single-codon AA contributes 0 to %MinMax (no NaN).
- INV-08: every cluster contains ≥ `minRareCodons` rare codons.

## Phase 7 Work Queue (clusters) — all ✅ Done

24 cluster cases (MM1–MM10, C1–C14) implemented and executed: **Failed: 0, Passed: 24**.
Remaining: 0.

## Status note

Per the validation protocol, the ROOT `ALGORITHMS_CHECKLIST_V2.md` Status for CODON-RARE-001
is reset to ☐ pending independent re-validation of this new capability.

---

# Review 2026-09 (campaign B02) — contract changes and new locks

Sources opened: python-codon-tables 0.1.18 (PyPI wheel; Kazusa-derived `e_coli_316407`,
`s_cerevisiae_4932`, `h_sapiens_9606` relative frequencies — all 192 preset values identical);
Clark-lab CHARMING repository (raw.githubusercontent.com/wrightgs/CHARMING: `CHARMING.py`
`calculateMinMax`, `README.md`, `ScerCUB.txt`, `EcolCUB.txt`); Kane (1995) rare-codon list via
WebSearch snippet (AGG, AGA, CUA, AUA, CGA, CCC).

| ID | Case | Expected (reference) | Test |
|----|------|----------------------|------|
| R1 | `AUGNNNAGAcgaRYU`, E. coli, 0.15 | (6,AGA,R,0.04), (9,CGA,R,0.06); NNN/RYU skipped | `FindRareCodons_AmbiguousTriplets_SkippedFramePreserved` |
| R2 | Kane 1995 codons `AGGAGACUAAUACGACCC`, default | all six flagged (0.02, 0.04, 0.04, 0.07, 0.06, 0.12) | `FindRareCodons_Kane1995EColiRareCodons_AllFlaggedAtDefault` |
| R3 | `UAGUAA`, E. coli, 0.10 | (0,UAG,*,0.07) | `FindRareCodons_RareStopCodon_ReportedAsStop` |
| R4 | %MinMax, yeast per-thousand, `CUGAGA`, w=2 | 33.99209486166008 (CHARMING: 33.99) | `CalculateMinMaxProfile_PerThousandTable_MatchesClarkLabReference` |
| R5 | %MinMax, yeast per-thousand, 20-codon gene, w=18 | 12.6879134095, 17.0285714286, 15.9910628790 (CHARMING: 12.69, 17.03, 15.99) | same |
| R6 | %MinMax, yeast relative preset, `CUGAGA`, w=2 | 58.620689655172406 (equal residue weighting) | `CalculateMinMaxProfile_RelativeFractionPreset_WeightsResiduesEqually` |
| R7 | %MinMax, E. coli, `UAGUAG`, w=2 | −100 (stop family is synonymous) | `CalculateMinMaxProfile_StopFamily_IsSynonymousFamily` |
| R8 | Clusters, `NNN`×7 / `AGA`×4+`NNN`×3 | none / (0,6,4) | `FindRareCodonClusters_AmbiguousTriplets_AreNotPausePositions` |

Randomised cross-check (400 sequences incl. lower case, RNA/DNA, IUPAC N/R/Y, partial codons;
4 tables incl. per-thousand yeast; thresholds 0–1; windows 1–12): 830 rare hits, 3989 %MinMax
windows, 101 clusters — all identical to the Python reference (1e-9).

---

# Review 2026-10 (FIN-B02 WP7) — Sherlocc threshold (A1-1), per-thousand presets (A1-3)

Sources opened: Sherlocc 1.0 (`git clone github.com/mtthchrtr/sherlocc`, commit b3b4d15; `sherlocc.pl`
l. 803 `if($avrg <= $avrg_lim){ #We have a pause`, README "equal or lower than the specified
threshold"); Kazusa 4932 per thousand (CHARMING `ScerCUB.txt`); Kazusa 9606 counts (PyPI
codon-optimizer 0.2.0 `Homo_sapiens_codon_frequency.csv`, Σ = 40 662 582); CHARMING
`calculateMinMax` + Python port as %MinMax oracle. Kazusa 316407 per thousand: not obtainable.

| ID | Case | Expected (reference) | Test |
|----|------|----------------------|------|
| R9 | Clusters, E. coli, `UCU UCC UCG AGU GUA GGG CGC` (six at 0.15) | default (≤, Sherlocc): (0,6,6); `inclusiveThreshold:false`: none | `FindRareCodonClusters_DefaultIsSherloccInclusive_…`, `…_StrictOption_ReproducesPreviousBehaviour` |
| R10 | `FindRareCodons("UCUAGA")`, E. coli, 0.15 | default: (3,AGA,R,0.04); inclusive: + (0,UCU,S,0.15) | `FindRareCodons_DefaultStaysStrict_InclusiveOptionAddsCodonsAtThreshold` |
| R11 | Preset per-thousand columns | EColiK12 null (superseded 2026-10-14 by SRC RO1, R17); Yeast 64 values (UUU 26.1, Σ 1000.5); Human UUU 17.56646934028931 (Σ 1000) | `Presets_PerThousand_SourcedEntries` |
| R12 | Per-thousand ratios vs stored fractions | Human all 64 within 0.005; Yeast non-stop within 0.005, stops within 0.013 | `Presets_PerThousandRatios_MatchStoredFractions` |
| R13 | %MinMax Yeast preset (default) `CUGAGA` w2 / 20-codon gene w18 | 33.99209486166008 / 12.68791340950089, 17.028571428571425, 15.991062879029688 (CHARMING 33.99 / 12.69, 17.03, 15.99) | `CalculateMinMaxProfile_YeastPreset_DefaultMatchesClarkLabReference` |
| R14 | %MinMax Human preset gene w18 / `UUACGU` w2 | 28.796581754809075, 23.264123403923747, 16.195449171123403 / −96.40404591634571 (CHARMING 28.8, 23.26, 16.2 / −96.4); `CodonFrequencies` basis 19.650655021834133, 13.878080415045439, 13.803548247511891 / −94.54545454545456 | `CalculateMinMaxProfile_HumanPreset_…` |
| R15 | %MinMax E. coli preset `CUGAGA` w2 | 36.470588235294116 under both bases (no per-thousand) — superseded 2026-10-14 by SRC RO1 (R19): default 58.10113519091846, `CodonFrequencies` 36.470588235294116 | `CalculateMinMaxProfile_EColiPreset_DefaultPerThousand_FractionBasisUnchanged` |
| R16 | Custom table with `PerThousand` init; undefined basis | 33.99209486166008; `ArgumentOutOfRangeException` | `…_CustomTableWithPerThousand_UsesIt`, `…_UndefinedBasis_Throws` |

R6 (yeast fraction value 58.620689655172406) is now requested with `MinMaxUsageBasis.CodonFrequencies`
(value unchanged). MCP: `find_rare_codons` `inclusive_threshold` (default false),
`find_rare_codon_clusters` `inclusive_threshold` (default true), `codon_minmax_profile` `usage_basis`
(default `PerThousandWhenAvailable`) — `SherloccThresholdPerThousandToolsTests.cs`.

---

# Review 2026-10 (FIN-B02 WP8) — Sherlocc conserved rare-codon clusters over an alignment (A1-2)

```csharp
public static SherloccResult FindConservedRareCodonClusters(
    IReadOnlyList<string> alignedCodingSequences, CodonUsageTable table,
    double threshold, SherloccOptions? options = null);
public static SherloccResult FindConservedRareCodonClusters(
    IReadOnlyList<string> alignedCodingSequences, IReadOnlyList<CodonUsageTable> tablesPerSequence,
    double threshold, SherloccOptions? options = null);
```

Source / oracle: Sherlocc 1.0 `sherlocc.pl` (commit b3b4d15) stages 2+3 run with perl 5.38 on the
alignments below (stage-1 records hand-written in the l. 369–374 format; l. 15/17 stage-1 reset removed
because stage 1 queries UniProt/ENA live); `files/cut` rows = Yeast (4932) / Human (9606) preset
per-thousand values. Unrounded column averages: independent Python port (equal to every Perl value).
Alignments are listed in `CodonOptimizer_SherloccConservedClusters_Tests.cs`.

| ID | Case | Expected (sherlocc.pl) | Test |
|----|------|------------------------|------|
| S1 | A: yeast, 3×20 columns, threshold 15 — column averages / usage counts | NA,NA,38.71176470588236 … 4.79047619047619 (col 10) … 41.32352941176471,NA,NA; counts 11,14,17,20,…,21,…,11 | `AlignmentA_ColumnAverages_MatchSherlocc` |
| S2 | A: pauses / clusters | pauses 8–11 (0-based); `clusters=9-12\|9.83928571428571` (4 windows merged), middle point 10, fraction 0.2 | `AlignmentA_OneMergedCluster_MatchesSherloccResults` |
| S3 | B: yeast 2×18, columns 6 & 12 gapped → column 9 has 10 < 10.5 valued cells | col 9 NA; pauses 5–8,10–12; `clusters=6-13\|5.265625` (NA column counted as 0), middle 9 | `AlignmentB_NonRepresentativeColumnInsideCluster_CountsAsZero` |
| S4 | C: per-sequence tables yeast/human/yeast, 3×34, threshold 16 | 14 pauses; `clusters=4-10\|9.20962052101093,21-27\|9.30926559372343`; `clusters_freq_avrg=9.25944305736718`; middle 6, 23 | `AlignmentC_PerSequenceTables_TwoClusters_MatchSherlocc` |
| S5 | Threshold at / just below column-8 average of A | `≤` (l. 803): pause, 1 cluster; below: 3 pauses, no cluster | `Threshold_IsInclusive_PerSherlocc` |
| S6 | Single-table vs per-sequence overload; RNA/lower case; `.` gaps | identical results | `SingleTableOverload_EqualsPerSequenceOverload_AndAcceptsRnaLowerCase` |
| S7 | Defaults | AveragingWindow 7, ClusterWindow 7, MinPauses 4, MinCoverage 0.75 (l. 772, 785, 838–848) | `Defaults_AreSherloccConstants` |
| S8 | E. coli preset (no per-thousand; since SRC RO1 2026-10-14: `EColiK12 with { PerThousand = null }`) | `ArgumentException`; `UseCodonFrequencies` → col 3 avg 0.03, pauses 2–4, no cluster | `TableWithoutPerThousand_Rejected_UnlessCodonFrequenciesRequested` |
| S9 | Invalid input (null, empty, length %3, unequal rows, table count, NaN, even window, coverage, MinPauses) | exceptions | `InvalidInput_Throws` |
| S10 | 4-column alignment | no representative column, no cluster window | `ShortAlignment_NoClusterWindow_NoClusters` |

Cluster P-values: `sherlocc.pl` has no P-value code; the paper's model could not be retrieved
(BLOCKED, see Evidence) — not implemented, no test.

### Review 2026-09 (F48 — genetic-code-aware API, audit A2-1)

`FindRareCodons` (both overloads) and `CalculateMinMaxProfile` (both overloads) gained an additive
optional `GeneticCode? geneticCode = null` (`null` = `GeneticCode.Standard`, outputs unchanged;
DNA Chisel `EnforceTranslation(genetic_table=…)` / `get_backtranslation_table`).

| # | Scenario | Expectation | Test |
|---|----------|-------------|------|
| G1 | `FindRareCodons`, default vs explicit Standard | identical tuples | `FindRareCodons_GeneticCode_ChangesOnlyTheAminoAcid` |
| G2 | `FindRareCodons`, table 2 | same positions/codons/frequencies (rarity is a per-codon table lookup); `UGA` reported as `W`, `AGA` as `*` | same test |
| G3 | `CalculateMinMaxProfile`, default vs explicit Standard; 4- vs 5-argument overload | identical profiles | `CalculateMinMaxProfile_GeneticCode_UsesThatCodesFamilies` |
| G4 | `CalculateMinMaxProfile`, table 2 | same window count, different values (Arg loses AGA/AGG; Met {AUA, AUG} and Trp {UGA, UGG} become two-codon families) | same test |

`FindRareCodonClusters` and `FindConservedRareCodonClusters` are genetic-code independent (a pause
is a per-codon usage lookup; no codon is translated and no family is formed) and therefore take no
`geneticCode` parameter — documented in the XML docs and in `Rare_Codon_Detection.md`.

# SRC RO1 (2026-10-14) — Kazusa 316407 per-thousand column on `EColiK12`

Source: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=316407 ("Escherichia coli W3110
[gbbct]: 4332 CDS's (1372057 codons)", fetched 2026-10-14); oracle: CHARMING `calculateMinMax` run
isolated + independent Python port (Evidence "SRC RO1"). F39 design kept: default %MinMax basis =
per thousand when present; Sherlocc default requires per thousand.

| ID | Case | Expected | Test |
|---|---|---|---|
| R17 | `EColiK12.PerThousand` | 64 values = Kazusa 316407 /1000 column verbatim (UUU 22.2 … GGG 11.0), Σ 1000.0 | `EColiPerThousand_EqualsKazusa316407Verbatim`, `Presets_PerThousand_SourcedEntries` |
| R18 | Per-thousand ratios vs stored fractions | all 64 within 0.0055 (one-decimal rounding) | `EColiPerThousandRatios_MatchStoredFractions` |
| R19 | %MinMax E. coli `CUGAGA` w2 | default 58.10113519091846 (CHARMING 58.1; was 36.470588235294116); `CodonFrequencies` 36.470588235294116 (unchanged) | `CalculateMinMaxProfile_EColiPreset_DefaultPerThousand_FractionBasisUnchanged`; MCP `CodonMinMaxProfile_PresetsUsePerThousandByDefault` |
| R20 | %MinMax E. coli 20-codon gene w18; `UUACGU` w2 | 7.417300184918832, 7.371763071733311, −6.565186442534101 (CHARMING 7.42, 7.37, −6.57); fractions −1.358411703239252, −0.9027081243731403, −8.181517058628705; `UUACGU` 16.202270381836932 (16.2) | `CalculateMinMaxProfile_EColiPreset_Gene20_MatchesClarkLabReference` |
| R21 | MM1/MM3/MM4 default basis (`AGAAGAAGA` w3, `CUGAGA` w2, `CGCCGCAGAAGAAGA` w3) | −88.91170431211498; 58.10113519091846; 48.280254777070056, −5.544147843942496, −88.91170431211498 (CHARMING −88.91; 58.1; 48.28, −5.54, −88.91). Hand-derived fraction values (−86.3636…, 36.4706…) kept under explicit `CodonFrequencies` | `CodonOptimizer_RareCodonClusters_Tests` MM1/MM3/MM4; MCP `CodonMinMaxProfile_MatchesLockedApiValues` |
| R22 | Sherlocc accepts `EColiK12` by default | AGA×10/AGG×10, threshold 15: columns 2–7 average 1.55, pauses 2–7, one cluster 2–7 (6 pauses, mean 1.55); 7-column alignment: col 3 1.55, pauses 2–4, no cluster | `EColiPreset_AcceptedByDefault_PerThousandUsage`; MCP `UseCodonFrequencies_ScoresFractions_EColiDefaultScoresPerThousand` |
| R23 | Metamorphic: no per-thousand ⇒ both bases equal | uses `EColiK12 with { PerThousand = null }` (invariant unchanged) | `MinMax_Basis_EqualWithoutPerThousand_AndBounded` |
