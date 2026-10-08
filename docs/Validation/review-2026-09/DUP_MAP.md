# Repo-wide duplication map (orchestrator sweep, 2026-09-28)

> Read-only static sweep; every 'probable bug' must be confirmed against sources + reference implementation before fixing (campaign rule 1). Batch owners fix items inside their ownership; the rest is Phase 2.

# Seqeron — duplicated algorithm implementations map

Scope: `src/Seqeron/Algorithms/**` and `src/Seqeron/Mcp/**` (bin/obj excluded). Read-only survey; paths are relative to `src/Seqeron/`.
Method: a per-name cross-file index of every declared method (2,248 declarations), then grep for characteristic constants and tables
(complement switches, codon tables, Kyte–Doolittle, NN ΔH/ΔS, Lanczos coefficients, `-p*log p`, `/2` medians, DP matrices),
followed by reading each candidate's body side by side.

Legend: **ID** = the same behaviour (safe to delegate). **DIFF** = the behaviour differs (the difference is described). **BUG?** = the difference is probably a defect.
**LEGIT** = the semantics really are different, so it is not a true duplicate.

Overall: the MCP layer delegates almost everywhere. Only one MCP tool re-implements logic (`primer_dimer`, §10). Most duplication is
inside `Algorithms/**`, where each analyzer carries private copies of statistics, genetic code, complement and GC helpers.

---

## 1. Statistical primitives (p-values, special functions, descriptive stats) — HIGHEST IMPACT

Canonical home: `Algorithms/Seqeron.Genomics.Infrastructure/StatisticsHelper.cs` has only `NormalCDF` and `Erf`. Every other primitive is copied privately.

### 1a. Two-sample t-test p-value — **DIFF / BUG?**
| Location | Method | Vis. | Behaviour |
|---|---|---|---|
| **Canonical:** `Annotation/TranscriptomeAnalyzer.cs:415` | `WelchTTestTwoSidedPValue` | private | Welch–Satterthwaite df plus the exact Student-t tail via `RegularizedIncompleteBeta` (:485) |
| `Annotation/TranscriptomeAnalyzer.cs:580` | `CalculateTTestPValue` | private | Welch statistic, but p = 2(1−Φ(t)) (normal approximation). **Anti-conservative for small n** (for example n=3 vs 3) |
| `Metagenomics/MetagenomicsAnalyzer.cs:1626` | `CalculateTTestPValue` | private | Same normal approximation. The degenerate-variance branch also differs: `var1==var2==0` gives 1/0, `se==0` gives 1 |

Impact: `TranscriptomeAnalyzer.AnalyzeDifferentialExpression` (:248, which uses the approximation) backs the MCP tool `differential_expression` (`Mcp/Seqeron.Mcp.Annotation/Tools/AnnotationTools.cs:1813`). That tool's description promises "Welch-style t-test p-values". The same class already has the exact test in `FindDifferentiallyExpressed` (:363).
The same approximation drives `MetagenomicsAnalyzer` differential abundance (:1619).
**Fix:** move `WelchTTestTwoSidedPValue` and `RegularizedIncompleteBeta`/`BetaContinuedFraction` into `StatisticsHelper`, then delete both `CalculateTTestPValue` copies.

### 1b. Hypergeometric / Fisher exact — **DIFF / BUG?**
| Location | Method | Vis. | Behaviour |
|---|---|---|---|
| **Canonical (one-sided ORA):** `Metagenomics/MetagenomicsAnalyzer.cs:1285` | `HypergeometricUpperTail` | public | Exact upper tail Σ P(X=i) in log space |
| **Canonical (two-sided 2×2):** `Oncology/OncologyAnalyzer.DriversArtifactsAnnotation.cs:524` (+`HypergeometricLogProbability` :563) | `FisherExactTwoSided` | private | Exact, log space, O(1)-per-term LogGamma |
| `Annotation/EpigeneticsAnalyzer.cs:793` (+`FisherExactProbability` :770 public) | `FisherExactTwoSided` | private | Exact and ID in result. Uses the O(k) `LogFactorial` loop (§1d), so it is slow at high coverage. Relative 1e-7 tie tolerance (Oncology uses 1e-7 in log space, which is equivalent) |
| `Annotation/TranscriptomeAnalyzer.cs:676` | `CalculateFisherPValue` | private | **Named "Fisher" but is a normal (z) approximation** of the hypergeometric: no continuity correction, P(X>obs) rather than P(X≥obs). Backs `PerformOverRepresentationAnalysis` and the MCP tool at `AnnotationTools.cs:1839` |

**Fix:** move `HypergeometricUpperTail` and `FisherExactTwoSided` into `StatisticsHelper`. The ORA path should call `HypergeometricUpperTail(overlap, background, pathwaySize, deGenes)`.

### 1c. Log-gamma — **DIFF (precision only)**
| Location | Vis. | Algorithm |
|---|---|---|
| `Annotation/TranscriptomeAnalyzer.cs:552` | private | Lanczos g=7, n=9 (≈1e-15) |
| `Oncology/OncologyAnalyzer.DriversArtifactsAnnotation.cs:590` | private | Same as Transcriptome (ID). Shared by all Oncology partials, for example `CtdnaMrdChip.LogChoose` |
| `Metagenomics/MetagenomicsAnalyzer.cs:1315` | private | Numerical Recipes gammln, 6 coefficients (≈2e-10); returns +∞ for x≤0 |
| `Population/PopulationGeneticsAnalyzer.cs:564` | private | Same NR copy as Metagenomics (ID) |

**Canonical:** the Lanczos-g7 version. Put `LogGamma`, `LogChoose` and `LogFactorial` in `StatisticsHelper`.

### 1d. LogChoose / LogFactorial
- `LogChoose`: `Metagenomics/MetagenomicsAnalyzer.cs:1306` (int, via LogGamma), `Oncology/OncologyAnalyzer.CtdnaMrdChip.cs:719` (double, via LogGamma), `Oncology/OncologyAnalyzer.DriversArtifactsAnnotation.cs:571` (int, via LogFactorial). **ID** up to LogGamma precision.
- `LogFactorial`: `Oncology/...DriversArtifactsAnnotation.cs:582` (LogGamma(n+1), O(1)) vs `Annotation/EpigeneticsAnalyzer.cs:831` (Σ log i, **O(k) per call**, called 9× per table × every table in the Fisher sum, so O(n²) per DMR test). **DIFF (perf).**

### 1e. Error function / normal CDF — **DIFF (accuracy)**
| Location | Vis. | Algorithm |
|---|---|---|
| `Infrastructure/StatisticsHelper.cs:19` `Erf` (+`NormalCDF` :11) | public | Abramowitz–Stegun 7.1.26, abs err ≈1.5e-7. `1-NormalCDF(z)` loses all precision for p < ~1e-7 |
| `Oncology/OncologyAnalyzer.CtdnaMrdChip.cs:1119` `Erf` (+`StandardNormalCdf` :1090) | private | Maclaurin series plus Lentz continued fraction, accurate to double precision |

**Fix:** promote the Oncology `Erf` into `StatisticsHelper` and add `NormalSf`/`Erfc` so two-sided tails do not use `1-CDF`. Callers: `TranscriptomeAnalyzer`:599/689, `MetagenomicsAnalyzer`:1645/1768, `PopulationGeneticsAnalyzer`:937.

### 1f. Benjamini–Hochberg — **ID** (same class, two copies)
`Annotation/TranscriptomeAnalyzer.cs:454` `BenjaminiHochbergAdjust` (unsorted input, returns original order) and :605 `BenjaminiHochberg` (requires pre-sorted input). Same numbers. **Fix:** keep `BenjaminiHochbergAdjust`, move it to `StatisticsHelper`, delete the other copy.

### 1g. Median / quantile — **ID**
- Median: `Oncology/OncologyAnalyzer.PhylogenyHeterogeneity.cs:588`, `Oncology/OncologyAnalyzer.SomaticCalling.cs:389`, inline in `IO/QualityScoreAnalyzer.cs:455-458`, `Chromosome/ChromosomeAnalyzer.cs:363-366`, `Annotation/StructuralVariantAnalyzer.cs:1051-1060` (median of non-zero values, a filtered variant). All use the same even/odd rule.
- R-type-7 quantile: `Oncology/OncologyAnalyzer.Signatures.cs:2402` `Percentile` and `Oncology/OncologyAnalyzer.CtdnaMrdChip.cs:1040` `Quantile`. **ID.**
- Population/sample variance is written inline about 10 times (for example `GcSkewCalculator.cs:418` `CalculateVariance`, `ChromosomeAnalyzer.cs:621`, `FastqParser.cs:421`, `MetagenomicsAnalyzer.cs:1060`, `ClinicalMisc.cs:752`, `TranscriptomeAnalyzer.cs:1264`). The divisor mixes N and N−1; each looks intentional for its context, but the choice is undocumented.
**Fix:** add `StatisticsHelper.Median`, `Quantile(type7)` and `Variance(sample: bool)`.

### 1h. Pearson correlation — **DIFF (edge cases)**
| Location | Vis. | Notes |
|---|---|---|
| **Canonical:** `Oncology/ImmuneAnalyzer.cs:1539` `ComputePearsonCorrelation` | internal | Guards against NaN from cancellation. But uses the one-pass Σx² formula, which is less stable |
| `Annotation/TranscriptomeAnalyzer.cs:1082` `CalculatePearsonCorrelation` | public | Two-pass (stable). Returns 0 on zero variance. Backs MCP `pearson_correlation` |
| `Metagenomics/MetagenomicsAnalyzer.cs:822` `PearsonCorrelation` | private | Two-pass, ID to Transcriptome |
| `Metagenomics/MetagenomicsAnalyzer.cs:1002` `TnfPearsonDistance` | private | Inline copy of the same, returns (1−r)/2 |

**Fix:** one `StatisticsHelper.Pearson` (two-pass, plus Immune's NaN guards). `TnfPearsonDistance` should call it.

### 1i. Dense linear algebra (Oncology) — **ID**
- Gaussian elimination: `Oncology/ImmuneAnalyzer.cs:1466` `SolveLinearSystem` (internal, singular if |pivot|<1e-15) vs `Oncology/OncologyAnalyzer.Signatures.cs:605` (private, singular only if pivot == 0.0). **DIFF** only on near-singular input.
- Lawson–Hanson NNLS: `Oncology/ImmuneAnalyzer.cs:1288` `SolveNnls` (+`SolvePassiveSetLeastSquares` :1416) vs `Oncology/OncologyAnalyzer.Signatures.cs:382` `SolveNonNegativeLeastSquares`. Two active-set NNLS solvers in one assembly.

**Fix:** a single internal `LinearAlgebra` helper in the Oncology assembly (or Infrastructure).

### 1j. Singletons that belong in the shared helper (not duplicated yet)
`PopulationGeneticsAnalyzer.cs:468` `ChiSquareCDF` / :486 `RegularizedGammaP`, and `MetagenomicsAnalyzer.cs:1698` `MannWhitneyU`. Moving them now prevents the next copy.

---

## 2. Genetic code / codon→amino-acid tables — HIGH

**Canonical:** `Core/GeneticCode.cs` (`CreateStandardCodonTable` :227, `Translate` :60, `IsStartCodon` :91, `IsStopCodon` :103, NCBI tables 1/2/3/11), plus `Core/Translator.cs`.

| Duplicate | Vis. | Notes |
|---|---|---|
| `MolTools/CodonOptimizer.cs:71` `StandardGeneticCode` (RNA keys, `string` values) + `TranslateCodon` :955 | private | ID to table 1. Single-codon-AA set derived from it (:128) |
| `MolTools/CodonUsageAnalyzer.cs:457` `CodonToAminoAcid` (DNA keys) | private | ID to table 1 |
| `IO/GenBankParser.cs:697` `CodonToAminoAcid` switch + `Translate` :686 | private | ID to table 1. **Ignores `/transl_table`**, so mito/plastid CDS translate wrongly. Should call `Translator.Translate(dna, GeneticCode.GetByTableNumber(n))` |
| `Analysis/GenomicAnalyzer.cs:412-413` `StartCodon="ATG"`, `StopCodons` + `IsStopCodon` :506 | private | Hard-coded standard stops |
| `Annotation/GenomeAnnotator.cs:62-73` `StartCodons {ATG,GTG,TTG}`, `StopCodons` | private | **DIFF**: bacterial start set, while `GeneticCode.Standard` starts are {AUG,UUG,CUG}. `FindOrfs` always uses these hard-coded sets even though `FindOrfsInFrame` receives a `GeneticCode` |

**Fix:** every consumer should take a `GeneticCode` (default `Standard`) and use `Translate` / `IsStartCodon` / `IsStopCodon` / `GetCodonsForAminoAcid`.

## 3. ORF finding — HIGH (three engines, three semantics)
| Location | Vis. | minLength unit | Starts | Nested ORFs | Open 3' ORF | Rev-strand coords |
|---|---|---|---|---|---|---|
| `Core/Translator.cs:74` `FindOrfs` (→:162) | public | amino acids | GeneticCode starts (ATG/TTG/CTG) | no (first start only) | reported | rev-comp frame coords, frame −1..−3, inclusive end |
| `Annotation/GenomeAnnotator.cs:83` `FindOrfs` (→:131) | public | amino acids (excludes stop) | ATG/GTG/TTG | yes (all starts) | only when `requireStartCodon=false` | mapped back to forward coords, half-open end |
| `Analysis/GenomicAnalyzer.cs:441` `FindOpenReadingFrames` (→:468) | public | **nucleotides incl. stop** | ATG only | yes (Rosalind) | no | rev-comp coords, `IsReverseComplement` flag |

MCP exposes two of them: `find_orfs` (Annotation) and `find_open_reading_frames` (Analysis). The same `minLength=100` means 300 nt in one tool and 100 nt in the other. **LEGIT** that the models differ (getorf vs Rosalind vs prokaryotic), but they share no scanning core and no code table.
**Fix:** one parametrised core (starts from `GeneticCode`, `nested`, `allowOpenEnd`, `lengthUnit`, coordinate convention) in `Translator`, with the other two as thin adapters. Document the unit on each MCP tool.

## 4. Codon counting / RSCU / CAI — HIGH
- **Codon counting (4 copies):**
  - `MolTools/CodonUsageAnalyzer.cs:16/25/33` `CountCodons`: ACGT-only, frame 0, upper-case.
  - `MolTools/CodonOptimizer.cs:892` `CalculateCodonUsage`: RNA keys, **no ambiguity filter**, depends on `SplitIntoCodons`.
  - `Annotation/GenomeAnnotator.cs:1080` `GetCodonUsage(string)`: ID to CodonUsageAnalyzer.
  - `Analysis/SequenceStatistics.cs:682` `CalculateCodonFrequencies`: frequencies plus `readingFrame`.

  MCP exposes all four (`count_codons`, `codon_usage`, `codon_frequencies`, plus `rscu`).
- **RSCU (2):** `MolTools/CodonUsageAnalyzer.cs:84` `CalculateRscuCore` vs `Annotation/GenomeAnnotator.cs:1137` `GetCodonUsage(IEnumerable<string>, GeneticCode)`. Same formula n·x/Σ. **DIFF:** CodonUsageAnalyzer also computes RSCU for the stop "family" and uses the standard code only. GenomeAnnotator excludes stops and honours `GeneticCode`, so it is the better canonical.
- **CAI (2):** `MolTools/CodonUsageAnalyzer.cs:119/138` `CalculateCai` (reference RSCU; excludes Met/Trp/stop always; **skips w=0 codons**) vs `MolTools/CodonOptimizer.cs:459` `CalculateCAI` (frequency table; Met/Trp are **included by default**, since `excludeSingleCodonAminoAcids=false`; w is clamped to 1e-6, so absent codons are heavily penalised, while CodonUsageAnalyzer silently drops them). **DIFF.** For the same gene the two public CAIs disagree (MCP `codon_adaptation_index` vs `cai_from_organism_table`). Neither applies Sharp & Li's 0.5 pseudo-count for w=0.

**Fix:** one codon-count core (ACGT filter, frame, `GeneticCode`) in `CodonUsageAnalyzer`. RSCU should come from GenomeAnnotator's code-aware version moved there. Pick one CAI policy and make the other a wrapper.

## 5. GC content — MEDIUM-HIGH (API unit hazard)
**Canonical:** `Core/SequenceExtensions.cs:24/61` `CalculateGcContent(Fast)` (percent) and :77/113 `CalculateGcFraction(Fast)` (fraction). Both count A,C,G,T,U, are case-insensitive and exclude other symbols from both counts. Also :177/216 Biopython ambiguity modes.

Methods named `CalculateGcContent` that delegate but **return different units**:
- percent: `MolTools/PrimerDesigner.cs:236` (public), `MolTools/CrisprDesigner.cs:969`, `Core/DnaSequence.cs:79`, `Core/RnaSequence.cs:92` `GcContent()`, MCP `gc_content`.
- fraction: `Annotation/MiRnaAnalyzer.cs:2709` (**public**), `Annotation/EpigeneticsAnalyzer.cs:354`, `Metagenomics/MetagenomicsAnalyzer.cs:1025`, `MolTools/CodonOptimizer.cs:960`, `MolTools/ProbeDesigner.cs:1377`.

So `PrimerDesigner.CalculateGcContent("GC")` returns 100 while `MiRnaAnalyzer.CalculateGcContent("GC")` returns 1. **Fix:** rename the fraction wrappers to `CalculateGcFraction`, or inline the canonical call.

Re-implementations (not delegating):
| Location | Vis. | DIFF |
|---|---|---|
| `Analysis/GcSkewCalculator.cs:387` `CalculateGcContent(seq, fraction)` | private | Counts only upper-case A/C/G/T. **U excluded**, so RNA input gives a wrong denominator. Callers pass upper-cased input, so case is fine |
| `Analysis/SequenceStatistics.cs:899` `CalculateGcContentProfile` (inline :915-930) | public | ACGTU. It duplicates `GcSkewCalculator.CalculateWindowedGcContentCore` (:364), and the two **disagree on RNA** |
| `MolTools/PrimerDesigner.cs:1747` `GcFraction` | private | ACGT only (fine for the DNA-only Owczarzy correction). Same as the canonical for DNA |
| `IO/FastqParser.cs:344` `CalculateStatistics` (GC inline :372-376) | public | **Denominator = totalBases incl. N**. Differs from the canonical |
| `Chromosome/GenomeAssemblyAnalyzer.cs:202/210` `CountGC`/`CountBases` | private | Denominator = all non-N chars (IUPAC and gaps included). DIFF |
| `MolTools/CodonOptimizer.cs:963` `GetCodonGcContent` | private | Per-codon; fine |
| `MolTools/CrisprDesigner.cs:557`, `MolTools/AzimuthRuleSet2.cs:179` | private | Fixed-window GC counts inside validated ACGT-only feature code. **LEGIT** |

## 6. Melting temperature & NN thermodynamics — MEDIUM-HIGH
**Basic Tm (Wallace / Marmur–Doty / salt):** constants are centralised in `Infrastructure/ThermoConstants.cs`, but the counting is copied 3×:
| Location | Vis. | DIFF |
|---|---|---|
| **Canonical:** `MolTools/PrimerDesigner.cs:193` `CalculateMeltingTemperature` | public | Counts valid ACGT. Wallace threshold uses the valid length. Clamps ≥0 |
| `Analysis/SequenceStatistics.cs:563` `CalculateMeltingTemperature(seq, useWallaceRule)` | public | Wallace threshold uses the **raw length** (N and gaps count), no ≥0 clamp. MCP `melting_temperature` uses this, while MCP `primer_melting_temperature` uses PrimerDesigner, so they can disagree on the same input |
| `MolTools/ProbeDesigner.cs:1380` `CalculateTm` | private | Case-sensitive counts. Wallace uses raw length. Long branch uses the **salt-adjusted** formula (81.5+16.6log…), not Marmur–Doty |

**Nearest-neighbour ΔH/ΔS tables (3 copies):**
| Location | Parameter set | Notes |
|---|---|---|
| **Canonical:** `MolTools/PrimerDesigner.cs:495` `NnUnifiedParams` + :603 `CalculateNearestNeighborThermodynamics`, :660 `CalculateMeltingTemperatureNN` | SantaLucia & Hicks 2004 (Biopython DNA_NN4; AA = −7.6/−21.3) | Owczarzy 2004/2008 salt. **Doc bug:** comments at :466 and :492 say "identical to SantaLucia 1998" / "SantaLucia (1998) unified", but the 1998 AA/TT values are −7.9/−22.2 (DNA_NN3) |
| `Analysis/SequenceStatistics.cs:429` `NearestNeighborParams` + :496 `CalculateThermodynamics` | Allawi & SantaLucia 1997 (DNA_NN3), terminal-init model, salt method 5, F=4 | DIFF: a different parameter set and salt model from PrimerDesigner's NN Tm. Both are public "NN Tm" |
| `MolTools/AzimuthRuleSet2.cs:239/262` `MeltingTemp`/`NearestNeighbor` | DNA_NN3, reproduces Biopython `Tm_NN` defaults | **LEGIT**: it must match the Azimuth reference features bit for bit |

(`NtthalDimer`/`NtthalHairpin` thal tables are a separate Primer3 port. **LEGIT**.)
**Fix:** `SequenceStatistics.CalculateThermodynamics` should call PrimerDesigner's NN core with an explicit `NnParameterSet` (NN3/NN4) and salt-mode argument. `ProbeDesigner.CalculateTm` and `SequenceStatistics.CalculateMeltingTemperature` should call `PrimerDesigner.CalculateMeltingTemperature` or state why they differ.

## 7. Protein / nucleotide physico-chemical properties — MEDIUM-HIGH
- **Molecular weight (protein):** `Core/ProteinSequence.cs:127` `MolecularWeight()` (2-dp masses, water 18.015, **rounded to 2 dp**) vs **canonical** `Analysis/SequenceStatistics.cs:181` `CalculateMolecularWeight` (Biopython 4-dp masses, water 18.0153, unrounded). **DIFF**, small but visible.
- **Isoelectric point:** `Core/ProteinSequence.cs:157` `IsoelectricPoint()` uses **Lehninger-type pKa** (N 9.69, C 2.34, Cys 8.3, His 6.0, Lys 10.5). **Canonical** `Analysis/SequenceStatistics.cs:~300` `CalculateIsoelectricPoint` uses EMBOSS (N 8.6, C 3.6, Cys 8.5, His 6.5, Lys 10.8). **DIFF**: the two public pI APIs give different answers for the same protein.
- **Kyte–Doolittle scale copied 4×:** `Analysis/SequenceStatistics.cs:353`, `Analysis/DisorderPredictor.cs:52`, `Analysis/ProteinMotifFinder.cs:802`, `Core/ProteinSequence.cs:~235` (local inside `Gravy()`). Values ID.
  - GRAVY: `SequenceStatistics.cs:372` `CalculateHydrophobicity` and `DisorderPredictor.cs:905` `CalculateHydropathy` are **ID**. `ProteinSequence.cs:231` `Gravy()` is DIFF (rounded to 3 dp).
  - Window profile: `SequenceStatistics.cs:398` `CalculateHydrophobicityProfile` divides by **windowSize** (unknown residues count as 0), while `ProteinMotifFinder.cs:815` `CalculateHydropathyProfile` divides by the **count of recognised residues** (Biopython). **DIFF** for X/B/Z.
- **Nucleotide MW:** **canonical** `Analysis/SequenceStatistics.cs:~218` `CalculateNucleotideMolecularWeight` (Biopython tables, separate DNA/RNA) vs `MolTools/ProbeDesigner.cs:1311` `CalculateMolecularWeight`. ProbeDesigner mixes DNA and RNA in one table, uses **U = 308.2** (that is dUMP; UMP is 324.18), gives unknown bases 330, and uses water 18.0. **BUG?** for RNA probes.

**Fix:** one internal `AminoAcidScales` class (KD, masses, EMBOSS pKa) in Core. `ProteinSequence` members should delegate to `SequenceStatistics`, or both should call the shared scales. ProbeDesigner should call `SequenceStatistics.CalculateNucleotideMolecularWeight(seq, isDna)`.

## 8. Complement / reverse complement — MEDIUM
**Canonical:** `Core/SequenceExtensions.cs:247` `GetComplementBase` (full IUPAC, case-insensitive, pass-through), :284 `GetRnaComplementBase`, :313/329 span versions. `Core/DnaSequence.cs:146` `GetReverseComplementString`.
| Duplicate | Vis. | DIFF |
|---|---|---|
| `Core/RnaSequence.cs:52` `Complement()`, :72 `ReverseComplement()` | public | ACGU only. IUPAC codes pass through **uncomplemented** (R stays R). Should use `GetRnaComplementBase` |
| `Core/ISequence.cs:131` `IupacDnaSequence._complements` (+`GetComplement` :~180, `GetReverseComplement` :~201) | public | Case-sensitive. Unknown characters become N (the canonical passes them through). Otherwise ID |
| `Core/ISequence.cs:404` `QualitySequence.GetComplement` | public | Uses the canonical, then maps "unchanged" to N. **BUG:** self-complementary codes S, W, N and gaps become N |
| `Annotation/MiRnaAnalyzer.cs:293` `GetReverseComplement` | public | RNA, T→A, unknown → N (backs MCP `rna_reverse_complement`) |
| `Annotation/MiRnaAnalyzer.cs:1490` `ReverseComplementRna` | private | Same as above but unknown characters pass through. Two copies in one class that differ from each other |
| `MolTools/PrimerDesigner.cs:1729` `Complement(string)` | private | ACGT, pass-through. ID to the canonical for upper-case DNA |
| `Metagenomics/MetagenomicsAnalyzer.cs:755` `ExtendWithReverseComplement` | private | Inline RC on an ACGT-filtered string. ID |
| `Oncology/OncologyAnalyzer.Signatures.cs:146` `Complement(char)`; `MolTools/CrisprDesigner.cs:879` `CfdComplement` | private | Throw on non-ACGT. **LEGIT** strict validators (SBS96 pyrimidine folding, CFD lookup) |

**Fix:** RnaSequence and the MiRnaAnalyzer copies should use `GetRnaComplementBase`. `IupacDnaSequence` should use `GetComplementBase`. Fix the `QualitySequence` N-mapping.

## 9. RNA base-pair predicates — MEDIUM (MCP-visible divergence)
| Location | Vis. | Accepts T? | Notes |
|---|---|---|---|
| **Canonical:** `Analysis/RnaSecondaryStructure.cs:418` `CanPair` (lookup :400), `GetBasePairType` :428 | public | **No** (G·T → false) | backs MCP analysis `can_pair` (`AnalysisTools.cs:1334`) |
| `Annotation/MiRnaAnalyzer.cs:322` `CanPair`, :336 `IsWobblePair` (+`NormalizeBase` :347 T→U) | public | **Yes** | backs MCP annotation `can_pair` (`AnnotationTools.cs:1103`) |
| `Analysis/RnaSecondaryStructure.cs:1763` `PairType` | private | no | third copy inside the same file (MFE inner loop, byte codes) |
| `Analysis/RnaSecondaryStructure.cs:829` `IsAUorGU` | private | no | **ID** to `Core/Turner2004Parameters.cs:118` `IsTerminalPenaltyPair` (MiRna already delegates, :1734) |
| `MolTools/CodonOptimizer.cs:636` `AreComplementary(char,char)` | private | no | WC only, **no G·U wobble**, for mRNA secondary-structure avoidance. **DIFF** (understates mRNA hairpins) |

Two MCP tools both named `can_pair` return different results for inputs containing T.
**Fix:** `MiRnaAnalyzer.CanPair`/`IsWobblePair` should be `RnaSecondaryStructure.CanPair(NormalizeBase(a), NormalizeBase(b))`. Replace `IsAUorGU` with the Turner helper. CodonOptimizer should use `RnaSecondaryStructure.CanPair`.

## 10. DNA complementarity / primer-dimer — MEDIUM (the only MCP re-implementation)
- `MolTools/PrimerDesigner.cs:394` `HasPrimerDimer` (public) and `IsComplementary` :2002 (private).
- `Mcp/Seqeron.Mcp.MolTools/Tools/MolToolsTools.cs:120-139` (`primer_dimer` tool) **re-implements** the complementarity-count loop and has its own local `IsComplementary`, "mirroring the inner loop in HasPrimerDimer". **ID** today, but the copy will drift. **Fix:** expose `PrimerDesigner.CountPrimerDimerComplementarity(p1,p2)` (or return it from `HasPrimerDimer`) and have the MCP tool call it.
- **BUG?** Both copies compare `end1 = suffix(primer1)` against `end2 = prefix(revcomp(primer2))` using *complementarity*. `end2` is already the reverse complement, so a real 3′–3′ dimer requires `end1 == end2`. The current test only detects parallel complementarity, so it flags poly-A vs poly-A (which cannot pair) as a dimer and misses GGCC/GGCC 3′ ends. The unit tests encode the current behaviour (`PrimerDesignerTests.cs:106`, `PrimerDesigner_PrimerDesign_Tests.cs:332`). Verify against Primer3 before changing it. The thermodynamic `CalculateDimerThermodynamicsNtthal` (:1633) is the correct engine.
- `PrimerDesigner.cs:1825` `IsSelfComplementary` and `:1991` `AreComplementary(string,string)` are local and correct.
- Self-complementarity heuristics: `MolTools/CrisprDesigner.cs:977` (offset-sum / L²) vs `MolTools/ProbeDesigner.cs:1420` (position-wise seq==rc fraction). Same name, **LEGIT-different ad-hoc metrics**. Both could be replaced by `PrimerDesigner.CalculateHairpin/SelfDimer*` ΔG.

## 11. IUPAC matching / encoding — MEDIUM
**Canonical:** `Core/IupacHelper.cs:16` `MatchesIupac` (15 codes, throws on invalid).
| Duplicate | Vis. | DIFF |
|---|---|---|
| `Analysis/MotifFinder.cs:47` `IupacCodes` dict (used by `FindDegenerateMotif(DnaSequence)` :~97, lookup :113) **and** a switch inside the string/cancellation overload :179-194 | private | Two copies in one class. Both ID to `IupacHelper` |
| `Chromosome/ChromosomeAnalyzer.cs:935` `MatchesIupac(seq,start,pattern)` | private | Handles only R/Y. Adequate for the fixed CENP-B consensus, but a latent trap. Should loop over `IupacHelper.MatchesIupac` |
| `Core/ISequence.cs:152` `_expansions` / `ExpandCode` :215 / `CodesMatch` :270 | public | Symmetric code-vs-code overlap (**LEGIT** semantics). **BUG?** `CodesMatch` does not upper-case its arguments, so MCP `iupac_match` (`SequenceTools.cs:648`) returns false for `("r","a")` |
| `CrisprDesigner.cs:150`, `RestrictionAnalyzer.cs:249` `MatchesIupac` | private | Delegates (fine) |
| Base-set → IUPAC code: `Core/ISequence.cs:225` `GetIupacCode(IEnumerable<char>)` vs `Analysis/MotifFinder.cs:462` `GetIupacCode(counts,total)` | public / private | MotifFinder adds a threshold step, then uses its own string→code switch. Should call `IupacDnaSequence.GetIupacCode(present)` |

## 12. Shannon entropy — MEDIUM (two MCP tools, two answers)
| Location | Vis. | Alphabet / notes |
|---|---|---|
| `Analysis/SequenceComplexity.cs:83` (→:89) `CalculateShannonEntropy` | public | **ACGT only** (U, N and IUPAC ignored), so RNA input gives 0 bits. MCP `complexity_shannon` (`SequenceTools.cs:~466`) |
| `Analysis/SequenceStatistics.cs:726` `CalculateShannonEntropy` | public | **any letter** (N, U and amino acids counted). MCP `shannon_entropy` (`SequenceTools.cs:~313`) |
| `Analysis/DisorderPredictor.cs:367` | private | Window over A–Z, protein. **LEGIT** (protein) |
| `Analysis/ProteinMotifFinder.cs:1194` `CalculateSegComplexity` | private | SEG window, protein. **LEGIT** |
| k-mer entropy: `Analysis/KmerAnalyzer.cs:325` `CalculateKmerEntropy` vs `Analysis/SequenceComplexity.cs:153/160` `CalculateKmerEntropy` | public / public | **ID** (both upper-case, overlapping k-mers). MCP exposes both (`kmer_entropy`, `complexity_kmer_entropy`) |
| ln-based diversity inline: `Metagenomics/MetagenomicsAnalyzer.cs:506` `CalculateShannonIndex`, `MetagenomicsAnalyzer.cs:1370` (functional diversity), `Oncology/OncologyAnalyzer.PhylogenyHeterogeneity.cs:563` | private | Natural log. **LEGIT** (ecological H′), but all three should share one `Entropy(counts, base)` helper |

**Fix:** one `Entropy(IEnumerable<int> counts, double logBase)` helper. `SequenceComplexity` gets an explicit alphabet parameter (ACGT vs ACGU). `SequenceComplexity.CalculateKmerEntropy` should delegate to `KmerAnalyzer`.

## 13. Linguistic complexity — MEDIUM (same name, three formulas)
- `Analysis/SequenceComplexity.cs:29` (→:34): Σobserved/Σpossible over k=1..10 (Trifonov product form).
- `Analysis/SequenceStatistics.cs:761`: **mean of the per-k ratios**, k=1..6.
- `Chromosome/GenomeAssemblyAnalyzer.cs:1052` (private): single k=4, skips N-containing k-mers.

All three are called "linguistic complexity" but return different numbers. **Fix:** pick SequenceComplexity as canonical. SequenceStatistics should delegate (or be renamed "mean vocabulary usage"). GenomeAssemblyAnalyzer should call it with `maxWordLength=4`, or be renamed.

## 14. Phred ↔ error probability & encoding detection — MEDIUM
| Location | Vis. | p≤0 | p≥1 | rounding |
|---|---|---|---|---|
| `Core/ISequence.cs:481/487` (QualitySequence) | public | 93 | 0 | **truncation** (byte cast) |
| `IO/FastqParser.cs:206/214` | public | 93 | **negative for p>1** | round |
| `IO/QualityScoreAnalyzer.cs:359/367` | public | **60** | 0 | round |

All three `PhredToErrorProbability` are ID. **Fix:** canonical = QualityScoreAnalyzer's clamping plus FastqParser's 93 cap. The other two delegate.

Encoding detection:
| Location | Vis. | Rule |
|---|---|---|
| `IO/FastqParser.cs:143` | public | **Order-dependent**: the first char <'@' gives Phred33, the first char >'I'(73) gives Phred64. Illumina 1.8 Q41 ('J'=74) is **misdetected as Phred64** when it appears before any char <64. **BUG?** |
| `IO/QualityScoreAnalyzer.cs:262` | public | min<59 → 33; min≥64 && max>74 → 64 |
| `IO/QualityScoreAnalyzer.cs:307` (multi-read, with confidence) | public | **Canonical**: min<64 → 33; max>74 → 64 |

Stats: `IO/FastqParser.cs:344` `CalculateStatistics` overlaps with `IO/QualityScoreAnalyzer.cs:377/404` (Q20/Q30/mean). **Fix:** FastqParser should delegate encoding detection and quality statistics to QualityScoreAnalyzer.

## 15. N50 — MEDIUM (bug)
- **Canonical:** `Chromosome/GenomeAssemblyAnalyzer.cs:233/272/290` `CalculateNx`/`CalculateN50`. Integer-exact `cum*100 ≥ total*x` (QUAST), long totals.
- `Alignment/SequenceAssembler.cs:570` `CalculateStats`: `halfLength = totalLength/2` (integer division) and an `int` total. **BUG:** for odd totals it stops one contig early. For lengths {4,3,2} it gives N50=4 (canonical: 3). It also overflows above 2^31 bp. Backs MCP `assembly_stats` / `assemble_*` (`AlignmentTools.cs:332`), while the Chromosome MCP `assembly_statistics` uses the canonical.
- **Fix:** `SequenceAssembler.CalculateStats` should call `GenomeAssemblyAnalyzer.CalculateN50` (or move Nx to a shared place, because Alignment does not reference Chromosome).

## 16. Hamming / identity / k-mer similarity — LOW-MEDIUM
- Hamming: `Core/SequenceExtensions.cs:388` (span) and `Alignment/ApproximateMatcher.cs:180` (string) are **ID** (case-insensitive, throw on unequal length). `MolTools/CrisprDesigner.cs:400` `CountMismatches` is case-sensitive and returns `int.MaxValue` on unequal length. `ApproximateMatcher.cs:405` `HammingDistanceFast` adds the length difference (**LEGIT** different). **Fix:** ApproximateMatcher should delegate to the span version, and CRISPR should call it after its ACGT validation.
- Ungapped identity:
  - `Alignment/SequenceAssembler.cs:243` `CalculateIdentity`: case-insensitive, 0 if lengths differ.
  - `Metagenomics/PanGenomeAnalyzer.cs:311` `CalculateSequenceIdentity`: case-sensitive, denominator is the shorter length.
  - `Alignment/SequenceAligner.cs:2115` `PercentIdentity`: gapped, percent, rounded (**LEGIT**, T-Coffee weight).
  - **DIFF** between the first two.
- k-mer Jaccard:
  - `Analysis/GenomicAnalyzer.cs:384` `CalculateSimilarity`: percent; `GetKmers` :522 is case-sensitive, but DnaSequence input is already upper-case.
  - `Analysis/ComparativeGenomics.cs:516` `CalculateSequenceSimilarity`: fraction, k=5, plus coverage.
  - Same Jaccard core. **Fix:** a shared `KmerAnalyzer.JaccardSimilarity(a,b,k)`.
- Levenshtein exists once (`ApproximateMatcher.cs:203`). NW/SW/semi-global DP loops exist only in `SequenceAligner` (`:109`, `:315`, `:423`, `:1390`, `:2208`). `AnchorBasedAligner` and `MiRnaAnalyzer.AlignMiRnaToTarget` use their own models (**LEGIT**). **No pairwise-DP duplication found.**

## 17. Repeats & palindromes — LOW-MEDIUM
- Palindromes: `Analysis/GenomicAnalyzer.cs:179` `FindPalindromes` is a **verbatim copy** of `Analysis/RepeatFinder.cs:973/988` `FindPalindromesCore` (same even-length loop). GenomicAnalyzer lacks the validation: an odd `minLength` silently scans only odd lengths and returns nothing, and there is no null check. **Fix:** delegate to RepeatFinder.
- Tandem repeats (4 engines):
  - `Analysis/GenomicAnalyzer.cs:111/124`: exact; reports non-primitive units such as "ATAT" alongside "AT"; not left-maximal.
  - `Annotation/GenomeAnnotator.cs:884`: private; primitive and left-maximal, deduplicated.
  - `Chromosome/GenomeAssemblyAnalyzer.cs:745`: 80% purity, approximate.
  - `Analysis/RepeatFinder.cs:22-82` `FindMicrosatellites`.
  - GenomeAnnotator's exact core is the most correct. GenomicAnalyzer should adopt it. The purity-based search is **LEGIT** different.
- Inverted repeats:
  - `Analysis/RepeatFinder.cs:697`: DNA, loop ≥ minLoop.
  - `Annotation/GenomeAnnotator.cs:958`: private, DNA, gap ≥ 0.
  - `Analysis/RnaSecondaryStructure.cs:2191`: RNA pairing rules. **LEGIT**.
  - GenomeAnnotator should call RepeatFinder with `minLoopLength=0`.

## 18. Inline k-mer extraction/counting — LOW
`KmerAnalyzer.CountKmers` (:14) and `SequenceExtensions.CountKmersSpan` (:349) are canonical. The same `Substring(i,k)` + dictionary/set loop is re-written in:
- `Analysis/MotifFinder.cs:544,625`
- `Chromosome/GenomeAssemblyAnalyzer.cs:589,684,708,878,910,928,1062`
- `Chromosome/ChromosomeAnalyzer.cs:588`
- `Alignment/SequenceAssembler.cs:1140`
- `Alignment/ApproximateMatcher.cs:342`
- `Annotation/GenomeAnnotator.cs:798`
- `Analysis/ComparativeGenomics.cs:527,1291`
- `Analysis/SequenceComplexity.cs:169,375`
- `Analysis/SequenceStatistics.cs:776`
- `Metagenomics/MetagenomicsAnalyzer.cs:223,380` (canonical k-mers, ACGT filter)

Mostly **ID** apart from ACGT/N filtering. Replace with `KmerAnalyzer.CountKmers` / a `DistinctKmers` helper where no filter is needed.

## 19. Transition / transversion — LOW
- `Annotation/VariantCaller.cs:185` `ClassifyMutation` (purine test) + :204 `CalculateTiTvRatio` returns **0** when there are no transversions.
- `IO/VcfParser.cs:646` `CalculateTiTvRatio` (pair set) returns **null** in that case. VcfParser counts any non-transition SNP (including N) as a transversion; VariantCaller treats N as a pyrimidine.
- `Phylogenetics/PhylogeneticAnalyzer.cs:272` `IsTransition`.

**Fix:** a shared `Nucleotide.IsTransition(a,b)` in Core. Agree on one no-transversion sentinel.

## 20. Synteny & rearrangements — LOW (conceptual duplicates, different models)
`Chromosome/ChromosomeAnalyzer.cs:1341` `FindSyntenyBlocks` (coordinate-gap heuristic) + :1419 `DetectRearrangements` vs `Analysis/ComparativeGenomics.cs:126` `FindSyntenicBlocks` (MCScanX DP) + :579 `DetectRearrangements` (breakpoint model). MCP exposes both pairs. **Mostly LEGIT** (different inputs and models), but users get two synteny answers. Consider making the Chromosome versions adapters over ComparativeGenomics.

## 21. k-means-style clustering — LOW
`Metagenomics/MetagenomicsAnalyzer.cs:916` `KMeansCluster` (multi-dim), `Oncology/OncologyAnalyzer.Clonality.cs:376` `ClusterCcfValues` (+:452/:481, 1-D), `Annotation/TranscriptomeAnalyzer.cs:1138` `ClusterGenesByExpression` ("k-means-like", correlation, single pass). Three Lloyd loops. **Fix:** a shared `KMeans(double[][] points, k, seed)`.

## 22. Trivial helpers — LOW
`Sigmoid` (`MiRnaAnalyzer.cs:2432`, `MhcflurryAffinityPredictor.cs:484`). `PercentScale = 100` constants (`GcSkewCalculator.cs:385`, `GenomicAnalyzer.cs:357`, `SequenceStatistics.cs:872`). `NormalizeBase` (MiRna T→U vs Oncology validate-ACGT, **LEGIT**-different).

---

## Legitimately different (not true duplicates — keep separate, maybe rename)
- `AzimuthRuleSet2.MeltingTemp`/`NearestNeighbor` (DNA_NN3 by design, to reproduce the reference model); `NtthalDimer`/`NtthalHairpin` (Primer3 thal port).
- `MiRnaAnalyzer.CalculateHairpinEnergy` (:1691). This is a simplified perfect-stem hairpin energy built on the shared `Turner2004Parameters`. It is **not** the full `RnaSecondaryStructure.CalculateHairpinLoopEnergy` (:663): no special loops, UU/GA/GG bonuses or all-C penalty. That is acceptable as a quick filter; the MFE path (`FindPreMiRnaHairpinsByMfe`) already uses `RnaSecondaryStructure`.
- `Oncology Complement` / `CrisprDesigner.CfdComplement` (strict ACGT validators).
- Protein entropy windows (DisorderPredictor, ProteinMotifFinder SEG) vs DNA Shannon entropy; ln-based ecological Shannon vs bit-based sequence entropy.
- BED `MergeOverlapping` (`IO/BedParser.cs:328`, half-open, merges book-ended intervals) vs GFF (`IO/GffParser.cs:509`, closed, merges adjacent intervals with +1). Correct for their coordinate systems. `BedParser.cs:57` `Overlaps` (half-open) vs `TranscriptomeAnalyzer.cs:923` (closed). Correct.
- `ApproximateMatcher.HammingDistanceFast` (length-penalising), `SequenceAligner.PercentIdentity` (gapped T-Coffee weight).
- `CodesMatch` (code-vs-code) vs `MatchesIupac` (base-vs-code).
- ORF engines differ in model (§3), but they should share one core.

## MCP layer summary
- Re-implements logic: **only** `Mcp/Seqeron.Mcp.MolTools/Tools/MolToolsTools.cs:120-139` (`primer_dimer` complementarity count).
- Minor redundant recount: `Mcp/Seqeron.Mcp.Sequence/Tools/SequenceTools.cs:365-367` (`gc_content` recounts G/C and valid bases next to the canonical call; consistent with it).
- Exposes divergent duplicates side by side, so fixing §1a, §3, §4, §6, §9, §12 and §15 changes MCP output:
  - `can_pair` ×2
  - `shannon_entropy` vs `complexity_shannon`
  - `kmer_entropy` vs `complexity_kmer_entropy`
  - `find_orfs` vs `find_open_reading_frames`
  - `melting_temperature` vs `primer_melting_temperature`
  - `assembly_stats` vs `assembly_statistics` (N50)
  - four codon-usage tools
  - `differential_expression` (normal-approximation p)
  - `over_representation_analysis` (normal-approximation "Fisher")

## Suggested consolidation order
1. `StatisticsHelper`: exact Welch t, hypergeometric/Fisher, accurate Erf/Erfc, LogGamma/LogChoose, BH, Median/Quantile, Pearson. This fixes the DE and ORA p-values.
2. N50 in `SequenceAssembler` → `CalculateNx`.
3. Genetic code: CodonOptimizer, CodonUsageAnalyzer, GenBankParser (honour `/transl_table`), GenomeAnnotator, GenomicAnalyzer → `GeneticCode`.
4. Codon count / RSCU / CAI unification; ORF core unification.
5. Protein scales / pI / MW. Nucleotide MW in ProbeDesigner.
6. RNA `CanPair` unification; complement helpers (RnaSequence, MiRna, QualitySequence bug, IupacDnaSequence).
7. GC-content naming (percent vs fraction). FASTQ encoding detection and Phred conversions → QualityScoreAnalyzer.
8. The rest (entropy, linguistic complexity, repeats, Jaccard, Ti/Tv, k-means, inline k-mers).
