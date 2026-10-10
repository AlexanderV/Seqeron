# seqeron-oncology — curated Method-ID map (C#-API only)

There is **no MCP oncology server** and Oncology is **not** in `docs/skills/domain-map.json`, so there
is **no generated tool slice** — this hand-curated index is the map. Every entry is a real
`Seqeron.Genomics.Oncology` symbol verified against source; paths are from repo root `../../../../`.

`OncologyAnalyzer` is one `static partial` class split over `OncologyAnalyzer.*.cs` partial files, so
the **file** column names the partial that holds the method (line numbers drift — locate a method with
`grep -n "public static .* <Name>(" src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/*.cs`).
All three classes share namespace **`Seqeron.Genomics.Oncology`**: `OncologyAnalyzer`,
`ImmuneAnalyzer` (`ImmuneAnalyzer.cs`), `MhcflurryAffinityPredictor` (`MhcflurryAffinityPredictor.cs`).
Below, unqualified members are on `OncologyAnalyzer`; `src/.../` = `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/`;
a file written `CopyNumberPloidy` means `OncologyAnalyzer.CopyNumberPloidy.cs`.

## Purity / ploidy / allele-specific copy number
| Method ID | file | Algorithm doc |
|---|---|---|
| `FitPurityPloidy` / `TryFitPurityPloidy` / `EvaluatePurityPloidy` (+ `AscatSexModel` overloads) | `CopyNumberPloidy` | `Allele_Specific_Copy_Number_Derivation.md`, `Tumor_Purity_Estimation.md`, `Tumor_Ploidy_Estimation.md` |
| `FitPurityPloidyFromAspcf` / `TryFitPurityPloidyFromAspcf` / `EvaluatePurityPloidyFromAspcf` | `CopyNumberPloidy` | `Allele_Specific_Copy_Number_Derivation.md` |
| `EstimatePurity` | `SomaticCalling` | `Tumor_Purity_Estimation.md` |
| `EstimatePurityFromVAF` / `EstimatePurityFromVaf` (scalar) | `SomaticCalling` | `Tumor_Purity_Estimation.md` |
| `AnalyzePurityPeaks` / `AnalyzeComplexKaryotypePeaks` / `AnalyzeSubclonalPurityPeaks` | `CopyNumberPloidy` | `Tumor_Purity_Estimation.md` |
| `FitBinomialMixture` | `CopyNumberPloidy` | `Tumor_Purity_Estimation.md` |
| `AdjustVAFForPurity` | `SomaticCalling` | `Variant_Allele_Frequency.md` |
| `EstimatePloidy` (bp-weighted; `probeCounts` overload) | `CopyNumberPloidy` | `Tumor_Ploidy_Estimation.md` |
| `DetectWholeGenomeDoubling` (sample span; `ReferenceGenome` overload) / `DetectWholeGenomeDoublingFromSuppliedLength` | `CopyNumberPloidy` | `Tumor_Ploidy_Estimation.md` |
| `ComputeAscatGenomeMetrics` / `GetAutosomeLengths` / `GetAutosomalGenomeLength` | `CopyNumberPloidy` | `Tumor_Ploidy_Estimation.md` |
| `SegmentAlleleSpecificAspcf` (summary; germline-aware; `AscatMaleXGenotyping` overloads) | `CopyNumberPloidy` | `Allele_Specific_Copy_Number_Derivation.md` |
| `SegmentAlleleSpecificAsMultiPcf` | `CopyNumberPloidy` | `Allele_Specific_Copy_Number_Derivation.md` |
| `SegmentAlleleSpecific` (= `SegmentAlleleSpecificAspcf(loci, 70)`; thresholds ignored) | `CopyNumberPloidy` | `Allele_Specific_Copy_Number_Derivation.md` |
| `DeriveMultiplicity` | `CopyNumberPloidy` | `Allele_Specific_Copy_Number_Derivation.md` |
| `FitSubclonalCopyNumber` / `FitSubclonalCopyNumberWithSnpTest` / `FitSubclonalCopyNumberWithBootstrap` | `CopyNumberPloidy` | `Allele_Specific_Copy_Number_Derivation.md` |
| `SegmentPhasedBaf` / `BuildBattenbergSegments` | `CopyNumberPloidy` | `Allele_Specific_Copy_Number_Derivation.md` |
| `MergeBattenbergSegments` / `MaskHighCopyNumberSegments` / `CallBattenbergSubclones` | `CopyNumberPloidy` | `Allele_Specific_Copy_Number_Derivation.md` |

## Copy-number calling / classification
| Method ID | file | Algorithm doc |
|---|---|---|
| `Log2RatioToCopyNumber` (+ `purity` overload) | `CopyNumberPloidy` | `Copy_Number_Alteration_Classification.md` |
| `CallCopyNumber` / `ClassifyCopyNumber` / `ClassifyCopyNumbers` (+ `purity` overloads) | `CopyNumberPloidy` | `Copy_Number_Alteration_Classification.md` |
| `DetectFocalAmplifications` / `IsFocalAmplification` | `CopyNumberPloidy` | `Focal_Amplification_Detection.md` |
| `DeconstructZiggurat` / `DetectFocalAmplificationEvents` | `CopyNumberPloidy` | `Focal_Amplification_Detection.md` |
| `IdentifyAmplifiedOncogenes` (arm labels; `GeneLocus` locus-overlap overload) | `CopyNumberPloidy` | `Focal_Amplification_Detection.md` |
| `GetChromosomeArmLengths` | `CopyNumberPloidy` | `Focal_Amplification_Detection.md` |
| `DetectHomozygousDeletions` / `IsHomozygousDeletion` (+ `purity` overloads) | `CopyNumberPloidy` | `Homozygous_Deletion_Detection.md` |
| `IdentifyDeletedTumorSuppressors` (arm labels; `GeneLocus` locus-overlap overload) | `CopyNumberPloidy` | `Homozygous_Deletion_Detection.md` |
| `CountCopyNumberStateOscillations` | `ClinicalMisc` | `Complex_Rearrangement_Classification.md` |

## Mutational signatures (SBS96)
| Method ID | file | Algorithm doc |
|---|---|---|
| `ClassifySbsContext` | `Signatures` | `SBS96_Trinucleotide_Context_Catalog.md` |
| `EnumerateSbs96Channels` | `Signatures` | `SBS96_Trinucleotide_Context_Catalog.md` |
| `FitSignatures` | `Signatures` | `Mutational_Signature_Fitting.md` |
| `ExtractSignatures` (2 overloads) | `Signatures` | `Mutational_Signature_Extraction_NMF.md` |
| `MatchToReferenceSignatures` | `Signatures` | `Mutational_Signature_Fitting.md` |
| `BootstrapExposures` | `Signatures` | `Mutational_Signature_Exposure_Bootstrap.md` |
| `GetMutationalProcess` / `ClassifyMutationalProcess` | `Signatures` | `Mutational_Process_Classification.md` |
| `CalculateSignatureScore` | `Signatures` | `Mutational_Process_Classification.md` |

## Gene fusions
| Method ID | file | Algorithm doc |
|---|---|---|
| `DetectFusions` | `Fusions` | `Fusion_Gene_Detection.md` |
| `IsInFrame` | `Fusions` | `Fusion_Gene_Detection.md` |
| `ComputeTotalSupport` | `Fusions` | `Fusion_Gene_Detection.md` |
| `GetFusionAnnotation` | `Fusions` | `Fusion_Gene_Detection.md` |
| `MatchKnownFusions` | `Fusions` | `Known_Fusion_Database_Lookup.md` |
| `AnalyzeBreakpoint` | `Fusions` | `Fusion_Breakpoint_Analysis.md` |
| `PredictFusionProtein` | `Fusions` | `Fusion_Breakpoint_Analysis.md` |
| `TestBreakpointClustering` | `ClinicalMisc` | `Complex_Rearrangement_Classification.md` |

## TMB / MSI
| Method ID | file | Algorithm doc |
|---|---|---|
| `CalculateTMB` (2 overloads) | `TmbMsi` | `Tumor_Mutational_Burden.md` |
| `ClassifyTMB` | `TmbMsi` | `Tumor_Mutational_Burden.md` |
| `CalculateMSIScore` | `TmbMsi` | `Microsatellite_Instability_Detection.md` |
| `ClassifyMSIStatus` | `TmbMsi` | `Microsatellite_Instability_Detection.md` |
| `ClassifyBethesdaPanel` | `TmbMsi` | `Microsatellite_Instability_Detection.md` |
| `DetectMSI` | `TmbMsi` | `Microsatellite_Instability_Detection.md` |

## Neoantigen / MHC binding
| Method ID | file | Notes |
|---|---|---|
| `GenerateNeoantigenPeptides` | `Neoantigen` | `Neoantigen_Peptide_Generation.md`; mutant peptide windows |
| `IsValidPeptideLength` | `Neoantigen` | length gate per `MhcClass` |
| `ClassifyMhcBinding` | `Neoantigen` | `MHC_Peptide_Binding_Classification.md`; classifies a **supplied** IC50 — ungated |
| `PredictBindingHalfLifeBimas` | `Neoantigen` | **GUARDED ONCO-MHC-001** (`LimitationPolicy.Enforce` in the method); caller supplies matrix |
| `PredictIc50Smm` | `Neoantigen` | **GUARDED ONCO-MHC-001** (`LimitationPolicy.Enforce` in the method); caller supplies matrix |
| `MhcflurryAffinityPredictor.LoadWeightPack` | `MhcflurryAffinityPredictor.cs` | caller supplies the ANN weight stream |
| `MhcflurryAffinityPredictor.PredictIc50` | `MhcflurryAffinityPredictor.cs` | ANN IC50 (nM) |
| `MhcflurryAffinityPredictor.PredictAndClassify` | `MhcflurryAffinityPredictor.cs` | IC50 + `BindingStrength` |

## Immune infiltration / deconvolution (GUARDED ONCO-IMMUNE-001)
| Method ID | file | Notes |
|---|---|---|
| `ImmuneAnalyzer.EstimateInfiltration` | `ImmuneAnalyzer.cs` | immune/stromal signature scores |
| `ImmuneAnalyzer.EstimateTumorPurity` | `ImmuneAnalyzer.cs` | **GUARDED** — `Enforce("ONCO-IMMUNE-001")`; ESTIMATE transform |
| `ImmuneAnalyzer.DeconvoluteImmuneCells` | `ImmuneAnalyzer.cs` | **GUARDED** |
| `ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr` | `ImmuneAnalyzer.cs` | **GUARDED** |
| `ImmuneAnalyzer.LoadSignatureMatrix` / `LoadBundledAbisSignatureMatrix` | `ImmuneAnalyzer.cs` | caller matrix / bundled ABIS (not exact LM22) |

## Clonal structure / CCF / phylogeny / heterogeneity
| Method ID | file | Algorithm doc |
|---|---|---|
| `ClassifyClonality` / `IdentifyClonalMutations` | `Clonality` | `Clonal_Subclonal_Classification.md` |
| `EstimateCcf` | `Clonality` | `Cancer_Cell_Fraction_Estimation.md` |
| `ClusterCcfValues` (fixed k; `(min, max)` BIC overload) | `Clonality` | `Cancer_Cell_Fraction_Estimation.md` |
| `ReconstructPhylogeny` / `TryReconstructPhylogeny` | `PhylogenyHeterogeneity` | `Tumor_Phylogeny_Reconstruction.md` |
| `ReconstructPhylogenyFromClusterSummaries` / `Try…` (+ `CcfClusterSummary.FromMembers`) | `PhylogenyHeterogeneity` | `Tumor_Phylogeny_Reconstruction.md` |
| `IdentifyTrunkMutations` / `IdentifyBranchMutations` | `PhylogenyHeterogeneity` | `Tumor_Phylogeny_Reconstruction.md` |
| `CalculateITH` (MATH; maftools-filter overload) / `CalculateCloneShannonDiversity` | `PhylogenyHeterogeneity` | `Tumor_Heterogeneity_Analysis.md` |
| `InferSubclones` / `AnalyzeHeterogeneity` | `PhylogenyHeterogeneity` | `Tumor_Heterogeneity_Analysis.md` |

## Drivers / HRD / LOH / HLA / ctDNA
| Method ID | file | Algorithm doc |
|---|---|---|
| `ClassifyGene` | `DriversArtifactsAnnotation` | `Driver_Mutation_Detection.md` |
| `ScoreDriverPotential` / `IdentifyDriverMutations` | `DriversArtifactsAnnotation` | `Driver_Mutation_Detection.md` |
| `DetectHRD` / `CalculateHrdLohScore` | `HrdLoh` | `HRD_Score.md` |
| `DetectLOH` / `CalculateLOHFraction` | `HrdLoh` | `Loss_Of_Heterozygosity.md` |
| `ParseHlaAllele` / `TryParseHlaAllele` / `DetectHlaLoh` | `ClinicalMisc` | `HLA_Nomenclature_And_Allele_Specific_LOH.md` |
| `IntegratedMutantAlleleFractionV2` | `CtdnaMrdChip` | `CtDNA_Analysis.md`, `MRD_Detection.md` |

## Recent additions (validation review B24, `docs/Validation/review-2026-09/B24.md`)
All additive (existing overloads keep their behaviour unless noted); F-numbers index the review's fix log.
- CNVkit purity overloads of `Log2RatioToCopyNumber` / `CallCopyNumber` / `ClassifyCopyNumber(s)` / `IsHomozygousDeletion` / `DetectHomozygousDeletions` (F29).
- `ComputeAscatGenomeMetrics` (F30); `EstimatePloidy(segments, probeCounts)` (F31); `DetectWholeGenomeDoubling(segments)` default denominator = facets-suite sample span, reference length via `(segments, ReferenceGenome)` (F32, deliberate default change).
- `AnalyzePurityPeaks` (F33/F34), `AnalyzeComplexKaryotypePeaks` / `AnalyzeSubclonalPurityPeaks` (F62), `FitBinomialMixture` + `PurityPeakOptions.FitMixturePeaks` / `BootstrapCount` (F63).
- `SegmentAlleleSpecific` now delegates to ASPCF (F35); germline-aware `SegmentAlleleSpecificAspcf(loci, germlineHeterozygous)` (F36) + `FitPurityPloidyFromAspcf` (F37); `AscatSexModel` overloads (F38); male X re-genotyping `AscatMaleXGenotyping` (F58); NA probes (F59).
- Battenberg: `FitSubclonalCopyNumberWithSnpTest` (F39), `SegmentPhasedBaf` / `BuildBattenbergSegments` (F40), `FitSubclonalCopyNumberWithBootstrap` (F41), `MergeBattenbergSegments` / `MaskHighCopyNumberSegments` (F60), `CallBattenbergSubclones` (F61).
- `SegmentAlleleSpecificAsMultiPcf` (ASCAT `ascat.asmultipcf`, F53/F59).
- LICHeE: `ReconstructPhylogenyFromClusterSummaries` + `CcfClusterSummary.FromMembers` (F42/F43); default ε `DefaultPhylogenyTolerance` = 0.1 (F44).
- `ClusterCcfValues(ccf, minClusters, maxClusters)` — BIC-selected k (F47).
- GISTIC2: `GeneLocus` locus-overlap gene mapping (F49), `DeconstructZiggurat` / `DetectFocalAmplificationEvents` (F50–F52, F54), `GetChromosomeArmLengths` (F55).
- `CalculateITH(vafs, vafCutOff, minMutations = 5)` maftools filters (F56); `CalculateCloneShannonDiversity` (F57).

## Notes for callers
- **Namespace is `Seqeron.Genomics.Oncology`** (single, unlike the split `Core`/`Alignment`/… layout
  documented in [`seqeron-dev`](../../seqeron-dev/SKILL.md)). `MhcClass` / `BindingStrength` /
  `PurityPloidyFit` / `AscatSexModel` / `CcfClusterSummary` etc. are nested types on `OncologyAnalyzer` — see source.
- **Two guarded units, both `MinimumMode = Moderate`** (run by default; throw only under `Strict`):
  ONCO-MHC-001 (`PredictIc50Smm`, `PredictBindingHalfLifeBimas`) and ONCO-IMMUNE-001
  (`ImmuneAnalyzer.EstimateTumorPurity` / `Deconvolute*`). ONCO-ASCAT-001 / ONCO-PURITY-001 have no
  runtime guard and no `LIMITATIONS.md` entry (their caveats live in the algorithm docs). Full policy
  mechanics + `SeqeronLimitationException` fields: [`seqeron-dev`](../../seqeron-dev/SKILL.md);
  envelope + STOP rule: [`bio-rigor`](../../bio-rigor/SKILL.md).
- **Symbols not exhaustively verified here:** the nested result/record types (e.g. `PurityPloidyFit`,
  `FusionCall`, `SignatureFitResult`) and their properties were not each grepped — read the method's
  XML doc in source before relying on a specific property name.
