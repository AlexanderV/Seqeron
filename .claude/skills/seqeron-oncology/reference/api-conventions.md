# seqeron-oncology — curated Method-ID map (C#-API only)

There is **no MCP oncology server** and Oncology is **not** in `docs/skills/domain-map.json`, so there
is **no generated tool slice** — this hand-curated index is the map. Every entry is a real
`Seqeron.Genomics.Oncology` symbol verified against source; paths are from repo root `../../../../`.
Confirm anything with
`grep -n "public static" src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`.

All three classes are `static` and share namespace **`Seqeron.Genomics.Oncology`**:
`OncologyAnalyzer` (`OncologyAnalyzer.cs:30`), `ImmuneAnalyzer` (`ImmuneAnalyzer.cs:23`),
`MhcflurryAffinityPredictor` (`MhcflurryAffinityPredictor.cs:45`). Below, unqualified members are on
`OncologyAnalyzer`; `src/.../` = `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/`.

## Purity / ploidy / allele-specific copy number
| Method ID | file:line | Algorithm doc |
|---|---|---|
| `FitPurityPloidy` | `OncologyAnalyzer.cs:7411` | `Tumor_Purity_Estimation.md`, `Tumor_Ploidy_Estimation.md`, `Allele_Specific_Copy_Number_Derivation.md` |
| `EstimatePurity` | `:515` | `Tumor_Purity_Estimation.md` |
| `EstimatePurityFromVAF` | `:451` | `Tumor_Purity_Estimation.md` |
| `EstimatePurityFromVaf` (scalar) | `:480` | `Tumor_Purity_Estimation.md` |
| `AdjustVAFForPurity` | `:414` | `Variant_Allele_Frequency.md` |
| `EstimatePloidy` | `:6554` | `Tumor_Ploidy_Estimation.md` |
| `SegmentAlleleSpecific` | `:7282` | `Allele_Specific_Copy_Number_Derivation.md` |
| `SegmentAlleleSpecificAspcf` | `:7641` | `Allele_Specific_Copy_Number_Derivation.md` |
| `DeriveMultiplicity` | `:7579` | `Allele_Specific_Copy_Number_Derivation.md` |
| `FitSubclonalCopyNumber` | `:7830` | `Copy_Number_Alteration_Classification.md` |
| `DetectWholeGenomeDoublingFromSuppliedLength` | `:6659` | `Tumor_Ploidy_Estimation.md` |

## Copy-number calling / classification
| Method ID | file:line | Algorithm doc |
|---|---|---|
| `Log2RatioToCopyNumber` | `:5900` | `Copy_Number_Alteration_Classification.md` |
| `CallCopyNumber` | `:5928` | `Copy_Number_Alteration_Classification.md` |
| `ClassifyCopyNumber` / `ClassifyCopyNumbers` | `:5971` / `:5995` | `Copy_Number_Alteration_Classification.md` |
| `DetectFocalAmplifications` | `:6157` | `Focal_Amplification_Detection.md` |
| `DetectHomozygousDeletions` | `:6297` | `Homozygous_Deletion_Detection.md` |
| `CountCopyNumberStateOscillations` | `:11823` | `Complex_Rearrangement_Classification.md` |

## Mutational signatures (SBS96)
| Method ID | file:line | Algorithm doc |
|---|---|---|
| `ClassifySbsContext` | `:2659` | `SBS96_Trinucleotide_Context_Catalog.md` |
| `EnumerateSbs96Channels` | `:2696` | `SBS96_Trinucleotide_Context_Catalog.md` |
| `FitSignatures` | `:2933` | `Mutational_Signature_Fitting.md` |
| `ExtractSignatures` (2 overloads) | `:3431`, `:3485` | `Mutational_Signature_Extraction_NMF.md` |
| `MatchToReferenceSignatures` | `:4609` | `Mutational_Signature_Fitting.md` |
| `BootstrapExposures` | `:4774` | `Mutational_Signature_Exposure_Bootstrap.md` |
| `GetMutationalProcess` / `ClassifyMutationalProcess` | `:5121` / `:5189` | `Mutational_Process_Classification.md` |
| `CalculateSignatureScore` | `:12130` | `Mutational_Process_Classification.md` |

## Gene fusions
| Method ID | file:line | Algorithm doc |
|---|---|---|
| `DetectFusions` | `:5421` | `Fusion_Gene_Detection.md` |
| `IsInFrame` | `:5382` | `Fusion_Gene_Detection.md` |
| `ComputeTotalSupport` | `:5365` | `Fusion_Gene_Detection.md` |
| `GetFusionAnnotation` | `:5521` | `Fusion_Gene_Detection.md` |
| `MatchKnownFusions` | `:5558` | `Known_Fusion_Database_Lookup.md` |
| `AnalyzeBreakpoint` | `:5703` | `Fusion_Breakpoint_Analysis.md` |
| `PredictFusionProtein` | `:5750` | `Fusion_Breakpoint_Analysis.md` |
| `TestBreakpointClustering` | `:11850` | `Complex_Rearrangement_Classification.md` |

## TMB / MSI
| Method ID | file:line | Algorithm doc |
|---|---|---|
| `CalculateTMB` (2 overloads) | `:1508`, `:1540` | `Tumor_Mutational_Burden.md` |
| `ClassifyTMB` | `:1557` | `Tumor_Mutational_Burden.md` |
| `CalculateMSIScore` | `:1634` | `Microsatellite_Instability_Detection.md` |
| `ClassifyMSIStatus` | `:1664` | `Microsatellite_Instability_Detection.md` |
| `ClassifyBethesdaPanel` | `:1688` | `Microsatellite_Instability_Detection.md` |
| `DetectMSI` | `:1725` | `Microsatellite_Instability_Detection.md` |

## Neoantigen / MHC binding
| Method ID | file:line | Notes |
|---|---|---|
| `GenerateNeoantigenPeptides` | `:8026` | `Neoantigen_Peptide_Generation.md`; mutant peptide windows |
| `IsValidPeptideLength` | `:8273` | length gate per `MhcClass` |
| `ClassifyMhcBinding` | `:8292` | `MHC_Peptide_Binding_Classification.md`; classifies a **supplied** IC50 — ungated |
| `PredictBindingHalfLifeBimas` | `:8483` | **GUARDED ONCO-MHC-001** — `Enforce` at `:8490`; caller supplies matrix |
| `PredictIc50Smm` | `:8516` | **GUARDED ONCO-MHC-001** — `Enforce` at `:8523`; caller supplies matrix |
| `MhcflurryAffinityPredictor.LoadWeightPack` | `MhcflurryAffinityPredictor.cs:516` | caller supplies the ANN weight stream |
| `MhcflurryAffinityPredictor.PredictIc50` | `MhcflurryAffinityPredictor.cs:619` | ANN IC50 (nM) |
| `MhcflurryAffinityPredictor.PredictAndClassify` | `MhcflurryAffinityPredictor.cs:673` | IC50 + `BindingStrength` |

## Immune infiltration / deconvolution (GUARDED ONCO-IMMUNE-001)
| Method ID | file:line | Notes |
|---|---|---|
| `ImmuneAnalyzer.EstimateInfiltration` | `ImmuneAnalyzer.cs:402` | immune/stromal signature scores |
| `ImmuneAnalyzer.EstimateTumorPurity` | `ImmuneAnalyzer.cs:473` | **GUARDED** — `Enforce("ONCO-IMMUNE-001")` at `:477`; ESTIMATE transform |
| `ImmuneAnalyzer.DeconvoluteImmuneCells` | `ImmuneAnalyzer.cs:516` | **GUARDED** — `Enforce` at `:525` |
| `ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr` | `ImmuneAnalyzer.cs:658` | **GUARDED** — `Enforce` at `:680` |
| `ImmuneAnalyzer.LoadSignatureMatrix` / `LoadBundledAbisSignatureMatrix` | `ImmuneAnalyzer.cs:813` / `:898` | caller matrix / bundled ABIS (not exact LM22) |

## Clonal structure / phylogeny / heterogeneity
| Method ID | file:line | Algorithm doc |
|---|---|---|
| `ClassifyClonality` | `:6789` | `Clonal_Subclonal_Classification.md` |
| `EstimateCcf` | `:7014` | `Cancer_Cell_Fraction_Estimation.md` |
| `ClusterCcfValues` | `:7061` | `Tumor_Heterogeneity_Analysis.md` |
| `InferSubclones` | `:11009` | `Tumor_Heterogeneity_Analysis.md` |
| `AnalyzeHeterogeneity` | `:11042` | `Tumor_Heterogeneity_Analysis.md` |
| `ReconstructPhylogeny` | `:10617` | `Tumor_Phylogeny_Reconstruction.md` |
| `IdentifyTrunkMutations` / `IdentifyBranchMutations` | `:10727` / `:10760` | `Tumor_Phylogeny_Reconstruction.md` |

## Drivers / HRD / LOH / HLA / ctDNA
| Method ID | file:line | Algorithm doc |
|---|---|---|
| `ClassifyGene` | `:738` | `Driver_Mutation_Detection.md` |
| `ScoreDriverPotential` / `IdentifyDriverMutations` | `:794` / `:830` | `Driver_Mutation_Detection.md` |
| `DetectHRD` / `CalculateHrdLohScore` | `:1896` / `:2017` | `HRD_Score.md` |
| `DetectLOH` / `CalculateLOHFraction` | `:1979` / `:2033` | `Loss_Of_Heterozygosity.md` |
| `ParseHlaAllele` / `TryParseHlaAllele` / `DetectHlaLoh` | `:11281` / `:11350` / `:11384` | `HLA_Nomenclature_And_Allele_Specific_LOH.md` |
| `IntegratedMutantAlleleFractionV2` | `:9087` | `CtDNA_Analysis.md`, `MRD_Detection.md` |

## Notes for callers
- **Namespace is `Seqeron.Genomics.Oncology`** (single, unlike the split `Core`/`Alignment`/… layout
  documented in [`seqeron-dev`](../../seqeron-dev/SKILL.md)). `MhcClass` / `BindingStrength` /
  `PurityPloidyFit` etc. are nested types on `OncologyAnalyzer` — see source.
- **Two guarded units, both `MinimumMode = Moderate`** (run by default; throw only under `Strict`):
  ONCO-MHC-001 (`PredictIc50Smm`, `PredictBindingHalfLifeBimas`) and ONCO-IMMUNE-001
  (`ImmuneAnalyzer.EstimateTumorPurity` / `Deconvolute*`). ONCO-ASCAT-001 / ONCO-PURITY-001 are
  documented-only (no runtime guard). Full policy mechanics + `SeqeronLimitationException` fields:
  [`seqeron-dev`](../../seqeron-dev/SKILL.md); envelope + STOP rule: [`bio-rigor`](../../bio-rigor/SKILL.md).
- **Symbols not exhaustively verified here:** the nested result/record types (e.g. `PurityPloidyFit`,
  `FusionCall`, `SignatureFitResult`) and their properties were not each grepped — read the source
  region around the listed line before relying on a specific property name.
