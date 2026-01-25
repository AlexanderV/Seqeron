# Полный аудит методов SuffixTree.Genomics

## Итого: 277 public static методов

### По файлам:

| # | Файл | Метод | Строка |
|---|------|-------|--------|
| 1 | ApproximateMatcher.cs | HammingDistance | 163 |
| 2 | ApproximateMatcher.cs | EditDistance | 186 |
| 3 | ApproximateMatcher.cs | CountApproximateOccurrences | 283 |
| 4 | BedParser.cs | GetTotalBlockLength | 469 |
| 5 | BedParser.cs | CalculateStatistics | 510 |
| 6 | BedParser.cs | WriteToFile | 549 |
| 7 | BedParser.cs | WriteToStream | 558 |
| 8 | BedParser.cs | ExtractSequence | 655 |
| 9 | CancellableOperations.cs | GlobalAlign | 355 |
| 10 | ChromosomeAnalyzer.cs | AnalyzeKaryotype | 136 |
| 11 | ChromosomeAnalyzer.cs | AnalyzeTelomeres | 235 |
| 12 | ChromosomeAnalyzer.cs | EstimateTelomereLengthFromTSRatio | 330 |
| 13 | ChromosomeAnalyzer.cs | AnalyzeCentromere | 346 |
| 14 | ChromosomeAnalyzer.cs | CalculateArmRatio | 844 |
| 15 | ChromosomeAnalyzer.cs | ClassifyChromosomeByArmRatio | 858 |
| 16 | ChromosomeAnalyzer.cs | EstimateCellDivisionsFromTelomereLength | 876 |
| 17 | CodonOptimizer.cs | OptimizeSequence | 213 |
| 18 | CodonOptimizer.cs | CalculateCAI | 388 |
| 19 | CodonOptimizer.cs | RemoveRestrictionSites | 436 |
| 20 | CodonOptimizer.cs | ReduceSecondaryStructure | 489 |
| 21 | CodonOptimizer.cs | CompareCodonUsage | 616 |
| 22 | CodonOptimizer.cs | CreateCodonTableFromSequence | 673 |
| 23 | CodonUsageAnalyzer.cs | CalculateCai (DnaSequence) | 120 |
| 24 | CodonUsageAnalyzer.cs | CalculateCai (string) | 131 |
| 25 | CodonUsageAnalyzer.cs | CalculateEnc (DnaSequence) | 324 |
| 26 | CodonUsageAnalyzer.cs | CalculateEnc (string) | 333 |
| 27 | CodonUsageAnalyzer.cs | GetStatistics (DnaSequence) | 410 |
| 28 | CodonUsageAnalyzer.cs | GetStatistics (string) | 419 |
| 29 | ComparativeGenomics.cs | CompareGenomes | 424 |
| 30 | ComparativeGenomics.cs | CalculateReversalDistance | 468 |
| 31 | ComparativeGenomics.cs | CalculateANI | 585 |
| 32 | CrisprDesigner.cs | GetSystem | 18 |
| 33 | CrisprDesigner.cs | EvaluateGuideRna | 200 |
| 34 | CrisprDesigner.cs | CalculateSpecificityScore | 371 |
| 35 | DisorderPredictor.cs | PredictDisorder | 130 |
| 36 | DisorderPredictor.cs | GetDisorderPropensity | 558 |
| 37 | DisorderPredictor.cs | IsDisorderPromoting | 566 |
| 38 | DnaSequence.cs | TryCreate | 129 |
| 39 | DnaSequence.cs | GetReverseComplementString | 149 |
| 40 | DnaSequence.cs | TryGetReverseComplement | 200 |
| 41 | EmblParser.cs | ParseLocation | 471 |
| 42 | EmblParser.cs | ExtractSequence | 600 |
| 43 | EpigeneticsAnalyzer.cs | CalculateCpGObservedExpected | 191 |
| 44 | EpigeneticsAnalyzer.cs | SimulateBisulfiteConversion | 283 |
| 45 | EpigeneticsAnalyzer.cs | GenerateMethylationProfile | 385 |
| 46 | EpigeneticsAnalyzer.cs | PredictChromatinState | 519 |
| 47 | EpigeneticsAnalyzer.cs | CalculateEpigeneticAge | 721 |
| 48 | FastaParser.cs | ToFasta | 78 |
| 49 | FastaParser.cs | WriteFile | 98 |
| 50 | FastqParser.cs | DetectEncoding | 148 |
| 51 | FastqParser.cs | EncodeQualityScores | 189 |
| 52 | FastqParser.cs | PhredToErrorProbability | 206 |
| 53 | FastqParser.cs | ErrorProbabilityToPhred | 214 |
| 54 | FastqParser.cs | TrimByQuality | 264 |
| 55 | FastqParser.cs | TrimAdapter | 297 |
| 56 | FastqParser.cs | CalculateStatistics | 344 |
| 57 | FastqParser.cs | WriteToFile | 437 |
| 58 | FastqParser.cs | WriteToStream | 446 |
| 59 | FastqParser.cs | ToFastqString | 474 |
| 60 | FeatureLocationHelper.cs | ExtractSequence (GenBank) | 60 |
| 61 | FeatureLocationHelper.cs | ExtractSequence (EMBL) | 66 |
| 62 | GcSkewCalculator.cs | CalculateGcSkew (DnaSequence) | 21 |
| 63 | GcSkewCalculator.cs | CalculateGcSkew (string) | 30 |
| 64 | GcSkewCalculator.cs | CalculateAtSkew (DnaSequence) | 167 |
| 65 | GcSkewCalculator.cs | CalculateAtSkew (string) | 176 |
| 66 | GcSkewCalculator.cs | PredictReplicationOrigin | 203 |
| 67 | GcSkewCalculator.cs | AnalyzeGcContent | 238 |
| 68 | GenBankParser.cs | ParseLocation | 471 |
| 69 | GenBankParser.cs | ExtractSequence | 556 |
| 70 | GeneticCode.cs | GetByTableNumber | 134 |
| 71 | GenomeAnnotator.cs | CalculateCodingPotential | 460 |
| 72 | GenomeAssemblyAnalyzer.cs | CalculateStatistics | 121 |
| 73 | GenomeAssemblyAnalyzer.cs | CalculateNx | 222 |
| 74 | GenomeAssemblyAnalyzer.cs | CalculateAuN | 270 |
| 75 | GenomeAssemblyAnalyzer.cs | AssessCompleteness | 469 |
| 76 | GenomeAssemblyAnalyzer.cs | CompareAssemblies | 791 |
| 77 | GenomicAnalyzer.cs | FindLongestRepeat | 20 |
| 78 | GenomicAnalyzer.cs | FindLongestCommonRegion | 178 |
| 79 | GenomicAnalyzer.cs | CalculateSimilarity | 238 |
| 80 | GffParser.cs | CalculateStatistics | 377 |
| 81 | GffParser.cs | WriteToFile | 417 |
| 82 | GffParser.cs | WriteToStream | 426 |
| 83 | GffParser.cs | ExtractSequence | 488 |
| 84 | ISequence.cs | GetIupacCode | 228 |
| 85 | ISequence.cs | CodesMatch | 273 |
| 86 | ISequence.cs | PhredToErrorProbability | 484 |
| 87 | ISequence.cs | ErrorProbabilityToPhred | 490 |
| 88 | IupacHelper.cs | MatchesIupac | 14 |
| 89 | KmerAnalyzer.cs | KmerDistance | 165 |
| 90 | KmerAnalyzer.cs | CalculateKmerEntropy | 243 |
| 91 | KmerAnalyzer.cs | AnalyzeKmers | 363 |
| 92 | MetagenomicsAnalyzer.cs | GenerateTaxonomicProfile | 229 |
| 93 | MetagenomicsAnalyzer.cs | CalculateAlphaDiversity | 288 |
| 94 | MetagenomicsAnalyzer.cs | CalculateBetaDiversity | 353 |
| 95 | MiRnaAnalyzer.cs | GetSeedSequence | 85 |
| 96 | MiRnaAnalyzer.cs | CreateMiRna | 96 |
| 97 | MiRnaAnalyzer.cs | GetReverseComplement | 226 |
| 98 | MiRnaAnalyzer.cs | CanPair | 254 |
| 99 | MiRnaAnalyzer.cs | IsWobblePair | 267 |
| 100 | MiRnaAnalyzer.cs | AlignMiRnaToTarget | 282 |
| 101 | MiRnaAnalyzer.cs | CalculateSiteAccessibility | 536 |
| 102 | MiRnaAnalyzer.cs | CalculateGcContent | 615 |
| 103 | MotifFinder.cs | CreatePwm | 145 |
| 104 | MotifFinder.cs | GenerateConsensus | 257 |
| 105 | PanGenomeAnalyzer.cs | ConstructPanGenome | 83 |
| 106 | PanGenomeAnalyzer.cs | FitHeapsLaw | 378 |
| 107 | PanGenomeAnalyzer.cs | CreateCoreGenomeAlignment | 487 |
| 108 | PhylogeneticAnalyzer.cs | BuildTree | 77 |
| 109 | PhylogeneticAnalyzer.cs | CalculatePairwiseDistance | 134 |
| 110 | PhylogeneticAnalyzer.cs | ToNewick | 430 |
| 111 | PhylogeneticAnalyzer.cs | ParseNewick | 469 |
| 112 | PhylogeneticAnalyzer.cs | CalculateTreeLength | 589 |
| 113 | PhylogeneticAnalyzer.cs | GetTreeDepth | 603 |
| 114 | PhylogeneticAnalyzer.cs | RobinsonFouldsDistance | 616 |
| 115 | PhylogeneticAnalyzer.cs | PatristicDistance | 682 |
| 116 | PopulationGeneticsAnalyzer.cs | CalculateMAF | 127 |
| 117 | PopulationGeneticsAnalyzer.cs | CalculateNucleotideDiversity | 168 |
| 118 | PopulationGeneticsAnalyzer.cs | CalculateWattersonTheta | 202 |
| 119 | PopulationGeneticsAnalyzer.cs | CalculateTajimasD | 220 |
| 120 | PopulationGeneticsAnalyzer.cs | CalculateDiversityStatistics | 261 |
| 121 | PopulationGeneticsAnalyzer.cs | TestHardyWeinberg | 352 |
| 122 | PopulationGeneticsAnalyzer.cs | CalculateFst | 453 |
| 123 | PopulationGeneticsAnalyzer.cs | CalculateFStatistics | 512 |
| 124 | PopulationGeneticsAnalyzer.cs | CalculateLD | 560 |
| 125 | PopulationGeneticsAnalyzer.cs | CalculateIHS | 684 |
| 126 | PopulationGeneticsAnalyzer.cs | CalculateInbreedingFromROH | 870 |
| 127 | PrimerDesigner.cs | DesignPrimers | 39 |
| 128 | PrimerDesigner.cs | EvaluatePrimer | 119 |
| 129 | PrimerDesigner.cs | CalculateMeltingTemperature | 193 |
| 130 | PrimerDesigner.cs | CalculateMeltingTemperatureWithSalt | 219 |
| 131 | PrimerDesigner.cs | CalculateGcContent | 229 |
| 132 | PrimerDesigner.cs | FindLongestHomopolymer | 235 |
| 133 | PrimerDesigner.cs | FindLongestDinucleotideRepeat | 262 |
| 134 | PrimerDesigner.cs | HasHairpinPotential | 296 |
| 135 | PrimerDesigner.cs | HasPrimerDimer | 387 |
| 136 | PrimerDesigner.cs | Calculate3PrimeStability | 414 |
| 137 | ProbeDesigner.cs | DesignTilingProbes | 411 |
| 138 | ProbeDesigner.cs | ValidateProbe | 552 |
| 139 | ProbeDesigner.cs | CheckSpecificity | 604 |
| 140 | ProbeDesigner.cs | CalculateMolecularWeight | 647 |
| 141 | ProbeDesigner.cs | CalculateExtinctionCoefficient | 673 |
| 142 | ProbeDesigner.cs | CalculateConcentration | 699 |
| 143 | ProteinMotifFinder.cs | ConvertPrositeToRegex | 201 |
| 144 | ProteinSequence.cs | TryCreate | 357 |
| 145 | QualityScoreAnalyzer.cs | CharToPhred | 66 |
| 146 | QualityScoreAnalyzer.cs | PhredToChar | 75 |
| 147 | QualityScoreAnalyzer.cs | PhredToQualityString | 99 |
| 148 | QualityScoreAnalyzer.cs | DetectEncoding | 107 |
| 149 | QualityScoreAnalyzer.cs | PhredToErrorProbability | 128 |
| 150 | QualityScoreAnalyzer.cs | ErrorProbabilityToPhred | 136 |
| 151 | QualityScoreAnalyzer.cs | CalculateStatistics (1) | 146 |
| 152 | QualityScoreAnalyzer.cs | CalculateStatistics (2) | 173 |
| 153 | QualityScoreAnalyzer.cs | QualityTrim | 254 |
| 154 | QualityScoreAnalyzer.cs | SlidingWindowTrim | 308 |
| 155 | QualityScoreAnalyzer.cs | CalculateExpectedErrors | 402 |
| 156 | QualityScoreAnalyzer.cs | MaskLowQualityBases | 413 |
| 157 | QualityScoreAnalyzer.cs | CalculateConsensusQuality | 569 |
| 158 | RepeatFinder.cs | GetTandemRepeatSummary | 358 |
| 159 | ReportGenerator.cs | CreateBuilder | 77 |
| 160 | ReportGenerator.cs | Generate | 174 |
| 161 | ReportGenerator.cs | SaveToFile | 190 |
| 162 | ReportGenerator.cs | CreateSequenceAnalysisReport | 817 |
| 163 | ReportGenerator.cs | CreateComparisonReport | 854 |
| 164 | RestrictionAnalyzer.cs | GetDigestSummary | 301 |
| 165 | RestrictionAnalyzer.cs | CreateMap | 322 |
| 166 | RestrictionAnalyzer.cs | AreCompatible | 375 |
| 167 | RnaSecondaryStructure.cs | CanPair | 142 |
| 168 | RnaSecondaryStructure.cs | GetComplement | 174 |
| 169 | RnaSecondaryStructure.cs | CalculateStemEnergy | 286 |
| 170 | RnaSecondaryStructure.cs | CalculateHairpinLoopEnergy | 309 |
| 171 | RnaSecondaryStructure.cs | CalculateMinimumFreeEnergy | 365 |
| 172 | RnaSecondaryStructure.cs | PredictStructure | 428 |
| 173 | RnaSecondaryStructure.cs | ValidateDotBracket | 620 |
| 174 | RnaSecondaryStructure.cs | CalculateStructureProbability | 684 |
| 175 | RnaSecondaryStructure.cs | GenerateRandomRna | 698 |
| 176 | RnaSequence.cs | FromDna | 147 |
| 177 | RnaSequence.cs | TryCreate | 176 |
| 178 | SequenceAligner.cs | GlobalAlign (1) | 55 |
| 179 | SequenceAligner.cs | GlobalAlign (2) | 69 |
| 180 | SequenceAligner.cs | GlobalAlign (3) | 93 |
| 181 | SequenceAligner.cs | GlobalAlign (4) | 106 |
| 182 | SequenceAligner.cs | LocalAlign (1) | 164 |
| 183 | SequenceAligner.cs | LocalAlign (2) | 178 |
| 184 | SequenceAligner.cs | SemiGlobalAlign | 286 |
| 185 | SequenceAligner.cs | CalculateStatistics | 406 |
| 186 | SequenceAligner.cs | FormatAlignment | 446 |
| 187 | SequenceAligner.cs | MultipleAlign | 495 |
| 188 | SequenceAssembler.cs | AssembleOLC | 48 |
| 189 | SequenceAssembler.cs | AssembleDeBruijn | 75 |
| 190 | SequenceAssembler.cs | CalculateIdentity | 175 |
| 191 | SequenceAssembler.cs | CalculateStats | 388 |
| 192 | SequenceAssembler.cs | MergeContigs | 423 |
| 193 | SequenceAssembler.cs | ComputeConsensus | 534 |
| 194 | SequenceComplexity.cs | CalculateLinguisticComplexity (DnaSequence) | 22 |
| 195 | SequenceComplexity.cs | CalculateLinguisticComplexity (string) | 33 |
| 196 | SequenceComplexity.cs | CalculateShannonEntropy (DnaSequence) | 78 |
| 197 | SequenceComplexity.cs | CalculateShannonEntropy (string) | 87 |
| 198 | SequenceComplexity.cs | CalculateKmerEntropy | 128 |
| 199 | SequenceComplexity.cs | CalculateDustScore (DnaSequence) | 296 |
| 200 | SequenceComplexity.cs | CalculateDustScore (string) | 305 |
| 201 | SequenceComplexity.cs | MaskLowComplexity | 346 |
| 202 | SequenceComplexity.cs | EstimateCompressionRatio (DnaSequence) | 391 |
| 203 | SequenceComplexity.cs | EstimateCompressionRatio (string) | 400 |
| 204 | SequenceExtensions.cs | CalculateGcContent | 21 |
| 205 | SequenceExtensions.cs | CalculateGcContentFast | 41 |
| 206 | SequenceExtensions.cs | CalculateGcFraction | 50 |
| 207 | SequenceExtensions.cs | CalculateGcFractionFast | 69 |
| 208 | SequenceExtensions.cs | GetComplementBase | 83 |
| 209 | SequenceExtensions.cs | GetRnaComplementBase | 98 |
| 210 | SequenceExtensions.cs | TryGetComplement | 115 |
| 211 | SequenceExtensions.cs | TryGetReverseComplement | 131 |
| 212 | SequenceExtensions.cs | EnumerateKmers | 174 |
| 213 | SequenceExtensions.cs | HammingDistance | 187 |
| 214 | SequenceExtensions.cs | IsValidDna | 210 |
| 215 | SequenceExtensions.cs | IsValidRna | 225 |
| 216 | SequenceIO.cs | ToGenBank | 328 |
| 217 | SequenceIO.cs | ToBed | 545 |
| 218 | SequenceIO.cs | ToGff | 662 |
| 219 | SequenceIO.cs | IsPaired | 749 |
| 220 | SequenceIO.cs | IsProperPair | 750 |
| 221 | SequenceIO.cs | IsUnmapped | 751 |
| 222 | SequenceIO.cs | IsMateUnmapped | 752 |
| 223 | SequenceIO.cs | IsReverse | 753 |
| 224 | SequenceIO.cs | IsMateReverse | 754 |
| 225 | SequenceIO.cs | IsRead1 | 755 |
| 226 | SequenceIO.cs | IsRead2 | 756 |
| 227 | SequenceIO.cs | IsSecondary | 757 |
| 228 | SequenceIO.cs | IsQcFail | 758 |
| 229 | SequenceIO.cs | IsDuplicate | 759 |
| 230 | SequenceIO.cs | IsSupplementary | 760 |
| 231 | SequenceIO.cs | ToPhylip | 913 |
| 232 | SequenceStatistics.cs | CalculateNucleotideComposition | 48 |
| 233 | SequenceStatistics.cs | CalculateAminoAcidComposition | 98 |
| 234 | SequenceStatistics.cs | CalculateMolecularWeight | 159 |
| 235 | SequenceStatistics.cs | CalculateNucleotideMolecularWeight | 180 |
| 236 | SequenceStatistics.cs | CalculateIsoelectricPoint | 228 |
| 237 | SequenceStatistics.cs | CalculateHydrophobicity | 306 |
| 238 | SequenceStatistics.cs | CalculateThermodynamics | 381 |
| 239 | SequenceStatistics.cs | CalculateMeltingTemperature | 441 |
| 240 | SequenceStatistics.cs | CalculateShannonEntropy | 580 |
| 241 | SequenceStatistics.cs | CalculateLinguisticComplexity | 615 |
| 242 | SequenceStatistics.cs | SummarizeNucleotideSequence | 775 |
| 243 | SpliceSitePredictor.cs | PredictGeneStructure | 478 |
| 244 | SpliceSitePredictor.cs | CalculateMaxEntScore | 713 |
| 245 | SpliceSitePredictor.cs | IsWithinCodingRegion | 756 |
| 246 | StatisticsHelper.cs | NormalCDF | 13 |
| 247 | StatisticsHelper.cs | Erf | 21 |
| 248 | ThermoConstants.cs | CalculateWallaceTm | 87 |
| 249 | ThermoConstants.cs | CalculateMarmurDotyTm | 96 |
| 250 | ThermoConstants.cs | CalculateSaltAdjustedTm | 109 |
| 251 | ThermoConstants.cs | CalculateSaltCorrection | 121 |
| 252 | TranscriptomeAnalyzer.cs | CalculateEnrichmentScore | 392 |
| 253 | TranscriptomeAnalyzer.cs | CalculatePearsonCorrelation | 573 |
| 254 | Translator.cs | Translate (DnaSequence) | 20 |
| 255 | Translator.cs | Translate (RnaSequence) | 37 |
| 256 | Translator.cs | Translate (string) | 54 |
| 257 | VariantAnnotator.cs | ClassifyVariant | 194 |
| 258 | VariantAnnotator.cs | NormalizeVariant | 226 |
| 259 | VariantAnnotator.cs | GetImpactLevel | 507 |
| 260 | VariantAnnotator.cs | PredictPathogenicity | 728 |
| 261 | VariantAnnotator.cs | ParseVcfVariant | 1193 |
| 262 | VariantAnnotator.cs | FormatAsVcfInfo | 1208 |
| 263 | VariantCaller.cs | ClassifyMutation | 178 |
| 264 | VariantCaller.cs | CalculateTiTvRatio | 197 |
| 265 | VariantCaller.cs | CalculateStatistics | 223 |
| 266 | VariantCaller.cs | PredictEffect | 263 |
| 267 | VcfParser.cs | ClassifyVariant | 376 |
| 268 | VcfParser.cs | IsSNP | 407 |
| 269 | VcfParser.cs | IsIndel | 412 |
| 270 | VcfParser.cs | GetVariantLength | 421 |
| 271 | VcfParser.cs | IsHomRef | 529 |
| 272 | VcfParser.cs | IsHomAlt | 539 |
| 273 | VcfParser.cs | IsHet | 554 |
| 274 | VcfParser.cs | CalculateStatistics | 631 |
| 275 | VcfParser.cs | WriteToFile | 702 |
| 276 | VcfParser.cs | WriteToStream | 715 |
| 277 | VcfParser.cs | HasInfoFlag | 835 |

---

## Сводка по классам

| Класс | Количество методов |
|-------|-------------------|
| ApproximateMatcher | 3 |
| BedParser | 5 |
| CancellableOperations | 1 |
| ChromosomeAnalyzer | 7 |
| CodonOptimizer | 6 |
| CodonUsageAnalyzer | 6 |
| ComparativeGenomics | 3 |
| CrisprDesigner | 3 |
| DisorderPredictor | 3 |
| DnaSequence | 3 |
| EmblParser | 2 |
| EpigeneticsAnalyzer | 5 |
| FastaParser | 2 |
| FastqParser | 10 |
| FeatureLocationHelper | 2 |
| GcSkewCalculator | 6 |
| GenBankParser | 2 |
| GeneticCode | 1 |
| GenomeAnnotator | 1 |
| GenomeAssemblyAnalyzer | 5 |
| GenomicAnalyzer | 3 |
| GffParser | 4 |
| ISequence | 4 |
| IupacHelper | 1 |
| KmerAnalyzer | 3 |
| MetagenomicsAnalyzer | 3 |
| MiRnaAnalyzer | 8 |
| MotifFinder | 2 |
| PanGenomeAnalyzer | 3 |
| PhylogeneticAnalyzer | 8 |
| PopulationGeneticsAnalyzer | 11 |
| PrimerDesigner | 10 |
| ProbeDesigner | 6 |
| ProteinMotifFinder | 1 |
| ProteinSequence | 1 |
| QualityScoreAnalyzer | 13 |
| RepeatFinder | 1 |
| ReportGenerator | 5 |
| RestrictionAnalyzer | 3 |
| RnaSecondaryStructure | 9 |
| RnaSequence | 2 |
| SequenceAligner | 10 |
| SequenceAssembler | 6 |
| SequenceComplexity | 10 |
| SequenceExtensions | 12 |
| SequenceIO | 16 |
| SequenceStatistics | 11 |
| SpliceSitePredictor | 3 |
| StatisticsHelper | 2 |
| ThermoConstants | 4 |
| TranscriptomeAnalyzer | 2 |
| Translator | 3 |
| VariantAnnotator | 6 |
| VariantCaller | 4 |
| VcfParser | 11 |
| **ИТОГО** | **277** |

---

## Примечания

1. **Перегрузки** считаются отдельно (например, GlobalAlign имеет 4 перегрузки = 4 метода)
2. **Extension methods** из SequenceExtensions включены (12 методов)
3. **StatisticsHelper** - внутренний helper, может быть исключен из MCP (2 метода)
4. **SAM flags** в SequenceIO - 12 однострочных методов для работы с SAM флагами

## Дата аудита: 2026-01-23
