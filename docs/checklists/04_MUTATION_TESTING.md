# Checklist 04: Mutation Testing (Stryker)

**Priority:** P1  
**Framework:** Stryker.NET  
**Date:** 2026-03-19  
**Total algorithms:** 258

---

## Description

Mutation testing вимірює ефективність тестового набору: Stryker вносить мутації у вихідний код (заміна операторів, видалення умов тощо) і перевіряє, що тести їх виявляють («вбивають»). Мутанти, що пережили всі тести (survived), вказують на прогалини в покритті. Не вимагає написання нового коду для запуску — лише конфігурацію та аналіз.

**Поточне покриття:** `MutationKillerTests.cs` таргетує MotifFinder.cs та RepeatFinder.cs. `stryker-config.json` існує. Результати в `StrykerOutput/`.

**Процес для кожного модуля:**
1. Запустити `dotnet stryker` для вихідного файлу
2. Проаналізувати survived мутантів
3. Написати прицільні killer-тести для survivors зі score < 80%

**Target:** mutation score ≥ 80% для кожного модуля.

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | Source File | Test File(s) | Target Score |
|---|--------|-----------|------|------------|-------------|:---:|
| 1 | ☑ | SEQ-GC-001 | Composition | SequenceExtensions.cs | SequenceExtensions_CalculateGcContent_Tests.cs | ≥ 80% |
| 2 | ☑ | SEQ-COMP-001 | Composition | SequenceExtensions.cs | SequenceExtensions_Complement_Tests.cs | ≥ 80% |
| 3 | ☑ | SEQ-REVCOMP-001 | Composition | SequenceExtensions.cs | SequenceExtensions_ReverseComplement_Tests.cs | ≥ 80% |
| 4 | ☑ | SEQ-VALID-001 | Composition | SequenceExtensions.cs | SequenceExtensions_SequenceValidation_Tests.cs | ≥ 80% |
| 5 | ☑ | SEQ-COMPLEX-001 | Composition | SequenceComplexity.cs | SequenceComplexityTests.cs | ≥ 80% |
| 6 | ☑ | SEQ-ENTROPY-001 | Composition | SequenceComplexity.cs | SequenceComplexityTests.cs | ≥ 80% |
| 7 | ☑ | SEQ-GCSKEW-001 | Composition | GcSkewCalculator.cs | GcSkewCalculatorTests.cs | ≥ 80% |
| 8 | ☑ | PAT-EXACT-001 | Matching | MotifFinder.cs | FindAllOccurrencesTests.cs, ContainsTests.cs, CountOccurrencesTests.cs | ≥ 80% |
| 9 | ☑ | PAT-APPROX-001 | Matching | ApproximateMatcher.cs | ApproximateMatcher_HammingDistance_Tests.cs | ≥ 80% |
| 10 | ☑ | PAT-APPROX-002 | Matching | ApproximateMatcher.cs | ApproximateMatcher_EditDistance_Tests.cs | ≥ 80% |
| 11 | ☑ | PAT-IUPAC-001 | Matching | MotifFinder.cs, IupacHelper.cs | IupacMotifMatchingTests.cs, MutationKillerTests.cs | ≥ 80% |
| 12 | ☑ | PAT-PWM-001 | Matching | MotifFinder.cs | MotifFinder_PWM_Tests.cs, MutationKillerTests.cs | ≥ 80% |
| 13 | ☑ | REP-STR-001 | Repeats | RepeatFinder.cs | RepeatFinder_Microsatellite_Tests.cs | ≥ 80% |
| 14 | ☑ | REP-TANDEM-001 | Repeats | RepeatFinder.cs, GenomicAnalyzer.cs | GenomicAnalyzer_TandemRepeat_Tests.cs, RepeatFinderTests.cs | ≥ 80% |
| 15 | ☑ | REP-INV-001 | Repeats | RepeatFinder.cs | RepeatFinder_InvertedRepeat_Tests.cs | ≥ 80% |
| 16 | ☑ | REP-DIRECT-001 | Repeats | RepeatFinder.cs | RepeatFinder_DirectRepeat_Tests.cs | ≥ 80% |
| 17 | ☑ | REP-PALIN-001 | Repeats | RepeatFinder.cs | RepeatFinder_Palindrome_Tests.cs | ≥ 80% |
| 18 | ☑ | CRISPR-PAM-001 | MolTools | CrisprDesigner.cs | CrisprDesigner_PAM_Tests.cs | ≥ 80% |
| 19 | ☑ | CRISPR-GUIDE-001 | MolTools | CrisprDesigner.cs | CrisprDesigner_GuideRNA_Tests.cs | ≥ 80% |
| 20 | ☑ | CRISPR-OFF-001 | MolTools | CrisprDesigner.cs | CrisprDesigner_OffTarget_Tests.cs | ≥ 80% |
| 21 | ☑ | PRIMER-TM-001 | MolTools | PrimerDesigner.cs | PrimerDesigner_MeltingTemperature_Tests.cs | ≥ 80% |
| 22 | ☑ | PRIMER-DESIGN-001 | MolTools | PrimerDesigner.cs | PrimerDesigner_PrimerDesign_Tests.cs | ≥ 80% |
| 23 | ☑ | PRIMER-STRUCT-001 | MolTools | PrimerDesigner.cs | PrimerDesigner_PrimerStructure_Tests.cs | ≥ 80% |
| 24 | ☑ | PROBE-DESIGN-001 | MolTools | ProbeDesigner.cs | ProbeDesigner_ProbeDesign_Tests.cs | ≥ 80% |
| 25 | ☑ | PROBE-VALID-001 | MolTools | ProbeDesigner.cs | ProbeDesigner_ProbeValidation_Tests.cs | ≥ 80% |
| 26 | ☑ | RESTR-FIND-001 | MolTools | RestrictionAnalyzer.cs | RestrictionAnalyzer_FindSites_Tests.cs | ≥ 80% |
| 27 | ☑ | RESTR-DIGEST-001 | MolTools | RestrictionAnalyzer.cs | RestrictionAnalyzer_Digest_Tests.cs | ≥ 80% |
| 28 | ☑ | ANNOT-ORF-001 | Annotation | GenomeAnnotator.cs | GenomeAnnotator_ORF_Tests.cs | ≥ 80% |
| 29 | ☑ | ANNOT-GENE-001 | Annotation | GenomeAnnotator.cs | GenomeAnnotator_Gene_Tests.cs | ≥ 80% |
| 30 | ☑ | ANNOT-PROM-001 | Annotation | GenomeAnnotator.cs | GenomeAnnotator_PromoterMotif_Tests.cs | ≥ 80% |
| 31 | ☑ | ANNOT-GFF-001 | Annotation | GenomeAnnotator.cs | GenomeAnnotator_GFF3_Tests.cs | ≥ 80% |
| 32 | ☑ | KMER-COUNT-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_CountKmers_Tests.cs | ≥ 80% |
| 33 | ☑ | KMER-FREQ-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_Frequency_Tests.cs | ≥ 80% |
| 34 | ☑ | KMER-FIND-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_Find_Tests.cs | ≥ 80% |
| 35 | ☑ | ALIGN-GLOBAL-001 | Alignment | SequenceAligner.cs | SequenceAligner_GlobalAlign_Tests.cs | ≥ 80% |
| 36 | ☑ | ALIGN-LOCAL-001 | Alignment | SequenceAligner.cs | SequenceAligner_LocalAlign_Tests.cs | ≥ 80% |
| 37 | ☑ | ALIGN-SEMI-001 | Alignment | SequenceAligner.cs | SequenceAligner_SemiGlobalAlign_Tests.cs | ≥ 80% |
| 38 | ☑ | ALIGN-MULTI-001 | Alignment | SequenceAligner.cs | SequenceAligner_MultipleAlign_Tests.cs | ≥ 80% |
| 39 | ☑ | PHYLO-DIST-001 | Phylogenetic | PhylogeneticAnalyzer.cs | PhylogeneticAnalyzer_DistanceMatrix_Tests.cs | ≥ 80% |
| 40 | ☑ | PHYLO-TREE-001 | Phylogenetic | PhylogeneticAnalyzer.cs | PhylogeneticAnalyzer_TreeConstruction_Tests.cs | ≥ 80% |
| 41 | ☑ | PHYLO-NEWICK-001 | Phylogenetic | PhylogeneticAnalyzer.cs | PhylogeneticAnalyzer_NewickIO_Tests.cs | ≥ 80% |
| 42 | ☑ | PHYLO-COMP-001 | Phylogenetic | PhylogeneticAnalyzer.cs | PhylogeneticAnalyzer_TreeComparison_Tests.cs | ≥ 80% |
| 43 | ☑ | POP-FREQ-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs | ≥ 80% |
| 44 | ☑ | POP-DIV-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_Diversity_Tests.cs | ≥ 80% |
| 45 | ☑ | POP-HW-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs | ≥ 80% |
| 46 | ☑ | POP-FST-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_FStatistics_Tests.cs | ≥ 80% |
| 47 | ☑ | POP-LD-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs | ≥ 80% |
| 48 | ☑ | CHROM-TELO-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_Telomere_Tests.cs | ≥ 80% |
| 49 | ☑ | CHROM-CENT-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_Centromere_Tests.cs | ≥ 80% |
| 50 | ☑ | CHROM-KARYO-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_Karyotype_Tests.cs | ≥ 80% |
| 51 | ☑ | CHROM-ANEU-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_Aneuploidy_Tests.cs | ≥ 80% |
| 52 | ☑ | CHROM-SYNT-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_Synteny_Tests.cs | ≥ 80% |
| 53 | ☑ | META-CLASS-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs | ≥ 80% |
| 54 | ☑ | META-PROF-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_TaxonomicProfile_Tests.cs | ≥ 80% |
| 55 | ☑ | META-ALPHA-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_AlphaDiversity_Tests.cs | ≥ 80% |
| 56 | ☑ | META-BETA-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_BetaDiversity_Tests.cs | ≥ 80% |
| 57 | ☑ | META-BIN-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_GenomeBinning_Tests.cs | ≥ 80% |
| 58 | ☑ | CODON-OPT-001 | Codon | CodonOptimizer.cs | CodonOptimizer_OptimizeSequence_Tests.cs | ≥ 80% |
| 59 | ☑ | CODON-CAI-001 | Codon | CodonOptimizer.cs | CodonOptimizer_CAI_Tests.cs | ≥ 80% |
| 60 | ☑ | CODON-RARE-001 | Codon | CodonOptimizer.cs | CodonOptimizer_FindRareCodons_Tests.cs | ≥ 80% |
| 61 | ☑ | CODON-USAGE-001 | Codon | CodonOptimizer.cs | CodonOptimizer_CodonUsage_Tests.cs | ≥ 80% |
| 62 | ☑ | TRANS-CODON-001 | Translation | GeneticCode.cs | GeneticCodeTests.cs | ≥ 80% |
| 63 | ☑ | TRANS-PROT-001 | Translation | Translator.cs | TranslatorTests.cs | ≥ 80% |
| 64 | ☑ | PARSE-FASTA-001 | FileIO | FastaParser.cs | FastaParserTests.cs | ≥ 80% |
| 65 | ☑ | PARSE-FASTQ-001 | FileIO | FastqParser.cs | FastqParserTests.cs | ≥ 80% |
| 66 | ☑ | PARSE-BED-001 | FileIO | BedParser.cs | BedParserTests.cs | ≥ 80% |
| 67 | ☑ | PARSE-VCF-001 | FileIO | VcfParser.cs | VcfParserTests.cs | ≥ 80% |
| 68 | ☑ | PARSE-GFF-001 | FileIO | GffParser.cs | GffParserTests.cs | ≥ 80% |
| 69 | ☑ | PARSE-GENBANK-001 | FileIO | GenBankParser.cs | GenBankParserTests.cs | ≥ 80% |
| 70 | ☑ | PARSE-EMBL-001 | FileIO | EmblParser.cs | EmblParserTests.cs | ≥ 80% |
| 71 | ☑ | RNA-STRUCT-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructureTests.cs | ≥ 80% |
| 72 | ☑ | RNA-STEMLOOP-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructureTests.cs | ≥ 80% |
| 73 | ☑ | RNA-ENERGY-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructureTests.cs | ≥ 80% |
| 74 | ☐ | MIRNA-SEED-001 | MiRNA | MiRnaAnalyzer.cs | MiRnaAnalyzer_SeedAnalysis_Tests.cs | ≥ 80% |
| 75 | ☐ | MIRNA-TARGET-001 | MiRNA | MiRnaAnalyzer.cs | MiRnaAnalyzer_TargetPrediction_Tests.cs | ≥ 80% |
| 76 | ☐ | MIRNA-PRECURSOR-001 | MiRNA | MiRnaAnalyzer.cs | MiRnaAnalyzer_PreMiRna_Tests.cs | ≥ 80% |
| 77 | ☑ | SPLICE-DONOR-001 | Splicing | SpliceSitePredictor.cs | SpliceSitePredictor_DonorSite_Tests.cs | ≥ 80% |
| 78 | ☑ | SPLICE-ACCEPTOR-001 | Splicing | SpliceSitePredictor.cs | SpliceSitePredictor_AcceptorSite_Tests.cs | ≥ 80% |
| 79 | ☑ | SPLICE-PREDICT-001 | Splicing | SpliceSitePredictor.cs | SpliceSitePredictor_GeneStructure_Tests.cs | ≥ 80% |
| 80 | ☑ | DISORDER-PRED-001 | ProteinPred | DisorderPredictor.cs | DisorderPredictor_DisorderPrediction_Tests.cs | ≥ 80% |
| 81 | ☐ | DISORDER-REGION-001 | ProteinPred | DisorderPredictor.cs | DisorderPredictor_DisorderedRegion_Tests.cs | ≥ 80% |
| 82 | ☑ | PROTMOTIF-FIND-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_MotifSearch_Tests.cs | ≥ 80% |
| 83 | ☑ | PROTMOTIF-PROSITE-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_PrositePattern_Tests.cs | ≥ 80% |
| 84 | ☑ | PROTMOTIF-DOMAIN-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_DomainPrediction_Tests.cs | ≥ 80% |
| 85 | ☑ | EPIGEN-CPG-001 | Epigenetics | EpigeneticsAnalyzer.cs | EpigeneticsAnalyzer_CpGDetection_Tests.cs | ≥ 80% |
| 86 | ☐ | ONCO-IMMUNE-001 | Oncology | ImmuneAnalyzer.cs | ImmuneAnalyzer_ImmuneInfiltration_Tests.cs | ≥ 80% |
| 87 | ☑ | ONCO-SOMATIC-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_CallSomaticMutations_Tests.cs | ≥ 80% |
| 88 | ☑ | ONCO-VAF-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_CalculateVAF_Tests.cs | ≥ 80% |
| 89 | ☑ | ONCO-DRIVER-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_IdentifyDriverMutations_Tests.cs | ≥ 80% |
| 90 | ☑ | ONCO-ARTIFACT-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_FilterArtifacts_Tests.cs | ≥ 80% |
| 91 | ☑ | ONCO-ANNOT-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_AnnotateCancerVariants_Tests.cs | ≥ 80% |
| 92 | ☑ | ONCO-TMB-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_CalculateTMB_Tests.cs | ≥ 80% |
| 93 | ☑ | ONCO-MSI-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_DetectMSI_Tests.cs | ≥ 80% |
| 94 | ☐ | ONCO-HRD-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_CalculateHRDScore_Tests.cs | ≥ 80% |
| 95 | ☐ | ONCO-LOH-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_DetectLOH_Tests.cs | ≥ 80% |
| 96 | ☑ | ONCO-SIG-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_ClassifySbsContext_Tests.cs | ≥ 80% |
| 97 | ☐ | ONCO-SIG-002 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_FitSignatures_Tests.cs | ≥ 80% |
| 98 | ☐ | ONCO-SIG-003 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_BootstrapExposures_Tests.cs | ≥ 80% |
| 99 | ☑ | ONCO-SIG-004 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_ClassifyMutationalProcess_Tests.cs | ≥ 80% |
| 100 | ☑ | ONCO-FUSION-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_DetectFusions_Tests.cs | ≥ 80% |
| 101 | ☑ | ONCO-FUSION-002 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_MatchKnownFusions_Tests.cs | ≥ 80% |
| 102 | ☑ | ONCO-FUSION-003 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_AnalyzeBreakpoint_Tests.cs | ≥ 80% |
| 103 | ☑ | ONCO-CNA-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_CopyNumberClassification_Tests.cs | ≥ 80% |
| 104 | ☑ | ONCO-CNA-002 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_DetectFocalAmplifications_Tests.cs | ≥ 80% |
| 105 | ☑ | ONCO-CNA-003 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_DetectHomozygousDeletions_Tests.cs | ≥ 80% |
| 106 | ☐ | ONCO-PURITY-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_EstimatePurity_Tests.cs | ≥ 80% |
| 107 | ☐ | ONCO-PLOIDY-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_EstimatePloidy_Tests.cs | ≥ 80% |
| 108 | ☑ | ONCO-CLONAL-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_ClassifyClonality_Tests.cs | ≥ 80% |
| 109 | ☑ | ONCO-NEO-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_GenerateNeoantigenPeptides_Tests.cs | ≥ 80% |
| 110 | ☐ | ONCO-MHC-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_ClassifyMhcBinding_Tests.cs | ≥ 80% |
| 111 | ☑ | ONCO-CTDNA-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_CtDnaAnalysis_Tests.cs | ≥ 80% |
| 112 | ☐ | ONCO-MRD-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_DetectMRD_Tests.cs | ≥ 80% |
| 113 | ☐ | ONCO-CHIP-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_FilterCHIP_Tests.cs | ≥ 80% |
| 114 | ☑ | ONCO-PHYLO-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_ReconstructPhylogeny_Tests.cs | ≥ 80% |
| 115 | ☑ | ONCO-CCF-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_EstimateCcf_Tests.cs | ≥ 80% |
| 116 | ☑ | ONCO-HETERO-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_AnalyzeHeterogeneity_Tests.cs | ≥ 80% |
| 117 | ☑ | ONCO-HLA-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_HlaAnalysis_Tests.cs | ≥ 80% |
| 118 | ☑ | ONCO-ACTION-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_AssessActionability_Tests.cs | ≥ 80% |
| 119 | ☑ | ONCO-SV-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_ClassifyComplexRearrangement_Tests.cs | ≥ 80% |
| 120 | ☑ | ONCO-EXPR-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_IdentifyOutlierGenes_Tests.cs | ≥ 80% |
| 121 | ☑ | SEQ-COMPOSITION-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_CalculateNucleotideComposition_Tests.cs | ≥ 80% |
| 122 | ☑ | SEQ-DINUC-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_CalculateDinucleotide_Tests.cs | ≥ 80% |
| 123 | ☑ | SEQ-HYDRO-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_CalculateHydrophobicity_Tests.cs | ≥ 80% |
| 124 | ☑ | SEQ-MW-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_CalculateMolecularWeight_Tests.cs | ≥ 80% |
| 125 | ☑ | SEQ-PI-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_CalculateIsoelectricPoint_Tests.cs | ≥ 80% |
| 126 | ☑ | SEQ-SECSTRUCT-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_PredictSecondaryStructure_Tests.cs | ≥ 80% |
| 127 | ☐ | SEQ-STATS-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_CalculateNucleotideComposition_Tests.cs | ≥ 80% |
| 128 | ☑ | SEQ-SUMMARY-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_SummarizeNucleotideSequence_Tests.cs | ≥ 80% |
| 129 | ☑ | SEQ-THERMO-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_CalculateThermodynamics_Tests.cs | ≥ 80% |
| 130 | ☑ | SEQ-TM-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_CalculateThermodynamics_Tests.cs | ≥ 80% |
| 131 | ☐ | COMPGEN-ANI-001 | Comparative | ComparativeGenomics.cs | ComparativeGenomics_CalculateANI_Tests.cs | ≥ 80% |
| 132 | ☑ | COMPGEN-CLUSTER-001 | Comparative | ComparativeGenomics.cs | ComparativeGenomics_FindConservedClusters_Tests.cs | ≥ 80% |
| 133 | ☑ | COMPGEN-COMPARE-001 | Comparative | ComparativeGenomics.cs | ComparativeGenomics_CompareGenomes_Tests.cs | ≥ 80% |
| 134 | ☑ | COMPGEN-DOTPLOT-001 | Comparative | ComparativeGenomics.cs | ComparativeGenomics_GenerateDotPlot_Tests.cs | ≥ 80% |
| 135 | ☑ | COMPGEN-ORTHO-001 | Comparative | ComparativeGenomics.cs | ComparativeGenomics_FindOrthologs_Tests.cs | ≥ 80% |
| 136 | ☑ | COMPGEN-RBH-001 | Comparative | ComparativeGenomics.cs | ComparativeGenomics_FindReciprocalBestHits_Tests.cs | ≥ 80% |
| 137 | ☑ | COMPGEN-REARR-001 | Comparative | ComparativeGenomics.cs | ComparativeGenomics_DetectRearrangements_Tests.cs | ≥ 80% |
| 138 | ☑ | COMPGEN-REVERSAL-001 | Comparative | ComparativeGenomics.cs | ComparativeGenomics_CalculateReversalDistance_Tests.cs | ≥ 80% |
| 139 | ☑ | COMPGEN-SYNTENY-001 | Comparative | ComparativeGenomics.cs | ComparativeGenomics_FindSyntenicBlocks_Tests.cs | ≥ 80% |
| 140 | ☑ | ASSEMBLY-CONSENSUS-001 | Assembly | SequenceAssembler.cs | SequenceAssembler_ComputeConsensus_Tests.cs | ≥ 80% |
| 141 | ☑ | ASSEMBLY-CORRECT-001 | Assembly | SequenceAssembler.cs | SequenceAssembler_ErrorCorrectReads_Tests.cs | ≥ 80% |
| 142 | ☑ | ASSEMBLY-COVER-001 | Assembly | SequenceAssembler.cs | SequenceAssembler_CalculateCoverage_Tests.cs | ≥ 80% |
| 143 | ☑ | ASSEMBLY-DBG-001 | Assembly | SequenceAssembler.cs | SequenceAssembler_AssembleDeBruijn_Tests.cs | ≥ 80% |
| 144 | ☑ | ASSEMBLY-MERGE-001 | Assembly | SequenceAssembler.cs | SequenceAssembler_MergeContigs_Tests.cs | ≥ 80% |
| 145 | ☑ | ASSEMBLY-OLC-001 | Assembly | SequenceAssembler.cs | SequenceAssembler_AssembleOLC_Tests.cs | ≥ 80% |
| 146 | ☑ | ASSEMBLY-SCAFFOLD-001 | Assembly | SequenceAssembler.cs | SequenceAssembler_Scaffold_Tests.cs | ≥ 80% |
| 147 | ☑ | ASSEMBLY-STATS-001 | Assembly | GenomeAssemblyAnalyzer.cs | GenomeAssemblyAnalyzer_AssemblyStatistics_Tests.cs | ≥ 80% |
| 148 | ☑ | ASSEMBLY-TRIM-001 | Assembly | SequenceAssembler.cs | SequenceAssembler_QualityTrimReads_Tests.cs | ≥ 80% |
| 149 | ☑ | RNA-DOTBRACKET-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructure_ParseDotBracket_Tests.cs | ≥ 80% |
| 150 | ☑ | RNA-HAIRPIN-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructure_HairpinEnergy_Tests.cs | ≥ 80% |
| 151 | ☑ | RNA-INVERT-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructure_FindInvertedRepeats_Tests.cs | ≥ 80% |
| 152 | ☑ | RNA-MFE-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructure_MinimumFreeEnergy_Tests.cs | ≥ 80% |
| 153 | ☑ | RNA-PAIR-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructure_CanPair_Tests.cs | ≥ 80% |
| 154 | ☑ | RNA-PARTITION-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructure_PartitionFunction_Tests.cs | ≥ 80% |
| 155 | ☑ | RNA-PSEUDOKNOT-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructure_DetectPseudoknots_Tests.cs | ≥ 80% |
| 156 | ☑ | KMER-ASYNC-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_CountKmersAsync_Tests.cs | ≥ 80% |
| 157 | ☑ | KMER-BOTH-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_CountKmersBothStrands_Tests.cs | ≥ 80% |
| 158 | ☑ | KMER-DIST-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_KmerDistance_Tests.cs | ≥ 80% |
| 159 | ☑ | KMER-GENERATE-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_GenerateAllKmers_Tests.cs | ≥ 80% |
| 160 | ☑ | KMER-POSITIONS-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_FindKmerPositions_Tests.cs | ≥ 80% |
| 161 | ☑ | KMER-STATS-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_AnalyzeKmers_Tests.cs | ≥ 80% |
| 162 | ☑ | KMER-UNIQUE-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_FindUniqueAndMinCount_Tests.cs | ≥ 80% |
| 163 | ☑ | PROTMOTIF-CC-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_PredictCoiledCoils_Tests.cs | ≥ 80% |
| 164 | ☑ | PROTMOTIF-COMMON-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_FindCommonMotifs_Tests.cs | ≥ 80% |
| 165 | ☑ | PROTMOTIF-LC-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_FindLowComplexityRegions_Tests.cs | ≥ 80% |
| 166 | ☑ | PROTMOTIF-PATTERN-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_FindMotifByPattern_Tests.cs | ≥ 80% |
| 167 | ☑ | PROTMOTIF-SP-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_PredictSignalPeptide_Tests.cs | ≥ 80% |
| 168 | ☑ | PROTMOTIF-TM-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_PredictTransmembraneHelices_Tests.cs | ≥ 80% |
| 169 | ☑ | MOTIF-CONS-001 | Matching | MotifFinder.cs | MotifFinder_CreateConsensusFromAlignment_Tests.cs | ≥ 80% |
| 170 | ☑ | MOTIF-DISCOVER-001 | Matching | MotifFinder.cs | MotifFinder_DiscoverMotifs_Tests.cs | ≥ 80% |
| 171 | ☑ | MOTIF-GENERATE-001 | Matching | MotifFinder.cs | MotifFinder_GenerateConsensus_Tests.cs | ≥ 80% |
| 172 | ☑ | MOTIF-REGULATORY-001 | Matching | MotifFinder.cs | MotifFinder_FindRegulatoryElements_Tests.cs | ≥ 80% |
| 173 | ☑ | MOTIF-SHARED-001 | Matching | MotifFinder.cs | MotifFinder_FindSharedMotifs_Tests.cs | ≥ 80% |
| 174 | ☑ | PAT-APPROX-003 | Matching | ApproximateMatcher.cs | ApproximateMatcher_FindBestMatch_Tests.cs | ≥ 80% |
| 175 | ☑ | GENOMIC-COMMON-001 | Analysis | GenomicAnalyzer.cs | GenomicAnalyzer_FindCommonRegion_Tests.cs | ≥ 80% |
| 176 | ☑ | GENOMIC-MOTIFS-001 | Analysis | GenomicAnalyzer.cs | GenomicAnalyzer_FindKnownMotifs_Tests.cs | ≥ 80% |
| 177 | ☑ | GENOMIC-ORF-001 | Analysis | GenomicAnalyzer.cs | GenomicAnalyzer_FindOpenReadingFrames_Tests.cs | ≥ 80% |
| 178 | ☑ | GENOMIC-REPEAT-001 | Analysis | GenomicAnalyzer.cs | GenomicAnalyzer_FindRepeats_Tests.cs | ≥ 80% |
| 179 | ☑ | GENOMIC-SIMILARITY-001 | Analysis | GenomicAnalyzer.cs | GenomicAnalyzer_CalculateSimilarity_Tests.cs | ≥ 80% |
| 180 | ☑ | GENOMIC-TANDEM-001 | Analysis | GenomicAnalyzer.cs | GenomicAnalyzer_TandemRepeat_Tests.cs | ≥ 80% |
| 181 | ☐ | EPIGEN-AGE-001 | Epigenetics | EpigeneticsAnalyzer.cs | EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests.cs | ≥ 80% |
| 182 | ☑ | EPIGEN-BISULF-001 | Epigenetics | EpigeneticsAnalyzer.cs | EpigeneticsAnalyzer_Bisulfite_Tests.cs | ≥ 80% |
| 183 | ☑ | EPIGEN-CHROM-001 | Epigenetics | EpigeneticsAnalyzer.cs | EpigeneticsAnalyzer_ChromatinState_Tests.cs | ≥ 80% |
| 184 | ☑ | EPIGEN-DMR-001 | Epigenetics | EpigeneticsAnalyzer.cs | EpigeneticsAnalyzer_DMR_Tests.cs | ≥ 80% |
| 185 | ☑ | EPIGEN-METHYL-001 | Epigenetics | EpigeneticsAnalyzer.cs | EpigeneticsAnalyzer_Methylation_Tests.cs | ≥ 80% |
| 186 | ☑ | VARIANT-ANNOT-001 | Variants | VariantAnnotator.cs | VariantAnnotator_FunctionalImpact_Tests.cs | ≥ 80% |
| 187 | ☐ | VARIANT-CALL-001 | Variants | VariantCaller.cs | VariantCaller_CallVariants_Tests.cs | ≥ 80% |
| 188 | ☑ | VARIANT-INDEL-001 | Variants | VariantCaller.cs | VariantCaller_FindIndels_Tests.cs | ≥ 80% |
| 189 | ☑ | VARIANT-SNP-001 | Variants | VariantCaller.cs | VariantCaller_FindSnps_Tests.cs | ≥ 80% |
| 190 | ☑ | PANGEN-CLUSTER-001 | PanGenome | PanGenomeAnalyzer.cs | PanGenomeAnalyzer_ClusterGenes_Tests.cs | ≥ 80% |
| 191 | ☑ | PANGEN-CORE-001 | PanGenome | PanGenomeAnalyzer.cs | PanGenomeAnalyzer_ConstructPanGenome_Tests.cs | ≥ 80% |
| 192 | ☑ | PANGEN-HEAP-001 | PanGenome | PanGenomeAnalyzer.cs | PanGenomeAnalyzer_FitHeapsLaw_Tests.cs | ≥ 80% |
| 193 | ☑ | PANGEN-MARKER-001 | PanGenome | PanGenomeAnalyzer.cs | PanGenomeAnalyzer_SelectPhylogeneticMarkers_Tests.cs | ≥ 80% |
| 194 | ☑ | META-FUNC-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_PredictFunctions_Tests.cs | ≥ 80% |
| 195 | ☑ | META-PATHWAY-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_FindPathwayEnrichment_Tests.cs | ≥ 80% |
| 196 | ☑ | META-RESIST-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_FindAntibioticResistanceGenes_Tests.cs | ≥ 80% |
| 197 | ☑ | META-TAXA-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_FindSignificantTaxa_Tests.cs | ≥ 80% |
| 198 | ☑ | TRANS-DIFF-001 | Transcriptome | TranscriptomeAnalyzer.cs | TranscriptomeAnalyzer_DifferentialExpression_Tests.cs | ≥ 80% |
| 199 | ☑ | TRANS-EXPR-001 | Transcriptome | TranscriptomeAnalyzer.cs | TranscriptomeAnalyzer_ExpressionQuantification_Tests.cs | ≥ 80% |
| 200 | ☑ | TRANS-SPLICE-001 | Transcriptome | TranscriptomeAnalyzer.cs | TranscriptomeAnalyzer_AlternativeSplicing_Tests.cs | ≥ 80% |
| 201 | ☑ | SV-BREAKPOINT-001 | StructuralVar | StructuralVariantAnalyzer.cs | StructuralVariantAnalyzer_FindBreakpoints_Tests.cs | ≥ 80% |
| 202 | ☑ | SV-CNV-001 | StructuralVar | StructuralVariantAnalyzer.cs | StructuralVariantAnalyzer_DetectCNV_Tests.cs | ≥ 80% |
| 203 | ☑ | SV-DETECT-001 | StructuralVar | StructuralVariantAnalyzer.cs | StructuralVariantAnalyzer_DetectSVs_Tests.cs | ≥ 80% |
| 204 | ☑ | DISORDER-LC-001 | ProteinPred | DisorderPredictor.cs | DisorderPredictor_LowComplexity_Tests.cs | ≥ 80% |
| 205 | ☑ | DISORDER-MORF-001 | ProteinPred | DisorderPredictor.cs | DisorderPredictor_MoRF_Tests.cs | ≥ 80% |
| 206 | ☑ | DISORDER-PROPENSITY-001 | ProteinPred | DisorderPredictor.cs | DisorderPredictor_GetDisorderPropensity_Tests.cs | ≥ 80% |
| 207 | ☑ | POP-ANCESTRY-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_EstimateAncestry_Tests.cs | ≥ 80% |
| 208 | ☑ | POP-ROH-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_FindROH_Tests.cs | ≥ 80% |
| 209 | ☑ | POP-SELECT-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_SelectionSignature_Tests.cs | ≥ 80% |
| 210 | ☑ | SEQ-ATSKEW-001 | Composition | GcSkewCalculator.cs | GcSkewCalculator_CalculateAtSkew_Tests.cs | ≥ 80% |
| 211 | ☑ | SEQ-REPLICATION-001 | Composition | GcSkewCalculator.cs | GcSkewCalculator_PredictReplicationOrigin_Tests.cs | ≥ 80% |
| 212 | ☑ | SEQ-RNACOMP-001 | Composition | SequenceExtensions.cs | SequenceExtensions_GetRnaComplementBase_Tests.cs | ≥ 80% |
| 213 | ☑ | CODON-ENC-001 | Codon | CodonUsageAnalyzer.cs | CodonUsageAnalyzer_CalculateEnc_Tests.cs | ≥ 80% |
| 214 | ☑ | CODON-RSCU-001 | Codon | CodonUsageAnalyzer.cs | CodonUsageAnalyzer_CalculateRscu_Tests.cs | ≥ 80% |
| 215 | ☑ | CODON-STATS-001 | Codon | CodonUsageAnalyzer.cs | CodonUsageAnalyzer_GetStatistics_Tests.cs | ≥ 80% |
| 216 | ☑ | ANNOT-CODING-001 | Annotation | GenomeAnnotator.cs | GenomeAnnotator_CalculateCodingPotential_Tests.cs | ≥ 80% |
| 217 | ☑ | ANNOT-CODONUSAGE-001 | Annotation | GenomeAnnotator.cs | GenomeAnnotator_GetCodonUsage_Tests.cs | ≥ 80% |
| 218 | ☑ | ANNOT-REPEAT-001 | Annotation | GenomeAnnotator.cs | GenomeAnnotator_FindRepetitiveElements_Tests.cs | ≥ 80% |
| 219 | ☑ | QUALITY-PHRED-001 | Quality | QualityScoreAnalyzer.cs | QualityScoreAnalyzer_ParseQualityString_Tests.cs | ≥ 80% |
| 220 | ☑ | QUALITY-STATS-001 | Quality | QualityScoreAnalyzer.cs | QualityScoreAnalyzer_CalculateStatistics_Tests.cs | ≥ 80% |
| 221 | ☑ | PHYLO-BOOT-001 | Phylogenetic | PhylogeneticAnalyzer.cs | PhylogeneticAnalyzer_Bootstrap_Tests.cs | ≥ 80% |
| 222 | ☑ | PHYLO-STATS-001 | Phylogenetic | PhylogeneticAnalyzer.cs | PhylogeneticAnalyzer_TreeStatistics_Tests.cs | ≥ 80% |
| 223 | ☑ | TRANS-SIXFRAME-001 | Translation | Translator.cs | Translator_SixFrames_Tests.cs | ≥ 80% |
| 224 | ☑ | RESTR-FILTER-001 | MolTools | RestrictionAnalyzer.cs | RestrictionAnalyzer_Filter_Tests.cs | ≥ 80% |
| 225 | ☑ | MIRNA-PAIR-001 | MiRNA | MiRnaAnalyzer.cs | MiRnaAnalyzer_AlignMiRnaToTarget_Tests.cs | ≥ 80% |
| 226 | ☑ | ALIGN-STATS-001 | Alignment | SequenceAligner.cs | SequenceAligner_CalculateStatistics_Tests.cs | ≥ 80% |
| 227 | ☑ | SEQ-CODON-FREQ-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_CalculateCodonFrequencies_Tests.cs | ≥ 80% |
| 228 | ☑ | SEQ-COMPLEX-COMPRESS-001 | Complexity | SequenceComplexity.cs | SequenceComplexity_EstimateCompressionRatio_Tests.cs | ≥ 80% |
| 229 | ☑ | SEQ-COMPLEX-DUST-001 | Complexity | SequenceComplexity.cs | SequenceComplexity_CalculateDustScore_Tests.cs | ≥ 80% |
| 230 | ☑ | SEQ-COMPLEX-KMER-001 | Complexity | SequenceComplexity.cs | SequenceComplexity_CalculateKmerEntropy_Tests.cs | ≥ 80% |
| 231 | ☑ | SEQ-COMPLEX-WINDOW-001 | Complexity | SequenceComplexity.cs | SequenceComplexity_CalculateWindowedComplexity_Tests.cs | ≥ 80% |
| 232 | ☑ | SEQ-ENTROPY-PROFILE-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_CalculateEntropyProfile_Tests.cs | ≥ 80% |
| 233 | ☑ | SEQ-GC-ANALYSIS-001 | Composition | GcSkewCalculator.cs | GcSkewCalculator_AnalyzeGcContent_Tests.cs | ≥ 80% |
| 234 | ☐ | SEQ-GC-PROFILE-001 | Statistics | SequenceStatistics.cs | SequenceStatistics_CalculateGcContentProfile_Tests.cs | ≥ 80% |
| 235 | ☐ | ONCO-ASCAT-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_AscatDerivation_Tests.cs | ≥ 80% |
| 236 | ☑ | RNA-PKPREDICT-001 | Analysis | RnaSecondaryStructure.cs | RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs | ≥ 80% |
| 237 | ☑ | RNA-PKRECURSIVE-001 | Analysis | RnaSecondaryStructure.cs | RnaSecondaryStructure_PredictStructurePseudoknotRecursive_Tests.cs | ≥ 80% |
| 238 | ☑ | RNA-ACCESS-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructure_UnpairedProbabilities_Tests.cs | ≥ 80% |
| 239 | ☐ | PROTMOTIF-HMM-001 | ProteinMotif | Plan7ProfileHmm.cs / ProteinMotifFinder.cs | ProteinMotifFinder_FindDomainsByHmm_Tests.cs | ≥ 80% |
| 240 | ☑ | PRIMER-NNTM-001 | MolTools | PrimerDesigner.cs | PrimerDesigner_NearestNeighborTm_Tests.cs | ≥ 80% |
| 241 | ☑ | PRIMER-HAIRPIN-001 | MolTools | PrimerDesigner.cs | PrimerDesigner_HairpinTm_Tests.cs / PrimerDesigner_HairpinSpecialLoop_Tests.cs | ≥ 80% |
| 242 | ☑ | PRIMER-DIMER-001 | MolTools | PrimerDesigner.cs / NtthalDimer.cs | PrimerDesigner_DimerTm_Tests.cs | ≥ 80% |
| 243 | ☑ | PROBE-LNATM-001 | MolTools | ProbeDesigner.cs | ProbeDesigner_LnaTm_Tests.cs | ≥ 80% |
| 244 | ☑ | PROBE-EVALUE-001 | MolTools | ProbeDesigner.cs | ProbeDesigner_ProbeValidation_Tests.cs | ≥ 80% |
| 245 | ☐ | MHC-NN-001 | Oncology | MhcflurryAffinityPredictor.cs | MhcflurryAffinityPredictor_PredictIc50_Tests.cs | ≥ 80% |
| 246 | ☐ | MHC-MATRIX-001 | Oncology | OncologyAnalyzer.cs | OncologyAnalyzer_ClassifyMhcBinding_Tests.cs | ≥ 80% |
| 247 | ☐ | IMMUNE-NUSVR-001 | Oncology | ImmuneAnalyzer.cs | ImmuneAnalyzer_ImmuneInfiltration_Tests.cs | ≥ 80% |
| 248 | ☑ | META-CHECKM-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_MarkerGeneQuality_Tests.cs | ≥ 80% |
| 249 | ☑ | META-TETRA-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_TetranucleotideZScore_Tests.cs | ≥ 80% |
| 250 | ☑ | SPLICE-MAXENT3-001 | Splicing | SpliceSitePredictor.cs | SpliceSitePredictor_AcceptorSite_Tests.cs | ≥ 80% |
| 251 | ☑ | SPLICE-MAXENT5-001 | Splicing | SpliceSitePredictor.cs | SpliceSitePredictor_DonorSite_Tests.cs | ≥ 80% |
| 252 | ☐ | MIRNA-CONTEXT-001 | MiRNA | MiRnaAnalyzer.cs | MiRnaAnalyzer_TargetPrediction_Tests.cs | ≥ 80% |
| 253 | ☐ | MIRNA-PCT-001 | MiRNA | MiRnaAnalyzer.cs | MiRnaAnalyzer_TargetPrediction_Tests.cs | ≥ 80% |
| 254 | ☐ | MIRNA-CLASSIFY-001 | MiRNA | MiRnaAnalyzer.cs | MiRnaAnalyzer_PreMiRna_Tests.cs | ≥ 80% |
| 255 | ☐ | MIRNA-CLEAVAGE-001 | MiRNA | MiRnaAnalyzer.cs | MiRnaAnalyzer_PreMiRna_Tests.cs | ≥ 80% |
| 256 | ☑ | REP-APPROX-001 | Repeats | RepeatFinder.cs | RepeatFinder_ApproximateTandemRepeats_Tests.cs | ≥ 80% |
| 257 | ☑ | CHROM-ALPHASAT-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_AlphaSatellite_Tests.cs | ≥ 80% |
| 258 | ☑ | CHROM-HOR-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_HigherOrderRepeat_Tests.cs | ≥ 80% |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 255 |
| ☑ Complete (run + killers written) | 230 |
| ☐ Not started | 28 |
| Unique source files to mutate | ~25 |
| Target mutation score per file | ≥ 80% |
