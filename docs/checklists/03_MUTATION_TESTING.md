# Checklist 03: Mutation Testing (Stryker)

**Priority:** P1  
**Framework:** Stryker.NET  
**Date:** 2026-03-19  
**Total algorithms:** 86

---

## Description

Mutation testing измеряет эффективность тестового набора: Stryker вносит мутации в исходный код (замена операторов, удаление условий и т.д.) и проверяет, что тесты их обнаруживают («убивают»). Мутанты, пережившие все тесты (survived), указывают на пробелы в покрытии. Не требует написания нового кода для запуска — только конфигурацию и анализ.

**Текущее покрытие:** `MutationKillerTests.cs` таргетирует MotifFinder.cs и RepeatFinder.cs. `stryker-config.json` существует. Результаты в `StrykerOutput/`.

**Процесс для каждого модуля:**
1. Запустить `dotnet stryker` для исходного файла
2. Проанализировать survived мутантов
3. Написать прицельные killer-тесты для survivors с score < 80%

**Target:** mutation score ≥ 80% для каждого модуля.

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | Source File | Test File(s) | Target Score |
|---|--------|-----------|------|------------|-------------|:---:|
| 1 | ☐ | SEQ-GC-001 | Composition | SequenceExtensions.cs | SequenceExtensions_CalculateGcContent_Tests.cs | ≥ 80% |
| 2 | ☐ | SEQ-COMP-001 | Composition | SequenceExtensions.cs | SequenceExtensions_Complement_Tests.cs | ≥ 80% |
| 3 | ☐ | SEQ-REVCOMP-001 | Composition | SequenceExtensions.cs | SequenceExtensions_ReverseComplement_Tests.cs | ≥ 80% |
| 4 | ☐ | SEQ-VALID-001 | Composition | SequenceExtensions.cs | SequenceExtensions_SequenceValidation_Tests.cs | ≥ 80% |
| 5 | ☐ | SEQ-COMPLEX-001 | Composition | SequenceComplexity.cs | SequenceComplexityTests.cs | ≥ 80% |
| 6 | ☐ | SEQ-ENTROPY-001 | Composition | SequenceComplexity.cs | SequenceComplexityTests.cs | ≥ 80% |
| 7 | ☐ | SEQ-GCSKEW-001 | Composition | GcSkewCalculator.cs | GcSkewCalculatorTests.cs | ≥ 80% |
| 8 | ☑ | PAT-EXACT-001 | Matching | MotifFinder.cs | FindAllOccurrencesTests.cs, ContainsTests.cs, CountOccurrencesTests.cs | ≥ 80% |
| 9 | ☐ | PAT-APPROX-001 | Matching | ApproximateMatcher.cs | ApproximateMatcher_HammingDistance_Tests.cs | ≥ 80% |
| 10 | ☐ | PAT-APPROX-002 | Matching | ApproximateMatcher.cs | ApproximateMatcher_EditDistance_Tests.cs | ≥ 80% |
| 11 | ☑ | PAT-IUPAC-001 | Matching | MotifFinder.cs, IupacHelper.cs | IupacMotifMatchingTests.cs, MutationKillerTests.cs | ≥ 80% |
| 12 | ☑ | PAT-PWM-001 | Matching | MotifFinder.cs | MotifFinder_PWM_Tests.cs, MutationKillerTests.cs | ≥ 80% |
| 13 | ☑ | REP-STR-001 | Repeats | RepeatFinder.cs | RepeatFinder_Microsatellite_Tests.cs | ≥ 80% |
| 14 | ☑ | REP-TANDEM-001 | Repeats | RepeatFinder.cs, GenomicAnalyzer.cs | GenomicAnalyzer_TandemRepeat_Tests.cs, RepeatFinderTests.cs | ≥ 80% |
| 15 | ☑ | REP-INV-001 | Repeats | RepeatFinder.cs | RepeatFinder_InvertedRepeat_Tests.cs | ≥ 80% |
| 16 | ☑ | REP-DIRECT-001 | Repeats | RepeatFinder.cs | RepeatFinder_DirectRepeat_Tests.cs | ≥ 80% |
| 17 | ☑ | REP-PALIN-001 | Repeats | RepeatFinder.cs | RepeatFinder_Palindrome_Tests.cs | ≥ 80% |
| 18 | ☐ | CRISPR-PAM-001 | MolTools | CrisprDesigner.cs | CrisprDesigner_PAM_Tests.cs | ≥ 80% |
| 19 | ☐ | CRISPR-GUIDE-001 | MolTools | CrisprDesigner.cs | CrisprDesigner_GuideRNA_Tests.cs | ≥ 80% |
| 20 | ☐ | CRISPR-OFF-001 | MolTools | CrisprDesigner.cs | CrisprDesigner_OffTarget_Tests.cs | ≥ 80% |
| 21 | ☐ | PRIMER-TM-001 | MolTools | PrimerDesigner.cs | PrimerDesigner_MeltingTemperature_Tests.cs | ≥ 80% |
| 22 | ☐ | PRIMER-DESIGN-001 | MolTools | PrimerDesigner.cs | PrimerDesigner_PrimerDesign_Tests.cs | ≥ 80% |
| 23 | ☐ | PRIMER-STRUCT-001 | MolTools | PrimerDesigner.cs | PrimerDesigner_PrimerStructure_Tests.cs | ≥ 80% |
| 24 | ☐ | PROBE-DESIGN-001 | MolTools | ProbeDesigner.cs | ProbeDesigner_ProbeDesign_Tests.cs | ≥ 80% |
| 25 | ☐ | PROBE-VALID-001 | MolTools | ProbeDesigner.cs | ProbeDesigner_ProbeValidation_Tests.cs | ≥ 80% |
| 26 | ☐ | RESTR-FIND-001 | MolTools | RestrictionAnalyzer.cs | RestrictionAnalyzer_FindSites_Tests.cs | ≥ 80% |
| 27 | ☐ | RESTR-DIGEST-001 | MolTools | RestrictionAnalyzer.cs | RestrictionAnalyzer_Digest_Tests.cs | ≥ 80% |
| 28 | ☐ | ANNOT-ORF-001 | Annotation | GenomeAnnotator.cs | GenomeAnnotator_ORF_Tests.cs | ≥ 80% |
| 29 | ☐ | ANNOT-GENE-001 | Annotation | GenomeAnnotator.cs | GenomeAnnotator_Gene_Tests.cs | ≥ 80% |
| 30 | ☐ | ANNOT-PROM-001 | Annotation | GenomeAnnotator.cs | GenomeAnnotator_PromoterMotif_Tests.cs | ≥ 80% |
| 31 | ☐ | ANNOT-GFF-001 | Annotation | GenomeAnnotator.cs | GenomeAnnotator_GFF3_Tests.cs | ≥ 80% |
| 32 | ☐ | KMER-COUNT-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_CountKmers_Tests.cs | ≥ 80% |
| 33 | ☐ | KMER-FREQ-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_Frequency_Tests.cs | ≥ 80% |
| 34 | ☐ | KMER-FIND-001 | K-mer | KmerAnalyzer.cs | KmerAnalyzer_Find_Tests.cs | ≥ 80% |
| 35 | ☐ | ALIGN-GLOBAL-001 | Alignment | SequenceAligner.cs | SequenceAligner_GlobalAlign_Tests.cs | ≥ 80% |
| 36 | ☐ | ALIGN-LOCAL-001 | Alignment | SequenceAligner.cs | SequenceAligner_LocalAlign_Tests.cs | ≥ 80% |
| 37 | ☐ | ALIGN-SEMI-001 | Alignment | SequenceAligner.cs | SequenceAligner_SemiGlobalAlign_Tests.cs | ≥ 80% |
| 38 | ☐ | ALIGN-MULTI-001 | Alignment | SequenceAligner.cs | SequenceAligner_MultipleAlign_Tests.cs | ≥ 80% |
| 39 | ☐ | PHYLO-DIST-001 | Phylogenetic | PhylogeneticAnalyzer.cs | PhylogeneticAnalyzer_DistanceMatrix_Tests.cs | ≥ 80% |
| 40 | ☐ | PHYLO-TREE-001 | Phylogenetic | PhylogeneticAnalyzer.cs | PhylogeneticAnalyzer_TreeConstruction_Tests.cs | ≥ 80% |
| 41 | ☐ | PHYLO-NEWICK-001 | Phylogenetic | PhylogeneticAnalyzer.cs | PhylogeneticAnalyzer_NewickIO_Tests.cs | ≥ 80% |
| 42 | ☐ | PHYLO-COMP-001 | Phylogenetic | PhylogeneticAnalyzer.cs | PhylogeneticAnalyzer_TreeComparison_Tests.cs | ≥ 80% |
| 43 | ☐ | POP-FREQ-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs | ≥ 80% |
| 44 | ☐ | POP-DIV-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_Diversity_Tests.cs | ≥ 80% |
| 45 | ☐ | POP-HW-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs | ≥ 80% |
| 46 | ☐ | POP-FST-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_FStatistics_Tests.cs | ≥ 80% |
| 47 | ☐ | POP-LD-001 | PopGen | PopulationGeneticsAnalyzer.cs | PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs | ≥ 80% |
| 48 | ☐ | CHROM-TELO-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_Telomere_Tests.cs | ≥ 80% |
| 49 | ☐ | CHROM-CENT-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_Centromere_Tests.cs | ≥ 80% |
| 50 | ☐ | CHROM-KARYO-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_Karyotype_Tests.cs | ≥ 80% |
| 51 | ☐ | CHROM-ANEU-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_Aneuploidy_Tests.cs | ≥ 80% |
| 52 | ☐ | CHROM-SYNT-001 | Chromosome | ChromosomeAnalyzer.cs | ChromosomeAnalyzer_Synteny_Tests.cs | ≥ 80% |
| 53 | ☐ | META-CLASS-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs | ≥ 80% |
| 54 | ☐ | META-PROF-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_TaxonomicProfile_Tests.cs | ≥ 80% |
| 55 | ☐ | META-ALPHA-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_AlphaDiversity_Tests.cs | ≥ 80% |
| 56 | ☐ | META-BETA-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_BetaDiversity_Tests.cs | ≥ 80% |
| 57 | ☐ | META-BIN-001 | Metagenomics | MetagenomicsAnalyzer.cs | MetagenomicsAnalyzer_GenomeBinning_Tests.cs | ≥ 80% |
| 58 | ☐ | CODON-OPT-001 | Codon | CodonOptimizer.cs | CodonOptimizer_OptimizeSequence_Tests.cs | ≥ 80% |
| 59 | ☐ | CODON-CAI-001 | Codon | CodonOptimizer.cs | CodonOptimizer_CAI_Tests.cs | ≥ 80% |
| 60 | ☐ | CODON-RARE-001 | Codon | CodonOptimizer.cs | CodonOptimizer_FindRareCodons_Tests.cs | ≥ 80% |
| 61 | ☐ | CODON-USAGE-001 | Codon | CodonOptimizer.cs | CodonOptimizer_CodonUsage_Tests.cs | ≥ 80% |
| 62 | ☐ | TRANS-CODON-001 | Translation | GeneticCode.cs | GeneticCodeTests.cs | ≥ 80% |
| 63 | ☐ | TRANS-PROT-001 | Translation | Translator.cs | TranslatorTests.cs | ≥ 80% |
| 64 | ☐ | PARSE-FASTA-001 | FileIO | FastaParser.cs | FastaParserTests.cs | ≥ 80% |
| 65 | ☐ | PARSE-FASTQ-001 | FileIO | FastqParser.cs | FastqParserTests.cs | ≥ 80% |
| 66 | ☐ | PARSE-BED-001 | FileIO | BedParser.cs | BedParserTests.cs | ≥ 80% |
| 67 | ☐ | PARSE-VCF-001 | FileIO | VcfParser.cs | VcfParserTests.cs | ≥ 80% |
| 68 | ☐ | PARSE-GFF-001 | FileIO | GffParser.cs | GffParserTests.cs | ≥ 80% |
| 69 | ☐ | PARSE-GENBANK-001 | FileIO | GenBankParser.cs | GenBankParserTests.cs | ≥ 80% |
| 70 | ☐ | PARSE-EMBL-001 | FileIO | EmblParser.cs | EmblParserTests.cs | ≥ 80% |
| 71 | ☐ | RNA-STRUCT-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructureTests.cs | ≥ 80% |
| 72 | ☐ | RNA-STEMLOOP-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructureTests.cs | ≥ 80% |
| 73 | ☐ | RNA-ENERGY-001 | RnaStructure | RnaSecondaryStructure.cs | RnaSecondaryStructureTests.cs | ≥ 80% |
| 74 | ☐ | MIRNA-SEED-001 | MiRNA | MiRnaAnalyzer.cs | MiRnaAnalyzer_SeedAnalysis_Tests.cs | ≥ 80% |
| 75 | ☐ | MIRNA-TARGET-001 | MiRNA | MiRnaAnalyzer.cs | MiRnaAnalyzer_TargetPrediction_Tests.cs | ≥ 80% |
| 76 | ☐ | MIRNA-PRECURSOR-001 | MiRNA | MiRnaAnalyzer.cs | MiRnaAnalyzer_PreMiRna_Tests.cs | ≥ 80% |
| 77 | ☐ | SPLICE-DONOR-001 | Splicing | SpliceSitePredictor.cs | SpliceSitePredictor_DonorSite_Tests.cs | ≥ 80% |
| 78 | ☐ | SPLICE-ACCEPTOR-001 | Splicing | SpliceSitePredictor.cs | SpliceSitePredictor_AcceptorSite_Tests.cs | ≥ 80% |
| 79 | ☐ | SPLICE-PREDICT-001 | Splicing | SpliceSitePredictor.cs | SpliceSitePredictor_GeneStructure_Tests.cs | ≥ 80% |
| 80 | ☐ | DISORDER-PRED-001 | ProteinPred | DisorderPredictor.cs | DisorderPredictor_DisorderPrediction_Tests.cs | ≥ 80% |
| 81 | ☐ | DISORDER-REGION-001 | ProteinPred | DisorderPredictor.cs | DisorderPredictor_DisorderedRegion_Tests.cs | ≥ 80% |
| 82 | ☐ | PROTMOTIF-FIND-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_MotifSearch_Tests.cs | ≥ 80% |
| 83 | ☐ | PROTMOTIF-PROSITE-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_PrositePattern_Tests.cs | ≥ 80% |
| 84 | ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | ProteinMotifFinder.cs | ProteinMotifFinder_DomainPrediction_Tests.cs | ≥ 80% |
| 85 | ☐ | EPIGEN-CPG-001 | Epigenetics | EpigeneticsAnalyzer.cs | EpigeneticsAnalyzer_CpGDetection_Tests.cs | ≥ 80% |
| 86 | ☐ | ONCO-IMMUNE-001 | Oncology | ImmuneAnalyzer.cs | ImmuneAnalyzer_ImmuneInfiltration_Tests.cs | ≥ 80% |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 86 |
| ☑ Complete (run + killers written) | 8 |
| ☐ Not started | 78 |
| Unique source files to mutate | ~25 |
| Target mutation score per file | ≥ 80% |
