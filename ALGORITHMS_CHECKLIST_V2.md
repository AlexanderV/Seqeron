# Algorithms Checklist v2.0

**Date:** 2026-02-12
**Version:** 2.5
**Library:** Seqeron.Genomics

---

## Quick Reference

| Metric | Value |
|--------|-------|
| **Total Test Units** | 255 |
| **Completed** | 254 |
| **In Progress** | 0 |
| **Blocked** | 0 |
| **Not Started** | 1 |

---

## Processing Registry

| Status | Test Unit ID | Area | Methods | Evidence | TestSpec | Test File(s) |
|--------|--------------|------|---------|----------|----------|--------------|
| ☑ | SEQ-GC-001 | Composition | 5 | Wikipedia, Biopython | [SEQ-GC-001.md](TestSpecs/SEQ-GC-001.md) | SequenceExtensions_CalculateGcContent_Tests.cs |
| ☑ | SEQ-COMP-001 | Composition | 3 | Wikipedia, Biopython | [SEQ-COMP-001.md](TestSpecs/SEQ-COMP-001.md) | SequenceExtensions_Complement_Tests.cs |
| ☑ | SEQ-REVCOMP-001 | Composition | 4 | Wikipedia, Biopython | [SEQ-REVCOMP-001.md](TestSpecs/SEQ-REVCOMP-001.md) | SequenceExtensions_ReverseComplement_Tests.cs |
| ☑ | SEQ-VALID-001 | Composition | 4 | Wikipedia, IUPAC 1970, Bioinformatics.org | [SEQ-VALID-001.md](TestSpecs/SEQ-VALID-001.md) | SequenceExtensions_SequenceValidation_Tests.cs |
| ☑ | SEQ-COMPLEX-001 | Composition | 4 | Wikipedia, Troyanskaya (2002), Orlov (2004) | [SEQ-COMPLEX-001.md](TestSpecs/SEQ-COMPLEX-001.md) | SequenceComplexityTests.cs |
| ☑ | SEQ-ENTROPY-001 | Composition | 2 | Wikipedia (Entropy, Sequence logo, K-mer), Shannon (1948) | [SEQ-ENTROPY-001.md](TestSpecs/SEQ-ENTROPY-001.md) | SequenceComplexityTests.cs |
| ☑ | SEQ-GCSKEW-001 | Composition | 4 | Wikipedia, Lobry (1996), Grigoriev (1998) | [SEQ-GCSKEW-001.md](TestSpecs/SEQ-GCSKEW-001.md) | GcSkewCalculatorTests.cs |
| ☑ | PAT-EXACT-001 | Matching | 4 | Wikipedia, Gusfield (1997), Rosalind | [PAT-EXACT-001.md](TestSpecs/PAT-EXACT-001.md) | FindAllOccurrencesTests.cs, ContainsTests.cs, CountOccurrencesTests.cs |
| ☑ | PAT-APPROX-001 | Matching | 2 | Wikipedia (Hamming), Rosalind (HAMM), Gusfield (1997), Navarro (2001) | [PAT-APPROX-001.md](TestSpecs/PAT-APPROX-001.md) | ApproximateMatcher_HammingDistance_Tests.cs |
| ☑ | PAT-APPROX-002 | Matching | 2 | Wikipedia (Levenshtein, Edit Distance), Rosetta Code, Navarro (2001) | [PAT-APPROX-002.md](TestSpecs/PAT-APPROX-002.md) | ApproximateMatcher_EditDistance_Tests.cs |
| ☑ | PAT-IUPAC-001 | Matching | 2 | Wikipedia (Nucleic acid notation), IUPAC-IUB 1970, Bioinformatics.org | [PAT-IUPAC-001.md](TestSpecs/PAT-IUPAC-001.md) | IupacMotifMatchingTests.cs |
| ☑ | PAT-PWM-001 | Matching | 2 | Wikipedia (PWM), Kel et al. (2003), Rosalind (CONS), Nishida (2008) | [PAT-PWM-001.md](TestSpecs/PAT-PWM-001.md) | MotifFinder_PWM_Tests.cs |
| ☑ | REP-STR-001 | Repeats | 4 | Wikipedia (Microsatellite, Trinucleotide repeat disorder), Richard et al. (2008), Benson (1999) TRF | [REP-STR-001.md](TestSpecs/REP-STR-001.md) | RepeatFinder_Microsatellite_Tests.cs |
| ☑ | REP-TANDEM-001 | Repeats | 2 | Wikipedia (Tandem repeat, Microsatellite), Richard et al. (2008) | [REP-TANDEM-001.md](TestSpecs/REP-TANDEM-001.md) | GenomicAnalyzer_TandemRepeat_Tests.cs, RepeatFinderTests.cs |
| ☑ | REP-INV-001 | Repeats | 1 | Wikipedia (Inverted repeat, Stem-loop, Palindromic sequence), EMBOSS einverted, Pearson (1996), Bissler (1998) | [REP-INV-001.md](TestSpecs/REP-INV-001.md) | RepeatFinder_InvertedRepeat_Tests.cs |
| ☑ | REP-DIRECT-001 | Repeats | 1 | Wikipedia (Direct repeat, Repeated sequence), Ussery (2009), Richard (2021) | [REP-DIRECT-001.md](TestSpecs/REP-DIRECT-001.md) | RepeatFinder_DirectRepeat_Tests.cs |
| ☑ | REP-PALIN-001 | Repeats | 2 | Wikipedia (Palindromic sequence, Restriction enzyme), Rosalind REVP | [REP-PALIN-001.md](TestSpecs/REP-PALIN-001.md) | RepeatFinder_Palindrome_Tests.cs |
| ☑ | CRISPR-PAM-001 | MolTools | 2 | Wikipedia (Protospacer adjacent motif, CRISPR), Jinek et al. (2012), Zetsche et al. (2015) | [CRISPR-PAM-001.md](TestSpecs/CRISPR-PAM-001.md) | CrisprDesigner_PAM_Tests.cs |
| ☑ | CRISPR-GUIDE-001 | MolTools | 2 | Wikipedia (Guide RNA, CRISPR gene editing, PAM), Addgene CRISPR Guide | [CRISPR-GUIDE-001.md](TestSpecs/CRISPR-GUIDE-001.md) | CrisprDesigner_GuideRNA_Tests.cs |
| ☑ | CRISPR-OFF-001 | MolTools | 2 | Wikipedia (Off-target genome editing), Hsu et al. (2013), Fu et al. (2013) | [CRISPR-OFF-001.md](TestSpecs/CRISPR-OFF-001.md) | CrisprDesigner_OffTarget_Tests.cs |
| ☑ | PRIMER-TM-001 | MolTools | 2 | Wikipedia (Nucleic acid thermodynamics), Marmur & Doty (1962), SantaLucia (1998), Owczarzy (2004) | [PRIMER-TM-001.md](TestSpecs/PRIMER-TM-001.md) | PrimerDesigner_MeltingTemperature_Tests.cs |
| ☑ | PRIMER-DESIGN-001 | MolTools | 3 | Wikipedia, Addgene, Primer3 | [PRIMER-DESIGN-001.md](TestSpecs/PRIMER-DESIGN-001.md) | PrimerDesigner_PrimerDesign_Tests.cs |
| ☑ | PRIMER-STRUCT-001 | MolTools | 3 | Wikipedia (Primer, Primer dimer, Stem-loop, Nucleic acid thermodynamics), Primer3 Manual, SantaLucia (1998) | [PRIMER-STRUCT-001.md](TestSpecs/PRIMER-STRUCT-001.md) | PrimerDesigner_PrimerStructure_Tests.cs |
| ☑ | PROBE-DESIGN-001 | MolTools | 3 | Wikipedia (Nucleic acid thermodynamics, Hybridization probe, FISH, DNA microarray), SantaLucia (1998), Applied Biosystems/Thermo Fisher/PREMIER Biosoft (TaqMan rules) | [PROBE-DESIGN-001.md](TestSpecs/PROBE-DESIGN-001.md) | ProbeDesigner_ProbeDesign_Tests.cs, ProbeDesigner_TaqMan_Tests.cs |
| ☑ | PROBE-VALID-001 | MolTools | 4 | Wikipedia (Hybridization probe, DNA microarray, Off-target genome editing, BLAST), Smith & Waterman (1981), Altschul et al. (1990), Kane et al. (2000), Amann & Ludwig (2000) | [PROBE-VALID-001.md](TestSpecs/PROBE-VALID-001.md) | ProbeDesigner_ProbeValidation_Tests.cs |
| ☑ | RESTR-FIND-001 | MolTools | 2 | Wikipedia (Restriction enzyme, Restriction site, EcoRI), Roberts (1976), REBASE | [RESTR-FIND-001.md](TestSpecs/RESTR-FIND-001.md) | RestrictionAnalyzer_FindSites_Tests.cs |
| ☑ | RESTR-DIGEST-001 | MolTools | 2 | Wikipedia (Restriction digest, Restriction enzyme, Restriction map), Addgene Protocol, Roberts (1976), REBASE | [RESTR-DIGEST-001.md](TestSpecs/RESTR-DIGEST-001.md) | RestrictionAnalyzer_Digest_Tests.cs |
| ☑ | ANNOT-ORF-001 | Annotation | 3 | Wikipedia (ORF), Rosalind ORF, NCBI ORF Finder, Deonier (2005), Claverie (1997) | [ANNOT-ORF-001.md](TestSpecs/ANNOT-ORF-001.md) | GenomeAnnotator_ORF_Tests.cs |
| ☑ | ANNOT-GENE-001 | Annotation | 2 | Wikipedia (Gene prediction, Shine-Dalgarno, RBS), Shine & Dalgarno (1975), Chen (1994), Laursen (2005) | [ANNOT-GENE-001.md](TestSpecs/ANNOT-GENE-001.md) | GenomeAnnotator_Gene_Tests.cs |
| ☑ | ANNOT-PROM-001 | Annotation | 1 | Wikipedia (Promoter, Pribnow box), Pribnow (1975), Harley & Reynolds (1987) | [ANNOT-PROM-001.md](TestSpecs/ANNOT-PROM-001.md) | GenomeAnnotator_PromoterMotif_Tests.cs |
| ☑ | ANNOT-GFF-001 | Annotation | 2 | Sequence Ontology GFF3 Spec v1.26, Wikipedia (GFF), RFC 3986 | [ANNOT-GFF-001.md](TestSpecs/ANNOT-GFF-001.md) | GenomeAnnotator_GFF3_Tests.cs |
| ☑ | KMER-COUNT-001 | K-mer | 3 | Wikipedia (K-mer), Rosalind (KMER, BA1E) | [KMER-COUNT-001.md](TestSpecs/KMER-COUNT-001.md) | KmerAnalyzer_CountKmers_Tests.cs |
| ☑ | KMER-FREQ-001 | K-mer | 3 | Wikipedia (K-mer, Entropy), Shannon (1948), Rosalind KMER | [KMER-FREQ-001.md](TestSpecs/KMER-FREQ-001.md) | KmerAnalyzer_Frequency_Tests.cs |
| ☑ | KMER-FIND-001 | K-mer | 3 | Wikipedia (K-mer), Rosalind BA1B (frequent words), Rosalind BA1E (clump finding) | [KMER-FIND-001.md](TestSpecs/KMER-FIND-001.md) | KmerAnalyzer_Find_Tests.cs |
| ☑ | ALIGN-GLOBAL-001 | Alignment | 1 | Wikipedia (Needleman–Wunsch, Sequence alignment) | [ALIGN-GLOBAL-001.md](TestSpecs/ALIGN-GLOBAL-001.md) | SequenceAligner_GlobalAlign_Tests.cs, PerformanceExtensionsTests.cs |
| ☑ | ALIGN-LOCAL-001 | Alignment | 1 | Wikipedia (Smith–Waterman, Sequence alignment) | [ALIGN-LOCAL-001.md](TestSpecs/ALIGN-LOCAL-001.md) | SequenceAligner_LocalAlign_Tests.cs |
| ☑ | ALIGN-SEMI-001 | Alignment | 1 | Wikipedia (Sequence alignment, Needleman–Wunsch, Smith–Waterman) | [ALIGN-SEMI-001.md](TestSpecs/ALIGN-SEMI-001.md) | SequenceAligner_SemiGlobalAlign_Tests.cs |
| ☑ | ALIGN-MULTI-001 | Alignment | 1 | Wikipedia (Multiple sequence alignment, Clustal) | [ALIGN-MULTI-001.md](TestSpecs/ALIGN-MULTI-001.md) | SequenceAligner_MultipleAlign_Tests.cs |
| ☑ | PHYLO-DIST-001 | Phylogenetic | 2 | Wikipedia (Models of DNA evolution, Distance matrices in phylogeny, Jukes-Cantor, Kimura 2-parameter) | [PHYLO-DIST-001.md](TestSpecs/PHYLO-DIST-001.md) | PhylogeneticAnalyzer_DistanceMatrix_Tests.cs |
| ☑ | PHYLO-TREE-001 | Phylogenetic | 1 | Wikipedia (UPGMA, Neighbor joining, Phylogenetic tree), Saitou & Nei (1987), Sokal & Michener (1958) | [PHYLO-TREE-001.md](TestSpecs/PHYLO-TREE-001.md) | PhylogeneticAnalyzer_TreeConstruction_Tests.cs |
| ☑ | PHYLO-NEWICK-001 | Phylogenetic | 2 | Wikipedia (Newick format), PHYLIP (Felsenstein), Olsen (1990) | [PHYLO-NEWICK-001.md](TestSpecs/PHYLO-NEWICK-001.md) | PhylogeneticAnalyzer_NewickIO_Tests.cs |
| ☑ | PHYLO-COMP-001 | Phylogenetic | 3 | Wikipedia (Robinson-Foulds, MRCA, Phylogenetic tree), Robinson & Foulds (1981) | [PHYLO-COMP-001.md](TestSpecs/PHYLO-COMP-001.md) | PhylogeneticAnalyzer_TreeComparison_Tests.cs |
| ☑ | POP-FREQ-001 | PopGen | 3 | Wikipedia (Allele frequency, Minor allele frequency, Genotype frequency), Gillespie (2004) | [POP-FREQ-001.md](TestSpecs/POP-FREQ-001.md) | PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs |
| ☑ | POP-DIV-001 | PopGen | 4 | Wikipedia (Nucleotide diversity, Watterson estimator, Tajima's D), Nei & Li (1979), Watterson (1975), Tajima (1989) | [POP-DIV-001.md](TestSpecs/POP-DIV-001.md) | PopulationGeneticsAnalyzer_Diversity_Tests.cs |
| ☑ | POP-HW-001 | PopGen | 1 | Wikipedia (Hardy-Weinberg principle, Chi-squared test), Hardy (1908), Weinberg (1908), Emigh (1980) | [POP-HW-001.md](TestSpecs/POP-HW-001.md) | PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs |
| ☑ | POP-FST-001 | PopGen | 2 | Wikipedia (Fixation index, F-statistics), Wright (1965), Weir & Cockerham (1984) | [POP-FST-001.md](TestSpecs/POP-FST-001.md) | PopulationGeneticsAnalyzer_FStatistics_Tests.cs |
| ☑ | POP-LD-001 | PopGen | 2 | Wikipedia (Linkage disequilibrium, Haplotype block), Lewontin (1964), Hill & Robertson (1968), Gabriel (2002) | [POP-LD-001.md](TestSpecs/POP-LD-001.md) | PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs |
| ☑ | CHROM-TELO-001 | Chromosome | 2 | Wikipedia (Telomere), Meyne (1989), Cawthon (2002) | [CHROM-TELO-001.md](TestSpecs/CHROM-TELO-001.md) | ChromosomeAnalyzer_Telomere_Tests.cs |
| ☐ | CHROM-CENT-001 | Chromosome | 3 | Willard (1985), Waye & Willard (1987), Masumoto et al. (1989), Levan (1964), McNulty & Sullivan (2018), Shepelev et al. (2009) | [CHROM-CENT-001.md](TestSpecs/CHROM-CENT-001.md) | ChromosomeAnalyzer_Centromere_Tests.cs, ChromosomeAnalyzer_AlphaSatellite_Tests.cs, ChromosomeAnalyzer_HigherOrderRepeat_Tests.cs, ChromosomeAnalyzer_SuprachromosomalFamily_Tests.cs |
| ☑ | CHROM-KARYO-001 | Chromosome | 2 | Wikipedia (Karyotype, Ploidy, Aneuploidy) | [CHROM-KARYO-001.md](TestSpecs/CHROM-KARYO-001.md) | ChromosomeAnalyzer_Karyotype_Tests.cs |
| ☑ | CHROM-ANEU-001 | Chromosome | 2 | Wikipedia (Aneuploidy, Copy Number Variation), Griffiths et al. (2000) | [CHROM-ANEU-001.md](TestSpecs/CHROM-ANEU-001.md) | ChromosomeAnalyzer_Aneuploidy_Tests.cs |
| ☑ | CHROM-SYNT-001 | Chromosome | 2 | Wikipedia (Synteny, Comparative genomics, Chromosomal rearrangement), Wang et al. (2012), Goel et al. (2019) | [CHROM-SYNT-001.md](TestSpecs/CHROM-SYNT-001.md) | ChromosomeAnalyzer_Synteny_Tests.cs |
| ☑ | META-CLASS-001 | Metagenomics | 2 | Wikipedia (Metagenomics), Kraken CCB JHU, Wood & Salzberg (2014) | [META-CLASS-001.md](TestSpecs/META-CLASS-001.md) | MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs |
| ☑ | META-PROF-001 | Metagenomics | 1 | Wikipedia (Metagenomics, Relative Abundance), Shannon (1948), Simpson (1949), Segata et al. (2012) | [META-PROF-001.md](TestSpecs/META-PROF-001.md) | MetagenomicsAnalyzer_TaxonomicProfile_Tests.cs |
| ☑ | META-ALPHA-001 | Metagenomics | 1 | Wikipedia (Diversity Index, Alpha Diversity), Shannon (1948), Simpson (1949), Hill (1973), Chao (1984) | [META-ALPHA-001.md](TestSpecs/META-ALPHA-001.md) | MetagenomicsAnalyzer_AlphaDiversity_Tests.cs |
| ☑ | META-BETA-001 | Metagenomics | 1 | Wikipedia (Beta diversity, Bray-Curtis dissimilarity, Jaccard index), Whittaker (1960), Bray & Curtis (1957), Jaccard (1901) | [META-BETA-001.md](TestSpecs/META-BETA-001.md) | MetagenomicsAnalyzer_BetaDiversity_Tests.cs |
| ☑ | META-BIN-001 | Metagenomics | 1 | Wikipedia (Binning metagenomics), Teeling (2004), Parks et al. (2014), Maguire et al. (2020) | [META-BIN-001.md](TestSpecs/META-BIN-001.md) | MetagenomicsAnalyzer_GenomeBinning_Tests.cs |
| ☑ | CODON-OPT-001 | Codon | 1 | Wikipedia (Codon usage bias, CAI), Sharp & Li (1987), Plotkin & Kudla (2011) | [CODON-OPT-001.md](TestSpecs/CODON-OPT-001.md) | CodonOptimizer_OptimizeSequence_Tests.cs |
| ☑ | CODON-CAI-001 | Codon | 1 | Wikipedia (Codon Adaptation Index), Sharp & Li (1987), Jansen et al. (2003) | [CODON-CAI-001.md](TestSpecs/CODON-CAI-001.md) | CodonOptimizer_CAI_Tests.cs |
| ☑ | CODON-RARE-001 | Codon | 1 | Kazusa, Shu et al. (2006), Sharp & Li (1987), Clarke & Clark (2008), Chartier et al. (2012) | [CODON-RARE-001.md](TestSpecs/CODON-RARE-001.md) | CodonOptimizer_FindRareCodons_Tests.cs, CodonOptimizer_RareCodonClusters_Tests.cs |
| ☑ | CODON-USAGE-001 | Codon | 2 | Wikipedia (Codon usage bias), Kazusa, Sharp & Li (1987) | [CODON-USAGE-001.md](TestSpecs/CODON-USAGE-001.md) | CodonOptimizer_CodonUsage_Tests.cs |
| ☑ | TRANS-CODON-001 | Translation | 3 | TRANS-CODON-001-Evidence.md | [TRANS-CODON-001.md](TestSpecs/TRANS-CODON-001.md) | GeneticCodeTests.cs |
| ☑ | TRANS-PROT-001 | Translation | 1 | Wikipedia (Translation, Reading frame, ORF), NCBI Genetic Codes | [TRANS-PROT-001.md](TestSpecs/TRANS-PROT-001.md) | TranslatorTests.cs |
| ☑ | PARSE-FASTA-001 | FileIO | 4 | Wikipedia (FASTA format), NCBI BLAST Help, Lipman & Pearson (1985) | [PARSE-FASTA-001.md](TestSpecs/PARSE-FASTA-001.md) | FastaParserTests.cs |
| ☑ | PARSE-FASTQ-001 | FileIO | 4 | Wikipedia (FASTQ format), Cock et al. (2009), NCBI SRA File Format Guide | [PARSE-FASTQ-001.md](TestSpecs/PARSE-FASTQ-001.md) | FastqParserTests.cs |
| ☑ | PARSE-BED-001 | FileIO | 6 | Wikipedia (BED format), UCSC Genome Browser FAQ, Quinlan & Hall (2010) | [PARSE-BED-001.md](TestSpecs/PARSE-BED-001.md) | BedParserTests.cs |
| ☑ | PARSE-VCF-001 | FileIO | 4 | Wikipedia (Variant Call Format), Danecek et al. (2011), SAMtools HTS-specs VCFv4.3 | [PARSE-VCF-001.md](TestSpecs/PARSE-VCF-001.md) | VcfParserTests.cs |
| ☑ | PARSE-GFF-001 | FileIO | 3 | Wikipedia (GFF), UCSC Genome Browser FAQ, Sequence Ontology GFF3 Spec v1.26 | [PARSE-GFF-001.md](TestSpecs/PARSE-GFF-001.md) | GffParserTests.cs |
| ☑ | PARSE-GENBANK-001 | FileIO | 3 | NCBI Sample Record, Wikipedia (GenBank), INSDC Feature Table Definition | [PARSE-GENBANK-001.md](TestSpecs/PARSE-GENBANK-001.md) | GenBankParserTests.cs |
| ☑ | PARSE-EMBL-001 | FileIO | 2 | EBI EMBL User Manual, INSDC Feature Table Definition v11.3 | [PARSE-EMBL-001.md](TestSpecs/PARSE-EMBL-001.md) | EmblParserTests.cs |
| ☑ | RNA-STRUCT-001 | RnaStructure | 4 | Wikipedia (Nucleic acid structure prediction, Nussinov algorithm), Nussinov (1980), Zuker (1981), Turner (2004) | [RNA-STRUCT-001.md](TestSpecs/RNA-STRUCT-001.md) | RnaSecondaryStructureTests.cs |
| ☑ | RNA-STEMLOOP-001 | RnaStructure | 3 | Wikipedia (Stem-loop, Tetraloop, Pseudoknot), Woese (1990), Heus & Pardi (1991) | [RNA-STEMLOOP-001.md](TestSpecs/RNA-STEMLOOP-001.md) | RnaSecondaryStructureTests.cs |
| ☑ | RNA-ENERGY-001 | RnaStructure | 2 | Wikipedia (RNA folding, Nearest neighbor parameters), Turner (2004), NNDB | [RNA-ENERGY-001.md](TestSpecs/RNA-ENERGY-001.md) | RnaSecondaryStructureTests.cs |
| ☑ | MIRNA-SEED-001 | MiRNA | 3 | miRBase, TargetScan (Bartel Lab), Bartel (2009), Lewis (2005) | [MIRNA-SEED-001.md](TestSpecs/MIRNA-SEED-001.md) | MiRnaAnalyzer_SeedAnalysis_Tests.cs |
| ☑ | MIRNA-TARGET-001 | MiRNA | 2 | Bartel (2009) Cell 136:215-233, Lewis et al. (2005), Grimson et al. (2007), Agarwal et al. (2015) eLife 4:e05005, TargetScan 8.0 | [MIRNA-TARGET-001.md](TestSpecs/MIRNA-TARGET-001.md) | MiRnaAnalyzer_TargetPrediction_Tests.cs |
| ☑ | MIRNA-PRECURSOR-001 | MiRNA | 2 | [Evidence](docs/Evidence/MIRNA-PRECURSOR-001-Evidence.md) | [TestSpec](tests/TestSpecs/MIRNA-PRECURSOR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs) |
| ☑ | SPLICE-DONOR-001 | Splicing | 2 | [Evidence](docs/Evidence/SPLICE-DONOR-001-Evidence.md) | [TestSpec](tests/TestSpecs/SPLICE-DONOR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs) |
| ☑ | SPLICE-ACCEPTOR-001 | Splicing | 3 | [Evidence](docs/Evidence/SPLICE-ACCEPTOR-001-Evidence.md) | [TestSpec](tests/TestSpecs/SPLICE-ACCEPTOR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_AcceptorSite_Tests.cs) |
| ☑ | SPLICE-PREDICT-001 | Splicing | 3 | [Evidence](docs/Evidence/SPLICE-PREDICT-001-Evidence.md) | [TestSpec](tests/TestSpecs/SPLICE-PREDICT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_GeneStructure_Tests.cs) |
| ☑ | DISORDER-PRED-001 | ProteinPred | 3 | [Evidence](docs/Evidence/DISORDER-PRED-001-Evidence.md) | [TestSpec](tests/TestSpecs/DISORDER-PRED-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs) |
| ☑ | DISORDER-REGION-001 | ProteinPred | 2 | [Evidence](docs/Evidence/DISORDER-REGION-001-Evidence.md) | [TestSpec](tests/TestSpecs/DISORDER-REGION-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs) |
| ☑ | PROTMOTIF-FIND-001 | ProteinMotif | 3 | [Evidence](docs/Evidence/PROTMOTIF-FIND-001-Evidence.md) | [TestSpec](tests/TestSpecs/PROTMOTIF-FIND-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_MotifSearch_Tests.cs) |
| ☑ | PROTMOTIF-PROSITE-001 | ProteinMotif | 2 | [Evidence](docs/Evidence/PROTMOTIF-PROSITE-001-Evidence.md) | [TestSpec](tests/TestSpecs/PROTMOTIF-PROSITE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_PrositePattern_Tests.cs) |
| ☑ | PROTMOTIF-DOMAIN-001 | ProteinMotif | 2 | [Evidence](docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md) | [TestSpec](tests/TestSpecs/PROTMOTIF-DOMAIN-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_DomainPrediction_Tests.cs) |
| ☑ | EPIGEN-CPG-001 | Epigenetics | 3 | [Evidence](docs/Evidence/EPIGEN-CPG-001-Evidence.md) | [TestSpec](tests/TestSpecs/EPIGEN-CPG-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_CpGDetection_Tests.cs) |
| ☑ | EPIGEN-METHYL-001 | Epigenetics | 3 | [Evidence](docs/Evidence/EPIGEN-METHYL-001-Evidence.md) | [TestSpec](tests/TestSpecs/EPIGEN-METHYL-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_Methylation_Tests.cs) |
| ☑ | EPIGEN-DMR-001 | Epigenetics | 2 | [Evidence](docs/Evidence/EPIGEN-DMR-001-Evidence.md) | [TestSpec](tests/TestSpecs/EPIGEN-DMR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_DMR_Tests.cs) |
| ☑ | VARIANT-CALL-001 | Variants | 3 | [Evidence](docs/Evidence/VARIANT-CALL-001-Evidence.md) | [TestSpec](tests/TestSpecs/VARIANT-CALL-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/VariantCaller_CallVariants_Tests.cs) |
| ☑ | VARIANT-SNP-001 | Variants | 2 | [Evidence](docs/Evidence/VARIANT-SNP-001-Evidence.md) | [TestSpec](tests/TestSpecs/VARIANT-SNP-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/VariantCaller_FindSnps_Tests.cs) |
| ☑ | VARIANT-INDEL-001 | Variants | 2 | [Evidence](docs/Evidence/VARIANT-INDEL-001-Evidence.md) | [TestSpec](tests/TestSpecs/VARIANT-INDEL-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/VariantCaller_FindIndels_Tests.cs) |
| ☑ | VARIANT-ANNOT-001 | Variants | 2 | [Evidence](docs/Evidence/VARIANT-ANNOT-001-Evidence.md) | [TestSpec](tests/TestSpecs/VARIANT-ANNOT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/VariantAnnotator_FunctionalImpact_Tests.cs) |
| ☑ | SV-DETECT-001 | StructuralVar | 3 | [Evidence](docs/Evidence/SV-DETECT-001-Evidence.md) | [TestSpec](tests/TestSpecs/SV-DETECT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/StructuralVariantAnalyzer_DetectSVs_Tests.cs) |
| ☑ | SV-BREAKPOINT-001 | StructuralVar | 2 | [Evidence](docs/Evidence/SV-BREAKPOINT-001-Evidence.md) | [TestSpec](tests/TestSpecs/SV-BREAKPOINT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/StructuralVariantAnalyzer_FindBreakpoints_Tests.cs) |
| ☑ | SV-CNV-001 | StructuralVar | 2 | [Evidence](docs/Evidence/SV-CNV-001-Evidence.md) | [TestSpec](tests/TestSpecs/SV-CNV-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/StructuralVariantAnalyzer_DetectCNV_Tests.cs) |
| ☑ | ASSEMBLY-OLC-001 | Assembly | 2 | [Evidence](docs/Evidence/ASSEMBLY-OLC-001-Evidence.md) | [TestSpec](tests/TestSpecs/ASSEMBLY-OLC-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_AssembleOLC_Tests.cs) |
| ☑ | ASSEMBLY-DBG-001 | Assembly | 2 | [Evidence](docs/Evidence/ASSEMBLY-DBG-001-Evidence.md) | [TestSpec](tests/TestSpecs/ASSEMBLY-DBG-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_AssembleDeBruijn_Tests.cs) |
| ☑ | ASSEMBLY-STATS-001 | Assembly | 4 | [Evidence](docs/Evidence/ASSEMBLY-STATS-001-Evidence.md) | [TestSpec](tests/TestSpecs/ASSEMBLY-STATS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GenomeAssemblyAnalyzer_AssemblyStatistics_Tests.cs) |
| ☑ | TRANS-EXPR-001 | Transcriptome | 3 | [Evidence](docs/Evidence/TRANS-EXPR-001-Evidence.md) | [TestSpec](tests/TestSpecs/TRANS-EXPR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/TranscriptomeAnalyzer_ExpressionQuantification_Tests.cs) |
| ☑ | TRANS-DIFF-001 | Transcriptome | 2 | [Evidence](docs/Evidence/TRANS-DIFF-001-Evidence.md) | [TestSpec](tests/TestSpecs/TRANS-DIFF-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/TranscriptomeAnalyzer_DifferentialExpression_Tests.cs) |
| ☑ | TRANS-SPLICE-001 | Transcriptome | 2 | [Evidence](docs/Evidence/TRANS-SPLICE-001-Evidence.md) | [TestSpec](tests/TestSpecs/TRANS-SPLICE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/TranscriptomeAnalyzer_AlternativeSplicing_Tests.cs) |
| ☑ | COMPGEN-SYNTENY-001 | Comparative | 2 | [Evidence](docs/Evidence/COMPGEN-SYNTENY-001-Evidence.md) | [TestSpec](tests/TestSpecs/COMPGEN-SYNTENY-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindSyntenicBlocks_Tests.cs) |
| ☑ | COMPGEN-ORTHO-001 | Comparative | 2 | [Evidence](docs/Evidence/COMPGEN-ORTHO-001-Evidence.md) | [TestSpec](tests/TestSpecs/COMPGEN-ORTHO-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindOrthologs_Tests.cs) |
| ☑ | COMPGEN-REARR-001 | Comparative | 2 | [Evidence](docs/Evidence/COMPGEN-REARR-001-Evidence.md) | [TestSpec](tests/TestSpecs/COMPGEN-REARR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_DetectRearrangements_Tests.cs) |
| ☑ | PANGEN-CORE-001 | PanGenome | 2 | [Evidence](docs/Evidence/PANGEN-CORE-001-Evidence.md) | [TestSpec](tests/TestSpecs/PANGEN-CORE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_ConstructPanGenome_Tests.cs) |
| ☑ | PANGEN-CLUSTER-001 | PanGenome | 2 | [Evidence](docs/Evidence/PANGEN-CLUSTER-001-Evidence.md) | [TestSpec](tests/TestSpecs/PANGEN-CLUSTER-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_ClusterGenes_Tests.cs) |
| ☑ | QUALITY-PHRED-001 | Quality | 3 | [Evidence](docs/Evidence/QUALITY-PHRED-001-Evidence.md) | [TestSpec](tests/TestSpecs/QUALITY-PHRED-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/QualityScoreAnalyzer_ParseQualityString_Tests.cs) |
| ☑ | QUALITY-STATS-001 | Quality | 2 | [Evidence](docs/Evidence/QUALITY-STATS-001-Evidence.md) | [TestSpec](tests/TestSpecs/QUALITY-STATS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/QualityScoreAnalyzer_CalculateStatistics_Tests.cs) |
| ☑ | SEQ-STATS-001 | Statistics | 3 | [Evidence](docs/Evidence/SEQ-STATS-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-STATS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateNucleotideComposition_Tests.cs) |
| ☑ | SEQ-MW-001 | Statistics | 2 | [Evidence](docs/Evidence/SEQ-MW-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-MW-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateMolecularWeight_Tests.cs) |
| ☑ | SEQ-PI-001 | Statistics | 1 | [Evidence](docs/Evidence/SEQ-PI-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-PI-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateIsoelectricPoint_Tests.cs) |
| ☑ | SEQ-HYDRO-001 | Statistics | 2 | [Evidence](docs/Evidence/SEQ-HYDRO-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-HYDRO-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateHydrophobicity_Tests.cs) |
| ☑ | SEQ-THERMO-001 | Statistics | 2 | [Evidence](docs/Evidence/SEQ-THERMO-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-THERMO-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateThermodynamics_Tests.cs) |
| ☑ | SEQ-DINUC-001 | Statistics | 3 | [Evidence](docs/Evidence/SEQ-DINUC-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-DINUC-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateDinucleotide_Tests.cs) |
| ☑ | SEQ-SECSTRUCT-001 | Statistics | 1 | [Evidence](docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-SECSTRUCT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_PredictSecondaryStructure_Tests.cs) |
| ☑ | CODON-RSCU-001 | Codon | 1 | [Evidence](docs/Evidence/CODON-RSCU-001-Evidence.md) | [TestSpec](tests/TestSpecs/CODON-RSCU-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/CodonUsageAnalyzer_CalculateRscu_Tests.cs) |
| ☑ | CODON-ENC-001 | Codon | 1 | [Evidence](docs/Evidence/CODON-ENC-001-Evidence.md) | [TestSpec](tests/TestSpecs/CODON-ENC-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/CodonUsageAnalyzer_CalculateEnc_Tests.cs) |
| ☑ | CODON-STATS-001 | Codon | 1 | [Evidence](docs/Evidence/CODON-STATS-001-Evidence.md) | [TestSpec](tests/TestSpecs/CODON-STATS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/CodonUsageAnalyzer_GetStatistics_Tests.cs) |
| ☑ | TRANS-SIXFRAME-001 | Translation | 1 | [Evidence](docs/Evidence/TRANS-SIXFRAME-001-Evidence.md) | [TestSpec](tests/TestSpecs/TRANS-SIXFRAME-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/Translator_SixFrames_Tests.cs) |
| ☑ | ASSEMBLY-MERGE-001 | Assembly | 1 | [Evidence](docs/Evidence/ASSEMBLY-MERGE-001-Evidence.md) | [TestSpec](tests/TestSpecs/ASSEMBLY-MERGE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_MergeContigs_Tests.cs) |
| ☑ | ASSEMBLY-SCAFFOLD-001 | Assembly | 1 | [Evidence](docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md) | [TestSpec](tests/TestSpecs/ASSEMBLY-SCAFFOLD-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_Scaffold_Tests.cs) |
| ☑ | ASSEMBLY-COVER-001 | Assembly | 1 | [Evidence](docs/Evidence/ASSEMBLY-COVER-001-Evidence.md) | [TestSpec](tests/TestSpecs/ASSEMBLY-COVER-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_CalculateCoverage_Tests.cs) |
| ☑ | ASSEMBLY-CONSENSUS-001 | Assembly | 1 | [Evidence](docs/Evidence/ASSEMBLY-CONSENSUS-001-Evidence.md) | [TestSpec](tests/TestSpecs/ASSEMBLY-CONSENSUS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_ComputeConsensus_Tests.cs) |
| ☑ | ASSEMBLY-TRIM-001 | Assembly | 1 | [Evidence](docs/Evidence/ASSEMBLY-TRIM-001-Evidence.md) | [TestSpec](tests/TestSpecs/ASSEMBLY-TRIM-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_QualityTrimReads_Tests.cs) |
| ☑ | ASSEMBLY-CORRECT-001 | Assembly | 1 | [Evidence](docs/Evidence/ASSEMBLY-CORRECT-001-Evidence.md) | [TestSpec](tests/TestSpecs/ASSEMBLY-CORRECT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_ErrorCorrectReads_Tests.cs) |
| ☑ | PAT-APPROX-003 | Matching | 3 | [Evidence](docs/Evidence/PAT-APPROX-003-Evidence.md) | [TestSpec](tests/TestSpecs/PAT-APPROX-003.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ApproximateMatcher_FindBestMatch_Tests.cs) |
| ☑ | ALIGN-STATS-001 | Alignment | 2 | [Evidence](docs/Evidence/ALIGN-STATS-001-Evidence.md) | [TestSpec](tests/TestSpecs/ALIGN-STATS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceAligner_CalculateStatistics_Tests.cs) |
| ☑ | EPIGEN-BISULF-001 | Epigenetics | 2 | [Evidence](docs/Evidence/EPIGEN-BISULF-001-Evidence.md) | [TestSpec](tests/TestSpecs/EPIGEN-BISULF-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_Bisulfite_Tests.cs) |
| ☑ | EPIGEN-CHROM-001 | Epigenetics | 3 | [Evidence](docs/Evidence/EPIGEN-CHROM-001-Evidence.md) | [TestSpec](tests/TestSpecs/EPIGEN-CHROM-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_ChromatinState_Tests.cs) |
| ☑ | EPIGEN-AGE-001 | Epigenetics | 1 | [Evidence](docs/Evidence/EPIGEN-AGE-001-Evidence.md) | [TestSpec](tests/TestSpecs/EPIGEN-AGE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests.cs) |
| ☑ | MIRNA-PAIR-001 | MiRNA | 3 | [Evidence](docs/Evidence/MIRNA-PAIR-001-Evidence.md) | [TestSpec](tests/TestSpecs/MIRNA-PAIR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_AlignMiRnaToTarget_Tests.cs) |
| ☑ | PANGEN-HEAP-001 | PanGenome | 1 | [Evidence](docs/Evidence/PANGEN-HEAP-001-Evidence.md) | [TestSpec](tests/TestSpecs/PANGEN-HEAP-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_FitHeapsLaw_Tests.cs) |
| ☑ | PANGEN-MARKER-001 | PanGenome | 2 | [Evidence](docs/Evidence/PANGEN-MARKER-001-Evidence.md) | [TestSpec](tests/TestSpecs/PANGEN-MARKER-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_SelectPhylogeneticMarkers_Tests.cs) |
| ☑ | POP-SELECT-001 | PopGen | 2 | [Evidence](docs/Evidence/POP-SELECT-001-Evidence.md) | [TestSpec](tests/TestSpecs/POP-SELECT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_SelectionSignature_Tests.cs) |
| ☑ | POP-ANCESTRY-001 | PopGen | 1 | [Evidence](docs/Evidence/POP-ANCESTRY-001-Evidence.md) | [TestSpec](tests/TestSpecs/POP-ANCESTRY-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_EstimateAncestry_Tests.cs) |
| ☑ | POP-ROH-001 | PopGen | 2 | [Evidence](docs/Evidence/POP-ROH-001-Evidence.md) | [TestSpec](tests/TestSpecs/POP-ROH-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_FindROH_Tests.cs) |
| ☑ | META-FUNC-001 | Metagenomics | 2 | [Evidence](docs/Evidence/META-FUNC-001-Evidence.md) | [TestSpec](tests/TestSpecs/META-FUNC-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_PredictFunctions_Tests.cs) |
| ☑ | META-RESIST-001 | Metagenomics | 1 | [Evidence](docs/Evidence/META-RESIST-001-Evidence.md) | [TestSpec](tests/TestSpecs/META-RESIST-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_FindAntibioticResistanceGenes_Tests.cs) |
| ☑ | META-PATHWAY-001 | Metagenomics | 1 | [Evidence](docs/Evidence/META-PATHWAY-001-Evidence.md) | [TestSpec](tests/TestSpecs/META-PATHWAY-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_FindPathwayEnrichment_Tests.cs) |
| ☑ | META-TAXA-001 | Metagenomics | 2 | [Evidence](docs/Evidence/META-TAXA-001-Evidence.md) | [TestSpec](tests/TestSpecs/META-TAXA-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_FindSignificantTaxa_Tests.cs) |
| ☑ | PHYLO-BOOT-001 | Phylogenetic | 1 | [Evidence](docs/Evidence/PHYLO-BOOT-001-Evidence.md) | [TestSpec](tests/TestSpecs/PHYLO-BOOT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_Bootstrap_Tests.cs) |
| ☑ | PHYLO-STATS-001 | Phylogenetic | 3 | [Evidence](docs/Evidence/PHYLO-STATS-001-Evidence.md) | [TestSpec](tests/TestSpecs/PHYLO-STATS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_TreeStatistics_Tests.cs) |
| ☑ | ANNOT-CODING-001 | Annotation | 1 | [Evidence](docs/Evidence/ANNOT-CODING-001-Evidence.md) | [TestSpec](tests/TestSpecs/ANNOT-CODING-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotator_CalculateCodingPotential_Tests.cs) |
| ☑ | ANNOT-REPEAT-001 | Annotation | 2 | [Evidence](docs/Evidence/ANNOT-REPEAT-001-Evidence.md) | [TestSpec](tests/TestSpecs/ANNOT-REPEAT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotator_FindRepetitiveElements_Tests.cs) |
| ☑ | ANNOT-CODONUSAGE-001 | Annotation | 1 | [Evidence](docs/Evidence/ANNOT-CODONUSAGE-001-Evidence.md) | [TestSpec](tests/TestSpecs/ANNOT-CODONUSAGE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotator_GetCodonUsage_Tests.cs) |
| ☑ | RESTR-FILTER-001 | MolTools | 3 | [Evidence](docs/Evidence/RESTR-FILTER-001-Evidence.md) | [TestSpec](tests/TestSpecs/RESTR-FILTER-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/RestrictionAnalyzer_Filter_Tests.cs) |
| ☑ | KMER-DIST-001 | K-mer | 1 | [Evidence](docs/Evidence/KMER-DIST-001-Evidence.md) | [TestSpec](tests/TestSpecs/KMER-DIST-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_KmerDistance_Tests.cs) |
| ☑ | MOTIF-CONS-001 | Matching | 1 | [Evidence](docs/Evidence/MOTIF-CONS-001-Evidence.md) | [TestSpec](tests/TestSpecs/MOTIF-CONS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_CreateConsensusFromAlignment_Tests.cs) |
| ☑ | GENOMIC-REPEAT-001 | Analysis | 2 | [Evidence](docs/Evidence/GENOMIC-REPEAT-001-Evidence.md) | [TestSpec](tests/TestSpecs/GENOMIC-REPEAT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindRepeats_Tests.cs) |
| ☑ | GENOMIC-COMMON-001 | Analysis | 2 | [Evidence](docs/Evidence/GENOMIC-COMMON-001-Evidence.md) | [TestSpec](tests/TestSpecs/GENOMIC-COMMON-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindCommonRegion_Tests.cs) |
| ☑ | GENOMIC-MOTIFS-001 | Analysis | 1 | [Evidence](docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md) | [TestSpec](tests/TestSpecs/GENOMIC-MOTIFS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindKnownMotifs_Tests.cs) |
| ☑ | SEQ-RNACOMP-001 | Composition | 1 | [Evidence](docs/Evidence/SEQ-RNACOMP-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-RNACOMP-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceExtensions_GetRnaComplementBase_Tests.cs) |
| ☑ | PROTMOTIF-PATTERN-001 | ProteinMotif | 4 | [Evidence](docs/Evidence/PROTMOTIF-PATTERN-001-Evidence.md) | [TestSpec](tests/TestSpecs/PROTMOTIF-PATTERN-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindMotifByPattern_Tests.cs) |
| ☑ | PROTMOTIF-SP-001 | ProteinMotif | 1 | [Evidence](docs/Evidence/PROTMOTIF-SP-001-Evidence.md) | [TestSpec](tests/TestSpecs/PROTMOTIF-SP-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_PredictSignalPeptide_Tests.cs) |
| ☑ | PROTMOTIF-TM-001 | ProteinMotif | 1 | [Evidence](docs/Evidence/PROTMOTIF-TM-001-Evidence.md) | [TestSpec](tests/TestSpecs/PROTMOTIF-TM-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_PredictTransmembraneHelices_Tests.cs) |
| ☑ | PROTMOTIF-CC-001 | ProteinMotif | 1 | [Evidence](docs/Evidence/PROTMOTIF-CC-001-Evidence.md) | [TestSpec](tests/TestSpecs/PROTMOTIF-CC-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_PredictCoiledCoils_Tests.cs) |
| ☑ | PROTMOTIF-LC-001 | ProteinMotif | 1 | [Evidence](docs/Evidence/PROTMOTIF-LC-001-Evidence.md) | [TestSpec](tests/TestSpecs/PROTMOTIF-LC-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindLowComplexityRegions_Tests.cs) |
| ☑ | PROTMOTIF-COMMON-001 | ProteinMotif | 2 | [Evidence](docs/Evidence/PROTMOTIF-COMMON-001-Evidence.md) | [TestSpec](tests/TestSpecs/PROTMOTIF-COMMON-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindCommonMotifs_Tests.cs) |
| ☑ | RNA-PAIR-001 | RnaStructure | 3 | [Evidence](docs/Evidence/RNA-PAIR-001-Evidence.md) | [TestSpec](tests/TestSpecs/RNA-PAIR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_CanPair_Tests.cs) |
| ☑ | RNA-HAIRPIN-001 | RnaStructure | 2 | [Evidence](docs/Evidence/RNA-HAIRPIN-001-Evidence.md) | [TestSpec](tests/TestSpecs/RNA-HAIRPIN-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_HairpinEnergy_Tests.cs) |
| ☑ | RNA-MFE-001 | RnaStructure | 2 | [Evidence](docs/Evidence/RNA-MFE-001-Evidence.md) | [TestSpec](tests/TestSpecs/RNA-MFE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_MinimumFreeEnergy_Tests.cs) |
| ☑ | RNA-PSEUDOKNOT-001 | RnaStructure | 1 | [Evidence](docs/Evidence/RNA-PSEUDOKNOT-001-Evidence.md) | [TestSpec](tests/TestSpecs/RNA-PSEUDOKNOT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_DetectPseudoknots_Tests.cs) |
| ☑ | RNA-DOTBRACKET-001 | RnaStructure | 2 | [Evidence](docs/Evidence/RNA-DOTBRACKET-001-Evidence.md) | [TestSpec](tests/TestSpecs/RNA-DOTBRACKET-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_ParseDotBracket_Tests.cs) |
| ☑ | RNA-INVERT-001 | RnaStructure | 1 | [Evidence](docs/Evidence/RNA-INVERT-001-Evidence.md) | [TestSpec](tests/TestSpecs/RNA-INVERT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_FindInvertedRepeats_Tests.cs) |
| ☑ | RNA-PARTITION-001 | RnaStructure | 4 | [Evidence](docs/Evidence/RNA-PARTITION-001-Evidence.md) | [TestSpec](tests/TestSpecs/RNA-PARTITION-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_PartitionFunction_Tests.cs) |
| ☑ | SEQ-COMPLEX-KMER-001 | Complexity | 1 | [Evidence](docs/Evidence/SEQ-COMPLEX-KMER-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-COMPLEX-KMER-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexity_CalculateKmerEntropy_Tests.cs) |
| ☑ | SEQ-COMPLEX-WINDOW-001 | Complexity | 1 | [Evidence](docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-COMPLEX-WINDOW-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexity_CalculateWindowedComplexity_Tests.cs) |
| ☑ | SEQ-COMPLEX-DUST-001 | Complexity | 2 | [Evidence](docs/Evidence/SEQ-COMPLEX-DUST-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-COMPLEX-DUST-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexity_CalculateDustScore_Tests.cs) |
| ☑ | SEQ-COMPLEX-COMPRESS-001 | Complexity | 1 | [Evidence](docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-COMPLEX-COMPRESS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexity_EstimateCompressionRatio_Tests.cs) |
| ☑ | COMPGEN-RBH-001 | Comparative | 1 | [Evidence](docs/Evidence/COMPGEN-RBH-001-Evidence.md) | [TestSpec](tests/TestSpecs/COMPGEN-RBH-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindReciprocalBestHits_Tests.cs) |
| ☑ | COMPGEN-COMPARE-001 | Comparative | 1 | [Evidence](docs/Evidence/COMPGEN-COMPARE-001-Evidence.md) | [TestSpec](tests/TestSpecs/COMPGEN-COMPARE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CompareGenomes_Tests.cs) |
| ☑ | COMPGEN-REVERSAL-001 | Comparative | 1 | [Evidence](docs/Evidence/COMPGEN-REVERSAL-001-Evidence.md) | [TestSpec](tests/TestSpecs/COMPGEN-REVERSAL-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CalculateReversalDistance_Tests.cs) |
| ☑ | COMPGEN-CLUSTER-001 | Comparative | 1 | [Evidence](docs/Evidence/COMPGEN-CLUSTER-001-Evidence.md) | [TestSpec](tests/TestSpecs/COMPGEN-CLUSTER-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindConservedClusters_Tests.cs) |
| ☑ | COMPGEN-ANI-001 | Comparative | 1 | [Evidence](docs/Evidence/COMPGEN-ANI-001-Evidence.md) | [TestSpec](tests/TestSpecs/COMPGEN-ANI-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CalculateANI_Tests.cs) |
| ☑ | COMPGEN-DOTPLOT-001 | Comparative | 1 | [Evidence](docs/Evidence/COMPGEN-DOTPLOT-001-Evidence.md) | [TestSpec](tests/TestSpecs/COMPGEN-DOTPLOT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_GenerateDotPlot_Tests.cs) |
| ☑ | MOTIF-DISCOVER-001 | Matching | 1 | [Evidence](docs/Evidence/MOTIF-DISCOVER-001-Evidence.md) | [TestSpec](tests/TestSpecs/MOTIF-DISCOVER-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_DiscoverMotifs_Tests.cs) |
| ☑ | MOTIF-SHARED-001 | Matching | 1 | [Evidence](docs/Evidence/MOTIF-SHARED-001-Evidence.md) | [TestSpec](tests/TestSpecs/MOTIF-SHARED-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_FindSharedMotifs_Tests.cs) |
| ☑ | MOTIF-REGULATORY-001 | Matching | 1 | [Evidence](docs/Evidence/MOTIF-REGULATORY-001-Evidence.md) | [TestSpec](tests/TestSpecs/MOTIF-REGULATORY-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_FindRegulatoryElements_Tests.cs) |
| ☑ | MOTIF-GENERATE-001 | Matching | 1 | [Evidence](docs/Evidence/MOTIF-GENERATE-001-Evidence.md) | [TestSpec](tests/TestSpecs/MOTIF-GENERATE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_GenerateConsensus_Tests.cs) |
| ☑ | KMER-ASYNC-001 | K-mer | 1 | [Evidence](docs/Evidence/KMER-ASYNC-001-Evidence.md) | [TestSpec](tests/TestSpecs/KMER-ASYNC-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_CountKmersAsync_Tests.cs) |
| ☑ | KMER-UNIQUE-001 | K-mer | 1 | [Evidence](docs/Evidence/KMER-UNIQUE-001-Evidence.md) | [TestSpec](tests/TestSpecs/KMER-UNIQUE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_FindUniqueAndMinCount_Tests.cs) |
| ☑ | KMER-GENERATE-001 | K-mer | 1 | [Evidence](docs/Evidence/KMER-GENERATE-001-Evidence.md) | [TestSpec](tests/TestSpecs/KMER-GENERATE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_GenerateAllKmers_Tests.cs) |
| ☑ | KMER-BOTH-001 | K-mer | 1 | [Evidence](docs/Evidence/KMER-BOTH-001-Evidence.md) | [TestSpec](tests/TestSpecs/KMER-BOTH-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_CountKmersBothStrands_Tests.cs) |
| ☑ | KMER-STATS-001 | K-mer | 1 | [Evidence](docs/Evidence/KMER-STATS-001-Evidence.md) | [TestSpec](tests/TestSpecs/KMER-STATS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_AnalyzeKmers_Tests.cs) |
| ☑ | KMER-POSITIONS-001 | K-mer | 1 | [Evidence](docs/Evidence/KMER-POSITIONS-001-Evidence.md) | [TestSpec](tests/TestSpecs/KMER-POSITIONS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_FindKmerPositions_Tests.cs) |
| ☑ | SEQ-ATSKEW-001 | Composition | 1 | [Evidence](docs/Evidence/SEQ-ATSKEW-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-ATSKEW-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculator_CalculateAtSkew_Tests.cs) |
| ☑ | SEQ-REPLICATION-001 | Composition | 1 | [Evidence](docs/Evidence/SEQ-REPLICATION-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-REPLICATION-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculator_PredictReplicationOrigin_Tests.cs) |
| ☑ | SEQ-GC-ANALYSIS-001 | Composition | 1 | [Evidence](docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-GC-ANALYSIS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculator_AnalyzeGcContent_Tests.cs) |
| ☑ | DISORDER-MORF-001 | ProteinPred | 1 | [Evidence](docs/Evidence/DISORDER-MORF-001-Evidence.md) | [TestSpec](tests/TestSpecs/DISORDER-MORF-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_MoRF_Tests.cs) |
| ☑ | DISORDER-PROPENSITY-001 | ProteinPred | 3 | [Evidence](docs/Evidence/DISORDER-PROPENSITY-001-Evidence.md) | [TestSpec](tests/TestSpecs/DISORDER-PROPENSITY-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_GetDisorderPropensity_Tests.cs) |
| ☑ | DISORDER-LC-001 | ProteinPred | 1 | [Evidence](docs/Evidence/DISORDER-LC-001-Evidence.md) | [TestSpec](tests/TestSpecs/DISORDER-LC-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_LowComplexity_Tests.cs) |
| ☑ | SEQ-COMPOSITION-001 | Statistics | 2 | [Evidence](docs/Evidence/SEQ-COMPOSITION-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-COMPOSITION-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateNucleotideComposition_Tests.cs) |
| ☑ | SEQ-TM-001 | Statistics | 2 | [Evidence](docs/Evidence/SEQ-TM-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-TM-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateThermodynamics_Tests.cs) |
| ☑ | SEQ-ENTROPY-PROFILE-001 | Statistics | 1 | docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md | tests/TestSpecs/SEQ-ENTROPY-PROFILE-001.md | tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateEntropyProfile_Tests.cs |
| ☑ | SEQ-GC-PROFILE-001 | Statistics | 1 | [Evidence](docs/Evidence/SEQ-GC-PROFILE-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-GC-PROFILE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateGcContentProfile_Tests.cs) |
| ☑ | SEQ-CODON-FREQ-001 | Statistics | 1 | [Evidence](docs/Evidence/SEQ-CODON-FREQ-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-CODON-FREQ-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateCodonFrequencies_Tests.cs) |
| ☑ | SEQ-SUMMARY-001 | Statistics | 1 | [Evidence](docs/Evidence/SEQ-SUMMARY-001-Evidence.md) | [TestSpec](tests/TestSpecs/SEQ-SUMMARY-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_SummarizeNucleotideSequence_Tests.cs) |
| ☑ | GENOMIC-TANDEM-001 | Analysis | 1 | [Evidence](docs/Evidence/GENOMIC-TANDEM-001-Evidence.md) | [TestSpec](tests/TestSpecs/GENOMIC-TANDEM-001.md) | GenomicAnalyzer_TandemRepeat_Tests.cs (consolidated into REP-TANDEM-001) |
| ☑ | GENOMIC-SIMILARITY-001 | Analysis | 1 | [Evidence](docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md) | [TestSpec](tests/TestSpecs/GENOMIC-SIMILARITY-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_CalculateSimilarity_Tests.cs) |
| ☑ | GENOMIC-ORF-001 | Analysis | 1 | [Evidence](docs/Evidence/GENOMIC-ORF-001-Evidence.md) | [TestSpec](tests/TestSpecs/GENOMIC-ORF-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindOpenReadingFrames_Tests.cs) |
| ☑ | ONCO-SOMATIC-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-SOMATIC-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-SOMATIC-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CallSomaticMutations_Tests.cs) |
| ☑ | ONCO-VAF-001 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-VAF-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-VAF-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateVAF_Tests.cs) |
| ☑ | ONCO-DRIVER-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-DRIVER-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-DRIVER-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_IdentifyDriverMutations_Tests.cs) |
| ☑ | ONCO-ARTIFACT-001 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-ARTIFACT-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-ARTIFACT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_FilterArtifacts_Tests.cs) |
| ☑ | ONCO-ANNOT-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-ANNOT-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-ANNOT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AnnotateCancerVariants_Tests.cs) |
| ☑ | ONCO-TMB-001 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-TMB-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-TMB-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateTMB_Tests.cs) |
| ☑ | ONCO-MSI-001 | Oncology | 3 | Niu et al. (2014) MSIsensor, niu-lab/msisensor2, Boland et al. (1998) | [ONCO-MSI-001.md](TestSpecs/ONCO-MSI-001.md) | OncologyAnalyzer_DetectMSI_Tests.cs |
| ☑ | ONCO-HRD-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-HRD-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-HRD-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CalculateHRDScore_Tests.cs) |
| ☑ | ONCO-LOH-001 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-LOH-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-LOH-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectLOH_Tests.cs) |
| ☑ | ONCO-SIG-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-SIG-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-SIG-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifySbsContext_Tests.cs) |
| ☑ | ONCO-SIG-002 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-SIG-002-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-SIG-002.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_FitSignatures_Tests.cs) |
| ☑ | ONCO-SIG-003 | Oncology | 1 | [Evidence](docs/Evidence/ONCO-SIG-003-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-SIG-003.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_BootstrapExposures_Tests.cs) |
| ☑ | ONCO-SIG-004 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-SIG-004-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-SIG-004.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyMutationalProcess_Tests.cs) |
| ☑ | ONCO-FUSION-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-FUSION-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-FUSION-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectFusions_Tests.cs) |
| ☑ | ONCO-FUSION-002 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-FUSION-002-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-FUSION-002.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_MatchKnownFusions_Tests.cs) |
| ☑ | ONCO-FUSION-003 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-FUSION-003-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-FUSION-003.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AnalyzeBreakpoint_Tests.cs) |
| ☑ | ONCO-CNA-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-CNA-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-CNA-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CopyNumberClassification_Tests.cs) |
| ☑ | ONCO-CNA-002 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-CNA-002-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-CNA-002.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectFocalAmplifications_Tests.cs) |
| ☑ | ONCO-CNA-003 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-CNA-003-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-CNA-003.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectHomozygousDeletions_Tests.cs) |
| ☑ | ONCO-PURITY-001 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-PURITY-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-PURITY-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_EstimatePurity_Tests.cs) |
| ☑ | ONCO-PLOIDY-001 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-PLOIDY-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-PLOIDY-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_EstimatePloidy_Tests.cs) |
| ☑ | ONCO-CLONAL-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-CLONAL-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-CLONAL-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyClonality_Tests.cs) |
| ☑ | ONCO-NEO-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-NEO-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-NEO-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_GenerateNeoantigenPeptides_Tests.cs) |
| ☑ | ONCO-MHC-001 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-MHC-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-MHC-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyMhcBinding_Tests.cs) |
| ☑ | ONCO-IMMUNE-001 | Oncology | 3 | 40 | [Evidence](docs/Evidence/ONCO-IMMUNE-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-IMMUNE-001.md) |
| ☑ | ONCO-CTDNA-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-CTDNA-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-CTDNA-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CtDnaAnalysis_Tests.cs) |
| ☑ | ONCO-MRD-001 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-MRD-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-MRD-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectMRD_Tests.cs) |
| ☑ | ONCO-CHIP-001 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-CHIP-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-CHIP-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_FilterCHIP_Tests.cs) |
| ☑ | ONCO-PHYLO-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-PHYLO-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-PHYLO-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ReconstructPhylogeny_Tests.cs) |
| ☑ | ONCO-CCF-001 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-CCF-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-CCF-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_EstimateCcf_Tests.cs) |
| ☑ | ONCO-HETERO-001 | Oncology | 2 | [Evidence](docs/Evidence/ONCO-HETERO-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-HETERO-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AnalyzeHeterogeneity_Tests.cs) |
| ☑ | ONCO-HLA-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-HLA-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-HLA-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_HlaAnalysis_Tests.cs) |
| ☑ | ONCO-ACTION-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-ACTION-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-ACTION-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AssessActionability_Tests.cs) |
| ☑ | ONCO-SV-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-SV-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-SV-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyComplexRearrangement_Tests.cs) |
| ☑ | ONCO-EXPR-001 | Oncology | 3 | [Evidence](docs/Evidence/ONCO-EXPR-001-Evidence.md) | [TestSpec](tests/TestSpecs/ONCO-EXPR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_IdentifyOutlierGenes_Tests.cs) |
| ☑ | RNA-ACCESS-001 | RnaStructure | 2 | McCaskill (1990), RNAplfold/Bernhart (2006) | [TestSpec](tests/TestSpecs/RNA-ACCESS-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_UnpairedProbabilities_Tests.cs) |
| ☑ | PROTMOTIF-HMM-001 | ProteinMotif | 4 | Eddy (1998, 2011) HMMER/Plan7, Pfam | [TestSpec](tests/TestSpecs/PROTMOTIF-HMM-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindDomainsByHmm_Tests.cs) |
| ☑ | PRIMER-NNTM-001 | MolTools | 2 | SantaLucia (1998), Allawi & SantaLucia (1997), Owczarzy (2004/2008), Bommarito (2000) | [TestSpec](tests/TestSpecs/PRIMER-NNTM-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_NearestNeighborTm_Tests.cs) |
| ☑ | PRIMER-HAIRPIN-001 | MolTools | 2 | SantaLucia (1998), SantaLucia & Hicks (2004), UNAFold DNA params | [TestSpec](tests/TestSpecs/PRIMER-HAIRPIN-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_HairpinTm_Tests.cs) |
| ☑ | PRIMER-DIMER-001 | MolTools | 3 | SantaLucia & Hicks (2004), primer3 ntthal | [TestSpec](tests/TestSpecs/PRIMER-DIMER-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_DimerTm_Tests.cs) |
| ☑ | PROBE-LNATM-001 | MolTools | 2 | McTigue (2004) LNA NN, Kutyavin (2000) MGB | [TestSpec](tests/TestSpecs/PROBE-LNATM-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_LnaTm_Tests.cs) |
| ☑ | PROBE-EVALUE-001 | MolTools | 2 | Karlin & Altschul (1990), Altschul et al. (1990) BLAST | [TestSpec](tests/TestSpecs/PROBE-EVALUE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_ProbeValidation_Tests.cs) |
| ☑ | MHC-NN-001 | Oncology | 3 | O'Donnell et al. (2018, 2020) MHCflurry, Apache-2.0 | [TestSpec](tests/TestSpecs/MHC-NN-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MhcflurryAffinityPredictor_PredictIc50_Tests.cs) |
| ☑ | MHC-MATRIX-001 | Oncology | 3 | Peters et al. (2005) SMM, Parker et al. (1994) BIMAS | [TestSpec](tests/TestSpecs/MHC-MATRIX-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyMhcBinding_Tests.cs) |
| ☑ | IMMUNE-NUSVR-001 | Oncology | 3 | Newman et al. (2015) CIBERSORT, Schölkopf et al. (2000) ν-SVR, Monaco et al. (2019) ABIS CC-BY | [TestSpec](tests/TestSpecs/IMMUNE-NUSVR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ImmuneAnalyzer_ImmuneInfiltration_Tests.cs) |
| ☑ | META-CHECKM-001 | Metagenomics | 4 | Parks et al. (2015) CheckM, Parks et al. (2018) GTDB, Pfam CC0 | [TestSpec](tests/TestSpecs/META-CHECKM-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_MarkerGeneQuality_Tests.cs) |
| ☑ | META-TETRA-001 | Metagenomics | 2 | Teeling et al. (2004) TETRA, Schbath (1995) | [TestSpec](tests/TestSpecs/META-TETRA-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_TetranucleotideZScore_Tests.cs) |
| ☑ | SPLICE-MAXENT3-001 | Splicing | 1 | Yeo & Burge (2004), maxentpy (MIT) | [TestSpec](tests/TestSpecs/SPLICE-MAXENT3-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_AcceptorSite_Tests.cs) |
| ☑ | SPLICE-MAXENT5-001 | Splicing | 1 | Yeo & Burge (2004), maxentpy (MIT) | [TestSpec](tests/TestSpecs/SPLICE-MAXENT5-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs) |
| ☑ | MIRNA-CONTEXT-001 | MiRNA | 2 | Agarwal et al. (2015) TargetScan context++ | [TestSpec](tests/TestSpecs/MIRNA-CONTEXT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_TargetPrediction_Tests.cs) |
| ☑ | MIRNA-PCT-001 | MiRNA | 1 | Friedman et al. (2009) PCT, TargetScan | [TestSpec](tests/TestSpecs/MIRNA-PCT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_TargetPrediction_Tests.cs) |
| ☑ | MIRNA-CLASSIFY-001 | MiRNA | 1 | Bonnet et al. (2004), miRBase (public domain), Zhang (2006) MFEI | [TestSpec](tests/TestSpecs/MIRNA-CLASSIFY-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs) |
| ☑ | MIRNA-CLEAVAGE-001 | MiRNA | 1 | Han et al. (2006), Park et al. (2011), Auyeung et al. (2013) | [TestSpec](tests/TestSpecs/MIRNA-CLEAVAGE-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs) |
| ☑ | REP-APPROX-001 | Repeats | 2 | Benson (1999) TRF | [TestSpec](tests/TestSpecs/REP-APPROX-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_ApproximateTandemRepeats_Tests.cs) |
| ☑ | CHROM-ALPHASAT-001 | Chromosome | 2 | Waye & Willard (1987), Henikoff et al. (2001), CENP-B box motif | [TestSpec](tests/TestSpecs/CHROM-ALPHASAT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_AlphaSatellite_Tests.cs) |
| ☑ | CHROM-HOR-001 | Chromosome | 1 | McNulty & Sullivan (2018), Alkan et al. (2007) | [TestSpec](tests/TestSpecs/CHROM-HOR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_HigherOrderRepeat_Tests.cs) |

**Statuses:** ☐ Not Started | ⏳ In Progress | ☑ Complete | ⛔ Blocked

---

## Definition of Done (DoD)

### Required Criteria

| # | Criterion | Artifact |
|---|-----------|----------|
| 1 | TestSpec created | `TestSpecs/{TestUnitID}.md` |
| 2 | Tests written | `*.Tests/{Class}_{Method}_Tests.cs` |
| 3 | Branch coverage ≥ 80% | Coverage report |
| 4 | Edge cases covered | null, empty, boundary, error |
| 5 | Tests pass | CI green |
| 6 | Evidence documented | PR/commit link in Registry |

### Test Quality Criteria

- [ ] Tests are independent (order-independent)
- [ ] Tests are deterministic (no random without seed)
- [ ] Naming: `Method_Scenario_ExpectedResult`
- [ ] Structure: Arrange-Act-Assert
- [ ] One assert per logical check

### For O(n²) and Higher Algorithms

- [ ] Property-based test for invariant
- [ ] Performance baseline recorded

---

## Test Units by Area

### 1. Sequence Composition (4 units)

#### SEQ-GC-001: GC Content Calculation

| Field | Value |
|-------|-------|
| **Canonical** | `SequenceExtensions.CalculateGcContent(ReadOnlySpan<char>)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ result ≤ 100 (percentage) or 0 ≤ result ≤ 1 (fraction) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateGcContent(ReadOnlySpan<char>)` | SequenceExtensions | Canonical |
| `CalculateGcFraction(ReadOnlySpan<char>)` | SequenceExtensions | Variant (0-1) |
| `CalculateGcContentFast(string)` | SequenceExtensions | Delegate |
| `CalculateGcFractionFast(string)` | SequenceExtensions | Delegate |
| `GcContent` (property) | DnaSequence | Delegate |

**Edge Cases:**
- [ ] Empty sequence → 0
- [ ] All G/C → 100% / 1.0
- [ ] All A/T → 0% / 0.0
- [ ] Mixed case input
- [ ] Non-ACGT characters (N, etc.)

---

#### SEQ-COMP-001: DNA Complement

| Field | Value |
|-------|-------|
| **Canonical** | `SequenceExtensions.GetComplementBase(char)` |
| **Complexity** | O(n) for sequence |
| **Invariant** | Complement(Complement(x)) = x |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetComplementBase(char)` | SequenceExtensions | Canonical |
| `TryGetComplement(ReadOnlySpan, Span)` | SequenceExtensions | Span API |
| `Complement()` | DnaSequence | Instance |

**Edge Cases:**
- [ ] A ↔ T, G ↔ C (both directions)
- [ ] Case insensitivity (a → T)
- [ ] RNA support (U → A)
- [ ] Unknown base → unchanged
- [ ] Destination too small → false

---

#### SEQ-REVCOMP-001: Reverse Complement

| Field | Value |
|-------|-------|
| **Canonical** | `SequenceExtensions.TryGetReverseComplement(ReadOnlySpan, Span)` |
| **Complexity** | O(n) |
| **Invariant** | RevComp(RevComp(x)) = x |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `TryGetReverseComplement(ReadOnlySpan, Span)` | SequenceExtensions | Canonical |
| `ReverseComplement()` | DnaSequence | Instance |
| `GetReverseComplementString(string)` | DnaSequence | Static helper |
| `TryWriteReverseComplement(Span)` | DnaSequence | Span API |

**Edge Cases:**
- [ ] Empty sequence
- [ ] Single nucleotide
- [ ] Palindrome (self-complementary)
- [ ] Destination too small → false

---

#### SEQ-VALID-001: Sequence Validation

| Field | Value |
|-------|-------|
| **Canonical** | `SequenceExtensions.IsValidDna(ReadOnlySpan<char>)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `IsValidDna(ReadOnlySpan<char>)` | SequenceExtensions | Canonical DNA |
| `IsValidRna(ReadOnlySpan<char>)` | SequenceExtensions | Canonical RNA |
| `TryCreate(string, out DnaSequence)` | DnaSequence | Factory |
| `IsValid` (property) | DnaSequence | Instance |

**Edge Cases:**
- [ ] Empty → true (or false? define!)
- [ ] All valid bases
- [ ] Single invalid character
- [ ] Lowercase valid
- [ ] Whitespace handling

---

### 2. Pattern Matching (5 units)

#### PAT-EXACT-001: Exact Pattern Search

| Field | Value |
|------|----------|
| **Canonical** | `SuffixTree.FindAllOccurrences(string)` |
| **Complexity** | O(m + k) where m=pattern length, k=occurrences |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindAllOccurrences(string)` | SuffixTree/DnaSequence | Canonical |
| `Contains(string)` | SuffixTree/DnaSequence | Existence check |
| `CountOccurrences(string)` | SuffixTree/DnaSequence | Count only |
| `FindMotif(DnaSequence, string)` | GenomicAnalyzer | Wrapper |
| `FindExactMotif(DnaSequence, string)` | MotifFinder | Wrapper |

**Edge Cases:**
- [ ] Pattern not found → empty
- [ ] Pattern = entire sequence
- [ ] Overlapping occurrences
- [ ] Empty pattern
- [ ] Pattern longer than sequence

---

#### PAT-APPROX-001: Approximate Matching (Hamming)

| Field | Value |
|------|----------|
| **Canonical** | `ApproximateMatcher.FindWithMismatches(...)` |
| **Complexity** | O(n × m) |
| **Invariant** | HammingDistance requires equal length |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindWithMismatches(seq, pattern, max)` | ApproximateMatcher | Canonical |
| `HammingDistance(s1, s2)` | ApproximateMatcher | Distance |
| `HammingDistance(span1, span2)` | SequenceExtensions | Span API |

**Edge Cases:**
- [ ] maxMismatches = 0 → exact match
- [ ] maxMismatches ≥ pattern length → all positions
- [ ] Unequal lengths for HammingDistance → exception

---

#### PAT-APPROX-002: Approximate Matching (Edit Distance)

| Field | Value |
|------|----------|
| **Canonical** | `ApproximateMatcher.FindWithEdits(...)` |
| **Complexity** | O(n × m²) |
| **Invariant** | EditDistance(a,b) = EditDistance(b,a) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindWithEdits(seq, pattern, maxEdits)` | ApproximateMatcher | Canonical |
| `EditDistance(s1, s2)` | ApproximateMatcher | Distance |

**Edge Cases:**
- [ ] Identical strings → 0
- [ ] One empty string → length of other
- [ ] Single character difference
- [ ] Insertion vs deletion vs substitution

---

#### PAT-IUPAC-001: IUPAC Degenerate Matching

| Field | Value |
|------|----------|
| **Canonical** | `MotifFinder.FindDegenerateMotif(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindDegenerateMotif(seq, motif)` | MotifFinder | Canonical |
| `MatchesIupac(nucleotide, iupacCode)` | IupacHelper | Helper |

**Edge Cases:**
- [ ] All IUPAC codes: R, Y, S, W, K, M, B, D, H, V, N
- [ ] Mixed standard + IUPAC
- [ ] N matches everything

---

#### PAT-PWM-001: Position Weight Matrix

| Field | Value |
|------|----------|
| **Canonical** | `MotifFinder.ScanWithPwm(...)` |
| **Complexity** | O(n × m) where m=PWM width |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CreatePwm(alignedSequences)` | MotifFinder | Construction |
| `ScanWithPwm(seq, pwm, threshold)` | MotifFinder | Scanning |

**Edge Cases:**
- [ ] Empty alignment
- [ ] Single sequence alignment
- [ ] Threshold at boundary
- [ ] All positions below threshold

---

### 3. Repeat Analysis (5 units)

#### REP-STR-001: Microsatellite Detection (STR)

| Field | Value |
|------|----------|
| **Canonical** | `RepeatFinder.FindMicrosatellites(...)` |
| **Complexity** | O(n × U × R) where U=maxUnitLength, R=maxRepeats |
| **Invariant** | Result positions are non-overlapping (or document overlap policy) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindMicrosatellites(DnaSequence, ...)` | RepeatFinder | Canonical |
| `FindMicrosatellites(string, ...)` | RepeatFinder | Overload |
| `FindMicrosatellites(..., CancellationToken)` | RepeatFinder | Cancellable |
| `FindMicrosatellites(..., IProgress)` | RepeatFinder | With progress |
| `FindApproximateTandemRepeats(...)` | RepeatFinder | Approximate (TRF, opt-in) |

**Edge Cases:**
- [ ] No repeats found
- [ ] Entire sequence is one repeat
- [ ] minRepeats = 2 (minimum)
- [ ] Unit length 1-6 (mono to hexa)
- [ ] Cancellation mid-operation
- [ ] Approximate: imperfect (one substitution) / interrupted tract reported as one repeat
- [ ] Approximate: single indel → exact percent-indels; Minscore threshold

---

#### REP-TANDEM-001: Tandem Repeat Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.FindTandemRepeats(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindTandemRepeats(seq, minUnitLength, minRepetitions)` | GenomicAnalyzer | Canonical |
| `GetTandemRepeatSummary(seq, minRepeats)` | RepeatFinder | Summary |

---

#### REP-INV-001: Inverted Repeat Detection

| Field | Value |
|------|----------|
| **Canonical** | `RepeatFinder.FindInvertedRepeats(...)` |
| **Complexity** | O(n² × L) where L=maxLoopLength |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindInvertedRepeats(seq, minArm, maxLoop)` | RepeatFinder | Canonical |

**Edge Cases:**
- [ ] Perfect palindrome (loop = 0)
- [ ] Maximum loop length
- [ ] Arm length at boundaries

---

#### REP-DIRECT-001: Direct Repeat Detection

| Field | Value |
|------|----------|
| **Canonical** | `RepeatFinder.FindDirectRepeats(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindDirectRepeats(seq, minLen, maxLen, minSpacing)` | RepeatFinder | Canonical |

---

#### REP-PALIN-001: Palindrome Detection

| Field | Value |
|------|----------|
| **Canonical** | `RepeatFinder.FindPalindromes(...)` |
| **Complexity** | O(n²) |
| **Invariant** | Palindrome = reverse complement equals self |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindPalindromes(seq, minLen, maxLen)` | RepeatFinder | Canonical |
| `FindPalindromes(seq, minLen, maxLen)` | GenomicAnalyzer | Alternate |

---

### 4. Molecular Tools (8 units)

#### CRISPR-PAM-001: PAM Site Detection

| Field | Value |
|------|----------|
| **Canonical** | `CrisprDesigner.FindPamSites(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindPamSites(seq, systemType)` | CrisprDesigner | Canonical |
| `GetSystem(systemType)` | CrisprDesigner | System info |

**CRISPR Systems:**
- [ ] SpCas9 (NGG)
- [ ] SpCas9-NAG
- [ ] SaCas9 (NNGRRT)
- [ ] Cas12a (TTTV)
- [ ] AsCas12a, LbCas12a
- [ ] CasX

---

#### CRISPR-GUIDE-001: Guide RNA Design

| Field | Value |
|------|----------|
| **Canonical** | `CrisprDesigner.DesignGuideRnas(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DesignGuideRnas(seq, start, end, type)` | CrisprDesigner | Canonical |
| `EvaluateGuideRna(guide, type, params)` | CrisprDesigner | Scoring |

---

#### CRISPR-OFF-001: Off-Target Analysis

| Field | Value |
|------|----------|
| **Canonical** | `CrisprDesigner.FindOffTargets(...)` |
| **Complexity** | O(n × m) with maxMismatches; may be higher |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindOffTargets(guide, genome, maxMismatches)` | CrisprDesigner | Canonical |
| `CalculateSpecificityScore(guide, genome, type)` | CrisprDesigner | Score |

---

#### PRIMER-TM-001: Melting Temperature

| Field | Value |
|------|----------|
| **Canonical** | `PrimerDesigner.CalculateMeltingTemperature(...)` |
| **Complexity** | O(n) |
| **Formula** | Wallace rule (<14bp), Marmur-Doty (≥14bp) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateMeltingTemperature(primer)` | PrimerDesigner | Canonical |
| `CalculateMeltingTemperatureWithSalt(primer, Na)` | PrimerDesigner | Salt corrected |

**Constants (ThermoConstants):**
- `WallaceMaxLength` = 14
- `CalculateWallaceTm(at, gc)` = 2×AT + 4×GC
- `CalculateMarmurDotyTm(gc, len)` = 64.9 + 41×(GC-16.4)/len
- `CalculateSaltCorrection(Na)` = 16.6 × log10(Na/1000)

**Edge Cases:**
- [ ] Empty primer → 0
- [ ] Short primer (<14) uses Wallace
- [ ] Long primer (≥14) uses Marmur-Doty
- [ ] Salt concentration variations

---

#### PRIMER-DESIGN-001: Primer Pair Design

| Field | Value |
|------|----------|
| **Canonical** | `PrimerDesigner.DesignPrimers(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DesignPrimers(template, start, end, params)` | PrimerDesigner | Canonical |
| `EvaluatePrimer(seq, pos, isForward, params)` | PrimerDesigner | Single primer |
| `GeneratePrimerCandidates(template, region)` | PrimerDesigner | All candidates |

**Parameters (PrimerParameters):**
- MinLength=18, MaxLength=25, OptimalLength=20
- MinGcContent=40, MaxGcContent=60
- MinTm=55, MaxTm=65, OptimalTm=60
- MaxHomopolymer=4, MaxDinucleotideRepeats=4

---

#### PRIMER-STRUCT-001: Primer Structure Analysis

| Field | Value |
|------|----------|
| **Canonical** | `PrimerDesigner.HasHairpinPotential(...)` |
| **Complexity** | O(m²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `HasHairpinPotential(seq, minStemLength)` | PrimerDesigner | Hairpin |
| `HasPrimerDimer(primer1, primer2, minComp)` | PrimerDesigner | Dimer |
| `Calculate3PrimeStability(seq)` | PrimerDesigner | ΔG calculation |
| `FindLongestHomopolymer(seq)` | PrimerDesigner | Structure |
| `FindLongestDinucleotideRepeat(seq)` | PrimerDesigner | Structure |

---

#### RESTR-FIND-001: Restriction Site Detection

| Field | Value |
|------|----------|
| **Canonical** | `RestrictionAnalyzer.FindSites(...)` |
| **Complexity** | O(n × k) where k=enzymes |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSites(seq, enzymeNames)` | RestrictionAnalyzer | Canonical |
| `FindAllSites(seq)` | RestrictionAnalyzer | All 40+ enzymes |
| `GetEnzyme(name)` | RestrictionAnalyzer | Lookup |

**Enzyme Database:** 40+ enzymes (EcoRI, BamHI, HindIII, NotI, etc.)

---

#### RESTR-DIGEST-001: Digest Simulation

| Field | Value |
|------|----------|
| **Canonical** | `RestrictionAnalyzer.Digest(...)` |
| **Complexity** | O(n + k log k) where k=cut sites |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Digest(seq, enzymeNames)` | RestrictionAnalyzer | Canonical |
| `GetDigestSummary(seq, enzymeNames)` | RestrictionAnalyzer | Summary |
| `CreateMap(seq, enzymeNames)` | RestrictionAnalyzer | Full map |

---

### 5. Genome Annotation (4 units)

#### ANNOT-ORF-001: ORF Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAnnotator.FindOrfs(...)` |
| **Complexity** | O(n) |
| **Invariant** | ORF starts with start codon, ends with stop codon |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindOrfs(dna, minLen, bothStrands, requireStart)` | GenomeAnnotator | Canonical |
| `FindLongestOrfsPerFrame(dna, bothStrands)` | GenomeAnnotator | Per-frame |
| `FindOpenReadingFrames(seq, minLen)` | GenomicAnalyzer | Alternate |

**Edge Cases:**
- [ ] No ORF found
- [ ] ORF extends to sequence end (no stop)
- [ ] All 6 reading frames
- [ ] Overlapping ORFs

---

#### ANNOT-GENE-001: Gene Prediction

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAnnotator.PredictGenes(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictGenes(dna, minOrfLen, prefix)` | GenomeAnnotator | Canonical |
| `FindRibosomeBindingSites(dna, window)` | GenomeAnnotator | RBS/SD |

---

#### ANNOT-PROM-001: Promoter Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAnnotator.FindPromoterMotifs(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindPromoterMotifs(dna)` | GenomeAnnotator | Canonical |

**Motifs:** -35 box (TTGACA), -10 box (TATAAT)

---

#### ANNOT-GFF-001: GFF3 I/O

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAnnotator.ParseGff3(...)` / `ToGff3(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ParseGff3(lines)` | GenomeAnnotator | Parse |
| `ToGff3(annotations, seqId)` | GenomeAnnotator | Export |

---

### 6. K-mer Analysis (3 units)

#### KMER-COUNT-001: K-mer Counting

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.CountKmers(...)` |
| **Complexity** | O(n) |
| **Invariant** | Sum of counts = n - k + 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CountKmers(seq, k)` | KmerAnalyzer | Canonical |
| `CountKmersSpan(seq, k)` | SequenceExtensions | Span API |
| `CountKmersBothStrands(dnaSeq, k)` | KmerAnalyzer | Both strands |

---

#### KMER-FREQ-001: K-mer Frequency Analysis

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.GetKmerSpectrum(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetKmerSpectrum(seq, k)` | KmerAnalyzer | Spectrum |
| `GetKmerFrequencies(seq, k)` | KmerAnalyzer | Normalized |
| `CalculateKmerEntropy(seq, k)` | KmerAnalyzer | Entropy |

---

#### KMER-FIND-001: K-mer Search

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.FindMostFrequentKmers(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindMostFrequentKmers(seq, k)` | KmerAnalyzer | Most frequent |
| `FindUniqueKmers(seq, k)` | KmerAnalyzer | Count = 1 |
| `FindClumps(seq, k, window, minOcc)` | KmerAnalyzer | Clumps |

---

### 7. Alignment (4 units)

#### ALIGN-GLOBAL-001: Global Alignment (Needleman-Wunsch)

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAligner.GlobalAlign(...)` |
| **Complexity** | O(n × m) |
| **Invariant** | Optimal global alignment score |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GlobalAlign(seq1, seq2, scoring)` | SequenceAligner | Canonical |

---

#### ALIGN-LOCAL-001: Local Alignment (Smith-Waterman)

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAligner.LocalAlign(...)` |
| **Complexity** | O(n × m) |
| **Invariant** | Score ≥ 0, finds best local match |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `LocalAlign(seq1, seq2, scoring)` | SequenceAligner | Canonical |

---

#### ALIGN-SEMI-001: Semi-Global Alignment

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAligner.SemiGlobalAlign(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `SemiGlobalAlign(seq1, seq2, scoring)` | SequenceAligner | Canonical |

---

#### ALIGN-MULTI-001: Multiple Sequence Alignment

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAligner.MultipleAlign(...)` |
| **Complexity** | O(n² × m) progressive alignment |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `MultipleAlign(sequences)` | SequenceAligner | Canonical |

---

### 8. Phylogenetics (4 units)

#### PHYLO-DIST-001: Distance Matrix

| Field | Value |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.CalculateDistanceMatrix(...)` |
| **Complexity** | O(n² × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateDistanceMatrix(seqs, method)` | PhylogeneticAnalyzer | Canonical |
| `CalculatePairwiseDistance(s1, s2, method)` | PhylogeneticAnalyzer | Single pair |

**Distance Methods:** p-distance, Jukes-Cantor, Kimura 2-parameter, Hamming

---

#### PHYLO-TREE-001: Tree Construction

| Field | Value |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.BuildTree(...)` |
| **Complexity** | O(n³) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `BuildTree(seqs, distMethod, treeMethod)` | PhylogeneticAnalyzer | Canonical |

**Tree Methods:** UPGMA, Neighbor-Joining

---

#### PHYLO-NEWICK-001: Newick I/O

| Field | Value |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.ToNewick(...)` / `ParseNewick(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ToNewick(treeNode)` | PhylogeneticAnalyzer | Export |
| `ParseNewick(newickString)` | PhylogeneticAnalyzer | Parse |

---

#### PHYLO-COMP-001: Tree Comparison

| Field | Value |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.RobinsonFouldsDistance(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `RobinsonFouldsDistance(tree1, tree2)` | PhylogeneticAnalyzer | Topology |
| `FindMRCA(root, taxon1, taxon2)` | PhylogeneticAnalyzer | MRCA |
| `PatristicDistance(root, t1, t2)` | PhylogeneticAnalyzer | Path length |

---

### 9. Population Genetics (5 units)

#### POP-FREQ-001: Allele Frequencies

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(...)` |
| **Complexity** | O(n) |
| **Invariant** | p + q = 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateAlleleFrequencies(hom_maj, het, hom_min)` | PopulationGeneticsAnalyzer | Canonical |
| `CalculateMAF(genotypes)` | PopulationGeneticsAnalyzer | MAF |
| `FilterByMAF(variants, min, max)` | PopulationGeneticsAnalyzer | Filter |

---

#### POP-DIV-001: Diversity Statistics

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateDiversityStatistics(...)` |
| **Complexity** | O(n² × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateNucleotideDiversity(seqs)` | PopulationGeneticsAnalyzer | π |
| `CalculateWattersonTheta(segSites, n)` | PopulationGeneticsAnalyzer | θ |
| `CalculateTajimasD(kHat, S, n)` | PopulationGeneticsAnalyzer | D |
| `CalculateObservedHeterozygosity(seqs)` | PopulationGeneticsAnalyzer | H_obs |
| `CalculateExpectedHeterozygosity(seqs)` | PopulationGeneticsAnalyzer | H_exp |
| `CalculateDiversityStatistics(seqs)` | PopulationGeneticsAnalyzer | All |

---

#### POP-HW-001: Hardy-Weinberg Test

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.TestHardyWeinberg(...)` |
| **Complexity** | O(1) per variant |
| **Invariant** | Expected: p², 2pq, q² |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `TestHardyWeinberg(variantId, counts)` | PopulationGeneticsAnalyzer | Canonical |

---

#### POP-FST-001: F-Statistics

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateFst(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ Fst ≤ 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateFst(pop1, pop2)` | PopulationGeneticsAnalyzer | Canonical |
| `CalculateFStatistics(variantData)` | PopulationGeneticsAnalyzer | Fis, Fit, Fst |

---

#### POP-LD-001: Linkage Disequilibrium

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateLD(...)` |
| **Complexity** | O(n) per pair |
| **Algorithm** | Squared Pearson correlation (r²), Lewontin D' |
| **Formula** | r² = Cov(X₁,X₂)² / (Var(X₁) × Var(X₂)); D = Cov/2; D' = \|D\|/D_max |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateLD(var1, var2, genotypes)` | PopulationGeneticsAnalyzer | D', r² |
| `FindHaplotypeBlocks(variants)` | PopulationGeneticsAnalyzer | Blocks |

**References:** Wikipedia (Linkage disequilibrium, LD for diploid frequencies), Hill & Robertson (1968), Lewontin (1964), Gabriel et al. (2002)

---

### 10. Chromosome Analysis (5 units)

#### CHROM-TELO-001: Telomere Analysis

| Field | Value |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.AnalyzeTelomeres(...)` |
| **Complexity** | O(n) |
| **Constant** | Human telomere repeat: TTAGGG |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeTelomeres(chrName, seq, repeat, ...)` | ChromosomeAnalyzer | Canonical |
| `EstimateTelomereLengthFromTSRatio(tsRatio)` | ChromosomeAnalyzer | qPCR estimate |

---

#### CHROM-CENT-001: Centromere Analysis

| Field | Value |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.AnalyzeCentromere(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeCentromere(chrName, seq, windowSize)` | ChromosomeAnalyzer | Canonical (generic repeat heuristic) |
| `DetectAlphaSatellite(seq)` | ChromosomeAnalyzer | Canonical (alpha-satellite-specific: 171-bp tandem + AT + CENP-B) |
| `FindCenpBBoxes(seq)` | ChromosomeAnalyzer | Canonical (CENP-B box `YTTCGTTGGAARCGGGA` scan) |
| `DetectHigherOrderRepeat(seq, monomerLength=171)` | ChromosomeAnalyzer | Canonical (alpha-satellite HOR structure: period, copy number, inter-/intra-HOR identity) |
| `AssignSuprachromosomalFamily(seq, reference=null)` | ChromosomeAnalyzer | Canonical (SF assignment vs bundled CC0 Dfam reference: SF3/SF4/SF5 + {SF1,SF2} via periodicity + A/B-box) |
| `LoadBundledAlphaSatelliteReference()` | ChromosomeAnalyzer | Bundled CC0 reference loader (Dfam ALR/ALRa/ALRb) |

---

#### CHROM-KARYO-001: Karyotype Analysis

| Field | Value |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.AnalyzeKaryotype(...)` |
| **Complexity** | O(k) where k=chromosomes |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeKaryotype(chromosomes, ploidy)` | ChromosomeAnalyzer | Canonical |
| `DetectPloidy(depths, expected)` | ChromosomeAnalyzer | Ploidy |

---

#### CHROM-ANEU-001: Aneuploidy Detection

| Field | Value |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.DetectAneuploidy(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectAneuploidy(depthData, medianDepth)` | ChromosomeAnalyzer | Canonical |
| `IdentifyWholeChromosomeAneuploidy(cnStates)` | ChromosomeAnalyzer | Classification |

---

#### CHROM-SYNT-001: Synteny Analysis

| Field | Value |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.FindSyntenyBlocks(...)` |
| **Complexity** | O(n log n) — requires verification |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSyntenyBlocks(orthologPairs, minGenes)` | ChromosomeAnalyzer | Canonical |
| `DetectRearrangements(syntenyBlocks)` | ChromosomeAnalyzer | Rearrangements |

---

### 11. Metagenomics (5 units)

#### META-CLASS-001: Taxonomic Classification

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.ClassifyReads(...)` |
| **Complexity** | O(n × m) where n=reads, m=read length |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ClassifyReads(reads, kmerDB, k)` | MetagenomicsAnalyzer | Canonical |
| `BuildKmerDatabase(refGenomes, k)` | MetagenomicsAnalyzer | DB construction |

---

#### META-PROF-001: Taxonomic Profile

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.GenerateTaxonomicProfile(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GenerateTaxonomicProfile(classifications)` | MetagenomicsAnalyzer | Canonical |

---

#### META-ALPHA-001: Alpha Diversity

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.CalculateAlphaDiversity(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateAlphaDiversity(abundances)` | MetagenomicsAnalyzer | Canonical |

**Indices:** Shannon, Simpson, Inverse Simpson, Chao1, Pielou's evenness

---

#### META-BETA-001: Beta Diversity

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.CalculateBetaDiversity(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateBetaDiversity(sample1, sample2)` | MetagenomicsAnalyzer | Canonical |

**Metrics:** Bray-Curtis, Jaccard

---

#### META-BIN-001: Genome Binning

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.BinContigs(...)` |
| **Complexity** | O(n × k × i) where k=bins, i=iterations — needs verification |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `BinContigs(contigs, numBins, minBinSize)` | MetagenomicsAnalyzer | Canonical |

---

### 12. Codon Optimization (4 units)

#### CODON-OPT-001: Sequence Optimization

| Field | Value |
|------|----------|
| **Canonical** | `CodonOptimizer.OptimizeSequence(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `OptimizeSequence(seq, organism, strategy)` | CodonOptimizer | Canonical |

**Strategies:** MaximizeCAI, BalancedOptimization, HarmonizeExpression, MinimizeSecondary, AvoidRareCodons

---

#### CODON-CAI-001: CAI Calculation

| Field | Value |
|------|----------|
| **Canonical** | `CodonOptimizer.CalculateCAI(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ CAI ≤ 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateCAI(codingSeq, codonTable, excludeSingleCodonAminoAcids=false)` | CodonOptimizer | Canonical |

---

#### CODON-RARE-001: Rare Codon Detection

| Field | Value |
|------|----------|
| **Canonical** | `CodonOptimizer.FindRareCodons(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindRareCodons(seq, codonTable, threshold)` | CodonOptimizer | Canonical |
| `CalculateMinMaxProfile(seq, codonTable, windowSize)` | CodonOptimizer | Canonical (%MinMax, Clarke & Clark 2008) |
| `FindRareCodonClusters(seq, codonTable, rareThreshold, windowSize, minRareCodons)` | CodonOptimizer | Canonical (Sherlocc RCC, Chartier et al. 2012) |

---

#### CODON-USAGE-001: Codon Usage Analysis

| Field | Value |
|------|----------|
| **Canonical** | `CodonOptimizer.CalculateCodonUsage(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateCodonUsage(seq)` | CodonOptimizer | Canonical |
| `CompareCodonUsage(seq1, seq2)` | CodonOptimizer | Comparison |

**Organism Tables:** E. coli K12, S. cerevisiae, H. sapiens

---

### 13. Translation (2 units)

#### TRANS-CODON-001: Codon Translation

| Field | Value |
|------|----------|
| **Canonical** | `GeneticCode.Translate(...)` |
| **Complexity** | O(1) per codon |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Translate(codon)` | GeneticCode | Canonical |
| `IsStartCodon(codon)` | GeneticCode | Check |
| `IsStopCodon(codon)` | GeneticCode | Check |
| `GetCodonsForAminoAcid(aa)` | GeneticCode | Reverse lookup |

**Genetic Codes:** Standard (1), Vertebrate Mito (2), Yeast Mito (3), Bacterial (11)

---

#### TRANS-PROT-001: Protein Translation

| Field | Value |
|------|----------|
| **Canonical** | `Translator.Translate(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Translate(dna, geneticCode)` | Translator | Canonical |

---

### 14. File I/O Parsers (7 units)

#### PARSE-FASTA-001: FASTA Parsing

| Field | Value |
|------|----------|
| **Canonical** | `FastaParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | FastaParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(fastaContent)` | FastaParser | Canonical |
| `ParseFile(filePath)` | FastaParser | File |
| `ParseFileAsync(filePath)` | FastaParser | Async |
| `ToFasta(entries, lineWidth)` | FastaParser | Export |
| `WriteFile(filePath, entries)` | FastaParser | Write |

**Edge Cases:**
- [ ] Empty content
- [ ] Multi-line sequences
- [ ] Special characters in headers
- [ ] Missing sequence after header

---

#### PARSE-FASTQ-001: FASTQ Parsing

| Field | Value |
|------|----------|
| **Canonical** | `FastqParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | FastqParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content, encoding)` | FastqParser | Canonical |
| `ParseFile(filePath, encoding)` | FastqParser | File |
| `CalculateStatistics(records)` | FastqParser | Stats |
| `FilterByQuality(records, minQ)` | FastqParser | Filter |

**Edge Cases:**
- [ ] Phred+33 vs Phred+64 encoding
- [ ] Malformed quality strings
- [ ] Empty records
- [ ] Auto-detect encoding

---

#### PARSE-BED-001: BED File Parsing

| Field | Value |
|------|----------|
| **Canonical** | `BedParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | BedParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content, format)` | BedParser | Canonical |
| `ParseFile(filePath, format)` | BedParser | File |
| `FilterByChrom(records, chrom)` | BedParser | Filter |
| `FilterByRegion(records, ...)` | BedParser | Filter |
| `MergeOverlapping(records)` | BedParser | Merge |
| `Intersect(records1, records2)` | BedParser | Set op |

**Edge Cases:**
- [ ] BED3 vs BED6 vs BED12 formats
- [ ] Block structures (BED12)
- [ ] Invalid coordinates

---

#### PARSE-VCF-001: VCF Parsing

| Field | Value |
|------|----------|
| **Canonical** | `VcfParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | VcfParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content)` | VcfParser | Canonical |
| `ParseFile(filePath)` | VcfParser | File |
| `ParseWithHeader(content)` | VcfParser | With header |
| `GetVariantType(record)` | VcfParser | Classification |

**Edge Cases:**
- [ ] Multi-allelic variants
- [ ] Missing values (.)
- [ ] Complex INFO fields
- [ ] Sample genotypes

---

#### PARSE-GFF-001: GFF/GTF Parsing

| Field | Value |
|------|----------|
| **Canonical** | `GffParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | GffParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content)` | GffParser | Canonical |
| `ParseFile(filePath)` | GffParser | File |
| `ToGff3(features)` | GffParser | Export |

---

#### PARSE-GENBANK-001: GenBank Parsing

| Field | Value |
|------|----------|
| **Canonical** | `GenBankParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | GenBankParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content)` | GenBankParser | Canonical |
| `ParseFile(filePath)` | GenBankParser | File |
| `ExtractFeatures(record)` | GenBankParser | Features |

---

#### PARSE-EMBL-001: EMBL Parsing

| Field | Value |
|------|----------|
| **Canonical** | `EmblParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | EmblParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content)` | EmblParser | Canonical |
| `ParseFile(filePath)` | EmblParser | File |

---

### 15. Sequence Complexity (3 units)

#### SEQ-COMPLEX-001: Linguistic Complexity

| Field | Value |
|------|----------|
| **Canonical** | `SequenceComplexity.CalculateLinguisticComplexity(...)` |
| **Complexity** | O(n × k) |
| **Invariant** | 0 ≤ result ≤ 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateLinguisticComplexity(seq, maxLen)` | SequenceComplexity | Canonical |
| `CalculateLinguisticComplexity(string, maxLen)` | SequenceComplexity | String |
| `FindLowComplexityRegions(seq, window, threshold)` | SequenceComplexity | Regions |
| `MaskLowComplexity(seq, threshold)` | SequenceComplexity | Masking |

---

#### SEQ-ENTROPY-001: Shannon Entropy

| Field | Value |
|------|----------|
| **Canonical** | `SequenceComplexity.CalculateShannonEntropy(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ result ≤ 2 for DNA |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateShannonEntropy(sequence)` | SequenceComplexity | Canonical |
| `CalculateKmerEntropy(seq, k)` | SequenceComplexity | K-mer based |

---

#### SEQ-GCSKEW-001: GC Skew

| Field | Value |
|------|----------|
| **Canonical** | `GcSkewCalculator.CalculateGcSkew(...)` |
| **Complexity** | O(n) |
| **Invariant** | -1 ≤ result ≤ 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateGcSkew(sequence)` | GcSkewCalculator | Canonical |
| `CalculateWindowedGcSkew(seq, window, step)` | GcSkewCalculator | Windowed |
| `CalculateCumulativeGcSkew(sequence)` | GcSkewCalculator | Cumulative |
| `PredictReplicationOrigin(sequence)` | GcSkewCalculator | Origin detection |

---

### 16. RNA Secondary Structure (3 units)

#### RNA-STRUCT-001: Secondary Structure Prediction

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.Predict(...)` |
| **Complexity** | O(n³) |
| **Class** | RnaSecondaryStructure |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Predict(sequence)` | RnaSecondaryStructure | Canonical |
| `PredictWithConstraints(seq, constraints)` | RnaSecondaryStructure | Constrained |
| `ToDotBracket(structure)` | RnaSecondaryStructure | Notation |
| `FromDotBracket(notation)` | RnaSecondaryStructure | Parse |

---

#### RNA-STEMLOOP-001: Stem-Loop Detection

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.FindStemLoops(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindStemLoops(sequence, minStem, maxLoop)` | RnaSecondaryStructure | Canonical |
| `FindHairpins(sequence, params)` | RnaSecondaryStructure | Hairpins |
| `FindPseudoknots(sequence)` | RnaSecondaryStructure | Pseudoknots |

---

#### RNA-ENERGY-001: Free Energy Calculation

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.CalculateFreeEnergy(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateFreeEnergy(structure)` | RnaSecondaryStructure | Canonical |
| `CalculateStackingEnergy(bp1, bp2)` | RnaSecondaryStructure | Stacking |

---

### 17. MicroRNA Analysis (3 units)

#### MIRNA-SEED-001: Seed Sequence Analysis

| Field | Value |
|------|----------|
| **Canonical** | `MiRnaAnalyzer.GetSeedSequence(...)` |
| **Complexity** | O(1) |
| **Class** | MiRnaAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetSeedSequence(miRnaSequence)` | MiRnaAnalyzer | Canonical |
| `CreateMiRna(name, sequence)` | MiRnaAnalyzer | Factory |
| `CompareSeedRegions(mirna1, mirna2)` | MiRnaAnalyzer | Compare |

---

#### MIRNA-TARGET-001: Target Site Prediction

| Field | Value |
|------|----------|
| **Canonical** | `MiRnaAnalyzer.FindTargetSites(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindTargetSites(mRna, miRna, minScore)` | MiRnaAnalyzer | Canonical |
| `ScoreTargetSite(site)` | MiRnaAnalyzer | Scoring |

**Site Types:**
- [x] 8mer, 7mer-m8, 7mer-A1, 6mer
- [x] Offset 6mer (positions 3-8)
- [ ] Supplementary pairing
- [x] ~~Centered sites~~ (removed per TargetScan 8.0 / McGeary 2019)

---

#### MIRNA-PRECURSOR-001: Pre-miRNA Detection

| Field | Value |
|------|----------|
| **Canonical** | `MiRnaAnalyzer.FindPreMiRnaHairpins(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindPreMiRnaHairpins(sequence)` | MiRnaAnalyzer | Canonical (default consecutive-pairing heuristic) |
| `FindPreMiRnaHairpinsByMfe(sequence, ...)` | MiRnaAnalyzer | Canonical (opt-in: real MFE structure via RNA-STRUCT-001 folder) |
| `AssessHairpinByMfe(candidate, minMfei, minLoopSize)` | MiRnaAnalyzer | Canonical (opt-in single-candidate MFE assessment) |
| `CalculateMfeIndex(freeEnergy, length, gcPercent)` | MiRnaAnalyzer | MFEI = AMFE/(G+C)% (Zhang 2006) |

---

### 18. Splice Site Prediction (3 units)

#### SPLICE-DONOR-001: Donor Site Detection

| Field | Value |
|------|----------|
| **Canonical** | `SpliceSitePredictor.FindDonorSites(...)` |
| **Complexity** | O(n) |
| **Class** | SpliceSitePredictor |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindDonorSites(sequence, minScore)` | SpliceSitePredictor | Canonical |
| `ScoreDonorSite(context)` | SpliceSitePredictor | Scoring |

---

#### SPLICE-ACCEPTOR-001: Acceptor Site Detection

| Field | Value |
|------|----------|
| **Canonical** | `SpliceSitePredictor.FindAcceptorSites(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindAcceptorSites(sequence, minScore)` | SpliceSitePredictor | Canonical |
| `ScoreAcceptorSite(context)` | SpliceSitePredictor | Scoring |

---

#### SPLICE-PREDICT-001: Gene Structure Prediction

| Field | Value |
|------|----------|
| **Canonical** | `SpliceSitePredictor.PredictGeneStructure(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictGeneStructure(sequence)` | SpliceSitePredictor | Canonical |
| `FindIntrons(sequence)` | SpliceSitePredictor | Introns |
| `FindExons(sequence)` | SpliceSitePredictor | Exons |

---

### 19. Protein Disorder Prediction (2 units)

#### DISORDER-PRED-001: Disorder Prediction

| Field | Value |
|------|----------|
| **Canonical** | `DisorderPredictor.PredictDisorder(...)` |
| **Complexity** | O(n) |
| **Class** | DisorderPredictor |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictDisorder(sequence, windowSize, threshold)` | DisorderPredictor | Canonical |
| `CalculateDisorderScore(window)` | DisorderPredictor | Score |
| `CalculateHydropathy(sequence)` | DisorderPredictor | Hydropathy |

---

#### DISORDER-REGION-001: Disordered Region Detection ☑

| Field | Value |
|------|----------|
| **Status** | ☑ Complete |
| **Canonical** | `DisorderPredictor.IdentifyDisorderedRegions(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `IdentifyDisorderedRegions(predictions, minLen)` | DisorderPredictor | Canonical |
| `ClassifyRegionType(region)` | DisorderPredictor | Classification |

---

### 20. Protein Motif Finding (3 units)

#### PROTMOTIF-FIND-001: Motif Search

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.FindMotifs(...)` |
| **Complexity** | O(n × m) |
| **Class** | ProteinMotifFinder |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindMotifs(sequence, patterns)` | ProteinMotifFinder | Canonical |
| `FindAllKnownMotifs(sequence)` | ProteinMotifFinder | All patterns |
| `ScanForPattern(sequence, pattern)` | ProteinMotifFinder | Single |

---

#### PROTMOTIF-PROSITE-001: PROSITE Pattern Matching ☑

| Field | Value |
|------|----------|
| **Status** | ☑ Complete |
| **Canonical** | `ProteinMotifFinder.FindMotifByProsite(...)` / `ConvertPrositeToRegex(...)` |
| **Complexity** | O(n × m) |
| **Evidence** | [Evidence](docs/Evidence/PROTMOTIF-PROSITE-001-Evidence.md) |
| **TestSpec** | [TestSpec](tests/TestSpecs/PROTMOTIF-PROSITE-001.md) |
| **Tests** | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_PrositePattern_Tests.cs) |
| **Algorithm Doc** | [Algorithm](docs/algorithms/ProteinMotif/PROSITE_Pattern_Matching.md) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindMotifByProsite(sequence, pattern, name)` | ProteinMotifFinder | Canonical |
| `ConvertPrositeToRegex(pattern)` | ProteinMotifFinder | Parse |

**Common Patterns:**
- [x] N-glycosylation (PS00001)
- [x] Phosphorylation sites (PS00004, PS00005)
- [x] Zinc fingers (PS00028)
- [x] EF-hand (PS00018)
- [x] Tachykinin / Pyrokinin with `[G>]` C-term brackets (PS00267, PS00539)

---

#### PROTMOTIF-DOMAIN-001: Domain Prediction

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.PredictDomains(...)` |
| **Complexity** | O(n × d) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictDomains(sequence)` | ProteinMotifFinder | Canonical |
| `PredictSignalPeptide(sequence)` | ProteinMotifFinder | Signal |

---

### 21. Epigenetics Analysis (3 units)

#### EPIGEN-CPG-001: CpG Site Detection

| Field | Value |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.FindCpGSites(...)` |
| **Complexity** | O(n) |
| **Class** | EpigeneticsAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindCpGSites(sequence)` | EpigeneticsAnalyzer | Canonical |
| `FindCpGIslands(sequence, params)` | EpigeneticsAnalyzer | Islands |
| `CalculateCpGObservedExpected(sequence)` | EpigeneticsAnalyzer | O/E ratio |

---

#### EPIGEN-METHYL-001: Methylation Analysis

| Field | Value |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.FindMethylationSites(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindMethylationSites(sequence)` | EpigeneticsAnalyzer | Canonical |
| `CalculateMethylationProfile(sites)` | EpigeneticsAnalyzer | Profile |
| `GetMethylationContext(site)` | EpigeneticsAnalyzer | Context (CpG/CHG/CHH) |

---

#### EPIGEN-DMR-001: Differentially Methylated Regions

| Field | Value |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.FindDMRs(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindDMRs(profile1, profile2, threshold)` | EpigeneticsAnalyzer | Canonical |
| `AnnotateDMRs(dmrs, annotations)` | EpigeneticsAnalyzer | Annotate |

---

### 22. Variant Calling (4 units)

#### VARIANT-CALL-001: Variant Detection

| Field | Value |
|------|----------|
| **Canonical** | `VariantCaller.CallVariants(...)` |
| **Complexity** | O(n × m) |
| **Class** | VariantCaller |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CallVariants(reference, query)` | VariantCaller | Canonical |
| `CallVariantsFromAlignment(aligned1, aligned2)` | VariantCaller | From alignment |
| `ClassifyVariant(variant)` | VariantCaller | Classification |

---

#### VARIANT-SNP-001: SNP Detection

| Field | Value |
|------|----------|
| **Canonical** | `VariantCaller.FindSnps(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSnps(reference, query)` | VariantCaller | Canonical |
| `FindSnpsDirect(ref, query)` | VariantCaller | Direct (no alignment) |

---

#### VARIANT-INDEL-001: Indel Detection

| Field | Value |
|------|----------|
| **Canonical** | `VariantCaller.FindInsertions/FindDeletions(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindInsertions(reference, query)` | VariantCaller | Insertions |
| `FindDeletions(reference, query)` | VariantCaller | Deletions |

---

#### VARIANT-ANNOT-001: Variant Annotation ☑

| Field | Value |
|------|----------|
| **Canonical** | `VariantAnnotator.Annotate(...)` |
| **Complexity** | O(v × g) |
| **Class** | VariantAnnotator |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Annotate(variants, annotations, referenceSequence?)` | VariantAnnotator | Canonical |
| `PredictFunctionalImpact(variant, transcript, referenceSequence)` | VariantAnnotator | Impact |

> Implemented per Ensembl VEP (McLaren et al. 2016) SO consequence terms + IMPACT/rank
> (Constants.pm) and VariationEffect.pm coding predicates; codons translated with the
> NCBI Standard code (table 1). `PredictFunctionalImpact` takes the overlapping
> transcript and a forward-strand reference window to translate ref/alt codons.

---

### 23. Structural Variant Analysis (3 units)

#### SV-DETECT-001: SV Detection ☑

| Field | Value |
|------|----------|
| **Canonical** | `StructuralVariantAnalyzer.DetectSVs(...)` |
| **Complexity** | O(n log n) |
| **Class** | StructuralVariantAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectSVs(readPairs, params)` | StructuralVariantAnalyzer | Canonical |
| `FindDiscordantPairs(readPairs, params)` | StructuralVariantAnalyzer | Discordant |
| `ClassifySV(pair, params)` | StructuralVariantAnalyzer | Classification |

**SV Types:**
- [x] Deletion, Duplication, Inversion
- [x] Insertion, Translocation
- [x] Complex rearrangements

---

#### SV-BREAKPOINT-001: Breakpoint Detection

| Field | Value |
|------|----------|
| **Canonical** | `StructuralVariantAnalyzer.FindBreakpoints(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindBreakpoints(splitReads)` | StructuralVariantAnalyzer | Canonical |
| `RefineBreakpoint(region, reads)` | StructuralVariantAnalyzer | Refinement |

---

#### SV-CNV-001: Copy Number Variation

| Field | Value |
|------|----------|
| **Canonical** | `StructuralVariantAnalyzer.DetectCNV(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectCNV(depthData, windowSize, referenceDepth, chromosome)` | StructuralVariantAnalyzer | Canonical |
| `SegmentCopyNumber(logRatios, chromosome)` | StructuralVariantAnalyzer | Segmentation |

**Status:** ☑ Complete — read-depth windowed log2-ratio → integer copy number (Yoon et al. 2009; CNVkit `CN = round(2·2^log2)`). Evidence: docs/Evidence/SV-CNV-001-Evidence.md; TestSpec: tests/TestSpecs/SV-CNV-001.md; Tests: StructuralVariantAnalyzer_DetectCNV_Tests.cs.

---

### 24. Sequence Assembly (3 units)

#### ASSEMBLY-OLC-001: Overlap-Layout-Consensus ☑

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAssembler.AssembleOLC(...)` |
| **Complexity** | O(n² × m) |
| **Class** | SequenceAssembler |
| **Evidence** | [ASSEMBLY-OLC-001-Evidence.md](docs/Evidence/ASSEMBLY-OLC-001-Evidence.md) |
| **TestSpec** | [ASSEMBLY-OLC-001.md](tests/TestSpecs/ASSEMBLY-OLC-001.md) |
| **Algorithm doc** | [Overlap_Layout_Consensus.md](docs/algorithms/Assembly/Overlap_Layout_Consensus.md) |
| **Tests** | [SequenceAssembler_AssembleOLC_Tests.cs](tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_AssembleOLC_Tests.cs) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AssembleOLC(reads, params)` | SequenceAssembler | Canonical |
| `FindAllOverlaps(reads, minOverlap, minId)` | SequenceAssembler | Overlaps |

**Edge Cases (derived from Evidence — no prior by-area list):**
- [x] null / empty read set → empty `AssemblyResult` (no exception)
- [x] no above-threshold overlap → each read is its own singleton contig (INV-05)
- [x] unambiguous chain → single merged superstring contig (INV-04)
- [x] overlap exactly = `minOverlap` accepted; below rejected (INV-02)
- [x] identity threshold gates approximate overlaps
- [x] case-insensitive overlaps
- [x] repeat longer than read length not collapsed below longest read (ASM-02)

**Definition of Done:**
- [x] TestSpec, tests, Evidence, algorithm doc created; Registry updated
- [x] All tests pass (`dotnet test`); deterministic; zero warnings in changed files
- [x] Edge cases covered (null, empty, boundary, invalid) for in-scope methods
- [x] O(n²) algorithm: property-based invariant test (INV-04) + performance baseline recorded (algorithm doc §7.2)

---

#### ASSEMBLY-DBG-001: De Bruijn Graph Assembly

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAssembler.AssembleDeBruijn(...)` |
| **Complexity** | O(n × k) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AssembleDeBruijn(reads, params)` | SequenceAssembler | Canonical |
| `BuildDeBruijnGraph(reads, k)` | SequenceAssembler | Graph construction |

---

#### ASSEMBLY-STATS-001: Assembly Statistics ☑

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAssemblyAnalyzer.CalculateStatistics(...)` |
| **Complexity** | O(n log n) |
| **Class** | GenomeAssemblyAnalyzer |

**Status:** ☑ Complete — N50/L50 = "smallest contig in the smallest set of largest contigs whose combined length is ≥ 50% of the assembly" (Miller, Koren & Sutton 2010, Genomics 95(6):315-327, §1.2); inclusive ≥ threshold% cumulative test verified against QUAST `quast_libs/N50.py`; auN = Σl²/Σl (Li 2020). Worked example {80,70,50,40,30,20}→N50=70,L50=2 (Wikipedia). Evidence: docs/Evidence/ASSEMBLY-STATS-001-Evidence.md; TestSpec: tests/TestSpecs/ASSEMBLY-STATS-001.md; Tests: GenomeAssemblyAnalyzer_AssemblyStatistics_Tests.cs.

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateStatistics(contigs)` | GenomeAssemblyAnalyzer | Canonical |
| `CalculateN50(contigs)` | GenomeAssemblyAnalyzer | N50 |
| `CalculateNx(contigs, threshold)` | GenomeAssemblyAnalyzer | Nx/Lx |
| `CalculateAuN(contigs)` | GenomeAssemblyAnalyzer | auN (Σl²/Σl) |
| `FindGaps(sequence)` | GenomeAssemblyAnalyzer | Gap detection |

**Edge Cases:**
- [x] Empty contig set → all-zero statistics / Nx=Lx=0 / auN=0
- [x] Single contig → N50 = its length, L50 = 1
- [x] Cumulative reaching exactly threshold% (inclusive boundary)
- [x] Leading / trailing / all-N gap runs
- [x] minGapLength filtering

---

### 25. Transcriptome Analysis (3 units)

#### TRANS-EXPR-001: Expression Quantification

| Field | Value |
|------|----------|
| **Canonical** | `TranscriptomeAnalyzer.CalculateTPM(...)` |
| **Complexity** | O(n) |
| **Class** | TranscriptomeAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateTPM(geneCounts)` | TranscriptomeAnalyzer | TPM |
| `CalculateFPKM(count, length, total)` | TranscriptomeAnalyzer | FPKM |
| `QuantileNormalize(samples)` | TranscriptomeAnalyzer | Normalization |

---

#### TRANS-DIFF-001: Differential Expression

| Field | Value |
|------|----------|
| **Canonical** | `TranscriptomeAnalyzer.FindDifferentiallyExpressed(...)` |
| **Complexity** | O(g × s) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindDifferentiallyExpressed(cond1, cond2, alpha)` | TranscriptomeAnalyzer | Canonical |
| `CalculateFoldChange(expr1, expr2)` | TranscriptomeAnalyzer | Fold change |

---

#### TRANS-SPLICE-001: Alternative Splicing ☑

| Field | Value |
|------|----------|
| **Canonical** | `TranscriptomeAnalyzer.DetectAlternativeSplicing(...)` |
| **Complexity** | O(n) (PSI O(1); classification O(g·k²·e)) |
| **Evidence** | [TRANS-SPLICE-001-Evidence.md](docs/Evidence/TRANS-SPLICE-001-Evidence.md) |
| **TestSpec** | [TRANS-SPLICE-001.md](tests/TestSpecs/TRANS-SPLICE-001.md) |
| **Tests** | [TranscriptomeAnalyzer_AlternativeSplicing_Tests.cs](tests/Seqeron/Seqeron.Genomics.Tests/TranscriptomeAnalyzer_AlternativeSplicing_Tests.cs) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectAlternativeSplicing(isoforms)` | TranscriptomeAnalyzer | Canonical |
| `CalculatePSI(inclusionReads, exclusionReads, ...)` | TranscriptomeAnalyzer | Canonical (Percent Spliced In) |

**Edge Cases:**
- [x] PSI of 0/0 (no reads) → NaN (undefined)
- [x] PSI = 1 (S=0) and PSI = 0 (I=0)
- [x] rMATS length-normalized PSI when both effective lengths > 0
- [x] negative read counts → ArgumentOutOfRangeException
- [x] fewer than two isoforms / identical isoforms → no event
- [x] null / empty isoform input → empty result
- [x] all five event classes (SE, RI, A5SS, A3SS, MXE) classified

---

### 26. Comparative Genomics (3 units)

#### COMPGEN-SYNTENY-001: Synteny Detection ☑

| Field | Value |
|------|----------|
| **Status** | ☑ Complete (MCScanX collinearity model; Wang et al. 2012, NAR 40(7):e49) |
| **Canonical** | `ComparativeGenomics.FindSyntenicBlocks(...)` |
| **Complexity** | O(n²) |
| **Class** | ComparativeGenomics |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSyntenicBlocks(genes1, genes2, orthologs)` | ComparativeGenomics | Canonical |
| `VisualizeSynteny(blocks)` | ComparativeGenomics | Visualization |

---

#### ☑ COMPGEN-ORTHO-001: Ortholog Identification

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.FindOrthologs(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindOrthologs(genes1, genes2, minIdentity, minCoverage)` | ComparativeGenomics | Canonical (Reciprocal Best Hits) |
| `FindParalogs(genes, minIdentity, minCoverage)` | ComparativeGenomics | Canonical (within-genome best hits / in-paralogs) |

**Note:** No prior by-area Edge Cases list existed; canonical methods, signatures, invariants, and edge cases were derived from the Evidence (Tatusov 1997; Moreno-Hagelsieb & Latimer 2008; Fitch 1970; Remm et al. 2001). `FindOrthologs` was corrected from a one-directional best hit to the reciprocal (RBH) criterion; `FindParalogs` was newly implemented.

---

#### COMPGEN-REARR-001: Genome Rearrangements

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.DetectRearrangements(...)` |
| **Complexity** | O(n log n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectRearrangements(blocks)` | ComparativeGenomics | Canonical |
| `ClassifyRearrangement(event)` | ComparativeGenomics | Classification |

---

### 27. Pan-Genome Analysis (2 units)

#### PANGEN-CORE-001: Core/Accessory Genome

| Field | Value |
|------|----------|
| **Canonical** | `PanGenomeAnalyzer.ConstructPanGenome(...)` |
| **Complexity** | O(g² × s) |
| **Class** | PanGenomeAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ConstructPanGenome(genomes, identityThreshold, coreFraction)` | PanGenomeAnalyzer | Canonical |
| `GetCoreGeneClusters(clusters, totalGenomes, threshold)` | PanGenomeAnalyzer | Core genes (implements `IdentifyCoreGenes` referent) |

---

#### PANGEN-CLUSTER-001: Gene Clustering

| Field | Value |
|------|----------|
| **Canonical** | `PanGenomeAnalyzer.ClusterGenes(...)` |
| **Complexity** | O(g² × s) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ClusterGenes(genomes, idThreshold)` | PanGenomeAnalyzer | Canonical |
| `CreatePresenceAbsenceMatrix(genomes, clusters)` | PanGenomeAnalyzer | Matrix |

> Note (PANGEN-CLUSTER-001): the matrix method ships as `CreatePresenceAbsenceMatrix` (sibling-consistent name), not `GeneratePresenceAbsenceMatrix`; API-naming only, see TestSpec §7. `ClusterGenes` was corrected from a k-mer Jaccard heuristic to CD-HIT global sequence identity (identical positions / shorter length) per Li & Godzik (2006).

---

### 28. Quality Score Analysis (2 units)

#### QUALITY-PHRED-001: Phred Score Handling

| Field | Value |
|------|----------|
| **Canonical** | `QualityScoreAnalyzer.ParseQualityString(...)` |
| **Complexity** | O(n) |
| **Class** | QualityScoreAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ParseQualityString(qualStr, encoding)` | QualityScoreAnalyzer | Canonical |
| `ToQualityString(scores, encoding)` | QualityScoreAnalyzer | Export |
| `ConvertEncoding(qualStr, from, to)` | QualityScoreAnalyzer | Convert |

---

#### QUALITY-STATS-001: Quality Statistics

| Field | Value |
|------|----------|
| **Canonical** | `QualityScoreAnalyzer.CalculateStatistics(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateStatistics(scores)` | QualityScoreAnalyzer | Canonical |
| `CalculateQ30Percentage(scores)` | QualityScoreAnalyzer | Q30 |

---

### 29. Probe Design (2 units)

#### PROBE-DESIGN-001: Hybridization Probe Design

| Field | Value |
|------|----------|
| **Canonical** | `ProbeDesigner.DesignProbes(...)` |
| **Complexity** | O(n²) |
| **Class** | ProbeDesigner |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DesignProbes(template, params)` | ProbeDesigner | Canonical |
| `DesignTilingProbes(template, overlap)` | ProbeDesigner | Tiling |
| `ScoreProbe(sequence, params)` | ProbeDesigner | Scoring |

---

#### PROBE-VALID-001: Probe Validation

| Field | Value |
|------|----------|
| **Canonical** | `ProbeDesigner.ValidateProbe(...)` |
| **Complexity** | O(n × g) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ValidateProbe(probe, genome)` | ProbeDesigner | Canonical |
| `CheckSpecificity(probe, database)` | ProbeDesigner | Specificity |

---

### 30. Sequence Statistics (7 units) (7 units)

#### SEQ-STATS-001: Sequence Composition Statistics

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateNucleotideComposition(...)` |
| **Complexity** | O(n) |
| **Class** | SequenceStatistics |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateNucleotideComposition(sequence)` | SequenceStatistics | Canonical |
| `CalculateAminoAcidComposition(sequence)` | SequenceStatistics | Protein |
| `SummarizeNucleotideSequence(sequence)` | SequenceStatistics | Summary |

---

#### SEQ-MW-001: Molecular Weight Calculation ☑

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateMolecularWeight(...)` |
| **Complexity** | O(n) |
| **Status** | ☑ Complete — Evidence/TestSpec/Tests linked in Processing Registry |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateMolecularWeight(proteinSequence)` | SequenceStatistics | Protein |
| `CalculateNucleotideMolecularWeight(sequence, isDna)` | SequenceStatistics | DNA/RNA |

---

#### SEQ-PI-001: Isoelectric Point Calculation

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateIsoelectricPoint(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ pI ≤ 14 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateIsoelectricPoint(proteinSequence)` | SequenceStatistics | Canonical |

---

#### SEQ-HYDRO-001: Hydrophobicity Analysis

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateHydrophobicity(...)` |
| **Complexity** | O(n) |
| **Scale** | Kyte-Doolittle |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateHydrophobicity(proteinSequence)` | SequenceStatistics | GRAVY |
| `CalculateHydrophobicityProfile(sequence, windowSize)` | SequenceStatistics | Profile |

---

#### SEQ-THERMO-001: Thermodynamic Properties

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateThermodynamics(...)` |
| **Complexity** | O(n) |
| **Method** | Nearest-neighbor |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateThermodynamics(dnaSequence, naConc, primerConc)` | SequenceStatistics | Canonical |
| `CalculateMeltingTemperature(dnaSequence, useWallaceRule)` | SequenceStatistics | Simple Tm |

---

#### SEQ-DINUC-001: Dinucleotide Analysis ☑

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateDinucleotideFrequencies(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateDinucleotideFrequencies(sequence)` | SequenceStatistics | Frequencies |
| `CalculateDinucleotideRatios(sequence)` | SequenceStatistics | O/E ratios |
| `CalculateCodonFrequencies(dnaSequence, readingFrame)` | SequenceStatistics | Codons |

---

#### SEQ-SECSTRUCT-001: Secondary Structure Prediction

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.PredictSecondaryStructure(...)` |
| **Complexity** | O(n) |
| **Method** | Chou-Fasman propensities |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictSecondaryStructure(proteinSequence, windowSize)` | SequenceStatistics | Canonical |

---

### 31. Codon Usage Analysis (3 units)

#### CODON-RSCU-001: Relative Synonymous Codon Usage

| Field | Value |
|------|----------|
| **Canonical** | `CodonUsageAnalyzer.CalculateRscu(...)` |
| **Complexity** | O(n) |
| **Class** | CodonUsageAnalyzer |
| **Invariant** | RSCU = 1 means no bias |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateRscu(sequence)` | CodonUsageAnalyzer | Canonical |
| `CountCodons(sequence)` | CodonUsageAnalyzer | Counting |

---

#### CODON-ENC-001: Effective Number of Codons

| Field | Value |
|------|----------|
| **Canonical** | `CodonUsageAnalyzer.CalculateEnc(...)` |
| **Complexity** | O(n) |
| **Invariant** | 20 ≤ ENC ≤ 61 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateEnc(sequence)` | CodonUsageAnalyzer | Canonical |

---

#### CODON-STATS-001: Codon Usage Statistics

| Field | Value |
|------|----------|
| **Canonical** | `CodonUsageAnalyzer.GetStatistics(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetStatistics(sequence)` | CodonUsageAnalyzer | Canonical |
| `CalculateCai(sequence, referenceRscu)` | CodonUsageAnalyzer | CAI |
| `EColiOptimalCodons` (property) | CodonUsageAnalyzer | Reference |
| `HumanOptimalCodons` (property) | CodonUsageAnalyzer | Reference |

---

### 32. Extended Translation (1 unit)

#### TRANS-SIXFRAME-001: Six-Frame Translation

| Field | Value |
|------|----------|
| **Canonical** | `Translator.TranslateSixFrames(...)` |
| **Complexity** | O(n) |
| **Class** | Translator |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `TranslateSixFrames(dna, geneticCode)` | Translator | Canonical |
| `FindOrfs(dna, geneticCode, minLength, searchBothStrands)` | Translator | ORF finding |

---

### 33. Extended Assembly (6 units)

#### ASSEMBLY-MERGE-001: Contig Merging

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAssembler.MergeContigs(...)` |
| **Complexity** | O(n) |
| **Class** | SequenceAssembler |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `MergeContigs(contig1, contig2, overlapLength)` | SequenceAssembler | Canonical |

---

#### ASSEMBLY-SCAFFOLD-001: Scaffolding

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAssembler.Scaffold(...)` |
| **Complexity** | O(n + k) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Scaffold(contigs, links, gapCharacter)` | SequenceAssembler | Canonical |

---

#### ASSEMBLY-COVER-001: Coverage Calculation

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAssembler.CalculateCoverage(...)` |
| **Complexity** | O(n × r) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateCoverage(reference, reads, minOverlap)` | SequenceAssembler | Canonical |

---

#### ASSEMBLY-CONSENSUS-001: Consensus Computation

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAssembler.ComputeConsensus(...)` |
| **Complexity** | O(n × r) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ComputeConsensus(alignedReads)` | SequenceAssembler | Canonical |

---

#### ASSEMBLY-TRIM-001: Quality Trimming

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAssembler.QualityTrimReads(...)` |
| **Complexity** | O(n × r) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `QualityTrimReads(reads, minQuality, minLength)` | SequenceAssembler | Canonical |

---

#### ASSEMBLY-CORRECT-001: Error Correction

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAssembler.ErrorCorrectReads(...)` |
| **Complexity** | O(n × r × k) |
| **Method** | K-mer frequency based |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ErrorCorrectReads(reads, kmerSize, minKmerFrequency)` | SequenceAssembler | Canonical |

---

### 34. Extended Approximate Matching (1 unit)

#### PAT-APPROX-003: Best Match and Frequency Analysis

| Field | Value |
|------|----------|
| **Canonical** | `ApproximateMatcher.FindBestMatch(...)` |
| **Complexity** | O(n × m) |
| **Class** | ApproximateMatcher |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindBestMatch(sequence, pattern)` | ApproximateMatcher | Best single match |
| `CountApproximateOccurrences(sequence, pattern, maxMismatches)` | ApproximateMatcher | Counting |
| `FindFrequentKmersWithMismatches(sequence, k, d)` | ApproximateMatcher | Frequency |

---

### 35. Alignment Statistics (1 unit)

#### ☑ ALIGN-STATS-001: Alignment Statistics and Formatting

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAligner.CalculateStatistics(...)` |
| **Complexity** | O(n) |
| **Class** | SequenceAligner |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateStatistics(alignment, scoring)` | SequenceAligner | Statistics |
| `FormatAlignment(alignment, lineWidth, scoring)` | SequenceAligner | Formatting |

**Behavioral note (evidence-corrected):** Identity/Similarity/Gaps follow the EMBOSS
needle definition — denominator is the alignment length **including gap columns**.
Similarity counts identical columns **plus** non-identical columns whose substitution
score is **positive** (EMBOSS/BLAST "positives"); it is NOT `(matches+mismatches)/length`.
For the DNA models exposed by this class (Mismatch < 0) Similarity equals Identity.
`FormatAlignment` uses the srspair markup legend (`|` identity, `:` similarity, space
gap/mismatch). The prior implementation's similarity formula was corrected in this unit.

- [x] Edge cases: null alignment, empty alignment, all-gap, perfect identity, non-positive lineWidth
- [x] Evidence-based exact values (EMBOSS HBA/HBB 43.6/60.4/6.0), invariants INV-1..INV-4
- [x] Tests pass, zero warnings in changed files

---

### 36. Extended Epigenetics (3 units)

#### EPIGEN-BISULF-001: Bisulfite Sequencing Analysis

| Field | Value |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.SimulateBisulfiteConversion(...)` |
| **Complexity** | O(n) |
| **Class** | EpigeneticsAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `SimulateBisulfiteConversion(sequence)` | EpigeneticsAnalyzer | Canonical |
| `CalculateMethylationFromBisulfite(bsSeq, refSeq)` | EpigeneticsAnalyzer | Calculation |
| `GenerateMethylationProfile(sites)` | EpigeneticsAnalyzer | Profile |

---

#### EPIGEN-CHROM-001: Chromatin State Prediction

| Field | Value |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.PredictChromatinState(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictChromatinState(histoneMarks)` | EpigeneticsAnalyzer | Canonical |
| `AnnotateHistoneModifications(peaks, annotations)` | EpigeneticsAnalyzer | Annotation |
| `FindAccessibleRegions(atacSignal, threshold)` | EpigeneticsAnalyzer | ATAC-seq |

---

#### EPIGEN-AGE-001: Epigenetic Age Estimation

| Field | Value |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.CalculateEpigeneticAge(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateEpigeneticAge(methylationAtClockCpGs, coefficients, intercept)` | EpigeneticsAnalyzer | Canonical |
| `HorvathAntiTransform(transformedAge)` | EpigeneticsAnalyzer | Canonical |
| `PredictImprintedGenes(...)` | EpigeneticsAnalyzer | Out of scope (no retrieved authoritative basis; not part of the age clock) |

> Signature corrected from the placeholder `(methylationProfile, clockType)` to the evidence-based
> `(methylationAtClockCpGs, coefficients, intercept)`: Horvath (2013) computes DNAm age as
> `anti.trafo(intercept + Σ coef·β)` with caller-supplied clock coefficients (the 353-CpG table is a
> large published table and is not bundled — fabricated coefficients are forbidden). See EPIGEN-AGE-001 Evidence/TestSpec.

---

### 37. Extended MiRNA Analysis (1 unit)

#### MIRNA-PAIR-001: MiRNA-Target Pairing Analysis

| Field | Value |
|------|----------|
| **Canonical** | `MiRnaAnalyzer.AlignMiRnaToTarget(...)` |
| **Complexity** | O(n × m) |
| **Class** | MiRnaAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AlignMiRnaToTarget(miRna, target)` | MiRnaAnalyzer | Canonical |
| `GetReverseComplement(sequence)` | MiRnaAnalyzer | Helper |
| `CanPair(base1, base2)` | MiRnaAnalyzer | Pairing check |
| `IsWobblePair(base1, base2)` | MiRnaAnalyzer | Wobble check |

**Status:** ☑ Complete — Watson-Crick (A-U, C-G) + G:U wobble (Crick 1966) antiparallel
reverse-complement duplex pairing. Invented free-energy constants (−2.0/−1.0/−0.5/+0.5)
replaced by NNDB Turner 2004 nearest-neighbor stacking sum; free-energy magnitude is
"Intentionally simplified" (only sign validated). See MIRNA-PAIR-001 Evidence/TestSpec.

---

### 38. Extended Pan-Genome Analysis (2 units)

#### PANGEN-HEAP-001: Pan-Genome Growth Model ☑

| Field | Value |
|------|----------|
| **Canonical** | `PanGenomeAnalyzer.FitHeapsLaw(...)` |
| **Complexity** | O(P·G·C) (P perms, G genomes, C clusters) |
| **Class** | PanGenomeAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FitHeapsLaw(matrix \| genomes)` | PanGenomeAnalyzer | Canonical |
| `CreatePresenceAbsenceMatrix(genomes, clusters)` | PanGenomeAnalyzer | Matrix |

**Behavior (external evidence overrides checklist):** Heaps' law is the *new-gene-discovery
decay* model `n(N) = K·N^(−alpha)` fit over permuted genome orderings (Tettelin et al. 2008;
micropan `heaps()`), NOT a cumulative pan-genome *size* growth `P=K·N^+gamma`. Open ⇔ alpha < 1,
closed ⇔ alpha > 1. The prior implementation (log-log fit of cumulative size with a positive
exponent, gene-id matching) was non-conforming and was rewritten to the micropan model.

**Edge cases:**
- [x] Fewer than 2 genomes / null / empty → degenerate fit (Intercept 0, predictor 0), no exception.
- [x] First genome contributes no "new" count (curve uses N = 2..G).
- [x] Binary presence: duplicate/over-counted presence collapses to 1.
- [x] alpha bounded to [0,2], Intercept to [0,10000] (micropan optim box).
- [x] Exact power-curve data → analytic (K, alpha) recovered.

---

#### PANGEN-MARKER-001: Phylogenetic Marker Selection

| Field | Value |
|------|----------|
| **Canonical** | `PanGenomeAnalyzer.SelectPhylogeneticMarkers(...)` |
| **Complexity** | O(n × g) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `SelectPhylogeneticMarkers(clusters, criteria)` | PanGenomeAnalyzer | Canonical |
| `GetCoreGeneClusters(clusters, threshold)` | PanGenomeAnalyzer | Core genes |
| `CreateCoreGenomeAlignment(coreGenes)` | PanGenomeAnalyzer | Alignment |
| `GetSingletonGenes(clusters)` | PanGenomeAnalyzer | Singletons |

---

### 39. Extended Population Genetics (3 units)

#### POP-SELECT-001: Selection Signature Detection ☑

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex)` |
| **Complexity** | O(n × h) |
| **Class** | PopulationGeneticsAnalyzer |
| **Status** | ☑ Complete — unstandardized iHS = ln(iHH_A/iHH_D) per Voight et al. (2006); EHH per Sabeti et al. (2002) / selscan. Edge cases (null, empty, monomorphic, invalid allele, out-of-range, bad bin/window) covered. |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateIHS(haplotypes, positions)` | PopulationGeneticsAnalyzer | iHS |
| `ScanForSelection(variants, windowSize)` | PopulationGeneticsAnalyzer | Genome-wide |

---

#### ☑ POP-ANCESTRY-001: Ancestry Estimation

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.EstimateAncestry(...)` |
| **Complexity** | O(I × iterations × J × K) (per-iteration O(JK), F fixed) |
| **Status** | ☑ Complete — supervised/projection ADMIXTURE EM (Alexander et al. 2009, Eq. 2/4/5) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `EstimateAncestry(individuals, referencePops, maxIterations)` | PopulationGeneticsAnalyzer | Canonical |

**Edge Cases:**
- [x] Empty individuals → empty result
- [x] Empty reference panels → empty result
- [x] Genotype length ≠ panel SNP count → individual skipped
- [x] Genotype value outside {0,1,2} → SNP treated as missing
- [x] All genotypes missing → uniform prior returned
- [x] Identical panels → uniform stays uniform

---

#### POP-ROH-001: Runs of Homozygosity

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.FindROH(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindROH(genotypes, minLength, minSnps)` | PopulationGeneticsAnalyzer | Canonical |
| `CalculateInbreedingFromROH(rohSegments, genomeLength)` | PopulationGeneticsAnalyzer | F_ROH |
| `CalculatePairwiseFst(pop1, pop2, variants)` | PopulationGeneticsAnalyzer | Fst matrix |

---

### 40. Extended Metagenomics (4 units)

#### META-FUNC-001: Functional Prediction

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.PredictFunctions(...)` |
| **Complexity** | O(n × g) |
| **Class** | MetagenomicsAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictFunctions(genes, database)` | MetagenomicsAnalyzer | Canonical |
| `FindPathwayEnrichment(functions, pathwayDb)` | MetagenomicsAnalyzer | Enrichment |

---

#### META-RESIST-001: Antibiotic Resistance Detection

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.FindAntibioticResistanceGenes(...)` |
| **Complexity** | O(n × d) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindAntibioticResistanceGenes(contigs, argDatabase)` | MetagenomicsAnalyzer | Canonical |

---

#### META-PATHWAY-001: Metabolic Pathway Analysis

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.FindPathwayEnrichment(...)` |
| **Complexity** | O(n × p) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindPathwayEnrichment(functions, pathwayDb)` | MetagenomicsAnalyzer | Canonical |

---

#### META-TAXA-001: Significant Taxa Detection

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.FindSignificantTaxa(...)` |
| **Complexity** | O(t × s log s) (t taxa, s samples; per-taxon Mann–Whitney sort) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSignificantTaxa(profiles, groups, pThreshold, useContinuityCorrection)` | MetagenomicsAnalyzer | Canonical |
| `MannWhitneyU(group1, group2, useContinuityCorrection)` | MetagenomicsAnalyzer | Canonical |

**Behavior:** per-taxon Wilcoxon rank-sum (Mann–Whitney U) test under the asymptotic normal approximation with midrank tie correction and a SciPy-default continuity correction; a taxon is significant when its two-tailed p-value < `pThreshold`. Source: Mann & Whitney (1947), SciPy `mannwhitneyu`, Xia & Sun (2017). No by-area edge-case list pre-existed; edge cases derived from Evidence (null/empty/single-group/invalid-label/all-tied/absent-taxon).

---

### 41. Extended Phylogenetics (2 units)

#### PHYLO-BOOT-001: Bootstrap Analysis

| Field | Value |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.Bootstrap(...)` |
| **Complexity** | O(b × n³) |
| **Class** | PhylogeneticAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Bootstrap(sequences, replicates, distanceMethod, treeMethod, seed)` | PhylogeneticAnalyzer | Canonical |

---

#### PHYLO-STATS-001: Tree Statistics

| Field | Value |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.GetLeaves(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetLeaves(tree)` | PhylogeneticAnalyzer | Leaves |
| `CalculateTreeLength(tree)` | PhylogeneticAnalyzer | Total length |
| `GetTreeDepth(tree)` | PhylogeneticAnalyzer | Max depth |

---

### 42. Extended Annotation (3 units)

#### ANNOT-CODING-001: Coding Potential Calculation

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAnnotator.CalculateCodingPotential(...)` |
| **Complexity** | O(n) |
| **Class** | GenomeAnnotator |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateCodingPotential(sequence)` | GenomeAnnotator | Canonical |

---

#### ANNOT-REPEAT-001: Repetitive Element Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAnnotator.FindRepetitiveElements(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindRepetitiveElements(sequence, minLength)` | GenomeAnnotator | Canonical |
| `ClassifyRepeat(sequence, repeatDb)` | GenomeAnnotator | Classification |

---

#### ANNOT-CODONUSAGE-001: Codon Usage in Annotations

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAnnotator.GetCodonUsage(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetCodonUsage(codingSequences)` | GenomeAnnotator | Canonical |

---

### 43. Extended Restriction Analysis (1 unit)

#### RESTR-FILTER-001: Enzyme Filtering

| Field | Value |
|------|----------|
| **Canonical** | `RestrictionAnalyzer.GetEnzymesByCutLength(...)` |
| **Complexity** | O(e) |
| **Class** | RestrictionAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetEnzymesByCutLength(minLength, maxLength)` | RestrictionAnalyzer | By length |
| `GetBluntCutters()` | RestrictionAnalyzer | Blunt ends |
| `GetStickyCutters()` | RestrictionAnalyzer | Sticky ends |

---

### 44. Extended K-mer Analysis (1 unit)

#### KMER-DIST-001: K-mer Distance

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.KmerDistance(...)` |
| **Complexity** | O(n + m) |
| **Class** | KmerAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `KmerDistance(seq1, seq2, k)` | KmerAnalyzer | Euclidean distance |

---

### 45. Extended Motif Analysis (1 unit)

#### MOTIF-CONS-001: Consensus from Alignment

| Field | Value |
|------|----------|
| **Canonical** | `MotifFinder.CreateConsensusFromAlignment(...)` |
| **Complexity** | O(n × m) |
| **Class** | MotifFinder |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CreateConsensusFromAlignment(alignedSequences)` | MotifFinder | Canonical |

---

### 46. Extended Genomic Analysis (3 units)

#### GENOMIC-REPEAT-001: Repeat Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.FindLongestRepeat(...)` |
| **Complexity** | O(n²) |
| **Class** | GenomicAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindLongestRepeat(sequence)` | GenomicAnalyzer | Longest |
| `FindRepeats(sequence, minLength)` | GenomicAnalyzer | All repeats |

---

#### GENOMIC-COMMON-001: Common Region Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.FindLongestCommonRegion(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindLongestCommonRegion(seq1, seq2)` | GenomicAnalyzer | Longest |
| `FindCommonRegions(seq1, seq2, minLength)` | GenomicAnalyzer | All regions |

---

#### GENOMIC-MOTIFS-001: Known Motif Search

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.FindKnownMotifs(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindKnownMotifs(sequence, motifSet)` | GenomicAnalyzer | Canonical |

---

### 47. RNA Complement (1 unit)

#### SEQ-RNACOMP-001: RNA-specific Complement

| Field | Value |
|------|----------|
| **Canonical** | `SequenceExtensions.GetRnaComplementBase(...)` |
| **Complexity** | O(1) per base |
| **Class** | SequenceExtensions |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetRnaComplementBase(char)` | SequenceExtensions | Canonical |

---

### 48. Extended Protein Motif Analysis (6 units)

#### PROTMOTIF-PATTERN-001: Pattern Matching Methods

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.FindMotifByPattern(...)` |
| **Complexity** | O(n × m) |
| **Class** | ProteinMotifFinder |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindMotifByPattern(sequence, pattern)` | ProteinMotifFinder | Regex pattern |
| `FindMotifByProsite(sequence, prositePattern)` | ProteinMotifFinder | PROSITE pattern |
| `ConvertPrositeToRegex(prositePattern)` | ProteinMotifFinder | Conversion |
| `FindDomains(sequence)` | ProteinMotifFinder | Domain detection |

---

#### PROTMOTIF-SP-001: Signal Peptide Prediction

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.PredictSignalPeptide(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictSignalPeptide(sequence)` | ProteinMotifFinder | Canonical |

---

#### PROTMOTIF-TM-001: Transmembrane Helix Prediction

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.PredictTransmembraneHelices(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictTransmembraneHelices(sequence)` | ProteinMotifFinder | Hydropathy analysis |

---

#### PROTMOTIF-CC-001: Coiled-Coil Prediction

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.PredictCoiledCoils(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictCoiledCoils(sequence)` | ProteinMotifFinder | Heptad repeat analysis |

---

#### PROTMOTIF-LC-001: Low Complexity Regions

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.FindLowComplexityRegions(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindLowComplexityRegions(sequence)` | ProteinMotifFinder | Compositionally biased |

---

#### PROTMOTIF-COMMON-001: Common Motif Finding ☑

| Field | Value |
|------|----------|
| **Status** | ☑ Complete |
| **Canonical** | `ProteinMotifFinder.FindCommonMotifs(...)` |
| **Complexity** | O(n × m) |
| **Evidence** | [Evidence](docs/Evidence/PROTMOTIF-COMMON-001-Evidence.md) |
| **TestSpec** | [TestSpec](tests/TestSpecs/PROTMOTIF-COMMON-001.md) |
| **Tests** | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindCommonMotifs_Tests.cs) |
| **Algorithm Doc** | [Algorithm](docs/algorithms/ProteinMotif/Common_Motif_Finding.md) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindCommonMotifs(sequence)` | ProteinMotifFinder | Canonical |

> Note: `FindAllKnownMotifs` was listed in the original method table but does not exist in
> the codebase and has no authoritative basis distinct from `FindCommonMotifs`; the canonical
> evidence-derived method is `FindCommonMotifs`. PROSITE patterns (PS00001/05/06/16/17, …)
> verified against https://prosite.expasy.org/ (2026-06-14).

**Edge Cases:**
- [x] null / empty input → empty result
- [x] Proline at `{P}` excluded position rejects N-glycosylation site (PS00001)
- [x] multiple distinct pattern types aggregated from one sequence
- [x] multiple occurrences of one pattern reported
- [x] case-insensitive matching

---

### 49. Extended RNA Secondary Structure (7 units)

#### RNA-PAIR-001: RNA Base Pairing

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.CanPair(...)` |
| **Complexity** | O(1) |
| **Class** | RnaSecondaryStructure |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CanPair(base1, base2)` | RnaSecondaryStructure | Canonical |
| `GetBasePairType(base1, base2)` | RnaSecondaryStructure | Pair type |
| `GetComplement(base)` | RnaSecondaryStructure | RNA complement |

**Status:** ☑ Complete — Watson-Crick A-U/G-C (Watson & Crick) + G:U wobble (Crick 1966)
classified as a distinct type; RNA complement per IUPAC-IUB (1970)/Biopython complement_rna.
CanPair/GetBasePairType are O(1) symmetric table lookups over the RNA alphabet {A,C,G,U}.
See RNA-PAIR-001 Evidence/TestSpec.

---

#### ☑ RNA-HAIRPIN-001: Hairpin Energy Calculation

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.CalculateHairpinLoopEnergy(...)` |
| **Complexity** | O(n) |
| **Status** | ☑ Complete — Turner 2004 NNDB; Evidence/TestSpec/Tests linked in Registry |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateHairpinLoopEnergy(loop)` | RnaSecondaryStructure | Hairpin loop |
| `CalculateStemEnergy(stem)` | RnaSecondaryStructure | Stem region |

---

#### RNA-MFE-001: Minimum Free Energy

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.CalculateMinimumFreeEnergy(...)` |
| **Complexity** | O(n³) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateMinimumFreeEnergy(sequence)` | RnaSecondaryStructure | MFE calculation |
| `PredictStructure(sequence)` | RnaSecondaryStructure | Structure prediction |

---

#### RNA-PSEUDOKNOT-001: Pseudoknot Detection

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.DetectPseudoknots(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectPseudoknots(basePairs)` | RnaSecondaryStructure | Canonical |

---

#### RNA-DOTBRACKET-001: Dot-Bracket Notation

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.ParseDotBracket(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ParseDotBracket(notation)` | RnaSecondaryStructure | Parse |
| `ValidateDotBracket(notation)` | RnaSecondaryStructure | Validation |

---

#### RNA-INVERT-001: RNA Inverted Repeats ☑

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.FindInvertedRepeats(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindInvertedRepeats(sequence)` | RnaSecondaryStructure | Potential stems |

---

#### RNA-PARTITION-001: Partition Function ☑

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.CalculatePartitionFunction(...)` |
| **Complexity** | O(n³) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculatePartitionFunction(sequence, basePairEnergy, temperature)` | RnaSecondaryStructure | McCaskill partition function Z + base-pair probabilities |
| `CalculateStructureProbability(structureEnergy, ensembleEnergy, temperature)` | RnaSecondaryStructure | Boltzmann structure probability `p = exp(−βE)/Z` |
| `GenerateRandomRna(length[, random], gcContent)` | RnaSecondaryStructure | Random RNA generation (seeded overload deterministic) |

> **Conflict note (external evidence wins):** the original checklist listed the canonical
> as `CalculateStructureProbability(sequence)`, but the implemented API takes energies, not
> a sequence. The substantive O(n³) partition function is implemented as
> `CalculatePartitionFunction(sequence)` per McCaskill (1990). Energy model is the simplified
> fixed-per-pair model (Freiburg teaching tool); see Evidence/TestSpec.
>
> Edge cases covered: - [x] null sequence - [x] empty sequence - [x] no admissible pair
> - [x] pair span below min loop - [x] probabilities in [0,1] - [x] Z ≥ 1 invariant
> - [x] monotonicity in E_bp - [x] non-positive temperature - [x] property-based invariants
> - [x] performance baseline (n=300)

---

### 50. Extended Sequence Complexity (4 units)

#### SEQ-COMPLEX-KMER-001: K-mer Entropy

| Field | Value |
|------|----------|
| **Canonical** | `SequenceComplexity.CalculateKmerEntropy(...)` |
| **Complexity** | O(n × k) |
| **Class** | SequenceComplexity |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateKmerEntropy(sequence, k)` | SequenceComplexity | Shannon entropy using k-mers |

---

#### SEQ-COMPLEX-WINDOW-001: Windowed Complexity

| Field | Value |
|------|----------|
| **Canonical** | `SequenceComplexity.CalculateWindowedComplexity(...)` |
| **Complexity** | O(n × w) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateWindowedComplexity(sequence, windowSize)` | SequenceComplexity | Sliding window |

---

#### SEQ-COMPLEX-DUST-001: DUST Score

| Field | Value |
|------|----------|
| **Canonical** | `SequenceComplexity.CalculateDustScore(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateDustScore(DnaSequence)` | SequenceComplexity | For DnaSequence |
| `CalculateDustScore(string)` | SequenceComplexity | For string |

---

#### SEQ-COMPLEX-COMPRESS-001: Compression Ratio

| Field | Value |
|------|----------|
| **Canonical** | `SequenceComplexity.EstimateCompressionRatio(...)` (= normalized Lempel–Ziv complexity) |
| **Complexity** | O(n) |

> Note: behavior is defined as the **Lempel–Ziv (1976) complexity** (number of exhaustive-history factors), not an unspecified compressor ratio. The prior heuristic was replaced; `EstimateCompressionRatio` now delegates to `CalculateNormalizedLempelZivComplexity`. See `docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md`.

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateLempelZivComplexity(sequence)` | SequenceComplexity | Raw LZ76 component count |
| `CalculateNormalizedLempelZivComplexity(sequence)` | SequenceComplexity | Length-normalized LZ complexity |
| `EstimateCompressionRatio(sequence)` | SequenceComplexity | Compression-based complexity (delegates to normalized LZ) |

---

### 51. Extended Comparative Genomics (6 units)

#### ☑ COMPGEN-RBH-001: Reciprocal Best Hits

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.FindReciprocalBestHits(...)` |
| **Complexity** | O(n × m) |
| **Class** | ComparativeGenomics |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindReciprocalBestHits(genes1, genes2, minIdentity, minCoverage)` | ComparativeGenomics | Canonical (RBH ortholog identification) |
| `FindOrthologs(...)` | ComparativeGenomics | Delegate (delegates to RBH) |

**Note:** No prior by-area Edge Cases list existed; canonical method, signature, invariants, and edge cases were derived from the Evidence (Moreno-Hagelsieb & Latimer 2008; Tatusov et al. 1997). The pre-existing `FindReciprocalBestHits` was **nonconforming** (no coverage gate, no deterministic tie-break, placeholder `Coverage=1.0`/`AlignmentLength=0`, `Identity` set to the score product); it was corrected to the canonical RBH (matching the corrected `FindOrthologs` from COMPGEN-ORTHO-001), and `FindOrthologs` now delegates to it so the two entry points cannot diverge. RBH here is the same algorithm as COMPGEN-ORTHO-001's RBH criterion (validate + fix + document + test this dedicated entry point).

**Edge Cases:**
- [x] Null genome list → `ArgumentNullException`
- [x] Empty genome → no pairs
- [x] Gene without sequence → skipped
- [x] One-directional best hit → excluded (reciprocity)
- [x] Sub-threshold pair → excluded (coverage/identity gate)
- [x] Actual hit metrics reported (not placeholders)
- [x] Deterministic / order-independent

---

#### COMPGEN-COMPARE-001: Genome Comparison

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.CompareGenomes(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CompareGenomes(genome1, genome2)` | ComparativeGenomics | Comprehensive comparison |

---

#### COMPGEN-REVERSAL-001: Reversal Distance

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.CalculateReversalDistance(...)` |
| **Complexity** | O(n log n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateReversalDistance(geneOrder1, geneOrder2)` | ComparativeGenomics | Canonical (unsigned breakpoint lower bound) |

**Note:** No prior by-area Edge Cases list existed; canonical method, invariants, and edge cases were derived from the Evidence (Bafna & Pevzner 1998 §2; Hunter College CompBio Lecture 16; Hübotter 2020). The method returns the **unsigned breakpoint lower bound** ⌈b/2⌉, NOT the exact signed Hannenhalli–Pevzner distance — documented as **Simplified**. The pre-existing implementation was correct as the unsigned bound but its XML doc wrongly claimed "signed permutations"; the doc was corrected and source-cited constants added. Actual complexity is **O(n)** (single pass + one hash map), not the O(n log n) listed above.

**Edge Cases:**
- [x] Identity (perm1 == perm2) → 0
- [x] Empty / single-element → 0
- [x] Fully reversed → 1 (b=2)
- [x] Unequal lengths → `ArgumentException`
- [x] Arbitrary (non 1..n) labels → relabelled to relative order
- [x] Symmetry d(α,β) = d(β,α)
- [x] Lower-bound property (result ≤ reversals actually applied)
- [x] Deterministic / order-independent

---

#### COMPGEN-CLUSTER-001: Conserved Gene Clusters ☑

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.FindConservedClusters(...)` |
| **Complexity** | O(n² × g) |
| **Model** | Common intervals of permutations (Uno & Yagiura 2000; Heber & Stoye 2001; def. per Bui-Xuan, Habib & Paul 2013) — a cluster is a set of ortholog groups contiguous in every genome |
| **Evidence** | [COMPGEN-CLUSTER-001-Evidence.md](docs/Evidence/COMPGEN-CLUSTER-001-Evidence.md) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindConservedClusters(genomes)` | ComparativeGenomics | Multi-genome clusters (common intervals) |

**Edge Cases:**
- [x] Fewer than two genomes → empty (family notion, K ≥ 2)
- [x] Set contiguous in some but not all genomes → excluded
- [x] Foreign ortholog group inside the window → breaks the cluster
- [x] Repeated group labels (paralogs) → any matching window counts
- [x] `minClusterSize` filters smaller clusters; null arguments → ArgumentNullException

---

#### COMPGEN-ANI-001: Average Nucleotide Identity

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.CalculateANI(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateANI(genome1, genome2)` | ComparativeGenomics | ANI calculation |

---

#### COMPGEN-DOTPLOT-001: Dot Plot Generation

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.GenerateDotPlot(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GenerateDotPlot(seq1, seq2)` | ComparativeGenomics | Dot plot visualization |

---

### 52. Extended Motif Finding (4 units)

#### MOTIF-DISCOVER-001: Motif Discovery

| Field | Value |
|------|----------|
| **Canonical** | `MotifFinder.DiscoverMotifs(...)` |
| **Complexity** | O(n × k) |
| **Class** | MotifFinder |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DiscoverMotifs(sequences, k)` | MotifFinder | Overrepresented k-mers |

---

#### MOTIF-SHARED-001: Shared Motifs

| Field | Value |
|------|----------|
| **Canonical** | `MotifFinder.FindSharedMotifs(...)` |
| **Complexity** | O(n × m × s) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSharedMotifs(sequences)` | MotifFinder | Common motifs |

---

#### MOTIF-REGULATORY-001: Regulatory Elements

| Field | Value |
|------|----------|
| **Canonical** | `MotifFinder.FindRegulatoryElements(...)` |
| **Complexity** | O(n × r) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindRegulatoryElements(sequence)` | MotifFinder | Known regulatory motifs |

---

#### MOTIF-GENERATE-001: Consensus Generation

| Field | Value |
|------|----------|
| **Status** | ☑ Complete |
| **Canonical** | `MotifFinder.GenerateConsensus(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GenerateConsensus(alignedSequences)` | MotifFinder | Consensus sequence |

---

### 53. Extended K-mer Analysis (6 units)

#### KMER-ASYNC-001: Asynchronous K-mer Counting

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.CountKmersAsync(...)` |
| **Complexity** | O(n) |
| **Class** | KmerAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CountKmersAsync(sequence, k)` | KmerAnalyzer | Async counting |

---

#### KMER-UNIQUE-001: Unique K-mers

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.FindUniqueKmers(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindUniqueKmers(sequence, k)` | KmerAnalyzer | Canonical |

---

#### KMER-GENERATE-001: K-mer Generation

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.GenerateAllKmers(...)` |
| **Complexity** | O(4^k) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GenerateAllKmers(k)` | KmerAnalyzer | All possible k-mers |

---

#### KMER-BOTH-001: Both Strand Analysis

| Field | Value |
|------|----------|
| **Status** | ☑ Complete |
| **Canonical** | `KmerAnalyzer.CountKmersBothStrands(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CountKmersBothStrands(sequence, k)` | KmerAnalyzer | Forward + reverse complement |

---

#### KMER-STATS-001: K-mer Statistics

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.AnalyzeKmers(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeKmers(sequence, k)` | KmerAnalyzer | Comprehensive statistics |

---

#### KMER-POSITIONS-001: K-mer Positions

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.FindKmerPositions(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindKmerPositions(sequence, kmer)` | KmerAnalyzer | Position finding |

---

### 54. Extended GC Skew Analysis (3 units)

#### SEQ-ATSKEW-001: AT Skew

| Field | Value |
|------|----------|
| **Canonical** | `GcSkewCalculator.CalculateAtSkew(...)` |
| **Complexity** | O(n) |
| **Class** | GcSkewCalculator |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateAtSkew(sequence)` | GcSkewCalculator | AT skew calculation |

---

#### SEQ-REPLICATION-001: Replication Origin Prediction

| Field | Value |
|------|----------|
| **Canonical** | `GcSkewCalculator.PredictReplicationOrigin(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictReplicationOrigin(sequence)` | GcSkewCalculator | Origin and terminus |

---

#### SEQ-GC-ANALYSIS-001: Comprehensive GC Analysis

| Field | Value |
|------|----------|
| **Canonical** | `GcSkewCalculator.AnalyzeGcContent(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeGcContent(sequence)` | GcSkewCalculator | GC skew, content, variability |

---

### 55. Extended Disorder Prediction (3 units)

#### DISORDER-MORF-001: MoRF Prediction

| Field | Value |
|------|----------|
| **Canonical** | `DisorderPredictor.PredictMoRFs(...)` |
| **Complexity** | O(n) |
| **Class** | DisorderPredictor |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictMoRFs(sequence)` | DisorderPredictor | Molecular Recognition Features |

---

#### DISORDER-PROPENSITY-001: Disorder Propensity

| Field | Value |
|------|----------|
| **Status** | ☑ Complete |
| **Canonical** | `DisorderPredictor.GetDisorderPropensity(...)` |
| **Complexity** | O(1) |
| **Evidence** | [DISORDER-PROPENSITY-001-Evidence.md](docs/Evidence/DISORDER-PROPENSITY-001-Evidence.md) |
| **TestSpec** | [DISORDER-PROPENSITY-001.md](tests/TestSpecs/DISORDER-PROPENSITY-001.md) |
| **Tests** | [DisorderPredictor_GetDisorderPropensity_Tests.cs](tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_GetDisorderPropensity_Tests.cs) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetDisorderPropensity(aminoAcid)` | DisorderPredictor | Propensity value (TOP-IDP, Campen 2008) |
| `IsDisorderPromoting(aminoAcid)` | DisorderPredictor | Boolean check (Dunker 2001) |
| `DisorderPromotingAminoAcids` | DisorderPredictor | Property (8 AA, Dunker 2001) |
| `OrderPromotingAminoAcids` | DisorderPredictor | Property (8 AA, Dunker 2001) |

---

#### DISORDER-LC-001: Low Complexity in Disorder

| Field | Value |
|------|----------|
| **Canonical** | `DisorderPredictor.PredictLowComplexityRegions(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictLowComplexityRegions(sequence)` | DisorderPredictor | Protein low complexity |

---

### 56. Extended Sequence Statistics (6 units)

#### SEQ-COMPOSITION-001: Sequence Composition ☑

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateNucleotideComposition(...)` |
| **Complexity** | O(n) |
| **Class** | SequenceStatistics |
| **Status** | ☑ Complete — duplicate of SEQ-STATS-001; resolved by consolidation (no new code/tests; same canonical fixture). See `tests/TestSpecs/SEQ-COMPOSITION-001.md` §7. |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateNucleotideComposition(sequence)` | SequenceStatistics | Nucleotide composition |
| `CalculateAminoAcidComposition(sequence)` | SequenceStatistics | Amino acid composition |

---

#### SEQ-TM-001: Melting Temperature ☑

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateMeltingTemperature(...)` |
| **Complexity** | O(n) |
| **Status** | ☑ Complete — duplicate of SEQ-THERMO-001; resolved by consolidation (no new code/tests; same canonical fixture). See `tests/TestSpecs/SEQ-TM-001.md` §7. |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateMeltingTemperature(sequence)` | SequenceStatistics | Wallace/GC formula |
| `CalculateThermodynamics(sequence)` | SequenceStatistics | Thermodynamic properties |

---

#### SEQ-ENTROPY-PROFILE-001: Entropy Profile

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateEntropyProfile(...)` |
| **Complexity** | O(n × w) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateEntropyProfile(sequence, windowSize)` | SequenceStatistics | Sliding window entropy |

---

#### SEQ-GC-PROFILE-001: GC Content Profile

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateGcContentProfile(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateGcContentProfile(sequence, windowSize)` | SequenceStatistics | GC in windows |

---

#### SEQ-CODON-FREQ-001: Codon Frequencies

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateCodonFrequencies(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateCodonFrequencies(sequence)` | SequenceStatistics | Codon usage |

---

#### SEQ-SUMMARY-001: Sequence Summary

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.SummarizeNucleotideSequence(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `SummarizeNucleotideSequence(sequence)` | SequenceStatistics | Comprehensive summary |

---

### 57. Extended Genomic Analysis (3 units)

#### GENOMIC-TANDEM-001: Tandem Repeat Detection ☑

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.FindTandemRepeats(...)` |
| **Complexity** | O(n²) |
| **Class** | GenomicAnalyzer |
| **Status** | ☑ Complete — duplicate of REP-TANDEM-001 (same method/class); resolved by consolidation (no new code/tests; same canonical fixture). See `tests/TestSpecs/GENOMIC-TANDEM-001.md` §7. |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindTandemRepeats(sequence)` | GenomicAnalyzer | Consecutive repeating units |

---

#### GENOMIC-SIMILARITY-001: Sequence Similarity ☑

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.CalculateSimilarity(...)` |
| **Complexity** | O(n + m) |
| **Metric** | k-mer Jaccard index `J=|A∩B|/|A∪B|` over distinct-k-mer sets, ×100 (Jaccard 1901; Ondov et al. 2016) |
| **Status** | ☑ Complete — see [Evidence](docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md), [TestSpec](tests/TestSpecs/GENOMIC-SIMILARITY-001.md). |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateSimilarity(seq1, seq2, k)` | GenomicAnalyzer | K-mer Jaccard similarity |

**Edge Cases:**
- [x] Identical sequences → 100.0 (J=1)
- [x] Disjoint k-mer sets → 0.0 (J=0)
- [x] Both empty / shorter than k → 0.0 (empty union; Jaccard undefined, impl convention)
- [x] Repeated k-mers counted once (distinct-set semantics)
- [x] Null sequence → ArgumentNullException; kmerSize < 1 → ArgumentOutOfRangeException

---

#### GENOMIC-ORF-001: ORF Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.FindOpenReadingFrames(...)` |
| **Complexity** | O(n²) worst case (per-ATG scan to next in-frame stop; O(n) typical) — corrected from O(n); see Evidence GENOMIC-ORF-001 |
| **Invariant** | Every ORF starts ATG, ends TAA/TAG/TGA, length % 3 == 0; six-frame search |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindOpenReadingFrames(sequence, minLength=100)` | GenomicAnalyzer | Canonical — every-ATG-to-first-in-frame-stop, both strands |

**Edge Cases:**
- [x] ATG with no in-frame stop → not reported
- [x] Nested ORFs sharing a stop → both reported (Rosalind canonical)
- [x] Reverse-complement-only ORF detected
- [x] minLength filtering (nucleotides, inclusive)
- [x] null sequence → ArgumentNullException

---

### 17. Oncology Genomics (35 units)

#### ONCO-SOMATIC-001: Somatic Mutation Calling

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.CallSomaticMutations(...)` |
| **Complexity** | O(n) |
| **Invariant** | Tumor VAF > threshold, absent in matched normal |
| **Depends on** | — (standalone, entry point) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CallSomaticMutations(tumorVcf, normalVcf)` | OncologyAnalyzer | Canonical |
| `FilterGermlineVariants(variants, normalSample)` | OncologyAnalyzer | Filter |
| `CalculateSomaticScore(variant)` | OncologyAnalyzer | Scoring |

**Edge Cases:**
- [x] Tumor-only mode (no matched normal)
- [x] Low tumor purity samples
- [x] Clonal hematopoiesis contamination

---

#### ONCO-VAF-001: Variant Allele Frequency Analysis ☑

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.CalculateVAF(...)` |
| **Complexity** | O(1) per variant |
| **Invariant** | 0 ≤ VAF ≤ 1 |
| **Depends on** | ONCO-SOMATIC-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateVAF(altReads, totalReads)` | OncologyAnalyzer | Canonical |
| `CalculateVAFConfidenceInterval(altReads, totalReads, confidence)` | OncologyAnalyzer | Canonical (Wilson 1927) |
| `AdjustVAFForPurity(vaf, purity, ploidy)` | OncologyAnalyzer | Correction |

**Edge Cases:**
- [x] totalReads = 0 (division by zero)
- [x] Multiallelic sites (multiple ALT alleles) — per-allele alt/total ratio is well-defined for each ALT
- [x] VAF > 1.0 due to alignment artifacts (altReads > totalReads → ArgumentOutOfRangeException)

---

#### ☑ ONCO-DRIVER-001: Driver Mutation Detection

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.IdentifyDriverMutations(...)` |
| **Complexity** | O(n × k) where k=known drivers |
| **Invariant** | driver_mutations ⊆ somatic_mutations |
| **Depends on** | ONCO-SOMATIC-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `IdentifyDriverMutations(variants)` | OncologyAnalyzer | Canonical |
| `ScoreDriverPotential(variant)` | OncologyAnalyzer | CADD/SIFT/PolyPhen |
| `MatchCancerHotspots(variant)` | OncologyAnalyzer | Database lookup |

**Databases:** COSMIC, OncoKB, ClinVar, Cancer Hotspots

---

#### ONCO-ARTIFACT-001: Sequencing Artifact Detection

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.FilterArtifacts(...)` |
| **Complexity** | O(n) |
| **Depends on** | ONCO-SOMATIC-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FilterArtifacts(variants, bamFile)` | OncologyAnalyzer | Canonical |
| `DetectOxoGArtifacts(variants)` | OncologyAnalyzer | 8-oxoG filter |

---

#### ONCO-ANNOT-001: Cancer-Specific Variant Annotation ☑

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.AnnotateCancerVariants(...)` / `ClassifyVariantTier(...)` |
| **Complexity** | O(n) per batch; O(1) per variant |
| **Standard** | AMP/ASCO/CAP 2017 four-tier system (Li et al. 2017, J Mol Diagn 19(1):4–23) |
| **Depends on** | ONCO-SOMATIC-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ClassifyVariantTier(variant)` | OncologyAnalyzer | Canonical |
| `AnnotateCancerVariants(variants)` | OncologyAnalyzer | Canonical |
| `GetCOSMICAnnotation(variant, catalog)` | OncologyAnalyzer | COSMIC lookup (caller-supplied catalog) |

**Edge Cases:**
- [x] Clinical evidence (Level A/B) takes priority over high population MAF (still Tier I)
- [x] Population MAF ≥ 1% with no clinical evidence ⇒ Tier IV (benign cutoff, Li 2017)
- [x] Rare variant, no evidence: Tier III if cancer association present, else Tier IV
- [x] MAF boundary exactly 0.01 (inclusive ≥) ⇒ Tier IV
- [x] Invalid MAF (negative / > 1 / NaN) throws; null inputs throw; empty batch ⇒ empty
- [x] COSMIC catalog miss ⇒ null (external content not fabricated)

---

#### ☑ ONCO-TMB-001: Tumor Mutational Burden

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.CalculateTMB(...)` |
| **Complexity** | O(n) |
| **Invariant** | TMB = mutations / Mb (coding region) |
| **Depends on** | ONCO-SOMATIC-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateTMB(mutationCount, targetRegionMb)` / `CalculateTMB(calls, targetRegionMb)` | OncologyAnalyzer | Canonical |
| `ClassifyTMB(tmb)` | OncologyAnalyzer | Low/High (FDA ≥10) |

**Thresholds:** TMB-High = TMB ≥ 10 mut/Mb (inclusive), else Low. Source: FDA pembrolizumab approval / FoundationOne CDx (Marcus et al. 2021, Clin Cancer Res 27(17):4685). NOTE: the previous "Low (<6/Mb), Intermediate (6–20/Mb), High (>20/Mb)" boundaries had no retrievable authoritative source and were replaced by the source-backed FDA ≥10 cutoff (evidence-first; see TestSpec §7).

**Edge Cases:**
- [x] targetRegionMb = 0 (division by zero — throws)
- [x] Panel < 0.5 Mb (value still computed; instability documented, not an error — Chalmers 2017)
- [x] Synonymous/germline variants excluded by upstream somatic caller (count is pre-filtered; SomaticCall overload counts only Somatic status)

---

#### ONCO-MSI-001: Microsatellite Instability Detection

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.DetectMSI(...)` |
| **Complexity** | O(n × m) where m=microsatellite loci |
| **Invariant** | 0 ≤ MSI_score ≤ 1 |
| **Depends on** | — (standalone) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectMSI(tumorBam, normalBam, loci)` | OncologyAnalyzer | Canonical |
| `CalculateMSIScore(instableLoci, totalLoci)` | OncologyAnalyzer | Score |
| `ClassifyMSIStatus(score)` | OncologyAnalyzer | MSS/MSI-L/MSI-H |

**Markers:** BAT25, BAT26, NR21, NR24, MONO27 (Bethesda panel)

**Edge Cases:**
- [x] Tumor-only mode (no matched normal) — fraction-of-unstable-loci definition is normal-independent; DetectMSI takes per-locus flags
- [x] Insufficient coverage at microsatellite loci — zero valid loci → CalculateMSIScore/DetectMSI throw (score undefined)
- [x] < 5 evaluable markers (unreliable classification) — ClassifyBethesdaPanel validates 0 ≤ unstable ≤ total; MSS/MSI-L ambiguity documented per Boland 1998

---

#### ONCO-HRD-001: Homologous Recombination Deficiency

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.CalculateHRDScore(...)` |
| **Complexity** | O(n) |
| **Invariant** | HRD = LOH + TAI + LST |
| **Depends on** | ONCO-LOH-001, ONCO-CNA-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateHRDScore(cnvSegments)` | OncologyAnalyzer | Canonical |
| `CalculateLOHScore(segments)` | OncologyAnalyzer | LOH component |
| `CalculateTAIScore(segments)` | OncologyAnalyzer | Telomeric AI |
| `CalculateLSTScore(segments)` | OncologyAnalyzer | Large-scale transitions |

**Clinical:** PARP inhibitor response prediction (BRCA1/2, PALB2)

**Edge Cases:**
- [x] Insufficient number of CNV segments for scoring — component counts supplied as input; sum is well-defined for any non-negative counts (scope per ONCO-HRD-001 NOTE: composite-sum + threshold)
- [x] Near-diploid tumors (low signal) — all-zero components → score 0 → HRD-negative (DetectHRD_NearDiploidZeroComponents_IsHrdNegative)

---

#### ONCO-LOH-001: Loss of Heterozygosity

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.DetectLOH(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ LOH_fraction ≤ 1 per chromosome |
| **Depends on** | ONCO-SOMATIC-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectLOH(IEnumerable<AlleleSpecificSegment>)` | OncologyAnalyzer | Canonical — qualifying HRD-LOH regions + score |
| `CalculateHrdLohScore(IEnumerable<AlleleSpecificSegment>)` | OncologyAnalyzer | HRD-LOH count (Abkevich 2012 / scarHRD) |
| `CalculateLOHFraction(IEnumerable<AlleleSpecificSegment>, chromosome)` | OncologyAnalyzer | Length-weighted per-chromosome LOH fraction |

> **Validation note (ONCO-LOH-001):** the HRD-LOH algorithm of Abkevich et al. (2012) / scarHRD / oncoscanR
> is defined over *allele-specific copy-number segments* (chromosome, start, end, major CN, minor CN), not
> raw `(tumorVcf, normalVcf)` text — raw segmentation/BAF is upstream (ONCO-CNA-001). The canonical signature
> is therefore `DetectLOH(IEnumerable<AlleleSpecificSegment>)`. (Earlier `DetectLOH(tumorVcf, normalVcf)` stub
> corrected 2026-06-16.)

---

#### ONCO-SIG-001: COSMIC Mutational Signature Extraction

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.ClassifySbsContext(...)` / `Build96ContextCatalog(...)` |
| **Complexity** | O(n) to build the 96-channel SBS catalog (n = #variants) |
| **Invariant** | Σ catalog counts = #classifiable SBS variants; 96 channels; pyrimidine reference base |
| **Depends on** | — (foundational 96-context catalog is self-contained) |

**Implemented (this unit):** the foundational, well-defined and authoritatively retrievable piece — the
96-class SBS trinucleotide-context catalog: each single-base substitution classified by the 6 pyrimidine
substitution types (C>A,C>G,C>T,T>A,T>C,T>G) × 4 5'-bases × 4 3'-bases, with reverse-complement folding of
purine-reference mutations onto the pyrimidine strand (Alexandrov et al. 2013, Nature 500:415-421; COSMIC
SBS96; SigProfilerMatrixGenerator, Bergstrom et al. 2019).

**Methods (implemented):**
| Method | Class | Type |
|-------|-------|-----|
| `ClassifySbsContext(fivePrime, ref, alt, threePrime)` | OncologyAnalyzer | Canonical |
| `Build96ContextCatalog(variants)` | OncologyAnalyzer | Canonical |
| `EnumerateSbs96Channels()` | OncologyAnalyzer | Canonical |

**Deferred to later units:** NMF signature extraction / COSMIC fitting (`ExtractSignatures`,
`FitToCosmicSignatures`, `DecomposeSpectrum`) require caller-supplied signature reference matrices and are NOT
implemented here (reference profiles must not be fabricated). Signatures: SBS (DBS, ID out of scope).

---

#### ONCO-SIG-002: Mutational Signature Fitting / Refitting ☑

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.FitSignatures(...)` / `OncologyAnalyzer.CosineSimilarity(...)` |
| **Complexity** | NNLS active-set: O(k³ + k²·n) per outer iteration, ≤ O(k) outer iterations (k signatures, n channels); cosine O(n) |
| **Invariant** | all exposures ≥ 0; cosine ∈ [0,1]; residual SSE ≤ ‖d‖²; normalized exposures sum to 1 |
| **Depends on** | ONCO-SIG-001 (96-context catalog supplies the observed `catalog` vector) |

**Scope note:** the by-area registry originally titled this unit "Trinucleotide Context Analysis", but the
trinucleotide-context catalog (`ClassifySbsContext` / `Build96ContextCatalog`) was already implemented under
ONCO-SIG-001. The genuinely next mutational-signature piece — and the one implemented here — is signature
**fitting/refitting**: decomposing an observed 96-channel catalog into a non-negative combination of
caller-supplied reference signatures via NNLS (min‖S·x−d‖², x≥0; Blokzijl et al. 2018; Lawson & Hanson 1974)
plus cosine similarity between catalogs/signatures (Blokzijl et al. 2018; iMutSig). Reference signature
profiles are caller-supplied (not fabricated). External evidence supersedes the original by-area title.

**Methods (implemented):**
| Method | Class | Type |
|-------|-------|-----|
| `CosineSimilarity(a, b)` | OncologyAnalyzer | Canonical |
| `FitSignatures(catalog, signatures)` | OncologyAnalyzer | Canonical |
| `ReconstructCatalog(signatures, exposures)` | OncologyAnalyzer | Canonical |

**Edge Cases:**
- [x] Zero observed catalog ⇒ all exposures 0, reconstruction 0
- [x] Zero-norm vector in cosine similarity (÷0) ⇒ 0.0
- [x] Unconstrained LS coefficient < 0 ⇒ clamped to 0, refit on remaining set
- [x] Null / empty / dimension-mismatched inputs ⇒ ArgumentNullException / ArgumentException

---

#### ONCO-SIG-003: Signature Exposure Estimation — Bootstrap Confidence Intervals ☑

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.BootstrapExposures(...)` |
| **Complexity** | O(R·(N + NNLS)) — R replicates, each a multinomial draw of N mutations + one NNLS refit |
| **Invariant** | all(exposure ≥ 0); lower ≤ upper; deterministic for a fixed seed |
| **Depends on** | ONCO-SIG-001 (catalog) and ONCO-SIG-002 (NNLS `FitSignatures` refit) |

**Scope note:** the by-area registry titled this unit "Signature Exposure Estimation" with a canonical
`EstimateExposures`, but point exposure estimation (NNLS) was already delivered by ONCO-SIG-002
(`FitSignatures`). The genuinely next mutational-signature piece — and the registry's second listed method
`BootstrapConfidenceIntervals(spectrum)` — is **bootstrap confidence intervals on exposures**, implemented
here as `OncologyAnalyzer.BootstrapExposures`: a parametric multinomial bootstrap (resample the integer
catalog as Multinomial(N, p), NNLS-refit each resample, take per-signature percentile intervals). Sources:
Senkin 2021 (MSA, BMC Bioinformatics 22:540), Huang et al. 2018 (Bioinformatics 34(2):330–337), sigminer
`sig_fit_bootstrap`, Efron 1979 percentile method, Hyndman & Fan 1996 type-7 quantile. External evidence
supersedes the original by-area title. Reference signature profiles are caller-supplied (not fabricated).

**Methods (implemented):**
| Method | Class | Type |
|-------|-------|-----|
| `BootstrapExposures(catalog, signatures, replicates, confidence, seed)` | OncologyAnalyzer | Canonical |

**Edge Cases:**
- [x] Zero-mutation catalog (N = 0) ⇒ every interval [0, 0], point/mean 0
- [x] Single non-zero channel ⇒ deterministic resample ⇒ lower = upper = mean = point estimate
- [x] Single replicate (R = 1) ⇒ percentile of one-element sample ⇒ lower = upper = mean
- [x] Null / empty / ragged / dimension-mismatched / negative-count inputs ⇒ ArgumentNullException / ArgumentException
- [x] replicates < 1 or confidence ∉ (0, 1) ⇒ ArgumentOutOfRangeException

---

#### ONCO-SIG-004: Mutational Process Classification ☑

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.ClassifyMutationalProcess(...)` / `OncologyAnalyzer.GetMutationalProcess(...)` |
| **Complexity** | O(k log k) (k signatures; ordering of ≤5 processes) |
| **Invariant** | each surviving contribution ∈ [0,1], Σ ≤ 1; active iff normalized contribution ≥ 0.06 (strict `<` excludes); dominant = max per-process sum; Σ exposure = 0 ⇒ no active processes |
| **Depends on** | ONCO-SIG-003 |

**Scope note:** the by-area registry placed this in a `MutationalSignatures` class, but the area's analyzer is
`OncologyAnalyzer` (ONCO-SIG-001..003 all live there), so the methods were added there. Classification follows
the deconstructSigs presence rule — exposures are converted to normalized relative contributions and a
signature is active only when its contribution ≥ 6% (Rosenthal et al. 2016, *Genome Biology* 17:31:
"any signature with Wᵢ < 6% is excluded"; reference `whichSignatures.R` `signature.cutoff = 0.06`) — and the
COSMIC SBS→aetiology map (https://cancer.sanger.ac.uk/signatures/sbs/; Alexandrov et al. 2020, *Nature* 578).
Reference signature profiles are caller-supplied (not fabricated); only exposures + COSMIC labels are consumed.

**Methods (implemented):**
| Method | Class | Type |
|-------|-------|-----|
| `ClassifyMutationalProcess(exposures, contributionCutoff)` | OncologyAnalyzer | Canonical |
| `GetMutationalProcess(signatureLabel)` | OncologyAnalyzer | Canonical |

**Processes:** Aging (SBS1/5), APOBEC (SBS2/13), Smoking (SBS4), UV (SBS7a–d), MMR deficiency (SBS6/15/20/26)

**Edge Cases:**
- [x] All-zero / empty exposures (Σ = 0) ⇒ no active processes, dominant = Unknown
- [x] Contribution exactly 0.06 ⇒ retained (strict `<` cutoff); just below ⇒ excluded
- [x] Multiple signatures per process (SBS2+SBS13, SBS6/15/20/26, SBS7a+7b) ⇒ summed
- [x] Unmapped / unknown-aetiology SBS label ⇒ contributes to no recognized process
- [x] Null exposures / null label ⇒ ArgumentNullException; negative exposure ⇒ ArgumentException; cutoff ∉ [0,1) ⇒ ArgumentOutOfRangeException

---

#### ☑ ONCO-FUSION-001: Fusion Gene Detection

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.DetectFusions(...)` |
| **Complexity** | O(n log n) |
| **Invariant** | gene5p ≠ gene3p |
| **Depends on** | — (standalone) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectFusions(candidates, thresholds?)` | OncologyAnalyzer | Canonical |
| `IsInFrame(fivePrimeCodingBases, threePrimeStartPhase)` | OncologyAnalyzer | Canonical |
| `ComputeTotalSupport(candidate)` | OncologyAnalyzer | Internal |
| `FindChimericReads(bamFile)` | (FusionDetector) | Read extraction — out of scope (raw-BAM step) |
| `ValidateFusion(fusion, refGenome)` | (FusionDetector) | Validation — out of scope |

Implemented on `OncologyAnalyzer` (this session's mandated class) using STAR-Fusion min-support thresholds
(MIN_JUNCTION_READS=1, MIN_SUM_FRAGS=2, MIN_SPANNING_FRAGS_ONLY=5) and Arriba total-support / reading-frame
definitions. `FindChimericReads`/`ValidateFusion` are raw-BAM/validation steps outside the count-based
canonical scope.

**Clinical Fusions:** BCR-ABL, EML4-ALK, ROS1, NTRK, FGFR, RET

**Edge Cases:**
- [x] Read-through transcripts (false positive fusions) — distinct-gene rule (gene5p ≠ gene3p) + support thresholds
- [x] Low supporting read count — MIN_SUM_FRAGS / MIN_SPANNING_FRAGS_ONLY thresholds reject low-evidence candidates
- [x] Reciprocal fusions (same partners, swapped orientation) — treated as distinct candidates by 5'/3' assignment

---

#### ☑ ONCO-FUSION-002: Known Fusion Database Lookup

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.MatchKnownFusions(...)` |
| **Complexity** | O(1) hash lookup (O(k) case-insensitive fallback) |
| **Invariant** | designation = 5'::3' (directional); A::B ≠ B::A |
| **Depends on** | ONCO-FUSION-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetFusionAnnotation(gene5p, gene3p)` | OncologyAnalyzer | Canonical |
| `MatchKnownFusions(fusion, knownFusions)` | OncologyAnalyzer | Canonical |

Implemented on `OncologyAnalyzer` (same class as ONCO-FUSION-001). Designation format and directional 5'→3'
keying follow HGNC nomenclature (Bruford et al. 2021, double colon `::`, 5' partner first). Known-fusion set
membership and annotations are **caller-supplied** — the library bundles no curated database (ChimerDB / COSMIC
Fusions / Mitelman content is the caller's responsibility), making this a Framework algorithm.

**Edge Cases:**
- [x] Reciprocal fusion not matched — directional designation (A::B ≠ B::A) per 5'-first rule
- [x] Case-varied symbols — case-insensitive (ordinal-ignore-case) matching
- [x] Null/empty partner or null known-fusion set — input validation (ArgumentException / ArgumentNullException)

---

#### ☑ ONCO-FUSION-003: Fusion Breakpoint Analysis

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.AnalyzeBreakpoint(...)` |
| **Complexity** | O(n) |
| **Depends on** | ONCO-FUSION-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeBreakpoint(fusion)` | OncologyAnalyzer | Canonical |
| `PredictFusionProtein(fusion, transcripts)` | OncologyAnalyzer | Canonical (protein product) |

Implemented on `OncologyAnalyzer` (same class as ONCO-FUSION-001/002). Breakpoint site categories follow
the Arriba `site`/`reading_frame` output schema (Uhrig et al. 2021); the fusion protein is the chimeric CDS
(5' CDS prefix ++ 3' CDS suffix) translated with the standard genetic code (NCBI Table 1) and truncated at
the first stop codon, exactly per AGFusion (`model.py`; Murphy & Elemento 2016). Frame call reuses the
ONCO-FUSION-001 codon-phase rule `(b − p) mod 3 == 0`. Partner CDS sequences and breakpoint offsets are
**caller-supplied** (no genome/transcript DB in repo) — a Framework algorithm.

**Edge Cases:**
- [x] Breakpoint not in CDS (UTR/intron/intergenic) — reading frame `NotPredicted` (Arriba `reading_frame = .`)
- [x] Premature stop codon at/after junction — peptide truncated at first stop; `HasPrematureStop` flag (AGFusion `protein[0:find("*")]`)
- [x] Out-of-frame junction — chimeric CDS trimmed to whole codons; 3' partner read in shifted frame (AGFusion)

---

#### ONCO-CNA-001: Cancer Copy Number Alteration Classification ☑

> **Scope note (resolved in implementation):** The well-defined, source-traceable
> foundational piece of CNA analysis is log2 copy-ratio → absolute/integer copy number →
> discrete CNA classification (deep deletion / loss / neutral / gain / amplification),
> implemented in `OncologyAnalyzer` per CNVkit `absolute_threshold` / GISTIC2.0. Full CBS
> segmentation (`SegmentCopyNumber`) already exists in `StructuralVariantAnalyzer`
> (SV-CNV-001); the original `CopyNumberAnalyzer` / `CallCopyNumberStates` names are
> superseded by the `OncologyAnalyzer.ClassifyCopyNumber` API below. Conflict noted in
> the TestSpec §7.

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.ClassifyCopyNumber(log2Ratio, thresholds, ploidy)` |
| **Complexity** | O(1) per region; O(n) per batch |
| **Invariant** | integer CN monotonic non-decreasing in log2; CN≥0; state↔CN mapping fixed |
| **Depends on** | — (standalone; reuses CNVkit n = 2·2^log2, cited shared with SV-CNV-001) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Log2RatioToCopyNumber(log2Ratio, ploidy)` | OncologyAnalyzer | Canonical (n = ploidy·2^log2) |
| `CallCopyNumber(log2Ratio, thresholds, ploidy)` | OncologyAnalyzer | Canonical (hard-threshold integer CN) |
| `ClassifyCopyNumber(log2Ratio, thresholds, ploidy)` | OncologyAnalyzer | Canonical (5-state classification) |
| `ClassifyCopyNumbers(log2Ratios, ...)` | OncologyAnalyzer | Delegate (batch) |

**Edge Cases:**
- [x] Noisy log-ratio data (low coverage) — neutral noise band (−0.25, 0.2]; NaN no-call → Neutral
- [x] Whole-genome doubling (WGD) — `ploidy` parameter exposed (default diploid; ASM-01)

---

#### ONCO-CNA-002: Focal Amplification Detection ☑

> **Scope note (resolved in implementation):** Placed in `OncologyAnalyzer` (consistent with
> ONCO-CNA-001, which superseded the `CopyNumberAnalyzer` names). A segment is a focal amplification
> when it is amplified (log2 > GISTIC2 `t_amp` = 0.1) AND focal (length < GISTIC2 `broad_len_cutoff`
> = 0.98 × chromosome-arm length), per Mermel et al. (2011) GISTIC2.0 length-based focal/arm-level
> split. Oncogene mapping uses NCBI Gene cytogenetic arms. Class-name conflict noted in TestSpec §7.

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.DetectFocalAmplifications(segments, thresholds?)` |
| **Complexity** | O(n) |
| **Depends on** | ONCO-CNA-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectFocalAmplifications(segments, thresholds?)` | OncologyAnalyzer | Canonical (length < 0.98·arm AND log2 > t_amp) |
| `IdentifyAmplifiedOncogenes(amplifications)` | OncologyAnalyzer | Canonical (arm → ERBB2/MYC/EGFR/CCND1/MDM2/CDK4) |
| `IsFocalAmplification(segment, thresholds)` | OncologyAnalyzer | Internal predicate |

**Oncogenes:** ERBB2/HER2 (17q), MYC (8q), EGFR (7p), CCND1 (11q), MDM2 (12q), CDK4 (12q)

**Edge Cases:**
- [x] Whole-arm event (≥ 98% of arm) excluded as arm-level/broad — GISTIC2 `broad_len_cutoff` 0.98
- [x] Low-amplitude gain (log2 ≤ t_amp 0.1) not called amplified — GISTIC2 `t_amp`
- [x] Boundary exactly 0.98 of arm is arm-level (strict focal test)
- [x] Null / empty input and invalid arm length / coordinates validated

---

#### ONCO-CNA-003: Homozygous Deletion Detection ☑

> **Scope note (resolved in implementation):** Placed in `OncologyAnalyzer` (consistent with
> ONCO-CNA-001/002, which superseded the `CopyNumberAnalyzer` names). A segment is a homozygous
> (deep) deletion when its hard-threshold integer copy number is 0 — total copy number 0, both
> alleles lost — i.e. the cBioPortal "−2" Deep Deletion / DeepDeletion state, per Cheng et al.
> (2017) "zero copies of both alleles" and CNVkit `absolute_threshold` integer-CN calling (reused
> from ONCO-CNA-001). Tumour-suppressor mapping uses NCBI Gene cytogenetic arms. Reuses the
> ONCO-CNA-002 `CopyNumberArmSegment`. Class-name conflict noted in TestSpec §7.

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.DetectHomozygousDeletions(segments, thresholds?, ploidy?)` |
| **Complexity** | O(n) |
| **Invariant** | reported iff integer CN = 0; single-copy (CN 1) loss never reported; order-preserving |
| **Depends on** | ONCO-CNA-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectHomozygousDeletions(segments, thresholds?, ploidy?)` | OncologyAnalyzer | Canonical (integer CN = 0) |
| `IsHomozygousDeletion(segment, thresholds?, ploidy?)` | OncologyAnalyzer | Internal predicate |
| `IdentifyDeletedTumorSuppressors(deletions)` | OncologyAnalyzer | Canonical (arm → TP53/RB1/CDKN2A/PTEN/BRCA1/BRCA2) |

**Tumor Suppressors:** TP53 (17p), RB1 (13q), CDKN2A (9p), PTEN (10q), BRCA1 (17q), BRCA2 (13q)

**Edge Cases:**
- [x] Single-copy / heterozygous loss (CN 1, cBioPortal −1) NOT reported as homozygous
- [x] Boundary log2 exactly at deletion cutoff (−1.1) is CN 0 (homozygous); just above is CN 1
- [x] NaN log2 no-call → neutral (CN = rounded ploidy), not homozygous
- [x] Null / empty input and invalid arm length / coordinates validated
- [x] Custom thresholds and ploidy shift the CN-0 boundary

---

#### ONCO-PURITY-001: Tumor Purity Estimation ☑

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.EstimatePurity(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ purity ≤ 1 |
| **Depends on** | ONCO-VAF-001, ONCO-CNA-001 |

**Methods:** (implemented on `OncologyAnalyzer`; closed-form inversion of the CNAqc expected-VAF relation v = mπ/[2(1−π)+π·n_tot])
| Method | Class | Type |
|-------|-------|-----|
| `EstimatePurity(IEnumerable<PurityVariant>)` | OncologyAnalyzer | Canonical (allele-specific) |
| `EstimatePurityFromVAF(IEnumerable<VariantObservation>)` | OncologyAnalyzer | Canonical (diploid het, ρ = 2·VAF) |
| `EstimatePurityFromVaf(double)` | OncologyAnalyzer | Delegate (single-VAF closed form) |

**Edge Cases:**
- [x] Purity < 0.1 (below detection limit)
- [x] No heterozygous SNPs / no informative variants (empty input → undefined)
- [x] High stromal contamination (low-purity / boundary handling)

---

#### ONCO-PLOIDY-001: Tumor Ploidy Estimation ☑

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.EstimatePloidy(...)` |
| **Complexity** | O(n) |
| **Invariant** | ploidy > 0; WGD ⇔ frac(major CN ≥ 2 by length) > 0.5 |
| **Depends on** | ONCO-CNA-001 / ONCO-LOH-001 (`AlleleSpecificSegment`) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `EstimatePloidy(IEnumerable<AlleleSpecificSegment>)` | OncologyAnalyzer | Canonical — ψ = Σ(CN·L)/Σ(L), CN = Major+Minor (Patchwork/ASCAT) |
| `DetectWholeGenomeDoubling(IEnumerable<AlleleSpecificSegment>)` | OncologyAnalyzer | Canonical — facets-suite rule: frac(major CN ≥ 2 by length) > 0.5 (Bielski 2018, PMID 30013179) |

**Edge Cases:**
- [x] Empty segment set → ArgumentException (weighted mean / fraction undefined)
- [x] Segment End ≤ Start (non-positive length) → ArgumentException
- [x] Negative copy number → ArgumentException
- [x] All-1:1 genome → ψ = 2.0, WGD = false (major CN, not total)
- [x] Exactly 50% major CN ≥ 2 → WGD = false (strict `>`)

> Note: the registry stub listed `DetectWholeGenomeDoubling(ploidy)` (scalar) and class `TumorAnalyzer`. The authoritative WGD definition (major CN ≥ 2 over >50% of the genome, Bielski 2018 / facets-suite) requires per-segment data, so the canonical method takes segments; the area's analyzer class is `OncologyAnalyzer` (matching all sibling ONCO units).

---

#### ONCO-CLONAL-001: Clonal vs Subclonal Classification ☑

> Note: the registry stub listed class `TumorAnalyzer` and signature `ClassifyClonality(variants, purity, ploidy)`. The authoritative Landau et al. (2013) CCF posterior model uses the **per-locus** absolute copy number q (not a genome-wide ploidy scalar), so the canonical method is `OncologyAnalyzer.ClassifyClonality(variants, purity)` with q carried per `ClonalityVariant` (mirrors the ONCO-WGD scalar→per-segment decision). The analyzer class is `OncologyAnalyzer` (all sibling ONCO units).

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.ClassifyClonality(...)` |
| **Complexity** | O(n) |
| **Invariant** | clonal_count + subclonal_count = total_variants |
| **Depends on** | ONCO-CCF-001 (CCF; consumed by `IdentifyClonalMutations`) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ClassifyClonality(variants, purity)` | OncologyAnalyzer | Canonical — clonal iff P(CCF>0.95) > 0.5 (Landau 2013) |
| `IdentifyClonalMutations(ccfValues)` | OncologyAnalyzer | Canonical — clonal iff CCF > 0.95 |

---

#### ONCO-NEO-001: Neoantigen Prediction (candidate peptide window generation)

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.GenerateNeoantigenPeptides(...)` |
| **Complexity** | O(Σ_k k) per missense mutation (constant for the fixed 8–11 class I range) |
| **Depends on** | ONCO-SOMATIC-001 |

> Scope (implemented): the well-defined windowing step — every 8–11-mer (MHC-I) window of the mutant
> protein that SPANS a somatic missense residue, paired with the wild-type agretope at the same
> coordinates (pVACtools / Hundal 2020; ProGeo-neo / Li 2020; TESLA / Wells 2020). Binding affinity /
> IC50 scoring is caller-supplied / out of scope (ONCO-MHC-001) — no MHC model is fabricated. The
> checklist placeholder class `NeoantigenPredictor` is superseded by `OncologyAnalyzer` (project layout).

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GenerateNeoantigenPeptides(wildTypeProtein, mutantResidue, mutationPosition, minLength, maxLength)` | OncologyAnalyzer | Canonical (8-11mer windowing + agretope pairing) |

**Edge Cases:**
- [x] Mutation at N-/C-terminus (truncated window count — only windows that fit while spanning the residue)
- [x] Requested peptide length exceeds protein length (length skipped; shorter lengths still returned)
- [x] Non-substitution (mutant residue == wild-type) rejected; out-of-range position / invalid length range rejected
- [ ] Frameshift / indel / fusion neopeptides (out of scope; separate translation step)
- [ ] Binding-affinity scoring (out of scope; ONCO-MHC-001)

---

#### ONCO-MHC-001: MHC-Peptide Binding Classification

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.ClassifyBindingAffinity(...)` / `ClassifyBindingRank(...)` / `IsValidPeptideLength(...)` |
| **Complexity** | O(1) per classification |
| **Invariant** | IC50 > 0; %Rank ∈ [0,100] |
| **Depends on** | ONCO-NEO-001 |

> Scope (implemented): the **classification** of a caller-supplied predicted affinity into binder
> categories using the standard IEDB / NetMHCpan-4.1 thresholds (Reynisson 2020; Sette 1994): IC50
> strong < 50 nM / weak < 500 nM; class I %Rank strong < 0.5% / weak < 2%; class II %Rank strong < 2% /
> weak < 10%; plus peptide-length validity (class I 8–11, class II 13–25). The peptide–MHC affinity / %Rank
> **prediction** is now also available as an **opt-in matrix-based predictor**: the BIMAS / Parker 1994
> product rule (`PredictBindingHalfLifeBimas`) and the SMM / Peters & Sette 2005 transform
> `IC50 = 50000^(1−score)` (`PredictIc50Smm` / `PredictAndClassifySmm`), which chain into the existing
> classifier. The trained coefficient **matrix is caller-supplied** via `LoadScoringMatrix` — no
> redistributable, cross-verifiable trained HLA matrix was obtainable (BIMAS CGI dead/unarchived; Parker
> 1994 paywalled; IEDB SMM non-commercial) — so only the published scoring rules are bundled; no model is
> fabricated. The **pan-allele NetMHCpan/MHCflurry neural** prediction remains out of scope. The checklist
> placeholder class `NeoantigenPredictor` is superseded by `OncologyAnalyzer` (project layout).

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ClassifyBindingAffinity(ic50Nm)` | OncologyAnalyzer | Canonical — IC50 → Strong/Weak/NonBinder (50/500 nM) |
| `ClassifyBindingRank(percentRank, mhcClass)` | OncologyAnalyzer | Canonical — %Rank → strength (class I 0.5/2; class II 2/10) |
| `IsValidPeptideLength(length, mhcClass)` | OncologyAnalyzer | Canonical — class I 8–11; class II 13–25 |
| `ClassifyMhcBinding(peptideLength, ic50Nm, mhcClass)` | OncologyAnalyzer | Delegate — length gate + affinity |
| `PredictIc50Smm(peptide, matrix)` | OncologyAnalyzer | Canonical — SMM sum → IC50 = 50000^(1−score) (opt-in; matrix caller-supplied) |
| `PredictBindingHalfLifeBimas(peptide, matrix)` | OncologyAnalyzer | Canonical — BIMAS product → half-life (opt-in; matrix caller-supplied) |
| `PredictAndClassifySmm(peptide, matrix)` | OncologyAnalyzer | Delegate — `PredictIc50Smm` + `ClassifyBindingAffinity` |
| `LoadScoringMatrix(lines)` | OncologyAnalyzer | Canonical — caller-supplied matrix loader |
| `PredictMHCBinding(peptide, hlaAllele)` | (NeoantigenPredictor) | Out of scope — pan-allele neural model |

**HLA Types:** HLA-A, HLA-B, HLA-C (Class I); HLA-DR/DQ/DP (Class II length range)

**Edge Cases:**
- [x] IC50 / %Rank boundary values classified by strict `<` (50/500 nM; 0.5/2/10 %Rank)
- [x] Peptide length outside class range rejected (class I 8–11; class II 13–25)
- [x] Invalid input rejected (IC50 ≤ 0 / non-finite; %Rank ∉ [0,100] / NaN)
- [ ] Affinity / %Rank prediction (out of scope; trained model is caller-supplied)

---

#### ONCO-IMMUNE-001: Immune Infiltration Estimation

| Field | Value |
|------|----------|
| **Canonical** | `ImmuneAnalyzer.EstimateInfiltration(...)` |
| **Complexity** | O(n × k) |
| **Depends on** | ONCO-EXPR-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `EstimateInfiltration(expressionProfile)` | ImmuneAnalyzer | Canonical |
| `DeconvoluteImmuneCells(expression)` | ImmuneAnalyzer | Cell type fractions |

---

#### ONCO-CTDNA-001: ctDNA Analysis

| Field | Value |
|------|----------|
| **Canonical** | `LiquidBiopsyAnalyzer.AnalyzeCtDNA(...)` |
| **Complexity** | O(n) |
| **Depends on** | ONCO-SOMATIC-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeCtDNA(plasmaVcf, tumorVcf)` | LiquidBiopsyAnalyzer | Canonical |
| `CalculateTumorFraction(variants)` | LiquidBiopsyAnalyzer | ctDNA fraction |
| `AnalyzeFragmentSizeDistribution(bamFile)` | LiquidBiopsyAnalyzer | Fragment analysis |

**Edge Cases:**
- [x] ctDNA fraction < 0.1% (below detection limit) — covered by Poisson LoD (λ < 1 ⇒ not detected); tests M6/S1–S5.
- [ ] High background of CHIP mutations — out of scope; handled by ONCO-CHIP-001 (see CtDNA_Analysis.md §5.3).
- [ ] No matched primary tumor VCF — out of scope; matched-tumour cross-referencing belongs to ONCO-MRD-001.

---

#### ONCO-MRD-001: Minimal Residual Disease Detection ☑

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.DetectMRD(...)` |
| **Complexity** | O(n × k) where k=tracked mutations |
| **Depends on** | ONCO-CTDNA-001 |

> Class note: the area's liquid-biopsy/ctDNA methods live on `OncologyAnalyzer` (no separate
> `LiquidBiopsyAnalyzer` exists), mirroring ONCO-CTDNA-001. Behaviour: tumour-informed panel-level call —
> MRD-positive when ≥2 tracked patient-specific variants are detected (Reinert 2019 / Signatera; PMC9265001
> Table 1). Reuses the ONCO-CTDNA-001 Poisson primitive `p = 1 − e^(−n·f·m)` and reports the INVAR
> integrated mutant allele fraction (Wan 2020).

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectMRD(tumorMarkers, positivityThreshold, minSupportingReads, genomeEquivalents)` | OncologyAnalyzer | Canonical |
| `TrackVariantsOverTime(timepoints, ...)` | OncologyAnalyzer | Longitudinal |

**Edge Cases:**
- [x] Exactly 1 variant detected — below the ≥2 rule ⇒ MRD-negative (tests M2; S2).
- [x] Empty / null marker panel — ArgumentException / ArgumentNullException (tests C1, C2).
- [x] Configurable positivity threshold and per-locus supporting-read cutoff (tests S1, S3).

---

#### ONCO-CHIP-001: Clonal Hematopoiesis Filtering ☑

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.FilterCHIP(...)` |
| **Complexity** | O(n) |
| **Depends on** | ONCO-CTDNA-001 |

> Conflict resolved: by-area table originally named `LiquidBiopsyAnalyzer`, but all ONCO units
> live in `OncologyAnalyzer` (the actual repository class); implemented there to match every sibling.

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FilterCHIP(variants, whiteBloodCellVariants, ...)` | OncologyAnalyzer | Canonical (matched-WBC subtraction, Razavi 2019) |
| `IdentifyCHIPVariants(variants, chipGenes?, minVaf?)` | OncologyAnalyzer | Canonical (gene + VAF≥0.02, Steensma 2015) |
| `IsCanonicalChipGene(gene, chipGenes?)` | OncologyAnalyzer | Internal (case-insensitive panel membership) |

**CHIP Genes (canonical default, caller-overridable):** DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D — Steensma 2015 Fig 2A / Genovese 2014.
**VAF threshold:** ≥ 2% (0.02), inclusive — Steensma et al. (2015) *Blood* 126(1):9–16.

**Edge Cases:**
- [x] CHIP gene at VAF ≥ 0.02 flagged; sub-2% not flagged (Steensma 2015).
- [x] Non-CHIP gene never flagged regardless of VAF.
- [x] Matched-WBC variant removed regardless of gene (Razavi 2019).
- [x] WBC-absent variant retained as tumour candidate.
- [x] null / empty / out-of-range inputs validated.

---

#### ONCO-PHYLO-001: Tumor Phylogeny Reconstruction

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.ReconstructPhylogeny(...)` (placed in OncologyAnalyzer, the area's existing analyzer) |
| **Complexity** | O(n² × k) |
| **Depends on** | ONCO-CCF-001 |
| **Status** | ☑ Complete — [Evidence](docs/Evidence/ONCO-PHYLO-001-Evidence.md) · [TestSpec](tests/TestSpecs/ONCO-PHYLO-001.md) · [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ReconstructPhylogeny_Tests.cs) |

**Method:** clonal tree built from CCF clusters via the lineage-precedence rule (ancestor CCF ≥ descendant CCF, Popic et al. 2015 *Genome Biology* 16:91 Eq. 2; Zheng et al. 2022 *Bioinformatics* PICTograph) and the sum rule (children CCF sum ≤ parent CCF, Popic 2015 Eq. 5). CCF clustering itself is ONCO-CCF-001.

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ReconstructPhylogeny(IReadOnlyList<CcfCluster>, tolerance)` | OncologyAnalyzer | Canonical |
| `IdentifyTrunkMutations(phylogeny)` | OncologyAnalyzer | Clonal mutations |
| `IdentifyBranchMutations(phylogeny)` | OncologyAnalyzer | Subclonal mutations |

**Edge Cases:**
- [x] Empty cluster set → root-only tree
- [x] Single cluster → trunk only, no branches
- [x] Sum-rule conflict (sibling CCFs exceed parent) → forced chain
- [x] Private (single-sample) clusters → branching via presence pattern
- [x] Null / NaN / out-of-[0,1] CCF / ragged sample counts / duplicate ids → exceptions

---

#### ONCO-CCF-001: Cancer Cell Fraction Estimation ☑

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.EstimateCcf(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ CCF ≤ 1 |
| **Depends on** | ONCO-VAF-001, ONCO-PURITY-001, ONCO-CNA-001 |
| **Status** | ☑ Complete — [Evidence](docs/Evidence/ONCO-CCF-001-Evidence.md) · [TestSpec](tests/TestSpecs/ONCO-CCF-001.md) · [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_EstimateCcf_Tests.cs) |

**Method:** point CCF estimate `CCF = VAF·(ρ·N_T + 2(1−ρ))/(ρ·m)` (McGranahan et al. 2016 *Science* 351:1463–1469; Tarabichi et al. 2021 *Nat. Methods* 18:144–155 Box 1; Zheng et al. 2022 *Bioinformatics* 38:3677), reported capped to [0,1] with raw exposed; CCF values clustered into clones/subclones by deterministic 1D Lloyd k-means (Lloyd 1982 *IEEE TIT* 28:129–137), clonal cluster = highest centroid (Tarabichi 2021). Distinct from ONCO-CLONAL-001 (Bayesian grid posterior). Implemented on `OncologyAnalyzer` (the area's analyzer class used by sibling units), not the registry's placeholder `TumorEvolutionAnalyzer`.

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `EstimateCcf(vaf, purity, tumorCopyNumber, multiplicity)` | OncologyAnalyzer | Canonical |
| `ClusterCcfValues(ccfValues, clusterCount)` | OncologyAnalyzer | Subclone inference |

**Edge Cases:**
- [x] Purity not determined (unknown purity) — purity ∉ (0,1] rejected (ArgumentOutOfRangeException)
- [x] Multi-copy loci (ambiguous CCF) — multiplicity m ∈ [1, N_T] required; tested at N_T=4, m=2

---

#### ONCO-HETERO-001: Tumor Heterogeneity Analysis ☑

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.AnalyzeHeterogeneity(...)` |
| **Complexity** | O(n log n) |
| **Invariant** | ITH_score ≥ 0 |
| **Depends on** | ONCO-CCF-001 |
| **Status** | ☑ Complete — [Evidence](docs/Evidence/ONCO-HETERO-001-Evidence.md) · [TestSpec](tests/TestSpecs/ONCO-HETERO-001.md) · [Tests](tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AnalyzeHeterogeneity_Tests.cs) |

**Note:** implemented in `OncologyAnalyzer` (the area's analyzer class; the registry's `TumorEvolutionAnalyzer` label is logical, mirroring ONCO-PHYLO-001/ONCO-CCF-001). ITH metrics: MATH = 100·1.4826·MAD(VAF)/median(VAF) (Mroz & Rocco 2013); Shannon H = −Σ pᵢ ln pᵢ over CCF-cluster fractions (Shannon 1948; Liu & Zhang 2017); subclone count = occupied CCF clusters; subclonal fraction = #(CCF < 0.95)/n (Landau 2013).

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeHeterogeneity(variants, ccfValues, clusterCount)` | OncologyAnalyzer | Canonical |
| `CalculateITH(ccfDistribution)` | OncologyAnalyzer | ITH score (MATH) |
| `InferSubclones(ccfClusters)` | OncologyAnalyzer | Subclone count |

---

#### ONCO-HLA-001: HLA nomenclature parsing + allele-specific HLA LOH (LOHHLA)

> **Scope note (external evidence wins):** the original by-area definition listed a full HLA genotyping
> caller (`HlaTyper.TypeHLA(bamFile)`). That caller requires a trained model / IPD-IMGT/HLA reference
> allele database whose behaviour is not retrievable as an exact specification; fabricating it would
> violate the evidence-first policy. This unit instead implements the two retrievable, formally-specified
> pieces and takes genotype/coverage as caller-supplied: (a) HLA allele nomenclature parsing/validation
> per the WHO HLA Nomenclature (Marsh et al. 2010; hla.alleles.org), and (b) allele-specific HLA LOH
> classification per LOHHLA (McGranahan et al. 2017, Cell).

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.DetectHlaLoh(...)`, `OncologyAnalyzer.ParseHlaAllele(...)` |
| **Complexity** | parse O(n); LOH classification O(1) |
| **Invariant** | HLA LOH ⇔ one allele CN < 0.5 ∧ allelic-imbalance p < 0.01 (McGranahan 2017) |
| **Depends on** | — (standalone) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ParseHlaAllele(name)` | OncologyAnalyzer | Canonical |
| `TryParseHlaAllele(name, out allele)` | OncologyAnalyzer | Delegate |
| `DetectHlaLoh(alleleCopyNumber)` | OncologyAnalyzer | Canonical |

**Databases:** IPD-IMGT/HLA (nomenclature only; not used as a runtime allele DB)

**Edge Cases:**
- [x] Boundary copy number (CN exactly 0.5 → retained)
- [x] Over-calling guard (low CN but allelic-imbalance p ≥ 0.01 → no LOH)
- [x] Homozygous loss (both alleles CN < 0.5 → not allele-specific LOH)
- [x] Malformed nomenclature (missing prefix / wrong field count / invalid suffix)

---

#### ONCO-ACTION-001: Clinical Actionability Assessment ☑

> Scoped to OncoKB therapeutic levels of evidence (Chakravarty 2017); AMP/ASCO/CAP tiering is covered by
> ONCO-ANNOT-001. Implemented in `OncologyAnalyzer` (not a separate `ClinicalInterpreter` class).

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.AssessActionability(...)` |
| **Complexity** | O(n × k) where k=knowledge base entries |
| **Invariant** | evidence tier ∈ {Tier I, Tier II, Tier III, Tier IV} |
| **Depends on** | ONCO-DRIVER-001, ONCO-ANNOT-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AssessActionability(variants)` | ClinicalInterpreter | Canonical |
| `GetTherapyRecommendations(variant)` | ClinicalInterpreter | Therapy lookup |
| `ClassifyEvidenceLevel(variant)` | ClinicalInterpreter | AMP/ASCO/CAP tiers |

**Databases:** OncoKB, CIViC, ClinVar

**Edge Cases:**
- [x] VUS (variants of uncertain significance) — no leveled association ⇒ NotActionable
- [x] Conflicting evidence across databases — highest level wins (max under combined order)
- [x] Off-label therapy recommendations — Level 3B (response in another indication)

---

#### ONCO-SV-001: Non-Fusion Structural Variant Detection

| Field | Value |
|------|----------|
| **Canonical** | `StructuralVariantDetector.DetectStructuralVariants(...)` |
| **Complexity** | O(n log n) |
| **Invariant** | SV type ∈ {DEL, DUP, INV, TRA, COMPLEX} |
| **Depends on** | — (standalone) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectStructuralVariants(bamFile)` | StructuralVariantDetector | Canonical |
| `ClassifySVType(sv)` | StructuralVariantDetector | Type classification |
| `FilterSomaticSVs(svs, normalBam)` | StructuralVariantDetector | Somatic filter |

**Edge Cases:**
- [x] Complex rearrangements (chromothripsis) — `ClassifyComplexRearrangement` (Korbel & Campbell 2013 hallmark criteria; Cortés-Ciriano 2020 thresholds)
- [ ] Low-coverage samples (insufficient split reads) — N/A: generic read-based SV detection is covered by `StructuralVariantAnalyzer` (out of scope for the oncology complex-rearrangement layer)
- [ ] Repetitive regions (ambiguous breakpoints) — N/A: out of scope (read-level breakpoint refinement is `StructuralVariantAnalyzer`)

> **Scope note (ONCO-SV-001):** Generic SV typing (DEL/DUP/INV/TRA from breakpoints) already exists in `StructuralVariantAnalyzer`. ONCO-SV-001 implements the oncology-specific complex-rearrangement layer — chromothripsis inference (`OncologyAnalyzer.ClassifyComplexRearrangement`, `CountCopyNumberStateOscillations`, `TestBreakpointClustering`) per Korbel & Campbell (2013) and Cortés-Ciriano et al. (2020). The two low-coverage / repetitive-region edge cases pertain to read-level detection and remain out of scope.

---

#### ONCO-EXPR-001: Tumor Gene Expression Quantification

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.IdentifyOutlierGenes(...)` (z-score outliers) + `CalculateSignatureScore(...)` |
| **Complexity** | O(g × n), g=gene count, n=reference cohort size |
| **Invariant** | z=(r−μ)/σ; outlier iff z>+t or z<−t (strict, t=2); signature a=Σz/√k |
| **Depends on** | — (standalone) |

**Methods (as implemented — oncology-specific expression layer on `OncologyAnalyzer`):**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateExpressionZScore(value, referenceCohort)` | OncologyAnalyzer | Canonical (z = (r−μ)/σ, sample SD n−1; cBioPortal) |
| `IdentifyOutlierGenes(sampleExpression, referenceCohorts, threshold)` | OncologyAnalyzer | Canonical (strict z>±2 over/under outliers) |
| `CalculateSignatureScore(memberZScores)` | OncologyAnalyzer | Canonical (combined z-score a=Σz/√k; Lee et al. 2008) |

Note: TPM/FPKM quantification (`CalculateTPM`/`CalculateFPKM`, TRANS-EXPR-001) and differential expression
(TRANS-DIFF-001) already exist in `TranscriptomeAnalyzer`; ONCO-EXPR-001 implements the oncology-specific
outlier/signature layer over caller-supplied reference cohorts and signatures.

**Edge Cases:**
- [x] Zero-variance / degenerate reference cohort (sd=0 → throws; mirrors NormalizeExpressionLevels.java)
- [x] No reference cohort for a queried gene / reference smaller than 2 samples (sample SD undefined → throws)
- [x] Batch effects / unrepresentative cohort (documented limitation: z-scores miscalibrated; caller responsibility — algorithm doc §6.2)

---

## Appendix A: Method Index

| Method | Test Unit ID |
|--------|--------------|
| `AssembleDeBruijn` | ASSEMBLY-DBG-001 |
| `AssembleOLC` | ASSEMBLY-OLC-001 |
| `CalculateAlphaDiversity` | META-ALPHA-001 |
| `CalculateAlleleFrequencies` | POP-FREQ-001 |
| `CalculateBetaDiversity` | META-BETA-001 |
| `CalculateCAI` | CODON-CAI-001 |
| `CalculateCodonUsage` | CODON-USAGE-001 |
| `CalculateCumulativeGcSkew` | SEQ-GCSKEW-001 |
| `CalculateDistanceMatrix` | PHYLO-DIST-001 |
| `CalculateDiversityStatistics` | POP-DIV-001 |
| `CalculateFoldChange` | TRANS-DIFF-001 |
| `CalculateFPKM` | TRANS-EXPR-001 |
| `CalculateFreeEnergy` | RNA-ENERGY-001 |
| `CalculateFst` | POP-FST-001 |
| `CalculateGcContent` | SEQ-GC-001 |
| `CalculateGcFraction` | SEQ-GC-001 |
| `CalculateGcSkew` | SEQ-GCSKEW-001 |
| `CalculateKmerEntropy` | KMER-FREQ-001 |
| `CalculateLD` | POP-LD-001 |
| `CalculateLinguisticComplexity` | SEQ-COMPLEX-001 |
| `CalculateMeltingTemperature` | PRIMER-TM-001 |
| `CalculateMethylationProfile` | EPIGEN-METHYL-001 |
| `CalculateN50` | ASSEMBLY-STATS-001 |
| `CalculateNucleotideDiversity` | POP-DIV-001 |
| `CalculateObservedExpectedCpG` | EPIGEN-CPG-001 |
| `CalculatePSI` | TRANS-SPLICE-001 |
| `CalculateQ30Percentage` | QUALITY-STATS-001 |
| `CalculateShannonEntropy` | SEQ-ENTROPY-001 |
| `CalculateSpecificityScore` | CRISPR-OFF-001 |
| `CalculateStatistics` | ASSEMBLY-STATS-001, QUALITY-STATS-001 |
| `CalculateTajimasD` | POP-DIV-001 |
| `CalculateTPM` | TRANS-EXPR-001 |
| `CalculateWattersonTheta` | POP-DIV-001 |
| `CalculateWindowedGcSkew` | SEQ-GCSKEW-001 |
| `CallVariants` | VARIANT-CALL-001 |
| `CallVariantsFromAlignment` | VARIANT-CALL-001 |
| `CheckSpecificity` | PROBE-VALID-001 |
| `ClassifyReads` | META-CLASS-001 |
| `ClusterGenes` | PANGEN-CLUSTER-001 |
| `Complement` | SEQ-COMP-001 |
| `ConstructPanGenome` | PANGEN-CORE-001 |
| `Contains` | PAT-EXACT-001 |
| `ConvertEncoding` | QUALITY-PHRED-001 |
| `CountKmers` | KMER-COUNT-001 |
| `CreateMiRna` | MIRNA-SEED-001 |
| `CreatePwm` | PAT-PWM-001 |
| `DesignGuideRnas` | CRISPR-GUIDE-001 |
| `DesignPrimers` | PRIMER-DESIGN-001 |
| `DesignProbes` | PROBE-DESIGN-001 |
| `DesignTilingProbes` | PROBE-DESIGN-001 |
| `DetectAlternativeSplicing` | TRANS-SPLICE-001 |
| `DetectAneuploidy` | CHROM-ANEU-001 |
| `DetectCNV` | SV-CNV-001 |
| `DetectPloidy` | CHROM-KARYO-001 |
| `DetectRearrangements` | COMPGEN-REARR-001 |
| `DetectSVs` | SV-DETECT-001 |
| `Digest` | RESTR-DIGEST-001 |
| `EditDistance` | PAT-APPROX-002 |
| `EvaluateGuideRna` | CRISPR-GUIDE-001 |
| `EvaluatePrimer` | PRIMER-DESIGN-001 |
| `FilterByQuality` | PARSE-FASTQ-001 |
| `FindAcceptorSites` | SPLICE-ACCEPTOR-001 |
| `FindAllOccurrences` | PAT-EXACT-001 |
| `FindAllOverlaps` | ASSEMBLY-OLC-001 |
| `FindBreakpoints` | SV-BREAKPOINT-001 |
| `FindClumps` | KMER-FIND-001 |
| `FindCpGIslands` | EPIGEN-CPG-001 |
| `FindCpGSites` | EPIGEN-CPG-001 |
| `FindDegenerateMotif` | PAT-IUPAC-001 |
| `FindDifferentiallyExpressed` | TRANS-DIFF-001 |
| `FindDirectRepeats` | REP-DIRECT-001 |
| `FindDiscordantPairs` | SV-DETECT-001 |
| `FindDMRs` | EPIGEN-DMR-001 |
| `FindDonorSites` | SPLICE-DONOR-001 |
| `FindExactMotif` | PAT-EXACT-001 |
| `FindExons` | SPLICE-PREDICT-001 |
| `FindHaplotypeBlocks` | POP-LD-001 |
| `FindIntrons` | SPLICE-PREDICT-001 |
| `FindInvertedRepeats` | REP-INV-001 |
| `FindLowComplexityRegions` | SEQ-COMPLEX-001 |
| `FindMethylationSites` | EPIGEN-METHYL-001 |
| `FindMicrosatellites` | REP-STR-001 |
| `FindMostFrequentKmers` | KMER-FIND-001 |
| `FindMotif` | PAT-EXACT-001 |
| `FindMotifs` | PROTMOTIF-FIND-001 |
| `FindOffTargets` | CRISPR-OFF-001 |
| `FindOrfs` | ANNOT-ORF-001 |
| `PredictReplicationOrigin` | SEQ-GCSKEW-001 |
| `FindOrthologs` | COMPGEN-ORTHO-001 |
| `FindPalindromes` | REP-PALIN-001 |
| `FindPamSites` | CRISPR-PAM-001 |
| `FindParalogs` | COMPGEN-ORTHO-001 |
| `FindPreMiRnaHairpins` | MIRNA-PRECURSOR-001 |
| `FindPromoterMotifs` | ANNOT-PROM-001 |
| `FindRareCodons` | CODON-RARE-001 |
| `FindRibosomeBindingSites` | ANNOT-GENE-001 |
| `FindSites` | RESTR-FIND-001 |
| `FindStemLoops` | RNA-STEMLOOP-001 |
| `FindSyntenicBlocks` | COMPGEN-SYNTENY-001 |
| `FindSyntenyBlocks` | CHROM-SYNT-001 |
| `FindTandemRepeats` | REP-TANDEM-001 |
| `FindTargetSites` | MIRNA-TARGET-001 |
| `FindUniqueKmers` | KMER-FIND-001 |
| `FindWithEdits` | PAT-APPROX-002 |
| `FindWithMismatches` | PAT-APPROX-001 |
| `GeneratePresenceAbsenceMatrix` | PANGEN-CLUSTER-001 |
| `GenerateTaxonomicProfile` | META-PROF-001 |
| `GetComplementBase` | SEQ-COMP-001 |
| `GetKmerFrequencies` | KMER-FREQ-001 |
| `GetKmerSpectrum` | KMER-FREQ-001 |
| `GetSeedSequence` | MIRNA-SEED-001 |
| `GlobalAlign` | ALIGN-GLOBAL-001 |
| `HammingDistance` | PAT-APPROX-001 |
| `HasHairpinPotential` | PRIMER-STRUCT-001 |
| `HasPrimerDimer` | PRIMER-STRUCT-001 |
| `GetCoreGeneClusters` (a.k.a. `IdentifyCoreGenes`) | PANGEN-CORE-001 |
| `IdentifyDisorderedRegions` | DISORDER-REGION-001 |
| `IsStartCodon` | TRANS-CODON-001 |
| `IsStopCodon` | TRANS-CODON-001 |
| `IsValidDna` | SEQ-VALID-001 |
| `IsValidRna` | SEQ-VALID-001 |
| `LocalAlign` | ALIGN-LOCAL-001 |
| `MaskLowComplexity` | SEQ-COMPLEX-001 |
| `MatchesIupac` | PAT-IUPAC-001 |
| `MatchPrositePattern` | PROTMOTIF-PROSITE-001 |
| `MergeOverlapping` | PARSE-BED-001 |
| `MultipleAlign` | ALIGN-MULTI-001 |
| `OptimizeSequence` | CODON-OPT-001 |
| `Parse` | PARSE-FASTA-001, PARSE-FASTQ-001, PARSE-BED-001, PARSE-VCF-001 |
| `ParseFile` | PARSE-FASTA-001, PARSE-FASTQ-001, PARSE-BED-001, PARSE-VCF-001 |
| `ParseGff3` | ANNOT-GFF-001 |
| `ParseNewick` | PHYLO-NEWICK-001 |
| `ParsePrositePattern` | PROTMOTIF-PROSITE-001 |
| `ParseQualityString` | QUALITY-PHRED-001 |
| `PatristicDistance` | PHYLO-COMP-001 |
| `Predict` | RNA-STRUCT-001 |
| `PredictDisorder` | DISORDER-PRED-001 |
| `PredictDomains` | PROTMOTIF-DOMAIN-001 |
| `PredictGenes` | ANNOT-GENE-001 |
| `PredictGeneStructure` | SPLICE-PREDICT-001 |
| `PredictSignalPeptide` | PROTMOTIF-DOMAIN-001 |
| `QuantileNormalize` | TRANS-EXPR-001 |
| `ReverseComplement` | SEQ-REVCOMP-001 |
| `RobinsonFouldsDistance` | PHYLO-COMP-001 |
| `ScanForPattern` | PROTMOTIF-FIND-001 |
| `ScanWithPwm` | PAT-PWM-001 |
| `ScoreDonorSite` | SPLICE-DONOR-001 |
| `ScoreAcceptorSite` | SPLICE-ACCEPTOR-001 |
| `ScoreProbe` | PROBE-DESIGN-001 |
| `ScoreTargetSite` | MIRNA-TARGET-001 |
| `SegmentCopyNumber` | SV-CNV-001 |
| `SemiGlobalAlign` | ALIGN-SEMI-001 |
| `TestHardyWeinberg` | POP-HW-001 |
| `ToDotBracket` | RNA-STRUCT-001 |
| `ToFasta` | PARSE-FASTA-001 |
| `ToGff3` | ANNOT-GFF-001 |
| `ToNewick` | PHYLO-NEWICK-001 |
| `ToQualityString` | QUALITY-PHRED-001 |
| `Translate` | TRANS-CODON-001, TRANS-PROT-001 |
| `TryGetComplement` | SEQ-COMP-001 |
| `TryGetReverseComplement` | SEQ-REVCOMP-001 |
| `ValidateProbe` | PROBE-VALID-001 |
| `CalculateNucleotideComposition` | SEQ-STATS-001 |
| `CalculateAminoAcidComposition` | SEQ-STATS-001 |
| `SummarizeNucleotideSequence` | SEQ-STATS-001 |
| `CalculateMolecularWeight` | SEQ-MW-001 |
| `CalculateNucleotideMolecularWeight` | SEQ-MW-001 |
| `CalculateIsoelectricPoint` | SEQ-PI-001 |
| `CalculateHydrophobicity` | SEQ-HYDRO-001 |
| `CalculateHydrophobicityProfile` | SEQ-HYDRO-001 |
| `CalculateThermodynamics` | SEQ-THERMO-001 |
| `CalculateDinucleotideFrequencies` | SEQ-DINUC-001 |
| `CalculateDinucleotideRatios` | SEQ-DINUC-001 |
| `CalculateCodonFrequencies` | SEQ-DINUC-001 |
| `PredictSecondaryStructure` | SEQ-SECSTRUCT-001 |
| `CalculateRscu` | CODON-RSCU-001 |
| `CountCodons` | CODON-RSCU-001 |
| `CalculateEnc` | CODON-ENC-001 |
| `GetStatistics` (CodonUsageAnalyzer) | CODON-STATS-001 |
| `CalculateCai` | CODON-STATS-001 |
| `EColiOptimalCodons` | CODON-STATS-001 |
| `HumanOptimalCodons` | CODON-STATS-001 |
| `TranslateSixFrames` | TRANS-SIXFRAME-001 |
| `FindOrfs` (Translator) | TRANS-SIXFRAME-001 |
| `MergeContigs` | ASSEMBLY-MERGE-001 |
| `Scaffold` | ASSEMBLY-SCAFFOLD-001 |
| `CalculateCoverage` | ASSEMBLY-COVER-001 |
| `ComputeConsensus` | ASSEMBLY-CONSENSUS-001 |
| `QualityTrimReads` | ASSEMBLY-TRIM-001 |
| `ErrorCorrectReads` | ASSEMBLY-CORRECT-001 |
| `FindBestMatch` | PAT-APPROX-003 |
| `CountApproximateOccurrences` | PAT-APPROX-003 |
| `FindFrequentKmersWithMismatches` | PAT-APPROX-003 |
| `CalculateStatistics` (SequenceAligner) | ALIGN-STATS-001 |
| `FormatAlignment` | ALIGN-STATS-001 |
| `SimulateBisulfiteConversion` | EPIGEN-BISULF-001 |
| `CalculateMethylationFromBisulfite` | EPIGEN-BISULF-001 |
| `GenerateMethylationProfile` | EPIGEN-BISULF-001 |
| `PredictChromatinState` | EPIGEN-CHROM-001 |
| `AnnotateHistoneModifications` | EPIGEN-CHROM-001 |
| `FindAccessibleRegions` | EPIGEN-CHROM-001 |
| `CalculateEpigeneticAge` | EPIGEN-AGE-001 |
| `PredictImprintedGenes` | EPIGEN-AGE-001 |
| `AlignMiRnaToTarget` | MIRNA-PAIR-001 |
| `GetReverseComplement` (MiRnaAnalyzer) | MIRNA-PAIR-001 |
| `CanPair` | MIRNA-PAIR-001 |
| `IsWobblePair` | MIRNA-PAIR-001 |
| `FitHeapsLaw` | PANGEN-HEAP-001 |
| `CreatePresenceAbsenceMatrix` | PANGEN-HEAP-001 |
| `SelectPhylogeneticMarkers` | PANGEN-MARKER-001 |
| `GetCoreGeneClusters` | PANGEN-MARKER-001 |
| `CreateCoreGenomeAlignment` | PANGEN-MARKER-001 |
| `GetSingletonGenes` | PANGEN-MARKER-001 |
| `CalculateIHS` | POP-SELECT-001 |
| `ScanForSelection` | POP-SELECT-001 |
| `EstimateAncestry` | POP-ANCESTRY-001 |
| `FindROH` | POP-ROH-001 |
| `CalculateInbreedingFromROH` | POP-ROH-001 |
| `CalculatePairwiseFst` | POP-ROH-001 |
| `PredictFunctions` | META-FUNC-001 |
| `FindAntibioticResistanceGenes` | META-RESIST-001 |
| `FindPathwayEnrichment` | META-PATHWAY-001 |
| `FindSignificantTaxa` | META-TAXA-001 |
| `MannWhitneyU` | META-TAXA-001 |
| `Bootstrap` | PHYLO-BOOT-001 |
| `GetLeaves` | PHYLO-STATS-001 |
| `CalculateTreeLength` | PHYLO-STATS-001 |
| `GetTreeDepth` | PHYLO-STATS-001 |
| `CalculateCodingPotential` | ANNOT-CODING-001 |
| `FindRepetitiveElements` | ANNOT-REPEAT-001 |
| `GetCodonUsage` | ANNOT-CODONUSAGE-001 |
| `GetEnzymesByCutLength` | RESTR-FILTER-001 |
| `GetBluntCutters` | RESTR-FILTER-001 |
| `GetStickyCutters` | RESTR-FILTER-001 |
| `KmerDistance` | KMER-DIST-001 |
| `CreateConsensusFromAlignment` | MOTIF-CONS-001 |
| `FindLongestRepeat` | GENOMIC-REPEAT-001 |
| `FindRepeats` | GENOMIC-REPEAT-001 |
| `FindLongestCommonRegion` | GENOMIC-COMMON-001 |
| `FindCommonRegions` | GENOMIC-COMMON-001 |
| `FindKnownMotifs` | GENOMIC-MOTIFS-001 |
| `GetRnaComplementBase` | SEQ-RNACOMP-001 |
| `FindMotifByPattern` | PROTMOTIF-PATTERN-001 |
| `FindMotifByProsite` | PROTMOTIF-PATTERN-001 |
| `ConvertPrositeToRegex` | PROTMOTIF-PATTERN-001 |
| `FindDomains` | PROTMOTIF-PATTERN-001 |
| `PredictTransmembraneHelices` | PROTMOTIF-TM-001 |
| `PredictCoiledCoils` | PROTMOTIF-CC-001 |
| `FindLowComplexityRegions` (ProteinMotifFinder) | PROTMOTIF-LC-001 |
| `FindCommonMotifs` | PROTMOTIF-COMMON-001 |
| `CanPair` (RnaSecondaryStructure) | RNA-PAIR-001 |
| `GetBasePairType` | RNA-PAIR-001 |
| `GetComplement` (RnaSecondaryStructure) | RNA-PAIR-001 |
| `CalculateStemEnergy` | RNA-HAIRPIN-001 |
| `CalculateHairpinLoopEnergy` | RNA-HAIRPIN-001 |
| `CalculateMinimumFreeEnergy` | RNA-MFE-001 |
| `PredictStructure` (RnaSecondaryStructure) | RNA-MFE-001 |
| `DetectPseudoknots` | RNA-PSEUDOKNOT-001 |
| `ParseDotBracket` | RNA-DOTBRACKET-001 |
| `ValidateDotBracket` | RNA-DOTBRACKET-001 |
| `FindInvertedRepeats` (RnaSecondaryStructure) | RNA-INVERT-001 |
| `CalculatePartitionFunction` | RNA-PARTITION-001 |
| `CalculateStructureProbability` | RNA-PARTITION-001 |
| `GenerateRandomRna` | RNA-PARTITION-001 |
| `CalculateKmerEntropy` (SequenceComplexity) | SEQ-COMPLEX-KMER-001 |
| `CalculateWindowedComplexity` | SEQ-COMPLEX-WINDOW-001 |
| `CalculateDustScore` | SEQ-COMPLEX-DUST-001 |
| `MaskLowComplexity` (SequenceComplexity) | SEQ-COMPLEX-DUST-001 |
| `EstimateCompressionRatio` | SEQ-COMPLEX-COMPRESS-001 |
| `FindReciprocalBestHits` | COMPGEN-RBH-001 |
| `CompareGenomes` | COMPGEN-COMPARE-001 |
| `CalculateReversalDistance` | COMPGEN-REVERSAL-001 |
| `FindConservedClusters` | COMPGEN-CLUSTER-001 |
| `CalculateANI` | COMPGEN-ANI-001 |
| `GenerateDotPlot` | COMPGEN-DOTPLOT-001 |
| `DiscoverMotifs` | MOTIF-DISCOVER-001 |
| `FindSharedMotifs` | MOTIF-SHARED-001 |
| `FindRegulatoryElements` | MOTIF-REGULATORY-001 |
| `GenerateConsensus` | MOTIF-GENERATE-001 |
| `CountKmersAsync` | KMER-ASYNC-001 |
| `CountKmersSpan` | KMER-ASYNC-001 |
| `FindUniqueKmers` | KMER-UNIQUE-001 |
| `FindKmersWithMinCount` | KMER-UNIQUE-001 |
| `GenerateAllKmers` | KMER-GENERATE-001 |
| `CountKmersBothStrands` | KMER-BOTH-001 |
| `AnalyzeKmers` | KMER-STATS-001 |
| `FindKmerPositions` | KMER-POSITIONS-001 |
| `CalculateAtSkew` | SEQ-ATSKEW-001 |
| `PredictReplicationOrigin` | SEQ-REPLICATION-001 |
| `AnalyzeGcContent` | SEQ-GC-ANALYSIS-001 |
| `PredictMoRFs` | DISORDER-MORF-001 |
| `GetDisorderPropensity` | DISORDER-PROPENSITY-001 |
| `IsDisorderPromoting` | DISORDER-PROPENSITY-001 |
| `DisorderPromotingAminoAcids` | DISORDER-PROPENSITY-001 |
| `OrderPromotingAminoAcids` | DISORDER-PROPENSITY-001 |
| `PredictLowComplexityRegions` (DisorderPredictor) | DISORDER-LC-001 |
| `CalculateNucleotideComposition` | SEQ-COMPOSITION-001 |
| `CalculateAminoAcidComposition` | SEQ-COMPOSITION-001 |
| `CalculateMeltingTemperature` (SequenceStatistics) | SEQ-TM-001 |
| `CalculateEntropyProfile` | SEQ-ENTROPY-PROFILE-001 |
| `CalculateGcContentProfile` | SEQ-GC-PROFILE-001 |
| `CalculateCodonFrequencies` (SequenceStatistics) | SEQ-CODON-FREQ-001 |
| `SummarizeNucleotideSequence` | SEQ-SUMMARY-001 |
| `FindTandemRepeats` (GenomicAnalyzer) | GENOMIC-TANDEM-001 |
| `CalculateSimilarity` | GENOMIC-SIMILARITY-001 |
| `FindOpenReadingFrames` | GENOMIC-ORF-001 |
| `AdjustVAFForPurity` | ONCO-VAF-001 |
| `AnalyzeBreakpoint` | ONCO-FUSION-003 |
| `AnalyzeCtDNA` | ONCO-CTDNA-001 |
| `AnalyzeFragmentSizeDistribution` | ONCO-CTDNA-001 |
| `AnalyzeHeterogeneity` (TumorEvolutionAnalyzer) | ONCO-HETERO-001 |
| `AnnotateCancerVariants` | ONCO-ANNOT-001 |
| `AssessActionability` (OncologyAnalyzer) | ONCO-ACTION-001 |
| `BootstrapConfidenceIntervals` (MutationalSignatures) | ONCO-SIG-003 |
| `Build96ChannelSpectrum` | ONCO-SIG-002 |
| `CalculateBindingAffinity` | ONCO-MHC-001 |
| `CalculateExpressionZScore` | ONCO-EXPR-001 |
| `CalculateHRDScore` | ONCO-HRD-001 |
| `CalculateITH` | ONCO-HETERO-001 |
| `CalculateLOHFraction` | ONCO-LOH-001 |
| `CalculateLOHScore` | ONCO-HRD-001 |
| `CalculateLSTScore` | ONCO-HRD-001 |
| `CalculateMSIScore` | ONCO-MSI-001 |
| `CalculateSignatureScore` | ONCO-EXPR-001 |
| `CalculateSomaticScore` | ONCO-SOMATIC-001 |
| `CalculateTAIScore` | ONCO-HRD-001 |
| `CalculateTMB` | ONCO-TMB-001 |
| `CalculateTumorFraction` | ONCO-CTDNA-001 |
| `CalculateVAF` | ONCO-VAF-001 |
| `CallCopyNumberStates` | ONCO-CNA-001 |
| `CallSomaticMutations` | ONCO-SOMATIC-001 |
| `ClassifyClonality` | ONCO-CLONAL-001 |
| `ClassifyEvidenceLevel` | ONCO-ACTION-001 |
| `ClassifyMSIStatus` | ONCO-MSI-001 |
| `ClassifyMutationalProcess` | ONCO-SIG-004 |
| `ClassifySVType` | ONCO-SV-001 |
| `ClassifyTMB` | ONCO-TMB-001 |
| `ClusterCCFValues` | ONCO-CCF-001 |
| `DecomposeSpectrum` | ONCO-SIG-001 |
| `DeconvoluteImmuneCells` | ONCO-IMMUNE-001 |
| `DetectFocalAmplifications` | ONCO-CNA-002 |
| `DetectFusions` | ONCO-FUSION-001 |
| `DetectHomozygousDeletions` | ONCO-CNA-003 |
| `DetectLOH` | ONCO-LOH-001 |
| `DetectMRD` | ONCO-MRD-001 |
| `DetectMSI` | ONCO-MSI-001 |
| `DetectOxoGArtifacts` | ONCO-ARTIFACT-001 |
| `DetectStructuralVariants` | ONCO-SV-001 |
| `DetectWholeGenomeDoubling` | ONCO-PLOIDY-001 |
| `EstimateCCF` | ONCO-CCF-001 |
| `EstimateExposures` | ONCO-SIG-003 |
| `EstimateInfiltration` | ONCO-IMMUNE-001 |
| `EstimatePloidy` | ONCO-PLOIDY-001 |
| `EstimatePurity` | ONCO-PURITY-001 |
| `EstimatePurityFromVAF` | ONCO-PURITY-001 |
| `ExtractSignatures` | ONCO-SIG-001 |
| `FilterArtifacts` | ONCO-ARTIFACT-001 |
| `FilterCHIP` | ONCO-CHIP-001 |
| `FilterGermlineVariants` | ONCO-SOMATIC-001 |
| `FilterSomaticSVs` | ONCO-SV-001 |
| `FindChimericReads` | ONCO-FUSION-001 |
| `FitToCosmicSignatures` | ONCO-SIG-001 |
| `GenerateMutantPeptides` | ONCO-NEO-001 |
| `GetCOSMICAnnotation` | ONCO-ANNOT-001 |
| `DetectHlaLoh` | ONCO-HLA-001 |
| `ParseHlaAllele` | ONCO-HLA-001 |
| `GetFusionAnnotation` | ONCO-FUSION-002 |
| `GetTherapyRecommendations` | ONCO-ACTION-001 |
| `GetTrinucleotideContext` | ONCO-SIG-002 |
| `IdentifyAmplifiedOncogenes` | ONCO-CNA-002 |
| `IdentifyBranchMutations` | ONCO-PHYLO-001 |
| `IdentifyCHIPVariants` | ONCO-CHIP-001 |
| `IdentifyClonalMutations` | ONCO-CLONAL-001 |
| `IdentifyDeletedTumorSuppressors` | ONCO-CNA-003 |
| `IdentifyDriverMutations` | ONCO-DRIVER-001 |
| `IdentifyOutlierGenes` | ONCO-EXPR-001 |
| `IdentifyTrunkMutations` | ONCO-PHYLO-001 |
| `InferSubclones` | ONCO-HETERO-001 |
| `MatchCancerHotspots` | ONCO-DRIVER-001 |
| `MatchKnownFusions` | ONCO-FUSION-002 |
| `PredictFusionProtein` | ONCO-FUSION-003 |
| `PredictMHCBinding` | ONCO-MHC-001 |
| `PredictNeoantigens` | ONCO-NEO-001 |
| `ReconstructPhylogeny` | ONCO-PHYLO-001 |
| `ScoreDriverPotential` | ONCO-DRIVER-001 |
| `ScoreNeoantigens` | ONCO-NEO-001 |
| `SegmentCopyNumber` (CopyNumberAnalyzer) | ONCO-CNA-001 |
| `TrackVariantsOverTime` | ONCO-MRD-001 |
| `TryParseHlaAllele` | ONCO-HLA-001 |
| `ValidateFusion` | ONCO-FUSION-001 |

---

## Appendix B: Complexity Notes

| Test Unit | Claimed | Verified | Notes |
|-----------|---------|----------|-------|
| REP-STR-001 | O(n²) | ⚠️ | Actually O(n × U × R), depends on parameters |
| REP-INV-001 | O(n²) | ⚠️ | O(n² × L) with maxLoopLength |
| CHROM-SYNT-001 | O(n log n) | ⚠️ | Has nested loops, needs verification |
| ALIGN-MULTI-001 | O(k² × m) | ✓ | Anchor-based star alignment; k=sequences, m=avg length |
| META-BIN-001 | O(n) | ⚠️ | K-means is O(n × k × i), verify iterations |
| CRISPR-OFF-001 | O(n × m) | ⚠️ | May be exponential with high mismatches |
| RNA-STRUCT-001 | O(n³) | ✓ | Standard Nussinov/Zuker |
| ASSEMBLY-OLC-001 | O(n² × m) | ⚠️ | Depends on overlap detection method |
| PANGEN-CORE-001 | O(g² × s) | ✓ | All-vs-all k-mer clustering dominates (verified PANGEN-CORE-001) |
| SV-CNV-001 | O(n) | ✓ | One pass over depth + O(w log w) median over w windows (w << n) |

---

## Appendix C: Canonical Implementations

### GC Content (SEQ-GC-001)

```
Canonical: SequenceExtensions.CalculateGcContent(ReadOnlySpan<char>)
         ↑
    ┌────┴────┐
    │         │
CalculateGcContentFast(string)   CalculateGcFraction(ReadOnlySpan)
    │                                      │
    │                              CalculateGcFractionFast(string)
    │                                      │
    ├──────────────────────────────────────┤
    │                                      │
PrimerDesigner.CalculateGcContent   MetagenomicsAnalyzer.CalculateGcContent
ChromosomeAnalyzer (internal)       DnaSequence.GcContent
```

**Test Strategy:** Test canonical `CalculateGcContent(ReadOnlySpan)` thoroughly. Delegates need only smoke tests verifying they call canonical correctly.

---

### DNA Complement (SEQ-COMP-001 + SEQ-REVCOMP-001)

```
Canonical: SequenceExtensions.GetComplementBase(char)
         ↑
    ┌────┴────────────────────┐
    │                         │
TryGetComplement(Span)   TryGetReverseComplement(Span)
    │                         │
DnaSequence.Complement   DnaSequence.ReverseComplement
                              │
                    GetReverseComplementString(string)
```

**Test Strategy:** Test `GetComplementBase` for all nucleotides. Test Span methods for sequence operations. Instance methods are wrappers.

---

## Appendix D: Class Coverage Summary

| Class | Test Units | Status |
|-------|------------|--------|
| ApproximateMatcher | PAT-APPROX-001, PAT-APPROX-002 | ✓ |
| BedParser | PARSE-BED-001 | ✓ |
| ChromosomeAnalyzer | CHROM-TELO-001 to CHROM-SYNT-001 | ✓ |
| CodonOptimizer | CODON-OPT-001 to CODON-USAGE-001 | ✓ |
| ComparativeGenomics | COMPGEN-SYNTENY-001 to COMPGEN-DOTPLOT-001 | ✓ |
| CrisprDesigner | CRISPR-PAM-001 to CRISPR-OFF-001 | ✓ |
| DisorderPredictor | DISORDER-PRED-001 to DISORDER-LC-001 | ✓ |
| EmblParser | PARSE-EMBL-001 | ✓ |
| EpigeneticsAnalyzer | EPIGEN-CPG-001 to EPIGEN-DMR-001 | ✓ |
| FastaParser | PARSE-FASTA-001 | ✓ |
| FastqParser | PARSE-FASTQ-001 | ✓ |
| GcSkewCalculator | SEQ-GCSKEW-001 to SEQ-GC-ANALYSIS-001 | ✓ |
| GenBankParser | PARSE-GENBANK-001 | ✓ |
| GenomeAnnotator | ANNOT-ORF-001 to ANNOT-GFF-001 | ✓ |
| GenomeAssemblyAnalyzer | ASSEMBLY-STATS-001 | ✓ |
| GenomicAnalyzer | REP-TANDEM-001 to GENOMIC-ORF-001 | ✓ |
| GffParser | PARSE-GFF-001 | ✓ |
| GeneticCode | TRANS-CODON-001 | ✓ |
| IupacHelper | PAT-IUPAC-001 | ✓ |
| KmerAnalyzer | KMER-COUNT-001 to KMER-POSITIONS-001 | ✓ |
| MetagenomicsAnalyzer | META-CLASS-001 to META-BIN-001 | ✓ |
| MiRnaAnalyzer | MIRNA-SEED-001 to MIRNA-PRECURSOR-001 | ✓ |
| MotifFinder | PAT-PWM-001 to MOTIF-GENERATE-001 | ✓ |
| PanGenomeAnalyzer | PANGEN-CORE-001, PANGEN-CLUSTER-001 | ✓ |
| PhylogeneticAnalyzer | PHYLO-DIST-001 to PHYLO-COMP-001 | ✓ |
| PopulationGeneticsAnalyzer | POP-FREQ-001 to POP-LD-001 | ✓ |
| PrimerDesigner | PRIMER-TM-001 to PRIMER-STRUCT-001 | ✓ |
| ProbeDesigner | PROBE-DESIGN-001, PROBE-VALID-001 | ✓ |
| ProteinMotifFinder | PROTMOTIF-FIND-001 to PROTMOTIF-COMMON-001 | ✓ |
| QualityScoreAnalyzer | QUALITY-PHRED-001, QUALITY-STATS-001 | ✓ |
| RepeatFinder | REP-STR-001 to REP-PALIN-001 | ✓ |
| RestrictionAnalyzer | RESTR-FIND-001, RESTR-DIGEST-001 | ✓ |
| RnaSecondaryStructure | RNA-STRUCT-001 to RNA-PARTITION-001 | ✓ |
| SequenceAligner | ALIGN-GLOBAL-001 to ALIGN-MULTI-001 | ✓ |
| SequenceAssembler | ASSEMBLY-OLC-001, ASSEMBLY-DBG-001 | ✓ |
| SequenceComplexity | SEQ-COMPLEX-001 to SEQ-COMPLEX-COMPRESS-001 | ✓ |
| SequenceExtensions | SEQ-GC-001 to SEQ-VALID-001 | ✓ |
| SpliceSitePredictor | SPLICE-DONOR-001 to SPLICE-PREDICT-001 | ✓ |
| StructuralVariantAnalyzer | SV-DETECT-001 to SV-CNV-001 | ✓ |
| TranscriptomeAnalyzer | TRANS-EXPR-001 to TRANS-SPLICE-001 | ✓ |
| Translator | TRANS-PROT-001 | ✓ |
| VariantAnnotator | VARIANT-ANNOT-001 | ✓ |
| VariantCaller | VARIANT-CALL-001 to VARIANT-INDEL-001 | ✓ |
| VcfParser | PARSE-VCF-001 | ✓ |
| SequenceStatistics | SEQ-STATS-001 to SEQ-SUMMARY-001 | ✓ |
| CodonUsageAnalyzer | CODON-RSCU-001 to CODON-STATS-001 | ✓ |
| OncologyAnalyzer | ONCO-SOMATIC-001 to ONCO-ANNOT-001, ONCO-TMB-001 to ONCO-LOH-001 | ☐ |
| MutationalSignatures | ONCO-SIG-001 to ONCO-SIG-004 | ☑ |
| OncologyAnalyzer (Fusion) | ONCO-FUSION-001 ☑, ONCO-FUSION-002 ☑, ONCO-FUSION-003 ☑ | ● |
| CopyNumberAnalyzer (→ OncologyAnalyzer) | ONCO-CNA-001 ☑, ONCO-CNA-002 ☑, ONCO-CNA-003 ☑ | ☑ |
| TumorAnalyzer | ONCO-PURITY-001 to ONCO-CLONAL-001 | ☐ |
| NeoantigenPredictor | ONCO-NEO-001 to ONCO-MHC-001 | ☐ |
| ImmuneAnalyzer | ONCO-IMMUNE-001 | ☑ |
| LiquidBiopsyAnalyzer (implemented in OncologyAnalyzer) | ONCO-CTDNA-001 to ONCO-CHIP-001 | ☑ |
| TumorEvolutionAnalyzer (ONCO-PHYLO-001 implemented in OncologyAnalyzer) | ONCO-PHYLO-001 to ONCO-HETERO-001 | ☐ (ONCO-PHYLO-001 ☑) |
| HLA analysis (implemented in OncologyAnalyzer) | ONCO-HLA-001 | ☑ |
| Clinical actionability (implemented in OncologyAnalyzer) | ONCO-ACTION-001 | ☑ |
| StructuralVariantDetector | ONCO-SV-001 | ☐ |
| TumorExpressionAnalyzer (implemented in OncologyAnalyzer) | ONCO-EXPR-001 | ☑ |

**Total Classes Covered: 44/57 (77%)**

---

*Generated: 2026-02-12*
*Checklist version: 2.5 (Oncology Genomics: consistency fixes, 4 new algorithms)*


---

## Test Units by Area — Campaign-Added Algorithms (21 units)

> Net-new algorithms implemented during the limitation-elimination campaign (see `docs/Validation/LIMITATIONS.md`).
> Each is `☐ Not Started` here pending **independent Stage A/B re-validation** to restore `☑` — the code ships and is
> covered by the listed test fixture, but the unit has not yet been re-validated under the project's two-stage protocol.

#### RNA-ACCESS-001: McCaskill Unpaired (Accessibility) Probabilities ☑

| Field | Value |
|-------|-------|
| **Canonical** | `CalculateUnpairedProbabilities` |
| **Area** | RnaStructure |
| **Status** | ☑ Complete (validated 2026-06-25 — Stage A/B PASS, CLEAN) |
| **Evidence** | McCaskill (1990), RNAplfold/Bernhart (2006) |
| **Source** | `RnaSecondaryStructure.cs` |

**Methods:** `CalculateUnpairedProbabilities`, `CalculateRegionUnpairedProbability`

**Tests:** `RnaSecondaryStructure_UnpairedProbabilities_Tests.cs`

---

#### PROTMOTIF-HMM-001: Plan7 Profile-HMM Domain Search ☑

| Field | Value |
|-------|-------|
| **Canonical** | `Plan7ProfileHmm.Viterbi/Forward` |
| **Area** | ProteinMotif |
| **Status** | ☑ Complete (validated 2026-06-25, CLEAN — see reports/PROTMOTIF-HMM-001.md) |
| **Evidence** | Eddy (1998, 2011) HMMER/Plan7, Pfam |
| **Source** | `Plan7ProfileHmm.cs / ProteinMotifFinder.cs` |

**Methods:** `Plan7ProfileHmm.Viterbi/Forward`, `FindDomainsByHmm`, `FindDomainEnvelopes`

**Tests:** `ProteinMotifFinder_FindDomainsByHmm_Tests.cs`

---

#### PRIMER-NNTM-001: Nearest-Neighbour Salt/Mismatch/Dangling-End Tm ☑

| Field | Value |
|-------|-------|
| **Canonical** | `CalculateMeltingTemperatureNN` |
| **Area** | MolTools |
| **Status** | ☑ Validated (Stage A ✅ / Stage B ✅ / CLEAN — 2026-06-25) |
| **Evidence** | SantaLucia (1998), Allawi & SantaLucia (1997), Owczarzy (2004/2008), Bommarito (2000) |
| **Source** | `PrimerDesigner.cs` |

**Methods:** `CalculateMeltingTemperatureNN`, `CalculateMeltingTemperatureNNMismatch`

**Tests:** `PrimerDesigner_NearestNeighborTm_Tests.cs`

---

#### PRIMER-HAIRPIN-001: DNA Hairpin Folder + Secondary-Structure Tm ☑

| Field | Value |
|-------|-------|
| **Canonical** | `FindMostStableHairpin` |
| **Area** | MolTools |
| **Status** | ☑ Complete — validated 2026-06-25 (Stage A ✅ / Stage B ✅ / CLEAN) |
| **Evidence** | SantaLucia (1998), SantaLucia & Hicks (2004), UNAFold DNA params |
| **Source** | `PrimerDesigner.cs` |

**Methods:** `FindMostStableHairpin`, `CalculateHairpinMeltingTemperature`

**Tests:** `PrimerDesigner_HairpinTm_Tests.cs`

---

#### PRIMER-DIMER-001: ntthal Self/Hetero-Dimer Tm ☑

| Field | Value |
|-------|-------|
| **Canonical** | `FindMostStableDimer` |
| **Area** | MolTools |
| **Status** | ☑ Complete (validated 2026-06-25; primer3-py 2.3.0 parity to machine precision) |
| **Evidence** | SantaLucia & Hicks (2004), primer3 ntthal |
| **Source** | `PrimerDesigner.cs / NtthalDimer.cs` |

**Methods:** `FindMostStableDimer`, `CalculateDimerMeltingTemperature`, `CalculateSelfDimerMeltingTemperature`

**Tests:** `PrimerDesigner_DimerTm_Tests.cs`

---

#### PROBE-LNATM-001: LNA-Adjusted NN Tm + MGB Probe Design ☑

| Field | Value |
|-------|-------|
| **Canonical** | `CalculateMeltingTemperatureNNLna` |
| **Area** | MolTools |
| **Status** | ☑ Complete (validated Stage A ✅ / Stage B ✅ / CLEAN — 2026-06-25) |
| **Evidence** | McTigue (2004) LNA NN, Kutyavin (2000) MGB |
| **Source** | `ProbeDesigner.cs` |

**Methods:** `CalculateMeltingTemperatureNNLna`, `EvaluateMgbProbeDesign`

**Tests:** `ProbeDesigner_LnaTm_Tests.cs`

---

#### PROBE-EVALUE-001: Karlin–Altschul Off-Target E-value ☑

| Field | Value |
|-------|-------|
| **Canonical** | `ComputeKarlinAltschul` |
| **Area** | MolTools |
| **Status** | ☑ Complete (Stage A/B validated 2026-06-25 — CLEAN) |
| **Evidence** | Karlin & Altschul (1990), Altschul et al. (1990) BLAST |
| **Source** | `ProbeDesigner.cs` |

**Methods:** `ComputeKarlinAltschul`, `ComputeLambdaNucleotide`

**Tests:** `ProbeDesigner_ProbeValidation_Tests.cs`

---

#### MHC-NN-001: MHCflurry Pan-Allele NN Binding Affinity ☑

| Field | Value |
|-------|-------|
| **Canonical** | `MhcflurryAffinityPredictor.PredictIc50` |
| **Area** | Oncology |
| **Status** | ☑ Complete (Stage A ✅ / Stage B ✅ / CLEAN — 2026-06-25) |
| **Evidence** | O'Donnell et al. (2018, 2020) MHCflurry, Apache-2.0 |
| **Source** | `MhcflurryAffinityPredictor.cs` |

**Methods:** `MhcflurryAffinityPredictor.PredictIc50`, ensemble geometric-mean combiner

**Tests:** `MhcflurryAffinityPredictor_PredictIc50_Tests.cs`

---

#### MHC-MATRIX-001: SMM / BIMAS Matrix pMHC Prediction ☑

| Field | Value |
|-------|-------|
| **Canonical** | `PredictIc50Smm` |
| **Area** | Oncology |
| **Status** | ☑ Complete (validated 2026-06-25; Stage A ✅ / Stage B ✅ / CLEAN) |
| **Evidence** | Peters et al. (2005) SMM, Parker et al. (1994) BIMAS |
| **Source** | `OncologyAnalyzer.cs` |

**Methods:** `PredictIc50Smm`, `PredictBindingHalfLifeBimas`, `PredictAndClassifySmm`

**Tests:** `OncologyAnalyzer_ClassifyMhcBinding_Tests.cs`

---

#### IMMUNE-NUSVR-001: CIBERSORT ν-SVR Immune Deconvolution (ABIS) ☑

| Field | Value |
|-------|-------|
| **Canonical** | `DeconvoluteImmuneCellsNuSvr` |
| **Area** | Oncology |
| **Status** | ☑ Complete (validated 2026-06-25; Stage A/B PASS, CLEAN; sklearn-verified < 1e-6) |
| **Evidence** | Newman et al. (2015) CIBERSORT, Schölkopf et al. (2000) ν-SVR, Monaco et al. (2019) ABIS CC-BY |
| **Source** | `ImmuneAnalyzer.cs` |

**Methods:** `DeconvoluteImmuneCellsNuSvr`, `LoadBundledAbisSignatureMatrix`, `LoadSignatureMatrix`

**Tests:** `ImmuneAnalyzer_ImmuneInfiltration_Tests.cs`

---

#### META-CHECKM-001: CheckM Marker-Gene Completeness/Contamination ☑

| Field | Value |
|-------|-------|
| **Canonical** | `EstimateBinQualityFromMarkers` |
| **Area** | Metagenomics |
| **Status** | ☑ Complete (Stage A/B re-validated 2026-06-25; CLEAN) |
| **Evidence** | Parks et al. (2015) CheckM, Parks et al. (2018) GTDB, Pfam CC0 |
| **Source** | `MetagenomicsAnalyzer.cs` |

**Methods:** `EstimateBinQualityFromMarkers`, `DetectMarkers`, `LoadBundledBacterial/ArchaealMarkerHmms`

**Tests:** `MetagenomicsAnalyzer_MarkerGeneQuality_Tests.cs`

---

#### META-TETRA-001: TETRA Tetranucleotide Z-score Signature ☑

| Field | Value |
|-------|-------|
| **Canonical** | `CalculateTetranucleotideZScores` |
| **Area** | Metagenomics |
| **Status** | ☐ Not Started (pending independent Stage A/B re-validation) |
| **Evidence** | Teeling et al. (2004) TETRA, Schbath (1995) |
| **Source** | `MetagenomicsAnalyzer.cs` |

**Methods:** `CalculateTetranucleotideZScores`, `TetranucleotideZScoreCorrelation`

**Tests:** `MetagenomicsAnalyzer_TetranucleotideZScore_Tests.cs`

---

#### SPLICE-MAXENT3-001: MaxEntScan score3 (3' Acceptor) ☑

| Field | Value |
|-------|-------|
| **Canonical** | `ScoreAcceptorMaxEnt` |
| **Area** | Splicing |
| **Status** | ☑ Complete (independently validated 2026-06-25 — reports/SPLICE-MAXENT3-001.md) |
| **Evidence** | Yeo & Burge (2004), maxentpy (MIT) |
| **Source** | `SpliceSitePredictor.cs` |

**Methods:** `ScoreAcceptorMaxEnt` (MaxEntScan score3)

**Tests:** `SpliceSitePredictor_AcceptorSite_Tests.cs`

---

#### SPLICE-MAXENT5-001: MaxEntScan score5 (5' Donor) ☑

| Field | Value |
|-------|-------|
| **Canonical** | `ScoreDonorMaxEnt` |
| **Area** | Splicing |
| **Status** | ☑ Complete (validated 2026-06-25, CLEAN) |
| **Evidence** | Yeo & Burge (2004), maxentpy (MIT) |
| **Source** | `SpliceSitePredictor.cs` |

**Methods:** `ScoreDonorMaxEnt` (MaxEntScan score5)

**Tests:** `SpliceSitePredictor_DonorSite_Tests.cs`

---

#### MIRNA-CONTEXT-001: TargetScan context++ Scoring ☑

| Field | Value |
|-------|-------|
| **Canonical** | `ScoreTargetSiteContextPlusPlus` |
| **Area** | MiRNA |
| **Status** | ☑ Validated (Stage A ✅ / Stage B ✅ / ✅ CLEAN, 2026-06-25) |
| **Evidence** | Agarwal et al. (2015) TargetScan context++ |
| **Source** | `MiRnaAnalyzer.cs` |

**Methods:** `ScoreTargetSiteContextPlusPlus` (+ SA accessibility wiring)

**Tests:** `MiRnaAnalyzer_TargetPrediction_Tests.cs`

---

#### MIRNA-PCT-001: TargetScan PCT (Branch-Length Conservation) ☑

| Field | Value |
|-------|-------|
| **Canonical** | `PCT (branch-length-score → logistic) from caller alignment+tree` |
| **Area** | MiRNA |
| **Status** | ☑ Complete (independently validated — see `docs/Validation/reports/MIRNA-PCT-001.md`) |
| **Evidence** | Friedman et al. (2009) PCT, TargetScan |
| **Source** | `MiRnaAnalyzer.cs` |

**Methods:** PCT (branch-length-score → logistic) from caller alignment+tree

**Tests:** `MiRnaAnalyzer_TargetPrediction_Tests.cs`

---

#### MIRNA-CLASSIFY-001: Pre-miRNA Structure-Feature Classifier ☑

| Field | Value |
|-------|-------|
| **Canonical** | `ClassifyPreMiRna` |
| **Area** | MiRNA |
| **Status** | ☑ Complete (independently validated 2026-06-25 — Stage A/B PASS, CLEAN) |
| **Evidence** | Bonnet et al. (2004), miRBase (public domain), Zhang (2006) MFEI |
| **Source** | `MiRnaAnalyzer.cs` |

**Methods:** `ClassifyPreMiRna` (logistic over MFE/AMFE/MFEI/GC/%paired)

**Tests:** `MiRnaAnalyzer_PreMiRna_Tests.cs`

---

#### MIRNA-CLEAVAGE-001: Drosha/Dicer Cleavage-Site Prediction ☑

| Field | Value |
|-------|-------|
| **Canonical** | `PredictDroshaDicerCleavage` |
| **Area** | MiRNA |
| **Status** | ☑ Complete (Stage A ✅ PASS / Stage B 🟡 PASS-WITH-NOTES — CLEAN; 3p linear-geometry approximation documented) |
| **Evidence** | Han et al. (2006), Park et al. (2011), Auyeung et al. (2013) |
| **Source** | `MiRnaAnalyzer.cs` |

**Methods:** `PredictDroshaDicerCleavage` (Han 11-bp + Park 22-nt + 2-nt overhang)

**Tests:** `MiRnaAnalyzer_PreMiRna_Tests.cs`

---

#### REP-APPROX-001: Approximate (TRF) Tandem-Repeat Detection ☑

| Field | Value |
|-------|-------|
| **Canonical** | `FindApproximateTandemRepeats` |
| **Area** | Repeats |
| **Status** | ☑ Complete (validated 2026-06-25; see docs/Validation/reports/REP-APPROX-001.md) |
| **Evidence** | Benson (1999) TRF |
| **Source** | `RepeatFinder.cs` |

**Methods:** `FindApproximateTandemRepeats`, `ComputeBernoulliStatistics`

**Tests:** `RepeatFinder_ApproximateTandemRepeats_Tests.cs`

---

#### CHROM-ALPHASAT-001: Alpha-Satellite Monomer Detection ☑

| Field | Value |
|-------|-------|
| **Canonical** | `DetectAlphaSatellite` |
| **Area** | Chromosome |
| **Status** | ☑ Complete (independently validated 2026-06-25; see docs/Validation/reports/CHROM-ALPHASAT-001.md) |
| **Evidence** | Waye & Willard (1987), Henikoff et al. (2001), CENP-B box motif |
| **Source** | `ChromosomeAnalyzer.cs` |

**Methods:** `DetectAlphaSatellite`, `FindCenpBBoxes`

**Tests:** `ChromosomeAnalyzer_AlphaSatellite_Tests.cs`

---

#### CHROM-HOR-001: Higher-Order Repeat (HOR) Detection ☑

| Field | Value |
|-------|-------|
| **Canonical** | `DetectHigherOrderRepeat` |
| **Area** | Chromosome |
| **Status** | ☑ Complete (independently validated 2026-06-25; report `docs/Validation/reports/CHROM-HOR-001.md`) |
| **Evidence** | McNulty & Sullivan (2018), Alkan et al. (2007) |
| **Source** | `ChromosomeAnalyzer.cs` |

**Methods:** `DetectHigherOrderRepeat`

**Tests:** `ChromosomeAnalyzer_HigherOrderRepeat_Tests.cs`

---
