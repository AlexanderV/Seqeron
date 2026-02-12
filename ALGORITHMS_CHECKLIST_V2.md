# Algorithms Checklist v2.0

**Date:** 2026-02-12
**Version:** 2.5
**Library:** Seqeron.Genomics

---

## Quick Reference

| Metric | Value |
|--------|-------|
| **Total Test Units** | 234 |
| **Completed** | 79 |
| **In Progress** | 0 |
| **Blocked** | 0 |
| **Not Started** | 155 |

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
| ☑ | REP-STR-001 | Repeats | 4 | Wikipedia (Microsatellite, Trinucleotide repeat disorder), Richard et al. (2008) | [REP-STR-001.md](TestSpecs/REP-STR-001.md) | RepeatFinder_Microsatellite_Tests.cs |
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
| ☑ | PROBE-DESIGN-001 | MolTools | 3 | Wikipedia (Nucleic acid thermodynamics, Hybridization probe, FISH, DNA microarray), SantaLucia (1998) | [PROBE-DESIGN-001.md](TestSpecs/PROBE-DESIGN-001.md) | ProbeDesigner_ProbeDesign_Tests.cs |
| ☑ | PROBE-VALID-001 | MolTools | 2 | Wikipedia (Hybridization probe, DNA microarray, Off-target genome editing, BLAST), Amann & Ludwig (2000) | [PROBE-VALID-001.md](TestSpecs/PROBE-VALID-001.md) | ProbeDesigner_ProbeValidation_Tests.cs |
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
| ☑ | CHROM-CENT-001 | Chromosome | 1 | Wikipedia (Centromere, Karyotype, Chromosome), Levan (1964) | [CHROM-CENT-001.md](TestSpecs/CHROM-CENT-001.md) | ChromosomeAnalyzer_Centromere_Tests.cs |
| ☑ | CHROM-KARYO-001 | Chromosome | 2 | Wikipedia (Karyotype, Ploidy, Aneuploidy) | [CHROM-KARYO-001.md](TestSpecs/CHROM-KARYO-001.md) | ChromosomeAnalyzer_Karyotype_Tests.cs |
| ☑ | CHROM-ANEU-001 | Chromosome | 2 | Wikipedia (Aneuploidy, Copy Number Variation), Griffiths et al. (2000) | [CHROM-ANEU-001.md](TestSpecs/CHROM-ANEU-001.md) | ChromosomeAnalyzer_Aneuploidy_Tests.cs |
| ☑ | CHROM-SYNT-001 | Chromosome | 2 | Wikipedia (Synteny, Comparative genomics, Chromosomal rearrangement), Wang et al. (2012), Goel et al. (2019) | [CHROM-SYNT-001.md](TestSpecs/CHROM-SYNT-001.md) | ChromosomeAnalyzer_Synteny_Tests.cs |
| ☑ | META-CLASS-001 | Metagenomics | 2 | Wikipedia (Metagenomics), Kraken CCB JHU, Wood & Salzberg (2014) | [META-CLASS-001.md](TestSpecs/META-CLASS-001.md) | MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs |
| ☑ | META-PROF-001 | Metagenomics | 1 | Wikipedia (Metagenomics, Relative Abundance), Shannon (1948), Simpson (1949), Segata et al. (2012) | [META-PROF-001.md](TestSpecs/META-PROF-001.md) | MetagenomicsAnalyzer_TaxonomicProfile_Tests.cs |
| ☑ | META-ALPHA-001 | Metagenomics | 1 | Wikipedia (Diversity Index, Alpha Diversity), Shannon (1948), Simpson (1949), Hill (1973), Chao (1984) | [META-ALPHA-001.md](TestSpecs/META-ALPHA-001.md) | MetagenomicsAnalyzer_AlphaDiversity_Tests.cs |
| ☑ | META-BETA-001 | Metagenomics | 1 | Wikipedia (Beta diversity, Bray-Curtis dissimilarity, Jaccard index), Whittaker (1960), Bray & Curtis (1957), Jaccard (1901) | [META-BETA-001.md](TestSpecs/META-BETA-001.md) | MetagenomicsAnalyzer_BetaDiversity_Tests.cs |
| ☑ | META-BIN-001 | Metagenomics | 1 | Wikipedia (Binning metagenomics), Teeling (2004), Parks et al. (2014), Maguire et al. (2020) | [META-BIN-001.md](TestSpecs/META-BIN-001.md) | MetagenomicsAnalyzer_GenomeBinning_Tests.cs |
| ☑ | CODON-OPT-001 | Codon | 1 | Wikipedia (Codon usage bias, CAI), Sharp & Li (1987), Plotkin & Kudla (2011) | [CODON-OPT-001.md](TestSpecs/CODON-OPT-001.md) | CodonOptimizer_OptimizeSequence_Tests.cs |
| ☑ | CODON-CAI-001 | Codon | 1 | Wikipedia (Codon Adaptation Index), Sharp & Li (1987) | [CODON-CAI-001.md](TestSpecs/CODON-CAI-001.md) | CodonOptimizer_CAI_Tests.cs |
| ☑ | CODON-RARE-001 | Codon | 1 | Wikipedia (Codon usage bias), Kazusa, Shu et al. (2006), Sharp & Li (1987) | [CODON-RARE-001.md](TestSpecs/CODON-RARE-001.md) | CodonOptimizer_FindRareCodons_Tests.cs |
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
| ☑ | MIRNA-TARGET-001 | MiRNA | 2 | Bartel (2009) Cell 136:215-233, Lewis et al. (2005), Grimson et al. (2007), TargetScan 8.0 | [MIRNA-TARGET-001.md](TestSpecs/MIRNA-TARGET-001.md) | MiRnaAnalyzer_TargetPrediction_Tests.cs |
| ☑ | MIRNA-PRECURSOR-001 | MiRNA | 2 | [Evidence](docs/Evidence/MIRNA-PRECURSOR-001-Evidence.md) | [TestSpec](tests/TestSpecs/MIRNA-PRECURSOR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs) |
| ☑ | SPLICE-DONOR-001 | Splicing | 2 | [Evidence](docs/Evidence/SPLICE-DONOR-001-Evidence.md) | [TestSpec](tests/TestSpecs/SPLICE-DONOR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs) |
| ☑ | SPLICE-ACCEPTOR-001 | Splicing | 2 | [Evidence](docs/Evidence/SPLICE-ACCEPTOR-001-Evidence.md) | [TestSpec](tests/TestSpecs/SPLICE-ACCEPTOR-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_AcceptorSite_Tests.cs) |
| ☑ | SPLICE-PREDICT-001 | Splicing | 3 | [Evidence](docs/Evidence/SPLICE-PREDICT-001-Evidence.md) | [TestSpec](tests/TestSpecs/SPLICE-PREDICT-001.md) | [Tests](tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_GeneStructure_Tests.cs) |
| ☐ | DISORDER-PRED-001 | ProteinPred | 3 | - | - | - |
| ☐ | DISORDER-REGION-001 | ProteinPred | 2 | - | - | - |
| ☐ | PROTMOTIF-FIND-001 | ProteinMotif | 3 | - | - | - |
| ☐ | PROTMOTIF-PROSITE-001 | ProteinMotif | 2 | - | - | - |
| ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | 2 | - | - | - |
| ☐ | EPIGEN-CPG-001 | Epigenetics | 3 | - | - | - |
| ☐ | EPIGEN-METHYL-001 | Epigenetics | 3 | - | - | - |
| ☐ | EPIGEN-DMR-001 | Epigenetics | 2 | - | - | - |
| ☐ | VARIANT-CALL-001 | Variants | 3 | - | - | - |
| ☐ | VARIANT-SNP-001 | Variants | 2 | - | - | - |
| ☐ | VARIANT-INDEL-001 | Variants | 2 | - | - | - |
| ☐ | VARIANT-ANNOT-001 | Variants | 2 | - | - | - |
| ☐ | SV-DETECT-001 | StructuralVar | 3 | - | - | - |
| ☐ | SV-BREAKPOINT-001 | StructuralVar | 2 | - | - | - |
| ☐ | SV-CNV-001 | StructuralVar | 2 | - | - | - |
| ☐ | ASSEMBLY-OLC-001 | Assembly | 2 | - | - | - |
| ☐ | ASSEMBLY-DBG-001 | Assembly | 2 | - | - | - |
| ☐ | ASSEMBLY-STATS-001 | Assembly | 4 | - | - | - |
| ☐ | TRANS-EXPR-001 | Transcriptome | 3 | - | - | - |
| ☐ | TRANS-DIFF-001 | Transcriptome | 2 | - | - | - |
| ☐ | TRANS-SPLICE-001 | Transcriptome | 2 | - | - | - |
| ☐ | COMPGEN-SYNTENY-001 | Comparative | 2 | - | - | - |
| ☐ | COMPGEN-ORTHO-001 | Comparative | 2 | - | - | - |
| ☐ | COMPGEN-REARR-001 | Comparative | 2 | - | - | - |
| ☐ | PANGEN-CORE-001 | PanGenome | 2 | - | - | - |
| ☐ | PANGEN-CLUSTER-001 | PanGenome | 2 | - | - | - |
| ☐ | QUALITY-PHRED-001 | Quality | 3 | - | - | - |
| ☐ | QUALITY-STATS-001 | Quality | 2 | - | - | - |
| ☐ | SEQ-STATS-001 | Statistics | 3 | - | - | - |
| ☐ | SEQ-MW-001 | Statistics | 2 | - | - | - |
| ☐ | SEQ-PI-001 | Statistics | 1 | - | - | - |
| ☐ | SEQ-HYDRO-001 | Statistics | 2 | - | - | - |
| ☐ | SEQ-THERMO-001 | Statistics | 2 | - | - | - |
| ☐ | SEQ-DINUC-001 | Statistics | 2 | - | - | - |
| ☐ | SEQ-SECSTRUCT-001 | Statistics | 1 | - | - | - |
| ☐ | CODON-RSCU-001 | Codon | 1 | - | - | - |
| ☐ | CODON-ENC-001 | Codon | 1 | - | - | - |
| ☐ | CODON-STATS-001 | Codon | 1 | - | - | - |
| ☐ | TRANS-SIXFRAME-001 | Translation | 1 | - | - | - |
| ☐ | ASSEMBLY-MERGE-001 | Assembly | 1 | - | - | - |
| ☐ | ASSEMBLY-SCAFFOLD-001 | Assembly | 1 | - | - | - |
| ☐ | ASSEMBLY-COVER-001 | Assembly | 1 | - | - | - |
| ☐ | ASSEMBLY-CONSENSUS-001 | Assembly | 1 | - | - | - |
| ☐ | ASSEMBLY-TRIM-001 | Assembly | 1 | - | - | - |
| ☐ | ASSEMBLY-CORRECT-001 | Assembly | 1 | - | - | - |
| ☐ | PAT-APPROX-003 | Matching | 3 | - | - | - |
| ☐ | ALIGN-STATS-001 | Alignment | 2 | - | - | - |
| ☐ | EPIGEN-BISULF-001 | Epigenetics | 2 | - | - | - |
| ☐ | EPIGEN-CHROM-001 | Epigenetics | 3 | - | - | - |
| ☐ | EPIGEN-AGE-001 | Epigenetics | 1 | - | - | - |
| ☐ | MIRNA-PAIR-001 | MiRNA | 3 | - | - | - |
| ☐ | PANGEN-HEAP-001 | PanGenome | 1 | - | - | - |
| ☐ | PANGEN-MARKER-001 | PanGenome | 2 | - | - | - |
| ☐ | POP-SELECT-001 | PopGen | 2 | - | - | - |
| ☐ | POP-ANCESTRY-001 | PopGen | 1 | - | - | - |
| ☐ | POP-ROH-001 | PopGen | 2 | - | - | - |
| ☐ | META-FUNC-001 | Metagenomics | 2 | - | - | - |
| ☐ | META-RESIST-001 | Metagenomics | 1 | - | - | - |
| ☐ | META-PATHWAY-001 | Metagenomics | 1 | - | - | - |
| ☐ | META-TAXA-001 | Metagenomics | 1 | - | - | - |
| ☐ | PHYLO-BOOT-001 | Phylogenetic | 1 | - | - | - |
| ☐ | PHYLO-STATS-001 | Phylogenetic | 3 | - | - | - |
| ☐ | ANNOT-CODING-001 | Annotation | 1 | - | - | - |
| ☐ | ANNOT-REPEAT-001 | Annotation | 2 | - | - | - |
| ☐ | ANNOT-CODONUSAGE-001 | Annotation | 1 | - | - | - |
| ☐ | RESTR-FILTER-001 | MolTools | 3 | - | - | - |
| ☐ | KMER-DIST-001 | K-mer | 1 | - | - | - |
| ☐ | MOTIF-CONS-001 | Matching | 1 | - | - | - |
| ☐ | GENOMIC-REPEAT-001 | Analysis | 2 | - | - | - |
| ☐ | GENOMIC-COMMON-001 | Analysis | 2 | - | - | - |
| ☐ | GENOMIC-MOTIFS-001 | Analysis | 1 | - | - | - |
| ☐ | SEQ-RNACOMP-001 | Composition | 1 | - | - | - |
| ☐ | PROTMOTIF-PATTERN-001 | ProteinMotif | 4 | - | - | - |
| ☐ | PROTMOTIF-SP-001 | ProteinMotif | 1 | - | - | - |
| ☐ | PROTMOTIF-TM-001 | ProteinMotif | 1 | - | - | - |
| ☐ | PROTMOTIF-CC-001 | ProteinMotif | 1 | - | - | - |
| ☐ | PROTMOTIF-LC-001 | ProteinMotif | 1 | - | - | - |
| ☐ | PROTMOTIF-COMMON-001 | ProteinMotif | 2 | - | - | - |
| ☐ | RNA-PAIR-001 | RnaStructure | 3 | - | - | - |
| ☐ | RNA-HAIRPIN-001 | RnaStructure | 2 | - | - | - |
| ☐ | RNA-MFE-001 | RnaStructure | 2 | - | - | - |
| ☐ | RNA-PSEUDOKNOT-001 | RnaStructure | 1 | - | - | - |
| ☐ | RNA-DOTBRACKET-001 | RnaStructure | 2 | - | - | - |
| ☐ | RNA-INVERT-001 | RnaStructure | 1 | - | - | - |
| ☐ | RNA-PARTITION-001 | RnaStructure | 2 | - | - | - |
| ☐ | SEQ-COMPLEX-KMER-001 | Complexity | 1 | - | - | - |
| ☐ | SEQ-COMPLEX-WINDOW-001 | Complexity | 1 | - | - | - |
| ☐ | SEQ-COMPLEX-DUST-001 | Complexity | 2 | - | - | - |
| ☐ | SEQ-COMPLEX-COMPRESS-001 | Complexity | 1 | - | - | - |
| ☐ | COMPGEN-RBH-001 | Comparative | 1 | - | - | - |
| ☐ | COMPGEN-COMPARE-001 | Comparative | 1 | - | - | - |
| ☐ | COMPGEN-REVERSAL-001 | Comparative | 1 | - | - | - |
| ☐ | COMPGEN-CLUSTER-001 | Comparative | 1 | - | - | - |
| ☐ | COMPGEN-ANI-001 | Comparative | 1 | - | - | - |
| ☐ | COMPGEN-DOTPLOT-001 | Comparative | 1 | - | - | - |
| ☐ | MOTIF-DISCOVER-001 | Matching | 1 | - | - | - |
| ☐ | MOTIF-SHARED-001 | Matching | 1 | - | - | - |
| ☐ | MOTIF-REGULATORY-001 | Matching | 1 | - | - | - |
| ☐ | MOTIF-GENERATE-001 | Matching | 1 | - | - | - |
| ☐ | KMER-ASYNC-001 | K-mer | 1 | - | - | - |
| ☐ | KMER-UNIQUE-001 | K-mer | 1 | - | - | - |
| ☐ | KMER-GENERATE-001 | K-mer | 1 | - | - | - |
| ☐ | KMER-BOTH-001 | K-mer | 1 | - | - | - |
| ☐ | KMER-STATS-001 | K-mer | 1 | - | - | - |
| ☐ | KMER-POSITIONS-001 | K-mer | 1 | - | - | - |
| ☐ | SEQ-ATSKEW-001 | Composition | 1 | - | - | - |
| ☐ | SEQ-REPLICATION-001 | Composition | 1 | - | - | - |
| ☐ | SEQ-GC-ANALYSIS-001 | Composition | 1 | - | - | - |
| ☐ | DISORDER-MORF-001 | ProteinPred | 1 | - | - | - |
| ☐ | DISORDER-PROPENSITY-001 | ProteinPred | 3 | - | - | - |
| ☐ | DISORDER-LC-001 | ProteinPred | 1 | - | - | - |
| ☐ | SEQ-COMPOSITION-001 | Statistics | 2 | - | - | - |
| ☐ | SEQ-TM-001 | Statistics | 2 | - | - | - |
| ☐ | SEQ-ENTROPY-PROFILE-001 | Statistics | 1 | - | - | - |
| ☐ | SEQ-GC-PROFILE-001 | Statistics | 1 | - | - | - |
| ☐ | SEQ-CODON-FREQ-001 | Statistics | 1 | - | - | - |
| ☐ | SEQ-SUMMARY-001 | Statistics | 1 | - | - | - |
| ☐ | GENOMIC-TANDEM-001 | Analysis | 1 | - | - | - |
| ☐ | GENOMIC-SIMILARITY-001 | Analysis | 1 | - | - | - |
| ☐ | GENOMIC-ORF-001 | Analysis | 1 | - | - | - |
| ☐ | ONCO-SOMATIC-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-VAF-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-DRIVER-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-ARTIFACT-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-ANNOT-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-TMB-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-MSI-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-HRD-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-LOH-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-SIG-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-SIG-002 | Oncology | 2 | - | - | - |
| ☐ | ONCO-SIG-003 | Oncology | 2 | - | - | - |
| ☐ | ONCO-SIG-004 | Oncology | 1 | - | - | - |
| ☐ | ONCO-FUSION-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-FUSION-002 | Oncology | 2 | - | - | - |
| ☐ | ONCO-FUSION-003 | Oncology | 2 | - | - | - |
| ☐ | ONCO-CNA-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-CNA-002 | Oncology | 2 | - | - | - |
| ☐ | ONCO-CNA-003 | Oncology | 2 | - | - | - |
| ☐ | ONCO-PURITY-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-PLOIDY-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-CLONAL-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-NEO-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-MHC-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-IMMUNE-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-CTDNA-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-MRD-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-CHIP-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-PHYLO-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-CCF-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-HETERO-001 | Oncology | 2 | - | - | - |
| ☐ | ONCO-HLA-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-ACTION-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-SV-001 | Oncology | 3 | - | - | - |
| ☐ | ONCO-EXPR-001 | Oncology | 3 | - | - | - |

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

**Edge Cases:**
- [ ] No repeats found
- [ ] Entire sequence is one repeat
- [ ] minRepeats = 2 (minimum)
- [ ] Unit length 1-6 (mono to hexa)
- [ ] Cancellation mid-operation

---

#### REP-TANDEM-001: Tandem Repeat Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.FindTandemRepeats(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindTandemRepeats(seq, minUnit, maxUnit, minReps)` | GenomicAnalyzer | Canonical |
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
| `CalculateTajimasD(pi, theta, S)` | PopulationGeneticsAnalyzer | D |
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
| **Algorithm** | Hill (1974) composite LD estimator |
| **Formula** | r² = Cov(X₁,X₂)² / (Var(X₁) × Var(X₂)) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateLD(var1, var2, genotypes)` | PopulationGeneticsAnalyzer | D', r² |
| `FindHaplotypeBlocks(variants)` | PopulationGeneticsAnalyzer | Blocks |

**Reference:** Hill WG (1974) "Estimation of linkage disequilibrium in randomly mating populations" Heredity 33:229-239

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
| `AnalyzeCentromere(chrName, seq, windowSize)` | ChromosomeAnalyzer | Canonical |

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
| `CalculateCAI(codingSeq, codonTable)` | CodonOptimizer | Canonical |

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
| `FindOriginOfReplication(sequence)` | GcSkewCalculator | Origin detection |

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
- [ ] 8mer, 7mer-m8, 7mer-A1, 6mer
- [ ] Supplementary pairing
- [ ] Centered sites

---

#### MIRNA-PRECURSOR-001: Pre-miRNA Detection

| Field | Value |
|------|----------|
| **Canonical** | `MiRnaAnalyzer.FindPreMiRnas(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindPreMiRnas(sequence)` | MiRnaAnalyzer | Canonical |
| `ValidateHairpin(structure)` | MiRnaAnalyzer | Validation |

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

#### DISORDER-REGION-001: Disordered Region Detection

| Field | Value |
|------|----------|
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

#### PROTMOTIF-PROSITE-001: PROSITE Pattern Matching

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.MatchPrositePattern(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `MatchPrositePattern(sequence, pattern)` | ProteinMotifFinder | Canonical |
| `ParsePrositePattern(pattern)` | ProteinMotifFinder | Parse |

**Common Patterns:**
- [ ] N-glycosylation (PS00001)
- [ ] Phosphorylation sites
- [ ] Zinc fingers
- [ ] Signal peptides

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
| `CalculateObservedExpectedCpG(sequence)` | EpigeneticsAnalyzer | O/E ratio |

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

#### VARIANT-ANNOT-001: Variant Annotation

| Field | Value |
|------|----------|
| **Canonical** | `VariantAnnotator.Annotate(...)` |
| **Complexity** | O(v × g) |
| **Class** | VariantAnnotator |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Annotate(variants, annotations)` | VariantAnnotator | Canonical |
| `PredictFunctionalImpact(variant)` | VariantAnnotator | Impact |

---

### 23. Structural Variant Analysis (3 units)

#### SV-DETECT-001: SV Detection

| Field | Value |
|------|----------|
| **Canonical** | `StructuralVariantAnalyzer.DetectSVs(...)` |
| **Complexity** | O(n log n) |
| **Class** | StructuralVariantAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectSVs(readPairs, splitReads)` | StructuralVariantAnalyzer | Canonical |
| `FindDiscordantPairs(readPairs, params)` | StructuralVariantAnalyzer | Discordant |
| `ClassifySV(sv)` | StructuralVariantAnalyzer | Classification |

**SV Types:**
- [ ] Deletion, Duplication, Inversion
- [ ] Insertion, Translocation
- [ ] Complex rearrangements

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
| `DetectCNV(depthData, window)` | StructuralVariantAnalyzer | Canonical |
| `SegmentCopyNumber(logRatios)` | StructuralVariantAnalyzer | Segmentation |

---

### 24. Sequence Assembly (3 units)

#### ASSEMBLY-OLC-001: Overlap-Layout-Consensus

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAssembler.AssembleOLC(...)` |
| **Complexity** | O(n² × m) |
| **Class** | SequenceAssembler |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AssembleOLC(reads, params)` | SequenceAssembler | Canonical |
| `FindAllOverlaps(reads, minOverlap, minId)` | SequenceAssembler | Overlaps |

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

#### ASSEMBLY-STATS-001: Assembly Statistics

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAssemblyAnalyzer.CalculateStatistics(...)` |
| **Complexity** | O(n log n) |
| **Class** | GenomeAssemblyAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateStatistics(contigs)` | GenomeAssemblyAnalyzer | Canonical |
| `CalculateN50(contigs)` | GenomeAssemblyAnalyzer | N50 |
| `CalculateNx(contigs, threshold)` | GenomeAssemblyAnalyzer | Nx/Lx |
| `FindGaps(sequence)` | GenomeAssemblyAnalyzer | Gap detection |

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

#### TRANS-SPLICE-001: Alternative Splicing

| Field | Value |
|------|----------|
| **Canonical** | `TranscriptomeAnalyzer.DetectAlternativeSplicing(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectAlternativeSplicing(isoforms)` | TranscriptomeAnalyzer | Canonical |
| `CalculatePSI(event, reads)` | TranscriptomeAnalyzer | Percent spliced in |

---

### 26. Comparative Genomics (3 units)

#### COMPGEN-SYNTENY-001: Synteny Detection

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.FindSyntenicBlocks(...)` |
| **Complexity** | O(n²) |
| **Class** | ComparativeGenomics |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSyntenicBlocks(genes1, genes2, orthologs)` | ComparativeGenomics | Canonical |
| `VisualizeSynteny(blocks)` | ComparativeGenomics | Visualization |

---

#### COMPGEN-ORTHO-001: Ortholog Identification

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.FindOrthologs(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindOrthologs(genes1, genes2, minIdentity)` | ComparativeGenomics | Canonical |
| `FindParalogs(genes, minIdentity)` | ComparativeGenomics | Paralogs |

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
| `ConstructPanGenome(genomes, idThreshold)` | PanGenomeAnalyzer | Canonical |
| `IdentifyCoreGenes(clusters, threshold)` | PanGenomeAnalyzer | Core genes |

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
| `GeneratePresenceAbsenceMatrix(clusters)` | PanGenomeAnalyzer | Matrix |

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

#### SEQ-MW-001: Molecular Weight Calculation

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateMolecularWeight(...)` |
| **Complexity** | O(n) |

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

#### SEQ-DINUC-001: Dinucleotide Analysis

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

#### ALIGN-STATS-001: Alignment Statistics and Formatting

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAligner.CalculateStatistics(...)` |
| **Complexity** | O(n) |
| **Class** | SequenceAligner |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateStatistics(alignment)` | SequenceAligner | Statistics |
| `FormatAlignment(alignment, lineWidth)` | SequenceAligner | Formatting |

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
| `CalculateEpigeneticAge(methylationProfile, clockType)` | EpigeneticsAnalyzer | Canonical |
| `PredictImprintedGenes(methylationProfile)` | EpigeneticsAnalyzer | Imprinting |

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

---

### 38. Extended Pan-Genome Analysis (2 units)

#### PANGEN-HEAP-001: Pan-Genome Growth Model

| Field | Value |
|------|----------|
| **Canonical** | `PanGenomeAnalyzer.FitHeapsLaw(...)` |
| **Complexity** | O(n) |
| **Class** | PanGenomeAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FitHeapsLaw(panGenomeData)` | PanGenomeAnalyzer | Canonical |
| `CreatePresenceAbsenceMatrix(clusters, genomes)` | PanGenomeAnalyzer | Matrix |

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

#### POP-SELECT-001: Selection Signature Detection

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateIHS(...)` |
| **Complexity** | O(n × h) |
| **Class** | PopulationGeneticsAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateIHS(haplotypes, positions)` | PopulationGeneticsAnalyzer | iHS |
| `ScanForSelection(variants, windowSize)` | PopulationGeneticsAnalyzer | Genome-wide |

---

#### POP-ANCESTRY-001: Ancestry Estimation

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.EstimateAncestry(...)` |
| **Complexity** | O(n × k) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `EstimateAncestry(genotypes, refPanels, k)` | PopulationGeneticsAnalyzer | Canonical |

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
| **Complexity** | O(n × t) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSignificantTaxa(profiles, groups, pThreshold)` | MetagenomicsAnalyzer | Canonical |

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
| `Bootstrap(sequences, nReplicates, treeMethod)` | PhylogeneticAnalyzer | Canonical |

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

#### PROTMOTIF-COMMON-001: Common Motif Finding

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.FindCommonMotifs(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindCommonMotifs(sequence)` | ProteinMotifFinder | Canonical |
| `FindAllKnownMotifs(sequence)` | ProteinMotifFinder | All patterns |

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

---

#### RNA-HAIRPIN-001: Hairpin Energy Calculation

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.CalculateHairpinLoopEnergy(...)` |
| **Complexity** | O(n) |

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

#### RNA-INVERT-001: RNA Inverted Repeats

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.FindInvertedRepeats(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindInvertedRepeats(sequence)` | RnaSecondaryStructure | Potential stems |

---

#### RNA-PARTITION-001: Partition Function

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.CalculateStructureProbability(...)` |
| **Complexity** | O(n³) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateStructureProbability(sequence)` | RnaSecondaryStructure | Probability |
| `GenerateRandomRna(length)` | RnaSecondaryStructure | Random generation |

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
| **Canonical** | `SequenceComplexity.EstimateCompressionRatio(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `EstimateCompressionRatio(sequence)` | SequenceComplexity | Compression-based complexity |

---

### 51. Extended Comparative Genomics (6 units)

#### COMPGEN-RBH-001: Reciprocal Best Hits

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.FindReciprocalBestHits(...)` |
| **Complexity** | O(n × m) |
| **Class** | ComparativeGenomics |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindReciprocalBestHits(genes1, genes2)` | ComparativeGenomics | RBH ortholog identification |

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
| `CalculateReversalDistance(geneOrder1, geneOrder2)` | ComparativeGenomics | Reversal distance |

---

#### COMPGEN-CLUSTER-001: Conserved Gene Clusters

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.FindConservedClusters(...)` |
| **Complexity** | O(n × g) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindConservedClusters(genomes)` | ComparativeGenomics | Multi-genome clusters |

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
| **Canonical** | `DisorderPredictor.GetDisorderPropensity(...)` |
| **Complexity** | O(1) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetDisorderPropensity(aminoAcid)` | DisorderPredictor | Propensity value |
| `IsDisorderPromoting(aminoAcid)` | DisorderPredictor | Boolean check |
| `DisorderPromotingAminoAcids` | DisorderPredictor | Property |

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

#### SEQ-COMPOSITION-001: Sequence Composition

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateNucleotideComposition(...)` |
| **Complexity** | O(n) |
| **Class** | SequenceStatistics |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateNucleotideComposition(sequence)` | SequenceStatistics | Nucleotide composition |
| `CalculateAminoAcidComposition(sequence)` | SequenceStatistics | Amino acid composition |

---

#### SEQ-TM-001: Melting Temperature

| Field | Value |
|------|----------|
| **Canonical** | `SequenceStatistics.CalculateMeltingTemperature(...)` |
| **Complexity** | O(n) |

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

#### GENOMIC-TANDEM-001: Tandem Repeat Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.FindTandemRepeats(...)` |
| **Complexity** | O(n²) |
| **Class** | GenomicAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindTandemRepeats(sequence)` | GenomicAnalyzer | Consecutive repeating units |

---

#### GENOMIC-SIMILARITY-001: Sequence Similarity

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.CalculateSimilarity(...)` |
| **Complexity** | O(n + m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateSimilarity(seq1, seq2, k)` | GenomicAnalyzer | K-mer based similarity |

---

#### GENOMIC-ORF-001: ORF Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.FindOpenReadingFrames(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindOpenReadingFrames(sequence)` | GenomicAnalyzer | Potential ORFs |

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
- [ ] Tumor-only mode (no matched normal)
- [ ] Low tumor purity samples
- [ ] Clonal hematopoiesis contamination

---

#### ONCO-VAF-001: Variant Allele Frequency Analysis

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
| `AdjustVAFForPurity(vaf, purity, ploidy)` | OncologyAnalyzer | Correction |

**Edge Cases:**
- [ ] totalReads = 0 (division by zero)
- [ ] Multiallelic sites (multiple ALT alleles)
- [ ] VAF > 1.0 due to alignment artifacts

---

#### ONCO-DRIVER-001: Driver Mutation Detection

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

#### ONCO-ANNOT-001: Cancer-Specific Variant Annotation

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.AnnotateCancerVariants(...)` |
| **Complexity** | O(n × k) |
| **Depends on** | ONCO-SOMATIC-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnnotateCancerVariants(variants)` | OncologyAnalyzer | Canonical |
| `GetCOSMICAnnotation(variant)` | OncologyAnalyzer | COSMIC lookup |

---

#### ONCO-TMB-001: Tumor Mutational Burden

| Field | Value |
|------|----------|
| **Canonical** | `OncologyAnalyzer.CalculateTMB(...)` |
| **Complexity** | O(n) |
| **Invariant** | TMB = mutations / Mb (coding region) |
| **Depends on** | ONCO-SOMATIC-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateTMB(variants, targetRegionSize)` | OncologyAnalyzer | Canonical |
| `ClassifyTMB(tmb)` | OncologyAnalyzer | Low/Medium/High |

**Thresholds:** Low (<6/Mb), Intermediate (6-20/Mb), High (>20/Mb)

**Edge Cases:**
- [ ] targetRegionSize = 0 (division by zero)
- [ ] Panel < 1 Mb (unstable TMB estimation)
- [ ] Synonymous-only variants (should be excluded)

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
- [ ] Tumor-only mode (no matched normal)
- [ ] Insufficient coverage at microsatellite loci
- [ ] < 5 evaluable markers (unreliable classification)

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
- [ ] Insufficient number of CNV segments for scoring
- [ ] Near-diploid tumors (low signal)

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
| `DetectLOH(tumorVcf, normalVcf)` | OncologyAnalyzer | Canonical |
| `CalculateLOHFraction(chromosome)` | OncologyAnalyzer | Per chromosome |

---

#### ONCO-SIG-001: COSMIC Mutational Signature Extraction

| Field | Value |
|------|----------|
| **Canonical** | `MutationalSignatures.ExtractSignatures(...)` |
| **Complexity** | O(n × k) where k=signature count |
| **Invariant** | sum(exposures) = total_mutations |
| **Depends on** | ONCO-SIG-002 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ExtractSignatures(mutations, refGenome)` | MutationalSignatures | Canonical |
| `FitToCosmicSignatures(spectrum)` | MutationalSignatures | COSMIC v3.3 |
| `DecomposeSpectrum(spectrum, signatures)` | MutationalSignatures | NMF decomposition |

**Signatures:** SBS, DBS, ID (COSMIC v3.3)

---

#### ONCO-SIG-002: Trinucleotide Context Analysis

| Field | Value |
|------|----------|
| **Canonical** | `MutationalSignatures.GetTrinucleotideContext(...)` |
| **Complexity** | O(n) |
| **Depends on** | ONCO-SOMATIC-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetTrinucleotideContext(variant, refGenome)` | MutationalSignatures | Canonical |
| `Build96ChannelSpectrum(variants)` | MutationalSignatures | SBS spectrum |

---

#### ONCO-SIG-003: Signature Exposure Estimation

| Field | Value |
|------|----------|
| **Canonical** | `MutationalSignatures.EstimateExposures(...)` |
| **Complexity** | O(n × k²) |
| **Invariant** | all(exposure ≥ 0) |
| **Depends on** | ONCO-SIG-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `EstimateExposures(spectrum, signatures)` | MutationalSignatures | Canonical |
| `BootstrapConfidenceIntervals(spectrum)` | MutationalSignatures | CI estimation |

---

#### ONCO-SIG-004: Mutational Process Classification

| Field | Value |
|------|----------|
| **Canonical** | `MutationalSignatures.ClassifyMutationalProcess(...)` |
| **Complexity** | O(k) |
| **Depends on** | ONCO-SIG-003 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ClassifyMutationalProcess(exposures)` | MutationalSignatures | Canonical |

**Processes:** Aging (SBS1/5), APOBEC (SBS2/13), Smoking (SBS4), UV (SBS7), MMR deficiency (SBS6/15/20/26)

---

#### ONCO-FUSION-001: Fusion Gene Detection

| Field | Value |
|------|----------|
| **Canonical** | `FusionDetector.DetectFusions(...)` |
| **Complexity** | O(n log n) |
| **Invariant** | gene5p ≠ gene3p |
| **Depends on** | — (standalone) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectFusions(readPairs, splitReads)` | FusionDetector | Canonical |
| `FindChimericReads(bamFile)` | FusionDetector | Read extraction |
| `ValidateFusion(fusion, refGenome)` | FusionDetector | Validation |

**Clinical Fusions:** BCR-ABL, EML4-ALK, ROS1, NTRK, FGFR, RET

**Edge Cases:**
- [ ] Read-through transcripts (false positive fusions)
- [ ] Low supporting read count
- [ ] Reciprocal fusions (same partners, swapped orientation)

---

#### ONCO-FUSION-002: Known Fusion Database Lookup

| Field | Value |
|------|----------|
| **Canonical** | `FusionDetector.MatchKnownFusions(...)` |
| **Complexity** | O(n × k) |
| **Depends on** | ONCO-FUSION-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `MatchKnownFusions(fusion)` | FusionDetector | Canonical |
| `GetFusionAnnotation(gene5p, gene3p)` | FusionDetector | Lookup |

**Databases:** ChimerDB, COSMIC Fusions, Mitelman

---

#### ONCO-FUSION-003: Fusion Breakpoint Analysis

| Field | Value |
|------|----------|
| **Canonical** | `FusionDetector.AnalyzeBreakpoint(...)` |
| **Complexity** | O(n) |
| **Depends on** | ONCO-FUSION-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeBreakpoint(fusion)` | FusionDetector | Canonical |
| `PredictFusionProtein(fusion, transcripts)` | FusionDetector | Protein product |

---

#### ONCO-CNA-001: Cancer Copy Number Segmentation

| Field | Value |
|------|----------|
| **Canonical** | `CopyNumberAnalyzer.SegmentCopyNumber(...)` |
| **Complexity** | O(n log n) |
| **Invariant** | segments cover entire genome without gaps |
| **Depends on** | — (standalone) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `SegmentCopyNumber(logRatios)` | CopyNumberAnalyzer | CBS algorithm |
| `CallCopyNumberStates(segments, purity)` | CopyNumberAnalyzer | State calling |

**Edge Cases:**
- [ ] Noisy log-ratio data (low coverage)
- [ ] Whole-genome doubling (WGD)

---

#### ONCO-CNA-002: Focal Amplification Detection

| Field | Value |
|------|----------|
| **Canonical** | `CopyNumberAnalyzer.DetectFocalAmplifications(...)` |
| **Complexity** | O(n) |
| **Depends on** | ONCO-CNA-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectFocalAmplifications(segments)` | CopyNumberAnalyzer | Canonical |
| `IdentifyAmplifiedOncogenes(amplifications)` | CopyNumberAnalyzer | Gene mapping |

**Oncogenes:** ERBB2/HER2, MYC, EGFR, CCND1, MDM2, CDK4

---

#### ONCO-CNA-003: Homozygous Deletion Detection

| Field | Value |
|------|----------|
| **Canonical** | `CopyNumberAnalyzer.DetectHomozygousDeletions(...)` |
| **Complexity** | O(n) |
| **Depends on** | ONCO-CNA-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectHomozygousDeletions(segments)` | CopyNumberAnalyzer | Canonical |
| `IdentifyDeletedTumorSuppressors(deletions)` | CopyNumberAnalyzer | Gene mapping |

**Tumor Suppressors:** TP53, RB1, CDKN2A, PTEN, BRCA1/2

---

#### ONCO-PURITY-001: Tumor Purity Estimation

| Field | Value |
|------|----------|
| **Canonical** | `TumorAnalyzer.EstimatePurity(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ purity ≤ 1 |
| **Depends on** | ONCO-VAF-001, ONCO-CNA-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `EstimatePurity(variants, cnvSegments)` | TumorAnalyzer | Canonical |
| `EstimatePurityFromVAF(variants)` | TumorAnalyzer | VAF-based |

**Edge Cases:**
- [ ] Purity < 0.1 (below detection limit)
- [ ] No heterozygous SNPs in matched normal
- [ ] High stromal contamination

---

#### ONCO-PLOIDY-001: Tumor Ploidy Estimation

| Field | Value |
|------|----------|
| **Canonical** | `TumorAnalyzer.EstimatePloidy(...)` |
| **Complexity** | O(n) |
| **Invariant** | ploidy > 0 |
| **Depends on** | ONCO-CNA-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `EstimatePloidy(cnvSegments)` | TumorAnalyzer | Canonical |
| `DetectWholeGenomeDoubling(ploidy)` | TumorAnalyzer | WGD detection |

---

#### ONCO-CLONAL-001: Clonal vs Subclonal Classification

| Field | Value |
|------|----------|
| **Canonical** | `TumorAnalyzer.ClassifyClonality(...)` |
| **Complexity** | O(n) |
| **Invariant** | clonal_count + subclonal_count = total_variants |
| **Depends on** | ONCO-CCF-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ClassifyClonality(variants, purity, ploidy)` | TumorAnalyzer | Canonical |
| `IdentifyClonalMutations(ccfValues)` | TumorAnalyzer | CCF ≈ 1.0 |

---

#### ONCO-NEO-001: Neoantigen Prediction

| Field | Value |
|------|----------|
| **Canonical** | `NeoantigenPredictor.PredictNeoantigens(...)` |
| **Complexity** | O(n × m) where m=HLA alleles |
| **Depends on** | ONCO-SOMATIC-001, ONCO-HLA-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictNeoantigens(variants, hlaType)` | NeoantigenPredictor | Canonical |
| `GenerateMutantPeptides(variant, lengths)` | NeoantigenPredictor | 8-11mer |
| `ScoreNeoantigens(peptides, hlaAlleles)` | NeoantigenPredictor | Binding affinity |

**Edge Cases:**
- [ ] Unknown HLA type (incomplete typing)
- [ ] Frameshift mutations (long peptide sequences)
- [ ] Variants in non-coding regions (no peptide generated)

---

#### ONCO-MHC-001: MHC-Peptide Binding Prediction

| Field | Value |
|------|----------|
| **Canonical** | `NeoantigenPredictor.PredictMHCBinding(...)` |
| **Complexity** | O(n × m) |
| **Invariant** | IC50 > 0 |
| **Depends on** | ONCO-NEO-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictMHCBinding(peptide, hlaAllele)` | NeoantigenPredictor | Canonical |
| `CalculateBindingAffinity(peptide, hla)` | NeoantigenPredictor | IC50 prediction |

**HLA Types:** HLA-A, HLA-B, HLA-C (Class I)

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
- [ ] ctDNA fraction < 0.1% (below detection limit)
- [ ] High background of CHIP mutations
- [ ] No matched primary tumor VCF

---

#### ONCO-MRD-001: Minimal Residual Disease Detection

| Field | Value |
|------|----------|
| **Canonical** | `LiquidBiopsyAnalyzer.DetectMRD(...)` |
| **Complexity** | O(n × k) where k=tracked mutations |
| **Depends on** | ONCO-CTDNA-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectMRD(plasmaVcf, tumorMarkers)` | LiquidBiopsyAnalyzer | Canonical |
| `TrackVariantsOverTime(timepoints)` | LiquidBiopsyAnalyzer | Longitudinal |

---

#### ONCO-CHIP-001: Clonal Hematopoiesis Filtering

| Field | Value |
|------|----------|
| **Canonical** | `LiquidBiopsyAnalyzer.FilterCHIP(...)` |
| **Complexity** | O(n) |
| **Depends on** | ONCO-CTDNA-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FilterCHIP(variants, whiteBloodCellVcf)` | LiquidBiopsyAnalyzer | Canonical |
| `IdentifyCHIPVariants(variants)` | LiquidBiopsyAnalyzer | Known CHIP genes |

**CHIP Genes:** DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1

---

#### ONCO-PHYLO-001: Tumor Phylogeny Reconstruction

| Field | Value |
|------|----------|
| **Canonical** | `TumorEvolutionAnalyzer.ReconstructPhylogeny(...)` |
| **Complexity** | O(n² × k) |
| **Depends on** | ONCO-CCF-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ReconstructPhylogeny(multiSampleVariants)` | TumorEvolutionAnalyzer | Canonical |
| `IdentifyTrunkMutations(phylogeny)` | TumorEvolutionAnalyzer | Clonal mutations |
| `IdentifyBranchMutations(phylogeny)` | TumorEvolutionAnalyzer | Subclonal mutations |

---

#### ONCO-CCF-001: Cancer Cell Fraction Estimation

| Field | Value |
|------|----------|
| **Canonical** | `TumorEvolutionAnalyzer.EstimateCCF(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ CCF ≤ 1 |
| **Depends on** | ONCO-VAF-001, ONCO-PURITY-001, ONCO-CNA-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `EstimateCCF(vaf, purity, localCopyNumber)` | TumorEvolutionAnalyzer | Canonical |
| `ClusterCCFValues(ccfValues)` | TumorEvolutionAnalyzer | Subclone inference |

**Edge Cases:**
- [ ] Purity not determined (unknown purity)
- [ ] Multi-copy loci (ambiguous CCF)

---

#### ONCO-HETERO-001: Tumor Heterogeneity Analysis

| Field | Value |
|------|----------|
| **Canonical** | `TumorEvolutionAnalyzer.AnalyzeHeterogeneity(...)` |
| **Complexity** | O(n log n) |
| **Invariant** | ITH_score ≥ 0 |
| **Depends on** | ONCO-CCF-001 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeHeterogeneity(variants, ccfValues)` | TumorEvolutionAnalyzer | Canonical |
| `CalculateITH(ccfDistribution)` | TumorEvolutionAnalyzer | ITH score |
| `InferSubclones(ccfClusters)` | TumorEvolutionAnalyzer | Subclone count |

---

#### ONCO-HLA-001: HLA Typing from NGS Data

| Field | Value |
|------|----------|
| **Canonical** | `HlaTyper.TypeHLA(...)` |
| **Complexity** | O(n × k) where k=HLA alleles in database |
| **Invariant** | result contains HLA-A, HLA-B, HLA-C (2 alleles each) |
| **Depends on** | — (standalone) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `TypeHLA(bamFile)` | HlaTyper | Canonical |
| `ResolveAmbiguities(candidates)` | HlaTyper | Disambiguation |
| `GetFourDigitType(hlaResult)` | HlaTyper | Resolution |

**Databases:** IPD-IMGT/HLA

**Edge Cases:**
- [ ] Homozygous alleles (single allele detected)
- [ ] Novel alleles not in database
- [ ] Low coverage at HLA loci

---

#### ONCO-ACTION-001: Clinical Actionability Assessment

| Field | Value |
|------|----------|
| **Canonical** | `ClinicalInterpreter.AssessActionability(...)` |
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
- [ ] VUS (variants of uncertain significance)
- [ ] Conflicting evidence across databases
- [ ] Off-label therapy recommendations

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
- [ ] Complex rearrangements (chromothripsis)
- [ ] Low-coverage samples (insufficient split reads)
- [ ] Repetitive regions (ambiguous breakpoints)

---

#### ONCO-EXPR-001: Tumor Gene Expression Quantification

| Field | Value |
|------|----------|
| **Canonical** | `TumorExpressionAnalyzer.QuantifyExpression(...)` |
| **Complexity** | O(n × g) where g=gene count |
| **Invariant** | all TPM values ≥ 0, sum(TPM) ≈ 1e6 |
| **Depends on** | — (standalone) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `QuantifyExpression(rnaBam, annotation)` | TumorExpressionAnalyzer | Canonical |
| `NormalizeExpression(counts)` | TumorExpressionAnalyzer | TPM/FPKM normalization |
| `IdentifyOutlierGenes(expression, reference)` | TumorExpressionAnalyzer | Z-score outliers |

**Edge Cases:**
- [ ] Low RNA quality (degraded, 3' bias)
- [ ] No matched normal expression profile
- [ ] Batch effects between samples

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
| `FindOriginOfReplication` | SEQ-GCSKEW-001 |
| `FindOrthologs` | COMPGEN-ORTHO-001 |
| `FindPalindromes` | REP-PALIN-001 |
| `FindPamSites` | CRISPR-PAM-001 |
| `FindParalogs` | COMPGEN-ORTHO-001 |
| `FindPreMiRnas` | MIRNA-PRECURSOR-001 |
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
| `IdentifyCoreGenes` | PANGEN-CORE-001 |
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
| `AssessActionability` | ONCO-ACTION-001 |
| `BootstrapConfidenceIntervals` (MutationalSignatures) | ONCO-SIG-003 |
| `Build96ChannelSpectrum` | ONCO-SIG-002 |
| `CalculateBindingAffinity` | ONCO-MHC-001 |
| `CalculateHRDScore` | ONCO-HRD-001 |
| `CalculateITH` | ONCO-HETERO-001 |
| `CalculateLOHFraction` | ONCO-LOH-001 |
| `CalculateLOHScore` | ONCO-HRD-001 |
| `CalculateLSTScore` | ONCO-HRD-001 |
| `CalculateMSIScore` | ONCO-MSI-001 |
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
| `GetFourDigitType` | ONCO-HLA-001 |
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
| `NormalizeExpression` | ONCO-EXPR-001 |
| `PredictFusionProtein` | ONCO-FUSION-003 |
| `PredictMHCBinding` | ONCO-MHC-001 |
| `PredictNeoantigens` | ONCO-NEO-001 |
| `QuantifyExpression` | ONCO-EXPR-001 |
| `ReconstructPhylogeny` | ONCO-PHYLO-001 |
| `ResolveAmbiguities` (HlaTyper) | ONCO-HLA-001 |
| `ScoreDriverPotential` | ONCO-DRIVER-001 |
| `ScoreNeoantigens` | ONCO-NEO-001 |
| `SegmentCopyNumber` (CopyNumberAnalyzer) | ONCO-CNA-001 |
| `TrackVariantsOverTime` | ONCO-MRD-001 |
| `TypeHLA` | ONCO-HLA-001 |
| `ValidateFusion` | ONCO-FUSION-001 |

---

## Appendix B: Complexity Notes

| Test Unit | Claimed | Verified | Notes |
|-----------|---------|----------|-------|
| REP-STR-001 | O(n²) | ⚠️ | Actually O(n × U × R), depends on parameters |
| REP-INV-001 | O(n²) | ⚠️ | O(n² × L) with maxLoopLength |
| CHROM-SYNT-001 | O(n log n) | ⚠️ | Has nested loops, needs verification |
| ALIGN-MULTI-001 | O(n² × m²) | ⚠️ | Progressive is typically O(n² × m) |
| META-BIN-001 | O(n) | ⚠️ | K-means is O(n × k × i), verify iterations |
| CRISPR-OFF-001 | O(n × m) | ⚠️ | May be exponential with high mismatches |
| RNA-STRUCT-001 | O(n³) | ✓ | Standard Nussinov/Zuker |
| ASSEMBLY-OLC-001 | O(n² × m) | ⚠️ | Depends on overlap detection method |
| PANGEN-CORE-001 | O(g² × s) | ⚠️ | All-vs-all comparison |
| SV-CNV-001 | O(n) | ⚠️ | Segmentation may be O(n log n) |

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
| MutationalSignatures | ONCO-SIG-001 to ONCO-SIG-004 | ☐ |
| FusionDetector | ONCO-FUSION-001 to ONCO-FUSION-003 | ☐ |
| CopyNumberAnalyzer | ONCO-CNA-001 to ONCO-CNA-003 | ☐ |
| TumorAnalyzer | ONCO-PURITY-001 to ONCO-CLONAL-001 | ☐ |
| NeoantigenPredictor | ONCO-NEO-001 to ONCO-MHC-001 | ☐ |
| ImmuneAnalyzer | ONCO-IMMUNE-001 | ☐ |
| LiquidBiopsyAnalyzer | ONCO-CTDNA-001 to ONCO-CHIP-001 | ☐ |
| TumorEvolutionAnalyzer | ONCO-PHYLO-001 to ONCO-HETERO-001 | ☐ |
| HlaTyper | ONCO-HLA-001 | ☐ |
| ClinicalInterpreter | ONCO-ACTION-001 | ☐ |
| StructuralVariantDetector | ONCO-SV-001 | ☐ |
| TumorExpressionAnalyzer | ONCO-EXPR-001 | ☐ |

**Total Classes Covered: 44/57 (77%)**

---

*Generated: 2026-02-12*
*Checklist version: 2.5 (Oncology Genomics: consistency fixes, 4 new algorithms)*
