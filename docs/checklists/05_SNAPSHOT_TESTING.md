# Checklist 05: Snapshot / Approval Testing (Verify)

**Priority:** P1  
**Framework:** Verify + VerifyNUnit  
**Date:** 2026-03-19  
**Total algorithms:** 258

---

## Description

Snapshot (approval) тести серіалізують повний вихід алгоритму та порівнюють із зафіксованим `.verified.txt` golden master файлом. Будь-яка зміна виходу вимагає явного схвалення розробника (оновлення снепшоту). Ловить ненавмисні регресії у складних структурованих виходах — вирівнювання, дерева, анотації, результати парсингу.

**Поточне покриття:** ~20 snapshot-файлів у `Snapshots/` покривають: Alignment, Annotation, Chromosome, Codon, CRISPR, Disorder, Epigenetics, FileIO, Metagenomics, MiRNA, MolTools, PatternMatching, Phylogenetic, Population, PrimerProbe, ProteinMotif, Repeat, Restriction, RNA, Splicing.

**Процес:**
1. Перший запуск: створює `.verified.txt` файл
2. Наступні запуски: порівнюють із `.verified.txt`, fail при diff
3. При зміні поведінки: ревʼю diff та accept нового снепшоту

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | What to Snapshot | Snapshot File |
|---|--------|-----------|------|-----------------|---------------|
| 1 | ☐ | SEQ-GC-001 | Composition | GC%, GcFraction для набору відомих послідовностей | CompositionSnapshotTests.cs (new) |
| 2 | ☐ | SEQ-COMP-001 | Composition | Complement для відомих послідовностей | CompositionSnapshotTests.cs (new) |
| 3 | ☐ | SEQ-REVCOMP-001 | Composition | Reverse-complement для відомих послідовностей | CompositionSnapshotTests.cs (new) |
| 4 | ☐ | SEQ-VALID-001 | Composition | Validation результати для набору valid/invalid входів | CompositionSnapshotTests.cs (new) |
| 5 | ☐ | SEQ-COMPLEX-001 | Composition | Complexity scores для відомих послідовностей | CompositionSnapshotTests.cs (new) |
| 6 | ☐ | SEQ-ENTROPY-001 | Composition | Shannon entropy для відомих розподілів | CompositionSnapshotTests.cs (new) |
| 7 | ☐ | SEQ-GCSKEW-001 | Composition | Cumulative GC skew масив | CompositionSnapshotTests.cs (new) |
| 8 | ☐ | PAT-EXACT-001 | Matching | Exact motif positions та count | PatternMatchingSnapshotTests.cs |
| 9 | ☐ | PAT-APPROX-001 | Matching | Hamming matches для відомої пари | PatternMatchingSnapshotTests.cs |
| 10 | ☐ | PAT-APPROX-002 | Matching | Edit distance matches для відомої пари | PatternMatchingSnapshotTests.cs |
| 11 | ☐ | PAT-IUPAC-001 | Matching | Degenerate motif matches | PatternMatchingSnapshotTests.cs |
| 12 | ☐ | PAT-PWM-001 | Matching | PWM creation та scan результати | PatternMatchingSnapshotTests.cs |
| 13 | ☐ | REP-STR-001 | Repeats | Microsatellite повний output | RepeatSnapshotTests.cs |
| 14 | ☐ | REP-TANDEM-001 | Repeats | Tandem repeat повний output | RepeatSnapshotTests.cs |
| 15 | ☐ | REP-INV-001 | Repeats | Inverted repeat повний output | RepeatSnapshotTests.cs |
| 16 | ☐ | REP-DIRECT-001 | Repeats | Direct repeat повний output | RepeatSnapshotTests.cs |
| 17 | ☐ | REP-PALIN-001 | Repeats | Palindrome повний output | RepeatSnapshotTests.cs |
| 18 | ☐ | CRISPR-PAM-001 | MolTools | PAM sites SpCas9 результати | CrisprSnapshotTests.cs |
| 19 | ☐ | CRISPR-GUIDE-001 | MolTools | Guide RNA design результати | CrisprSnapshotTests.cs |
| 20 | ☐ | CRISPR-OFF-001 | MolTools | Off-target analysis повний output | CrisprSnapshotTests.cs |
| 21 | ☐ | PRIMER-TM-001 | MolTools | Melting temperature результати | PrimerProbeSnapshotTests.cs |
| 22 | ☐ | PRIMER-DESIGN-001 | MolTools | Primer design результати | PrimerProbeSnapshotTests.cs |
| 23 | ☐ | PRIMER-STRUCT-001 | MolTools | Primer structure analysis (hairpin, dimer) | PrimerProbeSnapshotTests.cs |
| 24 | ☐ | PROBE-DESIGN-001 | MolTools | Probe design результати | PrimerProbeSnapshotTests.cs |
| 25 | ☐ | PROBE-VALID-001 | MolTools | Probe validation результати | PrimerProbeSnapshotTests.cs |
| 26 | ☐ | RESTR-FIND-001 | MolTools | Restriction site positions | RestrictionSnapshotTests.cs |
| 27 | ☐ | RESTR-DIGEST-001 | MolTools | Digest fragments та summary | RestrictionSnapshotTests.cs |
| 28 | ☐ | ANNOT-ORF-001 | Annotation | ORF positions та frames | AnnotationSnapshotTests.cs |
| 29 | ☐ | ANNOT-GENE-001 | Annotation | Gene prediction результати | AnnotationSnapshotTests.cs |
| 30 | ☐ | ANNOT-PROM-001 | Annotation | Promoter motif positions | AnnotationSnapshotTests.cs |
| 31 | ☐ | ANNOT-GFF-001 | Annotation | GFF3 serialization output | AnnotationSnapshotTests.cs |
| 32 | ☐ | KMER-COUNT-001 | K-mer | K-mer count table для відомої послідовності | KmerSnapshotTests.cs (new) |
| 33 | ☐ | KMER-FREQ-001 | K-mer | K-mer frequency distribution | KmerSnapshotTests.cs (new) |
| 34 | ☐ | KMER-FIND-001 | K-mer | Frequent k-mers повний output | KmerSnapshotTests.cs (new) |
| 35 | ☐ | ALIGN-GLOBAL-001 | Alignment | Global alignment: aligned seqs, score, type | AlignmentSnapshotTests.cs |
| 36 | ☐ | ALIGN-LOCAL-001 | Alignment | Local alignment: aligned seqs, score, positions | AlignmentSnapshotTests.cs |
| 37 | ☐ | ALIGN-SEMI-001 | Alignment | Semi-global alignment результати | AlignmentSnapshotTests.cs |
| 38 | ☐ | ALIGN-MULTI-001 | Alignment | Multiple alignment результати | AlignmentSnapshotTests.cs |
| 39 | ☐ | PHYLO-DIST-001 | Phylogenetic | Distance matrix Jukes-Cantor | PhylogeneticSnapshotTests.cs |
| 40 | ☐ | PHYLO-TREE-001 | Phylogenetic | UPGMA та NJ дерева | PhylogeneticSnapshotTests.cs |
| 41 | ☐ | PHYLO-NEWICK-001 | Phylogenetic | Newick string output | PhylogeneticSnapshotTests.cs |
| 42 | ☐ | PHYLO-COMP-001 | Phylogenetic | Tree comparison результати (RF distance) | PhylogeneticSnapshotTests.cs |
| 43 | ☐ | POP-FREQ-001 | PopGen | Allele frequencies | PopulationSnapshotTests.cs |
| 44 | ☐ | POP-DIV-001 | PopGen | Diversity statistics (π, θ, Tajima's D) | PopulationSnapshotTests.cs |
| 45 | ☐ | POP-HW-001 | PopGen | Hardy-Weinberg test результати | PopulationSnapshotTests.cs |
| 46 | ☐ | POP-FST-001 | PopGen | F-statistics повний output | PopulationSnapshotTests.cs |
| 47 | ☐ | POP-LD-001 | PopGen | Linkage disequilibrium результати | PopulationSnapshotTests.cs |
| 48 | ☐ | CHROM-TELO-001 | Chromosome | Telomere analysis результати | ChromosomeSnapshotTests.cs |
| 49 | ☐ | CHROM-CENT-001 | Chromosome | Centromere analysis результати | ChromosomeSnapshotTests.cs |
| 50 | ☐ | CHROM-KARYO-001 | Chromosome | Karyotype analysis результати | ChromosomeSnapshotTests.cs |
| 51 | ☐ | CHROM-ANEU-001 | Chromosome | Aneuploidy detection результати | ChromosomeSnapshotTests.cs |
| 52 | ☐ | CHROM-SYNT-001 | Chromosome | Synteny blocks результати | ChromosomeSnapshotTests.cs |
| 53 | ☐ | META-CLASS-001 | Metagenomics | Taxonomic classification результати | MetagenomicsSnapshotTests.cs |
| 54 | ☐ | META-PROF-001 | Metagenomics | Taxonomic profile | MetagenomicsSnapshotTests.cs |
| 55 | ☐ | META-ALPHA-001 | Metagenomics | Alpha diversity metrics | MetagenomicsSnapshotTests.cs |
| 56 | ☐ | META-BETA-001 | Metagenomics | Beta diversity matrix | MetagenomicsSnapshotTests.cs |
| 57 | ☐ | META-BIN-001 | Metagenomics | Genome binning результати | MetagenomicsSnapshotTests.cs |
| 58 | ☐ | CODON-OPT-001 | Codon | Optimized sequence | CodonSnapshotTests.cs |
| 59 | ☐ | CODON-CAI-001 | Codon | CAI score | CodonSnapshotTests.cs |
| 60 | ☐ | CODON-RARE-001 | Codon | Rare codons list | CodonSnapshotTests.cs |
| 61 | ☐ | CODON-USAGE-001 | Codon | Codon usage table | CodonSnapshotTests.cs |
| 62 | ☐ | TRANS-CODON-001 | Translation | Genetic code table serialization | TranslationSnapshotTests.cs (new) |
| 63 | ☐ | TRANS-PROT-001 | Translation | Translation result | TranslationSnapshotTests.cs (new) |
| 64 | ☐ | PARSE-FASTA-001 | FileIO | Parsed FASTA records | FileIOSnapshotTests.cs |
| 65 | ☐ | PARSE-FASTQ-001 | FileIO | Parsed FASTQ records та statistics | FileIOSnapshotTests.cs |
| 66 | ☐ | PARSE-BED-001 | FileIO | Parsed BED regions | FileIOSnapshotTests.cs |
| 67 | ☐ | PARSE-VCF-001 | FileIO | Parsed VCF variants | FileIOSnapshotTests.cs |
| 68 | ☐ | PARSE-GFF-001 | FileIO | Parsed GFF features та statistics | FileIOSnapshotTests.cs |
| 69 | ☐ | PARSE-GENBANK-001 | FileIO | Parsed GenBank record | FileIOSnapshotTests.cs |
| 70 | ☐ | PARSE-EMBL-001 | FileIO | Parsed EMBL record | FileIOSnapshotTests.cs |
| 71 | ☐ | RNA-STRUCT-001 | RnaStructure | Predicted structure (pairs, dot-bracket) | RnaSnapshotTests.cs |
| 72 | ☐ | RNA-STEMLOOP-001 | RnaStructure | Stem-loop detection результати | RnaSnapshotTests.cs |
| 73 | ☐ | RNA-ENERGY-001 | RnaStructure | MFE computation результати | RnaSnapshotTests.cs |
| 74 | ☐ | MIRNA-SEED-001 | MiRNA | miRNA creation та seed analysis | MiRnaSnapshotTests.cs |
| 75 | ☐ | MIRNA-TARGET-001 | MiRNA | Target site prediction результати | MiRnaSnapshotTests.cs |
| 76 | ☐ | MIRNA-PRECURSOR-001 | MiRNA | Pre-miRNA analysis output | MiRnaSnapshotTests.cs |
| 77 | ☐ | SPLICE-DONOR-001 | Splicing | Donor site detection результати | SplicingSnapshotTests.cs |
| 78 | ☐ | SPLICE-ACCEPTOR-001 | Splicing | Acceptor site detection результати | SplicingSnapshotTests.cs |
| 79 | ☐ | SPLICE-PREDICT-001 | Splicing | Gene structure prediction результати | SplicingSnapshotTests.cs |
| 80 | ☐ | DISORDER-PRED-001 | ProteinPred | Disorder prediction scores | DisorderSnapshotTests.cs |
| 81 | ☐ | DISORDER-REGION-001 | ProteinPred | Disordered regions output | DisorderSnapshotTests.cs |
| 82 | ☐ | PROTMOTIF-FIND-001 | ProteinMotif | Common motifs повний output | ProteinMotifSnapshotTests.cs |
| 83 | ☐ | PROTMOTIF-PROSITE-001 | ProteinMotif | PROSITE pattern matching output | ProteinMotifSnapshotTests.cs |
| 84 | ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | Domain prediction output | ProteinMotifSnapshotTests.cs |
| 85 | ☐ | EPIGEN-CPG-001 | Epigenetics | CpG island detection output | EpigeneticsSnapshotTests.cs |
| 86 | ☐ | ONCO-IMMUNE-001 | Oncology | Immune infiltration analysis output | OncologySnapshotTests.cs (new) |
| 87 | ☐ | ONCO-SOMATIC-001 | Oncology | Somatic mutation call set (full output) | OncologySnapshotTests.cs (new) |
| 88 | ☐ | ONCO-VAF-001 | Oncology | VAF values for known pileups | OncologySnapshotTests.cs (new) |
| 89 | ☐ | ONCO-DRIVER-001 | Oncology | Driver mutation list and scores | OncologySnapshotTests.cs (new) |
| 90 | ☐ | ONCO-ARTIFACT-001 | Oncology | Artifact-filtered variant set | OncologySnapshotTests.cs (new) |
| 91 | ☐ | ONCO-ANNOT-001 | Oncology | Cancer variant annotations | OncologySnapshotTests.cs (new) |
| 92 | ☐ | ONCO-TMB-001 | Oncology | TMB value and classification | OncologySnapshotTests.cs (new) |
| 93 | ☐ | ONCO-MSI-001 | Oncology | MSI status and unstable-loci output | OncologySnapshotTests.cs (new) |
| 94 | ☐ | ONCO-HRD-001 | Oncology | HRD score with component breakdown | OncologySnapshotTests.cs (new) |
| 95 | ☐ | ONCO-LOH-001 | Oncology | LOH regions output | OncologySnapshotTests.cs (new) |
| 96 | ☐ | ONCO-SIG-001 | Oncology | 96-channel SBS context matrix | OncologySnapshotTests.cs (new) |
| 97 | ☐ | ONCO-SIG-002 | Oncology | Signature exposures | OncologySnapshotTests.cs (new) |
| 98 | ☐ | ONCO-SIG-003 | Oncology | Bootstrap exposure confidence intervals | OncologySnapshotTests.cs (new) |
| 99 | ☐ | ONCO-SIG-004 | Oncology | Dominant mutational process classification | OncologySnapshotTests.cs (new) |
| 100 | ☐ | ONCO-FUSION-001 | Oncology | Detected fusions (full output) | OncologySnapshotTests.cs (new) |
| 101 | ☐ | ONCO-FUSION-002 | Oncology | Known-fusion matches | OncologySnapshotTests.cs (new) |
| 102 | ☐ | ONCO-FUSION-003 | Oncology | Breakpoint analysis output | OncologySnapshotTests.cs (new) |
| 103 | ☐ | ONCO-CNA-001 | Oncology | Copy-number segment classification | OncologySnapshotTests.cs (new) |
| 104 | ☐ | ONCO-CNA-002 | Oncology | Focal amplification calls | OncologySnapshotTests.cs (new) |
| 105 | ☐ | ONCO-CNA-003 | Oncology | Homozygous deletion calls | OncologySnapshotTests.cs (new) |
| 106 | ☐ | ONCO-PURITY-001 | Oncology | Tumor purity estimate | OncologySnapshotTests.cs (new) |
| 107 | ☐ | ONCO-PLOIDY-001 | Oncology | Ploidy estimate | OncologySnapshotTests.cs (new) |
| 108 | ☐ | ONCO-CLONAL-001 | Oncology | Clonality classification output | OncologySnapshotTests.cs (new) |
| 109 | ☐ | ONCO-NEO-001 | Oncology | Neoantigen peptide set | OncologySnapshotTests.cs (new) |
| 110 | ☐ | ONCO-MHC-001 | Oncology | MHC binding classification | OncologySnapshotTests.cs (new) |
| 111 | ☐ | ONCO-CTDNA-001 | Oncology | ctDNA analysis output | OncologySnapshotTests.cs (new) |
| 112 | ☐ | ONCO-MRD-001 | Oncology | MRD detection output | OncologySnapshotTests.cs (new) |
| 113 | ☐ | ONCO-CHIP-001 | Oncology | CHIP-filtered variants | OncologySnapshotTests.cs (new) |
| 114 | ☐ | ONCO-PHYLO-001 | Oncology | Tumor phylogeny (Newick + clones) | OncologySnapshotTests.cs (new) |
| 115 | ☐ | ONCO-CCF-001 | Oncology | CCF estimates | OncologySnapshotTests.cs (new) |
| 116 | ☐ | ONCO-HETERO-001 | Oncology | Heterogeneity (MATH) output | OncologySnapshotTests.cs (new) |
| 117 | ☐ | ONCO-HLA-001 | Oncology | HLA typing output | OncologySnapshotTests.cs (new) |
| 118 | ☐ | ONCO-ACTION-001 | Oncology | Actionability tiers | OncologySnapshotTests.cs (new) |
| 119 | ☐ | ONCO-SV-001 | Oncology | Complex rearrangement classification | OncologySnapshotTests.cs (new) |
| 120 | ☐ | ONCO-EXPR-001 | Oncology | Outlier gene list | OncologySnapshotTests.cs (new) |
| 121 | ☐ | SEQ-COMPOSITION-001 | Statistics | Nucleotide composition table | StatisticsSnapshotTests.cs (new) |
| 122 | ☐ | SEQ-DINUC-001 | Statistics | Dinucleotide frequency table | StatisticsSnapshotTests.cs (new) |
| 123 | ☐ | SEQ-HYDRO-001 | Statistics | Hydrophobicity profile | StatisticsSnapshotTests.cs (new) |
| 124 | ☐ | SEQ-MW-001 | Statistics | Molecular weight value | StatisticsSnapshotTests.cs (new) |
| 125 | ☐ | SEQ-PI-001 | Statistics | Isoelectric point value | StatisticsSnapshotTests.cs (new) |
| 126 | ☐ | SEQ-SECSTRUCT-001 | Statistics | Secondary structure assignment | StatisticsSnapshotTests.cs (new) |
| 127 | ☐ | SEQ-STATS-001 | Statistics | Composition statistics output | StatisticsSnapshotTests.cs (new) |
| 128 | ☐ | SEQ-SUMMARY-001 | Statistics | Sequence summary output | StatisticsSnapshotTests.cs (new) |
| 129 | ☐ | SEQ-THERMO-001 | Statistics | Thermodynamic parameters (ΔG/ΔH/ΔS) | StatisticsSnapshotTests.cs (new) |
| 130 | ☐ | SEQ-TM-001 | Statistics | Sequence-level Tm value | StatisticsSnapshotTests.cs (new) |
| 131 | ☐ | COMPGEN-ANI-001 | Comparative | ANI matrix output | ComparativeGenomicsSnapshotTests.cs (new) |
| 132 | ☐ | COMPGEN-CLUSTER-001 | Comparative | Conserved cluster output | ComparativeGenomicsSnapshotTests.cs (new) |
| 133 | ☐ | COMPGEN-COMPARE-001 | Comparative | Genome comparison report | ComparativeGenomicsSnapshotTests.cs (new) |
| 134 | ☐ | COMPGEN-DOTPLOT-001 | Comparative | Dot-plot coordinates | ComparativeGenomicsSnapshotTests.cs (new) |
| 135 | ☐ | COMPGEN-ORTHO-001 | Comparative | Ortholog pairs | ComparativeGenomicsSnapshotTests.cs (new) |
| 136 | ☐ | COMPGEN-RBH-001 | Comparative | Reciprocal-best-hit pairs | ComparativeGenomicsSnapshotTests.cs (new) |
| 137 | ☐ | COMPGEN-REARR-001 | Comparative | Rearrangement breakpoints | ComparativeGenomicsSnapshotTests.cs (new) |
| 138 | ☐ | COMPGEN-REVERSAL-001 | Comparative | Reversal distance value | ComparativeGenomicsSnapshotTests.cs (new) |
| 139 | ☐ | COMPGEN-SYNTENY-001 | Comparative | Syntenic blocks output | ComparativeGenomicsSnapshotTests.cs (new) |
| 140 | ☐ | ASSEMBLY-CONSENSUS-001 | Assembly | Consensus sequence | AssemblySnapshotTests.cs (new) |
| 141 | ☐ | ASSEMBLY-CORRECT-001 | Assembly | Corrected read set | AssemblySnapshotTests.cs (new) |
| 142 | ☐ | ASSEMBLY-COVER-001 | Assembly | Coverage profile | AssemblySnapshotTests.cs (new) |
| 143 | ☐ | ASSEMBLY-DBG-001 | Assembly | De Bruijn contigs | AssemblySnapshotTests.cs (new) |
| 144 | ☐ | ASSEMBLY-MERGE-001 | Assembly | Merged contigs | AssemblySnapshotTests.cs (new) |
| 145 | ☐ | ASSEMBLY-OLC-001 | Assembly | OLC contigs | AssemblySnapshotTests.cs (new) |
| 146 | ☐ | ASSEMBLY-SCAFFOLD-001 | Assembly | Scaffold layout | AssemblySnapshotTests.cs (new) |
| 147 | ☐ | ASSEMBLY-STATS-001 | Assembly | Assembly statistics (N50, L50, ...) | AssemblySnapshotTests.cs (new) |
| 148 | ☐ | ASSEMBLY-TRIM-001 | Assembly | Trimmed read set | AssemblySnapshotTests.cs (new) |
| 149 | ☐ | RNA-DOTBRACKET-001 | RnaStructure | Parsed pair list from dot-bracket | RnaSnapshotTests.cs |
| 150 | ☐ | RNA-HAIRPIN-001 | RnaStructure | Hairpin energy value | RnaSnapshotTests.cs |
| 151 | ☐ | RNA-INVERT-001 | RnaStructure | Inverted-repeat output | RnaSnapshotTests.cs |
| 152 | ☐ | RNA-MFE-001 | RnaStructure | MFE value and structure | RnaSnapshotTests.cs |
| 153 | ☐ | RNA-PAIR-001 | RnaStructure | Pairing matrix | RnaSnapshotTests.cs |
| 154 | ☐ | RNA-PARTITION-001 | RnaStructure | Partition function / pair probabilities | RnaSnapshotTests.cs |
| 155 | ☐ | RNA-PSEUDOKNOT-001 | RnaStructure | Detected pseudoknots | RnaSnapshotTests.cs |
| 156 | ☐ | KMER-ASYNC-001 | K-mer | Async k-mer counts | KmerSnapshotTests.cs (new) |
| 157 | ☐ | KMER-BOTH-001 | K-mer | Both-strand k-mer counts | KmerSnapshotTests.cs (new) |
| 158 | ☐ | KMER-DIST-001 | K-mer | K-mer distance value | KmerSnapshotTests.cs (new) |
| 159 | ☐ | KMER-GENERATE-001 | K-mer | Generated k-mer set | KmerSnapshotTests.cs (new) |
| 160 | ☐ | KMER-POSITIONS-001 | K-mer | K-mer position list | KmerSnapshotTests.cs (new) |
| 161 | ☐ | KMER-STATS-001 | K-mer | K-mer statistics | KmerSnapshotTests.cs (new) |
| 162 | ☐ | KMER-UNIQUE-001 | K-mer | Unique / min-count k-mers | KmerSnapshotTests.cs (new) |
| 163 | ☐ | PROTMOTIF-CC-001 | ProteinMotif | Coiled-coil prediction | ProteinMotifSnapshotTests.cs |
| 164 | ☐ | PROTMOTIF-COMMON-001 | ProteinMotif | Common motifs output | ProteinMotifSnapshotTests.cs |
| 165 | ☐ | PROTMOTIF-LC-001 | ProteinMotif | Low-complexity regions | ProteinMotifSnapshotTests.cs |
| 166 | ☐ | PROTMOTIF-PATTERN-001 | ProteinMotif | Pattern-match output | ProteinMotifSnapshotTests.cs |
| 167 | ☐ | PROTMOTIF-SP-001 | ProteinMotif | Signal-peptide prediction | ProteinMotifSnapshotTests.cs |
| 168 | ☐ | PROTMOTIF-TM-001 | ProteinMotif | Transmembrane helices | ProteinMotifSnapshotTests.cs |
| 169 | ☐ | MOTIF-CONS-001 | Matching | Consensus from alignment | PatternMatchingSnapshotTests.cs |
| 170 | ☐ | MOTIF-DISCOVER-001 | Matching | Discovered motifs | PatternMatchingSnapshotTests.cs |
| 171 | ☐ | MOTIF-GENERATE-001 | Matching | Generated consensus | PatternMatchingSnapshotTests.cs |
| 172 | ☐ | MOTIF-REGULATORY-001 | Matching | Regulatory elements | PatternMatchingSnapshotTests.cs |
| 173 | ☐ | MOTIF-SHARED-001 | Matching | Shared motifs | PatternMatchingSnapshotTests.cs |
| 174 | ☐ | PAT-APPROX-003 | Matching | Best-match output | PatternMatchingSnapshotTests.cs |
| 175 | ☐ | GENOMIC-COMMON-001 | Analysis | Common region output | GenomicAnalyzerSnapshotTests.cs (new) |
| 176 | ☐ | GENOMIC-MOTIFS-001 | Analysis | Known-motif hits | GenomicAnalyzerSnapshotTests.cs (new) |
| 177 | ☐ | GENOMIC-ORF-001 | Analysis | ORF list | GenomicAnalyzerSnapshotTests.cs (new) |
| 178 | ☐ | GENOMIC-REPEAT-001 | Analysis | Repeat list | GenomicAnalyzerSnapshotTests.cs (new) |
| 179 | ☐ | GENOMIC-SIMILARITY-001 | Analysis | Similarity value | GenomicAnalyzerSnapshotTests.cs (new) |
| 180 | ☐ | GENOMIC-TANDEM-001 | Analysis | Tandem-repeat output | GenomicAnalyzerSnapshotTests.cs (new) |
| 181 | ☐ | EPIGEN-AGE-001 | Epigenetics | Epigenetic age value | EpigeneticsSnapshotTests.cs |
| 182 | ☐ | EPIGEN-BISULF-001 | Epigenetics | Bisulfite-converted sequence | EpigeneticsSnapshotTests.cs |
| 183 | ☐ | EPIGEN-CHROM-001 | Epigenetics | Chromatin-state output | EpigeneticsSnapshotTests.cs |
| 184 | ☐ | EPIGEN-DMR-001 | Epigenetics | DMR list | EpigeneticsSnapshotTests.cs |
| 185 | ☐ | EPIGEN-METHYL-001 | Epigenetics | Methylation levels | EpigeneticsSnapshotTests.cs |
| 186 | ☐ | VARIANT-ANNOT-001 | Variants | Variant functional annotations | VariantSnapshotTests.cs (new) |
| 187 | ☐ | VARIANT-CALL-001 | Variants | Called variant set | VariantSnapshotTests.cs (new) |
| 188 | ☐ | VARIANT-INDEL-001 | Variants | Indel calls | VariantSnapshotTests.cs (new) |
| 189 | ☐ | VARIANT-SNP-001 | Variants | SNP calls | VariantSnapshotTests.cs (new) |
| 190 | ☐ | PANGEN-CLUSTER-001 | PanGenome | Gene clusters | PanGenomeSnapshotTests.cs (new) |
| 191 | ☐ | PANGEN-CORE-001 | PanGenome | Core / accessory partition | PanGenomeSnapshotTests.cs (new) |
| 192 | ☐ | PANGEN-HEAP-001 | PanGenome | Heaps' law fit | PanGenomeSnapshotTests.cs (new) |
| 193 | ☐ | PANGEN-MARKER-001 | PanGenome | Selected markers | PanGenomeSnapshotTests.cs (new) |
| 194 | ☐ | META-FUNC-001 | Metagenomics | Functional predictions | MetagenomicsSnapshotTests.cs |
| 195 | ☐ | META-PATHWAY-001 | Metagenomics | Pathway enrichment output | MetagenomicsSnapshotTests.cs |
| 196 | ☐ | META-RESIST-001 | Metagenomics | Resistance-gene hits | MetagenomicsSnapshotTests.cs |
| 197 | ☐ | META-TAXA-001 | Metagenomics | Significant taxa | MetagenomicsSnapshotTests.cs |
| 198 | ☐ | TRANS-DIFF-001 | Transcriptome | Differential expression table | TranscriptomeSnapshotTests.cs (new) |
| 199 | ☐ | TRANS-EXPR-001 | Transcriptome | Expression quantification | TranscriptomeSnapshotTests.cs (new) |
| 200 | ☐ | TRANS-SPLICE-001 | Transcriptome | Alternative splicing output | TranscriptomeSnapshotTests.cs (new) |
| 201 | ☐ | SV-BREAKPOINT-001 | StructuralVar | Breakpoint list | StructuralVariantSnapshotTests.cs (new) |
| 202 | ☐ | SV-CNV-001 | StructuralVar | CNV calls | StructuralVariantSnapshotTests.cs (new) |
| 203 | ☐ | SV-DETECT-001 | StructuralVar | SV call set | StructuralVariantSnapshotTests.cs (new) |
| 204 | ☐ | DISORDER-LC-001 | ProteinPred | Low-complexity regions | DisorderSnapshotTests.cs |
| 205 | ☐ | DISORDER-MORF-001 | ProteinPred | MoRF predictions | DisorderSnapshotTests.cs |
| 206 | ☐ | DISORDER-PROPENSITY-001 | ProteinPred | Disorder propensity profile | DisorderSnapshotTests.cs |
| 207 | ☐ | POP-ANCESTRY-001 | PopGen | Ancestry proportions | PopulationSnapshotTests.cs |
| 208 | ☐ | POP-ROH-001 | PopGen | ROH segments | PopulationSnapshotTests.cs |
| 209 | ☐ | POP-SELECT-001 | PopGen | Selection-signature output | PopulationSnapshotTests.cs |
| 210 | ☐ | SEQ-ATSKEW-001 | Composition | Cumulative AT-skew array | CompositionSnapshotTests.cs (new) |
| 211 | ☐ | SEQ-REPLICATION-001 | Composition | Predicted replication origin | CompositionSnapshotTests.cs (new) |
| 212 | ☐ | SEQ-RNACOMP-001 | Composition | RNA complement output | CompositionSnapshotTests.cs (new) |
| 213 | ☐ | CODON-ENC-001 | Codon | ENC value | CodonSnapshotTests.cs |
| 214 | ☐ | CODON-RSCU-001 | Codon | RSCU table | CodonSnapshotTests.cs |
| 215 | ☐ | CODON-STATS-001 | Codon | Codon usage statistics | CodonSnapshotTests.cs |
| 216 | ☐ | ANNOT-CODING-001 | Annotation | Coding-potential score | AnnotationSnapshotTests.cs |
| 217 | ☐ | ANNOT-CODONUSAGE-001 | Annotation | Codon usage output | AnnotationSnapshotTests.cs |
| 218 | ☐ | ANNOT-REPEAT-001 | Annotation | Repetitive elements | AnnotationSnapshotTests.cs |
| 219 | ☐ | QUALITY-PHRED-001 | Quality | Parsed Phred scores | FileIOSnapshotTests.cs |
| 220 | ☐ | QUALITY-STATS-001 | Quality | Quality statistics | FileIOSnapshotTests.cs |
| 221 | ☐ | PHYLO-BOOT-001 | Phylogenetic | Bootstrap support values | PhylogeneticSnapshotTests.cs |
| 222 | ☐ | PHYLO-STATS-001 | Phylogenetic | Tree statistics | PhylogeneticSnapshotTests.cs |
| 223 | ☐ | TRANS-SIXFRAME-001 | Translation | Six-frame translation | TranslationSnapshotTests.cs (new) |
| 224 | ☐ | RESTR-FILTER-001 | MolTools | Filtered restriction sites | RestrictionSnapshotTests.cs |
| 225 | ☐ | MIRNA-PAIR-001 | MiRNA | miRNA-target alignment | MiRnaSnapshotTests.cs |
| 226 | ☐ | ALIGN-STATS-001 | Alignment | Alignment statistics | AlignmentSnapshotTests.cs |
| 227 | ☐ | SEQ-CODON-FREQ-001 | Statistics | Codon frequency table for known CDS | SequenceStatisticsSnapshotTests.cs (new) |
| 228 | ☐ | SEQ-COMPLEX-COMPRESS-001 | Complexity | Compression ratio for known seqs | SequenceComplexitySnapshotTests.cs (new) |
| 229 | ☐ | SEQ-COMPLEX-DUST-001 | Complexity | DUST scores for known seqs | SequenceComplexitySnapshotTests.cs (new) |
| 230 | ☐ | SEQ-COMPLEX-KMER-001 | Complexity | K-mer entropy for known seqs | SequenceComplexitySnapshotTests.cs (new) |
| 231 | ☐ | SEQ-COMPLEX-WINDOW-001 | Complexity | Windowed complexity profile | SequenceComplexitySnapshotTests.cs (new) |
| 232 | ☐ | SEQ-ENTROPY-PROFILE-001 | Statistics | Windowed entropy profile | SequenceStatisticsSnapshotTests.cs (new) |
| 233 | ☐ | SEQ-GC-ANALYSIS-001 | Composition | GC content analysis report | CompositionSnapshotTests.cs |
| 234 | ☐ | SEQ-GC-PROFILE-001 | Statistics | Windowed GC% profile | SequenceStatisticsSnapshotTests.cs (new) |
| 235 | ☐ | ONCO-ASCAT-001 | Oncology | ASPCF segments + purity/ploidy fit for a reference profile | OncologySnapshotTests.cs (new) |
| 236 | ☐ | RNA-PKPREDICT-001 | Analysis | Dot-bracket [] for known H-type pseudoknots | RnaStructureSnapshotTests.cs (new) |
| 237 | ☐ | RNA-PKRECURSIVE-001 | Analysis | Nested/multiple-knot dot-bracket for designed RNAs | RnaStructureSnapshotTests.cs (new) |
| 238 | ☑ | RNA-ACCESS-001 | RnaStructure | per-base unpaired probability profile | RnaStructureSnapshotTests.cs (new) |
| 239 | ☑ | PROTMOTIF-HMM-001 | ProteinMotif | domain envelope list with bit scores | ProteinMotifSnapshotTests.cs (new) |
| 240 | ☑ | PRIMER-NNTM-001 | MolTools | Tm + ΔH/ΔS breakdown | PrimerSnapshotTests.cs (new) |
| 241 | ☑ | PRIMER-HAIRPIN-001 | MolTools | most-stable hairpin structure + Tm | PrimerSnapshotTests.cs (new) |
| 242 | ☑ | PRIMER-DIMER-001 | MolTools | most-stable dimer alignment + Tm | PrimerSnapshotTests.cs (new) |
| 243 | ☑ | PROBE-LNATM-001 | MolTools | LNA-adjusted Tm + MGB design report | ProbeSnapshotTests.cs (new) |
| 244 | ☑ | PROBE-EVALUE-001 | MolTools | λ, K, bit score, E-value | ProbeSnapshotTests.cs (new) |
| 245 | ☑ | MHC-NN-001 | Oncology | per-peptide IC50 + percentile | OncologySnapshotTests.cs (new) |
| 246 | ☑ | MHC-MATRIX-001 | Oncology | score + binding classification | OncologySnapshotTests.cs (new) |
| 247 | ☑ | IMMUNE-NUSVR-001 | Oncology | per-cell-type fraction vector | ImmuneSnapshotTests.cs (new) |
| 248 | ☑ | META-CHECKM-001 | Metagenomics | completeness/contamination + per-marker counts | MetagenomicsSnapshotTests.cs (new) |
| 249 | ☐ | META-TETRA-001 | Metagenomics | 256-dim tetranucleotide z-score vector | MetagenomicsSnapshotTests.cs (new) |
| 250 | ☐ | SPLICE-MAXENT3-001 | Splicing | score3 value per acceptor window | SpliceSnapshotTests.cs (new) |
| 251 | ☐ | SPLICE-MAXENT5-001 | Splicing | score5 value per donor window | SpliceSnapshotTests.cs (new) |
| 252 | ☐ | MIRNA-CONTEXT-001 | MiRNA | context++ feature vector + score | MiRnaSnapshotTests.cs (new) |
| 253 | ☐ | MIRNA-PCT-001 | MiRNA | branch-length score + PCT | MiRnaSnapshotTests.cs (new) |
| 254 | ☐ | MIRNA-CLASSIFY-001 | MiRNA | feature vector + native/background label | MiRnaSnapshotTests.cs (new) |
| 255 | ☐ | MIRNA-CLEAVAGE-001 | MiRNA | cleavage coordinates + mature sequence | MiRnaSnapshotTests.cs (new) |
| 256 | ☐ | REP-APPROX-001 | Repeats | repeat array (consensus, copies, match/indel%) | RepeatSnapshotTests.cs (new) |
| 257 | ☐ | CHROM-ALPHASAT-001 | Chromosome | monomer period, copies, CENP-B box positions | ChromosomeSnapshotTests.cs (new) |
| 258 | ☐ | CHROM-HOR-001 | Chromosome | HOR period, copy number, inter/intra identity | ChromosomeSnapshotTests.cs (new) |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 255 |
| ☑ Complete | 9 |
| ☐ Not started | 247 |
| New snapshot test files needed | 5 (Composition, Kmer, Translation, Oncology) |
| Existing snapshot files to extend | ~10 |
