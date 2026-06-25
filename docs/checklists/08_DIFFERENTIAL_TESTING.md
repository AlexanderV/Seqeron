# Checklist 08: Differential Testing

**Priority:** P2
**Date:** 2026-03-19
**Total algorithms:** 258

---

## Description

Differential testing порівнює виходи двох незалежних реалізацій одного алгоритму на ідентичних входах. Виявляє тонкі помилки реалізації. У біоінформатиці існують референсні реалізації (Biopython, EMBOSS, etc.) та альтернативні алгоритмічні стратегії.

**Поточне покриття:** `SafeVsUnsafeDifferentialTests` існує лише для SuffixTree. Нуль для Seqeron.Genomics.

**Стратегії:**
- **ALT** = Альтернативний алгоритм (інший підхід до тієї самої задачі)
- **BRUTE** = Brute-force реалізація для малих входів
- **REF** = Порівняння з обчисленим вручну результатом / референсною бібліотекою
- **DUAL** = Дві реалізації всередині проєкту (різні методи одного класу)

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | Impl A (Current) | Impl B (Reference) | Comparison |
|---|--------|-----------|------|------------------|-------------------|------------|
| 1 | ☐ | SEQ-GC-001 | Composition | Span-based counting | LINQ Count-based | Exact GC% |
| 2 | ☐ | SEQ-COMP-001 | Composition | Switch-based complement | Lookup table complement | Exact sequence |
| 3 | ☐ | SEQ-REVCOMP-001 | Composition | Optimized revcomp | Reverse + complement | Exact sequence |
| 4 | ☐ | SEQ-VALID-001 | Composition | Span loop validation | Regex-based validation | Same bool result |
| 5 | ☐ | SEQ-COMPLEX-001 | Composition | Linguistic complexity | Compression ratio | Correlated scores |
| 6 | ☐ | SEQ-ENTROPY-001 | Composition | Optimized Shannon | Naive histogram entropy | |Δ| < ε |
| 7 | ☐ | SEQ-GCSKEW-001 | Composition | Windowed calculation | Per-base cumulative | Same skew array |
| 8 | ☐ | PAT-EXACT-001 | Matching | Suffix tree / KMP | String.IndexOf loop | Same positions |
| 9 | ☐ | PAT-APPROX-001 | Matching | Optimized Hamming | Brute-force char-by-char | Exact distance |
| 10 | ☐ | PAT-APPROX-002 | Matching | DP edit distance | Recursive with memo | Exact distance |
| 11 | ☐ | PAT-IUPAC-001 | Matching | Switch-based IupacHelper | Dictionary-based matching | Same match results |
| 12 | ☐ | PAT-PWM-001 | Matching | Array-based PWM scan | Dictionary-based PWM scan | Same positions + scores |
| 13 | ☐ | REP-STR-001 | Repeats | RepeatFinder algorithm | Regex-based detection | Same microsatellites |
| 14 | ☐ | REP-TANDEM-001 | Repeats | GenomicAnalyzer | RepeatFinder | Same repeats |
| 15 | ☐ | REP-INV-001 | Repeats | RepeatFinder inverted | Brute-force revcomp search | Same arms |
| 16 | ☐ | REP-DIRECT-001 | Repeats | RepeatFinder direct | Substring search | Same positions |
| 17 | ☐ | REP-PALIN-001 | Repeats | RepeatFinder palindrome | Revcomp equality check | Same palindromes |
| 18 | ☐ | CRISPR-PAM-001 | MolTools | CrisprDesigner.FindPamSites | Regex search for PAM | Same positions |
| 19 | ☐ | CRISPR-GUIDE-001 | MolTools | CrisprDesigner.Design | Manual PAM extraction + scoring | Same top guides |
| 20 | ☐ | CRISPR-OFF-001 | MolTools | Off-target scoring | Brute Hamming search + score | Same off-target set |
| 21 | ☐ | PRIMER-TM-001 | MolTools | SantaLucia nearest-neighbor | Basic rule (4°C GC + 2°C AT) | Correlated Tm |
| 22 | ☐ | PRIMER-DESIGN-001 | MolTools | PrimerDesigner | Manual sliding window | Same candidate set |
| 23 | ☐ | PRIMER-STRUCT-001 | MolTools | Thermodynamic hairpin | Self-complement brute scan | Correlated ΔG |
| 24 | ☐ | PROBE-DESIGN-001 | MolTools | ProbeDesigner | Manual Tm-filtered scan | Same candidate set |
| 25 | ☐ | PROBE-VALID-001 | MolTools | ProbeDesigner.Validate | Manual Tm + specificity check | Same pass/fail |
| 26 | ☐ | RESTR-FIND-001 | MolTools | Pattern-based search | Regex-based search | Same positions |
| 27 | ☐ | RESTR-DIGEST-001 | MolTools | RestrictionAnalyzer.Digest | Manual split at positions | Same fragments |
| 28 | ☐ | ANNOT-ORF-001 | Annotation | GenomeAnnotator.FindOrfs | 6-frame ATG…stop scan | Same ORF set |
| 29 | ☐ | ANNOT-GENE-001 | Annotation | GenomeAnnotator.PredictGenes | ORF + RBS brute scan | Same gene set |
| 30 | ☐ | ANNOT-PROM-001 | Annotation | Promoter motif scoring | Regex -10/-35 box | Same positions |
| 31 | ☐ | ANNOT-GFF-001 | Annotation | GFF serializer | String.Format manual build | Same GFF output |
| 32 | ☐ | KMER-COUNT-001 | K-mer | HashMap counting | Sorted array counting | Exact counts |
| 33 | ☐ | KMER-FREQ-001 | K-mer | Optimized frequency | Count / total manual | Same frequencies |
| 34 | ☐ | KMER-FIND-001 | K-mer | KmerAnalyzer.Find | Sorted + top-N | Same k-mers |
| 35 | ☐ | ALIGN-GLOBAL-001 | Alignment | Needleman-Wunsch (current) | Brute-force O(2^n) short seqs | Score equality |
| 36 | ☐ | ALIGN-LOCAL-001 | Alignment | Smith-Waterman (current) | NW with zero-floor | Score equality |
| 37 | ☐ | ALIGN-SEMI-001 | Alignment | Semi-global variant | NW with free end gaps | Score equality |
| 38 | ☐ | ALIGN-MULTI-001 | Alignment | Progressive MSA | Star alignment | Same column count |
| 39 | ☐ | PHYLO-DIST-001 | Phylogenetic | Jukes-Cantor | Kimura 2-parameter | Both ≥ 0, d_K2P ≥ d_JC |
| 40 | ☐ | PHYLO-TREE-001 | Phylogenetic | UPGMA | NJ | Both valid trees |
| 41 | ☐ | PHYLO-NEWICK-001 | Phylogenetic | Serializer | Manual string construction | Same Newick |
| 42 | ☐ | PHYLO-COMP-001 | Phylogenetic | RF distance | Brute bipartition count | Same RF value |
| 43 | ☐ | POP-FREQ-001 | PopGen | Optimized frequency | Manual count / total | Same frequencies |
| 44 | ☐ | POP-DIV-001 | PopGen | π (pairwise) | Watterson θ | Both ≥ 0, ratio = Tajima's D |
| 45 | ☐ | POP-HW-001 | PopGen | Chi-squared test | Manual expected genotypes | Same chi², same p-value |
| 46 | ☐ | POP-FST-001 | PopGen | Weir-Cockerham Fst | Wright Fst | Correlated, both ∈ [0,1] |
| 47 | ☐ | POP-LD-001 | PopGen | Haplotype-based LD | Allele count LD | Same D', r² |
| 48 | ☐ | CHROM-TELO-001 | Chromosome | Repeat pattern search | Regex TTAGGG search | Same regions |
| 49 | ☐ | CHROM-CENT-001 | Chromosome | AT-richness window | GC-poverty window (1-GC) | Same position |
| 50 | ☐ | CHROM-KARYO-001 | Chromosome | Classification algorithm | Manual arm ratio calc | Same classification |
| 51 | ☐ | CHROM-ANEU-001 | Chromosome | Depth-based CN | Ratio-based CN | Correlated CN |
| 52 | ☐ | CHROM-SYNT-001 | Chromosome | Ortholog-based synteny | BLAST-based synteny | Same major blocks |
| 53 | ☐ | META-CLASS-001 | Metagenomics | K-mer classification | LCA classification | Correlated taxonomy |
| 54 | ☐ | META-PROF-001 | Metagenomics | Read-based profile | K-mer-based profile | Correlated abundances |
| 55 | ☐ | META-ALPHA-001 | Metagenomics | Shannon via log2 | Shannon via ln (normalized) | Proportional |
| 56 | ☐ | META-BETA-001 | Metagenomics | Bray-Curtis | Jaccard | Both valid distance metrics |
| 57 | ☐ | META-BIN-001 | Metagenomics | GC + coverage binning | Tetra-nucleotide freq binning | Consistent bins |
| 58 | ☐ | CODON-OPT-001 | Codon | Optimized (current) | Random synonymous | Both translate to same protein |
| 59 | ☐ | CODON-CAI-001 | Codon | Sharp-Li CAI | Manual w_i product / geometric mean | Same CAI |
| 60 | ☐ | CODON-RARE-001 | Codon | Threshold filter | Manual frequency lookup | Same rare set |
| 61 | ☐ | CODON-USAGE-001 | Codon | Optimized counting | Manual triplet scan | Same usage table |
| 62 | ☐ | TRANS-CODON-001 | Translation | GeneticCode class | Hardcoded codon table | Same mappings |
| 63 | ☐ | TRANS-PROT-001 | Translation | Translator.Translate | Manual triplet → AA loop | Same protein |
| 64 | ☐ | PARSE-FASTA-001 | FileIO | FastaParser | Regex line-by-line parser | Same records |
| 65 | ☐ | PARSE-FASTQ-001 | FileIO | FastqParser | 4-line block reader | Same records |
| 66 | ☐ | PARSE-BED-001 | FileIO | BedParser | Tab-split manual | Same regions |
| 67 | ☐ | PARSE-VCF-001 | FileIO | VcfParser | Tab-split manual | Same variants |
| 68 | ☐ | PARSE-GFF-001 | FileIO | GffParser | Tab-split manual | Same features |
| 69 | ☐ | PARSE-GENBANK-001 | FileIO | GenBankParser | State-machine parser | Same record |
| 70 | ☐ | PARSE-EMBL-001 | FileIO | EmblParser | State-machine parser | Same record |
| 71 | ☐ | RNA-STRUCT-001 | RnaStructure | Nussinov | Energy-based MFE | Nussinov pairs ≤ MFE |
| 72 | ☐ | RNA-STEMLOOP-001 | RnaStructure | Stem-loop finder | Complement search + loop check | Same stem-loops |
| 73 | ☐ | RNA-ENERGY-001 | RnaStructure | Nearest-neighbor energy | Manual stack energy sum | Same ΔG |
| 74 | ☐ | MIRNA-SEED-001 | MiRNA | MiRnaAnalyzer.Seed | Substring extraction | Same seed |
| 75 | ☐ | MIRNA-TARGET-001 | MiRNA | Target prediction | Complement search + scoring | Correlated scores |
| 76 | ☐ | MIRNA-PRECURSOR-001 | MiRNA | Precursor analysis | Structure fold + pattern match | Consistent classification |
| 77 | ☐ | SPLICE-DONOR-001 | Splicing | SpliceSitePredictor | PWM-based scoring | Correlated scores |
| 78 | ☐ | SPLICE-ACCEPTOR-001 | Splicing | SpliceSitePredictor | PWM-based scoring | Correlated scores |
| 79 | ☐ | SPLICE-PREDICT-001 | Splicing | Gene structure | Donor+acceptor pairing | Same introns |
| 80 | ☐ | DISORDER-PRED-001 | ProteinPred | DisorderPredictor | Amino acid propensity scan | Correlated scores |
| 81 | ☐ | DISORDER-REGION-001 | ProteinPred | Region extraction | Threshold-based from scores | Same regions |
| 82 | ☐ | PROTMOTIF-FIND-001 | ProteinMotif | ProteinMotifFinder | Regex pattern search | Same motifs |
| 83 | ☐ | PROTMOTIF-PROSITE-001 | ProteinMotif | PROSITE engine | Regex from PROSITE pattern | Same matches |
| 84 | ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | Domain predictor | HMM profile scan | Correlated domains |
| 85 | ☐ | EPIGEN-CPG-001 | Epigenetics | CpG island detector | Sliding window O/E calc | Same islands |
| 86 | ☐ | ONCO-IMMUNE-001 | Oncology | Immune infiltration | ssGSEA manual calc | Correlated scores |
| 87 | ☐ | ONCO-SOMATIC-001 | Oncology | Somatic caller | Mutect2-style manual calc | Concordant call set |
| 88 | ☐ | ONCO-VAF-001 | Oncology | VAF calc | alt/(ref+alt) by hand | Exact match |
| 89 | ☐ | ONCO-DRIVER-001 | Oncology | Driver scoring | OncoKB / dN-dS reference | Concordant driver ranking |
| 90 | ☐ | ONCO-ARTIFACT-001 | Oncology | Artifact filter | strand-bias formula | Same survivors |
| 91 | ☐ | ONCO-ANNOT-001 | Oncology | Variant annotation | VEP / SnpEff consequence | Same consequence class |
| 92 | ☐ | ONCO-TMB-001 | Oncology | TMB calc | count/Mb by hand | Exact match |
| 93 | ☐ | ONCO-MSI-001 | Oncology | MSI detection | MSIsensor logic | Concordant MSI status |
| 94 | ☐ | ONCO-HRD-001 | Oncology | HRD score | LOH+TAI+LST reference | Same composite |
| 95 | ☐ | ONCO-LOH-001 | Oncology | LOH detection | BAF-deviation reference | Same regions |
| 96 | ☐ | ONCO-SIG-001 | Oncology | SBS context | manual trinucleotide assignment | Identical 96-vector |
| 97 | ☐ | ONCO-SIG-002 | Oncology | Signature fit | NNLS reference solver | Correlated exposures |
| 98 | ☐ | ONCO-SIG-003 | Oncology | Bootstrap CI | resampling by hand | Overlapping CIs |
| 99 | ☐ | ONCO-SIG-004 | Oncology | Process classify | argmax reference | Same dominant process |
| 100 | ☐ | ONCO-FUSION-001 | Oncology | Fusion detection | STAR-Fusion-style logic | Concordant fusions |
| 101 | ☐ | ONCO-FUSION-002 | Oncology | Known-fusion match | DB lookup | Same matches |
| 102 | ☐ | ONCO-FUSION-003 | Oncology | Breakpoint frame | manual frame calc | Same in/out-of-frame |
| 103 | ☐ | ONCO-CNA-001 | Oncology | CN classification | log2→CN reference | Same CN class |
| 104 | ☐ | ONCO-CNA-002 | Oncology | Focal amplification | GISTIC-style threshold | Same focal calls |
| 105 | ☐ | ONCO-CNA-003 | Oncology | Homozygous deletion | CN≈0 reference | Same deletions |
| 106 | ☐ | ONCO-PURITY-001 | Oncology | Purity estimate | ASCAT-style reference | Correlated purity |
| 107 | ☐ | ONCO-PLOIDY-001 | Oncology | Ploidy estimate | weighted-CN reference | Correlated ploidy |
| 108 | ☐ | ONCO-CLONAL-001 | Oncology | Clonality | CCF-threshold reference | Same class |
| 109 | ☐ | ONCO-NEO-001 | Oncology | Neoantigen peptides | sliding-window reference | Same peptide set |
| 110 | ☐ | ONCO-MHC-001 | Oncology | MHC binding | NetMHC-style rank | Concordant binders |
| 111 | ☐ | ONCO-CTDNA-001 | Oncology | ctDNA fraction | tumor-read fraction by hand | Correlated fraction |
| 112 | ☐ | ONCO-MRD-001 | Oncology | MRD detection | tumor-informed reference | Same MRD status |
| 113 | ☐ | ONCO-CHIP-001 | Oncology | CHIP filter | gene+VAF reference | Same flags |
| 114 | ☐ | ONCO-PHYLO-001 | Oncology | Tumor phylogeny | parsimony reference | Concordant topology |
| 115 | ☐ | ONCO-CCF-001 | Oncology | CCF estimate | VAF·CN/purity by hand | Match within tolerance |
| 116 | ☐ | ONCO-HETERO-001 | Oncology | Heterogeneity | MATH formula | Exact MATH |
| 117 | ☐ | ONCO-HLA-001 | Oncology | HLA typing | OptiType-style reference | Concordant alleles |
| 118 | ☐ | ONCO-ACTION-001 | Oncology | Actionability | OncoKB tier reference | Same tier |
| 119 | ☐ | ONCO-SV-001 | Oncology | Complex SV class | breakpoint-pattern reference | Same class |
| 120 | ☐ | ONCO-EXPR-001 | Oncology | Outlier genes | z-score reference | Same outliers |
| 121 | ☐ | SEQ-COMPOSITION-001 | Statistics | Composition | manual count/length | Exact match |
| 122 | ☐ | SEQ-DINUC-001 | Statistics | Dinucleotide freq | manual sliding count | Exact match |
| 123 | ☐ | SEQ-HYDRO-001 | Statistics | Hydrophobicity | Kyte-Doolittle by hand | Exact match |
| 124 | ☐ | SEQ-MW-001 | Statistics | Molecular weight | residue-MW sum reference | Match within tolerance |
| 125 | ☐ | SEQ-PI-001 | Statistics | Isoelectric point | Henderson-Hasselbalch reference | Match within tolerance |
| 126 | ☐ | SEQ-SECSTRUCT-001 | Statistics | Secondary structure | Chou-Fasman reference | Same assignment |
| 127 | ☐ | SEQ-STATS-001 | Statistics | Composition stats | manual calc | Exact match |
| 128 | ☐ | SEQ-SUMMARY-001 | Statistics | Sequence summary | manual calc | Exact match |
| 129 | ☐ | SEQ-THERMO-001 | Statistics | Thermodynamics | nearest-neighbour reference | Match within tolerance |
| 130 | ☐ | SEQ-TM-001 | Statistics | Tm | Marmur-Doty/Wallace by hand | Match within tolerance |
| 131 | ☐ | COMPGEN-ANI-001 | Comparative | ANI | fastANI-style reference | Correlated ANI |
| 132 | ☐ | COMPGEN-CLUSTER-001 | Comparative | Conserved clusters | identity-threshold reference | Same clusters |
| 133 | ☐ | COMPGEN-COMPARE-001 | Comparative | Genome comparison | manual metric calc | Concordant metrics |
| 134 | ☐ | COMPGEN-DOTPLOT-001 | Comparative | Dot plot | k-mer match reference | Same coordinates |
| 135 | ☐ | COMPGEN-ORTHO-001 | Comparative | Orthologs | OrthoFinder-style reference | Concordant pairs |
| 136 | ☐ | COMPGEN-RBH-001 | Comparative | RBH | BLAST-RBH reference | Same pairs |
| 137 | ☐ | COMPGEN-REARR-001 | Comparative | Rearrangements | breakpoint reference | Same breakpoints |
| 138 | ☐ | COMPGEN-REVERSAL-001 | Comparative | Reversal distance | GRIMM-style reference | Exact distance |
| 139 | ☐ | COMPGEN-SYNTENY-001 | Comparative | Syntenic blocks | MCScanX-style reference | Concordant blocks |
| 140 | ☐ | ASSEMBLY-CONSENSUS-001 | Assembly | Consensus | majority-vote by hand | Exact match |
| 141 | ☐ | ASSEMBLY-CORRECT-001 | Assembly | Error correction | k-mer spectrum reference | Concordant reads |
| 142 | ☐ | ASSEMBLY-COVER-001 | Assembly | Coverage | total-bases/length by hand | Exact match |
| 143 | ☐ | ASSEMBLY-DBG-001 | Assembly | DBG assembly | Velvet-style reference | Concordant contigs |
| 144 | ☐ | ASSEMBLY-MERGE-001 | Assembly | Contig merge | overlap reference | Same merges |
| 145 | ☐ | ASSEMBLY-OLC-001 | Assembly | OLC assembly | overlap-graph reference | Concordant contigs |
| 146 | ☐ | ASSEMBLY-SCAFFOLD-001 | Assembly | Scaffolding | mate-pair reference | Same layout |
| 147 | ☐ | ASSEMBLY-STATS-001 | Assembly | Assembly stats | QUAST-style reference | Exact N50/L50 |
| 148 | ☐ | ASSEMBLY-TRIM-001 | Assembly | Quality trim | sliding-window reference | Same trimmed reads |
| 149 | ☐ | RNA-DOTBRACKET-001 | RnaStructure | Dot-bracket parse | manual pairing | Same pairs |
| 150 | ☐ | RNA-HAIRPIN-001 | RnaStructure | Hairpin energy | Turner-rule reference | Match within tolerance |
| 151 | ☐ | RNA-INVERT-001 | RnaStructure | Inverted repeats | revcomp scan reference | Same repeats |
| 152 | ☐ | RNA-MFE-001 | RnaStructure | MFE | Nussinov/Zuker reference | Match within tolerance |
| 153 | ☐ | RNA-PAIR-001 | RnaStructure | Base pairing | Watson-Crick+wobble table | Exact match |
| 154 | ☐ | RNA-PARTITION-001 | RnaStructure | Partition function | McCaskill reference | Correlated probabilities |
| 155 | ☐ | RNA-PSEUDOKNOT-001 | RnaStructure | Pseudoknots | crossing-pair reference | Same detections |
| 156 | ☐ | KMER-ASYNC-001 | K-mer | Async count | sync count | Exact match |
| 157 | ☐ | KMER-BOTH-001 | K-mer | Both-strand count | fwd+revcomp by hand | Exact match |
| 158 | ☐ | KMER-DIST-001 | K-mer | K-mer distance | manual metric | Exact match |
| 159 | ☐ | KMER-GENERATE-001 | K-mer | Generate k-mers | Cartesian product | Exact set |
| 160 | ☐ | KMER-POSITIONS-001 | K-mer | Positions | naive scan | Exact positions |
| 161 | ☐ | KMER-STATS-001 | K-mer | K-mer stats | manual calc | Exact match |
| 162 | ☐ | KMER-UNIQUE-001 | K-mer | Unique k-mers | count==1 reference | Exact set |
| 163 | ☐ | PROTMOTIF-CC-001 | ProteinMotif | Coiled-coil | COILS-style reference | Concordant scores |
| 164 | ☐ | PROTMOTIF-COMMON-001 | ProteinMotif | Common motifs | naive enumeration | Same motifs |
| 165 | ☐ | PROTMOTIF-LC-001 | ProteinMotif | Low-complexity | SEG-style reference | Same regions |
| 166 | ☐ | PROTMOTIF-PATTERN-001 | ProteinMotif | Pattern match | regex reference | Same matches |
| 167 | ☐ | PROTMOTIF-SP-001 | ProteinMotif | Signal peptide | SignalP-style reference | Concordant calls |
| 168 | ☐ | PROTMOTIF-TM-001 | ProteinMotif | TM helices | TMHMM-style reference | Concordant helices |
| 169 | ☐ | MOTIF-CONS-001 | Matching | Consensus | majority by hand | Exact consensus |
| 170 | ☐ | MOTIF-DISCOVER-001 | Matching | Motif discovery | enumeration reference | Same motifs |
| 171 | ☐ | MOTIF-GENERATE-001 | Matching | Generate consensus | majority reference | Exact consensus |
| 172 | ☐ | MOTIF-REGULATORY-001 | Matching | Regulatory elements | known-set reference | Same elements |
| 173 | ☐ | MOTIF-SHARED-001 | Matching | Shared motifs | intersection reference | Same set |
| 174 | ☐ | PAT-APPROX-003 | Matching | Best match | brute-force min distance | Exact distance |
| 175 | ☐ | GENOMIC-COMMON-001 | Analysis | Common region | LCS reference | Same region |
| 176 | ☐ | GENOMIC-MOTIFS-001 | Analysis | Known motifs | naive scan | Same hits |
| 177 | ☐ | GENOMIC-ORF-001 | Analysis | ORFs | NCBI ORFfinder logic | Same ORFs |
| 178 | ☐ | GENOMIC-REPEAT-001 | Analysis | Repeats | naive scan | Same repeats |
| 179 | ☐ | GENOMIC-SIMILARITY-001 | Analysis | Similarity | identity by hand | Match within tolerance |
| 180 | ☐ | GENOMIC-TANDEM-001 | Analysis | Tandem repeats | TRF-style reference | Concordant repeats |
| 181 | ☐ | EPIGEN-AGE-001 | Epigenetics | Epigenetic age | Horvath-clock reference | Match within tolerance |
| 182 | ☐ | EPIGEN-BISULF-001 | Epigenetics | Bisulfite | manual C→T conversion | Exact match |
| 183 | ☐ | EPIGEN-CHROM-001 | Epigenetics | Chromatin state | ChromHMM-style reference | Concordant states |
| 184 | ☐ | EPIGEN-DMR-001 | Epigenetics | DMR | t-test reference | Same regions |
| 185 | ☐ | EPIGEN-METHYL-001 | Epigenetics | Methylation level | methylated/total by hand | Exact match |
| 186 | ☐ | VARIANT-ANNOT-001 | Variants | Variant annotation | VEP-style reference | Same impact |
| 187 | ☐ | VARIANT-CALL-001 | Variants | Variant calling | bcftools-style logic | Concordant calls |
| 188 | ☐ | VARIANT-INDEL-001 | Variants | Indel calling | pileup reference | Same indels |
| 189 | ☐ | VARIANT-SNP-001 | Variants | SNP calling | pileup reference | Same SNPs |
| 190 | ☐ | PANGEN-CLUSTER-001 | PanGenome | Gene clusters | CD-HIT-style reference | Concordant clusters |
| 191 | ☐ | PANGEN-CORE-001 | PanGenome | Pan-genome | Roary-style reference | Same core/accessory |
| 192 | ☐ | PANGEN-HEAP-001 | PanGenome | Heaps' law | regression reference | Match within tolerance |
| 193 | ☐ | PANGEN-MARKER-001 | PanGenome | Markers | core-gene reference | Same markers |
| 194 | ☐ | META-FUNC-001 | Metagenomics | Functions | KEGG-style reference | Concordant assignments |
| 195 | ☐ | META-PATHWAY-001 | Metagenomics | Pathway enrichment | hypergeometric reference | Match within tolerance |
| 196 | ☐ | META-RESIST-001 | Metagenomics | Resistance genes | CARD-style reference | Same hits |
| 197 | ☐ | META-TAXA-001 | Metagenomics | Significant taxa | LEfSe-style reference | Concordant taxa |
| 198 | ☐ | TRANS-DIFF-001 | Transcriptome | Differential expr | DESeq-style reference | Correlated FC |
| 199 | ☐ | TRANS-EXPR-001 | Transcriptome | Expression | TPM by hand | Match within tolerance |
| 200 | ☐ | TRANS-SPLICE-001 | Transcriptome | Alt splicing | PSI by hand | Match within tolerance |
| 201 | ☐ | SV-BREAKPOINT-001 | StructuralVar | Breakpoints | split-read reference | Same breakpoints |
| 202 | ☐ | SV-CNV-001 | StructuralVar | CNV | read-depth reference | Concordant CN |
| 203 | ☐ | SV-DETECT-001 | StructuralVar | SV detection | Manta-style reference | Concordant SVs |
| 204 | ☐ | DISORDER-LC-001 | ProteinPred | Low-complexity | SEG-style reference | Same regions |
| 205 | ☐ | DISORDER-MORF-001 | ProteinPred | MoRF | ANCHOR-style reference | Concordant MoRFs |
| 206 | ☐ | DISORDER-PROPENSITY-001 | ProteinPred | Propensity | IUPred-style reference | Correlated scores |
| 207 | ☐ | POP-ANCESTRY-001 | PopGen | Ancestry | ADMIXTURE-style reference | Correlated proportions |
| 208 | ☐ | POP-ROH-001 | PopGen | ROH | PLINK-style reference | Same segments |
| 209 | ☐ | POP-SELECT-001 | PopGen | Selection | iHS/Fst reference | Correlated signal |
| 210 | ☐ | SEQ-ATSKEW-001 | Composition | AT skew | (A−T)/(A+T) by hand | Exact match |
| 211 | ☐ | SEQ-REPLICATION-001 | Composition | Replication origin | skew-minimum reference | Same index |
| 212 | ☐ | SEQ-RNACOMP-001 | Composition | RNA complement | base table | Exact match |
| 213 | ☐ | CODON-ENC-001 | Codon | ENC | Wright reference | Match within tolerance |
| 214 | ☐ | CODON-RSCU-001 | Codon | RSCU | manual calc | Exact match |
| 215 | ☐ | CODON-STATS-001 | Codon | Codon stats | manual calc | Exact match |
| 216 | ☐ | ANNOT-CODING-001 | Annotation | Coding potential | CPC-style reference | Correlated score |
| 217 | ☐ | ANNOT-CODONUSAGE-001 | Annotation | Codon usage | manual calc | Exact match |
| 218 | ☐ | ANNOT-REPEAT-001 | Annotation | Repetitive elements | naive scan | Same elements |
| 219 | ☐ | QUALITY-PHRED-001 | Quality | Phred parse | ASCII−offset by hand | Exact match |
| 220 | ☐ | QUALITY-STATS-001 | Quality | Quality stats | manual calc | Exact match |
| 221 | ☐ | PHYLO-BOOT-001 | Phylogenetic | Bootstrap | resampling reference | Concordant support |
| 222 | ☐ | PHYLO-STATS-001 | Phylogenetic | Tree stats | manual calc | Exact match |
| 223 | ☐ | TRANS-SIXFRAME-001 | Translation | Six-frame | manual translation | Exact frames |
| 224 | ☐ | RESTR-FILTER-001 | MolTools | Filter sites | criteria reference | Same survivors |
| 225 | ☐ | MIRNA-PAIR-001 | MiRNA | miRNA align | seed-pairing reference | Concordant alignment |
| 226 | ☐ | ALIGN-STATS-001 | Alignment | Alignment stats | manual calc | Exact match |
| 227 | ☐ | SEQ-CODON-FREQ-001 | Statistics | Triplet scan | LINQ GroupBy | Exact frequencies |
| 228 | ☐ | SEQ-COMPLEX-COMPRESS-001 | Complexity | Built-in compression estimate | GZip ratio | Ratio within tolerance |
| 229 | ☐ | SEQ-COMPLEX-DUST-001 | Complexity | DUST score | SDUST reference | Concordant score |
| 230 | ☐ | SEQ-COMPLEX-KMER-001 | Complexity | K-mer entropy | manual Shannon calc | Exact entropy |
| 231 | ☐ | SEQ-COMPLEX-WINDOW-001 | Complexity | Sliding window | naive per-window recompute | Identical profile |
| 232 | ☐ | SEQ-ENTROPY-PROFILE-001 | Statistics | Sliding window | naive per-window recompute | Identical profile |
| 233 | ☐ | SEQ-GC-ANALYSIS-001 | Composition | Windowed scan | LINQ Count-based | Exact GC% |
| 234 | ☐ | SEQ-GC-PROFILE-001 | Statistics | Sliding window | naive per-window recompute | Identical profile |
| 235 | ☐ | ONCO-ASCAT-001 | Oncology | ASPCF PCF DP | naive changepoint scan | same breakpoints |
| 236 | ☐ | RNA-PKPREDICT-001 | Analysis | pknotsRG canonical scan | brute-force H-type search | same structure |
| 237 | ☐ | RNA-PKRECURSIVE-001 | Analysis | recursive grammar | exhaustive nested search (small n) | same ΔG |
| 238 | ☑ | RNA-ACCESS-001 | RnaStructure | Boltzmann (McCaskill) DP | brute-force ensemble enumeration (small n) | equal P_unpaired ±1e-9 |
| 239 | ☑ | PROTMOTIF-HMM-001 | ProteinMotif | C# Plan7 (Viterbi/Forward/null2) | pyhmmer / hmmsearch | bit score ±1e-3, same envelopes |
| 240 | ☑ | PRIMER-NNTM-001 | MolTools | C# NN (unified SantaLucia) | primer3-py / Biopython MeltingTemp | Tm ±0.5°C |
| 241 | ☑ | PRIMER-HAIRPIN-001 | MolTools | C# ntthal hairpin (CalculateHairpinThermodynamicsNtthal) | primer3-py 2.3.0 calc_hairpin | ΔH exact, ΔS/Tm ≤1e-6 (machine precision) |
| 242 | ☑ | PRIMER-DIMER-001 | MolTools | C# ntthal DP | primer3-py calc_homodimer/heterodimer | ΔG/Tm to machine precision on contiguous optima |
| 243 | ☑ | PROBE-LNATM-001 | MolTools | C# LNA NN | MELTING 5 | Tm ±0.2°C |
| 244 | ☐ | PROBE-EVALUE-001 | MolTools | C# Karlin-Altschul | NCBI BLAST stats / published λ | λ≈1.374, E within tolerance |
| 245 | ☐ | MHC-NN-001 | Oncology | C# MHCflurry port | mhcflurry 2.1.5 (models_class1_pan) | IC50 < 0.03% |
| 246 | ☐ | MHC-MATRIX-001 | Oncology | C# SMM/BIMAS | published worked examples / IEDB (caller matrix) | exact on anchor cases |
| 247 | ☐ | IMMUNE-NUSVR-001 | Oncology | C# ν-SVR (SMO) | scikit-learn NuSVR | coefficients < 2e-3 |
| 248 | ☐ | META-CHECKM-001 | Metagenomics | C# CheckM formula | CheckM markerSets.py | completeness/contamination exact on synthetic bin |
| 249 | ☐ | META-TETRA-001 | Metagenomics | C# TETRA z-score | TETRA reference (Teeling) | z-vector ±1e-6 |
| 250 | ☐ | SPLICE-MAXENT3-001 | Splicing | C# score3 | MaxEntScan score3.pl | exact score |
| 251 | ☐ | SPLICE-MAXENT5-001 | Splicing | C# score5 | MaxEntScan score5.pl | exact score |
| 252 | ☐ | MIRNA-CONTEXT-001 | MiRNA | C# context++ | targetscan_70_context_scores.pl | computable subset byte-exact |
| 253 | ☐ | MIRNA-PCT-001 | MiRNA | C# PCT | Friedman 2009 logistic worked example | PCT within tolerance |
| 254 | ☐ | MIRNA-CLASSIFY-001 | MiRNA | C# classifier | held-out miRBase vs shuffled (AUC) | AUC ≈ 1.0 on held-out set |
| 255 | ☐ | MIRNA-CLEAVAGE-001 | MiRNA | C# cleavage rules | miRBase mature coordinates | mature 5'/3' exact |
| 256 | ☐ | REP-APPROX-001 | Repeats | C# approximate TRF | TRF (Benson) on benchmark repeats | consensus + match/indel% agree |
| 257 | ☐ | CHROM-ALPHASAT-001 | Chromosome | C# alpha-satellite detector | known centromeric reference arrays | period + CENP-B positions agree |
| 258 | ☐ | CHROM-HOR-001 | Chromosome | C# HOR detector | known HOR arrays (D-region) | period + copy number agree |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 255 |
| ☑ Complete | 4 |
| ☐ Not started | 252 |
| High-value pairs (ALT/BRUTE feasible) | ~25 |
| Medium-value pairs (REF comparison) | ~35 |
| Lower priority (DUAL re-impl needed) | ~26 |
