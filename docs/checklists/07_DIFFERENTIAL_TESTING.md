# Checklist 07: Differential Testing

**Priority:** P2
**Date:** 2026-03-19
**Total algorithms:** 86

---

## Description

Differential testing сравнивает выходы двух независимых реализаций одного алгоритма на идентичных входах. Обнаруживает тонкие ошибки реализации. В биоинформатике существуют референсные реализации (Biopython, EMBOSS, etc.) и альтернативные алгоритмические стратегии.

**Текущее покрытие:** `SafeVsUnsafeDifferentialTests` существует только для SuffixTree. Ноль для Seqeron.Genomics.

**Стратегии:**
- **ALT** = Альтернативный алгоритм (другой подход к той же задаче)
- **BRUTE** = Brute-force реализация для малых входов
- **REF** = Сравнение с вычисленным вручную результатом / референсной библиотекой
- **DUAL** = Две реализации внутри проекта (разные методы одного класса)

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

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 86 |
| ☑ Complete | 0 |
| ☐ Not started | 86 |
| High-value pairs (ALT/BRUTE feasible) | ~25 |
| Medium-value pairs (REF comparison) | ~35 |
| Lower priority (DUAL re-impl needed) | ~26 |
