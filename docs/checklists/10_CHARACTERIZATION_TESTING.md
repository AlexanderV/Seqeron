# Checklist 10: Characterization Testing

**Priority:** P3  
**Date:** 2026-03-19  
**Total algorithms:** 258

---

## Description

Characterization tests (golden master tests) фіксують поточну поведінку системи «as-is» перед рефакторингом. Не перевіряють коректність — перевіряють незмінність поведінки. Використовуються on-demand перед небезпечними рефакторингами.

**Поточне покриття:** 0. Фактично snapshot-тести в `Snapshots/` виконують схожу роль, але characterization tests специфічні для рефакторингу.

**Коли застосовувати:**
- Перед заміною алгоритму (нова імплементація)
- Перед оптимізацією (Span-based, SIMD, etc.)
- Перед виносом коду в окремий модуль
- Перед зміною API (параметри, типи значень, що повертаються)

**Процес:**
1. Генерація набору входів (corner cases + typical cases)
2. Запис поточних виходів у golden master
3. Рефакторинг коду
4. Запуск тестів — будь-яке розходження = fail
5. Ревʼю diff: intentional change → approve, regression → fix

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | What to Capture | Refactoring Trigger |
|---|--------|-----------|------|----------------|-------------------|
| 1 | ☐ | SEQ-GC-001 | Composition | GC% для 10+ edge-case seqs | Span-based optimization |
| 2 | ☐ | SEQ-COMP-001 | Composition | Complement для DNA/RNA/ambiguous | Switch → lookup table |
| 3 | ☐ | SEQ-REVCOMP-001 | Composition | RevComp для різних довжин | Algorithm optimization |
| 4 | ☐ | SEQ-VALID-001 | Composition | Validation results для valid/invalid | Regex → Span loop |
| 5 | ☐ | SEQ-COMPLEX-001 | Composition | Complexity scores | Algorithm change |
| 6 | ☐ | SEQ-ENTROPY-001 | Composition | Entropy для known distributions | Math library change |
| 7 | ☐ | SEQ-GCSKEW-001 | Composition | Skew arrays для test seqs | Optimization |
| 8 | ☐ | PAT-EXACT-001 | Matching | Positions для multi-pattern | Algorithm swap (KMP → Aho-Corasick) |
| 9 | ☐ | PAT-APPROX-001 | Matching | Hamming matches для threshold set | Distance function rewrite |
| 10 | ☐ | PAT-APPROX-002 | Matching | Edit distance matches для various maxDist | DP optimization |
| 11 | ☐ | PAT-IUPAC-001 | Matching | IUPAC expansion results | IupacHelper refactor |
| 12 | ☐ | PAT-PWM-001 | Matching | PWM scores для known matrices | Scoring optimization |
| 13 | ☐ | REP-STR-001 | Repeats | Microsatellites для test genomes | RepeatFinder rewrite |
| 14 | ☐ | REP-TANDEM-001 | Repeats | Tandem repeats повний output | Algorithm optimization |
| 15 | ☐ | REP-INV-001 | Repeats | Inverted repeat arms | Algorithm optimization |
| 16 | ☐ | REP-DIRECT-001 | Repeats | Direct repeat positions | Algorithm optimization |
| 17 | ☐ | REP-PALIN-001 | Repeats | Palindrome output | Algorithm optimization |
| 18 | ☐ | CRISPR-PAM-001 | MolTools | PAM site positions для SpCas9/Cas12a | Multi-PAM refactor |
| 19 | ☐ | CRISPR-GUIDE-001 | MolTools | Guide design output | Scoring model change |
| 20 | ☐ | CRISPR-OFF-001 | MolTools | Off-target results | Index structure change |
| 21 | ☐ | PRIMER-TM-001 | MolTools | Tm values для test primers | Thermodynamic model update |
| 22 | ☐ | PRIMER-DESIGN-001 | MolTools | Primer candidates output | Filter criteria change |
| 23 | ☐ | PRIMER-STRUCT-001 | MolTools | Structure analysis (hairpin, dimer) | Algorithm update |
| 24 | ☐ | PROBE-DESIGN-001 | MolTools | Probe candidates | Filter criteria change |
| 25 | ☐ | PROBE-VALID-001 | MolTools | Validation pass/fail + reasons | Validation rule change |
| 26 | ☐ | RESTR-FIND-001 | MolTools | Site positions для known enzymes | Pattern engine change |
| 27 | ☐ | RESTR-DIGEST-001 | MolTools | Fragment sizes + positions | Circular handling change |
| 28 | ☐ | ANNOT-ORF-001 | Annotation | ORF positions, frames, lengths | MinLength change |
| 29 | ☐ | ANNOT-GENE-001 | Annotation | Gene predictions | Scoring model change |
| 30 | ☐ | ANNOT-PROM-001 | Annotation | Promoter positions + scores | Weight matrix update |
| 31 | ☐ | ANNOT-GFF-001 | Annotation | GFF3 serialized strings | Format spec compliance fix |
| 32 | ☐ | KMER-COUNT-001 | K-mer | K-mer count tables для various k | HashMap → sorted array |
| 33 | ☐ | KMER-FREQ-001 | K-mer | Frequency distributions | Counting optimization |
| 34 | ☐ | KMER-FIND-001 | K-mer | Top k-mers | Sorting algorithm change |
| 35 | ☐ | ALIGN-GLOBAL-001 | Alignment | Aligned seqs + score для known pairs | DP matrix optimization |
| 36 | ☐ | ALIGN-LOCAL-001 | Alignment | Aligned seqs + score + positions | Traceback optimization |
| 37 | ☐ | ALIGN-SEMI-001 | Alignment | Semi-global output | Boundary condition change |
| 38 | ☐ | ALIGN-MULTI-001 | Alignment | MSA output | Guide tree algorithm change |
| 39 | ☐ | PHYLO-DIST-001 | Phylogenetic | Distance matrix values | Model change (JC → K2P) |
| 40 | ☐ | PHYLO-TREE-001 | Phylogenetic | Tree topology + branch lengths | UPGMA → NJ switch |
| 41 | ☐ | PHYLO-NEWICK-001 | Phylogenetic | Newick string serialization | Format compliance fix |
| 42 | ☐ | PHYLO-COMP-001 | Phylogenetic | RF distance values | Bipartition algorithm change |
| 43 | ☐ | POP-FREQ-001 | PopGen | Allele frequencies | Counting optimization |
| 44 | ☐ | POP-DIV-001 | PopGen | π, θ, Tajima's D | Formula correction |
| 45 | ☐ | POP-HW-001 | PopGen | Chi² values + p-values | Statistical method change |
| 46 | ☐ | POP-FST-001 | PopGen | Fst values | Weir-Cockerham update |
| 47 | ☐ | POP-LD-001 | PopGen | D', r² values | EM algorithm optimization |
| 48 | ☐ | CHROM-TELO-001 | Chromosome | Telomere regions | Pattern update |
| 49 | ☐ | CHROM-CENT-001 | Chromosome | Centromere position estimates | Window algorithm change |
| 50 | ☐ | CHROM-KARYO-001 | Chromosome | Karyotype classifications | Classification rule change |
| 51 | ☐ | CHROM-ANEU-001 | Chromosome | Copy number calls | Threshold change |
| 52 | ☐ | CHROM-SYNT-001 | Chromosome | Synteny block coordinates | Merge algorithm change |
| 53 | ☐ | META-CLASS-001 | Metagenomics | Taxonomic assignments | Database/k-mer size change |
| 54 | ☐ | META-PROF-001 | Metagenomics | Profile abundance values | Normalization change |
| 55 | ☐ | META-ALPHA-001 | Metagenomics | Shannon, Simpson values | Formula change |
| 56 | ☐ | META-BETA-001 | Metagenomics | Distance matrix values | Metric change |
| 57 | ☐ | META-BIN-001 | Metagenomics | Bin assignments | Clustering algorithm change |
| 58 | ☐ | CODON-OPT-001 | Codon | Optimized sequences | Usage table update |
| 59 | ☐ | CODON-CAI-001 | Codon | CAI scores | Reference table change |
| 60 | ☐ | CODON-RARE-001 | Codon | Rare codon lists | Threshold change |
| 61 | ☐ | CODON-USAGE-001 | Codon | Usage tables | Counting optimization |
| 62 | ☐ | TRANS-CODON-001 | Translation | Genetic code mappings | Alternative table support |
| 63 | ☐ | TRANS-PROT-001 | Translation | Translation results | Reading frame handling change |
| 64 | ☐ | PARSE-FASTA-001 | FileIO | Parsed records для sample files | Parser rewrite (streaming) |
| 65 | ☐ | PARSE-FASTQ-001 | FileIO | Parsed records + quality | Parser rewrite (streaming) |
| 66 | ☐ | PARSE-BED-001 | FileIO | Parsed BED regions | Parser rewrite |
| 67 | ☐ | PARSE-VCF-001 | FileIO | Parsed variants | VCF 4.3 compliance |
| 68 | ☐ | PARSE-GFF-001 | FileIO | Parsed features | GFF3 spec compliance |
| 69 | ☐ | PARSE-GENBANK-001 | FileIO | Parsed GenBank records | Feature location parser rewrite |
| 70 | ☐ | PARSE-EMBL-001 | FileIO | Parsed EMBL records | Format spec compliance |
| 71 | ☐ | RNA-STRUCT-001 | RnaStructure | Structure predictions | Algorithm upgrade (MFE) |
| 72 | ☐ | RNA-STEMLOOP-001 | RnaStructure | Stem-loop detection results | Detection criteria change |
| 73 | ☐ | RNA-ENERGY-001 | RnaStructure | MFE values | Energy parameter update |
| 74 | ☐ | MIRNA-SEED-001 | MiRNA | Seed extraction results | Seed definition change |
| 75 | ☐ | MIRNA-TARGET-001 | MiRNA | Target predictions + scores | Scoring model update |
| 76 | ☐ | MIRNA-PRECURSOR-001 | MiRNA | Pre-miRNA analysis | Classification criteria change |
| 77 | ☐ | SPLICE-DONOR-001 | Splicing | Donor site scores | PWM update |
| 78 | ☐ | SPLICE-ACCEPTOR-001 | Splicing | Acceptor site scores | PWM update |
| 79 | ☐ | SPLICE-PREDICT-001 | Splicing | Gene structure predictions | Prediction model change |
| 80 | ☐ | DISORDER-PRED-001 | ProteinPred | Disorder scores | Propensity scale update |
| 81 | ☐ | DISORDER-REGION-001 | ProteinPred | Disordered regions | Window/threshold change |
| 82 | ☐ | PROTMOTIF-FIND-001 | ProteinMotif | Motif positions + patterns | Pattern matching engine change |
| 83 | ☐ | PROTMOTIF-PROSITE-001 | ProteinMotif | PROSITE matches | Pattern database update |
| 84 | ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | Domain predictions | Profile/HMM update |
| 85 | ☐ | EPIGEN-CPG-001 | Epigenetics | CpG island regions | O/E criteria change |
| 86 | ☐ | ONCO-IMMUNE-001 | Oncology | Immune scores full output | Gene set / algorithm change |
| 87 | ☐ | ONCO-SOMATIC-001 | Oncology | Somatic call set golden master | Caller logic / threshold change |
| 88 | ☐ | ONCO-VAF-001 | Oncology | VAF outputs for fixed pileups | Counting / rounding change |
| 89 | ☐ | ONCO-DRIVER-001 | Oncology | Driver list + scores | Scoring / gene-list change |
| 90 | ☐ | ONCO-ARTIFACT-001 | Oncology | Surviving variant set | Filter threshold change |
| 91 | ☐ | ONCO-ANNOT-001 | Oncology | Annotation records | Consequence rules / DB change |
| 92 | ☐ | ONCO-TMB-001 | Oncology | TMB value + class | Panel-size / counting change |
| 93 | ☐ | ONCO-MSI-001 | Oncology | MSI status + loci | Threshold / locus-panel change |
| 94 | ☐ | ONCO-HRD-001 | Oncology | HRD score + components | Component formula change |
| 95 | ☐ | ONCO-LOH-001 | Oncology | LOH regions | BAF threshold change |
| 96 | ☐ | ONCO-SIG-001 | Oncology | 96-channel SBS matrix | Context-assignment change |
| 97 | ☐ | ONCO-SIG-002 | Oncology | Exposure vector | Fit algorithm change |
| 98 | ☐ | ONCO-SIG-003 | Oncology | Bootstrap CIs | Resampling / seed change |
| 99 | ☐ | ONCO-SIG-004 | Oncology | Process classification | Classification rule change |
| 100 | ☐ | ONCO-FUSION-001 | Oncology | Fusion call set | Detection logic change |
| 101 | ☐ | ONCO-FUSION-002 | Oncology | Known-fusion matches | DB change |
| 102 | ☐ | ONCO-FUSION-003 | Oncology | Breakpoint / frame output | Frame logic change |
| 103 | ☐ | ONCO-CNA-001 | Oncology | CN classification | log2→CN mapping change |
| 104 | ☐ | ONCO-CNA-002 | Oncology | Focal amplification calls | Threshold change |
| 105 | ☐ | ONCO-CNA-003 | Oncology | Homozygous deletion calls | Threshold change |
| 106 | ☐ | ONCO-PURITY-001 | Oncology | Purity estimate | Estimator change |
| 107 | ☐ | ONCO-PLOIDY-001 | Oncology | Ploidy estimate | Estimator change |
| 108 | ☐ | ONCO-CLONAL-001 | Oncology | Clonality classification | CCF threshold change |
| 109 | ☐ | ONCO-NEO-001 | Oncology | Neoantigen peptides | Window / length change |
| 110 | ☐ | ONCO-MHC-001 | Oncology | Binding classification | Affinity model change |
| 111 | ☐ | ONCO-CTDNA-001 | Oncology | ctDNA fraction + output | Estimator change |
| 112 | ☐ | ONCO-MRD-001 | Oncology | MRD status | Detection threshold change |
| 113 | ☐ | ONCO-CHIP-001 | Oncology | CHIP flags | Gene-list / VAF band change |
| 114 | ☐ | ONCO-PHYLO-001 | Oncology | Tumor phylogeny (Newick) | Reconstruction algorithm change |
| 115 | ☐ | ONCO-CCF-001 | Oncology | CCF estimates | Formula change |
| 116 | ☐ | ONCO-HETERO-001 | Oncology | MATH / heterogeneity output | Metric change |
| 117 | ☐ | ONCO-HLA-001 | Oncology | HLA alleles | Typing algorithm change |
| 118 | ☐ | ONCO-ACTION-001 | Oncology | Actionability tiers | Evidence DB / tier change |
| 119 | ☐ | ONCO-SV-001 | Oncology | Complex SV classification | Pattern rules change |
| 120 | ☐ | ONCO-EXPR-001 | Oncology | Outlier gene list | z-score threshold change |
| 121 | ☐ | SEQ-COMPOSITION-001 | Statistics | Composition vector | Counting / normalization change |
| 122 | ☐ | SEQ-DINUC-001 | Statistics | Dinucleotide table | Counting change |
| 123 | ☐ | SEQ-HYDRO-001 | Statistics | Hydrophobicity profile | Scale / window change |
| 124 | ☐ | SEQ-MW-001 | Statistics | Molecular weight | Residue-mass table change |
| 125 | ☐ | SEQ-PI-001 | Statistics | Isoelectric point | pKa table change |
| 126 | ☐ | SEQ-SECSTRUCT-001 | Statistics | Secondary structure assignment | Propensity table change |
| 127 | ☐ | SEQ-STATS-001 | Statistics | Composition statistics | Calculation change |
| 128 | ☐ | SEQ-SUMMARY-001 | Statistics | Sequence summary | Summary-field change |
| 129 | ☐ | SEQ-THERMO-001 | Statistics | Thermodynamic parameters | NN-parameter change |
| 130 | ☐ | SEQ-TM-001 | Statistics | Tm value | Formula change |
| 131 | ☐ | COMPGEN-ANI-001 | Comparative | ANI matrix | Identity / k-mer method change |
| 132 | ☐ | COMPGEN-CLUSTER-001 | Comparative | Conserved clusters | Threshold change |
| 133 | ☐ | COMPGEN-COMPARE-001 | Comparative | Comparison report | Metric change |
| 134 | ☐ | COMPGEN-DOTPLOT-001 | Comparative | Dot-plot coordinates | Word size change |
| 135 | ☐ | COMPGEN-ORTHO-001 | Comparative | Ortholog pairs | Detection algorithm change |
| 136 | ☐ | COMPGEN-RBH-001 | Comparative | RBH pairs | Scoring change |
| 137 | ☐ | COMPGEN-REARR-001 | Comparative | Rearrangement output | Breakpoint logic change |
| 138 | ☐ | COMPGEN-REVERSAL-001 | Comparative | Reversal distance | Algorithm change |
| 139 | ☐ | COMPGEN-SYNTENY-001 | Comparative | Syntenic blocks | minBlockSize / chaining change |
| 140 | ☐ | ASSEMBLY-CONSENSUS-001 | Assembly | Consensus sequence | Voting rule change |
| 141 | ☐ | ASSEMBLY-CORRECT-001 | Assembly | Corrected reads | Correction threshold change |
| 142 | ☐ | ASSEMBLY-COVER-001 | Assembly | Coverage profile | Calculation change |
| 143 | ☐ | ASSEMBLY-DBG-001 | Assembly | De Bruijn contigs | k / graph-cleaning change |
| 144 | ☐ | ASSEMBLY-MERGE-001 | Assembly | Merged contigs | Overlap rule change |
| 145 | ☐ | ASSEMBLY-OLC-001 | Assembly | OLC contigs | minOverlap change |
| 146 | ☐ | ASSEMBLY-SCAFFOLD-001 | Assembly | Scaffold layout | Linking rule change |
| 147 | ☐ | ASSEMBLY-STATS-001 | Assembly | Assembly statistics | Metric-definition change |
| 148 | ☐ | ASSEMBLY-TRIM-001 | Assembly | Trimmed reads | Quality cutoff / window change |
| 149 | ☐ | RNA-DOTBRACKET-001 | RnaStructure | Parsed pairs | Notation / parser change |
| 150 | ☐ | RNA-HAIRPIN-001 | RnaStructure | Hairpin energy | Energy model change |
| 151 | ☐ | RNA-INVERT-001 | RnaStructure | Inverted repeats | Detection logic change |
| 152 | ☐ | RNA-MFE-001 | RnaStructure | MFE structure | Energy parameter change |
| 153 | ☐ | RNA-PAIR-001 | RnaStructure | Pairing table | Pairing rule change |
| 154 | ☐ | RNA-PARTITION-001 | RnaStructure | Pair probabilities | Algorithm change |
| 155 | ☐ | RNA-PSEUDOKNOT-001 | RnaStructure | Pseudoknots | Detection logic change |
| 156 | ☐ | KMER-ASYNC-001 | K-mer | Async counts | Counting change |
| 157 | ☐ | KMER-BOTH-001 | K-mer | Both-strand counts | Counting change |
| 158 | ☐ | KMER-DIST-001 | K-mer | Distance value | Metric change |
| 159 | ☐ | KMER-GENERATE-001 | K-mer | Generated set | Generation change |
| 160 | ☐ | KMER-POSITIONS-001 | K-mer | Positions | Indexing change |
| 161 | ☐ | KMER-STATS-001 | K-mer | Statistics | Calculation change |
| 162 | ☐ | KMER-UNIQUE-001 | K-mer | Unique k-mers | Threshold change |
| 163 | ☐ | PROTMOTIF-CC-001 | ProteinMotif | Coiled-coil scores | Scoring change |
| 164 | ☐ | PROTMOTIF-COMMON-001 | ProteinMotif | Common motifs | Algorithm change |
| 165 | ☐ | PROTMOTIF-LC-001 | ProteinMotif | Low-complexity regions | Threshold change |
| 166 | ☐ | PROTMOTIF-PATTERN-001 | ProteinMotif | Pattern matches | Pattern engine change |
| 167 | ☐ | PROTMOTIF-SP-001 | ProteinMotif | Signal-peptide call | Model change |
| 168 | ☐ | PROTMOTIF-TM-001 | ProteinMotif | TM helices | Threshold change |
| 169 | ☐ | MOTIF-CONS-001 | Matching | Consensus | Voting change |
| 170 | ☐ | MOTIF-DISCOVER-001 | Matching | Discovered motifs | Algorithm change |
| 171 | ☐ | MOTIF-GENERATE-001 | Matching | Consensus | Tie-break change |
| 172 | ☐ | MOTIF-REGULATORY-001 | Matching | Regulatory elements | Reference-set change |
| 173 | ☐ | MOTIF-SHARED-001 | Matching | Shared motifs | Algorithm change |
| 174 | ☐ | PAT-APPROX-003 | Matching | Best match | Distance metric change |
| 175 | ☐ | GENOMIC-COMMON-001 | Analysis | Common region | Algorithm change |
| 176 | ☐ | GENOMIC-MOTIFS-001 | Analysis | Motif hits | Reference change |
| 177 | ☐ | GENOMIC-ORF-001 | Analysis | ORF list | ORF rule change |
| 178 | ☐ | GENOMIC-REPEAT-001 | Analysis | Repeats | Parameter change |
| 179 | ☐ | GENOMIC-SIMILARITY-001 | Analysis | Similarity | Metric change |
| 180 | ☐ | GENOMIC-TANDEM-001 | Analysis | Tandem repeats | Parameter change |
| 181 | ☐ | EPIGEN-AGE-001 | Epigenetics | Age value | Clock model change |
| 182 | ☐ | EPIGEN-BISULF-001 | Epigenetics | Converted sequence | Conversion rule change |
| 183 | ☐ | EPIGEN-CHROM-001 | Epigenetics | Chromatin states | Model change |
| 184 | ☐ | EPIGEN-DMR-001 | Epigenetics | DMRs | Threshold change |
| 185 | ☐ | EPIGEN-METHYL-001 | Epigenetics | Methylation | Calculation change |
| 186 | ☐ | VARIANT-ANNOT-001 | Variants | Annotations | Rule / DB change |
| 187 | ☐ | VARIANT-CALL-001 | Variants | Calls | Caller threshold change |
| 188 | ☐ | VARIANT-INDEL-001 | Variants | Indels | Detection change |
| 189 | ☐ | VARIANT-SNP-001 | Variants | SNPs | Detection change |
| 190 | ☐ | PANGEN-CLUSTER-001 | PanGenome | Clusters | Identity change |
| 191 | ☐ | PANGEN-CORE-001 | PanGenome | Core / accessory | Definition change |
| 192 | ☐ | PANGEN-HEAP-001 | PanGenome | Heaps fit | Regression change |
| 193 | ☐ | PANGEN-MARKER-001 | PanGenome | Markers | Selection change |
| 194 | ☐ | META-FUNC-001 | Metagenomics | Functions | DB change |
| 195 | ☐ | META-PATHWAY-001 | Metagenomics | Enrichment | Statistic change |
| 196 | ☐ | META-RESIST-001 | Metagenomics | Resistance hits | DB change |
| 197 | ☐ | META-TAXA-001 | Metagenomics | Significant taxa | Statistic change |
| 198 | ☐ | TRANS-DIFF-001 | Transcriptome | DE table | Model change |
| 199 | ☐ | TRANS-EXPR-001 | Transcriptome | Expression | Normalization change |
| 200 | ☐ | TRANS-SPLICE-001 | Transcriptome | Splicing | Algorithm change |
| 201 | ☐ | SV-BREAKPOINT-001 | StructuralVar | Breakpoints | Detection change |
| 202 | ☐ | SV-CNV-001 | StructuralVar | CNV calls | Threshold change |
| 203 | ☐ | SV-DETECT-001 | StructuralVar | SV calls | Detection change |
| 204 | ☐ | DISORDER-LC-001 | ProteinPred | Low-complexity regions | Threshold change |
| 205 | ☐ | DISORDER-MORF-001 | ProteinPred | MoRFs | Model change |
| 206 | ☐ | DISORDER-PROPENSITY-001 | ProteinPred | Propensity | Scale change |
| 207 | ☐ | POP-ANCESTRY-001 | PopGen | Ancestry | Model change |
| 208 | ☐ | POP-ROH-001 | PopGen | ROH segments | Parameter change |
| 209 | ☐ | POP-SELECT-001 | PopGen | Selection signal | Statistic change |
| 210 | ☐ | SEQ-ATSKEW-001 | Composition | AT-skew array | Calculation change |
| 211 | ☐ | SEQ-REPLICATION-001 | Composition | Origin | Method change |
| 212 | ☐ | SEQ-RNACOMP-001 | Composition | RNA complement | Mapping change |
| 213 | ☐ | CODON-ENC-001 | Codon | ENC | Formula change |
| 214 | ☐ | CODON-RSCU-001 | Codon | RSCU | Calculation change |
| 215 | ☐ | CODON-STATS-001 | Codon | Statistics | Calculation change |
| 216 | ☐ | ANNOT-CODING-001 | Annotation | Coding score | Model change |
| 217 | ☐ | ANNOT-CODONUSAGE-001 | Annotation | Codon usage | Calculation change |
| 218 | ☐ | ANNOT-REPEAT-001 | Annotation | Repetitive elements | Parameter change |
| 219 | ☐ | QUALITY-PHRED-001 | Quality | Phred scores | Offset change |
| 220 | ☐ | QUALITY-STATS-001 | Quality | Quality stats | Calculation change |
| 221 | ☐ | PHYLO-BOOT-001 | Phylogenetic | Support values | Resampling change |
| 222 | ☐ | PHYLO-STATS-001 | Phylogenetic | Tree stats | Metric change |
| 223 | ☐ | TRANS-SIXFRAME-001 | Translation | Six frames | Table change |
| 224 | ☐ | RESTR-FILTER-001 | MolTools | Filtered sites | Criteria change |
| 225 | ☐ | MIRNA-PAIR-001 | MiRNA | Alignment | Pairing rule change |
| 226 | ☐ | ALIGN-STATS-001 | Alignment | Statistics | Calculation change |
| 227 | ☐ | SEQ-CODON-FREQ-001 | Statistics | Codon frequencies for 10+ CDS | Counting optimization |
| 228 | ☐ | SEQ-COMPLEX-COMPRESS-001 | Complexity | Compression ratios for edge seqs | Compression algorithm change |
| 229 | ☐ | SEQ-COMPLEX-DUST-001 | Complexity | DUST scores for edge seqs | Scoring change |
| 230 | ☐ | SEQ-COMPLEX-KMER-001 | Complexity | K-mer entropy for edge seqs | Entropy formula change |
| 231 | ☐ | SEQ-COMPLEX-WINDOW-001 | Complexity | Windowed complexity profile | Window-stepping optimization |
| 232 | ☐ | SEQ-ENTROPY-PROFILE-001 | Statistics | Entropy profile for edge seqs | Window-stepping optimization |
| 233 | ☐ | SEQ-GC-ANALYSIS-001 | Composition | GC analysis for edge seqs | Span-based optimization |
| 234 | ☐ | SEQ-GC-PROFILE-001 | Statistics | GC profile for edge seqs | Window-stepping optimization |
| 235 | ☐ | ONCO-ASCAT-001 | Oncology | purity/ploidy for planted-truth profiles | grid optimization |
| 236 | ☐ | RNA-PKPREDICT-001 | Analysis | PK structures for known H-type RNAs | O(n³) scan |
| 237 | ☐ | RNA-PKRECURSIVE-001 | Analysis | nested-knot structures for designed RNAs | recursive DP |
| 238 | ☑ | RNA-ACCESS-001 | RnaStructure | unpaired-probability vectors for reference RNAs | partition-function model change |
| 239 | ☑ | PROTMOTIF-HMM-001 | ProteinMotif | domain hits (coords, bit score, E-value) for SH3/PDZ/WD40 | HMMER pipeline parity change |
| 240 | ☑ | PRIMER-NNTM-001 | MolTools | Tm for reference oligos at fixed salt | NN parameter set change |
| 241 | ☑ | PRIMER-HAIRPIN-001 | MolTools | best-hairpin ΔG/Tm for reference primers | loop-init table change |
| 242 | ☑ | PRIMER-DIMER-001 | MolTools | best-dimer ΔG/Tm for primer pairs | ntthal alignment extension change |
| 243 | ☑ | PROBE-LNATM-001 | MolTools | LNA Tm + MGB verdict for reference probes | LNA increment table change |
| 244 | ☑ | PROBE-EVALUE-001 | MolTools | E-value/bit score for reference HSPs | scoring-system change |
| 245 | ☑ | MHC-NN-001 | Oncology | IC50 predictions for benchmark peptide/allele pairs | weight-pack update |
| 246 | ☑ | MHC-MATRIX-001 | Oncology | IC50/half-life for reference peptides | matrix reload |
| 247 | ☑ | IMMUNE-NUSVR-001 | Oncology | cell-fraction estimates for reference mixtures | signature-matrix change |
| 248 | ☑ | META-CHECKM-001 | Metagenomics | completeness/contamination for synthetic + real bins | marker-set update |
| 249 | ☐ | META-TETRA-001 | Metagenomics | z-vector + pairwise correlations for reference contigs | expected-frequency model change |
| 250 | ☐ | SPLICE-MAXENT3-001 | Splicing | acceptor MaxEnt scores for reference 3' sites | maxent table change |
| 251 | ☐ | SPLICE-MAXENT5-001 | Splicing | donor MaxEnt scores for reference 5' sites | maxent table change |
| 252 | ☐ | MIRNA-CONTEXT-001 | MiRNA | context++ scores for reference sites | coefficient table change |
| 253 | ☐ | MIRNA-PCT-001 | MiRNA | PCT for reference conserved sites | sigmoid-parameter change |
| 254 | ☐ | MIRNA-CLASSIFY-001 | MiRNA | precursor probabilities for miRBase positives | classifier retraining |
| 255 | ☐ | MIRNA-CLEAVAGE-001 | MiRNA | Drosha/Dicer cut sites for reference precursors | measuring-rule change |
| 256 | ☐ | REP-APPROX-001 | Repeats | consensus + match/indel% for reference repeats | TRF scoring change |
| 257 | ☐ | CHROM-ALPHASAT-001 | Chromosome | monomer period + CENP-B boxes for reference arrays | motif/threshold change |
| 258 | ☐ | CHROM-HOR-001 | Chromosome | HOR period/copy#/identity for reference arrays | identity-threshold change |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 255 |
| ☑ Complete | 9 |
| ☐ Not started | 247 |
| Applies on-demand (before refactoring) | All 234 |
| High refactoring risk (complex algorithms) | ~20 (Alignment, Phylogenetic, RNA, Annotation) |
| Medium refactoring risk | ~40 |
| Lower refactoring risk (simple calculation) | ~26 |
