# Checklist 09: Characterization Testing

**Priority:** P3  
**Date:** 2026-03-19  
**Total algorithms:** 86

---

## Description

Characterization tests (golden master tests) фиксируют текущее поведение системы «as-is» перед рефакторингом. Не проверяют корректность — проверяют неизменность поведения. Используются on-demand перед опасными рефакторами.

**Текущее покрытие:** 0. Фактически snapshot-тесты в `Snapshots/` выполняют похожую роль, но characterization tests специфичны для рефакторинга.

**Когда применять:**
- Перед заменой алгоритма (новая имплементация)
- Перед оптимизацией (Span-based, SIMD, etc.)
- Перед выносом кода в отдельный модуль
- Перед изменением API (параметры, типы возвращаемых значений)

**Процесс:**
1. Генерация набора входов (corner cases + typical cases)
2. Запись текущих выходов в golden master
3. Рефакторинг кода
4. Запуск тестов — любое расхождение = fail
5. Ревью diff: intentional change → approve, regression → fix

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | What to Capture | Refactoring Trigger |
|---|--------|-----------|------|----------------|-------------------|
| 1 | ☐ | SEQ-GC-001 | Composition | GC% для 10+ edge-case seqs | Span-based optimization |
| 2 | ☐ | SEQ-COMP-001 | Composition | Complement для DNA/RNA/ambiguous | Switch → lookup table |
| 3 | ☐ | SEQ-REVCOMP-001 | Composition | RevComp для различных длин | Algorithm optimization |
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
| 14 | ☐ | REP-TANDEM-001 | Repeats | Tandem repeats полный output | Algorithm optimization |
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

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 120 |
| ☑ Complete | 0 |
| ☐ Not started | 120 |
| Applies on-demand (before refactoring) | All 120 |
| High refactoring risk (complex algorithms) | ~20 (Alignment, Phylogenetic, RNA, Annotation) |
| Medium refactoring risk | ~40 |
| Lower refactoring risk (simple calculation) | ~26 |
