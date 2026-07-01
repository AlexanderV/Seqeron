# Checklist 07: Architecture Testing (ArchUnitNET)

**Priority:** P2  
**Framework:** ArchUnitNET  
**Date:** 2026-06-27  
**Total algorithms:** 258

---

## Description

Architecture тести запобігають архітектурному дрейфу, перевіряючи залежності між модулями, naming conventions та структурні правила на рівні IL. Правила застосовуються до модулів (проєктів), а не до окремих алгоритмів.

**Поточне покриття:** `Architecture/ArchitectureTests.cs` — 22 правила (всі ☑): межі шарів Core/IO до всіх 13 модулів, відсутність циклічних залежностей між модулями, заборона System.IO у Core, незмінність Result/DTO-типів, а також placement/naming-інваріанти (парсери лише в IO, Core без алгоритмічних класів, namespace=assembly). Модульні правила не залежать від кількості алгоритмів — нові алгоритми автоматично покриваються правилами свого модуля.

**Модульна структура:**
- `Seqeron.Genomics.Core` — базові типи, послідовності, розширення
- `Seqeron.Genomics.Analysis` — аналізатори (K-mer, Repeat, Motif, Disorder, GC-skew, Complexity, Comparative, etc.)
- `Seqeron.Genomics.Alignment` — вирівнювання (NW, SW, MSA) та збірка (assembly: DBG/OLC/scaffold)
- `Seqeron.Genomics.IO` — парсери файлів (FASTA, FASTQ, BED, VCF, etc.) та якість (Phred)
- `Seqeron.Genomics.Annotation` — анотація геному (ORF, gene, promoter), miRNA, splicing, epigenetics, variants, SV, transcriptome
- `Seqeron.Genomics.MolTools` — молекулярні інструменти (CRISPR, primer, probe, restriction, codon optimization)
- `Seqeron.Genomics.Phylogenetics` — філогенетика
- `Seqeron.Genomics.Population` — популяційна генетика
- `Seqeron.Genomics.Chromosome` — хромосомний аналіз
- `Seqeron.Genomics.Metagenomics` — метагеноміка та пангеном
- `Seqeron.Genomics.Oncology` — онкологія
- `Seqeron.Genomics.Infrastructure` — інфраструктура
- `Seqeron.Genomics.Reports` — генерація звітів

---

## Checklist — Module Dependency Rules

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Rule | Description |
|---|--------|------|-------------|
| 1 | ☑ | Core !→ Analysis | Core не залежить від Analysis |
| 2 | ☑ | Core !→ IO | Core не залежить від IO |
| 3 | ☑ | Core !→ Alignment | Core не залежить від Alignment |
| 4 | ☑ | IO !→ Analysis | IO не залежить від Analysis |
| 5 | ☑ | Analyzers static | Predictor/Finder класи статичні |
| 6 | ☑ | Core !→ Annotation | Core не залежить від Annotation |
| 7 | ☑ | Core !→ MolTools | Core не залежить від MolTools |
| 8 | ☑ | Core !→ Phylogenetics | Core не залежить від Phylogenetics |
| 9 | ☑ | Core !→ Population | Core не залежить від Population |
| 10 | ☑ | Core !→ Chromosome | Core не залежить від Chromosome |
| 11 | ☑ | Core !→ Metagenomics | Core не залежить від Metagenomics |
| 12 | ☑ | Core !→ Oncology | Core не залежить від Oncology |
| 13 | ☑ | Core !→ Reports | Core не залежить від Reports |
| 14 | ☑ | IO !→ Alignment | IO не залежить від Alignment |
| 15 | ☑ | IO !→ Phylogenetics | IO не залежить від Phylogenetics |
| 16 | ☑ | IO !→ Oncology | IO не залежить від Oncology |
| 17 | ☑ | No circular deps | Немає циклічних залежностей між модулями |
| 18 | ☑ | No System.IO in Core | Core не використовує System.IO напряму |
| 19 | ☑ | DTOs immutable | Result/DTO типи — records або лише getters |
| 20 | ☑ | Parsers in IO | `*Parser` класи лежать лише в модулі IO |
| 21 | ☑ | Core primitives-only | Core не містить алгоритмічних класів (`*Parser/*Predictor/*Finder/*Designer/*Analyzer/*Optimizer/*Caller/*Aligner/*Assembler`) |
| 22 | ☑ | Namespace = assembly | namespace типу збігається з його модулем (проєкт = namespace) |

---

## Checklist — Algorithm-to-Module Mapping

Кожен алгоритм верифікується модульними правилами для свого модуля. Модуль визначено за фактичним розташуванням класу-реалізації у `src/Seqeron/Algorithms/` станом на 2026-06-27.

| # | Test Unit | Area | Module | Layer Rules Applied |
|---|-----------|------|--------|-------------------|
| 1 | SEQ-GC-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 2 | SEQ-COMP-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 3 | SEQ-REVCOMP-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 4 | SEQ-VALID-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 5 | SEQ-COMPLEX-001 | Composition | Analysis | Core !→ Analysis |
| 6 | SEQ-ENTROPY-001 | Composition | Analysis | Core !→ Analysis |
| 7 | SEQ-GCSKEW-001 | Composition | Analysis | Core !→ Analysis |
| 8 | PAT-EXACT-001 | Matching | Analysis | Core !→ Analysis |
| 9 | PAT-APPROX-001 | Matching | Alignment | Core !→ Alignment, IO !→ Alignment |
| 10 | PAT-APPROX-002 | Matching | Alignment | Core !→ Alignment, IO !→ Alignment |
| 11 | PAT-IUPAC-001 | Matching | Analysis | Core !→ Analysis |
| 12 | PAT-PWM-001 | Matching | Analysis | Core !→ Analysis |
| 13 | REP-STR-001 | Repeats | Analysis | Core !→ Analysis |
| 14 | REP-TANDEM-001 | Repeats | Analysis | Core !→ Analysis |
| 15 | REP-INV-001 | Repeats | Analysis | Core !→ Analysis |
| 16 | REP-DIRECT-001 | Repeats | Analysis | Core !→ Analysis |
| 17 | REP-PALIN-001 | Repeats | Analysis | Core !→ Analysis |
| 18 | CRISPR-PAM-001 | MolTools | MolTools | Core !→ MolTools |
| 19 | CRISPR-GUIDE-001 | MolTools | MolTools | Core !→ MolTools |
| 20 | CRISPR-OFF-001 | MolTools | MolTools | Core !→ MolTools |
| 21 | PRIMER-TM-001 | MolTools | MolTools | Core !→ MolTools |
| 22 | PRIMER-DESIGN-001 | MolTools | MolTools | Core !→ MolTools |
| 23 | PRIMER-STRUCT-001 | MolTools | MolTools | Core !→ MolTools |
| 24 | PROBE-DESIGN-001 | MolTools | MolTools | Core !→ MolTools |
| 25 | PROBE-VALID-001 | MolTools | MolTools | Core !→ MolTools |
| 26 | RESTR-FIND-001 | MolTools | MolTools | Core !→ MolTools |
| 27 | RESTR-DIGEST-001 | MolTools | MolTools | Core !→ MolTools |
| 28 | ANNOT-ORF-001 | Annotation | Annotation | Core !→ Annotation |
| 29 | ANNOT-GENE-001 | Annotation | Annotation | Core !→ Annotation |
| 30 | ANNOT-PROM-001 | Annotation | Annotation | Core !→ Annotation |
| 31 | ANNOT-GFF-001 | Annotation | Annotation | Core !→ Annotation |
| 32 | KMER-COUNT-001 | K-mer | Analysis | Core !→ Analysis |
| 33 | KMER-FREQ-001 | K-mer | Analysis | Core !→ Analysis |
| 34 | KMER-FIND-001 | K-mer | Analysis | Core !→ Analysis |
| 35 | ALIGN-GLOBAL-001 | Alignment | Alignment | Core !→ Alignment, IO !→ Alignment |
| 36 | ALIGN-LOCAL-001 | Alignment | Alignment | Core !→ Alignment, IO !→ Alignment |
| 37 | ALIGN-SEMI-001 | Alignment | Alignment | Core !→ Alignment, IO !→ Alignment |
| 38 | ALIGN-MULTI-001 | Alignment | Alignment | Core !→ Alignment, IO !→ Alignment |
| 39 | PHYLO-DIST-001 | Phylogenetic | Phylogenetics | Core !→ Phylogenetics, IO !→ Phylogenetics |
| 40 | PHYLO-TREE-001 | Phylogenetic | Phylogenetics | Core !→ Phylogenetics, IO !→ Phylogenetics |
| 41 | PHYLO-NEWICK-001 | Phylogenetic | Phylogenetics | Core !→ Phylogenetics, IO !→ Phylogenetics |
| 42 | PHYLO-COMP-001 | Phylogenetic | Phylogenetics | Core !→ Phylogenetics, IO !→ Phylogenetics |
| 43 | POP-FREQ-001 | PopGen | Population | Core !→ Population |
| 44 | POP-DIV-001 | PopGen | Population | Core !→ Population |
| 45 | POP-HW-001 | PopGen | Population | Core !→ Population |
| 46 | POP-FST-001 | PopGen | Population | Core !→ Population |
| 47 | POP-LD-001 | PopGen | Population | Core !→ Population |
| 48 | CHROM-TELO-001 | Chromosome | Chromosome | Core !→ Chromosome |
| 49 | CHROM-CENT-001 | Chromosome | Chromosome | Core !→ Chromosome |
| 50 | CHROM-KARYO-001 | Chromosome | Chromosome | Core !→ Chromosome |
| 51 | CHROM-ANEU-001 | Chromosome | Chromosome | Core !→ Chromosome |
| 52 | CHROM-SYNT-001 | Chromosome | Chromosome | Core !→ Chromosome |
| 53 | META-CLASS-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 54 | META-PROF-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 55 | META-ALPHA-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 56 | META-BETA-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 57 | META-BIN-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 58 | CODON-OPT-001 | Codon | MolTools | Core !→ MolTools |
| 59 | CODON-CAI-001 | Codon | MolTools | Core !→ MolTools |
| 60 | CODON-RARE-001 | Codon | MolTools | Core !→ MolTools |
| 61 | CODON-USAGE-001 | Codon | MolTools | Core !→ MolTools |
| 62 | TRANS-CODON-001 | Translation | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 63 | TRANS-PROT-001 | Translation | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 64 | PARSE-FASTA-001 | FileIO | IO | IO !→ Analysis, IO !→ Alignment |
| 65 | PARSE-FASTQ-001 | FileIO | IO | IO !→ Analysis, IO !→ Alignment |
| 66 | PARSE-BED-001 | FileIO | IO | IO !→ Analysis, IO !→ Alignment |
| 67 | PARSE-VCF-001 | FileIO | IO | IO !→ Analysis, IO !→ Alignment |
| 68 | PARSE-GFF-001 | FileIO | IO | IO !→ Analysis, IO !→ Alignment |
| 69 | PARSE-GENBANK-001 | FileIO | IO | IO !→ Analysis, IO !→ Alignment |
| 70 | PARSE-EMBL-001 | FileIO | IO | IO !→ Analysis, IO !→ Alignment |
| 71 | RNA-STRUCT-001 | RnaStructure | Analysis | Core !→ Analysis |
| 72 | RNA-STEMLOOP-001 | RnaStructure | Analysis | Core !→ Analysis |
| 73 | RNA-ENERGY-001 | RnaStructure | Analysis | Core !→ Analysis |
| 74 | MIRNA-SEED-001 | MiRNA | Annotation | Core !→ Annotation |
| 75 | MIRNA-TARGET-001 | MiRNA | Annotation | Core !→ Annotation |
| 76 | MIRNA-PRECURSOR-001 | MiRNA | Annotation | Core !→ Annotation |
| 77 | SPLICE-DONOR-001 | Splicing | Annotation | Core !→ Annotation |
| 78 | SPLICE-ACCEPTOR-001 | Splicing | Annotation | Core !→ Annotation |
| 79 | SPLICE-PREDICT-001 | Splicing | Annotation | Core !→ Annotation |
| 80 | DISORDER-PRED-001 | ProteinPred | Analysis | Core !→ Analysis |
| 81 | DISORDER-REGION-001 | ProteinPred | Analysis | Core !→ Analysis |
| 82 | PROTMOTIF-FIND-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 83 | PROTMOTIF-PROSITE-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 84 | PROTMOTIF-DOMAIN-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 85 | EPIGEN-CPG-001 | Epigenetics | Annotation | Core !→ Annotation |
| 86 | ONCO-IMMUNE-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 87 | ONCO-SOMATIC-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 88 | ONCO-VAF-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 89 | ONCO-DRIVER-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 90 | ONCO-ARTIFACT-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 91 | ONCO-ANNOT-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 92 | ONCO-TMB-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 93 | ONCO-MSI-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 94 | ONCO-HRD-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 95 | ONCO-LOH-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 96 | ONCO-SIG-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 97 | ONCO-SIG-002 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 98 | ONCO-SIG-003 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 99 | ONCO-SIG-004 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 100 | ONCO-FUSION-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 101 | ONCO-FUSION-002 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 102 | ONCO-FUSION-003 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 103 | ONCO-CNA-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 104 | ONCO-CNA-002 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 105 | ONCO-CNA-003 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 106 | ONCO-PURITY-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 107 | ONCO-PLOIDY-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 108 | ONCO-CLONAL-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 109 | ONCO-NEO-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 110 | ONCO-MHC-001 | Oncology | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 111 | ONCO-CTDNA-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 112 | ONCO-MRD-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 113 | ONCO-CHIP-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 114 | ONCO-PHYLO-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 115 | ONCO-CCF-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 116 | ONCO-HETERO-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 117 | ONCO-HLA-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 118 | ONCO-ACTION-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 119 | ONCO-SV-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 120 | ONCO-EXPR-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 121 | SEQ-COMPOSITION-001 | Statistics | Analysis | Core !→ Analysis |
| 122 | SEQ-DINUC-001 | Statistics | Analysis | Core !→ Analysis |
| 123 | SEQ-HYDRO-001 | Statistics | Analysis | Core !→ Analysis |
| 124 | SEQ-MW-001 | Statistics | Analysis | Core !→ Analysis |
| 125 | SEQ-PI-001 | Statistics | Analysis | Core !→ Analysis |
| 126 | SEQ-SECSTRUCT-001 | Statistics | Analysis | Core !→ Analysis |
| 127 | SEQ-STATS-001 | Statistics | Analysis | Core !→ Analysis |
| 128 | SEQ-SUMMARY-001 | Statistics | Analysis | Core !→ Analysis |
| 129 | SEQ-THERMO-001 | Statistics | Analysis | Core !→ Analysis |
| 130 | SEQ-TM-001 | Statistics | Analysis | Core !→ Analysis |
| 131 | COMPGEN-ANI-001 | Comparative | Analysis | Core !→ Analysis |
| 132 | COMPGEN-CLUSTER-001 | Comparative | Analysis | Core !→ Analysis |
| 133 | COMPGEN-COMPARE-001 | Comparative | Analysis | Core !→ Analysis |
| 134 | COMPGEN-DOTPLOT-001 | Comparative | Analysis | Core !→ Analysis |
| 135 | COMPGEN-ORTHO-001 | Comparative | Analysis | Core !→ Analysis |
| 136 | COMPGEN-RBH-001 | Comparative | Analysis | Core !→ Analysis |
| 137 | COMPGEN-REARR-001 | Comparative | Analysis | Core !→ Analysis |
| 138 | COMPGEN-REVERSAL-001 | Comparative | Analysis | Core !→ Analysis |
| 139 | COMPGEN-SYNTENY-001 | Comparative | Analysis | Core !→ Analysis |
| 140 | ASSEMBLY-CONSENSUS-001 | Assembly | Alignment | Core !→ Alignment, IO !→ Alignment |
| 141 | ASSEMBLY-CORRECT-001 | Assembly | Alignment | Core !→ Alignment, IO !→ Alignment |
| 142 | ASSEMBLY-COVER-001 | Assembly | Alignment | Core !→ Alignment, IO !→ Alignment |
| 143 | ASSEMBLY-DBG-001 | Assembly | Alignment | Core !→ Alignment, IO !→ Alignment |
| 144 | ASSEMBLY-MERGE-001 | Assembly | Alignment | Core !→ Alignment, IO !→ Alignment |
| 145 | ASSEMBLY-OLC-001 | Assembly | Alignment | Core !→ Alignment, IO !→ Alignment |
| 146 | ASSEMBLY-SCAFFOLD-001 | Assembly | Alignment | Core !→ Alignment, IO !→ Alignment |
| 147 | ASSEMBLY-STATS-001 | Assembly | Chromosome | Core !→ Chromosome |
| 148 | ASSEMBLY-TRIM-001 | Assembly | Alignment | Core !→ Alignment, IO !→ Alignment |
| 149 | RNA-DOTBRACKET-001 | RnaStructure | Analysis | Core !→ Analysis |
| 150 | RNA-HAIRPIN-001 | RnaStructure | Analysis | Core !→ Analysis |
| 151 | RNA-INVERT-001 | RnaStructure | Analysis | Core !→ Analysis |
| 152 | RNA-MFE-001 | RnaStructure | Analysis | Core !→ Analysis |
| 153 | RNA-PAIR-001 | RnaStructure | Annotation | Core !→ Annotation |
| 154 | RNA-PARTITION-001 | RnaStructure | Analysis | Core !→ Analysis |
| 155 | RNA-PSEUDOKNOT-001 | RnaStructure | Analysis | Core !→ Analysis |
| 156 | KMER-ASYNC-001 | K-mer | Analysis | Core !→ Analysis |
| 157 | KMER-BOTH-001 | K-mer | Analysis | Core !→ Analysis |
| 158 | KMER-DIST-001 | K-mer | Analysis | Core !→ Analysis |
| 159 | KMER-GENERATE-001 | K-mer | Analysis | Core !→ Analysis |
| 160 | KMER-POSITIONS-001 | K-mer | Analysis | Core !→ Analysis |
| 161 | KMER-STATS-001 | K-mer | Analysis | Core !→ Analysis |
| 162 | KMER-UNIQUE-001 | K-mer | Analysis | Core !→ Analysis |
| 163 | PROTMOTIF-CC-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 164 | PROTMOTIF-COMMON-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 165 | PROTMOTIF-LC-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 166 | PROTMOTIF-PATTERN-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 167 | PROTMOTIF-SP-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 168 | PROTMOTIF-TM-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 169 | MOTIF-CONS-001 | Matching | Analysis | Core !→ Analysis |
| 170 | MOTIF-DISCOVER-001 | Matching | Analysis | Core !→ Analysis |
| 171 | MOTIF-GENERATE-001 | Matching | Analysis | Core !→ Analysis |
| 172 | MOTIF-REGULATORY-001 | Matching | Analysis | Core !→ Analysis |
| 173 | MOTIF-SHARED-001 | Matching | Analysis | Core !→ Analysis |
| 174 | PAT-APPROX-003 | Matching | Alignment | Core !→ Alignment, IO !→ Alignment |
| 175 | GENOMIC-COMMON-001 | Analysis | Analysis | Core !→ Analysis |
| 176 | GENOMIC-MOTIFS-001 | Analysis | Analysis | Core !→ Analysis |
| 177 | GENOMIC-ORF-001 | Analysis | Analysis | Core !→ Analysis |
| 178 | GENOMIC-REPEAT-001 | Analysis | Analysis | Core !→ Analysis |
| 179 | GENOMIC-SIMILARITY-001 | Analysis | Analysis | Core !→ Analysis |
| 180 | GENOMIC-TANDEM-001 | Analysis | Analysis | Core !→ Analysis |
| 181 | EPIGEN-AGE-001 | Epigenetics | Annotation | Core !→ Annotation |
| 182 | EPIGEN-BISULF-001 | Epigenetics | Annotation | Core !→ Annotation |
| 183 | EPIGEN-CHROM-001 | Epigenetics | Annotation | Core !→ Annotation |
| 184 | EPIGEN-DMR-001 | Epigenetics | Annotation | Core !→ Annotation |
| 185 | EPIGEN-METHYL-001 | Epigenetics | Annotation | Core !→ Annotation |
| 186 | VARIANT-ANNOT-001 | Variants | Annotation | Core !→ Annotation |
| 187 | VARIANT-CALL-001 | Variants | Annotation | Core !→ Annotation |
| 188 | VARIANT-INDEL-001 | Variants | Annotation | Core !→ Annotation |
| 189 | VARIANT-SNP-001 | Variants | Annotation | Core !→ Annotation |
| 190 | PANGEN-CLUSTER-001 | PanGenome | Metagenomics | Core !→ Metagenomics |
| 191 | PANGEN-CORE-001 | PanGenome | Metagenomics | Core !→ Metagenomics |
| 192 | PANGEN-HEAP-001 | PanGenome | Metagenomics | Core !→ Metagenomics |
| 193 | PANGEN-MARKER-001 | PanGenome | Metagenomics | Core !→ Metagenomics |
| 194 | META-FUNC-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 195 | META-PATHWAY-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 196 | META-RESIST-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 197 | META-TAXA-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 198 | TRANS-DIFF-001 | Transcriptome | Annotation | Core !→ Annotation |
| 199 | TRANS-EXPR-001 | Transcriptome | Annotation | Core !→ Annotation |
| 200 | TRANS-SPLICE-001 | Transcriptome | Annotation | Core !→ Annotation |
| 201 | SV-BREAKPOINT-001 | StructuralVar | Annotation | Core !→ Annotation |
| 202 | SV-CNV-001 | StructuralVar | Annotation | Core !→ Annotation |
| 203 | SV-DETECT-001 | StructuralVar | Annotation | Core !→ Annotation |
| 204 | DISORDER-LC-001 | ProteinPred | Analysis | Core !→ Analysis |
| 205 | DISORDER-MORF-001 | ProteinPred | Analysis | Core !→ Analysis |
| 206 | DISORDER-PROPENSITY-001 | ProteinPred | Analysis | Core !→ Analysis |
| 207 | POP-ANCESTRY-001 | PopGen | Population | Core !→ Population |
| 208 | POP-ROH-001 | PopGen | Population | Core !→ Population |
| 209 | POP-SELECT-001 | PopGen | Population | Core !→ Population |
| 210 | SEQ-ATSKEW-001 | Composition | Analysis | Core !→ Analysis |
| 211 | SEQ-REPLICATION-001 | Composition | Analysis | Core !→ Analysis |
| 212 | SEQ-RNACOMP-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 213 | CODON-ENC-001 | Codon | MolTools | Core !→ MolTools |
| 214 | CODON-RSCU-001 | Codon | MolTools | Core !→ MolTools |
| 215 | CODON-STATS-001 | Codon | MolTools | Core !→ MolTools |
| 216 | ANNOT-CODING-001 | Annotation | Annotation | Core !→ Annotation |
| 217 | ANNOT-CODONUSAGE-001 | Annotation | Annotation | Core !→ Annotation |
| 218 | ANNOT-REPEAT-001 | Annotation | Annotation | Core !→ Annotation |
| 219 | QUALITY-PHRED-001 | Quality | IO | IO !→ Analysis, IO !→ Alignment |
| 220 | QUALITY-STATS-001 | Quality | IO | IO !→ Analysis, IO !→ Alignment |
| 221 | PHYLO-BOOT-001 | Phylogenetic | Phylogenetics | Core !→ Phylogenetics, IO !→ Phylogenetics |
| 222 | PHYLO-STATS-001 | Phylogenetic | Phylogenetics | Core !→ Phylogenetics, IO !→ Phylogenetics |
| 223 | TRANS-SIXFRAME-001 | Translation | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 224 | RESTR-FILTER-001 | MolTools | MolTools | Core !→ MolTools |
| 225 | MIRNA-PAIR-001 | MiRNA | Annotation | Core !→ Annotation |
| 226 | ALIGN-STATS-001 | Alignment | Alignment | Core !→ Alignment, IO !→ Alignment |
| 227 | SEQ-CODON-FREQ-001 | Statistics | Analysis | Core !→ Analysis |
| 228 | SEQ-COMPLEX-COMPRESS-001 | Complexity | Analysis | Core !→ Analysis |
| 229 | SEQ-COMPLEX-DUST-001 | Complexity | Analysis | Core !→ Analysis |
| 230 | SEQ-COMPLEX-KMER-001 | Complexity | Analysis | Core !→ Analysis |
| 231 | SEQ-COMPLEX-WINDOW-001 | Complexity | Analysis | Core !→ Analysis |
| 232 | SEQ-ENTROPY-PROFILE-001 | Statistics | Analysis | Core !→ Analysis |
| 233 | SEQ-GC-ANALYSIS-001 | Composition | Analysis | Core !→ Analysis |
| 234 | SEQ-GC-PROFILE-001 | Statistics | Analysis | Core !→ Analysis |
| 235 | ONCO-ASCAT-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 236 | RNA-PKPREDICT-001 | Analysis | Analysis | Core !→ Analysis |
| 237 | RNA-PKRECURSIVE-001 | Analysis | Analysis | Core !→ Analysis |
| 238 | RNA-ACCESS-001 | RnaStructure | Analysis | Core !→ Analysis |
| 239 | PROTMOTIF-HMM-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 240 | PRIMER-NNTM-001 | MolTools | MolTools | Core !→ MolTools |
| 241 | PRIMER-HAIRPIN-001 | MolTools | MolTools | Core !→ MolTools |
| 242 | PRIMER-DIMER-001 | MolTools | MolTools | Core !→ MolTools |
| 243 | PROBE-LNATM-001 | MolTools | MolTools | Core !→ MolTools |
| 244 | PROBE-EVALUE-001 | MolTools | MolTools | Core !→ MolTools |
| 245 | MHC-NN-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 246 | MHC-MATRIX-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 247 | IMMUNE-NUSVR-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 248 | META-CHECKM-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 249 | META-TETRA-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 250 | SPLICE-MAXENT3-001 | Splicing | Annotation | Core !→ Annotation |
| 251 | SPLICE-MAXENT5-001 | Splicing | Annotation | Core !→ Annotation |
| 252 | MIRNA-CONTEXT-001 | MiRNA | Annotation | Core !→ Annotation |
| 253 | MIRNA-PCT-001 | MiRNA | Annotation | Core !→ Annotation |
| 254 | MIRNA-CLASSIFY-001 | MiRNA | Annotation | Core !→ Annotation |
| 255 | MIRNA-CLEAVAGE-001 | MiRNA | Annotation | Core !→ Annotation |
| 256 | REP-APPROX-001 | Repeats | Analysis | Core !→ Analysis |
| 257 | CHROM-ALPHASAT-001 | Chromosome | Chromosome | Core !→ Chromosome |
| 258 | CHROM-HOR-001 | Chromosome | Chromosome | Core !→ Chromosome |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 258 |
| Module-level rules | 22 (all ☑) |
| Modules covered | 13 |
| Rules complete ☑ | 22 |
| Rules remaining ☐ | 0 |

### Algorithms per module

| Module | Algorithms |
|--------|:----------:|
| Core | 9 |
| Analysis | 89 |
| Alignment | 16 |
| IO | 9 |
| Annotation | 37 |
| MolTools | 23 |
| Phylogenetics | 6 |
| Population | 8 |
| Chromosome | 8 |
| Metagenomics | 15 |
| Oncology | 38 |
| **Total** | **258** |
