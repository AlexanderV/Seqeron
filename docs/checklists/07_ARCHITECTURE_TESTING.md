# Checklist 07: Architecture Testing (ArchUnitNET)

**Priority:** P2  
**Framework:** ArchUnitNET  
**Date:** 2026-03-19  
**Total algorithms:** 118

---

## Description

Architecture тести запобігають архітектурному дрейфу, перевіряючи залежності між модулями, naming conventions та структурні правила на рівні IL. Правила застосовуються до модулів (проєктів), а не до окремих алгоритмів.

**Поточне покриття:** `Architecture/ArchitectureTests.cs` — 19 правил (всі ☑): межі шарів Core/IO до всіх 13 модулів, відсутність циклічних залежностей між модулями, заборона System.IO у Core, незмінність Result/DTO-типів.

**Модульна структура:**
- `Seqeron.Genomics.Core` — базові типи, послідовності, розширення
- `Seqeron.Genomics.Analysis` — аналізатори (K-mer, Repeat, miRNA, Disorder, etc.)
- `Seqeron.Genomics.Alignment` — вирівнювання (NW, SW, MSA)
- `Seqeron.Genomics.IO` — парсери файлів (FASTA, FASTQ, BED, VCF, etc.)
- `Seqeron.Genomics.Annotation` — анотація геному (ORF, gene, promoter)
- `Seqeron.Genomics.MolTools` — молекулярні інструменти (CRISPR, primer, restriction)
- `Seqeron.Genomics.Phylogenetics` — філогенетика
- `Seqeron.Genomics.Population` — популяційна генетика
- `Seqeron.Genomics.Chromosome` — хромосомний аналіз
- `Seqeron.Genomics.Metagenomics` — метагеноміка
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

---

## Checklist — Algorithm-to-Module Mapping

Кожен алгоритм верифікується модульними правилами для свого модуля.

| # | Test Unit | Area | Module | Layer Rules Applied |
|---|-----------|------|--------|-------------------|
| 1 | SEQ-GC-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 2 | SEQ-COMP-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 3 | SEQ-REVCOMP-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 4 | SEQ-VALID-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 5 | SEQ-COMPLEX-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 6 | SEQ-ENTROPY-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 7 | SEQ-GCSKEW-001 | Composition | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 8 | PAT-EXACT-001 | Matching | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 9 | PAT-APPROX-001 | Matching | Alignment | IO !→ Alignment |
| 10 | PAT-APPROX-002 | Matching | Alignment | IO !→ Alignment |
| 11 | PAT-IUPAC-001 | Matching | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 12 | PAT-PWM-001 | Matching | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 13 | REP-STR-001 | Repeats | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 14 | REP-TANDEM-001 | Repeats | Core/Analysis | Core !→ Analysis (GenomicAnalyzer in Analysis) |
| 15 | REP-INV-001 | Repeats | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 16 | REP-DIRECT-001 | Repeats | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 17 | REP-PALIN-001 | Repeats | Core | Core !→ Analysis, Core !→ IO, No System.IO |
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
| 58 | CODON-OPT-001 | Codon | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 59 | CODON-CAI-001 | Codon | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 60 | CODON-RARE-001 | Codon | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 61 | CODON-USAGE-001 | Codon | Core | Core !→ Analysis, Core !→ IO, No System.IO |
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
| 74 | MIRNA-SEED-001 | MiRNA | Analysis | Core !→ Analysis |
| 75 | MIRNA-TARGET-001 | MiRNA | Analysis | Core !→ Analysis |
| 76 | MIRNA-PRECURSOR-001 | MiRNA | Analysis | Core !→ Analysis |
| 77 | SPLICE-DONOR-001 | Splicing | Analysis | Core !→ Analysis |
| 78 | SPLICE-ACCEPTOR-001 | Splicing | Analysis | Core !→ Analysis |
| 79 | SPLICE-PREDICT-001 | Splicing | Analysis | Core !→ Analysis |
| 80 | DISORDER-PRED-001 | ProteinPred | Analysis | Core !→ Analysis |
| 81 | DISORDER-REGION-001 | ProteinPred | Analysis | Core !→ Analysis |
| 82 | PROTMOTIF-FIND-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 83 | PROTMOTIF-PROSITE-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 84 | PROTMOTIF-DOMAIN-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 85 | EPIGEN-CPG-001 | Epigenetics | Analysis | Core !→ Analysis |
| 86 | ONCO-IMMUNE-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 87 | SEQ-CODON-FREQ-001 | Statistics | Analysis | Core !→ Analysis |
| 88 | SEQ-COMPLEX-COMPRESS-001 | Complexity | Analysis | Core !→ Analysis |
| 89 | SEQ-COMPLEX-DUST-001 | Complexity | Analysis | Core !→ Analysis |
| 90 | SEQ-COMPLEX-KMER-001 | Complexity | Analysis | Core !→ Analysis |
| 91 | SEQ-COMPLEX-WINDOW-001 | Complexity | Analysis | Core !→ Analysis |
| 92 | SEQ-ENTROPY-PROFILE-001 | Statistics | Analysis | Core !→ Analysis |
| 93 | SEQ-GC-ANALYSIS-001 | Composition | Analysis | Core !→ Analysis |
| 94 | SEQ-GC-PROFILE-001 | Statistics | Analysis | Core !→ Analysis |
| 95 | ONCO-ASCAT-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 96 | RNA-PKPREDICT-001 | RnaStructure | Analysis | Core !→ Analysis |
| 97 | RNA-PKRECURSIVE-001 | RnaStructure | Analysis | Core !→ Analysis |
| 98 | RNA-ACCESS-001 | RnaStructure | Analysis | Core !→ Analysis |
| 99 | PROTMOTIF-HMM-001 | ProteinMotif | Analysis | Core !→ Analysis |
| 100 | PRIMER-NNTM-001 | MolTools | MolTools | Core !→ MolTools |
| 101 | PRIMER-HAIRPIN-001 | MolTools | MolTools | Core !→ MolTools |
| 102 | PRIMER-DIMER-001 | MolTools | MolTools | Core !→ MolTools |
| 103 | PROBE-LNATM-001 | MolTools | MolTools | Core !→ MolTools |
| 104 | PROBE-EVALUE-001 | MolTools | MolTools | Core !→ MolTools |
| 105 | MHC-NN-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 106 | MHC-MATRIX-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 107 | IMMUNE-NUSVR-001 | Oncology | Oncology | Core !→ Oncology, IO !→ Oncology |
| 108 | META-CHECKM-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 109 | META-TETRA-001 | Metagenomics | Metagenomics | Core !→ Metagenomics |
| 110 | SPLICE-MAXENT3-001 | Splicing | Analysis | Core !→ Analysis |
| 111 | SPLICE-MAXENT5-001 | Splicing | Analysis | Core !→ Analysis |
| 112 | MIRNA-CONTEXT-001 | MiRNA | Analysis | Core !→ Analysis |
| 113 | MIRNA-PCT-001 | MiRNA | Analysis | Core !→ Analysis |
| 114 | MIRNA-CLASSIFY-001 | MiRNA | Analysis | Core !→ Analysis |
| 115 | MIRNA-CLEAVAGE-001 | MiRNA | Analysis | Core !→ Analysis |
| 116 | REP-APPROX-001 | Repeats | Core | Core !→ Analysis, Core !→ IO, No System.IO |
| 117 | CHROM-ALPHASAT-001 | Chromosome | Chromosome | Core !→ Chromosome |
| 118 | CHROM-HOR-001 | Chromosome | Chromosome | Core !→ Chromosome |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 118 |
| Module-level rules | 19 (all ☑) |
| Modules covered | 13 |
| Rules complete ☑ | 19 |
| Rules remaining ☐ | 0 |
