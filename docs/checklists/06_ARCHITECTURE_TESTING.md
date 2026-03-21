# Checklist 06: Architecture Testing (ArchUnitNET)

**Priority:** P2  
**Framework:** ArchUnitNET  
**Date:** 2026-03-19  
**Total algorithms:** 86

---

## Description

Architecture тесты предотвращают архитектурный дрейф, проверяя зависимости между модулями, naming conventions и структурные правила на уровне IL. Правила применяются к модулям (проектам), а не к индивидуальным алгоритмам.

**Текущее покрытие:** `Architecture/ArchitectureTests.cs` — 5 правил: Core !→ Analysis, Core !→ IO, Core !→ Alignment, IO !→ Analysis, статические анализаторы.

**Модульная структура:**
- `Seqeron.Genomics.Core` — базовые типы, последовательности, расширения
- `Seqeron.Genomics.Analysis` — анализаторы (K-mer, Repeat, miRNA, Disorder, etc.)
- `Seqeron.Genomics.Alignment` — выравнивание (NW, SW, MSA)
- `Seqeron.Genomics.IO` — парсеры файлов (FASTA, FASTQ, BED, VCF, etc.)
- `Seqeron.Genomics.Annotation` — аннотация генома (ORF, gene, promoter)
- `Seqeron.Genomics.MolTools` — молекулярные инструменты (CRISPR, primer, restriction)
- `Seqeron.Genomics.Phylogenetics` — филогенетика
- `Seqeron.Genomics.Population` — популяционная генетика
- `Seqeron.Genomics.Chromosome` — хромосомный анализ
- `Seqeron.Genomics.Metagenomics` — метагеномика
- `Seqeron.Genomics.Oncology` — онкология
- `Seqeron.Genomics.Infrastructure` — инфраструктура
- `Seqeron.Genomics.Reports` — генерация отчётов

---

## Checklist — Module Dependency Rules

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Rule | Description |
|---|--------|------|-------------|
| 1 | ☑ | Core !→ Analysis | Core не зависит от Analysis |
| 2 | ☑ | Core !→ IO | Core не зависит от IO |
| 3 | ☑ | Core !→ Alignment | Core не зависит от Alignment |
| 4 | ☑ | IO !→ Analysis | IO не зависит от Analysis |
| 5 | ☑ | Analyzers static | Predictor/Finder классы статические |
| 6 | ☐ | Core !→ Annotation | Core не зависит от Annotation |
| 7 | ☐ | Core !→ MolTools | Core не зависит от MolTools |
| 8 | ☐ | Core !→ Phylogenetics | Core не зависит от Phylogenetics |
| 9 | ☐ | Core !→ Population | Core не зависит от Population |
| 10 | ☐ | Core !→ Chromosome | Core не зависит от Chromosome |
| 11 | ☐ | Core !→ Metagenomics | Core не зависит от Metagenomics |
| 12 | ☐ | Core !→ Oncology | Core не зависит от Oncology |
| 13 | ☐ | Core !→ Reports | Core не зависит от Reports |
| 14 | ☐ | IO !→ Alignment | IO не зависит от Alignment |
| 15 | ☐ | IO !→ Phylogenetics | IO не зависит от Phylogenetics |
| 16 | ☐ | IO !→ Oncology | IO не зависит от Oncology |
| 17 | ☐ | No circular deps | Нет циклических зависимостей между модулями |
| 18 | ☐ | No System.IO in Core | Core не использует System.IO напрямую |
| 19 | ☐ | DTOs immutable | Result/DTO типы — records или только getters |

---

## Checklist — Algorithm-to-Module Mapping

Каждый алгоритм верифицируется модульными правилами для своего модуля.

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

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 86 |
| Module-level rules | 5 existing + 14 new = 19 |
| Modules covered | 13 |
| Existing rules ☑ | 5 |
| New rules ☐ | 14 |
