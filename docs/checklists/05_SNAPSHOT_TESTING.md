# Checklist 05: Snapshot / Approval Testing (Verify)

**Priority:** P1  
**Framework:** Verify + VerifyNUnit  
**Date:** 2026-03-19  
**Total algorithms:** 86

---

## Description

Snapshot (approval) тесты сериализуют полный выход алгоритма и сравнивают с зафиксированным `.verified.txt` golden master файлом. Любое изменение выхода требует явного одобрения разработчика (обновления снэпшота). Ловит непреднамеренные регрессии в сложных структурированных выходах — выравнивания, деревья, аннотации, результаты парсинга.

**Текущее покрытие:** ~20 snapshot-файлов в `Snapshots/` покрывают: Alignment, Annotation, Chromosome, Codon, CRISPR, Disorder, Epigenetics, FileIO, Metagenomics, MiRNA, MolTools, PatternMatching, Phylogenetic, Population, PrimerProbe, ProteinMotif, Repeat, Restriction, RNA, Splicing.

**Процесс:**
1. Первый запуск: создаёт `.verified.txt` файл
2. Последующие запуски: сравнивают с `.verified.txt`, fail при diff
3. При изменении поведения: ревью diff и accept нового снэпшота

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | What to Snapshot | Snapshot File |
|---|--------|-----------|------|-----------------|---------------|
| 1 | ☐ | SEQ-GC-001 | Composition | GC%, GcFraction для набора известных последовательностей | CompositionSnapshotTests.cs (new) |
| 2 | ☐ | SEQ-COMP-001 | Composition | Complement для известных последовательностей | CompositionSnapshotTests.cs (new) |
| 3 | ☐ | SEQ-REVCOMP-001 | Composition | Reverse-complement для известных последовательностей | CompositionSnapshotTests.cs (new) |
| 4 | ☐ | SEQ-VALID-001 | Composition | Validation результаты для набора valid/invalid входов | CompositionSnapshotTests.cs (new) |
| 5 | ☐ | SEQ-COMPLEX-001 | Composition | Complexity scores для известных последовательностей | CompositionSnapshotTests.cs (new) |
| 6 | ☐ | SEQ-ENTROPY-001 | Composition | Shannon entropy для известных распределений | CompositionSnapshotTests.cs (new) |
| 7 | ☐ | SEQ-GCSKEW-001 | Composition | Cumulative GC skew массив | CompositionSnapshotTests.cs (new) |
| 8 | ☑ | PAT-EXACT-001 | Matching | Exact motif positions и count | PatternMatchingSnapshotTests.cs |
| 9 | ☐ | PAT-APPROX-001 | Matching | Hamming matches для известной пары | PatternMatchingSnapshotTests.cs |
| 10 | ☐ | PAT-APPROX-002 | Matching | Edit distance matches для известной пары | PatternMatchingSnapshotTests.cs |
| 11 | ☑ | PAT-IUPAC-001 | Matching | Degenerate motif matches | PatternMatchingSnapshotTests.cs |
| 12 | ☑ | PAT-PWM-001 | Matching | PWM creation и scan результаты | PatternMatchingSnapshotTests.cs |
| 13 | ☑ | REP-STR-001 | Repeats | Microsatellite полный output | RepeatSnapshotTests.cs |
| 14 | ☐ | REP-TANDEM-001 | Repeats | Tandem repeat полный output | RepeatSnapshotTests.cs |
| 15 | ☑ | REP-INV-001 | Repeats | Inverted repeat полный output | RepeatSnapshotTests.cs |
| 16 | ☑ | REP-DIRECT-001 | Repeats | Direct repeat полный output | RepeatSnapshotTests.cs |
| 17 | ☑ | REP-PALIN-001 | Repeats | Palindrome полный output | RepeatSnapshotTests.cs |
| 18 | ☑ | CRISPR-PAM-001 | MolTools | PAM sites SpCas9 результаты | CrisprSnapshotTests.cs |
| 19 | ☑ | CRISPR-GUIDE-001 | MolTools | Guide RNA design результаты | CrisprSnapshotTests.cs |
| 20 | ☐ | CRISPR-OFF-001 | MolTools | Off-target analysis полный output | CrisprSnapshotTests.cs |
| 21 | ☑ | PRIMER-TM-001 | MolTools | Melting temperature результаты | PrimerProbeSnapshotTests.cs |
| 22 | ☑ | PRIMER-DESIGN-001 | MolTools | Primer design результаты | PrimerProbeSnapshotTests.cs |
| 23 | ☐ | PRIMER-STRUCT-001 | MolTools | Primer structure analysis (hairpin, dimer) | PrimerProbeSnapshotTests.cs |
| 24 | ☑ | PROBE-DESIGN-001 | MolTools | Probe design результаты | PrimerProbeSnapshotTests.cs |
| 25 | ☑ | PROBE-VALID-001 | MolTools | Probe validation результаты | PrimerProbeSnapshotTests.cs |
| 26 | ☑ | RESTR-FIND-001 | MolTools | Restriction site positions | RestrictionSnapshotTests.cs |
| 27 | ☑ | RESTR-DIGEST-001 | MolTools | Digest fragments и summary | RestrictionSnapshotTests.cs |
| 28 | ☑ | ANNOT-ORF-001 | Annotation | ORF positions и frames | AnnotationSnapshotTests.cs |
| 29 | ☑ | ANNOT-GENE-001 | Annotation | Gene prediction результаты | AnnotationSnapshotTests.cs |
| 30 | ☑ | ANNOT-PROM-001 | Annotation | Promoter motif positions | AnnotationSnapshotTests.cs |
| 31 | ☐ | ANNOT-GFF-001 | Annotation | GFF3 serialization output | AnnotationSnapshotTests.cs |
| 32 | ☐ | KMER-COUNT-001 | K-mer | K-mer count table для известной последовательности | KmerSnapshotTests.cs (new) |
| 33 | ☐ | KMER-FREQ-001 | K-mer | K-mer frequency distribution | KmerSnapshotTests.cs (new) |
| 34 | ☐ | KMER-FIND-001 | K-mer | Frequent k-mers полный output | KmerSnapshotTests.cs (new) |
| 35 | ☑ | ALIGN-GLOBAL-001 | Alignment | Global alignment: aligned seqs, score, type | AlignmentSnapshotTests.cs |
| 36 | ☑ | ALIGN-LOCAL-001 | Alignment | Local alignment: aligned seqs, score, positions | AlignmentSnapshotTests.cs |
| 37 | ☐ | ALIGN-SEMI-001 | Alignment | Semi-global alignment результаты | AlignmentSnapshotTests.cs |
| 38 | ☐ | ALIGN-MULTI-001 | Alignment | Multiple alignment результаты | AlignmentSnapshotTests.cs |
| 39 | ☑ | PHYLO-DIST-001 | Phylogenetic | Distance matrix Jukes-Cantor | PhylogeneticSnapshotTests.cs |
| 40 | ☑ | PHYLO-TREE-001 | Phylogenetic | UPGMA и NJ деревья | PhylogeneticSnapshotTests.cs |
| 41 | ☐ | PHYLO-NEWICK-001 | Phylogenetic | Newick string output | PhylogeneticSnapshotTests.cs |
| 42 | ☐ | PHYLO-COMP-001 | Phylogenetic | Tree comparison результаты (RF distance) | PhylogeneticSnapshotTests.cs |
| 43 | ☑ | POP-FREQ-001 | PopGen | Allele frequencies | PopulationSnapshotTests.cs |
| 44 | ☑ | POP-DIV-001 | PopGen | Diversity statistics (π, θ, Tajima's D) | PopulationSnapshotTests.cs |
| 45 | ☑ | POP-HW-001 | PopGen | Hardy-Weinberg test результаты | PopulationSnapshotTests.cs |
| 46 | ☐ | POP-FST-001 | PopGen | F-statistics полный output | PopulationSnapshotTests.cs |
| 47 | ☑ | POP-LD-001 | PopGen | Linkage disequilibrium результаты | PopulationSnapshotTests.cs |
| 48 | ☑ | CHROM-TELO-001 | Chromosome | Telomere analysis результаты | ChromosomeSnapshotTests.cs |
| 49 | ☑ | CHROM-CENT-001 | Chromosome | Centromere analysis результаты | ChromosomeSnapshotTests.cs |
| 50 | ☑ | CHROM-KARYO-001 | Chromosome | Karyotype analysis результаты | ChromosomeSnapshotTests.cs |
| 51 | ☑ | CHROM-ANEU-001 | Chromosome | Aneuploidy detection результаты | ChromosomeSnapshotTests.cs |
| 52 | ☑ | CHROM-SYNT-001 | Chromosome | Synteny blocks результаты | ChromosomeSnapshotTests.cs |
| 53 | ☑ | META-CLASS-001 | Metagenomics | Taxonomic classification результаты | MetagenomicsSnapshotTests.cs |
| 54 | ☑ | META-PROF-001 | Metagenomics | Taxonomic profile | MetagenomicsSnapshotTests.cs |
| 55 | ☑ | META-ALPHA-001 | Metagenomics | Alpha diversity metrics | MetagenomicsSnapshotTests.cs |
| 56 | ☑ | META-BETA-001 | Metagenomics | Beta diversity matrix | MetagenomicsSnapshotTests.cs |
| 57 | ☑ | META-BIN-001 | Metagenomics | Genome binning результаты | MetagenomicsSnapshotTests.cs |
| 58 | ☑ | CODON-OPT-001 | Codon | Optimized sequence | CodonSnapshotTests.cs |
| 59 | ☑ | CODON-CAI-001 | Codon | CAI score | CodonSnapshotTests.cs |
| 60 | ☑ | CODON-RARE-001 | Codon | Rare codons list | CodonSnapshotTests.cs |
| 61 | ☑ | CODON-USAGE-001 | Codon | Codon usage table | CodonSnapshotTests.cs |
| 62 | ☐ | TRANS-CODON-001 | Translation | Genetic code table serialization | TranslationSnapshotTests.cs (new) |
| 63 | ☐ | TRANS-PROT-001 | Translation | Translation result | TranslationSnapshotTests.cs (new) |
| 64 | ☐ | PARSE-FASTA-001 | FileIO | Parsed FASTA records | FileIOSnapshotTests.cs |
| 65 | ☑ | PARSE-FASTQ-001 | FileIO | Parsed FASTQ records и statistics | FileIOSnapshotTests.cs |
| 66 | ☑ | PARSE-BED-001 | FileIO | Parsed BED regions | FileIOSnapshotTests.cs |
| 67 | ☑ | PARSE-VCF-001 | FileIO | Parsed VCF variants | FileIOSnapshotTests.cs |
| 68 | ☑ | PARSE-GFF-001 | FileIO | Parsed GFF features и statistics | FileIOSnapshotTests.cs |
| 69 | ☐ | PARSE-GENBANK-001 | FileIO | Parsed GenBank record | FileIOSnapshotTests.cs |
| 70 | ☐ | PARSE-EMBL-001 | FileIO | Parsed EMBL record | FileIOSnapshotTests.cs |
| 71 | ☑ | RNA-STRUCT-001 | RnaStructure | Predicted structure (pairs, dot-bracket) | RnaSnapshotTests.cs |
| 72 | ☑ | RNA-STEMLOOP-001 | RnaStructure | Stem-loop detection результаты | RnaSnapshotTests.cs |
| 73 | ☑ | RNA-ENERGY-001 | RnaStructure | MFE computation результаты | RnaSnapshotTests.cs |
| 74 | ☑ | MIRNA-SEED-001 | MiRNA | miRNA creation и seed analysis | MiRnaSnapshotTests.cs |
| 75 | ☑ | MIRNA-TARGET-001 | MiRNA | Target site prediction результаты | MiRnaSnapshotTests.cs |
| 76 | ☐ | MIRNA-PRECURSOR-001 | MiRNA | Pre-miRNA analysis output | MiRnaSnapshotTests.cs |
| 77 | ☑ | SPLICE-DONOR-001 | Splicing | Donor site detection результаты | SplicingSnapshotTests.cs |
| 78 | ☑ | SPLICE-ACCEPTOR-001 | Splicing | Acceptor site detection результаты | SplicingSnapshotTests.cs |
| 79 | ☑ | SPLICE-PREDICT-001 | Splicing | Gene structure prediction результаты | SplicingSnapshotTests.cs |
| 80 | ☑ | DISORDER-PRED-001 | ProteinPred | Disorder prediction scores | DisorderSnapshotTests.cs |
| 81 | ☐ | DISORDER-REGION-001 | ProteinPred | Disordered regions output | DisorderSnapshotTests.cs |
| 82 | ☑ | PROTMOTIF-FIND-001 | ProteinMotif | Common motifs полный output | ProteinMotifSnapshotTests.cs |
| 83 | ☑ | PROTMOTIF-PROSITE-001 | ProteinMotif | PROSITE pattern matching output | ProteinMotifSnapshotTests.cs |
| 84 | ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | Domain prediction output | ProteinMotifSnapshotTests.cs |
| 85 | ☑ | EPIGEN-CPG-001 | Epigenetics | CpG island detection output | EpigeneticsSnapshotTests.cs |
| 86 | ☐ | ONCO-IMMUNE-001 | Oncology | Immune infiltration analysis output | OncologySnapshotTests.cs (new) |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 86 |
| ☑ Complete | 55 |
| ☐ Not started | 31 |
| New snapshot test files needed | 5 (Composition, Kmer, Translation, Oncology) |
| Existing snapshot files to extend | ~10 |
