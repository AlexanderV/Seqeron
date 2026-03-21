# Checklist 10: Combinatorial / Pairwise Testing

**Priority:** P3  
**Date:** 2026-03-19  
**Total algorithms:** 86

---

## Description

Combinatorial (pairwise) тестирование генерирует минимальный набор тестовых случаев, покрывающий все пары (или t-кортежи) параметров. Эффективно для алгоритмов с несколькими конфигурационными параметрами, где полный перебор непрактичен.

**Текущее покрытие:** 0. Нет pairwise/combinatorial тестов.

**Инструменты:**
- Pairwise: PICT (Microsoft), AllPairs, NUnit `[Combinatorial]`
- NUnit атрибуты: `[Values]`, `[Range]`, `[Combinatorial]`, `[Pairwise]`

**Когда применять:**
- Алгоритм принимает ≥3 параметра
- Параметры имеют дискретные значения или диапазоны
- Полный перебор: > 100 комбинаций

**Оценка:** Low = ≤2 параметра; Med = 3 параметра; High = ≥4 параметра

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | Parameters (dimensions) | Full Combos | Pairwise Est. | Priority |
|---|--------|-----------|------|------------------------|:-----------:|:-------------:|:--------:|
| 1 | ☐ | SEQ-GC-001 | Composition | seqType(2) × seqLen(3) | 6 | 6 | Low |
| 2 | ☐ | SEQ-COMP-001 | Composition | seqType(3: DNA/RNA/Ambig) × case(2) | 6 | 6 | Low |
| 3 | ☐ | SEQ-REVCOMP-001 | Composition | seqType(3) × seqLen(3) | 9 | 9 | Low |
| 4 | ☐ | SEQ-VALID-001 | Composition | alphabet(4: DNA/RNA/Prot/Amb) × strict(2) × seqLen(3) | 24 | 12 | Med |
| 5 | ☐ | SEQ-COMPLEX-001 | Composition | method(2) × windowSize(3) × seqLen(3) | 18 | 9 | Med |
| 6 | ☐ | SEQ-ENTROPY-001 | Composition | base(2: 2/e) × seqLen(3) | 6 | 6 | Low |
| 7 | ☐ | SEQ-GCSKEW-001 | Composition | windowSize(3) × step(3) × seqLen(3) | 27 | 12 | Med |
| 8 | ☐ | PAT-EXACT-001 | Matching | patLen(3) × seqLen(3) × algorithm(2: ST/KMP) | 18 | 9 | Med |
| 9 | ☐ | PAT-APPROX-001 | Matching | patLen(3) × maxDist(3) × seqLen(3) | 27 | 12 | Med |
| 10 | ☐ | PAT-APPROX-002 | Matching | patLen(3) × maxEdits(3) × seqLen(3) × substCost(2) | 54 | 15 | High |
| 11 | ☐ | PAT-IUPAC-001 | Matching | patLen(3) × ambiguityLevel(3) × seqLen(3) | 27 | 12 | Med |
| 12 | ☐ | PAT-PWM-001 | Matching | motifLen(3) × threshold(3) × seqLen(3) × pseudocount(2) | 54 | 15 | High |
| 13 | ☐ | REP-STR-001 | Repeats | unitLen(3) × minReps(3) × seqLen(3) | 27 | 12 | Med |
| 14 | ☐ | REP-TANDEM-001 | Repeats | minUnitLen(3) × maxUnitLen(3) × minReps(3) × seqLen(3) | 81 | 15 | High |
| 15 | ☐ | REP-INV-001 | Repeats | minArmLen(3) × maxGap(3) × seqLen(3) | 27 | 12 | Med |
| 16 | ☐ | REP-DIRECT-001 | Repeats | minLen(3) × maxGap(3) × seqLen(3) | 27 | 12 | Med |
| 17 | ☐ | REP-PALIN-001 | Repeats | minLen(3) × maxLen(3) × seqLen(3) | 27 | 12 | Med |
| 18 | ☐ | CRISPR-PAM-001 | MolTools | pamType(3: SpCas9/Cas12a/Custom) × strand(3: +/-/both) × seqLen(3) | 27 | 12 | Med |
| 19 | ☐ | CRISPR-GUIDE-001 | MolTools | pamType(3) × guideLen(3) × maxOff(3) × scoringMethod(2) | 54 | 15 | High |
| 20 | ☐ | CRISPR-OFF-001 | MolTools | maxMismatch(3) × seedLen(3) × scoringMethod(2) | 18 | 9 | Med |
| 21 | ☐ | PRIMER-TM-001 | MolTools | method(2: basic/SantaLucia) × saltConc(3) × primerLen(3) | 18 | 9 | Med |
| 22 | ☐ | PRIMER-DESIGN-001 | MolTools | minLen(3) × maxLen(3) × gcRange(3) × tmRange(3) | 81 | 15 | High |
| 23 | ☐ | PRIMER-STRUCT-001 | MolTools | primerLen(3) × saltConc(3) × tempC(3) | 27 | 12 | Med |
| 24 | ☐ | PROBE-DESIGN-001 | MolTools | minLen(3) × maxLen(3) × tmRange(3) × gcRange(3) | 81 | 15 | High |
| 25 | ☐ | PROBE-VALID-001 | MolTools | gcRange(3) × tmRange(3) × selfCompMax(3) | 27 | 12 | Med |
| 26 | ☐ | RESTR-FIND-001 | MolTools | enzyme(4) × topology(2: linear/circular) | 8 | 8 | Low |
| 27 | ☐ | RESTR-DIGEST-001 | MolTools | enzyme(4) × topology(2) × fragments(3) | 24 | 12 | Med |
| 28 | ☐ | ANNOT-ORF-001 | Annotation | minLen(3) × frame(3: fwd/rev/both) × startCodon(2: ATG/alt) | 18 | 9 | Med |
| 29 | ☐ | ANNOT-GENE-001 | Annotation | minOrfLen(3) × rbsWindow(3) × scoring(2) | 18 | 9 | Med |
| 30 | ☐ | ANNOT-PROM-001 | Annotation | threshold(3) × windowSize(3) × motifSet(2) | 18 | 9 | Med |
| 31 | ☐ | ANNOT-GFF-001 | Annotation | featureType(3) × strand(3) × phase(3) | 27 | 12 | Med |
| 32 | ☐ | KMER-COUNT-001 | K-mer | k(4: 3,5,7,11) × seqLen(3) | 12 | 12 | Low |
| 33 | ☐ | KMER-FREQ-001 | K-mer | k(4) × normalization(2) | 8 | 8 | Low |
| 34 | ☐ | KMER-FIND-001 | K-mer | k(4) × topN(3) × minFreq(3) | 36 | 12 | Med |
| 35 | ☐ | ALIGN-GLOBAL-001 | Alignment | matchScore(3) × mismatchPen(3) × gapPen(3) × seqLen(3) | 81 | 15 | High |
| 36 | ☐ | ALIGN-LOCAL-001 | Alignment | matchScore(3) × mismatchPen(3) × gapPen(3) × seqLen(3) | 81 | 15 | High |
| 37 | ☐ | ALIGN-SEMI-001 | Alignment | matchScore(3) × mismatchPen(3) × gapPen(3) | 27 | 12 | Med |
| 38 | ☐ | ALIGN-MULTI-001 | Alignment | nSeqs(3) × seqLen(3) × gapPen(3) × guideTree(2) | 54 | 15 | High |
| 39 | ☐ | PHYLO-DIST-001 | Phylogenetic | model(2: JC/K2P) × nSeqs(3) × seqLen(3) | 18 | 9 | Med |
| 40 | ☐ | PHYLO-TREE-001 | Phylogenetic | method(2: UPGMA/NJ) × nSeqs(3) × seqLen(3) | 18 | 9 | Med |
| 41 | ☐ | PHYLO-NEWICK-001 | Phylogenetic | nLeaves(3) × branchLengths(2: yes/no) | 6 | 6 | Low |
| 42 | ☐ | PHYLO-COMP-001 | Phylogenetic | nLeaves(3) × topology(3: identical/similar/different) | 9 | 9 | Low |
| 43 | ☐ | POP-FREQ-001 | PopGen | nAlleles(3) × nSamples(3) | 9 | 9 | Low |
| 44 | ☐ | POP-DIV-001 | PopGen | nSeqs(3) × seqLen(3) × method(3: π/θ/D) | 27 | 12 | Med |
| 45 | ☐ | POP-HW-001 | PopGen | nGenotypes(3) × nAlleles(2) × correction(2) | 12 | 8 | Low |
| 46 | ☐ | POP-FST-001 | PopGen | nPops(3) × nSamples(3) × method(2: Wright/WC) | 18 | 9 | Med |
| 47 | ☐ | POP-LD-001 | PopGen | nLoci(3) × nSamples(3) × metric(2: D'/r²) | 18 | 9 | Med |
| 48 | ☐ | CHROM-TELO-001 | Chromosome | repeatMotif(2) × minRepeats(3) × seqLen(3) | 18 | 9 | Med |
| 49 | ☐ | CHROM-CENT-001 | Chromosome | windowSize(3) × threshold(3) × seqLen(3) | 27 | 12 | Med |
| 50 | ☐ | CHROM-KARYO-001 | Chromosome | nChrom(3) × armRatio(3) | 9 | 9 | Low |
| 51 | ☐ | CHROM-ANEU-001 | Chromosome | nChrom(3) × depth(3) × threshold(3) | 27 | 12 | Med |
| 52 | ☐ | CHROM-SYNT-001 | Chromosome | nGenes(3) × minBlockSize(3) × nChroms(3) | 27 | 12 | Med |
| 53 | ☐ | META-CLASS-001 | Metagenomics | kmerSize(3) × database(2) × readLen(3) | 18 | 9 | Med |
| 54 | ☐ | META-PROF-001 | Metagenomics | nReads(3) × nTaxa(3) × normalization(2) | 18 | 9 | Med |
| 55 | ☐ | META-ALPHA-001 | Metagenomics | index(3: Shannon/Simpson/Chao1) × nSpecies(3) × evenness(3) | 27 | 12 | Med |
| 56 | ☐ | META-BETA-001 | Metagenomics | metric(3: BC/Jaccard/UniFrac) × nSamples(3) × nSpecies(3) | 27 | 12 | Med |
| 57 | ☐ | META-BIN-001 | Metagenomics | nContigs(3) × features(3: GC/tetra/coverage) × nBins(3) | 27 | 12 | Med |
| 58 | ☐ | CODON-OPT-001 | Codon | organism(3) × strategy(2: freq/random) × seqLen(3) | 18 | 9 | Med |
| 59 | ☐ | CODON-CAI-001 | Codon | referenceTable(3) × seqLen(3) | 9 | 9 | Low |
| 60 | ☐ | CODON-RARE-001 | Codon | threshold(3) × referenceTable(3) × seqLen(3) | 27 | 12 | Med |
| 61 | ☐ | CODON-USAGE-001 | Codon | geneticCode(3) × seqLen(3) | 9 | 9 | Low |
| 62 | ☐ | TRANS-CODON-001 | Translation | tableId(4: 1,2,3,5) × includeAlt(2) | 8 | 8 | Low |
| 63 | ☐ | TRANS-PROT-001 | Translation | frame(3) × tableId(4) × stopHandling(2: stop/readthrough) | 24 | 12 | Med |
| 64 | ☐ | PARSE-FASTA-001 | FileIO | multiRecord(2) × lineWrap(3) × headerStyle(3) | 18 | 9 | Med |
| 65 | ☐ | PARSE-FASTQ-001 | FileIO | qualEncoding(3: Phred33/64/Solexa) × multiRecord(2) × seqLen(3) | 18 | 9 | Med |
| 66 | ☐ | PARSE-BED-001 | FileIO | nFields(3: 3/6/12) × hasHeader(2) × encoding(2) | 12 | 8 | Low |
| 67 | ☐ | PARSE-VCF-001 | FileIO | nSamples(3) × nVariants(3) × vcfVersion(2: 4.1/4.3) × genotypeFormat(3) | 54 | 15 | High |
| 68 | ☐ | PARSE-GFF-001 | FileIO | version(2: GFF3/GTF) × nFeatures(3) × hasSubfeatures(2) | 12 | 8 | Low |
| 69 | ☐ | PARSE-GENBANK-001 | FileIO | nFeatures(3) × hasSequence(2) × division(3) | 18 | 9 | Med |
| 70 | ☐ | PARSE-EMBL-001 | FileIO | nFeatures(3) × hasSequence(2) × division(3) | 18 | 9 | Med |
| 71 | ☐ | RNA-STRUCT-001 | RnaStructure | algorithm(2: Nussinov/MFE) × seqLen(3) × temperature(3) | 18 | 9 | Med |
| 72 | ☐ | RNA-STEMLOOP-001 | RnaStructure | minStem(3) × minLoop(3) × maxLoop(3) | 27 | 12 | Med |
| 73 | ☐ | RNA-ENERGY-001 | RnaStructure | temperature(3) × saltConc(3) × seqLen(3) | 27 | 12 | Med |
| 74 | ☐ | MIRNA-SEED-001 | MiRNA | seedType(3: 6mer/7mer-m8/7mer-A1) × miRnaLen(3) | 9 | 9 | Low |
| 75 | ☐ | MIRNA-TARGET-001 | MiRNA | seedType(3) × utrLen(3) × scoringMethod(2) | 18 | 9 | Med |
| 76 | ☐ | MIRNA-PRECURSOR-001 | MiRNA | precLen(3) × minStem(3) × maxLoop(3) | 27 | 12 | Med |
| 77 | ☐ | SPLICE-DONOR-001 | Splicing | windowSize(3) × threshold(3) × seqLen(3) | 27 | 12 | Med |
| 78 | ☐ | SPLICE-ACCEPTOR-001 | Splicing | windowSize(3) × threshold(3) × seqLen(3) | 27 | 12 | Med |
| 79 | ☐ | SPLICE-PREDICT-001 | Splicing | minIntron(3) × maxIntron(3) × minExon(3) × scoring(2) | 54 | 15 | High |
| 80 | ☐ | DISORDER-PRED-001 | ProteinPred | windowSize(3) × propScale(2) × seqLen(3) | 18 | 9 | Med |
| 81 | ☐ | DISORDER-REGION-001 | ProteinPred | threshold(3) × minLen(3) × mergeGap(3) | 27 | 12 | Med |
| 82 | ☐ | PROTMOTIF-FIND-001 | ProteinMotif | minMotifLen(3) × maxMotifLen(3) × seqLen(3) | 27 | 12 | Med |
| 83 | ☐ | PROTMOTIF-PROSITE-001 | ProteinMotif | nPatterns(3) × matchMode(2: exact/fuzzy) | 6 | 6 | Low |
| 84 | ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | eValueThreshold(3) × minDomainLen(3) × nProfiles(3) | 27 | 12 | Med |
| 85 | ☐ | EPIGEN-CPG-001 | Epigenetics | windowSize(3) × minOE(3) × minGC(3) × minLen(3) | 81 | 15 | High |
| 86 | ☐ | ONCO-IMMUNE-001 | Oncology | geneSetSize(3) × normalization(3) × nPermutations(3) | 27 | 12 | Med |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 86 |
| ☑ Complete | 0 |
| ☐ Not started | 86 |
| High priority (≥4 params, >50 full combos) | 15 |
| Medium priority (3 params) | 52 |
| Low priority (≤2 params) | 19 |
| Total pairwise test cases (estimated) | ~900 |
