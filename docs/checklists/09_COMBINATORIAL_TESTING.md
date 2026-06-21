# Checklist 09: Combinatorial / Pairwise Testing

**Priority:** P3  
**Date:** 2026-03-19  
**Total algorithms:** 234

---

## Description

Combinatorial (pairwise) тестування генерує мінімальний набір тестових випадків, що покриває всі пари (або t-кортежі) параметрів. Ефективно для алгоритмів із кількома конфігураційними параметрами, де повний перебір непрактичний.

**Поточне покриття:** 0. Немає pairwise/combinatorial тестів.

**Інструменти:**
- Pairwise: PICT (Microsoft), AllPairs, NUnit `[Combinatorial]`
- NUnit атрибути: `[Values]`, `[Range]`, `[Combinatorial]`, `[Pairwise]`

**Коли застосовувати:**
- Алгоритм приймає ≥3 параметри
- Параметри мають дискретні значення або діапазони
- Повний перебір: > 100 комбінацій

**Оцінка:** Low = ≤2 параметри; Med = 3 параметри; High = ≥4 параметри

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete | ✗ Not Applicable

| # | Status | Test Unit | Area | Parameters (dimensions) | Full Combos | Pairwise Est. | Priority |
|---|--------|-----------|------|------------------------|:-----------:|:-------------:|:--------:|
| 1 | ✗ | SEQ-GC-001 | Composition | seqType(2) × seqLen(3) | 6 | 6 | Low |
| 2 | ✗ | SEQ-COMP-001 | Composition | seqType(3: DNA/RNA/Ambig) × case(2) | 6 | 6 | Low |
| 3 | ✗ | SEQ-REVCOMP-001 | Composition | seqType(3) × seqLen(3) | 9 | 9 | Low |
| 4 | ☑ | SEQ-VALID-001 | Composition | alphabet(4: DNA/RNA/Prot/Amb) × strict(2) × seqLen(3) | 24 | 12 | Med |
| 5 | ☑ | SEQ-COMPLEX-001 | Composition | method(2) × windowSize(3) × seqLen(3) | 18 | 9 | Med |
| 6 | ✗ | SEQ-ENTROPY-001 | Composition | base(2: 2/e) × seqLen(3) | 6 | 6 | Low |
| 7 | ☑ | SEQ-GCSKEW-001 | Composition | windowSize(3) × step(3) × seqLen(3) | 27 | 12 | Med |
| 8 | ☑ | PAT-EXACT-001 | Matching | patLen(3) × seqLen(3) × algorithm(2: ST/KMP) | 18 | 9 | Med |
| 9 | ☑ | PAT-APPROX-001 | Matching | patLen(3) × maxDist(3) × seqLen(3) | 27 | 12 | Med |
| 10 | ☑ | PAT-APPROX-002 | Matching | patLen(3) × maxEdits(3) × seqLen(3) × substCost(2) | 54 | 15 | High |
| 11 | ☑ | PAT-IUPAC-001 | Matching | patLen(3) × ambiguityLevel(3) × seqLen(3) | 27 | 12 | Med |
| 12 | ☑ | PAT-PWM-001 | Matching | motifLen(3) × threshold(3) × seqLen(3) × pseudocount(2) | 54 | 15 | High |
| 13 | ☑ | REP-STR-001 | Repeats | unitLen(3) × minReps(3) × seqLen(3) | 27 | 12 | Med |
| 14 | ☑ | REP-TANDEM-001 | Repeats | minUnitLen(3) × maxUnitLen(3) × minReps(3) × seqLen(3) | 81 | 15 | High |
| 15 | ☑ | REP-INV-001 | Repeats | minArmLen(3) × maxGap(3) × seqLen(3) | 27 | 12 | Med |
| 16 | ☑ | REP-DIRECT-001 | Repeats | minLen(3) × maxGap(3) × seqLen(3) | 27 | 12 | Med |
| 17 | ☑ | REP-PALIN-001 | Repeats | minLen(3) × maxLen(3) × seqLen(3) | 27 | 12 | Med |
| 18 | ☑ | CRISPR-PAM-001 | MolTools | pamType(3: SpCas9/Cas12a/Custom) × strand(3: +/-/both) × seqLen(3) | 27 | 12 | Med |
| 19 | ☑ | CRISPR-GUIDE-001 | MolTools | pamType(3) × guideLen(3) × maxOff(3) × scoringMethod(2) | 54 | 15 | High |
| 20 | ☑ | CRISPR-OFF-001 | MolTools | maxMismatch(3) × seedLen(3) × scoringMethod(2) | 18 | 9 | Med |
| 21 | ☑ | PRIMER-TM-001 | MolTools | method(2: basic/SantaLucia) × saltConc(3) × primerLen(3) | 18 | 9 | Med |
| 22 | ☑ | PRIMER-DESIGN-001 | MolTools | minLen(3) × maxLen(3) × gcRange(3) × tmRange(3) | 81 | 15 | High |
| 23 | ☑ | PRIMER-STRUCT-001 | MolTools | primerLen(3) × saltConc(3) × tempC(3) | 27 | 12 | Med |
| 24 | ☑ | PROBE-DESIGN-001 | MolTools | minLen(3) × maxLen(3) × tmRange(3) × gcRange(3) | 81 | 15 | High |
| 25 | ☑ | PROBE-VALID-001 | MolTools | gcRange(3) × tmRange(3) × selfCompMax(3) | 27 | 12 | Med |
| 26 | ✗ | RESTR-FIND-001 | MolTools | enzyme(4) × topology(2: linear/circular) | 8 | 8 | Low |
| 27 | ☑ | RESTR-DIGEST-001 | MolTools | enzyme(4) × topology(2) × fragments(3) | 24 | 12 | Med |
| 28 | ☑ | ANNOT-ORF-001 | Annotation | minLen(3) × frame(3: fwd/rev/both) × startCodon(2: ATG/alt) | 18 | 9 | Med |
| 29 | ☑ | ANNOT-GENE-001 | Annotation | minOrfLen(3) × rbsWindow(3) × scoring(2) | 18 | 9 | Med |
| 30 | ☑ | ANNOT-PROM-001 | Annotation | threshold(3) × windowSize(3) × motifSet(2) | 18 | 9 | Med |
| 31 | ☑ | ANNOT-GFF-001 | Annotation | featureType(3) × strand(3) × phase(3) | 27 | 12 | Med |
| 32 | ✗ | KMER-COUNT-001 | K-mer | k(4: 3,5,7,11) × seqLen(3) | 12 | 12 | Low |
| 33 | ✗ | KMER-FREQ-001 | K-mer | k(4) × normalization(2) | 8 | 8 | Low |
| 34 | ☑ | KMER-FIND-001 | K-mer | k(4) × topN(3) × minFreq(3) | 36 | 12 | Med |
| 35 | ☑ | ALIGN-GLOBAL-001 | Alignment | matchScore(3) × mismatchPen(3) × gapPen(3) × seqLen(3) | 81 | 15 | High |
| 36 | ☑ | ALIGN-LOCAL-001 | Alignment | matchScore(3) × mismatchPen(3) × gapPen(3) × seqLen(3) | 81 | 15 | High |
| 37 | ☑ | ALIGN-SEMI-001 | Alignment | matchScore(3) × mismatchPen(3) × gapPen(3) | 27 | 12 | Med |
| 38 | ☑ | ALIGN-MULTI-001 | Alignment | nSeqs(3) × seqLen(3) × gapPen(3) × guideTree(2) | 54 | 15 | High |
| 39 | ☑ | PHYLO-DIST-001 | Phylogenetic | model(2: JC/K2P) × nSeqs(3) × seqLen(3) | 18 | 9 | Med |
| 40 | ☑ | PHYLO-TREE-001 | Phylogenetic | method(2: UPGMA/NJ) × nSeqs(3) × seqLen(3) | 18 | 9 | Med |
| 41 | ✗ | PHYLO-NEWICK-001 | Phylogenetic | nLeaves(3) × branchLengths(2: yes/no) | 6 | 6 | Low |
| 42 | ✗ | PHYLO-COMP-001 | Phylogenetic | nLeaves(3) × topology(3: identical/similar/different) | 9 | 9 | Low |
| 43 | ✗ | POP-FREQ-001 | PopGen | nAlleles(3) × nSamples(3) | 9 | 9 | Low |
| 44 | ☑ | POP-DIV-001 | PopGen | nSeqs(3) × seqLen(3) × method(3: π/θ/D) | 27 | 12 | Med |
| 45 | ✗ | POP-HW-001 | PopGen | nGenotypes(3) × nAlleles(2) × correction(2) | 12 | 8 | Low |
| 46 | ☑ | POP-FST-001 | PopGen | nPops(3) × nSamples(3) × method(2: Wright/WC) | 18 | 9 | Med |
| 47 | ☑ | POP-LD-001 | PopGen | nLoci(3) × nSamples(3) × metric(2: D'/r²) | 18 | 9 | Med |
| 48 | ☑ | CHROM-TELO-001 | Chromosome | repeatMotif(2) × minRepeats(3) × seqLen(3) | 18 | 9 | Med |
| 49 | ☑ | CHROM-CENT-001 | Chromosome | windowSize(3) × threshold(3) × seqLen(3) | 27 | 12 | Med |
| 50 | ✗ | CHROM-KARYO-001 | Chromosome | nChrom(3) × armRatio(3) | 9 | 9 | Low |
| 51 | ☑ | CHROM-ANEU-001 | Chromosome | nChrom(3) × depth(3) × threshold(3) | 27 | 12 | Med |
| 52 | ☑ | CHROM-SYNT-001 | Chromosome | nGenes(3) × minBlockSize(3) × nChroms(3) | 27 | 12 | Med |
| 53 | ☑ | META-CLASS-001 | Metagenomics | kmerSize(3) × database(2) × readLen(3) | 18 | 9 | Med |
| 54 | ☑ | META-PROF-001 | Metagenomics | nReads(3) × nTaxa(3) × normalization(2) | 18 | 9 | Med |
| 55 | ☑ | META-ALPHA-001 | Metagenomics | index(3: Shannon/Simpson/Chao1) × nSpecies(3) × evenness(3) | 27 | 12 | Med |
| 56 | ☑ | META-BETA-001 | Metagenomics | metric(3: BC/Jaccard/UniFrac) × nSamples(3) × nSpecies(3) | 27 | 12 | Med |
| 57 | ☑ | META-BIN-001 | Metagenomics | nContigs(3) × features(3: GC/tetra/coverage) × nBins(3) | 27 | 12 | Med |
| 58 | ☑ | CODON-OPT-001 | Codon | organism(3) × strategy(2: freq/random) × seqLen(3) | 18 | 9 | Med |
| 59 | ✗ | CODON-CAI-001 | Codon | referenceTable(3) × seqLen(3) | 9 | 9 | Low |
| 60 | ☑ | CODON-RARE-001 | Codon | threshold(3) × referenceTable(3) × seqLen(3) | 27 | 12 | Med |
| 61 | ✗ | CODON-USAGE-001 | Codon | geneticCode(3) × seqLen(3) | 9 | 9 | Low |
| 62 | ✗ | TRANS-CODON-001 | Translation | tableId(4: 1,2,3,5) × includeAlt(2) | 8 | 8 | Low |
| 63 | ☑ | TRANS-PROT-001 | Translation | frame(3) × tableId(4) × stopHandling(2: stop/readthrough) | 24 | 12 | Med |
| 64 | ☑ | PARSE-FASTA-001 | FileIO | multiRecord(2) × lineWrap(3) × headerStyle(3) | 18 | 9 | Med |
| 65 | ☑ | PARSE-FASTQ-001 | FileIO | qualEncoding(3: Phred33/64/Solexa) × multiRecord(2) × seqLen(3) | 18 | 9 | Med |
| 66 | ✗ | PARSE-BED-001 | FileIO | nFields(3: 3/6/12) × hasHeader(2) × encoding(2) | 12 | 8 | Low |
| 67 | ☑ | PARSE-VCF-001 | FileIO | nSamples(3) × nVariants(3) × vcfVersion(2: 4.1/4.3) × genotypeFormat(3) | 54 | 15 | High |
| 68 | ✗ | PARSE-GFF-001 | FileIO | version(2: GFF3/GTF) × nFeatures(3) × hasSubfeatures(2) | 12 | 8 | Low |
| 69 | ☑ | PARSE-GENBANK-001 | FileIO | nFeatures(3) × hasSequence(2) × division(3) | 18 | 9 | Med |
| 70 | ☑ | PARSE-EMBL-001 | FileIO | nFeatures(3) × hasSequence(2) × division(3) | 18 | 9 | Med |
| 71 | ☑ | RNA-STRUCT-001 | RnaStructure | algorithm(2: Nussinov/MFE) × seqLen(3) × temperature(3) | 18 | 9 | Med |
| 72 | ☑ | RNA-STEMLOOP-001 | RnaStructure | minStem(3) × minLoop(3) × maxLoop(3) | 27 | 12 | Med |
| 73 | ☑ | RNA-ENERGY-001 | RnaStructure | temperature(3) × saltConc(3) × seqLen(3) | 27 | 12 | Med |
| 74 | ✗ | MIRNA-SEED-001 | MiRNA | seedType(3: 6mer/7mer-m8/7mer-A1) × miRnaLen(3) | 9 | 9 | Low |
| 75 | ☑ | MIRNA-TARGET-001 | MiRNA | seedType(3) × utrLen(3) × scoringMethod(2) | 18 | 9 | Med |
| 76 | ☑ | MIRNA-PRECURSOR-001 | MiRNA | precLen(3) × minStem(3) × maxLoop(3) | 27 | 12 | Med |
| 77 | ☑ | SPLICE-DONOR-001 | Splicing | windowSize(3) × threshold(3) × seqLen(3) | 27 | 12 | Med |
| 78 | ☑ | SPLICE-ACCEPTOR-001 | Splicing | windowSize(3) × threshold(3) × seqLen(3) | 27 | 12 | Med |
| 79 | ☑ | SPLICE-PREDICT-001 | Splicing | minIntron(3) × maxIntron(3) × minExon(3) × scoring(2) | 54 | 15 | High |
| 80 | ☑ | DISORDER-PRED-001 | ProteinPred | windowSize(3) × propScale(2) × seqLen(3) | 18 | 9 | Med |
| 81 | ☑ | DISORDER-REGION-001 | ProteinPred | threshold(3) × minLen(3) × mergeGap(3) | 27 | 12 | Med |
| 82 | ☑ | PROTMOTIF-FIND-001 | ProteinMotif | minMotifLen(3) × maxMotifLen(3) × seqLen(3) | 27 | 12 | Med |
| 83 | ✗ | PROTMOTIF-PROSITE-001 | ProteinMotif | nPatterns(3) × matchMode(2: exact/fuzzy) | 6 | 6 | Low |
| 84 | ☑ | PROTMOTIF-DOMAIN-001 | ProteinMotif | eValueThreshold(3) × minDomainLen(3) × nProfiles(3) | 27 | 12 | Med |
| 85 | ☑ | EPIGEN-CPG-001 | Epigenetics | windowSize(3) × minOE(3) × minGC(3) × minLen(3) | 81 | 15 | High |
| 86 | ☑ | ONCO-IMMUNE-001 | Oncology | geneSetSize(3) × normalization(3) × nPermutations(3) | 27 | 12 | Med |
| 87 | ☑ | ONCO-SOMATIC-001 | Oncology | minVaf(3) × minDepth(3) × strandBias(2) | 18 | 9 | High |
| 88 | ✗ | ONCO-VAF-001 | Oncology | depth(3) × altCount(3) | 9 | 9 | Low |
| 89 | ☑ | ONCO-DRIVER-001 | Oncology | recurrence(3) × hotspot(2) × geneList(2) | 12 | 8 | Med |
| 90 | ☑ | ONCO-ARTIFACT-001 | Oncology | strandBias(3) × baseQual(3) × position(2) | 18 | 9 | Med |
| 91 | ✗ | ONCO-ANNOT-001 | Oncology | consequenceType(3) × region(3) | 9 | 9 | Low |
| 92 | ☑ | ONCO-TMB-001 | Oncology | panelSize(3) × mutationCount(3) × includeSilent(2) | 18 | 9 | Med |
| 93 | ☑ | ONCO-MSI-001 | Oncology | nLoci(3) × instabilityThreshold(3) | 9 | 9 | Med |
| 94 | ☑ | ONCO-HRD-001 | Oncology | LOH(3) × TAI(3) × LST(3) | 27 | 12 | Med |
| 95 | ✗ | ONCO-LOH-001 | Oncology | bafThreshold(3) × segLen(3) | 9 | 9 | Low |
| 96 | ✗ | ONCO-SIG-001 | Oncology | nVariants(3) × contextSource(2) | 6 | 6 | Low |
| 97 | ☑ | ONCO-SIG-002 | Oncology | nSignatures(3) × nMutations(3) × solver(2) | 18 | 9 | High |
| 98 | ☑ | ONCO-SIG-003 | Oncology | nBootstrap(3) × seed(2) × nMutations(3) | 18 | 9 | Med |
| 99 | ✗ | ONCO-SIG-004 | Oncology | nProcesses(3) × confidenceThreshold(3) | 9 | 9 | Low |
| 100 | ☑ | ONCO-FUSION-001 | Oncology | splitReads(3) × spanningReads(3) × minMapQ(2) | 18 | 9 | High |
| 101 | ✗ | ONCO-FUSION-002 | Oncology | dbSize(3) × matchMode(2) | 6 | 6 | Low |
| 102 | ✗ | ONCO-FUSION-003 | Oncology | breakpointPos(3) × frame(2) | 6 | 6 | Low |
| 103 | ☑ | ONCO-CNA-001 | Oncology | log2Range(3) × binSize(3) × ploidy(2) | 18 | 9 | Med |
| 104 | ☑ | ONCO-CNA-002 | Oncology | cnThreshold(3) × segLen(3) | 9 | 9 | Med |
| 105 | ☑ | ONCO-CNA-003 | Oncology | cnThreshold(3) × segLen(3) | 9 | 9 | Med |
| 106 | ☑ | ONCO-PURITY-001 | Oncology | nVariants(3) × vafDist(3) × cnModel(2) | 18 | 9 | High |
| 107 | ☑ | ONCO-PLOIDY-001 | Oncology | nSegments(3) × cnDist(3) | 9 | 9 | Med |
| 108 | ✗ | ONCO-CLONAL-001 | Oncology | ccfThreshold(3) × nVariants(3) | 9 | 9 | Low |
| 109 | ☑ | ONCO-NEO-001 | Oncology | peptideLen(3) × mutationPos(3) × mutationType(2) | 18 | 9 | Med |
| 110 | ☑ | ONCO-MHC-001 | Oncology | allele(3) × peptideLen(3) × affinityThreshold(2) | 18 | 9 | Med |
| 111 | ☑ | ONCO-CTDNA-001 | Oncology | tumorFraction(3) × depth(3) × errorRate(2) | 18 | 9 | High |
| 112 | ☑ | ONCO-MRD-001 | Oncology | nTrackedVariants(3) × depth(3) × detectionThreshold(2) | 18 | 9 | High |
| 113 | ☑ | ONCO-CHIP-001 | Oncology | geneList(3) × vafBand(3) | 9 | 9 | Med |
| 114 | ☑ | ONCO-PHYLO-001 | Oncology | nClones(3) × nMutations(3) × method(2) | 18 | 9 | Med |
| 115 | ☑ | ONCO-CCF-001 | Oncology | vaf(3) × copyNumber(3) × purity(3) | 27 | 12 | High |
| 116 | ☑ | ONCO-HETERO-001 | Oncology | nVariants(3) × vafSpread(3) | 9 | 9 | Med |
| 117 | ☑ | ONCO-HLA-001 | Oncology | locus(3) × zygosity(2) × coverage(3) | 18 | 9 | Med |
| 118 | ☑ | ONCO-ACTION-001 | Oncology | evidenceLevel(4) × variantType(3) | 12 | 12 | Med |
| 119 | ☑ | ONCO-SV-001 | Oncology | nBreakpoints(3) × clustering(3) × svType(3) | 27 | 12 | Med |
| 120 | ☑ | ONCO-EXPR-001 | Oncology | nGenes(3) × zThreshold(3) × normalization(2) | 18 | 9 | Med |
| 121 | ✗ | SEQ-COMPOSITION-001 | Statistics | alphabet(3) × seqLen(3) | 9 | 9 | Low |
| 122 | ✗ | SEQ-DINUC-001 | Statistics | alphabet(3) × seqLen(3) | 9 | 9 | Low |
| 123 | ☑ | SEQ-HYDRO-001 | Statistics | scale(2) × windowSize(3) × seqLen(3) | 18 | 9 | Med |
| 124 | ✗ | SEQ-MW-001 | Statistics | moleculeType(2: DNA/protein) × seqLen(3) | 6 | 6 | Low |
| 125 | ✗ | SEQ-PI-001 | Statistics | pKaSet(2) × seqLen(3) | 6 | 6 | Low |
| 126 | ☑ | SEQ-SECSTRUCT-001 | Statistics | method(2) × windowSize(3) × seqLen(3) | 18 | 9 | Med |
| 127 | ✗ | SEQ-STATS-001 | Statistics | alphabet(3) × seqLen(3) | 9 | 9 | Low |
| 128 | ✗ | SEQ-SUMMARY-001 | Statistics | seqType(3) × seqLen(3) | 9 | 9 | Low |
| 129 | ☑ | SEQ-THERMO-001 | Statistics | saltConc(3) × seqLen(3) × gcContent(3) | 27 | 12 | Med |
| 130 | ☑ | SEQ-TM-001 | Statistics | method(2) × seqLen(3) × gcContent(3) | 18 | 9 | Med |
| 131 | ☑ | COMPGEN-ANI-001 | Comparative | kmerSize(3) × genomeLen(3) × divergence(3) | 27 | 12 | Med |
| 132 | ☑ | COMPGEN-CLUSTER-001 | Comparative | nGenomes(3) × identityThreshold(3) | 9 | 9 | Med |
| 133 | ✗ | COMPGEN-COMPARE-001 | Comparative | nGenomes(3) × metric(3) | 9 | 9 | Low |
| 134 | ☑ | COMPGEN-DOTPLOT-001 | Comparative | wordSize(3) × seqLen(3) × strand(2) | 18 | 9 | Med |
| 135 | ☐ | COMPGEN-ORTHO-001 | Comparative | nGenes(3) × identityThreshold(3) × eValue(2) | 18 | 9 | Med |
| 136 | ✗ | COMPGEN-RBH-001 | Comparative | nGenes(3) × scoreThreshold(3) | 9 | 9 | Low |
| 137 | ☐ | COMPGEN-REARR-001 | Comparative | nBlocks(3) × minBlockSize(3) | 9 | 9 | Med |
| 138 | ☐ | COMPGEN-REVERSAL-001 | Comparative | nGenes(3) × nReversals(3) | 9 | 9 | Med |
| 139 | ☐ | COMPGEN-SYNTENY-001 | Comparative | nAnchors(3) × minBlockSize(3) × maxGap(3) | 27 | 12 | Med |
| 140 | ☐ | ASSEMBLY-CONSENSUS-001 | Assembly | nReads(3) × coverage(3) × errorRate(2) | 18 | 9 | Med |
| 141 | ☐ | ASSEMBLY-CORRECT-001 | Assembly | k(3) × coverage(3) × errorRate(3) | 27 | 12 | High |
| 142 | ✗ | ASSEMBLY-COVER-001 | Assembly | nReads(3) × readLen(3) | 9 | 9 | Low |
| 143 | ☐ | ASSEMBLY-DBG-001 | Assembly | k(3) × coverage(3) × errorRate(3) | 27 | 12 | High |
| 144 | ☐ | ASSEMBLY-MERGE-001 | Assembly | nContigs(3) × minOverlap(3) | 9 | 9 | Med |
| 145 | ☐ | ASSEMBLY-OLC-001 | Assembly | nReads(3) × minOverlap(3) × errorRate(2) | 18 | 9 | High |
| 146 | ☐ | ASSEMBLY-SCAFFOLD-001 | Assembly | nContigs(3) × nLinks(3) × insertSize(2) | 18 | 9 | Med |
| 147 | ✗ | ASSEMBLY-STATS-001 | Assembly | nContigs(3) × lengthDist(3) | 9 | 9 | Low |
| 148 | ☐ | ASSEMBLY-TRIM-001 | Assembly | qualityCutoff(3) × windowSize(3) × readLen(3) | 27 | 12 | Med |
| 149 | ✗ | RNA-DOTBRACKET-001 | RnaStructure | notation(2) × seqLen(3) | 6 | 6 | Low |
| 150 | ☐ | RNA-HAIRPIN-001 | RnaStructure | minLoop(3) × stemLen(3) | 9 | 9 | Med |
| 151 | ☐ | RNA-INVERT-001 | RnaStructure | minArm(3) × maxGap(3) × seqLen(3) | 27 | 12 | Med |
| 152 | ☐ | RNA-MFE-001 | RnaStructure | algorithm(2) × seqLen(3) × temperature(3) | 18 | 9 | Med |
| 153 | ✗ | RNA-PAIR-001 | RnaStructure | pairType(3) × strand(2) | 6 | 6 | Low |
| 154 | ☐ | RNA-PARTITION-001 | RnaStructure | seqLen(3) × temperature(3) | 9 | 9 | Med |
| 155 | ☐ | RNA-PSEUDOKNOT-001 | RnaStructure | seqLen(3) × maxKnots(3) | 9 | 9 | Med |
| 156 | ☐ | KMER-ASYNC-001 | K-mer | k(4) × seqLen(3) × parallelism(2) | 24 | 12 | Med |
| 157 | ✗ | KMER-BOTH-001 | K-mer | k(4) × seqLen(3) | 12 | 12 | Low |
| 158 | ☐ | KMER-DIST-001 | K-mer | k(4) × metric(2) × seqLen(3) | 24 | 12 | Med |
| 159 | ✗ | KMER-GENERATE-001 | K-mer | k(4) × alphabet(2) | 8 | 8 | Low |
| 160 | ✗ | KMER-POSITIONS-001 | K-mer | k(4) × seqLen(3) | 12 | 12 | Low |
| 161 | ✗ | KMER-STATS-001 | K-mer | k(4) × seqLen(3) | 12 | 12 | Low |
| 162 | ☐ | KMER-UNIQUE-001 | K-mer | k(4) × minCount(3) × seqLen(3) | 36 | 12 | Med |
| 163 | ☐ | PROTMOTIF-CC-001 | ProteinMotif | windowSize(3) × threshold(3) × seqLen(3) | 27 | 12 | Med |
| 164 | ☐ | PROTMOTIF-COMMON-001 | ProteinMotif | nSeqs(3) × motifLen(3) | 9 | 9 | Med |
| 165 | ☐ | PROTMOTIF-LC-001 | ProteinMotif | windowSize(3) × threshold(3) | 9 | 9 | Med |
| 166 | ✗ | PROTMOTIF-PATTERN-001 | ProteinMotif | nPatterns(3) × matchMode(2) | 6 | 6 | Low |
| 167 | ☐ | PROTMOTIF-SP-001 | ProteinMotif | windowSize(3) × threshold(3) | 9 | 9 | Med |
| 168 | ☐ | PROTMOTIF-TM-001 | ProteinMotif | windowSize(3) × threshold(3) × seqLen(3) | 27 | 12 | Med |
| 169 | ✗ | MOTIF-CONS-001 | Matching | nSeqs(3) × width(3) | 9 | 9 | Low |
| 170 | ☐ | MOTIF-DISCOVER-001 | Matching | k(3) × support(3) × seqLen(3) | 27 | 12 | Med |
| 171 | ✗ | MOTIF-GENERATE-001 | Matching | width(3) × pseudocount(2) | 6 | 6 | Low |
| 172 | ✗ | MOTIF-REGULATORY-001 | Matching | elementSet(2) × seqLen(3) | 6 | 6 | Low |
| 173 | ☐ | MOTIF-SHARED-001 | Matching | nSeqs(3) × k(3) | 9 | 9 | Med |
| 174 | ☐ | PAT-APPROX-003 | Matching | patLen(3) × maxDist(3) × seqLen(3) | 27 | 12 | Med |
| 175 | ☐ | GENOMIC-COMMON-001 | Analysis | nSeqs(3) × minLen(3) | 9 | 9 | Med |
| 176 | ✗ | GENOMIC-MOTIFS-001 | Analysis | motifSet(2) × seqLen(3) | 6 | 6 | Low |
| 177 | ☐ | GENOMIC-ORF-001 | Analysis | minLen(3) × frame(3) × startCodon(2) | 18 | 9 | Med |
| 178 | ☐ | GENOMIC-REPEAT-001 | Analysis | minLen(3) × seqLen(3) | 9 | 9 | Med |
| 179 | ✗ | GENOMIC-SIMILARITY-001 | Analysis | method(2) × seqLen(3) | 6 | 6 | Low |
| 180 | ☐ | GENOMIC-TANDEM-001 | Analysis | unitLen(3) × minReps(3) × seqLen(3) | 27 | 12 | Med |
| 181 | ☐ | EPIGEN-AGE-001 | Epigenetics | clock(2) × nSites(3) | 6 | 6 | Med |
| 182 | ✗ | EPIGEN-BISULF-001 | Epigenetics | methylationRate(3) × seqLen(3) | 9 | 9 | Low |
| 183 | ☐ | EPIGEN-CHROM-001 | Epigenetics | nMarks(3) × windowSize(3) | 9 | 9 | Med |
| 184 | ☐ | EPIGEN-DMR-001 | Epigenetics | threshold(3) × minSites(3) × nSamples(2) | 18 | 9 | Med |
| 185 | ✗ | EPIGEN-METHYL-001 | Epigenetics | context(3) × coverage(3) | 9 | 9 | Low |
| 186 | ☐ | VARIANT-ANNOT-001 | Variants | consequence(3) × region(3) | 9 | 9 | Med |
| 187 | ☐ | VARIANT-CALL-001 | Variants | minVaf(3) × minDepth(3) × minQual(3) | 27 | 12 | High |
| 188 | ☐ | VARIANT-INDEL-001 | Variants | indelLen(3) × depth(3) | 9 | 9 | Med |
| 189 | ☐ | VARIANT-SNP-001 | Variants | minVaf(3) × depth(3) × baseQual(3) | 27 | 12 | High |
| 190 | ☐ | PANGEN-CLUSTER-001 | PanGenome | nGenes(3) × identity(3) | 9 | 9 | Med |
| 191 | ☐ | PANGEN-CORE-001 | PanGenome | nGenomes(3) × coreThreshold(3) | 9 | 9 | Med |
| 192 | ✗ | PANGEN-HEAP-001 | PanGenome | nGenomes(3) × permutations(2) | 6 | 6 | Low |
| 193 | ☐ | PANGEN-MARKER-001 | PanGenome | nGenomes(3) × nMarkers(3) | 9 | 9 | Med |
| 194 | ☐ | META-FUNC-001 | Metagenomics | dbSize(3) × nGenes(3) | 9 | 9 | Med |
| 195 | ☐ | META-PATHWAY-001 | Metagenomics | nGenes(3) × pathwaySet(2) | 6 | 6 | Med |
| 196 | ☐ | META-RESIST-001 | Metagenomics | dbSize(3) × identity(3) | 9 | 9 | Med |
| 197 | ☐ | META-TAXA-001 | Metagenomics | nSamples(3) × nTaxa(3) × test(2) | 18 | 9 | Med |
| 198 | ☐ | TRANS-DIFF-001 | Transcriptome | nReplicates(3) × foldChange(3) × test(2) | 18 | 9 | High |
| 199 | ☐ | TRANS-EXPR-001 | Transcriptome | nReads(3) × normalization(3) | 9 | 9 | Med |
| 200 | ☐ | TRANS-SPLICE-001 | Transcriptome | nIsoforms(3) × junctionReads(3) | 9 | 9 | Med |
| 201 | ☐ | SV-BREAKPOINT-001 | StructuralVar | splitReads(3) × minMapQ(3) | 9 | 9 | Med |
| 202 | ☐ | SV-CNV-001 | StructuralVar | binSize(3) × coverageRatio(3) × ploidy(2) | 18 | 9 | High |
| 203 | ☐ | SV-DETECT-001 | StructuralVar | svType(3) × minSize(3) × support(3) | 27 | 12 | High |
| 204 | ☐ | DISORDER-LC-001 | ProteinPred | windowSize(3) × threshold(3) | 9 | 9 | Med |
| 205 | ☐ | DISORDER-MORF-001 | ProteinPred | windowSize(3) × threshold(3) × seqLen(3) | 27 | 12 | Med |
| 206 | ✗ | DISORDER-PROPENSITY-001 | ProteinPred | scale(2) × windowSize(3) | 6 | 6 | Low |
| 207 | ☐ | POP-ANCESTRY-001 | PopGen | nPops(3) × nMarkers(3) × method(2) | 18 | 9 | Med |
| 208 | ☐ | POP-ROH-001 | PopGen | minLen(3) × density(3) | 9 | 9 | Med |
| 209 | ☐ | POP-SELECT-001 | PopGen | statistic(3) × windowSize(3) × nPops(2) | 18 | 9 | Med |
| 210 | ☐ | SEQ-ATSKEW-001 | Composition | windowSize(3) × step(3) × seqLen(3) | 27 | 12 | Med |
| 211 | ☐ | SEQ-REPLICATION-001 | Composition | windowSize(3) × seqLen(3) | 9 | 9 | Med |
| 212 | ✗ | SEQ-RNACOMP-001 | Composition | base(4) × case(2) | 8 | 8 | Low |
| 213 | ☐ | CODON-ENC-001 | Codon | geneticCode(3) × seqLen(3) | 9 | 9 | Med |
| 214 | ✗ | CODON-RSCU-001 | Codon | geneticCode(3) × seqLen(3) | 9 | 9 | Low |
| 215 | ✗ | CODON-STATS-001 | Codon | geneticCode(3) × seqLen(3) | 9 | 9 | Low |
| 216 | ☐ | ANNOT-CODING-001 | Annotation | method(2) × seqLen(3) | 6 | 6 | Med |
| 217 | ✗ | ANNOT-CODONUSAGE-001 | Annotation | geneticCode(3) × seqLen(3) | 9 | 9 | Low |
| 218 | ☐ | ANNOT-REPEAT-001 | Annotation | minLen(3) × repeatType(3) | 9 | 9 | Med |
| 219 | ☐ | QUALITY-PHRED-001 | Quality | encoding(3) × seqLen(3) | 9 | 9 | Med |
| 220 | ✗ | QUALITY-STATS-001 | Quality | encoding(3) × seqLen(3) | 9 | 9 | Low |
| 221 | ☐ | PHYLO-BOOT-001 | Phylogenetic | nReplicates(3) × method(2) × nSeqs(3) | 18 | 9 | Med |
| 222 | ✗ | PHYLO-STATS-001 | Phylogenetic | nLeaves(3) × topology(3) | 9 | 9 | Low |
| 223 | ☐ | TRANS-SIXFRAME-001 | Translation | tableId(4) × seqLen(3) | 12 | 12 | Med |
| 224 | ☐ | RESTR-FILTER-001 | MolTools | enzyme(3) × criteria(3) | 9 | 9 | Med |
| 225 | ☐ | MIRNA-PAIR-001 | MiRNA | seedType(3) × utrLen(3) | 9 | 9 | Med |
| 226 | ✗ | ALIGN-STATS-001 | Alignment | alignType(3) × seqLen(3) | 9 | 9 | Low |
| 227 | ✗ | SEQ-CODON-FREQ-001 | Statistics | seqType(2) × seqLen(3) | 6 | 6 | Low |
| 228 | ✗ | SEQ-COMPLEX-COMPRESS-001 | Complexity | seqType(3) × seqLen(3) | 9 | 9 | Low |
| 229 | ☐ | SEQ-COMPLEX-DUST-001 | Complexity | seqType(3) × window(3) | 9 | 9 | Med |
| 230 | ☐ | SEQ-COMPLEX-KMER-001 | Complexity | seqType(3) × k(3) | 9 | 9 | Med |
| 231 | ☐ | SEQ-COMPLEX-WINDOW-001 | Complexity | seqType(3) × window(3) × step(2) | 18 | 9 | Med |
| 232 | ✗ | SEQ-ENTROPY-PROFILE-001 | Statistics | seqType(3) × window(3) | 9 | 9 | Low |
| 233 | ✗ | SEQ-GC-ANALYSIS-001 | Composition | seqType(2) × window(3) | 6 | 6 | Low |
| 234 | ✗ | SEQ-GC-PROFILE-001 | Statistics | seqType(2) × window(3) | 6 | 6 | Low |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 234 |
| ☑ Complete | 100 |
| ☐ Not started | 69 |
| ✗ Not applicable | 65 |
| High priority (≥4 params, >50 full combos) | 15 |
| Medium priority (3 params) | 52 |
| Low priority (≤2 params) | 19 |
| Total pairwise test cases (estimated) | ~900 |
