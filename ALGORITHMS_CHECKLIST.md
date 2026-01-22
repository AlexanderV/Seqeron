# Чек-лист алгоритмов SuffixTree.Genomics

**Дата:** 2026-01-22
**Версия:** 1.0

Этот документ содержит полный каталог алгоритмов для решения задач биоинформатики с использованием библиотеки SuffixTree.Genomics.

---

## Оглавление

1. [Поиск и сопоставление паттернов](#1-поиск-и-сопоставление-паттернов)
2. [Анализ повторений](#2-анализ-повторений)
3. [Дизайн молекулярных инструментов](#3-дизайн-молекулярных-инструментов)
4. [Аннотация генома](#4-аннотация-генома)
5. [Анализ хромосом](#5-анализ-хромосом)
6. [Филогенетический анализ](#6-филогенетический-анализ)
7. [Популяционная генетика](#7-популяционная-генетика)
8. [Метагеномный анализ](#8-метагеномный-анализ)
9. [Выравнивание последовательностей](#9-выравнивание-последовательностей)
10. [Анализ K-mer](#10-анализ-k-mer)
11. [Оптимизация кодонов](#11-оптимизация-кодонов)
12. [Вспомогательные алгоритмы](#12-вспомогательные-алгоритмы)

---

## 1. Поиск и сопоставление паттернов

### GenomicAnalyzer

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Найти вхождения мотива | `FindMotif(sequence, motif)` | O(m) |
| ☐ | Найти известные мотивы | `FindKnownMotifs(sequence, motifs)` | O(m × k) |
| ☐ | Найти общие регионы | `FindCommonRegions(seq1, seq2, minLength)` | O(n + m) |
| ☐ | Вычислить сходство | `CalculateSimilarity(seq1, seq2, kmerSize)` | O(n + m) |

### MotifFinder

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Точный поиск мотива | `FindExactMotif(sequence, motif)` | O(m) |
| ☐ | Поиск с IUPAC кодами | `FindDegenerateMotif(sequence, motif)` | O(n × m) |
| ☐ | Сканировать с PWM | `ScanWithPwm(sequence, pwm, threshold)` | O(n × m) |
| ☐ | Создать PWM | `CreatePwm(alignedSequences)` | O(n × m) |
| ☐ | Генерировать консенсус | `GenerateConsensus(alignedSequences)` | O(n × m) |
| ☐ | De novo поиск мотивов | `DiscoverMotifs(sequence, k, minCount)` | O(n) |
| ☐ | Найти общие мотивы | `FindSharedMotifs(sequences, k)` | O(n × k) |
| ☐ | Найти регуляторные элементы | `FindRegulatoryElements(sequence)` | O(n) |

**Известные мотивы:** TATA box, CAAT box, Kozak, Shine-Dalgarno, Poly(A), E-box, AP-1, NF-κB, CREB

### ApproximateMatcher

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Поиск с mismatches (Hamming) | `FindWithMismatches(seq, pattern, maxMismatches)` | O(n × m) |
| ☐ | Поиск с edit distance | `FindWithEdits(seq, pattern, maxEdits)` | O(n × m²) |
| ☐ | Вычислить Hamming distance | `HammingDistance(s1, s2)` | O(m) |
| ☐ | Вычислить Edit distance | `EditDistance(s1, s2)` | O(n × m) |
| ☐ | Найти лучшее совпадение | `FindBestMatch(sequence, pattern)` | O(n × m) |
| ☐ | Частые k-mer с mismatches | `FindFrequentKmersWithMismatches(seq, k, d)` | O(n × 4^k) |

---

## 2. Анализ повторений

### GenomicAnalyzer

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Найти самый длинный повтор | `FindLongestRepeat(sequence)` | O(n) |
| ☐ | Найти все повторы | `FindRepeats(sequence, minLength)` | O(n) |
| ☐ | Найти тандемные повторы | `FindTandemRepeats(seq, minUnit, maxUnit, minReps)` | O(n²) |
| ☐ | Найти палиндромы | `FindPalindromes(sequence, minLength, maxLength)` | O(n²) |

### RepeatFinder

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Найти микросателлиты (STR) | `FindMicrosatellites(seq, minUnit, maxUnit, minRepeats)` | O(n²) |
| ☐ | Найти инвертированные повторы | `FindInvertedRepeats(seq, minArm, maxLoop)` | O(n²) |
| ☐ | Найти прямые повторы | `FindDirectRepeats(seq, minLen, maxLen, minSpacing)` | O(n²) |
| ☐ | Найти палиндромы | `FindPalindromes(sequence, minLength, maxLength)` | O(n²) |
| ☐ | Сводка тандемных повторов | `GetTandemRepeatSummary(sequence, minRepeats)` | O(n²) |

**Типы микросателлитов:** Моно-, ди-, три-, тетра-, пента-, гексануклеотидные повторы

---

## 3. Дизайн молекулярных инструментов

### CrisprDesigner

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Получить CRISPR систему | `GetSystem(systemType)` | O(1) |
| ☐ | Найти PAM сайты | `FindPamSites(sequence, systemType)` | O(n) |
| ☐ | Дизайн guide RNA | `DesignGuideRnas(seq, start, end, systemType)` | O(n) |
| ☐ | Оценить guide RNA | `EvaluateGuideRna(guide, systemType, params)` | O(m) |
| ☐ | Найти off-targets | `FindOffTargets(guide, genome, maxMismatches)` | O(n × m) |
| ☐ | Вычислить специфичность | `CalculateSpecificityScore(guide, genome, type)` | O(n) |

**CRISPR системы:** SpCas9 (NGG), SpCas9-NAG, SaCas9 (NNGRRT), Cas12a (TTTV), AsCas12a, LbCas12a, CasX

### PrimerDesigner

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Дизайн пары праймеров | `DesignPrimers(template, targetStart, targetEnd)` | O(n²) |
| ☐ | Оценить праймер | `EvaluatePrimer(sequence, position, isForward)` | O(m) |
| ☐ | Вычислить Tm (Wallace) | `CalculateMeltingTemperature(primer)` | O(m) |
| ☐ | Вычислить Tm (с солью) | `CalculateMeltingTemperatureWithSalt(primer, Na)` | O(m) |
| ☐ | Вычислить GC content | `CalculateGcContent(sequence)` | O(m) |
| ☐ | Найти гомополимер | `FindLongestHomopolymer(sequence)` | O(m) |
| ☐ | Найти повтор динуклеотида | `FindLongestDinucleotideRepeat(sequence)` | O(m) |
| ☐ | Проверить hairpin | `HasHairpinPotential(sequence, minStemLength)` | O(m²) |
| ☐ | Проверить primer-dimer | `HasPrimerDimer(primer1, primer2)` | O(m) |
| ☐ | Вычислить 3' стабильность | `Calculate3PrimeStability(sequence)` | O(1) |
| ☐ | Генерировать кандидаты | `GeneratePrimerCandidates(template, region)` | O(n²) |

**Параметры:** MinLength=18, MaxLength=25, OptimalTm=60°C, GC=40-60%

### RestrictionAnalyzer

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Получить фермент | `GetEnzyme(enzymeName)` | O(1) |
| ☐ | Найти сайты рестрикции | `FindSites(sequence, enzymeNames)` | O(n × k) |
| ☐ | Найти все сайты | `FindAllSites(sequence)` | O(n × 40) |
| ☐ | Симулировать digest | `Digest(sequence, enzymeNames)` | O(n) |
| ☐ | Сводка digest | `GetDigestSummary(sequence, enzymeNames)` | O(n) |
| ☐ | Создать рестрикционную карту | `CreateMap(sequence, enzymeNames)` | O(n) |
| ☐ | Найти совместимые ферменты | `FindCompatibleEnzymes()` | O(k²) |
| ☐ | Получить blunt cutters | `GetBluntCutters()` | O(k) |
| ☐ | Получить sticky cutters | `GetStickyCutters()` | O(k) |

**База данных:** 40+ ферментов (EcoRI, BamHI, HindIII, NotI, XhoI, SalI, PstI, KpnI, SacI, и др.)

### ProbeDesigner

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Дизайн hybridization probe | `DesignProbe(sequence, target)` | O(n) |
| ☐ | Оценить probe | `EvaluateProbe(probe, params)` | O(m) |

---

## 4. Аннотация генома

### GenomeAnnotator

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Найти ORF | `FindOrfs(dna, minLength, bothStrands, requireStart)` | O(n) |
| ☐ | Найти самый длинный ORF/frame | `FindLongestOrfsPerFrame(dna, bothStrands)` | O(n) |
| ☐ | Найти RBS (Shine-Dalgarno) | `FindRibosomeBindingSites(dna, upstreamWindow)` | O(n) |
| ☐ | Предсказать гены | `PredictGenes(dna, minOrfLength, prefix)` | O(n) |
| ☐ | Парсить GFF3 | `ParseGff3(lines)` | O(n) |
| ☐ | Экспортировать в GFF3 | `ToGff3(annotations, seqId)` | O(k) |
| ☐ | Найти промоторные мотивы | `FindPromoterMotifs(dna)` | O(n) |
| ☐ | Вычислить coding potential | `CalculateCodingPotential(sequence)` | O(n) |
| ☐ | Найти repetitive elements | `FindRepetitiveElements(dna, minLen, minCopies)` | O(n²) |
| ☐ | Получить codon usage | `GetCodonUsage(dna)` | O(n) |

**Промоторные элементы:** -35 box (TTGACA), -10 box (TATAAT/Pribnow box)

### GenomicAnalyzer

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Найти ORF | `FindOpenReadingFrames(sequence, minLength)` | O(n) |

---

## 5. Анализ хромосом

### ChromosomeAnalyzer

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Анализ кариотипа | `AnalyzeKaryotype(chromosomes, ploidyLevel)` | O(k) |
| ☐ | Определить плоидность | `DetectPloidy(depths, expectedDiploidDepth)` | O(n) |
| ☐ | Анализ теломер | `AnalyzeTelomeres(chrName, seq, repeat, searchLen)` | O(n) |
| ☐ | Оценить длину теломер (qPCR) | `EstimateTelomereLengthFromTSRatio(tsRatio)` | O(1) |
| ☐ | Анализ центромеры | `AnalyzeCentromere(chrName, seq, windowSize)` | O(n) |
| ☐ | Предсказать G-bands | `PredictGBands(chrName, seq, bandSize)` | O(n) |
| ☐ | Найти гетерохроматин | `FindHeterochromatinRegions(seq, windowSize)` | O(n) |
| ☐ | Найти блоки синтении | `FindSyntenyBlocks(orthologPairs, minGenes)` | O(n log n) |
| ☐ | Обнаружить перестройки | `DetectRearrangements(syntenyBlocks)` | O(n) |
| ☐ | Обнаружить анеуплоидию | `DetectAneuploidy(depthData, medianDepth)` | O(n) |
| ☐ | Идентифицировать хромосомную анеуплоидию | `IdentifyWholeChromosomeAneuploidy(cnStates)` | O(n) |
| ☐ | Вычислить arm ratio | `CalculateArmRatio(centromerePos, chrLength)` | O(1) |
| ☐ | Классифицировать хромосому | `ClassifyChromosomeByArmRatio(armRatio)` | O(1) |

**Типы хромосом:** Metacentric, Submetacentric, Acrocentric, Telocentric

---

## 6. Филогенетический анализ

### PhylogeneticAnalyzer

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Построить дерево | `BuildTree(seqs, distanceMethod, treeMethod)` | O(n³) |
| ☐ | Вычислить матрицу расстояний | `CalculateDistanceMatrix(seqs, method)` | O(n² × m) |
| ☐ | Попарное расстояние | `CalculatePairwiseDistance(seq1, seq2, method)` | O(m) |
| ☐ | Конвертировать в Newick | `ToNewick(treeNode)` | O(n) |
| ☐ | Парсить Newick | `ParseNewick(newickString)` | O(n) |
| ☐ | Получить листья | `GetLeaves(rootNode)` | O(n) |
| ☐ | Вычислить длину дерева | `CalculateTreeLength(rootNode)` | O(n) |
| ☐ | Получить глубину дерева | `GetTreeDepth(rootNode)` | O(n) |
| ☐ | Robinson-Foulds distance | `RobinsonFouldsDistance(tree1, tree2)` | O(n) |
| ☐ | Найти MRCA | `FindMRCA(root, taxon1, taxon2)` | O(n) |
| ☐ | Patristic distance | `PatristicDistance(root, taxon1, taxon2)` | O(n) |
| ☐ | Bootstrap анализ | `Bootstrap(sequences, replicates)` | O(r × n³) |

**Методы построения дерева:** UPGMA, Neighbor-Joining
**Модели расстояния:** p-distance, Jukes-Cantor, Kimura 2-parameter, Hamming

---

## 7. Популяционная генетика

### PopulationGeneticsAnalyzer

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Вычислить частоту аллелей | `CalculateAlleleFrequencies(hom_maj, het, hom_min)` | O(1) |
| ☐ | Вычислить MAF | `CalculateMAF(genotypes)` | O(n) |
| ☐ | Фильтровать по MAF | `FilterByMAF(variants, minMAF, maxMAF)` | O(n) |
| ☐ | Nucleotide diversity (π) | `CalculateNucleotideDiversity(sequences)` | O(n² × m) |
| ☐ | Watterson theta | `CalculateWattersonTheta(segSites, sampleSize)` | O(1) |
| ☐ | Tajima's D | `CalculateTajimasD(pi, theta, segSites)` | O(1) |
| ☐ | Статистика разнообразия | `CalculateDiversityStatistics(sequences)` | O(n² × m) |
| ☐ | Тест Hardy-Weinberg | `TestHardyWeinberg(variantId, observedCounts)` | O(1) |
| ☐ | Вычислить Fst | `CalculateFst(pop1, pop2)` | O(n) |
| ☐ | F-statistics (Fis, Fit, Fst) | `CalculateFStatistics(variantData)` | O(n) |
| ☐ | Linkage Disequilibrium | `CalculateLD(var1, var2, genotypes)` | O(n) |
| ☐ | Найти haplotype blocks | `FindHaplotypeBlocks(variants)` | O(n²) |
| ☐ | Integrated Haplotype Score | `CalculateIHS(ehh0, ehh1, positions)` | O(n) |
| ☐ | Сканировать на отбор | `ScanForSelection(regions)` | O(n) |
| ☐ | Оценить ancestry | `EstimateAncestry(individuals, refPops)` | O(n × k) |
| ☐ | Найти ROH | `FindROH(genotypes)` | O(n) |
| ☐ | Inbreeding из ROH | `CalculateInbreedingFromROH(rohSegments)` | O(k) |

---

## 8. Метагеномный анализ

### MetagenomicsAnalyzer

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Классифицировать reads | `ClassifyReads(reads, kmerDatabase, k)` | O(n × m) |
| ☐ | Построить k-mer базу | `BuildKmerDatabase(refGenomes, k)` | O(n) |
| ☐ | Генерировать таксономический профиль | `GenerateTaxonomicProfile(classifications)` | O(n) |
| ☐ | Alpha diversity | `CalculateAlphaDiversity(abundances)` | O(n) |
| ☐ | Beta diversity | `CalculateBetaDiversity(sample1, sample2)` | O(n) |
| ☐ | Биннинг контигов | `BinContigs(contigs, numBins, minBinSize)` | O(n) |
| ☐ | Предсказать функции | `PredictFunctions(proteins, funcDatabase)` | O(n × k) |
| ☐ | Функциональное разнообразие | `CalculateFunctionalDiversity(annotations)` | O(n) |
| ☐ | Найти гены резистентности | `FindResistanceGenes(genes, resistanceDB)` | O(n × k) |
| ☐ | Differential abundance | `DifferentialAbundance(cond1, cond2, pThreshold)` | O(n × m) |

**Индексы разнообразия:** Shannon, Simpson, Inverse Simpson, Chao1, Pielou's evenness
**Beta diversity:** Bray-Curtis, Jaccard, UniFrac

---

## 9. Выравнивание последовательностей

### SequenceAligner

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Global alignment (NW) | `GlobalAlign(seq1, seq2, scoring)` | O(n × m) |
| ☐ | Local alignment (SW) | `LocalAlign(seq1, seq2, scoring)` | O(n × m) |
| ☐ | Semi-global alignment | `SemiGlobalAlign(seq1, seq2, scoring)` | O(n × m) |
| ☐ | Статистика выравнивания | `CalculateStatistics(alignment)` | O(n) |
| ☐ | Форматировать выравнивание | `FormatAlignment(alignment, lineWidth)` | O(n) |
| ☐ | Множественное выравнивание | `MultipleAlign(sequences)` | O(n² × m²) |

**Scoring matrices:** SimpleDna, BlastDna, HighIdentityDna

---

## 10. Анализ K-mer

### KmerAnalyzer

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Подсчитать k-mer | `CountKmers(sequence, k)` | O(n) |
| ☐ | K-mer spectrum | `GetKmerSpectrum(sequence, k)` | O(n) |
| ☐ | Самые частые k-mer | `FindMostFrequentKmers(sequence, k)` | O(n) |
| ☐ | Нормализованные частоты | `GetKmerFrequencies(sequence, k)` | O(n) |
| ☐ | K-mer distance | `KmerDistance(seq1, seq2, k)` | O(n + m) |
| ☐ | Уникальные k-mer | `FindUniqueKmers(sequence, k)` | O(n) |
| ☐ | K-mer с min count | `FindKmersWithMinCount(seq, k, minCount)` | O(n) |
| ☐ | Сгенерировать все k-mer | `GenerateAllKmers(k, alphabet)` | O(4^k) |
| ☐ | K-mer entropy | `CalculateKmerEntropy(sequence, k)` | O(n) |
| ☐ | Найти clumps | `FindClumps(seq, k, windowSize, minOccurrences)` | O(n) |
| ☐ | Позиции k-mer | `FindKmerPositions(sequence, kmer)` | O(n) |
| ☐ | K-mer на обоих strand | `CountKmersBothStrands(dnaSeq, k)` | O(n) |
| ☐ | Комплексный анализ k-mer | `AnalyzeKmers(sequence, k)` | O(n) |

---

## 11. Оптимизация кодонов

### CodonOptimizer

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Оптимизировать последовательность | `OptimizeSequence(seq, organism, strategy)` | O(n) |
| ☐ | Вычислить CAI | `CalculateCAI(codingSeq, codonTable)` | O(n) |
| ☐ | Удалить сайты рестрикции | `RemoveRestrictionSites(seq, sites)` | O(n × k) |
| ☐ | Снизить secondary structure | `ReduceSecondaryStructure(seq, codonTable)` | O(n²) |
| ☐ | Найти редкие кодоны | `FindRareCodons(seq, codonTable, threshold)` | O(n) |
| ☐ | Вычислить codon usage | `CalculateCodonUsage(sequence)` | O(n) |
| ☐ | Сравнить codon usage | `CompareCodonUsage(seq1, seq2)` | O(n) |
| ☐ | Создать codon table | `CreateCodonTableFromSequence(refSeq, name)` | O(n) |

**Стратегии:** MaximizeCAI, BalancedOptimization, HarmonizeExpression, MinimizeSecondary, AvoidRareCodons
**Организмы:** E. coli K12, S. cerevisiae, H. sapiens

### CodonUsageAnalyzer

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Анализ использования кодонов | `AnalyzeCodonUsage(sequence)` | O(n) |
| ☐ | RSCU (Relative Synonymous Codon Usage) | `CalculateRSCU(sequence)` | O(n) |

---

## 12. Вспомогательные алгоритмы

### SequenceExtensions

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | GC content (%) | `CalculateGcContent(span)` | O(n) |
| ☐ | GC content (fraction) | `CalculateGcFraction(span)` | O(n) |
| ☐ | GC content (string, %) | `CalculateGcContentFast(string)` | O(n) |
| ☐ | GC content (string, fraction) | `CalculateGcFractionFast(string)` | O(n) |
| ☐ | Complement | `TryGetComplement(source, dest)` | O(n) |
| ☐ | Reverse complement | `TryGetReverseComplement(source, dest)` | O(n) |
| ☐ | Подсчёт k-mer (Span) | `CountKmersSpan(sequence, k)` | O(n) |
| ☐ | Enumerate k-mer | `EnumerateKmers(sequence, k)` | O(n) |
| ☐ | Hamming distance | `HammingDistance(s1, s2)` | O(n) |
| ☐ | Валидация DNA | `IsValidDna(sequence)` | O(n) |
| ☐ | Валидация RNA | `IsValidRna(sequence)` | O(n) |

### DnaSequence

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Complement | `Complement()` | O(n) |
| ☐ | Reverse complement | `ReverseComplement()` | O(n) |
| ☐ | GC content | `GcContent` | O(n) |
| ☐ | Transcribe to RNA | `Transcribe()` | O(n) |
| ☐ | Substring | `Substring(start, length)` | O(m) |
| ☐ | Contains | `Contains(pattern)` | O(m) |
| ☐ | Find all occurrences | `FindAllOccurrences(pattern)` | O(m + k) |
| ☐ | Count occurrences | `CountOccurrences(pattern)` | O(m) |

### GeneticCode

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Translate codon | `Translate(codon)` | O(1) |
| ☐ | Is start codon | `IsStartCodon(codon)` | O(1) |
| ☐ | Is stop codon | `IsStopCodon(codon)` | O(1) |
| ☐ | Get codons for amino acid | `GetCodonsForAminoAcid(aa)` | O(64) |
| ☐ | Get by table number | `GetByTableNumber(tableNumber)` | O(1) |

**Таблицы:** Standard (1), Vertebrate Mitochondrial (2), Yeast Mitochondrial (3), Bacterial/Plastid (11)

### Translator

| ☐ | Задача | Метод | Сложность |
|---|--------|-------|-----------|
| ☐ | Translate DNA to protein | `Translate(dna, geneticCode)` | O(n) |

---

## Дополнительные анализаторы

### Другие классы

| Класс | Описание |
|-------|----------|
| `GcSkewCalculator` | Анализ асимметрии GC |
| `SequenceComplexity` | Вычисление сложности последовательности |
| `QualityScoreAnalyzer` | Анализ PHRED quality scores |
| `VariantCaller` | Определение вариантов |
| `VariantAnnotator` | Аннотация вариантов |
| `StructuralVariantAnalyzer` | Анализ структурных вариантов |
| `ComparativeGenomics` | Сравнительный анализ геномов |
| `GenomeAssemblyAnalyzer` | Анализ качества сборки |
| `PanGenomeAnalyzer` | Анализ pan-genome |
| `RnaSecondaryStructure` | Предсказание вторичной структуры RNA |
| `ProteinMotifFinder` | Поиск мотивов в белках |
| `SpliceSitePredictor` | Предсказание сайтов сплайсинга |
| `MiRnaAnalyzer` | Анализ микроRNA |
| `EpigeneticsAnalyzer` | Эпигенетический анализ |
| `DisorderPredictor` | Предсказание неупорядоченных областей |
| `TranscriptomeAnalyzer` | Анализ транскриптома |

---

## Форматы файлов

### Парсеры

| ☐ | Формат | Класс | Операции |
|---|--------|-------|----------|
| ☐ | FASTA | `FastaParser` | Parse, ParseFile, ParseFileAsync |
| ☐ | FASTQ | `FastqParser` | Parse, ParseFile |
| ☐ | GenBank | `GenBankParser` | Parse, ParseFile |
| ☐ | EMBL | `EmblParser` | Parse, ParseFile |
| ☐ | GFF/GTF | `GffParser` | Parse, ParseFile |
| ☐ | BED | `BedParser` | Parse, ParseFile |
| ☐ | VCF | `VcfParser` | Parse, ParseFile |

### SequenceIO

| ☐ | Задача | Метод |
|---|--------|-------|
| ☐ | Записать FASTA | `WriteFasta(entries, path)` |
| ☐ | Записать FASTQ | `WriteFastq(entries, path)` |

---

## Статистика библиотеки

| Метрика | Значение |
|---------|----------|
| Анализаторов | 21+ |
| Публичных методов | 150+ |
| Ферментов рестрикции | 40+ |
| CRISPR систем | 7 |
| Генетических кодов | 4 |
| Парсеров форматов | 7 |

---

## Использование чек-листа

1. **Отметьте задачу** ☑ когда она решена
2. **Выберите подходящий метод** из соответствующей категории
3. **Обратите внимание на сложность** для больших данных
4. **Используйте Span-based методы** для высокой производительности
5. **Используйте CancellationToken** для длительных операций

---

*Сгенерировано автоматически на основе анализа SuffixTree.Genomics*
