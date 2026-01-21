# SuffixTree.Genomics - План реалізації

## Архітектурні принципи

### Clean Architecture
```
SuffixTree.Genomics/
├── Core/                    # Доменні моделі та інтерфейси
│   ├── Sequences/           # DnaSequence, RnaSequence, ProteinSequence
│   ├── Interfaces/          # ISequence, ISequenceAnalyzer, IAligner
│   └── Results/             # Result types для алгоритмів
├── Analysis/                # Аналітичні алгоритми
│   ├── Repeats/             # Пошук повторів
│   ├── Motifs/              # Пошук мотивів
│   ├── Structure/           # Вторинна структура
│   └── Statistics/          # Статистичний аналіз
├── Alignment/               # Вирівнювання послідовностей
├── Comparison/              # Порівняльна геноміка
├── Applications/            # Практичні застосування (праймери, CRISPR)
└── IO/                      # Парсери форматів (FASTA, GenBank, etc.)
```

### Coding Standards
- **Immutable результати**: всі Result types - readonly struct
- **Fluent API**: для налаштування параметрів алгоритмів
- **Lazy evaluation**: IEnumerable для великих результатів
- **Memory efficient**: Span/Memory для роботи з послідовностями
- **Cancellation support**: CancellationToken для довгих операцій
- **Nullable annotations**: повна підтримка nullable reference types

---

## Фаза 1: Core Foundation (Week 1-2)

### 1.1 Базові типи послідовностей

#### Завдання
- [ ] `ISequence` - базовий інтерфейс
- [ ] `RnaSequence` - РНК послідовність (A, C, G, U)
- [ ] `ProteinSequence` - амінокислотна послідовність (20 AA)
- [ ] `IupacDnaSequence` - ДНК з IUPAC кодами (N, R, Y, etc.)
- [ ] `QualitySequence` - послідовність з якістю (FASTQ)

#### Файли
```
Core/Sequences/
├── ISequence.cs
├── SequenceBase.cs
├── DnaSequence.cs (refactor existing)
├── RnaSequence.cs
├── ProteinSequence.cs
├── IupacDnaSequence.cs
└── QualitySequence.cs
```

#### Тести (25 тестів)
- [ ] RnaSequenceTests.cs (10 тестів)
- [ ] ProteinSequenceTests.cs (10 тестів)
- [ ] IupacDnaSequenceTests.cs (5 тестів)

---

### 1.2 Генетичний код та трансляція

#### Завдання
- [ ] `GeneticCode` - таблиця кодонів (Standard, Mitochondrial, etc.)
- [ ] `Translator` - трансляція ДНК/РНК → білок
- [ ] `CodonTable` - частоти кодонів для організмів

#### Файли
```
Core/Translation/
├── GeneticCode.cs
├── CodonTable.cs
├── Translator.cs
└── TranslationResult.cs
```

#### Тести (15 тестів)
- [ ] GeneticCodeTests.cs (5 тестів)
- [ ] TranslatorTests.cs (10 тестів)

---

## Фаза 2: DNA Analysis (Week 3-4)

### 2.1 Exact Pattern Matching (вже частково є)

#### Завдання
- [ ] Рефакторинг `FindMotif` → `PatternMatcher`
- [ ] Підтримка IUPAC wildcard patterns
- [ ] Batch pattern matching (множинний пошук)

#### Файли
```
Analysis/Motifs/
├── PatternMatcher.cs
├── IupacPatternMatcher.cs
├── PatternMatchResult.cs
└── BatchPatternMatcher.cs
```

#### Тести (20 тестів)
- [ ] PatternMatcherTests.cs (10 тестів)
- [ ] IupacPatternMatcherTests.cs (10 тестів)

---

### 2.2 Approximate Matching (з мутаціями)

#### Завдання
- [ ] `ApproximateMatcher` - пошук з k помилками
- [ ] Підтримка різних типів помилок:
  - Substitutions (заміни)
  - Insertions (вставки)
  - Deletions (видалення)
- [ ] Edit distance (Levenshtein)
- [ ] Hamming distance (тільки заміни)

#### Файли
```
Analysis/Motifs/
├── ApproximateMatcher.cs
├── EditDistance.cs
├── ApproximateMatchResult.cs
└── MismatchType.cs
```

#### Тести (25 тестів)
- [ ] ApproximateMatcherTests.cs (15 тестів)
- [ ] EditDistanceTests.cs (10 тестів)

---

### 2.3 Repeat Analysis (розширення існуючого)

#### Завдання
- [ ] `SupermaximalRepeatFinder` - максимальні повтори
- [ ] `InvertedRepeatFinder` - інвертовані повтори
- [ ] `MicrosatelliteFinder` - STR (Short Tandem Repeats)
- [ ] `MinisatelliteFinder` - VNTR

#### Файли
```
Analysis/Repeats/
├── RepeatFinder.cs (refactor existing)
├── SupermaximalRepeatFinder.cs
├── InvertedRepeatFinder.cs
├── MicrosatelliteFinder.cs
├── MinisatelliteFinder.cs
└── RepeatClassifier.cs
```

#### Тести (30 тестів)
- [ ] SupermaximalRepeatFinderTests.cs (10 тестів)
- [ ] InvertedRepeatFinderTests.cs (10 тестів)
- [ ] MicrosatelliteFinderTests.cs (10 тестів)

---

### 2.4 Unique Substring Analysis

#### Завдання
- [ ] `ShortestUniqueSubstring` - для праймерів
- [ ] `MinimalUniqueSubstrings` - всі мінімальні унікальні
- [ ] `UniqueKmerFinder` - унікальні k-mers

#### Файли
```
Analysis/Unique/
├── ShortestUniqueSubstringFinder.cs
├── MinimalUniqueSubstringsFinder.cs
├── UniqueKmerFinder.cs
└── UniquenessResult.cs
```

#### Тести (20 тестів)
- [ ] ShortestUniqueSubstringTests.cs (10 тестів)
- [ ] UniqueKmerFinderTests.cs (10 тестів)

---

## Фаза 3: RNA Analysis (Week 5-6)

### 3.1 RNA Secondary Structure

#### Завдання
- [ ] `StemLoopFinder` - пошук stem-loop структур
- [ ] `HairpinFinder` - пошук шпильок
- [ ] `PseudoknotDetector` - детекція псевдовузлів
- [ ] `FreeEnergyCalculator` - ΔG розрахунок

#### Файли
```
Analysis/Structure/
├── Rna/
│   ├── StemLoopFinder.cs
│   ├── HairpinFinder.cs
│   ├── PseudoknotDetector.cs
│   ├── FreeEnergyCalculator.cs
│   ├── SecondaryStructure.cs
│   └── BasePair.cs
```

#### Тести (25 тестів)
- [ ] StemLoopFinderTests.cs (10 тестів)
- [ ] HairpinFinderTests.cs (10 тестів)
- [ ] FreeEnergyCalculatorTests.cs (5 тестів)

---

### 3.2 Splicing Analysis

#### Завдання
- [ ] `SpliceSiteFinder` - донорні/акцепторні сайти
- [ ] `IntronExonPredictor` - передбачення інтронів/екзонів
- [ ] `AlternativeSplicingDetector` - альтернативний сплайсинг

#### Файли
```
Analysis/Splicing/
├── SpliceSiteFinder.cs
├── SpliceSite.cs
├── IntronExonPredictor.cs
└── AlternativeSplicingDetector.cs
```

#### Тести (20 тестів)
- [ ] SpliceSiteFinderTests.cs (10 тестів)
- [ ] IntronExonPredictorTests.cs (10 тестів)

---

### 3.3 miRNA Analysis

#### Завдання
- [ ] `MiRnaTargetPredictor` - передбачення мішеней miRNA
- [ ] `SeedMatcher` - пошук seed region matches
- [ ] `MiRnaScorer` - скоринг взаємодії

#### Файли
```
Analysis/MiRna/
├── MiRnaTargetPredictor.cs
├── SeedMatcher.cs
├── MiRnaScorer.cs
├── TargetSite.cs
└── MiRnaDatabase.cs
```

#### Тести (15 тестів)
- [ ] MiRnaTargetPredictorTests.cs (10 тестів)
- [ ] SeedMatcherTests.cs (5 тестів)

---

## Фаза 4: Protein Analysis (Week 7-8)

### 4.1 Protein Motifs

#### Завдання
- [ ] `ProteinMotifFinder` - PROSITE-подібні патерни
- [ ] `PrositeParser` - парсер PROSITE формату
- [ ] `RegexMotifMatcher` - regex для білків

#### Файли
```
Analysis/Protein/
├── ProteinMotifFinder.cs
├── PrositeParser.cs
├── PrositePattern.cs
├── RegexMotifMatcher.cs
└── ProteinMotifResult.cs
```

#### Тести (20 тестів)
- [ ] ProteinMotifFinderTests.cs (10 тестів)
- [ ] PrositeParserTests.cs (10 тестів)

---

### 4.2 Protein Structure Prediction

#### Завдання
- [ ] `SignalPeptidePredictor` - сигнальні пептиди
- [ ] `TransmembranePredictor` - TM домени
- [ ] `CoiledCoilPredictor` - coiled-coil структури
- [ ] `DisorderPredictor` - неструктуровані регіони

#### Файли
```
Analysis/Protein/
├── SignalPeptidePredictor.cs
├── TransmembranePredictor.cs
├── CoiledCoilPredictor.cs
├── DisorderPredictor.cs
└── ProteinRegion.cs
```

#### Тести (25 тестів)
- [ ] SignalPeptidePredictorTests.cs (7 тестів)
- [ ] TransmembranePredictorTests.cs (8 тестів)
- [ ] CoiledCoilPredictorTests.cs (5 тестів)
- [ ] DisorderPredictorTests.cs (5 тестів)

---

### 4.3 Protein Properties

#### Завдання
- [ ] `ProteinPropertiesCalculator`:
  - Molecular weight
  - Isoelectric point (pI)
  - GRAVY (hydropathicity)
  - Instability index
  - Aliphatic index
- [ ] `AminoAcidComposition` - склад амінокислот

#### Файли
```
Analysis/Protein/
├── ProteinPropertiesCalculator.cs
├── AminoAcidComposition.cs
├── AminoAcidProperties.cs
└── ProteinProperties.cs
```

#### Тести (15 тестів)
- [ ] ProteinPropertiesCalculatorTests.cs (10 тестів)
- [ ] AminoAcidCompositionTests.cs (5 тестів)

---

## Фаза 5: Sequence Alignment (Week 9-10)

### 5.1 Pairwise Alignment

#### Завдання
- [ ] `GlobalAligner` - Needleman-Wunsch
- [ ] `LocalAligner` - Smith-Waterman
- [ ] `SemiGlobalAligner` - overlap alignment
- [ ] `ScoringMatrix` - BLOSUM, PAM матриці

#### Файли
```
Alignment/
├── Pairwise/
│   ├── GlobalAligner.cs
│   ├── LocalAligner.cs
│   ├── SemiGlobalAligner.cs
│   ├── AlignmentResult.cs
│   └── AlignedSequence.cs
├── Scoring/
│   ├── IScoringMatrix.cs
│   ├── NucleotideScoringMatrix.cs
│   ├── BlosumMatrix.cs
│   ├── PamMatrix.cs
│   └── ScoringMatrixLoader.cs
```

#### Тести (30 тестів)
- [ ] GlobalAlignerTests.cs (10 тестів)
- [ ] LocalAlignerTests.cs (10 тестів)
- [ ] ScoringMatrixTests.cs (10 тестів)

---

### 5.2 Multiple Sequence Alignment

#### Завдання
- [ ] `ProgressiveAligner` - progressive MSA
- [ ] `ConsensusBuilder` - побудова консенсусу
- [ ] `AlignmentScorer` - оцінка якості MSA

#### Файли
```
Alignment/
├── Multiple/
│   ├── ProgressiveAligner.cs
│   ├── GuideTree.cs
│   ├── ConsensusBuilder.cs
│   ├── AlignmentScorer.cs
│   └── MultipleAlignment.cs
```

#### Тести (20 тестів)
- [ ] ProgressiveAlignerTests.cs (10 тестів)
- [ ] ConsensusBuilderTests.cs (10 тестів)

---

### 5.3 Suffix Tree-based Alignment

#### Завдання
- [ ] `MummerAligner` - MUM-based alignment
- [ ] `MaximalUniqueMatchFinder` - MUM finder
- [ ] `AnchorChainer` - chaining anchors

#### Файли
```
Alignment/
├── SuffixTreeBased/
│   ├── MummerAligner.cs
│   ├── MaximalUniqueMatchFinder.cs
│   ├── AnchorChainer.cs
│   └── MummerResult.cs
```

#### Тести (20 тестів)
- [ ] MummerAlignerTests.cs (10 тестів)
- [ ] MaximalUniqueMatchFinderTests.cs (10 тестів)

---

## Фаза 6: Comparative Genomics (Week 11-12)

### 6.1 Genome Comparison

#### Завдання
- [ ] `SnpDetector` - SNP detection
- [ ] `IndelDetector` - insertion/deletion detection
- [ ] `SyntenyFinder` - synteny blocks
- [ ] `GeneDuplicationFinder` - gene duplications

#### Файли
```
Comparison/
├── SnpDetector.cs
├── IndelDetector.cs
├── SyntenyFinder.cs
├── GeneDuplicationFinder.cs
├── Snp.cs
├── Indel.cs
└── SyntenyBlock.cs
```

#### Тести (25 тестів)
- [ ] SnpDetectorTests.cs (8 тестів)
- [ ] IndelDetectorTests.cs (7 тестів)
- [ ] SyntenyFinderTests.cs (10 тестів)

---

### 6.2 Pan-Genome Analysis

#### Завдання
- [ ] `CoreGenomeFinder` - core genome
- [ ] `PanGenomeBuilder` - pan genome
- [ ] `AccessoryGeneFinder` - accessory genes

#### Файли
```
Comparison/
├── PanGenome/
│   ├── CoreGenomeFinder.cs
│   ├── PanGenomeBuilder.cs
│   ├── AccessoryGeneFinder.cs
│   └── PanGenomeResult.cs
```

#### Тести (15 тестів)
- [ ] CoreGenomeFinderTests.cs (8 тестів)
- [ ] PanGenomeBuilderTests.cs (7 тестів)

---

## Фаза 7: Statistics (Week 13)

### 7.1 Sequence Statistics

#### Завдання
- [ ] `KmerAnalyzer` - k-mer spectrum
- [ ] `SequenceComplexity` - linguistic complexity
- [ ] `EntropyCalculator` - Shannon entropy
- [ ] `GcSkewCalculator` - GC skew analysis

#### Файли
```
Analysis/Statistics/
├── KmerAnalyzer.cs
├── KmerSpectrum.cs
├── SequenceComplexity.cs
├── EntropyCalculator.cs
├── GcSkewCalculator.cs
└── StatisticsResult.cs
```

#### Тести (25 тестів)
- [ ] KmerAnalyzerTests.cs (10 тестів)
- [ ] SequenceComplexityTests.cs (5 тестів)
- [ ] EntropyCalculatorTests.cs (5 тестів)
- [ ] GcSkewCalculatorTests.cs (5 тестів)

---

### 7.2 Codon Analysis

#### Завдання
- [ ] `CodonUsageAnalyzer` - codon usage bias
- [ ] `CaiCalculator` - Codon Adaptation Index
- [ ] `RscuCalculator` - Relative Synonymous Codon Usage

#### Файли
```
Analysis/Statistics/
├── CodonUsageAnalyzer.cs
├── CaiCalculator.cs
├── RscuCalculator.cs
└── CodonUsageTable.cs
```

#### Тести (15 тестів)
- [ ] CodonUsageAnalyzerTests.cs (8 тестів)
- [ ] CaiCalculatorTests.cs (7 тестів)

---

## Фаза 8: Practical Applications (Week 14-15)

### 8.1 Primer Design

#### Завдання
- [ ] `PrimerDesigner` - дизайн PCR праймерів
- [ ] `PrimerValidator` - валідація праймерів
- [ ] `MeltingTemperatureCalculator` - Tm розрахунок
- [ ] `PrimerDimerChecker` - перевірка димерів

#### Файли
```
Applications/Primers/
├── PrimerDesigner.cs
├── PrimerValidator.cs
├── MeltingTemperatureCalculator.cs
├── PrimerDimerChecker.cs
├── Primer.cs
├── PrimerPair.cs
└── PrimerDesignOptions.cs
```

#### Тести (30 тестів)
- [ ] PrimerDesignerTests.cs (10 тестів)
- [ ] PrimerValidatorTests.cs (10 тестів)
- [ ] MeltingTemperatureCalculatorTests.cs (10 тестів)

---

### 8.2 CRISPR Guide RNA

#### Завдання
- [ ] `GuideRnaDesigner` - дизайн gRNA
- [ ] `PamFinder` - пошук PAM sequences
- [ ] `OffTargetPredictor` - передбачення off-targets
- [ ] `GuideRnaScorer` - скоринг gRNA

#### Файли
```
Applications/Crispr/
├── GuideRnaDesigner.cs
├── PamFinder.cs
├── OffTargetPredictor.cs
├── GuideRnaScorer.cs
├── GuideRna.cs
├── PamSequence.cs
└── CrisprSystem.cs
```

#### Тести (25 тестів)
- [ ] GuideRnaDesignerTests.cs (10 тестів)
- [ ] PamFinderTests.cs (5 тестів)
- [ ] OffTargetPredictorTests.cs (10 тестів)

---

### 8.3 Restriction Analysis

#### Завдання
- [ ] `RestrictionMapper` - restriction map
- [ ] `RestrictionEnzymeDatabase` - база ензимів
- [ ] `DigestSimulator` - симуляція рестрикції
- [ ] `FragmentAnalyzer` - аналіз фрагментів

#### Файли
```
Applications/Restriction/
├── RestrictionMapper.cs
├── RestrictionEnzymeDatabase.cs
├── RestrictionEnzyme.cs
├── DigestSimulator.cs
├── FragmentAnalyzer.cs
└── RestrictionSite.cs
```

#### Тести (20 тестів)
- [ ] RestrictionMapperTests.cs (10 тестів)
- [ ] DigestSimulatorTests.cs (10 тестів)

---

### 8.4 Probe Design

#### Завдання
- [ ] `ProbeDesigner` - дизайн гібридизаційних зондів
- [ ] `ProbeSpecificityChecker` - перевірка специфічності
- [ ] `OligoAnalyzer` - аналіз олігонуклеотидів

#### Файли
```
Applications/Probes/
├── ProbeDesigner.cs
├── ProbeSpecificityChecker.cs
├── OligoAnalyzer.cs
├── Probe.cs
└── ProbeDesignOptions.cs
```

#### Тести (15 тестів)
- [ ] ProbeDesignerTests.cs (10 тестів)
- [ ] OligoAnalyzerTests.cs (5 тестів)

---

## Фаза 9: File Formats (Week 16)

### 9.1 Input Formats

#### Завдання
- [ ] `GenBankParser` - GenBank format
- [ ] `EmblParser` - EMBL format
- [ ] `FastqParser` - FASTQ with quality
- [ ] `GffParser` - GFF/GTF annotations
- [ ] `BedParser` - BED format
- [ ] `VcfParser` - VCF (variants)

#### Файли
```
IO/
├── Parsers/
│   ├── GenBankParser.cs
│   ├── EmblParser.cs
│   ├── FastqParser.cs
│   ├── GffParser.cs
│   ├── BedParser.cs
│   └── VcfParser.cs
├── Models/
│   ├── GenBankRecord.cs
│   ├── FastqRecord.cs
│   ├── GffFeature.cs
│   ├── BedRegion.cs
│   └── VcfVariant.cs
```

#### Тести (35 тестів)
- [ ] GenBankParserTests.cs (8 тестів)
- [ ] FastqParserTests.cs (7 тестів)
- [ ] GffParserTests.cs (7 тестів)
- [ ] BedParserTests.cs (6 тестів)
- [ ] VcfParserTests.cs (7 тестів)

---

### 9.2 Output Formats

#### Завдання
- [ ] `FastaWriter` - FASTA output (refactor)
- [ ] `GenBankWriter` - GenBank output
- [ ] `GffWriter` - GFF output
- [ ] `ReportGenerator` - HTML/JSON reports

#### Файли
```
IO/
├── Writers/
│   ├── FastaWriter.cs
│   ├── GenBankWriter.cs
│   ├── GffWriter.cs
│   └── ReportGenerator.cs
```

#### Тести (15 тестів)
- [ ] WriterTests.cs (15 тестів)

---

## Підсумок

### Загальна статистика

| Фаза | Тиждень | Алгоритми | Тести |
|------|---------|-----------|-------|
| 1. Core Foundation | 1-2 | 8 | 40 |
| 2. DNA Analysis | 3-4 | 12 | 95 |
| 3. RNA Analysis | 5-6 | 10 | 60 |
| 4. Protein Analysis | 7-8 | 12 | 60 |
| 5. Alignment | 9-10 | 9 | 70 |
| 6. Comparative Genomics | 11-12 | 7 | 40 |
| 7. Statistics | 13 | 7 | 40 |
| 8. Applications | 14-15 | 14 | 90 |
| 9. File Formats | 16 | 10 | 50 |
| **TOTAL** | **16 weeks** | **89 algorithms** | **545 tests** |

### Пріоритети (MVP)

**High Priority (реалізувати першими):**
1. ✅ DnaSequence (done)
2. ✅ GenomicAnalyzer basics (done)
3. ✅ FastaParser (done)
4. RnaSequence, ProteinSequence
5. ApproximateMatcher
6. KmerAnalyzer
7. PrimerDesigner

**Medium Priority:**
8. RNA Secondary Structure
9. Pairwise Alignment
10. Restriction Analysis
11. CRISPR gRNA

**Lower Priority:**
12. Multiple Alignment
13. Pan-Genome
14. All file formats

---

## Команди для виконання

```bash
# Створення структури
mkdir -p SuffixTree.Genomics/Core/Sequences
mkdir -p SuffixTree.Genomics/Core/Translation
mkdir -p SuffixTree.Genomics/Analysis/Repeats
mkdir -p SuffixTree.Genomics/Analysis/Motifs
# ... etc

# Запуск тестів
dotnet test --filter "FullyQualifiedName~Genomics"

# Запуск з покриттям
dotnet test --collect:"XPlat Code Coverage"
```

---

## Definition of Done (для кожного алгоритму)

- [ ] Реалізація з XML документацією
- [ ] Unit тести (мін. 5 на алгоритм)
- [ ] Integration тести з реальними даними
- [ ] Performance benchmarks
- [ ] README документація
- [ ] Приклади використання
