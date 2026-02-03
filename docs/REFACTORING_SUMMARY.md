# Финальная сводка пакетов Seqeron.Genomics

## 📊 Распределение файлов по пакетам

| # | Пакет | Файлов | Строк | Уровень | Зависит от |
|---|-------|--------|-------|---------|------------|
| 0 | Seqeron.Genomics.Infrastructure | 3 | ~192 | 0 | SuffixTree |
| 1 | Seqeron.Genomics.Core | 8 | ~1,832 | 1 | Infrastructure, SuffixTree |
| 2 | Seqeron.Genomics.IO | 10 | ~4,872 | 2 | Core |
| 3 | Seqeron.Genomics.Alignment | 4 | ~1,815 | 2 | Core, Infrastructure |
| 4 | Seqeron.Genomics.Analysis | 11 | ~5,461 | 3 | Core, Alignment |
| 5 | Seqeron.Genomics.Annotation | 8 | ~5,590 | 3 | Core, IO, Alignment, Analysis |
| 6 | Seqeron.Genomics.Phylogenetics | 1 | ~654 | 3 | Core, Alignment |
| 7 | Seqeron.Genomics.Population | 1 | ~853 | 3 | Core |
| 8 | Seqeron.Genomics.Metagenomics | 1 | ~582 | 3 | Core, Analysis |
| 9 | Seqeron.Genomics.MolTools | 9 | ~4,980 | 3 | Core, Analysis, Infrastructure |
| 10 | Seqeron.Genomics.Reports | 1 | ~749 | 4 | Все выше |
| 11 | Seqeron.Genomics (мета) | 1 | ~50 | 5 | Все пакеты |
| **ИТОГО** | **12 пакетов** | **58** | **~27,630** | - | - |

**Примечание:** +2 файла за счёт:
- AlignmentTypes.cs создаётся в Infrastructure (вынос из SequenceAligner.cs)
- GlobalUsings.cs в мета-пакете

---

## 📁 Полный список файлов

### Seqeron.Genomics.Infrastructure (3 файла) — НОВЫЙ
```
StatisticsHelper.cs      ← перенос из Genomics
ThermoConstants.cs       ← перенос из Genomics
AlignmentTypes.cs        ← создать: вынести типы из SequenceAligner.cs
```

**Типы в AlignmentTypes.cs:**
- `ScoringMatrix` (record)
- `AlignmentType` (enum) 
- `AlignmentResult` (record)
- `AlignmentStatistics` (record struct)
- `MultipleAlignmentResult` (record)

### Seqeron.Genomics.Core (8 файлов)
```
ISequence.cs
DnaSequence.cs
RnaSequence.cs
ProteinSequence.cs
SequenceExtensions.cs
GeneticCode.cs
Translator.cs
IupacHelper.cs
```

### Seqeron.Genomics.IO (10 файлов)
```
FastaParser.cs
FastqParser.cs
GenBankParser.cs
EmblParser.cs
GffParser.cs
BedParser.cs
VcfParser.cs
SequenceIO.cs
QualityScoreAnalyzer.cs
FeatureLocationHelper.cs
```

### Seqeron.Genomics.Alignment (4 файла)
```
SequenceAligner.cs
ApproximateMatcher.cs
SequenceAssembler.cs
CancellableOperations.cs
```

### Seqeron.Genomics.Analysis (11 файлов)
```
KmerAnalyzer.cs
SequenceStatistics.cs
SequenceComplexity.cs
GcSkewCalculator.cs
GenomicAnalyzer.cs
MotifFinder.cs
ProteinMotifFinder.cs
RepeatFinder.cs
RnaSecondaryStructure.cs
DisorderPredictor.cs
ComparativeGenomics.cs
```

### Seqeron.Genomics.Annotation (8 файлов)
```
GenomeAnnotator.cs
VariantCaller.cs
VariantAnnotator.cs
StructuralVariantAnalyzer.cs
SpliceSitePredictor.cs
TranscriptomeAnalyzer.cs
EpigeneticsAnalyzer.cs
GenomeAssemblyAnalyzer.cs
```

### Seqeron.Genomics.Phylogenetics (1 файл)
```
PhylogeneticAnalyzer.cs
```

### Seqeron.Genomics.Population (1 файл)
```
PopulationGeneticsAnalyzer.cs
```

### Seqeron.Genomics.Metagenomics (1 файл)
```
MetagenomicsAnalyzer.cs
```

### Seqeron.Genomics.MolTools (9 файлов)
```
CrisprDesigner.cs
PrimerDesigner.cs
ProbeDesigner.cs
RestrictionAnalyzer.cs
CodonOptimizer.cs
CodonUsageAnalyzer.cs
MiRnaAnalyzer.cs
ChromosomeAnalyzer.cs
PanGenomeAnalyzer.cs
```

### Seqeron.Genomics.Reports (1 файл)
```
ReportGenerator.cs
```

---

## ✅ Проверка: Все 56 файлов распределены

```
Infrastructure: 2 файла (перенос) + 1 файл (создать AlignmentTypes)
Core:           8 файлов
IO:            10 файлов
Alignment:      4 файла (включая CancellableOperations)
Analysis:      11 файлов
Annotation:     8 файлов
Phylogenetics:  1 файл
Population:     1 файл
Metagenomics:   1 файл
MolTools:       9 файлов
Reports:        1 файл
─────────────────────────
ИТОГО:         56 файлов (перенос) + 1 файл (создать) = 57 файлов ✓
```

**Примечание:** AlignmentTypes.cs — это новый файл, куда выносятся типы из SequenceAligner.cs.
SequenceAligner.cs остаётся в Alignment, но становится меньше (~60 строк выносятся).

---

## 🔗 Граф зависимостей (финальный)

```
                    ┌───────────────────────────────────────────┐
                    │         Seqeron.Genomics (мета-пакет)     │
                    │   TypeForwarding + GlobalUsings           │
                    └─────────────────────┬─────────────────────┘
                                          │
                    ┌─────────────────────┴─────────────────────┐
                    │                                           │
                    ▼                                           ▼
         ┌─────────────────┐                         ┌─────────────────┐
         │     Reports     │                         │   Annotation    │
         │  (1 файл)       │                         │   (8 файлов)    │
         └────────┬────────┘                         └────────┬────────┘
                  │                                           │
                  │            ┌──────────────────────────────┤
                  │            │              │               │
                  ▼            ▼              ▼               ▼
         ┌─────────────┐ ┌─────────────┐ ┌──────────────┐ ┌──────────┐
         │  MolTools   │ │Phylogenetics│ │ Metagenomics │ │Population│
         │ (9 файлов)  │ │  (1 файл)   │ │   (1 файл)   │ │ (1 файл) │
         └──────┬──────┘ └──────┬──────┘ └──────┬───────┘ └────┬─────┘
                │               │               │              │
                └───────────────┼───────────────┘              │
                                │                              │
                                ▼                              │
                       ┌────────────────┐                      │
                       │    Analysis    │◄─────────────────────┘
                       │  (11 файлов)   │
                       └───────┬────────┘
                               │
              ┌────────────────┼────────────────┐
              │                │                │
              ▼                ▼                │
     ┌────────────────┐ ┌────────────┐          │
     │   Alignment    │ │     IO     │          │
     │   (4 файла)    │ │ (10 файлов)│          │
     └───────┬────────┘ └─────┬──────┘          │
             │                │                 │
             └────────┬───────┘                 │
                      │                         │
                      ▼                         │
             ┌────────────────┐                 │
             │      Core      │◄────────────────┘
             │   (8 файлов)   │
             └───────┬────────┘
                     │
                     ▼
          ┌─────────────────────┐
          │   Infrastructure    │  ← НОВЫЙ (базовый уровень)
          │    (3 файла)        │
          └──────────┬──────────┘
                     │
                     ▼
             ┌────────────────┐
             │   SuffixTree   │
             │   (внешний)    │
             └────────────────┘
```

---

## 📋 Чеклист перед началом рефакторинга

### Предварительные проверки
```powershell
# 1. Подсчитать файлы (должно быть 56)
(Get-ChildItem -Path "src/Seqeron/Seqeron.Genomics" -Filter "*.cs" | 
    Where-Object { $_.Name -ne "obj" }).Count

# 2. Сборка проекта
dotnet build Seqeron.sln

# 3. Запуск тестов и сохранение baseline
dotnet test --logger "trx;LogFileName=baseline.trx"

# 4. Подсчёт количества тестов
dotnet test --list-tests | Measure-Object -Line
```

### После рефакторинга
```powershell
# 1. Подсчитать файлы во всех пакетах (должно быть 56)
$packages = @("Core", "IO", "Alignment", "Analysis", "Annotation", 
              "Phylogenetics", "Population", "Metagenomics", "MolTools", "Reports")
$total = 0
foreach ($pkg in $packages) {
    $count = (Get-ChildItem -Path "src/Seqeron/Seqeron.Genomics.$pkg" -Filter "*.cs" -ErrorAction SilentlyContinue).Count
    Write-Host "$pkg : $count"
    $total += $count
}
Write-Host "ИТОГО: $total (должно быть 56)"

# 2. Сборка
dotnet build Seqeron.sln

# 3. Тесты
dotnet test

# 4. Сравнение количества тестов с baseline
dotnet test --list-tests | Measure-Object -Line
```

---

## ⚠️ Критические точки внимания

### 1. CancellableOperations
- **Проблема:** Использует типы из SequenceAligner (ScoringMatrix, AlignmentResult)
- **Решение:** Перенесён в Alignment вместе с SequenceAligner

### 2. ReportGenerator  
- **Проблема:** Зависит от типов из многих пакетов
- **Решение:** Выделен в отдельный пакет Reports на верхнем уровне

### 3. VcfParser ↔ VariantAnnotator
- **Проблема:** VariantAnnotator использует типы из VcfParser
- **Решение:** VcfParser остаётся в IO, Annotation зависит от IO

### 4. ChromosomeAnalyzer
- **Проблема:** Использует RepeatFinder и GcSkewCalculator
- **Решение:** MolTools зависит от Analysis

---

## 🎯 Ожидаемый результат

После рефакторинга:
- ✅ 11 отдельных NuGet-совместимых пакетов
- ✅ Чёткие границы ответственности
- ✅ Однонаправленный граф зависимостей (без циклов)
- ✅ Все 56 файлов сохранены
- ✅ Все тесты проходят
- ✅ Обратная совместимость через мета-пакет
