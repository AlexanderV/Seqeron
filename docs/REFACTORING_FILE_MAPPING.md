# Точная маппировка файлов: Seqeron.Genomics → Модульные пакеты

## Легенда
- ✅ = Файл готов к переносу
- 🔗 = Имеет зависимости от других файлов
- 📦 = Целевой пакет

---

## 📦 Seqeron.Genomics.Core (10 файлов)

| # | Файл | Строк | Зависит от | Используется в |
|---|------|-------|------------|----------------|
| 1 | ISequence.cs | 427 | - | DnaSequence, RnaSequence, ProteinSequence |
| 2 | DnaSequence.cs | 178 | ISequence, SequenceExtensions, SuffixTree | 14+ файлов |
| 3 | RnaSequence.cs | 168 | SequenceExtensions, SuffixTree | Translator |
| 4 | ProteinSequence.cs | 346 | SuffixTree | Translator, GenomeAnnotator |
| 5 | SequenceExtensions.cs | 219 | - | DnaSequence, RnaSequence, многие другие |
| 6 | GeneticCode.cs | 260 | - | Translator, GenomeAnnotator, CodonOptimizer |
| 7 | Translator.cs | 200 | GeneticCode, DnaSequence, RnaSequence, ProteinSequence | GenomeAnnotator |
| 8 | IupacHelper.cs | 34 | - | MotifFinder |
| 9 | ThermoConstants.cs | 102 | - | PrimerDesigner |
| 10 | StatisticsHelper.cs | 30 | - | Многие анализаторы |

### Namespace изменения для Core:
```csharp
// Было:
namespace Seqeron.Genomics;

// Станет:
namespace Seqeron.Genomics.Core;
```

---

## 📦 Seqeron.Genomics.IO (10 файлов)

| # | Файл | Строк | Зависит от Core | Используется в |
|---|------|-------|-----------------|----------------|
| 1 | FastaParser.cs | 156 | DnaSequence | - |
| 2 | FastqParser.cs | 449 | - | - |
| 3 | GenBankParser.cs | 546 | - | EmblParser, SequenceIO |
| 4 | EmblParser.cs | 524 | GenBankParser (records) | - |
| 5 | GffParser.cs | 471 | - | - |
| 6 | BedParser.cs | 573 | - | - |
| 7 | VcfParser.cs | 703 | - | VariantAnnotator |
| 8 | SequenceIO.cs | 864 | FastaParser, GenBankParser... | - |
| 9 | QualityScoreAnalyzer.cs | 523 | - | - |
| 10 | FeatureLocationHelper.cs | 63 | - | GenBankParser |

**Перенесено в другие пакеты:**
- ReportGenerator.cs → Seqeron.Genomics.Reports (отдельный пакет)
- CancellableOperations.cs → Seqeron.Genomics.Core (сквозная функциональность)

### Namespace изменения для IO:
```csharp
// Было:
namespace Seqeron.Genomics;

// Станет:
namespace Seqeron.Genomics.IO;

// Добавить using:
using Seqeron.Genomics.Core;
```

---

## 📦 Seqeron.Genomics.Alignment (3 файла)

| # | Файл | Строк | Зависит от Core | Используется в |
|---|------|-------|-----------------|----------------|
| 1 | SequenceAligner.cs | 523 | DnaSequence | VariantCaller, PhylogeneticAnalyzer |
| 2 | ApproximateMatcher.cs | 354 | - | - |
| 3 | SequenceAssembler.cs | 554 | DnaSequence | GenomeAssemblyAnalyzer |

### Namespace изменения для Alignment:
```csharp
namespace Seqeron.Genomics.Alignment;

using Seqeron.Genomics.Core;
```

---

## 📦 Seqeron.Genomics.Analysis (11 файлов)

| # | Файл | Строк | Зависит от | Используется в |
|---|------|-------|------------|----------------|
| 1 | KmerAnalyzer.cs | 351 | DnaSequence | Metagenomics, Analysis |
| 2 | SequenceStatistics.cs | 677 | DnaSequence, ProteinSequence | Многие |
| 3 | SequenceComplexity.cs | 383 | DnaSequence | - |
| 4 | GcSkewCalculator.cs | 296 | - | ChromosomeAnalyzer |
| 5 | GenomicAnalyzer.cs | 367 | DnaSequence | - |
| 6 | MotifFinder.cs | 527 | DnaSequence, IupacHelper | - |
| 7 | ProteinMotifFinder.cs | 728 | ProteinSequence | - |
| 8 | RepeatFinder.cs | 469 | DnaSequence | ChromosomeAnalyzer |
| 9 | RnaSecondaryStructure.cs | 606 | - | - |
| 10 | DisorderPredictor.cs | 504 | - | - |
| 11 | ComparativeGenomics.cs | 553 | DnaSequence | - |

### Namespace изменения для Analysis:
```csharp
namespace Seqeron.Genomics.Analysis;

using Seqeron.Genomics.Core;
using Seqeron.Genomics.Alignment; // если нужно
```

---

## 📦 Seqeron.Genomics.Annotation (8 файлов)

| # | Файл | Строк | Зависит от | Используется в |
|---|------|-------|------------|----------------|
| 1 | GenomeAnnotator.cs | 521 | GeneticCode, Translator, DnaSequence | - |
| 2 | VariantCaller.cs | 364 | DnaSequence, SequenceAligner | VariantAnnotator |
| 3 | VariantAnnotator.cs | 1070 | VcfParser, VariantCaller | - |
| 4 | StructuralVariantAnalyzer.cs | 668 | - | - |
| 5 | SpliceSitePredictor.cs | 673 | - | - |
| 6 | TranscriptomeAnalyzer.cs | 645 | - | - |
| 7 | EpigeneticsAnalyzer.cs | 640 | - | - |
| 8 | GenomeAssemblyAnalyzer.cs | 1009 | DnaSequence, SequenceAssembler | - |

### Namespace изменения для Annotation:
```csharp
namespace Seqeron.Genomics.Annotation;

using Seqeron.Genomics.Core;
using Seqeron.Genomics.IO;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;
```

---

## 📦 Seqeron.Genomics.Phylogenetics (1 файл)

| # | Файл | Строк | Зависит от | Используется в |
|---|------|-------|------------|----------------|
| 1 | PhylogeneticAnalyzer.cs | 654 | SequenceAligner | - |

### Namespace изменения:
```csharp
namespace Seqeron.Genomics.Phylogenetics;

using Seqeron.Genomics.Core;
using Seqeron.Genomics.Alignment;
```

---

## 📦 Seqeron.Genomics.Population (1 файл)

| # | Файл | Строк | Зависит от | Используется в |
|---|------|-------|------------|----------------|
| 1 | PopulationGeneticsAnalyzer.cs | 853 | - | - |

### Namespace изменения:
```csharp
namespace Seqeron.Genomics.Population;

using Seqeron.Genomics.Core;
```

---

## 📦 Seqeron.Genomics.Metagenomics (1 файл)

| # | Файл | Строк | Зависит от | Используется в |
|---|------|-------|------------|----------------|
| 1 | MetagenomicsAnalyzer.cs | 582 | KmerAnalyzer | - |

### Namespace изменения:
```csharp
namespace Seqeron.Genomics.Metagenomics;

using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;
```

---

## 📦 Seqeron.Genomics.MolTools (9 файлов)

| # | Файл | Строк | Зависит от | Используется в |
|---|------|-------|------------|----------------|
| 1 | CrisprDesigner.cs | 490 | DnaSequence | - |
| 2 | PrimerDesigner.cs | 492 | DnaSequence, ThermoConstants | - |
| 3 | ProbeDesigner.cs | 717 | DnaSequence | - |
| 4 | RestrictionAnalyzer.cs | 437 | DnaSequence | - |
| 5 | CodonOptimizer.cs | 584 | GeneticCode, DnaSequence, ProteinSequence | - |
| 6 | CodonUsageAnalyzer.cs | 493 | GeneticCode | - |
| 7 | MiRnaAnalyzer.cs | 542 | - | - |
| 8 | ChromosomeAnalyzer.cs | 746 | DnaSequence, RepeatFinder, GcSkewCalculator | - |
| 9 | PanGenomeAnalyzer.cs | 479 | - | - |

### Namespace изменения:
```csharp
namespace Seqeron.Genomics.MolTools;

using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis; // для ChromosomeAnalyzer
```

---

## 📊 Сводная таблица

| Пакет | Файлов | Строк | Уровень |
|-------|--------|-------|---------|
| Core | 10 | 1,964 | 0 |
| IO | 12 | 6,005 | 1 |
| Alignment | 3 | 1,431 | 1 |
| Analysis | 11 | 5,461 | 2 |
| Annotation | 8 | 5,590 | 2 |
| Phylogenetics | 1 | 654 | 2 |
| Population | 1 | 853 | 2 |
| Metagenomics | 1 | 582 | 2 |
| MolTools | 9 | 4,980 | 2 |
| **ИТОГО** | **56** | **27,520** | - |

---

## 🔍 Проверочный чеклист

### Перед началом
- [ ] Убедиться: 56 файлов .cs в Seqeron.Genomics
- [ ] Записать: количество тестов (baseline)
- [ ] Записать: покрытие кода (baseline)

### После завершения
- [ ] Проверить: 56 файлов .cs распределены по пакетам
- [ ] Проверить: количество тестов не изменилось
- [ ] Проверить: все тесты проходят
- [ ] Проверить: покрытие не снизилось

---

## 📝 Шаблон .csproj для каждого пакета

### Seqeron.Genomics.Core.csproj
```xml
<Project Sdk="Microsoft.NET.Sdk">
  <PropertyGroup>
    <TargetFramework>net8.0</TargetFramework>
    <ImplicitUsings>enable</ImplicitUsings>
    <Nullable>enable</Nullable>
    <RootNamespace>Seqeron.Genomics.Core</RootNamespace>
  </PropertyGroup>

  <ItemGroup>
    <ProjectReference Include="..\..\SuffixTree\SuffixTree\SuffixTree.csproj" />
  </ItemGroup>
</Project>
```

### Seqeron.Genomics.IO.csproj
```xml
<Project Sdk="Microsoft.NET.Sdk">
  <PropertyGroup>
    <TargetFramework>net8.0</TargetFramework>
    <ImplicitUsings>enable</ImplicitUsings>
    <Nullable>enable</Nullable>
    <RootNamespace>Seqeron.Genomics.IO</RootNamespace>
  </PropertyGroup>

  <ItemGroup>
    <ProjectReference Include="..\Seqeron.Genomics.Core\Seqeron.Genomics.Core.csproj" />
  </ItemGroup>
</Project>
```

### Seqeron.Genomics.Alignment.csproj
```xml
<Project Sdk="Microsoft.NET.Sdk">
  <PropertyGroup>
    <TargetFramework>net8.0</TargetFramework>
    <ImplicitUsings>enable</ImplicitUsings>
    <Nullable>enable</Nullable>
    <RootNamespace>Seqeron.Genomics.Alignment</RootNamespace>
  </PropertyGroup>

  <ItemGroup>
    <ProjectReference Include="..\Seqeron.Genomics.Core\Seqeron.Genomics.Core.csproj" />
  </ItemGroup>
</Project>
```

### Seqeron.Genomics.Analysis.csproj
```xml
<Project Sdk="Microsoft.NET.Sdk">
  <PropertyGroup>
    <TargetFramework>net8.0</TargetFramework>
    <ImplicitUsings>enable</ImplicitUsings>
    <Nullable>enable</Nullable>
    <RootNamespace>Seqeron.Genomics.Analysis</RootNamespace>
  </PropertyGroup>

  <ItemGroup>
    <ProjectReference Include="..\Seqeron.Genomics.Core\Seqeron.Genomics.Core.csproj" />
    <ProjectReference Include="..\Seqeron.Genomics.Alignment\Seqeron.Genomics.Alignment.csproj" />
  </ItemGroup>
</Project>
```

### (Аналогично для остальных пакетов...)

---

## 🔄 Type Forwarding для обратной совместимости

В мета-пакете Seqeron.Genomics создать файл `TypeForwards.cs`:

```csharp
// Seqeron.Genomics/TypeForwards.cs
using System.Runtime.CompilerServices;

// Core types
[assembly: TypeForwardedTo(typeof(Seqeron.Genomics.Core.DnaSequence))]
[assembly: TypeForwardedTo(typeof(Seqeron.Genomics.Core.RnaSequence))]
[assembly: TypeForwardedTo(typeof(Seqeron.Genomics.Core.ProteinSequence))]
[assembly: TypeForwardedTo(typeof(Seqeron.Genomics.Core.ISequence))]
// ... и так далее для всех публичных типов
```

Это позволит существующему коду, использующему `Seqeron.Genomics.DnaSequence`, продолжать работать.
