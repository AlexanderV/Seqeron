# Отчёт: Анализ Clean Code для SuffixTree.Genomics

**Дата:** 2026-01-22
**Версия:** 1.0

## Общая оценка

Код **хорошего качества**, но имеет ряд нарушений принципов Clean Code. Основные проблемы: чрезмерно длинные методы, дублирование кода, магические числа и смешение уровней абстракции.

---

## 🔴 КРИТИЧЕСКИЕ нарушения

### 1. God Methods — методы делают слишком много

**Файлы:** `ChromosomeAnalyzer.cs:623-696`, `MetagenomicsAnalyzer.cs:421-478`

**Проблема:** Методы превышают 30 строк и делают несколько вещей:

```csharp
// ❌ FindSyntenyBlocks - 74 строки, множество ответственностей
public static IEnumerable<SyntenyBlock> FindSyntenyBlocks(...)
{
    var pairs = orthologPairs.ToList();
    // 1. Валидация
    // 2. Группировка
    // 3. Сортировка
    // 4. Поиск collinear runs
    // 5. Создание результата
    // ... 74 строки кода
}
```

**Рекомендация:** Разбить на методы по 5-20 строк:
```csharp
public static IEnumerable<SyntenyBlock> FindSyntenyBlocks(...)
{
    var pairs = ValidateAndPrepare(orthologPairs);
    var chromPairs = GroupByChromosomePairs(pairs);

    foreach (var group in chromPairs)
    {
        foreach (var block in FindCollinearBlocks(group, minGenes, maxGap))
            yield return block;
    }
}
```

---

### 2. Дублирование логики GC Content (DRY violation)

**Файлы:** Минимум 6 разных реализаций:
- `SequenceExtensions.cs:21-34`
- `PrimerDesigner.cs:235-243`
- `MetagenomicsAnalyzer.cs:480-486`
- `ChromosomeAnalyzer.cs:522-527`
- `DnaSequence.cs:96-102`
- `GenomeAnnotator.cs:489-496`

**Проблема:** Идентичная логика вычисления GC content повторяется:

```csharp
// Версия 1
int gc = sequence.Count(c => c == 'G' || c == 'C' || c == 'g' || c == 'c');
return (double)gc / sequence.Length;

// Версия 2
int gcCount = seq.Count(c => c == 'G' || c == 'C');
return (double)gcCount / seq.Length * 100;

// Версия 3...
```

**Рекомендация:** Единый helper-метод в `SequenceExtensions`:
```csharp
public static double CalculateGcContent(this ReadOnlySpan<char> sequence, bool asPercentage = true)
{
    // Единственная реализация
}
```

---

### 3. Дублирование логики Complement/ReverseComplement

**Файлы:**
- `SequenceExtensions.cs:53-94`
- `DnaSequence.cs:54-89`
- `DnaSequence.cs:169-188`
- `GenomeAnnotator.cs:68-80`

**Проблема:** Mapping `A↔T, G↔C` повторяется в 4+ местах:

```csharp
// GenomeAnnotator
private static readonly Dictionary<char, char> ComplementMap = new()
{
    ['A'] = 'T', ['T'] = 'A', ['C'] = 'G', ['G'] = 'C', ...
};

// SequenceExtensions
destination[i] = sequence[i] switch
{
    'A' or 'a' => 'T',
    'T' or 't' => 'A', ...
};

// DnaSequence
sb.Append(c switch { 'A' => 'T', 'T' => 'A', ... });
```

**Рекомендация:** Единый helper:
```csharp
public static class NucleotideHelper
{
    public static char GetComplement(char nucleotide) => nucleotide switch { ... };
}
```

---

### 4. Magic Numbers — магические константы

**Файлы:** Множество файлов

**Проблема:** Числа без объяснения смысла:

```csharp
// ChromosomeAnalyzer.cs:220
ploidy = Math.Max(1, Math.Min(8, ploidy)); // Почему 8?

// PrimerDesigner.cs:211
double tm = 64.9 + 41.0 * (gcCount - 16.4) / seq.Length; // Откуда эти числа?

// GenomeAnnotator.cs:501
return validRatio * 0.7 + Math.Abs(gc3ratio - 0.5) * 0.6; // Что это?

// ChromosomeAnalyzer.cs:310
if (similarity >= 0.7) // Magic threshold
```

**Рекомендация:** Именованные константы:
```csharp
private const int MaxReasonablePloidyLevel = 8;
private const double WallaceRuleTmConstant = 64.9;
private const double WallaceRuleGcFactor = 41.0;
private const double TelomereMinSimilarityThreshold = 0.7;
```

---

## 🟠 ВЫСОКИЕ нарушения

### 5. Слишком много параметров в методах (> 3)

**Файлы:**

```csharp
// RepeatFinder.cs:48-59 - 6 параметров
public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
    DnaSequence sequence,
    int minUnitLength,
    int maxUnitLength,
    int minRepeats,
    CancellationToken cancellationToken,
    IProgress<double>? progress = null)

// ChromosomeAnalyzer.cs:235-242 - 6 параметров
public static TelomereResult AnalyzeTelomeres(
    string chromosomeName,
    string sequence,
    string telomereRepeat = "TTAGGG",
    int searchLength = 10000,
    int minTelomereLength = 500,
    int criticalLength = 3000)

// GenomeAnnotator.cs:90-94 - 4 параметра
public static IEnumerable<OpenReadingFrame> FindOrfs(
    string dnaSequence,
    int minLength = 100,
    bool searchBothStrands = true,
    bool requireStartCodon = true)
```

**Рекомендация:** Использовать Parameter Objects:
```csharp
public readonly record struct TelomereAnalysisOptions(
    string TelomereRepeat = "TTAGGG",
    int SearchLength = 10000,
    int MinTelomereLength = 500,
    int CriticalLength = 3000);

public static TelomereResult AnalyzeTelomeres(
    string chromosomeName,
    string sequence,
    TelomereAnalysisOptions? options = null)
```

---

### 6. Нарушение Command-Query Separation

**Файлы:** `MetagenomicsAnalyzer.cs:575-606`

**Проблема:** Метод возвращает данные И модифицирует состояние (вычисляет несколько метрик одновременно):

```csharp
// ❌ Возвращает tuple с 3 разными вещами
public static (double FunctionalRichness, double FunctionalDiversity, IReadOnlyDictionary<string, int> PathwayCounts)
    CalculateFunctionalDiversity(IEnumerable<FunctionalAnnotation> annotations)
```

**Рекомендация:** Разбить на отдельные методы:
```csharp
public static double CalculateFunctionalRichness(...)
public static double CalculateFunctionalDiversity(...)
public static IReadOnlyDictionary<string, int> GetPathwayCounts(...)
```

---

### 7. Глубокая вложенность (> 2 уровня)

**Файлы:** `GenomeAnnotator.cs:432-455`, `ChromosomeAnalyzer.cs:637-695`

**Проблема:** 3-4 уровня вложенности:

```csharp
// GenomeAnnotator.cs:432-455
foreach (string motif in minus35Motifs)  // Level 1
{
    for (int i = 0; i <= seq.Length - motif.Length; i++)  // Level 2
    {
        if (seq.Substring(i, motif.Length) == motif)  // Level 3
        {
            // ...
        }
    }
}
```

**Рекомендация:** Early return, extract method:
```csharp
foreach (string motif in minus35Motifs)
{
    foreach (var match in FindMotifMatches(seq, motif))
        yield return ("-35 box", match.Position, match.Score);
}
```

---

### 8. Смешение уровней абстракции

**Файлы:** `PrimerDesigner.cs:39-114`

**Проблема:** Высокоуровневая логика смешана с низкоуровневыми деталями:

```csharp
public static PrimerPairResult DesignPrimers(...)
{
    // Высокий уровень: дизайн праймеров
    // Средний уровень: поиск кандидатов
    // Низкий уровень: substring, loop индексы

    for (int start = forwardSearchStart; start < forwardSearchEnd; start++)  // Низкий
    {
        for (int len = param.MinLength; len <= param.MaxLength && start + len <= targetStart; len++)
        {
            var candidate = EvaluatePrimer(template.Sequence.Substring(start, len), ...);  // Смесь
        }
    }
}
```

**Рекомендация:** Отдельные методы для каждого уровня:
```csharp
public static PrimerPairResult DesignPrimers(...)
{
    var forwardCandidates = FindForwardPrimerCandidates(template, targetStart, param);
    var reverseCandidates = FindReversePrimerCandidates(template, targetEnd, param);
    return SelectBestPrimerPair(forwardCandidates, reverseCandidates);
}
```

---

## 🟡 СРЕДНИЕ нарушения

### 9. Несогласованное именование (Inconsistent Naming)

**Проблемы:**

| Файл | Проблема | Пример |
|------|----------|--------|
| MetagenomicsAnalyzer | Сокращения | `IncrementCount` vs `CalculateBrayCurtis` |
| ChromosomeAnalyzer | Аббревиатуры | `gcValues`, `centMid`, `qArmLength` |
| PrimerDesigner | Hungarian notation | `seq`, `param`, `len` |
| Various | Inconsistent casing | `_enzymes`, `_codonTable`, `_sequence` |

**Рекомендация:** Единый стиль — полные слова:
```csharp
// ❌
int centMid = (centStart.Value + centEnd.Value) / 2;

// ✅
int centromereMiddlePosition = (centromereStart.Value + centromereEnd.Value) / 2;
```

---

### 10. Boolean параметры без контекста

**Файлы:**

```csharp
// GenomeAnnotator.cs:91-92
bool searchBothStrands = true,
bool requireStartCodon = true

// PrimerDesigner.cs:406
bool forward = true

// ChromosomeAnalyzer.cs:281
bool fromEnd
```

**Рекомендация:** Enum или named arguments:
```csharp
public enum StrandSearchMode { Forward, Both }
public enum CodonRequirement { RequireStart, AllowAny }

public static IEnumerable<OpenReadingFrame> FindOrfs(
    string dnaSequence,
    int minLength = 100,
    StrandSearchMode strandMode = StrandSearchMode.Both,
    CodonRequirement codonMode = CodonRequirement.RequireStart)
```

---

### 11. Избыточные комментарии (Noise Comments)

**Файлы:** Многие файлы

```csharp
// ❌ Комментарий повторяет код
/// <summary>
/// Gets the name of this genetic code.
/// </summary>
public string Name { get; }

// ❌ Очевидный комментарий
// Search forward strand
foreach (var orf in FindOrfsInStrand(dnaSequence, ...))

// ❌ TODO без тикета
// Simplified - would use pre-computed coding/non-coding hexamer tables
```

**Рекомендация:** Удалить очевидные комментарии, оставить только WHY:
```csharp
// Hexamer bias indicates coding potential - AT-rich regions typically non-coding
// Reference: Fickett 1982 "Recognition of protein coding regions"
```

---

### 12. Длинные файлы (> 500 строк)

| Файл | Строк | Рекомендация |
|------|-------|--------------|
| ChromosomeAnalyzer.cs | 892 | Разбить на: TelomereAnalyzer, CentromereAnalyzer, SyntenyAnalyzer |
| MetagenomicsAnalyzer.cs | 706 | Разбить на: TaxonomicClassifier, DiversityCalculator, GenomeBinner |
| PrimerDesigner.cs | 518 | Выделить: TmCalculator, HairpinDetector |
| RestrictionAnalyzer.cs | 516 | Выделить: EnzymeDatabase, DigestSimulator |

---

## 🟢 НИЗКИЕ нарушения / рекомендации

### 13. Использование `var` непоследовательно

```csharp
// Иногда используется
var pairs = orthologPairs.ToList();

// Иногда явный тип
Dictionary<string, int> counts = new Dictionary<string, int>();
```

**Рекомендация:** Использовать `var` когда тип очевиден справа.

---

### 14. Отсутствие Null Object Pattern

```csharp
// Возврат null при ошибке
if (taxonCounts.Count == 0)
{
    yield return new TaxonomicClassification(
        id, "Unclassified", "", "", "", "", "", "", 0, 0, totalKmers);
}
```

**Рекомендация:** Использовать `TaxonomicClassification.Unclassified` как Null Object.

---

### 15. String concatenation в циклах

```csharp
// GenomeAnnotator.cs:571
string repeatSeq = sequence.Substring(start, end - start);
yield return (start, end, "inverted_repeat", arm1 + "..." + arm2);
```

**Рекомендация:** StringBuilder или string interpolation.

---

## Сводная таблица

| # | Нарушение | Уровень | Принцип | Файлы затронуты |
|---|-----------|---------|---------|-----------------|
| 1 | God Methods (>30 строк) | 🔴 Критический | Functions | ChromosomeAnalyzer, MetagenomicsAnalyzer |
| 2 | Дублирование GC Content | 🔴 Критический | DRY | 6+ файлов |
| 3 | Дублирование Complement | 🔴 Критический | DRY | 4+ файлов |
| 4 | Magic Numbers | 🔴 Критический | Meaningful Names | Везде |
| 5 | >3 параметра в методах | 🟠 Высокий | Functions | RepeatFinder, ChromosomeAnalyzer |
| 6 | Command-Query violation | 🟠 Высокий | Functions | MetagenomicsAnalyzer |
| 7 | Глубокая вложенность | 🟠 Высокий | Formatting | GenomeAnnotator, ChromosomeAnalyzer |
| 8 | Смешение абстракций | 🟠 Высокий | Functions | PrimerDesigner |
| 9 | Inconsistent Naming | 🟡 Средний | Naming | Везде |
| 10 | Boolean параметры | 🟡 Средний | Naming | GenomeAnnotator, PrimerDesigner |
| 11 | Noise Comments | 🟡 Средний | Comments | Везде |
| 12 | Файлы >500 строк | 🟡 Средний | Formatting | 4 файла |
| 13 | Непоследовательный var | 🟢 Низкий | Formatting | Везде |
| 14 | Нет Null Object | 🟢 Низкий | Error Handling | MetagenomicsAnalyzer |
| 15 | String concat в циклах | 🟢 Низкий | Performance | GenomeAnnotator |

---

## Положительные аспекты ✅

1. **XML Documentation** — хорошие XML docs на публичных API
2. **Immutable Records** — использование `readonly record struct` для result types
3. **Nullable Reference Types** — корректное использование `?` и null checks
4. **ArgumentNullException.ThrowIfNull** — современный guard clause
5. **Expression-bodied members** — краткий синтаксис где уместно
6. **Pattern matching** — хорошее использование `switch expressions`
7. **Span<T>** — zero-allocation API в SequenceExtensions
8. **Regions** — логическая группировка кода (#region)

---

## Метрики качества (оценка)

| Метрика | Текущее | Цель |
|---------|---------|------|
| Средний размер метода | ~25 строк | <20 строк |
| Макс. параметров | 6 | ≤3 |
| Макс. вложенность | 4 уровня | ≤2 уровня |
| Дублирование кода | ~15% | <5% |
| Размер файлов | до 892 строк | <500 строк |

---

## Рекомендуемый план рефакторинга

### Фаза 1 (Критические)
1. Создать `NucleotideHelper` с единой логикой complement/GC
2. Разбить God Methods на 5-20 строк
3. Извлечь магические числа в именованные константы

### Фаза 2 (Высокие)
4. Создать Parameter Objects для методов с >3 параметрами
5. Разделить Command и Query
6. Уменьшить вложенность через early return

### Фаза 3 (Средние)
7. Разбить большие файлы (>500 строк)
8. Унифицировать именование
9. Заменить boolean параметры на enum

---

## Ссылки

- [Clean Code by Robert C. Martin](https://www.oreilly.com/library/view/clean-code-a/9780136083238/)
- [Refactoring by Martin Fowler](https://refactoring.com/)
