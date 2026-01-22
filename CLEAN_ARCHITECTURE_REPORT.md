# Отчёт: Анализ Clean Architecture для SuffixTree.Genomics

**Дата:** 2026-01-22
**Версия:** 1.0

## Общая оценка

Библиотека **частично соответствует** принципам Clean Architecture. Есть сильные стороны (разделение на проекты, immutable types, хорошая документация), но также присутствуют нарушения, требующие внимания.

---

## 🔴 КРИТИЧЕСКИЕ нарушения

### 1. Anemic Domain Model + God Static Classes

**Файлы:** `GenomicAnalyzer.cs`, `CrisprDesigner.cs`, `RepeatFinder.cs`, и др. (45+ классов)

**Проблема:** Вся бизнес-логика размещена в static классах, а доменные модели (`DnaSequence`, `RnaSequence`) — это просто контейнеры данных.

```csharp
// ❌ СЕЙЧАС: Anemic Domain Model
public sealed class DnaSequence
{
    private readonly string _sequence;
    public string Sequence => _sequence;
    public int Length => _sequence.Length;
    // Только базовые операции, основная логика - в static классах
}

// Static класс с бизнес-логикой
public static class GenomicAnalyzer
{
    public static RepeatInfo FindLongestRepeat(DnaSequence sequence) { ... }
    public static IEnumerable<TandemRepeat> FindTandemRepeats(DnaSequence sequence) { ... }
}
```

**Рекомендация:**
```csharp
// ✅ Rich Domain Model
public sealed class DnaSequence
{
    public RepeatInfo FindLongestRepeat() => _repeatAnalyzer.FindLongest(this);
    public IEnumerable<TandemRepeat> FindTandemRepeats() => _repeatAnalyzer.Find(this);
}

// Или через Extension методы для обратной совместимости
public static class DnaSequenceExtensions
{
    public static RepeatInfo FindLongestRepeat(this DnaSequence seq) { ... }
}
```

**Влияние:** Нарушает SRP (Single Responsibility) и OCP (Open/Closed). Затрудняет тестирование и расширение.

---

### 2. DnaSequence не реализует ISequence

**Файлы:** `DnaSequence.cs`, `ISequence.cs`

**Проблема:** Существует интерфейс `ISequence`, но основной класс `DnaSequence` его **не реализует**:

```csharp
// ISequence существует
public interface ISequence
{
    string Sequence { get; }
    int Length { get; }
    SequenceType Type { get; }
    IReadOnlySet<char> Alphabet { get; }
    // ...
}

// ❌ DnaSequence НЕ реализует ISequence
public sealed class DnaSequence  // <-- Нет : ISequence
{
    public string Sequence => _sequence;
    // ...
}
```

**Рекомендация:** `DnaSequence` должен реализовать `ISequence` для полиморфизма и DIP (Dependency Inversion).

**Влияние:** Нарушает LSP (Liskov Substitution) и ISP (Interface Segregation). Невозможно работать с разными типами последовательностей единообразно.

---

### 3. Result Types в одном файле с логикой

**Файлы:** `GenomicAnalyzer.cs:329-428`, `RepeatFinder.cs:418-498`, `CrisprDesigner.cs:499-593`

**Проблема:** Result types (`RepeatInfo`, `TandemRepeat`, `MicrosatelliteResult`, `CrisprSystem`, и др.) объявлены в одном файле с логикой анализаторов.

**Рекомендация:** Вынести в отдельные файлы/папку:
```
SuffixTree.Genomics/
├── Models/
│   ├── RepeatInfo.cs
│   ├── TandemRepeat.cs
│   ├── MicrosatelliteResult.cs
│   └── CrisprTypes.cs
```

**Влияние:** Нарушает SRP. Затрудняет навигацию и переиспользование.

---

## 🟠 ВЫСОКИЕ нарушения

### 4. Смешение Infrastructure и Domain в FastaParser

**Файл:** `FastaParser.cs`

**Проблема:** Parser напрямую создаёт `DnaSequence` — жёсткая связь с конкретной доменной моделью:

```csharp
// ❌ Жёсткая связь с DnaSequence
private static FastaEntry CreateEntry(string header, string sequence)
{
    return new FastaEntry(id, description, new DnaSequence(sequence));  // Hardcoded
}
```

**Рекомендация:**
```csharp
// ✅ Использовать фабрику или generics
public static class FastaParser<TSequence> where TSequence : ISequence
{
    public static IEnumerable<FastaEntry<TSequence>> Parse(
        string content,
        Func<string, TSequence> sequenceFactory) { ... }
}
```

**Влияние:** Нарушает OCP и DIP. Невозможно парсить в RnaSequence или ProteinSequence без дублирования кода.

---

### 5. Отсутствие интерфейсов для анализаторов

**Файлы:** Все анализаторы (45+ классов)

**Проблема:** Анализаторы — это static классы без интерфейсов:

```csharp
// ❌ Static класс, нет интерфейса
public static class CrisprDesigner
{
    public static IEnumerable<PamSite> FindPamSites(...) { ... }
}
```

**Рекомендация:**
```csharp
// ✅ Интерфейс + реализация
public interface ICrisprDesigner
{
    IEnumerable<PamSite> FindPamSites(ISequence sequence, CrisprSystemType type);
}

public class CrisprDesigner : ICrisprDesigner
{
    public IEnumerable<PamSite> FindPamSites(...) { ... }
}
```

**Влияние:** Нарушает DIP и затрудняет тестирование (невозможно mock-ить static классы).

---

### 6. Дублирование логики между DnaSequence и static helpers

**Файлы:** `DnaSequence.cs:54-89`, `DnaSequence.cs:169-188`

**Проблема:** Логика Complement/ReverseComplement дублируется:

```csharp
// В DnaSequence
public DnaSequence ReverseComplement() { ... }

// И как static helper
public static string GetReverseComplementString(string sequence) { ... }
```

**Рекомендация:** Оставить только метод экземпляра, static helper сделать internal или удалить.

---

## 🟡 СРЕДНИЕ нарушения

### 7. Отсутствие слоя Application (Use Cases)

**Проблема:** Нет явного разделения между:
- **Domain Logic** (правила биоинформатики)
- **Application Logic** (оркестрация use cases)

Все анализаторы смешивают обе ответственности.

**Рекомендация:** Создать слой Use Cases:
```csharp
// Application/UseCases/AnalyzeGenomeUseCase.cs
public class AnalyzeGenomeUseCase
{
    public GenomeAnalysisResult Execute(GenomeAnalysisRequest request)
    {
        var sequence = _sequenceRepository.Get(request.SequenceId);
        var repeats = _repeatFinder.FindMicrosatellites(sequence);
        var crisprSites = _crisprDesigner.FindPamSites(sequence);
        return new GenomeAnalysisResult(repeats, crisprSites);
    }
}
```

---

### 8. Нет Repository Pattern для последовательностей

**Проблема:** Отсутствует абстракция доступа к данным последовательностей.

**Рекомендация:**
```csharp
// Domain/Repositories/ISequenceRepository.cs
public interface ISequenceRepository
{
    Task<ISequence> GetByIdAsync(string id);
    Task<IEnumerable<ISequence>> SearchAsync(SequenceQuery query);
}
```

---

### 9. Enum в Domain без поведения

**Файлы:** `RepeatFinder.cs:421-430`, `CrisprDesigner.cs:502-518`

**Проблема:** Enum-ы используются как примитивы:

```csharp
public enum CrisprSystemType { SpCas9, SaCas9, Cas12a, ... }

// Switch по enum в каждом методе
public static CrisprSystem GetSystem(CrisprSystemType type) => type switch { ... };
```

**Рекомендация:** Использовать Smart Enum / Value Object pattern:
```csharp
public sealed class CrisprSystemType : Enumeration
{
    public static readonly CrisprSystemType SpCas9 = new(1, "SpCas9", "NGG", 20);
    public string PamSequence { get; }
    public int GuideLength { get; }
}
```

---

## 🟢 НИЗКИЕ нарушения / рекомендации

### 10. Отсутствие Value Objects для примитивов

**Пример:** Позиции, длины, scores — это примитивы int/double.

**Рекомендация:**
```csharp
public readonly record struct Position(int Value)
{
    public static Position operator +(Position a, Position b) => new(a.Value + b.Value);
}

public readonly record struct GcContent(double Percentage)
{
    public bool IsOptimal => Percentage >= 40 && Percentage <= 60;
}
```

---

### 11. Нет Domain Events

**Проблема:** При обнаружении важных паттернов (CRISPR site, repeat) нет механизма событий.

**Рекомендация:**
```csharp
public interface IDomainEvent { DateTime OccurredOn { get; } }
public record CrisprSiteFoundEvent(PamSite Site) : IDomainEvent;
```

---

### 12. Отсутствие Specification Pattern

**Проблема:** Критерии фильтрации зашиты в методах.

**Рекомендация:**
```csharp
public interface ISpecification<T>
{
    bool IsSatisfiedBy(T entity);
}

public class MinLengthSpecification : ISpecification<MicrosatelliteResult>
{
    public MinLengthSpecification(int minLength) { ... }
    public bool IsSatisfiedBy(MicrosatelliteResult result) => result.TotalLength >= _minLength;
}
```

---

## Сводная таблица

| # | Нарушение | Уровень | SOLID | Clean Architecture |
|---|-----------|---------|-------|-------------------|
| 1 | Anemic Domain + God Classes | 🔴 Критический | SRP, OCP | Domain Layer |
| 2 | DnaSequence ≠ ISequence | 🔴 Критический | LSP, DIP | Domain Layer |
| 3 | Result Types в одном файле | 🔴 Критический | SRP | Domain Layer |
| 4 | Parser ↔ DnaSequence coupling | 🟠 Высокий | OCP, DIP | Infrastructure |
| 5 | Static анализаторы без интерфейсов | 🟠 Высокий | DIP | Application Layer |
| 6 | Дублирование логики | 🟠 Высокий | DRY | Domain Layer |
| 7 | Нет Application Layer | 🟡 Средний | - | Application Layer |
| 8 | Нет Repository Pattern | 🟡 Средний | DIP | Infrastructure |
| 9 | Enum без поведения | 🟡 Средний | OCP | Domain Layer |
| 10 | Примитивы вместо Value Objects | 🟢 Низкий | - | Domain Layer |
| 11 | Нет Domain Events | 🟢 Низкий | - | Domain Layer |
| 12 | Нет Specification Pattern | 🟢 Низкий | OCP | Domain Layer |

---

## Положительные аспекты ✅

1. **Immutable Result Types** — использование `readonly record struct`
2. **Хорошая документация** — XML docs на всех публичных API
3. **Разделение проектов** — Core отделён от Genomics
4. **Nullable reference types** — включены
5. **Zero-allocation overloads** — Span<T> перегрузки для производительности
6. **Код качества** — TreatWarningsAsErrors, Code Analysis

---

## Рекомендуемый план рефакторинга

### Фаза 1 (Критические)
1. DnaSequence реализует ISequence
2. Вынести Result Types в отдельные файлы
3. Перенести логику из static классов в domain models (или Extension methods)

### Фаза 2 (Высокие)
4. Создать интерфейсы для анализаторов
5. Сделать parsers generic/configurable
6. Устранить дублирование

### Фаза 3 (Средние/Низкие)
7. Добавить Application Layer с Use Cases
8. Реализовать Repository pattern
9. Smart Enums и Value Objects

---

## Архитектурная диаграмма (текущая)

```
┌─────────────────────────────────────────────────────────────┐
│                    Console / Benchmarks                      │
│              (Presentation/Demo Layer)                       │
└───────────────────────┬─────────────────────────────────────┘
                        │
┌───────────────────────▼─────────────────────────────────────┐
│               SuffixTree.Genomics                            │
│          (Application/Use Case Layer)                        │
│  45+ Static Analyzer Classes (Domain Logic)                 │
├─────────────────────────────────────────────────────────────┤
│ • Pattern Matching (MotifFinder, ApproximateMatcher)        │
│ • Repeat Analysis (RepeatFinder, GenomeAssemblyAnalyzer)    │
│ • Population/Comparative (PhylogeneticAnalyzer, etc.)       │
│ • Codon Analysis (CodonOptimizer, CodonUsageAnalyzer)       │
│ • Molecular Tools (CrisprDesigner, RestrictionAnalyzer)     │
│ • Format Parsers (FASTA, FASTQ, GenBank, GFF, BED, VCF)     │
└───────────────────────┬─────────────────────────────────────┘
                        │
┌───────────────────────▼─────────────────────────────────────┐
│          SuffixTree.Genomics Domain Models                   │
│              (Domain Layer)                                  │
├─────────────────────────────────────────────────────────────┤
│ • ISequence (базовый интерфейс)                             │
│ • DnaSequence, RnaSequence, ProteinSequence                 │
│ • GeneticCode (таблицы кодонов)                             │
│ • Result types (MicrosatelliteResult, RepeatInfo, etc.)     │
└───────────────────────┬─────────────────────────────────────┘
                        │
┌───────────────────────▼─────────────────────────────────────┐
│                   SuffixTree Core                            │
│         (Data Structure / Algorithms)                        │
├─────────────────────────────────────────────────────────────┤
│ • ISuffixTree interface                                      │
│ • SuffixTree implementation (Ukkonen's O(n) algorithm)      │
│ • SuffixTreeNode (hybrid children storage)                  │
│ • Search algorithms (O(m), O(m+k))                          │
│ • Pattern algorithms (LRS, LCS)                             │
└─────────────────────────────────────────────────────────────┘
```

---

## Ссылки

- [Clean Architecture by Robert C. Martin](https://blog.cleancoder.com/uncle-bob/2012/08/13/the-clean-architecture.html)
- [Domain-Driven Design by Eric Evans](https://www.domainlanguage.com/ddd/)
- [SOLID Principles](https://en.wikipedia.org/wiki/SOLID)
