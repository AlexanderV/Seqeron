# Seqeron.Genomics Refactoring — Quick Reference

## 🎯 Цель
Разделить монолитный `Seqeron.Genomics` (56 файлов) на 12 модульных пакетов.

## 📦 Структура пакетов

```
Уровень 0:  Infrastructure (3 файла)     — общие типы, константы
                    ↓
Уровень 1:  Core (8 файлов)              — DnaSequence, RnaSequence, ProteinSequence
                    ↓
Уровень 2:  IO (10)  |  Alignment (4)    — парсеры | выравнивание
                    ↓
Уровень 3:  Analysis (11) | Annotation (8) | Phylogenetics (1) | Population (1) | Metagenomics (1) | MolTools (9)
                    ↓
Уровень 4:  Reports (1)                  — генерация отчётов
                    ↓
Уровень 5:  Seqeron.Genomics (мета)      — обратная совместимость
```

## 📊 Сводка

| Пакет | Файлов | Ключевые типы |
|-------|--------|---------------|
| Infrastructure | 3 | ScoringMatrix, AlignmentResult, StatisticsHelper, ThermoConstants |
| Core | 8 | DnaSequence, RnaSequence, ProteinSequence, GeneticCode, Translator |
| IO | 10 | FastaParser, GenBankParser, VcfParser, BedParser, GffParser |
| Alignment | 4 | SequenceAligner, ApproximateMatcher, CancellableOperations |
| Analysis | 11 | KmerAnalyzer, MotifFinder, RepeatFinder, SequenceComplexity |
| Annotation | 8 | GenomeAnnotator, VariantCaller, VariantAnnotator |
| Phylogenetics | 1 | PhylogeneticAnalyzer |
| Population | 1 | PopulationGeneticsAnalyzer |
| Metagenomics | 1 | MetagenomicsAnalyzer |
| MolTools | 9 | CrisprDesigner, PrimerDesigner, RestrictionAnalyzer |
| Reports | 1 | ReportGenerator |

## ⚡ Ключевые решения

1. **Infrastructure на нижнем уровне** — содержит общие типы (ScoringMatrix, AlignmentResult), чтобы избежать циклических зависимостей

2. **CancellableOperations в Alignment** — содержит реализации алгоритмов, а не абстракции

3. **Reports на верхнем уровне** — зависит от многих пакетов

4. **Обратная совместимость** — мета-пакет с TypeForwarding

## 📋 Фазы выполнения

| Фаза | Действие | Файлов |
|------|----------|--------|
| 0 | Подготовка (git branch, baseline tests) | - |
| 1 | Создать Infrastructure | 3 |
| 2 | Создать Core | 8 |
| 3 | Создать IO | 10 |
| 4 | Создать Alignment | 4 |
| 5 | Создать Analysis | 11 |
| 6 | Создать остальные пакеты | 21 |
| 7 | Создать мета-пакет | 1 |
| 8 | Мигрировать тесты (using) | - |
| 9 | Финальное тестирование | - |

## 🔧 Команды проверки

```powershell
# Перед началом
dotnet build
dotnet test

# После каждой фазы
dotnet build Seqeron.sln

# Финальная проверка
dotnet test --collect:"XPlat Code Coverage"
```

## 📁 Документация

- [Детальный план](REFACTORING_PLAN.md)
- [Маппинг файлов](REFACTORING_FILE_MAPPING.md) 
- [Infrastructure детали](REFACTORING_INFRASTRUCTURE.md)
- [Полная сводка](REFACTORING_SUMMARY.md)
