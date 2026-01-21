# C# Clean Code Examples

Practical C# code examples demonstrating Clean Code principles.

## Structure

### Core Examples

| File | Principle | Domain | Description |
|------|-----------|--------|-------------|
| [01-naming.md](01-naming.md) | Meaningful Names | Gaming | Variable, method, class naming conventions |
| [02-functions.md](02-functions.md) | Functions | IoT | Small functions, SRP, parameters |
| [03-error-handling.md](03-error-handling.md) | Error Handling | Banking | Result pattern, exceptions, IExceptionHandler |
| [04-objects-vs-data.md](04-objects-vs-data.md) | Objects vs Data | Healthcare | Rich models vs DTOs, Law of Demeter |
| [05-refactoring.md](05-refactoring.md) | Refactoring | Education | Code smells and refactoring patterns |
| [06-testing.md](06-testing.md) | Unit Tests | Social Media | F.I.R.S.T., AAA pattern, test naming |
| [07-modern-csharp.md](07-modern-csharp.md) | Modern C# | Mixed | Records, NRT, pattern matching, C# 11-12 |

### Advanced Examples

| File | Focus | Description |
|------|-------|-------------|
| [08-complete-example.md](08-complete-example.md) | **All Principles** | Complete Banking feature with all principles |
| [09-refactoring-journey.md](09-refactoring-journey.md) | **Refactoring** | Step-by-step code transformation (Reporting) |
| [10-architecture-integration.md](10-architecture-integration.md) | **Architecture** | Clean Code + Clean Architecture integration |
| [11-metrics-and-tools.md](11-metrics-and-tools.md) | **Metrics** | Code quality measurement and tooling |
| [12-decision-trees.md](12-decision-trees.md) | **Visual Guides** | Decision trees for clean code choices |

## Domain Diversity

Examples use **7 different business domains** to demonstrate universal applicability:

| Domain | Files | Key Entities |
|--------|-------|--------------|
| **Gaming** | 01-naming | Player, Game, Score, Level, Achievement |
| **IoT** | 02-functions | Device, Sensor, Reading, Alert, Telemetry |
| **Banking** | 03-error, 08-complete | Account, Money, Transfer, Transaction |
| **Healthcare** | 04-objects | Patient, Appointment, MedicalRecord, Doctor |
| **Education** | 05-refactoring | Student, Course, Grade, Enrollment |
| **Social Media** | 06-testing | User, Post, Comment, Like, Follow |
| **Reporting** | 09-journey | ReportGenerator, DataSource, Template |
| **Insurance** | 10-architecture | Policy, Claim, PolicyHolder |

## Requirements

- **.NET 8+** (or .NET 6+ for most examples)
- **C# 12** (C# 11 minimum for some features)
- **Nullable Reference Types enabled**

```xml
<!-- .csproj -->
<PropertyGroup>
    <TargetFramework>net8.0</TargetFramework>
    <Nullable>enable</Nullable>
    <ImplicitUsings>enable</ImplicitUsings>
</PropertyGroup>
```

## How to Use

### Learning Path
1. Start with [01-naming.md](01-naming.md) — foundation of readable code
2. Progress through [02-functions.md](02-functions.md) — building blocks
3. Master [03-error-handling.md](03-error-handling.md) — robust code
4. Understand [04-objects-vs-data.md](04-objects-vs-data.md) — design decisions
5. Apply [07-modern-csharp.md](07-modern-csharp.md) — leverage latest features

### During Code Review
Reference specific examples when discussing code quality issues.

### Refactoring
Use [05-refactoring.md](05-refactoring.md) for step-by-step improvement patterns.

## Example Format

Each file uses consistent format:

```markdown
### ❌ BAD
[Anti-pattern code]

### ✅ GOOD  
[Clean code solution]
```

## Relationship to Principles

These examples are **C#-specific implementations** of the language-agnostic principles in [principles/](../../principles/):

| Example | Principle File |
|---------|----------------|
| 01-naming.md | [01-meaningful-names.md](../../principles/01-meaningful-names.md) |
| 02-functions.md | [02-functions.md](../../principles/02-functions.md) |
| 03-error-handling.md | [06-error-handling.md](../../principles/06-error-handling.md) |
| 04-objects-vs-data.md | [05-objects-and-data-structures.md](../../principles/05-objects-and-data-structures.md) |
| 06-testing.md | [08-unit-tests.md](../../principles/08-unit-tests.md) |

## Related

- [Clean Code Principles](../../PRINCIPLES.md) — Theory and concepts
- [Clean Code Checklist](../../CHECKLIST.md) — Code review checklist
- [Clean Architecture Examples](../../../clean-architecture/examples/csharp/) — Architectural patterns
