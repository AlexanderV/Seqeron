---
name: clean-code
version: 2.0.0
description: Guide for writing clean, maintainable, and readable code following Clean Code principles. Use when reviewing code quality, refactoring, naming variables/functions, improving readability, reducing complexity, or when user mentions code smells, technical debt, maintainability, or asks for code quality improvements.
allowed-tools: Read, Grep, Glob, Edit, Write
triggers: [
  "clean code", "code quality", "readable code", "maintainable code",
  "refactor", "refactoring", "code smell", "code smells",
  "naming", "variable names", "function names", "meaningful names",
  "comments", "commenting", "documentation",
  "function size", "method length", "too long", "complex function",
  "cyclomatic complexity", "cognitive complexity",
  "DRY", "don't repeat yourself", "duplication", "duplicate code",
  "KISS", "keep it simple",
  "YAGNI", "you aren't gonna need it",
  "magic numbers", "magic strings", "constants",
  "nested if", "nested loops", "deeply nested",
  "single responsibility", "does too much",
  "code review", "review code quality",
  "technical debt", "legacy code", "unreadable code",
  "boy scout rule", "leave it cleaner",
  "async await", "concurrency", "thread safety"
]
---

# Clean Code Development Guide

This skill helps you write clean, maintainable, and readable code following Robert C. Martin's (Uncle Bob) Clean Code principles.

> "Clean code is simple and direct. Clean code reads like well-written prose." — Grady Booch

## Quick Links

| Document | Description |
|----------|-------------|
| [CHECKLIST.md](CHECKLIST.md) | Code review checklist |
| [PRINCIPLES.md](PRINCIPLES.md) | Index of all principles |
| [Examples](examples/csharp/) | C# code examples |
| [Metrics & Tools](examples/csharp/11-metrics-and-tools.md) | Code quality measurement |
| [Decision Trees](examples/csharp/12-decision-trees.md) | Visual decision guides |

## When to Use This Skill

- Reviewing code for quality issues
- Refactoring existing code for better readability
- Improving naming conventions
- Simplifying complex functions
- Reducing code duplication
- Eliminating code smells
- Writing self-documenting code
- Reducing technical debt

## Relationship with Clean Architecture

| Aspect | Clean Code | Clean Architecture |
|--------|------------|-------------------|
| Focus | Micro-level | Macro-level |
| Scope | Functions, classes, naming | Layers, dependencies, modules |
| Goal | Readable, maintainable code | Flexible, testable systems |

**Use both together:** Clean Architecture provides structure; Clean Code ensures each piece is well-written.

See: [Clean Architecture Skill](../clean-architecture/SKILL.md)

---

## Core Principles Overview

### 1. [Meaningful Names](principles/01-meaningful-names.md)

Names should reveal intent, avoid disinformation, and be pronounceable.

| Rule | Bad | Good |
|------|-----|------|
| Reveal intent | `int d` | `int elapsedDays` |
| Avoid abbreviations | `string custNm` | `string customerName` |
| Class names = nouns | `class Process` | `class OrderProcessor` |
| Method names = verbs | `void Order()` | `void PlaceOrder()` |
| Booleans = predicates | `bool flag` | `bool isActive` |

### 2. [Functions](principles/02-functions.md)

Functions should be small, do one thing, and have minimal parameters.

| Metric | Target |
|--------|--------|
| Lines | 5-20 (max 30) |
| Parameters | 0-2 (max 3) |
| Levels of abstraction | 1 per function |
| Nesting depth | Max 2 levels |

**Key rules:**
- Do ONE thing
- No side effects
- Command-Query Separation
- Prefer exceptions to error codes

### 3. [Comments](principles/03-comments.md)

Good code is self-documenting. Comments should explain **WHY**, not **WHAT**.

**✅ Acceptable comments:**
- Legal/copyright notices
- Intent explanation
- Warning of consequences
- TODO with ticket references
- Public API documentation

**❌ Delete these:**
- Redundant comments
- Commented-out code
- Journal/changelog comments
- Closing brace comments

### 4. [Formatting](principles/04-formatting.md)

Code formatting is about communication.

| Metric | Target | Max |
|--------|--------|-----|
| File size | ~200 lines | 500 |
| Line width | 80-100 chars | 120 |

**Rules:**
- Vertical: Related code close together
- Horizontal: Consistent indentation
- Team: Agree on standards and enforce with tools

### 5. [Objects vs Data Structures](principles/05-objects-and-data-structures.md)

| Aspect | Objects | Data Structures |
|--------|---------|-----------------|
| Data | Hidden | Exposed |
| Behavior | Rich methods | None |
| Use case | Domain entities | DTOs |

**Law of Demeter:** Don't talk to strangers. Avoid `a.GetB().GetC().DoSomething()`.

### 6. [Error Handling](principles/06-error-handling.md)

Errors should be handled gracefully without obscuring logic.

**Key rules:**
- Use exceptions, not error codes
- Provide context in exceptions
- Don't return null → use Null Object, Option, or throw
- Don't pass null
- Wrap third-party exceptions

### 7. [Boundaries](principles/07-boundaries.md)

Manage third-party code with wrappers and adapters.

- Wrap third-party APIs in adapters
- Use interfaces you control
- Write learning tests for third-party code
- Convert to domain types at boundaries

### 8. [Unit Tests](principles/08-unit-tests.md)

Test code is just as important as production code.

**F.I.R.S.T.:**
- **F**ast
- **I**ndependent
- **R**epeatable
- **S**elf-validating
- **T**imely

**Key rules:**
- One concept per test
- Arrange-Act-Assert pattern
- Test edge cases
- Don't test private methods

### 9. [Classes](principles/09-classes.md)

Classes should be small and have single responsibility.

| Metric | Target |
|--------|--------|
| Responsibilities | 1 |
| Public methods | 5-10 max |
| Lines | 100-200 typical |
| Dependencies | 3-5 max |

**Key rules:**
- Single Responsibility Principle
- High cohesion
- Depend on abstractions
- Avoid "Manager", "Handler", "Processor" names

### 10. [Emergence](principles/10-emergence.md)

Kent Beck's four rules of simple design (in priority order):

1. **Runs all tests** — System must be verifiable
2. **No duplication** — DRY principle
3. **Expresses intent** — Self-documenting code
4. **Minimal elements** — No unnecessary complexity

### 11. [Concurrency](principles/11-concurrency.md)

Writing clean concurrent code requires special care.

**Key rules:**
- Separate concurrency from business logic
- Minimize shared mutable state
- Use immutable objects
- Async/await all the way
- Support cancellation
- Avoid deadlocks with consistent lock ordering

### 12. [Smells & Heuristics](principles/12-smells-and-heuristics.md)

Comprehensive catalog of code smells from Clean Code Appendix B.

Categories: Comments (C), Environment (E), Functions (F), General (G), Names (N), Tests (T)

---

## Quick Decision Tree

```
Is this code hard to understand?
├─ Bad naming? → Rename to reveal intent
├─ Too long? → Extract smaller pieces
├─ Too complex? → Simplify logic, reduce nesting
├─ Duplicated? → Extract to shared method
└─ Needs comments? → Refactor until self-documenting

Is this function doing too much?
├─ Multiple abstraction levels? → Extract methods
├─ Multiple responsibilities? → Split into focused functions
├─ Too many parameters? → Create parameter object
└─ Side effects? → Separate commands from queries

Is this class too large?
├─ Multiple reasons to change? → Split by responsibility (SRP)
├─ Too many dependencies? → Review architecture
└─ Hard to test? → Apply Dependency Inversion
```

---

## Code Quality Metrics

| Metric | Target | Tool |
|--------|--------|------|
| Cyclomatic Complexity | < 10/method | Roslyn Analyzers |
| Code Coverage | > 80% business logic | dotnet test |
| Maintainability Index | > 70 | SonarQube |
| Line Length | < 120 chars | EditorConfig |

---

## The Boy Scout Rule

> "Leave the code cleaner than you found it."

When touching code, always:
- Fix bad names
- Extract long methods
- Remove duplication
- Delete commented code
- Simplify complex expressions

---

## Language-Specific Examples

### C# / .NET
- [Naming Conventions](examples/csharp/01-naming.md)
- [Functions](examples/csharp/02-functions.md)
- [Error Handling](examples/csharp/03-error-handling.md)
- [Objects vs Data](examples/csharp/04-objects-vs-data.md)
- [Refactoring](examples/csharp/05-refactoring.md)
- [Unit Testing](examples/csharp/06-testing.md)
- [Modern C# Practices](examples/csharp/07-modern-csharp.md)
- [Complete Feature Example](examples/csharp/08-complete-example.md) ⭐
- [Refactoring Journey](examples/csharp/09-refactoring-journey.md) — Step-by-step transformation

---

## Working with This Skill

When I use this skill, I will:

1. **Analyze** code for quality issues
2. **Identify** code smells and anti-patterns
3. **Suggest** concrete refactorings
4. **Apply** clean code principles
5. **Ensure** code remains testable and maintainable
6. **Reference** Clean Architecture when structural issues are found

---

## Reference Materials

- [PRINCIPLES.md](PRINCIPLES.md) — Index of all detailed principles
- [CHECKLIST.md](CHECKLIST.md) — Code review checklist
- [Clean Architecture](../clean-architecture/SKILL.md) — Architectural patterns

**Books:**
- *Clean Code* by Robert C. Martin
- *Refactoring* by Martin Fowler
- *Code Complete* by Steve McConnell
