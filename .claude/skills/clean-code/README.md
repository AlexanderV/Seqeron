# Clean Code Skill

A comprehensive guide for writing clean, maintainable, and readable code following Robert C. Martin's Clean Code principles.

## Overview

This skill helps developers write code that is:
- **Readable** - Easy to understand at a glance
- **Maintainable** - Easy to change without breaking
- **Testable** - Easy to verify correctness
- **Simple** - No unnecessary complexity

## Quick Start

1. **Read** [SKILL.md](SKILL.md) for an overview and quick reference
2. **Review** [CHECKLIST.md](CHECKLIST.md) when doing code reviews
3. **Study** [PRINCIPLES.md](PRINCIPLES.md) for detailed explanations
4. **Explore** language-specific examples in [examples/](examples/)

## Documentation Structure

```
clean-code/
├── SKILL.md                   # Main skill file with overview and triggers
├── PRINCIPLES.md              # Detailed principles explanation
├── CHECKLIST.md               # Code review checklist
├── README.md                  # This file
└── examples/
    └── csharp/
        ├── 01-naming.md       # Naming conventions
        ├── 02-functions.md    # Function best practices
        ├── 03-error-handling.md
        ├── 04-objects-vs-data.md
        ├── 05-refactoring.md
        └── 06-testing.md
```

## Core Principles

1. **Meaningful Names** - Names reveal intent
2. **Small Functions** - Do one thing well (5-20 lines)
3. **Comments Are a Failure** - Code should be self-documenting
4. **DRY** - Don't Repeat Yourself
5. **Error Handling** - Use exceptions, not error codes
6. **Formatting** - Consistent and clean
7. **Objects vs Data** - Clear separation

## Relationship with Clean Architecture

**Clean Code** and [Clean Architecture](../clean-architecture/README.md) work together:

- **Clean Code**: Micro-level (functions, classes, naming)
- **Clean Architecture**: Macro-level (layers, dependencies, modules)

Use both together for maintainable software.

## When to Use This Skill

- Reviewing code for quality
- Refactoring legacy code
- Improving naming conventions
- Simplifying complex functions
- Eliminating code duplication
- Writing self-documenting code

## Code Review Process

1. Run automated tools (linters, analyzers)
2. Use [CHECKLIST.md](CHECKLIST.md) systematically
3. Focus on high-priority issues first
4. Apply **Boy Scout Rule**: Leave code cleaner than you found it

## Language Support

### Currently Available
- **C#** - Complete examples in [examples/csharp/](examples/csharp/)

### Planned
- JavaScript/TypeScript
- Python
- Java
- Go

## Tools and Resources

### C# / .NET
- **Roslyn Analyzers** - Code analysis
- **StyleCop** - Style enforcement
- **SonarQube** - Quality metrics
- **ReSharper** - Refactoring tools

### Books
- **Clean Code** by Robert C. Martin
- **Refactoring** by Martin Fowler
- **Code Complete** by Steve McConnell

## Metrics

Track these metrics to measure improvement:

| Metric | Target |
|--------|--------|
| Cyclomatic Complexity | < 10 per method |
| Lines per Method | < 20 |
| Lines per Class | < 300 |
| Code Coverage | > 80% for business logic |
| Maintainability Index | > 70 |

## Quick Decision Tree

```
Code hard to understand?
├─ Bad names? → Rename
├─ Too long? → Extract methods
├─ Too complex? → Simplify
├─ Duplicated? → DRY
└─ Needs comments? → Refactor until self-documenting

Function doing too much?
├─ Multiple abstractions? → Extract per level
├─ Multiple responsibilities? → Split
├─ Too many parameters? → Parameter object
└─ Side effects? → Separate commands/queries

Class too large?
├─ Multiple responsibilities? → Split (SRP)
├─ Too many dependencies? → See Clean Architecture
└─ Hard to test? → Dependency Inversion
```

## Examples

### Before (Bad)
```csharp
int d; // elapsed time in days
public void proc(Order o) {
    if (o.s == 2) {
        // 50 lines of mixed concerns...
    }
}
```

### After (Good)
```csharp
int elapsedTimeInDays;
public void ProcessApprovedOrder(Order order) {
    ValidateOrder(order);
    CalculateTotal(order);
    SaveOrder(order);
}
```

## Contributing

To add examples for a new language:

1. Create directory: `examples/{language}/`
2. Follow the same structure as C# examples
3. Adapt examples to language idioms
4. Update this README

## Getting Help

- Review [SKILL.md](SKILL.md) for triggers and overview
- Check [PRINCIPLES.md](PRINCIPLES.md) for detailed explanations
- Use [CHECKLIST.md](CHECKLIST.md) during code reviews
- Explore [examples/csharp/](examples/csharp/) for patterns

## Related Skills

- [Clean Architecture](../clean-architecture/README.md) - Architectural patterns and layer design

---

> "Any fool can write code that a computer can understand. Good programmers write code that humans can understand."
>
> — Martin Fowler

> "Clean code is simple and direct. Clean code reads like well-written prose."
>
> — Grady Booch

**Remember:** Leave the code cleaner than you found it!
