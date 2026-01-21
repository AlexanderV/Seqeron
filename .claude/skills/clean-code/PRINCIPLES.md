# Clean Code Principles

This document serves as an index to the detailed Clean Code principles. Each principle has its own dedicated file with comprehensive explanations and examples.

## Principles Index

| # | Principle | Description |
|---|-----------|-------------|
| 1 | [Meaningful Names](principles/01-meaningful-names.md) | Names that reveal intent and avoid disinformation |
| 2 | [Functions](principles/02-functions.md) | Small functions that do one thing well |
| 3 | [Comments](principles/03-comments.md) | When to comment and when to refactor instead |
| 4 | [Formatting](principles/04-formatting.md) | Vertical and horizontal code organization |
| 5 | [Objects and Data Structures](principles/05-objects-and-data-structures.md) | The fundamental dichotomy in OOP |
| 6 | [Error Handling](principles/06-error-handling.md) | Exceptions, null handling, and error context |
| 7 | [Boundaries](principles/07-boundaries.md) | Managing third-party code and adapters |
| 8 | [Unit Tests](principles/08-unit-tests.md) | F.I.R.S.T. principles and clean test design |
| 9 | [Classes](principles/09-classes.md) | Small classes with single responsibility |
| 10 | [Emergence](principles/10-emergence.md) | Kent Beck's four rules of simple design |
| 11 | [Concurrency](principles/11-concurrency.md) | Thread safety, async/await, immutability |
| 12 | [Smells & Heuristics](principles/12-smells-and-heuristics.md) | Comprehensive code smell catalog |

## Quick Summary

### The Foundation

1. **Meaningful Names** — Code should read like well-written prose
2. **Functions** — Small, focused, one level of abstraction
3. **Comments** — Code should be self-documenting

### Organization

4. **Formatting** — Communication through structure
5. **Objects vs Data** — Know when to use each
6. **Classes** — Single responsibility, high cohesion

### Quality

7. **Error Handling** — Handle errors without obscuring logic
8. **Boundaries** — Isolate third-party code
9. **Unit Tests** — Tests are first-class citizens

### Advanced

10. **Emergence** — Simple design rules
11. **Concurrency** — Clean concurrent code
12. **Smells** — Recognize and fix code problems

## The Four Rules of Simple Design

From [Emergence](principles/10-emergence.md) (in priority order):

1. **Runs all tests** — System must be verifiable
2. **No duplication** — DRY principle
3. **Expresses intent** — Self-documenting code
4. **Minimal elements** — No unnecessary complexity

## The Boy Scout Rule

> "Leave the code cleaner than you found it."

## See Also

- [SKILL.md](SKILL.md) — Main skill guide with quick reference
- [CHECKLIST.md](CHECKLIST.md) — Code review checklist
- [Examples](examples/csharp/) — Language-specific code examples
- [Clean Architecture](../clean-architecture/SKILL.md) — Architectural patterns

## Books

- *Clean Code* by Robert C. Martin
- *Refactoring* by Martin Fowler
- *Code Complete* by Steve McConnell
