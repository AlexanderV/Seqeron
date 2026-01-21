# CQRS Pattern

> **Navigation:** [← Back to Principles](../PRINCIPLES.md) | [Code Examples →](../examples/csharp/PATTERNS-EXAMPLES.md)

## Definition

Command Query Responsibility Segregation - separate read and write operations.

## Core Principle

- **Commands:** Change state, return void (or status)
- **Queries:** Return data, never change state

## Benefits

- Optimized read/write models
- Scalability (scale reads and writes independently)
- Clarity (intent is explicit)

## Advanced Concept

Write Model can be normalized (3NF), Read Model can be denormalized for performance.

## Related Principles

- [Interface Segregation Principle](02-solid-principles.md#i---interface-segregation-principle-isp) - CQRS is an application of ISP to operations
- [Domain-Driven Design](03-domain-driven-design.md) - Commands and Queries work with domain objects
