# The Dependency Rule

> **Navigation:** [← Back to Principles](../PRINCIPLES.md) | [Code Examples →](../examples/csharp/PRINCIPLES-EXAMPLES.md)

## Definition
Source code dependencies must point only inward, toward higher-level policies.

## Explanation
- **Inner circles** = high-level policy (business rules, entities)
- **Outer circles** = low-level details (UI, database, frameworks)
- Nothing in an inner circle can know anything about something in an outer circle
- Data formats from outer circle must not be used by inner circle

## Why It Matters
- Business logic becomes independent of frameworks
- Database can be swapped without changing business rules
- UI can change without touching core logic
- External services can be replaced easily

## Key Insight
The dependency arrow points FROM Infrastructure TO Domain (through the interface), not the other way around.

## Related Principles

- [SOLID: Dependency Inversion Principle](02-solid-principles.md#d---dependency-inversion-principle-dip) - Concrete implementation of this rule
- [Hexagonal Architecture](04-hexagonal-architecture.md) - Applying Dependency Rule through ports and adapters
- [Onion Architecture](05-onion-architecture.md) - Another way to visualize this rule
