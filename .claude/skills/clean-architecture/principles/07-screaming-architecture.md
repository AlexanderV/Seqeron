# Screaming Architecture

> **Navigation:** [← Back to Principles](../PRINCIPLES.md) | [Code Examples →](../examples/csharp/SKILL-EXAMPLES.md)

## Definition

The architecture should "scream" the intent of the application, not the frameworks used.

## Principle

Looking at the project structure should immediately tell you what the application does, not what frameworks it uses.

## Key Insight

The folder structure reflects use cases and domain concepts, not technical patterns (Controllers, Models, Views).

**Bad:** MyApp/Controllers/Models/Views - tells you it's MVC, not what it does.
**Good:** MyApp/Orders/Inventory/Shipping/Payments - tells you it's e-commerce.

## Related Principles

- [Domain-Driven Design](03-domain-driven-design.md) - Ubiquitous Language applies to project structure
- [The Dependency Rule](01-dependency-rule.md) - Structure reflects architecture layers
