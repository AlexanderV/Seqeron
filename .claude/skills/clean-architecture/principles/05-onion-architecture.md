# Onion Architecture

> **Navigation:** [← Back to Principles](../PRINCIPLES.md) | [Code Examples →](../examples/csharp/SKILL-EXAMPLES.md)

## Core Concept

Similar to Clean Architecture - layers of abstraction with dependencies pointing inward.

## Layers (Inside to Outside)

### 1. Domain Model (Center)
- Entities
- Value Objects
- Domain Events
- NO dependencies

### 2. Domain Services
- Domain logic spanning multiple entities
- Depends only on Domain Model

### 3. Application Services
- Use Cases
- Orchestrates domain objects
- Defines repository interfaces
- Depends on Domain Model and Domain Services

### 4. Infrastructure (Outer)
- Implements repository interfaces
- Database, file system, external APIs
- Depends on Application Services (through interfaces)

**Key Difference from Clean Architecture:** Onion explicitly separates Domain Services from Application Services.

## Related Principles

- [The Dependency Rule](01-dependency-rule.md) - Foundation of Onion Architecture
- [Hexagonal Architecture](04-hexagonal-architecture.md) - Alternative visualization of the same concept
- [Domain-Driven Design](03-domain-driven-design.md) - Explicit separation of Domain Services and Application Services
