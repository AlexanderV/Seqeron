# Clean Architecture Principles - Detailed Reference

This document provides comprehensive explanations of all architectural principles used in Clean Architecture.

> **For Code Examples:** See [SOLID & DDD Examples](examples/csharp/PRINCIPLES-EXAMPLES.md) for practical implementations.

## Clean Architecture Principles

### 1. [The Dependency Rule](principles/01-dependency-rule.md)
Source code dependencies must point only inward, toward higher-level policies.

### 2. [SOLID Principles](principles/02-solid-principles.md)
Five fundamental principles of object-oriented design:
- **S** - Single Responsibility Principle
- **O** - Open/Closed Principle
- **L** - Liskov Substitution Principle
- **I** - Interface Segregation Principle
- **D** - Dependency Inversion Principle

### 3. [Domain-Driven Design (DDD)](principles/03-domain-driven-design.md)
Tactical patterns for modeling domain logic: Entities, Value Objects, Aggregates, Domain Services, Repositories, Factories, Specifications.

### 4. [Hexagonal Architecture (Ports & Adapters)](principles/04-hexagonal-architecture.md)
Isolating business logic from the external world through ports (interfaces) and adapters (implementations).

### 5. [Onion Architecture](principles/05-onion-architecture.md)
Abstraction layers with dependencies pointing inward: Domain Model → Domain Services → Application Services → Infrastructure.

### 6. [CQRS Pattern](principles/06-cqrs-pattern.md)
Command Query Responsibility Segregation: Commands change state, Queries return data.

### 7. [Screaming Architecture](principles/07-screaming-architecture.md)
Architecture should "scream" about the application's purpose, not the frameworks used.

### 8. [Package Principles](principles/08-package-principles.md)
Principles for organizing packages and modules: cohesion, coupling, stability, and abstraction.

### 9. [Putting It All Together](principles/09-putting-it-together.md)
Integration of all principles: decision matrix and application execution flow.

### 10. [Vertical Slices (Feature-Based Structure)](principles/10-vertical-slices.md)
Organize code by business features rather than technical layers. Each slice contains all layers for one use case, enabling high cohesion, team scalability, and natural CQRS alignment.

---

## Summary

**The Goal:** Software that is:
- **Testable** - Business logic can be tested without database/UI
- **Independent of Frameworks** - Frameworks are tools, not architecture
- **Independent of UI** - UI can change without changing business rules
- **Independent of Database** - Can swap Oracle for MongoDB
- **Independent of External Services** - Business rules don't know about external APIs

**The Method:**
1. **Apply Dependency Rule** - Dependencies point inward
2. **Follow SOLID** - Ensure clean, maintainable code
3. **Use DDD Tactically** - Entities, Value Objects, Aggregates model the domain
4. **Define Clear Boundaries** - Ports & Adapters, Bounded Contexts
5. **Separate Concerns** - CQRS, layers, modules
6. **Let Architecture Scream** - Structure reflects intent

**For working code examples:**
- [SOLID & DDD Examples](examples/csharp/PRINCIPLES-EXAMPLES.md) - Practical implementations
- [Design Pattern Examples](examples/csharp/PATTERNS-EXAMPLES.md) - Pattern catalog with code
