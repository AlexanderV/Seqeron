# Putting It All Together

> **Navigation:** [← Back to Principles](../PRINCIPLES.md)

## Decision Matrix

| Concern | Principle/Pattern | Key Question |
|---------|------------------|--------------|
| Layer dependency | Dependency Rule | Does this point inward? |
| Class responsibility | SRP | How many reasons to change? |
| Adding features | OCP | Can I extend without modifying? |
| Substitutability | LSP | Can I swap implementations? |
| Interface design | ISP | Am I forcing unused methods? |
| Dependencies | DIP | Am I depending on abstractions? |
| Business concepts | DDD | Does this reflect domain language? |
| External systems | Hexagonal | Did I define a port? |
| Read vs Write | CQRS | Is this changing or querying? |
| Project structure | Screaming | Does structure show intent? |

## Application Flow

1. User Request
2. Driving Adapter (Controller)
3. Driving Port (IOrderService interface)
4. Application Service (Use Case)
5. Domain Model (Order.placeOrder())
6. Driven Port (IOrderRepository interface)
7. Driven Adapter (SqlOrderRepository)
8. Database

All dependencies point toward Domain Model (steps 1-5).
Infrastructure (steps 6-8) depends on abstractions defined in steps 4-5.

## Related Principles

- [The Dependency Rule](01-dependency-rule.md) - Foundational rule of the entire architecture
- [SOLID Principles](02-solid-principles.md) - Applied at all levels
- [Domain-Driven Design](03-domain-driven-design.md) - Tactical patterns for the domain layer
- [Hexagonal Architecture](04-hexagonal-architecture.md) - Ports and adapters for isolation
- [CQRS Pattern](06-cqrs-pattern.md) - Separation of commands and queries
