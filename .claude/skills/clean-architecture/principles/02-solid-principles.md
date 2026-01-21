# SOLID Principles

> **Navigation:** [← Back to Principles](../PRINCIPLES.md) | [Code Examples →](../examples/csharp/PRINCIPLES-EXAMPLES.md)

## S - Single Responsibility Principle (SRP)

**Definition:** A class should have one, and only one, reason to change.

**Explanation:**
- A "reason to change" represents an actor or stakeholder
- If a class serves multiple actors, it violates SRP
- Separating responsibilities prevents unintended side effects

**Key Question:** Can you describe the class without using "and" or "or"?

## O - Open/Closed Principle (OCP)

**Definition:** Software entities should be open for extension, but closed for modification.

**Explanation:**
- Add new functionality by adding new code, not changing existing code
- Use abstraction and polymorphism
- Prevents bugs in tested code

**Implementation:** Use interfaces and Strategy pattern to allow extension without modification.

## L - Liskov Substitution Principle (LSP)

**Definition:** Subtypes must be substitutable for their base types without altering program correctness.

**Explanation:**
- If S is a subtype of T, then objects of type T can be replaced with objects of type S
- Subclasses must honor the contract of the base class
- Violation leads to unexpected behavior

**Key Rule:** Derived classes cannot strengthen preconditions or weaken postconditions.

## I - Interface Segregation Principle (ISP)

**Definition:** No client should be forced to depend on methods it does not use.

**Explanation:**
- Many small, specific interfaces are better than one large, general interface
- Prevents "fat" interfaces
- Clients only depend on what they actually need

**Guideline:** If implementing an interface requires throwing NotImplementedException, the interface is too broad.

## D - Dependency Inversion Principle (DIP)

**Definition:**
- High-level modules should not depend on low-level modules. Both should depend on abstractions.
- Abstractions should not depend on details. Details should depend on abstractions.

**Explanation:**
- "High-level" = business logic, use cases
- "Low-level" = technical details, infrastructure
- Both depend on interfaces (abstractions)

**Key Insight:** The dependency arrow points FROM the concrete class TO the interface, inverting the traditional dependency direction.

## Related Principles

- [The Dependency Rule](01-dependency-rule.md) - DIP is a concrete implementation of this rule
- [Domain-Driven Design](03-domain-driven-design.md) - SOLID principles apply when creating domain objects
- [Hexagonal Architecture](04-hexagonal-architecture.md) - Uses DIP to create ports and adapters
