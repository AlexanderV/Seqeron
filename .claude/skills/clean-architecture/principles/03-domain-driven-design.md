# Domain-Driven Design (DDD)

> **Navigation:** [← Back to Principles](../PRINCIPLES.md) | [Code Examples →](../examples/csharp/PRINCIPLES-EXAMPLES.md)

## Core Concepts

### 1. Ubiquitous Language
**Definition:** A common language shared by developers and domain experts.

**Principles:**
- Use the same terms in code that domain experts use
- Model should reflect business terminology
- Reduces translation overhead and misunderstandings

**Application:** Class names, method names, and variable names should match business vocabulary exactly.

### 2. Bounded Context
**Definition:** Explicit boundary within which a domain model is defined.

**Principles:**
- Different contexts can have different models for the same concept
- Clear boundaries prevent model pollution
- Each microservice typically represents one bounded context

**Key Rule:** A "Customer" in Sales context may have different properties than "Customer" in Shipping context.

### 3. Entities
**Definition:** Objects with distinct identity that persists over time.

**Characteristics:**
- Has unique identifier
- Equality based on identity, not attributes
- Mutable (state can change)
- Contains business logic

**Rule:** Identity remains constant even if all attributes change.

### 4. Value Objects
**Definition:** Objects defined by their attributes, not identity.

**Characteristics:**
- No unique identifier
- Immutable
- Equality based on all attributes
- Can be shared

**Common Examples:**
- Money (amount + currency)
- Email (validated string)
- Address (street, city, zip)
- DateRange (start + end)
- PhoneNumber (validated string)

### 5. Aggregates
**Definition:** Cluster of entities and value objects with a root entity.

**Principles:**
- Aggregate Root is the only entry point
- External objects can only reference the root
- Enforces invariants within the boundary
- Transaction boundary

**Rule:** All changes to aggregate members must go through the root.

### 6. Domain Services
**Definition:** Operations that don't naturally belong to an entity or value object.

**When to Use:**
- Operation involves multiple aggregates
- Operation is stateless
- Represents domain logic, not application logic

**Example:** Money transfer involves TWO accounts, so it belongs in a service, not an entity.

### 7. Domain Events
**Definition:** Something significant that happened in the domain.

**Characteristics:**
- Immutable
- Past tense naming (OrderPlaced, PaymentReceived)
- Contains all data needed by subscribers
- Triggers side effects

**Usage:** Decouple side effects from main business flow.

### 8. Repositories
**Definition:** Abstraction for accessing and persisting aggregates.

**Principles:**
- Interface defined in Domain layer
- Implementation in Infrastructure layer
- Works with aggregate roots only
- Hides persistence details

**Key Rule:** Repository provides collection-like interface - get, save, delete, find.

### 9. Factories
**Definition:** Encapsulates complex object creation logic.

**When to Use:**
- Creating aggregates with complex rules
- Object creation involves multiple steps
- Invariants must be validated during construction

**Benefit:** Centralizes creation logic and enforces invariants.

### 10. Specifications
**Definition:** Encapsulates business rules for querying and validation.

**Benefits:**
- Reusable business rules
- Composable with AND/OR/NOT
- Testable in isolation

**Pattern:** Specifications can be combined to create complex business rules.

## Related Principles

- [SOLID Principles](02-solid-principles.md) - SOLID applies when creating domain objects
- [The Dependency Rule](01-dependency-rule.md) - Domain resides in the inner circle
- [Hexagonal Architecture](04-hexagonal-architecture.md) - DDD tactical patterns are used in the Application Core
- [Onion Architecture](05-onion-architecture.md) - Explicitly separates Domain Model and Domain Services
