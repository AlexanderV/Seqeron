# Design Patterns for Clean Architecture

This document catalogs essential design patterns used in Clean Architecture applications.

> **For Code Examples:** See [Design Pattern Examples](examples/csharp/PATTERNS-EXAMPLES.md) for complete implementations with MediatR, EF Core, and other frameworks.

## Table of Contents

1. [Creational Patterns](#creational-patterns)
2. [Structural Patterns](#structural-patterns)
3. [Behavioral Patterns](#behavioral-patterns)
4. [Architectural Patterns](#architectural-patterns)

---

## Creational Patterns

### Factory Pattern

**Purpose:** Encapsulate complex object creation logic.

**When to Use:**
- Creating aggregates with multiple invariants to validate
- Object construction requires business rules
- Hide complex initialization from clients

**Key Benefits:**
- Centralizes complex creation logic
- Enforces invariants during construction
- Provides named constructors for clarity (e.g., createGiftOrder, createBulkOrder)

---

### Builder Pattern

**Purpose:** Construct complex objects step-by-step.

**When to Use:**
- Object has many optional parameters
- Construction process has multiple steps
- Need to enforce construction order

**Advantages:**
- Fluent interface (method chaining)
- Immutable objects after build()
- Clear intention through named methods

---

### Dependency Injection Pattern

**Purpose:** Invert dependencies - depend on abstractions, not concretions.

**When to Use:**
- Always (in Clean Architecture)
- Decoupling layers
- Enabling testability

**Implementation:**
- Domain/Application defines interface
- Infrastructure provides implementation
- Composition Root wires dependencies

---

## Structural Patterns

### Repository Pattern

**Purpose:** Abstract data access and provide collection-like interface for aggregates.

**When to Use:**
- Always, for aggregate roots
- Hide persistence implementation details
- Provide domain-friendly query interface

**Key Principles:**
- Interface in Domain, implementation in Infrastructure
- Work with aggregate roots only
- Return domain entities, not database entities
- Hide all persistence details (SQL, ORM, etc.)

---

### Unit of Work Pattern

**Purpose:** Maintain a list of objects affected by a transaction and coordinate changes.

**When to Use:**
- Managing multiple repository operations as a single transaction
- Ensuring consistency across aggregate changes
- Optimizing database calls (batch updates)

**Benefits:**
- Transactional consistency across multiple repositories
- Single point to commit or rollback all changes
- Reduces database round-trips

---

### Adapter Pattern (Hexagonal Architecture)

**Purpose:** Convert one interface to another, allowing incompatible interfaces to work together.

**When to Use:**
- Integrating external services
- Translating between layers
- Implementing ports in Hexagonal Architecture

**Structure:**
- Application defines port (interface)
- Infrastructure provides adapter (implementation)
- Adapter translates between domain types and external service types

**Key Insight:** The application defines the interface. External services adapt to application needs.

---

### Mapper Pattern (DTO Mapping)

**Purpose:** Convert between domain models and DTOs/ViewModels.

**When to Use:**
- Presenting domain entities to external clients
- Receiving external data and converting to domain
- Preventing domain model leakage

**Direction:**
- toDto(): Domain → DTO (outbound)
- toDomain(): DTO → Domain (inbound)

**Benefit:** Keeps domain model independent of presentation concerns.

---

## Behavioral Patterns

### Strategy Pattern

**Purpose:** Encapsulate interchangeable algorithms and make them substitutable.

**When to Use:**
- Multiple ways to perform an operation
- Avoiding conditional logic for algorithm selection
- Runtime algorithm selection

**Application in Clean Architecture:**
- Discount calculation strategies
- Payment processing strategies
- Shipping cost calculation strategies

**Benefit:** Open/Closed Principle - add new strategies without modifying existing code.

---

### Command Pattern (CQRS)

**Purpose:** Encapsulate a request as an object, enabling queuing, logging, and undo operations.

**When to Use:**
- CQRS architecture
- Separating request from execution
- Queuing operations
- Implementing undo/redo

**Components:**
- **Command:** Immutable data object (CreateOrderCommand)
- **CommandHandler:** Executes the command logic
- **Result:** Return status or identifier

**Benefit:** Clear separation between intent (command) and execution (handler).

---

### Specification Pattern

**Purpose:** Encapsulate business rules as composable, reusable objects.

**When to Use:**
- Complex business rules
- Composing rules with AND/OR/NOT
- Reusing rules across queries and validation

**Operations:**
- isSatisfiedBy(entity): boolean
- and(other): Specification
- or(other): Specification
- not(): Specification

**Power:** Compose simple specifications into complex business rules.

---

### Observer Pattern (Domain Events)

**Purpose:** Define one-to-many dependency so that when one object changes state, all dependents are notified.

**When to Use:**
- Domain events
- Decoupling side effects from main flow
- Event-driven architecture

**Flow:**
1. Entity raises domain event
2. Application service collects events
3. Event dispatcher notifies all handlers
4. Handlers execute side effects (email, inventory reservation, etc.)

**Benefit:** Main flow remains clean, side effects decoupled.

---

### Chain of Responsibility Pattern

**Purpose:** Pass request along a chain of handlers until one handles it.

**When to Use:**
- Validation pipelines
- Middleware
- Processing requests through multiple steps

**Examples:**
- Validation → Authorization → Business Logic
- Request logging → Authentication → Rate limiting → Request processing

**Structure:** Each handler either processes the request or passes it to the next handler.

---

## Architectural Patterns

### CQRS (Command Query Responsibility Segregation)

**Purpose:** Separate read and write operations using different models.

**When to Use:**
- Complex domains with different read/write requirements
- High-performance read operations needed
- Event sourcing architecture

**Structure:**
- **Write Side:** Commands → Handlers → Domain Model → Write Database
- **Read Side:** Queries → Handlers → Read Model → Read Database
- **Synchronization:** Events update read model from write model

**Benefit:**
- Optimize each side independently
- Scale reads and writes separately
- Simple, focused models

---

### Event Sourcing

**Purpose:** Store state as a sequence of events instead of current state.

**When to Use:**
- Audit trail required
- Temporal queries (what was the state at time X?)
- Event-driven architecture
- Complex domains with many state changes

**Concept:**
- Events are immutable, append-only
- Current state reconstructed by replaying events
- Complete history preserved

**Benefits:**
- Perfect audit trail
- Time travel (reconstruct state at any point)
- Event-driven integration

**Trade-offs:**
- Increased complexity
- Eventually consistent reads
- Event schema evolution challenges

---

### Saga Pattern

**Purpose:** Manage distributed transactions across multiple aggregates/services.

**When to Use:**
- Distributed systems
- Long-running business processes
- Coordinating multiple aggregates
- Microservices architecture

**Types:**

**1. Choreography (Event-Driven):**
- Each service publishes events
- Other services subscribe and react
- Decentralized coordination

**2. Orchestration (Centralized):**
- Saga coordinator manages the process
- Explicit compensation logic
- Centralized control

**Compensation:**
- Each step has compensating transaction
- On failure, execute compensations in reverse order

---

## Pattern Selection Guide

| Concern | Pattern | Use When |
|---------|---------|----------|
| Creating complex objects | Factory | Validation/business rules during construction |
| Step-by-step construction | Builder | Many optional parameters |
| Decoupling layers | Dependency Injection | Always |
| Data access | Repository | Working with aggregates |
| Transactions | Unit of Work | Multiple repos in one transaction |
| External integration | Adapter | Integrating third-party services |
| DTO conversion | Mapper | Presenting domain to external |
| Algorithm selection | Strategy | Multiple ways to perform operation |
| Encapsulating requests | Command | CQRS, undo/redo, queuing |
| Business rules | Specification | Composable, reusable rules |
| Side effects | Observer | Domain events |
| Request processing | Chain of Responsibility | Validation pipelines, middleware |
| Read/Write separation | CQRS | Different read/write requirements |
| Audit trail | Event Sourcing | Need full history |
| Distributed transactions | Saga | Microservices, long-running processes |

---

## Summary

**Core Principles:**
1. **Patterns serve architecture** - Use patterns to implement Clean Architecture principles
2. **Don't over-pattern** - Only use patterns that solve actual problems
3. **Combine patterns** - Many patterns work together (Repository + Unit of Work + Mapper)
4. **Understand trade-offs** - Each pattern has costs and benefits

**For working implementations with code:**
- [Design Pattern Examples](examples/csharp/PATTERNS-EXAMPLES.md) - All patterns with complete C# code
- [SOLID & DDD Examples](examples/csharp/PRINCIPLES-EXAMPLES.md) - Foundation patterns
