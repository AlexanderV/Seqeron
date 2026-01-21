# Clean Architecture Review Checklist

Use this checklist to review code for architectural compliance. Each section should be verified before considering the architecture sound.

> **For Code Examples:** This checklist contains only review criteria. For implementation examples, see:
> - [Complete Feature Example](examples/csharp/CHECKLIST-EXAMPLES.md) - Full vertical slice with tests
> - [Project Structure Examples](examples/csharp/SKILL-EXAMPLES.md) - Layer organization
> - [SOLID & DDD Examples](examples/csharp/PRINCIPLES-EXAMPLES.md) - Principles in practice
> - [Design Pattern Examples](examples/csharp/PATTERNS-EXAMPLES.md) - Pattern implementations

## Quick Reference

- ✅ Compliant
- ⚠️ Needs improvement
- ❌ Violation

---

## 1. Layer Dependencies

### 1.1 Dependency Rule Compliance

- [ ] Domain layer has NO dependencies on other layers
- [ ] Domain layer has NO framework dependencies (except language primitives)
- [ ] Application layer depends ONLY on Domain
- [ ] Infrastructure depends on Domain and Application (through interfaces)
- [ ] Presentation depends on Application (not directly on Infrastructure)
- [ ] No circular dependencies between layers
- [ ] Dependencies point inward (toward Domain)

### 1.2 Interface Segregation

- [ ] Infrastructure implements interfaces defined in Domain/Application
- [ ] No leaking of infrastructure details into Domain
- [ ] Interfaces are defined by their consumers (not implementers)

**Red Flags:**
- Domain code importing Infrastructure namespaces
- Application code directly instantiating database classes
- Concrete infrastructure types in Domain/Application signatures

---

## 2. Domain Layer

### 2.1 Entities

- [ ] Entities have unique identity
- [ ] Entities protect their invariants
- [ ] Entities expose behavior, not just properties
- [ ] Entity state changes go through methods (not public setters)
- [ ] Entities validate themselves
- [ ] No public parameterless constructors on entities
- [ ] No ORM/persistence attributes on domain entities

**Red Flags:**
- Public setters on entity properties
- Entities with only getters/setters (anemic model)
- Database annotations (@Entity, [Table], etc.) on domain entities
- Validation logic outside entities

### 2.2 Value Objects

- [ ] Value objects are immutable
- [ ] Value objects have no identity
- [ ] Equality based on all properties
- [ ] Value objects validate in constructor
- [ ] No public setters on value objects
- [ ] Value objects are side-effect free

**Red Flags:**
- Mutable properties on value objects
- Value objects with identity (ID field)
- Validation logic outside value object

### 2.3 Aggregates

- [ ] Clear aggregate boundaries defined
- [ ] Only aggregate root is accessible from outside
- [ ] Internal entities accessed only through root
- [ ] Aggregate enforces invariants across all contained entities
- [ ] Transaction boundary matches aggregate boundary
- [ ] Aggregate root has repository, internal entities don't

**Red Flags:**
- Direct access to internal entities
- Repositories for non-root entities
- Business rules spanning multiple aggregates in single transaction

### 2.4 Domain Services

- [ ] Domain services are stateless
- [ ] Domain services contain domain logic only (not application orchestration)
- [ ] Domain services used when logic doesn't belong to single entity
- [ ] Domain services operate on domain types
- [ ] No infrastructure dependencies in domain services

**Red Flags:**
- Domain service with state
- Domain service orchestrating multiple use cases
- Domain service depending on repositories

### 2.5 Domain Events

- [ ] Domain events are immutable
- [ ] Domain events use past tense naming (OrderPlaced, not PlaceOrder)
- [ ] Domain events contain all necessary data
- [ ] Entities raise events for significant state changes
- [ ] Events are value objects

**Red Flags:**
- Mutable event properties
- Present tense event names
- Events requiring additional lookups to be useful

---

## 3. Application Layer

### 3.1 Use Cases / Application Services

- [ ] One use case per class/function
- [ ] Use cases orchestrate domain objects
- [ ] Use cases are stateless
- [ ] Use cases don't contain business logic
- [ ] Use cases depend on interfaces (not concrete implementations)
- [ ] Use cases handle application concerns (transactions, events, etc.)

**Red Flags:**
- Multiple use cases in one class
- Business logic in use case
- Use case directly instantiating infrastructure

### 3.2 DTOs (Data Transfer Objects)

- [ ] DTOs are simple data structures
- [ ] DTOs have no business logic
- [ ] DTOs used for cross-boundary communication
- [ ] DTOs don't expose domain entities
- [ ] Separate DTOs for commands and queries

**Red Flags:**
- DTOs with methods beyond simple getters
- Domain entities used as DTOs
- Single DTO for multiple use cases

### 3.3 Interfaces (Ports)

- [ ] Interfaces defined in Application/Domain, not Infrastructure
- [ ] Interfaces designed for use case needs (not database structure)
- [ ] Interfaces return domain types (not infrastructure types)
- [ ] Interface segregation principle followed (many small > one large)

**Red Flags:**
- Interfaces in Infrastructure layer
- "Fat" interfaces with many methods
- Interfaces exposing ORM-specific methods

---

## 4. Infrastructure Layer

### 4.1 Repositories

- [ ] Repository implements interface from Domain
- [ ] Repository works with aggregate roots only
- [ ] Repository returns domain entities (not DB entities)
- [ ] Repository hides all persistence details
- [ ] Repository maps between domain and persistence models
- [ ] Repository methods are domain-centric (not SQL-centric)

**Red Flags:**
- Repository defined in Infrastructure
- Repository returning database entities
- Domain code aware of SQL/ORM
- Repository for non-aggregate-root

### 4.2 Adapters

- [ ] Adapters implement interfaces from Application/Domain
- [ ] Adapters translate between domain and external world
- [ ] Adapters handle external service failures gracefully
- [ ] Adapters map external data to domain types

**Red Flags:**
- External service types leaking to Domain/Application
- Domain code aware of third-party APIs
- No error handling for external services

### 4.3 Persistence

- [ ] Separate persistence models from domain models
- [ ] Mapping logic contained in Infrastructure
- [ ] Database schema doesn't dictate domain model
- [ ] No lazy loading in domain entities
- [ ] Transactions managed at application boundary (not in domain)

**Red Flags:**
- ORM entities = domain entities
- Domain model constrained by database
- Lazy-loaded properties in domain
- Transaction management in domain

---

## 5. Presentation Layer

### 5.1 Controllers / Handlers

- [ ] Controllers are thin (delegate to Application layer)
- [ ] Controllers don't contain business logic
- [ ] Controllers map between HTTP/UI and Application DTOs
- [ ] Controllers don't reference Infrastructure directly
- [ ] Controllers handle only presentation concerns

**Red Flags:**
- Business logic in controllers
- Controllers instantiating repositories
- Controllers with complex logic

### 5.2 ViewModels / Response DTOs

- [ ] ViewModels are UI-specific
- [ ] ViewModels don't leak domain entities
- [ ] Separate mapping from domain DTOs to ViewModels

**Red Flags:**
- Domain entities exposed in API
- Single DTO for both application and presentation

---

## 6. SOLID Principles

### 6.1 Single Responsibility Principle (SRP)

- [ ] Each class has one reason to change
- [ ] Classes serve single actor/stakeholder
- [ ] No "god" classes with multiple responsibilities

**Check Each Class:**
- Can you describe it without using "and" or "or"?
- Does it have exactly one reason to change?

### 6.2 Open/Closed Principle (OCP)

- [ ] Behavior extended through new classes (not modifying existing)
- [ ] Use interfaces/abstractions for extension points
- [ ] Avoid switch/if-else for type checking

**Red Flags:**
- Adding new feature requires modifying existing class
- Switch statements on types
- No abstraction for varying behavior

### 6.3 Liskov Substitution Principle (LSP)

- [ ] Subtypes substitutable for base types
- [ ] Derived classes honor base class contract
- [ ] No strengthening preconditions or weakening postconditions

**Red Flags:**
- Derived class throws exception not in base
- Derived class has stricter validation
- Derived class breaks base class assumptions

### 6.4 Interface Segregation Principle (ISP)

- [ ] Clients not forced to depend on unused methods
- [ ] Many small, specific interfaces
- [ ] No "fat" interfaces

**Red Flags:**
- Interface with many methods
- Implementations throwing NotImplementedException
- Clients depending on more than they use

### 6.5 Dependency Inversion Principle (DIP)

- [ ] High-level modules don't depend on low-level modules
- [ ] Both depend on abstractions
- [ ] Abstractions don't depend on details
- [ ] Dependency injection used throughout

**Red Flags:**
- Direct instantiation of concrete classes
- "new" keyword for infrastructure in domain/application
- No interfaces for dependencies

---

## 7. Domain-Driven Design (DDD)

### 7.1 Ubiquitous Language

- [ ] Code uses same terms as domain experts
- [ ] No translation between business terms and code
- [ ] Class/method names reflect business concepts

**Check:**
- Do variable/class names match business vocabulary?
- Would a domain expert understand the code?

### 7.2 Bounded Contexts

- [ ] Clear context boundaries identified
- [ ] No model pollution across contexts
- [ ] Each context has own model for shared concepts
- [ ] Anti-corruption layer for external contexts

**Red Flags:**
- Same entity class shared across contexts
- Context boundaries unclear
- Contexts directly calling each other

### 7.3 Strategic Design

- [ ] Context map documents relationships
- [ ] Upstream/downstream relationships identified
- [ ] Integration patterns chosen (Shared Kernel, ACL, etc.)

---

## 8. Testing

### 8.1 Domain Testing

- [ ] Domain logic has high unit test coverage
- [ ] Domain tests require no mocking
- [ ] Domain tests are fast
- [ ] Domain tests test business rules

**Red Flags:**
- Domain tests requiring database
- Domain tests with many mocks
- No tests for business rules

### 8.2 Application Testing

- [ ] Use cases tested with mocked dependencies
- [ ] Application tests verify orchestration
- [ ] Application tests don't test business logic (that's in domain)

### 8.3 Integration Testing

- [ ] Infrastructure tested against real dependencies
- [ ] Repository tests use real database (test db)
- [ ] Integration tests verify mappings
- [ ] Integration tests slow, kept separate from unit tests

### 8.4 Test Independence

- [ ] Domain can be tested without infrastructure
- [ ] Business logic tests don't require database/network
- [ ] Tests don't depend on execution order

**Red Flags:**
- Can't test without full environment
- Business logic tests needing database
- Tests failing when run in different order

---

## 9. Common Anti-Patterns

### 9.1 Anemic Domain Model

- [ ] **NOT PRESENT:** Entities are just data containers
- [ ] **NOT PRESENT:** Business logic in services instead of entities
- [ ] **PRESENT:** Behavior and data together in entities

### 9.2 God Objects

- [ ] **NOT PRESENT:** Classes with too many responsibilities
- [ ] **NOT PRESENT:** Single service doing everything
- [ ] **PRESENT:** Single Responsibility followed

### 9.3 Leaky Abstractions

- [ ] **NOT PRESENT:** Infrastructure details in Domain/Application
- [ ] **NOT PRESENT:** Domain entities exposed in API
- [ ] **PRESENT:** Clean separation through interfaces

### 9.4 Tight Coupling

- [ ] **NOT PRESENT:** Direct dependencies between layers
- [ ] **NOT PRESENT:** Classes instantiating their dependencies
- [ ] **PRESENT:** Dependency injection throughout

### 9.5 Transaction Script

- [ ] **NOT PRESENT:** Procedural code with no objects
- [ ] **NOT PRESENT:** Long methods with sequential steps
- [ ] **PRESENT:** Object-oriented domain model

---

## 10. Architecture Quality Metrics

### 10.1 Dependencies

- [ ] Dependency graph flows inward
- [ ] No circular dependencies
- [ ] Stable dependencies principle followed
- [ ] Stable abstractions principle followed

**Tools:**
- Dependency analyzers (NDepend, SonarQube, etc.)
- Architecture testing frameworks

### 10.2 Cohesion

- [ ] High cohesion within modules
- [ ] Related code grouped together
- [ ] Clear module boundaries

### 10.3 Coupling

- [ ] Low coupling between modules
- [ ] Coupling through interfaces only
- [ ] No implementation coupling

---

## 11. Vertical Slices (if applicable)

Use this section when Application Layer is organized by features.

### 11.1 Layer Separation

- [ ] Vertical Slices are in **Application Layer only**
- [ ] Domain Layer remains horizontal (shared across features)
- [ ] Presentation Layer is separate from Application
- [ ] Infrastructure Layer is horizontal (implements ports)

**Red Flags:**
- Endpoints in same folder as Handlers
- Domain entities inside slice folders
- Business logic in Handlers (should be in Domain)

### 11.2 Slice Structure

- [ ] Each slice represents one use case
- [ ] Slice contains: Command/Query, Handler, Validator, Response
- [ ] Endpoints are in Presentation layer, separate folder
- [ ] API DTOs are separate from Application Commands

**Correct Structure:**
```
Application/Features/Orders/CreateOrder/
├── CreateOrderCommand.cs      ← Application
├── CreateOrderHandler.cs      ← Application
├── CreateOrderValidator.cs    ← Application
└── CreateOrderResponse.cs     ← Application

Presentation/Orders/
├── CreateOrderEndpoint.cs     ← Presentation
└── Contracts/
    └── CreateOrderRequest.cs  ← API DTO
```

### 11.3 Dependency Rule in Slices

- [ ] Handlers depend on Domain Layer (not vice versa)
- [ ] Handlers use Ports (interfaces) defined in Application
- [ ] Handlers don't contain business logic (delegate to Domain)
- [ ] Domain Events dispatched AFTER persistence

**Red Flags:**
- Handler implementing business invariants
- Domain entities depending on Handler types
- Direct database calls in Handlers (should use repository)

### 11.4 Cross-Cutting Concerns

- [ ] Pipeline Behaviors handle logging, validation, transactions
- [ ] Behaviors registered in correct order (outermost first)
- [ ] Validation Behavior uses FluentValidation validators
- [ ] Exception handling is centralized

**Behavior Order:**
1. Logging (outermost)
2. Validation
3. Transaction
4. Handler (innermost)

### 11.5 API DTO Mapping

- [ ] Presentation layer maps API DTOs → Application Commands
- [ ] Application layer maps Domain → Response DTOs
- [ ] No direct use of Application Commands in HTTP binding
- [ ] API contracts are versioned separately

**Red Flags:**
- `CreateOrderCommand` in `[FromBody]` parameter
- Application Command structure dictated by API needs
- No mapping layer between Presentation and Application

### 11.6 Shared Code

- [ ] Feature-specific shared code in `Features/{Feature}/Shared/`
- [ ] Cross-feature shared code in `Application/Shared/`
- [ ] Domain logic never duplicated (stays in Domain Layer)
- [ ] Duplication extracted when 3+ occurrences

**Red Flags:**
- Business logic duplicated across slices
- Domain services inside slice folders
- No Shared folder structure

---

## 12. Code Organization

### 12.1 Project Structure

- [ ] Clear layer separation in project structure
- [ ] Domain project has no external dependencies
- [ ] Infrastructure project isolated
- [ ] Project names reflect architecture (not frameworks)

**Folder Structure:**
- [ ] Organized by feature/domain (not technical layers)
- [ ] Screaming architecture (intent clear from structure)
- [ ] Easy to find code by business concept

### 12.2 Naming Conventions

- [ ] Consistent naming across codebase
- [ ] Names reflect ubiquitous language
- [ ] No technical jargon in domain

---

## 13. Final Validation

### 13.1 The "Can We?" Questions

- [ ] Can we test business logic without database?
- [ ] Can we swap database without changing domain?
- [ ] Can we swap UI without changing domain?
- [ ] Can we swap framework without changing domain?
- [ ] Can we understand business from code structure?

### 13.2 The "Does It?" Questions

- [ ] Does architecture scream the domain (not the framework)?
- [ ] Does changing requirements only affect expected layers?
- [ ] Does code match how domain experts talk?
- [ ] Does adding a feature mean adding code (not changing)?

---

## Summary Score

**Critical Violations (❌):** _____
- Any violation here requires immediate attention

**Improvements Needed (⚠️):** _____
- Should be addressed in next refactoring

**Compliant (✅):** _____
- Percentage of checklist items passing

**Overall Assessment:**
- [ ] ✅ Architecture is sound
- [ ] ⚠️ Architecture needs improvements but is acceptable
- [ ] ❌ Architecture has critical violations requiring refactoring

---

## Code Examples

For working implementation examples:

**C# / .NET:**
- [Complete Feature Example](examples/csharp/CHECKLIST-EXAMPLES.md) - Full vertical slice
- [Vertical Slices Examples](examples/csharp/VERTICAL-SLICES-EXAMPLES.md) - Feature-based structure with proper layer separation
- [Project Structure Examples](examples/csharp/SKILL-EXAMPLES.md) - Solution organization
- [SOLID & DDD Examples](examples/csharp/PRINCIPLES-EXAMPLES.md) - Domain modeling
- [Design Pattern Examples](examples/csharp/PATTERNS-EXAMPLES.md) - Pattern implementations

## Resources

- [PRINCIPLES.md](PRINCIPLES.md) - Detailed explanation of principles
- [PATTERNS.md](PATTERNS.md) - Design patterns catalog
- [SKILL.md](SKILL.md) - Quick reference guide
- [Clean Code Skill](../clean-code/SKILL.md) - Code-level quality (naming, functions, comments)
