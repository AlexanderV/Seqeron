---
name: clean-architecture
version: 1.0.0
description: Guide for building software with Clean Architecture, DDD, SOLID principles, and hexagonal/onion patterns. Use when designing new features, refactoring code, reviewing architecture, or when user mentions clean code, architecture patterns, layers, separation of concerns, maintainability, testability, or asks for architectural guidance.
allowed-tools: Read, Grep, Glob, Edit, Write
triggers: [
  "clean architecture", "Uncle Bob", "Robert Martin",
  "dependency rule", "dependency inversion", "inversion of control", "composition root", "separation of concerns",
  "SOLID", "single responsibility", "open closed", "liskov", "interface segregation",
  "DDD", "domain-driven design", "domain driven",
  "aggregate root", "value object", "domain event", "bounded context",
  "domain service", "anemic domain", "anemic model", "rich domain model",
  "anti-corruption layer",
  "hexagonal", "ports and adapters", "onion architecture",
  "screaming architecture", "vertical slice",
  "presentation layer", "application layer", "domain layer", "infrastructure layer",
  "repository pattern", "unit of work", "specification pattern",
  "CQRS", "command query",
  "architecture review"
]
---

# Clean Architecture Development Guide

This skill helps you design and implement software following Clean Architecture principles, Domain-Driven Design (DDD), SOLID principles, and related architectural patterns (Hexagonal, Onion).

## When to Use This Skill

- Designing new features or modules
- Refactoring existing code to improve architecture
- Reviewing code for architectural issues
- Planning microservices or modular monoliths
- Setting up project structure
- Making architectural decisions

## Core Principles (The Foundation)

### 1. The Dependency Rule
**Dependencies must point inward only.**

```
┌─────────────────────────────────────┐
│   Infrastructure & Frameworks       │  ← Outer layer (most concrete)
│  ┌───────────────────────────────┐  │
│  │   Interface Adapters          │  │
│  │  ┌─────────────────────────┐  │  │
│  │  │   Application/Use Cases │  │  │
│  │  │  ┌───────────────────┐  │  │  │
│  │  │  │   Domain/Entities │  │  │  │  ← Inner layer (most abstract)
│  │  │  │   Business Rules  │  │  │  │
│  │  │  └───────────────────┘  │  │  │
│  │  └─────────────────────────┘  │  │
│  └───────────────────────────────┘  │
└─────────────────────────────────────┘
     Dependencies point inward ──→
```

**Key Rules:**
- Inner layers NEVER depend on outer layers
- Outer layers depend on inner layers through abstractions (interfaces)
- Business logic knows nothing about databases, UI, frameworks

### 2. SOLID Principles

- **S**ingle Responsibility: One class, one reason to change
- **O**pen/Closed: Open for extension, closed for modification
- **L**iskov Substitution: Subtypes must be substitutable for base types
- **I**nterface Segregation: Many small interfaces better than one large
- **D**ependency Inversion: Depend on abstractions, not concretions

### 3. Separation of Concerns

Each layer has a clear, single responsibility:
- **Domain**: Business rules and entities
- **Application**: Use cases and orchestration
- **Infrastructure**: Technical implementation details
- **Presentation**: User interface and APIs

## Architecture Implementation Steps

### Step 1: Identify the Domain

Before writing code, understand the business:

1. **Define Bounded Contexts** (DDD)
   - What are the logical boundaries?
   - What is the Ubiquitous Language?

2. **Identify Entities**
   - What are the core business objects?
   - What are their invariants (rules that must always be true)?

3. **Find Value Objects**
   - Immutable objects defined by their attributes
   - Examples: Money, Email, Address

4. **Discover Aggregates**
   - Clusters of entities with a root
   - Enforce consistency boundaries

### Step 2: Design the Layers

You have two main approaches for organizing the Application layer:

#### Option A: Traditional Horizontal Layers

```
Domain/
├── Entities/           # Business objects with identity
├── ValueObjects/       # Immutable objects
├── Aggregates/         # Aggregate roots
├── DomainServices/     # Domain logic that doesn't belong to entities
├── DomainEvents/       # Events that happened in the domain
├── Repositories/       # Interfaces (implementations in Infrastructure)
└── Specifications/     # Business rules encapsulation

Application/
├── UseCases/           # Application-specific business rules
├── DTOs/               # Data Transfer Objects
├── Interfaces/         # Ports (for adapters to implement)
├── Services/           # Application services
├── Commands/           # CQRS commands
├── Queries/            # CQRS queries
└── Validators/         # Input validation
```

#### Option B: Vertical Slices (Feature-Based) — Recommended for 3+ Developers

Organize Application layer by **features**, not technical concerns. **Important:** Domain Layer remains horizontal, Presentation Layer is separate.

```
src/
├── Domain/                          # HORIZONTAL - shared
│   └── Orders/
│       ├── Order.cs                 # Aggregate Root
│       └── OrderLine.cs             # Entity
│
├── Application/                     # VERTICAL SLICES
│   └── Features/
│       └── Orders/
│           ├── CreateOrder/
│           │   ├── CreateOrderCommand.cs
│           │   ├── CreateOrderHandler.cs
│           │   ├── CreateOrderValidator.cs
│           │   └── CreateOrderResponse.cs
│           └── Shared/
│               └── IOrderRepository.cs  # Port
│
├── Presentation/                    # SEPARATE from Application
│   └── Orders/
│       ├── CreateOrderEndpoint.cs
│       └── Contracts/
│           └── CreateOrderRequest.cs    # API DTO
│
└── Infrastructure/                  # HORIZONTAL
    └── Persistence/
        └── SqlOrderRepository.cs    # Adapter
```

**Vertical Slices Benefits:**
| Benefit | Description |
|---------|-------------|
| High Cohesion | All Application code for one feature in one folder |
| Layer Separation | Domain, Presentation, Infrastructure remain horizontal |
| Team Scalability | Multiple developers work without conflicts |
| Easy Navigation | Find everything in one place |
| Natural CQRS | Each slice is Command or Query |

**When to Use Vertical Slices:**
- Teams with 3+ developers
- Frequently changing requirements
- Projects using CQRS/MediatR
- Planning microservice extraction

**Important:** Domain Layer remains horizontal (shared), and Presentation Layer is separate from Application slices. See [Vertical Slices Principle](principles/10-vertical-slices.md) and [Code Examples](examples/csharp/VERTICAL-SLICES-EXAMPLES.md).

#### Layer 1: Domain (Core)
```
Domain/
├── Entities/           # Business objects with identity
├── ValueObjects/       # Immutable objects
├── Aggregates/         # Aggregate roots
├── DomainServices/     # Domain logic that doesn't belong to entities
├── DomainEvents/       # Events that happened in the domain
├── Repositories/       # Interfaces (implementations in Infrastructure)
└── Specifications/     # Business rules encapsulation
```

**Rules:**
- No dependencies on outer layers
- No framework code
- Pure business logic
- Framework-agnostic

#### Layer 2: Application (Use Cases)

**Traditional Structure:**
```
Application/
├── UseCases/           # Application-specific business rules
├── DTOs/               # Data Transfer Objects
├── Interfaces/         # Ports (for adapters to implement)
├── Services/           # Application services
├── Commands/           # CQRS commands
├── Queries/            # CQRS queries
└── Validators/         # Input validation
```

**Rules:**
- Orchestrates domain objects
- Implements use cases
- Depends only on Domain layer
- Defines interfaces for infrastructure

#### Layer 3: Infrastructure (Technical Details)
```
Infrastructure/
├── Persistence/        # Database implementations
│   ├── Repositories/   # Repository implementations
│   ├── EntityConfigurations/
│   └── Migrations/
├── ExternalServices/   # API clients, message queues
├── Identity/           # Authentication/Authorization
├── FileStorage/        # File system, cloud storage
└── Logging/            # Logging implementations
```

**Rules:**
- Implements interfaces from Application/Domain
- Contains all technical details
- Depends on Application and Domain through interfaces

#### Layer 4: Presentation (API/UI)
```
Presentation/
├── Controllers/        # API endpoints
├── ViewModels/         # UI-specific models
├── Middleware/         # HTTP pipeline
├── Filters/            # Cross-cutting concerns
└── Mappers/            # DTO to ViewModel mapping
```

**Rules:**
- Thin layer - delegates to Application
- No business logic
- Request/Response handling only

### Step 3: Apply Design Patterns

Use these patterns appropriately:

**Creational:**
- **Factory**: Create complex aggregates
- **Builder**: Construct complex objects step-by-step
- **Dependency Injection**: Invert dependencies

**Structural:**
- **Adapter**: Convert interfaces (Hexagonal Architecture)
- **Repository**: Abstract data access
- **Unit of Work**: Manage transactions

**Behavioral:**
- **Strategy**: Encapsulate algorithms
- **Command**: Encapsulate requests (CQRS)
- **Observer**: Domain events
- **Specification**: Business rules composition

### Step 4: Implement Dependency Inversion

The outer layer (Infrastructure) provides the concrete implementation at runtime through dependency injection.

**Violates DIP:** High-level code directly instantiates low-level concrete classes.
**Follows DIP:** High-level code depends on abstractions (interfaces), concrete implementations injected at runtime.

### Step 5: Testing Strategy

**Test Pyramid:**
```
        ┌─────┐
       /  E2E  \          ← Few, test critical paths
      /_________\
     /           \
    / Integration \        ← Moderate, test layer boundaries
   /_______________\
  /                 \
 /   Unit Tests      \       ← Many, test business logic
/_____________________\
```

**What to Test:**
- **Domain**: Unit test all business rules (most important!)
- **Application**: Unit test use cases, mock dependencies
- **Infrastructure**: Integration tests with real dependencies
- **Presentation**: Integration/E2E tests for critical flows

## Architecture Review Checklist

Before considering implementation complete, review using [CHECKLIST.md](CHECKLIST.md).

## Common Anti-Patterns to Avoid

### 1. Anemic Domain Model

**❌ BAD** - Business logic in service, entity is just data:
```csharp
public class Order
{
    public Guid Id { get; set; }
    public decimal Total { get; set; }
    public List<OrderItem> Items { get; set; }
}

public class OrderService
{
    public void AddItem(Order order, Product product, int quantity)
    {
        // Business logic in SERVICE instead of ENTITY
        var item = new OrderItem { ProductId = product.Id, Quantity = quantity };
        order.Items.Add(item);
        order.Total = order.Items.Sum(i => i.Subtotal);
    }
}
```

**✅ GOOD** - Business logic in entity:
```csharp
public class Order
{
    public OrderId Id { get; private set; }
    public Money Total { get; private set; }
    private readonly List<OrderItem> _items = new();

    public void AddItem(Product product, int quantity)
    {
        // Business logic in ENTITY
        if (Status != OrderStatus.Draft)
            throw new DomainException("Cannot add items to non-draft order");

        _items.Add(OrderItem.Create(product, quantity));
        RecalculateTotal();
    }
}
```

### 2. Leaking Domain to Presentation

**❌ BAD** - Domain entity exposed in API:
```csharp
[HttpGet]
public Order GetOrder(Guid id)
{
    return _orderRepository.GetById(id); // Exposes domain entity!
}
```

**✅ GOOD** - DTO isolates presentation from domain:
```csharp
[HttpGet]
public OrderDto GetOrder(Guid id)
{
    var order = _orderRepository.GetById(id);
    return _mapper.Map<OrderDto>(order); // Returns DTO
}
```

### 3. Infrastructure in Domain

**❌ BAD** - ORM annotations pollute domain:
```csharp
[Table("Orders")]
public class Order
{
    [Key]
    [Column("order_id")]
    public Guid Id { get; set; }
}
```

**✅ GOOD** - Configuration in Infrastructure:
```csharp
// Domain - pure C#
public class Order
{
    public OrderId Id { get; private set; }
}

// Infrastructure configuration
public class OrderConfiguration : IEntityTypeConfiguration<Order>
{
    public void Configure(EntityTypeBuilder<Order> builder)
    {
        builder.ToTable("Orders");
        builder.Property(o => o.Id).HasColumnName("order_id");
    }
}
```

### 4. God Services

**❌ BAD** - One service does everything:
```csharp
public class OrderService
{
    public void CreateOrder() { }
    public void CancelOrder() { }
    public void ShipOrder() { }
    public void ProcessPayment() { }
    public void SendEmail() { }
    public void GenerateInvoice() { }
    // ... 20 more methods
}
```

**✅ GOOD** - One handler per use case:
```csharp
public class CreateOrderCommandHandler : IRequestHandler<CreateOrderCommand>
public class CancelOrderCommandHandler : IRequestHandler<CancelOrderCommand>
public class ShipOrderCommandHandler : IRequestHandler<ShipOrderCommand>
```

### 5. Transaction Script

**❌ BAD** - Procedural, database-centric:
```csharp
public void PlaceOrder(Guid customerId, List<ItemDto> items)
{
    var order = new Order();
    _db.Orders.Add(order);
    _db.SaveChanges();

    foreach (var item in items)
    {
        var orderItem = new OrderItem { OrderId = order.Id, ProductId = item.ProductId };
        _db.OrderItems.Add(orderItem);
    }
    _db.SaveChanges();

    var total = _db.OrderItems.Where(i => i.OrderId == order.Id).Sum(i => i.Price);
    order.Total = total;
    _db.SaveChanges();
}
```

**✅ GOOD** - Object-oriented, domain-centric:
```csharp
public async Task<OrderId> Handle(PlaceOrderCommand request)
{
    // Create domain object
    var order = Order.Create(customerId, address);

    // Use domain logic
    foreach (var item in request.Items)
        order.AddItem(product, item.Quantity);

    // Persist through repository
    await _orderRepository.AddAsync(order);
    await _unitOfWork.SaveChangesAsync();

    return order.Id;
}
```

## Decision Framework

When making architectural decisions, ask:

1. **Does this follow the Dependency Rule?**
   - Inner layers independent of outer?

2. **Is business logic in the domain?**
   - Not scattered in services or controllers?

3. **Are dependencies inverted?**
   - Using interfaces, not concrete types?

4. **Is it testable?**
   - Can I test without database/network?

5. **Is it simple?**
   - Not over-engineered for current needs?

## Reference Materials

- **[PRINCIPLES.md](PRINCIPLES.md)** - Detailed explanation of all architectural principles
- **[PATTERNS.md](PATTERNS.md)** - Design patterns catalog
- **[CHECKLIST.md](CHECKLIST.md)** - Complete architecture review checklist
- **[Clean Code Skill](../clean-code/SKILL.md)** - Code-level quality (naming, functions, error handling)

## Quick Decision Tree

```
Need to add new feature?
├─ Using Vertical Slices?
│  ├─ Yes → Create new folder: Features/{Domain}/{FeatureName}/
│  │        └─ Add: Command/Query, Handler, Validator, Response, Endpoint
│  └─ No → Traditional layers:
│     ├─ Is it business logic?
│     │  ├─ Yes → Add to Domain layer
│     │  └─ No → Is it orchestration?
│     │     ├─ Yes → Add to Application layer
│     │     └─ No → Is it technical?
│     │        ├─ Yes → Add to Infrastructure
│     │        └─ No → Add to Presentation
│
How to structure Application layer?
├─ Team size 3+? → Vertical Slices (feature-based)
├─ Frequent changes? → Vertical Slices
├─ Simple CRUD? → Traditional layers may suffice
└─ Complex domain? → Hybrid (DDD domain + Vertical Slices application)
│
Need to access data?
├─ Define interface in Domain (or Features/{Domain}/Shared/)
├─ Implement in Infrastructure
└─ Inject into Application/Handler
│
Need external service?
├─ Define port (interface) in Application
├─ Create adapter in Infrastructure
└─ Use dependency injection
│
Found business rule?
├─ Belongs to one entity? → Method on entity
├─ Spans entities? → Domain Service
└─ Application-specific? → Use Case / Handler
│
Cross-cutting concern (logging, validation)?
├─ Using MediatR? → Pipeline Behavior in Shared/Behaviors/
└─ Traditional? → Middleware or Filters
```

## Working with This Skill

When I use this skill, I will:

1. **Analyze** the current architecture
2. **Identify** violations of principles
3. **Suggest** improvements aligned with Clean Architecture
4. **Implement** changes following the patterns
5. **Validate** with the checklist
6. **Ensure** proper testing strategy

## Code Examples

For complete runnable code examples:

**C# / .NET Examples:**
- [Project Structure Examples](examples/csharp/SKILL-EXAMPLES.md) - Solution structure, basic layer examples
- [SOLID & DDD Examples](examples/csharp/PRINCIPLES-EXAMPLES.md) - Entities, Value Objects, Domain Events
- [Design Pattern Examples](examples/csharp/PATTERNS-EXAMPLES.md) - Repository, Unit of Work, CQRS, Observers
- [Complete Feature Example](examples/csharp/CHECKLIST-EXAMPLES.md) - Full vertical slice example
- [Error Handling Examples](examples/csharp/ERROR-HANDLING-EXAMPLES.md) - Exception strategy, Result pattern

## Getting Started

To apply Clean Architecture to your code:

1. Review existing code structure
2. Identify current layer violations
3. Create proper layer structure
4. Move code to appropriate layers
5. Apply dependency inversion
6. Add tests for business logic
7. Validate with checklist

Let's build maintainable, testable, and scalable software!
