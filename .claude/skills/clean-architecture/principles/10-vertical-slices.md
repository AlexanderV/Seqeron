# Vertical Slices (Feature-Based Structure)

> **Navigation:** [← Back to Principles](../PRINCIPLES.md) | [Code Examples →](../examples/csharp/VERTICAL-SLICES-EXAMPLES.md)

## Definition

Vertical Slices is an architectural pattern where code is organized by **business features** rather than technical layers. Each "slice" contains the **Application Layer components** needed to implement one use case, while respecting Clean Architecture layer boundaries.

## Core Principle

**One use case = One folder containing Application Layer code for that feature.**

Vertical Slices organize the **Application Layer only**. The Domain Layer remains horizontal (shared across features), and the Presentation Layer stays separate.

## Visual Comparison

```
┌─────────────────────────────────────────────────────────┐
│              HORIZONTAL (Traditional Layers)            │
├─────────────────────────────────────────────────────────┤
│  Controllers/     Services/        Repositories/        │
│  ────────────     ─────────        ─────────────        │
│  UserController   UserService      UserRepository       │
│  OrderController  OrderService     OrderRepository      │
│  PaymentCtrl      PaymentService   PaymentRepository    │
│                                                         │
│  To implement CreateOrder: touch 5+ files in 5 folders  │
└─────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────┐
│               VERTICAL SLICES (Feature-Based)           │
├─────────────────────────────────────────────────────────┤
│   Application/Features/Orders/CreateOrder/              │
│   ├── CreateOrderCommand.cs      (Application Layer)    │
│   ├── CreateOrderHandler.cs      (Application Layer)    │
│   └── CreateOrderValidator.cs    (Application Layer)    │
│                                                         │
│   Presentation/Orders/                                  │
│   └── CreateOrderEndpoint.cs     (Presentation Layer)   │
│                                                         │
│   To implement CreateOrder: touch 1-2 folders           │
└─────────────────────────────────────────────────────────┘
```

## Integration with Clean Architecture Layers

Vertical Slices organize the **Application Layer**, while Clean Architecture provides **layer separation**:

```
┌─────────────────────────────────────────────────────────┐
│ PRESENTATION LAYER                                      │
│ └── Endpoints/Controllers (per feature)                 │
│     └── Maps API DTOs → Application Commands            │
├─────────────────────────────────────────────────────────┤
│ APPLICATION LAYER (VERTICAL SLICES)                     │
│ ├── Features/Orders/CreateOrder/                        │
│ │   ├── CreateOrderCommand.cs   ← Request DTO           │
│ │   ├── CreateOrderHandler.cs   ← Orchestrates Domain   │
│ │   ├── CreateOrderValidator.cs ← Input validation      │
│ │   └── CreateOrderResponse.cs  ← Response DTO          │
│ └── Features/Orders/Shared/                             │
│     └── IOrderRepository.cs     ← Port (interface)      │
├─────────────────────────────────────────────────────────┤
│ DOMAIN LAYER (HORIZONTAL - DDD)                         │
│ ├── Orders/                                             │
│ │   ├── Order.cs                ← Aggregate Root        │
│ │   ├── OrderLine.cs            ← Entity                │
│ │   ├── OrderStatus.cs          ← Value Object          │
│ │   └── OrderCreatedEvent.cs    ← Domain Event          │
│ ├── Customers/                                          │
│ └── Shared/                     ← Base classes          │
├─────────────────────────────────────────────────────────┤
│ INFRASTRUCTURE LAYER                                    │
│ └── Persistence/Orders/                                 │
│     └── SqlOrderRepository.cs   ← Implements Port       │
└─────────────────────────────────────────────────────────┘
```

**Key Points:**
- Domain Layer is **shared horizontally** across all features
- Application Layer is organized by **vertical slices**
- Presentation Layer has **separate endpoints** that map to slices
- Infrastructure Layer **implements ports** defined in Application

## Anatomy of a Vertical Slice

Each slice contains **Application Layer components only**:

| Component | Layer | Purpose |
|-----------|-------|---------|
| **Command/Query** | Application | Request DTO with intent |
| **Handler** | Application | Orchestrates domain objects |
| **Validator** | Application | Input validation rules |
| **Response** | Application | Response DTO |
| **Endpoint** | Presentation | HTTP mapping (separate folder) |

## Correct Folder Structure

```
src/
├── Domain/                              # HORIZONTAL - shared across features
│   ├── Orders/
│   │   ├── Order.cs                     # Aggregate Root
│   │   ├── OrderLine.cs                 # Entity
│   │   ├── OrderStatus.cs               # Value Object
│   │   └── Events/
│   │       └── OrderCreatedEvent.cs     # Domain Event
│   ├── Customers/
│   │   └── Customer.cs
│   └── Shared/
│       ├── AggregateRoot.cs
│       ├── Entity.cs
│       └── IDomainEvent.cs
│
├── Application/                         # VERTICAL SLICES
│   ├── Features/
│   │   ├── Orders/
│   │   │   ├── CreateOrder/
│   │   │   │   ├── CreateOrderCommand.cs
│   │   │   │   ├── CreateOrderHandler.cs
│   │   │   │   ├── CreateOrderValidator.cs
│   │   │   │   └── CreateOrderResponse.cs
│   │   │   ├── GetOrderById/
│   │   │   │   ├── GetOrderByIdQuery.cs
│   │   │   │   ├── GetOrderByIdHandler.cs
│   │   │   │   └── GetOrderByIdResponse.cs
│   │   │   └── Shared/                  # Feature-specific ports
│   │   │       └── IOrderRepository.cs
│   │   │
│   │   └── Payments/
│   │       ├── ProcessPayment/
│   │       └── Shared/
│   │
│   ├── Shared/                          # Cross-feature Application concerns
│   │   ├── Behaviors/                   # MediatR Pipeline Behaviors
│   │   │   ├── LoggingBehavior.cs
│   │   │   ├── ValidationBehavior.cs
│   │   │   └── TransactionBehavior.cs
│   │   └── Interfaces/
│   │       ├── IUnitOfWork.cs
│   │       └── IDomainEventDispatcher.cs
│   │
│   └── DependencyInjection.cs           # Application layer DI registration
│
├── Infrastructure/                      # HORIZONTAL - implements ports
│   ├── Persistence/
│   │   ├── Orders/
│   │   │   └── SqlOrderRepository.cs
│   │   ├── AppDbContext.cs
│   │   └── UnitOfWork.cs
│   ├── Messaging/
│   │   └── DomainEventDispatcher.cs
│   └── DependencyInjection.cs
│
├── Presentation/                        # SEPARATE from Application
│   ├── Orders/
│   │   ├── CreateOrderEndpoint.cs
│   │   ├── GetOrderByIdEndpoint.cs
│   │   └── Contracts/                   # API-specific DTOs
│   │       ├── CreateOrderRequest.cs
│   │       └── OrderResponse.cs
│   ├── Payments/
│   └── DependencyInjection.cs
│
└── Program.cs
```

## Key Benefits

| Benefit | Description |
|---------|-------------|
| **High Cohesion** | All Application Layer code for a feature in one location |
| **Low Coupling** | Features are independent of each other |
| **Easy Navigation** | Find all use case code in one folder |
| **Team Scalability** | Multiple developers work on different features without conflicts |
| **Natural CQRS** | Each slice is either a Command or Query |
| **Clear Boundaries** | Respects Clean Architecture layer separation |

## Important Clarifications

### Code Duplication Policy

**Minimize duplication through Shared folders, but prefer duplication over wrong abstraction.**

- **Domain logic**: Never duplicate. Keep in Domain Layer.
- **Application orchestration**: May have similar patterns, extract to base handlers if 3+ duplicates.
- **Validation rules**: Shared rules go to Shared/Validators.
- **DTOs**: Each slice has its own DTOs (not shared).

### Microservice Extraction

Vertical Slices **ease** microservice extraction but don't solve:

- Bounded Context boundaries (requires DDD analysis)
- Database separation (requires data partitioning strategy)
- Distributed transactions (requires Saga/Outbox patterns)
- Service communication (requires API contracts/events)

Vertical Slices help because feature code is already grouped, but migration still requires careful planning.

## Comparison with Clean Architecture

| Aspect | Clean Architecture | Vertical Slices + Clean Architecture |
|--------|-------------------|--------------------------------------|
| **Domain Layer** | Horizontal | Horizontal (unchanged) |
| **Application Layer** | By technical concern | By business feature |
| **Presentation Layer** | Horizontal | Parallel to features |
| **Navigation** | Jump between folders | Stay in feature folder |
| **Adding Feature** | Modify multiple folders | Add new slice folder |
| **Dependency Rule** | Enforced | Enforced (same) |

## When to Use Vertical Slices

### ✅ Good Fit

- CRUD applications with many use cases
- Teams with 3+ developers
- Frequently changing requirements
- Projects already using CQRS/MediatR
- Modular monoliths
- Clear use case boundaries

### ❌ Not Recommended

- Simple services with 2-3 endpoints
- Solo developer projects
- Libraries or SDKs
- Systems where features heavily share logic
- Legacy systems without clear use case boundaries

## Trade-offs and Mitigations

| Trade-off | Mitigation |
|-----------|------------|
| Similar code in slices | Extract to Shared when 3+ duplicates appear |
| Harder to see "all handlers" | IDE search, naming conventions |
| Defining feature boundaries | 1 API endpoint = 1 Feature (rule of thumb) |
| Team learning curve | Start with new features, gradual adoption |
| Domain logic location | Always in Domain Layer, never in handlers |

## Decision Tree

```
Adding new functionality?
│
├── Is it business logic (invariants, calculations)?
│   └── YES → Add to Domain Layer (Entity/Aggregate/Service)
│
├── Is it a new use case (API endpoint)?
│   └── YES → Create new slice in Application/Features/
│
├── Is it shared across 3+ features?
│   └── YES → Put in Application/Shared/ or Domain/
│
└── Is it technical (DB, messaging, caching)?
    └── YES → Add to Infrastructure Layer
```

## Related Principles

- [CQRS Pattern](06-cqrs-pattern.md) - Vertical Slices naturally implement CQRS
- [Screaming Architecture](07-screaming-architecture.md) - Features scream business intent
- [Hexagonal Architecture](04-hexagonal-architecture.md) - Ports defined in Application, Adapters in Infrastructure
- [Single Responsibility](02-solid-principles.md#s---single-responsibility-principle-srp) - Each handler has one responsibility
- [Dependency Rule](01-dependency-rule.md) - Dependencies still point inward

## Code Examples

For complete implementation examples with Domain Layer integration, see:
- [Vertical Slices Examples](../examples/csharp/VERTICAL-SLICES-EXAMPLES.md) - Full slice with proper layer separation
