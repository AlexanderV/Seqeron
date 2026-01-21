# Vertical Slices - Code Examples

> **Navigation:** [← Back to Vertical Slices Principle](../../principles/10-vertical-slices.md) | [Other Examples →](./PATTERNS-EXAMPLES.md)

This file contains complete code examples demonstrating Vertical Slices with proper Clean Architecture layer separation.

---

## Table of Contents

1. [Complete Slice Example: CreateOrder](#1-complete-slice-example-createorder)
2. [Domain Layer (Horizontal)](#2-domain-layer-horizontal)
3. [Application Layer (Vertical Slice)](#3-application-layer-vertical-slice)
4. [Presentation Layer (Separate)](#4-presentation-layer-separate)
5. [Infrastructure Layer](#5-infrastructure-layer)
6. [Cross-Cutting Concerns (Behaviors)](#6-cross-cutting-concerns-behaviors)
7. [Dependency Injection Setup](#7-dependency-injection-setup)
8. [Testing a Vertical Slice](#8-testing-a-vertical-slice)
9. [Query Example: GetOrderById](#9-query-example-getorderbyid)

---

## 1. Complete Slice Example: CreateOrder

This example shows how all layers work together for a single use case.

### Layer Diagram

```
┌─────────────────────────────────────────────────────────────┐
│ Presentation: CreateOrderEndpoint                           │
│   └── Maps CreateOrderRequest → CreateOrderCommand          │
├─────────────────────────────────────────────────────────────┤
│ Application: CreateOrder/ (Vertical Slice)                  │
│   ├── CreateOrderCommand     ← Application DTO              │
│   ├── CreateOrderHandler     ← Orchestrates Domain          │
│   ├── CreateOrderValidator   ← FluentValidation             │
│   └── CreateOrderResponse    ← Response DTO                 │
├─────────────────────────────────────────────────────────────┤
│ Domain: Orders/                                             │
│   ├── Order (Aggregate Root) ← Business logic               │
│   ├── OrderLine (Entity)     ← Part of Aggregate            │
│   ├── OrderStatus (Value Object)                            │
│   └── OrderCreatedEvent      ← Domain Event                 │
├─────────────────────────────────────────────────────────────┤
│ Infrastructure:                                             │
│   ├── SqlOrderRepository     ← Implements IOrderRepository  │
│   └── DomainEventDispatcher  ← Publishes Domain Events      │
└─────────────────────────────────────────────────────────────┘
```

---

## 2. Domain Layer (Horizontal)

Domain Layer contains business logic and is **shared across all features**.

### 2.1 Base Classes

```csharp
// Domain/Shared/Entity.cs
namespace Domain.Shared;

public abstract class Entity<TId> where TId : notnull
{
    public TId Id { get; protected set; } = default!;
    
    public override bool Equals(object? obj)
    {
        if (obj is not Entity<TId> other) return false;
        if (ReferenceEquals(this, other)) return true;
        if (GetType() != other.GetType()) return false;
        return Id.Equals(other.Id);
    }
    
    public override int GetHashCode() => Id.GetHashCode();
}
```

```csharp
// Domain/Shared/AggregateRoot.cs
namespace Domain.Shared;

public abstract class AggregateRoot<TId> : Entity<TId> where TId : notnull
{
    private readonly List<IDomainEvent> _domainEvents = new();
    
    public IReadOnlyList<IDomainEvent> DomainEvents => _domainEvents.AsReadOnly();
    
    protected void AddDomainEvent(IDomainEvent domainEvent)
    {
        _domainEvents.Add(domainEvent);
    }
    
    public void ClearDomainEvents()
    {
        _domainEvents.Clear();
    }
}
```

```csharp
// Domain/Shared/IDomainEvent.cs
namespace Domain.Shared;

public interface IDomainEvent
{
    DateTime OccurredAt { get; }
}
```

```csharp
// Domain/Shared/Result.cs
namespace Domain.Shared;

public class Result<T>
{
    public bool IsSuccess { get; }
    public bool IsFailure => !IsSuccess;
    public T? Value { get; }
    public string Error { get; }
    
    private Result(bool isSuccess, T? value, string error)
    {
        IsSuccess = isSuccess;
        Value = value;
        Error = error;
    }
    
    public static Result<T> Success(T value) => new(true, value, string.Empty);
    public static Result<T> Failure(string error) => new(false, default, error);
    
    public TResult Match<TResult>(
        Func<T, TResult> onSuccess,
        Func<string, TResult> onFailure) =>
        IsSuccess ? onSuccess(Value!) : onFailure(Error);
}
```

### 2.2 Value Objects

```csharp
// Domain/Orders/OrderId.cs
namespace Domain.Orders;

public readonly record struct OrderId(Guid Value)
{
    public static OrderId New() => new(Guid.NewGuid());
    public static OrderId From(Guid value) => new(value);
    public override string ToString() => Value.ToString();
}
```

```csharp
// Domain/Orders/OrderStatus.cs
namespace Domain.Orders;

public enum OrderStatus
{
    Pending,
    Confirmed,
    Shipped,
    Delivered,
    Cancelled
}
```

```csharp
// Domain/Shared/Money.cs
namespace Domain.Shared;

public readonly record struct Money
{
    public decimal Amount { get; }
    public string Currency { get; }
    
    private Money(decimal amount, string currency)
    {
        if (amount < 0)
            throw new ArgumentException("Amount cannot be negative", nameof(amount));
        if (string.IsNullOrWhiteSpace(currency))
            throw new ArgumentException("Currency is required", nameof(currency));
            
        Amount = amount;
        Currency = currency.ToUpperInvariant();
    }
    
    public static Money Create(decimal amount, string currency = "USD") => 
        new(amount, currency);
    
    public static Money Zero(string currency = "USD") => 
        new(0, currency);
    
    public Money Add(Money other)
    {
        if (Currency != other.Currency)
            throw new InvalidOperationException("Cannot add money with different currencies");
        return new Money(Amount + other.Amount, Currency);
    }
    
    public Money Multiply(int quantity) => 
        new(Amount * quantity, Currency);
}
```

### 2.3 Order Aggregate Root

```csharp
// Domain/Orders/Order.cs
namespace Domain.Orders;

using Domain.Shared;
using Domain.Orders.Events;

public class Order : AggregateRoot<OrderId>
{
    private readonly List<OrderLine> _items = new();
    
    public CustomerId CustomerId { get; private set; }
    public OrderStatus Status { get; private set; }
    public Money TotalAmount { get; private set; }
    public DateTime CreatedAt { get; private set; }
    public IReadOnlyList<OrderLine> Items => _items.AsReadOnly();
    
    // Private constructor - enforces factory method usage
    private Order() { }
    
    /// <summary>
    /// Factory method with business invariants
    /// </summary>
    public static Result<Order> Create(
        CustomerId customerId, 
        IEnumerable<(ProductId ProductId, int Quantity, Money UnitPrice)> items)
    {
        // Business Invariant: Customer is required
        if (customerId == default)
            return Result<Order>.Failure("Customer ID is required");
        
        var itemsList = items.ToList();
        
        // Business Invariant: Order must have at least one item
        if (!itemsList.Any())
            return Result<Order>.Failure("Order must contain at least one item");
        
        // Business Invariant: All quantities must be positive
        if (itemsList.Any(i => i.Quantity <= 0))
            return Result<Order>.Failure("All item quantities must be positive");
        
        var order = new Order
        {
            Id = OrderId.New(),
            CustomerId = customerId,
            Status = OrderStatus.Pending,
            CreatedAt = DateTime.UtcNow,
            TotalAmount = Money.Zero()
        };
        
        // Add items through domain logic
        foreach (var (productId, quantity, unitPrice) in itemsList)
        {
            var lineResult = OrderLine.Create(productId, quantity, unitPrice);
            if (lineResult.IsFailure)
                return Result<Order>.Failure(lineResult.Error);
                
            order._items.Add(lineResult.Value!);
        }
        
        // Calculate total
        order.RecalculateTotal();
        
        // Raise domain event
        order.AddDomainEvent(new OrderCreatedEvent(order.Id, order.CustomerId, order.TotalAmount));
        
        return Result<Order>.Success(order);
    }
    
    /// <summary>
    /// Confirm order - business operation with invariants
    /// </summary>
    public Result<Order> Confirm()
    {
        // Business Invariant: Only pending orders can be confirmed
        if (Status != OrderStatus.Pending)
            return Result<Order>.Failure($"Cannot confirm order in {Status} status");
        
        Status = OrderStatus.Confirmed;
        AddDomainEvent(new OrderConfirmedEvent(Id));
        
        return Result<Order>.Success(this);
    }
    
    /// <summary>
    /// Cancel order - business operation with invariants
    /// </summary>
    public Result<Order> Cancel(string reason)
    {
        // Business Invariant: Cannot cancel delivered or already cancelled orders
        if (Status == OrderStatus.Delivered)
            return Result<Order>.Failure("Cannot cancel delivered order");
        if (Status == OrderStatus.Cancelled)
            return Result<Order>.Failure("Order is already cancelled");
        
        Status = OrderStatus.Cancelled;
        AddDomainEvent(new OrderCancelledEvent(Id, reason));
        
        return Result<Order>.Success(this);
    }
    
    private void RecalculateTotal()
    {
        TotalAmount = _items.Aggregate(
            Money.Zero(),
            (total, line) => total.Add(line.Subtotal));
    }
}
```

### 2.4 OrderLine Entity

```csharp
// Domain/Orders/OrderLine.cs
namespace Domain.Orders;

using Domain.Shared;

public class OrderLine : Entity<Guid>
{
    public ProductId ProductId { get; private set; }
    public int Quantity { get; private set; }
    public Money UnitPrice { get; private set; }
    public Money Subtotal => UnitPrice.Multiply(Quantity);
    
    private OrderLine() { }
    
    public static Result<OrderLine> Create(ProductId productId, int quantity, Money unitPrice)
    {
        if (productId == default)
            return Result<OrderLine>.Failure("Product ID is required");
        if (quantity <= 0)
            return Result<OrderLine>.Failure("Quantity must be positive");
        if (unitPrice.Amount <= 0)
            return Result<OrderLine>.Failure("Unit price must be positive");
        
        return Result<OrderLine>.Success(new OrderLine
        {
            Id = Guid.NewGuid(),
            ProductId = productId,
            Quantity = quantity,
            UnitPrice = unitPrice
        });
    }
}
```

### 2.5 Domain Events

```csharp
// Domain/Orders/Events/OrderCreatedEvent.cs
namespace Domain.Orders.Events;

using Domain.Shared;

public sealed record OrderCreatedEvent(
    OrderId OrderId,
    CustomerId CustomerId,
    Money TotalAmount) : IDomainEvent
{
    public DateTime OccurredAt { get; } = DateTime.UtcNow;
}
```

```csharp
// Domain/Orders/Events/OrderConfirmedEvent.cs
namespace Domain.Orders.Events;

using Domain.Shared;

public sealed record OrderConfirmedEvent(OrderId OrderId) : IDomainEvent
{
    public DateTime OccurredAt { get; } = DateTime.UtcNow;
}
```

```csharp
// Domain/Orders/Events/OrderCancelledEvent.cs
namespace Domain.Orders.Events;

using Domain.Shared;

public sealed record OrderCancelledEvent(
    OrderId OrderId, 
    string Reason) : IDomainEvent
{
    public DateTime OccurredAt { get; } = DateTime.UtcNow;
}
```

---

## 3. Application Layer (Vertical Slice)

### 3.1 Repository Interface (Port)

```csharp
// Application/Features/Orders/Shared/IOrderRepository.cs
namespace Application.Features.Orders.Shared;

using Domain.Orders;

/// <summary>
/// Port for Order persistence - defined in Application, implemented in Infrastructure
/// </summary>
public interface IOrderRepository
{
    Task<Order?> GetByIdAsync(OrderId id, CancellationToken ct = default);
    Task AddAsync(Order order, CancellationToken ct = default);
    Task UpdateAsync(Order order, CancellationToken ct = default);
}
```

### 3.2 CreateOrder Command

```csharp
// Application/Features/Orders/CreateOrder/CreateOrderCommand.cs
namespace Application.Features.Orders.CreateOrder;

using MediatR;
using Domain.Shared;

/// <summary>
/// Application-level Command DTO
/// Contains data needed to execute the use case
/// </summary>
public sealed record CreateOrderCommand(
    Guid CustomerId,
    List<CreateOrderItemDto> Items
) : IRequest<Result<CreateOrderResponse>>;

public sealed record CreateOrderItemDto(
    Guid ProductId,
    int Quantity,
    decimal UnitPrice,
    string Currency = "USD"
);
```

### 3.3 CreateOrder Response

```csharp
// Application/Features/Orders/CreateOrder/CreateOrderResponse.cs
namespace Application.Features.Orders.CreateOrder;

/// <summary>
/// Application-level Response DTO
/// </summary>
public sealed record CreateOrderResponse(
    Guid OrderId,
    decimal TotalAmount,
    string Currency,
    DateTime CreatedAt
);
```

### 3.4 CreateOrder Validator

```csharp
// Application/Features/Orders/CreateOrder/CreateOrderValidator.cs
namespace Application.Features.Orders.CreateOrder;

using FluentValidation;

/// <summary>
/// Input validation - runs BEFORE handler via Pipeline Behavior
/// Validates Application-level concerns (presence, format)
/// Business invariants are checked in Domain Layer
/// </summary>
public sealed class CreateOrderValidator : AbstractValidator<CreateOrderCommand>
{
    public CreateOrderValidator()
    {
        RuleFor(x => x.CustomerId)
            .NotEmpty()
            .WithMessage("Customer ID is required");
        
        RuleFor(x => x.Items)
            .NotEmpty()
            .WithMessage("Order must contain at least one item");
        
        RuleForEach(x => x.Items).ChildRules(item =>
        {
            item.RuleFor(x => x.ProductId)
                .NotEmpty()
                .WithMessage("Product ID is required");
            
            item.RuleFor(x => x.Quantity)
                .GreaterThan(0)
                .WithMessage("Quantity must be positive");
            
            item.RuleFor(x => x.UnitPrice)
                .GreaterThan(0)
                .WithMessage("Unit price must be positive");
            
            item.RuleFor(x => x.Currency)
                .NotEmpty()
                .Length(3)
                .WithMessage("Currency must be a 3-letter code");
        });
    }
}
```

### 3.5 CreateOrder Handler

```csharp
// Application/Features/Orders/CreateOrder/CreateOrderHandler.cs
namespace Application.Features.Orders.CreateOrder;

using MediatR;
using Domain.Orders;
using Domain.Shared;
using Application.Features.Orders.Shared;
using Application.Shared.Interfaces;

/// <summary>
/// Application Service / Use Case Handler
/// 
/// Responsibilities:
/// 1. Load/verify external data (customer existence)
/// 2. Call Domain Layer for business logic
/// 3. Persist the aggregate
/// 4. Dispatch domain events
/// 5. Return response
/// 
/// Does NOT contain business logic - that's in Domain Layer
/// </summary>
public sealed class CreateOrderHandler 
    : IRequestHandler<CreateOrderCommand, Result<CreateOrderResponse>>
{
    private readonly IOrderRepository _orderRepository;
    private readonly ICustomerRepository _customerRepository;
    private readonly IUnitOfWork _unitOfWork;
    private readonly IDomainEventDispatcher _eventDispatcher;
    
    public CreateOrderHandler(
        IOrderRepository orderRepository,
        ICustomerRepository customerRepository,
        IUnitOfWork unitOfWork,
        IDomainEventDispatcher eventDispatcher)
    {
        _orderRepository = orderRepository;
        _customerRepository = customerRepository;
        _unitOfWork = unitOfWork;
        _eventDispatcher = eventDispatcher;
    }
    
    public async Task<Result<CreateOrderResponse>> Handle(
        CreateOrderCommand request, 
        CancellationToken ct)
    {
        // 1. Verify customer exists (Application concern)
        var customerExists = await _customerRepository.ExistsAsync(
            CustomerId.From(request.CustomerId), ct);
        
        if (!customerExists)
            return Result<CreateOrderResponse>.Failure("Customer not found");
        
        // 2. Prepare domain data
        var orderItems = request.Items.Select(item => (
            ProductId: ProductId.From(item.ProductId),
            Quantity: item.Quantity,
            UnitPrice: Money.Create(item.UnitPrice, item.Currency)
        ));
        
        // 3. Call Domain Layer - business logic is HERE
        var orderResult = Order.Create(
            CustomerId.From(request.CustomerId),
            orderItems
        );
        
        // 4. Handle domain validation failure
        if (orderResult.IsFailure)
            return Result<CreateOrderResponse>.Failure(orderResult.Error);
        
        var order = orderResult.Value!;
        
        // 5. Persist aggregate
        await _orderRepository.AddAsync(order, ct);
        await _unitOfWork.SaveChangesAsync(ct);
        
        // 6. Dispatch domain events AFTER successful persistence
        await _eventDispatcher.DispatchAsync(order.DomainEvents, ct);
        order.ClearDomainEvents();
        
        // 7. Return response
        return Result<CreateOrderResponse>.Success(new CreateOrderResponse(
            order.Id.Value,
            order.TotalAmount.Amount,
            order.TotalAmount.Currency,
            order.CreatedAt
        ));
    }
}
```

### 3.6 Application Shared Interfaces

```csharp
// Application/Shared/Interfaces/IUnitOfWork.cs
namespace Application.Shared.Interfaces;

public interface IUnitOfWork
{
    Task<int> SaveChangesAsync(CancellationToken ct = default);
}
```

```csharp
// Application/Shared/Interfaces/IDomainEventDispatcher.cs
namespace Application.Shared.Interfaces;

using Domain.Shared;

public interface IDomainEventDispatcher
{
    Task DispatchAsync(IEnumerable<IDomainEvent> domainEvents, CancellationToken ct = default);
}
```

```csharp
// Application/Shared/Interfaces/ICustomerRepository.cs
namespace Application.Shared.Interfaces;

using Domain.Orders;

public interface ICustomerRepository
{
    Task<bool> ExistsAsync(CustomerId id, CancellationToken ct = default);
}
```

---

## 4. Presentation Layer (Separate)

### 4.1 API Request DTO

```csharp
// Presentation/Orders/Contracts/CreateOrderRequest.cs
namespace Presentation.Orders.Contracts;

/// <summary>
/// API-level DTO - external contract
/// Separate from Application Command to allow API evolution
/// </summary>
public sealed record CreateOrderRequest(
    Guid CustomerId,
    List<OrderItemRequest> Items
);

public sealed record OrderItemRequest(
    Guid ProductId,
    int Quantity,
    decimal Price,
    string Currency = "USD"
);
```

### 4.2 API Response DTO

```csharp
// Presentation/Orders/Contracts/OrderResponse.cs
namespace Presentation.Orders.Contracts;

/// <summary>
/// API-level Response - external contract
/// </summary>
public sealed record OrderResponse(
    Guid Id,
    decimal Total,
    string Currency,
    string Status,
    DateTime CreatedAt
);
```

### 4.3 Endpoint with Proper Mapping

```csharp
// Presentation/Orders/CreateOrderEndpoint.cs
namespace Presentation.Orders;

using Microsoft.AspNetCore.Http;
using Microsoft.AspNetCore.Routing;
using Microsoft.AspNetCore.Builder;
using MediatR;
using Presentation.Orders.Contracts;
using Application.Features.Orders.CreateOrder;

/// <summary>
/// Presentation Layer Endpoint
/// 
/// Responsibilities:
/// 1. Receive HTTP request
/// 2. Map API DTO to Application Command
/// 3. Send to MediatR
/// 4. Map Application Response to API Response
/// 5. Return HTTP response
/// </summary>
public static class CreateOrderEndpoint
{
    public static void Map(IEndpointRouteBuilder app)
    {
        app.MapPost("/api/orders", Handle)
            .WithTags("Orders")
            .WithName("CreateOrder")
            .WithDescription("Creates a new order for a customer")
            .Produces<OrderResponse>(StatusCodes.Status201Created)
            .Produces(StatusCodes.Status400BadRequest)
            .Produces(StatusCodes.Status404NotFound);
    }
    
    private static async Task<IResult> Handle(
        CreateOrderRequest request,     // ← API DTO (external contract)
        IMediator mediator,
        CancellationToken ct)
    {
        // 1. Map API DTO → Application Command
        var command = MapToCommand(request);
        
        // 2. Send to Application Layer
        var result = await mediator.Send(command, ct);
        
        // 3. Map result to HTTP response
        return result.Match(
            success => Results.Created(
                $"/api/orders/{success.OrderId}",
                MapToResponse(success)),
            error => error.Contains("not found") 
                ? Results.NotFound(error) 
                : Results.BadRequest(error)
        );
    }
    
    /// <summary>
    /// Maps API DTO to Application Command
    /// Isolates external contract from internal structure
    /// </summary>
    private static CreateOrderCommand MapToCommand(CreateOrderRequest request) =>
        new(
            request.CustomerId,
            request.Items.Select(i => new CreateOrderItemDto(
                i.ProductId,
                i.Quantity,
                i.Price,
                i.Currency
            )).ToList()
        );
    
    /// <summary>
    /// Maps Application Response to API Response
    /// </summary>
    private static OrderResponse MapToResponse(CreateOrderResponse response) =>
        new(
            response.OrderId,
            response.TotalAmount,
            response.Currency,
            "Pending",
            response.CreatedAt
        );
}
```

### 4.4 Endpoint Registration

```csharp
// Presentation/DependencyInjection.cs
namespace Presentation;

using Microsoft.AspNetCore.Builder;
using Microsoft.AspNetCore.Routing;

public static class DependencyInjection
{
    public static IEndpointRouteBuilder MapPresentationEndpoints(this IEndpointRouteBuilder app)
    {
        // Orders
        Orders.CreateOrderEndpoint.Map(app);
        Orders.GetOrderByIdEndpoint.Map(app);
        
        // Payments
        // Payments.ProcessPaymentEndpoint.Map(app);
        
        return app;
    }
}
```

---

## 5. Infrastructure Layer

### 5.1 Repository Implementation

```csharp
// Infrastructure/Persistence/Orders/SqlOrderRepository.cs
namespace Infrastructure.Persistence.Orders;

using Microsoft.EntityFrameworkCore;
using Domain.Orders;
using Application.Features.Orders.Shared;

/// <summary>
/// Adapter - implements Port defined in Application Layer
/// </summary>
public sealed class SqlOrderRepository : IOrderRepository
{
    private readonly AppDbContext _context;
    
    public SqlOrderRepository(AppDbContext context)
    {
        _context = context;
    }
    
    public async Task<Order?> GetByIdAsync(OrderId id, CancellationToken ct = default)
    {
        return await _context.Orders
            .Include(o => o.Items)
            .FirstOrDefaultAsync(o => o.Id == id, ct);
    }
    
    public async Task AddAsync(Order order, CancellationToken ct = default)
    {
        await _context.Orders.AddAsync(order, ct);
    }
    
    public Task UpdateAsync(Order order, CancellationToken ct = default)
    {
        _context.Orders.Update(order);
        return Task.CompletedTask;
    }
}
```

### 5.2 Domain Event Dispatcher

```csharp
// Infrastructure/Messaging/DomainEventDispatcher.cs
namespace Infrastructure.Messaging;

using MediatR;
using Domain.Shared;
using Application.Shared.Interfaces;

public sealed class DomainEventDispatcher : IDomainEventDispatcher
{
    private readonly IPublisher _publisher;
    
    public DomainEventDispatcher(IPublisher publisher)
    {
        _publisher = publisher;
    }
    
    public async Task DispatchAsync(
        IEnumerable<IDomainEvent> domainEvents, 
        CancellationToken ct = default)
    {
        foreach (var domainEvent in domainEvents)
        {
            // Wrap domain event in MediatR notification
            var notification = new DomainEventNotification(domainEvent);
            await _publisher.Publish(notification, ct);
        }
    }
}

/// <summary>
/// Wrapper to publish domain events via MediatR
/// </summary>
public sealed record DomainEventNotification(IDomainEvent DomainEvent) : INotification;
```

### 5.3 Unit of Work

```csharp
// Infrastructure/Persistence/UnitOfWork.cs
namespace Infrastructure.Persistence;

using Application.Shared.Interfaces;

public sealed class UnitOfWork : IUnitOfWork
{
    private readonly AppDbContext _context;
    
    public UnitOfWork(AppDbContext context)
    {
        _context = context;
    }
    
    public async Task<int> SaveChangesAsync(CancellationToken ct = default)
    {
        return await _context.SaveChangesAsync(ct);
    }
}
```

---

## 6. Cross-Cutting Concerns (Behaviors)

### 6.1 Validation Behavior

```csharp
// Application/Shared/Behaviors/ValidationBehavior.cs
namespace Application.Shared.Behaviors;

using MediatR;
using FluentValidation;

/// <summary>
/// Runs FluentValidation validators BEFORE handler
/// Registered as open generic in DI
/// </summary>
public sealed class ValidationBehavior<TRequest, TResponse> 
    : IPipelineBehavior<TRequest, TResponse>
    where TRequest : IRequest<TResponse>
{
    private readonly IEnumerable<IValidator<TRequest>> _validators;
    
    public ValidationBehavior(IEnumerable<IValidator<TRequest>> validators)
    {
        _validators = validators;
    }
    
    public async Task<TResponse> Handle(
        TRequest request,
        RequestHandlerDelegate<TResponse> next,
        CancellationToken ct)
    {
        if (!_validators.Any())
            return await next();
        
        var context = new ValidationContext<TRequest>(request);
        
        var validationResults = await Task.WhenAll(
            _validators.Select(v => v.ValidateAsync(context, ct)));
        
        var failures = validationResults
            .SelectMany(r => r.Errors)
            .Where(f => f != null)
            .ToList();
        
        if (failures.Count != 0)
            throw new ValidationException(failures);
        
        return await next();
    }
}
```

### 6.2 Logging Behavior

```csharp
// Application/Shared/Behaviors/LoggingBehavior.cs
namespace Application.Shared.Behaviors;

using MediatR;
using Microsoft.Extensions.Logging;
using System.Diagnostics;

/// <summary>
/// Logs request handling with timing
/// </summary>
public sealed class LoggingBehavior<TRequest, TResponse> 
    : IPipelineBehavior<TRequest, TResponse>
    where TRequest : IRequest<TResponse>
{
    private readonly ILogger<LoggingBehavior<TRequest, TResponse>> _logger;
    
    public LoggingBehavior(ILogger<LoggingBehavior<TRequest, TResponse>> logger)
    {
        _logger = logger;
    }
    
    public async Task<TResponse> Handle(
        TRequest request,
        RequestHandlerDelegate<TResponse> next,
        CancellationToken ct)
    {
        var requestName = typeof(TRequest).Name;
        
        _logger.LogInformation("Handling {RequestName}", requestName);
        
        var stopwatch = Stopwatch.StartNew();
        
        try
        {
            var response = await next();
            
            stopwatch.Stop();
            _logger.LogInformation(
                "Handled {RequestName} in {ElapsedMs}ms",
                requestName, 
                stopwatch.ElapsedMilliseconds);
            
            return response;
        }
        catch (Exception ex)
        {
            stopwatch.Stop();
            _logger.LogError(
                ex,
                "Error handling {RequestName} after {ElapsedMs}ms",
                requestName,
                stopwatch.ElapsedMilliseconds);
            throw;
        }
    }
}
```

### 6.3 Transaction Behavior

```csharp
// Application/Shared/Behaviors/TransactionBehavior.cs
namespace Application.Shared.Behaviors;

using MediatR;
using Application.Shared.Interfaces;

/// <summary>
/// Wraps handler execution in a transaction
/// Only for commands (not queries)
/// </summary>
public sealed class TransactionBehavior<TRequest, TResponse> 
    : IPipelineBehavior<TRequest, TResponse>
    where TRequest : IRequest<TResponse>
{
    private readonly IUnitOfWork _unitOfWork;
    
    public TransactionBehavior(IUnitOfWork unitOfWork)
    {
        _unitOfWork = unitOfWork;
    }
    
    public async Task<TResponse> Handle(
        TRequest request,
        RequestHandlerDelegate<TResponse> next,
        CancellationToken ct)
    {
        // Skip for queries (naming convention)
        if (typeof(TRequest).Name.EndsWith("Query"))
            return await next();
        
        var response = await next();
        
        // Commit transaction
        await _unitOfWork.SaveChangesAsync(ct);
        
        return response;
    }
}
```

---

## 7. Dependency Injection Setup

### 7.1 Application Layer DI

```csharp
// Application/DependencyInjection.cs
namespace Application;

using Microsoft.Extensions.DependencyInjection;
using MediatR;
using FluentValidation;
using Application.Shared.Behaviors;

public static class DependencyInjection
{
    public static IServiceCollection AddApplication(this IServiceCollection services)
    {
        var assembly = typeof(DependencyInjection).Assembly;
        
        // MediatR
        services.AddMediatR(cfg => cfg.RegisterServicesFromAssembly(assembly));
        
        // FluentValidation
        services.AddValidatorsFromAssembly(assembly);
        
        // Pipeline Behaviors (order matters!)
        // 1. Logging (outermost)
        services.AddTransient(typeof(IPipelineBehavior<,>), typeof(LoggingBehavior<,>));
        // 2. Validation
        services.AddTransient(typeof(IPipelineBehavior<,>), typeof(ValidationBehavior<,>));
        // 3. Transaction (innermost before handler)
        services.AddTransient(typeof(IPipelineBehavior<,>), typeof(TransactionBehavior<,>));
        
        return services;
    }
}
```

### 7.2 Infrastructure Layer DI

```csharp
// Infrastructure/DependencyInjection.cs
namespace Infrastructure;

using Microsoft.Extensions.DependencyInjection;
using Microsoft.EntityFrameworkCore;
using Application.Features.Orders.Shared;
using Application.Shared.Interfaces;
using Infrastructure.Persistence;
using Infrastructure.Persistence.Orders;
using Infrastructure.Messaging;

public static class DependencyInjection
{
    public static IServiceCollection AddInfrastructure(
        this IServiceCollection services,
        string connectionString)
    {
        // DbContext
        services.AddDbContext<AppDbContext>(options =>
            options.UseSqlServer(connectionString));
        
        // Unit of Work
        services.AddScoped<IUnitOfWork, UnitOfWork>();
        
        // Repositories
        services.AddScoped<IOrderRepository, SqlOrderRepository>();
        services.AddScoped<ICustomerRepository, SqlCustomerRepository>();
        
        // Domain Event Dispatcher
        services.AddScoped<IDomainEventDispatcher, DomainEventDispatcher>();
        
        return services;
    }
}
```

### 7.3 Program.cs

```csharp
// Program.cs
using Application;
using Infrastructure;
using Presentation;

var builder = WebApplication.CreateBuilder(args);

// Add layers
builder.Services.AddApplication();
builder.Services.AddInfrastructure(
    builder.Configuration.GetConnectionString("Default")!);

var app = builder.Build();

// Map endpoints
app.MapPresentationEndpoints();

app.Run();
```

---

## 8. Testing a Vertical Slice

### 8.1 Unit Test for Handler

```csharp
// Tests/Application/Features/Orders/CreateOrderHandlerTests.cs
namespace Tests.Application.Features.Orders;

using Moq;
using Xunit;
using Domain.Orders;
using Domain.Shared;
using Application.Features.Orders.CreateOrder;
using Application.Features.Orders.Shared;
using Application.Shared.Interfaces;

public class CreateOrderHandlerTests
{
    private readonly Mock<IOrderRepository> _orderRepository;
    private readonly Mock<ICustomerRepository> _customerRepository;
    private readonly Mock<IUnitOfWork> _unitOfWork;
    private readonly Mock<IDomainEventDispatcher> _eventDispatcher;
    private readonly CreateOrderHandler _handler;
    
    public CreateOrderHandlerTests()
    {
        _orderRepository = new Mock<IOrderRepository>();
        _customerRepository = new Mock<ICustomerRepository>();
        _unitOfWork = new Mock<IUnitOfWork>();
        _eventDispatcher = new Mock<IDomainEventDispatcher>();
        
        _handler = new CreateOrderHandler(
            _orderRepository.Object,
            _customerRepository.Object,
            _unitOfWork.Object,
            _eventDispatcher.Object);
    }
    
    [Fact]
    public async Task Handle_ValidCommand_ReturnsSuccess()
    {
        // Arrange
        var customerId = Guid.NewGuid();
        var command = new CreateOrderCommand(
            customerId,
            new List<CreateOrderItemDto>
            {
                new(Guid.NewGuid(), 2, 100.00m, "USD")
            });
        
        _customerRepository
            .Setup(x => x.ExistsAsync(It.IsAny<CustomerId>(), It.IsAny<CancellationToken>()))
            .ReturnsAsync(true);
        
        // Act
        var result = await _handler.Handle(command, CancellationToken.None);
        
        // Assert
        Assert.True(result.IsSuccess);
        Assert.NotEqual(Guid.Empty, result.Value!.OrderId);
        Assert.Equal(200.00m, result.Value.TotalAmount);
        
        _orderRepository.Verify(x => x.AddAsync(It.IsAny<Order>(), It.IsAny<CancellationToken>()), Times.Once);
        _unitOfWork.Verify(x => x.SaveChangesAsync(It.IsAny<CancellationToken>()), Times.Once);
        _eventDispatcher.Verify(x => x.DispatchAsync(It.IsAny<IEnumerable<IDomainEvent>>(), It.IsAny<CancellationToken>()), Times.Once);
    }
    
    [Fact]
    public async Task Handle_CustomerNotFound_ReturnsFailure()
    {
        // Arrange
        var command = new CreateOrderCommand(
            Guid.NewGuid(),
            new List<CreateOrderItemDto>
            {
                new(Guid.NewGuid(), 1, 50.00m, "USD")
            });
        
        _customerRepository
            .Setup(x => x.ExistsAsync(It.IsAny<CustomerId>(), It.IsAny<CancellationToken>()))
            .ReturnsAsync(false);
        
        // Act
        var result = await _handler.Handle(command, CancellationToken.None);
        
        // Assert
        Assert.True(result.IsFailure);
        Assert.Equal("Customer not found", result.Error);
        
        _orderRepository.Verify(x => x.AddAsync(It.IsAny<Order>(), It.IsAny<CancellationToken>()), Times.Never);
    }
    
    [Fact]
    public async Task Handle_EmptyItems_ReturnsFailure()
    {
        // Arrange
        var command = new CreateOrderCommand(
            Guid.NewGuid(),
            new List<CreateOrderItemDto>());  // Empty items
        
        _customerRepository
            .Setup(x => x.ExistsAsync(It.IsAny<CustomerId>(), It.IsAny<CancellationToken>()))
            .ReturnsAsync(true);
        
        // Act
        var result = await _handler.Handle(command, CancellationToken.None);
        
        // Assert
        Assert.True(result.IsFailure);
        Assert.Contains("at least one item", result.Error);
    }
}
```

### 8.2 Domain Unit Test

```csharp
// Tests/Domain/Orders/OrderTests.cs
namespace Tests.Domain.Orders;

using Xunit;
using Domain.Orders;
using Domain.Shared;

public class OrderTests
{
    [Fact]
    public void Create_ValidData_ReturnsSuccess()
    {
        // Arrange
        var customerId = CustomerId.From(Guid.NewGuid());
        var items = new[]
        {
            (ProductId.From(Guid.NewGuid()), 2, Money.Create(100.00m)),
            (ProductId.From(Guid.NewGuid()), 1, Money.Create(50.00m))
        };
        
        // Act
        var result = Order.Create(customerId, items);
        
        // Assert
        Assert.True(result.IsSuccess);
        Assert.Equal(OrderStatus.Pending, result.Value!.Status);
        Assert.Equal(250.00m, result.Value.TotalAmount.Amount);
        Assert.Equal(2, result.Value.Items.Count);
        Assert.Single(result.Value.DomainEvents);
    }
    
    [Fact]
    public void Create_EmptyItems_ReturnsFailure()
    {
        // Arrange
        var customerId = CustomerId.From(Guid.NewGuid());
        var items = Array.Empty<(ProductId, int, Money)>();
        
        // Act
        var result = Order.Create(customerId, items);
        
        // Assert
        Assert.True(result.IsFailure);
        Assert.Contains("at least one item", result.Error);
    }
    
    [Fact]
    public void Confirm_PendingOrder_ReturnsSuccess()
    {
        // Arrange
        var order = CreateValidOrder();
        
        // Act
        var result = order.Confirm();
        
        // Assert
        Assert.True(result.IsSuccess);
        Assert.Equal(OrderStatus.Confirmed, result.Value!.Status);
    }
    
    [Fact]
    public void Cancel_DeliveredOrder_ReturnsFailure()
    {
        // Arrange
        var order = CreateValidOrder();
        order.Confirm();
        // Simulate delivered status (in real code, use Ship() and Deliver() methods)
        
        // Act & Assert - would fail on delivered order
        var cancelResult = order.Cancel("Changed my mind");
        Assert.True(cancelResult.IsSuccess); // Can cancel confirmed order
    }
    
    private static Order CreateValidOrder()
    {
        var result = Order.Create(
            CustomerId.From(Guid.NewGuid()),
            new[] { (ProductId.From(Guid.NewGuid()), 1, Money.Create(100.00m)) });
        return result.Value!;
    }
}
```

---

## 9. Query Example: GetOrderById

This section demonstrates a **Query** slice (read-only operation) to contrast with the **Command** example above.

### 9.1 Query vs Command Differences

| Aspect | Command (CreateOrder) | Query (GetOrderById) |
|--------|----------------------|---------------------|
| **Intent** | Change state | Read state |
| **Return** | Result + new entity ID | Result + DTO |
| **Side Effects** | Creates entity, raises events | None |
| **Validation** | Complex business rules | Simple parameter validation |
| **Transaction** | Yes (write to DB) | No (read-only) |
| **Caching** | Not applicable | Can be cached |
| **Domain Logic** | Calls Domain methods | Maps Domain to DTO |

### 9.2 GetOrderById Query

```csharp
// Application/Features/Orders/GetOrderById/GetOrderByIdQuery.cs
namespace Application.Features.Orders.GetOrderById;

using MediatR;
using Domain.Shared;

/// <summary>
/// Query DTO - read-only request
/// </summary>
public sealed record GetOrderByIdQuery(Guid OrderId)
    : IRequest<Result<GetOrderByIdResponse>>;
```

### 9.3 GetOrderById Response

```csharp
// Application/Features/Orders/GetOrderById/GetOrderByIdResponse.cs
namespace Application.Features.Orders.GetOrderById;

/// <summary>
/// Read model - optimized for display
/// Can be denormalized for performance
/// </summary>
public sealed record GetOrderByIdResponse(
    Guid Id,
    Guid CustomerId,
    string CustomerName,
    string Status,
    decimal TotalAmount,
    string Currency,
    DateTime CreatedAt,
    List<OrderItemResponse> Items
);

public sealed record OrderItemResponse(
    Guid ProductId,
    string ProductName,
    int Quantity,
    decimal UnitPrice,
    decimal Subtotal
);
```

### 9.4 GetOrderById Validator

```csharp
// Application/Features/Orders/GetOrderById/GetOrderByIdValidator.cs
namespace Application.Features.Orders.GetOrderById;

using FluentValidation;

/// <summary>
/// Simpler validation for queries - only format checks
/// No business rule validation (that's for commands)
/// </summary>
public sealed class GetOrderByIdValidator : AbstractValidator<GetOrderByIdQuery>
{
    public GetOrderByIdValidator()
    {
        RuleFor(x => x.OrderId)
            .NotEmpty()
            .WithMessage("Order ID is required");
    }
}
```

### 9.5 GetOrderById Handler

```csharp
// Application/Features/Orders/GetOrderById/GetOrderByIdHandler.cs
namespace Application.Features.Orders.GetOrderById;

using MediatR;
using Domain.Orders;
using Domain.Shared;
using Application.Features.Orders.Shared;

/// <summary>
/// Query Handler - read-only operation
///
/// Responsibilities:
/// 1. Load aggregate from repository
/// 2. Map Domain Entity to Response DTO
/// 3. Return response
///
/// No business logic, no state changes, no events
/// </summary>
public sealed class GetOrderByIdHandler
    : IRequestHandler<GetOrderByIdQuery, Result<GetOrderByIdResponse>>
{
    private readonly IOrderRepository _orderRepository;
    private readonly ICustomerRepository _customerRepository;

    public GetOrderByIdHandler(
        IOrderRepository orderRepository,
        ICustomerRepository customerRepository)
    {
        _orderRepository = orderRepository;
        _customerRepository = customerRepository;
    }

    public async Task<Result<GetOrderByIdResponse>> Handle(
        GetOrderByIdQuery request,
        CancellationToken ct)
    {
        // 1. Load aggregate (read-only)
        var order = await _orderRepository.GetByIdAsync(
            OrderId.From(request.OrderId),
            ct);

        if (order == null)
            return Result<GetOrderByIdResponse>.Failure("Order not found");

        // 2. Load related data (optional - for denormalization)
        var customer = await _customerRepository.GetByIdAsync(order.CustomerId, ct);

        // 3. Map Domain to Response DTO
        var response = MapToResponse(order, customer?.Name ?? "Unknown");

        // 4. Return response (no persistence, no events)
        return Result<GetOrderByIdResponse>.Success(response);
    }

    /// <summary>
    /// Maps Domain Entity to Response DTO
    /// This is where we denormalize for read performance
    /// </summary>
    private static GetOrderByIdResponse MapToResponse(Order order, string customerName)
    {
        return new GetOrderByIdResponse(
            order.Id.Value,
            order.CustomerId.Value,
            customerName,  // ← Denormalized from Customer aggregate
            order.Status.ToString(),
            order.TotalAmount.Amount,
            order.TotalAmount.Currency,
            order.CreatedAt,
            order.Items.Select(item => new OrderItemResponse(
                item.ProductId.Value,
                "Product Name",  // ← In real app, load from Product aggregate or read model
                item.Quantity,
                item.UnitPrice.Amount,
                item.Subtotal.Amount
            )).ToList()
        );
    }
}
```

### 9.6 GetOrderById Endpoint

```csharp
// Presentation/Orders/GetOrderByIdEndpoint.cs
namespace Presentation.Orders;

using Microsoft.AspNetCore.Http;
using Microsoft.AspNetCore.Routing;
using Microsoft.AspNetCore.Builder;
using MediatR;
using Application.Features.Orders.GetOrderById;

/// <summary>
/// Read endpoint - GET request
/// </summary>
public static class GetOrderByIdEndpoint
{
    public static void Map(IEndpointRouteBuilder app)
    {
        app.MapGet("/api/orders/{id:guid}", Handle)
            .WithTags("Orders")
            .WithName("GetOrderById")
            .WithDescription("Retrieves an order by ID")
            .Produces<GetOrderByIdResponse>(StatusCodes.Status200OK)
            .Produces(StatusCodes.Status404NotFound);
    }

    private static async Task<IResult> Handle(
        Guid id,  // ← From route parameter
        IMediator mediator,
        CancellationToken ct)
    {
        // 1. Create query from route parameter
        var query = new GetOrderByIdQuery(id);

        // 2. Send to Application Layer
        var result = await mediator.Send(query, ct);

        // 3. Map result to HTTP response
        return result.Match(
            success => Results.Ok(success),
            error => Results.NotFound(new { error })
        );
    }
}
```

### 9.7 Key Differences from Command

**No Validation of Business Rules:**
```csharp
// ✅ Query: Only validate format
RuleFor(x => x.OrderId).NotEmpty();

// ❌ Query: Don't validate business rules like "order status"
// Business rule validation is only for Commands
```

**No Domain Logic Call:**
```csharp
// ✅ Command: Calls Domain method
var orderResult = Order.Create(customerId, items);  // Business logic

// ✅ Query: Just loads and maps
var order = await _repository.GetByIdAsync(id);  // No business logic
var response = MapToResponse(order);
```

**No Transaction:**
```csharp
// Application/Shared/Behaviors/TransactionBehavior.cs
public async Task<TResponse> Handle(...)
{
    // Skip for queries (naming convention)
    if (typeof(TRequest).Name.EndsWith("Query"))
        return await next();  // ← Queries bypass transaction

    // Commands get transaction
    var response = await next();
    await _unitOfWork.SaveChangesAsync(ct);
    return response;
}
```

**No Domain Events:**
```csharp
// ✅ Command: Dispatches domain events
await _eventDispatcher.DispatchAsync(order.DomainEvents, ct);

// ✅ Query: No events (read-only)
// Just returns data
```

### 9.8 Query Testing

```csharp
// Tests/Application/Features/Orders/GetOrderByIdHandlerTests.cs
namespace Tests.Application.Features.Orders;

using Moq;
using Xunit;
using Domain.Orders;
using Domain.Shared;
using Application.Features.Orders.GetOrderById;
using Application.Features.Orders.Shared;

public class GetOrderByIdHandlerTests
{
    private readonly Mock<IOrderRepository> _orderRepository;
    private readonly Mock<ICustomerRepository> _customerRepository;
    private readonly GetOrderByIdHandler _handler;

    public GetOrderByIdHandlerTests()
    {
        _orderRepository = new Mock<IOrderRepository>();
        _customerRepository = new Mock<ICustomerRepository>();

        _handler = new GetOrderByIdHandler(
            _orderRepository.Object,
            _customerRepository.Object);
    }

    [Fact]
    public async Task Handle_ExistingOrder_ReturnsSuccess()
    {
        // Arrange
        var orderId = OrderId.New();
        var order = CreateTestOrder(orderId);

        _orderRepository
            .Setup(x => x.GetByIdAsync(orderId, It.IsAny<CancellationToken>()))
            .ReturnsAsync(order);

        var query = new GetOrderByIdQuery(orderId.Value);

        // Act
        var result = await _handler.Handle(query, CancellationToken.None);

        // Assert
        Assert.True(result.IsSuccess);
        Assert.Equal(orderId.Value, result.Value!.Id);
        Assert.Equal(order.Status.ToString(), result.Value.Status);
        Assert.Equal(order.TotalAmount.Amount, result.Value.TotalAmount);
    }

    [Fact]
    public async Task Handle_NonExistentOrder_ReturnsFailure()
    {
        // Arrange
        var orderId = OrderId.New();

        _orderRepository
            .Setup(x => x.GetByIdAsync(orderId, It.IsAny<CancellationToken>()))
            .ReturnsAsync((Order?)null);  // Order not found

        var query = new GetOrderByIdQuery(orderId.Value);

        // Act
        var result = await _handler.Handle(query, CancellationToken.None);

        // Assert
        Assert.True(result.IsFailure);
        Assert.Equal("Order not found", result.Error);
    }

    private static Order CreateTestOrder(OrderId orderId)
    {
        var result = Order.Create(
            CustomerId.From(Guid.NewGuid()),
            new[] { (ProductId.From(Guid.NewGuid()), 1, Money.Create(100.00m)) });

        // Use reflection to set Id for testing (in real code, Id is set in Create method)
        var order = result.Value!;
        typeof(Order)
            .GetProperty(nameof(Order.Id))!
            .SetValue(order, orderId);

        return order;
    }
}
```

### 9.9 CQRS Pattern in Action

This example demonstrates the **Command Query Responsibility Segregation** pattern:

```
┌─────────────────────────────────────────────────────────┐
│                    WRITE SIDE (Commands)                │
├─────────────────────────────────────────────────────────┤
│ CreateOrderCommand                                      │
│ ├── Complex validation (business rules)                │
│ ├── Calls Domain methods (Order.Create)                │
│ ├── Writes to database                                 │
│ ├── Raises Domain Events                               │
│ └── Returns Result<OrderId>                            │
└─────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────┐
│                    READ SIDE (Queries)                  │
├─────────────────────────────────────────────────────────┤
│ GetOrderByIdQuery                                       │
│ ├── Simple validation (format only)                    │
│ ├── Loads from repository (read-only)                  │
│ ├── Maps to denormalized DTO                           │
│ ├── No Domain Events                                   │
│ └── Returns Result<GetOrderByIdResponse>               │
└─────────────────────────────────────────────────────────┘
```

**Benefits of Separation:**
- **Optimized models**: Write model normalized, read model denormalized
- **Independent scaling**: Scale reads and writes separately
- **Clarity**: Intent is explicit (changing vs querying)
- **Performance**: Queries can bypass validation, transactions, event dispatching

---

## Summary

This example demonstrates:

1. **Proper Layer Separation**: Domain, Application, Presentation, and Infrastructure in separate projects/folders
2. **Vertical Slices in Application Layer**: Features organized by use case, not technical concern
3. **Rich Domain Model**: Business logic in Aggregate Root with invariants
4. **API DTO Mapping**: Presentation layer maps API DTOs to Application Commands
5. **Ports and Adapters**: Interfaces in Application, implementations in Infrastructure
6. **Pipeline Behaviors**: Cross-cutting concerns (logging, validation, transactions) handled uniformly
7. **Testability**: Each layer can be tested in isolation
8. **CQRS Pattern**: Clear separation between Commands (write) and Queries (read) with different responsibilities and optimizations
