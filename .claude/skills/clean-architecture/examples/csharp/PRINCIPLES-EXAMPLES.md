# Clean Architecture - SOLID & DDD Principles Examples (C#)

Examples demonstrating SOLID principles and Domain-Driven Design tactical patterns in C#/.NET.

> These examples complement [PRINCIPLES.md](../../PRINCIPLES.md)

## SOLID Principles in Practice

### Single Responsibility Principle (SRP)

Each class has one reason to change:

```csharp
// ✅ GOOD - Single responsibility: Handle one command
public class CreateOrderCommandHandler : IRequestHandler<CreateOrderCommand, OrderId>
{
    private readonly IOrderRepository _orderRepository;
    private readonly IProductRepository _productRepository;
    private readonly IUnitOfWork _unitOfWork;

    public async Task<OrderId> Handle(CreateOrderCommand request, CancellationToken cancellationToken)
    {
        var customerId = CustomerId.From(request.CustomerId);
        var shippingAddress = ShippingAddress.From(request.ShippingAddress);
        var order = Order.Create(customerId, shippingAddress);

        foreach (var itemDto in request.Items)
        {
            var product = await _productRepository.GetByIdAsync(
                ProductId.From(itemDto.ProductId),
                cancellationToken);

            if (product == null)
                throw new NotFoundException($"Product {itemDto.ProductId} not found");

            order.AddItem(product, itemDto.Quantity);
        }

        order.PlaceOrder();

        await _orderRepository.AddAsync(order, cancellationToken);
        await _unitOfWork.SaveChangesAsync(cancellationToken);

        return order.Id;
    }
}

// ✅ GOOD - Single responsibility: Cancel one order
public class CancelOrderCommandHandler : IRequestHandler<CancelOrderCommand>
{
    private readonly IOrderRepository _orderRepository;
    private readonly IUnitOfWork _unitOfWork;

    public async Task Handle(CancelOrderCommand request, CancellationToken cancellationToken)
    {
        var order = await _orderRepository.GetByIdAsync(
            OrderId.From(request.OrderId),
            cancellationToken);

        if (order == null)
            throw new NotFoundException($"Order {request.OrderId} not found");

        order.Cancel(new CancellationReason(request.Reason));

        await _orderRepository.UpdateAsync(order, cancellationToken);
        await _unitOfWork.SaveChangesAsync(cancellationToken);
    }
}
```

### Open/Closed Principle (OCP)

Open for extension, closed for modification:

```csharp
// ✅ GOOD - Can add new specifications without modifying existing code
public interface ISpecification<T>
{
    bool IsSatisfiedBy(T entity);
    Expression<Func<T, bool>> ToExpression();
}

public class OrderIsPlacedSpecification : ISpecification<Order>
{
    public bool IsSatisfiedBy(Order order)
    {
        return order.Status == OrderStatus.Placed;
    }

    public Expression<Func<Order, bool>> ToExpression()
    {
        return order => order.Status == OrderStatus.Placed;
    }
}

public class OrderTotalExceedsSpecification : ISpecification<Order>
{
    private readonly Money _threshold;

    public OrderTotalExceedsSpecification(Money threshold)
    {
        _threshold = threshold;
    }

    public bool IsSatisfiedBy(Order order)
    {
        return order.Total.IsGreaterThan(_threshold);
    }

    public Expression<Func<Order, bool>> ToExpression()
    {
        var amount = _threshold.Amount;
        return order => order.Total.Amount > amount;
    }
}

// Can compose specifications without modifying existing ones
public class AndSpecification<T> : ISpecification<T>
{
    private readonly ISpecification<T> _left;
    private readonly ISpecification<T> _right;

    public AndSpecification(ISpecification<T> left, ISpecification<T> right)
    {
        _left = left;
        _right = right;
    }

    public bool IsSatisfiedBy(T entity)
    {
        return _left.IsSatisfiedBy(entity) && _right.IsSatisfiedBy(entity);
    }

    public Expression<Func<T, bool>> ToExpression()
    {
        var leftExpr = _left.ToExpression();
        var rightExpr = _right.ToExpression();

        var parameter = Expression.Parameter(typeof(T));
        var combined = Expression.AndAlso(
            Expression.Invoke(leftExpr, parameter),
            Expression.Invoke(rightExpr, parameter));

        return Expression.Lambda<Func<T, bool>>(combined, parameter);
    }
}
```

### Liskov Substitution Principle (LSP)

Subtypes must be substitutable for base types:

```csharp
// ✅ GOOD - All implementations honor the contract
public interface IOrderRepository
{
    Task<Order?> GetByIdAsync(OrderId id, CancellationToken cancellationToken = default);
    Task AddAsync(Order order, CancellationToken cancellationToken = default);
    Task UpdateAsync(Order order, CancellationToken cancellationToken = default);
}

// SQL implementation
public class SqlOrderRepository : IOrderRepository
{
    private readonly ApplicationDbContext _context;

    public async Task<Order?> GetByIdAsync(OrderId id, CancellationToken cancellationToken = default)
    {
        return await _context.Orders
            .Include(o => o.Items)
            .FirstOrDefaultAsync(o => o.Id == id, cancellationToken);
    }

    public async Task AddAsync(Order order, CancellationToken cancellationToken = default)
    {
        await _context.Orders.AddAsync(order, cancellationToken);
    }

    public Task UpdateAsync(Order order, CancellationToken cancellationToken = default)
    {
        _context.Orders.Update(order);
        return Task.CompletedTask;
    }
}

// In-memory implementation (for testing) - honors same contract
public class InMemoryOrderRepository : IOrderRepository
{
    private readonly Dictionary<OrderId, Order> _orders = new();

    public Task<Order?> GetByIdAsync(OrderId id, CancellationToken cancellationToken = default)
    {
        _orders.TryGetValue(id, out var order);
        return Task.FromResult(order);
    }

    public Task AddAsync(Order order, CancellationToken cancellationToken = default)
    {
        _orders[order.Id] = order;
        return Task.CompletedTask;
    }

    public Task UpdateAsync(Order order, CancellationToken cancellationToken = default)
    {
        _orders[order.Id] = order;
        return Task.CompletedTask;
    }
}
```

### Interface Segregation Principle (ISP)

Many small interfaces better than one large:

```csharp
// ✅ GOOD - Segregated interfaces, clients only depend on what they need

// Order read operations
public interface IOrderQueryRepository
{
    Task<Order?> GetByIdAsync(OrderId id, CancellationToken ct = default);
    Task<IEnumerable<Order>> GetByCustomerAsync(CustomerId customerId, CancellationToken ct = default);
}

// Order write operations
public interface IOrderCommandRepository
{
    Task AddAsync(Order order, CancellationToken ct = default);
    Task UpdateAsync(Order order, CancellationToken ct = default);
    Task DeleteAsync(OrderId id, CancellationToken ct = default);
}

// Query handlers only need read interface
public class GetOrderByIdQueryHandler : IRequestHandler<GetOrderByIdQuery, OrderDto>
{
    private readonly IOrderQueryRepository _repository;

    // Only depends on query operations
    public GetOrderByIdQueryHandler(IOrderQueryRepository repository)
    {
        _repository = repository;
    }
}

// Command handlers need write interface
public class CreateOrderCommandHandler : IRequestHandler<CreateOrderCommand, OrderId>
{
    private readonly IOrderCommandRepository _repository;

    // Only depends on command operations
    public CreateOrderCommandHandler(IOrderCommandRepository repository)
    {
        _repository = repository;
    }
}
```

### Dependency Inversion Principle (DIP)

Depend on abstractions, not concretions:

```csharp
// ✅ GOOD - High-level policy (Application) defines the interface
// Domain/Repositories/IOrderRepository.cs
namespace YourApp.Domain.Repositories;

public interface IOrderRepository
{
    Task<Order?> GetByIdAsync(OrderId id, CancellationToken cancellationToken = default);
    Task AddAsync(Order order, CancellationToken cancellationToken = default);
}

// ✅ GOOD - Low-level detail (Infrastructure) implements the interface
// Infrastructure/Persistence/Repositories/OrderRepository.cs
namespace YourApp.Infrastructure.Persistence.Repositories;

public class OrderRepository : IOrderRepository
{
    private readonly ApplicationDbContext _context;

    public OrderRepository(ApplicationDbContext context)
    {
        _context = context;
    }

    public async Task<Order?> GetByIdAsync(OrderId id, CancellationToken cancellationToken = default)
    {
        return await _context.Orders
            .Include(o => o.Items)
            .FirstOrDefaultAsync(o => o.Id == id, cancellationToken);
    }

    public async Task AddAsync(Order order, CancellationToken cancellationToken = default)
    {
        await _context.Orders.AddAsync(order, cancellationToken);
    }
}

// ✅ GOOD - Dependency injection wires them up
// API/Program.cs
builder.Services.AddScoped<IOrderRepository, OrderRepository>();
```

---

## Domain-Driven Design (DDD) Tactical Patterns

### Entities

Objects with unique identity. See [Complete Feature Example](CHECKLIST-EXAMPLES.md#3-domain-layer) for full implementation.

**Key characteristics:**
```csharp
public class Order
{
    // ✅ Has unique identity
    public OrderId Id { get; private set; }

    // ✅ Private setters protect invariants
    public OrderStatus Status { get; private set; }
    public Money Total { get; private set; }

    // ✅ Encapsulated collections
    private readonly List<OrderItem> _items = new();
    public IReadOnlyList<OrderItem> Items => _items.AsReadOnly();

    // ✅ Private constructor enforces factory usage
    private Order() { }

    // ✅ Factory method with business rules
    public static Order Create(CustomerId customerId, ShippingAddress address)
    {
        // Validation and initialization
        var order = new Order { /* ... */ };
        order._domainEvents.Add(new OrderCreatedEvent(/*...*/));
        return order;
    }

    // ✅ Business logic in methods, not setters
    public void AddItem(Product product, int quantity)
    {
        // Enforce invariants
        if (Status != OrderStatus.Draft)
            throw new DomainException("Cannot add items to non-draft order");

        // Business logic
        _items.Add(OrderItem.Create(product, quantity));
        RecalculateTotal(); // Maintain invariant: Total = sum(items)
    }

    // ✅ Raises domain events for significant changes
    public void PlaceOrder()
    {
        if (!_items.Any())
            throw new DomainException("Cannot place empty order");

        Status = OrderStatus.Placed;
        _domainEvents.Add(new OrderPlacedEvent(Id, CustomerId, Total));
    }
}
```

### Value Objects

Immutable objects defined by their attributes:

```csharp
// Domain/ValueObjects/Money.cs
namespace YourApp.Domain.ValueObjects;

public class Money : IEquatable<Money>, IComparable<Money>
{
    public decimal Amount { get; }
    public Currency Currency { get; }

    private Money() { } // For EF Core

    public Money(decimal amount, Currency currency)
    {
        if (amount < 0)
            throw new ArgumentException("Amount cannot be negative", nameof(amount));

        // Use MidpointRounding.AwayFromZero for financial calculations to avoid precision loss
        Amount = decimal.Round(amount, 2, MidpointRounding.AwayFromZero);
        Currency = currency ?? throw new ArgumentNullException(nameof(currency));
    }

    public static Money Zero(Currency currency) => new Money(0, currency);

    // Immutable operations
    public Money Add(Money other)
    {
        if (other == null)
            throw new ArgumentNullException(nameof(other));
        if (Currency != other.Currency)
            throw new InvalidOperationException($"Cannot add {other.Currency} to {Currency}");

        return new Money(Amount + other.Amount, Currency);
    }

    public Money Subtract(Money other)
    {
        if (other == null)
            throw new ArgumentNullException(nameof(other));
        if (Currency != other.Currency)
            throw new InvalidOperationException($"Cannot subtract {other.Currency} from {Currency}");

        return new Money(Amount - other.Amount, Currency);
    }

    public Money Multiply(decimal multiplier)
    {
        return new Money(Amount * multiplier, Currency);
    }

    public bool IsGreaterThan(Money other)
    {
        if (other == null)
            throw new ArgumentNullException(nameof(other));
        if (Currency != other.Currency)
            throw new InvalidOperationException("Cannot compare different currencies");

        return Amount > other.Amount;
    }

    // Equality based on all properties
    public bool Equals(Money other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true; // Performance optimization
        return Amount == other.Amount && Currency == other.Currency;
    }

    public override bool Equals(object obj) => obj is Money money && Equals(money);

    public override int GetHashCode() => HashCode.Combine(Amount, Currency);

    public int CompareTo(Money other)
    {
        if (other == null) return 1;
        if (Currency != other.Currency)
            throw new InvalidOperationException("Cannot compare different currencies");

        return Amount.CompareTo(other.Amount);
    }

    public override string ToString() => $"{Amount:N2} {Currency}";

    public static bool operator ==(Money left, Money right) =>
        left?.Equals(right) ?? right is null;

    public static bool operator !=(Money left, Money right) => !(left == right);

    public static bool operator >(Money left, Money right) =>
        left?.CompareTo(right) > 0;

    public static bool operator <(Money left, Money right) =>
        left?.CompareTo(right) < 0;
}

// Domain/ValueObjects/Email.cs
public class Email : IEquatable<Email>
{
    public string Address { get; }

    private Email() { } // For EF Core

    public Email(string address)
    {
        if (string.IsNullOrWhiteSpace(address))
            throw new ArgumentException("Email cannot be empty", nameof(address));

        if (!IsValid(address))
            throw new ArgumentException("Invalid email format", nameof(address));

        Address = address.ToLowerInvariant();
    }

    // Cached compiled regex for better performance (thread-safe)
    private static readonly System.Text.RegularExpressions.Regex EmailRegex = new(
        @"^[^@\s]+@[^@\s]+\.[^@\s]+$",
        System.Text.RegularExpressions.RegexOptions.IgnoreCase | 
        System.Text.RegularExpressions.RegexOptions.Compiled);

    private static bool IsValid(string email)
    {
        if (string.IsNullOrWhiteSpace(email))
            return false;

        try
        {
            // Basic validation using MailAddress
            var addr = new System.Net.Mail.MailAddress(email);

            // Additional validation: check format more strictly (RFC 5322 simplified)
            return addr.Address == email && EmailRegex.IsMatch(email);
        }
        catch
        {
            return false;
        }
    }

    public bool Equals(Email other)
    {
        if (other is null) return false;
        return Address == other.Address;
    }

    public override bool Equals(object obj) => obj is Email email && Equals(email);
    public override int GetHashCode() => Address.GetHashCode();
    public override string ToString() => Address;

    public static implicit operator string(Email email) => email.Address;
}
```

### Strongly-Typed IDs

Type-safe identifiers using records:

```csharp
// Domain/ValueObjects/OrderId.cs
namespace YourApp.Domain.ValueObjects;

public record OrderId(Guid Value)
{
    public static OrderId NewId() => new(Guid.NewGuid());
    public static OrderId From(Guid value) => new(value);
    public static OrderId From(string value) => new(Guid.Parse(value));

    public override string ToString() => Value.ToString();
}

public record CustomerId(Guid Value)
{
    public static CustomerId NewId() => new(Guid.NewGuid());
    public static CustomerId From(Guid value) => new(value);
    public override string ToString() => Value.ToString();
}

public record ProductId(Guid Value)
{
    public static ProductId NewId() => new(Guid.NewGuid());
    public static ProductId From(Guid value) => new(value);
    public override string ToString() => Value.ToString();
}
```

### Domain Events

Immutable, past-tense events:

```csharp
// Domain/DomainEvents/IDomainEvent.cs
namespace YourApp.Domain.DomainEvents;

public interface IDomainEvent
{
    DateTime OccurredAt { get; }
}

// Domain/DomainEvents/OrderPlacedEvent.cs
public class OrderPlacedEvent : IDomainEvent
{
    public OrderId OrderId { get; }
    public CustomerId CustomerId { get; }
    public Money Total { get; }
    public DateTime OccurredAt { get; }

    public OrderPlacedEvent(OrderId orderId, CustomerId customerId, Money total)
    {
        OrderId = orderId;
        CustomerId = customerId;
        Total = total;
        OccurredAt = DateTime.UtcNow;
    }
}

public class OrderCancelledEvent : IDomainEvent
{
    public OrderId OrderId { get; }
    public CancellationReason Reason { get; }
    public DateTime OccurredAt { get; }

    public OrderCancelledEvent(OrderId orderId, CancellationReason reason)
    {
        OrderId = orderId;
        Reason = reason;
        OccurredAt = DateTime.UtcNow;
    }
}
```

### Repositories

Collection-like interface for aggregates:

```csharp
// Domain/Repositories/IOrderRepository.cs
namespace YourApp.Domain.Repositories;

public interface IOrderRepository
{
    Task<Order?> GetByIdAsync(OrderId id, CancellationToken cancellationToken = default);
    Task<IEnumerable<Order>> GetByCustomerAsync(CustomerId customerId, CancellationToken cancellationToken = default);
    Task<Order?> GetByOrderNumberAsync(string orderNumber, CancellationToken cancellationToken = default);
    Task AddAsync(Order order, CancellationToken cancellationToken = default);
    Task UpdateAsync(Order order, CancellationToken cancellationToken = default);
    Task DeleteAsync(OrderId id, CancellationToken cancellationToken = default);
    Task<bool> ExistsAsync(OrderId id, CancellationToken cancellationToken = default);
}
```

### Specifications

Composable business rules:

```csharp
// Domain/Specifications/ISpecification.cs
namespace YourApp.Domain.Specifications;

public interface ISpecification<T>
{
    bool IsSatisfiedBy(T entity);
    Expression<Func<T, bool>> ToExpression();
}

// Domain/Specifications/OrderSpecifications.cs
public class OrderIsPlacedSpecification : ISpecification<Order>
{
    public bool IsSatisfiedBy(Order order)
    {
        return order.Status == OrderStatus.Placed;
    }

    public Expression<Func<Order, bool>> ToExpression()
    {
        return order => order.Status == OrderStatus.Placed;
    }
}

public class OrderTotalExceedsSpecification : ISpecification<Order>
{
    private readonly Money _threshold;

    public OrderTotalExceedsSpecification(Money threshold)
    {
        _threshold = threshold;
    }

    public bool IsSatisfiedBy(Order order)
    {
        return order.Total.IsGreaterThan(_threshold);
    }

    public Expression<Func<Order, bool>> ToExpression()
    {
        var amount = _threshold.Amount;
        return order => order.Total.Amount > amount;
    }
}

// Example: Composing specifications with AND/OR
public class CompositeSpecificationExample
{
    public void Example()
    {
        // Create individual specifications
        var isPlaced = new OrderIsPlacedSpecification();
        var exceedsThreshold = new OrderTotalExceedsSpecification(new Money(100, Currency.USD));

        // Compose: placed AND exceeds threshold
        var highValuePlacedOrders = new AndSpecification<Order>(isPlaced, exceedsThreshold);

        // Use in repository query
        var orders = await _orderRepository.FindAsync(
            highValuePlacedOrders.ToExpression(),
            cancellationToken);

        // Or use in-memory
        var matchingOrders = allOrders.Where(o => highValuePlacedOrders.IsSatisfiedBy(o));
    }
}
```

### Ubiquitous Language

Code using same terms as domain experts:

```csharp
// ✅ GOOD - Domain language is clear
public class Order
{
    public void PlaceOrder() { /* ... */ }
    public void Cancel(CancellationReason reason) { /* ... */ }
    public void Ship(TrackingNumber trackingNumber) { /* ... */ }
}

public enum OrderStatus
{
    Draft,
    Placed,
    Paid,
    Shipped,
    Delivered,
    Cancelled
}

// ❌ BAD - Technical jargon, not domain language
public class Order
{
    public void SetStatus(int statusCode) { /* ... */ }
    public void UpdateRecord() { /* ... */ }
    public void Persist() { /* ... */ }
}
```

---

## Key Takeaways

1. **SOLID makes code maintainable** - Each principle addresses specific design problems
2. **DDD tactical patterns model the domain** - Entities, Value Objects, Events reflect business concepts
3. **Value Objects are immutable** - Validation in constructor, equality by value
4. **Entities have identity** - Can change over time but maintain same ID
5. **Domain Events decouple side effects** - Main flow stays clean
6. **Specifications are composable** - Reusable business rules
7. **Repositories abstract persistence** - Domain doesn't know about database
8. **Ubiquitous Language bridges gaps** - Code speaks business language
