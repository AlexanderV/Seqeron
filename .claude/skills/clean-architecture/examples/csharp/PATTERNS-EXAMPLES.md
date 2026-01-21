# Clean Architecture - Design Patterns Examples (C#)

Examples demonstrating design patterns commonly used in Clean Architecture applications with C#/.NET.

> These examples complement [PATTERNS.md](../../PATTERNS.md)

## Creational Patterns

### Factory Pattern

Encapsulate complex object creation:

```csharp
// Domain/Entities/Order.cs - Static Factory Method
public class Order
{
    // Private constructor prevents direct instantiation
    private Order() { }

    // Static factory method with clear intent
    public static Order Create(CustomerId customerId, ShippingAddress shippingAddress)
    {
        // Validate before creating
        if (customerId == null)
            throw new ArgumentNullException(nameof(customerId));
        if (shippingAddress == null)
            throw new ArgumentNullException(nameof(shippingAddress));

        var order = new Order
        {
            Id = OrderId.NewId(),
            CustomerId = customerId,
            Status = OrderStatus.Draft,
            ShippingAddress = shippingAddress,
            Total = Money.Zero(Currency.USD),
            CreatedAt = DateTime.UtcNow
        };

        // Raise domain event after successful creation
        order._domainEvents.Add(new OrderCreatedEvent(order.Id, order.CustomerId));
        return order;
    }

    // Named factory methods for different scenarios
    public static Order CreateGiftOrder(CustomerId customerId, ShippingAddress shippingAddress, string giftMessage)
    {
        var order = Create(customerId, shippingAddress);
        order.MarkAsGift(giftMessage);
        return order;
    }
}

// Domain/ValueObjects/Money.cs - Factory Methods
public class Money
{
    public static Money Zero(Currency currency) => new Money(0, currency);
    public static Money FromDecimal(decimal amount, Currency currency) => new Money(amount, currency);
}
```

### Builder Pattern

Construct complex objects step-by-step:

```csharp
// Domain/Builders/OrderBuilder.cs
namespace YourApp.Domain.Builders;

public class OrderBuilder
{
    private CustomerId _customerId;
    private ShippingAddress _shippingAddress;
    private readonly List<(Product product, int quantity)> _items = new();
    private string _giftMessage;

    public OrderBuilder ForCustomer(CustomerId customerId)
    {
        _customerId = customerId;
        return this;
    }

    public OrderBuilder ShipTo(ShippingAddress address)
    {
        _shippingAddress = address;
        return this;
    }

    public OrderBuilder AddItem(Product product, int quantity)
    {
        _items.Add((product, quantity));
        return this;
    }

    public OrderBuilder AsGift(string message)
    {
        _giftMessage = message;
        return this;
    }

    public Order Build()
    {
        if (_customerId == null)
            throw new InvalidOperationException("Customer is required");
        if (_shippingAddress == null)
            throw new InvalidOperationException("Shipping address is required");

        var order = Order.Create(_customerId, _shippingAddress);

        foreach (var (product, quantity) in _items)
        {
            order.AddItem(product, quantity);
        }

        if (!string.IsNullOrEmpty(_giftMessage))
        {
            order.MarkAsGift(_giftMessage);
        }

        return order;
    }
}

// Usage:
var order = new OrderBuilder()
    .ForCustomer(customerId)
    .ShipTo(address)
    .AddItem(product1, 2)
    .AddItem(product2, 1)
    .AsGift("Happy Birthday!")
    .Build();
```

### Dependency Injection Pattern

Invert dependencies throughout the application:

```csharp
// Domain defines interface
namespace YourApp.Domain.Repositories;

public interface IOrderRepository
{
    Task<Order?> GetByIdAsync(OrderId id, CancellationToken ct = default);
    Task AddAsync(Order order, CancellationToken ct = default);
}

// Infrastructure implements interface
namespace YourApp.Infrastructure.Persistence.Repositories;

public class OrderRepository : IOrderRepository
{
    private readonly ApplicationDbContext _context;

    public OrderRepository(ApplicationDbContext context)
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
}

// Application uses interface (not implementation)
namespace YourApp.Application.Commands;

public class CreateOrderCommandHandler : IRequestHandler<CreateOrderCommand, OrderId>
{
    private readonly IOrderRepository _orderRepository; // Depends on abstraction

    public CreateOrderCommandHandler(IOrderRepository orderRepository)
    {
        _orderRepository = orderRepository;
    }

    public async Task<OrderId> Handle(CreateOrderCommand request, CancellationToken ct)
    {
        var order = Order.Create(/*...*/);
        await _orderRepository.AddAsync(order, ct);
        return order.Id;
    }
}

// Composition root wires everything up
// API/Program.cs
builder.Services.AddScoped<IOrderRepository, OrderRepository>();
```

---

## Structural Patterns

### Repository Pattern

Abstract data access with collection-like interface:

```csharp
// Domain/Repositories/IOrderRepository.cs - Interface in Domain
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

// Infrastructure/Persistence/Repositories/OrderRepository.cs - Implementation in Infrastructure
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
                .ThenInclude(i => i.Product)
            .FirstOrDefaultAsync(o => o.Id == id, cancellationToken);
    }

    public async Task<IEnumerable<Order>> GetByCustomerAsync(
        CustomerId customerId,
        CancellationToken cancellationToken = default)
    {
        return await _context.Orders
            .Where(o => o.CustomerId == customerId)
            .Include(o => o.Items)
            .OrderByDescending(o => o.CreatedAt)
            .ToListAsync(cancellationToken);
    }

    public async Task AddAsync(Order order, CancellationToken cancellationToken = default)
    {
        await _context.Orders.AddAsync(order, cancellationToken);
    }

    public async Task UpdateAsync(Order order, CancellationToken cancellationToken = default)
    {
        var exists = await _context.Orders.AnyAsync(o => o.Id == order.Id, cancellationToken);
        if (!exists)
            throw new NotFoundException($"Order {order.Id} not found");

        _context.Orders.Update(order);
    }

    public async Task DeleteAsync(OrderId id, CancellationToken cancellationToken = default)
    {
        var order = await GetByIdAsync(id, cancellationToken);
        if (order != null)
        {
            _context.Orders.Remove(order);
        }
    }

    public async Task<bool> ExistsAsync(OrderId id, CancellationToken cancellationToken = default)
    {
        return await _context.Orders.AnyAsync(o => o.Id == id, cancellationToken);
    }
}
```

### Unit of Work Pattern

Coordinate changes across multiple repositories:

```csharp
// Application/Interfaces/IUnitOfWork.cs - Interface in Application
namespace YourApp.Application.Interfaces;

public interface IUnitOfWork
{
    Task<int> SaveChangesAsync(CancellationToken cancellationToken = default);
}

// Infrastructure/Persistence/ApplicationDbContext.cs - Implementation in Infrastructure
namespace YourApp.Infrastructure.Persistence;

// Domain/Entities/IAggregateRoot.cs
public interface IAggregateRoot
{
    IReadOnlyList<IDomainEvent> DomainEvents { get; }
    void ClearDomainEvents();
}

public class ApplicationDbContext : DbContext, IUnitOfWork
{
    private readonly IDomainEventDispatcher _domainEventDispatcher;

    public ApplicationDbContext(
        DbContextOptions<ApplicationDbContext> options,
        IDomainEventDispatcher domainEventDispatcher) : base(options)
    {
        _domainEventDispatcher = domainEventDispatcher;
    }

    public DbSet<Order> Orders { get; set; }
    public DbSet<Customer> Customers { get; set; }
    public DbSet<Product> Products { get; set; }

    public override async Task<int> SaveChangesAsync(CancellationToken cancellationToken = default)
    {
        // Collect domain events before saving
        var domainEvents = ChangeTracker.Entries<IAggregateRoot>()
            .Select(x => x.Entity)
            .SelectMany(x => x.DomainEvents)
            .ToList();

        // Save changes
        var result = await base.SaveChangesAsync(cancellationToken);

        // Dispatch events after successful save
        foreach (var domainEvent in domainEvents)
        {
            await _domainEventDispatcher.DispatchAsync(domainEvent, cancellationToken);
        }

        // Clear events after dispatch
        foreach (var entry in ChangeTracker.Entries<IAggregateRoot>())
        {
            entry.Entity.ClearDomainEvents();
        }

        return result;
    }
}

// Usage in command handler
public class CreateOrderCommandHandler : IRequestHandler<CreateOrderCommand, OrderId>
{
    private readonly IOrderRepository _orderRepository;
    private readonly IProductRepository _productRepository;
    private readonly IUnitOfWork _unitOfWork;

    public async Task<OrderId> Handle(CreateOrderCommand request, CancellationToken cancellationToken)
    {
        var order = Order.Create(/*...*/);
        await _orderRepository.AddAsync(order, cancellationToken);

        // Commit all changes in one transaction
        await _unitOfWork.SaveChangesAsync(cancellationToken);

        return order.Id;
    }
}
```

### Adapter Pattern (Hexagonal Architecture)

Adapt external services to application interfaces:

```csharp
// Application/Interfaces/IEmailService.cs - Port (interface) defined in Application
namespace YourApp.Application.Interfaces;

public interface IEmailService
{
    Task SendOrderConfirmationAsync(OrderId orderId, CancellationToken cancellationToken = default);
    Task SendShippingNotificationAsync(OrderId orderId, TrackingNumber trackingNumber, CancellationToken cancellationToken = default);
}

// Infrastructure/ExternalServices/SendGridEmailService.cs - Adapter in Infrastructure
namespace YourApp.Infrastructure.ExternalServices;

// Application/Exceptions/EmailServiceException.cs
public class EmailServiceException : Exception
{
    public EmailServiceException(string message) : base(message) { }
    public EmailServiceException(string message, Exception innerException) : base(message, innerException) { }
}

public class SendGridEmailService : IEmailService
{
    private readonly ISendGridClient _client;
    private readonly ILogger<SendGridEmailService> _logger;
    private readonly EmailSettings _settings;

    public SendGridEmailService(
        ISendGridClient client,
        IOptions<EmailSettings> settings,
        ILogger<SendGridEmailService> logger)
    {
        _client = client;
        _settings = settings.Value;
        _logger = logger;
    }

    public async Task SendOrderConfirmationAsync(OrderId orderId, CancellationToken cancellationToken)
    {
        try
        {
            var message = new SendGridMessage
            {
                From = new EmailAddress(_settings.FromEmail, _settings.FromName),
                Subject = $"Order Confirmation - {orderId}",
                HtmlContent = GenerateOrderConfirmationHtml(orderId)
            };

            message.AddTo(_settings.ToEmail);

            var response = await _client.SendEmailAsync(message, cancellationToken);

            if (!response.IsSuccessStatusCode)
            {
                _logger.LogError("Failed to send email. Status: {StatusCode}", response.StatusCode);
                throw new EmailServiceException("Failed to send confirmation email");
            }

            _logger.LogInformation("Order confirmation email sent for {OrderId}", orderId);
        }
        catch (SendGridException ex)
        {
            _logger.LogError(ex, "SendGrid error sending confirmation email for {OrderId}", orderId);
            throw new EmailServiceException("Email service error", ex);
        }
        catch (Exception ex)
        {
            _logger.LogError(ex, "Error sending confirmation email for {OrderId}", orderId);
            throw new EmailServiceException("Failed to send email", ex);
        }
    }

    private string GenerateOrderConfirmationHtml(OrderId orderId)
    {
        return $"<h1>Your order {orderId} has been confirmed!</h1>";
    }
}
```

### Mapper Pattern

Map between domain entities and DTOs:

```csharp
// Application/DTOs/OrderDto.cs
namespace YourApp.Application.DTOs;

public class OrderDto
{
    public Guid Id { get; set; }
    public Guid CustomerId { get; set; }
    public required string Status { get; set; }
    public decimal TotalAmount { get; set; }
    public required string TotalCurrency { get; set; }
    public DateTime CreatedAt { get; set; }
    public required List<OrderItemDto> Items { get; set; }
}

public class OrderItemDto
{
    public Guid ProductId { get; set; }
    public required string ProductName { get; set; }
    public int Quantity { get; set; }
    public decimal UnitPrice { get; set; }
    public decimal Subtotal { get; set; }
}

// Application/Mappers/OrderProfile.cs - AutoMapper profile
namespace YourApp.Application.Mappers;

public class OrderProfile : Profile
{
    public OrderProfile()
    {
        // Domain to DTO
        CreateMap<Order, OrderDto>()
            .ForMember(dest => dest.Id, opt => opt.MapFrom(src => src.Id.Value))
            .ForMember(dest => dest.CustomerId, opt => opt.MapFrom(src => src.CustomerId.Value))
            .ForMember(dest => dest.Status, opt => opt.MapFrom(src => src.Status.ToString()))
            .ForMember(dest => dest.TotalAmount, opt => opt.MapFrom(src => src.Total.Amount))
            .ForMember(dest => dest.TotalCurrency, opt => opt.MapFrom(src => src.Total.Currency.Code));

        CreateMap<OrderItem, OrderItemDto>()
            .ForMember(dest => dest.ProductId, opt => opt.MapFrom(src => src.ProductId.Value))
            .ForMember(dest => dest.ProductName, opt => opt.MapFrom(src => src.Product.Name));
    }
}

// Usage in query handler
public class GetOrderByIdQueryHandler : IRequestHandler<GetOrderByIdQuery, OrderDto>
{
    private readonly IOrderRepository _orderRepository;
    private readonly IMapper _mapper;

    public async Task<OrderDto> Handle(GetOrderByIdQuery request, CancellationToken cancellationToken)
    {
        var order = await _orderRepository.GetByIdAsync(
            OrderId.From(request.OrderId),
            cancellationToken);

        if (order == null)
            throw new NotFoundException($"Order {request.OrderId} not found");

        return _mapper.Map<OrderDto>(order);
    }
}
```

---

## Behavioral Patterns

### Strategy Pattern

Encapsulate interchangeable algorithms:

```csharp
// Domain/Services/IDiscountStrategy.cs
namespace YourApp.Domain.Services;

public interface IDiscountStrategy
{
    Money CalculateDiscount(Order order);
    bool IsApplicable(Order order);
}

// Domain/Services/PercentageDiscountStrategy.cs
public class PercentageDiscountStrategy : IDiscountStrategy
{
    private readonly decimal _percentage;

    public PercentageDiscountStrategy(decimal percentage)
    {
        _percentage = percentage;
    }

    public Money CalculateDiscount(Order order)
    {
        var discountAmount = order.Total.Amount * (_percentage / 100);
        return new Money(discountAmount, order.Total.Currency);
    }

    public bool IsApplicable(Order order)
    {
        return order.Total.Amount > 0;
    }
}

// Domain/Services/FixedAmountDiscountStrategy.cs
public class FixedAmountDiscountStrategy : IDiscountStrategy
{
    private readonly Money _discountAmount;
    private readonly Money _minimumOrderAmount;

    public FixedAmountDiscountStrategy(Money discountAmount, Money minimumOrderAmount)
    {
        _discountAmount = discountAmount;
        _minimumOrderAmount = minimumOrderAmount;
    }

    public Money CalculateDiscount(Order order)
    {
        return _discountAmount;
    }

    public bool IsApplicable(Order order)
    {
        return order.Total.IsGreaterThan(_minimumOrderAmount);
    }
}

// Domain/Services/DiscountService.cs
public class DiscountService
{
    public Money ApplyDiscount(Order order, IDiscountStrategy strategy)
    {
        // Guard clauses
        if (order == null)
            throw new ArgumentNullException(nameof(order));
        if (strategy == null)
            throw new ArgumentNullException(nameof(strategy));

        if (!strategy.IsApplicable(order))
            return Money.Zero(order.Total.Currency);

        return strategy.CalculateDiscount(order);
    }
}

// Usage:
var percentageDiscount = new PercentageDiscountStrategy(10);
var fixedDiscount = new FixedAmountDiscountStrategy(
    new Money(5, Currency.USD),
    new Money(50, Currency.USD));

var discount = discountService.ApplyDiscount(order, percentageDiscount);
```

### Command Pattern (CQRS)

Encapsulate requests as objects:

```csharp
// Application/Commands/CreateOrderCommand.cs
namespace YourApp.Application.Commands;

public record CreateOrderCommand : IRequest<OrderId>
{
    public Guid CustomerId { get; init; }
    public required List<OrderItemDto> Items { get; init; }
    public required AddressDto ShippingAddress { get; init; }
}

// Application/Commands/CreateOrderCommandHandler.cs
public class CreateOrderCommandHandler : IRequestHandler<CreateOrderCommand, OrderId>
{
    private readonly IOrderRepository _orderRepository;
    private readonly IProductRepository _productRepository;
    private readonly IUnitOfWork _unitOfWork;

    public CreateOrderCommandHandler(
        IOrderRepository orderRepository,
        IProductRepository productRepository,
        IUnitOfWork unitOfWork)
    {
        _orderRepository = orderRepository;
        _productRepository = productRepository;
        _unitOfWork = unitOfWork;
    }

    public async Task<OrderId> Handle(CreateOrderCommand request, CancellationToken cancellationToken)
    {
        // Create aggregate
        var customerId = CustomerId.From(request.CustomerId);
        var shippingAddress = ShippingAddress.From(request.ShippingAddress);
        var order = Order.Create(customerId, shippingAddress);

        // Add items
        foreach (var itemDto in request.Items)
        {
            var product = await _productRepository.GetByIdAsync(
                ProductId.From(itemDto.ProductId),
                cancellationToken);

            if (product == null)
                throw new NotFoundException($"Product {itemDto.ProductId} not found");

            order.AddItem(product, itemDto.Quantity);
        }

        // Place order
        order.PlaceOrder();

        // Persist
        await _orderRepository.AddAsync(order, cancellationToken);
        await _unitOfWork.SaveChangesAsync(cancellationToken);

        return order.Id;
    }
}

// Application/Queries/GetOrderByIdQuery.cs
public record GetOrderByIdQuery(Guid OrderId) : IRequest<OrderDto>;

public class GetOrderByIdQueryHandler : IRequestHandler<GetOrderByIdQuery, OrderDto>
{
    private readonly IOrderRepository _orderRepository;
    private readonly IMapper _mapper;

    public async Task<OrderDto> Handle(GetOrderByIdQuery request, CancellationToken cancellationToken)
    {
        var order = await _orderRepository.GetByIdAsync(
            OrderId.From(request.OrderId),
            cancellationToken);

        if (order == null)
            throw new NotFoundException($"Order {request.OrderId} not found");

        return _mapper.Map<OrderDto>(order);
    }
}
```

### Observer Pattern (Domain Events)

Decouple side effects from main flow:

```csharp
// Domain/DomainEvents/OrderPlacedEvent.cs
namespace YourApp.Domain.DomainEvents;

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

// Application/EventHandlers/OrderPlacedEventHandler.cs
namespace YourApp.Application.EventHandlers;

public class OrderPlacedEventHandler : INotificationHandler<OrderPlacedEvent>
{
    private readonly IEmailService _emailService;
    private readonly ILoyaltyService _loyaltyService;
    private readonly ILogger<OrderPlacedEventHandler> _logger;

    public async Task Handle(OrderPlacedEvent notification, CancellationToken cancellationToken)
    {
        _logger.LogInformation("Handling OrderPlacedEvent for order {OrderId}", notification.OrderId);

        // Send confirmation email
        await _emailService.SendOrderConfirmationAsync(
            notification.OrderId,
            cancellationToken);

        // Award loyalty points
        await _loyaltyService.AwardPointsAsync(
            notification.CustomerId,
            notification.Total,
            cancellationToken);

        _logger.LogInformation("OrderPlacedEvent handled successfully");
    }
}

// Domain/Entities/Order.cs - Raising events
public class Order
{
    private readonly List<IDomainEvent> _domainEvents = new();
    public IReadOnlyList<IDomainEvent> DomainEvents => _domainEvents.AsReadOnly();

    public void PlaceOrder()
    {
        if (Status != OrderStatus.Draft)
            throw new DomainException("Order is already placed");

        Status = OrderStatus.Placed;
        _domainEvents.Add(new OrderPlacedEvent(Id, CustomerId, Total));
    }

    public void ClearDomainEvents() => _domainEvents.Clear();
}
```

### Chain of Responsibility Pattern (Pipeline Behaviors)

Process requests through multiple handlers:

```csharp
// Application/Behaviors/ValidationBehavior.cs
namespace YourApp.Application.Behaviors;

public class ValidationBehavior<TRequest, TResponse> : IPipelineBehavior<TRequest, TResponse>
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
        CancellationToken cancellationToken)
    {
        if (!_validators.Any())
            return await next();

        var context = new ValidationContext<TRequest>(request);

        var validationResults = await Task.WhenAll(
            _validators.Select(v => v.ValidateAsync(context, cancellationToken)));

        var failures = validationResults
            .SelectMany(r => r.Errors)
            .Where(f => f != null)
            .ToList();

        if (failures.Any())
            throw new ValidationException(failures);

        return await next();
    }
}

// Application/Behaviors/LoggingBehavior.cs
public class LoggingBehavior<TRequest, TResponse> : IPipelineBehavior<TRequest, TResponse>
    where TRequest : IRequest<TResponse>
{
    private readonly ILogger<LoggingBehavior<TRequest, TResponse>> _logger;

    public async Task<TResponse> Handle(
        TRequest request,
        RequestHandlerDelegate<TResponse> next,
        CancellationToken cancellationToken)
    {
        var requestName = typeof(TRequest).Name;

        _logger.LogInformation("Handling {RequestName}", requestName);

        var stopwatch = Stopwatch.StartNew();

        try
        {
            var response = await next();

            stopwatch.Stop();

            _logger.LogInformation(
                "Handled {RequestName} in {ElapsedMilliseconds}ms",
                requestName,
                stopwatch.ElapsedMilliseconds);

            return response;
        }
        catch (Exception ex)
        {
            stopwatch.Stop();

            _logger.LogError(
                ex,
                "Error handling {RequestName} after {ElapsedMilliseconds}ms",
                requestName,
                stopwatch.ElapsedMilliseconds);

            throw;
        }
    }
}

// API/Program.cs - Register behaviors
builder.Services.AddMediatR(cfg =>
{
    cfg.RegisterServicesFromAssembly(typeof(CreateOrderCommand).Assembly);
    cfg.AddOpenBehavior(typeof(ValidationBehavior<,>));
    cfg.AddOpenBehavior(typeof(LoggingBehavior<,>));
});

// Execution flow: Request → ValidationBehavior → LoggingBehavior → CommandHandler
```

---

## Key Takeaways

1. **Factory Pattern** - Encapsulates creation logic, enforces invariants
2. **Builder Pattern** - Fluent API for complex object construction
3. **Dependency Injection** - Inverts dependencies, enables testability
4. **Repository Pattern** - Abstracts data access, provides collection-like interface
5. **Unit of Work** - Coordinates transactions across repositories
6. **Adapter Pattern** - Translates between application and external services
7. **Mapper Pattern** - Converts between domain and DTO representations
8. **Strategy Pattern** - Encapsulates algorithms, enables runtime selection
9. **Command Pattern** - Separates request from execution (CQRS)
10. **Observer Pattern** - Decouples side effects through domain events
11. **Chain of Responsibility** - Processes requests through pipeline of handlers
