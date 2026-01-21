# Clean Architecture - Architecture Review Examples (C#)

Complete feature example with tests demonstrating all architectural principles working together.

> These examples complement [CHECKLIST.md](../../CHECKLIST.md)

## Complete Feature: Place Order

This example shows a complete vertical slice from API to database, demonstrating proper Clean Architecture implementation.

### 1. Presentation Layer (API)

**Request Model** - API-specific DTO with validation:

```csharp
// API/Models/PlaceOrderRequest.cs
namespace YourApp.API.Models;

public class PlaceOrderRequest
{
    [Required(ErrorMessage = "Customer ID is required")]
    public Guid CustomerId { get; set; }

    [Required(ErrorMessage = "Items are required")]
    [MinLength(1, ErrorMessage = "Order must have at least one item")]
    public required List<OrderItemRequest> Items { get; set; }

    [Required(ErrorMessage = "Shipping address is required")]
    public required AddressDto ShippingAddress { get; set; }
}

public class OrderItemRequest
{
    [Required(ErrorMessage = "Product ID is required")]
    public Guid ProductId { get; set; }

    [Range(1, int.MaxValue, ErrorMessage = "Quantity must be greater than 0")]
    public int Quantity { get; set; }
}
```

**Controller** - Thin layer, delegates to application:

```csharp
// API/Controllers/OrdersController.cs
using System.ComponentModel.DataAnnotations;
using Microsoft.AspNetCore.Mvc;
using MediatR;

namespace YourApp.API.Controllers;

[ApiController]
[Route("api/[controller]")]
public class OrdersController : ControllerBase
{
    private readonly IMediator _mediator;
    private readonly ILogger<OrdersController> _logger;

    public OrdersController(IMediator mediator, ILogger<OrdersController> logger)
    {
        _mediator = mediator;
        _logger = logger;
    }

    [HttpPost]
    [ProducesResponseType(typeof(Guid), StatusCodes.Status201Created)]
    [ProducesResponseType(StatusCodes.Status400BadRequest)]
    public async Task<IActionResult> PlaceOrder(
        [FromBody] PlaceOrderRequest request,
        CancellationToken cancellationToken)
    {
        _logger.LogInformation("Placing order for customer {CustomerId}", request.CustomerId);

        var command = new PlaceOrderCommand
        {
            CustomerId = request.CustomerId,
            Items = request.Items.Select(i => new OrderItemDto
            {
                ProductId = i.ProductId,
                Quantity = i.Quantity
            }).ToList(),
            ShippingAddress = request.ShippingAddress
        };

        var orderId = await _mediator.Send(command, cancellationToken);

        return CreatedAtAction(nameof(GetOrder), new { id = orderId }, orderId);
    }

    [HttpGet("{id}")]
    [ProducesResponseType(typeof(OrderDto), StatusCodes.Status200OK)]
    [ProducesResponseType(StatusCodes.Status404NotFound)]
    public async Task<IActionResult> GetOrder(Guid id, CancellationToken cancellationToken)
    {
        var query = new GetOrderByIdQuery(id);
        var order = await _mediator.Send(query, cancellationToken);

        return Ok(order);
    }
}
```

### 2. Application Layer

**Command** - Immutable request object:

```csharp
// Application/Commands/PlaceOrderCommand.cs
namespace YourApp.Application.Commands;

public record PlaceOrderCommand : IRequest<Guid>
{
    public Guid CustomerId { get; init; }
    public required List<OrderItemDto> Items { get; init; }
    public required AddressDto ShippingAddress { get; init; }
}

public class OrderItemDto
{
    public Guid ProductId { get; set; }
    public int Quantity { get; set; }
}

public class AddressDto
{
    public required string Street { get; set; }
    public required string City { get; set; }
    public required string PostalCode { get; set; }
    public required string Country { get; set; }
}
```

**Validation** - FluentValidation:

```csharp
// Application/Validators/PlaceOrderCommandValidator.cs
namespace YourApp.Application.Validators;

public class PlaceOrderCommandValidator : AbstractValidator<PlaceOrderCommand>
{
    public PlaceOrderCommandValidator()
    {
        RuleFor(x => x.CustomerId)
            .NotEmpty()
            .WithMessage("Customer ID is required");

        RuleFor(x => x.Items)
            .NotEmpty()
            .WithMessage("Order must have at least one item");

        RuleForEach(x => x.Items)
            .SetValidator(new OrderItemDtoValidator());

        RuleFor(x => x.ShippingAddress)
            .NotNull()
            .WithMessage("Shipping address is required")
            .SetValidator(new AddressDtoValidator());
    }
}

public class OrderItemDtoValidator : AbstractValidator<OrderItemDto>
{
    public OrderItemDtoValidator()
    {
        RuleFor(x => x.ProductId)
            .NotEmpty()
            .WithMessage("Product ID is required");

        RuleFor(x => x.Quantity)
            .GreaterThan(0)
            .WithMessage("Quantity must be greater than 0");
    }
}
```

**Command Handler** - Orchestrates domain objects:

```csharp
// Application/Commands/PlaceOrderCommandHandler.cs
using Microsoft.EntityFrameworkCore;

namespace YourApp.Application.Commands;

// Application/Exceptions/ConcurrencyException.cs
public class ConcurrencyException : Exception
{
    public ConcurrencyException(string message) : base(message) { }
    public ConcurrencyException(string message, Exception innerException) : base(message, innerException) { }
}

// Application/Exceptions/NotFoundException.cs
public class NotFoundException : Exception
{
    public NotFoundException(string message) : base(message) { }
}

public class PlaceOrderCommandHandler : IRequestHandler<PlaceOrderCommand, Guid>
{
    private readonly IOrderRepository _orderRepository;
    private readonly IProductRepository _productRepository;
    private readonly IUnitOfWork _unitOfWork;
    private readonly ILogger<PlaceOrderCommandHandler> _logger;

    public PlaceOrderCommandHandler(
        IOrderRepository orderRepository,
        IProductRepository productRepository,
        IUnitOfWork unitOfWork,
        ILogger<PlaceOrderCommandHandler> logger)
    {
        _orderRepository = orderRepository;
        _productRepository = productRepository;
        _unitOfWork = unitOfWork;
        _logger = logger;
    }

    public async Task<Guid> Handle(PlaceOrderCommand request, CancellationToken cancellationToken)
    {
        _logger.LogInformation("Creating order for customer {CustomerId}", request.CustomerId);

        // 1. Create aggregate using domain factory
        var customerId = CustomerId.From(request.CustomerId);
        var shippingAddress = ShippingAddress.From(request.ShippingAddress);
        var order = Order.Create(customerId, shippingAddress);

        // 2. Add items - domain logic enforces business rules
        foreach (var itemDto in request.Items)
        {
            var product = await _productRepository.GetByIdAsync(
                ProductId.From(itemDto.ProductId),
                cancellationToken);

            if (product == null)
                throw new NotFoundException($"Product {itemDto.ProductId} not found");

            // Domain logic validates and calculates
            order.AddItem(product, itemDto.Quantity);
        }

        // 3. Place order - domain logic validates and raises event
        order.PlaceOrder();

        // 4. Persist using repository abstraction
        await _orderRepository.AddAsync(order, cancellationToken);

        try
        {
            await _unitOfWork.SaveChangesAsync(cancellationToken);
        }
        catch (DbUpdateConcurrencyException ex)
        {
            _logger.LogError(ex, "Concurrency conflict while creating order");
            throw new ConcurrencyException("The order was modified by another process. Please try again.");
        }

        _logger.LogInformation("Order {OrderId} created successfully", order.Id);

        // 5. Return identifier
        return order.Id.Value;
    }
}
```

**Domain Event Handler** - Side effects decoupled from main flow:

```csharp
// Application/EventHandlers/OrderPlacedEventHandler.cs
using System.Text.Json;

namespace YourApp.Application.EventHandlers;

// Application/Interfaces/IOutboxRepository.cs
public interface IOutboxRepository
{
    Task AddAsync(OutboxMessage message, CancellationToken cancellationToken = default);
}

// Domain/Entities/OutboxMessage.cs
public class OutboxMessage
{
    public Guid Id { get; set; } = Guid.NewGuid();
    public required string EventType { get; set; }
    public required string Payload { get; set; }
    public string? Error { get; set; }
    public int RetryCount { get; set; }
    public DateTime CreatedAt { get; set; } = DateTime.UtcNow;
    public DateTime? ProcessedAt { get; set; }
}

public class OrderPlacedEventHandler : INotificationHandler<OrderPlacedEvent>
{
    private readonly IEmailService _emailService;
    private readonly ILoyaltyService _loyaltyService;
    private readonly IOutboxRepository _outboxRepository;
    private readonly ILogger<OrderPlacedEventHandler> _logger;

    public OrderPlacedEventHandler(
        IEmailService emailService,
        ILoyaltyService loyaltyService,
        IOutboxRepository outboxRepository,
        ILogger<OrderPlacedEventHandler> logger)
    {
        _emailService = emailService;
        _loyaltyService = loyaltyService;
        _outboxRepository = outboxRepository;
        _logger = logger;
    }

    public async Task Handle(OrderPlacedEvent notification, CancellationToken cancellationToken)
    {
        _logger.LogInformation("Handling OrderPlacedEvent for order {OrderId}", notification.OrderId);

        try
        {
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
        catch (Exception ex)
        {
            _logger.LogError(ex, "Error handling OrderPlacedEvent for order {OrderId}", notification.OrderId);

            // Store failed event in outbox for retry
            await _outboxRepository.AddAsync(new OutboxMessage
            {
                EventType = nameof(OrderPlacedEvent),
                Payload = JsonSerializer.Serialize(notification),
                Error = ex.Message,
                RetryCount = 0
            }, cancellationToken);

            // Don't throw - side effects should not fail the main operation
        }
    }
}
```

### 3. Domain Layer

**Aggregate Root** - Rich behavior, encapsulated state:

```csharp
// Domain/Entities/Order.cs
using System.ComponentModel.DataAnnotations;

namespace YourApp.Domain.Entities;

// Domain/Exceptions/DomainException.cs
public class DomainException : Exception
{
    public DomainException(string message) : base(message) { }
    public DomainException(string message, Exception innerException) : base(message, innerException) { }
}

public class Order
{
    // Private setters protect invariants
    public OrderId Id { get; private set; }
    public CustomerId CustomerId { get; private set; }
    public OrderStatus Status { get; private set; }
    public Money Total { get; private set; }
    public DateTime CreatedAt { get; private set; }
    public ShippingAddress ShippingAddress { get; private set; }

    // Optimistic concurrency control
    [Timestamp]
    public byte[] RowVersion { get; private set; }

    // Encapsulated collection - readonly from outside
    private readonly List<OrderItem> _items = new();
    public IReadOnlyList<OrderItem> Items => _items.AsReadOnly();

    // Domain events
    private readonly List<IDomainEvent> _domainEvents = new();
    public IReadOnlyList<IDomainEvent> DomainEvents => _domainEvents.AsReadOnly();

    // Private constructor prevents invalid construction
    private Order() { }

    // Static factory enforces creation rules
    public static Order Create(CustomerId customerId, ShippingAddress shippingAddress)
    {
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

        order._domainEvents.Add(new OrderCreatedEvent(order.Id, order.CustomerId));
        return order;
    }

    // Business logic method - validates and enforces rules
    public void AddItem(Product product, int quantity)
    {
        // Invariant: Can only add items to draft orders
        if (Status != OrderStatus.Draft)
            throw new DomainException("Cannot add items to non-draft order");

        // Invariant: Quantity must be positive
        if (quantity <= 0)
            throw new DomainException("Quantity must be positive");

        // Business logic: Combine duplicate items
        var existingItem = _items.FirstOrDefault(i => i.ProductId == product.Id);
        if (existingItem != null)
        {
            existingItem.IncreaseQuantity(quantity);
        }
        else
        {
            _items.Add(OrderItem.Create(product, quantity));
        }

        // Maintain invariant: Total always reflects items
        RecalculateTotal();
    }

    // Business logic method - validates and raises event
    public void PlaceOrder()
    {
        // Invariant: Can only place draft orders
        if (Status != OrderStatus.Draft)
            throw new DomainException("Order is already placed");

        // Invariant: Cannot place empty orders
        if (!_items.Any())
            throw new DomainException("Cannot place empty order");

        Status = OrderStatus.Placed;

        // Raise domain event for side effects
        _domainEvents.Add(new OrderPlacedEvent(Id, CustomerId, Total));
    }

    // Private method encapsulates calculation logic
    private void RecalculateTotal()
    {
        // Use the order's currency, not hardcoded USD
        var currency = Total?.Currency ?? Currency.USD;
        Total = _items.Aggregate(
            Money.Zero(currency),
            (sum, item) => sum.Add(item.Subtotal)
        );
    }

    public void ClearDomainEvents() => _domainEvents.Clear();
}
```

**Domain Event** - Immutable, past tense:

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
        OrderId = orderId ?? throw new ArgumentNullException(nameof(orderId));
        CustomerId = customerId ?? throw new ArgumentNullException(nameof(customerId));
        Total = total ?? throw new ArgumentNullException(nameof(total));
        OccurredAt = DateTime.UtcNow;
    }
}
```

### 4. Infrastructure Layer

**Repository Implementation** - Hides persistence details:

```csharp
// Infrastructure/Persistence/Repositories/OrderRepository.cs
namespace YourApp.Infrastructure.Persistence.Repositories;

public class OrderRepository : IOrderRepository
{
    private readonly ApplicationDbContext _context;
    private readonly ILogger<OrderRepository> _logger;

    public OrderRepository(ApplicationDbContext context, ILogger<OrderRepository> logger)
    {
        _context = context;
        _logger = logger;
    }

    public async Task<Order?> GetByIdAsync(OrderId id, CancellationToken cancellationToken = default)
    {
        _logger.LogDebug("Getting order {OrderId}", id);

        return await _context.Orders
            .AsNoTracking() // Read-only query optimization
            .AsSplitQuery() // Avoid cartesian explosion with multiple includes
            .Include(o => o.Items)
                .ThenInclude(i => i.Product)
            .FirstOrDefaultAsync(o => o.Id == id, cancellationToken);
    }

    public async Task AddAsync(Order order, CancellationToken cancellationToken = default)
    {
        _logger.LogDebug("Adding order {OrderId}", order.Id);

        await _context.Orders.AddAsync(order, cancellationToken);
    }

    public Task UpdateAsync(Order order, CancellationToken cancellationToken = default)
    {
        _logger.LogDebug("Updating order {OrderId}", order.Id);

        _context.Orders.Update(order);
        return Task.CompletedTask;
    }
}
```

---

## Testing Examples

### Domain Tests (Unit Tests)

Pure unit tests - no mocks, no database:

```csharp
// Domain.Tests/OrderTests.cs
namespace YourApp.Domain.Tests;

public class OrderTests
{
    [Fact]
    public void Create_ValidInputs_ShouldCreateOrderInDraftStatus()
    {
        // Arrange
        var customerId = CustomerId.NewId();
        var address = CreateTestAddress();

        // Act
        var order = Order.Create(customerId, address);

        // Assert
        order.Should().NotBeNull();
        order.Status.Should().Be(OrderStatus.Draft);
        order.CustomerId.Should().Be(customerId);
        order.Total.Should().Be(Money.Zero(Currency.USD));
        order.Items.Should().BeEmpty();
    }

    [Fact]
    public void AddItem_ValidProduct_ShouldAddItemAndRecalculateTotal()
    {
        // Arrange
        var order = CreateTestOrder();
        var product = CreateTestProduct(name: "Widget", price: 10.00m);

        // Act
        order.AddItem(product, quantity: 2);

        // Assert
        order.Items.Should().HaveCount(1);
        order.Items[0].ProductId.Should().Be(product.Id);
        order.Items[0].Quantity.Should().Be(2);
        order.Total.Amount.Should().Be(20.00m);
    }

    [Fact]
    public void AddItem_DuplicateProduct_ShouldIncreaseQuantity()
    {
        // Arrange
        var order = CreateTestOrder();
        var product = CreateTestProduct(price: 10.00m);

        // Act
        order.AddItem(product, quantity: 2);
        order.AddItem(product, quantity: 3);

        // Assert
        order.Items.Should().HaveCount(1);
        order.Items[0].Quantity.Should().Be(5);
        order.Total.Amount.Should().Be(50.00m);
    }

    [Fact]
    public void AddItem_ToPlacedOrder_ShouldThrowDomainException()
    {
        // Arrange
        var order = CreateTestOrder();
        var product = CreateTestProduct();
        order.AddItem(product, 1);
        order.PlaceOrder(); // Now order is placed

        // Act
        Action act = () => order.AddItem(product, 1);

        // Assert
        act.Should().Throw<DomainException>()
            .WithMessage("Cannot add items to non-draft order");
    }

    [Fact]
    public void PlaceOrder_WithNoItems_ShouldThrowDomainException()
    {
        // Arrange
        var order = CreateTestOrder();

        // Act
        Action act = () => order.PlaceOrder();

        // Assert
        act.Should().Throw<DomainException>()
            .WithMessage("Cannot place empty order");
    }

    [Fact]
    public void PlaceOrder_ValidOrder_ShouldChangeStatusAndRaiseEvent()
    {
        // Arrange
        var order = CreateTestOrder();
        var product = CreateTestProduct();
        order.AddItem(product, 1);

        // Act
        order.PlaceOrder();

        // Assert
        order.Status.Should().Be(OrderStatus.Placed);
        order.DomainEvents.Should().ContainSingle()
            .Which.Should().BeOfType<OrderPlacedEvent>();
    }

    [Fact]
    public void PlaceOrder_AlreadyPlaced_ShouldThrowDomainException()
    {
        // Arrange
        var order = CreateTestOrder();
        var product = CreateTestProduct();
        order.AddItem(product, 1);
        order.PlaceOrder();

        // Act
        Action act = () => order.PlaceOrder();

        // Assert
        act.Should().Throw<DomainException>()
            .WithMessage("Order is already placed");
    }

    // Test helpers
    private Order CreateTestOrder()
    {
        return Order.Create(CustomerId.NewId(), CreateTestAddress());
    }

    private ShippingAddress CreateTestAddress()
    {
        return new ShippingAddress(
            "123 Main St",
            "Anytown",
            "12345",
            "USA");
    }

    private Product CreateTestProduct(string name = "Test Product", decimal price = 10.00m)
    {
        return Product.Create(name, new Money(price, Currency.USD));
    }
}
```

### Application Tests (Integration Tests with Mocks)

Test orchestration logic:

```csharp
// Application.Tests/PlaceOrderCommandHandlerTests.cs
namespace YourApp.Application.Tests;

public class PlaceOrderCommandHandlerTests
{
    private readonly Mock<IOrderRepository> _orderRepositoryMock;
    private readonly Mock<IProductRepository> _productRepositoryMock;
    private readonly Mock<IUnitOfWork> _unitOfWorkMock;
    private readonly Mock<ILogger<PlaceOrderCommandHandler>> _loggerMock;
    private readonly PlaceOrderCommandHandler _handler;

    public PlaceOrderCommandHandlerTests()
    {
        _orderRepositoryMock = new Mock<IOrderRepository>();
        _productRepositoryMock = new Mock<IProductRepository>();
        _unitOfWorkMock = new Mock<IUnitOfWork>();
        _loggerMock = new Mock<ILogger<PlaceOrderCommandHandler>>();

        _handler = new PlaceOrderCommandHandler(
            _orderRepositoryMock.Object,
            _productRepositoryMock.Object,
            _unitOfWorkMock.Object,
            _loggerMock.Object);
    }

    [Fact]
    public async Task Handle_ValidCommand_ShouldCreateAndPersistOrder()
    {
        // Arrange
        var productId = Guid.NewGuid();
        var command = new PlaceOrderCommand
        {
            CustomerId = Guid.NewGuid(),
            Items = new List<OrderItemDto>
            {
                new() { ProductId = productId, Quantity = 2 }
            },
            ShippingAddress = CreateTestAddressDto()
        };

        var product = CreateTestProduct(productId);
        _productRepositoryMock
            .Setup(x => x.GetByIdAsync(It.IsAny<ProductId>(), It.IsAny<CancellationToken>()))
            .ReturnsAsync(product);

        // Act
        var result = await _handler.Handle(command, CancellationToken.None);

        // Assert
        result.Should().NotBeEmpty();

        _orderRepositoryMock.Verify(
            x => x.AddAsync(
                It.Is<Order>(o => o.CustomerId.Value == command.CustomerId),
                It.IsAny<CancellationToken>()),
            Times.Once);

        _unitOfWorkMock.Verify(
            x => x.SaveChangesAsync(It.IsAny<CancellationToken>()),
            Times.Once);
    }

    [Fact]
    public async Task Handle_ProductNotFound_ShouldThrowNotFoundException()
    {
        // Arrange
        var command = new PlaceOrderCommand
        {
            CustomerId = Guid.NewGuid(),
            Items = new List<OrderItemDto>
            {
                new() { ProductId = Guid.NewGuid(), Quantity = 1 }
            },
            ShippingAddress = CreateTestAddressDto()
        };

        _productRepositoryMock
            .Setup(x => x.GetByIdAsync(It.IsAny<ProductId>(), It.IsAny<CancellationToken>()))
            .ReturnsAsync((Product)null);

        // Act
        Func<Task> act = async () => await _handler.Handle(command, CancellationToken.None);

        // Assert
        await act.Should().ThrowAsync<NotFoundException>();

        _orderRepositoryMock.Verify(
            x => x.AddAsync(It.IsAny<Order>(), It.IsAny<CancellationToken>()),
            Times.Never);
    }

    private AddressDto CreateTestAddressDto()
    {
        return new AddressDto
        {
            Street = "123 Main St",
            City = "Anytown",
            PostalCode = "12345",
            Country = "USA"
        };
    }

    private Product CreateTestProduct(Guid id, decimal price = 10.00m)
    {
        return Product.Create("Test Product", new Money(price, Currency.USD));
    }
}
```

### Infrastructure Tests (Integration Tests with Real Database)

Test against real dependencies:

```csharp
// Infrastructure.Tests/OrderRepositoryTests.cs
namespace YourApp.Infrastructure.Tests;

public class OrderRepositoryTests : IClassFixture<DatabaseFixture>
{
    private readonly ApplicationDbContext _context;
    private readonly OrderRepository _repository;

    public OrderRepositoryTests(DatabaseFixture fixture)
    {
        _context = fixture.CreateContext();
        _repository = new OrderRepository(_context, Mock.Of<ILogger<OrderRepository>>());
    }

    [Fact]
    public async Task GetByIdAsync_ExistingOrder_ShouldReturnOrderWithItems()
    {
        // Arrange
        var order = CreateTestOrder();
        await _context.Orders.AddAsync(order);
        await _context.SaveChangesAsync();

        // Act
        var result = await _repository.GetByIdAsync(order.Id);

        // Assert
        result.Should().NotBeNull();
        result.Id.Should().Be(order.Id);
        result.Items.Should().HaveCount(order.Items.Count);
    }

    [Fact]
    public async Task AddAsync_ValidOrder_ShouldPersistToDatabase()
    {
        // Arrange
        var order = CreateTestOrder();

        // Act
        await _repository.AddAsync(order);
        await _context.SaveChangesAsync();

        // Assert
        var saved = await _context.Orders.FindAsync(order.Id);
        saved.Should().NotBeNull();
        saved.Id.Should().Be(order.Id);
    }

    private Order CreateTestOrder()
    {
        var customerId = CustomerId.NewId();
        var address = new ShippingAddress("123 Main St", "Anytown", "12345", "USA");
        return Order.Create(customerId, address);
    }
}
```

---

## Key Architectural Achievements

This example demonstrates:

1. ✅ **Dependency Rule** - Dependencies point inward (Infrastructure → Application → Domain)
2. ✅ **Rich Domain Model** - Business logic in entities, not anemic data containers
3. ✅ **SOLID Principles** - SRP, OCP, DIP all applied
4. ✅ **DDD Tactical Patterns** - Entities, Value Objects, Aggregates, Domain Events, Repositories
5. ✅ **CQRS** - Separate commands and queries
6. ✅ **Repository Pattern** - Interface in Domain, implementation in Infrastructure
7. ✅ **Unit of Work** - Transactional consistency
8. ✅ **Domain Events** - Side effects decoupled from main flow
9. ✅ **Validation** - FluentValidation in Application layer
10. ✅ **Testing** - Domain (pure unit), Application (mocked), Infrastructure (real DB), API (integration)
11. ✅ **Logging** - Cross-cutting concern handled properly
12. ✅ **Exception Handling** - Domain exceptions for business rule violations

---

## Advanced Testing Examples

### API Integration Tests

Test the full HTTP pipeline:

```csharp
// API.Tests/OrdersControllerIntegrationTests.cs
namespace YourApp.API.Tests;

public class OrdersControllerIntegrationTests : IClassFixture<WebApplicationFactory<Program>>
{
    private readonly WebApplicationFactory<Program> _factory;
    private readonly HttpClient _client;

    public OrdersControllerIntegrationTests(WebApplicationFactory<Program> factory)
    {
        _factory = factory;
        _client = factory.CreateClient();
    }

    [Fact]
    public async Task PlaceOrder_ValidRequest_ReturnsCreatedWithOrderId()
    {
        // Arrange
        var request = new
        {
            CustomerId = Guid.NewGuid(),
            Items = new[]
            {
                new { ProductId = Guid.NewGuid(), Quantity = 2 }
            },
            ShippingAddress = new
            {
                Street = "123 Main St",
                City = "TestCity",
                PostalCode = "12345",
                Country = "USA"
            }
        };

        // Act
        var response = await _client.PostAsJsonAsync("/api/orders", request);

        // Assert
        response.StatusCode.Should().Be(HttpStatusCode.Created);
        var orderId = await response.Content.ReadFromJsonAsync<Guid>();
        orderId.Should().NotBeEmpty();
    }

    [Fact]
    public async Task PlaceOrder_EmptyItems_ReturnsBadRequest()
    {
        // Arrange
        var request = new
        {
            CustomerId = Guid.NewGuid(),
            Items = Array.Empty<object>(),
            ShippingAddress = new { Street = "123 Main St", City = "TestCity", PostalCode = "12345", Country = "USA" }
        };

        // Act
        var response = await _client.PostAsJsonAsync("/api/orders", request);

        // Assert
        response.StatusCode.Should().Be(HttpStatusCode.BadRequest);
    }
}
```

### Testing Domain Events

Verify events are raised and handled:

```csharp
// Domain.Tests/OrderDomainEventsTests.cs
public class OrderDomainEventsTests
{
    [Fact]
    public void PlaceOrder_ShouldRaiseOrderPlacedEvent()
    {
        // Arrange
        var order = CreateTestOrder();
        var product = CreateTestProduct();
        order.AddItem(product, 1);

        // Act
        order.PlaceOrder();

        // Assert
        var domainEvent = order.DomainEvents.Should().ContainSingle()
            .Which.Should().BeOfType<OrderPlacedEvent>().Subject;

        domainEvent.OrderId.Should().Be(order.Id);
        domainEvent.CustomerId.Should().Be(order.CustomerId);
        domainEvent.Total.Should().Be(order.Total);
    }

    [Fact]
    public void Cancel_ShouldRaiseOrderCancelledEvent()
    {
        // Arrange
        var order = CreatePlacedOrder();
        var reason = new CancellationReason("Customer request");

        // Act
        order.Cancel(reason);

        // Assert
        order.DomainEvents.Should().Contain(e =>
            e is OrderCancelledEvent cancelledEvent &&
            cancelledEvent.OrderId == order.Id &&
            cancelledEvent.Reason == reason);
    }
}

// Application.Tests/OrderPlacedEventHandlerTests.cs
public class OrderPlacedEventHandlerTests
{
    [Fact]
    public async Task Handle_OrderPlacedEvent_ShouldSendEmail()
    {
        // Arrange
        var emailServiceMock = new Mock<IEmailService>();
        var handler = new OrderPlacedEventHandler(
            emailServiceMock.Object,
            Mock.Of<ILoyaltyService>(),
            Mock.Of<ILogger<OrderPlacedEventHandler>>());

        var @event = new OrderPlacedEvent(
            OrderId.NewId(),
            CustomerId.NewId(),
            new Money(100, Currency.USD));

        // Act
        await handler.Handle(@event, CancellationToken.None);

        // Assert
        emailServiceMock.Verify(
            x => x.SendOrderConfirmationAsync(
                @event.OrderId,
                It.IsAny<CancellationToken>()),
            Times.Once);
    }
}
```

### Test Doubles Guide

When to use Mock vs Stub vs Fake:

```csharp
// MOCK - Verify behavior (interactions)
// Use when: Testing that dependencies are called correctly
[Fact]
public async Task CreateOrder_ShouldCallRepository()
{
    var repositoryMock = new Mock<IOrderRepository>();

    var handler = new CreateOrderCommandHandler(
        repositoryMock.Object,
        Mock.Of<IProductRepository>(),
        Mock.Of<IUnitOfWork>());

    await handler.Handle(new CreateOrderCommand { /*...*/ }, CancellationToken.None);

    // Verify the interaction happened
    repositoryMock.Verify(
        x => x.AddAsync(It.IsAny<Order>(), It.IsAny<CancellationToken>()),
        Times.Once);
}

// STUB - Provide predetermined responses
// Use when: Need dependencies to return specific data
[Fact]
public async Task CreateOrder_ProductNotFound_ShouldThrow()
{
    var productRepositoryStub = new Mock<IProductRepository>();
    // Stub returns null (product not found)
    productRepositoryStub
        .Setup(x => x.GetByIdAsync(It.IsAny<ProductId>(), It.IsAny<CancellationToken>()))
        .ReturnsAsync((Product)null);

    var handler = new CreateOrderCommandHandler(
        Mock.Of<IOrderRepository>(),
        productRepositoryStub.Object,
        Mock.Of<IUnitOfWork>());

    // Act & Assert
    await FluentActions
        .Invoking(() => handler.Handle(new CreateOrderCommand { /*...*/ }, CancellationToken.None))
        .Should().ThrowAsync<NotFoundException>();
}

// FAKE - Working implementation with shortcuts
// Use when: Need realistic behavior but not production dependencies
public class InMemoryOrderRepository : IOrderRepository
{
    private readonly Dictionary<OrderId, Order> _orders = new();

    public Task<Order?> GetByIdAsync(OrderId id, CancellationToken ct = default)
    {
        _orders.TryGetValue(id, out var order);
        return Task.FromResult(order);
    }

    public Task AddAsync(Order order, CancellationToken ct = default)
    {
        _orders[order.Id] = order;
        return Task.CompletedTask;
    }
}

[Fact]
public async Task CreateOrder_UsingFakeRepository_ShouldPersist()
{
    // Fake repository - real behavior, in-memory storage
    var fakeRepository = new InMemoryOrderRepository();

    var handler = new CreateOrderCommandHandler(
        fakeRepository,
        Mock.Of<IProductRepository>(),
        Mock.Of<IUnitOfWork>());

    var command = new CreateOrderCommand { /*...*/ };
    var orderId = await handler.Handle(command, CancellationToken.None);

    // Verify using the same fake repository
    var savedOrder = await fakeRepository.GetByIdAsync(orderId);
    savedOrder.Should().NotBeNull();
}
```

### Test Pyramid Summary

```
        ┌─────────────┐
       /   E2E Tests   \      ← Few (5-10) - Full HTTP→DB flow
      /_________________\      - Use WebApplicationFactory
     /                   \     - Slow, high value
    / Integration Tests   \   ← Some (20-30) - API + mocked deps
   /_______________________\   - Use Mocks for external services
  /                         \  - Medium speed
 /   Application Tests       \ ← Many (50-100) - Use case orchestration
/_____________________________\ - All dependencies mocked
/                             \  - Fast
/       Domain Tests           \ ← Most (100+) - Business logic
/_______________________________\ - Zero mocks, pure unit tests
                                 - Fastest, highest value
```

**Key Testing Principles:**
1. **Domain tests** - No mocks, pure business logic
2. **Application tests** - Mock repositories/services (behavior verification)
3. **Infrastructure tests** - Real database (test container/in-memory)
4. **API tests** - Full pipeline with WebApplicationFactory
5. **Use Mocks** for behavior verification (was X called?)
6. **Use Stubs** for state-based testing (return specific data)
7. **Use Fakes** for complex dependencies needing realistic behavior
