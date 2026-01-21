# Clean Architecture - Project Structure Examples (C#)

Examples demonstrating the overall structure and layer organization for Clean Architecture in C#/.NET.

> These examples complement [SKILL.md](../../SKILL.md)

## Project Structure

### Recommended Solution Structure

```
YourApp.sln
├── src/
│   ├── YourApp.Domain/
│   │   ├── Entities/
│   │   ├── ValueObjects/
│   │   ├── Aggregates/
│   │   ├── DomainServices/
│   │   ├── DomainEvents/
│   │   ├── Repositories/           # Interfaces only
│   │   ├── Specifications/
│   │   └── Exceptions/
│   ├── YourApp.Application/
│   │   ├── Commands/
│   │   ├── Queries/
│   │   ├── DTOs/
│   │   ├── Interfaces/             # Ports
│   │   ├── Behaviors/
│   │   ├── Validators/
│   │   └── Mappers/
│   ├── YourApp.Infrastructure/
│   │   ├── Persistence/
│   │   │   ├── Configurations/
│   │   │   ├── Repositories/      # Implementations
│   │   │   ├── Migrations/
│   │   │   └── ApplicationDbContext.cs
│   │   ├── ExternalServices/
│   │   ├── Identity/
│   │   └── FileStorage/
│   └── YourApp.API/
│       ├── Controllers/
│       ├── Middleware/
│       ├── Filters/
│       └── Program.cs
└── tests/
    ├── YourApp.Domain.Tests/
    ├── YourApp.Application.Tests/
    ├── YourApp.Infrastructure.Tests/
    └── YourApp.API.Tests/
```

## Domain Layer Example

**File:** `YourApp.Domain/Entities/Order.cs`

```csharp
namespace YourApp.Domain.Entities;

public class Order
{
    public OrderId Id { get; private set; }
    public CustomerId CustomerId { get; private set; }
    public OrderStatus Status { get; private set; }
    public Money Total { get; private set; }

    private readonly List<OrderItem> _items = new();
    public IReadOnlyList<OrderItem> Items => _items.AsReadOnly();

    private Order() { }

    public static Order Create(CustomerId customerId, ShippingAddress address)
    {
        return new Order
        {
            Id = OrderId.NewId(),
            CustomerId = customerId,
            Status = OrderStatus.Draft,
            Total = Money.Zero("USD")
        };
    }

    public void AddItem(Product product, int quantity)
    {
        if (Status != OrderStatus.Draft)
            throw new DomainException("Cannot add items to non-draft order");

        _items.Add(OrderItem.Create(product, quantity));
        RecalculateTotal();
    }

    private void RecalculateTotal()
    {
        Total = _items.Aggregate(
            Money.Zero("USD"),
            (sum, item) => sum.Add(item.Subtotal)
        );
    }
}
```

**File:** `YourApp.Domain/Repositories/IOrderRepository.cs`

```csharp
namespace YourApp.Domain.Repositories;

// Interface in Domain, implementation in Infrastructure
public interface IOrderRepository
{
    Task<Order?> GetByIdAsync(OrderId id, CancellationToken cancellationToken = default);
    Task<List<Order>> GetByCustomerIdAsync(CustomerId customerId, CancellationToken cancellationToken = default);
    Task AddAsync(Order order, CancellationToken cancellationToken = default);
    Task UpdateAsync(Order order, CancellationToken cancellationToken = default);
    Task DeleteAsync(OrderId id, CancellationToken cancellationToken = default);
}
```

## Application Layer Example

**File:** `YourApp.Application/Commands/CreateOrderCommand.cs`

```csharp
namespace YourApp.Application.Commands;

public record CreateOrderCommand(
    Guid CustomerId,
    List<OrderItemDto> Items,
    AddressDto ShippingAddress
) : IRequest<Guid>;

public record OrderItemDto(Guid ProductId, int Quantity);
public record AddressDto(string Street, string City, string PostalCode);
```

**File:** `YourApp.Application/Commands/CreateOrderCommandHandler.cs`

```csharp
namespace YourApp.Application.Commands;

public class CreateOrderCommandHandler : IRequestHandler<CreateOrderCommand, Guid>
{
    private readonly IOrderRepository _orderRepository;
    private readonly IProductRepository _productRepository;

    public CreateOrderCommandHandler(
        IOrderRepository orderRepository,
        IProductRepository productRepository)
    {
        _orderRepository = orderRepository;
        _productRepository = productRepository;
    }

    public async Task<Guid> Handle(CreateOrderCommand request, CancellationToken cancellationToken)
    {
        // Create domain entity
        var order = Order.Create(
            CustomerId.From(request.CustomerId),
            ShippingAddress.From(request.ShippingAddress)
        );

        // Add items using domain logic
        foreach (var itemDto in request.Items)
        {
            var product = await _productRepository.GetByIdAsync(
                ProductId.From(itemDto.ProductId),
                cancellationToken
            );

            if (product == null)
                throw new NotFoundException($"Product {itemDto.ProductId} not found");

            order.AddItem(product, itemDto.Quantity);
        }

        // Persist through repository
        await _orderRepository.AddAsync(order, cancellationToken);

        return order.Id.Value;
    }
}
```

## Infrastructure Layer Example

**File:** `YourApp.Infrastructure/Persistence/Repositories/OrderRepository.cs`

```csharp
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
        // Load aggregate with all internal entities
        return await _context.Orders
            .Include(o => o.Items)
            .FirstOrDefaultAsync(o => o.Id == id, cancellationToken);
    }

    public async Task AddAsync(Order order, CancellationToken cancellationToken = default)
    {
        await _context.Orders.AddAsync(order, cancellationToken);
        await _context.SaveChangesAsync(cancellationToken);
    }
}
```

**File:** `YourApp.Infrastructure/Persistence/ApplicationDbContext.cs`

```csharp
namespace YourApp.Infrastructure.Persistence;

public class ApplicationDbContext : DbContext
{
    public DbSet<Order> Orders => Set<Order>();
    public DbSet<Product> Products => Set<Product>();

    public ApplicationDbContext(DbContextOptions<ApplicationDbContext> options)
        : base(options)
    {
    }

    protected override void OnModelCreating(ModelBuilder modelBuilder)
    {
        base.OnModelCreating(modelBuilder);

        // Apply all configurations from assembly
        modelBuilder.ApplyConfigurationsFromAssembly(typeof(ApplicationDbContext).Assembly);
    }
}
```

## Presentation Layer Example

**File:** `YourApp.API/Controllers/OrdersController.cs`

```csharp
namespace YourApp.API.Controllers;

[ApiController]
[Route("api/[controller]")]
public class OrdersController : ControllerBase
{
    private readonly IMediator _mediator;

    public OrdersController(IMediator mediator)
    {
        _mediator = mediator;
    }

    [HttpPost]
    public async Task<ActionResult<Guid>> CreateOrder(
        [FromBody] CreateOrderRequest request,
        CancellationToken cancellationToken)
    {
        var command = new CreateOrderCommand(
            request.CustomerId,
            request.Items,
            request.ShippingAddress
        );

        var orderId = await _mediator.Send(command, cancellationToken);

        return CreatedAtAction(
            nameof(GetOrder),
            new { id = orderId },
            orderId
        );
    }

    [HttpGet("{id}")]
    public async Task<ActionResult<OrderDto>> GetOrder(
        Guid id,
        CancellationToken cancellationToken)
    {
        var query = new GetOrderQuery(id);
        var order = await _mediator.Send(query, cancellationToken);

        return order == null ? NotFound() : Ok(order);
    }
}
```

## Dependency Injection Setup

**File:** `YourApp.API/Program.cs`

```csharp
var builder = WebApplication.CreateBuilder(args);

// Add services to the container
builder.Services.AddControllers();

// Infrastructure
builder.Services.AddDbContext<ApplicationDbContext>(options =>
    options.UseSqlServer(builder.Configuration.GetConnectionString("DefaultConnection")));

// Register repositories
builder.Services.AddScoped<IOrderRepository, OrderRepository>();
builder.Services.AddScoped<IProductRepository, ProductRepository>();

// MediatR for CQRS
builder.Services.AddMediatR(cfg =>
    cfg.RegisterServicesFromAssembly(typeof(CreateOrderCommand).Assembly));

// FluentValidation
builder.Services.AddValidatorsFromAssembly(typeof(CreateOrderCommand).Assembly);

var app = builder.Build();

app.UseHttpsRedirection();
app.UseAuthorization();
app.MapControllers();

app.Run();
```

## Package Dependencies

### Domain Project
```xml
<Project Sdk="Microsoft.NET.Sdk">
  <PropertyGroup>
    <TargetFramework>net8.0</TargetFramework>
  </PropertyGroup>
  <!-- NO DEPENDENCIES - Pure C# -->
</Project>
```

### Application Project
```xml
<ItemGroup>
  <PackageReference Include="MediatR" Version="12.2.0" />
  <PackageReference Include="FluentValidation" Version="11.9.0" />
</ItemGroup>
<ItemGroup>
  <ProjectReference Include="..\YourApp.Domain\YourApp.Domain.csproj" />
</ItemGroup>
```

### Infrastructure Project
```xml
<ItemGroup>
  <PackageReference Include="Microsoft.EntityFrameworkCore.SqlServer" Version="8.0.0" />
  <PackageReference Include="Microsoft.EntityFrameworkCore.Design" Version="8.0.0" />
</ItemGroup>
<ItemGroup>
  <ProjectReference Include="..\YourApp.Application\YourApp.Application.csproj" />
  <ProjectReference Include="..\YourApp.Domain\YourApp.Domain.csproj" />
</ItemGroup>
```

### API Project
```xml
<ItemGroup>
  <PackageReference Include="Microsoft.AspNetCore.OpenApi" Version="8.0.0" />
  <PackageReference Include="Swashbuckle.AspNetCore" Version="6.5.0" />
</ItemGroup>
<ItemGroup>
  <ProjectReference Include="..\YourApp.Application\YourApp.Application.csproj" />
  <ProjectReference Include="..\YourApp.Infrastructure\YourApp.Infrastructure.csproj" />
</ItemGroup>
```

## Key Takeaways

1. **Domain has NO dependencies** - Pure C# code
2. **Application depends ONLY on Domain** - Plus MediatR for CQRS
3. **Infrastructure depends on Application and Domain** - Implements interfaces
4. **API depends on Application and Infrastructure** - Composition root for DI

This structure ensures proper dependency flow: inward toward the Domain!
