````markdown
# Clean Architecture - Error Handling Examples (C#)

Comprehensive error handling strategy for Clean Architecture applications.

> These examples complement [SKILL.md](../../SKILL.md) and [CHECKLIST.md](../../CHECKLIST.md)

## Error Handling Strategy Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                      Error Flow                                  │
├─────────────────────────────────────────────────────────────────┤
│  Domain Layer        →  Application Layer  →  Presentation      │
│  ───────────────        ─────────────────     ─────────────     │
│  DomainException        ValidationException   HTTP 400/404/500  │
│  BusinessRuleViolation  NotFoundException     ProblemDetails    │
│  InvariantViolation     ConflictException     Error Response    │
└─────────────────────────────────────────────────────────────────┘
```

---

## 1. Domain Layer Exceptions

Domain exceptions represent **business rule violations**. They should be:
- Specific to business concepts
- Self-documenting
- Never catch infrastructure errors

### Base Domain Exception

```csharp
// Domain/Exceptions/DomainException.cs
namespace YourApp.Domain.Exceptions;

/// <summary>
/// Base exception for all domain/business rule violations.
/// </summary>
public class DomainException : Exception
{
    public string Code { get; }

    public DomainException(string message) : base(message)
    {
        Code = "DOMAIN_ERROR";
    }

    public DomainException(string code, string message) : base(message)
    {
        Code = code;
    }

    public DomainException(string message, Exception innerException) 
        : base(message, innerException)
    {
        Code = "DOMAIN_ERROR";
    }
}
```

### Specific Domain Exceptions

```csharp
// Domain/Exceptions/OrderExceptions.cs
namespace YourApp.Domain.Exceptions;

/// <summary>
/// Thrown when attempting to modify a non-draft order.
/// </summary>
public class OrderNotModifiableException : DomainException
{
    public OrderId OrderId { get; }
    public OrderStatus CurrentStatus { get; }

    public OrderNotModifiableException(OrderId orderId, OrderStatus status)
        : base("ORDER_NOT_MODIFIABLE", 
               $"Order {orderId} cannot be modified. Current status: {status}")
    {
        OrderId = orderId;
        CurrentStatus = status;
    }
}

/// <summary>
/// Thrown when order cannot be placed due to business rules.
/// </summary>
public class OrderCannotBePlacedException : DomainException
{
    public OrderId OrderId { get; }
    public string Reason { get; }

    public OrderCannotBePlacedException(OrderId orderId, string reason)
        : base("ORDER_CANNOT_BE_PLACED", 
               $"Order {orderId} cannot be placed: {reason}")
    {
        OrderId = orderId;
        Reason = reason;
    }
}

/// <summary>
/// Thrown when order item quantity is invalid.
/// </summary>
public class InvalidQuantityException : DomainException
{
    public int Quantity { get; }

    public InvalidQuantityException(int quantity)
        : base("INVALID_QUANTITY", 
               $"Quantity must be greater than 0. Provided: {quantity}")
    {
        Quantity = quantity;
    }
}
```

### Value Object Validation Exceptions

```csharp
// Domain/Exceptions/ValueObjectExceptions.cs
namespace YourApp.Domain.Exceptions;

/// <summary>
/// Thrown when email format is invalid.
/// </summary>
public class InvalidEmailException : DomainException
{
    public string ProvidedEmail { get; }

    public InvalidEmailException(string email)
        : base("INVALID_EMAIL", $"'{email}' is not a valid email address")
    {
        ProvidedEmail = email;
    }
}

/// <summary>
/// Thrown when money amount is negative.
/// </summary>
public class NegativeMoneyException : DomainException
{
    public decimal Amount { get; }

    public NegativeMoneyException(decimal amount)
        : base("NEGATIVE_MONEY", $"Money amount cannot be negative: {amount}")
    {
        Amount = amount;
    }
}

/// <summary>
/// Thrown when currency mismatch occurs in money operations.
/// </summary>
public class CurrencyMismatchException : DomainException
{
    public string Currency1 { get; }
    public string Currency2 { get; }

    public CurrencyMismatchException(string currency1, string currency2)
        : base("CURRENCY_MISMATCH", 
               $"Cannot perform operation between {currency1} and {currency2}")
    {
        Currency1 = currency1;
        Currency2 = currency2;
    }
}
```

### Using Domain Exceptions in Entities

```csharp
// Domain/Entities/Order.cs
namespace YourApp.Domain.Entities;

public class Order
{
    public void AddItem(Product product, int quantity)
    {
        // Validate invariant - throws domain exception
        if (Status != OrderStatus.Draft)
            throw new OrderNotModifiableException(Id, Status);

        if (quantity <= 0)
            throw new InvalidQuantityException(quantity);

        if (product.IsDiscontinued)
            throw new DomainException("PRODUCT_DISCONTINUED",
                $"Product {product.Name} is discontinued and cannot be ordered");

        _items.Add(OrderItem.Create(product, quantity));
        RecalculateTotal();
    }

    public void PlaceOrder()
    {
        if (!_items.Any())
            throw new OrderCannotBePlacedException(Id, "Order has no items");

        if (Total.Amount <= 0)
            throw new OrderCannotBePlacedException(Id, "Order total must be greater than zero");

        Status = OrderStatus.Placed;
        PlacedAt = DateTime.UtcNow;
        _domainEvents.Add(new OrderPlacedEvent(Id, CustomerId, Total));
    }
}
```

---

## 2. Application Layer Exceptions

Application exceptions handle **use case failures** - validation, not found, conflicts.

### Application Exception Types

```csharp
// Application/Exceptions/ApplicationException.cs
namespace YourApp.Application.Exceptions;

/// <summary>
/// Base exception for application layer errors.
/// </summary>
public abstract class ApplicationException : Exception
{
    public string Code { get; }

    protected ApplicationException(string code, string message) : base(message)
    {
        Code = code;
    }

    protected ApplicationException(string code, string message, Exception inner) 
        : base(message, inner)
    {
        Code = code;
    }
}

/// <summary>
/// Thrown when a requested entity is not found.
/// </summary>
public class NotFoundException : ApplicationException
{
    public string EntityType { get; }
    public string EntityId { get; }

    public NotFoundException(string entityType, string entityId)
        : base("NOT_FOUND", $"{entityType} with ID '{entityId}' was not found")
    {
        EntityType = entityType;
        EntityId = entityId;
    }

    public static NotFoundException For<T>(object id) =>
        new(typeof(T).Name, id.ToString() ?? "unknown");
}

/// <summary>
/// Thrown when input validation fails.
/// </summary>
public class ValidationException : ApplicationException
{
    public IDictionary<string, string[]> Errors { get; }

    public ValidationException(IDictionary<string, string[]> errors)
        : base("VALIDATION_FAILED", "One or more validation errors occurred")
    {
        Errors = errors;
    }

    public ValidationException(string field, string error)
        : base("VALIDATION_FAILED", error)
    {
        Errors = new Dictionary<string, string[]>
        {
            { field, new[] { error } }
        };
    }
}

/// <summary>
/// Thrown when operation conflicts with current state.
/// </summary>
public class ConflictException : ApplicationException
{
    public ConflictException(string message)
        : base("CONFLICT", message)
    {
    }
}

/// <summary>
/// Thrown when concurrent modification is detected.
/// </summary>
public class ConcurrencyException : ApplicationException
{
    public ConcurrencyException(string message)
        : base("CONCURRENCY_CONFLICT", message)
    {
    }

    public ConcurrencyException(string message, Exception inner)
        : base("CONCURRENCY_CONFLICT", message, inner)
    {
    }
}

/// <summary>
/// Thrown when user is not authorized for operation.
/// </summary>
public class ForbiddenException : ApplicationException
{
    public ForbiddenException(string message = "You are not authorized to perform this action")
        : base("FORBIDDEN", message)
    {
    }
}
```

### Using Application Exceptions in Handlers

```csharp
// Application/Commands/PlaceOrderCommandHandler.cs
namespace YourApp.Application.Commands;

public class PlaceOrderCommandHandler : IRequestHandler<PlaceOrderCommand, Guid>
{
    private readonly IOrderRepository _orderRepository;
    private readonly IProductRepository _productRepository;
    private readonly IUnitOfWork _unitOfWork;
    private readonly ILogger<PlaceOrderCommandHandler> _logger;

    public async Task<Guid> Handle(PlaceOrderCommand request, CancellationToken ct)
    {
        // 1. Validate input - throw ValidationException
        if (request.Items == null || !request.Items.Any())
            throw new ValidationException("Items", "At least one item is required");

        // 2. Create domain objects
        var customerId = CustomerId.From(request.CustomerId);
        var shippingAddress = ShippingAddress.From(request.ShippingAddress);
        var order = Order.Create(customerId, shippingAddress);

        // 3. Add items - check existence
        foreach (var itemDto in request.Items)
        {
            var product = await _productRepository.GetByIdAsync(
                ProductId.From(itemDto.ProductId), ct);

            // Throw NotFoundException if product doesn't exist
            if (product == null)
                throw NotFoundException.For<Product>(itemDto.ProductId);

            // Domain exception may be thrown here (e.g., product discontinued)
            order.AddItem(product, itemDto.Quantity);
        }

        // 4. Place order - domain exception may be thrown
        order.PlaceOrder();

        // 5. Persist
        await _orderRepository.AddAsync(order, ct);

        try
        {
            await _unitOfWork.SaveChangesAsync(ct);
        }
        catch (DbUpdateConcurrencyException ex)
        {
            _logger.LogError(ex, "Concurrency conflict for order {OrderId}", order.Id);
            throw new ConcurrencyException(
                "The data was modified by another process. Please retry.", ex);
        }

        return order.Id.Value;
    }
}
```

---

## 3. FluentValidation Integration

Use FluentValidation for command/query validation before reaching handlers.

### Validators

```csharp
// Application/Validators/PlaceOrderCommandValidator.cs
namespace YourApp.Application.Validators;

public class PlaceOrderCommandValidator : AbstractValidator<PlaceOrderCommand>
{
    public PlaceOrderCommandValidator()
    {
        RuleFor(x => x.CustomerId)
            .NotEmpty()
            .WithMessage("Customer ID is required")
            .WithErrorCode("CUSTOMER_REQUIRED");

        RuleFor(x => x.Items)
            .NotEmpty()
            .WithMessage("Order must have at least one item")
            .WithErrorCode("ITEMS_REQUIRED");

        RuleForEach(x => x.Items)
            .SetValidator(new OrderItemDtoValidator());

        RuleFor(x => x.ShippingAddress)
            .NotNull()
            .WithMessage("Shipping address is required")
            .SetValidator(new AddressDtoValidator()!);
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
            .WithMessage("Quantity must be greater than 0")
            .LessThanOrEqualTo(100)
            .WithMessage("Quantity cannot exceed 100 items");
    }
}

public class AddressDtoValidator : AbstractValidator<AddressDto>
{
    public AddressDtoValidator()
    {
        RuleFor(x => x.Street)
            .NotEmpty().WithMessage("Street is required")
            .MaximumLength(200).WithMessage("Street cannot exceed 200 characters");

        RuleFor(x => x.City)
            .NotEmpty().WithMessage("City is required")
            .MaximumLength(100).WithMessage("City cannot exceed 100 characters");

        RuleFor(x => x.PostalCode)
            .NotEmpty().WithMessage("Postal code is required")
            .Matches(@"^\d{5}(-\d{4})?$").WithMessage("Invalid postal code format");

        RuleFor(x => x.Country)
            .NotEmpty().WithMessage("Country is required")
            .Length(2).WithMessage("Country must be 2-letter ISO code");
    }
}
```

### Validation Behavior (MediatR Pipeline)

```csharp
// Application/Behaviors/ValidationBehavior.cs
namespace YourApp.Application.Behaviors;

public class ValidationBehavior<TRequest, TResponse> 
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

        if (failures.Count != 0)
        {
            var errors = failures
                .GroupBy(f => f.PropertyName)
                .ToDictionary(
                    g => g.Key,
                    g => g.Select(f => f.ErrorMessage).ToArray());

            throw new ValidationException(errors);
        }

        return await next();
    }
}
```

---

## 4. Presentation Layer - Global Exception Handler

### Exception Handling Middleware

```csharp
// API/Middleware/ExceptionHandlingMiddleware.cs
namespace YourApp.API.Middleware;

public class ExceptionHandlingMiddleware
{
    private readonly RequestDelegate _next;
    private readonly ILogger<ExceptionHandlingMiddleware> _logger;

    public ExceptionHandlingMiddleware(
        RequestDelegate next,
        ILogger<ExceptionHandlingMiddleware> logger)
    {
        _next = next;
        _logger = logger;
    }

    public async Task InvokeAsync(HttpContext context)
    {
        try
        {
            await _next(context);
        }
        catch (Exception ex)
        {
            await HandleExceptionAsync(context, ex);
        }
    }

    private async Task HandleExceptionAsync(HttpContext context, Exception exception)
    {
        var (statusCode, problemDetails) = exception switch
        {
            // Application Exceptions
            ValidationException ex => (
                StatusCodes.Status400BadRequest,
                CreateValidationProblem(ex)),

            NotFoundException ex => (
                StatusCodes.Status404NotFound,
                CreateProblem(ex.Code, ex.Message, StatusCodes.Status404NotFound)),

            ConflictException ex => (
                StatusCodes.Status409Conflict,
                CreateProblem(ex.Code, ex.Message, StatusCodes.Status409Conflict)),

            ConcurrencyException ex => (
                StatusCodes.Status409Conflict,
                CreateProblem(ex.Code, ex.Message, StatusCodes.Status409Conflict)),

            ForbiddenException ex => (
                StatusCodes.Status403Forbidden,
                CreateProblem(ex.Code, ex.Message, StatusCodes.Status403Forbidden)),

            // Domain Exceptions - map to 400 Bad Request (business rule violation)
            DomainException ex => (
                StatusCodes.Status400BadRequest,
                CreateProblem(ex.Code, ex.Message, StatusCodes.Status400BadRequest)),

            // Unexpected errors - 500
            _ => (
                StatusCodes.Status500InternalServerError,
                CreateProblem("INTERNAL_ERROR", "An unexpected error occurred", 
                    StatusCodes.Status500InternalServerError))
        };

        // Log based on severity
        if (statusCode >= 500)
        {
            _logger.LogError(exception, "Unhandled exception: {Message}", exception.Message);
        }
        else
        {
            _logger.LogWarning("Handled exception: {Type} - {Message}", 
                exception.GetType().Name, exception.Message);
        }

        context.Response.StatusCode = statusCode;
        context.Response.ContentType = "application/problem+json";

        await context.Response.WriteAsJsonAsync(problemDetails);
    }

    private static ProblemDetails CreateProblem(string code, string detail, int status)
    {
        return new ProblemDetails
        {
            Status = status,
            Title = GetTitleForStatus(status),
            Detail = detail,
            Extensions = { ["code"] = code }
        };
    }

    private static ValidationProblemDetails CreateValidationProblem(ValidationException ex)
    {
        return new ValidationProblemDetails(ex.Errors)
        {
            Status = StatusCodes.Status400BadRequest,
            Title = "Validation Failed",
            Detail = ex.Message,
            Extensions = { ["code"] = ex.Code }
        };
    }

    private static string GetTitleForStatus(int status) => status switch
    {
        400 => "Bad Request",
        403 => "Forbidden",
        404 => "Not Found",
        409 => "Conflict",
        500 => "Internal Server Error",
        _ => "Error"
    };
}
```

### Register Middleware

```csharp
// API/Program.cs
var builder = WebApplication.CreateBuilder(args);

// ... other services

var app = builder.Build();

// Add exception handling middleware FIRST
app.UseMiddleware<ExceptionHandlingMiddleware>();

app.UseHttpsRedirection();
app.UseAuthorization();
app.MapControllers();

app.Run();
```

---

## 5. Result Pattern (Alternative to Exceptions)

For expected failures, consider the Result pattern instead of exceptions.

### Result Types

```csharp
// Domain/Common/Result.cs
namespace YourApp.Domain.Common;

public class Result
{
    public bool IsSuccess { get; }
    public bool IsFailure => !IsSuccess;
    public Error Error { get; }

    protected Result(bool isSuccess, Error error)
    {
        IsSuccess = isSuccess;
        Error = error;
    }

    public static Result Success() => new(true, Error.None);
    public static Result Failure(Error error) => new(false, error);
    public static Result<T> Success<T>(T value) => new(value, true, Error.None);
    public static Result<T> Failure<T>(Error error) => new(default!, false, error);
}

public class Result<T> : Result
{
    public T Value { get; }

    protected internal Result(T value, bool isSuccess, Error error) 
        : base(isSuccess, error)
    {
        Value = value;
    }

    public static implicit operator Result<T>(T value) => Success(value);
}

public record Error(string Code, string Message)
{
    public static readonly Error None = new(string.Empty, string.Empty);
}
```

### Domain Errors

```csharp
// Domain/Errors/OrderErrors.cs
namespace YourApp.Domain.Errors;

public static class OrderErrors
{
    public static Error NotFound(OrderId id) => 
        new("Order.NotFound", $"Order {id} was not found");

    public static Error NotModifiable(OrderStatus status) => 
        new("Order.NotModifiable", $"Order cannot be modified in {status} status");

    public static Error EmptyOrder => 
        new("Order.Empty", "Order must have at least one item");

    public static Error ProductNotFound(ProductId id) => 
        new("Order.ProductNotFound", $"Product {id} was not found");
}
```

### Using Result in Handlers

```csharp
// Application/Commands/PlaceOrderCommandHandler.cs
public class PlaceOrderCommandHandler 
    : IRequestHandler<PlaceOrderCommand, Result<Guid>>
{
    public async Task<Result<Guid>> Handle(
        PlaceOrderCommand request, 
        CancellationToken ct)
    {
        var order = Order.Create(customerId, address);

        foreach (var itemDto in request.Items)
        {
            var product = await _productRepository.GetByIdAsync(
                ProductId.From(itemDto.ProductId), ct);

            if (product == null)
                return Result.Failure<Guid>(
                    OrderErrors.ProductNotFound(ProductId.From(itemDto.ProductId)));

            var addResult = order.AddItem(product, itemDto.Quantity);
            if (addResult.IsFailure)
                return Result.Failure<Guid>(addResult.Error);
        }

        var placeResult = order.PlaceOrder();
        if (placeResult.IsFailure)
            return Result.Failure<Guid>(placeResult.Error);

        await _orderRepository.AddAsync(order, ct);
        await _unitOfWork.SaveChangesAsync(ct);

        return order.Id.Value;
    }
}
```

### Controller with Result Pattern

```csharp
// API/Controllers/OrdersController.cs
[HttpPost]
public async Task<IActionResult> PlaceOrder(
    PlaceOrderRequest request,
    CancellationToken ct)
{
    var result = await _mediator.Send(new PlaceOrderCommand(...), ct);

    if (result.IsFailure)
    {
        return result.Error.Code switch
        {
            "Order.NotFound" => NotFound(result.Error),
            "Order.ProductNotFound" => NotFound(result.Error),
            _ => BadRequest(result.Error)
        };
    }

    return CreatedAtAction(nameof(GetOrder), new { id = result.Value }, result.Value);
}
```

---

## 6. Error Handling Best Practices

### ✅ DO

```csharp
// ✅ Use specific exception types
throw new OrderNotModifiableException(orderId, status);

// ✅ Include relevant context
throw new NotFoundException("Product", productId.ToString());

// ✅ Validate early (fail fast)
if (request.Items.Count == 0)
    throw new ValidationException("Items", "At least one item required");

// ✅ Let domain exceptions bubble up naturally
order.PlaceOrder(); // May throw DomainException

// ✅ Log with appropriate level
_logger.LogWarning("Order {OrderId} not found", orderId);
_logger.LogError(ex, "Database error while saving order");
```

### ❌ DON'T

```csharp
// ❌ Don't catch and rethrow without adding value
try { ... }
catch (Exception ex) { throw ex; } // Loses stack trace!

// ❌ Don't swallow exceptions
try { ... }
catch (Exception) { } // Silent failure

// ❌ Don't use generic exceptions
throw new Exception("Something went wrong"); // Not specific

// ❌ Don't expose internal details to clients
throw new Exception($"SQL Error: {sqlException.Message}"); // Security risk

// ❌ Don't mix exception types across layers
// Domain shouldn't throw ApplicationException
// Application shouldn't throw DbException
```

---

## Summary

| Layer | Exception Type | HTTP Status | Example |
|-------|---------------|-------------|---------|
| Domain | `DomainException` | 400 | Business rule violated |
| Domain | `InvalidEmailException` | 400 | Value object validation |
| Application | `ValidationException` | 400 | Input validation |
| Application | `NotFoundException` | 404 | Entity not found |
| Application | `ConflictException` | 409 | State conflict |
| Application | `ConcurrencyException` | 409 | Optimistic locking |
| Application | `ForbiddenException` | 403 | Authorization |
| Infrastructure | Caught & wrapped | 500 | Database errors |

**Key Principle:** Each layer handles its own concerns. Domain exceptions are business problems. Application exceptions are use case problems. Infrastructure exceptions should be caught and wrapped.

````
