# C# Error Handling Examples

C#-specific error handling patterns and idioms. For general principles, see [Error Handling Principle](../../principles/06-error-handling.md).

## C# 11+ Guard Clauses

```csharp
// Modern ArgumentNullException
public void ProcessOrder(Order order, Customer customer)
{
    ArgumentNullException.ThrowIfNull(order);
    ArgumentNullException.ThrowIfNull(customer);
    ArgumentException.ThrowIfNullOrWhiteSpace(order.Description);
    
    // Process...
}
```

## Result Pattern with Modern C#

Full implementation using records and pattern matching:

```csharp
// Result type using records (C# 10+)
public abstract record Result<T>
{
    public sealed record Success(T Value) : Result<T>;
    public sealed record Failure(Error Error) : Result<T>;
    
    public bool IsSuccess => this is Success;
    public bool IsFailure => this is Failure;
}

public record Error(string Code, string Message)
{
    public static Error NotFound(string entity, object id) => 
        new("NOT_FOUND", $"{entity} with ID {id} not found");
    
    public static Error Validation(string message) => 
        new("VALIDATION", message);
    
    public static Error Conflict(string message) => 
        new("CONFLICT", message);
}

// Extension methods for fluent API
public static class ResultExtensions
{
    public static TResult Match<T, TResult>(
        this Result<T> result,
        Func<T, TResult> onSuccess,
        Func<Error, TResult> onFailure) => result switch
    {
        Result<T>.Success s => onSuccess(s.Value),
        Result<T>.Failure f => onFailure(f.Error),
        _ => throw new InvalidOperationException()
    };
    
    public static Result<TResult> Map<T, TResult>(
        this Result<T> result,
        Func<T, TResult> mapper) => result switch
    {
        Result<T>.Success s => new Result<TResult>.Success(mapper(s.Value)),
        Result<T>.Failure f => new Result<TResult>.Failure(f.Error),
        _ => throw new InvalidOperationException()
    };
    
    public static async Task<Result<TResult>> MapAsync<T, TResult>(
        this Result<T> result,
        Func<T, Task<TResult>> mapper) => result switch
    {
        Result<T>.Success s => new Result<TResult>.Success(await mapper(s.Value)),
        Result<T>.Failure f => new Result<TResult>.Failure(f.Error),
        _ => throw new InvalidOperationException()
    };
}
```

### Usage in Application Layer

```csharp
public class TransferMoneyHandler(
    IAccountRepository repository,
    IUnitOfWork unitOfWork)
{
    public async Task<Result<TransferResult>> Handle(TransferCommand cmd)
    {
        var fromAccount = await repository.FindAsync(cmd.FromAccountId);
        if (fromAccount is null)
            return new Result<TransferResult>.Failure(
                Error.NotFound("Account", cmd.FromAccountId));
        
        var toAccount = await repository.FindAsync(cmd.ToAccountId);
        if (toAccount is null)
            return new Result<TransferResult>.Failure(
                Error.NotFound("Account", cmd.ToAccountId));
        
        if (fromAccount.Balance < cmd.Amount)
            return new Result<TransferResult>.Failure(
                Error.Validation("Insufficient funds"));
        
        fromAccount.Withdraw(cmd.Amount);
        toAccount.Deposit(cmd.Amount);
        
        await unitOfWork.SaveChangesAsync();
        
        return new Result<TransferResult>.Success(
            new TransferResult(fromAccount.Id, toAccount.Id, cmd.Amount));
    }
}
```

### Minimal API Integration

```csharp
app.MapPost("/api/transfers", async (TransferCommand cmd, TransferMoneyHandler handler) =>
{
    var result = await handler.Handle(cmd);
    
    return result.Match<IResult>(
        onSuccess: r => Results.Ok(r),
        onFailure: e => e.Code switch
        {
            "NOT_FOUND" => Results.NotFound(e.Message),
            "VALIDATION" => Results.BadRequest(e.Message),
            "CONFLICT" => Results.Conflict(e.Message),
            _ => Results.Problem(e.Message)
        });
});
```

## Exception Hierarchy for Domain

```csharp
// Base domain exception
public abstract class DomainException : Exception
{
    public string Code { get; }
    
    protected DomainException(string code, string message) : base(message)
    {
        Code = code;
    }
}

// Specific exceptions with context (C# 12 primary constructors)
public class EntityNotFoundException<TId>(string entityName, TId id)
    : DomainException("NOT_FOUND", $"{entityName} with ID {id} not found")
{
    public string EntityName { get; } = entityName;
    public TId EntityId { get; } = id;
}

public class InsufficientFundsException(decimal available, decimal requested)
    : DomainException("INSUFFICIENT_FUNDS", 
        $"Insufficient funds. Available: {available:C}, Requested: {requested:C}")
{
    public decimal Available { get; } = available;
    public decimal Requested { get; } = requested;
}

public class BusinessRuleViolationException(string rule, string details)
    : DomainException("BUSINESS_RULE", $"Business rule '{rule}' violated: {details}")
{
    public string Rule { get; } = rule;
}
```

## Global Exception Handling (ASP.NET Core 8+)

```csharp
// IExceptionHandler (new in .NET 8)
public class GlobalExceptionHandler(ILogger<GlobalExceptionHandler> logger) 
    : IExceptionHandler
{
    public async ValueTask<bool> TryHandleAsync(
        HttpContext context,
        Exception exception,
        CancellationToken cancellationToken)
    {
        var (statusCode, response) = exception switch
        {
            EntityNotFoundException<Guid> e => (
                StatusCodes.Status404NotFound,
                new ProblemDetails { Title = "Not Found", Detail = e.Message }),
            
            BusinessRuleViolationException e => (
                StatusCodes.Status422UnprocessableEntity,
                new ProblemDetails { Title = "Business Rule Violation", Detail = e.Message }),
            
            ValidationException e => (
                StatusCodes.Status400BadRequest,
                new ProblemDetails { Title = "Validation Error", Detail = e.Message }),
            
            _ => (
                StatusCodes.Status500InternalServerError,
                new ProblemDetails { Title = "Server Error", Detail = "An unexpected error occurred" })
        };
        
        logger.LogError(exception, "Exception handled: {Message}", exception.Message);
        
        context.Response.StatusCode = statusCode;
        await context.Response.WriteAsJsonAsync(response, cancellationToken);
        
        return true;
    }
}

// Registration in Program.cs
builder.Services.AddExceptionHandler<GlobalExceptionHandler>();
app.UseExceptionHandler();
```

## FluentValidation + MediatR Pipeline

```csharp
// Validator
public class CreateOrderValidator : AbstractValidator<CreateOrderCommand>
{
    public CreateOrderValidator()
    {
        RuleFor(x => x.CustomerId)
            .NotEmpty()
            .WithMessage("Customer ID is required");
        
        RuleFor(x => x.Items)
            .NotEmpty()
            .WithMessage("Order must have at least one item");
        
        RuleForEach(x => x.Items).ChildRules(item =>
        {
            item.RuleFor(i => i.Quantity)
                .GreaterThan(0)
                .WithMessage("Quantity must be positive");
        });
    }
}

// MediatR Pipeline Behavior
public class ValidationBehavior<TRequest, TResponse>(
    IEnumerable<IValidator<TRequest>> validators)
    : IPipelineBehavior<TRequest, TResponse>
    where TRequest : IRequest<TResponse>
{
    public async Task<TResponse> Handle(
        TRequest request,
        RequestHandlerDelegate<TResponse> next,
        CancellationToken cancellationToken)
    {
        if (!validators.Any())
            return await next();
        
        var context = new ValidationContext<TRequest>(request);
        
        var failures = validators
            .Select(v => v.Validate(context))
            .SelectMany(r => r.Errors)
            .Where(f => f != null)
            .ToList();
        
        if (failures.Count != 0)
            throw new ValidationException(failures);
        
        return await next();
    }
}
```

## OneOf for Discriminated Unions

Using [OneOf](https://github.com/mcintyre321/OneOf) library:

```csharp
// Install: dotnet add package OneOf

public class GetOrderHandler
{
    public async Task<OneOf<Order, NotFound, Forbidden>> Handle(Guid orderId, Guid userId)
    {
        var order = await _repository.FindAsync(orderId);
        
        if (order is null)
            return new NotFound();
        
        if (order.CustomerId != userId)
            return new Forbidden();
        
        return order;
    }
}

// Usage in controller
public async Task<IActionResult> Get(Guid id)
{
    var result = await _handler.Handle(id, User.GetId());
    
    return result.Match<IActionResult>(
        order => Ok(order),
        notFound => NotFound(),
        forbidden => Forbid());
}
```

---

## Summary: C# Error Handling Checklist

- [ ] Use `ArgumentNullException.ThrowIfNull()` for guard clauses
- [ ] Prefer Result pattern over exceptions for expected failures
- [ ] Use exceptions only for exceptional/unexpected situations
- [ ] Create domain-specific exception hierarchy
- [ ] Implement `IExceptionHandler` in ASP.NET Core 8+
- [ ] Use FluentValidation for input validation
- [ ] Consider OneOf for explicit discriminated unions
- [ ] Log exceptions with structured logging

## Related

- [Error Handling Principle](../../principles/06-error-handling.md) — Theory and concepts
- [Complete Example](08-complete-example.md) — Error handling in context
- [Clean Architecture Error Handling](../../../clean-architecture/examples/csharp/ERROR-HANDLING-EXAMPLES.md) — Architectural patterns
