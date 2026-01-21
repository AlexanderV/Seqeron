# Clean Code + Clean Architecture Integration

> **How Clean Code principles apply within Clean Architecture layers**

This guide shows how to apply Clean Code principles at each architectural layer, demonstrating the synergy between tactical code quality and strategic architecture.

## Overview: Where Clean Code Meets Clean Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    CLEAN ARCHITECTURE LAYERS                     │
├─────────────────────────────────────────────────────────────────┤
│  Presentation    │  Application     │  Domain        │  Infra   │
│  (Controllers)   │  (Use Cases)     │  (Entities)    │  (DB)    │
├─────────────────────────────────────────────────────────────────┤
│                     CLEAN CODE PRINCIPLES                        │
├─────────────────────────────────────────────────────────────────┤
│  Naming ──────────────────────────────────────────────────────► │
│  Small Functions ─────────────────────────────────────────────► │
│  SRP ─────────────────────────────────────────────────────────► │
│  Error Handling ──────────────────────────────────────────────► │
│  Comments/Docs ───────────────────────────────────────────────► │
└─────────────────────────────────────────────────────────────────┘
```

## Layer-Specific Clean Code Guidelines

### Domain Layer

The domain layer benefits most from Clean Code—it represents your business logic.

**Naming: Use Ubiquitous Language**
```csharp
// ❌ Generic names
public class User { }
public class Item { }
public decimal Amount { get; }

// ✅ Domain-specific names (Insurance Domain)
public class PolicyHolder { }
public class InsuranceClaim { }
public Money CoverageLimit { get; }
```

**Rich Domain Models (Not Anemic)**
```csharp
// ❌ Anemic model - logic lives elsewhere
public class Policy
{
    public Guid Id { get; set; }
    public decimal Premium { get; set; }
    public PolicyStatus Status { get; set; }
    public DateTime? CancelledAt { get; set; }
}

// ✅ Rich domain model - behavior with data
public class Policy
{
    public PolicyId Id { get; }
    public Money Premium { get; private set; }
    public PolicyStatus Status { get; private set; }

    public Result Cancel(CancellationReason reason, IClock clock)
    {
        if (Status == PolicyStatus.Cancelled)
            return Result.Failure("Policy already cancelled");

        if (Status == PolicyStatus.Expired)
            return Result.Failure("Cannot cancel expired policy");

        Status = PolicyStatus.Cancelled;
        CancelledAt = clock.UtcNow;
        AddDomainEvent(new PolicyCancelledEvent(Id, reason));

        return Result.Success();
    }

    public Money CalculateRefund(IClock clock)
    {
        if (Status != PolicyStatus.Cancelled)
            return Money.Zero;

        var remainingDays = (EndDate - clock.UtcNow).Days;
        var totalDays = (EndDate - StartDate).Days;
        var refundRatio = (decimal)remainingDays / totalDays;

        return Premium * refundRatio;
    }
}
```

**Value Objects for Domain Concepts**
```csharp
public record PolicyNumber
{
    public string Value { get; }

    public PolicyNumber(string value)
    {
        if (!Regex.IsMatch(value, @"^POL-\d{4}-\d{6}$"))
            throw new InvalidPolicyNumberException(value);
        Value = value;
    }

    public int Year => int.Parse(Value.Substring(4, 4));
}

public record Money
{
    public decimal Amount { get; }
    public Currency Currency { get; }

    public static Money Zero => new(0, Currency.USD);

    public Money(decimal amount, Currency currency)
    {
        if (amount < 0)
            throw new NegativeMoneyException(amount);
        Amount = amount;
        Currency = currency;
    }

    public static Money operator *(Money money, decimal multiplier) =>
        new(money.Amount * multiplier, money.Currency);
}
```

### Application Layer (Use Cases)

**Single Responsibility per Use Case**
```csharp
// ❌ God use case doing too much
public class PolicyService
{
    public async Task CreatePolicy(...) { }
    public async Task CancelPolicy(...) { }
    public async Task RenewPolicy(...) { }
    public async Task ProcessClaim(...) { }
    public async Task GenerateReport(...) { }
}

// ✅ One class per use case (Vertical Slice)
public class CancelPolicyCommand : IRequest<Result>
{
    public required PolicyId PolicyId { get; init; }
    public required CancellationReason Reason { get; init; }
}

public class CancelPolicyHandler : IRequestHandler<CancelPolicyCommand, Result>
{
    private readonly IPolicyRepository _policies;
    private readonly IClock _clock;
    private readonly IUnitOfWork _unitOfWork;

    public CancelPolicyHandler(
        IPolicyRepository policies,
        IClock clock,
        IUnitOfWork unitOfWork)
    {
        _policies = policies;
        _clock = clock;
        _unitOfWork = unitOfWork;
    }

    public async Task<Result> Handle(
        CancelPolicyCommand command,
        CancellationToken cancellationToken)
    {
        var policy = await _policies.GetByIdAsync(command.PolicyId, cancellationToken);

        if (policy is null)
            return Result.Failure($"Policy {command.PolicyId} not found");

        var result = policy.Cancel(command.Reason, _clock);

        if (result.IsFailure)
            return result;

        await _unitOfWork.SaveChangesAsync(cancellationToken);

        return Result.Success();
    }
}
```

**Clean Validation with FluentValidation**
```csharp
public class CancelPolicyValidator : AbstractValidator<CancelPolicyCommand>
{
    public CancelPolicyValidator()
    {
        RuleFor(x => x.PolicyId)
            .NotEmpty()
            .WithMessage("Policy ID is required");

        RuleFor(x => x.Reason)
            .IsInEnum()
            .WithMessage("Valid cancellation reason is required");
    }
}
```

### Infrastructure Layer

**Repository Pattern with Clean Naming**
```csharp
// ❌ Generic repository with unclear methods
public interface IRepository<T>
{
    T Get(Guid id);
    IEnumerable<T> GetAll();
    void Add(T entity);
}

// ✅ Domain-specific repository
public interface IPolicyRepository
{
    Task<Policy?> GetByIdAsync(PolicyId id, CancellationToken ct = default);
    Task<Policy?> GetByNumberAsync(PolicyNumber number, CancellationToken ct = default);
    Task<IReadOnlyList<Policy>> GetActiveByHolderAsync(PolicyHolderId holderId, CancellationToken ct = default);
    Task<IReadOnlyList<Policy>> GetExpiringWithinDaysAsync(int days, CancellationToken ct = default);
    void Add(Policy policy);
    void Remove(Policy policy);
}
```

**Clean EF Core Configurations**
```csharp
public class PolicyConfiguration : IEntityTypeConfiguration<Policy>
{
    public void Configure(EntityTypeBuilder<Policy> builder)
    {
        builder.ToTable("Policies");

        builder.HasKey(p => p.Id);

        builder.Property(p => p.Id)
            .HasConversion(
                id => id.Value,
                value => new PolicyId(value));

        builder.Property(p => p.Number)
            .HasConversion(
                n => n.Value,
                value => new PolicyNumber(value))
            .HasMaxLength(15);

        builder.OwnsOne(p => p.Premium, money =>
        {
            money.Property(m => m.Amount).HasColumnName("PremiumAmount");
            money.Property(m => m.Currency).HasColumnName("PremiumCurrency");
        });

        builder.HasIndex(p => p.Number).IsUnique();
    }
}
```

### Presentation Layer

**Clean Controller Actions**
```csharp
// ❌ Fat controller with business logic
[HttpPost("cancel")]
public async Task<IActionResult> Cancel([FromBody] CancelRequest request)
{
    var policy = await _context.Policies.FindAsync(request.PolicyId);
    if (policy == null) return NotFound();

    if (policy.Status == "Cancelled") return BadRequest("Already cancelled");

    policy.Status = "Cancelled";
    policy.CancelledAt = DateTime.Now;

    await _context.SaveChangesAsync();
    await _emailService.SendCancellationEmail(policy);

    return Ok();
}

// ✅ Thin controller delegating to use case
[HttpPost("{policyId}/cancel")]
[ProducesResponseType(StatusCodes.Status200OK)]
[ProducesResponseType(StatusCodes.Status400BadRequest)]
[ProducesResponseType(StatusCodes.Status404NotFound)]
public async Task<IActionResult> Cancel(
    [FromRoute] Guid policyId,
    [FromBody] CancelPolicyRequest request,
    CancellationToken cancellationToken)
{
    var command = new CancelPolicyCommand
    {
        PolicyId = new PolicyId(policyId),
        Reason = request.Reason
    };

    var result = await _mediator.Send(command, cancellationToken);

    return result.IsSuccess
        ? Ok()
        : BadRequest(result.Error);
}
```

## Cross-Cutting Concerns

### Error Handling Across Layers

```csharp
// Domain exceptions (throw from domain)
public class PolicyNotFoundException : DomainException
{
    public PolicyNotFoundException(PolicyId id)
        : base($"Policy {id} was not found") { }
}

public class PolicyAlreadyCancelledException : DomainException
{
    public PolicyAlreadyCancelledException(PolicyId id)
        : base($"Policy {id} is already cancelled") { }
}

// Global exception handler (Presentation layer)
public class GlobalExceptionHandler : IExceptionHandler
{
    public async ValueTask<bool> TryHandleAsync(
        HttpContext context,
        Exception exception,
        CancellationToken cancellationToken)
    {
        var (statusCode, response) = exception switch
        {
            DomainException e => (400, new ProblemDetails
            {
                Status = 400,
                Title = "Business Rule Violation",
                Detail = e.Message
            }),
            NotFoundException e => (404, new ProblemDetails
            {
                Status = 404,
                Title = "Resource Not Found",
                Detail = e.Message
            }),
            _ => (500, new ProblemDetails
            {
                Status = 500,
                Title = "Internal Server Error",
                Detail = "An unexpected error occurred"
            })
        };

        context.Response.StatusCode = statusCode;
        await context.Response.WriteAsJsonAsync(response, cancellationToken);
        return true;
    }
}
```

### Logging with Meaningful Context

```csharp
public class CancelPolicyHandler : IRequestHandler<CancelPolicyCommand, Result>
{
    private readonly ILogger<CancelPolicyHandler> _logger;

    public async Task<Result> Handle(CancelPolicyCommand command, CancellationToken ct)
    {
        _logger.LogInformation(
            "Cancelling policy {PolicyId} with reason {Reason}",
            command.PolicyId,
            command.Reason);

        var policy = await _policies.GetByIdAsync(command.PolicyId, ct);

        if (policy is null)
        {
            _logger.LogWarning("Policy {PolicyId} not found for cancellation", command.PolicyId);
            return Result.Failure($"Policy {command.PolicyId} not found");
        }

        var result = policy.Cancel(command.Reason, _clock);

        if (result.IsFailure)
        {
            _logger.LogWarning(
                "Failed to cancel policy {PolicyId}: {Error}",
                command.PolicyId,
                result.Error);
            return result;
        }

        await _unitOfWork.SaveChangesAsync(ct);

        _logger.LogInformation(
            "Policy {PolicyId} cancelled successfully. Refund: {RefundAmount}",
            command.PolicyId,
            policy.CalculateRefund(_clock));

        return Result.Success();
    }
}
```

## Summary: Clean Code in Each Layer

| Layer | Key Clean Code Principles |
|-------|--------------------------|
| **Domain** | Ubiquitous language, Rich models, Value objects, No dependencies |
| **Application** | Single responsibility per use case, Clean validation, Result pattern |
| **Infrastructure** | Domain-specific repositories, Clean configurations, Adapter pattern |
| **Presentation** | Thin controllers, Clean DTOs, Global exception handling |

## Related Resources

- [Clean Architecture Skill](../../../clean-architecture/SKILL.md)
- [Clean Architecture Checklist](../../../clean-architecture/CHECKLIST.md)
- [Complete Banking Example](08-complete-example.md)
- [Refactoring Journey](09-refactoring-journey.md)
