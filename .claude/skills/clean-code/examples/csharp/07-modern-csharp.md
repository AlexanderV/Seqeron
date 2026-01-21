# Modern C# Practices for Clean Code

> Modern C# features enable writing cleaner, safer, and more expressive code.

## Contents

1. [Records](#records)
2. [Nullable Reference Types (NRT)](#nullable-reference-types)
3. [Pattern Matching](#pattern-matching)
4. [Primary Constructors](#primary-constructors)
5. [Required Properties](#required-properties)
6. [Global Usings & File-Scoped Namespaces](#global-usings--file-scoped-namespaces)
7. [Collection Expressions](#collection-expressions)
8. [Raw String Literals](#raw-string-literals)

---

## Records

Records provide immutability and value equality out of the box.

### ❌ BAD — Traditional Approach

```csharp
public class CustomerDto
{
    public string FirstName { get; set; }
    public string LastName { get; set; }
    public string Email { get; set; }
    
    // Must write manually
    public override bool Equals(object? obj)
    {
        if (obj is not CustomerDto other) return false;
        return FirstName == other.FirstName 
            && LastName == other.LastName 
            && Email == other.Email;
    }
    
    public override int GetHashCode()
    {
        return HashCode.Combine(FirstName, LastName, Email);
    }
    
    public override string ToString()
    {
        return $"CustomerDto {{ FirstName = {FirstName}, LastName = {LastName}, Email = {Email} }}";
    }
}
```

### ✅ GOOD — Records

```csharp
// Value equality, immutability, ToString — all built-in
public record CustomerDto(string FirstName, string LastName, string Email);

// Mutation via with-expression
var updated = original with { Email = "new@email.com" };
```

### ✅ Record Struct for High-Performance Scenarios

```csharp
// No heap allocation
public readonly record struct Point(int X, int Y);

public readonly record struct Money(decimal Amount, string Currency)
{
    public Money Add(Money other)
    {
        if (Currency != other.Currency)
            throw new InvalidOperationException("Cannot add different currencies");
        return this with { Amount = Amount + other.Amount };
    }
}
```

### ✅ Records for Domain Events

```csharp
public abstract record DomainEvent
{
    public Guid EventId { get; } = Guid.NewGuid();
    public DateTimeOffset OccurredAt { get; } = DateTimeOffset.UtcNow;
}

public record OrderCreated(Guid OrderId, Guid CustomerId, decimal TotalAmount) : DomainEvent;
public record OrderShipped(Guid OrderId, string TrackingNumber, DateTimeOffset ShippedAt) : DomainEvent;
public record OrderCancelled(Guid OrderId, string Reason) : DomainEvent;
```

---

## Nullable Reference Types

NRT helps prevent NullReferenceException at compile time.

### Enabling NRT

```xml
<!-- In .csproj -->
<PropertyGroup>
    <Nullable>enable</Nullable>
</PropertyGroup>
```

### ❌ BAD — Ignoring Warnings

```csharp
public class UserService
{
    private readonly IUserRepository _repository;
    
    // Warning: possible null reference
    public User GetUser(int id)
    {
        return _repository.FindById(id); // May return null!
    }
    
    public string GetUserEmail(int id)
    {
        var user = GetUser(id);
        return user.Email; // NullReferenceException!
    }
}
```

### ✅ GOOD — Explicit Nullable Types

```csharp
public class UserService
{
    private readonly IUserRepository _repository;
    
    public User? FindUser(int id)
    {
        return _repository.FindById(id);
    }
    
    // Returns non-null or throws exception
    public User GetUser(int id)
    {
        return _repository.FindById(id) 
            ?? throw new NotFoundException($"User {id} not found");
    }
    
    public string? GetUserEmail(int id)
    {
        var user = FindUser(id);
        return user?.Email;
    }
}
```

### ✅ Using Attributes

```csharp
using System.Diagnostics.CodeAnalysis;

public class CacheService
{
    private readonly Dictionary<string, object> _cache = new();
    
    // Compiler understands that on true, result is not null
    public bool TryGet<T>(string key, [NotNullWhen(true)] out T? value)
    {
        if (_cache.TryGetValue(key, out var obj) && obj is T typed)
        {
            value = typed;
            return true;
        }
        value = default;
        return false;
    }
    
    // Return value is never null
    [return: NotNull]
    public T GetOrCreate<T>(string key, Func<T> factory) where T : notnull
    {
        if (!_cache.TryGetValue(key, out var obj) || obj is not T typed)
        {
            typed = factory();
            _cache[key] = typed;
        }
        return typed;
    }
    
    // After calling, parameter is guaranteed non-null
    public void EnsureNotNull([NotNull] string? parameter)
    {
        if (string.IsNullOrEmpty(parameter))
            throw new ArgumentException("Parameter required");
    }
}
```

---

## Pattern Matching

Pattern matching makes code more declarative and readable.

### ❌ BAD — Traditional Switch

```csharp
public decimal CalculateDiscount(Customer customer)
{
    if (customer == null)
        return 0;
    
    if (customer.IsPremium)
    {
        if (customer.YearsActive > 5)
            return 0.25m;
        else
            return 0.15m;
    }
    else if (customer.YearsActive > 3)
    {
        return 0.10m;
    }
    else
    {
        return 0.05m;
    }
}
```

### ✅ GOOD — Switch Expression

```csharp
public decimal CalculateDiscount(Customer? customer) => customer switch
{
    null => 0,
    { IsPremium: true, YearsActive: > 5 } => 0.25m,
    { IsPremium: true } => 0.15m,
    { YearsActive: > 3 } => 0.10m,
    _ => 0.05m
};
```

### ✅ Pattern Matching for Type Checks

```csharp
public string FormatPayment(IPayment payment) => payment switch
{
    CreditCardPayment { IsExpired: true } => "Credit card expired",
    CreditCardPayment card => $"Card ending in {card.LastFourDigits}",
    BankTransferPayment transfer => $"Transfer from {transfer.BankName}",
    CryptoPayment { Currency: "BTC" } crypto => $"{crypto.Amount} Bitcoin",
    CryptoPayment crypto => $"{crypto.Amount} {crypto.Currency}",
    null => throw new ArgumentNullException(nameof(payment)),
    _ => throw new NotSupportedException($"Unknown payment type: {payment.GetType().Name}")
};
```

### ✅ List Patterns (C# 11+)

```csharp
public string AnalyzeArgs(string[] args) => args switch
{
    [] => "No arguments provided",
    ["help" or "-h" or "--help"] => "Showing help...",
    ["version" or "-v"] => "Version 1.0.0",
    ["run", var file] => $"Running {file}",
    ["run", var file, .. var rest] => $"Running {file} with {rest.Length} additional args",
    [var single] => $"Single argument: {single}",
    [_, _, ..] => $"Multiple arguments: {args.Length}"
};
```

### ✅ Relational Patterns

```csharp
public string GetTemperatureDescription(double celsius) => celsius switch
{
    < -20 => "Extreme cold",
    >= -20 and < 0 => "Freezing",
    >= 0 and < 10 => "Cold",
    >= 10 and < 20 => "Cool",
    >= 20 and < 30 => "Comfortable",
    >= 30 and < 40 => "Hot",
    >= 40 => "Extreme heat"
};
```

---

## Primary Constructors

Primary constructors (C# 12) reduce boilerplate for DI.

### ❌ BAD — Traditional Style

```csharp
public class OrderService
{
    private readonly IOrderRepository _orderRepository;
    private readonly IPaymentService _paymentService;
    private readonly INotificationService _notificationService;
    private readonly ILogger<OrderService> _logger;
    
    public OrderService(
        IOrderRepository orderRepository,
        IPaymentService paymentService,
        INotificationService notificationService,
        ILogger<OrderService> logger)
    {
        _orderRepository = orderRepository;
        _paymentService = paymentService;
        _notificationService = notificationService;
        _logger = logger;
    }
    
    public async Task ProcessOrder(Order order)
    {
        _logger.LogInformation("Processing order {OrderId}", order.Id);
        await _paymentService.Charge(order);
        await _orderRepository.Save(order);
        await _notificationService.SendConfirmation(order);
    }
}
```

### ✅ GOOD — Primary Constructors

```csharp
public class OrderService(
    IOrderRepository orderRepository,
    IPaymentService paymentService,
    INotificationService notificationService,
    ILogger<OrderService> logger)
{
    public async Task ProcessOrder(Order order)
    {
        logger.LogInformation("Processing order {OrderId}", order.Id);
        await paymentService.Charge(order);
        await orderRepository.Save(order);
        await notificationService.SendConfirmation(order);
    }
}
```

### ✅ Primary Constructors with Validation

```csharp
public class EmailAddress(string value)
{
    public string Value { get; } = IsValid(value) 
        ? value 
        : throw new ArgumentException($"Invalid email: {value}");
    
    private static bool IsValid(string email) => 
        !string.IsNullOrWhiteSpace(email) && email.Contains('@');
    
    public override string ToString() => Value;
}
```

---

## Required Properties

Required modifier (C# 11) guarantees initialization.

### ❌ BAD — Nullable with Hope

```csharp
public class CreateUserRequest
{
    public string? Name { get; set; }  // Can forget
    public string? Email { get; set; } // Can forget
}

// Usage
var request = new CreateUserRequest { Name = "John" }; // Email forgotten!
```

### ✅ GOOD — Required Properties

```csharp
public class CreateUserRequest
{
    public required string Name { get; init; }
    public required string Email { get; init; }
    public string? Phone { get; init; } // Optional
}

// Compiler requires all required properties
var request = new CreateUserRequest 
{ 
    Name = "John", 
    Email = "john@example.com" 
};
```

### ✅ Combining with Records

```csharp
public record CreateOrderRequest
{
    public required Guid CustomerId { get; init; }
    public required List<OrderItemRequest> Items { get; init; }
    public string? Notes { get; init; }
    public ShippingAddress? ShippingAddress { get; init; }
}

public record OrderItemRequest
{
    public required Guid ProductId { get; init; }
    public required int Quantity { get; init; }
}
```

---

## Global Usings & File-Scoped Namespaces

Reducing visual noise in files.

### ❌ BAD — Repeating Usings

```csharp
// In every file:
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Microsoft.Extensions.Logging;
using Microsoft.EntityFrameworkCore;

namespace MyCompany.MyProject.Domain.Entities
{
    public class Order
    {
        // ...
    }
}
```

### ✅ GOOD — GlobalUsings.cs

```csharp
// GlobalUsings.cs
global using System;
global using System.Collections.Generic;
global using System.Linq;
global using System.Threading;
global using System.Threading.Tasks;
global using Microsoft.Extensions.Logging;

// Or in .csproj
// <ItemGroup>
//   <Using Include="System.Collections.Generic" />
//   <Using Include="Microsoft.Extensions.Logging" />
// </ItemGroup>
```

### ✅ File-Scoped Namespace

```csharp
// One line instead of braces
namespace MyCompany.MyProject.Domain.Entities;

public class Order
{
    public Guid Id { get; private set; }
    public required Guid CustomerId { get; init; }
    public List<OrderLine> Lines { get; } = [];
}

public class OrderLine
{
    public required Guid ProductId { get; init; }
    public required int Quantity { get; init; }
}
```

---

## Collection Expressions

Simplified collection syntax (C# 12).

### ❌ BAD — Old Syntax

```csharp
int[] numbers = new int[] { 1, 2, 3, 4, 5 };
List<string> names = new List<string> { "Alice", "Bob" };
Dictionary<string, int> scores = new Dictionary<string, int>
{
    { "Alice", 100 },
    { "Bob", 85 }
};
```

### ✅ GOOD — Collection Expressions

```csharp
int[] numbers = [1, 2, 3, 4, 5];
List<string> names = ["Alice", "Bob"];

// Spread operator
int[] first = [1, 2, 3];
int[] second = [4, 5, 6];
int[] combined = [..first, ..second]; // [1, 2, 3, 4, 5, 6]

// In method arguments
DoSomething([1, 2, 3]);
```

### ✅ Collection Expressions for Immutable Collections

```csharp
using System.Collections.Immutable;

ImmutableArray<int> immutableArray = [1, 2, 3];
ImmutableList<string> immutableList = ["a", "b", "c"];

// Empty collection
List<Order> emptyOrders = [];
```

---

## Raw String Literals

Multi-line strings without escape sequences (C# 11).

### ❌ BAD — Escape Hell

```csharp
string json = "{\n  \"name\": \"John\",\n  \"email\": \"john@example.com\"\n}";

string sql = @"
SELECT *
FROM Users
WHERE Email = @Email
    AND Status = 'Active'";
```

### ✅ GOOD — Raw String Literals

```csharp
string json = """
    {
        "name": "John",
        "email": "john@example.com"
    }
    """;

string sql = """
    SELECT *
    FROM Users
    WHERE Email = @Email
        AND Status = 'Active'
    """;
```

### ✅ Raw Strings with Interpolation

```csharp
string name = "John";
string email = "john@example.com";

// Use $$ to indicate that {} is interpolation
// (number of $ matches number of { for interpolation)
string json = $$"""
    {
        "name": "{{name}}",
        "email": "{{email}}"
    }
    """;
```

---

## Summary: Modern C# Checklist

### Records
- [ ] DTOs use `record` instead of `class`
- [ ] Domain Events are records
- [ ] High-performance scenarios use `readonly record struct`

### Nullable Reference Types
- [ ] `<Nullable>enable</Nullable>` in project
- [ ] All nullable parameters marked with `?`
- [ ] Using `[NotNull]`, `[NotNullWhen]` attributes

### Pattern Matching
- [ ] Switch expressions instead of if-else chains
- [ ] Property patterns for complex conditions
- [ ] List patterns for collections

### Primary Constructors
- [ ] Services use primary constructors
- [ ] No duplicate field assignments

### Required Properties
- [ ] Required properties instead of constructors for DTOs
- [ ] init-only for immutability

### Syntax Improvements
- [ ] File-scoped namespaces
- [ ] Collection expressions `[]`
- [ ] Raw string literals for JSON/SQL

---

## Related

- [Meaningful Names](../../principles/01-meaningful-names.md) — Naming conventions
- [Objects and Data Structures](../../principles/05-objects-and-data-structures.md) — Records for DTOs
- [Error Handling](../../principles/06-error-handling.md) — NRT for null safety
- [Concurrency](../../principles/11-concurrency.md) — Async patterns with modern C#
