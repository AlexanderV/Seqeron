# Principle 6: Error Handling

> "Error handling is important, but if it obscures logic, it's wrong."
> — Robert C. Martin

## Overview

Error handling is one of the things that programmers do to ensure that code works properly. Things can go wrong, and when they do, we as programmers are responsible for making sure that our code does what it needs to do.

## Key Rules

### 6.1 Use Exceptions Rather Than Return Codes

**❌ BAD - Error Codes:**
```csharp
public int SaveCustomer(Customer customer)
{
    if (customer == null)
        return ERROR_NULL_CUSTOMER;

    if (!IsValid(customer))
        return ERROR_INVALID_CUSTOMER;

    if (!_db.IsConnected)
        return ERROR_DATABASE_UNAVAILABLE;

    // ...
    return SUCCESS;
}

// Caller must check every call
int result = SaveCustomer(customer);
if (result == SUCCESS)
{
    result = SendConfirmation(customer);
    if (result == SUCCESS)
    {
        // ...
    }
}
```

**✅ GOOD - Exceptions:**
```csharp
public void SaveCustomer(Customer customer)
{
    ArgumentNullException.ThrowIfNull(customer);

    if (!IsValid(customer))
        throw new ValidationException("Customer is not valid");

    // Save logic...
}

// Caller focuses on happy path
try
{
    SaveCustomer(customer);
    SendConfirmation(customer);
    UpdateStatistics();
}
catch (ValidationException ex)
{
    HandleValidationError(ex);
}
catch (DatabaseException ex)
{
    HandleDatabaseError(ex);
}
```

### 6.2 Write Your Try-Catch-Finally Statement First

When writing code that might throw exceptions, start with try-catch-finally first to define the scope.

```csharp
public List<Customer> GetActiveCustomers()
{
    try
    {
        return FetchActiveCustomersFromDatabase();
    }
    catch (DatabaseException ex)
    {
        _logger.LogError(ex, "Failed to fetch active customers");
        throw new CustomerServiceException("Unable to retrieve customers", ex);
    }
    finally
    {
        _metrics.RecordDatabaseCall();
    }
}
```

### 6.3 Provide Context with Exceptions

Create informative error messages that help debugging.

**❌ BAD:**
```csharp
throw new Exception("Error");
throw new Exception("Not found");
throw new InvalidOperationException("Invalid");
```

**✅ GOOD:**
```csharp
throw new CustomerNotFoundException(
    $"Customer with ID {customerId} was not found in the database");

throw new OrderValidationException(
    $"Order {orderId} cannot be processed: {string.Join(", ", errors)}");

throw new PaymentFailedException(
    $"Payment of {amount:C} for order {orderId} failed with gateway response: {gatewayResponse}");
```

### 6.4 Define Exception Classes in Terms of Caller's Needs

Create domain-specific exceptions and wrap third-party exceptions.

**❌ BAD - Multiple catch blocks for third-party exceptions:**
```csharp
try
{
    _paymentGateway.Charge(card, amount);
}
catch (GatewayTimeoutException ex) { /* handle */ }
catch (GatewayAuthException ex) { /* handle */ }
catch (GatewayNetworkException ex) { /* handle */ }
catch (GatewayValidationException ex) { /* handle */ }
```

**✅ GOOD - Wrap in adapter with domain exception:**
```csharp
public class PaymentGatewayAdapter : IPaymentGateway
{
    private readonly ThirdPartyPaymentGateway _gateway;

    public void Charge(CreditCard card, Money amount)
    {
        try
        {
            _gateway.Process(card.Number, amount.Value);
        }
        catch (Exception ex)
        {
            throw new PaymentProcessingException(
                $"Payment failed for amount {amount}", ex);
        }
    }
}

// Caller - single exception to handle
try
{
    _paymentGateway.Charge(card, amount);
}
catch (PaymentProcessingException ex)
{
    HandlePaymentFailure(ex);
}
```

### 6.5 Define a Domain Exception Hierarchy

```csharp
// Base domain exception
public abstract class DomainException : Exception
{
    protected DomainException(string message) : base(message) { }
    protected DomainException(string message, Exception inner) : base(message, inner) { }
}

// Entity not found
public class EntityNotFoundException : DomainException
{
    public string EntityType { get; }
    public object EntityId { get; }

    public EntityNotFoundException(string entityType, object entityId)
        : base($"{entityType} with ID {entityId} was not found")
    {
        EntityType = entityType;
        EntityId = entityId;
    }
}

public class CustomerNotFoundException : EntityNotFoundException
{
    public CustomerNotFoundException(int customerId)
        : base("Customer", customerId) { }
}

// Validation errors
public class ValidationException : DomainException
{
    public IReadOnlyList<ValidationError> Errors { get; }

    public ValidationException(IEnumerable<ValidationError> errors)
        : base("Validation failed")
    {
        Errors = errors.ToList();
    }
}

// Business rule violations
public class BusinessRuleViolationException : DomainException
{
    public string RuleName { get; }

    public BusinessRuleViolationException(string ruleName, string message)
        : base(message)
    {
        RuleName = ruleName;
    }
}
```

### 6.6 Don't Return Null

Returning null leads to NullReferenceException throughout the codebase.

**❌ BAD:**
```csharp
public Customer GetCustomer(int id)
{
    return _customers.Find(id); // Returns null if not found
}

// Caller forgets to check
var customer = GetCustomer(123);
var name = customer.Name; // NullReferenceException!
```

**✅ GOOD - Option 1: Throw Exception:**
```csharp
public Customer GetCustomer(int id)
{
    var customer = _customers.Find(id);
    if (customer is null)
        throw new CustomerNotFoundException(id);
    return customer;
}
```

**✅ GOOD - Option 2: Null Object Pattern:**
```csharp
public Customer GetCustomer(int id)
{
    var customer = _customers.Find(id);
    return customer ?? Customer.Unknown;
}

public class Customer
{
    public static readonly Customer Unknown = new UnknownCustomer();

    public virtual bool IsUnknown => false;
    public virtual string Name { get; }
}

private class UnknownCustomer : Customer
{
    public override bool IsUnknown => true;
    public override string Name => "Unknown Customer";
}
```

**✅ GOOD - Option 3: Optional/Maybe Pattern:**
```csharp
public Option<Customer> FindCustomer(int id)
{
    var customer = _customers.Find(id);
    return customer is not null
        ? Option.Some(customer)
        : Option.None<Customer>();
}

// Usage - forces caller to handle both cases
var result = FindCustomer(123);
result.Match(
    some: customer => ProcessOrder(customer),
    none: () => NotifyCustomerNotFound()
);
```

**✅ GOOD - Option 4: Result Pattern:**
```csharp
public Result<Customer> GetCustomer(int id)
{
    var customer = _customers.Find(id);
    return customer is not null
        ? Result.Success(customer)
        : Result.Failure<Customer>($"Customer {id} not found");
}

// Usage
var result = GetCustomer(123);
if (result.IsSuccess)
{
    ProcessOrder(result.Value);
}
else
{
    _logger.LogWarning(result.Error);
}
```

### 6.7 Don't Pass Null

Passing null as a parameter is even worse than returning it.

**❌ BAD:**
```csharp
public void ProcessOrder(Customer customer, Order order)
{
    // What if customer or order is null?
    // Defensive programming everywhere!
    if (customer == null || order == null)
        return; // Silent failure
}
```

**✅ GOOD:**
```csharp
public void ProcessOrder(Customer customer, Order order)
{
    ArgumentNullException.ThrowIfNull(customer);
    ArgumentNullException.ThrowIfNull(order);

    // Now we can trust the parameters
}
```

### 6.8 Use Nullable Reference Types (C# 8+)

Enable nullable reference types for compile-time null checking.

```csharp
#nullable enable

public class OrderService
{
    // Non-nullable - guaranteed not null
    private readonly IOrderRepository _repository;

    // Nullable - explicitly may be null
    public Customer? FindCustomer(int id)
    {
        return _customers.FirstOrDefault(c => c.Id == id);
    }

    // Non-nullable return - must not return null
    public Customer GetCustomer(int id)
    {
        return FindCustomer(id)
            ?? throw new CustomerNotFoundException(id);
    }
}
```

### 6.9 Extract Error Handling

Don't let error handling obscure business logic.

**❌ BAD:**
```csharp
public void ProcessOrder(Order order)
{
    try
    {
        var customer = _customerRepository.GetById(order.CustomerId);
        try
        {
            ValidateOrder(order, customer);
            try
            {
                _paymentService.Charge(order);
                try
                {
                    _inventoryService.Reserve(order);
                    _orderRepository.Save(order);
                }
                catch (InventoryException ex)
                {
                    _paymentService.Refund(order);
                    throw;
                }
            }
            catch (PaymentException ex)
            {
                _logger.LogError(ex, "Payment failed");
                throw;
            }
        }
        catch (ValidationException ex)
        {
            _logger.LogWarning(ex, "Order validation failed");
            throw;
        }
    }
    catch (CustomerNotFoundException ex)
    {
        _logger.LogError(ex, "Customer not found");
        throw;
    }
}
```

**✅ GOOD:**
```csharp
public void ProcessOrder(Order order)
{
    var customer = GetCustomerOrThrow(order.CustomerId);
    ValidateOrderOrThrow(order, customer);
    ChargeAndReserve(order);
    _orderRepository.Save(order);
}

private Customer GetCustomerOrThrow(int customerId)
{
    try
    {
        return _customerRepository.GetById(customerId);
    }
    catch (CustomerNotFoundException ex)
    {
        _logger.LogError(ex, "Customer {Id} not found", customerId);
        throw;
    }
}

private void ChargeAndReserve(Order order)
{
    _paymentService.Charge(order);
    try
    {
        _inventoryService.Reserve(order);
    }
    catch
    {
        _paymentService.Refund(order);
        throw;
    }
}
```

## Quick Checklist

- [ ] Uses exceptions, not error codes
- [ ] Specific exception types (not generic Exception)
- [ ] Exceptions provide context
- [ ] No returning null - use exceptions, Null Object, or Option/Result
- [ ] No passing null as parameters
- [ ] Nullable reference types enabled
- [ ] Third-party exceptions wrapped in domain exceptions
- [ ] Try-catch blocks don't obscure logic
- [ ] Domain exception hierarchy defined

## See Also

- [C# Error Handling Examples](../examples/csharp/03-error-handling.md) - Detailed C# examples
- [Functions](02-functions.md) - Extract Try/Catch blocks
- [Boundaries](07-boundaries.md) - Wrapping third-party exceptions
