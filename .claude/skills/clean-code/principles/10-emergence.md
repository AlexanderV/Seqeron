# Principle 10: Emergence

> "Following the practice of simple design can and does encourage and enable developers to adhere to good principles and patterns."
> — Robert C. Martin

## Overview

Kent Beck's four rules of Simple Design help create well-designed software through emergent behavior. A design is "simple" if it follows these rules (in order of importance).

## The Four Rules of Simple Design

### Rule 1: Runs All the Tests

A system that cannot be verified should not be deployed. Tests are the primary driver of design quality.

**Why it matters:**
- Testable systems tend to have small, single-purpose classes
- Tight coupling makes testing difficult → forces better design
- Tests provide safety net for refactoring

```csharp
// Testable code = better design
public class OrderProcessor
{
    private readonly IOrderRepository _repository;
    private readonly IPaymentGateway _payment;
    private readonly INotificationService _notifications;

    public OrderProcessor(
        IOrderRepository repository,
        IPaymentGateway payment,
        INotificationService notifications)
    {
        _repository = repository;
        _payment = payment;
        _notifications = notifications;
    }

    public async Task ProcessAsync(Order order)
    {
        await _payment.ChargeAsync(order.Total);
        await _repository.SaveAsync(order);
        await _notifications.SendConfirmationAsync(order);
    }
}

// All dependencies are injectable → easy to test
[Fact]
public async Task Process_ShouldChargePayment()
{
    var mockPayment = new Mock<IPaymentGateway>();
    var processor = new OrderProcessor(
        Mock.Of<IOrderRepository>(),
        mockPayment.Object,
        Mock.Of<INotificationService>());

    await processor.ProcessAsync(new Order { Total = 100m });

    mockPayment.Verify(p => p.ChargeAsync(100m), Times.Once);
}
```

### Rule 2: Contains No Duplication (DRY)

Duplication is the primary enemy of a well-designed system.

**Types of duplication:**
1. **Exact duplication** - Copy-pasted code
2. **Structural duplication** - Similar algorithms
3. **Conceptual duplication** - Same concept implemented differently

**❌ BAD - Duplicated validation:**
```csharp
public class UserService
{
    public void CreateUser(User user)
    {
        if (string.IsNullOrWhiteSpace(user.Name))
            throw new ValidationException("Name is required");
        if (string.IsNullOrWhiteSpace(user.Email))
            throw new ValidationException("Email is required");
        if (!IsValidEmail(user.Email))
            throw new ValidationException("Email is invalid");
        // Create...
    }

    public void UpdateUser(User user)
    {
        if (string.IsNullOrWhiteSpace(user.Name))
            throw new ValidationException("Name is required");
        if (string.IsNullOrWhiteSpace(user.Email))
            throw new ValidationException("Email is required");
        if (!IsValidEmail(user.Email))
            throw new ValidationException("Email is invalid");
        // Update...
    }
}
```

**✅ GOOD - Extracted validation:**
```csharp
public class UserValidator
{
    public void Validate(User user)
    {
        ValidateRequired(user.Name, nameof(user.Name));
        ValidateRequired(user.Email, nameof(user.Email));
        ValidateEmail(user.Email);
    }

    private void ValidateRequired(string value, string fieldName)
    {
        if (string.IsNullOrWhiteSpace(value))
            throw new ValidationException($"{fieldName} is required");
    }

    private void ValidateEmail(string email)
    {
        if (!Regex.IsMatch(email, @"^[^@]+@[^@]+\.[^@]+$"))
            throw new ValidationException("Email is invalid");
    }
}

public class UserService
{
    private readonly UserValidator _validator;

    public void CreateUser(User user)
    {
        _validator.Validate(user);
        // Create...
    }

    public void UpdateUser(User user)
    {
        _validator.Validate(user);
        // Update...
    }
}
```

**Template Method for structural duplication:**
```csharp
// Similar algorithms with slight variations
public abstract class ReportGenerator
{
    public Report Generate(ReportData data)
    {
        var report = new Report();
        report.Header = CreateHeader(data);
        report.Body = CreateBody(data);    // Varies
        report.Footer = CreateFooter(data);
        return report;
    }

    protected virtual string CreateHeader(ReportData data)
        => $"Report: {data.Title}";

    protected abstract string CreateBody(ReportData data);

    protected virtual string CreateFooter(ReportData data)
        => $"Generated: {DateTime.Now}";
}

public class SalesReportGenerator : ReportGenerator
{
    protected override string CreateBody(ReportData data)
    {
        return FormatSalesData(data.Sales);
    }
}

public class InventoryReportGenerator : ReportGenerator
{
    protected override string CreateBody(ReportData data)
    {
        return FormatInventoryData(data.Inventory);
    }
}
```

### Rule 3: Expresses the Intent of the Programmer

Code should clearly communicate what it does and why.

**Techniques for expressiveness:**

**1. Good Names:**
```csharp
// Clear intent from names alone
public class VipDiscountPolicy : IDiscountPolicy
{
    public decimal Calculate(Customer customer, Order order)
    {
        if (customer.IsVip && order.IsEligibleForVipDiscount)
        {
            return order.Subtotal * VipDiscountRate;
        }
        return 0m;
    }
}
```

**2. Small Functions:**
```csharp
public async Task ProcessOrderAsync(Order order)
{
    await ValidateOrder(order);
    await ReserveInventory(order);
    await ChargePayment(order);
    await SendConfirmation(order);
}
```

**3. Standard Patterns:**
```csharp
// Using well-known patterns makes intent clear
public class OrderFactory  // Everyone knows what a Factory does
{
    public Order CreateFromCart(ShoppingCart cart)
    {
        // ...
    }
}

public class CustomerRepository  // Repository pattern
{
    public Task<Customer> GetByIdAsync(CustomerId id);
}

public class PriceCalculatorStrategy  // Strategy pattern
{
    public decimal Calculate(Order order);
}
```

**4. Well-Written Tests:**
```csharp
// Tests document behavior
public class DiscountCalculatorTests
{
    [Fact]
    public void Calculate_ForVipCustomer_Returns10PercentDiscount()
    {
        var customer = new Customer { IsVip = true };
        var order = new Order { Subtotal = 100m };
        var calculator = new DiscountCalculator();

        var discount = calculator.Calculate(customer, order);

        Assert.Equal(10m, discount);
    }

    [Fact]
    public void Calculate_ForRegularCustomer_ReturnsNoDiscount()
    {
        var customer = new Customer { IsVip = false };
        var order = new Order { Subtotal = 100m };
        var calculator = new DiscountCalculator();

        var discount = calculator.Calculate(customer, order);

        Assert.Equal(0m, discount);
    }
}
```

### Rule 4: Minimizes the Number of Classes and Methods

While following the other rules, keep the count of classes and methods low.

**Balance is key:**
- Don't create classes for every tiny concept
- Don't create methods for one-liners
- Avoid dogmatic adherence to rules

**❌ Over-engineered:**
```csharp
// Too many tiny classes
public interface INameValidator { }
public class NameValidator : INameValidator { }
public interface IEmailValidator { }
public class EmailValidator : IEmailValidator { }
public interface IPhoneValidator { }
public class PhoneValidator : IPhoneValidator { }
public interface IAddressValidator { }
public class AddressValidator : IAddressValidator { }
public class UserValidatorFactory { }
public class UserValidatorContext { }
```

**✅ Balanced:**
```csharp
public class UserValidator
{
    public ValidationResult Validate(User user)
    {
        var errors = new List<ValidationError>();

        if (!IsValidName(user.Name))
            errors.Add(new ValidationError("Name", "Invalid name"));

        if (!IsValidEmail(user.Email))
            errors.Add(new ValidationError("Email", "Invalid email"));

        return new ValidationResult(errors);
    }

    private bool IsValidName(string name) => !string.IsNullOrWhiteSpace(name);
    private bool IsValidEmail(string email) => email.Contains("@");
}
```

## Applying the Rules

### Priority Order

1. **First:** Make it work (all tests pass)
2. **Second:** Remove duplication
3. **Third:** Make it expressive
4. **Fourth:** Minimize elements

### Refactoring Workflow

```
Write Test → Make it Pass → Refactor
    ↑                            │
    └────────────────────────────┘
```

**During refactoring:**
1. Does it still pass all tests? (Rule 1)
2. Is there any duplication? (Rule 2)
3. Is the intent clear? (Rule 3)
4. Are there unnecessary elements? (Rule 4)

## Simple Design in Practice

### Before: Complex, Duplicated Code
```csharp
public class OrderService
{
    public decimal CalculateTotalWithTax(Order order)
    {
        decimal total = 0;
        foreach (var item in order.Items)
        {
            total += item.Price * item.Quantity;
        }

        if (order.Customer.IsVip)
        {
            total = total * 0.9m; // 10% discount
        }

        if (order.ShippingAddress.State == "NY")
        {
            total = total * 1.08m; // NY tax
        }
        else if (order.ShippingAddress.State == "CA")
        {
            total = total * 1.0725m; // CA tax
        }

        return total;
    }

    public decimal CalculateTotalWithTaxAndShipping(Order order)
    {
        decimal total = 0;
        foreach (var item in order.Items)
        {
            total += item.Price * item.Quantity;
        }

        if (order.Customer.IsVip)
        {
            total = total * 0.9m; // 10% discount - DUPLICATED!
        }

        // Same tax logic DUPLICATED!
        if (order.ShippingAddress.State == "NY")
        {
            total = total * 1.08m;
        }
        else if (order.ShippingAddress.State == "CA")
        {
            total = total * 1.0725m;
        }

        total += CalculateShipping(order);
        return total;
    }
}
```

### After: Simple Design Applied
```csharp
// Rule 3: Expressive - clear interfaces
public interface IDiscountPolicy
{
    decimal Apply(decimal amount, Customer customer);
}

public interface ITaxCalculator
{
    decimal Calculate(decimal amount, Address address);
}

// Rule 2: No duplication - single responsibility
public class VipDiscountPolicy : IDiscountPolicy
{
    private const decimal VipDiscountRate = 0.10m;

    public decimal Apply(decimal amount, Customer customer)
    {
        return customer.IsVip ? amount * (1 - VipDiscountRate) : amount;
    }
}

public class StateTaxCalculator : ITaxCalculator
{
    private static readonly Dictionary<string, decimal> TaxRates = new()
    {
        ["NY"] = 0.08m,
        ["CA"] = 0.0725m
    };

    public decimal Calculate(decimal amount, Address address)
    {
        var rate = TaxRates.GetValueOrDefault(address.State, 0m);
        return amount * rate;
    }
}

// Rule 4: Minimal - just what's needed
public class OrderCalculator
{
    private readonly IDiscountPolicy _discount;
    private readonly ITaxCalculator _tax;

    public OrderCalculator(IDiscountPolicy discount, ITaxCalculator tax)
    {
        _discount = discount;
        _tax = tax;
    }

    public decimal CalculateTotal(Order order)
    {
        var subtotal = order.Items.Sum(i => i.Price * i.Quantity);
        var discounted = _discount.Apply(subtotal, order.Customer);
        var tax = _tax.Calculate(discounted, order.ShippingAddress);
        return discounted + tax;
    }
}

// Rule 1: Testable
[Fact]
public void CalculateTotal_WithVipCustomer_AppliesDiscount()
{
    var discount = new VipDiscountPolicy();
    var tax = new StateTaxCalculator();
    var calculator = new OrderCalculator(discount, tax);
    var order = CreateTestOrder(isVip: true, state: "NY", subtotal: 100m);

    var total = calculator.CalculateTotal(order);

    Assert.Equal(97.2m, total); // 100 - 10% + 8% tax
}
```

## Quick Checklist

- [ ] All tests pass (system is verifiable)
- [ ] No duplicated code or concepts
- [ ] Code is expressive and self-documenting
- [ ] Minimal classes and methods needed
- [ ] Uses well-known patterns appropriately
- [ ] Tests serve as documentation

## See Also

- [Unit Tests](08-unit-tests.md) - Writing clean tests
- [Functions](02-functions.md) - Small, expressive functions
- [Classes](09-classes.md) - Small, focused classes
