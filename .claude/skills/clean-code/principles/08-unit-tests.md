# Principle 8: Unit Tests

> "Test code is just as important as production code."
> — Robert C. Martin

## Overview

Tests are what keep our code flexible, maintainable, and reusable. Without tests, every change is a possible bug. Clean tests follow the same principles as clean production code.

## The Three Laws of TDD

1. **First Law:** You may not write production code until you have written a failing unit test.
2. **Second Law:** You may not write more of a unit test than is sufficient to fail (compilation failures count as failures).
3. **Third Law:** You may not write more production code than is sufficient to pass the currently failing test.

## F.I.R.S.T. Principles

Clean tests follow these principles:

| Principle | Description |
|-----------|-------------|
| **F**ast | Tests should run quickly (milliseconds) |
| **I**ndependent | Tests should not depend on each other |
| **R**epeatable | Tests should work in any environment |
| **S**elf-Validating | Tests should have boolean output (pass/fail) |
| **T**imely | Tests should be written just before production code |

## Key Rules

### 8.1 One Concept per Test

Each test should verify one behavior or concept.

**❌ BAD - Multiple concepts:**
```csharp
[Fact]
public void TestCustomer()
{
    var customer = new Customer("John", "john@example.com");

    // Testing creation
    Assert.Equal("John", customer.Name);
    Assert.Equal("john@example.com", customer.Email);

    // Testing validation
    Assert.True(customer.IsValid());

    // Testing behavior
    customer.PlaceOrder(new Order());
    Assert.Equal(1, customer.OrderCount);
}
```

**✅ GOOD - One concept per test:**
```csharp
public class CustomerTests
{
    [Fact]
    public void Constructor_ShouldSetName()
    {
        var customer = new Customer("John", "john@example.com");

        Assert.Equal("John", customer.Name);
    }

    [Fact]
    public void Constructor_ShouldSetEmail()
    {
        var customer = new Customer("John", "john@example.com");

        Assert.Equal("john@example.com", customer.Email);
    }

    [Fact]
    public void IsValid_ShouldReturnTrue_WhenDataIsComplete()
    {
        var customer = new Customer("John", "john@example.com");

        Assert.True(customer.IsValid());
    }

    [Fact]
    public void PlaceOrder_ShouldIncrementOrderCount()
    {
        var customer = new Customer("John", "john@example.com");

        customer.PlaceOrder(new Order());

        Assert.Equal(1, customer.OrderCount);
    }
}
```

### 8.2 Arrange-Act-Assert (AAA) Pattern

Structure tests clearly with three sections.

```csharp
[Fact]
public void CalculateTotal_ShouldApplyDiscount_WhenCustomerIsVip()
{
    // Arrange
    var customer = new Customer("John") { IsVip = true };
    var order = new Order();
    order.AddItem(new Product("Widget", 100m), quantity: 2);
    var calculator = new OrderCalculator();

    // Act
    var total = calculator.CalculateTotal(order, customer);

    // Assert
    Assert.Equal(180m, total); // 10% VIP discount
}
```

### 8.3 Clean Test Naming

Test names should describe the scenario and expected outcome.

**Format:** `MethodName_Scenario_ExpectedBehavior`

```csharp
// Good test names
public void GetCustomer_WithValidId_ReturnsCustomer()
public void GetCustomer_WithInvalidId_ThrowsNotFoundException()
public void PlaceOrder_WhenInventoryInsufficient_ThrowsException()
public void CalculateDiscount_ForNewCustomer_ReturnsZero()
public void SendEmail_WithEmptyRecipient_ThrowsArgumentException()
```

### 8.4 Given-When-Then for BDD Style

Alternative naming using behavior-driven format.

```csharp
public class OrderProcessingTests
{
    [Fact]
    public void GivenVipCustomer_WhenPlacingOrder_ThenDiscountIsApplied()
    {
        // Given
        var customer = CustomerBuilder.Vip().Build();
        var order = OrderBuilder.WithItems(2).Build();

        // When
        var result = _processor.Process(order, customer);

        // Then
        Assert.Equal(0.10m, result.DiscountApplied);
    }
}
```

### 8.5 Use Test Fixtures and Builders

Reduce test setup duplication with builders.

**❌ BAD - Repeated setup:**
```csharp
[Fact]
public void Test1()
{
    var customer = new Customer();
    customer.Name = "John";
    customer.Email = "john@example.com";
    customer.IsVip = true;
    customer.Address = new Address { City = "NYC", Country = "USA" };
    // Test...
}

[Fact]
public void Test2()
{
    var customer = new Customer();
    customer.Name = "Jane";
    customer.Email = "jane@example.com";
    customer.IsVip = true;
    customer.Address = new Address { City = "LA", Country = "USA" };
    // Test...
}
```

**✅ GOOD - Use builders:**
```csharp
public class CustomerBuilder
{
    private string _name = "Default Name";
    private string _email = "default@example.com";
    private bool _isVip = false;

    public static CustomerBuilder Create() => new();
    public static CustomerBuilder Vip() => new CustomerBuilder().AsVip();

    public CustomerBuilder WithName(string name)
    {
        _name = name;
        return this;
    }

    public CustomerBuilder WithEmail(string email)
    {
        _email = email;
        return this;
    }

    public CustomerBuilder AsVip()
    {
        _isVip = true;
        return this;
    }

    public Customer Build()
    {
        return new Customer(_name, _email) { IsVip = _isVip };
    }
}

// Usage in tests
[Fact]
public void Test_WithVipCustomer()
{
    var customer = CustomerBuilder.Vip()
        .WithName("John")
        .Build();
    // Test...
}
```

### 8.6 Test Only Public Interface

Test behavior, not implementation details.

**❌ BAD - Testing private methods:**
```csharp
[Fact]
public void ValidateEmail_ShouldReturnTrue_ForValidFormat()
{
    // Using reflection to test private method
    var method = typeof(Customer).GetMethod("ValidateEmail",
        BindingFlags.NonPublic | BindingFlags.Instance);
    // ...
}
```

**✅ GOOD - Test through public interface:**
```csharp
[Fact]
public void Constructor_WithInvalidEmail_ThrowsArgumentException()
{
    Assert.Throws<ArgumentException>(() =>
        new Customer("John", "invalid-email"));
}
```

### 8.7 Don't Test External Dependencies

Mock external dependencies to test your code, not theirs.

```csharp
public class OrderServiceTests
{
    private readonly Mock<IOrderRepository> _mockRepository;
    private readonly Mock<IEmailService> _mockEmailService;
    private readonly OrderService _sut; // System Under Test

    public OrderServiceTests()
    {
        _mockRepository = new Mock<IOrderRepository>();
        _mockEmailService = new Mock<IEmailService>();
        _sut = new OrderService(_mockRepository.Object, _mockEmailService.Object);
    }

    [Fact]
    public async Task PlaceOrder_ShouldSaveOrder()
    {
        // Arrange
        var order = new Order();
        _mockRepository
            .Setup(r => r.SaveAsync(It.IsAny<Order>()))
            .ReturnsAsync(order);

        // Act
        await _sut.PlaceOrder(order);

        // Assert
        _mockRepository.Verify(r => r.SaveAsync(order), Times.Once);
    }

    [Fact]
    public async Task PlaceOrder_ShouldSendConfirmationEmail()
    {
        // Arrange
        var order = new Order { CustomerEmail = "test@example.com" };

        // Act
        await _sut.PlaceOrder(order);

        // Assert
        _mockEmailService.Verify(
            e => e.SendAsync("test@example.com", It.IsAny<string>()),
            Times.Once);
    }
}
```

### 8.8 Use Assertions Wisely

Make assertions clear and specific.

**❌ BAD - Vague assertions:**
```csharp
Assert.NotNull(result);
Assert.True(result.Count > 0);
```

**✅ GOOD - Specific assertions:**
```csharp
Assert.NotNull(result);
Assert.Equal(3, result.Count);
Assert.Contains(result, c => c.Name == "Expected Customer");
```

### 8.9 Test Edge Cases and Error Conditions

```csharp
public class CalculatorTests
{
    [Fact]
    public void Divide_ByZero_ThrowsDivideByZeroException()
    {
        var calc = new Calculator();

        Assert.Throws<DivideByZeroException>(() => calc.Divide(10, 0));
    }

    [Fact]
    public void Add_WithMaxValues_HandlesOverflow()
    {
        var calc = new Calculator();

        Assert.Throws<OverflowException>(() =>
            calc.Add(int.MaxValue, 1));
    }

    [Theory]
    [InlineData(null)]
    [InlineData("")]
    [InlineData("   ")]
    public void ProcessName_WithEmptyInput_ThrowsArgumentException(string input)
    {
        var processor = new NameProcessor();

        Assert.Throws<ArgumentException>(() => processor.Process(input));
    }
}
```

### 8.10 Parameterized Tests

Use theory tests for multiple inputs.

```csharp
public class ValidationTests
{
    [Theory]
    [InlineData("john@example.com", true)]
    [InlineData("jane@company.org", true)]
    [InlineData("invalid", false)]
    [InlineData("missing@", false)]
    [InlineData("@nodomain.com", false)]
    public void IsValidEmail_ShouldValidateCorrectly(string email, bool expected)
    {
        var validator = new EmailValidator();

        var result = validator.IsValid(email);

        Assert.Equal(expected, result);
    }

    [Theory]
    [MemberData(nameof(DiscountTestData))]
    public void CalculateDiscount_ShouldReturnCorrectAmount(
        Customer customer, decimal orderTotal, decimal expectedDiscount)
    {
        var calculator = new DiscountCalculator();

        var discount = calculator.Calculate(customer, orderTotal);

        Assert.Equal(expectedDiscount, discount);
    }

    public static IEnumerable<object[]> DiscountTestData()
    {
        yield return new object[] { new Customer { IsVip = false }, 100m, 0m };
        yield return new object[] { new Customer { IsVip = true }, 100m, 10m };
        yield return new object[] { new Customer { IsVip = true }, 1000m, 150m }; // Higher discount for large orders
    }
}
```

## Quick Checklist

- [ ] Tests are fast (run in milliseconds)
- [ ] Tests are independent (no shared state)
- [ ] Tests are repeatable (work in any environment)
- [ ] One concept per test
- [ ] AAA pattern followed (Arrange-Act-Assert)
- [ ] Clear, descriptive test names
- [ ] External dependencies mocked
- [ ] Edge cases and errors tested
- [ ] No testing of private methods
- [ ] Test builders reduce duplication

## See Also

- [C# Testing Examples](../examples/csharp/06-testing.md) - Detailed C# examples
- [Functions](02-functions.md) - Testable function design
- [Objects and Data Structures](05-objects-and-data-structures.md) - Mockable interfaces
