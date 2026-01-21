# Principle 9: Classes

> "The first rule of classes is that they should be small. The second rule of classes is that they should be smaller than that."
> — Robert C. Martin

## Overview

Classes are the building blocks of object-oriented systems. Clean classes follow the same principles as clean functions: small size, single responsibility, and high cohesion.

## Key Rules

### 9.1 Class Organization

Follow a consistent order within classes:

```csharp
public class OrderProcessor
{
    // 1. Constants (public, then private)
    public const int MaxOrderItems = 100;
    private const decimal TaxRate = 0.13m;

    // 2. Static fields
    private static readonly ILogger _logger = LoggerFactory.Create<OrderProcessor>();

    // 3. Instance fields
    private readonly IOrderRepository _repository;
    private readonly IEmailService _emailService;

    // 4. Constructor(s)
    public OrderProcessor(IOrderRepository repository, IEmailService emailService)
    {
        _repository = repository;
        _emailService = emailService;
    }

    // 5. Public properties
    public int ProcessedCount { get; private set; }

    // 6. Public methods (high-level first)
    public async Task ProcessAsync(Order order)
    {
        ValidateOrder(order);
        await SaveOrderAsync(order);
        await NotifyCustomerAsync(order);
        ProcessedCount++;
    }

    // 7. Private methods (in order of use by public methods)
    private void ValidateOrder(Order order)
    {
        // ...
    }

    private async Task SaveOrderAsync(Order order)
    {
        // ...
    }

    private async Task NotifyCustomerAsync(Order order)
    {
        // ...
    }
}
```

### 9.2 Classes Should Be Small

**Measure by responsibilities, not lines of code.**

**❌ BAD - God Class:**
```csharp
public class SuperDashboard
{
    // UI rendering
    public void Render() { }
    public void UpdateLayout() { }

    // User management
    public void CreateUser() { }
    public void DeleteUser() { }
    public void UpdateUserRole() { }

    // Reporting
    public void GenerateReport() { }
    public void ExportToPdf() { }
    public void ExportToExcel() { }

    // Database operations
    public void SaveToDatabase() { }
    public void LoadFromDatabase() { }

    // Configuration
    public void LoadSettings() { }
    public void SaveSettings() { }

    // Notifications
    public void SendNotification() { }
    public void ScheduleReminder() { }
}
```

**✅ GOOD - Focused classes:**
```csharp
public class DashboardRenderer
{
    public void Render() { }
    public void UpdateLayout() { }
}

public class UserManager
{
    public void CreateUser() { }
    public void DeleteUser() { }
    public void UpdateUserRole() { }
}

public class ReportGenerator
{
    public void Generate() { }
}

public class ReportExporter
{
    public void ExportToPdf(Report report) { }
    public void ExportToExcel(Report report) { }
}

public class NotificationService
{
    public void Send(Notification notification) { }
    public void Schedule(Reminder reminder) { }
}
```

### 9.3 Single Responsibility Principle (SRP)

> **A class should have one, and only one, reason to change.**

**Ask yourself:** "What is this class's responsibility?" If you can't describe it in one sentence without using "and" or "or", it probably has too many responsibilities.

**❌ BAD - Multiple responsibilities:**
```csharp
public class Employee
{
    public string Name { get; set; }
    public decimal Salary { get; set; }

    // Responsibility 1: Calculate pay
    public decimal CalculatePay() { /* ... */ }

    // Responsibility 2: Generate report
    public string GenerateReport() { /* ... */ }

    // Responsibility 3: Save to database
    public void Save() { /* ... */ }
}
```

**✅ GOOD - Single responsibility:**
```csharp
// Responsibility: Represent employee data and behavior
public class Employee
{
    public string Name { get; }
    public Salary Salary { get; }

    public Employee(string name, Salary salary)
    {
        Name = name;
        Salary = salary;
    }
}

// Responsibility: Calculate pay
public class PayrollCalculator
{
    public Money CalculatePay(Employee employee, PayPeriod period)
    {
        // ...
    }
}

// Responsibility: Generate employee reports
public class EmployeeReportGenerator
{
    public Report Generate(Employee employee)
    {
        // ...
    }
}

// Responsibility: Persist employees
public interface IEmployeeRepository
{
    Task SaveAsync(Employee employee);
}
```

### 9.4 Cohesion

**High cohesion:** Methods use most of the instance variables.
**Low cohesion:** Methods use few of the instance variables.

**✅ GOOD - High cohesion:**
```csharp
public class Stack<T>
{
    private readonly List<T> _elements = new();
    private int _topOfStack = -1;

    public int Size => _topOfStack + 1;
    public bool IsEmpty => _topOfStack < 0;

    public void Push(T element)
    {
        _topOfStack++;
        if (_topOfStack >= _elements.Count)
            _elements.Add(element);
        else
            _elements[_topOfStack] = element;
    }

    public T Pop()
    {
        if (IsEmpty)
            throw new InvalidOperationException("Stack is empty");

        return _elements[_topOfStack--];
    }

    public T Peek()
    {
        if (IsEmpty)
            throw new InvalidOperationException("Stack is empty");

        return _elements[_topOfStack];
    }
}
```

All methods use `_elements` and `_topOfStack` - high cohesion!

**❌ BAD - Low cohesion (split this class!):**
```csharp
public class UserUtilities
{
    private readonly IDatabase _database;
    private readonly IEmailClient _emailClient;
    private readonly ILogger _logger;

    // Only uses _database
    public User GetUser(int id) => _database.Find<User>(id);

    // Only uses _emailClient
    public void SendEmail(string to, string body) => _emailClient.Send(to, body);

    // Only uses _logger
    public void LogActivity(string activity) => _logger.Info(activity);
}
```

### 9.5 Maintaining Cohesion Results in Many Small Classes

When a class loses cohesion, split it!

```csharp
// Original class losing cohesion
public class PrintPrimes
{
    public void PrintPrimesInRange(int min, int max)
    {
        for (int i = min; i <= max; i++)
        {
            if (IsPrime(i))
            {
                FormatAndPrint(i);
            }
        }
    }

    private bool IsPrime(int n) { /* ... */ }
    private void FormatAndPrint(int n) { /* ... */ }
}

// Split into cohesive classes
public class PrimeChecker
{
    public bool IsPrime(int n)
    {
        if (n < 2) return false;
        for (int i = 2; i <= Math.Sqrt(n); i++)
        {
            if (n % i == 0) return false;
        }
        return true;
    }
}

public class NumberPrinter
{
    public void Print(int number, string format = "{0}")
    {
        Console.WriteLine(string.Format(format, number));
    }
}

public class PrimeRangePrinter
{
    private readonly PrimeChecker _checker;
    private readonly NumberPrinter _printer;

    public PrimeRangePrinter(PrimeChecker checker, NumberPrinter printer)
    {
        _checker = checker;
        _printer = printer;
    }

    public void PrintPrimesInRange(int min, int max)
    {
        for (int i = min; i <= max; i++)
        {
            if (_checker.IsPrime(i))
            {
                _printer.Print(i);
            }
        }
    }
}
```

### 9.6 Organizing for Change (Open/Closed Principle)

Classes should be open for extension, closed for modification.

**❌ BAD - Must modify class to add features:**
```csharp
public class Sql
{
    public string Create(string table, Column[] columns) { /* ... */ }
    public string Insert(string table, object[] values) { /* ... */ }
    public string Select(string table, string[] columns, Criteria criteria) { /* ... */ }
    public string Update(string table, Column[] columns, Criteria criteria) { /* ... */ }
    public string Delete(string table, Criteria criteria) { /* ... */ }

    // Adding new operation (e.g., UPSERT) requires modifying this class!
}
```

**✅ GOOD - Extend without modifying:**
```csharp
public abstract class SqlStatement
{
    public abstract string Generate();
}

public class CreateStatement : SqlStatement
{
    private readonly string _table;
    private readonly Column[] _columns;

    public CreateStatement(string table, Column[] columns)
    {
        _table = table;
        _columns = columns;
    }

    public override string Generate()
    {
        var columnDefs = string.Join(", ", _columns.Select(c => c.Definition));
        return $"CREATE TABLE {_table} ({columnDefs})";
    }
}

public class SelectStatement : SqlStatement
{
    // Similar structure...
}

// Adding UPSERT - just add new class, no modifications to existing code!
public class UpsertStatement : SqlStatement
{
    public override string Generate()
    {
        // ...
    }
}
```

### 9.7 Isolating from Change (Dependency Inversion)

Depend on abstractions, not concrete implementations.

**❌ BAD - Depends on concrete class:**
```csharp
public class Portfolio
{
    private readonly TokyoStockExchange _exchange;

    public Portfolio()
    {
        _exchange = new TokyoStockExchange();
    }

    public decimal Value()
    {
        // Can't test without hitting real exchange!
        return _holdings.Sum(h => h.Shares * _exchange.GetPrice(h.Symbol));
    }
}
```

**✅ GOOD - Depends on abstraction:**
```csharp
public interface IStockExchange
{
    decimal GetPrice(string symbol);
}

public class Portfolio
{
    private readonly IStockExchange _exchange;
    private readonly List<Holding> _holdings;

    public Portfolio(IStockExchange exchange)
    {
        _exchange = exchange;
    }

    public decimal Value()
    {
        return _holdings.Sum(h => h.Shares * _exchange.GetPrice(h.Symbol));
    }
}

// Easy to test with mock
public class FakeStockExchange : IStockExchange
{
    private readonly Dictionary<string, decimal> _prices = new();

    public void SetPrice(string symbol, decimal price) => _prices[symbol] = price;
    public decimal GetPrice(string symbol) => _prices[symbol];
}
```

### 9.8 Avoid Class Names with "Manager", "Handler", "Processor"

These names often indicate a class that does too much.

**❌ Suspicious names:**
- `UserManager`
- `OrderHandler`
- `DataProcessor`
- `ServiceController`
- `HelperUtility`

**✅ Better alternatives:**
- `UserAuthenticator`, `UserRegistration`, `UserProfileEditor`
- `OrderValidator`, `OrderCalculator`, `OrderShipper`
- `DataParser`, `DataTransformer`, `DataExporter`

## Quick Checklist

- [ ] Classes are small (measured by responsibilities)
- [ ] Single Responsibility - one reason to change
- [ ] High cohesion - methods use most instance variables
- [ ] Consistent organization (constants → fields → constructor → methods)
- [ ] Depends on abstractions, not concretions
- [ ] Open for extension, closed for modification
- [ ] No "God classes" doing everything
- [ ] Class names are specific nouns, not vague "Manager/Handler"

## Size Guidelines

| Metric | Guideline |
|--------|-----------|
| Responsibilities | 1 |
| Public methods | 5-10 max |
| Lines of code | ~100-200 typical, 300 max |
| Dependencies | 3-5 max |
| Instance variables | 5-7 max |

## See Also

- [Functions](02-functions.md) - Small functions within classes
- [Objects and Data Structures](05-objects-and-data-structures.md) - Class vs data structure
- [Clean Architecture - SOLID](../../clean-architecture/principles/02-solid-principles.md) - SOLID principles in depth
