# Principle 11: Concurrency

> "Objects are abstractions of processing. Threads are abstractions of schedule."
> — James O. Coplien

## Overview

Concurrent programming is challenging. Writing clean concurrent code is even more challenging. This principle covers the essential practices for writing maintainable, correct concurrent code.

## Why Concurrency is Hard

Concurrency introduces non-determinism. The same code can produce different results depending on thread timing.

**Common problems:**
- Race conditions
- Deadlocks
- Starvation
- Memory visibility issues
- Hard-to-reproduce bugs

## Key Rules

### 11.1 Single Responsibility for Concurrency

Separate concurrency code from other code.

**❌ BAD - Mixed concerns:**
```csharp
public class OrderProcessor
{
    private readonly object _lock = new();
    private int _processedCount;

    public async Task ProcessAsync(Order order)
    {
        // Business logic mixed with concurrency
        lock (_lock)
        {
            ValidateOrder(order);
            CalculateTotal(order);
            _processedCount++;
        }

        await SaveOrderAsync(order);
    }
}
```

**✅ GOOD - Separated concerns:**
```csharp
// Business logic - no concurrency knowledge
public class OrderProcessor
{
    public void Process(Order order)
    {
        ValidateOrder(order);
        CalculateTotal(order);
    }
}

// Concurrency handled separately
public class ConcurrentOrderService
{
    private readonly OrderProcessor _processor;
    private readonly SemaphoreSlim _semaphore = new(10); // Max 10 concurrent
    private int _processedCount;

    public async Task ProcessAsync(Order order)
    {
        await _semaphore.WaitAsync();
        try
        {
            _processor.Process(order);
            Interlocked.Increment(ref _processedCount);
            await SaveOrderAsync(order);
        }
        finally
        {
            _semaphore.Release();
        }
    }
}
```

### 11.2 Limit Scope of Shared Data

Minimize the amount of shared mutable state.

**❌ BAD - Too much shared state:**
```csharp
public class ShoppingCart
{
    public List<CartItem> Items { get; } = new();
    public decimal Total { get; set; }
    public Customer Customer { get; set; }
    public DateTime LastModified { get; set; }
    // All mutable, all potentially shared!
}
```

**✅ GOOD - Encapsulated, minimal shared state:**
```csharp
public class ShoppingCart
{
    private readonly object _lock = new();
    private readonly List<CartItem> _items = new();

    public IReadOnlyList<CartItem> Items
    {
        get
        {
            lock (_lock)
            {
                return _items.ToList(); // Return copy
            }
        }
    }

    public void AddItem(CartItem item)
    {
        lock (_lock)
        {
            _items.Add(item);
        }
    }

    public decimal CalculateTotal()
    {
        lock (_lock)
        {
            return _items.Sum(i => i.Price * i.Quantity);
        }
    }
}
```

### 11.3 Use Immutable Objects

Immutable objects are inherently thread-safe.

**✅ GOOD - Immutable types:**
```csharp
// Immutable record
public record OrderSummary(
    OrderId Id,
    CustomerId CustomerId,
    IReadOnlyList<OrderItem> Items,
    Money Total,
    DateTime CreatedAt
);

// Immutable class
public sealed class Money
{
    public decimal Amount { get; }
    public Currency Currency { get; }

    public Money(decimal amount, Currency currency)
    {
        Amount = amount;
        Currency = currency;
    }

    // Returns new instance instead of mutating
    public Money Add(Money other)
    {
        if (Currency != other.Currency)
            throw new InvalidOperationException("Currency mismatch");

        return new Money(Amount + other.Amount, Currency);
    }

    public Money Multiply(decimal factor)
    {
        return new Money(Amount * factor, Currency);
    }
}

// Usage - no synchronization needed!
var price = new Money(100, Currency.USD);
var taxedPrice = price.Multiply(1.08m); // New object, original unchanged
```

### 11.4 Use Thread-Safe Collections

.NET provides thread-safe collections in `System.Collections.Concurrent`.

```csharp
// Thread-safe alternatives
ConcurrentDictionary<string, Customer> _customers = new();
ConcurrentQueue<Order> _orderQueue = new();
ConcurrentBag<LogEntry> _logs = new();
BlockingCollection<WorkItem> _workItems = new();

// Channel for producer-consumer (modern approach)
Channel<Order> _orderChannel = Channel.CreateUnbounded<Order>();

// Example: Producer-Consumer with Channel
public class OrderProcessor
{
    private readonly Channel<Order> _channel = Channel.CreateBounded<Order>(100);

    public async Task EnqueueAsync(Order order)
    {
        await _channel.Writer.WriteAsync(order);
    }

    public async Task ProcessQueueAsync(CancellationToken ct)
    {
        await foreach (var order in _channel.Reader.ReadAllAsync(ct))
        {
            await ProcessOrderAsync(order);
        }
    }
}
```

### 11.5 Async/Await Best Practices

Modern C# concurrency uses async/await.

**❌ BAD - Blocking async code:**
```csharp
public void ProcessOrder(Order order)
{
    // Blocks thread, can cause deadlocks!
    var customer = GetCustomerAsync(order.CustomerId).Result;
    var result = ProcessPaymentAsync(order).GetAwaiter().GetResult();
}
```

**✅ GOOD - Async all the way:**
```csharp
public async Task ProcessOrderAsync(Order order)
{
    var customer = await GetCustomerAsync(order.CustomerId);
    var result = await ProcessPaymentAsync(order);
}
```

**✅ GOOD - ConfigureAwait for library code:**
```csharp
// In library code, use ConfigureAwait(false) to avoid deadlocks
public async Task<Customer> GetCustomerAsync(int id)
{
    var data = await _httpClient.GetAsync($"/api/customers/{id}")
        .ConfigureAwait(false);

    return await data.Content.ReadFromJsonAsync<Customer>()
        .ConfigureAwait(false);
}
```

**✅ GOOD - Parallel execution when possible:**
```csharp
public async Task<OrderSummary> GetOrderSummaryAsync(OrderId orderId)
{
    // Execute independent operations in parallel
    var orderTask = _orderRepository.GetByIdAsync(orderId);
    var customerTask = _customerRepository.GetByOrderIdAsync(orderId);
    var paymentsTask = _paymentRepository.GetByOrderIdAsync(orderId);

    await Task.WhenAll(orderTask, customerTask, paymentsTask);

    return new OrderSummary(
        orderTask.Result,
        customerTask.Result,
        paymentsTask.Result
    );
}
```

### 11.6 Avoid Lock Contention

Design to minimize time spent waiting for locks.

**❌ BAD - Coarse-grained locking:**
```csharp
public class OrderService
{
    private readonly object _lock = new();

    public void ProcessOrder(Order order)
    {
        lock (_lock) // Everything locked!
        {
            ValidateOrder(order);           // Could be lock-free
            CalculateTotal(order);          // Could be lock-free
            SaveToDatabase(order);          // Long operation while locked!
            SendConfirmationEmail(order);   // Network I/O while locked!
        }
    }
}
```

**✅ GOOD - Fine-grained locking:**
```csharp
public class OrderService
{
    private readonly object _counterLock = new();
    private int _processedCount;

    public async Task ProcessOrderAsync(Order order)
    {
        // Lock-free operations
        ValidateOrder(order);
        CalculateTotal(order);

        // Async I/O without locks
        await SaveToDatabaseAsync(order);
        await SendConfirmationEmailAsync(order);

        // Only lock what needs synchronization
        lock (_counterLock)
        {
            _processedCount++;
        }
    }
}
```

### 11.7 Use Interlocked for Simple Operations

For simple atomic operations, use `Interlocked` instead of locks.

```csharp
public class Counter
{
    private int _count;

    // ❌ BAD - Lock for simple increment
    public void IncrementWithLock()
    {
        lock (this) { _count++; }
    }

    // ✅ GOOD - Interlocked
    public void Increment()
    {
        Interlocked.Increment(ref _count);
    }

    public int GetValue()
    {
        return Interlocked.CompareExchange(ref _count, 0, 0);
    }

    // Atomic compare-and-swap
    public bool TryUpdateIfEquals(int expected, int newValue)
    {
        return Interlocked.CompareExchange(ref _count, newValue, expected) == expected;
    }
}
```

### 11.8 Cancellation Support

Always support cancellation for async operations.

```csharp
public class OrderProcessor
{
    public async Task ProcessBatchAsync(
        IEnumerable<Order> orders,
        CancellationToken cancellationToken = default)
    {
        foreach (var order in orders)
        {
            // Check cancellation regularly
            cancellationToken.ThrowIfCancellationRequested();

            await ProcessOrderAsync(order, cancellationToken);
        }
    }

    public async Task ProcessOrderAsync(
        Order order,
        CancellationToken cancellationToken = default)
    {
        // Pass token to all async operations
        var customer = await _customerService
            .GetByIdAsync(order.CustomerId, cancellationToken);

        await _paymentService
            .ChargeAsync(order.Total, cancellationToken);

        await _repository
            .SaveAsync(order, cancellationToken);
    }
}

// Usage with timeout
using var cts = new CancellationTokenSource(TimeSpan.FromSeconds(30));
await processor.ProcessBatchAsync(orders, cts.Token);
```

### 11.9 Avoid Deadlocks

**Deadlock conditions:**
1. Mutual exclusion
2. Hold and wait
3. No preemption
4. Circular wait

**❌ BAD - Potential deadlock:**
```csharp
public class Transfer
{
    public void TransferMoney(Account from, Account to, decimal amount)
    {
        lock (from)
        {
            lock (to)  // Thread 1: lock A then B
            {           // Thread 2: lock B then A → DEADLOCK!
                from.Withdraw(amount);
                to.Deposit(amount);
            }
        }
    }
}
```

**✅ GOOD - Consistent lock ordering:**
```csharp
public class Transfer
{
    public void TransferMoney(Account from, Account to, decimal amount)
    {
        // Always lock in consistent order (by ID)
        var first = from.Id < to.Id ? from : to;
        var second = from.Id < to.Id ? to : from;

        lock (first)
        {
            lock (second)
            {
                from.Withdraw(amount);
                to.Deposit(amount);
            }
        }
    }
}
```

**✅ BETTER - Use async with SemaphoreSlim:**
```csharp
public class Account
{
    private readonly SemaphoreSlim _semaphore = new(1, 1);
    private decimal _balance;

    public async Task<bool> TryWithdrawAsync(decimal amount, CancellationToken ct)
    {
        if (!await _semaphore.WaitAsync(TimeSpan.FromSeconds(5), ct))
            return false;

        try
        {
            if (_balance < amount)
                return false;

            _balance -= amount;
            return true;
        }
        finally
        {
            _semaphore.Release();
        }
    }
}
```

### 11.10 Testing Concurrent Code

Concurrent code requires special testing techniques.

```csharp
public class ConcurrentCounterTests
{
    [Fact]
    public async Task Counter_ShouldBeThreadSafe()
    {
        // Arrange
        var counter = new ThreadSafeCounter();
        const int iterations = 1000;
        const int parallelism = 10;

        // Act - multiple threads incrementing
        var tasks = Enumerable.Range(0, parallelism)
            .Select(_ => Task.Run(() =>
            {
                for (int i = 0; i < iterations; i++)
                {
                    counter.Increment();
                }
            }));

        await Task.WhenAll(tasks);

        // Assert
        Assert.Equal(iterations * parallelism, counter.Value);
    }

    [Fact]
    public async Task Service_ShouldHandleCancellation()
    {
        // Arrange
        var service = new LongRunningService();
        using var cts = new CancellationTokenSource();

        // Act
        var task = service.ProcessAsync(cts.Token);
        await Task.Delay(100);
        cts.Cancel();

        // Assert
        await Assert.ThrowsAsync<OperationCanceledException>(() => task);
    }
}
```

## Quick Checklist

- [ ] Concurrency code separated from business logic
- [ ] Shared mutable state minimized
- [ ] Immutable objects used where possible
- [ ] Thread-safe collections used (ConcurrentDictionary, Channel)
- [ ] Async/await used correctly (no .Result or .Wait())
- [ ] ConfigureAwait(false) in library code
- [ ] Cancellation tokens supported
- [ ] Fine-grained locking (short critical sections)
- [ ] Consistent lock ordering to prevent deadlocks
- [ ] Concurrent code tested for thread safety

## Common Patterns

| Pattern | Use Case | Example |
|---------|----------|---------|
| Immutability | Shared data | `record`, `readonly struct` |
| Producer-Consumer | Work queues | `Channel<T>`, `BlockingCollection<T>` |
| Reader-Writer | Read-heavy workloads | `ReaderWriterLockSlim` |
| Parallel Processing | CPU-bound work | `Parallel.ForEach`, `PLINQ` |
| Async I/O | I/O-bound work | `async/await` |

## See Also

- [Error Handling](06-error-handling.md) - Exception handling in async code
- [Classes](09-classes.md) - Thread-safe class design
- [Unit Tests](08-unit-tests.md) - Testing concurrent code
