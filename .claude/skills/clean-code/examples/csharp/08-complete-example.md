# Complete Feature Example: All Principles Together

This example demonstrates how Clean Code principles work together in a real feature implementation.

## Scenario: Transfer Money Between Accounts

A banking feature that transfers money between accounts, demonstrating:

- ✅ Meaningful Names
- ✅ Small Functions
- ✅ Error Handling
- ✅ Objects vs Data Structures
- ✅ Clean Classes
- ✅ Unit Tests

---

## Domain Layer

### Value Objects

```csharp
// Meaningful Names: Money instead of decimal
// Objects vs Data: Encapsulates business rules
public readonly record struct Money(decimal Amount, string Currency)
{
    public static Money Zero(string currency) => new(0, currency);
    
    public Money Add(Money other)
    {
        EnsureSameCurrency(other);
        return this with { Amount = Amount + other.Amount };
    }
    
    public Money Subtract(Money other)
    {
        EnsureSameCurrency(other);
        if (Amount < other.Amount)
            throw new InsufficientFundsException(this, other);
        return this with { Amount = Amount - other.Amount };
    }
    
    // Small Function: Single purpose
    private void EnsureSameCurrency(Money other)
    {
        if (Currency != other.Currency)
            throw new CurrencyMismatchException(Currency, other.Currency);
    }
    
    public override string ToString() => $"{Amount:N2} {Currency}";
}

// Meaningful Names: AccountId not Guid
public readonly record struct AccountId(Guid Value)
{
    public static AccountId New() => new(Guid.NewGuid());
    public override string ToString() => Value.ToString()[..8];
}
```

### Domain Entity

```csharp
// Clean Class: High cohesion, single responsibility
public class Account
{
    public AccountId Id { get; }
    public string HolderName { get; }
    public Money Balance { get; private set; }
    
    private readonly List<Transaction> _transactions = [];
    public IReadOnlyList<Transaction> Transactions => _transactions.AsReadOnly();
    
    // Private constructor: Factory method controls creation
    private Account(AccountId id, string holderName, Money initialBalance)
    {
        Id = id;
        HolderName = holderName;
        Balance = initialBalance;
    }
    
    // Meaningful Names: Create not Constructor
    public static Account Create(string holderName, string currency)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(holderName);
        return new Account(AccountId.New(), holderName, Money.Zero(currency));
    }
    
    // Small Function: Does one thing
    public void Deposit(Money amount, string description)
    {
        ValidatePositiveAmount(amount);
        Balance = Balance.Add(amount);
        RecordTransaction(TransactionType.Deposit, amount, description);
    }
    
    public void Withdraw(Money amount, string description)
    {
        ValidatePositiveAmount(amount);
        Balance = Balance.Subtract(amount); // Throws if insufficient
        RecordTransaction(TransactionType.Withdrawal, amount, description);
    }
    
    // Small Functions: Extracted validation
    private static void ValidatePositiveAmount(Money amount)
    {
        if (amount.Amount <= 0)
            throw new InvalidAmountException("Amount must be positive");
    }
    
    private void RecordTransaction(TransactionType type, Money amount, string description)
    {
        _transactions.Add(new Transaction(type, amount, description, DateTimeOffset.UtcNow));
    }
}
```

### Domain Service

```csharp
// Clean Class: Orchestrates operations spanning multiple aggregates
public class MoneyTransferService
{
    // Small Function: Clear, single purpose
    public TransferResult Transfer(
        Account fromAccount,
        Account toAccount,
        Money amount,
        string description)
    {
        // Guard clauses: Fail fast
        ValidateTransferRequest(fromAccount, toAccount, amount);
        
        try
        {
            ExecuteTransfer(fromAccount, toAccount, amount, description);
            return TransferResult.Success(fromAccount.Id, toAccount.Id, amount);
        }
        catch (InsufficientFundsException ex)
        {
            return TransferResult.InsufficientFunds(ex.Available, ex.Requested);
        }
    }
    
    // Extracted: Validation logic
    private static void ValidateTransferRequest(Account from, Account to, Money amount)
    {
        ArgumentNullException.ThrowIfNull(from);
        ArgumentNullException.ThrowIfNull(to);
        
        if (from.Id == to.Id)
            throw new SameAccountTransferException();
            
        if (amount.Amount <= 0)
            throw new InvalidAmountException("Transfer amount must be positive");
    }
    
    // Extracted: Core transfer logic
    private static void ExecuteTransfer(
        Account from, Account to, Money amount, string description)
    {
        var transferDescription = $"Transfer: {description}";
        from.Withdraw(amount, $"To {to.HolderName}: {transferDescription}");
        to.Deposit(amount, $"From {from.HolderName}: {transferDescription}");
    }
}
```

---

## Application Layer

### Command (Data Structure)

```csharp
// Objects vs Data: This IS a data structure (DTO)
public record TransferMoneyCommand(
    Guid FromAccountId,
    Guid ToAccountId,
    decimal Amount,
    string Currency,
    string Description);
```

### Result Pattern (Error Handling)

```csharp
// Error Handling: Result instead of exceptions for expected failures
public record TransferResult
{
    public bool IsSuccess { get; }
    public string? ErrorMessage { get; }
    public TransferDetails? Details { get; }
    
    private TransferResult(bool success, string? error, TransferDetails? details)
    {
        IsSuccess = success;
        ErrorMessage = error;
        Details = details;
    }
    
    public static TransferResult Success(AccountId from, AccountId to, Money amount) =>
        new(true, null, new TransferDetails(from, to, amount));
    
    public static TransferResult InsufficientFunds(Money available, Money requested) =>
        new(false, $"Insufficient funds. Available: {available}, Requested: {requested}", null);
    
    public static TransferResult AccountNotFound(Guid accountId) =>
        new(false, $"Account {accountId} not found", null);
}

public record TransferDetails(AccountId FromAccount, AccountId ToAccount, Money Amount);
```

### Command Handler

```csharp
// Clean Class: Single responsibility - orchestrate use case
public class TransferMoneyHandler(
    IAccountRepository accountRepository,
    MoneyTransferService transferService,
    IUnitOfWork unitOfWork,
    ILogger<TransferMoneyHandler> logger)
{
    public async Task<TransferResult> Handle(
        TransferMoneyCommand command,
        CancellationToken cancellationToken = default)
    {
        // Meaningful Names: Variables describe what they hold
        var fromAccountId = new AccountId(command.FromAccountId);
        var toAccountId = new AccountId(command.ToAccountId);
        var transferAmount = new Money(command.Amount, command.Currency);
        
        // Error Handling: Early return pattern
        var fromAccount = await accountRepository.FindByIdAsync(fromAccountId, cancellationToken);
        if (fromAccount is null)
            return TransferResult.AccountNotFound(command.FromAccountId);
        
        var toAccount = await accountRepository.FindByIdAsync(toAccountId, cancellationToken);
        if (toAccount is null)
            return TransferResult.AccountNotFound(command.ToAccountId);
        
        // Delegate to domain service
        var result = transferService.Transfer(
            fromAccount, toAccount, transferAmount, command.Description);
        
        if (result.IsSuccess)
        {
            await unitOfWork.SaveChangesAsync(cancellationToken);
            LogSuccessfulTransfer(result.Details!);
        }
        
        return result;
    }
    
    // Extracted: Logging as separate concern
    private void LogSuccessfulTransfer(TransferDetails details)
    {
        logger.LogInformation(
            "Transfer completed: {Amount} from {From} to {To}",
            details.Amount,
            details.FromAccount,
            details.ToAccount);
    }
}
```

---

## Custom Exceptions

```csharp
// Error Handling: Domain-specific exceptions with context

public class InsufficientFundsException(Money available, Money requested)
    : DomainException($"Insufficient funds. Available: {available}, Requested: {requested}")
{
    public Money Available { get; } = available;
    public Money Requested { get; } = requested;
}

public class CurrencyMismatchException(string expected, string actual)
    : DomainException($"Currency mismatch. Expected: {expected}, Actual: {actual}")
{
    public string Expected { get; } = expected;
    public string Actual { get; } = actual;
}

public class SameAccountTransferException()
    : DomainException("Cannot transfer to the same account");

public class InvalidAmountException(string message)
    : DomainException(message);

public abstract class DomainException(string message) : Exception(message);
```

---

## Unit Tests

```csharp
// Unit Tests: F.I.R.S.T. principles, AAA pattern

public class MoneyTests
{
    [Fact]
    public void Add_SameCurrency_ReturnsSum()
    {
        // Arrange
        var money1 = new Money(100, "USD");
        var money2 = new Money(50, "USD");
        
        // Act
        var result = money1.Add(money2);
        
        // Assert
        result.Amount.Should().Be(150);
        result.Currency.Should().Be("USD");
    }
    
    [Fact]
    public void Subtract_InsufficientFunds_ThrowsException()
    {
        // Arrange
        var balance = new Money(50, "USD");
        var withdrawal = new Money(100, "USD");
        
        // Act
        var act = () => balance.Subtract(withdrawal);
        
        // Assert
        act.Should().Throw<InsufficientFundsException>()
            .Which.Available.Should().Be(balance);
    }
}

public class AccountTests
{
    [Fact]
    public void Deposit_ValidAmount_IncreasesBalance()
    {
        // Arrange
        var account = Account.Create("John Doe", "USD");
        var depositAmount = new Money(100, "USD");
        
        // Act
        account.Deposit(depositAmount, "Initial deposit");
        
        // Assert
        account.Balance.Amount.Should().Be(100);
        account.Transactions.Should().HaveCount(1);
    }
    
    [Fact]
    public void Withdraw_MoreThanBalance_ThrowsInsufficientFunds()
    {
        // Arrange
        var account = Account.Create("John Doe", "USD");
        account.Deposit(new Money(50, "USD"), "Deposit");
        
        // Act
        var act = () => account.Withdraw(new Money(100, "USD"), "Withdrawal");
        
        // Assert
        act.Should().Throw<InsufficientFundsException>();
    }
}

public class TransferMoneyHandlerTests
{
    private readonly Mock<IAccountRepository> _repositoryMock = new();
    private readonly Mock<IUnitOfWork> _unitOfWorkMock = new();
    private readonly Mock<ILogger<TransferMoneyHandler>> _loggerMock = new();
    private readonly MoneyTransferService _transferService = new();
    
    [Fact]
    public async Task Handle_ValidTransfer_ReturnsSuccess()
    {
        // Arrange
        var fromAccount = Account.Create("Alice", "USD");
        fromAccount.Deposit(new Money(500, "USD"), "Initial");
        
        var toAccount = Account.Create("Bob", "USD");
        
        _repositoryMock.Setup(r => r.FindByIdAsync(fromAccount.Id, It.IsAny<CancellationToken>()))
            .ReturnsAsync(fromAccount);
        _repositoryMock.Setup(r => r.FindByIdAsync(toAccount.Id, It.IsAny<CancellationToken>()))
            .ReturnsAsync(toAccount);
        
        var handler = CreateHandler();
        var command = new TransferMoneyCommand(
            fromAccount.Id.Value, toAccount.Id.Value, 100, "USD", "Test transfer");
        
        // Act
        var result = await handler.Handle(command);
        
        // Assert
        result.IsSuccess.Should().BeTrue();
        fromAccount.Balance.Amount.Should().Be(400);
        toAccount.Balance.Amount.Should().Be(100);
    }
    
    [Fact]
    public async Task Handle_AccountNotFound_ReturnsError()
    {
        // Arrange
        _repositoryMock.Setup(r => r.FindByIdAsync(It.IsAny<AccountId>(), It.IsAny<CancellationToken>()))
            .ReturnsAsync((Account?)null);
        
        var handler = CreateHandler();
        var command = new TransferMoneyCommand(Guid.NewGuid(), Guid.NewGuid(), 100, "USD", "Test");
        
        // Act
        var result = await handler.Handle(command);
        
        // Assert
        result.IsSuccess.Should().BeFalse();
        result.ErrorMessage.Should().Contain("not found");
    }
    
    private TransferMoneyHandler CreateHandler() =>
        new(_repositoryMock.Object, _transferService, _unitOfWorkMock.Object, _loggerMock.Object);
}
```

---

## Principles Summary

| Principle | How Applied |
|-----------|-------------|
| **Meaningful Names** | `Money`, `AccountId`, `TransferResult` instead of primitives |
| **Small Functions** | `ValidateTransferRequest`, `ExecuteTransfer`, `EnsureSameCurrency` |
| **Comments** | None needed — code is self-documenting |
| **Error Handling** | Result pattern for expected failures, exceptions for programming errors |
| **Objects vs Data** | `Account` (object with behavior) vs `TransferMoneyCommand` (data structure) |
| **Classes** | High cohesion: `Account` handles account logic, `MoneyTransferService` handles transfers |
| **Unit Tests** | F.I.R.S.T., AAA pattern, meaningful test names |
| **Modern C#** | Records, primary constructors, pattern matching |

---

## Related

- [Meaningful Names](../../principles/01-meaningful-names.md)
- [Functions](../../principles/02-functions.md)
- [Error Handling](../../principles/06-error-handling.md)
- [Objects and Data Structures](../../principles/05-objects-and-data-structures.md)
- [Classes](../../principles/09-classes.md)
- [Unit Tests](../../principles/08-unit-tests.md)
