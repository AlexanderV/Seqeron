# Principle 7: Boundaries

> "We manage third-party boundaries by having very few places in the code that refer to them."
> — Robert C. Martin

## Overview

We seldom control all the software in our systems. We use third-party packages, open source, internal components from other teams. Managing these boundaries cleanly is essential for maintainable systems.

## Key Rules

### 7.1 Using Third-Party Code

Third-party interfaces often provide more than we need and expose us to changes we don't control.

**Problem:** Direct usage scatters third-party calls throughout codebase.

**❌ BAD - Scattered third-party usage:**
```csharp
public class SensorMonitor
{
    private readonly ThirdPartySensorAPI _api = new();

    public void Monitor()
    {
        var data = _api.FetchRawData();
        var processed = _api.ProcessWith(data, new APISettings
        {
            Format = "JSON",
            Timeout = 5000,
            RetryCount = 3
        });
        // Third-party types and settings everywhere
    }
}

public class SensorReporter
{
    private readonly ThirdPartySensorAPI _api = new();

    public void GenerateReport()
    {
        var data = _api.FetchRawData();
        // Same third-party calls in multiple places
    }
}
```

**✅ GOOD - Wrap in adapter:**
```csharp
// Our interface - we control it
public interface ISensorReader
{
    SensorData Read();
    SensorReport GenerateReport(DateTime from, DateTime to);
}

// Our domain types
public record SensorData(double Temperature, double Humidity, DateTime Timestamp);
public record SensorReport(IReadOnlyList<SensorData> Readings, SensorStatistics Stats);

// Adapter - single place for third-party calls
public class ThirdPartySensorAdapter : ISensorReader
{
    private readonly ThirdPartySensorAPI _api;
    private readonly SensorSettings _settings;

    public ThirdPartySensorAdapter(SensorSettings settings)
    {
        _api = new ThirdPartySensorAPI();
        _settings = settings;
    }

    public SensorData Read()
    {
        try
        {
            var raw = _api.FetchRawData();
            return MapToSensorData(raw);
        }
        catch (APIException ex)
        {
            throw new SensorReadException("Failed to read sensor", ex);
        }
    }

    private SensorData MapToSensorData(RawAPIData raw)
    {
        return new SensorData(
            raw.TempValue / 100.0,  // Convert from API format
            raw.HumidValue / 100.0,
            DateTime.Parse(raw.TimestampStr)
        );
    }
}
```

### 7.2 Benefits of Wrapping

1. **Isolation** - Changes to third-party API affect only the wrapper
2. **Testability** - Mock the interface, not the third-party
3. **Domain language** - Use your terms, not theirs
4. **Exception handling** - Convert to domain exceptions
5. **Feature limitation** - Expose only what you need

### 7.3 Learning Tests

Write tests to explore and understand third-party APIs.

```csharp
public class PaymentGatewayLearningTests
{
    [Fact]
    public void LearningTest_Gateway_AcceptsTestCard()
    {
        var gateway = new ThirdPartyPaymentGateway(TestCredentials);

        var result = gateway.Charge("4111111111111111", 100);

        Assert.True(result.Success);
        Assert.NotNull(result.TransactionId);
    }

    [Fact]
    public void LearningTest_Gateway_RejectsInvalidCard()
    {
        var gateway = new ThirdPartyPaymentGateway(TestCredentials);

        var result = gateway.Charge("0000000000000000", 100);

        Assert.False(result.Success);
        Assert.Equal("INVALID_CARD", result.ErrorCode);
    }

    [Fact]
    public void LearningTest_Gateway_ThrowsOnTimeout()
    {
        var gateway = new ThirdPartyPaymentGateway(TestCredentials);
        gateway.SetTimeout(1); // 1ms - will timeout

        Assert.Throws<GatewayTimeoutException>(() =>
            gateway.Charge("4111111111111111", 100));
    }
}
```

**Benefits of learning tests:**
- Understand API behavior before using it
- Verify our assumptions
- Detect breaking changes when upgrading packages
- Document expected behavior

### 7.4 Using Code That Doesn't Exist Yet

When waiting for another team's API, define your own interface first.

```csharp
// We need notification service, but API isn't ready yet

// Step 1: Define the interface WE need
public interface INotificationService
{
    Task SendAsync(Notification notification);
    Task<NotificationStatus> GetStatusAsync(string notificationId);
}

public record Notification(
    string RecipientId,
    string Channel,  // email, sms, push
    string Subject,
    string Body
);

// Step 2: Create a fake for development/testing
public class FakeNotificationService : INotificationService
{
    private readonly List<Notification> _sent = new();

    public Task SendAsync(Notification notification)
    {
        _sent.Add(notification);
        Console.WriteLine($"[FAKE] Sent {notification.Channel}: {notification.Subject}");
        return Task.CompletedTask;
    }

    public Task<NotificationStatus> GetStatusAsync(string notificationId)
    {
        return Task.FromResult(NotificationStatus.Delivered);
    }

    // For testing
    public IReadOnlyList<Notification> SentNotifications => _sent;
}

// Step 3: When real API is ready, create adapter
public class RealNotificationServiceAdapter : INotificationService
{
    private readonly ExternalNotificationAPI _api;

    public async Task SendAsync(Notification notification)
    {
        var request = MapToApiRequest(notification);
        await _api.SendNotificationAsync(request);
    }

    // ...
}
```

### 7.5 Clean Boundaries

**Rules for managing boundaries:**

1. **Minimize boundary points** - Few places should know about third-party code
2. **Use interfaces we control** - Program to our abstractions
3. **Wrap third-party code** - Adapters isolate external dependencies
4. **Don't pass boundary types** - Convert to domain types at the boundary
5. **Test the boundary** - Learning tests verify behavior

### 7.6 Boundary Patterns

**Adapter Pattern:**
```csharp
public class StripePaymentAdapter : IPaymentGateway
{
    private readonly StripeClient _stripe;

    public async Task<PaymentResult> ChargeAsync(Money amount, CreditCard card)
    {
        var options = new ChargeCreateOptions
        {
            Amount = (long)(amount.Value * 100), // Stripe uses cents
            Currency = amount.Currency.Code.ToLower(),
            Source = card.Token
        };

        var charge = await _stripe.ChargeService.CreateAsync(options);
        return MapToPaymentResult(charge);
    }
}
```

**Facade Pattern:**
```csharp
public class CloudStorageFacade : IFileStorage
{
    private readonly S3Client _s3;
    private readonly AzureBlobClient _azure;
    private readonly StorageConfig _config;

    public async Task<FileInfo> UploadAsync(Stream content, string fileName)
    {
        return _config.Provider switch
        {
            "aws" => await UploadToS3(content, fileName),
            "azure" => await UploadToAzure(content, fileName),
            _ => throw new ConfigurationException($"Unknown provider: {_config.Provider}")
        };
    }
}
```

**Anti-Corruption Layer:**
```csharp
public class LegacyOrderAdapter : IOrderService
{
    private readonly LegacyOrderSystem _legacy;

    public Order GetOrder(OrderId id)
    {
        // Legacy system uses different format
        var legacyOrder = _legacy.FindOrder(id.Value.ToString("N"));
        return TranslateToDomainOrder(legacyOrder);
    }

    private Order TranslateToDomainOrder(LegacyOrder legacy)
    {
        // Translate legacy concepts to our domain
        return new Order(
            new OrderId(Guid.Parse(legacy.ORDER_NUM)),
            TranslateCustomer(legacy.CUST_INFO),
            TranslateItems(legacy.LINE_ITEMS),
            TranslateStatus(legacy.STAT_CODE)
        );
    }
}
```

## Quick Checklist

- [ ] Third-party code wrapped in adapters
- [ ] Programming to our interfaces, not theirs
- [ ] Boundary types not passed through system
- [ ] Domain exceptions wrap third-party exceptions
- [ ] Learning tests written for third-party APIs
- [ ] Minimal boundary points in codebase
- [ ] Fake implementations available for development

## See Also

- [Clean Architecture - Infrastructure Layer](../../clean-architecture/SKILL.md) - Adapter implementation
- [Error Handling](06-error-handling.md) - Wrapping exceptions
- [Objects and Data Structures](05-objects-and-data-structures.md) - Interface design
