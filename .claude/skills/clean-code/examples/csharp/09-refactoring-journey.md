# Refactoring Journey: Step-by-Step Transformation

A complete refactoring example showing the evolution of code through multiple iterations.

## The Scenario

We have a `ReportGenerator` class that has grown messy over time. We'll transform it step by step, applying Clean Code principles at each stage.

---

## Step 0: Original Code (All Problems)

```csharp
public class ReportGenerator
{
    public string Generate(int type, List<object> data, bool pdf, string template)
    {
        string result = "";
        
        // Check if data exists
        if (data == null || data.Count == 0)
        {
            return "No data";
        }
        
        // Process based on type
        if (type == 1) // Sales report
        {
            double total = 0;
            for (int i = 0; i < data.Count; i++)
            {
                var item = data[i];
                if (item is Dictionary<string, object> d)
                {
                    if (d.ContainsKey("amount"))
                    {
                        var amt = d["amount"];
                        if (amt != null)
                        {
                            try
                            {
                                total += Convert.ToDouble(amt);
                            }
                            catch
                            {
                                // Skip invalid
                            }
                        }
                    }
                }
            }
            
            result = "Sales Report\n";
            result += "Total: $" + total.ToString("0.00") + "\n";
            result += "Items: " + data.Count.ToString() + "\n";
            
            // Add details
            for (int i = 0; i < data.Count; i++)
            {
                var item = data[i];
                if (item is Dictionary<string, object> d)
                {
                    result += "- ";
                    if (d.ContainsKey("name")) result += d["name"].ToString();
                    result += ": $";
                    if (d.ContainsKey("amount")) result += Convert.ToDouble(d["amount"]).ToString("0.00");
                    result += "\n";
                }
            }
        }
        else if (type == 2) // Inventory report
        {
            int totalItems = 0;
            int lowStock = 0;
            
            for (int i = 0; i < data.Count; i++)
            {
                var item = data[i];
                if (item is Dictionary<string, object> d)
                {
                    if (d.ContainsKey("quantity"))
                    {
                        int qty = Convert.ToInt32(d["quantity"]);
                        totalItems += qty;
                        if (qty < 10) lowStock++;
                    }
                }
            }
            
            result = "Inventory Report\n";
            result += "Total Items: " + totalItems.ToString() + "\n";
            result += "Low Stock Items: " + lowStock.ToString() + "\n";
        }
        else
        {
            result = "Unknown report type";
        }
        
        // Convert to PDF if needed
        if (pdf)
        {
            // PDF conversion logic...
            result = "<PDF>" + result + "</PDF>";
        }
        
        return result;
    }
}
```

### Problems Identified:

| Problem | Principle Violated |
|---------|-------------------|
| Magic numbers (1, 2, 10) | Meaningful Names |
| God method (80+ lines) | Functions |
| `type` parameter as int | Meaningful Names |
| Empty catch block | Error Handling |
| String concatenation | Modern C# |
| `object` and `Dictionary` | Objects vs Data |
| Nested if statements | Functions |
| No separation of concerns | Classes |

---

## Step 1: Meaningful Names

Replace magic numbers with enums, improve variable names:

```csharp
public enum ReportType
{
    Sales,
    Inventory
}

public class ReportGenerator
{
    private const int LowStockThreshold = 10;
    
    public string Generate(
        ReportType reportType, 
        List<object> reportData, 
        bool convertToPdf, 
        string templateName)
    {
        string reportContent = "";
        
        if (reportData == null || reportData.Count == 0)
        {
            return "No data";
        }
        
        if (reportType == ReportType.Sales)
        {
            double salesTotal = 0;
            for (int itemIndex = 0; itemIndex < reportData.Count; itemIndex++)
            {
                var currentItem = reportData[itemIndex];
                if (currentItem is Dictionary<string, object> itemData)
                {
                    if (itemData.ContainsKey("amount"))
                    {
                        var itemAmount = itemData["amount"];
                        if (itemAmount != null)
                        {
                            try
                            {
                                salesTotal += Convert.ToDouble(itemAmount);
                            }
                            catch
                            {
                                // Skip invalid amounts
                            }
                        }
                    }
                }
            }
            // ... rest similar
        }
        // ...
    }
}
```

**What Changed:**
- `type` → `ReportType` enum
- `1, 2` → `ReportType.Sales`, `ReportType.Inventory`
- `10` → `LowStockThreshold` constant
- `i` → `itemIndex`
- `d` → `itemData`
- `amt` → `itemAmount`

---

## Step 2: Strongly Typed Data Models

Replace `Dictionary<string, object>` with proper classes:

```csharp
// Domain models
public record SalesItem(string Name, decimal Amount);
public record InventoryItem(string Sku, string Name, int Quantity);

// Report data wrappers
public record SalesReportData(IReadOnlyList<SalesItem> Items);
public record InventoryReportData(IReadOnlyList<InventoryItem> Items);

public class ReportGenerator
{
    private const int LowStockThreshold = 10;
    
    public string GenerateSalesReport(SalesReportData data, bool convertToPdf)
    {
        if (data.Items.Count == 0)
            return "No data";
        
        decimal salesTotal = data.Items.Sum(item => item.Amount);
        
        var reportContent = new StringBuilder();
        reportContent.AppendLine("Sales Report");
        reportContent.AppendLine($"Total: {salesTotal:C}");
        reportContent.AppendLine($"Items: {data.Items.Count}");
        
        foreach (var item in data.Items)
        {
            reportContent.AppendLine($"- {item.Name}: {item.Amount:C}");
        }
        
        return FormatOutput(reportContent.ToString(), convertToPdf);
    }
    
    public string GenerateInventoryReport(InventoryReportData data, bool convertToPdf)
    {
        if (data.Items.Count == 0)
            return "No data";
        
        int totalQuantity = data.Items.Sum(item => item.Quantity);
        int lowStockCount = data.Items.Count(item => item.Quantity < LowStockThreshold);
        
        var reportContent = new StringBuilder();
        reportContent.AppendLine("Inventory Report");
        reportContent.AppendLine($"Total Items: {totalQuantity}");
        reportContent.AppendLine($"Low Stock Items: {lowStockCount}");
        
        return FormatOutput(reportContent.ToString(), convertToPdf);
    }
    
    private string FormatOutput(string content, bool convertToPdf)
    {
        return convertToPdf ? $"<PDF>{content}</PDF>" : content;
    }
}
```

**What Changed:**
- `List<object>` → `SalesReportData`, `InventoryReportData`
- `Dictionary<string, object>` → `SalesItem`, `InventoryItem` records
- Split into separate methods per report type
- Extracted `FormatOutput` method
- Used `StringBuilder` instead of string concatenation
- Used LINQ instead of manual loops
- Used string interpolation with formatting

---

## Step 3: Extract Functions & Single Responsibility

Each method does one thing:

```csharp
public class SalesReportGenerator
{
    public string Generate(SalesReportData data, ReportFormat format)
    {
        if (data.Items.Count == 0)
            return CreateEmptyReport();
        
        var summary = CalculateSummary(data);
        var content = BuildReportContent(summary, data.Items);
        
        return FormatReport(content, format);
    }
    
    private static string CreateEmptyReport() => "No data available";
    
    private static SalesSummary CalculateSummary(SalesReportData data)
    {
        return new SalesSummary(
            TotalAmount: data.Items.Sum(i => i.Amount),
            ItemCount: data.Items.Count);
    }
    
    private static string BuildReportContent(SalesSummary summary, IEnumerable<SalesItem> items)
    {
        var builder = new StringBuilder();
        
        AppendHeader(builder);
        AppendSummary(builder, summary);
        AppendItemDetails(builder, items);
        
        return builder.ToString();
    }
    
    private static void AppendHeader(StringBuilder builder)
    {
        builder.AppendLine("Sales Report");
        builder.AppendLine(new string('=', 40));
    }
    
    private static void AppendSummary(StringBuilder builder, SalesSummary summary)
    {
        builder.AppendLine($"Total: {summary.TotalAmount:C}");
        builder.AppendLine($"Items: {summary.ItemCount}");
        builder.AppendLine();
    }
    
    private static void AppendItemDetails(StringBuilder builder, IEnumerable<SalesItem> items)
    {
        foreach (var item in items)
        {
            builder.AppendLine($"- {item.Name}: {item.Amount:C}");
        }
    }
    
    private static string FormatReport(string content, ReportFormat format)
    {
        return format switch
        {
            ReportFormat.Pdf => PdfConverter.Convert(content),
            ReportFormat.Html => HtmlConverter.Convert(content),
            ReportFormat.Text => content,
            _ => content
        };
    }
}

public enum ReportFormat { Text, Pdf, Html }
public record SalesSummary(decimal TotalAmount, int ItemCount);
```

**What Changed:**
- Split into multiple small functions
- Each function does ONE thing
- Function names are descriptive
- Separated calculation from formatting
- `bool convertToPdf` → `ReportFormat` enum
- No function longer than 10-15 lines

---

## Step 4: Proper Error Handling

Handle edge cases and provide context:

```csharp
public class SalesReportGenerator(ILogger<SalesReportGenerator> logger)
{
    public Result<ReportOutput> Generate(SalesReportData data, ReportFormat format)
    {
        try
        {
            return GenerateInternal(data, format);
        }
        catch (Exception ex)
        {
            logger.LogError(ex, "Failed to generate sales report");
            return Result<ReportOutput>.Failure(
                Error.Internal("Failed to generate report. Please try again."));
        }
    }
    
    private Result<ReportOutput> GenerateInternal(SalesReportData data, ReportFormat format)
    {
        ArgumentNullException.ThrowIfNull(data);
        
        if (data.Items.Count == 0)
            return Result<ReportOutput>.Success(ReportOutput.Empty("No sales data available"));
        
        var invalidItems = data.Items.Where(i => i.Amount < 0).ToList();
        if (invalidItems.Any())
        {
            logger.LogWarning("Found {Count} items with negative amounts", invalidItems.Count);
            return Result<ReportOutput>.Failure(
                Error.Validation($"Found {invalidItems.Count} items with invalid amounts"));
        }
        
        var summary = CalculateSummary(data);
        var content = BuildReportContent(summary, data.Items);
        var formatted = FormatReport(content, format);
        
        return Result<ReportOutput>.Success(new ReportOutput(formatted, format));
    }
    
    // ... other methods
}

public record ReportOutput(string Content, ReportFormat Format)
{
    public static ReportOutput Empty(string message) => new(message, ReportFormat.Text);
}
```

**What Changed:**
- Added Result pattern for error handling
- Guard clauses with meaningful exceptions
- Validation of input data
- Logging for debugging
- No empty catch blocks

---

## Step 5: Final Clean Version

Complete refactored solution with all principles:

```csharp
// Domain Models
public record SalesItem(string Name, decimal Amount)
{
    public bool IsValid => Amount >= 0;
}

public record SalesReportData(IReadOnlyList<SalesItem> Items)
{
    public bool IsEmpty => Items.Count == 0;
    public IEnumerable<SalesItem> ValidItems => Items.Where(i => i.IsValid);
    public IEnumerable<SalesItem> InvalidItems => Items.Where(i => !i.IsValid);
}

public record SalesSummary(decimal TotalAmount, int ItemCount, decimal AverageAmount);

public record ReportOutput(string Content, ReportFormat Format, DateTimeOffset GeneratedAt)
{
    public static ReportOutput Empty(string message) => 
        new(message, ReportFormat.Text, DateTimeOffset.UtcNow);
}

// Interface for testability
public interface IReportGenerator<TData>
{
    Result<ReportOutput> Generate(TData data, ReportFormat format);
}

// Implementation
public sealed class SalesReportGenerator(
    IReportFormatter formatter,
    ILogger<SalesReportGenerator> logger) : IReportGenerator<SalesReportData>
{
    public Result<ReportOutput> Generate(SalesReportData data, ReportFormat format)
    {
        ArgumentNullException.ThrowIfNull(data);
        
        var validationResult = ValidateData(data);
        if (validationResult.IsFailure)
            return validationResult;
        
        if (data.IsEmpty)
            return Result<ReportOutput>.Success(ReportOutput.Empty("No sales data available"));
        
        var summary = CalculateSummary(data);
        var content = BuildContent(summary, data.ValidItems);
        var output = CreateOutput(content, format);
        
        logger.LogInformation(
            "Generated sales report: {ItemCount} items, {Total:C}", 
            summary.ItemCount, 
            summary.TotalAmount);
        
        return Result<ReportOutput>.Success(output);
    }
    
    private Result<ReportOutput> ValidateData(SalesReportData data)
    {
        var invalidItems = data.InvalidItems.ToList();
        
        if (invalidItems.Count == 0)
            return Result<ReportOutput>.Success(null!); // Will be replaced
        
        logger.LogWarning(
            "Validation failed: {Count} invalid items found", 
            invalidItems.Count);
        
        return Result<ReportOutput>.Failure(
            Error.Validation($"Found {invalidItems.Count} items with invalid amounts"));
    }
    
    private static SalesSummary CalculateSummary(SalesReportData data)
    {
        var validItems = data.ValidItems.ToList();
        var total = validItems.Sum(i => i.Amount);
        var count = validItems.Count;
        var average = count > 0 ? total / count : 0;
        
        return new SalesSummary(total, count, average);
    }
    
    private static string BuildContent(SalesSummary summary, IEnumerable<SalesItem> items)
    {
        var builder = new StringBuilder();
        
        builder.AppendLine("Sales Report");
        builder.AppendLine(new string('═', 50));
        builder.AppendLine();
        builder.AppendLine($"Total Revenue: {summary.TotalAmount:C}");
        builder.AppendLine($"Total Items:   {summary.ItemCount}");
        builder.AppendLine($"Average Sale:  {summary.AverageAmount:C}");
        builder.AppendLine();
        builder.AppendLine("Details:");
        builder.AppendLine(new string('─', 50));
        
        foreach (var item in items)
        {
            builder.AppendLine($"  • {item.Name,-30} {item.Amount,12:C}");
        }
        
        builder.AppendLine(new string('═', 50));
        
        return builder.ToString();
    }
    
    private ReportOutput CreateOutput(string content, ReportFormat format)
    {
        var formatted = formatter.Format(content, format);
        return new ReportOutput(formatted, format, DateTimeOffset.UtcNow);
    }
}

// Tests
public class SalesReportGeneratorTests
{
    private readonly SalesReportGenerator _sut;
    private readonly Mock<IReportFormatter> _formatterMock = new();
    private readonly Mock<ILogger<SalesReportGenerator>> _loggerMock = new();
    
    public SalesReportGeneratorTests()
    {
        _formatterMock
            .Setup(f => f.Format(It.IsAny<string>(), It.IsAny<ReportFormat>()))
            .Returns<string, ReportFormat>((c, _) => c);
        
        _sut = new SalesReportGenerator(_formatterMock.Object, _loggerMock.Object);
    }
    
    [Fact]
    public void Generate_WithValidData_ReturnsSuccessWithCorrectTotal()
    {
        // Arrange
        var data = new SalesReportData([
            new SalesItem("Widget", 100m),
            new SalesItem("Gadget", 200m)
        ]);
        
        // Act
        var result = _sut.Generate(data, ReportFormat.Text);
        
        // Assert
        result.IsSuccess.Should().BeTrue();
        result.Value.Content.Should().Contain("$300.00");
    }
    
    [Fact]
    public void Generate_WithEmptyData_ReturnsEmptyReport()
    {
        // Arrange
        var data = new SalesReportData([]);
        
        // Act
        var result = _sut.Generate(data, ReportFormat.Text);
        
        // Assert
        result.IsSuccess.Should().BeTrue();
        result.Value.Content.Should().Contain("No sales data");
    }
    
    [Fact]
    public void Generate_WithInvalidAmounts_ReturnsValidationError()
    {
        // Arrange
        var data = new SalesReportData([
            new SalesItem("Valid", 100m),
            new SalesItem("Invalid", -50m)
        ]);
        
        // Act
        var result = _sut.Generate(data, ReportFormat.Text);
        
        // Assert
        result.IsFailure.Should().BeTrue();
        result.Error.Code.Should().Be("VALIDATION");
    }
}
```

---

## Summary: Transformation Steps

| Step | Focus | Key Changes |
|------|-------|-------------|
| 0 → 1 | **Names** | Magic numbers → enums/constants, clear variable names |
| 1 → 2 | **Types** | `object/Dictionary` → records, LINQ, StringBuilder |
| 2 → 3 | **Functions** | One method → many small focused functions |
| 3 → 4 | **Errors** | Result pattern, validation, logging |
| 4 → 5 | **Polish** | Interfaces, DI, tests, formatting |

## Metrics Comparison

| Metric | Before | After |
|--------|--------|-------|
| Lines in main method | 80+ | 15 |
| Longest function | 80+ | 12 |
| Cyclomatic complexity | ~15 | 2-3 |
| Test coverage | 0% | 100% |
| Type safety | ❌ | ✅ |

## Related

- [Meaningful Names](../../principles/01-meaningful-names.md)
- [Functions](../../principles/02-functions.md)
- [Error Handling](../../principles/06-error-handling.md)
- [Complete Example](08-complete-example.md)
