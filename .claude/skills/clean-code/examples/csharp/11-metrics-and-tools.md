# Code Quality Metrics and Tools

> **Measuring and improving code quality in C#/.NET projects**

## Quick Start

```bash
# Install .NET tools globally
dotnet tool install -g dotnet-format
dotnet tool install -g JetBrains.ReSharper.GlobalTools

# Run code analysis
dotnet build /p:TreatWarningsAsErrors=true
dotnet format --verify-no-changes

# Check with SonarQube (if configured)
dotnet sonarscanner begin /k:"project-key"
dotnet build
dotnet sonarscanner end
```

---

## Key Code Quality Metrics

### 1. Cyclomatic Complexity

**What it measures:** Number of linearly independent paths through code.

**Target:** ≤ 10 per method, ≤ 20 for complex algorithms

**How to check:**
```xml
<!-- .editorconfig -->
dotnet_diagnostic.CA1502.severity = warning
```

```csharp
// ❌ High complexity (CC = 12)
public decimal CalculatePrice(Order order)
{
    decimal price = 0;
    if (order.Type == OrderType.Standard)
    {
        if (order.Quantity > 100) price = order.Quantity * 0.8m;
        else if (order.Quantity > 50) price = order.Quantity * 0.9m;
        else price = order.Quantity * 1.0m;
    }
    else if (order.Type == OrderType.Premium)
    {
        if (order.Customer.IsVip) price = order.Quantity * 0.7m;
        else price = order.Quantity * 0.85m;
    }
    // ... more branches
    return price;
}

// ✅ Low complexity (CC = 2 per method)
public decimal CalculatePrice(Order order)
{
    var strategy = _pricingStrategyFactory.Create(order.Type);
    return strategy.Calculate(order);
}
```

### 2. Cognitive Complexity

**What it measures:** How difficult code is to understand (more intuitive than CC).

**Target:** ≤ 15 per method

**Key factors:**
- Nesting depth (each level adds more)
- Breaks in linear flow (loops, conditionals)
- Recursion

### 3. Lines of Code (LOC)

**Targets:**
| Scope | Target | Warning |
|-------|--------|---------|
| Method | ≤ 20 | > 30 |
| Class | ≤ 200 | > 300 |
| File | ≤ 400 | > 500 |

### 4. Maintainability Index

**Scale:** 0-100 (higher is better)

| Score | Rating | Action |
|-------|--------|--------|
| 80-100 | High | Good maintainability |
| 60-79 | Moderate | Consider refactoring |
| 40-59 | Low | Plan refactoring |
| 0-39 | Very Low | Refactor immediately |

### 5. Code Coverage

**Targets:**
| Type | Minimum | Recommended |
|------|---------|-------------|
| Line coverage | 60% | 80% |
| Branch coverage | 50% | 70% |
| Mutation score | 40% | 60% |

---

## Tools Setup

### 1. Roslyn Analyzers (Built-in)

```xml
<!-- Directory.Build.props -->
<Project>
  <PropertyGroup>
    <AnalysisLevel>latest-recommended</AnalysisLevel>
    <EnableNETAnalyzers>true</EnableNETAnalyzers>
    <TreatWarningsAsErrors>true</TreatWarningsAsErrors>
    <Nullable>enable</Nullable>
  </PropertyGroup>

  <ItemGroup>
    <!-- Additional analyzers -->
    <PackageReference Include="Microsoft.CodeAnalysis.NetAnalyzers" Version="8.*" />
    <PackageReference Include="SonarAnalyzer.CSharp" Version="9.*" />
    <PackageReference Include="Roslynator.Analyzers" Version="4.*" />
    <PackageReference Include="StyleCop.Analyzers" Version="1.*" />
  </ItemGroup>
</Project>
```

### 2. EditorConfig

```ini
# .editorconfig
root = true

[*.cs]
# Indentation and spacing
indent_style = space
indent_size = 4
tab_width = 4

# New line preferences
end_of_line = crlf
insert_final_newline = true

# Organize usings
dotnet_sort_system_directives_first = true
dotnet_separate_import_directive_groups = true

# Naming conventions
dotnet_naming_rule.private_fields_should_be_camel_case.severity = error
dotnet_naming_rule.private_fields_should_be_camel_case.symbols = private_fields
dotnet_naming_rule.private_fields_should_be_camel_case.style = camel_case_underscore

dotnet_naming_symbols.private_fields.applicable_kinds = field
dotnet_naming_symbols.private_fields.applicable_accessibilities = private

dotnet_naming_style.camel_case_underscore.required_prefix = _
dotnet_naming_style.camel_case_underscore.capitalization = camel_case

# Code quality rules
dotnet_diagnostic.CA1062.severity = error   # Validate arguments
dotnet_diagnostic.CA1502.severity = warning # Cyclomatic complexity
dotnet_diagnostic.CA1505.severity = warning # Maintainability
dotnet_diagnostic.CA1506.severity = warning # Coupling

# Clean Code rules
dotnet_diagnostic.CA1707.severity = error   # No underscores in names
dotnet_diagnostic.CA1708.severity = error   # Identifiers should differ
dotnet_diagnostic.CA1710.severity = warning # Collections suffix
dotnet_diagnostic.CA1715.severity = error   # Interface prefix I
dotnet_diagnostic.CA1716.severity = warning # Reserved keywords

# Null safety
dotnet_diagnostic.CS8600.severity = error   # Null literal conversion
dotnet_diagnostic.CS8602.severity = error   # Dereference null
dotnet_diagnostic.CS8603.severity = error   # Null reference return
```

### 3. SonarQube/SonarCloud

```yaml
# sonar-project.properties
sonar.projectKey=my-project
sonar.organization=my-org
sonar.sources=src
sonar.tests=tests
sonar.exclusions=**/Migrations/**,**/obj/**
sonar.cs.opencover.reportsPaths=**/coverage.opencover.xml
sonar.cs.vstest.reportsPaths=**/*.trx

# Quality gate thresholds
sonar.qualitygate.wait=true
```

**CI Pipeline (GitHub Actions):**
```yaml
name: Code Quality

on: [push, pull_request]

jobs:
  sonarcloud:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with:
          fetch-depth: 0

      - name: Setup .NET
        uses: actions/setup-dotnet@v4
        with:
          dotnet-version: '8.0.x'

      - name: Install SonarScanner
        run: dotnet tool install --global dotnet-sonarscanner

      - name: Build and analyze
        env:
          SONAR_TOKEN: ${{ secrets.SONAR_TOKEN }}
        run: |
          dotnet sonarscanner begin \
            /k:"project-key" \
            /o:"org" \
            /d:sonar.token="${SONAR_TOKEN}" \
            /d:sonar.host.url="https://sonarcloud.io" \
            /d:sonar.cs.opencover.reportsPaths="**/coverage.opencover.xml"
          
          dotnet build --no-incremental
          dotnet test --collect:"XPlat Code Coverage" -- DataCollectionRunSettings.DataCollectors.DataCollector.Configuration.Format=opencover
          
          dotnet sonarscanner end /d:sonar.token="${SONAR_TOKEN}"
```

### 4. Code Coverage with Coverlet

```xml
<!-- Test project .csproj -->
<ItemGroup>
  <PackageReference Include="coverlet.collector" Version="6.*" />
  <PackageReference Include="coverlet.msbuild" Version="6.*" />
</ItemGroup>
```

```bash
# Generate coverage report
dotnet test /p:CollectCoverage=true \
            /p:CoverletOutputFormat=opencover \
            /p:CoverletOutput=./coverage/

# Generate HTML report
dotnet tool install -g dotnet-reportgenerator-globaltool
reportgenerator -reports:./coverage/coverage.opencover.xml \
                -targetdir:./coverage/report \
                -reporttypes:Html
```

---

## Automated Quality Gates

### Pre-Commit Hooks (Husky.NET)

```bash
# Install Husky
dotnet new tool-manifest
dotnet tool install Husky
dotnet husky install
```

```json
// .husky/task-runner.json
{
  "tasks": [
    {
      "name": "format",
      "command": "dotnet",
      "args": ["format", "--verify-no-changes"],
      "include": ["**/*.cs"]
    },
    {
      "name": "build",
      "command": "dotnet",
      "args": ["build", "-c", "Release", "/warnaserror"]
    }
  ]
}
```

### CI Quality Gate

```yaml
# .github/workflows/quality.yml
name: Quality Gate

on:
  pull_request:
    branches: [main]

jobs:
  quality:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Setup .NET
        uses: actions/setup-dotnet@v4
        with:
          dotnet-version: '8.0.x'

      - name: Restore
        run: dotnet restore

      - name: Check formatting
        run: dotnet format --verify-no-changes

      - name: Build with warnings as errors
        run: dotnet build -c Release /warnaserror

      - name: Run tests with coverage
        run: |
          dotnet test -c Release \
            --collect:"XPlat Code Coverage" \
            --results-directory ./coverage \
            /p:Threshold=80 \
            /p:ThresholdType=line

      - name: Check coverage threshold
        run: |
          # Fail if coverage below 80%
          coverage=$(grep -oP 'line-rate="\K[^"]+' coverage/**/coverage.cobertura.xml | head -1)
          if (( $(echo "$coverage < 0.80" | bc -l) )); then
            echo "Coverage $coverage is below 80% threshold"
            exit 1
          fi
```

---

## Metric Visualization

### Visual Studio Metrics

1. **Analyze** → **Calculate Code Metrics** → **For Solution**
2. Review: Maintainability Index, Cyclomatic Complexity, Depth of Inheritance, Class Coupling

### JetBrains dotCover/ReSharper

```bash
# Command line coverage
dotnet tool install -g JetBrains.dotCover.GlobalTool
dotnet dotcover test --dcOutput=coverage.html --dcReportType=HTML
```

### NDepend (Advanced)

Key metrics available:
- Technical Debt estimation
- Code coverage visualization
- Dependency graphs
- Trend analysis over time

---

## Metrics Cheat Sheet

| Metric | Good | Warning | Bad |
|--------|------|---------|-----|
| Cyclomatic Complexity | ≤ 10 | 11-20 | > 20 |
| Cognitive Complexity | ≤ 15 | 16-25 | > 25 |
| Method LOC | ≤ 20 | 21-40 | > 40 |
| Class LOC | ≤ 200 | 201-400 | > 400 |
| Maintainability Index | ≥ 80 | 60-79 | < 60 |
| Code Coverage | ≥ 80% | 60-79% | < 60% |
| Duplicated Lines | ≤ 3% | 3-5% | > 5% |
| Technical Debt Ratio | ≤ 5% | 5-10% | > 10% |

---

## Quick Commands Reference

```bash
# Format all code
dotnet format

# Check formatting without changes
dotnet format --verify-no-changes

# Build with all warnings as errors
dotnet build -c Release /warnaserror

# Run tests with coverage
dotnet test --collect:"XPlat Code Coverage"

# Generate coverage report
reportgenerator -reports:**/coverage.cobertura.xml -targetdir:coveragereport

# Run Roslyn analyzers
dotnet build /p:RunAnalyzersDuringBuild=true

# Check for outdated packages
dotnet list package --outdated
```

---

## Related Resources

- [Clean Code Principles](../../PRINCIPLES.md)
- [Code Smells Catalog](../../principles/12-smells-and-heuristics.md)
- [Refactoring Examples](05-refactoring.md)
