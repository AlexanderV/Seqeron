# Clean Code Review Checklist

Use this checklist to systematically review code quality. Mark items as you review them.

## 📝 Naming

- [ ] All names reveal intent
- [ ] No abbreviations or encodings (no Hungarian notation)
- [ ] Names are pronounceable and searchable
- [ ] Class names are nouns (Customer, Order, not Manager, Processor)
- [ ] Method names are verbs (GetCustomer, SaveOrder)
- [ ] Boolean names are predicates (IsActive, HasPermission, CanExecute)
- [ ] One consistent word per concept (not Get/Fetch/Retrieve for same thing)
- [ ] No magic numbers or strings - replaced with named constants
- [ ] Meaningful distinctions (not data1/data2 or ProductInfo/ProductData)
- [ ] Context is clear (AddressState, not just State)

**Common Naming Issues:**
- [ ] Variables named 'data', 'info', 'object', 'thing'
- [ ] Single letter variables outside loops
- [ ] Type information in names (stringName, intCount)
- [ ] Unnecessary prefixes (m_, _, the, my)

---

## ⚡ Functions

### Size and Structure
- [ ] Functions are small (5-20 lines max)
- [ ] Each function does ONE thing
- [ ] Functions have one level of abstraction
- [ ] Functions read top-to-bottom (stepdown rule)
- [ ] No nested functions deeper than 2 levels

### Parameters
- [ ] 0-2 parameters (ideal)
- [ ] No more than 3 parameters
- [ ] No flag/boolean parameters (use separate methods instead)
- [ ] Complex parameters extracted to objects

### Behavior
- [ ] No side effects (function does what name says)
- [ ] Command-Query Separation (functions either DO or ANSWER, not both)
- [ ] No hidden couplings or temporal dependencies

### Error Handling
- [ ] Uses exceptions instead of error codes
- [ ] Try-catch blocks extracted to separate functions
- [ ] Error handling is one thing (not mixed with business logic)

### Duplication
- [ ] No duplicate code blocks
- [ ] Similar logic extracted to shared functions
- [ ] Repeated validation consolidated

**Common Function Issues:**
- [ ] Functions with "And", "Or", "Process", "Handle" in name (doing too much)
- [ ] Functions longer than screen height
- [ ] Multiple return statements in complex logic
- [ ] Functions that both query and modify state

---

## 💬 Comments

- [ ] No redundant comments explaining obvious code
- [ ] No commented-out code (delete it, use version control)
- [ ] No journal/change log comments
- [ ] No position markers or separators
- [ ] No closing brace comments
- [ ] Comments explain WHY, not WHAT
- [ ] Complex logic extracted to named methods instead of comments

**Acceptable Comments:**
- [ ] Legal/copyright notices (if required)
- [ ] Intent explanation for non-obvious decisions
- [ ] Warning of consequences
- [ ] TODO comments (sparingly, with ticket references)
- [ ] Public API documentation (XML docs for libraries)

**Replace Comments With:**
- [ ] Descriptive function names
- [ ] Descriptive variable names
- [ ] Extract method refactoring

---

## 📐 Formatting

### Vertical Formatting
- [ ] Files are small (< 500 lines)
- [ ] Related code is close together
- [ ] Concepts separated by blank lines
- [ ] Caller above callee (top-to-bottom reading)
- [ ] Variables declared close to usage

### Horizontal Formatting
- [ ] Lines are short (< 120 characters)
- [ ] Proper indentation shows hierarchy
- [ ] Whitespace used to show relationships
- [ ] No horizontal alignment (IDE handles it)

### Consistency
- [ ] Team coding standards followed
- [ ] Consistent formatting throughout file
- [ ] EditorConfig or formatting rules applied
- [ ] No mixture of tabs and spaces

---

## 🏗️ Classes

### Size and Responsibility
- [ ] Classes are small (measured by responsibilities, not lines)
- [ ] Single Responsibility Principle - one reason to change
- [ ] High cohesion - methods use most instance variables
- [ ] Class name clearly describes single purpose

### Organization
- [ ] Public constants first
- [ ] Private fields next
- [ ] Constructors after fields
- [ ] Public methods after constructors
- [ ] Private helper methods after public methods

### Dependencies
- [ ] Minimal dependencies on other classes
- [ ] Depends on abstractions, not concretions (Dependency Inversion)
- [ ] No circular dependencies

**Common Class Issues:**
- [ ] "God" classes doing everything
- [ ] Classes with names ending in "Manager", "Handler", "Processor", "Utility"
- [ ] Classes with both high-level and low-level operations
- [ ] Low cohesion (methods don't use instance variables)

---

## 🎯 Objects vs Data Structures

- [ ] Clear distinction between objects and data structures
- [ ] Objects hide data, expose behavior
- [ ] Data structures expose data, no behavior (DTOs)
- [ ] No hybrid structures (objects with getters/setters but no real methods)
- [ ] Law of Demeter followed (no train wrecks: `a.GetB().GetC().DoSomething()`)

### For Domain Objects
- [ ] Private fields with public behavior methods
- [ ] Business logic inside objects, not in services
- [ ] Rich domain model, not anemic

### For DTOs
- [ ] Public properties, no behavior
- [ ] Used only for data transfer between layers
- [ ] Separate from domain entities

---

## ⚠️ Error Handling

- [ ] Uses exceptions, not error codes or status flags
- [ ] Specific exception types (not generic Exception)
- [ ] Exceptions provide context (error message includes relevant details)
- [ ] No returning null - throw exception or use Null Object/Option pattern
- [ ] No passing null as parameters
- [ ] Null checks with ArgumentNullException.ThrowIfNull()
- [ ] Third-party exceptions wrapped in custom exceptions
- [ ] Try-catch blocks don't obscure logic

**Exception Strategy:**
- [ ] Domain exceptions defined (CustomerNotFoundException, etc.)
- [ ] Exception hierarchy makes sense
- [ ] Validation exceptions separate from domain exceptions
- [ ] Infrastructure exceptions wrapped

---

## 🧪 Testability

- [ ] Code is testable (can write unit tests without database/network)
- [ ] Dependencies can be mocked/stubbed
- [ ] Functions are pure (no side effects) where possible
- [ ] No static methods that can't be tested
- [ ] No direct instantiation of dependencies (uses DI)

### Test Quality
- [ ] Tests exist for all public methods
- [ ] Tests are readable and well-named
- [ ] One assertion per test (or one concept)
- [ ] Tests are fast, independent, repeatable
- [ ] Tests follow Arrange-Act-Assert pattern
- [ ] No production code in test files

---

## 🔄 DRY - Don't Repeat Yourself

- [ ] No copy-pasted code blocks
- [ ] Similar algorithms consolidated
- [ ] Repeated validation extracted to methods
- [ ] Common operations in utility/helper methods
- [ ] Business rules have single source of truth

**Look for:**
- [ ] Same validation in multiple places
- [ ] Similar calculations with slight variations
- [ ] Duplicate error handling
- [ ] Repeated null checks

---

## 🎨 Code Smells

### Rigidity
- [ ] Easy to change code without breaking other parts
- [ ] No "I can't touch this" areas

### Fragility
- [ ] Changes in one place don't break distant, unrelated parts
- [ ] Tests catch regressions

### Immobility
- [ ] Code can be reused in other contexts
- [ ] No tight coupling to specific implementations

### Needless Complexity
- [ ] No speculative generality ("might need it someday")
- [ ] YAGNI principle followed
- [ ] Simple solutions preferred over clever ones

### Needless Repetition
- [ ] DRY principle followed throughout

### Opacity
- [ ] Code is self-explanatory
- [ ] Intent is clear from reading the code

---

## ⚡ Concurrency

See: [Concurrency Principle](principles/11-concurrency.md)

### Async/Await Basics
- [ ] Async methods named with `Async` suffix
- [ ] `ConfigureAwait(false)` used in library code
- [ ] No `async void` (except event handlers)
- [ ] No `.Result` or `.Wait()` blocking calls
- [ ] Cancellation tokens supported for long operations

### Thread Safety
- [ ] Shared state properly synchronized
- [ ] Immutable objects preferred for shared data
- [ ] No race conditions in concurrent code
- [ ] Lock scope as small as possible
- [ ] `ConcurrentDictionary` used instead of `Dictionary` with locks

### Common Pitfalls Avoided
- [ ] No `Task.Run` in ASP.NET Core request handlers
- [ ] ValueTask not double-awaited
- [ ] IDisposable properly disposed in async methods
- [ ] HttpClient reused (not created per request)
- [ ] Database connections not held during long async operations

---

## 🏛️ Architecture Compatibility

This section ensures clean code aligns with [Clean Architecture](../clean-architecture/SKILL.md) principles:

- [ ] Business logic in Domain layer, not scattered in services
- [ ] Domain entities are rich (not anemic)
- [ ] No infrastructure concerns in Domain layer
- [ ] Dependencies point inward (Dependency Rule)
- [ ] DTOs separate from domain entities
- [ ] Repository interfaces in Domain/Application, implementations in Infrastructure

**Cross-Reference:**
- If architectural issues found, use [Clean Architecture Checklist](../clean-architecture/CHECKLIST.md)
- For layer violations, refer to [Dependency Rule](../clean-architecture/principles/01-dependency-rule.md)
- For entity design, refer to [SOLID Principles](../clean-architecture/principles/02-solid-principles.md)
- For feature organization, refer to [Vertical Slices](../clean-architecture/principles/10-vertical-slices.md)

---

## 🔍 Review Priorities

Review in this order:

### Priority 1: Critical Issues (Fix Immediately)
- Security vulnerabilities
- Null reference exceptions
- Resource leaks
- Thread safety issues
- Breaking changes to public APIs

### Priority 2: Major Issues (Fix Soon)
- God classes (violate SRP)
- Long methods (> 30 lines)
- High cyclomatic complexity (> 10)
- Major code duplication
- Missing error handling
- Anemic domain models

### Priority 3: Medium Issues (Fix Next Sprint)
- Poor naming
- Functions with too many parameters
- Low cohesion classes
- Minor duplication
- Lack of tests

### Priority 4: Minor Issues (Refactor When Touching Code)
- Missing comments for public APIs
- Inconsistent formatting
- Minor code smells
- Optimization opportunities

---

## 📊 Metrics

For detailed metrics, tools setup, and CI integration, see:
→ [Metrics & Tools Guide](examples/csharp/11-metrics-and-tools.md)

**Quick Reference:**

| Metric | Target | Red Flag |
|--------|--------|----------|
| Cyclomatic Complexity | < 10 | > 15 |
| Method LOC | < 20 | > 30 |
| Class LOC | < 200 | > 300 |
| Code Coverage | > 80% | < 60% |
| Maintainability Index | > 70 | < 50 |

---

## 🛠️ Tools

Use these tools to automate checks:

### C# / .NET
- **Roslyn Analyzers** - Built-in code analysis
- **StyleCop** - Style and consistency
- **SonarQube / SonarLint** - Code quality and security
- **ReSharper** - Refactoring and analysis
- **Code Metrics** - Complexity analysis
- **FxCop** - Framework design guidelines

### General
- **EditorConfig** - Consistent formatting
- **Git hooks** - Pre-commit quality checks
- **CI/CD** - Automated quality gates

---

## ✅ Quick Review Template

Copy this for quick reviews:

```
### File: [FileName]

**Naming:** ✅ / ⚠️ / ❌
- Issues:

**Functions:** ✅ / ⚠️ / ❌
- Issues:

**Comments:** ✅ / ⚠️ / ❌
- Issues:

**Error Handling:** ✅ / ⚠️ / ❌
- Issues:

**DRY:** ✅ / ⚠️ / ❌
- Issues:

**Testability:** ✅ / ⚠️ / ❌
- Issues:

**Overall Score:** ✅ Good / ⚠️ Needs Work / ❌ Major Issues

**Action Items:**
1.
2.
3.
```

---

## 🎯 Red Flags

Stop and refactor immediately if you see:

- [ ] Method longer than 50 lines
- [ ] Class longer than 500 lines
- [ ] Cyclomatic complexity > 15
- [ ] More than 3 levels of nesting
- [ ] Duplicate code blocks > 5 lines
- [ ] Method with > 5 parameters
- [ ] Returning null from public methods
- [ ] Empty catch blocks
- [ ] God classes
- [ ] Commented-out code throughout

---

## 📚 References

- **Clean Code** by Robert C. Martin
- **Refactoring** by Martin Fowler
- **Code Complete** by Steve McConnell
- [Clean Architecture Skill](../clean-architecture/SKILL.md) for architectural checks
- [PRINCIPLES.md](PRINCIPLES.md) for detailed explanations

---

## 💡 Remember

> "Indeed, the ratio of time spent reading versus writing is well over 10 to 1. We are constantly reading old code as part of the effort to write new code. Making it easy to read makes it easier to write."
>
> — Robert C. Martin

**The Boy Scout Rule:** Always leave the code cleaner than you found it.
