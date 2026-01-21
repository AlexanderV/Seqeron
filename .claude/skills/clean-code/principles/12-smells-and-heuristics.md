# Principle 12: Smells and Heuristics

> "A code smell is a surface indication that usually corresponds to a deeper problem in the system."
> — Martin Fowler

## Overview

This is a comprehensive catalog of code smells and heuristics from Clean Code Appendix B, organized by category. Use this as a reference when reviewing code or deciding what to refactor.

---

## Comments (C)

### C1: Inappropriate Information
Comments should not contain information better held in source control, issue tracking, or other systems.

**❌ BAD:**
```csharp
// Modified by John on 2020-03-15
// Bug fix for issue #1234
```

### C2: Obsolete Comment
Comments that are old, irrelevant, or incorrect are worse than no comments.

**❌ BAD:**
```csharp
// Returns the customer name
public Customer GetCustomer(int id) { } // Actually returns whole Customer!
```

### C3: Redundant Comment
A comment is redundant if it describes something that adequately describes itself.

**❌ BAD:**
```csharp
i++; // Increment i
```

### C4: Poorly Written Comment
Comments should be brief, precise, and grammatically correct.

### C5: Commented-Out Code
Delete it. Version control remembers everything.

**❌ BAD:**
```csharp
// var old = GetOldValue();
// if (old != null) Process(old);
```

---

## Environment (E)

### E1: Build Requires More Than One Step
Building should be a single trivial operation.

**✅ GOOD:**
```bash
dotnet build
```

### E2: Tests Require More Than One Step
Running all tests should be a single trivial operation.

**✅ GOOD:**
```bash
dotnet test
```

---

## Functions (F)

### F1: Too Many Arguments
Functions should have few arguments (0-2 ideal, 3 max).

**❌ BAD:**
```csharp
void CreateUser(string firstName, string lastName, string email,
    string phone, DateTime birthDate, string address, bool isActive);
```

**✅ GOOD:**
```csharp
void CreateUser(UserRegistration registration);
```

### F2: Output Arguments
Avoid using output arguments; return objects instead.

**❌ BAD:**
```csharp
void GetCustomer(int id, out Customer customer);
```

**✅ GOOD:**
```csharp
Customer GetCustomer(int id);
```

### F3: Flag Arguments
Boolean arguments indicate the function does more than one thing.

**❌ BAD:**
```csharp
void Save(Customer customer, bool validate);
```

**✅ GOOD:**
```csharp
void Save(Customer customer);
void SaveWithValidation(Customer customer);
```

### F4: Dead Function
Functions that are never called should be deleted.

---

## General (G)

### G1: Multiple Languages in One Source File
Avoid mixing languages (C# with embedded SQL, HTML, JavaScript, etc.)

**❌ BAD:**
```csharp
public string GetHtml()
{
    return "<html><body><script>alert('hello');</script></body></html>";
}
```

### G2: Obvious Behavior Is Unimplemented
Functions should implement the behaviors that their names suggest.

**❌ BAD:**
```csharp
public Day GetDayOfWeek(string name)
{
    // Doesn't handle "Thursday" - surprising!
}
```

### G3: Incorrect Behavior at the Boundaries
Always test boundary conditions.

### G4: Overridden Safeties
Don't disable compiler warnings, tests, or safety checks.

**❌ BAD:**
```csharp
#pragma warning disable
```

### G5: Duplication
DRY - Don't Repeat Yourself. Every piece of knowledge must have a single, unambiguous representation.

### G6: Code at Wrong Level of Abstraction
High-level concepts should not be mixed with low-level details.

**❌ BAD:**
```csharp
public abstract class Animal
{
    public abstract void Speak();
    public abstract double GetWingSpan(); // Not all animals have wings!
}
```

### G7: Base Classes Depending on Their Derivatives
Base classes should not know about their derived classes.

### G8: Too Much Information
Well-defined modules have small interfaces that allow you to do a lot with little.

### G9: Dead Code
Code that isn't executed should be deleted.

### G10: Vertical Separation
Variables and functions should be defined close to where they are used.

### G11: Inconsistency
If you do something a certain way, do all similar things the same way.

**❌ BAD:**
```csharp
GetCustomer();
FetchOrder();
RetrieveProduct();
```

**✅ GOOD:**
```csharp
GetCustomer();
GetOrder();
GetProduct();
```

### G12: Clutter
Keep code clean of useless things: unused variables, functions, comments.

### G13: Artificial Coupling
Things that don't depend on each other should not be coupled.

### G14: Feature Envy
Methods that use more features from another class than its own should move there.

**❌ BAD:**
```csharp
public decimal CalculateDiscount(Customer customer)
{
    if (customer.IsVip && customer.YearsActive > 5 && customer.OrderCount > 100)
    {
        return customer.TotalSpent * 0.1m;
    }
}
```

**✅ GOOD:**
```csharp
// Move to Customer class
public class Customer
{
    public decimal CalculateDiscount()
    {
        if (IsVip && YearsActive > 5 && OrderCount > 100)
            return TotalSpent * 0.1m;
        return 0;
    }
}
```

### G15: Selector Arguments
Avoid passing flags that select behavior; use polymorphism.

### G16: Obscured Intent
Code should be as expressive as possible.

**❌ BAD:**
```csharp
return m_vIK > 0 ? m_vIK * m_uLLM : m_vIK;
```

**✅ GOOD:**
```csharp
return inventoryCount > 0
    ? inventoryCount * unitPrice
    : inventoryCount;
```

### G17: Misplaced Responsibility
Put code where it belongs.

### G18: Inappropriate Static
Prefer non-static methods unless there's no reasonable object to invoke them on.

### G19: Use Explanatory Variables
Make calculations clear by using intermediate variables with good names.

**✅ GOOD:**
```csharp
var isEligibleForDiscount = customer.IsVip && order.Total > 100;
var discountPercentage = isEligibleForDiscount ? 0.1m : 0m;
var discountAmount = order.Total * discountPercentage;
```

### G20: Function Names Should Say What They Do

**❌ BAD:**
```csharp
public Date Add(int days);  // Adds to this object or returns new?
```

**✅ GOOD:**
```csharp
public Date PlusDays(int days);  // Clearly returns new Date
```

### G21: Understand the Algorithm
Before you write code, make sure you understand the algorithm.

### G22: Make Logical Dependencies Physical
If one module depends on another, make it explicit.

### G23: Prefer Polymorphism to If/Else or Switch/Case
Use polymorphism for type-based behavior.

### G24: Follow Standard Conventions
Follow your team's coding standards.

### G25: Replace Magic Numbers with Named Constants

**❌ BAD:**
```csharp
if (age > 18) { }
```

**✅ GOOD:**
```csharp
const int LegalAdultAge = 18;
if (age > LegalAdultAge) { }
```

### G26: Be Precise
Don't be lazy about decisions.

**❌ BAD:**
```csharp
return customers.FirstOrDefault();  // What if null?
```

**✅ GOOD:**
```csharp
return customers.FirstOrDefault()
    ?? throw new CustomerNotFoundException();
```

### G27: Structure over Convention
Enforce design decisions with structure, not just documentation.

### G28: Encapsulate Conditionals

**❌ BAD:**
```csharp
if (timer.HasExpired && !timer.IsRecurrent)
```

**✅ GOOD:**
```csharp
if (ShouldBeDeleted(timer))
```

### G29: Avoid Negative Conditionals

**❌ BAD:**
```csharp
if (!buffer.IsNotEmpty())
```

**✅ GOOD:**
```csharp
if (buffer.IsEmpty())
```

### G30: Functions Should Do One Thing

### G31: Hidden Temporal Couplings
Make temporal coupling explicit.

**❌ BAD:**
```csharp
public void Initialize()
{
    A();  // Must call in this order
    B();  // but nothing enforces it
    C();
}
```

**✅ GOOD:**
```csharp
public void Initialize()
{
    var a = A();
    var b = B(a);  // Dependencies make order explicit
    C(b);
}
```

### G32: Don't Be Arbitrary
Have a reason for the way you structure your code.

### G33: Encapsulate Boundary Conditions

**❌ BAD:**
```csharp
if (level + 1 < tags.Length)
{
    parts = Parse(level + 1);
}
```

**✅ GOOD:**
```csharp
int nextLevel = level + 1;
if (nextLevel < tags.Length)
{
    parts = Parse(nextLevel);
}
```

### G34: Functions Should Descend Only One Level of Abstraction

### G35: Keep Configurable Data at High Levels
Constants and configuration should be at the top of an application.

### G36: Avoid Transitive Navigation (Law of Demeter)

**❌ BAD:**
```csharp
order.Customer.Address.City
```

**✅ GOOD:**
```csharp
order.GetShippingCity()
```

---

## Java/C# (J) - Language-Specific

### J1: Avoid Long Import Lists by Using Wildcards
In C#: prefer targeted `using` statements.

### J2: Don't Inherit Constants
Use static classes or dedicated constant files.

### J3: Constants versus Enums
Prefer enums for related constants.

**❌ BAD:**
```csharp
public const int StatusPending = 0;
public const int StatusApproved = 1;
public const int StatusRejected = 2;
```

**✅ GOOD:**
```csharp
public enum OrderStatus
{
    Pending,
    Approved,
    Rejected
}
```

---

## Names (N)

### N1: Choose Descriptive Names

**❌ BAD:**
```csharp
int d;
void Proc();
```

**✅ GOOD:**
```csharp
int elapsedDays;
void ProcessPayment();
```

### N2: Choose Names at the Appropriate Level of Abstraction

### N3: Use Standard Nomenclature Where Possible
Use pattern names: Factory, Builder, Visitor, Strategy, etc.

### N4: Unambiguous Names
Names should not require context to understand.

**❌ BAD:**
```csharp
void DoRename();  // Rename what?
```

**✅ GOOD:**
```csharp
void RenamePageTitle();
```

### N5: Use Long Names for Long Scopes
Short names are OK for small scopes.

### N6: Avoid Encodings
No Hungarian notation, no type prefixes.

### N7: Names Should Describe Side Effects

**❌ BAD:**
```csharp
public ObjectOutputStream GetOS()
{
    if (_os == null)
        _os = new ObjectOutputStream();  // Creates if null!
    return _os;
}
```

**✅ GOOD:**
```csharp
public ObjectOutputStream GetOrCreateOS()
```

---

## Tests (T)

### T1: Insufficient Tests
Test everything that could possibly break.

### T2: Use a Coverage Tool
Coverage tools help find untested code.

### T3: Don't Skip Trivial Tests
They are easy to write and document assumptions.

### T4: An Ignored Test Is a Question about an Ambiguity
If requirements are unclear, tests can document questions.

### T5: Test Boundary Conditions
Edge cases are where bugs hide.

### T6: Exhaustively Test Near Bugs
If you find a bug, test around it thoroughly.

### T7: Patterns of Failure Are Revealing
Analyze which tests fail together to understand bugs.

### T8: Test Coverage Patterns Can Be Revealing
Look at what code isn't covered by tests.

### T9: Tests Should Be Fast
Slow tests don't get run.

---

## Quick Reference Table

| Code | Category | Smell |
|------|----------|-------|
| C1 | Comments | Inappropriate Information |
| C2 | Comments | Obsolete Comment |
| C3 | Comments | Redundant Comment |
| C4 | Comments | Poorly Written Comment |
| C5 | Comments | Commented-Out Code |
| E1 | Environment | Build Requires More Than One Step |
| E2 | Environment | Tests Require More Than One Step |
| F1 | Functions | Too Many Arguments |
| F2 | Functions | Output Arguments |
| F3 | Functions | Flag Arguments |
| F4 | Functions | Dead Function |
| G1-G36 | General | Various |
| N1-N7 | Names | Various |
| T1-T9 | Tests | Various |

---

## Refactoring Priorities

When reviewing code, prioritize fixing smells in this order:

1. **Critical:** G4 (Overridden Safeties), G9 (Dead Code), F4 (Dead Function)
2. **High:** G5 (Duplication), G14 (Feature Envy), G6 (Wrong Abstraction Level)
3. **Medium:** N1 (Descriptive Names), G11 (Inconsistency), F1 (Too Many Arguments)
4. **Low:** C1-C5 (Comments), G12 (Clutter)

## See Also

- [Functions](02-functions.md) - Function-related smells in depth
- [Classes](09-classes.md) - Class-related smells
- [Meaningful Names](01-meaningful-names.md) - Naming smells in depth
- [Comments](03-comments.md) - Comment smells in depth
