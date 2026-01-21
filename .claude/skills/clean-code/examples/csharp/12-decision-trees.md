# Clean Code Decision Trees

> **Visual guides for making clean code decisions**

## 1. Naming Decision Tree

```
                    ┌─────────────────────┐
                    │  Need to name       │
                    │  something?         │
                    └──────────┬──────────┘
                               │
            ┌──────────────────┼──────────────────┐
            ▼                  ▼                  ▼
       ┌─────────┐        ┌─────────┐        ┌─────────┐
       │ Variable│        │ Method  │        │  Class  │
       └────┬────┘        └────┬────┘        └────┬────┘
            │                  │                  │
            ▼                  ▼                  ▼
    ┌───────────────┐  ┌───────────────┐  ┌───────────────┐
    │ Use noun      │  │ Use verb      │  │ Use noun      │
    │ (what it is)  │  │ (what it does)│  │ (what it      │
    │               │  │               │  │  represents)  │
    └───────┬───────┘  └───────┬───────┘  └───────┬───────┘
            │                  │                  │
            ▼                  ▼                  ▼
    ┌───────────────┐  ┌───────────────┐  ┌───────────────┐
    │ Is it bool?   │  │ Returns data? │  │ Single        │
    └───┬───────┬───┘  └───┬───────┬───┘  │ responsibility│
      Yes       No       Yes       No     └───────┬───────┘
        │       │         │        │              │
        ▼       ▼         ▼        ▼              ▼
    ┌───────┐ ┌───────┐ ┌───────┐ ┌───────┐  ┌─────────┐
    │is/has/│ │Describe│ │Get/   │ │Action │  │Name the │
    │can/   │ │content │ │Find/  │ │verb:  │  │concept: │
    │should │ │        │ │Calc   │ │Save/  │  │Customer │
    └───────┘ └───────┘ └───────┘ │Delete │  │Service  │
                                  └───────┘  └─────────┘

    Examples:
    ─────────
    bool isActive         GetCustomer()        CustomerService
    bool hasPermission    FindByEmail()        OrderProcessor
    bool canExecute       CalculateTotal()     PaymentGateway
    decimal totalAmount   Save(), Delete()     EmailValidator
```

## 2. Function Size Decision Tree

```
                    ┌─────────────────────┐
                    │  How long is your   │
                    │  function?          │
                    └──────────┬──────────┘
                               │
        ┌──────────────────────┼──────────────────────┐
        ▼                      ▼                      ▼
   ┌─────────┐            ┌─────────┐            ┌─────────┐
   │ ≤10 LOC │            │ 11-30   │            │  >30    │
   │         │            │  LOC    │            │  LOC    │
   └────┬────┘            └────┬────┘            └────┬────┘
        │                      │                      │
        ▼                      ▼                      ▼
   ┌─────────┐            ┌─────────┐            ┌─────────┐
   │  ✅ OK  │            │ Review: │            │❌ REFACTOR│
   │         │            │ Can be  │            │  NOW!   │
   │         │            │ split?  │            │         │
   └─────────┘            └────┬────┘            └────┬────┘
                               │                      │
                    ┌──────────┴──────────┐          │
                    ▼                     ▼          │
               ┌─────────┐           ┌─────────┐     │
               │  Yes    │           │   No    │     │
               │         │           │ (rare)  │     │
               └────┬────┘           └────┬────┘     │
                    │                     │          │
                    ▼                     ▼          │
               ┌─────────┐           ┌─────────┐     │
               │ Extract │           │Document │     │
               │ Method  │           │ why     │     │
               └─────────┘           └─────────┘     │
                                                     │
                    ┌────────────────────────────────┘
                    ▼
        ┌─────────────────────────────────────┐
        │  Refactoring Steps:                 │
        │  1. Identify distinct operations    │
        │  2. Extract each to named method    │
        │  3. Main function becomes outline   │
        │  4. Each helper ≤10 lines           │
        └─────────────────────────────────────┘
```

## 3. Error Handling Decision Tree

```
                    ┌─────────────────────┐
                    │  Operation can fail │
                    └──────────┬──────────┘
                               │
        ┌──────────────────────┼──────────────────────┐
        ▼                      ▼                      ▼
   ┌──────────┐           ┌──────────┐           ┌──────────┐
   │Programmer│           │ Expected │           │Unexpected│
   │  Error   │           │ Business │           │ System   │
   │(bug)     │           │  Case    │           │  Error   │
   └────┬─────┘           └────┬─────┘           └────┬─────┘
        │                      │                      │
        ▼                      ▼                      ▼
   ┌──────────┐           ┌──────────┐           ┌──────────┐
   │ArgumentNul│          │ Result<T>│           │ Exception│
   │Exception, │          │ OneOf<>  │           │ (let it  │
   │InvalidOp  │          │ Either   │           │  bubble) │
   └──────────┘           └────┬─────┘           └──────────┘
                               │
            ┌──────────────────┼──────────────────┐
            ▼                  ▼                  ▼
       ┌─────────┐        ┌─────────┐        ┌─────────┐
       │ Simple  │        │Multiple │        │ Railway │
       │ Success/│        │ Error   │        │ Oriented│
       │ Failure │        │ Types   │        │ Program.│
       └────┬────┘        └────┬────┘        └────┬────┘
            │                  │                  │
            ▼                  ▼                  ▼
       ┌─────────┐        ┌─────────┐        ┌─────────┐
       │Result<T>│        │ OneOf   │        │ Then()  │
       │.Success │        │ <T,E1,  │        │ chain   │
       │.Failure │        │  E2,E3> │        │ methods │
       └─────────┘        └─────────┘        └─────────┘


    Code Examples:
    ──────────────

    // Result pattern (simple)
    public Result<Customer> GetCustomer(CustomerId id)
    {
        var customer = _repo.Find(id);
        return customer is null
            ? Result.Failure<Customer>("Not found")
            : Result.Success(customer);
    }

    // OneOf (multiple outcomes)
    public OneOf<Success, NotFound, ValidationError> CreateOrder(...)
    {
        if (!IsValid(request)) 
            return new ValidationError(errors);
        // ...
    }
```

## 4. Class Design Decision Tree

```
                    ┌─────────────────────┐
                    │  Creating a class?  │
                    └──────────┬──────────┘
                               │
                               ▼
                    ┌─────────────────────┐
                    │ Does it have single │
                    │ responsibility?     │
                    └──────────┬──────────┘
                               │
              ┌────────────────┴────────────────┐
              ▼                                 ▼
         ┌─────────┐                       ┌─────────┐
         │   Yes   │                       │   No    │
         └────┬────┘                       └────┬────┘
              │                                 │
              ▼                                 ▼
    ┌─────────────────┐              ┌─────────────────────┐
    │ Is it cohesive? │              │ Split into multiple │
    │ (all members    │              │ classes, each with  │
    │  work together) │              │ one responsibility  │
    └────────┬────────┘              └─────────────────────┘
             │
   ┌─────────┴─────────┐
   ▼                   ▼
┌──────┐           ┌──────┐
│ Yes  │           │  No  │
└──┬───┘           └──┬───┘
   │                  │
   ▼                  ▼
┌──────────┐    ┌────────────────┐
│    ✅    │    │ Extract class  │
│   GOOD   │    │ for unrelated  │
│          │    │ members        │
└──────────┘    └────────────────┘


    Cohesion Check:
    ───────────────
    
    ❌ LOW COHESION               ✅ HIGH COHESION
    
    class UserManager {           class UserAuthenticator {
      void Login()                  void Login()
      void Register()               void Logout()
      void SendEmail()              bool ValidateToken()
      void GenerateReport()       }
      void BackupData()           
    }                             class UserRegistration {
                                    void Register()
                                    void VerifyEmail()
                                  }
```

## 5. Refactoring Priority Decision Tree

```
                    ┌─────────────────────┐
                    │  Found code smell?  │
                    └──────────┬──────────┘
                               │
                               ▼
                    ┌─────────────────────┐
                    │ Is it blocking      │
                    │ current work?       │
                    └──────────┬──────────┘
                               │
              ┌────────────────┴────────────────┐
              ▼                                 ▼
         ┌─────────┐                       ┌─────────┐
         │   Yes   │                       │   No    │
         └────┬────┘                       └────┬────┘
              │                                 │
              ▼                                 ▼
    ┌─────────────────┐              ┌─────────────────┐
    │ FIX NOW         │              │ How severe?     │
    │ (Boy Scout Rule)│              │                 │
    └─────────────────┘              └────────┬────────┘
                                              │
                     ┌────────────────────────┼────────────────────────┐
                     ▼                        ▼                        ▼
                ┌─────────┐              ┌─────────┐              ┌─────────┐
                │Critical │              │ Medium  │              │  Minor  │
                │(security│              │(readabil│              │ (style) │
                │ bugs)   │              │  ity)   │              │         │
                └────┬────┘              └────┬────┘              └────┬────┘
                     │                        │                        │
                     ▼                        ▼                        ▼
                ┌─────────┐              ┌─────────┐              ┌─────────┐
                │ Create  │              │ Add to  │              │ Fix if  │
                │ urgent  │              │ backlog │              │ touching│
                │ ticket  │              │         │              │ the file│
                └─────────┘              └─────────┘              └─────────┘


    Severity Examples:
    ──────────────────
    
    CRITICAL                MEDIUM                 MINOR
    ────────                ──────                 ─────
    • Security holes        • Long methods         • Naming style
    • Data corruption       • Deep nesting         • Extra whitespace
    • Memory leaks          • Duplicate code       • Import order
    • Race conditions       • Magic numbers        • Comment format
```

## 6. Testing Strategy Decision Tree

```
                    ┌─────────────────────┐
                    │  What to test?      │
                    └──────────┬──────────┘
                               │
        ┌──────────────────────┼──────────────────────┐
        ▼                      ▼                      ▼
   ┌──────────┐           ┌──────────┐           ┌──────────┐
   │ Business │           │Integration│          │   UI/    │
   │  Logic   │           │  Points   │          │   E2E    │
   └────┬─────┘           └────┬─────┘           └────┬─────┘
        │                      │                      │
        ▼                      ▼                      ▼
   ┌──────────┐           ┌──────────┐           ┌──────────┐
   │UNIT TESTS│           │INTEGRATION│          │ E2E TESTS│
   │  (many)  │           │  TESTS    │          │  (few)   │
   │   Fast   │           │ (moderate)│          │  Slow    │
   └──────────┘           └──────────┘           └──────────┘


         Testing Pyramid:
         ────────────────
         
              /\            E2E (5%)
             /  \           - Critical user journeys
            /    \          - Smoke tests
           /──────\
          /        \        Integration (20%)
         /          \       - API endpoints
        /            \      - Database operations
       /──────────────\     - External services
      /                \
     /                  \   Unit Tests (75%)
    /                    \  - Business rules
   /______________________\ - Pure functions
                           - Domain logic


    Test Type Selection:
    ───────────────────
    
    ┌─────────────────────┐
    │ Does it have        │
    │ external deps?      │
    └──────────┬──────────┘
               │
      ┌────────┴────────┐
      ▼                 ▼
    ┌───┐             ┌───┐
    │Yes│             │No │
    └─┬─┘             └─┬─┘
      │                 │
      ▼                 ▼
    ┌─────────────┐   ┌─────────────┐
    │ Can mock?   │   │ UNIT TEST   │
    └──────┬──────┘   │ (pure logic)│
           │          └─────────────┘
      ┌────┴────┐
      ▼         ▼
    ┌───┐     ┌───┐
    │Yes│     │No │
    └─┬─┘     └─┬─┘
      │         │
      ▼         ▼
   ┌───────┐ ┌────────────┐
   │ UNIT  │ │ INTEGRATION│
   │ TEST  │ │ TEST       │
   │ +mock │ │ (real deps)│
   └───────┘ └────────────┘
```

## 7. When to Comment Decision Tree

```
                    ┌─────────────────────┐
                    │  Should I add a     │
                    │  comment?           │
                    └──────────┬──────────┘
                               │
                               ▼
                    ┌─────────────────────┐
                    │ Can I express it    │
                    │ in code instead?    │
                    └──────────┬──────────┘
                               │
              ┌────────────────┴────────────────┐
              ▼                                 ▼
         ┌─────────┐                       ┌─────────┐
         │   Yes   │                       │   No    │
         └────┬────┘                       └────┬────┘
              │                                 │
              ▼                                 ▼
    ┌─────────────────┐              ┌─────────────────┐
    │ DO THAT INSTEAD │              │ What type?      │
    │                 │              │                 │
    │ • Better names  │              └────────┬────────┘
    │ • Extract method│                       │
    │ • Add type      │    ┌──────────────────┼──────────────────┐
    └─────────────────┘    ▼                  ▼                  ▼
                      ┌─────────┐        ┌─────────┐        ┌─────────┐
                      │  WHY    │        │  WHAT   │        │  HOW    │
                      │(intent) │        │(summary)│        │(explain)│
                      └────┬────┘        └────┬────┘        └────┬────┘
                           │                  │                  │
                           ▼                  ▼                  ▼
                      ┌─────────┐        ┌─────────┐        ┌─────────┐
                      │    ✅   │        │ Usually │        │    ❌   │
                      │  GOOD   │        │   OK    │        │  AVOID  │
                      │         │        │         │        │(code    │
                      │ // Due  │        │ // XML  │        │ should  │
                      │ // to   │        │ // docs │        │ speak)  │
                      │ // RFC  │        │ /// <sum│        │         │
                      │ // 1234 │        │ /// mary│        │         │
                      └─────────┘        └─────────┘        └─────────┘


    Good vs Bad Comments:
    ─────────────────────
    
    ✅ GOOD (WHY)                    ❌ BAD (WHAT)
    ─────────────                    ────────────
    // Compensate for legacy         // Loop through items
    // API's 1-based indexing        foreach (var item in items)
    index = apiIndex - 1;
                                     // Increment counter
    // Timeout required per          counter++;
    // SLA agreement (RFC-2023-45)
    await Task.Delay(500);           // Check if null
                                     if (x == null)
```

## 8. Extract Method Decision Tree

```
                    ┌─────────────────────┐
                    │ Should I extract    │
                    │ this code block?    │
                    └──────────┬──────────┘
                               │
                               ▼
                    ┌─────────────────────┐
                    │ Is it:              │
                    │ • Reused elsewhere? │
                    │ • Has clear name?   │
                    │ • Different abstrac.│
                    └──────────┬──────────┘
                               │
              ┌────────────────┴────────────────┐
              ▼                                 ▼
    ┌─────────────────┐              ┌─────────────────┐
    │  ANY = YES      │              │  ALL = NO       │
    └────────┬────────┘              └────────┬────────┘
             │                                │
             ▼                                ▼
    ┌─────────────────┐              ┌─────────────────┐
    │ EXTRACT METHOD  │              │ Keep inline     │
    │                 │              │ (but consider   │
    │ Name should     │              │  other reasons) │
    │ describe WHAT   │              └─────────────────┘
    │ not HOW         │
    └─────────────────┘


    Extract Method Checklist:
    ─────────────────────────
    
    Before:
    ┌─────────────────────────────────────┐
    │ public void Process()               │
    │ {                                   │
    │     // Validate input               │
    │     if (x == null) throw ...;       │
    │     if (y < 0) throw ...;           │
    │                                     │
    │     // Calculate result             │
    │     var a = x * 2;                  │
    │     var b = y + 10;                 │
    │     var result = a + b;             │
    │                                     │
    │     // Save to database             │
    │     _db.Save(result);               │
    │ }                                   │
    └─────────────────────────────────────┘
    
    After:
    ┌─────────────────────────────────────┐
    │ public void Process()               │
    │ {                                   │
    │     ValidateInput();                │
    │     var result = CalculateResult(); │
    │     SaveResult(result);             │
    │ }                                   │
    └─────────────────────────────────────┘
```

---

## Quick Reference Card

```
┌─────────────────────────────────────────────────────────────┐
│                 CLEAN CODE QUICK DECISIONS                   │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  NAMING           │  FUNCTIONS        │  ERROR HANDLING     │
│  ────────         │  ─────────        │  ──────────────     │
│  Class → Noun     │  ≤20 LOC          │  Expected → Result  │
│  Method → Verb    │  One thing        │  Bug → Exception    │
│  Bool → is/has    │  ≤3 params        │  System → Throw     │
│  Const → MEANING  │  No side effects  │                     │
│                                                             │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  CLASSES          │  COMMENTS         │  REFACTORING        │
│  ───────          │  ────────         │  ───────────        │
│  Single purpose   │  WHY → Good       │  Tests first        │
│  High cohesion    │  WHAT → Maybe     │  Small steps        │
│  ≤200 LOC         │  HOW → Avoid      │  Commit often       │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

---

## Related Resources

- [Clean Code Principles](../../PRINCIPLES.md)
- [Refactoring Examples](05-refactoring.md)
- [Error Handling](03-error-handling.md)
- [Testing Examples](06-testing.md)
