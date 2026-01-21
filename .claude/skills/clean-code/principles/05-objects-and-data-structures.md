# Principle 5: Objects and Data Structures

> "Objects hide their data behind abstractions and expose functions that operate on that data. Data structures expose their data and have no meaningful functions."
> — Robert C. Martin

## Overview

There is a fundamental dichotomy between objects and data structures. Understanding this distinction is crucial for clean design. They are complementary but opposite approaches, each with different strengths and trade-offs.

## Why It Matters

- **Design clarity**: Understanding when to use each leads to better designs
- **Flexibility**: Choose the right approach for the problem
- **Maintainability**: Clear distinction makes code easier to modify
- **Encapsulation**: Proper use of objects protects invariants
- **Communication**: Structure reveals intent

## The Fundamental Difference

**Core distinction:**

| Aspect | Objects | Data Structures |
|--------|---------|-----------------|
| Data | Hidden (private) | Exposed (public) |
| Behavior | Rich methods | No behavior |
| Purpose | Encapsulate business logic | Transfer or hold data |
| Example | Domain entities | DTOs, Records, Config |
| Adding new types | Easy | Hard (change all functions) |
| Adding new functions | Hard (change all classes) | Easy |

**Key insight:** This isn't a quality judgment - neither is better. They solve different problems and you need both.

## Data Abstraction

### 5.1 Hide Implementation

**Principle:** Abstract the data, don't just add getters/setters. True abstraction hides whether data is computed or stored, and in what form.

**Why it matters:**
- Implementation can change without affecting clients
- Allows coordinate system changes
- Enables caching and lazy computation
- Protects invariants
- Reduces coupling

**Good abstraction characteristics:**
- Exposes operations, not data
- Implementation-agnostic interface
- Clients don't know internal representation
- Methods express intent, not mechanics

### 5.2 Getters/Setters Are Not Abstraction

**Common misconception:** Adding getters/setters makes data private.

**Reality:** Getters/setters are just syntactic sugar for public fields if they expose internal representation.

**True abstraction:**
- Hides how data is stored
- Provides operations that make sense in the problem domain
- Protects invariants
- May compute values rather than store them

## Data/Object Anti-Symmetry

### 5.3 Procedural Code (Data Structures)

**Characteristics:**
- Data structures with no behavior
- Functions operate on data structures
- Easy to add new functions without changing data structures
- Hard to add new data structures (requires changing all functions)

**When procedural is better:**
- New functions are added frequently
- Types are stable
- Algorithms need to operate on multiple types
- Performance-critical code
- Simple data transformations

**Use cases:**
- Utility functions
- Data processing pipelines
- Report generation
- Algorithm implementations

### 5.4 Object-Oriented Code (Objects)

**Characteristics:**
- Objects with behavior
- Data is hidden
- Easy to add new classes without changing existing functions
- Hard to add new functions (requires changing all classes)

**When OO is better:**
- New types are added frequently
- Behavior varies by type
- Need polymorphism
- Domain modeling
- Complex business logic

**Use cases:**
- Domain entities
- Business logic
- Plugin architectures
- Framework extension points

### 5.5 When to Use Each

**Choose objects when:**
- New types will be added frequently
- Behavior is complex and varies by type
- Need polymorphism and substitutability
- Modeling a business domain
- Protecting invariants is important
- Encapsulation benefits outweigh flexibility

**Choose data structures when:**
- New functions will be added frequently
- Types are stable and well-known
- Transferring data between layers
- Serialization/deserialization needed
- API contracts and DTOs
- Performance is critical
- Simple data aggregation

## The Law of Demeter

### 5.6 "Only Talk to Your Friends"

**Principle:** A method should only call methods on:
1. Its own class
2. Objects it creates
3. Objects passed as parameters
4. Objects held in instance variables

**Do NOT call methods on objects returned by other calls.**

**Why it matters:**
- Reduces coupling
- Makes changes easier
- Hides implementation details
- Prevents ripple effects from changes
- Improves testability

**Violation indicators:**
- Multiple dots in a chain (except fluent interfaces)
- Reaching through objects to get to data
- Knowledge of internal structure of other objects

### 5.7 Train Wrecks

**Problem:** Chains of calls that reach through multiple objects.

**Why bad:**
- Tight coupling to internal structure
- Changes in structure break many clients
- Hard to test
- Violates encapsulation
- Exposes implementation

**Fix strategies:**
1. **Tell, Don't Ask**: Give objects commands instead of querying for data
2. **Extract methods**: Add methods to intermediate objects
3. **Re-examine design**: Maybe objects don't have right responsibilities

**Exception:** Fluent interfaces (builders, LINQ) where chaining is intentional and returns same type.

### 5.8 Hybrids: The Worst of Both Worlds

**Problem:** Structures that are half object, half data structure.

**Why hybrids are bad:**
- Hard to add new functions (like objects)
- Hard to add new types (like data structures)
- Gets worst of both worlds
- Confusing intent
- Unclear responsibilities

**Symptoms:**
- Public getters/setters AND business methods
- Mix of procedural and OO style
- Unclear whether it's data or behavior
- Active Record pattern

**Fix:**
- Separate data from behavior
- Make pure objects OR pure data structures
- Use DTOs for data transfer
- Use domain objects for behavior

## Data Transfer Objects (DTOs)

### 5.9 When DTOs Are Appropriate

**Purpose:** DTOs are data structures with no behavior - and that's perfectly OK for their purpose!

**Valid DTO uses:**
- API request/response models
- Data transfer between layers
- Configuration objects
- Database query results (before mapping to domain)
- Message queue payloads
- External system integration

**DTO characteristics:**
- Public properties
- No business logic
- Simple validation OK (required, format)
- Can have computed properties for convenience
- Serializable
- Immutable preferred (use records)

### 5.10 DTO Rules

**Do:**
- Use for data transfer between layers
- Make properties public (it's a data structure!)
- Use records in C# 9+ for immutability
- Keep them simple
- Have simple computed/derived properties
- Use data annotations for validation rules

**Don't:**
- Add business logic
- Mix with domain entities
- Use for internal domain modeling
- Make them do behavior
- Share between layers if representations differ

### 5.11 Active Record Anti-Pattern

**Problem:** Data structures with save/find methods mix data and persistence.

**Why it's bad:**
- Violates Single Responsibility
- Mixes domain logic with data access
- Hard to test domain logic
- Tight coupling to database
- Difficult to change persistence strategy
- Violates Dependency Inversion

**Better approach:**
- Separate domain entities from persistence
- Use Repository pattern
- Keep domain objects persistence-ignorant
- Use separate data access layer
- Test domain logic without database

## Quick Checklist

- [ ] Clear distinction between objects and data structures
- [ ] Objects hide data, expose behavior
- [ ] Data structures expose data, have no (or minimal) behavior
- [ ] No hybrid structures (pick one!)
- [ ] Law of Demeter followed (no train wrecks)
- [ ] DTOs used only for data transfer
- [ ] Domain entities have rich, meaningful behavior
- [ ] No Active Record pattern
- [ ] Public properties only on DTOs, not domain objects
- [ ] Abstractions hide implementation details

## When to Apply

**When designing a class, ask:**
1. Is this primarily data or behavior?
2. Will I add more types or more operations?
3. Should clients know the internal structure?
4. Does this need to protect invariants?
5. Is this for data transfer or domain logic?

**Red flags:**
- Object with mostly getters/setters and little behavior
- Data structure with complex methods
- Needing to change many classes to add a function
- Needing to change many functions to add a type
- Train wreck chains of method calls
- Public fields with business logic methods

**Refactoring triggers:**
- Hybrid structures (has public data AND business logic)
- Law of Demeter violations
- Active Record pattern in domain layer
- DTOs with business logic
- Domain objects used as DTOs

## See Examples

For detailed code examples demonstrating these principles:
- [C# Objects and Data Structures Examples](../examples/csharp/05-objects-data-structures.md) - Comprehensive examples

## See Also

- [Classes](09-classes.md) - Class design and organization
- [Error Handling](06-error-handling.md) - Null Object pattern
- [Functions](02-functions.md) - Tell, Don't Ask principle
- [Clean Architecture](../../clean-architecture/SKILL.md) - Layer separation and domain modeling
