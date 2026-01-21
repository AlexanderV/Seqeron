# Principle 1: Meaningful Names

> "The name of a variable, function, or class, should answer all the big questions. It should tell you why it exists, what it does, and how it is used."
> — Robert C. Martin

## Overview

Naming is one of the most important aspects of clean code. Good names make code self-documenting and reduce the need for comments. Every name in your code should clearly communicate its purpose, making the code easier to read, understand, and maintain.

## Why It Matters

- **Reduces cognitive load**: Developers spend more time reading code than writing it
- **Makes code self-documenting**: Good names eliminate the need for explanatory comments
- **Prevents errors**: Clear names reduce misunderstandings and bugs
- **Improves maintainability**: Future developers (including yourself) understand intent immediately
- **Facilitates code reviews**: Reviewers can focus on logic rather than deciphering names

## Key Rules

### 1.1 Use Intention-Revealing Names

**Principle:** A name should tell you why it exists, what it does, and how it is used without requiring comments.

**Rules:**
- The name should answer: Why does this exist? What does it do? How is it used?
- If a name requires a comment to explain it, the name doesn't reveal its intent
- Be explicit rather than clever
- Don't use generic names like "data", "info", "thing", "list", "array"

### 1.2 Avoid Disinformation

**Principle:** Don't use names that vary in small ways or mean something different than intended.

**Rules:**
- Don't refer to a grouping as "List" unless it's actually a List type
- Avoid names that vary only slightly (e.g., XYZControllerForEfficientHandlingOfStrings vs XYZControllerForEfficientStorageOfStrings)
- Don't use names that look like platform names or variable types
- Avoid lowercase 'l' and uppercase 'O' as they look like 1 and 0

### 1.3 Make Meaningful Distinctions

**Principle:** Don't add noise words or numbers to make names unique.

**Rules:**
- Avoid number-series naming (a1, a2, ... aN)
- Eliminate noise words: Info, Data, Manager, Processor (unless they add real meaning)
- Don't use synonyms without clear distinction (ProductInfo vs ProductData)
- Make each name meaningfully different from every other name

### 1.4 Use Pronounceable Names

**Principle:** Names should be easy to pronounce and discuss.

**Rules:**
- If you can't pronounce it, you can't discuss it intelligently
- Avoid cryptic abbreviations (genymdhms for generation year-month-day-hour-minute-second)
- Use complete words rather than abbreviations when possible
- This makes code reviews and discussions much easier

### 1.5 Use Searchable Names

**Principle:** Single-letter names and numeric constants are hard to locate across a codebase.

**Rules:**
- Avoid single-letter names except for local loop variables in small scopes
- Replace magic numbers with named constants
- The length of a name should correspond to the size of its scope
- For variables used in multiple places, use meaningful, searchable names

### 1.6 Avoid Encodings

**Principle:** Don't encode type or scope information into names.

**Rules:**
- No Hungarian notation (strName, iCount, bFlag)
- No member prefixes (m_description, _name)
- Modern IDEs provide type information without needing it in names
- Interfaces: Prefer ICustomer over CustomerInterface (language-dependent)
- Don't encode scope (g_variable for global, l_variable for local)

### 1.7 Class Names

**Principle:** Classes should have noun or noun phrase names.

**Rules:**
- Use nouns: Customer, WikiPage, Account, AddressParser
- Avoid verbs: Manager, Processor, Data (unless they have specific meaning)
- Don't use generic names: Data, Info, Object
- Be specific: CustomerAccount instead of Account if context matters

### 1.8 Method Names

**Principle:** Methods should have verb or verb phrase names.

**Rules:**
- Use verbs: Save, Delete, Create, Update
- Accessors/mutators/predicates: Get, Set, Is, Has, Can
- Factory methods: CreateFrom, BuildWith, NewInstance
- Boolean methods: IsValid, HasPermission, CanExecute, ShouldRetry
- Follow language conventions (C# properties vs Get/Set methods)

### 1.9 Pick One Word per Concept

**Principle:** Be consistent across your codebase - use the same word for the same abstract concept.

**Rules:**
- Choose one word per concept and stick with it throughout the codebase
- Don't mix: Fetch, Retrieve, Get for the same operation
- Don't mix: Controller, Manager, Driver for the same concept
- Create a project lexicon if needed
- Consistency reduces cognitive load

### 1.10 Don't Pun

**Principle:** Avoid using the same word for two different purposes.

**Rules:**
- Don't use "Add" for both mathematical addition and collection insertion
- Don't use "Get" for both simple retrieval and computed values
- Use "Insert" or "Append" instead of "Add" when adding to collections if "Add" means arithmetic
- Keep the semantics of each word consistent

### 1.11 Use Solution Domain Names

**Principle:** Use computer science terms, algorithm names, pattern names, math terms.

**Rules:**
- Use technical terms when appropriate: Queue, Stack, Factory, Visitor
- Algorithm names: QuickSort, BinarySearch
- Pattern names: Singleton, Observer, Strategy
- Developers will understand these terms
- This is preferred over problem domain names for technical concepts

### 1.12 Use Problem Domain Names

**Principle:** When solution domain names don't apply, use names from the problem domain.

**Rules:**
- Use terms that domain experts would recognize
- Medical system: Patient, Diagnosis, Treatment
- Financial system: Account, Transaction, Balance
- E-commerce: Cart, Order, Inventory
- At least maintainers can ask domain experts what it means

### 1.13 Add Meaningful Context

**Principle:** When a variable name needs context, add prefixes or encapsulate in a class.

**Rules:**
- Add context through prefixes when needed: addressState, orderTotal
- Better: Encapsulate related variables in a class (Address with State property)
- Don't add gratuitous context (prefix everything with "Application")
- Shorter names are better, but only if they're clear
- Context should be just enough, not more

### 1.14 Don't Add Gratuitous Context

**Principle:** Don't add needless context that makes names unnecessarily long.

**Rules:**
- Don't prefix every class with the application name
- AccountAddress is better than GSDAccountAddress (Gas Station Deluxe)
- Only add context when it adds clarity
- Keep names as short as possible while remaining clear

## Quick Checklist

- [ ] Names reveal intent without comments
- [ ] No abbreviations or encodings
- [ ] Names are pronounceable and searchable
- [ ] Class names are nouns
- [ ] Method names are verbs
- [ ] Boolean names are predicates (IsActive, HasPermission, CanExecute)
- [ ] One consistent word per concept
- [ ] No magic numbers or strings - use named constants
- [ ] Context is clear but not gratuitous
- [ ] No disinformation or misleading names

## Common Naming Anti-Patterns

| Anti-Pattern | Description | Fix |
|--------------|-------------|-----|
| Single letter | `d`, `x`, `i` (outside loop scope) | Use descriptive names: `elapsedDays`, `xCoordinate` |
| Noise words | `ProductInfo`, `ProductData`, `TheProduct` | Be specific: `Product`, `ProductSpecification` |
| Hungarian notation | `strName`, `iCount`, `bFlag` | Remove prefixes: `name`, `count`, `isEnabled` |
| Meaningless | `data`, `info`, `thing`, `stuff` | Use specific name for context |
| Inconsistent verbs | Mix of `Get`, `Fetch`, `Retrieve` | Pick one verb and use consistently |
| Encoded types | `customerList` (but it's an array) | `customers` or `customerArray` |
| Ambiguous | `process()`, `handle()`, `manage()` | Be specific: `validateOrder()`, `formatAddress()` |

## When to Apply

**Always apply these rules when:**
- Creating new variables, functions, or classes
- Refactoring existing code
- During code reviews
- When you find yourself writing a comment to explain a name

**Pay extra attention when:**
- Names will be used in multiple places (high scope)
- Names represent important domain concepts
- Names represent public APIs or interfaces
- Working in a team (consistency matters more)

## See Examples

For detailed code examples demonstrating these principles:
- [C# Naming Examples](../examples/csharp/01-naming.md) - Comprehensive examples with before/after comparisons

## See Also

- [Functions](02-functions.md) - Naming conventions for functions
- [Classes](09-classes.md) - Naming conventions for classes
- [Comments](03-comments.md) - How good names reduce need for comments
