# Principle 2: Functions

> "Functions should do one thing. They should do it well. They should do it only."
> — Robert C. Martin

## Overview

Functions are the first line of organization in any program. Making them clean is essential for readable, maintainable code. Well-crafted functions make code easier to read, understand, and modify.

## Why It Matters

- **Improves readability**: Small, focused functions are easier to understand
- **Enhances testability**: Single-purpose functions are easier to test
- **Reduces bugs**: Less complexity means fewer places for bugs to hide
- **Facilitates reuse**: Functions that do one thing can be reused in multiple contexts
- **Simplifies debugging**: Smaller functions make it easier to isolate problems
- **Enables better naming**: Functions that do one thing are easier to name clearly

## Key Rules

### 2.1 Small!

**Principle:** Functions should be small. Really small. 5-20 lines is ideal.

**Why:**
- Easier to understand at a glance
- Easier to test
- Easier to debug
- Less likely to have hidden bugs
- Encourages proper abstraction

**Rules:**
- Target 5-20 lines per function
- If a function is longer, look for extraction opportunities
- Blocks within if, else, while statements should typically be one line (a function call)
- Indent level should not be greater than one or two

### 2.2 Do One Thing

**Principle:** A function should do one thing, do it well, and do it only.

**How to tell if a function does one thing:**
- Can you extract another function from it with a name that is not merely a restatement of its implementation?
- Are all statements in the function at the same level of abstraction?
- Can you describe what the function does without using the word "and" or "or"?

**Rules:**
- A function should have one reason to change
- If you can extract a meaningful function from another, the original was doing more than one thing
- Sections within functions indicate it's doing multiple things

### 2.3 One Level of Abstraction per Function

**Principle:** Don't mix high-level and low-level operations in the same function.

**Why:**
- Mixing levels makes code confusing
- Hard to tell what's essential vs. detail
- Encourages more detail mixing as code evolves

**Rules:**
- High-level: Business operations (ProcessOrder, GenerateReport)
- Medium-level: Algorithm steps (ValidateData, CalculateTotal)
- Low-level: Implementation details (StringManipulation, FileIO)
- Keep all statements in a function at the same level

### 2.4 Reading Code from Top to Bottom (Stepdown Rule)

**Principle:** Code should read like a top-down narrative. Each function should be followed by those at the next level of abstraction.

**Structure:**
- To do X, we do A, B, and C
- To do A, we do A1 and A2
- To do A1, we do...

**Rules:**
- High-level functions call medium-level functions
- Medium-level functions call low-level functions
- Read from top to bottom with decreasing abstraction levels
- This creates a "stepdown" effect

### 2.5 Switch Statements

**Principle:** Avoid switch statements. Use polymorphism instead.

**Why switch statements are problematic:**
- They violate Single Responsibility (does N things)
- They violate Open/Closed (must modify to add cases)
- They duplicate (same structure repeated elsewhere)
- Hard to test all branches

**Rules:**
- Bury switches in low-level classes, never repeated
- Use polymorphism to hide switches
- Factory pattern for creating polymorphic objects
- Consider Strategy or State pattern
- If switch is unavoidable, make it create polymorphic objects

### 2.6 Use Descriptive Names

**Principle:** Long descriptive names are better than short enigmatic names or long descriptive comments.

**Rules:**
- Don't be afraid of long names (15-20 chars is fine)
- Name should describe what the function does
- Use verbs or verb phrases
- Be consistent with naming patterns
- Spend time choosing good names

**Examples of good patterns:**
- includeSetupAndTeardownPages
- isTestable
- validateOrderDetails
- processPaymentTransaction

### 2.7 Function Arguments

**Principle:** Minimize the number of function arguments. Zero is ideal, one or two is good, three should be avoided, more than three requires special justification.

**Argument count guidelines:**

| Count | Name | Recommendation | Reason |
|-------|------|----------------|--------|
| 0 | Niladic | Ideal | Easiest to understand and test |
| 1 | Monadic | Good | Clear purpose, easy to understand |
| 2 | Dyadic | Acceptable | More mental effort required |
| 3 | Triadic | Avoid | Significantly harder to understand |
| 3+ | Polyadic | Strongly avoid | Use parameter objects |

**When you need many arguments:**
- Create a parameter object or class
- Group related parameters into objects
- Consider if the function is doing too much

**Common monadic forms:**
- Asking a question about the argument: fileExists(file)
- Transforming the argument: fileOpen(name)
- Event: passwordAttemptFailedNTimes(attempts)

### 2.8 Flag Arguments Are Ugly

**Principle:** Don't pass boolean flags to functions - it loudly proclaims the function does more than one thing.

**Why flags are bad:**
- They complicate the signature
- They violate Single Responsibility
- They make the function do at least two things (one for true, one for false)

**Rules:**
- Split the function into two functions
- Name each function after what it does
- If flags are unavoidable, use enums instead of booleans for clarity

### 2.9 No Side Effects

**Principle:** Functions should not have hidden side effects. They should do what their name says and nothing else.

**Common side effects:**
- Modifying global state
- Modifying passed arguments
- Performing I/O operations not indicated by name
- Starting sessions or connections
- Modifying class member variables unexpectedly

**Rules:**
- A function should either change the state of an object, or return information about it, but not both
- Make side effects explicit in the function name if unavoidable
- Prefer pure functions when possible
- If a function changes state, make it obvious from the name

### 2.10 Command Query Separation (CQS)

**Principle:** Functions should either do something (command) or answer something (query), but not both.

**Commands:**
- Change the state of the system
- Return void (or status for success)
- Examples: Save, Delete, Update, Send

**Queries:**
- Return information
- Don't change state
- Examples: Get, Find, Calculate, Is, Has

**Rules:**
- Separate commands from queries
- Don't return values from commands (except success status)
- Don't change state in queries
- This makes code more predictable and testable

### 2.11 Prefer Exceptions to Returning Error Codes

**Principle:** Use exceptions for error handling rather than returning error codes.

**Why error codes are bad:**
- Lead to deeply nested structures
- Caller must handle errors immediately
- Violate Command Query Separation
- Error codes can be ignored

**Why exceptions are better:**
- Separate error handling from happy path
- Can't be silently ignored
- Can be handled at appropriate level
- Make happy path clearer

**Rules:**
- Use exceptions for exceptional conditions
- Use return values for normal control flow
- Create specific exception types
- Document exceptions in function comments

### 2.12 Extract Try/Catch Blocks

**Principle:** Error handling is one thing. Functions that handle errors should do nothing else.

**Why:**
- Try/catch blocks are ugly and confuse structure
- Error handling is "one thing"
- Separates error handling from business logic

**Rules:**
- Extract the body of try block into its own function
- Extract the body of catch block into its own function
- The try/catch function does only error handling
- The extracted function does the actual work

### 2.13 Don't Repeat Yourself (DRY)

**Principle:** Duplication is the root of all evil in software. Every piece of knowledge should have a single representation.

**Why duplication is bad:**
- Maintenance nightmare (change in N places)
- Increases likelihood of bugs
- Makes code harder to understand
- Increases size unnecessarily

**How to avoid duplication:**
- Extract common code into functions
- Use inheritance or composition
- Create utility functions
- Use configuration over code
- Apply design patterns (Template Method, Strategy)

**Rules:**
- If you copy-paste code, you're doing it wrong
- Look for similar code and extract commonality
- Three strikes rule: On third duplication, refactor

### 2.14 Structured Programming

**Principle:** Every function should have one entry and one exit (Dijkstra's rules).

**Rules (for large functions):**
- One return statement
- No break or continue in loops
- No goto statements (ever)

**Modern interpretation:**
- For small functions, multiple returns are often clearer
- Early returns can avoid deep nesting
- Break/continue in loops can be readable if function is small
- The benefit of structure increases with function size

## Quick Checklist

- [ ] Functions are small (5-20 lines)
- [ ] Each function does ONE thing only
- [ ] One level of abstraction per function
- [ ] Descriptive, intention-revealing names
- [ ] 0-2 parameters (max 3, use parameter objects for more)
- [ ] No flag/boolean parameters
- [ ] No side effects or hidden behaviors
- [ ] Command-Query Separation followed
- [ ] Uses exceptions, not error codes
- [ ] Try/catch blocks extracted to separate functions
- [ ] No duplicate code blocks (DRY)
- [ ] Functions read top-to-bottom at descending levels

## When to Apply

**Always apply these rules when:**
- Writing new functions
- Refactoring existing code
- Code reviews
- When functions become difficult to test

**Prioritize refactoring when:**
- Function exceeds 20 lines
- Function name doesn't clearly describe what it does
- Function has more than 3 parameters
- You need comments to explain what the function does
- Function mixes abstraction levels
- You find duplicated code

## See Examples

For detailed code examples demonstrating these principles:
- [C# Functions Examples](../examples/csharp/02-functions.md) - Comprehensive function refactoring examples

## See Also

- [Meaningful Names](01-meaningful-names.md) - Naming conventions for functions
- [Comments](03-comments.md) - When functions need comments
- [Error Handling](06-error-handling.md) - Exception handling in depth
- [Classes](09-classes.md) - Similar principles for classes
- [Emergence](10-emergence.md) - DRY and simple design
