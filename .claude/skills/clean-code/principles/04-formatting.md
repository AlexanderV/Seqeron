# Principle 4: Formatting

> "Code formatting is about communication, and communication is the professional developer's first order of business."
> — Robert C. Martin

## Overview

Your style and discipline survives even after your code has changed. The formatting you set today influences the readability of code for years to come. Good formatting makes code easier to read, understand, and maintain.

## Why It Matters

- **First impressions count**: Poorly formatted code suggests carelessness
- **Long-term impact**: Formatting decisions affect readability for years
- **Team productivity**: Consistent formatting reduces cognitive load
- **Code reviews**: Well-formatted code is easier to review
- **Maintainability**: Clear structure makes code easier to modify

## Vertical Formatting

### 4.1 File Size

**Principle:** Small files are easier to understand than large files.

**Guidelines:**

| Metric | Target | Maximum |
|--------|--------|---------|
| File size | ~200 lines | 500 lines |
| Class size | ~100 lines | 300 lines |

**Important note:** These are guidelines, not hard rules. A well-organized 600-line file can be cleaner than a poorly organized 200-line file. Focus on clarity and organization first.

**When files grow too large:**
- Extract classes or modules
- Split into multiple related files
- Identify natural boundaries for separation

### 4.2 The Newspaper Metaphor

**Principle:** A source file should read like a newspaper article.

**Structure:**
- **Name (headline)**: Simple but explanatory - tells you if you're in the right place
- **Top (summary)**: High-level concepts and algorithms - the "what" and "why"
- **Middle**: Supporting details
- **Bottom (details)**: Low-level functions and implementation details

**Why this works:**
- Readers can decide quickly if they need to read further
- Important concepts come first
- Natural top-to-bottom flow
- Details are available but don't obscure the big picture

### 4.3 Vertical Openness Between Concepts

**Principle:** Separate different concepts with blank lines to create visual groups.

**Rules:**
- Blank line between methods
- Blank line between logical sections within a method
- Blank line after variable declarations before logic
- No blank lines within tightly coupled code
- Blank line between imports/usings and code

**Why:**
- Creates visual paragraphs
- Shows conceptual boundaries
- Makes scanning code easier
- Improves readability

### 4.4 Vertical Density

**Principle:** Related code should appear close together with minimal blank lines between.

**Rules:**
- Related variable declarations should be grouped
- No unnecessary blank lines within a concept
- No comments separating related fields
- Keep tightly coupled code together

**Why:**
- Shows that code is related
- Reduces scrolling to see related items
- Makes relationships obvious

### 4.5 Vertical Distance

**Principle:** Concepts that are closely related should be close vertically.

**Key rules:**

**Variable declarations:**
- Declare close to usage
- Loop variables in loop header
- Local variables at top of short functions
- For long functions, declare near first usage

**Instance variables:**
- Declare at the top of the class
- All together in one place
- Well-known location (convention)

**Dependent functions:**
- Caller should be above callee
- Creates top-down flow
- Natural reading order
- Called functions close to caller

**Conceptual affinity:**
- Related functions grouped together
- Similar operations near each other
- Common pattern: group by feature not by type

### 4.6 Vertical Ordering

**Principle:** The most important concepts should come first, with increasing detail as you read down.

**Standard order for classes:**
1. Public constants
2. Private static variables
3. Private instance variables
4. Constructors
5. Public methods (high-level operations)
6. Private methods (low-level details, in order called)

**Why this order:**
- Important/public information first
- Top-down narrative flow
- Matches how we think about the class
- Details hidden until needed

## Horizontal Formatting

### 4.7 Line Length

**Guidelines:**

| Metric | Target | Maximum |
|--------|--------|---------|
| Line width | 80-100 chars | 120 chars |

**Why shorter lines:**
- No horizontal scrolling needed
- Easier to read (eye tracking research)
- Allows side-by-side comparison
- Works well on smaller screens
- Multiple files can be viewed side-by-side

**When lines are too long:**
- Extract variables or methods
- Break into multiple lines
- Simplify complex expressions

### 4.8 Horizontal Openness and Density

**Principle:** Use whitespace to show relationships and precedence.

**Rules:**
- Spaces around binary operators (=, +, -, *, /)
- No space between function name and opening parenthesis
- Space after commas in parameter lists
- Spaces can show operator precedence
- No space inside parentheses (usually)

**Why:**
- Shows related vs separate concepts
- Makes precedence clear
- Improves readability
- Reduces errors from misreading

### 4.9 Horizontal Alignment

**Principle:** Don't align variable declarations or assignments - it draws attention to the wrong thing.

**Why alignment is bad:**
- Draws eye to variable names, not types or meanings
- Creates busywork maintaining alignment
- Hides the true structure
- Breaks when viewed in different fonts
- If list needs alignment, it's probably too long

**Fix:** Keep declarations unaligned and simple. If the list is too long to read, the class is too large.

### 4.10 Indentation

**Principle:** Indentation shows the hierarchical structure of code. It makes scope visually obvious.

**Rules:**
- Consistent indentation (spaces or tabs, not mixed)
- Each scope level indented one level deeper
- Never break indentation for short statements
- File level: 0 indent
- Class level: 1 indent
- Method level: 2 indents
- Block level: 3+ indents

**Why:**
- Makes scope obvious at a glance
- Shows nesting depth
- Makes structure clear
- Consistency reduces cognitive load

**Common mistake:** Collapsing short statements onto one line breaks visual structure and should be avoided.

## Team Rules

### 4.11 Consistent Standards

**The most important rule:** Every team should agree on formatting standards and follow them consistently.

**Why team standards matter:**
- Eliminates formatting debates
- Reduces diff noise in version control
- Makes all code look familiar
- Removes personal preferences from code reviews
- Reduces cognitive load when reading others' code

**How to enforce:**
- Document standards in team wiki
- Use automated formatters
- Add formatter checks to CI/CD
- Run formatter on pre-commit hooks
- Make it automatic, not manual

**Tools for consistency:**
- **EditorConfig** - Cross-IDE formatting rules
- **Roslyn Analyzers / StyleCop** - C# code style enforcement
- **Prettier** - JavaScript/TypeScript
- **Black / Ruff** - Python
- **RustFmt** - Rust
- **gofmt** - Go

### 4.12 Automation Over Documentation

**Principle:** Automate formatting enforcement rather than relying on documentation and discipline.

**Benefits:**
- No debates during code review
- Consistent formatting everywhere
- New team members follow standards automatically
- Reduces mental overhead
- Eliminates formatting as a code review concern

**Implementation:**
- Configure IDE/editor to format on save
- Add pre-commit hooks to format code
- Run formatter checks in CI/CD
- Fail builds on formatting violations
- Make it impossible to commit badly formatted code

## Quick Checklist

### Vertical
- [ ] Files are small (< 500 lines, ideally < 200)
- [ ] Related code is close together
- [ ] Concepts separated by blank lines
- [ ] Caller above callee (top-to-bottom reading)
- [ ] Variables declared close to usage
- [ ] Instance variables at class top
- [ ] High-level concepts before low-level details

### Horizontal
- [ ] Lines are short (< 120 characters, ideally < 100)
- [ ] Proper indentation shows hierarchy
- [ ] Whitespace used to show relationships
- [ ] No horizontal alignment of declarations
- [ ] Consistent spacing around operators
- [ ] No broken indentation for short statements

### Team
- [ ] Team coding standards documented
- [ ] EditorConfig or formatting tool configured
- [ ] Automated formatting in pre-commit hooks
- [ ] Formatting checks in CI/CD
- [ ] No mixture of tabs and spaces
- [ ] All team members use same formatter settings

## When to Apply

**Always:**
- Follow team standards consistently
- Run automated formatter before committing
- Keep files and lines short
- Use vertical and horizontal whitespace intentionally

**During code review:**
- Check that formatting is consistent
- Ensure automated formatter was run
- Don't debate formatting if standards exist
- Focus on logic, not style

**When joining a project:**
- Learn the existing formatting standards
- Configure your IDE to match
- Follow the established patterns
- Don't reformat existing code unless specifically working on it

## See Examples

For detailed code examples demonstrating these principles:
- [C# Formatting Examples](../examples/csharp/04-formatting.md) - Before/after formatting examples

## See Also

- [Classes](09-classes.md) - Class organization and structure
- [Functions](02-functions.md) - Function organization
- [Meaningful Names](01-meaningful-names.md) - How names affect readability
