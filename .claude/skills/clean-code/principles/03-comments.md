# Principle 3: Comments

> "Don't comment bad code - rewrite it."
> — Brian W. Kernighan and P. J. Plaugher

## Overview

Comments are, at best, a necessary evil. The proper use of comments is to compensate for our failure to express ourselves in code. Every comment represents a failure to make the code self-explanatory. Good code mostly documents itself through clear names and structure.

## Why It Matters

- **Code changes, comments don't**: Comments quickly become outdated lies
- **Comments indicate failure**: Need for comments suggests code isn't clear enough
- **Maintenance burden**: Comments need maintenance just like code
- **Source of confusion**: Incorrect comments are worse than no comments
- **Distraction**: Too many comments make code harder to read

## The Truth About Comments

**Key insights:**
- Comments do not make up for bad code
- Code changes frequently, comments often don't, making them lies
- The only truly good comment is the comment you found a way not to write
- Comments are a last resort, not a first option
- Good names and clean code structure eliminate most need for comments

**The ideal:**
- Spend time making code expressive rather than writing comments
- Use comments only when code alone cannot express intent
- Keep comments minimal and focused on "why" not "what"

## Good Comments

### 3.1 Legal Comments

**Purpose:** Copyright and license statements required at the beginning of files.

**When to use:**
- Copyright notices
- License information
- Legal disclaimers
- Regulatory compliance

**Guidelines:**
- Keep them short and reference external documents
- Place at the top of the file
- Use standard headers from legal team

### 3.2 Informative Comments

**Purpose:** Provide information that cannot be expressed in code alone.

**Valid uses:**
- Explain complex regular expressions
- Document format specifications
- Describe abstract return values
- Clarify intent of abstract methods

**Guidelines:**
- Better: Move this information into function names or constants
- Use only when code truly cannot express the information
- Be concise and accurate

### 3.3 Explanation of Intent

**Purpose:** Explain WHY a decision was made, not WHAT the code does.

**Valid uses:**
- Explain trade-offs in design decisions
- Document why a particular algorithm was chosen
- Clarify business rules that seem odd
- Explain performance considerations
- Document non-obvious optimizations

**Guidelines:**
- Focus on intent and reasoning
- Include data or measurements when relevant
- Link to tickets or documents for more context
- Never explain what the code does (that should be obvious)

### 3.4 Warning of Consequences

**Purpose:** Warn other developers about consequences they might not anticipate.

**Valid uses:**
- Performance warnings (long-running operations)
- Memory usage warnings
- Thread safety issues
- External dependencies or side effects
- Breaking changes or deprecated patterns

**Guidelines:**
- Be specific about the consequence
- Suggest alternatives if available
- Use clear, attention-grabbing language
- Consider using attributes or annotations instead

### 3.5 TODO Comments

**Purpose:** Mark work that needs to be done in the future.

**Guidelines:**
- Use sparingly and track them
- Include ticket references (JIRA-1234)
- Add context about why it's a TODO
- Set a date or milestone when possible
- Review and clean up regularly
- Don't use as excuse for bad code
- Better: Create actual tickets and remove TODOs

**Format:**
- TODO: [Ticket-ID] Brief description
- TODO: Remove after X event/date
- TODO: Waiting on Y dependency

### 3.6 Amplification

**Purpose:** Emphasize importance of something that might seem inconsequential but is critical.

**Valid uses:**
- Highlight subtle bugs that can occur
- Emphasize critical security considerations
- Call attention to order dependencies
- Warn about non-obvious side effects

**Guidelines:**
- Use rarely and only for truly critical items
- Explain the consequence of ignoring it
- Consider if design can make this unnecessary

### 3.7 Public API Documentation

**Purpose:** Document public APIs in libraries and frameworks for external consumers.

**Valid uses:**
- XML documentation for library methods
- Parameter descriptions and constraints
- Return value descriptions
- Exception documentation
- Usage examples for complex APIs

**Guidelines:**
- Required for public APIs
- Keep synchronized with code
- Focus on contract and usage
- Include examples for complex scenarios
- Document preconditions and postconditions
- List all exceptions that can be thrown

## Bad Comments

### 3.8 Mumbling / Redundant Comments

**Problem:** Comments that simply restate what the code already clearly says.

**Why it's bad:**
- Adds no information
- Creates maintenance burden
- Clutters code
- Insults the reader's intelligence

**Examples of redundancy:**
- Getter/setter descriptions that just repeat the property name
- "Default constructor" for parameterless constructors
- "The X" for a field named X
- Method comments that restate the method name

**Fix:** Delete the comment and improve names if needed.

### 3.9 Misleading Comments

**Problem:** Comments that don't accurately describe what the code actually does.

**Why it's bad:**
- Worse than no comment at all
- Leads to bugs and confusion
- Results from code changing after comment was written
- Wastes developer time

**How to avoid:**
- Delete comments rather than updating them
- Make code self-explanatory instead
- Use tests to document behavior
- Keep comments minimal to reduce drift

### 3.10 Mandated Comments

**Problem:** Comments required by policy but add no value (e.g., every function must have a Javadoc).

**Why it's bad:**
- Creates noise and clutter
- Developers write meaningless comments to satisfy rules
- Makes finding useful comments harder
- Creates false sense of documentation

**Fix:**
- Require comments only for public APIs
- Focus on code quality, not comment quantity
- Use linting tools to enforce meaningful standards

### 3.11 Journal Comments

**Problem:** Change logs maintained in source files listing every modification.

**Why it's bad:**
- Version control systems do this better
- Clutters the top of files
- Never gets updated properly
- Irrelevant to current code

**Fix:**
- Delete journal comments entirely
- Use git log to see history
- Use git blame to see who changed what
- Link commits to tickets for context

### 3.12 Noise Comments

**Problem:** Comments that provide no information beyond what's already obvious.

**Why it's bad:**
- Pure clutter
- Train developers to ignore all comments
- Waste time reading them
- Suggest the coder was just filling space

**Common patterns:**
- "The X" for a variable named X
- "Constructor" for constructors
- "Set the X" for setters
- Restating each line of obvious code

**Fix:** Delete all noise comments immediately.

### 3.13 Position Markers

**Problem:** Banner comments or section separators in code.

**Why it's bad:**
- If you need markers, your file is too long
- Adds visual noise
- Becomes meaningless if overused
- Better solved by file organization

**Examples:**
- "======== PRIVATE METHODS ========"
- "// Actions /////////////////"
- Banner boxes around sections

**Fix:**
- Split large files into smaller ones
- Use proper file organization
- Let IDEs show structure
- Use whitespace for separation

### 3.14 Closing Brace Comments

**Problem:** Comments after closing braces indicating what they close.

**Why it's bad:**
- If you need them, your functions are too long
- Modern IDEs show matching braces
- Adds clutter
- Gets out of sync with code

**Fix:**
- Keep functions small (no need for markers)
- Use IDE features for brace matching
- Extract nested blocks to functions

### 3.15 Commented-Out Code

**Problem:** Code that's been commented out rather than deleted.

**Why it's bad:**
- Others won't delete it (might be important?)
- Accumulates over time
- Makes code harder to read
- Version control keeps history

**Fix:**
- DELETE IT IMMEDIATELY
- Trust version control
- If unsure, git has the history
- If it's not important enough to keep, delete it
- If it is important, keep it working

## Replace Comments with Better Code

**Strategy:** Instead of writing comments, improve the code to be self-explanatory.

### Technique 1: Extract Method with Descriptive Name

**Instead of:** Comment explaining complex logic
**Do this:** Extract logic to a well-named method

**Benefit:** Code becomes self-documenting, logic becomes reusable

### Technique 2: Extract Variable with Meaningful Name

**Instead of:** Comment explaining complex condition
**Do this:** Extract condition to a well-named boolean variable

**Benefit:** Condition becomes readable and testable

### Technique 3: Use Intention-Revealing Names

**Instead of:** Comment explaining what variable holds
**Do this:** Rename variable to reveal its purpose

**Benefit:** Eliminates need for comment entirely

### Technique 4: Simplify Complex Expressions

**Instead of:** Comment explaining complex boolean expression
**Do this:** Break into multiple named conditions or extract method

**Benefit:** Each piece of logic is independently understandable

## Quick Checklist

- [ ] No redundant comments explaining obvious code
- [ ] No commented-out code (delete it immediately!)
- [ ] No journal/change log comments
- [ ] No position markers, banners, or separators
- [ ] No closing brace comments
- [ ] Comments explain WHY, not WHAT
- [ ] Complex logic extracted to named methods instead of comments
- [ ] No misleading or outdated comments
- [ ] Public API documentation is complete and accurate
- [ ] TODOs include ticket references and are tracked

## When to Apply

**Ask these questions before writing a comment:**

1. Can I rename something to make this comment unnecessary?
2. Can I extract a method to make this clear?
3. Can I simplify the code to make it obvious?
4. Is this comment explaining "what" instead of "why"?
5. Will this comment stay accurate as code changes?

**Only write a comment if:**
- You're documenting a public API
- You're explaining a non-obvious business rule
- You're warning about consequences
- You're explaining "why" not "what"
- You've exhausted code-based alternatives

## The Golden Rule

**"If you must write a comment, ask yourself: Can I rename something or extract a method to make this comment unnecessary?"**

Usually, the answer is yes.

## See Examples

For detailed code examples demonstrating these principles:
- [C# Comments Examples](../examples/csharp/03-comments.md) - Before/after refactoring examples

## See Also

- [Meaningful Names](01-meaningful-names.md) - Making code self-documenting through names
- [Functions](02-functions.md) - Extracting well-named functions instead of comments
- [Formatting](04-formatting.md) - Using whitespace instead of position markers
