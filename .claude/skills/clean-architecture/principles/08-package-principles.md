# Package Principles

> **Navigation:** [← Back to Principles](../PRINCIPLES.md)

## Cohesion Principles (What to Put in a Package)

### 1. Reuse/Release Equivalence Principle (REP)
- Only release things that make sense to reuse together
- Versioning applies to whole package

### 2. Common Closure Principle (CCP)
- Classes that change together should be packaged together
- Minimizes releases

### 3. Common Reuse Principle (CRP)
- Don't force users to depend on things they don't use
- Split packages if clients only need part of it

## Coupling Principles (Package Dependencies)

### 1. Acyclic Dependencies Principle (ADP)
- No circular dependencies between packages
- Dependency graph must be a DAG (Directed Acyclic Graph)

### 2. Stable Dependencies Principle (SDP)
- Depend in the direction of stability
- Unstable packages depend on stable packages
- Not the other way around

### 3. Stable Abstractions Principle (SAP)
- Stable packages should be abstract
- Unstable packages should be concrete
- Matches dependency flow

## Related Principles

- [The Dependency Rule](01-dependency-rule.md) - Package Principles apply this rule to modules
- [SOLID Principles](02-solid-principles.md) - Package Principles extend SOLID to the module level
- [Screaming Architecture](07-screaming-architecture.md) - Proper package organization reflects architectural intent
