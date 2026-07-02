# Hexagonal Architecture (Ports & Adapters)

> **Navigation:** [← Back to Principles](../PRINCIPLES.md) | [Code Examples →](../examples/csharp/PATTERNS-EXAMPLES.md)

## Core Concept

Isolate the application core from external concerns using ports (interfaces) and adapters (implementations).

## Components

**Application Core:**
- Contains business logic
- Defines ports (interfaces)
- No knowledge of external world

**Ports:**
- Primary/Driving Ports: The application's own API — defined by the app and CALLED BY driving adapters (controllers, CLI) to drive the application inward.
- Secondary/Driven Ports: Interfaces the app defines and CALLS OUT through — implemented by driven adapters (database, email, external APIs).

**Adapters:**
- Primary/Driving Adapters: Call the application (REST API, GraphQL, CLI)
- Secondary/Driven Adapters: Called by application (Database, Email, External APIs)

## Key Insight

The application defines BOTH ports. Adapters conform to the application's needs, not the other way around.

## Related Principles

- [The Dependency Rule](01-dependency-rule.md) - Ports and adapters implement this rule
- [Dependency Inversion Principle](02-solid-principles.md#d---dependency-inversion-principle-dip) - DIP is the foundation of ports and adapters
- [Onion Architecture](05-onion-architecture.md) - Another visualization of the same concept
- [Domain-Driven Design](03-domain-driven-design.md) - Repositories are an example of Driven Ports
