# Clean Architecture Skill

Comprehensive guidance for building software with Clean Architecture, DDD, SOLID principles, and hexagonal/onion patterns.

## File Structure

```
clean-architecture/
├── SKILL.md              # Main skill file (start here!)
├── PRINCIPLES.md         # Detailed principles reference
├── PATTERNS.md           # Design patterns catalog
├── CHECKLIST.md          # Architecture review checklist
├── README.md             # This file
└── examples/
    └── csharp/           # C# / .NET Examples (extensible to other languages)
        ├── SKILL-EXAMPLES.md       # Project structure & layers ⭐
        ├── PRINCIPLES-EXAMPLES.md  # SOLID & DDD code examples ⭐
        ├── PATTERNS-EXAMPLES.md    # Design pattern implementations ⭐
        ├── CHECKLIST-EXAMPLES.md   # Complete feature with tests ⭐
        └── ERROR-HANDLING-EXAMPLES.md  # Exception strategy & Result pattern ⭐
```

> **Note:** Currently examples are provided for C#/.NET only. The structure is designed to be extensible - contributions for other languages (TypeScript, Java, Python, Go) are welcome!

## Quick Start

| Goal | Start Here | Then |
|------|------------|------|
| 🆕 New Project | [SKILL.md](SKILL.md) | [Project Structure](examples/csharp/SKILL-EXAMPLES.md) |
| 🔧 Refactoring | [CHECKLIST.md](CHECKLIST.md) | [PRINCIPLES.md](PRINCIPLES.md) |
| 👀 Code Review | [CHECKLIST.md](CHECKLIST.md) | [Anti-patterns in SKILL.md](SKILL.md) |
| 📐 Domain Modeling | [PRINCIPLES.md](PRINCIPLES.md) | [DDD Examples](examples/csharp/PRINCIPLES-EXAMPLES.md) |
| 🧩 Choose Pattern | [PATTERNS.md](PATTERNS.md) | [Pattern Examples](examples/csharp/PATTERNS-EXAMPLES.md) |
| ⚠️ Error Handling | [Error Handling](examples/csharp/ERROR-HANDLING-EXAMPLES.md) | - |
| ✅ Full Feature | [CHECKLIST-EXAMPLES.md](examples/csharp/CHECKLIST-EXAMPLES.md) | - |

## Key Principles

- **Dependency Rule:** Dependencies point inward only (Presentation → Application → Domain)
- **Theory** (top-level .md files) - Language-agnostic concepts
- **Practice** (examples/) - Production-ready code

## Related Skills

- **[Clean Code](../clean-code/README.md)** - Code-level quality: naming, functions, comments, error handling

## Extending This Skill

To add examples for another language:

1. Create folder: `examples/<language>/` (e.g., `examples/typescript/`)
2. Mirror the C# structure with language-specific implementations
3. Update this README to include the new language

**Start here:** [SKILL.md](SKILL.md) | **Validate:** [CHECKLIST.md](CHECKLIST.md)

Happy architecting! 🏗️
