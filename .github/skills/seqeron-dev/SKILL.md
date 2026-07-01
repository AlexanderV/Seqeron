---
name: seqeron-dev
description: Use the Seqeron.Genomics C# API directly (not via MCP) with the library's real conventions and guardrails. Triggers — "call Seqeron from C#", "write code using Seqeron.Genomics", "use the Seqeron library / API in code", "why does this throw LimitationPolicy / SeqeronLimitationException", "how do I set the limitation mode", "TryCreate vs the DnaSequence constructor", "bootstrap the policy in a test", "which C# method backs this MCP tool / Method ID", "guarded unit / MinimumMode in code". Developer-profile counterpart to bio-rigor: the CODE mechanics of construction, the 3-tier LimitationPolicy, the guarded-unit exception, and mapping MCP Method IDs to real calls.
allowed-tools: Read, Grep, Glob
---

# seqeron-dev — using the Seqeron.Genomics C# API

For the **developer profile**: calling `Seqeron.Genomics` from C# code, not through MCP. This skill
is the code-side counterpart of `bio-rigor`. It does **not** restate rigor rules and it does **not**
cover *developing* the library (tests/architecture/refactoring live in `clean-architecture`,
`clean-code`, and the testing campaigns) — it covers the **API mechanics** a caller needs.

Delegate scientific rigor (provenance, cross-checking, the STOP rule, the alpha/clinical caveat) to
[`bio-rigor`](../bio-rigor/SKILL.md). This skill tells you *how to write the code*.

## Namespaces (verified in source)

There is **no single `Seqeron.Genomics` namespace**; types live in sub-namespaces. The README
Quickstart writes `using Seqeron.Genomics;` as shorthand, but the real declarations are:

| Type | Real namespace | Source |
|---|---|---|
| `DnaSequence`, `RnaSequence`, `ProteinSequence` | `Seqeron.Genomics.Core` | `src/.../Seqeron.Genomics.Core/DnaSequence.cs:6` |
| `LimitationPolicy`, `LimitationMode`, `SeqeronLimitationException`, `LimitationCatalog` | `Seqeron.Genomics.Core` | `src/.../Seqeron.Genomics.Core/LimitationPolicy.cs:6` |
| `SequenceAligner` | `Seqeron.Genomics.Alignment` | `src/.../Seqeron.Genomics.Alignment/SequenceAligner.cs:9` |
| `AlignmentResult`, `AlignmentStatistics`, `MultipleAlignmentResult` | `Seqeron.Genomics.Infrastructure` | `src/.../Seqeron.Genomics.Infrastructure/AlignmentTypes.cs:3` |
| `PrimerDesigner` | `Seqeron.Genomics.MolTools` | `src/.../Seqeron.Genomics.MolTools/PrimerDesigner.cs:7` |

(Paths relative to repo root `../../../` from here.) When unsure which namespace/type an operation
lives in, use [`seqeron-discovery`](../seqeron-discovery/SKILL.md).

## Constructing sequences: throwing ctor vs `TryCreate`

Two paths, both real (`DnaSequence.cs`):

```csharp
using Seqeron.Genomics.Core;

// Throwing constructor — throws on invalid input. Use when input is trusted.
var dna = new DnaSequence("AAAGAATTCAAA");          // DnaSequence.cs:22
Console.WriteLine(dna.GcContent());                 // DnaSequence.cs:82

// Validation-friendly path — no exception; returns false on invalid input.
if (!DnaSequence.TryCreate("ACGTNN", out var seq))  // DnaSequence.cs:129
    Console.WriteLine("Invalid DNA sequence");
else
    Console.WriteLine(seq.GcContent());
```

Prefer `TryCreate` for any untrusted / user-supplied input (the "validation-friendly path" in
[`README.md`](../../../README.md)). Use the constructor only when a failure *should* throw.

## The 3-tier LimitationPolicy (read/set in code)

`LimitationPolicy` is a **static class** with an `AsyncLocal` scoped override
(`LimitationPolicy.cs:23`). Modes are the enum `LimitationMode` — ordered
`Strict < Moderate < Permissive`, **default `Moderate`** (`LimitationPolicy.cs:33`):

```csharp
using Seqeron.Genomics.Core;

LimitationPolicy.DefaultMode = LimitationMode.Moderate;   // process-wide (LimitationPolicy.cs:33)
var effective = LimitationPolicy.CurrentMode;             // scoped-or-default (LimitationPolicy.cs:37)

using (LimitationPolicy.Use(LimitationMode.Permissive))   // scoped region (LimitationPolicy.cs:54)
{
    // guarded branch reachable only inside this scope
}
// convenience: LimitationPolicy.UsePermissive() / .UseStrict()  (LimitationPolicy.cs:57,60)
```

> Note: it is `LimitationMode.Permissive` (the **enum**), not `LimitationPolicy.Mode.Permissive`.

## The guarded-unit exception

When a task hits a guarded branch and the effective mode is **more restrictive** than that unit's
`MinimumMode`, the call throws `SeqeronLimitationException` (`LimitationPolicy.cs:333`, a
`sealed class : InvalidOperationException`). It carries `LimitationId`, `Category`, `MinimumMode`,
`Branch`, `RelatedTo`, `Workaround`, `ReportPath`, and a composed message naming the alternative.

**Do not** widen the mode just to force a number — that is the `bio-rigor` STOP rule. In production,
catch it and report the limitation + the documented workaround:

```csharp
try { /* guarded call */ }
catch (SeqeronLimitationException ex)
{
    // ex.LimitationId, ex.MinimumMode, ex.Workaround, ex.ReportPath — report, don't silently retry.
}
```

The 9 guarded units + their `MinimumMode` and the STOP rule: see
[`reference/limitation-policy.md`](reference/limitation-policy.md). Source of truth for scope:
[`docs/Validation/LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md).

## Test bootstrap (Permissive)

Tests that intentionally exercise a guarded best-effort branch must bootstrap `Permissive`, or the
default `Moderate`/`Strict` throws and the test fails for the wrong reason. The real pattern is a
`[ModuleInitializer]` that sets the assembly default
(`tests/.../Seqeron.Genomics.Tests/_LimitationPolicyTestBootstrap.cs:14-15`); production code stays
at `Moderate`. Exact snippet + the scoped-override alternative:
[`reference/limitation-policy.md`](reference/limitation-policy.md).

## MCP Method IDs → C# calls

Every per-tool MCP doc has a **`Method ID`** row (`Type.Method`) and a **source link** to the exact
`.cs#Lnnn`. That is the bridge from an MCP tool to the C# call:

- `docs/mcp/tools/moltools/primer_melting_temperature.md` → Method ID
  `PrimerDesigner.CalculateMeltingTemperature`, source `PrimerDesigner.cs#L197`.

To find the C# entry point for any tool, read its tool doc's `Method ID` + source link, or ask
[`seqeron-discovery`](../seqeron-discovery/SKILL.md). Algorithm contracts/invariants:
[`docs/algorithms/`](../../../docs/algorithms/). Construction, the Method-ID lookup recipe, and 2–3
worked pipelines: [`reference/api-conventions.md`](reference/api-conventions.md).
