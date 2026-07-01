# LimitationPolicy — the real C# API (developer reference)

The runtime envelope guard. This file documents the **code mechanics**; for *why* an envelope exists
and the rigor STOP rule, see [`bio-rigor`](../../bio-rigor/SKILL.md) and its
[`reference/envelope.md`](../../bio-rigor/reference/envelope.md). Source of truth for scope:
[`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md).

All symbols below are verified against
[`src/Seqeron/Algorithms/Seqeron.Genomics.Core/LimitationPolicy.cs`](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/LimitationPolicy.cs)
(one file: policy + enum + catalog + exception).

## Types & namespace

Namespace: **`Seqeron.Genomics.Core`** (`LimitationPolicy.cs:6`).

| Symbol | Kind | Source |
|---|---|---|
| `LimitationPolicy` | `static class` | `LimitationPolicy.cs:23` |
| `LimitationMode` | `enum { Strict, Moderate, Permissive }` | `LimitationPolicy.cs:106` |
| `SeqeronLimitationException` | `sealed class : InvalidOperationException` | `LimitationPolicy.cs:333` |
| `LimitationCatalog` | `static class` (id → `LimitationInfo`) | `LimitationPolicy.cs:160` |
| `LimitationInfo` | `sealed record` (`Id`, `MinimumMode`, `Branch`, `Workaround`, …) | `LimitationPolicy.cs:144` |

## Reading / setting the mode

```csharp
using Seqeron.Genomics.Core;

// Modes are ordered Strict < Moderate < Permissive; default is Moderate.
LimitationPolicy.DefaultMode = LimitationMode.Moderate;   // process-wide (cs:33)

LimitationMode m = LimitationPolicy.CurrentMode;          // scoped override ?? default (cs:37)
bool strict    = LimitationPolicy.IsStrict;               // cs:40
bool allowed   = LimitationPolicy.IsAllowed("META-BIN-001"); // cs:47 — vs that unit's MinimumMode

// Scoped, nestable, AsyncLocal — restores previous mode on Dispose (cs:54,79-97):
using (LimitationPolicy.Use(LimitationMode.Permissive)) { /* ... */ }
using (LimitationPolicy.UsePermissive())               { /* ... */ }   // cs:57
using (LimitationPolicy.UseStrict())                   { /* ... */ }   // cs:60
```

- `DefaultMode` is a process-wide singleton; `Use(...)` pushes an **async-local** override for the
  current flow only. `CurrentMode` returns the scoped value if set, else `DefaultMode`.
- **It is the enum value `LimitationMode.Permissive`**, not `LimitationPolicy.Mode.Permissive`
  (the enum is a sibling type, not nested). The `bio-rigor/reference/envelope.md` snippet writes
  `LimitationPolicy.Mode.Permissive` — that form does **not** compile; use `LimitationMode`.

## Catching the guarded-unit exception

A guarded branch calls `LimitationPolicy.Enforce(limitationId)` (`cs:72`); it throws when
`CurrentMode < info.MinimumMode` (`cs:75-76`). The exception exposes actionable fields (`cs:335-354`):

```csharp
try
{
    // e.g. a guarded MHC / disorder / miRNA / bin-quality call
}
catch (SeqeronLimitationException ex)
{
    // ex.LimitationId  — e.g. "ONCO-MHC-001"
    // ex.MinimumMode   — lowest LimitationMode that would allow it
    // ex.Branch        — the specific guarded branch
    // ex.Workaround    — how to get the ideal result (supply matrix / use ref tool / …)
    // ex.ReportPath    — docs/Validation/reports/<UNIT>.md
    // ex.Message       — composed, human-readable, names the alternative
}
```

**STOP rule (do not widen to force output):** report the limitation and its `Workaround`; only run at
a wider mode if the caller has explicitly accepted the documented caveat. See `bio-rigor`.

`LimitationCatalog.Get(id)` throws `KeyNotFoundException` for a non-guarded id (`cs:169-175`);
`LimitationCatalog.Entries` (`cs:165`) enumerates all guarded entries.

## The 9 guarded units + MinimumMode (pointer, not a copy)

The authoritative, current list is [`LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md)
(mirrored by `LimitationCatalog`, `LimitationPolicy.cs:177-319`). As of last review they group as:

- **MinimumMode `Permissive`** — non-ideal output, blocked in Strict **and** Moderate:
  `PARSE-FASTQ-001`, `CHROM-CENT-001`, `DISORDER-REGION-001`, `MIRNA-TARGET-001`, `MIRNA-CLEAVAGE-001`.
- **MinimumMode `Moderate`** — correct-but-incomplete / narrower-contract, blocked only in Strict:
  `ONCO-MHC-001`, `ONCO-IMMUNE-001`, `META-BIN-001`, `PROBE-DESIGN-001`.

Each catalog entry carries an exhaustive `Workaround` string (e.g. *supply a scoring matrix to
`PredictIc50Smm`*, *pass lineage markers to `EstimateBinQualityFromMarkers`*, *use `PredictDisorderRegions`
for the ideal result*). Two rows (`RNA-STRUCT-001`, `MIRNA-PRECURSOR-001`) are **documented-only with
no runtime guard** — envelope-awareness there means reporting the caveat, not catching an exception.
**Always read `LIMITATIONS.md`** for the live list rather than trusting this grouping.

## Test bootstrap (Permissive)

Real pattern — a per-assembly `[ModuleInitializer]` that sets the default so canonical/discipline
fixtures can exercise the best-effort branches
([`_LimitationPolicyTestBootstrap.cs:12-16`](../../../../tests/Seqeron/Seqeron.Genomics.Tests/_LimitationPolicyTestBootstrap.cs)):

```csharp
using System.Runtime.CompilerServices;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

internal static class LimitationPolicyTestBootstrap
{
    [ModuleInitializer]
    internal static void Init() => LimitationPolicy.DefaultMode = LimitationMode.Permissive; // :15
}
```

The same one-liner appears in `Seqeron.Mcp.Parsers.Tests` and `Seqeron.Mcp.Analysis.Tests`
bootstraps. Tests that must verify the **production** Strict/Moderate behaviour re-tighten locally
with a scoped override — e.g.
[`LimitationPolicy_Strict_Tests.cs:32`](../../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/LimitationPolicy_Strict_Tests.cs):

```csharp
using (LimitationPolicy.UseStrict())   // exercise the default production guard within a Permissive assembly
{
    // assert SeqeronLimitationException is thrown
}
```

This is **test-only**: production code stays at `Moderate` and respects the STOP rule.
