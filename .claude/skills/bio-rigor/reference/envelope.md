# Operating envelope — how to check it before you compute

Seqeron algorithms are validated for a **stated scope**, not for everything. Many are faithful but
**simplified or subset** realisations of fuller published methods. Before reporting a result, confirm
the task is inside the unit's validated envelope. Do **not** force output outside it.

**Single source of truth:** [`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md)
(honest scope boundaries) and, per unit, `docs/Validation/reports/{UNIT}.md`. Do not restate the table
here — read it there so it can never drift.

## The runtime guard: `LimitationPolicy`

`Seqeron.Genomics.Core.LimitationPolicy` has **three modes**, least → most permissive:

- **`Strict`** — only the ideal *and complete* result; throws on every guarded branch.
- **`Moderate`** (**default**) — throws on **non-ideal-output** branches; allows
  **correct-but-incomplete / narrower-contract** branches.
- **`Permissive`** — allows everything (historical best-effort).

A guarded call throws `SeqeronLimitationException` when the effective mode is **more restrictive**
than that limitation's **minimum access mode**. The exception names the limitation, what it relates
to, and how to obtain the result another way.

### The STOP rule

If the task hits a guarded branch below its MinimumMode:

1. **STOP.** Do not silently widen the mode to `Permissive` just to get a number.
2. **Report** the limitation: name the unit, why it is out of envelope, and the named alternative
   (from the exception message / `LIMITATIONS.md`).
3. Only proceed at a wider mode if the **user explicitly accepts** the documented caveat — and then
   record the chosen mode in the provenance `Envelope:` line.

## The 9 guarded units (MinimumMode) — pointer, not a copy

`LIMITATIONS.md` lists them against `LimitationCatalog`. As of last review they group as:

- **MinimumMode `Permissive`** (non-ideal output, blocked in Strict & Moderate):
  `PARSE-FASTQ-001`, `CHROM-CENT-001`, `DISORDER-REGION-001`, `MIRNA-TARGET-001`, `MIRNA-CLEAVAGE-001`.
- **MinimumMode `Moderate`** (correct-but-incomplete, blocked only in Strict):
  `ONCO-MHC-001`, `ONCO-IMMUNE-001`, `META-BIN-001`, `PROBE-DESIGN-001`.

Two rows are **documented-only, no runtime guard** (`RNA-STRUCT-001`, `MIRNA-PRECURSOR-001`) — the
shortfall is undetectable per call, so envelope-awareness here means reporting the caveat, not
catching an exception. **Always read `LIMITATIONS.md` for the authoritative, current list.**

## Setting / scoping the mode (C# API)

```csharp
LimitationPolicy.DefaultMode = LimitationPolicy.Mode.Moderate;   // process-wide default
using (LimitationPolicy.Use(LimitationPolicy.Mode.Permissive))   // scoped region
{
    // guarded call that the caller has explicitly accepted the caveat for
}
```

### Test bootstrap (Permissive)

Tests that intentionally exercise a guarded branch must **bootstrap `Permissive`**, otherwise the
default `Moderate` throws and the test fails for the wrong reason. Wrap the guarded call in
`using (LimitationPolicy.Use(LimitationPolicy.Mode.Permissive)) { … }` (or set the default in the
fixture) so the branch is reachable. This is test-only; production code stays at `Moderate` and
respects the STOP rule above.

## MCP-client mode

MCP tools that wrap a guarded unit surface the same limitation as a structured error. Treat that
error exactly like the STOP rule: report the limitation and the alternative; do not retry with a
looser setting unless the user has accepted the documented caveat.
