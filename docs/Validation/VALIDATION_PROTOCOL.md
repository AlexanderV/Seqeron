# Biological Algorithm Validation Protocol

**Library:** Seqeron.Genomics
**Started:** 2026-06-12
**Scope (phase 1):** the 86 *implemented* test units (status ☑ in `ALGORITHMS_CHECKLIST_V2.md`).
**Tracker:** [VALIDATION_LEDGER.md](VALIDATION_LEDGER.md)
**Reports:** per-unit reports were consolidated into the ledger + findings register and archived in git (`git show cb113ce:docs/Validation/reports/{UnitID}.md`).

---

## Purpose

This is a critical library. Existing artifacts (TestSpecs, Evidence docs, unit tests) assert
correctness, but they were authored alongside the code and may share the same blind spots.
This campaign **independently re-validates** each algorithm against authoritative external
sources, in a fresh context per algorithm, in two ordered stages:

1. **Stage A — Validate the description** (the biology/maths is correct *in the abstract*).
2. **Stage B — Validate the implementation** (the code faithfully realises the validated description).

Stage B is only meaningful once Stage A passes. If the description is wrong, fixing the code
to match it would be wrong too. Always do A before B.

---

## One algorithm = one session (fresh context)

Each test unit is validated in its own session so that no reasoning from a previous unit
contaminates the next. Use the prompt template below to start each session. Record the
outcome in the ledger and write a report.

### Per-session prompt template

> We are validating the Seqeron.Genomics library, one biological algorithm per session.
> Validate test unit **`{UNIT-ID}`** following `docs/Validation/VALIDATION_PROTOCOL.md`.
>
> Do **Stage A (description)** first, then **Stage B (implementation)**. Use authoritative
> external sources (peer-reviewed papers, reference textbooks, Wikipedia with its cited
> primary refs, reference implementations such as Biopython/EMBOSS/samtools, and public
> datasets / Rosalind problems). Cross-check at least one independent reference implementation
> or hand-computed dataset.
>
> Produce a report at `docs/Validation/reports/{UNIT-ID}.md` using the report template in the
> protocol and give a final verdict (PASS / PASS-WITH-NOTES / FAIL) for each stage.
>
> **Completion criteria — a session ends in exactly one of two states:**
> 1. **FIXED / CLEAN** — either no defect was found, or every defect found was *completely
>    fixed* in this session: code corrected, tests added/updated to lock the sourced values,
>    and `dotnet build` + the unit's tests pass. The algorithm is fully functional.
> 2. **LIMITED** — a defect or gap exists that you could *not* completely fix in this session.
>    End by recording, in the report, precisely **why** it is limited and **what** is missing
>    for it to be fully functional (root cause, the correct behaviour per source, why the fix
>    is out of reach here — e.g. needs a dataset, a model, an upstream API, or a larger
>    redesign). Never leave a half-fix or "fix to match a wrong spec".
>
> Do not silently change code to match a description that Stage A proved wrong — fix the
> description too. Report the chosen end-state clearly so the orchestrator can set the status.

---

## Environment (build & test)

The repo targets **net10.0**. The SDK is installed at `~/.dotnet` (no system PATH). Prefix
commands or export PATH:

```bash
export PATH="$HOME/.dotnet:$PATH" DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
dotnet build tests/Seqeron/Seqeron.Genomics.Tests/Seqeron.Genomics.Tests.csproj -c Debug
dotnet test  tests/Seqeron/Seqeron.Genomics.Tests/Seqeron.Genomics.Tests.csproj -c Debug --no-build \
  --filter FullyQualifiedName~<TestClassName>
```

Green baseline (2026-06-12): `Seqeron.Genomics.Tests` = 4484 passed, 0 failed (after the 3 perfect
fixes: RNA IUPAC complement, Newick multifurcation rejection, IUPAC positive-set matching). Any
session that touches code must leave this project building and all its tests passing (or explain LIMITED).

## Stage A — Validate the description

Inputs: `tests/TestSpecs/{UnitID}.md`, `docs/Evidence/{UnitID}-Evidence.md` (if present),
the unit block in `ALGORITHMS_CHECKLIST_V2.md`.

Checklist:

- [ ] **Source quality** — sources cited are authoritative and actually say what is claimed
      (open them; do not trust the citation label). Prefer primary literature over Wikipedia.
- [ ] **Formula correctness** — every formula matches the cited source exactly (symbols,
      normalisation, units, log base, edge conventions).
- [ ] **Definitions & conventions** — coordinate base (0/1), inclusive/exclusive ends, strand,
      percentage vs fraction, ambiguity-code handling — all explicit and standard.
- [ ] **Edge-case semantics** — empty / single element / all-same / boundary / invalid input
      have a *defined, sourced* expected behaviour (not "implementation-defined").
- [ ] **Independent cross-check** — at least one of: a reference tool's output, a worked
      example from a paper, a Rosalind dataset, or a hand computation. Record exact numbers.
- [ ] **Invariants** — listed invariants are genuinely true mathematical properties.

Stage A verdict: PASS / PASS-WITH-NOTES (minor, documented divergences) / FAIL (description
is biologically or mathematically wrong → stop, log, do not proceed to B until resolved).

## Stage B — Validate the implementation

Inputs: the canonical method(s) and all variants/delegates listed in the unit; the test file(s).

Checklist:

- [ ] **Code realises the formula** — read the implementation; confirm it computes the
      validated formula, not an approximation, for the right inputs.
- [ ] **Edge cases in code** — each Stage-A edge case is handled as specified (trace or test).
- [ ] **Cross-verification values** — recompute the unit's cross-check table against the actual
      code (run the tests or trace by hand); every value matches the external reference.
- [ ] **Variant consistency** — delegates/`*Fast`/instance properties agree with the canonical.
- [ ] **Numerical robustness** — no precision loss, overflow, or div-by-zero on stated ranges.
- [ ] **Tests are real** — assertions check exact sourced values (not just "no throw" or
      tautologies), are deterministic, and actually cover the Stage-A edge cases.

Stage B verdict: PASS / PASS-WITH-NOTES / FAIL (code diverges from validated description → log
defect with minimal repro; fix only after approval).

---

## Report template (`docs/Validation/reports/{UnitID}.md`)

```markdown
# Validation Report: {UNIT-ID} — {Algorithm name}

- **Validated:** {date}   **Area:** {area}
- **Canonical method(s):** {…}
- **Stage A verdict:** PASS / PASS-WITH-NOTES / FAIL
- **Stage B verdict:** PASS / PASS-WITH-NOTES / FAIL

## Stage A — Description
- Sources opened & what they confirm
- Formula check (cite source line/equation)
- Edge-case semantics check
- Independent cross-check (numbers)
- Findings / divergences

## Stage B — Implementation
- Code path reviewed (file:line)
- Formula realised correctly? (evidence)
- Cross-verification table recomputed vs code
- Variant/delegate consistency
- Test quality audit
- Findings / defects

## Verdict & follow-ups
- Final verdict + any logged defects/issues
```

---

## Verdict legend (ledger)

Per-stage verdict (Stage A / Stage B columns):

| Symbol | Meaning |
|--------|---------|
| ✅ | PASS — independently confirmed correct |
| 🟡 | PASS-WITH-NOTES — correct, minor documented divergence |
| ❌ | FAIL — defect found (see report) |
| ⬜ | not yet validated |

End-state (State column) — every finished session is exactly one of:

| Symbol | Meaning |
|--------|---------|
| ✅ CLEAN | no defect, or every defect completely fixed; algorithm fully functional |
| 🔧 LIMITED | defect/gap could not be fully fixed; report explains why & what is missing |
| ⬜ | not yet processed |
