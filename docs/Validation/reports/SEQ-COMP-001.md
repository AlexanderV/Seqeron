# Validation Report: SEQ-COMP-001 — DNA/RNA Complement (Watson–Crick base complement)

- **Validated:** 2026-06-24   **Area:** Composition
- **Canonical method(s):** `SequenceExtensions.GetComplementBase(char)`, `SequenceExtensions.GetRnaComplementBase(char)`, `SequenceExtensions.TryGetComplement(ReadOnlySpan<char>, Span<char>)`; delegates `DnaSequence.Complement()`, `RnaSequence.Complement()`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN

> Note: the orchestration prompt labelled this unit "nucleotide composition / base counts". That is incorrect — "Composition" is only the *area* name. SEQ-COMP-001 is the **DNA/RNA single-base complement** unit (Watson–Crick + IUPAC complement), as defined in `tests/TestSpecs/SEQ-COMP-001.md` and `ALGORITHMS_CHECKLIST_V2.md`. Validated the unit as actually defined.

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Nucleic acid notation** (https://en.wikipedia.org/wiki/Nucleic_acid_notation), fetched live. Full IUPAC/NC-IUB 1984 complement table retrieved:
  A→T, C→G, G→C, T→A, U→A, W→W, S→S, M→K, K→M, R→Y, Y→R, B→V, D→H, H→D, V→B, N→N, gap (-)→gap. This matches the TestSpec §2 Source 2 table and the code switch **exactly**.
- **Biopython Bio.Seq API** (https://biopython.org/docs/latest/api/Bio.Seq.html), fetched live. Docstring examples confirmed:
  - `complement("CGA")` → `"GCT"`
  - `complement("CGAUT")` → `"GCTAA"` (U treated as T in DNA context)
  - `complement("ACGTUacgtuXYZxyz")` → `"TGCAAtgcaaXRZxrz"` (case-preserving; X/Z pass through; Y→R)
  - `complement_rna("ACGTUacgtuXYZxyz")` → `"UGCAAugcaaXRZxrz"` (T treated as U; U output)

### Formula check
Complement is a per-base bijection defined by Watson–Crick pairing (A↔T/U, G↔C) extended to IUPAC ambiguity codes by complementing the underlying base *set* (NC-IUB 1984). Every entry in the code's `switch` matches the Wikipedia table cell-for-cell. The RNA path is identical except the base alphabet emits U instead of T (A→U, U→A) — matching Biopython `complement_rna`.

### Edge-case semantics (all defined & sourced)
- Empty sequence → success, no output (standard API contract).
- Lowercase input → uppercase output (deliberate divergence D1 from Biopython's case-preservation; justified by `DnaSequence`/`RnaSequence` uppercase-normalisation contract).
- U in DNA context → A; T in RNA context → A (Biopython-sourced).
- Non-IUPAC chars (gap `-`, `.`, `?`, `*`, X, Z) pass through unchanged (Biopython + IUPAC gap rule).
- Destination too small → `false` (Try-pattern contract).

### Independent cross-checks (exact numbers)
- DNA `complement("ACTG-NH")` hand-traced: A→T C→G T→A G→C -→- N→N H→D = **`TGAC-ND`** — matches IUPAC table (H complement = D) and the TestSpec cross-verification row.
- DNA full string with our uppercase convention: `"ACGTUacgtuXYZxyz"` → recognised bases uppercased, unknown x/z verbatim = **`TGCAATGCAAXRZxRz`** (Biopython `TGCAAtgcaaXRZxrz` modulo case). Test asserts this exactly.
- RNA `complement_rna("CGAUT")` → **`GCUAA`**; full IUPAC string `"ACGURYSWKMBDHVN"` → **`UGCAYRSWMKVHDBN`** — hand-verified per-char.

### Findings / divergences
- D1 (case → uppercase) is documented and intentional; not a defect.
- D2 (IUPAC ambiguity codes) previously a divergence, now RESOLVED (MUST-11) — code complements all 11 codes per NC-IUB 1984.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:251` (`GetComplementBase`), `:288` (`GetRnaComplementBase`), `:317` (`TryGetComplement`), `:333` (`TryGetReverseComplement`).
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs:54` (`Complement`) — delegates to `GetComplementBase`.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/RnaSequence.cs:54` (`Complement`) — inline A/U/C/G switch.
- Recently changed in commit `6e900e92` (opt-in Biopython/VCF compatibility modes); re-validated against live sources.

### Formula realised correctly?
Yes. The `switch` expressions are exhaustive over the IUPAC alphabet and emit precisely the table values; the `_ => nucleotide` default implements pass-through for non-IUPAC characters. `TryGetComplement` length-guards then maps each position with `GetComplementBase` (O(n), O(1) per base).

### Cross-verification table recomputed vs code
| Input | Function | Code output | Reference | Match |
|-------|----------|-------------|-----------|-------|
| `ACGT` | GetComplementBase | `TGCA` | Biopython | ✅ |
| `CGAUT` | GetComplementBase | `GCTAA` | Biopython | ✅ |
| `ACTG-NH` | GetComplementBase | `TGAC-ND` | IUPAC + Biopython | ✅ |
| `ACGTUacgtuXYZxyz` | GetComplementBase | `TGCAATGCAAXRZxRz` | Biopython (mod case D1) | ✅ |
| `CGAUT` | GetRnaComplementBase | `GCUAA` | Biopython complement_rna | ✅ |
| `ACGURYSWKMBDHVN` | GetRnaComplementBase | `UGCAYRSWMKVHDBN` | IUPAC RNA | ✅ |

All values reproduced by the live test run (445 complement tests, 0 failed).

### Variant/delegate consistency
- `DnaSequence.Complement()` / `ReverseComplement()` / `GetReverseComplementString` all route through `GetComplementBase` — consistent by construction.
- `RnaSequence.Complement()` uses a private A/U/C/G switch. `RnaSequence.ValidateSequence` rejects anything outside {A,C,G,U} at construction, so ambiguity codes/T can never reach it; over its *reachable* input domain the inline switch is total and agrees with `GetRnaComplementBase`. No inconsistency.

### Numerical robustness
No arithmetic; pure char mapping. Length guard prevents buffer overrun; empty input is a no-op. No precision/overflow concerns.

### Test quality audit
`SequenceExtensions_Complement_Tests.cs` asserts exact sourced values (Watson–Crick, all 11 IUPAC codes both cases, involution, T→A RNA context, gap/unknown pass-through, buffer-size edges incl. destination-untouched checks, and 5 explicit Biopython cross-verification cases). Assertions are deterministic and value-based, not "no-throw" tautologies. Covers every Stage-A edge case.

## Verdict & follow-ups
Both stages PASS. Implementation faithfully realises the IUPAC/Biopython-sourced complement description; all cross-checks reproduced exactly against live external sources. No code changes required. **State: CLEAN.** No follow-ups.
