# Validation Report: SEQ-SUMMARY-001 — Sequence Summary

- **Validated:** 2026-06-16   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.SummarizeNucleotideSequence(string?)` → `SequenceSummary` record (aggregates Length, GcContent, Entropy, Complexity, MeltingTemperature, Composition)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Unit nature

This is an **AGGREGATOR**, not a duplicate. `SummarizeNucleotideSequence` performs no new
computation; it copies the results of four already-validated sibling canonical methods into one
record:

| Summary field | Source method | Sibling unit |
|---|---|---|
| `Length` / `GcContent` / `Composition[A,T,G,C,U,N]` | `CalculateNucleotideComposition` | SEQ-COMPOSITION-001 / SEQ-STATS-001 |
| `Entropy` | `CalculateShannonEntropy` | SEQ-ENTROPY-PROFILE-001 |
| `Complexity` | `CalculateLinguisticComplexity` | (linguistic-complexity sibling) |
| `MeltingTemperature` | `CalculateMeltingTemperature(seq, useWallaceRule: len<14)` | SEQ-TM-001 / SEQ-THERMO-001 |

The unit's own correctness obligation is therefore **field-wise consistency**: each field must equal
its canonical method's value on the same input. The underlying per-metric correctness is the scope of
the sibling units. Validation here focuses on (a) the aggregated field definitions being sound against
external sources, and (b) the aggregation faithfully delegating with no rounding/alphabet/branch drift.

## Stage A — Description

### Sources opened this session (with extracted numbers)

1. **Biopython `Bio.SeqUtils.MeltingTemp`** (WebFetch, 2026-06-16):
   `https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py`
   - Wallace docstring verbatim: **"Tm = 4 degC * (G + C) + 2 degC * (A+T)"**; worked example
     `Tm_Wallace('ACGTTGCAATGCCGTA')` → **48.0**.
   - GC general form verbatim: **"Tm = A + B(%GC) - C/N + salt correction - D(%mismatch)"**;
     Marmur–Doty 1962 valueset **A=69.3, B=0.41, C=650** → "Tm = 69.3 + 0.41(%GC) - 650/N".
2. **Wikipedia "Entropy (information theory)"** (WebFetch, 2026-06-16):
   `https://en.wikipedia.org/wiki/Entropy_(information_theory)`
   - **H(X) := − ∑ p(x) log p(x)**; "Base 2 gives the unit of bits"; "The maximal entropy of an event
     with n different outcomes is **log_b(n)**: attained by the uniform distribution." → 4 equal
     symbols ⇒ H = log₂4 = 2.0 bits.
3. **Wikipedia "Nucleic acid thermodynamics"** (WebFetch, 2026-06-16) — confirmed the article does
   *not* state the repo's `64.9 + 41·(GC−16.4)/N` variant or the `<14` threshold (those belong to
   the SEQ-TM-001 sibling). Recorded as a Stage-A note, not a defect of the aggregator.

### Formula check

- **GcContent** = (G+C)/counted-bases, case-insensitive, 0 on empty — matches Biopython `gc_fraction`
  semantics (Evidence Source 1). ✔
- **Entropy** = −Σ p·log₂p in bits — matches Shannon 1948 / Wikipedia. ✔
- **Wallace Tm** = 2(A+T)+4(G+C) — matches Biopython docstring exactly; Biopython's own worked example
  `ACGTTGCAATGCCGTA` → 48.0 reproduced by the repo formula (independent hand-calc this session). ✔
- **Marmur–Doty branch** = `64.9 + 41·(GC−16.4)/N` — a widely-used Tm variant but **not** the exact
  Biopython Marmur–Doty valueset (69.3/0.41/650 with %GC). This is the SEQ-TM-001 unit's documented
  convention; SEQ-SUMMARY-001 only delegates and asserts equality with `CalculateMeltingTemperature`.
  → **PASS-WITH-NOTES** (the variant choice is sourced/owned by SEQ-TM-001, out of scope here).
- **Complexity** = mean of per-word-size vocabulary-usage ratios (not the Trifonov **product**).
  Documented divergence owned by the linguistic-complexity sibling, faithfully reported here.

### Edge-case semantics

- empty / null → degenerate summary (all zeros), guarded; null normalised to `""`. ✔
- case-insensitive (each per-metric method uppercases internally). ✔
- RNA `U` and ambiguous `N` counted as distinct composition keys; `T=0` for RNA. ✔
- Tm threshold strict `<14` (length 14 → GC branch). ✔

### Independent cross-check (hand-computed this session, Python)

| Input | Field | Expected (sourced) |
|---|---|---|
| ATGCATGC | GcContent | 0.5 (4/8) |
| ATGCATGC | Entropy | 2.0 (log₂4) |
| ATGCATGC | Wallace Tm | 24.0 (2·4+4·4) |
| ATGCATGC | Complexity | 0.8396825396825397 (mean of U₁..U₆ = 1, 4/7, 2/3, 4/5, 1, 1) |
| ATGCATGCATGCATGC | Marmur–Doty Tm | 43.375 (64.9+41·(8−16.4)/16) |
| ATGC | Wallace Tm | 12.0 (2·2+4·2) |
| ACGTTGCAATGCCGTA | Wallace Tm | 48.0 (Biopython docstring example — reproduced) |

All independently derived from the cited formulas, not from repo output.

### Stage A findings

PASS-WITH-NOTES — two sourced, sibling-owned divergences (Marmur–Doty variant constants; mean vs
Trifonov-product complexity). Neither is an aggregation defect: the summary's contract is field
equality with the canonical methods, which holds.

## Stage B — Implementation

### Code path reviewed

- `src/.../Seqeron.Genomics.Analysis/SequenceStatistics.cs:990-1020` — `SummarizeNucleotideSequence`:
  null→`""`, calls the 4 canonical methods, builds the 6-key composition dict, assembles the record.
- Composition `:49-94`; Entropy `:732-762` (Math.Log2, bits); Linguistic complexity `:767-795`
  (mean of `observed/min(4^k, n-k+1)` over k=1..6); Tm `:569-590`; constants in
  `Seqeron.Genomics.Infrastructure/ThermoConstants.cs` (`WallaceMaxLength=14`, `CalculateWallaceTm`,
  `CalculateMarmurDotyTm` = 64.9/41/16.4).

### Formula realised correctly?

Yes. Each field is a direct copy of the canonical method's return on the same `seq`; the Tm flag is
`seq.Length < ThermoConstants.WallaceMaxLength` (14), matching the documented branch. No rounding,
no alphabet remap, no recomputation in the aggregator.

### Cross-verification table recomputed vs code

`dotnet test --filter SummarizeNucleotideSequence` → 12 passed, 0 failed. Every value in the Stage-A
cross-check table is asserted by the tests and passes (GcContent 0.5, Entropy 2.0, Wallace 24.0/12.0,
Marmur–Doty 43.375, Complexity 0.8396825396825397, RNA composition A=2/U=1/G=1/C=1/N=2/T=0).

### Variant / delegate consistency

M2 (16-mer) and S2 (lowercase) assert every field `== ` its canonical method — the aggregation
contract — and pass. No `*Fast`/instance variants exist for this method.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** M1/M3/M4/M5/M6 lock exact externally-derived values; M2/S2 assert the
  aggregation contract (field == canonical) which is the unit's actual obligation. Good.
- **Gap found & fixed this session:** the `Complexity` field had **no exact-value assertion** (M1 did
  not assert it; only M2/S2 checked it against the code itself). A deliberately-wrong complexity copy
  would have survived. **Fixed (test-only, 0 code change):** added an exact lock
  `summary.Complexity == 0.8396825396825397` to M1, derived independently from the documented
  vocabulary-usage-mean definition (hand-computed this session, not read from output).
- **No green-washing:** C1's `GreaterThan(0)/LessThan(1)` are correct for an *invariant/bound* test
  (INV-07 is a bound, not an exact value); no tolerance widened, nothing skipped.
- **Coverage:** all 6 fields, both Tm branches, RNA/N, empty, null, case-insensitivity, bounds — all
  exercised.
- **Honest green:** FULL unfiltered suite **6618 passed, 0 failed** (1 pre-existing benchmark skip,
  unrelated); `dotnet build` 0 errors (4 pre-existing warnings in unrelated test files).

### Stage B findings

PASS. One test-strengthening (A39, no code defect): added the missing exact Complexity lock.

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES (sibling-owned Marmur–Doty-variant & mean-complexity divergences; the
  aggregator's field-equality contract is sound and externally cross-checked).
- **Stage B:** PASS.
- **End-state:** ✅ CLEAN — aggregator delegates faithfully; tests strengthened with an exact
  externally-sourced Complexity value; full suite green.
- **Follow-ups:** none for this unit. The two notes are tracked by the SEQ-TM-001 and
  linguistic-complexity sibling units, not defects of the aggregation.
