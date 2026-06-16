# Validation Report: SEQ-ATSKEW-001 — AT Skew

- **Validated:** 2026-06-16   **Area:** Composition
- **Canonical method(s):** `GcSkewCalculator.CalculateAtSkew(string)` (canonical); `GcSkewCalculator.CalculateAtSkew(DnaSequence)` (delegate); private `CalculateAtSkewCore(string)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Scope note

The unit computes a **single scalar** AT skew = (A − T)/(A + T) over a whole sequence/window.
The session prompt mentioned a "per-window/cumulative AT-skew profile"; per the TestSpec and the
algorithm description this is **explicitly out of scope** for this unit — windowed/cumulative
profiles and origin/terminus prediction are provided by the GC-skew methods in the same class
(`CalculateWindowedGcSkew`, `CalculateCumulativeGcSkew`, `PredictReplicationOrigin`). The scalar
statistic is exactly what Lobry (1996) and Charneski (2011) define, so the scoping is sound.

## Stage A — Description

### Sources opened this session (URLs + extracted numbers)

1. **Wikipedia "GC skew"** — https://en.wikipedia.org/wiki/GC_skew (fetched 2026-06-16).
   Verbatim formulas: GC skew = **(G − C)/(G + C)**, AT skew = **(A − T)/(A + T)**.
   Range (verbatim): "The nucleotide composition skew spectra ranges from −1, which corresponds
   to G = 0 or A = 0, to +1, which corresponds to T = 0 or C = 0." ⇒ for AT skew, **−1 ⇔ A = 0**,
   **+1 ⇔ T = 0**. Primary attribution: Lobry, J. R. (1996), Mol Biol Evol 13(5):660–665.

2. **Charneski et al. (2011), PLoS Genetics** —
   https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002283 (fetched 2026-06-16).
   Verbatim: "excess of A over T in the leading strand, or positive AT skew given as **(A–T)/(A+T)**."
   Formula confirmed directly from a peer-reviewed primary source.

3. **Biopython `Bio.SeqUtils.GC_skew` source** —
   https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py (fetched 2026-06-16).
   Case handling (verbatim): `g = s.count("G") + s.count("g")`, `c = s.count("C") + s.count("c")` ⇒
   case-insensitive. Zero-denominator (verbatim): `try: skew = (g - c)/(g + c) except ZeroDivisionError:
   skew = 0.0`. Docstring (verbatim): "Does NOT look at any ambiguous nucleotides." These justify
   INV-03 (A+T=0 ⇒ 0), INV-04 (case-insensitive), INV-05 (non-A/T ignored) **by analogy** to the
   shipped GC-skew code (Biopython ships no AT-skew-specific line).

### Formula check

`(A − T)/(A + T)` matches Charneski (2011) verbatim and Wikipedia citing Lobry (1996) — exact, no
divergence in symbols, normalisation, or sign convention.

### Edge-case semantics check

- Bounds INV-01: |A − T| ≤ A + T ⇒ result ∈ [−1, +1]; +1 ⇔ T = 0, −1 ⇔ A = 0 — matches Wikipedia/Lobry.
- INV-02: A = T (A+T>0) ⇒ numerator 0 ⇒ 0 — arithmetic.
- INV-03: A + T = 0 ⇒ 0 (no NaN/throw) — Biopython ZeroDivisionError → 0.0 convention.
- INV-04 / INV-05: case-insensitive; non-A/T ignored — Biopython.
- Null/empty string ⇒ 0; null DnaSequence ⇒ ArgumentNullException — documented validation, reasonable.

### Independent cross-check (hand-computed from the sourced formula)

| Input | A | T | (A−T)/(A+T) | Expected |
|-------|---|---|-------------|----------|
| `AAAA` | 4 | 0 | 4/4 | **1.0** |
| `TTTT` | 0 | 4 | −4/4 | **−1.0** |
| `ATAT` | 2 | 2 | 0/4 | **0.0** |
| `AAAT` | 3 | 1 | 2/4 | **0.5** |
| `ATTT` | 1 | 3 | −2/4 | **−0.5** |
| `GGCC` | 0 | 0 | 0/0→0 | **0.0** |
| `AAATGGGCCC` | 3 | 1 | 2/4 | **0.5** |
| `AAATTGCGCAATA` | 6 | 3 | 3/9 | **0.3333… (1/3)** |

All values are arithmetic consequences of the sourced (A−T)/(A+T) definition; none trace to repo output.

### Findings / divergences

None. The description, formula, range, invariants, and edge-case semantics are all confirmed against
authoritative sources. The one documented assumption (lowercase/non-A/T handling borrowed from the
directly analogous Biopython `GC_skew`) is reasonable and explicitly disclosed; the formula itself is
fully primary-sourced. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs:174-206`.
- `CalculateAtSkew(string)`: null/empty ⇒ 0; else `ToUpperInvariant()` then core.
- `CalculateAtSkew(DnaSequence)`: `ArgumentNullException.ThrowIfNull`; forwards normalized `.Sequence`.
- `CalculateAtSkewCore`: counts 'A' and 'T', `total = aCount + tCount`, returns
  `total > 0 ? (double)(aCount - tCount)/total : 0`.

### Formula realised correctly?

Yes — `(double)(aCount − tCount)/total` is exactly (A − T)/(A + T). The `total > 0` guard implements
the Biopython zero-denominator → 0 convention (INV-03). `ToUpperInvariant` gives INV-04. Counting only
'A'/'T' gives INV-05. The `(double)` cast prevents integer division. No precision/overflow risk on
realistic input.

### Cross-verification table recomputed vs code

Every row of the Stage-A table is asserted by a passing test (M1–M7, S1, plus the strengthened S5/C1).
The 1/3 mixed case was recomputed independently (A=6, T=3 ⇒ 0.3333333333333333) and matches the code.

### Variant/delegate consistency

`DnaSequence` overload forwards to the same core; C1 confirms it equals the string overload **and**
(after strengthening) equals the exact sourced value 0.5 on `AAATGGGCCC`.

### Test quality audit (HARD gate)

- **Sourced, not code-echoed:** M1–M7 assert exact values traceable to the formula/range sources, with
  `Is.EqualTo(x).Within(1e-10)` — a wrong-sign or off-by-one implementation would fail them.
- **No greenwashing:** no exact-value case uses Greater/AtLeast/Contains. Two weaknesses were found and
  **fixed this session**:
  - **S5** previously asserted only the range bound (`>= -1`, `<= 1`) — that would pass against many
    wrong implementations. Rewritten to also assert the exact sourced value **1/3** for
    `AAATTGCGCAATA` (renamed `CalculateAtSkew_MixedSequence_ReturnsExactValueWithinBounds`).
  - **C1** previously only asserted `viaDna == viaString` — both could be identically wrong. Strengthened
    to also assert the exact sourced value **0.5**.
- **Coverage:** both public overloads exercised; all five invariants (INV-01..05), all Stage-A
  edge/error cases (pure A, pure T, balanced, ±asymmetric, no-A/T zero-denominator, G/C-ignored,
  lowercase, null string, empty string, null DnaSequence, range bound, delegate equivalence) covered.
- **Honest green:** FULL unfiltered suite **Failed: 0, Passed: 6607**; changed test file builds
  warning-free (the 4 build warnings are pre-existing in an unrelated file, `ApproximateMatcher_EditDistance_Tests.cs`).

### Findings / defects

No implementation defect. Two **test-quality defects** (S5 range-only, C1 echo-only) found and fully
fixed in this session by locking exact sourced values. **Stage B: PASS.**

## Verdict & follow-ups

- **Stage A: PASS** — formula, range, invariants, edge semantics confirmed against Charneski (2011),
  Wikipedia/Lobry (1996), and Biopython.
- **Stage B: PASS** — code realises (A−T)/(A+T) exactly with correct zero-denominator/case/symbol
  handling; tests now lock exact sourced values for every case.
- **Test-quality gate:** PASS (after fixing S5 and C1; full suite green, warning-free build of changed files).
- **End-state: CLEAN.** No code change needed; the two test weaknesses were completely fixed.
