# Validation Report: MOTIF-DISCOVER-001 — Motif Discovery via Overrepresented k-mers (O/E enrichment)

- **Validated:** 2026-06-16   **Area:** Matching / Motif Discovery
- **Canonical method(s):** `MotifFinder.DiscoverMotifs(DnaSequence sequence, int k = 6, int minCount = 2)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (one test-quality defect found and fixed in-session)

## Scope note

The unit *name* in the prompt ("de novo motif discovery — enumerative / median-string / greedy /
Gibbs") describes the classic `(l,d)` planted-motif family. The **registered** unit is **not** that:
per the Registry row, TestSpec, and Evidence, MOTIF-DISCOVER-001 is **overrepresented-k-mer
discovery in a single DNA sequence** using the observed/expected (O/E) ratio under a zero-order
i.i.d. uniform background. The classic motif-search algorithms (median string, greedy, Gibbs) are
not part of this unit (consensus/PWM-based units are separate: MOTIF-CONS-001, GENOMIC-MOTIFS-001,
etc.). Validation therefore targets the O/E-enrichment method actually registered.

## Stage A — Description

### Sources opened & what they confirm

1. **Compeau & Pevzner, *Bioinformatics Algorithms* (Ch. 2), reproduced at the wikiselev wiki**
   (https://github.com/wikiselev/bioinformatics-algorithms/wiki/Kmer-expected-number-of-occurrences-in-a-DNA-string,
   retrieved this session via WebFetch). Confirmed verbatim:
   - Background: "selecting each nucleotide (A, C, G, T) with the same probability (0.25)" — zero-order i.i.d. uniform.
   - Probability formula: `Pr(N,A,Pattern,t) ≈ ( N − t·(k−1) | t ) / A^(t·k)` (binomial over A^(t·k)).
   - R reference impl: `Pr <- function(N,A,k,t) choose(N - t*(k-1), t)/A^(t*k)`.
   - Worked example: `Pr(1000, 4, 9, 1)·500 ≈ 1.9`.
   - Self-overlap caveat: "only an approximation because it assumes that Pattern can not overlap
     with itself" — affects the *probability* statistic, **not** the deterministic count or O/E denominator.
2. **Independent web corroboration** (WebSearch, this session): "the expected number of random
   matches to a particular k-mer is (N−k+1)/4^k", with (N−k+1) = number of start positions and 4^k
   = number of possible k-mers; "calculate the expected versus observed ratios … and identify
   overrepresented sequences." Matches the Evidence/TestSpec claim independently of the textbook.

### Formula check

- Expected count of a *specific* k-mer = the `t = 1` case of the textbook formula:
  `choose(N − (k−1), 1) / 4^k = (N − k + 1) / 4^k`. Verified numerically:
  `choose(992,1) = 992 = 1000 − 9 + 1`. ✓
- O/E enrichment = observed count / E. A value > 1 means overrepresented (Sources 1 & 2). ✓
- Worked example reproduced: `Pr(1000,4,9,1)·500 = 1.89208… ≈ 1.9` (matches the wiki "≈ 1.9"). ✓

### Edge-case semantics check

- `k > N`: `N − k + 1 ≤ 0` ⇒ no length-k windows ⇒ empty result. Sourced (Evidence corner case). ✓
- Null sequence → `ArgumentNullException`; `k < 1` → `ArgumentOutOfRangeException` — method contract. ✓
- `minCount` is a *presentation filter* (Evidence Assumption 1): it never alters a returned record's
  Count/Positions/Enrichment, only which rows appear. Correct and sourced. ✓
- INV-04 (`E > 0` for any returned k-mer): true because a k-mer can only be counted when at least one
  window exists, i.e. `N − k + 1 ≥ 1`, so `E = (N−k+1)/4^k > 0`. ✓

### Independent cross-check (hand-computed, this session)

| Dataset | N | k | windows | E | k-mer | count | enrichment (exact) |
|---------|---|---|---------|---|-------|-------|--------------------|
| ATGCATGCATGC | 12 | 4 | 9 | 9/256 = 0.03515625 | ATGC @ {0,4,8} | 3 | 3/(9/256) = **768/9 ≈ 85.3333** |
| AAAAAAAAAA | 10 | 3 | 8 | 8/64 = 0.125 | AAA @ {0..7} | 8 | 8/0.125 = **64.0** |
| AACCGGAACCGG | 12 | 6 | 7 | 7/4096 | AACCGG @ {0,6} | 2 | 2/(7/4096) = **8192/7 ≈ 1170.286** |

All three independently recomputed from the external formula (Python, this session) — match the
Evidence and the test expectations.

### Findings / divergences

None. Description is mathematically correct and every non-trivial number traces to the retrieved
external sources. The "intentionally simplified" notes (zero-order background only; no p-value/E-value)
are honest and correct.

**Stage A verdict: PASS.**

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs:515-558` (`DiscoverMotifs`).

- One-pass window enumeration `for (i = 0; i <= seq.Length - k; i++)` builds a `Dictionary<string,
  List<int>>` of 0-based start positions (overlap allowed). ✓ INV-01.
- `windowCount = seq.Length - k + 1.0`; `expectedCount = windowCount / Math.Pow(4, k)`;
  `enrichment = positions.Count / expectedCount` — i.e. exactly `Count / ((N−k+1)/4^k)`. ✓ INV-02.
  No `max(E, 0.1)` floor (the previously-removed untraceable clamp); E is unclamped. ✓
- `positions.Count >= minCount` filter applied before yield. ✓ INV-03.
- Guards: `ArgumentNullException.ThrowIfNull(sequence)`; `if (k < 1) throw ArgumentOutOfRangeException`. ✓
- `k > N`: loop bound negative ⇒ no iterations ⇒ empty dictionary ⇒ empty result. ✓

### Cross-verification table recomputed vs code

The three hand-computed datasets above are reproduced exactly by the code (tests M1–M4, C1 pass with
`Within(1e-10)`/`1e-9`). 768/9, 64.0, and 8192/7 all verified against the code path.

### Variant/delegate consistency

`DiscoverMotifs` has a single canonical signature for this unit. The sibling `FindSharedMotifs`
(cross-sequence quorum) is a *different* registered unit (MOTIF-SHARED-001) and is out of scope here.

### Test quality audit (test-quality gate)

Tests file: `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_DiscoverMotifs_Tests.cs`.

| Test | Assessment |
|------|------------|
| M1 count=3 | Exact sourced value. ✓ |
| M2 positions {0,4,8} | Exact ordered set. ✓ |
| M3 enrichment 768/9 | Exact sourced value, tight tolerance. ✓ |
| M4 count=8 + enrichment 64.0 | Exact sourced values. ✓ |
| **M5 minCount filter** | **DEFECT (found & fixed).** Original used `"ATGCAAAA"`, k=4, where *every* 4-mer occurs exactly once, so with minCount=2 the result is **empty** and `motifs.All(m => m.Count >= 2)` passed **vacuously** — it would pass against an implementation that returns nothing or that ignores the filter entirely. Green-washing per the gate. |
| S1 null | Exact exception type. ✓ |
| S2 k<1 | Exact exception type. ✓ |
| S3 k>N empty | Real assertion. ✓ |
| C1 no-floor 8192/7 | Exact value; pins against clamp re-introduction. ✓ |

**Fix applied (this session):** M5 rewritten to use `"ACGTACGTAA"`, k=4. By hand enumeration the
4-mer multiset is `ACGT{0,4}`, `CGTA{1,5}` (count 2) and `GTAC{2}`, `TACG{3}`, `GTAA{6}` (count 1).
The rewritten test now asserts (a) INV-03, (b) the returned set is **exactly** `{ACGT, CGTA}`
(inclusion side), and (c) the three singletons are **excluded**. It now fails both a return-nothing
impl and a filter-ignoring impl. The expected k-mer multiset was derived independently by hand window
enumeration, not from code output.

### Findings / defects

- **D1 (test green-washing, FIXED-NOW):** M5 was a vacuously-passing filter test. Fixed in-session
  with a non-vacuous, hand-derived dataset. Full suite re-run green.

No implementation defect: the code faithfully realises the validated O/E formula and all edge cases.

**Stage B verdict: PASS-WITH-NOTES** (code correct; one test defect found and fully fixed).

## Verdict & follow-ups

- **Stage A: PASS.** **Stage B: PASS-WITH-NOTES.**
- **Test-quality gate: PASS after fix.** All assertions now use exact externally-sourced values; the
  one green-washed test (M5) was rewritten to a non-vacuous, hand-derived dataset. No assertion
  weakened, no tolerance widened, no test skipped. Full unfiltered suite: **6606 passed, 0 failed**;
  build 0 errors (4 pre-existing warnings in unrelated test files, none in MotifFinder).
- **End-state: ✅ CLEAN** — the single defect (vacuous M5 filter test) was completely fixed this
  session; the algorithm is fully functional and its tests lock the sourced values.
