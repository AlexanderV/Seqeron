# Validation Report: SEQ-COMPLEX-KMER-001 — K-mer Entropy

- **Validated:** 2026-06-16   **Area:** Complexity (Extended Sequence Complexity)
- **Canonical method(s):** `SequenceComplexity.CalculateKmerEntropy(DnaSequence, int k = 2)`; delegate `CalculateKmerEntropy(string, int k = 2)`; private core `CalculateKmerEntropyCore(string, int)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Scope note (prompt vs. registry)

The session prompt described this unit as *"linguistic complexity = observed distinct k-mers / maximum
possible distinct k-mers"* (Trifonov / Gabrielian–Bolshoy). That is a **different** algorithm. The repo
defines SEQ-COMPLEX-KMER-001 consistently across the Registry (`ALGORITHMS_CHECKLIST_V2.md` §50:
"K-mer Entropy", canonical `SequenceComplexity.CalculateKmerEntropy`), the TestSpec, the Evidence, the
algorithm doc, and the code as **Shannon entropy (in bits) of the overlapping k-mer frequency
distribution**. (Linguistic complexity is a separate `SequenceComplexity.CalculateLinguisticComplexity`
method / different unit.) This report validates the unit **as the repo defines it** (K-mer Entropy);
the prompt's "linguistic complexity" framing does not apply here.

## Stage A — Description

### Sources opened & what they confirm (re-derived this session)

- **Shannon entropy formula & bounds** — standard, textbook. `H = −Σ pᵢ log₂ pᵢ`; `0 ≤ H ≤ log_b(n)`;
  `H = 0` for a deterministic distribution (one pᵢ = 1); `H = log_b(n)` for the uniform distribution over
  `n` outcomes (Shannon 1948; Wikipedia "Entropy (information theory)"). These are well-established
  mathematical facts; I re-derived the four invariants and the worked numbers independently below.
- **Overlapping k-mers convention** — a length-L sequence has `N = L − k + 1` overlapping k-mers
  (sliding window, step 1); `pᵢ = nᵢ / N`. This is the longdust (Li 2025) k-mer-frequency formulation
  the Evidence cites, and is the universal overlapping-k-mer count.
- **Log base 2 → bits**; single-nucleotide maximum `= log₂ 4 = 2 bits`.

### Formula check

The documented model `H = −Σ pᵢ log₂ pᵢ`, `pᵢ = nᵢ/N`, `N = L−k+1`, in bits, matches the cited
Shannon/longdust formulation exactly: correct symbols, correct normalisation (Σpᵢ=1), correct log base
(2 → bits), overlapping-window count. No divergence.

### Edge-case semantics

- `L < k` (no k-mers) → 0 (entropy of empty multiset; documented assumption, consistent with sibling
  methods). Mathematically sound.
- Single distinct k-mer (deterministic) → H = 0 (Shannon certainty).
- All-distinct k-mers (uniform) → H = log₂(N) (Shannon uniform maximum).
- `k < 1` → `ArgumentOutOfRangeException`; null `DnaSequence` → `ArgumentNullException`; null/empty
  string → 0 — these are API-contract choices (not entropy literature), reasonable and consistent with
  siblings. They do not affect entropy values for valid input.

### Independent cross-check (numbers — computed this session, NOT from the implementation)

Recomputed every expected value with an independent Python reference
(`Counter` over overlapping k-mers, `-Σ (n/N)·math.log2(n/N)`):

| Input | k | counts | N | H (independent) | spec value |
|-------|---|--------|---|-----------------|------------|
| `ACGT` | 1 | A,C,G,T each 1 | 4 | 2.0 | 2.0 ✓ |
| `ACGT` | 2 | AC,CG,GT each 1 | 3 | 1.584962500721156 = log₂3 | 1.5849625007211562 ✓ |
| `ATATAT` | 2 | AT=3,TA=2 | 5 | 0.9709505944546686 | 0.9709505944546686 ✓ |
| `AAAA` | 2 | AA=3 | 3 | 0.0 | 0.0 ✓ |
| `AAACGT` | 2 | AA=2,AC=1,CG=1,GT=1 | 5 | 1.921928094887362 = log₂5−0.4 | 1.9219280948873623 ✓ |
| `AC` | 5 | (L<k) | — | 0.0 | 0.0 ✓ |

Bounds invariant `0 ≤ H ≤ log₂N` independently verified for all four C1 cases
(`ACGTACGTAA`/2 → 2.281≤3.170; `AAAAAAAA`/3 → 0; `ACGTACGTACGT`/4 → 1.975≤3.170; `GCGCGCGCGCAT`/1 →
1.650≤3.585).

### Findings / divergences

None. **Stage A: PASS.** The description is mathematically correct and matches authoritative
definitions; every numeric expectation traces to an independent computation done this session.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs:140–189`.

- `CalculateKmerEntropy(DnaSequence,k)` (l.140): null-check, `k<1` guard, delegates core with
  `sequence.Sequence`.
- `CalculateKmerEntropy(string,k)` (l.157): **own** `k<1` guard, null/empty → 0, upper-cases, delegates.
- `CalculateKmerEntropyCore(string,k)` (l.164): `L<k`→0; single linear scan building a
  `Dictionary<string,int>` of overlapping k-mer counts (`total = L−k+1`); `H = −Σ (count/total)·log₂(count/total)`.

### Formula realised correctly?

Yes — verbatim `H = −Σ pᵢ log₂ pᵢ`, `pᵢ = count/total`, `total = L−k+1`, `Math.Log2` (base 2 → bits).
Overlapping windows enumerated by `seq.Substring(i,k)` for `i = 0..L−k`. For p=1, `Math.Log2(1)=0` so a
single-k-mer sequence yields exactly 0.0. No precision loss, overflow, or div-by-zero on stated ranges.

### Cross-verification table recomputed vs code

All 15 canonical tests pass with `.Within(1e-10)` against the independently-sourced values in the Stage A
table (run this session): M1 2.0, M2 log₂3, M3 0.9709505944546686, M4 0.0, M5 log₂5−0.4, M6 0.0, plus the
four C1 bounds rows. Code output == independent reference for every case.

### Variant/delegate consistency

String overload delegates to the same core after `ToUpperInvariant()`; S1 confirms
`string("ATATAT") == DnaSequence("ATATAT")` exact value, S2 confirms lowercase == uppercase. Consistent.

### Test quality audit

Original canonical file (12 named cases / 15 with C1 rows) used **exact sourced values** for all
M-tests (no Greater/AtLeast/Contains on the happy path), exact throws for M7/M8, and an appropriate
bounds **property** assertion for the C1 invariant. Sound. Two **coverage gaps** found:

1. The **string overload's own `k < 1` guard** (l.159) was untested — M7 only exercises the
   *DnaSequence* overload's guard (a distinct code path).
2. The **string-overload `L < k` (non-empty, too short)** core branch was untested — M6/S3 cover only
   the DnaSequence path and null/empty string, never a non-empty short string via the string overload.

A wrong string-overload guard (e.g. missing/mis-ordered `k<1` check) would have passed the whole suite.

### Findings / defects

No code defect. Two **test-coverage** gaps (test-only fix, 0 code change):
- Added **S4** `CalculateKmerEntropy_StringOverload_InvalidK_Throws` (`"ACGT"`, k=0 →
  `ArgumentOutOfRangeException`) — locks the string overload's own guard.
- Added **S5** `CalculateKmerEntropy_StringOverload_SequenceShorterThanK_ReturnsZero` (`"AC"`, k=5 → 0.0)
  — exercises the `L<k` core branch via the string path; value independently confirmed (empty-multiset
  entropy = 0).

Fixture 15 → 17 assertions.

## Verdict & follow-ups

- **Stage A: PASS** — description mathematically correct, matches Shannon/longdust definitions,
  all numbers independently reproduced this session.
- **Stage B: PASS** — code realises the formula verbatim; all values match the independent reference;
  variants consistent.
- **Test-quality gate: PASS** (after fix) — exact sourced values throughout, no green-washing, all
  public overloads and Stage-A branches now covered (incl. both overloads' `k<1` guards and `L<k`
  paths). Full unfiltered suite **6598 passed / 0 failed** (1 benchmark skipped by design); changed test
  file builds warning-free.
- **End-state: CLEAN.** No code defect; the two test-coverage gaps were completely fixed in-session.
- **Note:** the session prompt's "linguistic complexity = distinct/max" framing describes a different
  algorithm; this unit is K-mer Entropy, validated as the repo defines it.
