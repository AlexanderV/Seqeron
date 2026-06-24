# Validation Report: SEQ-ENTROPY-001 — Shannon Entropy

- **Validated:** 2026-06-24   **Area:** Composition
- **Canonical method(s):** `SequenceComplexity.CalculateShannonEntropy(DnaSequence)`,
  `SequenceComplexity.CalculateShannonEntropy(string)`,
  `SequenceComplexity.CalculateKmerEntropy(DnaSequence, k)`
  (wrapper smoke: `SequenceStatistics.CalculateShannonEntropy(string)` — separate impl, all-letters alphabet)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Entropy (information theory)** (fetched live, not from memory). Confirmed verbatim:
  - Formula: `H(X) := −Σ_{x∈X} p(x) log p(x)`.
  - 0·log0 convention: "In the case of p(x)=0 …, the value of the corresponding summand
    0 log_b(0) is taken to be 0, which is consistent with the limit lim_{p→0+} p log(p) = 0."
  - Maximum: "The maximal entropy of an event with n different outcomes is log_b(n):
    it is attained by the uniform probability distribution."
  - Units: "Base 2 gives the unit of bits (or 'shannons')."
- **Shannon (1948), A Mathematical Theory of Communication** — original definition of H and the
  0·log0 = 0 limit convention (cited in spec; consistent with the Wikipedia derivation above).
- DNA specialisation (n = 4): H_max = log2(4) = **2 bits**.

### Formula check
Spec formula `H = −Σ p_i log2 p_i` (bits, base 2) matches the authoritative source **exactly**:
log base 2 → bits, 0·log0 = 0 convention, max = log_b(n) at uniform. Invariants INV-ENT-001..005
are all genuine mathematical properties (range 0..2 for DNA; max at uniform = log2(4) = 2;
homopolymer = 0; k-mer bound log2(4^k) = 2k; empty sum = 0 via the expansibility axiom).

### Edge-case semantics check
Empty → 0, single symbol/homopolymer → 0, uniform → max, non-DNA chars excluded from both
numerator and denominator (alphabet fixed to {A,T,G,C}), null guard, invalid k guard,
length < k → 0 — all defined and sourced in the spec.

### Independent cross-check (hand-recomputed, exact, this session)
| Input | Distribution | H computed |
|-------|--------------|------------|
| `ATGCATGCATGCATGC` (uniform 4) | 4,4,4,4 | **2.0** |
| `AAAAAAA` (homopolymer) | 7 | **0.0** |
| `ATATATAT` (50/50) | 4,4 | **1.0** |
| `ATGATGATG` (3 uniform) | 3,3,3 | **1.5849625…** = log2(3) |
| `ATGCNN` (N ignored) | 1,1,1,1 | **2.0** |
| k-mer `ATCG` k=2 (AT,TC,CG) | 1,1,1 | **1.5849625…** = log2(3) |
| k-mer `ATGCATGCATGCATGC` k=2 (AT4,TG4,GC4,CA3 /15) | 4,4,4,3 | **1.9898981…** |

All values reproduce the spec's cross-verification table exactly.

### Findings / divergences
None. The spec's k=2 `ATGCATGCATGCATGC` table value reads **≈ 1.98990** (a prior session had
already corrected an earlier 1.98082 typo); my independent recomputation gives 1.9898981…, which
rounds to 1.98990 — confirmed correct. No description changes needed this session.

---

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs`
  - `CalculateShannonEntropyCore` (l.93–120): fixed alphabet {A,T,G,C}, counts only those chars,
    `total==0 → 0`, `count>0` guard implements 0·log0 = 0, `entropy -= p*Math.Log2(p)` → **bits**.
  - `CalculateKmerEntropyCore` (l.164–189): overlapping windows (`Substring(i,k)`, step 1),
    `length<k → 0`; since `length≥k ⇒ total≥1`, the k-mer loop never divides by zero.
  - Guards: `ArgumentNullException.ThrowIfNull` (l.80,142); `k<1 → ArgumentOutOfRangeException`
    (l.143,159). String overloads short-circuit null/empty → 0 (l.89,160).
  - String overload uppercases via `ToUpperInvariant()` (l.90,161); DnaSequence overload feeds
    `sequence.Sequence` (already canonical uppercase) — both reach the same Core.

### Formula realised correctly?
Yes. `entropy -= p * Math.Log2(p)` over `p = count/total` is exactly `−Σ p log2 p` in bits, with
0·log0 handled by the `count>0` / `total==0` guards — both match the Stage-A validated description.

### Cross-verification table recomputed vs code
All Stage-A values are asserted by the canonical tests (exact `2.0`/`0.0`/`1.0`, `log2(3)` within
1e-10, the k=2 formula within 1e-10). The 51 `SequenceComplexityTests` pass, including all entropy cases.

### Variant/delegate consistency
String vs DnaSequence overload bitwise-equal (S02 test). `SequenceStatistics` wrapper is a distinct
implementation (counts all letters, not a fixed DNA alphabet) and is correctly scoped "smoke only";
its 3 smoke tests (uniform→2, homopolymer→0, empty→0) pass.

### Numerical robustness
`p ∈ (0,1]` always; `Math.Log2` is exact at powers of two; no overflow/div-by-zero on stated ranges
(homopolymer returns IEEE `-0.0`, value-equal to `0.0` and accepted by `Is.EqualTo(0.0)`).

### Test quality audit
Assertions check exact sourced values (not "no throw"): exact `2.0`/`0.0`/`1.0`, `log2(3)` and k-mer
formula within 1e-10, range invariant 0..2. Edge cases (empty, null, k<1, length<k, homopolymer,
non-DNA) all covered. Deterministic.

### Findings / defects
None in code or tests.

---

## Verdict & follow-ups
- **Stage A: PASS** — formula, log base (bits), 0·log0 convention, and max-at-uniform all match the
  authoritative Wikipedia/Shannon sources opened this session; all cross-check numbers reproduced by hand.
- **Stage B: PASS** — code faithfully realises `H = −Σ p log2 p` in bits with correct 0·log0 and
  edge handling; every worked example reproduces.
- **End state: CLEAN.** Build succeeds (0 warn/err); `SequenceComplexityTests` 51/51 pass.
- Files changed: none (no source, test, or spec change required this session).
</content>
</invoke>
