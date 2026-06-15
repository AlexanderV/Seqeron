# Validation Report: KMER-DIST-001 — K-mer Euclidean Distance

- **Validated:** 2026-06-15   **Area:** K-mer / Alignment-free sequence comparison
- **Canonical method(s):** `KmerAnalyzer.KmerDistance(string seq1, string seq2, int k)`
  (supported by `KmerAnalyzer.GetKmerFrequencies` and `KmerAnalyzer.CountKmers`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (retrieved live)

| # | Source | URL | What it confirmed |
|---|--------|-----|-------------------|
| 1 | Zielezinski et al. (2017), Genome Biology 18:186, Fig. 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5627421/ | Verbatim worked example: x="ATGTGTG", y="CATGTG", k=3; W_X3={ATG,TGT,GTG}, W_Y3={CAT,ATG,TGT,GTG}; union W3={ATG,CAT,GTG,TGT}; count vectors c_X3=(1,0,2,2), c_Y3=(1,1,1,1); "distance function … typically using Euclidean distance"; "two identical sequences will result in a distance of 0." |
| 2 | Lau et al. (2022), NAR Genom Bioinform | https://pmc.ncbi.nlm.nih.gov/articles/PMC9442500/ | Verbatim: "The k-mer frequencies are derived from the counts by dividing each k-mer count by the total number of k-mers in the sequence (i.e. the sequence length minus the k-mer length)." and "various distance measures can be used … include the Euclidian, Manhattan, Canberra or Chebyshev distances." |
| 3 | Wikipedia — K-mer | https://en.wikipedia.org/wiki/K-mer | "a sequence of length L will have L − k + 1 k-mers" (number of overlapping windows). |
| 4 | scikit-bio `Sequence.kmer_frequencies` | https://scikit.bio/docs/latest/generated/skbio.sequence.Sequence.kmer_frequencies.html | Independent reference implementation: relative frequency = count ÷ total number of k-mers (sum of counts). Worked example 'ACACATTTATTA', k=3, non-overlap: counts {ACA:1,CAT:1,TTA:2} → relative {ACA:0.25,CAT:0.25,TTA:0.5}. |

### Formula check

- Word-vector model + union over both sequences with 0 for absent words — matches Source 1 Fig. 1 exactly (the spec's union, word order, and count vectors are reproduced verbatim from the source).
- Frequency normalization f_s(w) = count ÷ (number of k-mer windows) — matches Source 2 (count ÷ total k-mers) and Source 4 (count ÷ sum of counts).
- Euclidean distance d = √(Σ (f_x(w) − f_y(w))²) over the union — matches Source 1 ("Euclidean distance") applied to relative-frequency vectors (Source 2; spec Source 4 Boden 2014).
- INV-1…INV-4 are genuine mathematical properties of the Euclidean metric on frequency vectors; INV-4 (√2 for disjoint single-k-mer sequences) is confirmed by hand and by the independent Python recomputation below.

### Note (the "minor" in PASS-WITH-NOTES)

Source 2 (Lau et al.) parenthetically equates "the total number of k-mers in the sequence" with
"the sequence length minus the k-mer length" (= L − k). The established and correct count of
overlapping k-mer windows is **L − k + 1** (Source 3, Wikipedia; consistently across the
literature). This is an off-by-one imprecision in the Lau prose. The Seqeron implementation does
**not** inherit that error: it divides each count by `counts.Values.Sum()`, which equals the true
window count L − k + 1, and this is exactly what reproduces the Fig. 1 count vectors from Source 1
and the scikit-bio convention from Source 4. So the description is correct as implemented; the only
divergence is in one cited source's wording, which is documented here and does not affect any
numeric expectation. Logged BY-DESIGN in FINDINGS_REGISTER (the implementation uses the correct
standard divisor).

### Edge-case semantics

- Identical ⇒ 0 (Source 1). Absent word ⇒ 0 component (Source 1, c_X has 0 for CAT). Frequencies
  require L ≥ k (Source 2). L < k / empty / null are not defined by any source; the implementation's
  "empty = zero vector" is a sound, explicitly-flagged ASSUMPTION (A2), and is the natural extension.

### Independent cross-check (numbers recomputed from scratch this session, not from the repo code)

A standalone Python reimplementation (Counter-based k-mer counting, count ÷ sum, Euclidean over the
key union) reproduced every spec value:

| Case | Inputs | Independent result | Spec expectation |
|------|--------|--------------------|------------------|
| M1 | "ATGTGTG","CATGTG",k=3 | 0.33166247903554 | √0.11 = 0.33166247903553997 ✓ |
| M2 | identical, k=3 | 0.0 | 0.0 ✓ |
| M3 | "AAAA","AAAT",k=1 | 0.3535533905932738 | √0.125 ✓ |
| M5/S1 | "AAAA","TTTT",k=2 / "CCCCC","GGGGG",k=2 | 1.4142135623730951 | √2 ✓ |
| S2 | "ACGT","AAAAAA",k=5 | 1.0 | 1.0 ✓ |
| S3 | "atgtgtg","CATGTG",k=3 | 0.33166247903554 | √0.11 ✓ |
| C2 | "","",k=3 | 0.0 | 0.0 ✓ |

The frequency variant (spec's choice) is the one endorsed by Sources 2 and 4; the count-based
variant from Source 1's raw vectors gives √3 ≈ 1.7320508 and is correctly documented as the
alternative, not the implemented behaviour.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs`
- `KmerDistance` (lines ~213–240): validates k>0, builds freq1/freq2 via `GetKmerFrequencies`,
  unions the key sets, sums squared per-word differences (missing key ⇒ 0), returns `Math.Sqrt`.
- `GetKmerFrequencies` (~175–195): `count ÷ Σcounts`; returns empty dict when Σ=0.
- `CountKmers` (~20–41): null/empty ⇒ empty; `k≤0` ⇒ ArgumentOutOfRangeException; `k>len` ⇒ empty;
  upper-cases input.

### Formula realised correctly?

Yes. The divisor `counts.Values.Sum()` equals the overlapping window count L − k + 1 for every
sequence whose positions all form k-mers (always true for contiguous k-mers). The union-of-keys
sparse computation is mathematically identical to the dense 4^k vector since unobserved words are 0
on both sides. The double-validation of `k≤0` (in both `KmerDistance` and `CountKmers`) is harmless.

### Cross-verification vs code (full suite executed)

All 11 tests in the canonical fixture pass; their expected constants equal the independently
recomputed values in the Stage-A table. No code value was used to derive an expectation.

### Variant/delegate consistency

`KmerDistance` is the single public entry point for this unit; no `*Fast`/instance variant exists.
`GetKmerFrequencies`/`CountKmers` are exercised transitively and behave as specified.

### Test quality audit (HARD gate)

- **Sourced, not echoed:** every assertion is `Is.EqualTo(<exact sourced constant>).Within(1e-10)`
  — √0.11, √0.125, √2, 1.0, 0.0 — each traceable to Source 1/2 + the independent recompute. A
  deliberately-wrong implementation (e.g. count-based → √3, or wrong divisor) would fail these. No
  Greater/AtLeast/Contains used as a primary assertion (M5's `GreaterThanOrEqualTo(0)` is a
  secondary INV-03 check alongside the exact √2 assertion).
- **No green-washing:** the three pre-existing permissive tests (tolerance-only identity,
  `GreaterThan(0)`, `LessThan` ordering) were already removed from `KmerAnalyzerTests.cs` and
  replaced with exact-value tests in the dedicated file. No tolerance widened, nothing skipped.
- **Coverage:** all Stage-A branches/edges covered — worked example, identity, k=1 derivation,
  symmetry, disjoint √2, short-input (L<k), case-insensitivity, invalid-k throw, both-empty.
- **Gap closed this session:** the contract documents `null` inputs as allowed, but no test covered
  null. Added **C3 `KmerDistance_NullInput_TreatedAsZeroVector`** (null vs "AAAAA" k=3 ⇒ 1.0;
  null vs null ⇒ 0.0), with both expectations independently recomputed in Python (1.0 and 0.0).
- **Honest green:** full unfiltered suite = **Failed: 0, Passed: 6570** (was 6569; +1 for C3).
  Build: 0 errors; target test file builds warning-free (4 pre-existing warnings are in unrelated
  ApproximateMatcher files).

### Findings / defects

No code defect. One source-wording note (Stage-A, BY-DESIGN). One test-coverage gap (null) found
and fixed in-session.

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — formula, worked example, and invariants confirmed against primary
  literature; lone note is an off-by-one in one source's prose that the implementation does not share.
- **Stage B:** PASS — code faithfully realises the validated frequency-variant Euclidean distance;
  tests now cover all branches with exact sourced values including null.
- **Test-quality gate:** PASS.
- **End-state:** ✅ CLEAN (no defect; the one coverage gap completely fixed and locked with a sourced test).
