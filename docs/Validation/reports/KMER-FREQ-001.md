# Validation Report: KMER-FREQ-001 — K-mer Frequency Analysis

- **Validated:** 2026-06-24   **Area:** K-mer Analysis
- **Canonical method(s):** `KmerAnalyzer.GetKmerSpectrum(string, int)`, `KmerAnalyzer.GetKmerFrequencies(string, int)`, `KmerAnalyzer.CalculateKmerEntropy(string, int)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "K-mer"** (fetched 2026-06-24): a sequence of length L contains **"L − k + 1"** k-mers; there exist **"n^k"** possible k-mers (n = alphabet size, four for DNA); the k-mer spectrum shows **"the multiplicity of each k-mer in a sequence versus the number of k-mers with that multiplicity"** — i.e. a frequency-of-frequencies distribution. Exactly the spectrum the TestSpec/code describe.
- **Rosalind "k-Mer Composition" (KMER)** (fetched 2026-06-24): k-mer composition = "an array encoding the number of times that each possible k-mer occurs in the string," extracted by a sliding window (each start position except the last k−1). Confirms the count-over-sliding-window model and the L−k+1 window count (worked example CTTCGAAAGTTT → 9 4-mers, L=12, 12−4+1=9).
- **Shannon entropy** convention (relied on from established theory and prior report): H = −Σ p log₂ p, base 2 → bits; 0·log0 = 0; max = log₂(n) for n equiprobable outcomes; 0 under certainty.

### Formula check
- Relative frequency of a k-mer = count / (total k-mers) = count / (L − k + 1). Since Σ counts over all distinct k-mers = L − k + 1 (every window position contributes exactly one occurrence), the frequencies form a probability distribution summing to **1.0**. No pseudocount is applied (none warranted for a relative-frequency spectrum).
- k-mer Shannon entropy: H = −Σ f·log₂ f over the observed k-mers; base 2 (bits); 0·log0 = 0.

### Edge-case semantics
- Empty sequence or k > L → L − k + 1 ≤ 0 → no k-mers → empty frequency/spectrum maps; entropy = 0. Defined and consistent.
- Single distinct k-mer (homopolymer or k = L) → frequency 1.0, spectrum {count:1}, entropy 0.
- Uniform distribution of n distinct k-mers → entropy = log₂(n) (maximum).

### Independent cross-check (hand computation, reproduced in Python)
- **"ACGTACGT", k=2** (L=8 ⇒ 7 windows): counts AC=2, CG=2, GT=2, TA=1; Σ freq = 7/7 = **1.0 ✓**; H = **1.950212 bits**.
- **"AAACGT", k=2**: AA=2/5=0.4, AC=CG=GT=1/5=0.2 → Σ=1.0; H = log₂(5) − 0.4 = **1.921928 bits** (matches test S3 exactly).
- **"ACGT", k=1**: uniform 4 bases → H = log₂(4) = **2.0 bits**.
- **"ACGT", k=2**: 3 distinct 2-mers → H = log₂(3) = **1.584963 bits**.
- **"AAAA", k=2**: single k-mer → H = **0**.
- Spectra: "ACGTACGT" k=4 → **{1:3, 2:1}**; "AACGAACG" k=2 → **{1:1, 2:3}**; "AAAAA" k=2 → **{4:1}**.
- Spectrum total invariant: "ACGTACGT" k=3 → {1:2, 2:2}, Σ(mult×count)=6 = L−k+1 ✓.
- Edge: freq("",2)=∅, freq("ACG",4)=∅, ent("ACG",5)=0 ✓.

### Findings / divergences
None. Description matches authoritative sources exactly.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs`
- `CountKmers` (L20–42): sliding window i = 0 … L−k, `ToUpperInvariant` normalization → counts. Empty/null → empty dict; k≤0 → throws; k > L → empty dict.
- `GetKmerFrequencies` (L177–189): `total = counts.Values.Sum()` (= L−k+1); freq = count/total; total==0 → empty dict.
- `GetKmerSpectrum` (L136–148): inverts counts → {multiplicity : number_of_kmers_with_that_multiplicity}.
- `CalculateKmerEntropy` (L331–346): H = −Σ f·Math.Log2(f), guarded by `if (freq > 0)` (implements 0log0=0); empty frequencies → 0.

### Formula realised correctly?
- Denominator = Σ counts = L − k + 1 (each window contributes exactly one occurrence) → sum-to-1 invariant holds. ✓
- Entropy uses `Math.Log2` (bits) with the `freq > 0` guard for the 0log0 convention. ✓
- Edge cases (empty, k>L) return empty maps / 0 entropy as specified. ✓

### Cross-verification table recomputed vs code (32 passing tests + Python reference)
| Input | Quantity | Expected (source/hand) | Code result |
|-------|----------|------------------------|-------------|
| ACGTACGT, k=2 | Σ frequencies | 1.0 | 1.0 ✓ |
| AAACGT, k=2 | AA / AC,CG,GT | 0.4 / 0.2 | match ✓ |
| ACGTACGT, k=4 | spectrum | {1:3, 2:1} | match ✓ |
| AACGAACG, k=2 | spectrum | {1:1, 2:3} | match ✓ |
| AAAAA, k=2 | spectrum | {4:1} | match ✓ |
| AAAA, k=2 | entropy | 0.0 | 0.0 ✓ |
| ACGT, k=1 | entropy | log₂4 = 2.0 | 2.0 ✓ |
| ACGT, k=2 | entropy | log₂3 ≈ 1.58496 | match ✓ |
| AAACGT, k=2 | entropy | log₂5 − 0.4 ≈ 1.92193 | match ✓ |

### Variant/delegate consistency
`CalculateKmerEntropy` delegates to `GetKmerFrequencies`; `KmerDistance` and `AnalyzeKmers` build on `GetKmerFrequencies`/`CountKmers` — no divergent reimplementation. The `CountKmers` overloads (cancellation, DnaSequence, Span) share the same window/normalization logic.

### Test quality audit
`tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_Frequency_Tests.cs` (32 tests, all passing). Assertions check exact sourced values (frequencies, spectra, entropies) plus invariants (sum-to-1, Σ(mult×count)=L−k+1, 0 ≤ H ≤ log₂n) and edge cases (empty, k>L) and case-insensitivity (full key/value comparison). Deterministic and source-backed; no tautology/"no-throw" weakness.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS**, **State: CLEAN.**
- No code changes required. Filtered run: 32 passed / 0 failed (`KmerAnalyzer_Frequency_Tests`). Implementation unchanged since prior validation (cb113ce); independently re-confirmed against Wikipedia + Rosalind and hand/Python computation.
