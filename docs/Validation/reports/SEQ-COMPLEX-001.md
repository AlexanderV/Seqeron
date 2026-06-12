# Validation Report: SEQ-COMPLEX-001 — Sequence Complexity Metrics

- **Validated:** 2026-06-12   **Area:** Sequence Composition
- **Canonical method(s):** `SequenceComplexity.CalculateLinguisticComplexity` (+ Shannon/k-mer entropy, windowed complexity, low-complexity regions, DUST score, masking, compression ratio)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

Source file: `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs`
Test file: `tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexityTests.cs` (91 tests)

---

## Stage A — Description

### Sources opened
- **Wikipedia, "Linguistic sequence complexity"** (WebFetch). Confirms vocabulary usage
  `U_i = actual vocabulary size / maximum possible vocabulary size`, with max possible =
  `min(4^i, N − i + 1)`. The Wikipedia article presents the **product** form
  `C = U₁·U₂·…·U_w`. It gives the worked example `ACGGGAAGCTGATTCCA` (U₂=14/16, U₃=15/15,
  U₄=14/14) and the dinucleotide repeat `ACACACACACACACACA` (U₁=1/2, U₂=2/16, U₃=2/15).
- **Troyanskaya et al. (2002), Bioinformatics 18(5):679–688** (WebSearch / Oxford Academic abstract).
  Confirms the **summation** form: "summing up the values of combinatorial complexity over all
  values l from 1 to N … and dividing by maximal dictionary size" — i.e.
  `LC = Σ_i V_obs(i) / Σ_i V_max(i)`, computed via suffix trees in linear time.
- **Morgulis et al. (2006), J Comput Biol 13(5):1028–1040** (WebSearch; PubMed / sdust repo).
  Confirms the symmetric DUST scoring function `S = (1/(L−1)) · Σ_t c_t(c_t−1)/2`, where the
  c_t are triplet counts in the window and `L−1` is (number of triplets − 1). The function
  grows linearly with window length for windows longer than 4^k.
- **Shannon (1948) / Wikipedia "Entropy (information theory)"**: `H = −Σ p_i log₂ p_i`,
  max for a 4-letter alphabet = log₂4 = 2.

### Formula check
| Measure | Spec formula | Source | Status |
|---|---|---|---|
| Linguistic complexity | `Σ V_obs / Σ V_max`, `V_max(i)=min(4^i, N−i+1)` | Troyanskaya 2002 | ✅ matches |
| Shannon entropy | `H = −Σ p_i log₂ p_i` | Shannon 1948 | ✅ matches |
| k-mer entropy | Shannon H over k-mer frequencies | Shannon 1948 | ✅ matches |
| DUST score | `Σ c_t(c_t−1)/2 / (w−1)` | Morgulis 2006 | ✅ matches |

### Note (the only Stage-A caveat)
Two **distinct, both-published** "linguistic complexity" definitions exist: the Wikipedia
**product** form `C = ∏ U_i` and the Troyanskaya **summation** form `LC = Σobs/Σmax`. The
TestSpec explicitly selects and cites the **Troyanskaya summation** form, and the code
implements exactly that. This is internally consistent and correctly sourced; it is recorded
here only so the divergence from the Wikipedia *product* form is not mistaken for an error.
(The spec's "Wikipedia example" tests still reproduce because they apply the summation form to
Wikipedia's example sequences — the U_i component values agree; only the final aggregation
operator differs between the two definitions.)

### Edge-case semantics
Empty → 0; single nucleotide LC → 1.0 (V_max=min(4,1)=1, V_obs=1); homopolymer → minimal
LC and 0 entropy; maximally diverse (`ATGC`) → LC=1.0, H=2.0. All defined and sourced.

### Independent cross-check (hand/Python recomputation — all exact)
| Input | Measure | Recomputed | Spec/test expects |
|---|---|---|---|
| `ATGCTAGCATGCAATG` | LC (mw10) | obs91/max103 = 0.883495 | 91/103 ✅ |
| `AAAAAAAAAAAAAAAA` | LC (mw10) | obs10/max103 = 0.097087 | 10/103 ✅ |
| `A` | LC | 1/1 = 1.0 | 1.0 ✅ |
| `ATGC` | LC | 10/10 = 1.0 | 1.0 ✅ |
| `ACGGGAAGCTGATTCCA` | LC (mw4) | obs47/max49 = 0.959184 | 47/49 ✅ |
| `ACACACACACACACACA` | LC (mw10) | obs20/max112 = 0.178571 | 5/28 ✅ |
| `AAAAAAAAAAAAAAAAAA` | DUST | 120/15 = 8.0 | 8.0 ✅ |
| `ATGCTAGCATGCTAGC` | DUST | 6/13 = 0.461538 | 6/13 ✅ |
| `AAAAAAA` | DUST | 10/4 = 2.5 | 2.5 ✅ |
| `ATGCTAGC…AATGC` (N=30) | Compression | 112/216 = 0.518519 | 14/27 ✅ |
| 31×`A` | Compression | 10/224 = 0.044643 | 5/112 ✅ |
| `ATCG` | k-mer H (k2) | 1.5849625 | log₂3 ✅ |

**Stage A verdict: PASS-WITH-NOTES** — every formula matches its cited authoritative source
exactly; sole note is the documented product-vs-summation LC distinction, which the spec
selects and cites correctly.

---

## Stage B — Implementation

### Code path reviewed
`SequenceComplexity.cs`:
- LC: `CalculateLinguisticComplexityCore` (lines 39–66) — loops word length 1..min(maxWord,N),
  counts distinct substrings (HashSet), `maxPossible = min(4^wordLen, N−wordLen+1)`, returns
  `observedTotal/possibleTotal`. Matches Troyanskaya summation exactly.
- Shannon: lines 93–120 — frequencies over fixed `{A,T,G,C}` alphabet (non-ATGC ignored from
  numerator and denominator), `−Σ p log₂ p`. Matches.
- k-mer entropy: lines 136–161 — `H` over k-mer counts; `< k` returns 0. Matches.
- DUST: `CalculateDustScoreCore` (lines 311–336) — triplet counts, `Σ c(c−1)/2`, normalized by
  `total−1` (= w−1). Matches Morgulis. `len < wordSize` → 0; `total ≤ 1` → 0.
- Windowed / regions / masking / compression: read and traced; consistent with the above cores.

### Formula realised correctly?
Yes — exact, not approximate. The independent Python reimplementation above mirrors the C#
control flow line-for-line and produces identical numbers to both the code and the test
assertions.

### Cross-verification recomputed vs code
The 91 `~Complexity` tests assert the exact sourced values in the table above (`Within(1e-10)`
or exact equality for power-of-two-clean cases). All recomputed values match.

### Variant/delegate consistency
String overloads upper-case then call the same `*Core` as the `DnaSequence` overloads
(`StringOverload_Matches…` tests pass — bitwise identical). Windowed complexity / masking /
region finding reuse `CalculateShannonEntropyCore` and `CalculateDustScoreCore`, so variants
are consistent by construction.

### Numerical robustness
Div-by-zero guarded (`possibleTotal>0`, `total>1`, `total==0` early returns). `count*(count-1)/2.0`
uses double; `Math.Pow(4,wordLen)` cast to long is safe within the bounded word lengths used.
No overflow on stated ranges.

### Test quality audit
Assertions check exact externally-sourced values (Troyanskaya/Wikipedia/Morgulis/Shannon),
not tautologies; guard clauses (null, k<1, maxWord<1, windowSize<1, stepSize<1) tested; range
invariants (LC∈[0,1], H∈[0,2], k-mer H∈[0,log₂4^k], compression∈[0,1]) tested; deterministic.

### Findings / defects
None. No code change required.

---

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** (LC product-vs-summation distinction documented; spec selects the
  cited Troyanskaya summation form correctly).
- **Stage B: PASS** — implementation faithfully realises every validated formula; all worked
  examples reproduce exactly.
- **End-state: CLEAN** — no defect found; no code changed.
- **Tests:** `--filter FullyQualifiedName~Complexity` → 91 passed, 0 failed. Full suite:
  4461 passed, 0 failed.
