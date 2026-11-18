# Tandem Repeat Detection (Genomic Analysis)

| Field | Value |
|-------|-------|
| Algorithm Group | Genomic Analysis |
| Test Unit ID | GENOMIC-TANDEM-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production (exact tandem repeats = maximal runs with a primitive unit; approximate repeats are `RepeatFinder.FindApproximateTandemRepeats`, TRF model) |
| Last Reviewed | 2026-11-18 (B09 review campaign, F35–F38) |

## 1. Overview

Tandem repeat detection identifies contiguous DNA segments where a repeat unit (pattern) occurs
consecutively two or more times, e.g. `ATTCG ATTCG ATTCG` for the unit `ATTCG` repeated three
times [1][2]. `GenomicAnalyzer.FindTandemRepeats` is the **exact** detector: for every unit length it
reports each maximal run of a primitive unit once, with the unit, the 0-based start, the number of
complete copies, the total length and the full repeated substring [3][4][5].

> **Relation to REP-TANDEM-001 / microsatellites.** `FindTandemRepeats` is a facade over the
> canonical `RepeatFinder.FindMicrosatellites` (REP-STR-001: one linear scan per unit length; same
> maximal-run + primitive-unit convention as MISA and pytrf/Krait), with an additional run-enumeration
> engine for wide unit-length ranges (§4). The only differences to `FindMicrosatellites` are the
> parameter defaults (unit length ≥ 2, ≥ 2 copies, no upper bound on the unit length) and the result
> type (`TandemRepeat`). See also [Repeat_Analysis/Tandem_Repeat_Detection.md](../Repeat_Analysis/Tandem_Repeat_Detection.md).
>
> **Approximate (interrupted / noisy) tandem repeats are implemented in-repo, not delegated to external
> tools:** `RepeatFinder.FindApproximateTandemRepeats` (Tandem Repeats Finder model of Benson 1999,
> parity with TRF 4.10.0; REP-APPROX-001) and the MCP tools `find_approximate_tandem_repeats` and
> `mask_approximate_tandem_repeats`. `GenomicAnalyzer` deliberately does not wrap them: the exact detector
> and the TRF detector have different result types (copies/consensus/percent-matches/score vs. unit/copies)
> and parameter sets, and the approximate detector is a self-contained canonical implementation in
> `RepeatFinder` (see §5.3).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A tandem repeat is the product of a tandem-duplication event: a pattern of nucleotides repeated in
directly adjacent copies [1][2]. Tandem repeats constitute about 8% of the human genome and are
implicated in more than 50 lethal human diseases [2]. By unit length they are classified as
microsatellites/STRs (short units), minisatellites/VNTRs (10–60 nt) and macrosatellites (≈1,000
nt) [2].

### 2.2 Core Model

Per Benson (1999), "a tandem repeat in DNA is two or more contiguous … copies of a pattern of
nucleotides" [1] (Benson's copies are *approximate*; this detector accepts exact copies only). For a
sequence $S$ of length $n$ and a unit length (period) $p$:

* An interval $S[a..e)$ has period $p$ when $S[x] = S[x+p]$ for all $a \le x < e-p$.
* It is a **maximal repetition (run)** when it cannot be extended to the left or right with the same
  period: $a = 0$ or $S[a-1] \ne S[a-1+p]$, and $e = n$ or $S[e] \ne S[e-p]$ (Kolpakov & Kucherov 1999 [3];
  Bannai et al. 2017, Definition 1 [4]).
* Its **copy number** is $\lfloor (e-a)/p \rfloor$ (complete copies; a trailing partial copy is not
  counted) and it is reported once, at its left end $a$, with unit $U = S[a..a+p)$.
* The unit must be **primitive** (not a power $w^m$, $m \ge 2$, of a shorter word): `ATATATAT` is `AT`×4,
  not additionally `ATAT`×2. This is the microsatellite convention — MISA's `misa.pl` rejects "false type
  motifs [e.g. (TT)6 or (ACAC)5]" [5] — and it is exactly the Kolpakov–Kucherov notion of a run: if the
  unit of a period-$p$ interval of length $\ge 2p$ is primitive, $p$ is the *smallest* period of the
  interval (Fine–Wilf), and conversely a run with smallest period $p$ has a primitive unit of length $p$.

The repeat is reported when $p \in [\text{minUnitLength}, \text{maxUnitLength}]$ and
$\lfloor (e-a)/p \rfloor \ge \text{minRepetitions}$ ($\ge 2$). Reported set (the *runs* of $S$ with
admissible period and copy number):

$$
\{(a, p, k) : S[a..a+kp) = U^k,\; U \text{ primitive},\; k = \lfloor (e-a)/p\rfloor \ge \text{minRepetitions},\; (a,e,p)\text{ a run}\}
$$

A string of length $n$ has fewer than $n$ runs (the "runs theorem", Bannai et al. 2017 [4]), so the
output is always $O(n)$ long (for `minUnitLength = 1`).

### 2.3 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported repeat has Repetitions ≥ `minRepetitions` ≥ 2 (k ≥ 2) and $\text{minUnitLength} \le |U| \le \text{maxUnitLength}$. | a tandem repeat requires two or more contiguous copies [1]. |
| INV-02 | `TotalLength = Unit.Length × Repetitions`. | total length = period × copy number [1]. |
| INV-03 | `Position + Unit.Length × Repetitions ≤ sequence.Length`. | only complete copies inside the run are counted. |
| INV-04 | The unit is primitive; the run is maximal (not extendable left/right with period $|U|$). | Kolpakov–Kucherov run definition [3][4]; MISA primitive motifs [5]. |
| INV-05 | Output order: unit length ascending, then start ascending; each (unit length, start) at most once. | per-unit-length left-to-right scan; the run enumeration is sorted by (period, start). |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` | required | DNA sequence to scan (upper-case A/C/G/T by construction; case is irrelevant). | `null` → `ArgumentNullException`. |
| `minUnitLength` | `int` | `2` | Minimum repeat-unit length (period). | ≥ 1, else `ArgumentOutOfRangeException` (eager). |
| `minRepetitions` | `int` | `2` | Minimum number of complete consecutive copies (k). | ≥ 2 per the definition [1], else `ArgumentOutOfRangeException` (eager). |
| `maxUnitLength` | `int` | `int.MaxValue` (3-argument overload) | Maximum unit length, e.g. 6 for microsatellites, 60 for minisatellites. | ≥ `minUnitLength`, else `ArgumentOutOfRangeException` (eager). A unit longer than `n / minRepetitions` can never repeat, so "no bound" costs nothing. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Tandem repeats | `IEnumerable<TandemRepeat>` | Each hit: `Unit` (primitive, the motif phase at the run's left end), `Position` (0-based start), `Repetitions` (complete copies), `TotalLength`, `FullSequence`. Order: unit length ascending, then start. |

### 3.3 Preconditions and Validation

All three parameter checks and the null check run eagerly (at the call, not on enumeration). An empty
sequence, or a sequence shorter than `minUnitLength × minRepetitions`, yields no hits. Positions are
0-based. The MCP tool `find_tandem_repeats` exposes the same three parameters (`maxUnitLength`
optional) and delegates to this method.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate; let `K = min(maxUnitLength, n / minRepetitions) − minUnitLength + 1` be the number of admissible unit lengths.
2. **`K ≤ 256` — per-unit-length scan** (`RepeatFinder.FindMicrosatellites(sequence, minUnitLength, maxUnitLength, minRepetitions)`): for each unit length `p` visit only left-maximal positions (`S[i−1] ≠ S[i−1+p]`), extend to the right-maximal end, report the run once when the unit `S[i..i+p)` is primitive and `⌊len/p⌋ ≥ minRepetitions`. O(n) comparisons per unit length.
3. **`K > 256` — run enumeration** (`MaximalRepetitions.FindRuns`, internal): all runs of the string via the Lyndon-array algorithm of the runs theorem [4] (§4.2), then keep those with `p` in range and `⌊len/p⌋ ≥ minRepetitions`, sorted by (period, start). The result is identical to step 2 (same set, same order); the cost does not depend on `K`.

### 4.2 Run enumeration (Bannai, I, Inenaga, Nakashima, Takeda, Tsuruta 2017) [4]

For each of the two alphabet orders $\prec_0$ (natural) and $\prec_1$ (reversed), reading $S$ as $S\$$ with $\$$ smaller than every letter under $\prec_0$ (hence larger under $\prec_1$): $\lambda(i)$ = the longest Lyndon word starting at $i$ = $S[i..j)$ with $j$ the next position whose suffix is smaller (one monotonic-stack pass; suffix comparisons by longest-common-extension queries). Every run has, for one of the two orders, an L-root $S[i..i+p)$ equal to $\lambda(i)$ with $i > a$ (Lemma 7 [4]); for each candidate $(i, p = |\lambda(i)|)$ the period is extended right with $\text{LCE}(i, i+p)$ and left with the longest common suffix of the prefixes ending at $i-1$ and $i+p-1$, and the candidate is a run when it spans $\ge 2p$ (§4.2 [4]; a Lyndon word is primitive, so the smallest period of the interval is $p$). LCE queries are O(1) from the suffix array (prefix doubling, `SequenceComplexity.BuildSuffixArray`), the Kasai LCP array and a block sparse-table RMQ of the text and of its reverse. Duplicates (a run found from both orders / several roots) are removed.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| per-unit-length scan | `O(n · K)` | `O(1)` plus output | ≈ 17 ms per unit length at 1 Mbp (Debug build); 1 Mbp random, unit lengths 2–60, ≥ 5 copies: ≈ 1 s. |
| run enumeration | `O(n log² n)` (suffix array by prefix doubling; the scan itself is O(n) LCE queries) | `O(n)` | 1 Mbp random, unbounded unit length: ≈ 3 s (Debug); 1 Mbp homopolymer / period-8 string (worst case for the doubling): ≈ 6 s. Break-even ≈ 175 unit lengths at 1 Mbp; the switch is set to 256. |
| former substring scan (removed) | `O(n² · m)` | `O(1)` | hours on 1 Mbp. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs), [MaximalRepetitions.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MaximalRepetitions.cs), canonical scan in [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs).

- `GenomicAnalyzer.FindTandemRepeats(DnaSequence, int minUnitLength = 2, int minRepetitions = 2)`: exact detector, no upper unit length.
- `GenomicAnalyzer.FindTandemRepeats(DnaSequence, int, int, int maxUnitLength)`: restricted unit-length range.
- MCP: `find_tandem_repeats` (Analysis server), parameters `sequence`, `minUnitLength`, `minRepetitions`, optional `maxUnitLength`.

### 5.2 Current Behavior

The facade validates the arguments, then either runs the canonical per-unit-length scan or the run
enumeration (§4.1). Before the 2026-11 review the method was a separate O(n²·m) substring scan that
(a) reported non-primitive units alongside their primitive root (`ATATATAT`, minUnit 1, minReps 2 →
`AT`×4@0 **and** `ATAT`×2@0), (b) was cubic in practice and (c) had no upper bound on the unit length.
On the 1156-row oracle (§6.3) the old scan differed from the definition on 276 rows, in every case by
extra non-primitive units (487 extra rows; on 20 000 further random strings 26 747 extras, never a
missing or non-maximal primitive hit). The old doc/DUP_MAP remark "not left-maximal" does not hold for
the primitive hits: they were already maximal runs.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact detection of two-or-more contiguous copies of a unit (k ≥ 2) [1][2], reported as maximal repetitions with a primitive unit [3][4][5].
- Unit, 0-based start, complete-copy count, total length, deterministic order.
- Unit-length bounds `[minUnitLength, maxUnitLength]` and the copy-number threshold.

**Not part of this method (implemented elsewhere in the repository):**

- Approximate / interrupted tandem repeats with alignment scoring: `RepeatFinder.FindApproximateTandemRepeats` (TRF 4.10.0 model, 99.8 % row parity on the B04 oracle; REP-APPROX-001) and MCP `find_approximate_tandem_repeats` / `mask_approximate_tandem_repeats`. No pass-through is added to `GenomicAnalyzer`: it would duplicate a canonical public API with a different result type.
- Microsatellite statistics and MISA/Krait classification: `RepeatFinder.GetTandemRepeatSummary` and `FindMicrosatellites` overloads (REP-STR-001 / REP-TANDEM-001).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Exact copies only (Benson's TRF copies are approximate). | Scope | Repeats with substitutions/indels are not reported here. | by design | Use `FindApproximateTandemRepeats`. |
| 2 | The reported unit is the motif phase at the run's left end (e.g. `GCAGCAGCAGCAGG` → `GCA`×4@0, not `CAG`×4@1). | Convention | Rotation is not canonicalised. | by design | MISA / pytrf convention [5]; rotations of one run are never reported twice. |
| 3 | Unit lengths ≤ 256 distinct values use the per-unit-length scan, wider ranges the run enumeration. | Implementation | Same results (differential-tested), different cost. | by design | Threshold = measured break-even. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence / shorter than `minUnitLength × minRepetitions` | no results | `minRepetitions` copies cannot fit. |
| No tandem present | empty enumerable | no run reaches the copy threshold. |
| Entire sequence is one tandem | one result spanning the repeated region | the run is the whole sequence (`⌊n/p⌋` copies). |
| Non-primitive multiples (`ATATATAT`, unit 4) | not reported (`AT`×4 only) | primitive-unit rule [5]. |
| Partial last copy (`ATATATA`) | `AT`×3@0 | only complete copies are counted. |
| Rotations (`TA` inside `ATATATA`) | not reported again | one report per maximal run. |
| Runs of different unit lengths over overlapping regions (`AAAGAAAGAAAG`, minUnit 1) | each reported (`A`×3@0/4/8, `AAAG`×3@0) | different periods are different runs. |
| `minUnitLength > maxUnitLength`, `minUnitLength < 1`, `minRepetitions < 2`, `null` | `ArgumentOutOfRangeException` / `ArgumentNullException`, eagerly | contract §3.3. |
| lower case input | `DnaSequence` normalises to upper case | — |

### 6.2 Limitations

Exact detector only (see §5.3 for the approximate one). There is no raw-string overload (a `DnaSequence`
is A/C/G/T only; `RepeatFinder.FindMicrosatellites(string, …)` accepts N/IUPAC text). Results are
materialised eagerly in the run-enumeration regime (O(n) memory).

### 6.3 Oracles

* **Python, three independent methods asserted equal** (`scripts/analysis/generate_tandem_oracle.py`): (A) brute-force Kolpakov–Kucherov runs (every interval, its smallest period by definition, maximality by extension), (B) boolean-array maximal runs of `S[x] = S[x+p]` + primitivity test, (C) regex `(?=((.{p})\2{k-1,}))` at left-maximal starts + primitivity; 1156 `(sequence, minUnit, minReps, maxUnit)` rows → `TestData/Analysis/tandem_oracle.tsv`. Both C# engines match all rows.
* **In-test brute force** (C#, every interval's smallest period): 1500 random strings × bounded/unbounded, plus 62 000 exhaustive strings (all `{A,C}` ≤ 14, all `{A,C,G}` ≤ 9) for the run enumeration, plus 60 long (2.5 kb) low-entropy strings engine-vs-engine.
* EMBOSS `etandem` 6.6.0 was run as an informational third reference: it scores diagonals with its own threshold model (score = matches, `-threshold` global across sizes), reports consensus rotations differently and is not an exact-run enumerator (agreement 297/400 random cases at fixed `-minrepeat = -maxrepeat`, differences all in rotation phase / sub-threshold diagonals), so it is not used as an oracle.

## 7. Examples and Related Material

### 7.1 Worked Examples

```csharp
var dna = new DnaSequence("ATTCGATTCGATTCG");
var hits = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 5, minRepetitions: 2).ToList();
// hits[0].Unit == "ATTCG", hits[0].Position == 0, hits[0].Repetitions == 3 (total length 15)

var at = GenomicAnalyzer.FindTandemRepeats(new DnaSequence("ATATATAT"), 1, 2).ToList();
// AT x4 @0 only (not also ATAT x2)

var micro = GenomicAnalyzer.FindTandemRepeats(genome, minUnitLength: 1, minRepetitions: 3, maxUnitLength: 6);
// microsatellite-sized units only
```

The first is the Wikipedia worked example: unit `ATTCG` (period 5), copy number 3 [2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomicAnalyzer_TandemRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_TandemRepeat_Tests.cs) (sourced worked examples, `INV-01`–`INV-03`), [GenomicAnalyzer_TandemRepeat_Oracle_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_TandemRepeat_Oracle_Tests.cs) (oracles, both engines, 1 Mbp), fuzz [GenomicTandemFuzzTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Fuzzing/GenomicTandemFuzzTests.cs)
- Evidence: [GENOMIC-TANDEM-001-Evidence.md](../../../docs/Evidence/GENOMIC-TANDEM-001-Evidence.md); report [GENOMIC-TANDEM-001.md](../../Validation/reports/GENOMIC-TANDEM-001.md)
- Related algorithms: [Repeat_Analysis/Tandem_Repeat_Detection.md](../Repeat_Analysis/Tandem_Repeat_Detection.md) (REP-TANDEM-001, summary helpers), [Repeat_Analysis/Microsatellite_Detection.md](../Repeat_Analysis/Microsatellite_Detection.md) (the canonical scan), approximate repeats: `RepeatFinder.FindApproximateTandemRepeats` (REP-APPROX-001)

## 8. References

1. Benson, G. 1999. Tandem repeats finder: a program to analyze DNA sequences. *Nucleic Acids Research* 27(2):573–580. https://doi.org/10.1093/nar/27.2.573
2. Wikipedia contributors. Tandem repeat. https://en.wikipedia.org/wiki/Tandem_repeat (accessed 2026-06-14)
3. Kolpakov, R., Kucherov, G. 1999. Finding maximal repetitions in a word in linear time. *Proc. 40th IEEE FOCS*, 596–604.
4. Bannai, H., I, T., Inenaga, S., Nakashima, Y., Takeda, M., Tsuruta, K. 2017. The "Runs" Theorem. *SIAM J. Comput.* 46(5):1501–1514. https://arxiv.org/abs/1406.0263 (opened 2026-11-18: Definition 1, Lemmas 3, 6, 7, 8, §3, §4.1–4.2)
5. Thiel, T., Michalek, W., Varshney, R.K., Graner, A. 2003. Exploiting EST databases for the development and characterization of gene-derived SSR-markers in barley. *Theor. Appl. Genet.* 106:411–422; `misa.pl` (mirror https://raw.githubusercontent.com/cfljam/SSR_marker_design/master/misa.pl, opened 2026-11-18: line 110 "reject false type motifs [e.g. (TT)6 or (ACAC)5]").
