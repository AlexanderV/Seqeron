# Tandem Repeat Detection (Genomic Analysis)

| Field | Value |
|-------|-------|
| Algorithm Group | Genomic Analysis |
| Test Unit ID | GENOMIC-TANDEM-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Tandem repeat detection identifies contiguous DNA segments where a repeat unit (pattern) occurs
consecutively two or more times, e.g. `ATTCG ATTCG ATTCG` for the unit `ATTCG` repeated three
times [1][2]. `GenomicAnalyzer.FindTandemRepeats` is an exact detector: it reports the unit, the
0-based start, the repetition count, the total length and the full repeated substring for each
exact tandem block. The computation is exact (not approximate) and uses direct string comparison
over candidate unit lengths and start positions [1].

> **Consolidation note.** This unit's canonical method is identical to the one already delivered
> under Test Unit **REP-TANDEM-001** (`GenomicAnalyzer.FindTandemRepeats`). GENOMIC-TANDEM-001 is a
> duplicate Registry entry (same method, same class, same O(n²) complexity); it is resolved by
> consolidation, not re-implementation. See the REP-TANDEM-001 algorithm doc
> ([Tandem_Repeat_Detection.md](../Repeat_Analysis/Tandem_Repeat_Detection.md)) for the full
> treatment, including the `RepeatFinder.GetTandemRepeatSummary` aggregation helper.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A tandem repeat is the product of a tandem-duplication event: a pattern of nucleotides repeated in
directly adjacent copies [1][2]. Tandem repeats constitute about 8% of the human genome and are
implicated in more than 50 lethal human diseases [2]. By unit length they are classified as
microsatellites/STRs (short units), minisatellites/VNTRs (10–60 nt) and macrosatellites (≈1,000
nt) [2].

### 2.2 Core Model

Per Benson (1999), "a tandem repeat in DNA is two or more contiguous … copies of a pattern of
nucleotides" [1]. For a sequence $S$, start $p$, unit (pattern) $U$ of length (period) $|U|$, and
copy number $k$, a tandem repeat is present when:

$$
S[p..p + k|U|) = U^k \quad \text{with} \quad k \ge 2
$$

The **period** is $|U|$ (for an exact repeat) and the **copy number** is $k$ [1]. Wikipedia gives
the worked example `ATTCG ATTCG ATTCG`: unit `ATTCG`, period 5, copy number 3 [2].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported repeat has Repetitions ≥ `minRepetitions` ≥ 2 (k ≥ 2). | a tandem repeat requires two or more contiguous copies [1]. |
| INV-02 | `TotalLength = Unit.Length × Repetitions`. | total length = period × copy number [1]. |
| INV-03 | `Position + Unit.Length × Repetitions ≤ sequence.Length`. | copies are contiguous and counted within sequence bounds [1][2]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` | required | DNA sequence to scan. | dereferenced directly; caller supplies non-null. |
| `minUnitLength` | `int` | `2` | Minimum candidate repeat-unit (period) length. | positive; not explicitly guarded. |
| `minRepetitions` | `int` | `2` | Minimum number of contiguous copies (k). | ≥ 2 per the definition [1]; not explicitly guarded. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Tandem repeats | `IEnumerable<TandemRepeat>` | Each hit: unit, 0-based start, repetition count, total length, full repeated substring. |

### 3.3 Preconditions and Validation

`FindTandemRepeats` reads `sequence.Sequence` directly and does not validate threshold values; it
relies on the caller to provide a non-null `DnaSequence` and sensible thresholds. An empty
sequence yields no hits. Positions are 0-based. Matching is exact (case is normalized by
`DnaSequence` on construction).

## 4. Algorithm

### 4.1 High-Level Steps

1. For each unit length from `minUnitLength` to `sequence.Length / minRepetitions`, choose a candidate period.
2. For each start where at least `minRepetitions` copies could fit, extract the candidate unit.
3. Count consecutive exact copies by advancing in unit-length increments until the pattern breaks.
4. If the count meets the threshold, yield a `TandemRepeat` and advance the scan to the end of that block within the current unit-length pass.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `GenomicAnalyzer.FindTandemRepeats` | `O(n² × m)` | `O(1)` plus yielded output | `m` = effective max unit length explored; brute-force substring comparison. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `GenomicAnalyzer.FindTandemRepeats(DnaSequence, int, int)`: exact detector for consecutive tandem repeats.

### 5.2 Current Behavior

`FindTandemRepeats` uses a brute-force scan over candidate unit lengths and start positions,
compares units with direct substring equality, and skips forward after each hit within the current
unit-length pass. **Suffix-tree decision:** although Wikipedia notes tandem repeats can be detected
with suffix trees/arrays [2], the repository suffix tree (`SuffixTree.FindAllOccurrences`) is
designed for enumerating occurrences of a *given* pattern, not for discovering unknown repeat
units; using it here would not change the formally required output and would not simplify the
unit-length/period enumeration this detector performs. The detector therefore uses the
brute-force scan (the same code path validated under REP-TANDEM-001); the suffix tree was evaluated
and not adopted.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact detection of two-or-more contiguous copies of a unit (k ≥ 2) [1][2].
- Reporting of unit (pattern), 0-based start, copy number (repetitions) and total repeated length [1].

**Intentionally simplified:**

- Exact-copy detection only, instead of Benson's *approximate* tandem copies [1]; **consequence:**
  repeats with internal mismatches/indels are not reported, but exact tandems match the definition.
- Brute-force comparison instead of suffix-tree/suffix-array optimization; **consequence:** runtime
  grows on long sequences; best suited to moderate lengths [1][2].

**Not implemented:**

- Approximate / interrupted tandem-repeat scoring; **users should rely on:** external tools
  (e.g. Tandem Repeats Finder [1]) for noisy repeat families.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `FindTandemRepeats` assumes valid threshold inputs rather than validating them. | Assumption | Nonsensical thresholds can yield undefined/exception behavior. | accepted | Definition requires k ≥ 2; floor not enforced in code. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | no results | the scan loops do not execute when no full block fits. |
| No tandem present | empty enumerable | no candidate unit reaches the copy threshold. |
| Entire sequence is one tandem | one result spanning the repeated region | the counting loop runs until the pattern breaks. |
| Same region, multiple periods | reported once per unit-length interpretation meeting the threshold | the definition is satisfied for each valid period [1][2]. |

### 6.2 Limitations

The detector is exact and does not score approximate, interrupted, or noisy tandem families, and
does not canonicalize across competing period interpretations.

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
var dna = new DnaSequence("ATTCGATTCGATTCG");
var hits = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 5, minRepetitions: 2).ToList();
// hits[0].Unit == "ATTCG", hits[0].Position == 0, hits[0].Repetitions == 3 (total length 15)
```

This is the Wikipedia worked example: unit `ATTCG` (period 5), copy number 3 [2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomicAnalyzer_TandemRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_TandemRepeat_Tests.cs) — covers `INV-01`–`INV-03` (shared canonical fixture with REP-TANDEM-001)
- Evidence: [GENOMIC-TANDEM-001-Evidence.md](../../../docs/Evidence/GENOMIC-TANDEM-001-Evidence.md)
- Related algorithms: [Tandem_Repeat_Detection.md](../Repeat_Analysis/Tandem_Repeat_Detection.md) (REP-TANDEM-001 — same method)

## 8. References

1. Benson, G. 1999. Tandem repeats finder: a program to analyze DNA sequences. *Nucleic Acids Research* 27(2):573–580. https://doi.org/10.1093/nar/27.2.573
2. Wikipedia contributors. Tandem repeat. https://en.wikipedia.org/wiki/Tandem_repeat (accessed 2026-06-14)
