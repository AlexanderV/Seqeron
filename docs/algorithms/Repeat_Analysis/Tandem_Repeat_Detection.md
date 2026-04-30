# Tandem Repeat Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-TANDEM-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Tandem repeat detection identifies contiguous DNA segments where a repeat unit occurs consecutively two or more times, such as `ATTCGATTCGATTCG` for the unit `ATTCG` repeated three times [1]. The repository exposes an exact tandem detector in `GenomicAnalyzer.FindTandemRepeats` and a separate aggregation helper, `RepeatFinder.GetTandemRepeatSummary`, that summarizes only microsatellite-sized tandem repeats with unit lengths from 1 to 6 bp [1][2]. The canonical detector uses direct string comparison and skips to the end of each detected tandem block within the current unit-length pass, though the same region can still be reported under a different unit-length interpretation. Tandem repeats are biologically important because they occupy a substantial fraction of the human genome and include disease-associated repeat expansions and microsatellite-instability markers [1][2].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A tandem repeat is a sequence of the form $U^k$, where a unit $U$ is repeated consecutively $k$ times with no gap between copies [1]. Common terminology and size classes are:

| Term | Definition | Typical Unit Length |
|------|------------|---------------------|
| Microsatellite (STR) | Short tandem repeat | 1-6 bp, sometimes extended to 1-10 bp [2][3] |
| Minisatellite (VNTR) | Variable-number tandem repeat | 10-60 bp [1] |
| Macrosatellite | Large tandem repeat array | approximately 1,000+ bp [1] |

The repository summary logic uses the standard mono-, di-, tri-, tetra-, penta-, and hexanucleotide categories [2]. Tandem repeats are implicated in trinucleotide-repeat disorders, microsatellite instability in cancer, and broader genome variability; the legacy reference set also notes that tandem repeats account for roughly 8% of the human genome and are linked to more than 50 human diseases [1][2]. Their dominant mutation mechanism is replication slippage, with microsatellite mutation rates reported around one slippage event per 1,000 generations [2][3].

### 2.2 Core Model

For a sequence $S$, position $p$, repeat unit $U$, and repetition count $k$, a tandem repeat is present when:

$$
S[p..p + k|U|) = U^k \quad \text{with} \quad k \ge 2
$$

The canonical detector searches candidate unit lengths and starting positions, counts consecutive copies of each candidate unit, and emits a result when the repetition count reaches the configured threshold. The summary helper applies the same idea indirectly through microsatellite detection and aggregates counts, bases covered, repeat-type totals, longest repeat, and most frequent repeat unit.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported tandem repeat contains at least `minRepetitions` consecutive copies of a unit whose length is at least `minUnitLength`. | The detector yields only after counting consecutive equal substrings. |
| INV-02 | `TotalLength = Unit.Length × Repetitions` for every `TandemRepeat`. | Total length is defined by the unit size and repetition count. |
| INV-03 | `Position + Unit.Length × Repetitions <= sequence.Length`. | The counting loop stops when the next full unit would exceed sequence bounds. |
| INV-04 | The summary percentages and totals are derived only from reported microsatellites. | `GetTandemRepeatSummary` delegates to `FindMicrosatellites(sequence, 1, 6, minRepeats)`. |
| INV-05 | Within a fixed candidate unit length, later starts inside a detected tandem block are skipped. | After yielding a result, the implementation advances the start index to the end of the detected tandem block for that unit-length pass. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` | required | DNA sequence to analyze. | `FindTandemRepeats` dereferences `sequence.Sequence` directly; `GetTandemRepeatSummary` throws on `null`. |
| `minUnitLength` | `int` | `2` | Minimum candidate repeat-unit length for `FindTandemRepeats`. | The algorithm assumes a positive value; no explicit guard is implemented in `GenomicAnalyzer.FindTandemRepeats`. |
| `minRepetitions` | `int` | `2` | Minimum number of consecutive unit copies for `FindTandemRepeats`. | The algorithm assumes at least two repetitions; no explicit guard is implemented in `GenomicAnalyzer.FindTandemRepeats`. |
| `minRepeats` | `int` | `3` | Minimum repeat count used by `GetTandemRepeatSummary`. | Passed directly to `FindMicrosatellites(sequence, 1, 6, minRepeats)`, which rejects values below 2. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Tandem repeats | `IEnumerable<TandemRepeat>` | Exact tandem-repeat hits with unit, 0-based start position, repetition count, total length, and full repeated sequence. |
| Tandem summary | `TandemRepeatSummary` | Aggregate summary over microsatellite-sized tandem repeats, including total repeat count, total repeat bases, percentage of sequence, longest repeat, most frequent repeat unit, and dedicated per-class counts for mono-, di-, tri-, and tetranucleotide repeats. |

### 3.3 Preconditions and Validation

`GetTandemRepeatSummary` throws `ArgumentNullException` when `sequence` is `null` because it delegates to `FindMicrosatellites`. `FindTandemRepeats` does not perform explicit argument validation; it reads `sequence.Sequence` immediately and therefore relies on the caller to provide a non-null `DnaSequence` and sensible threshold values. Empty sequences produce no tandem-repeat hits, and an empty sequence summarized through `GetTandemRepeatSummary` returns zero totals and `0` percent coverage.

## 4. Algorithm

### 4.1 High-Level Steps

1. For each unit length from `minUnitLength` to `sequence.Length / minRepetitions`, choose a candidate repeat size.
2. For each start position where at least `minRepetitions` copies could fit, extract the candidate unit.
3. Count consecutive occurrences of that unit by advancing in `unitLength` increments until the pattern breaks.
4. If the repetition count meets the threshold, yield a `TandemRepeat` and advance the scan to the end of that tandem block within the current unit-length pass.
5. For summary mode, detect microsatellites with unit sizes 1-6 and aggregate counts, bases covered, repeat-type totals, longest repeat, and most frequent unit.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `GenomicAnalyzer.FindTandemRepeats` | `O(n^2 × m)` | `O(1)` plus yielded output | `m` is the effective maximum unit length explored by the nested loops and substring comparisons. |
| `RepeatFinder.GetTandemRepeatSummary` | `O(n × U × R)` | `O(k)` | Delegates to microsatellite detection, where `U` is the searched unit-length range, `R` is the average repeat count, and `k` is the number of microsatellites retained. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation locations:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs), [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `GenomicAnalyzer.FindTandemRepeats(DnaSequence, int, int)`: Canonical exact detector for consecutive tandem repeats.
- `RepeatFinder.GetTandemRepeatSummary(DnaSequence, int)`: Summary helper that aggregates microsatellite-sized tandem repeats.

### 5.2 Current Behavior

`GenomicAnalyzer.FindTandemRepeats` uses a brute-force scan over candidate unit lengths and positions, compares units with direct substring equality, and skips forward after each hit within the current unit-length pass. This suppresses later starts inside the same detected block for that unit length, but it does not prevent the same region from being reported again under a different unit-length interpretation. It does not normalize case or validate parameters before iterating. `RepeatFinder.GetTandemRepeatSummary` does validate `sequence`, then delegates to `FindMicrosatellites(sequence, 1, 6, minRepeats)`, meaning the summary covers only tandem repeats with 1-6 bp units and inherits microsatellite overlap suppression and redundant-unit filtering from that implementation. The summary record tracks total tandem-repeat counts across that full 1-6 bp range, but dedicated per-class fields stop at tetranucleotide repeats.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact detection of directly adjacent repeated units in DNA sequences [1][2].
- Standard short tandem repeat classification by unit length from mononucleotide through hexanucleotide classes [2][3].
- Reporting of tandem-repeat start position, unit, repetition count, and total repeated length for each exact hit.

**Intentionally simplified:**

- Brute-force direct substring comparison instead of suffix-tree, suffix-array, or Tandem Repeats Finder style optimization; **consequence:** runtime grows rapidly on long sequences and the implementation is best suited to moderate sequence lengths [1][4].
- `GetTandemRepeatSummary` is restricted to microsatellite-sized units from 1 to 6 bp; **consequence:** longer minisatellite and macrosatellite tandems are excluded from the summary even though `FindTandemRepeats` can detect longer exact units.
- `GetTandemRepeatSummary` aggregates 1-6 bp microsatellites into totals, but its dedicated count fields stop at tetranucleotide repeats; **consequence:** penta- and hexanucleotide repeats contribute to total counts and bases without receiving their own named output fields.

**Not implemented:**

- Optimized large-genome tandem-repeat indexing or approximate repeat scoring; **users should rely on:** `GenomicAnalyzer.FindTandemRepeats` for exact tandem blocks and `RepeatFinder.GetTandemRepeatSummary` only for microsatellite aggregation.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `FindTandemRepeats` assumes valid threshold inputs instead of validating them explicitly. | Assumption | Callers can trigger undefined or exception-driven behavior with nonsensical values such as very small unit lengths or repetition counts. | accepted | The legacy doc states that `minRepetitions` should be at least 2, but that floor is not enforced in code. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | `FindTandemRepeats` yields no results; summary returns zero totals. | The scan loops do not execute when no full repeat block can fit. |
| No repeats found | Empty tandem-repeat enumerable or zero-count summary. | No candidate unit reaches the repetition threshold. |
| Entire sequence is one tandem block | One result spanning the full repeated region. | The counting loop continues until the pattern breaks or the sequence ends. |
| Overlapping tandem patterns | Later overlapping starts inside a detected block are skipped within the current unit-length pass, but the same region can still be reported under a different unit length. | The canonical detector advances `start` to the end of the detected tandem only inside the active unit-length loop. |
| `unitLength > sequence.Length / minRepetitions` | No result for that unit size. | The outer loop bounds prevent impossible repeat sizes from being checked. |

### 6.2 Limitations

The canonical detector is exact and does not score approximate tandem repeats, interrupted repeats, or noisy repeat families. The summary helper is narrower than the canonical detector because it only considers 1-6 bp units. The detector also does not canonicalize across competing unit-length interpretations, so the same genomic region can appear more than once when different repeat-unit sizes satisfy the threshold. There is also no raw-string overload for `FindTandemRepeats`, so callers must provide a `DnaSequence` and handle any normalization before calling the algorithm.

## 7. Examples and Related Material

### 7.2 Related Use Cases

- Forensic DNA profiling: STR markers are standard forensic markers, with tetra- and pentanucleotide loci commonly preferred for robust genotyping [2].
- Paternity and kinship analysis: High STR polymorphism makes tandem repeats useful for relationship inference [2].
- Population genetics: Tandem-repeat variability supports diversity and lineage studies [2][3].
- Cancer diagnostics: Microsatellite instability is a repeat-based marker used in oncology workflows [1][2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomicAnalyzer_TandemRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_TandemRepeat_Tests.cs)
- Test spec: [REP-TANDEM-001.md](../../../tests/TestSpecs/REP-TANDEM-001.md)
- Related property tests: [RepeatFinderProperties.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Properties/RepeatFinderProperties.cs)
- Related metamorphic tests: [MetamorphicTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/MetamorphicTests.cs)

### 7.4 Change History

| Date | Version | Author | Changes |
|------|---------|--------|---------|
| 2026-01-22 | 1.0 | Algorithm QA | Initial documentation |

## 8. References

1. Wikipedia. 2026. Tandem repeat. Wikipedia. https://en.wikipedia.org/wiki/Tandem_repeat
2. Wikipedia. 2026. Microsatellite. Wikipedia. https://en.wikipedia.org/wiki/Microsatellite
3. Richard GF, Kerrest A, Dujon B. 2008. Comparative genomics and molecular dynamics of DNA repeats in eukaryotes. Microbiology and Molecular Biology Reviews. 72(4):686-727.
4. Benson G. 1999. Tandem Repeats Finder: a program to analyze DNA sequences. Nucleic Acids Research. 27(2):573-580.
