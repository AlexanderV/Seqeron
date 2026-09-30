# Tandem Repeat Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-TANDEM-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production (`RepeatFinder.GetTandemRepeatSummary`: MISA/Krait statistics incl. MISA per-unit-size thresholds and repeat-type classes, misa.pl-verified); `GenomicAnalyzer.FindTandemRepeats` exact brute-force detector (see §5.4) |
| Last Reviewed | 2026-09-30 |

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
| INV-06 | The six per-class counts (mono … hexa) sum to `TotalRepeats`; `0 ≤ PercentageOfSequence ≤ 100`; covered bases ≤ `TotalRepeatBases`. | Every reported unit has length 1–6; coverage is the union of spans, the base total is their sum [5][6]. |
| INV-05 | Within a fixed candidate unit length, later starts inside a detected tandem block are skipped. | After yielding a result, the implementation advances the start index to the end of the detected tandem block for that unit-length pass. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` | required | DNA sequence to analyze. | `FindTandemRepeats` dereferences `sequence.Sequence` directly; `GetTandemRepeatSummary` throws on `null`. |
| `minUnitLength` | `int` | `2` | Minimum candidate repeat-unit length for `FindTandemRepeats`. | Values below `1` throw `ArgumentOutOfRangeException` (eager). |
| `minRepetitions` | `int` | `2` | Minimum number of consecutive unit copies for `FindTandemRepeats`. | Values below `2` throw `ArgumentOutOfRangeException` (eager). |
| `minRepeats` | `int` | `3` | Minimum repeat count used by `GetTandemRepeatSummary(DnaSequence, int)` for every unit length 1–6. | Passed to `FindMicrosatellites(sequence, 1, 6, minRepeats)`, which rejects values below 2. |
| `minRepeatsByUnitLength` | `IReadOnlyDictionary<int,int>` | — | Per-unit-size minimum copies for `GetTandemRepeatSummary(DnaSequence, IReadOnlyDictionary<int,int>)`, e.g. `RepeatFinder.MisaDefaultMinRepeats` = MISA `misa.ini` `1-10 2-6 3-5 4-5 5-5 6-5` [5]; unit lengths absent from the map are not searched. | Non-null, non-empty; unit lengths 1–6; values ≥ 2 (else `ArgumentNullException` / `ArgumentException` / `ArgumentOutOfRangeException`). |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Tandem repeats | `IEnumerable<TandemRepeat>` | Exact tandem-repeat hits with unit, 0-based start position, repetition count, total length, and full repeated sequence. |
| Repeat-type classes | `IReadOnlyDictionary<string,int>` | `RepeatFinder.GetCanonicalMotifFrequencies(ssrs)`: SSR count per MISA class "considering sequence complementary" — `X/Y` with X, Y the smallest rotations of the motif and of its reverse complement, smaller first (AC/CA/GT/TG → `AC/GT`) [5]; `GetStandardMotifFrequencies(ssrs, level)`: per Krait standard motif (`StandardMotif`, order A < T < C < G, levels 0–4; level 2 = rotations + reverse complement, e.g. ACAT → `ATAC`) [6]. |
| Tandem summary | `TandemRepeatSummary` | Aggregate summary over the perfect microsatellites (1–6 bp units) reported by `FindMicrosatellites`: `TotalRepeats`; per-class counts for mono-, di-, tri-, tetra-, penta- and hexanucleotide repeats (MISA "Distribution to different repeat type classes", Krait Mono…Hexa [5][6]); `TotalRepeatBases` = sum of repeat lengths (Krait "Length (bp)" = `SUM(length)` [6]; overlapping runs of different unit lengths each count in full); `PercentageOfSequence` = bases covered by the union of repeat spans / length × 100; `LongestRepeat` (largest `TotalLength`, ties → shorter unit then leftmost; `null` when none); `MostFrequentUnit` (reported unit string — motif phase at run start, not rotation/strand-canonicalized, as in MISA's "Frequency of identified SSR motifs" [5]; ties → first in unit-length/position order; `null` when none). |

### 3.3 Preconditions and Validation

`GetTandemRepeatSummary` throws `ArgumentNullException` when `sequence` is `null`. `FindTandemRepeats` throws `ArgumentNullException` for a `null` sequence and `ArgumentOutOfRangeException` for `minUnitLength < 1` or `minRepetitions < 2` (eager). Empty sequences produce no tandem-repeat hits, and an empty sequence summarized through `GetTandemRepeatSummary` returns zero totals, `0` percent coverage and `null` `LongestRepeat` / `MostFrequentUnit`. `GetTandemRepeatSummary` throws `ArgumentOutOfRangeException` for `minRepeats < 2` (eager validation in `FindMicrosatellites`).

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
- `RepeatFinder.GetTandemRepeatSummary(DnaSequence, IReadOnlyDictionary<int,int>)`: Same summary with MISA-style per-unit-size thresholds (`RepeatFinder.MisaDefaultMinRepeats`).
- `RepeatFinder.GetCanonicalMotifClass(string)` / `GetCanonicalMotifFrequencies(IEnumerable<MicrosatelliteResult>)`: MISA repeat-type classes (rotation + reverse complement).
- `RepeatFinder.GetStandardMotif(string, int level = 2)` / `GetStandardMotifFrequencies(IEnumerable<MicrosatelliteResult>, int level = 2)`: Krait standard motifs.
- Compound SSRs (MISA types `c` / `c*`): `RepeatFinder.FindCompoundMicrosatellites` / `AssembleCompoundMicrosatellites` — see [Microsatellite_Detection.md](Microsatellite_Detection.md) §5.

### 5.2 Current Behavior

`GenomicAnalyzer.FindTandemRepeats` uses a brute-force scan over candidate unit lengths and positions, compares units with direct substring equality, and skips forward after each hit within the current unit-length pass. This suppresses later starts inside the same detected block for that unit length, but it does not prevent the same region from being reported again under a different unit-length interpretation. It validates its parameters eagerly (2026-09, B09) but reports non-primitive units (e.g. `ATATATAT` → `AT×4@0` and `ATAT×2@0`); delegation to `RepeatFinder.FindMicrosatellites` is a recorded cross-batch request (B09 owns `GenomicAnalyzer.cs`). `RepeatFinder.GetTandemRepeatSummary` does validate `sequence`, then delegates to `FindMicrosatellites(sequence, 1, 6, minRepeats)`, meaning the summary covers only tandem repeats with 1-6 bp units and inherits that implementation's conventions (each maximal primitive ACGT run reported once per unit length; runs of different unit lengths may overlap). The summary record has a dedicated count field for each of the six classes, and the six counts sum to `TotalRepeats`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact detection of directly adjacent repeated units in DNA sequences [1][2].
- Standard short tandem repeat classification by unit length from mononucleotide through hexanucleotide classes [2][3].
- Reporting of tandem-repeat start position, unit, repetition count, and total repeated length for each exact hit.

**Intentionally simplified:**

- Brute-force direct substring comparison instead of suffix-tree, suffix-array, or Tandem Repeats Finder style optimization; **consequence:** runtime grows rapidly on long sequences and the implementation is best suited to moderate sequence lengths [1][4].
- `GetTandemRepeatSummary` is restricted to microsatellite-sized units from 1 to 6 bp; **consequence:** longer minisatellite and macrosatellite tandems are excluded from the summary even though `FindTandemRepeats` can detect longer exact units.

**Implemented 2026-09-30 (MISA / Krait parity):**

- Per-unit-size thresholds (`GetTandemRepeatSummary(DnaSequence, IReadOnlyDictionary<int,int>)`, `MisaDefaultMinRepeats` = `1-10 2-6 3-5 4-5 5-5 6-5`) [5]. **Cross-check:** SSR lists = brute-force maximal primitive runs with per-size thresholds on 6 048 sequences (6 `misa.ini` configurations, 72 974 SSRs), 0 mismatches; vs a real `perl misa.pl` run the per-class totals differ only through the documented REP-STR-001 run conventions (e.g. default ini: misa 6 639 SSRs, this library 6 641).
- MISA table "Frequency of classified repeat types (considering sequence complementary)" (`GetCanonicalMotifFrequencies`) [5] and Krait standard motifs (`GetStandardMotifFrequencies`, levels 0–4) [6]. **Cross-check:** class name = misa.pl `.statistics` row for all 5 356 primitive motifs of 1–6 bp (one misa.pl run each); Krait `StandardMotif.standard()` for all 5 460 motifs × 5 levels, 0 mismatches; the class table built from misa.pl's own SSR list equals misa.pl's table in all 6 configurations.
- MISA compound SSRs: see Microsatellite_Detection.md §5.
- `MostFrequentUnit` remains the raw reported unit (MISA "Frequency of identified SSR motifs"); the canonicalized counts are the separate class table above.

**Not implemented:**

- Suffix-structure indexing inside `GenomicAnalyzer.FindTandemRepeats` (B09-owned brute-force scan; the delegation request routes it to the O(n)-per-unit-length `RepeatFinder.FindMicrosatellites`). Approximate tandem repeats are provided by `RepeatFinder.FindApproximateTandemRepeats` (TRF model, REP-APPROX-001).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `GenomicAnalyzer.FindTandemRepeats` reports non-primitive units and one hit per unit-length interpretation (no primitivity / maximal-run rule), O(n²·m) substring scan. | Deviation | `ATATATAT` → `AT×4` and `ATAT×2`. | cross-batch request (B09 owns `GenomicAnalyzer.cs`): delegate to `RepeatFinder.FindMicrosatellites` | Parameter validation was added (2026-09); the former "no validation" assumption no longer applies. |

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

`GenomicAnalyzer.FindTandemRepeats` is exact and does not score approximate tandem repeats, interrupted repeats, or noisy repeat families (use `RepeatFinder.FindApproximateTandemRepeats`, TRF model). The summary helper is narrower than the canonical detector because it only considers 1-6 bp units. The detector also does not canonicalize across competing unit-length interpretations, so the same genomic region can appear more than once when different repeat-unit sizes satisfy the threshold. There is also no raw-string overload for `FindTandemRepeats`, so callers must provide a `DnaSequence` and handle any normalization before calling the algorithm.

## 7. Examples and Related Material

### 7.2 Related Use Cases

- Forensic DNA profiling: STR markers are standard forensic markers, with tetra- and pentanucleotide loci commonly preferred for robust genotyping [2].
- Paternity and kinship analysis: High STR polymorphism makes tandem repeats useful for relationship inference [2].
- Population genetics: Tandem-repeat variability supports diversity and lineage studies [2][3].
- Cancer diagnostics: Microsatellite instability is a repeat-based marker used in oncology workflows [1][2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomicAnalyzer_TandemRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_TandemRepeat_Tests.cs), [RepeatFinder_MisaCompound_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_MisaCompound_Tests.cs)
- Test spec: [REP-TANDEM-001.md](../../../tests/TestSpecs/REP-TANDEM-001.md)
- Related property tests: [RepeatFinderProperties.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Properties/RepeatFinderProperties.cs)
- Related metamorphic tests: [MetamorphicTests.cs](../../../tests/SuffixTree/SuffixTree.Tests/Algorithms/MetamorphicTests.cs)

### 7.4 Change History

| Date | Version | Author | Changes |
|------|---------|--------|---------|
| 2026-01-22 | 1.0 | Algorithm QA | Initial documentation |
| 2026-09-29 | 1.1 | Review 2026-09 (B04) | Summary: penta/hexa classes, null `LongestRepeat`, sourced field conventions (MISA / Krait) |
| 2026-09-30 | 1.2 | Review 2026-09 (B04 audit WP3) | MISA per-unit-size thresholds, MISA repeat-type classes, Krait standard motifs (misa.pl / Krait verified); status → Production |

## 8. References

1. Wikipedia. 2026. Tandem repeat. Wikipedia. https://en.wikipedia.org/wiki/Tandem_repeat
2. Wikipedia. 2026. Microsatellite. Wikipedia. https://en.wikipedia.org/wiki/Microsatellite
3. Richard GF, Kerrest A, Dujon B. 2008. Comparative genomics and molecular dynamics of DNA repeats in eukaryotes. Microbiology and Molecular Biology Reviews. 72(4):686-727.
4. Benson G. 1999. Tandem Repeats Finder: a program to analyze DNA sequences. Nucleic Acids Research. 27(2):573-580.
5. Thiel T, Michalek W, Varshney RK, Graner A. 2003. Exploiting EST databases for the development and characterization of gene-derived SSR-markers in barley. Theoretical and Applied Genetics 106:411-422. MISA `misa.pl` v1.0 source (`.statistics` output sections), opened via a raw GitHub mirror 2026-09-29.
6. Du L, Zhang C, Liu Q, Zhang X, Yue B. 2018. Krait: an ultrafast tool for genome-wide survey of microsatellites and primer design. Bioinformatics 34(4):681-683. `src/statistics.py` (lmdu/krait, raw.githubusercontent.com, opened 2026-09-29).
