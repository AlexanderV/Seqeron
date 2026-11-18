# Tandem Repeat Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-TANDEM-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production (`RepeatFinder.GetTandemRepeatSummary`: MISA/Krait statistics incl. MISA per-unit-size thresholds and repeat-type classes, misa.pl-verified); `GenomicAnalyzer.FindTandemRepeats` exact maximal-run detector (facade over `RepeatFinder.FindMicrosatellites` + a run-enumeration engine for wide unit ranges, B09 F35–F36; see GENOMIC-TANDEM-001) |
| Last Reviewed | 2026-09-30 |

## 1. Overview

Tandem repeat detection identifies contiguous DNA segments where a repeat unit occurs consecutively two or more times, such as `ATTCGATTCGATTCG` for the unit `ATTCG` repeated three times [1]. The repository exposes an exact tandem detector in `GenomicAnalyzer.FindTandemRepeats` and a separate aggregation helper, `RepeatFinder.GetTandemRepeatSummary`, that summarizes only microsatellite-sized tandem repeats with unit lengths from 1 to 6 bp [1][2]. The canonical detector reports, for each unit length, every maximal run with a primitive unit once (a facade over `RepeatFinder.FindMicrosatellites`; see [Genomic_Analysis/Tandem_Repeat_Detection.md](../Genomic_Analysis/Tandem_Repeat_Detection.md)), though the same region can still be reported under a different, non-multiple unit-length interpretation. Tandem repeats are biologically important because they occupy a substantial fraction of the human genome and include disease-associated repeat expansions and microsatellite-instability markers [1][2].

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
| INV-01 | Every reported tandem repeat contains at least `minRepetitions` consecutive copies of a unit whose length is at least `minUnitLength` (and at most `maxUnitLength`). | Only maximal runs with `⌊length / p⌋ ≥ minRepetitions` are reported. |
| INV-02 | `TotalLength = Unit.Length × Repetitions` for every `TandemRepeat`. | Total length is defined by the unit size and repetition count. |
| INV-03 | `Position + Unit.Length × Repetitions <= sequence.Length`. | Only complete copies inside the run are counted. |
| INV-04 | The summary percentages and totals are derived only from reported microsatellites. | `GetTandemRepeatSummary` delegates to `FindMicrosatellites(sequence, 1, 6, minRepeats)`. |
| INV-06 | `CountsByUnitLength` (one entry per searched unit size) sums to `TotalRepeats`; its entries 1–6 equal the six named class counts, which therefore sum to `TotalRepeats` whenever no size > 6 is searched (always for the uniform overloads); `0 ≤ PercentageOfSequence ≤ 100`; covered bases ≤ `TotalRepeatBases`. | Each reported SSR is counted once under its unit length (misa.pl `$count_class{size}++`); coverage is the union of spans, the base total is their sum [5][6]. |
| INV-05 | Within a unit length, reported runs are maximal and do not overlap (the next run starts after `e − p`); the unit is primitive. | Left-maximal start + right-maximal end + primitivity test (Kolpakov & Kucherov 1999). |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to analyze (string overloads of `GetTandemRepeatSummary` accept N/IUPAC, case-insensitive). | `FindTandemRepeats` throws `ArgumentNullException` on `null`; `GetTandemRepeatSummary(DnaSequence, …)` throws on `null`; the string overloads return the empty summary for `null`/empty. |
| `minUnitLength` | `int` | `2` | Minimum candidate repeat-unit length for `FindTandemRepeats`. | Values below `1` throw `ArgumentOutOfRangeException` (eager). |
| `minRepetitions` | `int` | `2` | Minimum number of consecutive unit copies for `FindTandemRepeats`. | Values below `2` throw `ArgumentOutOfRangeException` (eager). |
| `minRepeats` | `int` | `3` | Minimum repeat count used by `GetTandemRepeatSummary(DnaSequence, int)` for every unit length 1–6. | Passed to `FindMicrosatellites(sequence, 1, 6, minRepeats)`, which rejects values below 2. |
| `minRepeatsByUnitLength` | `IReadOnlyDictionary<int,int>` | — | Per-unit-size minimum copies for `GetTandemRepeatSummary(DnaSequence, IReadOnlyDictionary<int,int>)`, e.g. `RepeatFinder.MisaDefaultMinRepeats` = MISA `misa.ini` `1-10 2-6 3-5 4-5 5-5 6-5` [5]; unit lengths absent from the map are not searched; any unit length ≥ 1 (misa.pl accepts any size in its `def` line, e.g. `… 7-5 8-5`); `RepeatFinder.ParseMisaDefinition(string)` builds the map from misa.ini `def` syntax. | Non-null, non-empty; unit lengths ≥ 1; values ≥ 2 (else `ArgumentNullException` / `ArgumentException` / `ArgumentOutOfRangeException`). |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Tandem repeats | `IEnumerable<TandemRepeat>` | Exact tandem-repeat hits with unit, 0-based start position, repetition count, total length, and full repeated sequence. |
| Repeat-type classes | `IReadOnlyDictionary<string,int>` | `RepeatFinder.GetCanonicalMotifFrequencies(ssrs)`: SSR count per MISA class "considering sequence complementary" — `X/Y` with X, Y the smallest rotations of the motif and of its reverse complement, smaller first (AC/CA/GT/TG → `AC/GT`) [5]; `GetStandardMotifFrequencies(ssrs, level)`: per Krait standard motif (`StandardMotif`, order A < T < C < G, levels 0–4; level 2 = rotations + reverse complement, e.g. ACAT → `ATAC`) [6]. |
| Tandem summary | `TandemRepeatSummary` | Aggregate summary over the perfect microsatellites reported by `FindMicrosatellites` (unit sizes 1–6 for the uniform overloads, the map's sizes — any size — for the per-unit-length overloads): `TotalRepeats`; `CountsByUnitLength` = count per searched unit size (MISA "Distribution to different repeat type classes", every size incl. > 6; 0 for searched sizes without SSRs, which misa.pl omits) and the named per-class counts for mono-, di-, tri-, tetra-, penta- and hexanucleotide repeats (Krait Mono…Hexa [5][6]); `TotalRepeatBases` = sum of repeat lengths (Krait "Length (bp)" = `SUM(length)` [6]; overlapping runs of different unit lengths each count in full); `PercentageOfSequence` = bases covered by the union of repeat spans / length × 100; `LongestRepeat` (largest `TotalLength`, ties → shorter unit then leftmost; `null` when none); `MostFrequentUnit` (reported unit string — motif phase at run start, not rotation/strand-canonicalized, as in MISA's "Frequency of identified SSR motifs" [5]; ties → first in unit-length/position order; `null` when none). |

### 3.3 Preconditions and Validation

`GetTandemRepeatSummary` throws `ArgumentNullException` when `sequence` is `null`. `FindTandemRepeats` throws `ArgumentNullException` for a `null` sequence and `ArgumentOutOfRangeException` for `minUnitLength < 1`, `minRepetitions < 2` or (4-argument overload) `maxUnitLength < minUnitLength` (eager). Empty sequences produce no tandem-repeat hits, and an empty sequence summarized through `GetTandemRepeatSummary` returns zero totals, `0` percent coverage and `null` `LongestRepeat` / `MostFrequentUnit`. `GetTandemRepeatSummary` throws `ArgumentOutOfRangeException` for `minRepeats < 2` (eager validation in `FindMicrosatellites`).

## 4. Algorithm

### 4.1 High-Level Steps

1. For each unit length `p` from `minUnitLength` to `min(maxUnitLength, sequence.Length / minRepetitions)`, visit only the left-maximal run starts (`S[i-1] ≠ S[i-1+p]`) and extend each to its right-maximal end (`RepeatFinder.FindMicrosatellites`, O(n) comparisons per unit length).
2. Report the run once, at its left end, with `⌊length / p⌋` complete copies, when the unit `S[i..i+p)` is primitive and the copies reach `minRepetitions`.
3. When more than 256 unit lengths are admissible (e.g. the default "no upper bound" on a long sequence), enumerate all maximal repetitions once instead (`MaximalRepetitions`: Lyndon-array runs algorithm, Bannai et al. 2017, O(n log² n)) and keep those in range; the result is identical.
5. For summary mode, detect microsatellites with unit sizes 1-6 (or the sizes of a per-unit-length MISA definition) and aggregate counts per unit size, bases covered, repeat-type totals, longest repeat, and most frequent unit.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `GenomicAnalyzer.FindTandemRepeats` | `O(n × K)` for `K ≤ 256` admissible unit lengths, else `O(n log² n)` | `O(1)` plus output, resp. `O(n)` | 1 Mbp: ≈ 1 s per ten unit lengths via the scan; unbounded unit length via runs, see the GENOMIC-TANDEM-001 evidence. Was `O(n² × m)` substring comparisons. |
| `RepeatFinder.GetTandemRepeatSummary` | `O(n × U × R)` | `O(k)` | Delegates to microsatellite detection, where `U` is the searched unit-length range, `R` is the average repeat count, and `k` is the number of microsatellites retained. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation locations:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs), [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `GenomicAnalyzer.FindTandemRepeats(DnaSequence, int, int)` and `(…, int maxUnitLength)`: canonical exact detector for consecutive tandem repeats (maximal runs, primitive units; B09 F35–F36).
- `RepeatFinder.GetTandemRepeatSummary(DnaSequence, int)`: Summary helper that aggregates microsatellite-sized tandem repeats.
- `RepeatFinder.GetTandemRepeatSummary(DnaSequence, IReadOnlyDictionary<int,int>)`: Same summary with MISA-style per-unit-size thresholds (`RepeatFinder.MisaDefaultMinRepeats`). Any unit size (e.g. a misa.ini `def` line with 7–10 bp units, parsed by `RepeatFinder.ParseMisaDefinition`); per-size counts in `CountsByUnitLength` (B04 WP16, F61).
- `RepeatFinder.GetTandemRepeatSummary(string, int)`, `(string, IReadOnlyDictionary<int,int>)`, `(string, IReadOnlyDictionary<int,int>, MicrosatelliteScanMode)`: raw-string, N/IUPAC-tolerant counterparts (case-insensitive; only A/C/G/T form units, so N/IUPAC never belong to an SSR — MISA `[acgt]`); `PercentageOfSequence` uses the full length incl. N (misa.pl `length $seq`); `null`/empty → empty summary. With `MisaRegex` the totals and per-unit-size counts equal misa.pl `.statistics` (6 000 N-containing sequences × 6 definitions, 0 mismatches; 12 000 per-sequence runs identical).
- `RepeatFinder.GetCanonicalMotifClass(string)` / `GetCanonicalMotifFrequencies(IEnumerable<MicrosatelliteResult>)`: MISA repeat-type classes (rotation + reverse complement).
- `RepeatFinder.GetStandardMotif(string, int level = 2)` / `GetStandardMotifFrequencies(IEnumerable<MicrosatelliteResult>, int level = 2)`: Krait standard motifs.
- Compound SSRs (MISA types `c` / `c*`): `RepeatFinder.FindCompoundMicrosatellites` / `AssembleCompoundMicrosatellites` — see [Microsatellite_Detection.md](Microsatellite_Detection.md) §5.

### 5.2 Current Behavior

`GenomicAnalyzer.FindTandemRepeats` validates its parameters eagerly and delegates to `RepeatFinder.FindMicrosatellites(sequence, minUnitLength, maxUnitLength, minRepetitions)` (maximal run, primitive unit, one report per run and unit length; `ATATATAT` → `AT×4@0` only), switching to the run enumeration `MaximalRepetitions` when more than 256 unit lengths are admissible (B09 F35–F36, oracle-verified: 1156 Python rows and 1500 + 1200 brute-force strings). `RepeatFinder.GetTandemRepeatSummary` does validate `sequence`, then delegates to `FindMicrosatellites(sequence, 1, 6, minRepeats)`, meaning the summary covers only tandem repeats with 1-6 bp units and inherits that implementation's conventions (each maximal primitive ACGT run reported once per unit length; runs of different unit lengths may overlap). The summary record has a dedicated count field for each of the six classes, and the six counts sum to `TotalRepeats`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact detection of directly adjacent repeated units in DNA sequences [1][2].
- Standard short tandem repeat classification by unit length from mononucleotide through hexanucleotide classes [2][3].
- Reporting of tandem-repeat start position, unit, repetition count, and total repeated length for each exact hit.

**Intentionally simplified:**

- ~~`GetTandemRepeatSummary` is restricted to microsatellite-sized units from 1 to 6 bp~~ — **resolved 2026-10-01 (B04 audit WP16, F61):** the uniform-threshold overloads keep the classical 1–6 bp STR window, but the per-unit-length (MISA definition) overloads accept any unit size as misa.pl does (`def` = free list of `size-min` pairs, one regex per size; `.statistics` distribution row per size), with `TandemRepeatSummary.CountsByUnitLength` carrying every size (the six named fields stay 1–6) and MCP `misaDefinition` on `find_microsatellites` / `tandem_repeat_summary`. **Cross-check:** real `perl misa.pl` on 6 000 sequences × 7 definitions with unit sizes 7–12 (434 008 SSRs, 107 762 with units > 6; 79 690 `.misa` rows incl. compounds): SSR lists, rows and per-sequence per-size counts 0 mismatches (ties by SSR number, F48); `.statistics` totals and every distribution row identical; 2 100 / 2 100 per-sequence stock runs identical.

**Implemented 2026-09-30 (MISA / Krait parity):**

- Per-unit-size thresholds (`GetTandemRepeatSummary(DnaSequence, IReadOnlyDictionary<int,int>)`, `MisaDefaultMinRepeats` = `1-10 2-6 3-5 4-5 5-5 6-5`) [5]. **Cross-check:** SSR lists = brute-force maximal primitive runs with per-size thresholds on 6 048 sequences (6 `misa.ini` configurations, 72 974 SSRs), 0 mismatches; vs a real `perl misa.pl` run the per-class totals differ only through the documented REP-STR-001 run conventions (e.g. default ini: misa 6 639 SSRs, this library 6 641).
- MISA table "Frequency of classified repeat types (considering sequence complementary)" (`GetCanonicalMotifFrequencies`) [5] and Krait standard motifs (`GetStandardMotifFrequencies`, levels 0–4) [6]. **Cross-check:** class name = misa.pl `.statistics` row for all 5 356 primitive motifs of 1–6 bp (one misa.pl run each); Krait `StandardMotif.standard()` for all 5 460 motifs × 5 levels, 0 mismatches; the class table built from misa.pl's own SSR list equals misa.pl's table in all 6 configurations.
- MISA compound SSRs: see Microsatellite_Detection.md §5.
- `MostFrequentUnit` remains the raw reported unit (MISA "Frequency of identified SSR motifs"); the canonicalized counts are the separate class table above.

**Not implemented:**

- Approximate (interrupted) tandem repeats are not part of `GenomicAnalyzer.FindTandemRepeats`; they are `RepeatFinder.FindApproximateTandemRepeats` (TRF 4.10.0 model, REP-APPROX-001) and the MCP tools `find_approximate_tandem_repeats` / `mask_approximate_tandem_repeats`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | ~~`GenomicAnalyzer.FindTandemRepeats` reported non-primitive units (`ATATATAT` → `AT×4` and `ATAT×2`), O(n²·m).~~ | Deviation | — | **resolved 2026-11-18 (B09 F35–F36)**: facade over `RepeatFinder.FindMicrosatellites` + run enumeration | See [Genomic_Analysis/Tandem_Repeat_Detection.md](../Genomic_Analysis/Tandem_Repeat_Detection.md). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | `FindTandemRepeats` yields no results; summary returns zero totals. | No unit length admits `minRepetitions` copies. |
| No repeats found | Empty tandem-repeat enumerable or zero-count summary. | No candidate unit reaches the repetition threshold. |
| Entire sequence is one tandem block | One result spanning the full repeated region. | The maximal run is the whole sequence. |
| Overlapping tandem patterns | A run is reported once per unit length (rotations are not repeated); runs of different unit lengths may overlap and are each reported. Multiples of the primitive unit are not reported (`ATATATAT`: `AT×4` only). | Maximal-run + primitive-unit rule (Kolpakov & Kucherov 1999; MISA). |
| `unitLength > sequence.Length / minRepetitions` | No result for that unit size. | `minRepetitions` copies cannot fit. |

### 6.2 Limitations

`GenomicAnalyzer.FindTandemRepeats` is exact and does not score approximate tandem repeats, interrupted repeats, or noisy repeat families (use `RepeatFinder.FindApproximateTandemRepeats`, TRF model, or the MCP tool `find_approximate_tandem_repeats`). The uniform-threshold summary considers 1-6 bp units only; for longer units pass a per-unit-length (MISA) definition with the wanted sizes. Runs of different (non-multiple) unit lengths over the same region are each reported. There is no raw-string overload for `FindTandemRepeats` (a `DnaSequence` is A/C/G/T only); `RepeatFinder.FindMicrosatellites(string, …)` accepts N/IUPAC text.

## 7. Examples and Related Material

### 7.2 Related Use Cases

- Forensic DNA profiling: STR markers are standard forensic markers, with tetra- and pentanucleotide loci commonly preferred for robust genotyping [2].
- Paternity and kinship analysis: High STR polymorphism makes tandem repeats useful for relationship inference [2].
- Population genetics: Tandem-repeat variability supports diversity and lineage studies [2][3].
- Cancer diagnostics: Microsatellite instability is a repeat-based marker used in oncology workflows [1][2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomicAnalyzer_TandemRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_TandemRepeat_Tests.cs), [GenomicAnalyzer_TandemRepeat_Oracle_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzer_TandemRepeat_Oracle_Tests.cs), [RepeatFinder_MisaCompound_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_MisaCompound_Tests.cs)
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
