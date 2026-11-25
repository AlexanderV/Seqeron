# Repetitive Element Detection and Classification

| Field | Value |
|-------|-------|
| Algorithm Group | Annotation |
| Test Unit ID | ANNOT-REPEAT-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Complete (exact detection via RepeatFinder; ClassifyRepeat exact + Smith-Waterman-Gotoh options) |
| Last Reviewed | 2026-11-25 (B11 review: F14–F17) |

## 1. Overview

Detects repetitive elements in a DNA sequence and assigns a repeat class. `FindRepetitiveElements` reports two structural classes: **tandem repeats** (a motif repeated head-to-tail in two or more directly adjacent copies) and **inverted repeats** (a sequence followed downstream by its reverse complement, the stem-loop form WGW^R). `ClassifyRepeat` assigns a RepeatMasker-style class label by matching a query against a supplied repeat-element library, falling back to a simple-repeat class by motif size. Detection is exact (no mismatches) and delegates to the canonical `RepeatFinder` (maximal runs; EMBOSS palindrome maximal stems); classification is library-driven (exact both-strand containment, or scored Smith-Waterman-Gotoh homology) and therefore as complete as the library supplied.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Repetitive DNA makes up a large fraction of eukaryotic genomes. Tandem repeats (microsatellites/STRs with 1–6 bp motifs, minisatellites/VNTRs with ~10–60 bp motifs) arise by replication slippage; inverted repeats form hairpin/cruciform secondary structures and are recognition sites for many regulatory and restriction enzymes [1][2]. Interspersed repeats (SINE, LINE, LTR, DNA transposons) are classified by homology to curated libraries such as Repbase [4].

### 2.2 Core Model

**Tandem repeat.** "A pattern of one or more nucleotides is repeated and the repetitions are directly adjacent to each other" (head-to-tail), with a minimum of two copies [1][5]. An exact tandem repeat of period `p` is a **maximal run** (maximal repetition, Kolpakov & Kucherov 1999 [7]): the longest interval `S[a..e)` with `S[x] = S[x+p]` for all `a ≤ x < e−p`, not extendable left or right. It is reported once, at its left end `a`, with unit `S[a..a+p)` and `c = ⌊(e−a)/p⌋` complete copies, spanning `c·p` bases (MISA / pytrf convention); rotations of the same run (TATATA inside CATATATA) are the same run. The unit must be **primitive** (non-primitive units such as "AA" = "A"×2 are reported at the primitive period) and consist of A/C/G/T only (MISA `[acgt]`; runs of N are assembly gaps, not repeats).

**Inverted repeat.** An inverted repeat is "a single stranded sequence of nucleotides followed downstream by its reverse complement" [2]. Formally it has the form **WGW̄ᴿ**, where `W` is the left arm, `G` is the gap/loop/spacer with `|G| ≥ 0`, and `W̄ᴿ` is the reverse complement of `W` (the right arm) [3]. When `|G| = 0` the element is a reverse-complement palindrome [2][3].

Among all exact stems, only those **not contained in both arms in another stem** are reported (EMBOSS `palindrome -nummismatches 0 -overlap Y` [6]); sub-stems and slipped re-pairings inside a longer stem are the same structure. Default limits are the EMBOSS `palindrome` / IUPACpal defaults: arm ≤ 100 (`-maxpallen` / `-M`), loop ≤ 100 (`-gaplimit` / `-g`) [3][6]; a stem longer than the arm cap is not reported and not split (EMBOSS `-maxpallen` semantics).

**Repeat classification.** RepeatMasker assigns a class by the best homology match (Smith-Waterman-Gotoh via cross_match/RMBlast, scored against a Repbase/Dfam library, on both strands — output orientation "+" or "C") above a minimum score (`-cutoff`); classes are SINE, LINE, LTR, DNA, Satellite, Simple_repeat, Low_complexity, RNA, and Unclassified/Unknown [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A reported tandem repeat's `sequence` equals `dnaSequence[start..end]` (0-based, end-exclusive) and its length is an integer multiple of the unit length | reported span is `start..start+copies·|unit|`, sequence taken as that slice [1] |
| INV-02 | A reported tandem repeat has ≥ `minCopies` (≥ 2) directly adjacent identical primitive copies | copy count is incremented only on exact adjacent unit match; arrays below `minCopies` are dropped [1] |
| INV-03 | For every reported inverted repeat, the right arm equals the reverse complement of the left arm (A·T / C·G pairs only) and the stem is maximal | `RepeatFinder.FindInvertedRepeats` pairs via `SequenceExtensions.GetComplementBase`, ACGT only, `palindrome_AInB` filter [2][3][6] |
| INV-04 | `ClassifyRepeat` returns a label from the library or the fixed fallback set {Simple_repeat, Unknown} | return paths are library value, `Simple_repeat`, or `Unknown` [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| dnaSequence | string | required | DNA sequence to scan | non-null; case-insensitive (upper-cased internally) |
| minRepeatLength | int | 10 | tandem: minimum total span (bp); inverted: minimum arm length (raised to 2, RepeatFinder's minimum stem) | ≥ 1 |
| minCopies | int | 2 | minimum adjacent copies for a tandem repeat | ≥ 2 |
| sequence (ClassifyRepeat) | string | required | repeat sequence to classify | non-null; case-insensitive |
| repeatDb (ClassifyRepeat) | IReadOnlyDictionary<string,string> | required | library: element sequence → class label | non-null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| start | int | 0-based inclusive start of the element |
| end | int | exclusive end of the element |
| type | string | `"tandem_repeat"` or `"inverted_repeat"` |
| sequence | string | tandem: the array slice; inverted: `arm1arm2` (gap 0) or `arm1[gap]arm2` |
| (ClassifyRepeat) | string | library class label, or `"Simple_repeat"` / `"Unknown"` |

### 3.3 Preconditions and Validation

Coordinates are 0-based with exclusive end. The alphabet is DNA; input is upper-cased. `FindRepetitiveElements` throws `ArgumentNullException` for null `dnaSequence`, `ArgumentOutOfRangeException` for `minCopies < 2` or `minRepeatLength < 1`; an empty sequence yields an empty result. `ClassifyRepeat` throws `ArgumentNullException` for null `sequence` or `repeatDb`; an empty sequence returns `Simple_repeat`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; upper-case the sequence.
2. **Tandem repeats:** `RepeatFinder.FindMicrosatellites(seq, 1, 60, minCopies)` (one O(n) scan per period: visit run starts only, extend to the right-maximal end, report the primitive ACGT unit with its complete copies); keep runs with span ≥ `minRepeatLength`.
3. **Inverted repeats:** `RepeatFinder.FindInvertedRepeats(seq, max(minRepeatLength, 2), 100, minLoopLength: 0, maxArmLength: max(100, arm minimum))` — maximal exact stems (EMBOSS palindrome model); the composite is `leftArm + rightArm` or `leftArm[loop]rightArm`.
4. **Classification (exact):** the class of the longest library element occurring in the query on either strand (element or its reverse complement ⊆ query); **(scored):** the class of the library element with the best Smith-Waterman-Gotoh score (`SequenceAligner.LocalAlignAffine`, both strands) ≥ `minScore`, ties to the longer element; else `Simple_repeat` if the whole query is ≥ 2 complete copies of a primitive ACGT motif of 1–6 bp (`RepeatFinder.FindMicrosatellites`), else `Unknown`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Tandem unit length capped at 60 bp (covers STR 1–6 bp and minisatellite 10–60 bp) [1]; inverted arm capped at 100 bp and loop at 100 nt (EMBOSS palindrome `-maxpallen` / `-gaplimit` and IUPACpal `-M` / `-g` defaults [3][6]).
- Class vocabulary (Simple_repeat, Unknown, plus library labels SINE/LINE/LTR/DNA/Satellite) from RepeatMasker output classes [4].
- Primitivity and ACGT-only unit tests are RepeatFinder's (`IsReportableUnit`); no local copy.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindRepetitiveElements (tandem) | O(60·n) | O(r) | n = length, r = reported runs; one linear scan per period |
| FindRepetitiveElements (inverted) | O(n · G · A) | O(r) | G = loop limit 100, A = stem length (maximality check) |
| ClassifyRepeat (exact) | O(d · (m + q)) | O(m) | d = library size, m = element length, q = query length (substring scan, 2 strands) |
| ClassifyRepeat (scored) | O(d · m · q) | O(m · q) | Smith-Waterman-Gotoh per element and strand |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomeAnnotator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs)

- `GenomeAnnotator.FindRepetitiveElements(string, int, int)`: detects tandem and inverted repeats.
- `GenomeAnnotator.ClassifyRepeat(string, IReadOnlyDictionary<string,string>)`: assigns a repeat class via exact both-strand library containment with simple-repeat / Unknown fallback.
- `GenomeAnnotator.ClassifyRepeat(string, IReadOnlyDictionary<string,string>, ScoringMatrix, int minScore)`: scored Smith-Waterman-Gotoh homology (both strands, best score ≥ `minScore`).

### 5.2 Current Behavior

Detection is exact-match only (no mismatches/IUPAC degeneracy in arms). Tandem runs are reported once at their left end via their **primitive** ACGT unit, so `AAAAAA` is reported once as motif `A` and `CATATATA` once as `ATATAT` [1,7). Inverted repeats are the maximal stems only (`AGGATCCT` is one arm-4 stem, not also `AGG[AT]CCT` / `GGATCC`); the composite `sequence` field embeds the loop in square brackets when it is non-empty. Cross-checked (B11 F14/F15) on 1,860 random/planted cases: tandem 6,677 runs vs an independent maximal-run oracle, 0 differences; inverted 106,801 stems vs EMBOSS 6.6.0 `palindrome` on ACGT-only cases and 98,864 stems vs a brute-force maximal-stem oracle (incl. N), 0 differences.

**Suffix tree evaluation (search-reuse policy).** The repository [SuffixTree](../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.cs) supports O(m) occurrence queries for a *known* pattern after O(n) construction. It was evaluated and **not used** for `FindRepetitiveElements`: repeat detection here is not a fixed-pattern occurrence query — it must *discover* unknown primitive units and check head-to-tail adjacency (tandem) or pair a left arm with the occurrence of its reverse complement within a bounded gap (inverted). These are structural-discovery problems, not "find occurrences of pattern X in text Y", so a direct suffix-tree lookup does not yield the required output without substantial additional machinery; a bounded O(n²) scan is the correct fit for the in-scope input sizes. For `ClassifyRepeat`, library matching is conceptually a homology (scored alignment) problem (RepeatMasker uses Smith-Waterman-Gotoh [4]); the implemented exact-substring containment is a simplification (§5.3) for which the suffix tree is also unnecessary at the small library/query sizes intended.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Tandem repeat = head-to-tail, ≥ 2 directly adjacent copies, primitive unit reported [1].
- Inverted repeat detected as WGW̄ᴿ with right arm = exact reverse complement of left arm, `|G| ≥ 0` (gap 0 ⇒ palindrome) [2][3].
- Repeat-class vocabulary (SINE/LINE/LTR/DNA/Satellite/Simple_repeat/Unknown) from RepeatMasker output classes [4].

**Scope (by design, other entry points exist):**

- `ClassifyRepeat(sequence, repeatDb)` matches by exact containment; scored homology is the `ScoringMatrix` overload (Smith-Waterman-Gotoh, both strands, `minScore` cutoff). RepeatMasker's own GC-dependent transition/transversion matrices and complexity-adjusted scores are not reproduced (`ScoringMatrix` is a match/mismatch scheme); the scheme and cutoff are the caller's choice [4].
- Detection is perfect-match only; approximate tandem repeats are `RepeatFinder.FindApproximateTandemRepeats` (TRF), k-mismatch inverted repeats `RepeatFinder.FindInvertedRepeats(maxMismatches)` / `FindInvertedRepeatsScored` (EMBOSS einverted).

**Not implemented:**

- Curated Repbase/Dfam library, low-complexity masking, and interspersed-element family annotation; **users should rely on:** external tools (RepeatMasker) supplying a real library to `ClassifyRepeat`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Exact overload = containment (both strands); scored overload = SW-Gotoh with a match/mismatch scheme | Scope | RM matrices / complexity adjustment not reproduced | accepted | §5.3; Evidence |
| 2 | Inverted-repeat arm minimum is max(minRepeatLength, 2) | Deviation | minRepeatLength 1 reports no arm-1 stems (EMBOSS `-minpallen` allows 1) | B04 request R15 (QUEUE.md) | RepeatFinder validates `minArmLength ≥ 2` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null dnaSequence | `ArgumentNullException` | input contract |
| empty dnaSequence | empty result | nothing to scan |
| minCopies < 2 | `ArgumentOutOfRangeException` | tandem repeat needs ≥ 2 copies [1] |
| single motif occurrence | not reported as tandem repeat | requires "two or more" [1] |
| non-primitive unit (AA) | collapsed to primitive (A) | avoids double-counting [1] |
| inverted repeat gap 0 | reported (palindrome) | `|G| = 0` allowed [2][3] |
| rotation of a tandem run (CATATATA) | one run [1,7) `ATATAT` | maximal run [7] |
| run of N | no tandem / inverted repeat | ACGT-only units and pairs [6] |
| sub-stem of a longer stem | not reported | EMBOSS palindrome `-overlap Y` [6] |
| stem arm > 100 (minRepeatLength ≤ 100) | not reported, not split | EMBOSS `-maxpallen` [6] |
| query holds revcomp of a library element | that element's class | RepeatMasker "C" orientation [4] |
| ClassifyRepeat no library match, simple motif | `Simple_repeat` | STR 1–6 bp [1] |
| ClassifyRepeat no match, non-simple | `Unknown` | RepeatMasker Unclassified [4] |

### 6.2 Limitations

No mismatch/IUPAC tolerance in `FindRepetitiveElements` (use the RepeatFinder approximate APIs); classification quality is bounded by the supplied library and scoring scheme.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var elements = GenomeAnnotator.FindRepetitiveElements("ATTCGATTCGATTCG", minRepeatLength: 5, minCopies: 2);
// -> (0, 15, "tandem_repeat", "ATTCGATTCGATTCG")  // unit "ATTCG" x3

var db = new Dictionary<string, string> { ["GGCCGGGCGCGGTGGCTCAC"] = "SINE/Alu" };
var cls = GenomeAnnotator.ClassifyRepeat("...GGCCGGGCGCGGTGGCTCAC...", db); // "SINE/Alu"
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomeAnnotator_FindRepetitiveElements_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/GenomeAnnotator_FindRepetitiveElements_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ANNOT-REPEAT-001-Evidence.md](../../../docs/Evidence/ANNOT-REPEAT-001-Evidence.md)
- Related algorithms: [Microsatellite detection](../Repeat_Analysis/Microsatellite_Detection.md) (`RepeatFinder.FindMicrosatellites`, `FindInvertedRepeats`)

## 8. References

1. Wikipedia contributors. 2026. *Tandem repeat*. Wikipedia. https://en.wikipedia.org/wiki/Tandem_repeat (cites Duitama J et al. 2014, *Nucleic Acids Res* 42(9):5728–5741, https://doi.org/10.1093/nar/gku212).
2. Wikipedia contributors. 2026. *Inverted repeat*. Wikipedia. https://en.wikipedia.org/wiki/Inverted_repeat (cites Ye C, Ji G, Liang C. 2014. *PLoS ONE* 9(11):e113349, https://doi.org/10.1371/journal.pone.0113349).
3. Hampson SE, Pissis SP, et al. 2021. *IUPACpal: efficient identification of inverted repeats in IUPAC-encoded DNA sequences*. BMC Bioinformatics 22:51. https://doi.org/10.1186/s12859-021-03983-2.
4. Smit AFA, Hubley R, Green P. *RepeatMasker Open-4* documentation. https://www.repeatmasker.org/webrepeatmaskerhelp.html.
5. Benson G. 1999. Tandem repeats finder. *Nucleic Acids Res* 27(2):573–580.
6. Rice P, Longden I, Bleasby A. 2000. EMBOSS. *Trends Genet* 16:276–277 — `palindrome` (M. Faller), `palindrome.acd` defaults (minpallen 10, maxpallen 100, gaplimit 100); EMBOSS 6.6.0 binary run locally.
7. Kolpakov R, Kucherov G. 1999. Finding maximal repetitions in a word in linear time. *FOCS 1999*, 596–604.
