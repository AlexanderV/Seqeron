# Repetitive Element Detection and Classification

| Field | Value |
|-------|-------|
| Algorithm Group | Annotation |
| Test Unit ID | ANNOT-REPEAT-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Detects repetitive elements in a DNA sequence and assigns a repeat class. `FindRepetitiveElements` reports two structural classes: **tandem repeats** (a motif repeated head-to-tail in two or more directly adjacent copies) and **inverted repeats** (a sequence followed downstream by its reverse complement, the stem-loop form WGW^R). `ClassifyRepeat` assigns a RepeatMasker-style class label by matching a query against a supplied repeat-element library, falling back to a simple-repeat class by motif size. Detection is exact (no mismatches); classification is library-driven and therefore as complete as the library supplied.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Repetitive DNA makes up a large fraction of eukaryotic genomes. Tandem repeats (microsatellites/STRs with 1–6 bp motifs, minisatellites/VNTRs with ~10–60 bp motifs) arise by replication slippage; inverted repeats form hairpin/cruciform secondary structures and are recognition sites for many regulatory and restriction enzymes [1][2]. Interspersed repeats (SINE, LINE, LTR, DNA transposons) are classified by homology to curated libraries such as Repbase [4].

### 2.2 Core Model

**Tandem repeat.** "A pattern of one or more nucleotides is repeated and the repetitions are directly adjacent to each other" (head-to-tail), with a minimum of two copies [1]. A tandem array of unit `u` with `c` copies spans `c·|u|` bases. The reported unit is the **primitive** period of the array (the shortest `u` such that the array is `u` repeated); this avoids double-counting non-primitive units such as "AA" = "A"×2 [1].

**Inverted repeat.** An inverted repeat is "a single stranded sequence of nucleotides followed downstream by its reverse complement" [2]. Formally it has the form **WGW̄ᴿ**, where `W` is the left arm, `G` is the gap/loop/spacer with `|G| ≥ 0`, and `W̄ᴿ` is the reverse complement of `W` (the right arm) [3]. When `|G| = 0` the element is a reverse-complement palindrome [2][3].

**Repeat classification.** RepeatMasker assigns a class by the best homology match (Smith-Waterman-Gotoh, scored against a Repbase library) above a minimum score; classes are SINE, LINE, LTR, DNA, Satellite, Simple_repeat, Low_complexity, RNA, and Unclassified/Unknown [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A reported tandem repeat's `sequence` equals `dnaSequence[start..end]` (0-based, end-exclusive) and its length is an integer multiple of the unit length | reported span is `start..start+copies·|unit|`, sequence taken as that slice [1] |
| INV-02 | A reported tandem repeat has ≥ `minCopies` (≥ 2) directly adjacent identical primitive copies | copy count is incremented only on exact adjacent unit match; arrays below `minCopies` are dropped [1] |
| INV-03 | For every reported inverted repeat, the right arm equals the reverse complement of the left arm | a right arm is accepted only when it ordinal-equals `revcomp(leftArm)` [2][3] |
| INV-04 | `ClassifyRepeat` returns a label from the library or the fixed fallback set {Simple_repeat, Unknown} | return paths are library value, `Simple_repeat`, or `Unknown` [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| dnaSequence | string | required | DNA sequence to scan | non-null; case-insensitive (upper-cased internally) |
| minRepeatLength | int | 10 | minimum total span (bp) of a reported element | ≥ 1 |
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
2. **Tandem repeats:** for each primitive unit length `1..min(60, n/minCopies)` and each start, extend head-to-tail while the next unit matches; report each left-maximal array (≥ `minCopies`, span ≥ `minRepeatLength`) once, de-duplicated by span.
3. **Inverted repeats:** for each left-arm start and arm length (`minArm..min(100, remaining/2)`), compute `revcomp(arm)` and scan the gap window `[i+len, i+len+maxGap]` for a right arm equal to it; report each `(start,end)` once.
4. **Classification:** assign the class of the longest library element that the query contains or is contained by; else `Simple_repeat` if the query is a 1–6 bp primitive tandem repeat, else `Unknown`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Tandem unit length capped at 60 bp (covers STR 1–6 bp and minisatellite 10–60 bp) [1]; inverted arm capped at 100 bp; default max inverted gap 50 nt (hairpin loop range; matches the repo `RepeatFinder` default).
- Class vocabulary (Simple_repeat, Unknown, plus library labels SINE/LINE/LTR/DNA/Satellite) from RepeatMasker output classes [4].
- Primitivity test: a unit of length `n` is non-primitive iff it equals its own period `p` (for some divisor `p < n`) repeated `n/p` times.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindRepetitiveElements (tandem) | O(n²) | O(r) | n = length, r = reported arrays; unit-length and start loops |
| FindRepetitiveElements (inverted) | O(n² · L) | O(r) | L = max arm length; arm × start × gap window |
| ClassifyRepeat | O(d · (m + q)) | O(1) | d = library size, m = element length, q = query length (substring scan) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomeAnnotator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs)

- `GenomeAnnotator.FindRepetitiveElements(string, int, int)`: detects tandem and inverted repeats.
- `GenomeAnnotator.ClassifyRepeat(string, IReadOnlyDictionary<string,string>)`: assigns a repeat class via library match with simple-repeat / Unknown fallback.

### 5.2 Current Behavior

Detection is exact-match only (no mismatches/IUPAC degeneracy in arms). Tandem arrays are reported via their **primitive** unit and de-duplicated by span, so `AAAAAA` is reported once as motif `A`, not also as `AA`/`AAA`. Inverted repeats are reported once per `(start,end)` span; the composite `sequence` field embeds the gap in square brackets when the gap is non-empty.

**Suffix tree evaluation (search-reuse policy).** The repository [SuffixTree](../../../src/SuffixTree/Algorithms/SuffixTree/SuffixTree.cs) supports O(m) occurrence queries for a *known* pattern after O(n) construction. It was evaluated and **not used** for `FindRepetitiveElements`: repeat detection here is not a fixed-pattern occurrence query — it must *discover* unknown primitive units and check head-to-tail adjacency (tandem) or pair a left arm with the occurrence of its reverse complement within a bounded gap (inverted). These are structural-discovery problems, not "find occurrences of pattern X in text Y", so a direct suffix-tree lookup does not yield the required output without substantial additional machinery; a bounded O(n²) scan is the correct fit for the in-scope input sizes. For `ClassifyRepeat`, library matching is conceptually a homology (scored alignment) problem (RepeatMasker uses Smith-Waterman-Gotoh [4]); the implemented exact-substring containment is a simplification (§5.3) for which the suffix tree is also unnecessary at the small library/query sizes intended.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Tandem repeat = head-to-tail, ≥ 2 directly adjacent copies, primitive unit reported [1].
- Inverted repeat detected as WGW̄ᴿ with right arm = exact reverse complement of left arm, `|G| ≥ 0` (gap 0 ⇒ palindrome) [2][3].
- Repeat-class vocabulary (SINE/LINE/LTR/DNA/Satellite/Simple_repeat/Unknown) from RepeatMasker output classes [4].

**Intentionally simplified:**

- `ClassifyRepeat` library matching: exact substring containment instead of Smith-Waterman-Gotoh scoring; **consequence:** a query that is only *similar* (not identical/contained) to a library element is not matched and falls back to Simple_repeat/Unknown, whereas RepeatMasker would report it via an alignment score [4].
- Detection is perfect-match only; **consequence:** imperfect tandem arrays and imperfect (k-mismatch) inverted repeats defined by IUPACpal are not reported [3].

**Not implemented:**

- Curated Repbase/Dfam library, low-complexity masking, and interspersed-element family annotation; **users should rely on:** external tools (RepeatMasker) supplying a real library to `ClassifyRepeat`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Exact-substring library match vs. scored homology | Assumption | similar-but-not-identical repeats are not classified by library | accepted | §5.3; documented in Evidence Assumption 1 |

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
| ClassifyRepeat no library match, simple motif | `Simple_repeat` | STR 1–6 bp [1] |
| ClassifyRepeat no match, non-simple | `Unknown` | RepeatMasker Unclassified [4] |

### 6.2 Limitations

No mismatch/IUPAC tolerance; classification quality is bounded by the supplied library and the exact-substring relaxation; O(n²) scan is intended for short-to-moderate sequences, not whole chromosomes.

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
- Related algorithms: [Microsatellite detection](../Repeat_Analysis/Microsatellite_Detection.md)

## 8. References

1. Wikipedia contributors. 2026. *Tandem repeat*. Wikipedia. https://en.wikipedia.org/wiki/Tandem_repeat (cites Duitama J et al. 2014, *Nucleic Acids Res* 42(9):5728–5741, https://doi.org/10.1093/nar/gku212).
2. Wikipedia contributors. 2026. *Inverted repeat*. Wikipedia. https://en.wikipedia.org/wiki/Inverted_repeat (cites Ye C, Ji G, Liang C. 2014. *PLoS ONE* 9(11):e113349, https://doi.org/10.1371/journal.pone.0113349).
3. Hampson SE, Pissis SP, et al. 2021. *IUPACpal: efficient identification of inverted repeats in IUPAC-encoded DNA sequences*. BMC Bioinformatics 22:51. https://doi.org/10.1186/s12859-021-03983-2.
4. Smit AFA, Hubley R, Green P. *RepeatMasker Open-4* documentation. https://www.repeatmasker.org/webrepeatmaskerhelp.html.
