# RNA Inverted Repeats (Potential Stem Regions)

| Field | Value |
|-------|-------|
| Algorithm Group | RNA Secondary Structure |
| Test Unit ID | RNA-INVERT-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Finds inverted repeats in a nucleotide sequence: a left arm `W` followed by an intervening loop and a right arm that is the **reverse complement** of `W`. Inverted repeats are the sequence signature of potential RNA hairpin stems — the two complementary arms can base-pair antiparallel to form the stem, leaving the intervening region as the loop [1][3]. The implementation reports **perfect** (zero-mismatch), ungapped stems only; it is an exact combinatorial search over a bounded loop window, not the scored dynamic-programming variant of EMBOSS einverted [3].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

An inverted repeat (IR) is "a single stranded sequence of nucleotides followed downstream by its reverse complement"; the intervening sequence may be any length, including zero, in which case the composite is a palindrome [2]. In RNA, the two complementary arms fold back on each other to form the double-stranded stem of a hairpin, with the spacer becoming the single-stranded loop [3].

### 2.2 Core Model

An IR is a string of the form `W G W̄ᴿ`, where `W` is the left arm, `G` the gap/loop with `|G| ≥ 0`, and `W̄ᴿ` the right arm equal to the reverse complement of `W` [1]. A gapped IR is *within k mismatches* when it can be written `W G W̄ᴿ` with `δH(W, W̄ᴿ) ≤ k`, where `δH` is the Hamming distance [1]. This implementation realizes the **perfect** case `k = 0`: every position of the left arm pairs antiparallel with the corresponding position of the right arm under the strict Watson-Crick/IUPAC complement (A⟷U, C⟷G, etc.) [1][5]. Concretely, for an arm of length `m` with left start `s1` and right start `s2`, the pairing is `complement(seq[s2 + m - 1 - k]) == seq[s1 + k]` for all `k ∈ [0, m)`.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Arms pair perfectly (zero mismatches, no internal gaps) under strict Watson-Crick/IUPAC complement (no G-U wobble). | Imperfect/wobble stems are not reported; output is a strict subset of einverted's scored hits. |
| ASM-02 | The loop length `|G|` is bounded by `[minSpacing, maxSpacing]`. | Stems with loops outside the window are not reported. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Right arm = reverse complement of left arm: `complement(seq[Start2+Length-1-k]) == seq[Start1+k]` for all `k`. | Enforced by the outward antiparallel extension loop [1]. |
| INV-02 | Both arms have equal length: `End1-Start1+1 == Length == End2-Start2+1`. | Arms are extended in lockstep. |
| INV-03 | Loop length `Start2 - End1 - 1 ∈ [minSpacing, maxSpacing]`. | Loop boundary scan restricts `q ∈ [p+1+minSpacing, p+1+maxSpacing]`. |
| INV-04 | `Length ≥ minLength`. | Hits with `len < minLength` are discarded. |
| INV-05 | `Start1 ≤ End1 < Start2 ≤ End2` (arms disjoint, left precedes right). | `q > p` and arms cannot cross into the loop. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | Nucleotide sequence (DNA or RNA), case-insensitive | null/empty → empty result |
| minLength | int | 4 | Minimum arm length | `< 1` → empty result |
| minSpacing | int | 3 | Minimum loop length | `< 0` → empty result |
| maxSpacing | int | 100 | Maximum loop length | `< minSpacing` → empty result |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Start1, End1 | int | 0-based inclusive bounds of the left arm `W` |
| Start2, End2 | int | 0-based inclusive bounds of the right arm `W̄ᴿ` |
| Length | int | Common arm length (`End1-Start1+1`) |

Returned lazily as `IEnumerable<(int Start1, int End1, int Start2, int End2, int Length)>`.

### 3.3 Preconditions and Validation

0-based, inclusive coordinates. Input is upper-cased internally (case-insensitive). Accepts DNA or RNA; complement uses `GetRnaComplementBase` (A→U, U/T→A, G→C, C→G, IUPAC degenerate codes mapped, non-IUPAC passed through). No exceptions are thrown: null/empty input, out-of-range parameters, or sequences shorter than `2·minLength + minSpacing` yield an empty sequence. Reporting is non-overlapping: one maximal stem is computed per loop boundary, then repeats are selected longest-first and any repeat overlapping an already-selected arm is skipped.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; reject null/empty, bad parameters, and too-short sequences (empty result).
2. For each loop boundary defined by the last index of a candidate left arm `p` and the first index of a candidate right arm `q` with `q - p - 1 ∈ [minSpacing, maxSpacing]`:
   - Extend the stem **outward** from the innermost pair `(seq[p], seq[q])`: while `complement(seq[q+k]) == seq[p-k]`, grow the arm.
   - Keep the longest valid arm (≥ minLength) for this `p`.
3. Collect the best (maximal) stem per loop boundary, sort by length descending (ties by position), and emit non-overlapping repeats longest-first.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Complement table: `GetRnaComplementBase` (strict Watson-Crick + IUPAC) [5]. No scoring matrix is used — matching is exact (perfect arms). Selection: candidates are sorted by arm length descending, then by left start, then by right start; non-overlapping repeats are emitted in that order. Per loop boundary, among equal-length arms the smaller loop (smaller `q`) is kept by scanning `q` ascending and replacing only on strictly greater length.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindInvertedRepeats | O(n² · m) | O(1) extra | n = sequence length, m = max arm length (≤ maxSpacing-bounded). Loop-boundary pairs are O(n · maxSpacing); each extension is O(m). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.FindInvertedRepeats(string, int, int, int)`: returns maximal, non-overlapping perfect inverted repeats.

### 5.2 Current Behavior

The detector scans loop boundaries `(p, q)` and extends the stem outward, which directly encodes antiparallel reverse-complement pairing (the bug in the previous version computed the right-arm offset from `minLength` rather than the matched arm length, so it could not locate arms longer than `minLength`; this is corrected). Reporting is greedy and non-overlapping.

**Search reuse (suffix tree):** the repository `SuffixTree` was evaluated and **not** used. Inverted-repeat detection is antiparallel reverse-complement extension over a bounded loop window, not exact-substring occurrence enumeration; a suffix tree's strength (O(m) lookup of a fixed pattern after O(n) construction) does not map onto "for every center, extend complementary arms while bounding the loop". The established tool for this task, EMBOSS einverted, uses local dynamic programming rather than a suffix tree [3]. Correctness is unaffected by this choice — the formal `WGW̄ᴿ` output is identical.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The `W G W̄ᴿ` inverted-repeat pattern with the right arm = reverse complement of the left arm [1][2].
- Perfect (k = 0) antiparallel pairing under strict Watson-Crick/IUPAC complement [1][5].
- Bounded loop length `|G| ∈ [minSpacing, maxSpacing]` and minimum arm length [1].

**Intentionally simplified:**

- Mismatches and internal gaps: only `k = 0` perfect, ungapped arms are reported; **consequence:** imperfect stems (bulges, mismatches, G-U wobble) that EMBOSS einverted would score are not returned [1][3].
- Reporting is greedy/non-overlapping rather than enumerating every sub-arm; **consequence:** nested or overlapping shorter repeats inside a longer reported one are not separately listed.

**Not implemented:**

- Scored local-alignment IR detection (match/mismatch/gap penalties, threshold score); **users should rely on:** EMBOSS einverted for imperfect/scored inverted repeats [3].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Right-arm offset computed from `minLength` (previous code) | Deviation | Long arms mislocated / missed | fixed | Rewritten to extend from the matched arm length (ASM-01). |
| 2 | Perfect-arm-only restriction | Assumption | No imperfect stems | accepted | ASM-01; documented as Simplified. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty sequence | empty result, no throw | sibling-method convention (yield break) |
| sequence shorter than `2·minLength + minSpacing` | empty result | cannot hold two arms + loop |
| parallel direct repeat (left == right 5'→3') | empty result | right arm must be the *reverse* complement, not identical [1] |
| loop length outside `[minSpacing, maxSpacing]` | repeat not reported | INV-03 |
| arm shorter than minLength | repeat not reported | INV-04 |

### 6.2 Limitations

Detects perfect, ungapped stems only — no mismatches, bulges, or G-U wobble. Greedy non-overlapping reporting may omit alternative overlapping stems. Loop is bounded by `maxSpacing`; very long-range pairings beyond the window are not found. For imperfect/scored inverted repeats use EMBOSS einverted or IUPACpal.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// UUACG ... loop AAAAAA ... CGUAA  (CGUAA = reverse complement of UUACG)
var irs = RnaSecondaryStructure.FindInvertedRepeats("UUACGAAAAAACGUAA").ToList();
// irs[0] == (Start1: 0, End1: 4, Start2: 11, End2: 15, Length: 5)
```

**Numerical walk-through:** loop boundary `p = 4` (G), `q = 11` (C); loop length `11-4-1 = 6 ∈ [3,100]`. Extend outward: (G,C)✓ (C,G)✓ (A,U)✓ (U,A)✓ (U,A)✓ → arm length 5; left arm `[0,4]` = UUACG, right arm `[11,15]` = CGUAA = reverse complement of UUACG.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_FindInvertedRepeats_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructure_FindInvertedRepeats_Tests.cs) — covers INV-01..INV-05
- Evidence: [RNA-INVERT-001-Evidence.md](../../../docs/Evidence/RNA-INVERT-001-Evidence.md)

## 8. References

1. Alamro H, Alzamel M, Iliopoulos CS, Pissis SP, Watts S. 2021. IUPACpal: efficient identification of inverted repeats in IUPAC-encoded DNA sequences. BMC Bioinformatics 22:51. https://doi.org/10.1186/s12859-021-03983-2
2. Ussery DW, Wassenaar TM, Borini S. 2008. Computing for Comparative Microbial Genomics: Bioinformatics for Microbiologists. Springer. https://doi.org/10.1007/978-1-84800-255-5 (definition via Wikipedia "Inverted repeat": https://en.wikipedia.org/wiki/Inverted_repeat)
3. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite (einverted). Trends in Genetics 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2 (manual: https://emboss.bioinformatics.nl/cgi-bin/emboss/help/einverted)
5. Biopython Bio.Seq complement_rna (Watson-Crick/IUPAC complement basis). https://biopython.org/docs/latest/api/Bio.Seq.html
