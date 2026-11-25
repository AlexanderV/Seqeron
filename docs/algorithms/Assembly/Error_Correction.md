# K-mer Spectrum Read Error Correction

| Field | Value |
|-------|-------|
| Algorithm Group | Assembly |
| Test Unit ID | ASSEMBLY-CORRECT-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Complete (Musket 1.1 port; see §5.4 deviations) |
| Last Reviewed | 2026-11-25 |

## 1. Overview

K-mer spectrum error correction repairs substitution errors in sequencing reads by exploiting the fact that, at sufficient depth, a correct k-mer recurs many times across reads while a k-mer carrying a sequencing error is rare. The algorithm classifies k-mers as *trusted* (multiplicity ≥ a coverage cut-off) or *untrusted*, then changes bases that are covered only by untrusted k-mers [1][2]. Seqeron implements Musket [1] (port of musket-1.1 `ParaCorrect`): stage 1 is the conservative two-sided correction (unique alternative, ≤1 error per k-mer, ambiguity ⇒ unchanged); the opt-in full corrector adds Musket's aggressive one-sided stage and voting stage and counts canonical (both-strand) k-mers. It is a deterministic heuristic, not a maximum-likelihood corrector.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

High-throughput sequencing reads contain base-call errors (predominantly substitutions for Illumina). For a genome sequenced to coverage *c*, a k-mer from the true genome appears with multiplicity centred on a high value, whereas a k-mer produced by a base-call error appears rarely. Plotting k-mer multiplicities yields a histogram with a high-coverage mode (genomic k-mers) separated by a valley from a low-coverage mode (error k-mers); the valley defines the coverage cut-off [1][2].

### 2.2 Core Model

Let *S* be the multiset of all length-*k* substrings (k-mers) over all reads, and let `mult(x)` be the multiplicity of k-mer *x* in *S*. Given a cut-off *t*:

- A k-mer *x* is **trusted** (solid) iff `mult(x) ≥ t`, and **untrusted** (weak) otherwise [1][3].
- A read position *i* is **trusted** iff some k-mer covering *i* is trusted; bases inside trusted k-mers are assumed correct and are not modified [1][2].
- For an untrusted position *i*, the **two-sided correction** seeks a base *b* that makes the **leftmost and the rightmost** k-mer covering *i* trusted (Musket: "evaluates both the leftmost and the rightmost k-mers that cover position i"; musket-1.1 `paraCorrectCoreTwoSides`). The correction is applied only if such *b* is **unique**; if zero or more than one candidate qualifies, the base is left unchanged (ambiguity) [1]. Each iteration applies only the first such correction (left to right) and recomputes the trusted bases.
- With `CountReverseComplement`, `mult(x)` = occurrences of *x* plus occurrences of its reverse complement (Musket counts the canonical representative; a palindromic k-mer once per occurrence).

Corrections are single-base substitutions; the rule conservatively assumes at most one substitution error per k-mer [1]. Quake frames the same operation as a search for a set of single-base edits that makes all k-mers in the error region trusted [2].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Errors are substitutions (no indels) | Indel errors are not corrected and may shift k-mers, defeating correction [1][2] |
| ASM-02 | At most one error per k-mer | Multiple close errors in one k-mer may be uncorrectable [1] |
| ASM-03 | The cut-off separates genomic from error k-mers | Too-low/high cut-off mislabels k-mers (a solid k-mer may carry an error; a weak k-mer may be error-free) [2][3] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output read count = input read count | Each read maps to exactly one corrected read [1][2] |
| INV-02 | Each output read length = its input length | Only substitutions are applied (no indels) [1][2] |
| INV-03 | A position covered by a trusted k-mer is never modified | Trusted-base rule [1][2] |
| INV-04 | In the two-sided stage a base changes only to the unique base making its leftmost and rightmost covering k-mers trusted | Two-sided unique-alternative rule [1] |
| INV-05 | Deterministic output | Spectrum built once; Musket candidate order and libstdc++ region-heap order reproduced |

INV-03/INV-04 hold for the two-sided stage (the three-argument overload, `MaxErrorsPerKmer = 0`); Musket's voting stage (stage 3) may change a trusted base.

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| reads | `IReadOnlyList<string>` | required | Reads to correct (ASCII case-insensitive) | non-null; no null element |
| kmerSize | `int` | 21 | k-mer length *k* (Musket `option.cpp`: default `-k` 21) | ≥ 1 |
| minKmerFrequency | `int` | 3 | Trusted cut-off *t* (Musket `-minmulti`); Musket's default (`-minmulti -1`) estimates it from the histogram valley, `main.cpp` starts the valley scan at index 2 and increments once more, so the estimate is always ≥ 3 and 3 is the smallest value Musket itself can choose | ≥ 1 |
| options | `ErrorCorrectionOptions` | (3-arg overload: `TwoSided` = canonical spectrum, stage 1 only) | `CountReverseComplement` (Musket: true), `MaxIterations` (`-maxiter`, Musket 2, ≥ 1), `MaxErrorsPerKmer` (`-maxerr`, Musket 4, 0 = stage 1 only, ≥ 0); `ErrorCorrectionOptions.Musket` = Musket 1.1 defaults; `TwoSidedForwardStrand` = explicit forward-strand-only opt-out | non-null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `IReadOnlyList<string>` | Corrected reads, upper-cased, in input order; count and per-read length preserved |

### 3.3 Preconditions and Validation

Null `reads`/`options` → `ArgumentNullException`; null read element → `ArgumentException`; `kmerSize < 1`, `minKmerFrequency < 1` (a cut-off of 0 would trust absent k-mers), `MaxIterations < 1`, `MaxErrorsPerKmer < 0` → `ArgumentOutOfRangeException`. Reads are ASCII upper-cased. A read shorter than *k* contributes no k-mers and is returned unchanged (upper-cased; Musket skips it). The k-mer spectrum (Core `SequenceExtensions.CountKmersSpan`) is built once from the original reads and is not updated during correction. Substitution candidates are A, C, G, T; non-ACGT symbols are never produced.

## 4. Algorithm

### 4.1 High-Level Steps

1. Build the k-mer spectrum: count every overlapping length-*k* substring across all reads (Core `CountKmersSpan`); with `CountReverseComplement` merge each k-mer with its reverse complement (Core `TryGetReverseComplement`).
2. For each read (musket-1.1 `ParaCorrect`):
   - **Stage 1** (≤ `MaxIterations` iterations): if every base is trusted, stop. Otherwise scan left to right; at the first untrusted position *i* whose leftmost and rightmost covering k-mers have exactly one common trusted alternative, substitute it and start the next iteration; if none exists, leave stage 1.
   - **Stage 2** (for e = 1..`MaxErrorsPerKmer`, ≤ `MaxIterations` iterations each, regions alternately longest-first / shortest-first): extend every maximal run of trusted k-mers to the right and to the left. An untrusted k-mer is fixed at its new end base by the unique trusted alternative, else by the first alternative of strictly highest multiplicity whose following (preceding) `MaxErrorsPerKmer − e + 1` k-mers are trusted; stop at the first unfixable k-mer or at a trusted base; more than e corrections within k bases are discarded.
   - **Stage 3** (after each e): each untrusted k-mer votes for every single substitution that makes it trusted; if the maximum vote is ≥ 3, each base whose maximum is reached by exactly one alternative is replaced. An error-free read stops all stages.
3. Return the corrected reads in input order.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Trusted cut-off** *t* (`minKmerFrequency`): `mult(x) ≥ t ⇒ trusted` [1].
- **Candidate alphabet:** `{A, C, G, T}`, enumerated as Musket does: on the representative strand (`km < twin ? km : twin` under Musket's packed 2-bit `memcmp` order), cyclically after the original base. The order only decides stage-2 ties of equal multiplicity.
- **Spectrum:** a hash map `k-mer → multiplicity` (ordinal string keys, span lookups).
- **Stage-3 minimum vote:** 3 (musket-1.1 `fixErroneousBase(.., 3, ..)`).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Build spectrum | O(N·L·k) | O(distinct k-mers · k) | N reads of length L; substring of length k per position |
| Stage 1 (one iteration) | O(L · k) | O(L + k) | solids + two probes of ≤3 alternatives per untrusted position; ≤ L iterations to convergence |
| Stages 2–3 | O(L · k²) per iteration | O(L) | look-ahead ≤ `MaxErrorsPerKmer` k-mers; voting 3k alternatives per untrusted k-mer |
| Total (3-arg) | O(N·L·k + Σ iterations · L · k) | O(distinct k-mers · k) | |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.ErrorCorrectReads(reads, kmerSize, minKmerFrequency)`: two-sided stage to convergence over a canonical spectrum (= `ErrorCorrectionOptions.TwoSided`); `ErrorCorrectionOptions.TwoSidedForwardStrand` is the explicit forward-only opt-out.
- `SequenceAssembler.ErrorCorrectReads(reads, kmerSize, minKmerFrequency, ErrorCorrectionOptions)`: Musket corrector; `ErrorCorrectionOptions.Musket` = Musket 1.1 defaults.
- Private: `BuildKmerSpectrum` (Core `CountKmersSpan`), `MusketCorrector` (`TwoSided` = `paraCorrectCoreTwoSides`, `OneSided` = `paraCorrectCore`, `SelectAllMutations`, `Successor`, `Predecessor`, `VoteAndFix`), `StlHeap` (libstdc++ heap order).
- MCP: `error_correct_reads` (`musket`, `countReverseComplement`, `maxIterations`, `maxErrorsPerKmer`).

### 5.2 Current Behavior

The spectrum is computed once from the original reads and held fixed (as Musket). Stage 1 applies one correction per iteration and recomputes the trusted bases from the corrected read. Verified against a line-by-line Python transcription of musket-1.1 on 3000 random cases (Evidence). Search reuse: the repository suffix tree was **not** used — this is a frequency-table (hash-map) computation, not exact substring occurrence enumeration in one text; a suffix tree offers no advantage for counting overlapping k-mers across many short reads.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Trusted/untrusted classification by multiplicity cut-off [1][2].
- Trusted-base rule: a base covered by any trusted k-mer is not modified [1][2].
- Two-sided correction: change a base to the unique alternative making its leftmost and rightmost covering k-mers trusted, one correction per iteration [1].
- Ambiguity rule: >1 valid alternative ⇒ leave unchanged [1].
- Canonical k-mer counting, aggressive one-sided stage and voting stage (opt-in, `ErrorCorrectionOptions.Musket`) [1].
- Single-base substitution model (no indels) [1][2].

**Intentionally simplified:**

- Coverage cut-off: parameterised instead of auto-selected from the histogram valley; **consequence:** the caller must supply a data-appropriate cut-off rather than relying on automatic selection [1][2].
- k selection: fixed Musket default k = 21 (as Musket without `-k`); **consequence:** caller should pass *k* suited to the genome size and depth [2].

**Not implemented:**

- Quality-weighted (q-mer) counting; **users should rely on:** no current alternative in-repo [2].
- Read trimming/discarding of uncorrectable reads (Musket `-maxtrim`, default 0 = off); **users should rely on:** `SequenceAssembler.QualityTrimReads` for quality-based trimming [2].
- Musket multiple-k mode (`-multik`) and automatic cut-off estimation from the histogram valley.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Cut-off and k are parameters, not auto-selected | Assumption | Caller must choose them from data | accepted | ASM-03; defaults documented as non-behavioral for tested contract |
| 2 | 3-arg overload counts canonical k-mers as Musket (default flipped from forward-only, audit A9-1) | Resolved | forward-only counting miscorrects reads from the opposite strand (test `ErrorCorrectReads_ReverseStrandReads_*`) | resolved | explicit opt-out `ErrorCorrectionOptions.TwoSidedForwardStrand` / `CountReverseComplement = false` |
| 3 | Reads shorter than 2k − 1: stage 1 probes the true leftmost covering k-mer | Deviation from musket-1.1 code (follows the paper) | musket-1.1 advances its left k-mer from position 0 and probes a spliced k-mer | accepted | oracle `fix_quirks` |
| 4 | Stage-2 left look-ahead tests the k-mers really preceding the corrected k-mer | Deviation from musket-1.1 code | musket-1.1 `predecessor(seq, ipos - 1, ..)` skips base ipos − 1 | accepted | test `ErrorCorrectReads_AggressiveStage_LookAheadUsesTruePrecedingKmers` |
| 5 | Non-ACGT symbols kept distinct; cut-off 1 honoured; output fully upper-cased | Deviation | Musket encodes non-ACGT as A, drops singleton k-mers, keeps input case of uncorrected bases | accepted | |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null reads / options | `ArgumentNullException` | Contract validation |
| Null read element | `ArgumentException` | Contract validation |
| `kmerSize < 1`, `minKmerFrequency < 1`, `MaxIterations < 1`, `MaxErrorsPerKmer < 0` | `ArgumentOutOfRangeException` | |
| Read shorter than *k* | Returned unchanged (upper-cased) | No k-mers cover any position |
| All reads identical/error-free | Unchanged | Every position covered by a trusted k-mer (INV-03) |
| Ambiguous position | Unchanged | >1 valid alternative (INV-04) |
| No valid alternative | Unchanged | No substitution makes covering k-mers trusted |

### 6.2 Limitations

Corrects substitutions only; indel errors are unsupported and may degrade nearby correction (ASM-01). Stage 1 assumes ≤1 error per k-mer (ASM-02); stages 2–3 (opt-in) handle up to `MaxErrorsPerKmer`. Effectiveness depends on a cut-off that genuinely separates genomic from error k-mers; at low coverage or in repeats the cut-off mislabels k-mers (ASM-03). It is a greedy heuristic, not an optimal maximum-likelihood corrector (Quake's quality-aware likelihood search is not implemented).

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through (k=3, cut-off=2):** reads `ACGTACGT` ×3 (true) and `ACGTTCGT` ×1 (substitution A→T at index 4).

Spectrum: `ACG=7, CGT=8, GTA=3, TAC=3, GTT=1, TTC=1, TCG=1`. Trusted (≥2): ACG, CGT, GTA, TAC. In `ACGTTCGT`, only index 4 is covered solely by untrusted k-mers (GTT, TTC, TCG). Leftmost covering k-mer GTT@2 (index 2): `A`→GTA(3) trusted, GTC/GTG absent; rightmost TCG@4 (index 0): `A`→ACG(7) trusted, CCG/GCG absent. The unique common base is `A`, so index 4 is corrected, yielding `ACGTACGT`. All four outputs are `ACGTACGT`.

**API usage example:**

```csharp
var reads = new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTTCGT" };
IReadOnlyList<string> corrected =
    SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);
// corrected[3] == "ACGTACGT"
```

**Performance baseline (O(N·L·k²), above O(n²) — recorded here per Definition of Done):** measured in the test suite, the property test `ErrorCorrectReads_RandomReads_PreservesLengthAndTrustedPositions` corrects 60 reads of length 40 (k=7) and the fixture completes in milliseconds. The dominant cost is the one-time spectrum build (O(N·L·k)); per-read correction adds O(L·k²) per read.

### 7.2 Applications and Use Cases

- **Pre-assembly cleanup:** running k-spectrum correction before OLC/de Bruijn assembly reduces spurious branches and tips caused by substitution errors, improving contig contiguity [1][2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAssembler_ErrorCorrectReads_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAssembler_ErrorCorrectReads_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [ASSEMBLY-CORRECT-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-CORRECT-001-Evidence.md)
- Related algorithms: [Quality_Trimming](Quality_Trimming.md), [De_Bruijn_Graph_Assembly](De_Bruijn_Graph_Assembly.md)

## 8. References

1. Liu Y, Schmidt B, Maskell DL. 2013. Musket: a multistage k-mer spectrum-based error corrector for Illumina sequence data. *Bioinformatics* 29(3):308-315. https://doi.org/10.1093/bioinformatics/bts690
2. Kelley DR, Schatz MC, Salzberg SL. 2010. Quake: quality-aware detection and correction of sequencing errors. *Genome Biology* 11:R116. https://doi.org/10.1186/gb-2010-11-11-r116
3. Song L, Florea L. 2018. Mining statistically-solid k-mers for accurate NGS error correction. PMC6311904. https://pmc.ncbi.nlm.nih.gov/articles/PMC6311904/
