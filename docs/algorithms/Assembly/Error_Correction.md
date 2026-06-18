# K-mer Spectrum Read Error Correction

| Field | Value |
|-------|-------|
| Algorithm Group | Assembly |
| Test Unit ID | ASSEMBLY-CORRECT-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

K-mer spectrum error correction repairs substitution errors in sequencing reads by exploiting the fact that, at sufficient depth, a correct k-mer recurs many times across reads while a k-mer carrying a sequencing error is rare. The algorithm classifies k-mers as *trusted* (multiplicity ≥ a coverage cut-off) or *untrusted*, then changes bases that are covered only by untrusted k-mers to the unique alternative base that makes their covering k-mers trusted [1][2]. It is a deterministic heuristic (not exact): it corrects single substitutions under the conservative assumption of at most one error per k-mer and leaves a base unchanged whenever the correction is ambiguous [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

High-throughput sequencing reads contain base-call errors (predominantly substitutions for Illumina). For a genome sequenced to coverage *c*, a k-mer from the true genome appears with multiplicity centred on a high value, whereas a k-mer produced by a base-call error appears rarely. Plotting k-mer multiplicities yields a histogram with a high-coverage mode (genomic k-mers) separated by a valley from a low-coverage mode (error k-mers); the valley defines the coverage cut-off [1][2].

### 2.2 Core Model

Let *S* be the multiset of all length-*k* substrings (k-mers) over all reads, and let `mult(x)` be the multiplicity of k-mer *x* in *S*. Given a cut-off *t*:

- A k-mer *x* is **trusted** (solid) iff `mult(x) ≥ t`, and **untrusted** (weak) otherwise [1][3].
- A read position *i* is **trusted** iff some k-mer covering *i* is trusted; bases inside trusted k-mers are assumed correct and are not modified [1][2].
- For an untrusted position *i*, the **two-sided correction** seeks a base *b* such that *every* k-mer covering *i* becomes trusted after the substitution. The correction is applied only if such *b* is **unique**; if zero or more than one candidate qualifies, the base is left unchanged (ambiguity) [1].

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
| INV-04 | A base changes only to the unique base making all covering k-mers trusted | Two-sided unique-alternative rule [1] |
| INV-05 | Deterministic output | Candidates tested in fixed A,C,G,T order; spectrum built once |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| reads | `IReadOnlyList<string>` | required | Reads to correct (compared case-insensitively) | non-null; null elements skipped in spectrum |
| kmerSize | `int` | 15 | k-mer length *k* | ≥ 1 |
| minKmerFrequency | `int` | 2 | Trusted cut-off *t*: k-mers with multiplicity ≥ t are trusted | ≥ 1 in practice |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `IReadOnlyList<string>` | Corrected reads, upper-cased, in input order; count and per-read length preserved |

### 3.3 Preconditions and Validation

Null `reads` → `ArgumentNullException`; `kmerSize < 1` → `ArgumentOutOfRangeException`. Reads are upper-cased internally (case-insensitive). A read shorter than *k* contributes no k-mers and is returned unchanged (upper-cased). The k-mer spectrum is built once from the original reads and is not updated during correction. Substitution candidates are the four DNA bases A, C, G, T; non-ACGT symbols are not used as replacements.

## 4. Algorithm

### 4.1 High-Level Steps

1. Build the k-mer spectrum: count every overlapping length-*k* substring across all reads (case-insensitive).
2. For each read, for each position *i*:
   - If *i* is covered by any trusted k-mer (multiplicity ≥ *t*), leave the base.
   - Otherwise, for each alternative base, substitute and test whether *all* k-mers covering *i* become trusted; count the qualifying alternatives.
   - Apply the substitution only if exactly one alternative qualifies; otherwise restore the original base.
3. Return the corrected reads in input order.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Trusted cut-off** *t* (`minKmerFrequency`): `mult(x) ≥ t ⇒ trusted` [1].
- **Candidate alphabet:** `{A, C, G, T}` in fixed order (determinism) [1][2].
- **Spectrum:** a hash map `k-mer → multiplicity` (ordinal string keys).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Build spectrum | O(N·L·k) | O(distinct k-mers · k) | N reads of length L; substring of length k per position |
| Correct one read | O(L · 4 · k · k) | O(k) | per position: 3 candidates × ≤k covering k-mers × O(k) k-mer read |
| Total | O(N·L·k²) | O(distinct k-mers · k) | matches Registry O(n·r·k) order (extra k from k-mer construction) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.ErrorCorrectReads(reads, kmerSize, minKmerFrequency)`: public entry point; builds the spectrum and corrects each read.
- `BuildKmerSpectrum`, `CorrectRead`, `IsPositionTrusted`, `AllCoveringKmersTrusted`, `CoveringKmerStarts`, `KmerAt`: private helpers realizing the two-sided rule.

### 5.2 Current Behavior

The spectrum is computed once from the original reads and held fixed during correction (a correction does not feed back into k-mer counts within the same pass), matching a single Musket correction stage. Positions are processed left to right; an applied correction at position *i* is visible to covering-k-mer tests for later positions in the same read. Non-ACGT residues are never produced (candidate set is ACGT). Search reuse: the repository suffix tree was **not** used — this is a frequency-table (hash-map) computation, not exact substring occurrence enumeration in one text; a suffix tree offers no advantage for counting overlapping k-mers across many short reads.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Trusted/untrusted classification by multiplicity cut-off [1][2].
- Trusted-base rule: a base covered by any trusted k-mer is not modified [1][2].
- Two-sided correction: change a base to the unique alternative making all covering k-mers trusted [1].
- Ambiguity rule: >1 valid alternative ⇒ leave unchanged [1].
- Single-base substitution model (no indels) [1][2].

**Intentionally simplified:**

- Coverage cut-off: parameterised instead of auto-selected from the histogram valley; **consequence:** the caller must supply a data-appropriate cut-off rather than relying on automatic selection [1][2].
- k selection: fixed default; **consequence:** caller should pass *k* suited to the data depth [2].

**Not implemented:**

- Quality-weighted (q-mer) counting; **users should rely on:** no current alternative in-repo [2].
- Read trimming/discarding of uncorrectable reads; **users should rely on:** `SequenceAssembler.QualityTrimReads` for quality-based trimming [2].
- Multi-stage iteration and aggressive/two-sided combined passes of full Musket; **users should rely on:** repeated calls if needed [1].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Cut-off and k are parameters, not auto-selected | Assumption | Caller must choose them from data | accepted | ASM-03; defaults documented as non-behavioral for tested contract |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null reads | `ArgumentNullException` | Contract validation |
| `kmerSize < 1` | `ArgumentOutOfRangeException` | k-mer length must be positive |
| Read shorter than *k* | Returned unchanged (upper-cased) | No k-mers cover any position |
| All reads identical/error-free | Unchanged | Every position covered by a trusted k-mer (INV-03) |
| Ambiguous position | Unchanged | >1 valid alternative (INV-04) |
| No valid alternative | Unchanged | No substitution makes covering k-mers trusted |

### 6.2 Limitations

Corrects substitutions only; indel errors are unsupported and may degrade nearby correction (ASM-01). Assumes ≤1 error per k-mer (ASM-02). Effectiveness depends on a cut-off that genuinely separates genomic from error k-mers; at low coverage or in repeats the cut-off mislabels k-mers (ASM-03). It is a single-pass heuristic, not an optimal maximum-likelihood corrector.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through (k=3, cut-off=2):** reads `ACGTACGT` ×3 (true) and `ACGTTCGT` ×1 (substitution A→T at index 4).

Spectrum: `ACG=7, CGT=8, GTA=3, TAC=3, GTT=1, TTC=1, TCG=1`. Trusted (≥2): ACG, CGT, GTA, TAC. In `ACGTTCGT`, only index 4 is covered solely by untrusted k-mers (GTT, TTC, TCG). Testing substitutions at index 4: `A` makes the covering windows GTA(3), TAC(3), ACG(7) — all trusted; `C`→GTC (absent) fails; `G`→GTG (absent) fails. The unique valid base is `A`, so index 4 is corrected, yielding `ACGTACGT`. All four outputs are `ACGTACGT`.

**API usage example:**

```csharp
var reads = new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTTCGT" };
IReadOnlyList<string> corrected =
    SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);
// corrected[3] == "ACGTACGT"
```

**Performance baseline (O(N·L·k²), above O(n²) — recorded here per Definition of Done):** measured in the test suite, the property test `ErrorCorrectReads_RandomReads_PreservesLengthAndTrustedPositions` corrects 60 reads of length 40 (k=7) and the full 13-test fixture completes in ~7 ms on the CI host. The dominant cost is the one-time spectrum build (O(N·L·k)); per-read correction adds O(L·k²) per read.

### 7.2 Applications and Use Cases

- **Pre-assembly cleanup:** running k-spectrum correction before OLC/de Bruijn assembly reduces spurious branches and tips caused by substitution errors, improving contig contiguity [1][2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAssembler_ErrorCorrectReads_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_ErrorCorrectReads_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [ASSEMBLY-CORRECT-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-CORRECT-001-Evidence.md)
- Related algorithms: [Quality_Trimming](Quality_Trimming.md), [De_Bruijn_Graph_Assembly](De_Bruijn_Graph_Assembly.md)

## 8. References

1. Liu Y, Schmidt B, Maskell DL. 2013. Musket: a multistage k-mer spectrum-based error corrector for Illumina sequence data. *Bioinformatics* 29(3):308-315. https://doi.org/10.1093/bioinformatics/bts690
2. Kelley DR, Schatz MC, Salzberg SL. 2010. Quake: quality-aware detection and correction of sequencing errors. *Genome Biology* 11:R116. https://doi.org/10.1186/gb-2010-11-11-r116
3. Song L, Florea L. 2018. Mining statistically-solid k-mers for accurate NGS error correction. PMC6311904. https://pmc.ncbi.nlm.nih.gov/articles/PMC6311904/
