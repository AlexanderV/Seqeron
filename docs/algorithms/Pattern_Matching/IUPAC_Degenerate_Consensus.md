# IUPAC-Degenerate Consensus Generation

| Field | Value |
|-------|-------|
| Algorithm Group | Pattern Matching / Matching |
| Test Unit ID | MOTIF-GENERATE-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-09-30 (review-2026-09 B05, F13) |

## 1. Overview

Given a set of equal-length aligned DNA sequences, this algorithm produces a single consensus string in which each column is summarised by an IUPAC nucleotide symbol. Unlike a plain most-frequent ("plurality") consensus, ambiguous columns are encoded with IUPAC degeneracy codes (R, Y, B, N, …) so the consensus retains the set of bases that occur with appreciable frequency at that position [1][2]. A base is included in a column's code only if its frequency exceeds a fixed threshold; the surviving base set is then mapped to its IUPAC symbol. The library also provides the published Cavener (1987) rule set (`GenerateCavenerConsensus`, as in TRANSFAC and Biopython `degenerate_consensus`) [5][6]. The computation is exact and deterministic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A multiple alignment of related sequences can be collapsed into one representative string. When several bases coexist at a column, a single most-frequent base discards information; the IUPAC degenerate code instead names the *set* of bases present, which is the standard way to report position variability, mixed probes, and consensus motifs [1].

### 2.2 Core Model

For each column *j* of *n* aligned sequences, count the occurrences of each standard base. Retain the set `B_j = { b : count(b) > θ·n }`, where θ = 0.25 is this library's per-base design threshold; remaining (low-frequency) bases are dropped. This belongs to the threshold-consensus family but is **not** DECIPHER's rule (DECIPHER drops the least frequent characters while their *cumulative* fraction stays below `threshold`, default 0.05) [4]. If no base passes, the column emits the IUPAC symbol of all bases tied at the maximum count — DECIPHER: "degeneracy codes are always used in cases where multiple characters are equally abundant" [4] — and `N` when the column contains no A/C/G/T.

**Cavener (1987) rules** (`GenerateCavenerConsensus`) [5][6]: with counts sorted c1 ≥ c2 ≥ c3 ≥ c4 (ties in A,C,G,T order): single base if c1 > c2+c3+c4 and c1 > 2·c2; else two-base code of the top two if c1+c2 > 75 %; else three-base code of the top three if c4 = 0; else `N`. The column emits `IUPAC(B_j)`, the single symbol that the NC-IUB 1984 nomenclature assigns to that base set [1]:

| Base set | Symbol | Base set | Symbol |
|----------|--------|----------|--------|
| {A} | A | {A,G} | R |
| {C} | C | {C,T} | Y |
| {G} | G | {C,G} | S |
| {T} | T | {A,T} | W |
| {G,T} | K | {A,C} | M |
| {C,G,T} | B | {A,G,T} | D |
| {A,C,T} | H | {A,C,G} | V |
| {A,C,G,T} | N | | |

The mapping is bijective over the 15 non-empty subsets of {A,C,G,T} [1][2][3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output length = length of the first input sequence | one symbol emitted per column [4] |
| INV-02 | A unanimous column yields that standard base | singleton set → standard base [1] |
| INV-03 | Every output character is one of the 15 IUPAC symbols | image of the NC-IUB mapping [1][2] |
| INV-04 | The symbol for a passing base set equals NC-IUB's symbol for that set | bijective table [1][2][3] |
| INV-05 | A base with count ≤ θ·n is excluded (strict `>`) | threshold filtering precedes encoding [4]; θ = 0.25 here |

### 2.5 Comparison with Related Methods

| Aspect | IUPAC-degenerate consensus (this) | Most-frequent consensus (`CreateConsensusFromAlignment`) |
|--------|-----------------------------------|----------------------------------------------------------|
| Ambiguous column | IUPAC degeneracy code (R/Y/…/N) | single most-frequent base |
| Output alphabet | 15 IUPAC symbols | {A,C,G,T} |
| Threshold | frequency cut (θ = 0.25) | none (plurality) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequences | `IEnumerable<string>` | required | Aligned DNA sequences | equal length; upper/lower-case; non-ACGT characters not counted (still count in n) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `string` | Consensus over the 15 IUPAC symbols; one symbol per column |

### 3.3 Preconditions and Validation

Null `sequences` throws `ArgumentNullException`; a null element or rows of unequal length throw `ArgumentException` (F13; Biopython `MultipleSeqAlignment`: "Sequences must all be the same length"). An empty collection returns `""`. Input is upper-cased before counting (case-insensitive). Only A/C/G/T are counted; other characters at a position are ignored. Indexing is 0-based per column; the column count equals the first sequence's length.

## 4. Algorithm

### 4.1 High-Level Steps

1. Upper-case all sequences; if none, return `""`.
2. For each column *j* (0 .. firstLength−1), tally counts of A, C, G, T.
3. Compute `threshold = n × 0.25`; keep bases whose count is strictly greater than the threshold.
4. If at least one base passes, map the surviving base set to its IUPAC symbol; otherwise map the set of bases tied at the maximum count (`N` if no A/C/G/T at all).
5. Append the symbol; return the assembled string.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **IUPAC set→symbol table** — NC-IUB 1984 / Cornish-Bowden [1], corroborated by UCSC [2] and Wikipedia Table 1 [3]; realised by the canonical `IupacDnaSequence.GetIupacCode(IEnumerable<char>)` (Core/ISequence.cs), which `MotifFinder.GetIupacCode` calls after its threshold step.
- **Inclusion threshold** — `IupacInclusionThreshold = 0.25`; a base must occur in strictly more than a quarter of the sequences. Design constant (kept for API/MCP compatibility); not DECIPHER's cumulative rule [4] nor Cavener's [5].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| GenerateConsensus | O(n × m) | O(m) | n sequences, m columns; constant 4-base tally and O(1) symbol lookup per column |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs)

- `MotifFinder.GenerateConsensus(IEnumerable<string>)`: builds the IUPAC-degenerate consensus.
- `MotifFinder.GetIupacCode(...)` (private): applies the >25 % inclusion threshold, then maps the passing base set to its NC-IUB symbol via canonical `IupacDnaSequence.GetIupacCode`.
- `MotifFinder.GenerateCavenerConsensus(IEnumerable<string>)`: Cavener 1987 rules on the shared private `BuildCountMatrix` (rejects null/unequal/non-ACGT rows); set→symbol via `IupacDnaSequence.GetIupacCode`. Locked against Biopython 1.88 `degenerate_consensus` (tutorial WACVC/GBGTW/CV + 12 random alignments).

### 5.2 Current Behavior

All rows must have the same length. Bases at exactly the threshold are excluded (strict `>`). When no base passes the threshold (e.g. four equally-frequent bases each at 25 %, or a gap-diluted column), the column is the IUPAC code of the bases tied at the maximum count — a four-equal column yields `N` (was `A` before F13, 2026-09); a column with no A/C/G/T yields `N` (was `A`). This is a single linear scan over the sequences per column; the repository suffix tree is not applicable (no substring search or occurrence enumeration is involved).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The NC-IUB 1984 set→symbol mapping for all 15 non-empty base sets (A,C,G,T,R,Y,S,W,K,M,B,D,H,V,N) [1][2][3].
- Threshold-then-encode mechanism: low-frequency bases removed, surviving set encoded with an IUPAC code [4].

**Intentionally simplified:**

- Threshold value: fixed at θ = 0.25 (strict `>`); **consequence:** users cannot tune the frequency cut (DECIPHER's default is 0.05); minority bases at ≤25 % are silently dropped.
- Empty-information column (only gaps/N) → `N` (Biopython's Cavener code would give `V` for an all-zero column, an artefact of its sort; DECIPHER would give `-`).

**Not implemented:**

- Gap (`-`) handling and U/RNA columns; **users should rely on:** pre-normalising input to DNA without gaps.
- Configurable threshold / weighted consensus; **users should rely on:** `GenerateCavenerConsensus` for the published rule, or an external tool (DECIPHER, EMBOSS `cons`) for parameterised thresholds.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | θ = 0.25 strict-`>` inclusion threshold | Assumption | Determines which minority bases enter the IUPAC code | accepted | Documented design constant; threshold-consensus family is authoritative [4]; symbol per passing set is fully source-backed [1] |
| 2 | No-pass fallback → IUPAC code of the bases tied at the maximum; no A/C/G/T → N | Sourced (F13) | Four-equal column → N | fixed 2026-09 | DECIPHER equal-abundance rule [4]; Biopython `degenerate_consensus` = N [6] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty collection | `""` | nothing to summarise |
| Null collection | `ArgumentNullException` | guard contract |
| Unanimous column | the standard base | singleton set → base [1] |
| Base at exactly 25 % | excluded | strict `>` boundary (INV-05) |
| Four equal bases (each 25 %) | `N` | equally abundant → degeneracy code [4][6] |
| Column with only gaps / N | `N` | no information |
| Null element / unequal lengths | `ArgumentException` | not an alignment |
| Lowercase input | same as upper-cased | case-insensitive normalisation |

### 6.2 Limitations

Gaps, IUPAC-degenerate input symbols, and RNA (U) are not counted. The threshold is fixed and not exposed; weighted or quality-aware consensus is out of scope. The 25 % per-base threshold is a design constant, not a published rule; use `GenerateCavenerConsensus` for the Cavener/TRANSFAC/Biopython result.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Column 0 has {A,G}; columns 1-3 unanimous → "RTGC"
string consensus = MotifFinder.GenerateConsensus(new[] { "ATGC", "GTGC" });
```

**Numerical walk-through:** for `["C","G","T"]` (n=3): threshold = 3×0.25 = 0.75; each base count = 1 > 0.75, so the set is {C,G,T} → IUPAC symbol `B`.

### 7.2 Applications and Use Cases

- **Motif / binding-site summarisation:** an IUPAC consensus (e.g. `TGASTCA` for an AP-1-like site) captures positional ambiguity that a single-base consensus would lose [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MotifFinder_GenerateConsensus_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/MotifFinder_GenerateConsensus_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [MOTIF-GENERATE-001-Evidence.md](../../../docs/Evidence/MOTIF-GENERATE-001-Evidence.md)
- Related algorithms: [Consensus_From_Alignment](./Consensus_From_Alignment.md)

## 8. References

1. Cornish-Bowden A. 1985. Nomenclature for incompletely specified bases in nucleic acid sequences: recommendations 1984. Nucleic Acids Research 13(9):3021–3030. https://doi.org/10.1093/nar/13.9.3021
2. UCSC Genome Browser. IUPAC ambiguity codes. https://genome.ucsc.edu/goldenPath/help/iupac.html
3. Wikipedia. Nucleic acid notation (Table 1, citing NC-IUB 1984). https://en.wikipedia.org/wiki/Nucleic_acid_notation
4. Wright E.S. DECIPHER `ConsensusSequence` (Bioconductor). https://rdrr.io/bioc/DECIPHER/man/ConsensusSequence.html
5. Cavener D.R. 1987. Comparison of the consensus sequence flanking translational start sites in Drosophila and vertebrates. Nucleic Acids Research 15(4):1353–1361.
6. Biopython 1.88 `Bio.motifs.matrix.GenericPositionMatrix.degenerate_consensus` (installed source) and Tutorial `chapter_motifs.rst` (WACVC / GBGTW / CV examples). https://raw.githubusercontent.com/biopython/biopython/master/Doc/Tutorial/chapter_motifs.rst
