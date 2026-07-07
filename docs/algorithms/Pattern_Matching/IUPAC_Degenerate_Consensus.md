# IUPAC-Degenerate Consensus Generation

| Field | Value |
|-------|-------|
| Algorithm Group | Pattern Matching / Matching |
| Test Unit ID | MOTIF-GENERATE-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Given a set of equal-length aligned DNA sequences, this algorithm produces a single consensus string in which each column is summarised by an IUPAC nucleotide symbol. Unlike a plain most-frequent ("plurality") consensus, ambiguous columns are encoded with IUPAC degeneracy codes (R, Y, B, N, …) so the consensus retains the set of bases that occur with appreciable frequency at that position [1][2]. A base is included in a column's code only if its frequency exceeds a fixed threshold; the surviving base set is then mapped to its IUPAC symbol [4]. The computation is exact and deterministic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A multiple alignment of related sequences can be collapsed into one representative string. When several bases coexist at a column, a single most-frequent base discards information; the IUPAC degenerate code instead names the *set* of bases present, which is the standard way to report position variability, mixed probes, and consensus motifs [1].

### 2.2 Core Model

For each column *j* of *n* aligned sequences, count the occurrences of each standard base. Retain the set `B_j = { b : count(b) > θ·n }`, where θ is a frequency threshold; remaining (low-frequency) bases are dropped [4]. The column emits `IUPAC(B_j)`, the single symbol that the NC-IUB 1984 nomenclature assigns to that base set [1]:

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
| sequences | `IEnumerable<string>` | required | Aligned DNA sequences | upper/lower-case A/C/G/T; column length taken from first sequence |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `string` | Consensus over the 15 IUPAC symbols; one symbol per column |

### 3.3 Preconditions and Validation

Null `sequences` throws `ArgumentNullException`. An empty collection returns `""`. Input is upper-cased before counting (case-insensitive). Only A/C/G/T are counted; other characters at a position are ignored. Indexing is 0-based per column; the column count equals the first sequence's length.

## 4. Algorithm

### 4.1 High-Level Steps

1. Upper-case all sequences; if none, return `""`.
2. For each column *j* (0 .. firstLength−1), tally counts of A, C, G, T.
3. Compute `threshold = n × 0.25`; keep bases whose count is strictly greater than the threshold.
4. If at least one base passes, map the surviving base set to its IUPAC symbol; otherwise emit the single most-frequent base (alphabetical tie-break).
5. Append the symbol; return the assembled string.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **IUPAC set→symbol table** — NC-IUB 1984 / Cornish-Bowden [1], corroborated by UCSC [2] and Wikipedia Table 1 [3]; realised as the `switch` in `GetIupacCode`.
- **Inclusion threshold** — `IupacInclusionThreshold = 0.25`; a base must occur in strictly more than a quarter of the sequences. The threshold-consensus mechanism is from DECIPHER `ConsensusSequence` [4]; the 0.25 value and strict `>` boundary are this implementation's documented design constant.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| GenerateConsensus | O(n × m) | O(m) | n sequences, m columns; constant 4-base tally and O(1) symbol lookup per column |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs)

- `MotifFinder.GenerateConsensus(IEnumerable<string>)`: builds the IUPAC-degenerate consensus.
- `MotifFinder.GetIupacCode(...)` (private): maps a column's passing base set to the IUPAC symbol via the threshold and NC-IUB table.

### 5.2 Current Behavior

Column length is taken from the first sequence; shorter sequences simply contribute fewer counts at trailing columns. Bases at exactly the threshold are excluded (strict `>`). When no base passes the threshold (e.g. four equally-frequent bases each at 25 %), the implementation falls back to the single most-frequent base, with ties resolved alphabetically — so a four-equal column yields `A`, never `N`. This is a single linear scan over the sequences per column; the repository suffix tree is not applicable (no substring search or occurrence enumeration is involved).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The NC-IUB 1984 set→symbol mapping for all 15 non-empty base sets (A,C,G,T,R,Y,S,W,K,M,B,D,H,V,N) [1][2][3].
- Threshold-then-encode mechanism: low-frequency bases removed, surviving set encoded with an IUPAC code [4].

**Intentionally simplified:**

- Threshold value: fixed at θ = 0.25 (strict `>`); **consequence:** users cannot tune the frequency cut (DECIPHER's default is 0.05); minority bases at ≤25 % are silently dropped.
- No-pass fallback emits the most-frequent base rather than a degeneracy code; **consequence:** a uniformly ambiguous column (four equal bases) returns a single base, not `N`.

**Not implemented:**

- Gap (`-`) handling and U/RNA columns; **users should rely on:** pre-normalising input to DNA without gaps.
- Configurable threshold / weighted consensus; **users should rely on:** an external tool (DECIPHER, EMBOSS `cons`) for parameterised thresholds.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | θ = 0.25 strict-`>` inclusion threshold | Assumption | Determines which minority bases enter the IUPAC code | accepted | Documented design constant; threshold-consensus family is authoritative [4]; symbol per passing set is fully source-backed [1] |
| 2 | No-pass fallback → most-frequent base (alphabetical tie) | Assumption | Four-equal column → single base not N | accepted | Implementation contract; no authoritative spec for this corner |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty collection | `""` | nothing to summarise |
| Null collection | `ArgumentNullException` | guard contract |
| Unanimous column | the standard base | singleton set → base [1] |
| Base at exactly 25 % | excluded | strict `>` boundary (INV-05) |
| Four equal bases (each 25 %) | most-frequent base (`A` on tie) | no base passes; fallback §5.2 |
| Lowercase input | same as upper-cased | case-insensitive normalisation |

### 6.2 Limitations

Gaps, IUPAC-degenerate input symbols, and RNA (U) are not counted. The threshold is fixed and not exposed; weighted or quality-aware consensus is out of scope. For four-way ambiguity the output is a single base rather than `N` due to the strict threshold and fallback.

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
