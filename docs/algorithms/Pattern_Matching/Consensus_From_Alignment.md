# Consensus Sequence from a Multiple Alignment

| Field | Value |
|-------|-------|
| Algorithm Group | Matching / Motif analysis |
| Test Unit ID | MOTIF-CONS-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-10-01 (F28: Auto type, ragged padding, `?` under -snucleotide) |

## 1. Overview

Given a set of equal-length aligned DNA sequences, this algorithm computes the consensus string: at each column it emits the most frequently occurring nucleotide [1][2]. It is the classical column-wise majority consensus used to summarise a multiple alignment as a single representative sequence. The computation is exact and deterministic; ties between equally-frequent bases are resolved in a fixed alphabetical order [4]. It differs from the IUPAC-degenerate consensus (`GenerateConsensus`), which encodes ambiguous columns with IUPAC codes rather than picking a single base.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A multiple sequence alignment places homologous sequences in a common coordinate frame so that each column contains positionally-corresponding residues. Summarising the alignment as one "typical" sequence — the consensus — is a standard step in motif description and primer/probe design [1].

### 2.2 Core Model

For aligned strings of length *n*, build a 4×*n* profile matrix *P* where *P*[b, j] is the number of times base *b* ∈ {A, C, G, T} occurs in column *j* [2]. The consensus *c* is defined position-wise: the *j*th symbol of *c* is "the symbol having the maximum value in the *j*-th column of the profile matrix" [2], i.e. the most frequent residue at that position [1]. When several symbols share the maximum count there may be more than one valid consensus [2]; this implementation fixes the choice to the alphabetically-earliest tied base [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output length = common input length | one symbol emitted per column [2] |
| INV-02 | Each output symbol is a base attaining the maximum count in its column | definition of consensus [2] |
| INV-03 | Ties resolve to the alphabetically-earliest base (A<C<G<T) | fixed scan order over the alphabet [4] |
| INV-04 | Identical inputs ⇒ output equals that sequence | every column is unanimous [1][2] |
| INV-05 | Deterministic for any valid input | INV-03 removes the only ambiguity [4] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | Most-frequent consensus (this) | IUPAC-degenerate consensus (`GenerateConsensus`) |
|--------|-------------------------------|--------------------------------------------------|
| Ambiguous column | single most-frequent base | IUPAC ambiguity code (e.g. R for A/G) |
| Output alphabet | A, C, G, T | A, C, G, T + IUPAC codes |

Two further alignment-consensus algorithms of reference tools are implemented alongside (B05 follow-up, §5.4):

| Aspect | EMBOSS `cons` (`GenerateEmbossConsensus`) | Biopython `dumb_consensus` (`GenerateDumbConsensus`) |
|--------|-------------------------------------------|------------------------------------------------------|
| Column score | substitution matrix (EDNAFULL / EBLOSUM62), weighted | raw residue counts |
| Gate | positive-match weight ≥ plurality (default half total weight); optional identity count | unique maximum with fraction ≥ threshold (default 0.7) |
| No consensus | `N` / `X`, lower case when positive matches ≤ setcase | the `ambiguous` symbol (default `X`) |
| Gaps / protein | yes / yes | yes / yes (any alphabet, case-sensitive) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| alignedSequences | `IEnumerable<string>` | required | Aligned DNA sequences | Equal length; alphabet {A,C,G,T}; case-insensitive |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `string` | Consensus; one uppercase base per alignment column. Empty string for an empty collection |

### 3.3 Preconditions and Validation

Null collection → `ArgumentNullException`. Null element → `ArgumentException`. Empty collection → `""`. Sequences of unequal length → `ArgumentException`. Any character outside {A,C,G,T} (after uppercasing) → `ArgumentException`. Input is uppercased before processing (case-insensitive); indexing is 0-based.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input; uppercase all sequences.
2. For each column, count occurrences of A, C, G, T (the profile column) [2].
3. Select the base with the maximum count, scanning A→C→G→T so ties resolve alphabetically [4].
4. Append the selected base to the consensus.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

Alphabet/order table: `{'A','C','G','T'}` — also the tie-break order [4]. No scoring matrix or plurality threshold is applied (contrast EMBOSS `cons`, which gates output on a weighted plurality value defaulting to half the total sequence weight [3]).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Consensus | O(n × m) | O(n) | n = column count, m = sequence count; O(1) extra per column for the 4-element profile |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs)

- `MotifFinder.CreateConsensusFromAlignment(IEnumerable<string>)`: column-wise most-frequent consensus with alphabetical tie-break.
- `MotifFinder.GenerateEmbossConsensus(IEnumerable<string>, ConsensusResidueType = Nucleotide, float? plurality = null, int identity = 0, float? setcase = null, IReadOnlyList<float>? weights = null)` ([MotifFinder.AlignmentConsensus.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.AlignmentConsensus.cs)): EMBOSS 6.6.0 `cons` [3][6].
- `MotifFinder.GenerateEmbossConsensus(IEnumerable<string>, bool padRaggedRows, ConsensusResidueType = Nucleotide, …)` (F28): the same with EMBOSS sequence-set padding — `padRaggedRows = true` appends `-` to rows shorter than the longest (`ajSeqsetFill`, applied by ACD to every `aligned: "Y"` input of `cons`); `ConsensusResidueType.Auto` (either overload) types the set like `cons` without `-snucleotide`/`-sprotein`.
- `MotifFinder.GenerateDumbConsensus(IEnumerable<string>, double threshold = 0.7, char ambiguous = 'X', bool requireMultiple = false)`: Biopython `SummaryInfo.dumb_consensus` [7].

### 5.2 Current Behavior

The column counts come from the private `BuildCountMatrix` (4 × L profile matrix, rows A, C, G, T; characters uppercased with `char.ToUpperInvariant`), the same helper `CreatePwm` uses, so the consensus and the PWM are derived from one count profile. The maximum is found with a strict `>` comparison while iterating rows in alphabetical order, so the first (alphabetically-earliest) maximum wins on a tie — the same scan as Biopython `Bio.motifs` `GenericPositionMatrix.consensus` (`if count > maximum` over the alphabet ACGT) [5]. **Search reuse:** the repository suffix tree was evaluated and is N/A — this is a column-wise tally over aligned positions, not a substring/occurrence search, so no pattern matching is involved.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Most-frequent residue per column = profile-matrix column maximum [1][2].
- Deterministic alphabetical tie-break (A<C<G<T) [4].

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none) — the EMBOSS `cons` weighted plurality consensus, formerly listed here, is implemented as `GenerateEmbossConsensus` (§5.4).

### 5.4 EMBOSS `cons` and Biopython `dumb_consensus` (B05 follow-up, 2026-09-30)

**EMBOSS `cons`** — line-by-line port of `embConsCalc` (EMBOSS 6.6.0 `nucleus/embcons.c`, driver `emboss/cons.c`, defaults `emboss/acd/cons.acd`) [6]. Per column, with weights w (default 1) and matrix M (EDNAFULL for nucleotides, EBLOSUM62 for proteins, both embedded verbatim from EMBOSS `data/`):

1. score(i) = Σ_{j≠i, both non-gap} M(rᵢ, rⱼ)·wⱼ (single precision, as in C);
2. candidate = first row of maximal score, except that a gap incumbent is displaced by a later row of equal score;
3. positive matches(r) = Σ wⱼ over rows j (r itself included) with M(r, rⱼ) > 0;
4. emit the candidate if positive matches ≥ plurality, else `N` (nucleotide) / `X` (protein); lower-case if positive matches ≤ setcase (this also lower-cases `N`/`X`); plurality and setcase default to half the total weight;
5. if identity > 0 and fewer than identity rows carry the residue with the most positive matches (ties → more identical weight), emit upper-case `N`/`X`.

Characters absent from the matrix (gaps, `*` in DNA, `J`/`O`/`U` in protein) have code 0: no score, no positive matches (they are only emitted when plurality ≤ 0). Input normalisation as the EMBOSS reader with an explicit type (`cons -snucleotide`/`-sprotein`): upper-casing, `.`/`~` → `-`, `X` → `N` for nucleotides (`ajSeqSetNuc`, applied before the `gapany` `?` → `X` conversion in `ajSeqTypeCheckIn`), so `?` is an unscored `X` (code 0) in both types; an `X` emitted into a nucleotide consensus is written as `N` (`cons.c` writes the result through `ajSeqSetNuc`). (F28 fix: `?` was previously scored as `N` in the explicit nucleotide type, which `cons -snucleotide` does not do — e.g. one-column `USC?TK`, `-plurality 4.51 -setcase 1.1` → `n`, was `N`.)

**Residue type `Auto` (F28)** — what `cons` does without a type flag (EMBOSS 6.6.0 `ajSeqsetFromList`, `ajSeqType`, `ajSeqIsNuc`/`ajSeqTypeGapnucS`, `ajSeqsetIsNuc`/`ajSeqsetIsProt`, `acdprotein`): every sequence is typed when read — nucleotide iff all its characters are in `ACGTU` + `BDHKMNRSVWXY?` + `.~-` (case-insensitive; empty → nucleotide), otherwise protein; a nucleotide-typed row reads `?` → `X` → `N` and `X` → `N`, a protein row keeps `X` (unscored in EDNAFULL). The *set* takes the first sequence's type, which selects the matrix (`$(acdprotein)`) and the no-consensus symbol (`ajSeqsetIsNuc`).

**Ragged rows (F28)** — `cons` never sees ragged rows: ACD pads every aligned sequence set with trailing gaps (`ajSeqsetFill`, `ajax/acd/ajacd.c`) before `embConsCalc`, so `cons` on `ACGTAC`, `ACG`, `AC` prints `ACGnnn`. `padRaggedRows = true` does the same; the original overload keeps rejecting unequal lengths.

**Remaining deliberate difference:** with the explicit `Protein` type the no-consensus symbol is `X`; `cons -sprotein` still decides `N` vs `X` (and writes emitted `X` as `N`) from the composition of the *first* sequence (`ajSeqsetIsNuc` tests the first sequence's characters even when the set type is `P`), so a protein alignment whose first row contains only nucleotide IUPAC letters prints `n` there where this type prints `x`. `Auto` reproduces `cons` without flags exactly.

**Cross-check (exact string equality) vs the EMBOSS 6.6.0 `cons` binary** (Ubuntu `emboss 6.6.0+dfsg-12ubuntu2`, `em_cons`): 780 alignments — 8 classic × 5 parameter sets + 700 seeded random (seeds 20260930 and 7; DNA and protein, 2–12 rows, 1–50 columns, gaps incl. `.`/`~`, lower case, IUPAC/B/Z/X, MSF weights 0.25–3 on ~30 %, random `-plurality`/`-identity`/`-setcase`): 780/780 identical, 14 of them after the documented first-sequence `N`/`X` difference. 110 are locked in `MotifFinder_AlignmentConsensus_Tests` with their command lines.

**F28 re-run (2026-10-01, Ubuntu noble `emboss` 6.6.0 `/usr/lib/emboss/cons`)**: `Auto` + `padRaggedRows` vs `cons` without type flags — 3,700/3,700 identical (seeds 20261001 ×1200 and 99 ×1000: DNA, IUPAC incl. `U`/`X`/`?`, protein ± `BZX*`, mixed rows, nucleotide-looking first row + protein rows, 50 % ragged, `-plurality` incl. 0 and negative, `-identity`, `-setcase`; seeds 20260930 ×800 and 4242 ×700 equal-length). Explicit types vs `cons -snucleotide`/`-sprotein` on the same 1,500 equal-length alignments: nucleotide 592/592 after the `?` fix (before it, 33/315 of seed 20260930 differed, all with `?`); protein 865/908 — the 43 differences are exactly the first-row-looks-nucleotide cases above (where `Auto` = `cons` without flags). 70 Auto cases + typing probes locked in `MotifFinder_EmbossConsensusAutoPad_Tests`.

**Biopython `dumb_consensus`** — port of `Bio/Align/AlignInfo.py` `SummaryInfo.dumb_consensus` (Biopython 1.85; deprecated since 1.82 and absent from the installed 1.88) [7]: per column count residues other than `-` and `.` (case-sensitive); emit the residue if it is the unique maximum and max/non-gap ≥ threshold, else `ambiguous`; with `requireMultiple` a column with exactly one non-gap residue is ambiguous. Cross-check vs Biopython 1.85: 708 alignments (4 classic + 704 seeded random, DNA/RNA/protein, gaps `-`/`.`, lower case, thresholds 0–1, require_multiple) → 708/708 identical; 34 locked.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty collection | `""` | nothing to summarise; mirrors `GenerateConsensus` |
| Single sequence | returns it (uppercased) | each column's only base is its maximum [2] |
| Tie column (A,G) | alphabetically-earliest (A) | tie-break rule [4] |
| Lowercase input | normalised to uppercase | case-insensitive contract |

### 6.2 Limitations

`CreateConsensusFromAlignment` operates on the DNA alphabet {A,C,G,T} only (no IUPAC ambiguity input, no gaps, no protein) — use `GenerateEmbossConsensus` or `GenerateDumbConsensus` for gapped / protein alignments. Requires pre-aligned equal-length input; it does not perform alignment. No confidence/plurality threshold is applied — every column yields a base.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var aligned = new[]
{
    "ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC",
    "TTGGAACT", "ATGCCATT", "ATGGCACT"
};
string consensus = MotifFinder.CreateConsensusFromAlignment(aligned); // "ATGCAACT"
```

**Numerical walk-through:** For the Rosalind CONS sample [2], the profile is A=`5 1 0 0 5 5 0 0`, C=`0 0 1 4 2 0 6 1`, G=`1 1 6 3 0 1 0 0`, T=`1 5 0 0 0 1 1 6`. Taking the column maxima (5→A, 5→T, 6→G, 4→C, 5→A, 5→A, 6→C, 6→T) gives `ATGCAACT`.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MotifFinder_CreateConsensusFromAlignment_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/MotifFinder_CreateConsensusFromAlignment_Tests.cs) — covers `INV-01`–`INV-05`
- Tests: [MotifFinder_AlignmentConsensus_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/MotifFinder_AlignmentConsensus_Tests.cs) (EMBOSS / Biopython-locked), `Properties/AlignmentConsensusProperties.cs`, `Metamorphic/AlignmentConsensusMetamorphicTests.cs`, `Fuzzing/AlignmentConsensusFuzzTests.cs`
- Evidence: [MOTIF-CONS-001-Evidence.md](../../../docs/Evidence/MOTIF-CONS-001-Evidence.md)

## 8. References

1. Wikipedia contributors. 2026. Consensus sequence. Wikipedia. https://en.wikipedia.org/wiki/Consensus_sequence
2. Rosalind. Consensus and Profile (CONS). https://rosalind.info/problems/cons/
3. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite. Trends in Genetics 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2 (program docs: https://www.bioinformatics.nl/cgi-bin/emboss/help/cons)
4. Los Alamos HIV Sequence Database. Advanced Consensus Maker — explanation. https://hfv.lanl.gov/content/sequence/CONSENSUS/AdvConExplain.html
5. Biopython 1.88, `Bio/motifs/matrix.py` — `GenericPositionMatrix.consensus` (installed package source; review 2026-09). Cock et al. 2009, Bioinformatics 25(11):1422.
6. EMBOSS 6.6.0 source: `nucleus/embcons.c` (`embConsCalc`, Tim Carver 2001), `emboss/cons.c`, `emboss/acd/cons.acd`, `ajax/core/ajseqtype.c`, `ajax/core/ajseq.c` (raw.githubusercontent.com/kimrutherford/EMBOSS, master); matrices `data/EDNAFULL`, `data/EBLOSUM62` (Ubuntu `emboss-data 6.6.0+dfsg-12ubuntu2`).
7. Biopython 1.85, `Bio/Align/AlignInfo.py` — `SummaryInfo.dumb_consensus` (raw.githubusercontent.com/biopython/biopython/biopython-185; PyPI wheel biopython==1.85 used as the reference).
