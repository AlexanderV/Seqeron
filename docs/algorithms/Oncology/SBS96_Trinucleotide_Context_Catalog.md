# SBS-96 Trinucleotide Context Catalog

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Mutational Signatures |
| Test Unit ID | ONCO-SIG-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

The SBS-96 trinucleotide-context catalog is the foundational data structure of single-base-substitution (SBS)
mutational-signature analysis. Each somatic single-base substitution is classified by its substitution type
together with the bases immediately 5' and 3' of the mutated base, yielding 96 categories (channels). Because
a substitution and its complement on the opposite strand are biologically equivalent, every mutation is folded
onto the pyrimidine strand — the substitution is referred to by the pyrimidine (C or T) of the mutated
Watson-Crick base pair [1][2][3]. The result is an exact, specification-driven classification: each variant
maps deterministically to exactly one of the 96 channels, and a multiset of variants is tallied into the
96-channel spectrum that downstream signature decomposition (NMF / fitting) consumes.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Somatic mutations in a cancer genome are the cumulative imprint of the mutational processes that operated in
the cell lineage. Many processes leave a characteristic bias in *which* base changes occur and in *what
sequence context*. Alexandrov et al. (2013) analysed 4,938,362 mutations from 7,042 cancers and showed that
representing each substitution by its trinucleotide context (the mutated base plus its immediate neighbours)
exposes these process-specific patterns [1]. The 96-channel representation is the standard substrate for the
COSMIC reference signatures [2].

### 2.2 Core Model

A single-base substitution at a reference position is described by the reference (mutated) base, the alternate
base, and the two flanking reference bases. The classification has three factors [1][2][3]:

- **Substitution type**, referred to by the pyrimidine of the mutated Watson-Crick base pair, so there are six:
  **C>A, C>G, C>T, T>A, T>C, T>G** [2].
- **5' base** ∈ {A, C, G, T}.
- **3' base** ∈ {A, C, G, T}.

This gives 6 × 4 × 4 = **96** mutation types [1][2]. The mutated base is centred in the trinucleotide; a channel
is written `5'[REF>ALT]3'` (e.g. `A[C>A]A`).

**Pyrimidine-strand folding.** When the reference (mutated) base is a purine (A or G), the mutation is not one
of the six pyrimidine substitutions and must be folded onto the pyrimidine strand by taking the reverse
complement of both the trinucleotide context and the substitution [3]. Using complement map A↔T, C↔G [4]:

```
fold(5', REF, ALT, 3')  when REF ∈ {A,G}:
    REF' = complement(REF)        # purine → pyrimidine
    ALT' = complement(ALT)
    5'' = complement(3')          # reverse + complement of the context
    3'' = complement(5')
    channel = 5''[REF'>ALT']3''
```

For example a G>T at context 5'-T G A-3' folds to `T[C>A]A` (reverse-complement of TGA is TCA; G>T becomes
C>A) [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every channel's reference base is a pyrimidine (C or T) | purine-reference inputs are reverse-complemented before counting [3] |
| INV-02 | The channel space is exactly 96 distinct labels | 6 substitutions × 4 5'-bases × 4 3'-bases [1][2] |
| INV-03 | Σ catalog counts = number of classifiable input variants | the catalog is a partition; each variant increments exactly one channel |
| INV-04 | Folding is reverse-complement: a purine-ref variant and its pyrimidine-strand form map to the same channel | strand-equivalence of complementary substitutions [3] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | SBS-96 | SBS-6 |
|--------|--------|-------|
| Factors | substitution × 5' × 3' | substitution only |
| Channels | 96 | 6 |
| Context sensitivity | yes (trinucleotide) | none |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| fivePrime | char | required | reference base 5' of the mutation | A/C/G/T, case-insensitive |
| referenceBase | char | required | mutated (reference) base | A/C/G/T, case-insensitive |
| alternateBase | char | required | substituted base | A/C/G/T, case-insensitive, ≠ reference |
| threePrime | char | required | reference base 3' of the mutation | A/C/G/T, case-insensitive |
| variants | IEnumerable<(char,char,char,char)> | required | SBS variants (5', ref, alt, 3') | each a valid substitution |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| ClassifySbsContext | string | channel label `5'[REF>ALT]3'` with a pyrimidine reference base |
| EnumerateSbs96Channels | IReadOnlyList<string> | the 96 canonical channel labels, deterministic order |
| Build96ContextCatalog | IReadOnlyDictionary<string,int> | all 96 channels → variant count (zeros included) |

### 3.3 Preconditions and Validation

Bases are upper-cased and must be A/C/G/T; any other character raises `ArgumentException`. For a single
substitution, reference == alternate raises `ArgumentException` (not a mutation). `Build96ContextCatalog(null)`
raises `ArgumentNullException`. There is no coordinate/indexing convention beyond the centred-mutation
trinucleotide; folding to the pyrimidine strand is automatic. Only single-base substitutions are in scope —
indels, doublet (DBS), and multi-base substitutions belong to other catalogues and are not handled here.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate and upper-case the four bases; reject ref == alt.
2. If the reference base is a purine (A/G): reverse-complement the context and the substitution onto the
   pyrimidine strand.
3. Emit the channel label `5'[REF>ALT]3'`.
4. For a catalog: initialise all 96 channels to 0, classify each variant, increment its channel.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

- Six pyrimidine substitutions: C>A, C>G, C>T, T>A, T>C, T>G [2].
- Complement map: A↔T, C↔G [4].
- Channel enumeration order: substitution-major, then 5' (A,C,G,T), then 3' (A,C,G,T) — a deterministic
  presentation convention; per-variant classification is independent of it.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| ClassifySbsContext | O(1) | O(1) | fixed-size base operations |
| EnumerateSbs96Channels | O(96) | O(96) | constant |
| Build96ContextCatalog | O(n) | O(1) extra (96-entry dict) | n = number of variants |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.ClassifySbsContext(char, char, char, char)`: folds one SBS to its 96-channel label.
- `OncologyAnalyzer.EnumerateSbs96Channels()`: returns the 96 canonical channel labels.
- `OncologyAnalyzer.Build96ContextCatalog(...)`: tallies SBS variants into the 96-channel spectrum.

### 5.2 Current Behavior

All 96 channels are always present in the catalog (zero-count channels included) so the spectrum has a fixed
shape suitable as a 96-dimensional vector for downstream decomposition. The dictionary is keyed by `Ordinal`
string comparison. No substring search is performed (classification is a constant-time base computation), so
the repository suffix tree is not applicable to this unit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Six pyrimidine substitution subtypes C>A, C>G, C>T, T>A, T>C, T>G [2].
- 96 = 6 × 4 × 4 trinucleotide channels with the mutated base centred [1][2].
- Reverse-complement folding of purine-reference mutations onto the pyrimidine strand [3].
- Watson-Crick complement map A↔T, C↔G [4].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Signature reference profiles / decomposition (NMF, COSMIC fitting): out of scope; **users should rely on**
  caller-supplied signature matrices and ONCO-SIG-002..004 (not yet implemented).
- Genome trinucleotide-frequency normalisation of the spectrum: out of scope; **users should rely on** a
  separate normalisation step when comparing to COSMIC profiles.
- DBS / ID / SBS-1536 catalogues: out of scope; **users should rely on** dedicated catalogues.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Purine reference base (A/G) | reverse-complemented to pyrimidine channel | strand equivalence [3] |
| Empty variant collection | all 96 channels present, each count 0 | partition of the empty set |
| Lower-case bases | classified identically to upper-case | bases upper-cased on input |
| ref == alt | ArgumentException | not a substitution |
| Non-ACGT base | ArgumentException | no defined trinucleotide context |

### 6.2 Limitations

Single-base substitutions only; the caller must supply the flanking reference bases (this unit does not read a
reference genome — extracting the trinucleotide context from a genome is the responsibility of ONCO-SIG-002).
The catalog is unnormalised counts, not genome-frequency-normalised proportions, so it should be normalised
before direct comparison against COSMIC reference profiles.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
// A G>T at 5'-T G A-3' folds to the pyrimidine strand: T[C>A]A
string channel = OncologyAnalyzer.ClassifySbsContext('T', 'G', 'T', 'A'); // "T[C>A]A"

var catalog = OncologyAnalyzer.Build96ContextCatalog(new[]
{
    ('A', 'C', 'A', 'A'),  // A[C>A]A  (pyrimidine, unchanged)
    ('T', 'G', 'T', 'A'),  // folds to T[C>A]A
});
// catalog["A[C>A]A"] == 1, catalog["T[C>A]A"] == 1, all other 94 channels == 0
```

**Numerical / biological walk-through:**

Fold G>T at 5'-T G A-3': complement each base of TGA → ACT, reverse → TCA; centre G→C, alt T→A ⇒ `T[C>A]A`.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_ClassifySbsContext_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_ClassifySbsContext_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-SIG-001-Evidence.md](../../../docs/Evidence/ONCO-SIG-001-Evidence.md)

## 8. References

1. Alexandrov, L.B., Nik-Zainal, S., Wedge, D.C., et al. 2013. Signatures of mutational processes in human cancer. Nature 500(7463):415-421. https://www.nature.com/articles/nature12477 (DOI: 10.1038/nature12477)
2. COSMIC Mutational Signatures — SBS96. Wellcome Sanger Institute. https://cancer.sanger.ac.uk/signatures/sbs/sbs96/
3. Bergstrom, E.N., Huang, M.N., Mahto, U., et al. 2019. SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. BMC Genomics 20:685. https://pmc.ncbi.nlm.nih.gov/articles/PMC6717374/ (DOI: 10.1186/s12864-019-6041-2)
4. Complementarity (molecular biology). Wikipedia. https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)
