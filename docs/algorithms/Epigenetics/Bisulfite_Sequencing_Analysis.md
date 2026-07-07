# Bisulfite Sequencing Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Epigenetics |
| Test Unit ID | EPIGEN-BISULF-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Bisulfite sequencing analysis turns the chemistry of sodium-bisulfite treatment into
base-level DNA methylation calls. Sodium bisulfite converts unmethylated cytosine to
uracil (read as thymine) while leaving 5-methylcytosine unchanged, so methylation becomes
visible as a C↔T difference between a read and the reference [1]. This unit provides three
specification-driven operations: an in-silico conversion simulator, a per-CpG methylation
caller from aligned reads, and a profile aggregator that summarizes per-context methylation.
The operations are deterministic and exact (no probabilistic modeling): each output is a
direct function of the input bases and the cited formulas.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

DNA methylation in animals occurs predominantly at the 5-position of cytosine in the CpG
dinucleotide; plants additionally methylate CHG and CHH contexts (H = A, C, or T). Frommer
et al. (1992) introduced bisulfite genomic sequencing, which "yields a positive display of
5-methylcytosine residues": after treatment, a remaining cytosine marks a methylated
position, and a thymine marks an originally-unmethylated cytosine [1].

### 2.2 Core Model

**Conversion (Frommer 1992) [1].** For each base b at position i of a single strand:
- b is an unprotected cytosine → thymine (uracil, read as T);
- b is a protected (5-methyl) cytosine → cytosine (unchanged);
- otherwise → b unchanged.

**Methylation call (Krueger & Andrews 2011, Bismark) [2].** At a reference CpG cytosine,
a read base C = methylated, a read base T = unmethylated; the methylation level (fraction)
is the Bismark percentage expressed in [0,1]:

```
level = methylated / (methylated + unmethylated)      [3]
```

**Per-context profile (Schultz 2012) [4].** For a set of cytosines of one context, the
weighted methylation level is the read-pooled fraction:

```
weighted level = Σ(methylated reads) / Σ(methylated + unmethylated reads)
              = Σ(levelᵢ · coverageᵢ) / Σ(coverageᵢ)                       [4]
```

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Conversion output length = input length | Per-base substitution only [1] |
| INV-02 | Protected cytosines stay C; unprotected C→T (case preserved) | Frommer conversion rule [1] |
| INV-03 | Non-cytosine bases are unchanged | Only cytosine reacts with bisulfite [1] |
| INV-04 | Each reported level ∈ [0,1], coverage ≥ 1 | level = meth/(meth+unmeth), denominator > 0 [3] |
| INV-05 | Per-context profile level = Σ(level·coverage)/Σ(coverage) | Weighted methylation definition [4] |
| INV-06 | Zero-coverage CpG sites are excluded from calling output | Percentage undefined when denominator = 0 [3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | Strand to convert | Null/empty → `""`; case preserved |
| `methylatedPositions` | `IReadOnlySet<int>?` | `null` | 0-based protected cytosine indices | Null = none protected |
| `referenceSequence` | `string` | required | Reference for calling | Null/empty → no sites; uppercased internally |
| `bisulfiteReads` | `IEnumerable<(string ReadSequence, int StartPosition)>` | required | Aligned reads | 0-based start; bases outside reference ignored |
| `sites` | `IEnumerable<MethylationSite>` | required | Measured sites for the profile | Empty → zero profile |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (conversion) | `string` | Converted strand, same length as input |
| `MethylationSite` | record | `Position` (0-based CpG C), `Type` = CpG, `Context`, `MethylationLevel` ∈ [0,1], `Coverage` |
| `MethylationProfile` | record | Global/CpG/CHG/CHH weighted levels, total/methylated CpG counts, per-position levels |

### 3.3 Preconditions and Validation

Indexing is 0-based. Sequences are case-insensitive for context/calling (uppercased
internally); the conversion simulator preserves case. Accepted alphabet is DNA {A,C,G,T};
a CpG requires both the C and the following G (the last reference base cannot start a CpG).
A read base at a reference-C position that is neither C nor T is not a valid bisulfite call
and is ignored. Null or empty `sequence`/`referenceSequence` produce empty results; an empty
`sites` enumerable produces an all-zero `MethylationProfile`.

## 4. Algorithm

### 4.1 High-Level Steps

1. **Conversion:** iterate bases; emit T (or t) for unprotected cytosines, keep protected
   cytosines and all non-cytosines.
2. **Calling:** locate reference CpG sites; for each aligned read, at every covered CpG
   tally C (methylated) and T (unmethylated); emit level = meth/(meth+unmeth) and coverage
   for sites with coverage > 0.
3. **Profile:** split sites by context; compute the read-pooled weighted level per context
   and globally.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Methylation call symbols (Bismark) [3]: `z/Z` CpG, `x/X` CHG, `h/H` CHH; lowercase =
  unmethylated, uppercase = methylated. This unit reports the numeric level rather than the
  symbol but uses the same C=methylated / T=unmethylated rule.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Conversion | O(n) | O(n) | One pass over the strand |
| Calling | O(R·L + n) | O(s) | R reads of length L, n reference length, s CpG sites |
| Profile | O(k) | O(k) | k input sites |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [EpigeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs)

- `EpigeneticsAnalyzer.SimulateBisulfiteConversion(sequence, methylatedPositions)`: in-silico conversion of one strand.
- `EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(referenceSequence, bisulfiteReads)`: per-CpG calling.
- `EpigeneticsAnalyzer.GenerateMethylationProfile(sites)`: per-context weighted aggregation.

### 5.2 Current Behavior

Conversion preserves input case and converts only the supplied strand. Calling reuses
`FindCpGSites` to enumerate reference CpG positions, then counts C/T read bases at each;
sites with no coverage are omitted. The profile uses the weighted methylation level and
falls back to an unweighted mean of per-site levels only when total coverage is zero
(e.g. sequence-only sites), so such sites are not silently dropped.

**Search reuse:** the only "find occurrences" step is locating CpG dinucleotides, a single
linear scan over the reference (`FindCpGSites`). The repository suffix tree is not used: a
suffix tree pays off for many queries of varying patterns against one text, whereas this is
one fixed two-character pattern in one pass — O(n) scanning is already optimal and simpler.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Unmethylated C → T, 5mC stays C, non-C unchanged (Frommer 1992) [1].
- Read C = methylated, read T = unmethylated at reference CpG (Krueger & Andrews 2011) [2].
- level = methylated / (methylated + unmethylated) (Bismark) [3].
- Per-context weighted methylation level Σ(meth)/Σ(meth+unmeth) (Schultz 2012) [4].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Complementary-strand conversion/merging; bisulfite non-conversion rate correction; CHG/CHH
  calling from reads (the caller targets CpG sites); **users should rely on:** dedicated
  pipelines (Bismark, methylKit) for full BS-seq processing.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Registry signature `CalculateMethylationFromBisulfite(bsSeq, refSeq)` | Deviation | Per-site coverage needs per-read input | accepted | Implemented as `(referenceSequence, bisulfiteReads)`; a single converted string carries no read multiplicity. See TestSpec §7. |
| 2 | Single-strand conversion (no complement merge) | Assumption | Output reflects supplied strand only | accepted | Frommer protocol is strand-specific [1]. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null/empty sequence | `""` | Validation contract |
| No cytosines | Returned unchanged | Only C reacts [1] |
| Lowercase c | → t (case preserved) | Conversion contract |
| CpG with no covering read | Site excluded | Percentage undefined [3] |
| Read base A/G at CpG C | Call ignored | Not a valid bisulfite call [2] |
| Read past reference end | Extra bases ignored | Out-of-reference [2] |
| Empty site list | Zero profile | Validation contract |

### 6.2 Limitations

CpG-only calling (no CHG/CHH from reads); single-strand conversion; no incomplete-conversion
modeling, sequencing-error filtering, or statistical significance testing. For full-genome
BS-seq, use a dedicated aligner/caller. This is a correctness-focused reference for the core
conversion and calling formulas, not a production read aligner.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Conversion: unmethylated cytosines become T; protected C@1 stays C.
string converted = EpigeneticsAnalyzer.SimulateBisulfiteConversion(
    "ACGTCGAA", new HashSet<int> { 1 });           // "ACGTTGAA"

// Calling: one methylated (C) and one unmethylated (T) read at CpG index 1 -> level 0.5.
var sites = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(
    "ACGTACGT",
    new[] { ("C", 1), ("T", 1) });                 // level 0.5, coverage 2
```

**Numerical walk-through (weighted profile, Schultz 2012 [4]):** two CpG sites with
levels 1.0 (coverage 10) and 0.0 (coverage 30) give weighted CpG methylation
(1.0·10 + 0.0·30)/(10+30) = 10/40 = 0.25 — not the unweighted mean 0.5.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [EpigeneticsAnalyzer_Bisulfite_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/EpigeneticsAnalyzer_Bisulfite_Tests.cs) — covers INV-01..INV-06
- Evidence: [EPIGEN-BISULF-001-Evidence.md](../../../docs/Evidence/EPIGEN-BISULF-001-Evidence.md)
- Related algorithms: [Methylation_Analysis](./Methylation_Analysis.md), [CpG_Site_Detection](./CpG_Site_Detection.md)

## 8. References

1. Frommer M, McDonald LE, Millar DS, Collis CM, Watt F, Grigg GW, Molloy PL, Paul CL. 1992. A genomic sequencing protocol that yields a positive display of 5-methylcytosine residues in individual DNA strands. Proc Natl Acad Sci USA 89(5):1827–1831. https://doi.org/10.1073/pnas.89.5.1827
2. Krueger F, Andrews SR. 2011. Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. Bioinformatics 27(11):1571–1572. https://doi.org/10.1093/bioinformatics/btr167
3. Babraham Bioinformatics. 2016. Bismark Bisulfite Mapper – User Guide v0.15.0. https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide_v0.15.0.pdf
4. Schultz MD, Schmitz RJ, Ecker JR. 2012. 'Leveling' the playing field for analyses of single-base resolution DNA methylomes. Trends Genet. 28(12):583–585. https://doi.org/10.1016/j.tig.2012.10.012
