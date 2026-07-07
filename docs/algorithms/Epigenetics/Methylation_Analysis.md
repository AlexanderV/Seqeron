# Methylation Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Epigenetics |
| Test Unit ID | EPIGEN-METHYL-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Methylation analysis classifies each cytosine in a DNA sequence into its methylation sequence context — CpG, CHG, or CHH — and aggregates per-cytosine methylation measurements into a per-context profile [2][3]. Context classification is an exact, deterministic lookup over the cytosine and its one or two downstream bases [1][2]. The profile reports weighted methylation levels per context, defined as the ratio of methylated reads to total reads (Schultz et al. 2012) [4]. The implementation exposes a single-site context classifier, a whole-sequence site enumerator, and a profile aggregator in one analyzer class.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In mammals DNA methylation occurs almost exclusively at CpG dinucleotides; in plant and stem-cell genomes methylation also occurs at non-CpG cytosines in the CHG and CHH contexts [3]. The three contexts are distinguished by the bases immediately 3' of the cytosine [2]. CpG and CHG are *symmetric* (the context recurs on the complementary strand), whereas CHH is *asymmetric* [2].

### 2.2 Core Model

For a cytosine at 0-based position *i* in a 5'→3' sequence, with downstream bases *b₁ = s[i+1]* and *b₂ = s[i+2]*:

- **CpG** ⟺ *b₁* = G [2][3].
- **CHG** ⟺ *b₁* ∈ H and *b₂* = G [1][2].
- **CHH** ⟺ *b₁* ∈ H and *b₂* ∈ H [2].

where **H = {A, C, T}** is the IUPAC ambiguity code for "not G" (Cornish-Bowden 1985) [1].

The per-context **weighted methylation level** (Schultz et al. 2012) is [4]:

$$
\mathrm{WML}_{\text{context}} = \frac{\sum_{c}\ \text{methylated reads}(c)}{\sum_{c}\ \text{total reads}(c)}
= \frac{\sum_c \mathrm{level}(c)\cdot \mathrm{coverage}(c)}{\sum_c \mathrm{coverage}(c)}
$$

summed over all cytosines *c* of that context. The per-cytosine methylation level is `methylated reads / total reads` ∈ [0, 1] [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A cytosine is CpG ⟺ the next base is G | CpG definition [2][3] |
| INV-02 | CHG/CHH require the H position(s) ∈ {A,C,T}; a G is never H | IUPAC H excludes G [1] |
| INV-03 | Reported site positions are 0-based and index a cytosine | enumeration contract [2] |
| INV-04 | Each per-site methylation level ∈ [0,1] | definition of a read fraction [4] |
| INV-05 | Per-context profile level = Σ(level·coverage)/Σ(coverage); equals the mean of fractions only under equal coverage | weighted methylation level [4] |
| INV-06 | Classification is case-insensitive and deterministic | bases upper-cased before comparison |

### 2.5 Comparison with Related Methods

| Aspect | Weighted methylation level (this profile) | Unweighted mean of fractions |
|--------|-------------------------------------------|------------------------------|
| Coverage weighting | Yes (Σmeth/Σtotal) [4] | No |
| Effect of low-coverage outliers | Down-weighted | Equal weight, can dominate |
| Equal-coverage case | Identical to the mean | Identical to weighted |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | DNA sequence | case-insensitive; bases compared as {A,C,G,T}; 0-based |
| index | int | required | Position to classify (`GetMethylationContext`) | any int; out-of-range → null |
| sites | IEnumerable\<MethylationSite\> | required | Measured sites (`GenerateMethylationProfile`) | each carries Type, Level∈[0,1], Coverage≥0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| GetMethylationContext | MethylationType? | CpG/CHG/CHH, or null when not a classifiable cytosine |
| FindMethylationSites | IEnumerable\<MethylationSite\> | one site per classifiable C: Position (0-based), Type, Context string, Level=0, Coverage=0 |
| MethylationProfile | record struct | Global/CpG/CHG/CHH weighted levels, TotalCpGSites, MethylatedCpGSites, MethylationByPosition |

### 3.3 Preconditions and Validation

Null/empty sequence → empty enumeration (sites) or null (context). Index < 0 or ≥ length → null. A base that is not a cytosine, an incomplete downstream window, or a context base outside {A,C,G,T} → null/unclassified (the cytosine is skipped). Case is normalized to uppercase. Coordinates are 0-based. Sequence-only sites carry methylation level 0 and coverage 0 because no bisulfite read evidence is present; an empty site list yields an all-zero profile.

## 4. Algorithm

### 4.1 High-Level Steps

1. For a target cytosine, read the next base; if G → CpG.
2. Otherwise require the next base to be H (A/C/T); if not, unclassified.
3. Read the third base; if G → CHG, else if H → CHH, else unclassified.
4. `FindMethylationSites` applies steps 1–3 to every position, yielding one site per classifiable C.
5. `GenerateMethylationProfile` partitions sites by context and computes the weighted methylation level Σ(level·coverage)/Σ(coverage) per context, plus CpG counts.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

IUPAC ambiguity code **H = A, C, or T** (Cornish-Bowden 1985) [1]; the classifier rejects G in any H position. The "methylated CpG" count uses a descriptive fractional cutoff of 0.5 (assumption; does not affect continuous levels). When total coverage of a context is 0 (sequence-only sites), the profile falls back to the unweighted mean of per-site levels so the sites are not dropped.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| GetMethylationContext | O(1) | O(1) | constant-base lookup |
| FindMethylationSites | O(n) | O(1) per yielded site | single linear scan |
| GenerateMethylationProfile | O(m) | O(m) | m = number of sites; partition + sums |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [EpigeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs)

- `EpigeneticsAnalyzer.GetMethylationContext(string, int)`: classifies one cytosine into CpG/CHG/CHH or null.
- `EpigeneticsAnalyzer.FindMethylationSites(string)`: enumerates all classifiable cytosine sites with context.
- `EpigeneticsAnalyzer.GenerateMethylationProfile(IEnumerable<MethylationSite>)`: per-context weighted methylation profile (Registry name: `CalculateMethylationProfile`).

### 5.2 Current Behavior

`FindMethylationSites` delegates classification to `GetMethylationContext`, guaranteeing the two entry points agree. Site `Context` strings are the up-to-3-base window upper-cased; terminal CpG yields a 2-base context. Profile levels are weighted by coverage (Schultz 2012); the zero-coverage fallback path supports sequence-only sites. This unit does not involve substring search/occurrence enumeration, so the repository suffix tree is **not used** — classification is a constant-window lookup per position, for which a single O(n) scan is optimal.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- CpG/CHG/CHH classification with H = A/C/T (G excluded) [1][2][3].
- Weighted methylation level Σ(meth)/Σ(total) per context [4].
- 0-based cytosine-anchored site positions [2].

**Intentionally simplified:**

- MethylatedCpGSites count: uses a 0.5 fractional cutoff instead of the per-site binomial test [4]; **consequence:** the count field is descriptive, not a statistically tested call; continuous methylation levels are unaffected.
- Sequence-only sites: methylation level/coverage default to 0; **consequence:** profiles built solely from `FindMethylationSites` (no reads) report level 0 until bisulfite data is supplied.

**Not implemented:**

- Strand-aware (top/bottom) symmetric-context pairing and read alignment; **users should rely on:** `CalculateMethylationFromBisulfite` for read-derived levels, or external tools (Bismark) for full WGBS calling.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Profile method named `GenerateMethylationProfile` (Registry: `CalculateMethylationProfile`) | Deviation | naming only; same behavior | accepted | kept for API stability (used by MCP tools/prior tests) |
| 2 | 0.5 cutoff for MethylatedCpGSites count | Assumption | count field only | accepted | Schultz (2012) recommends binomial test |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Non-cytosine index | null | only cytosines have methylation context [2] |
| Terminal `CG` | CpG | CpG needs only two bases [2] |
| Terminal `C` + one H | null (CHG/CHH need a third base) | window incomplete [2] |
| G in H position (`CGG`) | CpG (next is G) | H excludes G [1] |
| Non-ACGT in context (`CNG`) | null | undetermined context [1] |
| Empty/null sequence | empty / null | validation contract |
| Empty site list | all-zero profile | validation contract |

### 6.2 Limitations

Single-strand classification only: it does not pair symmetric CpG/CHG contexts across strands, does not perform read alignment, and treats degenerate/ambiguous (non-ACGT) bases as unclassifiable rather than expanding them. Methylation levels are only as meaningful as the coverage supplied; sequence-only input yields zero-level placeholder sites.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Classify a single cytosine and enumerate all sites.
var ctx = EpigeneticsAnalyzer.GetMethylationContext("CGACAGCAA", 3); // CHG (C,A,G)
var sites = EpigeneticsAnalyzer.FindMethylationSites("CGACAGCAA").ToList();
// CpG@0, CHG@3, CHH@6

var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(new[]
{
    new EpigeneticsAnalyzer.MethylationSite(0, EpigeneticsAnalyzer.MethylationType.CpG, "CG", 0.8, 10),
    new EpigeneticsAnalyzer.MethylationSite(5, EpigeneticsAnalyzer.MethylationType.CpG, "CG", 0.2, 10),
});
// profile.CpGMethylation == (8 + 2) / (10 + 10) == 0.5  (weighted, Schultz 2012)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [EpigeneticsAnalyzer_Methylation_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/EpigeneticsAnalyzer_Methylation_Tests.cs) — covers `INV-01`–`INV-06`
- Evidence: [EPIGEN-METHYL-001-Evidence.md](../../../docs/Evidence/EPIGEN-METHYL-001-Evidence.md)
- Related algorithms: [CpG_Site_Detection](./CpG_Site_Detection.md)

## 8. References

1. Cornish-Bowden A. 1985. Nomenclature for incompletely specified bases in nucleic acid sequences (IUPAC-IUB). Nucleic Acids Research 13(9):3021–3030. https://doi.org/10.1093/nar/13.9.3021
2. Krueger F, Andrews SR. 2011. Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. Bioinformatics 27(11):1571–1572. https://doi.org/10.1093/bioinformatics/btr167
3. Lister R, Pelizzola M, Dowen RH, et al. 2009. Human DNA methylomes at base resolution show widespread epigenomic differences. Nature 462(7271):315–322. https://doi.org/10.1038/nature08514
4. Schultz MD, Schmitz RJ, Ecker JR. 2012. 'Leveling' the playing field for analyses of single-base resolution DNA methylomes. Trends in Genetics 28(12):583–585. https://doi.org/10.1016/j.tig.2012.10.012
