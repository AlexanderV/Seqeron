# Alpha-Satellite Higher-Order Repeat (HOR) Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Chromosome Analysis |
| Test Unit ID | CHROM-CENT-001 |
| Related Projects | Seqeron.Genomics.Chromosome, Seqeron.Genomics.Alignment |
| Implementation Status | Production |
| Last Reviewed | 2026-06-24 |

## 1. Overview

Human (and primate) centromeric alpha-satellite DNA is organised hierarchically: a block of *N*
diverged ~171 bp monomers forms a **higher-order repeat (HOR) unit**, and that unit is itself tandemly
repeated — hundreds to thousands of times — into near-identical, chromosome-specific arrays.[1][2]
`DetectHigherOrderRepeat` recovers this structure from an alpha-satellite array: it splits the array
into ~171 bp monomers, measures monomer-vs-monomer identity with the library aligner, and reports the
HOR period (monomers per unit), the HOR unit length (bp), the HOR copy number, and the mean inter-HOR
versus intra-HOR identity. It is a deterministic, alignment-based detector with source-derived
identity thresholds and no external trained data. It is an **opt-in, additive** capability that does
not change `AnalyzeCentromere`, `DetectAlphaSatellite`, or the Levan classification.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Alpha satellite is the most abundant centromeric satellite DNA in primates; its fundamental unit is a
~171 bp monomer.[1][2] Two organisational levels matter:[1][2][3]

- **Monomer level (intra-HOR):** the distinct monomers *within one HOR unit* are only **50–70%
  identical** (they "differ in sequence by 10–40%"; "the degree of divergence between any two monomers
  within each HOR copy … rang[es] from approximately 20% to 40%").[1][3]
- **HOR level (inter-HOR):** copies of the *whole HOR unit* are highly homogeneous — "HOR within a
  given array are **97–100% identical**", with "minimal divergence between HOR copies, typically
  **less than 5%**", i.e. "mutual sequence divergence of **<5%**".[1][2][4]

A HOR made of *n* monomers is called an "*n*mer HOR".[2]

### 2.2 Core Model

Index the array's monomers `m[0], m[1], …, m[M-1]` (each ~171 bp). The **HOR period** `k` is the unit
size such that monomers a full period apart are HOR copies of the same position and therefore
near-identical, while monomers within one unit are divergent. The source defines the period
operationally:[1]

> "the HOR unit length is determined by where the next monomer shows nearly total sequence identity to
> the first monomer in the HOR."

Formally, the HOR period is the **smallest** `k ≥ 1` for which the *k*-periodic monomer pairs are
inter-HOR identical consistently across the array:

```text
period = min { k ≥ 1 :  identity(m[i], m[i+k]) ≥ θ_inter  for (≥ ρ of) all valid i }
```

with the inter-HOR identity threshold `θ_inter = 95%` (the <5%-divergence bound)[2][4] and a
consistency fraction `ρ = 0.90` (arrays are highly homogeneous, so the periodicity must hold across
essentially the whole array, not at a single pair).[1] Once `period = k` is found:

- HOR unit length = `k × 171` bp;
- HOR copy number = `⌊M / k⌋`;
- mean inter-HOR identity = mean of `identity(m[i], m[i+k])`;
- mean intra-HOR identity = mean pairwise identity among the `k` distinct monomers of the first unit
  (defined only when `k ≥ 2`).

`period = 1` means adjacent monomers are already near-identical — a **homogeneous single-monomer
(1-mer) array**, which is *not* a multi-monomer HOR.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-HOR-01 | `MonomersPerUnit ≥ 1`; `HasHigherOrderStructure ⇔ MonomersPerUnit ≥ 2`. | The search starts at `k = 1`; a multi-monomer HOR requires a period of at least 2. |
| INV-HOR-02 | `HorUnitLengthBp = MonomersPerUnit × monomerLength`. | Computed directly from the detected period. |
| INV-HOR-03 | `HorCopyNumber = ⌊MonomerCount / MonomersPerUnit⌋`. | Number of complete units tiling the analysed monomers. |
| INV-HOR-04 | For a detected HOR, `MeanInterHorIdentity ≥ 95%` and (when defined) `MeanInterHorIdentity > MeanIntraHorIdentity`. | The period is accepted only at ≥95% inter-HOR identity, while intra-HOR monomers are the divergent 50–70% set — the hallmark of HOR organisation.[1] |
| INV-HOR-05 | Deterministic: identical input (case-insensitive) yields identical output. | Pure function of the sequence; the aligner is deterministic. |

### 2.5 Comparison with Related Methods

| Aspect | `DetectHigherOrderRepeat` | `DetectAlphaSatellite` |
|--------|---------------------------|------------------------|
| Level detected | HOR hierarchy (multi-monomer unit + copy number) | monomer-level signal (171 bp period + AT-richness + CENP-B box) |
| Output | period, unit length, copy number, inter/intra identity | IsAlphaSatellite, periodicity, best period, AT content, CENP-B count |
| Identity measure | gapped global alignment of whole monomers | base-level self-similarity at the monomer period |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | Alpha-satellite array. | `null`/empty or fewer than two full monomers → no-structure result. Case-insensitive (uppercased internally). |
| `monomerLength` | `int` | `171` | Monomer length used to split the array. | Must be ≥ 1 (else `ArgumentOutOfRangeException`). 171 bp = the alphoid monomer.[1] |

### 3.2 Output / Return Value

`HorResult` (read-only record struct):

| Field | Type | Description |
|-------|------|-------------|
| `HasHigherOrderStructure` | `bool` | True iff a HOR period ≥ 2 was detected. |
| `MonomersPerUnit` | `int` | HOR period (monomers per unit); 1 when no multi-monomer HOR. |
| `HorUnitLengthBp` | `int` | `MonomersPerUnit × monomerLength`. |
| `HorCopyNumber` | `int` | `⌊MonomerCount / MonomersPerUnit⌋`. |
| `MonomerCount` | `int` | Number of full monomers the array was split into. |
| `MeanInterHorIdentity` | `double` | Mean percent identity between same-position monomers in different copies; `NaN` when undefined. |
| `MeanIntraHorIdentity` | `double` | Mean percent identity between distinct monomers within one unit; `NaN` when the unit has a single monomer. |

### 3.3 Preconditions and Validation

`monomerLength < 1` throws `ArgumentOutOfRangeException`. `null`/empty sequence returns
`(false, 1, monomerLength, 0, 0, NaN, NaN)`. A sequence with fewer than two full monomers returns a
no-structure result with `MonomerCount` set accordingly. The input is uppercased before alignment; a
trailing partial monomer (sequence length not a multiple of `monomerLength`) is ignored. Coordinates
are monomer indices (0-based); identities are percentages in `[0, 100]`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `monomerLength`; handle `null`/empty and `< 2` monomers up front.
2. Split the uppercased array into `M = ⌊len / monomerLength⌋` consecutive monomers (drop the tail).
3. Define `identity(i, j)` = percent identity of monomers `i, j` via `SequenceAligner.GlobalAlign`
   (Needleman-Wunsch) + `SequenceAligner.CalculateStatistics` (EMBOSS-style identity), memoised.
4. For `k = 1 … ⌊M/2⌋`, accept the first `k` where ≥ 90% of the *k*-periodic monomer pairs are ≥ 95%
   identical. That `k` is the HOR period.
5. If no `k` qualifies, report no structure (period 1, copy number = monomer count).
6. Otherwise compute unit length, copy number, mean inter-HOR identity (i vs i+period), and mean
   intra-HOR identity (distinct monomers of the first unit); set `HasHigherOrderStructure = (period ≥ 2)`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Constant | Value | Origin |
|----------|-------|--------|
| `AlphaSatelliteMonomerLength` | 171 bp | Willard 1985; review PMC6121732[1] |
| `InterHorMinIdentityPercent` | 95.0 | <5% inter-HOR divergence (PMC11050224[2]; Alkan 2007[4]; "97–100% identical" PMC6121732[1]) |
| `HorPeriodConsistencyFraction` | 0.90 | Arrays are "homogenous" / "repeated hundreds to thousands of times"[1] |

Monomer identities are cached in a dictionary keyed by the ordered pair `(min, max)` so each pair is
aligned at most once.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DetectHigherOrderRepeat` | `O(M² · L²)` worst case | `O(M² + L²)` | `M` = monomer count, `L` = monomer length (~171). Each of up to `O(M²)` monomer pairs is a Needleman-Wunsch alignment costing `O(L²)`; identities are memoised, so each pair is aligned once. |

A property/invariant test (INV-HOR-04: inter > intra) and the per-pair memoisation cover the O(n²)
guidance; for typical multi-kb arrays the run completes in milliseconds (see §7.3 tests).

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ChromosomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs)

- `ChromosomeAnalyzer.DetectHigherOrderRepeat(string, int)`: splits the array, builds the memoised
  monomer-identity function, finds the smallest consistent HOR period, and returns a `HorResult`.
- `ChromosomeAnalyzer.HorResult`: the result record struct.

### 5.2 Current Behavior

The detector reuses the repository aligner (`SequenceAligner.GlobalAlign` + `CalculateStatistics`)
rather than a custom Hamming distance, so monomers that differ in length (167–171 bp in real arrays)[4]
are compared by gapped global alignment. The search returns the **smallest** qualifying period, so a
dimeric array is reported as period 2 (not 4 or 6). Period 1 is reported (with
`HasHigherOrderStructure = false`) for homogeneous single-monomer arrays. A trailing partial monomer is
silently dropped. The method does not re-test the alphoid signature; callers gate on
`DetectAlphaSatellite` first if needed.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- HOR period = smallest monomer block size at which the next copy shows near-total identity to the
  first monomer of the unit — the source's operational period definition.[1]
- Inter-HOR acceptance at ≥ 95% identity (< 5% divergence) and intra-HOR monomers in the divergent
  50–70% band; reporting period, unit length, copy number, and both identities.[1][2][4]

**Intentionally simplified:**

- Period consistency is a fixed 0.90 fraction at a fixed 95% identity bar; **consequence:** arrays with
  heavy monomer turnover or cascading/variant HORs may need the thresholds exposed/tuned (both are
  named parameters in code with sourced defaults).

**Not implemented:**

- Suprachromosomal-family / specific α-satellite family (J1/J2/W/…) assignment; **users should rely
  on:** curated reference HOR libraries / dedicated centromere tooling (data-blocked — see §6.2).
- Cascading-HOR decomposition (nested HOR-of-HORs);[2] **users should rely on:** dedicated HOR tools
  (HORmon, alpha-CENTAURI).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Inter-HOR threshold fixed at 95% identity, consistency at 0.90. | Assumption | Borderline-homogeneous arrays could be classified differently under a stricter/looser bound. | accepted | Sourced defaults (<5% divergence); exposed as named constants for future parameterisation. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `null` / empty sequence | `(false, 1, monomerLength, 0, 0, NaN, NaN)` | No monomers to analyse. |
| Fewer than two full monomers | No-structure result with `MonomerCount` set | A single monomer cannot show periodicity. |
| Homogeneous single-monomer array | `MonomersPerUnit = 1`, `HasHigherOrderStructure = false`, `MeanInterHorIdentity = 100` | 1-mer array is not a multi-monomer HOR. |
| Monomeric, mutually divergent array (no pair ≥95%) | `MonomersPerUnit = 1`, `HasHigherOrderStructure = false`, `MeanInterHorIdentity = NaN` | No period clears the inter-HOR bar. |
| Trailing partial monomer | Ignored; period/copy number from full monomers only | The tail is not a complete monomer. |
| `monomerLength < 1` | `ArgumentOutOfRangeException` | Invalid split width. |
| Lowercase input | Same result as uppercase | Sequence is uppercased before alignment. |

### 6.2 Limitations

The detector recovers HOR *structure* (period, copy number, inter-/intra-HOR identity) but does **not**
assign a **suprachromosomal family or specific α-satellite family label** (J1/J2/W/…): that requires
curated chromosome-specific reference HOR libraries, which are external trained/curated data not
embedded in the library (data-blocked). It also does not decompose cascading / nested HORs.[2] Inputs
should already be alpha-satellite arrays; the method does not itself verify the alphoid signature.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// array = a 3-monomer HOR unit (A,B,C ~58% mutually identical) repeated 5× with exact copies.
var hor = ChromosomeAnalyzer.DetectHigherOrderRepeat(array);
// hor.HasHigherOrderStructure == true
// hor.MonomersPerUnit       == 3      (3-monomer HOR unit)
// hor.HorCopyNumber         == 5      (15 monomers / 3)
// hor.HorUnitLengthBp       == 513    (3 × 171)
// hor.MeanInterHorIdentity  == 100.0  (exact copies)
// hor.MeanIntraHorIdentity  ~= 57.9   (distinct monomers, 50–70% band)  < inter-HOR
```

### 7.2 Applications and Use Cases

- **Centromere array characterisation:** distinguishing a true multi-monomer HOR array from a
  homogeneous monomeric repeat, and estimating the HOR copy number of an assembled centromere.[1]

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ChromosomeAnalyzer_HigherOrderRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_HigherOrderRepeat_Tests.cs) — covers `INV-HOR-01`..`INV-HOR-05`
- Evidence: [CHROM-CENT-001-Evidence.md](../../../docs/Evidence/CHROM-CENT-001-Evidence.md)
- Related algorithms: [Centromere_Analysis.md](Centromere_Analysis.md)

## 8. References

1. McNulty SM, Sullivan BA. 2018. Alpha satellite DNA biology: finding function in the recesses of the genome. Chromosome Research 26:115–138. https://pmc.ncbi.nlm.nih.gov/articles/PMC6121732/
2. Rosandić M, Paar V, et al. 2024. Novel Concept of Alpha Satellite Cascading Higher-Order Repeats (HORs) … in the T2T-CHM13 Assembly of Human Chromosome 15. https://pmc.ncbi.nlm.nih.gov/articles/PMC11050224/
3. Warburton PE, Willard HF. 1990. Genomic analysis of sequence variation in tandemly repeated DNA. Evidence for localized homogeneous sequence domains within arrays of alpha-satellite DNA. Journal of Molecular Biology 216(1):3–16. https://doi.org/10.1016/S0022-2836(05)80056-7
4. Paar V, Pavin N, Rosandić M, et al. (ColorHOR); Alkan C, et al. 2007. Genome-wide characterization of centromeric satellites. Bioinformatics 21(7):846–852. https://academic.oup.com/bioinformatics/article/21/7/846/268781
