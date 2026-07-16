---
type: concept
title: "Selection-signature detection (iHS / EHH haplotype scan)"
tags: [population-genetics, algorithm]
mcp_tools:
  - integrated_haplotype_score
sources:
  - docs/Evidence/POP-SELECT-001-Evidence.md
  - docs/algorithms/Population_Genetics/Integrated_Haplotype_Score.md
source_commit: 954cc1e1b3751da20446ced26b4265a0c543dd91
created: 2026-07-10
updated: 2026-07-16
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pop-select-001-evidence
      evidence: "Test Unit ID: POP-SELECT-001 ... Methods: CalculateEhh, CalculateIHS, StandardizeIHS, ScanForSelection (population-genetics unit)."
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:allele-genotype-frequencies
      source: pop-select-001-evidence
      evidence: "iHS is computed only for SNPs with ancestral state and MAF > 5%, and StandardizeIHS normalizes within derived-allele-frequency bins — it consumes the per-locus allele/derived frequencies."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:linkage-disequilibrium
      source: pop-select-001-evidence
      evidence: "EHH is the probability that two core-allele chromosomes are homozygous over a surrounding region — extended haplotype homozygosity, the same multi-locus haplotype-association structure LD quantifies."
      confidence: medium
      status: current
---

# Selection-signature detection (iHS / EHH haplotype scan)

Detect the genomic footprint of **recent positive natural selection** from **haplotype structure**:
an allele driven up in frequency by selection carries an unusually **long, low-diversity haplotype**
around it (a "selective sweep"), because there has not yet been time for recombination to break it
down. This is a population-genetics `POP-*` unit (**POP-SELECT-001**) in the family anchored by
[[ancestry-estimation-admixture]]. It is genuinely distinct from its POP siblings: it is a
**haplotype-length / decay** statistic, not a per-locus count ([[allele-genotype-frequencies]]), a
within-sample diversity summary ([[genetic-diversity-statistics]]), a between-population
differentiation scalar ([[population-differentiation-fst]]), a single-locus goodness-of-fit test
([[hardy-weinberg-equilibrium-test]]), a pairwise inter-locus association ([[linkage-disequilibrium]]),
or a within-genome autozygosity scan ([[runs-of-homozygosity-inbreeding]]). It **consumes** the
derived/ancestral allele frequencies of [[allele-genotype-frequencies]] (MAF filter + frequency-bin
standardization). Validated under **POP-SELECT-001**; the literature-traced record is
[[pop-select-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## The pipeline: EHH → iHH → iHS → scan

**1. EHH (Extended Haplotype Homozygosity)** — Sabeti et al. (2002). At a genomic offset from a
**core SNP**, EHH is the probability that two randomly chosen chromosomes **carrying the same core
allele** are identical (homozygous) over the whole stretch between the core and that offset. In the
combinatorial form (selscan Eq. 3, algebraically identical in rehh):

```
EHH_c(x) = Σ_{h ∈ H_c(x)} C(n_h, 2) / C(n_c, 2)
```

where `n_c` is the number of chromosomes carrying core allele `c` and `n_h` is the count of each
distinct extended haplotype among them. EHH starts at 1.0 at the core and **decays toward 0** with
distance as recombination and mutation diversify the flanks. A sweep keeps EHH high for far longer on
the selected allele.

**2. iHH (integrated EHH)** — the **area under the EHH-vs-position curve**, computed separately for
the ancestral (`iHH_A`) and derived (`iHH_D`) core allele and **summed in both directions** from the
core. Integration uses the **trapezoidal rule** (EHH values at successive markers joined by straight
lines) and is **truncated at the first marker where EHH drops below 0.05** (the `limehh` cutoff).

**3. iHS (integrated Haplotype Score)** — the log-ratio of the two integrated curves:

```
unstandardized iHS = ln(iHH_A / iHH_D)      # Voight et al. 2006: ancestral numerator
```

Then **standardized** to an approximately standard-normal score within bins of matching
**derived-allele frequency**, so scores at different frequencies are directly comparable:

```
iHS = ( ln(iHH_A/iHH_D) − E_p[ln(iHH_A/iHH_D)] ) / SD_p[ln(iHH_A/iHH_D)]
```

**4. Genome scan** (`ScanForSelection`) — slide a window (standard **50 SNPs**) and report the
**proportion of SNPs with |iHS| > 2**; high-proportion windows are candidate selection regions.

## Sign convention — declare it or get it backwards

The iHS sign is **not universal**. This unit follows **Voight (2006): `ln(iHH_A / iHH_D)`** (ancestral
over derived), so:

- **Large negative iHS ⇒ long haplotype on the *derived* allele** (the common "derived allele under
  selection" signal).
- **Large positive iHS ⇒ long haplotype on the *ancestral* allele.**
- **Balanced EHH decay ⇒ `iHH_A/iHH_D ≈ 1` ⇒ iHS ≈ 0** (no signal).

The *selscan* reference implementation uses the **opposite** numerator (`ln(iHH_1/iHH_0)`,
derived/ancestral) — same magnitude, flipped sign. Because `ln(a/b) = −ln(b/a)`, an implementation
**must state which convention it uses**; the Evidence file records this as a documented convention
difference, not a contradiction.

## Worked oracles

- **rehh SNP F1205400:** `iHH_A = 284429.9`, `iHH_D = 2057107.4` ⇒ ln(A/D) = **−1.978569274**.
- **Constructed panel** (3 identical derived haplotypes vs 3 all-distinct ancestral): `EHH_D = 1.0`
  at each flank, `EHH_A = 0.0` (truncated immediately), `iHH_D = 40.0`, `iHH_A = 10.0` ⇒
  iHS = ln(0.25) = **−1.386294361**, derived freq 0.5 — a long derived haplotype giving the expected
  negative iHS.
- **EHH unit values:** `11,11,11,10` (n_c=4) ⇒ 3/6 = **0.5**; `00,00,01,01` ⇒ 2/6 = **0.3333**;
  single haplotype ⇒ **1.0**; three all-distinct ⇒ **0.0**.

## Invariants and edge cases

- **EHH ∈ [0, 1]**: 1.0 for a single chromosome (trivially homozygous), 0.0 when all haplotypes are
  distinct; empty core-allele sample ⇒ 0.
- **Integration** always truncates at EHH < 0.05 and is symmetric (both directions from the core).
- **Standardization:** single-element frequency bin ⇒ standardized iHS = 0 (no spread); the SD uses
  the sample (N−1) estimator; frequency binning defaults to 20 bins of width 0.05 (`freqbin`),
  overridable.
- **`ScanForSelection` window value = (count of |iHS| > 2) / window size.**
- **Throws / no value:** monomorphic core (throws), null inputs, empty or inconsistent-length
  haplotypes, invalid core allele, out-of-range `coreIndex`; SNPs lacking ancestral state or with
  MAF ≤ 5% get no iHS.
- **Property:** `ln(iHH_A/iHH_D) = −ln(iHH_D/iHH_A)` (Voight ↔ selscan sign symmetry).

## Implementation notes

Implemented in `PopulationGeneticsAnalyzer` (`Seqeron.Genomics.Population`) as four entry points:
`CalculateEhh(IReadOnlyList<string>)` (EHH of one core-allele subset), `CalculateIHS(haplotypes,
positions, coreIndex)` → `IhsResult` (unstandardized Voight iHS plus `IhhAncestral`/`IhhDerived`/
`DerivedAlleleFrequency`), `StandardizeIHS(scores, binCount=20)`, and `ScanForSelection(
standardizedScores, windowSize=50)` → `SelectionScanWindow` records
(`WindowIndex, SnpCount, ExtremeCount, ProportionExtreme`).

- **EHH counts literal whole-window haplotypes by hashing substrings**, so any alphabet is accepted
  for the flanking markers — only the *core* character is constrained to `'0'`/`'1'`. This is **not**
  a substring-search problem (it counts exact whole-window multiplicities, not occurrences of a query
  pattern in a text), so the repository suffix tree deliberately does **not** apply.
- **Cost:** `CalculateEhh` is O(n·L) (n chromosomes, L window length); `CalculateIHS` is O(n·h)
  (h markers) because EHH is **recomputed per outward marker** rather than incrementally — exact but
  not tuned for genome-scale panels. `StandardizeIHS` and `ScanForSelection` are single passes, O(N).
- A separate pre-existing overload `CalculateIHS(ehh0, ehh1, positions)` (EHH curves supplied
  directly) and a region-threshold `ScanForSelection` variant back the MCP layer and are out of scope
  for this unit.
- **Versus site-frequency-spectrum scans:** iHS keys on haplotype-length *asymmetry* around one
  allele within a single population (best for incomplete/ongoing sweeps at intermediate frequency),
  whereas Tajima's D / Fay & Wu's H key on SFS distortion (completed sweeps / frequency skew) —
  see [[genetic-diversity-statistics]].

## Scope

Faithful implementation of the EHH/iHH/iHS haplotype-scan family (exact match to Voight et al. 2006,
Sabeti et al. 2002, selscan Szpiech & Hernandez 2014, and the rehh package). It computes single-core
EHH, the trapezoidal integrated iHH with the 0.05 cutoff, the unstandardized Voight iHS, its
frequency-binned standardization, and a windowed |iHS| > 2 scan. It does **not** implement the
cross-population variants (XP-EHH), Fst-outlier or Tajima's-D-based scans, genetic-map interpolation
of marker spacing, or the >200 kb gap / chromosome-end data-curation masking. The two documented
assumptions (N−1 SD estimator; default bin width) affect only standardized magnitude, never the sign,
ordering, or the canonical unstandardized iHS. No source contradictions; Open Questions: none.
