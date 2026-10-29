# Off-Target Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | CRISPR-OFF-001 |
| Related Projects | N/A |
| Implementation Status | Implemented (published substitution models, bulge enumeration) + one declared heuristic |
| Last Reviewed | 2026-10-25 (D19 bulge enumeration) |

## 1. Overview

Off-target analysis estimates potential unintended cleavage sites for CRISPR guide RNAs by scanning a genome or reference sequence for near matches adjacent to valid PAMs. The established `FindOffTargets` path uses guide-length targets with substitutions and feeds its hits to the existing specificity scorers. `FindOffTargetsWithBulges` separately enumerates PAM-supported DNA/RNA bulge alignments with variable-length protospacers. It reports alignments, not activity scores.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

CRISPR off-target activity arises when a guide RNA binds and cleaves genomic loci that are similar, but not identical, to the intended target. The original document highlights the importance of the PAM-proximal seed region, cites observed off-target activity with 3-5 mismatches, and notes that PAM-proximal mismatches are especially important for specificity. Sources: Hsu et al. (2013), Fu et al. (2013), Wikipedia (Off-target genome editing).

### 2.2 Core Model

The repository's off-target search applies three criteria:

1. The candidate site must have a valid PAM for the selected CRISPR system.
2. The PAM-adjacent target must have the correct guide length.
3. The guide-to-target mismatch count must satisfy `0 < mismatches <= maxMismatches`.

The sibling `FindOffTargetsWithBulges` follows [Cas-OFFinder 3.0](https://github.com/snugel/cas-offinder/blob/develop/README.md): it scans PAM-adjacent windows from `guideLength + maxDnaBulge` through `guideLength - maxRnaBulge`. A DNA bulge is extra genomic bases and gaps the aligned guide; an RNA bulge is extra guide bases and gaps the aligned target. The gap is contiguous and has size 1 or 2. The admissible gap offsets are those of Cas-OFFinder 3.0.0b3 (`cas-offinder.cpp` `parseInput`), in guide coordinates for a guide of length `L`: a DNA bulge may sit at `1..L` when the PAM follows the target (SpCas9, SaCas9) and at `0..L-1` when it precedes it (Cas12a, CasX) — the remaining offset, on the side of the tool's `N` padding, rebuilds the ungapped query, which the tool reports as type `X` and which here is the ungapped alignment; an RNA bulge removes a run starting at `0..L-size-1`, so the guide base at the end away from the padding is never removable. The only deliberate difference is at the sequence ends: Cas-OFFinder pads its pattern with `maxDnaBulge` wildcard columns and therefore cannot match a PAM whose padded window crosses the end of a sequence; every alignment inside the supplied sequence is reported here. Mismatches are counted only where both aligned strings have bases. Exact ungapped matches are excluded by the existing no-on-target-locus convention — or, when `onTargetStart`/`onTargetEnd` are given, only when their protospacer+PAM footprint lies inside that interval (CRISPOR `annotateOfftargets`, as in the specificity aggregates); zero-mismatch bulged alignments are reported, so a perfect site also yields the 0-mismatch RNA-bulge rows that remove a leading guide base, as the reference does. Each hit carries its PAM (`BulgedOffTargetSite.Pam`), which the reference prints as the last columns of its aligned strings. The method retains the same PAM-orientation support as `FindPamSites`, returns the PAM coordinate plus the forward-strand protospacer start, and scans a supplied sequence directly rather than a genome-scale index. Cas-OFFinder's original [paper](https://doi.org/10.1093/bioinformatics/btu048) describes the mismatch-search tool; the 3.0 implementation/README supplies the bulge behavior.

All three specificity scorers (`CalculateHeuristicSpecificityScore`, `CalculateMitSpecificityScore`, `CalculateCfdSpecificityScore`) remain substitution-only. Neither the MIT nor CFD model defines a weight for a bulged position, and the heuristic has no sourced bulge weight. Bulged hits are never sent to them. The exact Cas-OFFinder oracle rows and reproducible input are in [CRISPR-OFF-001-Evidence.md](../../Evidence/CRISPR-OFF-001-Evidence.md).

The off-target score is computed by summing mismatch penalties by position:

| Region | Weight |
|--------|--------|
| Seed region mismatch | `5` |
| Non-seed mismatch | `2` |

The seed region is the last 12 bp for PAM-after-target systems such as Cas9 and the first 12 bp for PAM-before-target systems such as Cas12a. `CalculateHeuristicSpecificityScore(...)` finds off-targets with up to 4 mismatches, sums their penalties, and returns `max(0, 100 - totalPenalty)`. The word "heuristic" is in the method name on purpose: the weights `5`/`2`, the 12-nt seed width, the linear sum and the subtraction from 100 are unsourced (Hsu 2013 motivates only the direction of the seed weighting). The historical shorter name `CalculateSpecificityScore` is kept as a backward-compatible alias that forwards to it, and the MCP tool name `crispr_specificity_score` is unchanged.

> **`CalculateSpecificityScore` is this library's own heuristic, not a published model.** Hsu et al. (2013) motivates only the *direction* of the seed weighting; the weights `5`/`2`, the 12-nt seed width, the linear sum and the subtraction from 100 are unsourced. It is kept as a coarse triage number and declared as such in its XML doc and in the `crispr_specificity_score` MCP tool description.

### 2.3 Published models (reference-verified)

Three published off-target models live in the same class and are the ones to use for any
publication-grade number. All are verified numerically against CRISPOR
(`raw.githubusercontent.com/maximilianh/crisporWebsite/master/crispor.py`); see
[CRISPR-OFF-001-Evidence.md](../../Evidence/CRISPR-OFF-001-Evidence.md) for the opened line numbers
and every reference value.

| Model | Method | Form |
|---|---|---|
| MIT / Hsu 2013 single hit | `CalculateMitHitScore(guide20, offTarget20)` | the published position-weight product, ×100 |
| SaCas9, Tycko et al. 2018 | `CalculateSaCas9HitScore(guide21, offTarget21)` | product of trained per-position mismatch coefficients, ×100 |
| Doench 2016 CFD single hit | `CalculateCfdScore(guide20, offTarget20, offTargetPam)` | mismatch-matrix product × PAM-activity score |

Both guide-level aggregates follow `annotateOfftargets`:

| Aggregate | Method | Form |
|---|---|---|
| MIT / Hsu guide score | `CalculateMitSpecificityScore(...)` | `100 / (100 + Σ hits) × 100`, unrounded |
| CFD guide score | `CalculateCfdSpecificityScore(..., CfdAggregation.MitStyle)` | `100 / (100 + Σ raw CFD) × 100` (the reference default) |
| (all four genome aggregates) | optional `onTargetStart` / `onTargetEnd` | CRISPOR's on-target interval skip; see *Intentionally simplified* |
| CFD guide score over the alternative PAMs | `CalculateCfdSpecificityScoreWithAlternativePams(..., minAlternativePamHitScore = 1.0)` | CRISPOR's wider NAG/NGA search with `calcCfdScore` on every hit; same hit set and MIT-score cutoff as the MIT method, no 0.2 factor (B08 round 3 G1) |
| CFD guide score, `doCfdFix` | `CalculateCfdSpecificityScore(..., CfdAggregation.ParkinsonNormalized)` | `100 / Σ raw CFD`; undefined at Σ = 0 |

The genome-scanning MIT aggregate picks the per-hit model by system exactly as the reference does: a
SaCas9 PAM with a 21-nt guide uses the Tycko 2018 model, any other PAM-after-target system uses the
Hsu model. PAM-before-target systems (Cas12a-class, CasX) are **refused** by both aggregates, because
the reference reports no MIT or CFD off-target score for them (`mitScore = cfdScore = 0.0`,
`guideScore = guideCfdScore = -1`). CFD additionally requires a 20-nt protospacer, so SaCas9 is
refused there too — the reference emits its CFD sentinel (`cfdScore = -1`) for SaCas9 hits.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Exact matches are not returned as off-targets | `FindOffTargets(...)` yields only when `mismatches > 0` |
| INV-02 | Returned mismatch counts are bounded by the requested `maxMismatches` | The source filters on `mismatches <= maxMismatches` |
| INV-03 | Specificity scores are clamped to the range `[0, 100]` | `CalculateSpecificityScore(...)` applies `Math.Max(0, 100 - totalPenalty)`; the MIT and `MitStyle` CFD aggregates are `100/(100+Σ)×100` with `Σ ≥ 0`. `CfdAggregation.ParkinsonNormalized` is deliberately **not** bounded by 100 (the reference's `calcCfdGuideScore` is not either) |
| INV-04 | `MismatchPositions.Count == Mismatches` for every returned site | both are derived from one case-folding comparison (`GetMismatchPositions`), which also drives the off-target penalty, so the three cannot disagree |
| INV-05 | Enumeration order is deterministic | forward-strand PAM sites in ascending position, then reverse-strand sites |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `guideSequence` | `string` | required | Guide RNA sequence to analyze | Must match the guide length of the selected CRISPR system |
| `genome` | `DnaSequence` | required | Genome or reference sequence to search | Null input throws `ArgumentNullException` |
| `maxMismatches` | `int` | `3` | Maximum number of mismatches allowed in `FindOffTargets(...)` | Must be between `0` and `5` inclusive |
| `maxDnaBulge`, `maxRnaBulge` | `int` | `2` each | Largest single contiguous DNA/RNA gap in `FindOffTargetsWithBulges(...)` | Each must be between `0` and `2` inclusive |
| `systemType` | `CrisprSystemType` | `SpCas9` | CRISPR system definition used for PAM and guide length | Must resolve through `GetSystem(...)` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | Position of the off-target PAM site |
| `Sequence` | `string` | Target sequence adjacent to the PAM |
| `Mismatches` | `int` | Guide-to-target mismatch count |
| `MismatchPositions` | `IReadOnlyList<int>` | Zero-based mismatch positions |
| `IsForwardStrand` | `bool` | Strand orientation of the off-target site |
| `OffTargetScore` | `double` | Position-weighted mismatch penalty |
| `specificity` | `double` | Overall specificity score in the range `0-100` |

### 3.3 Preconditions and Validation

**Guide alphabet (B08 round 3 G2).** Every raw-string guide entry point (`FindOffTargets`, `FindOffTargetsWithPam`, the heuristic / MIT / CFD genome aggregates, the MIT aggregate over hand-supplied hits, `CalculateMitSpecificityScoreWithAlternativePams`, `CalculateMitHitScore`, `GetGrafMotif`) applies the rule of `EvaluateGuideRna` (F70): ASCII-only upper-casing, `U` read as `T` (an RNA-spelled guide gives the same result as its DNA spelling), any other character (`N`, IUPAC codes, `-`, space, `ſ`) is an `ArgumentException` naming the parameter, raised after the null/empty, genome, `maxMismatches` and length guards. Before, an `N`/junk guide mismatched every site and the aggregates returned a false perfect 100. `CalculateCfdScore` and `CalculateSaCas9HitScore` keep their own stricter A/C/G/T check (no `U` reading; the CFD reference screens `[^ATCG]`).

`FindOffTargets(...)` throws `ArgumentNullException` for null guide or genome input, `ArgumentException` for an empty guide (the null/empty split of `ArgumentException.ThrowIfNullOrEmpty`, audit round 4 H19 / round 5 A5-1) and `ArgumentOutOfRangeException` when `maxMismatches` is outside `0-5`. It also throws `ArgumentException` when the guide length does not match the expected guide length for the selected CRISPR system. `CalculateSpecificityScore(...)` always evaluates off-targets using a fixed mismatch cap of 4. Both genome-scanning published aggregates enforce the same `0-5` bound on their own `maxMismatches` parameter, and both throw `ArgumentException` for a PAM-before-target system; `CalculateCfdSpecificityScore(...)` additionally throws `ArgumentException` for any system whose guide length is not 20, and `ArgumentOutOfRangeException` for `CfdAggregation.ParkinsonNormalized` when the CFD sum is 0.

## 4. Algorithm

### 4.1 High-Level Steps

1. Resolve the CRISPR system into its PAM pattern, guide length, and orientation.
2. Find all PAM sites in the genome on both strands.
3. Extract the guide-length target sequence for each PAM site.
4. Count mismatches between the guide and each target.
5. Yield sites whose mismatch count is positive and within the configured threshold.
6. Score each site by assigning a higher penalty to mismatches inside the 12-base seed region.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

CRISPR system variations documented for off-target analysis:

| System | PAM | Guide Length | PAM Position | Seed Region |
|--------|-----|--------------|--------------|-------------|
| SpCas9 | NGG | 20 bp | After target | Last 12 bp |
| SaCas9 | NNGRRT | 21 bp | After target | Last 12 bp |
| Cas12a | TTTV | 23 bp | Before target | First 12 bp |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindOffTargets` | `O(n × m)` | `O(k)` | `n` is genome length, `m` is guide length, `k` is number of off-targets |
| `FindOffTargetsWithBulges` | `O(n × m² × (d + r + 1))` | `O(n + k)` | `d`/`r` are the bounded DNA/RNA gap ceilings; every gap offset is compared across a guide-length window, and the reverse complement is materialized |
| `CalculateSpecificityScore` | `O(n × m)` | `O(k)` | Reuses `FindOffTargets(...)` with a mismatch cap of 4 |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CrisprDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs)

- `CrisprDesigner.FindOffTargets(string, DnaSequence, int, CrisprSystemType)`: Finds PAM-constrained near matches for a guide sequence.
- `CrisprDesigner.FindOffTargetsWithBulges(string, DnaSequence, int, int, int, CrisprSystemType, int?, int?)`: Finds PAM-constrained substitution and single-bulge alignments (optional on-target interval as above); returns `BulgedOffTargetSite` with the aligned strings, gap type/offset, PAM coordinate and sequence, and protospacer start. Throws `ArgumentNullException` (null guide, null genome), `ArgumentOutOfRangeException` (`maxMismatches` outside 0-5, a bulge bound outside 0-2, a negative or empty interval) and `ArgumentException` (empty guide, guide length, guide alphabet, one interval bound only). No activity score is implied.
- `CrisprDesigner.CalculateSpecificityScore(string, DnaSequence, CrisprSystemType)`: Converts the heuristic off-target penalties into a `0-100` score (declared heuristic — see §2.2).
- `CrisprDesigner.CalculateMitHitScore(string, string)` / `CalculateSaCas9HitScore(string, string)` / `CalculateCfdScore(string, string, string)`: published single-hit models.
- `CrisprDesigner.CalculateMitSpecificityScore(...)` / `CalculateCfdSpecificityScore(...)` (an aggregate-only overload and a genome-scanning overload each): published guide-level aggregates; `CalculateMitSpecificityScoreWithAlternativePams(...)` / `CalculateCfdSpecificityScoreWithAlternativePams(...)`: the same over CRISPOR's wider alternative-PAM search.

### 5.2 Current Behavior

The current implementation enumerates PAM sites first and then compares the guide against each PAM-adjacent target sequence. Mismatch positions are collected explicitly, and exact matches are excluded from the off-target results. `CalculateSpecificityScore(...)` always uses up to 4 mismatches, regardless of the `FindOffTargets(...)` default of 3. Seed mismatches add `5` points and non-seed mismatches add `2` points to the off-target score.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- PAM-constrained candidate enumeration.
- A PAM-proximal seed region with stronger mismatch penalties.
- A specificity score that decreases as off-target penalties accumulate.

**Intentionally simplified:**

- The `FindOffTargets` / `CalculateSpecificityScore` heuristic uses only two positional weights (`5` in the seed region and `2` outside it) and is **unsourced**; **consequence:** it does not represent base-specific or context-specific activity — for a published number call `CalculateMitSpecificityScore` (MIT/Hsu 2013, with Tycko 2018 for SaCas9) or `CalculateCfdSpecificityScore` (Doench 2016 CFD).
- The implementation uses strict PAM matching for the chosen system, one system per scan; **consequence:** off-targets at other PAMs are ignored unless that PAM is modeled as its own system (for example `SpCas9_NAG`). CRISPOR's alternative-PAM ×0.2 MIT downweight (`crispor.py` L1947-1958) is therefore not reachable from the single-system aggregate `CalculateMitSpecificityScore(string, DnaSequence, int, CrisprSystemType)` — it IS implemented, on the multi-PAM entry points `CalculateMitSpecificityScore(string, IEnumerable<PamAnnotatedOffTarget>, CrisprSystemType)` (with the public `FindOffTargetsWithPam`) and `CalculateMitSpecificityScoreWithAlternativePams`, which scans the main PAM plus the alternative PAMs of the reference's `offtargetPams` table that this library offers as systems (`NGG → NAG, NGA` via `SpCas9_NAG` / `SpCas9_NGA`; `NNGRRT → NNGRRN` via `SaCas9_NNGRRN`; `TTTN` has no PAM-after-target system and Cas12a is refused anyway). A site matching both the main and a wider alternative pattern (every `NNGRRT` site is also `NNGRRN`) is one off-target and counted once; the 0.2 weight is gated on a main `NGG`, so SaCas9 `NNGRRN` hits enter unweighted, as in the reference. Only the MIT score is downweighted, as in the reference. The CFD aggregate has its own multi-PAM entry point, `CalculateCfdSpecificityScoreWithAlternativePams` (B08 round 3 G1): `annotateOfftargets` computes `calcCfdScore(guideSeq, otSeq)` for every hit of the wider search, with no 0.2 factor and the same `ALTPAMMINSCORE` pre-filter on the **MIT** hit score (L4104-4122); the single-system `CalculateCfdSpecificityScore` keeps scanning the one PAM (a guide whose only off-target sits on `NAG` scores 100 there, 99.815 with the alternative scan). See Evidence, divergence 2.
- The three specificity scorers still consider base substitutions only; **consequence:** `FindOffTargetsWithBulges` reports gap-containing sites without a score. The reference's SaCas9 model refuses gaps explicitly, and CFD is defined only for a substitution-only 20-nt:20-nt comparison.
- Only the 21-nt SaCas9 coefficient table is transcribed (the reference ships 20-24); **consequence:** `CalculateSaCas9HitScore` refuses other lengths, which no CRISPR system in this library can produce.
- **On-target and repeats (B08 round 4 H4/H5).** `FindOffTargets` drops every exact (0-mismatch) match - its contract is the 1..max-mismatch enumeration. The four CRISPOR-faithful aggregates (`CalculateMitSpecificityScore(string, DnaSequence, ...)`, `CalculateMitSpecificityScoreWithAlternativePams`, `CalculateCfdSpecificityScore(string, DnaSequence, ...)`, `CalculateCfdSpecificityScoreWithAlternativePams`; MCP `crispr_mit_specificity_score` / `crispr_cfd_specificity_score`) follow `annotateOfftargets` (crispor.py L1917-1928): without an interval the **first** exact match in scan order (forward-strand PAMs ascending, then reverse-strand) is the guide's own on-target - CRISPOR's "auto-ontarget" mode - and every further perfect match is scored as a full off-target (MIT hit score 100, CFD its PAM score; a guide with 2 / 5 / 50 / 500 perfect sites scores 50 / 20 / 2 / 0.2); with `onTargetStart`/`onTargetEnd` (0-based, end-exclusive, forward-strand) the on-target is every exact match whose protospacer+PAM footprint lies inside the interval, either strand. **Consequence:** the aggregates assume the first perfect site is the intended one; name the interval when it is not (a genome that lacks the guide's own site loses its first perfect off-target). With an alternative-PAM aggregate the first exact match of the main-PAM scan is the on-target and exact alternative-PAM matches are kept. **Repeats:** `maxOccurrences` (default 60000 = CRISPOR `MAXOCC`) - a guide with more alignments (guide-length windows, either strand, any PAM, within `maxMismatches` substitutions = BWA's `X0 + X1`) scores 0 under every aggregate, as the reference's `repCount` does; its `isRep` flag is the same condition. The reference additionally empties the off-target list for such a guide; the enumerations here have no score to zero and still list their hits. Oracle: CRISPOR's `annotateOfftargets` executed, and BWA 0.7.17 for the count (Evidence, divergences 10-11).

**Intentionally simplified (alternative-PAM cutoff):**

- The reference's alternative-PAM score cutoff (`ALTPAMMINSCORE = 1.0`, `--minAltPamScore`, `crispor.py` L318 / L4108-4122: a hit whose PAM is in the alternative list is dropped when its unweighted `calcHitScore` is `< 1.0`) is the optional `minAlternativePamHitScore` of `CalculateMitSpecificityScoreWithAlternativePams`, **defaulting to CRISPOR's `1.0`** (B08 round 3 G6; the method was created in this campaign to reproduce CRISPOR's wider search, and the default was `0` before — pass `0` for no filtering). As in the reference (alternative list tested first), a SaCas9 cutoff also applies to its `NNGRRT` main-PAM hits. The cutoff must be a number that is not negative: `NaN` or a negative value is an `ArgumentOutOfRangeException` (`minAlternativePamHitScore`; zero = no filtering, `+Infinity` = every alternative-PAM hit filtered out). It is validated with the other arguments, **before** the repeat-guide rule can return 0, in both alternative-PAM aggregates, so whether it is refused never depends on the genome's content (B08 A6-1: round 5 had moved the repeat rule in front of the check, and a repeat guide came back 0 for `-1`). The same holds for the CFD aggregates' `aggregation`: an undefined `CfdAggregation` value (an enum arrives over MCP as any number) is an `ArgumentOutOfRangeException` before the repeat rule, in both CFD genome aggregates (B08 A7-1; it was detected only in the last line of the aggregation, after a repeat guide had returned 0).

**Not implemented:**

- Chromatin accessibility and an experimentally calibrated activity score for bulged alignments. The new scan finds gap-containing candidate sites, but activity prediction for those sites needs external/experimental evidence. All three published sequence-level substitution-only models the reference uses (MIT/Hsu, Tycko 2018 SaCas9, Doench 2016 CFD) and both of its guide-level aggregates are implemented here.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| No off-targets found | Specificity score is `100` | Maximum specificity per implementation |
| `maxMismatches < 0` or `> 5` | Throws `ArgumentOutOfRangeException` | Explicit source guard |
| Guide length mismatches the system | Throws `ArgumentException` | The method enforces system-specific guide length |
| Exact ungapped match | Not returned as an off-target | The no-locus convention excludes it (inside the on-target interval only, when one is given); exact bulged alignments are reported |

### 6.2 Limitations

`CalculateSpecificityScore` is a simplified, unsourced position-weighted sequence model. The published models (`CalculateMitSpecificityScore`, `CalculateCfdSpecificityScore`) are reference-verified but still sequence-only: none of them includes chromatin state or indel/bulge alignments, and all rely on the CRISPR system's explicit PAM definition for site discovery. The substitution-only genome scan is also a naïve sweep over every PAM site, not CRISPOR's BWA-based search, so it is meant for a supplied reference sequence rather than a whole genome. The bulge-aware scan has a larger cost because it evaluates every allowed gap offset; it was not part of the following benchmark. **Measured size envelope (B08 round 3 G14, 2026-10-21):** random ACGT sequence, 20-nt guide, one Debug build (`net10.0`, JIT optimiser disabled, so a Release build is faster), single thread on a 4-vCPU Intel Xeon @ 2.10 GHz cloud VM (16 GB RAM; wall clock, median of 3 runs after one warm-up run). Time grows linearly with the sequence length (about 0.15-0.16 s per Mb per PAM scan):

| Call | 1 Mb | 5 Mb | 10 Mb | 20 Mb |
|---|---|---|---|---|
| `FindOffTargets` (3 mismatches) | 0.13 s | 0.81 s | 1.56 s | 3.13 s |
| `CalculateMitSpecificityScore` (4 mismatches, one PAM) | 0.12 s | 0.77 s | 1.55 s | 3.29 s |
| `CalculateCfdSpecificityScore` (4 mismatches, one PAM) | - | 0.62 s | - | - |
| `CalculateMitSpecificityScoreWithAlternativePams` (NGG + NAG + NGA scans) | - | 2.13 s | - | - |
| `CalculateCfdSpecificityScoreWithAlternativePams` (NGG + NAG + NGA scans) | - | 2.38 s | - | - |

The alternative-PAM aggregates run three scans, hence about three times the single-scan time. These are single data points on one machine, not a guarantee: absolute numbers depend on hardware, build configuration and hit density (random sequence has almost no hits within 4 mismatches; a repeat-rich genome produces more off-targets to score). Whole-genome use (Gb scale, about 150-500 s per scan at this rate) is outside this method's design - use an indexed aligner as CRISPOR does.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
string guide = "ACGTACGTACGTACGTACGT";
var genome = new DnaSequence("ACGT...long genome...ACGT");

var offTargets = CrisprDesigner.FindOffTargets(
    guide,
    genome,
    maxMismatches: 3,
    CrisprSystemType.SpCas9);

// This library's own heuristic (coarse triage only).
double heuristic = CrisprDesigner.CalculateSpecificityScore(
    guide, genome, CrisprSystemType.SpCas9);

// The published aggregates — prefer these.
double mit = CrisprDesigner.CalculateMitSpecificityScore(
    guide, genome, maxMismatches: 4, CrisprSystemType.SpCas9);
double cfd = CrisprDesigner.CalculateCfdSpecificityScore(
    guide, genome, maxMismatches: 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle);
```

## 8. References

1. Hsu PD, Scott DA, Weinstein JA, et al. (2013). DNA targeting specificity of RNA-guided Cas9 nucleases. *Nature Biotechnology*, 31(9):827-832. doi:10.1038/nbt.2647
2. Fu Y, Foden JA, Khayter C, et al. (2013). High-frequency off-target mutagenesis induced by CRISPR-Cas nucleases in human cells. *Nature Biotechnology*, 31(9):822-826. doi:10.1038/nbt.2623
3. Wikipedia. Off-target genome editing. https://en.wikipedia.org/wiki/Off-target_genome_editing
4. Doench JG, Fusi N, Sullender M, et al. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nature Biotechnology*, 34(2):184-191. doi:10.1038/nbt.3437 (CFD)
5. Tycko J, Barrera LA, Huston NC, et al. (2018). Pairwise library screen systematically interrogates *Staphylococcus aureus* Cas9 specificity in human cells. *Nature Communications*, 9:2962. doi:10.1038/s41467-018-05391-2 (SaCas9 model)
6. Concordet J-P, Haeussler M (2018). CRISPOR: intuitive guide selection for CRISPR/Cas9 genome editing experiments and screens. *Nucleic Acids Research*, 46(W1):W242-W245. doi:10.1093/nar/gky354 (reference implementation)
