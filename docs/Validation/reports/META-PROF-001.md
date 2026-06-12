# Validation Report: META-PROF-001 — Taxonomic Profiling / Relative Abundance

- **Validated:** 2026-06-12   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.GenerateTaxonomicProfile(IEnumerable<TaxonomicClassification>)`
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs:237-288`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_TaxonomicProfile_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Relative species abundance"** (https://en.wikipedia.org/wiki/Relative_species_abundance):
  "Relative abundance is the percent composition of an organism of a particular kind relative to
  the total number of organisms in the area." → relative abundance of a taxon = (count of that
  taxon) / (total count of all taxa), the values being proportional parts of the whole and thus
  summing to 1.0 (100%). Confirms the defined computation.
- **Wikipedia — "Metagenomics"** / **Segata et al. (2012), MetaPhlAn** (Nature Methods,
  doi:10.1038/nmeth.2066): MetaPhlAn estimates *organismal relative abundances* per taxonomic
  rank; higher ranks are sums of their leaf (species) abundances. Profiles are normalized over
  the *assigned/classified* fraction (reads/markers that map to a clade); "unclassified" reads
  form a residual that is not part of the per-clade relative-abundance denominator. This matches
  the spec's explicit statement (M3/M4): unclassified reads are excluded from both
  `ClassifiedReads` and the abundance denominator.

### Formula check (matches cited sources exactly)
| Metric | Formula | Source |
|--------|---------|--------|
| Relative abundance | `count(taxon) / Σ count(all classified taxa)` | Wikipedia, MetaPhlAn |
| Shannon diversity   | `H = −Σ pᵢ·ln(pᵢ)` (natural log → nats) | Shannon (1948) |
| Simpson concentration | `λ = Σ pᵢ²` (original Simpson, not 1−λ or 1/λ) | Simpson (1949) |

### Denominator semantics (spec states this explicitly)
Denominator = **classified reads only** (reads whose Kingdom is non-empty and ≠ "Unclassified").
Unclassified / empty-kingdom reads are counted in `TotalReads` but excluded from `ClassifiedReads`
and from every per-rank abundance denominator. This is the MetaPhlAn convention and is asserted
by spec rows M3, M4, M6.

### Edge-case semantics (all defined & sourced)
- Empty input → `TotalReads=0`, `ClassifiedReads=0`, empty maps, Shannon=0, Simpson=0
  (empty-sum convention: Σ over zero terms = 0).
- Single classified read → abundance = 1/1 = 1.0; Shannon=0, Simpson=1.0.
- All unclassified → `ClassifiedReads=0`, empty abundance maps.
- Missing rank values → empty strings filtered from rank maps.

### Independent cross-check (hand computation)
- Worked example from prompt: A=50, B=30, C=20, total=100 → 0.50, 0.30, 0.20; **Σ = 1.0** ✓.
  Code computes `value/total` per taxon with `total = classifiedReads = 100`, reproducing this exactly.
- 3 equal taxa → each 1/3 ≈ 0.3333, Σ = 1.0 ✓ (test `EqualCounts_ProduceEqualAbundances`).
- 3 uniform species → H = −3·(⅓·ln ⅓) = ln 3 ≈ 1.0986 ✓.
- Species counts [2,1,1] → p=[0.5,0.25,0.25], λ = 0.25+0.0625+0.0625 = 0.375 ✓.
- Skew p=[0.9,0.1] → H = −(0.9·ln0.9 + 0.1·ln0.1) ≈ 0.325; λ = 0.81+0.01 = 0.82 ✓.

### Findings / divergences
None. Description is biologically and mathematically correct.

## Stage B — Implementation

### Code path reviewed
`GenerateTaxonomicProfile` (`MetagenomicsAnalyzer.cs:237-280`), helper `IncrementCount` (282-288),
and the shared `CalculateShannonIndex` (327-339) / `CalculateSimpsonIndex` (341-353).

### Formula realised correctly? (evidence)
- `classified = classList.Where(c => c.Kingdom != "Unclassified" && !IsNullOrEmpty(c.Kingdom))`
  → exactly the classified-only filter (line 242); `classifiedReads = classified.Count` (243).
- `double total = classifiedReads > 0 ? classifiedReads : 1;` (258) — denominator is classified reads;
  the `:1` guard only fires when there are 0 classified reads, in which case the count dictionaries
  are empty so no division actually occurs. No div-by-zero, no spurious entries.
- Per-rank abundance = `count / total` (260-266); empty rank keys filtered for phylum/genus/species
  (and excluded at increment time via `IncrementCount`). Matches `count/Σcount`.
- Shannon/Simpson are computed from species-level abundance **values**; both helpers internally
  re-normalize by `sum` of the passed values (lines 329/343), so they yield correct entropy /
  concentration even if species abundances sum to <1 (e.g., some classified reads lack species),
  and correctly reduce to the textbook p=value/Σ form. Natural log (`Math.Log`) = nats; Simpson is
  λ=Σpᵢ² returned directly.

### Cross-verification table recomputed vs code (all pass)
| Case | Expected | Test | Result |
|------|----------|------|--------|
| Empty input | TotalReads=0, Classified=0, Shannon=0, Simpson=0 | M1 | ✓ |
| Single read | abundance 1.0 all ranks | M2 | ✓ |
| Σ kingdom & species = 1.0 | ≈1.0 | M8 | ✓ |
| 3 equal taxa | each 1/3 | M14 | ✓ |
| Unclassified excluded from count | Total=3, Classified=1 | M3 | ✓ |
| Unclassified excluded from denom | Bacteria=1.0 (not 0.5) | M4 | ✓ |
| 3 uniform species | Shannon=ln3 | M9 | ✓ |
| [2,1,1] | Simpson=0.375 | M10 | ✓ |
| single species | Shannon=0, Simpson=1.0 | M11/M12 | ✓ |
| 4 uniform | Shannon=ln4, Simpson=0.25 | S1/S2 | ✓ |
| skew [0.9,0.1] | Shannon≈0.325, Simpson=0.82 | S4 | ✓ |
| all 4 ranks | each Σ=1.0 | S3 | ✓ |

### Variant/delegate consistency
Single canonical method; the diversity helpers are shared with `CalculateAlphaDiversity` and
behave identically (re-normalize by sum). No divergent variant.

### Numerical robustness
`double` arithmetic; denominator guarded; no overflow on stated O(n) input. Within(0.001) tolerances
are appropriate for the floating-point values asserted.

### Test quality audit
18 spec tests + 2 edge cases (20 total) assert **exact sourced values** (ln 3, ln 4, 0.375, 0.82,
1/3, 1.0, 0/0) rather than tautologies; edge cases (empty, single, all-unclassified, missing ranks,
empty kingdom) are covered. Deterministic.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS** — relative abundance = count/classified-total summing to 1.0; unclassified
  excluded from the denominator; Shannon (nats) and Simpson (λ=Σpᵢ²) match primary sources.
- **Stage B: PASS** — code faithfully realizes the validated formulas; worked example A/B/C=50/30/20
  reproduces 0.50/0.30/0.20 (Σ=1.0); all 20 tests pass.
- **State: CLEAN** — no defects found; no code changes required.
- Tests: `~TaxonomicProfile` filter = 20 passed; full `Seqeron.Genomics.Tests` = 4484 passed,
  0 failed (baseline preserved).
