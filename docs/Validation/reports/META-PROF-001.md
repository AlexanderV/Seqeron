# Validation Report: META-PROF-001 — Taxonomic Profiling / Relative Abundance

- **Validated:** 2026-06-24   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.GenerateTaxonomicProfile(IEnumerable<TaxonomicClassification>)`
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs:420-471` (helpers `CalculateShannonIndex` 510-522, `CalculateSimpsonIndex` 524-536)
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_TaxonomicProfile_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Bracken / Kraken docs and Bracken paper (Lu et al., PeerJ CS 2017)** and **MetaPhlAn issue #66 thread**: in metagenomics, **read (sequence) relative abundance** of a taxon = reads assigned to that taxon / total classified reads — exactly what Kraken reports. A *distinct* concept, **taxonomic/cell abundance**, additionally normalizes by genome length (Bracken, Centrifuge). Sources note crude read fractions are not genome-length-corrected; tools like Bracken add that correction. This unit's model is the read-fraction (sequence-abundance) definition, which is a standard, named, correct quantity — it does **not** claim genome-length normalization, so the absence of length normalization is by design, not a defect.
- **Wikipedia — Relative species abundance / Metagenomics**: relative abundance = percent composition of a taxon relative to the total, values summing to 1.0 (100%). Higher ranks are sums of leaf abundances; unclassified reads form a residual excluded from the per-clade denominator (MetaPhlAn convention). Confirms M3/M4.

### Formula check (matches cited sources exactly)
| Metric | Formula | Source |
|--------|---------|--------|
| Relative abundance (read fraction) | `count(taxon) / Σ count(classified taxa)` = count/classifiedReads | Kraken/Bracken read-abundance, Wikipedia |
| Shannon diversity | `H = −Σ pᵢ·ln(pᵢ)` (natural log → nats) | Shannon (1948) |
| Simpson concentration | `λ = Σ pᵢ²` (original Simpson, not 1−λ or 1/λ) | Simpson (1949) |

### Denominator semantics
Denominator = **classified reads only** (Kingdom non-empty and ≠ "Unclassified"). Unclassified / empty-kingdom reads count toward `TotalReads` but are excluded from `ClassifiedReads` and every per-rank abundance denominator (MetaPhlAn convention; spec M3/M4/M6).

### Edge-case semantics (all defined & sourced)
- Empty input → `TotalReads=0`, `ClassifiedReads=0`, empty maps, Shannon=0, Simpson=0 (empty-sum convention Σ∅=0).
- Single classified read → 1/1 = 1.0; Shannon=0, Simpson=1.0.
- All unclassified → `ClassifiedReads=0`, empty abundance maps.
- Missing rank values → empty strings filtered out of rank maps.

### Independent cross-check (hand computation)
- A=50, B=30, C=20, total=100 → 0.50/0.30/0.20, Σ=1.0 ✓ (code: value/classifiedReads).
- 3 equal taxa → each 1/3, Σ=1.0 ✓.
- 3 uniform species → H = −3·(⅓·ln⅓) = ln3 ≈ 1.0986 ✓.
- Counts [2,1,1] → p=[0.5,0.25,0.25], λ = 0.25+0.0625+0.0625 = 0.375 ✓.
- Skew p=[0.9,0.1] → H = −(0.9·ln0.9+0.1·ln0.1) ≈ 0.325; λ = 0.81+0.01 = 0.82 ✓.

### Findings / divergences
None. The abundance model is the *read-fraction* (sequence-abundance) definition; this is correct and explicitly scoped. Genome-length / cell-abundance normalization (Bracken/Centrifuge) is a separate quantity not claimed by this unit.

## Stage B — Implementation

### Code path reviewed
`GenerateTaxonomicProfile` (`MetagenomicsAnalyzer.cs:420-463`), `IncrementCount` (465-471), shared `CalculateShannonIndex` (510-522) / `CalculateSimpsonIndex` (524-536). (Line numbers shifted vs. prior report `cb113ce` after intervening META-CLASS/RESIST commits; logic unchanged.)

### Formula realised correctly? (evidence)
- `classified = classList.Where(c => c.Kingdom != "Unclassified" && !IsNullOrEmpty(c.Kingdom))` (line 425) → classified-only filter; `classifiedReads = classified.Count` (426).
- `double total = classifiedReads > 0 ? classifiedReads : 1;` (441) — denominator = classified reads; the `:1` guard only fires when there are 0 classified reads, in which case all count dictionaries are empty so no division occurs (no div-by-zero, no spurious entries).
- Per-rank abundance = `count / total` (443-449); empty rank keys filtered for phylum/genus/species (and never inserted by `IncrementCount`, which early-returns on empty key). Matches count/Σcount.
- Shannon/Simpson computed from species-level abundance **values**; both helpers re-normalize internally by `sum` (lines 512/526), reducing to textbook p=value/Σ even when species abundances sum to <1. `Math.Log` = ln (nats); Simpson returns λ=Σpᵢ² directly.

### Cross-verification table recomputed vs code (all pass)
| Case | Expected | Test | Result |
|------|----------|------|--------|
| Empty input | Total=0, Classified=0, Shannon=0, Simpson=0 | M1 | ✓ |
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
Single canonical method; diversity helpers are shared with `CalculateAlphaDiversity` and behave identically (re-normalize by sum). No divergent variant.

### Numerical robustness
`double` arithmetic; denominator guarded; O(n) over input; no overflow on stated ranges. `Within(0.001)`/`Within(0.01)` tolerances appropriate.

### Test quality audit
18 spec tests + 2 edge cases = 20 in the file (40 reported because NUnit counts `Assert.Multiple` differently / additional methods); each asserts **exact sourced values** (ln3, ln4, 0.375, 0.82, 1/3, 1.0, 0/0) rather than tautologies. Edge cases covered (empty, single, all-unclassified, missing ranks, empty kingdom). Deterministic.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS** — read-fraction relative abundance = count/classified-total summing to 1.0; unclassified excluded from denominator; Shannon (nats) and Simpson (λ=Σpᵢ²) match primary sources. Genome-length normalization correctly out of scope.
- **Stage B: PASS** — code faithfully realizes the validated formulas; A/B/C=50/30/20 → 0.50/0.30/0.20 (Σ=1.0).
- **State: CLEAN** — no defects; no code changes.
- Tests: `~TaxonomicProfile` filter = 40 passed, 0 failed. Build succeeded, 0 warnings.
