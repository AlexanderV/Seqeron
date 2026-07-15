---
type: source
title: "Evidence: META-PROF-001 (taxonomic profile — community abundance profiling)"
tags: [validation, metagenomics]
doc_path: docs/Evidence/META-PROF-001-Evidence.md
sources:
  - docs/Evidence/META-PROF-001-Evidence.md
source_commit: 02f28f4a5999dc16f47afa3db37eca0cb1eda2ee
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: META-PROF-001

The validation-evidence artifact for test unit **META-PROF-001** — **taxonomic profile
generation**, the aggregation of per-read taxonomic classifications into per-taxon
relative-abundance estimates by `MetagenomicsAnalyzer.GenerateTaxonomicProfile`. The
community-abundance **profiling** step downstream of the per-read
[[taxonomic-classification|META-CLASS-001 classification]] unit — it consumes an
`IEnumerable<TaxonomicClassification>`. One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the method is synthesized in its own
concept [[taxonomic-profile]]; [[test-unit-registry]] tracks the unit. See
`docs/Evidence/META-PROF-001-Evidence.md`.

## What this file records

- **Online sources (mutually consistent):**
  - **Wikipedia — Metagenomics** — defines taxonomic profiling; "taxonomic profiling aggregates
    classification results into abundance distributions"; MetaPhlAn / AMPHORA marker-based relative
    abundance; metagenomes compared by normalized species/genus abundance.
  - **Wikipedia — Relative abundance distribution** — abundance normalized as fractions summing to
    1.0 (SAD = species abundance distribution).
  - **MetaPhlAn documentation** — canonical profiler; profiles at kingdom → phylum → … → species;
    relative abundance as percentages/fractions.
  - **Segata et al. (2012), Nature Methods (DOI 10.1038/nmeth.2066)** — MetaPhlAn primary reference
    (unique clade-specific marker genes).

## Algorithm (from the Evidence file)

Input `IEnumerable<TaxonomicClassification>` → `TaxonomicProfile` record (impl
`MetagenomicsAnalyzer.cs` lines 229–280) with:

- **Abundance maps at four ranks** (kingdom, phylum, genus, species): relative abundance =
  `count(taxon) / Σcount(all classified taxa)`.
- **Diversity metrics** `Shannon = −Σpᵢln(pᵢ)` (natural log = nats) and `Simpson = Σpᵢ²`
  (concentration index λ), computed from **species-level** abundances.
- **Read counts** `TotalReads`, `ClassifiedReads`.

Behaviour notes: excludes `Unclassified` from abundance denominators; filters empty rank strings;
abundances normalized to sum 1.0 before diversity computation.

## Source-verified invariants and oracles

**Invariants:** Σ(abundance) ≈ 1.0 per rank; `ClassifiedReads ≤ TotalReads`;
`ClassifiedReads = Σ(counts at any rank)`; `Shannon ≥ 0`; `0 ≤ Simpson ≤ 1.0`.

**Documented test cases / edge cases:** empty classification list → `TotalReads = 0`,
`ClassifiedReads = 0`, empty abundances; all Unclassified → `ClassifiedReads = 0`, no entries;
single taxon → abundance 1.0, Shannon 0, Simpson 1.0; uniform over `n` → Shannon `ln(n)`,
Simpson `1/n`; skewed → low Shannon, high Simpson; missing genus/species → empty strings filtered.

**Strengthened oracles (2026-03-09 coverage pass):** `Shannon = ln(3)` for 3 uniform species (M9);
`Simpson = 0.375` for the non-uniform `[2,1,1]` distribution (M10); `TotalReads = 3 /
ClassifiedReads = 2` (M7); all 4 ranks sum to 1.0 for fully-populated reads (S3). Result: 18 tests,
all covered.

## Deviations

**Verified design decisions (all mathematical facts, not assumptions):** Shannon uses natural log
(nats, ecology standard, `Math.Log`); Simpson is the concentration index `λ = Σpᵢ²` (distinct from
Gini-Simpson `1−λ` and inverse Simpson `1/λ`); empty input yields `Shannon = 0, Simpson = 0`
(empty-sum convention). No literature deviations; no source contradictions.
