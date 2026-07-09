---
type: concept
title: "Methylation context classification (CpG/CHG/CHH)"
tags: [epigenetics, algorithm]
sources:
  - docs/Evidence/EPIGEN-METHYL-001-Evidence.md
  - docs/algorithms/Epigenetics/Methylation_Analysis.md
source_commit: f59358fda234a26aa74779bcfa3844f8fbdae7f8
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: epigen-methyl-001-evidence
      evidence: "Test Unit ID: EPIGEN-METHYL-001 ... Algorithm: Methylation Site Detection, Sequence-Context Classification (CpG/CHG/CHH), and Methylation Profile ... Algorithm Group: Epigenetics"
      confidence: high
      status: current
---

# Methylation context classification (CpG/CHG/CHH)

Classifying each **cytosine** of a DNA sequence into its **methylation sequence context** — **CpG**,
**CHG**, or **CHH** — from the base(s) immediately 3' of it, and enumerating one methylation *site* per
classifiable cytosine. This is the **sixth (final) ingested unit of the Epigenetics family** and a
**sequence-only** classifier: distinct from every sibling because it partitions *non-CpG* cytosines
(the CHG/CHH contexts that [[bisulfite-methylation-calling]] explicitly does **not** call from reads).
Its siblings: [[bisulfite-methylation-calling]] *measures* methylation from bisulfite reads (and shares
this unit's `GenerateMethylationProfile` aggregator), [[cpg-island-detection]] scores CpG *density*
(O/E ratio + Gardiner-Garden islands), [[epigenetic-age-horvath-clock]] scores age from measured
β-values, [[chromatin-state-prediction]] works on histone marks, and
[[differentially-methylated-regions]] compares two samples' β-values. Validated under test unit
**EPIGEN-METHYL-001**; the record is [[epigen-methyl-001-evidence]], [[test-unit-registry]] tracks the
unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

DNA methylation occurs almost exclusively at **CpG** in mammals, but in plant and stem-cell genomes
also at non-CpG cytosines in the **CHG** and **CHH** contexts (Lister et al. 2009). The three contexts
are distinguished purely by the bases 3' of the cytosine (Krueger & Andrews 2011). The unit exposes
three deterministic operations on `EpigeneticsAnalyzer`, all comparing bases as upper-case {A,C,G,T}:

1. `GetMethylationContext(sequence, index)` — classify one cytosine → CpG/CHG/CHH or null.
2. `FindMethylationSites(sequence)` — enumerate every classifiable cytosine site.
3. `GenerateMethylationProfile(sites)` — per-context weighted profile (shared with
   [[bisulfite-methylation-calling]]; detailed there).

## 1. The trinucleotide context rule (IUPAC H = not G)

For a cytosine at 0-based position `i`, with downstream bases `b₁ = s[i+1]` and `b₂ = s[i+2]`:

```
CpG  ⟺  b₁ = G
CHG  ⟺  b₁ ∈ H  and  b₂ = G
CHH  ⟺  b₁ ∈ H  and  b₂ ∈ H
```

where **H = {A, C, T}** is the IUPAC ambiguity code for **"not G"** (Cornish-Bowden 1985). The single
sharp edge is that **a G is never H**: `CGG` classifies its first C as **CpG** (next base is G), *not*
CHG — the context is decided by the exact base, G taking precedence. CpG and CHG are **symmetric** (the
context recurs on the complementary strand); **CHH is asymmetric** and has no symmetric counterpart.
Classification is `O(1)`, case-insensitive, deterministic. Oracle: `CGACAGCAA` → C@0 `CG…`=CpG,
C@3 `CAG`=CHG (H=A, then G), C@6 `CAA`=CHH (H=A, H=A).

## 2. Methylation-site enumeration (sequence-only)

`FindMethylationSites` applies the rule at every position via `GetMethylationContext` (so the two entry
points can never disagree), yielding one `MethylationSite` per classifiable C: `Position` (0-based),
`Type` (CpG/CHG/CHH), `Context` (the up-to-3-base upper-cased window; a terminal CpG carries a 2-base
context), and — because there is **no bisulfite read evidence** — `MethylationLevel = 0` and
`Coverage = 0`. This zero is a **representational default** (the site is *potentially* methylatable),
not a claim that the level is zero (Assumption 1). Real levels arrive from
`CalculateMethylationFromBisulfite` ([[bisulfite-methylation-calling]]) or caller-supplied sites.
`O(n)` single scan; the repository suffix tree is deliberately **not** used (a constant-window lookup
per position is already optimal).

## 3. Per-context weighted profile (Schultz 2012, shared)

`GenerateMethylationProfile` partitions sites by context and computes the **weighted methylation
level** `Σ(level·coverage)/Σ(coverage)` per context (CG/CHG/CHH separately), plus `TotalCpGSites`,
`MethylatedCpGSites`, and `MethylationByPosition`. This is the **same aggregator** documented in
[[bisulfite-methylation-calling#3-per-context-weighted-profile-schultz-2012|bisulfite-methylation-calling]]
— weighted (read-pooled) fraction, not the unweighted mean; when a context's total coverage is 0
(sequence-only sites) it **falls back to the unweighted mean** so such sites are not dropped. Worked
oracle: sites (0.8, cov 10) and (0.2, cov 10) → weighted CpG = (8+2)/(10+10) = **0.5** (equals the
mean here only because coverage is equal). The `MethylatedCpGSites` count uses a **descriptive 0.5
fractional cutoff** (Assumption 2); Schultz (2012) recommends a binomial test, so the count is not a
statistically tested call — the continuous levels are unaffected.

## Invariants and edge cases

- **INV:** a cytosine is CpG ⟺ the next base is G; CHG/CHH require the H position(s) ∈ {A,C,T} (G is
  never H); reported positions are 0-based and index a cytosine.
- **INV:** each per-site level ∈ [0,1]; per-context profile level = Σ(level·coverage)/Σcoverage, equal
  to the mean of fractions only under equal coverage.
- **Terminal `CG`** is still CpG (two bases suffice); a **terminal `C` + one H** is unclassified
  (CHG/CHH need the third base).
- A **non-ACGT** base in an H or third position (`CNG`) → undetermined context → unclassified.
- Null/empty sequence → empty enumeration (or null context); index < 0 or ≥ length → null; empty site
  list → all-zero profile.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for context classification, site
enumeration, and the shared weighted profile. **Single-strand only:** it does not pair symmetric
CpG/CHG contexts across strands, performs no read alignment, and treats degenerate bases as
unclassifiable rather than expanding them. Methylation levels are only as meaningful as the coverage
supplied — sequence-only input yields zero-level placeholder sites. **Biological realism** (Lister
2009): non-CG (CHG+CHH) methylation is a stem-cell/plant phenomenon (~25% in H1 ES cells) and is
essentially absent in differentiated somatic cells (IMR90 = 99.98% CpG). For read-derived non-CpG
calling or full-genome WGBS, use dedicated pipelines (Bismark). No source contradictions —
Cornish-Bowden 1985, Krueger & Andrews 2011, Lister 2009, and Schultz 2012 are mutually consistent.
</content>
