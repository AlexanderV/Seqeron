---
type: concept
title: "Allele & genotype frequencies (MAF, MAF filtering)"
tags: [population-genetics, algorithm]
mcp_tools:
  - allele_frequencies
  - filter_variants_by_maf
sources:
  - docs/Evidence/POP-FREQ-001-Evidence.md
source_commit: fec2c72b4f77c252586394fe43424909b13d98d6
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pop-freq-001-evidence
      evidence: "Test Unit ID: POP-FREQ-001 ... Algorithm: Allele Frequency Calculation (Major/Minor Frequencies, MAF, Filtering)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:ancestry-estimation-admixture
      source: pop-freq-001-evidence
      evidence: "Ancestry estimation (POP-ANCESTRY-001) consumes reference-panel allele frequencies F; allele/genotype frequency estimation is the foundational primitive that produces such per-population allele frequencies from genotype counts."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:genetic-diversity-statistics
      source: pop-freq-001-evidence
      evidence: "Per-site allele frequencies p_i feed the gene-diversity / heterozygosity term H = 1 − Σ p_i² used by POP-DIV-001; allele-frequency estimation is the shared upstream primitive."
      confidence: high
      status: current
---

# Allele & genotype frequencies (MAF, MAF filtering)

The **foundational population-genetics primitive**: turn per-locus **genotype counts** into
**allele frequencies**, derive the **minor allele frequency (MAF)**, and **filter variants by MAF**.
This is a population-genetics `POP-*` unit (**POP-FREQ-001**) and the numeric substrate the rest of
the family builds on — the reference allele frequencies **F** that
[[ancestry-estimation-admixture]] treats as fixed, and the per-site frequencies `p_i` behind the
gene-diversity/heterozygosity term of [[genetic-diversity-statistics]], and the per-population
allele frequencies that [[population-differentiation-fst]] compares to measure Fst, are exactly
these quantities. Validated under test unit **POP-FREQ-001**; the literature-traced record is
[[pop-freq-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## Allele frequency from genotype counts

For a **biallelic** locus with major allele A and minor allele B, given diploid genotype counts
`n_AA`, `n_AB`, `n_BB` (Wikipedia *Allele frequency* / *Genotype frequency*):

```
total_alleles = 2·(n_AA + n_AB + n_BB)
major_alleles = 2·n_AA + n_AB
minor_alleles = 2·n_BB + n_AB
p = major_alleles / total_alleles          # freq(A)
q = minor_alleles / total_alleles          # freq(B)
```

Equivalently `p = f(AA) + ½f(AB)`, `q = f(BB) + ½f(AB)`. **Invariant `p + q = 1`** (biallelic),
with `0 ≤ p,q ≤ 1` and `major_alleles + minor_alleles = total_alleles`. Worked oracle (Wikipedia
four-o'clock example): 49 AA, 42 Aa, 9 aa → `freq(a) = (42 + 2·9)/200 = 0.30`,
`freq(A) = (2·49 + 42)/200 = 0.70`.

## MAF from a genotype vector

Under the standard **VCF/PLINK dosage encoding** (0 = homozygous reference, 1 = heterozygous,
2 = homozygous alternate), the alt-allele count per individual **is** the genotype value:

```
alt_freq = Σ genotypes / (2 · n_individuals)
MAF      = min(alt_freq, 1 − alt_freq)
```

**Invariant `0 ≤ MAF ≤ 0.5`** with **symmetry** `MAF = min(f, 1−f)` — MAF is which of the two
alleles is rarer, not which is "alt". So `alt_freq = 0.3 → MAF = 0.3` but `alt_freq = 0.7 → MAF = 0.3`
too (ref is then the minor allele). Monomorphic loci (all-0 or all-2 genotypes) → **MAF = 0**;
perfect 50/50 → **MAF = 0.5** (the maximum).

## MAF filtering

Keep variant `v` iff `minMAF ≤ MAF(v) ≤ maxMAF`. The thresholds are conventional, not derived:
the **HapMap** project targeted SNPs with MAF ≥ 0.05, and the **common vs rare** boundary is
drawn at **MAF = 0.05** (rare variants, MAF < 0.05, are enriched in coding regions). A lower gate
of MAF < 0.01 is a frequent QC filter for very rare / likely-genotyping-error variants.

## Invariants and edge cases

- **Simplex / accounting:** `p + q = 1`; `major_alleles + minor_alleles = total_alleles`.
- **MAF bounds:** `0 ≤ MAF ≤ 0.5`, symmetric under ref/alt swap; MAF = 0 is a valid fixed
  (monomorphic) allele; MAF = 0.5 is the balanced maximum.
- **Empty / zero input:** zero samples → frequencies `(0, 0)` handled gracefully; empty genotype
  vector → MAF = 0; empty filter input → empty result.
- **Negative counts rejected:** genotype counts are non-negative by definition, so a negative count
  throws `ArgumentOutOfRangeException` (input validation, not silent undefined behaviour).
- **Threshold consistency:** at an exact MAF threshold, include/exclude decisions must be applied
  consistently (the filter is an inclusive `[minMAF, maxMAF]` band).

## Scope

Faithful implementation of the textbook biallelic frequency formulae (exact match to Wikipedia
*Allele frequency* / *Minor allele frequency* / *Genotype frequency*, Gillespie 2004). It is a
counting/normalization primitive: it does **not** test Hardy–Weinberg equilibrium (that is the
sibling unit [[hardy-weinberg-equilibrium-test]], POP-HW-001, which consumes these frequencies),
handle multiallelic (>2) sites, phase haplotypes, or impute missing genotypes. The allele and
two-locus haplotype frequencies it produces are also the inputs to [[linkage-disequilibrium]]
(POP-LD-001), which measures non-random association *between* loci. No source contradictions —
the algorithm is fully determined by the sources (Open Questions: none).
