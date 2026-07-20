---
type: gotcha
title: "find_snps/find_indels is alignment-based, not a pileup genotype caller — and indels aren't left-normalized"
tags: [variants, gotcha]
mcp_tools:
  - find_snps
  - find_indels
  - titv_ratio
sources:
  - docs/algorithms/Variants/Variant_Detection.md
  - docs/algorithms/Variants/Indel_Detection.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# find_snps/find_indels is alignment-based, not a pileup genotype caller — and indels aren't left-normalized

**The trap.** These call variants from a **single query-vs-reference global alignment**, not from a
read pileup. There is **no depth model, no base-quality weighting, and no diploid genotype
assignment** — the output is the set of differences between two sequences. Separately, reported
indels are **not left-aligned / parsimony-normalized** (Tan 2015): each indel is reported at the
column `GlobalAlign` happened to produce.

**Why it bites.**

- **It is not a genotyper.** Feeding a consensus vs reference gives sequence differences, not
  genotypes — do not expect GT/DP/GQ, allele fractions, or heterozygous calls from read data.
- **Indel positions can disagree with a normalized VCF.** In low-complexity / tandem-repeat
  regions the alignment is non-unique, so the reported **position** may differ from the canonical
  left-aligned convention (e.g. a tandem-repeat `AATGA→A` left-shifts to `GATGA→G` — **same count
  and type, different coordinate**). Comparing positions against a normalized VCF without
  normalizing both sides will spuriously mismatch. This affects position only, **not** variant
  counts or types.
- The VCF **padding base** and **1-based POS** conventions apply only to the serialized
  `ToVcfLines` output; the in-memory `Variant` uses `"-"` and a **0-based** position.

**What to rely on instead.** Left-align + parsimony-normalize both sides before comparing
coordinates across callers; for read-based genotyping use a pileup caller. Full model:
[[germline-variant-calling-snp-indel]]. Pipeline view: [[comparative-genomics-pipeline-silent-traps]].
