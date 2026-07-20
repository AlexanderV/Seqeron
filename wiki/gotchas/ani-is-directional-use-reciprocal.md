---
type: gotcha
title: "calculate_ani is directional — ANI(A,B) ≠ ANI(B,A); use reciprocal ANI for a symmetric value"
tags: [comparative-genomics, gotcha]
mcp_tools:
  - calculate_ani
sources:
  - docs/algorithms/Comparative_Genomics/Average_Nucleotide_Identity.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# calculate_ani is directional — ANI(A,B) ≠ ANI(B,A); use reciprocal ANI for a symmetric value

**The trap.** `calculate_ani` fragments **only the query** genome and searches each fragment against
the reference, so **ANI(query → reference) need not equal ANI(reference → query)** (pyani notes
"non-symmetrical result matrices"). Swapping the two arguments can return a different number.

**Why it bites.**

- **Asymmetric pairwise matrices.** Building an ANI matrix by calling `calculate_ani` once per pair
  yields an asymmetric matrix, and the **~95 % species cutoff** you apply depends on which genome
  was the query — two directions can straddle the boundary and flip a same-species call.
- **Non-conserved fragments are discarded, not scored 0.** ANI is a mean over **conserved**
  fragments only; a genome that shares one strong region reads high even if most of it is unrelated.
- **A reference shorter than one fragment (1020 nt) → ANI = 0**, because no fragment gets a
  full-length placement — that is "not measurable", not "0 % similar".

**What to rely on instead.** Use **`CalculateReciprocalAni`** (the mean of both directions) for an
order-independent value whenever symmetry matters — species delimitation, clustering, a distance
matrix. Single-direction ANI is well defined, but record which genome was the query. Full model:
[[average-nucleotide-identity]]. Pipeline view: [[comparative-genomics-pipeline-silent-traps]].
