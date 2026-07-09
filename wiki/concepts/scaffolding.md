---
type: concept
title: "Scaffolding (contig ordering and N-gap linking)"
tags: [assembly, algorithm]
sources:
  - docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md
  - docs/algorithms/Extended_Assembly/Scaffolding.md
source_commit: 9ce49bade5c11e63eebbf8c06dd642662321d5a2
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: assembly-scaffold-001-evidence
      evidence: "Test Unit ID: ASSEMBLY-SCAFFOLD-001 ... Scaffolding (joining ordered contigs into scaffolds with N-gaps)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:contig-merge-overlap-collapse
      source: assembly-scaffold-001-evidence
      evidence: "It is possible that the distance estimate is negative, indicating that the two contigs should in fact overlap. If such an overlap is indeed found in the contig overlap graph, the two contigs are merged — the negative-gap case is where scaffolding hands off to overlap merging"
      confidence: medium
      status: current
---

# Scaffolding (contig ordering and N-gap linking)

**Scaffolding** joins an ordered series of contigs into a **scaffold** — a longer sequence in
which adjacent contigs are separated by **gaps of estimated length, represented as runs of the
character `N`**. It is the step *after* contig construction ([[de-bruijn-graph-assembly|DBG]] /
[[overlap-layout-consensus-assembly|OLC]] / [[contig-merge-overlap-collapse|merge]]): the contigs
are already built and ordered/oriented (here supplied pre-ordered as links); scaffolding lays them
out on a common coordinate frame with sized gaps. This is the anchor for the assembly **SCAFFOLD**
family, validated under test unit **ASSEMBLY-SCAFFOLD-001**. The literature-traced validation record
is [[assembly-scaffold-001-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern.

## The scaffold-construction rule (source-traced)

Traced verbatim to Jackman et al. (ABySS 2.0, *Genome Research* 2017) and corroborated by the
NCBI AGP Specification v2.1, Sahlin et al. (2012) and Pop et al. (Bambus, 2004):

- **Scaffold = concatenation of ordered contigs interspersed with `N`-runs.** "The sequences of the
  vertices in a path are concatenated, interspersed with gaps represented by a run of the character
  N, whose length corresponds to the estimate of the distance between those two contigs."
- **Gap length = distance estimate.** A positive gap estimate `g` between two contigs emits exactly
  `g` gap characters. Scaffold length = Σ|contig| + Σ gap.
- **Distance estimation is upstream.** A maximum-likelihood estimator derives the inter-contig
  distance from paired-read alignments; **this unit consumes a supplied gap estimate** and does not
  compute it (cf. Sahlin et al. 2012, whose whole subject is improving that estimator).
- **Contigs form a path; each appears once.** A scaffold is a path of *distinct* contigs (ABySS /
  Bambus path model); a link to an already-placed contig is skipped, and an unlinked contig becomes
  its own single-contig scaffold. Bambus, the canonical greedy scaffolder, "joins together contigs
  with the most links first" — here links are supplied pre-ordered, so ordering is not re-derived.

## Non-positive gap → the AGP unknown-size default (100 N)

The gap estimate can be **zero or negative**, and a negative estimate is a real, expected input
class — Sahlin et al.: "The negative gap case frequently occurs since a de Bruijn-based assembler
splits its contigs at a given node ... that leaves an overlap (negative gap) of one k-mer length."
A negative estimate means the two contigs should in fact **overlap** (and if the overlap is *found*
in the overlap graph, ABySS merges them — see [[contig-merge-overlap-collapse]]).

This unit does **not** perform overlap resolution, so a non-positive estimate is emitted as a gap of
**unknown size**. The NCBI AGP Specification v2.1 fixes the placeholder:

- "Gap lengths must be positive. Negative gaps and gap lines with zero length are not valid."
- "For negative gaps, or gaps of unknown size, use `U` as the component_type and **100** as the gap
  size, since 100 is the GenBank/EMBL/DDBJ standard for gaps of unknown size."

So a **zero or negative** estimate emits exactly **100 `N`** — the constant itself is source-backed
(the INSDC unknown-gap length); only the *decision to fall back to it* rather than resolve the
overlap is the scoping assumption (see below). AGP distinguishes gap component type `N`
(specified size) from `U` (unknown size); the negative/zero case is the `U` (=100) branch.

## Worked oracles (published)

| Inputs | Gap handling | Expected scaffold | Length |
|--------|--------------|-------------------|--------|
| contigs `["ACGT","TTGG","CCAA"]`, links `[(0,1,3),(1,2,2)]` | positive `g` = `g`×`N` | `ACGT`+`NNN`+`TTGG`+`NN`+`CCAA` = `ACGTNNNTTGGNNCCAA` | 17 = 4+3+4+2+4 |
| contigs `["AAAA","TTTT"]`, links `[(0,1,-5)]` | negative → unknown-size default | `AAAA`+(`N`×100)+`TTTT` | 108 = 4+100+4 |

Single positive-gap input → **one** scaffold. A custom gap character is used verbatim in place of `N`
(the fill character is parameterized — the source only fixes the *default* as `N`).

## Assumptions (from the artifact)

- **Unresolved-overlap placeholder uses the AGP unknown-gap length (100).** ABySS *merges* contigs
  when a negative estimate's overlap is *found*; this unit does no overlap resolution, so a
  non-positive estimate is emitted as an unknown-size gap of **100 `N`**. The 100 constant is
  source-backed (GenBank/EMBL/DDBJ standard, NCBI AGP v2.1); the *scoping* decision (fall back
  rather than resolve) is the assumption. **No numeric value is invented.**

## Relation to the other assembly formulations

Scaffolding sits **downstream** of contig assembly and is orthogonal to the overlap-vs-k-mer
contrast between [[overlap-layout-consensus-assembly|OLC]] and [[de-bruijn-graph-assembly|DBG]]: it
takes finished contigs and *lays them out with sized gaps* rather than reconstructing sequence. Its
one point of contact with the rest of the pipeline is the **negative-gap → overlap** case, where it
hands off to the [[contig-merge-overlap-collapse|suffix–prefix merge]] primitive (if the overlap can
be resolved). The per-column [[consensus-sequence|consensus]] and [[coverage-depth-calculation|
coverage]] steps operate on contigs, not on the inter-contig `N`-gaps that scaffolding introduces.

No contradictions among the sources — ABySS 2.0, the NCBI AGP spec, Sahlin et al. and Bambus all
give the same "ordered contigs + sized `N`-gaps" scaffold model; the AGP 100-N unknown-size default
and the ABySS negative-gap = overlap rule are complementary (they cover the same non-positive case
from the file-format and the assembler side respectively).
