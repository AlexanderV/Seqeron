---
type: concept
title: "Average Nucleotide Identity (ANI)"
tags: [comparative-genomics, algorithm]
mcp_tools:
  - calculate_ani
sources:
  - docs/Evidence/COMPGEN-ANI-001-Evidence.md
  - docs/algorithms/Comparative_Genomics/Average_Nucleotide_Identity.md
  - docs/Validation/reports/COMPGEN-ANI-001.md
source_commit: 205b259dc3168dfda72a89caf5103f39ac5e1ce9
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: compgen-ani-001-evidence
      evidence: "Test Unit ID: COMPGEN-ANI-001 ... Algorithm: Average Nucleotide Identity (ANI), ANIb definition (Goris et al. 2007)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:synteny-and-rearrangement-detection
      source: compgen-ani-001-evidence
      evidence: "ANI (COMPGEN-ANI-001) and synteny (COMPGEN-SYNTENY) are sibling Comparative-genomics units; ANI is the genome-similarity metric, synteny the gene-order-conservation view of the same two-genome comparison"
      confidence: medium
      status: current
---

# Average Nucleotide Identity (ANI)

**ANI** is the genome-wide sequence-similarity metric that operationalizes the prokaryotic
species boundary. It is the first ingested unit of the **Comparative-genomics** family
(`COMPGEN-*`), a sibling of the shared synteny anchor
[[synteny-and-rearrangement-detection]]: where synteny asks whether *gene order* is
conserved between two genomes, ANI asks how *nucleotide-identical* their conserved regions
are. Validated under test unit **COMPGEN-ANI-001**; the validation record is
[[compgen-ani-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern. The independent two-stage
re-validation [[compgen-ani-001-report]] closed the verdict at **Stage A/B PASS · CLEAN**
(20/20 tests, zero code change), confirming the gapped-placement and reciprocal-ANI additions
resolved the prior report's two PASS-WITH-NOTES items.

Seqeron exposes `CalculateAni` (single-direction) and `CalculateReciprocalAni`
(symmetric, both-directions mean). A third COMPGEN sibling,
[[conserved-gene-clusters-common-intervals]], compares the two genomes at the
gene-*content* level (sets of genes contiguous in both) rather than nucleotide identity.
The end-to-end pipeline [[genome-comparison-core-dispensable]] orchestrates these views into a
single core/dispensable gene partition plus a syntenic-gene fraction. The
[[dot-plot-word-match]] sibling keeps the same two-sequence comparison as a **visual** word-match
matrix (diagonals = similarity) rather than reducing it to a single identity number.

## The ANIb algorithm (Goris et al. 2007)

ANI here means the **ANIb** definition — the BLAST-based procedure of Goris et al. 2007,
which reproduces the ~1 kb fragmentation of the classic DNA–DNA hybridization (DDH)
experiment in silico:

1. **Fragment** the query genome into consecutive **1020 nt** fragments (the cut-off chosen
   to mirror the ~1 kb sheared DNA of DDH). A trailing partial fragment is ignored.
2. **Search** each fragment against the whole reference genome (BLASTN in the canonical
   pipeline); keep the **best match** per fragment.
3. **Filter**: a fragment qualifies only when its best match shows **> 30 % identity**
   (recalculated over the entire fragment length) **over an alignable region of ≥ 70 %** of
   the fragment length. Non-conserved fragments are **discarded**, not scored as zero.
4. **Average**: ANI = the **mean per-fragment identity** across the qualifying fragments.

Two definitional subtleties, both source-critical:

- **Identity is recalculated over the whole query fragment**, not just the locally aligned
  sub-region: per-fragment identity = identical aligned bases / query-fragment length
  (pyani `ani_pid = ani_alnids / qlen`). A short but perfect local hit still scores low.
- **Coverage is likewise over the query length**: `ani_coverage = ani_alnlen / qlen`
  (ungapped aligned length / fragment length), gated at ≥ 0.70.

## Gapped placement (indels)

The canonical placement is **gapped** — pyani's BLASTN runs with `-xdrop_gap_final 150` and
computes `ani_alnlen = blast_alnlen - blast_gaps`. An **ungapped** placement underestimates
identity across indel-containing homologous regions; the gapped path (`gapped: true`)
recovers it. Worked oracle: query fragment `AAAACCCC` vs reference `AAAATCCCC` (one inserted
`T`) → ungapped best window 7/8 = 0.875, but gapped `AAAA-CCCC` gives 8 identical columns →
identity 1.0, coverage 1.0.

## Reciprocal (symmetric) ANI

ANI(query → reference) need not equal ANI(reference → query), because only the **query** is
fragmented (pyani notes "non-symmetrical result matrices"). Goris perform **reverse
searching** "to provide reciprocal values"; the symmetric ANIb value is the **mean of the
two directions**. `CalculateReciprocalAni` implements this mean — order-independent by
construction, so ANI(A,B) = ANI(B,A). The single-direction value is well defined on its own.

## Species boundary

ANI is the modern replacement for wet-lab DDH: **ANI ≈ 95 %** corresponds to the 70 % DDH
standard for delineating a prokaryotic species (Goris et al. 2007; Konstantinidis & Tiedje
2005 place the equivalence at ≈ 94–95 %).

## Parameters and invariants

| Parameter | Meaning | Default |
|-----------|---------|---------|
| `fragmentLength` | consecutive fragment size | 1020 nt |
| `minIdentity` | per-fragment identity cut-off | 0.30 |
| `minAlignableFraction` | per-fragment coverage cut-off | 0.70 |
| `gapped` | gapped Smith-Waterman placement | true (recovers indels) |

- Result is a fraction in **[0, 1]** (invariant of the identity definition).
- Identical genomes → ANI = 1.0 (every fragment is a perfect substring).
- Number of fragments = ⌊len / fragmentLength⌋; trailing partial fragment ignored.
- **Null / empty inputs → 0**; **non-positive `fragmentLength` → `ArgumentOutOfRangeException`**.

## Corner cases

- **Non-conserved fragments discarded** — a fragment whose best match is ≤ 30 % identity or
  < 70 % alignable is excluded from the mean, so ANI is computed over conserved fragments
  only, not the whole genome. (Oracle: query `AAAACGTC` vs ref `AAAAAAAA`, fragLen 4 →
  `CGTC` scores 0 and is excluded → ANI = 1.0 from the surviving `AAAA` fragment.)
- **Reference shorter than the fragment** → no full-length placement → alignable fraction 0
  < 0.70 → no qualifying fragment → ANI = 0.
- **Query shorter than the fragment** → no fragment fits → ANI = 0.

## Reference tools and engine choice

The definition traces to **Goris et al. 2007** (IJSEM, the ANIb method), **Konstantinidis &
Tiedje 2005** (PNAS, the species-threshold framing), and the **pyani** ANIb reference
implementation (Pritchard et al. 2016). Seqeron's one documented **decision** (not a
correctness gap): the gapped path uses the library's own Smith-Waterman aligner
(`SequenceAligner.LocalAlign`, BLAST DNA scoring) rather than the NCBI BLASTN engine.
Smith-Waterman is full dynamic programming — more sensitive than BLAST's heuristic seeding,
not less — and identity/coverage use the identical recalculated-over-fragment definition;
numeric ANI on real genomes may differ slightly from NCBI-BLASTN pipelines because the
alignment engine differs, but the indel-handling behaviour is correct.

**Sharp edge:** [[ani-is-directional-use-reciprocal]] — `calculate_ani` is **directional** (ANI(A,B) != ANI(B,A)); use reciprocal ANI for a symmetric value.
