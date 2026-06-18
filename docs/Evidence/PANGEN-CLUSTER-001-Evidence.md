# Evidence Artifact: PANGEN-CLUSTER-001

**Test Unit ID:** PANGEN-CLUSTER-001
**Algorithm:** Gene Clustering (homolog / ortholog grouping by sequence identity)
**Date Collected:** 2026-06-13

---

## Online Sources

### CD-HIT User's Guide (Li & Godzik) — identity definition and `-c` option

**URL:** https://vcru.wisc.edu/simonlab/bioinformatics/programs/cd-hit/cdhit-user-guide.pdf
**Accessed:** 2026-06-13 (PDF fetched with WebFetch, text extracted locally via zlib stream decode)
**Authority rank:** 3 (reference implementation documentation)

**Key Extracted Points:**

1. **`-c` cutoff and default:** Verbatim: "`-c` sequence identity threshold, default 0.9 this is the default cd-hit's 'global sequence identity' calculated as: number of identical amino acids in alignment divided by the full length of the shorter sequence".
2. **Global identity (default, `-G 1`):** identity = (identical residues in alignment) / (full length of the **shorter** sequence). This is the default mode.
3. **Local identity (`-G 0`):** Verbatim: "use local sequence identity, calculated as: number of identical amino acids in alignment divided by the length of the alignment". (Not the default; recorded for contrast.)
4. **Cluster structure:** Verbatim: "CD-HIT clusters proteins into clusters that meet a user-defined similarity threshold, usually a sequence identity. Each cluster has one representative sequence."
5. **Per-member identity reported:** in the `.clstr` output, "`%` is the identity between this sequence and the representative" — confirming members are compared to the cluster representative, not all-pairs.

### CD-HIT Algorithm wiki — greedy incremental clustering procedure

**URL:** https://github.com/weizhongli/cdhit/wiki/1.-Algorithm
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 3 (reference implementation documentation)

**Key Extracted Points:**

1. **Greedy incremental procedure (verbatim):** "CD-HIT is a greedy incremental clustering approach. The basic CD-HIT algorithm sorts the input sequences from long to short, and processes them sequentially from the longest to the shortest. The first sequence is automatically classified as the first cluster representative sequence. Then each query sequence of the remaining sequences is compared to the representative sequences found before it, and is classified as redundant or representative based on whether it is similar to one of the existing representative sequences."
2. **First-match (fast) assignment:** "In default manner (fast mode), a query is grouped into the first representative without comparing to other representatives." — i.e. the query joins the first representative meeting the threshold (greedy, not best-hit).

### Roary: rapid large-scale prokaryote pan genome analysis (Page et al., 2015) — pan-genome clustering context and default identity

**URL:** https://academic.oup.com/bioinformatics/article/31/22/3691/240757
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 1 (peer-reviewed, Bioinformatics)

**Key Extracted Points:**

1. **Workflow:** coding regions are converted to proteins, "filtered to remove partial sequences and iteratively pre-clustered with CD-HIT," then an "all-against-all comparison ... with BLASTP on the reduced sequences with a user defined percentage sequence identity (default 95%)," then clustered with MCL.
2. **Default identity for pan-genome ortholog grouping:** 95% (`-i` option), confirming sequence-identity thresholding (not k-mer similarity) is the basis for grouping homologous genes into clusters.
3. **Cluster = ortholog/homolog group:** clustering produces homologous gene groups; conserved gene neighbourhood is later used to split paralogs from true orthologs (out of scope for a greedy identity clusterer).

### EMBOSS needle manual — percent-identity convention (corroboration)

**URL:** https://galaxy-iuc.github.io/emboss-5.0-docs/needle.html
**Accessed:** 2026-06-13 (WebSearch summary of the EMBOSS needle manual)
**Authority rank:** 3 (EMBOSS reference tool documentation)

**Key Extracted Points:**

1. **Percent identity convention:** identity is "the number of corresponding positions in the alignment showing an identical ... nucleotide in both sequences" over the alignment length; the numerator is identical aligned positions. This corroborates the numerator (count of identical positions) used by CD-HIT; the denominator differs by convention (alignment length vs shorter-sequence length). CD-HIT's shorter-sequence denominator is the one tied to the clustering model used here.

---

## Documented Corner Cases and Failure Modes

### From CD-HIT User's Guide / Algorithm wiki

1. **Representative is the longest sequence:** because input is sorted long→short and the first unassigned sequence becomes a representative, the cluster representative is the longest member.
2. **Greedy first-match assignment:** a sequence joins the first representative that meets the threshold; with overlapping representatives the assignment depends on processing order (longest-first), so the result is deterministic for a fixed input.
3. **Threshold is inclusive:** a sequence is grouped when identity meets/exceeds the cutoff (`>=`), so `idThreshold = 1.0` clusters only exact-identity (over the shorter length) sequences.

### From general percent-identity definitions (EMBOSS)

1. **Empty vs non-empty:** percent identity counts identical positions; with no residues to align there are no identical positions. Two empty sequences have zero differences over zero positions (defined here as identical); one empty + one non-empty has no shared residues → 0.

---

## Test Datasets

### Dataset: Hand-derived greedy clustering (CD-HIT model)

**Source:** Derived from the CD-HIT greedy incremental algorithm (Li & Godzik, 2006; CD-HIT Algorithm wiki).

Sequences (length-sorted long→short), threshold = 0.8, ungapped global identity = identical positions / shorter length:

| Gene | Sequence | Length |
|------|----------|--------|
| R (rep) | `AAAAAAAAAA` (10) | 10 |
| Q1 | `AAAAAAAAAT` (10, 9/10 = 0.9) | 10 |
| Q2 | `AAAAAAAAAAAA` (12; vs R 10/10 over shorter=10 = 1.0) | 12 |
| Q3 | `CCCCCCCCCC` (10; vs R 0/10 = 0.0) | 10 |

Expected with threshold 0.8: Q2 (longest, becomes rep), then R joins Q2 (1.0), Q1 joins (vs Q2 rep: 9/10=0.9 ≥ 0.8), Q3 starts a new cluster → **2 clusters**: {Q2,R,Q1} and {Q3}.

### Dataset: Global identity worked values (CD-HIT shorter-length denominator)

**Source:** CD-HIT global identity formula (Li & Godzik, 2006).

| seq1 | seq2 | identical positions (over shorter) | shorter len | identity |
|------|------|-----------------------------------|-------------|----------|
| `ATGCATGC` | `ATGCATGC` | 8 | 8 | 1.0 |
| `ATGCATGC` | `ATGCATGG` | 7 | 8 | 0.875 |
| `ATGCATGC` | `ATGCATGCAAAA` | 8 | 8 | 1.0 |
| `ATGC` | `CGTA` | 0 | 4 | 0.0 |

---

## Assumptions

1. **ASSUMPTION: Ungapped alignment** — CD-HIT computes identity over an *alignment* (with banded short-word filtering). This unit computes identity over an ungapped positional comparison of the shared prefix, divided by the shorter length. This matches CD-HIT exactly for sequences differing only by substitutions or by a length difference (no internal indels). For sequences requiring internal gaps to align, identity may be underestimated. Justification: implementing CD-HIT's full banded aligner is out of scope for a single unit; the denominator (shorter length) and numerator (identical positions) conventions are taken verbatim from the source. Recorded as a deviation in the algorithm doc §5.3 / §5.4.
2. **ASSUMPTION: Paralog splitting not performed** — Roary splits homolog groups into orthologs using gene neighbourhood; this clusterer produces homolog groups only (no synteny step). Justification: the canonical method `ClusterGenes` is a sequence-identity clusterer, matching CD-HIT pre-clustering, not the full Roary pipeline.

---

## Recommendations for Test Coverage

1. **MUST Test:** Global identity = identical positions / shorter length for identical, substitution-differing, length-differing, and disjoint sequences — Evidence: CD-HIT User's Guide `-c`/`-G 1`.
2. **MUST Test:** Greedy clustering groups sequences ≥ threshold into one cluster and separates sequences below threshold — Evidence: CD-HIT Algorithm wiki greedy procedure.
3. **MUST Test:** Threshold 1.0 clusters only exact-identity sequences; lowering the threshold merges near-identical sequences — Evidence: CD-HIT `-c` inclusive cutoff.
4. **MUST Test:** Representative is the longest member; identical members across genomes give GenomeCount = number of distinct genomes — Evidence: CD-HIT long→short sort.
5. **SHOULD Test:** Empty/null genomes → no clusters; null inner gene list skipped — Rationale: documented input-validation contract.
6. **SHOULD Test:** Singleton cluster AverageIdentity = 1.0 — Rationale: a sequence is 100% identical to itself.
7. **COULD Test:** Determinism across repeated calls — Rationale: greedy order is fixed by stable long→short sort.

---

## References

1. Li W, Godzik A. (2006). Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics 22(13):1658–1659. https://doi.org/10.1093/bioinformatics/btl158
2. CD-HIT User's Guide (Li lab). CD-HIT documentation. https://vcru.wisc.edu/simonlab/bioinformatics/programs/cd-hit/cdhit-user-guide.pdf
3. CD-HIT Algorithm wiki. https://github.com/weizhongli/cdhit/wiki/1.-Algorithm
4. Page AJ, Cummins CA, Hunt M, et al. (2015). Roary: rapid large-scale prokaryote pan genome analysis. Bioinformatics 31(22):3691–3693. https://doi.org/10.1093/bioinformatics/btv421
5. EMBOSS needle manual (percent-identity convention). https://galaxy-iuc.github.io/emboss-5.0-docs/needle.html

---

## Change History

- **2026-06-13**: Initial documentation.
