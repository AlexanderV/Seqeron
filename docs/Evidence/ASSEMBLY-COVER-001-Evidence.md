# Evidence Artifact: ASSEMBLY-COVER-001

**Test Unit ID:** ASSEMBLY-COVER-001
**Algorithm:** Coverage (Depth) Calculation — per-base sequencing depth over a reference
**Date Collected:** 2026-06-13

---

## Online Sources

### Illumina — "Sequencing Coverage for NGS Experiments"

**URL:** https://sapac.illumina.com/science/technology/next-generation-sequencing/plan-experiments/coverage.html
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 2 (official vendor specification / standard usage)

**Key Extracted Points:**

1. **Coverage definition:** Retrieved text: "Next-generation sequencing (NGS) coverage describes the average number of reads that align to, or 'cover,' known reference bases." Coverage / sequencing depth = average number of reads covering each reference base.
2. **Average-coverage formula:** The page presents the Lander-Waterman formula **C = LN / G**, where **C** = coverage, **L** = read length, **N** = number of reads, **G** = haploid genome length. It is "a method for computing genome coverage."

### Daniel E. Cook — "Calculate Depth and Breadth of Coverage From a bam File"

**URL:** https://www.danielecook.com/calculate-depth-and-breadth-of-coverage-from-a-bam-file/
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (reference-implementation / established workflow over samtools output)

**Key Extracted Points:**

1. **Depth of coverage:** Retrieved text states depth and coverage "refer to the same thing, the average number of reads aligned to an individual base." Computed as **Sum of Depths / genome size** (or contig length per chromosome).
2. **Per-base depth:** measured via `samtools depth` as "the number of reads aligned over a given base."
3. **Breadth of coverage:** "the percentage of the genome with aligned reads to it," computed as **Bases Mapped / genome size**, a fraction representing the proportion of bases with at least one aligned read.

### Metagenomics Wiki — "SAMtools: get breadth of coverage"

**URL:** https://www.metagenomics.wiki/tools/samtools/breadth-of-coverage
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (established bioinformatics tooling documentation)

**Key Extracted Points:**

1. **Per-base depth:** "the number of reads mapping to a specific reference position."
2. **Breadth of coverage:** "percentage/fraction of reference bases covered by at least one read" (example: 32,876 covered bases / 45,678 total = 0.719×).
3. **Average depth:** "total bases mapped / reference length" = "sum of per-base depths / genome size."

### Predicting the Number of Bases to Attain Sufficient Coverage (Daley et al.) — PMC7398442

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7398442/
**Accessed:** 2026-06-13 (fetched via WebFetch, after 301 redirect from ncbi.nlm.nih.gov)
**Authority rank:** 1 (peer-reviewed, restates the Lander-Waterman model)

**Key Extracted Points:**

1. **Lander-Waterman assumption:** Retrieved text: "the basic statistical assumption of the Lander/Waterman model is that reads are generated uniformly at random from the genome. Under this assumption, for each site in the genome, the number of reads that cover this site can be approximated by a Poisson distribution, where the Poisson rate λ is estimated by average sequencing depth."
2. **Sufficient coverage:** "a site as having sufficient coverage if the site is covered by at least r reads."

### Lander & Waterman (1988) Poisson model (via search-retrieved restatement)

**URL:** Search results for "Lander Waterman 1988 Genomics genomic mapping by fingerprinting random clones Poisson coverage" and "probability that a base is not covered e^-c gaps" (WebSearch; primary: Genomics 3:231-239).
**Accessed:** 2026-06-13 (WebSearch summary; primary not directly fetchable as full text in this session)
**Authority rank:** 1 (primary paper) / 4 (the retrieved restatement)

**Key Extracted Points:**

1. **Gap probability:** Retrieved text: "The probability that any given base in a sequence is not sequenced can be calculated by the equation **P = e^−m**, where m is the fold coverage." Worked numbers: 1× → e^−1 = 0.37 (≈37% uncovered); 5× → e^−5 = 0.0067.
2. **Breadth complement:** "The complementary probability (**1 − e^−c**) represents the probability that a base is covered at least once."

---

## Documented Corner Cases and Failure Modes

### From samtools / Metagenomics Wiki / Daniel Cook

1. **Position with zero overlapping reads:** depth at that reference position is 0 (it does not contribute to breadth-covered bases).
2. **Read extending past reference end:** only the overlapping portion of a read contributes to per-base depth (clipping at the reference boundary).
3. **Uncovered genome → breadth 0, average depth 0:** when no read aligns, every per-base depth is 0; average depth = 0 and breadth = 0.

### From Lander-Waterman / PMC7398442

1. **Highly variable coverage:** the uniform-Poisson expectation (mean ≈ C) holds only under uniform random read placement; observed per-base depth can deviate. This is a modeling caveat, not an arithmetic rule — the per-base depth array itself is exact regardless.

---

## Test Datasets

### Dataset: Hand-constructed exact-placement example (per-base depth)

**Source:** Definition from samtools / Metagenomics Wiki (per-base depth = number of reads covering each position); deterministic by construction.

Reference `ACGTTGCAAT` (length 10). Reads are distinct 5-mer substrings, each matching uniquely, so placement is unambiguous; per-base depth = count of reads spanning each position.

| Parameter | Value |
|-----------|-------|
| Reference | `ACGTTGCAAT` (len 10) |
| Read 1 | `ACGTT` aligns at pos 0 → covers [0,5) |
| Read 2 | `TTGCA` aligns at pos 3 → covers [3,8) |
| Read 3 | `GCAAT` aligns at pos 5 → covers [5,10) |
| Expected depth array | `[1,1,1,2,2,2,2,2,1,1]` |
| Total bases mapped (Σ depth) | 15 |
| Average depth (Σ/G) | 15/10 = 1.5 |
| Breadth (covered ≥1× / G) | 10/10 = 1.0 |

### Dataset: Lander-Waterman gap/breadth sanity (Poisson)

**Source:** Lander & Waterman (1988); P(uncovered) = e^−c; breadth = 1 − e^−c.

| Parameter | Value |
|-----------|-------|
| c = 1× | P(uncovered) = e^−1 ≈ 0.3679; breadth ≈ 0.6321 |
| c = 5× | P(uncovered) = e^−5 ≈ 0.006738 |

(Used as a property/derivation check, not asserted against the deterministic per-base depth array which is exact.)

---

## Assumptions

1. **ASSUMPTION: Mapping model for read placement** — The unit signature is `CalculateCoverage(reference, reads, minOverlap)`; the sources define depth given an *alignment*, but do not prescribe how an aligner places a read. The repository uses an ungapped best-match scan requiring ≥ `minOverlap` matching characters (`FindBestAlignment`). This placement rule is implementation-level (it determines *where* a read maps, not the depth-counting arithmetic). The depth-counting arithmetic (per-base = number of placed reads spanning the position; clip at reference end) is fully source-defined. Tests therefore use exact-match reads where placement is unambiguous, isolating the source-defined counting rule.

---

## Recommendations for Test Coverage

1. **MUST Test:** Per-base depth array equals the count of reads spanning each position for unambiguous exact-match placements. — Evidence: samtools / Metagenomics Wiki ("number of reads mapping to a specific reference position").
2. **MUST Test:** Read extending past reference end contributes only its overlapping portion (boundary clipping). — Evidence: per-base counting at reference boundary (implementation contract derived from depth definition).
3. **MUST Test:** Unmatched read (below `minOverlap`) contributes 0 to all positions. — Evidence: a read that does not align does not add depth (samtools depth semantics).
4. **MUST Test:** Empty reads list → all-zero depth array of reference length. — Evidence: no aligned reads → every per-base depth 0 (Daniel Cook / Metagenomics Wiki).
5. **SHOULD Test:** Average depth = Σ(per-base depth) / reference length on the worked dataset (1.5). — Rationale: Daniel Cook "Sum of Depths / genome size."
6. **SHOULD Test:** Breadth = (# positions with depth ≥ 1) / reference length on the worked dataset (1.0). — Rationale: Metagenomics Wiki breadth definition.
7. **COULD Test:** Case-insensitive matching (lowercase read vs uppercase reference still maps). — Rationale: implementation normalizes case; verifies depth is placement-driven, not case-driven.

---

## References

1. Illumina, Inc. (n.d.). Sequencing Coverage for NGS Experiments. https://sapac.illumina.com/science/technology/next-generation-sequencing/plan-experiments/coverage.html
2. Cook, D.E. (n.d.). Calculate Depth and Breadth of Coverage From a bam File. https://www.danielecook.com/calculate-depth-and-breadth-of-coverage-from-a-bam-file/
3. Metagenomics Wiki. (n.d.). SAMtools: get breadth of coverage. https://www.metagenomics.wiki/tools/samtools/breadth-of-coverage
4. Daley, T. et al. (2020). Predicting the Number of Bases to Attain Sufficient Coverage in High-Throughput Sequencing Experiments. PMC7398442. https://pmc.ncbi.nlm.nih.gov/articles/PMC7398442/
5. Lander, E.S., Waterman, M.S. (1988). Genomic mapping by fingerprinting random clones: a mathematical analysis. Genomics 3(2):231-239. https://doi.org/10.1016/0888-7543(88)90007-9 (primary; full text not fetched this session — Poisson gap/breadth formulas taken from retrieved restatements in refs 1 & 4 and WebSearch summaries).

---

## Change History

- **2026-06-13**: Initial documentation.
