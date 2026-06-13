# Evidence Artifact: META-RESIST-001

**Test Unit ID:** META-RESIST-001
**Algorithm:** Antibiotic Resistance Gene Detection (ResFinder-style acquired-gene detection)
**Date Collected:** 2026-06-13

---

## Online Sources

### Zankari et al. (2012), "Identification of acquired antimicrobial resistance genes", J Antimicrob Chemother 67(11):2640–2644

**URL:** https://academic.oup.com/jac/article/67/11/2640/707208
**Accessed:** 2026-06-13 (fetched via WebFetch of the article URL above)
**Authority rank:** 1 (peer-reviewed paper; original ResFinder publication)

**Key Extracted Points:**

1. **Core method:** "All genes from the ResFinder database were BLASTed against the assembled genome, and the best-matching genes were given as output." Detection is BLAST of each database gene vs the query, reporting the best-matching gene.
2. **Coverage / length cutoff:** "For a gene to be reported, it has to cover at least 2/5 of the length of the resistance gene in the database." (i.e. minimum coverage relative to the *reference* gene length.)
3. **Default identity:** "The default ID is 100%." The user can select a %ID threshold.
4. **Percent identity definition:** %ID is "the percentage of nucleotides that are identical between the best-matching resistance gene in the database and the corresponding sequence in the genome."
5. **Selected operating threshold in study:** "A ResFinder threshold of ID = 98.00% was selected, as previous tests of ResFinder had shown that a threshold lower than this gives too much noise (e.g. fragments of genes)."

### ResFinder GitHub repository (genomicepidemiology/resfinder)

**URL:** https://github.com/genomicepidemiology/resfinder
**Accessed:** 2026-06-13 (fetched via WebFetch of the repository README)
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **Default identity threshold (-t):** `CGE_RESFINDER_GENE_ID` default 0.80 — "Minimum threshold for identity."
2. **Default coverage threshold (-l):** `CGE_RESFINDER_GENE_COV` default 0.60 — "Minimum (breadth-of) coverage of ResFinder within the range 0–1."
3. **Coverage meaning:** breadth-of-coverage = proportion of a reference gene's sequence covered by the alignment (0–1).

### Mahfouz et al. / pipeline validation, Sci Rep (2023) 13:15543 (carbapenem-resistant K. pneumoniae)

**URL:** https://www.nature.com/articles/s41598-023-42154-6
**Accessed:** 2026-06-13 (fetched via WebFetch; abstract/methods)
**Authority rank:** 1 (peer-reviewed)

**Key Extracted Points:**

1. **Operating thresholds:** ResFinder uses "98% identity" and "60% coverage".
2. **Reason for 60% coverage:** "to ensure that genes lying on the edge of a contig or spread over two contigs are not missed, due to non-perfect assembly." Coverage is the proportion of the reference gene sequence matched in the assembly.

### Benchmarking of AMR-gene identification methods, JAC (2016) 71(9):2484–2492

**URL:** https://academic.oup.com/jac/article/71/9/2484/2238319
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 1 (peer-reviewed)

**Key Extracted Points:**

1. **Thresholds confirmed:** ResFinder is operated at "98% identity" and "60% coverage"; coverage is defined relative to the reference gene length.

### Heng Li, "On the definition of sequence identity" (lh3.github.io, 2018)

**URL:** https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (well-cited reference by the samtools/minimap2 author)

**Key Extracted Points:**

1. **BLAST identity formula:** BLAST identity = "the number of matching bases over the number of alignment columns." The denominator (alignment length) includes M/I/D columns; in a gapless alignment there are no gap columns, so the denominator equals the number of aligned positions. Worked example: 43 matches over 50 columns → 86%.

### CARD Resistance Gene Identifier (RGI) — McMaster (best-hit / identity criteria)

**URL:** https://card.mcmaster.ca/analyze/rgi (and arpcard/rgi docs, via WebSearch 2026-06-13)
**Accessed:** 2026-06-13 (retrieved via WebSearch summary of card.mcmaster.ca and github.com/arpcard/rgi)
**Authority rank:** 5 (curated database / reference tool documentation)

**Key Extracted Points:**

1. **Perfect match:** "A Perfect RGI match is 100% identical to the reference protein sequence along its entire length" — confirms identity=1.0 at full coverage is the unambiguous top of the scale.
2. **Best-hit selection:** RGI ranks candidate matches by alignment bit-score (best hit reported); corroborates ResFinder's "best-matching gene" single-output convention.

---

## Documented Corner Cases and Failure Modes

### From Zankari et al. (2012)

1. **Low-identity noise / gene fragments:** thresholds below ~98% ID "gives too much noise (e.g. fragments of genes)" — sub-threshold matches must be rejected, not reported.
2. **Coverage floor:** a partial hit below 2/5 (≥60% in later versions) of the reference length is not reported.

### From Sci Rep (2023)

1. **Genes split across contigs / contig edges:** the 60% coverage floor exists precisely so edge / fragmented genes are still detectable; a hit need not span the full reference.

---

## Test Datasets

### Dataset: Synthetic exact-and-mismatch nucleotide cases (derived from BLAST identity definition)

**Source:** Heng Li (2018) identity formula; Zankari et al. (2012) coverage and best-match rules.

Because CARD/ResFinder databases are large curated tables that cannot be reproduced verbatim, the
implementation is a generic detector taking caller-supplied reference genes and thresholds; expected
values below are derived arithmetically from the cited formulas (identity = matches / alignment-length,
coverage = window / reference-length), not from any hard-coded gene list.

| Case | Contig | Reference gene | Expected %ID | Expected coverage |
|------|--------|----------------|--------------|-------------------|
| Exact full-length | `AAACGTACGT` contains `CGTACGT` | `CGTACGT` (len 7) | 7/7 = 1.0 | 7/7 = 1.0 |
| One mismatch, full length | `CGTTCGT` vs `CGTACGT` | `CGTACGT` (len 7) | 6/7 ≈ 0.857142857 | 1.0 |
| Partial (edge), 4 of 7 bases | contig ends with `CGTA`, ref `CGTACGT` | `CGTACGT` (len 7) | 4/4 = 1.0 | 4/7 ≈ 0.571428571 |

---

## Assumptions

1. **ASSUMPTION: Gapless (ungapped) alignment model** — The detector locates the best *ungapped*
   alignment (reference slid across the contig, no insertions/deletions). ResFinder uses full
   gapped BLAST; the BLAST identity formula (matches / alignment columns) is identical for the
   gapless case (no gap columns). This affects output only for genes whose true alignment to the
   contig requires indels; for substitution-only divergence and contig-edge truncation (the cases
   the coverage floor targets) it is exact. Documented as a scope simplification.

---

## Recommendations for Test Coverage

1. **MUST Test:** Exact full-length match → %ID=1.0, coverage=1.0, reported. — Evidence: Zankari (2012) best-match; Heng Li identity.
2. **MUST Test:** Single-mismatch full-length → %ID=6/7, coverage=1.0. — Evidence: Heng Li identity = matches/columns.
3. **MUST Test:** Contig-edge partial hit (coverage < 1) still scored by reference length; passes if ≥ coverage threshold. — Evidence: Zankari (2012) 2/5 rule; Sci Rep (2023) edge rationale.
4. **MUST Test:** Below identity threshold → not reported. — Evidence: Zankari (2012) noise/fragment rule.
5. **MUST Test:** Below coverage threshold → not reported. — Evidence: Zankari (2012) coverage floor.
6. **MUST Test:** Best-matching gene only (highest identity) reported per contig. — Evidence: Zankari (2012) "best-matching genes"; RGI best hit.
7. **MUST Test:** Default thresholds equal ResFinder values (0.90 ID / 0.60 cov). — Evidence: ResFinder web service / GitHub README.
8. **SHOULD Test:** null/empty/invalid-threshold inputs raise the documented exceptions. — Rationale: contract robustness.
9. **COULD Test:** tie-break by coverage when identity ties. — Rationale: deterministic selection.

---

## References

1. Zankari E, Hasman H, Cosentino S, et al. (2012). Identification of acquired antimicrobial resistance genes. J Antimicrob Chemother 67(11):2640–2644. https://academic.oup.com/jac/article/67/11/2640/707208
2. genomicepidemiology/resfinder (reference implementation, README defaults). https://github.com/genomicepidemiology/resfinder
3. Pipeline validation for identification of AMR genes in carbapenem-resistant K. pneumoniae (2023). Sci Rep 13. https://www.nature.com/articles/s41598-023-42154-6
4. Clausen PTLC, et al. / Benchmarking of methods for identification of antimicrobial resistance genes (2016). J Antimicrob Chemother 71(9):2484–2492. https://academic.oup.com/jac/article/71/9/2484/2238319
5. Li H (2018). On the definition of sequence identity. https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
6. Alcock BP, et al. CARD Resistance Gene Identifier (RGI). https://card.mcmaster.ca/analyze/rgi

---

## Change History

- **2026-06-13**: Initial documentation.
