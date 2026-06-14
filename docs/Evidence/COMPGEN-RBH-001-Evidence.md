# Evidence Artifact: COMPGEN-RBH-001

**Test Unit ID:** COMPGEN-RBH-001
**Algorithm:** Reciprocal Best Hits (RBH / bidirectional best hits) for ortholog identification
**Date Collected:** 2026-06-14

---

## Online Sources

### Moreno-Hagelsieb & Latimer (2008) — operational RBH definition and thresholds

**URL:** https://academic.oup.com/bioinformatics/article/24/3/319/252715 — DOI https://doi.org/10.1093/bioinformatics/btm585
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed, *Bioinformatics*)
**Retrieved how:** WebSearch "Moreno-Hagelsieb Latimer 2008 reciprocal best hits orthologs BLAST bit-score coverage definition", then WebFetch of the OUP article URL. The four quotations below are copied verbatim from the fetched article text.

**Key Extracted Points:**

1. **RBH definition:** *"two genes residing in two different genomes are deemed orthologs if their protein products find each other as the best hit in the opposite genome."* → orthology = mutual (reciprocal) best hit; a one-directional best hit is not enough.
2. **Best-hit ranking:** hits are *"sorted the BLASTP hits from highest to lowest bit-score, then, if the bit-scores were identical, from smallest to highest E-values."* → best hit = maximum similarity score, with a deterministic tie-break.
3. **Coverage threshold:** *"we also required coverage of at least 50% of any of the protein sequences in the alignments."* → qualifying gate of ≥ 50% coverage (the fetched text says 50%, not 60%; a search-engine summary that said "60%" was contradicted by the article body and is rejected).
4. **E-value threshold:** *"maximum E-value threshold of 1×10−6."* → a significance gate; in an alignment-free implementation it maps to a minimum-similarity gate.

### Tatusov, Koonin & Lipman (1997) — COG, symmetrical (reciprocal) best hits

**URL:** https://www.ncbi.nlm.nih.gov/books/NBK21090/ (NCBI Handbook, "The Clusters of Orthologous Groups (COGs) Database", describing the Tatusov et al. 1997 method); citation: *Science* 278(5338):631–637, https://doi.org/10.1126/science.278.5338.631 (PubMed https://pubmed.ncbi.nlm.nih.gov/9381173/)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed *Science*; method described by the official NCBI Handbook)
**Retrieved how:** WebSearch "Tatusov Koonin Lipman 1997 genomic perspective protein families COG symmetrical best hit BLAST"; the *Science* full-text PDF (fire.biol.wwu.edu) returned HTTP 404, so the construction procedure was taken from the WebFetch of the NCBI Handbook page NBK21090, copied below.

**Key Extracted Points:**

1. **Best hits between genomes:** proteins from different genomes that are *"more similar to each other than they are to any other proteins from the same genomes are most likely to form an orthologous set."* → the genome-specific best hit (BeT) in each direction is the building block.
2. **Mutual best-hit triangles:** COGs are built by *"Detect triangles of mutually consistent, genome-specific best hits (BeTs)"* — i.e. best-hit relationships that hold reciprocally/consistently across genomes; the pairwise special case is the reciprocal best hit.
3. **COG composition:** *"Each COG is a group of three or more proteins that are inferred to be orthologs, i.e., they are direct evolutionary counterparts."*

---

## Documented Corner Cases and Failure Modes

### From Moreno-Hagelsieb & Latimer (2008)

1. **Tie in best hit:** equal bit-scores are broken by smaller E-value; an implementation must apply a deterministic tie-break so the "best hit" is unique and the matching is order-independent.
2. **Coverage filter:** a high-scoring local match over a short region must be rejected unless it covers ≥ 50% of a sequence (avoids spurious orthologs from shared domains).

### From the RBH definition (Tatusov 1997; Moreno-Hagelsieb 2008)

1. **Non-reciprocity:** if gene A's best hit is B but B's best hit is C ≠ A, then A–B is **not** an ortholog pair (one-directional best hit ≠ ortholog). This is the defect class the unit must guard against.
2. **No qualifying hit:** a gene with no above-threshold hit in the other genome yields no pair.
3. **Empty genome:** with 0 sequence-bearing genes in either genome there can be no between-genome pair.

---

## Test Datasets

### Dataset: Reciprocity worked example (RBH definition, Moreno-Hagelsieb 2008)

**Source:** Moreno-Hagelsieb & Latimer (2008), RBH definition (best hit each way). Sequences chosen so the 5-mer-content similarity ranking is unambiguous.

| Gene (genome) | Sequence | Best hit in other genome |
|---------------|----------|--------------------------|
| a1 (G1) | `ACGTACGTACGTAC` | b1 |
| a2 (G1) | `TTTTGGGGCCCCAA` | b2 |
| b1 (G2) | `ACGTACGTACGTAC` | a1 |
| b2 (G2) | `TTTTGGGGCCCCAA` | a2 |

Expected RBH: {a1↔b1, a2↔b2} (each pair is identical sequences ⇒ maximal similarity each way; identity = coverage = 1.0).

### Dataset: Non-reciprocity (one-directional best hit is NOT an ortholog)

**Source:** Tatusov 1997 symmetrical-best-hit requirement; Moreno-Hagelsieb 2008.

| Gene (genome) | Sequence |
|---------------|----------|
| a1 (G1) | `ACGTACGTACGTAC` |
| b1 (G2) | `ACGTACGTACGTAC` (identical to a1) |
| b2 (G2) | `ACGTACGTACGTACGTACGT` (a1's superstring; shares all of a1's 5-mers) |

a1's best hit is b1 (Jaccard 1.0) over b2 (Jaccard < 1.0). b1's best hit is a1; b2's best hit is a1, but a1's best hit is b1, not b2. Expected RBH: {a1↔b1} only; b2 excluded (not reciprocal). RBH count = 1.

### Dataset: Coverage / minimum-identity gate

**Source:** Moreno-Hagelsieb (2008) ≥ 50% coverage and significance gate.

A pair whose similarity score falls below `minIdentity`, or whose coverage falls below `minCoverage`, does not qualify as a hit and therefore yields no RBH pair even when it is the mutual top candidate. Setting `minIdentity = 1.0` (above any non-identical Jaccard) rejects all non-identical pairs.

---

## Assumptions

1. **ASSUMPTION: Alignment-free similarity score as the best-hit ranking metric.** Moreno-Hagelsieb (2008) ranks best hits by BLAST bit-score. The `ComparativeGenomics` class has no alignment-based bit-score available (the `Seqeron.Genomics.Alignment` project is not referenced), so candidates are ranked by a 5-mer Jaccard similarity, a monotone alignment-free similarity (cf. Mash, Ondov et al. 2016). This affects **which** pair wins only for near-identical ties; the RBH reciprocity rule, the deterministic tie-break, the ≥ 50% coverage gate, and the minimum-similarity gate (the correctness-critical parts) are source-backed. Identical sequences score 1.0, so the metric is order-preserving for the datasets above. The ≥ 50% coverage gate is mapped to "shared k-mers ≥ 50% of the smaller k-mer set."

---

## Recommendations for Test Coverage

1. **MUST Test:** Mutual best hits are returned as RBH pairs (reciprocity). — Evidence: Moreno-Hagelsieb & Latimer (2008) RBH definition.
2. **MUST Test:** A one-directional best hit (non-reciprocal) is NOT returned; count = 1 for the superstring dataset. — Evidence: Tatusov (1997) symmetrical best hits; Moreno-Hagelsieb (2008).
3. **MUST Test:** Returned pair carries the actual hit identity and coverage (identical sequences ⇒ identity = coverage = 1.0), not a hardcoded coverage. — Evidence: Moreno-Hagelsieb (2008) best-hit metrics.
4. **MUST Test:** Sub-threshold pairs are rejected by the minimum-identity / coverage gate. — Evidence: Moreno-Hagelsieb (2008) ≥ 50% coverage and significance gate.
5. **MUST Test:** Empty genome → no RBH pairs. — Evidence: definition requires a between-genome pair.
6. **MUST Test:** Null inputs throw `ArgumentNullException`. — Rationale: repository contract (sibling methods).
7. **SHOULD Test:** RBH yields a matching — each gene appears at most once across the returned pairs. — Rationale: best-hit reciprocity is a partial matching.
8. **SHOULD Test:** Genes with empty/absent sequence are skipped. — Rationale: similarity undefined without sequence.
9. **COULD Test:** Determinism — same input yields identical pair set regardless of input order. — Rationale: deterministic tie-break / order-independence.

---

## References

1. Moreno-Hagelsieb G, Latimer K. 2008. Choosing BLAST options for better detection of orthologs as reciprocal best hits. *Bioinformatics* 24(3):319–324. https://doi.org/10.1093/bioinformatics/btm585 — https://academic.oup.com/bioinformatics/article/24/3/319/252715
2. Tatusov RL, Koonin EV, Lipman DJ. 1997. A genomic perspective on protein families. *Science* 278(5338):631–637. https://doi.org/10.1126/science.278.5338.631 — https://pubmed.ncbi.nlm.nih.gov/9381173/ ; method described in NCBI Handbook, "The Clusters of Orthologous Groups (COGs) Database", https://www.ncbi.nlm.nih.gov/books/NBK21090/
3. Ondov BD, et al. 2016. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol.* 17:132. https://doi.org/10.1186/s13059-016-0997-x (basis for alignment-free k-mer similarity ranking; see Assumption 1)

---

## Change History

- **2026-06-14**: Initial documentation.
