# Evidence Artifact: COMPGEN-ORTHO-001

**Test Unit ID:** COMPGEN-ORTHO-001
**Algorithm:** Ortholog identification by Reciprocal Best Hits (RBH); paralog (in-paralog) identification by within-genome best hits
**Date Collected:** 2026-06-13

---

## Online Sources

### Fitch (1970) — original orthology / paralogy definitions (via PMC obituary that quotes the paper verbatim)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3178060/ (Koonin EV, 2011. "Walter Fitch and the orthology paradigm", quoting Fitch 1970, *Systematic Zoology* 19:99–106)
**Accessed:** 2026-06-13
**Authority rank:** 1 (primary definition, peer-reviewed)
**Retrieved how:** WebSearch "Fitch 1970 orthology paralogy definition…" then WebFetch of the PMC article; quotations below are copied from the fetched text.

**Key Extracted Points:**

1. **Orthology (speciation):** Quoting Fitch 1970: *"Where the homology is the result of speciation so that the history of the gene reflects the history of the species (for example a hemoglobin in man and mouse) the genes should be called orthologous (ortho = exact)."*
2. **Paralogy (gene duplication):** Quoting Fitch 1970: *"Where the homology is the result of gene duplication so that both copies have descended side by side during the history of an organism (for example, α and β hemoglobin) the genes should be called paralogous (para = in parallel)."*
3. **Consequence:** Orthologs split by speciation (between genomes); paralogs split by duplication (within a genome / lineage).

### Tatusov, Koonin & Lipman (1997) — COG, symmetrical best hits

**URL:** https://www.science.org/doi/10.1126/science.278.5338.631 (citation confirmed via https://www.scirp.org/reference/referencespapers?referenceid=2329029); PubMed https://pubmed.ncbi.nlm.nih.gov/9381173/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed, *Science*)
**Retrieved how:** WebSearch "Tatusov 1997 Genomic perspective protein families COG reciprocal best hit"; full-text PDF (fire.biol.wwu.edu) returned 404 and science.org returned 403, so the method statement below is taken from the search-result summary plus the scirp citation page that were fetched; the citation (Science 278:631–637, 1997, DOI 10.1126/science.278.5338.631) was confirmed from the fetched scirp page.

**Key Extracted Points:**

1. **Method:** COGs are built from "reciprocal best BLAST hits" / symmetrical best hits: each COG consists of orthologous proteins (or orthologous sets of paralogs) detected by genome-to-genome best-hit comparison across lineages.
2. **Orthologs + paralogs grouped:** "Each COG consists of individual orthologous proteins or orthologous sets of paralogs from at least three lineages" — i.e. recent paralogs are co-clustered with the ortholog.

### Moreno-Hagelsieb & Latimer (2008) — operational RBH definition and thresholds

**URL:** https://academic.oup.com/bioinformatics/article/24/3/319/252715 — DOI https://doi.org/10.1093/bioinformatics/btm585
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed, *Bioinformatics*)
**Retrieved how:** WebSearch then WebFetch of the OUP article URL; facts below copied from the fetched text.

**Key Extracted Points:**

1. **RBH definition:** *"two genes residing in two different genomes are deemed orthologs if their protein products find each other as the best hit in the opposite genome."*
2. **Best-hit ranking:** hits are sorted "from highest to lowest bit-score, then, if the bit-scores were identical, from smallest to highest E-values"; the first hit is the best hit. → best hit = the candidate with the **maximum similarity score** (ties broken deterministically).
3. **Coverage threshold:** require "coverage of at least 50% of any of the protein sequences in the alignments."
4. **E-value threshold:** maximum E-value 1×10⁻⁶ (a significance gate; in an alignment-free implementation this maps to a minimum-similarity gate).

### Remm, Storm & Sonnhammer (2001) — InParanoid: seed orthologs by RBH, in-paralog rule

**URL:** https://www.sciencedirect.com/science/article/abs/pii/S0022283600951970 (J. Mol. Biol. 314:1041–1052, 2001); review corroboration https://pmc.ncbi.nlm.nih.gov/articles/PMC5674930
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed, *J. Mol. Biol.*)
**Retrieved how:** WebSearch "Remm Storm Sonnhammer 2001 …" (abstract fetched in results) and WebFetch of the Frontiers review PMC5674930 for corroboration; facts below copied from those fetched texts.

**Key Extracted Points:**

1. **Seed:** an ortholog cluster is "seeded by a reciprocally best-matching orthologous pair" (a between-genome RBH).
2. **In-paralog rule (OrthoMCL/InParanoid):** recent paralogs (in-paralogs) are "within-species BLAST hits that are reciprocally better than between-species hits" — i.e. a within-genome pair whose mutual similarity exceeds either gene's similarity to the seed ortholog in the other genome. In-paralogs arose by duplication **after** speciation; outparalogs (before speciation) are excluded.

---

## Documented Corner Cases and Failure Modes

### From Moreno-Hagelsieb & Latimer (2008)

1. **Tie in best hit:** equal bit-scores are broken by smaller E-value; an implementation must apply a deterministic tie-break so the "best hit" is unique.
2. **Coverage filter:** a high-scoring local match over a short region must be rejected unless it covers ≥ 50% of a sequence (avoids spurious orthologs from shared domains).

### From the RBH definition (Tatusov 1997; Moreno-Hagelsieb 2008)

1. **Non-reciprocity:** if gene A's best hit is B but B's best hit is C≠A, then A–B is **not** an ortholog pair (one-directional best hit ≠ ortholog). This is the defect class the unit must guard against.
2. **No qualifying hit:** a gene with no above-threshold hit in the other genome yields no ortholog pair.

### From Fitch (1970) / Remm et al. (2001)

1. **Empty / single-gene genome:** with < 1 gene in either genome there can be no between-genome pair; with < 2 genes in a genome there can be no within-genome paralog pair.
2. **Outparalogs:** within-genome pairs that predate speciation are not in-paralogs; the within-genome best-hit rule alone identifies the closest within-genome relatives (recent paralogs) and is the documented operational proxy.

---

## Test Datasets

### Dataset: Reciprocity worked example (derived from the RBH definition, Moreno-Hagelsieb 2008)

**Source:** Moreno-Hagelsieb & Latimer (2008), RBH definition (best hit each way).

Sequences chosen so similarity ranking is unambiguous by shared 5-mer content.

| Gene (genome) | Sequence | Best hit in other genome |
|---------------|----------|--------------------------|
| a1 (G1) | `ACGTACGTACGTAC` | b1 |
| a2 (G1) | `TTTTGGGGCCCCAA` | b2 |
| b1 (G2) | `ACGTACGTACGTAC` | a1 |
| b2 (G2) | `TTTTGGGGCCCCAA` | a2 |

Expected orthologs (RBH): {a1↔b1, a2↔b2}. a1↔b1 and a2↔b2 are mutual best hits (identical sequences ⇒ maximal similarity each way).

### Dataset: Non-reciprocity (one-directional best hit is NOT an ortholog)

**Source:** Tatusov 1997 symmetrical-best-hit requirement; Moreno-Hagelsieb 2008.

| Gene (genome) | Sequence |
|---------------|----------|
| a1 (G1) | `ACGTACGTACGTAC` |
| b1 (G2) | `ACGTACGTACGTAC` (identical to a1) |
| b2 (G2) | `ACGTACGTACGTACGTACGT` (a1's superstring; shares all of a1's 5-mers) |

a1's best hit is b1 (Jaccard 1.0) over b2 (Jaccard < 1.0). b1's best hit is a1. b2's best hit is a1 too, but a1's best hit is b1, not b2. Expected orthologs: {a1↔b1} only; b2 is excluded (not reciprocal). RBH count = 1.

### Dataset: In-paralog (within-genome best hit)

**Source:** Fitch (1970) paralogy = within-genome duplication; Remm et al. (2001) in-paralog rule.

| Gene (genome) | Sequence |
|---------------|----------|
| p1 (G1) | `GGGGCCCCAAAATT` |
| p2 (G1) | `GGGGCCCCAAAATT` (duplicate of p1, identical) |
| q1 (G1) | `ACGTACGTACGTAC` (unrelated) |

Expected paralogs within G1: {p1↔p2} (mutual within-genome best hits, similarity 1.0); q1 has no qualifying within-genome partner.

---

## Assumptions

1. **ASSUMPTION: Alignment-free similarity score as the best-hit ranking metric.** Moreno-Hagelsieb (2008) ranks best hits by BLAST bit-score. The repository `ComparativeGenomics` class has no alignment-based bit-score available in its project (the `Seqeron.Genomics.Alignment` project is not referenced). The implementation therefore ranks candidates by a k-mer (5-mer) Jaccard similarity, which is a monotone alignment-free similarity (cf. Mash, Ondov et al. 2016). This affects **which** pair wins ties only when sequences are near-identical; the RBH reciprocity rule, coverage gate, and threshold semantics (the correctness-critical parts) are source-backed. The metric is order-preserving for the test datasets above (identical sequences score 1.0). The 50%-coverage gate is mapped from Moreno-Hagelsieb to "shared k-mers ≥ 50% of the smaller k-mer set."

---

## Recommendations for Test Coverage

1. **MUST Test:** Mutual best hits are returned as ortholog pairs (RBH). — Evidence: Moreno-Hagelsieb & Latimer (2008) RBH definition.
2. **MUST Test:** A one-directional best hit (non-reciprocal) is NOT returned. — Evidence: Tatusov (1997) symmetrical best hits; Moreno-Hagelsieb (2008).
3. **MUST Test:** Within-genome mutual best hits are returned as paralog pairs; unrelated gene excluded. — Evidence: Fitch (1970); Remm et al. (2001).
4. **MUST Test:** Coverage / minimum-identity threshold rejects sub-threshold pairs. — Evidence: Moreno-Hagelsieb (2008) ≥ 50% coverage.
5. **MUST Test:** Empty genome → no orthologs; single-gene genome → no paralogs. — Evidence: definitions require a pair.
6. **MUST Test:** Null inputs throw `ArgumentNullException`. — Rationale: repository contract (sibling methods).
7. **SHOULD Test:** Each ortholog pair is symmetric / no gene paired twice in one direction. — Rationale: RBH yields a partial matching.
8. **SHOULD Test:** Genes with empty/absent sequence are skipped. — Rationale: similarity undefined without sequence.
9. **COULD Test:** Determinism — same input yields identical pair set across runs. — Rationale: order-independence requirement.

---

## References

1. Fitch WM. 1970. Distinguishing homologous from analogous proteins. *Systematic Zoology* 19(2):99–113 (definitions p.99–106). Quoted in: Koonin EV. 2011. Walter Fitch and the orthology paradigm. *Brief Bioinform* — https://pmc.ncbi.nlm.nih.gov/articles/PMC3178060/
2. Tatusov RL, Koonin EV, Lipman DJ. 1997. A genomic perspective on protein families. *Science* 278(5338):631–637. https://doi.org/10.1126/science.278.5338.631 — https://pubmed.ncbi.nlm.nih.gov/9381173/
3. Moreno-Hagelsieb G, Latimer K. 2008. Choosing BLAST options for better detection of orthologs as reciprocal best hits. *Bioinformatics* 24(3):319–324. https://doi.org/10.1093/bioinformatics/btm585
4. Remm M, Storm CEV, Sonnhammer ELL. 2001. Automatic clustering of orthologs and in-paralogs from pairwise species comparisons. *J. Mol. Biol.* 314(5):1041–1052. https://doi.org/10.1006/jmbi.2000.5197 — https://www.sciencedirect.com/science/article/abs/pii/S0022283600951970
5. Li L, Stoeckert CJ, Roos DS. 2003. OrthoMCL: identification of ortholog groups for eukaryotic genomes. *Genome Res.* 13(9):2178–2189. https://doi.org/10.1101/gr.1224503 (in-paralog corroboration via review https://pmc.ncbi.nlm.nih.gov/articles/PMC5674930)
6. Ondov BD, et al. 2016. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol.* 17:132. https://doi.org/10.1186/s13059-016-0997-x (basis for alignment-free k-mer similarity ranking; see Assumption 1)

---

## Change History

- **2026-06-13**: Initial documentation.
