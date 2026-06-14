# Evidence Artifact: COMPGEN-COMPARE-001

**Test Unit ID:** COMPGEN-COMPARE-001
**Algorithm:** Comprehensive two-genome comparison — core/dispensable (conserved vs genome-specific) gene partition and the syntenic-gene fraction
**Date Collected:** 2026-06-14

---

## Online Sources

### Tettelin et al. (2005) — pan-genome: core genome vs dispensable (genome-specific) genome

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC1216834/ — DOI https://doi.org/10.1073/pnas.0506758102 (PNAS 102(39):13950–13955)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed, *PNAS*; the paper that coined "pan-genome", "core genome", "dispensable genome")
**Retrieved how:** WebSearch "Tettelin 2005 PNAS genome analysis Streptococcus agalactiae pan-genome core dispensable doi 10.1073", then WebFetch of the PMC full-text PMC1216834. The quotations below are taken from the fetched article text.

**Key Extracted Points:**

1. **Core genome:** the species pan-genome comprises *"a core genome containing genes present in all strains"* — i.e. a gene shared by every genome in the compared set is a **core / conserved** gene.
2. **Dispensable genome:** *"a dispensable genome composed of genes absent from one or more strains and genes that are unique to each strain."* → a gene missing from one or more genomes (in the two-genome case: present in only one) is **dispensable / genome-specific**.
3. **Strain-specific (unique) genes:** *"genes that are unique to each strain"* form part of the dispensable genome; in the pairwise case these are exactly the genes of one genome with no ortholog in the other.
4. **Conservation gate (when is a gene "present"/conserved):** *"A gene was considered conserved if at least one of these three methods produced an alignment with a minimum of 50% sequence conservation over 50% of the protein/gene length."* → presence/conservation is decided by an alignment gate of ≥50% identity over ≥50% length, not by exact identity.

### Reciprocal Best Hits — operational ortholog (shared-gene) criterion (Moreno-Hagelsieb & Latimer 2008; Tatusov et al. 1997)

**URL:** https://academic.oup.com/bioinformatics/article/24/3/319/252715 — DOI https://doi.org/10.1093/bioinformatics/btm585
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed, *Bioinformatics*)
**Retrieved how:** Reused from this session's lineage of comparative-genomics units (COMPGEN-RBH-001 Evidence, retrieved 2026-06-14 by WebSearch + WebFetch of the OUP article). Quotation re-confirmed against that artifact.

**Key Extracted Points:**

1. **Shared gene = reciprocal best hit:** *"two genes residing in two different genomes are deemed orthologs if their protein products find each other as the best hit in the opposite genome."* → in `CompareGenomes` the **conserved (core) genes** are exactly the reciprocal-best-hit ortholog pairs; everything else in each genome is genome-specific. This is the validated sub-unit COMPGEN-RBH-001.
2. **Coverage gate:** *"we also required coverage of at least 50% of any of the protein sequences in the alignments."* → operationalises Tettelin's "over 50% of the … length".

### Synteny — fraction of syntenic (collinear) genes as a conservation metric

**URL:** https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/synteny (ScienceDirect "Synteny — an overview", citing primary studies) and https://en.wikipedia.org/wiki/Synteny
**Accessed:** 2026-06-14
**Authority rank:** 4 (overview citing primaries) / 4 (Wikipedia citing primaries)
**Retrieved how:** WebSearch "synteny conservation index fraction of genes in syntenic blocks proportion conserved gene order metric"; the returned ScienceDirect/Wikipedia summaries report the metric below.

**Key Extracted Points:**

1. **Syntenic blocks:** *"conservation of blocks of order within two sets of chromosomes … referred to as syntenic blocks"* (Wikipedia, Synteny). A syntenic block is a run of collinear orthologs.
2. **Fraction of syntenic genes:** *"The fraction of syntenic genes is a metric used to measure synteny conservation"* (e.g. 50.2% for the coral genus *Acropora*). → an overall-synteny score is naturally the fraction of genes that fall in syntenic blocks; `CompareGenomes` reports `OverallSynteny = (genes in syntenic blocks) / min(|genome1|, |genome2|)`, clamped to ≤1.0.
3. **Block construction (MCScanX):** the syntenic blocks themselves come from the MCScanX collinearity model (validated sub-unit COMPGEN-SYNTENY-001), which reports non-overlapping chains scoring ≥250 (≥5 collinear anchors at MatchScore 50). Wang et al. 2012, *NAR* 40(7):e49, https://doi.org/10.1093/nar/gkr1293.

---

## Documented Corner Cases and Failure Modes

### From Tettelin et al. (2005)

1. **All genes shared → all core:** if every gene of each genome has an ortholog in the other, the dispensable/genome-specific count is 0 for both genomes (all genes are core).
2. **No genes shared → all dispensable:** if no ortholog pairs exist, the core (conserved) count is 0 and every gene of each genome is genome-specific.
3. **One shared, one unique each:** a gene present in both genomes is core (conserved); a gene present in only one genome is that genome's dispensable/specific gene.

### From the RBH criterion (Moreno-Hagelsieb 2008)

1. **One-directional best hit is not a shared gene:** a gene whose best hit is not reciprocated does not count toward the conserved/core set, hence remains genome-specific.
2. **Empty genome:** with no sequence-bearing genes in a genome there are no ortholog pairs; conserved = 0.

### From synteny (MCScanX collinearity)

1. **Too few collinear anchors:** fewer than 5 collinear orthologs form no reported syntenic block, so `OverallSynteny` can be 0 even when conserved orthologs exist (a documented behavior of the block threshold; see Assumption 2).

---

## Test Datasets

### Dataset: Two-genome core/dispensable partition (Tettelin 2005 model)

**Source:** Tettelin et al. (2005) core/dispensable definitions; ortholog detection by RBH (Moreno-Hagelsieb 2008). Sequences chosen so the 5-mer-content RBH similarity ranking is unambiguous (≥0.3 Jaccard, ≥0.5 coverage for the intended pairs; well below for the unique genes).

| Scenario | genome1 genes | genome2 genes | shared (orthologous) pairs | Expected Conserved (core) | Expected Specific1 | Expected Specific2 |
|----------|---------------|---------------|----------------------------|---------------------------|--------------------|--------------------|
| One shared, one unique each | a1(=S), b1(=U1) | c2(=S), d2(=U2) | a1↔c2 | 1 | 1 | 1 |
| Disjoint content | x1(=U1), y1(=W1) | x2(=U2), y2(=W2) | (none) | 0 | 2 | 2 |
| Identical content (5 collinear + 1 unique each) | 5 shared S₀…S₄, u1 | 5 shared S₀…S₄, u2 | 5 collinear | 5 | 1 | 1 |

Here S, S₀…S₄ are distinct ≥60-nt sequences each shared identically by both genomes; U1, U2, W1, W2 are mutually dissimilar 60-nt sequences (each genome-specific).

### Dataset: Syntenic-gene fraction (identical-content scenario)

**Source:** "fraction of syntenic genes" metric (ScienceDirect Synteny overview); block from MCScanX (Wang et al. 2012).

| Parameter | Value |
|-----------|-------|
| genome1 / genome2 sizes | 6 / 6 |
| collinear orthologs (one block) | 5 |
| `OverallSynteny` | 5 / min(6,6) = 5/6 = 0.8333… |
| Rearrangements (identity permutation) | 0 breakpoints |

---

## Assumptions

1. **ASSUMPTION: Alignment-free similarity in place of the Tettelin 50%/50% alignment gate.** The conserved-gene set is found by RBH, and RBH similarity here is 5-mer-content Jaccard (identity ≥0.3) with k-mer coverage ≥0.5, not a Needleman–Wunsch/BLAST alignment. This maps the Tettelin "≥50% conservation over ≥50% length" and the Moreno-Hagelsieb "≥50% coverage, E≤1e-6" gates onto alignment-free space (inherited verbatim from the validated sub-unit COMPGEN-RBH-001, Assumption 1). It does not change the *partition logic* tested here (core = reciprocal pairs, specific = the rest); identical sequences pass the gate and disjoint sequences fail it, which is all the partition tests rely on.

2. **ASSUMPTION: Minimum syntenic block size = 5 collinear anchors for `OverallSynteny`.** `CompareGenomes` reports `OverallSynteny` as the fraction of genes inside MCScanX syntenic blocks, and MCScanX only reports chains scoring ≥250 (≥5 collinear anchors). Hence `OverallSynteny` can be 0 even with a few conserved orthologs. This threshold is the MCScanX default (Wang et al. 2012; validated in COMPGEN-SYNTENY-001), not invented here.

---

## Recommendations for Test Coverage

1. **MUST Test:** one-shared-one-unique partition → Conserved=1, Specific1=1, Specific2=1 — Evidence: Tettelin (2005) core/dispensable definition.
2. **MUST Test:** disjoint content → Conserved=0, Specific1=|g1|, Specific2=|g2| (all dispensable) — Evidence: Tettelin (2005) "unique to each strain".
3. **MUST Test:** identical content with ≥5 collinear orthologs → Conserved=all-shared, Specific=non-shared, `OverallSynteny` = sharedInBlock/min(|g1|,|g2|), no rearrangements — Evidence: Tettelin (core) + fraction-of-syntenic-genes metric (synteny) + MCScanX block threshold.
4. **MUST Test:** empty genomes → Conserved=0, Specific1=0, Specific2=0, OverallSynteny=0 — Evidence: corner case (no ortholog pairs possible).
5. **SHOULD Test:** orthologs detected but fewer than 5 collinear → OverallSynteny=0 while Conserved>0 — Rationale: documents Assumption 2 boundary.
6. **SHOULD Test:** conserved partition is symmetric (swapping genome1/genome2 swaps Specific1/Specific2, same Conserved) — Rationale: invariant from the symmetric RBH matching.
7. **COULD Test:** `OverallSynteny` is clamped to ≤1.0 — Rationale: documented clamp; degenerate inputs.

---

## References

1. Tettelin H, Masignani V, Cieslewicz MJ, et al. (2005). Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*: Implications for the microbial "pan-genome". *PNAS* 102(39):13950–13955. https://doi.org/10.1073/pnas.0506758102 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC1216834/)
2. Moreno-Hagelsieb G, Latimer K (2008). Choosing BLAST options for better detection of orthologs as reciprocal best hits. *Bioinformatics* 24(3):319–324. https://doi.org/10.1093/bioinformatics/btm585
3. Tatusov RL, Koonin EV, Lipman DJ (1997). A genomic perspective on protein families. *Science* 278(5338):631–637. https://doi.org/10.1126/science.278.5338.631
4. Wang Y, Tang H, DeBarry JD, et al. (2012). MCScanX: a toolkit for detection and evolutionary analysis of gene synteny and collinearity. *Nucleic Acids Research* 40(7):e49. https://doi.org/10.1093/nar/gkr1293
5. Synteny — an overview. ScienceDirect Topics. https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/synteny (accessed 2026-06-14); Wikipedia, "Synteny", https://en.wikipedia.org/wiki/Synteny (accessed 2026-06-14; used for the cited primary metric only).

---

## Change History

- **2026-06-14**: Initial documentation.
