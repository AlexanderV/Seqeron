# Evidence Artifact: TRANS-SPLICE-001

**Test Unit ID:** TRANS-SPLICE-001
**Algorithm:** Alternative Splicing — Event Classification and Percent-Spliced-In (PSI / Ψ)
**Date Collected:** 2026-06-13

---

## Online Sources

### Wang et al. (2008) — Alternative isoform regulation in human tissue transcriptomes (Nature)

**URL:** https://pubmed.ncbi.nlm.nih.gov/18978772/ (record); https://www.nature.com/articles/nature07509 (article)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary paper)
**Retrieved via:** WebSearch query *"Wang 2008 Nature Alternative isoform regulation in human tissue transcriptomes PMC full text exon skipping mutually exclusive frequency"*, then WebFetch of the Nature article URL (returned a redirect to an IdP authorize URL; the bibliographic facts and event-type list were taken from the WebSearch result block citing the article).

**Key Extracted Points:**

1. **Five canonical AS event classes:** "Five types of alternative splicing … which are differential exon inclusion (exon skipping), intron retention, alternative 5′ and 3′ exon splicing, and mutually exclusive splicing (Wang et al, 2008)." (extracted from WebSearch result summarizing/citing Wang 2008).
2. **PSI as inclusion fraction:** the fraction of transcripts that include a cassette exon is estimated from the relative density of inclusion-supporting reads vs exclusion-supporting reads; alternative splicing is pervasive (92–94% of multi-exon human genes).
3. **Citation:** Wang ET, Sandberg R, Luo S, Khrebtukova I, Zhang L, Mayr C, Kingsmore SF, Schroth GP, Burge CB. Nature. 2008;456(7221):470–476.

---

### Pruitt context — Challenges in estimating percent inclusion of alternatively spliced junctions (BMC Bioinformatics, PMC3330053)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3330053/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary paper)
**Retrieved via:** WebSearch query *"Percent Spliced In PSI formula inclusion exclusion reads Katz MISO alternative splicing quantification"*, then WebFetch of the PMC article (followed the 301 redirect to the `pmc.ncbi.nlm.nih.gov` host).

**Key Extracted Points:**

1. **Core PSI definition (verbatim, paraphrased from fetched text):** "PSI is defined as the expression of isoforms containing the alternatively spliced exon as a fraction of the total expression for both alternatively and constitutively spliced isoforms." → Ψ = inclusion / (inclusion + exclusion).
2. **Gaussian-model form (verbatim from fetched text):** `μ̃ = γᵢ / (γᵢ + γₑ)`, where γᵢ is normalized expression of the inclusion junction and γₑ of the exclusion junction.
3. **Read classification:** reads aligning to the alternative exon or to its junctions with adjacent constitutive exons support inclusion; reads aligning to the junction between the two adjacent constitutive exons support exclusion.

---

### Shen et al. (2014) — rMATS (PNAS, PMC4280593)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4280593/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary paper) / 3 (rMATS is the reference implementation)
**Retrieved via:** WebSearch query *"rMATS Shen 2014 PNAS replicate multivariate analysis transcript splicing PSI inclusion isoform length normalization formula PMC"*, then WebFetch of the PMC article.

**Key Extracted Points:**

1. **Length-normalized PSI (verbatim from fetched text):** `ψ̂ = (I/lᵢ) / (I/lᵢ + S/lₛ)`, where I = inclusion read count, S = skipping read count, lᵢ = effective length of the inclusion isoform, lₛ = effective length of the skipping isoform.
2. **Binomial model (verbatim from fetched text):** `I | ψ ∼ Binomial(n = I+S, p = f(ψ) = lᵢψ / (lᵢψ + lₛ(1−ψ)))`; f normalizes ψ by the effective isoform lengths.
3. **Effective length definition (verbatim):** "the number of unique isoform-specific read positions."
4. **PSI meaning:** "the ratio of isoforms (supported by junction-spanning reads) that include or exclude a specific RNA segment."
5. **Citation:** Shen S, Park JW, Lu ZX, Lin L, Henry MD, Wu YN, Zhou Q, Xing Y. PNAS. 2014;111(51):E5593–E5601.

---

### rMATS project page (reference implementation) + SUPPA2 (Trincado et al. 2018)

**URL:** https://rmats.sourceforge.io/ ; SUPPA via https://pubmed.ncbi.nlm.nih.gov/29571299/ and WebSearch summary
**Accessed:** 2026-06-13
**Authority rank:** 3 (established reference implementations)
**Retrieved via:** WebFetch of `https://rmats.sourceforge.io/`; WebSearch query *"SUPPA2 PSI percent spliced in formula skipped exon inclusion isoform exclusion junction definition transcript abundance"*.

**Key Extracted Points:**

1. **Event-type vocabulary (rMATS):** rMATS detects all major AS patterns; output files name event classes SE (skipped exon), A5SS, A3SS, MXE, RI — matching the Wang 2008 five-class taxonomy.
1a. **A5SS / A3SS coordinate convention (verified 2026-06-15 vs https://github.com/Xinglab/rmats-turbo README, columns `longExonStart_0base,longExonEnd,shortES,shortEE,flankingES,flankingEE`):** on the **+ strand**, for **A5SS** the long/short exon forms differ at their **downstream (END/3′) boundary** while sharing the upstream (START/5′) boundary (alternative **donor** = 5′ splice site = 3′ end of the upstream exon); for **A3SS** they differ at their **upstream (START/5′) boundary** while sharing the downstream (END/3′) boundary (alternative **acceptor** = 3′ splice site = 5′ start of the downstream exon). Confirmed independently by molecular-biology splice-site definitions (5′ splice site = donor at the exon→intron boundary; 3′ splice site = acceptor at the intron→exon boundary) and by NAR 34(21):6305 (alternative 5′SS = exon extension/truncation at the downstream boundary of the upstream exon).
2. **SUPPA PSI (verbatim from fetched WebSearch summary):** "PSI is calculated as the number of reads covering the included exon divided by the sum of reads covering the skipped exon and reads covering the included exon" — i.e. Ψ = inclusion/(inclusion+exclusion); inclusion = junctions between the upstream exon and the exon of interest, exclusion = junction between upstream and downstream exons skipping the exon of interest.

---

## Documented Corner Cases and Failure Modes

### From PMC3330053 (BMC Bioinformatics) and Shen 2014 (rMATS)

1. **Zero total reads (0/0):** PSI is undefined when inclusion + exclusion = 0. The Gaussian model in PMC3330053 adds a pseudo-count (`1/|P|`) to the normalized junction expression to avoid division by zero; in practice an event with no supporting reads yields no PSI estimate.
2. **No exclusion reads:** when only inclusion reads exist, PSI = 1 (fully included). When only skipping reads exist, PSI = 0 (fully excluded). Direct consequence of Ψ = I/(I+S).
3. **Length normalization changes the estimate:** rMATS divides each read count by the isoform's effective length before forming the ratio, so the unnormalized Ψ = I/(I+S) and the normalized Ψ = (I/lᵢ)/(I/lᵢ + S/lₛ) differ whenever lᵢ ≠ lₛ.

### From Wang 2008

1. **Event classification requires ≥2 isoforms:** an alternative-splicing event is defined relative to two transcript isoforms of the same gene that differ in exon structure; a single isoform cannot define an event.

---

## Test Datasets

### Dataset: Skipped-exon PSI (read-count ratio) — derived from Ψ = I/(I+S)

**Source:** PMC3330053 core definition; SUPPA PSI definition.

| Parameter | Value |
|-----------|-------|
| Inclusion reads I | 80 |
| Skipping reads S | 20 |
| Expected PSI (unnormalized) | 80 / (80+20) = 0.80 |

### Dataset: Length-normalized PSI — derived from rMATS ψ̂ = (I/lᵢ)/(I/lᵢ + S/lₛ)

**Source:** Shen et al. (2014) rMATS, PMC4280593.

| Parameter | Value |
|-----------|-------|
| Inclusion reads I | 80 |
| Skipping reads S | 20 |
| Inclusion effective length lᵢ | 200 |
| Skipping effective length lₛ | 100 |
| I/lᵢ | 0.40 |
| S/lₛ | 0.20 |
| Expected PSI (normalized) | 0.40 / (0.40+0.20) = 0.666666… |

### Dataset: Event classification (two isoforms of one gene)

**Source:** Wang et al. (2008) five-class taxonomy.

| Isoform A exons | Isoform B exons | Differs by | Expected class |
|-----------------|-----------------|-----------|----------------|
| (1,100),(200,300),(400,500) | (1,100),(400,500) | middle exon (200,300) present in A, absent in B | SkippedExon |
| (1,100),(200,300) | (1,100),(101,199),(200,300) | intron (101,199) retained as exon-body in B | RetainedIntron |
| (1,100),(200,300) | (1,100),(200,350) | shared 5′ start 200, different 3′ end (300 vs 350) ⇒ alternative **donor** (5′ splice site = exon END) | AlternativeFivePrimeSS (A5SS) |
| (1,100),(200,300) | (1,100),(150,300) | shared 3′ end 300, different 5′ start (200 vs 150) ⇒ alternative **acceptor** (3′ splice site = exon START) | AlternativeThreePrimeSS (A3SS) |
| (1,100),(200,300),(500,600) | (1,100),(350,400),(500,600) | exactly one of two alternative middle exons used | MutuallyExclusiveExons |

---

## Assumptions

1. **ASSUMPTION: Length normalization is opt-in.** The two cited definitions (unnormalized Ψ = I/(I+S) and rMATS length-normalized Ψ) are both authoritative. `CalculatePSI` implements the unnormalized read-count ratio by default (the definition shared by Wang 2008, PMC3330053, and SUPPA) and exposes effective lengths as optional parameters that, when both > 0, switch to the rMATS form. This is an API-shape choice; both numerical behaviors are source-backed.
2. **ASSUMPTION: Strand is forward / coordinates are ascending.** The event classifier compares exon coordinate lists assuming exons are ordered 5′→3′ on the same strand. Wang 2008 defines events on a per-gene basis; strand handling is an input-normalization concern outside the formula.

---

## Recommendations for Test Coverage

1. **MUST Test:** `CalculatePSI` returns I/(I+S) for the (80,20) dataset → 0.80. — Evidence: PMC3330053; SUPPA.
2. **MUST Test:** `CalculatePSI` with effective lengths returns the rMATS normalized value 0.6666… for (I=80,S=20,lᵢ=200,lₛ=100). — Evidence: Shen 2014 rMATS (PMC4280593).
3. **MUST Test:** PSI = 1 when S = 0; PSI = 0 when I = 0. — Evidence: Ψ = I/(I+S) (PMC3330053).
4. **MUST Test:** PSI of 0/0 (no reads) is the documented undefined case → NaN. — Evidence: 0/0 undefined (PMC3330053 pseudo-count rationale).
5. **MUST Test:** `DetectAlternativeSplicing` classifies a skipped-exon isoform pair as SkippedExon. — Evidence: Wang 2008.
6. **MUST Test:** classification of RetainedIntron, AlternativeThreePrimeSS, AlternativeFivePrimeSS, MutuallyExclusiveExons. — Evidence: Wang 2008.
7. **SHOULD Test:** invariant 0 ≤ PSI ≤ 1 for any non-negative inputs. — Rationale: ratio of a part to a whole.
8. **SHOULD Test:** fewer than two isoforms for a gene yields no event. — Rationale: an event needs two isoforms (Wang 2008).
9. **COULD Test:** identical isoforms produce no event. — Rationale: no structural difference ⇒ no AS.

---

## References

1. Wang ET, Sandberg R, Luo S, Khrebtukova I, Zhang L, Mayr C, Kingsmore SF, Schroth GP, Burge CB. (2008). Alternative isoform regulation in human tissue transcriptomes. Nature 456(7221):470–476. https://doi.org/10.1038/nature07509 (PubMed: https://pubmed.ncbi.nlm.nih.gov/18978772/)
2. (BMC Bioinformatics, 2012). Challenges in estimating percent inclusion of alternatively spliced junctions from RNA-seq data. BMC Bioinformatics 13(Suppl 6):S11. https://pmc.ncbi.nlm.nih.gov/articles/PMC3330053/
3. Shen S, Park JW, Lu ZX, Lin L, Henry MD, Wu YN, Zhou Q, Xing Y. (2014). rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. PNAS 111(51):E5593–E5601. https://doi.org/10.1073/pnas.1419161111 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC4280593/)
4. rMATS project documentation. https://rmats.sourceforge.io/ (accessed 2026-06-13)
5. Trincado JL, Entizne JC, Hysenaj G, Singh B, Skalic M, Elliott DJ, Eyras E. (2018). SUPPA2: fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions. (PubMed: https://pubmed.ncbi.nlm.nih.gov/29571299/)

---

## Change History

- **2026-06-13**: Initial documentation.
