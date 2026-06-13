# Evidence Artifact: SV-DETECT-001

**Test Unit ID:** SV-DETECT-001
**Algorithm:** Structural Variant Detection from Paired-End Mapping (PEM) signatures
**Date Collected:** 2026-06-13

---

## Online Sources

### Medvedev P, Stanciu M, Brudno M (2009) — "Computational methods for discovering structural variation with next-generation sequencing", Nature Methods 6(11s):S13–S20

**URL:** http://www.cs.toronto.edu/~brudno/nmeth_review09.pdf (author copy of DOI:10.1038/nmeth.1374)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed review)
**Retrieved how:** WebSearch query `computational methods discovering structural variation paired-end read signatures deletion insertion inversion span orientation Medvedev Stanciu Brudno 2009`; WebFetch of the PDF; PDF text extracted with `pypdf` from the saved binary and read verbatim.

**Key Extracted Points:**

1. **Basic deletion signature:** "A mate pair that spans an isolated deletion event maps to the corresponding regions of the reference, but the mapped distance is greater than the insert size." → deletion ⇒ observed span LARGER than expected insert size.
2. **Basic insertion signature:** "Conversely, if the event is an insertion, then the distance is smaller." → insertion ⇒ observed span SMALLER than expected insert size.
3. **Basic inversion signature:** "A mate pair that spans either (but not both) of its breakpoints will map to the reference with the orientation of the read, lying within the inversion, flipped. Two such mate pairs ... form the 'basic inversion' signature." → inversion ⇒ one mate's orientation flipped relative to the concordant pair.
4. **Linking signature / translocation:** "A mate pair spanning the donor's breakpoint will map with a distance much greater than the insert size. ... Other linking signatures can connect regions that are arbitrarily distant or even on different chromosomes." → mates on different chromosomes form a linking (translocation) signature.
5. **Signature-then-cluster paradigm:** methods "distinguish themselves by the signatures they can detect, and ... the way they cluster or window these signatures" — discordant pairs supporting the same event are clustered.

### Chen K et al. (2009) — BreakDancer, distribution README (BreakDancer-1.3.6, genome/breakdancer)

**URL:** https://raw.githubusercontent.com/genome/breakdancer/master/README
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation documentation; primary paper DOI:10.1038/nmeth.1363, Nat Methods 6:677–681)
**Retrieved how:** WebSearch `BreakDancer discordant read pairs structural variant detection insert size standard deviation Chen 2009`; `curl` confirmed the file exists; WebFetch of the raw README.

**Key Extracted Points:**

1. **Anomaly cutoff (-c):** "cutoff in unit of standard deviation"; default 3 standard deviations. "upper threshold = mean + (std × -c value) and lower threshold = mean − (std × -c value)." → a pair is discordant by span if insertSize < mean − c·sd OR insertSize > mean + c·sd.
2. **Minimum support (-r):** "minimum number of read pairs required to establish a connection"; default 2. → an SV call requires ≥ 2 supporting read pairs by default.
3. **Detection inputs:** anomalous read pairs are identified by "unexpected separation distances or orientation."
4. **SV type codes:** DEL = deletions, INS = insertions, INV = inversions, ITX = intra-chromosomal translocations, CTX = inter-chromosomal translocations.

### BreakDancer protocol (Fan X, Chen K et al.) — Curr Protoc Bioinformatics, PMC3661775

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3661775/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed protocol)
**Retrieved how:** WebFetch (followed 301 redirect from ncbi.nlm.nih.gov to pmc.ncbi.nlm.nih.gov).

**Key Extracted Points:**

1. **Classification basis:** read pairs are "independently classified into six types ... based on 1) the separation distance and alignment orientation between the paired reads, 2) the user-specified threshold, and 3) the empirical insert size distribution." → the six classes are normal + DEL/INS/INV/ITX/CTX.
2. **Threshold:** experiments use separation thresholds of "3 s.d." or "4 s.d." from the mean insert size, confirming the standard-deviation cutoff convention.

### Concordant Illumina paired-end orientation (FR) — cureffi.org, citing Heng Li / BWA & SAM proper-pair flag

**URL:** https://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/
**Accessed:** 2026-06-13
**Authority rank:** 4 (secondary article quoting the BWA author and the SAM FLAG 0x02 convention; primary = SAM/BAM spec proper-pair bit and BWA behavior)
**Retrieved how:** WebSearch `paired-end read concordant orientation forward reverse FR proper pair Illumina mate +/- strand structural variant`; WebFetch of the article.

**Key Extracted Points:**

1. **Concordant FR orientation:** "one read should align to the forward strand, and the other should align to the reverse strand, ... so that they are pointed towards one another." → concordant pair has one mate on + and the mate on −.
2. **Abnormal orientation:** "If they instead align RF, FF or RR, that's a problem and often indicates the reads aligned incorrectly ... or ... a real inversion or translocation exists." → same-orientation (FF/RR) supports inversion.
3. **Proper-pair flag:** "BWA will only set the 'proper pair flag' to 1 for Illumina reads aligned FR." (SAM FLAG 0x02). → FR is the formal concordant orientation.

---

## Documented Corner Cases and Failure Modes

### From Medvedev, Stanciu & Brudno (2009)

1. **Insertion larger than insert size:** "the basic insertion signature does not appear when the size of the insertion is greater than the insert size of the sequenced fragment, and it does not indicate the inserted sequence itself." → insertions above the fragment length are invisible to PEM span; inserted sequence is not recovered from the span signature.

### From BreakDancer README

1. **Below-support clusters:** clusters with fewer than the minimum supporting read pairs (default 2) are not reported as SVs.

---

## Test Datasets

### Dataset: Synthetic PEM signatures (derived from the cited signature rules)

**Source:** Medvedev, Stanciu & Brudno (2009) signature definitions + BreakDancer cutoff (mean=400, sd=50, c=3 ⇒ bounds [250, 550]).

| Parameter | Value |
|-----------|-------|
| Expected (mean) insert size | 400 |
| Insert size s.d. | 50 |
| Cutoff c (standard deviations) | 3 |
| Lower span bound (mean − 3·sd) | 250 |
| Upper span bound (mean + 3·sd) | 550 |
| Concordant pair (chr1==chr2, FR, span 400) | NOT discordant |
| Large-span pair (span 5000, FR, same chr) | Deletion signature (span > 550) |
| Small-span pair (span 100, FR, same chr) | Insertion signature (span < 250) |
| Same-orientation pair (FF, same chr, span 400) | Inversion signature |
| Inter-chromosomal pair (chr1 ≠ chr2) | Translocation signature |

---

## Assumptions

1. **ASSUMPTION: Translocation vs. Inversion precedence for inter-chromosomal same-orientation pairs.** When mates map to different chromosomes, the inter-chromosomal (translocation/CTX) signature is reported regardless of relative orientation. The cited sources define inter-chromosomal mapping as a linking/translocation signature (Medvedev et al. 2009; BreakDancer CTX) and define inversion (INV) only for intra-chromosomal flipped pairs; chromosome difference is therefore evaluated first. Justification: a flipped orientation is undefined as "inversion" across chromosomes, where the event is by definition a translocation.

---

## Recommendations for Test Coverage

1. **MUST Test:** Span greater than (mean + c·sd) on the same chromosome with concordant FR orientation classifies as Deletion. — Evidence: Medvedev et al. 2009 ("mapped distance is greater than the insert size"); BreakDancer DEL.
2. **MUST Test:** Span smaller than (mean − c·sd) on the same chromosome classifies as Insertion. — Evidence: Medvedev et al. 2009 ("the distance is smaller"); BreakDancer INS.
3. **MUST Test:** Same-orientation (FF or RR) intra-chromosomal pair classifies as Inversion. — Evidence: Medvedev et al. 2009 (flipped orientation); cureffi/BWA (FF/RR ⇒ inversion).
4. **MUST Test:** Different-chromosome pair classifies as Translocation. — Evidence: Medvedev et al. 2009 (linking signature across chromosomes); BreakDancer CTX.
5. **MUST Test:** Concordant pair (same chr, FR, span within [mean−c·sd, mean+c·sd]) is NOT discordant and yields no SV. — Evidence: BreakDancer "normal" class; FR proper pair.
6. **MUST Test:** A cluster with fewer than min-support pairs yields no SV; ≥ min-support yields one SV. — Evidence: BreakDancer -r default 2.
7. **SHOULD Test:** Cutoff bounds are exact: span exactly at mean±c·sd is concordant; one unit beyond is discordant. — Rationale: boundary correctness of the s.d. cutoff.
8. **COULD Test:** Empty input yields empty output. — Rationale: defined trivial behavior.

---

## References

1. Medvedev P, Stanciu M, Brudno M. 2009. Computational methods for discovering structural variation with next-generation sequencing. Nature Methods 6(11s):S13–S20. https://doi.org/10.1038/nmeth.1374 (author PDF: http://www.cs.toronto.edu/~brudno/nmeth_review09.pdf)
2. Chen K, Wallis JW, McLellan MD, et al. 2009. BreakDancer: an algorithm for high-resolution mapping of genomic structural variation. Nature Methods 6:677–681. https://doi.org/10.1038/nmeth.1363 (distribution README: https://raw.githubusercontent.com/genome/breakdancer/master/README)
3. Fan X, Abbott TE, Larson D, Chen K. 2014. BreakDancer: Identification of Genomic Structural Variation from Paired-End Read Mapping. Curr Protoc Bioinformatics 45:15.6.1–15.6.11. https://pmc.ncbi.nlm.nih.gov/articles/PMC3661775/
4. Kennedy A. 2012. Forward and reverse reads in paired-end sequencing (citing the SAM proper-pair FLAG 0x02 and BWA FR convention). https://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/

---

## Change History

- **2026-06-13**: Initial documentation.
