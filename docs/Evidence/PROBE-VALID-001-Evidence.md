# Evidence Artifact: PROBE-VALID-001

**Test Unit ID:** PROBE-VALID-001
**Algorithm:** Hybridization Probe Validation — Gapped (Smith–Waterman) Off-Target Scan
**Date Collected:** 2026-06-24

---

## Online Sources

### Smith & Waterman (1981) — local alignment recurrence

**URL:** https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm (citing Smith TF, Waterman MS (1981) *J Mol Biol* 147(1):195–197)
**Accessed:** 2026-06-24
**Authority rank:** 1 (primary peer-reviewed paper, via Wikipedia citing the primary)

**Key Extracted Points:**

1. **Local-alignment recurrence:** `H(i,j) = max{ H(i-1,j-1) + s(a_i,b_j),  max_{k≥1} H(i-k,j) - W_k,  max_{l≥1} H(i,j-l) - W_l,  0 }`. The fourth alternative — the **zero floor** — is the defining feature: "negative scoring matrix cells are set to zero."
2. **What makes it LOCAL (vs Needleman–Wunsch global):** "any element has a score lower than zero ... will then be set to zero to eliminate influence from previous alignment." Traceback "begins with the highest score, ends when 0 is encountered" — so it returns the best-scoring *subsequence* match rather than an end-to-end alignment.
3. **Gaps / indels:** the `max` over gap lengths `k`/`l` lets the alignment introduce insertions and deletions, which an ungapped fixed-window comparison cannot.

### Altschul et al. (1990) — BLAST (Basic Local Alignment Search Tool)

**URL:** https://en.wikipedia.org/wiki/BLAST_(biotechnology) (citing Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990) *J Mol Biol* 215(3):403–410)
**Accessed:** 2026-06-24
**Authority rank:** 1 (primary peer-reviewed paper, via Wikipedia citing the primary)

**Key Extracted Points:**

1. **Purpose:** a sequence-similarity search tool — compare a query against a library/database to find matches above thresholds. This is the canonical "off-target / homology search" rationale.
2. **Gapped vs ungapped:** the original BLAST produced only ungapped alignments; gapped BLAST (BLAST2) "produce[s] a single alignment with gaps ... allowing detection of insertions and deletions that the original version missed." This is exactly the improvement of a gapped local-alignment off-target scan over a pure ungapped Hamming scan.
3. **Seed-and-extend:** BLAST is a *heuristic* — it locates short seed matches, then extends them into high-scoring segment pairs (HSPs). A full Smith–Waterman scan is exact (not seeded); the "BLAST-grade" property realized here is *gapped local alignment* (indel-aware), not the seed heuristic or a genome-scale index.
4. **Identity thresholds:** the BLAST page does not give a percent-identity cutoff; significance is by E-value/HSP score. Hence the off-target identity threshold is sourced separately (Kane et al. 2000) and exposed as a caller parameter.

### Altschul et al. (1997) — Gapped BLAST and PSI-BLAST

**URL:** https://academic.oup.com/nar/article/25/17/3389/1061651 (Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ (1997) *Nucleic Acids Research* 25(17):3389–3402)
**Accessed:** 2026-06-24
**Authority rank:** 1

**Key Extracted Points:**

1. **Gapped local alignment as the second generation:** the 1997 paper formalizes producing a single gapped alignment, confirming that gap (indel) handling is the recognized improvement over the 1990 ungapped HSPs.

### Kane et al. (2000) — 50-mer oligonucleotide microarray specificity

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC113865/ (Kane MD, Jatkoe TA, Stumpf CR, Lu J, Thomas JD, Madore SJ (2000) *Nucleic Acids Research* 28(22):4552–4557)
**Accessed:** 2026-06-24
**Authority rank:** 1 (peer-reviewed paper) / 5 (curated specificity guideline)

**Key Extracted Points:**

1. **Off-target identity threshold:** "for a given oligonucleotide probe any 'non-target' transcripts (cDNAs) **>75% similar** over the 50 base target may show cross-hybridization." → 0.75 identity over the probe length is the empirically-grounded default above which a hit is called an off-target.
2. **Gene-specific rule:** "oligonucleotide probes with **<75% overall sequence similarity** with non-target sequences and <14 contiguous complementary base pairs are gene-specific" — the complement of the off-target call.
3. **Contiguous-stretch caveat:** "if the 50 base target region is marginally similar, it must not include a stretch of complementary sequence >15 contiguous bases."

---

## Documented Corner Cases and Failure Modes

### From Smith & Waterman (1981)

1. **Local trimming:** the zero floor and "stop traceback at 0" mean a trailing mismatched tail is excluded from the reported local alignment — identity is over the matched core, so a hit with a mismatched end is reported with the core only (e.g. a single trailing mismatch reduces both identity and coverage by the trimmed length).

### From Altschul et al. (1990/1997)

1. **Ungapped scans miss indels:** an off-target reachable only through an insertion or deletion is invisible to a fixed-length ungapped comparison; a single indel shifts the reading frame so every downstream position mismatches.

### From Kane et al. (2000)

1. **Threshold sensitivity:** lowering the identity cutoff admits more (lower-similarity) off-targets; the cutoff is a caller parameter, default 0.75.

---

## Test Datasets

### Dataset: Indel off-target reachable only by a single insertion (hand-derived)

**Source:** derived from the Smith–Waterman recurrence (Smith & Waterman 1981) with BLAST DNA scoring (+2/−3, gap −2, as in `SequenceAligner.BlastDna`); cross-checked by an independent Python re-implementation of the recurrence.

| Parameter | Value |
|-----------|-------|
| Probe | `ACGTACGTACGT` (12 nt) |
| Reference | `NNNNN` + probe + `NNNNNNNNNN` + `ACGTACTGTACGT` + `NNNNN` |
| On-target site | start 5, exact 12/12 match, no gaps → identity 1.0, coverage 1.0 |
| Off-target site | `ACGTACTGTACGT` (probe with a `T` inserted after position 6) at start 27 |
| Off-target gapped alignment | probe `ACGTAC-GTACGT` vs ref `ACGTACTGTACGT` → 12/12 identical aligned columns, one gap |
| Off-target identity | 12 / 12 = 1.0 (HasGaps = true) |
| Ungapped Hamming scan (maxMismatches = 3) | finds ONLY the on-target at 5; every fixed 12-window over the indel region has ≥ 6 mismatches → misses the off-target |
| Gapped scan | finds both: 1 on-target (no gap) + 1 off-target (with gap) |

### Dataset: Indel + mismatch off-target, identity < 1 (hand-derived)

**Source:** same recurrence/scoring as above.

| Parameter | Value |
|-----------|-------|
| Probe | `ACGTACGTACGT` (12 nt) |
| Off-target region | `ACGTACTGTACTT` (insert `T` after position 6; trailing `GTACGT`→`GTACTT`) |
| Local alignment | probe `ACGTAC-GTAC` vs ref `ACGTACTGTAC` — SW trims the mismatched `TT` tail (zero-floor) |
| Identical aligned columns | 10 |
| Identity | 10 / 12 = 0.8333… |
| Threshold behaviour | with minIdentity = 0.75 → called as an off-target; with minIdentity = 0.90 → rejected |

---

## Assumptions

1. **ASSUMPTION: On-target = first perfect ungapped full-coverage exact match.** The literature defines specificity as on-target signal vs non-intended (off-target) signal, but does not prescribe an algorithmic on/off label for pooled references. We classify the single perfect (identity 1.0, coverage 1.0, no gaps) exact match as the intended on-target; additional perfect repeats and all imperfect/indel hits are off-targets. This is the most defensible operationalization (the intended hybridization site is the exact complement) and is exposed transparently in the result record. It is API/labelling, not a sourced numeric constant.

---

## Recommendations for Test Coverage

1. **MUST Test:** an off-target reachable only via a single indel is found by the gapped scan but missed by the ungapped `ValidateProbe` Hamming scan (maxMismatches = 3) — Evidence: Altschul et al. 1990 (gapped vs ungapped), hand-derived dataset.
2. **MUST Test:** the perfect on-target exact match is classified as on-target and excluded from the off-target count — Evidence: on/off separation; Kane et al. 2000 (intended vs non-intended signal).
3. **MUST Test:** exact identity and coverage values on a hand-derived gapped alignment (1.0 indel hit; 0.8333 indel+mismatch hit) — Evidence: Smith & Waterman 1981 recurrence.
4. **SHOULD Test:** the identity threshold gates hits (0.8333 hit admitted at 0.75, rejected at 0.90) — Rationale: Kane et al. 2000 threshold, caller-configurable.
5. **COULD Test:** null/empty probe and null references guard behaviour — Rationale: API robustness.

---

## References

1. Smith TF, Waterman MS (1981). Identification of common molecular subsequences. *Journal of Molecular Biology* 147(1):195–197. https://doi.org/10.1016/0022-2836(81)90087-5 (recurrence retrieved via https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
2. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990). Basic local alignment search tool. *Journal of Molecular Biology* 215(3):403–410. https://doi.org/10.1016/S0022-2836(05)80360-2 (retrieved via https://en.wikipedia.org/wiki/BLAST_(biotechnology))
3. Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Research* 25(17):3389–3402. https://doi.org/10.1093/nar/25.17.3389
4. Kane MD, Jatkoe TA, Stumpf CR, Lu J, Thomas JD, Madore SJ (2000). Assessment of the sensitivity and specificity of oligonucleotide (50mer) microarrays. *Nucleic Acids Research* 28(22):4552–4557. https://pmc.ncbi.nlm.nih.gov/articles/PMC113865/

---

## Change History

- **2026-06-24**: Initial Evidence for the gapped (Smith–Waterman) off-target scan + on/off-target separation (limitation fix). The prior ungapped-Hamming validation evidence is preserved in the TestSpec/algorithm doc.
