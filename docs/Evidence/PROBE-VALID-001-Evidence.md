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

### Karlin & Altschul (1990) / Altschul et al. (1990) — E-value, λ, K, bit score

**URL:** https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html ("The Statistics of Sequence Similarity Scores", Altschul) and http://www.cs.cmu.edu/~durand/03-711/2011/Lectures/Blast-informationContent-2011.pdf (Durand, "BLAST (Karlin–Altschul) Statistics", CMU 03-711, explicitly citing Karlin & Altschul 1990 *PNAS* 87:2264 and Altschul et al. 1990 *J Mol Biol* 215:403)
**Accessed:** 2026-06-24
**Authority rank:** 1 (primary peer-reviewed papers, via NCBI's own tutorial and a university lecture quoting the primaries verbatim)

**Key Extracted Points (verbatim):**

1. **E-value:** "the expected number of HSPs with score at least *S* is given by the formula" **E = K·m·n·e^(−λS)** (Altschul tutorial; Durand renders it `E = Kmn e^−λS`). m = query length, n = database length, S = raw score.
2. **λ defining equation:** "The parameter λ is specified by the equation **1 = Σ p_i p_j e^(−λ S[ij])**" (Durand) — i.e. λ is the unique positive root of Σ_{i,j} p_i p_j e^(λ s_ij) = 1 (p = residue background frequencies, s_ij = score matrix).
3. **Bit score:** "By normalizing a raw score using the formula" **S' = (λS − ln K)/ln 2** "one attains a 'bit score' *S'*, which has a standard set of units" (Altschul tutorial; Durand renders the same `S' = (λS − ln K)/ln 2`).
4. **E from bit score:** "The *E*-value corresponding to a given bit score is simply" **E = m·n·2^(−S')** (both sources).
5. **Scoring-scheme precondition:** "the expected score for aligning a random pair of … is required to be negative. Were this not the case, long alignments would tend to have high score independently of whether the segments aligned were related, and the statistical theory would break down" (Altschul tutorial). The complementary requirement — at least one positive score so the positive root exists — is the standard statement of the same theory ("for valid scoring matrices (ones where at least one positive score exists), λ will have a unique positive solution").
6. **K:** "K is a constant that depends on S[i,j] and can be computed from the theory for any scoring function" (Durand). It is computed exactly as NCBI BLAST+ computes it (`BlastKarlinLHtoK`, see the NCBI BLAST+ source entry below); a caller-supplied K still overrides it.
7. **λ ≈ 1.37, K ≈ 0.711 cross-check (+1/−3, uniform 0.25):** NCBI blastn reports Lambda ≈ 1.37 and K ≈ 0.711 for match=+1/mismatch=−3 (https://www.biostars.org/p/9596760/ shows a blastn run with "matrix:1 -3 … Lambda: 1.37, K: 0.711"). Solving 0.25·e^(λ·1) + 0.75·e^(λ·(−3)) = 1 independently gives **λ = 1.3740631** (re-derived in this session by bisection), matching the published 1.37/1.374. The expected per-pair score is 0.25·1 + 0.75·(−3) = **−2.0 < 0** and a positive score (+1) exists, so both preconditions hold.

### NCBI BLAST+ source (blast_stat.c / blast_setup.c / blast_hits.c / ncbi_math.c) + blastn 2.12.0+

**URL:** https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/blast/core/blast_stat.c (and `blast_setup.c`, `blast_hits.c`, `ncbi_math.c`; `include/algo/blast/core/blast_stat.h`), retrieved 2026-10-01. Oracle binary: NCBI blastn 2.12.0+ (Debian `ncbi-blast+`).

**Authority rank:** 2 (reference implementation of the published method)

**Key Extracted Points:**

1. Ungapped λ, H, K (`Blast_KarlinBlkUngappedCalc`): λ by safeguarded Newton (`Blast_KarlinLambdaNR`), `H = λ·Σ s·p_s·e^{λs}` (`BlastKarlinLtoH`), K by `BlastKarlinLHtoK`: on the lattice reduced by δ = gcd, `K = (p₋₁ − p₁)²/p₋₁` when low = −1 and high = 1; `K = (H/λ)(1 − e^{−λ})` when high = 1; `K = (μ²/(H/λ))(1 − e^{−λ})` when low = −1; else `K = −exp(−2·Σ_j inner_j/j) / ((H/λ)·expm1(−λ))` with `inner_j = Σ_{i<0} P(i,j)e^{λi} + Σ_{i≥0} P(i,j)` over gapless alignments of j pairs (sum limit 10⁻⁴, ≤ 100 terms). Comment example in blast_stat.c: scores −2/0/3 with probabilities 0.7/0.1/0.2 → λ = 0.330, K = 0.154 (Python port: 0.32995, 0.15399).
2. blastn ungapped uses the standard (uniform 0.25) nucleotide composition: a GC-rich query still prints 1.37/0.711/1.31 for +1/−3.
3. Gapped λ/K/H/α/β are the `blastn_values_<reward>_<penalty>` tables (row {open, extend, λ, K, H, α, β, θ}); a leading {0,0} row is the non-affine (megablast greedy) entry (`s_SplitArrayOf8`) — e.g. 2/−3 {0,0} → 0.55/0.21 is NOT the ungapped value; gap costs ≥ (gap_open_max, gap_extend_max) copy the ungapped block; reward/penalty with gcd d > 1 use the reduced table with gap costs × d, λ and α ÷ d; 2/−3, 2/−5, 2/−7, 3/−4 set round_down (E-value from `score & ~1`).
4. α/β: table values; otherwise α = λ_ungapped/H, β = −2 for 1/−1 and 2/−3 else 0 (`s_GetUngappedBeta`). Length adjustment `BLAST_ComputeLengthAdjustment(K, logK, α/λ, β, m, n, N)`; effective search space `(m − ℓ)·max(1, n − N·ℓ)` (`BLAST_CalcEffLengths`); `E = searchsp·exp(−λS + ln K)`; bit score `(S·λ − ln K)/ln 2` on the raw score.
5. NCBI quirk: for a scheme whose scores share a divisor d > 1, `BlastKarlinLHtoK` indexes the probability array by the reduced offset from the unreduced lowest score, so blastn prints K = 1.17 for +4/−6 although K must be scale-invariant (+2/−3: 0.408). Seqeron computes K on the reduced lattice by default (= 0.408 for +4/−6); `KarlinKMethod.NcbiBlast` reproduces the BLAST+ indexing (default in `ComputeBlastnStatistics` / `GetBlastnGappedKarlinParameters`). blastn 2.12.0+ `-task blastn -ungapped` footers (λ K H): 4/−6 0.317 1.17 0.912; 6/−9 and 8/−12 K 1.17; 4/−10 0.340 1.06 1.24; 6/−4 0.136 1.63 0.222; 6/−10 0.216 1.03 0.997; 8/−10 0.151 1.07 0.753; 10/−8 0.0958 1.31 0.357; 4/−2 0.132 0.0532 0.0722; 2/−2 0.549 0.333 0.549; 2/−4 0.666 0.621 1.12; 4/−6 gap 12/8 "Gapped 0.317 1.17 0.912". Python port (verbatim indexing): K 1.1666856431064105, 1.0551429674627688, 1.633485999782881, 1.030863998562112, 1.0682924577828306, 1.3071030460951232, 0.05322292075469216. End-to-end (200-nt query, 3100-nt subject): 4/−6 ungapped space 573996, S 400 → 182 bits, E 6.03e-50 (port 6.034606696310044e-50); 6/−4 ungapped space 425600, S 602 → 117 bits, E 2.47e-30 (port 2.471093453808567e-30). (B07 F54)

### Dataset: NCBI blastn 2.12.0+ statistics (oracle for `ComputeUngappedKarlinParameters`, `GetBlastnGappedKarlinParameters`, `ComputeBlastnStatistics`)

| Case | blastn output | Python port of blast_stat.c |
|---|---|---|
| ungapped 1/−3 | λ 1.37, K 0.711, H 1.31 | 1.3740631224599753, 0.7106027952162398, 1.3072466039090012 |
| ungapped 2/−3 | 0.634, 0.408, 0.912 | 0.6337314430979075, 0.4081456625463167, 0.9124383922742278 |
| ungapped 1/−2 | 1.33, 0.621, 1.12 | 1.3327057628202603, 0.6209911172603866, 1.1240918464926624 |
| gapped 2/−3 5/2 | 0.625, 0.410, 0.780 | table |
| gapped 1/−3 2/2; 1/−2 2/2; 2/−3 4/4 | 1.37/0.700/1.20; 1.33/0.620/1.10; 0.630/0.420/0.840 | table |
| m 40, n 3079, 2/−3 5/2 | eff. space 88972; S 80 → 73.4 bits, 7e-18; S 15 → 5.8 | ℓ 11, 88972, 7.035793990394873e-18, 73.42105622960482; S 15 → 14 → 5.780434617461434 |
| same, N = 5, n = 6040 | 167440; S 80 → 1.32e-17 | ℓ 12, 167440, 1.3240944856266213e-17 |
| m 40, n 3079, 1/−3 2/2 | 98272; S 31 → 2e-14 | ℓ 8, 98272, 2.4719585736391905e-14 |
| m 40, n 3079, ungapped 2/−3 | 95170; S 80 → 74.4 bits, 4e-18 | ℓ 9, 95170, 3.725887650102598e-18 |
| probe 40 nt vs 279-nt subject (1 mismatch + 1-nt deletion), 2/−3 5/2 | score 66, 60.8 bits, 4.33e-15, space 8672 | Biopython local score 66.0; 4.327686048582086e-15, 60.79747462182638 |

### Kane et al. (2000) — 50-mer oligonucleotide microarray specificity

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC113865/ (Kane MD, Jatkoe TA, Stumpf CR, Lu J, Thomas JD, Madore SJ (2000) *Nucleic Acids Research* 28(22):4552–4557)
**Accessed:** 2026-06-24
**Authority rank:** 1 (peer-reviewed paper) / 5 (curated specificity guideline)

**Key Extracted Points:**

1. **Off-target identity threshold:** "for a given oligonucleotide probe any 'non-target' transcripts (cDNAs) **>75% similar** over the 50 base target may show cross-hybridization." → 0.75 identity over the probe length is the empirically-grounded default above which a hit is called an off-target.
2. **Gene-specific rule:** "oligonucleotide probes with **<75% overall sequence similarity** with non-target sequences and <14 contiguous complementary base pairs are gene-specific" — the complement of the off-target call.
3. **Contiguous-stretch caveat:** "if the 50 base target region is marginally similar, it must not include a stretch of complementary sequence >15 contiguous bases."
4. **Implemented decision rule (B07 review, 2026-10-01):** a non-target strand cross-hybridizes when identity (identical columns of the best local alignment ÷ probe length) **> 0.75** or the longest identical stretch **> 15 nt** — the reading used by later microarray-design pipelines ("a probe is likely to cross-hybridize with a nontarget if overall sequence identity is > 75% or if there is a contiguous match > 15 bp", Satya RV, Zavaljevski N, Kumar K, Reifman J (2008) BMC Bioinformatics 9:185; restated by Chen & Sharp (2002) Oliz, BMC Bioinformatics 3:27). Sources opened: WebSearch snippets of the Kane abstract (academic.oup.com, PubMed 11071945), Satya et al. 2008 and Oliz (publisher/PMC pages blocked). The paper's conclusion "<14 contiguous ... gene-specific" leaves 14–15 nt as a grey zone; `maxContiguousMatch` is configurable.

### Primer3 hybridization-probe (internal-oligo) self-structure screen

**Source:** primer3 `libprimer3.cc` (`o_args`: PRIMER_INTERNAL_MAX_SELF_ANY_TH = PRIMER_INTERNAL_MAX_SELF_END_TH = PRIMER_INTERNAL_MAX_HAIRPIN_TH = 47 °C; internal-oligo conditions 50 nM DNA, 50 mM monovalent, 0 Mg²⁺, 0 dNTP; `oligo_compl_thermod`, `oligo_hairpin`), opened in PROBE-DESIGN-001 (B07 F18). Reference: primer3-py 2.3.1 `calc_homodimer` / `calc_end_stability` / `calc_hairpin`.

### OligoArray 2.0 (Rouillard, Zuker & Gulari 2003, NAR 31:3057) — duplex-Tm specificity

WebSearch extract: specificity is computed from the thermodynamics of hybridization of the probe with every BLAST hit; "if there is no possible cross-hybridization with a Tm above the specificity threshold set by the user, the oligonucleotide is considered to be specific". → `CrossHybridizationAssessment.DuplexTm` (ntthal THAL_ANY Tm of the probe with the complementary strand of the aligned site) + optional `maxDuplexTm`.

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

### Dataset: Karlin–Altschul λ / bit-score / E-value worked example (hand-derived)

**Source:** Karlin & Altschul (1990); Altschul et al. (1990) formulas above; λ solved in-session, scalars cross-checked in Python.

| Parameter | Value |
|-----------|-------|
| Scoring scheme | match = +1, mismatch = −3, base freq = 0.25 (NCBI blastn +1/−3) |
| λ (root of 0.25·e^λ + 0.75·e^(−3λ) = 1) | 1.3740631224599755 (≈ published 1.37) |
| K (caller-supplied in this example; computed default 0.7106027952162398) | 0.711 |
| Raw score S | 30 |
| Query length m | 20 |
| Database length n | 1000 |
| Bit score S' = (λ·30 − ln 0.711)/ln 2 | 59.962700114285006 |
| E = K·m·n·e^(−λ·30) = m·n·2^(−S') | 1.7801583686083893e−14 |
| Monotonicity | E(S=31) = 4.5052e−15 < E(S=30) (decreases with score) |
| Linear in m·n | E(n=2000) = 2 × E(n=1000) |

### Dataset: Kane criteria (Biopython 1.88 PairwiseAligner local, match 2 / mismatch −3 / open −7 / extend −2 = BLAST+ blastn 2/−3/5/2)

Probe `TATGCCTCCGGTACATCAACTACAGTTAGCCTTAAGAGAAAAATCCCAAA` (random.seed 2000).

| Non-target | Strand | Score | Identical / 50 | Longest contiguous | Kane |
|---|---|---:|---:|---:|---|
| A (substitution every 5th base, random flanks) | fwd | 53 | 40 (0.80) | 6 | identity |
| B (probe[10..28) embedded) | fwd | 37 | 28 (0.56) | 18 | contiguous |
| C (unrelated) | fwd / rc | 10 / 16 | 5 / 8 | 5 / 8 | — |
| D (revcomp(probe) embedded) | rc | 100 | 50 (1.00) | 50 | both |
| E (probe[0..15) embedded) | fwd | 30 | 15 (0.30) | 15 | — (15 is not > 15) |

Site duplex Tm (primer3-py `calc_heterodimer(probe, revcomp(site))`, mv 50, dv 0, dntp 0, dna 50): A fwd site 20..68 → 36.11423712379826 °C; D rc site 15..64 → 66.04038852959525 °C.

Random cross-check (this review): 420 probe/non-target pairs (20–70-nt probes, 0–9000-nt non-targets incl. mutated/indel/reverse-complement copies and chunk-boundary cases), 840 strands — alignment score and longest contiguous match identical to Biopython / DP LCS on all 840; reported identity always one of Biopython's co-optimal alignments' identities (48 strands have co-optimal alignments with different identities); duplex Tm identical to primer3-py `calc_heterodimer` on all 612 ≤ 60-nt cases (max |Δ| = 0).

### Dataset: ntthal self-structure (primer3-py 2.3.1, mv 50, dv 0, dntp 0, dna 50)

| Probe | calc_homodimer Tm | calc_end_stability Tm | calc_hairpin Tm | fold-back fraction |
|---|---:|---:|---:|---:|
| GCGCGCGCGCGCGCGCGCGC | 78.85652531616256 | 78.85652531616256 | 87.30265612393043 | 1.00 |
| ACGTACGTACGTACGTACGTACGT | 59.1857717189107 | 59.1857717189107 | 67.29188756961071 | 1.00 |
| CTAGAAATGCTGTCGGGACTTCTAC | −6.43 (→ 0) | −99.94 (→ 0) | 0 | 0.64 |
| GCGCGCGCGC | 52.763 | 52.763 | 55.851 | 1.00 |

At mv 100, dv 2, dntp 0.2, dna 250: ACGTACGTACGTACGTACGTACGT → 69.17069845823409 / 69.17069845823409 / 74.99462150250321.

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
6. **MUST Test:** the computed λ for +1/−3 with uniform 0.25 frequencies equals the published 1.374 to ≤1e-6 — Evidence: Karlin & Altschul 1990; NCBI blastn cross-check. A wrong solver fails this.
7. **MUST Test:** bit score and E-value for a hand-derived (S, m, n) match S' = (λS − ln K)/ln 2 and E = K·m·n·e^(−λS); E = m·n·2^(−S') equals the K·m·n·e^(−λS) form — Evidence: Altschul et al. 1990.
8. **MUST Test:** E decreases as S increases and scales linearly with m·n — Evidence: E = K·m·n·e^(−λS).

---

## References

1. Smith TF, Waterman MS (1981). Identification of common molecular subsequences. *Journal of Molecular Biology* 147(1):195–197. https://doi.org/10.1016/0022-2836(81)90087-5 (recurrence retrieved via https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
2. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990). Basic local alignment search tool. *Journal of Molecular Biology* 215(3):403–410. https://doi.org/10.1016/S0022-2836(05)80360-2 (retrieved via https://en.wikipedia.org/wiki/BLAST_(biotechnology))
3. Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Research* 25(17):3389–3402. https://doi.org/10.1093/nar/25.17.3389
4. Kane MD, Jatkoe TA, Stumpf CR, Lu J, Thomas JD, Madore SJ (2000). Assessment of the sensitivity and specificity of oligonucleotide (50mer) microarrays. *Nucleic Acids Research* 28(22):4552–4557. https://pmc.ncbi.nlm.nih.gov/articles/PMC113865/
5. Karlin S, Altschul SF (1990). Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes. *PNAS* 87(6):2264–2268. https://doi.org/10.1073/pnas.87.6.2264 (formulas retrieved via the NCBI tutorial and CMU 03-711 lecture cited above)
6. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990). Basic local alignment search tool. *J Mol Biol* 215(3):403–410. https://doi.org/10.1016/S0022-2836(05)80360-2 (E-value/bit-score statements retrieved via https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html and http://www.cs.cmu.edu/~durand/03-711/2011/Lectures/Blast-informationContent-2011.pdf)

---

## Change History

- **2026-10-02** (B07 audit round 2, A6): the fallback self-dimer criterion (> 60 nt, non-ACGT, `Heuristic`) is Primer3's alignment-mode internal-oligo `oligo_compl` (dpal `self_any` DPAL_LOCAL / `self_end` DPAL_GLOBAL_END, PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END = 12.00; libprimer3.cc / dpal.c from raw.githubusercontent.com/primer3-org/primer3). Datasets: dpal.c compiled with Primer3's `align()` — CTAGAAATGCTGTCGGGACTTCTAC 9.00 / 7.00, (ACGT)16 64.00 / 64.00, 80-nt stem-loop GGATCACAG…GATCC 60.00, random 80-mer CCCTGAGTCC…GGTTCA 7.00 / 1.00; 40 000/40 000 random values identical. The fold-back fraction is no longer a criterion.
- **2026-06-24**: Initial Evidence for the gapped (Smith–Waterman) off-target scan + on/off-target separation (limitation fix). The prior ungapped-Hamming validation evidence is preserved in the TestSpec/algorithm doc.
- **2026-10-08** (B07 audit round 3, A3-12/A3-14): no published source for `SpecificityScore` = 1/N or `maxMismatches` = 3 (searched: Kane et al. 2000, OligoArray 2.0, Li & Stormo 2001 abstract, Primer3 `libprimer3.cc` `o_args.max_repeat_compl` = PRIMER_INTERNAL_MAX_LIBRARY_MISHYB 12.00) → documented as library conventions; the off-target issue now counts only sites meeting the Kane criteria on their ungapped diagonal. Dataset (Python brute-force Hamming oracle): `GATCCGACGCTA` vs `TTTTTGATCCGACGCTATTTTTGCTCCTACGGTATTTTT` → hits 2 (0, 3 mismatches: 0.75, run 3) → 1 Kane site; + `GGGGGGAACCGAGGCTAGGGGG` (2 mismatches, 0.8333) → 2. 40-mer `TATGCCGTACAGTTTTAAGATAGAGCGAAAGCGCAGACAA`, radius 12: exact, clustered (0.70, run 28), spread (0.70, run 5) → 2 Kane sites. MCP conditions: primer3-py 2.3.1 GCGC…(20) at mv 100 / dv 2 / dntp 0.2 / dna 250 → 88.57933031369095 / 88.57933031369095 / 93.4845818217363.
- **2026-10-01** (B07 PROBE-VALID-001 review): Kane contiguous-stretch criterion + strict > 75 % identity (`AssessCrossHybridization`, both strands), Primer3 ntthal self-structure screen in `ValidateProbe`, OligoArray-style site duplex Tm, `CheckSpecificity` both-strand option; datasets above.
- **2026-06-24**: Added the Karlin–Altschul E-value / bit-score / λ evidence (sources 5–6), the +1/−3 λ≈1.374 cross-check, and the worked-example dataset, for the opt-in `ComputeLambdaNucleotide` / `ComputeKarlinAltschul` statistics.
- **2026-10-01** (B07 PROBE-EVALUE-001 review): K computed as NCBI BLAST+ computes it (no longer the +1/−3 constant for every scheme); NCBI BLAST+ gapped tables, length adjustment and blastn E-values; NCBI source + blastn 2.12.0+ datasets above.
