# Evidence Artifact: COMPGEN-ANI-001

**Test Unit ID:** COMPGEN-ANI-001
**Algorithm:** Average Nucleotide Identity (ANI), ANIb definition (Goris et al. 2007)
**Date Collected:** 2026-06-14 (refreshed 2026-06-23 for gapped placement + reciprocal ANI)

---

## Online Sources

### Goris et al. 2007 — DNA-DNA hybridization values and their relationship to whole-genome sequence similarities

**URL:** https://www.microbiologyresearch.org/content/journal/ijsem/10.1099/ijsem.0.000760 (OrthoANI, which reproduces and cites the Goris method); primary citation Int J Syst Evol Microbiol 57:81-91, DOI:10.1099/ijs.0.64483-0
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary paper, IJSEM)

**Retrieval:** Web search query `Goris 2007 ANI "1020 nt" OR "1020 bp" fragments BLASTN "two genomic sequences was calculated" identity coverage 30% 70% methods`, which returned the verbatim Methods passages of Goris et al. 2007 (reproduced in the secondary literature and pyani docs).

**Key Extracted Points:**

1. **Fragmentation:** The query genome "was cut into consecutive 1020 nt fragments, and the 1020 nt cut-off was used to correspond with the fragmentation of the genomic DNA to approximately 1 kb fragments during the DDH experiments." (verbatim from retrieved text)
2. **Search:** "The 1020 nt fragments were then used to search against the whole genomic sequence of the other genome in the pair (the reference) by using the BLASTN algorithm; the best BLASTN match was saved for further analysis."
3. **BLAST settings (for completeness):** "X=150 ... q=-1 ... F=F ... the rest of the parameters were used at the default settings."
4. **ANI definition / cut-offs (verbatim):** "ANI between the query genome and the reference genome was calculated as the mean identity of all BLASTN matches that showed more than 30 % overall sequence identity (recalculated to an identity along the entire sequence) over an alignable region of at least 70 % of their length."
5. **Species boundary:** "ANI values of approximately 95% correspond to the 70% DNA-DNA hybridization standard for defining a species." (retrieved in the Goris-2007 search summary)

### Konstantinidis & Tiedje 2005 — Genomic insights that advance the species definition for prokaryotes

**URL:** https://www.pnas.org/doi/abs/10.1073/pnas.0409727102 (PNAS 102(7):2567-2572, DOI:10.1073/pnas.0409727102)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary paper, PNAS)

**Retrieval:** Web search query `Konstantinidis Tiedje 2005 genomic insights species definition prokaryotes ANI 95% nucleotide identity`.

**Key Extracted Points:**

1. **ANI as a relatedness measure:** "the average nucleotide identity (ANI) of the shared genes between two strains is a robust means to compare genetic relatedness among strains" (retrieved summary).
2. **Species threshold:** "ANI values of ≈94% corresponded to the traditional 70% DNA–DNA reassociation standard of the current species definition."

### pyani (Pritchard et al.) — reference implementation of ANIb

**URL:** https://pyani.readthedocs.io/en/latest/api/pyani.anib.html and https://github.com/widdowquinn/pyani/blob/master/README.md
**Accessed:** 2026-06-14
**Authority rank:** 3 (established reference implementation)

**Retrieval:** WebFetch of the pyani README and anib API page.

**Key Extracted Points:**

1. **Fragment length:** pyani's ANIb fragments the input sequences into "1020nt fragments" and aligns them with BLASTN+.
2. **Cut-offs (confirming Goris):** "Coverage filter: Alignments must span at least 70% of the query fragment length"; "Identity filter: Matches must show at least 30% nucleotide identity 'recalculated to an identity along the entire sequence' (i.e., calculated over the full query fragment length, not just the aligned region)."
3. **Averaging:** ANIb is "the percentage nucleotide identity of the matching regions, as an average for all matching regions" — i.e. the mean per-fragment identity across qualifying fragments.

### Refresh 2026-06-23 — gapped placement + reciprocal ANI (this session's retrieval)

**Retrieval A:** WebSearch query `Goris 2007 ANI BLASTN 1020 nt fragments "30%" identity "70%" alignable region mean identity reciprocal both directions average`. Returned verbatim Goris Methods text plus the reciprocal statement.

**Extracted (verbatim from the returned Goris Methods passages):**
- "The ANI between the query and the reference genome was calculated as the mean identity of all BLASTn matches that showed more than 30% overall sequence identity (recalculated to an identity along the entire sequence) over an alignable region of at least 70%."
- **Reciprocal:** "reverse searching, i.e., in which the reference genome is used as the query, was also performed to provide reciprocal values." → confirms the bidirectional computation; the symmetric ANIb value is the mean of the two directions.

**Retrieval B:** WebFetch of `https://pyani.readthedocs.io/en/latest/_modules/pyani/anib.html` (pyani ANIb source, reference implementation, authority rank 3).

**Extracted (from the fetched source/page):**
- Fragment length: query "cut into consecutive 1020 nt fragments".
- **Gapped:** the blastn command includes `-xdrop_gap_final 150`, and the code computes `ani_alnlen = blast_alnlen - blast_gaps` — i.e. BLASTN is **gapped** and gaps are explicitly accounted for.
- **Identity recalculated over query length:** `ani_pid = ani_alnids / qlen` (identical aligned bases / query fragment length).
- **Coverage over query length:** `ani_coverage = ani_alnlen / qlen` (ungapped aligned length / query fragment length).
- **Both directions:** "Reverse searching, i.e. in which the reference genome is used as the query, was also performed to provide reciprocal values"; the code "generates paired jobs in both directions".

These confirm: (1) the canonical method is **gapped**; (2) identity AND coverage are recalculated over the full query-fragment length; (3) the search is performed in **both directions** (the reciprocal/symmetric ANIb value is the mean of the two).

---

## Documented Corner Cases and Failure Modes

### From Goris et al. 2007 / pyani

1. **Non-conserved fragments are discarded:** fragments whose best match falls below 30 % identity or below 70 % alignable length do NOT contribute to the mean. ANI is computed only over conserved (alignable) fragments, not over the whole genome.
2. **Identity recalculated over the whole fragment:** per-fragment identity is matching bases over the fragment length, not over only the locally aligned sub-region. A short, very high-identity local hit still has low recalculated identity.
3. **Asymmetry / reciprocal:** ANI(query→reference) need not equal ANI(reference→query) because the query is the genome that is fragmented (pyani notes "non-symmetrical result matrices"). Goris perform "reverse searching ... to provide reciprocal values"; the symmetric ANIb value is the mean of the two directions. `CalculateReciprocalAni` implements this mean; the single-direction value is well defined on its own.
4. **Gapped placement (indels):** the best BLASTN match is gapped (pyani `-xdrop_gap_final 150`, `ani_alnlen = blast_alnlen - blast_gaps`). An ungapped placement underestimates identity for indel-containing homologous regions; the gapped path (`gapped: true`) recovers it.

---

## Test Datasets

### Dataset: Synthetic exact-arithmetic fragments

**Source:** Derived directly from the Goris et al. 2007 formula (mean of per-fragment identities = matching bases / fragment length, with the >30 % identity and ≥70 % alignable cut-offs). All values are computed by hand and independently re-derived in Python; they do not depend on the implementation.

| Parameter | Value |
|-----------|-------|
| Reference R | `AAAACCCCGGGGTTTT` (16 nt) |
| Fragment length | 4 |
| Query = R (identical) | ANI = 1.0 (four substrings, each identity 1.0) |
| Query `AAAACCCCGGGGTTTA` (1 mismatch in last frag) | last frag `TTTA` best = 3/4 = 0.75 → ANI = (1+1+1+0.75)/4 = 0.9375 |
| Query `AAAACCCCGGGGAATT` (last frag `AATT`) | best vs `TTTT` = 2/4 = 0.5 → ANI = (1+1+1+0.5)/4 = 0.875 |
| Identity cut-off: query `AAAACGTC`, ref `AAAAAAAA` | frag1 `AAAA`=1.0 kept; frag2 `CGTC`=0/4=0.0 NOT >0.30 → excluded → ANI = 1.0 |
| Alignable cut-off: query `AAAA`, ref `AA` (ref < frag) | no full-length placement, alignable fraction 0 < 0.70 → no qualifying frag → ANI = 0 |
| Query shorter than fragment | no fragment of length 4 fits → ANI = 0 |

### Dataset: Gapped placement recovers an indel (hand-derived alignment)

**Source:** Derived by hand from the Goris/pyani definition (identity = identical aligned columns / query-fragment length, gapped Smith-Waterman placement). Independently re-derivable.

| Parameter | Value |
|-----------|-------|
| Query fragment | `AAAACCCC` (8 nt) |
| Reference | `AAAATCCCC` (9 nt; one `T` inserted after `AAAA`) |
| Ungapped best offset | best full-length 8-mer window matches 7/8 bases → identity 0.875 |
| Gapped alignment | `AAAA-CCCC` vs `AAAATCCCC`: 8 identical columns, 1 gap → identity 8/8 = 1.0; coverage 8/8 = 1.0 (≥ 0.70) → qualifies |
| Expected | gapped ANI = 1.0 > ungapped ANI = 0.875 |

### Dataset: Reciprocal (bidirectional) ANI

**Source:** Goris "reverse searching ... to provide reciprocal values"; reciprocal ANI = mean of the two single-direction ANIs.

| Parameter | Value |
|-----------|-------|
| Identical genomes (A=B) | both directions 1.0 → reciprocal = 1.0 |
| Symmetry | reciprocal ANI(A,B) = reciprocal ANI(B,A) by construction (mean is order-independent) |
| A `AAAACGTC`, B `AAAAAAAA`, fragLen 4 | A→B = 1.0 (CGTC excluded), B→A = 1.0 (both AAAA frags match) → reciprocal = (1.0+1.0)/2 = 1.0 |

---

## Assumptions

**No unresolved correctness-affecting assumptions remain.** The earlier "ungapped placement" assumption is RESOLVED by this fix: gapped Smith-Waterman placement is now implemented (`gapped: true`) per Goris/pyani, and the reciprocal mean is implemented (`CalculateReciprocalAni`). The remaining decision is an engine choice, not a correctness gap:

1. **DECISION (not an assumption): the gapped path uses the library's own Smith-Waterman aligner (`SequenceAligner.LocalAlign`, BLAST DNA scoring) rather than the NCBI BLASTN engine.** Smith-Waterman is full dynamic programming (more sensitive than BLAST's heuristic seeding, not less); identity and coverage are computed by the same recalculated-over-fragment definition (pyani `ani_pid = ani_alnids/qlen`, `ani_coverage = ani_alnlen/qlen`). Numeric ANI on real genomes may differ slightly from NCBI-BLASTN pipelines because the alignment engine differs; the indel-handling behaviour is correct. Documented in algorithm doc §5.3.

---

## Recommendations for Test Coverage

1. **MUST Test:** Identical genomes give ANI = 1.0 (each fragment is a perfect substring). — Evidence: Goris 2007 mean-identity formula.
2. **MUST Test:** One substituted base in one fragment lowers that fragment's identity to (frag−1)/frag and ANI to the exact mean (worked example 0.9375). — Evidence: Goris 2007 "recalculated to an identity along the entire sequence".
3. **MUST Test:** A fragment whose best match is ≤30 % identity is excluded from the mean. — Evidence: Goris 2007 ">more than 30 % overall sequence identity".
4. **MUST Test:** A fragment that cannot align over ≥70 % of its length (reference shorter than fragment) is excluded. — Evidence: Goris 2007 "over an alignable region of at least 70 % of their length".
5. **SHOULD Test:** Consecutive non-overlapping fragmentation (number of fragments = ⌊len/fragmentLength⌋; trailing partial fragment ignored). — Rationale: Goris 2007 "consecutive 1020 nt fragments".
6. **SHOULD Test:** Result is a fraction in [0, 1]. — Rationale: invariant from the identity definition.
7. **COULD Test:** Custom fragmentLength / minIdentity / minAlignableFraction parameters behave per definition. — Rationale: parameter exposure.
8. **MUST Test:** null / empty inputs → 0; non-positive fragmentLength → ArgumentOutOfRangeException. — Rationale: documented validation.
9. **MUST Test (gapped):** identical genomes → 1.0 gapped; an indel fragment where gapped identity (1.0) > ungapped (0.875). — Evidence: gapped dataset above; Goris/pyani gapped BLASTN.
10. **MUST Test (gapped cut-offs):** gapped identity and alignable-fraction cut-offs still gate (low-identity and low-coverage fragments excluded). — Evidence: Goris 30 %/70 % cut-offs over the query length.
11. **MUST Test (reciprocal):** reciprocal ANI of identical genomes = 1.0; reciprocal symmetry ANI(A,B)=ANI(B,A); reciprocal = mean of both directions; null/empty → 0; non-positive fragmentLength throws. — Evidence: Goris reverse searching; reciprocal dataset above.

---

## References

1. Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, Tiedje JM (2007). DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. Int J Syst Evol Microbiol 57(1):81-91. https://doi.org/10.1099/ijs.0.64483-0
2. Konstantinidis KT, Tiedje JM (2005). Genomic insights that advance the species definition for prokaryotes. Proc Natl Acad Sci USA 102(7):2567-2572. https://doi.org/10.1073/pnas.0409727102
3. Lee I, Ouk Kim Y, Park S-C, Chun J (2016). OrthoANI: An improved algorithm and software for calculating average nucleotide identity. Int J Syst Evol Microbiol 66(2):1100-1103. https://doi.org/10.1099/ijsem.0.000760
4. Pritchard L, Glover RH, Humphris S, Elphinstone JG, Toth IK (2016). Genomics and taxonomy in diagnostics for food security (pyani ANIb). Anal Methods 8:12-24. pyani docs: https://pyani.readthedocs.io/en/latest/api/pyani.anib.html

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-06-23**: Limitation fix — added gapped Smith-Waterman fragment placement (`gapped` option) and reciprocal (bidirectional) ANI (`CalculateReciprocalAni`). Refreshed Goris/pyani retrieval (gapped BLASTN, identity/coverage over query length, reverse searching). Resolved the ungapped-placement assumption.
