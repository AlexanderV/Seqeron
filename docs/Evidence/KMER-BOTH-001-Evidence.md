# Evidence Artifact: KMER-BOTH-001

**Test Unit ID:** KMER-BOTH-001
**Algorithm:** K-mer counting over both strands of double-stranded DNA (forward + reverse-complement)
**Date Collected:** 2026-06-14

---

## Online Sources

### kPAL — Methodology (k-mer profile analysis library)

**URL:** https://kpal.readthedocs.io/en/latest/method.html
**Accessed:** 2026-06-14 (WebFetch of the page)
**Authority rank:** 3 (reference implementation / project documentation; the package is described in a peer-reviewed paper — Anvar et al. 2014)

**Key Extracted Points:**

1. **Both-strand balancing (verbatim):** "kPAL can forcefully balance the k-mer profiles (if desired) by adding the values of each k-mer to its reverse complement." This is the additive both-strand operation: the both-strand count of a k-mer is the sum of its own count and the count of its reverse complement.
2. **Purpose:** balancing "enforce[s] balance between sequence information from the minus or plus strand" — i.e. it makes the profile strand-symmetric, which is the both-strand view of double-stranded DNA.

### Anvar et al. (2014) — Determining the quality and complexity of NGS data (kPAL paper)

**URL:** https://link.springer.com/article/10.1186/s13059-014-0555-3 (search-result summary; full text behind Springer IDP redirect, summary retrieved via WebSearch)
**Accessed:** 2026-06-14 (WebSearch summary; primary URL redirected to authentication)
**Authority rank:** 1 (peer-reviewed, Genome Biology 2014, 15:555)

**Key Extracted Points:**

1. **Balance = sum of k-mer and reverse complement:** "The balance function uses a sum of k-mers and their reverse complements to enforce balance between sequence information from the minus or plus strand." This grounds the additive both-strand semantics in the peer-reviewed source.

### Shporer et al. (2016) — Inversion symmetry of DNA k-mer counts (generalized second Chargaff rule)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5006273/
**Accessed:** 2026-06-14 (WebFetch)
**Authority rank:** 1 (peer-reviewed; PMC full text)

**Key Extracted Points:**

1. **Inversion symmetry (verbatim):** "the counts of any string of nucleotides of length k on a single chromosomal strand equal the counts of its inverse (reverse-complement) k-mer."
2. **Strand reading (verbatim):** "the number of times a string of nucleotides of length k is observed on a strand, when read from 5' to 3', is almost equal to the number of times it is observed on the other strand when the latter is read from its 5' end to 3' end." This justifies the identity: occurrences of w on the reverse-complement strand (read 5'→3') = occurrences of RC(w) on the forward strand. Hence both-strand count[w] = forward[w] + forward[RC(w)].

### Marçais & Kingsford (2011) — Jellyfish: a fast k-mer counter (canonical k-mer contrast)

**URL:** https://academic.oup.com/bioinformatics/article/27/6/764/234905 ; manual https://vcru.wisc.edu/simonlab/bioinformatics/programs/jellyfish/jellyfish-manual-1.1.pdf
**Accessed:** 2026-06-14 (WebFetch of the Oxford article page and the manual PDF)
**Authority rank:** 1 (peer-reviewed; Bioinformatics 27(6):764–770)

**Key Extracted Points:**

1. **Definition of the k-mer counting problem (verbatim):** Jellyfish counts "the number of occurrences of every k-mer (substring of length k) in a long string." This is the single-strand counting primitive (`CountKmers`) on which both-strand counting builds.
2. **Canonical option contrast:** Jellyfish offers a canonical (`-C`) mode that collapses a k-mer and its reverse complement onto one representative. KMER-BOTH-001 is NOT canonical collapsing — it is the additive (kPAL "balance") both-strand profile that keeps a key per observed k-mer. (Recorded to distinguish the two strand-aware semantics; the canonical mode wording itself was not extractable from the man-page PDF stream.)

### Mash issue #45 / Ondov et al. — canonical k-mer definition (contrast reference)

**URL:** https://github.com/marbl/Mash/issues/45
**Accessed:** 2026-06-14 (WebFetch)
**Authority rank:** 3 (project maintainer explanation citing the Mash paper)

**Key Extracted Points:**

1. **Canonical definition (verbatim from the Mash paper):** "only the lexicographically smaller of the forward and reverse complement representations of a k-mer is hashed." Confirms the canonical (collapsing) approach, which is explicitly NOT the additive both-strand approach implemented here.

### BioInfoLogics — k-mer counting, part I (Clavijo, 2018)

**URL:** https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
**Accessed:** 2026-06-14 (WebFetch)
**Authority rank:** 4 (tutorial citing primaries; used only for the strand rationale)

**Key Extracted Points:**

1. **Strand rationale (verbatim):** "DNA['s] double-stranded nature means researchers encounter sequences on either strand." Canonical counting "means if GTCGAT appears it will be counted as ATCGAC." Used to confirm the biological motivation for strand-aware counting (the both-strand sum here is the non-collapsing variant).

---

## Documented Corner Cases and Failure Modes

### From Marçais & Kingsford (2011) / Shporer et al. (2016)

1. **Window count:** A length-L sequence has L − k + 1 overlapping k-mers; each strand contributes L − k + 1 windows, so the both-strand grand total over all keys is 2·(L − k + 1). (Derived from the k-mer-occurrence definition above.)
2. **Palindromic k-mers (reverse-complement palindromes):** A k-mer equal to its own reverse complement (e.g. AT, GC, ACGT) receives forward and reverse-complement contributions on the SAME key, so its both-strand count = forward[w] + forward[w]. (Direct consequence of inversion symmetry with RC(w)=w.)
3. **k > L:** L − k + 1 ≤ 0 ⇒ no windows ⇒ empty result on each strand ⇒ empty combined result. (k-mer-occurrence definition.)

---

## Test Datasets

### Dataset: Both-strand worked example ATGGC, k=2

**Source:** Derived from kPAL balance (Anvar et al. 2014) + inversion symmetry (Shporer et al. 2016): count[w] = forward[w] + forward[RC(w)].

| Quantity | Value |
|----------|-------|
| Sequence | ATGGC (L=5) |
| k | 2 |
| Forward 2-mers | AT, TG, GG, GC |
| RC(ATGGC) | GCCAT |
| RC-strand 2-mers | GC, CC, CA, AT |
| Combined counts | AT:2, TG:1, GG:1, GC:2, CC:1, CA:1 |
| Grand total | 8 = 2·(5−2+1) |

### Dataset: Palindromic homopolymer/2-mer example ACGT, k=2

**Source:** Derived; ACGT is a reverse-complement palindrome (RC(ACGT)=ACGT).

| Quantity | Value |
|----------|-------|
| Sequence | ACGT |
| k | 2 |
| Forward 2-mers | AC, CG, GT |
| RC-strand 2-mers | AC, CG, GT |
| Combined counts | AC:2, CG:2, GT:2 |

### Dataset: Non-palindromic homopolymer AAA, k=2

**Source:** Derived; RC(AAA)=TTT.

| Quantity | Value |
|----------|-------|
| Sequence | AAA |
| k | 2 |
| Combined counts | AA:2, TT:2 |

---

## Assumptions

1. **ASSUMPTION: Empty/short input returns empty dictionary** — No authoritative source explicitly defines the both-strand result for an empty sequence or k > L; resolved consistently with the k-mer-occurrence definition (L − k + 1 ≤ 0 ⇒ no windows) and with sibling `CountKmers` behavior in this repository. Non-correctness-affecting beyond the empty-result boundary already implied by the window formula.
2. **ASSUMPTION: k ≤ 0 throws `ArgumentOutOfRangeException`** — Sources define k as a positive substring length but do not prescribe an exception type; resolved to match the sibling `CountKmers`/`GenerateAllKmers` contract in this repository (API-shape only, not correctness-affecting on valid input).

---

## Recommendations for Test Coverage

1. **MUST Test:** Worked example ATGGC, k=2 ⇒ exact dictionary {AT:2,TG:1,GG:1,GC:2,CC:1,CA:1}. — Evidence: kPAL balance (Anvar et al. 2014) + inversion symmetry (Shporer et al. 2016).
2. **MUST Test:** Palindromic ACGT, k=2 ⇒ {AC:2,CG:2,GT:2} (each key doubled). — Evidence: inversion symmetry with RC(w)=w.
3. **MUST Test:** Non-palindromic AAA, k=2 ⇒ {AA:2,TT:2}. — Evidence: forward[w]+forward[RC(w)].
4. **MUST Test:** Grand total = 2·(L − k + 1). — Evidence: k-mer window count (Marçais & Kingsford 2011).
5. **MUST Test:** Identity count[w] = forward[w] + forward[RC(w)] for every key. — Evidence: Shporer et al. 2016.
6. **SHOULD Test:** Case-insensitivity (lowercase input == uppercase). — Rationale: sibling `CountKmers` upper-cases; both-strand must agree.
7. **SHOULD Test:** k = L (single window per strand). — Rationale: boundary of the window formula.
8. **SHOULD Test:** DnaSequence overload delegates to the string overload (smoke). — Rationale: delegate wrapper.
9. **COULD Test:** Empty / null sequence ⇒ empty; k > L ⇒ empty; k ≤ 0 ⇒ throws. — Rationale: documented edge/failure modes.

---

## References

1. Anvar SY, et al. (2014). Determining the quality and complexity of next-generation sequencing data without a reference genome. Genome Biology, 15:555. https://doi.org/10.1186/s13059-014-0555-3
2. kPAL documentation — Methodology. https://kpal.readthedocs.io/en/latest/method.html (accessed 2026-06-14)
3. Shporer S, Chor B, Rosset S, Horn D (2016). Inversion symmetry of DNA k-mer counts: validity and deviations. BMC Genomics. https://pmc.ncbi.nlm.nih.gov/articles/PMC5006273/
4. Marçais G, Kingsford C (2011). A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Bioinformatics, 27(6):764–770. https://doi.org/10.1093/bioinformatics/btr011
5. Ondov BD, et al. — Mash; maintainer explanation of canonical k-mers, GitHub issue #45. https://github.com/marbl/Mash/issues/45 (accessed 2026-06-14)
6. Clavijo BJ (2018). BioInfoLogics — k-mer counting, part I: Introduction. https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/ (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
