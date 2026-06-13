# Evidence Artifact: ANNOT-REPEAT-001

**Test Unit ID:** ANNOT-REPEAT-001
**Algorithm:** Repetitive Element Detection and Classification (tandem repeats, inverted repeats, repeat-class assignment)
**Date Collected:** 2026-06-13

---

## Online Sources

### Wikipedia — Tandem repeat (citing primary sources Duitama et al. 2014; MeSH; Jorda et al. 2010)

**URL:** https://en.wikipedia.org/wiki/Tandem_repeat
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primaries)

**Key Extracted Points:**

1. **Definition (verbatim):** "A pattern of one or more nucleotides is repeated and the repetitions are directly adjacent to each other" (example given: ATTCG ATTCG ATTCG). The defining structural feature is a **head-to-tail consecutive arrangement** — repetitions occur immediately adjacent with no intervening sequence.
2. **Microsatellite / STR size:** the repeated unit is short ("fewer than 10 nucleotides repeated per unit"); microsatellites/STRs are commonly defined as motifs 1–6 bp.
3. **Minisatellite size:** "Between 10 and 60 nucleotides are repeated."
4. **Minimum copies:** a tandem repeat is defined by "two or more" directly adjacent repetitions (a single occurrence of a motif is not a tandem repeat).

### Wikipedia — Inverted repeat (citing Ussery et al. 2008; Ye et al. 2014)

**URL:** https://en.wikipedia.org/wiki/Inverted_repeat
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primaries)

**Key Extracted Points:**

1. **Definition (verbatim):** an inverted repeat is "a single stranded sequence of nucleotides followed downstream by its reverse complement."
2. **Spacer/loop:** "The intervening sequence of nucleotides between the initial sequence and the reverse complement can be any length including zero."
3. **Palindrome relationship:** "When the intervening length is zero, the composite sequence is a palindromic sequence." Example: `5'---TTACGnnnnnnCGTAA---3'`.

### IUPACpal — Hampson et al. (2021), BMC Bioinformatics (PMC7866733)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7866733/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Formal IR definition (verbatim):** "An *inverted repeat* (IR) is a string that can be expressed in the form WW̄ᴿ for some string *W*." Here W is the left arm and W̄ᴿ is the reverse complement of W (the right arm).
2. **Gapped IR:** the string is a gapped inverted repeat when expressible as **WGW̄ᴿ** for strings W and G where |G| ≥ 0 (G is the gap/loop/spacer).
3. **Imperfect IR:** a gapped inverted repeat within k mismatches when expressible as WGW̄ᴿ with Hamming distance δ_H(W, W̄ᴿ) ≤ k. A **perfect** inverted repeat is ungapped (G = ε) with zero mismatches.
4. **Detection parameters:** users specify minimum arm length, maximum arm length, maximum gap, and maximum mismatches.

### RepeatMasker documentation (Smit, Hubley & Green; Repbase library)

**URL:** https://www.repeatmasker.org/webrepeatmaskerhelp.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (established reference tool / curated database)

**Key Extracted Points:**

1. **Repeat classes (output categories):** SINE, LINE, LTR, DNA (DNA transposons), Satellite, Simple_repeat (microsatellites), Low_complexity, Small RNA, and Unclassified/Unknown.
2. **Classification method (verbatim sense):** "Sequence comparisons in RepeatMasker are performed by the program cross_match … Smith-Waterman-Gotoh"; it "lists all best matches (above a set minimum score) between the query sequence and any of the sequences in the repeat database." Classification is by **homology / best match to a known repeat element in the library**, not exact identity.
3. **Database-driven:** the screened repeat databases are based on Repbase Update (G.I.R.I.); the class of the best-matching library entry is assigned to the query.

---

## Documented Corner Cases and Failure Modes

### From Wikipedia — Tandem repeat

1. **Single copy is not a repeat:** a motif appearing once (copies = 1) is not a tandem repeat; the minimum is two adjacent copies.
2. **Non-primitive units:** a length-2 unit "AA" is really the mononucleotide A repeated; reporting it as a distinct di-nucleotide repeat double-counts. The primitive (shortest) period should be preferred.

### From IUPACpal (Hampson et al. 2021)

1. **Zero gap = palindrome:** when |G| = 0 the IR coincides with a reverse-complement palindrome (even-length).
2. **Odd-length palindrome:** with a 1-nt loop the centre base is unpaired; arms must still be reverse complements.

### From RepeatMasker

1. **No match above threshold:** if no library entry matches above the minimum score, the query is Unclassified/Unknown rather than forced into a class.

---

## Test Datasets

### Dataset: Tandem repeat worked example (Wikipedia definition)

**Source:** Wikipedia "Tandem repeat" (https://en.wikipedia.org/wiki/Tandem_repeat)

| Parameter | Value |
|-----------|-------|
| Sequence | `ATTCGATTCGATTCG` |
| Repeat unit | `ATTCG` (5 bp) |
| Copies | 3 (head-to-tail) |
| Span | start 0, end 15 (exclusive), length 15 |

### Dataset: Inverted repeat / palindrome (Wikipedia + IUPACpal WW̄ᴿ form)

**Source:** Wikipedia "Inverted repeat"; Hampson et al. (2021) PMC7866733

| Parameter | Value |
|-----------|-------|
| Example (gap 0) | `GAATTC` → arm `GAA`, revcomp arm `TTC`; WW̄ᴿ palindrome |
| Example (gap 6) | `TTACGnnnnnnCGTAA` → left arm `TTACG`, gap `nnnnnn`, right arm `CGTAA` = revcomp(`TTACG`) |
| revcomp(`TTACG`) | `CGTAA` (3'→5' complement reversed) |

### Dataset: Repeat classification by motif size (microsatellite / Simple_repeat)

**Source:** Wikipedia "Tandem repeat" (STR 1–6 bp); RepeatMasker classes

| Query motif | Size | RepeatMasker-style class |
|-------------|------|--------------------------|
| `A` | 1 | Simple_repeat (mononucleotide) |
| `CA` | 2 | Simple_repeat (dinucleotide) |
| `CAG` | 3 | Simple_repeat (trinucleotide) |
| 7+ bp tandem | >6 | Satellite / minisatellite (not STR) |

---

## Assumptions

1. **ASSUMPTION: `ClassifyRepeat` library matching is exact-substring containment, not Smith-Waterman homology.** RepeatMasker uses Smith-Waterman-Gotoh against Repbase. A full local-alignment + curated Repbase library is out of scope for one unit; the implemented `ClassifyRepeat(sequence, repeatDb)` assigns the class of the library entry that the query exactly contains / is contained by, and falls back to motif-size simple-repeat classification when no library entry matches. The repeat-class vocabulary (SINE/LINE/LTR/DNA/Satellite/Simple_repeat/Unknown) is source-backed; only the *matching* relaxation (exact substring vs. scored alignment) is assumed. This is documented as a Framework/Simplified limitation, not an invented constant.

---

## Recommendations for Test Coverage

1. **MUST Test:** tandem repeat `ATTCGATTCGATTCG` is found as unit `ATTCG`, 3 copies, span [0,15) — Evidence: Wikipedia Tandem repeat worked example.
2. **MUST Test:** single occurrence of a motif (copies = 1) is NOT reported as a tandem repeat — Evidence: Wikipedia ("two or more").
3. **MUST Test:** inverted repeat `GAATTC` detected with arms `GAA`/`TTC` (revcomp), gap 0 — Evidence: Wikipedia + IUPACpal WW̄ᴿ.
4. **MUST Test:** gapped inverted repeat `TTACGAAAAAACGTAA` detected: left `TTACG`, right `CGTAA`=revcomp, gap 6 — Evidence: Wikipedia example, IUPACpal WGW̄ᴿ.
5. **MUST Test:** `ClassifyRepeat` returns the library class of an exactly-matching entry (e.g. an `AluY` SINE) — Evidence: RepeatMasker best-match classification.
6. **MUST Test:** `ClassifyRepeat` with no library match falls back to Simple_repeat by motif size / Unknown — Evidence: RepeatMasker Unclassified.
7. **SHOULD Test:** primitive-unit preference — `AAAAAA` reported as mononucleotide `A`, not `AA`/`AAA` — Rationale: avoids double counting (Wikipedia non-primitive corner case).
8. **SHOULD Test:** null / empty / too-short input handling — Rationale: documented input-validation contract; matches sibling analyzer methods.
9. **COULD Test:** O(n²) invariant — every reported tandem repeat's reported `sequence` equals `sequence[start..end]` and is an integer number of unit copies — Rationale: structural invariant of the definition.

---

## References

1. Wikipedia contributors. 2026. *Tandem repeat*. Wikipedia. https://en.wikipedia.org/wiki/Tandem_repeat (cites Duitama J et al. 2014, *Nucleic Acids Res* 42(9):5728–5741, https://doi.org/10.1093/nar/gku212; U.S. NLM MeSH "Tandem Repeat Sequences").
2. Wikipedia contributors. 2026. *Inverted repeat*. Wikipedia. https://en.wikipedia.org/wiki/Inverted_repeat (cites Ussery DW, Wassenaar TM, Borini S. 2008. *Computing for Comparative Microbial Genomics*, Springer; Ye C, Ji G, Liang C. 2014. detectIR. *PLoS ONE* 9(11):e113349, https://doi.org/10.1371/journal.pone.0113349).
3. Hampson SE, Pissis SP, et al. 2021. *IUPACpal: efficient identification of inverted repeats in IUPAC-encoded DNA sequences*. BMC Bioinformatics 22:51. https://doi.org/10.1186/s12859-021-03983-2 (PMC7866733: https://pmc.ncbi.nlm.nih.gov/articles/PMC7866733/).
4. Smit AFA, Hubley R, Green P. *RepeatMasker Open-4*. RepeatMasker documentation. https://www.repeatmasker.org/webrepeatmaskerhelp.html (Repbase Update, G.I.R.I.).

---

## Change History

- **2026-06-13**: Initial documentation.
