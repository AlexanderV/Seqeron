# Evidence Artifact: MOTIF-CONS-001

**Test Unit ID:** MOTIF-CONS-001
**Algorithm:** Consensus Sequence from a Multiple Alignment (plurality / most-frequent residue)
**Date Collected:** 2026-06-13

---

## Online Sources

### Wikipedia — Consensus sequence (cited primaries used)

**URL:** https://en.wikipedia.org/wiki/Consensus_sequence
**Retrieved:** WebFetch of the URL above on 2026-06-13.
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primary sources Schneider & Stephens 1990; Pierce 2002)

**Key Extracted Points:**

1. **Formal definition:** A consensus sequence is "the calculated sequence of most frequent residues, either nucleotide or amino acid, found at each position in a sequence alignment." Each position's residue is determined by frequency analysis across the aligned sequences — the most frequent base/amino acid at each position becomes the consensus character.
2. **Limitation (not tie-breaking):** A consensus "reduce[s] variability to a single residue per position"; the article does not specify a tie-breaking procedure, deferring richer representations to sequence logos (Schneider & Stephens 1990).

### Rosalind — "Consensus and Profile" (CONS) problem

**URL:** https://rosalind.info/problems/cons/
**Retrieved:** WebFetch of the URL above on 2026-06-13.
**Accessed:** 2026-06-13
**Authority rank:** 5 (curated bioinformatics learning database with verified datasets)

**Key Extracted Points:**

1. **Profile matrix:** a 4×n matrix P where P[1,j] is the number of times 'A' occurs in the jth position of the aligned strings (and likewise C, G, T).
2. **Consensus rule:** "A consensus string is formed by taking the most common symbol at each position; the jth symbol of c corresponds to the symbol having the maximum value in the j-th column of the profile matrix."
3. **Tie behaviour:** "there may be more than one most common symbol, leading to multiple possible consensus strings" — any valid consensus may be returned. (Determinism is therefore implementation-defined; see Assumptions.)
4. **Worked example — sample input (7 DNA strings of length 8):** `ATCCAGCT`, `GGGCAACT`, `ATGGATCT`, `AAGCAACC`, `TTGGAACT`, `ATGCCATT`, `ATGGCACT`.
5. **Worked example — sample output consensus:** `ATGCAACT`. Profile: A = `5 1 0 0 5 5 0 0`; C = `0 0 1 4 2 0 6 1`; G = `1 1 6 3 0 1 0 0`; T = `1 5 0 0 0 1 1 6`.

### EMBOSS `cons` — consensus from a multiple alignment

**URL:** https://www.bioinformatics.nl/cgi-bin/emboss/help/cons
**Retrieved:** WebFetch of the URL above on 2026-06-13.
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation; EMBOSS suite)

**Key Extracted Points:**

1. **Purpose:** "cons calculates a consensus sequence from a multiple sequence alignment."
2. **Selection rule:** the highest-scoring residue in a column becomes the consensus residue if its weighted positive-match count exceeds the "plurality" value, otherwise there is no consensus there.
3. **Plurality default:** "Half the total sequence weighting" — i.e. a residue must be supported by more than half the (equally weighted) sequences to call a confident consensus.
4. **No-consensus output:** at positions lacking consensus, "an 'n' (nucleotide sequence alignment) or an 'x' (protein sequence alignment) character is written to the consensus sequence."

### Geneious / LANL HIV database — tie-breaking conventions

**URL:** https://hfv.lanl.gov/content/sequence/CONSENSUS/AdvConExplain.html
**Retrieved:** WebFetch of the URL above on 2026-06-13.
**Accessed:** 2026-06-13
**Authority rank:** 5 (curated database; Los Alamos HIV sequence database documentation)

**Key Extracted Points:**

1. **Most-frequent rule:** the consensus at each column is the most frequently occurring character in that column.
2. **Tie-breaking options:** ties may be broken (a) with the correct IUPAC ambiguity code (nucleotides only), (b) by a specified residue order, or (c) by an ambiguity symbol ('?'). The web-search summary of the same family of tools (Geneious manual) also documents an explicit **alphabetical** tie-break: "In the event of a tie, the residue letter occurring earlier in the alphabet was chosen."

---

## Documented Corner Cases and Failure Modes

### From Rosalind (CONS)

1. **Equal-length precondition:** the profile/consensus is defined over a "collection of equal-length DNA strings." Inputs of unequal length are outside the stated definition.
2. **Multiple consensus strings:** when two symbols tie for maximum count in a column, more than one consensus string is valid; a deterministic implementation must pick one rule (see Assumptions).

### From EMBOSS `cons`

1. **No-consensus position:** when no residue reaches the plurality threshold, EMBOSS writes 'n' (nucleotide). For a pure most-frequent consensus (no threshold) there is always a most-frequent symbol, so this only arises under a plurality cutoff.

---

## Test Datasets

### Dataset: Rosalind CONS sample

**Source:** Rosalind, "Consensus and Profile" (https://rosalind.info/problems/cons/), accessed 2026-06-13.

| Parameter | Value |
|-----------|-------|
| Input strings | ATCCAGCT, GGGCAACT, ATGGATCT, AAGCAACC, TTGGAACT, ATGCCATT, ATGGCACT |
| String count | 7 |
| String length | 8 |
| Profile A | 5 1 0 0 5 5 0 0 |
| Profile C | 0 0 1 4 2 0 6 1 |
| Profile G | 1 1 6 3 0 1 0 0 |
| Profile T | 1 5 0 0 0 1 1 6 |
| Consensus | ATGCAACT |

### Dataset: Alphabetical tie-break (derived)

**Source:** Geneious/LANL alphabetical tie-break rule (above); column composition chosen so two symbols tie for the maximum count.

| Parameter | Value |
|-----------|-------|
| Input strings | AT, GT (column 1: one A, one G — tie; column 2: T,T) |
| Expected consensus | AT |
| Reasoning | Column 1 has A and G tied at count 1; alphabetical tie-break selects A (earlier letter). Column 2 is unanimous T. |

---

## Assumptions

1. **ASSUMPTION: Alphabetical tie-break (A<C<G<T).** Rosalind explicitly permits any most-common symbol on a tie; EMBOSS uses scoring/plurality; the Geneious/LANL family documents an explicit alphabetical tie-break. To make the method deterministic (a library requirement) we adopt the alphabetical-order tie-break documented by Geneious/LANL. This is correctness-affecting only on tied columns; on the Rosalind worked example there are no ties affecting the published consensus, so conformance to the rank-5 dataset is unaffected.
2. **ASSUMPTION: Pure most-frequent consensus, no plurality threshold.** The Registry canonical signature `CreateConsensusFromAlignment(alignedSequences)` takes no threshold parameter, matching the Rosalind/Wikipedia "most common symbol" definition rather than EMBOSS's parameterised plurality. Threshold-based no-consensus ('n'/'x') output is therefore out of scope for this method (the area already exposes IUPAC-degenerate consensus via `GenerateConsensus`).

---

## Recommendations for Test Coverage

1. **MUST Test:** Rosalind CONS sample → consensus `ATGCAACT`. — Evidence: Rosalind CONS sample input/output.
2. **MUST Test:** identical sequences → that exact sequence (every column unanimous). — Evidence: most-frequent rule (Wikipedia/Rosalind).
3. **MUST Test:** alphabetical tie-break: column with A,G tied → A. — Evidence: Geneious/LANL alphabetical tie-break.
4. **MUST Test:** single sequence → returns it unchanged (each column's only symbol is the most frequent). — Evidence: most-frequent rule.
5. **MUST Test:** case-insensitivity (lowercase input normalised) — Rationale: consistent with sibling MotifFinder methods (`ToUpperInvariant`).
6. **SHOULD Test:** null input → ArgumentNullException; empty collection → empty string. — Rationale: library failure-mode convention.
7. **SHOULD Test:** unequal-length sequences → ArgumentException. — Rationale: Rosalind equal-length precondition.
8. **COULD Test:** non-ACGT character → ArgumentException. — Rationale: alphabet validation consistent with `CreatePwm`.

---

## References

1. Wikipedia contributors. 2026. Consensus sequence. Wikipedia. https://en.wikipedia.org/wiki/Consensus_sequence (citing Schneider TD, Stephens RM. 1990. Sequence Logos. Nucleic Acids Res 18(20):6097–6100, https://doi.org/10.1093/nar/18.20.6097).
2. Rosalind. Consensus and Profile (CONS). https://rosalind.info/problems/cons/ (accessed 2026-06-13).
3. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite. Trends in Genetics 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2. Program documentation: https://www.bioinformatics.nl/cgi-bin/emboss/help/cons (accessed 2026-06-13).
4. Los Alamos HIV Sequence Database. Advanced Consensus Maker — explanation. https://hfv.lanl.gov/content/sequence/CONSENSUS/AdvConExplain.html (accessed 2026-06-13).

---

## Change History

- **2026-06-13**: Initial documentation.
