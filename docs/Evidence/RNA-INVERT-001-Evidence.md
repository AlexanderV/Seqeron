# Evidence Artifact: RNA-INVERT-001

**Test Unit ID:** RNA-INVERT-001
**Algorithm:** RNA Inverted Repeats (potential stem regions)
**Date Collected:** 2026-06-14

---

## Online Sources

### IUPACpal: efficient identification of inverted repeats in IUPAC-encoded DNA sequences (Alamro, Alzamel, Iliopoulos, Pissis, Watts; BMC Bioinformatics 2021)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7866733/
**Accessed:** 2026-06-14 (fetched the PMC full-text page; query "inverted repeat definition reverse complement EMBOSS einverted stem loop")
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Core definition:** "An inverted repeat (IR) is a single stranded sequence of nucleotides with a subsequent downstream sequence consisting of its reverse complement."
2. **Formal model (gapped IR):** "a string that can be expressed in the form WGW̄ᴿ for some pair of strings W and G where |G|≥0" — W is the left arm, G the gap/spacer/loop (length ≥ 0), W̄ᴿ the right arm equal to the reverse complement of W.
3. **Mismatch model:** "a gapped inverted repeat within k mismatches when it can be expressed in the form WGW̄ᴿ with δH(W,W̄ᴿ)≤k", where δH is the Hamming distance (count of differing positions). A perfect IR is the case k = 0.
4. **Watson-Crick basis:** "the natural choice of a complement function over the alphabet Σ={A,C,G,T} is such that A⟷T and C⟷G" (for RNA, A⟷U and C⟷G).
5. **Tool parameters:** IUPACpal accepts minimum length, maximum length, maximum gap size, and maximum mismatches allowed.
6. **Worked example:** the sequence `CT-CGCAGTCACCG-GA` is an IR with left arm W = `CT`, an internal gap, and right arm = reverse complement of the left arm, with a single mismatch toward the tail ends.

### Wikipedia: Inverted repeat (definition cites Ussery et al. 2008)

**URL:** https://en.wikipedia.org/wiki/Inverted_repeat
**Accessed:** 2026-06-14 (fetched the article; primary citation is Ussery, Wassenaar & Borini 2008)
**Authority rank:** 4 (Wikipedia citing a primary textbook source — Ussery et al. 2008)

**Key Extracted Points:**

1. **Definition:** an inverted repeat is "a single stranded sequence of nucleotides followed downstream by its reverse complement"; the intervening sequence may be any length including zero.
2. **Worked example:** `5'---TTACGnnnnnnCGTAA---3'` is an inverted repeat, where `nnnnnn` is the (any) intervening loop. The right arm `CGTAA` is the reverse complement of the left arm `TTACG` (TTACG → complement AATGC → reverse CGTAA).
3. **Palindrome relationship:** with zero intervening nucleotides the composite `5'TTACGCGTAA3'` is a palindromic sequence (equal to its own reverse complement).

### EMBOSS einverted manual (EMBOSS application documentation)

**URL:** https://emboss.bioinformatics.nl/cgi-bin/emboss/help/einverted
**Accessed:** 2026-06-14 (fetched the application help page)
**Authority rank:** 3 (reference implementation / established tool, EMBOSS suite)

**Key Extracted Points:**

1. **Purpose:** "einverted finds inverted repeats (stem loops) in nucleotide sequences. It identifies regions of local alignment of the input sequence and its reverse complement that exceed a threshold score." This confirms inverted repeats are equivalent to potential stem-loops, found by comparing a sequence to its reverse complement.
2. **Stem-loop = inverted repeat:** the two arms are complementary; the intervening region is the loop/bulge. Mismatches and gaps correspond to bulges in the stem-loop.
3. **Scoring (default, EMBOSS):** match +3, mismatch −4, gap penalty 12, minimum reportable score 50. (Not used here — the repository implementation reports perfect, ungapped arms only; see Assumptions.)

---

## Documented Corner Cases and Failure Modes

### From IUPACpal (PMC7866733)

1. **Zero gap (palindrome):** |G| ≥ 0; when |G| = 0 the IR is a palindrome (`WW̄ᴿ`). This repository's RNA variant requires a loop (minSpacing ≥ 1 by default 3), so true palindromes with no loop are out of scope of the default call.
2. **Minimum/maximum arm length:** arms shorter than the minimum length are not reported; arms are extended to their maximal complementary run.

### From Wikipedia / Ussery et al. (2008)

1. **Intervening length:** the loop may be any length including zero; for RNA hairpins a non-zero loop is biologically required (a hairpin loop typically needs ≥ 3 unpaired nucleotides).

---

## Test Datasets

### Dataset: Wikipedia inverted-repeat example (perfect arms)

**Source:** Wikipedia "Inverted repeat" (citing Ussery, Wassenaar & Borini 2008), example `5'---TTACGnnnnnnCGTAA---3'`.

Adapted to RNA with a concrete 6-nt loop (`AAAAAA`): `UUACGAAAAAACGUAA` (length 16).

| Parameter | Value |
|-----------|-------|
| Sequence | `UUACGAAAAAACGUAA` |
| Left arm W | `UUACG` (positions 0–4) |
| Loop G | `AAAAAA` (positions 5–10, length 6) |
| Right arm W̄ᴿ | `CGUAA` (positions 11–15) |
| Arm length | 5 |
| Right arm = revcomp(left arm)? | Yes: UUACG → complement AAUGC → reverse CGUAA |

### Dataset: Palindromic IR with minimal RNA loop

**Source:** Derived directly from the WGW̄ᴿ definition (IUPACpal) with a 3-nt loop.

| Parameter | Value |
|-----------|-------|
| Sequence | `GGCCAAAGGCC` (length 11) |
| Left arm W | `GGCC` (positions 0–3) |
| Loop G | `AAA` (positions 4–6, length 3) |
| Right arm W̄ᴿ | `GGCC` (positions 7–10) |
| Right arm = revcomp(left arm)? | Yes: GGCC → complement CCGG → reverse GGCC |

---

## Assumptions

1. **ASSUMPTION: Perfect (zero-mismatch), ungapped arms only.** The repository's `GetComplement`/`GetRnaComplementBase` is strict Watson-Crick + IUPAC (no internal gaps in the arm, no mismatches). The implementation therefore reports the perfect-IR case k = 0 of the IUPACpal model (WGW̄ᴿ with δH = 0). EMBOSS einverted additionally scores mismatches/gaps; that scored/DP variant is out of scope and noted as Not Implemented. The k = 0 case is itself a documented special case (Hamming distance 0) in the IUPACpal model, so this is a scope restriction, not an invented behavior.
2. **ASSUMPTION: Loop bounds via minSpacing/maxSpacing.** The public parameters `minSpacing`/`maxSpacing` bound |G| (loop length). This matches IUPACpal's "maximum gap size" and the einverted loop concept; defaults (minSpacing 3, maxSpacing 100) are repository API conventions (3-nt minimal hairpin loop) and do not affect the formal correctness of a reported IR.
3. **ASSUMPTION: Maximal arm + non-overlapping greedy reporting.** For each accepted left start the longest perfect antiparallel arm is reported, and reported repeats do not overlap. This mirrors einverted's maximal-local-alignment reporting and avoids reporting every sub-arm; it does not change which positions are valid complementary pairs.

---

## Recommendations for Test Coverage

1. **MUST Test:** A known IR (`UUACGAAAAAACGUAA`) returns arms at the exact positions with the right arm = reverse complement of the left arm (antiparallel pairing). — Evidence: Wikipedia/Ussery example; IUPACpal WGW̄ᴿ.
2. **MUST Test:** A palindromic IR with a 3-nt loop (`GGCCAAAGGCC`) returns left 0–3, right 7–10, length 4. — Evidence: IUPACpal WGW̄ᴿ definition.
3. **MUST Test:** Antiparallel (not parallel) pairing — a sequence whose arms match in parallel but not antiparallel must NOT be reported. — Evidence: IUPACpal "reverse complement" (reverse, not just complement).
4. **MUST Test:** No inverted repeat present → empty result. — Evidence: definition (no WGW̄ᴿ decomposition exists).
5. **SHOULD Test:** Arm shorter than minLength is not reported; arm exactly minLength is. — Rationale: minimum-length parameter boundary (IUPACpal).
6. **SHOULD Test:** Loop length below minSpacing / above maxSpacing excludes the IR. — Rationale: gap-size bounds (IUPACpal).
7. **COULD Test:** Maximal-arm extension — a perfect arm longer than minLength is reported at its full length, not truncated. — Rationale: einverted maximal local alignment.

---

## References

1. Alamro H, Alzamel M, Iliopoulos CS, Pissis SP, Watts S. 2021. IUPACpal: efficient identification of inverted repeats in IUPAC-encoded DNA sequences. BMC Bioinformatics 22:51. https://doi.org/10.1186/s12859-021-03983-2 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC7866733/)
2. Ussery DW, Wassenaar TM, Borini S. 2008. Computing for Comparative Microbial Genomics: Bioinformatics for Microbiologists. Springer. https://doi.org/10.1007/978-1-84800-255-5 (definition accessed via Wikipedia "Inverted repeat": https://en.wikipedia.org/wiki/Inverted_repeat)
3. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite (einverted application). Trends in Genetics 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2 (einverted manual: https://emboss.bioinformatics.nl/cgi-bin/emboss/help/einverted)

---

## Change History

- **2026-06-14**: Initial documentation.
