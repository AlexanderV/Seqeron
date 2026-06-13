# Evidence Artifact: TRANS-SIXFRAME-001

**Test Unit ID:** TRANS-SIXFRAME-001
**Algorithm:** Six-Frame Translation and Open Reading Frame (ORF) finding
**Date Collected:** 2026-06-13

---

## Online Sources

### EMBOSS transeq — application documentation (reference implementation, frame numbering)

**URL:** https://emboss.sourceforge.net/apps/cvs/emboss/apps/transeq.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation in an established library — EMBOSS)

**Key Extracted Points:**

1. **Reverse-frame definition (phase-locked convention):** "The reverse frame '-1' is defined as the translation you get when you use the reverse-complement of the sequence with the same codon phase as the codon in frame '1'." Frames -2 and -3 use the phase of frames 2 and 3 respectively.
2. **Worked example (phase-locked):** For sequence `ACTGG`, frame -1 translation is `S` (phase of frame 1: reverse-complement read as `...AGT`), frame -2 is `QX`, frame -3 is `PV`.
3. **Documented alternative convention (the one used here):** The doc explicitly records a second, widely used convention: "The alternative way of generating the reverse translation frames used by some people is that frame -1 is made by taking the frame '1' of the reverse complement. There is no correspondence between the codons used in frame 1 and -1 …" — i.e., frame -1 = offset 0 of the reverse complement.
4. **Frame option values:** Permitted `-frame` values are `1, 2, 3, F` (forward three), `-1, -2, -3, R` (reverse three), and `6` (all six frames) — confirming exactly six frames.

### Biopython `Bio.SeqUtils.six_frame_translations` — reference implementation (governing convention)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation — Biopython)

**Key Extracted Points:**

1. **Reverse strand computed once:** `anti = reverse_complement(seq)` (`reverse_complement_rna` for RNA).
2. **Reverse-frame offsets (independent / "alternative" convention):** the reverse frames are filled as `frames[-(i + 1)] = translate(anti[i : i + fragment_length], genetic_code)[::-1]` for `i in range(3)`. Thus frame -1 = translate(reverse-complement, offset 0), frame -2 = offset 1, frame -3 = offset 2. The trailing `[::-1]` only reverses the *string for positional display alignment* under the forward sequence; the residue multiset/content of each reverse-frame protein is the translation of the reverse complement read 5'→3' at that offset.
3. **Forward frames:** `frames[i + 1] = translate(seq[i : i + fragment_length], genetic_code)` for `i in range(3)`, i.e. forward frames are offsets 0, 1, 2 of the input.
4. **Six frames keyed ±1..±3:** the dictionary uses keys `+1,+2,+3` and `-1,-2,-3`.

### NCBI — The Genetic Codes, Standard Code (transl_table=1)

**URL:** https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
**Accessed:** 2026-06-13
**Authority rank:** 2 (official NCBI specification / standard)

**Key Extracted Points:**

1. **Codon table 1 rows (verbatim):**
   `AAs    = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG`
   `Starts = ---M------**--*----M---------------M----------------------------`
   `Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG`
   `Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG`
   `Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`
2. **Start codons (the three `M` positions in the Starts line):** TTG, CTG, ATG.
3. **Stop codons:** TAA, TAG, TGA.
4. **Initiator translated as Met:** "The initiator codon - whether it is AUG, CTG, TTG or something else, - is by default translated as methionine (Met, M)." (Note: this is a display convention for the *initiator* position; see Assumptions for how the repository handles it.)

### EMBOSS getorf — application documentation (ORF definition)

**URL:** https://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation — EMBOSS)

**Key Extracted Points:**

1. **ORF definitions:** getorf recognises an ORF as either "a region that is free of STOP codons" or "a region that begins with a START codon and ends with a STOP codon."
2. **`-find` option values:** 0 = translation between STOP codons; 1 = translation between START and STOP codons; 2/3 = corresponding nucleic sequences; 4/5 = flanking nucleotides. Option **1 (START→STOP)** is the model implemented here.
3. **Minimum size:** default `-minsize` is 30 bases (= 10 amino acids); it measures nucleotide length.
4. **Both strands & reverse numbering:** by default both forward and reverse strands are searched; reverse-strand ORFs are labelled `(REVERSE SENSE)`.

### Wikipedia — Reading frame (uses cited primaries: Lodish 2007; Pierce 2012; Badger & Olsen 1999; Anderson et al. 1981)

**URL:** https://en.wikipedia.org/wiki/Reading_frame
**Accessed:** 2026-06-13
**Authority rank:** 4 (encyclopedic, citing primary sources)

**Key Extracted Points:**

1. **Six frames:** "any given sequence of DNA can therefore be read in six different ways: Three reading frames in one direction (starting at different nucleotides) and three in the opposite direction."
2. **Reverse strand direction:** the three additional frames "may be read from the other, complementary strand in the 5′→3′ direction along this strand"; the 5′→3′ direction on the second strand corresponds to 3′→5′ on the first — i.e. the reverse-complement read forward.

---

## Documented Corner Cases and Failure Modes

### From EMBOSS transeq / Biopython

1. **Incomplete trailing codon:** translation stops at the last complete codon; trailing 1 or 2 nucleotides that cannot form a full codon are ignored (`fragment_length = 3 * ((length - i) // 3)` in Biopython).
2. **Reverse-frame numbering ambiguity:** two conventions exist (phase-locked vs. independent offset). The repository follows the Biopython "independent offset" convention (frame -1 = reverse-complement offset 0). This is the documented "alternative" in EMBOSS transeq.

### From EMBOSS getorf

1. **No START codon present:** under `-find 1` (START→STOP) no ORF is emitted for a region lacking a START codon.
2. **ORF running off the end of the sequence:** a region beginning at a START codon with no downstream STOP before the sequence end is an incomplete ORF.
3. **minsize filtering:** ORFs shorter than the minimum are discarded.

### From NCBI Standard Code

1. **IUPAC-ambiguous codons:** codons containing ambiguity codes (N, R, Y, …) are not in the 64-codon table.

---

## Test Datasets

### Dataset: Six-frame translation of `ACTGG` (EMBOSS transeq worked example)

**Source:** EMBOSS transeq documentation (worked example), cross-checked with Biopython convention.

| Frame | Phase-locked (EMBOSS default) | Independent offset (Biopython / this repo) |
|-------|-------------------------------|--------------------------------------------|
| +1 | T (ACT,GG→T) | T |
| -1 | S | P (reverse-complement `CCAGT` offset 0 → CCA=P) |

Reverse complement of `ACTGG` = `CCAGT`. This repo returns frame -1 = `P`; EMBOSS default returns `S`. The conflict is documented (different conventions).

### Dataset: Six-frame translation of a 39-nt sequence (derived from cited algorithm)

**Source:** Standard genetic code (NCBI table 1) applied per Biopython six-frame algorithm.

Input DNA: `ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG` (length 39).
Reverse complement: `CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT`.

| Frame | Protein |
|-------|---------|
| +1 | `MAIVMGR*KGAR*` |
| +2 | `WPL*WAAERVPD` |
| +3 | `GHCNGPLKGCPI` |
| -1 | `LSGTLSAAHYNGH` |
| -2 | `YRAPFQRPITMA` |
| -3 | `IGHPFSGPLQWP` |

### Dataset: ORF finding (forward strand, START→STOP)

**Source:** EMBOSS getorf `-find 1` model; positions derived from the START→STOP definition.

Input DNA: `GGGATGAAACCCTAAGGG`. ATG begins at index 3; stop `TAA` occupies indices 12–14.

| Field | Value |
|-------|-------|
| StartPosition (0-based, start codon first base) | 3 |
| EndPosition (0-based, stop codon last base, inclusive) | 14 |
| Frame | 1 |
| Protein (start residue included, stop excluded) | `MKP` |
| AminoAcidLength | 3 |
| NucleotideLength (End − Start + 1) | 12 |

---

## Assumptions

1. **ASSUMPTION: Reverse-frame numbering convention.** Two documented conventions exist (EMBOSS phase-locked vs. Biopython independent-offset). The repository follows the **Biopython** convention (frame -k = reverse-complement offset k−1), which is the dominant reference-implementation behaviour and is explicitly listed as an accepted alternative in the EMBOSS transeq documentation. This is a convention choice, not an invented value; both produce correct biology, only the −1/−2/−3 labels differ.
2. **ASSUMPTION: Stop codons rendered as `*`; ambiguous IUPAC codons rendered as `X`.** The `*` for stop is universal (NCBI). Rendering ambiguous codons as `X` (unknown amino acid) follows the IUPAC single-letter "any amino acid" code; it is the established behaviour of the existing `GeneticCode.Translate` and is not exercised as a six-frame-specific MUST.
3. **ASSUMPTION: ORF length filter is in amino acids.** getorf's `-minsize` is in nucleotides; the repository's `FindOrfs(minLength)` parameter counts amino acids (protein length). This is an API-shape choice documented in the contract; behaviour is well-defined for any value.

---

## Recommendations for Test Coverage

1. **MUST Test:** `TranslateSixFrames` returns exactly six frames keyed +1,+2,+3,−1,−2,−3 with the exact protein strings of the 39-nt dataset — Evidence: Biopython six_frame_translations algorithm + NCBI table 1.
2. **MUST Test:** Forward frames +1/+2/+3 equal `Translate` at offsets 0/1/2 — Evidence: Biopython forward-frame loop.
3. **MUST Test:** Reverse frames −1/−2/−3 equal translation of the reverse complement at offsets 0/1/2 — Evidence: Biopython reverse-frame loop.
4. **MUST Test:** Incomplete trailing codon is ignored (length not a multiple of 3) — Evidence: `fragment_length` truncation in Biopython.
5. **MUST Test:** `FindOrfs` returns START→STOP ORF with exact StartPosition/EndPosition/Frame/Protein for the ORF dataset — Evidence: getorf `-find 1` model.
6. **MUST Test:** No START codon ⇒ no ORF; ORF below `minLength` filtered — Evidence: getorf `-find 1` + `-minsize`.
7. **MUST Test:** Null input throws `ArgumentNullException`; empty sequence yields six empty frames / no ORFs — Evidence: implementation contract; truncation rule.
8. **SHOULD Test:** Both-strand ORF search finds reverse-strand ORFs with negative frame labels — Rationale: getorf searches both strands by default.
9. **COULD Test:** Alternative start codon (TTG/CTG) initiates an ORF and is translated by its actual residue — Rationale: NCBI table 1 lists TTG/CTG as starts.

---

## References

1. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite — transeq application documentation. https://emboss.sourceforge.net/apps/cvs/emboss/apps/transeq.html
2. Rice P, Longden I, Bleasby A. 2000. EMBOSS — getorf application documentation. https://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html
3. Cock PJA et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25(11):1422–1423. Source: Bio/SeqUtils/__init__.py (six_frame_translations). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
4. NCBI. The Genetic Codes (transl_table=1, The Standard Code). https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
5. Wikipedia contributors. Reading frame. (cites Lodish 2007; Pierce 2012; Badger & Olsen 1999; Anderson et al. 1981). https://en.wikipedia.org/wiki/Reading_frame

---

## Change History

- **2026-06-13**: Initial documentation.
