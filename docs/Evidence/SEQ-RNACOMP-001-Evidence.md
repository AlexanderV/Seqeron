# Evidence Artifact: SEQ-RNACOMP-001

**Test Unit ID:** SEQ-RNACOMP-001
**Algorithm:** RNA-specific Complement (per-base, IUPAC-complete)
**Date Collected:** 2026-06-13

---

## Online Sources

### Biopython — Bio/Data/IUPACData.py (`ambiguous_rna_complement`)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/Data/IUPACData.py
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation in an established library — Biopython)
**Retrieved via:** WebFetch of the raw file URL above; prompt requested verbatim `ambiguous_rna_complement` and `ambiguous_dna_complement` dicts.

**Key Extracted Points:**

1. **RNA complement table (verbatim):** `ambiguous_rna_complement = {"A":"U", "C":"G", "G":"C", "U":"A", "M":"K", "R":"Y", "W":"W", "S":"S", "Y":"R", "K":"M", "V":"B", "H":"D", "D":"H", "B":"V", "X":"X", "N":"N"}`.
2. **DNA complement table (verbatim, for contrast):** `ambiguous_dna_complement = {"A":"T", "C":"G", "G":"C", "T":"A", "M":"K", "R":"Y", "W":"W", "S":"S", "Y":"R", "K":"M", "V":"B", "H":"D", "D":"H", "B":"V", "X":"X", "N":"N"}`. The RNA table is identical except the base alphabet swaps T→U.
3. **Self-complementary codes:** W→W, S→S, X→X, N→N. Reciprocal pairs: A↔U, C↔G, R↔Y, M↔K, D↔H, B↔V.

### Biopython — Bio/Seq.py (`complement_rna` / `_rna_complement_table`)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py
**Accessed:** 2026-06-13
**Authority rank:** 3 (Biopython reference implementation)
**Retrieved via:** WebFetch of the raw file URL above; prompt requested how `complement_rna` builds its table and how T/U are handled.

**Key Extracted Points:**

1. **T treated as U:** the RNA complement table is built as `ambiguous_rna_complement = dict(IUPACData.ambiguous_rna_complement)` then `ambiguous_rna_complement["T"] = ambiguous_rna_complement["U"]`. Since `ambiguous_rna_complement["U"] == "A"`, this means **T → A** in the RNA complement.
2. **Case handling:** `_maketrans` builds a translation table covering "lower case and upper case sequences", i.e. Biopython preserves the case of the input letter (lowercase in → lowercase out).

### Biopython — Bio.Seq module documentation (`complement_rna`, `reverse_complement_rna`)

**URL:** https://biopython.org/docs/1.79/api/Bio.Seq.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (Biopython official API documentation)
**Retrieved via:** WebFetch of the doc page; prompt requested the docstring worked example for `complement_rna` and T handling.

**Key Extracted Points:**

1. **Worked example (forward complement):** `complement("ACG")` → `'TGC'` (DNA), `complement_rna("ACG")` → `'UGC'`.
2. **T handling, verbatim:** "Any T in the sequence is treated as a U."
3. **Full-alphabet worked example (from `reverse_complement_rna`):** `reverse_complement_rna("ACGTUacgtuXYZxyz")` → `'zrxZRXaacguAACGU'`. Un-reversing this string (computed locally with Python `''.join(reversed(...))`) gives the per-base forward complement of `"ACGTUacgtuXYZxyz"` = `"UGCAAugcaaXRZxrz"` (Biopython preserves case; X passes through; Y→R, Z passes through). Mapping per input char: A→U, C→G, G→C, T→A, U→A, X→X, Y→R, Z→Z.

### IUPAC nucleotide codes table (bioinformatics.org SMS)

**URL:** https://www.bioinformatics.org/sms/iupac.html
**Accessed:** 2026-06-13
**Authority rank:** 5 (curated bioinformatics resource summarizing the IUPAC/NC-IUB standard)
**Retrieved via:** WebFetch of the page; prompt requested the code→bases→complement table.

**Key Extracted Points:**

1. **Code → complement (verbatim table):** A→T, C→G, G→C, T→A, U→A, R→Y, Y→R, S→S, W→W, K→M, M→K, B→V, D→H, H→D, V→B, N→N. Gap (`.`/`-`) → gap.
2. **Bases represented:** R=A|G, Y=C|T, S=G|C, W=A|T, K=G|T, M=A|C, B=C|G|T, D=A|G|T, H=A|C|T, V=A|C|G, N=any. (For RNA, U substitutes for T throughout.)

### NC-IUB 1984 nomenclature (primary source — bibliographic confirmation)

**URL:** https://academic.oup.com/nar/article/13/9/3021/2381659 ; https://pubmed.ncbi.nlm.nih.gov/2582368/
**Accessed:** 2026-06-13
**Authority rank:** 2 (official IUPAC-IUB / NC-IUB standard)
**Retrieved via:** WebSearch for "Cornish-Bowden 1985 Nomenclature for incompletely specified bases"; WebFetch of the NAR article page (only metadata is web-accessible; the full table is behind a PDF wall).

**Key Extracted Points:**

1. **Bibliographic record confirmed online:** Cornish-Bowden A. (1985) "Nomenclature for incompletely specified bases in nucleic acid sequences: recommendations 1984." *Nucleic Acids Research* 13(9):3021–3030. DOI 10.1093/nar/13.9.3021. This is the originating standard for the ambiguity codes; the table content used here is taken from the retrievable sources above (bioinformatics.org SMS and Biopython IUPACData), which encode this standard.

---

## Documented Corner Cases and Failure Modes

### From Biopython (Seq.py / IUPACData.py)

1. **T in an RNA context:** treated as U → complement is A (not preserved as T). Documented verbatim: "Any T in the sequence is treated as a U."
2. **X and unknown letters:** `X→X` is an explicit map; characters not in the table (e.g., Z, gaps, digits) are not translated and pass through unchanged.

### From bioinformatics.org IUPAC table

1. **Gap characters (`.`/`-`):** complement is the gap itself (pass-through), not an error.

---

## Test Datasets

### Dataset: Biopython full-alphabet worked example

**Source:** Biopython 1.79 `Bio.Seq` docs — `reverse_complement_rna("ACGTUacgtuXYZxyz")` → `'zrxZRXaacguAACGU'`; un-reversed forward complement = `"UGCAAugcaaXRZxrz"`.

| Input char | Biopython forward complement | This repo `GetRnaComplementBase` |
|-----------|------------------------------|----------------------------------|
| A | U | U |
| C | G | G |
| G | C | C |
| T | A | A |
| U | A | A |
| a | u | U (uppercased per repo convention) |
| c | g | G |
| g | c | C |
| t | a | A |
| u | a | A |
| X | X | X (pass-through) |
| Y | R | R |
| Z | Z | Z (pass-through) |
| x | x | X (uppercased recognized base; x is not a base → pass-through `x`) |
| y | r | R |
| z | z | Z (pass-through) |

### Dataset: IUPAC ambiguity-code complements (RNA alphabet)

**Source:** bioinformatics.org SMS IUPAC table + Biopython `ambiguous_rna_complement`.

| Code | RNA bases | Complement |
|------|-----------|-----------|
| R | A,G | Y |
| Y | C,U | R |
| S | G,C | S |
| W | A,U | W |
| K | G,U | M |
| M | A,C | K |
| B | C,G,U | V |
| D | A,G,U | H |
| H | A,C,U | D |
| V | A,C,G | B |
| N | any | N |

---

## Assumptions

1. **ASSUMPTION: Case normalization to uppercase for recognized bases.** Biopython preserves the case of recognized letters (e.g., `a → u`), whereas this repository's `GetRnaComplementBase` (and its DNA sibling `GetComplementBase`, SEQ-COMP-001) always returns uppercase for recognized bases. This is an established, intentional repository convention ("DnaSequence/RnaSequence normalize to uppercase", per SEQ-COMP-001 MUST-02 test) and is documented in the method's XML remarks. It is the only behavioral divergence from Biopython and does not affect the complement identity (which base pairs with which) — only letter casing. Tests assert the repo (uppercase) convention; the divergence is recorded in the algorithm doc §5.4.

---

## Recommendations for Test Coverage

1. **MUST Test:** Standard RNA pairing A↔U, C↔G (both directions). — Evidence: bioinformatics.org IUPAC table; Biopython `ambiguous_rna_complement`.
2. **MUST Test:** T → A (thymine treated as uracil in RNA complement). — Evidence: Biopython Seq.py (`ambiguous_rna_complement["T"]=...["U"]`); docs "Any T is treated as a U".
3. **MUST Test:** All IUPAC ambiguity codes (R,Y,S,W,K,M,B,D,H,V,N) → correct RNA complement. — Evidence: Biopython `ambiguous_rna_complement`; bioinformatics.org table.
4. **MUST Test:** Lowercase recognized input → uppercase complement (repo convention). — Evidence: SEQ-COMP-001 MUST-02; this unit's Assumption.
5. **MUST Test:** Full-alphabet worked example string per-base. — Evidence: Biopython `reverse_complement_rna("ACGTUacgtuXYZxyz")`.
6. **SHOULD Test:** Non-IUPAC characters (gap `-`, `.`, digit, Z) pass through unchanged. — Rationale: documented pass-through behavior; gaps map to gaps.
7. **SHOULD Test:** Involution on RNA-emitting alphabet: complement(complement(x)) returns to the same letter for self/reciprocal pairs (within the U-alphabet; T is absorbed into U). — Rationale: complement is an involution on the canonical RNA bases and ambiguity codes.
8. **COULD Test:** Distinction from DNA complement — `GetRnaComplementBase('A')` = 'U' vs `GetComplementBase('A')` = 'T'. — Rationale: confirms RNA-specific behavior is not the DNA path.

---

## References

1. Cornish-Bowden A. (1985). Nomenclature for incompletely specified bases in nucleic acid sequences: recommendations 1984. *Nucleic Acids Research* 13(9):3021–3030. https://doi.org/10.1093/nar/13.9.3021 (PubMed: https://pubmed.ncbi.nlm.nih.gov/2582368/)
2. Biopython contributors. (2026, accessed). `Bio/Data/IUPACData.py` — `ambiguous_rna_complement`, `ambiguous_dna_complement`. https://raw.githubusercontent.com/biopython/biopython/master/Bio/Data/IUPACData.py
3. Biopython contributors. (2026, accessed). `Bio/Seq.py` — `complement_rna`, `_rna_complement_table`. https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py
4. Biopython 1.79. `Bio.Seq` module documentation — `complement_rna`, `reverse_complement_rna` worked examples. https://biopython.org/docs/1.79/api/Bio.Seq.html
5. Bioinformatics.org Sequence Manipulation Suite. IUPAC codes table (code → bases → complement). https://www.bioinformatics.org/sms/iupac.html

---

## Change History

- **2026-06-13**: Initial documentation.
