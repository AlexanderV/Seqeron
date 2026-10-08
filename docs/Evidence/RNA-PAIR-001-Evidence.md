# Evidence Artifact: RNA-PAIR-001

**Test Unit ID:** RNA-PAIR-001
**Algorithm:** RNA Base Pairing (CanPair / GetBasePairType / GetComplement)
**Date Collected:** 2026-06-14

---

## Online Sources

### Crick FHC (1966) — Codon–anticodon pairing: the wobble hypothesis

**URL:** https://en.wikipedia.org/wiki/Wobble_base_pair (Wikipedia article citing the primary, used to obtain the primary citation and the G–U rule); primary: Crick FHC (1966) *J Mol Biol* 19:548–555, DOI 10.1016/S0022-2836(66)80022-0
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary) / 4 (Wikipedia → cited primary)

**Retrieval:** WebSearch query "G:U wobble base pair Crick 1966 codon anticodon wobble hypothesis RNA" then WebFetch of https://en.wikipedia.org/wiki/Wobble_base_pair; a second WebSearch "Crick wobble hypothesis 1966 wobble rules table G pairs U C anticodon U pairs A G" confirmed the pairing rules.

**Key Extracted Points:**

1. **Primary citation (verbatim from fetched page):** Crick FH, "Codon–anticodon pairing: the wobble hypothesis", *Journal of Molecular Biology*, August 1966, vol. 19, pages 548–555, DOI 10.1016/S0022-2836(66)80022-0.
2. **G–U is a wobble pair:** "The four main wobble base pairs are guanine–uracil (G–U), hypoxanthine–uracil (I–U), hypoxanthine–adenine (I–A), and hypoxanthine–cytosine (I–C)." Only G–U involves the four standard RNA bases (A,C,G,U); the inosine (I) pairs are not over the standard alphabet.
3. **Wobble pairing rules (from confirming search of the Crick 1966 rules):** first (wobble) base U recognizes A or G; first base G recognizes U or C. Therefore over the standard alphabet, G pairs with C (Watson-Crick) and U (wobble); U pairs with A (Watson-Crick) and G (wobble).
4. **Distinct from Watson-Crick:** a wobble pair "does not follow Watson–Crick base pair rules", yet "The thermodynamic stability of a wobble base pair is comparable to that of a Watson–Crick base pair" — i.e. G–U is a real, distinct pair type, not a Watson-Crick pair.

### Wikipedia — Base pair (canonical Watson-Crick pairs)

**URL:** https://en.wikipedia.org/wiki/Base_pair
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia citing primaries)

**Retrieval:** WebSearch "Watson-Crick base pairs RNA A-U G-C canonical hydrogen bonds definition" then WebFetch of https://en.wikipedia.org/wiki/Base_pair.

**Key Extracted Points:**

1. **Canonical RNA pairs:** the canonical Watson-Crick pairs in RNA are Adenine–Uracil (A•U) and Guanine–Cytosine (G•C); in RNA thymine is replaced by uracil.
2. **Hydrogen bonds:** A•U has 2 hydrogen bonds; G•C has 3 hydrogen bonds.
3. **Pairing is reciprocal / symmetric:** a base pair is between two complementary bases; A•U is the same pair as U•A, and G•C the same as C•G.

### Wikipedia — Nucleic acid notation (IUPAC codes & complements)

**URL:** https://en.wikipedia.org/wiki/Nucleic_acid_notation
**Accessed:** 2026-06-14
**Authority rank:** 2 (official IUPAC-IUB nomenclature, via Wikipedia) — primary: IUPAC-IUB Commission on Biochemical Nomenclature (1970), *Biochemistry* 9(20):4022–4027

**Retrieval:** WebFetch of https://en.wikipedia.org/wiki/Nucleic_acid_notation.

**Key Extracted Points:**

1. **Base complements (verbatim table):** A↔T, C↔G, G↔C, U↔A; W↔W, S↔S, M↔K, K↔M, R↔Y, Y↔R, B↔V, D↔H, H↔D, V↔B, N↔N.
2. **RNA complement:** identical to the DNA complement table except the complement of A is U (not T) in RNA output; T is treated as U.

### Biopython — Bio.Seq.complement_rna (reference implementation)

**URL:** https://biopython.org/docs/latest/api/Bio.Seq.html
**Accessed:** 2026-06-14
**Authority rank:** 3 (established reference library)

**Retrieval:** WebFetch of https://biopython.org/docs/latest/api/Bio.Seq.html.

**Key Extracted Points:**

1. **complement_rna mapping (from documented example):** `complement_rna(Seq("CGAUT"))` returns `Seq("GCUAA")`, i.e. C→G, G→C, A→U, U→A, and T→A (T treated as U, whose complement is A).
2. This matches the repository helper `SequenceExtensions.GetRnaComplementBase` (A→U, U→A, G→C, C→G, T→A).

### ViennaRNA 2.7 — default pair matrix and nucleotide encoding (reference implementation)

**URL:** https://raw.githubusercontent.com/ViennaRNA/ViennaRNA/master/src/ViennaRNA/model.c and .../src/ViennaRNA/sequences/alphabet.c
**Accessed:** 2026-09-28
**Authority rank:** 3 (established reference implementation; Lorenz et al. 2011, Algorithms Mol Biol 6:26)

**Key Extracted Points:**

1. `BP_ENCODING_DEFAULT` (model.c): over `_ A C G U X K I`, the non-zero pair types are CG=1, GC=2, GU=3, UG=4, AU=5, UA=6 — exactly Watson-Crick + G·U wobble.
2. `vrna_nucleotide_encode` (alphabet.c): `c = toupper(c)`; lookup in `"_ACGUTXKI"`; `if (code > 5) code = 0;` `if (code > 4) code--; /* make T and U equivalent */` — case-insensitive, T ≡ U, every other character → 0 (never pairs).
3. Numerical confirmation (RNA Python 2.7.2, `md = RNA.md()`, encodings from `fold_compound(...).sequence_encoding`): over the alphabet "ACGUTacgutNnXxRYI-." exactly 40 ordered pairs have `md.pair != 0` (e.g. AT:5, GT:3, TA:6, TG:4, at:5); `RNA.fold("GGGATAAAAATATCCC")` = `RNA.fold("GGGAUAAAAAUAUCCC")` = `((((((....))))))`, −6.10 (dangles=0).

---

## Documented Corner Cases and Failure Modes

### From the sources

1. **Order independence (Base pair article):** pairing is symmetric; `CanPair(x,y)` must equal `CanPair(y,x)` and `GetBasePairType(x,y)` must equal `GetBasePairType(y,x)`.
2. **G–U is Wobble, not Watson-Crick (Wobble base pair article):** `GetBasePairType('G','U')` must report a distinct Wobble type, never WatsonCrick.
3. **Non-pairs return false / null:** combinations other than A-U, U-A, G-C, C-G, G-U, U-G do not form pairs (e.g., A-A, A-G, A-C, C-U, G-G, C-C).
4. **DNA T (Biopython, IUPAC, ViennaRNA):** T is treated as U for RNA complement (`GetComplement('T')` = A) and for pairing (`CanPair('A','T')` = true, WatsonCrick; `CanPair('G','T')` = true, Wobble), matching ViennaRNA's encoding. (Revised 2026-09: the earlier "T does not pair" contract had no source and contradicted the reference and the library's own folding engines.)

---

## Test Datasets

### Dataset: Canonical RNA base-pair truth table

**Source:** Crick (1966) DOI 10.1016/S0022-2836(66)80022-0; Wikipedia Base pair / Wobble base pair.

| base1 | base2 | CanPair | GetBasePairType |
|-------|-------|---------|-----------------|
| A | U | true | WatsonCrick |
| U | A | true | WatsonCrick |
| G | C | true | WatsonCrick |
| C | G | true | WatsonCrick |
| G | U | true | Wobble |
| U | G | true | Wobble |
| A | A | false | null |
| A | G | false | null |
| A | C | false | null |
| C | U | false | null |
| G | G | false | null |
| C | C | false | null |

### Dataset: RNA complement (IUPAC / Biopython)

**Source:** IUPAC-IUB (1970) Biochemistry 9(20):4022–4027; Biopython complement_rna.

| base | GetComplement |
|------|---------------|
| A | U |
| U | A |
| G | C |
| C | G |
| T | A |
| N | N |
| R | Y |
| Y | R |

---

## Assumptions

1. **ASSUMPTION: Case-insensitive input** — Sources define pairing over uppercase RNA bases. The implementation upper-cases inputs before lookup; the sources do not address case, but lower/upper case denotes the same nucleotide, so this is a non-correctness-affecting normalization (does not change which nucleotide is meant).

---

## Recommendations for Test Coverage

1. **MUST Test:** All four Watson-Crick pairs (A-U, U-A, G-C, C-G) return CanPair=true and GetBasePairType=WatsonCrick — Evidence: Wikipedia Base pair (canonical A•U, G•C).
2. **MUST Test:** G-U and U-G return CanPair=true and GetBasePairType=Wobble (distinct from WatsonCrick) — Evidence: Crick (1966); Wikipedia Wobble base pair.
3. **MUST Test:** Non-pairing combinations (A-A, A-G, A-C, C-U, G-G, C-C) return CanPair=false and GetBasePairType=null — Evidence: pairing rules (A only with U, C only with G).
4. **MUST Test:** GetComplement maps A→U, U→A, G→C, C→G, T→A, N→N — Evidence: IUPAC-IUB (1970), Biopython complement_rna.
5. **MUST Test:** Symmetry — CanPair(x,y)==CanPair(y,x) and GetBasePairType(x,y)==GetBasePairType(y,x) for all pairs — Evidence: base pairing is reciprocal (Wikipedia Base pair).
6. **SHOULD Test:** Case-insensitivity (lowercase inputs behave identically) — Rationale: normalization contract.
7. **SHOULD Test:** GetComplement('T') = A (T treated as U); CanPair does not pair T (RNA alphabet only) — Rationale: Biopython complement_rna for complement; Crick/WC pairing is over {A,C,G,U}.
8. **COULD Test:** Out-of-alphabet / non-ASCII chars return false/null (no out-of-range exception) — Rationale: robustness, not in source.

---

## References

1. Crick, F.H.C. (1966). Codon–anticodon pairing: the wobble hypothesis. *Journal of Molecular Biology* 19(2):548–555. https://doi.org/10.1016/S0022-2836(66)80022-0
2. Wikipedia. Base pair. https://en.wikipedia.org/wiki/Base_pair (accessed 2026-06-14)
3. Wikipedia. Wobble base pair. https://en.wikipedia.org/wiki/Wobble_base_pair (accessed 2026-06-14)
4. IUPAC-IUB Commission on Biochemical Nomenclature (1970). Abbreviations and symbols for nucleic acids, polynucleotides, and their constituents. *Biochemistry* 9(20):4022–4027. (via https://en.wikipedia.org/wiki/Nucleic_acid_notation, accessed 2026-06-14)
5. Biopython. Bio.Seq.complement_rna. https://biopython.org/docs/latest/api/Bio.Seq.html (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-09-28**: Added ViennaRNA reference (pair matrix + T≡U encoding); T now pairs as U.
