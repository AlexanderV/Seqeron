# Evidence Artifact: RNA-DOTBRACKET-001

**Test Unit ID:** RNA-DOTBRACKET-001
**Algorithm:** Dot-Bracket (extended WUSS) Notation — parsing and validation
**Date Collected:** 2026-06-14

---

## Online Sources

### ViennaRNA Package — RNA Structure Notations (Read the Docs)

**URL:** https://viennarna.readthedocs.io/en/latest/io/rna_structures.html
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation documentation — ViennaRNA Package)

**Key Extracted Points:**

1. **Base pairs and unpaired:** "Base pairs by matching pairs of parenthesis `()` and unpaired nucleotides by dots `.`". The structure must use balanced, nested parentheses where each opening bracket has a corresponding closing bracket in the correct order.
2. **Extended (pseudoknot) brackets:** the extended version "may use additional pairs of brackets, such as `<>`, `{}`, and `[]`, and matching pairs of uppercase/lowercase letters." Different pairs of brackets are NOT required to be nested with each other, which is what allows pseudoknots (crossing helices).
3. **Equivalence of bracket families:** three retrieved equivalent encodings of two crossing helices of size 4 — `<<<<[[[[....>>>>]]]]`, `((((AAAA....))))aaaa`, `AAAA{{{{....aaaa}}}}` — show each bracket family / letter pair is an independent pairing system and that the **uppercase** letter is the 5' (opening) partner and the matching **lowercase** letter is the 3' (closing) partner (e.g. `AAAA...aaaa`).

---

### ViennaRNA Package — Dot-Bracket Notation of Secondary Structures

**URL:** https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/utils/struct/dotbracket.html
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation documentation)

**Key Extracted Points:**

1. **Basic notation:** "matching pairs of parenthesis `()` and unpaired nucleotides by dots `.`"; example `((((....))))`.
2. **Extended notation:** permits additional pairs of brackets `<>`, `{}`, `[]` and uppercase/lowercase letter pairs; "different bracket pairs are not required to be nested" — used for pseudoknots.
3. **Equivalent examples (verbatim):** `<<<<[[[[....>>>>]]]]`, `((((AAAA....))))aaaa`, `AAAA{{{{....aaaa}}}}` are three equivalent representations of the same crossing structure.

---

### ViennaRNA Package — WUSS notation / `vrna_db_from_WUSS()`

**URL:** https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/utils/struct/wuss.html (and the WUSS overview returned by web search of the same documentation set)
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation documentation)

**Key Extracted Points:**

1. **Nested base pairs:** "Nested base pairs are annotated by matching pairs of the symbols `<>`, `()`, `{}`, and `[]`."
2. **Pseudoknots:** "annotation of pseudo-knots using pairs of upper-case/lower-case letters ... matching nested pairs of any alphabet character, such as Aa, Bb, etc."
3. **Infernal interpretation:** "any matching nested pair of `()`, `[]`, `{}` symbols indicates a base pair; the exact choice of symbol has no meaning, so long as the left and right partners match up." `vrna_db_from_WUSS()` "flattens all brackets, and treats pseudo-knots annotated by matching pairs of upper/lowercase letters as unpaired nucleotides" — confirming each family is an independent pairing system.

---

### Rfam Documentation — Glossary (WUSS format)

**URL:** https://docs.rfam.org/en/latest/glossary.html
**Accessed:** 2026-06-14
**Authority rank:** 5 (curated database documentation, EBI/Rfam)

**Key Extracted Points:**

1. **WUSS symbols:** `<>` represent "basepairs in simple stem loops"; `()`, `[]`, `{}` denote "basepairs enclosing multifurcations"; i.e. all four bracket families annotate base pairs.
2. **Unpaired markers:** `-` marks "internal loops and bulges"; `,` "single strand between helices"; `:` "single stranded residues external to any secondary structure"; `.` "insertions relative to the consensus." All non-bracket WUSS symbols denote single-stranded (unpaired) residues.

---

## Documented Corner Cases and Failure Modes

### From ViennaRNA / WUSS documentation

1. **Crossing bracket families (pseudoknots):** `([)]` and the size-4 examples above are valid because `()` and `[]` are matched independently; a single shared stack would mis-pair them. A correct parser must keep one stack per family.
2. **Mismatched partners:** "so long as the left and right partners match up" implies a closing symbol must match an opening of the *same* family. `(]` is therefore NOT well-formed even though one open and one close are present.
3. **Letter direction:** the equivalent example `((((AAAA....))))aaaa` fixes uppercase = opening (5') and the matching lowercase = closing (3').

### From Rfam glossary

4. **Non-bracket WUSS symbols are unpaired:** `-`, `,`, `:`, `.` are all single-stranded; a parser/validator must treat them as unpaired (ignore), not as errors.

---

## Test Datasets

### Dataset: ViennaRNA equivalent crossing structures

**Source:** ViennaRNA RNA Structure Notations (retrieved 2026-06-14)

| Input (dot-bracket) | Meaning | Base pairs (0-based, 5'→3') |
|---------------------|---------|------------------------------|
| `((((....))))` | one hairpin, 4 bp | (0,11),(1,10),(2,9),(3,8) |
| `<<<<[[[[....>>>>]]]]` | two crossing helices, 4 bp each (families nest LIFO) | `<` family: (0,15),(1,14),(2,13),(3,12); `[` family: (4,19),(5,18),(6,17),(7,16) |
| `((((AAAA....))))aaaa` | same structure via letters | `(` family: (0,15),(1,14),(2,13),(3,12); `A` family: (4,19),(5,18),(6,17),(7,16) |

### Dataset: Well-formedness (validation)

**Source:** ViennaRNA (balanced/nested requirement) + WUSS ("partners must match up")

| Input | Valid? | Reason |
|-------|--------|--------|
| `(((...)))` | true | balanced, nested |
| `(([[]]))` | true | nested mixed families, each balanced |
| `([)]` | true | crossing families, each family balanced |
| `(((...)` | false | unclosed `(` |
| `...)` | false | `)` with no open partner |
| `)(` | false | `)` before any `(` |
| `(]` | false | mismatched families: `(` unclosed and `]` unopened |

---

## Assumptions

<!-- Every assumption MUST have bold ASSUMPTION prefix. -->

1. **ASSUMPTION: Best-effort parse of malformed input** — `ParseDotBracket` on a malformed string (e.g. a stray `)`) yields only the pairs it can match and silently drops unmatched closers, rather than throwing. The sources define behavior for *well-formed* notation; the documented contract is that callers test well-formedness with `ValidateDotBracket` first. This affects only malformed input and is recorded as a contract decision, not a numeric/scoring value.
2. **ASSUMPTION: Empty/null string** — empty (and null) notation is treated as a valid, pair-free structure (`ValidateDotBracket("") == true`, `ParseDotBracket("")` empty). ViennaRNA does not define empty input; the empty balanced string is unambiguously balanced under the balanced-bracket definition.

---

## Recommendations for Test Coverage

1. **MUST Test:** Simple nested hairpin `((((....))))` parses to the 4 expected (i,j) pairs with exact positions and outermost-with-outermost nesting. — Evidence: ViennaRNA basic notation + example `((((....))))`.
2. **MUST Test:** Crossing families parse independently — `([)]` → (0,2) and (1,3); `<<<<[[[[....>>>>]]]]` and `((((AAAA....))))aaaa` produce the two retrieved 4-bp helices. — Evidence: ViennaRNA equivalent examples; WUSS independent families.
3. **MUST Test:** Validation accepts balanced/nested and crossing-family strings and rejects unclosed, unopened, reversed, and mismatched-family strings (`(]` false). — Evidence: ViennaRNA balanced requirement + WUSS "partners must match up".
4. **SHOULD Test:** Uppercase/lowercase letter pairs parse with uppercase as 5' opener. — Rationale: retrieved example `AAAA...aaaa`.
5. **SHOULD Test:** Non-bracket WUSS symbols (`-`, `,`, `:`) and dots are treated as unpaired. — Rationale: Rfam glossary single-stranded symbols.
6. **COULD Test:** Empty/null handling. — Rationale: documented contract decision (Assumptions).

---

## References

1. Lorenz R, Bernhart SH, Höner zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL (2011). ViennaRNA Package 2.0. Algorithms for Molecular Biology 6:26. RNA Structure Notations documentation: https://viennarna.readthedocs.io/en/latest/io/rna_structures.html (accessed 2026-06-14)
2. ViennaRNA Package. Dot-Bracket Notation of Secondary Structures. https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/utils/struct/dotbracket.html (accessed 2026-06-14)
3. ViennaRNA Package. Washington University Secondary Structure (WUSS) notation. https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/utils/struct/wuss.html (accessed 2026-06-14)
4. Nawrocki EP, Eddy SR (2013). Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics 29(22):2933-2935. WUSS notation as documented in ViennaRNA/Infernal (matching nested pairs of (), [], {} indicate base pairs). https://doi.org/10.1093/bioinformatics/btt509
5. Rfam Documentation — Glossary (WUSS format). https://docs.rfam.org/en/latest/glossary.html (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
