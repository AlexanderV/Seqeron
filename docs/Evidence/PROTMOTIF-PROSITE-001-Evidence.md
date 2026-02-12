# Evidence Artifact: PROTMOTIF-PROSITE-001

**Test Unit ID:** PROTMOTIF-PROSITE-001
**Algorithm:** PROSITE Pattern Matching
**Date Collected:** 2026-02-12

---

## Online Sources

### PROSITE User Manual — PA Line Specification

**URL:** https://prosite.expasy.org/prosuser.html
**Accessed:** 2026-02-12
**Authority rank:** 2 (Official specification)

**Key Extracted Points:**

1. **Standard IUPAC codes:** Patterns use one-letter amino acid codes (A–Y).
2. **Any amino acid:** The symbol `x` matches any amino acid.
3. **Character class:** `[ALT]` means Ala or Leu or Thr.
4. **Exclusion class:** `{AM}` means any amino acid except Ala and Met.
5. **Separator:** Each element is separated by `-`.
6. **Repetition:** An element followed by `(n)` or `(n,m)` indicates repetition.
   - `x(3)` = `x-x-x`; `x(2,4)` = 2 to 4 of any; `A(3)` = `A-A-A`.
   - **Range `(n,m)` only valid with `x`**, not with specific amino acids.
   - Fixed `(n)` valid with any element.
7. **N-terminus anchor:** `<` at start restricts pattern to N-terminus.
8. **C-terminus anchor:** `>` at end restricts pattern to C-terminus.
9. **C-terminus inside brackets:** In rare cases (PS00267, PS00539), `>` can appear
   inside square brackets: `[G>]` means G or end-of-sequence.
10. **Period:** A period `.` ends the pattern in data files (PA lines).

### ScanProsite Documentation

**URL:** https://prosite.expasy.org/scanprosite/scanprosite_doc.html
**Accessed:** 2026-02-12
**Authority rank:** 2 (Official specification)

**Key Extracted Points:**

1. **Extended syntax:** If pattern contains no ambiguous residues, separators `-`
   can be omitted (e.g., `MASKE` = `M-A-S-K-E`).
2. **Match modes:** greedy (default), overlap, include — three parameters control
   matching behavior.
3. **Repetition note:** `A(2,4)` is NOT a valid pattern element (ranges only for `x`).
4. **Pattern validation examples:** `[AC]-x-V-x(4)-{ED}` and `<A-x-[ST](2)-x(0,1)-V`.

### PROSITE Entry PS00001 (N-glycosylation)

**URL:** https://prosite.expasy.org/PS00001
**Accessed:** 2026-02-12
**Authority rank:** 2 (Official database entry)

**Key Extracted Points:**

1. **Pattern:** `N-{P}-[ST]-{P}.`
2. **Accession:** PS00001
3. **Entry name:** ASN_GLYCOSYLATION
4. **Site:** carbohydrate at position 1

### PROSITE Entry PS00028 (Zinc Finger C2H2)

**URL:** https://prosite.expasy.org/PS00028
**Accessed:** 2026-02-12
**Authority rank:** 2 (Official database entry)

**Key Extracted Points:**

1. **Pattern:** `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H.`
2. **Accession:** PS00028
3. **Entry name:** ZINC_FINGER_C2H2_1
4. **Max repeat:** 36 (up to 36 zinc finger domains per protein)

### ScanProsite Verified Dataset: Human Transferrin (P02787)

**URL:** https://prosite.expasy.org/cgi-bin/prosite/scanprosite/PSScan.cgi?seq=P02787&sig=PS00001&output=json
**Accessed:** 2026-02-12
**Authority rank:** 2 (Official tool output)

**Key Extracted Points:**

1. **Matches found:** 2 N-glycosylation sites (PS00001) in human transferrin.
2. **Match 1:** positions 432–435 (1-based, ScanProsite numbering).
3. **Match 2:** positions 630–633 (1-based, ScanProsite numbering).
4. **UniProt accession:** P02787, TRFE_HUMAN.

### Hulo et al. (2007) — PROSITE Database Publication

**URL:** https://doi.org/10.1093/nar/gkm977
**Accessed:** 2026-02-12
**Authority rank:** 1 (Peer-reviewed paper)

**Key Extracted Points:**

1. PROSITE is the authoritative collection of protein motif patterns.
2. Patterns use a formal notation documented in the PROSITE User Manual.

### De Castro et al. (2006) — ScanProsite

**URL:** https://doi.org/10.1093/nar/gkl124
**Accessed:** 2026-02-12
**Authority rank:** 1 (Peer-reviewed paper)

**Key Extracted Points:**

1. ScanProsite is the reference implementation for scanning PROSITE patterns.
2. Overlapping matches are supported via match mode parameters.

---

## Documented Corner Cases and Failure Modes

### From PROSITE User Manual

1. **Empty pattern:** No specification; implementation should return empty.
2. **Trailing period:** Patterns in data files end with `.` — must be ignored
   during conversion.
3. **C-terminus inside brackets:** `[G>]` in rare patterns (PS00267, PS00539) means
   G-or-end-of-sequence. Implementation converts `[G>]` to `(?:G|$)` regex alternation.
4. **Range only for x:** `A(2,4)` is invalid; only `x(2,4)` is valid. Fixed
   repetition `A(3)` is valid.

### From ScanProsite Documentation

1. **Greedy matching:** Default mode extends variable-length elements maximally.
2. **Overlap handling:** Default allows partially overlapping matches but not
   included (entirely contained) matches.

---

## Test Datasets

### Dataset: PROSITE Syntax Conversion Verification

**Source:** PROSITE User Manual PA line specification, official PROSITE entries

| Input Pattern | Expected Regex | Source |
|---------------|---------------|--------|
| `R-G-D` | `RGD` | PS00016 |
| `A-x-G` | `A.G` | User Manual (x = any) |
| `x(3)-A` | `.{3}A` | User Manual (repetition) |
| `A-x(2,4)-G` | `A.{2,4}G` | User Manual (range) |
| `[ST]-x-[RK]` | `[ST].[RK]` | PS00005 |
| `N-{P}-[ST]-{P}` | `N[^P][ST][^P]` | PS00001 |
| `<M-x-K` | `^M.K` | User Manual (N-terminus) |
| `A-x-G>` | `A.G$` | User Manual (C-terminus) |
| `[RK](2)-x-[ST]` | `[RK]{2}.[ST]` | PS00004 |
| `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` | `C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H` | PS00028 |
| `G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}` | `G[^EDRKHPFYW].{2}[STAGCN][^P]` | PS00008 |
| `D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]` | `D[^W][DNS][^ILVFYW][DENSTG][DNQGHRK][^GP][LIVMC][DENQSTAGC].{2}[DE][LIVMFYW]` | PS00018 |
| (empty) | (empty) | Edge case |
| `R-G-D.` | `RGD` | User Manual (trailing period) |
| `F-[IVFY]-G-[LM]-M-[G>].` | `F[IVFY]G[LM]M(?:G\|$)` | PS00267 (C-term inside brackets) |
| `F-[GSTV]-P-R-L-[G>].` | `F[GSTV]PRL(?:G\|$)` | PS00539 (C-term inside brackets) |

### Dataset: PROSITE Pattern Matching Verification

**Source:** Manual derivation from pattern definitions

| Sequence | Pattern | Expected Matches (0-based start) | Reasoning |
|----------|---------|----------------------------------|-----------|
| `AAARGDAAA` | `R-G-D` | [3] | RGD at index 3 |
| `AANASAAANGTAAA` | `N-{P}-[ST]-{P}` | [2, 8] | NAS at 2 (A≠P,S∈[ST],A≠P); NGT at 8 (G≠P,T∈[ST],A≠P) |
| `AANPSAAA` | `N-{P}-[ST]-{P}` | [] | P at pos 3 violates {P} |
| `AASARKAAA` | `[ST]-x-[RK]` | [2] | SAR: S∈[ST], A=any, R∈[RK] |
| (empty) | `R-G-D` | [] | Empty sequence |
| `AAARGDAAA` | (empty) | [] | Empty pattern |

### Dataset: Human Transferrin N-glycosylation (Real Published Data)

**Source:** ScanProsite scan of P02787 against PS00001 (JSON output)

| Parameter | Value |
|-----------|-------|
| Protein | Human Transferrin (P02787, TRFE_HUMAN) |
| Pattern | PS00001: `N-{P}-[ST]-{P}` |
| Match count | 2 |
| Match 1 position | 432–435 (1-based) / 431–434 (0-based) |
| Match 2 position | 630–633 (1-based) / 631–632 (0-based) |

---

## Recommendations for Test Coverage

1. **MUST Test:** ConvertPrositeToRegex converts simple literal pattern (PS00016 RGD) — Evidence: PROSITE PS00016
2. **MUST Test:** ConvertPrositeToRegex converts `x` to `.` — Evidence: PROSITE User Manual
3. **MUST Test:** ConvertPrositeToRegex converts `x(n)` exact repetition — Evidence: PROSITE User Manual
4. **MUST Test:** ConvertPrositeToRegex converts `x(n,m)` range repetition — Evidence: PROSITE User Manual
5. **MUST Test:** ConvertPrositeToRegex converts character class `[ABC]` — Evidence: PROSITE User Manual
6. **MUST Test:** ConvertPrositeToRegex converts exclusion class `{ABC}` — Evidence: PROSITE User Manual
7. **MUST Test:** ConvertPrositeToRegex handles `<` N-terminus anchor — Evidence: PROSITE User Manual
8. **MUST Test:** ConvertPrositeToRegex handles `>` C-terminus anchor — Evidence: PROSITE User Manual
9. **MUST Test:** ConvertPrositeToRegex handles element repetition `A(n)` — Evidence: PS00004 `[RK](2)-x-[ST]`
10. **MUST Test:** ConvertPrositeToRegex converts complex real pattern (PS00028) — Evidence: PROSITE PS00028
11. **MUST Test:** ConvertPrositeToRegex converts complex PS00008 pattern — Evidence: PROSITE PS00008
12. **MUST Test:** ConvertPrositeToRegex converts complex PS00018 pattern — Evidence: PROSITE PS00018
13. **MUST Test:** ConvertPrositeToRegex handles empty input — Evidence: Edge case
14. **MUST Test:** ConvertPrositeToRegex strips trailing period — Evidence: PROSITE User Manual
15. **MUST Test:** FindMotifByProsite matches PS00001 at correct positions — Evidence: Manual derivation
16. **MUST Test:** FindMotifByProsite matches PS00016 (RGD) at correct position — Evidence: PS00016
17. **MUST Test:** FindMotifByProsite returns empty when no match — Evidence: Edge case
18. **MUST Test:** FindMotifByProsite handles empty sequence — Evidence: Edge case
19. **MUST Test:** FindMotifByProsite handles empty pattern — Evidence: Edge case
20. **MUST Test:** FindMotifByProsite is case-insensitive — Evidence: Implementation spec
21. **MUST Test:** FindMotifByProsite finds multiple matches — Evidence: PS00001 can match multiple times
22. **MUST Test:** ConvertPrositeToRegex handles `[G>]` C-terminus inside brackets (PS00267) — Evidence: PROSITE User Manual §IV.E
23. **MUST Test:** ConvertPrositeToRegex handles `[G>]` C-terminus inside brackets (PS00539) — Evidence: PROSITE User Manual §IV.E
24. **MUST Test:** FindMotifByProsite matches `[G>]` pattern via G branch — Evidence: PS00267
25. **MUST Test:** FindMotifByProsite matches `[G>]` pattern via C-terminus branch — Evidence: PS00267
26. **MUST Test:** FindMotifByProsite rejects `[G>]` mid-sequence without G — Evidence: PS00267
27. **SHOULD Test:** FindMotifByProsite with real protein (Human Transferrin PS00001) — Evidence: ScanProsite output
28. **SHOULD Test:** N-terminal anchored pattern only matches at start — Evidence: PROSITE User Manual
29. **SHOULD Test:** C-terminal anchored pattern only matches at end — Evidence: PROSITE User Manual

---

## References

1. PROSITE User Manual. SIB Swiss Institute of Bioinformatics. https://prosite.expasy.org/prosuser.html (Accessed 2026-02-12)
2. ScanProsite Documentation. https://prosite.expasy.org/scanprosite/scanprosite_doc.html (Accessed 2026-02-12)
3. Hulo N, Bairoch A, Bulliard V, et al. (2007). "The 20 years of PROSITE." Nucleic Acids Res. 36(Database issue):D245-9. https://doi.org/10.1093/nar/gkm977
4. De Castro E, Sigrist CJA, Gattiker A, et al. (2006). "ScanProsite: detection of PROSITE signature matches and ProRule-associated functional and structural residues in proteins." Nucleic Acids Res. 34:W362-365. https://doi.org/10.1093/nar/gkl124
5. PROSITE PS00001 (ASN_GLYCOSYLATION). https://prosite.expasy.org/PS00001 (Accessed 2026-02-12)
6. PROSITE PS00028 (ZINC_FINGER_C2H2_1). https://prosite.expasy.org/PS00028 (Accessed 2026-02-12)

---

## Change History

- **2026-02-12**: Initial documentation.
- **2026-02-13**: Implemented `[G>]` C-terminus inside brackets (PS00267, PS00539). Removed all assumptions.
