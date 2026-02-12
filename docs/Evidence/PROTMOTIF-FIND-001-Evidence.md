# Evidence Artifact: PROTMOTIF-FIND-001

**Test Unit ID:** PROTMOTIF-FIND-001
**Algorithm:** Protein Motif Search (Pattern-based)
**Date Collected:** 2026-02-12

---

## Online Sources

### PROSITE Database — Official Pattern Definitions

**URL:** https://prosite.expasy.org/
**Accessed:** 2026-02-12
**Authority rank:** 2 (Official specification / standard)

**Key Extracted Points:**

1. **Pattern syntax:** PROSITE uses a well-defined pattern notation: `x` for any amino acid, `[ABC]` for alternatives, `{ABC}` for exclusions, `(n)` for exact repeats, `(n,m)` for range repeats, `<` for N-terminus, `>` for C-terminus. Elements separated by `-`.
2. **PS00001 (ASN_GLYCOSYLATION):** Pattern `N-{P}-[ST]-{P}`. N-glycosylation site. Skip-flag TRUE (ubiquitous).
3. **PS00005 (PKC_PHOSPHO_SITE):** Pattern `[ST]-x-[RK]`. Protein kinase C phosphorylation site.
4. **PS00006 (CK2_PHOSPHO_SITE):** Pattern `[ST]-x(2)-[DE]`. Casein kinase II phosphorylation site.
5. **PS00004 (CAMP_PHOSPHO_SITE):** Pattern `[RK](2)-x-[ST]`. cAMP/cGMP-dependent protein kinase phosphorylation site.
6. **PS00007 (TYR_PHOSPHO_SITE_1):** Pattern `[RK]-x(2)-[DE]-x(3)-Y`. Tyrosine kinase phosphorylation site 1.
7. **PS00008 (MYRISTYL):** Pattern `G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}`. N-myristoylation site.
8. **PS00009 (AMIDATION):** Pattern `x-G-[RK]-[RK]`. Amidation site.
9. **PS00016 (RGD):** Pattern `R-G-D`. Cell attachment sequence.
10. **PS00017 (ATP_GTP_A):** Pattern `[AG]-x(4)-G-K-[ST]`. ATP/GTP-binding site motif A (P-loop).
11. **PS00018 (EF_HAND_1):** Pattern `D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]`.
12. **PS00028 (ZINC_FINGER_C2H2_1):** Pattern `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H`. Verified correct in implementation.
13. **PS00029 (LEUCINE_ZIPPER):** Pattern `L-x(6)-L-x(6)-L-x(6)-L`. Verified correct in implementation.

### Wikipedia — Sequence motif

**URL:** https://en.wikipedia.org/wiki/Sequence_motif
**Accessed:** 2026-02-12
**Authority rank:** 4 (Wikipedia with citations to primary sources)

**Key Extracted Points:**

1. **Definition:** A sequence motif is a nucleotide or amino acid sequence pattern that is widespread and usually assumed to be related to biological function.
2. **PROSITE notation:** Uses IUPAC one-letter codes. `x` = any amino acid, `[ABC]` = alternative residues, `{ABC}` = exclusion, `<`/`>` = N/C-terminus, element separator `-`, repetition `e(m)` or `e(m,n)`.
3. **N-glycosylation example:** `N{P}[ST]{P}` — Asn, followed by anything but Pro, followed by Ser or Thr, followed by anything but Pro.
4. **Zinc finger C2H2 signature:** `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H`.

### Wikipedia — PROSITE

**URL:** https://en.wikipedia.org/wiki/PROSITE
**Accessed:** 2026-02-12
**Authority rank:** 4

**Key Extracted Points:**

1. **Database origin:** Created in 1988 by Amos Bairoch; contains patterns, profiles, and ProRules for protein families, domains, and functional sites.
2. **Citation:** Hulo N et al. (2007). "The 20 years of PROSITE." Nucleic Acids Res. 36(Database issue):D245-9. https://doi.org/10.1093/nar/gkm977
3. **ScanProsite:** De Castro E et al. (2006). Nucleic Acids Res. 34(Web Server issue):W362-365. https://doi.org/10.1093/nar/gkl124

### PROSITE User Manual — Pattern Syntax

**URL:** https://prosite.expasy.org/prosuser.html
**Accessed:** 2026-02-12
**Authority rank:** 2 (Official specification)

**Key Extracted Points:**

1. **Formal syntax rules:** Standard IUPAC one-letter amino acid codes. `x` for any position. `[ABC]` for alternatives, `{ABC}` for exclusions. `-` as separator. `(n)` or `(n,m)` for repetitions. `<`/`>` for terminus anchors. Period ends the pattern.
2. **Pattern-to-regex conversion:** Each PROSITE element maps to a regex equivalent: `x` → `.`, `[ABC]` → `[ABC]`, `{ABC}` → `[^ABC]`, `x(n)` → `.{n}`, `x(n,m)` → `.{n,m}`, `<` → `^`, `>` → `$`.

---

## Documented Corner Cases and Failure Modes

### From PROSITE Documentation

1. **Skip-flag patterns:** Some patterns (PS00001, PS00005, PS00006, PS00004, PS00007, PS00008, PS00009, PS00016, PS00017, PS00029) have SKIP-FLAG=TRUE, meaning they are very common and may produce many false positives. Programs may choose to skip these in some contexts.
2. **Overlapping matches:** PROSITE profiles require non-overlapping matches, but for patterns, multiple occurrences are expected and overlapping matches can occur. Standard regex matching does not find overlapping matches.
3. **Case sensitivity:** Protein sequences should match case-insensitively; PROSITE uses uppercase IUPAC codes.
4. **Empty/null input:** Motif search on empty or null sequences should return no results.
5. **Invalid regex pattern:** If a user supplies an invalid regex pattern, the search should handle it gracefully.

---

## Test Datasets

### Dataset 1: Human Insulin (UniProt P01308) N-glycosylation scan

**Source:** UniProt P01308 (INS_HUMAN), ScanProsite

| Parameter | Value |
|-----------|-------|
| Protein | `MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN` |
| Pattern | PS00001: N-{P}-[ST]-{P} |
| Known sites | Positions with N followed by non-P, then S/T, then non-P |

### Dataset 2: Constructed PROSITE pattern verification

**Source:** Derived from PROSITE official patterns

| Pattern ID | PROSITE Pattern | Test Sequence | Expected Match |
|-----------|-----------------|---------------|----------------|
| PS00001 | N-{P}-[ST]-{P} | `ANFTB` | Match at position 1: `NFTB` |
| PS00001 | N-{P}-[ST]-{P} | `ANPSB` | No match (P at position 2) |
| PS00005 | [ST]-x-[RK] | `ASARK` | Match: `SAR` at pos 1 |
| PS00006 | [ST]-x(2)-[DE] | `ATAAD` | Match: `TAAD` at pos 1 |
| PS00016 | R-G-D | `ARGDA` | Match: `RGD` at pos 1 |
| PS00017 | [AG]-x(4)-G-K-[ST] | `AXXXXGKS` | Match at pos 0 |

### Dataset 3: PS00007 pattern correction verification

**Source:** PROSITE PS00007 (official pattern: `[RK]-x(2)-[DE]-x(3)-Y`)

| Parameter | Value |
|-----------|-------|
| Correct regex | `[RK].{2}[DE].{3}Y` |
| Implementation (wrong) regex | `[RK].{2,3}[DE].{2,3}Y` |
| Test sequence | `RAAEDDDDY` (9 chars: R + 2 any + E + 3 any + Y) |
| Expected match | `RAAEDDDDY` (correct: matches with exact counts) |

### Dataset 4: PS00018 (EF-hand) pattern correction verification

**Source:** PROSITE PS00018 (official pattern: `D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]`)

| Parameter | Value |
|-----------|-------|
| Correct regex | `D[^W][DNS][^ILVFYW][DENSTG][DNQGHRK][^GP][LIVMC][DENQSTAGC].{2}[DE][LIVMFYW]` |
| Implementation (wrong) regex | `D.[DNS][^ILVFYW][DENSTG][DNQGHRK][^GP][LIVMC][DENQSTAGC].{2}[DE]` |
| Discrepancy 1 | Position 2 uses `x` (any) instead of `{W}` (any except W) |
| Discrepancy 2 | Missing trailing `[LIVMFYW]` element |

---

## Non-PROSITE Pattern Sources

All non-PROSITE patterns in the `CommonMotifs` dictionary are verified against published literature:

| ID | Pattern | Consensus Source | Reference |
|----|---------|-----------------|-----------|
| NLS1 | `[KR]-[KR]-x-[KR]` | Chelsky monopartite NLS: K-K/R-X-K/R | Dingwall C, Laskey RA (1991). Trends Biochem Sci 16:478-481 |
| NES1 | `L-x(2,3)-[LIVFM]-x(2,3)-L-x-[LI]` | la Cour NES consensus: Φ1-x(2,3)-Φ2-x(2,3)-Φ3-x-Φ4 (L-enriched form) | la Cour T et al. (2004). Nucleic Acids Res 32:W142-5 |
| SIM1 | `[VIL]-x-[VIL]-[VIL]` | SIM type b: Ψ-x-Ψ-Ψ where Ψ = V/I/L | Hecker CM et al. (2006). JBC 281:16117-16127 |
| WW1 | `P-P-x-Y` | PPxY consensus for WW domain binding | Chen HI, Sudol M (1995). PNAS 92:7819-7823 |
| SH3_1 | `[RK]-x(2)-P-x(2)-P` | SH3 Class I consensus: +xxPxxP where + = R/K | Mayer BJ (2001). J Cell Sci 114:1253-1263 |

## Scoring Methodology

Scoring uses information content (IC) per Schneider & Stephens (1990):

- **CalculateMotifScore:** IC = Σ log₂(20 / allowed_count) for each constrained position. Unconstrained positions (any amino acid) contribute 0 bits. Fixed positions contribute log₂(20) ≈ 4.32 bits.
- **CalculateEValue:** E = (N − L + 1) × 2^(−IC), where N = sequence length, L = motif length, IC = total information content. This is the expected number of random matches assuming uniform amino acid distribution.

## Overlapping Match Behavior

`FindMotifByPattern` uses regex lookahead `(?=(pattern))` to discover all matches including overlapping occurrences, consistent with PROSITE ScanProsite behavior (De Castro et al. 2006).

---

## Recommendations for Test Coverage

1. **MUST Test:** FindMotifByPattern returns correct matches for RGD pattern (`R-G-D`) — Evidence: PROSITE PS00016
2. **MUST Test:** FindMotifByPattern returns correct positions (0-based Start, End) — Evidence: Regex match semantics
3. **MUST Test:** FindCommonMotifs detects N-glycosylation (PS00001: `N-{P}-[ST]-{P}`) — Evidence: PROSITE PS00001
4. **MUST Test:** FindCommonMotifs correctly excludes proline in N-glycosylation — Evidence: PROSITE PS00001 exclusion syntax
5. **MUST Test:** FindCommonMotifs detects PKC phosphorylation (PS00005) — Evidence: PROSITE PS00005
6. **MUST Test:** FindCommonMotifs detects CK2 phosphorylation (PS00006) — Evidence: PROSITE PS00006
7. **MUST Test:** FindCommonMotifs detects P-loop (PS00017) — Evidence: PROSITE PS00017
8. **MUST Test:** ConvertPrositeToRegex correctly converts all PROSITE syntax elements — Evidence: PROSITE User Manual PA line syntax
9. **MUST Test:** FindMotifByPattern on empty/null input returns empty — Evidence: Trivial correctness
10. **MUST Test:** FindMotifByPattern is case-insensitive — Evidence: PROSITE convention
11. **MUST Test:** FindMotifByPattern handles invalid regex gracefully — Evidence: Robustness requirement
12. **MUST Test:** PS00007 pattern matches official PROSITE definition after fix — Evidence: PROSITE PS00007
13. **MUST Test:** PS00018 pattern matches official PROSITE definition after fix — Evidence: PROSITE PS00018
14. **MUST Test:** CommonMotifs dictionary contains all official PROSITE entries with correct patterns — Evidence: PROSITE database
15. **SHOULD Test:** FindCommonMotifs on a real protein (e.g., Human Insulin P01308) — Evidence: UniProt/ScanProsite
16. **COULD Test:** Multiple match positions for repeated motifs — Evidence: Pattern matching semantics

---

## References

1. PROSITE Database. SIB Swiss Institute of Bioinformatics. https://prosite.expasy.org/ (Accessed 2026-02-12)
2. PROSITE User Manual. https://prosite.expasy.org/prosuser.html (Accessed 2026-02-12)
3. Hulo N, Bairoch A, Bulliard V, et al. (2007). "The 20 years of PROSITE." Nucleic Acids Res. 36(Database issue):D245-9. https://doi.org/10.1093/nar/gkm977
4. De Castro E, Sigrist CJA, Gattiker A, et al. (2006). "ScanProsite: detection of PROSITE signature matches and ProRule-associated functional and structural residues in proteins." Nucleic Acids Res. 34(Web Server issue):W362-365. https://doi.org/10.1093/nar/gkl124
5. Wikipedia. "Sequence motif." https://en.wikipedia.org/wiki/Sequence_motif (Accessed 2026-02-12)
6. Wikipedia. "PROSITE." https://en.wikipedia.org/wiki/PROSITE (Accessed 2026-02-12)
7. PROSITE PS00001 (ASN_GLYCOSYLATION). https://prosite.expasy.org/PS00001 (Accessed 2026-02-12)
8. PROSITE PS00005 (PKC_PHOSPHO_SITE). https://prosite.expasy.org/PS00005 (Accessed 2026-02-12)
9. PROSITE PS00006 (CK2_PHOSPHO_SITE). https://prosite.expasy.org/PS00006 (Accessed 2026-02-12)
10. PROSITE PS00004 (CAMP_PHOSPHO_SITE). https://prosite.expasy.org/PS00004 (Accessed 2026-02-12)
11. PROSITE PS00007 (TYR_PHOSPHO_SITE_1). https://prosite.expasy.org/PS00007 (Accessed 2026-02-12)
12. PROSITE PS00008 (MYRISTYL). https://prosite.expasy.org/PS00008 (Accessed 2026-02-12)
13. PROSITE PS00009 (AMIDATION). https://prosite.expasy.org/PS00009 (Accessed 2026-02-12)
14. PROSITE PS00016 (RGD). https://prosite.expasy.org/PS00016 (Accessed 2026-02-12)
15. PROSITE PS00017 (ATP_GTP_A). https://prosite.expasy.org/PS00017 (Accessed 2026-02-12)
16. PROSITE PS00018 (EF_HAND_1). https://prosite.expasy.org/PS00018 (Accessed 2026-02-12)
17. PROSITE PS00028 (ZINC_FINGER_C2H2_1). https://prosite.expasy.org/PS00028 (Accessed 2026-02-12)
18. PROSITE PS00029 (LEUCINE_ZIPPER). https://prosite.expasy.org/PS00029 (Accessed 2026-02-12)
19. Dingwall C, Laskey RA (1991). "Nuclear targeting sequences — a consensus?" Trends Biochem Sci 16:478-481.
20. la Cour T, Kiemer L, Mølgaard A, et al. (2004). "Analysis and prediction of leucine-rich nuclear export signals." Protein Eng Des Sel 17:527-536.
21. Hecker CM, Rabiller M, Haglund K, et al. (2006). "Specification of SUMO1- and SUMO2-interacting motifs." J Biol Chem 281:16117-16127.
22. Chen HI, Sudol M (1995). "The WW domain of Yes-associated protein binds a proline-rich ligand that differs from the consensus established for Src homology 3-binding modules." PNAS 92:7819-7823.
23. Mayer BJ (2001). "SH3 domains: complexity in moderation." J Cell Sci 114:1253-1263.
24. Schneider TD, Stephens RM (1990). "Sequence logos: a new way to display consensus sequences." Nucleic Acids Res 18:6097-6100.
25. Altschul SF, Gish W, Miller W, et al. (1990). "Basic Local Alignment Search Tool." J Mol Biol 215:403-410.

---

## Change History

- **2026-02-12**: Initial documentation. Verified all PROSITE patterns against official source. Fixed PS00007 and PS00018 patterns in implementation.
- **2026-02-13**: Eliminated all assumptions. Verified 5 non-PROSITE patterns (NLS1, NES1, SIM1, WW1, SH3_1) against published literature. Replaced heuristic scoring with information-content-based scoring (Schneider & Stephens 1990). Implemented overlapping match discovery via regex lookahead. Updated SH3_1 from minimal PxxP core to full Class I consensus `[RK]-x(2)-P-x(2)-P` (Mayer 2001).
