# Evidence Artifact: PROTMOTIF-COMMON-001

**Test Unit ID:** PROTMOTIF-COMMON-001
**Algorithm:** Common Motif Finding (`ProteinMotifFinder.FindCommonMotifs`)
**Date Collected:** 2026-06-14

---

## Online Sources

### PROSITE entry PS00001 — N-glycosylation site (ASN_GLYCOSYLATION)

**URL:** https://prosite.expasy.org/PS00001
**Accessed:** 2026-06-14
**Authority rank:** 2 (official database / specification)

**Key Extracted Points:**

1. **Entry name:** Retrieved page reports entry name `ASN_GLYCOSYLATION`.
2. **Pattern (verbatim):** `N-{P}-[ST]-{P}.`
3. **Syntax interpretation:** `N` = Asparagine (required first position); `{P}` = any amino acid except Proline; `[ST]` = Serine or Threonine; `-` is the element separator.

### PROSITE entry PS00005 — Protein kinase C phosphorylation site (PKC_PHOSPHO_SITE)

**URL:** https://prosite.expasy.org/PS00005
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **Entry name:** `PKC_PHOSPHO_SITE`.
2. **Pattern (verbatim):** `[ST]-x-[RK]` — Ser/Thr, then any amino acid, then Arg/Lys.

### PROSITE entry PS00006 — Casein kinase II phosphorylation site (CK2_PHOSPHO_SITE)

**URL:** https://prosite.expasy.org/PS00006
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **Entry name:** `CK2_PHOSPHO_SITE`.
2. **Pattern (verbatim):** `[ST]-x(2)-[DE]` — Ser/Thr, then two arbitrary residues, then Asp/Glu.

### PROSITE entry PS00016 — Cell attachment sequence / RGD (RGD)

**URL:** https://prosite.expasy.org/PS00016
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **Entry name:** `RGD`.
2. **Pattern (verbatim):** `R-G-D.` — Arg, Gly, Asp in that fixed order.

### PROSITE entry PS00017 — ATP/GTP-binding site motif A / P-loop (ATP_GTP_A)

**URL:** https://prosite.expasy.org/PS00017
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **Entry name:** `ATP_GTP_A`.
2. **Pattern (verbatim):** `[AG]-x(4)-G-K-[ST]` — Ala/Gly, four arbitrary residues, then G, K, then Ser/Thr.

### PROSITE pattern-syntax rules (ScanProsite documentation / PROSITE User Manual)

**URL:** https://prosite.expasy.org/scanprosite/scanprosite_doc.html
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **Letters:** "The standard IUPAC one letter code for the amino acids is used in PROSITE."
2. **Wildcard:** "The symbol 'x' is used for a position where any amino acid is accepted."
3. **Allowed set:** "Ambiguities are indicated by listing the acceptable amino acids for a given position, between square brackets '[ ]'."
4. **Excluded set:** "Ambiguities are also indicated by listing between a pair of curly brackets '{ }' the amino acids that are not accepted."
5. **Separator:** "Each element in a pattern is separated from its neighbor by a '-'."
6. **Fixed repetition:** "x(3) corresponds to x-x-x".
7. **Variable repetition:** "x(2,4) corresponds to x-x or x-x-x or x-x-x-x".
8. **Overlap reporting (default):** Default ScanProsite behavior is "greedy, overlaps, no includes," meaning overlapping matches are reported unless one is entirely contained within another.
9. **Position reporting:** ScanProsite reports hits as `[hit start]-[hit stop]` using 1-based, inclusive coordinates.

---

## Documented Corner Cases and Failure Modes

### From PROSITE PS00001

1. **Proline exclusion:** `{P}` at positions 2 and 4 means an `N-P-[ST]-x` or `N-x-[ST]-P` window is NOT an N-glycosylation match — proline at an excluded position rejects the site.

### From ScanProsite documentation

1. **Overlapping occurrences:** Two matches of the same or different patterns that overlap are both reported (default "overlaps, no includes"); only a match fully contained inside another is suppressed.

---

## Test Datasets

### Dataset: Synthetic windows constructed to satisfy/violate each PROSITE pattern

**Source:** PROSITE patterns above; positions derived from 0-based string indexing as used by `ProteinMotifFinder` (`MotifMatch.Start`/`End` are 0-based, inclusive — repository convention; PROSITE reports 1-based).

| Sequence | Pattern (entry) | Expected match (0-based Start..End, Sequence) |
|----------|-----------------|------------------------------------------------|
| `AAAANFTAAAA` | PS00001 `N[^P][ST][^P]` | 4..7 `NFTA` |
| `AAAANPSAAAAANPTAAA` | PS00001 | no match (Pro at excluded pos 2) |
| `AAAAASARKAAA` | PS00005 `[ST].[RK]` | 5..7 `SAR` |
| `AAAASAAEASDEDAAA` | PS00006 `[ST].{2}[DE]` | 4..7 `SAAE`; 9..12 `SDED` |
| `AAAAAGXXXXGKSAAAA` | PS00017 `[AG].{4}GK[ST]` | 5..12 `GXXXXGKS` |
| `AARGDKK` | PS00016 `RGD` | 2..4 `RGD` |
| `RGDRGD` | PS00016 `RGD` (overlap test) | 0..2 `RGD`; 3..5 `RGD` (two non-overlapping) |

---

## Assumptions

1. **ASSUMPTION: 0-based inclusive coordinates** — PROSITE/ScanProsite reports 1-based coordinates; the repository `MotifMatch` records 0-based `Start`/`End` (matching sibling units PROTMOTIF-FIND-001 and PROTMOTIF-PATTERN-001). This is an API-shape convention, not a correctness-affecting parameter: the matched substring content and relative positions are identical; only the coordinate origin differs. Tests assert 0-based positions consistent with the repository convention.

---

## Recommendations for Test Coverage

1. **MUST Test:** `FindCommonMotifs` finds each canonical PROSITE motif (PS00001, PS00005, PS00006, PS00016, PS00017) at the exact 0-based position with the exact matched substring — Evidence: PROSITE entries above.
2. **MUST Test:** Proline exclusion `{P}` for PS00001 rejects `N-P-[ST]` windows — Evidence: PS00001 `{P}`.
3. **MUST Test:** Multiple distinct pattern types are aggregated from one sequence (whole-dictionary scan) — Evidence: `FindCommonMotifs` scans every entry in `CommonMotifs`.
4. **MUST Test:** Overlapping occurrences of a pattern are both reported — Evidence: ScanProsite default "overlaps".
5. **SHOULD Test:** Every returned `MotifMatch.Sequence` equals `protein.Substring(Start, End-Start+1)` (substring invariant) — Rationale: position/sequence consistency.
6. **SHOULD Test:** `MotifName`/`Pattern` (accession) on each match come from the matching `CommonMotifs` dictionary entry — Rationale: aggregation must carry the entry's identity.
7. **COULD Test:** Determinism (same input → identical ordered result) — Rationale: regex scan is deterministic.
8. **MUST Test:** null/empty input returns empty — Evidence: trivial guard.

---

## References

1. Sigrist CJA, de Castro E, Cerutti L, Cuche BA, Hulo N, Bridge A, Bougueleret L, Xenarios I. 2013. New and continuing developments at PROSITE. Nucleic Acids Research 41(D1):D344–D347. https://doi.org/10.1093/nar/gks1067
2. ExPASy PROSITE database. PS00001, PS00005, PS00006, PS00016, PS00017 entries. https://prosite.expasy.org/ (accessed 2026-06-14).
3. ExPASy ScanProsite documentation (pattern syntax, overlap and coordinate reporting). https://prosite.expasy.org/scanprosite/scanprosite_doc.html (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation.
