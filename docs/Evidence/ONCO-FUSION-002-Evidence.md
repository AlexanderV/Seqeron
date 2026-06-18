# Evidence Artifact: ONCO-FUSION-002

**Test Unit ID:** ONCO-FUSION-002
**Algorithm:** Known Fusion Database Lookup (HGNC gene-fusion designation + caller-supplied known-fusion match)
**Date Collected:** 2026-06-14

---

## Online Sources

### Bruford et al. (2021) — HUGO Gene Nomenclature Committee (HGNC) recommendations for the designation of gene fusions

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/
**Accessed:** 2026-06-14
**Authority rank:** 2 (official nomenclature standard, HGNC) — also peer-reviewed (Leukemia)

**Retrieval:** WebSearch query `HGVS gene fusion nomenclature double colon "::" separator recommendation`
located the Leukemia/PMC article; WebFetch of the PMC full text returned the verbatim quotes below.

**Key Extracted Points:**

1. **Separator:** "HGNC recommends that a new separator—a double colon (::)—be used in describing gene
   fusions, e.g., *BCR*::*ABL1*." The double colon was chosen as a unique, instantly recognizable separator,
   replacing the previously ambiguous hyphen (-) and forward slash (/).
2. **Partner order (directional):** "the 5′ partner gene should always be listed first in the description of a
   fusion gene, i.e., before the double colon, irrespective of chromosomal location or the orientation of the
   gene." → The designation is directional: the 5′ gene is written before `::`, the 3′ gene after.
3. **Approved symbols:** "genes involved in fusions should be designated by their HGNC approved gene symbols."
4. **Read-through transcripts use a hyphen, not `::`:** "The HGNC has approved the use of the hyphen
   separator … for denoting readthrough transcripts, e.g., *INS-IGF2*." The double colon is reserved for true
   fusion genes; read-throughs keep the hyphen.
5. **Worked example (5′/3′ assignment):** "in the *BCR*::*ABL1* fusion gene—the outcome of the translocation
   t(9;22)(q34.1;q11.2)—the *BCR* gene in chromosome 22 is the 5′ gene, the *ABL1* gene from chromosome 9 is
   the 3′ gene." → `GetFusionAnnotation("BCR","ABL1")` must yield `BCR::ABL1`.

### Recommendations for future extensions to the HGNC gene fusion nomenclature (Leukemia 2021, response/consortium letter)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC8632684/
**Accessed:** 2026-06-14
**Authority rank:** 2 (consortium consensus endorsing the HGNC standard)

**Retrieval:** WebSearch query `Bruford 2021 HGNC guidelines gene fusion nomenclature 5' 3' partner ordering`;
WebFetch of the PMC text.

**Key Extracted Points:**

1. **Endorsement of `::`:** "The primary recommendation from this manuscript is for a double-colon ('::')
   delimiter for indicating fusions." / "Our consortium thus supports the use of double-colon for fusion
   representation." Confirms `::` as the cross-resource standard separator.

---

## Documented Corner Cases and Failure Modes

### From Bruford et al. (2021)

1. **Direction matters:** because the 5′ gene is always written first, the designations `A::B` and `B::A`
   describe two *different* fusions (e.g. a reciprocal fusion). A directional lookup must not treat them as
   equal. (Derived from point 2: order is fixed by 5′/3′ role, not alphabetical.)
2. **Hyphen ≠ double colon:** `INS-IGF2` (read-through, hyphen) must not be confused with a fusion
   designation; this unit produces `::` designations for true fusions only.

---

## Test Datasets

### Dataset: HGNC worked example (BCR::ABL1)

**Source:** Bruford et al. (2021), PMC8550944 (point 5 above).

| Parameter | Value |
|-----------|-------|
| 5′ gene | BCR (chromosome 22) |
| 3′ gene | ABL1 (chromosome 9) |
| Canonical designation | `BCR::ABL1` |
| Reciprocal (different fusion) | `ABL1::BCR` |

### Dataset: Caller-supplied known-fusion list (illustrative; lists are NOT bundled with the library)

**Source:** designations are formatted per Bruford et al. (2021); membership is caller-supplied
(Mitelman / COSMIC / ChimerDB are NOT fabricated or bundled — see Assumption 1).

| Designation | Annotation (caller-supplied label) |
|-------------|------------------------------------|
| `BCR::ABL1` | Chronic myeloid leukemia driver |
| `EML4::ALK` | NSCLC driver, ALK TKI target |

---

## Assumptions

1. **ASSUMPTION: Known-fusion membership is caller-supplied** — The library does not bundle Mitelman, COSMIC,
   or ChimerDB content (licensing + curation are out of scope; per the unit mandate). `MatchKnownFusions`
   takes the known-fusion set as a parameter. Only the *designation format* and the *directional 5′/3′ keying*
   are evidence-defined (Bruford et al. 2021); the set contents are the caller's responsibility. This makes
   the unit a Framework algorithm: format/keying are source-backed, data is supplied.
2. **ASSUMPTION: Symbol case-insensitivity** — HGNC approved symbols are case-defined (uppercase Latin), but
   real inputs vary in case. Lookup compares symbols case-insensitively (ordinal-ignore-case) while preserving
   directionality. This affects matching only, not the formal `::`/order rule. (Not contradicted by the source;
   flagged as an assumption because Bruford et al. specify approved symbols without addressing case folding.)

---

## Recommendations for Test Coverage

1. **MUST Test:** `GetFusionAnnotation("BCR","ABL1") == "BCR::ABL1"` (5′ first, `::` separator) — Evidence: Bruford et al. (2021) points 1,2,5.
2. **MUST Test:** Direction matters — `GetFusionAnnotation("ABL1","BCR") == "ABL1::BCR"` ≠ `BCR::ABL1` — Evidence: point 2.
3. **MUST Test:** `MatchKnownFusions` returns the annotation for a designation present in the caller's set, keyed by `5′::3′` — Evidence: points 1,2.
4. **MUST Test:** `MatchKnownFusions` returns no-match when only the reciprocal `B::A` is in the set (directional) — Evidence: corner case 1.
5. **SHOULD Test:** case-insensitive symbol match (`bcr`/`abl1`) — Rationale: Assumption 2.
6. **MUST Test:** null/empty symbol and null arguments rejected — Rationale: input-validation contract (mirrors sibling methods).
7. **COULD Test:** annotation string round-trips through `MatchKnownFusions` on a `FusionCall` produced by `DetectFusions` — Rationale: integration with ONCO-FUSION-001.

---

## References

1. Bruford EA, Antonescu CR, Carroll AJ, et al. (2021). HUGO Gene Nomenclature Committee (HGNC) recommendations for the designation of gene fusions. *Leukemia* 35(11):3040–3043. https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/ (DOI: 10.1038/s41375-021-01436-6)
2. Various authors (2021). Recommendations for future extensions to the HGNC gene fusion nomenclature. *Leukemia* 35(11):3044–3045. https://pmc.ncbi.nlm.nih.gov/articles/PMC8632684/ (DOI: 10.1038/s41375-021-01438-4)

---

## Change History

- **2026-06-14**: Initial documentation.
