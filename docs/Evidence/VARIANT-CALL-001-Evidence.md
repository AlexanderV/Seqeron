# Evidence Artifact: VARIANT-CALL-001

**Test Unit ID:** VARIANT-CALL-001
**Algorithm:** Variant Detection (SNP / Insertion / Deletion calling from a reference↔query comparison, with transition/transversion classification)
**Date Collected:** 2026-06-13

---

## Online Sources

### The Variant Call Format Specification VCFv4.3 (samtools/hts-specs)

**URL:** https://samtools.github.io/hts-specs/VCFv4.3.pdf
**Accessed:** 2026-06-13
**Authority rank:** 2 (official specification / standard)

**How retrieved:** `WebSearch` for "VCF specification v4.3 variant representation REF ALT POS samtools hts-specs" returned the PDF URL; `WebFetch` of the PDF saved the binary locally; the text was decompressed from the PDF FlateDecode streams and the field-definition section read directly.

**Key Extracted Points:**

1. **POS field (1-based):** "POS - position: The reference position, with the 1st base having position 1." (Section 1.4.1 / 1.6.1, field 2.)
2. **REF field:** "REF - reference base(s): Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted. The value in the POS field refers to the position of the first base in the String." (field 4)
3. **Indel padding base:** "For simple insertions and deletions in which either the REF or one of the ALT alleles would otherwise be null/empty, the REF and ALT Strings must include the base before the event (which must be reflected in the POS field), unless the event occurs at position 1 on the contig in which case it must include the base after the event; this padding base is not required (although it is permitted) for e.g. complex substitutions or other events where all alleles have at least one base represented in their Strings." (field 4)
4. **ALT field alphabet:** "ALT - alternate base(s): … the ALT field must … match the regular expression `^([ACGTNacgtn]+|\*|\.)$`." (field 5)
5. **Worked example (the canonical VCF example, §1.1):** record `20  1234567  microsat1  GTC  G,GTCT  50  PASS …` — the spec text states this is "a microsatellite with two alternative alleles, one a deletion of 2 bases (TC), and the other an insertion of one base (T)." This shows the padding-base convention: REF `GTC` vs ALT `G` is a 2-base deletion (TC), and REF `GTC` vs ALT `GTCT` is a 1-base insertion (T), both anchored at the preceding base `G`.
6. **Simple SNP example (§1.1):** record `20  14370  rs6054257  G  A  29  PASS …` — described as "a good simple SNP": single-base REF `G`, single-base ALT `A`, no padding base.

### Wikipedia — Transition (genetics)

**URL:** https://en.wikipedia.org/wiki/Transition_(genetics)
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primary sources; primaries used below)

**How retrieved:** `WebSearch` for "transition transversion definition purine pyrimidine SNP molecular evolution"; `WebFetch` of the article URL.

**Key Extracted Points:**

1. **Definition:** A transition is "a point mutation that changes a purine nucleotide to another purine (A ↔ G) or a pyrimidine nucleotide to another pyrimidine (C ↔ T)."
2. **SNP frequency:** "approximately two out of three single nucleotide polymorphisms (SNPs) are transitions." Cited to Collins & Jukes (1994) (primary, below).
3. **Bias:** "although there are twice as many possible transversions, transitions occur more often in genomes, a pattern known as transition/transversion bias." Cited to Ebersberger et al. (2002).

### Wikipedia — Transversion

**URL:** https://en.wikipedia.org/wiki/Transversion
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primary source Futuyma, *Evolution*, 2013)

**How retrieved:** `WebFetch` of the article URL.

**Key Extracted Points:**

1. **Definition:** A transversion is "a point mutation in DNA in which a single (two ring) purine (A or G) is changed for a (one ring) pyrimidine (T or C), or vice versa."
2. **Enumerated transversions:** "The transversions are: A↔C, A↔T, G↔C, and G↔T." (4 possible transversions vs 2 possible transitions per base.)

### Collins DW, Jukes TH (1994) — Rates of transition and transversion

**URL:** https://pubmed.ncbi.nlm.nih.gov/8034311/  (DOI: 10.1006/geno.1994.1192)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary)

**How retrieved:** `WebSearch` for the title; `WebFetch` of the PubMed record.

**Key Extracted Points:**

1. **Citation confirmed:** Collins DW, Jukes TH (1994). "Rates of transition and transversion in coding sequences since the human-rodent divergence." *Genomics* 20(3):386–396. DOI 10.1006/geno.1994.1192. PMID 8034311.
2. **Transitional bias:** abstract states "The rates of transitional and transversional silent substitutions in fourfold degenerate sites are estimated as 1.71 × 10⁻⁹ and 1.22 × 10⁻⁹ site⁻¹ year⁻¹, respectively" — i.e. transitions exceed transversions (the basis for Ti/Tv > 0.5 expectations).

### Tan A, Abecasis GR, Kang HM (2015) — Unified representation of genetic variants

**URL:** https://academic.oup.com/bioinformatics/article/31/13/2202/196142  (DOI: 10.1093/bioinformatics/btv112)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary)

**How retrieved:** `WebSearch` for "left alignment normalization indel variant calling … Tan Abecasis 2015"; `WebFetch` of the Oxford Academic article.

**Key Extracted Points:**

1. **Normalization:** "A VCF entry is *normalized* if and only if it is left aligned and parsimonious."
2. **Left alignment:** "A VCF entry is *left aligned* if and only if its base position is smallest among all potential VCF entries having the same allele length and representing the same variant."
3. **Parsimony:** "A VCF entry is *parsimonious* if and only if the entry has the shortest allele length among all VCF entries representing the same variant."
4. **Consequence for this unit:** an indel detected from an arbitrary alignment is not guaranteed to be in the canonical (normalized) representation; the same biological indel can be written multiple ways. This bounds the correctness claim of the alignment-based caller (see Assumptions).

### Danecek P, et al. (2011) — The variant call format and VCFtools

**URL:** https://pubmed.ncbi.nlm.nih.gov/21653522/  (DOI: 10.1093/bioinformatics/btr330)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary; the VCF reference paper)

**How retrieved:** `WebSearch` for "Danecek 2011 variant call format VCFtools Bioinformatics DOI".

**Key Extracted Points:**

1. **Citation confirmed:** Danecek P, Auton A, Abecasis G, et al. (2011). "The variant call format and VCFtools." *Bioinformatics* 27(15):2156–2158. DOI 10.1093/bioinformatics/btr330. PMID 21653522.
2. **Scope:** "VCF is a generic format for storing DNA polymorphism data such as SNPs, insertions, deletions and structural variants" — confirms the three variant classes this unit must distinguish (SNP / insertion / deletion).

---

## Documented Corner Cases and Failure Modes

### From VCFv4.3 spec

1. **Empty allele forbidden in REF/ALT:** a pure insertion (REF would be empty) or pure deletion (ALT would be empty) MUST carry a padding base. The repository's column-based caller instead emits a per-column event with the `"-"` gap sentinel; this is an internal representation, not VCF output (see Assumptions ASM-01).
2. **Case insensitivity:** REF/ALT bases are case-insensitive (A,C,G,T,N). The classification of transition/transversion must therefore be case-insensitive.

### From Tan et al. (2015)

1. **Non-unique indel placement:** in a low-complexity / repeated region, the alignment that produces the indel column is not unique, so the reported indel position is implementation-dependent unless normalized (left-aligned + parsimonious).

---

## Test Datasets

### Dataset: VCF v4.3 §1.1 canonical example records

**Source:** VCFv4.3 specification §1.1 (samtools/hts-specs).

| Record | POS | REF | ALT | Meaning (per spec text) |
|--------|-----|-----|-----|-------------------------|
| simple SNP | 14370 | G | A | single-base substitution G→A (a transition) |
| microsatellite (del) | 1234567 | GTC | G | deletion of 2 bases (TC) |
| microsatellite (ins) | 1234567 | GTC | GTCT | insertion of one base (T) |

### Dataset: Transition / Transversion classification table

**Source:** Wikipedia Transition (genetics) + Transversion, citing Futuyma (2013) and Collins & Jukes (1994).

| REF→ALT | Class | Reason |
|---------|-------|--------|
| A→G | Transition | purine ↔ purine |
| G→A | Transition | purine ↔ purine |
| C→T | Transition | pyrimidine ↔ pyrimidine |
| T→C | Transition | pyrimidine ↔ pyrimidine |
| A→C | Transversion | purine ↔ pyrimidine |
| A→T | Transversion | purine ↔ pyrimidine |
| G→C | Transversion | purine ↔ pyrimidine |
| G→T | Transversion | purine ↔ pyrimidine |
| C→A | Transversion | pyrimidine ↔ purine |
| T→A | Transversion | pyrimidine ↔ purine |

### Dataset: Ti/Tv ratio derivation

**Source:** definition of Ti/Tv = (#transitions)/(#transversions), built from the classes above.

| Variants | #Ti | #Tv | Ti/Tv |
|----------|-----|-----|-------|
| {A→G, A→C} | 1 | 1 | 1.0 |
| {A→G, C→T} | 2 | 0 | 0 by repo convention (no transversions; ratio undefined → 0) |

---

## Assumptions

1. **ASSUMPTION: Internal gap-sentinel representation for indels.** The repository's `CallVariantsFromAlignment` reports a per-column indel using the `"-"` gap character for the absent allele and a 0-based `Position`, rather than the VCF padded-allele, 1-based representation. The VCF spec (field 4) mandates a padding base and 1-based POS only for the *serialized VCF* (which `ToVcfLines` produces, out of scope here). The in-memory `Variant` model is an implementation choice not governed by a source; it is internally consistent and is the existing contract of sibling methods. — Not changed; documented.
2. **ASSUMPTION: Indels are not left-aligned / parsimony-normalized.** Per Tan et al. (2015) the canonical representation requires left-alignment and parsimony; the alignment-based caller reports the indel at the column produced by `SequenceAligner.GlobalAlign` without a normalization pass. This is correctness-affecting only for *position* in repeated regions, not for *variant counts/types*. Tests therefore assert counts/types/alleles on unambiguous inputs and assert position only where the alignment is unique. — Documented as a limitation, not a defect of the detection logic.
3. **ASSUMPTION: Ti/Tv with zero transversions returns 0.** The mathematically-undefined case (#Tv = 0) is mapped to 0 by the existing contract rather than throwing or returning +∞. No source mandates a specific sentinel; tested as the documented contract.

---

## Recommendations for Test Coverage

1. **MUST Test:** identical sequences produce zero variants — Evidence: Danecek 2011 (a variant is a difference from reference).
2. **MUST Test:** a single substitution column yields exactly one SNP with correct REF/ALT/position — Evidence: VCFv4.3 simple-SNP example.
3. **MUST Test:** a one-base insertion column (`-`/base) yields exactly one Insertion with correct alleles and position bookkeeping — Evidence: VCFv4.3 microsatellite insertion example.
4. **MUST Test:** a one-base deletion column (base/`-`) yields exactly one Deletion — Evidence: VCFv4.3 microsatellite deletion example.
5. **MUST Test:** transition/transversion classification for all four transition pairs and representative transversions, case-insensitively — Evidence: Wikipedia Transition/Transversion (Futuyma 2013).
6. **MUST Test:** Ti/Tv ratio numerator/denominator from classified SNPs; Ti/Tv=0 when no transversions — Evidence: definition + Collins & Jukes (1994) bias.
7. **SHOULD Test:** mismatched aligned lengths throw `ArgumentException`; empty input yields empty — Rationale: documented contract / input validation.
8. **SHOULD Test:** non-SNP variant classifies as `Other` — Rationale: classification is defined only for SNPs.
9. **COULD Test:** property — every emitted variant's position is within reference bounds; SNP REF≠ALT always — Rationale: structural invariant (O(n×m) algorithm property test).

---

## References

1. Danecek P, Auton A, Abecasis G, et al. (2011). The variant call format and VCFtools. *Bioinformatics* 27(15):2156–2158. https://doi.org/10.1093/bioinformatics/btr330 (PMID 21653522)
2. The SAM/BCF/VCF working group. (2024). The Variant Call Format Specification VCFv4.3. samtools/hts-specs. https://samtools.github.io/hts-specs/VCFv4.3.pdf
3. Collins DW, Jukes TH (1994). Rates of transition and transversion in coding sequences since the human-rodent divergence. *Genomics* 20(3):386–396. https://doi.org/10.1006/geno.1994.1192 (PMID 8034311)
4. Tan A, Abecasis GR, Kang HM (2015). Unified representation of genetic variants. *Bioinformatics* 31(13):2202–2204. https://doi.org/10.1093/bioinformatics/btv112
5. Wikipedia. Transition (genetics). https://en.wikipedia.org/wiki/Transition_(genetics) (accessed 2026-06-13; primary: Collins & Jukes 1994)
6. Wikipedia. Transversion. https://en.wikipedia.org/wiki/Transversion (accessed 2026-06-13; primary: Futuyma DJ, *Evolution*, 3rd ed., 2013, ISBN 978-1605351155)

---

## Change History

- **2026-06-13**: Initial documentation.
