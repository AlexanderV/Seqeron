# Evidence Artifact: VARIANT-INDEL-001

**Test Unit ID:** VARIANT-INDEL-001
**Algorithm:** Indel Detection (identification of insertions and deletions between a reference and a query sequence — alignment-based `FindInsertions` / `FindDeletions`, filtering the aligned-column variant caller)
**Date Collected:** 2026-06-13

---

## Online Sources

### The Variant Call Format Specification VCFv4.3 (samtools/hts-specs)

**URL:** https://raw.githubusercontent.com/samtools/hts-specs/master/VCFv4.3.tex
**Accessed:** 2026-06-13
**Authority rank:** 2 (official specification / standard)

**How retrieved:** `WebFetch` of the master `VCFv4.3.tex` source asking for the verbatim indel examples (single-base insertion of A after position 3; single-base deletion of C at position 3), the padding-base rule, and the POS / REF field definitions.

**Key Extracted Points:**

1. **Single-base insertion example (verbatim):** "A single base insertion of A after position 3 becomes REF=C, ALT=CA". An insertion is encoded as REF = the single anchor base, ALT = the anchor base followed by the inserted base(s); ALT is longer than REF.
2. **Single-base deletion example (verbatim):** "A single base deletion of C at position 3 becomes REF=TC, ALT=T". A deletion is encoded as REF = anchor base followed by the deleted base(s), ALT = the single anchor base; REF is longer than ALT.
3. **Padding-base rule (verbatim):** "For simple insertions and deletions in which either the REF or one of the ALT alleles would otherwise be null/empty, the REF and ALT Strings must include the base before the event (which must be reflected in the POS field), unless the event occurs at position 1 on the contig in which case it must include the base after the event; this padding base is not required (although it is permitted) for e.g. complex substitutions or other events where all alleles have at least one base represented in their Strings." (REF field, field 4)
4. **POS field (1-based) (verbatim):** "The reference position, with the 1st base having position 1." — POS for an indel refers to the position of the (anchor) base, the first base of the REF String.
5. **REF field alphabet (verbatim):** "Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted. The value in the POS field refers to the position of the first base in the String."
6. **Canonical microsatellite example (§1.1):** record `20  1234567  microsat1  GTC  G,GTCT  50  PASS` — "a deletion of 2 bases (TC), and the other an insertion of one base (T)", both anchored at the preceding base `G`. Confirms: deletion ⇒ REF longer than ALT; insertion ⇒ ALT longer than REF.

### Tan A, Abecasis GR, Kang HM (2015) — Unified representation of genetic variants

**URL:** https://academic.oup.com/bioinformatics/article/31/13/2202/196142  (DOI: 10.1093/bioinformatics/btv112)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary)

**How retrieved:** `WebSearch` for "Tan Abecasis Kang 2015 Unified representation genetic variants Bioinformatics normalization left align parsimony" located the Oxford Academic article (DOI 10.1093/bioinformatics/btv112; PMID 25701572); `WebFetch` of the article asked for the normalized/left-aligned/parsimonious definitions and Algorithm 1.

**Key Extracted Points:**

1. **Normalized (verbatim):** "A VCF entry is *normalized* if and only if it is left aligned and parsimonious."
2. **Left aligned (verbatim):** "A VCF entry is *left aligned* if and only if its base position is smallest among all potential VCF entries having the same allele length and representing the same variant."
3. **Parsimonious (verbatim):** "A VCF entry is *parsimonious* if and only if the entry has the shortest allele length among all VCF entries representing the same variant."
4. **Algorithm 1 (verbatim):** `1: do; 2: if all alleles end with same nucleotide then; 3: truncate the rightmost nucleotide of each allele; 4: if any allele is length zero then; 5: extend all alleles by 1 nucleotide to the left; 6: while changes made in the VCF entry in the loop; 7: while all alleles start with same nucleotide and length ≥2 do; 8: truncate the leftmost nucleotide of each allele; 9: end while; 10: return the VCF entry`.
5. **Consequence for this unit:** the same biological indel can be written multiple ways; an indel emitted from an arbitrary alignment column is not guaranteed to be in the canonical normalized representation. This bounds the *position* claim of the alignment-based caller in repeated/low-complexity regions (see Assumptions). It does not affect indel *counts* or *types* on unambiguous inputs.

### minimal_representation — reference implementation of variant trimming (E. Minikel)

**URL:** https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation; implements Tan et al. 2015 Algorithm 1)

**How retrieved:** `WebSearch` surfaced `github.com/ericminikel/minimal_representation/blob/master/normalize.py`; `WebFetch` of the raw file asked for the trimming algorithm and worked examples.

**Key Extracted Points:**

1. **Suffix-then-prefix trimming (verbatim structure):** the `normalize` routine first removes shared *suffix* bases (`while ref[-1] == alt[-1]`, re-padding to the left when an allele empties and decrementing `pos`), then removes shared *prefix* bases (`while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]`, incrementing `pos`). This is the operational form of Tan et al. Algorithm 1.
2. **Worked example — CFTR p.F508del (verbatim test case):** `('7', 117199646, 'CTT', '-')` normalizes to `('7', 117199644, 'ATCT', 'A')` — a 3-base deletion, padded to a left anchor `A`, REF longer than ALT.
3. **Worked example — BRCA2 (verbatim test case):** `('13', 32914438, 'T', '-')` normalizes to `('13', 32914437, 'GT', 'G')` — a 1-base deletion, padded to a left anchor `G`.
4. **Consequence for this unit:** confirms the directional length invariant — a deletion has `len(REF) > len(ALT)`; an insertion has `len(ALT) > len(REF)` — and that pure (empty-allele) indels are padded with a left anchor base. The repository's in-memory model instead uses the `-` gap sentinel for the absent allele (see Assumptions).

### PharmCAT — Variant Normalization worked examples

**URL:** https://pharmcat.clinpgx.org/using/Variant-Normalization/
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference tool documentation; applies Tan et al. 2015 normalization)

**How retrieved:** `WebSearch` for VCF normalization worked examples surfaced the PharmCAT page; `WebFetch` asked for before/after POS/REF/ALT values.

**Key Extracted Points:**

1. **Left-alignment in a tandem repeat (verbatim):** before `POS 97740414, REF AATGA, ALT A` → after `POS 97740410, REF GATGA, ALT G`. Demonstrates that an un-normalized deletion can be shifted left to a smaller POS within a repeat; the *count* (one deletion) and *type* are unchanged, only the *position*/allele spelling change.
2. **Consequence for this unit:** corroborates Tan et al. (2015) ASM — indel position is alignment/normalization dependent in repeats. Tests assert position only where the alignment is unique; counts and types are asserted generally.

---

## Documented Corner Cases and Failure Modes

### From VCFv4.3 spec

1. **Empty allele forbidden in serialized VCF; padding base required:** a pure insertion (REF empty) or pure deletion (ALT empty) MUST carry a left anchor base in serialized VCF (field 4). The repository's column-based caller represents the absent allele with the `"-"` gap sentinel in the in-memory `Variant`; the padded form is produced only by `ToVcfLines` (out of scope). Documented as ASM-01.
2. **Directional length:** insertion ⇒ ALT longer than REF; deletion ⇒ REF longer than ALT (microsatellite example and the single-base examples). In the in-memory per-column model this manifests as: insertion column has `ReferenceAllele == "-"` and a one-base `AlternateAllele`; deletion column has a one-base `ReferenceAllele` and `AlternateAllele == "-"`.

### From Tan et al. (2015)

1. **Non-unique indel placement:** in a repeated / low-complexity region the alignment producing the indel column is not unique, so the reported indel position is implementation-dependent unless normalized (left-aligned + parsimonious). Affects *position* only, not *counts* or *types* (ASM-02).

---

## Test Datasets

### Dataset: VCF v4.3 single-base indel examples

**Source:** VCFv4.3 specification (REF field text, samtools/hts-specs).

| Event | POS (1-based) | REF | ALT | Length relation |
|-------|---------------|-----|-----|-----------------|
| single-base insertion of A after position 3 | 3 | C | CA | ALT longer than REF |
| single-base deletion of C at position 3 | 3 | TC | T | REF longer than ALT |
| microsatellite deletion (2 bp TC) | 1234567 | GTC | G | REF longer than ALT |
| microsatellite insertion (1 bp T) | 1234567 | GTC | GTCT | ALT longer than REF |

In the repository's per-column in-memory model the same events are represented column-wise: an inserted base is a column with `ReferenceAllele = "-"`, `AlternateAllele = <inserted base>`, `Type = Insertion`; a deleted base is a column with `ReferenceAllele = <deleted base>`, `AlternateAllele = "-"`, `Type = Deletion`.

### Dataset: Alignment-derived indel columns (unambiguous, non-repeat inputs)

**Source:** Needleman–Wunsch global alignment of reference vs query (`SequenceAligner.GlobalAlign`), with the VCF directional/length semantics above applied to each gap column. Inputs are chosen so the optimal alignment is unique (no internal repeat that would allow left-shifting), making position deterministic.

| Reference | Query | Event | Aligned column | Detected (Type, RefAllele, AltAllele) |
|-----------|-------|-------|----------------|----------------------------------------|
| `ATGCAT` | `ATGTCAT` | 1-base insertion of `T` after ref index 2 | ref `ATG-CAT` / qry `ATGTCAT` | (Insertion, `-`, `T`) |
| `ATGTCAT` | `ATGCAT` | 1-base deletion of `T` at ref index 3 | ref `ATGTCAT` / qry `ATG-CAT` | (Deletion, `T`, `-`) |
| `ATGC` | `ATGC` | no event | identity | (none) |

### Dataset: minimal_representation reference trimming (directional length invariant)

**Source:** ericminikel/minimal_representation `normalize.py` test cases (implements Tan et al. 2015 Algorithm 1).

| Input (chrom, pos, ref, alt) | Normalized (chrom, pos, ref, alt) | Invariant verified |
|------------------------------|-----------------------------------|--------------------|
| `(7, 117199646, CTT, -)` | `(7, 117199644, ATCT, A)` | deletion ⇒ len(REF) > len(ALT) |
| `(13, 32914438, T, -)` | `(13, 32914437, GT, G)` | deletion ⇒ len(REF) > len(ALT) |

---

## Assumptions

1. **ASSUMPTION: Internal gap-sentinel representation for indels (ASM-01).** The repository's `CallVariantsFromAlignment` reports each indel column using the `"-"` gap character for the absent allele and a 0-based `Position`, rather than the VCF padded-allele, 1-based representation. The VCF spec mandates a left padding base and 1-based POS only for *serialized VCF* (produced by `ToVcfLines`, out of scope here). The in-memory `Variant` model is an implementation choice not governed by a source; it is internally consistent and is the established contract of sibling methods (VARIANT-CALL-001). This is API/representation shape, not a correctness-affecting algorithm parameter (it changes no detection decision: counts and types are unchanged). — Documented; unchanged.
2. **ASSUMPTION: Indels are not left-aligned / parsimony-normalized (ASM-02).** Per Tan et al. (2015) the canonical representation requires left-alignment and parsimony; the alignment-based caller reports the indel at the column produced by `SequenceAligner.GlobalAlign` without a normalization pass. This is correctness-affecting only for *position* in repeated regions, not for *variant counts/types*. Tests assert counts/types/alleles generally and assert exact position only where the global alignment is provably unique (no repeat permitting a left shift). — Documented as a bounded limitation, not a defect of the detection logic.
3. **ASSUMPTION: 0-based `Position` in the in-memory `Variant`.** VCF serialization is 1-based (POS), but the in-memory `Variant.Position` is 0-based; 1-based POS is produced only by `ToVcfLines`. Matches the sibling contract (VARIANT-CALL-001) and is internally consistent. Not a source-governed value for the in-memory model.

---

## Recommendations for Test Coverage

1. **MUST Test:** a one-base insertion (a query base with no reference counterpart) yields exactly one Insertion with `ReferenceAllele == "-"`, the inserted base as `AlternateAllele`, and the correct reference position — Evidence: VCFv4.3 single-base insertion example (REF=C, ALT=CA; ALT longer than REF).
2. **MUST Test:** a one-base deletion (a reference base with no query counterpart) yields exactly one Deletion with the deleted base as `ReferenceAllele`, `AlternateAllele == "-"`, and the correct reference position — Evidence: VCFv4.3 single-base deletion example (REF=TC, ALT=T; REF longer than ALT).
3. **MUST Test:** `FindInsertions` returns insertions only (no deletions, no SNPs) and `FindDeletions` returns deletions only — Evidence: VCFv4.3 (insertion, deletion and SNP are distinct variant classes); these methods are filters over `CallVariants`.
4. **MUST Test:** a multi-base insertion / multi-base deletion is reported as the corresponding per-base indel columns at consecutive positions — Evidence: VCFv4.3 microsatellite example (2-bp deletion TC) + minimal_representation (multi-base indels).
5. **MUST Test:** identical sequences produce no insertions and no deletions — Evidence: an indel is a length-changing difference from the reference (VCFv4.3); identical sequences have none.
6. **SHOULD Test:** null reference / null query throw `ArgumentNullException` (propagated from `CallVariants`); a sequence with only substitutions yields no indels — Rationale: documented input-validation contract; insertion/deletion are distinct from substitution.
7. **SHOULD Test:** the directional length invariant — every Insertion has `AlternateAllele` longer than (non-gap) `ReferenceAllele` and every Deletion has `ReferenceAllele` longer than (non-gap) `AlternateAllele`, i.e. insertion ⇒ ref is the gap sentinel, deletion ⇒ alt is the gap sentinel — Evidence: VCFv4.3 directional length; minimal_representation.
8. **COULD Test:** property — over an input where a contiguous block of `k` bases is inserted (resp. deleted) in an otherwise unique alignment, `FindInsertions` (resp. `FindDeletions`) reports exactly `k` indel columns, all of the correct type — Evidence: VCFv4.3 (each absent/extra base is one indel event); ASM-02 bounds this to unique alignments.

---

## References

1. The SAM/BCF/VCF working group. (2024). The Variant Call Format Specification VCFv4.3. samtools/hts-specs. https://raw.githubusercontent.com/samtools/hts-specs/master/VCFv4.3.tex (also https://samtools.github.io/hts-specs/VCFv4.3.pdf)
2. Tan A, Abecasis GR, Kang HM (2015). Unified representation of genetic variants. *Bioinformatics* 31(13):2202–2204. https://doi.org/10.1093/bioinformatics/btv112 (PMID 25701572)
3. Minikel E. minimal_representation: find the minimal representation of a variant (implements Tan et al. 2015 normalization). https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py (accessed 2026-06-13)
4. PharmCAT. Variant Normalization. https://pharmcat.clinpgx.org/using/Variant-Normalization/ (accessed 2026-06-13; applies Tan et al. 2015)

---

## Change History

- **2026-06-13**: Initial documentation.
