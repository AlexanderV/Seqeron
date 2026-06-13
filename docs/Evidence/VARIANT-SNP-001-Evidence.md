# Evidence Artifact: VARIANT-SNP-001

**Test Unit ID:** VARIANT-SNP-001
**Algorithm:** SNP Detection (single-nucleotide substitution identification between a reference and a query sequence — alignment-based `FindSnps` and positional/Hamming-style `FindSnpsDirect`, with transition/transversion classification)
**Date Collected:** 2026-06-13

---

## Online Sources

### The Variant Call Format Specification VCFv4.3 (samtools/hts-specs)

**URL:** https://raw.githubusercontent.com/samtools/hts-specs/master/VCFv4.3.tex
**Accessed:** 2026-06-13
**Authority rank:** 2 (official specification / standard)

**How retrieved:** `WebSearch` for "VCF specification SNP single nucleotide polymorphism definition substitution REF ALT samtools hts-specs" returned the hts-specs repository; `WebFetch` of the master `VCFv4.3.tex` source asking for the simple-SNP example record, the POS definition, the REF definition, and the SNP/indel storage sentence.

**Key Extracted Points:**

1. **Simple SNP example (§1.1):** the record `20  14370  rs6054257  G  A  29  PASS  NS=3;DP=14;AF=0.5;DB;H2` is described as "a good simple SNP" — a single-base REF `G` substituted by a single-base ALT `A`, no padding base. (Retrieved verbatim from §1.1 of the fetched `.tex`.)
2. **POS field (1-based):** "The reference position, with the 1st base having position 1." (field definition, fetched verbatim)
3. **REF field alphabet / case:** "Each base must be one of A,C,G,T,N (case insensitive)." (REF field definition, fetched verbatim) — therefore SNP base comparison and transition/transversion classification must be case-insensitive.
4. **SNP is a single-base substitution:** the spec encodes a SNP as a single-base REF allele and a single-base ALT allele differing at one position (the §1.1 `G`→`A` example). A position where REF == ALT is not a variant.

### Wikipedia — Transversion (uses primary Futuyma 2013)

**URL:** https://en.wikipedia.org/wiki/Transversion
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primary source Futuyma, *Evolution*, 2013)

**How retrieved:** `WebFetch` of the article URL asking for the verbatim transversion definition and the explicit enumeration of the four transversions.

**Key Extracted Points:**

1. **Definition (verbatim):** a transversion is "a point mutation in DNA in which a single (two ring) purine (A or G) is changed for a (one ring) pyrimidine (T or C), or vice versa."
2. **The four transversions (verbatim list):** A↔C, A↔T, G↔C, G↔T.
3. **Counting:** four transversions are possible per base vs. one transition per base; transitions nonetheless occur more frequently (transition/transversion bias).

### Hamming Distance as a Concept in DNA Molecular Recognition (PMC5410656)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5410656/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary)

**How retrieved:** `WebSearch` for "Hamming distance positional substitution two equal length sequences point differences mismatch bioinformatics"; `WebFetch` of the PMC article asking for the verbatim Hamming-distance definition.

**Key Extracted Points:**

1. **Definition (verbatim):** "the number of positions that two codewords of the same length differ is the Hamming distance." (article cites Bierbrauer J., *Introduction to Coding Theory*, 2nd ed., Taylor & Francis, 2016.)
2. **Consequence for this unit:** positional SNP detection (`FindSnpsDirect`) over two same-length sequences is exactly the enumeration of the Hamming-mismatch positions; each mismatched index is one substitution (one SNP). No gap handling is involved — substitutions only.

### Collins DW, Jukes TH (1994) — Rates of transition and transversion

**URL:** https://pubmed.ncbi.nlm.nih.gov/8034311/  (DOI: 10.1006/geno.1994.1192)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary)

**How retrieved:** Cross-referenced from the sibling VARIANT-CALL-001 evidence; the PubMed record establishes that transitional silent substitution rates (1.71×10⁻⁹ yr⁻¹) exceed transversional rates (1.22×10⁻⁹ yr⁻¹) at fourfold-degenerate sites — the empirical basis for Ti/Tv > 0.5. Used here only to ground that transition vs transversion is a real, asymmetric classification, not for any expected numeric test value.

---

## Documented Corner Cases and Failure Modes

### From VCFv4.3 spec

1. **REF == ALT is not a variant:** a position whose query base equals the reference base is a match, not a SNP (a SNP is a *substitution*). The detector must skip equal columns.
2. **Case insensitivity:** REF/ALT bases are A,C,G,T,N case-insensitive, so `a`/`g` must classify identically to `A`/`G`.

### From the Hamming-distance definition (PMC5410656)

1. **Equal-length precondition:** the Hamming distance is defined only for strings "of the same length". When the two inputs differ in length, positional comparison is well-defined only over the common prefix; the trailing region of the longer sequence is not a substitution and is out of scope for `FindSnpsDirect` (indels are handled by VARIANT-INDEL-001, not here).

---

## Test Datasets

### Dataset: VCF v4.3 §1.1 canonical simple-SNP record

**Source:** VCFv4.3 specification §1.1 (samtools/hts-specs).

| Record | POS (1-based) | REF | ALT | Meaning (per spec text) |
|--------|---------------|-----|-----|-------------------------|
| simple SNP | 14370 | G | A | single-base substitution G→A ("a good simple SNP") |

In the repository's 0-based per-position model, a SNP at a single mismatched index `i` carries `Position = i`, `ReferenceAllele = ref[i]`, `AlternateAllele = query[i]`.

### Dataset: Positional substitution detection (Hamming mismatches)

**Source:** Hamming-distance definition (PMC5410656): mismatch positions between two equal-length strings.

| Reference | Query | Mismatch positions (0-based) | SNP alleles (ref→alt) |
|-----------|-------|------------------------------|------------------------|
| `ATGC`    | `ATTC`| {2}                          | G→T |
| `AAAA`    | `TGTA`| {0,1,2}                      | A→T, A→G, A→T |
| `ATGC`    | `ATGC`| {} (Hamming distance 0)      | (none) |

### Dataset: Transition / Transversion classification of SNP alleles

**Source:** Wikipedia Transversion (Futuyma 2013) + transition definition (purine↔purine / pyrimidine↔pyrimidine).

| REF→ALT | Class | Reason |
|---------|-------|--------|
| A→G | Transition | purine ↔ purine |
| C→T | Transition | pyrimidine ↔ pyrimidine |
| A→C | Transversion | purine ↔ pyrimidine |
| G→T | Transversion | purine ↔ pyrimidine |

---

## Assumptions

1. **ASSUMPTION: Unequal-length inputs to `FindSnpsDirect` compare only the common prefix.** The Hamming distance is defined for equal-length strings only (PMC5410656). The repository contract iterates `min(reference.Length, query.Length)` and reports substitutions over the common prefix; the trailing region of the longer input is not reported as a SNP (it is indel territory, handled by VARIANT-INDEL-001). No source mandates substitution semantics beyond the common length, so this prefix behavior is the documented contract, not a defect. Tests assert it explicitly.
2. **ASSUMPTION: 0-based `Position` in the in-memory `Variant`.** The VCF serialization is 1-based (spec field POS), but the in-memory `Variant.Position` reported by the detector is 0-based; VCF 1-based POS is produced only by `ToVcfLines` (out of scope here). This matches the existing sibling contract (VARIANT-CALL-001) and is internally consistent. Not a source-governed value for the in-memory model.

---

## Recommendations for Test Coverage

1. **MUST Test:** identical sequences produce zero SNPs (Hamming distance 0) — Evidence: PMC5410656 (Hamming distance = number of differing positions).
2. **MUST Test:** a single substitution yields exactly one SNP with exact 0-based position and ref/alt alleles — Evidence: VCFv4.3 §1.1 simple-SNP record; PMC5410656.
3. **MUST Test:** multiple substitutions yield one SNP per mismatched position with exact positions and alleles — Evidence: PMC5410656 (each mismatch is one substitution).
4. **MUST Test:** `FindSnps` (alignment-based) reports SNPs only and no indels for a substitution-only input — Evidence: VCFv4.3 (SNP = substitution; insertions/deletions are distinct classes).
5. **MUST Test:** `FindSnpsDirect` over unequal-length inputs compares only the common prefix — Evidence: PMC5410656 equal-length precondition (ASSUMPTION-1).
6. **SHOULD Test:** null inputs to `FindSnps` throw `ArgumentNullException`; empty inputs to `FindSnpsDirect` yield empty — Rationale: documented input-validation contract.
7. **SHOULD Test:** every variant emitted by `FindSnps`/`FindSnpsDirect` is of type SNP with distinct ref/alt alleles — Rationale: structural invariant of SNP detection (a SNP is a substitution, ref ≠ alt).
8. **COULD Test:** property — the SNP count from `FindSnpsDirect` over two equal-length sequences equals their Hamming distance — Evidence: PMC5410656 (Hamming distance = count of differing positions).

---

## References

1. The SAM/BCF/VCF working group. (2024). The Variant Call Format Specification VCFv4.3. samtools/hts-specs. https://raw.githubusercontent.com/samtools/hts-specs/master/VCFv4.3.tex (also https://samtools.github.io/hts-specs/VCFv4.3.pdf)
2. Wikipedia. Transversion. https://en.wikipedia.org/wiki/Transversion (accessed 2026-06-13; primary: Futuyma DJ, *Evolution*, 3rd ed., 2013, ISBN 978-1605351155)
3. Acharya T, et al. (2017). Hamming Distance as a Concept in DNA Molecular Recognition. *ACS Omega*. PMC5410656. https://pmc.ncbi.nlm.nih.gov/articles/PMC5410656/ (Hamming distance definition; cites Bierbrauer J., *Introduction to Coding Theory*, 2nd ed., Taylor & Francis, 2016)
4. Collins DW, Jukes TH (1994). Rates of transition and transversion in coding sequences since the human-rodent divergence. *Genomics* 20(3):386–396. https://doi.org/10.1006/geno.1994.1192 (PMID 8034311)

---

## Change History

- **2026-06-13**: Initial documentation.
