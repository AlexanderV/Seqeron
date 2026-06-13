# Evidence Artifact: VARIANT-ANNOT-001

**Test Unit ID:** VARIANT-ANNOT-001
**Algorithm:** Variant Annotation — functional impact / consequence prediction (VEP / Sequence Ontology)
**Date Collected:** 2026-06-13

---

## Online Sources

### Ensembl VEP — `Bio::EnsEMBL::Variation::Utils::Constants` (OverlapConsequence table)

**URL:** https://raw.githubusercontent.com/Ensembl/ensembl-variation/release/110/modules/Bio/EnsEMBL/Variation/Utils/Constants.pm
**Accessed:** 2026-06-13 (retrieved with `curl`; HTTP 200; identical content on `main` and `release/111`)
**Authority rank:** 3 (reference implementation — ensembl-variation, release/110)

**Key Extracted Points:**

1. **Consequence → IMPACT → rank table (verbatim from the file):** each consequence is an `OverlapConsequence` with `SO_term`, `SO_accession`, `impact`, and `rank`. The ordered table (rank, SO_term, impact, SO_accession):
   - 1 `transcript_ablation` HIGH SO:0001893
   - 2 `splice_acceptor_variant` HIGH SO:0001574
   - 3 `splice_donor_variant` HIGH SO:0001575
   - 4 `stop_gained` HIGH SO:0001587
   - 5 `frameshift_variant` HIGH SO:0001589
   - 6 `stop_lost` HIGH SO:0001578
   - 7 `start_lost` HIGH SO:0002012
   - 8 `transcript_amplification` HIGH SO:0001889
   - 9 `feature_elongation` HIGH SO:0001907
   - 10 `feature_truncation` HIGH SO:0001906
   - 11 `inframe_insertion` MODERATE SO:0001821
   - 12 `inframe_deletion` MODERATE SO:0001822
   - 13 `missense_variant` MODERATE SO:0001583
   - 14 `protein_altering_variant` MODERATE SO:0001818
   - 16 `splice_region_variant` LOW SO:0001630
   - 19 `incomplete_terminal_codon_variant` LOW SO:0001626
   - 20 `start_retained_variant` LOW SO:0002019
   - 21 `stop_retained_variant` LOW SO:0001567
   - 22 `synonymous_variant` LOW SO:0001819
   - 23 `coding_sequence_variant` MODIFIER SO:0001580
   - 25 `5_prime_UTR_variant` MODIFIER SO:0001623
   - 26 `3_prime_UTR_variant` MODIFIER SO:0001624
   - 27 `non_coding_transcript_exon_variant` MODIFIER SO:0001792
   - 28 `intron_variant` MODIFIER SO:0001627
   - 30 `non_coding_transcript_variant` MODIFIER SO:0001619
   - 32 `upstream_gene_variant` MODIFIER SO:0001631
   - 33 `downstream_gene_variant` MODIFIER SO:0001632
   - 40 `intergenic_variant` MODIFIER SO:0001628
2. **SO term descriptions (verbatim `description` fields):**
   - `missense_variant`: "A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved".
   - `stop_gained`: "A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript".
   - `synonymous_variant`: "A sequence variant where there is no resulting change to the encoded amino acid".
   - `splice_acceptor_variant`: "A splice variant that changes the 2 base region at the 3' end of an intron".
3. **Impact is a property of the consequence term, not computed separately:** the `impact` field is stored on each `OverlapConsequence` (e.g. `frameshift_variant` → HIGH, `missense_variant` → MODERATE, `synonymous_variant` → LOW, `intron_variant` → MODIFIER).

### Ensembl VEP — `Bio::EnsEMBL::Variation::Utils::VariationEffect` (consequence predicates)

**URL:** https://raw.githubusercontent.com/Ensembl/ensembl-variation/release/110/modules/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm
**Accessed:** 2026-06-13 (retrieved with `curl`; HTTP 200)
**Authority rank:** 3 (reference implementation — algorithm logic)

**Key Extracted Points:**

1. **`synonymous_variant` predicate (verbatim):** `return ( ($alt_pep eq $ref_pep) and (not stop_retained(@_) and ($alt_pep !~ /X/) and ($ref_pep !~ /X/)) );` — alt peptide equals ref peptide, neither is unknown (`X`), and it is not a retained stop.
2. **`missense_variant` predicate (verbatim):** returns false if `start_lost`, `stop_lost`, `stop_gained`, or `partial_codon`; otherwise `return ( $ref_pep ne $alt_pep ) && ( length($ref_pep) == length($alt_pep) );` — ref and alt peptides differ and are the same length.
3. **`stop_gained` predicate (verbatim):** `$cache->{stop_gained} = ( ($alt_pep =~ /\*/) and ($ref_pep !~ /\*/) );` — alt peptide contains a stop (`*`) that the ref peptide does not.
4. **`stop_lost` predicate (verbatim, SNV branch):** `$cache->{stop_lost} = ( ($alt_pep !~ /\*/) and ($ref_pep =~ /\*/) );` — ref peptide is a stop, alt peptide is not.
5. **`frameshift` predicate (verbatim, sequence-variant branch):** `my $var_len = $bvfo->cds_end - $bvfo->cds_start + 1; my $allele_len = $bvfoa->seq_length; return abs( $allele_len - $var_len ) % 3;` — a frameshift iff the length difference between alt and ref alleles is not a multiple of 3.
6. **inframe insertion vs deletion:** `inframe_insertion` requires `length($alt_codon) > length($ref_codon)`; `inframe_deletion` requires `length($alt_codon) < length($ref_codon)` (both with the indel length divisible by 3, i.e. not frameshift).
7. **start_lost gating:** `missense_variant` and `inframe_insertion` explicitly `return 0 if start_lost(@_)`, establishing the precedence order in which `start_lost` (a HIGH/rank-7 term) wins over coding substitution terms.

### McLaren W, et al. (2016) — The Ensembl Variant Effect Predictor

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4893825/
**Accessed:** 2026-06-13 (WebFetch; open-access PMC mirror of Genome Biology 17:122)
**Authority rank:** 1 (peer-reviewed publication)

**Key Extracted Points:**

1. **SO terms (verbatim):** "Variant consequences are described using a standardized set of variant annotation terms which were defined in collaboration with the Sequence Ontology (SO)."
2. **Most-specific / hierarchy (verbatim):** "Hierarchy in the predicate system preserves the tree structure of the SO such that only the most specific term that applies under any given parent term is assigned."
3. **Impact classification:** VEP provides an "Impact classification" over its Sequence-Ontology consequence types (Table 1).

### NCBI Genetic Codes — Standard code (transl_table 1)

**URL:** https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 2 (official NCBI standard)

**Key Extracted Points:**

1. **Standard code AAs (verbatim):** `FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG` indexed by Base1/Base2/Base3 = `TTTT…GGGG` / `TTTTCCCC…GGGG` / `TCAGTCAG…`. Stop codons (`*`) are TAA, TAG, TGA; the primary start codon (`M` in the Starts string) is ATG. This is the table used to translate reference and alternate codons to compare peptides (drives synonymous/missense/stop_gained/stop_lost classification).

---

## Documented Corner Cases and Failure Modes

### From VEP `VariationEffect.pm`

1. **Unknown peptide (`X`) excludes synonymous:** the synonymous predicate explicitly excludes any allele whose translation contains `X` (an untranslatable/ambiguous codon), so ambiguous codons are not reported as synonymous.
2. **stop_gained overrides missense:** `missense_variant` returns false when `stop_gained` is true; a substitution that introduces a premature stop is `stop_gained`, never `missense`.
3. **start_lost overrides coding substitution:** `missense_variant`/`inframe_insertion` return false when `start_lost` is true.
4. **Frameshift is length-based:** any coding indel whose `|alt − ref|` length is not divisible by 3 is `frameshift_variant`, regardless of sequence content.

### From NCBI gc.prt

1. **Stop in reference codon:** a substitution at a reference stop codon (TAA/TAG/TGA) that yields a sense codon is `stop_lost`.

---

## Test Datasets

### Dataset: Standard-code codon substitutions (NCBI transl_table 1)

**Source:** NCBI gc.prt, Standard code; VEP `VariationEffect.pm` peptide predicates.

| Variant (codon ref→alt) | Ref AA → Alt AA | Expected consequence | Expected IMPACT |
|--------------------------|------------------|----------------------|------------------|
| GAA → GTA (Glu→Val)      | E → V            | missense_variant     | MODERATE         |
| TTA → TTG (Leu→Leu)      | L → L            | synonymous_variant   | LOW              |
| CAA → TAA (Gln→Stop)     | Q → *            | stop_gained          | HIGH             |
| TAA → CAA (Stop→Gln)     | * → Q            | stop_lost            | HIGH             |
| ATG → ATC (Met→Ile, start)| M → I (at CDS start) | start_lost      | HIGH             |

### Dataset: Coding indels (length rule)

**Source:** VEP `VariationEffect.pm` frameshift / inframe predicates.

| Variant (ref→alt) | Δlen | Expected consequence | Expected IMPACT |
|-------------------|------|----------------------|------------------|
| AC → A            | −1   | frameshift_variant   | HIGH             |
| A → ATTT (≡ +3)   | +3   | inframe_insertion    | MODERATE         |
| ATTT → A (≡ −3)   | −3   | inframe_deletion     | MODERATE         |

---

## Assumptions

1. **ASSUMPTION: Standard genetic code (table 1).** VEP uses transcript-specific codon tables; this unit translates with the NCBI Standard code (table 1) only. Non-standard organism codes are out of scope and would change `*`-codon classification. Justified: human nuclear transcripts (the dominant VEP use case) use table 1; the Standard code is the documented default.
2. **ASSUMPTION: single-codon SNV peptide comparison.** Reference/alternate peptides are compared over the codon(s) directly overlapped by the variant, matching the VEP predicate which compares `ref_pep`/`alt_pep` for the affected codon(s); long-range effects (NMD, downstream re-initiation) are not modelled.

---

## Recommendations for Test Coverage

1. **MUST Test:** missense vs synonymous vs stop_gained vs stop_lost vs start_lost via real codon translation with the exact expected SO term and IMPACT — Evidence: VEP `VariationEffect.pm` predicates + NCBI Standard code.
2. **MUST Test:** frameshift vs inframe insertion/deletion by indel length modulo 3 — Evidence: VEP frameshift/inframe predicates.
3. **MUST Test:** IMPACT mapping for every term equals the `impact` field in `Constants.pm` (HIGH/MODERATE/LOW/MODIFIER) — Evidence: `Constants.pm`.
4. **MUST Test:** most-severe selection across transcripts/consequences returns the lowest-rank term — Evidence: `Constants.pm` ranks + McLaren 2016.
5. **SHOULD Test:** ambiguous codon (`X`) is not reported synonymous — Rationale: explicit VEP exclusion.
6. **COULD Test:** null/empty inputs surface as argument validation — Rationale: library robustness (not specified by VEP).

---

## References

1. McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GRS, Thormann A, Flicek P, Cunningham F (2016). The Ensembl Variant Effect Predictor. *Genome Biology* 17:122. https://doi.org/10.1186/s13059-016-0974-4 (open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC4893825/)
2. Ensembl ensembl-variation (release/110). `Bio::EnsEMBL::Variation::Utils::Constants` — OverlapConsequence impact/rank table. https://raw.githubusercontent.com/Ensembl/ensembl-variation/release/110/modules/Bio/EnsEMBL/Variation/Utils/Constants.pm
3. Ensembl ensembl-variation (release/110). `Bio::EnsEMBL::Variation::Utils::VariationEffect` — consequence predicate logic. https://raw.githubusercontent.com/Ensembl/ensembl-variation/release/110/modules/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm
4. NCBI. The Genetic Codes (gc.prt), Standard code (transl_table 1). https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
5. Eilbeck K, et al. (2005). The Sequence Ontology: a tool for the unification of genome annotations. *Genome Biology* 6:R44. https://doi.org/10.1186/gb-2005-6-5-r44 (SO accessions cited in source 2)

---

## Change History

- **2026-06-13**: Initial documentation (VARIANT-ANNOT-001).
