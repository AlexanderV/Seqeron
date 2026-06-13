# Test Specification: VARIANT-ANNOT-001

**Test Unit ID:** VARIANT-ANNOT-001
**Area:** Variants
**Algorithm:** Variant Annotation — functional impact / consequence prediction (VEP / Sequence Ontology)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | McLaren W, et al. (2016) The Ensembl Variant Effect Predictor, Genome Biology 17:122 | 1 | https://doi.org/10.1186/s13059-016-0974-4 (https://pmc.ncbi.nlm.nih.gov/articles/PMC4893825/) | 2026-06-13 |
| 2 | Ensembl ensembl-variation rel/110, `Utils/Constants.pm` (OverlapConsequence impact/rank) | 3 | https://raw.githubusercontent.com/Ensembl/ensembl-variation/release/110/modules/Bio/EnsEMBL/Variation/Utils/Constants.pm | 2026-06-13 |
| 3 | Ensembl ensembl-variation rel/110, `Utils/VariationEffect.pm` (consequence predicates) | 3 | https://raw.githubusercontent.com/Ensembl/ensembl-variation/release/110/modules/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm | 2026-06-13 |
| 4 | NCBI The Genetic Codes (gc.prt), Standard code (transl_table 1) | 2 | https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt | 2026-06-13 |

### 1.2 Key Evidence Points

1. Each consequence carries an `impact` (HIGH/MODERATE/LOW/MODIFIER) and an integer `rank`; severity ordering is by rank (1 = most severe) — source 2.
2. `synonymous_variant` iff alt peptide == ref peptide and neither is unknown (`X`); `missense_variant` iff peptides differ, same length, and not start/stop_lost/stop_gained — source 3.
3. `stop_gained` iff alt peptide contains `*` absent from ref; `stop_lost` iff ref peptide is `*` and alt is not — source 3.
4. `frameshift_variant` iff coding indel length difference is not a multiple of 3; otherwise `inframe_insertion` (alt longer) or `inframe_deletion` (alt shorter) — source 3.
5. IMPACT values: frameshift/stop_gained/stop_lost/start_lost/splice_acceptor/splice_donor = HIGH; missense/inframe_insertion/inframe_deletion = MODERATE; synonymous/splice_region = LOW; intron/UTR/up/downstream/intergenic = MODIFIER — source 2.
6. VEP reports the most specific / most severe consequence and uses SO terms — source 1.
7. Codons translated with NCBI Standard code (table 1): `FFLLSSSS…GGGG`; stops TAA/TAG/TGA; start ATG — source 4.

### 1.3 Documented Corner Cases

- Ambiguous/untranslatable codon (peptide `X`) is excluded from `synonymous_variant` (source 3).
- A substitution introducing a premature stop is `stop_gained`, never `missense` (source 3, predicate gating).
- `start_lost` overrides coding substitution terms at the start codon (source 3, predicate gating).

### 1.4 Known Failure Modes / Pitfalls

1. Culture-sensitive numeric formatting (e.g. `0.010` rendered as `0,010` under non-invariant cultures) breaks VCF INFO field output — must format with `CultureInfo.InvariantCulture` (VCF 4.x fields are locale-independent ASCII). Source: VCF spec convention; observed in `FormatAsVcfInfo`.
2. Placeholder peptide values (`refAa="X"`, `altAa="Y"`) produce non-evidence consequences — replaced by real codon translation (source 3, 4).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictFunctionalImpact(Variant, Transcript, string)` | VariantAnnotator | Canonical | Codon-translation consequence + impact + HGVS-style change |
| `Annotate(IEnumerable<Variant>, IEnumerable<Transcript>, string?)` | VariantAnnotator | Canonical | Most-severe annotation per variant |
| `GetImpactLevel(ConsequenceType)` | VariantAnnotator | Internal→promoted | IMPACT mapping per Constants.pm; tested directly |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | IMPACT(consequence) equals the `impact` field in Constants.pm for every term | Yes | Source 2 |
| INV-2 | A coding SNV is synonymous ⇔ ref peptide == alt peptide (Standard code) | Yes | Sources 3,4 |
| INV-3 | A coding indel is frameshift ⇔ \|alt−ref\| length not divisible by 3 | Yes | Source 3 |
| INV-4 | `Annotate` returns the lowest-rank (most severe) consequence per variant | Yes | Sources 1,2 |
| INV-5 | `PredictFunctionalImpact` consequence severity rank ≥ HIGH > MODERATE > LOW > MODIFIER ordering matches Constants.pm | Yes | Source 2 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Missense SNV | GAA(Glu)→GTA(Val) at codon 1 | MissenseVariant, Moderate, aa "p.E1V" | Src 3,4 |
| M2 | Synonymous SNV | TTA(Leu)→TTG(Leu) | SynonymousVariant, Low | Src 3,4 |
| M3 | Stop gained | CAA(Gln)→TAA(Stop) | StopGained, High | Src 3,4 |
| M4 | Stop lost | TAA(Stop)→CAA(Gln) at stop codon | StopLost, High | Src 3,4 |
| M5 | Start lost | ATG→ATC at CDS start codon | StartLost, High | Src 3,4 |
| M6 | Frameshift deletion | ref "AC"→"A" (Δ−1) in coding | FrameshiftVariant, High | Src 3 |
| M7 | Inframe insertion | ref "A"→"ATTT" (Δ+3) in coding | InframeInsertion, Moderate | Src 3 |
| M8 | Inframe deletion | ref "ATTT"→"A" (Δ−3) in coding | InframeDeletion, Moderate | Src 3 |
| M9 | IMPACT mapping HIGH | GetImpactLevel(StopGained/Frameshift/SpliceDonor/SpliceAcceptor/StopLost/StartLost) | High | Src 2 |
| M10 | IMPACT mapping MODERATE | GetImpactLevel(Missense/InframeInsertion/InframeDeletion) | Moderate | Src 2 |
| M11 | IMPACT mapping LOW | GetImpactLevel(Synonymous/SpliceRegion) | Low | Src 2 |
| M12 | IMPACT mapping MODIFIER | GetImpactLevel(Intron/5'UTR/3'UTR/Up/Down/Intergenic) | Modifier | Src 2 |
| M13 | Most-severe selection | Variant overlapping two transcripts: one missense, one synonymous | Annotate returns the missense (lower rank) annotation | Src 1,2 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Ambiguous codon not synonymous | ref codon contains N → peptide X | Not SynonymousVariant | Src 3 exclusion |
| S2 | Intergenic | variant far from any transcript | IntergenicVariant, Modifier | Annotate path |
| S3 | InvariantCulture VCF formatting | FormatAsVcfInfo under de-DE culture | "SIFT=0.010" with a dot | VCF locale-independence |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Null variants in Annotate | Annotate(null, …) | ArgumentNullException | Robustness |
| C2 | Empty reference sequence | PredictFunctionalImpact with empty refSeq | ArgumentException | Robustness |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Existing tests: `tests/Seqeron/Seqeron.Genomics.Tests/VariantAnnotatorTests.cs` — covers `ClassifyVariant`, `NormalizeVariant`, `AnnotateVariant` consequence routing, `GetImpactLevel`, pathogenicity, conservation, regulatory, VCF. None translate codons (consequence determination uses placeholder peptides). The pathogenicity/conservation/SIFT/PolyPhen methods are NOT in this unit's scope and contain invented constants (out of scope; left untouched).
- No `VariantAnnotator_*_Tests.cs` canonical file existed for this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 Missense (real translation) | ❌ Missing | old test asserts term but on placeholder logic |
| M2 Synonymous | ❌ Missing | old code never returns Synonymous |
| M3 Stop gained | ❌ Missing | old code never returns StopGained |
| M4 Stop lost | ❌ Missing | — |
| M5 Start lost | ❌ Missing | — |
| M6 Frameshift | ❌ Missing | covered shape-only on placeholder |
| M7 Inframe insertion | ❌ Missing | — |
| M8 Inframe deletion | ❌ Missing | — |
| M9–M12 IMPACT mapping | ⚠ Weak | old GetImpactLevel tests partial, no MODIFIER UTR rows |
| M13 Most-severe selection | ❌ Missing | — |
| S1 Ambiguous codon | ❌ Missing | — |
| S2 Intergenic | ✅ Covered (old) | re-covered via Annotate |
| S3 InvariantCulture VCF | ❌ Missing | pre-existing failing test is culture bug |
| C1 Null Annotate | ❌ Missing | — |
| C2 Empty refSeq | ❌ Missing | — |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/VariantAnnotator_FunctionalImpact_Tests.cs` — all VARIANT-ANNOT-001 cases (PredictFunctionalImpact, Annotate, GetImpactLevel, InvariantCulture VCF).
- **Remove:** nothing. `VariantAnnotatorTests.cs` stays (covers out-of-scope methods); its weak `GetImpactLevel` tests are superseded by exact rows here but left in place (not in this unit's scope to delete).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `VariantAnnotator_FunctionalImpact_Tests.cs` | Canonical VARIANT-ANNOT-001 | 18 |
| `VariantAnnotatorTests.cs` | Pre-existing, out-of-scope methods | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ | Implemented | ✅ Done |
| 2 | M2 | ❌ | Implemented | ✅ Done |
| 3 | M3 | ❌ | Implemented | ✅ Done |
| 4 | M4 | ❌ | Implemented | ✅ Done |
| 5 | M5 | ❌ | Implemented | ✅ Done |
| 6 | M6 | ❌ | Implemented | ✅ Done |
| 7 | M7 | ❌ | Implemented | ✅ Done |
| 8 | M8 | ❌ | Implemented | ✅ Done |
| 9 | M9 | ⚠ | Rewritten exact | ✅ Done |
| 10 | M10 | ⚠ | Rewritten exact | ✅ Done |
| 11 | M11 | ⚠ | Rewritten exact | ✅ Done |
| 12 | M12 | ⚠ | Rewritten exact (incl. UTR) | ✅ Done |
| 13 | M13 | ❌ | Implemented | ✅ Done |
| 14 | S1 | ❌ | Implemented | ✅ Done |
| 15 | S2 | ✅ | Re-covered via Annotate | ✅ Done |
| 16 | S3 | ❌ | Implemented (culture fix) | ✅ Done |
| 17 | C1 | ❌ | Implemented | ✅ Done |
| 18 | C2 | ❌ | Implemented | ✅ Done |

**Total items:** 18
**✅ Done:** 18 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | PredictFunctionalImpact_MissenseSnv_ReturnsMissenseModerate |
| M2 | ✅ | PredictFunctionalImpact_SynonymousSnv_ReturnsSynonymousLow |
| M3 | ✅ | PredictFunctionalImpact_PrematureStop_ReturnsStopGainedHigh |
| M4 | ✅ | PredictFunctionalImpact_StopLoss_ReturnsStopLostHigh |
| M5 | ✅ | PredictFunctionalImpact_StartCodonChange_ReturnsStartLostHigh |
| M6 | ✅ | PredictFunctionalImpact_CodingDeletion_ReturnsFrameshiftHigh |
| M7 | ✅ | PredictFunctionalImpact_InframeInsertion_ReturnsInframeInsertionModerate |
| M8 | ✅ | PredictFunctionalImpact_InframeDeletion_ReturnsInframeDeletionModerate |
| M9 | ✅ | GetImpactLevel_HighTerms_ReturnsHigh |
| M10 | ✅ | GetImpactLevel_ModerateTerms_ReturnsModerate |
| M11 | ✅ | GetImpactLevel_LowTerms_ReturnsLow |
| M12 | ✅ | GetImpactLevel_ModifierTerms_ReturnsModifier |
| M13 | ✅ | Annotate_TwoTranscripts_ReturnsMostSevere |
| S1 | ✅ | PredictFunctionalImpact_AmbiguousCodon_NotSynonymous |
| S2 | ✅ | Annotate_VariantFarFromTranscripts_ReturnsIntergenic |
| S3 | ✅ | FormatAsVcfInfo_NonInvariantCulture_UsesDotDecimal |
| C1 | ✅ | Annotate_NullVariants_Throws |
| C2 | ✅ | PredictFunctionalImpact_EmptyReferenceSequence_Throws |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | Standard genetic code (NCBI table 1) for codon translation | M1–M5, S1 |
| A2 | Peptide comparison over the codon(s) directly overlapped by the variant | M1–M5 |

---

## 7. Open Questions / Decisions

1. The out-of-scope `PredictPathogenicity`, SIFT/PolyPhen, and conservation methods in `VariantAnnotator` contain invented constants. They are NOT part of VARIANT-ANNOT-001 (canonical = Annotate + functional impact) and are left untouched; they should be addressed under a dedicated future Test Unit.
