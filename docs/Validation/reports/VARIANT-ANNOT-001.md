# Validation Report: VARIANT-ANNOT-001 — Variant Annotation (functional impact / SO consequences)

- **Validated:** 2026-06-15   **Area:** Variants
- **Canonical method(s):** `VariantAnnotator.PredictFunctionalImpact`, `VariantAnnotator.Annotate`, `VariantAnnotator.GetImpactLevel`, `VariantAnnotator.GetConsequenceRank`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent of repo artifacts)

1. **Ensembl ensembl-variation rel/110 `Utils/Constants.pm`** — OverlapConsequence impact/rank table.
   Retrieved via WebFetch from
   `https://raw.githubusercontent.com/Ensembl/ensembl-variation/release/110/modules/Bio/EnsEMBL/Variation/Utils/Constants.pm`.
2. **Ensembl ensembl-variation rel/110 `Utils/VariationEffect.pm`** — consequence predicates.
   Retrieved via WebFetch from the rel/110 raw URL.
3. **NCBI The Genetic Codes (gc.prt), Standard code (transl_table 1)** — retrieved via WebFetch from
   `https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt`, then decoded the AAs/Starts/Base1-3 strings by
   hand (Python) to recover exact codon→AA and start/stop assignments.
4. McLaren et al. (2016) *The Ensembl Variant Effect Predictor*, Genome Biology 17:122 (most-severe
   reporting + SO terms) — relied on the previously-cited PMC open-access copy.

### Formula / logic check (against retrieved source text)

- **IMPACT ⟷ rank table.** Every rank/impact in the implementation's `ConsequenceRank` dictionary and
  `GetImpactLevel` switch was compared element-by-element against the Constants.pm table fetched this
  session. All shared terms match exactly: transcript_ablation=1, splice_acceptor=2, splice_donor=3,
  stop_gained=4, frameshift=5, stop_lost=6, start_lost=7, inframe_insertion=11, inframe_deletion=12,
  missense=13, splice_region=16, synonymous=22, coding_sequence=23, 5'UTR=25, 3'UTR=26, intron=28,
  upstream=32, downstream=33, intergenic=40 (all HIGH/MODERATE/LOW/MODIFIER as in source).
- **Consequence predicates** (verbatim Perl retrieved from VariationEffect.pm):
  - `synonymous_variant`: `alt_pep eq ref_pep` AND not stop_retained AND neither peptide matches `/X/`. ✓
  - `missense_variant`: returns 0 if start_lost/stop_lost/stop_gained/partial_codon; else `ref ne alt && len==len`. ✓
  - `stop_gained`: `alt_pep =~ /\*/ and ref_pep !~ /\*/`. ✓
  - `stop_lost`: `alt_pep !~ /\*/ and ref_pep =~ /\*/`. ✓
  - `frameshift`: `abs(allele_len − var_len) % 3` (non-zero ⇒ frameshift). ✓
  - `inframe_insertion`/`inframe_deletion`: length(alt_codon) > / < length(ref_codon), divisible by 3, gated by start_lost. ✓
  - `start_lost`: requires `_overlaps_start_codon` and `translation_start == 1` and the alt peptide no
    longer matches the reference start peptide. The implementation models this as
    `codonNumber==1 && IsStartCodon(refCodon) && !IsStartCodon(altCodon)` — a faithful simplification for
    the SNV case.
- **Standard genetic code.** Decoded NCBI table 1 by hand:
  GAA→E, GTA→V, TTA→L, TTG→L, CAA→Q, TAA→`*`, ATG→M, ATC→I. Start codons (M in Starts string) = **TTG, CTG, ATG**;
  stop codons = TAA, TAG, TGA. The implementation's `GeneticCode.Standard` start set `{AUG, UUG, CUG}` and
  stop set `{UAA, UAG, UGA}` match this exactly. (Note: a WebFetch summary erroneously claimed "only ATG";
  the hand-decode of the Starts string disproves that and confirms the code is right.)

### Edge-case semantics

- Ambiguous codon (peptide `X`) excluded from synonymous — sourced from the synonymous predicate's `!~ /X/` guard. ✓
- stop_gained / start_lost precedence over missense — sourced from the missense predicate's leading `return 0 if …` guards. ✓
- Frameshift is length-based regardless of sequence — sourced from the `% 3` predicate. ✓

### Independent cross-check (exact numbers)

Hand computation against the retrieved NCBI table for every coding test variant (CDS starts at genomic 100,
window `ATGCAAGAATTATAANAAGGG`):

| Variant | cdsPos | codon# | ref codon→AA | alt codon→AA | Consequence | Impact |
|---|---|---|---|---|---|---|
| 107 A>T | 8 | 3 | GAA→E | GTA→V | missense | MODERATE |
| 111 A>G | 12 | 4 | TTA→L | TTG→L | synonymous | LOW |
| 103 C>T | 4 | 2 | CAA→Q | TAA→`*` | stop_gained | HIGH |
| 112 T>C | 13 | 5 | TAA→`*` | CAA→Q | stop_lost | HIGH |
| 102 G>C | 3 | 1 | ATG (start) | ATC (not start) | start_lost | HIGH |
| 117 A>G | 18 | 6 | NAA→X | NAG→X | coding_sequence (not synonymous) | MODIFIER |
| 106 AC>A | — | — | Δ−1 | — | frameshift | HIGH |
| 106 A>ATTT | — | — | Δ+3 | — | inframe_insertion | MODERATE |
| 106 ATTT>A | — | — | Δ−3 | — | inframe_deletion | MODERATE |

All match the implementation's output and the sourced predicates/codon table.

### Findings / divergences

None that affect correctness. The `start_lost` model is a documented simplification (Assumption A2) but is
correct for the single-codon SNV scope of this unit and is consistent with VEP's `translation_start == 1`
requirement.

## Stage B — Implementation

### Code path reviewed

- `src/.../Seqeron.Genomics.Annotation/VariantAnnotator.cs`:
  - `PredictFunctionalImpact` (L586) → `DetermineConsequence` (L343) → `DetermineCodingConsequence` (L462)
    + `RefineCodingSubstitution` (L624) for SNV/MNV codon translation.
  - `Annotate` (L701) → most-severe selection by `GetConsequenceRank` (L570) with ordinal transcript tie-break.
  - `GetImpactLevel` (L761), `ConsequenceRank` table (L516), `ExtractCodon` (L684), `CalculateCdsPosition` (L812).
  - `FormatAsVcfInfo` (L1462) — numeric INFO fields formatted with `CultureInfo.InvariantCulture`.
- `src/.../Seqeron.Genomics.Core/GeneticCode.cs` — Standard table verified codon-by-codon against NCBI.

### Formula realised correctly?

Yes. `RefineCodingSubstitution` extracts the affected forward-strand codon, substitutes the alt base(s),
translates ref/alt with `GeneticCode.Standard`, and applies the VEP predicate ordering
(start_lost → stop_gained → stop_lost → X-exclusion → synonymous → missense). Indels are classified purely
by `(alt.Length − ref.Length) % 3` with sign deciding inframe insertion vs deletion — matching the frameshift
predicate. The most-severe selection orders by the exact Constants.pm rank.

### Cross-verification table

The Stage-A hand-computed table was reproduced exactly by the unit tests (all 21 functional-impact tests pass).

### Variant/delegate consistency

`Annotate` re-refines via `RefineAnnotation` → `PredictFunctionalImpact`, so the most-severe path and the
single-transcript path share the same codon logic. `GetConsequenceRank` ↔ `GetImpactLevel` tiers are
consistent (verified by the new ordering test).

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** every expected consequence/impact/AA traces to the NCBI table + VEP
  predicates retrieved this session, not to the implementation output.
- **Defects fixed in this session (test strengthening, no green-washing):**
  1. **S1** previously asserted only `Is.Not.EqualTo(SynonymousVariant)` — a weak assertion where the term is
     deterministic. Strengthened to assert the exact sourced term `CodingSequenceVariant` (Constants.pm rank 23)
     and its MODIFIER impact.
  2. Added `GetConsequenceRank_MatchesConstantsPmRanks` — direct exact-value coverage of the public
     `GetConsequenceRank` method against the Constants.pm ranks (was only exercised implicitly by M13).
  3. Added `GetConsequenceRank_OrderingConsistentWithImpactTiers` — locks INV-5 (HIGH<MODERATE<LOW<MODIFIER rank ordering).
- **Coverage:** all Stage-A MUST/SHOULD/COULD cases (M1–M13, S1–S3, C1–C2) exercised; every public method in
  scope (`PredictFunctionalImpact`, `Annotate`, `GetImpactLevel`, `GetConsequenceRank`, `FormatAsVcfInfo`) is
  now directly tested. Error cases (null variants/annotations, empty reference) covered.
- **Honest green:** full unfiltered suite = **6484 passed, 0 failed, 1 skipped** (`MFE_Benchmark_AllScenarios`,
  a pre-existing benchmark skip). `dotnet build` 0 errors; the 4 build warnings are pre-existing NUnit2007
  warnings in unrelated files (`ApproximateMatcher_*`), none in files touched here.

**Test-quality gate: PASS** (two weak/missing-coverage defects found and fully fixed from sourced values).

### Findings / defects

No implementation defect. Two test-quality defects (weak S1 assertion; missing `GetConsequenceRank` coverage)
found and fixed in-session.

## Verdict & follow-ups

- **Stage A: PASS.** Description and formulas independently confirmed against Constants.pm, VariationEffect.pm,
  and NCBI transl_table 1.
- **Stage B: PASS.** Code faithfully realises the validated predicates and codon table; tests now lock exact
  sourced values across the full public surface.
- **End-state: CLEAN.**
- Out of scope (unchanged, flagged for a future unit): `PredictPathogenicity`, SIFT/PolyPhen, and conservation
  methods contain invented heuristic constants — not part of VARIANT-ANNOT-001.
