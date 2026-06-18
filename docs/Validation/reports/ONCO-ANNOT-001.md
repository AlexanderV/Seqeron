# Validation Report: ONCO-ANNOT-001 — Cancer-Specific Variant Annotation (AMP/ASCO/CAP 2017 four-tier classification)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.ClassifyVariantTier`, `OncologyAnalyzer.AnnotateCancerVariants`, `OncologyAnalyzer.GetCOSMICAnnotation`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session
1. **Li MM et al. (2017), "Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer" (AMP/ASCO/CAP), J Mol Diagn 19(1):4–23, DOI 10.1016/j.jmoldx.2016.10.002.** Primary consensus guideline. Full text retrieved this session two ways:
   - PMC HTML: https://pmc.ncbi.nlm.nih.gov/articles/PMC5707196/ (WebFetch — confirmed tier↔level mapping, 1% cutoff, Tier III/IV distinction).
   - PDF (MCW-hosted authoritative copy): https://ocpe.mcw.edu/sites/default/files/course/2024-03/AMP-ASCO-CAP%20guidelines%20-%20somatic%20variants.pdf — downloaded and converted with `pdftotext -layout`; Figure 2, Tables 6/7, and the Population Databases section read verbatim.
2. **Tate JG et al. (2019), COSMIC, Nucleic Acids Res 47(D1):D941–D947** — confirms COSMIC is an external curated database (millions of mutations), not reproducible in-library; a lookup must be against caller-supplied records.
3. **COSMIC search for BRAF V600E** (WebSearch + sanger.ac.uk; the DB itself is login-gated). COSV56056643 is the established COSMIC genomic-mutation ID for BRAF V600E.

### Formula / decision-rule check (against the retrieved primary source)
The unit is a qualitative classification, not a numeric formula. Verbatim from the retrieved sources:

- **Abstract (l.42–44):** "tier I, variants with strong clinical significance; tier II, variants with potential clinical significance; tier III, variants of unknown clinical significance; and tier IV, variants deemed benign or likely benign." ✔ matches the four `VariantTier` values.
- **Body (l.467–468):** "...(level A and B evidence); tier III, variants with unknown clinical significance; and tier IV, variants that are benign or likely benign (Figure 2)." Combined with Figure 2 columns showing **Level A / Level B under "Strong Clinical Significance" (Tier I)** and **Level C / Level D under "Potential Clinical Significance" (Tier II)**. ✔ matches INV-2 / the code's first two branches exactly.
- **Population cutoff (l.296–298, verbatim):** "In the absence of paired normal tissue, the work group recommends using 1% (0.01) as a primary cutoff." ✔ matches `BenignPopulationMafThreshold = 0.01`.
- **Table 7 (Tier IV), Population database row (l.718):** "MAF [≥] 1% in the general population; or high MAF in some ethnic populations" (the "≥" glyph is lost in extraction but the 1% primary cutoff is the benign threshold and is **inclusive**). ✔ matches the `>=` comparison and the M9 boundary expectation.
- **Figure 2 Tier III box (l.489–504, verbatim):** "Not observed at a significant allele frequency in the general or specific subpopulation databases..." AND "No convincing published evidence of cancer association."
- **Figure 2 Tier IV box (l.491–504, verbatim):** "Observed at significant allele frequency in the general or specific subpopulation databases" / "No existing published evidence of cancer association."
- **Tier III narrative (l.751–758):** Tier III "may include somatic variants in cancer genes reported in the same or different cancer types with unknown clinical significance and variants in cancer genes that have not been reported in any cancers... These variants should not have been observed at significant allele frequencies in the general population."

### Edge-case semantics
- MAF exactly 0.01 ⇒ benign (Tier IV): inclusive cutoff is sourced (1% primary cutoff; Table 7 "≥ 1%"). ✔
- MAF just below cutoff + cancer association ⇒ Tier III. ✔ (Tier III is rare, in cancer genes, no clinical evidence)
- Rare + no cancer association ⇒ Tier IV. Sourced from the Figure 2 Tier IV box "No existing published evidence of cancer association." ✔
- Invalid MAF (NaN, <0, >1) and null inputs throwing are API-contract decisions, not guideline content — reasonable.

### Independent cross-check (numbers / verbatim text)
The Figure 2 box text and the 1% cutoff sentence were extracted verbatim from the PDF (lines quoted above), independent of the repo's Evidence doc. The tier↔level mapping was independently confirmed by the PMC HTML fetch and by a third-party (jmdjournal/CAP) description. All match the description.

### Findings / divergences (the NOTES on Stage A)
**Tier III vs Tier IV discriminator is a documented interpretation, not a verbatim algorithm.** Read literally, **both** Figure 2 boxes (III and IV) say there is *no* convincing published cancer association; the boxes differ on population frequency (III = "not observed at significant frequency", IV = "observed at significant frequency"). The guideline gives no machine-decidable rule; it is a qualitative narrative ("cancer genes...with unknown clinical significance"). The implementation operationalizes this as: *no evidence level + (MAF ≥ 0.01 **OR** no cancer association) ⇒ Tier IV; else Tier III* — i.e. it uses a `HasCancerAssociation` flag (cancer-gene/somatic-DB context) to lift a rare variant out of "benign" into "unknown." This is a **defensible, explicitly-documented reading** (Evidence ASSUMPTION 2; spec §5.4 #2) that correctly captures the population-frequency axis and the spirit of Tier III (rare, cancer-gene context, undetermined significance) vs Tier IV (common, or nothing pointing to cancer). It is not the only possible mapping, hence **PASS-WITH-NOTES** rather than PASS. No biological/mathematical error.

The single numeric constant (1% / 0.01, inclusive) is verbatim-correct. The level→tier mapping is verbatim-correct.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `ClassifyVariantTier` (l.1383–1413): validation (NaN/<0/>1 ⇒ `ArgumentOutOfRangeException`); A/B⇒Tier I; C/D⇒Tier II; else `MAF >= 0.01 || !HasCancerAssociation` ⇒ Tier IV; else Tier III.
- `AnnotateCancerVariants` (l.1424–1436): null-guard; one annotation per input, input order preserved.
- `GetCOSMICAnnotation` (l.1455–1461): null-guard; `TryGetValue((Gene,ProteinChange))`, returns value or null.
- `BenignPopulationMafThreshold = 0.01` (l.1267).

### Formula realised correctly?
Yes. The cascade is exactly the validated decision rule. The `>=` comparison realises the inclusive 1% cutoff; evidence-level branches precede the frequency rule, so a Level A/B biomarker stays Tier I even at high MAF (guideline categorizes by evidence level first). The COSMIC lookup is an exact ordinal-tuple dictionary lookup — does not fabricate annotations (consistent with COSMIC being external, Tate 2019).

### Cross-verification table recomputed vs code
| Input (level, MAF, assoc) | Sourced expected tier | Code output | Source |
|---|---|---|---|
| A, 0.0, true | Tier I | Tier I | Fig 2 (Level A) |
| B, 0.0, true | Tier I | Tier I | Fig 2 (Level B) |
| C, 0.0, true | Tier II | Tier II | Fig 2 (Level C) |
| D, 0.0, true | Tier II | Tier II | Fig 2 (Level D) |
| A, 0.30, false | Tier I | Tier I | Fig 2 (evidence-first) |
| None, 0.25, true | Tier IV | Tier IV | Table 7 (MAF ≥ 1%) |
| None, 0.0001, false | Tier IV | Tier IV | Fig 2 Tier IV box (no cancer assoc.) |
| None, 0.0001, true | Tier III | Tier III | Fig 2 Tier III box / Table 6 |
| None, 0.01, true | Tier IV | Tier IV | 1% inclusive cutoff |
| None, 0.0099, true | Tier III | Tier III | below cutoff, has assoc. |
| (BRAF,p.V600E) in catalog | "COSV56056643" | "COSV56056643" | caller catalog round-trip |
| (TP53,p.R175H) not in catalog | null | null | external lookup |

All twelve match.

### Variant/delegate consistency
`AnnotateCancerVariants` delegates per-variant to `ClassifyVariantTier`; the M11 batch (A,C,common,rare) returns [I, II, IV, III] in order with count 4 — consistent with the canonical.

### Test quality audit (HARD gate)
File: `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AnnotateCancerVariants_Tests.cs` (20 tests).
- **Exact sourced values, no code-echo:** every M-test asserts an exact `VariantTier` enum value (`Is.EqualTo`), never Greater/Contains/range. The values trace to Figure 2 / Tables 6–7 quoted above, not to code output.
- **Discriminating inputs (not tautologies):**
  - M8 (Level A + MAF 0.30 + no assoc.) would be Tier IV without the evidence-first override — genuinely kills a wrong "frequency-first" implementation.
  - M5 (MAF 0.25 + assoc.=true) isolates *high MAF* as the sole Tier-IV driver (association present, so it can't be the no-association branch).
  - M9 (MAF exactly 0.01) kills a strict-`>` off-by-one; M10 (0.0099) is the matching below-boundary case.
  - M6 (rare + no assoc.) vs M7 (rare + assoc.) isolate the `HasCancerAssociation` branch.
- **All branches/edge cases covered:** Tier I (A,B), Tier II (C,D), Tier IV (high MAF; rare+no assoc.), Tier III (rare+assoc.), boundary 0.01/0.0099, batch order/count, COSMIC hit/miss, NaN/<0/>1, null variants, null catalog, empty batch, totality property (C1).
- **No green-washing:** no skipped/ignored/commented-out tests; no weakened assertions; no widened tolerances. C1 asserts only `Enum.IsDefined` but that is its stated purpose (INV-1 totality) and the exact values are pinned by M1–M10.
- **Minor coverage note (not a defect):** the batch test M11 exercises only Levels A and C; Levels B and D in a batch are covered transitively (`AnnotateCancerVariants` is a thin delegate) and directly by M2/M4 on `ClassifyVariantTier`. No fix required.

### Honest-green result
- `dotnet build` (Debug): **0 errors**, 4 warnings — all pre-existing NUnit2007 in the unrelated `ApproximateMatcher_EditDistance_Tests.cs`; this unit's file is warning-free.
- Full unfiltered `dotnet test` (Debug, --no-build): **Failed: 0, Passed: 6628, Skipped: 0.**

### Findings / defects
None. Code faithfully realises the validated description; tests are sourced and discriminating.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES — level→tier mapping and the 1% inclusive cutoff are verbatim-correct; the Tier III/IV discriminator (`HasCancerAssociation`) is a documented, defensible operationalization of a qualitative guideline, not a verbatim rule (the literal Figure 2 boxes distinguish III/IV by population frequency, with both saying "no convincing cancer association").
- **Stage B:** PASS — implementation matches; tests pin exact sourced values and cover all branches/edges.
- **Test-quality gate:** PASS.
- **End-state:** ✅ CLEAN. No code or test change required; algorithm fully functional.
