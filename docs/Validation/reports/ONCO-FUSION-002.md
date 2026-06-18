# Validation Report: ONCO-FUSION-002 — Known Fusion Database Lookup (HGNC designation + directional match)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.GetFusionAnnotation(gene5p, gene3p)`, `OncologyAnalyzer.MatchKnownFusions(FusionCall, IReadOnlyDictionary<string,string>)` (+ `FusionDesignationSeparator`, `KnownFusionMatch`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **Bruford et al. (2021), HGNC recommendations for the designation of gene fusions, *Leukemia* 35(11):3040–3043** — WebFetched `https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/`. Verbatim:
   - Separator: "a double colon (::)—be used in describing gene fusions, e.g. *BCR*::*ABL1*".
   - 5′-first rule: "the 5′ partner gene should always be listed first in the description of a fusion gene, i.e., before the double colon, irrespective of chromosomal location or the orientation of the gene".
   - BCR::ABL1 worked example: "in the *BCR*::*ABL1* fusion gene—the outcome of the translocation t(9;22)(q34.1;q11.2)—the *BCR* gene in chromosome 22 is the 5′ gene, the *ABL1* gene from chromosome 9 is the 3′ gene".
   - Read-throughs: "The HGNC has approved the use of the hyphen separator…for denoting readthrough transcripts, e.g. *INS-IGF2*".
2. **Independent cross-check — Wikipedia "Philadelphia chromosome"** (WebFetched): confirms BCR on chr22, ABL1 on chr9, and the standard designation **BCR::ABL1** (BCR = 5′ upstream, ABL1 = 3′ downstream).
3. **Independent cross-check — EML4-ALK** (WebSearch → Soda et al. 2007, *Nature* `nature05945`, corroborated by Mano 2008 and others): a chr2p inversion "juxtaposes the 5′ end of the EML4 gene with the 3′ end of the ALK gene" → **EML4 = 5′, ALK = 3′**, so the tests' `EML4::ALK` keying is biologically correct.

### Formula check
- Designation = `gene5p + "::" + gene3p` — matches the verbatim 5′-first + double-colon rule (Source 1). ✓
- Directionality (A::B ≠ B::A for A≠B) is a direct logical consequence of "5′ partner always listed first…irrespective of orientation" (the order is fixed by biological role, not alphabetical). ✓
- Read-throughs (hyphen) are explicitly out of scope; the unit emits `::` for true fusions only — consistent with Source 1. ✓

### Edge-case semantics
- Reciprocal-only set → no match (directional). Sourced from the 5′-first rule. ✓
- Null/empty/whitespace partner → `ArgumentException`; null map → `ArgumentNullException` — input-validation contract (standard, not a biology claim). ✓
- Case-insensitive matching (INV-4) is explicitly flagged as an **ASSUMPTION** in Evidence/TestSpec — Bruford et al. specify HGNC approved (case-defined) symbols but do not address case folding. The assumption is reasonable, honestly disclosed, and affects matching only (designation string is preserved verbatim). ✓

### Findings / divergences
None. Description, Evidence, and TestSpec quote the primary source faithfully; every load-bearing biological fact (BCR/ABL1 and EML4/ALK 5′/3′ assignment) was independently re-confirmed from a second source this session.

## Stage B — Implementation

### Code path reviewed
- `OncologyAnalyzer.cs:3580` `FusionDesignationSeparator = "::"`.
- `OncologyAnalyzer.cs:3603–3616` `GetFusionAnnotation`: validates both symbols with `IsNullOrWhiteSpace` (throws `ArgumentException`), returns `gene5p + "::" + gene3p` (verbatim concat → preserves case).
- `OncologyAnalyzer.cs:3640–3665` `MatchKnownFusions`: `ArgumentNullException.ThrowIfNull(knownFusions)`; builds the directional key via `GetFusionAnnotation` (so it inherits the partner validation); `TryGetValue` fast path (honours the dict's own comparer), then a case-insensitive (`OrdinalIgnoreCase`) fallback scan; returns `KnownFusionMatch(designation, IsKnown, annotation)`.

### Formula realised correctly?
Yes. The designation is exactly `5′::3′`; the lookup is exactly directional (the key is the queried fusion's own designation, never the reciprocal). The fallback scan only ever matches the same directional string case-insensitively — it never relaxes direction.

### Cross-verification table (recomputed vs code, values traced to sources)

| Input | Expected (sourced) | Source |
|-------|--------------------|--------|
| `GetFusionAnnotation("BCR","ABL1")` | `"BCR::ABL1"` | Bruford 2021 worked example |
| `GetFusionAnnotation("ABL1","BCR")` | `"ABL1::BCR"` (≠ forward) | 5′-first rule (directional) |
| `GetFusionAnnotation("bcr","abl1")` | `"bcr::abl1"` (case preserved) | verbatim concat |
| Match `EML4→ALK` vs `{EML4::ALK}` | IsKnown=true, annotation returned | Soda 2007 (EML4 5′, ALK 3′) + Bruford keying |
| Match `EML4→ALK` vs `{ALK::EML4}` only | IsKnown=false | directionality |
| Match `eml4→alk` vs `{EML4::ALK}` (case-sensitive dict) | IsKnown=true (fallback scan) | INV-4 assumption |

All verified green by running the suite.

### Variant/delegate consistency
`MatchKnownFusions` builds its key through `GetFusionAnnotation`, so the two methods cannot diverge on format or validation. Fast-path (`OrdinalIgnoreCase` dict, C1) and fallback-scan (default `Ordinal` dict, S1/M3) both reach the same result.

### Test quality audit (TEST-QUALITY GATE)
- **Sourced, not code-echoes:** designation assertions use exact `Is.EqualTo("BCR::ABL1")` etc., traceable to the Bruford worked example, not to code output. The directional negative (M4) and case-insensitive (S1) tests would fail a deliberately-wrong (unordered/case-sensitive) implementation. ✓
- **No green-washing:** no weakened assertions, no widened tolerances, no skips. ✓
- **Coverage gap found and fixed (Stage B defect per the gate):** the validation contract ("null/empty/whitespace partner throw") was only half-covered — M6 tested null-5′ only, M7 empty-3′ only. Added 4 tests with no code change:
  - `GetFusionAnnotation_EmptyFivePrime_Throws`
  - `GetFusionAnnotation_NullThreePrime_Throws`
  - `GetFusionAnnotation_WhitespaceSymbols_Throw` (space on 5′, tab on 3′)
  - `MatchKnownFusions_EmptyPartner_Throws` (propagation of `ArgumentException` through the key build)
- **Honest green:** full unfiltered `dotnet test` = **6648 passed, 0 failed** (was 6644; +4). `dotnet build` = 0 errors. The 4 build warnings are pre-existing NUnit2007 warnings in the unrelated `ApproximateMatcher_EditDistance_Tests.cs`; the changed test file (`OncologyAnalyzer_MatchKnownFusions_Tests.cs`) compiles warning-free.

### Findings / defects
No algorithm or description defect. One test-coverage gap (partial validation-branch coverage), fully fixed this session (logged as FINDINGS_REGISTER A48).

## Verdict & follow-ups
- **Stage A: PASS** — designation format, 5′-first directionality, and the BCR::ABL1 / EML4::ALK assignments independently re-confirmed against Bruford et al. 2021 + Wikipedia + Soda et al. 2007.
- **Stage B: PASS** — code realises the validated spec exactly; the directional lookup and case-insensitive fallback are correct; tests strengthened to cover all documented validation branches.
- **End-state: CLEAN** — fully functional; full suite green (6648/0); working tree committed in one `validate(ONCO-FUSION-002)` commit.
- **Test-quality gate: PASS** (after adding 4 missing edge/error-case tests; no green-washing).
