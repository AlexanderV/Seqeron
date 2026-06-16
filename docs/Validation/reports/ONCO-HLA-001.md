# Validation Report: ONCO-HLA-001 — HLA Allele Nomenclature Parsing & Allele-Specific HLA LOH (LOHHLA)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.ParseHlaAllele(string)`, `OncologyAnalyzer.TryParseHlaAllele(string?, out HlaAllele)` (delegate), `OncologyAnalyzer.DetectHlaLoh(HlaAlleleCopyNumber)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

This unit is distinct from the Phase-1 sibling ONCO-IMMUNE-001 (immune infiltration). It provides
(1) a deterministic validator/parser of the WHO HLA nomenclature grammar and (2) a deterministic
threshold classifier of allele-specific HLA LOH from caller-supplied per-allele copy number and an
allelic-imbalance p value. It deliberately does **not** genotype HLA from reads.

## Stage A — Description

### Sources opened this session (retrieved, not trusted from the repo)

| # | Source | URL | What it confirms |
|---|--------|-----|------------------|
| 1 | WHO/IPD-IMGT/HLA — "Naming Alleles" | https://hla.alleles.org/pages/nomenclature/naming_alleles/ | Name structure, field meanings, two-field minimum, suffix set |
| 2 | IPD-IMGT/HLA — "Assigning Suffixes" | https://hla.alleles.org/pages/nomenclature/assigning_suffixes/ | Exact definition of each suffix N/L/S/C/A/Q; this is the complete suffix set |
| 3 | McGranahan et al. 2017, Cell 171(6):1259–1271 (PMC5720478) — LOHHLA | https://pmc.ncbi.nlm.nih.gov/articles/PMC5720478/ | CN < 0.5 loss threshold; allelic imbalance p < 0.01 paired t-test |

### Formula / grammar check (verbatim from sources)

- **Name structure** (Source 1, verbatim): an allele name is "a unique number corresponding to up to
  **four sets of digits separated by colons**"; gene separated from fields by `*`. → grammar
  `HLA-<Gene>*F1:F2[:F3[:F4]][suffix]` is correct.
- **Field meanings** (Source 1, verbatim): F1 "describe[s] the type, which often corresponds to the
  serological antigen"; F2 lists "the subtypes" (specific protein); F3 "synonymous nucleotide
  substitutions … within the coding sequence"; F4 "polymorphisms in the introns, or in the 5´ or 3´
  untranslated regions". → matches the doc/code field semantics exactly.
- **Two-field minimum** (Source 1, verbatim): "All alleles receive **at least a four digit name,
  which corresponds to the first two sets of digits**." → minimum two fields. ✔
- **Four-field maximum** (Source 1, "up to four sets of digits"). → maximum four. ✔
- **Expression suffixes** (Source 2, verbatim, all six confirmed and stated as the standard set):
  N = "not expressed" (Null); L = "Low cell surface expression"; S = "soluble, Secreted molecule …
  not present on the cell surface"; C = "present in the Cytoplasm and not on the cell surface";
  A = "Aberrant expression where there is some doubt as to whether a protein is actually expressed";
  Q = "questionable expression". → code's enum + `TryMapSuffix`/`SuffixLetter` map is exactly correct.
  (Note: Source 2 records that as of March 2026 no alleles have been *named* with C or A, but both
  are valid, defined suffix letters — accepting them is correct.)
- **LOH loss threshold** (Source 3, verbatim): "**A copy number < 0.5, is classified as subject to
  loss, and thereby indicative of LOH.**" → strict `< 0.5`. ✔
- **Allelic-imbalance guard** (Source 3, verbatim): "**Allelic imbalance is determined if p < 0.01
  using the paired Student's t-Test** between the two distributions." → strict `< 0.01`. ✔

### Edge-case semantics check

- Single field / five fields / missing `HLA-` prefix / bad suffix letter → invalid, all directly
  sourced (two-field min, four-field max, suffix set). ✔
- CN = 0.5 retained; p = 0.01 not significant → strict inequalities are the verbatim source wording. ✔
- Both alleles < 0.5 → labelled `HomozygousLoss`, not allele-specific LOH. The source defines a *lost*
  allele as one with CN < 0.5 and does not address two simultaneously-lost homologs; the repo's
  ASM-01 assumption (homozygous deletion ≠ allele-specific LOH) is biologically sound and explicitly
  flagged as an assumption (only the label changes, the two thresholds are unchanged). Accepted as a
  documented, reasonable divergence — it does not contradict any source.

### Independent cross-check (numbers)

Hand-recomputed all six LOH cases against the validated rule (LOH ⇔ exactly one allele CN < 0.5 ∧
p < 0.01); every expected value traces to the verbatim thresholds, not to code output:

| Case | CN1 | CN2 | p | <0.5? | p<0.01? | Expected |
|------|-----|-----|---|-------|---------|----------|
| M9  | 1.8  | 0.30 | 0.001  | A2 | yes | LOH, lost Allele2 |
| M10 | 0.10 | 1.50 | 0.0005 | A1 | yes | LOH, lost Allele1 |
| M11 | 1.10 | 0.90 | 0.30   | none | no | no LOH |
| M12 | 1.60 | 0.40 | 0.05   | A2 | no  | no LOH (over-call guard) |
| M13 | 1.50 | 0.50 | 0.001  | none (0.50 not < 0.5) | yes | no LOH |
| M14 | 1.70 | 0.40 | 0.01   | A2 | no (0.01 not < 0.01) | no LOH |
| C1  | 0.20 | 0.30 | 0.001  | both | yes | HomozygousLoss (assumption) |

Nomenclature examples independently checked against Source 1 grammar (`HLA-A*02:01` valid 2-field;
`HLA-A*24:02:01:02L` valid 4-field + L; `HLA-A*02` invalid; `A*02:01` invalid; `HLA-A*02:01:01:01:01`
invalid; `HLA-A*02:01X` invalid). All consistent.

### Findings / divergences (Stage A)

None material. ASM-01 (homozygous loss labelling) is a sound, documented assumption. The doc's prose
gloss of suffix C/A ("Cytoplasm", "Aberrant") matches Source 2 verbatim. **Stage A = PASS.**

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `ParseHlaAllele` L6779–6835; suffix map `TryMapSuffix` L6749–6762, `SuffixLetter` L6738–6747; `HlaAllele.Name` L6704–6711.
- `TryParseHlaAllele` L6845–6863.
- `DetectHlaLoh` L6879–6914; constants `HlaLohCopyNumberThreshold = 0.5` L6632, `HlaLohAllelicImbalancePValueThreshold = 0.01` L6641, field bounds `2..4` L6618/6625.

### Formula realised correctly?

- Parser: requires `HLA-` prefix (case-insensitive), splits gene at first `*`, strips a single
  trailing letter only if `char.IsLetter`, validates 2–4 colon fields each a non-empty digit group,
  upper-cases the gene, preserves field digit strings verbatim (leading zeros kept). Matches the
  validated grammar exactly. Strict comparisons in `DetectHlaLoh` (`< 0.5`, `< 0.01`) match the
  verbatim thresholds. Both-lost → `HomozygousLoss` per ASM-01.
- All six suffixes round-trip (verified by new parametrised test): N/L/S/C/A/Q ↔ enum ↔ letter.

### Cross-verification table recomputed vs code

All seven LOH rows above reproduced by running the suite (Failed: 0). Parser examples all behave as
the grammar dictates. No divergence between sourced expectation and code output.

### Variant / delegate consistency

`TryParseHlaAllele` delegates to `ParseHlaAllele` and converts FormatException/ArgumentException to
`false` + default; null → `false` (no throw). Consistent with the canonical (tested S1/S2/S2b).

### Test quality audit (HARD gate)

- **Sourced, not code-echoed:** every expected value is the externally-verified value (grammar
  validity, suffix meaning, 0.5/0.01 thresholds). A deliberately-wrong implementation would fail.
- **No green-washing:** assertions are exact (`Is.EqualTo`, `Throws.TypeOf`), no widened tolerances,
  no skips, no weakened comparisons.
- **Coverage gaps found and fixed this session:** the original 21 tests exercised only **one** of the
  six expression suffixes (`L`) and **none** of the `*`-separator error branches. I added:
  - `ParseHlaAllele_EachExpressionSuffix_MapsAndRoundTrips` — all six suffixes N/L/S/C/A/Q map to the
    correct WHO status and round-trip in `Name` (locks the full suffix map against Source 2).
  - `ParseHlaAllele_LowercaseSuffix_MapsCaseInsensitively` — case-insensitive suffix match.
  - `ParseHlaAllele_MissingOrEmptyGeneFieldSeparator_ThrowsFormatException` — `HLA-A02:01` (no `*`),
    `HLA-*02:01` (empty gene), `HLA-A*` (empty field block) all FormatException.
  Net: +10 test cases; the LOH classifier and validation paths were already fully covered.
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6691` (was 6681); `dotnet build` 0
  errors, 0 new warnings (the 4 build warnings are pre-existing in
  `ApproximateMatcher_EditDistance_Tests.cs`, unrelated to this unit).

**Test-quality gate: PASS** (after adding the missing suffix and separator coverage).

### Findings / defects (Stage B)

No implementation defect. Sole change is additive test coverage (no code change). **Stage B = PASS.**

## Verdict & follow-ups

- **Stage A: PASS** — grammar, field semantics, suffix set, and both LOHHLA thresholds confirmed
  verbatim against primary sources retrieved this session.
- **Stage B: PASS** — code faithfully realises the validated description; all values reproduce.
- **End-state: ✅ CLEAN** — no defect; the test-coverage gap (5 of 6 suffixes and the `*`-separator
  branches untested) was completely fixed this session; full suite green.
- No open defects to log.
