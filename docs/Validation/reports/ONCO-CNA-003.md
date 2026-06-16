# Validation Report: ONCO-CNA-003 — Homozygous (Deep) Deletion Detection

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.DetectHomozygousDeletions(segments, thresholds?, ploidy?)`,
  `OncologyAnalyzer.IsHomozygousDeletion(segment, thresholds?, ploidy?)`,
  `OncologyAnalyzer.IdentifyDeletedTumorSuppressors(deletions)`
  (built on the shared `CallCopyNumber` integer-CN caller and `CopyNumberArmSegment` record)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

> Note on the session prompt: the prompt speculated this "copy-number step 3" might be LST / allele-specific
> CN / genome-instability. It is **not** — the artifacts (Registry row 242, TestSpec, Evidence, algorithm
> doc) consistently define ONCO-CNA-003 as **Homozygous (Deep) Deletion Detection**: a segment is a
> homozygous deletion iff its hard-threshold integer copy number is 0 (total CN 0 = both alleles lost),
> plus an arm→tumour-suppressor-gene panel mapping. Validation proceeded against the *registered*
> definition.

## Stage A — Description

### Sources opened this session (all retrieved live, 2026-06-16)

| Source | URL | What it confirms |
|--------|-----|------------------|
| Cheng et al. 2017, Nat Commun 8:1221 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5663922/ | Homozygous deletion = "regions having **zero copies of both alleles** in the tumour cells" (total CN 0). Distinguished from hemizygous (single-copy) loss where one allele remains. Recurrent TSG targets: CDKN2A (108 events), RB1 (23), PTEN (16), TP53 (4), BRCA2. |
| cBioPortal FAQ | https://docs.cbioportal.org/user-guide/faq/ | Verbatim: "−2 or Deep Deletion indicates a deep loss, possibly a homozygous deletion"; "−1 or Shallow Deletion … possibly a heterozygous deletion"; "0 is diploid"; "1 or Gain"; "2 or Amplification". Calls are putative/unreviewed. |
| CNVkit `cnvlib/call.py` `absolute_threshold` | https://raw.githubusercontent.com/etal/cnvkit/master/cnvlib/call.py | Docstring: "Integer values are assigned for log2 ratio values less than each given threshold value in sequence, counting up from zero." Body uses **`if row.log2 <= thresh:`** (inclusive `<=`); above last threshold `int(np.ceil(...))`; **NaN log2 → neutral copy number** (with a warning). |
| CNVkit pipeline docs | https://cnvkit.readthedocs.io/en/stable/pipeline.html | `call --thresholds`: "log2 ratios up to the first threshold value are assigned a copy number 0 …". Documented **default thresholds = -1.1, -0.4, 0.3, 0.7** (see Note 1). |
| NCBI Gene (×6) | gene/7157, /5925, /1029, /5728, /672, /675 | TP53 **17p13.1**, RB1 **13q14.2**, CDKN2A **9p21.3**, PTEN **10q23.31**, BRCA1 **17q21.31**, BRCA2 **13q13.1**. |

### Formula / definition check

- **Homozygous = total CN 0**: matches Cheng et al. verbatim. ✓
- **CN 0 ⇔ cBioPortal −2 (Deep Deletion)**: matches FAQ verbatim. −1 (shallow/heterozygous) is explicitly **not** homozygous. ✓
- **Integer-CN calling**: CN = index of the first cutoff the log2 value is `<=` (inclusive), counting from 0;
  above the last cutoff `ceil(ploidy·2^log2)`; NaN → neutral (rounded ploidy). The repo's `CallCopyNumber`
  reproduces this exactly (verified against the CNVkit source body, including the `<=` operator). ✓
- **Arm mapping**: all six cytobands confirmed against NCBI Gene; arm letters (p/q) correct; 13q correctly
  carries **both** RB1 and BRCA2. ✓

### Edge-case semantics

- log2 exactly at deletion cutoff (−1.1) → CN 0 (inclusive `<=`, sourced from CNVkit). ✓
- log2 just above → CN 1, not homozygous. ✓
- NaN log2 → neutral no-call (rounded ploidy), not homozygous. ✓ (matches CNVkit "replacing with neutral").
- Single-copy loss (CN 1) excluded. ✓

### Independent cross-check (hand-computed against the verified CNVkit rule, default repo thresholds −1.1,−0.25,0.2,0.7, ploidy 2)

| log2 | CN trace | CN | Homozygous? |
|------|----------|----|-------------|
| −2.0 | ≤ −1.1 | 0 | yes |
| −1.5 | ≤ −1.1 | 0 | yes |
| −1.1 | ≤ −1.1 (boundary) | 0 | yes |
| −1.0999 | > −1.1; ≤ −0.25 | 1 | no |
| −0.5 | > −1.1; ≤ −0.25 | 1 | no |
| 0.0 | ≤ 0.2 | 2 | no |
| 0.5 | ≤ 0.7 | 3 | no |
| 1.0 | > 0.7 → ceil(2·2)=4 | 4 | no |
| NaN | no-call → round(2)=2 | 2 | no |
| NaN (ploidy 3) | no-call → round(3)=3 | 3 | no |
| −0.5 (custom cutoff −0.4) | ≤ −0.4 | 0 | yes |

All of these match the test fixture's expected outcomes.

### Findings / divergences

- **NOTE 1 (non-blocking, carried from ONCO-CNA-001):** the published CNVkit `call --thresholds` default is
  **(-1.1, -0.4, 0.3, 0.7)**; the repo's `DefaultCopyNumberThresholds` is **(-1.1, -0.25, 0.2, 0.7)**. For
  *this* unit only the deletion cutoff governs CN 0, and **-1.1 is identical** in both, so no homozygous-deletion
  result is affected. This is a copy-number-classification matter owned by ONCO-CNA-001, not a defect of
  ONCO-CNA-003. Recorded as a Stage-A note.
- **NOTE 2 (by design):** homozygous status is inferred from **total** integer CN 0, not allele-specific CN,
  so copy-neutral LOH / allele-specific zero with a retained allele is not distinguished. Honestly disclosed in
  the algorithm doc §6.2; consistent with cBioPortal discrete calls. Not a defect.
- **NOTE 3 (interpretation):** cBioPortal hedges "−2 … **possibly** a homozygous deletion" (putative calls).
  The doc/spec faithfully reproduce this caveat; the deterministic CN-0 rule is the correct discrete operationalization.

Stage A is correct in the abstract; the only divergence (Note 1) is non-load-bearing for this unit → **PASS-WITH-NOTES**.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `CallCopyNumber` (lines ~4010-4038): inclusive `<=` loop counting from 0, NaN→neutral, ceil above last — matches CNVkit source verbatim.
- `IsHomozygousDeletion` (lines ~4355-4362): `ValidateArmSegment` then `CallCopyNumber(...) == 0`. ✓
- `DetectHomozygousDeletions` (lines ~4379-4396): null-guard + order-preserving filter on `IsHomozygousDeletion`. ✓
- `IdentifyDeletedTumorSuppressors` (lines ~4408-4432): distinct deleted arms (HashSet, OrdinalIgnoreCase) → panel emitted in fixed order, each gene once. ✓
- `TumorSuppressorArms` panel (lines ~4440-4448): TP53 17p, RB1 13q, CDKN2A 9p, PTEN 10q, BRCA1 17q, BRCA2 13q — matches NCBI Gene; 13q→{RB1,BRCA2} in panel order. ✓

### Formula realised correctly?

Yes. The integer-CN calling, the CN==0 predicate, the order-preserving filter, and the arm→gene mapping all
realise the validated description. Hand-traced cross-check table (above) matches code behaviour and the tests.

### Variant / delegate consistency

`DetectHomozygousDeletions` is defined purely in terms of `IsHomozygousDeletion`, which is defined in terms of
`CallCopyNumber` (the shared ONCO-CNA-001 caller). No divergent reimplementation. ✓

### Numerical robustness

No precision loss in the integer-CN comparison; NaN handled explicitly; ploidy/threshold validation delegated to
`CallCopyNumber`/`ValidateThresholds` (throws on ploidy≤0, non-4/non-ascending thresholds). ✓

### Test quality audit (HARD gate)

Pre-existing fixture: `OncologyAnalyzer_DetectHomozygousDeletions_Tests.cs`, 22 tests. Audited against the
externally-sourced values, not the code:

- **Sourced, not echoes:** every assertion uses exact `Is.EqualTo` / exact counts / ordered sequence equality.
  Expected values trace to the sources above (CN-0 boundary −1.1, mapping cytobands, NaN no-call). No
  `Greater/AtLeast/Contains`, no ranges where an exact value is known, no widened tolerances, no skips. ✓
- **Coverage gaps found (documented behaviour with no test) — FIXED this session:**
  1. **C3** in the TestSpec says invalid = "ArmLength ≤ 0 **or** End ≤ Start", but only `End ≤ Start` was tested.
     Added `DetectHomozygousDeletions_NonPositiveArmLength_Throws`.
  2. Contract §3.3 documents `ArgumentOutOfRangeException` for non-positive ploidy — untested. Added
     `DetectHomozygousDeletions_NonPositivePloidy_Throws`.
  3. Contract §3.3 documents `ArgumentException` for thresholds not four strictly-ascending values — untested.
     Added `DetectHomozygousDeletions_InvalidThresholds_Throws` (non-ascending + wrong-count).
  4. `IsHomozygousDeletion` validation and boundary paths were only happy-path tested. Added
     `IsHomozygousDeletion_DeletionCutoffBoundary_TrueThenFalse` (−1.1 true, −1.0999 false) and
     `IsHomozygousDeletion_InvalidSegment_Throws`.
- No assertion was weakened; no expected value was changed to match output.

Fixture: 22 → **28** tests, all green. **Gate: PASS.**

### Findings / defects

- **No algorithm defect.** Code faithfully realises the validated description.
- **Test-coverage gap (fixed):** 5 documented error/boundary behaviours had no test (above) — added, locked to
  sourced/contract values. Logged as FINDINGS_REGISTER A20.

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** (Note 1 default-thresholds divergence is owned by ONCO-CNA-001 and does not affect
  CN-0 calling; Notes 2–3 are disclosed by-design simplifications/caveats).
- **Stage B: PASS** (no code defect; test-coverage gap fixed).
- **End-state: CLEAN.** `dotnet build` 0 errors; full unfiltered suite **6660 passed, 0 failed** (1 pre-existing
  unrelated skip). This unit's fixture 22→28.
- Follow-up: none required for this unit. The CNVkit default-threshold value (Note 1) remains an ONCO-CNA-001
  concern, already noted there.
