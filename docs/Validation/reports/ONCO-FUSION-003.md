# Validation Report: ONCO-FUSION-003 — Fusion Breakpoint Analysis (reading-frame consequence + fusion protein prediction)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.AnalyzeBreakpoint(FusionBreakpoint)`, `OncologyAnalyzer.PredictFusionProtein(FusionBreakpoint, (string,string))`; reuses `IsInFrame` (ONCO-FUSION-001).
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** ✅ CLEAN (defects found were fully fixed this session)

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)
- **Arriba output spec** (https://github.com/suhrig/arriba/wiki/05-Output-files, WebFetch 2026-06-16): `reading_frame ∈ {in-frame, out-of-frame, stop-codon, .}`, and it "states whether the gene at the 3' end of the fusion is fused in-frame or out-of-frame" — i.e. the call is about the **3' gene's reading frame**. `site1/site2 ∈ {5'UTR, 3'UTR, UTR, CDS, exon, intron, intergenic}`. Confirms the repo's `BreakpointSite` enum, `BreakpointFrameStatus` enum, and `NotPredicted` (= `.`) semantics.
- **AGFusion model.py** (https://raw.githubusercontent.com/murphycj/AGFusion/master/agfusion/model.py, WebFetch 2026-06-16): verbatim — `cds_5prime = transcript1.coding_sequence[0:junction5]`; `cds_3prime = transcript2.coding_sequence[junction3:]`; `seq = cds_5prime + cds_3prime`; frame rule `if len5%3==0 and len3%3==0: in-frame elif round(frac5+frac3,2)==1.0: in-frame (with mutation) else: out-of-frame`; out-of-frame trims `cds[0:3*(len//3)]` before `translate()`; both paths truncate at `protein_seq[0:find("*")]`. Confirmed that `junction3` is the **coding-base offset** where the 3' suffix begins (computed by `_fetch_transcript_cds`), and the frame test uses `len(cds_3prime)%3`.
- **Codon table** (https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables, WebFetch 2026-06-16): ATG→M, AAA→K, GAT→D, GGT→G, AAG→K; TAA/TAG/TGA are stop. All hand-computed peptides confirmed.
- **In-frame definition** (search + Wikipedia "Reading frame"): in-frame fusion = downstream ORF joined so the reading frame is preserved and the 3' protein domains (e.g. kinase domain) are produced intact.

### Formula check
- Site classification (Arriba): correct.
- Frame rule: repo uses `(fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0` (Genomics England exon-phase rule, ONCO-FUSION-001). This is the **Arriba two-way** model — "is the 3' gene read in its native frame" — which is the correct functional notion of in-frame. Verified hand-derivation: 5' contributes `b` bases (chimeric phase `b mod 3`); 3' suffix's first base is at native phase `p`; native-frame reading requires `b ≡ p (mod 3)`. Matches `IsInFrame`.
- Protein prediction: chimeric = 5' prefix ++ 3' suffix, translate, truncate at first stop, trim to whole codons — matches AGFusion exactly.

### Independent cross-check (numbers)
Python recomputation of the worked vectors vs both AGFusion's three-way rule and the repo's `IsInFrame`:

| Vector | chimeric | AGFusion(len5,len3) | repo IsInFrame(b,p) | peptide |
|--------|----------|---------------------|---------------------|---------|
| M8 `ATGAAA`+`GATGGT` (b6,p0) | `ATGAAAGATGGT` | in-frame | InFrame | MKDG |
| M9 `ATGAAA`+`GATTAAGGT` (b6,p0) | `ATGAAAGATTAAGGT` | in-frame | InFrame | MKD (stop at TAA) |
| M10 `ATGA`+`AAGGT` (b4,p0) | `ATGAAAGGT` | **in-frame (with mutation)** | **OutOfFrame** | MKG |
| M11 `ATGAA`+`GATGGT` (b5,p0) | `ATGAAGATGGT` | out-of-frame | OutOfFrame | MKM (trim to 9) |

### Findings / divergences (Stage A)
- **DEFECT-A1 (fixed):** the Evidence dataset table (line 107) and TestSpec M10 row labelled the M10 vector (`ATGA`+`AAGGT`, phase 0) as **in-frame**, citing AGFusion's "in-frame (with mutation)" complement rule (`4%3=1, 5%3=2; 1+2=3`). But the implementation models the **Arriba two-way** in/out call (3' gene native frame), under which `(4−0)%3=1 ≠ 0` → **out-of-frame**, which is the biologically correct functional call (the 3' gene is read frameshifted, so its domains are not preserved). The Evidence conflated AGFusion's *three-way ORF-continuity* class with the repo's *two-way 3'-gene-frame* model. Corrected the Evidence table, added a MODEL note explaining the Arriba-vs-AGFusion distinction, and fixed the TestSpec M10 row (+ added M10b for the genuine in-frame mid-codon case).

The underlying biology/maths is correct; only the M10 narrative label was wrong → PASS-WITH-NOTES.

## Stage B — Implementation

### Code path reviewed
`OncologyAnalyzer.cs` lines 3464–3480 (`IsInFrame`), 3785–3804 (`AnalyzeBreakpoint`), 3832–3890 (`PredictFusionProtein`), 3676–3766 (enums/records).

### Formula realised correctly?
- `AnalyzeBreakpoint`: frame call only when both sites `Cds`, else `NotPredicted`; otherwise `IsInFrame` → In/OutOfFrame. Correct (Arriba `.`).
- `PredictFusionProtein`: `cds5 = fivePrimeCds[0:junction5]`, `cds3 = threePrimeCds[junction3:]`, concat, translate via validated `GeneticCode.Standard`, break at first stop (`HasPrematureStop`), trim to whole codons (`translatableLength = len − len%3`). Matches AGFusion. All hand-computed peptides (MKDG, MKD, MKG, MKM, MK) reproduced.
- Validation: negative/over-range offsets → `ArgumentOutOfRangeException`; null CDS → `ArgumentNullException`; phase outside {0,1,2} → `ArgumentOutOfRangeException` (via `IsInFrame`).

### Variant/delegate consistency
`AnalyzeBreakpoint` and `PredictFusionProtein` both route the frame call through `IsInFrame`; consistent.

### Test quality audit (HARD gate)
- **DEFECT-B1 (fixed):** test `PredictFusionProtein_MidCodonJunction_InFramePeptide` (M10) was **green-washed** — its name and comment claimed "in-frame (1+2 complement)" but it deliberately **omitted the `Effect` assertion**, the one load-bearing claim, precisely because asserting `InFrame` would fail (impl returns `OutOfFrame`) and asserting `OutOfFrame` would contradict its own narrative. A test that dodges the disputed assertion to stay green is a defect. **Fix:** renamed to `PredictFusionProtein_MidCodonJunctionPhaseMismatch_OutOfFrame`, now asserts `Effect == OutOfFrame` with sourced reasoning; added a new honest in-frame mid-codon case `PredictFusionProtein_MidCodonJunctionPhaseMatch_InFrame` (b4, p1, suffix `AAGGT` from `TAAGGT`, `(4−1)%3=0` → InFrame, MKG) so the in-frame-with-non-zero-phase protein path is covered without weakening anything.
- Other tests assert exact sourced values (peptides, enum values, exceptions); cover both methods, all four `BreakpointFrameStatus` outcomes (InFrame/OutOfFrame/NotPredicted via AnalyzeBreakpoint; premature-stop flag via PredictFusionProtein), phases 0/1/2, premature stop, empty 3' suffix, null/out-of-range/invalid-phase. No tolerances, no skips, no inequality-where-exact.
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6649`; `dotnet build` 0 errors (4 pre-existing warnings in unrelated files, none in changed files).

### Findings / defects
Both defects (A1 narrative, B1 green-washed test) were completely fixed this session. The algorithm code was already correct (Arriba model); no code change was needed — the fixes were to the wrong narrative and the dodging test.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (M10 label corrected; Arriba-vs-AGFusion model documented).
- **Stage B:** PASS-WITH-NOTES (green-washed M10 test fixed; in-frame mid-codon coverage added).
- **End-state:** ✅ CLEAN — algorithm fully functional; all defects fixed; full suite green.
- **Test-quality gate:** PASS (after fix). One green-washed test found and corrected; coverage extended; full unfiltered suite green.
