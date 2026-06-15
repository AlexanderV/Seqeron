# Validation Report: TRANS-SPLICE-001 — Alternative Splicing (Event Classification + Percent Spliced In)

- **Validated:** 2026-06-15   **Area:** Transcriptome
- **Canonical method(s):** `TranscriptomeAnalyzer.CalculatePSI(inclusionReads, exclusionReads, inclusionEffectiveLength?, exclusionEffectiveLength?)`; `TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms)` (+ private `ClassifyIsoformPair`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** FAIL (defect found and **fixed** in-session)
- **End-state:** ✅ CLEAN

---

## Stage A — Description

### Sources opened this session (live)
- **rMATS-turbo README** (reference implementation) — `https://github.com/Xinglab/rmats-turbo` (raw README): A5SS/A3SS output columns `longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE`. **On the + strand, for A5SS the long/short exon forms differ at their downstream (END/3′) boundary while sharing the upstream (START/5′) boundary; for A3SS they differ at their upstream (START/5′) boundary while sharing the downstream (END/3′) boundary.**
- **NAR 34(21):6305** — `https://academic.oup.com/nar/article/34/21/6305/3100610`: alternative 5′/3′SS = "exon extension/truncation"; "donor/acceptor splice site, i.e. 5′/3′SS"; 5′ donor = downstream boundary of the upstream exon, 3′ acceptor = upstream boundary of the downstream exon.
- **Molecular-biology splice-site definitions** (WebSearch, multiple concurring sources): 5′ splice site = **donor**, at the exon→intron boundary (GT), i.e. the **3′ END of the upstream exon**; 3′ splice site = **acceptor**, at the intron→exon boundary (AG), i.e. the **5′ START of the downstream exon**.
- **SUPPA2 / biostars / Outrigger** (WebSearch): PSI = inclusion reads / (inclusion + skipping reads).
- Carried from Evidence doc (Phase-1 fetches): PMC3330053 (Ψ = γᵢ/(γᵢ+γₑ)), Shen 2014 rMATS PMC4280593 (ψ̂ = (I/lᵢ)/(I/lᵢ + S/lₛ)), Wang 2008 five-class taxonomy.

### Formula check
- **Unnormalized PSI** Ψ = I/(I+S) — confirmed verbatim by PMC3330053 and SUPPA2 definition. ✅
- **rMATS length-normalized PSI** ψ̂ = (I/lᵢ)/(I/lᵢ + S/lₛ) — confirmed verbatim by Shen 2014 (PMC4280593). ✅
- **Five-class taxonomy** SE, RI, A5SS, A3SS, MXE — confirmed (Wang 2008; rMATS event vocabulary). ✅

### Edge-case semantics
- 0/0 → undefined → NaN (PMC3330053 pseudo-count rationale). ✅
- S=0,I>0 → 1; I=0,S>0 → 0 (direct from Ψ=I/(I+S)). ✅
- < 2 isoforms per gene → no event (Wang 2008). ✅
- Identical isoforms → no event. ✅ (no structural difference)

### Independent cross-check (numbers)
- M1: 80/(80+20) = **0.80**. ✅
- M2: (80/200)/((80/200)+(20/100)) = 0.40/0.60 = **2/3 = 0.6666…**. ✅
- A5SS/A3SS assignment (numeric, + strand, ascending coords):
  - exons sharing START, differing END → alternative donor → **A5SS** (rMATS: differ at END boundary).
  - exons sharing END, differing START → alternative acceptor → **A3SS** (rMATS: differ at START boundary).

### Findings / divergences (Stage A)
- **PASS-WITH-NOTES.** The biology/maths of the description are correct *except* the A5SS/A3SS labels in the **TestSpec (M8/M9)** and the **Evidence dataset table** were **swapped** relative to the standard convention. The Evidence row for "AlternativeFivePrimeSS" even used a *different* exon pair (first exon `(1,100)` vs `(1,150)`) than the test it claims to back. These description defects were corrected this session (TestSpec M8/M9 + §5.6 table; Evidence dataset table; added the rMATS coordinate-convention evidence point). The PSI formulas, SE/RI/MXE definitions, edge cases and invariants were all confirmed correct.

---

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs`
  - `CalculatePSI` (lines ~761–786)
  - `DetectAlternativeSplicing` (lines ~804–831) → `ClassifyIsoformPair` (~838–891) and helpers `SpansIntron`/`Spans`/`Overlaps`/`MakeEvent`.

### Formula realised correctly?
- **`CalculatePSI`**: ✅ realises both forms exactly. Negative reads → `ArgumentOutOfRangeException`; both lengths > 0 → rMATS normalized form (NaN when both rates 0); otherwise Ψ=I/(I+S) (NaN when total 0). All matched the sourced numbers.
- **`DetectAlternativeSplicing` / `ClassifyIsoformPair`**: SE, RI, MXE branches correct. **A5SS/A3SS branches were SWAPPED** — the **defect**:
  - Old: `Start==Start && End!=End → AlternativeThreePrimeSS`; `End==End && Start!=Start → AlternativeFivePrimeSS`.
  - This is the reverse of the rMATS/biology convention (shared-start/different-end is an alternative **donor** = **A5SS**; shared-end/different-start is an alternative **acceptor** = **A3SS**).
  - **Fix applied** (lines ~875–888): the two labels were swapped and the rationale + rMATS source documented in a code comment.

### Cross-verification table recomputed vs (fixed) code
| Case | Input pair | Sourced expectation | Fixed-code output |
|------|-----------|---------------------|-------------------|
| M1 | I=80,S=20 | 0.80 | 0.80 ✅ |
| M2 | I=80,S=20,lᵢ=200,lₛ=100 | 0.6666… | 0.6666… ✅ |
| M3 | I=50,S=0 | 1.0 | 1.0 ✅ |
| M4 | I=0,S=40 | 0.0 | 0.0 ✅ |
| M5 | I=0,S=0 | NaN | NaN ✅ |
| M6 | A=(1,100),(200,300),(400,500); B=(1,100),(400,500) | SkippedExon | SkippedExon ✅ |
| M7 | A=(1,100),(200,300); B=(1,300) | RetainedIntron | RetainedIntron ✅ |
| M8 | A=…(200,300); B=…(200,350) (shared START, diff END) | **A5SS** | **AlternativeFivePrimeSS** ✅ (was A3SS before fix) |
| M9 | A=…(200,300); B=…(150,300) (shared END, diff START) | **A3SS** | **AlternativeThreePrimeSS** ✅ (was A5SS before fix) |
| M10 | two non-overlapping middle exons | MutuallyExclusiveExons | MutuallyExclusiveExons ✅ |

### Variant/delegate consistency
- No other caller in `src/` references the `AlternativeFivePrimeSS`/`AlternativeThreePrimeSS` labels (grep clean), so the fix is self-contained; the unrelated `FindSkippedExonEvents` PSI path is unaffected.

### Test quality audit (HARD gate)
- **Before:** M8/M9 tests, the TestSpec and the Evidence all encoded the *swapped* labels — classic code-echoes that would pass against the wrong implementation. This is itself a Stage-B defect.
- **After (this session):**
  - M8/M9 renamed and re-asserted to the **sourced** event types, with the splice-site rationale and rMATS source in the comment; **added exact Start/End span assertions** (M8 → [200,350]; M9 → [150,300]) to also lock the `MakeEvent` coordinate logic.
  - PSI tests assert exact sourced values within 1e-10 (M2 uses `2.0/3.0`, not a code echo); NaN, bounds, negative-throw, partial-length fallback all covered.
  - All 5 event classes + single-isoform, identical-isoform, null/empty, gene-id attribution covered.
- **Honest green:** FULL unfiltered suite = **Failed: 0, Passed: 6501**. `dotnet build` = **0 errors**; the 4 build warnings are pre-existing `NUnit2007` in `ApproximateMatcher_EditDistance_Tests.cs` (untouched) — the files changed this session build warning-free. **Gate: PASS.**

### Findings / defects (Stage B)
- **Defect A15 (real algorithm bug):** swapped A5SS/A3SS labels in `ClassifyIsoformPair`. **Fixed** in-session (code + tests + TestSpec + Evidence), logged in `FINDINGS_REGISTER.md` §A15.

---

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — formulas/taxonomy/edge-cases correct; swapped A5SS/A3SS labels in TestSpec+Evidence corrected this session.
- **Stage B: FAIL → fixed** — the implementation swapped the A5SS/A3SS labels; corrected, tests strengthened to sourced values + coordinate assertions.
- **End-state: ✅ CLEAN** — code corrected, tests lock the sourced behaviour, `dotnet build` 0 errors, full suite 6501 passed / 0 failed.
- No deferred follow-ups.
