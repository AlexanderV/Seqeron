---
type: source
title: "Validation report: TRANS-SPLICE-001 (alternative splicing — event classification + Percent-Spliced-In)"
tags: [validation, transcriptome, governance]
doc_path: docs/Validation/reports/TRANS-SPLICE-001.md
sources:
  - docs/Validation/reports/TRANS-SPLICE-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: TRANS-SPLICE-001

The two-stage **validation write-up** for test unit **TRANS-SPLICE-001** (Alternative
Splicing — event classification + Percent-Spliced-In / PSI), validated 2026-06-15 (Area:
Transcriptome). This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's **verdict** on both the algorithm
description and the shipped code. The two-stage methodology is the [[validation-protocol]];
the splicing algorithm itself is summarized in [[alternative-splicing-psi]]. Distinct from
the pre-implementation [[trans-splice-001-evidence]] artifact.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: FAIL → fixed · State: ✅ CLEAN.** A real algorithm
bug — the **A5SS/A3SS event labels were swapped** in `ClassifyIsoformPair` — was found and
**corrected** in-session (code + tests + TestSpec + Evidence), logged FINDINGS_REGISTER
§A15. The PSI formulas, SE/RI/MXE definitions, edge cases and invariants were all confirmed
correct. Full unfiltered suite **6501 passed / 0 failed**, `dotnet build` 0 errors (the 4
warnings are pre-existing NUnit2007 in an untouched file); changed files build warning-free.

Canonical methods: `TranscriptomeAnalyzer.CalculatePSI(inclusionReads, exclusionReads,
inclusionEffectiveLength?, exclusionEffectiveLength?)`; `DetectAlternativeSplicing(isoforms)`
(+ private `ClassifyIsoformPair`).

## Stage A — description (algorithm faithfulness)

- Sources opened this session: **rMATS-turbo README** (A5SS/A3SS coordinate columns — on the
  + strand, A5SS forms differ at their downstream END/3′ boundary sharing the upstream
  START/5′; A3SS differ at their upstream START/5′ boundary sharing the downstream END/3′);
  **NAR 34(21):6305** and molecular-biology splice-site definitions (5′ = donor = 3′ END of
  the upstream exon; 3′ = acceptor = 5′ START of the downstream exon); **SUPPA2 / biostars /
  Outrigger** (PSI = inclusion / (inclusion + skipping)); carried from Evidence: PMC3330053
  (Ψ = γᵢ/(γᵢ+γₑ)), Shen 2014 rMATS PMC4280593 (ψ̂ = (I/lᵢ)/(I/lᵢ + S/lₛ)), Wang 2008
  five-class taxonomy.
- Formulas confirmed verbatim: unnormalized **Ψ = I/(I+S)** (PMC3330053, SUPPA2); rMATS
  **length-normalized ψ̂ = (I/lᵢ)/(I/lᵢ + S/lₛ)** (Shen 2014); five-class taxonomy SE, RI,
  A5SS, A3SS, MXE (Wang 2008). Numeric anchors: M1 80/(80+20)=0.80; M2
  (80/200)/((80/200)+(20/100))=2/3=0.6666…
- Edge cases defined: 0/0 → NaN; S=0,I>0 → 1; I=0,S>0 → 0; <2 isoforms or identical isoforms
  → no event.
- **PASS-WITH-NOTES:** biology/maths correct **except** the A5SS/A3SS labels in the TestSpec
  (M8/M9) and the Evidence dataset table were **swapped** relative to the standard convention
  (the Evidence "AlternativeFivePrimeSS" row even used a different exon pair than the test it
  backs). These description defects were corrected this session (TestSpec M8/M9 + §5.6 table;
  Evidence dataset table; added the rMATS coordinate-convention evidence point).

## Stage B — implementation (code review + cross-check)

- Code path: `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs`
  — `CalculatePSI` (~761–786), `DetectAlternativeSplicing` (~804–831) → `ClassifyIsoformPair`
  (~838–891) with helpers `SpansIntron`/`Spans`/`Overlaps`/`MakeEvent`.
- `CalculatePSI` realises both PSI forms exactly: negative reads → `ArgumentOutOfRangeException`;
  both lengths > 0 → rMATS normalized form (NaN when both rates 0); otherwise Ψ=I/(I+S) (NaN
  when total 0). SE/RI/MXE branches correct.
- **Defect A15 (the FAIL):** the **A5SS/A3SS branches were swapped** — old code mapped
  `Start==Start && End!=End → AlternativeThreePrimeSS` and `End==End && Start!=Start →
  AlternativeFivePrimeSS`, the reverse of the rMATS/biology convention (shared-start /
  different-end = alternative **donor** = **A5SS**; shared-end / different-start = alternative
  **acceptor** = **A3SS**). **Fix applied** (~875–888): the two labels swapped, with the
  rationale + rMATS source in a code comment. No other caller in `src/` references these labels
  (grep clean), so the fix is self-contained; the `FindSkippedExonEvents` PSI path is unaffected.
- Cross-verification after fix (M1–M10): PSI M1 0.80, M2 0.6666…, M3 1.0, M4 0.0, M5 NaN;
  M6 SkippedExon, M7 RetainedIntron, **M8 (shared START, diff END) → AlternativeFivePrimeSS
  (was A3SS before fix)**, **M9 (shared END, diff START) → AlternativeThreePrimeSS (was A5SS
  before fix)**, M10 MutuallyExclusiveExons — all match sourced expectations.
- Test-quality audit (HARD gate): before the fix, M8/M9 tests, the TestSpec and the Evidence
  all encoded the *swapped* labels — code-echoes that would pass against the wrong
  implementation (itself a Stage-B defect). After: M8/M9 renamed and re-asserted to the
  **sourced** event types with the splice-site rationale + rMATS source in the comment, and
  **exact Start/End span assertions added** (M8 → [200,350]; M9 → [150,300]) to also lock the
  `MakeEvent` coordinate logic. PSI tests assert exact sourced values within 1e-10 (M2 uses
  `2.0/3.0`, not a code echo); all 5 event classes + single/identical isoform, null/empty,
  gene-id attribution covered. **Gate: PASS.**
- **Stage B: FAIL → fixed** — swapped A5SS/A3SS labels corrected; tests strengthened to
  sourced values + coordinate assertions.

## Findings

- **State ✅ CLEAN.** Defect **A15** (swapped A5SS/A3SS labels in `ClassifyIsoformPair`)
  **fixed** in-session across code, tests, TestSpec and Evidence; logged in
  `FINDINGS_REGISTER.md` §A15. No deferred follow-ups.

See the report at `docs/Validation/reports/TRANS-SPLICE-001.md`.
