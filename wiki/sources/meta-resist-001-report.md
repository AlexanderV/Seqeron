---
type: source
title: "Validation report: META-RESIST-001 (antibiotic-resistance gene detection, ResFinder-style identity/coverage best-match)"
tags: [validation, metagenomics, governance]
doc_path: docs/Validation/reports/META-RESIST-001.md
sources:
  - docs/Validation/reports/META-RESIST-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: META-RESIST-001

The two-stage **validation write-up** for test unit **META-RESIST-001** (Antibiotic
Resistance Gene Detection, ResFinder-style), validated 2026-06-15. This is the *report*
artifact that feeds one row of the [[validation-ledger]]; it records the validator's
**verdict** on both the algorithm description and the shipped code. The identity /
coverage / dual-threshold best-match model is summarized in
[[antibiotic-resistance-gene-detection]]; the two-stage methodology is the
[[validation-protocol]]. Distinct from the pre-implementation
[[meta-resist-001-evidence]] artifact.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS (after in-session fix) · End-state: ✅ CLEAN.**
A real tie-break defect (`BestUngappedMatch`) was found and completely fixed in-session,
producing a genuine false-negative under default thresholds before the fix. Full
unfiltered suite **6556 passed / 0 failed** (1 pre-existing benchmark skip); `dotnet
build` 0 errors, no new warnings. Test quality gate PASS.

## Stage A — description (algorithm faithfulness)

- Canonical method:
  `MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, referenceGenes,
  identityThreshold, coverageThreshold)`; internal `BestUngappedMatch(contig, reference)`.
- Sources opened: **Heng Li (2018)** "On the definition of sequence identity" (BLAST
  identity = matches / alignment columns; gapless ⇒ denominator = window length; worked
  43/50 = 86%); **Zankari et al. (2012) ResFinder** (best-matching gene single output;
  ≥ 2/5 reference-length coverage floor; default ID 100%); **ResFinder GitHub** (CLI
  defaults `-t 0.80` / `-l 0.60`; coverage = breadth-of-coverage vs reference);
  **web-search corroboration** (90% identity / 60% coverage operating point); **CARD RGI**
  ("Perfect" match = 100% over full reference length, best hit by bit-score).
- Formula check: percent identity = matches / window (gapless) ✔; coverage = window /
  reference length ✔; dual-threshold reporting (identity ≥ id AND coverage ≥ cov) ✔;
  best-matching gene per contig ✔.
- Edge cases: empty contig skipped, empty reference ignored, null →
  `ArgumentNullException`, threshold ∉ [0,1] → `ArgumentOutOfRangeException`.
- **Verdict PASS-WITH-NOTES** (see Findings for N-A1 threshold-version note and N-A2
  description tie-break inconsistency fixed this session).

## Stage B — implementation (code review + cross-check)

- Code path:
  `Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs:1084-1191`
  (`FindAntibioticResistanceGenes` + `BestUngappedMatch`). Legacy `FindResistanceGenes`
  (motif stub, line 1196) is a separate MCP method, out of scope.
- **Defect F-RESIST-001 (fixed) — tie-break direction in `BestUngappedMatch`:** the
  original objective maximized match count but on ties preferred the *longer* window. A
  longer window is only reached by padding with mismatching flanks, which lowers identity
  — not how BLAST reports an HSP. Outcome-changing probe: `contig=GGGACGTACG`,
  `ref=ACGTACGTAC` — the true perfect 7-base suffix HSP (identity 1.0, coverage 7/10 =
  0.70) **passes** 0.90/0.60, but the old code picked the padded 9-wide window (identity
  0.778) which **fails** the identity threshold → gene **wrongly missed (false negative)**.
- **Fix:** tie-break now prefers the **shorter (higher-identity)** window
  (`matches == bestMatches && bestWindow != 0 && window < bestWindow`), with
  `bestWindow == 0` marking the unset state. Doc §4.2 corrected.
- Cross-verification (independent Python model, fixed code): M1 exact 1.0/1.0 ✔; M2
  1-mismatch 6/7≈0.857/1.0 ✔; M3 edge-4 1.0/0.571 ✔; M4 partial-5 1.0/0.714 ✔; M5 low-id
  0.714/1.0 ✔; M6 best-of-2 geneA wins ✔; C1 tie→coverage geneFull wins ✔; edge-perfect
  probe 1.0/0.70 reported ✔. All values trace to Li (2018) / Zankari (2012), not code.

## Findings

- **N-A1 (NOTE):** default identity threshold **0.90** — a documented, version-specific
  web-service operating point (2012 paper default 100%; GitHub CLI default 0.80).
  Recorded as a sourced note, not a defect.
- **N-A2 (description fixed):** doc §4.2 documented the tie-break as "ties → longer
  window," but the Evidence M3 value is only producible with a shorter/higher-identity
  tie-break; the documented rule is biologically wrong. Description corrected.
- **F-RESIST-001 (Stage B, fixed):** `BestUngappedMatch` tie-break preferred the longer
  mismatch-padded window (real false negative at default thresholds); fixed to prefer the
  shorter/higher-identity window; M3 strengthened (was `Is.Empty` only, now locks
  identity 1.0 / coverage 4/7); added **M3b** regression guard
  (`FindAntibioticResistanceGenes_EdgePerfectHsp_PreferredOverPaddedWindow`) that fails
  against the old code. Honest green: full suite 6556/0. The gapless-vs-gapped-BLAST
  simplification remains a documented, accepted scope limitation (ASM-01), not a defect.
