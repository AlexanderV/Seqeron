---
type: source
title: "Validation report: ANNOT-CODING-001 (coding potential — CPAT hexamer usage-bias score)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/ANNOT-CODING-001.md
sources:
  - docs/Validation/reports/ANNOT-CODING-001.md
source_commit: e748206486a14ab05fe3c14e312e74cd77874af2
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ANNOT-CODING-001

The two-stage **validation write-up** for test unit **ANNOT-CODING-001** (Coding Potential
Calculation — the CPAT hexamer usage-bias score), validated 2026-06-15. This is the *report*
artifact that feeds one row of the [[validation-ledger]]; it records the validator's **verdict**
on both the algorithm description and the shipped code. The score itself is summarized in
[[coding-potential-hexamer-score]]; the two-stage methodology is the [[validation-protocol]].
Distinct from the pre-implementation [[annot-coding-001-evidence]] artifact.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS (after in-session fix) · State: ✅ CLEAN.** A genuine
code defect was found **and completely fixed** in-session. Full unfiltered suite **6561 passed /
0 failed**, build 0 errors (4 pre-existing unrelated warnings). Test-quality gate PASS.

## The defect (both-zero hexamer counted, not skipped)

The canonical method is `GenomeAnnotator.CalculateCodingPotential(...)` at
`GenomeAnnotator.cs:571-618`. Canonical CPAT `kmer_ratio` (`liguowang/cpat`
`src/cpmodule/FrameKmer.py`, and the identical `WGLab/lncScore` copy — both re-opened this
session 2026-06-15) has an explicit `elif coding[k]==0 and noncoding[k]==0: continue` branch:
a hexamer absent from **both** tables is skipped and **not** counted.

The shipped code (line ~610) had **no** both-zero branch — it fell through and unconditionally ran
`scoredHexamers++`, so a both-zero hexamer was counted as a scored **0**, inflating the denominator
and diluting the mean whenever such a hexamer occurs. **Fix:** added
`else if (coding == 0 && noncoding == 0) continue;` plus a trailing `else continue;` before the
counter, verbatim to the canonical reference. This mirror-defect also lived in three docs and one
test (below) and was corrected across all of them.

## Stage A — description (algorithm faithfulness)

- Formula confirmed against the CPAT paper (Wang et al. 2013, *NAR* 41(6):e74) and canonical
  `kmer_ratio` (frame-0 loop): score = **mean of `ln(coding[k]/noncoding[k])`** over in-frame
  hexamers (word 6, step 3, frame 0, full-length words only), `math.log` = natural log; positive
  ⇒ coding, negative ⇒ noncoding.
- Four contribution branches sourced verbatim: both>0 ⇒ `+= ln(coding/noncoding)`, count++;
  coding-only ⇒ `+= 1`, count++; noncoding-only ⇒ `-= 1`, count++; **both-zero ⇒ `continue`,
  NOT counted**; missing-from-either-table ⇒ `continue`; `len(seq) < word_size` ⇒ 0.
- Independent hand cross-check (5 rows): `ATGAAACCC` with C{ATGAAA:8,AAACCC:2}/N{ATGAAA:2,AAACCC:4}
  ⇒ (ln4 + ln0.5)/2 = **0.34657359027997264**; ln4/1 = **1.3862943611198906**; coding-only ⇒
  **1.0**; noncoding-only ⇒ **−1.0**; **both-zero case: only the scorable hexamer counts** ⇒ ln4
  = **1.3862943611198906** (the discriminating value that proves the fix).
- **DEFECT (description):** algorithm-doc §2.2, the Evidence `kmer_ratio` transcription, and TestSpec
  C1 all omitted the both-zero `continue` and described it as "contributes 0 but is counted."
  Corrected this session in all three docs.
- ASSUMPTION-1 retained: no-scorable-hexamer returns **0** (port) vs reference **−1** — both are
  non-scores; defensible documented port choice.

## Stage B — implementation (code review + cross-check)

- After the fix, the formula is realised correctly: in-frame extraction (offset 0, step `stepSize`,
  full-length `wordSize`), `Math.Log` (= ln), the four branches, skip-if-missing, mean over scored
  count, `len < wordSize → 0`, `count == 0 → 0` (ASSUMPTION-1).
- All five cross-verification rows reproduced by the suite (6561 passed). C1 now asserts the sourced
  **1.3862943611198906** — it previously asserted the code-echoing **0.6931471805599453**.
- Variant consistency: single canonical method; MCP `AnnotationTools.cs:144` delegates with no
  behavioural dependency on the fixed branch.
- Test-quality audit (HARD gate) PASS: C1 was a **code-echo** (asserted the impl's wrong both-zero
  behaviour) and was rewritten to the canonical sourced value; exact `Is.EqualTo(...).Within(1e-10)`
  throughout; no green-washing, no skips; all Stage-A branches covered (both-positive, ±1
  pseudo-scores, both-zero skip, missing-key skip, in-frame stepping, too-short/empty→0, no-scorable
  →0, case-insensitivity, four validation throws). The out-of-contract negative-value `else: continue`
  branch is left untested by design (inputs constrained non-negative).

## Findings

- **FIND (fixed):** both-zero hexamer counted instead of skipped — one code defect + one code-echo
  test (C1) + three docs (algorithm doc §2.2, Evidence, TestSpec). All corrected in-session; logged
  in [[findings-register|FINDINGS_REGISTER]]. End-state ✅ CLEAN.
- This is a **test-quality catch**: the pre-existing green test encoded the buggy implementation, so
  the suite was passing on a wrong value until the validator re-anchored C1 to the canonical source.
