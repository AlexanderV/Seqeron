# Validation Report: PROBE-VALID-001 — Hybridization Probe Validation

- **Validated:** 2026-06-12   **Area:** MolTools
- **Canonical method(s):** `ProbeDesigner.ValidateProbe(probeSequence, referenceSequences, maxMismatches=3, selfComplementarityThreshold=0.3)`, `ProbeDesigner.CheckSpecificity(probeSequence, genomeIndex)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_ProbeValidation_Tests.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: Hybridization probe** — Confirms specificity is governed by **stringency** (temperature/salt): high stringency permits only highly-similar duplexes; low stringency tolerates mismatched (cross-hybridizing) duplexes. Confirms the qualitative model that cross-hybridization is a sequence-similarity (mismatch) phenomenon. Notes a probe may still hybridize to an unknown/undesired target — i.e., specificity is never absolute.
- **Wikipedia: DNA microarray** — Cross-hybridization between probe and partially-complementary off-target sequences is a known artifact in high-density arrays; specificity is assessed by signal from intended vs. non-intended targets.
- **Wikipedia: Off-target genome editing** — Confirms **CRISPR/Cas9 tolerates 3–5 bp mismatches per 20 nt guide**; more mismatches generally lower off-target activity, but position (seed region near PAM) and GC content modulate tolerance. Supports `maxMismatches = 3` as the defensible **lower bound** of mismatch tolerance.
- **Amann & Ludwig (2000), *FEMS Microbiol. Rev.* 24:555–565** — Foundational rRNA-targeted probe specificity reference; specificity = probe binds intended target and discriminates against non-targets primarily by number/position of mismatches and resulting duplex stability (Tm).
- **PLOS One 2018 (oligoprobe on/off-target trade-offs)** and **NAR 36:2395 (microarray specificity)** — Confirm specificity is a *trade-off*: both on-target and cross-hybridization signals rise with binding energy; specificity is the ratio of on-target to off-target binding. Real specificity scoring uses binding-energy/Tm, k-mer genome occurrence, and self-folding — i.e., it is thermodynamic, not a flat mismatch count.

### Validation-logic check (the defined description)
- "Specific = matches intended target but does NOT hybridize elsewhere" — correct in the abstract.
- "Cross-hyb risk assessed by sequence similarity (mismatch count) with a threshold" — correct; the implementation uses **Hamming distance ≤ maxMismatches** over the full probe length (ungapped, fixed-length window). This is a defensible *heuristic* approximation of the similarity criterion.
- Thresholds are sourced and defensible: `maxMismatches=3` (CRISPR lower bound, Hsu 2013 / Fu 2013 via Wikipedia); `selfComplementarityThreshold=0.3` (random-DNA baseline ≈0.25, +20% headroom; per-application defaults qPCR 0.25 → FISH/Southern 0.4).

### Edge-case semantics (Stage-A defined)
- Unique probe → specific (1 hit → score 1.0). Defined, sensible.
- Exact off-target match → flagged (N hits → score 1/N, issue emitted when N>1). Defined.
- Near-threshold identity → governed by `maxMismatches`. Defined.
- Empty off-target set / 0 matches → score 0.0. Defined (see note below).

### Findings / divergences (why PASS-WITH-NOTES, not PASS)
1. **Heuristic, not BLAST-grade / not thermodynamic.** The TestSpec name "BLAST" (in the checklist source list) and "Tm of off-target duplex" (in the validation brief) are NOT realized. Off-target search is a **naive ungapped Hamming-distance substring scan** with a flat mismatch cap — no gaps/indels, no seed weighting, no positional mismatch model, no E-value, and **no Tm of the off-target duplex**. This is a *declared heuristic*; the code comments and TestSpec attribute the invariants to "Implementation," which is honest. The specificity claim must NOT be over-stated as BLAST-grade — flagged here per protocol.
2. **Specificity score conflates on-target and off-target hits.** `OffTargetHits` counts ALL matches across the pooled reference sequences, including the intended target. The score `1/totalHits` therefore treats the on-target match itself as reducing specificity, and `0 matches → 0.0` means "matches nothing" (useless probe), not "no off-targets." This is an internally-consistent declared convention (Invariants #4/#5/#6 sourced to "Implementation"), but it is NOT the literature ratio of on-target:off-target signal. It is a usable relative-uniqueness heuristic, correctly documented.

These are documentation/scope notes, not biological errors. The abstract logic (similarity-based cross-hyb detection with a sourced mismatch threshold) is correct.

---

## Stage B — Implementation

### Code path reviewed
- `ValidateProbe` — `ProbeDesigner.cs:491-564`
- `CheckSpecificity` — `ProbeDesigner.cs:569-587`
- `FindApproximateMatches` (Hamming-distance window scan) — `ProbeDesigner.cs:789-804`
- `CalculateSelfComplementarity` — `ProbeDesigner.cs:721-733`
- `HasSecondaryStructurePotential` — `ProbeDesigner.cs:735-761`

### Realises the validated (heuristic) description? — Yes
- Off-target detection: `FindApproximateMatches` slides the probe across each `ToUpperInvariant()` reference, counts mismatches with early-exit at `maxMismatches`, yields each start with `mismatches ≤ maxMismatches`. Exactly the Hamming-≤-k heuristic validated in Stage A.
- Specificity decision: 0 hits→0.0, 1 hit→1.0, N hits→1/N (`:548-553`). Matches Invariants #4/#5/#6 and `CheckSpecificity` (suffix-tree exact-match variant, `:579-586`) — consistent.
- Self-complementarity, secondary-structure, and issue/`IsValid` logic match the TestSpec.

### Cross-verification (recomputed by hand vs. code)
| Case | Input | Hand calc | Code/test result | Match |
|------|-------|-----------|------------------|-------|
| Exact off-target (M4) | `AAAAAAAAAA` in 34×A | 34−10+1 = 25 hits, 1/25 = 0.04 | offTargetHits=25, spec=0.04 | ✅ |
| Many-mismatch passes (S4) | `ACGTACGTACGTACGT` vs `TTTTTACGAACGAACGAACGTTTTT` | pos5 sub `ACGAACGAACGAACGT` = 3 mismatches; strict→0, approx(3)→1 | strict 0, approx 1, spec 1.0 | ✅ (Python-verified) |
| Unique (M3) | `ATCGATCG…` in single-match ref | 1 hit → 1.0 | offTargetHits=1, spec=1.0 | ✅ |
| Empty off-target set (M2) | StandardProbe, no refs | 0 hits → 0.0 | offTargetHits=0, spec=0.0 | ✅ |
| Empty probe (M1) | "" | early return, spec 0.0, invalid | spec=0.0, IsValid=false | ✅ |
| No-match probe (ZeroHits) | poly-T vs ACG-repeat | 0 hits → 0.0 | offTargetHits=0, spec=0.0 | ✅ |

### Variant/delegate consistency
`CheckSpecificity` (suffix-tree, exact) and `ValidateProbe` (Hamming, approximate) share the same 0/1/N→0.0/1.0/(1/N) scoring rule — consistent. `DesignProbes(genomeIndex,…)` calls `CheckSpecificity` and filters non-unique probes when `requireUnique`.

### Test quality audit
19 ProbeValidation tests, all exact-value assertions (not "no-throw" tautologies): exact hit counts (25, 16, 12, 3, 1, 0), exact specificity (0.04, 1/3, 1.0, 0.0), exact selfComp (1.0), issue-string contents, IsValid booleans. Edge cases (empty probe, empty refs, null refs→throws, mixed case, long ref, multi ref, approximate vs strict) all covered. Numerically robust (no div-by-zero: 0-hit branch short-circuits before `1/N`).

### Findings / defects
None. The code faithfully and correctly realises its declared heuristic specificity model. No precision, overflow, or edge-case defects found.

---

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** — abstract similarity-based cross-hyb logic is biologically correct and thresholds are sourced/defensible; specificity is a documented **heuristic (Hamming-≤-k substring scan + 1/N uniqueness score)**, NOT BLAST-grade and NOT Tm/thermodynamic, and `OffTargetHits` pools on- and off-target matches. Correctly attributed to "Implementation" in the spec, so not overclaimed.
- **Stage B: PASS** — implementation matches the validated description exactly; both worked examples recompute correctly; all edge cases handled.
- **End-state: CLEAN** — no defect; no code changes required. Build green, 19/19 ProbeValidation tests pass, full suite 4461/4461 pass.
- **Follow-up (non-blocking):** if BLAST-grade or Tm-aware off-target scoring is ever claimed in user-facing docs, upgrade `FindApproximateMatches` to a gapped/seed-and-extend search with duplex-Tm filtering, and split `OffTargetHits` into intended-target vs. off-target counts. Until then, keep the "heuristic" framing.
