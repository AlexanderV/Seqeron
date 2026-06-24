# Validation Report: PROBE-VALID-001 — Hybridization Probe Validation

- **Validated:** 2026-06-24   **Area:** MolTools
- **Canonical method(s):** `ProbeDesigner.ValidateProbe(probeSequence, referenceSequences, maxMismatches=3, selfComplementarityThreshold=0.3)`, `ProbeDesigner.CheckSpecificity(probeSequence, genomeIndex)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_ProbeValidation_Tests.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

This is an independent re-validation (fresh context). The code and tests are unchanged
since the prior validation (last source commit `6f2ea3ef`, a coverage-classification audit;
core logic identical). All cross-check numbers below were re-derived by hand (Python),
not copied from the prior report.

---

## Stage A — Description

### What the unit actually checks (per TestSpec + code)
`ValidateProbe` returns a `ProbeValidation` record with five signals:
1. **OffTargetHits** — count of approximate matches (Hamming distance ≤ `maxMismatches`,
   ungapped fixed-length window) of the probe across the pooled reference sequences.
2. **SpecificityScore** — derived from hit count: `0 → 0.0`, `1 → 1.0`, `N>1 → 1/N`.
3. **SelfComplementarity** — fraction of positions where `seq[i] == revComp(seq)[i]`, in [0,1].
4. **HasSecondaryStructure** — inverted-repeat (hairpin) scan: stem ≥4 nt, loop gap = 3,
   stem complementarity ≥ 80%.
5. **Issues / IsValid** — human-readable issues plus an `IsValid` verdict.

The "valid" verdict (`ProbeDesigner.cs:555`):
`isValid = issues.Count == 0 || (offTargetHits <= 1 && selfComp <= 0.4)`.

### Sources opened & what they confirm
- **Wikipedia: Hybridization probe** — specificity is governed by **stringency**
  (temperature/salt); high stringency permits only highly-similar duplexes, low stringency
  tolerates mismatched (cross-hybridizing) duplexes. Confirms the qualitative model that
  cross-hybridization is a sequence-similarity (mismatch) phenomenon; a probe may still
  hybridize to an unknown target (specificity is never absolute).
- **Wikipedia: DNA microarray** — cross-hybridization to partially-complementary off-targets
  is a known artifact in high-density arrays; specificity = intended vs. non-intended signal.
- **Wikipedia: Off-target genome editing** — CRISPR/Cas9 tolerates **3–5 bp mismatches per
  20 nt guide** (Hsu 2013, Fu 2013), with seed-region and GC modulation. Supports
  `maxMismatches = 3` as the defensible **lower bound** of mismatch tolerance.
- **Wikipedia: Nucleic acid thermodynamics** — mismatches destabilize duplexes; real
  specificity scoring is thermodynamic (binding-energy/Tm), of which a flat mismatch cap is
  a heuristic surrogate. Random-DNA self-complementarity baseline ≈ 0.25 supports the 0.3
  threshold (+20% headroom; per-application defaults qPCR 0.25 → FISH/Southern 0.4).

These align with PROBE-DESIGN-001's source set (same Wikipedia references + SantaLucia/Breslauer
for the Tm side, which is not exercised by *validation*).

### Validation-logic check (abstract correctness)
- "Specific = matches intended target, does not hybridize elsewhere" — correct in the abstract.
- "Cross-hyb risk via sequence similarity (mismatch count) with a threshold" — correct; the
  implementation uses **Hamming distance ≤ maxMismatches** over the full probe length
  (ungapped, fixed-length window). A defensible *heuristic* of the similarity criterion.
- Self-complementarity and hairpin checks are standard probe-QC concerns; thresholds sourced.

### Edge-case semantics (Stage-A defined)
- Empty probe → early return: spec 0.0, hits 0, selfComp 0.0, IsValid false, issue emitted. Defined.
- Empty / no references → 0 hits → spec 0.0. Defined.
- Null references → `ArgumentNullException`. Defined.
- Unique probe (1 hit) → spec 1.0. Exact off-target (N hits) → spec 1/N, issue when N>1. Defined.
- Near-threshold identity governed by `maxMismatches`. Defined.

### Findings / divergences (why PASS-WITH-NOTES)
1. **Heuristic, not BLAST-grade / not thermodynamic.** The test-file doc-comment lists
   "Wikipedia: BLAST (approximate matching algorithms)" and the task brief mentions
   "Tm of off-target duplex." Neither is realized: off-target search is a **naive ungapped
   Hamming-distance substring scan** with a flat mismatch cap — no gaps/indels, no seed
   weighting, no positional mismatch model, no E-value, no off-target duplex Tm. The
   invariants are honestly attributed to "Implementation" in the TestSpec, so this is a
   *declared heuristic*, not an over-claim — but it must not be presented as BLAST-grade.
2. **SpecificityScore conflates on-target and off-target hits.** `OffTargetHits` counts ALL
   matches across the pooled references, including the intended target; `1/totalHits` therefore
   treats the on-target match itself as reducing specificity, and `0 matches → 0.0` means
   "matches nothing" (useless probe), not "no off-targets." Internally consistent and
   documented (Invariants #4/#5/#6 sourced to "Implementation"), but it is a relative-uniqueness
   heuristic, not the literature on-target:off-target signal ratio.
3. **`IsValid` ignores `hasStructure` in its second clause.** A probe with detected secondary
   structure but ≤1 hit and selfComp ≤ 0.4 yields `IsValid = true` while still emitting a
   structure issue. This is a deliberate "soft" QC design (structure is advisory, not
   disqualifying); consistent with the tested cases but worth noting as a documented choice.

These are documentation/scope notes, not biological errors. The abstract logic
(similarity-based cross-hyb detection with sourced thresholds, plus standard self-structure QC)
is correct.

---

## Stage B — Implementation

### Code path reviewed
- `ValidateProbe` — `ProbeDesigner.cs:491-564`
- `CheckSpecificity` (suffix-tree, exact) — `ProbeDesigner.cs:569-587`
- `FindApproximateMatches` (Hamming-window scan, early-exit) — `ProbeDesigner.cs:789-804`
- `CalculateSelfComplementarity` — `ProbeDesigner.cs:721-733`
- `HasSecondaryStructurePotential` — `ProbeDesigner.cs:735-761`

### Realises the validated (heuristic) description? — Yes
- Off-target: `FindApproximateMatches` slides probe across each `ToUpperInvariant()` reference,
  counts mismatches with early-exit at `maxMismatches`, yields each start with `mm ≤ maxMismatches`.
  Exactly the Hamming-≤-k heuristic validated in Stage A.
- Specificity: `0→0.0, 1→1.0, N→1/N` (`:548-553`); `CheckSpecificity` applies the same rule
  over exact suffix-tree occurrences (`:579-586`) — consistent.
- SelfComplementarity, hairpin, and issue/`IsValid` logic match the TestSpec.

### Cross-verification (re-derived by hand in Python vs. code/tests)
| Case | Input | Hand calc | Test asserts | Match |
|------|-------|-----------|--------------|-------|
| M4 exact off-target | `A×10` in `A×34` | 34−10+1 = **25** hits, 1/25 = 0.04 | hits=25, spec=0.04 | ✅ |
| S2 problematic | `A×10` in `A×25` | 25−10+1 = **16** hits | hits=16, "16 potential off-target sites" | ✅ |
| S3 multi-problem | `GCGCGCGCGC` in `(GC)×16` | even positions only = **12** hits; selfComp **1.0** | hits=12, selfComp=1.0, 2 issues, IsValid=false | ✅ |
| S4 approx | `ACGTACGTACGTACGT` vs `TTTTTACGAACGAACGAACGTTTTT` | strict(0)→**0**, approx(3)→**1** | strict 0, approx 1, spec 1.0 | ✅ |
| M8 palindrome | `(GC)×10` | selfComp = **1.0** | selfComp=1.0, "Self-complementarity" issue | ✅ |
| M3 unique | `UniqueProbe` in single-match ref | **1** hit → 1.0 | hits=1, spec=1.0 | ✅ |
| C2 multi-ref | `UniqueProbe` in 3 refs | **3** hits → 1/3 | hits=3, spec=1/3 | ✅ |

S3 note: default `maxMismatches=3`, but odd-position offsets of a GC-repeat give 10 mismatches
(every base shifts), so only the 12 even positions match — confirmed independent of mm=3.

### Variant/delegate consistency
`CheckSpecificity` (suffix-tree, exact) and `ValidateProbe` (Hamming, approximate) share the
identical `0/1/N → 0.0/1.0/(1/N)` scoring rule. `DesignProbes(genomeIndex,…)` calls
`CheckSpecificity` and filters non-unique probes when `requireUnique` (integration test confirms).

### Numerical robustness
No div-by-zero (the `0-hit` branch short-circuits before `1/N`); selfComp denominator is the
non-empty probe length (empty probe handled by the early return). Hits accumulate as `int`;
no overflow on stated ranges. Early-exit in the mismatch loop bounds inner cost.

### Test quality audit
19 tests, all **exact-value** assertions (not "no-throw" tautologies): exact hit counts
(25, 16, 12, 3, 1, 0), exact specificity (0.04, 1/3, 1.0, 0.0), exact selfComp (1.0), exact
issue strings, IsValid booleans. Edge cases (empty probe, empty refs, null→throws, mixed case,
20k-char ref, multi-ref, approximate vs strict, hairpin) all covered. Deterministic.

### Findings / defects
None. The code faithfully realises its declared heuristic specificity + self-structure QC model.

---

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** — abstract similarity-based cross-hyb logic and self-structure QC
  are correct, thresholds sourced/defensible; the default `ValidateProbe` specificity is a
  documented **heuristic** (Hamming-≤-k substring scan + 1/N uniqueness) whose `OffTargetHits`
  pools on- and off-target matches.
- **Stage B: PASS** — implementation matches the validated description exactly; all seven
  hand-recomputed cross-checks agree; all edge cases handled; tests are real and exact.
- **End-state: CLEAN** — no defect; no code change. Build green, 19/19 ProbeValidation tests pass.

---

## 2026-06-24 — Limitation fix (re-validation pending)

The follow-up below was implemented. An opt-in **gapped (Smith–Waterman) off-target scan** was added
that REUSES the library's validated `SequenceAligner.LocalAlign` (`BlastDna` scoring) and **separates
on-target from off-target hits**, addressing both PASS-WITH-NOTES findings:

- **New method:** `ProbeDesigner.ScanOffTargetsGapped(string probe, IEnumerable<string> references,
  double minIdentity = 0.75, ScoringMatrix? scoring = null)` → `GappedSpecificityResult`
  (`OnTargetHits`, `OffTargetHits`, `OffTargetCount`, `IsSpecific`), with per-hit `GappedProbeHit`
  (`ReferenceIndex`, `Start`, `End`, `Identity`, `Coverage`, `HasGaps`, `AlignedProbe`, `AlignedReference`).
- **Indel-aware:** finds off-targets reachable only via an insertion/deletion that the ungapped
  Hamming scan misses (Smith & Waterman 1981; Altschul et al. 1990). Identity threshold 0.75 over the
  probe length from Kane et al. (2000), caller-configurable.
- **On/off separation:** the first perfect ungapped full-coverage exact match is the intended
  on-target and is excluded from `OffTargetHits`; imperfect/indel hits (and extra perfect repeats)
  are off-targets. Fixes the pooling of the default `ValidateProbe.OffTargetHits` (which is unchanged).
- **Default behaviour unchanged:** `ValidateProbe` (ungapped Hamming) and `CheckSpecificity` are untouched.
- **Tests:** 9 new evidence-based tests (MG1–MG4, SG1–SG2 + 3 guards) in
  `ProbeDesigner_ProbeValidation_Tests.cs`, all with exact hand-derived values (indel hit identity 1.0;
  indel+mismatch hit identity 10/12 = 0.8333; on-target at start 5/end 16; off-target at start 27; the
  Hamming scan misses the indel site). Full fixture 28/28 pass.
- **Residual:** `ScanOffTargetsGapped` is an exhaustive sliding Smith–Waterman scan, NOT a seeded
  BLAST k-mer index over a whole genome, and has no duplex-Tm / E-value model.
- **Status reset to ☐** in the root registry for independent re-validation of the changed unit.

(Optional `hasStructure`-in-`IsValid` change was NOT made; structure remains advisory by design.)
