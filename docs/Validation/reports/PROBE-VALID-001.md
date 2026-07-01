# Validation Report: PROBE-VALID-001 — Hybridization Probe Validation

- **Validated:** 2026-06-24   **Area:** MolTools
- **Canonical method(s):** `ProbeDesigner.ValidateProbe(probeSequence, referenceSequences, maxMismatches=3, selfComplementarityThreshold=0.3)`, `ProbeDesigner.CheckSpecificity(probeSequence, genomeIndex)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_ProbeValidation_Tests.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **Latest re-validation:** 2026-06-25 (see bottom section) — ✅ CLEAN, +1 locking test (SG3), no code change.

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

---

## 2026-06-24 — Limitation fix: Karlin–Altschul E-value / bit-score (re-validation pending)

The follow-up below was implemented. The citable statistical piece of "BLAST-grade" off-target search —
the **Karlin–Altschul E-value and bit score** — was added as opt-in methods; the scan, `ValidateProbe`,
and all defaults are unchanged.

- **New methods:**
  - `ProbeDesigner.ComputeLambdaNucleotide(int match, int mismatch, double baseFrequency = 0.25)` →
    `double` λ, the unique positive root of Σ p_i p_j e^{λ s_ij} = 1, solved numerically by bisection.
  - `ProbeDesigner.ComputeKarlinAltschul(double rawScore, int queryLength, long databaseLength, ScoringMatrix? scoring = null, double k = 0.711, double baseFrequency = 0.25)`
    → `KarlinAltschulStatistics(RawScore, Lambda, K, BitScore, EValue, QueryLength, DatabaseLength)`.
- **Formulas (verbatim, Karlin & Altschul 1990 PNAS 87:2264; Altschul et al. 1990 J Mol Biol 215:403):**
  E = K·m·n·e^{−λS}; S' = (λS − ln K)/ln 2; E = m·n·2^{−S'}; λ = root of Σ p_i p_j e^{λ s_ij} = 1.
  Preconditions enforced: at least one positive score and negative expected per-pair score.
- **λ derived, not hard-coded:** the +1/−3 uniform-0.25 scheme solves to λ = 1.3740631 (asserted ≈ 1.374,
  matching the published NCBI blastn value). K is exposed as a caller parameter (its closed form needs the
  Karlin–Altschul score-lattice machinery), defaulted to the published nucleotide value 0.711.
- **Tests:** 9 new evidence-based tests (KA1–KA7) in `ProbeDesigner_ProbeValidation_Tests.cs`: λ ≈ 1.374
  to 1e-6, the root satisfies its defining equation, a hand-derived (S=30, m=20, n=1000) bit score
  59.9627 and E = 1.7802e-14, the two E-value forms agree, E decreases with score and scales linearly
  with m·n, plus precondition/argument guards. Full fixture 37/37 pass.
- **`ScanOffTargetsGapped`, `ValidateProbe`, `CheckSpecificity` and every default unchanged.**
- **Residual:** genome-scale **performance** only — no precomputed seeded k-mer index over a whole-genome
  database. The exhaustive sliding Smith–Waterman scan already finds every hit a seed would (correctness
  is complete); the seed index is purely a speed optimization. No capability-level gap remains.
- **Status stays ☐** in the root registry for independent re-validation of the changed unit.

---

## 2026-06-25 — FRESH RE-VALIDATION (independent, this session)

Independent re-validation of the changed unit in a fresh context. Scope confirmed against the
registry row + report: **PROBE-VALID-001 = `ValidateProbe` + `CheckSpecificity` + the campaign-added
`ScanOffTargetsGapped` (gapped Smith–Waterman off-target scan) + on/off-target separation**
(records `GappedProbeHit`, `GappedSpecificityResult`; helper `ScanReferenceGapped`). The
Karlin–Altschul methods (`ComputeLambdaNucleotide`, `ComputeKarlinAltschul`, tests KA1–KA12) are
the **separate** unit PROBE-EVALUE-001 (already validated CLEAN) and are out of scope here; they
were not touched. All numbers below were re-derived this session with an independent Python
Smith–Waterman, not copied from the prior report or the repo's tests.

- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN (no code change; one locking test added)

### Stage A — Description (external first sources, retrieved this session)
- **Kane MD et al. (2000) Nucleic Acids Res 28(22):4552–4557** — confirmed from the primary
  abstract/results (academic.oup.com/nar/article/28/22/4552): two specificity criteria —
  (1) *"any 'non-target' transcripts (cDNAs) >75% similar over the 50 base target may show
  cross-hybridization"* (a specific probe needs **<75 % overall identity** to non-targets); and
  (2) a 100 %-identical stretch with a non-target *"should be limited to ∼15 bases"* (>15 contiguous
  identical bases cross-hybridizes; ~15 bp ≈ 1 % signal, 20 bp ≈ 4 %). The code's
  `DefaultOffTargetMinIdentity = 0.75` is exactly the Kane criterion-1 boundary, sourced and
  caller-overridable.
- **Smith TF & Waterman MS (1981) J Mol Biol 147:195–197** and **Altschul SF et al. (1990) J Mol
  Biol 215:403–410** — standard, uncontroversial grounding for the local-alignment / gapped
  off-target scan that detects indel-reachable off-targets the ungapped Hamming scan misses.
- **Findings / divergences (why PASS-WITH-NOTES):**
  1. `ScanOffTargetsGapped` realizes **Kane criterion 1** (the 75 %-identity threshold) faithfully.
     **Kane criterion 2** (the explicit *>15 contiguous identical bases* rule) is **not** a separate
     code path — it is subsumed by the identity-over-probe-length heuristic and is honestly
     documented as such (it is not over-claimed as a distinct contiguous-stretch detector). A
     declared scope boundary, not a biological error.
  2. The default `ValidateProbe.OffTargetHits` remains the ungapped Hamming-≤-k heuristic that
     **pools** on- and off-target matches (documented, unchanged); `ScanOffTargetsGapped` is the
     corrected, on/off-separating, indel-aware path. Both are internally consistent and documented.

### Stage B — Implementation (code path + recomputed cross-checks)
- Code reviewed: `ScanOffTargetsGapped` (`ProbeDesigner.cs:1000-1052`), `ScanReferenceGapped`
  (`:1059-1134`), `ValidateProbe` (`:853-926`), `CheckSpecificity` (`:931-949`),
  `FindApproximateMatches` (`:1488-1503`). The scan slides `SequenceAligner.LocalAlign` (BlastDna
  +2/−3, gap −2) across each reference, computes `identity = identical/probeLen`,
  `coverage = ungapped/probeLen`, flags gaps, gates on `minIdentity`, collapses overlapping windows
  to one best-per-site hit, and separates the first perfect ungapped full-coverage exact match as the
  intended on-target.
- **Independent cross-check — hand-constructed probe + target set** (probe `ACGTGGCATTACGGCATTCA`,
  20 nt; references = on-target / 0.85-mismatched / low-identity / indel-only). Hand SW vs the actual
  compiled library agree **exactly**:

  | Target | Hand SW | Code (`ScanOffTargetsGapped`, minId 0.75) | Match |
  |--------|---------|-------------------------------------------|-------|
  | (a) exact on-target (ref0 [5,24]) | id 20/20=1.0, cov 1.0, no gap | OnTarget [5,24] id 1.0 cov 1.0 gap=F | ✅ |
  | (b) high-identity off-target (ref1) | `…CATTAC…`/`…CATAAC…` 17/20=**0.85**, no gap | OffTarget id 0.85 gap=F | ✅ |
  | (c) low-identity (ref2) | best local block 4/20=0.20 (<0.75) | **not called** | ✅ |
  | (d) indel-only off-target (ref3 [5,25]) | `ACGTGGCATT-ACGGCATTCA`/`…TTAACGG…` 20/20=**1.0**, 1 gap | OffTarget id 1.0 gap=T | ✅ |

  Result: OnTargets=1, OffTargetCount=2, IsSpecific=false; the on-target is excluded from off-targets;
  low-identity is correctly rejected; the indel site is flagged only via the gap.
- **Indel cross-checks re-derived** (matching MG3/MG4): probe `ACGTACGTACGT` (12 nt) —
  `ACGTAC-GTACGT` vs `ACGTACTGTACGT` → 12/12 = **identity 1.0** (MG3); SW zero-floor trims the
  mismatched tail of `ACGTACTGTACTT` → `ACGTAC-GTAC` vs `ACGTACTGTAC` → 10/12 = **0.8333** (MG4).
  Positions verified: on-target [5,16], indel site at 27; the Hamming scan (maxMismatches=3) sees
  **only** the on-target (every indel window has ≥6 mismatches) — confirming the gapped scan finds
  what the ungapped misses.
- **Variant consistency:** `CheckSpecificity` (suffix-tree exact) and `ValidateProbe` (Hamming)
  share the `0/1/N → 0.0/1.0/(1/N)` rule; `ScanOffTargetsGapped` is the separate, additive,
  on/off-separating path. `DesignProbes(genomeIndex,…)` filters via `CheckSpecificity` (integration
  test passes).
- **Test-quality audit:** assertions trace to Kane 2000 (0.75 threshold; SG1 boundary 0.75 vs 0.90),
  to hand-derived SW values (MG1–MG4 identities 1.0 and 10/12; on-target [5,16]; off-target at 27),
  and to the on/off separation contract (MG2, SG2). Guards: null probe, null refs, empty probe.
  **Coverage gap closed this session:** added **SG3**
  (`ScanOffTargetsGapped_OnTargetPlusHighIdentityAndIndelOffTargets_SeparatedAndLowIdentityExcluded`)
  locking the full Kane separation — exact on-target + a 0.85 ungapped off-target + a low-identity
  exclusion + an indel off-target in one scan, all values hand-derived. No green-washing; every
  expected value is externally sourced or hand-computed.

### Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — the 75 %-identity Kane criterion and the SW/BLAST gapped-scan +
  on/off-separation logic are correct and sourced; the explicit >15-contiguous-base Kane criterion-2
  is a documented heuristic subsumption, and the default `ValidateProbe` pooling is by-design.
- **Stage B: PASS** — code realizes the validated description exactly; all hand-recomputed
  cross-checks (independent SW) agree to the digit; all edge cases handled; tests are real, exact,
  and now cover the on-target / high-identity-off-target / low-identity-exclusion / indel-off-target
  separation. Added 1 locking test (SG3).
- **End-state: ✅ CLEAN.** No code change; the genome-scale seeded BLAST index is a documented
  performance boundary (out of scope, not a correctness gap). Full unfiltered suite green:
  **Seqeron.Genomics.Tests 18783 passed / 0 failed** (was 18782; +1 SG3); all assemblies Failed: 0;
  0 warnings on the changed test file. ProbeValidation fixture 43/43.
