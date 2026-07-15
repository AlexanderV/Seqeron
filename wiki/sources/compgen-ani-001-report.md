---
type: source
title: "Validation report: COMPGEN-ANI-001 (Average Nucleotide Identity — ANIb, ComparativeGenomics.CalculateANI / CalculateReciprocalAni)"
tags: [validation, comparative-genomics, governance]
doc_path: docs/Validation/reports/COMPGEN-ANI-001.md
sources:
  - docs/Validation/reports/COMPGEN-ANI-001.md
source_commit: 205b259dc3168dfda72a89caf5103f39ac5e1ce9
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: COMPGEN-ANI-001

The two-stage **validation write-up** for test unit **COMPGEN-ANI-001** — Average Nucleotide
Identity under the **ANIb** definition (Goris et al. 2007), validated 2026-06-24. This is the
*report* artifact that feeds one row of the [[validation-ledger]]; it records the validator's
independent **verdict** on both the algorithm description (Stage A) and the shipped code
(Stage B). The wider campaign is [[validation-and-testing]]. The algorithm itself — parameters,
invariants, worked oracles, corner cases — is synthesized in the concept
[[average-nucleotide-identity]], and [[test-unit-registry]] tracks the unit. Distinct from
[[compgen-ani-001-evidence]] — the pre-implementation evidence artifact sourced from
`docs/Evidence/` — this page is the independent two-stage re-validation verdict.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: CLEAN.** No defect found; no code or test change
required. This re-validation covers the limitations-campaign change in commit `69c51fa0`, which
added **gapped fragment alignment** (reusing `SequenceAligner.LocalAlign`) and **reciprocal
(two-way) ANI**. The prior report (2026-06-16) had graded Stage B **PASS-WITH-NOTES** because the
implementation was ungapped / single-direction and the `minAlignableFraction` parameter was
structurally inert; **both notes are now resolved**. The filtered ANI class ran **20/20 pass**;
the broader `~Comparative` filter (fuzz / metamorphic / algebraic / combinatorial / property
families) ran **480/480 pass**; build 0 Error(s) / 0 Warning(s).

## Canonical methods & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:1076–1228`:

- `CalculateANI(query, reference, fragmentLength, minIdentity, minAlignableFraction, gapped)`
  — null/empty → 0; `fragmentLength ≤ 0` → `ArgumentOutOfRangeException`; fragments the query
  consecutively (`start += fragmentLength`, trailing partial dropped), places each fragment
  (gapped or ungapped), keeps fragments with `identity > minIdentity && alignableFraction ≥
  minAlignableFraction`, and returns the mean of the kept identities (else 0).
- `CalculateReciprocalAni(genomeA, genomeB, …)` (`:1126–1144`) — returns `(CalculateANI(A,B,…) +
  CalculateANI(B,A,…)) / 2.0`; order-independent by construction.
- private `BestGappedFragmentMatch` (`:1194–1225`) — full Smith-Waterman via
  `SequenceAligner.LocalAlign` with `BlastDna` scoring; counts `identicalColumns` (gaps skipped)
  and `ungappedColumns`, giving `identity = identicalColumns/fragLen` and `alignableFraction =
  ungappedColumns/fragLen` — exactly pyani's `ani_pid = ani_alnids/qlen` and `ani_coverage =
  (blast_alnlen − blast_gaps)/qlen`.
- private `BestUngappedFragmentMatch` (`:1155–1182`) — best-offset matching-base scan (unchanged);
  coverage 1.0 when ref ≥ fragment, else (0,0).
- Tests: `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CalculateANI_Tests.cs` (20 tests).

## Stage A — description (algorithm faithfulness)

Confirmed against three retrieved sources: **Goris et al. 2007** (IJSEM 57:81–91,
DOI 10.1099/ijs.0.64483-0) — query cut into **consecutive 1020 nt fragments**, each BLASTN-searched
against the reference, ANI = "mean identity of all BLASTN matches that showed **more than 30 %**
overall sequence identity (recalculated to an identity along the entire sequence) over an
alignable region of **at least 70 %** of their length", **reverse searching** performed "to
provide reciprocal values", BLASTN gapped by default so the best match recovers indels, and
ANI ≈ 95 % ↔ 70 % DDH species boundary; **Konstantinidis & Tiedje 2005** (PNAS 102:2567–2572) —
the ≈94–96 % species boundary; **pyani** ANIb source (`anib.html`) — the exact per-fragment
conventions (`ani_pid = ani_alnids/qlen`, `ani_coverage = ani_alnlen/qlen`, filter
`cov > 0.7 & pid > 0.3`, gapped BLASTN `-xdrop_gap_final 150`, `-max_target_seqs 1`).

The documented model (identity and coverage both recalculated over the full query-fragment length;
strict `> 30 %` identity; `≥ 70 %` coverage; reverse-search averaging; defaults `L=1020`, `0.30`,
`0.70` traced to Goris) matches Goris and pyani **exactly**. Edge-case semantics all confirmed:
trailing partial fragment dropped; non-conserved fragments discarded (not scored as 0); empty /
no-qualifying-fragment → 0; direction-dependence acknowledged and symmetry restored only by
reciprocal averaging; gaps excluded from both identity and coverage. **No divergences.** The only
out-of-scope item — the OrthoANI / FastANI orthologous-best-pair variant — is **honestly disclosed
as not implemented**, not a defect.

## Stage B — implementation

Every changed-path value was hand-derived and reproduced by the code:

- **G2 indel recovery** — `AAAACCCC` (8 nt) vs `AAAATCCCC` (one inserted T): ungapped best window
  7/8 = **0.875**; gapped `AAAA-CCCC` scores `8·(+2) + 1 gap(−5−2) = +9`, beating any ungapped
  local, giving 8 identical columns → **gapped ANI 1.0 > ungapped 0.875** (INV-6).
- **G3 coverage cut-off** — `AAAACCCC` vs `AA`: best local covers ≤ 2 query columns → `cov = 2/8 =
  0.25 < 0.70` → excluded → **0.0** (this exercises the previously inert `minAlignableFraction`).
- **G4 identity cut-off (gapped)** — `AAAACGTC` vs `AAAAAAAA` fragLen 4: `CGTC` → 0 identical → `id
  = 0` not `> 0.30` → excluded; surviving `AAAA` → **1.0**.
- **R3 reciprocal = mean of directions** — A=`AAAACGTC`, B=`AAAAAAAA`, fragLen 4: A→B = 1.0, B→A =
  1.0 → **(1.0+1.0)/2 = 1.0** (INV-7, INV-8). R1/R2 confirm identical→1.0 and order-independent
  symmetry.

The full cross-verification table (G1–G4, R1–R3, plus M1–M8/S1/S2/C1 from the prior ungapped
report) matches the hand-derived values. Numerical robustness holds — identity/coverage bounded by
`Math.Min(…, 1.0)`, denominators (`fragLen`) positive by validation, empty / ref-shorter → (0,0).
Two public methods with no divergent logic (`CalculateReciprocalAni` delegates per direction; the
MCP wrapper forwards to `CalculateANI`; no `*Fast` variants).

**Test-quality audit (HARD gate): PASS.** Every deterministic value is hand-derived from
Goris/pyani, not code-echoed; G2 contrasts gapped 1.0 vs ungapped 0.875 (a real indel discriminator);
R3 ties reciprocal to the mean of two directions; M4 would expose an average-all-fragments bug
(would return 0.5). Exact `Within(1e-10)` equality on exact fractions; the sole `Is.InRange` (S1)
is the range invariant and the sole `Is.GreaterThan` (G2) is an additive assertion atop exact
equalities — no green-washing, no skips/ignores/widened tolerances. Branch coverage spans gapped
vs ungapped, both cut-offs, fragmentation/trailing, null/empty, invalid `fragmentLength`, and
reciprocal symmetry & mean.

## Findings

- **No code defect and no test change (State CLEAN).** Every worked example and cross-check
  reproduced exactly.
- **Both prior PASS-WITH-NOTES items resolved:** (1) `minAlignableFraction` is now genuinely active
  under the gapped path (G3 exercises coverage 0.25 < 0.70 exclusion); (2) the ungapped /
  single-direction simplification is replaced by gapped placement + reciprocal averaging matching
  the full ANIb definition.
- **Documented out-of-scope (not a defect):** the OrthoANI / FastANI orthologous-best-pair variant
  remains an honestly disclosed extension. The engine choice — Seqeron's own Smith-Waterman
  (`SequenceAligner.LocalAlign`, BLAST DNA scoring) rather than NCBI BLASTN — is a design decision
  carried on [[average-nucleotide-identity]], not a correctness gap.
</content>
</invoke>
