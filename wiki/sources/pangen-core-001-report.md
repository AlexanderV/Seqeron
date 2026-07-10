---
type: source
title: "Validation report: PANGEN-CORE-001 (pan-genome core/accessory/unique partition + genomic fluidity + open/closed)"
tags: [validation, pan-genome, governance]
doc_path: docs/Validation/reports/PANGEN-CORE-001.md
sources:
  - docs/Validation/reports/PANGEN-CORE-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: PANGEN-CORE-001

The two-stage **validation write-up** for test unit **PANGEN-CORE-001** (pan-genome
construction — core / accessory / unique partition by cluster occupancy, genomic fluidity,
and open/closed classification), validated 2026-06-15. This is the *report* artifact that
feeds one row of the [[validation-ledger]]; it records the validator's **verdict** on both the
algorithm description and the shipped code. The two-stage methodology is the
[[validation-protocol]]; the algorithm itself is summarized in
[[pan-genome-core-accessory-partition]]. Distinct from the pre-implementation
[[pangen-core-001-evidence]] artifact.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS-WITH-NOTES · State: CLEAN.** One genuine defect —
the **core-membership formula** (floor `floor(coreFraction·N)` instead of the sourced fractional
`occupancy/N ≥ coreFraction`) — was found and **completely fixed** across code, tests, and all
three description artifacts. The fluidity and openness formulas were already correct. **18 tests
pass**; full unfiltered suite **6509 passed / 0 failed**, build 0 errors (4 pre-existing warnings
in unrelated files). The order-dependent Heaps α fit is a disclosed simplification, not a defect.

## Stage A — description (algorithm faithfulness)

- Canonical methods: `PanGenomeAnalyzer.ConstructPanGenome`, `GetCoreGeneClusters`, private
  `CalculateGenomeFluidity`, `DeterminePanGenomeType` / `EstimateHeapsDecayExponent`.
- Sources opened this session (external, first-source): **Kislyuk et al. 2011** (PMC3030549) —
  fluidity `φ = (2/[N(N−1)])·Σ_{k<l}(U_k+U_l)/(M_k+M_l)`, range 0..1, confirmed **verbatim**;
  **micropan `heaps()`** — openness rule "alpha<1.0 open, alpha>1.0 closed" verbatim; **Page et al.
  2015 (Roary)** — core = "in **at least 99% of samples**", `-cd` default 99, standard tiers
  core 99–100 / soft-core 95–99 / shell 15–95 / cloud <15.
- Formula check: fluidity matches Kislyuk exactly ✅; openness matches micropan/Tettelin ✅.
- **Core membership: DEFECT.** The description (algorithm doc, Evidence, TestSpec INV-03) stated
  core ⟺ `occupancy ≥ floor(coreFraction·N)` — not the Roary fractional rule
  `occupancy/N ≥ coreFraction`. They diverge for small N: at N=3, coreFraction=0.99,
  `floor(2.97)=2` wrongly makes a 2/3 (66.7%) cluster core. `floor(coreFraction·N)` is an unsourced
  artifact. **DEF-A1 fixed** in the algorithm doc, Evidence, and TestSpec (INV-03, key evidence, S4
  row, added S4b).
- Independent cross-check (numbers): hand-derived fluidity A/B/C → `φ = 10/18 = 0.5̄` ✅; Heaps α
  recomputed in Python — curve (1,1,1) ⇒ slope 0 ⇒ α 0 < 1 ⇒ **Open**, curve (4,2,1) ⇒ slope
  −1.9809 ⇒ α 1.9809 > 1 ⇒ **Closed** ✅.

## Stage B — implementation (code review + cross-check)

- Code path: `PanGenomeAnalyzer.cs` — `ConstructPanGenome` L103–168, `GetCoreGeneClusters` L789+,
  `CalculateGenomeFluidity` L371–421, `DeterminePanGenomeType` / `EstimateHeapsDecayExponent`
  L440–510.
- Fluidity: code computes `unique = |s1\s2| + |s2\s1|`, `total = |s1| + |s2|`, averaged over pairs
  with `total>0` = `2/(N(N−1))·Σ` — matches Kislyuk exactly ✅. Openness: `alpha = −slope` of
  `log(newClusters)` vs `log(k)` for k≥2, open iff α<1 — rule correct ✅ (the single dictionary-order
  regression, not the permutation-averaged micropan fit, is a disclosed §5.2 simplification making
  borderline calls order-dependent — PASS-WITH-NOTES).
- **Core membership: DEFECT (fixed).** Code had `coreThreshold = (int)(totalGenomes*coreFraction)`
  (floor) with `occupancy >= coreThreshold` in both `ConstructPanGenome` and `GetCoreGeneClusters`.
  Replaced with shared `IsCoreOccupancy`: `occupancy/N ≥ coreFraction − ε` (ε=1e-9 absorbs
  0.99·100 = 98.999… round-off). Now matches Roary "at least 99% of samples". **DEF-B1 fixed**;
  tests S4 corrected, S4b added.
- Cross-verification table recomputed vs code (post-fix): M1 partition (cf 1.0) core 1 / accessory 1
  / unique 3 ✅; M3 fluidity 10/18 ✅ within 1e-10; M4 identical → 0 ✅; M5 disjoint → 1 ✅; M7 Open ✅;
  M8 Closed ✅; M9 core thr 1.0 → {c1} ✅; **S4 cf 0.99 N=3 → core 1 (3/3 only)** ✅ (was core 2 before
  fix); S4b cf 0.99 N=100 → 99/100 core, 98/100 not ✅.
- Test-quality audit (HARD gate) PASS: M3/M4/M5/M7/M8/M9/M1/S4/S4b assert exact externally-derived
  values (a deliberately-wrong floor implementation now **fails** S4/S4b); no green-washing — the
  pre-existing M6 `Is.InRange` and C2 property are retained as supplementary invariants alongside
  the exact asserts, S4 strengthened from "2 core (floor)" to exact sourced "1 core" + float-boundary
  S4b; all public surface exercised; 18 tests pass.

## Findings

- **One genuine defect (floor vs fractional core membership), found and completely fixed** across
  code (`IsCoreOccupancy` shared by `ConstructPanGenome` L122 and `GetCoreGeneClusters` L794), tests
  (S4 corrected, S4b added), and all three description artifacts. Full suite green (6509/0).
- **Note (not a defect):** openness α uses a single dictionary-order regression (order-dependent),
  a disclosed spec §5.2 simplification distinct from the permutation-averaged public `FitHeapsLaw`.
- **Note (documented assumption):** empty-pair fluidity (`M_k+M_l=0`) contributes 0 as a neutral
  element; only arises for empty genomes.
- MCP `MetagenomicsTools.ConstructPanGenome` / `GetCoreGenes` delegate unchanged; benefit from the
  fix.

See the full report at `docs/Validation/reports/PANGEN-CORE-001.md`.
