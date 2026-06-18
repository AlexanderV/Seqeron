# Validation Report: PANGEN-CORE-001 — Pan-genome construction (core/accessory/unique), genome fluidity, open/closed

- **Validated:** 2026-06-15   **Area:** PanGenome (Metagenomics)
- **Canonical method(s):** `PanGenomeAnalyzer.ConstructPanGenome`, `GetCoreGeneClusters`, private `CalculateGenomeFluidity`, `DeterminePanGenomeType`/`EstimateHeapsDecayExponent`
- **Stage A verdict:** PASS-WITH-NOTES (one description defect found and fixed: core-membership formula)
- **Stage B verdict:** PASS-WITH-NOTES (one code defect found and fixed: floor core threshold; documented order-dependence of openness)
- **End-state:** CLEAN (defect completely fixed; full suite green)

## Stage A — Description

### Sources opened this session (not the repo's citations)
- **Kislyuk et al. 2011, BMC Genomics 12:32** — PMC3030549 (full text fetched). Confirmed the genomic fluidity equation **verbatim**: `φ = (2/[N(N−1)]) · Σ_{k<l} (U_k+U_l)/(M_k+M_l)`; N genomes; U_k/U_l = gene families unique to k/l; M_k/M_l = total families in k/l; range 0..1; φ=0 ⇒ share all genes, φ=1 ⇒ share none; "fluidity of 0.1 ⇒ 10% unique, 90% shared".
- **micropan `heaps()` (CRAN refman)** — fetched. Confirmed openness rule **verbatim**: "If `alpha<1.0` the pan-genome is open, if `alpha>1.0` it is closed."
- **Page et al. 2015 (Roary), PMC4817141** — fetched. Confirmed core definition **verbatim**: "Core is defined as a gene being in **at least 99% of samples**, which allows for some assembly errors in very large datasets." `-cd` default 99.
- **Roary README + bioinformatics.org search** — confirmed `-cd` = "percentage of isolates a gene must be in to be core [99]" and identity `-i` default 95 (BLASTP). The standard 4-tier scheme: core 99–100%, soft-core 95–99%, shell 15–95%, cloud <15% (a percentage/fraction test).

### Formula check
- Fluidity: matches Kislyuk exactly. ✅
- Openness: alpha<1 open / alpha>1 closed matches micropan/Tettelin. ✅
- **Core membership: DEFECT.** The description (algorithm doc, Evidence, TestSpec INV-03) stated core ⟺ `occupancy ≥ floor(coreFraction · N)`. This is **not** the Roary definition. Roary's rule is fractional — "in **at least 99% of samples**", i.e. `occupancy / N ≥ coreFraction`. The two diverge for small N: at N=3, coreFraction=0.99, floor(2.97)=2 wrongly makes a 2/3 (66.7%) cluster "core", whereas Roary requires ≥99% (3/3). `floor(coreFraction·N)` is an unsourced artifact.

### Edge-case semantics
- Empty/null → empty result (sourced as "no data"). ✅
- N=1 → occupancy 1/1 = 100% ≥ coreFraction ⇒ all core; no pairs ⇒ fluidity 0. ✅ (consistent under the corrected fractional rule too).
- N<3 → openness not estimable ⇒ Closed default (sourced, micropan needs several genomes). ✅
- Fluidity empty-pair (M_k+M_l=0) contributes 0 — reasonable neutral-element convention (not explicitly in Kislyuk; only arises for empty genomes). 🟡 documented assumption.

### Independent cross-check (numbers)
- **Hand-derived fluidity** A={c1,c2,c3}, B={c1,c2,c4}, C={c1,c5,c6}: pairs 2/6, 4/6, 4/6 ⇒ φ = (1/3)(10/6) = **10/18 = 0.5̄**. Recomputed by hand; matches M3. ✅
- **Heaps alpha** recomputed in Python: M7 new-gene curve (1,1,1) at k=2,3,4 ⇒ slope 0 ⇒ alpha 0 < 1 ⇒ **Open**. M8 curve (4,2,1) ⇒ slope −1.9809 ⇒ alpha 1.9809 > 1 ⇒ **Closed**. Both match the asserted exact `PanGenomeType`. ✅

### Findings / divergences (Stage A)
- **DEF-A1 (fixed):** core-membership formula was `floor(coreFraction·N)` instead of the sourced fractional rule `occupancy/N ≥ coreFraction`. Fixed in algorithm doc, Evidence, and TestSpec (INV-03, key-evidence point, S4 row, added S4b).

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs`
- `ConstructPanGenome` L103–168 (partition + fluidity + type)
- `GetCoreGeneClusters` L789+
- `CalculateGenomeFluidity` L371–421 (Kislyuk)
- `DeterminePanGenomeType`/`EstimateHeapsDecayExponent` L440–510 (log-log regression on dictionary order)

### Formula realised correctly?
- **Fluidity:** code computes `unique = |s1\s2| + |s2\s1|`, `total = |s1| + |s2|`, averaged over pairs with `total>0`, divided by pair count = `2/(N(N−1))·Σ`. Matches Kislyuk exactly. ✅
- **Openness:** alpha = −slope of log(newClusters) vs log(k) for k≥2; open iff alpha<1. Rule correct. The α *fit* is a single dictionary-order regression, not the permutation-averaged micropan nonlinear fit used by the separate public `FitHeapsLaw`. This is disclosed in the spec §5.2 and means borderline calls are genome-order-dependent. 🟡 PASS-WITH-NOTES.
- **Core membership: DEFECT (fixed).** Code had `coreThreshold = (int)(totalGenomes * coreFraction)` (floor) and `occupancy >= coreThreshold` in both `ConstructPanGenome` and `GetCoreGeneClusters`. Replaced with `IsCoreOccupancy`: `occupancy / N ≥ coreFraction − ε` (ε=1e-9 absorbs 0.99·100 = 98.999… round-off). Now matches Roary "at least 99% of samples".

### Cross-verification table recomputed vs code
| Case | Expected (sourced) | Code after fix |
|------|--------------------|----------------|
| M1 partition (cf 1.0) | core 1, accessory 1, unique 3, total 5 | ✅ |
| M3 fluidity | 10/18 | ✅ within 1e-10 |
| M4 identical | 0 | ✅ | 
| M5 disjoint | 1 | ✅ |
| M7 open | Open | ✅ |
| M8 closed | Closed | ✅ |
| M9 core thr 1.0 {3,2,1}/3 | {c1} only | ✅ |
| S4 cf 0.99 N=3 | core 1 (3/3 only), accessory 1, unique 1 | ✅ (was core 2 before fix) |
| S4b cf 0.99 N=100 | 99/100 core, 98/100 not | ✅ |

### Variant/delegate consistency
- `GetCoreGeneClusters` and `ConstructPanGenome` now share `IsCoreOccupancy` — consistent. ✅
- MCP `MetagenomicsTools.ConstructPanGenome`/`GetCoreGenes` delegate unchanged; benefit from the fix.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** M3 (10/18), M4 (0), M5 (1), M7/M8 (exact Open/Closed verified by independent Python regression), M9, M1, S4, S4b all assert exact externally-derived values. A deliberately-wrong floor implementation now **fails** S4/S4b (it would call 2/3 core / 98/100 core).
- **No green-washing:** the pre-existing M6 (`Is.InRange`) and C2 (property) are retained as supplementary invariant checks alongside the exact M3/M4/M5; no exact assertion was weakened. S4 strengthened from "2 core (floor)" to exact sourced "1 core" + added float-boundary S4b.
- **Coverage:** all public surface for this unit exercised — `ConstructPanGenome` (partition/fluidity/type/fractions), `GetCoreGeneClusters` (incl. empty + float boundary), null/empty/single-genome edges. 18 tests, all pass.
- **Honest green:** FULL unfiltered suite **6509 passed, 0 failed**; build **0 errors** (4 pre-existing warnings in unrelated files). ✅
- **Gate result: PASS.**

### Findings / defects (Stage B)
- **DEF-B1 (fixed):** floor core threshold (`(int)(N·coreFraction)`) in `ConstructPanGenome` L122 and `GetCoreGeneClusters` L794 → replaced with fractional `IsCoreOccupancy`. Tests S4 corrected, S4b added.
- **Note:** openness α uses single dictionary-order regression (order-dependent), disclosed simplification — not a defect.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES — core formula corrected to the sourced fractional rule; fluidity & openness already correct.
- **Stage B:** PASS-WITH-NOTES — core-threshold floor defect fixed in code + tests + description; order-dependent α fit documented.
- **End-state:** CLEAN. One genuine defect (floor vs fractional core membership) found and completely fixed across code, tests, and all three description artifacts; full suite green.
