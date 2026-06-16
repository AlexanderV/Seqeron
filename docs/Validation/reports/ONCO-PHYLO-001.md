# Validation Report: ONCO-PHYLO-001 — Tumor Phylogeny Reconstruction (clonal tree from CCF clusters)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.ReconstructPhylogeny(IReadOnlyList<CcfCluster>, double tolerance)`, `OncologyAnalyzer.IdentifyTrunkMutations(ClonalPhylogeny)`, `OncologyAnalyzer.IdentifyBranchMutations(ClonalPhylogeny)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm

1. **Popic V et al. (2015), LICHeE, *Genome Biology* 16:91** — retrieved the full-text PDF from
   https://arxiv.org/pdf/1412.8574 this session, extracted text with `pypdf`, and quoted the
   constraint passages directly (not from the repo Evidence doc). Confirms:
   - The three perfect-phylogeny ordering constraints, **verbatim**: "(1) a mutation present in a
     given set of samples cannot be a successor of a mutation that is present in a smaller subset of
     these samples …; (2) a given mutation cannot have a VAF higher than that of its predecessor
     mutation (except due to CNVs) …; (3) the sum of the VAFs of mutations disjointly present in
     distinct subclones cannot exceed the VAF of a common predecessor mutation present in these
     subclones."
   - **Edge rule, Eq. 2 (verbatim):** "an edge (u → v) is added only if the nodes satisfy the
     following two constraints ∀ i ∈ samples … (1) `u.VAFi ≥ v.VAFi − ϵuv` and (2) if `u.VAFi = 0,
     v.VAFi = 0`, where ϵuv is the VAF noise error margin."
   - **Sum rule, Eq. 5 (verbatim):** "∀ nodes u∈T : ∀ i∈ samples : `Σ_{v s.t. (u→v)∈T} v.VAFi ≤
     u.VAFi + ϵ`. That is, the sum of the VAF centroids of all the children must not exceed the
     centroid of the parent. … we use inequality here since our method does not require all the true
     lineage branches to have been observed."
   - **Trunk = structural concept:** LICHeE's "trunk branch" is the first edge from the germline
     root and may contain private/partially-shared mutations ("The trunk branch … includes all shared
     mutations as well as one private and two partially shared mutations"), confirming a *topological*
     (root→first-branch-point) trunk definition rather than a strict CCF=1 criterion.

2. **Zheng L et al. (2022), PICTograph, *Bioinformatics* 38(15):3677–3683** — fetched the PMC
   full text (https://pmc.ncbi.nlm.nih.gov/articles/PMC9344857/) this session. Confirms the CCF form
   of both rules:
   - **Sum condition (verbatim):** "Trees … where the sum of the CCFs of the descendants of a parent
     node exceeded the parent node's CCF by more than ε₂ were excluded." (default ε₂ = 0.2).
   - **Lineage precedence (verbatim):** "the descendant cluster must be present in the same sample or
     a subset of samples in which the ancestral cluster is present and … the CCF of the descendant
     cluster are at most ε₁ greater than that of the ancestral cluster in each sample." (default
     ε₁ = 0.1).

3. **Trunk/branch biology cross-check** — general literature (e.g. PMC5558210, PMC7867630) confirms
   truncal = clonal mutations in essentially every tumor cell; branch = subclonal (subset of cells).
   The implementation uses LICHeE's *structural* trunk-branch notion (see Notes).

### Formula check
- Lineage precedence in the spec/doc (`u.CCF[i] ≥ v.CCF[i] − ε` and `u=0 ⇒ v=0`) matches LICHeE
  Eq. 2 exactly (VAF→CCF substitution is the standard CCF analogue, confirmed by PICTograph).
- Sum rule (`Σ_children v.CCF[i] ≤ u.CCF[i] + ε`) matches LICHeE Eq. 5 and PICTograph's sum condition
  exactly, including the inequality (not equality) rationale and the per-node, per-sample quantifier.

### Edge-case semantics check
- Inequality (not equality) tolerates unobserved branches/noise — sourced (Eq. 5 prose).
- Presence pattern (`parent absent ⇒ child absent`) — sourced (Eq. 2(2); PICTograph subset rule).
- Private/single-sample clusters are under-constrained → deterministic tie-break required — sourced
  (LICHeE Methods/Results; the deepest-valid-ancestor + id tie-break is a *documented assumption*,
  flagged as ASSUMPTION in spec §6 and doc ASM-03, and it only selects among already-valid trees).
- Empty/single-cluster/duplicate-id/NaN/out-of-range — defined in the contract.

### Independent cross-check (numbers)
Re-derived the expected topologies directly from the verbatim equations with a standalone Python
checker (no repo code):
- **M2** (A=[1,1], B=[0.6,0], C=[0,0.7]): lineage check gives A→B valid, A→C valid, **B→C invalid,
  C→B invalid** (presence pattern). Therefore B and C are *forced* siblings under A — not a tie-break
  artifact. Sum under A: s1 0.6≤1, s2 0.7≤1. Matches the test.
- **M3** (A=1.0, B=0.6, C=0.6 single sample): all pairwise lineage valid, but B+C=1.2 > 1.0 violates
  the sum rule, so they **cannot both** be A's children → a chain is forced. Matches the test.
- **M1** (1.0,0.6,0.3): each child ≤ previous, deepest-valid-ancestor → caterpillar root→A→B→C.

### Findings / divergences
None at the description level. All formulae and edge semantics trace to retrieved primary sources.

**Stage A verdict: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`
- `ReconstructPhylogeny` (6115–6215), `IdentifyTrunkMutations` (6225–6250),
  `IdentifyBranchMutations` (6258–6276), `SatisfiesLineagePrecedence` (6282–6299),
  `FitsSumRule` (6305–6316), `ValidateAndGetSampleCount` (6357–6400).

### Formula realised correctly? (evidence)
- `SatisfiesLineagePrecedence`: rejects when `parentCcf[i] < childCcf[i] - tolerance` (i.e. requires
  `parent ≥ child − ε`) and when `parent ≤ 0 && child > 0` (presence pattern). Exactly LICHeE Eq. 2.
- `FitsSumRule` enforced via a per-node, per-sample remaining-budget that is debited by each attached
  child; admission requires `budget[i] ≥ childCcf[i] − ε`. This is precisely Eq. 5
  (`Σ_children ≤ parent + ε`) realised incrementally. Verified against the green suite.
- Synthetic root has CCF=1 in every sample (`RootCcf=1.0`) and a unique id below the min cluster id;
  spanning-tree invariant (every cluster gets exactly one parent) holds.

### Cross-verification table recomputed vs code
| Case | Independent-source expectation | Code output (suite) | Match |
|------|-------------------------------|---------------------|-------|
| M1 chain (1.0,0.6,0.3) | root→A→B→C | same | ✅ |
| M2 branch ([1,1],[.6,0],[0,.7]) | root→A; A→B; A→C (B,C forced siblings) | same | ✅ |
| M3 equal subclones (1.0,0.6,0.6) | sum-rule forces root→A→B→C | same | ✅ |
| M8 single | root→A only | same | ✅ |
| S2 ε=0.1 ([.5,.5],[.55,0]) | B nests under A (0.5 ≥ 0.55−0.1) | same | ✅ |
| S3 ε=0 same input | B attaches to root (0.5 < 0.55 rejects A) | same | ✅ |
| INV-1 / INV-2 | hold on M1/M2 and 50 random valid inputs | hold | ✅ |

### Variant/delegate consistency
`IdentifyBranchMutations` = `Clusters \ IdentifyTrunkMutations` (disjoint partition, INV-4). Trunk is
the single-child path from the root. Consistent with LICHeE's structural trunk-branch (see Notes).

### Test quality audit (HARD gate)
- **Sourced expectations, not code echoes:** M1–M3, M8, S2, S3 assert exact parent ids/edges that are
  *forced* by the verbatim Eq. 2 / Eq. 5 (re-derived independently this session), not copied from
  output. A deliberately-wrong implementation (e.g. dropping the sum rule, or flipping ancestor≥desc)
  would fail M3/M2/S3. Confirmed by hand re-derivation. ✓
- **No green-washing:** assertions use exact `Is.EqualTo` on parents and edge sets; invariant tests
  (M6, M7) appropriately use `GreaterThanOrEqualTo`/`LessThanOrEqualTo` because they verify the
  *inequality invariants themselves* (this is the correct, sourced form for INV-1/INV-2, not a
  weakening). No skipped/ignored/commented tests; no widened tolerances (ε=0 strict path tested in
  S3). ✓
- **Coverage:** all three public methods exercised; both ε-branches (S2/S3); empty (S1), single (M8),
  chain (M1), branching (M2), sum-rule chain (M3), determinism (C1), and the full validation matrix
  (null, NaN, out-of-range, ragged, duplicate id, negative tolerance) with the precise exception
  subtypes. Stage-A edge cases all covered. ✓
- **Honest green:** full unfiltered suite = **6677 passed, 0 failed, 0 skipped-fail**; `dotnet build`
  0 errors (4 pre-existing warnings in unrelated `ApproximateMatcher_EditDistance_Tests.cs`, not
  touched). ✓

**Test-quality gate: PASS.**

### Findings / defects
No defects. One documented divergence (Note 1) drives the PASS-WITH-NOTES verdict.

- **Note 1 (BY-DESIGN, sourced):** `IdentifyTrunkMutations` defines the trunk *topologically* — the
  maximal single-child path from the root to the first branch point — rather than by the strict
  biological criterion "CCF≈1 in all samples / present in every cell." Consequence: on a pure chain
  (M1: 1.0→0.6→0.3) it returns trunk={A,B,C}, even though only A is strictly clonal. This is
  consistent with LICHeE's own structural "trunk branch" (which may contain private/partially-shared
  mutations) and is explicitly documented in the spec (§trunk) and algorithm doc; it is internally
  consistent and sourced, so it is not a defect. Recorded for transparency.
- **Note 2 (BY-DESIGN):** default ε=0 is stricter than the source defaults (LICHeE ϵ; PICTograph
  ε₁=0.1, ε₂=0.2); exposed via `tolerance`. Setting ε>0 only widens admissibility. Documented.
- **Note 3 (BY-DESIGN):** deepest-valid-ancestor + ascending-id tie-break is an assumption used only
  to return one deterministic tree among the constraint-valid set; it never admits an invalid edge
  (INV-1/INV-2 always hold). Documented (ASM-03).

## Verdict & follow-ups
- **Stage A: PASS** — all formulae and edge semantics confirmed against the retrieved primary sources
  (LICHeE Eq. 2 / Eq. 5 verbatim; PICTograph CCF-form sum condition + lineage precedence).
- **Stage B: PASS-WITH-NOTES** — code faithfully realises Eq. 2 and Eq. 5; tests are real, sourced,
  and complete; the only divergence (structural trunk definition, Note 1) is sourced and documented.
- **End-state: CLEAN** — no defect required fixing; full suite green (6677/0).
- No code or test changes were necessary.
