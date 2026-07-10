---
type: source
title: "Validation report: COMPGEN-CLUSTER-001 (Conserved gene clusters — common intervals, ComparativeGenomics.FindConservedClusters)"
tags: [validation, comparative-genomics, governance]
doc_path: docs/Validation/reports/COMPGEN-CLUSTER-001.md
sources:
  - docs/Validation/reports/COMPGEN-CLUSTER-001.md
source_commit: 665dc3361ce2789ca8ede9ad2e88ea718c20310e
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: COMPGEN-CLUSTER-001

The two-stage **validation write-up** for test unit **COMPGEN-CLUSTER-001** — Conserved Gene
Clusters under the **common-interval** model of permutations (a set of ortholog-group labels that
forms a contiguous *interval* in every genome, order- and strand-free inside the window), validated
2026-06-16. This is the *report* artifact that feeds one row of the [[validation-ledger]]; it records
the validator's independent **verdict** on both the algorithm description (Stage A) and the shipped
code (Stage B). The wider campaign is [[validation-and-testing]]. The algorithm itself — parameters,
invariants, worked oracles, corner cases — is synthesized in the concept
[[conserved-gene-clusters-common-intervals]], and [[test-unit-registry]] tracks the unit. Distinct
from [[compgen-cluster-001-evidence]] — the pre-implementation evidence artifact sourced from
`docs/Evidence/` — this page is the independent two-stage re-validation verdict.

## Verdict

**Stage A: PASS · Stage B: PASS-WITH-NOTES · End state: CLEAN.** No implementation defect and no
code change. Stage B is PASS-WITH-NOTES purely because **three weak test assertions were
strengthened in-session** (a test-quality defect, not a code defect): permissive
`Contains.Item` / `Does.Not.Contain` checks — where the exact set is known and brute-forceable —
were rewritten to exact sourced sets, closing a green-wash gap. The full unfiltered suite ran
**6605 passed, 0 failed** (1 benchmark skipped as designed); `dotnet build` 0 errors (4 pre-existing
warnings, all in unrelated untouched files). Test-quality gate: **PASS** (after the M3/M5/S3 fixes).

## Scope clarification

The orchestrator prompt described this unit generically as "ortholog/gene clustering … COG/OrthoMCL
… connected-components / MCL." The unit's Registry row, TestSpec and Evidence all unambiguously
define COMPGEN-CLUSTER-001 as the **common-interval gene-cluster model** — the single public method
under test is `FindConservedClusters`. COG/OrthoMCL ortholog *grouping* is a different concept,
covered by [[ortholog-detection-reciprocal-best-hits]] and its units (COMPGEN-ORTHO-001,
COMPGEN-RBH-001). This report therefore validates the common-interval model, which is what the code
actually computes.

## Canonical method & source under test

Canonical method `ComparativeGenomics.FindConservedClusters(genomes, orthologGroups,
minClusterSize=3, maxGap=2)`. Code path
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:914–1021`
(`FindConservedClusters` + helper `IsIntervalOf`):

- Maps each genome to its ordered ortholog-group label sequence (genes with no group → a
  window-breaking sentinel).
- **Candidate generation (`:948–965`):** the group-set of every foreign-free contiguous window of
  genome 0. Complete because every common interval is in particular an interval of genome 0, hence
  equals `Set(T[i..j])` of some genome-0 window, and that window cannot contain a foreign group.
- **Filter (`:967–978`):** keep candidates that are `IsIntervalOf` *every* genome. `IsIntervalOf`
  (`:1000–1021`) finds a window whose label set equals the target exactly (breaks on any foreign
  element) — the `Set(T[i..j])` sequence-interval test, matching Didier et al. exactly, paralogs
  included.
- **Size threshold** `effectiveMin = max(minClusterSize, 2)` (`:929`, `MinCommonIntervalSize=2`)
  enforces the `i < j` ≥ 2 rule; **< 2 genomes → empty** (`:926`); determinism by size then
  lexicographic joined labels (`:980–986`); null guards `ThrowIfNull` on both args (`:920–921`).
- Tests: `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindConservedClusters_Tests.cs`
  (12 tests). Old weak/duplicate tests were already removed from `ComparativeGenomicsTests.cs`.

## Stage A — description (algorithm faithfulness)

Confirmed against four sources, all fetched fresh this session as PDFs and extracted with
`pdftotext`: **Bui-Xuan, Habib & Paul 2013** (arXiv:1304.5140) — interval `[i,j]` defined only for
`1 ≤ i < j ≤ n` (size ≥ 2, singletons excluded); Definition 1 (attributed to Uno & Yagiura): a
common interval is a set that is an interval of **each** `P_k`; **Didier, Schmidt, Stoye & Tsur
2013** (arXiv:1310.4290) — interval of a sequence `I = Set(T[i..j])`, common interval of two
sequences = interval of both (handles paralogs/duplicates); **Uno & Yagiura 2000** and **Heber &
Stoye 2001** — provenance of the K=2 originating model and the k-permutation generalisation (full
text behind Springer auth; bibliographic + definitional facts confirmed via the two citing arXiv
papers).

The documented model — a conserved cluster = a common interval = a label set occupying a contiguous
foreign-free window in *every* genome, size ≥ 2, whole set always trivially common — matches
Source 1 Def. 1 and Source 2 §Prelim/Def. 1 **exactly**. Edge-case semantics all confirmed: **< 2
genomes** vacuous → empty; **repeated labels (paralogs)** — any matching window in each genome
counts; **a foreign group inside the window** breaks the interval. The `maxGap` parameter is
documented (Assumption 1) as **API-shape only** — the validated contract is the strict, gap-free
common-interval model. **No divergences.**

**Independent cross-check** (brute force over all subsets, contiguity = positions span exactly
`|S|−1`, none a code echo): the two-permutation golden vector `P1=Id₇`, `P2=(7 2 1 3 6 4 5)` min 2 →
`{1,2}`, `{1,2,3}`, `{3,4,5,6}`, `{4,5}`, `{4,5,6}`, `{1..6}`, `{1..7}` (7 sets) **= paper Example
1**; `{2,3}` not an interval of `P2` (positions 1,3); Didier `T,S`: `{1,2}` No / `{1,2,3,4}` Yes
(matches Source 2 Example 1); plus M3/M5/S2/S3/C2 sets all brute-forced.

## Stage B — implementation

The formula is realised correctly (candidate completeness + interval filter + size threshold + < 2
→ empty + determinism + null guards, all traced above). The strengthened exact-set tests now assert
the independently brute-forced sets, and the full suite passes, so the code produces **exactly**
those sets for every dataset (not a superset/subset) — the strongest available confirmation.
Single public method; no `*Fast`/instance variant; the MCP wrapper (`AnalysisTools.cs`) delegates
to the same method.

**Test-quality audit (HARD gate).** M1/M4/C2/S2 lock the *full* exact set against the paper's
Example 1 or an independent brute force; P1 is a property test (every returned cluster contiguous in
every genome) over a seeded random input, guarding over-reporting. Three green-wash-risk assertions
were rewritten this session against the known/computable exact sets:

- **M3** (`RepeatedLabels_WindowStillMatches`): `Contains.Item({a,b,c})` →
  `Is.EquivalentTo({a,b,c},{a,b,c,x},{b,c,x})`.
- **M5** (`ForeignGroupInsideWindow_BreaksCluster`): `Does.Not.Contain({a,b,c})` → `Is.Empty`
  (no size-≥3 common interval).
- **S3** (`ConservedInTwoButSplitInThird`): `Does.Not.Contain({1,2})` + `Contains.Item({2,3,4})` →
  `Is.EquivalentTo({2,3},{3,4},{2,3,4})`.

These now fail against an over- or under-reporting implementation. Coverage exercises the single
public method and every Stage-A branch: golden vector (M1), one-genome-only negative (M2),
paralog/repeated-label path (M3), size filter (M4), foreign-group window break (M5), < 2 genomes
empty (S1), k=3 all-conserved (S2), k=3 split-in-one (S3), determinism (C1), trivial-whole-set-only
(C2), null args (C3), property/INV-01 (P1).

## Findings

- **No code defect and no code change (End state CLEAN).** Description, formulas, definitions, edge
  cases and worked example all match the cited primary literature verbatim; every expected value
  traces to the paper's Example 1 or an independent brute force run this session.
- **Test-quality defect fixed in-session:** three weak assertions (M3/M5/S3) rewritten from
  permissive `Contains`/`Does.Not.Contain` to exact brute-forced sets, verified green — this is the
  sole reason Stage B is PASS-WITH-NOTES rather than PASS.
- **Documented assumption (not a defect):** the public method keeps a `maxGap` parameter but the
  validated behaviour is the strict, gap-free common-interval model; `maxGap` does not relax it, and
  the gene-teams gapped extension (Bergeron, Corteel & Raffinot 2002) is not implemented — a design
  decision carried on [[conserved-gene-clusters-common-intervals]].
