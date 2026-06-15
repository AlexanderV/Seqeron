# Validation Report: PANGEN-CLUSTER-001 — Gene Clustering (homolog grouping by global sequence identity, CD-HIT greedy model)

- **Validated:** 2026-06-15   **Area:** PanGenome
- **Canonical method(s):** `PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold)`; `CalculateSequenceIdentity` (private, tested indirectly); `CreatePresenceAbsenceMatrix` (delegate)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent retrieval)

| Source | URL | What it confirms (verbatim) |
|--------|-----|------------------------------|
| CD-HIT Algorithm wiki | https://github.com/weizhongli/cdhit/wiki/1.-Algorithm | "sorts the input sequences from long to short, and processes them sequentially from the longest to the shortest. The first sequence is automatically classified as the first cluster representative sequence. Then each query sequence … is compared to the representative sequences found before it, and is classified as redundant or representative…" and "In default manner (fast mode), a query is grouped into the first representative without comparing to other representatives." |
| CD-HIT User's Guide (wiki source) | https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki | Global sequence identity = "number of identical amino acids in alignment divided by the full length of the shorter sequence"; "-c sequence identity threshold, default 0.9"; "-G use global sequence identity, default 1 … if set to 0, then use local sequence identity, calculated as: number of identical amino acids in alignment divided by the length of the alignment" |
| CD-HIT User's Guide (vcru mirror) | https://vcru.wisc.edu/simonlab/bioinformatics/programs/cd-hit/cdhit-user-guide.pdf | corroborates the `-c`/`-G` definitions above (cross-mirror) |

All three retrieved this session; the numbers below trace to them, not to the repo's own artifacts.

### Formula check
- **Global identity** = identical residues (ungapped positional alignment) / length of the **shorter** sequence. Matches CD-HIT `-G 1` default verbatim. Algorithm doc §2.2 eq. (1) is correct.
- **Greedy incremental clustering**: long→short sort, longest = first representative, query joins the **first** representative meeting the cutoff (fast mode), else new representative. Matches the wiki verbatim.
- **Inclusive cutoff `>=`, default 0.9**: matches `-c` default.
- **Representative = longest member** (consequence of long→short ordering): correct.

### Edge-case semantics check
- null genomes → empty (no throw): documented contract, sibling-consistent with `ConstructPanGenome`. Not a CD-HIT property but an explicit local API contract; acceptable.
- empty/empty → identity 1.0; empty/non-empty → 0.0: a defined convention (0 differences over 0 positions vs no shared residues), consistent with percent-identity numerator semantics. Documented.
- threshold 1.0 → only exact-over-shorter-length sequences cluster: follows from `>=` and the shorter-length denominator. Correct.
- singleton AverageIdentity = 1.0 (self-identity): correct.

### Independent cross-check (hand-computed against the sourced formula)
- `ATGCATGC` vs `ATGCATGC` → 8/8 = **1.0**
- `ATGCATGC` vs `ATGCATGG` → 7/8 = **0.875** (single substitution at last position)
- `ATGCATGC` vs `ATGCATGCAAAA` → prefix 8/8 over shorter(8) = **1.0**
- `ATGC` vs `CGTA` → 0/4 = **0.0**
- **M6 greedy trace** (threshold 0.8): stable long→short = Q2(12),R(10),Q1(10),Q3(10). Q2→rep; R: id(R,Q2)=10/10=1.0 join; Q1: id(Q1,Q2)=9/10=0.9 join; Q3: id(Q3,Q2)=0/10=0.0 → new. **2 clusters {Q2,R,Q1},{Q3}**, big-cluster consensus = `AAAAAAAAAAAA`. Matches the implementation and tests.

### Findings / divergences
- ASM-01 (ungapped, no internal indels) and ASM-02 (homolog not ortholog; no paralog/synteny split) are clearly documented as intentional simplifications relative to full CD-HIT/Roary. Acceptable scope; recorded.
- Minor doc-label drift: TestSpec §4 row **M4** says "2 seqs at 0.9 identity, threshold 0.9 → 1 cluster", but the implemented M4 test uses the 0.875/0.875 boundary (a stronger exact-equality inclusive-cutoff test). Both are valid; no defect. Noted, not blocking.

**Stage A verdict: PASS** — every formula, default, and edge convention matches the authoritative CD-HIT sources retrieved this session.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs`
- `ClusterGenes` lines 210–306 (flatten, stable long→short sort, greedy first-match, cluster assembly, all-pairs mean identity).
- `CalculateSequenceIdentity` lines 314–343 (shorter-length denominator, ungapped prefix scan, empty/empty=1.0, empty/non-empty=0.0).
- `CreatePresenceAbsenceMatrix` lines 352–383.

### Formula realised correctly?
- Shorter-length denominator: `int shorter = Math.Min(len1,len2); … identical/shorter`. Correct (`-G 1`).
- Ungapped prefix scan over `shorter` positions. Correct for substitution/length-only differences (matches ASM-01).
- Greedy first-match: `for c … if (identity >= identityThreshold){ joined=c; break; }`. Correct (fast mode, first representative wins).
- Long→short via `OrderByDescending(length)` — `OrderByDescending` is a stable sort in .NET, so ties keep input order → deterministic, longest is representative. Correct.
- AverageIdentity: singletons 1.0; multi-member = mean over all C(n,2) pairwise identities. Correct.

### Cross-verification table recomputed vs code (full suite run)
All 18 canonical tests pass (17 original + 1 added) with the hand-computed sourced values above; identity 0.875, 1.0, 0.0 and the M6 2-cluster partition reproduce exactly. New M8b 3-member mean = (1.0+0.9+0.9)/3 = 2.8/3 ≈ 0.9333 confirmed against code.

### Variant/delegate consistency
`CreatePresenceAbsenceMatrix` rows derive from cluster `GeneIds`; S5 confirms one row per genome with correct present-gene counts. Consistent.

### Test quality audit (HARD gate)
- **Sourced, not code-echo:** M1–M9, S-cases use exact values traceable to CD-HIT (identity 0.875/1.0/0.0, inclusive 0.875 boundary, shorter-length 1.0 at threshold 1.0, 2-cluster greedy partition). M5 discriminates the shorter-length denominator from the alignment-length convention; M4 pins the inclusive `>=` boundary. These would fail a deliberately-wrong implementation.
- **Two defects fixed this session (test-quality):**
  1. `PanGenomeAnalyzerTests.ClusterGenes_CalculatesAverageIdentity` asserted `AverageIdentity Is.GreaterThan(0.9)` for two **identical** sequences whose exact value is **1.0** — a green-washing lower bound on the unit under test. Tightened to `Is.EqualTo(1.0).Within(1e-10)` with sourced rationale.
  2. The multi-pair AverageIdentity averaging path (cluster size ≥ 3) had no exact-value assertion (M8 only covers a single 2-member pair). Added `ClusterGenes_ThreeMembers_AverageIdentityIsMeanOfAllPairs` pinning the hand-computed 2.8/3.
- **No green-washing introduced:** no widened tolerances, no skipped/ignored tests, no expected value adjusted to match output.
- **Coverage:** identity branches (identical, substitution, length-diff, disjoint, empty/empty, empty/non-empty), inclusive boundary, threshold granularity, greedy multi-cluster, representative = longest, partition invariant, GenomeCount, presence/absence delegate, determinism, [0,1] property, singleton, null/empty/no-genes — all exercised.
- **Honest green:** full unfiltered suite **Passed 6510, Failed 0, Skipped 1** (the pre-existing `MFE_Benchmark_AllScenarios`); changed files build warning-free (the 4 remaining warnings are pre-existing and unrelated).

### Findings / defects
- No implementation defect. Two test-quality defects (one green-washed assertion, one missing formula-path) found and **completely fixed** in this session.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. End-state: CLEAN.**
- Test-quality gate: PASS (after fixing the two test defects above).
- Follow-up (non-blocking): align TestSpec §4 row M4's label (says "0.9 identity") with the implemented 0.875/0.875 boundary test for documentation tidiness.
