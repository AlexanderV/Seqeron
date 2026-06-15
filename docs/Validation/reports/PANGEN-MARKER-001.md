# Validation Report: PANGEN-MARKER-001 — Phylogenetic Marker Selection

- **Validated:** 2026-06-15   **Area:** PanGenome
- **Canonical method(s):** `PanGenomeAnalyzer.CountParsimonyInformativeSites(IReadOnlyList<string>)`,
  `PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, coreClusters, totalGenomes, maxMarkers)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session
1. **Wikipedia "Informative site"** (cites Zvelebil & Baum 2008, *Understanding Bioinformatics*) —
   https://en.wikipedia.org/wiki/Informative_site — retrieved verbatim:
   "The informative site is a position in the relevant set of aligned sequences at which there are
   **at least two different character states and each of those states occurs in at least two of the
   sequences**." Uninformative = invariable (monomorphic) sites and singleton sites (a difference in
   only one sequence). Worked diagram: sites 2 and 4 informative; site 1 invariable, site 3 singleton.
2. **Memorial University parsimony notes** — https://www.mun.ca/biology/scarr/2900_Parsimony_Analysis.htm —
   confirms: "a parsimony-informative site must have at least two character states, each appearing in at
   least two taxa." Worked positions: informative position has state a in two taxa + state g in two taxa;
   singleton (one taxon differs) and all-distinct (4 states) positions are uninformative.
3. **IQ-TREE FAQ** (via WebSearch) — "Parsimony informative sites are patterns that have at least two
   different characters … and each character should occur in at least two species."
4. **panX (Ding, Baumdicker, Neher 2018)** and **Roary (Page et al. 2015)** + Roary docs — already cited
   in the Evidence doc for the single-copy-core / 99%-core / paralog-filtering rules (not re-fetched this
   session; the Evidence doc records the verbatim quotes and they match the standard panX/Roary practice).

### Formula check
- **Parsimony-informative column criterion** (Zvelebil 2008 / IQ-TREE / Memorial U.): a column is PI iff
  there exist **≥ 2 distinct states, each occurring in ≥ 2 sequences**. Monomorphic and singleton columns
  are excluded. All three independent sources agree exactly. The repo description (`§2.2`, Evidence §1.2.5)
  copies this verbatim. ✔
- **Disambiguation of the "each of those states" wording.** The literal phrase could be misread as "every
  state present must occur ≥ 2 times," which would wrongly mark a column like A,A,C,C,G (A:2,C:2,G:1) as
  uninformative. All three sources resolve this: the operative requirement is "**at least two** states each
  in ≥ 2 taxa" (the mutation-count argument: a singleton state cannot change the relative parsimony cost of
  competing topologies, so its presence is irrelevant once two supported states exist). The repo
  implementation uses exactly this reading (`statesWithSupport >= 2`), which is correct. ✔
- **Marker selection rule** (panX/Roary): keep single-copy core clusters — present in all genomes with
  exactly one gene each (panX "all strains represented exactly once"; Roary paralog filtering) — extract
  variable positions, keep clusters with ≥ 1 PI site, rank by descending PIS. The description matches. ✔

### Edge-case semantics check
All documented and sourced: monomorphic → 0; singleton → 0; 4-singletons → 0; two-states-each-≥2 → 1;
< 2 rows → 0; unequal-length members → 0 (Assumption 1, no in-repo aligner); paralog / non-core / fully
conserved clusters excluded. INV-1 bound 0 ≤ PIS ≤ length. All consistent with sources.

### Independent cross-check (numbers)
Hand-computed against the Zvelebil definition:

| Column | States | PI? | Source agreement |
|--------|--------|-----|------------------|
| A,A,A,A | A:4 | No (mono) | Wikipedia site 1 |
| A,A,A,C | A:3,C:1 | No (singleton) | Wikipedia site 3 |
| A,A,C,C | A:2,C:2 | **Yes** | Wikipedia site 2 / Memorial U. |
| A,C,G,T | all 1 | No | Memorial U. (all-distinct) |
| A,A,C,C,G | A:2,C:2,G:1 | **Yes** | resolved above (singleton irrelevant) |
| A,A,A,C,C | A:3,C:2 | **Yes** | two supported states |

M1 worked alignment (AAAAA/AAACA/AACCG/ACCTG) → cols 3 & 5 PI → **PIS = 2**. Confirmed by hand.

### Findings / divergences
None. The description, formula, edge cases and the subtle "two-supported-states" disambiguation are all
correct and source-backed. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs:934-976`
  (`CountParsimonyInformativeSites`) and `:997-1043` (`SelectPhylogeneticMarkers`).

### Formula realised correctly? (evidence)
- PIS: tallies per-column state counts, skips columns with < 2 distinct states (monomorphic), then counts
  states with count ≥ 2 (`MinSequencesPerState = 2`) and marks the column informative iff ≥ 2 such states
  (`statesWithSupport >= MinDistinctStates`). This is exactly the sourced "≥ 2 states each in ≥ 2 taxa"
  criterion, including the correct handling of a third singleton state. ✔
- Guards: null / < 2 rows / first-row length ≤ 0 / any row null-or-unequal-length → 0. Matches Assumption 1. ✔
- Selection: builds geneId→sequence map; per cluster requires `GenomeCount == totalGenomes` AND
  `GeneIds.Count == totalGenomes` (single-copy core: no missing, no paralog); recovers members; counts PIS;
  keeps PIS ≥ 1; orders by descending PIS, ties by ordinal ClusterId; `Take(maxMarkers)`. Matches panX/Roary. ✔
- Degenerate-input guard: null genomes/clusters, totalGenomes ≤ 0, maxMarkers ≤ 0 → empty. ✔

### Cross-verification table recomputed vs code
Ran the full test suite; values match hand computation:
- M1 alignment → 2; mono → 0; singleton → 0; A,A,C,C → 1; A,C,G,T → 0;
- **NEW** A,A,C,C,G → 1; A,A,A,C,C → 1; all-informative 3 cols → 3 (upper bound);
- selection: only single-copy-core informative cluster selected; paralog/non-core/conserved excluded;
  descending-PIS ranking (2 before 1) overrides ordinal id; maxMarkers cap = 1; equal-PIS tie by id;
  unequal-length members → not selected; null/empty → empty.

### Variant/delegate consistency
Two public static methods; `SelectPhylogeneticMarkers` delegates PIS to `CountParsimonyInformativeSites`.
Consistent; no `*Fast` variant.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** all assertions use exact `Is.EqualTo` to values derived from the Zvelebil /
  Memorial U. / IQ-TREE definitions, not from running the code. ✔
- **No green-washing:** exact equality everywhere; no weakened/widened/skipped assertions; ordering tests
  use full `Is.EqualTo(new[]{…})` sequence equality. ✔
- **Coverage gap found and FIXED:** the original 17 tests did **not** cover the disambiguating branch —
  a column with two supported states *plus* a third singleton state (A,A,C,C,G), nor unequal support counts
  (A,A,A,C,C), nor the INV-1 upper bound (PIS = length). A deliberately-wrong implementation reading
  "every state must occur ≥ 2 times" would pass all 17 original tests yet fail biology. Added three tests
  to lock the sourced interpretation:
  - `CountParsimonyInformativeSites_TwoSupportedStatesPlusSingleton_ReturnsOne` (A,A,C,C,G → 1)
  - `CountParsimonyInformativeSites_UnequalSupportCounts_ReturnsOne` (A,A,A,C,C → 1)
  - `CountParsimonyInformativeSites_AllColumnsInformative_EqualsAlignmentLength` (→ 3, INV-1 upper bound)
- **All Stage-A branches exercised:** mono / singleton / minimal-PI / 4-singleton / mixed-with-singleton /
  unequal-support / all-informative / < 2 rows / null / empty / unequal-length; selection: single-copy-core
  select, paralog exclude, non-core exclude, conserved exclude, descending-PIS rank, maxMarkers cap,
  equal-PIS tie-break, null/empty/totalGenomes≤0/maxMarkers≤0. ✔
- **Honest green:** FULL unfiltered suite = **Failed: 0, Passed: 6548** (was 6545; +3 new). Build 0 errors,
  no new warnings (4 pre-existing warnings in unrelated files).

### Findings / defects
No implementation defect. One **test-coverage gap** (missing disambiguating branch) — fixed in-session by
adding three sourced exact-value tests. **Stage B: PASS.**

## Verdict & follow-ups
- **Stage A: PASS** — description, formula, and the subtle "≥ 2 supported states" disambiguation are
  correct and confirmed against three independent authoritative sources.
- **Stage B: PASS** — code faithfully realises the validated criterion; tests strengthened to lock the
  disambiguating branch and the INV-1 upper bound.
- **End-state: CLEAN.** Test-quality gate: PASS (gap found and fully fixed; full suite green).
- No open defects.
