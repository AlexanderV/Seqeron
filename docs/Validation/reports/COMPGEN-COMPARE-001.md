# Validation Report: COMPGEN-COMPARE-001 — Comprehensive Two-Genome Comparison (core/dispensable partition + syntenic-gene fraction)

- **Validated:** 2026-06-16   **Area:** Comparative Genomics
- **Canonical method(s):** `ComparativeGenomics.CompareGenomes(genome1Genes, genome2Genes, minOrthologIdentity = 0.3, minSyntenicBlockSize = 3)` (`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:765`); delegates to validated sub-units `FindReciprocalBestHits` (COMPGEN-RBH-001), `FindSyntenicBlocks` (COMPGEN-SYNTENY-001), `DetectRearrangements` (COMPGEN-REARR-001).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

`CompareGenomes` is an **aggregator**. The unit-specific logic it defines (not delegated) is exactly three things, each independently sourced this session:

1. **Core / dispensable partition (Tettelin et al. 2005, PNAS 102(39):13950–13955).**
   - Retrieved this session via WebFetch of PMC1216834. Verbatim quote confirmed: the pan-genome *"includes a core genome containing genes present in all strains and a dispensable genome composed of genes absent from one or more strains and genes that are unique to each strain."*
   - Also confirmed verbatim: *"A gene was considered conserved if at least one of these three methods produced an alignment with a minimum of 50% sequence conservation over 50% of the protein/gene length."*
   - Mapping to two genomes (the pairwise case): "present in all" = present in both → **core**; "unique to each strain" / "absent from one" = present in only one → **dispensable / genome-specific**. This mapping is the standard pairwise specialisation and is the partition the code implements (`ConservedGenes = |orthologs|`; `GenomeSpecificGenes_i = |genome_i| − (genes of genome_i in an ortholog pair)|`).

2. **Shared gene = reciprocal best hit (Moreno-Hagelsieb & Latimer 2008, Bioinformatics 24(3):319–324).**
   - Retrieved this session via WebFetch of the OUP article. Verbatim: *"two genes residing in two different genomes are deemed orthologs if their protein products find each other as the best hit in the opposite genome,"* with coverage gate *"coverage of at least 50% of any of the protein sequences in the alignments."*
   - This is the operational orthology criterion; the conserved set = RBH pairs. (RBH itself is validated in COMPGEN-RBH-001.)

3. **OverallSynteny = "fraction of syntenic genes" (synteny conservation metric).**
   - Retrieved this session via WebSearch. Confirmed as a real comparative-genomics metric: the proportion of genes that are syntenic (e.g. coral *Acropora* 50.2%; Rhizobiales ~80% of orthologs syntenic; Fugu/human 40–50%). The code defines it as `Σ block.GeneCount / min(|g1|,|g2|)`, clamped to ≤1, with blocks coming from MCScanX (score ≥250 ≈ ≥5 anchors, Wang et al. 2012) — validated in COMPGEN-SYNTENY-001.

**Formula check.** `ConservedGenes = orthologs.Count` (INV-01 ✓); `GenomeSpecificGenes_i = |genome_i| − |orthologGenes_i|` with `core + specific_i = |genome_i|` (INV-02) — holds because RBH is a matching so each gene appears in ≤1 pair (the distinct `Gene1Id`/`Gene2Id` guarantee from COMPGEN-RBH-001); `OverallSynteny ∈ [0,1]` (INV-03, lower bound natural, upper bound by `Math.Min(1.0, …)`); swap symmetry (INV-04) from reciprocal RBH. All four invariants are genuine mathematical properties given the sourced definitions.

**Edge-case semantics.** All-shared → both specific counts 0 (Tettelin); disjoint → conserved 0, all dispensable; empty genomes → all zeros (no pairs possible); <5 collinear orthologs → conserved>0 but synteny 0 (MCScanX threshold). All sourced/derived.

**Independent cross-check (numbers).**
- Partition (M1/M2/C1) is pure set arithmetic over the sourced definitions: 1 shared + 1 unique each → core 1, specific 1/1; disjoint → 0, 2/2; all shared → 2, 0/0. Independent of code.
- Synteny (M3): one block of 5 collinear anchors (MCScanX score 5×50 = 250 ≥ 250 → reported, GeneCount 5); `min(6,6) = 6`; `OverallSynteny = 5/6 = 0.8333…`. Hand-traced, independent of the implementation's return value.

**Findings / divergences (Stage A).** One documented, sourced simplification (BY-DESIGN, inherited from COMPGEN-RBH-001): the Tettelin "50% conservation over 50% length" alignment gate is approximated by alignment-free 5-mer Jaccard (id ≥0.3, cov ≥0.5). This does not change the partition logic — identical sequences pass, disjoint fail — which is all the partition tests rely on. Not a defect. → Stage A **PASS**.

## Stage B — Implementation

**Code path reviewed.** `ComparativeGenomics.cs:765–810`.
- `orthologs = FindReciprocalBestHits(...)` → `ConservedGenes = orthologs.Count` (INV-01 by construction). ✓
- `specific_i = genome_i.Count(g => !orthologGenes_i.Contains(g.Id))` over `HashSet` of pair gene ids → INV-02 holds because RBH yields distinct gene ids per side. ✓
- `OverallSynteny = Σ block.GeneCount / min(|g1|,|g2|)` when blocks exist and smaller>0, else 0, then `Math.Min(1.0, …)` (INV-03). ✓
- `ArgumentNullException.ThrowIfNull` on both lists; empty lists → all-zero result. ✓

**Cross-verification table recomputed vs code (full suite executed this session).**

| Case | Sourced expectation | Code output | Match |
|------|--------------------|-------------|-------|
| M1 one-shared/one-unique | Conserved 1, Spec 1/1, Orthologs 1 | same | ✓ |
| M2 disjoint | Conserved 0, Spec 2/2, Orthologs ∅ | same | ✓ |
| C1 all shared | Conserved 2, Spec 0/0 | same | ✓ |
| M3 5 collinear + 1 unique | Conserved 5, Spec 1/1, 1 block of 5, Synteny 5/6, 0 rearr | same | ✓ |
| S1 3 collinear (<5) | Conserved 3, no block, Synteny 0 | same | ✓ |
| S2 swap symmetry | Conserved 2 fixed, Spec 2/1 ↔ 1/2 | same | ✓ |
| M4 empty | all 0, all collections ∅ | same | ✓ |
| Null | ArgumentNullException ×2 | same | ✓ |

**Variant/delegate consistency.** Orthologs/blocks/rearrangements are produced by the already-validated sub-methods; the aggregator only adds the partition + synteny-fraction arithmetic, which is verified above.

**Test quality audit (HARD gate).**
- **Sourced, not code echoes:** every expected value is an exact constant tied to a retrieved source (Tettelin core/dispensable, RBH, fraction-of-syntenic-genes). M1's `Conserved=1 / Spec 1/1` would *fail* against a one-directional or no-ortholog (wrong) implementation, so it rejects the documented failure mode. M3's `5/6` is recomputed from the formula independently of code.
- **No green-washing:** all assertions are exact (`Is.EqualTo`, `Within(1e-10)` only for the rational 5/6 and the 0 floats), no `Greater/AtLeast/Contains`, no ranges, no skip/ignore, no widened tolerance. (The TestSpec confirms the prior permissive `GreaterThan`/`GreaterThanOrEqualTo` tests were removed.)
- **Coverage:** all 4 Stage-A MUST cases (M1–M4), both SHOULD (S1 block threshold, S2 symmetry), the COULD all-shared (C1), the null-contract path, plus all four invariants INV-01…INV-04 are exercised. 8 tests, all on the public `CompareGenomes`.
- **Honest green:** FULL unfiltered suite ran this session = **6605 passed, 0 failed, 1 skipped** (the skip is an unrelated MFE benchmark); `dotnet build` 0 errors (4 pre-existing warnings in unrelated files, none in changed files — no files changed).

**Findings / defects (Stage B).** None. Minor, non-defect coverage notes (BY-DESIGN, no change made): (a) the `minOrthologIdentity`/`minSyntenicBlockSize` non-default overload arguments are not exercised here — they merely forward to sub-methods validated in COMPGEN-RBH-001/SYNTENY-001, and the aggregator-specific logic is fully covered; (b) the `Math.Min(1.0,…)` upper clamp (INV-03 upper bound) is a defensive guard not separately triggered (MCScanX blocks are non-overlapping so `Σ GeneCount ≤ min(|g1|,|g2|)` structurally); both are "nice-to-have" per the Evidence recommendations, not MUST. → Stage B **PASS**.

## Verdict & follow-ups

- **Stage A: PASS. Stage B: PASS. End-state: CLEAN.** No defect found; no code/spec/test change required. The description is sourced and correct; the implementation faithfully realises it; the tests assert exact sourced values, cover all MUST/SHOULD/COULD cases and all four invariants, and pass in the full unfiltered suite (6605, Failed: 0).
- **Test-quality gate: PASS** (sourced expectations, no green-washing, full Stage-A coverage, honest green).
