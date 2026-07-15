---
type: source
title: "Validation report: COMPGEN-COMPARE-001 (comprehensive two-genome comparison — core/dispensable partition + syntenic-gene fraction, ComparativeGenomics.CompareGenomes)"
tags: [validation, comparative-genomics, governance]
doc_path: docs/Validation/reports/COMPGEN-COMPARE-001.md
sources:
  - docs/Validation/reports/COMPGEN-COMPARE-001.md
source_commit: 654fe3363f991b6b57549179139dcaa83c31f491
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: COMPGEN-COMPARE-001

The two-stage **validation write-up** for test unit **COMPGEN-COMPARE-001** — the end-to-end
two-genome comparison pipeline `ComparativeGenomics.CompareGenomes` (core/dispensable gene partition
+ overall syntenic-gene fraction), validated 2026-06-16. This is the *report* artifact that feeds one
row of the [[validation-ledger]]; it records the validator's independent **verdict** on both the
algorithm description (Stage A) and the shipped code (Stage B), and the wider campaign is
[[validation-and-testing]]. The algorithm, its invariants, the documented oracles and the edge cases
are synthesized in the concept [[genome-comparison-core-dispensable]];
[[test-unit-registry]] defines the unit. Distinct from
[[compgen-compare-001-evidence]] — the pre-implementation evidence artifact sourced from
`docs/Evidence/` — this page is the independent two-stage re-validation verdict.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: CLEAN.** No defect found; no code / spec / test change
required. The full unfiltered suite ran **6605 passed, 0 failed, 1 skipped** (the skip is an
unrelated MFE benchmark); `dotnet build` 0 errors (4 pre-existing warnings in unrelated files, none
in changed files — no files changed). The 8 unit tests assert exact sourced values, cover all
MUST/SHOULD/COULD cases and all four invariants.

## Canonical method & source under test

`ComparativeGenomics.CompareGenomes(genome1Genes, genome2Genes, minOrthologIdentity = 0.3,
minSyntenicBlockSize = 3)` — `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:765`
(code path reviewed `:765–810`). `CompareGenomes` is an **aggregator**: it delegates to already-validated
sub-units `FindReciprocalBestHits` (COMPGEN-RBH-001), `FindSyntenicBlocks` (COMPGEN-SYNTENY-001) and
`DetectRearrangements` (COMPGEN-REARR-001), and adds exactly three pieces of unit-specific logic,
each independently sourced this session.

## Stage A — description (algorithm faithfulness)

The three non-delegated pieces of logic and their sources:

1. **Core / dispensable partition** — Tettelin et al. 2005, *PNAS* 102(39):13950–13955 (WebFetch of
   PMC1216834). Verbatim: the pan-genome *"includes a core genome containing genes present in all
   strains and a dispensable genome composed of genes absent from one or more strains and genes that
   are unique to each strain."* Pairwise specialisation: present in both → **core**
   (`ConservedGenes = |orthologs|`); present in one only → **dispensable / genome-specific**
   (`GenomeSpecificGenes_i = |genome_i| − |ortholog genes of genome_i|`).
2. **Shared gene = reciprocal best hit** — Moreno-Hagelsieb & Latimer 2008, *Bioinformatics*
   24(3):319–324 (WebFetch of the OUP article). Verbatim: *"two genes residing in two different
   genomes are deemed orthologs if their protein products find each other as the best hit in the
   opposite genome,"* coverage gate *"at least 50% of any of the protein sequences."* The conserved
   set = RBH pairs (RBH itself validated in COMPGEN-RBH-001).
3. **OverallSynteny = fraction of syntenic genes** — confirmed via WebSearch as a real
   comparative-genomics metric (coral *Acropora* 50.2%; Rhizobiales ~80% of orthologs syntenic;
   Fugu/human 40–50%). Code defines `OverallSynteny = Σ block.GeneCount / min(|g1|,|g2|)`, clamped
   ≤ 1, with blocks from MCScanX (score ≥ 250 ≈ ≥ 5 anchors, Wang et al. 2012, COMPGEN-SYNTENY-001).

**Four invariants confirmed genuine.** INV-01 `Conserved = |orthologs|`; INV-02 `core + specific_i =
|genome_i|` (holds because RBH is a matching — each gene in ≤ 1 pair, the distinct `Gene1Id`/`Gene2Id`
guarantee from COMPGEN-RBH-001); INV-03 `OverallSynteny ∈ [0,1]` (lower bound natural, upper bound by
`Math.Min(1.0, …)`); INV-04 swap symmetry (from reciprocal RBH).

**Edge-case semantics** (all sourced/derived): all-shared → both specific counts 0 (Tettelin);
disjoint → conserved 0, all dispensable; empty genomes → all zeros; < 5 collinear orthologs →
conserved > 0 but synteny 0 (MCScanX threshold).

**Stage A finding (BY-DESIGN, not a defect).** One documented, sourced simplification inherited from
COMPGEN-RBH-001: the Tettelin "50% conservation over 50% length" alignment gate is approximated by
alignment-free 5-mer Jaccard (id ≥ 0.3, cov ≥ 0.5). It does not change the partition logic
(identical sequences pass, disjoint fail). → **Stage A PASS.**

## Stage B — implementation

Cross-verification table recomputed vs code, full suite executed this session — every case matched:

| Case | Sourced expectation | Match |
|------|--------------------|-------|
| M1 one-shared/one-unique | Conserved 1, Spec 1/1, Orthologs 1 | ✓ |
| M2 disjoint | Conserved 0, Spec 2/2, Orthologs ∅ | ✓ |
| C1 all shared | Conserved 2, Spec 0/0 | ✓ |
| M3 5 collinear + 1 unique | Conserved 5, Spec 1/1, 1 block of 5, Synteny 5/6, 0 rearr | ✓ |
| S1 3 collinear (< 5) | Conserved 3, no block, Synteny 0 | ✓ |
| S2 swap symmetry | Conserved 2 fixed, Spec 2/1 ↔ 1/2 | ✓ |
| M4 empty | all 0, all collections ∅ | ✓ |
| Null | ArgumentNullException ×2 | ✓ |

The M3 synteny value `5/6 = 0.8333…` was hand-traced independently: one block of 5 collinear anchors
(MCScanX score 5×50 = 250 ≥ 250 → reported, GeneCount 5); `min(6,6) = 6`.

**Test-quality audit (HARD gate) — PASS.** Every expected value is an exact constant tied to a
retrieved source, not a code echo; M1's `Conserved=1 / Spec 1/1` would *fail* a one-directional /
no-ortholog wrong implementation, and M3's `5/6` is recomputed independently. No green-washing: all
assertions exact (`Is.EqualTo`, `Within(1e-10)` only for the rational 5/6 and 0 floats), no
`Greater/AtLeast/Contains`, no ranges, no skips, no widened tolerance (the prior permissive
`GreaterThan`/`GreaterThanOrEqualTo` tests were removed). Coverage: all 4 MUST (M1–M4), both SHOULD
(S1 block threshold, S2 symmetry), the COULD all-shared (C1), the null contract, plus all four
invariants — 8 tests, all on the public `CompareGenomes`.

**Stage B notes (BY-DESIGN, no change made):** (a) the `minOrthologIdentity` / `minSyntenicBlockSize`
non-default overload arguments are not exercised — they merely forward to the sub-methods validated in
COMPGEN-RBH-001/SYNTENY-001; (b) the `Math.Min(1.0,…)` upper clamp is a defensive guard not separately
triggered (MCScanX blocks are non-overlapping, so `Σ GeneCount ≤ min(|g1|,|g2|)` structurally). Both are
"nice-to-have", not MUST. → **Stage B PASS.**

## Findings

- **No code defect and no test change (End state CLEAN).** Description is sourced and correct; the
  implementation faithfully realises it; tests assert exact sourced values with full MUST/SHOULD/COULD
  and all-four-invariant coverage, and pass in the full unfiltered suite (6605, Failed 0).
- **Test-quality gate: PASS** (sourced expectations, no green-washing, honest green).
- **No follow-ups.** The two Stage-B coverage notes (non-default overload args; the upper-clamp guard)
  are documented nice-to-haves, not defects.
</content>
</invoke>
