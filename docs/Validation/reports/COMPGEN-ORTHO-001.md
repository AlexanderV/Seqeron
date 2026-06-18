# Validation Report: COMPGEN-ORTHO-001 — Ortholog Identification (Reciprocal Best Hits) and Paralog Identification

- **Validated:** 2026-06-15   **Area:** Comparative Genomics
- **Canonical method(s):** `ComparativeGenomics.FindOrthologs`, `ComparativeGenomics.FindParalogs`, `ComparativeGenomics.FindReciprocalBestHits` (+ private `FindBestHit`, `CalculateSequenceSimilarity`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** ✅ CLEAN (no code defect; tests strengthened in-session)

## Stage A — Description

### Sources opened this session (independent of the repo)

| Source | Retrieved | What it confirms |
|--------|-----------|------------------|
| Wikipedia *Sequence homology*, quoting Fitch 1970 | WebFetch | Ortholog = homology by **speciation**; paralog = homology by **gene duplication**. Verbatim Fitch quote matches the repo Evidence quote. |
| Moreno-Hagelsieb & Latimer 2008, *Bioinformatics* 24(3):319–324 (OUP article) | WebFetch | RBH def: *"two genes residing in two different genomes are deemed orthologs if their protein products find each other as the best hit in the opposite genome."* Coverage gate: *"coverage of at least 50% of any of the protein sequences."* Ranking: sort by bit-score desc, then E-value asc (deterministic tie-break). E-value gate 1×10⁻⁶. |
| WebSearch — RBH/BBH definition (multiple hits incl. OUP, bioRxiv, widdowquinn RBH tutorial) | WebSearch | RBH/BBH = mutual best hits → inferred orthologs; orthology is symmetric so orthologs are reciprocal best matches. A one-directional best hit is NOT sufficient. |
| WebSearch — InParanoid / Remm et al. 2001 | WebSearch | Ortholog clusters **seeded by a two-way (reciprocal) best match**; in-paralogs added relative to that seed; outparalogs excluded. |

### Formula / definition check

- **Orthology vs paralogy (Fitch 1970):** matches the source verbatim. ✅
- **RBH ortholog criterion:** `BH(g1,G2)=g2 ∧ BH(g2,G1)=g1` matches Moreno-Hagelsieb & Latimer (2008) exactly. ✅
- **Symmetry / non-reciprocity rejection (Tatusov 1997 symmetrical best hits):** correct — a one-directional best hit is excluded. ✅
- **Coverage ≥ 50% gate:** present and source-backed. ✅
- **Deterministic best hit with tie-break:** source ranks by bit-score then E-value; repo ranks by score then coverage then ordinal id — both are deterministic total orders, consistent with the source's intent (a unique best hit). ✅

### Edge-case semantics

Empty/single-gene genome → no pair; gene without sequence skipped; sub-threshold excluded; one-directional excluded — all defined and sourced (definitions require a pair; coverage/identity gates per Moreno-Hagelsieb 2008).

### Independent cross-check (hand-computed, this session)

5-mer Jaccard recomputed in Python for every test dataset:

| Pair | shared | union | identity (Jaccard) | coverage (shared/min) |
|------|--------|-------|--------------------|------------------------|
| ACGTrepeat vs itself | 4 | 4 | **1.000** | 1.000 |
| ACGTrepeat vs AcgtPlusAA | 4 | 6 | **0.667** | 1.000 |
| ACGTrepeat vs TtBlock | 0 | 16 | 0.000 | 0.000 |
| GcBlock vs itself | 12 | 12 | 1.000 | 1.000 |
| GcBlock vs ACGTrepeat | 0 | 16 | 0.000 | 0.000 |
| **TtBlock vs GcBlock** | **8** | **16** | **0.500** | **0.667** |

All identity values asserted by the tests (1.0, 0.667 ranking, 0.0 rejection) reproduce exactly.

### Findings / divergences (Stage A)

- **NOTE A1 (documented simplification, not a defect):** best-hit ranking uses alignment-free 5-mer Jaccard similarity, not a BLAST bit-score. Honestly recorded as Evidence Assumption 1 and in the algorithm doc §5.3. The correctness-critical reciprocity, coverage gate, and threshold semantics are source-backed.
- **NOTE A2 (in-paralog proxy):** Remm et al. (2001) define in-paralogs *relative to a between-genome seed ortholog* (within-species pair more similar to each other than to the seed). The repo's `FindParalogs` uses the simpler "within-genome mutual best hit" proxy and explicitly documents that it does not discriminate in- vs out-paralogs (doc §5.3 "Not implemented", §6.2). This is a faithfully-labelled simplification, not a mislabelled exact method.
- **NOTE A3 (metric property):** for this set-based metric, `coverage = shared/min(|A|,|B|) ≥ identity = shared/|union|` always. A brute-force search (200k random DNA pairs) found **no** input that passes the identity gate (≥0.3) while failing the coverage gate (<0.5). The 50% coverage gate is therefore largely subsumed by the identity gate for realistic short DNA inputs. The gate is present and correct per source; this is a consequence of the alignment-free metric (Assumption 1), not a divergence from the description.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs`:
- `FindOrthologs` (l.334) → delegates to `FindReciprocalBestHits` (l.465) — single source of truth, cannot diverge. ✅
- `FindReciprocalBestHits`: null guards (l.471–472), skips genes without sequence (l.474–475), empty-genome short-circuit (l.477–478), best-hit each direction (l.482–496), keeps only reciprocal pairs (l.500–513). ✅
- `FindParalogs` (l.356): null guard, `<2` short-circuit, within-genome best hit excluding self, mutual-best-hit + unordered dedup. ✅
- `FindBestHit` (l.410): qualification gate `identity≥minIdentity ∧ coverage≥minCoverage`; deterministic argmax (score, coverage, ordinal id). ✅
- `CalculateSequenceSimilarity` (l.518): k=5; identity = Jaccard; coverage = shared/min(|A|,|B|); seqs <5 nt → (0,0,0). ✅

### Formula realised correctly?

Yes. Reciprocity is enforced by intersecting the two directional best-hit maps; the non-reciprocal `FindOrthologs` defect noted in the spec history has been corrected (verified by M2 test and by hand: a1↔b1 kept at identity 1.0, b2 excluded). Hand-recomputed RBH matchings for M1 and S1 match the code output exactly.

### Cross-verification table recomputed vs code

| Test | Hand result | Code result | Match |
|------|-------------|-------------|-------|
| M1 (ACGT,Tt ×2) | {a1↔b1, a2↔b2}, id 1.0, cov 1.0 | same | ✅ |
| M2 (a1; b1=ACGT, b2=AcgtPlusAA) | {a1↔b1} id 1.0; b2 excluded | same | ✅ |
| M3 (ACGT vs Tt) | 0 pairs (Jaccard 0) | 0 | ✅ |
| S1 (3 distinct pairs incl. Tt/Gc) | 3 pairs, perfect matching | 3 | ✅ |
| M6 (p1=p2=Gc, q1=ACGT) | {p1↔p2} id 1.0; q1 unpaired | same | ✅ |

### Variant/delegate consistency

`FindOrthologs` ≡ `FindReciprocalBestHits` by construction (delegation); locked by new test S5. ✅

### Test quality audit (HARD gate)

- **Sourced expectations:** identity values (1.0 / 0.667 ranking / 0.0 rejection) and the reciprocity/exclusion outcomes trace to the externally-retrieved Moreno-Hagelsieb (2008) RBH definition and the hand-computed Jaccard table — not code echoes. ✅
- **No green-washing:** exact equalities used (`EqualTo(1.0).Within(1e-10)`, exact counts, exact id mappings); no weakened assertions, widened tolerances, or skipped tests. ✅
- **Coverage of logic:** all three public methods now exercised, all Stage-A branches (mutual hit, non-reciprocal, sub-threshold, empty, single-gene, no-sequence, null) and both error cases covered.
- **Honest green:** full unfiltered suite **6506 passed, 0 failed, 0 skipped-of-relevance** (`MFE_Benchmark` is an unrelated pre-skipped benchmark); `dotnet build` 0 errors. Pre-existing NUnit2007 warnings are in `ApproximateMatcher_EditDistance_Tests.cs`, untouched by this unit.

### Findings / defects (Stage B)

No code defect found. Three **test-quality improvements** applied in-session (locking sourced values, not changing behaviour):

1. **Fixed a factually-wrong comment** claiming `TtBlock`/`GcBlock` are "unrelated, Jaccard 0.0" — they actually share 8/12 5-mers (Jaccard 0.5). Corrected to the hand-computed value; explained why the RBH matching is still unaffected (the 1.0 self-match dominates).
2. **Added missing coverage-field assertions** to M1 and M6 (coverage = 1.0 for identical sequences) — the documented `OrthologPair.Coverage` output was previously unasserted anywhere.
3. **Added direct tests** for the public `FindReciprocalBestHits` entry point (S5 equivalence-to-`FindOrthologs`; S6 null guards) — previously only exercised indirectly via delegation. Strengthened M2 to assert the kept pair's identity is the 1.0 reciprocal hit (guards against picking the 0.667 hit).

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** — all definitions/formulas confirmed against primary sources; notes A1–A3 are honestly-documented alignment-free simplifications, not errors.
- **Stage B: PASS-WITH-NOTES** — code faithfully realises the validated RBH/in-paralog description; reciprocity defect already corrected; tests strengthened to lock sourced values and cover the third public method.
- **End-state: ✅ CLEAN.** No open defect. Follow-up (non-blocking): if an alignment/bit-score becomes available to `Seqeron.Genomics.Analysis`, replace the Jaccard ranking metric so the 50% coverage gate becomes an independently-binding constraint (note A3).
