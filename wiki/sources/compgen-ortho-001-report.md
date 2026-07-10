---
type: source
title: "Validation report: COMPGEN-ORTHO-001 (ortholog detection by Reciprocal Best Hits + in-paralog identification, ComparativeGenomics.FindOrthologs / FindParalogs / FindReciprocalBestHits)"
tags: [validation, comparative-genomics, governance]
doc_path: docs/Validation/reports/COMPGEN-ORTHO-001.md
sources:
  - docs/Validation/reports/COMPGEN-ORTHO-001.md
source_commit: 0752d91e99ea903a09b9b514deb4d3bbd8688f68
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: COMPGEN-ORTHO-001

The two-stage **validation write-up** for test unit **COMPGEN-ORTHO-001** — **ortholog identification
by Reciprocal Best Hits (RBH)** and **paralog (in-paralog) identification by within-genome best
hits**, validated 2026-06-15. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's independent **verdict** on both the algorithm
description (Stage A) and the shipped code (Stage B). The wider campaign is
[[validation-and-testing]]. The algorithm itself — the RBH rule, the in-paralog proxy, gates,
invariants, worked oracles, and corner cases — is synthesized in the concept
[[ortholog-detection-reciprocal-best-hits]], and [[test-unit-registry]] tracks the unit. Distinct
from [[compgen-ortho-001-evidence]] — the pre-implementation evidence artifact sourced from
`docs/Evidence/` — this page is the independent two-stage re-validation verdict. Its RBH-only sibling
unit is [[compgen-rbh-001-evidence|COMPGEN-RBH-001]].

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS-WITH-NOTES · End state: ✅ CLEAN.** No code defect; the
`PASS-WITH-NOTES` grades are honestly-documented alignment-free simplifications (notes A1–A3, below),
not errors. The full unfiltered suite ran **6506 passed, 0 failed, 0 skipped-of-relevance**
(`MFE_Benchmark` is an unrelated pre-skipped benchmark); `dotnet build` 0 errors. Three
**test-quality improvements** were applied in-session — locking sourced values, not changing
behaviour.

## Canonical methods & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs`:

- `FindOrthologs` (`:334`) → **delegates** to `FindReciprocalBestHits` (`:465`) — single source of
  truth, so the two cannot diverge (locked by new test S5).
- `FindReciprocalBestHits` (`:465`) — null guards (`:471–472`); skips genes without sequence
  (`:474–475`); empty-genome short-circuit (`:477–478`); best hit each direction (`:482–496`); keeps
  only reciprocal pairs by intersecting the two directional best-hit maps (`:500–513`).
- `FindParalogs` (`:356`) — null guard, `<2` short-circuit, within-genome best hit excluding self,
  mutual-best-hit + unordered dedup (the in-paralog proxy).
- private `FindBestHit` (`:410`) — qualification gate `identity ≥ minIdentity ∧ coverage ≥
  minCoverage`; deterministic argmax over (score, coverage, ordinal id).
- private `CalculateSequenceSimilarity` (`:518`) — k=5; identity = Jaccard = shared/|union|;
  coverage = shared/min(|A|,|B|); sequences < 5 nt → (0,0,0).
- Tests: `tests/…/ComparativeGenomics_*` — RBH/ortholog + paralog cases (M1–M6, S1, S5, S6).

## Stage A — description (algorithm faithfulness)

Confirmed against independently-retrieved primary sources: **Fitch 1970** (via Wikipedia *Sequence
homology*, verbatim quote matched) — ortholog = homology by **speciation**, paralog = homology by
**gene duplication**; **Moreno-Hagelsieb & Latimer 2008** (*Bioinformatics* 24(3):319–324, OUP) —
the RBH definition (*"two genes … are deemed orthologs if their protein products find each other as
the best hit in the opposite genome"*), the **≥ 50 % coverage** gate, ranking by bit-score desc then
E-value asc, and the **1×10⁻⁶** E-value gate; **Tatusov 1997** symmetrical best hits (COGs) — a
one-directional best hit is excluded; **Remm et al. 2001** (InParanoid) — clusters seeded by a
two-way best match, in-paralogs added relative to the seed, out-paralogs excluded.

The RBH criterion `BH(g1,G2)=g2 ∧ BH(g2,G1)=g1`, the symmetry / non-reciprocity rejection, the
≥ 50 % coverage gate, and the deterministic best hit with tie-break all match the sources exactly.
An independent **5-mer Jaccard recomputation in Python** for every test dataset reproduced all
asserted identities (1.0 self-match, 0.667 ranking, 0.0 rejection; and notably **TtBlock vs GcBlock =
0.500 Jaccard / 0.667 coverage**, not 0.0).

**Stage A notes (documented simplifications, not defects):**

1. **A1 — alignment-free ranking.** Best-hit ranking uses **5-mer Jaccard similarity**, not a BLAST
   bit-score (recorded as Evidence Assumption 1 / algorithm doc §5.3). The correctness-critical
   reciprocity, coverage gate, and threshold semantics remain source-backed.
2. **A2 — in-paralog proxy.** Remm et al. define in-paralogs *relative to a between-genome seed
   ortholog*; `FindParalogs` uses the simpler **within-genome mutual best hit** and explicitly
   documents that it does not discriminate in- vs out-paralogs (doc §5.3 "Not implemented", §6.2).
   A faithfully-labelled simplification, not a mislabelled exact method.
3. **A3 — coverage gate largely subsumed.** For this set-based metric `coverage = shared/min ≥
   identity = shared/|union|` always; a 200k-random-pair brute force found **no** input passing the
   identity gate (≥ 0.3) while failing the coverage gate (< 0.5). The gate is present and correct per
   source; its independence is a casualty of the alignment-free metric (A1), not a divergence.

## Stage B — implementation

Reciprocity is enforced by intersecting the two directional best-hit maps; the historical
non-reciprocal `FindOrthologs` defect is **already corrected** (verified by test M2 and by hand:
a1↔b1 kept at identity 1.0, b2 excluded). Every worked example reproduced exactly: **M1** (ACGT,Tt
×2) → {a1↔b1, a2↔b2} id 1.0 / cov 1.0; **M2** → {a1↔b1} id 1.0, b2 excluded; **M3** (ACGT vs Tt) →
0 pairs (Jaccard 0); **S1** (3 distinct pairs incl. Tt/Gc) → 3 pairs, perfect matching; **M6**
(p1=p2=Gc, q1=ACGT) → paralogs {p1↔p2} id 1.0, q1 unpaired. `FindOrthologs ≡ FindReciprocalBestHits`
by construction (delegation), locked by S5. No `*Fast`/variant methods.

**Test-quality audit (HARD gate): PASS.** Identity values (1.0 / 0.667 ranking / 0.0 rejection) and
the reciprocity/exclusion outcomes trace to Moreno-Hagelsieb 2008 + the hand-computed Jaccard table,
not code echoes; exact equalities (`EqualTo(1.0).Within(1e-10)`, exact counts and id mappings), no
weakened assertions or skipped tests; all three public methods and all Stage-A branches (mutual hit,
non-reciprocal, sub-threshold, empty, single-gene, no-sequence, null) plus both error cases covered.

**Three test-quality improvements applied in-session** (locking sourced values, not changing
behaviour): (1) fixed a factually-wrong comment claiming TtBlock/GcBlock are "unrelated, Jaccard 0.0"
— they share 8/12 5-mers (Jaccard 0.5); corrected to the hand-computed value and explained why the
RBH matching is unaffected (the 1.0 self-match dominates); (2) added the previously-unasserted
`OrthologPair.Coverage` field assertions (= 1.0 for identical sequences) to M1 and M6; (3) added
direct tests for the public `FindReciprocalBestHits` entry point (S5 equivalence-to-`FindOrthologs`,
S6 null guards) and strengthened M2 to assert the kept pair's identity is the 1.0 reciprocal hit
(guards against picking the 0.667 hit).

## Findings

- **No code defect (State CLEAN).** The historical non-reciprocity defect was already fixed; every
  worked example and hand cross-check reproduced exactly. Notes A1–A3 are honestly-documented
  alignment-free simplifications, not errors.
- **Follow-up (non-blocking):** if an alignment / bit-score becomes available to
  `Seqeron.Genomics.Analysis`, replace the Jaccard ranking metric so the 50 % coverage gate becomes
  an independently-binding constraint (note A3).
