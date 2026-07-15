---
type: source
title: "Validation report: COMPGEN-RBH-001 (reciprocal best hits / RBH ortholog detection, ComparativeGenomics.FindReciprocalBestHits)"
tags: [validation, comparative-genomics, governance]
doc_path: docs/Validation/reports/COMPGEN-RBH-001.md
sources:
  - docs/Validation/reports/COMPGEN-RBH-001.md
source_commit: 00c5ea423943b75acdee59bd9ce88cc801bb2a37
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: COMPGEN-RBH-001

The two-stage **validation write-up** for test unit **COMPGEN-RBH-001** — **reciprocal best hits
(RBH / BBH)**, the core between-genome **ortholog-detection primitive**
`ComparativeGenomics.FindReciprocalBestHits`, validated 2026-06-16. This is the *report* artifact that
feeds one row of the [[validation-ledger]]; it records the validator's independent **verdict** on both
the algorithm description (Stage A) and the shipped code (Stage B). The wider campaign is
[[validation-and-testing]] and [[test-unit-registry]] tracks the unit. The algorithm itself — the RBH
reciprocity rule, the two gates, invariants, worked oracles, and corner cases — is synthesized in the
concept [[ortholog-detection-reciprocal-best-hits]] (RBH **is** its core method). This page is the
independent re-validation verdict; keep it distinct from **(a)** [[compgen-rbh-001-evidence]] — the
pre-implementation Evidence artifact sourced from `docs/Evidence/` — and **(b)**
[[compgen-ortho-001-report]], the sibling report for the broader **RBH + in-paralog** unit
COMPGEN-ORTHO-001. COMPGEN-RBH-001 is the **RBH-only slice**: the between-genome ortholog rule without
the within-genome in-paralog rule.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS · End state: ✅ CLEAN.** No code defect and no code change.
Stage A is with-notes only for the honestly-documented alignment-free simplification (5-mer Jaccard
ranking substitutes for a BLAST bit-score, note below); the correctness-critical parts are all
source-backed. The only actionable work this session was **closing two test-coverage gaps** (a missing
`minCoverage`-gate test and a missing `< k` short-sequence-edge test) — a test-surface gap, not an
implementation defect. Full unfiltered suite ran **6605 passed, 0 failed, 1 skipped** (the pre-existing
unrelated `MFE_Benchmark`); `dotnet build` 0 warnings, 0 errors.

## Canonical methods & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs`:

- `FindReciprocalBestHits` (`:465–516`) — null-checks both lists; drops genes without a sequence;
  computes best hit each direction; emits a pair **iff** g1's best hit is g2 **AND** g2's best hit
  back is g1 (the reciprocity conjunction).
- private `FindBestHit` (`:410–440`) — argmax identity; tie-break by larger coverage then ordinal id →
  unique deterministic winner; qualifies only when `identity ≥ minIdentity ∧ coverage ≥ minCoverage`
  (both gates).
- private `CalculateSequenceSimilarity` (`:518–549`) — k=5; identity = Jaccard = shared/|union|;
  coverage = shared/min(|A|,|B| k-mer counts); alignLen = min sequence length; sequences `< k` →
  (0,0,0).
- `FindOrthologs` (`:334–341`) — **pure delegation** to `FindReciprocalBestHits`; no divergence
  possible (smoke test asserts exact tuple equality).

## Stage A — description (algorithm faithfulness)

Confirmed against independently-retrieved primary sources this session: **Moreno-Hagelsieb & Latimer
2008** (*Bioinformatics* 24(3):319–324) — the operational RBH definition confirmed verbatim from the
PubMed abstract (*"orthologs are assumed if two genes each in a different genome find each other as the
best hit in the other genome"*); the `≥50%` coverage / `1×10⁻⁶` E-value body quotes are paywalled but
the load-bearing *definition* is independently confirmed. **WebSearch survey** of the RBH/BBH
literature (Genome Biol. Evol. "Bidirectional Best Hits Miss Many Orthologs"; the Best-Match-Graph
arXiv series 1803.10989 / 1903.07920 / 2001.00958) confirms the **symmetric requirement** — a
one-directional best hit is insufficient. **Tatusov, Koonin & Lipman 1997** via NCBI Handbook NBK21090
confirms COGs are built from *"triangles of mutually consistent, genome-specific best hits (BeTs)"* —
the pairwise special case of a mutually-consistent BeT is exactly the reciprocal best hit.

The formula `RBH(a,b) ⇔ bestHit(a→G2)=b ∧ bestHit(b→G1)=a` (with best hit = maximum qualifying
similarity + deterministic tie-break), the reciprocity conjunction, the significance/coverage gates,
and the deterministic tie-break are all source-backed. An **independent Python 5-mer Jaccard
recomputation** for every dataset reproduced the asserted values: self-match 1.0 (id/cov 1.0, alignLen
14); superstring case 4/6 = **0.667** (cov 1.0); and the new coverage-gate case
`AAAAACCCCCGGGGG` vs `AAAAACCCCCTTTTT` → shared 6 / union 16 = **0.375**, coverage 6/11 = **0.5455**.

**Stage A note (documented simplification, not a defect).** The ranking metric is a **5-mer Jaccard
similarity, not a BLAST bit-score** (Assumption ASM-01; backed by the alignment-free family, Mash /
Ondov 2016). This affects only *which* near-identical candidate wins a tie; the correctness-critical
parts — the reciprocity rule, the deterministic tie-break, the `≥50%` coverage gate, and the
minimum-similarity gate — are source-backed, and identical sequences score 1.0 so the ranking is
order-preserving on every dataset in this unit.

## Stage B — implementation

Best-hit selection, the reciprocity conjunction, both gates, the deterministic tie-break, and the
**actual-metric reporting** (identity/coverage/alignLen, not placeholders) all match Stage A, verified
by code trace and 13 passing tests. Cross-verification table recomputed vs code: **A≡A** → id 1.0 /
cov 1.0 / alignLen 14 (M1/M3 pass); **A vs A-superstring** → Jaccard 0.667, cov 1.0, and M2 **excludes**
the non-reciprocal b2 (a one-directional impl would wrongly emit b2→a1, so M2 genuinely fails a
non-reciprocal implementation) while M4 rejects at `minIdentity = 1.0`; **AcPrefix vs AcPrefixAlt** →
id 0.375 / cov 0.5455, M7 keeps it at the default `minCoverage 0.5` but rejects at 0.6; **len 4 (< k=5)**
→ similarity 0, S4 yields no pair. `FindOrthologs ≡ FindReciprocalBestHits` by construction.

**Test-quality audit (HARD gate): PASS.** M1/M3 assert exact Jaccard/coverage/alignLen from
independent hand computation; M2 asserts count==1 and b2 absent; M4 isolates the identity gate; **M7
(added this session)** isolates the coverage gate with exact 6/16 and 6/11 values; **S4 (added)** locks
the `< k` similarity-0 edge. No `Greater`/`AtLeast`/range assertions used where an exact value is known;
no weakened assertions, no widened tolerances, no skipped tests, no expected value adjusted to output.
Every public method exercised (RBH + `FindOrthologs` delegate); all Stage-A branches/edges covered —
reciprocity (M2, S2), both gates (M4 identity, **M7 coverage**), empty (M5), null both args (M6),
no-sequence skip (S3), short-seq edge (**S4**), matching (S1), determinism (C1).

## Findings

- **No code defect (State CLEAN).** The implementation already handled the coverage gate and the `< k`
  short-sequence edge correctly; the gap was in the **test surface only**. Two edge-case tests were
  added and logged as a FIXED-NOW finding: **M7** (coverage-gate rejection with exact 6/16 = 0.375 and
  6/11 = 0.5455 values) and **S4** (short-sequence similarity-0).
- **Follow-up (non-blocking, shared with the concept):** if an alignment / bit-score becomes available
  to `Seqeron.Genomics.Analysis`, replace the Jaccard ranking metric so the 50% coverage gate binds as
  an independent constraint (the alignment-free note above).
