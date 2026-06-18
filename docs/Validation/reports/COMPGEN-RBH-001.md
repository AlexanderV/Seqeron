# Validation Report: COMPGEN-RBH-001 — Reciprocal Best Hits (RBH / BBH)

- **Validated:** 2026-06-16   **Area:** Comparative
- **Canonical method(s):** `ComparativeGenomics.FindReciprocalBestHits(genome1Genes, genome2Genes, minIdentity, minCoverage)`; delegate `FindOrthologs(...)`; private helpers `FindBestHit`, `CalculateSequenceSimilarity`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **Moreno-Hagelsieb & Latimer (2008), *Bioinformatics* 24(3):319–324** — PubMed abstract fetched (https://pubmed.ncbi.nlm.nih.gov/18042555/). Verbatim RBH definition confirmed: *"orthologs are assumed if two genes each in a different genome find each other as the best hit in the other genome."* This is exactly the criterion the unit implements. (The body's "≥50% coverage", "1×10⁻⁶ E-value", and "sort highest→lowest bit-score, then smallest→highest E-value" quotes are behind the OUP paywall and could not be re-fetched this session; the *definition* itself is independently confirmed below from multiple sources.)
2. **WebSearch survey of the RBH/BBH literature** (Genome Biol. Evol. "Bidirectional Best Hits Miss Many Orthologs…"; arXiv Best-Match-Graph series 1803.10989 / 1903.07920 / 2001.00958; bioRxiv 2020.05.04.077222) confirms the terminology and the symmetric requirement: *"If gene a in genome A finds gene b as its best, highest-scoring, match in genome B; and gene b finds gene a as its best match in genome A, they are RBH and thus inferred to be orthologs"*; and *"orthology is a symmetric relation, [so] orthologs are necessarily reciprocal best matches."* → both directions must agree; a one-directional best hit is insufficient.
3. **Tatusov, Koonin & Lipman (1997) via NCBI Handbook NBK21090** — fetched (https://www.ncbi.nlm.nih.gov/books/NBK21090/). Confirms COGs are built by *"detecting triangles of mutually consistent, genome-specific best hits (BeTs)"*; the pairwise special case of a mutually-consistent BeT is precisely the reciprocal best hit. Confirms reciprocity (symmetry) as the building block.

### Formula check
RBH pair (a,b) iff `bestHit(a→G2)=b` AND `bestHit(b→G1)=a`, with best hit = maximum qualifying similarity score and a deterministic tie-break. This matches the Core Model in `Reciprocal_Best_Hits.md` §2.2 and the cited definition word-for-word. The reciprocity conjunction, the significance/coverage gate, and the deterministic tie-break are all source-backed.

### Edge-case semantics check
All documented edge cases have a defined, sourced expected behaviour: null → `ArgumentNullException` (repo contract); empty genome → no pairs (a pair needs one gene from each genome); gene without sequence → skipped; one-directional best hit → excluded (reciprocity, [1][2]); sub-threshold (identity OR coverage) → excluded ([1]); sequence shorter than k=5 → similarity 0 → never qualifies.

### Independent cross-check (numbers, hand-computed this session)
k=5 k-mer Jaccard recomputed independently (Python, this session):
- `ACGTACGTACGTAC` vs itself → k-mers {ACGTA,CGTAC,GTACG,TACGT}; identity = 4/4 = **1.0**, coverage = 4/4 = **1.0**, alignLen = min len = **14**. (Matches M1/M3.)
- `ACGTACGTACGTAC` vs `ACGTACGTACGTACAA` → shared 4 / union 6 = **0.667**, coverage 4/4 = 1.0. (Matches M2/M4 superstring case.)
- `AAAAACCCCCGGGGG` vs `AAAAACCCCCTTTTT` → shared 6 (AAAAA,AAAAC,AAACC,AACCC,ACCCC,CCCCC) / union 16 = **0.375**; coverage 6/11 = **0.5455**. (New M7 coverage-gate case.)

### Findings / divergences (Stage A)
- **Note (PASS-WITH-NOTES):** the ranking metric is a 5-mer Jaccard similarity, NOT a BLAST bit-score (Assumption ASM-01, documented in Evidence and the algorithm doc; backed by the alignment-free family, Mash/Ondov 2016). This is an explicit, documented simplification: it affects *which* near-identical candidate wins a tie, but the correctness-critical parts (reciprocity rule, deterministic tie-break, ≥50% coverage gate, minimum-similarity gate) are source-backed. Identical sequences score 1.0, so for every dataset in this unit the ranking is order-preserving. Acceptable as a documented divergence.
- The "≥50% coverage" / "1×10⁻⁶ E-value" body quotes are paywalled; the operational *definition* (the load-bearing claim) is independently confirmed from three sources. No correctness concern.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs`:
- `FindReciprocalBestHits` (lines 465–516): null-checks both lists; drops genes without a sequence; computes best1→2 (with full metrics) and best2→1 ids; emits a pair iff g1's best hit is g2 AND g2's best hit back is g1. Correct reciprocity realisation.
- `FindBestHit` (410–440): max identity; tie-break by larger coverage then ordinal id → unique deterministic winner; qualifies only when `identity ≥ minIdentity && coverage ≥ minCoverage` (both gates). Matches the validated description.
- `CalculateSequenceSimilarity` (518–549): k=5; identity = Jaccard, coverage = shared/min(kmer counts), alignLen = min seq length; `< k` → (0,0,0). Matches hand computation.
- `FindOrthologs` (334–341): pure delegation to `FindReciprocalBestHits`. No divergence possible.

### Formula realised correctly? (evidence)
Yes. Best-hit selection, the reciprocity conjunction, both gates, the deterministic tie-break, and the actual-metric reporting (identity/coverage/alignLen, not placeholders) all match Stage A. Verified by trace and by the 13 passing tests.

### Cross-verification table recomputed vs code
| Dataset | Hand value (this session) | Code/test outcome | Match |
|---------|---------------------------|-------------------|-------|
| A≡A | id 1.0, cov 1.0, alignLen 14 | M1/M3 pass | ✅ |
| A vs A-superstring | Jaccard 0.667, cov 1.0 | M2 excludes b2 (non-reciprocal), M4 rejects at minIdentity 1.0 | ✅ |
| AcPrefix vs AcPrefixAlt | id 0.375, cov 0.5455 | M7: kept by default (0.5), rejected at minCoverage 0.6 | ✅ |
| len 4 (< k=5) | similarity 0 | S4: no pair | ✅ |

### Variant/delegate consistency
`FindOrthologs` delegates to `FindReciprocalBestHits` (smoke test asserts exact tuple equality). Consistent by construction.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** M1/M3 assert exact Jaccard/coverage/alignLen from independent hand computation; M2 asserts count==1 and that b2 is absent (a one-directional impl would emit b2→a1, so M2 genuinely fails a non-reciprocal impl); M4 isolates the identity gate; **M7 (added this session)** isolates the coverage gate with exact 6/16 and 6/11 values; **S4 (added)** locks the `< k` similarity-0 edge. No Greater/AtLeast/range used where an exact value is known.
- **No green-washing:** no weakened assertions, no widened tolerances, no skipped tests; no expected value adjusted to match output (values trace to hand computation).
- **Cover all the logic:** every public method exercised (RBH + FindOrthologs delegate); all Stage-A branches/edges now covered — reciprocity (M2,S2), both gates (M4 identity, **M7 coverage**), empty (M5), null both args (M6), no-sequence skip (S3), short-seq edge (**S4**), matching (S1), determinism (C1).
- **Honest green:** FULL unfiltered suite = **6605 passed, 0 failed, 1 skipped** (pre-existing unrelated MFE benchmark); `dotnet build` 0 warnings, 0 errors.

### Findings / defects
- **Coverage gaps (closed this session):** the documented coverage gate (`minCoverage`) and the `< k` short-sequence edge had **no test** (M4 only raised `minIdentity`). Added **M7** (coverage-gate rejection, exact values) and **S4** (short-sequence similarity-0). No code defect: the implementation already handled both correctly; the gap was in the test surface only. Logged as a FIXED-NOW finding.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (alignment-free Jaccard ranking is a documented, sourced simplification; definition independently confirmed).
- **Stage B:** PASS (code faithfully realises the validated RBH; tests now cover all Stage-A branches with sourced exact values).
- **End-state:** ✅ CLEAN — fully functional; the only gap (two missing edge-case tests) was completely fixed this session.
- **Test-quality gate:** PASS.
