# Validation Report: ALIGN-MULTI-001 — Multiple Sequence Alignment (Progressive guide-tree, re-validation)

- **Validated:** 2026-06-17   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.MultipleAlignProgressive(IEnumerable<DnaSequence>, ScoringMatrix?)` (NEW, this re-validation); `SequenceAligner.MultipleAlign(...)` (star, unchanged)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

This is a **Phase-3 independent re-validation** of the enhancement in commit `7d0b2cb`, which added a
guide-tree progressive (Feng-Doolittle / Clustal-style) aligner beside the existing star MSA. The
validator session is independent of the implementer session; the implementer's code AND tests were
treated as untrusted and verified against externally-sourced theory and against hand derivations,
**not** against the repo's existing assertions.

---

## Stage A — Description

### Sources opened this session (URLs recorded)

| Source | URL | What it confirms |
|--------|-----|------------------|
| Wikipedia "Multiple sequence alignment" §Progressive alignment | https://en.wikipedia.org/wiki/Multiple_sequence_alignment | Guide tree built by "an efficient clustering method such as neighbor-joining or UPGMA"; "combining pairwise alignments beginning with the most similar pair and progressing to the most distantly related"; once-a-gap: "once a sequence has been aligned into the MSA, its alignment is not considered further" → early errors "propagated through to the final result". |
| Wikipedia "UPGMA" | https://en.wikipedia.org/wiki/UPGMA | Smallest-distance merge; proportional size-weighted update `d((A∪B),X) = (|A|·d(A,X) + |B|·d(B,X)) / (|A|+|B|)`; "unweighted" = equal weight to original distances. |
| Feng & Doolittle (1987), J Mol Evol 25:351-360, PubMed 3118049 | (citation; method) | Progressive alignment along a guide tree; the "once a gap, always a gap" preservation rule across subsequent fusions. |

### Formula / convention check

- **Distance:** d = 1 − fractional identity, identity = identical (non-gap, equal) columns / pairwise-NW
  alignment length. Standard Clustal-style identity-to-distance conversion. ✔
- **Guide tree:** UPGMA, smallest distance merged first, proportional cluster averaging — matches the
  Wikipedia UPGMA formula verbatim (size-weighted). ✔
- **Sequence order:** most-similar-first via post-order traversal of the UPGMA tree. ✔
- **Once a gap, always a gap:** existing columns copied verbatim; merges only insert whole new all-gap
  columns; an inserted gap is never later filled with a residue. Matches Feng-Doolittle. ✔
- **Profile-profile score:** sum-of-pairs over all cross-profile residue pairs (average); gap-gap
  neutral, residue-gap = GapExtend. Standard SP profile scoring. ✔
- **SimpleDna scoring:** Match +1, Mismatch −1, GapExtend −1 (confirmed at source lines 22-26). ✔

### Edge-case semantics

Empty → `Empty`; single → verbatim, SP 0; null → `ArgumentNullException`; two empty seqs → distance 0;
one empty vs non-empty → distance 1. All defined and sourced/standard. ✔

**Stage A verdict: PASS.** The abstract method (distance → UPGMA → progressive profile NW with
once-a-gap) is biologically and mathematically correct and matches the sources opened this session.

---

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs`:
- `MultipleAlignProgressive` L1143-1197 — pipeline + reprojection to input order.
- `PairwiseIdentityDistance` L1205-1227 — d = 1 − identical/length.
- `BuildProgressiveGuideTree` L1248-1318 — UPGMA, lowest-index tie-break, proportional averaging.
- `AlignProfileSubtree` L1335-1349 — post-order recursion.
- `AlignProfiles` L1358-1425 — column-wise NW + traceback; gaps inserted as new columns only.
- `AppendReprojectedRows` L1432-1441 — existing columns copied verbatim, −1 → fresh '-'.
- `ProfileColumnScore` / `ProfileColumnVsGap` L1444-1483 — SP profile scoring.
- Shared `GlobalAlignCore` L220, `BuildConsensus` L1019, `ComputeSumOfPairsScore` L1051.

### Independent hand-derivations (verified against the CODE, not the repo's tests)

**Case 1 — M07 `["ACGT","ACGT","AGT"]` (distances → merge order → exact columns):**
- Pairwise NW: d(0,1)=0 (identical); NW("ACGT","AGT") = "ACGT"/"A-GT", 3 identical of 4 → d=0.25.
- UPGMA: smallest distance 0 ⇒ merge (0,1) → profile [ACGT/ACGT]; then add leaf "AGT".
- Profile-profile NW places the single gap opposite the 'C' column ⇒ third row "A-GT".
- **Exact rows:** `ACGT / ACGT / A-GT`. SP = col0 AAA(+3) + col1 {C,C,−}(+1−1−1=−1) + col2 GGG(+3) +
  col3 TTT(+3) = **8**. Independent probe output: rows `ACGT/ACGT/A-GT`, SP 8. ✔

**Case 2 — discriminating progressive ≠ star `["AAGAA","AACAA","GGTGG","GGTGG"]`:**
- Distances: d(2,3)=0, d(0,1)=0.2, all cross-cluster ≥0.8. UPGMA merges (0,1) then (2,3) then the two
  clusters; each cluster is internally gap-free (equal length, no indel improves the diagonal), and the
  profile-profile step keeps the ungapped diagonal optimal.
- **Progressive (probe):** `AAGAA / AACAA / GGTGG / GGTGG`, length 5, gap-free, SP **−12**.
- **Star (probe):** `AAG-AA / -AACAA / -GGTGG / -GGTGG`, length 6, SP **−13**.
- The two methods genuinely differ, and the progressive answer is the better (a gap can only cost here),
  so it is the correct/optimal alignment. The discriminating test's exact-row assertion **fails if run
  against the star aligner** (proven by S01 which locks the star's length-6 gapped rows). ✔

**Case 3 — profile-level once-a-gap `["ACGT","ACGT","AGT","AGT"]` (validator-added M12):**
- UPGMA: d(0,1)=0 and d(2,3)=0 minimal; lowest-index tie-break merges (0,1) first, then (2,3), then the
  two profiles. P1=[ACGT/ACGT] vs P2=[AGT/AGT] inserts one new all-gap column in P2 opposite 'C'.
- **Exact rows:** `ACGT / ACGT / A-GT / A-GT`, SP = col0 AAAA(+6) + col1 {C,C,−,−}(−3) + col2 GGGG(+6) +
  col3 TTTT(+6) = **15**. Probe confirmed rows + SP 15. Gap-removal recovers all four inputs ⇒ no
  inserted gap was ever residue-filled, at the profile level too. ✔

### MSA correctness traps checked (independent probes)

- **Equal-length rows:** ✔ (all probes).
- **Gap-removal recovers each input:** ✔ (M07, M12, INSERT `GATTACA×2 + GATTTACA`, GAPPY 5-seq).
- **No all-gap column:** ✔ (GAPPY 5-seq, M10).
- **Deterministic tie-breaks:** ✔ (same input twice → byte-identical rows).
- **Once a gap, always a gap (leaf & profile level):** ✔ — gap columns carried verbatim; inserted gaps
  never filled.

### Test-quality audit (HARD gate)

- Implementer tests M01-M11 + S01: 12 strict tests. Exact aligned strings where derivable (M04-M08),
  exact SP scores hand-stated, the discriminating M08 + the additivity S01 together lock progressive
  *and* star byte-for-byte. The discriminator is real: asserting the gap-free length-5 progressive rows
  fails against the star aligner.
- **Validator addition — M12** `MultipleAlignProgressive_TwoGappedProfilesMerge_GapsCarriedVerbatim_ExactColumns`:
  closes the one genuine gap (no prior test pinned exact columns when *two already-gapped profiles* merge
  — the core profile-level once-a-gap guarantee). Exact rows + SP 15, independently hand-derived.
- **Mutation check:** mutating `AppendReprojectedRows` to fill an inserted (−1) column with a residue
  instead of '-' (a once-a-gap violation) → **4 tests fail** including M12. The suite is not green-washed.

### Variant / additivity consistency

Star `MultipleAlign` signature and output are byte-for-byte unchanged (S01 locks
`AAG-AA/-AACAA/-GGTGG/-GGTGG`; full suite green). Both aligners reuse `GlobalAlignCore`,
`BuildConsensus`, `ComputeSumOfPairsScore` consistently.

**Stage B verdict: PASS.** The code faithfully realises the validated description; all hand-derived
cases match the actual code output; tests are strict, deterministic, mutation-checked, and not
tautological. No defect found.

---

## Verdict & follow-ups

- **Stage A: PASS · Stage B: PASS · End-state: ✅ CLEAN.**
- No defect; nothing logged in FINDINGS_REGISTER (the only change is a strict added test, M12).
- Full unfiltered suite: **6772 passed, 0 failed**. Build 0 errors.
- The progressive aligner is correct (not garbage / non-optimal on the cases derived); it is NOT LIMITED.

### Sources (this session)
1. Wikipedia "Multiple sequence alignment" — https://en.wikipedia.org/wiki/Multiple_sequence_alignment
2. Wikipedia "UPGMA" — https://en.wikipedia.org/wiki/UPGMA
3. Feng DF, Doolittle RF (1987) "Progressive sequence alignment as a prerequisite to correct
   phylogenetic trees." J Mol Evol 25(4):351-360. PubMed 3118049.
