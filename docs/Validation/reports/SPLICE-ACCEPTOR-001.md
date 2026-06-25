# Validation Report: SPLICE-ACCEPTOR-001 — Acceptor Site Detection (3' splice site)

- **Validated:** 2026-06-25 (fresh re-validation)   **Area:** Splicing
- **Canonical method(s) (THIS unit's surface):**
  `SpliceSitePredictor.FindAcceptorSites(sequence, minScore, includeNonCanonical)`;
  `SpliceSitePredictor.FindAcceptorBranchPoint(sequence, acceptorAgPosition, minScore)`;
  internal `ScoreAcceptorSite`, `ScoreU12AcceptorSite`, `ScoreBranchPointConsensus`,
  `ComputeBranchPointPptFraction`.
- **Out of scope (separate CLEAN unit SPLICE-MAXENT3-001):** `ScoreAcceptorMaxEnt` /
  the MaxEntScan `score3ss` path and its ME1–ME10 tests — referenced, not re-litigated here.
- **Stage A verdict:** ✅ PASS (re-confirmed; prior PASS-WITH-NOTES notes retained as design notes)
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN — no defect; two test-coverage gaps closed this session (BP9, C3).

## Re-validation context (2026-06-25)

The 2026-06-24 limitation-elimination campaign added explicit branch-point detection
(`FindAcceptorBranchPoint`, Gao et al. 2008 `yUnAy`) to this existing unit and split the
MaxEntScan `score3ss` path into the separate already-CLEAN unit SPLICE-MAXENT3-001, so this
unit was reset to ⬜ pending. This report re-validates the OWN canonical surface (acceptor
consensus + branch point + polypyrimidine tract) fresh against external first sources retrieved
this session. All numbers below were hand-recomputed this session, not echoed from the code/tests.

## Stage A — Description

### Sources opened this session & what they confirm
- **Gao K, Masuda A, Matsuura T, Ohno K (2008) "Human branch point consensus sequence is yUnAy",
  Nucleic Acids Res 36(7):2257–2267** (DOI 10.1093/nar/gkn073), fetched 2026-06-25 via
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2367711/ — confirms VERBATIM:
  - Consensus **`yUnAy`** (branch adenosine at position 0, motif positions −3..+1).
  - Per-position conservation (n = 181 confirmed branch points): **−3 y = 79.0%**, **−2 U = 74.6%**,
    **0 A = 92.3%**, **+1 y = 75.1%** (−1 = n, uninformative).
  - Distance from the 3' splice site: range −50..−5, median −26, mean −27.7 ± 7.6;
    **83% of branch points at −34..−21**.
  - **Polypyrimidine tract spans 4–24 nt downstream of the branch point** (U preferred, +4..+12).
  These match the source-cited constants in the code (`BpPos3PyrimidineConservation=0.790`,
  `…Pos2Uracil=0.746`, `…BranchAdenosine=0.923`, `…Pos1Pyrimidine=0.751`; PPT span 4–24).
- **Wikipedia — RNA splicing** (fetched 2026-06-25): acceptor (3' splice site) consensus
  **`Y-rich-N-C-A-G-[cut]-G`**, intron ends with an "almost invariant AG"; GT-AG rule confirmed;
  branch sequence **`Y-U-R-A-C`** 20–50 nt upstream of the acceptor; polypyrimidine tract (C/U-rich)
  between the branch point and the AG. Consistent with the acceptor PWM `(Y)nNCAG|G`
  (−2 A=1.00, −1 G=1.00, −3 C=0.70, +0 G=0.50, pyrimidine-enriched upstream).

### Feature / coordinate check
- Invariant **AG** at the intron 3' end recognised via `upper[i]=='A' && upper[i+1]=='G'`.
- **Acceptor (AG) vs donor (GU/GC/AU)** — no dinucleotide confusion; U12 acceptor uses AC.
- **Polypyrimidine tract** in `ScoreAcceptorSite`: C/U count in window `[pos−15, pos−3)`,
  contribution `pptScore/12*2` — matches the documented C/U-rich PPT upstream of the AG.
- **Branch-point window 18–40 nt** (`FindAcceptorBranchPoint`): the conservative envelope
  bracketing the Gao −34..−21 core (and within the spliceosome's outer ~18–50 nt envelope).
- **Coordinate convention:** 0-based scan; reported `Position = i+1` = index of the terminal **G**
  of AG (last intronic nt). `FindAcceptorBranchPoint` measures distance from the branch A to that G.

### Independent cross-check (numbers, hand-computed this session)
- **M1** `UUUUUUUUUUUUUUUUCAGGG`: AG at i=17 → Position **18**. PPT window `[2,14)` = 12×U →
  contribution **2.0**; PWM log2 over 8 offsets = **12.21362** raw, count=8;
  normalized `(12.21362/9 + 2)/4` = **0.8392671** → rounds to **0.8393** (test M1 ✓).
- **BP1** `UU CUUAC U×22 AG GGG`: AG at index 29–30 (Position **30**); branch A at index **5**;
  distance = 30 − 5 = **25** (within 18–40); motif `CUUAC`; perfect yUnAy →
  `(0.790+0.746+0.923+0.751)/3.210 = ` **1.0**; PPT (all U) fraction **1.0** (test BP1 ✓).
- **BP4** `AUUAC` (A at −3, purine, loses the 0.790 term): `(0+0.746+0.923+0.751)/3.210` = **0.753894**
  (test BP4 ✓).
- **BP3** all-pyrimidine (no branch A): best `(0.790+0.746+0+0.751)/3.210` = **0.712461** < 0.8 →
  not found at minScore 0.8 (test BP3 ✓).
- **BP9** (added this session) multiple candidates: a perfect `CUUAC` at distance 30 vs a weaker
  `AUUAC` at distance 20 — the detector reports the higher score (1.0) at distance 30, proving
  best-score (not nearest) selection. Hand-verified.

### Findings / divergences (Stage A notes — design notes, not defects)
- **NOTE:** the default acceptor *score* (`ScoreAcceptorSite`) omits a branch-point term; the
  branch point is scored separately (`FindAcceptorBranchPoint`). Biologically sound (a distinct
  element) and not claimed by the default scorer.
- **NOTE:** the `(score/(count+1)+2)/4` normalisation mapping to [0,1] is a documented heuristic,
  not a literature formula; the spec asserts behavioural properties (range, monotonicity, ordering)
  which the tests verify.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs`
- `FindAcceptorSites` (lines ~264–305): guard len<20; scan from i=15; AG/AC recognition;
  `Position=i+1`; motif; U12 gated by `includeNonCanonical`.
- `ScoreAcceptorSite` (~436–471): PPT count `[pos−15,pos−3)` + PWM at `position+2+offset`;
  normalisation to [0,1].
- `ScoreU12AcceptorSite` (~572–616): YCCAC consensus + PPT; /3.5 normalisation.
- `FindAcceptorBranchPoint` (~331–378): scans distance 18→40 (nearest-first); guards
  `acceptorAgPosition<=0`/`>=Length` and `branchA-3<0 || branchA+1>=Length`; keeps the
  strictly-highest-scoring candidate (`score <= best.Score` skip).
- `ScoreBranchPointConsensus` (~480–497): conservation-weighted yUnAy fraction in [0,1].
- `ComputeBranchPointPptFraction` (~505–529): C/U fraction in the 4–24 nt downstream tract,
  clamped to the bases between the branch A and the A of the AG.

### Formula realised correctly? (evidence)
Yes. PWM offset `position + 2 + offset` aligns −2→A, −1→G of the AG, +0→first exonic base.
M1 reproduced to 7 dp (0.8392671). Branch-point conservation constants equal the Gao 2008
table values; `ScoreBranchPointConsensus` divides matched conservation by the maximum (perfect
= 1.0). Distance `agEnd − branchA` and the 18–40 window match Gao's core. PPT fraction over the
4–24 nt downstream span matches Gao's tract definition.

### Cross-verification table recomputed vs code
| Case | Expected (source/hand) | Recomputed/Observed | Match |
|------|------------------------|---------------------|-------|
| M1 position / score | 18 / 0.8393 | 18 / 0.8392671 | ✅ |
| M2 no AG | empty | empty | ✅ |
| M3/M4 guards | empty | empty (null/empty, len<20) | ✅ |
| M5 strong>weak PPT | strong>0.7, weak<0.7 | passes | ✅ |
| M10 two AG | pos 18 & 22, first higher | passes | ✅ |
| S1 U12 YCCAC | pos 17, score 1.0 | 3.5/3.5 = 1.0 | ✅ |
| S3 position | 16 (AG at i=15) | passes | ✅ |
| C1 AG before scan start | empty | empty (scan starts i=15) | ✅ |
| C3 non-ACGT context | 1 site, pos 18 | passes (no throw) | ✅ |
| BP1 canonical | A@5, dist 25, CUUAC, 1.0, PPT 1.0 | matches | ✅ |
| BP4 purine@−3 | 0.753894 | 0.753894 | ✅ |
| BP3 no A @minScore 0.8 | not found (best 0.712461) | not found | ✅ |
| BP5 near edge | 18 found / 17 not | passes | ✅ |
| BP6 far edge | 40 found / 41 not | passes | ✅ |
| BP9 multiple candidates | best (1.0) at dist 30 wins | A@10, dist 30, 1.0 | ✅ |

### Variant/delegate consistency
DNA(T)↔RNA(U) equivalence and case-insensitivity verified for both `FindAcceptorSites`
(M8/M9) and `FindAcceptorBranchPoint` (BP8). U12 AC path gated by `includeNonCanonical` (S1/S2).
Donor path uses GU/GC/AU — no acceptor/donor swap.

### Test quality audit
`SpliceSitePredictor_AcceptorSite_Tests.cs` — the in-scope acceptor + branch-point cases
(M1–M10, S1–S5, C1–C3, BP1–BP9) assert exact sourced positions/scores (0.8393, 1.0, 0.753894,
0.712461 boundary), exact distances, motifs, [0,1] ranges, ordering, threshold subset, window
edges, multiple-candidate selection, non-ACGT robustness, and DNA≡RNA/case. Expected values
trace to Gao 2008 / the acceptor consensus / hand-computation — not code echoes. Deterministic.
(The ME1–ME10 MaxEnt cases in the same file belong to SPLICE-MAXENT3-001 and were not re-litigated.)
Full solution suite: **18792 passed / 0 failed**; the in-scope filter: **37 passed / 0 failed**.

### Gaps closed this session
- **BP9** added — best-candidate (highest score, not nearest) selection in the branch-point
  window; the prior suite had no multiple-candidate test.
- **C3** added — non-ACGT characters in the acceptor context are tolerated (no throw, AG still found).

### Findings / defects
None. No off-by-one, no dinucleotide confusion, no missing spec-claimed feature.

## Verdict & follow-ups
- **Stage A:** ✅ PASS — acceptor consensus, branch-point `yUnAy` motif + conservation, 18–40 nt
  window, and 4–24 nt PPT span independently confirmed against Gao et al. (2008) and Wikipedia
  (RNA splicing) fetched this session.
- **Stage B:** ✅ PASS — code faithfully realises the validated description; M1/BP1/BP4/BP3/BP9
  reproduced by hand; two coverage gaps (BP9, C3) closed.
- **State:** ✅ CLEAN — no defect found; only test additions, no production-code change.
