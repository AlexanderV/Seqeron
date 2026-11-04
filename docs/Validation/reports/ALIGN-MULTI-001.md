# Validation Report: ALIGN-MULTI-001 — Multiple Sequence Alignment (star + progressive + consistency, full re-validation)

- **Validated:** 2026-06-24   **Area:** Alignment
- **Canonical method(s):**
  - `SequenceAligner.MultipleAlign(IEnumerable<DnaSequence>, ScoringMatrix?)` — star / center-star, SP score
  - `SequenceAligner.MultipleAlignProgressive(IEnumerable<DnaSequence>, ScoringMatrix?)` — UPGMA guide-tree progressive (Feng-Doolittle / Clustal-style)
  - `SequenceAligner.MultipleAlignConsistency(IEnumerable<DnaSequence>, ScoringMatrix?)` — T-Coffee primary library + triplet extension + library-scored progressive DP
  - (sibling `MultipleAlignIterative` also exercised for additivity; covered by its own addendum)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

This re-validation independently re-confirms **all three** MSA variants named in the unit (star,
progressive, consistency), with emphasis on the newest variant `MultipleAlignConsistency`
(commit `5b7a9f37`) not covered by the prior 2026-06-17 report. The implementer's code and tests were
treated as untrusted: theory was checked against the actual papers (the T-Coffee PDF was opened and
read page-by-page), and all cross-check numbers were recomputed against the live code via a temporary
probe (since removed), not read out of the repo's own assertions.

---

## Stage A — Description

### Sources opened this session (URLs recorded)

| Source | URL | What it confirms |
|--------|-----|------------------|
| Wikipedia "UPGMA" | https://en.wikipedia.org/wiki/UPGMA | Proportional size-weighted update **d((A∪B),X) = (|A|·d(A,X) + |B|·d(B,X)) / (|A|+|B|)**; "unweighted" = equal weight to the original distances. Matches the code's averaging verbatim. |
| Notredame, Higgins & Heringa (2000) T-Coffee, J Mol Biol 302:205-217 (full PDF, pp.208-210 read) | https://web.stanford.edu/class/gene211/pdfs/Notredame-Tcoffee.pdf | GARFIELD worked example (Fig 2 + p.209 text): direct A(G)/B(G) primary weight **88** = percent identity of the A–B pair; through C, W1=W(A(G),C(G))=**77**, W2=W(C(G),B(G))=**100**, triplet = min(77,100)=**77**; "added to the previous one to give a total weight of **165** (i.e. 77+88)"; through D = 0 (uninformative); "the weight … will be the **sum of all the weights** gathered through … all the triplets." p.210: progressive DP is Gotoh (1982) with "gap-opening penalties and gap-extension penalties set to zero"; group-vs-group: "the **average** library scores in each column of existing alignment are taken." |
| Wikipedia "Multiple sequence alignment" / "Clustal" (prior session, re-checked) | https://en.wikipedia.org/wiki/Multiple_sequence_alignment | MSA definition (length L ≥ max nᵢ, no all-gap column, reversibility); SP = "sum of all the pairs of characters at each position"; progressive = guide tree (NJ/UPGMA) + most-similar-first + once-a-gap. |
| Feng & Doolittle (1987), J Mol Evol 25:351-360 (citation; method) | PubMed 3118049 | Progressive alignment along a guide tree; "once a gap, always a gap" preservation across fusions. |

### Formula / convention check

- **Star:** center selection (max total similarity to others) → pairwise NW to center → gap reconciliation
  into one MSA coordinate space → majority consensus → SP score. SP convention: match/mismatch from
  matrix, residue-gap = GapExtend, **gap-gap = 0** (standard, documented). ✔
- **Progressive distance:** d = 1 − fractional identity (identical non-gap columns / pairwise-NW length).
  Standard Clustal-style conversion. ✔
- **UPGMA:** smallest distance merged first, lowest-index tie-break, **proportional size-weighted**
  averaging — matches the Wikipedia formula symbol-for-symbol. ✔
- **Once a gap, always a gap:** existing columns copied verbatim; merges insert only whole all-gap
  columns; inserted gaps never residue-filled (leaf and profile level). Matches Feng-Doolittle. ✔
- **Consistency primary library:** each aligned residue pair weighted by pairwise **percent identity**
  (integer, ×100), global (NW) + local (SW) combined by **signal addition** (sum of weights for
  duplicate pairs). Matches T-Coffee p.207. ✔
- **Consistency extension:** extended = direct primary + Σ over intermediates Sₖ of **min(W1,W2)**;
  uninformative triplets contribute 0 ⇒ extension never lowers a weight. Matches T-Coffee p.209
  (GARFIELD 88 → 165). ✔
- **Consistency progressive DP:** column score = library weight of cross-profile residue pairs,
  **zero gap penalty**; once-a-gap enforced. Matches T-Coffee p.210. ✔

### Edge-case semantics

Empty → `Empty`; single → verbatim, SP 0; null → `ArgumentNullException`; two empty → distance 0;
empty vs non-empty → distance 1; k=2 consistency → extended = primary (no triplets) = pairwise global;
identical inputs → gap-free exact MSA. All defined and sourced/standard, consistent across all variants. ✔

### Documented minor divergence (carried, not a defect)

The T-Coffee paper says group-vs-group columns use the **average** library score; the code **sums** the
cross-profile pair weights (`ColumnLibraryScore`). The evidence doc records this as an assumption: the
sum is the same objective up to a per-merge constant scaling, so the DP's argmax (the chosen alignment)
is unaffected for fixed group sizes within a single merge. The reported `TotalScore` is deliberately the
SP score (for cross-aligner comparability), not the consistency objective — also documented. Sound; PASS.

**Stage A verdict: PASS.** All three abstract methods (center-star SP; distance→UPGMA→profile NW with
once-a-gap; T-Coffee primary-library + min-triplet extension + library-scored zero-gap progressive DP)
match the sources opened this session.

---

## Stage B — Implementation

### Code path reviewed (`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs`)

- **Star:** `MultipleAlign` L702; `SelectCenterSequence` L791; gap reconciliation L902-1044;
  `BuildConsensus` L1046 (majority, gap participates, nucleotide-over-gap tie); `ComputeSumOfPairsScore`
  L1078 (gap-gap=0, residue-gap=GapExtend). ✔
- **Progressive:** `MultipleAlignProgressive` L1170; `PairwiseIdentityDistance` L1232 (1−id);
  `BuildProgressiveGuideTree` L1275 (UPGMA, proportional averaging L1331 = `(dIK·sizes[i]+dJK·sizes[j])/newSize`,
  lowest-index tie-break); `AlignProfileSubtree`/`AlignProfiles` L1362/1385 (column NW; gaps as new
  columns only). ✔
- **Consistency:** `MultipleAlignConsistency` L1918; `BuildExtendedLibrary` L1990 (primary global+local
  signal-add L2004-2010; triplet extension `extended[key] += min(w1,w2)` over intersected neighbours
  L2036-2046); `PercentIdentity` L2116 (`round(100·identical/len)`); `ColumnLibraryScore` L2299 (sum of
  cross-profile extended weights, gap=0); `AlignConsistencyProfiles` L2203 (zero-gap NW over columns,
  once-a-gap via column-map append). ✔

### Independent hand/probe derivations (verified against the CODE, then probe removed)

**Star vs progressive discriminator `["AAGAA","AACAA","GGTGG","GGTGG"]`** — probe output:
- STAR: `AAG-AA / -AACAA / -GGTGG / -GGTGG`, len 6, SP **−13**.
- PROG: `AAGAA / AACAA / GGTGG / GGTGG`, len 5, gap-free, SP **−12**.
- Methods genuinely differ; progressive is the better (a gap can only cost here). ✔

**Progressive `["ACGT","ACGT","AGT"]`** — probe: rows `ACGT / ACGT / A-GT`, SP **8**
(col0 AAA +3, col1 {C,C,−} −1, col2 GGG +3, col3 TTT +3 = 8). Hand-derivation matches code. ✔

**Consistency `["ACGT","ACGT","AGT"]`** — probe: degaps to inputs, rows `ACGT / ACGT / A-GT`, SP **8**;
consistency reaches the same MSA as progressive on this easy case. ✔

**Consistency library, distinct sequences `["ACGT","ACGT","ACGA"]`** (probe via internal `GetLibraryWeights`):
- Pair (S0.0='A', S1.0='A'): primary = **200** (global gapless 100%id=100 + local SW whole-seq 100%id=100,
  signal-added), extended = **375** = 200 + min-triplet 175 through S2.0 (where S2's global pos0 'A' weight 75
  + local "ACG" 100 = 175). Independently reproduces the GARFIELD relation `extended = primary + Σ min(W1,W2)`
  and `extended > primary` for a consistency-supported pair. ✔
- Pair (S0.3='T', S1.3='T') primary 200, extended 275 — strictly less third-party support (S2.3='A'
  mismatches and the local segment "ACG" excludes pos3). Supported pair > less-supported pair. ✔

These independently confirm the paper's worked numbers (88→165 ⇒ here 200→375) over the DNA alphabet —
the integer relation `extended = primary + Σ min-triplet` is alphabet-independent.

### MSA correctness traps (probes + tests)

Equal-length rows ✔; degap recovers each input ✔; no all-gap column ✔; deterministic (byte-identical on
repeat) ✔; once-a-gap at leaf and profile level ✔; k=2 consistency = pairwise global (length match) ✔.

### Test-quality audit

- **Star/progressive** tests (`SequenceAligner_MultipleAlign_Tests.cs`, `…Progressive_Tests.cs`): exact
  rows + hand-computed SP (SP=8, SP=24, the discriminator), strict, mutation-checked per the prior report.
- **Consistency** tests (`SequenceAligner_MultipleAlignConsistency_Tests.cs`, 12 tests): TM04 pins the
  GARFIELD relation with exact integers (200→400 for identical seqs; my distinct-seq probe gives 200→375,
  same relation); the consistency objective (TM08) is recomputed independently via `GetLibraryWeights`,
  not echoed from the DP; TM05 supported>unsupported; TM06 identical→gapless; TM07/TS02 validity
  invariants; TM09 k=2=global; TM10 determinism; TS01 sibling additivity. Assertions check exact
  sourced values, not "no-throw". Not green-washed.
- All 96 MultipleAlign-family tests pass (star + progressive + iterative + consistency + properties +
  benchmark).

### Variant / additivity consistency

All four aligners reuse `GlobalAlignCore`/`LocalAlignCore`, `BuildProgressiveGuideTree`, `BuildConsensus`,
`ComputeSumOfPairsScore`. The consistency aligner adds no perturbation: star/progressive/iterative remain
byte-for-byte unchanged (TS01 + full green suite).

**Stage B verdict: PASS.** The code faithfully realises the three validated descriptions; every
hand/probe-derived value matches the live code; tests are strict, deterministic, and not tautological.

---

## Verdict & follow-ups

- **Stage A: PASS · Stage B: PASS · End-state: ✅ CLEAN.**
- **No code changed.** Only a temporary probe test was added and removed.
- One documented minor divergence (group columns summed vs the paper's "average"; same argmax; SP reported
  as TotalScore for comparability) — already recorded as an assumption in the evidence doc; not a defect.
- **Full unfiltered suite: 18208 passed, 0 failed** (build 0 warnings / 0 errors).
- All three MSA variants are correct and optimal on the derived discriminating cases; NOT LIMITED.

### Sources (this session)
1. Notredame C, Higgins DG, Heringa J (2000) "T-Coffee." J Mol Biol 302:205-217 — full PDF pp.208-210 read.
2. Wikipedia "UPGMA" — https://en.wikipedia.org/wiki/UPGMA
3. Wikipedia "Multiple sequence alignment" / "Clustal".
4. Feng DF, Doolittle RF (1987) J Mol Evol 25:351-360. PubMed 3118049.


---

# Re-validation 2026-11-04 (campaign 2026-09, batch B09) — supersedes the verdicts above

- **Stage A verdict:** PASS-WITH-NOTES · **Stage B verdict:** PASS-WITH-NOTES after fixes (four defects found and fixed, see below) · **End-state:** ✅ FIXED (declared deviations remain in the algorithm doc §5.3).

## Stage A — description

Sources opened this time: the T-Coffee paper PDF (tcoffee.org, pp.207-210), the T-Coffee sources `evaluate.c` / `aln_convertion_util.c`, ClustalW 1.83 `pairalign.c`, Abouelhoda & Ohlebusch 2005 (PDF), Gusfield 1993 and Higgins & Sharp 1988 (bibliographic records via WebSearch). Findings against the previous description: (1) the "star" method was described as the Gusfield/Wikipedia center-star but its pairwise stage was a suffix-tree anchored heuristic, so the claim "pairwise alignments of a center-star" did not hold; (2) the Clustal-style distance claimed in the progressive description is not the ClustalW distance (ClustalW divides identities by the *shorter sequence length*, not by the alignment length); (3) T-Coffee extension was described as applying to "each pair of residues", implemented only for pairs of the primary library; (4) the T-Coffee weight "percent identity" was taken over the full alignment length (including gap columns) whereas the reference implementation counts aligned residue pairs only.

## Stage B — implementation, defects and fixes

| # | Defect | Evidence | Fix |
|---|---|---|---|
| F4 | `MultipleAlign` (k = 2) was not an optimal pairwise alignment: anchors constrain the path | 4/2000 random related pairs; e.g. SP 9 vs Biopython 24, 22 vs 27, −35 vs −3 | exact center-star (Gusfield): canonical NW for every pair, center = best total pair score; anchored variant kept as `MultipleAlignAnchored` |
| F5 | `AnchorBasedAligner`: overlapping anchors discarded by a greedy filter, LIS maximised the *count* of anchors, one-sided gaps scored affinely while DP gaps were linear (inconsistent `Score`), > 5000-base segments aligned by a diagonal heuristic, comments claimed "all MEMs" | unit tests (chain of overlapping anchors 10 vs 16 bases; Score ≠ re-scored rows) | weighted colinear chain with overlap trimming, linear scoring everywhere, exact Hirschberg for big gaps, comments corrected |
| F6 | `SelectCenterSequence` counted 4-mers with a private hash | B06 request | Core `CountKmersSpan` (identical counts) |
| F7 | T-Coffee extension only for primary pairs and percent identity over gap columns | 170/300 rows differ from the independent T-Coffee reference | on-demand extension for every residue pair; weight = ⌊100·identical/aligned⌋ |
| F8 | Progressive distance not the ClustalW distance; UPGMA ties decided by floating-point noise | `pairalign.c`; 1/300 random sets differ between exact-rational and double ties | identities / min length; ties within 1e-9 |

Cross-verification: independent Python implementations of all five methods reproduce the rows of 300 random sets exactly (see Evidence addendum); MCP tools added for every method; GapOpen is documented as unused by the linear-gap MSA family.

## Verdict

PASS-WITH-NOTES / PASS-WITH-NOTES, FIXED. Full fast tier: see `docs/Validation/review-2026-09/B09.md`.
