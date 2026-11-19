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

## B09 current state (2026-11-19) - supersedes the original verdict lines above and below

- **Stage A:** PASS-WITH-NOTES. **Stage B:** PASS-WITH-NOTES. **State:** FIXED, with DOABLE items still open (see Open items); not CLEAN.
- **Fixes (B09.md, F-numbers):** F4-F10 (first pass: exact k = 2 star, anchored aligner, T-Coffee extension, ClustalW distance, MCP tools), F48, F51-F56, F58, F60-F76 (audit round 1: exact ClustalW 2.1 pairwise / multiple / protein stages, MUSCLE 3.8.1551, T-Coffee 13.41, affine MSA, NJ guide tree, exact sum-of-pairs MSA, protein / matrix MSA); audit round 2: F81-F84 (input contract, size guards, parallel exceptions), F89 (MUSCLE round 3), F92 / F101 (ClustalW `-TYPE` guess, `-CLUSTERING=UPGMA`: 600 / 600 user-run sets), F93 (T-Coffee paper-path library options), F100 (MCP stage tools), F102 (`-ITERATION` / `-NUMITER`: 700 / 700), F103 (`-PROFILE` / `-SEQUENCES`: 1 200 / 1 200).
- **Oracle numbers:** ClustalW 2.1: 3 600 pairs (full DP), 8 257 pairs (quick), 1 000 DNA sets, protein pair stages 4 204 + 3 774, 240 + 800 matrix-file sets, 340 sets x 4 end-gap modes, 360 matrix sets; T-Coffee 13.41: 800 progressive sets, 820 libraries + 1 480 trees, 430 RNA sets, 4 715 option sets (4 712 equal); ClustalW 2.1 UPGMA 600 / 600 user-run sets, `-ITERATION` 700 / 700, `-PROFILE` / `-SEQUENCES` 1 200 / 1 200; MUSCLE 3.8.1551 (round 1 420 / 420 and round 2 630 / 630 user-run Docker sets; round 3 60 / 60 run locally with the Debian muscle3 package + a fresh sample of 1 800 / 1 800); exact MSA = two independent Python programs + exhaustive enumeration.
- **Open items:** A2-14b / c / d (ClustalW guide-tree files `-NEWTREE` / `-USETREE`, `-OUTORDER=ALIGNED`, secondary-structure gap-penalty masks; UPGMA with `-ITERATION` awaits a user-run `-O0` binary). (A1-10c a3 done by A2-1, F89: 60 / 60 locally run sets of round 3 equal, `-matrix` for nucleotides implemented, fresh sample 1 800 / 1 800); A1-13d (rest): fasta heuristics / `fast_pair` / `-quicktree` ports (DOABLE); A1-13d BLOCKED methods needing external programs / databases (blast, 3D-coffee, sap) - proof in B09.md. The two 2026-06 'CLEAN' statements below describe the first pass only.
- The unit's TestSpec (`tests/TestSpecs/ALIGN-MULTI-001.md`), Evidence (`docs/Evidence/ALIGN-MULTI-001-Evidence.md`) and algorithm document carry the same state; where the original text below says CLEAN / no defects / no code changes, that was the state of the first pass only.

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
- **Progressive distance:** d = 1 − identities / min(len1, len2) (the ClustalW `pairalign.c` distance; corrected from the first-pass "1 − fractional identity (identical columns / alignment length)" in F8, audit A1-17: the tool description, comments and docs now state the exact formula).
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


## Audit r1 addendum (A1-15, B09 F51)

`MultipleAlignProgressive/Iterative/Consistency` accept `MsaOptions.Distance`; `ClustalWAffine` takes the guide-tree distance from ClustalW 2.1's own pairwise stage (`SequenceAligner.ClustalWPairwiseAlign`, an exact port of `FullPairwiseAlign`: affine local Smith-Waterman-Gotoh + Myers-Miller, IUB ×100, gaps 15/6.66). 3 600 pairs equal `clustalw -ALIGN` "Score:"; progressive/iterative rows equal an independent Python implementation on 140 sets. Default stays `Linear`. Tests: `SequenceAligner_ClustalWDistance_Tests` (17) + 2 MCP tests.


## Audit r1 addendum (A1-10, B09 F52)

Additive `MsaOptions.GapModel = MsaGapModel.Affine` on all five multiple aligners (MCP `gapModel=affine` honours `gapOpen`): exact natural affine sum-of-pairs objective (`SequenceAligner.SumOfPairsScore`), exact affine pairwise stage for the star methods, Gotoh three-state affine profile DP for the progressive / iterative methods (heuristic for k ≥ 3, equal to `GlobalAlignAffine` for k = 2), affine pair library for the consistency method. Oracles: independent Python (objective, exact optimum of 260 tiny sets, star/progressive/iterative/consistency rows). Tests: `SequenceAligner_MultipleAlignAffine_Tests` (16) + 2 MCP tests. Default stays linear.


## Audit r1 addendum (A1-12 / A1-16, B09 F53 / F54)

A1-12: the private UPGMA copy is gone (`Phylogenetics` never used `Alignment`; the reference was removed and `Alignment` now calls the canonical `PhylogeneticAnalyzer` builders through `TreeBuildOptions`); additive `MsaOptions.GuideTree` (`Upgma` default / `NeighborJoining` = ClustalW's tree, MCP `guideTree=upgma|nj`). The NJ tree and the progressive merge order equal the installed `clustalw` 2.1 (`.dnd` and "Group n" lines) on 420 DNA sets; the canonical NJ/UPGMA equal Biopython 1.88 NJ / scipy average linkage on 400 random matrices. A1-16: `ChainAnchors` is a range-maximum sparse DP (O(a log² a)) identical to the O(a²) scan on > 1 900 random, MEM-like and adversarial sets; the dead internal `MultipleAlignClassic` baseline is removed (the explicit benchmark fixture compares `MultipleAlign` with `MultipleAlignAnchored`).


## Audit r1 addendum (A1-15b / A1-10b, B09 F55 / F56)

`MsaDistanceMode.ClustalWQuick` (+ `MsaOptions.QuickPairwise`, public `SequenceAligner.ClustalWQuickPairwiseAlign`): exact port of ClustalW 2.1's fast k-tuple pairwise stage (`-QUICKTREE`), printed `Score:` equal on 8 257 pairs of the installed binary (k-tuple 1-4, window / topdiags / pairgap / absolute-score variants, IUB-N, partial-alignment cases); MCP `distance = linear | clustalw | clustalw-quick` with `ktuple / window / topDiagonals / pairGap`. New `SequenceAligner.MultipleAlignClustalW(sequences, ClustalWOptions?)` (+ MCP `multiple_align_clustalw`): step-by-step port of ClustalW 2.1's multiple-alignment stage for DNA (NJ guide tree with ClustalW's arithmetic and rooting, tree-derived sequence weights, tree similarities, delayed divergent sequences, position-specific gap penalties, shifted IUB matrix with `-TRANSWEIGHT`, Myers-Miller profile alignment); every row equal to `clustalw -TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` on 1 000 committed sets (+ 2 500 development sets) covering `-GAPOPEN -GAPEXT -TRANSWEIGHT -MAXDIV -NOWEIGHTS -NEGATIVE -QUICKTREE -KTUPLE -WINDOW -TOPDIAGS -PAIRGAP -PWGAPOPEN -PWGAPEXT`. Not ported (recorded as audit items): MUSCLE profile score (A1-10c), protein stage and protein-only penalties (A1-14).


## Audit r1 addendum (A1-13, B09 F58)

`SequenceAligner.MultipleAlignTCoffee` (+ `TCoffeeLibrary`, `TCoffeeTreeNode`, `MsaOptions.ConsistencyEngine = TCoffeeLinked`, MCP `consistencyEngine`): exact port of the T-Coffee 13.41 program's progressive stage (library relaxation + sparse linked pair-wise DP, gap open = gap extend = 0) for an explicit primary library and guide tree: 800 / 800 committed random DNA sets equal `t_coffee` (relax 1 and 0; libraries from this repository and from the program). The earlier explanation (F7: `-extend_mode`, extension rescale) was wrong - `-extend_mode` has no effect with this DP. Default behaviour unchanged (`PaperExtension`). Remaining (A1-13b, DOABLE): the program's own primary libraries (`slow_pair` Myers-Miller + `dna_idmat`, `lalign_id_pair` SIM top-10, default `proba_pair`) and its k-tuple NJ guide tree; end to end today 21 / 150 sets equal `t_coffee -method slow_pair,lalign_id_pair`.


## Audit r1 addendum (A1-11, B09 F60)

"Not implemented: exact optimal MSA" is lifted for small k: `MultipleAlignExact` (k = 2 = `GlobalAlign*`; k >= 3 exact k-dimensional DP with natural-affine gap-history states and the Carrillo-Lipman bound, complexity guard `MsaOptions.ExactMaxStateCells`, <= 8 sequences linear / 6 affine; MCP `multiple_align_exact`). Equal to two independent Python exact programs, exhaustive enumeration and Biopython (k = 2) on 804 / 804 sets; the heuristic aligners never exceed it and reach it on 86-92 % of the small related sets measured (algorithm doc §4.7).


## Audit r1 addendum (A1-14, B09 F61)

Protein / substitution-matrix multiple alignment: every aligner except `MultipleAlignClustalW` has `ProteinSequence` / `string` + `MatrixScoring` overloads (BLOSUM62, PAM250, EDNAFULL; linear and affine; fractional gap scores) that share all dynamic programs, merge, profile, refinement, library, lattice and anchor code with the DNA overloads through the exact-integer `PairScoring` (the scalar scoring is its scale-1 case; all pre-existing DNA tests unchanged); `MatrixSumOfPairsScore`, `MultipleAlignmentResult.ExactTotalScore`, `ExactMsaResult.ScoreScale`; MCP `type = dna | protein`, `matrix` on `multiple_align*`. Rows equal an independent Python implementation on 238 sets, the exact optimum equals the lattice DP / Biopython / enumeration on 308 sets, 320 rows of the `t_coffee` progressive stage on protein equal, see Evidence. Not covered: ClustalW's protein profile stage (A1-10c(b)).


## Audit r1 addendum (A1-15b part b, B09 F62)

ClustalW 2.1's PROTEIN pairwise stages are exact ports: `ClustalWPairwiseAlign(a, b, ClustalWProteinPairwiseOptions)` (matrix GONNET250 / BLOSUM30 / PAM350 / identity x intScale, length-dependent gap opening `2 * matrixAvg * (PWGAPOPEN + ln MIN(len))`, `PWGAPEXT * intScale` in single precision) and `ClustalWQuickPairwiseAlign(..., ClustalWSequenceType.Protein, ...)` (k-tuple 1-2 over 25 residue codes, defaults 1 / 5 / 5 / 3); `SubstitutionMatrix.ClustalWGonnet250 / ClustalWBlosum30 / ClustalWPam350 / ClustalWIdentity` (canonical references: EMBOSS EBLOSUM30, Biopython GONNET1992); `MsaOptions.ClustalWProtein` feeds `MsaDistanceMode.ClustalWAffine | ClustalWQuick` for protein matrices; MCP `pwMatrix`, `pwGapOpen`, `pwGapExtend`. Printed `Score:` of the installed `clustalw -TYPE=PROTEIN` equal on 4 204 full-stage and 3 774 quick-stage pairs, protein MSA rows equal an independent Python implementation on 60 sets x 3 methods. Not covered: ClustalW's protein profile stage / `MultipleAlignClustalW` for protein (A1-10c(b)).

## Audit r1 addendum (A1-10c part b, B09 F63)

ClustalW 2.1's PROTEIN multiple-alignment stage is a step-by-step port: `SequenceAligner.MultipleAlignClustalW(IEnumerable<ProteinSequence> | IEnumerable<string>, ClustalWProteinOptions, MatrixScoring?)` = the DNA stage (F56) with the protein branches of `profileAlign` / `getProfileAlignMatrix` / `calcGapCoeff`: the matrix series chosen by the percent identity of the two profiles (GONNET80 / 250 / 120 / 350 / 160, BLOSUM80 / 62 / 45 / 30, PAM20 / 60 / 120 / 350, ID; `-NEGATIVE` = one matrix), the matrix- and length-scaled gap opening (`scale * matrixAvg * (float)(GAPOPEN / (logdiff * logmin))`), the residue-specific (Pascarella & Argos), hydrophilic-run and variable-region (`-NOVGAP`) gap penalties and the gap-separation distance (`-GAPDIST`), the protein pairwise stages of F62, and the options `-MATRIX -PWMATRIX -GAPOPEN -GAPEXT -PWGAPOPEN -PWGAPEXT -MAXDIV -NOWEIGHTS -NEGATIVE -QUICKTREE -KTUPLE -WINDOW -TOPDIAGS -PAIRGAP -NOPGAP -NOHGAP -NOVGAP -GAPDIST -HGAPRESIDUES`. Additive: `ClustalWProteinOptions`, `ClustalWProfileMatrix`, `SubstitutionMatrix.ClustalWBlosum40 / 45 / 80 / Pam20 / 60 / 120 / Gonnet40 / 80 / 120 / 160 / 350` (the BLOSUM / PAM tables are EMBOSS data files cell for cell), MCP `multiple_align_clustalw type=protein` (+ `matrix pwMatrix noPreferredGaps noHydrophilicGaps variableRegionGaps gapDist hydrophilicResidues scoreMatrix`); the DNA behaviour and every existing test are unchanged.

Oracle: 1 000 committed protein sets (2-25 sequences, <= 365 residues, 20 option variants, B/Z/X/O/U, duplicates, tiny sets): every row of `clustalw -TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` identical (two `-NOVGAP` sets with the forced `gdist` below), two development panels (700 + 3 000 sets, other seeds) equal except 5 residual `-NOVGAP` sets. Residual and its source: with `-NOVGAP` and a two-sequence profile of > 60 % identity ClustalW reads the uninitialised member `ProfileBase::gdist` (undefined behaviour); the value is whatever the recycled heap chunk held (0, the user's `-GAPDIST`, or a stale residue code). The port emulates the usual value and `ClustalWProteinOptions.UninitialisedGapDistances` reproduces a specific run; 7 of 4 700 sets (0.15 %) are affected, all with `-NOVGAP`.

Not done (registered, A1-10d): user matrix files `-MATRIX=file` / `-PWMATRIX=file` (series files with limits), the interactive-menu end-gap options (`endGapPenalties`, `useEndGaps = true`; the command-line `-ENDGAPS` is a verified no-op), for which a Docker oracle package is prepared.


## Audit r1 addendum (A1-10d file part, B09 F64)

ClustalW's user-defined matrices are ported: `ClustalWUserMatrix.Parse` / `ClustalWMatrixFile.Parse` read the files of `-MATRIX=file` (one matrix or a `CLUSTAL_SERIES` file, at most 10 entries, scale `0.5 + (pcid - llimit) / ((ulimit - llimit)·2)`) and `-PWMATRIX=file` exactly as `SubMatrix::readUserMatrix` / `readMatrixSeries` do (NCBI / EMBOSS and unlabelled layouts, per-row x10 floating-point scale, lower triangle only, `letters + 1` rows, quirks of unknown header letters), through the same `getMatrix` port that serves the built-in tables (`ClustalWProteinOptions.UserMatrix / UserMatrixSeries`, `ClustalWProteinPairwiseOptions.UserMatrix`; `ClustalWUserMatrix.FromSubstitutionMatrix` accepts a canonical matrix); MCP `matrixFile` / `matrixFiles` / `pwMatrixFile`. Oracle: 240 committed sets (rows and 2 000+ pairwise `Score:` lines of the installed `clustalw` equal; 800-set development run: only the two documented `-NOVGAP` residuals differ).

Closed by the addenda below: A1-10d part 2 (end-gap options, F67) and part 3 (nucleotide matrices, F68).

## Audit r1 addendum (A1-10d part 2, B09 F67)

`ClustalWOptions` / `ClustalWProteinOptions` gained `UseEndGaps` (`UserParameters::useEndGaps`) and `EndGapPenalties` (`endGapPenalties`), the program's menu-only end-gap options (the command-line `-ENDGAPS` is a no-op); MCP `multiple_align_clustalw` `useEndGaps` / `endGapPenalties`. The program's silent abort (a negative first-pass profile score with end-gap penalties: no output file, exit status 0) is an explicit contract, `ClustalWAbortedException`. Oracle = the user-run Docker package `handoff/A1-10d` (patched `-ENDGAPS` handler; control = installed binary 340 / 340): 340 sets x 4 modes equal (306 + 306 aligned rows, 34 + 34 aborts, 340 + 340 control / useendgaps rows).

## Audit r1 addendum (A1-10d part 3, B09 F68)

ClustalW's nucleotide matrices are ported: `-DNAMATRIX` / `-PWDNAMATRIX` = IUB (default), CLUSTALW (`clustalvdnamt`, profile gap scale 0.66, pairwise 0.6667 / 0.751) or a user file (the F64 reader), `-TRANSWEIGHT` through the matrix' own Xref (profile) / the residue codes (pairwise); `ClustalWOptions.DnaMatrix / PairwiseDnaMatrix / UserDnaMatrix / PairwiseUserDnaMatrix / UserDnaMatrixUse`, `ClustalWDnaPairwiseOptions` + `SequenceAligner.ClustalWPairwiseAlign(.., options)`, MCP `dnaMatrix` / `pwDnaMatrix` / `dnaMatrixFile` / `pwDnaMatrixFile` / `userDnaMatrixUse`. Finding: `clustalw -DNAMATRIX=file -PWDNAMATRIX=file` reads and validates the file and then registers it under the PROTEIN slot, so the option has no effect on a DNA alignment (a defect of the program); the interactive menu applies the file - both behaviours are ported (`UserDnaMatrixUse.CommandLine` / `Applied`). Oracle: 360 sets (125 command-line runs, 235 interactive-menu runs): every row, 3 485 pairwise scores and 1 393 group-score lines equal; undefined behaviour with headers of fewer than 13 residues (writes outside the C++ matrix: 8 of 37 menu runs crashed) is documented, not modelled.

## Audit r1 addendum (A1-10c part a + a2, B09 F65 / F69 / F70)

`SequenceAligner.MultipleAlignMuscle` (+ `MuscleOptions`, MCP `multiple_align_muscle`) is a step-by-step port of MUSCLE **3.8.1551** (k-mer distance, UPGMB tree, ClustalW weights, log-expectation / SPN profile-profile DP with position-specific gap penalties, Kimura-distance tree refinement of changed subtrees; §4.9 of the algorithm doc). The installed `muscle` is 5.x (a different program), so the oracle is the MUSCLE v3.8.1551 binary built from the same source by the **Docker hand-off package `handoff/A1-10c-muscle`, run by the user: rows and output order equal on 420 / 420 sets (F69)**, no defect found. Audit part a2 (F70, F76) adds, implemented from the program's source and **verified against the binary on the second oracle round** (`handoff/A1-10c-muscle-r2`, run by the user: 630 recorded sets of 638, rows in the written order, trees / distance / weight files, `SP=` / `Score=` and the program's refusals **equal on 630 / 630 after six defects of the port found by the oracle were fixed**: the `kbit20_3` bit count, `-brenner` after the tree refinement, the leaf-index / sequence-id weight lookup of `-usetree`, negative zero in the files, `-maketree` ignoring `-distance2`, a one-sequence profile re-rooted; `SequenceAligner_MuscleOracleRound2_Tests`; 8 `-matrix blosum62nostar` sets could not be recorded by the round-2 harness and are re-run in `handoff/A1-10c-muscle-r3`, with `-matrix` for nucleotides and corner cases - awaiting the user's run): iterations >= 3 (`RefineHoriz` / `RefineVert` with anchors), the objective scores `-objscore sp/ps/dp/xp/spf/spm`, the weightings `none/henikoff/henikoffpb/gsc/clustalw/threeway`, `-termgaps`, the distances `kmer20_3/kmer20_4/kbit20_3/pwkimura/pwscoredist` and `pctid_kimura/pctid_log/edit/scoredist`, neighbour joining and `-root1/-root2`, `-brenner`, `-matrix`, `-usetree`, the files `-tree1/-tree2/-distmx1/-distmx2/-computeweights`, `-clusteronly`, and the command lines `-profile`, `-refine`, `-spscore`, `-ppscore`, `-maketree`, `-refinew`, `-sw` (MCP `muscle_msa_command`); the options the 3.8.1551 source compiles out (`-diags*`, `-hydro*`, `-gapopen2/-gapextend2/-gapambig`) are accepted with no effect. They are also tested against the published definitions and hand traces of the source (`SequenceAligner_MuscleOptions_Tests`, 41 tests).

## Audit r1 addendum (A1-13b, B09 F66): the T-Coffee 13.41 program's own libraries, guide tree and default method

`SequenceAligner.MultipleAlignTCoffeeProgram` (+ MCP `multiple_align_tcoffee`) reproduces `t_coffee` 13.41 end to end: its `proba_pair` default (ProbCons pair-HMM), `slow_pair` (Myers-Miller with `idmat`, GOP -10 / GEP -1, or BLOSUM62), `lalign_id_pair` (Huang-Miller SIM top 10, residues mapped to the sequence by a second Myers-Miller alignment), the way the program merges several methods (a forked merge that returns list positions and therefore loses and duplicates entries), the k-tuple neighbour-joining guide tree (slow / fast NJ, rooting shared with the ClustalW port) and the progressive stage of A1-13. Verified against the installed program: 820 / 820 primary libraries (pairs, weights, counts), 1 480 / 1 480 guide trees and row sets (DNA and protein; plain, `slow_pair,lalign_id_pair`, single methods, reverse order, three methods), 100 / 101-sequence sets (fast NJ). F58's description of the primary library was corrected (`dna_idmat` is +10 / 0; the "gop -11" fit was `g + h k`). Open: RNA `proba_pair` (A1-13c), other `-method`s / tree modes / `-seq_weight` (A1-13d). Details: Evidence file, algorithm doc §4.10.


## Audit r1 addendum (A1-13c / A1-13d, B09 F71-F75): RNA `proba_pair` and the command-line options of the T-Coffee 13.41 program

`MultipleAlignTCoffeeProgram` now reproduces, row for row against the installed `t_coffee` 13.41 on 4 715 committed sets: RNA input (separate pair-HMM and `get_tot_prob`), `-relax_lib N`, `-filter_lib`, `-gapopen`, `-gapext`, `-cosmetic_penalty`, `-nomatch`, `-do_normalise`, `-tree_mode upgma`, `-distance_matrix_mode ktup / fast / idscore / no` with the `-distance_matrix_sim_mode` measures, `-seq_weight t_coffee` / weight files, `-weight no`, `-matrix`, `method@EP@GOP|GEP|MATRIX@value`, `clustalw_pair`, and the name sort of `-inorder aligned` (the "name-prefix quirk") with `-outorder`. Options are properties of `TCoffeeProgramOptions` or a flag string (`FromCommandLine`, MCP `options`, `names`). Open (DOABLE, B09.md A1-13d): the fasta-heuristic methods and `-quicktree`, `clean_slow_pair`, `hash_pair`, `biphasic_pair`, `exon*`, `promo_pair`, `blastr_pair`, `lalign_rs_s_pair`, `subop*`, `-dp_mode` variants, the `slow` / `cscore` distance modes. Blocked with proof: BLAST / PDB / structure methods.

## Audit r2 addendum (A2-38 ... A2-42, B09 F81-F86): input contract, size guards, parallel exceptions, cancellation

One input contract for every multiple-alignment entry point: null element -> `ArgumentNullException` (it was `NullReferenceException` on the ClustalW / MUSCLE DNA overloads and the Seqeron aligners' `DnaRows`), malformed text of a user tree / T-Coffee value -> `ArgumentException`; the three program ports share one empty-member message (the programs refuse empty input), the Seqeron aligners keep the all-gap row. T-Coffee `proba_pair` refuses pairs above 2^28 pair-HMM floats and the MUSCLE profile DP above 1e9 cells with a clear `ArgumentException` (both programs have no limit of their own and wrap / abort). The star aligners rethrow the original exception of a parallel pair stage (never `AggregateException`) and take `MsaOptions.MaxDegreeOfParallelism`. Every long-running aligner is cancellable (ambient token polled at row / iteration boundaries; `MsaOptions` / `ClustalWOptions` / `ClustalWProteinOptions` / `MuscleOptions` / `TCoffeeProgramOptions.CancellationToken`, `SequenceAligner.WithCancellation`; MCP tools take the request token). Tests: `SequenceAligner_InputContract_Tests`, `SequenceAligner_Cancellation_Tests`; TestSpec / Evidence sections of this unit.
