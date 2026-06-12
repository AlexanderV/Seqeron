# Validation Report: CHROM-SYNT-001 — Synteny Analysis (collinear blocks & rearrangement detection)

- **Validated:** 2026-06-12   **Area:** Chromosome Analysis
- **Canonical method(s):** `ChromosomeAnalyzer.FindSyntenyBlocks(orthologPairs, minGenes=3, maxGap=10)`, `ChromosomeAnalyzer.DetectRearrangements(syntenyBlocks)`
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs` (lines 95–116 records, 642–822 algorithms)
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_Synteny_Tests.cs` (20 tests)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES

---

## Stage A — Description

### Sources opened & what they confirm

- **Wikipedia "Synteny"** (https://en.wikipedia.org/wiki/Synteny): modern synteny = *collinearity*, "conservation of blocks of order within two sets of chromosomes." Syntenic blocks = regions preserving the precise order of homologous genes from a common ancestor. (Note: the article does **not** explicitly discuss strand/inverted blocks — see divergence below.)
- **Wikipedia "Chromosomal rearrangement"**: confirms the four structural types used by `DetectRearrangements`:
  - **Inversion** — a segment reversed/flipped within a chromosome.
  - **Translocation** — segment transferred between non-homologous chromosomes.
  - **Deletion** — segment removed (loss of material).
  - **Duplication** — segment copied (extra material).
- **Wang et al. (2012) MCScanX** (Nucleic Acids Res. 40(7):e49): synteny detection = dynamic-programming chaining of *pair-wise collinear anchor genes*; rewards adjacent anchor pairs, penalises distance. Key parameters: `-s MATCH_SIZE` (min genes to call synteny, default 5), `-m MAX_GAPS` (max gaps allowed, default 20), `-u UNIT_DIST`. This corroborates the spec's two parameters: **minimum block size (`minGenes`)** and **gap tolerance (`maxGap`)**. Inverted (reverse-collinear, anti-diagonal) blocks are standard MCScanX/SyRI output and represent inversions.
- **Goel et al. (2019) SyRI** (Genome Biology 20:277): classifies syntenic vs rearranged regions (inversions, translocations, duplications) between assemblies — corroborates forward vs inverted distinction and the rearrangement taxonomy.

### Block definition (validated)
A **synteny block** is a maximal run of homologous marker pairs that are **collinear** — i.e. consecutive markers preserve a consistent ordering on both genomes. A run is **forward ('+')** when target positions increase with reference positions, and **inverted ('-')** when target positions decrease. Parameters: **min anchors** (`minGenes`) and **allowed gap** (`maxGap`). A gap larger than the threshold, or an order reversal, breaks the run.

### Worked examples (hand-computed)
Reference order genome1 = 1,2,3,4; let "genome2" be the target order.

| genome2 | Expected | Reasoning |
|---------|----------|-----------|
| 1,2,3,4 | **1 forward block of 4 ('+')** | strictly increasing target → collinear forward |
| 4,3,2,1 | **1 inverted block of 4 ('-')** | strictly decreasing target → reverse-collinear |
| 1,2,5,3,4 (markers 1,2,5 increasing, then 3,4 after a reversal) | **1 forward block of 3** {1,2,5}; trailing {3,4} is a 2-anchor run dropped below minGenes=3 | the 5→3 step reverses direction (target decreases) → block break |

These match the spec's M1/M2 and the gap/break semantics in §4 of the TestSpec and §2.2 of the Evidence doc.

### Edge-case semantics (all defined & sourced)
empty → empty; below minGenes → empty; fully collinear → 1 forward block; fully inverted → 1 inverted block; gap > maxGap → split; multiple chromosome pairs → separate blocks. `SequenceIdentity = NaN` because identity is not computable from coordinate-only input (requires BLAST/alignment per MCScanX) — sound.

### Stage A findings / divergences (notes)
1. **Strand/inverted blocks not in the Wikipedia "Synteny" article.** Forward/inverted ('+'/'-') classification is nonetheless standard and authoritative via MCScanX dot-plots and SyRI inversion calls. Definition is correct; the citation label "Wikipedia (Synteny)" for the inverted-block test (M2) is weakly sourced — MCScanX/SyRI is the proper authority. Minor, documented.
2. The spec's collinearity rule is the simpler "consecutive consistent ordering," not MCScanX's full DP/LIS chaining. Both are valid block definitions; the spec's is a legitimate simplification.

**Stage A verdict: PASS-WITH-NOTES** — biology/definitions are correct; the inverted-block concept is correctly attributable to MCScanX/SyRI rather than the (silent-on-strand) Wikipedia Synteny page.

---

## Stage B — Implementation

### Code path reviewed
- `FindSyntenyBlocks` (lines 642–715): groups pairs by (Chr1,Chr2); sorts by `Start1`; walks consecutive pairs; `currentForward = curr.Start2 > prev.End2`; first step fixes block direction `isForward`; a pair stays collinear iff `currentForward == isForward` **and** `gap1 ≤ maxGap·1e6` **and** `gap2 ≤ maxGap·1e6`; on break or last index, emits a block when `geneCount ≥ minGenes`; block strand = `isForward ? '+' : '-'`; Species2 coords use Min/Max of first/last to stay valid under inversion; `SequenceIdentity = NaN`.
- `DetectRearrangements` (lines 720–822): sorts blocks by (Chr1, Start1); for adjacent same-Chr1 pairs: different Chr2 → **Translocation**; same Chr2 + different strand → **Inversion**; same Chr2 + same strand + asymmetric gap (`gap1 > 0 && gap2 ≥ 0 && gap1 > 2·gap2`) → **Deletion** (Size = gap1 − max(0,gap2)); separate O(n²) pass: overlapping Species1 regions mapping to different Species2 targets → **Duplication**.

### Formula realised correctly? (traced vs code)
- **M1 forward 4-gene** → 1 block, '+', GeneCount 4, Species1 1000–8000, Species2 1000–8000. ✔ (traced)
- **M2 reverse 4-gene** → 1 block, '-', GeneCount 4, Species2 Min/Max = 2000/9000. ✔ (traced)
- **M5 exactly minGenes (3)** → 1 block, GeneCount 3, '+', 1000–6000. ✔
- **M6 two chr pairs** → 2 blocks, 3 genes each. ✔
- **M16 gap split** (maxGap=2 → 2 MB; genes 3→4 gap = 2,994,000 bp > 2 MB) → 2 blocks of 3 (1000–6000 and 3,000,000–8,000,000). ✔ (traced step-by-step)
- **M9 Inversion** → Pos1=50000, Pos2=60000, Size=10000. ✔
- **M10 Translocation** → Chr2="chrB", Pos1=50000, Pos2=1000, Size=null. ✔
- **M14 Deletion** → gap1=100000, gap2=5000, 100000 > 2·5000 → Size = 100000−5000 = 95000. ✔
- **M15 Duplication** → overlap 20000–50000, different Species2 → Size = 30000. ✔
- **Interleaved 1,2,5,3,4** → 1 forward block of 3 {1,2,5}; trailing {3,4} dropped (below minGenes). ✔ (traced)
- Edge: empty→empty; single marker→empty; break exactly at last index does not lose a valid block (traced: 4-forward + 1-breaker → 1 block of 4, breaker dropped as singleton). ✔

### Variant/delegate consistency
No `*Fast`/delegate variants; two independent canonical methods. Consistent.

### Numerical robustness
`maxGap·1e6` computed as `int` — fits within `int.MaxValue` (~2.1e9) for `maxGap` up to ~2147; default 10 is safe. `gap2` uses `Math.Abs` to guard inverted runs. No div-by-zero.

### Test quality audit
20 tests; assertions are exact hand-calculated values (counts, coords, strand, sizes), deterministic, covering forward/inverted/gap-split/min-genes/multi-chromosome and all four rearrangement types plus invariants (strand ∈ {+,−}, Start ≤ End, NaN identity, valid type strings, Position1 set). Strong.

### Stage B findings / defects (notes, not blocking)
1. **Greedy pairwise chaining, not DP/LIS** (as in Stage A note 2): markers that remain *locally* monotonic but globally interleaved (e.g. m5 before m3,m4) are absorbed into a block instead of being broken out. This is faithful to the spec's "consecutive consistent ordering" definition but diverges from MCScanX's scoring DP. No test asserts the contrary; not a defect against the validated spec.
2. **`DetectRearrangements` adjacency-only / heuristic thresholds** (e.g. Deletion uses `gap1 > 2·gap2`): these are reasonable coordinate-only heuristics, sourced to the structural definitions, and exactly match the asserted test values. They are not the full SyRI classification but are internally consistent and tested.

**Stage B verdict: PASS-WITH-NOTES** — code faithfully realises the validated block/rearrangement definitions; every worked example and cross-check value reproduced exactly. Divergences are documented simplifications, not defects.

---

## Verdict & follow-ups

- **State: CLEAN** — no defect found; implementation matches the Stage-A-validated description; no code changes required.
- Tests: `--filter FullyQualifiedName~Synteny` → 20 passed, 0 failed. Full suite → **4484 passed, 0 failed** (baseline preserved).
- Follow-up (optional, non-blocking): the spec/Evidence could re-attribute the inverted-block (M2) and gap-split (M16) tests to MCScanX/SyRI rather than the Wikipedia "Synteny" page, which is silent on strand. Documentation-only.
