# Validation Report: COMPGEN-ANI-001 — Average Nucleotide Identity (ANI)

- **Validated:** 2026-06-24   **Area:** Comparative Genomics
- **Canonical method(s):**
  - `ComparativeGenomics.CalculateANI(query, reference, fragmentLength, minIdentity, minAlignableFraction, gapped)`
  - `ComparativeGenomics.CalculateReciprocalAni(genomeA, genomeB, fragmentLength, minIdentity, minAlignableFraction, gapped)`
  - private `BestUngappedFragmentMatch`, `BestGappedFragmentMatch` (Smith-Waterman via `SequenceAligner.LocalAlign`, `BlastDna` scoring)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** CLEAN

This re-validation covers the limitations-campaign change in commit `69c51fa0`: gapped fragment
alignment (reusing `SequenceAligner`) and reciprocal (two-way) ANI. The prior report (2026-06-16)
graded Stage B PASS-WITH-NOTES because the implementation was ungapped/single-direction and the
`minAlignableFraction` parameter was structurally inert. Both notes are now resolved.

## Stage A — Description

### Sources opened & what they confirm (retrieved this session, 2026-06-24)

1. **Goris et al. 2007** (IJSEM 57:81-91, DOI:10.1099/ijs.0.64483-0) — method text confirmed verbatim
   via WebSearch (reproduced in pyani docs / secondary literature):
   - Query genome cut into **consecutive 1020 nt fragments** (to mirror ~1 kb DDH fragmentation).
   - Each fragment searched against the reference whole genome by **BLASTN**; best match kept.
   - **ANI = "mean identity of all BLASTN matches that showed more than 30 % overall sequence
     identity (recalculated to an identity along the entire sequence) over an alignable region of
     at least 70 % of their length."**
   - **"Reverse searching, in which the reference genome is used as the query, was also performed
     to provide reciprocal values."** → the symmetric ANIb value averages the two directions.
   - BLASTN is **gapped** (default), so the best match recovers indels.
   - "ANI values of approximately 95 % correspond to the 70 % DNA-DNA hybridization standard."
2. **Konstantinidis & Tiedje 2005** (PNAS 102:2567-2572) — the ≈94 % ANI / 70 % DDH correspondence;
   the modern ≈95–96 % species boundary. Confirmed (prior session + this session search context).
3. **pyani ANIb source module** (`pyani/anib.html`, WebFetch this session) — the reference
   implementation's exact per-fragment conventions:
   - `ani_alnlen = blast_alnlen − blast_gaps` (ungapped aligned columns),
   - `ani_alnids = ani_alnlen − blast_mismatch` (identical ungapped columns),
   - `ani_pid = ani_alnids / qlen` → **identity recalculated over the full query-fragment length**,
   - `ani_coverage = ani_alnlen / qlen` → coverage = ungapped aligned length / fragment length,
   - filter: `(ani_coverage > 0.7) & (ani_pid > 0.3)`, gapped BLASTN with `-xdrop_gap_final 150`,
     `-max_target_seqs 1`.

### Formula check
Documented model: per fragment, `id = identical(ungapped) columns / L`, `cov = ungapped aligned
columns / L`; qualify iff `id > 0.30 ∧ cov ≥ 0.70`; ANI = mean of qualifying `id`. Reciprocal ANI =
`(ANI(A→B) + ANI(B→A)) / 2`. Both match Goris 2007 and pyani exactly (identity over full fragment
length; strict `>30 %`; `≥70 %` coverage; reverse-search averaging). Defaults `L=1020`, `0.30`,
`0.70` trace to Goris 2007. ✔

### Edge-case semantics
- Trailing partial fragment dropped (consecutive non-overlapping). ✔ (Goris "consecutive 1020 nt").
- Non-conserved fragments (below identity or coverage cut-off) discarded, not counted as 0. ✔
- Empty / no-qualifying-fragment → 0 (mean over empty set reported as 0). ✔
- Direction-dependence acknowledged; symmetry restored only by reciprocal averaging. ✔
- Gapped placement: identity/coverage both over fragment length; gaps excluded from both. ✔ (pyani).

### Independent cross-check (numbers)
Hand-derived the changed-path values:
- **G2 indel recovery.** Fragment `AAAACCCC` (8 nt) vs ref `AAAATCCCC` (one inserted T).
  - Ungapped: any full-length 8-mer offset loses column registration past the insert → best 7/8 =
    **0.875**.
  - Gapped (`BlastDna`: match +2, mismatch −3, gapOpen −5, gapExtend −2): alignment `AAAA-CCCC` vs
    `AAAATCCCC` scores 8·(+2) + 1 gap(−5−2) = 16 − 7 = **+9**, beating any ungapped local (max 7
    matched columns → ≤ 14, and a broken alignment scores less). All 8 query columns identical →
    `identicalColumns/8 = 1.0`, `ungappedColumns/8 = 1.0` ≥ 0.70 → qualifies → **gapped ANI = 1.0 >
    ungapped 0.875**. Confirms INV-6.
- **G3 coverage cut-off.** `AAAACCCC` (8) vs `AA` (2): best local covers ≤ 2 query columns →
  `cov = 2/8 = 0.25 < 0.70` → excluded → **0.0**. ✔
- **G4 identity cut-off (gapped).** `CGTC` vs all-A → 0 identical columns → `id = 0.0`, not `> 0.30`
  → excluded; only `AAAA` (1.0) kept → **1.0**. ✔
- **R3 reciprocal = mean of directions.** A=`AAAACGTC`, B=`AAAAAAAA`, fragLen 4: A→B = 1.0 (CGTC
  excluded), B→A = 1.0 (both `AAAA` match A's lead) → **(1.0+1.0)/2 = 1.0**. ✔
- **R1/R2.** Identical → mean(1,1)=1.0; symmetry holds for any pair since the mean of the two
  directions is order-independent. ✔ (INV-7, INV-8)

### Findings / divergences (Stage A)
None. The description now matches the full ANIb definition (gapped BLASTN + reciprocal averaging)
per Goris 2007 / pyani. The remaining out-of-scope item (OrthoANI/FastANI orthologous-best-pair
variant) is honestly disclosed as not implemented. → **Stage A PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:1076-1228`.

- **Validation** (`:1084-1087`): null/empty → 0; `fragmentLength ≤ 0` → `ArgumentOutOfRangeException`. ✔
- **Fragmentation** (`:1094`): `for (start=0; start+fragmentLength <= query.Length; start += fragmentLength)`
  — consecutive, non-overlapping, trailing partial dropped. ✔ (INV-4)
- **Per-fragment placement** (`:1096-1098`): `gapped` switch selects `BestGappedFragmentMatch` (full
  Smith-Waterman) or `BestUngappedFragmentMatch` (best-offset scan). ✔
- **Filter** (`:1101`): `identity > minIdentity && alignableFraction >= minAlignableFraction`
  — strict `>` on identity (Goris "more than 30 %"), `≥` on coverage; mean of kept identities, else 0. ✔
- **`BestGappedFragmentMatch`** (`:1194-1225`): aligns fragment vs reference with `BlastDna`; counts
  `identicalColumns` (both bases equal, gaps skipped) and `ungappedColumns` (no gap either side);
  `identity = identicalColumns/fragLen`, `alignableFraction = ungappedColumns/fragLen` — exactly
  pyani `ani_pid = ani_alnids/qlen` and `ani_coverage = (blast_alnlen − blast_gaps)/qlen`. ✔
- **`CalculateReciprocalAni`** (`:1126-1144`): same validation; returns `(forward + reverse)/2.0`
  with `forward = CalculateANI(A,B,…)`, `reverse = CalculateANI(B,A,…)`. ✔ (INV-7, INV-8)
- **`BestUngappedFragmentMatch`** (`:1155-1182`): unchanged best-offset matching-base count;
  identity = bestMatches/fragLen, coverage 1.0 (ref ≥ fragLen) else (0,0). ✔

### Formula realised correctly?
Yes. The gapped path computes the pyani-defined identity and coverage over the fragment length; the
reciprocal method computes the Goris reverse-search mean. The `minAlignableFraction` parameter is now
**genuinely active** under the gapped path (coverage can be any value in [0,1]), resolving the prior
report's "structurally inert" note (G3 exercises coverage 0.25 < 0.70). ✔

### Cross-verification table recomputed vs code
| Case | Inputs | Hand value | Test result |
|------|--------|-----------|-------------|
| G1 gapped identical | R,R fragLen 4 | 1.0 | 1.0 ✔ |
| G2 gapped indel | `AAAACCCC`/`AAAATCCCC` fragLen 8 | gapped 1.0 > ungapped 0.875 | 1.0 / 0.875 ✔ |
| G3 gapped coverage cut-off | `AAAACCCC`/`AA` | 0.0 | 0.0 ✔ |
| G4 gapped identity cut-off | `AAAACGTC`/`AAAAAAAA` fragLen 4 | 1.0 | 1.0 ✔ |
| R1 reciprocal identical | R,R | 1.0 | 1.0 ✔ |
| R2 reciprocal symmetry | `…TTTT`/`…TTTA` | ANI(A,B)=ANI(B,A) | equal ✔ |
| R3 reciprocal = mean | `AAAACGTC`/`AAAAAAAA` | (1.0+1.0)/2 = 1.0 | 1.0 ✔ |
| M1–M8,S1,S2,C1 (ungapped) | per prior report | unchanged | all ✔ |

All match the independently hand-derived values.

### Variant/delegate consistency
Two public methods; `CalculateReciprocalAni` delegates to `CalculateANI` per direction (no divergent
logic). The MCP wrapper (`AnalysisTools.cs`) forwards to `CalculateANI`. No `*Fast` variants. ✔

### Numerical robustness
Identity/coverage are bounded by `Math.Min(…, 1.0)`; division denominators (`fragLen`) are positive
(validated). Empty alignment / ref-shorter cases return `(0,0)`. No overflow/precision concern on
stated ranges. ✔

### Test quality audit (HARD gate)
File: `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CalculateANI_Tests.cs` (20 tests).
- **Sourced, not code-echoed:** every deterministic value is hand-derived from Goris/pyani; G2
  contrasts gapped 1.0 vs ungapped 0.875 (a real indel-recovery discriminator); R3 ties reciprocal
  to the mean of two directions; M4 would expose an average-all-fragments bug (would return 0.5). ✔
- **No green-washing:** exact `Within(1e-10)` equality on exact fractions; the only `Is.InRange` (S1)
  is the range invariant; `Is.GreaterThan` (G2) is an additive correctness assertion atop exact
  equalities, not a weakened oracle. No skips/ignores/widened tolerances. ✔
- **Branch coverage:** gapped vs ungapped, identity cut-off (M4/G4/C1), coverage cut-off (M5/G3),
  fragmentation/trailing (M6), null/empty (M7/R4), invalid fragmentLength (M8/R5), reciprocal
  symmetry & mean (R1/R2/R3), range invariant incl. default-1020 run (S1). ✔
- **Honest green:** filtered class 20/20 pass; broader `~Comparative` filter (fuzz/metamorphic/
  algebraic/combinatorial/property ANI families) 480/480 pass; build 0 Error(s)/0 Warning(s). ✔

**Gate result: PASS.**

### Findings / defects (Stage B)
None. The two prior PASS-WITH-NOTES items are resolved:
1. `minAlignableFraction` is now active (gapped coverage varies; G3 exercises exclusion).
2. The ungapped/single-direction simplification is replaced by gapped placement + reciprocal
   averaging matching the full ANIb definition.

## Verdict & follow-ups
- **Stage A: PASS** — description matches Goris 2007 (1020 nt fragments, >30 %/≥70 % cut-offs, gapped
  BLASTN, reverse-search reciprocal) and pyani identity/coverage conventions verbatim; species
  boundary confirmed.
- **Stage B: PASS** — gapped and reciprocal code faithfully realise the validated formulas; all
  cross-check values reproduce; tests are exact, sourced, discriminating, and branch-complete.
- **End-state: CLEAN** — no defect found; no code or test changes required. OrthoANI/FastANI
  orthologous-best-pair variant remains an honestly disclosed out-of-scope extension, not a defect.
