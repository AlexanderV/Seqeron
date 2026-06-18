# Validation Report: PHYLO-BOOT-001 — Phylogenetic Bootstrap Analysis (Felsenstein's bootstrap proportions)

- **Validated:** 2026-06-15   **Area:** Phylogenetics
- **Canonical method(s):** `PhylogeneticAnalyzer.Bootstrap(sequences, replicates, distanceMethod, treeMethod, seed)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (no algorithm defect; test-coverage gaps fixed in-session)

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **Felsenstein, J. (1985) — Confidence Limits on Phylogenies** (Evolution 39(4):783–791),
   via OSTI bibliographic record `https://www.osti.gov/biblio/6044842` (Wiley DOI page is paywalled).
   Confirmed verbatim:
   - resampling is **"resampling points from one's own data, with replacement, to create … bootstrap
     samples of the same size as the original data"** → each replicate has the same number of sites.
   - **"keep all of the original species while sampling characters with replacement"** → taxa (rows)
     are kept; characters/columns (sites) are resampled with replacement.
   - **"If a group shows up 95% of the time or more, the evidence for it is taken to be statistically
     significant"** → thresholds are interpretive, not part of the computation.
   - support P for a group = fraction of bootstrap samples in which the group appears; 100% → P=1.

2. **Lemoine et al. (2018) — Renewing Felsenstein's Phylogenetic Bootstrap** (Nature 556:452–456),
   PMC `https://pmc.ncbi.nlm.nih.gov/articles/PMC6030568/`. Restates the formal FBP procedure:
   "resample, with replacement, the sites of the alignment to obtain pseudo-alignments of the same
   length"; "measure the support of every branch in the **reference tree** as the proportion of
   pseudo-trees containing that branch"; per-replicate scoring is **binary** (a branch matches or
   not). Confirms the reference tree fixes the entity set.

3. **Biopython `Bio.Phylo.Consensus`** master source, fetched this session from
   `https://raw.githubusercontent.com/biopython/biopython/master/Bio/Phylo/Consensus.py`.
   Confirmed verbatim code:
   - `bootstrap` / `bootstrap_trees`: `for j in range(length): col = random.randint(0, length - 1)`
     then append column `msa[:, col:col+1]` → exactly `length` columns resampled with replacement,
     same-length pseudo-alignment.
   - `get_support`: scores only `target_tree.find_clades(terminal=False)` (non-terminal clades);
     a match is by terminal/leaf bitstring (`_clade_to_bitstr`); confidence `(t + 1) * 100.0 / size`.
     The `(t+1)` is the in-loop increment pattern: after a clade matches in `m` total trees the final
     stored confidence is `m * 100 / size` (= matches/size as a percentage). So the repo's
     `count/replicates ∈ [0,1]` equals Biopython's final value ÷ 100 — consistent.

### Formula check
Description §2.2 `support(c) = #{ r : c ∈ C(Tᵣ) } / B` over the non-trivial clades `C(T₀)` of the
reference tree matches Lemoine "proportion of pseudo-trees containing that branch" and Biopython's
final confidence exactly (units: proportion vs percentage — documented Assumption #2). ✔

### Edge-case semantics
- `null` → `ArgumentNullException`; `<2` taxa → `ArgumentException`; `replicates<1` → `ArgumentException`;
  unequal lengths → `ArgumentException` (from `BuildTree`). All are defined and sourced (a tree needs
  ≥2 taxa; a denominator needs ≥1; bootstrap requires an *alignment*).
- all-identical sequences → every reported clade support 1.0 (zero-distance matrix is invariant under
  any column resample → one topology every replicate). Sourced to INV-5 / Felsenstein 100%→P=1.
- two well-separated invariant groups → support 1.0 (distances invariant under resampling). ✔

### Independent cross-check (numbers)
- M1 (and NJ M7): with A=B=`A…`, C=D=`G…`, every column is identical w.r.t. the topology, so the
  distance matrix is identical for every resample → UPGMA **and** NJ recover {A,B},{C,D} on every
  replicate → support exactly **1.0**. Confirmed by an out-of-tree probe program against the built
  DLL: `A|B=1, C|D=1` for UPGMA (seed 42) and for NJ (seeds 1, 42), and `M6 A|B=1`.
- Resampling-axis sanity: a 4-column probe (`A=AAAA, B=AAAG, C=GGGA, D=GGGG`, p-distance, NJ, 200 reps)
  produced *fractional* seed-dependent support (`A|B` ≈ 0.925–0.975 across seeds 1/7/42/123), i.e. the
  topology genuinely varies across replicates as column resampling predicts — not a degenerate echo.

### Findings / divergences (Stage A)
- Two documented modeling choices (Assumptions, accepted): rooted-clade scoring (matches Biopython
  `get_support`) rather than unrooted bipartitions; support as a proportion in [0,1] (×100 = published
  percentage). Both are units/representation, not correctness, choices. **PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs`
- `Bootstrap` lines 894–956: validates (null/`<2`/`<1`); builds reference tree once and collects
  `GetClades(refTree.Root)`; per replicate draws `alignmentLength` column indices via one seeded
  `Random(seed)` and assembles a same-length pseudo-alignment over **all taxa**; rebuilds the tree;
  increments each reference clade present in the replicate; returns `count/replicates`.
- `GetClades`/`CollectClades` lines 764–794: clade key = sorted, `|`-joined leaf names; non-trivial =
  `>1 and < totalLeaves` leaves → exactly Biopython's `find_clades(terminal=False)` semantics.
- `CalculatePairwiseDistance` / `BuildUPGMA` / `BuildNeighborJoining` reused unchanged.

### Formula realised correctly?
Yes. Same-length column resampling with replacement (lines 928–938) matches Felsenstein/Biopython;
support = matching-replicate count ÷ replicates over the reference tree's non-trivial clades (lines
945–955) matches Lemoine/Biopython. Single seeded RNG → deterministic per seed (INV-4).

### Cross-verification recomputed vs code (out-of-tree probe + tests)
| Case | Source expectation | Code output |
|------|--------------------|-------------|
| Two-group UPGMA, seed 42 | {A,B}=1.0, {C,D}=1.0 (invariant distances) | `A|B=1, C|D=1` ✔ |
| Two-group NJ, seeds 1 & 42 | same (procedure tree-agnostic) | `A|B=1, C|D=1` ✔ |
| All-identical 3 taxa | every clade 1.0 | `A|B=1` ✔ |
| null / single / reps=0 / reps=−5 / unequal len | `ArgumentNullException` / `ArgumentException`×4 | matches ✔ |
| reps=1 (boundary) | valid proportions (here 1.0) | `A|B=1, C|D=1` ✔ |

### Variant/delegate consistency
`treeMethod` has two branches; both (`UPGMA`, `NeighborJoining`) now produce the sourced 1.0 result on
the invariant dataset. `distanceMethod` is orthogonal to the bootstrap loop (default JukesCantor used by
most tests; PDistance probed). No `*Fast`/instance variants for this method.

### Test quality audit (gate)
- **Pre-existing fixture (10 tests):** M1 exact 1.0 (sourced, not a code echo — derivation is
  RNG-independent), M2 [0,1] (INV-1), M3 quantization k/replicates (INV-2), M4 keys=reference clades
  (INV-3), M5 determinism (INV-4), M6 all-identical→1.0 (INV-5), S1/S2/S3 validation, S4 multi-seed
  validity. All assertions carry messages; no `ignore`/skip; no widened tolerances. The TestSpec §5.4
  claimed 11 but the file actually contained 10 — corrected in the spec.
- **Coverage gaps found (no algorithm defect):** (a) the `NeighborJoining` branch of `treeMethod` was
  never exercised; (b) the *documented* Stage-A edge case "unequal-length sequences throw" (doc §6.1 /
  contract) had no test; (c) the negative-replicate boundary was untested.
- **Fixed in-session:** added **M7** (NeighborJoining two-group → 1.0, source-anchored), **S5**
  (unequal-length → `ArgumentException`), **S3b** (negative replicates → `ArgumentException`). Fixture
  10→13. No assertion weakened, no value adjusted to match output, no skips.
- **Honest green:** full unfiltered suite **6561 passed / 0 failed / 0 skipped\*** (\*one pre-existing
  unrelated MFE benchmark is intentionally skipped); `dotnet build` 0 errors; the changed test file
  builds warning-free (the 4 build warnings are pre-existing NUnit-analyzer notes in unrelated
  `ApproximateMatcher`/other test files, untouched here).

### Findings / defects
No algorithm defect. Test-coverage gap only → classified FIXED-NOW (added M7/S5/S3b; corrected TestSpec
count + new rows).

## Verdict & follow-ups
- **Stage A: PASS** — description, formula, edge semantics and invariants all match Felsenstein (1985),
  Lemoine (2018) and Biopython `Consensus.py` retrieved this session; documented units/rooting choices
  are sourced and accepted.
- **Stage B: PASS-WITH-NOTES** — implementation faithfully realises the validated procedure; the only
  findings were untested branches/edge cases, fixed in-session with sourced assertions.
- **End-state: CLEAN** — no defect; coverage gaps closed; build + full suite green.
- **Test-quality gate: PASS** (after the M7/S5/S3b additions).

## Update (2026-06-17) — N-ary (multifurcating) tree model (C2)

After the `PhyloNode` N-ary refactor (C2), bootstrap clade collection traverses `Children` and is
correct over multifurcations. Neighbor-Joining now yields its natural **unrooted trifurcation** at the
centre (Saitou & Nei 1987), so a 4-taxon two-group NJ reference tree is `((A,B),C,D)`: its single
non-trivial rooted clade is `{A,B}` (support 1.0 in every replicate), while `{C,D}` crosses the
unrooted trifurcation and is no longer a rooted clade. The M7 NJ bootstrap test was corrected to this
sourced reality (previously it expected both `{A,B}` and `{C,D}` under the old binary root). No
assertion was weakened; the corrected expectation is hand-derived from the NJ trifurcation. See
FINDINGS_REGISTER §C2.

## Independent re-validation (2026-06-17, Phase-3 row #12) — N-ary clade collection

Fresh-context independent re-validation triggered by the N-ary (multifurcating) refactor (commit
c4f0190): `PhyloNode` now holds a `Children` list; bootstrap clade collection (`CollectClades`) and
support counting were rewritten, and the NJ-bootstrap test was updated for the trifurcating root.
The bootstrap *algorithm* is unchanged. Per the governing rule, expected values were hand-derived
from the bootstrap definition, not read from code output.

### Stage A — re-confirmed PASS
The bootstrap definition matches the three in-repo first sources (Felsenstein 1985 OSTI 6044842;
Lemoine 2018 PMC6030568; Biopython `Bio.Phylo.Consensus`): resample columns with replacement to the
**same length**, keep **all taxa**, support = **proportion in [0,1]** of replicate trees containing
each non-trivial reference-tree clade matched by leaf-name set; 100% agreement → 1.0. Pinned
invariants: support ∈ [0,1]; quantized to k/replicates; **1.0 when every replicate agrees**; resample
preserves length + taxa; deterministic per seed.

### Stage B — re-confirmed PASS; one test-coverage gap hardened
`CollectClades` (`PhylogeneticAnalyzer.cs:1011`) `foreach`-iterates **every** `node.Children` when
gathering a node's subtree leaf set, so all children — including the 3rd+ of a multifurcation — enter
the clade sets for both the reference and every replicate tree. Confirmed by a standalone probe against
the compiled assembly.

**N-ary clade-collection check (hand-derived):**
- 4-taxon two-group → NJ root `((A,B),C,D)` is a genuine 3-child trifurcation (probe:
  `root.Children.Count == 3`); the only non-trivial rooted clade is `{A,B}` → support 1.0.
- 6-taxon three-pair `((A,B),(C,D),(E,F))` → 3rd root child is the **internal** node `(E,F)`; clades =
  `{A|B, C|D, E|F}`, each support 1.0 (invariant distances → one topology per replicate).

**HARD mutation gate (run this session via standalone probe):** mutating collection to
`Children.Take(2)` drops the 3rd child — on the 6-taxon tree it **loses `{E,F}`** and emits a spurious
`{A,B,C,D}` (`{A|B, A|B|C|D, C|D}`); on the 4-taxon tree it emits a spurious `{A,B,C}`
(`{A|B, A|B|C}`). The pre-existing **M7** (4-taxon NJ) detects the drop only via an *extra* key, never
via a *lost* real clade.

**Fix (test-only hardening, no algorithm defect):** added **M8**
`Bootstrap_NeighborJoining_TrifurcationWithThreeInternalChildren_AllCladesScored` — the 6-taxon
trifurcation case, asserting keys ≡ `{A|B, C|D, E|F}`, each support exactly 1.0, and **no** spurious
`A|B|C|D` union. M8 is a lost-clade guard (it fails under the dropped-3rd-child mutation), complementing
M7's extra-key guard. All expected values hand-derived; no assertion weakened, no skip. Bootstrap
fixture 13→14 tests.

### End-state: CLEAN
No code or description defect. Test-quality gate PASS (M1 1.0-on-identical and resample-length/taxa
invariants exact; determinism via fixed seed; M8 fails if a 3rd child is dropped). Full unfiltered
suite **6780 passed / 0 failed**; build 0 errors. See FINDINGS_REGISTER §A57 and ledger Phase-3 #12.
