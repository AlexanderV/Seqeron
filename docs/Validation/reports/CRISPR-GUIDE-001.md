# Validation Report: CRISPR-GUIDE-001 — On-target guide efficacy (Doench 2014 "Rule Set 1")

- **Validated:** 2026-06-17 (Phase-3 independent re-validation of the Doench-2014 enhancement, commit 129c2ca)   **Area:** MolTools
- **Canonical method(s):** `CrisprDesigner.CalculateOnTargetDoench2014(string context30Mer)` (`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs:521`); coefficient table `DoenchParams` + `DoenchIntercept`/`DoenchGcLow`/`DoenchGcHigh` (lines 473–502).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN — no defect found; no code/test change required.

This re-validation is scoped to the **Doench 2014 "Rule Set 1" on-target model** added in commit
129c2ca. The pre-existing honest-heuristic `CalculateGuideScore`/`DesignGuides` paths and the
MIT/Hsu off-target additions (CRISPR-OFF-001) were not the subject of this session.

## Stage A — Description

### Sources opened this session (independent of the repo)
- **Doench, Hartenian, Graham, et al.**, "Rational design of highly active sgRNAs for CRISPR-Cas9–mediated gene inactivation." *Nat Biotechnol* 32:1262–1267 (2014), PMID 25184501 — the model is a logistic-regression linear model over a fixed 30-nt context.
- **Reference implementation `doenchScore.py`** (CRISPOR, Haeussler et al. 2016, *Genome Biol* 17:148; github.com/maximilianh/crisporWebsite). **Re-downloaded the raw file this session** to `/tmp/doenchScore.py` (54 lines) and read it directly — did NOT rely on the repo's copied constants.

### Model layout (confirmed against the reference)
- **30-mer window:** `[4 nt 5′ context] + [20 nt protospacer] + [3 nt PAM] + [3 nt 3′ context] = 30`. The reference computes the GC term over `seq[4:24]` (= the 20-nt protospacer) and indexes all features by raw 0-based offset into the 30-mer (`subSeq = seq[pos:pos+len(modelSeq)]`). Confirmed.
- **Intercept:** `0.59763615` (reference line 30). ✔
- **GC-count term:** `score += abs(10 - gcCount) * gcWeight` with `gcWeight = gcLow (-0.2026259)` if `gcCount <= 10` else `gcHigh (-0.1665878)` (reference lines 31–32, 38–43). The form is `|GC − 10|` (absolute value) and the boundary is `<= 10 → gcLow`. Confirmed: the repo uses `Math.Abs(10 - gcCount)` and `gcCount <= 10 ? DoenchGcLow : DoenchGcHigh` (CrisprDesigner.cs:556–557) — identical.
- **Output:** `1.0/(1.0+exp(-score))` → probability in (0,1). The repo multiplies by 100 to expose a 0–100 scale. Confirmed and documented.

### Coefficient table — full independent diff
The repo's 70-entry `DoenchParams` table was parsed from the C# source and compared **tuple-by-tuple, in order** against the reference `params` list read from `/tmp/doenchScore.py`:

```
ref count: 70   repo count: 70
EXACT MATCH (ordered): True
In ref not repo: ∅
In repo not ref: ∅
```

There is **no transcription error**: every `(pos, seq, weight)` is byte-for-byte identical, including the reference's own intentional quirk where the dinucleotide tuples `(24,'AG'/'CG'/'TG')` reuse the single-nucleotide weights `(24,'A'/'C'/'T')` and `(26,'GT',0.11787758)` equals `(27,'T',0.11787758)`. Spot-checked weights explicitly against the reference file: `(22,'T',-0.8770074)`, `(23,'C',-0.8762358)`, `(11,'GG',-1.5169074)`, `(20,'GG',-0.7822076)`, `(21,'TC',-1.029693)`, `(28,'GG',-0.69774)`, `(29,'G',0.38634258)`, `(6,'C',-0.7411813)`, `(17,'G',-0.6780964)` — all match.

### Independent cross-check (worked examples reproduced from scratch)
Implemented the reference formula in Python from the coefficients (no trust of the repo's test literals) and ran the reference's own two self-test 30-mers:

| 30-mer | My computation ×1 | Reference published value | |
|---|---|---|---|
| `TATAGCTGCGATCTGAGGTAGGGAGGGACC` | 0.713089346955 | 0.713089368437 | Δ 2.2e-08 |
| `TCCGCACCTGTCACGGTCGGGGCTTGGCGC` | 0.018983843193 | 0.0189838463593 | Δ 3.2e-09 |

The tiny deltas are float-print precision in the reference's quoted literal; the model reproduces exactly.

### Edge-case semantics
- **Wrong length:** the reference assumes a 30-mer (no guard). The repo adds an explicit `length != 30 → ArgumentException` guard — a stricter, defensible contract for a typed API. Sourced as the model's documented input requirement.
- **Lowercase:** the repo upper-cases before scoring (`ToUpperInvariant`), so lowercase is accepted and gives the identical value. Reasonable extension.
- **Non-ACGT:** repo throws `ArgumentException`. The reference would silently produce a degenerate score; throwing is the safer typed contract.
- **0–100 vs 0–1:** the repo deliberately scales reference×100; documented in code and tests.

### Stage A divergences (all documented, none a defect)
- The repo adds an **NGG PAM check** (offsets 25–26 must be `GG`) that the reference does NOT enforce. This is a defensible SpCas9-specificity guard and is documented; it does not alter scoring for valid inputs. (Note: it restricts inputs the bare reference would still score; it is a guard, not a model change.)

**Stage A verdict: PASS** — layout, intercept, GC term form/boundary, output transform, and the full coefficient table all match the primary reference exactly; worked examples reproduced independently.

## Stage B — Implementation

- **Code path reviewed:** `CrisprDesigner.CalculateOnTargetDoench2014` (CrisprDesigner.cs:521–568) and `DoenchParams`/constants (473–502).
- **Formula realised correctly:** `score = intercept`; GC over `seq[4..24)` with `Math.Abs(10-gc)*gcWeight`, `gc<=10→gcLow`; feature loop adds `weight` when `string.CompareOrdinal(seq, pos, modelSeq, 0, len)==0` (= `seq[pos:pos+len]==modelSeq`); `1/(1+e^-score)*100`. Matches the reference line-for-line.
- **Cross-verification vs the actual C#:** the 11 Doench tests pass (`--filter ~Doench2014` → Passed 11, Failed 0). M-001/M-002 lock the reference worked examples ×100 (71.3089368437 / 1.89838463593, tol 1e-4). M-003 all-A-with-PAM expected 4.4338168085 — independently reproduced from the coefficients in Python: `4.4338168085440035`. ✔
- **Test-quality audit (HARD gate):** PASS, no green-washing.
  - M-001/M-002 are exact sourced cross-checks, not tautologies. Mutation test: perturbing the intercept by +0.001 shifts the score by ~0.020 (≫ tol 1e-4) — a wrong intercept/weight **would** fail these tests. (Implementer's own commit notes a 0.59763615→0.69763615 flip kills both; reconfirmed the sensitivity here.)
  - Edge tests genuinely cover wrong length (29/31/short), null/empty (`ArgumentNullException`), non-ACGT (`N`), lowercase equivalence, range [0,100] + ordering, and non-NGG PAM rejection.
- **Numerical robustness:** sigmoid is bounded, no overflow/precision concern over the coefficient range.

**Stage B verdict: PASS** — the code faithfully realises the validated reference model; tests assert exact externally-sourced values and a wrong constant fails.

## Verdict & follow-ups
- **End-state: ✅ CLEAN.** No transcription error in the 70-coefficient table, intercept, or GC terms; both reference worked examples reproduced independently; tests are real (exact sourced values, mutation-sensitive). No code or test change required this session.
- Residual (as of this 2026-06-17 Doench-2014 session): Doench Rule Set 2 / Azimuth (GBT model, no coefficient table) and CFD (binary pickle) remain unimplemented. **Both since cleared — CFD on 2026-06-17 (CRISPR-OFF-001) and Rule Set 2 / Azimuth on 2026-06-18 (see the "Rule Set 2 / Azimuth" section appended below). C7 fully resolved.**
- **Full unfiltered suite:** 6772 passed, 0 failed; build 0 errors (4 pre-existing NUnit2007 warnings in unrelated `ApproximateMatcher_EditDistance_Tests.cs`).

---

# Validation Report: CRISPR-GUIDE-001 (continued) — Doench 2016 "Rule Set 2" / Azimuth on-target score

- **Validated/Implemented:** 2026-06-18   **Area:** MolTools
- **Canonical methods added:** `CrisprDesigner.CalculateOnTargetRuleSet2(string context30Mer)` (sequence-only) and
  `CalculateOnTargetRuleSet2(string context30Mer, int aminoAcidCutPosition, double percentPeptide)` (gene-context),
  returning the Azimuth score in ~[0,1]. Engine: internal `AzimuthRuleSet2`; trained models embedded as
  `Resources/azimuth_rs2_nopos.bin` / `azimuth_rs2_full.bin`; reproducible extractor
  `scripts/azimuth/extract_azimuth_model.py`.
- **Stage A verdict:** PASS · **Stage B verdict:** PASS · **End-state:** ✅ CLEAN — **clears the LAST C7 residual.**

## Scope

Rule Set 2 is a trained scikit-learn `GradientBoostingRegressor` (100 trees, depth 3), **not** a coefficient table —
the reason it was the sole remaining C7 residual. It is reproduced faithfully from Microsoft Research's Azimuth
trained pickles (BSD-3-Clause), with **no scikit-learn dependency** at build or runtime. Governing rule unchanged:
every expected value is re-derived from an external, independently verified source, never read off the C# code.

## Stage A — Description

### Model recovery (sklearn-free)
The two Azimuth pickles (`V3_model_nopos.pickle`, `V3_model_full.pickle`) were decoded with a custom `Unpickler`
that reads the raw `sklearn.tree._tree.Tree` node arrays without instantiating scikit-learn. Recovered: 100 trees,
init (training-target mean) `0.5023237009327475`, learning-rate `0.1`, `n_features` 627 (nopos) / 630 (full),
1498 / 1496 total nodes. Prediction = `init + Σ_trees learning_rate · leaf_value`.

### Featurization reproduced EXACTLY (three-way verification)
1. **Featurizer == upstream:** a verbatim py3 port of `azimuth/features/featurization.py` run with **real
   Biopython** produces feature vectors element-for-element identical to our implementation (max |Δ| 1e-13) for
   the worst-case and representative guides — order-1/2 position-dependent + position-independent nucleotide
   one-hot/counts, GC features, NGGX, and the 4 melting-temperature features.
2. **Tm == real Biopython:** the nearest-neighbor melting temperature (DNA_NN3, salt-correction method 5,
   dnac1=dnac2=25 nM, Na=50 mM — the exact parameters azimuth passes to `MeltingTemp.Tm_NN`) matches the real
   installed Biopython `Tm_NN` to 4 decimals on the whole 30-mer **and** the short AT-rich sub-segments
   (e.g. `ATTTT` → −52.8234, `AGTTT` → −45.1363).
3. **Column order == CPython-2.7 dict order:** azimuth concatenates feature blocks in `dict.keys()` order, which at
   training time was the deterministic CPython-2.7 (64-bit, no hash randomization) iteration order. We reproduce it
   from first principles and validate the simulator against **documented** behaviour — `{'a','b','c'}` → `['a','c','b']`,
   `{'one','two','three'}` → `['three','two','one']` — and the py2 string hash against known values
   (`hash('a')=12416037344`, `hash('abc')=1453079729188098211`, `hash('foo')=−4177197833195190597`). The recovered
   orders are additionally consistent with the model's own split-threshold fingerprints (the 4 Tm columns carry
   continuous thresholds in [−54, 76]; count columns carry half-integer thresholds; one-hot columns only 0.5).

### Model traversal verified against scikit-learn itself
The extracted trees were reconstructed into **scikit-learn 1.6.1** `Tree` objects (node dtype migrated) and its own
`Tree.predict` was run on our feature matrix: **bit-identical** to our hand-written traversal (max |Δ| = 0.0) across
all 947 guides. Extraction + traversal are therefore provably correct.

### KEY FINDING — the upstream fixture is stale, not a defect in our code
Upstream ships `azimuth/tests/1000guides.csv` (`truth nopos` / `truth pos`). Our faithful pipeline reproduces it on
only **585/947** (nopos) and **637/947** (full) rows; the rest differ by ≤0.04. This is **not** an implementation
error: (a) scikit-learn's own `Tree.predict` agrees with us bit-for-bit on the same features; (b) our featurizer is
1e-13-identical to a verbatim upstream featurization with real Biopython; (c) no feature ordering reproduces all 947
(exhaustively searched). The upstream test file itself warns it "can fail due to randomness ... feature reordering" —
the fixture drifted from the shipped pickles. **Consequently the authoritative oracle is our verified reference (≡
`azimuth.model_comparison.predict` for the shipped model), and the upstream fixture is used only as independent
third-party corroboration on the agreeing subset.** A spurious "max-agreement" full-model order (651) that fit the
stale fixture slightly better was rejected in favour of the py2-correct order (637) after pinning the CPython-2.7
hash/probe behaviour (one binary `gc_below_10` split distinguishes them).

### Independent cross-checks (from the verified reference + upstream)
| Check | Result |
|---|---|
| C# == verified reference, all 947 nopos guides | max |Δ| < 1e-5 ✓ |
| C# == verified reference, all 947 full guides | max |Δ| < 1e-5 ✓ |
| C# == upstream `truth nopos` on the 585 agreeing rows | all < 1e-3 ✓ |
| C# == upstream `truth pos` on the 637 agreeing rows | all < 1e-3 ✓ |
| Score range over all guides | [0,1] ✓ |

**Stage A verdict: PASS** — model recovered sklearn-free; featurization 1e-13-identical to upstream incl. real
Biopython Tm; column order matches documented CPython-2.7 dict behaviour; traversal bit-identical to scikit-learn;
the upstream fixture's partial disagreement proven to be its own drift, not our defect.

## Stage B — Implementation

- **Code path:** `AzimuthRuleSet2.cs` (model reader via `MemoryMarshal` over the embedded blob, featurizer, GBRT
  traversal) and the two thin public wrappers in `CrisprDesigner.cs` (Rule Set 2 region). Additive; no existing
  method/signature/test changed.
- **Binary format:** little-endian header (magic `ARS2`, version, flags, treeCount, nodeCount, nFeatures, init f64,
  learning-rate f64) + `treeStart[]` + AoS nodes (24 bytes: f64 threshold-or-value, i32 left/right/feature/pad).
  Leaf values are pre-scaled by the learning rate; f64 thresholds preserve the pickle's exact split points. Nodes are
  read zero-copy and copied once into a managed array; models load lazily and are thread-safe. Round-trips the
  reference to < 5e-7 (limited only by the 6-decimal oracle CSV).
- **Featurizer realised correctly:** writes each block at its py2-dict offset (nopos 627: gc_count | pd2 | pd1 |
  gc_above | pi1 | pi2 | Tm | gc_below | NGGX; full 630: same with AA-cut, pct-peptide, pct<50% interleaved). Tm via
  the DNA_NN3 nearest-neighbor model. ACGT-only, NGG-PAM-at-25/26, length-30 contract enforced.
- **Tests:** new fixture `CrisprDesigner_RuleSet2_Tests` (13 tests): the two model-vs-reference sweeps over all 947
  guides, the two upstream-agreeing-subset corroborations (count locked at 585 / 637), unit-interval, case-insensitive,
  determinism, full-vs-nopos differs, and the full validation/error set (null/empty, wrong length ×2, non-ACGT,
  missing NGG PAM). Every numeric oracle is the verified reference or upstream — none read off the C# code.
- **Defects:** none.

**Stage B verdict: PASS.**

## Verdict & follow-ups

- **Stage A: PASS · Stage B: PASS · End-state: ✅ CLEAN.**
- **C7 is now fully resolved** — all four CRISPR scoring models (Doench 2014 Rule Set 1, MIT/Hsu 2013, CFD 2016, and
  Doench 2016 Rule Set 2 / Azimuth) are implemented and validated. No CRISPR scoring residual remains.
- Reproducibility: `python3 scripts/azimuth/extract_azimuth_model.py` re-downloads the pinned Azimuth pickles and
  regenerates both `.bin` blobs and the embedded oracle CSVs.
- Test spec: the four published CRISPR scoring models (Rule Set 1, Rule Set 2 / Azimuth, MIT/Hsu, CFD) are
  consolidated in `tests/TestSpecs/CRISPR-SCORE-001.md`.
- Full unfiltered Genomics suite: **6825 passed, 0 failed**; build 0 errors, 0 warnings.
