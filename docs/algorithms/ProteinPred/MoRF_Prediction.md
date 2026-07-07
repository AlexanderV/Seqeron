# MoRF (Molecular Recognition Feature) Prediction

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinPred |
| Test Unit ID | DISORDER-MORF-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

A Molecular Recognition Feature (MoRF) is a short protein segment, located **within** a longer
intrinsically disordered region, that undergoes a disorder-to-order transition upon binding a partner
protein [1][2]. In a per-residue disorder prediction profile a MoRF appears as a downward "dip" — a short
stretch of relative **order** flanked by disorder [3]. `PredictMoRFs` detects these dips: maximal runs of
residues whose disorder score falls below the 0.5 order/disorder threshold [3], whose length lies in the
10–70 residue band [1], and which are flanked by predicted disorder on both sides [1][3]. It is an
evidence-driven heuristic annotator over amino-acid composition, not a trained machine-learning predictor.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Intrinsically disordered regions (IDRs) lack stable tertiary structure in isolation. A subset of IDRs
contains short interaction-prone elements that fold upon binding. Mohan et al. (2006) named these
Molecular Recognition Features (MoRFs) and classified them by bound-state secondary structure into
α-MoRFs (α-helix), β-MoRFs (β-strand), and ι-MoRFs (irregular) [1]. They are "relatively short
(10-70 residues), loosely structured protein regions" within "longer, largely disordered sequences" [1].

### 2.2 Core Model

Let `d(i) ∈ [0,1]` be the per-residue disorder score (higher = more disordered), computed by
`PredictDisorder` from the normalized TOP-IDP scale [4]. A residue is predicted **ordered** when
`d(i) < 0.5` and **disordered** when `d(i) ≥ 0.5` [3]. A MoRF (dip) is a maximal interval `[s, e]` with:

1. `d(i) < 0.5` for all `i ∈ [s, e]` (the ordered dip) [3];
2. `10 ≤ (e − s + 1) ≤ 70` (Mohan length band) [1];
3. `d(s−1) ≥ 0.5` and `d(e+1) ≥ 0.5` (embedded within disorder; flanked on both sides) [1][3].

The reported score is the dip depth below the threshold, normalized to `[0,1]`:

`score = (0.5 − mean_{i∈[s,e]} d(i)) / 0.5`, clamped to `[0,1]`.

Because every `d(i) < 0.5` inside the dip, `mean d(i) ∈ [0, 0.5)` and `score ∈ (0, 1]` — a deeper dip
(more ordered) scores higher. The 0.5 divisor is the maximum possible depth, so the normalization is a
direct derivation from the 0.5 threshold [3], not a tuned constant.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 ≤ Start ≤ End < len(sequence)` | Coordinates taken from valid residue indices |
| INV-02 | `minLength ≤ (End − Start + 1) ≤ maxLength` | Length filter on each dip [1] |
| INV-03 | Each MoRF is flanked by `d ≥ 0.5` on both immediate sides (within disorder) | Flank check [1][3] |
| INV-04 | Reported MoRFs are non-overlapping and ordered by `Start` | Maximal disjoint runs scanned left→right |
| INV-05 | `0 ≤ Score ≤ 1`, monotone in dip depth | Clamped normalization by the 0.5 threshold [3] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | This dip heuristic | Trained predictors (MoRFchibi, ANCHOR2) |
|--------|--------------------|------------------------------------------|
| Input features | Single composition scale (TOP-IDP) | ML over multiple features / energy models |
| Output | Dip coordinates + depth score | Per-residue MoRF propensity |
| Accuracy | Coarse first-pass | Substantially higher [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | Protein sequence | Case-insensitive; non-standard chars contribute 0 propensity |
| minLength | int | 10 | Minimum MoRF length | Residues; Mohan band lower bound [1] |
| maxLength | int | 70 | Maximum MoRF length | Residues; Mohan band upper bound [1] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Start | int | 0-based inclusive start of the dip |
| End | int | 0-based inclusive end of the dip |
| Score | double | Dip depth normalized to [0,1] (deeper order = higher) |

### 3.3 Preconditions and Validation

`null` or empty input returns an empty sequence. Input is upper-cased internally. Indexing is 0-based,
inclusive on both ends. Residues outside the 20 standard amino acids contribute 0 disorder propensity via
`PredictDisorder`. No exceptions are thrown for valid string input.

## 4. Algorithm

### 4.1 High-Level Steps

1. Compute per-residue disorder scores via `PredictDisorder` (normalized TOP-IDP) [4].
2. Scan left→right for maximal runs where `d(i) < 0.5` (ordered dips) [3].
3. Discard runs whose length is outside `[minLength, maxLength]` [1].
4. Discard runs not flanked by `d ≥ 0.5` on both immediate sides (not embedded in disorder) [1][3].
5. Emit each surviving dip with score `(0.5 − mean d) / 0.5`, ordered by start.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

| Constant | Value | Source |
|----------|-------|--------|
| Order/disorder threshold | 0.5 | Cheng/Oldfield PMC2570644 [3] |
| MoRF length band | 10–70 residues | Mohan et al. 2006 [1] |
| Per-residue scores | normalized TOP-IDP scale | Campen et al. 2008 [4] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `PredictMoRFs` | O(n·w) | O(n) | `n` = length, `w` = disorder window (21); the dip scan itself is O(n) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)

- `DisorderPredictor.PredictMoRFs(string, int, int)`: detects dip-in-disorder MoRFs.
- `DisorderPredictor.PredictDisorder(...)`: supplies the per-residue disorder profile consumed above.

### 5.2 Current Behavior

The dip scan operates on the per-residue scores already produced by `PredictDisorder`; it does not
re-window. The MoRF dip is defined against the **0.5** order/disorder threshold from the MoRF literature
[3], which is independent of the TOP-IDP decision cutoff (0.542) `PredictDisorder` uses for general IDR
calling — the two thresholds serve different purposes. No substring/pattern search is involved, so the
repository suffix tree is **not applicable** to this unit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- MoRF = short region of order ("dip", `d < 0.5`) within a longer region of disorder [3].
- Length band 10–70 residues [1].
- Embedding requirement: flanked by predicted disorder on both sides [1][3].

**Intentionally simplified:**

- Exact dip flank/run-length detection parameters: Oldfield et al. (2005) defines precise numeric
  parameters in a paywalled Methods section that could not be retrieved; **consequence:** the flank rule
  uses "≥1 disordered residue on each immediate side" rather than a published fixed flank length.
- Score model: dip depth below the 0.5 threshold rather than a trained propensity; **consequence:** scores
  rank dips by relative order, not calibrated binding probability.

**Not implemented:**

- MoRF subtype classification (α/β/ι); **users should rely on:** structural data or dedicated predictors
  (MoRFchibi, ANCHOR2).
- Trained machine-learning MoRF propensity; **users should rely on:** MoRFchibi / ANCHOR2 for
  publication-grade prediction [1].

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Flank length detail | Assumption | Boundary residue inclusion at dip edges | accepted | Oldfield 2005 Methods paywalled; threshold 0.5 and 10–70 band are source-traceable |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Fully ordered sequence | No MoRFs | No surrounding disorder to embed a dip [3] |
| Fully disordered sequence | No MoRFs | No ordered dip exists [3] |
| Dip < 10 or > 70 residues | Not reported | Outside Mohan length band [1] |
| Dip at sequence terminus | Not reported | Not flanked by disorder on both sides [1][3] |
| null / empty / very short | Empty result | Cannot contain a 10-residue embedded dip |

### 6.2 Limitations

Single-scale composition heuristic: it cannot match trained predictors, does not assign MoRF subtype, and
its boundaries depend on the disorder-window smoothing in `PredictDisorder`. Suitable for coarse first-pass
annotation and education, not clinical or publication-grade MoRF calling.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
// 25 disordered P + 30 ordered L + 25 disordered P → one MoRF (dip) inside the IDR.
string seq = new string('P', 25) + new string('L', 30) + new string('P', 25);
var morfs = DisorderPredictor.PredictMoRFs(seq).ToList();
// morfs[0] == (Start: 29, End: 50, Score: ~0.2759)
```

**Numerical / biological walk-through:**

For the sequence above, the smoothed disorder profile dips below 0.5 over residues 29–50 (length 22, within
10–70). Mean disorder over the dip ≈ 0.362033, so the score = (0.5 − 0.362033) / 0.5 ≈ 0.275934. Replacing
the L core with the more order-promoting I (TOP-IDP −0.486) deepens the dip (mean ≈ 0.300196,
score ≈ 0.399608), confirming score monotonicity in dip depth (INV-05).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [DisorderPredictor_MoRF_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/DisorderPredictor_MoRF_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [DISORDER-MORF-001-Evidence.md](../../../docs/Evidence/DISORDER-MORF-001-Evidence.md)
- Related algorithms: [Disorder_Prediction](./Disorder_Prediction.md)

## 8. References

1. Mohan A, Oldfield CJ, Radivojac P, Vacic V, Cortese MS, Dunker AK, Uversky VN. 2006. Analysis of molecular recognition features (MoRFs). J Mol Biol 362(5):1043–1059. https://doi.org/10.1016/j.jmb.2006.07.087 (PMID 16935303)
2. Wikipedia contributors. Molecular recognition feature. https://en.wikipedia.org/wiki/Molecular_recognition_feature (accessed 2026-06-14)
3. Cheng Y, Oldfield CJ, Meng J, Romero P, Uversky VN, Dunker AK. Mining α-helix-forming molecular recognition features with cross species sequence alignments. Biochemistry. https://pmc.ncbi.nlm.nih.gov/articles/PMC2570644/
4. Campen A, Williams RM, Brown CJ, Meng J, Uversky VN, Dunker AK. 2008. TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder. Protein Pept Lett 15(9):956–963. https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/
5. Oldfield CJ, Cheng Y, Cortese MS, Romero P, Uversky VN, Dunker AK. 2005. Coupled folding and binding with α-helix-forming molecular recognition elements. Biochemistry 44(37):12454–12470. https://pubmed.ncbi.nlm.nih.gov/16156658/
