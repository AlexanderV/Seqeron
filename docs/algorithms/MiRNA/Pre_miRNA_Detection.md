# Pre-miRNA Hairpin Detection

| Field | Value |
|-------|-------|
| Algorithm Group | MiRNA |
| Test Unit ID | MIRNA-PRECURSOR-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Production (default MFE fold, F40) + opt-in consecutive-pairing screen + opt-in trained classifier |
| Last Reviewed | 2026-11-25 (B13 review: F8–F10, F40) |

## 1. Overview

Pre-miRNA detection scans RNA sequence windows for stem-loop candidates that satisfy basic precursor-miRNA length, stem, and loop constraints [1][2][4]. The repository implements this in `MiRnaAnalyzer.FindPreMiRnaHairpins` by enumerating candidate windows, requiring uninterrupted complementary pairing from both ends inward, extracting a 5' mature strand and 3' star strand, and assigning the Turner 2004 free energy of the reported stem-loop computed by the canonical `RnaSecondaryStructure` loop model. The search is deterministic and exhaustive within the caller's length bounds, but it is intentionally stricter than real pre-miRNA structure because internal bulges, asymmetric loops, and arm switching are not modeled. This makes the implementation suitable for exact educational hairpin screening rather than production-grade precursor discovery on natural miRBase sequences [1][4][5].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Pre-miRNAs are stem-loop RNA intermediates in miRNA biogenesis, typically on the order of roughly 60-120 nt, from which Dicer processes a mature miRNA duplex of about 22 nt [1][2]. A valid precursor must fold into a hairpin with a stem and terminal loop, and the mature miRNA resides on one arm of that structure [1][2]. Thermodynamic stability of RNA hairpins is commonly modeled with nearest-neighbor parameters, including base-stacking energies, loop penalties, and terminal mismatch terms [3].

### 2.2 Core Model

For a candidate RNA subsequence $S$ of length $n$, the repository accepts a precursor candidate only when it can be decomposed into:

$$
S = \text{5' stem} + \text{loop} + \text{3' stem}
$$

with uninterrupted pairing between mirrored positions from the two ends inward (extended as far as the bases pair, but never closing a loop shorter than the 3-nt NNDB minimum) under the pairing set:

$$
\{A\!-\!U,\; U\!-\!A,\; G\!-\!C,\; C\!-\!G,\; G\!-\!U,\; U\!-\!G\}
$$

The acceptance criteria of the opt-in screen (`useConsecutiveScreen: true`) are:

- stem length $\ge 18$ bp,
- loop length $= n - 2 \times \text{stem length}$ with $3 \le \text{loop length} \le 25$,
- candidate length constrained by the caller, defaulting to 55-120 nt.

The free energy is that of the reported structure under the full Turner 2004 model [3], computed by the
canonical stem-loop evaluator (`RnaSecondaryStructure.FindStemLoops` with the loop fixed →
`CalculateStemEnergy` + `CalculateHairpinLoopEnergy`):

$$
\Delta G = \sum \Delta G_{\text{stack}} + 0.45 \times n_{\text{AU/GU helix ends}} + \Delta G_{\text{hairpin}}
$$

where $\Delta G_{\text{hairpin}}$ is the NNDB hairpin-loop term: the experimental total for special
tri/tetra/hexaloops, otherwise initiation(n) + terminal mismatch + UU/GA and GG first-mismatch bonuses +
special GU closure + all-C penalty (triloops: initiation + all-C penalty only). The AU/GU end penalty is
applied to both helix ends, including the pair closing the loop (NNDB 2004 hairpin example 5.3). For a
perfect hairpin this structure is the MFE structure and the energy equals
`RnaSecondaryStructure.CalculateMinimumFreeEnergy`.

#### 2.2.1 MFE-structure-based detection (the default since F40)

`FindPreMiRnaHairpins` (default, `useConsecutiveScreen = false`) is `FindPreMiRnaHairpinsByMfe` with each accepted
window mapped onto the `PreMiRna` record (below); `AssessHairpinByMfe` / `FindPreMiRnaHairpinsByMfe` replace the consecutive-pairing scan with the
**real** minimum-free-energy secondary structure produced by the validated Zuker–Stiegler folder
(`RnaSecondaryStructure.CalculateMfeStructure`, RNA-STRUCT-001, Turner 2004 NN model) [3][6]. The
candidate is folded once; its hairpin features are read from the ACTUAL MFE dot-bracket:

- the structure must be a **single dominant hairpin** — exactly one terminal (apical) loop, a nested
  stem that may contain small internal loops/bulges, and **no multibranch** (no `(` after a `)`) —
  consistent with the fold-back single-arm duplex of the annotation criteria [2][8];
- **stem base pairs $\ge 16$** (number of pairs in the MFE structure) — Ambros et al. (2003) [2];
- **terminal loop $\in [3, 25]$ nt** — Bartel (2004) [1];
- **MFEI $\ge 0.85$** — Zhang et al. (2006) [7], where

$$
\text{AMFE} = \frac{100 \cdot |\Delta G^\circ|}{n}, \qquad
\text{MFEI} = \frac{\text{AMFE}}{(G\!+\!C)\%}.
$$

Because genuine pre-miRNAs fold to a markedly more stable structure than random sequences [6], this
path detects natural miRBase precursors (e.g. hsa-mir-21, hsa-let-7a-1) that the consecutive-pairing
heuristic rejects. $\Delta G^\circ$ is taken from the engine's `CalculateMinimumFreeEnergy`
verbatim; MFEI uses $|\Delta G^\circ|$ so the published "MFEI > 0.85" threshold applies to the
folder's negative $\Delta G^\circ$.

**Mapping onto `PreMiRna` (default mode).** `Structure` / `FreeEnergy` are the canonical MFE dot-bracket and
$\Delta G^\circ$ of the window. A window carries no basal-junction knowledge, so the Drosha ruler of
`PredictDroshaDicerCleavage` cannot be applied; the window is read as the Drosha product (the pre-miRNA):
the **mature** (5p) product starts at the 5′ base of the outermost pair of the MFE stem and is
`matureLength` (22) nt long (Dicer 5′ counting, Park et al. 2011), capped at the 5′ arm; the **star** is the
opposite-arm product of the Dicer 2-nt-3′-overhang rule — the same miRDeep2 `find_star` code as
`PredictDroshaDicerCleavage` (Friedländer et al. 2008): start = partner(mature end − 2 − off) − off, end =
partner(mature start) + 2, clipped to the window. These are arm-dominant stem ends, not an annotated
guide/passenger call; for processing-accurate coordinates of a known pri-miRNA use `PredictDroshaDicerCleavage`
with the junction. Oracle (ViennaRNA 2.8 `RNA.fold`, 5 miRBase pri-miRNA hairpins, all 554 accepted windows):
identical dot-bracket in every window, |ΔG° − Vienna| ≤ 0.53 kcal/mol (best windows: mir-21 −35.17/−34.70,
let-7a-1 −34.31/−34.20, mir-16-1 −37.05/−37.00, mir-155 −29.76/−29.50, mir-30a −37.37/−37.20).
**Cost:** every window is folded (≈ 20 ms for 55–120 nt, O(#windows · W³)); the scan is parallelised over window
starts with deterministic ordered output. `useConsecutiveScreen: true` restores the fast screen.

#### 2.2.2 Opt-in Drosha/Dicer cleavage-site prediction (default unchanged)

`PredictDroshaDicerCleavage(sequence, basalJunction[, dotBracket, dicerRule])` predicts the excision
coordinates of a pri-miRNA hairpin by applying the PUBLISHED measuring rules to its folded stem (the
caller's structure, e.g. RNAfold output, or the canonical MFE fold) — it does **not** use a trained
classifier:

- **Drosha (basal-junction ruler):** Drosha cleaves $\sim 11$ bp ($\approx$ one helical turn) from the
  basal ssRNA–dsRNA junction [9] ("an 11 bp basal stem was preferred" for all four pri-miRNAs selected in [12]). The count is in **base pairs** of the unbranched
  stem starting at the junction pair (bulged nucleotides are skipped): the 5' end of the pre-miRNA is the
  5' base of the 12th pair, its 3' end is that base's partner $+2$ (2-nt 3' overhang [9][11]).
- **Dicer:** 5' counting (default) — the 5p product is $22$ nt from the Drosha 5' end [10]; 3' counting —
  the 3p product is $22$ nt from the Drosha 3' end [10]. The opposite product end follows the 2-nt
  3'-overhang pairing rule of miRDeep2 `find_star` [17] (start = partner(end − 2), compensated for bulged
  anchors; end = partner(start) + 2) — the same rule Auyeung et al. used to infer Microprocessor sites
  from RNAfold structures [12].
- **Motifs (confidence flags, not required) [12]:** positions are numbered from the Drosha sites
  (negative = upstream of the 5p cut, positive = downstream of the 3p cut). `HasBasalUgMotif` = U at −14
  and G at −13; `HasCnncMotif` = a C-N-N-C starting at +16..+18. The apical UGU/GUG motif is not flagged:
  the paper gives only "enrichment in this region of the loop" (window in Fig. S7F; P22–P24 in Fig. 7A vs
  UGUG at P24–P27 for mir-30a, which RNAfold pairs in the apical stem), so there is no reproducible rule.
- Gu et al. 2012 "loop counting" is a precision condition (Dicer cuts accurately when the cut is 2 nt from
  a loop/bulge), not a rule that places a cut; it is not a counting option.

Cross-check (miRBase hairpins + hg38 flanks from UCSC, RNAfold 2.8 structure, junction = the pair 11 bp
below the annotated 5p start; offsets predicted − annotated, 5p start/end | 3p start/end, 5' counting):
hsa-mir-21 0/−2 | +1/−1; hsa-let-7a-1 +1/0 | −1/0; hsa-mir-16-1 0/0 | 0/+1; hsa-mir-30a 0/0 | +2/0;
hsa-mir-155 0/−2 | +2/0. miRBase ends are read-defined (several 5p entries are 23–24-nt isoforms). The
11-bp ruler itself: the lower stem of the folded pri-miRNA (from the annotated 5p start down to the first
long single-stranded segment) is 8–13 bp, i.e. ~11 bp as published. The basal UG and CNNC flags are set
for mir-16-1 and mir-30a, the two of these five that Auyeung et al. name as CNNC examples.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | (Screen mode only, `useConsecutiveScreen: true`.) Candidate stems are approximated by uninterrupted end-to-end pairing with no internal bulges or mismatches. The default mode does not make this assumption (MFE fold, F40). | Natural pre-miRNAs with internal loops or asymmetric bulges (e.g. hsa-mir-21, hsa-let-7a-1) are missed by the screen. |
| ASM-02 | (Screen mode.) The reported energy is the full Turner 2004 energy of the reported (consecutive-stem) structure, not the MFE of the window; the default mode reports the MFE. | For windows whose MFE fold differs (bulges, internal loops), the energy is that of the consecutive stem; use the opt-in MFE path for the optimal fold. |
| ASM-03 | The mature product is extracted from the 5' arm (default mode: from the outermost pair, 22 nt, star by the 2-nt-overhang find_star rule; screen mode: first/last `matureLength` nt of the window, blunt). | Valid precursors whose dominant mature product comes from the opposite arm are represented asymmetrically. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported precursor has a balanced dot-bracket structure string of the same length as `Sequence`. | The structure is the canonical rendering of the stem-loop (first `stemLength` positions `(`, loop `.`, last `stemLength` `)`). |
| INV-02 | Every reported precursor has stem length `>= 18`. | Candidates with shorter uninterrupted stems are rejected before record creation. |
| INV-03 | Every reported precursor has loop length between `3` and `25` nt. | `AnalyzeHairpin` explicitly rejects loops outside that range. |
| INV-04 | `Start` and `End` are zero-based inclusive coordinates in the scanned input sequence. | `FindPreMiRnaHairpins` constructs `End = Start + length - 1`. |
| INV-05 | `MatureSequence.Length = StarSequence.Length = min(matureLength, stemLength)`. | Both arms are extracted using the same bounded length. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | RNA or DNA sequence to scan for precursor candidates. | `null`, empty, or shorter-than-`minHairpinLength` input yields no results; `T` is normalized to `U`. |
| `minHairpinLength` | `int` | `55` | Minimum candidate window length considered during scanning (the only length floor). | Must be ≥ 1 (else `ArgumentOutOfRangeException`, thrown eagerly); values larger than the input length produce no results. |
| `maxHairpinLength` | `int` | `120` | Maximum candidate window length considered during scanning. | Values below `minHairpinLength` make the scan empty. |
| `matureLength` | `int` | `22` | Maximum length extracted for mature and star strands. | Must be ≥ 1 (else `ArgumentOutOfRangeException`, thrown eagerly). |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Start` | `int` | Zero-based inclusive start coordinate of the candidate in the original input sequence. |
| `End` | `int` | Zero-based inclusive end coordinate of the candidate in the original input sequence. |
| `Sequence` | `string` | Uppercase RNA subsequence that satisfied the precursor tests. |
| `MatureSequence` | `string` | Extracted mature strand from the 5' arm. |
| `StarSequence` | `string` | Extracted star strand from the mirrored 3' arm. |
| `Structure` | `string` | Dot-bracket representation with contiguous stem parentheses and loop dots. |
| `FreeEnergy` | `double` | Turner 2004 ΔG° (kcal/mol) of `Structure`, from the canonical `RnaSecondaryStructure` stem-loop evaluator. |

### 3.3 Preconditions and Validation

`FindPreMiRnaHairpins` is non-throwing for `null`, empty, or too-short input: those cases simply yield no candidates. `minHairpinLength < 1` or `matureLength < 1` throws `ArgumentOutOfRangeException` at call time (eager validation). Input is normalized to uppercase RNA by replacing `T` with `U`. Candidate windows are rejected if uninterrupted complementary pairing from the ends inward yields a stem shorter than 18 bp, or if the resulting loop exceeds 25 nt (the stem never closes a loop below 3 nt). The window length is bounded only by the caller's `minHairpinLength`/`maxHairpinLength` (the former hard-coded 55-nt floor was removed). Coordinates are zero-based and inclusive. `FindPreMiRnaHairpinsByMfe` validates `minHairpinLength ≥ 1` the same way.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input to uppercase RNA.
2. Slide a start position across the sequence and enumerate every candidate window whose length is between `minHairpinLength` and `maxHairpinLength`.
3. For each candidate window, count uninterrupted complementary base pairs from the two ends inward (Watson-Crick and G:U wobble, `RnaSecondaryStructure.CanPair`), stopping at the first non-pair or when one more pair would leave a loop shorter than 3 nt.
4. Reject the candidate unless the resulting uninterrupted stem is at least 18 bp and the remaining loop is at most 25 nt.
5. Extract the mature strand from the 5' arm and the star strand from the mirrored 3' arm.
6. Obtain the dot-bracket and the Turner 2004 energy from the canonical stem-loop (`RnaSecondaryStructure.FindStemLoops` with the loop fixed).
7. Emit a `PreMiRna` record with genomic coordinates relative to the scanned input sequence.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Hairpin stability is delegated to the canonical Turner 2004 evaluator in `RnaSecondaryStructure`
(RNA-STRUCT-001): stacking (incl. the GGUC/CUGG special stack), 0.45 kcal/mol AU/GU penalty per helix end,
and the full hairpin-loop term (special loops, terminal mismatch, first-mismatch bonuses, special GU
closure, all-C penalty). Reference cross-check (B13 F8): ViennaRNA 2.8 `eval_structure` agrees to
≤ 0.24 kcal/mol on the test fixtures — the residual is ViennaRNA's 0.1-kcal-rounded parameter set
(stacks, TerminalAU = 0.50) — and the hairpin-loop terms are identical (e.g. UUCG tetraloop closed by
C-G: 3.70 in both).

The pairing test accepts both Watson-Crick and G:U wobble pairs, but the uninterrupted stem count stops at the first non-pairing mirrored position.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Exhaustive precursor scan | `O(nH^2)` | `O(H)` | `n` is input length and `H` is the maximum scanned hairpin length, because up to `O(nH)` windows are enumerated and each candidate can cost `O(H)` to analyze. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MiRnaAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs)

- `MiRnaAnalyzer.FindPreMiRnaHairpins(string, int, int, int, bool useConsecutiveScreen = false, double minMfei = 0.85)`: Public precursor scan over all candidate windows (default MFE fold via `FindPreMiRnaHairpinsByMfe` + `ToPreMiRna`; `useConsecutiveScreen: true` = consecutive-pairing screen).
- `MiRnaAnalyzer.AnalyzeHairpin(string, int)`: Private validator of the screen mode that checks stem continuity, loop size, arm extraction, and structure generation.
- `MiRnaAnalyzer.FindPreMiRnaHairpinsByMfe(string, int, int, double, int)`: Opt-in precursor scan that folds each candidate window with the RNA-STRUCT-001 MFE engine.
- `MiRnaAnalyzer.AssessHairpinByMfe(string, double, int)`: Opt-in single-candidate assessment from the real MFE structure (single hairpin + stem≥16 + loop 3–25 + MFEI≥0.85).
- `MiRnaAnalyzer.CalculateMfeIndex(double, int, double)`: MFEI = AMFE/(G+C)%, AMFE = 100·|ΔG°|/length (Zhang 2006).
- `MiRnaAnalyzer.PredictDroshaDicerCleavage(string, int[, string?, DicerCountingRule])`: Opt-in cleavage-site prediction on the folded stem — Drosha cut (11 bp from the basal junction), Dicer cut (22 nt, 5' or 3' counting), 5p + 3p products on opposite arms with 2-nt 3' overhangs, basal UG and CNNC flags (Han 2006 / Park 2011 / miRDeep2 find_star / Auyeung 2013).
- `MiRnaAnalyzer.ExtractPreMiRnaFeatures(string, int)`: Opt-in extraction of the published structure/sequence features (MFE, AMFE, MFEI, GC, %paired = 2·dP, dP = |bp|/L [16], stem bp, loop, length) from the real MFE structure (RNA-STRUCT-001).
- `MiRnaAnalyzer.ClassifyPreMiRna(string, double, int)`: Opt-in **trained** logistic-regression natural-vs-background classifier (an **in-repo model, not a published classifier**); returns the features, P(natural), and the boolean call.
- `MiRnaAnalyzer.ClassifyPreMiRnaTripletSvm(string, string?, int, int)`: Opt-in **published** classifier — Triplet-SVM (Xue et al. 2005, listed in ref. 13) with the authors' own trained libsvm RBF model (163 human pre-miRNAs vs 168 pseudo hairpins, 141 support vectors, γ = 16, ρ = 2.74269); ports the package's stem-loop retrieval (≥ 22 paired bases per arm) and 32 triplet-element coding, evaluated with the bundled support vectors (`MiRnaAnalyzer.TripletSvmModel.cs`, source URL + SHA-256 in its header). Reproduces the package's own 40-example run (37 positive/3 negative) and the Perl + libsvm output on Seqeron MFE folds; the only difference from the package is the folder (this library's MFE instead of RNAfold 1.x).
- `MiRnaAnalyzer.ScorePreMiRnaFeatures(PreMiRnaFeatures)`: Opt-in evaluation of the bundled trained model on a feature vector.
- `MiRnaAnalyzer.DinucleotideShuffle(string, Random)`: Altschul-Erickson (1985) dinucleotide-preserving shuffle (uniform over all sequences with the same doublet counts and ends), last-exit arborescence drawn by Wilson's walk as in uShuffle (Jiang et al. 2008) [15]; draw-for-draw identical to `ushuffle.c` — the background-set generator (Bonnet 2004 / randfold `-d` convention).
- `MiRnaAnalyzer.CalculateRandfoldPValue(string, int, Random, int)`: Opt-in randfold test (Bonnet 2004 [6]; `randfold.c` counting): P = (1 + #{shuffle ΔG° ≤ ΔG°})/(N + 1) over N dinucleotide shuffles.

### 5.2 Current Behavior

By default (F40) every window is folded with the canonical MFE engine and accepted by the Ambros/Bartel/Zhang criteria (§2.2.1). In `useConsecutiveScreen: true` mode the implementation searches candidates exhaustively rather than folding. Stem detection is based on uninterrupted end-to-end pairing, so the counter stops at the first mirrored mismatch and does not resume across bulges or internal loops; it is not truncated otherwise (the former unsourced `min(n/2 − 5, 35)` cap, which forced every reported loop to ≥ 10 nt and rejected perfect stems > 35 bp, was removed in B13 F9). `MatureSequence` is always taken from the 5' arm beginning at the candidate start, and `StarSequence` is the mirrored 3' arm segment it base-pairs with (a blunt duplex label, not Dicer's 2-nt 3'-overhang product — see `PredictDroshaDicerCleavage` for processing coordinates). Structure and energy come from the canonical `RnaSecondaryStructure` stem-loop (B13 F8; previously an inline stacking + initiation + mismatch sum without special loops or first-mismatch bonuses).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- RNA hairpin validation using Watson-Crick and G:U wobble pairing [1][4].
- Precursor-style constraints on hairpin length, stem length, and loop size [1][2][4].
- Full Turner 2004 energy of the reported stem-loop via the canonical `RnaSecondaryStructure` evaluator [3].
- **Opt-in:** real MFE-structure-based hairpin assessment that folds the candidate with the
  validated Zuker–Stiegler engine (RNA-STRUCT-001) [3][6] and applies the Ambros (2003) ≥16-bp
  stem [2], Bartel (2004) 3–25-nt loop [1], and Zhang (2006) MFEI ≥ 0.85 [7] criteria to the actual
  MFE structure — tolerant of internal bulges/loops and detecting natural miRBase precursors.
- **Opt-in:** Drosha/Dicer cleavage-site prediction (`PredictDroshaDicerCleavage`) applying the published
  measuring rules to the folded stem — Drosha 11 bp (pairs) from the basal junction [9][12], Dicer ~22 nt
  5'- or 3'-counting [10], the RNase III 2-nt 3' overhang on the opposite arm [9][11][17], and the
  Auyeung (2013) basal UG and CNNC confidence motifs [12].
- **Opt-in:** a **trained** structure/sequence-feature natural-vs-background classifier
  (`ClassifyPreMiRna`). Logistic regression [13] over published features — MFE [6], AMFE/MFEI [7],
  GC%, and %paired (= 2·dP, base-pairing propensity [16]) — ridge-penalised maximum likelihood (λ = 1e-3)
  solved by Newton–Raphson/IRLS [13] (= scikit-learn `LogisticRegression(C = 1/(λ·45))`), with bundled
  coefficients that the internal `TrainPreMiRnaClassifier` re-derives from the repository (a unit test
  checks it). Positives: 13 public-domain miRBase precursors [5]; negatives: 4 Altschul-Erickson
  dinucleotide shuffles [14][15] per positive (Random(20060101)) — the standard pre-miRNA-classifier
  background convention [6][8]. Held-out (last 20 of 65): accuracy 0.90, AUC 0.984375. Default
  threshold 0.5 = minimum-error cutoff under the training prior (1 : 4).
  70/30 split. **No GPL miRDeep2 code or weights are used** (its published method was consulted only).

**Intentionally simplified:**

- Stem discovery requires uninterrupted end-pairing from the candidate ends inward; **consequence:** biologically valid precursors with internal bulges, asymmetric loops, or local mismatches are frequently missed.
- Mature and star sequences are extracted symmetrically from fixed arm positions; **consequence:** arm switching, offset mature products, and heterogeneous isomiR processing are not represented.
- Candidate discovery is exhaustive over fixed windows rather than full secondary-structure optimization; **consequence:** overlapping windows can generate multiple related candidates and the search remains purely heuristic.

**Not implemented:**

- The **read-stacking** (small-RNA-seq pileup) signal of miRDeep2 — the probabilistic score component
  that requires an alignment of the caller's sequencing reads to the candidate locus (read-end
  consistency with Drosha/Dicer processing). This needs the caller's reads and cannot be derived from
  sequence/structure alone; **users should rely on:** miRDeep2 with their own small-RNA-seq reads for
  that read-supported signal. (The structure/sequence-feature classifier is now bundled via
  `ClassifyPreMiRna`; bulge-tolerant folding via the opt-in MFE-fold path; cleavage coordinates via
  `PredictDroshaDicerCleavage`.) Pseudoknotted precursors remain out of scope (energy-model floor,
  RNA-STRUCT-001).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Mature-strand extraction is fixed to the 5' arm. | Assumption | Precursors whose dominant mature product derives from the opposite arm are represented asymmetrically. | accepted | See ASM-03. |
| 2 | Uninterrupted stem pairing (screen mode only) is stricter than natural pre-miRNA structure. | Deviation | Real miRBase precursors with bulges are rejected by the screen; the default MFE mode detects them. | resolved by default (F40) | See ASM-01. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `null` input | Returns no candidates. | The public method exits early for null or empty input. |
| DNA input | `T` is normalized to `U` before scanning. | The method treats DNA input as RNA-like sequence content. |
| Input shorter than `minHairpinLength` | Returns no candidates. | No candidate window can be formed. |
| Candidate with stem `< 18` | Rejected. | The validator requires a precursor-scale stem. |
| Candidate with loop `< 3` or `> 25` | Rejected. | The validator enforces the current loop-size bounds. |
| Perfect synthetic hairpin with uninterrupted complementarity | Accepted if other length constraints are satisfied. | It matches the implementation's strongest-case topology. |

### 6.2 Limitations

The opt-in screen mode of `FindPreMiRnaHairpins` (`useConsecutiveScreen: true`; the default before F40) is deliberately stricter than natural pre-miRNA folding: it
does not tolerate internal bulges, asymmetric loops, pseudoknots, or arm ambiguity, and therefore
can miss bona fide miRBase precursors that do not exhibit uninterrupted end-to-end pairing. The
**default** MFE mode (`FindPreMiRnaHairpinsByMfe` / `AssessHairpinByMfe`) removes that bulge-intolerance by
folding each candidate with the validated Zuker–Stiegler MFE engine and reading the hairpin from the
real MFE structure (detecting e.g. hsa-mir-21 and hsa-let-7a-1). The opt-in
`PredictDroshaDicerCleavage` adds Drosha/Dicer cleavage-site (mature/star excision-coordinate)
prediction from the published measuring rules (Han 2006 / Park 2011). The opt-in `ClassifyPreMiRna`
adds a **trained** structure/sequence-feature natural-vs-background classifier (logistic regression
over MFE/AMFE/MFEI/GC/%paired, trained from public-domain miRBase positives vs di-shuffled negatives;
held-out AUC = 1.0). **Residual scope:** only the **read-stacking** (small-RNA-seq pileup) signal of
miRDeep2 — which requires the caller's sequencing reads — remains out of scope; no GPL miRDeep2 code
is used.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
using Seqeron.Genomics.Annotation;

string sequence =
	"GCAUAGCUAGCUAGCUAGCUAGCUA" +
	"GAAAUUU" +
	"UAGCUAGCUAGCUAGCUAGCUAUGC";

var candidates = MiRnaAnalyzer.FindPreMiRnaHairpins(sequence).ToList();

// candidates[0].Structure contains a simple dot-bracket hairpin
// candidates[0].MatureSequence.Length <= 22

// Opt-in: assess a candidate from its REAL MFE structure (RNA-STRUCT-001 folder).
var mfe = MiRnaAnalyzer.AssessHairpinByMfe(sequence);
// mfe.FreeEnergy == RnaSecondaryStructure.CalculateMinimumFreeEnergy(sequence) == -48.48
// mfe.StemBasePairs == 27, mfe.TerminalLoopSize == 3, mfe.Mfei ≈ 1.9392 (≥ 0.85)

// hsa-mir-21 (MI0000077): heuristic returns no candidates, but the MFE fold detects it:
var detected = MiRnaAnalyzer.AssessHairpinByMfe(
    "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA");
// detected.FreeEnergy == -35.13, detected.StemBasePairs == 32, detected.Mfei ≈ 1.0037

// Opt-in: predict Drosha/Dicer cleavage coordinates from the published measuring rules.
// pri-miR-21 = miRBase MI0000077 with 30 nt of hg38 flank each side; RNAfold 2.8 structure; junction 23.
var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(
    "CAUCGUGACAUCUCCAUGGCUGUACCACCUUGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACAUUUUGGUAUCUUUCAUCUGACCAUCCAUAU", 23,
    "..(((((.......)))))..((((((...((((((((((((((((.(((((.(((((.((((.((...))))))))))).)))))))))))))))))))))...)))))).....................")!.Value;
// cut.DroshaCut5Prime == 37 (5' base of the 12th stem pair, Han 2006); cut.DroshaCut3Prime == 96 (partner 94 + 2)
// cut.MatureSequence == "UAGCUUAUCAGACUGAUGUUGA" (22 nt, Park 2011; miRBase hsa-miR-21-5p)
// cut.StarStart..StarEnd == 76..96 (3' arm; miRBase hsa-miR-21-3p = 75..97)

// Opt-in: trained structure/sequence-feature natural-vs-background classifier.
var call = MiRnaAnalyzer.ClassifyPreMiRna(
    "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA")!.Value;
// call.IsNatural == true, call.NaturalProbability ≈ 0.99999 (real hsa-mir-21 precursor)

// A dinucleotide-shuffled (composition-matched) version scores as background:
var bg = MiRnaAnalyzer.ClassifyPreMiRna(
    MiRnaAnalyzer.DinucleotideShuffle(
        "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA",
        new Random(999)))!.Value;
// bg.IsNatural == false, bg.NaturalProbability ≈ 0.0005
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MiRnaAnalyzer_PreMiRna_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/MiRnaAnalyzer_PreMiRna_Tests.cs)
- Test spec: [MIRNA-PRECURSOR-001.md](../../../tests/TestSpecs/MIRNA-PRECURSOR-001.md)
- Evidence: [MIRNA-PRECURSOR-001-Evidence.md](../../../docs/Evidence/MIRNA-PRECURSOR-001-Evidence.md)
- Related property tests: [MiRnaProperties.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Properties/MiRnaProperties.cs)
- Related snapshots: [MiRnaSnapshotTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Snapshots/MiRnaSnapshotTests.cs)

## 8. References

1. Bartel DP. 2004. MicroRNAs: genomics, biogenesis, mechanism, and function. Cell. 116(2):281-297.
2. Ambros V, Bartel B, Bartel DP, Burge CB, Carrington JC, Chen X, Dreyfuss G, Eddy SR, Griffiths-Jones S, Marshall M, et al. 2003. A uniform system for microRNA annotation. RNA. 9(3):277-279.
3. Turner DH, Mathews DH. 2010. NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure. Nucleic Acids Research. 38(Database issue):D280-D282.
4. Bartel DP. 2009. MicroRNAs: target recognition and regulatory functions. Cell. 136(2):215-233.
5. Kozomara A, Birgaoanu M, Griffiths-Jones S. 2019. miRBase: from microRNA sequences to function. Nucleic Acids Research. 47(D1):D155-D162.
6. Bonnet E, Wuyts J, Rouzé P, Van de Peer Y. 2004. Evidence that microRNA precursors, unlike other non-coding RNAs, have lower folding free energies than random sequences. Bioinformatics. 20(17):2911-2917. doi:10.1093/bioinformatics/bth374.
7. Zhang BH, Pan XP, Cox SB, Cobb GP, Anderson TA. 2006. Evidence that miRNAs are different from other RNAs. Cellular and Molecular Life Sciences. 63(2):246-254. (AMFE = 100·MFE/length; MFEI = AMFE/(G+C)%; pre-miRNA MFEI > 0.85.)
8. Meyers BC, Axtell MJ, Bartel B, Bartel DP, Baulcombe D, Bowman JL, et al. 2008. Criteria for annotation of plant MicroRNAs. The Plant Cell. 20(12):3186-3190. doi:10.1105/tpc.108.064311.
9. Han J, Lee Y, Yeom KH, Nam JW, Heo I, Rhee JK, Sohn SY, Cho Y, Zhang BT, Kim VN. 2006. Molecular basis for the recognition of primary microRNAs by the Drosha-DGCR8 complex. Cell. 125(5):887-901. doi:10.1016/j.cell.2006.03.043. ("The cleavage site is determined mainly by the distance (approximately 11 bp) from the stem-ssRNA junction.")
10. Park JE, Heo I, Tian Y, Simanshu DK, Chang H, Jee D, Patel DJ, Kim VN. 2011. Dicer recognizes the 5' end of RNA for efficient and accurate processing. Nature. 475(7355):201-205. doi:10.1038/nature10198. ("the cleavage site determined mainly by the distance (∼22 nucleotides) from the 5' end (5' counting rule).")
11. Lee Y, Ahn C, Han J, Choi H, Kim J, Yim J, Lee J, Provost P, Rådmark O, Kim S, Kim VN. 2003. The nuclear RNase III Drosha initiates microRNA processing. Nature. 425(6956):415-419. doi:10.1038/nature01957. (RNase III staggered cleavage leaves a 2-nt 3' overhang.)
12. Auyeung VC, Ulitsky I, McGeary SE, Bartel DP. 2013. Beyond secondary structure: primary-sequence determinants license pri-miRNA hairpins for processing. Cell. 152(4):844-858. doi:10.1016/j.cell.2013.01.031. (Basal UG, apical UGU(G), and CNNC motifs; CNNC positioned 16-18 nt from the Drosha cut.)
13. Hastie T, Tibshirani R, Friedman J. 2009. The Elements of Statistical Learning, 2nd ed., §4.4.1 (logistic regression; sigmoid link; log-likelihood maximised by Newton–Raphson / IRLS, eqs. 4.21–4.28). doi:10.1007/978-0-387-84858-7. Also: Batuwita R, Palade V. 2009. microPred: effective classification of pre-miRNAs for human miRNA gene prediction. Bioinformatics. 25(8):989-995. doi:10.1093/bioinformatics/btp107 (MFE/MFEI + base-pair-composition / %paired structural feature family) and Xue C, Li F, He T, Liu GP, Li Y, Zhang X. 2005. Classification of real and pseudo microRNA precursors using local structure-sequence features and support vector machine. BMC Bioinformatics. 6:310. doi:10.1186/1471-2105-6-310.
14. Altschul SF, Erickson BW. 1985. Significance of nucleotide sequence alignments: a method for random sequence permutation that preserves dinucleotide and codon usage. Molecular Biology and Evolution. 2(6):526-538. (Eulerian-walk dinucleotide-preserving shuffle; implemented as the background-set generator. See also Jiang M, et al. 2008. uShuffle. BMC Bioinformatics. 9:192. doi:10.1186/1471-2105-9-192.)
15. Jiang M, Anderson J, Gillespie J, Mayne M. 2008. uShuffle: a useful tool for shuffling biological sequences while preserving the k-let counts. BMC Bioinformatics. 9:192. (source `ushuffle.c` 1.1.2, `shuffle2`: Wilson random arborescence.)
16. Ng KLS, Mishra SK. 2007. De novo SVM classification of precursor microRNAs from genomic pseudo hairpins using global and intrinsic folding measures. Bioinformatics. 23(11):1321-1330. (dP = |bp|/L, Schultes et al. 1999.)
17. Friedländer MR, Chen W, Adamidi C, Maaskola J, Einspanier R, Knespel S, Rajewsky N. 2008. Discovering microRNAs from deep sequencing data using miRDeep. Nat Biotechnol. 26(4):407-415. doi:10.1038/nbt1394. (`find_star` in rajewsky-lab/mirdeep2 `src/miRDeep2_core_algorithm.pl`: star from mature via base-pair partners with a 2-nt 3' overhang.)
