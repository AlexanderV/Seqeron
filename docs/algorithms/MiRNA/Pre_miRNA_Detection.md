# Pre-miRNA Hairpin Detection

| Field | Value |
|-------|-------|
| Algorithm Group | MiRNA |
| Test Unit ID | MIRNA-PRECURSOR-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified (default heuristic) + Production (opt-in MFE fold) |
| Last Reviewed | 2026-06-24 |

## 1. Overview

Pre-miRNA detection scans RNA sequence windows for stem-loop candidates that satisfy basic precursor-miRNA length, stem, and loop constraints [1][2][4]. The repository implements this in `MiRnaAnalyzer.FindPreMiRnaHairpins` by enumerating candidate windows, requiring uninterrupted complementary pairing from both ends inward, extracting a 5' mature strand and 3' star strand, and assigning a Turner-parameter-based hairpin free energy. The search is deterministic and exhaustive within the caller's length bounds, but it is intentionally stricter than real pre-miRNA structure because internal bulges, asymmetric loops, and arm switching are not modeled. This makes the implementation suitable for exact educational hairpin screening rather than production-grade precursor discovery on natural miRBase sequences [1][4][5].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Pre-miRNAs are stem-loop RNA intermediates in miRNA biogenesis, typically on the order of roughly 60-120 nt, from which Dicer processes a mature miRNA duplex of about 22 nt [1][2]. A valid precursor must fold into a hairpin with a stem and terminal loop, and the mature miRNA resides on one arm of that structure [1][2]. Thermodynamic stability of RNA hairpins is commonly modeled with nearest-neighbor parameters, including base-stacking energies, loop penalties, and terminal mismatch terms [3].

### 2.2 Core Model

For a candidate RNA subsequence $S$ of length $n$, the repository accepts a precursor candidate only when it can be decomposed into:

$$
S = \text{5' stem} + \text{loop} + \text{3' stem}
$$

with uninterrupted pairing between mirrored positions from the two ends inward under the pairing set:

$$
\{A\!-\!U,\; U\!-\!A,\; G\!-\!C,\; C\!-\!G,\; G\!-\!U,\; U\!-\!G\}
$$

The current acceptance criteria are:

- stem length $\ge 18$ bp,
- loop length $= n - 2 \times \text{stem length}$ with $3 \le \text{loop length} \le 25$,
- candidate length constrained by the caller, defaulting to 55-120 nt.

Hairpin energy is then estimated with Turner-style nearest-neighbor terms:

$$
\Delta G = \sum \Delta G_{\text{stack}} + \Delta G_{\text{loop}} + \Delta G_{\text{terminal mismatch}} + 0.45 \times n_{\text{terminal AU/GU}}
$$

using the repository's Turner 2004 lookup tables [3].

#### 2.2.1 Opt-in MFE-structure-based assessment (default unchanged)

`AssessHairpinByMfe` / `FindPreMiRnaHairpinsByMfe` replace the consecutive-pairing scan with the
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

#### 2.2.2 Opt-in Drosha/Dicer cleavage-site prediction (default unchanged)

`PredictDroshaDicerCleavage(sequence, basalJunction)` predicts the excision coordinates of a
pri-/pre-miRNA hairpin from the PUBLISHED measuring ("ruler") rules — it does **not** use a trained
classifier:

- **Drosha (basal-junction ruler):** Drosha cleaves $\sim 11$ bp ($\approx$ one helical turn) from the
  basal ssRNA–dsRNA junction [9]. The 5' cut is the 5' end of the 5p mature:
  $\text{DroshaCut5'} = \text{basalJunction} + 11$.
- **Dicer (5'-counting ruler):** Dicer cleaves $\sim 22$ nt from the Drosha-generated 5' end (the 5'
  counting rule) [10], fixing the mature length at $22$ nt:
  $\text{mature} = [\text{DroshaCut5'},\ \text{DroshaCut5'} + 21]$.
- **2-nt 3' overhang:** each RNase III cut (Drosha, Dicer) leaves a 2-nt 3' overhang [9][11]; the 3p
  (miRNA\*) span is the same $\sim 22$ nt with its Drosha-generated 3' end 2 nt 3' of the 5p mature end.
- **CNNC motif (optional confidence):** a C-N-N-C 16–18 nt 3' of the Drosha cut [12] sets
  `HasCnncMotif`; it is reported, not required.

Cross-checked against miRBase hsa-miR-21-5p (MIMAT0000076, `UAGCUUAUCAGACUGAUGUUGA`): feeding a
pri-miRNA whose 11-bp lower stem places the $+11$ Drosha cut at the annotated 5p start reproduces the
miRBase mature exactly.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Candidate stems are approximated by uninterrupted end-to-end pairing with no internal bulges or mismatches. | Natural pre-miRNAs with internal loops or asymmetric bulges can be missed even if they are biologically valid. |
| ASM-02 | Turner nearest-neighbor parameters are used as a stability proxy for the accepted hairpin topology. | Energies are meaningful only to the extent that the simplified topology resembles a real RNA hairpin at standard conditions [3]. |
| ASM-03 | The mature product is extracted from the 5' arm and the star sequence from the mirrored 3' arm. | Valid precursors whose dominant mature product comes from the opposite arm are represented asymmetrically. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported precursor has a balanced dot-bracket structure string of the same length as `Sequence`. | `AnalyzeHairpin` fills the first `stemLength` positions with `(`, the last `stemLength` with `)`, and the loop with `.`. |
| INV-02 | Every reported precursor has stem length `>= 18`. | Candidates with shorter uninterrupted stems are rejected before record creation. |
| INV-03 | Every reported precursor has loop length between `3` and `25` nt. | `AnalyzeHairpin` explicitly rejects loops outside that range. |
| INV-04 | `Start` and `End` are zero-based inclusive coordinates in the scanned input sequence. | `FindPreMiRnaHairpins` constructs `End = Start + length - 1`. |
| INV-05 | `MatureSequence.Length = StarSequence.Length = min(matureLength, stemLength)`. | Both arms are extracted using the same bounded length. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | RNA or DNA sequence to scan for precursor candidates. | `null`, empty, or shorter-than-`minHairpinLength` input yields no results; `T` is normalized to `U`. |
| `minHairpinLength` | `int` | `55` | Minimum candidate window length considered during scanning. | Not explicitly validated; values larger than the input length produce no results. |
| `maxHairpinLength` | `int` | `120` | Maximum candidate window length considered during scanning. | Not explicitly validated; values below `minHairpinLength` make the inner scan empty. |
| `matureLength` | `int` | `22` | Maximum length extracted for mature and star strands. | Expected to be positive; explicit validation is not performed. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Start` | `int` | Zero-based inclusive start coordinate of the candidate in the original input sequence. |
| `End` | `int` | Zero-based inclusive end coordinate of the candidate in the original input sequence. |
| `Sequence` | `string` | Uppercase RNA subsequence that satisfied the precursor tests. |
| `MatureSequence` | `string` | Extracted mature strand from the 5' arm. |
| `StarSequence` | `string` | Extracted star strand from the mirrored 3' arm. |
| `Structure` | `string` | Dot-bracket representation with contiguous stem parentheses and loop dots. |
| `FreeEnergy` | `double` | Turner-table-based heuristic energy for the accepted hairpin topology. |

### 3.3 Preconditions and Validation

`FindPreMiRnaHairpins` is non-throwing for `null`, empty, or too-short input: those cases simply yield no candidates. Input is normalized to uppercase RNA by replacing `T` with `U`. Candidate windows are rejected if they are shorter than 55 nt inside `AnalyzeHairpin`, if uninterrupted complementary pairing from the ends inward yields a stem shorter than 18 bp, or if the resulting loop falls outside 3-25 nt. Coordinates are zero-based and inclusive. Length parameters are not explicitly validated beyond their effect on the scan loops.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input to uppercase RNA.
2. Slide a start position across the sequence and enumerate every candidate window whose length is between `minHairpinLength` and `maxHairpinLength`.
3. For each candidate window, count uninterrupted complementary base pairs from the two ends inward, allowing Watson-Crick and G:U wobble pairs.
4. Reject the candidate unless the resulting uninterrupted stem is at least 18 bp and the remaining loop is between 3 and 25 nt.
5. Extract the mature strand from the 5' arm and the star strand from the mirrored 3' arm.
6. Build a dot-bracket structure and compute free energy from the Turner lookup tables.
7. Emit a `PreMiRna` record with genomic coordinates relative to the scanned input sequence.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Hairpin stability uses repository-resident Turner 2004 tables for:

- nearest-neighbor stacking energies,
- hairpin loop initiation energies,
- terminal mismatch energies,
- a fixed `0.45 kcal/mol` AU/GU terminal penalty.

The pairing test accepts both Watson-Crick and G:U wobble pairs, but the uninterrupted stem count stops at the first non-pairing mirrored position.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Exhaustive precursor scan | `O(nH^2)` | `O(H)` | `n` is input length and `H` is the maximum scanned hairpin length, because up to `O(nH)` windows are enumerated and each candidate can cost `O(H)` to analyze. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MiRnaAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs)

- `MiRnaAnalyzer.FindPreMiRnaHairpins(string, int, int, int)`: Public precursor scan over all candidate windows (default consecutive-pairing heuristic).
- `MiRnaAnalyzer.AnalyzeHairpin(string, int)`: Private validator that checks stem continuity, loop size, arm extraction, and structure generation.
- `MiRnaAnalyzer.FindPreMiRnaHairpinsByMfe(string, int, int, double, int)`: Opt-in precursor scan that folds each candidate window with the RNA-STRUCT-001 MFE engine.
- `MiRnaAnalyzer.AssessHairpinByMfe(string, double, int)`: Opt-in single-candidate assessment from the real MFE structure (single hairpin + stem≥16 + loop 3–25 + MFEI≥0.85).
- `MiRnaAnalyzer.CalculateMfeIndex(double, int, double)`: MFEI = AMFE/(G+C)%, AMFE = 100·|ΔG°|/length (Zhang 2006).
- `MiRnaAnalyzer.PredictDroshaDicerCleavage(string, int)`: Opt-in cleavage-site prediction — Drosha cut (~11 bp from the basal junction), Dicer cut / 22-nt mature length, mature (5p) + star (3p) spans, 2-nt 3' overhang, optional CNNC flag (Han 2006 / Park 2011 / Auyeung 2013).

### 5.2 Current Behavior

The implementation searches candidates exhaustively rather than folding the full transcript once. Stem detection is based on uninterrupted end-to-end pairing, so the counter stops at the first mirrored mismatch and does not resume across bulges or internal loops. `MatureSequence` is always taken from the 5' arm beginning at the candidate start, and `StarSequence` is taken from the mirrored 3' arm of equal extracted length. Hairpin energy is not the older linear approximation described by the legacy document; it is computed from Turner 2004 stacking, loop, terminal mismatch, and AU/GU penalty tables stored directly in `MiRnaAnalyzer`. Although the loop-energy table contains extrapolation logic for loops larger than 30 nt, accepted precursors currently never reach that branch because `AnalyzeHairpin` rejects loops above 25 nt.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- RNA hairpin validation using Watson-Crick and G:U wobble pairing [1][4].
- Precursor-style constraints on hairpin length, stem length, and loop size [1][2][4].
- Turner-style nearest-neighbor energy components from repository lookup tables [3].
- **Opt-in:** real MFE-structure-based hairpin assessment that folds the candidate with the
  validated Zuker–Stiegler engine (RNA-STRUCT-001) [3][6] and applies the Ambros (2003) ≥16-bp
  stem [2], Bartel (2004) 3–25-nt loop [1], and Zhang (2006) MFEI ≥ 0.85 [7] criteria to the actual
  MFE structure — tolerant of internal bulges/loops and detecting natural miRBase precursors.
- **Opt-in:** Drosha/Dicer cleavage-site prediction (`PredictDroshaDicerCleavage`) using the published
  measuring rules verbatim — Drosha ~11 bp from the basal junction [9], Dicer ~22 nt 5'-counting [10],
  the RNase III 2-nt 3' overhang [9][11], and the Auyeung (2013) CNNC confidence motif [12].

**Intentionally simplified:**

- Stem discovery requires uninterrupted end-pairing from the candidate ends inward; **consequence:** biologically valid precursors with internal bulges, asymmetric loops, or local mismatches are frequently missed.
- Mature and star sequences are extracted symmetrically from fixed arm positions; **consequence:** arm switching, offset mature products, and heterogeneous isomiR processing are not represented.
- Candidate discovery is exhaustive over fixed windows rather than full secondary-structure optimization; **consequence:** overlapping windows can generate multiple related candidates and the search remains purely heuristic.

**Not implemented:**

- A competitive **trained** natural-vs-background pre-miRNA classifier (a fitted probabilistic model
  such as miRDeep2 that scores genuine precursors against genomic background hairpins using
  read-stacking signatures) and pseudoknotted precursors; **users should rely on:** miRDeep2 / miRBase
  tooling for decision-grade, model-based precursor discovery. (Bulge-tolerant folding is covered by
  the opt-in MFE-fold path, and cleavage-site coordinates by the opt-in
  `PredictDroshaDicerCleavage` measuring-rule path above; only the trained classifier remains.)

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Mature-strand extraction is fixed to the 5' arm. | Assumption | Precursors whose dominant mature product derives from the opposite arm are represented asymmetrically. | accepted | See ASM-03. |
| 2 | Uninterrupted stem pairing is stricter than natural pre-miRNA structure. | Deviation | Real miRBase precursors with bulges or internal mismatches can be rejected. | accepted | See ASM-01. |

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

The **default** `FindPreMiRnaHairpins` is deliberately stricter than natural pre-miRNA folding: it
does not tolerate internal bulges, asymmetric loops, pseudoknots, or arm ambiguity, and therefore
can miss bona fide miRBase precursors that do not exhibit uninterrupted end-to-end pairing. The
**opt-in** `FindPreMiRnaHairpinsByMfe` / `AssessHairpinByMfe` remove that bulge-intolerance by
folding each candidate with the validated Zuker–Stiegler MFE engine and reading the hairpin from the
real MFE structure (detecting e.g. hsa-mir-21 and hsa-let-7a-1). The opt-in
`PredictDroshaDicerCleavage` adds Drosha/Dicer cleavage-site (mature/star excision-coordinate)
prediction from the published measuring rules (Han 2006 / Park 2011). **Residual scope:** only a
competitive **trained** natural-vs-background precursor classifier (e.g. miRDeep2 — a fitted
probabilistic model, data-blocked) remains out of scope.

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
// Pri-miRNA = 11-nt lower stem + miR-21 stem region; junction at index 0.
var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(
    "CCCCCCCCCCC" + "UAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGU", 0)!.Value;
// cut.DroshaCut5Prime == 11 (junction + 11 bp, Han 2006)
// cut.MatureSequence == "UAGCUUAUCAGACUGAUGUUGA" (22 nt; == miRBase hsa-miR-21-5p)
// cut.ThreePrimeOverhang == 2 (RNase III 2-nt 3' overhang)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MiRnaAnalyzer_PreMiRna_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs)
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
