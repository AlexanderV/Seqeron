# Probe Validation

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PROBE-VALID-001 |
| Related Projects | N/A |
| Implementation Status | Complete (Primer3 thermodynamic self-structure screen; Kane et al. 2000 cross-hybridization criteria) |
| Last Reviewed | 2026-10-08 |

## 1. Overview

Probe validation assesses whether a hybridization probe is likely to bind specifically to its intended targets with limited cross-hybridization. In this repository, `ValidateProbe` combines (1) substitution-tolerant fixed-length window matching against reference sequences, (2) Primer3's thermodynamic hybridization-probe self-structure screen — ntthal self-dimer, 3′ self-dimer and hairpin Tm ≤ 47 °C, the same screen `DesignProbes` uses — and (3) optionally the Kane et al. (2000) cross-hybridization criteria against known non-target sequences (`AssessCrossHybridization`: overall identity > 75 % or a contiguous identical stretch > 15 nt) [7]. The result is a validation record with a uniqueness score, the ntthal Tm values, per-non-target assessments, an issue list and `IsValid` = no issue.

An opt-in **gapped** off-target scan (`ScanOffTargetsGapped`) supplements the default ungapped scan: it reuses the library's validated Smith–Waterman local aligner [1] to find off-target sites reachable through insertions or deletions — the indel-aware "BLAST-grade" improvement over a pure ungapped Hamming scan [2] — and it separates the single intended on-target match from genuine off-target hits (correcting the on/off-target pooling of `ValidateProbe`'s `OffTargetHits`). An opt-in **Karlin–Altschul statistics** layer (`ComputeKarlinAltschul` / `ComputeLambdaNucleotide`) quantifies the statistical significance of an off-target hit's raw alignment score as a bit score and an E-value [8][9]. `CheckSpecificity` counts exact suffix-tree occurrences, optionally on both strands; `ValidateProbe`'s ungapped reference scan likewise searches the given strand by default and, opt-in (`bothStrands`), the probe's reverse complement too.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Cross-hybridization occurs when a probe binds sequences other than its intended target, and this is a central design concern for FISH, DNA microarrays, qPCR, and related assays. The original document also notes that mismatch tolerance, assay stringency, and probe length affect specificity, and that self-complementarity and low-complexity content can increase the risk of non-specific behavior. Sources: Wikipedia (Hybridization probe, DNA microarray, BLAST), Altschul et al. (1990), Amann & Ludwig (2000).

Off-target detection by *local alignment* follows the Smith–Waterman recurrence [1]: `H(i,j) = max{ H(i-1,j-1) + s(a_i,b_j), H(i-1,j) - W, H(i,j-1) - W, 0 }`, whose zero floor returns the best-scoring local subsequence match and whose gap terms admit insertions/deletions. The gapped local alignment is the recognized improvement over ungapped matching for similarity/homology search [2]. **Kane et al. (2000) criteria [7].** Kane et al. hybridized 50-mer oligonucleotide arrays with non-target cDNAs of graded similarity: "any 'non-target' transcripts (cDNAs) >75% similar over the 50 base target may show cross-hybridization", and a non-target must not share "a stretch of complementary sequence >15 contiguous bases" (abstract; later probe-design pipelines summarise the rule as "overall sequence identity > 75 % or a contiguous match > 15 bp" [10]). `AssessCrossHybridization` implements both criteria per non-target strand: identity = identical columns of the best local alignment (Smith–Waterman–Gotoh with BLAST+ blastn scoring +2/−3, gap existence 5 / extension 2 — the canonical `SequenceAligner.LocalAlignAffine`) ÷ probe length, flagged when strictly > `maxIdentity` (0.75); longest contiguous match = longest common substring (canonical suffix tree), flagged when strictly > `maxContiguousMatch` (15). Both strands are assessed by default (a double-stranded non-target offers both). `ScanOffTargetsGapped` uses the 0.75 value as its site-reporting threshold (identity ≥ `minIdentity`).

**Self-structure (Primer3) [11].** Primer3 screens a hybridization probe (internal oligo) with ntthal: self-dimer (ANY) Tm ≤ PRIMER_INTERNAL_MAX_SELF_ANY_TH, 3′ self-dimer (END1) Tm ≤ PRIMER_INTERNAL_MAX_SELF_END_TH, hairpin Tm ≤ PRIMER_INTERNAL_MAX_HAIRPIN_TH, all 47 °C by default, at the internal-oligo conditions (50 nM, 50 mM monovalent, 0 Mg²⁺, 0 dNTP). `ValidateProbe` applies this screen (via `PrimerDesigner.CalculatePrimer3OligoStructure`, bit-exact to primer3-py) to A/C/G/T probes of ≤ `ProbeParameters.ThermodynamicScreenMaxLength` nt (thal.h `THAL_MAX_ALIGN`, default 60 as in Primer3; a larger value is an opt-in that runs the unchanged ntthal recursions on longer probes — thal.c compiled with `-DTHAL_MAX_ALIGN=…` parity, audit round 3 A3-9 / F55); longer probes and non-ACGT probes fall back to Primer3's alignment-mode internal-oligo self-dimer screen (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0, `oligo_compl`: dpal `self_any` / `self_end` > PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END = 12.00, `ProbeParameters.MaxSelfAny/MaxSelfEnd`; no length limit; bit-exact to compiled dpal.c) and the sequence-only stem-loop screen `PrimerDesigner.HasHairpinPotential` (an exactly complementary stem of ≥ 4 bp — library convention, F49 — closing a loop of ≥ 3 nt, thal.c `min_hrpn_loop` = 3; audit round 5, A5-2: formerly a private scan of 3-nt loops only with an unsourced ≥ 80 % stem match, which missed a 74-nt probe's 10-bp stem with a 4-nt loop). It flags almost every probe of > 60 nt; raise `ThermodynamicScreenMaxLength` for a hairpin-Tm decision. `ProbeValidation.SelfAny` / `SelfEnd` report the alignment-mode values for every probe; `SelfComplementarity` (position-wise fold-back fraction) is an informational library metric only.

### 2.2 Core Model

The validation workflow normalizes the probe to uppercase, counts approximate matches across all supplied reference sequences, and judges every hit by the **Kane et al. (2000) criteria on its ungapped diagonal** [7]: identity = (L − mismatches) / L over the probe length > `maxNonTargetIdentity` (0.75), or a run of identical positions > `maxContiguousMatch` (15 nt) — the same two measures `AssessCrossHybridization` takes from the best local alignment / longest common substring. More than one such site (`CrossHybridizingHits` > 1) records the off-target issue: the probe can bind somewhere besides its intended site. With the default radius of 3 mismatches every hit of a probe ≥ 13 nt meets the identity criterion ((L − 3) / L > 0.75 ⇔ L > 12), so there the decision equals the former "more than one hit" rule; for shorter probes (or a larger radius) a site at ≤ 75 % identity without a > 15-nt run is no longer an issue (audit round 3, A3-12).

**Strands (opt-in `bothStrands`, audit round 4, A4-2).** A double-stranded reference offers the probe a binding site wherever its reverse complement occurs: BLAST+ blastn searches both strands by default (`-strand both`, `blastn -help`: "Query strand(s) to search against database/subject. Default = `both'"), and OligoArray 2.0 checks each candidate with BLAST against the genome [12]. With `bothStrands = true` each reference is scanned for the probe and for its reverse complement (canonical `DnaSequence.GetReverseComplementString`, the same canonical `ApproximateMatcher.FindWithMismatches` scan); the sites are the **union of the hit positions** per reference, and a position hit in both orientations — every hit of a reverse-palindromic probe, or a reverse-palindromic site (Hamming(S, P) = Hamming(rc S, rc P) = Hamming(S, rc P) when S = rc S) — is one site, cross-hybridizing when either orientation meets a Kane criterion. For exact matching this is precisely `CheckSpecificity(…, bothStrands: true)` (a probe and its reverse complement can start at the same position only when the probe is a reverse palindrome, which that method counts once). The default (`false`) keeps the given-strand scan (backward compatible); the genome-index `DesignProbes` overload has the same opt-in, passed to `CheckSpecificity`.

**Library conventions (no published source).** The search radius `maxMismatches = 3` (the figure is the CRISPR/Cas9 guide mismatch tolerance, not a hybridization criterion; to reach every ungapped site above 75 % identity pass ⌈L/4⌉ − 1) and the hit count mapped to a **library-defined uniqueness score** (the share of the probe's N candidate binding sites taken by one site). No published specificity metric of the form 1/N was found: Kane et al. (2000) judge a non-target by identity / contiguous identity [7], OligoArray 2.0 by the duplex stability of each BLAST hit [12], Li & Stormo (2001) by sequence uniqueness and hybridization free energy (abstract only), Primer3 by the dpal alignment score against a mishybridization library (`o_args.max_repeat_compl` = PRIMER_INTERNAL_MAX_LIBRARY_MISHYB = 12.00, `libprimer3.cc`) [11]. The score is therefore reported only; it does not enter `IsValid`:

$$
specificity =
\begin{cases}
0.0, & offTargetHits = 0 \\
1.0, & offTargetHits = 1 \\
1.0 / offTargetHits, & offTargetHits > 1
\end{cases}
$$

It then applies the Primer3 thermodynamic self-structure screen (fallback: Primer3 alignment-mode self_any / self_end > 12.00 and a sequence-level inverted-repeat screen), the optional Kane assessment of the non-targets, and sets `IsValid` = no issue recorded (every issue rests on a sourced criterion: Kane sites, Primer3 limits, Kane non-targets, optional OligoArray duplex Tm). The fallback inverted-repeat hairpin flag (stem ≥ 4 bp, a library convention) is reported in `HasSecondaryStructure` / `Warnings` but is not an issue (audit round 5, A5-5): it flags 98–100 % of random > 60-nt probes, so it made almost every FISH / Northern / Southern probe invalid.

**Karlin–Altschul statistics of an off-target hit (opt-in).** For a hit's raw alignment score `S` against a search space of query length `m` and database length `n`, the statistical significance follows the Karlin–Altschul framework [8][9]:

$$E = K \cdot m \cdot n \cdot e^{-\lambda S}, \qquad S' = \frac{\lambda S - \ln K}{\ln 2}, \qquad E = m \cdot n \cdot 2^{-S'}$$

where the scale parameter `λ` is the unique positive root of the defining equation [8]

$$\sum_{i,j} p_i\, p_j\, e^{\lambda s_{ij}} = 1$$

with `p_i` the background base frequencies and `s_ij` the score matrix. For four equiprobable bases (`p_i = 0.25`) and a match/mismatch scheme this reduces to `0.25·e^{λ·match} + 0.75·e^{λ·mismatch} = 1`; for the BLAST `+1/−3` scheme it solves to `λ ≈ 1.374` (matching the value NCBI blastn reports). The theory requires a scoring scheme with **negative expected per-pair score** and **at least one positive score** [9]; both are checked. `K` is computed for the same scheme and composition with the Karlin–Altschul lattice formula as NCBI BLAST+ computes it [13] (`ComputeUngappedKarlinParameters`: `BlastKarlinLtoH` `H = λ·Σ s·p_s·e^{λs}`, `BlastKarlinLHtoK` closed forms when the gcd-reduced lattice has lowest score −1 or highest +1, otherwise the convergent series over gapless score distributions, sum limit 10⁻⁴, ≤ 100 terms): `+1/−3 → K = 0.7106`, `+2/−3 → 0.4081` (blastn prints `0.711`, `0.408`); a caller-supplied `K` still overrides it. For a scheme whose scores share a divisor δ > 1 the default (`KarlinKMethod.ReducedLattice`) is scale-invariant (`K(+4/−6) = K(+2/−3)`), as the theory requires; `KarlinKMethod.NcbiBlast` reproduces BLAST+ exactly — `BlastKarlinLHtoK` reduces low/high/λ by δ but its series reads `probArrayStartLow[j]` from the unreduced score array — giving the K blastn prints (`+4/−6 → 1.17`, `+4/−10 → 1.06`, `+6/−4 → 1.63`; closed-form cases such as `+2/−2`, `+2/−4`, `+4/−2` agree under both).

**BLAST+ blastn statistics (opt-in, `ComputeBlastnStatistics`).** For a *gapped* (affine) alignment score — e.g. the `AlignmentScore` of `AssessCrossHybridization`, produced by the canonical `SequenceAligner.LocalAlignAffine` — λ and K have no closed form; BLAST+ uses simulated values tabulated per reward/penalty/gap cost (`blastn_values_*` in `blast_stat.c` [13]; `+2/−3`, gap 5/2 — the blastn-task default and `SequenceAligner.BlastDna` — → `λ 0.625, K 0.41, H 0.78, α 0.8, β −2`; gap costs at or beyond the table's "infinite" domain use the ungapped values; scores of a gcd > 1 scheme rescale λ and α by the gcd; `2/−3`, `2/−5`, `2/−7`, `3/−4` round odd scores down to even for the E-value). The search space is corrected for edge effects: the length adjustment ℓ is the integer approximation of the fixed point of `ℓ = β + (α/λ)(ln K + ln((m − ℓ)(n − Nℓ)))` (`BLAST_ComputeLengthAdjustment`; ungapped `α/λ = 1/H`), and `E = K·(m − ℓ)·max(1, n − Nℓ)·e^{−λS}` with `N` database sequences, bit score `S' = (λS − ln K)/ln 2` on the raw score. Its ungapped K (ungapped statistics, and the ungapped block copied in the infinite gap-cost domain) is BLAST+'s own (`KarlinKMethod.NcbiBlast`, default; `ReducedLattice` opt-out): blastn `+4/−6 -ungapped`, 200 × 3100 nt → effective search space 573996, score 400 → 182 bits, E 6.03e-50. The overload `ComputeBlastnStatistics(probe, subject)` aligns with `LocalAlignAffine` and evaluates the bl2seq search space (m = probe length, n = subject length, N = 1).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0.0 <= SpecificityScore <= 1.0` | The score is explicitly mapped from hit counts to `0`, `1`, or `1/hits` |
| INV-02 | `SelfComplementarity` is non-negative and at most `1.0` | It is computed as a fraction of aligned positions |
| INV-03 | `OffTargetHits >= 0` | Hit counts are accumulated from match enumeration |
| INV-04 | `OffTargetHits == 1` implies `SpecificityScore == 1.0` | That mapping is explicit in source |
| INV-05 | Each gapped hit has `0.0 <= Identity <= 1.0` and `0.0 <= Coverage <= 1.0` | Both are counts of aligned columns divided by probe length |
| INV-06 | The intended on-target (first perfect ungapped full-coverage exact match) is excluded from `OffTargetHits` | On/off separation in `ScanOffTargetsGapped` |
| INV-07 | An indel-only off-target has `HasGaps == true` and is found by `ScanOffTargetsGapped` but not by the ungapped `ValidateProbe` scan | Gapped local alignment admits indels [1][2] |
| INV-08 | `ComputeLambdaNucleotide` returns the unique positive root of `Σ p_i p_j e^{λ s_ij} = 1`; for `+1/−3`, `p=0.25` it equals `1.374` to numerical tolerance | Bisection on the strictly-crossing Karlin–Altschul equation [8] |
| INV-09 | The two E-value forms agree: `K·m·n·e^{−λS} == m·n·2^{−S'}`; E decreases as `S` increases and scales linearly with `m·n` | Algebraic identity of the Karlin–Altschul formulas [8][9] |
| INV-10 | `ComputeUngappedKarlinParameters`: `K > 0`, `H > 0`; K invariant and λ inversely scaled when all scores are multiplied by an integer | Karlin & Altschul (1990) lattice theory [8][13] |
| INV-11 | `ComputeBlastnStatistics`: `0 ≤ ℓ`, `K(m − ℓ)(n − Nℓ) > max(m, n)` whenever `ℓ > 0`; E non-increasing in S; `EValueScore` is S or S rounded down to even | `BLAST_ComputeLengthAdjustment`, `Blast_HSPListGetEvalues` [13] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `probeSequence` | `string` | required | Probe sequence to validate | Null input throws `ArgumentNullException`; empty string yields a structured invalid result |
| `referenceSequences` | `IEnumerable<string>` | required | Reference sequences scanned for approximate matches | Null input throws `ArgumentNullException` |
| `maxMismatches` | `int` | `3` | Search radius of the ungapped site scan (library convention) | ≥ 0; found sites are judged by the Kane criteria |
| `selfComplementarityThreshold` | `double` | `0.3` | Legacy fold-back-fraction limit | Kept for source compatibility; no screen uses it (fallback self-dimer limits: `ProbeParameters.MaxSelfAny/MaxSelfEnd` = 12.00) |
| `conditions` | `ProbeParameters?` | Primer3 probe conditions (the `Defaults.Microarray` settings at 50 nM / 50 mM / 0 / 0) | Salt / dNTP / oligo concentrations (ntthal screen and non-target duplex Tm), `StructureScreen`, `MaxStructureTm` (47 °C), `MaxSelfAny` / `MaxSelfEnd` (12.00), `ThermodynamicScreenMaxLength` (THAL_MAX_ALIGN of the screen and the duplex Tm, default 60) | Primer3 `_pr_data_control`: monovalent and oligo > 0, Mg²⁺ and dNTP ≥ 0; THAL_MAX_ALIGN 60–10 000; else `ArgumentOutOfRangeException` |
| `nonTargetSequences` | `IEnumerable<string>?` | `null` | Known non-targets for the Kane assessment | Optional |
| `maxNonTargetIdentity` / `maxContiguousMatch` | `double` / `int` | `0.75` / `15` | Kane thresholds (strict `>`) for the reference sites and the non-targets | [7]; identity ∈ [0, 1], contiguous ≥ 0, else `ArgumentOutOfRangeException` |
| `maxDuplexTm` | `double?` | `null` | Optional off-target site duplex-Tm threshold (°C) | [12] |
| `genomeIndex` | `ISuffixTree` | required for `CheckSpecificity(...)` | Pre-built suffix tree for exact hit counting | Used only by the suffix-tree specificity helper |
| `bothStrands` (`CheckSpecificity`) | `bool` | `false` | Also count reverse-complement occurrences (other strand of a ds genome; palindromes once) | Optional |
| `bothStrands` (`ValidateProbe`) | `bool` | `false` | Also scan the references for the probe's reverse complement; a position hit in both orientations (palindromic probe or site) is one site, Kane-flagged when either orientation is | Optional (A4-2); the non-target assessment is two-stranded regardless |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `IsValid` | `bool` | True when no issue was recorded |
| `SpecificityScore` | `double` | Library uniqueness score in `0.0-1.0` (0 or 1/N; reported only, not used by `IsValid`) |
| `OffTargetHits` | `int` | Total approximate hits across all reference sequences |
| `CrossHybridizingHits` | `int` | Hits meeting a Kane criterion on their ungapped diagonal; > 1 records the off-target issue |
| `SelfComplementarity` | `double` | Fold-back fraction (informational library metric) |
| `SelfAny` / `SelfEnd` | `double?` | Primer3 alignment-mode internal-oligo self_any / self_end (dpal), every non-empty probe; fallback self-dimer criterion (> 12.00) |
| `HasSecondaryStructure` | `bool` | ntthal hairpin Tm > `MaxStructureTm` (an issue); fallback: `PrimerDesigner.HasHairpinPotential`, exact stem ≥ 4 bp + loop ≥ 3 nt (a warning, not an issue) |
| `Issues` | `IReadOnlyList<string>` | Recorded validation issues (sourced criteria; any issue → `IsValid` false) |
| `Warnings` | `IReadOnlyList<string>` | Reported findings that do not decide `IsValid`: the fallback stem-loop flag ("Potential secondary structure formation"; A5-5) |
| `ThermodynamicScreen` | `bool` | True when the Primer3 ntthal screen applied |
| `SelfDimerTm` / `SelfEndDimerTm` / `HairpinTm` | `double?` | ntthal Tm (°C; Primer3 reports negative Tm as 0); null for the fallback |
| `CrossHybridization` | `IReadOnlyList<CrossHybridizationAssessment>` | Per non-target strand: identity, identical columns, score, longest contiguous match, the two Kane flags |

### 3.3 Preconditions and Validation

`ValidateProbe(...)` uppercases the input probe before analysis. Null probe or reference collections raise `ArgumentNullException`. An empty probe sequence returns a structured invalid result with `SpecificityScore = 0.0`, `OffTargetHits = 0`, and an `"Empty probe sequence"` issue. `CheckSpecificity(...)` uppercases the probe before querying the suffix tree.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the probe sequence to uppercase.
2. Search every reference sequence for approximate matches within the mismatch tolerance (with `bothStrands` also for the probe's reverse complement; union of the hit positions per reference).
3. Accumulate the total number of hits across all references; count the hits meeting a Kane criterion (identity (L − d)/L > 0.75 or a > 15-nt identical run); more than one → issue.
4. Self-structure: ntthal self-dimer / 3′ self-dimer / hairpin Tm vs `MaxStructureTm` (≤ 60-nt ACGT, thermodynamic screen) or the sequence-only screens.
5. Optional Kane assessment of each non-target (both strands).
6. Map hit count to the uniqueness score; `IsValid` = no issue.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Validation defaults preserved from the original document and source:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `maxMismatches` | `3` | Search radius of the ungapped site scan (library convention) |
| Reference-site decision | > 1 site with identity > 0.75 or a > 15-nt identical run | Kane et al. (2000) [7] |
| `selfComplementarityThreshold` | `0.3` | Legacy fold-back-fraction limit; kept for compatibility, unused by the screen |
| `MaxSelfAny` / `MaxSelfEnd` | `12.00` | Fallback self-dimer limits (PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END) |
| `MaxStructureTm` | `47 °C` | Primer3 PRIMER_INTERNAL_MAX_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH |
| Kane identity / contiguous | `> 0.75` / `> 15 nt` | Kane et al. (2000) [7] |
| Secondary-structure check | enabled | `ValidateProbe(...)` always checks the hairpin |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ValidateProbe` | `O(n × g × m)` | `O(1)` auxiliary | The original document describes dependence on probe length `n`, reference count `g`, and reference lengths `m` |
| `CheckSpecificity` | `O(m)` | `O(1)` | Suffix-tree exact hit counting for probe length `m` |
| `AssessCrossHybridization` | `O(Σ n_k · m)` | `O(m · min(n_k, 4096))` | Affine local alignment per non-target strand (chunked above 4096 nt with an overlap ≥ the longest positive-scoring alignment span, so the score is exact) + streaming suffix-tree LCS |
| `ScanOffTargetsGapped` | `O(g × n × m)` | `O(n × m)` per window | Sliding Smith–Waterman: for each of `g` reference start positions, a local alignment of probe length `n` against a window of width `m = n + 2`. Reuses `SequenceAligner.LocalAlign`. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProbeDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs)

- `ProbeDesigner.ValidateProbe(string, IEnumerable<string>, int, double, ProbeParameters?, IEnumerable<string>?, double, int, double?, bool)`: ungapped approximate-match count, Primer3 ntthal self-structure screen (fallback screens), optional Kane assessment.
- `ProbeDesigner.AssessCrossHybridization(string, IEnumerable<string>, double, int, bool, ScoringMatrix?)`: Kane et al. (2000) criteria per non-target strand.
- `ProbeDesigner.CheckSpecificity(string, ISuffixTree, bool)`: exact-hit uniqueness from suffix-tree occurrence counts (optionally both strands).
- `ProbeDesigner.ScanOffTargetsGapped(string, IEnumerable<string>, double, ScoringMatrix?)`: Opt-in gapped (Smith–Waterman) off-target scan. Returns a `GappedSpecificityResult` separating `OnTargetHits` (the perfect ungapped full-coverage exact match) from `OffTargetHits` (imperfect/indel hits ≥ `minIdentity`, default 0.75). Reuses `SequenceAligner.LocalAlign` for the indel-aware alignment.
- `ProbeDesigner.ComputeLambdaNucleotide(int, int, double)`: Opt-in. Solves `Σ p_i p_j e^{λ s_ij} = 1` numerically (bisection) for a match/mismatch scheme under uniform base frequencies (`p(match) = 4p²`, `p` in (0, 0.5)); returns the Karlin–Altschul `λ` [8][9]. Throws when the scheme has no positive score or a non-negative expected score, or for `p` outside (0, 0.5).
- `ProbeDesigner.ComputeUngappedKarlinParameters(int, int, double)` / `(int, int, IReadOnlyList<double>)`: ungapped `λ`, `K`, `H`, `α = λ/H`, `β` (BLAST+ `s_GetUngappedBeta`) under uniform or explicit A/C/G/T composition [8][13].
- `ProbeDesigner.GetBlastnGappedKarlinParameters(int reward, int penalty, int gapOpen, int gapExtend)`: BLAST+ gapped `λ, K, H, α, β`, round-down flag (`blastn_values_*` tables, infinite domain → ungapped); unsupported combinations throw `ArgumentException` [13].
- `ProbeDesigner.ComputeLengthAdjustment(double k, double alphaOverLambda, double beta, int m, long n, int N = 1)`: BLAST+ `BLAST_ComputeLengthAdjustment` [13].
- `ProbeDesigner.ComputeBlastnStatistics(int rawScore, int m, long n, int N = 1, ScoringMatrix? = BlastDna, bool gapped = true)` and `ComputeBlastnStatistics(string probe, string subject, ScoringMatrix?)`: `BlastnStatistics` (raw and E-value score, parameters, bit score, E-value, ℓ, effective search space) as blastn computes them; the string overload aligns with `SequenceAligner.LocalAlignAffine` [13].
- `ProbeDesigner.ComputeKarlinAltschul(double, int, long, ScoringMatrix?, double?, double)`: Opt-in. Returns a `KarlinAltschulStatistics` (`RawScore`, `Lambda`, `K`, `BitScore`, `EValue`, `QueryLength`, `DatabaseLength`) for a hit's raw score over a search space `m·n`, using `E = K·m·n·e^{−λS}` and `S' = (λS − ln K)/ln 2` [8][9]. `K` defaults to the value computed for the scoring scheme (`ComputeUngappedKarlinParameters`); ungapped statistics on the raw search space `m·n` (for gapped scores and the edge-effect correction use `ComputeBlastnStatistics`).

### 5.2 Current Behavior

The current validator treats an empty probe as invalid rather than throwing. It records an issue when more than one hit across all references meets a Kane criterion on its ungapped diagonal, when the ntthal self-dimer / 3′ self-dimer Tm exceeds `MaxStructureTm` (fallback: Primer3 alignment-mode self_any / self_end > 12.00), when the ntthal hairpin Tm exceeds it, and for every non-target strand meeting a Kane criterion. On the fallback path the inverted-repeat stem screen (`PrimerDesigner.HasHairpinPotential`) sets `HasSecondaryStructure` and a `Warnings` entry only (library convention, not a validity criterion; audit round 5, A5-5 — a correction of the validity rule, before which it was an issue); raise `ThermodynamicScreenMaxLength` (F55 opt-in) for a sourced hairpin decision on > 60-nt probes. `ValidateProbe`'s approximate matching is an ungapped fixed-length sliding scan with mismatch tolerance (its `OffTargetHits` pools the on-target with off-targets). `IsValid` is `true` exactly when no issue was recorded (the former lenient rule "≤ 1 hit and fold-back fraction ≤ 0.4 ⇒ valid despite issues" had no source and was removed). `CheckSpecificity(...)` uses exact suffix-tree hits (optionally on both strands). With `bothStrands` the reference scan also covers the probe's reverse complement (default: the given strand only, as before).

`ScanOffTargetsGapped(...)` is the opt-in gapped alternative. For each reference it slides a window of length `probeLen + 2` and runs `SequenceAligner.LocalAlign` (Smith–Waterman, `BlastDna` scoring by default), computes identity = identical aligned columns / probe length and coverage = ungapped columns / probe length, keeps every site whose identity ≥ `minIdentity`, and collapses overlapping window detections to one best (highest-identity, then highest-coverage, then leftmost) hit per site via greedy non-overlapping selection. The first perfect ungapped full-coverage exact match (identity = coverage = 1.0, no gaps) is classified as the intended on-target and excluded from the off-target count; all imperfect or indel-containing hits — and any additional perfect repeats — are off-targets.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Ungapped (Hamming) approximate-match cross-hybridization screening (`ValidateProbe`), each site judged by the Kane et al. (2000) identity / contiguous-identity criteria [7]; opt-in on both strands of the references (`bothStrands`, as blastn `-strand both`; palindromic sites once).
- Primer3 thermodynamic hybridization-probe self-structure screen (ntthal self-dimer / 3′ self-dimer / hairpin Tm ≤ 47 °C at the internal-oligo conditions; primer3-py parity) [11].
- Kane et al. (2000) cross-hybridization criteria: identity > 75 % or contiguous identical stretch > 15 nt, both strands (`AssessCrossHybridization`, optional in `ValidateProbe`) [7].
- Thermodynamic stability of each off-target site: ntthal duplex Tm of the probe with the strand complementary to the aligned site at the stated conditions (`CrossHybridizationAssessment.DuplexTm`, primer3-py `calc_heterodimer` parity — computed when the probe or the site is ≤ the conditions' `ThermodynamicScreenMaxLength`, THAL_MAX_ALIGN default 60, opt-in up to 10 000 = thal.c compiled with `-DTHAL_MAX_ALIGN`, A3-27; the duplex-Tm cross-hybridization check of OligoArray 2.0 [12]); reported, and flagged only against an optional caller threshold `maxDuplexTm` (OligoArray: a cross-hybridization with Tm above the user's specificity threshold makes the probe non-specific; the stringency is assay-specific, so there is no default).
- Gapped (Smith–Waterman) local-alignment off-target scan with indel handling [1][2] and on/off-target separation (`ScanOffTargetsGapped`).
- Sequence-only self-complementarity / stem screens as the documented fallback (> 60 nt: thal.c THAL_MAX_ALIGN; non-ACGT; `Heuristic`); the hairpin part is the canonical `PrimerDesigner.HasHairpinPotential` (stem ≥ 4 bp library convention, loop ≥ 3 nt thal.c; A5-2), reported as a warning that does not decide `IsValid` (A5-5).
- Uniqueness scoring as `0` or `1 / hits` (`ValidateProbe`, `CheckSpecificity`) and the default search radius of 3 mismatches — library conventions, not published metrics; the score does not enter `IsValid`.
- Off-target identity threshold (default 0.75 over the probe length) per Kane et al. (2000) [7].
- Karlin–Altschul E-value, bit score, and the `λ` defining equation (`ComputeKarlinAltschul` / `ComputeLambdaNucleotide`), with the negative-expected-score and at-least-one-positive-score preconditions [8][9]; `K` and `H` computed as BLAST+ computes them; BLAST+ gapped blastn parameters, edge-effect length adjustment and effective search space, even-score round-down (`ComputeBlastnStatistics`, equal to NCBI blastn 2.12.0+ E-values) [13].

**Intentionally simplified:**

- `ValidateProbe`'s uniqueness score collapses all multi-hit outcomes to `1 / hits`; **consequence:** it distinguishes hit multiplicity only; mismatch severity is assessed by the Kane identity/contiguity criteria (`AssessCrossHybridization`).
- `ValidateProbe`'s approximate matching is substitution-only and fixed-length, and its `OffTargetHits` pools the on-target match with off-targets; **consequence:** for indel-aware detection and on/off separation use `ScanOffTargetsGapped` instead.
- Suffix-tree specificity uses exact hits only; **consequence:** approximate off-targets are only modeled through `ValidateProbe(...)`/`ScanOffTargetsGapped(...)`, not through `CheckSpecificity(...)`.
- On/off-target labelling: the first perfect ungapped full-coverage exact match is taken as the intended on-target; **consequence:** when several identical perfect sites exist, the first is on-target and the rest are off-targets.

**Not implemented:**

- A seeded BLAST k-mer index over a whole genome; `ScanOffTargetsGapped` is an exhaustive sliding Smith–Waterman scan (O(g · n·m)), not a genome-scale seed-and-extend index [2]; **users should rely on:** an external seeded aligner for genome-scale off-target *performance* (the exhaustive scan already finds every hit a seed would, so this is a speed, not a correctness, gap).


## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty probe sequence | Returns a structured invalid result | Explicit special case in source |
| Unique probe hit | `SpecificityScore = 1.0` | Exact mapping in implementation |
| No probe hits | `SpecificityScore = 0.0` | Probe does not match the references |
| Multiple hits | `SpecificityScore = 1.0 / hits` | Cross-hybridization penalty |
| `bothStrands`, reverse-palindromic probe or site | One site per reference position (not two) | Same rule as `CheckSpecificity(bothStrands)` |

### 6.2 Limitations

The implementation is a screening tool. The opt-in `ScanOffTargetsGapped` adds indel-aware (gapped) off-target detection and on/off-target separation, and `ComputeKarlinAltschul` adds the Karlin–Altschul bit-score / E-value significance of a hit; but off-target search remains an exhaustive sliding Smith–Waterman scan, not a seeded BLAST k-mer index over a whole genome (a genome-scale *performance* technique — the exhaustive scan already finds every hit a seed would). Off-target sites assessed by `AssessCrossHybridization` carry the ntthal duplex Tm (ACGT; computed when the probe or the site is ≤ the conditions' `ThermodynamicScreenMaxLength`, default 60 = Primer3 THAL_MAX_ALIGN, opt-in up to 10 000 — F55/F56); mismatch-position weighting and assay stringency are not modelled, and the suffix-tree helper only captures exact-hit uniqueness.

## 8. References

1. Smith TF, Waterman MS (1981) - Identification of common molecular subsequences, J. Mol. Biol. 147(1):195–197. https://doi.org/10.1016/0022-2836(81)90087-5 (recurrence via https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
2. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990) - Basic local alignment search tool, J. Mol. Biol. 215(3):403–410. https://doi.org/10.1016/S0022-2836(05)80360-2 (via https://en.wikipedia.org/wiki/BLAST_(biotechnology))
3. Wikipedia: Hybridization probe - https://en.wikipedia.org/wiki/Hybridization_probe
4. Wikipedia: DNA microarray - https://en.wikipedia.org/wiki/DNA_microarray
5. Wikipedia: Off-target genome editing - https://en.wikipedia.org/wiki/Off-target_genome_editing
6. Amann R, Ludwig W (2000) - Ribosomal RNA-targeted nucleic acid probes for studies in microbial ecology, FEMS Microbiology Reviews.
7. Kane MD, Jatkoe TA, Stumpf CR, Lu J, Thomas JD, Madore SJ (2000) - Assessment of the sensitivity and specificity of oligonucleotide (50mer) microarrays, Nucleic Acids Research 28(22):4552–4557. https://pmc.ncbi.nlm.nih.gov/articles/PMC113865/
8. Karlin S, Altschul SF (1990) - Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes, PNAS 87(6):2264–2268. https://doi.org/10.1073/pnas.87.6.2264 (formulas via https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html and http://www.cs.cmu.edu/~durand/03-711/2011/Lectures/Blast-informationContent-2011.pdf)
9. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ (1990) - Basic local alignment search tool (E-value / bit-score statistics), J. Mol. Biol. 215(3):403–410. https://doi.org/10.1016/S0022-2836(05)80360-2 (statements via the NCBI tutorial "The Statistics of Sequence Similarity Scores")
10. Satya RV, Zavaljevski N, Kumar K, Reifman J (2008) - A high-throughput pipeline for designing microarray-based pathogen diagnostic assays, BMC Bioinformatics 9:185 (summary of the Kane criteria: "overall sequence identity is > 75% or ... a contiguous match > 15 bp"); Chen H, Sharp BM (2002) Oliz, BMC Bioinformatics 3:27 (Kane criteria restated). Both via WebSearch snippets (publisher pages blocked).
11. Untergasser A et al. (2012) Primer3 — new capabilities and interfaces, Nucleic Acids Res 40:e115; primer3 `libprimer3.cc` (`o_args` internal-oligo defaults, `oligo_compl_thermod`, `oligo_hairpin`); primer3-py 2.3.1 `calc_homodimer` / `calc_end_stability` / `calc_hairpin`.
12. Rouillard JM, Zuker M, Gulari E (2003) - OligoArray 2.0: design of oligonucleotide probes for DNA microarrays using a thermodynamic approach, Nucleic Acids Res 31(12):3057–3062 (cross-hybridization by the duplex stability of the probe with each BLAST hit).
13. NCBI C++ Toolkit BLAST+ core: `src/algo/blast/core/blast_stat.c` (`BlastScoreFreqCalc`, `Blast_KarlinBlkUngappedCalc`, `BlastKarlinLtoH`, `BlastKarlinLHtoK`, `blastn_values_*`, `Blast_KarlinBlkNuclGappedCalc`, `Blast_GetNuclAlphaBeta`, `BLAST_ComputeLengthAdjustment`), `blast_setup.c` (`BLAST_CalcEffLengths`), `blast_hits.c` (`Blast_HSPListGetEvalues`, `Blast_HSPListGetBitScores`). https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/blast/core/blast_stat.c (retrieved 2026-10-01); NCBI blastn 2.12.0+ used as the numerical oracle.
