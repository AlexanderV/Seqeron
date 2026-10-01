# Probe Validation

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PROBE-VALID-001 |
| Related Projects | N/A |
| Implementation Status | Complete (Primer3 thermodynamic self-structure screen; Kane et al. 2000 cross-hybridization criteria) |
| Last Reviewed | 2026-10-01 |

## 1. Overview

Probe validation assesses whether a hybridization probe is likely to bind specifically to its intended targets with limited cross-hybridization. In this repository, `ValidateProbe` combines (1) substitution-tolerant fixed-length window matching against reference sequences, (2) Primer3's thermodynamic hybridization-probe self-structure screen — ntthal self-dimer, 3′ self-dimer and hairpin Tm ≤ 47 °C, the same screen `DesignProbes` uses — and (3) optionally the Kane et al. (2000) cross-hybridization criteria against known non-target sequences (`AssessCrossHybridization`: overall identity > 75 % or a contiguous identical stretch > 15 nt) [7]. The result is a validation record with a uniqueness score, the ntthal Tm values, per-non-target assessments, an issue list and `IsValid` = no issue.

An opt-in **gapped** off-target scan (`ScanOffTargetsGapped`) supplements the default ungapped scan: it reuses the library's validated Smith–Waterman local aligner [1] to find off-target sites reachable through insertions or deletions — the indel-aware "BLAST-grade" improvement over a pure ungapped Hamming scan [2] — and it separates the single intended on-target match from genuine off-target hits (correcting the on/off-target pooling of `ValidateProbe`'s `OffTargetHits`). An opt-in **Karlin–Altschul statistics** layer (`ComputeKarlinAltschul` / `ComputeLambdaNucleotide`) quantifies the statistical significance of an off-target hit's raw alignment score as a bit score and an E-value [8][9]. `CheckSpecificity` counts exact suffix-tree occurrences, optionally on both strands.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Cross-hybridization occurs when a probe binds sequences other than its intended target, and this is a central design concern for FISH, DNA microarrays, qPCR, and related assays. The original document also notes that mismatch tolerance, assay stringency, and probe length affect specificity, and that self-complementarity and low-complexity content can increase the risk of non-specific behavior. Sources: Wikipedia (Hybridization probe, DNA microarray, BLAST), Altschul et al. (1990), Amann & Ludwig (2000).

Off-target detection by *local alignment* follows the Smith–Waterman recurrence [1]: `H(i,j) = max{ H(i-1,j-1) + s(a_i,b_j), H(i-1,j) - W, H(i,j-1) - W, 0 }`, whose zero floor returns the best-scoring local subsequence match and whose gap terms admit insertions/deletions. The gapped local alignment is the recognized improvement over ungapped matching for similarity/homology search [2]. **Kane et al. (2000) criteria [7].** Kane et al. hybridized 50-mer oligonucleotide arrays with non-target cDNAs of graded similarity: "any 'non-target' transcripts (cDNAs) >75% similar over the 50 base target may show cross-hybridization", and a non-target must not share "a stretch of complementary sequence >15 contiguous bases" (abstract; later probe-design pipelines summarise the rule as "overall sequence identity > 75 % or a contiguous match > 15 bp" [10]). `AssessCrossHybridization` implements both criteria per non-target strand: identity = identical columns of the best local alignment (Smith–Waterman–Gotoh with BLAST+ blastn scoring +2/−3, gap existence 5 / extension 2 — the canonical `SequenceAligner.LocalAlignAffine`) ÷ probe length, flagged when strictly > `maxIdentity` (0.75); longest contiguous match = longest common substring (canonical suffix tree), flagged when strictly > `maxContiguousMatch` (15). Both strands are assessed by default (a double-stranded non-target offers both). `ScanOffTargetsGapped` uses the 0.75 value as its site-reporting threshold (identity ≥ `minIdentity`).

**Self-structure (Primer3) [11].** Primer3 screens a hybridization probe (internal oligo) with ntthal: self-dimer (ANY) Tm ≤ PRIMER_INTERNAL_MAX_SELF_ANY_TH, 3′ self-dimer (END1) Tm ≤ PRIMER_INTERNAL_MAX_SELF_END_TH, hairpin Tm ≤ PRIMER_INTERNAL_MAX_HAIRPIN_TH, all 47 °C by default, at the internal-oligo conditions (50 nM, 50 mM monovalent, 0 Mg²⁺, 0 dNTP). `ValidateProbe` applies this screen (via `PrimerDesigner.CalculatePrimer3OligoStructure`, bit-exact to primer3-py) to ≤ 60-nt A/C/G/T probes; longer probes (thal.c `THAL_MAX_ALIGN` = 60) and non-ACGT probes fall back to the sequence-only fold-back fraction and inverted-repeat screens.

### 2.2 Core Model

The validation workflow normalizes the probe to uppercase, counts approximate matches across all supplied reference sequences, and maps the hit count to a **library-defined uniqueness score** (the share of the probe's N candidate binding sites taken by one site; not a published specificity metric — the sourced cross-hybridization decision is the Kane assessment):

$$
specificity =
\begin{cases}
0.0, & offTargetHits = 0 \\
1.0, & offTargetHits = 1 \\
1.0 / offTargetHits, & offTargetHits > 1
\end{cases}
$$

It then applies the Primer3 thermodynamic self-structure screen (fallback: fold-back fraction against the reverse complement and a sequence-level inverted-repeat screen), the optional Kane assessment, and sets `IsValid` = no issue recorded.

**Karlin–Altschul statistics of an off-target hit (opt-in).** For a hit's raw alignment score `S` against a search space of query length `m` and database length `n`, the statistical significance follows the Karlin–Altschul framework [8][9]:

$$E = K \cdot m \cdot n \cdot e^{-\lambda S}, \qquad S' = \frac{\lambda S - \ln K}{\ln 2}, \qquad E = m \cdot n \cdot 2^{-S'}$$

where the scale parameter `λ` is the unique positive root of the defining equation [8]

$$\sum_{i,j} p_i\, p_j\, e^{\lambda s_{ij}} = 1$$

with `p_i` the background base frequencies and `s_ij` the score matrix. For four equiprobable bases (`p_i = 0.25`) and a match/mismatch scheme this reduces to `0.25·e^{λ·match} + 0.75·e^{λ·mismatch} = 1`; for the BLAST `+1/−3` scheme it solves to `λ ≈ 1.374` (matching the value NCBI blastn reports). The theory requires a scoring scheme with **negative expected per-pair score** and **at least one positive score** [9]; both are checked. `K` (whose full closed form needs the Karlin–Altschul score-lattice machinery) is exposed as a caller parameter, defaulted to the published nucleotide value `0.711`.

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

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `probeSequence` | `string` | required | Probe sequence to validate | Null input throws `ArgumentNullException`; empty string yields a structured invalid result |
| `referenceSequences` | `IEnumerable<string>` | required | Reference sequences scanned for approximate matches | Null input throws `ArgumentNullException` |
| `maxMismatches` | `int` | `3` | Maximum mismatch tolerance for approximate matching | Passed through to the internal approximate-match search |
| `selfComplementarityThreshold` | `double` | `0.3` | Fold-back-fraction limit of the sequence-only fallback screen | Used only when the thermodynamic screen does not apply |
| `conditions` | `ProbeParameters?` | `Defaults.Microarray` | Salt / oligo concentrations, `StructureScreen`, `MaxStructureTm` (47 °C) | Primer3 internal-oligo defaults |
| `nonTargetSequences` | `IEnumerable<string>?` | `null` | Known non-targets for the Kane assessment | Optional |
| `maxNonTargetIdentity` / `maxContiguousMatch` | `double` / `int` | `0.75` / `15` | Kane thresholds (strict `>`) | [7] |
| `maxDuplexTm` | `double?` | `null` | Optional off-target site duplex-Tm threshold (°C) | [12] |
| `genomeIndex` | `ISuffixTree` | required for `CheckSpecificity(...)` | Pre-built suffix tree for exact hit counting | Used only by the suffix-tree specificity helper |
| `bothStrands` (`CheckSpecificity`) | `bool` | `false` | Also count reverse-complement occurrences (other strand of a ds genome; palindromes once) | Optional |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `IsValid` | `bool` | True when no issue was recorded |
| `SpecificityScore` | `double` | Specificity value in the range `0.0-1.0` |
| `OffTargetHits` | `int` | Total approximate hits across all reference sequences |
| `SelfComplementarity` | `double` | Fold-back fraction (criterion only of the fallback screen) |
| `HasSecondaryStructure` | `bool` | ntthal hairpin Tm > `MaxStructureTm` (fallback: inverted-repeat stem) |
| `Issues` | `IReadOnlyList<string>` | Recorded validation issues |
| `ThermodynamicScreen` | `bool` | True when the Primer3 ntthal screen applied |
| `SelfDimerTm` / `SelfEndDimerTm` / `HairpinTm` | `double?` | ntthal Tm (°C; Primer3 reports negative Tm as 0); null for the fallback |
| `CrossHybridization` | `IReadOnlyList<CrossHybridizationAssessment>` | Per non-target strand: identity, identical columns, score, longest contiguous match, the two Kane flags |

### 3.3 Preconditions and Validation

`ValidateProbe(...)` uppercases the input probe before analysis. Null probe or reference collections raise `ArgumentNullException`. An empty probe sequence returns a structured invalid result with `SpecificityScore = 0.0`, `OffTargetHits = 0`, and an `"Empty probe sequence"` issue. `CheckSpecificity(...)` uppercases the probe before querying the suffix tree.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the probe sequence to uppercase.
2. Search every reference sequence for approximate matches within the mismatch tolerance.
3. Accumulate the total number of hits across all references.
4. Self-structure: ntthal self-dimer / 3′ self-dimer / hairpin Tm vs `MaxStructureTm` (≤ 60-nt ACGT, thermodynamic screen) or the sequence-only screens.
5. Optional Kane assessment of each non-target (both strands).
6. Map hit count to the uniqueness score; `IsValid` = no issue.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Validation defaults preserved from the original document and source:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `maxMismatches` | `3` | Approximate-match tolerance |
| `selfComplementarityThreshold` | `0.3` | Fallback-screen fold-back-fraction limit |
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

- `ProbeDesigner.ValidateProbe(string, IEnumerable<string>, int, double, ProbeParameters?, IEnumerable<string>?, double, int)`: ungapped approximate-match count, Primer3 ntthal self-structure screen (fallback screens), optional Kane assessment.
- `ProbeDesigner.AssessCrossHybridization(string, IEnumerable<string>, double, int, bool, ScoringMatrix?)`: Kane et al. (2000) criteria per non-target strand.
- `ProbeDesigner.CheckSpecificity(string, ISuffixTree, bool)`: exact-hit uniqueness from suffix-tree occurrence counts (optionally both strands).
- `ProbeDesigner.ScanOffTargetsGapped(string, IEnumerable<string>, double, ScoringMatrix?)`: Opt-in gapped (Smith–Waterman) off-target scan. Returns a `GappedSpecificityResult` separating `OnTargetHits` (the perfect ungapped full-coverage exact match) from `OffTargetHits` (imperfect/indel hits ≥ `minIdentity`, default 0.75). Reuses `SequenceAligner.LocalAlign` for the indel-aware alignment.
- `ProbeDesigner.ComputeLambdaNucleotide(int, int, double)`: Opt-in. Solves `Σ p_i p_j e^{λ s_ij} = 1` numerically (bisection) for a match/mismatch scheme under uniform base frequencies; returns the Karlin–Altschul `λ` [8][9]. Throws when the scheme has no positive score or a non-negative expected score.
- `ProbeDesigner.ComputeKarlinAltschul(double, int, long, ScoringMatrix?, double, double)`: Opt-in. Returns a `KarlinAltschulStatistics` (`RawScore`, `Lambda`, `K`, `BitScore`, `EValue`, `QueryLength`, `DatabaseLength`) for a hit's raw score over a search space `m·n`, using `E = K·m·n·e^{−λS}` and `S' = (λS − ln K)/ln 2` [8][9]. `K` defaults to the published nucleotide value `0.711`.

### 5.2 Current Behavior

The current validator treats an empty probe as invalid rather than throwing. It records an issue when more than one hit is found across all references, when the ntthal self-dimer / 3′ self-dimer Tm exceeds `MaxStructureTm` (fallback: fold-back fraction > threshold), when the ntthal hairpin Tm exceeds it (fallback: inverted-repeat stem), and for every non-target strand meeting a Kane criterion. `ValidateProbe`'s approximate matching is an ungapped fixed-length sliding scan with mismatch tolerance (its `OffTargetHits` pools the on-target with off-targets). `IsValid` is `true` exactly when no issue was recorded (the former lenient rule "≤ 1 hit and fold-back fraction ≤ 0.4 ⇒ valid despite issues" had no source and was removed). `CheckSpecificity(...)` uses exact suffix-tree hits (optionally on both strands).

`ScanOffTargetsGapped(...)` is the opt-in gapped alternative. For each reference it slides a window of length `probeLen + 2` and runs `SequenceAligner.LocalAlign` (Smith–Waterman, `BlastDna` scoring by default), computes identity = identical aligned columns / probe length and coverage = ungapped columns / probe length, keeps every site whose identity ≥ `minIdentity`, and collapses overlapping window detections to one best (highest-identity, then highest-coverage, then leftmost) hit per site via greedy non-overlapping selection. The first perfect ungapped full-coverage exact match (identity = coverage = 1.0, no gaps) is classified as the intended on-target and excluded from the off-target count; all imperfect or indel-containing hits — and any additional perfect repeats — are off-targets.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Ungapped (Hamming) approximate-match cross-hybridization screening (`ValidateProbe`).
- Primer3 thermodynamic hybridization-probe self-structure screen (ntthal self-dimer / 3′ self-dimer / hairpin Tm ≤ 47 °C at the internal-oligo conditions; primer3-py parity) [11].
- Kane et al. (2000) cross-hybridization criteria: identity > 75 % or contiguous identical stretch > 15 nt, both strands (`AssessCrossHybridization`, optional in `ValidateProbe`) [7].
- Thermodynamic stability of each off-target site: ntthal duplex Tm of the probe with the strand complementary to the aligned site at the stated conditions (`CrossHybridizationAssessment.DuplexTm`, primer3-py `calc_heterodimer` parity; the duplex-Tm cross-hybridization check of OligoArray 2.0 [12]); reported, and flagged only against an optional caller threshold `maxDuplexTm` (OligoArray: a cross-hybridization with Tm above the user's specificity threshold makes the probe non-specific; the stringency is assay-specific, so there is no default).
- Gapped (Smith–Waterman) local-alignment off-target scan with indel handling [1][2] and on/off-target separation (`ScanOffTargetsGapped`).
- Sequence-only self-complementarity / stem screens as the documented fallback (> 60 nt: thal.c THAL_MAX_ALIGN; non-ACGT; `Heuristic`).
- Uniqueness scoring as `0` or `1 / hits` (`ValidateProbe`, `CheckSpecificity`) — library-defined, not a published metric.
- Off-target identity threshold (default 0.75 over the probe length) per Kane et al. (2000) [7].
- Karlin–Altschul E-value, bit score, and the `λ` defining equation (`ComputeKarlinAltschul` / `ComputeLambdaNucleotide`), with the negative-expected-score and at-least-one-positive-score preconditions [8][9].

**Intentionally simplified:**

- `ValidateProbe`'s uniqueness score collapses all multi-hit outcomes to `1 / hits`; **consequence:** it distinguishes hit multiplicity only; mismatch severity is assessed by the Kane identity/contiguity criteria (`AssessCrossHybridization`).
- `ValidateProbe`'s approximate matching is substitution-only and fixed-length, and its `OffTargetHits` pools the on-target match with off-targets; **consequence:** for indel-aware detection and on/off separation use `ScanOffTargetsGapped` instead.
- Suffix-tree specificity uses exact hits only; **consequence:** approximate off-targets are only modeled through `ValidateProbe(...)`/`ScanOffTargetsGapped(...)`, not through `CheckSpecificity(...)`.
- On/off-target labelling: the first perfect ungapped full-coverage exact match is taken as the intended on-target; **consequence:** when several identical perfect sites exist, the first is on-target and the rest are off-targets.
- The Karlin–Altschul `K` parameter is a caller-supplied value (default the published nucleotide `0.711`) rather than computed from its full closed form; **consequence:** `λ`, the bit score, and the E-value's score-dependence are exact, but `K`'s value depends on the supplied constant (use the matching published `K` for a non-default scoring scheme).

**Not implemented:**

- A seeded BLAST k-mer index over a whole genome; `ScanOffTargetsGapped` is an exhaustive sliding Smith–Waterman scan (O(g · n·m)), not a genome-scale seed-and-extend index [2]; **users should rely on:** an external seeded aligner for genome-scale off-target *performance* (the exhaustive scan already finds every hit a seed would, so this is a speed, not a correctness, gap).
- The Karlin–Altschul `K` closed form (score-probability lattice / geometric-spacing machinery of [8]); **users should rely on:** the caller-supplied `K` parameter (published values per scoring scheme).


## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty probe sequence | Returns a structured invalid result | Explicit special case in source |
| Unique probe hit | `SpecificityScore = 1.0` | Exact mapping in implementation |
| No probe hits | `SpecificityScore = 0.0` | Probe does not match the references |
| Multiple hits | `SpecificityScore = 1.0 / hits` | Cross-hybridization penalty |

### 6.2 Limitations

The implementation is a screening tool. The opt-in `ScanOffTargetsGapped` adds indel-aware (gapped) off-target detection and on/off-target separation, and `ComputeKarlinAltschul` adds the Karlin–Altschul bit-score / E-value significance of a hit; but off-target search remains an exhaustive sliding Smith–Waterman scan, not a seeded BLAST k-mer index over a whole genome (a genome-scale *performance* technique — the exhaustive scan already finds every hit a seed would). Off-target sites assessed by `AssessCrossHybridization` carry the ntthal duplex Tm (≤ 60-nt ACGT probes); mismatch-position weighting and assay stringency are not modelled, and the suffix-tree helper only captures exact-hit uniqueness.

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
