# Hybridization Probe Design

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | PROBE-DESIGN-001 |
| Related Projects | N/A |
| Implementation Status | Complete (Primer3 probe picker + heuristic ranking designer) |
| Last Reviewed | 2026-10-01 |

## 1. Overview

Hybridization probe design generates oligonucleotide probes that detect specific nucleic acid sequences through complementary binding. In this repository, probe design supports multiple applications including FISH, DNA microarrays, Northern blots, qPCR, and Southern blots. The documented implementation is a heuristic candidate-generation and ranking workflow driven by GC content, melting temperature, self-complementarity, secondary structure, simple repeats, and application-specific parameter ranges.

A **Primer3 hybridization-probe picker** (`DesignProbesPrimer3`) reproduces Primer3's `PRIMER_TASK=pick_hyb_probe_only` (internal-oligo) selection exactly (§4.4).

An **opt-in TaqMan (5'-nuclease hydrolysis probe) rule set** is also provided (`EvaluateTaqManProbe`, `SelectTaqManStrand`) for the qPCR/TaqMan chemistry. It is separate from and does not alter the generic designer, which remains the default.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Hybridization probes are used in sequence-detection assays such as FISH, DNA microarrays, Northern blots, and qPCR. Probe quality depends on a balance among target complementarity, melting temperature, GC content, and avoidance of self-structure or repetitive motifs that can reduce specificity. The original document ties these considerations to nucleic-acid thermodynamics and to application-specific probe lengths and Tm targets. Sources: Wikipedia (Hybridization probe, Fluorescence in situ hybridization, DNA microarray, Nucleic acid thermodynamics), SantaLucia (1998), Breslauer et al. (1986).

### 2.2 Core Model

The candidate score starts at `1.0` and is reduced by a fixed set of penalties:

| Factor | Penalty |
|--------|---------|
| GC content outside the configured range | `-0.3` |
| Tm outside the configured range | `-0.3` |
| Homopolymer run above the configured maximum | `-0.2` |
| Self-complementarity above the configured maximum | `-0.2` |
| Secondary-structure potential | `-0.15` |
| Simple repeats | `-0.1` |
| Starts with `G/C` | `-0.02` |
| Ends with `G/C` | `-0.02` |

**Probe Tm.** Every probe Tm (designer, tiling, molecular beacon, TaqMan check, `AnalyzeOligo`) is Primer3's `seqtm` (`PrimerDesigner.CalculateMeltingTemperaturePrimer3`, identical to primer3-py `calc_tm(..., max_nn_length=36)`, the MAX_NN_TM_LENGTH libprimer3 uses when picking oligos): SantaLucia (1998) nearest-neighbour ΔH/ΔS with the SantaLucia entropy salt correction and the von Ahsen Mg²⁺→Na⁺-equivalent for ≤ 36 nt, and Primer3 `long_seq_tm` = 81.5 + 16.6·log10([Na⁺]eq) + 0.41·%GC − 600/N above 36 nt, at the hybridization conditions in `ProbeParameters` (default = Primer3 internal-oligo defaults: 50 nM probe, 50 mM monovalent, 0 Mg²⁺, 0 dNTP). A probe with a non-ACGT base has no computable Tm (reported 0 with a warning and the Tm penalty; `AnalyzeOligo`/TaqMan report NaN).

**Self-structure.** With `StructureScreen = Thermodynamic` (default) a probe of ≤ 60 nt of A/C/G/T is screened like a Primer3 hybridization probe: ntthal self-dimer (ANY), 3′ self-dimer (END1) and hairpin Tm (`PrimerDesigner.CalculatePrimer3OligoStructure`, bit-exact to primer3-py) must not exceed `MaxStructureTm` (PRIMER_INTERNAL_MAX_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH = 47 °C); the self-dimer limit drives the "self-complementarity" penalty and the hairpin limit the "secondary structure" penalty. Longer probes (thal.c THAL_MAX_ALIGN = 60), non-ACGT probes and `StructureScreen = Heuristic` use the fallback screens: the self-dimer criterion is Primer3's alignment-mode internal-oligo screen (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0, `oligo_compl`): dpal `self_any` (`PrimerDesigner.CalculatePrimerSelfAnyComplementarity`, DPAL_LOCAL) and `self_end` (`CalculatePrimerSelfEndComplementarity`, DPAL_GLOBAL_END), which have no length limit, must not exceed `MaxSelfAny` / `MaxSelfEnd` (PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END = 12.00); the hairpin criterion stays the sequence-only inverted-repeat stem ≥ 4 bp with a 3-nt loop and ≥ 80 % matches (a thermodynamic hairpin for > 60-nt DNA needs a DNA-parameter MFE fold — cross-batch request to B12). Note: Primer3 itself picks internal oligos of ≤ 36 nt, so 12.00 is calibrated for short oligos; random 200–500-nt probes exceed it by chance in roughly half of the cases (a −0.2 ranking penalty in `DesignProbes`) — raise `MaxSelfAny`/`MaxSelfEnd` for long FISH/blot probes. `MaxSelfComplementarity` (fold-back fraction limit) is no longer used by any screen (kept for source compatibility).

**Simple repeats** are di-/trinucleotide microsatellites of ≥ 4 complete copies (canonical `RepeatFinder.FindMicrosatellites`, MISA convention: primitive units, A/C/G/T only).

#### 2.2.1 TaqMan (5'-nuclease hydrolysis probe) rules — opt-in

For TaqMan qPCR chemistry, the published Applied Biosystems / Thermo Fisher guidelines [7][8][9] add chemistry-specific constraints that the generic designer does not enforce. `EvaluateTaqManProbe` checks each as a boolean; `SelectTaqManStrand` chooses the better strand.

| Rule | Threshold | Source |
|------|-----------|--------|
| No G at the 5' end | first base ≠ `G` (a 5' G adjacent to the reporter dye quenches reporter fluorescence even after cleavage) | [7][9] |
| More Cs than Gs | `count(C) > count(G)` | [7] |
| No run of ≥4 consecutive Gs | max G-run `< 4` | [7] |
| G+C content | `30%–80%` | [7] |
| Probe length | `18–22 nt` (default; configurable) | [7] |
| Probe Tm above primer Tm | probe Tm `≥ primerTm + 10 °C` | [7][8] |

If a 5' G cannot be avoided on the sense strand, the probe is designed on the complement (antisense) strand [8]; `SelectTaqManStrand` returns the reverse-complement strand when it better satisfies the rules.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The base ranking pass retains only raw-score-positive candidates before optional specificity rescaling | `DesignProbesOptimized(...)` rejects candidates when `score <= 0`, but the suffix-tree overload can later rescale shortlisted scores |
| INV-02 | Probe GC content is computed as a fraction of length | The source uses `gcCount / length` |
| INV-03 | `CheckSpecificity(...)` returns `0` for no hits, `1` for a unique hit, and `1 / hits` otherwise | That mapping is explicit in source |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `targetSequence` | `string` | required | Sequence from which probe candidates are generated | Uppercased before processing |
| `parameters` | `ProbeParameters?` | `Defaults.Microarray` | Application-specific probe design limits | Includes length, Tm, GC, homopolymer, and self-complementarity thresholds |
| `maxProbes` | `int` | `10` | Maximum number of returned probes | Applied after ranking |
| `genomeIndex` | `ISuffixTree` | required for specificity overload | Pre-built suffix tree for genome-wide uniqueness filtering | Used only by the overload with specificity checking |
| `requireUnique` | `bool` | `true` | Whether non-unique probes are excluded when `genomeIndex` is provided | Filters candidates with specificity `< 1.0` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Sequence` | `string` | Probe sequence |
| `Start` / `End` | `int` | Probe coordinates in the source sequence |
| `Tm` | `double` | Melting-temperature estimate |
| `GcContent` | `double` | GC fraction |
| `Score` | `double` | Heuristic quality score |
| `Type` | `ProbeType` | Probe category such as `Standard`, `Tiling`, `Antisense`, `LNA`, or `MolecularBeacon` |
| `Warnings` | `IReadOnlyList<string>` | Quality warnings recorded during evaluation |

### 3.3 Preconditions and Validation

`DesignProbes(...)` returns no probes when the target sequence is null, empty, or shorter than the configured minimum length. All sequences are converted to uppercase before processing. The suffix-tree overload first builds a larger raw-score shortlist than the final requested count, then either filters that shortlist for uniqueness or scales shortlisted scores by the specificity value returned by `CheckSpecificity(...)`; when `requireUnique` is `false`, a shortlisted probe can therefore remain in the output with final score `0` if specificity is `0`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input sequence to uppercase.
2. Precompute GC prefix sums for the target sequence.
3. Enumerate all candidate windows within the configured length range.
4. Evaluate each candidate for GC content, Tm, homopolymers, self-complementarity, secondary structure, repeats, and terminal G/C penalties.
5. Keep only candidates with positive raw scores, sort them by score, and return the top results.
6. In the genome-index overload, apply suffix-tree specificity filtering or post-shortlist score adjustment before yielding the final probes.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Application-specific defaults from `ProbeDesigner.Defaults`:

| Application | Length (bp) | Tm (°C) | GC Content | Notes |
|-------------|-------------|---------|------------|-------|
| Microarray | 50-60 | 75-85 | 0.40-0.60 | Default parameter set |
| FISH | 200-500 | 70-90 | 0.35-0.65 | Higher self-complementarity tolerance |
| Northern Blot | 100-300 | 65-80 | 0.40-0.60 | Intermediate probe sizes |
| qPCR | 20-30 | 68-70 | 0.30-0.80 | TaqMan probe Tm/GC per Applied Biosystems Primer Express guidelines [10] |
| Southern Blot | 150-500 | 65-75 | 0.35-0.65 | Long-probe setting |

The implementation uses prefix sums for O(1) GC lookup and a suffix-tree overload for O(m) uniqueness checking of each candidate probe. The thermodynamic self-structure screens only lower a score, so they are evaluated lazily in descending base-score order with a branch-and-bound stop (identical result to an exhaustive evaluation; ties by length, then start).

The Microarray, FISH, Northern and Southern Tm windows are library conventions (assay- and tool-specific; e.g. OligoArray designs use other conditions) — set `MinTm`/`MaxTm` and the hybridization conditions for your protocol. Lengths: microarray long oligos 50–60 nt (Kane et al. 2000 50-mers; Agilent 60-mers).

### 4.4 Primer3 hybridization-probe picker (`DesignProbesPrimer3`)

Exact port of Primer3 `make_internal_oligo_list` → `pick_primer_range` → `calc_and_check_oligo_features` → `p_obj_fn` (internal-oligo branch) in the default thermodynamic mode, with `Primer3ProbeSettings` = the `PRIMER_INTERNAL_*` defaults (size 18/20/27, Tm 57/60/63 °C, GC 20–80 %, max poly-X 5, Ns 0, self-any/self-end/hairpin Tm ≤ 47 °C, 50 mM monovalent, 0 Mg²⁺, 0 dNTP, 50 nM). Accepted windows are ranked by penalty = |Tm − OptTm| + |len − OptSize| (+ PRIMER_INTERNAL_WT_GC_PERCENT_GT/_LT × the GC deviation from PRIMER_INTERNAL_OPT_GC_PERCENT when `WeightGcPercentGt/Lt` ≠ 0, optimum `OptGcPercent`, null = 50 %; illegal conditions — salt or DNA concentration ≤ 0, negative divalent/dNTP — throw as Primer3 `_pr_data_control`; audit round 3, A3-1, primer3-py parity in `Primer_Design.md` §5.3) and ordered as `primer_rec_comp` (penalty ↑, start ↓, length ↑); `numReturn` = PRIMER_NUM_RETURN (5; < 1 throws `ArgumentOutOfRangeException`, Primer3 `_pr_data_control` "PRIMER_NUM_RETURN < 1" — audit round 3, A3-15). Verified against primer3-py 2.3.1 `design_primers` on 950 random templates/settings (positions, Tm, penalty, SELF_ANY_TH, SELF_END_TH, HAIRPIN_TH identical, max |Δ| = 0).

Mishybridization library (PRIMER_INTERNAL_MISHYB_LIBRARY; audit round 3, A3-3 part 2): `Primer3ProbeSettings.MishybLibrary` (a `PrimerMisprimingLibrary`, primer3-py `mishyb_lib`: name → sequence, optional `*weight`, IUPAC codes, Primer3 adds each entry's reverse complement), `MaxLibraryMishyb` (PRIMER_INTERNAL_MAX_LIBRARY_MISHYB, 12), `WeightLibraryMishyb` (PRIMER_INTERNAL_WT_LIBRARY_MISHYB, 0), `LibraryAmbiguityCodesConsensus` (PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS, 0). Every window passing the structure limits is scored by `libprimer3.c` `oligo_repeat_library_mispriming` (OT_INTL): $w_i = weight_i \cdot align(probe, y_i)$, dpal unanchored local alignment (`local` / `local_ambig`: +1 match, −1 mismatch, −0.25 vs N, −2 per single-base gap, max gap 1; an entry < 3 nt scores its length; `PrimerDesigner.CalculateLibraryMishyb`, the same dpal port as the primer mispriming library); a window with any $w_i$ > (short) `MaxLibraryMishyb` is rejected and — Primer3 OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ is a "five-prime problem" — ends the 5′ extension of that 3′ end; the penalty gains `WeightLibraryMishyb` × score (`p_obj_fn`). `Primer3Probe.LibraryMishyb` / `LibraryMishybName` report Primer3's `repeat_sim.max` score and entry (PRIMER_INTERNAL_n_LIBRARY_MISHYB; primer3-py key `…_LIBRARY_MISPRIMING`). A weight without a library and, in alignment mode, a limit > 32767 throw (`_pr_data_control`). The same core serves `PrimerDesigner.DesignPrimers` with PRIMER_PICK_INTERNAL_OLIGO (`Primer_Design.md` §2.2 item 7). Cross-check: primer3-py 2.3.1 `design_primers(mishyb_lib=…)` with PRIMER_TASK = pick_hyb_probe_only: 814/814 random 40–220-nt templates (7875 probes; random 1–8-entry libraries — template fragments, mutations, IUPAC codes, 1–2-nt and N-rich entries, `*w` weights — PRIMER_INTERNAL_MAX_LIBRARY_MISHYB ∈ {6, 8, 10, 12, 12.9, 15, 20}, PRIMER_INTERNAL_WT_LIBRARY_MISHYB ∈ {0, 0.1, 0.5, 1}, both PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT modes, PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS ∈ {0, 1}, PRIMER_NUM_RETURN 5–50; probes rejected by the library in 699 templates) identical on position, PRIMER_INTERNAL_k_PENALTY (|Δ| ≤ 1e-9) and the library score + entry name.

Fraction bound (PRIMER_ANNEALING_TEMP; audit round 3, A3-5 part 1): `Primer3ProbeSettings.AnnealingTemperature` (default −10 = off; ≤ 100), `MinBound` / `MaxBound` / `OptBound` (PRIMER_INTERNAL_MIN/MAX/OPT_BOUND, −10 / 110 / 97 %), `WeightBoundGt` / `WeightBoundLt` (PRIMER_INTERNAL_WT_BOUND_GT/_LT, 0). With $T_a$ > 0 every window gets Primer3's `oligotm` fraction bound at $T_a$ (`PrimerDesigner.CalculateFractionBoundPrimer3` at the probe conditions; `Primer3Probe.Bound` = PRIMER_INTERNAL_n_BOUND) and is rejected outside [MinBound, MaxBound] right after the Tm check (not a five-prime problem). The bound terms of the internal-oligo `p_obj_fn` are not gated by $T_a$: without it (or for a window > 36 nt) the bound is Primer3's OLIGOTM_ERROR −999999.9999, so `WeightBoundLt` w adds w·(OptBound + 999999.9999), reproduced. OptBound outside [MinBound, MaxBound] or $T_a$ > 100 → `ArgumentOutOfRangeException` (`_pr_data_control`). For the internal oligo of a primer pair the annealing temperature is `PrimerParameters.AnnealingTemperature` (one global Primer3 setting). Cross-check: 85 pick_hyb_probe_only lists (432 probes) within the F44 harness identical to primer3-py 2.3.1 (position, penalty, bound).

Sequence quality (SEQUENCE_QUALITY; audit round 3, A3-5 part 2a): `DesignProbesPrimer3(template, settings, numReturn, sequenceQuality)` takes one integer quality per template base; every window gets `Primer3Probe.MinSequenceQuality` (PRIMER_INTERNAL_n_MIN_SEQ_QUALITY = min(PRIMER_QUALITY_RANGE_MAX, qualities over the window), `PrimerDesigner.CalculateSequenceQualityPrimer3`), is rejected when it is below `Primer3ProbeSettings.MinQuality` (PRIMER_INTERNAL_MIN_QUALITY, checked after the GC check; a five-prime problem) and gets `WeightSequenceQuality` × (`QualityRangeMax` − min quality) (PRIMER_INTERNAL_WT_SEQ_QUAL, the last internal-oligo `p_obj_fn` term). PRIMER_INTERNAL_WT_END_QUAL (`WeightEndQuality`) has no effect (never read by Primer3 2.3.1). `_pr_data_control`: quality length ≠ template, MinQuality ≠ 0 without quality or outside [`QualityRangeMin`, `QualityRangeMax`] (0 / 100), a value outside the range, a weight without quality → `ArgumentException`. For the internal oligo of a primer pair the quality comes from `PrimerPairOptions.SequenceQuality` and the range from `PrimerParameters.QualityRangeMin/Max`. Cross-check: the F45 harness (primer3-py 2.3.1, probe-only and PRIMER_PICK_INTERNAL_OLIGO cases, position / penalty / MIN_SEQ_QUALITY identical).

### 4.5 Oligo properties

- Molecular weight: canonical `SequenceStatistics.CalculateNucleotideMolecularWeight` (= Biopython `molecular_weight`, single-stranded; `CalculateMolecularWeight(seq)` infers RNA when the oligo has U and no T).
- ε260: `CalculateExtinctionCoefficient` = sum of mononucleotide values (no hypochromicity); `CalculateExtinctionCoefficientNearestNeighbor` = nearest-neighbour model ε = Σ ε(NᵢNᵢ₊₁) − Σ ε(internal Nᵢ) with the Cantor, Warshaw & Shapiro (1970) DNA / Warshaw & Tinoco (1966) RNA tables [11].
- Molecular beacon: loop Tm window [T+7, T+10] °C and stem-loop (ntthal hairpin) Tm ≥ T+7 °C for a detection temperature T (Tyagi & Kramer 1996 design rules [12]); without T the window is 55–65 °C.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DesignProbes` | `O(n × m)` | `O(k)` | `n` is sequence length, `m` is the scanned length range, and `k` is the number of retained candidates |
| `CheckSpecificity` | `O(m)` | `O(1)` | Uses suffix-tree lookups per probe |
| `DesignTilingProbes` | `O(n)` over fixed-length windows | `O(k)` | Produces overlapping probes for coverage |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProbeDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs), [ThermoConstants.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/ThermoConstants.cs)

- `ProbeDesigner.DesignProbes(string, ProbeParameters?, int)`: Main probe-generation and ranking routine.
- `ProbeDesigner.DesignProbes(string, ISuffixTree, ProbeParameters?, int, bool)`: Uniqueness-aware overload using a suffix tree.
- `ProbeDesigner.DesignTilingProbes(...)`: Generates overlapping tiling probes for coverage.
- `ProbeDesigner.CheckSpecificity(string, ISuffixTree)`: Maps suffix-tree hit counts to a specificity score.
- `ProbeDesigner.EvaluateTaqManProbe(string, double?, int, int)`: Opt-in TaqMan rule check; returns a `TaqManProbeEvaluation` with one boolean per rule and a `PassesAll` conjunction.
- `ProbeDesigner.SelectTaqManStrand(string, double?)`: Chooses the sense strand or its reverse complement, whichever better satisfies the TaqMan rules (no 5'-G, more C than G first).
- `ProbeDesigner.DesignProbesPrimer3(string, Primer3ProbeSettings?, int, IReadOnlyList<int>?)`: Primer3 `pick_hyb_probe_only` picker (§4.4), incl. the mishybridization library.
- `PrimerDesigner.CalculateLibraryMishyb(string, PrimerMisprimingLibrary, bool)`: Primer3 internal-oligo library mishybridization score and entry.

### 5.2 Current Behavior

The implementation evaluates candidates with prefix-sum GC optimization and begins from a raw-score-positive shortlist. Probe sequences are uppercased before evaluation. The suffix-tree overload does not rescore the full candidate universe: it rechecks only a larger raw-score shortlist, currently `maxProbes * 5`, can enforce uniqueness on that shortlist, or scales shortlisted scores by `1 / hitCount`; if `requireUnique` is `false` and specificity is `0`, a shortlisted probe can remain in the final output with score `0`, so final ordering still inherits the initial raw-score pass. `DesignTilingProbes(...)` includes suboptimal probes when needed for coverage and reports coverage, mean Tm, and Tm range. The source also defines probe types `Standard`, `Tiling`, `Antisense`, `LNA`, and `MolecularBeacon`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Application-specific probe-length, GC, and Tm windows.
- Heuristic penalties for GC, Tm, self-complementarity, secondary structure, and repeats.
- Genome-index-based uniqueness checking through a suffix tree.
- Opt-in TaqMan rules (no 5'-G, more C than G, no ≥4-G run, GC 30–80%, length 18–22 nt, probe Tm ≥ primer Tm + 10 °C) and strand selection. [7][8][9]

**Intentionally simplified:**

- `DesignProbes` ranks with a fixed additive penalty score (the published Primer3 objective is available as `DesignProbesPrimer3`); **consequence:** scores rank candidates but are not hybridization probabilities.
- Probes > 60 nt (Northern/Southern/FISH presets) use the sequence-only self-structure screens: ntthal (thal.c THAL_MAX_ALIGN = 60) and Primer3 do not handle longer oligos.
- Genome-index specificity is applied only after an initial raw-score shortlist is formed; **consequence:** uniqueness-aware results are specificity-filtered or specificity-scaled subsets of the top raw-score candidates rather than a full-candidate rerank.

**Not implemented:**

- Database-style alignment or experimentally calibrated hybridization prediction; **users should rely on:** external probe-validation workflows when those are required.
- MGB (minor-groove binder), LNA, and dual-quencher probe chemistries; **users should rely on:** the relevant chemistry's own design tool for those. The TaqMan rules implemented here target standard (single reporter/quencher) hydrolysis probes.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Probe Tm = Primer3 `seqtm` at Primer3 internal-oligo conditions | Assumption | Primer Express' own Tm scale (on which the 68–70 °C TaqMan window is defined) is proprietary; Primer3 is the open reference | accepted | Conditions are configurable in `ProbeParameters` / TaqMan arguments |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null or empty target sequence | Returns no probes | Explicit early return in source |
| Sequence shorter than `MinLength` | Returns no probes | No valid candidate window exists |
| Candidate with score `<= 0` | Rejected | The evaluator returns `null` for non-positive scores |
| Specificity check with no genome hits | Returns `0` | The probe does not match the indexed genome |
| `DesignProbesPrimer3` with `numReturn` < 1 | `ArgumentOutOfRangeException` | Primer3 `_pr_data_control` "PRIMER_NUM_RETURN < 1" |
| `DesignProbesPrimer3` with `WeightLibraryMishyb` ≠ 0 and no mishyb library | `ArgumentException` | Primer3 "Internal oligo mispriming score is part of objective function while mishyb library is not defined" |

### 6.2 Limitations

The current implementation uses heuristic penalties, simple self-structure detection, and shared Tm helpers instead of a full thermodynamic or database-backed specificity model. It is suitable for fast candidate generation and filtering, but not for high-confidence experimental validation by itself.

## 8. References

1. Wikipedia. "Nucleic acid thermodynamics." https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
2. Wikipedia. "Hybridization probe." https://en.wikipedia.org/wiki/Hybridization_probe
3. Wikipedia. "Fluorescence in situ hybridization." https://en.wikipedia.org/wiki/Fluorescence_in_situ_hybridization
4. Wikipedia. "DNA microarray." https://en.wikipedia.org/wiki/DNA_microarray
5. SantaLucia, J. (1998). "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." PNAS 95(4):1460-5.
6. Breslauer, K.J. et al. (1986). "Predicting DNA Duplex Stability from the Base Sequence." PNAS 83:3746-3750.
7. PREMIER Biosoft. "TaqMan® probe design tips." http://www.premierbiosoft.com/tech_notes/TaqMan.html (accessed 2026-06-24).
8. Applied Biosystems / Thermo Fisher Scientific. "Designing a TaqMan Gene Expression Assay." https://www.thermofisher.com/us/en/home/life-science/pcr/real-time-pcr/real-time-pcr-learning-center/gene-expression-analysis-real-time-pcr-information/designing-taqman-gene-expression-assay.html (accessed 2026-06-24).
9. ScienceDirect Topics. "TaqMan — an overview." https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/taqman (accessed 2026-06-24).
10. Applied Biosystems. Primer Express / "Designing TaqMan® MGB Probe and Primer Sets" guidelines: probe Tm 68–70 °C, no 5′ G, more C than G, avoid ≥ 4 G runs, G+C 30–80 %, probe ≥ 13 nt (via WebSearch extracts, 2026-10-01).
11. Cantor CR, Warshaw MM, Shapiro H (1970) Biopolymers 9:1059–1077; Warshaw MM, Tinoco I (1966) J Mol Biol 20:29–38 (ε260 nearest-neighbour tables as tabulated by OligoCalc / ATDBio; values confirmed via WebSearch extracts).
12. Tyagi S, Kramer FR (1996) Nat Biotechnol 14:303–308; Marras SAE et al. molecular-beacon design rules (stem 5–7 bp; probe and stem Tm 7–10 °C above the detection temperature).
13. Rozen S, Skaletsky H (2000) Primer3 on the WWW; Untergasser A et al. (2012) Nucleic Acids Res 40:e115; primer3 `libprimer3.cc` (raw.githubusercontent.com/primer3-org/primer3).
